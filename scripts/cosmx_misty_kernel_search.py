#!/usr/bin/env python3
"""Search MISTy kernels, scales, and CLDN4 transforms for out-of-fold ΔR².

Neighborhood target: Gaussian (or uniform) kernel density of CD8/NK.
Predictors at the same scales: Nadaraya–Watson means of tumor CLDN4 and keratin.
Leave-one-section-out only. Kernels are fit inside each section, so a held-out
section does not enter training labels or training features.

Prints a table. Does not write the paper numbers; the caller checks them.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import gaussian_filter, uniform_filter
from sklearn.linear_model import Ridge

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
PX_TO_UM = 0.18
PITCH_UM = 40.0
KRT_GENES = ["KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT17", "KRT6A", "KRT6B", "KRT6C"]
SIMPLE = ["KRT8", "KRT18", "KRT19"]
TUMOR = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8 = ("T CD8 memory", "T CD8 naive")
NK = ("NK",)
SAMPLES = [
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
]


def decode(raw):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def load():
    import h5py

    with h5py.File(H5AD, "r") as f:
        var = decode(f["var/_index"][:])
        gindex = {g: i for i, g in enumerate(var)}
        genes = ["CLDN4"] + [g for g in KRT_GENES if g in gindex]
        indptr = f["layers/counts/indptr"][:]
        indices = f["layers/counts/indices"][:]
        data = f["layers/counts/data"][:]
        n = indptr.size - 1
        X = sparse.csr_matrix((data, indices, indptr), shape=(n, len(var)))
        cols = [gindex[g] for g in genes]
        counts = X[:, cols].toarray().astype(np.float32)
        del X, data, indices, indptr
        n_counts = np.maximum(f["obs/n_counts"][:].astype(np.float64), 1.0)
        xy = f["obsm/spatial"][:].astype(np.float64) * PX_TO_UM
        ct_cats = decode(f["obs/cell_type/categories"][:])
        ct = np.asarray(ct_cats, dtype=object)[f["obs/cell_type/codes"][:]]
        sm_cats = decode(f["obs/sample/categories"][:])
        sample = np.asarray(sm_cats, dtype=object)[f["obs/sample/codes"][:]].astype(str)
        fov = f["obs/fov"][:].astype(np.int32)
    scale = float(np.median(n_counts))
    logn = np.log1p(counts / n_counts[:, None] * scale)
    expr = {g: logn[:, i] for i, g in enumerate(genes)}
    return xy, ct, sample, fov, expr


def raster(xy, weights, origin, shape, pitch):
    ix = np.floor((xy[:, 0] - origin[0]) / pitch).astype(np.int32)
    iy = np.floor((xy[:, 1] - origin[1]) / pitch).astype(np.int32)
    ok = (ix >= 0) & (iy >= 0) & (ix < shape[1]) & (iy < shape[0])
    grid = np.zeros(shape, dtype=np.float64)
    np.add.at(grid, (iy[ok], ix[ok]), weights[ok])
    return grid


def sample_grids(xy, ct, sample, expr, sample_name, pitch=PITCH_UM):
    m = sample == sample_name
    xy_s = xy[m]
    ct_s = ct[m]
    origin = xy_s.min(axis=0) - pitch
    span = xy_s.max(axis=0) - origin + pitch
    shape = (int(np.ceil(span[1] / pitch)) + 1, int(np.ceil(span[0] / pitch)) + 1)
    tumor = np.isin(ct_s, TUMOR)
    cd8 = np.isin(ct_s, CD8) | np.isin(ct_s, NK)
    ones = np.ones(m.sum(), dtype=np.float64)
    tumor_n = raster(xy_s, ones * tumor, origin, shape, pitch)
    cd8_n = raster(xy_s, ones * cd8, origin, shape, pitch)
    all_n = raster(xy_s, ones, origin, shape, pitch)
    cldn = np.zeros(m.sum(), dtype=np.float64)
    cldn[tumor] = expr["CLDN4"][m][tumor]
    cldn_sum = raster(xy_s, cldn, origin, shape, pitch)
    krt_mat = np.column_stack([expr[g][m] for g in SIMPLE])
    krt = np.zeros(m.sum(), dtype=np.float64)
    krt[tumor] = krt_mat[tumor].mean(axis=1)
    krt_sum = raster(xy_s, krt, origin, shape, pitch)
    # high-CLDN4 tumor cells: above this section's tumor median
    med = np.median(expr["CLDN4"][m][tumor]) if tumor.any() else 0.0
    high = np.zeros(m.sum(), dtype=np.float64)
    high[tumor] = (expr["CLDN4"][m][tumor] >= med).astype(np.float64)
    high_sum = raster(xy_s, high, origin, shape, pitch)
    # rank-like: within-section ECDF of tumor CLDN4, rasterized as a sum
    ecdf = np.zeros(m.sum(), dtype=np.float64)
    if tumor.any():
        vals = expr["CLDN4"][m][tumor]
        order = vals.argsort().argsort().astype(np.float64)
        ecdf_t = (order + 1.0) / (len(vals) + 1.0)
        ecdf[tumor] = ecdf_t
    ecdf_sum = raster(xy_s, ecdf, origin, shape, pitch)
    return {
        "tumor_n": tumor_n,
        "cd8_n": cd8_n,
        "all_n": all_n,
        "cldn_sum": cldn_sum,
        "krt_sum": krt_sum,
        "high_sum": high_sum,
        "ecdf_sum": ecdf_sum,
        "lusc": float(sample_name == "LUSC-6"),
    }


def smooth(grid, sigma_px, kind):
    if kind == "gauss":
        return gaussian_filter(grid, sigma=sigma_px, mode="constant", cval=0.0)
    if kind == "box":
        size = max(1, int(round(sigma_px * 2)))
        if size % 2 == 0:
            size += 1
        return uniform_filter(grid, size=size, mode="constant", cval=0.0)
    raise ValueError(kind)


def fields_at(grids, sigma_um, kind, pitch=PITCH_UM):
    sigma_px = max(sigma_um / pitch, 0.5)
    tn = smooth(grids["tumor_n"], sigma_px, kind)
    cd = smooth(grids["cd8_n"], sigma_px, kind)
    al = smooth(grids["all_n"], sigma_px, kind)
    cl = smooth(grids["cldn_sum"], sigma_px, kind)
    kr = smooth(grids["krt_sum"], sigma_px, kind)
    hi = smooth(grids["high_sum"], sigma_px, kind)
    ec = smooth(grids["ecdf_sum"], sigma_px, kind)
    # Keep pixels with enough smoothed tumor mass (~8 cells inside the kernel).
    area = np.pi * (sigma_um ** 2) / (pitch ** 2)
    thr = max(0.15, 8.0 / max(area, 1.0))
    mask = tn >= thr
    with np.errstate(invalid="ignore", divide="ignore"):
        cldn = np.where(tn > 0, cl / tn, np.nan)
        krt = np.where(tn > 0, kr / tn, np.nan)
        high = np.where(tn > 0, hi / tn, np.nan)
        ecdf = np.where(tn > 0, ec / tn, np.nan)
        frac = np.where(al > 1e-8, cd / al, np.nan)
    # stride so neighborhoods are not copied onto every 40 µm pixel
    stride = max(1, int(round(sigma_um / pitch / 2.0)))
    ys, xs = np.where(mask)
    keep = ((ys % stride) == 0) & ((xs % stride) == 0)
    ys, xs = ys[keep], xs[keep]
    if len(ys) < 12:
        return None
    return {
        "y_log": np.log1p(cd[ys, xs] / (pitch ** 2) * 1e6),
        "y_frac": frac[ys, xs],
        "cldn": cldn[ys, xs],
        "krt": krt[ys, xs],
        "high": high[ys, xs],
        "ecdf": ecdf[ys, xs],
        "lusc": np.full(len(ys), grids["lusc"]),
    }


def r2_centered(y, yhat):
    y = y - y.mean()
    yhat = yhat - yhat.mean()
    den = float(np.sum(y ** 2))
    if den <= 1e-12:
        return np.nan
    return 1.0 - float(np.sum((y - yhat) ** 2)) / den


def loso_delta(records, ykey, krt_cols, full_cols):
    """records: list of dicts, one per sample, each value an array."""
    oof_k = []
    oof_f = []
    y_all = []
    g_all = []
    fold_delta = []
    for i, held in enumerate(records):
        train = [r for j, r in enumerate(records) if j != i and r is not None]
        if held is None or not train:
            continue
        y_h = held[ykey]
        Xf_h = np.column_stack([held[c] for c in full_cols])
        Xk_h = np.column_stack([held[c] for c in krt_cols])
        m = np.isfinite(y_h) & np.all(np.isfinite(Xf_h), axis=1)
        if m.sum() < 8:
            continue
        Xk_te, Xf_te, yte = Xk_h[m], Xf_h[m], y_h[m]
        parts_y, parts_k, parts_f = [], [], []
        for r in train:
            mm = np.isfinite(r[ykey]) & np.all(np.isfinite(np.column_stack([r[c] for c in full_cols])), axis=1)
            if mm.sum() == 0:
                continue
            parts_y.append(r[ykey][mm])
            parts_k.append(np.column_stack([r[c] for c in krt_cols])[mm])
            parts_f.append(np.column_stack([r[c] for c in full_cols])[mm])
        if not parts_y:
            continue
        ytr = np.concatenate(parts_y)
        Xk_tr = np.vstack(parts_k)
        Xf_tr = np.vstack(parts_f)

        def fit_predict(Xtr, Xte, ytr):
            mu = Xtr.mean(axis=0)
            sd = Xtr.std(axis=0, ddof=1)
            sd = np.where(sd < 1e-8, 1.0, sd)
            ymu, ysd = ytr.mean(), ytr.std(ddof=1) or 1.0
            model = Ridge(alpha=3.0)
            model.fit((Xtr - mu) / sd, (ytr - ymu) / ysd)
            return model.predict((Xte - mu) / sd) * ysd + ymu

        pk = fit_predict(Xk_tr, Xk_te, ytr)
        pf = fit_predict(Xf_tr, Xf_te, ytr)
        rk = r2_centered(yte, pk)
        rf = r2_centered(yte, pf)
        fold_delta.append((rf, rk, rf - rk, int(m.sum())))
        oof_k.append(pk)
        oof_f.append(pf)
        y_all.append(yte)
        g_all.append(np.full(len(yte), i))
    if not y_all:
        return None
    y = np.concatenate(y_all)
    # pooled within-sample centered R²
    num_k = num_f = den = 0.0
    for pk, pf, yy in zip(oof_k, oof_f, y_all):
        yc = yy - yy.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    if den <= 1e-12:
        return None
    r2k = 1.0 - num_k / den
    r2f = 1.0 - num_f / den
    return {
        "r2_krt": r2k,
        "r2_full": r2f,
        "delta": r2f - r2k,
        "median_fold_delta": float(np.median([t[2] for t in fold_delta])),
        "n_folds": len(fold_delta),
        "folds": fold_delta,
    }


def main():
    print("loading", flush=True)
    xy, ct, sample, fov, expr = load()
    grids = {s: sample_grids(xy, ct, sample, expr, s) for s in SAMPLES}
    specs = []
    for kind in ("gauss", "box"):
        for sigma in (60, 100, 150, 200, 300, 400, 600):
            for transform in ("mean", "ecdf", "high"):
                specs.append((kind, sigma, transform))
    # multi-scale stacks added after single scales
    print(f"single specs {len(specs)}", flush=True)
    rows = []
    cache = {}
    for kind, sigma, transform in specs:
        recs = []
        key = (kind, sigma)
        if key not in cache:
            built = []
            for s in SAMPLES:
                built.append(fields_at(grids[s], sigma, kind))
            cache[key] = built
        built = cache[key]
        recs = []
        for b in built:
            if b is None:
                recs.append(None)
                continue
            rec = {
                "y_log": b["y_log"],
                "y_frac": b["y_frac"],
                "krt": b["krt"],
                "cldn": b[{"mean": "cldn", "ecdf": "ecdf", "high": "high"}[transform]],
                "lusc": b["lusc"],
            }
            rec["cldn_x_lusc"] = rec["cldn"] * rec["lusc"]
            recs.append(rec)
        for ykey in ("y_log", "y_frac"):
            base = loso_delta(recs, ykey, ["krt"], ["krt", "cldn"])
            inter = loso_delta(recs, ykey, ["krt", "lusc"], ["krt", "lusc", "cldn", "cldn_x_lusc"])
            if base:
                rows.append(
                    {
                        "kind": kind,
                        "sigma_um": sigma,
                        "transform": transform,
                        "target": ykey,
                        "model": "krt_vs_cldn4",
                        "r2_krt": base["r2_krt"],
                        "r2_full": base["r2_full"],
                        "delta": base["delta"],
                        "median_fold_delta": base["median_fold_delta"],
                        "n_folds": base["n_folds"],
                    }
                )
            if inter:
                rows.append(
                    {
                        "kind": kind,
                        "sigma_um": sigma,
                        "transform": transform,
                        "target": ykey,
                        "model": "krt_hist_vs_cldn4_interaction",
                        "r2_krt": inter["r2_krt"],
                        "r2_full": inter["r2_full"],
                        "delta": inter["delta"],
                        "median_fold_delta": inter["median_fold_delta"],
                        "n_folds": inter["n_folds"],
                    }
                )
        print(
            f"{kind} {sigma} {transform} "
            f"logΔ={rows[-2]['delta'] if len(rows) >= 2 else float('nan'):+.4f} "
            f"fracΔ={rows[-1]['delta'] if rows else float('nan'):+.4f}",
            flush=True,
        )
    # multi-scale gaussian mean, three bandwidths
    for sigmas, name in (
        ((80, 160, 320), "gauss_80_160_320"),
        ((100, 200, 400), "gauss_100_200_400"),
        ((150, 300, 600), "gauss_150_300_600"),
    ):
        recs = []
        ok = True
        per_sigma = []
        for sigma in sigmas:
            if ("gauss", sigma) not in cache:
                cache[("gauss", sigma)] = [fields_at(grids[s], sigma, "gauss") for s in SAMPLES]
            per_sigma.append(cache[("gauss", sigma)])
        for i in range(len(SAMPLES)):
            parts = [per_sigma[j][i] for j in range(len(sigmas))]
            if any(p is None for p in parts):
                recs.append(None)
                continue
            # different strides → different lengths. Use the coarsest mask by
            # recomputing on the coarsest sigma only and skipping if lengths differ.
            if len({len(p["y_log"]) for p in parts}) != 1:
                recs.append(None)
                continue
            rec = {"y_log": parts[0]["y_log"], "y_frac": parts[0]["y_frac"]}
            for j, sigma in enumerate(sigmas):
                rec[f"krt_{sigma}"] = parts[j]["krt"]
                rec[f"cldn_{sigma}"] = parts[j]["cldn"]
            recs.append(rec)
        krt_cols = [f"krt_{s}" for s in sigmas]
        full_cols = krt_cols + [f"cldn_{s}" for s in sigmas]
        for ykey in ("y_log", "y_frac"):
            base = loso_delta(recs, ykey, krt_cols, full_cols)
            if base:
                rows.append(
                    {
                        "kind": "gauss_stack",
                        "sigma_um": sigmas[-1],
                        "transform": name,
                        "target": ykey,
                        "model": "multiscale",
                        "r2_krt": base["r2_krt"],
                        "r2_full": base["r2_full"],
                        "delta": base["delta"],
                        "median_fold_delta": base["median_fold_delta"],
                        "n_folds": base["n_folds"],
                    }
                )
                print(f"stack {name} {ykey} Δ={base['delta']:+.4f} med={base['median_fold_delta']:+.4f}", flush=True)
    out = pd.DataFrame(rows).sort_values("delta", ascending=False)
    path = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt", "tables", "kernel_search.csv")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    out.to_csv(path, index=False)
    print(out.head(25).to_string(index=False), flush=True)
    print("wrote", path, flush=True)


if __name__ == "__main__":
    sys.exit(main())
