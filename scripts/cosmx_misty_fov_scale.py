#!/usr/bin/env python3
"""FOV-scale and within-section ceiling for CLDN4 beyond keratin.

Units are whole FOVs (no cross-FOV kernel, so FOV holdout does not leak).
Also reports in-sample partial R² of the Gaussian neighborhood field, which is
the ceiling an out-of-fold model can approach.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
PX_TO_UM = 0.18
KRT = ["KRT8", "KRT18", "KRT19"]
TUMOR = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8NK = ("T CD8 memory", "T CD8 naive", "NK")
SAMPLES = [
    "LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3", "LUSC-6",
    "LUAD-9 R1", "LUAD-9 R2", "LUAD-12", "LUAD-13",
]


def decode(raw):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def load():
    import h5py

    with h5py.File(H5AD, "r") as f:
        var = decode(f["var/_index"][:])
        gindex = {g: i for i, g in enumerate(var)}
        genes = ["CLDN4"] + KRT
        indptr = f["layers/counts/indptr"][:]
        indices = f["layers/counts/indices"][:]
        data = f["layers/counts/data"][:]
        n = indptr.size - 1
        X = sparse.csr_matrix((data, indices, indptr), shape=(n, len(var)))
        counts = X[:, [gindex[g] for g in genes]].toarray().astype(np.float32)
        del X
        n_counts = np.maximum(f["obs/n_counts"][:].astype(np.float64), 1.0)
        xy = f["obsm/spatial"][:].astype(np.float64) * PX_TO_UM
        ct = np.asarray(decode(f["obs/cell_type/categories"][:]), dtype=object)[f["obs/cell_type/codes"][:]]
        sample = np.asarray(decode(f["obs/sample/categories"][:]), dtype=object)[f["obs/sample/codes"][:]].astype(str)
        fov = f["obs/fov"][:].astype(np.int32)
    logn = np.log1p(counts / n_counts[:, None] * float(np.median(n_counts)))
    return xy, ct, sample, fov, logn


def r2(y, yhat):
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    den = float(np.sum((y - y.mean()) ** 2))
    if den < 1e-12:
        return np.nan
    return 1.0 - float(np.sum((y - yhat) ** 2)) / den


def fit_predict(Xtr, ytr, Xte):
    mu = Xtr.mean(0)
    sd = Xtr.std(0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    ymu = ytr.mean()
    ysd = ytr.std(ddof=1) or 1.0
    model = Ridge(alpha=1.0)
    model.fit((Xtr - mu) / sd, (ytr - ymu) / ysd)
    return model.predict((Xte - mu) / sd) * ysd + ymu


def pooled_r2(parts):
    """parts: list of (y, yhat_k, yhat_f) centered inside each part."""
    num_k = num_f = den = 0.0
    for y, pk, pf in parts:
        yc = y - y.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    if den < 1e-12:
        return np.nan, np.nan, np.nan
    rk = 1 - num_k / den
    rf = 1 - num_f / den
    return rk, rf, rf - rk


def fov_table(xy, ct, sample, fov, logn):
    tumor = np.isin(ct, TUMOR)
    cd8 = np.isin(ct, CD8NK)
    cldn = logn[:, 0]
    krt = logn[:, 1:].mean(axis=1)
    rows = []
    for s in SAMPLES:
        sm = sample == s
        for fv in np.unique(fov[sm]):
            m = sm & (fov == fv)
            t = m & tumor
            if t.sum() < 30 or m.sum() < 80:
                continue
            area = 0.98 * 0.65  # mm^2, CosMx FOV
            n_cd8 = int((m & cd8).sum())
            c = cldn[t]
            # within-FOV high fraction uses the section median, computed below
            rows.append(
                {
                    "sample": s,
                    "fov": int(fv),
                    "n_tumor": int(t.sum()),
                    "n_cd8": n_cd8,
                    "y_log": np.log1p(n_cd8 / area),
                    "y_frac": n_cd8 / m.sum(),
                    "krt": float(krt[t].mean()),
                    "cldn": float(c.mean()),
                    "cldn_sd": float(c.std()),
                    "cldn_p90": float(np.quantile(c, 0.9)),
                }
            )
    df = pd.DataFrame(rows)
    # section-wise high fraction and within-section rank of FOV mean CLDN4
    df["high"] = np.nan
    df["ecdf"] = np.nan
    # recompute high from cells: need section median. Approximate with FOV means is weaker.
    # Do it properly below by a second pass stored in the row via cell data.
    return df, tumor, cd8, cldn


def add_cell_transforms(df, ct, sample, fov, logn):
    tumor = np.isin(ct, TUMOR)
    cldn = logn[:, 0]
    rows_high = []
    for s in SAMPLES:
        sm = (sample == s) & tumor
        med = np.median(cldn[sm])
        q75 = np.quantile(cldn[sm], 0.75)
        for fv, g in df[df["sample"] == s].groupby("fov"):
            m = sm & (fov == fv)
            rows_high.append((s, int(fv), float((cldn[m] >= med).mean()), float((cldn[m] >= q75).mean())))
    h = pd.DataFrame(rows_high, columns=["sample", "fov", "frac_ge_median", "frac_ge_q75"])
    return df.merge(h, on=["sample", "fov"])


def cv_fov(df, ycol, krt_cols, full_cols, groups):
    y = df[ycol].to_numpy(float)
    Xk = df[krt_cols].to_numpy(float)
    Xf = df[full_cols].to_numpy(float)
    g = df[groups].to_numpy()
    # if groups are sample, LOSO; if fov_key, GroupKFold
    uniq = pd.unique(g)
    if len(uniq) <= 8 and groups == "sample":
        splits = [(g != u, g == u) for u in uniq]
    else:
        n_splits = min(8, len(uniq))
        splits = []
        for tr, te in GroupKFold(n_splits=n_splits).split(df, y, g):
            mask_tr = np.zeros(len(df), dtype=bool)
            mask_te = np.zeros(len(df), dtype=bool)
            mask_tr[tr] = True
            mask_te[te] = True
            splits.append((mask_tr, mask_te))
    parts = []
    fold_delta = []
    for tr, te in splits:
        if tr.sum() < 10 or te.sum() < 4:
            continue
        pk = fit_predict(Xk[tr], y[tr], Xk[te])
        pf = fit_predict(Xf[tr], y[tr], Xf[te])
        # center within each test sample so a section mean is not required
        yy, pkc, pfc = [], [], []
        for s in pd.unique(df.loc[te, "sample"]):
            m = te & (df["sample"].to_numpy() == s)
            if m.sum() < 3:
                continue
            yy.append(y[m])
            pkc.append(pk[m.nonzero()[0] if False else np.where(m)[0] * 0])  # placeholder
        # simpler: group test indices
        te_idx = np.where(te)[0]
        pred_map_k = {i: p for i, p in zip(te_idx, pk)}
        pred_map_f = {i: p for i, p in zip(te_idx, pf)}
        deltas = []
        for s in pd.unique(df.loc[te, "sample"]):
            idx = te_idx[df.loc[te_idx, "sample"].to_numpy() == s]
            if len(idx) < 3:
                continue
            ys = y[idx]
            pks = np.array([pred_map_k[i] for i in idx])
            pfs = np.array([pred_map_f[i] for i in idx])
            parts.append((ys, pks, pfs))
            deltas.append(r2(ys, pfs) - r2(ys, pks))
        if deltas:
            fold_delta.append(float(np.mean(deltas)))
    rk, rf, delta = pooled_r2(parts)
    return rk, rf, delta, float(np.median(fold_delta)) if fold_delta else np.nan


def within_sample_fov_holdout(df, ycol, krt_cols, full_cols, seed=25976224):
    """Inside each section, hold out half the FOVs. Pool the held-out R²."""
    rng = np.random.default_rng(seed)
    parts = []
    per = []
    for s, g in df.groupby("sample"):
        idx = g.index.to_numpy()
        if len(idx) < 8:
            continue
        perm = rng.permutation(idx)
        mid = len(perm) // 2
        # two folds
        for te_idx, tr_idx in ((perm[:mid], perm[mid:]), (perm[mid:], perm[:mid])):
            ytr = df.loc[tr_idx, ycol].to_numpy(float)
            yte = df.loc[te_idx, ycol].to_numpy(float)
            Xk_tr = df.loc[tr_idx, krt_cols].to_numpy(float)
            Xk_te = df.loc[te_idx, krt_cols].to_numpy(float)
            Xf_tr = df.loc[tr_idx, full_cols].to_numpy(float)
            Xf_te = df.loc[te_idx, full_cols].to_numpy(float)
            pk = fit_predict(Xk_tr, ytr, Xk_te)
            pf = fit_predict(Xf_tr, ytr, Xf_te)
            parts.append((yte, pk, pf))
            per.append((s, r2(yte, pf) - r2(yte, pk), len(te_idx)))
    rk, rf, delta = pooled_r2(parts)
    return rk, rf, delta, per


def in_sample_partial(df, ycol, krt_cols, full_cols):
    rows = []
    for s, g in df.groupby("sample"):
        if len(g) < 8:
            continue
        y = g[ycol].to_numpy(float)
        Xk = np.column_stack([np.ones(len(g)), g[krt_cols].to_numpy(float)])
        Xf = np.column_stack([np.ones(len(g)), g[full_cols].to_numpy(float)])
        bk, *_ = np.linalg.lstsq(Xk, y, rcond=None)
        bf, *_ = np.linalg.lstsq(Xf, y, rcond=None)
        rk = r2(y, Xk @ bk)
        rf = r2(y, Xf @ bf)
        rows.append((s, rk, rf, rf - rk, len(g)))
    return rows


def main():
    print("loading", flush=True)
    xy, ct, sample, fov, logn = load()
    df, *_ = fov_table(xy, ct, sample, fov, logn)
    df = add_cell_transforms(df, ct, sample, fov, logn)
    df["fov_key"] = df["sample"] + "::" + df["fov"].astype(str)
    df["lusc"] = (df["sample"] == "LUSC-6").astype(float)
    df["cldn_x_lusc"] = df["cldn"] * df["lusc"]
    df["high_x_lusc"] = df["frac_ge_q75"] * df["lusc"]
    print("fovs", len(df), flush=True)
    specs = [
        ("y_log", ["krt"], ["krt", "cldn"], "mean"),
        ("y_log", ["krt"], ["krt", "frac_ge_q75"], "q75"),
        ("y_log", ["krt"], ["krt", "frac_ge_median"], "ge_median"),
        ("y_log", ["krt"], ["krt", "cldn_p90"], "p90"),
        ("y_frac", ["krt"], ["krt", "cldn"], "frac_mean"),
        ("y_frac", ["krt"], ["krt", "frac_ge_q75"], "frac_q75"),
        ("y_log", ["krt", "lusc"], ["krt", "lusc", "cldn", "cldn_x_lusc"], "interact_mean"),
        ("y_log", ["krt", "lusc"], ["krt", "lusc", "frac_ge_q75", "high_x_lusc"], "interact_q75"),
        ("y_frac", ["krt", "lusc"], ["krt", "lusc", "frac_ge_q75", "high_x_lusc"], "interact_q75_frac"),
    ]
    rows = []
    for ycol, kcols, fcols, name in specs:
        rk, rf, delta, med = cv_fov(df, ycol, kcols, fcols, "sample")
        wrk, wrf, wdelta, per = within_sample_fov_holdout(df, ycol, kcols, fcols)
        ins = in_sample_partial(df, ycol, kcols, fcols)
        med_in = float(np.median([t[3] for t in ins]))
        rows.append(
            {
                "spec": name,
                "loso_r2_krt": rk,
                "loso_r2_full": rf,
                "loso_delta": delta,
                "loso_median_section_delta": med,
                "within_r2_krt": wrk,
                "within_r2_full": wrf,
                "within_delta": wdelta,
                "insample_median_delta": med_in,
            }
        )
        print(f"{name}: LOSO Δ={delta:+.3f} (full {rf:.3f}) within-FOV-holdout Δ={wdelta:+.3f} (full {wrf:.3f}) in-sample median Δ={med_in:+.3f}", flush=True)
        print("  in-sample", [(a, round(d, 3)) for a, b, c, d, n in ins], flush=True)
    out = pd.DataFrame(rows)
    path = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt", "tables", "fov_scale_search.csv")
    out.to_csv(path, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
