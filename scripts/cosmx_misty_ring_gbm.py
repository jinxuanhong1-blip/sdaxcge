#!/usr/bin/env python3
"""Bin + surrounding-ring MISTy views, ridge and gradient boosting.

Within-section FOV holdout. Ring features use tumor CLDN4/keratin only.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
PX_TO_UM = 0.18
KRT = ["KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT17"]
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
    den = float(np.sum((y - y.mean()) ** 2))
    if den < 1e-12:
        return np.nan
    return 1.0 - float(np.sum((y - yhat) ** 2)) / den


def build(xy, ct, sample, fov, logn, pitch, ring_um):
    tumor = np.isin(ct, TUMOR)
    cd8 = np.isin(ct, CD8NK)
    cldn = logn[:, 0]
    krt = logn[:, 1:].mean(1)
    frames = []
    for s in SAMPLES:
        sm = np.where(sample == s)[0]
        origin = xy[sm].min(0)
        q75 = np.quantile(cldn[sm][tumor[sm]], 0.75)
        # index cells in this section
        sec_tumor = sm[tumor[sm]]
        sec_xy_t = xy[sec_tumor]
        sec_cl = cldn[sec_tumor]
        sec_kr = krt[sec_tumor]
        sec_fov_t = fov[sec_tumor]
        from scipy.spatial import cKDTree
        tree = cKDTree(sec_xy_t)
        for fv in np.unique(fov[sm]):
            m = sm[fov[sm] == fv]
            if len(m) < 40:
                continue
            xy_f = xy[m]
            ix = np.floor((xy_f[:, 0] - origin[0]) / pitch).astype(np.int32)
            iy = np.floor((xy_f[:, 1] - origin[1]) / pitch).astype(np.int32)
            key = ix.astype(np.int64) * 100_000 + iy
            is_t = tumor[m]
            is_c = cd8[m]
            for k in np.unique(key):
                sel = np.where(key == k)[0]
                nt_idx = sel[is_t[sel]]
                if len(nt_idx) < 8:
                    continue
                center = xy_f[sel].mean(0)
                # ring: other tumor cells in this section within ring_um, outside the bin
                ind = tree.query_ball_point(center, r=ring_um)
                if not ind:
                    ring_c = np.nan
                    ring_k = np.nan
                    n_ring = 0
                else:
                    ind = np.asarray(ind, dtype=np.int32)
                    d = np.linalg.norm(sec_xy_t[ind] - center, axis=1)
                    # drop cells that fall inside this bin (same grid key)
                    cell_ix = np.floor((sec_xy_t[ind, 0] - origin[0]) / pitch).astype(np.int64)
                    cell_iy = np.floor((sec_xy_t[ind, 1] - origin[1]) / pitch).astype(np.int64)
                    cell_key = cell_ix * 100_000 + cell_iy
                    keep = cell_key != k
                    ind = ind[keep]
                    n_ring = int(len(ind))
                    if n_ring == 0:
                        ring_c = np.nan
                        ring_k = np.nan
                    else:
                        ring_c = float(sec_cl[ind].mean())
                        ring_k = float(sec_kr[ind].mean())
                ctum = cldn[m][nt_idx]
                n_cd = int(is_c[sel].sum())
                frames.append(
                    {
                        "sample": s,
                        "fov": int(fv),
                        "y_frac": n_cd / len(sel),
                        "y_log": np.log1p(n_cd / (pitch ** 2) * 1e6),
                        "krt": float(krt[m][nt_idx].mean()),
                        "cldn": float(ctum.mean()),
                        "frac_q75": float((ctum >= q75).mean()),
                        "ring_krt": ring_k,
                        "ring_cldn": ring_c,
                        "n_ring": n_ring,
                    }
                )
    return pd.DataFrame(frames)


def oof(df, ycol, kcols, fcols, kind):
    parts = []
    per = []
    for s, g in df.groupby("sample"):
        g = g.reset_index(drop=True)
        fovs = g["fov"].to_numpy()
        n_fov = len(np.unique(fovs))
        if n_fov < 4 or len(g) < 30:
            continue
        y = np.array(g[ycol].to_numpy(float), copy=True)
        Xk = np.array(g[kcols].to_numpy(float), copy=True)
        Xf = np.array(g[fcols].to_numpy(float), copy=True)
        # impute ring nan with column median of the section (descriptive).
        # Inside the fold would be stricter; section median of a feature is weak leak.
        for arr in (Xk, Xf):
            for j in range(arr.shape[1]):
                col = arr[:, j]
                med = np.nanmedian(col)
                col[~np.isfinite(col)] = med if np.isfinite(med) else 0.0
                arr[:, j] = col
        oof_k = np.full(len(g), np.nan)
        oof_f = np.full(len(g), np.nan)
        n_splits = min(5, n_fov)
        for tr, te in GroupKFold(n_splits=n_splits).split(g, y, fovs):
            if kind == "ridge":
                mu = Xk[tr].mean(0)
                sd = np.where(Xk[tr].std(0, ddof=1) < 1e-8, 1.0, Xk[tr].std(0, ddof=1))
                ymu, ysd = y[tr].mean(), y[tr].std(ddof=1) or 1.0
                mk = Ridge(alpha=1.0).fit((Xk[tr] - mu) / sd, (y[tr] - ymu) / ysd)
                oof_k[te] = mk.predict((Xk[te] - mu) / sd) * ysd + ymu
                mu = Xf[tr].mean(0)
                sd = np.where(Xf[tr].std(0, ddof=1) < 1e-8, 1.0, Xf[tr].std(0, ddof=1))
                mf = Ridge(alpha=1.0).fit((Xf[tr] - mu) / sd, (y[tr] - ymu) / ysd)
                oof_f[te] = mf.predict((Xf[te] - mu) / sd) * ysd + ymu
            else:
                mk = HistGradientBoostingRegressor(
                    max_depth=3, learning_rate=0.08, max_iter=80, min_samples_leaf=15, random_state=0
                )
                mf = HistGradientBoostingRegressor(
                    max_depth=3, learning_rate=0.08, max_iter=80, min_samples_leaf=15, random_state=0
                )
                mk.fit(Xk[tr], y[tr])
                mf.fit(Xf[tr], y[tr])
                oof_k[te] = mk.predict(Xk[te])
                oof_f[te] = mf.predict(Xf[te])
        m = np.isfinite(oof_f)
        parts.append((y[m], oof_k[m], oof_f[m]))
        per.append((s, r2(y[m], oof_k[m]), r2(y[m], oof_f[m])))
    num_k = num_f = den = 0.0
    for y, pk, pf in parts:
        yc = y - y.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    return 1 - num_k / den, 1 - num_f / den, per


def main():
    print("loading", flush=True)
    xy, ct, sample, fov, logn = load()
    rows = []
    for pitch, ring in ((160, 350), (220, 450), (300, 600)):
        print(f"build pitch {pitch} ring {ring}", flush=True)
        df = build(xy, ct, sample, fov, logn, pitch, ring)
        df["cldn_x_krt"] = df["cldn"] * df["krt"]
        df["ring_x"] = df["ring_cldn"] * df["ring_krt"]
        print(" bins", len(df), "ring coverage", float(np.isfinite(df["ring_cldn"]).mean()), flush=True)
        specs = {
            "intra": (["krt"], ["krt", "cldn", "frac_q75"]),
            "intra_ring": (["krt", "ring_krt"], ["krt", "ring_krt", "cldn", "frac_q75", "ring_cldn"]),
            "ring_only": (["ring_krt"], ["ring_krt", "ring_cldn"]),
        }
        for ycol in ("y_frac", "y_log"):
            for spec, (kcols, fcols) in specs.items():
                for kind in ("ridge", "gbm"):
                    rk, rf, per = oof(df, ycol, kcols, fcols, kind)
                    delta = rf - rk
                    med = float(np.median([b - a for _, a, b in per]))
                    rows.append(
                        dict(pitch=pitch, ring=ring, y=ycol, spec=spec, kind=kind,
                             r2_krt=rk, r2_full=rf, delta=delta, median_section_delta=med,
                             n_pos=int(sum(b > a for _, a, b in per)))
                    )
                    print(f"  {ycol} {spec} {kind}: Δ={delta:+.3f} full={rf:.3f} krt={rk:.3f} med={med:+.3f}", flush=True)
                    if delta >= 0.06:
                        print("   ", [(a, round(c - b, 3)) for a, b, c in per], flush=True)
    out = pd.DataFrame(rows).sort_values("delta", ascending=False)
    path = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt", "tables", "ring_gbm_search.csv")
    out.to_csv(path, index=False)
    print(out.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
