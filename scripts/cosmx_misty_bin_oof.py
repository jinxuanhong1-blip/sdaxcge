#!/usr/bin/env python3
"""Within-section FOV-holdout at the bin scale.

Each bin sits inside one FOV. CD8/NK and CLDN4/keratin are computed only from
cells in that bin, so holding out a FOV does not leak across the fold.
The model for a section is fit on that section's other FOVs.
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


def fit_predict(Xtr, ytr, Xte, alpha=1.0):
    mu = Xtr.mean(0)
    sd = Xtr.std(0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    ymu = float(ytr.mean())
    ysd = float(ytr.std(ddof=1) or 1.0)
    model = Ridge(alpha=alpha)
    model.fit((Xtr - mu) / sd, (ytr - ymu) / ysd)
    return model.predict((Xte - mu) / sd) * ysd + ymu


def make_bins(xy, ct, sample, fov, logn, pitch, min_tumor=8):
    tumor = np.isin(ct, TUMOR)
    cd8 = np.isin(ct, CD8NK)
    cldn = logn[:, 0]
    krt = logn[:, 1:].mean(1)
    frames = []
    for s in SAMPLES:
        sm = sample == s
        origin = xy[sm].min(0)
        # section quantile of tumor CLDN4, used only as a fixed transform label.
        # Recomputed inside folds for the selected model; here for the search table.
        tmask = sm & tumor
        q75 = np.quantile(cldn[tmask], 0.75)
        med = np.median(cldn[tmask])
        for fv in np.unique(fov[sm]):
            m = sm & (fov == fv)
            if m.sum() < 40:
                continue
            xy_f = xy[m]
            ix = np.floor((xy_f[:, 0] - origin[0]) / pitch).astype(np.int32)
            iy = np.floor((xy_f[:, 1] - origin[1]) / pitch).astype(np.int32)
            key = ix * 10_000 + iy
            ct_f = ct[m]
            cldn_f = cldn[m]
            krt_f = krt[m]
            is_t = np.isin(ct_f, TUMOR)
            is_c = np.isin(ct_f, CD8NK)
            for k in np.unique(key):
                sel = key == k
                nt = int(is_t[sel].sum())
                if nt < min_tumor:
                    continue
                n_all = int(sel.sum())
                n_cd = int(is_c[sel].sum())
                ctum = cldn_f[sel & is_t]
                frames.append(
                    {
                        "sample": s,
                        "fov": int(fv),
                        "n_tumor": nt,
                        "n_all": n_all,
                        "n_cd8": n_cd,
                        "y_log": np.log1p(n_cd / (pitch ** 2) * 1e6),
                        "y_frac": n_cd / n_all,
                        "krt": float(krt_f[sel & is_t].mean()),
                        "cldn": float(ctum.mean()),
                        "cldn2": float(np.mean(ctum ** 2)),
                        "frac_q75": float((ctum >= q75).mean()),
                        "frac_med": float((ctum >= med).mean()),
                    }
                )
    return pd.DataFrame(frames)


def oof_within_section(df, ycol, kcols, fcols):
    parts = []
    per = []
    for s, g in df.groupby("sample"):
        g = g.reset_index(drop=True)
        fovs = g["fov"].to_numpy()
        n_fov = len(np.unique(fovs))
        if n_fov < 4 or len(g) < 40:
            continue
        n_splits = min(5, n_fov)
        y = g[ycol].to_numpy(float)
        Xk = g[kcols].to_numpy(float)
        Xf = g[fcols].to_numpy(float)
        oof_k = np.full(len(g), np.nan)
        oof_f = np.full(len(g), np.nan)
        for tr, te in GroupKFold(n_splits=n_splits).split(g, y, fovs):
            oof_k[te] = fit_predict(Xk[tr], y[tr], Xk[te])
            oof_f[te] = fit_predict(Xf[tr], y[tr], Xf[te])
        m = np.isfinite(oof_f)
        parts.append((y[m], oof_k[m], oof_f[m]))
        per.append((s, r2(y[m], oof_k[m]), r2(y[m], oof_f[m]), int(m.sum())))
    num_k = num_f = den = 0.0
    for y, pk, pf in parts:
        yc = y - y.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    rk = 1 - num_k / den
    rf = 1 - num_f / den
    return rk, rf, rf - rk, per


def main():
    print("loading", flush=True)
    xy, ct, sample, fov, logn = load()
    rows = []
    for pitch in (80, 120, 160, 220, 300):
        df = make_bins(xy, ct, sample, fov, logn, pitch)
        print(f"pitch {pitch} bins {len(df)}", flush=True)
        specs = [
            ("y_log", ["krt"], ["krt", "cldn"], "mean"),
            ("y_log", ["krt"], ["krt", "frac_q75"], "q75"),
            ("y_log", ["krt"], ["krt", "frac_med"], "ge_med"),
            ("y_log", ["krt"], ["krt", "cldn", "cldn2"], "quad"),
            ("y_log", ["krt"], ["krt", "cldn", "frac_q75"], "mean_q75"),
            ("y_frac", ["krt"], ["krt", "cldn"], "frac_mean"),
            ("y_frac", ["krt"], ["krt", "frac_q75"], "frac_q75"),
            ("y_frac", ["krt"], ["krt", "cldn", "frac_q75"], "frac_both"),
        ]
        for ycol, kcols, fcols, name in specs:
            rk, rf, delta, per = oof_within_section(df, ycol, kcols, fcols)
            rows.append(
                {
                    "pitch_um": pitch,
                    "spec": name,
                    "r2_krt": rk,
                    "r2_full": rf,
                    "delta": delta,
                    "median_section_delta": float(np.median([p[2] - p[1] for p in per])),
                    "n_sections_delta_pos": int(sum((p[2] - p[1]) > 0 for p in per)),
                }
            )
            print(
                f"  {name}: Δ={delta:+.3f} full={rf:.3f} krt={rk:.3f} "
                f"med={rows[-1]['median_section_delta']:+.3f} "
                f"{rows[-1]['n_sections_delta_pos']}/8",
                flush=True,
            )
            if delta > 0.08:
                print("   ", [(p[0], round(p[2] - p[1], 3), round(p[2], 3)) for p in per], flush=True)
    out = pd.DataFrame(rows).sort_values("delta", ascending=False)
    path = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt", "tables", "bin_oof_search.csv")
    out.to_csv(path, index=False)
    print(out.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
