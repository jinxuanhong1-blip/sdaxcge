#!/usr/bin/env python3
"""Neighborhood-scale MISTy: CLDN4 beyond keratin, FOV holdout.

Cell-level log1p counts gave an out-of-fold ΔR² near 0.02 because most tumor
cells have no CD8/NK neighbor. This script moves the same question to 220 µm
bins inside each FOV (the scale selected in the kernel/bin search).

Target: CD8+NK fraction of cells in the bin.
Baseline: mean keratin (KRT8/18/19/7/5/17) in the bin.
Added: mean CLDN4, and the fraction of tumor cells at or above the training
FOVs' 75th percentile of CLDN4.
Folds: 5-fold GroupKFold on FOV, fit separately in each section.
A bin never uses cells from another FOV, so the holdout does not leak.

Also reports leave-one-section-out on the same features, and a within-section
permutation null.
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
PX_TO_UM = 0.18
PITCH = 220.0
MIN_TUMOR = 8
KRT = ["KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT17"]
TUMOR = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8NK = ("T CD8 memory", "T CD8 naive", "NK")
SAMPLES = [
    "LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3", "LUSC-6",
    "LUAD-9 R1", "LUAD-9 R2", "LUAD-12", "LUAD-13",
]
PATIENT = {
    "LUAD-5 R1": "Lung5", "LUAD-5 R2": "Lung5", "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6", "LUAD-9 R1": "Lung9", "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12", "LUAD-13": "Lung13",
}
SEED = 25976224


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


def fit_predict(Xtr, ytr, Xte, alpha=1.0):
    mu = Xtr.mean(0)
    sd = Xtr.std(0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    ymu = float(ytr.mean())
    ysd = float(ytr.std(ddof=1) or 1.0)
    model = Ridge(alpha=alpha)
    model.fit((Xtr - mu) / sd, (ytr - ymu) / ysd)
    return model.predict((Xte - mu) / sd) * ysd + ymu, model


def make_bins(xy, ct, sample, fov, logn):
    tumor = np.isin(ct, TUMOR)
    cd8 = np.isin(ct, CD8NK)
    cldn = logn[:, 0]
    krt = logn[:, 1:].mean(1)
    rows = []
    cells = []
    for s in SAMPLES:
        sm = sample == s
        origin = xy[sm].min(0)
        for fv in np.unique(fov[sm]):
            m = sm & (fov == fv)
            if int(m.sum()) < 40:
                continue
            xy_f = xy[m]
            ix = np.floor((xy_f[:, 0] - origin[0]) / PITCH).astype(np.int64)
            iy = np.floor((xy_f[:, 1] - origin[1]) / PITCH).astype(np.int64)
            key = ix * 100_000 + iy
            is_t = tumor[m]
            is_c = cd8[m]
            cldn_f = cldn[m]
            krt_f = krt[m]
            for k in np.unique(key):
                sel = key == k
                nt = sel & is_t
                if int(nt.sum()) < MIN_TUMOR:
                    continue
                vals = cldn_f[nt].astype(np.float64)
                rows.append(
                    {
                        "sample": s,
                        "patient": PATIENT[s],
                        "fov": int(fv),
                        "n_tumor": int(nt.sum()),
                        "n_all": int(sel.sum()),
                        "n_cd8": int((sel & is_c).sum()),
                        "y": (sel & is_c).sum() / sel.sum(),
                        "krt": float(krt_f[nt].mean()),
                        "cldn": float(vals.mean()),
                    }
                )
                cells.append(vals)
    return pd.DataFrame(rows), cells


def design(df, cells, train_idx, test_idx):
    """Training-fold CLDN4 75th percentile. Returns Xk, Xf, y for train and test."""
    thr = np.quantile(np.concatenate([cells[i] for i in train_idx]), 0.75)
    frac = np.array([float((cells[i] >= thr).mean()) for i in range(len(cells))])
    y = df["y"].to_numpy(float)
    krt = df["krt"].to_numpy(float)
    cldn = df["cldn"].to_numpy(float)
    Xk = krt[:, None]
    Xf = np.column_stack([krt, cldn, frac])
    return Xk[train_idx], y[train_idx], Xf[train_idx], Xk[test_idx], y[test_idx], Xf[test_idx]


def cv_section(df, cells, sample_name):
    g = df[df["sample"] == sample_name].reset_index(drop=True)
    # cells aligned to full df; subset
    idx_full = np.where(df["sample"].to_numpy() == sample_name)[0]
    local_cells = [cells[i] for i in idx_full]
    fovs = g["fov"].to_numpy()
    n_fov = len(np.unique(fovs))
    n_splits = min(5, n_fov)
    y = g["y"].to_numpy(float)
    oof_k = np.full(len(g), np.nan)
    oof_f = np.full(len(g), np.nan)
    for tr, te in GroupKFold(n_splits=n_splits).split(g, y, fovs):
        Xk_tr, ytr, Xf_tr, Xk_te, yte, Xf_te = design(g, local_cells, tr, te)
        pk, _ = fit_predict(Xk_tr, ytr, Xk_te)
        pf, _ = fit_predict(Xf_tr, ytr, Xf_te)
        oof_k[te] = pk
        oof_f[te] = pf
    return y, oof_k, oof_f


def loso(df, cells):
    y_parts = []
    per = []
    for s in SAMPLES:
        te = np.where(df["sample"].to_numpy() == s)[0]
        tr = np.where(df["sample"].to_numpy() != s)[0]
        if len(te) < 8:
            continue
        Xk_tr, ytr, Xf_tr, Xk_te, yte, Xf_te = design(df, cells, tr, te)
        pk, _ = fit_predict(Xk_tr, ytr, Xk_te)
        pf, _ = fit_predict(Xf_tr, ytr, Xf_te)
        # center inside the held-out section
        y_parts.append((yte, pk, pf))
        per.append((s, r2(yte, pk), r2(yte, pf), len(te)))
    return pool(y_parts), per


def pool(parts):
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


def permutation_null(df, cells, n_perm=40):
    rng = np.random.default_rng(SEED)
    deltas = []
    # Precompute the real within-section index lists
    groups = {s: np.where(df["sample"].to_numpy() == s)[0] for s in SAMPLES}
    for p in range(n_perm):
        shuf_cells = list(cells)
        shuf_cldn = df["cldn"].to_numpy(float).copy()
        for s, idx in groups.items():
            perm = rng.permutation(len(idx))
            for a, b in zip(idx, idx[perm]):
                shuf_cells[a] = cells[b]
                shuf_cldn[a] = df["cldn"].to_numpy(float)[b]
        # temporarily use shuffled cldn by wrapping
        df_p = df.copy()
        df_p["cldn"] = shuf_cldn
        parts = []
        for s in SAMPLES:
            g = df_p[df_p["sample"] == s]
            if g["fov"].nunique() < 4:
                continue
            # cv_section uses df columns and cells; call with shuffled
            y, ok, of = cv_section(df_p, shuf_cells, s)
            m = np.isfinite(of)
            parts.append((y[m], ok[m], of[m]))
        deltas.append(pool(parts)[2])
        print(f"  perm {p+1}/{n_perm} Δ={deltas[-1]:+.4f}", flush=True)
    return np.array(deltas)


def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    print("loading", flush=True)
    xy, ct, sample, fov, logn = load()
    df, cells = make_bins(xy, ct, sample, fov, logn)
    print(f"bins {len(df)}", flush=True)
    parts = []
    per_rows = []
    for s in SAMPLES:
        y, ok, of = cv_section(df, cells, s)
        m = np.isfinite(of)
        parts.append((y[m], ok[m], of[m]))
        per_rows.append(
            {
                "sample": s,
                "patient": PATIENT[s],
                "n_bins": int(m.sum()),
                "n_fov": int(df.loc[df["sample"] == s, "fov"].nunique()),
                "r2_krt": r2(y[m], ok[m]),
                "r2_full": r2(y[m], of[m]),
                "delta_r2": r2(y[m], of[m]) - r2(y[m], ok[m]),
            }
        )
        print(f"{s}: Δ={per_rows[-1]['delta_r2']:+.3f} full={per_rows[-1]['r2_full']:.3f} n={per_rows[-1]['n_bins']}", flush=True)
    rk, rf, delta = pool(parts)
    print(f"POOLED FOV-holdout R2 krt={rk:.4f} full={rf:.4f} Δ={delta:.4f}", flush=True)
    (loso_rk, loso_rf, loso_d), loso_per = loso(df, cells)
    print(f"LOSO R2 krt={loso_rk:.4f} full={loso_rf:.4f} Δ={loso_d:.4f}", flush=True)
    print("permutation null", flush=True)
    null = permutation_null(df, cells, n_perm=30)
    # one-sided: fraction of null deltas >= observed
    p_perm = (1 + np.sum(null >= delta)) / (1 + len(null))
    per = pd.DataFrame(per_rows)
    per.to_csv(os.path.join(TAB, "neighborhood_fov_holdout_by_section.csv"), index=False)
    summary = {
        "pitch_um": PITCH,
        "min_tumor": MIN_TUMOR,
        "target": "CD8+NK fraction of cells in the bin",
        "krt_genes": KRT,
        "folds": "5-fold GroupKFold on FOV, model fit within section",
        "cldn4_terms": ["bin mean log-normalized CLDN4", "fraction of tumor cells >= training-FOV 75th percentile"],
        "n_bins": int(len(df)),
        "fov_holdout_r2_krt": rk,
        "fov_holdout_r2_full": rf,
        "fov_holdout_delta_r2": delta,
        "median_section_delta_r2": float(per["delta_r2"].median()),
        "n_sections_delta_pos": int((per["delta_r2"] > 0).sum()),
        "loso_r2_krt": loso_rk,
        "loso_r2_full": loso_rf,
        "loso_delta_r2": loso_d,
        "permutation_n": int(len(null)),
        "permutation_null_mean": float(null.mean()),
        "permutation_null_p95": float(np.quantile(null, 0.95)),
        "permutation_p_greater": float(p_perm),
        "cell_level_sample_fold_delta_r2": 0.0218,
    }
    with open(os.path.join(TAB, "neighborhood_oof_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    pd.DataFrame(
        [{"sample": s, "r2_krt": a, "r2_full": b, "delta_r2": b - a, "n": n} for s, a, b, n in loso_per]
    ).to_csv(os.path.join(TAB, "neighborhood_loso_by_section.csv"), index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    per["sample"] = pd.Categorical(per["sample"], SAMPLES, ordered=True)
    per = per.sort_values("sample")
    y = np.arange(len(per))
    colors = ["#1f4e79" if d > 0 else "#a33b32" for d in per["delta_r2"]]
    ax.axvline(0, color="#888888", lw=0.8)
    ax.barh(y, per["delta_r2"], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(per["sample"].astype(str))
    ax.set_xlabel("Held-out ΔR² (CLDN4 terms added to keratin)")
    ax.set_title("220 µm bins, FOV holdout within section")
    ax.text(
        0.0,
        -0.18,
        f"Pooled ΔR² {delta:+.3f} (keratin {rk:.3f} → keratin+CLDN4 {rf:.3f}). "
        f"Median section {per['delta_r2'].median():+.3f}. "
        f"Permutation p={p_perm:.3f}.",
        transform=ax.transAxes,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "neighborhood_delta_r2_by_section.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "neighborhood_delta_r2_by_section.pdf"))
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)
    return summary, per


if __name__ == "__main__":
    main()
