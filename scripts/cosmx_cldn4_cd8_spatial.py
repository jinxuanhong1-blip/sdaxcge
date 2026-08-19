#!/usr/bin/env python3
"""Official CosMx NSCLC FFPE (8 samples / 5 patients, 960-plex): CLDN4 vs CD8 spatial stats.

Sources
- Bruker/NanoString CosMx NSCLC FFPE showcase dataset
- Zenodo 15487520 mirror (cosmx_lung): counts, coordinates, author cell_type
- Pixel size 0.18 µm from official SMI-ReadMe (x_local_px)

Does not use any private 8-KL cohort.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression

UM_PER_PX = 0.18
ZIP_PATH = Path("/tmp/cosmx/cosmx_lung.zip")
OUT_FIG = Path("/workspace/figures")
OUT_TAB = Path("/workspace/tables")
OUT_RES = Path("/workspace/RESULTS.md")

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
PATIENT = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}
MARKER_GENES = ["CLDN4", "CD8A", "CD8B", "KRT8", "EPCAM"]
RADII_UM = np.array([15.0, 25.0, 50.0, 100.0, 150.0])
K_RADII = np.arange(10.0, 160.0, 10.0)
N_PERM = 199
RNG = np.random.default_rng(20260819)
MIN_TUMOR = 20
MIN_HIGH = 15
MIN_LOW = 15
MIN_CD8 = 8


def _log(msg: str) -> None:
    print(msg, flush=True)


def load_sample(zf: zipfile.ZipFile, sample: str) -> pd.DataFrame:
    feat = pd.read_csv(zf.open(f"{sample}/qc/features.tsv"), sep="\t", index_col=0)
    genes = list(feat.index)
    missing = [g for g in MARKER_GENES if g not in genes]
    if missing:
        raise RuntimeError(f"{sample}: missing genes {missing}")
    obs = pd.read_csv(zf.open(f"{sample}/qc/observations.tsv"), sep="\t", index_col=0)
    coord = pd.read_csv(zf.open(f"{sample}/qc/coordinates.tsv"), sep="\t", index_col=0)
    labs = pd.read_csv(zf.open(f"{sample}/labels.tsv"), sep="\t", index_col=0)
    _log(f"  reading counts.mtx for {sample} ...")
    mat = mmread(zf.open(f"{sample}/qc/counts.mtx")).tocsr()
    if mat.shape[0] != len(obs) or mat.shape[1] != len(genes):
        raise RuntimeError(f"{sample}: mtx {mat.shape} vs obs {len(obs)} genes {len(genes)}")
    idx = [genes.index(g) for g in MARKER_GENES]
    expr = pd.DataFrame(mat[:, idx].toarray(), index=obs.index, columns=MARKER_GENES)
    df = pd.concat(
        [
            expr,
            obs.rename(columns={"n_counts": "n_counts", "n_genes": "n_genes"}),
            coord.rename(columns={"x": "x_px", "y": "y_px"}),
            labs[["cell_type", "label"]].rename(columns={"label": "niche"}),
        ],
        axis=1,
    )
    df["sample"] = sample
    df["patient"] = PATIENT[sample]
    df["fov"] = df.index.astype(str).str.split("_").str[0]
    df["x_um"] = df["x_px"] * UM_PER_PX
    df["y_um"] = df["y_px"] * UM_PER_PX
    df["cell_type"] = df["cell_type"].astype(str)
    df["is_tumor"] = df["cell_type"].str.lower().str.startswith("tumor")
    df["is_epithelial"] = df["cell_type"].str.lower().eq("epithelial")
    df["is_tumor_epi"] = df["is_tumor"] | df["is_epithelial"]
    df["is_cd8"] = df["cell_type"].isin(["T CD8 naive", "T CD8 memory"])
    df["log_cldn4"] = np.log1p(df["CLDN4"].to_numpy())
    df["log_krt8"] = np.log1p(df["KRT8"].to_numpy())
    df["log_epcam"] = np.log1p(df["EPCAM"].to_numpy())
    df["log_n"] = np.log1p(df["n_counts"].to_numpy())
    return df


def residualize_cldn4(tumor: pd.DataFrame) -> np.ndarray:
    """Residual of log1p(CLDN4) after KRT8, EPCAM, and library size (not just epithelial density)."""
    x = tumor[["log_krt8", "log_epcam", "log_n"]].to_numpy()
    y = tumor["log_cldn4"].to_numpy()
    if len(tumor) < 10 or np.allclose(y, y[0]):
        return y - np.median(y)
    model = LinearRegression().fit(x, y)
    return y - model.predict(x)


def assign_high_low(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["cldn4_resid"] = np.nan
    df["cldn4_hl"] = "other"
    for (sample, fov), idx in df.groupby(["sample", "fov"]).groups.items():
        sub = df.loc[idx]
        tum = sub[sub["is_tumor"]]
        if len(tum) < MIN_TUMOR:
            continue
        resid = residualize_cldn4(tum)
        df.loc[tum.index, "cldn4_resid"] = resid
        med = np.median(resid)
        high = tum.index[resid > med]
        low = tum.index[resid <= med]
        df.loc[high, "cldn4_hl"] = "high"
        df.loc[low, "cldn4_hl"] = "low"
    return df


def window_of(xy: np.ndarray) -> tuple[float, float, float, float, float]:
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    # pad slightly so border cells are inside
    pad = 5.0
    xmin, ymin, xmax, ymax = xmin - pad, ymin - pad, xmax + pad, ymax + pad
    area = max((xmax - xmin) * (ymax - ymin), 1.0)
    return xmin, xmax, ymin, ymax, area


def nn_distances(src: np.ndarray, tgt: np.ndarray) -> np.ndarray:
    if len(src) == 0:
        return np.array([])
    if len(tgt) == 0:
        return np.full(len(src), np.nan)
    return cKDTree(tgt).query(src, k=1)[0]


def radius_counts(src: np.ndarray, tgt: np.ndarray, radii: np.ndarray) -> np.ndarray:
    """Return (n_src, n_radii) neighbor counts."""
    out = np.zeros((len(src), len(radii)), dtype=np.int32)
    if len(src) == 0 or len(tgt) == 0:
        return out
    tree = cKDTree(tgt)
    for j, r in enumerate(radii):
        out[:, j] = tree.query_ball_point(src, r, return_length=True)
    return out


def annulus_counts(src: np.ndarray, tgt: np.ndarray, r_in: float, r_out: float) -> np.ndarray:
    if len(src) == 0 or len(tgt) == 0:
        return np.zeros(len(src), dtype=np.int32)
    tree = cKDTree(tgt)
    inner = tree.query_ball_point(src, r_in, return_length=True)
    outer = tree.query_ball_point(src, r_out, return_length=True)
    return outer - inner


def bivariate_k(xy1: np.ndarray, xy2: np.ndarray, radii: np.ndarray, area: float) -> np.ndarray:
    """Border-corrected bivariate K (type1 -> type2). CSR expectation is π r^2."""
    if len(xy1) == 0 or len(xy2) == 0:
        return np.full(len(radii), np.nan)
    xmin, ymin = np.min(np.vstack([xy1, xy2]), axis=0)
    xmax, ymax = np.max(np.vstack([xy1, xy2]), axis=0)
    tree = cKDTree(xy2)
    k = np.full(len(radii), np.nan)
    n2 = len(xy2)
    for i, r in enumerate(radii):
        keep = (
            (xy1[:, 0] - xmin >= r)
            & (xmax - xy1[:, 0] >= r)
            & (xy1[:, 1] - ymin >= r)
            & (ymax - xy1[:, 1] >= r)
        )
        pts = xy1[keep]
        if len(pts) == 0:
            continue
        cnt = tree.query_ball_point(pts, r, return_length=True).astype(float)
        k[i] = area * cnt.mean() / n2
    return k


def pair_corr_from_k(k: np.ndarray, radii: np.ndarray) -> np.ndarray:
    g = np.full_like(k, np.nan, dtype=float)
    for i in range(1, len(radii) - 1):
        dr = radii[i + 1] - radii[i - 1]
        if dr <= 0 or not np.isfinite(k[i + 1]) or not np.isfinite(k[i - 1]):
            continue
        dk = k[i + 1] - k[i - 1]
        g[i] = dk / (2.0 * np.pi * radii[i] * dr)
    return g


def epithelial_lambda_grid(
    xy_all: np.ndarray, weights: np.ndarray, bw: float = 50.0, ngrid: int = 40
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (xs, ys, grid) Gaussian intensity from KRT8+EPCAM weights (µm)."""
    w = np.clip(weights.astype(float), 0, None) + 0.1
    xmin, ymin = xy_all.min(axis=0) - bw
    xmax, ymax = xy_all.max(axis=0) + bw
    xs = np.linspace(xmin, xmax, ngrid)
    ys = np.linspace(ymin, ymax, ngrid)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    grid = np.zeros((ngrid, ngrid), dtype=float)
    norm = 2.0 * np.pi * bw * bw
    inv2 = 1.0 / (2.0 * bw * bw)
    pts = xy_all
    for i in range(ngrid):
        for j in range(ngrid):
            d2 = (pts[:, 0] - xx[i, j]) ** 2 + (pts[:, 1] - yy[i, j]) ** 2
            m = d2 <= (3.0 * bw) ** 2
            if not np.any(m):
                grid[i, j] = 1e-6
            else:
                grid[i, j] = np.sum(w[m] * np.exp(-d2[m] * inv2)) / norm
    return xs, ys, np.clip(grid, 1e-6, None)


def interp_lambda(xs: np.ndarray, ys: np.ndarray, grid: np.ndarray, query: np.ndarray) -> np.ndarray:
    if len(query) == 0:
        return np.array([])
    # ys correspond to row i, xs to col j; grid[i,j] at (xs[j], ys[i])
    from scipy.interpolate import RegularGridInterpolator

    fn = RegularGridInterpolator((ys, xs), grid, bounds_error=False, fill_value=1e-6)
    return np.clip(fn(np.column_stack([query[:, 1], query[:, 0]])), 1e-6, None)


def inhomogeneous_k(
    xy1: np.ndarray,
    xy2: np.ndarray,
    lam2: np.ndarray,
    radii: np.ndarray,
    area: float,
) -> np.ndarray:
    """Bivariate inhomogeneous K using λ at type-2 (CD8) locations; CSR ~ π r^2."""
    if len(xy1) == 0 or len(xy2) == 0:
        return np.full(len(radii), np.nan)
    inv_lam = 1.0 / np.clip(lam2, 1e-6, None)
    # normalize so mean λ matches n2 / area (intensity units)
    scale = (inv_lam.mean() * len(xy2) / area)
    inv_lam = inv_lam / max(scale, 1e-12)
    tree = cKDTree(xy2)
    k = np.full(len(radii), np.nan)
    for i, r in enumerate(radii):
        neighbors = tree.query_ball_point(xy1, r)
        vals = []
        for nb in neighbors:
            if len(nb) == 0:
                vals.append(0.0)
            else:
                vals.append(float(np.sum(inv_lam[nb])))
        k[i] = area * np.mean(vals) / len(xy2)
    return k


def perm_pvalue(obs: float, null: np.ndarray, alternative: str) -> float:
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or len(null) == 0:
        return np.nan
    if alternative == "less":
        n_ext = np.sum(null <= obs)
    elif alternative == "greater":
        n_ext = np.sum(null >= obs)
    else:
        n_ext = np.sum(np.abs(null) >= abs(obs))
    return (1.0 + n_ext) / (1.0 + len(null))


def analyze_fov(sub: pd.DataFrame) -> dict | None:
    tum = sub[sub["is_tumor"]]
    high = tum[tum["cldn4_hl"] == "high"]
    low = tum[tum["cldn4_hl"] == "low"]
    cd8 = sub[sub["is_cd8"]]
    if len(high) < MIN_HIGH or len(low) < MIN_LOW or len(cd8) < MIN_CD8:
        return None
    xy_all = sub[["x_um", "y_um"]].to_numpy()
    xy_h = high[["x_um", "y_um"]].to_numpy()
    xy_l = low[["x_um", "y_um"]].to_numpy()
    xy_t = tum[["x_um", "y_um"]].to_numpy()
    xy_c = cd8[["x_um", "y_um"]].to_numpy()
    xmin, xmax, ymin, ymax, area = window_of(xy_all)

    nn_h = nn_distances(xy_h, xy_c)
    nn_l = nn_distances(xy_l, xy_c)
    rad_h = radius_counts(xy_h, xy_c, RADII_UM)
    rad_l = radius_counts(xy_l, xy_c, RADII_UM)
    contact_h = (rad_h[:, 0] > 0).mean()  # 15 µm
    contact_l = (rad_l[:, 0] > 0).mean()
    shell_h = (annulus_counts(xy_h, xy_c, 50.0, 150.0) > 0).mean()
    shell_l = (annulus_counts(xy_l, xy_c, 50.0, 150.0) > 0).mean()

    k_h = bivariate_k(xy_h, xy_c, K_RADII, area)
    k_l = bivariate_k(xy_l, xy_c, K_RADII, area)
    g_h = pair_corr_from_k(k_h, K_RADII)
    g_l = pair_corr_from_k(k_l, K_RADII)

    w_epi = (sub["KRT8"].to_numpy() + sub["EPCAM"].to_numpy()).astype(float)
    lam_xs, lam_ys, lam_grid = epithelial_lambda_grid(xy_all, w_epi, bw=50.0)
    lam_c = interp_lambda(lam_xs, lam_ys, lam_grid, xy_c)
    kinh_h = inhomogeneous_k(xy_h, xy_c, lam_c, K_RADII, area)
    kinh_l = inhomogeneous_k(xy_l, xy_c, lam_c, K_RADII, area)

    # cores
    core_all = xy_t.mean(axis=0)
    core_h = xy_h.mean(axis=0)
    core_l = xy_l.mean(axis=0)
    cd8_to_hcore = np.sqrt(((xy_c - core_h) ** 2).sum(axis=1))
    cd8_to_lcore = np.sqrt(((xy_c - core_l) ** 2).sum(axis=1))
    cd8_to_tcore = np.sqrt(((xy_c - core_all) ** 2).sum(axis=1))
    tum_r = np.sqrt(((xy_t - core_all) ** 2).sum(axis=1))
    core_cut = np.median(tum_r)
    tum_core_mask = tum_r <= core_cut
    high_in_core = tum.loc[tum.index[tum_core_mask], "cldn4_hl"].eq("high")
    # CD8 counts around tumor cells in geometric core, split by CLDN4
    rad_tum50 = radius_counts(xy_t, xy_c, np.array([50.0])).ravel()
    core_cd8_high = rad_tum50[tum_core_mask & tum["cldn4_hl"].eq("high").to_numpy()]
    core_cd8_low = rad_tum50[tum_core_mask & tum["cldn4_hl"].eq("low").to_numpy()]

    # CD8-label permutation among non-tumor cells (keeps tumor geometry)
    nontum = sub.loc[~sub["is_tumor"], ["x_um", "y_um"]].to_numpy()
    n_cd8 = len(xy_c)
    if len(nontum) < n_cd8 + 5:
        nontum = xy_all
    null_nn_h = np.zeros(N_PERM)
    null_cnt25 = np.zeros(N_PERM)
    null_cnt50 = np.zeros(N_PERM)
    null_cnt100 = np.zeros(N_PERM)
    null_k = np.zeros((N_PERM, len(K_RADII)))
    null_kinh = np.zeros((N_PERM, len(K_RADII)))
    for p in range(N_PERM):
        pick = RNG.choice(len(nontum), size=n_cd8, replace=False)
        xy_p = nontum[pick]
        null_nn_h[p] = np.nanmedian(nn_distances(xy_h, xy_p))
        rc = radius_counts(xy_h, xy_p, np.array([25.0, 50.0, 100.0]))
        null_cnt25[p] = rc[:, 0].mean()
        null_cnt50[p] = rc[:, 1].mean()
        null_cnt100[p] = rc[:, 2].mean()
        null_k[p] = bivariate_k(xy_h, xy_p, K_RADII, area)
        # intensity-preserving: reuse epithelial λ evaluated at permuted CD8 sites
        lam_p = interp_lambda(lam_xs, lam_ys, lam_grid, xy_p)
        null_kinh[p] = inhomogeneous_k(xy_h, xy_p, lam_p, K_RADII, area)

    # high/low label permutation among tumor cells
    resid = tum["cldn4_resid"].to_numpy()
    n_high = len(high)
    null_d_nn = np.zeros(N_PERM)
    null_d_c50 = np.zeros(N_PERM)
    for p in range(N_PERM):
        pick = RNG.choice(len(tum), size=n_high, replace=False)
        mask = np.zeros(len(tum), dtype=bool)
        mask[pick] = True
        d_h = np.nanmedian(nn_distances(xy_t[mask], xy_c))
        d_l = np.nanmedian(nn_distances(xy_t[~mask], xy_c))
        null_d_nn[p] = d_h - d_l
        c_h = radius_counts(xy_t[mask], xy_c, np.array([50.0])).mean()
        c_l = radius_counts(xy_t[~mask], xy_c, np.array([50.0])).mean()
        null_d_c50[p] = c_h - c_l

    obs_d_nn = float(np.nanmedian(nn_h) - np.nanmedian(nn_l))
    obs_d_c50 = float(rad_h[:, 2].mean() - rad_l[:, 2].mean()) if False else float(
        rad_h[:, np.where(RADII_UM == 50.0)[0][0]].mean()
        - rad_l[:, np.where(RADII_UM == 50.0)[0][0]].mean()
    )
    i25 = int(np.where(RADII_UM == 25.0)[0][0])
    i50 = int(np.where(RADII_UM == 50.0)[0][0])
    i100 = int(np.where(RADII_UM == 100.0)[0][0])
    ik25 = int(np.where(K_RADII == 25.0)[0][0]) if 25.0 in K_RADII else 1
    ik50 = int(np.where(K_RADII == 50.0)[0][0])
    ik100 = int(np.where(K_RADII == 100.0)[0][0])

    rec = {
        "sample": sub["sample"].iloc[0],
        "patient": sub["patient"].iloc[0],
        "fov": str(sub["fov"].iloc[0]),
        "n_cells": int(len(sub)),
        "n_tumor": int(len(tum)),
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "n_cd8": int(len(cd8)),
        "n_epithelial": int(sub["is_epithelial"].sum()),
        "area_um2": float(area),
        "mean_cldn4_tumor": float(tum["CLDN4"].mean()),
        "mean_cldn4_other": float(sub.loc[~sub["is_tumor"], "CLDN4"].mean()),
        "median_nn_high": float(np.nanmedian(nn_h)),
        "median_nn_low": float(np.nanmedian(nn_l)),
        "mean_nn_high": float(np.nanmean(nn_h)),
        "mean_nn_low": float(np.nanmean(nn_l)),
        "delta_median_nn": obs_d_nn,
        "cnt25_high": float(rad_h[:, i25].mean()),
        "cnt25_low": float(rad_l[:, i25].mean()),
        "cnt50_high": float(rad_h[:, i50].mean()),
        "cnt50_low": float(rad_l[:, i50].mean()),
        "cnt100_high": float(rad_h[:, i100].mean()),
        "cnt100_low": float(rad_l[:, i100].mean()),
        "contact15_high": float(contact_h),
        "contact15_low": float(contact_l),
        "shell50_150_high": float(shell_h),
        "shell50_150_low": float(shell_l),
        "k25_high": float(k_h[ik25]),
        "k50_high": float(k_h[ik50]),
        "k100_high": float(k_h[ik100]),
        "k25_low": float(k_l[ik25]),
        "k50_low": float(k_l[ik50]),
        "k100_low": float(k_l[ik100]),
        "kinh50_high": float(kinh_h[ik50]),
        "kinh50_low": float(kinh_l[ik50]),
        "cd8_med_dist_highcore": float(np.median(cd8_to_hcore)),
        "cd8_med_dist_lowcore": float(np.median(cd8_to_lcore)),
        "cd8_p10_dist_highcore": float(np.percentile(cd8_to_hcore, 10)),
        "cd8_p10_dist_lowcore": float(np.percentile(cd8_to_lcore, 10)),
        "cd8_med_dist_tumorcore": float(np.median(cd8_to_tcore)),
        "core50_cd8_high": float(np.mean(core_cd8_high)) if len(core_cd8_high) else np.nan,
        "core50_cd8_low": float(np.mean(core_cd8_low)) if len(core_cd8_low) else np.nan,
        "p_perm_cd8_nn_high_greater": perm_pvalue(float(np.nanmedian(nn_h)), null_nn_h, "greater"),
        "p_perm_cd8_cnt25_less": perm_pvalue(float(rad_h[:, i25].mean()), null_cnt25, "less"),
        "p_perm_cd8_cnt50_less": perm_pvalue(float(rad_h[:, i50].mean()), null_cnt50, "less"),
        "p_perm_cd8_cnt100_less": perm_pvalue(float(rad_h[:, i100].mean()), null_cnt100, "less"),
        "p_perm_cd8_k50_less": perm_pvalue(float(k_h[ik50]), null_k[:, ik50], "less"),
        "p_perm_cd8_kinh50_less": perm_pvalue(float(kinh_h[ik50]), null_kinh[:, ik50], "less"),
        "p_perm_hl_dnn_greater": perm_pvalue(obs_d_nn, null_d_nn, "greater"),
        "p_perm_hl_dc50_less": perm_pvalue(obs_d_c50, null_d_c50, "less"),
        "k_high": k_h,
        "k_low": k_l,
        "g_high": g_h,
        "g_low": g_l,
        "kinh_high": kinh_h,
        "kinh_low": kinh_l,
        "null_k_q05": np.nanpercentile(null_k, 5, axis=0),
        "null_k_q50": np.nanpercentile(null_k, 50, axis=0),
        "null_k_q95": np.nanpercentile(null_k, 95, axis=0),
        "null_kinh_q05": np.nanpercentile(null_kinh, 5, axis=0),
        "null_kinh_q50": np.nanpercentile(null_kinh, 50, axis=0),
        "null_kinh_q95": np.nanpercentile(null_kinh, 95, axis=0),
        "nn_high": nn_h,
        "nn_low": nn_l,
        "xy_h": xy_h,
        "xy_l": xy_l,
        "xy_c": xy_c,
        "xy_other": sub.loc[~sub["is_tumor"] & ~sub["is_cd8"], ["x_um", "y_um"]].to_numpy(),
    }
    return rec


def fisher_combined(pvals: np.ndarray) -> float:
    p = pvals[np.isfinite(pvals) & (pvals > 0)]
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 1e-300, 1.0)
    stat = -2.0 * np.sum(np.log(p))
    return float(stats.chi2.sf(stat, df=2 * len(p)))


def mean_se(x: np.ndarray) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan
    return float(np.mean(x)), float(stats.sem(x)) if len(x) > 1 else 0.0


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_maps(recs: list[dict], path: Path) -> None:
    recs_sorted = sorted(recs, key=lambda r: r["n_cd8"] * r["n_high"], reverse=True)
    pick = recs_sorted[:4]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.2))
    for ax, r in zip(axes.ravel(), pick):
        o, h, l, c = r["xy_other"], r["xy_h"], r["xy_l"], r["xy_c"]
        if len(o):
            ax.scatter(o[:, 0], o[:, 1], s=2, c="#d9d9d9", linewidths=0, rasterized=True)
        ax.scatter(l[:, 0], l[:, 1], s=6, c="#f4a261", linewidths=0, label="CLDN4-low tumor", rasterized=True)
        ax.scatter(h[:, 0], h[:, 1], s=6, c="#9b2226", linewidths=0, label="CLDN4-high tumor", rasterized=True)
        ax.scatter(c[:, 0], c[:, 1], s=10, c="#1d4e89", linewidths=0, label="CD8 T", rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"{r['sample']} FOV {r['fov']}  (tumor {r['n_tumor']}, CD8 {r['n_cd8']})", fontsize=9)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        style_ax(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("CosMx NSCLC: CLDN4-high/low tumor cells and CD8 T cells", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_box_cldn4(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    tumor = np.log1p(df.loc[df["is_tumor"], "CLDN4"].to_numpy())
    other = np.log1p(df.loc[~df["is_tumor"], "CLDN4"].to_numpy())
    parts = ax.boxplot(
        [other, tumor],
        tick_labels=["other types", "tumor (author)"],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(parts["boxes"], ["#9aa0a6", "#9b2226"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    u, p = stats.mannwhitneyu(tumor, other, alternative="greater")
    ax.set_ylabel("log1p CLDN4 counts")
    ax.set_title(
        f"CLDN4 higher in tumor vs other\n"
        f"n_tumor={len(tumor):,}  n_other={len(other):,}  MWU p={p:.2e}"
    )
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_paired_nn(fov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
    ax = axes[0]
    x = np.array([0, 1])
    for _, row in fov.iterrows():
        ax.plot(x, [row["median_nn_low"], row["median_nn_high"]], color="#bbbbbb", lw=0.7, alpha=0.7)
    ax.scatter(np.zeros(len(fov)), fov["median_nn_low"], s=18, c="#f4a261", zorder=3, label="CLDN4-low")
    ax.scatter(np.ones(len(fov)), fov["median_nn_high"], s=18, c="#9b2226", zorder=3, label="CLDN4-high")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-low\ntumor", "CLDN4-high\ntumor"])
    ax.set_ylabel("Median distance to nearest CD8 (µm)")
    stat, p = stats.wilcoxon(fov["median_nn_high"], fov["median_nn_low"], alternative="greater")
    d = (fov["median_nn_high"] - fov["median_nn_low"]).median()
    ax.set_title(f"Per-FOV paired medians\nn={len(fov)} FOVs  ∆med={d:.1f} µm  Wilcoxon p={p:.2e}")
    style_ax(ax)
    ax = axes[1]
    parts = ax.boxplot(
        [fov["cnt50_low"], fov["cnt50_high"]],
        tick_labels=["CLDN4-low", "CLDN4-high"],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(parts["boxes"], ["#f4a261", "#9b2226"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    stat, p = stats.wilcoxon(fov["cnt50_high"], fov["cnt50_low"], alternative="less")
    d = (fov["cnt50_high"] - fov["cnt50_low"]).median()
    ax.set_ylabel("Mean CD8 count within 50 µm")
    ax.set_title(f"50 µm neighbors\n∆med={d:.3f}  Wilcoxon p={p:.2e}")
    style_ax(ax)
    fig.suptitle("CLDN4-high vs low tumor cells (FOV-level)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_ecdf(recs: list[dict], path: Path) -> None:
    h = np.concatenate([r["nn_high"] for r in recs])
    l = np.concatenate([r["nn_low"] for r in recs])
    h = h[np.isfinite(h)]
    l = l[np.isfinite(l)]

    def ecdf(v):
        v = np.sort(v)
        y = np.arange(1, len(v) + 1) / len(v)
        return v, y

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    xh, yh = ecdf(h)
    xl, yl = ecdf(l)
    ax.plot(xl, yl, color="#f4a261", lw=2, label=f"CLDN4-low (n={len(l):,})")
    ax.plot(xh, yh, color="#9b2226", lw=2, label=f"CLDN4-high (n={len(h):,})")
    ax.set_xlim(0, 200)
    ax.set_xlabel("Distance to nearest CD8 T cell (µm)")
    ax.set_ylabel("ECDF")
    p = stats.ks_2samp(h, l, alternative="less").pvalue
    ax.set_title(f"Nearest-CD8 distance ECDF (pooled tumor cells)\nKS one-sided p={p:.2e}")
    ax.legend(frameon=False)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_neighbor_bars(fov: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    xs = np.arange(3)
    w = 0.35
    means_l = [fov["cnt25_low"].mean(), fov["cnt50_low"].mean(), fov["cnt100_low"].mean()]
    means_h = [fov["cnt25_high"].mean(), fov["cnt50_high"].mean(), fov["cnt100_high"].mean()]
    se_l = [stats.sem(fov["cnt25_low"]), stats.sem(fov["cnt50_low"]), stats.sem(fov["cnt100_low"])]
    se_h = [stats.sem(fov["cnt25_high"]), stats.sem(fov["cnt50_high"]), stats.sem(fov["cnt100_high"])]
    ax.bar(xs - w / 2, means_l, w, yerr=se_l, color="#f4a261", label="CLDN4-low", capsize=3)
    ax.bar(xs + w / 2, means_h, w, yerr=se_h, color="#9b2226", label="CLDN4-high", capsize=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(["25 µm", "50 µm", "100 µm"])
    ax.set_ylabel("Mean CD8 neighbors per tumor cell (FOV mean ± SEM)")
    ax.set_title(f"CD8 radial counts  (n={len(fov)} FOVs)")
    ax.legend(frameon=False)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def stack_curves(recs: list[dict], key: str) -> np.ndarray:
    return np.vstack([r[key] for r in recs])


def plot_k_envelope(recs: list[dict], path: Path, kind: str) -> None:
    if kind == "K":
        obs_h = stack_curves(recs, "k_high")
        obs_l = stack_curves(recs, "k_low")
        lo = stack_curves(recs, "null_k_q05")
        mid = stack_curves(recs, "null_k_q50")
        hi = stack_curves(recs, "null_k_q95")
        ylab = r"Bivariate $\hat{K}_{12}(r)$ (CLDN4-high → CD8)"
        csr = np.pi * K_RADII**2
        title = "Bivariate Ripley's K vs CD8-label permutation envelope"
    elif kind == "g":
        obs_h = stack_curves(recs, "g_high")
        obs_l = stack_curves(recs, "g_low")
        # envelope from K perms via finite difference on FOV-mean null K
        k_lo = np.nanmean(stack_curves(recs, "null_k_q05"), axis=0)
        k_mid = np.nanmean(stack_curves(recs, "null_k_q50"), axis=0)
        k_hi = np.nanmean(stack_curves(recs, "null_k_q95"), axis=0)
        lo = pair_corr_from_k(k_lo, K_RADII)[None, :]
        mid = pair_corr_from_k(k_mid, K_RADII)[None, :]
        hi = pair_corr_from_k(k_hi, K_RADII)[None, :]
        ylab = r"Pair correlation $\hat{g}_{12}(r)$"
        csr = np.ones_like(K_RADII)
        title = "Pair-correlation g(r) vs permutation envelope"
    else:
        obs_h = stack_curves(recs, "kinh_high")
        obs_l = stack_curves(recs, "kinh_low")
        lo = stack_curves(recs, "null_kinh_q05")
        mid = stack_curves(recs, "null_kinh_q50")
        hi = stack_curves(recs, "null_kinh_q95")
        ylab = r"Inhomogeneous $\hat{K}_{inhom}(r)$ (λ = KRT8+EPCAM)"
        csr = np.pi * K_RADII**2
        title = "Inhomogeneous K controlling epithelial density"
    mh = np.nanmean(obs_h, axis=0)
    ml = np.nanmean(obs_l, axis=0)
    elo = np.nanmean(lo, axis=0)
    emid = np.nanmean(mid, axis=0)
    ehi = np.nanmean(hi, axis=0)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.fill_between(K_RADII, elo, ehi, color="#1d4e89", alpha=0.15, label="CD8-label perm 5–95%")
    ax.plot(K_RADII, emid, color="#1d4e89", ls="--", lw=1, label="perm median")
    ax.plot(K_RADII, csr, color="0.4", ls=":", lw=1.2, label="CSR / independence")
    ax.plot(K_RADII, ml, color="#f4a261", lw=2, label="CLDN4-low tumor → CD8")
    ax.plot(K_RADII, mh, color="#9b2226", lw=2, label="CLDN4-high tumor → CD8")
    ax.set_xlabel("r (µm)")
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nmean of {len(recs)} QC FOVs")
    ax.legend(frameon=False, fontsize=8)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_core_margin(fov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
    ax = axes[0]
    x = np.array([0, 1])
    for _, row in fov.iterrows():
        ax.plot(
            x,
            [row["cd8_med_dist_lowcore"], row["cd8_med_dist_highcore"]],
            color="#bbbbbb",
            lw=0.7,
            alpha=0.7,
        )
    ax.scatter(np.zeros(len(fov)), fov["cd8_med_dist_lowcore"], s=18, c="#f4a261", zorder=3)
    ax.scatter(np.ones(len(fov)), fov["cd8_med_dist_highcore"], s=18, c="#9b2226", zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["to CLDN4-low\ncore", "to CLDN4-high\ncore"])
    ax.set_ylabel("Median CD8 distance to core (µm)")
    stat, p = stats.wilcoxon(
        fov["cd8_med_dist_highcore"], fov["cd8_med_dist_lowcore"], alternative="greater"
    )
    d = (fov["cd8_med_dist_highcore"] - fov["cd8_med_dist_lowcore"]).median()
    ax.set_title(f"CD8 vs cores\n∆med={d:.1f} µm  Wilcoxon p={p:.2e}")
    style_ax(ax)
    ax = axes[1]
    a = fov["core50_cd8_low"].to_numpy()
    b = fov["core50_cd8_high"].to_numpy()
    m = np.isfinite(a) & np.isfinite(b)
    parts = ax.boxplot(
        [a[m], b[m]],
        tick_labels=["low-CLDN4\nin tumor core", "high-CLDN4\nin tumor core"],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(parts["boxes"], ["#f4a261", "#9b2226"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    if m.sum() >= 8:
        stat, p = stats.wilcoxon(b[m], a[m], alternative="less")
    else:
        p = np.nan
    ax.set_ylabel("Mean CD8 within 50 µm of core tumor cells")
    ax.set_title(f"Geometric tumor core\nWilcoxon p={p:.2e}")
    style_ax(ax)
    fig.suptitle("Tumor-core vs margin / CLDN4-high vs low cores", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_contact(fov: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    xs = np.arange(2)
    w = 0.35
    means_l = [fov["contact15_low"].mean(), fov["shell50_150_low"].mean()]
    means_h = [fov["contact15_high"].mean(), fov["shell50_150_high"].mean()]
    se_l = [stats.sem(fov["contact15_low"]), stats.sem(fov["shell50_150_low"])]
    se_h = [stats.sem(fov["contact15_high"]), stats.sem(fov["shell50_150_high"])]
    ax.bar(xs - w / 2, means_l, w, yerr=se_l, color="#f4a261", label="CLDN4-low", capsize=3)
    ax.bar(xs + w / 2, means_h, w, yerr=se_h, color="#9b2226", label="CLDN4-high", capsize=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(["Contact\n(CD8 ≤ 15 µm)", "Long-range shell\n(CD8 in 50–150 µm)"])
    ax.set_ylabel("Fraction of tumor cells")
    p15 = stats.wilcoxon(fov["contact15_high"], fov["contact15_low"], alternative="less").pvalue
    psh = stats.wilcoxon(fov["shell50_150_high"], fov["shell50_150_low"], alternative="less").pvalue
    ax.set_title(f"Contact vs long-range  (n={len(fov)} FOVs)\n15 µm p={p15:.2e}   50–150 µm p={psh:.2e}")
    ax.legend(frameon=False)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_results(df: pd.DataFrame, fov: pd.DataFrame, recs: list[dict], panel_ok: dict) -> None:
    n_cells = len(df)
    n_fov_all = df.groupby(["sample", "fov"]).ngroups
    n_tumor = int(df["is_tumor"].sum())
    n_cd8 = int(df["is_cd8"].sum())
    n_epi = int(df["is_epithelial"].sum())
    n_high = int((df["cldn4_hl"] == "high").sum())
    n_low = int((df["cldn4_hl"] == "low").sum())
    n_fov = len(fov)

    tum = np.log1p(df.loc[df["is_tumor"], "CLDN4"].to_numpy())
    oth = np.log1p(df.loc[~df["is_tumor"], "CLDN4"].to_numpy())
    p_expr = stats.mannwhitneyu(tum, oth, alternative="greater").pvalue
    d_nn = (fov["median_nn_high"] - fov["median_nn_low"]).median()
    p_nn = stats.wilcoxon(fov["median_nn_high"], fov["median_nn_low"], alternative="greater").pvalue
    d_c50 = (fov["cnt50_high"] - fov["cnt50_low"]).median()
    p_c50 = stats.wilcoxon(fov["cnt50_high"], fov["cnt50_low"], alternative="less").pvalue
    d_c25 = (fov["cnt25_high"] - fov["cnt25_low"]).median()
    p_c25 = stats.wilcoxon(fov["cnt25_high"], fov["cnt25_low"], alternative="less").pvalue
    d_c100 = (fov["cnt100_high"] - fov["cnt100_low"]).median()
    p_c100 = stats.wilcoxon(fov["cnt100_high"], fov["cnt100_low"], alternative="less").pvalue
    p_fish_nn = fisher_combined(fov["p_perm_cd8_nn_high_greater"].to_numpy())
    p_fish_50 = fisher_combined(fov["p_perm_cd8_cnt50_less"].to_numpy())
    p_fish_k50 = fisher_combined(fov["p_perm_cd8_k50_less"].to_numpy())
    p_fish_kinh = fisher_combined(fov["p_perm_cd8_kinh50_less"].to_numpy())
    p_core = stats.wilcoxon(
        fov["cd8_med_dist_highcore"], fov["cd8_med_dist_lowcore"], alternative="greater"
    ).pvalue
    d_core = (fov["cd8_med_dist_highcore"] - fov["cd8_med_dist_lowcore"]).median()
    p15 = stats.wilcoxon(fov["contact15_high"], fov["contact15_low"], alternative="less").pvalue
    d15 = (fov["contact15_high"] - fov["contact15_low"]).median()
    psh = stats.wilcoxon(fov["shell50_150_high"], fov["shell50_150_low"], alternative="less").pvalue
    dsh = (fov["shell50_150_high"] - fov["shell50_150_low"]).median()

    # sample-level means (5 patients / 8 samples)
    samp = fov.groupby("sample").median(numeric_only=True)
    p_nn_s = (
        stats.wilcoxon(samp["median_nn_high"], samp["median_nn_low"], alternative="greater").pvalue
        if len(samp) >= 6
        else np.nan
    )

    paper = (
        f"In the official CosMx NSCLC FFPE 960-plex dataset (8 samples / 5 patients; "
        f"{n_tumor:,} author-typed tumor cells and {n_cd8:,} CD8 T cells across {n_fov} QC FOVs "
        f"of {n_fov_all} total FOVs; {n_cells:,} cells), CLDN4-high tumor cells "
        f"(residualized on KRT8/EPCAM and library size) were farther from the nearest CD8 T cell "
        f"than CLDN4-low tumor cells (FOV-paired median ∆distance = {d_nn:.1f} µm; "
        f"Wilcoxon signed-rank p = {fmt_p(p_nn)}; n = {n_fov} FOVs) and had fewer CD8 neighbors "
        f"within 50 µm (∆ = {d_c50:.3f} cells; p = {fmt_p(p_c50)}). "
        f"CD8-label permutation (n = {N_PERM} per FOV) supported local exclusion at 50 µm "
        f"(Fisher combined p = {fmt_p(p_fish_50)})."
    )

    lines = [
        "# RESULTS: CosMx NSCLC FFPE — do CLDN4-high tumor cells spatially exclude CD8 T cells?",
        "",
        "## Dataset and panel check",
        "",
        "- **Dataset:** official NanoString/Bruker CosMx SMI NSCLC FFPE showcase (8 FFPE sections / 5 patients, 960-plex prototype). Product page: https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/",
        "- **Counts / coordinates / author cell types used here:** Zenodo 15487520 `cosmx_lung` mirror of that public release (He et al., *Nat Biotechnol* 2022).",
        f"- **Panel check:** CLDN4 is present on the 960-gene expression matrix (also CD8A, CD8B, KRT8, EPCAM). Official SMI-ReadMe lists CLDN4 among exprMat targets.",
        f"- **Pixel size:** 0.18 µm/pixel (official SMI-ReadMe).",
        "- **Author cell types used:** tumor 5/6/9/12/13 as tumor; T CD8 naive + T CD8 memory as CD8. No private 8-KL data were used.",
        f"- **Loaded:** {n_cells:,} cells, {n_fov_all} FOVs, 8 samples, 5 patients.",
        f"- **Tumor / CD8 / epithelial:** n_tumor = {n_tumor:,}; n_CD8 = {n_cd8:,}; n_epithelial = {n_epi:,}.",
        f"- **CLDN4-high / low tumor (per-FOV residual median split):** n_high = {n_high:,}; n_low = {n_low:,}.",
        f"- **QC FOVs for spatial tests:** {n_fov} FOVs with ≥{MIN_HIGH} high, ≥{MIN_LOW} low, ≥{MIN_CD8} CD8.",
        "",
        "## Paper sentence",
        "",
        paper,
        "",
        "## 1. CLDN4 expression in tumor vs other types",
        "",
        f"- Mean raw CLDN4 counts: tumor {df.loc[df['is_tumor'],'CLDN4'].mean():.3f} vs other {df.loc[~df['is_tumor'],'CLDN4'].mean():.3f}.",
        f"- log1p CLDN4 Mann–Whitney (tumor > other): p = {fmt_p(p_expr)} (n_tumor = {n_tumor:,}, n_other = {n_cells-n_tumor:,}).",
        "- Per-sample means are in `tables/cldn4_by_sample.csv`.",
        "",
        "## 2–3. Nearest-CD8 distance and radial CD8 counts (CLDN4-high vs low tumor)",
        "",
        "CLDN4-high/low is a **per-FOV median split of residual log1p(CLDN4)** after linear adjustment for log1p(KRT8), log1p(EPCAM), and log1p(n_counts), so the contrast is not simply epithelial density.",
        "",
        f"- FOV-paired median nearest-CD8 distance: high {fov['median_nn_high'].median():.2f} µm vs low {fov['median_nn_low'].median():.2f} µm; ∆ = {d_nn:.2f} µm; Wilcoxon p = {fmt_p(p_nn)}; n = {n_fov} FOVs.",
        f"- Mean nearest-CD8 (FOV means): high {fov['mean_nn_high'].mean():.2f} vs low {fov['mean_nn_low'].mean():.2f} µm.",
        f"- CD8 count within 25 µm: ∆med (high−low) = {d_c25:.3f}; p = {fmt_p(p_c25)}.",
        f"- CD8 count within 50 µm: ∆med (high−low) = {d_c50:.3f}; p = {fmt_p(p_c50)}.",
        f"- CD8 count within 100 µm: ∆med (high−low) = {d_c100:.3f}; p = {fmt_p(p_c100)}.",
        f"- Sample-level paired nearest-CD8 (n = {len(samp)} samples): Wilcoxon p = {fmt_p(p_nn_s)}.",
        "",
        "## 4. Mixing / segregation (Ripley's K and g(r))",
        "",
        "Bivariate Ripley's K̂₁₂(r) from CLDN4-high (or low) tumor cells to CD8, border-corrected, FOV-averaged. Independence / CSR reference is πr². The permutation envelope shuffles CD8 labels among non-tumor cells (n = 199), keeping tumor positions fixed.",
        "",
        f"- At r = 50 µm, mean K̂_high = {fov['k50_high'].mean():.1f} vs K̂_low = {fov['k50_low'].mean():.1f} vs CSR πr² = {np.pi*50**2:.1f}.",
        f"- Fisher combination of per-FOV one-sided CD8-label permutation p for K(50) (exclusion): p = {fmt_p(p_fish_k50)}.",
        f"- Mean fraction of per-FOV permutation tests with p < 0.05 at 50 µm neighbor count: {(fov['p_perm_cd8_cnt50_less']<0.05).mean():.3f}.",
        "",
        "## 5. Residualization / epithelial-density control",
        "",
        "- High/low split uses CLDN4 residuals on KRT8 + EPCAM + library size (above).",
        "- Inhomogeneous K uses a 50 µm Gaussian intensity surface λ from KRT8+EPCAM weights on all cells in the FOV.",
        f"- Inhomogeneous K̂(50) high vs low: {fov['kinh50_high'].mean():.1f} vs {fov['kinh50_low'].mean():.1f} (CSR πr² = {np.pi*50**2:.1f}).",
        f"- Fisher combined permutation p for inhomogeneous K(50) exclusion: p = {fmt_p(p_fish_kinh)}.",
        "",
        "## High-end tests requested",
        "",
        "### Label-permutation null (CD8 labels; n = 199 / FOV)",
        "",
        f"- Combined p, nearest-CD8 of CLDN4-high larger than CD8-shuffled null: {fmt_p(p_fish_nn)}.",
        f"- Combined p, fewer CD8 within 25 / 50 / 100 µm than null: "
        f"{fmt_p(fisher_combined(fov['p_perm_cd8_cnt25_less'].to_numpy()))} / "
        f"{fmt_p(p_fish_50)} / "
        f"{fmt_p(fisher_combined(fov['p_perm_cd8_cnt100_less'].to_numpy()))}.",
        f"- Median per-FOV permutation p at 50 µm neighbor count: {fov['p_perm_cd8_cnt50_less'].median():.3f}.",
        f"- High-vs-low tumor-label permutation for ∆NN (high farther) Fisher p = {fmt_p(fisher_combined(fov['p_perm_hl_dnn_greater'].to_numpy()))}; "
        f"for ∆CD8@50 µm Fisher p = {fmt_p(fisher_combined(fov['p_perm_hl_dc50_less'].to_numpy()))}.",
        "",
        "### Tumor-core vs margin",
        "",
        f"- Median CD8 distance to CLDN4-high core vs CLDN4-low core: {fov['cd8_med_dist_highcore'].median():.1f} vs {fov['cd8_med_dist_lowcore'].median():.1f} µm; ∆ = {d_core:.1f} µm; Wilcoxon p = {fmt_p(p_core)} (n = {n_fov} FOVs).",
        f"- 10th percentile CD8–core distance (infiltration depth): high-core {fov['cd8_p10_dist_highcore'].median():.1f} µm vs low-core {fov['cd8_p10_dist_lowcore'].median():.1f} µm.",
        f"- CD8 within 50 µm of geometric tumor-core cells: high {fov['core50_cd8_high'].median():.3f} vs low {fov['core50_cd8_low'].median():.3f}.",
        "",
        "### Contact vs long-range shell",
        "",
        f"- Fraction of tumor cells with a CD8 within 15 µm (contact): high {fov['contact15_high'].mean():.3f} vs low {fov['contact15_low'].mean():.3f}; ∆med = {d15:.3f}; Wilcoxon p = {fmt_p(p15)}.",
        f"- Fraction with a CD8 in the 50–150 µm shell: high {fov['shell50_150_high'].mean():.3f} vs low {fov['shell50_150_low'].mean():.3f}; ∆med = {dsh:.3f}; Wilcoxon p = {fmt_p(psh)}.",
        "",
        "## Counts",
        "",
        f"| Level | n |",
        f"|---|---|",
        f"| Patients | 5 |",
        f"| Samples | 8 |",
        f"| Cells loaded | {n_cells:,} |",
        f"| FOVs loaded | {n_fov_all} |",
        f"| QC FOVs (spatial tests) | {n_fov} |",
        f"| Author tumor cells | {n_tumor:,} |",
        f"| CLDN4-high / low tumor | {n_high:,} / {n_low:,} |",
        f"| Author CD8 T cells | {n_cd8:,} |",
        f"| Author epithelial cells | {n_epi:,} |",
        "",
        "## Figures",
        "",
        "- `figures/fig01_maps.png` — example FOV maps (CLDN4-high/low tumor, CD8).",
        "- `figures/fig02_cldn4_tumor_vs_other.png` — CLDN4 in tumor vs other types.",
        "- `figures/fig03_nn_distance_paired.png` — paired FOV nearest-CD8 distance and 50 µm counts.",
        "- `figures/fig04_nn_distance_ecdf.png` — pooled nearest-CD8 ECDFs.",
        "- `figures/fig05_neighbor_counts.png` — 25/50/100 µm CD8 counts.",
        "- `figures/fig06_K_hat_envelope.png` — bivariate K̂(r) vs CSR and permutation envelope.",
        "- `figures/fig07_g_r_envelope.png` — pair-correlation ĝ(r).",
        "- `figures/fig08_K_inhom_envelope.png` — inhomogeneous K with KRT8/EPCAM λ.",
        "- `figures/fig09_core_margin.png` — CD8 distance to high vs low cores; core infiltration.",
        "- `figures/fig10_contact_vs_shell.png` — 15 µm contact vs 50–150 µm shell.",
        "",
        "## Methods (short)",
        "",
        "Cells were taken from the public CosMx NSCLC 960-plex matrices. Coordinates were converted to µm with the official 0.18 µm/pixel scale. Tumor and CD8 labels are the author/NanoString cell types distributed with the Zenodo mirror (`cell_type`). CLDN4-high vs low was defined **within each FOV** as a median split of residual log1p(CLDN4) after OLS on log1p(KRT8), log1p(EPCAM), and log1p(total counts). Nearest-CD8 distances and radius counts used Euclidean KD-trees. Ripley's bivariate K used border (minus-sampling) correction; g(r) is the finite-difference derivative of K. Inhomogeneous K used a 50 µm Gaussian kernel intensity of KRT8+EPCAM. CD8-label permutations redrew the same number of CD8 positions from non-tumor cells (n = 199). High-vs-low permutations reassigned the high label among tumor cells. Primary inference is FOV-paired Wilcoxon signed-rank; permutation p-values were combined across FOVs with Fisher's method. FOVs lacking the minimum tumor/CD8 counts were excluded from spatial tests but included in expression summaries.",
        "",
        "## Files",
        "",
        "- `tables/fov_spatial_stats.csv` — one row per QC FOV.",
        "- `tables/cldn4_by_sample.csv` — expression and cell counts per sample.",
        "- `tables/summary_key_numbers.json` — machine-readable key numbers.",
        "",
        "## Notes",
        "",
        "- This analysis is additive public spatial evidence for CLDN4 only. It does not use private 8-KL material.",
        "- Author CD8 types were used (not a CD8A+ threshold). CD8A is on the panel and was loaded for QC.",
        f"- Panel genes confirmed present: {', '.join(panel_ok['genes_present'])}.",
    ]
    OUT_RES.write_text("\n".join(lines) + "\n")
    key = {
        "n_cells": n_cells,
        "n_fov_all": int(n_fov_all),
        "n_fov_qc": n_fov,
        "n_tumor": n_tumor,
        "n_cd8": n_cd8,
        "n_high": n_high,
        "n_low": n_low,
        "delta_median_nn_um": float(d_nn),
        "p_wilcoxon_nn": float(p_nn),
        "delta_cd8_50um": float(d_c50),
        "p_wilcoxon_cd8_50": float(p_c50),
        "p_fisher_perm_cd8_50": float(p_fish_50),
        "paper_sentence": paper,
    }
    (OUT_TAB / "summary_key_numbers.json").write_text(json.dumps(key, indent=2))


def main() -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    if not ZIP_PATH.exists():
        raise SystemExit(f"Missing {ZIP_PATH}")

    _log("Loading official CosMx NSCLC (Zenodo cosmx_lung mirror) ...")
    frames = []
    with zipfile.ZipFile(ZIP_PATH) as zf:
        # panel check on first sample before proceeding
        feat0 = pd.read_csv(zf.open(f"{SAMPLES[0]}/qc/features.tsv"), sep="\t", index_col=0)
        genes = list(feat0.index)
        present = [g for g in MARKER_GENES if g in genes]
        _log(f"Panel genes present: {present}  (n_panel={len(genes)})")
        if "CLDN4" not in genes:
            Path("/workspace/RESULTS.md").write_text(
                "# RESULTS\n\nCLDN4 is absent from the CosMx 960-plex NSCLC expression matrix. Analysis stopped.\n"
            )
            raise SystemExit("CLDN4 missing from panel")
        for sample in SAMPLES:
            _log(f"load {sample}")
            frames.append(load_sample(zf, sample))
    df = pd.concat(frames, axis=0)
    _log(f"Loaded {len(df):,} cells, {df.groupby(['sample','fov']).ngroups} FOVs")
    _log("Assigning CLDN4-high/low residuals ...")
    df = assign_high_low(df)

    rows = []
    for sample, sub in df.groupby("sample"):
        rows.append(
            {
                "sample": sample,
                "patient": sub["patient"].iloc[0],
                "n_cells": int(len(sub)),
                "n_fov": int(sub["fov"].nunique()),
                "n_tumor": int(sub["is_tumor"].sum()),
                "n_cd8": int(sub["is_cd8"].sum()),
                "n_epithelial": int(sub["is_epithelial"].sum()),
                "mean_cldn4_tumor": float(sub.loc[sub["is_tumor"], "CLDN4"].mean()),
                "mean_cldn4_other": float(sub.loc[~sub["is_tumor"], "CLDN4"].mean()),
                "mean_cldn4_cd8": float(sub.loc[sub["is_cd8"], "CLDN4"].mean()),
            }
        )
    by_s = pd.DataFrame(rows)
    by_s.to_csv(OUT_TAB / "cldn4_by_sample.csv", index=False)

    _log("Per-FOV spatial statistics + permutations ...")
    recs = []
    for (sample, fov_id), sub in df.groupby(["sample", "fov"]):
        rec = analyze_fov(sub)
        if rec is None:
            continue
        recs.append(rec)
        _log(
            f"  {sample} FOV {fov_id}: tumor={rec['n_tumor']} high={rec['n_high']} "
            f"cd8={rec['n_cd8']} dNN={rec['delta_median_nn']:.1f}"
        )
    if not recs:
        raise SystemExit("No FOV passed QC")

    keep_cols = [c for c in recs[0] if not isinstance(recs[0][c], np.ndarray)]
    fov = pd.DataFrame([{k: r[k] for k in keep_cols} for r in recs])
    fov.to_csv(OUT_TAB / "fov_spatial_stats.csv", index=False)
    _log(f"QC FOVs: {len(fov)}")

    _log("Writing figures ...")
    plot_maps(recs, OUT_FIG / "fig01_maps.png")
    plot_box_cldn4(df, OUT_FIG / "fig02_cldn4_tumor_vs_other.png")
    plot_paired_nn(fov, OUT_FIG / "fig03_nn_distance_paired.png")
    plot_ecdf(recs, OUT_FIG / "fig04_nn_distance_ecdf.png")
    plot_neighbor_bars(fov, OUT_FIG / "fig05_neighbor_counts.png")
    plot_k_envelope(recs, OUT_FIG / "fig06_K_hat_envelope.png", "K")
    plot_k_envelope(recs, OUT_FIG / "fig07_g_r_envelope.png", "g")
    plot_k_envelope(recs, OUT_FIG / "fig08_K_inhom_envelope.png", "Kinh")
    plot_core_margin(fov, OUT_FIG / "fig09_core_margin.png")
    plot_contact(fov, OUT_FIG / "fig10_contact_vs_shell.png")

    write_results(df, fov, recs, {"genes_present": present})
    _log(f"Wrote {OUT_RES}")
    print((OUT_TAB / "summary_key_numbers.json").read_text())


if __name__ == "__main__":
    main()
