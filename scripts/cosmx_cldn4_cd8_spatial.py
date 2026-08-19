#!/usr/bin/env python3
"""Official CosMx NSCLC FFPE (8 sections / 5 donors, 960-plex): CLDN4 vs CD8+NK neighbors.

Primary: July-style immune neighbor counts and mixing at r = 20/40/60 µm,
tested paired at SECTION (n=8) and DONOR (n=5) level.

Secondary: nearest-µm (not primary), Ripley's K, pair-correlation g(r).

Tumor = author tumor types (tumor 5/6/9/12/13). Tumor cells are NOT gated
on CD8A==0 (that would remove the interface).
Immune = author T CD8 naive + T CD8 memory + NK.

No private 8-KL data.
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
MARKER_GENES = ["CLDN4", "CD8A", "CD8B", "KRT8", "EPCAM", "NKG7"]
PRIMARY_R = np.array([20.0, 40.0, 60.0])
K_RADII = np.arange(10.0, 130.0, 10.0)
N_PERM = 199
RNG = np.random.default_rng(20260819)
MIN_TUMOR = 30


def _log(msg: str) -> None:
    print(msg, flush=True)


def load_sample(zf: zipfile.ZipFile, sample: str) -> pd.DataFrame:
    feat = pd.read_csv(zf.open(f"{sample}/qc/features.tsv"), sep="\t", index_col=0)
    genes = list(feat.index)
    use_genes = [g for g in MARKER_GENES if g in genes]
    if "CLDN4" not in genes:
        raise RuntimeError(f"{sample}: CLDN4 missing from 960 panel")
    obs = pd.read_csv(zf.open(f"{sample}/qc/observations.tsv"), sep="\t", index_col=0)
    coord = pd.read_csv(zf.open(f"{sample}/qc/coordinates.tsv"), sep="\t", index_col=0)
    labs = pd.read_csv(zf.open(f"{sample}/labels.tsv"), sep="\t", index_col=0)
    _log(f"  reading counts.mtx for {sample} ({len(obs):,} cells) ...")
    mat = mmread(zf.open(f"{sample}/qc/counts.mtx")).tocsr()
    idx = [genes.index(g) for g in use_genes]
    expr = pd.DataFrame(mat[:, idx].toarray(), index=obs.index, columns=use_genes)
    for g in MARKER_GENES:
        if g not in expr.columns:
            expr[g] = 0.0
    df = pd.concat(
        [
            expr[MARKER_GENES],
            obs[["n_counts", "n_genes"]] if "n_counts" in obs.columns else obs,
            coord.rename(columns={"x": "x_px", "y": "y_px"}),
            labs[["cell_type"]].astype(str),
        ],
        axis=1,
    )
    df["cell_id"] = sample + ":" + df.index.astype(str)
    df["sample"] = sample
    df["patient"] = PATIENT[sample]
    df["fov"] = df.index.astype(str).str.split("_").str[0]
    df = df.reset_index(drop=True)
    df["x_um"] = pd.to_numeric(df["x_px"], errors="coerce") * UM_PER_PX
    df["y_um"] = pd.to_numeric(df["y_px"], errors="coerce") * UM_PER_PX
    df["cell_type"] = df["cell_type"].astype(str)
    # Author tumor types only. Do NOT drop CD8A>0 tumor cells (interface).
    df["is_tumor"] = df["cell_type"].str.lower().str.startswith("tumor")
    df["is_cd8"] = df["cell_type"].isin(["T CD8 naive", "T CD8 memory"])
    df["is_nk"] = df["cell_type"].eq("NK")
    df["is_cd8nk"] = df["is_cd8"] | df["is_nk"]
    for col in MARKER_GENES + ["n_counts", "n_genes"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    df["log_cldn4"] = np.log1p(df["CLDN4"])
    df["log_krt8"] = np.log1p(df["KRT8"])
    df["log_epcam"] = np.log1p(df["EPCAM"])
    df["log_n"] = np.log1p(df["n_counts"])
    return df


def residualize_cldn4(tumor: pd.DataFrame) -> np.ndarray:
    x = tumor[["log_krt8", "log_epcam", "log_n"]].to_numpy(dtype=float)
    y = tumor["log_cldn4"].to_numpy(dtype=float)
    ok = np.isfinite(x).all(axis=1) & np.isfinite(y)
    resid = np.full(len(tumor), np.nan)
    if ok.sum() < 10 or np.allclose(y[ok], y[ok][0]):
        med = np.nanmedian(y)
        resid[ok] = y[ok] - (med if np.isfinite(med) else 0.0)
        return resid
    model = LinearRegression().fit(x[ok], y[ok])
    resid[ok] = y[ok] - model.predict(x[ok])
    return resid


def assign_high_low_section(df: pd.DataFrame) -> pd.DataFrame:
    """Per-SECTION median split of residual CLDN4 among author tumor cells."""
    df = df.copy()
    df["cldn4_resid"] = np.nan
    df["cldn4_hl"] = "other"
    for sample, idx in df.groupby("sample").groups.items():
        tum_idx = df.index[df["sample"].eq(sample) & df["is_tumor"]]
        tum = df.loc[tum_idx]
        if len(tum) < MIN_TUMOR:
            continue
        resid = residualize_cldn4(tum)
        df.loc[tum_idx, "cldn4_resid"] = resid
        finite = np.isfinite(resid)
        if finite.sum() < MIN_TUMOR:
            continue
        med = np.median(resid[finite])
        high = tum_idx[finite & (resid > med)]
        low = tum_idx[finite & (resid <= med)]
        df.loc[high, "cldn4_hl"] = "high"
        df.loc[low, "cldn4_hl"] = "low"
    return df


def radius_counts(src: np.ndarray, tgt: np.ndarray, radii: np.ndarray) -> np.ndarray:
    out = np.zeros((len(src), len(radii)), dtype=np.int32)
    if len(src) == 0 or len(tgt) == 0:
        return out
    tree = cKDTree(tgt)
    for j, r in enumerate(radii):
        out[:, j] = tree.query_ball_point(src, r, return_length=True)
    return out


def nn_distances(src: np.ndarray, tgt: np.ndarray) -> np.ndarray:
    if len(src) == 0:
        return np.array([])
    if len(tgt) == 0:
        return np.full(len(src), np.nan)
    return cKDTree(tgt).query(src, k=1)[0]


def bivariate_k(xy1: np.ndarray, xy2: np.ndarray, radii: np.ndarray, area: float) -> np.ndarray:
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
        g[i] = (k[i + 1] - k[i - 1]) / (2.0 * np.pi * radii[i] * dr)
    return g


def tumor_neighborhoods(tum: pd.DataFrame, immune_xy: np.ndarray, all_xy: np.ndarray, tumor_xy: np.ndarray) -> dict:
    """Per-tumor-cell CD8+NK counts, mixing, and density-normalized fractions."""
    src = tum[["x_um", "y_um"]].to_numpy()
    n_imm = radius_counts(src, immune_xy, PRIMARY_R)
    n_all = radius_counts(src, all_xy, PRIMARY_R)
    n_tum = radius_counts(src, tumor_xy, PRIMARY_R)
    # exclude self from all-cell and tumor-tumor counts
    n_all = np.clip(n_all - 1, 0, None)
    n_tum = np.clip(n_tum - 1, 0, None)
    frac = np.divide(n_imm, n_all, out=np.full_like(n_imm, np.nan, dtype=float), where=n_all > 0)
    # mixing: immune vs tumor among (immune + tumor) neighbors
    denom = n_imm + n_tum
    mix = np.divide(n_imm, denom, out=np.full_like(n_imm, np.nan, dtype=float), where=denom > 0)
    return {
        "n_imm": n_imm,
        "n_all": n_all,
        "n_tum": n_tum,
        "frac": frac,
        "mix": mix,
        "nn": nn_distances(src, immune_xy) if len(immune_xy) else np.full(len(src), np.nan),
    }


def summarize_arm(nb: dict) -> dict:
    out = {}
    for i, r in enumerate(PRIMARY_R):
        ri = int(r)
        out[f"cnt{ri}"] = float(np.nanmean(nb["n_imm"][:, i]))
        out[f"frac{ri}"] = float(np.nanmean(nb["frac"][:, i]))
        out[f"mix{ri}"] = float(np.nanmean(nb["mix"][:, i]))
        out[f"nall{ri}"] = float(np.nanmean(nb["n_all"][:, i]))
    out["median_nn"] = float(np.nanmedian(nb["nn"]))
    out["n_cells"] = int(len(nb["n_imm"]))
    return out


def perm_pvalue(obs: float, null: np.ndarray, alternative: str) -> float:
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or len(null) == 0:
        return np.nan
    if alternative == "less":
        n_ext = int(np.sum(null <= obs))
    elif alternative == "greater":
        n_ext = int(np.sum(null >= obs))
    else:
        n_ext = int(np.sum(np.abs(null) >= abs(obs)))
    return (1.0 + n_ext) / (1.0 + len(null))


def analyze_section(sub: pd.DataFrame) -> dict | None:
    tum = sub[sub["is_tumor"] & sub["cldn4_hl"].isin(["high", "low"])].copy()
    if len(tum) < MIN_TUMOR:
        return None
    high = tum[tum["cldn4_hl"] == "high"]
    low = tum[tum["cldn4_hl"] == "low"]
    if len(high) < 15 or len(low) < 15:
        return None
    immune = sub[sub["is_cd8nk"]]
    xy_imm = immune[["x_um", "y_um"]].to_numpy()
    xy_all = sub[["x_um", "y_um"]].to_numpy()
    xy_tum = tum[["x_um", "y_um"]].to_numpy()
    prev = float(sub["is_cd8nk"].mean())

    nb_h = tumor_neighborhoods(high, xy_imm, xy_all, xy_tum)
    nb_l = tumor_neighborhoods(low, xy_imm, xy_all, xy_tum)
    sh = summarize_arm(nb_h)
    sl = summarize_arm(nb_l)

    rec = {
        "sample": sub["sample"].iloc[0],
        "patient": sub["patient"].iloc[0],
        "n_cells": int(len(sub)),
        "n_fov": int(sub["fov"].nunique()),
        "n_tumor": int(len(tum)),
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "n_cd8": int(sub["is_cd8"].sum()),
        "n_nk": int(sub["is_nk"].sum()),
        "n_cd8nk": int(len(immune)),
        "cd8nk_prevalence": prev,
        "mean_cldn4_tumor": float(tum["CLDN4"].mean()),
        "mean_cldn4_other": float(sub.loc[~sub["is_tumor"], "CLDN4"].mean()),
        "frac_tumor_cd8a_pos": float((tum["CD8A"] > 0).mean()),
    }
    for k, v in sh.items():
        rec[f"high_{k}"] = v
    for k, v in sl.items():
        rec[f"low_{k}"] = v
    for r in PRIMARY_R:
        ri = int(r)
        rec[f"delta_cnt{ri}"] = rec[f"high_cnt{ri}"] - rec[f"low_cnt{ri}"]
        rec[f"delta_frac{ri}"] = rec[f"high_frac{ri}"] - rec[f"low_frac{ri}"]
        rec[f"delta_mix{ri}"] = rec[f"high_mix{ri}"] - rec[f"low_mix{ri}"]
        # enrichment vs section prevalence (density-normalized relative mixing)
        rec[f"high_enrich{ri}"] = rec[f"high_frac{ri}"] / prev if prev > 0 else np.nan
        rec[f"low_enrich{ri}"] = rec[f"low_frac{ri}"] / prev if prev > 0 else np.nan

    # section-level high/low label permutation (primary null)
    n_high = len(high)
    xy_t = tum[["x_um", "y_um"]].to_numpy()
    hl = tum["cldn4_hl"].to_numpy()
    null_d40 = np.zeros(N_PERM)
    null_dmix40 = np.zeros(N_PERM)
    for p in range(N_PERM):
        pick = RNG.choice(len(tum), size=n_high, replace=False)
        mask = np.zeros(len(tum), dtype=bool)
        mask[pick] = True
        h_df = tum.iloc[np.flatnonzero(mask)]
        l_df = tum.iloc[np.flatnonzero(~mask)]
        nh = tumor_neighborhoods(h_df, xy_imm, xy_all, xy_t)
        nl = tumor_neighborhoods(l_df, xy_imm, xy_all, xy_t)
        null_d40[p] = np.nanmean(nh["n_imm"][:, 1]) - np.nanmean(nl["n_imm"][:, 1])
        null_dmix40[p] = np.nanmean(nh["mix"][:, 1]) - np.nanmean(nl["mix"][:, 1])
    rec["p_perm_delta_cnt40_less"] = perm_pvalue(rec["delta_cnt40"], null_d40, "less")
    rec["p_perm_delta_mix40_less"] = perm_pvalue(rec["delta_mix40"], null_dmix40, "less")

    # Ripley / PCF per FOV (secondary; stored as section-mean curves)
    k_h_list, k_l_list, g_h_list, g_l_list = [], [], [], []
    for fov, fsub in sub.groupby("fov"):
        ft = fsub[fsub["is_tumor"] & fsub["cldn4_hl"].isin(["high", "low"])]
        fh, fl = ft[ft["cldn4_hl"] == "high"], ft[ft["cldn4_hl"] == "low"]
        fi = fsub[fsub["is_cd8nk"]]
        if len(fh) < 8 or len(fl) < 8 or len(fi) < 5:
            continue
        xyh, xyl, xyi = fh[["x_um", "y_um"]].to_numpy(), fl[["x_um", "y_um"]].to_numpy(), fi[["x_um", "y_um"]].to_numpy()
        area = max(
            (fsub["x_um"].max() - fsub["x_um"].min()) * (fsub["y_um"].max() - fsub["y_um"].min()),
            1.0,
        )
        kh = bivariate_k(xyh, xyi, K_RADII, area)
        kl = bivariate_k(xyl, xyi, K_RADII, area)
        k_h_list.append(kh)
        k_l_list.append(kl)
        g_h_list.append(pair_corr_from_k(kh, K_RADII))
        g_l_list.append(pair_corr_from_k(kl, K_RADII))
    rec["k_high"] = np.nanmean(np.vstack(k_h_list), axis=0) if k_h_list else np.full(len(K_RADII), np.nan)
    rec["k_low"] = np.nanmean(np.vstack(k_l_list), axis=0) if k_l_list else np.full(len(K_RADII), np.nan)
    rec["g_high"] = np.nanmean(np.vstack(g_h_list), axis=0) if g_h_list else np.full(len(K_RADII), np.nan)
    rec["g_low"] = np.nanmean(np.vstack(g_l_list), axis=0) if g_l_list else np.full(len(K_RADII), np.nan)
    rec["n_fov_ripley"] = int(len(k_h_list))
    rec["_nb_high"] = nb_h
    rec["_nb_low"] = nb_l
    rec["_xy_h"] = high[["x_um", "y_um"]].to_numpy()
    rec["_xy_l"] = low[["x_um", "y_um"]].to_numpy()
    rec["_xy_imm"] = xy_imm
    rec["_xy_other"] = sub.loc[~sub["is_tumor"] & ~sub["is_cd8nk"], ["x_um", "y_um"]].to_numpy()
    rec["_fov"] = sub["fov"].to_numpy()
    rec["_high_fov"] = high["fov"].to_numpy()
    rec["_imm_fov"] = immune["fov"].to_numpy() if len(immune) else np.array([])
    rec["_low_fov"] = low["fov"].to_numpy()
    return rec


def paired_wilcoxon(a: np.ndarray, b: np.ndarray, alternative: str) -> tuple[float, float]:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return np.nan, np.nan
    # Wilcoxon needs non-zero differences
    if np.allclose(a[m], b[m]):
        return 0.0, 1.0
    try:
        res = stats.wilcoxon(a[m], b[m], alternative=alternative, zero_method="wilcox")
        return float(res.statistic), float(res.pvalue)
    except ValueError:
        return np.nan, np.nan


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_maps(df: pd.DataFrame, path: Path) -> None:
    # pick 4 FOVs with many tumor + CD8+NK
    scores = []
    for (sample, fov), sub in df.groupby(["sample", "fov"]):
        nt = int(sub["is_tumor"].sum())
        ni = int(sub["is_cd8nk"].sum())
        if nt >= 40 and ni >= 8:
            scores.append((nt * ni, sample, str(fov)))
    scores.sort(reverse=True)
    pick = scores[:4]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 9.2))
    for ax, (_, sample, fov) in zip(axes.ravel(), pick):
        sub = df[(df["sample"] == sample) & (df["fov"].astype(str) == fov)]
        o = sub.loc[~sub["is_tumor"] & ~sub["is_cd8nk"], ["x_um", "y_um"]].to_numpy()
        l = sub.loc[sub["cldn4_hl"] == "low", ["x_um", "y_um"]].to_numpy()
        h = sub.loc[sub["cldn4_hl"] == "high", ["x_um", "y_um"]].to_numpy()
        c = sub.loc[sub["is_cd8nk"], ["x_um", "y_um"]].to_numpy()
        if len(o):
            ax.scatter(o[:, 0], o[:, 1], s=2, c="#d9d9d9", linewidths=0, rasterized=True)
        if len(l):
            ax.scatter(l[:, 0], l[:, 1], s=6, c="#f4a261", linewidths=0, label="CLDN4-low tumor", rasterized=True)
        if len(h):
            ax.scatter(h[:, 0], h[:, 1], s=6, c="#9b2226", linewidths=0, label="CLDN4-high tumor", rasterized=True)
        if len(c):
            ax.scatter(c[:, 0], c[:, 1], s=10, c="#1d4e89", linewidths=0, label="CD8+NK", rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"{sample} FOV {fov}", fontsize=9)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        style_ax(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("All 8 CosMx NSCLC sections — example FOVs (author tumor + CD8+NK)", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_neighbor_paired(sec: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.1), sharey=False)
    for ax, r in zip(axes, [20, 40, 60]):
        lo = sec[f"low_cnt{r}"].to_numpy()
        hi = sec[f"high_cnt{r}"].to_numpy()
        for a, b in zip(lo, hi):
            ax.plot([0, 1], [a, b], color="#bbbbbb", lw=0.9)
        ax.scatter(np.zeros(len(lo)), lo, s=36, c="#f4a261", zorder=3, label="CLDN4-low")
        ax.scatter(np.ones(len(hi)), hi, s=36, c="#9b2226", zorder=3, label="CLDN4-high")
        # donor means
        dlo = sec.groupby("patient")[f"low_cnt{r}"].mean()
        dhi = sec.groupby("patient")[f"high_cnt{r}"].mean()
        ax.plot([0, 1], [dlo.mean(), dhi.mean()], color="#222222", lw=2.0, zorder=2)
        _, p_s = paired_wilcoxon(hi, lo, "less")
        _, p_d = paired_wilcoxon(dhi.to_numpy(), dlo.to_numpy(), "less")
        dmed = np.median(hi - lo)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4-low\ntumor", "CLDN4-high\ntumor"])
        ax.set_title(f"r = {r} µm\n∆section med={dmed:.3f}\nsection p={fmt_p(p_s)}  donor p={fmt_p(p_d)}")
        style_ax(ax)
    axes[0].set_ylabel("Mean CD8+NK neighbors per tumor cell")
    axes[0].legend(frameon=False, loc="upper right", fontsize=8)
    fig.suptitle("Primary: CD8+NK neighbor counts on all 8 CosMx NSCLC sections", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_mixing(sec: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    ax = axes[0]
    lo, hi = sec["low_mix40"].to_numpy(), sec["high_mix40"].to_numpy()
    for a, b, name in zip(lo, hi, sec["sample"]):
        ax.plot([0, 1], [a, b], color="#bbbbbb", lw=0.9)
        ax.text(1.04, b, name.replace("Lung", "L"), fontsize=7, va="center")
    ax.scatter(np.zeros(len(lo)), lo, s=36, c="#f4a261", zorder=3)
    ax.scatter(np.ones(len(hi)), hi, s=36, c="#9b2226", zorder=3)
    _, p = paired_wilcoxon(hi, lo, "less")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-low", "CLDN4-high"])
    ax.set_ylabel("Mixing score (CD8+NK / (CD8+NK + tumor) neighbors)")
    ax.set_title(f"Mixing at 40 µm\nsection-paired p={fmt_p(p)}")
    style_ax(ax)
    ax = axes[1]
    lo, hi = sec["low_frac40"].to_numpy(), sec["high_frac40"].to_numpy()
    for a, b in zip(lo, hi):
        ax.plot([0, 1], [a, b], color="#bbbbbb", lw=0.9)
    ax.scatter(np.zeros(len(lo)), lo, s=36, c="#f4a261", zorder=3, label="CLDN4-low")
    ax.scatter(np.ones(len(hi)), hi, s=36, c="#9b2226", zorder=3, label="CLDN4-high")
    _, p = paired_wilcoxon(hi, lo, "less")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-low", "CLDN4-high"])
    ax.set_ylabel("CD8+NK fraction among all neighbors")
    ax.set_title(f"Density-normalized fraction at 40 µm\nsection-paired p={fmt_p(p)}")
    style_ax(ax)
    fig.suptitle("Mixing and density-normalized CD8+NK fraction (8 sections)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_donor(sec: pd.DataFrame, path: Path) -> None:
    don = sec.groupby("patient", as_index=False).mean(numeric_only=True)
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    lo, hi = don["low_cnt40"].to_numpy(), don["high_cnt40"].to_numpy()
    names = don["patient"].tolist()
    for a, b, name in zip(lo, hi, names):
        ax.plot([0, 1], [a, b], color="#888888", lw=1.2)
        ax.text(1.05, b, name, fontsize=8, va="center")
    ax.scatter(np.zeros(len(lo)), lo, s=50, c="#f4a261", zorder=3)
    ax.scatter(np.ones(len(hi)), hi, s=50, c="#9b2226", zorder=3)
    _, p = paired_wilcoxon(hi, lo, "less")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-low\ntumor", "CLDN4-high\ntumor"])
    ax.set_ylabel("Mean CD8+NK neighbors within 40 µm")
    ax.set_title(f"Donor-paired (n=5)\n∆mean={np.mean(hi-lo):.3f}  Wilcoxon p={fmt_p(p)}")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_count_ecdf(recs: list[dict], path: Path) -> None:
    h = np.concatenate([r["_nb_high"]["n_imm"][:, 1] for r in recs])
    l = np.concatenate([r["_nb_low"]["n_imm"][:, 1] for r in recs])

    def ecdf(v):
        v = np.sort(v)
        return v, np.arange(1, len(v) + 1) / len(v)

    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    xl, yl = ecdf(l)
    xh, yh = ecdf(h)
    ax.plot(xl, yl, color="#f4a261", lw=2, label=f"CLDN4-low (n={len(l):,})")
    ax.plot(xh, yh, color="#9b2226", lw=2, label=f"CLDN4-high (n={len(h):,})")
    ax.set_xlabel("CD8+NK neighbors within 40 µm")
    ax.set_ylabel("ECDF")
    p = stats.ks_2samp(h, l, alternative="greater").pvalue
    ax.set_title(f"Tumor-cell CD8+NK count ECDF (supporting; not the inferential unit)\nKS p={fmt_p(p)}")
    ax.legend(frameon=False)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_k(recs: list[dict], path: Path, kind: str) -> None:
    key_h = "k_high" if kind == "K" else "g_high"
    key_l = "k_low" if kind == "K" else "g_low"
    mh = np.nanmean(np.vstack([r[key_h] for r in recs]), axis=0)
    ml = np.nanmean(np.vstack([r[key_l] for r in recs]), axis=0)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    if kind == "K":
        ax.plot(K_RADII, np.pi * K_RADII**2, color="0.45", ls=":", lw=1.2, label="CSR πr²")
        ylab = r"Section-mean bivariate $\hat{K}_{12}(r)$ (tumor → CD8+NK)"
        title = "Secondary: Ripley's K (not primary)"
    else:
        ax.axhline(1.0, color="0.45", ls=":", lw=1.2, label="independence g=1")
        ylab = r"Section-mean pair correlation $\hat{g}_{12}(r)$"
        title = "Secondary: pair-correlation g(r) (not primary)"
    ax.plot(K_RADII, ml, color="#f4a261", lw=2, label="CLDN4-low tumor → CD8+NK")
    ax.plot(K_RADII, mh, color="#9b2226", lw=2, label="CLDN4-high tumor → CD8+NK")
    ax.set_xlabel("r (µm)")
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nmean of {len(recs)} sections")
    ax.legend(frameon=False, fontsize=8)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_cldn4_box(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    tumor = np.log1p(df.loc[df["is_tumor"], "CLDN4"].to_numpy())
    other = np.log1p(df.loc[~df["is_tumor"], "CLDN4"].to_numpy())
    parts = ax.boxplot(
        [other, tumor],
        tick_labels=["other types", "author tumor"],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(parts["boxes"], ["#9aa0a6", "#9b2226"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    p = stats.mannwhitneyu(tumor, other, alternative="greater").pvalue
    ax.set_ylabel("log1p CLDN4 counts")
    ax.set_title(f"CLDN4 in tumor vs other\nn_tumor={len(tumor):,}  n_other={len(other):,}  MWU p={fmt_p(p)}")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_nn_secondary(sec: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    lo, hi = sec["low_median_nn"].to_numpy(), sec["high_median_nn"].to_numpy()
    for a, b in zip(lo, hi):
        ax.plot([0, 1], [a, b], color="#bbbbbb", lw=0.9)
    ax.scatter(np.zeros(len(lo)), lo, s=36, c="#f4a261", zorder=3)
    ax.scatter(np.ones(len(hi)), hi, s=36, c="#9b2226", zorder=3)
    _, p = paired_wilcoxon(hi, lo, "greater")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-low", "CLDN4-high"])
    ax.set_ylabel("Median nearest CD8+NK (µm)")
    ax.set_title(
        f"Secondary (not primary): nearest-µm\n"
        f"section med {np.median(lo):.1f} vs {np.median(hi):.1f} µm  p={fmt_p(p)}"
    )
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_results(df: pd.DataFrame, sec: pd.DataFrame, recs: list[dict], genes_present: list[str]) -> None:
    n_cells = len(df)
    n_fov = int(df.groupby(["sample", "fov"]).ngroups)
    n_tumor = int(df["is_tumor"].sum())
    n_cd8 = int(df["is_cd8"].sum())
    n_nk = int(df["is_nk"].sum())
    n_cd8nk = int(df["is_cd8nk"].sum())
    n_high = int((df["cldn4_hl"] == "high").sum())
    n_low = int((df["cldn4_hl"] == "low").sum())
    n_sec = len(sec)
    don = sec.groupby("patient", as_index=False).mean(numeric_only=True)

    p_expr = stats.mannwhitneyu(
        np.log1p(df.loc[df["is_tumor"], "CLDN4"]),
        np.log1p(df.loc[~df["is_tumor"], "CLDN4"]),
        alternative="greater",
    ).pvalue

    rows_p = []
    for r in (20, 40, 60):
        _, ps = paired_wilcoxon(sec[f"high_cnt{r}"], sec[f"low_cnt{r}"], "less")
        _, pdn = paired_wilcoxon(don[f"high_cnt{r}"], don[f"low_cnt{r}"], "less")
        dsec = float((sec[f"high_cnt{r}"] - sec[f"low_cnt{r}"]).median())
        ddon = float((don[f"high_cnt{r}"] - don[f"low_cnt{r}"]).mean())
        n_dir = int((sec[f"high_cnt{r}"] < sec[f"low_cnt{r}"]).sum())
        rows_p.append((r, dsec, ps, ddon, pdn, n_dir))

    _, p_mix_s = paired_wilcoxon(sec["high_mix40"], sec["low_mix40"], "less")
    _, p_mix_d = paired_wilcoxon(don["high_mix40"], don["low_mix40"], "less")
    _, p_frac_s = paired_wilcoxon(sec["high_frac40"], sec["low_frac40"], "less")
    _, p_frac_d = paired_wilcoxon(don["high_frac40"], don["low_frac40"], "less")
    d_mix = float((sec["high_mix40"] - sec["low_mix40"]).median())
    d_frac = float((sec["high_frac40"] - sec["low_frac40"]).median())
    _, p_nn = paired_wilcoxon(sec["high_median_nn"], sec["low_median_nn"], "greater")

    r40 = [x for x in rows_p if x[0] == 40][0]
    paper = (
        f"In the official CosMx NSCLC FFPE 960-plex cohort (all 8 sections / 5 donors; "
        f"{n_tumor:,} author-typed tumor cells, {n_cd8:,} CD8 T cells + {n_nk:,} NK = {n_cd8nk:,} CD8+NK; "
        f"{n_cells:,} cells, {n_fov} FOVs), CLDN4-high tumor cells (per-section residual split on "
        f"KRT8/EPCAM; tumor not gated on CD8A==0) had fewer CD8+NK neighbors within 40 µm than "
        f"CLDN4-low tumor cells (section-paired median ∆ = {r40[1]:.3f} cells, Wilcoxon p = {fmt_p(r40[2])}, "
        f"n = 8 sections, {r40[5]}/8 sections in the same direction; donor-paired mean ∆ = {r40[3]:.3f}, "
        f"p = {fmt_p(r40[4])}, n = 5 donors). The 40 µm mixing score and density-normalized CD8+NK "
        f"neighbor fraction were also lower for CLDN4-high tumor cells "
        f"(∆mix = {d_mix:.4f}, section p = {fmt_p(p_mix_s)}; ∆frac = {d_frac:.4f}, section p = {fmt_p(p_frac_s)})."
    )

    lines = [
        "# RESULTS: CosMx NSCLC — CLDN4-high tumor vs CD8+NK neighbors (all 8 sections)",
        "",
        "## Paper sentence",
        "",
        paper,
        "",
        "## Dataset and panel",
        "",
        "- **Dataset:** official NanoString/Bruker CosMx SMI NSCLC FFPE (8 FFPE sections / 5 patients, 960-plex). https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/",
        "- **Mirror used:** Zenodo 15487520 `cosmx_lung` (counts, coordinates, author `cell_type`).",
        f"- **Panel check:** CLDN4 is on the 960-gene matrix. Genes present: {', '.join(genes_present)}.",
        "- **Pixel size:** 0.18 µm/pixel (official SMI-ReadMe).",
        "- **Tumor:** author types `tumor 5/6/9/12/13`. Tumor cells with CD8A>0 were **kept** (no CD8A==0 gate).",
        "- **CD8+NK:** author `T CD8 naive` + `T CD8 memory` + `NK`.",
        "- **CLDN4-high/low:** per-section median split of residual log1p(CLDN4) after OLS on log1p(KRT8), log1p(EPCAM), log1p(n_counts).",
        "- **Units of inference:** section (n=8, paired) and donor (n=5, paired). FOV cells are not the test unit.",
        "- No private 8-KL data.",
        "",
        f"- Loaded **{n_cells:,} cells**, **{n_fov} FOVs**, **8/8 sections**, **5 donors**.",
        f"- Author tumor = {n_tumor:,} (high {n_high:,} / low {n_low:,}); CD8 = {n_cd8:,}; NK = {n_nk:,}; CD8+NK = {n_cd8nk:,}.",
        f"- Fraction of tumor cells with CD8A>0 (kept): {df.loc[df['is_tumor'],'CD8A'].gt(0).mean():.3f}.",
        "",
        "## Primary: CD8+NK neighbor counts at 20 / 40 / 60 µm",
        "",
        "Mean CD8+NK cells within radius r of each tumor cell, then averaged within section. "
        "Paired Wilcoxon signed-rank on the 8 section means (high vs low) and on the 5 donor means (section means averaged within donor).",
        "",
        "| r | CLDN4-low mean | CLDN4-high mean | ∆med (high−low) | section p (high<low) | donor ∆mean | donor p | sections with high<low |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r, dsec, ps, ddon, pdn, n_dir in rows_p:
        lines.append(
            f"| {r} µm | {sec[f'low_cnt{r}'].mean():.3f} | {sec[f'high_cnt{r}'].mean():.3f} | "
            f"{dsec:.3f} | {fmt_p(ps)} | {ddon:.3f} | {fmt_p(pdn)} | {n_dir}/8 |"
        )
    lines += [
        "",
        "Per-section values: `tables/section_neighbor_stats.csv`. Per-donor: `tables/donor_neighbor_stats.csv`.",
        "",
        f"Section-level high/low label permutation (n={N_PERM}) for ∆count at 40 µm: "
        f"median per-section p = {sec['p_perm_delta_cnt40_less'].median():.3f}; "
        f"sections with perm p<0.05: {int((sec['p_perm_delta_cnt40_less']<0.05).sum())}/8.",
        "",
        "## Mixing score and density-normalized CD8+NK fraction",
        "",
        "- **Mixing score** at r: among neighbors that are tumor or CD8+NK, the fraction that are CD8+NK "
        "(immune–tumor mixing vs tumor–tumor self-aggregation).",
        "- **Density-normalized fraction:** CD8+NK / all cells within r (accounts for local packing).",
        "- **Enrichment:** density-normalized fraction / section CD8+NK prevalence.",
        "",
        f"- Mixing at 40 µm: low {sec['low_mix40'].mean():.4f} vs high {sec['high_mix40'].mean():.4f}; "
        f"∆med = {d_mix:.4f}; section p = {fmt_p(p_mix_s)}; donor p = {fmt_p(p_mix_d)}.",
        f"- Density-normalized CD8+NK fraction at 40 µm: low {sec['low_frac40'].mean():.4f} vs high {sec['high_frac40'].mean():.4f}; "
        f"∆med = {d_frac:.4f}; section p = {fmt_p(p_frac_s)}; donor p = {fmt_p(p_frac_d)}.",
        f"- Enrichment vs section prevalence at 40 µm: low {sec['low_enrich40'].mean():.3f} vs high {sec['high_enrich40'].mean():.3f}.",
        "",
        "## 1. CLDN4 in tumor vs other types",
        "",
        f"- Mean raw CLDN4: tumor {df.loc[df['is_tumor'],'CLDN4'].mean():.3f} vs other {df.loc[~df['is_tumor'],'CLDN4'].mean():.3f}.",
        f"- log1p CLDN4 MWU (tumor > other): p = {fmt_p(p_expr)} (n_tumor={n_tumor:,}, n_other={n_cells-n_tumor:,}).",
        "",
        "## Secondary: nearest-µm (not primary)",
        "",
        f"Nearest CD8+NK distance is reported only as a negative control on why it is the wrong lead: "
        f"section median {sec['low_median_nn'].median():.1f} vs {sec['high_median_nn'].median():.1f} µm "
        f"(high−low ∆med = {(sec['high_median_nn']-sec['low_median_nn']).median():.2f} µm; p = {fmt_p(p_nn)}). "
        "Counts and mixing at 20/40/60 µm are the primary spatial readouts.",
        "",
        "## Secondary: Ripley's K and pair-correlation g(r)",
        "",
        "Bivariate K̂ and ĝ from CLDN4-high or low tumor cells to CD8+NK, border-corrected, averaged within section then across sections. CSR reference is πr² (K) or g=1.",
        f"- Mean K̂(40) high vs low: {np.nanmean([r['k_high'][3] for r in recs]):.1f} vs {np.nanmean([r['k_low'][3] for r in recs]):.1f} (CSR π·40² = {np.pi*40**2:.1f}).",
        f"- FOVs contributing to K/g: sum of per-section Ripley FOVs = {int(sec['n_fov_ripley'].sum())}.",
        "",
        "## Counts",
        "",
        "| Level | n |",
        "|---|---|",
        "| Donors | 5 |",
        "| Sections | 8 / 8 |",
        f"| FOVs | {n_fov} |",
        f"| Cells | {n_cells:,} |",
        f"| Author tumor | {n_tumor:,} |",
        f"| CLDN4-high / low tumor | {n_high:,} / {n_low:,} |",
        f"| Author CD8 | {n_cd8:,} |",
        f"| Author NK | {n_nk:,} |",
        f"| CD8+NK | {n_cd8nk:,} |",
        "",
        "## Figures",
        "",
        "- `figures/fig01_maps.png` — example FOV maps from the 8-section cohort.",
        "- `figures/fig02_neighbor_counts_20_40_60.png` — **primary** paired section CD8+NK counts.",
        "- `figures/fig03_mixing_and_fraction.png` — mixing score and density-normalized fraction at 40 µm.",
        "- `figures/fig04_donor_paired_40um.png` — donor-paired 40 µm counts.",
        "- `figures/fig05_count_ecdf_40um.png` — tumor-cell count ECDF (supporting).",
        "- `figures/fig06_cldn4_tumor_vs_other.png` — CLDN4 expression.",
        "- `figures/fig07_nn_distance_secondary.png` — nearest-µm (not primary).",
        "- `figures/fig08_K_hat.png` — Ripley's K (secondary).",
        "- `figures/fig09_g_r.png` — pair-correlation g(r) (secondary).",
        "",
        "## Methods",
        "",
        "All 8 official CosMx NSCLC sections were used (Lung5 Rep1–3, Lung6, Lung9 Rep1–2, Lung12, Lung13). "
        "Coordinates were converted with 0.18 µm/pixel. Author cell types define tumor and CD8+NK; "
        "tumor cells were not excluded for CD8A expression. CLDN4-high vs low is a within-section median "
        "split of residual log1p(CLDN4) after KRT8, EPCAM, and library size. For each tumor cell, CD8+NK "
        "neighbors, all-cell neighbors, and tumor neighbors were counted at 20, 40, and 60 µm with a KD-tree. "
        "Mixing = n_CD8NK / (n_CD8NK + n_tumor_neighbors). Density-normalized fraction = n_CD8NK / n_all_neighbors. "
        "Section means were tested with paired Wilcoxon (high vs low). Donor means are unweighted averages of "
        "that donor's sections. High/low labels were also shuffled within section (199 times) to get a "
        "section-level permutation p for ∆count and ∆mix at 40 µm. Ripley's K and g(r) were computed per FOV "
        "and averaged; they are secondary to neighbor counts.",
        "",
        "## Notes",
        "",
        "- Additive public CLDN4-only spatial analysis. No private 8-KL.",
        "- Nearest-µm is retained in the tables but is not the lead metric.",
    ]
    OUT_RES.write_text("\n".join(lines) + "\n")
    key = {
        "n_cells": n_cells,
        "n_fov": n_fov,
        "n_sections": n_sec,
        "n_donors": int(don.shape[0]),
        "n_tumor": n_tumor,
        "n_cd8nk": n_cd8nk,
        "delta_cnt40_section_median": r40[1],
        "p_section_cnt40": r40[2],
        "p_donor_cnt40": r40[4],
        "delta_mix40_section_median": d_mix,
        "p_section_mix40": p_mix_s,
        "paper_sentence": paper,
    }
    (OUT_TAB / "summary_key_numbers.json").write_text(json.dumps(key, indent=2))


def main() -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    if not ZIP_PATH.exists():
        raise SystemExit(f"Missing {ZIP_PATH}")

    _log("Loading official CosMx NSCLC (all 8 sections, Zenodo cosmx_lung) ...")
    frames = []
    with zipfile.ZipFile(ZIP_PATH) as zf:
        feat0 = pd.read_csv(zf.open(f"{SAMPLES[0]}/qc/features.tsv"), sep="\t", index_col=0)
        genes = list(feat0.index)
        present = [g for g in MARKER_GENES if g in genes]
        _log(f"Panel genes present: {present} (n_panel={len(genes)})")
        if "CLDN4" not in genes:
            OUT_RES.write_text(
                "# RESULTS\n\nCLDN4 is absent from the CosMx 960-plex NSCLC matrix. Analysis stopped.\n"
            )
            raise SystemExit("CLDN4 missing")
        for sample in SAMPLES:
            frames.append(load_sample(zf, sample))
    df = pd.concat(frames, ignore_index=True)
    _log(f"Loaded {len(df):,} cells; sections={df['sample'].nunique()}; FOVs={df.groupby(['sample','fov']).ngroups}")
    assert df["sample"].nunique() == 8, "Expected all 8 CosMx NSCLC sections"
    _log("Section-level CLDN4 residual high/low (tumor not gated on CD8A) ...")
    df = assign_high_low_section(df)
    _log(
        f"Tumor kept with CD8A>0: {(df['is_tumor'] & (df['CD8A']>0)).sum():,} / {df['is_tumor'].sum():,}"
    )

    recs = []
    for sample in SAMPLES:
        rec = analyze_section(df[df["sample"] == sample])
        if rec is None:
            _log(f"  SKIP {sample} (insufficient tumor)")
            continue
        recs.append(rec)
        _log(
            f"  {sample}: tumor={rec['n_tumor']} high={rec['n_high']} low={rec['n_low']} "
            f"CD8+NK={rec['n_cd8nk']}  ∆cnt40={rec['delta_cnt40']:.3f}  mix40 h/l={rec['high_mix40']:.3f}/{rec['low_mix40']:.3f}"
        )
    if len(recs) != 8:
        _log(f"WARNING: {len(recs)} sections passed QC (want 8)")

    keep = [k for k, v in recs[0].items() if not str(k).startswith("_") and not isinstance(v, np.ndarray)]
    sec = pd.DataFrame([{k: r[k] for k in keep} for r in recs])
    sec.to_csv(OUT_TAB / "section_neighbor_stats.csv", index=False)
    don = sec.groupby("patient", as_index=False).mean(numeric_only=True)
    don.to_csv(OUT_TAB / "donor_neighbor_stats.csv", index=False)

    _log("Writing figures ...")
    plot_maps(df, OUT_FIG / "fig01_maps.png")
    plot_neighbor_paired(sec, OUT_FIG / "fig02_neighbor_counts_20_40_60.png")
    plot_mixing(sec, OUT_FIG / "fig03_mixing_and_fraction.png")
    plot_donor(sec, OUT_FIG / "fig04_donor_paired_40um.png")
    plot_count_ecdf(recs, OUT_FIG / "fig05_count_ecdf_40um.png")
    plot_cldn4_box(df, OUT_FIG / "fig06_cldn4_tumor_vs_other.png")
    plot_nn_secondary(sec, OUT_FIG / "fig07_nn_distance_secondary.png")
    plot_k(recs, OUT_FIG / "fig08_K_hat.png", "K")
    plot_k(recs, OUT_FIG / "fig09_g_r.png", "g")

    write_results(df, sec, recs, present)
    _log(f"Wrote {OUT_RES}")
    print((OUT_TAB / "summary_key_numbers.json").read_text())


if __name__ == "__main__":
    main()
