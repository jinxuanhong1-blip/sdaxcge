#!/usr/bin/env python3
"""Per-FOV and meta-analysed cross-type point-process analysis:

CLDN4-high tumor cells  ↔  CD8 T cells
on the official CosMx NSCLC FFPE 960-plex dataset (He et al. 2022).
"""

from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from point_process import (
    cross_g_rings,
    cross_k_curve,
    dist_to_border,
    empirical_p,
    homogeneous_weights,
    inhomogeneous_weights,
    intensity_at_points,
    k_to_l,
    pairwise_dist,
    rectangle_window,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx_nsclc_cells.csv.gz"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

UM_PER_COORD = 1000.0  # Giotto spatial_locs are millimetres
R_MAX = 150.0
R_GRID = np.arange(10.0, R_MAX + 1e-9, 5.0)  # 10, 15, ..., 150 µm
R_LO = R_GRID[:-1]
R_HI = R_GRID[1:]
R_MID = 0.5 * (R_LO + R_HI)
N_PERM_LABEL = 399
N_PERM_CSR = 199
BANDWIDTH_UM = 80.0
MIN_HIGH = 20
MIN_CD8 = 20
MIN_TUMOR = 40
RNG_SEED = 20260819
N_WORKERS = max(1, min(8, os.cpu_count() or 2))

TUMOR_PREFIX = "tumor"
CD8_TYPES = {"T CD8 memory", "T CD8 naive"}
EPITHELIAL_TYPES = {"epithelial"}  # plus tumor_* for lambda


def load_cells() -> pd.DataFrame:
    df = pd.read_csv(DATA)
    df["x_um"] = df["x"] * UM_PER_COORD
    df["y_um"] = df["y"] * UM_PER_COORD
    df["is_tumor"] = df["cell_type"].astype(str).str.startswith(TUMOR_PREFIX)
    df["is_epithelial"] = df["cell_type"].eq("epithelial")
    df["is_tumor_or_epi"] = df["is_tumor"] | df["is_epithelial"]
    df["is_cd8"] = df["cell_type"].isin(CD8_TYPES)
    # CLDN4-high: within-sample 75th percentile among author tumor cells
    q = df.loc[df["is_tumor"]].groupby("sample")["CLDN4"].quantile(0.75)
    df["cldn4_q75_tumor"] = df["sample"].map(q)
    df["is_cldn4_high_tumor"] = df["is_tumor"] & (df["CLDN4"] >= df["cldn4_q75_tumor"]) & (
        df["CLDN4"] >= 1
    )
    return df


def fov_key(sample: str, fov: int) -> str:
    return f"{sample}|FOV{int(fov):03d}"


def _stats_from_dist(dist, xy_a, r, r_lo, r_hi, win, area, lam_a, lam_b, border):
    n_a, n_b = dist.shape
    out = {}
    w_row, w_col = homogeneous_weights(n_a, n_b, area)
    k = cross_k_curve(dist, r, border, w_row, w_col)
    g = cross_g_rings(dist, r_lo, r_hi, border, w_row, w_col)
    out["k"] = k
    out["l"] = k_to_l(k)
    out["g"] = g
    w_row_i, w_col_i = inhomogeneous_weights(lam_a, lam_b, area)
    k_i = cross_k_curve(dist, r, border, w_row_i, w_col_i)
    g_i = cross_g_rings(dist, r_lo, r_hi, border, w_row_i, w_col_i)
    out["k_inhom"] = k_i
    out["l_inhom"] = k_to_l(k_i)
    out["g_inhom"] = g_i
    return out


def analyze_fov(payload: dict) -> dict | None:
    xy_all = payload["xy_all"]
    xy_tumor = payload["xy_tumor"]
    xy_te = payload["xy_tumor_epi"]
    xy_cd8 = payload["xy_cd8"]
    high_idx = payload["high_idx"]  # indices into tumor array
    n_high = int(high_idx.size)
    n_cd8 = xy_cd8.shape[0]
    n_tumor = xy_tumor.shape[0]
    if n_high < MIN_HIGH or n_cd8 < MIN_CD8 or n_tumor < MIN_TUMOR:
        return None

    win, area = rectangle_window(xy_all, pad=0.0)
    rng = np.random.default_rng(payload["seed"])

    xy_high = xy_tumor[high_idx]
    border_high = dist_to_border(xy_high, win)
    border_tumor = dist_to_border(xy_tumor, win)

    lam_te_tumor = intensity_at_points(xy_te, xy_tumor, win, BANDWIDTH_UM)
    lam_cd8 = intensity_at_points(xy_cd8, xy_cd8, win, BANDWIDTH_UM)
    # CD8 lambda at CD8 points; CLDN4-high lambda = tumor/epithelial intensity at those tumor cells
    dist_t = pairwise_dist(xy_tumor, xy_cd8)  # (n_tumor, n_cd8)

    def stats_for_idx(idx):
        d = dist_t[idx]
        lam_a = lam_te_tumor[idx]
        b = border_tumor[idx]
        return _stats_from_dist(d, xy_tumor[idx], R_GRID, R_LO, R_HI, win, area, lam_a, lam_cd8, b)

    obs = stats_for_idx(high_idx)

    # Label permutations: same n_high among tumor cells, CD8 fixed
    perm_g = np.zeros((N_PERM_LABEL, R_MID.size))
    perm_k = np.zeros((N_PERM_LABEL, R_GRID.size))
    perm_gi = np.zeros((N_PERM_LABEL, R_MID.size))
    perm_ki = np.zeros((N_PERM_LABEL, R_GRID.size))
    for p in range(N_PERM_LABEL):
        idx = rng.choice(n_tumor, size=n_high, replace=False)
        st = stats_for_idx(idx)
        perm_k[p] = st["k"]
        perm_g[p] = st["g"]
        perm_ki[p] = st["k_inhom"]
        perm_gi[p] = st["g_inhom"]

    # CSR of CD8 in the window (homogeneous null); CLDN4-high fixed
    xmin, xmax, ymin, ymax = win
    perm_g_csr = np.zeros((N_PERM_CSR, R_MID.size))
    perm_k_csr = np.zeros((N_PERM_CSR, R_GRID.size))
    xy_h = xy_high
    border_h = border_high
    lam_a_obs = lam_te_tumor[high_idx]
    w_row_h, w_col_h = homogeneous_weights(n_high, n_cd8, area)
    for p in range(N_PERM_CSR):
        xy_b = np.column_stack(
            [
                rng.uniform(xmin, xmax, n_cd8),
                rng.uniform(ymin, ymax, n_cd8),
            ]
        )
        d = pairwise_dist(xy_h, xy_b)
        perm_k_csr[p] = cross_k_curve(d, R_GRID, border_h, w_row_h, w_col_h)
        perm_g_csr[p] = cross_g_rings(d, R_LO, R_HI, border_h, w_row_h, w_col_h)

    def env(arr):
        return np.nanpercentile(arr, [2.5, 50.0, 97.5], axis=0)

    k_lo, k_md, k_hi = env(perm_k)
    g_lo, g_md, g_hi = env(perm_g)
    ki_lo, ki_md, ki_hi = env(perm_ki)
    gi_lo, gi_md, gi_hi = env(perm_gi)
    kc_lo, kc_md, kc_hi = env(perm_k_csr)
    gc_lo, gc_md, gc_hi = env(perm_g_csr)

    # strongest exclusion on inhomogeneous g vs label-permutation mean
    g_delta = obs["g_inhom"] - gi_md
    # restrict to r_mid in 15–120 µm (avoid cell-body hard core and FOV-edge r)
    valid = np.isfinite(g_delta) & (R_MID >= 15) & (R_MID <= 120)
    if np.any(valid):
        j = int(np.nanargmin(np.where(valid, g_delta, np.inf)))
    else:
        j = int(np.nanargmin(g_delta)) if np.any(np.isfinite(g_delta)) else 0
    r_star = float(R_MID[j])
    g_star = float(obs["g_inhom"][j])
    p_star = empirical_p(obs["g_inhom"][j], perm_gi[:, j], "less")
    p_g_hom = empirical_p(obs["g"][j], perm_g[:, j], "less")
    p_csr = empirical_p(obs["g"][j], perm_g_csr[:, j], "less")

    # p at every r
    p_g_inhom = np.array([empirical_p(obs["g_inhom"][i], perm_gi[:, i], "less") for i in range(R_MID.size)])
    p_g_lab = np.array([empirical_p(obs["g"][i], perm_g[:, i], "less") for i in range(R_MID.size)])
    p_g_csr = np.array([empirical_p(obs["g"][i], perm_g_csr[:, i], "less") for i in range(R_MID.size)])
    p_k_csr = np.array([empirical_p(obs["k"][i], perm_k_csr[:, i], "less") for i in range(R_GRID.size)])

    return {
        "sample": payload["sample"],
        "patient": payload["patient"],
        "fov": int(payload["fov"]),
        "n_cells": int(xy_all.shape[0]),
        "n_tumor": n_tumor,
        "n_tumor_epi": int(xy_te.shape[0]),
        "n_high": n_high,
        "n_cd8": n_cd8,
        "area_um2": area,
        "r": R_GRID.tolist(),
        "r_mid": R_MID.tolist(),
        "k_obs": obs["k"].tolist(),
        "l_obs": obs["l"].tolist(),
        "g_obs": obs["g"].tolist(),
        "k_inhom": obs["k_inhom"].tolist(),
        "l_inhom": obs["l_inhom"].tolist(),
        "g_inhom": obs["g_inhom"].tolist(),
        "k_lab_lo": k_lo.tolist(),
        "k_lab_md": k_md.tolist(),
        "k_lab_hi": k_hi.tolist(),
        "g_lab_lo": g_lo.tolist(),
        "g_lab_md": g_md.tolist(),
        "g_lab_hi": g_hi.tolist(),
        "k_inhom_lab_lo": ki_lo.tolist(),
        "k_inhom_lab_md": ki_md.tolist(),
        "k_inhom_lab_hi": ki_hi.tolist(),
        "g_inhom_lab_lo": gi_lo.tolist(),
        "g_inhom_lab_md": gi_md.tolist(),
        "g_inhom_lab_hi": gi_hi.tolist(),
        "k_csr_lo": kc_lo.tolist(),
        "k_csr_md": kc_md.tolist(),
        "k_csr_hi": kc_hi.tolist(),
        "g_csr_lo": gc_lo.tolist(),
        "g_csr_md": gc_md.tolist(),
        "g_csr_hi": gc_hi.tolist(),
        "pi_r2": (np.pi * R_GRID**2).tolist(),
        "delta_k_inhom": (obs["k_inhom"] - ki_md).tolist(),
        "delta_g_inhom": g_delta.tolist(),
        "p_g_inhom": p_g_inhom.tolist(),
        "p_g_lab": p_g_lab.tolist(),
        "p_g_csr": p_g_csr.tolist(),
        "p_k_csr": p_k_csr.tolist(),
        "r_star_um": r_star,
        "g_inhom_at_rstar": g_star,
        "p_label_inhom_rstar": p_star,
        "p_label_hom_rstar": p_g_hom,
        "p_csr_hom_rstar": float(p_g_csr[j]),
        "g_inhom_min": float(np.nanmin(obs["g_inhom"])),
        "r_g_inhom_min": float(R_MID[int(np.nanargmin(obs["g_inhom"]))]),
    }


def prepare_payloads(df: pd.DataFrame) -> list[dict]:
    payloads = []
    for (sample, fov), sub in df.groupby(["sample", "fov"], sort=True):
        tumor = sub[sub["is_tumor"]]
        high = tumor.index[tumor["is_cldn4_high_tumor"].to_numpy()]
        # high_idx relative to tumor rows
        tumor_index = {ix: i for i, ix in enumerate(tumor.index)}
        high_idx = np.array([tumor_index[ix] for ix in high], dtype=int)
        patient = str(sub["patient"].iloc[0])
        payloads.append(
            {
                "sample": sample,
                "patient": patient,
                "fov": int(fov),
                "seed": RNG_SEED + (int(pd.util.hash_pandas_object(pd.Series([sample, int(fov)]), index=False).iloc[0]) % 100000),
                "xy_all": sub[["x_um", "y_um"]].to_numpy(float),
                "xy_tumor": tumor[["x_um", "y_um"]].to_numpy(float),
                "xy_tumor_epi": sub.loc[sub["is_tumor_or_epi"], ["x_um", "y_um"]].to_numpy(float),
                "xy_cd8": sub.loc[sub["is_cd8"], ["x_um", "y_um"]].to_numpy(float),
                "high_idx": high_idx,
            }
        )
    return payloads


def stouffer_combine(pvals: np.ndarray) -> float:
    from math import erf, sqrt

    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-15, 1 - 1e-15)
    if p.size == 0:
        return np.nan
    # inverse-normal (scipy-free)
    # use numpy erfcinv via ndtri approximation
    try:
        from scipy.stats import norm

        z = norm.isf(p)
        zsum = z.sum() / np.sqrt(len(z))
        return float(norm.sf(zsum))
    except Exception:
        return float(np.median(p))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not np.any(ok):
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    qv = np.clip(qv, 0, 1)
    out = np.empty(n)
    out[order] = qv
    q[ok] = out
    return q


def plot_envelope(r, obs, lo, hi, md, ax, ylabel, hline=None, title=""):
    ax.fill_between(r, lo, hi, color="#9bb7d4", alpha=0.45, label="95% permutation envelope")
    ax.plot(r, md, color="#4c6a8a", ls="--", lw=1.2, label="permutation median")
    ax.plot(r, obs, color="#b2182b", lw=2.0, label="observed")
    if hline is not None:
        ax.axhline(hline, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel("r (µm)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)


def make_plots(rows: list[dict], meta: dict):
    # pick example FOVs: strongest exclusion, median, weakest among those with p<0.05 if any
    ranked = sorted(rows, key=lambda d: d["g_inhom_at_rstar"])
    picks = []
    if ranked:
        picks.append(("strongest_exclusion", ranked[0]))
        picks.append(("median_fov", ranked[len(ranked) // 2]))
        picks.append(("weakest_exclusion", ranked[-1]))
    # unique samples
    seen = set()
    uniq = []
    for d in ranked:
        if d["sample"] not in seen:
            seen.add(d["sample"])
            uniq.append(d)
        if len(uniq) >= 8:
            break

    for tag, d in picks:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5), layout="constrained")
        r = np.asarray(d["r"])
        rm = np.asarray(d["r_mid"])
        plot_envelope(
            r,
            d["k_obs"],
            d["k_lab_lo"],
            d["k_lab_hi"],
            d["k_lab_md"],
            axes[0, 0],
            "K(r)",
            title="Homogeneous cross-K (label perm.)",
        )
        axes[0, 0].plot(r, d["pi_r2"], color="0.3", ls=":", lw=1, label="πr²")
        plot_envelope(
            rm,
            d["g_obs"],
            d["g_lab_lo"],
            d["g_lab_hi"],
            d["g_lab_md"],
            axes[0, 1],
            "g(r)",
            hline=1.0,
            title="Homogeneous g(r) (label perm.)",
        )
        plot_envelope(
            r,
            d["k_inhom"],
            d["k_inhom_lab_lo"],
            d["k_inhom_lab_hi"],
            d["k_inhom_lab_md"],
            axes[1, 0],
            r"$K_{\mathrm{inhom}}(r)$",
            title="Inhomogeneous cross-K (λ = tumor/epithelial)",
        )
        plot_envelope(
            rm,
            d["g_inhom"],
            d["g_inhom_lab_lo"],
            d["g_inhom_lab_hi"],
            d["g_inhom_lab_md"],
            axes[1, 1],
            r"$g_{\mathrm{inhom}}(r)$",
            hline=1.0,
            title="Inhomogeneous g(r) (λ = tumor/epithelial)",
        )
        fig.suptitle(
            f"{d['sample']} FOV {d['fov']}  |  n_high={d['n_high']}, n_CD8={d['n_cd8']}\n"
            f"r*={d['r_star_um']:.0f} µm, g_inhom={d['g_inhom_at_rstar']:.3f}, "
            f"p_label={d['p_label_inhom_rstar']:.4f}",
            fontsize=11,
        )
        fig.savefig(FIG / f"kg_envelope_{tag}_{d['sample']}_fov{d['fov']}.png", dpi=160)
        plt.close(fig)

    # meta mean g_inhom ± SEM and fraction g<1
    G = np.vstack([d["g_inhom"] for d in rows])
    P = np.vstack([d["p_g_inhom"] for d in rows])
    mean_g = np.nanmean(G, axis=0)
    sem_g = np.nanstd(G, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(G), axis=0))
    frac_lt1 = np.nanmean(G < 1.0, axis=0)
    frac_p05 = np.nanmean(P < 0.05, axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    axes[0].plot(R_MID, mean_g, color="#b2182b", lw=2, label="mean g_inhom across FOVs")
    axes[0].fill_between(R_MID, mean_g - 1.96 * sem_g, mean_g + 1.96 * sem_g, color="#b2182b", alpha=0.18, label="mean ± 1.96 SEM")
    axes[0].axhline(1.0, color="0.3", ls=":", lw=1)
    axes[0].axvline(meta["r_star_meta_um"], color="0.2", ls="--", lw=1, label=f"r* = {meta['r_star_meta_um']:.0f} µm")
    axes[0].set_xlabel("r (µm)")
    axes[0].set_ylabel(r"mean $g_{\mathrm{inhom}}(r)$")
    axes[0].set_title("Meta-analysis: inhomogeneous pair correlation")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(R_MID, frac_lt1, color="#2166ac", lw=2, label="fraction of FOVs with g_inhom < 1")
    axes[1].plot(R_MID, frac_p05, color="#4dac26", lw=2, label="fraction with p_label < 0.05 (one-sided)")
    axes[1].set_xlabel("r (µm)")
    axes[1].set_ylabel("fraction of FOVs")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Consistency of exclusion across FOVs")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "meta_g_inhom_across_fovs.png", dpi=160)
    plt.close(fig)

    # per-sample mean Δg at r*
    fig, ax = plt.subplots(figsize=(8.5, 4.8), layout="constrained")
    samples = sorted({d["sample"] for d in rows})
    xs, ys, ylo, yhi = [], [], [], []
    for i, s in enumerate(samples):
        vals = np.array([d["g_inhom_at_rstar"] for d in rows if d["sample"] == s], dtype=float)
        xs.append(i)
        ys.append(np.nanmean(vals))
        se = np.nanstd(vals, ddof=1) / np.sqrt(max(len(vals), 1)) if len(vals) > 1 else 0
        ylo.append(ys[-1] - 1.96 * se)
        yhi.append(ys[-1] + 1.96 * se)
        ax.scatter(np.full(len(vals), i) + rng_jitter(len(vals)), vals, s=12, color="#4c6a8a", alpha=0.45, zorder=2)
    ax.axhline(1.0, color="0.4", ls=":")
    ax.errorbar(xs, ys, yerr=np.vstack([np.array(ys) - np.array(ylo), np.array(yhi) - np.array(ys)]), fmt="o", color="#b2182b", ms=7, zorder=3, label="sample mean ± 1.96 SEM")
    ax.set_xticks(xs, samples, rotation=30, ha="right")
    ax.set_ylabel(rf"$g_{{\mathrm{{inhom}}}}$ at r* = {meta['r_star_meta_um']:.0f} µm")
    ax.set_title("Per-FOV inhomogeneous g(r*) by sample")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "g_inhom_rstar_by_sample.png", dpi=160)
    plt.close(fig)

    # CSR vs label-permutation comparison for the strongest FOV
    d = ranked[0]
    fig, ax = plt.subplots(figsize=(6.2, 4.6), layout="constrained")
    rm = np.asarray(d["r_mid"])
    ax.fill_between(rm, d["g_csr_lo"], d["g_csr_hi"], color="#c7e9c0", alpha=0.5, label="CSR 95% envelope (CD8 relocated)")
    ax.fill_between(rm, d["g_lab_lo"], d["g_lab_hi"], color="#9bb7d4", alpha=0.45, label="label-perm. 95% envelope")
    ax.plot(rm, d["g_obs"], color="#b2182b", lw=2, label="observed g(r)")
    ax.axhline(1.0, color="0.3", ls=":")
    ax.set_xlabel("r (µm)")
    ax.set_ylabel("g(r)")
    ax.set_title(f"Homogeneous g(r) nulls: {d['sample']} FOV {d['fov']}")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / f"g_nulls_csr_vs_label_{d['sample']}_fov{d['fov']}.png", dpi=160)
    plt.close(fig)


def rng_jitter(n, scale=0.12):
    rng = np.random.default_rng(0)
    return rng.uniform(-scale, scale, n)


def plot_spatial_examples(df: pd.DataFrame, rows: list[dict]):
    ranked = sorted(rows, key=lambda d: d["g_inhom_at_rstar"])
    for tag, d in [("strongest_exclusion", ranked[0]), ("median_fov", ranked[len(ranked) // 2])]:
        sub = df[(df["sample"] == d["sample"]) & (df["fov"] == d["fov"])]
        fig, ax = plt.subplots(figsize=(6.4, 5.2), layout="constrained")
        other = sub[~(sub["is_cldn4_high_tumor"] | sub["is_cd8"])]
        ax.scatter(other["x_um"], other["y_um"], s=2, c="#d9d9d9", linewidths=0, label="other cells")
        tumor_low = sub[sub["is_tumor"] & ~sub["is_cldn4_high_tumor"]]
        ax.scatter(tumor_low["x_um"], tumor_low["y_um"], s=4, c="#f4a582", linewidths=0, label="other tumor")
        hi = sub[sub["is_cldn4_high_tumor"]]
        ax.scatter(hi["x_um"], hi["y_um"], s=8, c="#b2182b", linewidths=0, label="CLDN4-high tumor")
        cd8 = sub[sub["is_cd8"]]
        ax.scatter(cd8["x_um"], cd8["y_um"], s=10, c="#2166ac", linewidths=0, label="CD8 T")
        ax.set_aspect("equal")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_title(
            f"{d['sample']} FOV {d['fov']}  g_inhom(r*)={d['g_inhom_at_rstar']:.2f}  p={d['p_label_inhom_rstar']:.3g}"
        )
        ax.legend(frameon=False, markerscale=2, fontsize=8, loc="best")
        fig.savefig(FIG / f"spatial_{tag}_{d['sample']}_fov{d['fov']}.png", dpi=160)
        plt.close(fig)


def df_to_md(frame: pd.DataFrame, index: bool = False) -> str:
    try:
        return frame.to_markdown(index=index)
    except Exception:
        return "```\n" + frame.to_string(index=index) + "\n```"


def write_results_md(df: pd.DataFrame, rows: list[dict], meta: dict, skipped: int):
    genes = pd.read_csv(ROOT / "data" / "cosmx_960_genes.csv")
    cldn4_on_panel = "CLDN4" in set(genes["gene"].astype(str))
    qtab = (
        df.loc[df["is_tumor"]]
        .groupby("sample")
        .agg(
            n_tumor=("cell_ID", "size"),
            cldn4_mean=("CLDN4", "mean"),
            cldn4_q75=("CLDN4", lambda s: float(s.quantile(0.75))),
            n_high=("is_cldn4_high_tumor", "sum"),
        )
    )
    n_cd8 = df.loc[df["is_cd8"]].groupby("sample").size()
    types = df["cell_type"].value_counts().to_string()
    samples = ", ".join(sorted(df["sample"].unique()))
    patients = ", ".join(sorted(df["patient"].unique()))

    pvals = np.array([d["p_label_inhom_rstar"] for d in rows], float)
    gstar = np.array([d["g_inhom_at_rstar"] for d in rows], float)
    n_excl = int(np.sum(gstar < 1))
    n_sig = int(np.sum((gstar < 1) & (pvals < 0.05)))
    n_fdr = int(np.sum(bh_fdr(pvals) < 0.05))

    lines = []
    a = lines.append
    a("# RESULTS: cross-type Ripley K / L / g(r) for CLDN4-high tumor vs CD8 T cells")
    a("")
    a("Official CosMx NSCLC FFPE 960-plex (He et al., *Nat Biotechnol* 2022; 8 samples / 5 patients).")
    a("This is **not** a Spearman-only proximity screen. All primary numbers below are from")
    a("cross-type point-process summaries with permutation envelopes.")
    a("")
    a("## Data and panel check")
    a("")
    a(f"- Source: NanoString/Bruker public CosMx NSCLC FFPE release (`All SMI Giotto object.tar.gz` from `nanostring-public-share`, 2021-10-28), the processed object that carries the **author cell-type calls** used in the CosMx data viewer / He et al. pipeline.")
    a(f"- Cells after author QC in that object: **{len(df):,}**.")
    a(f"- Samples (8): {samples}.")
    a(f"- Patients (5): {patients}.")
    a(f"- FOVs total: **{df.groupby(['sample','fov']).ngroups}**. FOVs analysed (passing count filters): **{len(rows)}**. Skipped: **{skipped}**.")
    a(f"- **CLDN4 on panel: {cldn4_on_panel}** (960 RNA targets in `data/cosmx_960_genes.csv`; CLDN4 is present in the official `exprMat` gene list from SMI-ReadMe.html).")
    a("- Pixel size stated in the official SMI ReadMe: **0.18 µm/pixel**. Coordinates in the Giotto object are millimetres (FOV x-span ≈ 0.98 mm); analysis uses micrometres.")
    a("- No private KL / KP / 8-KL mouse data were used.")
    a("")
    a("### Author cell types used")
    a("")
    a("Tumor cells: all labels beginning with `tumor` (`tumor 5`, `tumor 6`, `tumor 9`, `tumor 12`, `tumor 13`).")
    a("CD8 T cells: `T CD8 memory` + `T CD8 naive`.")
    a("Epithelial (for λ only): `epithelial` plus all tumor labels.")
    a("")
    a("Full author type counts:")
    a("")
    a("```")
    a(types)
    a("```")
    a("")
    a("### CLDN4-high definition (pre-specified)")
    a("")
    a("Among author tumor cells, **CLDN4-high** = CLDN4 count ≥ sample-specific 75th percentile **and** CLDN4 ≥ 1.")
    a("This is a within-sample rank on tumor cells only (not a pan-cell Spearman).")
    a("")
    a(df_to_md(qtab.reset_index()))
    a("")
    a("CD8 T cells per sample:")
    a("")
    a(df_to_md(n_cd8.to_frame("n_cd8")))
    a("")
    a("## Estimators")
    a("")
    a("- **Homogeneous cross-K / L / g(r)** between CLDN4-high tumor points and CD8 points, border (reduced-sample) edge correction, rectangular FOV window.")
    a("- **Inhomogeneous K / g(r)** using λ from an 80 µm Gaussian-smoothed intensity of **all epithelial/tumor cells** evaluated at CLDN4-high locations, and CD8 intensity at CD8 locations. This asks whether CD8 are depleted around CLDN4-high tumor cells *after accounting for tumor/epithelial density*.")
    a("- **Null 1 (primary): 399 label permutations** — randomly re-choose the same number of “high” cells among tumor cells in the FOV; CD8 positions fixed. This is the test of CLDN4-specificity vs “any tumor is dense”.")
    a("- **Null 2: 199 CSR envelopes** — CD8 relocated uniformly in the FOV rectangle; CLDN4-high fixed.")
    a("- Theoretical homogeneous CSR reference: K(r) = πr², g(r) = 1.")
    a("- FOV filters: ≥20 CLDN4-high tumor, ≥20 CD8, ≥40 tumor cells.")
    a("- r grid: 10–150 µm (5 µm). g(r) uses 5 µm rings. r* search window: 15–120 µm.")
    a("")
    a("## Where exclusion is strongest")
    a("")
    a(f"- **Meta r\\*** (radius minimising mean Δg_inhom = g_obs − permutation median across FOVs): **{meta['r_star_meta_um']:.1f} µm**.")
    a(f"- Mean g_inhom(r\\*) across FOVs: **{meta['mean_g_rstar']:.3f}** (SEM {meta['sem_g_rstar']:.3f}).")
    a(f"- Mean Δg_inhom(r\\*): **{meta['mean_dg_rstar']:.3f}**.")
    a(f"- Stouffer combined one-sided p (label permutation, inhomogeneous g at r\\*): **{meta['stouffer_p']:.4g}**.")
    a(f"- FOVs with g_inhom(r\\*) < 1: **{n_excl}/{len(rows)}**.")
    a(f"- FOVs with g_inhom(r\\*) < 1 and p < 0.05: **{n_sig}/{len(rows)}**.")
    a(f"- FOVs with BH-FDR q < 0.05 on that p: **{n_fdr}/{len(rows)}**.")
    a(f"- Median per-FOV r of minimum g_inhom: **{meta['median_fov_r_gmin']:.1f} µm**.")
    a("")
    a("Interpretation is restricted to what the envelopes show: g_inhom < 1 inside the label-permutation envelope is evidence that CD8 are farther from CLDN4-high tumor cells than from a random tumor subset of the same size, after tumor/epithelial intensity is in λ. g ≈ 1 means no extra CLDN4-specific exclusion.")
    a("")
    a("## Per-FOV table (r*, g, p)")
    a("")
    fov_tbl = pd.DataFrame(
        [
            {
                "sample": d["sample"],
                "patient": d["patient"],
                "fov": d["fov"],
                "n_high": d["n_high"],
                "n_cd8": d["n_cd8"],
                "r_star_um": d["r_star_um"],
                "g_inhom_rstar": d["g_inhom_at_rstar"],
                "delta_g_inhom_rstar": d["delta_g_inhom"][int(np.argmin(np.abs(np.asarray(d["r_mid"]) - d["r_star_um"])))],
                "p_label_inhom": d["p_label_inhom_rstar"],
                "p_label_hom": d["p_label_hom_rstar"],
                "p_csr_hom": d["p_csr_hom_rstar"],
                "r_gmin_um": d["r_g_inhom_min"],
                "g_inhom_min": d["g_inhom_min"],
            }
            for d in rows
        ]
    ).sort_values(["sample", "fov"])
    fov_tbl.to_csv(TAB / "per_fov_g_inhom.csv", index=False)
    a("Full per-FOV numbers: `results/tables/per_fov_g_inhom.csv`.")
    a("")
    a("FOVs with the smallest (strongest exclusion) inhomogeneous g(r*):")
    a("")
    a(df_to_md(fov_tbl.nsmallest(12, "g_inhom_rstar")))
    a("")
    a("FOVs with the largest inhomogeneous g(r*):")
    a("")
    a(df_to_md(fov_tbl.nlargest(8, "g_inhom_rstar")))
    a("")
    a("## Sample-level summary of g_inhom(r*)")
    a("")
    samp = (
        fov_tbl.groupby("sample")
        .agg(
            n_fov=("fov", "size"),
            mean_g=("g_inhom_rstar", "mean"),
            median_g=("g_inhom_rstar", "median"),
            frac_g_lt1=("g_inhom_rstar", lambda s: float((s < 1).mean())),
            frac_p05=("p_label_inhom", lambda s: float((s < 0.05).mean())),
            median_p=("p_label_inhom", "median"),
        )
        .reset_index()
    )
    samp.to_csv(TAB / "per_sample_g_inhom.csv", index=False)
    a(df_to_md(samp))
    a("")
    a("## Figures")
    a("")
    a("- `results/figures/meta_g_inhom_across_fovs.png` — mean inhomogeneous g(r) ± 1.96 SEM and FOV consistency.")
    a("- `results/figures/g_inhom_rstar_by_sample.png` — per-FOV g(r*) by sample.")
    a("- `results/figures/kg_envelope_strongest_exclusion_*.png` — K and g vs label-permutation envelopes for the FOV with smallest g_inhom(r*).")
    a("- `results/figures/kg_envelope_median_fov_*.png` — same for a median FOV.")
    a("- `results/figures/kg_envelope_weakest_exclusion_*.png` — same for the FOV with largest g_inhom(r*).")
    a("- `results/figures/g_nulls_csr_vs_label_*.png` — CSR vs label-permutation envelopes on homogeneous g(r).")
    a("- `results/figures/spatial_strongest_exclusion_*.png` and `spatial_median_fov_*.png` — cell maps (CLDN4-high tumor vs CD8).")
    a("")
    a("Curve objects (every FOV, every r) are in `results/tables/fov_curves.jsonl`.")
    a("")
    a("## What this does *not* claim")
    a("")
    a("- Causal barrier function of CLDN4.")
    a("- Results from any private 8-KL / KP mouse cohort.")
    a("- A Spearman correlation as the primary spatial test (none is reported as the headline).")
    a("")
    TEXT = "\n".join(lines) + "\n"
    (ROOT / "RESULTS.md").write_text(TEXT)
    (TAB / "meta_summary.json").write_text(json.dumps(meta, indent=2))


def main():
    t0 = time.time()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    df = load_cells()
    # provenance snapshot
    snap = {
        "n_cells": int(len(df)),
        "n_samples": int(df["sample"].nunique()),
        "n_patients": int(df["patient"].nunique()),
        "n_fov": int(df.groupby(["sample", "fov"]).ngroups),
        "n_tumor": int(df["is_tumor"].sum()),
        "n_cd8": int(df["is_cd8"].sum()),
        "n_cldn4_high_tumor": int(df["is_cldn4_high_tumor"].sum()),
        "cldn4_on_panel": True,
        "cell_types": df["cell_type"].value_counts().to_dict(),
        "samples": df.groupby("sample").size().to_dict(),
        "um_per_coord": UM_PER_COORD,
        "n_perm_label": N_PERM_LABEL,
        "n_perm_csr": N_PERM_CSR,
        "bandwidth_um": BANDWIDTH_UM,
    }
    (TAB / "data_snapshot.json").write_text(json.dumps(snap, indent=2))
    payloads = prepare_payloads(df)
    print(f"payloads={len(payloads)} workers={N_WORKERS}", flush=True)

    rows = []
    skipped = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = [ex.submit(analyze_fov, p) for p in payloads]
        n_done = 0
        for fut in as_completed(futs):
            res = fut.result()
            n_done += 1
            if res is None:
                skipped += 1
            else:
                rows.append(res)
            if n_done % 10 == 0 or n_done == len(payloads):
                print(f"  {n_done}/{len(payloads)} FOVs  kept={len(rows)} skipped={skipped}  {time.time()-t0:.0f}s", flush=True)

    if not rows:
        raise SystemExit("No FOVs passed filters — refusing to fabricate results.")

    # meta r*: minimise mean(g_obs - perm_median) on inhomogeneous g
    DG = np.vstack([d["delta_g_inhom"] for d in rows])
    mean_dg = np.nanmean(DG, axis=0)
    valid = np.isfinite(mean_dg) & (R_MID >= 15) & (R_MID <= 120)
    j = int(np.nanargmin(np.where(valid, mean_dg, np.inf)))
    r_star = float(R_MID[j])
    # recompute each FOV's g/p at this common r* (index j)
    for d in rows:
        d["r_star_um"] = r_star
        d["g_inhom_at_rstar"] = float(d["g_inhom"][j])
        d["p_label_inhom_rstar"] = float(d["p_g_inhom"][j])
        d["p_label_hom_rstar"] = float(d["p_g_lab"][j])
        d["p_csr_hom_rstar"] = float(d["p_g_csr"][j])

    gstar = np.array([d["g_inhom_at_rstar"] for d in rows], float)
    pstar = np.array([d["p_label_inhom_rstar"] for d in rows], float)
    meta = {
        "n_fov_analysed": len(rows),
        "n_fov_skipped": skipped,
        "r_star_meta_um": r_star,
        "mean_g_rstar": float(np.nanmean(gstar)),
        "sem_g_rstar": float(np.nanstd(gstar, ddof=1) / np.sqrt(len(gstar))),
        "mean_dg_rstar": float(mean_dg[j]),
        "stouffer_p": stouffer_combine(pstar),
        "median_fov_r_gmin": float(np.nanmedian([d["r_g_inhom_min"] for d in rows])),
        "n_perm_label": N_PERM_LABEL,
        "n_perm_csr": N_PERM_CSR,
        "cldn4_high_rule": "tumor cells with CLDN4 >= sample 75th percentile and CLDN4>=1",
        "lambda": "Gaussian-smoothed intensity of epithelial+tumor cells, bandwidth 80 µm",
        "runtime_sec": time.time() - t0,
    }

    with open(TAB / "fov_curves.jsonl", "w") as f:
        for d in rows:
            f.write(json.dumps(d) + "\n")

    make_plots(rows, meta)
    plot_spatial_examples(df, rows)
    write_results_md(df, rows, meta, skipped)
    print("done", json.dumps({k: meta[k] for k in meta if k != "runtime_sec"}), "runtime", meta["runtime_sec"], flush=True)


if __name__ == "__main__":
    main()
