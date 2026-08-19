#!/usr/bin/env python3
"""CLDN4-only CosMx NSCLC metrics on Lung5_Rep2 (He et al. 2022 public set)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path("/workspace/data/cosmx/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images")
OUT = Path("/workspace/results")
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# CosMx SMI pixel size used by NanoString / Giotto loaders for this release
UM_PER_PX = 0.18

GENES = ["CLDN4", "CD8A", "CD3D", "CD3E", "NKG7", "KRT8", "KRT19", "EPCAM", "PTPRC"]


def spr(a, b, mask=None):
    if mask is None:
        mask = np.ones(len(a), dtype=bool)
    aa = np.asarray(a)[mask]
    bb = np.asarray(b)[mask]
    ok = np.isfinite(aa) & np.isfinite(bb)
    aa, bb = aa[ok], bb[ok]
    if len(aa) < 10:
        return {"rho": np.nan, "p": np.nan, "n": int(len(aa))}
    r, p = spearmanr(aa, bb)
    return {"rho": float(r), "p": float(p), "n": int(len(aa))}


def residualize(y, x):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(y) & np.isfinite(x)
    out = np.full_like(y, np.nan, dtype=float)
    if mask.sum() < 10:
        return out
    A = np.column_stack([np.ones(mask.sum()), x[mask]])
    coef, _, _, _ = np.linalg.lstsq(A, y[mask], rcond=None)
    out[mask] = y[mask] - (coef[0] + coef[1] * x[mask])
    return out


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def main():
    usecols = ["fov", "cell_ID"] + GENES
    expr = pd.read_csv(ROOT / "Lung5_Rep2_exprMat_file.csv", usecols=usecols)
    meta = pd.read_csv(
        ROOT / "Lung5_Rep2_metadata_file.csv",
        usecols=[
            "fov",
            "cell_ID",
            "CenterX_global_px",
            "CenterY_global_px",
            "Area",
            "Mean.PanCK",
            "Mean.CD45",
            "Mean.CD3",
        ],
    )
    df = meta.merge(expr, on=["fov", "cell_ID"], how="inner")
    gene_sum = df[GENES].sum(axis=1)
    # keep cells with at least a few of the tracked transcripts or any area
    keep = (df["Area"] > 0) & (gene_sum >= 1)
    df = df.loc[keep].copy()

    lib_proxy = df[GENES].sum(axis=1).clip(lower=1)
    # panel is 960-plex; use raw counts for Spearman (standard for targeted CosMx)
    # also log1p of per-cell total of selected genes as a sensitivity check
    for g in GENES:
        df[f"{g}_log"] = np.log1p(df[g].to_numpy(float))

    xy_um = df[["CenterX_global_px", "CenterY_global_px"]].to_numpy(float) * UM_PER_PX

    panck = df["Mean.PanCK"].to_numpy(float)
    cd3_prot = df["Mean.CD3"].to_numpy(float)
    # epithelial / tumor: PanCK above median among cells with PanCK>0, else top 40%
    pos = panck > 0
    panck_cut = np.quantile(panck[pos], 0.60) if pos.sum() else np.median(panck)
    is_epi = panck >= panck_cut

    cldn4 = df["CLDN4_log"].to_numpy(float)
    cd8a = df["CD8A_log"].to_numpy(float)
    cd3e = df["CD3E_log"].to_numpy(float)
    nkg7 = df["NKG7_log"].to_numpy(float)
    krt8 = df["KRT8_log"].to_numpy(float)

    cldn4_q75 = np.quantile(cldn4[is_epi], 0.75)
    cldn4_q25 = np.quantile(cldn4[is_epi], 0.25)
    is_hi = is_epi & (cldn4 >= cldn4_q75)
    is_lo = is_epi & (cldn4 <= cldn4_q25)

    # CD8 cells: CD8A count >= 1 (RNA). CosMx is sparse; count>=1 is the usual call.
    is_cd8 = df["CD8A"].to_numpy(float) >= 1

    same_all = {
        "CLDN4_vs_CD8A": spr(cldn4, cd8a),
        "CLDN4_vs_CD3E": spr(cldn4, cd3e),
        "CLDN4_vs_NKG7": spr(cldn4, nkg7),
    }
    same_epi = {
        "CLDN4_vs_CD8A": spr(cldn4, cd8a, is_epi),
        "CLDN4_vs_CD3E": spr(cldn4, cd3e, is_epi),
        "CLDN4_vs_NKG7": spr(cldn4, nkg7, is_epi),
    }
    resid = residualize(cldn4, krt8)
    resid_all = spr(resid, cd8a)
    resid_epi = spr(resid, cd8a, is_epi)

    tree = cKDTree(xy_um[is_cd8]) if is_cd8.sum() else None

    def nn(mask):
        if tree is None or mask.sum() == 0:
            return np.array([])
        d, _ = tree.query(xy_um[mask], k=1)
        return np.asarray(d, dtype=float)

    d_hi, d_lo = nn(is_hi), nn(is_lo)
    if len(d_hi) and len(d_lo):
        _, nn_p = mannwhitneyu(d_hi, d_lo, alternative="two-sided")
    else:
        nn_p = np.nan

    radii = [25, 50, 100]
    nhood = {}
    nhood_spr = {}
    if tree is not None:
        for r in radii:
            c_hi = np.asarray(tree.query_ball_point(xy_um[is_hi], r=r, return_length=True), float)
            c_lo = np.asarray(tree.query_ball_point(xy_um[is_lo], r=r, return_length=True), float)
            c_all = np.asarray(tree.query_ball_point(xy_um, r=r, return_length=True), float)
            if len(c_hi) and len(c_lo):
                _, p = mannwhitneyu(c_hi, c_lo, alternative="two-sided")
            else:
                p = np.nan
            nhood[str(r)] = {
                "radius_um": r,
                "cldn4_hi_mean": float(np.mean(c_hi)) if len(c_hi) else np.nan,
                "cldn4_hi_median": float(np.median(c_hi)) if len(c_hi) else np.nan,
                "cldn4_lo_mean": float(np.mean(c_lo)) if len(c_lo) else np.nan,
                "cldn4_lo_median": float(np.median(c_lo)) if len(c_lo) else np.nan,
                "mwu_p": float(p) if p == p else np.nan,
                "n_hi": int(len(c_hi)),
                "n_lo": int(len(c_lo)),
            }
            nhood_spr[str(r)] = {
                "all_cells": spr(cldn4, c_all),
                "epithelial": spr(cldn4, c_all, is_epi),
            }

    # per-FOV Spearman as a robustness slice
    fov_rows = []
    for fov, sub in df.groupby("fov"):
        if len(sub) < 80:
            continue
        r, p = spearmanr(np.log1p(sub["CLDN4"]), np.log1p(sub["CD8A"]))
        fov_rows.append({"fov": int(fov), "n": int(len(sub)), "rho": float(r), "p": float(p)})
    fov_df = pd.DataFrame(fov_rows)
    fov_df.to_csv(TAB / "cosmx_lung5rep2_per_fov_spearman.csv", index=False)

    summary = {
        "dataset": "NanoString CosMx SMI NSCLC FFPE Lung5_Rep2 (He et al. 2022 prototype 960-plex)",
        "accession": "NanoString-CosMx-NSCLC-Lung5_Rep2",
        "platform": "CosMx SMI 960-plex RNA + PanCK/CD45/CD3 protein",
        "n_cells": int(len(df)),
        "n_fov": int(df["fov"].nunique()),
        "um_per_pixel": UM_PER_PX,
        "genes_present": {g: True for g in GENES},
        "n_epithelial_panck": int(is_epi.sum()),
        "n_cldn4_high_epi": int(is_hi.sum()),
        "n_cldn4_low_epi": int(is_lo.sum()),
        "n_cd8_cells_CD8A_ge1": int(is_cd8.sum()),
        "panck_cut": float(panck_cut),
        "same_cell_spearman_all": same_all,
        "same_cell_spearman_epithelial": same_epi,
        "cldn4_krt8_residual_vs_CD8A_all": resid_all,
        "cldn4_krt8_residual_vs_CD8A_epithelial": resid_epi,
        "nearest_cd8_um": {
            "cldn4_hi_median": float(np.median(d_hi)) if len(d_hi) else np.nan,
            "cldn4_hi_mean": float(np.mean(d_hi)) if len(d_hi) else np.nan,
            "cldn4_lo_median": float(np.median(d_lo)) if len(d_lo) else np.nan,
            "cldn4_lo_mean": float(np.mean(d_lo)) if len(d_lo) else np.nan,
            "mwu_p": float(nn_p) if nn_p == nn_p else np.nan,
            "n_hi": int(len(d_hi)),
            "n_lo": int(len(d_lo)),
        },
        "nhood_cd8_count": nhood,
        "nhood_spearman_cldn4_vs_cd8count": nhood_spr,
        "per_fov_spearman_CLDN4_CD8A": {
            "n_fov": int(len(fov_df)),
            "median_rho": float(fov_df["rho"].median()) if len(fov_df) else np.nan,
            "mean_rho": float(fov_df["rho"].mean()) if len(fov_df) else np.nan,
            "n_negative": int((fov_df["rho"] < 0).sum()) if len(fov_df) else 0,
        },
        "definition": {
            "epithelial": "Mean.PanCK >= 60th percentile of PanCK-positive cells",
            "cldn4_high": "epithelial AND log1p(CLDN4) >= Q75 within epithelial",
            "cd8_cell": "CD8A RNA count >= 1",
            "not_used": "TACSTD2 intersection was not applied (CLDN4-only)",
        },
    }
    (TAB / "cosmx_lung5rep2_summary.json").write_text(json.dumps(summary, indent=2))

    # figures
    plt.rcParams.update({"font.size": 9, "figure.dpi": 140})

    # subsample for scatter/maps
    rng = np.random.default_rng(0)
    idx = rng.choice(len(df), size=min(40000, len(df)), replace=False)

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for ax, y, title in [(axes[0], cd8a, "CD8A"), (axes[1], cd3e, "CD3E"), (axes[2], nkg7, "NKG7")]:
        ax.scatter(cldn4[idx], y[idx], s=2, alpha=0.15, c="#335588", linewidths=0)
        r = same_all[f"CLDN4_vs_{title}"]
        ax.set_xlabel("CLDN4 (log1p count)")
        ax.set_ylabel(f"{title} (log1p count)")
        ax.set_title(f"ρ={r['rho']:.3f}  p={r['p']:.1e}  n={r['n']}")
    fig.suptitle("CosMx Lung5_Rep2 — same-cell Spearman (CLDN4-only)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "cosmx_lung5rep2_same_cell_spearman.png", bbox_inches="tight")
    fig.savefig(FIG / "cosmx_lung5rep2_same_cell_spearman.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    sc0 = axes[0].scatter(
        xy_um[idx, 0], -xy_um[idx, 1], c=cldn4[idx], s=1.5, cmap="magma", linewidths=0
    )
    axes[0].set_title("CLDN4")
    axes[0].set_aspect("equal")
    axes[0].axis("off")
    fig.colorbar(sc0, ax=axes[0], fraction=0.03, pad=0.02)
    sc1 = axes[1].scatter(
        xy_um[idx, 0], -xy_um[idx, 1], c=cd8a[idx], s=1.5, cmap="viridis", linewidths=0
    )
    axes[1].set_title("CD8A")
    axes[1].set_aspect("equal")
    axes[1].axis("off")
    fig.colorbar(sc1, ax=axes[1], fraction=0.03, pad=0.02)
    fig.suptitle("CosMx Lung5_Rep2 — spatial maps (µm, 40k-cell subsample)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "cosmx_lung5rep2_spatial_maps.png", bbox_inches="tight")
    fig.savefig(FIG / "cosmx_lung5rep2_spatial_maps.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    axes[0].boxplot([d_lo, d_hi], tick_labels=["CLDN4-low epi", "CLDN4-high epi"], showfliers=False)
    axes[0].set_ylabel("Nearest CD8A+ cell (µm)")
    axes[0].set_title(
        f"median {np.median(d_lo):.1f} vs {np.median(d_hi):.1f} µm; p={nn_p:.2e}"
        if len(d_hi) and len(d_lo)
        else "nearest CD8"
    )
    r50 = nhood.get("50", {})
    if tree is not None:
        c_lo = np.asarray(tree.query_ball_point(xy_um[is_lo], r=50, return_length=True), float)
        c_hi = np.asarray(tree.query_ball_point(xy_um[is_hi], r=50, return_length=True), float)
    else:
        c_lo = c_hi = np.array([0.0])
    axes[1].boxplot([c_lo, c_hi], tick_labels=["CLDN4-low epi", "CLDN4-high epi"], showfliers=False)
    axes[1].set_ylabel("CD8A+ cells within 50 µm")
    axes[1].set_title(
        f"mean {r50.get('cldn4_lo_mean', np.nan):.2f} vs {r50.get('cldn4_hi_mean', np.nan):.2f}; "
        f"p={r50.get('mwu_p', np.nan):.2e}"
    )
    fig.suptitle("CosMx Lung5_Rep2 — CLDN4-high vs low epithelial neighborhoods", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "cosmx_lung5rep2_neighborhood.png", bbox_inches="tight")
    fig.savefig(FIG / "cosmx_lung5rep2_neighborhood.pdf", bbox_inches="tight")
    plt.close(fig)

    if len(fov_df):
        fig, ax = plt.subplots(figsize=(6.2, 3.4))
        ax.axhline(0, color="#888", lw=0.8)
        ax.scatter(fov_df["fov"], fov_df["rho"], s=28, c="#224477")
        ax.set_xlabel("FOV")
        ax.set_ylabel("Spearman CLDN4 vs CD8A")
        ax.set_title(
            f"Per-FOV same-cell ρ (median {fov_df['rho'].median():.3f}; "
            f"{int((fov_df['rho']<0).sum())}/{len(fov_df)} negative)"
        )
        fig.tight_layout()
        fig.savefig(FIG / "cosmx_lung5rep2_per_fov_rho.png", bbox_inches="tight")
        fig.savefig(FIG / "cosmx_lung5rep2_per_fov_rho.pdf", bbox_inches="tight")
        plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
