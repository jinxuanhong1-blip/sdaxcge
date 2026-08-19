#!/usr/bin/env python3
"""CLDN4-only spatial immune-exclusion metrics (no TACSTD2 intersection)."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr

OUT = Path("/workspace/results")
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

EPI_GENES = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1", "ELF3", "SFN"]
# CLDN4 deliberately excluded from the epithelial score


def load_visium_h5(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        barcodes = [
            x.decode() if isinstance(x, bytes) else x for x in f["matrix/barcodes"][:]
        ]
        genes = [
            x.decode() if isinstance(x, bytes) else x for x in f["matrix/features/name"][:]
        ]
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(int(x) for x in f["matrix/shape"][:])
    # shape is (n_genes, n_spots)
    X = sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    return pd.Index(barcodes, name="barcode"), pd.Index(genes, name="gene"), X


def log_cp10k(X):
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib = np.maximum(lib, 1.0)
    return X.multiply(1e4 / lib[:, None]).log1p().tocsr()


def gene_vec(X_log, genes, name):
    if name not in genes:
        return None
    return np.asarray(X_log[:, genes.get_loc(name)].todense()).ravel()


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def residualize(y, x):
    """Simple OLS residual of y ~ x (with intercept)."""
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


def analyze_visium_lusc():
    root = Path("/workspace/data/visium_lusc")
    barcodes, genes, X = load_visium_h5(root / "filtered_feature_bc_matrix.h5")
    pos = pd.read_csv(root / "spatial/tissue_positions.csv")
    pos = pos.set_index("barcode").reindex(barcodes)
    scales = json.loads((root / "spatial/scalefactors_json.json").read_text())
    um_per_px = 55.0 / float(scales["spot_diameter_fullres"])

    in_tissue = pos["in_tissue"].fillna(0).astype(int).to_numpy() == 1
    lib = np.asarray(X.sum(axis=1)).ravel()
    keep = in_tissue & (lib >= 250)
    X = X[keep]
    barcodes = barcodes[keep]
    pos = pos.loc[barcodes]
    lib = lib[keep]
    X_log = log_cp10k(X)

    present = {
        g: (g in genes)
        for g in [
            "CLDN4",
            "CD8A",
            "CD8B",
            "CD3D",
            "CD3E",
            "NKG7",
            "KRT8",
            "EPCAM",
            "PTPRC",
        ]
    }

    cldn4 = gene_vec(X_log, genes, "CLDN4")
    cd8a = gene_vec(X_log, genes, "CD8A")
    cd3e = gene_vec(X_log, genes, "CD3E")
    nkg7 = gene_vec(X_log, genes, "NKG7")
    krt8 = gene_vec(X_log, genes, "KRT8")

    xy = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    xy_um = xy * um_per_px

    epi_cols = [g for g in EPI_GENES if g in genes]
    epi_score = np.zeros(X_log.shape[0])
    for g in epi_cols:
        epi_score += zscore(gene_vec(X_log, genes, g))
    epi_score /= max(len(epi_cols), 1)

    # epithelial spots: top 50% epithelial score
    epi_cut = np.median(epi_score)
    is_epi = epi_score >= epi_cut

    # CLDN4-high among epithelial: top quartile of CLDN4 within epi
    cldn4_epi_q75 = np.quantile(cldn4[is_epi], 0.75)
    is_cldn4_hi = is_epi & (cldn4 >= cldn4_epi_q75)
    is_cldn4_lo = is_epi & (cldn4 <= np.quantile(cldn4[is_epi], 0.25))

    # CD8+ spots: CD8A > 0 (log1p CP10K) and above 75th percentile of nonzero, fallback >0
    cd8_pos = cd8a > 0
    if cd8_pos.sum() >= 20:
        cd8_cut = np.quantile(cd8a[cd8_pos], 0.50)
        is_cd8 = cd8a >= max(cd8_cut, np.quantile(cd8a, 0.75))
    else:
        is_cd8 = cd8_pos

    # same-spot Spearman (all spots and epithelial-only)
    def spr(a, b, mask=None):
        if mask is None:
            mask = np.ones(len(a), dtype=bool)
        aa, bb = a[mask], b[mask]
        if len(aa) < 10:
            return {"rho": np.nan, "p": np.nan, "n": int(len(aa))}
        r, p = spearmanr(aa, bb)
        return {"rho": float(r), "p": float(p), "n": int(len(aa))}

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

    # KRT8 residualization
    cldn4_krt8_resid = residualize(cldn4, krt8)
    resid_all = spr(cldn4_krt8_resid, cd8a)
    resid_epi = spr(cldn4_krt8_resid, cd8a, is_epi)

    # nearest CD8 distance from CLDN4-high / low epithelial spots
    tree_cd8 = cKDTree(xy_um[is_cd8]) if is_cd8.sum() else None

    def nn_dist(mask):
        if tree_cd8 is None or mask.sum() == 0:
            return np.array([])
        d, _ = tree_cd8.query(xy_um[mask], k=1)
        return np.asarray(d, dtype=float)

    d_hi = nn_dist(is_cldn4_hi)
    d_lo = nn_dist(is_cldn4_lo)
    if len(d_hi) and len(d_lo):
        u_stat, u_p = mannwhitneyu(d_hi, d_lo, alternative="two-sided")
    else:
        u_stat, u_p = np.nan, np.nan

    # neighborhood CD8 count in radii (µm)
    radii = [100, 150, 250]
    nhood = {}
    if tree_cd8 is not None:
        for r in radii:
            cnt_hi = tree_cd8.query_ball_point(xy_um[is_cldn4_hi], r=r, return_length=True)
            cnt_lo = tree_cd8.query_ball_point(xy_um[is_cldn4_lo], r=r, return_length=True)
            cnt_hi = np.asarray(cnt_hi, dtype=float)
            cnt_lo = np.asarray(cnt_lo, dtype=float)
            if len(cnt_hi) and len(cnt_lo):
                _, p = mannwhitneyu(cnt_hi, cnt_lo, alternative="two-sided")
            else:
                p = np.nan
            nhood[str(r)] = {
                "radius_um": r,
                "cldn4_hi_mean": float(np.mean(cnt_hi)) if len(cnt_hi) else np.nan,
                "cldn4_hi_median": float(np.median(cnt_hi)) if len(cnt_hi) else np.nan,
                "cldn4_lo_mean": float(np.mean(cnt_lo)) if len(cnt_lo) else np.nan,
                "cldn4_lo_median": float(np.median(cnt_lo)) if len(cnt_lo) else np.nan,
                "mwu_p": float(p) if p == p else np.nan,
                "n_hi": int(len(cnt_hi)),
                "n_lo": int(len(cnt_lo)),
            }

    # also Spearman of CLDN4 vs neighborhood CD8 count among epithelial spots
    nhood_spearman = {}
    if tree_cd8 is not None:
        for r in radii:
            cnt_all = np.asarray(
                tree_cd8.query_ball_point(xy_um, r=r, return_length=True), dtype=float
            )
            nhood_spearman[str(r)] = {
                "all_spots": spr(cldn4, cnt_all),
                "epithelial": spr(cldn4, cnt_all, is_epi),
            }

    summary = {
        "dataset": "10x Visium CytAssist FFPE Human Lung Squamous Cell Carcinoma (demo)",
        "accession": "10x-CytAssist-FFPE-Human-Lung-SCC",
        "platform": "Visium CytAssist FFPE (whole transcriptome probe set v2.0)",
        "n_spots_in_tissue": int(keep.sum()),
        "n_genes": int(len(genes)),
        "genes_present": present,
        "um_per_pixel": float(um_per_px),
        "n_epithelial": int(is_epi.sum()),
        "n_cldn4_high_epi": int(is_cldn4_hi.sum()),
        "n_cldn4_low_epi": int(is_cldn4_lo.sum()),
        "n_cd8_spots": int(is_cd8.sum()),
        "cldn4_q75_in_epi": float(cldn4_epi_q75),
        "same_spot_spearman_all": same_all,
        "same_spot_spearman_epithelial": same_epi,
        "cldn4_krt8_residual_vs_CD8A_all": resid_all,
        "cldn4_krt8_residual_vs_CD8A_epithelial": resid_epi,
        "nearest_cd8_um": {
            "cldn4_hi_median": float(np.median(d_hi)) if len(d_hi) else np.nan,
            "cldn4_hi_mean": float(np.mean(d_hi)) if len(d_hi) else np.nan,
            "cldn4_lo_median": float(np.median(d_lo)) if len(d_lo) else np.nan,
            "cldn4_lo_mean": float(np.mean(d_lo)) if len(d_lo) else np.nan,
            "mwu_U": float(u_stat) if u_stat == u_stat else np.nan,
            "mwu_p": float(u_p) if u_p == u_p else np.nan,
            "n_hi": int(len(d_hi)),
            "n_lo": int(len(d_lo)),
        },
        "nhood_cd8_count": nhood,
        "nhood_spearman_cldn4_vs_cd8count": nhood_spearman,
        "definition": {
            "epithelial": "top 50% z-mean of EPCAM/KRT7/KRT8/KRT18/KRT19/CDH1/ELF3/SFN (CLDN4 excluded)",
            "cldn4_high": "epithelial AND CLDN4 >= Q75 within epithelial spots",
            "cldn4_low": "epithelial AND CLDN4 <= Q25 within epithelial spots",
            "cd8_spot": "CD8A >= max(median of nonzero, Q75 of all spots)",
            "not_used": "TACSTD2 intersection was not applied (CLDN4-only)",
        },
    }
    (TAB / "visium_lusc_summary.json").write_text(json.dumps(summary, indent=2))

    # per-spot table (compact)
    spot_df = pd.DataFrame(
        {
            "barcode": barcodes,
            "x_um": xy_um[:, 0],
            "y_um": xy_um[:, 1],
            "CLDN4": cldn4,
            "CD8A": cd8a,
            "CD3E": cd3e,
            "NKG7": nkg7,
            "KRT8": krt8,
            "CLDN4_KRT8_resid": cldn4_krt8_resid,
            "epi_score": epi_score,
            "is_epithelial": is_epi.astype(int),
            "is_cldn4_high_epi": is_cldn4_hi.astype(int),
            "is_cldn4_low_epi": is_cldn4_lo.astype(int),
            "is_cd8": is_cd8.astype(int),
        }
    )
    spot_df.to_csv(TAB / "visium_lusc_spots.csv.gz", index=False)

    # figures
    plt.rcParams.update({"font.size": 9, "figure.dpi": 140})

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for ax, y, title in [
        (axes[0], cd8a, "CD8A"),
        (axes[1], cd3e, "CD3E"),
        (axes[2], nkg7, "NKG7"),
    ]:
        ax.scatter(cldn4, y, s=4, alpha=0.25, c="#335588", linewidths=0)
        r = same_all[f"CLDN4_vs_{title}"]
        ax.set_xlabel("CLDN4 (log1p CP10K)")
        ax.set_ylabel(f"{title} (log1p CP10K)")
        ax.set_title(f"ρ={r['rho']:.3f}  p={r['p']:.1e}  n={r['n']}")
    fig.suptitle("10x Visium LUSC FFPE — same-spot Spearman (CLDN4-only)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "visium_lusc_same_spot_spearman.png", bbox_inches="tight")
    fig.savefig(FIG / "visium_lusc_same_spot_spearman.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    # spatial CLDN4
    sc0 = axes[0].scatter(
        xy_um[:, 0],
        -xy_um[:, 1],
        c=cldn4,
        s=6,
        cmap="magma",
        linewidths=0,
    )
    axes[0].set_title("CLDN4")
    axes[0].set_aspect("equal")
    axes[0].axis("off")
    fig.colorbar(sc0, ax=axes[0], fraction=0.046, pad=0.02)
    sc1 = axes[1].scatter(
        xy_um[:, 0],
        -xy_um[:, 1],
        c=cd8a,
        s=6,
        cmap="viridis",
        linewidths=0,
    )
    axes[1].set_title("CD8A")
    axes[1].set_aspect("equal")
    axes[1].axis("off")
    fig.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.02)
    fig.suptitle("10x Visium LUSC FFPE — spatial maps (µm)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "visium_lusc_spatial_maps.png", bbox_inches="tight")
    fig.savefig(FIG / "visium_lusc_spatial_maps.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    axes[0].boxplot(
        [d_lo, d_hi],
        tick_labels=["CLDN4-low epi", "CLDN4-high epi"],
        showfliers=False,
    )
    axes[0].set_ylabel("Nearest CD8+ spot (µm)")
    axes[0].set_title(
        f"median {np.median(d_lo):.0f} vs {np.median(d_hi):.0f} µm; MWU p={u_p:.2e}"
        if len(d_hi) and len(d_lo)
        else "nearest CD8"
    )
    r150 = nhood.get("150", {})
    axes[1].boxplot(
        [
            np.asarray(
                tree_cd8.query_ball_point(xy_um[is_cldn4_lo], r=150, return_length=True),
                dtype=float,
            )
            if tree_cd8 is not None
            else np.array([0]),
            np.asarray(
                tree_cd8.query_ball_point(xy_um[is_cldn4_hi], r=150, return_length=True),
                dtype=float,
            )
            if tree_cd8 is not None
            else np.array([0]),
        ],
        tick_labels=["CLDN4-low epi", "CLDN4-high epi"],
        showfliers=False,
    )
    axes[1].set_ylabel("CD8+ spots within 150 µm")
    axes[1].set_title(
        f"mean {r150.get('cldn4_lo_mean', np.nan):.2f} vs {r150.get('cldn4_hi_mean', np.nan):.2f}; "
        f"p={r150.get('mwu_p', np.nan):.2e}"
    )
    fig.suptitle("CLDN4-high vs CLDN4-low epithelial neighborhoods (CLDN4-only)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "visium_lusc_neighborhood.png", bbox_inches="tight")
    fig.savefig(FIG / "visium_lusc_neighborhood.pdf", bbox_inches="tight")
    plt.close(fig)

    return summary


if __name__ == "__main__":
    s = analyze_visium_lusc()
    print(json.dumps(s, indent=2))
