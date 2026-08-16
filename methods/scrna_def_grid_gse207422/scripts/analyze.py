#!/usr/bin/env python3
"""Exhaustive public-UMI definition grid for GSE207422 TACSTD2.

Hu et al. Genome Medicine 2023 (PMID 36869384), GEO GSE207422.
Author CopyKAT / epithelium RDS barcodes are not public. This script
does not reconstruct that object.

Grid axes (public UMI only, 12 post-treatment patients):
  malignant defs: Hu markers (mean and %pos), EPCAM+KRT (pos and mean),
                  CopyKAT-like / inferCNV-like (stromal ref), DRMref public
                  labels, DRMref-like 16-type rebuild. Author = not present.
  immune defs:    T+NK, CD8 only, CD8+NK, CXCL13+, cytotoxicity-high T
  TACSTD2 scores: mean log1p(CP10k), %pos (UMI>0), UCell module

Unit = patient/sample. A3 slide numbers are taken as given and are not
re-argued here. Combinations with Spearman ρ ≤ −0.35 or NMPR>MPR
direction are flagged in the highlight table / figure.
"""
from __future__ import annotations

import argparse
import gzip
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

WIN = 25
UCELL_MAXRANK = 1500
RHO_HIGHLIGHT = -0.35
CYTO_Q = 0.75

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
NORMAL_LUNG = [
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
    "CAPS",
]
EPCAM_KRT = ["EPCAM", "KRT8", "KRT18", "KRT19"]
CYTO_GENES = ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7"]
UCELL_MODULE = ["TACSTD2", "CLDN4", "EPCAM"]
DRMREF_LIKE = {
    "Malignant cells": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CEACAM5", "CEACAM6", "MUC1", "ELF3"],
    "CD8+ T cells": ["CD3D", "CD3E", "CD8A", "CD8B"],
    "CD4+ T cells": ["CD3D", "CD3E", "CD4", "IL7R"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B cells": ["CD79A", "MS4A1", "CD19"],
    "Plasma cells": ["MZB1", "JCHAIN", "SDC1"],
    "Mono/Macro": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophils": ["FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9"],
    "pDCs": ["LILRA4", "IL3RA"],
    "Mast cells": ["TPSAB1", "CPA3", "KIT"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial cells": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
EXTRA = [
    "TACSTD2",
    "CLDN4",
    "PTPRC",
    "CD8A",
    "CD8B",
    "CD4",
    "CXCL13",
    "CEACAM5",
    "CEACAM6",
    "MUC1",
    "ELF3",
    "IL7R",
    "S100A8",
    "S100A9",
    "LILRA4",
    "IL3RA",
    "KIT",
    "NCAM1",
]
DRMREF_TNK = {"CD8+ T cells", "CD4+ T cells", "NK cells"}

MAL_DEFS = [
    "hu_markers_mean",
    "hu_markers_pctpos",
    "epcam_krt_pos",
    "epcam_krt_mean",
    "copykat_stromal_p95",
    "infercnv_stromal_p95",
    "drmref_public",
    "drmref_like",
]
IMMUNE_DEFS = ["T_NK", "CD8_only", "CD8_NK", "CXCL13_pos", "cyto_high_T"]
SCORE_DEFS = ["mean_log1p", "pct_pos", "ucell_module"]

MAL_LABELS = {
    "hu_markers_mean": "Hu markers (mean)",
    "hu_markers_pctpos": "Hu markers (%pos)",
    "epcam_krt_pos": "EPCAM+KRT (pos)",
    "epcam_krt_mean": "EPCAM+KRT (mean)",
    "copykat_stromal_p95": "CopyKAT-like stromal p95",
    "infercnv_stromal_p95": "inferCNV-like stromal p95",
    "drmref_public": "DRMref public labels",
    "drmref_like": "DRMref-like rebuild",
    "author": "author CopyKAT (not public)",
}
IMM_LABELS = {
    "T_NK": "T+NK",
    "CD8_only": "CD8 only",
    "CD8_NK": "CD8+NK",
    "CXCL13_pos": "CXCL13+",
    "cyto_high_T": "cyto-high T",
}
SCORE_LABELS = {
    "mean_log1p": "mean log1p(CP10k)",
    "pct_pos": "%pos UMI>0",
    "ucell_module": "UCell module",
}


def wanted_markers() -> set[str]:
    genes: set[str] = set(EXTRA) | set(NORMAL_LUNG) | set(CYTO_GENES) | set(UCELL_MODULE)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    for vs in DRMREF_LIKE.values():
        genes.update(vs)
    return genes


def sample_from_barcode(bc: str) -> str:
    parts = bc.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else bc


def stream_pass1(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_genes = np.zeros(n, dtype=np.int32)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            n_genes += arr > 0
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  pass1 genes={n_streamed} stored={len(found)}", flush=True)
    print(f"pass1 done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_genes, n_streamed


def stream_ucell(matrix_path: Path, module_umi: dict[str, np.ndarray], max_rank: int = UCELL_MAXRANK):
    """UCell (Andreatta & Carmona 2021) for a small module, one extra stream.

    Rank of each module gene among all genes (competition rank = 1 + n_strictly_above).
    Genes with rank > maxRank are set to maxRank+1. Score is 1 - U/(n_sig * maxRank).
    """
    names = [g for g in module_umi]
    n = next(iter(module_umi.values())).size
    n_above = {g: np.zeros(n, dtype=np.int32) for g in names}
    n_streamed = 0
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n_cells = len(header) - 1
        if n_cells != n:
            raise ValueError(f"UCell header cells {n_cells} != {n}")
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            for g in names:
                n_above[g] += arr > module_umi[g]
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  ucell genes={n_streamed}", flush=True)
    ranks = []
    for g in names:
        rank = n_above[g].astype(np.float32) + 1.0
        rank = np.where(rank > max_rank, float(max_rank + 1), rank)
        ranks.append(rank)
    rank_mat = np.vstack(ranks)
    n_sig = rank_mat.shape[0]
    u_stat = rank_mat.sum(axis=0) - n_sig * (n_sig + 1) / 2.0
    score = 1.0 - u_stat / (n_sig * max_rank)
    score = np.clip(score, 0.0, 1.0).astype(np.float32)
    print(f"ucell done genes={n_streamed} module={names} maxRank={max_rank}", flush=True)
    return score, n_streamed


def stream_pass2(matrix_path: Path, keep_idx: np.ndarray, gene_set: set[str]):
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = keep_idx.size
    gene_names: list[str] = []
    blocks: list[np.ndarray] = []
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n = len(header) - 1
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep or gene not in gene_set:
                n_streamed += 1
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            gene_names.append(gene)
            blocks.append(arr[keep_idx].astype(np.float32, copy=False))
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  pass2 genes_scanned={n_streamed} stored={len(gene_names)}", flush=True)
    mat = np.vstack(blocks) if blocks else np.zeros((0, n_keep), dtype=np.float32)
    print(f"pass2 done stored_genes={len(gene_names)} keep_cells={n_keep}", flush=True)
    return gene_names, mat


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def pctpos_score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [(expr[g] > 0).astype(np.float32) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray, min_top: float = 0.12) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < min_top] = "Unassigned"
    return labels


def moving_average_2d(block: np.ndarray, k: int) -> np.ndarray:
    if block.shape[0] < 3:
        return block.copy()
    k = min(k, block.shape[0] if block.shape[0] % 2 == 1 else block.shape[0] - 1)
    if k < 3:
        return block.copy()
    if k % 2 == 0:
        k -= 1
    pad = k // 2
    padded = np.pad(block, ((pad, pad), (0, 0)), mode="edge")
    c = np.cumsum(np.vstack([np.zeros((1, block.shape[1]), dtype=block.dtype), padded]), axis=0)
    return (c[k:] - c[:-k]) / k


def cnv_scores(
    counts: np.ndarray,
    genes: list[str],
    gene_pos: pd.DataFrame,
    ref_mask: np.ndarray,
    win: int = WIN,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """CopyKAT-like (median ref, sum|smooth|) and inferCNV-like (mean ref, mean|smooth|)."""
    n_g, n_c = counts.shape
    expressed = (counts > 0).sum(axis=1)
    mean_umi = counts.mean(axis=1)
    keep_g = (expressed >= 20) & (mean_umi > 0.05)
    pos = gene_pos.set_index("gene")
    have_pos = np.array([g in pos.index for g in genes])
    keep_g = keep_g & have_pos
    info = {
        "n_genes_input": int(n_g),
        "n_genes_kept": int(keep_g.sum()),
        "n_ref_cells": int(ref_mask.sum()),
        "win": int(win),
    }
    nan = np.full(n_c, np.nan, dtype=np.float32)
    if keep_g.sum() < 200 or ref_mask.sum() < 20:
        info["ok"] = False
        return nan, nan, info
    g_idx = np.where(keep_g)[0]
    gnames = [genes[i] for i in g_idx]
    sub = counts[g_idx]
    lib = np.maximum(sub.sum(axis=0), 1.0)
    logcp = np.log1p(sub / lib * 1e4).astype(np.float32)
    ref_med = np.median(logcp[:, ref_mask], axis=1, keepdims=True)
    ref_mean = np.mean(logcp[:, ref_mask], axis=1, keepdims=True)
    rel_med = logcp - ref_med
    rel_mean = logcp - ref_mean
    chrom = np.array([str(pos.loc[g, "chrom"]) for g in gnames])
    start = np.array([int(pos.loc[g, "start"]) for g in gnames])
    order = np.lexsort((start, chrom))
    rel_med = rel_med[order]
    rel_mean = rel_mean[order]
    chrom = chrom[order]
    sm_med = np.zeros_like(rel_med)
    sm_mean = np.zeros_like(rel_mean)
    chroms_used = []
    for c in pd.unique(chrom):
        idx = np.where(chrom == c)[0]
        if len(idx) < 5:
            continue
        sm_med[idx] = moving_average_2d(rel_med[idx], win)
        sm_mean[idx] = moving_average_2d(rel_mean[idx], win)
        chroms_used.append((c, int(len(idx))))
    copykat = np.sum(np.abs(sm_med), axis=0).astype(np.float32)
    infercnv = np.mean(np.abs(sm_mean), axis=0).astype(np.float32)
    info["ok"] = True
    info["chromosomes"] = chroms_used
    info["copykat_ref_p95"] = float(np.quantile(copykat[ref_mask], 0.95))
    info["infercnv_ref_p95"] = float(np.quantile(infercnv[ref_mask], 0.95))
    return copykat, infercnv, info


def exact_wilcoxon(values: np.ndarray, is_a: np.ndarray) -> dict:
    values = np.asarray(values, float)
    is_a = np.asarray(is_a, bool)
    ok = np.isfinite(values)
    values, is_a = values[ok], is_a[ok]
    n = len(values)
    n_a = int(is_a.sum())
    n_b = n - n_a
    empty = {
        "n_a": n_a,
        "n_b": n_b,
        "U": None,
        "p_exact": None,
        "p_onesided_A_gt_B": None,
        "mean_a": None,
        "mean_b": None,
        "median_a": None,
        "median_b": None,
        "delta_a_minus_b": None,
        "n_perm": 0,
    }
    if n_a < 1 or n_b < 1:
        return empty
    obs_u = float(stats.mannwhitneyu(values[is_a], values[~is_a], alternative="two-sided").statistic)
    expected = n_a * n_b / 2.0
    obs_ext = abs(obs_u - expected)
    obs_mean_a = float(values[is_a].mean())
    obs_mean_b = float(values[~is_a].mean())
    count = 0
    count_ge = 0
    total = 0
    for combo in combinations(range(n), n_a):
        mask = np.zeros(n, dtype=bool)
        mask[list(combo)] = True
        u = float(stats.mannwhitneyu(values[mask], values[~mask], alternative="two-sided").statistic)
        if abs(u - expected) >= obs_ext - 1e-12:
            count += 1
        if values[mask].mean() >= obs_mean_a - 1e-15:
            count_ge += 1
        total += 1
    return {
        "n_a": n_a,
        "n_b": n_b,
        "U": obs_u,
        "p_exact": count / total,
        "p_onesided_A_gt_B": count_ge / total,
        "mean_a": obs_mean_a,
        "mean_b": obs_mean_b,
        "median_a": float(np.median(values[is_a])),
        "median_b": float(np.median(values[~is_a])),
        "delta_a_minus_b": obs_mean_a - obs_mean_b,
        "n_perm": total,
    }


def spearman_safe(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 4:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p)}


def load_meta(metadata: Path, xlsx: Path) -> pd.DataFrame:
    if metadata.exists():
        meta = pd.read_csv(metadata, sep="\t")
    else:
        meta = pd.read_excel(xlsx)
        meta.columns = [c.replace(" ", "_") for c in meta.columns]
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["response"] = meta["Pathologic_Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta["residual_tumor"] = pd.to_numeric(meta["Residual_Tumor"], errors="coerce")
    return meta


def draw_grid(spear: pd.DataFrame, nmpr: pd.DataFrame, out_png: Path) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(9.6, 13.2), gridspec_kw={"height_ratios": [1, 1, 1, 0.72]})
    mal_order = MAL_DEFS
    imm_order = IMMUNE_DEFS
    cmap = plt.cm.RdBu_r
    norm = plt.Normalize(vmin=-0.70, vmax=0.70)

    for ax, score in zip(axes[:3], SCORE_DEFS):
        grid = np.full((len(mal_order), len(imm_order)), np.nan)
        pgrid = np.full_like(grid, np.nan)
        ngrid = np.full_like(grid, np.nan)
        for i, mal in enumerate(mal_order):
            for j, imm in enumerate(imm_order):
                hit = spear[(spear["malignant_def"] == mal) & (spear["immune_def"] == imm) & (spear["score"] == score)]
                if len(hit):
                    grid[i, j] = hit.iloc[0]["rho"] if pd.notna(hit.iloc[0]["rho"]) else np.nan
                    pgrid[i, j] = hit.iloc[0]["p"] if pd.notna(hit.iloc[0]["p"]) else np.nan
                    ngrid[i, j] = hit.iloc[0]["n"] if pd.notna(hit.iloc[0]["n"]) else np.nan
        im = ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(len(imm_order)))
        ax.set_xticklabels([IMM_LABELS[k] for k in imm_order], fontsize=8)
        ax.set_yticks(range(len(mal_order)))
        ax.set_yticklabels([MAL_LABELS[k] for k in mal_order], fontsize=7.5)
        ax.set_title(f"TACSTD2 {SCORE_LABELS[score]} vs immune fraction  (post-tx patients)", fontsize=10)
        for i in range(len(mal_order)):
            for j in range(len(imm_order)):
                rho = grid[i, j]
                if not np.isfinite(rho):
                    ax.text(j, i, "NA", ha="center", va="center", fontsize=6.5, color="#555")
                    continue
                p = pgrid[i, j]
                n = ngrid[i, j]
                txt = f"{rho:+.2f}\np={p:.2f}\nn={int(n)}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.2, color="black")
                if rho <= RHO_HIGHLIGHT:
                    ax.add_patch(
                        mpatches.Rectangle(
                            (j - 0.48, i - 0.48),
                            0.96,
                            0.96,
                            fill=False,
                            edgecolor="#c9a227",
                            lw=2.2,
                        )
                    )
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ")

    ax = axes[3]
    grid = np.full((len(mal_order), len(SCORE_DEFS)), np.nan)
    pgrid = np.full_like(grid, np.nan)
    dgrid = np.full_like(grid, np.nan)
    for i, mal in enumerate(mal_order):
        for j, score in enumerate(SCORE_DEFS):
            hit = nmpr[(nmpr["malignant_def"] == mal) & (nmpr["score"] == score)]
            if len(hit):
                dgrid[i, j] = hit.iloc[0]["delta_NMPR_minus_MPR"] if pd.notna(hit.iloc[0]["delta_NMPR_minus_MPR"]) else np.nan
                pgrid[i, j] = hit.iloc[0]["p"] if pd.notna(hit.iloc[0]["p"]) else np.nan
                grid[i, j] = 1.0 if (np.isfinite(dgrid[i, j]) and dgrid[i, j] > 0) else (-1.0 if np.isfinite(dgrid[i, j]) else np.nan)
    im = ax.imshow(dgrid, cmap=cmap, aspect="auto")
    # autoscale is fine for mixed units; annotate instead of forcing a shared scale
    ax.set_xticks(range(len(SCORE_DEFS)))
    ax.set_xticklabels([SCORE_LABELS[k] for k in SCORE_DEFS], fontsize=8)
    ax.set_yticks(range(len(mal_order)))
    ax.set_yticklabels([MAL_LABELS[k] for k in mal_order], fontsize=7.5)
    ax.set_title("NMPR − MPR (exact Wilcoxon; gold box = NMPR>MPR)", fontsize=10)
    for i in range(len(mal_order)):
        for j in range(len(SCORE_DEFS)):
            d = dgrid[i, j]
            if not np.isfinite(d):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=6.5, color="#555")
                continue
            p = pgrid[i, j]
            ax.text(j, i, f"Δ={d:+.3f}\np={p:.2f}", ha="center", va="center", fontsize=6.5)
            if d > 0:
                ax.add_patch(
                    mpatches.Rectangle(
                        (j - 0.48, i - 0.48),
                        0.96,
                        0.96,
                        fill=False,
                        edgecolor="#c9a227",
                        lw=2.2,
                    )
                )
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Δ (NMPR−MPR)")
    fig.suptitle(
        "GSE207422 public UMI definition grid  ·  gold box: ρ≤−0.35 or NMPR>MPR  ·  author CopyKAT IDs not public",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    ap.add_argument("--xlsx", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--metadata", type=Path, default=Path("data/GSE207422/geo_scRNAseq_sample_metadata.tsv"))
    ap.add_argument("--gene-chr", type=Path, default=here / "assets" / "gene_chr.tsv")
    ap.add_argument("--drmref", type=Path, default=here / "assets" / "drmref_cell_annotation.tsv.gz")
    ap.add_argument("--outdir", type=Path, default=here)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    meta = load_meta(args.metadata, args.xlsx)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    gene_chr = pd.read_csv(args.gene_chr, sep="\t")
    gene_chr["chrom"] = gene_chr["chrom"].astype(str)
    cnv_gene_set = set(gene_chr["gene"])
    markers = wanted_markers()

    print("pass1: markers + library sizes", flush=True)
    cell_ids, expr, n_umi, n_genes, n_matrix_genes = stream_pass1(args.matrix, markers)
    n = len(cell_ids)
    cells = np.array(cell_ids)
    sample = np.array([sample_from_barcode(c) for c in cells])
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}

    print("ucell module", UCELL_MODULE, flush=True)
    module_umi = {g: expr[g] for g in UCELL_MODULE if g in expr}
    if len(module_umi) < 1:
        raise RuntimeError("TACSTD2 module genes missing from matrix")
    ucell, _ = stream_ucell(args.matrix, module_umi)

    scores_mean = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    scores_pos = {name: pctpos_score(expr, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    normal_score = module_score(log_cp, NORMAL_LUNG, n)
    normal_umi = np.zeros(n, dtype=np.float32)
    for g in NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]
    cd3_log = log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n)))
    lineage_mean = assign_lineage(scores_mean, cd3_log)
    lineage_pos = assign_lineage(scores_pos, cd3_log, min_top=0.15)

    drm_like_scores = {name: module_score(log_cp, genes, n) for name, genes in DRMREF_LIKE.items()}
    drm_like_mat = np.vstack([drm_like_scores[k] for k in DRMREF_LIKE])
    drm_like_best = np.argmax(drm_like_mat, axis=0)
    drm_like_top = drm_like_mat[drm_like_best, np.arange(n)]
    drm_like_lab = np.array(list(DRMREF_LIKE), dtype=object)[drm_like_best]
    drm_like_lab[drm_like_top < 0.12] = "Unassigned"
    # DRMref has no normal-epithelial class; subtract alveolar/club/ciliated so the
    # rebuild is not just "all epithelium".
    is_drm_like_mal = (drm_like_lab == "Malignant cells") & (normal_score < 0.25)

    tac = expr.get("TACSTD2", np.zeros(n, dtype=np.float32))
    cd3e = expr.get("CD3E", np.zeros(n, dtype=np.float32))
    cd8a = expr.get("CD8A", np.zeros(n, dtype=np.float32))
    cd8b = expr.get("CD8B", np.zeros(n, dtype=np.float32))
    nkg7 = expr.get("NKG7", np.zeros(n, dtype=np.float32))
    cxcl13 = expr.get("CXCL13", np.zeros(n, dtype=np.float32))
    epcam = expr.get("EPCAM", np.zeros(n, dtype=np.float32))
    krt8 = expr.get("KRT8", np.zeros(n, dtype=np.float32))
    krt18 = expr.get("KRT18", np.zeros(n, dtype=np.float32))
    krt19 = expr.get("KRT19", np.zeros(n, dtype=np.float32))
    ptprc = expr.get("PTPRC", np.zeros(n, dtype=np.float32))

    is_epi_mean = lineage_mean == "Epithelial"
    is_epi_pos = lineage_pos == "Epithelial"
    is_stromal = np.isin(lineage_mean, ["Fibroblast", "Endothelial"])
    is_t = lineage_mean == "T"
    is_nk = lineage_mean == "NK"
    epcam_krt_mod = module_score(log_cp, EPCAM_KRT, n)
    # mean-style EPCAM+KRT: high module among epithelial-range cells.
    # The 80th percentile of all cells is near zero (zero-inflated), so the
    # cutoff is the median of Hu-mean epithelial scores, with a floor of 0.40.
    if is_epi_mean.sum() >= 20:
        ek_cut = float(max(0.40, np.median(epcam_krt_mod[is_epi_mean])))
    else:
        ek_cut = 0.40
    is_ek_mean = (epcam_krt_mod >= ek_cut) & (ptprc < 1)
    is_ek_pos = (epcam >= 1) & ((krt8 >= 1) | (krt18 >= 1) | (krt19 >= 1))

    # Hu malignant: epithelial AND not a strong normal-lung program
    if is_epi_mean.sum() >= 50:
        nl_cut_mean = float(np.quantile(normal_score[is_epi_mean], 0.60))
    else:
        nl_cut_mean = 0.25
    is_hu_mean = is_epi_mean & (normal_score < nl_cut_mean)
    is_hu_pos = is_epi_pos & (normal_umi == 0)

    drm_type = np.array([""] * n, dtype=object)
    if args.drmref.exists():
        drm = pd.read_csv(args.drmref, sep="\t")
        dmap = dict(zip(drm["cell_barcode"].astype(str), drm["celltype"].astype(str)))
        drm_type = np.array([dmap.get(c, "") for c in cells], dtype=object)
    is_drm_mal = drm_type == "Malignant cells"

    cyto = module_score(log_cp, CYTO_GENES, n)
    if is_t.sum() >= 20:
        cyto_cut = float(np.quantile(cyto[is_t], CYTO_Q))
    else:
        cyto_cut = float(np.quantile(cyto, CYTO_Q))
    is_cd8 = is_t & ((cd8a >= 1) | (cd8b >= 1))
    is_tnk = is_t | is_nk
    is_cd8_nk = is_cd8 | is_nk
    is_cxcl13 = is_tnk & (cxcl13 >= 1)
    is_cyto_high_t = is_t & (cyto >= cyto_cut)

    print(
        f"lineage epi_mean={is_epi_mean.sum()} epi_pos={is_epi_pos.sum()} "
        f"stromal={is_stromal.sum()} T={is_t.sum()} NK={is_nk.sum()} "
        f"hu_mean={is_hu_mean.sum()} hu_pos={is_hu_pos.sum()} "
        f"ek_pos={is_ek_pos.sum()} ek_mean={is_ek_mean.sum()} "
        f"drmref={is_drm_mal.sum()} drm_like={is_drm_like_mal.sum()} "
        f"CD8={is_cd8.sum()} CXCL13+TNK={is_cxcl13.sum()} cytoT={is_cyto_high_t.sum()}",
        flush=True,
    )

    keep = is_epi_mean | is_stromal | is_ek_pos
    keep_idx = np.where(keep)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    print(f"pass2 keep_cells={keep_idx.size}", flush=True)
    cnv_genes, cnv_mat = stream_pass2(args.matrix, keep_idx, cnv_gene_set)
    ref_keep = is_stromal[keep_idx]
    copykat_s, infer_s, cnv_info = cnv_scores(cnv_mat, cnv_genes, gene_chr, ref_keep, WIN)
    print("CNV info", json.dumps({k: v for k, v in cnv_info.items() if k != "chromosomes"}), flush=True)

    copykat_full = np.full(n, np.nan, dtype=np.float32)
    infer_full = np.full(n, np.nan, dtype=np.float32)
    for i, j in keep_map.items():
        copykat_full[i] = copykat_s[j]
        infer_full[i] = infer_s[j]
    ck_p95 = float(cnv_info.get("copykat_ref_p95", np.nan))
    ic_p95 = float(cnv_info.get("infercnv_ref_p95", np.nan))
    is_copykat = is_epi_mean & np.isfinite(copykat_full) & (copykat_full > ck_p95)
    is_infercnv = is_epi_mean & np.isfinite(infer_full) & (infer_full > ic_p95)
    print(f"malignant copykat={is_copykat.sum()} infercnv={is_infercnv.sum()}", flush=True)

    mal_masks = {
        "hu_markers_mean": is_hu_mean,
        "hu_markers_pctpos": is_hu_pos,
        "epcam_krt_pos": is_ek_pos,
        "epcam_krt_mean": is_ek_mean,
        "copykat_stromal_p95": is_copykat,
        "infercnv_stromal_p95": is_infercnv,
        "drmref_public": is_drm_mal,
        "drmref_like": is_drm_like_mal,
    }
    imm_masks = {
        "T_NK": is_tnk,
        "CD8_only": is_cd8,
        "CD8_NK": is_cd8_nk,
        "CXCL13_pos": is_cxcl13,
        "cyto_high_T": is_cyto_high_t,
    }

    tac_log = np.log1p(tac * scale).astype(np.float32)
    cell_df = pd.DataFrame(
        {
            "barcode": cells,
            "sample": sample,
            "n_umi": n_umi,
            "TACSTD2": tac,
            "tacstd2_log1p_cp10k": tac_log,
            "ucell_module": ucell,
            "lineage_mean": lineage_mean,
            "drmref_celltype": drm_type,
            "drmref_like": drm_like_lab,
        }
    )
    for name, mask in mal_masks.items():
        cell_df[f"mal_{name}"] = mask.astype(np.int8)
    for name, mask in imm_masks.items():
        cell_df[f"imm_{name}"] = mask.astype(np.int8)
    cell_df.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    smap = meta.set_index("Sample")
    rows = []
    for samp, g in cell_df.groupby("sample", sort=True):
        if samp not in smap.index:
            continue
        rec = smap.loc[samp]
        n_s = len(g)
        row = {
            "sample": samp,
            "patient": rec["Patient"],
            "is_post": bool(rec["is_post"]),
            "response_raw": rec["Pathologic_Response"],
            "response": rec["response"],
            "residual_tumor": rec["residual_tumor"],
            "pathology": rec["Pathology"],
            "n_cells": n_s,
        }
        for name in IMMUNE_DEFS:
            row[f"n_{name}"] = int(g[f"imm_{name}"].sum())
            row[f"frac_{name}"] = float(g[f"imm_{name}"].mean())
        for name in MAL_DEFS:
            sub = g.loc[g[f"mal_{name}"].astype(bool)]
            row[f"{name}_n"] = int(len(sub))
            if len(sub) == 0:
                row[f"{name}_mean_log1p"] = np.nan
                row[f"{name}_pct_pos"] = np.nan
                row[f"{name}_ucell_module"] = np.nan
            else:
                row[f"{name}_mean_log1p"] = float(sub["tacstd2_log1p_cp10k"].mean())
                row[f"{name}_pct_pos"] = float((sub["TACSTD2"] > 0).mean())
                row[f"{name}_ucell_module"] = float(sub["ucell_module"].mean())
        rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(args.outdir / "per_sample.tsv", sep="\t", index=False)

    post = per[per["is_post"]].copy()
    spear_rows = []
    nmpr_rows = []
    for mal in MAL_DEFS:
        for score in SCORE_DEFS:
            col = f"{mal}_{score}"
            use = post[np.isfinite(post[col])].copy()
            is_nmpr = use["response"].to_numpy() == "NMPR"
            w = exact_wilcoxon(use[col].to_numpy(), is_nmpr)
            nmpr_rows.append(
                {
                    "family": "NMPR_vs_MPR",
                    "malignant_def": mal,
                    "malignant_label": MAL_LABELS[mal],
                    "score": score,
                    "score_label": SCORE_LABELS[score],
                    "n": int(len(use)),
                    "n_NMPR": int(is_nmpr.sum()) if len(use) else 0,
                    "n_MPR": int((~is_nmpr).sum()) if len(use) else 0,
                    "mean_NMPR": w["mean_a"],
                    "mean_MPR": w["mean_b"],
                    "median_NMPR": w["median_a"],
                    "median_MPR": w["median_b"],
                    "delta_NMPR_minus_MPR": w["delta_a_minus_b"],
                    "U": w["U"],
                    "p": w["p_exact"],
                    "p_onesided_NMPR_gt_MPR": w["p_onesided_A_gt_B"],
                    "nmpr_gt_mpr": bool(w["delta_a_minus_b"] is not None and w["delta_a_minus_b"] > 0),
                    "highlight": bool(w["delta_a_minus_b"] is not None and w["delta_a_minus_b"] > 0),
                    "note": (
                        f"exact C({w['n_a']+w['n_b']},{w['n_a']})={w['n_perm']}; "
                        "pCR P06 counted as MPR; author CopyKAT IDs not public"
                    ),
                }
            )
            for imm in IMMUNE_DEFS:
                tcol = f"frac_{imm}"
                u2 = use[np.isfinite(use[tcol])].copy()
                s = spearman_safe(u2[col], u2[tcol])
                rho = s["rho"]
                spear_rows.append(
                    {
                        "family": "vs_immune",
                        "malignant_def": mal,
                        "malignant_label": MAL_LABELS[mal],
                        "immune_def": imm,
                        "immune_label": IMM_LABELS[imm],
                        "score": score,
                        "score_label": SCORE_LABELS[score],
                        "n": s["n"],
                        "n_NMPR": int((u2["response"] == "NMPR").sum()) if len(u2) else 0,
                        "n_MPR": int((u2["response"] == "MPR").sum()) if len(u2) else 0,
                        "rho": rho,
                        "p": s["p"],
                        "highlight_rho_le_neg035": bool(rho is not None and rho <= RHO_HIGHLIGHT),
                        "note": "Spearman; unit=post-treatment patient; scipy.stats.spearmanr",
                    }
                )

    spear = pd.DataFrame(spear_rows)
    nmpr = pd.DataFrame(nmpr_rows)
    spear.to_csv(args.outdir / "grid_spearman.tsv", sep="\t", index=False)
    nmpr.to_csv(args.outdir / "grid_nmpr_mpr.tsv", sep="\t", index=False)

    hi_rho = spear[spear["highlight_rho_le_neg035"]].copy()
    hi_nmpr = nmpr[nmpr["highlight"]].copy()
    hi_rho.to_csv(args.outdir / "highlight_rho_le_neg035.tsv", sep="\t", index=False)
    hi_nmpr.to_csv(args.outdir / "highlight_nmpr_gt_mpr.tsv", sep="\t", index=False)

    counts = {
        "n_cells_public_matrix": int(n),
        "n_genes_matrix": int(n_matrix_genes),
        "user_claimed_cells": 90652,
        "author_copykat_ids_public": False,
        "n_epithelial_hu_mean": int(is_epi_mean.sum()),
        "n_epithelial_hu_pctpos": int(is_epi_pos.sum()),
        "n_stromal": int(is_stromal.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_CD8": int(is_cd8.sum()),
        "n_CXCL13_TNK": int(is_cxcl13.sum()),
        "n_cyto_high_T": int(is_cyto_high_t.sum()),
        "cyto_cut_log1p_cp10k": cyto_cut,
        "epcam_krt_mean_cut": ek_cut,
        "normal_lung_cut_hu_mean": nl_cut_mean,
        "ucell_module_genes": UCELL_MODULE,
        "ucell_maxRank": UCELL_MAXRANK,
        "malignant_n": {k: int(v.sum()) for k, v in mal_masks.items()},
        "cnv": {k: v for k, v in cnv_info.items() if k != "chromosomes"},
        "n_grid_spearman": int(len(spear)),
        "n_grid_nmpr": int(len(nmpr)),
        "n_highlight_rho_le_neg035": int(len(hi_rho)),
        "n_highlight_nmpr_gt_mpr": int(len(hi_nmpr)),
    }
    with open(args.outdir / "summary.json", "w") as fh:
        json.dump(counts, fh, indent=2, default=str)
    with open(args.outdir / "cnv_info.json", "w") as fh:
        json.dump(cnv_info, fh, indent=2, default=str)

    lin_counts = (
        cell_df.groupby(["sample", "lineage_mean"]).size().rename("n").reset_index().pivot_table(
            index="sample", columns="lineage_mean", values="n", fill_value=0
        )
    )
    lin_counts.to_csv(args.outdir / "lineage_counts_by_sample.tsv", sep="\t")

    draw_grid(spear, nmpr, args.outdir / "fig_def_grid.png")
    print("wrote", args.outdir, "spearman", len(spear), "nmpr", len(nmpr), flush=True)
    print("highlight ρ≤-0.35:", len(hi_rho), "NMPR>MPR:", len(hi_nmpr), flush=True)


if __name__ == "__main__":
    main()
