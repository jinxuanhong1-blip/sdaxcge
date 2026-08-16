#!/usr/bin/env python3
"""Shared lineage, CopyKAT-like CNV, and sample-level tests for the malignant-definition grid."""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

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
EXTRA = [
    "TACSTD2",
    "CLDN4",
    "PTPRC",
    "CD8A",
    "CD4",
    "NCAM1",
    "MKI67",
    "DST",
    "SERPINB9",
    "TOP2A",
    "PCNA",
    "KRT17",
    "KRT5",
    "TP63",
]
DRMREF_TNK = {"CD8+ T cells", "CD4+ T cells", "NK cells"}
WIN = 25
MIN_MAL = 5


def wanted_markers() -> set[str]:
    genes: set[str] = set(EXTRA) | set(NORMAL_LUNG)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    return genes


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
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
    labels[top < 0.12] = "Unassigned"
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


def copykat_like(
    counts: np.ndarray,
    genes: list[str],
    gene_pos: pd.DataFrame,
    ref_mask: np.ndarray,
    win: int = WIN,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Window-smoothed expression CNV vs a diploid reference (CopyKAT/inferCNV-style).

    counts: n_genes x n_cells raw UMI.
    Returns (score_sumabs, score_meanabs, info).
    """
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
    if keep_g.sum() < 200 or ref_mask.sum() < 20:
        info["ok"] = False
        z = np.full(n_c, np.nan, dtype=np.float32)
        return z, z, info
    g_idx = np.where(keep_g)[0]
    gnames = [genes[i] for i in g_idx]
    sub = counts[g_idx]
    lib = np.maximum(sub.sum(axis=0), 1.0)
    logcp = np.log1p(sub / lib * 1e4).astype(np.float32)
    ref_med = np.median(logcp[:, ref_mask], axis=1, keepdims=True)
    rel = logcp - ref_med
    chrom = np.array([str(pos.loc[g, "chrom"]) for g in gnames])
    start = np.array([int(pos.loc[g, "start"]) for g in gnames])
    order = np.lexsort((start, chrom))
    rel = rel[order]
    chrom = chrom[order]
    smoothed = np.zeros_like(rel)
    chroms_used = []
    for c in pd.unique(chrom):
        idx = np.where(chrom == c)[0]
        if len(idx) < 5:
            continue
        smoothed[idx] = moving_average_2d(rel[idx], win)
        chroms_used.append((c, int(len(idx))))
    score_sumabs = np.sum(np.abs(smoothed), axis=0).astype(np.float32)
    score_meanabs = np.mean(np.abs(smoothed), axis=0).astype(np.float32)
    info["ok"] = True
    info["chromosomes"] = chroms_used
    info["ref_score_median"] = float(np.median(score_sumabs[ref_mask]))
    info["ref_score_p95"] = float(np.quantile(score_sumabs[ref_mask], 0.95))
    return score_sumabs, score_meanabs, info


def kmeans1d_cut(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return float(np.nan)
    c1, c2 = np.quantile(x, [0.25, 0.75])
    for _ in range(40):
        d1 = np.abs(x - c1)
        d2 = np.abs(x - c2)
        m1 = x[d1 <= d2]
        m2 = x[d2 < d1]
        if m1.size == 0 or m2.size == 0:
            break
        c1, c2 = float(m1.mean()), float(m2.mean())
    lo, hi = (c1, c2) if c1 < c2 else (c2, c1)
    return float(0.5 * (lo + hi))


def exact_wilcoxon(values: np.ndarray, is_a: np.ndarray) -> dict:
    """Exact two-sided Wilcoxon by enumerating label assignments."""
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
        "p_onesided_A_gt_B_mean": None,
        "mean_a": None,
        "mean_b": None,
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
    # C(n, k) explodes; fall back to asymptotic MWU when the grid is large.
    from math import comb

    n_perm = comb(n, n_a)
    if n_perm > 30000:
        p_as = float(stats.mannwhitneyu(values[is_a], values[~is_a], alternative="two-sided").pvalue)
        p_one = float(stats.mannwhitneyu(values[is_a], values[~is_a], alternative="greater").pvalue)
        return {
            "n_a": n_a,
            "n_b": n_b,
            "U": obs_u,
            "p_exact": p_as,
            "p_onesided_A_gt_B_mean": p_one,
            "mean_a": obs_mean_a,
            "mean_b": obs_mean_b,
            "delta_a_minus_b": obs_mean_a - obs_mean_b,
            "n_perm": n_perm,
            "p_method": "asymptotic_mwu",
        }
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
        "p_onesided_A_gt_B_mean": count_ge / total,
        "mean_a": obs_mean_a,
        "mean_b": obs_mean_b,
        "delta_a_minus_b": obs_mean_a - obs_mean_b,
        "n_perm": total,
        "p_method": "exact_enumeration",
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


def direction_nmpr(delta: float | None) -> str:
    if delta is None or not np.isfinite(delta):
        return "NA"
    return "NMPR>MPR" if delta > 0 else ("NMPR<MPR" if delta < 0 else "tie")


def direction_rho(rho: float | None) -> str:
    if rho is None or not np.isfinite(rho):
        return "NA"
    return "negative" if rho < 0 else ("positive" if rho > 0 else "zero")


def recovers_user_direction(delta, rho) -> str:
    d = direction_nmpr(delta)
    r = direction_rho(rho)
    if d == "NA" or r == "NA":
        return "NA"
    if d == "NMPR>MPR" and r == "negative":
        return "both"
    if d == "NMPR>MPR":
        return "NMPR>MPR_only"
    if r == "negative":
        return "rho_neg_only"
    return "neither"


def score_compartment(sub: pd.DataFrame) -> dict:
    if len(sub) == 0:
        return {
            "n": 0,
            "mean_log1p_cp10k": np.nan,
            "mean_log1p_umi": np.nan,
            "pct_pos_ge1": np.nan,
            "pct_pos_ge2": np.nan,
            "pct_pos_ge3": np.nan,
            "pseudobulk_cpm": np.nan,
            "mean_cnv": np.nan,
        }
    lib = float(sub["n_umi"].sum()) if "n_umi" in sub.columns else np.nan
    tac = sub["TACSTD2"].to_numpy(float)
    out = {
        "n": int(len(sub)),
        "mean_log1p_cp10k": float(sub["tacstd2_log1p_cp10k"].mean()),
        "mean_log1p_umi": float(np.log1p(tac).mean()),
        "pct_pos_ge1": float((tac >= 1).mean()),
        "pct_pos_ge2": float((tac >= 2).mean()),
        "pct_pos_ge3": float((tac >= 3).mean()),
        "pseudobulk_cpm": float(tac.sum() / lib * 1e6) if lib and np.isfinite(lib) and lib > 0 else np.nan,
        "mean_cnv": float(sub["cnv_sumabs"].mean()) if "cnv_sumabs" in sub.columns and sub["cnv_sumabs"].notna().any() else np.nan,
    }
    return out
