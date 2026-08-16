#!/usr/bin/env python3
"""A3 wave-2: recompute malignant calls on public GSE207422 UMI.

Hu et al. Genome Medicine 2023 (PMID 36869384). Author CopyKAT barcode IDs
are not public. This script:

  1. Streams the GEO UMI matrix (92,330 cells).
  2. Assigns lineages from published canonical markers.
  3. Scores CopyKAT-like / inferCNV-like aneuploidy on epithelial cells
     using stromal (fibroblast + endothelial) as the diploid reference.
  4. Tests malignant-only TACSTD2 vs MPR and vs T/NK under several
     definitions, including %positive at UMI>=1/2/3 and residuals after
     regressing out epithelial fraction.

This is NOT the authors' CopyKAT object. It is a window-smoothed
expression-CNV score on the public UMI, documented as such.
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
    """Gene-axis moving average. block is n_genes x n_cells."""
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
    """Window-smoothed expression CNV vs a diploid reference.

    counts: n_genes x n_cells raw UMI (keep-cell subset).
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
    if n_a < 1 or n_b < 1:
        return {"n_a": n_a, "n_b": n_b, "U": None, "p_exact": None, "mean_a": None, "mean_b": None}
    obs_u = float(stats.mannwhitneyu(values[is_a], values[~is_a], alternative="two-sided").statistic)
    expected = n_a * n_b / 2.0

    def ext(u: float) -> float:
        return abs(u - expected)

    obs_ext = ext(obs_u)
    count = 0
    total = 0
    one_gt = 0
    obs_mean_a = float(values[is_a].mean())
    obs_mean_b = float(values[~is_a].mean())
    for combo in combinations(range(n), n_a):
        mask = np.zeros(n, dtype=bool)
        mask[list(combo)] = True
        u = float(stats.mannwhitneyu(values[mask], values[~mask], alternative="two-sided").statistic)
        if ext(u) >= obs_ext - 1e-12:
            count += 1
        if values[mask].mean() >= obs_mean_a - 1e-15:
            # one-sided: assignments with group-A mean as high or higher
            one_gt += 1
        total += 1
    # one-sided p for A > B uses rank U in the greater direction
    count_ge = 0
    for combo in combinations(range(n), n_a):
        mask = np.zeros(n, dtype=bool)
        mask[list(combo)] = True
        u = float(stats.mannwhitneyu(values[mask], values[~mask], alternative="greater").statistic)
        # use mean difference for one-sided to avoid U-tie confusion
        if values[mask].mean() >= obs_mean_a - 1e-15:
            count_ge += 1
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


def residualize(y, x) -> np.ndarray:
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = np.isfinite(y) & np.isfinite(x)
    out = np.full(y.shape, np.nan, dtype=float)
    if m.sum() < 3:
        return out
    xv = x[m]
    yv = y[m]
    A = np.vstack([np.ones(m.sum()), xv]).T
    coef, _, _, _ = np.linalg.lstsq(A, yv, rcond=None)
    out[m] = yv - (coef[0] + coef[1] * xv)
    return out


def rank_residualize(y, x) -> np.ndarray:
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = np.isfinite(y) & np.isfinite(x)
    out = np.full(y.shape, np.nan, dtype=float)
    if m.sum() < 3:
        return out
    yr = stats.rankdata(y[m])
    xr = stats.rankdata(x[m])
    A = np.vstack([np.ones(m.sum()), xr]).T
    coef, _, _, _ = np.linalg.lstsq(A, yr, rcond=None)
    out[m] = yr - (coef[0] + coef[1] * xr)
    return out


def match_nmpr(delta: float | None, p: float | None) -> str:
    if delta is None or p is None:
        return "NA"
    if delta > 0 and p < 0.05:
        return "MATCH"
    if delta > 0:
        return "PARTIAL_direction_NS"
    return "MISMATCH"


def match_rho(rho: float | None) -> str:
    if rho is None:
        return "NA"
    if -0.50 <= rho <= -0.40:
        return "MATCH"
    if -0.55 <= rho <= -0.35:
        return "NEAR"
    return "MISMATCH"


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    ap.add_argument("--metadata", type=Path, default=Path("data/GSE207422/geo_scRNAseq_sample_metadata.tsv"))
    ap.add_argument("--xlsx", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--gene-chr", type=Path, default=Path("scripts/rework/A3_wave2/gene_chr.tsv"))
    ap.add_argument("--drmref", type=Path, default=Path("scripts/rework/A3_wave2/drmref_cell_annotation.tsv.gz"))
    ap.add_argument("--outdir", type=Path, default=Path("results/rework/A3_wave2"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.metadata.exists():
        meta = pd.read_csv(args.metadata, sep="\t")
    else:
        meta = pd.read_excel(args.xlsx)
        meta.columns = [c.replace(" ", "_") for c in meta.columns]
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["response"] = meta["Pathologic_Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta["residual_tumor"] = pd.to_numeric(meta["Residual_Tumor"], errors="coerce")
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

    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    normal_score = module_score(log_cp, NORMAL_LUNG, n)
    lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n))))
    tac = expr.get("TACSTD2", np.zeros(n, dtype=np.float32))
    cd3e = expr.get("CD3E", np.zeros(n, dtype=np.float32))
    cd8a = expr.get("CD8A", np.zeros(n, dtype=np.float32))
    nkg7 = expr.get("NKG7", np.zeros(n, dtype=np.float32))

    is_epi = lineage == "Epithelial"
    is_stromal = np.isin(lineage, ["Fibroblast", "Endothelial"])
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    is_myeloid = lineage == "Myeloid"
    is_lineage_tnk = is_t | is_nk
    is_umi_tnk = (cd3e >= 1) | (cd8a >= 1) | (nkg7 >= 1)
    is_umi_t = cd3e >= 1
    is_umi_nk = (nkg7 >= 1) & (cd3e < 1)
    is_umi_strict_tnk = is_umi_t | is_umi_nk

    drm = None
    drm_type = np.array([""] * n, dtype=object)
    if args.drmref.exists():
        drm = pd.read_csv(args.drmref, sep="\t")
        dmap = dict(zip(drm["cell_barcode"].astype(str), drm["celltype"].astype(str)))
        drm_type = np.array([dmap.get(c, "") for c in cells], dtype=object)
    is_drm_mal = drm_type == "Malignant cells"
    is_drm_tnk = np.isin(drm_type, list(DRMREF_TNK))
    is_drm_annot = drm_type != ""

    print(
        f"lineage epi={is_epi.sum()} stromal={is_stromal.sum()} "
        f"T={is_t.sum()} NK={is_nk.sum()} myeloid={is_myeloid.sum()}",
        flush=True,
    )

    # Keep epithelial + all stromal + a cap of immune for alternate reference.
    rng = np.random.default_rng(0)
    immune_idx = np.where(is_t | is_myeloid)[0]
    if immune_idx.size > 2500:
        immune_idx = rng.choice(immune_idx, size=2500, replace=False)
    keep = np.zeros(n, dtype=bool)
    keep[is_epi | is_stromal] = True
    keep[immune_idx] = True
    keep_idx = np.where(keep)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    print(f"pass2 keep_cells={keep_idx.size} (epi+stromal+immune-ref)", flush=True)

    cnv_genes, cnv_mat = stream_pass2(args.matrix, keep_idx, cnv_gene_set)
    ref_keep = is_stromal[keep_idx]
    if ref_keep.sum() < 80:
        print(f"stromal ref n={ref_keep.sum()} < 80; adding T+myeloid subsample", flush=True)
        ref_keep = ref_keep | np.isin(keep_idx, immune_idx)
    score_sumabs, score_meanabs, cnv_info = copykat_like(cnv_mat, cnv_genes, gene_chr, ref_keep, WIN)
    print("CNV info", json.dumps({k: v for k, v in cnv_info.items() if k != "chromosomes"}), flush=True)

    cnv_full = np.full(n, np.nan, dtype=np.float32)
    cnv_mean_full = np.full(n, np.nan, dtype=np.float32)
    for i, j in keep_map.items():
        cnv_full[i] = score_sumabs[j]
        cnv_mean_full[i] = score_meanabs[j]

    ref_scores = score_sumabs[ref_keep]
    ref_scores = ref_scores[np.isfinite(ref_scores)]
    p95 = float(np.quantile(ref_scores, 0.95)) if ref_scores.size else np.nan
    p90 = float(np.quantile(ref_scores, 0.90)) if ref_scores.size else np.nan
    mu2 = float(ref_scores.mean() + 2 * ref_scores.std()) if ref_scores.size else np.nan
    gmm_cut = kmeans1d_cut(cnv_full[is_epi])
    cnv_info.update({"threshold_p95": p95, "threshold_p90": p90, "threshold_mean2sd": mu2, "threshold_gmm": gmm_cut})

    # Primary malignant: epithelial AND CNV > stromal/reference 95th percentile.
    is_mal_p95 = is_epi & np.isfinite(cnv_full) & (cnv_full > p95)
    is_mal_p90 = is_epi & np.isfinite(cnv_full) & (cnv_full > p90)
    is_mal_mu2 = is_epi & np.isfinite(cnv_full) & (cnv_full > mu2)
    is_mal_gmm = is_epi & np.isfinite(cnv_full) & (cnv_full > gmm_cut)
    # Sensitivity: CNV-high and not a strong normal-lung program.
    if is_epi.sum() >= 50:
        nl_cut = float(np.quantile(normal_score[is_epi], 0.75))
    else:
        nl_cut = 0.3
    is_mal_p95_nnl = is_mal_p95 & (normal_score < nl_cut)
    is_dip_p95 = is_epi & np.isfinite(cnv_full) & (cnv_full <= p95)

    print(
        f"malignant p95={is_mal_p95.sum()} p90={is_mal_p90.sum()} "
        f"mu2={is_mal_mu2.sum()} gmm={is_mal_gmm.sum()} p95_nnl={is_mal_p95_nnl.sum()} "
        f"diploid_epi={is_dip_p95.sum()}",
        flush=True,
    )

    cell_df = pd.DataFrame(
        {
            "barcode": cells,
            "sample": sample,
            "lineage": lineage,
            "n_umi": n_umi,
            "n_genes": n_genes,
            "TACSTD2": tac,
            "CD3E": cd3e,
            "CD8A": cd8a,
            "NKG7": nkg7,
            "tacstd2_log1p_cp10k": np.log1p(tac * scale),
            "cnv_sumabs": cnv_full,
            "cnv_meanabs": cnv_mean_full,
            "normal_lung_score": normal_score,
            "is_epithelial": is_epi.astype(int),
            "is_stromal": is_stromal.astype(int),
            "is_mal_p95": is_mal_p95.astype(int),
            "is_mal_p90": is_mal_p90.astype(int),
            "is_mal_mu2": is_mal_mu2.astype(int),
            "is_mal_gmm": is_mal_gmm.astype(int),
            "is_mal_p95_nnl": is_mal_p95_nnl.astype(int),
            "is_dip_p95": is_dip_p95.astype(int),
            "is_lineage_tnk": is_lineage_tnk.astype(int),
            "is_umi_tnk": is_umi_tnk.astype(int),
            "is_umi_strict_tnk": is_umi_strict_tnk.astype(int),
            "drmref_celltype": drm_type,
            "is_drm_mal": is_drm_mal.astype(int),
            "is_drm_tnk": is_drm_tnk.astype(int),
        }
    )
    # do not write 92k-row file uncompressed if huge; write a compact gzip
    cell_df.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    smap = meta.set_index("Sample")
    rows = []
    for samp, g in cell_df.groupby("sample", sort=True):
        if samp not in smap.index:
            continue
        rec = smap.loc[samp]
        n_s = len(g)
        n_epi = int(g["is_epithelial"].sum())
        row = {
            "sample": samp,
            "patient": rec["Patient"],
            "resource": rec["Resource"],
            "is_post": bool(rec["is_post"]),
            "response_raw": rec["Pathologic_Response"],
            "response": rec["response"],
            "residual_tumor": rec["residual_tumor"],
            "pathology": rec["Pathology"],
            "n_cells": n_s,
            "n_epithelial": n_epi,
            "frac_epithelial": n_epi / n_s if n_s else np.nan,
            "n_stromal": int(g["is_stromal"].sum()),
            "n_mal_p95": int(g["is_mal_p95"].sum()),
            "n_mal_p90": int(g["is_mal_p90"].sum()),
            "n_mal_mu2": int(g["is_mal_mu2"].sum()),
            "n_mal_gmm": int(g["is_mal_gmm"].sum()),
            "n_mal_p95_nnl": int(g["is_mal_p95_nnl"].sum()),
            "n_dip_p95": int(g["is_dip_p95"].sum()),
            "n_drm_mal": int(g["is_drm_mal"].sum()),
            "n_lineage_tnk": int(g["is_lineage_tnk"].sum()),
            "n_umi_tnk": int(g["is_umi_tnk"].sum()),
            "n_umi_strict_tnk": int(g["is_umi_strict_tnk"].sum()),
            "n_drm_tnk": int(g["is_drm_tnk"].sum()),
            "n_drm_annot": int((g["drmref_celltype"] != "").sum()),
            "frac_lineage_tnk": g["is_lineage_tnk"].mean(),
            "frac_umi_tnk": g["is_umi_tnk"].mean(),
            "frac_umi_strict_tnk": g["is_umi_strict_tnk"].mean(),
            "frac_drm_tnk": (
                g["is_drm_tnk"].sum() / (g["drmref_celltype"] != "").sum()
                if (g["drmref_celltype"] != "").sum()
                else np.nan
            ),
        }
        mal_defs = {
            "p95": g["is_mal_p95"].astype(bool),
            "p90": g["is_mal_p90"].astype(bool),
            "mu2": g["is_mal_mu2"].astype(bool),
            "gmm": g["is_mal_gmm"].astype(bool),
            "p95_nnl": g["is_mal_p95_nnl"].astype(bool),
            "drmref": g["is_drm_mal"].astype(bool),
            "all_epi": g["is_epithelial"].astype(bool),
            "dip_p95": g["is_dip_p95"].astype(bool),
        }
        for name, mask in mal_defs.items():
            sub = g.loc[mask]
            row[f"{name}_n"] = int(len(sub))
            if len(sub) == 0:
                row[f"{name}_mean_log1p_cp10k"] = np.nan
                row[f"{name}_mean_log1p_umi"] = np.nan
                row[f"{name}_pct_pos_ge1"] = np.nan
                row[f"{name}_pct_pos_ge2"] = np.nan
                row[f"{name}_pct_pos_ge3"] = np.nan
                row[f"{name}_mean_cnv"] = np.nan
                continue
            row[f"{name}_mean_log1p_cp10k"] = float(sub["tacstd2_log1p_cp10k"].mean())
            row[f"{name}_mean_log1p_umi"] = float(np.log1p(sub["TACSTD2"]).mean())
            row[f"{name}_pct_pos_ge1"] = float((sub["TACSTD2"] >= 1).mean())
            row[f"{name}_pct_pos_ge2"] = float((sub["TACSTD2"] >= 2).mean())
            row[f"{name}_pct_pos_ge3"] = float((sub["TACSTD2"] >= 3).mean())
            row[f"{name}_mean_cnv"] = float(sub["cnv_sumabs"].mean()) if sub["cnv_sumabs"].notna().any() else np.nan
        rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(args.outdir / "per_sample.tsv", sep="\t", index=False)

    post = per[per["is_post"]].copy()
    tests = []

    mal_metrics = [
        ("p95", "mean_log1p_cp10k"),
        ("p95", "pct_pos_ge1"),
        ("p95", "pct_pos_ge2"),
        ("p95", "pct_pos_ge3"),
        ("p90", "mean_log1p_cp10k"),
        ("p90", "pct_pos_ge1"),
        ("mu2", "mean_log1p_cp10k"),
        ("gmm", "mean_log1p_cp10k"),
        ("p95_nnl", "mean_log1p_cp10k"),
        ("drmref", "mean_log1p_cp10k"),
        ("drmref", "pct_pos_ge1"),
        ("all_epi", "mean_log1p_cp10k"),
        ("dip_p95", "mean_log1p_cp10k"),
    ]
    tnk_defs = [
        ("lineage_TNK", "frac_lineage_tnk"),
        ("umi_CD3E_CD8A_NKG7", "frac_umi_tnk"),
        ("umi_CD3E_or_NKG7notCD3E", "frac_umi_strict_tnk"),
        ("drmref_TNK", "frac_drm_tnk"),
    ]

    def add_group_test(df: pd.DataFrame, metric_col: str, mal_name: str, metric_name: str, min_n: int):
        use = df.copy()
        if min_n > 0:
            use = use[use[f"{mal_name}_n"] >= min_n]
        if use.empty:
            return
        is_nmpr = use["response"].to_numpy() == "NMPR"
        vals = use[metric_col].to_numpy()
        w = exact_wilcoxon(vals, is_nmpr)
        tests.append(
            {
                "family": "NMPR_vs_MPR",
                "malignant_def": mal_name,
                "metric": metric_name,
                "min_mal_cells": min_n,
                "n": int(len(use)),
                "n_NMPR": int(is_nmpr.sum()),
                "n_MPR": int((~is_nmpr).sum()),
                "mean_NMPR": w["mean_a"],
                "mean_MPR": w["mean_b"],
                "delta_NMPR_minus_MPR": w["delta_a_minus_b"],
                "U": w["U"],
                "p": w["p_exact"],
                "p_onesided_NMPR_gt_MPR": w["p_onesided_A_gt_B_mean"],
                "rho": None,
                "tnk_def": "",
                "match_NMPR_gt_MPR": match_nmpr(w["delta_a_minus_b"], w["p_exact"]),
                "match_rho_window": "",
                "note": f"exact enumeration C({w['n_a']+w['n_b']},{w['n_a']})={w['n_perm']}; pCR counted as MPR",
            }
        )

    for mal_name, metric_name in mal_metrics:
        col = f"{mal_name}_{metric_name}"
        add_group_test(post, col, mal_name, metric_name, 0)
        add_group_test(post, col, mal_name, metric_name, MIN_MAL)

    def add_corr(df: pd.DataFrame, metric_col: str, tnk_col: str, mal_name: str, metric_name: str, tnk_name: str, min_n: int, extra=""):
        use = df.copy()
        if min_n > 0:
            use = use[use[f"{mal_name}_n"] >= min_n]
        use = use[np.isfinite(use[metric_col]) & np.isfinite(use[tnk_col])]
        if len(use) < 4:
            return
        s = spearman_safe(use[metric_col], use[tnk_col])
        tests.append(
            {
                "family": "vs_TNK",
                "malignant_def": mal_name,
                "metric": metric_name,
                "min_mal_cells": min_n,
                "n": s["n"],
                "n_NMPR": int((use["response"] == "NMPR").sum()),
                "n_MPR": int((use["response"] == "MPR").sum()),
                "mean_NMPR": None,
                "mean_MPR": None,
                "delta_NMPR_minus_MPR": None,
                "U": None,
                "p": s["p"],
                "p_onesided_NMPR_gt_MPR": None,
                "rho": s["rho"],
                "tnk_def": tnk_name,
                "match_NMPR_gt_MPR": "",
                "match_rho_window": match_rho(s["rho"]),
                "note": extra,
            }
        )
        # residual on epithelial fraction
        res = residualize(use[metric_col].to_numpy(), use["frac_epithelial"].to_numpy())
        sr = spearman_safe(res, use[tnk_col].to_numpy())
        tests.append(
            {
                "family": "vs_TNK_residual_on_epi_frac",
                "malignant_def": mal_name,
                "metric": metric_name,
                "min_mal_cells": min_n,
                "n": sr["n"],
                "n_NMPR": int((use["response"] == "NMPR").sum()),
                "n_MPR": int((use["response"] == "MPR").sum()),
                "mean_NMPR": None,
                "mean_MPR": None,
                "delta_NMPR_minus_MPR": None,
                "U": None,
                "p": sr["p"],
                "p_onesided_NMPR_gt_MPR": None,
                "rho": sr["rho"],
                "tnk_def": tnk_name,
                "match_NMPR_gt_MPR": "",
                "match_rho_window": match_rho(sr["rho"]),
                "note": extra + "; residual = TACSTD2 ~ epithelial_fraction (OLS)",
            }
        )
        # partial Spearman via rank residuals
        ry = rank_residualize(use[metric_col].to_numpy(), use["frac_epithelial"].to_numpy())
        rx = rank_residualize(use[tnk_col].to_numpy(), use["frac_epithelial"].to_numpy())
        sp = spearman_safe(ry, rx)
        tests.append(
            {
                "family": "vs_TNK_partial_epi_frac",
                "malignant_def": mal_name,
                "metric": metric_name,
                "min_mal_cells": min_n,
                "n": sp["n"],
                "n_NMPR": int((use["response"] == "NMPR").sum()),
                "n_MPR": int((use["response"] == "MPR").sum()),
                "mean_NMPR": None,
                "mean_MPR": None,
                "delta_NMPR_minus_MPR": None,
                "U": None,
                "p": sp["p"],
                "p_onesided_NMPR_gt_MPR": None,
                "rho": sp["rho"],
                "tnk_def": tnk_name,
                "match_NMPR_gt_MPR": "",
                "match_rho_window": match_rho(sp["rho"]),
                "note": extra + "; partial Spearman controlling epithelial_fraction (rank residuals)",
            }
        )

    primary_corr_metrics = [
        ("p95", "mean_log1p_cp10k"),
        ("p95", "pct_pos_ge1"),
        ("p95", "pct_pos_ge2"),
        ("p95", "pct_pos_ge3"),
        ("p95_nnl", "mean_log1p_cp10k"),
        ("drmref", "mean_log1p_cp10k"),
        ("all_epi", "mean_log1p_cp10k"),
    ]
    for mal_name, metric_name in primary_corr_metrics:
        col = f"{mal_name}_{metric_name}"
        for tnk_name, tnk_col in tnk_defs:
            add_corr(post, col, tnk_col, mal_name, metric_name, tnk_name, 0)
            add_corr(post, col, tnk_col, mal_name, metric_name, tnk_name, MIN_MAL)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(args.outdir / "stats.tsv", sep="\t", index=False)

    # lineage counts
    lin_counts = (
        cell_df.groupby(["sample", "lineage"]).size().rename("n").reset_index().pivot_table(
            index="sample", columns="lineage", values="n", fill_value=0
        )
    )
    lin_counts.to_csv(args.outdir / "lineage_counts_by_sample.tsv", sep="\t")

    # figures
    fig, ax = plt.subplots(1, 1, figsize=(6.2, 4.2))
    epi_s = cnv_full[is_epi]
    ref_s = cnv_full[is_stromal]
    ax.hist(ref_s[np.isfinite(ref_s)], bins=40, alpha=0.6, label=f"stromal n={np.isfinite(ref_s).sum()}", density=True)
    ax.hist(epi_s[np.isfinite(epi_s)], bins=40, alpha=0.5, label=f"epithelial n={np.isfinite(epi_s).sum()}", density=True)
    if np.isfinite(p95):
        ax.axvline(p95, color="k", ls="--", lw=1, label=f"ref p95={p95:.1f}")
    ax.set_xlabel("CopyKAT-like CNV score (sum |smoothed|)")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    ax.set_title("GSE207422 epithelial vs stromal CNV scores")
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_cnv_score.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    pp = post[post["p95_n"] >= 0].copy()
    colors = {"MPR": "#2ca02c", "NMPR": "#d62728"}
    ax = axes[0]
    for resp, sub in pp.groupby("response"):
        ax.scatter(
            [resp] * len(sub),
            sub["p95_mean_log1p_cp10k"],
            c=colors.get(resp, "gray"),
            s=40,
            alpha=0.85,
        )
    ax.set_ylabel("malignant (CNV p95) mean log1p(CP10k) TACSTD2")
    ax.set_title("NMPR vs MPR")
    ax = axes[1]
    for resp, sub in pp.groupby("response"):
        ax.scatter(
            sub["frac_lineage_tnk"],
            sub["p95_mean_log1p_cp10k"],
            c=colors.get(resp, "gray"),
            s=40,
            label=resp,
        )
        for _, r in sub.iterrows():
            if np.isfinite(r["p95_mean_log1p_cp10k"]):
                ax.annotate(r["patient"], (r["frac_lineage_tnk"], r["p95_mean_log1p_cp10k"]), fontsize=7)
    ax.set_xlabel("T/NK fraction (lineage)")
    ax.set_ylabel("malignant TACSTD2")
    ax.legend(fontsize=8)
    ax.set_title("vs T/NK")
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_tacstd2_nmpr_tnk.png", dpi=140)
    plt.close(fig)

    # primary rows for JSON
    def pick(family, mal, metric, tnk="", min_n=0):
        q = tests_df[
            (tests_df["family"] == family)
            & (tests_df["malignant_def"] == mal)
            & (tests_df["metric"] == metric)
            & (tests_df["min_mal_cells"] == min_n)
        ]
        if tnk:
            q = q[q["tnk_def"] == tnk]
        return q.iloc[0].to_dict() if len(q) else {}

    primary = {
        "dataset": "GSE207422",
        "n_cells_public_matrix": int(n),
        "n_genes_matrix": int(n_matrix_genes),
        "user_claimed_cells": 90652,
        "n_epithelial": int(is_epi.sum()),
        "n_stromal": int(is_stromal.sum()),
        "n_malignant_p95": int(is_mal_p95.sum()),
        "n_malignant_drmref": int(is_drm_mal.sum()),
        "author_copykat_ids_public": False,
        "cnv": cnv_info,
        "primary_NMPR_vs_MPR": pick("NMPR_vs_MPR", "p95", "mean_log1p_cp10k", min_n=0),
        "primary_vs_TNK_lineage": pick("vs_TNK", "p95", "mean_log1p_cp10k", "lineage_TNK", 0),
        "primary_vs_TNK_umi": pick("vs_TNK", "p95", "mean_log1p_cp10k", "umi_CD3E_CD8A_NKG7", 0),
        "primary_residual_lineage": pick(
            "vs_TNK_residual_on_epi_frac", "p95", "mean_log1p_cp10k", "lineage_TNK", 0
        ),
    }
    with open(args.outdir / "summary.json", "w") as fh:
        json.dump(primary, fh, indent=2, default=str)

    with open(args.outdir / "cnv_info.json", "w") as fh:
        json.dump(cnv_info, fh, indent=2, default=str)

    print("wrote", args.outdir, "n_tests", len(tests_df), flush=True)


if __name__ == "__main__":
    main()
