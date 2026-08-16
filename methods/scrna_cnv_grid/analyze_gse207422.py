#!/usr/bin/env python3
"""GSE207422 malignant-definition grid: author (none public), marker, CNV, intersections.

Hu et al. Genome Medicine 2023 (PMID 36869384). Public GEO UMI only.
Author CopyKAT barcode IDs are not deposited. DRMref labels are third-party.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from grid_common import (  # noqa: E402
    DRMREF_TNK,
    EXTRA,
    LINEAGE_MARKERS,
    MIN_MAL,
    NORMAL_LUNG,
    WIN,
    assign_lineage,
    copykat_like,
    direction_nmpr,
    direction_rho,
    exact_wilcoxon,
    kmeans1d_cut,
    module_score,
    recovers_user_direction,
    score_compartment,
    spearman_safe,
    wanted_markers,
)


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
    ap.add_argument("--xlsx", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--gene-chr", type=Path, default=HERE / "gene_chr.tsv")
    ap.add_argument("--drmref", type=Path, default=HERE / "drmref_cell_annotation.tsv.gz")
    ap.add_argument("--outdir", type=Path, default=Path("results/scrna_cnv_grid/GSE207422"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

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
    normal_umi = np.zeros(n, dtype=np.float32)
    for g in NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]

    is_epi = lineage == "Epithelial"
    is_stromal = np.isin(lineage, ["Fibroblast", "Endothelial"])
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    is_lineage_tnk = is_t | is_nk
    is_umi_tnk = (cd3e >= 1) | (cd8a >= 1) | (nkg7 >= 1)

    drm_type = np.array([""] * n, dtype=object)
    if args.drmref.exists():
        drm = pd.read_csv(args.drmref, sep="\t")
        dmap = dict(zip(drm["cell_barcode"].astype(str), drm["celltype"].astype(str)))
        drm_type = np.array([dmap.get(c, "") for c in cells], dtype=object)
    is_drm_mal = drm_type == "Malignant cells"
    is_drm_tnk = np.isin(drm_type, list(DRMREF_TNK))

    print(
        f"lineage epi={is_epi.sum()} stromal={is_stromal.sum()} T={is_t.sum()} NK={is_nk.sum()} drm_mal={is_drm_mal.sum()}",
        flush=True,
    )

    rng = np.random.default_rng(0)
    immune_idx = np.where(is_t | (lineage == "Myeloid"))[0]
    if immune_idx.size > 2500:
        immune_idx = rng.choice(immune_idx, size=2500, replace=False)
    keep = np.zeros(n, dtype=bool)
    keep[is_epi | is_stromal] = True
    keep[immune_idx] = True
    keep_idx = np.where(keep)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    print(f"pass2 keep_cells={keep_idx.size}", flush=True)

    cnv_genes, cnv_mat = stream_pass2(args.matrix, keep_idx, cnv_gene_set)
    ref_keep = is_stromal[keep_idx]
    if ref_keep.sum() < 80:
        print(f"stromal ref n={ref_keep.sum()} < 80; adding T+myeloid subsample", flush=True)
        ref_keep = ref_keep | np.isin(keep_idx, immune_idx)
    score_sumabs, score_meanabs, cnv_info = copykat_like(cnv_mat, cnv_genes, gene_chr, ref_keep, WIN)
    print("CNV info", json.dumps({k: v for k, v in cnv_info.items() if k != "chromosomes"}), flush=True)

    cnv_full = np.full(n, np.nan, dtype=np.float32)
    for i, j in keep_map.items():
        cnv_full[i] = score_sumabs[j]

    ref_scores = score_sumabs[ref_keep]
    ref_scores = ref_scores[np.isfinite(ref_scores)]
    p95 = float(np.quantile(ref_scores, 0.95)) if ref_scores.size else np.nan
    p90 = float(np.quantile(ref_scores, 0.90)) if ref_scores.size else np.nan
    mu2 = float(ref_scores.mean() + 2 * ref_scores.std()) if ref_scores.size else np.nan
    gmm_cut = kmeans1d_cut(cnv_full[is_epi])
    cnv_info.update({"threshold_p95": p95, "threshold_p90": p90, "threshold_mean2sd": mu2, "threshold_gmm": gmm_cut})

    if is_epi.sum() >= 50:
        nl_cut = float(np.quantile(normal_score[is_epi], 0.75))
    else:
        nl_cut = 0.3

    is_marker_nnl_zero = is_epi & (normal_umi == 0)
    is_marker_nnl_p75 = is_epi & (normal_score < nl_cut)
    is_cnv_p95 = is_epi & np.isfinite(cnv_full) & (cnv_full > p95)
    is_cnv_p90 = is_epi & np.isfinite(cnv_full) & (cnv_full > p90)
    is_cnv_mu2 = is_epi & np.isfinite(cnv_full) & (cnv_full > mu2)
    is_cnv_gmm = is_epi & np.isfinite(cnv_full) & (cnv_full > gmm_cut)

    # Grid definitions (author CopyKAT not public).
    mal_flags = {
        "thirdparty_drmref": is_drm_mal,
        "marker_all_epi": is_epi,
        "marker_epi_nnl_zero": is_marker_nnl_zero,
        "marker_epi_nnl_p75": is_marker_nnl_p75,
        "cnv_p95": is_cnv_p95,
        "cnv_p90": is_cnv_p90,
        "cnv_mu2": is_cnv_mu2,
        "cnv_gmm": is_cnv_gmm,
        "intersect_markerNNL_cnvP95": is_marker_nnl_p75 & is_cnv_p95,
        "intersect_drmref_cnvP95": is_drm_mal & is_cnv_p95,
        "intersect_drmref_markerNNL": is_drm_mal & is_marker_nnl_p75,
        "intersect_drmref_markerNNL_cnvP95": is_drm_mal & is_marker_nnl_p75 & is_cnv_p95,
        "intersect_drmref_allEpi": is_drm_mal & is_epi,
    }
    print("malignant n:", {k: int(v.sum()) for k, v in mal_flags.items()}, flush=True)

    cell_df = pd.DataFrame(
        {
            "barcode": cells,
            "sample": sample,
            "lineage": lineage,
            "n_umi": n_umi,
            "n_genes": n_genes,
            "TACSTD2": tac,
            "tacstd2_log1p_cp10k": np.log1p(tac * scale),
            "cnv_sumabs": cnv_full,
            "normal_lung_score": normal_score,
            "drmref_celltype": drm_type,
            "is_lineage_tnk": is_lineage_tnk.astype(int),
            "is_umi_tnk": is_umi_tnk.astype(int),
            "is_drm_tnk": is_drm_tnk.astype(int),
        }
    )
    for name, mask in mal_flags.items():
        cell_df[f"mal_{name}"] = mask.astype(int)
    cell_df.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    smap = meta.set_index("Sample")
    rows = []
    for samp, g in cell_df.groupby("sample", sort=True):
        if samp not in smap.index:
            continue
        rec = smap.loc[samp]
        n_s = len(g)
        n_drm_annot = int((g["drmref_celltype"] != "").sum())
        row = {
            "dataset": "GSE207422",
            "sample": samp,
            "patient": rec["Patient"],
            "resource": rec["Resource"],
            "is_post": bool(rec["is_post"]),
            "response_raw": rec["Pathologic_Response"],
            "response": rec["response"],
            "residual_tumor": rec["residual_tumor"],
            "n_cells": n_s,
            "frac_lineage_tnk": float(g["is_lineage_tnk"].mean()),
            "frac_umi_tnk": float(g["is_umi_tnk"].mean()),
            "frac_drm_tnk": (float(g["is_drm_tnk"].sum() / n_drm_annot) if n_drm_annot else np.nan),
            "n_drm_annot": n_drm_annot,
        }
        for name in mal_flags:
            sc = score_compartment(g.loc[g[f"mal_{name}"].astype(bool)])
            for k, v in sc.items():
                row[f"{name}_{k}"] = v
        rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(args.outdir / "per_sample.tsv", sep="\t", index=False)

    post = per[per["is_post"]].copy()
    tnk_defs = [
        ("lineage_TNK", "frac_lineage_tnk"),
        ("umi_CD3E_CD8A_NKG7", "frac_umi_tnk"),
        ("drmref_TNK", "frac_drm_tnk"),
    ]
    metrics = ["mean_log1p_cp10k", "pct_pos_ge1", "pseudobulk_cpm"]
    tests = []

    def add_pair(df, mal_name, metric_name, tnk_name, tnk_col, min_n):
        use = df.copy()
        if min_n > 0:
            use = use[use[f"{mal_name}_n"] >= min_n]
        mcol = f"{mal_name}_{metric_name}"
        use = use[np.isfinite(use[mcol])]
        if use.empty:
            return
        is_nmpr = use["response"].to_numpy() == "NMPR"
        w = exact_wilcoxon(use[mcol].to_numpy(), is_nmpr)
        s = {"n": None, "rho": None, "p": None}
        use_c = use[np.isfinite(use[tnk_col])]
        if len(use_c) >= 4:
            s = spearman_safe(use_c[mcol], use_c[tnk_col])
        tests.append(
            {
                "dataset": "GSE207422",
                "cohort": "post_n12",
                "malignant_def": mal_name,
                "metric": metric_name,
                "tnk_def": tnk_name,
                "min_mal_cells": min_n,
                "n": int(len(use)),
                "n_corr": s["n"],
                "n_NMPR": int(is_nmpr.sum()),
                "n_MPR": int((~is_nmpr).sum()),
                "mean_NMPR": w["mean_a"],
                "mean_MPR": w["mean_b"],
                "delta_NMPR_minus_MPR": w["delta_a_minus_b"],
                "p_nmpr": w["p_exact"],
                "p_onesided_NMPR_gt_MPR": w["p_onesided_A_gt_B_mean"],
                "p_method": w.get("p_method"),
                "rho": s["rho"],
                "p_rho": s["p"],
                "dir_nmpr": direction_nmpr(w["delta_a_minus_b"]),
                "dir_rho": direction_rho(s["rho"]),
                "recovers_user_direction": recovers_user_direction(w["delta_a_minus_b"], s["rho"]),
                "note": "pCR counted as MPR; author CopyKAT IDs not public; drmref is third-party",
            }
        )

    for mal_name in mal_flags:
        for metric_name in metrics:
            for tnk_name, tnk_col in tnk_defs:
                add_pair(post, mal_name, metric_name, tnk_name, tnk_col, 0)
                add_pair(post, mal_name, metric_name, tnk_name, tnk_col, MIN_MAL)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(args.outdir / "grid_stats.tsv", sep="\t", index=False)

    # Primary slice: mean log1p CP10k, min_n=0, each T/NK def.
    primary = tests_df[
        (tests_df["metric"] == "mean_log1p_cp10k") & (tests_df["min_mal_cells"] == 0)
    ].copy()
    primary.to_csv(args.outdir / "direction_table.tsv", sep="\t", index=False)

    both = primary[primary["recovers_user_direction"] == "both"]
    print("defs recovering both directions (any T/NK):", sorted(both["malignant_def"].unique()), flush=True)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0))
    colors = {"MPR": "#2ca02c", "NMPR": "#d62728"}
    ax = axes[0]
    for resp, sub in post.groupby("response"):
        ax.scatter([resp] * len(sub), sub["cnv_p95_mean_log1p_cp10k"], c=colors.get(resp, "gray"), s=40)
    ax.set_ylabel("CNV-p95 malignant mean log1p(CP10k) TACSTD2")
    ax.set_title("GSE207422 NMPR vs MPR")
    ax = axes[1]
    for resp, sub in post.groupby("response"):
        ax.scatter(sub["frac_lineage_tnk"], sub["cnv_p95_mean_log1p_cp10k"], c=colors.get(resp, "gray"), s=40, label=resp)
        for _, r in sub.iterrows():
            if np.isfinite(r["cnv_p95_mean_log1p_cp10k"]):
                ax.annotate(r["patient"], (r["frac_lineage_tnk"], r["cnv_p95_mean_log1p_cp10k"]), fontsize=7)
    ax.set_xlabel("T/NK fraction (lineage)")
    ax.legend(fontsize=8)
    ax.set_title("vs T/NK")
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_grid_cnv_p95.png", dpi=140)
    plt.close(fig)

    # heatmap-like bar of direction recoveries
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    sub = primary[primary["tnk_def"] == "lineage_TNK"].copy()
    y = np.arange(len(sub))
    ax.barh(y, sub["delta_NMPR_minus_MPR"].fillna(0), color=["#d62728" if d > 0 else "#1f77b4" for d in sub["delta_NMPR_minus_MPR"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels(sub["malignant_def"], fontsize=8)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Δ mean log1p(CP10k) TACSTD2 (NMPR − MPR)")
    ax.set_title("GSE207422 definition grid — NMPR vs MPR")
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_grid_nmpr_delta.png", dpi=140)
    plt.close(fig)

    summary = {
        "dataset": "GSE207422",
        "n_cells": int(n),
        "n_genes_matrix": int(n_matrix_genes),
        "author_copykat_public": False,
        "thirdparty_drmref_malignant": int(is_drm_mal.sum()),
        "n_epithelial": int(is_epi.sum()),
        "n_cnv_p95": int(is_cnv_p95.sum()),
        "cnv": {k: v for k, v in cnv_info.items() if k != "chromosomes"},
        "n_tests": int(len(tests_df)),
        "defs_both_directions": sorted(both["malignant_def"].unique().tolist()),
        "primary_rows": int(len(primary)),
    }
    with open(args.outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    with open(args.outdir / "cnv_info.json", "w") as fh:
        json.dump(cnv_info, fh, indent=2, default=str)
    print("wrote", args.outdir, "n_tests", len(tests_df), flush=True)


if __name__ == "__main__":
    main()
