#!/usr/bin/env python3
"""GSE241934 malignant-definition grid: author Epi, marker epithelium, CNV, intersections.

Zhang et al. Cell Reports Medicine 2024 (NEOTIDE/CTONG2104; GEO GSE241934).
Public MTX + author major.cell.type and Pathological Response.
Author Epi (IIT n=1,699) is the deposited epithelial compartment after the paper's
CopyKAT filter (paper text: N=1,669 malignant epithelial).
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


def group_of(x) -> str:
    x = str(x)
    if x in {"non-MPR", "NMPR", "nonMPR"}:
        return "NMPR"
    if x in {"MPR", "pCR"}:
        return "MPR"
    return ""


def load_features(path: Path) -> list[str]:
    genes = []
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # prefer gene symbol in col1 if present
            if len(parts) >= 2 and parts[1] and parts[1] != parts[0]:
                genes.append(parts[1])
            else:
                genes.append(parts[0])
    return genes


def stream_mtx_markers(
    mtx_path: Path,
    n_cells: int,
    marker_rows: dict[int, str],
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    marker_arr = {g: np.zeros(n_cells, dtype=np.float32) for g in marker_rows.values()}
    n_umi = np.zeros(n_cells, dtype=np.float64)
    n_nz = 0
    with gzip.open(mtx_path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            _nr, _nc, _nnz = line.split()
            break
        for line in f:
            r_s, c_s, v_s = line.split()
            r = int(r_s) - 1
            c = int(c_s) - 1
            v = float(v_s)
            n_umi[c] += v
            n_nz += 1
            if r in marker_rows:
                marker_arr[marker_rows[r]][c] = v
            if n_nz % 8_000_000 == 0:
                print(f"  marker-pass nz={n_nz}", flush=True)
    print(f"marker-pass done nz={n_nz} genes={len(marker_arr)}", flush=True)
    return marker_arr, n_umi


def stream_mtx_cnv(
    mtx_path: Path,
    cnv_rows: dict[int, str],
    keep_cells: np.ndarray,
) -> tuple[list[str], np.ndarray]:
    keep_idx = np.where(keep_cells)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    n_keep = keep_idx.size
    cnv_names = list(dict.fromkeys(cnv_rows.values()))
    name_to_i = {g: i for i, g in enumerate(cnv_names)}
    cnv_mat = np.zeros((len(cnv_names), n_keep), dtype=np.float32)
    n_nz = 0
    with gzip.open(mtx_path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            _nr, _nc, _nnz = line.split()
            break
        for line in f:
            r_s, c_s, v_s = line.split()
            r = int(r_s) - 1
            c = int(c_s) - 1
            n_nz += 1
            if r in cnv_rows and c in keep_map:
                cnv_mat[name_to_i[cnv_rows[r]], keep_map[c]] = float(v_s)
            if n_nz % 8_000_000 == 0:
                print(f"  cnv-pass nz={n_nz}", flush=True)
    print(f"cnv-pass done nz={n_nz} genes={len(cnv_names)} keep={n_keep}", flush=True)
    return cnv_names, cnv_mat


def analyze_cohort(
    name: str,
    mtx: Path,
    features: Path,
    barcodes: Path,
    meta: Path,
    gene_chr: pd.DataFrame,
    outdir: Path,
    max_stromal_ref: int = 2500,
) -> pd.DataFrame:
    print(f"=== {name} ===", flush=True)
    md = pd.read_csv(meta, sep="\t", dtype=str, low_memory=False)
    bc = pd.read_csv(barcodes, sep="\t", header=None)
    genes = load_features(features)
    n_cells = len(bc)
    print(f"meta={md.shape} barcodes={n_cells} features={len(genes)}", flush=True)
    if len(md) != n_cells:
        raise RuntimeError(f"{name}: meta {len(md)} != barcodes {n_cells}")
    md = md.reset_index(drop=True)
    md["barcode"] = bc.iloc[:, 0].astype(str).to_numpy()

    major = md["major.cell.type"].astype(str)
    is_author_epi = major.eq("Epi").to_numpy()
    is_author_t = major.eq("T").to_numpy()
    is_author_nk = major.eq("NK").to_numpy()
    is_author_tnk = is_author_t | is_author_nk
    is_author_fibro = major.eq("Fibro").to_numpy()
    is_author_endo = major.eq("Endo").to_numpy()
    is_author_stromal = is_author_fibro | is_author_endo
    sample = md["sampleID"].astype(str).to_numpy()
    response = md["Pathological Response"].map(group_of).to_numpy()

    wanted = wanted_markers()
    cnv_set = set(gene_chr["gene"])
    marker_rows = {}
    cnv_rows = {}
    for i, g in enumerate(genes):
        if g in wanted:
            marker_rows[i] = g
        if g in cnv_set:
            cnv_rows[i] = g
    print(f"marker_rows={len(marker_rows)} cnv_rows={len(cnv_rows)} author_epi={is_author_epi.sum()} stromal={is_author_stromal.sum()}", flush=True)

    rng = np.random.default_rng(0)
    stromal_idx = np.where(is_author_stromal)[0]
    if stromal_idx.size > max_stromal_ref:
        stromal_keep = np.zeros(n_cells, dtype=bool)
        stromal_keep[rng.choice(stromal_idx, size=max_stromal_ref, replace=False)] = True
    else:
        stromal_keep = is_author_stromal.copy()

    print("marker pass", flush=True)
    marker_arr, n_umi = stream_mtx_markers(mtx, n_cells, marker_rows)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(marker_arr[g] * scale).astype(np.float32) for g in marker_arr}
    scores = {ln: module_score(log_cp, gs, n_cells) for ln, gs in LINEAGE_MARKERS.items()}
    normal_score = module_score(log_cp, NORMAL_LUNG, n_cells)
    lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n_cells))))
    tac = marker_arr.get("TACSTD2", np.zeros(n_cells, dtype=np.float32))
    cd3e = marker_arr.get("CD3E", np.zeros(n_cells, dtype=np.float32))
    cd8a = marker_arr.get("CD8A", np.zeros(n_cells, dtype=np.float32))
    nkg7 = marker_arr.get("NKG7", np.zeros(n_cells, dtype=np.float32))
    normal_umi = np.zeros(n_cells, dtype=np.float32)
    for g in NORMAL_LUNG:
        if g in marker_arr:
            normal_umi += marker_arr[g]

    is_marker_epi = lineage == "Epithelial"
    is_marker_stromal = np.isin(lineage, ["Fibroblast", "Endothelial"])
    is_lineage_tnk = np.isin(lineage, ["T", "NK"])
    is_umi_tnk = (cd3e >= 1) | (cd8a >= 1) | (nkg7 >= 1)
    print(
        f"marker epi={is_marker_epi.sum()} stromal={is_marker_stromal.sum()} "
        f"TNK_lineage={is_lineage_tnk.sum()} author_TNK={is_author_tnk.sum()}",
        flush=True,
    )

    keep_cells = is_author_epi | is_marker_epi | stromal_keep
    print(f"cnv keep_cells={keep_cells.sum()} (authorEpi+markerEpi+stromal_ref)", flush=True)
    print("cnv pass", flush=True)
    cnv_names, cnv_mat = stream_mtx_cnv(mtx, cnv_rows, keep_cells)
    keep_idx = np.where(keep_cells)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    # Prefer author stromal as diploid reference; fall back to marker stromal.
    ref_keep = stromal_keep[keep_idx]
    if ref_keep.sum() < 80:
        ref_keep = is_marker_stromal[keep_idx]
    score_sumabs, _, cnv_info = copykat_like(cnv_mat, cnv_names, gene_chr, ref_keep, WIN)
    print("CNV", json.dumps({k: v for k, v in cnv_info.items() if k != "chromosomes"}), flush=True)
    cnv_full = np.full(n_cells, np.nan, dtype=np.float32)
    for i, j in keep_map.items():
        cnv_full[i] = score_sumabs[j]

    ref_scores = score_sumabs[ref_keep]
    ref_scores = ref_scores[np.isfinite(ref_scores)]
    p95 = float(np.quantile(ref_scores, 0.95)) if ref_scores.size else np.nan
    p90 = float(np.quantile(ref_scores, 0.90)) if ref_scores.size else np.nan
    gmm_cut = kmeans1d_cut(cnv_full[is_author_epi | is_marker_epi])
    cnv_info.update({"threshold_p95": p95, "threshold_p90": p90, "threshold_gmm": gmm_cut, "cohort": name})

    epi_for_nl = is_author_epi | is_marker_epi
    nl_cut = float(np.quantile(normal_score[epi_for_nl], 0.75)) if epi_for_nl.sum() >= 50 else 0.3
    is_marker_nnl = is_marker_epi & (normal_score < nl_cut)
    is_marker_nnl_zero = is_marker_epi & (normal_umi == 0)
    is_cnv_p95_marker = is_marker_epi & np.isfinite(cnv_full) & (cnv_full > p95)
    is_cnv_p95_author = is_author_epi & np.isfinite(cnv_full) & (cnv_full > p95)
    is_cnv_p90_author = is_author_epi & np.isfinite(cnv_full) & (cnv_full > p90)
    is_cnv_gmm_author = is_author_epi & np.isfinite(cnv_full) & (cnv_full > gmm_cut)

    mal_flags = {
        "author_epi": is_author_epi,
        "marker_all_epi": is_marker_epi,
        "marker_epi_nnl_zero": is_marker_nnl_zero,
        "marker_epi_nnl_p75": is_marker_nnl,
        "cnv_p95_on_markerEpi": is_cnv_p95_marker,
        "cnv_p95_on_authorEpi": is_cnv_p95_author,
        "cnv_p90_on_authorEpi": is_cnv_p90_author,
        "cnv_gmm_on_authorEpi": is_cnv_gmm_author,
        "intersect_authorEpi_markerEpi": is_author_epi & is_marker_epi,
        "intersect_authorEpi_markerNNL": is_author_epi & is_marker_nnl,
        "intersect_authorEpi_cnvP95": is_cnv_p95_author,
        "intersect_markerNNL_cnvP95": is_marker_nnl & np.isfinite(cnv_full) & (cnv_full > p95),
        "intersect_author_markerNNL_cnvP95": is_author_epi & is_marker_nnl & np.isfinite(cnv_full) & (cnv_full > p95),
    }
    print("malignant n:", {k: int(v.sum()) for k, v in mal_flags.items()}, flush=True)

    cell_df = pd.DataFrame(
        {
            "barcode": md["barcode"],
            "sample": sample,
            "response": response,
            "lineage": lineage,
            "author_major": major,
            "n_umi": n_umi,
            "TACSTD2": tac,
            "tacstd2_log1p_cp10k": np.log1p(tac * scale),
            "cnv_sumabs": cnv_full,
            "normal_lung_score": normal_score,
            "is_author_tnk": is_author_tnk.astype(int),
            "is_lineage_tnk": is_lineage_tnk.astype(int),
            "is_umi_tnk": is_umi_tnk.astype(int),
        }
    )
    for nm, mask in mal_flags.items():
        cell_df[f"mal_{nm}"] = mask.astype(int)
    cell_df.to_csv(outdir / f"{name}_cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    rows = []
    for samp, g in cell_df.groupby("sample", sort=True):
        resp = str(g["response"].iloc[0])
        if resp not in {"MPR", "NMPR"}:
            continue
        n_s = len(g)
        row = {
            "dataset": "GSE241934",
            "cohort": name,
            "sample": samp,
            "patient": samp,
            "response": resp,
            "n_cells": n_s,
            "frac_author_tnk": float(g["is_author_tnk"].mean()),
            "frac_lineage_tnk": float(g["is_lineage_tnk"].mean()),
            "frac_umi_tnk": float(g["is_umi_tnk"].mean()),
        }
        for nm in mal_flags:
            sc = score_compartment(g.loc[g[f"mal_{nm}"].astype(bool)])
            for k, v in sc.items():
                row[f"{nm}_{k}"] = v
        rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(outdir / f"{name}_per_sample.tsv", sep="\t", index=False)

    tnk_defs = [
        ("author_TNK", "frac_author_tnk"),
        ("lineage_TNK", "frac_lineage_tnk"),
        ("umi_CD3E_CD8A_NKG7", "frac_umi_tnk"),
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
                "dataset": "GSE241934",
                "cohort": name,
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
                "note": "pCR counted as MPR; author_epi is deposited major.cell.type==Epi",
            }
        )

    for mal_name in mal_flags:
        for metric_name in metrics:
            for tnk_name, tnk_col in tnk_defs:
                add_pair(per, mal_name, metric_name, tnk_name, tnk_col, 0)
                add_pair(per, mal_name, metric_name, tnk_name, tnk_col, MIN_MAL)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(outdir / f"{name}_grid_stats.tsv", sep="\t", index=False)
    with open(outdir / f"{name}_cnv_info.json", "w") as fh:
        json.dump(cnv_info, fh, indent=2, default=str)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0))
    colors = {"MPR": "#2ca02c", "NMPR": "#d62728"}
    ax = axes[0]
    for resp, sub in per.groupby("response"):
        ax.scatter([resp] * len(sub), sub["author_epi_mean_log1p_cp10k"], c=colors.get(resp, "gray"), s=40)
    ax.set_ylabel("author Epi mean log1p(CP10k) TACSTD2")
    ax.set_title(f"GSE241934 {name} NMPR vs MPR")
    ax = axes[1]
    for resp, sub in per.groupby("response"):
        ax.scatter(sub["frac_author_tnk"], sub["author_epi_mean_log1p_cp10k"], c=colors.get(resp, "gray"), s=40, label=resp)
        for _, r in sub.iterrows():
            if np.isfinite(r["author_epi_mean_log1p_cp10k"]):
                ax.annotate(r["patient"], (r["frac_author_tnk"], r["author_epi_mean_log1p_cp10k"]), fontsize=7)
    ax.set_xlabel("T/NK fraction (author)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_fig_author_epi.png", dpi=140)
    plt.close(fig)

    print(f"wrote {name} n_samples={len(per)} n_tests={len(tests_df)}", flush=True)
    return tests_df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE241934"))
    ap.add_argument("--gene-chr", type=Path, default=HERE / "gene_chr.tsv")
    ap.add_argument("--outdir", type=Path, default=Path("results/scrna_cnv_grid/GSE241934"))
    ap.add_argument("--skip-real", action="store_true")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_chr = pd.read_csv(args.gene_chr, sep="\t")
    gene_chr["chrom"] = gene_chr["chrom"].astype(str)
    d = args.datadir

    frames = []
    frames.append(
        analyze_cohort(
            "IIT_EGFRmut",
            d / "GSE241934_IIT_Matrix.mtx.gz",
            d / "GSE241934_IIT_features.tsv.gz",
            d / "GSE241934_IIT_barcodes.tsv.gz",
            d / "GSE241934_IIT_Meta.txt.gz",
            gene_chr,
            args.outdir,
        )
    )
    if not args.skip_real and (d / "GSE241934_Real_Matrix.mtx.gz").exists():
        frames.append(
            analyze_cohort(
                "REAL_EGFRWT",
                d / "GSE241934_Real_Matrix.mtx.gz",
                d / "GSE241934_RWC_features.tsv.gz",
                d / "GSE241934_RWC_barcodes.tsv.gz",
                d / "GSE241934_Real_Meta.txt.gz",
                gene_chr,
                args.outdir,
                max_stromal_ref=2000,
            )
        )
    all_t = pd.concat(frames, ignore_index=True)
    all_t.to_csv(args.outdir / "grid_stats.tsv", sep="\t", index=False)
    primary = all_t[(all_t["metric"] == "mean_log1p_cp10k") & (all_t["min_mal_cells"] == 0)]
    primary.to_csv(args.outdir / "direction_table.tsv", sep="\t", index=False)
    both = primary[primary["recovers_user_direction"] == "both"]
    with open(args.outdir / "summary.json", "w") as fh:
        json.dump(
            {
                "dataset": "GSE241934",
                "defs_both_directions": sorted(both["malignant_def"].unique().tolist()),
                "n_tests": int(len(all_t)),
            },
            fh,
            indent=2,
        )
    print("GSE241934 defs recovering both directions:", sorted(both["malignant_def"].unique()), flush=True)


if __name__ == "__main__":
    main()
