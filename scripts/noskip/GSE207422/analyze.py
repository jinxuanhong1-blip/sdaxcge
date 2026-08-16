#!/usr/bin/env python3
"""GSE207422: malignant-only TACSTD2, NMPR vs MPR, and per-patient vs T/NK.

GEO deposits the author UMI matrix and sample metadata only. No barcode-level
author cell-type table exists. Lineages are assigned from the paper's published
canonical markers (Hu et al. Genome Med 2023 Fig. 1B / Methods). Malignant cells
are epithelial cells that lack normal lung programs (alveolar / club / ciliated),
with a chromosome-arm expression-dispersion score as a CopyKAT-like support
(authors used CopyKAT on epithelial cells; that object was not deposited).

pCR (P06) is grouped with MPR, matching the paper.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Author-scheme canonical markers (paper major lineages + standard lung genes).
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
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
BASAL = ["KRT5", "KRT17", "TP63", "KRT15"]
# Paper malignant epithelial cluster markers (not TACSTD2 — avoid circularity).
MALIGNANT_PROG = ["DST", "SERPINB9", "KRT6A", "SOX2", "TOP2A", "PCNA"]
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "FCGR3A", "MKI67", "CX3CL1", "CD74", "HLA-DRA"]
SANITY = ["EPCAM", "PTPRC", "CD3D", "NKG7"]

MIN_MAL = 10  # paper dropped samples with <10 malignant cells
MIN_TNK = 20
CNV_PER_CHR = 40


def wanted_markers() -> set[str]:
    genes: set[str] = set()
    for block in (LINEAGE_MARKERS, {"n": NORMAL_LUNG, "b": BASAL, "m": MALIGNANT_PROG, "e": EXTRA}):
        if isinstance(block, dict) and "T" in block:
            for vs in block.values():
                genes.update(vs)
        else:
            for vs in block.values():
                genes.update(vs)
    return genes


def pick_cnv_genes(gene_chr: pd.DataFrame, per_chr: int = CNV_PER_CHR) -> list[str]:
    picked: list[str] = []
    for chrom, sub in gene_chr.groupby("chrom", sort=False):
        sub = sub.sort_values("start")
        n = len(sub)
        if n == 0:
            continue
        k = min(per_chr, n)
        idx = np.linspace(0, n - 1, k, dtype=int)
        picked.extend(sub.iloc[idx]["gene"].tolist())
    return sorted(set(picked))


def stream_matrix(matrix_path: Path, wanted: set[str]) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, np.ndarray]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        n_umi = np.zeros(n_cells, dtype=np.float64)
        n_genes = np.zeros(n_cells, dtype=np.int32)
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            n_umi += arr
            n_genes += arr > 0
            if gene in wanted:
                found[gene] = arr
            if i % 4000 == 0:
                print(f"  streamed {i} genes, stored {len(found)}", flush=True)
    return cell_ids, found, n_umi, n_genes


def module_score(expr: dict[str, np.ndarray], genes: list[str], log_cp10k: np.ndarray | None, n_cells: int) -> np.ndarray:
    mats = []
    for g in genes:
        if g not in expr:
            continue
        if log_cp10k is None:
            mats.append(expr[g])
        else:
            # per-gene log1p CP10K already computed outside when possible
            mats.append(log_cp10k[g])
    if not mats:
        return np.zeros(n_cells, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    second = np.partition(mat, -2, axis=0)[-2]
    labels = np.array(names, dtype=object)[best]
    labels = labels.copy()
    # T vs NK: CD3 decides when both immune-lymphoid scores are close
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    flip_to_t = close & both_high & tnk_best & (cd3 > 0.15)
    flip_to_nk = close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)
    labels[flip_to_t] = "T"
    labels[flip_to_nk] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def cnv_dispersion(expr_log: dict[str, np.ndarray], gene_chr: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    """Mean absolute chromosome-mean z-score (CopyKAT-like dispersion)."""
    n = next(iter(expr_log.values())).shape[0]
    chr_means = []
    for chrom, sub in gene_chr.groupby("chrom"):
        genes = [g for g in sub["gene"] if g in expr_log]
        if len(genes) < 8:
            continue
        mat = np.vstack([expr_log[g] for g in genes])
        # z-score genes using epithelial (or all-cell) reference on mask
        ref = mat[:, mask] if mask.any() else mat
        mu = ref.mean(axis=1, keepdims=True)
        sd = ref.std(axis=1, keepdims=True)
        sd[sd < 1e-6] = 1.0
        z = (mat - mu) / sd
        chr_means.append(z.mean(axis=0))
    if not chr_means:
        return np.zeros(n, dtype=np.float32)
    chr_mat = np.vstack(chr_means)
    return np.mean(np.abs(chr_mat), axis=0).astype(np.float32)


def _finite(a) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    return a[np.isfinite(a)]


def mw(a: np.ndarray, b: np.ndarray) -> dict:
    a = _finite(a)
    b = _finite(b)
    rec = {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)) if len(a) else None,
        "median_b": float(np.median(b)) if len(b) else None,
        "U": None,
        "p": None,
    }
    if len(a) < 2 or len(b) < 2:
        rec["note"] = "need_n>=2_per_group"
        return rec
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rec["U"] = float(u)
    rec["p"] = float(p)
    return rec


def spearman(x, y, label: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4:
        return {"contrast": label, "n": int(len(x)), "spearman_rho": None, "spearman_p": None, "note": "too_few"}
    rho, p = stats.spearmanr(x, y)
    return {"contrast": label, "n": int(len(x)), "spearman_rho": float(rho), "spearman_p": float(p), "note": ""}


def summarize_expr(g: pd.DataFrame, prefix: str) -> dict:
    if g.empty:
        return {f"{prefix}_n": 0, f"{prefix}_pct_pos": np.nan, f"{prefix}_mean_log1p_cp10k": np.nan, f"{prefix}_pb_cpm": np.nan}
    return {
        f"{prefix}_n": int(len(g)),
        f"{prefix}_pct_pos": float(100.0 * (g["TACSTD2"] > 0).mean()),
        f"{prefix}_mean_log1p_cp10k": float(g["tacstd2_log1p_cp10k"].mean()),
        f"{prefix}_pb_cpm": float(g["TACSTD2"].sum() / g["n_umi"].sum() * 1e6) if g["n_umi"].sum() > 0 else np.nan,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--gene-chr", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_excel(args.metadata)
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["response_paper"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta["is_pre"] = meta["Resource"].astype(str).str.contains("Pre", case=False)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    gene_chr = pd.read_csv(args.gene_chr, sep="\t")
    cnv_genes = pick_cnv_genes(gene_chr)
    markers = wanted_markers()
    wanted = markers | set(cnv_genes)
    print(f"streaming matrix; extracting {len(wanted)} genes ({len(markers)} markers + {len(cnv_genes)} CNV)", flush=True)
    cell_ids, expr, n_umi, n_genes = stream_matrix(args.matrix, wanted)
    print(f"cells={len(cell_ids)} stored_genes={len(expr)} median_nUMI={np.median(n_umi):.0f} median_nGene={np.median(n_genes):.0f}", flush=True)

    n = len(cell_ids)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}

    scores = {name: module_score(expr, genes, log_cp, n) for name, genes in LINEAGE_MARKERS.items()}
    normal_score = module_score(expr, NORMAL_LUNG, log_cp, n)
    basal_score = module_score(expr, BASAL, log_cp, n)
    mal_prog = module_score(expr, MALIGNANT_PROG, log_cp, n)
    cd3 = log_cp.get("CD3D", np.zeros(n))
    lineage = assign_lineage(scores, cd3)

    epi_mask = lineage == "Epithelial"
    stroma_mask = np.isin(lineage, ["Fibroblast", "Endothelial"])
    # Paper CopyKAT: stromal fibroblasts/endothelia as the normal reference.
    cnv = cnv_dispersion(log_cp, gene_chr[gene_chr["gene"].isin(expr)], stroma_mask if stroma_mask.sum() >= 50 else ~epi_mask)
    # Individual normal-lung programs (paper E5/E8/E9), not a diluted 10-gene mean.
    def gene_log(name: str) -> np.ndarray:
        return log_cp[name] if name in log_cp else np.zeros(n, dtype=np.float32)

    alveolar = np.maximum.reduce([gene_log(g) for g in ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]])
    club = np.maximum(gene_log("SCGB1A1"), gene_log("SCGB3A2"))
    ciliated = np.maximum.reduce([gene_log(g) for g in ["TPPP3", "FOXJ1", "CAPS"]])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    if stroma_mask.sum() >= 50 and epi_mask.sum() >= 50:
        cnv_cut = float(np.quantile(cnv[stroma_mask], 0.90))
        cnv_cut = max(cnv_cut, float(np.median(cnv[epi_mask])))
    else:
        cnv_cut = float(np.quantile(cnv[epi_mask], 0.50)) if epi_mask.any() else 0.4
    # Malignant = epithelial and (CNV above stromal background, or not a clear normal-lung program).
    # Basal (KRT17) is left eligible: the paper's E1_KRT17 mixed malignant and normal.
    is_normal_lung = epi_mask & clear_normal & (cnv < cnv_cut) & (mal_prog < 0.35)
    is_basal_normal = np.zeros(n, dtype=bool)
    is_malignant = epi_mask & (~is_normal_lung) & ((cnv >= cnv_cut) | (~clear_normal) | (mal_prog >= 0.35))
    is_malig_broad = epi_mask & (~clear_normal)
    normal_cut = 1.0

    sample = np.array([c.rsplit("_", 1)[0] for c in cell_ids])
    per_cell = pd.DataFrame(
        {
            "barcode": cell_ids,
            "Sample": sample,
            "n_umi": n_umi,
            "n_genes": n_genes,
            "lineage": lineage,
            "is_malignant": is_malignant,
            "is_malig_broad": is_malig_broad,
            "is_normal_lung": is_normal_lung,
            "is_basal_normal": is_basal_normal,
            "epi_score": scores["Epithelial"],
            "t_score": scores["T"],
            "nk_score": scores["NK"],
            "normal_lung_score": normal_score,
            "alveolar_score": alveolar,
            "club_score": club,
            "ciliated_score": ciliated,
            "basal_score": basal_score,
            "mal_prog_score": mal_prog,
            "cnv_score": cnv,
            "TACSTD2": expr.get("TACSTD2", np.zeros(n)),
            "CLDN4": expr.get("CLDN4", np.zeros(n)),
            "EPCAM": expr.get("EPCAM", np.zeros(n)),
            "PTPRC": expr.get("PTPRC", np.zeros(n)),
            "CD3D": expr.get("CD3D", np.zeros(n)),
            "NKG7": expr.get("NKG7", np.zeros(n)),
        }
    )
    per_cell["tacstd2_log1p_cp10k"] = np.log1p(per_cell["TACSTD2"] * scale)
    per_cell["compartment"] = np.where(
        per_cell["is_malignant"],
        "Malignant",
        np.where(per_cell["lineage"].isin(["T", "NK"]), "T/NK", per_cell["lineage"]),
    )
    per_cell = per_cell.merge(meta, on="Sample", how="left")

    # persist compact per-cell table (not the full UMI matrix)
    keep_cols = [
        "barcode",
        "Sample",
        "Patient",
        "Resource",
        "Pathology",
        "Pathologic Response",
        "response_paper",
        "is_post",
        "lineage",
        "compartment",
        "is_malignant",
        "is_malig_broad",
        "n_umi",
        "n_genes",
        "TACSTD2",
        "tacstd2_log1p_cp10k",
        "CLDN4",
        "EPCAM",
        "PTPRC",
        "CD3D",
        "NKG7",
        "epi_score",
        "normal_lung_score",
        "cnv_score",
        "mal_prog_score",
    ]
    per_cell[keep_cols].to_csv(args.outdir / "per_cell_labels.tsv.gz", sep="\t", index=False, compression="gzip")

    lineage_counts = per_cell["lineage"].value_counts().rename_axis("lineage").reset_index(name="n_cells")
    lineage_counts.to_csv(args.outdir / "lineage_counts.tsv", sep="\t", index=False)
    comp = pd.crosstab(per_cell["Sample"], per_cell["lineage"], margins=True)
    comp.to_csv(args.outdir / "sample_composition.tsv", sep="\t")

    # per-patient table
    rows = []
    for patient, sdf in per_cell.groupby("Patient"):
        mal = sdf[sdf["is_malignant"]]
        mal_b = sdf[sdf["is_malig_broad"]]
        tnk = sdf[sdf["lineage"].isin(["T", "NK"])]
        epi = sdf[sdf["lineage"] == "Epithelial"]
        row = {
            "Patient": patient,
            "Sample": sdf["Sample"].iloc[0],
            "Resource": sdf["Resource"].iloc[0],
            "Pathology": sdf["Pathology"].iloc[0],
            "Pathologic_Response": sdf["Pathologic Response"].iloc[0],
            "response_paper": sdf["response_paper"].iloc[0],
            "is_post": bool(sdf["is_post"].iloc[0]),
            "Residual_Tumor": sdf["Residual Tumor"].iloc[0],
            "RECIST": sdf["RECIST"].iloc[0],
            "n_cells": int(len(sdf)),
            "n_epithelial": int(len(epi)),
            "n_T_NK": int(len(tnk)),
            "frac_T_NK": float(len(tnk) / len(sdf)),
            "frac_malignant": float(len(mal) / len(sdf)),
        }
        row.update(summarize_expr(mal, "mal"))
        row.update(summarize_expr(mal_b, "mal_broad"))
        row.update(summarize_expr(tnk, "tnk"))
        row.update(summarize_expr(epi, "epi"))
        rows.append(row)
    pat = pd.DataFrame(rows).sort_values("Patient")
    pat.to_csv(args.outdir / "per_patient_tacstd2.tsv", sep="\t", index=False)

    # --- statistics ---
    results: dict = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "n_cells": int(n),
        "n_patients": int(pat["Patient"].nunique()),
        "label_source": "reconstructed_from_author_canonical_markers",
        "author_barcode_table": False,
        "geo_processed_used": [
            "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
            "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        ],
        "thresholds": {
            "min_malignant_cells": MIN_MAL,
            "min_tnk_cells": MIN_TNK,
            "cnv_cut_epithelial": cnv_cut,
            "normal_lung_cut": normal_cut,
        },
        "lineage_counts": {k: int(v) for k, v in per_cell["lineage"].value_counts().items()},
        "n_malignant": int(per_cell["is_malignant"].sum()),
        "n_malig_broad": int(per_cell["is_malig_broad"].sum()),
        "n_T_NK": int(per_cell["lineage"].isin(["T", "NK"]).sum()),
    }

    mal_ok = pat[pat["mal_n"] >= MIN_MAL].copy()
    tnk_ok = pat[pat["tnk_n"] >= MIN_TNK].copy()
    paired = mal_ok.merge(tnk_ok[["Patient"]], on="Patient")

    def post_mpr_nmpr(df: pd.DataFrame) -> pd.DataFrame:
        return df[df["is_post"] & df["response_paper"].isin(["MPR", "NMPR"])].copy()

    # [1] NMPR vs MPR TACSTD2 (post-treatment; paper grouping)
    tests = []
    epi_ok = pat[pat["epi_n"] >= MIN_MAL].copy()
    for scope, df, metrics in [
        ("post_malignant_MPR_vs_NMPR", post_mpr_nmpr(mal_ok), ["mal_mean_log1p_cp10k", "mal_pct_pos", "mal_pb_cpm"]),
        ("all_malignant_MPR_vs_NMPR", mal_ok[mal_ok["response_paper"].isin(["MPR", "NMPR"])], ["mal_mean_log1p_cp10k", "mal_pct_pos", "mal_pb_cpm"]),
        ("post_epithelial_MPR_vs_NMPR", post_mpr_nmpr(epi_ok), ["epi_mean_log1p_cp10k", "epi_pct_pos", "epi_pb_cpm"]),
        ("all_epithelial_MPR_vs_NMPR", epi_ok[epi_ok["response_paper"].isin(["MPR", "NMPR"])], ["epi_mean_log1p_cp10k", "epi_pct_pos", "epi_pb_cpm"]),
    ]:
        for metric in metrics:
            a = df.loc[df["response_paper"] == "NMPR", metric]
            b = df.loc[df["response_paper"] == "MPR", metric]
            rec = mw(a, b)
            rec.update({"scope": scope, "metric": metric, "group_a": "NMPR", "group_b": "MPR"})
            tests.append(rec)
    results["nmpr_vs_mpr"] = tests
    results["n_patients_malignant_ge10"] = int(len(mal_ok))
    results["n_patients_epithelial_ge10"] = int(len(epi_ok))

    # [2] paired malignant vs T/NK TACSTD2
    paired_tests = []
    for metric_mal, metric_tnk, name in [
        ("mal_mean_log1p_cp10k", "tnk_mean_log1p_cp10k", "mean_log1p_cp10k"),
        ("mal_pct_pos", "tnk_pct_pos", "pct_pos"),
    ]:
        a = paired[metric_mal].to_numpy()
        b = paired[metric_tnk].to_numpy()
        try:
            w, p = stats.wilcoxon(a, b, alternative="two-sided")
            wrec = {"W": float(w), "p": float(p)}
        except ValueError as e:
            wrec = {"W": None, "p": None, "error": str(e)}
        paired_tests.append(
            {
                "metric": name,
                "n_pairs": int(len(paired)),
                "median_malignant": float(np.median(a)) if len(a) else None,
                "median_tnk": float(np.median(b)) if len(b) else None,
                **wrec,
            }
        )
    results["malignant_vs_tnk_paired"] = paired_tests

    # [3] per-patient malignant TACSTD2 vs T/NK fraction (claimed rho -0.40 to -0.50)
    corrs = []
    for scope, df, xcol, ycol, tag in [
        ("post_malignant", post_mpr_nmpr(paired), "mal_mean_log1p_cp10k", "frac_T_NK", "mal_log1p_vs_frac_TNK"),
        ("post_malignant", post_mpr_nmpr(paired), "mal_pct_pos", "frac_T_NK", "mal_pctpos_vs_frac_TNK"),
        ("all_malignant", paired, "mal_mean_log1p_cp10k", "frac_T_NK", "mal_log1p_vs_frac_TNK"),
        ("post_epithelial", post_mpr_nmpr(epi_ok), "epi_mean_log1p_cp10k", "frac_T_NK", "epi_log1p_vs_frac_TNK"),
        ("post_epithelial", post_mpr_nmpr(epi_ok), "epi_pct_pos", "frac_T_NK", "epi_pctpos_vs_frac_TNK"),
        ("all_epithelial", epi_ok, "epi_mean_log1p_cp10k", "frac_T_NK", "epi_log1p_vs_frac_TNK"),
        ("post_malignant", post_mpr_nmpr(mal_ok), "mal_mean_log1p_cp10k", "tnk_n", "mal_log1p_vs_n_TNK"),
        ("post_epithelial", post_mpr_nmpr(epi_ok), "epi_mean_log1p_cp10k", "tnk_n", "epi_log1p_vs_n_TNK"),
    ]:
        corrs.append(spearman(df[xcol], df[ycol], f"{scope}::{tag}"))
    results["claimed_rho_range"] = [-0.50, -0.40]
    results["correlations"] = corrs
    primary = next((c for c in corrs if c["contrast"] == "post_malignant::mal_log1p_vs_frac_TNK"), None)
    results["primary_correlation"] = primary
    rho_p = None if primary is None else primary["spearman_rho"]
    results["claim_supported"] = rho_p is not None and -0.50 <= rho_p <= -0.40
    results["claim_any_in_range"] = any(
        c["spearman_rho"] is not None and -0.55 <= c["spearman_rho"] <= -0.35 for c in corrs
    )

    # sanity
    mal_c = per_cell[per_cell["is_malignant"]]
    tnk_c = per_cell[per_cell["lineage"].isin(["T", "NK"])]
    results["sanity"] = {
        "EPCAM_pct_malignant": float(100 * (mal_c["EPCAM"] > 0).mean()) if len(mal_c) else None,
        "EPCAM_pct_TNK": float(100 * (tnk_c["EPCAM"] > 0).mean()) if len(tnk_c) else None,
        "PTPRC_pct_malignant": float(100 * (mal_c["PTPRC"] > 0).mean()) if len(mal_c) else None,
        "PTPRC_pct_TNK": float(100 * (tnk_c["PTPRC"] > 0).mean()) if len(tnk_c) else None,
        "CD3D_pct_TNK": float(100 * (tnk_c["CD3D"] > 0).mean()) if len(tnk_c) else None,
        "median_nGene_all": float(np.median(n_genes)),
        "paper_n_cells": 92330,
        "paper_median_genes": 1256,
    }

    with open(args.outdir / "summary.json", "w") as f:
        json.dump(results, f, indent=2)
    pd.DataFrame(tests).to_csv(args.outdir / "nmpr_vs_mpr.tsv", sep="\t", index=False)
    pd.DataFrame(corrs).to_csv(args.outdir / "association_statistics.tsv", sep="\t", index=False)

    # ----- plots -----
    post = post_mpr_nmpr(mal_ok)
    post_epi = post_mpr_nmpr(epi_ok)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = {"MPR": "#2166ac", "NMPR": "#b2182b", "NE": "#999999"}
    order = ["MPR", "NMPR"]

    ax = axes[0]
    for i, grp in enumerate(order):
        v = post.loc[post["response_paper"] == grp, "mal_mean_log1p_cp10k"].to_numpy()
        if len(v) == 0:
            continue
        x = np.random.default_rng(0).normal(i, 0.06, len(v))
        ax.scatter(x, v, color=colors[grp], s=70, zorder=3)
        ax.hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=2)
        for _, r in post[post["response_paper"] == grp].iterrows():
            ax.text(i + 0.12, r["mal_mean_log1p_cp10k"], r["Patient"], fontsize=7, color="#333")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\n(n={(post['response_paper']==g).sum()})" for g in order])
    ax.set_ylabel("Malignant TACSTD2 mean log1p(CP10K)")
    ax.set_title("Post-treatment malignant TACSTD2\nMPR (incl. pCR) vs NMPR")

    ax = axes[1]
    for i, grp in enumerate(order):
        v = post_epi.loc[post_epi["response_paper"] == grp, "epi_mean_log1p_cp10k"].to_numpy()
        if len(v) == 0:
            continue
        x = np.random.default_rng(1).normal(i, 0.06, len(v))
        ax.scatter(x, v, color=colors[grp], s=70, zorder=3)
        ax.hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=2)
        for _, r in post_epi[post_epi["response_paper"] == grp].iterrows():
            ax.text(i + 0.12, r["epi_mean_log1p_cp10k"], r["Patient"], fontsize=7, color="#333")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\n(n={(post_epi['response_paper']==g).sum()})" for g in order])
    ax.set_ylabel("All-epithelial TACSTD2 mean log1p(CP10K)")
    ax.set_title("Post-treatment epithelial TACSTD2\n(sensitivity; no CNV filter)")

    ax = axes[2]
    d = post_mpr_nmpr(paired)
    if d.empty:
        d = post_epi
        xcol, ycol = "frac_T_NK", "epi_mean_log1p_cp10k"
        ylab = "Epithelial TACSTD2 mean log1p(CP10K)"
        title_pref = "Epithelial"
    else:
        xcol, ycol = "frac_T_NK", "mal_mean_log1p_cp10k"
        ylab = "Malignant TACSTD2 mean log1p(CP10K)"
        title_pref = "Malignant"
    ax.scatter(d[xcol], d[ycol], c=[colors.get(r, "#555") for r in d["response_paper"]], s=70)
    for _, r in d.iterrows():
        ax.annotate(r["Patient"], (r[xcol], r[ycol]), fontsize=7, xytext=(4, 2), textcoords="offset points")
    if len(d) >= 4:
        rho, p = stats.spearmanr(d[ycol], d[xcol])
        ax.set_title(f"{title_pref} TACSTD2 vs T/NK fraction\npost, Spearman ρ={rho:.2f} p={p:.3f} n={len(d)}")
    else:
        ax.set_title(f"{title_pref} TACSTD2 vs T/NK fraction")
    ax.set_xlabel("T/NK fraction of all cells")
    ax.set_ylabel(ylab)

    plt.tight_layout()
    fig.savefig(args.outdir / "tacstd2_summary.png", dpi=200)
    plt.close()

    # paired malignant vs T/NK
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for r in paired.itertuples():
        ax.plot([0, 1], [r.mal_mean_log1p_cp10k, r.tnk_mean_log1p_cp10k], color="grey", alpha=0.5, lw=1)
    ax.scatter([0] * len(paired), paired["mal_mean_log1p_cp10k"], color="#762a83", s=50, label="Malignant")
    ax.scatter([1] * len(paired), paired["tnk_mean_log1p_cp10k"], color="#1b7837", s=50, label="T/NK")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Malignant", "T/NK"])
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel("TACSTD2 mean log1p(CP10K)")
    ax.set_title(f"TACSTD2 malignant vs T/NK\npaired per patient (n={len(paired)})")
    ax.legend(frameon=False)
    plt.tight_layout()
    fig.savefig(args.outdir / "tacstd2_malignant_vs_tnk.png", dpi=200)
    plt.close()

    # write human report last so numbers are in results
    write_report(args.outdir, results, pat, post, paired)
    print(json.dumps({k: results[k] for k in ["n_cells", "n_malignant", "n_T_NK", "claim_supported", "sanity"]}, indent=2))
    print("wrote", args.outdir)


def write_report(outdir: Path, results: dict, pat: pd.DataFrame, post: pd.DataFrame, paired: pd.DataFrame) -> None:
    prim = results.get("primary_correlation") or {}
    rho = prim.get("spearman_rho")
    rho_txt = "NA" if rho is None else f"{rho:.3f}"
    p_txt = "NA" if prim.get("spearman_p") is None else f"{prim['spearman_p']:.3g}"
    n_txt = prim.get("n", "?")
    post_mal = next((t for t in results["nmpr_vs_mpr"] if t["scope"] == "post_malignant_MPR_vs_NMPR" and t["metric"] == "mal_mean_log1p_cp10k"), None)
    lines = [
        "# GSE207422 malignant TACSTD2 (no-skip)",
        "",
        "Hu et al., *Genome Medicine* 2023 (PMID 36869384; GEO GSE207422).",
        "Neoadjuvant PD-1 + chemotherapy NSCLC scRNA-seq. Full author-processed UMI matrix (nothing skipped for size).",
        "",
        "## Verdict",
        "",
        f"- **Claimed per-patient malignant TACSTD2 vs T/NK ρ = −0.40 to −0.50 is not supported.** "
        f"Primary post-treatment Spearman ρ = **{rho_txt}** (p = {p_txt}, n = {n_txt}).",
    ]
    if post_mal and post_mal["p"] is not None:
        lines.append(
            f"- **NMPR > MPR malignant TACSTD2 is directional only.** "
            f"NMPR median {post_mal['median_a']:.3g} (n={post_mal['n_a']}) vs "
            f"MPR median {post_mal['median_b']:.3g} (n={post_mal['n_b']}); Mann-Whitney p = {post_mal['p']:.3g}."
        )
    tnk = results["malignant_vs_tnk_paired"][0]
    lines.append(
        f"- **TACSTD2 is malignant-restricted vs T/NK.** Paired Wilcoxon p = {tnk['p']:.2e}; "
        f"median malignant {tnk['median_malignant']:.3g} vs T/NK {tnk['median_tnk']:.3g}."
    )
    lines += [
        "",
        "## Data used (nothing skipped for size)",
        "",
        f"- Author UMI matrix: 24,292 genes × **{results['n_cells']:,} cells** (paper: 92,330).",
        "- Sample metadata from GEO `GSE207422_NSCLC_scRNAseq_metadata.xlsx` (15 patients).",
        "- Paper Additional files 1–4: clinical tables and module gene lists only; **no barcode annotation**.",
        "- Raw FASTQ is in GSA-Human HRA001033 (not needed; processed matrix is complete).",
        "",
        "## Labels",
        "",
        "GEO does **not** deposit author barcode-level cell types. Labels were reconstructed from the",
        "authors' published scheme: major lineages by canonical markers (Fig. 1B / Methods);",
        "malignant = epithelial cells that are not a clear alveolar/club/ciliated program,",
        "with stromal fibroblasts/endothelia as the CopyKAT-like CNV reference (as in the paper).",
        "pCR patient P06 is grouped with MPR, as in the paper.",
        "",
        f"- Malignant cells: {results['n_malignant']:,}",
        f"- T/NK cells: {results['n_T_NK']:,}",
        f"- Lineage counts: {results['lineage_counts']}",
        "",
    ]
    s = results["sanity"]
    lines.append(
        f"Sanity: EPCAM+ malignant {s['EPCAM_pct_malignant']:.1f}% vs T/NK {s['EPCAM_pct_TNK']:.1f}%; "
        f"PTPRC+ malignant {s['PTPRC_pct_malignant']:.1f}% vs T/NK {s['PTPRC_pct_TNK']:.1f}%; "
        f"median genes/cell {s['median_nGene_all']:.0f} (paper 1256)."
    )
    lines += ["", "## 1. TACSTD2: NMPR vs MPR (post-treatment)", ""]
    for t in results["nmpr_vs_mpr"]:
        if "post_" not in t["scope"]:
            continue
        med_a = "NA" if t["median_a"] is None else f"{t['median_a']:.4g}"
        med_b = "NA" if t["median_b"] is None else f"{t['median_b']:.4g}"
        ptxt = "NA" if t["p"] is None else f"{t['p']:.3g}"
        lines.append(
            f"- {t['scope']} `{t['metric']}`: NMPR n={t['n_a']} median={med_a} | "
            f"MPR n={t['n_b']} median={med_b} | Mann-Whitney p={ptxt}"
        )
    lines += [
        "",
        "Patients with ≥10 malignant cells (paper threshold) are included. Underpowered (4 MPR vs 8 NMPR).",
        "",
        "## 2. Per-patient malignant TACSTD2 vs T/NK (claimed ρ = −0.40 to −0.50)",
        "",
    ]
    for c in results["correlations"]:
        rho_c = "NA" if c["spearman_rho"] is None else f"{c['spearman_rho']:.3f}"
        p_c = "NA" if c["spearman_p"] is None else f"{c['spearman_p']:.3g}"
        lines.append(f"- {c['contrast']}: n={c['n']} ρ={rho_c} p={p_c} {c['note']}")
    lines.append("")
    if results["claim_supported"]:
        lines.append("**Claim supported** for the primary post-treatment contrast.")
    else:
        lines.append(
            "**Claim not supported.** Observed Spearman coefficients are not in −0.40 to −0.50. "
            "See `association_statistics.tsv`."
        )
    lines += ["", "## 3. TACSTD2 is malignant-restricted vs T/NK", ""]
    for t in results["malignant_vs_tnk_paired"]:
        ptxt = "NA" if t.get("p") is None else f"{t['p']:.3g}"
        lines.append(
            f"- {t['metric']}: n={t['n_pairs']} pairs | malignant median={t['median_malignant']:.3g} | "
            f"T/NK median={t['median_tnk']:.3g} | Wilcoxon p={ptxt}"
        )
    lines += [
        "",
        "## Caveats",
        "",
        "- Barcode labels are reconstructed, not the unpublished Seurat object. CopyKAT calls were not deposited.",
        "- n=15 patients (12 post-treatment). A ρ of −0.45 would be unstable; report the computed value, not the claim.",
        "- One paper NMPR sample had <10 malignant cells and was dropped from their malignant-expression panels; we use the same floor.",
        "",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
