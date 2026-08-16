#!/usr/bin/env python3
"""GSE207422 two-way extra: malignant TACSTD2 and CLDN4 in pre/post × MPR/NMPR.

Public GEO UMI matrix + sample xlsx only. Lineages from Hu et al. Genome Med 2023
canonical markers. Malignant = epithelial minus clear alveolar/club/ciliated
programs, with a CopyKAT-like chromosome-dispersion support (author CopyKAT
barcodes are not deposited). pCR is grouped with MPR.

Patient is the unit of n. Empty 2×2 cells are kept and named.
"""

from __future__ import annotations

import argparse
import gzip
import json
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
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
BASAL = ["KRT5", "KRT17", "TP63", "KRT15"]
MALIGNANT_PROG = ["DST", "SERPINB9", "KRT6A", "SOX2", "TOP2A", "PCNA"]
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "FCGR3A", "MKI67"]
MIN_MAL = 10
CNV_PER_CHR = 40
TIME_LEVELS = ["pre", "post"]
MPR_LEVELS = ["MPR", "NMPR"]
RECIST_LEVELS = ["PR", "SD"]


def wanted_markers() -> set[str]:
    genes: set[str] = set()
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    genes.update(BASAL)
    genes.update(MALIGNANT_PROG)
    genes.update(EXTRA)
    return genes


def pick_cnv_genes(gene_chr: pd.DataFrame, per_chr: int = CNV_PER_CHR) -> list[str]:
    picked: list[str] = []
    for _, sub in gene_chr.groupby("chrom", sort=False):
        sub = sub.sort_values("start")
        n = len(sub)
        if n == 0:
            continue
        k = min(per_chr, n)
        idx = np.linspace(0, n - 1, k, dtype=int)
        picked.extend(sub.iloc[idx]["gene"].tolist())
    return sorted(set(picked))


def stream_matrix(matrix_path: Path, wanted: set[str]):
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


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n_cells: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n_cells, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best]
    labels = labels.copy()
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
    n = next(iter(expr_log.values())).shape[0]
    chr_means = []
    for _, sub in gene_chr.groupby("chrom"):
        genes = [g for g in sub["gene"] if g in expr_log]
        if len(genes) < 8:
            continue
        mat = np.vstack([expr_log[g] for g in genes])
        ref = mat[:, mask] if mask.any() else mat
        mu = ref.mean(axis=1, keepdims=True)
        sd = ref.std(axis=1, keepdims=True)
        sd[sd < 1e-6] = 1.0
        z = (mat - mu) / sd
        chr_means.append(z.mean(axis=0))
    if not chr_means:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.abs(np.vstack(chr_means)), axis=0).astype(np.float32)


def _finite(a) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    return a[np.isfinite(a)]


def mw(a, b) -> dict:
    a = _finite(a)
    b = _finite(b)
    rec = {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)) if len(a) else None,
        "median_b": float(np.median(b)) if len(b) else None,
        "mean_a": float(np.mean(a)) if len(a) else None,
        "mean_b": float(np.mean(b)) if len(b) else None,
        "U": None,
        "p": None,
        "note": "",
    }
    if len(a) < 2 or len(b) < 2:
        rec["note"] = "need_n>=2_per_group"
        return rec
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rec["U"] = float(u)
    rec["p"] = float(p)
    return rec


def summarize_gene(g: pd.DataFrame, gene: str, prefix: str) -> dict:
    umi = g[gene]
    logv = g[f"{gene.lower()}_log1p_cp10k"]
    if g.empty:
        return {
            f"{prefix}_n": 0,
            f"{prefix}_pct_pos": np.nan,
            f"{prefix}_mean_log1p_cp10k": np.nan,
            f"{prefix}_pb_cpm": np.nan,
        }
    return {
        f"{prefix}_n": int(len(g)),
        f"{prefix}_pct_pos": float(100.0 * (umi > 0).mean()),
        f"{prefix}_mean_log1p_cp10k": float(logv.mean()),
        f"{prefix}_pb_cpm": float(umi.sum() / g["n_umi"].sum() * 1e6) if g["n_umi"].sum() > 0 else np.nan,
    }


def two_by_two(pat: pd.DataFrame, time_col: str, resp_col: str, time_levels, resp_levels, metric: str) -> pd.DataFrame:
    rows = []
    for t in time_levels:
        for r in resp_levels:
            sub = pat[(pat[time_col] == t) & (pat[resp_col] == r)]
            vals = _finite(sub[metric]) if metric in sub else np.array([])
            rows.append(
                {
                    "timepoint": t,
                    "response": r,
                    "n_patients": int(len(sub)),
                    "n_patients_with_metric": int(len(vals)),
                    "patients": ",".join(sub["Patient"].astype(str)) if len(sub) else "",
                    "n_malignant_cells": int(sub["mal_n"].sum()) if len(sub) else 0,
                    "mean": float(np.mean(vals)) if len(vals) else None,
                    "median": float(np.median(vals)) if len(vals) else None,
                    "empty": int(len(vals) == 0),
                }
            )
    return pd.DataFrame(rows)


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
    meta["mpr"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["timepoint"] = np.where(meta["Resource"].astype(str).str.contains("Pre", case=False), "pre", "post")
    meta["recist"] = meta["RECIST"].astype(str)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    gene_chr = pd.read_csv(args.gene_chr, sep="\t")
    cnv_genes = pick_cnv_genes(gene_chr)
    wanted = wanted_markers() | set(cnv_genes)
    print(f"streaming matrix; extracting {len(wanted)} genes", flush=True)
    cell_ids, expr, n_umi, n_genes = stream_matrix(args.matrix, wanted)
    print(f"cells={len(cell_ids)} stored={len(expr)} median_nUMI={np.median(n_umi):.0f}", flush=True)

    n = len(cell_ids)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    mal_prog = module_score(log_cp, MALIGNANT_PROG, n)
    cd3 = log_cp.get("CD3D", np.zeros(n, dtype=np.float32))
    lineage = assign_lineage(scores, cd3)

    epi_mask = lineage == "Epithelial"
    stroma_mask = np.isin(lineage, ["Fibroblast", "Endothelial"])
    cnv = cnv_dispersion(log_cp, gene_chr[gene_chr["gene"].isin(expr)], stroma_mask if stroma_mask.sum() >= 50 else ~epi_mask)

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
    is_normal_lung = epi_mask & clear_normal & (cnv < cnv_cut) & (mal_prog < 0.35)
    is_malignant = epi_mask & (~is_normal_lung) & ((cnv >= cnv_cut) | (~clear_normal) | (mal_prog >= 0.35))
    is_malig_broad = epi_mask & (~clear_normal)

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
            "TACSTD2": expr.get("TACSTD2", np.zeros(n)),
            "CLDN4": expr.get("CLDN4", np.zeros(n)),
        }
    )
    per_cell["tacstd2_log1p_cp10k"] = np.log1p(per_cell["TACSTD2"] * scale)
    per_cell["cldn4_log1p_cp10k"] = np.log1p(per_cell["CLDN4"] * scale)
    per_cell = per_cell.merge(meta, on="Sample", how="left")
    per_cell[["barcode", "Sample", "Patient", "timepoint", "mpr", "recist", "lineage", "is_malignant", "is_malig_broad", "n_umi", "TACSTD2", "CLDN4", "tacstd2_log1p_cp10k", "cldn4_log1p_cp10k"]].to_csv(
        args.outdir / "per_cell_labels.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    rows = []
    for patient, sdf in per_cell.groupby("Patient"):
        mal = sdf[sdf["is_malignant"]]
        mal_b = sdf[sdf["is_malig_broad"]]
        row = {
            "Patient": patient,
            "Sample": sdf["Sample"].iloc[0],
            "timepoint": sdf["timepoint"].iloc[0],
            "Pathologic_Response": sdf["Pathologic Response"].iloc[0],
            "mpr": sdf["mpr"].iloc[0],
            "recist": sdf["recist"].iloc[0],
            "Pathology": sdf["Pathology"].iloc[0],
            "Residual_Tumor": sdf["Residual Tumor"].iloc[0],
            "n_cells": int(len(sdf)),
            "n_epithelial": int((sdf["lineage"] == "Epithelial").sum()),
        }
        row.update(summarize_gene(mal, "TACSTD2", "mal_tacstd2"))
        row.update(summarize_gene(mal, "CLDN4", "mal_cldn4"))
        row.update(summarize_gene(mal_b, "TACSTD2", "malb_tacstd2"))
        row.update(summarize_gene(mal_b, "CLDN4", "malb_cldn4"))
        row["mal_n"] = row["mal_tacstd2_n"]
        row["malb_n"] = row["malb_tacstd2_n"]
        rows.append(row)
    pat = pd.DataFrame(rows).sort_values("Patient")
    pat.to_csv(args.outdir / "per_patient.tsv", sep="\t", index=False)

    mal_ok = pat[pat["mal_n"] >= MIN_MAL].copy()
    malb_ok = pat[pat["malb_n"] >= MIN_MAL].copy()

    metrics = [
        ("mal_tacstd2_mean_log1p_cp10k", "TACSTD2", "mean_log1p_cp10k", "cnv_or_program"),
        ("mal_tacstd2_pct_pos", "TACSTD2", "pct_pos", "cnv_or_program"),
        ("mal_cldn4_mean_log1p_cp10k", "CLDN4", "mean_log1p_cp10k", "cnv_or_program"),
        ("mal_cldn4_pct_pos", "CLDN4", "pct_pos", "cnv_or_program"),
        ("malb_tacstd2_mean_log1p_cp10k", "TACSTD2", "mean_log1p_cp10k", "epi_minus_normal"),
        ("malb_cldn4_mean_log1p_cp10k", "CLDN4", "mean_log1p_cp10k", "epi_minus_normal"),
    ]

    grid_frames = []
    for metric, gene, kind, mal_def in metrics:
        src = mal_ok if mal_def == "cnv_or_program" else malb_ok
        g = two_by_two(src, "timepoint", "mpr", TIME_LEVELS, MPR_LEVELS, metric)
        g.insert(0, "response_axis", "MPR")
        g.insert(0, "metric", metric)
        g.insert(0, "gene", gene)
        g.insert(0, "malignant_def", mal_def)
        grid_frames.append(g)
        g2 = two_by_two(src, "timepoint", "recist", TIME_LEVELS, RECIST_LEVELS, metric)
        g2.insert(0, "response_axis", "RECIST")
        g2.insert(0, "metric", metric)
        g2.insert(0, "gene", gene)
        g2.insert(0, "malignant_def", mal_def)
        grid_frames.append(g2)
    grid = pd.concat(grid_frames, ignore_index=True)
    grid.to_csv(args.outdir / "twoway_grid.tsv", sep="\t", index=False)

    tests = []

    def add_test(scope, df, col, a_mask, b_mask, group_a, group_b, gene, kind):
        rec = mw(df.loc[a_mask, col], df.loc[b_mask, col])
        rec.update(
            {
                "scope": scope,
                "metric": col,
                "gene": gene,
                "kind": kind,
                "group_a": group_a,
                "group_b": group_b,
                "patients_a": ",".join(df.loc[a_mask, "Patient"].astype(str)),
                "patients_b": ",".join(df.loc[b_mask, "Patient"].astype(str)),
            }
        )
        tests.append(rec)

    for df, tag in [(mal_ok, "cnv_or_program"), (malb_ok, "epi_minus_normal")]:
        mpr_df = df[df["mpr"].isin(MPR_LEVELS)]
        recist_df = df[df["recist"].isin(RECIST_LEVELS)]
        for gene, mean_col, pct_col in [
            ("TACSTD2", f"{'mal' if tag == 'cnv_or_program' else 'malb'}_tacstd2_mean_log1p_cp10k", f"{'mal' if tag == 'cnv_or_program' else 'malb'}_tacstd2_pct_pos"),
            ("CLDN4", f"{'mal' if tag == 'cnv_or_program' else 'malb'}_cldn4_mean_log1p_cp10k", f"{'mal' if tag == 'cnv_or_program' else 'malb'}_cldn4_pct_pos"),
        ]:
            for col, kind in [(mean_col, "mean_log1p_cp10k"), (pct_col, "pct_pos")]:
                add_test(
                    f"{tag}|post_NMPR_vs_MPR",
                    mpr_df,
                    col,
                    (mpr_df["timepoint"] == "post") & (mpr_df["mpr"] == "NMPR"),
                    (mpr_df["timepoint"] == "post") & (mpr_df["mpr"] == "MPR"),
                    "post_NMPR",
                    "post_MPR",
                    gene,
                    kind,
                )
                add_test(
                    f"{tag}|pre_NMPR_vs_MPR",
                    mpr_df,
                    col,
                    (mpr_df["timepoint"] == "pre") & (mpr_df["mpr"] == "NMPR"),
                    (mpr_df["timepoint"] == "pre") & (mpr_df["mpr"] == "MPR"),
                    "pre_NMPR",
                    "pre_MPR",
                    gene,
                    kind,
                )
                add_test(
                    f"{tag}|NMPR_post_vs_pre",
                    mpr_df,
                    col,
                    (mpr_df["mpr"] == "NMPR") & (mpr_df["timepoint"] == "post"),
                    (mpr_df["mpr"] == "NMPR") & (mpr_df["timepoint"] == "pre"),
                    "NMPR_post",
                    "NMPR_pre",
                    gene,
                    kind,
                )
                add_test(
                    f"{tag}|MPR_post_vs_pre",
                    mpr_df,
                    col,
                    (mpr_df["mpr"] == "MPR") & (mpr_df["timepoint"] == "post"),
                    (mpr_df["mpr"] == "MPR") & (mpr_df["timepoint"] == "pre"),
                    "MPR_post",
                    "MPR_pre",
                    gene,
                    kind,
                )
                add_test(
                    f"{tag}|post_SD_vs_PR",
                    recist_df,
                    col,
                    (recist_df["timepoint"] == "post") & (recist_df["recist"] == "SD"),
                    (recist_df["timepoint"] == "post") & (recist_df["recist"] == "PR"),
                    "post_SD",
                    "post_PR",
                    gene,
                    kind,
                )
                add_test(
                    f"{tag}|pre_SD_vs_PR",
                    recist_df,
                    col,
                    (recist_df["timepoint"] == "pre") & (recist_df["recist"] == "SD"),
                    (recist_df["timepoint"] == "pre") & (recist_df["recist"] == "PR"),
                    "pre_SD",
                    "pre_PR",
                    gene,
                    kind,
                )

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(args.outdir / "twoway_tests.tsv", sep="\t", index=False)

    # Primary figure: TACSTD2 + CLDN4 mean log1p, MPR 2×2
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=False)
    rng = np.random.default_rng(7)
    colors = {("pre", "MPR"): "#4C78A8", ("pre", "NMPR"): "#F58518", ("post", "MPR"): "#54A24B", ("post", "NMPR"): "#E45756"}
    for ax, gene, col in [
        (axes[0], "TACSTD2", "mal_tacstd2_mean_log1p_cp10k"),
        (axes[1], "CLDN4", "mal_cldn4_mean_log1p_cp10k"),
    ]:
        xs, labels = [], []
        i = 0
        for t in TIME_LEVELS:
            for r in MPR_LEVELS:
                sub = mal_ok[(mal_ok["timepoint"] == t) & (mal_ok["mpr"] == r)]
                vals = _finite(sub[col])
                x = np.full(len(vals), i, dtype=float)
                if len(vals):
                    x = x + rng.uniform(-0.08, 0.08, size=len(vals))
                    ax.scatter(x, vals, s=42, color=colors[(t, r)], zorder=3, edgecolors="k", linewidths=0.4)
                    ax.hlines(np.median(vals), i - 0.18, i + 0.18, colors="k", linewidth=1.2, zorder=4)
                else:
                    ax.text(i, 0.05, "empty", ha="center", va="bottom", fontsize=8, color="#666666")
                labels.append(f"{t}\n{r}\nn={len(vals)}")
                xs.append(i)
                i += 1
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"Malignant {gene}  (GSE207422)")
        ax.set_ylabel("patient mean log1p(CP10k)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_twoway_mpr.png", dpi=160)
    fig.savefig(args.outdir / "fig_twoway_mpr.pdf")
    plt.close(fig)

    missing = {
        "pre_MPR": "No public GSE207422 pre-treatment biopsy is labeled MPR. P01 is pre / pathologic NE (RECIST SD). P05 and P08 are pre / NMPR.",
        "paired_pre_post": "The 3 pre patients (P01, P05, P08) are not the same people as the 12 post patients. This is a cross-sectional 2×2, not a paired delta.",
        "second_series": "No second public lung scRNA series with both a timepoint axis and a response axis and a malignant/epithelial matrix was found to stack.",
    }

    primary = tests_df[tests_df["scope"] == "cnv_or_program|post_NMPR_vs_MPR"]
    summary = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "n_cells": int(n),
        "n_patients": int(pat["Patient"].nunique()),
        "n_malignant": int(per_cell["is_malignant"].sum()),
        "n_malig_broad": int(per_cell["is_malig_broad"].sum()),
        "cnv_cut": cnv_cut,
        "min_malignant_cells": MIN_MAL,
        "n_patients_malignant_ge10": int(len(mal_ok)),
        "genes_present": {"TACSTD2": "TACSTD2" in expr, "CLDN4": "CLDN4" in expr},
        "mpr_2x2_patients": {
            "pre_MPR": pat[(pat["timepoint"] == "pre") & (pat["mpr"] == "MPR")]["Patient"].tolist(),
            "pre_NMPR": pat[(pat["timepoint"] == "pre") & (pat["mpr"] == "NMPR")]["Patient"].tolist(),
            "post_MPR": pat[(pat["timepoint"] == "post") & (pat["mpr"] == "MPR")]["Patient"].tolist(),
            "post_NMPR": pat[(pat["timepoint"] == "post") & (pat["mpr"] == "NMPR")]["Patient"].tolist(),
            "pre_NE": pat[(pat["timepoint"] == "pre") & (pat["mpr"] == "NE")]["Patient"].tolist(),
        },
        "recist_2x2_patients": {
            "pre_PR": pat[(pat["timepoint"] == "pre") & (pat["recist"] == "PR")]["Patient"].tolist(),
            "pre_SD": pat[(pat["timepoint"] == "pre") & (pat["recist"] == "SD")]["Patient"].tolist(),
            "post_PR": pat[(pat["timepoint"] == "post") & (pat["recist"] == "PR")]["Patient"].tolist(),
            "post_SD": pat[(pat["timepoint"] == "post") & (pat["recist"] == "SD")]["Patient"].tolist(),
        },
        "missing": missing,
        "stacked_series": [],
        "primary_post_NMPR_vs_MPR": primary.to_dict(orient="records"),
        "lineage_counts": {k: int(v) for k, v in per_cell["lineage"].value_counts().items()},
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ["n_cells", "n_malignant", "n_patients_malignant_ge10", "mpr_2x2_patients", "missing"]}, indent=2))


if __name__ == "__main__":
    main()
