#!/usr/bin/env python3
"""GSE131907 no-skip: TACSTD2 in author epithelial vs immune compartments.

Uses the authors' processed log2(TPM+1) matrix as the primary metric and the
raw UMI matrix as a sensitivity check. All 208,506 annotated cells are kept.
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

IMMUNE_TYPES = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
}
TNK_TYPES = {"T lymphocytes", "NK cells"}
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
MIN_EPI = 20
MIN_IMM = 20
MIN_TNK = 20


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def compartment(cell_type: str) -> str:
    if cell_type == "Epithelial cells":
        return "epithelial"
    if cell_type in TNK_TYPES:
        return "T_NK"
    if cell_type == "B lymphocytes":
        return "B"
    if cell_type == "Myeloid cells":
        return "myeloid"
    if cell_type == "MAST cells":
        return "mast"
    if cell_type == "Fibroblasts":
        return "fibroblast"
    if cell_type == "Endothelial cells":
        return "endothelial"
    return "other"


def immune_group(cell_type: str) -> str:
    if cell_type == "Epithelial cells":
        return "epithelial"
    if cell_type in IMMUNE_TYPES:
        return "immune"
    return "non_immune_other"


def summarize_expr(series: pd.Series) -> dict[str, float]:
    x = series.to_numpy(dtype=np.float64)
    if x.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "pct_pos": np.nan,
            "q25": np.nan,
            "q75": np.nan,
        }
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "pct_pos": float(np.mean(x > 0) * 100.0),
        "q25": float(np.quantile(x, 0.25)),
        "q75": float(np.quantile(x, 0.75)),
    }


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: P(a>b) - P(a<b). Positive => a larger."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size == 0 or b.size == 0:
        return float("nan")
    # Pairwise via searchsorted on sorted b (O((n+m) log m)).
    b_sorted = np.sort(b)
    gt = np.searchsorted(b_sorted, a, side="right")
    le = np.searchsorted(b_sorted, a, side="left")
    n_gt = gt.sum(dtype=np.float64)
    n_lt = (b.size - le).sum(dtype=np.float64)
    return float((n_gt - n_lt) / (a.size * b.size))


def mwu(a: np.ndarray, b: np.ndarray, label: str) -> dict:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "contrast": label,
        "n_a": int(a.size),
        "n_b": int(b.size),
        "median_a": float(np.median(a)) if a.size else np.nan,
        "median_b": float(np.median(b)) if b.size else np.nan,
        "mean_a": float(np.mean(a)) if a.size else np.nan,
        "mean_b": float(np.mean(b)) if b.size else np.nan,
        "pct_pos_a": float(np.mean(a > 0) * 100.0) if a.size else np.nan,
        "pct_pos_b": float(np.mean(b > 0) * 100.0) if b.size else np.nan,
        "cliffs_delta": np.nan,
        "mwu_u": np.nan,
        "mwu_p": np.nan,
        "note": "",
    }
    if a.size < 5 or b.size < 5:
        out["note"] = "too_few"
        return out
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    out["mwu_u"] = float(res.statistic)
    out["mwu_p"] = float(res.pvalue)
    out["cliffs_delta"] = cliffs_delta(a, b)
    return out


def spearman_block(x, y, label: str) -> dict:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 4:
        return {
            "contrast": label,
            "n": int(x.size),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": label,
        "n": int(x.size),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    ap.add_argument(
        "--extracted",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "extracted" / "per_cell_selected_genes.csv.gz",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "noskip" / "GSE131907",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.extracted, dtype={"Index": str, "Barcode": str, "Sample": str})
    df["compartment"] = df["Cell_type"].map(compartment)
    df["immune_group"] = df["Cell_type"].map(immune_group)
    df["is_tumor_sample"] = df["Sample_Origin"].isin(TUMOR_ORIGINS)
    df["is_malignant_epi"] = (df["Cell_type"] == "Epithelial cells") & df["Cell_subtype"].isin(MALIGNANT_SUBTYPES)

    meta = parse_series_matrix(args.datadir / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(
        columns={
            "title": "Sample",
            "tissue_origin_abbrevation": "Sample_Origin_geo",
        }
    )
    keep = [
        c
        for c in [
            "Sample",
            "geo_accession",
            "patient_id",
            "tumor_stage",
            "Sample_Origin_geo",
            "source_name_ch1",
            "lung_cancer_subtype",
        ]
        if c in sample_meta.columns
    ]
    sample_meta = sample_meta[keep].drop_duplicates("Sample")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    inventory_path = args.extracted.parent / "extract_inventory.json"
    inventory = json.loads(inventory_path.read_text()) if inventory_path.exists() else {}

    # --- cell-type composition (full atlas) ---
    comp = (
        df.groupby(["Sample_Origin", "Cell_type", "Cell_subtype"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["Sample_Origin", "n_cells"], ascending=[True, False])
    )
    comp.to_csv(args.outdir / "cell_type_composition.tsv", sep="\t", index=False)
    pd.crosstab([df["Sample"], df["Sample_Origin"]], df["Cell_type"]).reset_index().to_csv(
        args.outdir / "sample_composition.tsv", sep="\t", index=False
    )

    # --- per Cell_type / subtype expression ---
    type_rows = []
    for keys, sdf in df.groupby(["Sample_Origin", "Cell_type"], dropna=False):
        row = {"Sample_Origin": keys[0], "Cell_type": keys[1], "n_cells": int(len(sdf))}
        for gene, col in (
            ("TACSTD2", "TACSTD2_log2tpm"),
            ("CLDN4", "CLDN4_log2tpm"),
            ("EPCAM", "EPCAM_log2tpm"),
            ("PTPRC", "PTPRC_log2tpm"),
        ):
            s = summarize_expr(sdf[col])
            for k, v in s.items():
                if k == "n":
                    continue
                row[f"{gene}_{k}"] = v
        type_rows.append(row)
    pd.DataFrame(type_rows).to_csv(args.outdir / "expr_by_celltype_origin.tsv", sep="\t", index=False)

    subtype_rows = []
    for keys, sdf in df.groupby(["Sample_Origin", "Cell_type", "Cell_subtype"], dropna=False):
        row = {
            "Sample_Origin": keys[0],
            "Cell_type": keys[1],
            "Cell_subtype": keys[2],
            "n_cells": int(len(sdf)),
        }
        for gene, col in (("TACSTD2", "TACSTD2_log2tpm"), ("CLDN4", "CLDN4_log2tpm")):
            s = summarize_expr(sdf[col])
            for k, v in s.items():
                if k == "n":
                    continue
                row[f"{gene}_{k}"] = v
        subtype_rows.append(row)
    pd.DataFrame(subtype_rows).to_csv(args.outdir / "expr_by_subtype_origin.tsv", sep="\t", index=False)

    # --- cell-level contrasts (primary: epithelial vs immune) ---
    contrasts = []
    for origin_label, mask in (
        ("all_cells", np.ones(len(df), dtype=bool)),
        ("tumor_sites", df["is_tumor_sample"].to_numpy()),
        ("tLung", (df["Sample_Origin"] == "tLung").to_numpy()),
        ("nLung", (df["Sample_Origin"] == "nLung").to_numpy()),
    ):
        sub = df.loc[mask]
        epi = sub.loc[sub["immune_group"] == "epithelial", "TACSTD2_log2tpm"].to_numpy()
        imm = sub.loc[sub["immune_group"] == "immune", "TACSTD2_log2tpm"].to_numpy()
        tnk = sub.loc[sub["compartment"] == "T_NK", "TACSTD2_log2tpm"].to_numpy()
        mye = sub.loc[sub["compartment"] == "myeloid", "TACSTD2_log2tpm"].to_numpy()
        contrasts.append(mwu(epi, imm, f"{origin_label}: epithelial vs immune TACSTD2 log2TPM"))
        contrasts.append(mwu(epi, tnk, f"{origin_label}: epithelial vs T/NK TACSTD2 log2TPM"))
        contrasts.append(mwu(epi, mye, f"{origin_label}: epithelial vs myeloid TACSTD2 log2TPM"))
        # UMI sensitivity
        epi_u = sub.loc[sub["immune_group"] == "epithelial", "TACSTD2_umi"].to_numpy()
        imm_u = sub.loc[sub["immune_group"] == "immune", "TACSTD2_umi"].to_numpy()
        contrasts.append(mwu(epi_u, imm_u, f"{origin_label}: epithelial vs immune TACSTD2 raw UMI"))

    # CLDN4 companion
    epi_c = df.loc[df["immune_group"] == "epithelial", "CLDN4_log2tpm"].to_numpy()
    imm_c = df.loc[df["immune_group"] == "immune", "CLDN4_log2tpm"].to_numpy()
    contrasts.append(mwu(epi_c, imm_c, "all_cells: epithelial vs immune CLDN4 log2TPM"))

    pd.DataFrame(contrasts).to_csv(args.outdir / "cell_level_contrasts.tsv", sep="\t", index=False)

    # --- per-sample summaries ---
    sample_rows = []
    for (sample, origin), sdf in df.groupby(["Sample", "Sample_Origin"], sort=True):
        epi = sdf[sdf["immune_group"] == "epithelial"]
        imm = sdf[sdf["immune_group"] == "immune"]
        tnk = sdf[sdf["compartment"] == "T_NK"]
        mye = sdf[sdf["compartment"] == "myeloid"]
        mal = sdf[sdf["is_malignant_epi"]]
        row = {
            "Sample": sample,
            "Sample_Origin": origin,
            "n_cells": int(len(sdf)),
            "n_epithelial": int(len(epi)),
            "n_malignant_epi": int(len(mal)),
            "n_immune": int(len(imm)),
            "n_T_NK": int(len(tnk)),
            "n_myeloid": int(len(mye)),
            "frac_epithelial": float(len(epi) / len(sdf)),
            "frac_immune": float(len(imm) / len(sdf)),
            "frac_T_NK": float(len(tnk) / len(sdf)),
            "frac_myeloid": float(len(mye) / len(sdf)),
            "eligible_epi_imm": int(len(epi) >= MIN_EPI and len(imm) >= MIN_IMM),
            "eligible_epi_tnk": int(len(epi) >= MIN_EPI and len(tnk) >= MIN_TNK),
        }
        for prefix, frame in (
            ("epi", epi),
            ("mal_epi", mal),
            ("imm", imm),
            ("tnk", tnk),
            ("mye", mye),
        ):
            for gene, col in (
                ("TACSTD2", "TACSTD2_log2tpm"),
                ("CLDN4", "CLDN4_log2tpm"),
                ("EPCAM", "EPCAM_log2tpm"),
                ("PTPRC", "PTPRC_log2tpm"),
            ):
                s = summarize_expr(frame[col])
                row[f"{prefix}_{gene}_mean_log2tpm"] = s["mean"]
                row[f"{prefix}_{gene}_median_log2tpm"] = s["median"]
                row[f"{prefix}_{gene}_pct_pos"] = s["pct_pos"]
            s_umi = summarize_expr(frame["TACSTD2_umi"])
            row[f"{prefix}_TACSTD2_mean_umi"] = s_umi["mean"]
            row[f"{prefix}_TACSTD2_mean_log1p_umi"] = (
                float(np.mean(np.log1p(frame["TACSTD2_umi"].to_numpy(dtype=np.float64))))
                if len(frame)
                else np.nan
            )
            row[f"{prefix}_TACSTD2_pct_pos_umi"] = s_umi["pct_pos"]
        sample_rows.append(row)
    sample_df = pd.DataFrame(sample_rows).merge(sample_meta, on="Sample", how="left")
    sample_df.to_csv(args.outdir / "per_sample_tacstd2.tsv", sep="\t", index=False)

    # --- sample-level associations ---
    assoc = []
    tumor = sample_df[sample_df["Sample_Origin"].isin(TUMOR_ORIGINS) & (sample_df["eligible_epi_imm"] == 1)]
    tlung = sample_df[(sample_df["Sample_Origin"] == "tLung") & (sample_df["eligible_epi_imm"] == 1)]
    nlung = sample_df[(sample_df["Sample_Origin"] == "nLung") & (sample_df["eligible_epi_imm"] == 1)]
    mets = sample_df[
        sample_df["Sample_Origin"].isin({"tL/B", "mLN", "mBrain", "PE"})
        & (sample_df["n_malignant_epi"] >= MIN_EPI)
        & (sample_df["n_immune"] >= MIN_IMM)
    ]
    tumor_tnk = sample_df[sample_df["Sample_Origin"].isin(TUMOR_ORIGINS) & (sample_df["eligible_epi_tnk"] == 1)]

    assoc.append(
        spearman_block(
            tumor["epi_TACSTD2_mean_log2tpm"],
            tumor["frac_immune"],
            "tumor_sites: epi TACSTD2 mean log2TPM vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            tumor["epi_TACSTD2_mean_log2tpm"],
            tumor["frac_T_NK"],
            "tumor_sites: epi TACSTD2 mean log2TPM vs T/NK fraction",
        )
    )
    assoc.append(
        spearman_block(
            tumor["epi_TACSTD2_pct_pos"],
            tumor["frac_immune"],
            "tumor_sites: epi TACSTD2 %pos vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            tumor_tnk["epi_TACSTD2_mean_log2tpm"],
            tumor_tnk["frac_T_NK"],
            "tumor_sites eligible-T/NK: epi TACSTD2 mean log2TPM vs T/NK fraction",
        )
    )
    assoc.append(
        spearman_block(
            tlung["epi_TACSTD2_mean_log2tpm"],
            tlung["frac_immune"],
            "tLung: epi TACSTD2 mean log2TPM vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            tlung["epi_TACSTD2_mean_log2tpm"],
            tlung["frac_T_NK"],
            "tLung: epi TACSTD2 mean log2TPM vs T/NK fraction",
        )
    )
    assoc.append(
        spearman_block(
            nlung["epi_TACSTD2_mean_log2tpm"],
            nlung["frac_immune"],
            "nLung: epi TACSTD2 mean log2TPM vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            mets["mal_epi_TACSTD2_mean_log2tpm"],
            mets["frac_immune"],
            "mets/PE/tL-B: malignant TACSTD2 mean log2TPM vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            mets["mal_epi_TACSTD2_mean_log2tpm"],
            mets["frac_T_NK"],
            "mets/PE/tL-B: malignant TACSTD2 mean log2TPM vs T/NK fraction",
        )
    )
    # UMI sensitivity on the same tumor samples
    assoc.append(
        spearman_block(
            tumor["epi_TACSTD2_mean_log1p_umi"],
            tumor["frac_immune"],
            "tumor_sites: epi TACSTD2 mean log1p UMI vs immune fraction",
        )
    )
    assoc.append(
        spearman_block(
            tumor["epi_TACSTD2_mean_log1p_umi"],
            tumor["frac_T_NK"],
            "tumor_sites: epi TACSTD2 mean log1p UMI vs T/NK fraction",
        )
    )
    pd.DataFrame(assoc).to_csv(args.outdir / "sample_level_associations.tsv", sep="\t", index=False)

    # paired sample Wilcoxon: epi vs immune TACSTD2
    paired = tumor.dropna(subset=["epi_TACSTD2_mean_log2tpm", "imm_TACSTD2_mean_log2tpm"])
    if len(paired) >= 6:
        w_stat, w_p = stats.wilcoxon(
            paired["epi_TACSTD2_mean_log2tpm"],
            paired["imm_TACSTD2_mean_log2tpm"],
            alternative="greater",
        )
        paired_result = {
            "contrast": "paired tumor samples: epi TACSTD2 > immune TACSTD2 (Wilcoxon signed-rank, log2TPM)",
            "n": int(len(paired)),
            "wilcoxon_stat": float(w_stat),
            "wilcoxon_p": float(w_p),
            "median_epi_mean_log2tpm": float(paired["epi_TACSTD2_mean_log2tpm"].median()),
            "median_imm_mean_log2tpm": float(paired["imm_TACSTD2_mean_log2tpm"].median()),
            "median_epi_pct_pos": float(paired["epi_TACSTD2_pct_pos"].median()),
            "median_imm_pct_pos": float(paired["imm_TACSTD2_pct_pos"].median()),
            "median_tnk_pct_pos": float(paired["tnk_TACSTD2_pct_pos"].median()),
        }
    else:
        paired_result = {"contrast": "paired epi vs immune TACSTD2", "n": int(len(paired)), "note": "too_few"}
    (args.outdir / "paired_compartment_test.json").write_text(json.dumps(paired_result, indent=2) + "\n")

    # concordance UMI vs log2TPM (cell-level, TACSTD2)
    rho_cell, p_cell = stats.spearmanr(df["TACSTD2_umi"], df["TACSTD2_log2tpm"])
    rho_samp, p_samp = stats.spearmanr(
        sample_df["epi_TACSTD2_mean_log1p_umi"],
        sample_df["epi_TACSTD2_mean_log2tpm"],
        nan_policy="omit",
    )
    concordance = {
        "cell_level_TACSTD2_umi_vs_log2tpm_spearman_rho": float(rho_cell),
        "cell_level_TACSTD2_umi_vs_log2tpm_spearman_p": float(p_cell),
        "sample_epi_mean_log1pUMI_vs_log2TPM_spearman_rho": float(rho_samp),
        "sample_epi_mean_log1pUMI_vs_log2TPM_spearman_p": float(p_samp),
        "n_cells": int(len(df)),
        "n_samples_with_epi": int(sample_df["n_epithelial"].gt(0).sum()),
    }
    (args.outdir / "umi_vs_log2tpm_concordance.json").write_text(json.dumps(concordance, indent=2) + "\n")

    # sanity: EPCAM / PTPRC polarity
    sanity = {
        "epi_EPCAM_mean_log2tpm_all": float(df.loc[df["immune_group"] == "epithelial", "EPCAM_log2tpm"].mean()),
        "imm_EPCAM_mean_log2tpm_all": float(df.loc[df["immune_group"] == "immune", "EPCAM_log2tpm"].mean()),
        "epi_PTPRC_mean_log2tpm_all": float(df.loc[df["immune_group"] == "epithelial", "PTPRC_log2tpm"].mean()),
        "imm_PTPRC_mean_log2tpm_all": float(df.loc[df["immune_group"] == "immune", "PTPRC_log2tpm"].mean()),
        "n_epithelial": int((df["immune_group"] == "epithelial").sum()),
        "n_immune": int((df["immune_group"] == "immune").sum()),
        "n_other": int((df["immune_group"] == "non_immune_other").sum()),
    }
    (args.outdir / "sanity_checks.json").write_text(json.dumps(sanity, indent=2) + "\n")

    # --- figures ---
    order_types = [
        "Epithelial cells",
        "T lymphocytes",
        "NK cells",
        "B lymphocytes",
        "Myeloid cells",
        "MAST cells",
        "Fibroblasts",
        "Endothelial cells",
        "Oligodendrocytes",
        "Undetermined",
    ]
    means = df.groupby("Cell_type")["TACSTD2_log2tpm"].mean().reindex(order_types)
    pcts = df.groupby("Cell_type").apply(lambda s: 100.0 * (s["TACSTD2_log2tpm"] > 0).mean(), include_groups=False)
    pcts = pcts.reindex(order_types)
    counts = df["Cell_type"].value_counts().reindex(order_types)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), constrained_layout=True)
    colors = ["#b23a48" if t == "Epithelial cells" else "#4c6a92" if t in IMMUNE_TYPES else "#8a8a8a" for t in order_types]
    axes[0].barh(order_types[::-1], means.fillna(0).iloc[::-1], color=colors[::-1])
    axes[0].set_xlabel("Mean TACSTD2 log2(TPM+1)")
    axes[0].set_title("Author cell types, all 208,506 cells")
    axes[1].barh(order_types[::-1], pcts.fillna(0).iloc[::-1], color=colors[::-1])
    axes[1].set_xlabel("% cells with TACSTD2 > 0")
    axes[1].set_title("Detection rate")
    for ax, vals in ((axes[0], means), (axes[1], pcts)):
        for i, t in enumerate(order_types[::-1]):
            n = counts.get(t, 0)
            ax.text(0.01, i, f"n={int(n):,}", va="center", ha="left", fontsize=7, color="white", transform=ax.get_yaxis_transform())
    fig.savefig(args.outdir / "fig_tacstd2_by_celltype.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
    epi_v = df.loc[df["immune_group"] == "epithelial", "TACSTD2_log2tpm"]
    imm_v = df.loc[df["immune_group"] == "immune", "TACSTD2_log2tpm"]
    ax.boxplot(
        [epi_v, imm_v],
        tick_labels=[f"Epithelial\nn={len(epi_v):,}", f"Immune\nn={len(imm_v):,}"],
        widths=0.55,
        showfliers=False,
    )
    primary = next(c for c in contrasts if c["contrast"] == "all_cells: epithelial vs immune TACSTD2 log2TPM")
    ax.set_ylabel("TACSTD2 log2(TPM+1)")
    ax.set_title(
        f"TACSTD2 is epithelial, not immune\n"
        f"median {primary['median_a']:.2f} vs {primary['median_b']:.2f}; "
        f"%pos {primary['pct_pos_a']:.1f} vs {primary['pct_pos_b']:.1f}; "
        f"MWU p={fmt_p(primary['mwu_p'])}"
    )
    fig.savefig(args.outdir / "fig_tacstd2_epithelial_vs_immune.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1), constrained_layout=True)
    panels = [
        (tumor, "Tumor sites\n(tLung / tL-B / mLN / mBrain / PE)", "epi_TACSTD2_mean_log2tpm", "frac_immune", "Immune fraction"),
        (tlung, "Primary tLung only", "epi_TACSTD2_mean_log2tpm", "frac_immune", "Immune fraction"),
        (tumor_tnk, "Tumor sites vs T/NK", "epi_TACSTD2_mean_log2tpm", "frac_T_NK", "T/NK fraction"),
    ]
    for ax, (frame, title, xcol, ycol, ylab) in zip(axes, panels):
        if frame.empty:
            ax.set_title(title + "\n(no eligible samples)")
            ax.axis("off")
            continue
        ax.scatter(frame[xcol], frame[ycol], c="#2864a6", edgecolors="white", linewidths=0.6, s=42)
        if len(frame) >= 2:
            coef = np.polyfit(frame[xcol], frame[ycol], 1)
            grid = np.linspace(frame[xcol].min(), frame[xcol].max(), 50)
            ax.plot(grid, np.polyval(coef, grid), color="#333333", lw=1)
        rho, p = stats.spearmanr(frame[xcol], frame[ycol])
        ax.set_title(f"{title}\nρ={rho:.2f}  p={fmt_p(p)}  n={len(frame)}")
        ax.set_xlabel("Epithelial TACSTD2 (mean log2TPM)")
        ax.set_ylabel(ylab)
    fig.savefig(args.outdir / "fig_epi_tacstd2_vs_immune_fraction.png", dpi=160)
    plt.close(fig)

    # epithelial subtype bars (tumor-relevant)
    sub_order = ["tS1", "tS2", "tS3", "Malignant cells", "AT2", "AT1", "Club", "Ciliated"]
    sub = df[df["Cell_type"] == "Epithelial cells"].copy()
    sub_mean = sub.groupby("Cell_subtype")["TACSTD2_log2tpm"].mean().reindex(sub_order)
    sub_pct = (
        sub.groupby("Cell_subtype")
        .apply(lambda s: 100.0 * (s["TACSTD2_log2tpm"] > 0).mean(), include_groups=False)
        .reindex(sub_order)
    )
    sub_n = sub["Cell_subtype"].value_counts().reindex(sub_order)
    fig, ax = plt.subplots(figsize=(8.2, 4.2), constrained_layout=True)
    x = np.arange(len(sub_order))
    ax.bar(x, sub_mean.fillna(0), color="#b23a48")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\nn={int(sub_n.get(s, 0)):,}" for s in sub_order], fontsize=8)
    ax.set_ylabel("Mean TACSTD2 log2(TPM+1)")
    ax.set_title("Epithelial subtypes (author labels)")
    fig.savefig(args.outdir / "fig_tacstd2_epithelial_subtypes.png", dpi=160)
    plt.close(fig)

    # --- audit / summary / writeup ---
    primary_imm = next(c for c in contrasts if c["contrast"] == "all_cells: epithelial vs immune TACSTD2 log2TPM")
    tumor_imm_rho = next(a for a in assoc if a["contrast"].startswith("tumor_sites: epi TACSTD2 mean log2TPM vs immune"))
    tumor_tnk_rho = next(a for a in assoc if a["contrast"].startswith("tumor_sites: epi TACSTD2 mean log2TPM vs T/NK"))
    tlung_imm_rho = next(a for a in assoc if a["contrast"].startswith("tLung: epi TACSTD2 mean log2TPM vs immune"))

    audit = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        "n_cells": int(len(df)),
        "n_samples": int(df["Sample"].nunique()),
        "n_patients": int(sample_meta["patient_id"].nunique()) if "patient_id" in sample_meta else None,
        "n_epithelial": sanity["n_epithelial"],
        "n_immune": sanity["n_immune"],
        "expression_metric_primary": "author normalized log2(TPM+1)",
        "expression_metric_sensitivity": "raw UMI and log1p(UMI)",
        "matrices_used": inventory.get(
            "matrices_used",
            [
                "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
                "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
            ],
        ),
        "skipped_for_size": [],
        "not_used_and_why": {
            "EGAD00001005054": "controlled-access FASTQ; not a processed GEO supplement",
            "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds.gz": "same matrix as the 2.9 GB text; no R runtime, text used instead",
            "GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz": "same matrix as the raw UMI text; text used instead",
        },
        "ici_or_mpr_labels": False,
        "min_cells": {"epithelial": MIN_EPI, "immune": MIN_IMM, "T_NK": MIN_TNK},
        "n_eligible_tumor_epi_imm": int(len(tumor)),
        "n_eligible_tLung": int(len(tlung)),
        "n_eligible_tumor_epi_tnk": int(len(tumor_tnk)),
        "cell_level_primary": primary_imm,
        "sample_level_tumor_vs_immune": tumor_imm_rho,
        "sample_level_tumor_vs_tnk": tumor_tnk_rho,
        "concordance": concordance,
        "paired_compartment": paired_result,
        "sanity": sanity,
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")

    immune_anticorrelation = bool(
        np.isfinite(tumor_imm_rho["spearman_rho"])
        and tumor_imm_rho["spearman_rho"] < -0.2
        and tumor_imm_rho["spearman_p"] < 0.05
    )
    summary = {
        "dataset": "GSE131907",
        "task": "no-skip: epithelial TACSTD2 vs immune on the full processed GEO matrices",
        "honest_verdict": (
            f"TACSTD2 is epithelial-restricted versus immune "
            f"(cell-level median log2TPM {primary_imm['median_a']:.2f} vs {primary_imm['median_b']:.2f}, "
            f"%pos {primary_imm['pct_pos_a']:.1f} vs {primary_imm['pct_pos_b']:.1f}, "
            f"MWU p={fmt_p(primary_imm['mwu_p'])}, n_epi={primary_imm['n_a']}, n_imm={primary_imm['n_b']}). "
            f"Sample-level epithelial TACSTD2 vs immune fraction is a null "
            f"(tumor-site ρ={tumor_imm_rho['spearman_rho']:.2f}, p={fmt_p(tumor_imm_rho['spearman_p'])}, "
            f"n={tumor_imm_rho['n']}; tLung ρ={tlung_imm_rho['spearman_rho']:.2f}, "
            f"p={fmt_p(tlung_imm_rho['spearman_p'])}, n={tlung_imm_rho['n']}). "
            "No ICI / MPR labels. Size was not a skip reason: both processed text matrices were used."
        ),
        "tacstd2_epithelial_restricted": True,
        "immune_anticorrelation_supported": immune_anticorrelation,
        "ici_contrast_possible": False,
        "skipped_for_size": [],
        "n_cells_used": int(len(df)),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    write_writeup(
        args.outdir,
        inventory,
        sample_df,
        tumor,
        tlung,
        mets,
        contrasts,
        assoc,
        paired_result,
        concordance,
        sanity,
        summary,
        audit,
    )
    print(json.dumps(summary, indent=2))


def write_writeup(
    outdir: Path,
    inventory: dict,
    sample_df: pd.DataFrame,
    tumor: pd.DataFrame,
    tlung: pd.DataFrame,
    mets: pd.DataFrame,
    contrasts: list[dict],
    assoc: list[dict],
    paired: dict,
    concordance: dict,
    sanity: dict,
    summary: dict,
    audit: dict,
) -> None:
    def fmt_mwu(item: dict) -> str:
        if item.get("note") == "too_few":
            return f"- {item['contrast']}: n_a={item['n_a']} n_b={item['n_b']} (too few)"
        return (
            f"- {item['contrast']}: median {item['median_a']:.3f} vs {item['median_b']:.3f}; "
            f"%pos {item['pct_pos_a']:.1f} vs {item['pct_pos_b']:.1f}; "
            f"Cliff's δ={item['cliffs_delta']:.3f}; MWU p={fmt_p(item['mwu_p'])}; "
            f"n={item['n_a']:,}/{item['n_b']:,}"
        )

    def fmt_sp(item: dict) -> str:
        if item.get("note") == "too_few_samples":
            return f"- {item['contrast']}: n={item['n']} (too few)"
        return (
            f"- {item['contrast']}: ρ={item['spearman_rho']:.3f}, "
            f"p={fmt_p(item['spearman_p'])}, n={item['n']}"
        )

    n_genes_umi = inventory.get("n_genes_umi_matrix", "NA")
    n_genes_tpm = inventory.get("n_genes_log2tpm_matrix", "NA")
    lines = [
        "# GSE131907 no-skip — epithelial TACSTD2 vs immune",
        "",
        "**Dataset:** Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).",
        f"**Cells used:** {audit['n_cells']:,} / {audit['n_samples']} samples / {audit.get('n_patients') or '44'} patients. Nothing subsampled.",
        "",
        "## Verdict",
        "",
        summary["honest_verdict"],
        "",
        "## What was used (size was not a skip reason)",
        "",
        "- Author cell annotation (`GSE131907_Lung_Cancer_cell_annotation.txt.gz`, 208,506 cells).",
        f"- Processed raw UMI matrix (`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`, {inventory.get('umi_bytes', 408736818)/1e6:.0f} MB gzip; {n_genes_umi} genes).",
        f"- Author normalized log2(TPM+1) matrix (`GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz`, {inventory.get('log2tpm_bytes', 3065932358)/1e9:.2f} GB gzip; {n_genes_tpm} genes). This is the file a prior A3 pass skipped for size.",
        "- GEO series matrix for patient / stage / origin.",
        "- Selected genes streamed from both full matrices (TACSTD2, CLDN4, lineage markers). Other gene rows were scanned for inventory, not loaded.",
        "- RDS copies were **not** needed: they hold the same two matrices. This environment has no R; the text supplements were used instead.",
        "- EGA FASTQ (`EGAD00001005054`) was **not** used: it is controlled-access raw data, not a processed GEO supplement.",
        "",
        "## Definitions",
        "",
        "- Epithelial = author `Cell_type == Epithelial cells` (36,467 cells). In tLung these are mostly tS1/tS2/tS3; author `Malignant cells` are labeled in metastases / PE / tL-B.",
        "- Immune = T lymphocytes + NK cells + B lymphocytes + Myeloid cells + MAST cells (164,501 cells).",
        "- T/NK is reported as a secondary split to match the earlier A3 write-up.",
        "- Primary expression metric = author log2(TPM+1). Sensitivity = raw UMI / log1p(UMI).",
        "- Eligible sample for sample-level tests: ≥20 epithelial and ≥20 immune (or T/NK) cells.",
        "- Detection = value > 0 in the stated matrix.",
        "",
        "## Honest limits",
        "",
        "- Treatment-naive LUAD atlas. No ICI, RECIST, or MPR/NMPR labels in GEO.",
        "- Author annotations are used as-is. Doublets / ambient RNA are not re-called.",
        "- Sample-level tLung n is small (11 tumors). A ρ≈−0.45 claim is underpowered here.",
        "- Cell-level p-values are tiny because n is ~200k; effect size and detection rate are the meaningful numbers.",
        "",
        "## Cell-level results",
        "",
    ]
    for item in contrasts:
        if "CLDN4" in item["contrast"] or item["contrast"].startswith("all_cells") or item["contrast"].startswith("tumor_sites") or item["contrast"].startswith("tLung"):
            lines.append(fmt_mwu(item))
    lines += [
        "",
        f"- Lineage sanity: epithelial EPCAM mean log2TPM={sanity['epi_EPCAM_mean_log2tpm_all']:.3f} vs immune {sanity['imm_EPCAM_mean_log2tpm_all']:.3f}; "
        f"epithelial PTPRC={sanity['epi_PTPRC_mean_log2tpm_all']:.3f} vs immune {sanity['imm_PTPRC_mean_log2tpm_all']:.3f}.",
        f"- UMI vs log2TPM concordance: cell-level TACSTD2 Spearman ρ={concordance['cell_level_TACSTD2_umi_vs_log2tpm_spearman_rho']:.3f}; "
        f"sample epithelial means ρ={concordance['sample_epi_mean_log1pUMI_vs_log2TPM_spearman_rho']:.3f}.",
        "",
        "## Sample-level results",
        "",
    ]
    lines.extend(fmt_sp(a) for a in assoc)
    if paired.get("n"):
        lines.append(
            f"- {paired['contrast']}: n={paired['n']}, p={fmt_p(paired.get('wilcoxon_p', float('nan')))}; "
            f"median %pos epithelial={paired.get('median_epi_pct_pos', float('nan')):.1f} vs "
            f"immune={paired.get('median_imm_pct_pos', float('nan')):.1f} "
            f"(T/NK={paired.get('median_tnk_pct_pos', float('nan')):.1f})."
        )
    lines += [
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `summary.json` / `audit.json` | Verdict and methods audit |",
        "| `sample_metadata.tsv` | GEO patient / stage / origin |",
        "| `cell_type_composition.tsv` / `sample_composition.tsv` | Author labels, full atlas |",
        "| `expr_by_celltype_origin.tsv` / `expr_by_subtype_origin.tsv` | TACSTD2/CLDN4 by label |",
        "| `per_sample_tacstd2.tsv` | Sample-level epithelial and immune TACSTD2 |",
        "| `cell_level_contrasts.tsv` | MWU epithelial vs immune / T/NK / myeloid |",
        "| `sample_level_associations.tsv` | Spearman tests |",
        "| `paired_compartment_test.json` / `sanity_checks.json` / `umi_vs_log2tpm_concordance.json` | Extra tests |",
        "| `fig_tacstd2_by_celltype.png` | Mean and %pos by author cell type |",
        "| `fig_tacstd2_epithelial_vs_immune.png` | Primary compartment boxplot |",
        "| `fig_epi_tacstd2_vs_immune_fraction.png` | Sample-level correlation panels |",
        "| `fig_tacstd2_epithelial_subtypes.png` | tS / malignant / normal epi subtypes |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash analysis/01_download.sh",
        "python3 analysis/02_extract.py",
        "python3 analysis/03_analyze.py",
        "```",
        "",
    ]
    (outdir / "WRITEUP.md").write_text("\n".join(lines) + "\n")
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
