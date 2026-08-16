#!/usr/bin/env python3
"""A3 analog on GSE131907: epithelial TACSTD2 vs T/NK from the processed UMI matrix."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SELECTED = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B",
    "NKG7",
    "GZMB",
]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MIN_EPI = 20
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


def stream_selected_genes(matrix_path: Path, wanted: set[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene not in wanted:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            found[gene] = arr
            if len(found) == len(wanted):
                break
    return cell_ids, found


def compartment(row: pd.Series) -> str:
    if row.Cell_type == "Epithelial cells":
        return "epithelial"
    if row.Cell_type in {"T lymphocytes", "NK cells"}:
        return "T_NK"
    if row.Cell_type == "B lymphocytes":
        return "B"
    if row.Cell_type == "Myeloid cells":
        return "myeloid"
    return "other"


def summarize_group(df: pd.DataFrame, genes: list[str]) -> dict:
    out = {"n_cells": int(len(df))}
    if df.empty:
        for gene in genes:
            out[f"{gene}_mean_umi"] = np.nan
            out[f"{gene}_mean_log1p"] = np.nan
            out[f"{gene}_pct_pos"] = np.nan
        return out
    for gene in genes:
        umi = df[gene].to_numpy()
        out[f"{gene}_mean_umi"] = float(np.mean(umi))
        out[f"{gene}_mean_log1p"] = float(np.mean(np.log1p(umi)))
        out[f"{gene}_pct_pos"] = float(np.mean(umi > 0) * 100.0)
    return out


def spearman_block(x: np.ndarray, y: np.ndarray, label: str) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return {
            "contrast": label,
            "n": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": label,
        "n": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "results" / "w200" / "A3_GSE131907",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(args.workdir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    meta = parse_series_matrix(args.workdir / "GSE131907_series_matrix.txt.gz")
    sample_meta = meta.rename(columns={"title": "Sample", "tissue_origin_abbrevation": "Sample_Origin_geo"})
    sample_cols = [
        c
        for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "Sample_Origin_geo", "source_name_ch1"]
        if c in sample_meta.columns
    ]
    sample_meta = sample_meta[sample_cols].drop_duplicates("Sample")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    comp = (
        ann.groupby(["Sample_Origin", "Cell_type", "Cell_subtype"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["Sample_Origin", "n_cells"], ascending=[True, False])
    )
    comp.to_csv(args.outdir / "cell_type_composition.tsv", sep="\t", index=False)

    sample_comp = pd.crosstab([ann["Sample"], ann["Sample_Origin"]], ann["Cell_type"]).reset_index()
    sample_comp.to_csv(args.outdir / "sample_composition.tsv", sep="\t", index=False)

    matrix_path = args.workdir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    print(f"streaming selected genes from {matrix_path}")
    cell_ids, expr = stream_selected_genes(matrix_path, set(SELECTED))
    missing = sorted(set(SELECTED) - set(expr))
    print("genes found", sorted(expr), "missing", missing)

    per_cell = ann.set_index("Index").reindex(cell_ids).reset_index()
    if per_cell["Sample"].isna().any():
        raise SystemExit("matrix cell IDs do not align with annotation Index")
    for gene, arr in expr.items():
        per_cell[gene] = arr
    per_cell["compartment"] = per_cell.apply(compartment, axis=1)
    per_cell["is_tumor_sample"] = per_cell["Sample_Origin"].isin(TUMOR_ORIGINS)
    per_cell["tumor_epithelial"] = (
        (per_cell["Cell_type"] == "Epithelial cells")
        & per_cell["Cell_subtype"].isin(["Malignant cells", "tS1", "tS2", "tS3"])
    )

    compact_cols = [
        "Index",
        "Sample",
        "Sample_Origin",
        "Cell_type",
        "Cell_subtype",
        "compartment",
        "tumor_epithelial",
    ] + [g for g in SELECTED if g in per_cell]
    per_cell[compact_cols].to_csv(args.workdir / "per_cell_selected_genes.csv.gz", index=False)

    genes_present = [g for g in SELECTED if g in per_cell]
    rows = []
    for (sample, origin), sdf in per_cell.groupby(["Sample", "Sample_Origin"], sort=True):
        epi = sdf[sdf["compartment"] == "epithelial"]
        tnk = sdf[sdf["compartment"] == "T_NK"]
        malig = sdf[sdf["tumor_epithelial"]]
        row = {
            "Sample": sample,
            "Sample_Origin": origin,
            "n_cells": int(len(sdf)),
            "n_epithelial": int(len(epi)),
            "n_tumor_epithelial": int(len(malig)),
            "n_T_NK": int(len(tnk)),
            "frac_epithelial": float(len(epi) / len(sdf)),
            "frac_T_NK": float(len(tnk) / len(sdf)),
            "eligible_epi_tnk": int(len(epi) >= MIN_EPI and len(tnk) >= MIN_TNK),
        }
        for prefix, frame in (("epi", epi), ("tumor_epi", malig), ("tnk", tnk)):
            stats_row = summarize_group(frame, genes_present)
            for key, value in stats_row.items():
                if key == "n_cells":
                    continue
                row[f"{prefix}_{key}"] = value
        rows.append(row)
    sample_df = pd.DataFrame(rows).merge(sample_meta, on="Sample", how="left")
    sample_df.to_csv(args.outdir / "per_sample_tacstd2.tsv", sep="\t", index=False)

    contrasts = []
    tumor = sample_df[sample_df["Sample_Origin"].isin(TUMOR_ORIGINS) & (sample_df["eligible_epi_tnk"] == 1)]
    tlung = sample_df[(sample_df["Sample_Origin"] == "tLung") & (sample_df["eligible_epi_tnk"] == 1)]
    mets = sample_df[
        sample_df["Sample_Origin"].isin({"tL/B", "mLN", "mBrain", "PE"})
        & (sample_df["n_tumor_epithelial"] >= MIN_EPI)
        & (sample_df["n_T_NK"] >= MIN_TNK)
    ]
    contrasts.append(
        spearman_block(
            tumor["epi_TACSTD2_mean_log1p"].to_numpy(),
            tumor["frac_T_NK"].to_numpy(),
            "tumor_sites: epi TACSTD2 mean_log1p vs T/NK fraction",
        )
    )
    contrasts.append(
        spearman_block(
            tumor["epi_TACSTD2_pct_pos"].to_numpy(),
            tumor["frac_T_NK"].to_numpy(),
            "tumor_sites: epi TACSTD2 %pos vs T/NK fraction",
        )
    )
    contrasts.append(
        spearman_block(
            tlung["epi_TACSTD2_mean_log1p"].to_numpy(),
            tlung["frac_T_NK"].to_numpy(),
            "tLung: epi TACSTD2 mean_log1p vs T/NK fraction",
        )
    )
    contrasts.append(
        spearman_block(
            mets["tumor_epi_TACSTD2_mean_log1p"].to_numpy(),
            mets["frac_T_NK"].to_numpy(),
            "metastasis/PE/tL-B: malignant TACSTD2 mean_log1p vs T/NK fraction",
        )
    )

    paired = tumor.dropna(subset=["epi_TACSTD2_mean_log1p", "tnk_TACSTD2_mean_log1p"])
    if len(paired) >= 6:
        w_stat, w_p = stats.wilcoxon(
            paired["epi_TACSTD2_mean_log1p"],
            paired["tnk_TACSTD2_mean_log1p"],
            alternative="greater",
        )
        paired_result = {
            "contrast": "paired tumor samples: epi TACSTD2 > T/NK TACSTD2 (Wilcoxon signed-rank)",
            "n": int(len(paired)),
            "wilcoxon_stat": float(w_stat),
            "wilcoxon_p": float(w_p),
            "median_epi_mean_log1p": float(paired["epi_TACSTD2_mean_log1p"].median()),
            "median_tnk_mean_log1p": float(paired["tnk_TACSTD2_mean_log1p"].median()),
            "median_epi_pct_pos": float(paired["epi_TACSTD2_pct_pos"].median()),
            "median_tnk_pct_pos": float(paired["tnk_TACSTD2_pct_pos"].median()),
        }
    else:
        paired_result = {"contrast": "paired epi vs T/NK TACSTD2", "n": int(len(paired)), "note": "too_few"}

    sanity = []
    if not tumor.empty:
        sanity.append(
            {
                "check": "EPCAM higher in epithelial than T/NK (tumor samples, mean of sample means)",
                "epi_EPCAM_mean_log1p": float(tumor["epi_EPCAM_mean_log1p"].mean()),
                "tnk_EPCAM_mean_log1p": float(tumor["tnk_EPCAM_mean_log1p"].mean()),
                "epi_PTPRC_mean_log1p": float(tumor["epi_PTPRC_mean_log1p"].mean()),
                "tnk_PTPRC_mean_log1p": float(tumor["tnk_PTPRC_mean_log1p"].mean()),
            }
        )

    stats_rows = []
    for item in contrasts:
        stats_rows.append(item)
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(args.outdir / "association_statistics.tsv", sep="\t", index=False)
    with (args.outdir / "paired_compartment_test.json").open("w") as fh:
        json.dump(paired_result, fh, indent=2)
    with (args.outdir / "sanity_checks.json").open("w") as fh:
        json.dump(sanity, fh, indent=2)

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), constrained_layout=True)
    panels = [
        (tumor, "Tumor sites (tLung/tL-B/mLN/mBrain/PE)", "epi_TACSTD2_mean_log1p"),
        (tlung, "Primary tLung only", "epi_TACSTD2_mean_log1p"),
        (mets, "Author malignant cells (mets/PE/tL-B)", "tumor_epi_TACSTD2_mean_log1p"),
    ]
    for ax, (frame, title, xcol) in zip(axes, panels):
        if frame.empty:
            ax.set_title(title + "\n(no eligible samples)")
            ax.axis("off")
            continue
        ax.scatter(frame[xcol], frame["frac_T_NK"], c="#2864a6", edgecolors="white", linewidths=0.6, s=42)
        if len(frame) >= 2:
            coef = np.polyfit(frame[xcol], frame["frac_T_NK"], 1)
            grid = np.linspace(frame[xcol].min(), frame[xcol].max(), 50)
            ax.plot(grid, np.polyval(coef, grid), color="#333333", lw=1)
        rho, p = stats.spearmanr(frame[xcol], frame["frac_T_NK"])
        ax.set_title(f"{title}\nρ={rho:.2f}  p={p:.3g}  n={len(frame)}")
        ax.set_xlabel("Epithelial TACSTD2 (mean log1p UMI)")
        ax.set_ylabel("T/NK fraction")
    fig.savefig(args.outdir / "fig_epi_tacstd2_vs_tnk_fraction.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    if not paired.empty:
        x = np.arange(2)
        meds = [paired["epi_TACSTD2_pct_pos"].median(), paired["tnk_TACSTD2_pct_pos"].median()]
        ax.boxplot(
            [paired["epi_TACSTD2_pct_pos"], paired["tnk_TACSTD2_pct_pos"]],
            tick_labels=["Epithelial", "T/NK"],
            widths=0.55,
        )
        ax.set_ylabel("% TACSTD2+ cells (UMI>0)")
        ax.set_title(
            f"Tumor samples: TACSTD2 is epithelial-restricted\n"
            f"median {meds[0]:.1f}% vs {meds[1]:.1f}%; Wilcoxon p={paired_result.get('wilcoxon_p', float('nan')):.2g}"
        )
    fig.savefig(args.outdir / "fig_tacstd2_epithelial_vs_tnk.png", dpi=160)
    plt.close(fig)

    audit = {
        "dataset": "GSE131907",
        "n_cells": int(len(per_cell)),
        "n_samples": int(per_cell["Sample"].nunique()),
        "genes_found": sorted(expr),
        "genes_missing": missing,
        "n_eligible_tumor_samples": int(len(tumor)),
        "n_eligible_tLung": int(len(tlung)),
        "n_eligible_malignant_mets": int(len(mets)),
        "mpr_or_ici_labels": False,
        "expression_metric": "mean log1p(raw UMI) and % cells with UMI>0; CP10k not computed (full-matrix totals skipped)",
        "min_cells": {"epithelial": MIN_EPI, "T_NK": MIN_TNK},
        "associations": contrasts,
        "paired_compartment": paired_result,
    }
    with (args.outdir / "audit.json").open("w") as fh:
        json.dump(audit, fh, indent=2)

    write_results_readme(args.outdir, sample_df, tumor, tlung, mets, contrasts, paired_result, missing)
    print(json.dumps(audit, indent=2))


def write_results_readme(
    outdir: Path,
    sample_df: pd.DataFrame,
    tumor: pd.DataFrame,
    tlung: pd.DataFrame,
    mets: pd.DataFrame,
    contrasts: list[dict],
    paired: dict,
    missing: list[str],
) -> None:
    def fmt(item: dict) -> str:
        if item.get("note") == "too_few_samples":
            return f"- {item['contrast']}: n={item['n']} (too few)"
        return (
            f"- {item['contrast']}: ρ={item['spearman_rho']:.3f}, "
            f"p={item['spearman_p']:.3g}, n={item['n']}"
        )

    n_tlung_all = int((sample_df.Sample_Origin == "tLung").sum())
    lines = [
        "# W200-A3 · GSE131907 — epithelial TACSTD2 vs T/NK",
        "",
        "**Task:** A3 analog (malignant/epithelial TACSTD2 vs immune) on GSE131907, using processed data only.",
        "",
        "## Verdict",
        "",
        "Processed data **were** in budget. The A3-style **epithelial TACSTD2 vs T/NK** correlation was run.",
        "The A3 **NMPR > MPR** half **cannot** be tested: this atlas has no immunotherapy or pathologic-response labels.",
        "Raw EGA FASTQ and the 2.86 GB log2TPM text matrix were **not** downloaded.",
        "",
        "## What was used",
        "",
        "- Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).",
        "- 208,506 cells / 58 samples / 44 LUAD patients; author `Cell_type` and `Cell_subtype`.",
        "- Expression: processed raw UMI matrix (0.38 GB gzip). Selected genes streamed; other rows not parsed.",
        "- Epithelial = `Cell_type == Epithelial cells`. In primary tLung these are tS1/tS2/tS3; author `Malignant cells` appear in metastases / PE / tL-B only.",
        "- T/NK = `T lymphocytes` + `NK cells`.",
        "- Eligible sample: ≥20 epithelial and ≥20 T/NK cells.",
        f"- Genes missing from the UMI matrix: {missing or 'none'}.",
        "",
        "## Honest limits vs claim A3",
        "",
        "- A3 source claim is GSE207422 (neoadjuvant ICI, MPR vs NMPR). GSE131907 is a treatment-naive LUAD atlas.",
        "- No MPR / RECIST / PD-1 labels in the GEO series matrix.",
        f"- Primary tLung n={n_tlung_all} samples (eligible n={len(tlung)}); underpowered for a ρ≈−0.45 claim.",
        "- Metric is mean log1p(raw UMI), not author log2(TPM+1), because the 2.86 GB normalized text file was skipped.",
        "",
        "## Results",
        "",
    ]
    lines.extend(fmt(c) for c in contrasts)
    if paired.get("n"):
        lines.append(
            f"- {paired['contrast']}: n={paired['n']}, p={paired.get('wilcoxon_p', float('nan')):.3g}; "
            f"median %pos epithelial={paired.get('median_epi_pct_pos', float('nan')):.1f} vs "
            f"T/NK={paired.get('median_tnk_pct_pos', float('nan')):.1f}."
        )
    lines += [
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `feasibility.json` / `file_manifest.tsv` | Size budget and skip decisions |",
        "| `sample_metadata.tsv` | GEO patient / stage / origin |",
        "| `cell_type_composition.tsv` / `sample_composition.tsv` | Author labels |",
        "| `per_sample_tacstd2.tsv` | Sample-level epithelial and T/NK TACSTD2 |",
        "| `association_statistics.tsv` | Spearman tests |",
        "| `paired_compartment_test.json` / `sanity_checks.json` / `audit.json` | Extra tests |",
        "| `fig_epi_tacstd2_vs_tnk_fraction.png` | A3-style correlation panels |",
        "| `fig_tacstd2_epithelial_vs_tnk.png` | Tumor-restricted TACSTD2 |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 scripts/w200/A3_GSE131907/download.py",
        "python3 scripts/w200/A3_GSE131907/analyze.py",
        "```",
        "",
    ]
    (outdir / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
