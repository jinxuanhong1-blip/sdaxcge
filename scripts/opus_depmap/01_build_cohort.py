#!/usr/bin/env python3
"""Build the lung cell-line cohort and cache the DepMap matrices as parquet.

Writes:
  results/opus_depmap/tables/cohort_lung_models.csv
  results/opus_depmap/tables/cohort_summary.csv
  notes/opus_depmap/01_cohort_qc.md
"""

from __future__ import annotations

import json
import sys

import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DATA,
    GENES_OF_INTEREST,
    NOTES,
    RELEASE,
    TABLES,
    ensure_dirs,
    load_expression,
    load_gene_dependency,
    load_gene_effect,
    load_models,
    lung_cohort,
)


def main() -> int:
    ensure_dirs()
    manifest = json.loads((DATA / "MANIFEST.json").read_text())

    models = load_models()
    lung = lung_cohort(models)

    print("loading gene effect ...", flush=True)
    gene_effect = load_gene_effect()
    print("loading gene dependency ...", flush=True)
    gene_dep = load_gene_dependency()
    print("loading expression ...", flush=True)
    expr = load_expression()

    lung_ge = lung.index.intersection(gene_effect.index)
    lung_ex = lung.index.intersection(expr.index)
    lung_both = lung_ge.intersection(lung_ex)

    cohort = lung.loc[
        :,
        [
            "CellLineName",
            "StrippedCellLineName",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "OncotreeCode",
            "LungGroup",
            "PrimaryOrMetastasis",
            "Sex",
            "Age",
        ],
    ].copy()
    cohort["HasCRISPR"] = cohort.index.isin(gene_effect.index)
    cohort["HasExpression"] = cohort.index.isin(expr.index)
    for gene in GENES_OF_INTEREST:
        cohort[f"{gene}_geneEffect"] = gene_effect[gene].reindex(cohort.index)
        cohort[f"{gene}_depProb"] = gene_dep[gene].reindex(cohort.index)
        cohort[f"{gene}_log2TPM1"] = expr[gene].reindex(cohort.index)
    cohort = cohort.sort_values(["LungGroup", "StrippedCellLineName"])
    cohort.to_csv(TABLES / "cohort_lung_models.csv")

    summary_rows = [
        {"metric": "release", "value": RELEASE},
        {"metric": "figshare_doi", "value": manifest["figshare_doi"]},
        {"metric": "models_total", "value": len(models)},
        {"metric": "lung_lineage_cancer_models", "value": len(lung)},
        {"metric": "lung_with_CRISPR", "value": len(lung_ge)},
        {"metric": "lung_with_expression", "value": len(lung_ex)},
        {"metric": "lung_with_CRISPR_and_expression", "value": len(lung_both)},
        {"metric": "all_models_with_CRISPR", "value": int(gene_effect.shape[0])},
        {"metric": "genes_in_CRISPR_matrix", "value": int(gene_effect.shape[1])},
        {"metric": "all_models_with_expression", "value": int(expr.shape[0])},
        {"metric": "genes_in_expression_matrix", "value": int(expr.shape[1])},
    ]
    for group, count in lung["LungGroup"].value_counts().items():
        summary_rows.append({"metric": f"lung_{group}_models", "value": int(count)})
        sub = lung.index[lung["LungGroup"] == group]
        summary_rows.append(
            {
                "metric": f"lung_{group}_with_CRISPR",
                "value": int(sub.isin(gene_effect.index).sum()),
            }
        )
        summary_rows.append(
            {
                "metric": f"lung_{group}_with_expression",
                "value": int(sub.isin(expr.index).sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES / "cohort_summary.csv", index=False)

    subtype_counts = (
        lung.loc[lung_ge, "OncotreeSubtype"].value_counts().rename("n_with_CRISPR")
    )

    lines = [
        "# 01 · Cohort and data QC",
        "",
        f"Release: **{RELEASE}** (figshare DOI `{manifest['figshare_doi']}`, "
        f"retrieved {manifest['retrieved_utc']}). All files are open release files; "
        "md5 checksums are recorded in `data/opus_depmap/MANIFEST.json`.",
        "",
        "## Cohort definition",
        "",
        "* Lung lineage = `Model.csv` rows with `OncotreeLineage == \"Lung\"`.",
        "* Non-cancerous lung models (`OncotreePrimaryDisease == \"Non-Cancerous\"`) are dropped.",
        "* `LungGroup` collapses Oncotree annotations into NSCLC / SCLC / Other lung "
        "(the latter is mostly SMARCA4-deficient undifferentiated thoracic tumours and "
        "non-SCLC neuroendocrine models).",
        "",
        "## Counts",
        "",
        "| metric | value |",
        "| --- | --- |",
    ]
    for row in summary_rows:
        lines.append(f"| {row['metric']} | {row['value']} |")
    lines += [
        "",
        "## Lung subtypes with CRISPR screens",
        "",
        "| OncotreeSubtype | n with CRISPR |",
        "| --- | --- |",
    ]
    for subtype, count in subtype_counts.items():
        lines.append(f"| {subtype} | {count} |")
    lines += [
        "",
        "## Measurement notes",
        "",
        "* `CRISPRGeneEffect.csv` is the integrated Chronos gene effect matrix, scaled so "
        "that the median common essential is -1 and the median non-essential is 0.",
        "* `CRISPRGeneDependency.csv` gives the posterior probability that a model is "
        "dependent on a gene; DepMap's convention is that probability > 0.5 counts as a "
        "dependent line.",
        "* Expression is `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, i.e. "
        "log2(TPM+1) from the GTEx-style RNA-seq pipeline, model-level.",
        "* Duplicate gene symbols (same symbol, different Entrez id) are de-duplicated by "
        "keeping the first occurrence; neither TACSTD2 nor CLDN4 is affected.",
    ]
    (NOTES / "01_cohort_qc.md").write_text("\n".join(lines) + "\n")

    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
