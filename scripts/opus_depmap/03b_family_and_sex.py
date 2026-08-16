#!/usr/bin/env python3
"""Targeted follow-ups that do not need a genome-wide scan.

1. Family / epithelial paralog co-dependency (CLDN3/7, EPCAM, CDH1, keratins).
2. Donor-sex check, because Y-linked genes appeared among the lung-only top hits
   of the genome-wide scan (those hits are not trusted until this is ruled out).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    GENES_OF_INTEREST,
    NOTES,
    TABLES,
    bh_fdr,
    ensure_dirs,
    load_gene_effect,
    load_models,
    lung_cohort,
)

FAMILIES = {
    "TACSTD2": ["EPCAM", "CLDN4", "CLDN3", "CLDN7", "KRT8", "KRT18", "CDH1"],
    "CLDN4": ["CLDN3", "CLDN7", "CLDN1", "CLDN6", "TACSTD2", "EPCAM", "CDH1"],
}
Y_GENES = ["UTY", "USP9Y", "DDX3Y", "KDM5D", "RPS4Y1", "EIF1AY", "ZFY"]


def main() -> int:
    ensure_dirs()
    models = load_models()
    lung = lung_cohort(models)
    gene_effect = load_gene_effect()
    lung_ids = lung.index.intersection(gene_effect.index)

    targeted = []
    for query, partners in FAMILIES.items():
        for cohort_mat, label in (
            (gene_effect, "pan-cancer"),
            (gene_effect.loc[lung_ids], "lung"),
        ):
            for partner in partners:
                if partner not in cohort_mat.columns or partner == query:
                    continue
                pair = cohort_mat[[query, partner]].dropna()
                if len(pair) < 20:
                    continue
                r, p = stats.pearsonr(pair[query], pair[partner])
                targeted.append(
                    {
                        "query": query,
                        "partner": partner,
                        "cohort": label,
                        "n": int(len(pair)),
                        "pearson_r": float(r),
                        "p": float(p),
                    }
                )
    targeted = pd.DataFrame(targeted)
    targeted["q_within_family_test"] = bh_fdr(targeted["p"].to_numpy())
    targeted.to_csv(TABLES / "codependency_targeted_family.csv", index=False)

    y_genes = [g for g in Y_GENES if g in gene_effect.columns]
    sex = models["Sex"].reindex(lung_ids)
    sex_rows = []
    for gene in GENES_OF_INTEREST:
        prof = gene_effect.loc[lung_ids, gene]
        male = prof[sex == "Male"].dropna()
        female = prof[sex == "Female"].dropna()
        y_rs = []
        for yg in y_genes:
            pair = gene_effect.loc[lung_ids, [gene, yg]].dropna()
            if len(pair) >= 20:
                y_rs.append(float(pair.corr().iloc[0, 1]))
        sex_rows.append(
            {
                "gene": gene,
                "n_male": int(len(male)),
                "n_female": int(len(female)),
                "mean_male": float(male.mean()) if len(male) else np.nan,
                "mean_female": float(female.mean()) if len(female) else np.nan,
                "mannwhitney_p": (
                    float(stats.mannwhitneyu(male, female, alternative="two-sided").pvalue)
                    if len(male) > 5 and len(female) > 5
                    else np.nan
                ),
                "mean_r_with_Y_linked_gene_effects": float(np.mean(y_rs)) if y_rs else np.nan,
                "n_Y_genes": len(y_rs),
                "y_genes_used": ";".join(y_genes),
            }
        )
    sex_frame = pd.DataFrame(sex_rows)
    sex_frame.to_csv(TABLES / "codependency_sex_confounder_check.csv", index=False)

    lines = [
        "# 03b · Targeted family co-dependency and donor-sex check",
        "",
        "These tests are pre-specified and do not use the genome-wide scan.",
        "",
        "## Family / epithelial partners",
        "",
        "| query | partner | cohort | n | Pearson r | p | q (within this table) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in targeted.sort_values(["query", "cohort", "p"]).itertuples():
        lines.append(
            f"| {r.query} | {r.partner} | {r.cohort} | {r.n} | {r.pearson_r:.3f} | "
            f"{r.p:.3g} | {r.q_within_family_test:.3g} |"
        )
    lines += [
        "",
        "## Donor sex vs lung gene effect",
        "",
        "| gene | n male | n female | mean male | mean female | Mann-Whitney p | mean r vs Y-linked gene effects |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in sex_frame.itertuples():
        lines.append(
            f"| {r.gene} | {r.n_male} | {r.n_female} | {r.mean_male:.4f} | {r.mean_female:.4f} | "
            f"{r.mannwhitney_p:.3g} | {r.mean_r_with_Y_linked_gene_effects:.3f} |"
        )
    (NOTES / "03b_family_and_sex.md").write_text("\n".join(lines) + "\n")
    print(targeted.to_string(index=False))
    print()
    print(sex_frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
