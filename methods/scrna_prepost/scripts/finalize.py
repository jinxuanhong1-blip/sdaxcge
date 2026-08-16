#!/usr/bin/env python3
"""Rebuild paper tables from an existing GSE207422 per-sample file + GSE337519."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze import (
    GENES,
    MIN_MAL,
    analyze_gse337519,
    mw,
    stouffer,
)

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "tables"
ROOT = HERE.parents[2]


def tests_from_samples(sample_df: pd.DataFrame) -> list[dict]:
    mal_ok = sample_df[sample_df["n_malignant"] >= MIN_MAL].copy()
    epi_ok = sample_df[sample_df["n_epithelial"] >= MIN_MAL].copy()
    contrasts = [
        ("epi_all_post_vs_pre", epi_ok, "epi", "unmatched 3 pre vs 12 post; epithelial (CopyKAT not deposited)"),
        (
            "epi_nmpr_post_vs_pre",
            epi_ok[epi_ok["timepoint"].eq("pre") | epi_ok["response_paper"].eq("NMPR")],
            "epi",
            "unmatched pre vs NMPR post; epithelial",
        ),
        (
            "mal_all_post_vs_pre",
            mal_ok,
            "mal",
            "unmatched; ≥10 malignant-like (epithelial minus alveolar/club/ciliated)",
        ),
        (
            "mal_nmpr_post_vs_pre",
            mal_ok[mal_ok["timepoint"].eq("pre") | mal_ok["response_paper"].eq("NMPR")],
            "mal",
            "unmatched pre vs NMPR post; ≥10 malignant-like",
        ),
    ]
    tests = []
    for contrast, sub, prefix, note in contrasts:
        pre = sub[sub["timepoint"] == "pre"]
        post = sub[sub["timepoint"] == "post"]
        for gene in GENES:
            for metric, suffix in [
                ("mean_log1p_cp10k", "mean_log1p"),
                ("pct_pos", "pct_pos"),
                ("pb_cpm", "pb_cpm"),
            ]:
                col = f"{prefix}_{gene}_{suffix}"
                if col not in sub.columns:
                    continue
                stat = mw(pre[col], post[col])
                stat.update(
                    {
                        "accession": "GSE207422",
                        "contrast": contrast,
                        "gene": gene,
                        "metric": metric,
                        "compartment": "epithelial" if prefix == "epi" else "malignant-like",
                        "pairing": "unmatched",
                        "n_pre_patients": int(pre["Patient"].nunique()),
                        "n_post_patients": int(post["Patient"].nunique()),
                        "pre_patients": ",".join(pre["Patient"].astype(str)),
                        "post_patients": ",".join(post["Patient"].astype(str)),
                        "design_note": note,
                    }
                )
                tests.append(stat)
    return tests, epi_ok


def main() -> None:
    sample_df = pd.read_csv(OUT / "gse207422_per_sample.tsv", sep="\t")
    tests, epi_ok = tests_from_samples(sample_df)
    pd.DataFrame(tests).to_csv(OUT / "gse207422_pre_vs_post.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    rng = np.random.default_rng(1)
    for ax, gene in zip(axes, GENES):
        col = f"epi_{gene}_mean_log1p"
        pre_v = epi_ok.loc[epi_ok["timepoint"] == "pre", col]
        post_v = epi_ok.loc[epi_ok["timepoint"] == "post", col]
        ax.boxplot([pre_v, post_v], tick_labels=["pre", "post"], widths=0.45)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(pre_v)), pre_v, c="#1f77b4", s=36, zorder=3)
        colors = epi_ok.loc[epi_ok["timepoint"] == "post", "response_paper"].map(
            {"MPR": "#2ca02c", "NMPR": "#d62728", "NE": "#7f7f7f"}
        )
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(post_v)), post_v, c=list(colors), s=36, zorder=3)
        stat = mw(pre_v, post_v)
        ax.set_title(f"{gene}\nMWU p={stat['p']:.3g}")
        ax.set_ylabel("epithelial mean log1p(CP10K)")
    fig.suptitle("GSE207422 unmatched pre vs post (sample unit, epithelial, 3 vs 12)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse207422_epithelial_pre_post.png", dpi=160)
    plt.close(fig)

    g337 = analyze_gse337519(ROOT / "data/scrna_prepost/GSE337519/raw", OUT)

    paper_rows = []
    for gene in GENES:
        for contrast, metric_label, note in [
            (
                "epi_all_post_vs_pre",
                "epithelial mean log1p(CP10K)",
                "3 pre vs 12 post; sample unit; CopyKAT not on GEO",
            ),
            (
                "epi_nmpr_post_vs_pre",
                "epithelial mean log1p(CP10K)",
                "3 pre vs 8 NMPR post (residual tumor)",
            ),
            (
                "mal_all_post_vs_pre",
                "malignant-like mean log1p(CP10K)",
                "≥10 cells after dropping alveolar/club/ciliated epithelial",
            ),
        ]:
            row = next(r for r in tests if r["gene"] == gene and r["contrast"] == contrast and r["metric"] == "mean_log1p_cp10k")
            paper_rows.append(
                {
                    "accession": "GSE207422",
                    "citation": "Hu et al. Genome Med 2023",
                    "regimen": "neoadjuvant PD-1 + chemo",
                    "pairing": "unmatched (different patients)",
                    "contrast": contrast,
                    "n_pre": row["n_pre"],
                    "n_post": row["n_post"],
                    "gene": gene,
                    "metric": metric_label,
                    "median_pre": row["median_pre"],
                    "median_post": row["median_post"],
                    "direction": row["direction"],
                    "test": row["test"],
                    "U": row["U"],
                    "p": row["p"],
                    "note": note,
                }
            )
    pd.DataFrame(paper_rows).to_csv(OUT / "paper_table_scrna_prepost.tsv", sep="\t", index=False)

    combined = []
    for gene in GENES:
        gene_rows = [r for r in paper_rows if r["gene"] == gene and r["contrast"] == "epi_all_post_vs_pre"]
        s = stouffer(gene_rows)
        s["gene"] = gene
        s["metric"] = "epithelial mean log1p(CP10K)"
        s["studies"] = "GSE207422"
        combined.append(s)
    pd.DataFrame(combined).to_csv(OUT / "combined_stouffer.tsv", sep="\t", index=False)

    dropped = sample_df.loc[sample_df["n_malignant"] < MIN_MAL, ["Patient", "timepoint", "response_paper", "n_epithelial", "n_malignant"]]
    summary = {
        "testable_tumor_cell_cohorts": ["GSE207422"],
        "primary": "epithelial 3 pre vs 12 post, unmatched MWU",
        "malignant_like_floor": MIN_MAL,
        "dropped_below_10_malignant": dropped.to_dict(orient="records"),
        "combined": combined,
        "gse337519": g337,
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print(pd.DataFrame(paper_rows).to_string(index=False))


if __name__ == "__main__":
    main()
