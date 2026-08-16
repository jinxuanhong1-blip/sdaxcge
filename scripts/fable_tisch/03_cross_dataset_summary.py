#!/usr/bin/env python3
"""Collect per-dataset stats_summary.json files into one cross-dataset table."""

from __future__ import annotations

import argparse
import glob
import json
import os

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", default="results/fable_tisch")
    args = ap.parse_args()

    rows = []
    for path in sorted(glob.glob(os.path.join(args.out_root, "*", "stats_summary.json"))):
        s = json.load(open(path))
        for gene, r in s["cell_level_malignant_vs_all_immune"].items():
            paired = s.get("per_patient_paired_wilcoxon", {}).get(gene, {})
            epi = s.get("epithelial_restriction", {}).get(gene, {})
            pos = s.get("positive_population_composition", {}).get(gene, {})
            ifr = s.get("immune_fraction_spearman", {}).get(gene, {})
            rows.append({
                "dataset": s["dataset"],
                "gene": gene,
                "malignant_definition": s.get("malignant_definition", ""),
                "n_malignant": s["n_malignant"],
                "n_immune": s["n_immune"],
                "mean_malignant": r["mean_malignant"],
                "mean_immune": r["mean_other"],
                "pct_pos_malignant": r["pct_pos_malignant"],
                "pct_pos_immune": r["pct_pos_other"],
                "log2FC_of_means": r["log2FC_of_means"],
                "auroc_malignant_vs_immune": r["auroc"],
                "cell_level_p": r["p_value"],
                "n_paired_patients": paired.get("n_patients"),
                "n_patients_malignant_higher": paired.get("n_patients_malignant_higher"),
                "paired_wilcoxon_p": paired.get("p_value"),
                "n_nonmalignant_epithelial": epi.get("n_nonmalignant_epithelial"),
                "mean_nonmalignant_epithelial": epi.get("mean_nonmalignant_epithelial"),
                "pct_pos_nonmalignant_epithelial": epi.get("pct_pos_nonmalignant_epithelial"),
                "log2FC_malignant_vs_nonmalig_epi": epi.get("log2FC_malignant_vs_nonmalig_epi"),
                "p_malignant_vs_nonmalig_epi": epi.get("p_malignant_vs_nonmalig_epi"),
                "pct_of_positive_in_malignant": pos.get("pct_of_positive_in_malignant"),
                "pct_of_positive_in_immune": pos.get("pct_of_positive_in_immune"),
                "spearman_rho_immune_frac_vs_malignant_mean": ifr.get("spearman_rho_immune_frac_vs_malignant_mean"),
                "spearman_p_immune_frac": ifr.get("p_value"),
            })
    df = pd.DataFrame(rows).sort_values(["gene", "dataset"])
    out_csv = os.path.join(args.out_root, "cross_dataset_summary.csv")
    df.to_csv(out_csv, index=False)
    with open(os.path.join(args.out_root, "cross_dataset_summary.md"), "w") as fh:
        fh.write(df.to_markdown(index=False, floatfmt=".3g"))
        fh.write("\n")
    print(df.to_string(index=False))
    print(f"\nwritten -> {out_csv}")

    # compact cross-dataset figure (TACSTD2 / CLDN4 only)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    sub = df[df.gene.isin(["TACSTD2", "CLDN4"])].copy()
    datasets = list(sub.dataset.unique())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(len(datasets))
    width = 0.35
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
        g = sub[sub.gene == gene].set_index("dataset").reindex(datasets)
        ax.bar(x - width / 2, g["pct_pos_malignant"], width, color="#c0392b", label="% pos malignant")
        ax.bar(x + width / 2, g["pct_pos_immune"], width, color="#2980b9", label="% pos immune")
        ax.set_xticks(x)
        ax.set_xticklabels(datasets, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel("% cells expressing")
        ax.set_ylim(0, 100)
        ax.set_title(gene)
        ax.legend(fontsize=8, frameon=False)
        for i, (au, lfc) in enumerate(zip(g["auroc_malignant_vs_immune"], g["log2FC_of_means"])):
            ax.text(i, 92, f"AUROC {au:.2f}\nlog2FC {lfc:.1f}", ha="center", va="top", fontsize=7)
    fig.suptitle("TACSTD2 / CLDN4: malignant vs immune (TISCH2 NSCLC objects)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_root, "cross_dataset_pct_pos.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
