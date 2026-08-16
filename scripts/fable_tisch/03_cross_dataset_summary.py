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
            rows.append({
                "dataset": s["dataset"],
                "gene": gene,
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
            })
    df = pd.DataFrame(rows).sort_values(["gene", "dataset"])
    out_csv = os.path.join(args.out_root, "cross_dataset_summary.csv")
    df.to_csv(out_csv, index=False)
    with open(os.path.join(args.out_root, "cross_dataset_summary.md"), "w") as fh:
        fh.write(df.to_markdown(index=False, floatfmt=".3g"))
        fh.write("\n")
    print(df.to_string(index=False))
    print(f"\nwritten -> {out_csv}")


if __name__ == "__main__":
    main()
