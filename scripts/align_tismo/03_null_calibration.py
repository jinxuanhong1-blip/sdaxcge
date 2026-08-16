#!/usr/bin/env python3
"""Calibrate the Tacstd2 result against a genome-wide null and reference genes.

The headline design compares each cohort's ICB arm to its own baseline arm. If
the TISMO ICB arms are globally shifted relative to baseline arms (different
tumour purity, immune content, or residual normalisation effects), then *many*
genes would show the same "most cohorts go up" pattern, and a small p-value for
Tacstd2 alone would say little about Tacstd2.

This runs the identical cohort-level paired statistic on:
  * a seeded random panel of genes -> empirical null,
  * hand-picked reference genes: ICB pharmacodynamic markers that should move
    (Cd274, Pdcd1, Ifng, Cd8a) and epithelial markers that share Tacstd2's
    tumour-cell-fraction dependence (Epcam, Krt8).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc
from importlib import import_module

replicate = import_module("02_replicate")

GENE = "Tacstd2"


def stat_for_frame(raw: pd.DataFrame, cancer_types: dict[str, str]) -> dict | None:
    cohorts = replicate.build_cohort_table(raw, cancer_types)
    if cohorts.empty or len(cohorts) < 5:
        return None
    s = replicate.paired_summary(cohorts, "x")
    return {
        "n_cohorts": s["n_cohorts"],
        "n_up": s["n_up"],
        "frac_up": s["frac_up"],
        "mean_baseline": s["mean_baseline"],
        "mean_icb": s["mean_icb"],
        "mean_delta": s["mean_delta"],
        "paired_t_p": s["paired_t_p"],
        "wilcoxon_p": s["wilcoxon_p"],
        "cohens_dz": s["cohens_dz"],
    }


def collect(paths: dict[str, Path], cancer_types: dict[str, str]) -> pd.DataFrame:
    rows = []
    for gene, path in sorted(paths.items()):
        try:
            raw = tc.load_expression_csv(path, gene=gene)
        except Exception as err:  # noqa: BLE001 - skip unreadable payloads, keep the panel
            print(f"  skipping {gene}: {err}")
            continue
        s = stat_for_frame(raw, cancer_types)
        if s is None:
            continue
        rows.append({"gene": gene, **s})
    return pd.DataFrame(rows)


def main() -> int:
    tc.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    cancer_types = replicate.load_cell_line_cancer_types()

    target = stat_for_frame(
        pd.read_csv(tc.DATA_DIR / f"vivo_expression_{GENE}.csv"), cancer_types
    )
    assert target is not None

    null_dir = tc.ensure_null_genes()
    null_paths = {p.stem: p for p in sorted(null_dir.glob("*.csv"))} if null_dir.exists() else {}
    if not null_paths:
        print("No null panel found; run 01_download.py --null-genes N first.")
        return 1

    null_df = collect(null_paths, cancer_types)
    null_df.to_csv(tc.TABLES_DIR / "null_panel_statistics.csv", index=False)

    ref_dir = tc.DATA_DIR / "reference_genes"
    ref_paths = {
        p.stem.replace("vivo_expression_", ""): p for p in sorted(ref_dir.glob("*.csv"))
    } if ref_dir.exists() else {}
    ref_df = collect(ref_paths, cancer_types)
    if not ref_df.empty:
        ref_df.to_csv(tc.TABLES_DIR / "reference_gene_statistics.csv", index=False)

    n = len(null_df)
    # One-sided empirical percentiles: how extreme is Tacstd2 among random genes?
    frac_up_ge = float((null_df["n_up"] >= target["n_up"]).mean())
    delta_ge = float((null_df["mean_delta"] >= target["mean_delta"]).mean())
    wilcox_le = float((null_df["wilcoxon_p"] <= target["wilcoxon_p"]).mean())

    calib = {
        "gene": GENE,
        "target": target,
        "null_panel_size": n,
        "null_median_n_up": float(null_df["n_up"].median()),
        "null_mean_n_up": float(null_df["n_up"].mean()),
        "null_frac_of_genes_with_majority_up": float((null_df["frac_up"] > 0.5).mean()),
        "null_median_mean_delta": float(null_df["mean_delta"].median()),
        "null_frac_mean_delta_positive": float((null_df["mean_delta"] > 0).mean()),
        "null_frac_wilcoxon_p_below_0.05": float((null_df["wilcoxon_p"] < 0.05).mean()),
        "null_frac_wilcoxon_p_below_target": wilcox_le,
        "empirical_p_n_up": frac_up_ge,
        "empirical_p_mean_delta": delta_ge,
        # Bonferroni-style reference point for a 21,235-gene transcriptome.
        "target_wilcoxon_p_times_transcriptome": float(target["wilcoxon_p"] * 21235),
    }
    (tc.TABLES_DIR / "null_calibration.json").write_text(json.dumps(calib, indent=2))

    print("=" * 78)
    print(f"Null calibration for {GENE} ({n} random genes)")
    print("=" * 78)
    print(json.dumps(calib, indent=2))
    print()
    if not ref_df.empty:
        print("Reference genes (same cohort-level paired statistic):")
        cols = ["gene", "n_cohorts", "n_up", "mean_baseline", "mean_icb", "mean_delta",
                "paired_t_p", "wilcoxon_p"]
        show = pd.concat([ref_df, pd.DataFrame([{"gene": GENE, **target}])], ignore_index=True)
        print(show[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
