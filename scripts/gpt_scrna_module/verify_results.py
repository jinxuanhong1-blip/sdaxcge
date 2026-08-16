#!/usr/bin/env python3
"""Independent invariants and statistic checks for generated result tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def residual(values: np.ndarray, purity: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), stats.rankdata(purity)])
    ranked = stats.rankdata(values)
    return ranked - design @ np.linalg.lstsq(design, ranked, rcond=None)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.results / "sample_scores.csv")
    reported = pd.read_csv(args.results / "association_statistics.csv")
    thesis = pd.read_csv(args.results / "thesis_results.csv")
    manifest = pd.read_csv(args.results / "source_manifest.csv")
    with (args.results / "audit.json").open() as handle:
        audit = json.load(handle)

    assert {"GSE131907", "GSE148071"}.issubset(set(samples["atlas"]))
    assert {"tacstd2", "cldn4", "tj_module", "immune_score", "immune_fraction"}.issubset(samples.columns)
    assert audit["GSE131907"]["matrix_cells"] == 208_506
    assert audit["GSE148071"]["matrix_cells"] == 89_887
    assert audit["GSE154826"]["status"] == "included_per_sample_geo_mtx"
    assert manifest.loc[manifest["used_for_statistics"], "bytes"].max() < 2_000_000_000
    assert not samples["tacstd2"].isna().all()
    assert not samples["cldn4"].isna().all()

    eligible = (
        (samples["n_epithelial"] >= 20)
        & (samples["n_T"] >= 20)
        & (samples["n_B"] >= 10)
        & (samples["n_myeloid"] >= 20)
    )
    for row in reported.itertuples(index=False):
        frame = samples[(samples["atlas"] == row.atlas) & eligible].dropna(
            subset=["epithelial_module", row.outcome]
        )
        if len(frame) < 4 or not np.isfinite(row.spearman_rho):
            continue
        observed = stats.spearmanr(frame["epithelial_module"], frame[row.outcome]).statistic
        assert len(frame) == row.n_samples
        assert np.isclose(observed, row.spearman_rho, atol=1e-12)
        assert 0 <= row.permutation_p <= 1
        assert row.bootstrap_ci_low <= row.spearman_rho <= row.bootstrap_ci_high

    for row in thesis.itertuples(index=False):
        assert row.match in {"match", "mismatch", "inconclusive"}
        assert row.metric in {"partial_spearman_rho", "OR", "log2FC", "spearman_rho"}
        if not np.isfinite(row.p):
            assert row.match == "inconclusive"
            continue
        assert 0 <= row.p <= 1
        if row.p >= 0.05:
            assert row.match == "inconclusive"
        frame = samples[(samples["atlas"] == row.atlas) & eligible]
        if row.metric == "spearman_rho" and row.outcome in frame.columns:
            usable = frame.dropna(subset=[row.exposure, row.outcome])
            observed = stats.spearmanr(usable[row.exposure], usable[row.outcome]).statistic
            assert np.isclose(observed, row.estimate, atol=1e-10)
        if row.metric == "log2FC":
            usable = frame.dropna(subset=[row.exposure, "immune_fraction"])
            high = usable[row.exposure] > usable[row.exposure].median()
            observed = np.log2(
                usable.loc[high, "immune_fraction"].mean()
                / usable.loc[~high, "immune_fraction"].mean()
            )
            assert np.isclose(observed, row.estimate, atol=1e-10)

    print("PASS: inputs, gene scores, exclusions, correlations, thesis table, and match labels verified")


if __name__ == "__main__":
    main()
