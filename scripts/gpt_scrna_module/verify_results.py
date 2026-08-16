#!/usr/bin/env python3
"""Independent invariants and statistic checks for generated result tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.results / "sample_scores.csv")
    reported = pd.read_csv(args.results / "association_statistics.csv")
    manifest = pd.read_csv(args.results / "source_manifest.csv")
    with (args.results / "audit.json").open() as handle:
        audit = json.load(handle)

    assert set(samples["atlas"]) == {"GSE131907", "GSE148071"}
    assert len(samples[samples["atlas"] == "GSE131907"]) == 21
    assert len(samples[samples["atlas"] == "GSE148071"]) == 42
    assert audit["GSE131907"]["matrix_cells"] == 208_506
    assert audit["GSE131907"]["tumor_bearing_specimens"] == 21
    assert audit["GSE148071"]["matrix_cells"] == 89_887
    assert audit["GSE154826"]["status"] == "excluded_before_analysis"
    assert audit["GSE154826"]["author_hca_rds_bytes"] > 2_000_000_000
    assert manifest.loc[manifest["used_for_statistics"], "bytes"].max() < 2_000_000_000

    eligible = (
        (samples["n_epithelial"] >= 20)
        & (samples["n_T"] >= 20)
        & (samples["n_B"] >= 10)
        & (samples["n_myeloid"] >= 20)
    )
    for row in reported.itertuples(index=False):
        frame = samples[(samples["atlas"] == row.atlas) & eligible]
        observed = stats.spearmanr(frame["epithelial_module"], frame[row.outcome]).statistic
        assert len(frame) == row.n_samples
        assert np.isclose(observed, row.spearman_rho, atol=1e-12)
        assert 0 <= row.permutation_p <= 1
        assert row.bootstrap_ci_low <= row.spearman_rho <= row.bootstrap_ci_high
        assert row.permutation_p <= row.bh_q_within_atlas <= 1

    assert not reported["bh_q_within_atlas"].lt(0.05).any()
    print("PASS: inputs, sample counts, exclusions, correlations, intervals, and FDR verified")


if __name__ == "__main__":
    main()
