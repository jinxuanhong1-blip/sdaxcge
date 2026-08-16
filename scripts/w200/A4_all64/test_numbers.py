#!/usr/bin/env python3
"""Lock the A4 recompute against the committed TISMO table."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path("results/w200/A4_all64")


def main() -> None:
    summary = json.loads((ROOT / "summary.json").read_text())
    models = pd.read_csv(ROOT / "per_model.tsv", sep="\t")
    split = pd.read_csv(ROOT / "cancer_type_split.tsv", sep="\t")
    raw = pd.read_csv(ROOT / "tacstd2_vivo_icb.csv")

    assert len(models) == 64
    assert raw["cell_line"].nunique() == 64
    assert raw["geneID"].nunique() == 1
    assert raw["geneID"].iloc[0] == "Tacstd2"

    med = summary["primary_median_sign_test"]
    assert med["n_up"] == 40 and med["n_down"] == 18 and med["n_tie"] == 6
    assert abs(med["sign_p_two_sided"] - 5.354859650800294e-3) < 1e-12

    mean = summary["mean_sign_test_user_count"]
    assert mean["n_up"] == 49 and mean["n_down"] == 15
    assert abs(mean["sign_p_two_sided"] - 2.436457216293327e-5) < 1e-15

    w = summary["wilcoxon_mean_deltas_user_p"]
    assert abs(w["p_two_sided"] - 5.8398482923347564e-5) < 1e-15

    mam = split[(split["cancer_group"] == "Mammary") & (split["location"] == "median")].iloc[0]
    assert int(mam["n_up"]) == 16 and int(mam["n_down"]) == 12 and int(mam["n_tie"]) == 1
    print("ok")


if __name__ == "__main__":
    main()
