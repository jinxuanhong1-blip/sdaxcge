#!/usr/bin/env python3
"""Fail if the Tacstd2 headline numbers drift from the replicated claim."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

EPS_MEAN = 0.02
EPS_P = 1e-6


def main() -> int:
    tab = pd.read_csv(tc.TABLES_DIR / "cohort_level_Tacstd2.csv")
    s = pd.read_csv(tc.TABLES_DIR / "paired_summaries_Tacstd2.csv")
    head = s[s["subset"] == "ALL cohorts (headline design)"].iloc[0]
    lung = s[s["subset"] == "Lung cohorts only"].iloc[0]
    calib = json.loads((tc.TABLES_DIR / "null_calibration.json").read_text())
    errors: list[str] = []

    def check(name: str, got, want, tol=0):
        if abs(float(got) - float(want)) > tol:
            errors.append(f"{name}: got {got}, expected {want} ± {tol}")

    check("n_cohorts", head.n_cohorts, 64)
    check("n_up", head.n_up, 49)
    check("n_down", head.n_down, 15)
    check("mean_baseline", head.mean_baseline, 1.051, EPS_MEAN)
    check("mean_icb", head.mean_icb, 1.317, EPS_MEAN)
    check("wilcoxon_p", head.wilcoxon_p, 5.84e-5, EPS_P)
    check("lung_n", lung.n_cohorts, 2)
    check("lung_cell_lines", lung.n_cell_lines, 1)
    check("lung_up", lung.n_up, 2)
    if set(tab.loc[tab.cancer_type.str.contains("Lung", case=False), "cell_line"]) != {"LLC"}:
        errors.append("lung cell lines are not exactly {LLC}")
    check("null_panel_size", calib["null_panel_size"], 500)
    if errors:
        print("HEADLINE ASSERTIONS FAILED")
        print("\n".join(f"  - {e}" for e in errors))
        return 1
    print("headline assertions OK: 64 cohorts, 49/64 up, 1.05→1.32, Wilcoxon p=5.84e-5; lung=2 LLC")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
