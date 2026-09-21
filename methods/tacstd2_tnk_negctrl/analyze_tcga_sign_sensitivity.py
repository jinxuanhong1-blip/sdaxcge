#!/usr/bin/env python3
"""Descriptive TCGA sensitivity: cohorts whose unadjusted TACSTD2–T/NK rho is negative.

Not the pre-specified test. The pre-specified meta is the seven-cohort funnel
in analyze_tcga.py, which also includes CESC (positive) and KIRC (about zero).

This script rebuilds the DerSimonian–Laird table from the saved patient scores.
The label-swap p-values are already in tcga_negcohort_head_to_head.tsv
(10,000 within-patient swaps, seed 1). Pass --perm to recompute them.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_tcga import CONTROLS, label_swap_p
from stats import bh_fdr, dl_meta, partial_spearman, spearman

ROOT = Path(__file__).resolve().parent
TAB = ROOT / "results" / "tables"
NEG = ["LUAD", "BRCA", "STAD", "BLCA", "PAAD"]


def meta_table(patient: pd.DataFrame) -> pd.DataFrame:
    rows = []
    frames = [patient[patient.cohort == c] for c in NEG]
    unadj_rho, unadj_n = [], []
    for frame in frames:
        x = frame["TACSTD2"].to_numpy(float)
        y = frame["tnk"].to_numpy(float)
        rho = spearman(x, y)
        unadj_rho.append(rho)
        unadj_n.append(len(frame))
    base = dl_meta(unadj_rho, unadj_n, 0)
    for gene in ["none", "CLDN4", "CLDN3", "CLDN7", "EPCAM", "MUC1", "KRT19"]:
        rhos, ns = [], []
        for frame in frames:
            x = frame["TACSTD2"].to_numpy(float)
            y = frame["tnk"].to_numpy(float)
            if gene == "none":
                rho = spearman(x, y)
                k = 0
            else:
                rho = partial_spearman(x, y, [frame[gene].to_numpy(float)])
                k = 1
            rhos.append(rho)
            ns.append(len(frame))
        meta = dl_meta(rhos, ns, k)
        rows.append(
            {
                "block": "funnel_negative_unadj",
                "cohorts": ",".join(NEG),
                "outcome": "tnk",
                "conditioner": gene,
                "rho_unadj": base["rho"],
                "p_unadj": base["p"],
                "rho_partial": meta["rho"],
                "p_partial": meta["p"],
                "I2": meta["I2"],
                "ci_lo": meta["ci_lo"],
                "ci_hi": meta["ci_hi"],
                "attenuation": meta["rho"] - base["rho"],
                "k": meta["k"],
                "N": meta["N"],
                "note": "Cohorts selected on a negative unadjusted rho. Not pre-specified.",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--perm", action="store_true")
    args = parser.parse_args()
    patient = pd.read_csv(TAB / "tcga_patient_scores.tsv.gz", sep="\t")
    meta_table(patient).to_csv(TAB / "tcga_negcohort_meta.tsv", sep="\t", index=False)
    if not args.perm:
        print("wrote meta; skipped label-swap (use --perm)")
        return
    frames = [patient[patient.cohort == c].reset_index(drop=True) for c in NEG]
    rng = np.random.default_rng(1)
    rows = []
    for control in CONTROLS:
        swap = label_swap_p(frames, "tnk", control, rng, 10000)
        rows.append(
            {
                "block": "funnel_negative_unadj",
                "outcome": "tnk",
                "control": control,
                "meta_partial_CLDN4": swap["partial_cldn4"],
                "meta_partial_control": swap["partial_control"],
                "delta_partial_CLDN4_minus_control": swap["delta"],
                "perm_p_two_sided": swap["p"],
                "n_perm": 10000,
                "cohorts": ",".join(NEG),
                "note": "Descriptive. Cohorts chosen because unadjusted TACSTD2–T/NK rho was negative. Not the pre-specified seven-cohort test.",
            }
        )
    q = bh_fdr([r["perm_p_two_sided"] for r in rows])
    for row, qq in zip(rows, q):
        row["q_bh"] = qq
    pd.DataFrame(rows).to_csv(TAB / "tcga_negcohort_head_to_head.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
