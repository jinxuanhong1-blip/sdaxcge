#!/usr/bin/env python3
"""Write FINDING + extra figures from saved LR/meta tables (no matrix reload)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze as A
from lib import fmt_num, fmt_p

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
GIVEN_RHO = -0.260
GIVEN_N = 95
GIVEN_P = 0.0195
MIN_K = 2
MIN_N = 15


def reportable(df: pd.DataFrame, n_col: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = df[(df["k_units"] >= MIN_K) & (df[n_col] >= MIN_N)].copy()
    return out.sort_values(["p_meta", n_col], ascending=[True, False])


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    coverage = pd.read_csv(RES / "patient_coverage.tsv", sep="\t")
    meta_q4 = pd.read_csv(RES / "lr_meta_q4q1.tsv", sep="\t")
    meta_med = pd.read_csv(RES / "lr_meta_median.tsv", sep="\t")
    meta_between = pd.read_csv(RES / "lr_meta_between_q4q1.tsv", sep="\t")
    pairs = pd.read_csv(RES / "patient_lr_long.tsv.gz", sep="\t")

    q4_rep = reportable(meta_q4, "n_patients")
    med_rep = reportable(meta_med, "n_patients")
    bet_rep = reportable(meta_between, "n_compared")
    q4_rep.to_csv(RES / "lr_meta_q4q1_reportable.tsv", sep="\t", index=False)
    med_rep.to_csv(RES / "lr_meta_median_reportable.tsv", sep="\t", index=False)
    bet_rep.to_csv(RES / "lr_meta_between_q4q1_reportable.tsv", sep="\t", index=False)

    skip = {
        "GSE207422": "DRMref barcodes not public; epithelial proxy on the same 12 locked samples",
        "GSE205335": "author malignant / T/NK",
        "GSE291670": "marker malignant; no author labels",
        "GSE253013": "9.3 GB RDS not downloaded; LR n=0",
        "GSE131907": "author malignant, locked n_mal≥20",
        "GSE325414": "author EpithelialCells_TumorCells",
    }
    A.plot_n(coverage, FIG / "fig_honest_n.png")
    A.plot_forest(q4_rep, FIG / "fig_forest_q4q1.png", "Within-patient Q4 vs Q1 outgoing ΔP (k≥2, n≥15)")
    A.plot_ligand_table(q4_rep, FIG / "fig_extra_ligand_table.png", "Extra: reportable outgoing patient-ΔP (Q4 vs Q1)")
    A.plot_forest(bet_rep, FIG / "fig_forest_between_q4q1.png", "Between-patient given-CLDN4 Q4 vs Q1 (k≥2, n≥15)")
    A.plot_forest(med_rep, FIG / "fig_forest_median.png", "Within-patient median-split outgoing ΔP (k≥2, n≥15)")
    if not q4_rep.empty:
        A.plot_patient_deltas(pairs, q4_rep.iloc[0]["interaction_name"], "q4q1", FIG / "fig_patient_delta_top.png")
        inhib = q4_rep[q4_rep["lr_class"].astype(str).str.contains("inhibitory|barrier|recruit")]
        if not inhib.empty:
            A.plot_ligand_table(
                inhib,
                FIG / "fig_extra_preclass.png",
                "Extra: barrier / inhibitory / recruit patient-ΔP",
            )

    A.write_finding(coverage, q4_rep, med_rep, bet_rep, skip, RES)

    summary = {
        "given_spearman": {"n": GIVEN_N, "rho": GIVEN_RHO, "p": GIVEN_P, "re_audited": False},
        "eligible_q4q1": int(coverage["eligible_q4q1"].sum()),
        "eligible_median": int(coverage["eligible_median"].sum()),
        "n_meta_q4_all": int(len(meta_q4)),
        "n_meta_q4_reportable": int(len(q4_rep)),
        "n_meta_between_reportable": int(len(bet_rep)),
        "top_q4": q4_rep.head(5)[["interaction_name", "n_patients", "mean_delta", "p_meta", "k_units"]].to_dict("records")
        if not q4_rep.empty
        else [],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
