#!/usr/bin/env python3
"""TACSTD2 / CLDN4 partial association and product-method ACME vs T/NK.

Input is the locked concordant-4 patient table (n = 65): malignant
percent detected and malignant mean log1p(UMI) for TACSTD2 and CLDN4,
and frac_tnk. The CLDN4 % positive meta-Spearman is a pipeline check
against the locked value ρ = −0.5311678045689989.

Primary pool is DerSimonian–Laird within the four cohorts. A single
Spearman on all 65 rows is written as a sensitivity only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from stats import (  # noqa: E402
    acme_bootstrap,
    attenuation_percent,
    dl_means,
    fisher_dl,
    partial_spearman,
    spearman,
    spearman_ci,
    verdict_row,
)

TABLE = ROOT / "data" / "concordant4_patient_scores.tsv"
OUT = ROOT / "results" / "tables"
LOCKED_CLDN4_PCT_RHO = -0.5311678045689989
READOUTS = (
    ("pct_pos", "pct_TACSTD2", "pct_CLDN4"),
    ("mean_log1p_umi", "mean_TACSTD2", "mean_CLDN4"),
)
COHORTS = ("GSE123902", "GSE131907", "GSE205335", "GSE189357")


def one_fit(x, m, y, seed: int) -> dict:
    raw_x, p_x = spearman(x, y)
    raw_m, p_m = spearman(m, y)
    px = partial_spearman(x, y, m)
    pm = partial_spearman(m, y, x)
    med = acme_bootstrap(x, m, y, n_boot=4000, seed=seed)
    lo, hi = spearman_ci(raw_x, len(x))
    att = attenuation_percent(raw_x, px["rho"])
    call = verdict_row(
        raw_x, lo, hi, px["rho"], px["p"],
        med["acme"], med["acme_lo"], med["acme_hi"], pm["rho"],
    )
    return {
        "n": int(len(x)),
        "rho_tacstd2": raw_x,
        "p_tacstd2": p_x,
        "rho_tacstd2_lo": lo,
        "rho_tacstd2_hi": hi,
        "rho_cldn4": raw_m,
        "p_cldn4": p_m,
        "rho_tacstd2_given_cldn4": px["rho"],
        "p_tacstd2_given_cldn4": px["p"],
        "rho_tacstd2_given_cldn4_rankresid": px["rank_residual_rho"],
        "rho_cldn4_given_tacstd2": pm["rho"],
        "p_cldn4_given_tacstd2": pm["p"],
        "rho_cldn4_given_tacstd2_rankresid": pm["rank_residual_rho"],
        "rho_tacstd2_cldn4": px["rho_xz"],
        "attenuation_pct": att,
        "acme": med["acme"],
        "acme_lo": med["acme_lo"],
        "acme_hi": med["acme_hi"],
        "acme_se": med["acme_se"],
        "acme_p": med["acme_p"],
        "ade": med["ade"],
        "total_linear": med["total"],
        "a_x_to_m": med["a"],
        "b_m_to_y": med["b"],
        "verdict": call,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    patient = pd.read_csv(TABLE, sep="\t")
    if patient["unit_id"].nunique() != 65 or set(patient["dataset"]) != set(COHORTS):
        raise SystemExit(f"patient table is not the locked n=65 concordant-4: {patient.groupby('dataset').size().to_dict()}")

    rows = []
    seed = 4639
    for ds in COHORTS:
        g = patient[patient.dataset == ds]
        for readout, xc, mc in READOUTS:
            seed += 1
            rec = one_fit(g[xc], g[mc], g["frac_tnk"], seed=seed)
            rec.update({"dataset": ds, "readout": readout, "level": "cohort"})
            rows.append(rec)

    cohort = pd.DataFrame(rows)
    # Pipeline check: CLDN4 %pos meta must reproduce the locked Spearman.
    pct = cohort[cohort.readout == "pct_pos"].set_index("dataset").loc[list(COHORTS)]
    check = fisher_dl(pct.rho_cldn4, pct.n, partial=False)
    if abs(check["rho"] - LOCKED_CLDN4_PCT_RHO) > 1e-12:
        raise SystemExit(f"CLDN4 %pos pipeline check failed: {check['rho']} != {LOCKED_CLDN4_PCT_RHO}")

    meta_rows = []
    for readout, _, _ in READOUTS:
        sub = cohort[cohort.readout == readout].set_index("dataset").loc[list(COHORTS)]
        raw = fisher_dl(sub.rho_tacstd2, sub.n, partial=False)
        cld = fisher_dl(sub.rho_cldn4, sub.n, partial=False)
        par_x = fisher_dl(sub.rho_tacstd2_given_cldn4, sub.n, partial=True)
        par_m = fisher_dl(sub.rho_cldn4_given_tacstd2, sub.n, partial=True)
        ac = dl_means(sub.acme, sub.acme_se)
        att = attenuation_percent(raw["rho"], par_x["rho"])
        call = verdict_row(
            raw["rho"], raw["ci_lo"], raw["ci_hi"], par_x["rho"], par_x["p"],
            ac["estimate"], ac["ci_lo"], ac["ci_hi"], par_m["rho"],
        )
        meta_rows.append({
            "dataset": "concordant-4",
            "readout": readout,
            "level": "dl_meta",
            "n": int(sub.n.sum()),
            "k": int(raw["k"]),
            "rho_tacstd2": raw["rho"],
            "p_tacstd2": raw["p"],
            "rho_tacstd2_lo": raw["ci_lo"],
            "rho_tacstd2_hi": raw["ci_hi"],
            "I2_tacstd2": raw["I2"],
            "rho_cldn4": cld["rho"],
            "p_cldn4": cld["p"],
            "rho_cldn4_lo": cld["ci_lo"],
            "rho_cldn4_hi": cld["ci_hi"],
            "I2_cldn4": cld["I2"],
            "rho_tacstd2_given_cldn4": par_x["rho"],
            "p_tacstd2_given_cldn4": par_x["p"],
            "rho_tacstd2_given_cldn4_lo": par_x["ci_lo"],
            "rho_tacstd2_given_cldn4_hi": par_x["ci_hi"],
            "I2_tacstd2_given_cldn4": par_x["I2"],
            "rho_cldn4_given_tacstd2": par_m["rho"],
            "p_cldn4_given_tacstd2": par_m["p"],
            "rho_cldn4_given_tacstd2_lo": par_m["ci_lo"],
            "rho_cldn4_given_tacstd2_hi": par_m["ci_hi"],
            "attenuation_pct": att,
            "acme": ac["estimate"],
            "acme_lo": ac["ci_lo"],
            "acme_hi": ac["ci_hi"],
            "acme_p": ac["p"],
            "I2_acme": ac["I2"],
            "verdict": call,
            "cldn4_pipeline_check_rho": check["rho"] if readout == "pct_pos" else np.nan,
        })

    # Sensitivity: ignore cohort and correlate all 65 patients.
    sens = []
    for readout, xc, mc in READOUTS:
        seed += 1
        rec = one_fit(patient[xc], patient[mc], patient["frac_tnk"], seed=seed)
        rec.update({"dataset": "concordant-4-pooled-rows", "readout": readout, "level": "sensitivity_pooled_rows"})
        sens.append(rec)

    # Leave-one-patient influence on cohort ACME for the two smallest cohorts.
    influence = []
    for ds in ("GSE123902", "GSE189357"):
        g = patient[patient.dataset == ds].reset_index(drop=True)
        for readout, xc, mc in READOUTS:
            full = one_fit(g[xc], g[mc], g["frac_tnk"], seed=7)
            for i in range(len(g)):
                keep = np.ones(len(g), dtype=bool)
                keep[i] = False
                held = one_fit(g.loc[keep, xc], g.loc[keep, mc], g.loc[keep, "frac_tnk"], seed=7)
                influence.append({
                    "dataset": ds,
                    "readout": readout,
                    "left_out": g.loc[i, "unit_id"],
                    "acme_full": full["acme"],
                    "acme_loo": held["acme"],
                    "rho_tacstd2_loo": held["rho_tacstd2"],
                    "rho_partial_loo": held["rho_tacstd2_given_cldn4"],
                    "delta_acme": held["acme"] - full["acme"],
                })

    cohort.to_csv(OUT / "concordant4_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(meta_rows).to_csv(OUT / "concordant4_meta.tsv", sep="\t", index=False)
    pd.DataFrame(sens).to_csv(OUT / "concordant4_pooled_rows.tsv", sep="\t", index=False)
    pd.DataFrame(influence).to_csv(OUT / "concordant4_leave_one_patient.tsv", sep="\t", index=False)
    (OUT / "concordant4_check.json").write_text(json.dumps({
        "locked_cldn4_pct_rho": LOCKED_CLDN4_PCT_RHO,
        "recomputed": check,
        "n": int(len(patient)),
        "cohort_n": patient.groupby("dataset").size().to_dict(),
    }, indent=2))
    print(pd.DataFrame(meta_rows)[["dataset", "readout", "rho_tacstd2", "rho_cldn4",
                                   "rho_tacstd2_given_cldn4", "rho_cldn4_given_tacstd2",
                                   "attenuation_pct", "acme", "verdict"]].to_string(index=False))
    print("pipeline check CLDN4 %pos", check["rho"])


if __name__ == "__main__":
    main()
