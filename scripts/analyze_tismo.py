#!/usr/bin/env python3
"""Tacstd2 in TISMO in vivo syngeneic mouse tumours: ICB arm vs matched control arm.

Important scope note carried into the outputs: TISMO in vivo ICB samples are *not*
within-animal longitudinal biopsies. Every ICB-treated sample in TISMO has Baseline=0 and
comes from a different mouse than the untreated/isotype tumours it is compared with, so a
"contrast" here is a treated arm versus a control arm of the same cell line in the same
study. It approximates "after ICB vs before/without ICB" at the model level only.
"""
import json
import os

import pandas as pd
import pyreadr
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "tismo")
OUT = os.path.join(ROOT, "results", "hunt_paired_up", "tismo")

GENE = "Tacstd2"
CONTROL_TOKENS = ("no_treatment", "isotype", "vehicle", "untreated", "control", "pbs", "igg")
# Expression is log2-scale; below this both arms are effectively at the floor and the sign
# of the difference is noise rather than biology.
EXPRESSED_MIN = 0.5


def is_control(treatment):
    s = str(treatment).lower()
    return any(tok in s for tok in CONTROL_TOKENS)


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = pd.read_csv(os.path.join(DATA, "TISMO_vivosample_annotations.csv"))
    expr = pyreadr.read_r(os.path.join(DATA, "TISMO_expressionvivo_profiles.RDS"))[None]

    audit = {
        "gene": GENE,
        "meta_samples": int(len(meta)),
        "expression_matrix_shape": list(expr.shape),
        "gene_present_in_matrix": bool(GENE in expr.index),
    }

    meta["has_expr"] = meta.SRX_ID.isin(expr.columns)
    audit["samples_with_expression"] = int(meta.has_expr.sum())
    audit["samples_missing_expression"] = int((~meta.has_expr).sum())
    meta = meta[meta.has_expr].copy()

    meta[GENE] = expr.loc[GENE, meta.SRX_ID].to_numpy()
    meta["is_control_arm"] = (meta.ICB == 0) & meta.Mouse_treatment.map(is_control)
    meta["is_icb_arm"] = meta.ICB == 1

    audit["icb_samples"] = int(meta.is_icb_arm.sum())
    audit["control_samples"] = int(meta.is_control_arm.sum())
    audit["non_icb_non_control_samples_excluded"] = int(
        ((meta.ICB == 0) & ~meta.is_control_arm).sum()
    )

    meta.sort_values(["Study_ID", "Cell_Line", "Mouse_treatment"]).to_csv(
        os.path.join(OUT, "tismo_invivo_samples_tacstd2.csv"), index=False
    )

    controls = meta[meta.is_control_arm]
    treated = meta[meta.is_icb_arm]

    rows = []
    arm_keys = ["Study_ID", "Cell_Line", "Mouse_treatment", "Timepoint"]
    for key, arm in treated.groupby(arm_keys, dropna=False):
        study, line, treatment, timepoint = key
        pool = controls[(controls.Study_ID == study) & (controls.Cell_Line == line)]
        match_level = "study+cell_line"
        if len(pool) and pd.notna(timepoint):
            same_tp = pool[pool.Timepoint == timepoint]
            if len(same_tp) >= 2:
                pool = same_tp
                match_level = "study+cell_line+timepoint"
        if pool.empty:
            rows.append(
                {
                    "study": study,
                    "cell_line": line,
                    "cancer_type": arm.Cancer_type.iloc[0],
                    "icb_arm": treatment,
                    "timepoint": timepoint,
                    "n_icb": len(arm),
                    "n_control": 0,
                    "status": "no_matched_control_arm",
                }
            )
            continue

        t_vals, c_vals = arm[GENE].to_numpy(), pool[GENE].to_numpy()
        t_mean, c_mean = float(t_vals.mean()), float(c_vals.mean())
        p = (
            float(stats.mannwhitneyu(t_vals, c_vals, alternative="two-sided").pvalue)
            if len(t_vals) >= 2 and len(c_vals) >= 2
            else float("nan")
        )
        rows.append(
            {
                "study": study,
                "cell_line": line,
                "cancer_type": arm.Cancer_type.iloc[0],
                "icb_arm": treatment,
                "timepoint": timepoint,
                "n_icb": len(arm),
                "n_control": len(pool),
                "control_match_level": match_level,
                "mean_log2_icb": round(t_mean, 4),
                "mean_log2_control": round(c_mean, 4),
                "delta_log2": round(t_mean - c_mean, 4),
                "mannwhitney_p": p,
                "expressed": bool(max(t_mean, c_mean) >= EXPRESSED_MIN),
                "direction": "up" if t_mean > c_mean else ("down" if t_mean < c_mean else "flat"),
                "status": "ok",
            }
        )

    res = pd.DataFrame(rows).sort_values(["study", "cell_line", "icb_arm"])
    res.to_csv(os.path.join(OUT, "tismo_icb_vs_control_contrasts.csv"), index=False)

    ok = res[res.status == "ok"]
    powered = ok[(ok.n_icb >= 2) & (ok.n_control >= 2)]
    expressed = powered[powered.expressed]

    def tally(df, label):
        up = int((df.delta_log2 > 0).sum())
        sig_up = int(((df.delta_log2 > 0) & (df.mannwhitney_p < 0.05)).sum())
        sig_dn = int(((df.delta_log2 < 0) & (df.mannwhitney_p < 0.05)).sum())
        binom = (
            float(stats.binomtest(up, len(df), 0.5).pvalue) if len(df) else float("nan")
        )
        return {
            "label": label,
            "n_contrasts": int(len(df)),
            "n_up": up,
            "n_down": int((df.delta_log2 < 0).sum()),
            "frac_up": round(up / len(df), 3) if len(df) else None,
            "n_up_p<0.05": sig_up,
            "n_down_p<0.05": sig_dn,
            "median_delta_log2": round(float(df.delta_log2.median()), 4) if len(df) else None,
            "sign_test_p": binom,
        }

    audit["contrasts_attempted"] = int(len(res))
    audit["contrasts_without_control"] = int((res.status != "ok").sum())
    audit["arm_size_min_icb"] = int(ok.n_icb.min())
    audit["arm_size_min_control"] = int(ok.n_control.min())
    audit["tallies"] = [
        tally(ok, "all contrasts with a matched control arm"),
        tally(powered, "contrasts with n>=2 per arm"),
        tally(expressed, f"n>=2 per arm AND {GENE} above floor (max arm mean >= {EXPRESSED_MIN} log2)"),
        tally(
            expressed[expressed.cancer_type.str.contains("Lung", na=False)],
            "lung carcinoma subset (expressed, n>=2 per arm)",
        ),
    ]

    # Sensitivity: the number of contrasts depends entirely on how an "arm" is defined, so
    # report the two other groupings a reader might reasonably pick.
    sens = []
    for label, keys in [
        ("arm = study + cell line + treatment (timepoints pooled)",
         ["Study_ID", "Cell_Line", "Mouse_treatment"]),
        ("arm = study + cell line + ICB drug class (treatments pooled)",
         ["Study_ID", "Cell_Line", "ICB_group"]),
        ("arm = study + cell line (all ICB pooled)", ["Study_ID", "Cell_Line"]),
    ]:
        deltas = []
        for key, arm in treated.groupby(keys, dropna=False):
            key = key if isinstance(key, tuple) else (key,)
            study, line = key[0], key[1]
            pool = controls[(controls.Study_ID == study) & (controls.Cell_Line == line)]
            if pool.empty:
                continue
            deltas.append(float(arm[GENE].mean() - pool[GENE].mean()))
        up = sum(d > 0 for d in deltas)
        sens.append(
            {
                "definition": label,
                "n_contrasts": len(deltas),
                "n_up": up,
                "frac_up": round(up / len(deltas), 3) if deltas else None,
                "sign_test_p": float(stats.binomtest(up, len(deltas), 0.5).pvalue) if deltas else None,
            }
        )
    audit["sensitivity_to_contrast_definition"] = sens

    with open(os.path.join(OUT, "tismo_summary.json"), "w") as fh:
        json.dump(audit, fh, indent=2)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
