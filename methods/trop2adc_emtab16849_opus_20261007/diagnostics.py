#!/usr/bin/env python3
"""Why the E-MTAB-16849 contrast is null: a failing pre-dose negative control
and panel deltas that co-move with arm-level offsets.

0h tumours were harvested before dosing, so any SG-vs-control difference at 0h
is a floor on arm-to-arm technical/compositional noise.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

pp = pd.read_csv(f"{RES}/emtab16849_score_delta_per_pair.csv")
pc = pd.read_csv(f"{RES}/emtab16849_percell_by_timepoint.csv")

rows = []
for panel in ["core_APM", "core_IFN", "STING_axis"]:
    x, y = pp["delta_sanity"], pp[f"delta_{panel}"]
    pr, pp_p = stats.pearsonr(x, y)
    sr, sp_p = stats.spearmanr(x, y)
    rows.append({"panel": panel, "n_pairs": len(x), "pearson_r_vs_sanity": pr,
                 "pearson_p": pp_p, "spearman_rho_vs_sanity": sr, "spearman_p": sp_p})
co = pd.DataFrame(rows)
co.to_csv(f"{RES}/emtab16849_offset_correlation.csv", index=False)

pre = pc[pc["timepoint"] == "0h"][["panel", "delta", "cohens_d"]].set_index("panel")
post = (pc[pc["timepoint"] != "0h"].groupby("panel")["delta"]
        .agg(["mean", "min", "max"]))
nc = pre.join(post, how="left").rename(columns={
    "delta": "delta_0h_PREDOSE", "cohens_d": "cohens_d_0h_PREDOSE",
    "mean": "mean_delta_postdose", "min": "min_delta_postdose",
    "max": "max_delta_postdose"})
nc["abs_0h_exceeds_mean_postdose"] = (nc["delta_0h_PREDOSE"].abs()
                                      > nc["mean_delta_postdose"].abs())
nc.to_csv(f"{RES}/emtab16849_predose_negative_control.csv")

print("Panel-score deltas co-move with the CLDN4/TACSTD2 'sanity' offset:")
print(co.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
print("\n0h PRE-DOSE negative control (cell level) vs post-dose timepoints:")
print(nc.to_string(float_format=lambda v: f"{v: .4f}"))
print("\nIf the 0h row is as large as the post-dose rows, the apparent effects")
print("are arm-level offsets, not SG pharmacology.")
