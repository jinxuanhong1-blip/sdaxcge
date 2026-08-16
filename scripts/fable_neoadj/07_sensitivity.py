"""Sensitivity / confounder checks on the primary GSE207422 baseline bulk set.

TACSTD2 is a known squamous-enriched gene and a bulk epithelial marker.
If MPR rates or tumor content differ by histology, the unadjusted
TACSTD2-vs-MPR AUC can be a histology or purity artifact. This script
does not change the primary null; it documents those alternatives.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from common import GENES, TAB, FIG, mwu_stats, spearman

BULK = os.path.join(TAB, "GSE207422_bulk_per_sample.csv")
G241 = os.path.join(TAB, "GSE241934_scrna_per_patient.csv")


def hist_group(p):
    p = str(p).strip().lower()
    if p.startswith("adeno"):
        return "Adeno"
    if p.startswith("squam"):
        return "Squamous"
    return "Other"


def residualize(y, x):
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = ~(np.isnan(y) | np.isnan(x))
    out = np.full_like(y, np.nan, dtype=float)
    if m.sum() < 3:
        return out
    slope, intercept, *_ = stats.linregress(x[m], y[m])
    out[m] = y[m] - (intercept + slope * x[m])
    return out


def main():
    df = pd.read_csv(BULK)
    df["hist"] = df["Pathology"].map(hist_group)
    rows = []

    # 1) gene vs histology (Squamous vs Adeno)
    for g in GENES + ["EPCAM"]:
        pos = df.loc[df["hist"] == "Squamous", g].values
        neg = df.loc[df["hist"] == "Adeno", g].values
        st = mwu_stats(pos, neg)
        row = {"dataset": "GSE207422_bulk", "gene": g, "comparison": "Squamous_vs_Adeno"}
        row.update(st)
        rows.append(row)

    # 2) MPR vs non-MPR within histology
    for h in ["Adeno", "Squamous"]:
        sub = df[df["hist"] == h]
        for g in GENES:
            pos = sub.loc[sub["MPR"] == 1, g].values
            neg = sub.loc[sub["MPR"] == 0, g].values
            st = mwu_stats(pos, neg)
            row = {"dataset": "GSE207422_bulk", "gene": g,
                   "comparison": f"MPR_vs_nonMPR_within_{h}"}
            row.update(st)
            rows.append(row)

    # 3) pCR vs non-pCR (all histologies)
    for g in GENES:
        pos = df.loc[df["pCR"] == 1, g].values
        neg = df.loc[df["pCR"] == 0, g].values
        st = mwu_stats(pos, neg)
        row = {"dataset": "GSE207422_bulk", "gene": g, "comparison": "pCR_vs_nonpCR"}
        row.update(st)
        rows.append(row)

    # 4) TACSTD2 vs EPCAM (tumor-content proxy) and EPCAM-residualized MPR test
    for g in GENES:
        sp = spearman(df[g].values, df["EPCAM"].values)
        rows.append({"dataset": "GSE207422_bulk", "gene": g,
                     "comparison": "Spearman_vs_EPCAM",
                     "n_pos": sp["n"], "n_neg": np.nan,
                     "median_pos": np.nan, "median_neg": np.nan,
                     "U": np.nan, "p": sp["p"], "auc": np.nan,
                     "cliffs_delta": sp["rho"]})
        resid = residualize(df[g].values, df["EPCAM"].values)
        pos = resid[df["MPR"] == 1]
        neg = resid[df["MPR"] == 0]
        st = mwu_stats(pos, neg)
        row = {"dataset": "GSE207422_bulk", "gene": g,
               "comparison": "MPR_vs_nonMPR_EPCAM_residualized"}
        row.update(st)
        rows.append(row)

    # 5) MPR rate by histology (Fisher)
    tab = pd.crosstab(df["hist"], df["MPR"])
    # Adeno vs Squamous only
    sub = df[df["hist"].isin(["Adeno", "Squamous"])]
    ct = pd.crosstab(sub["hist"], sub["MPR"])
    if ct.shape == (2, 2):
        oddsratio, p_fish = stats.fisher_exact(ct.values)
    else:
        oddsratio, p_fish = np.nan, np.nan
    rows.append({"dataset": "GSE207422_bulk", "gene": "MPR_rate",
                 "comparison": "Fisher_Adeno_vs_Squamous_MPR",
                 "n_pos": int((sub["hist"] == "Adeno").sum()),
                 "n_neg": int((sub["hist"] == "Squamous").sum()),
                 "median_pos": float(sub.loc[sub["hist"] == "Adeno", "MPR"].mean()),
                 "median_neg": float(sub.loc[sub["hist"] == "Squamous", "MPR"].mean()),
                 "U": float(oddsratio) if oddsratio == oddsratio else np.nan,
                 "p": float(p_fish) if p_fish == p_fish else np.nan,
                 "auc": np.nan, "cliffs_delta": np.nan})

    # 6) GSE241934 histology if present (cohort column + we can join from meta? per-patient csv has no hist)
    # skip if missing

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(TAB, "GSE207422_bulk_sensitivity.csv"), index=False)
    print(res.to_string())

    _plots(df)

    summary = {
        "adeno_n": int((df["hist"] == "Adeno").sum()),
        "squamous_n": int((df["hist"] == "Squamous").sum()),
        "other_n": int((df["hist"] == "Other").sum()),
        "adeno_mpr_rate": float(df.loc[df["hist"] == "Adeno", "MPR"].mean()),
        "squamous_mpr_rate": float(df.loc[df["hist"] == "Squamous", "MPR"].mean()),
        "fisher_p_hist_mpr": float(p_fish) if p_fish == p_fish else None,
        "note": "TACSTD2 tracks EPCAM (purity) and squamous histology; residualized MPR test remains null",
    }
    with open(os.path.join(TAB, "GSE207422_bulk_sensitivity_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _plots(df):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))
    # TACSTD2 by histology, colored by MPR
    ax = axes[0]
    order = ["Adeno", "Squamous", "Other"]
    for i, h in enumerate(order, 1):
        sub = df[df["hist"] == h]
        ax.boxplot([sub["TACSTD2"].values], positions=[i], widths=0.5, showfliers=False)
        rng = np.random.default_rng(5)
        colors = np.where(sub["MPR"] == 1, "#d1495b", "#8aa1b1")
        ax.scatter(rng.normal(i, 0.08, len(sub)), sub["TACSTD2"], c=colors, s=28, zorder=3)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(order)
    ax.set_ylabel("TACSTD2 log2 TPM")
    ax.set_title("GSE207422 bulk: TACSTD2 by histology\n(red=MPR, grey=non-MPR)")

    ax = axes[1]
    ax.scatter(df["EPCAM"], df["TACSTD2"],
               c=np.where(df["MPR"] == 1, "#d1495b", "#8aa1b1"), s=36)
    sp = spearman(df["TACSTD2"].values, df["EPCAM"].values)
    ax.set_xlabel("EPCAM log2 TPM (tumor-content proxy)")
    ax.set_ylabel("TACSTD2 log2 TPM")
    ax.set_title(f"TACSTD2 vs EPCAM  rho={sp['rho']:.2f} p={sp['p']:.3g}")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE207422_bulk_sensitivity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
