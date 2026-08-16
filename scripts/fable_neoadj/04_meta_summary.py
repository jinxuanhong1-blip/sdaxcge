"""Cross-dataset harmonized summary for TACSTD2 / CLDN4 vs MPR.

Builds one master table and a forest plot of the marker->MPR ranking AUC
(AUC = P(marker higher in responder than non-responder)) with bootstrap 95%
CIs. Baseline (pre-treatment) and post-treatment cohorts are kept visually
distinct because they answer different questions and are NOT pooled into a
single p-value.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import GENES, TAB, FIG, mwu_stats, stouffer

# (dataset key, csv, MPR col, gene->value col template, timing, label)
SOURCES = [
    ("GSE207422_bulk", "GSE207422_bulk_per_sample.csv", "MPR", "{g}",
     "baseline", "GSE207422 bulk (baseline biopsy)"),
    ("GSE207422_scRNA_malig_post", "GSE207422_scrna_malignant_per_sample.csv", "MPR", "{g}_mean_log1p",
     "post-treatment", "GSE207422 scRNA malig (post-tx EPCAM+/PTPRC-)"),
    ("GSE241934_scRNA_epi", "GSE241934_scrna_per_patient.csv", "MPR", "{g}_cp10k_log1p",
     "post-treatment", "GSE241934 scRNA epi (resected tumor)"),
    ("GSE205335_malig_NSCLC", "GSE205335_malignant_per_patient.csv", "responder", "{g}_cp10k_log1p",
     "advanced-ICI-RECIST", "GSE205335 malig leftover (advanced ICI, RECIST, NSCLC)"),
]


def boot_auc_ci(pos, neg, n=2000, seed=0):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n):
        p = rng.choice(pos, len(pos), replace=True)
        q = rng.choice(neg, len(neg), replace=True)
        aucs.append((p[:, None] > q[None, :]).mean() + 0.5 * (p[:, None] == q[None, :]).mean())
    return (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


def main():
    rows = []
    for key, csv, mprcol, tmpl, timing, label in SOURCES:
        df = pd.read_csv(os.path.join(TAB, csv))
        if "pass_qc" in df.columns:
            df = df[df["pass_qc"]]
        if key == "GSE207422_scRNA_malig_post" and "timing" in df.columns:
            df = df[df["timing"] == "post"]
        if key == "GSE205335_malig_NSCLC" and "subtype" in df.columns:
            df = df[df["subtype"].isin(["ADC", "SQ"])]
        df = df[df[mprcol].notna()]
        for g in GENES:
            col = tmpl.format(g=g)
            pos = df.loc[df[mprcol] == 1, col].values
            neg = df.loc[df[mprcol] == 0, col].values
            st = mwu_stats(pos, neg)
            lo, hi = boot_auc_ci(pos, neg)
            rows.append({
                "dataset": key, "label": label, "timing": timing, "gene": g,
                "n_MPR": st["n_pos"], "n_nonMPR": st["n_neg"],
                "auc": st["auc"], "auc_lo": lo, "auc_hi": hi,
                "cliffs_delta": st["cliffs_delta"], "p": st["p"],
                "direction": "higher_in_MPR" if (st["auc"] or 0.5) > 0.5 else "lower_in_MPR",
            })
    master = pd.DataFrame(rows)
    master.to_csv(os.path.join(TAB, "MASTER_auc_summary.csv"), index=False)
    print(master.to_string())

    # Baseline-only signed combination (currently only GSE207422 baseline)
    combos = {}
    for g in GENES:
        sub = master[(master["timing"] == "baseline") & (master["gene"] == g)]
        signs = np.sign(sub["auc"].values - 0.5)
        weights = (sub["n_MPR"] + sub["n_nonMPR"]).values.astype(float)
        combos[g] = stouffer(sub["p"].values, signs, weights)
    with open(os.path.join(TAB, "MASTER_baseline_combined.json"), "w") as fh:
        json.dump(combos, fh, indent=2)
    print("baseline combined:", json.dumps(combos, indent=2))

    _forest(master)


def _forest(master):
    fig, axes = plt.subplots(1, len(GENES), figsize=(7.2 * len(GENES), 4.6), sharex=True)
    tcolor = {"baseline": "#2a6f97", "post-treatment": "#bc4749",
              "advanced-ICI-RECIST": "#6a994e"}
    for ax, g in zip(axes, GENES):
        sub = master[master["gene"] == g].reset_index(drop=True)
        y = np.arange(len(sub))[::-1]
        for yi, (_, r) in zip(y, sub.iterrows()):
            c = tcolor[r["timing"]]
            ax.plot([r["auc_lo"], r["auc_hi"]], [yi, yi], color=c, lw=2)
            ax.plot(r["auc"], yi, "o", color=c, ms=8)
        ax.axvline(0.5, color="grey", ls="--", lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['label']}\n(MPR {r['n_MPR']} vs {r['n_nonMPR']})"
                            for _, r in sub.iterrows()], fontsize=8)
        ax.set_xlim(0, 1)
        ax.set_xlabel("AUC: marker higher in MPR  (0.5 = no assoc.)")
        ax.set_title(g)
    fig.suptitle("TACSTD2 / CLDN4 association with MPR across neoadjuvant lung PD-1 datasets", y=1.05)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "MASTER_forest_auc.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
