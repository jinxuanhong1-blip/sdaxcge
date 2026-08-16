"""Aggregate all cohorts into one master table + effect-size forest plot.

Provides the multi-dataset "pre / on / post" synthesis and applies
Benjamini-Hochberg FDR across the primary two-group comparisons.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import config as C
from common import savefig, PALETTE
import matplotlib.pyplot as plt


def main() -> None:
    sc = pd.read_csv(C.TABLES_DIR / "gse207422_scrna_stats.csv")
    bulk = pd.read_csv(C.TABLES_DIR / "baseline_bulk_stats.csv")
    paired = pd.read_csv(C.TABLES_DIR / "gse91061_paired_stats.csv")
    ar = pd.read_csv(C.TABLES_DIR / "gse248249_paired_stats.csv")
    mouse = pd.read_csv(C.TABLES_DIR / "gse246922_mouse_stats.csv")

    # unify the two-group comparisons (those with an AUC / log2fc)
    keep = ["dataset", "gene", "comparison", "n_hi", "n_lo",
            "log2fc_hi_vs_lo", "auc", "test", "statistic", "p_value"]
    two_group = pd.concat([
        sc[keep], bulk[keep],
        paired[paired["auc"].notna()][keep],
        ar[ar["auc"].notna()][keep] if "auc" in ar.columns else ar.iloc[0:0],
        mouse[keep],
    ], ignore_index=True)

    # tag a "timepoint axis" for interpretation
    def axis(row):
        c = row["comparison"]
        if "baseline" in c:
            return "PRE (baseline vs response)"
        if "on-pre delta" in c:
            return "ON (pre->on change vs response)"
        if "acquired resistance" in c:
            return "POST (acquired resistance)"
        if "ICB_resistant" in c or "IFNg" in c:
            return "MOUSE (lung ICB models)"
        if "post:" in c or "post_vs_pre" in c:
            return "POST (post-treatment)"
        return "other"
    two_group["timepoint_axis"] = two_group.apply(axis, axis=1)

    # BH-FDR across all response-contrast comparisons
    mask = two_group["p_value"].notna()
    two_group["p_fdr_bh"] = np.nan
    if mask.sum() > 1:
        two_group.loc[mask, "p_fdr_bh"] = multipletests(
            two_group.loc[mask, "p_value"], method="fdr_bh")[1]

    two_group.to_csv(C.TABLES_DIR / "master_summary.csv", index=False)
    print(two_group.to_string())

    # ---------- forest plot of effect sizes (AUC of hi vs lo) ----------
    resp_contrasts = two_group[two_group["comparison"].str.contains(
        "responder|DCB", regex=True)].copy()
    resp_contrasts = resp_contrasts.sort_values(["gene", "dataset"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True)
    for ax, gene in zip(axes, C.TARGET_GENES):
        d = resp_contrasts[resp_contrasts.gene == gene].reset_index(drop=True)
        y = np.arange(len(d))
        colors = [PALETTE["responder"] if a >= 0.5 else PALETTE["non_responder"] for a in d["auc"]]
        ax.scatter(d["auc"], y, s=90, color=colors, zorder=3, edgecolor="black", linewidth=0.6)
        for yi, (_, r) in zip(y, d.iterrows()):
            star = "*" if (r["p_value"] == r["p_value"] and r["p_value"] < 0.05) else ""
            ax.text(r["auc"], yi + 0.18, f"p={r['p_value']:.3f}{star}", fontsize=8, ha="center")
        ax.axvline(0.5, color="grey", ls="--", lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.dataset}\n{r.comparison.split(':')[0]}" for _, r in d.iterrows()],
                           fontsize=8)
        ax.set_xlim(0, 1)
        ax.set_xlabel("AUC  (>0.5: higher in responders)")
        ax.set_title(gene)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Effect sizes: TACSTD2 / CLDN4 vs ICI response across cohorts\n"
                 "(AUC=probability responder > non-responder; dashed line = no effect)",
                 fontsize=11)
    savefig(fig, C.FIGURES_DIR / "fig4_effectsize_forest.png")


if __name__ == "__main__":
    main()
