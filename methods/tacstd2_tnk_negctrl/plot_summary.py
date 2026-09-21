#!/usr/bin/env python3
"""Redraw the partial-correlation figure from the saved tables. No resampling."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
ORDER = ["none", "CLDN4", "CLDN7", "KRT19", "EPCAM", "MUC1", "CLDN3"]
LABELS = ["unadjusted", "CLDN4", "CLDN7", "KRT19", "EPCAM", "MUC1", "CLDN3"]


def _forest(ax, rows, title, xlabel):
    by = {r["conditioner"]: r for r in rows}
    ypos = list(range(len(ORDER)))[::-1]
    rhos = [by[g]["rho_partial"] for g in ORDER]
    lo = [by[g]["rho_partial"] - by[g]["ci_lo"] for g in ORDER]
    hi = [by[g]["ci_hi"] - by[g]["rho_partial"] for g in ORDER]
    colors = []
    for g in ORDER:
        if g == "none":
            colors.append("#4d4d4d")
        elif g == "CLDN4":
            colors.append("#b2182b")
        elif g == "CLDN7":
            colors.append("#e08214")
        else:
            colors.append("#2166ac")
    ax.errorbar(rhos, ypos, xerr=[lo, hi], fmt="none", ecolor="#666666", elinewidth=1, capsize=2)
    ax.scatter(rhos, ypos, c=colors, s=38, zorder=3)
    ax.axvline(0, color="#888888", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels(LABELS)
    ax.set_xlabel(xlabel)
    ax.set_title(title)


def _shift(ax, meta_rows, cohort_rows, title):
    genes = [g for g in ORDER if g != "none"]
    by = {r["conditioner"]: r for r in meta_rows}
    ypos = list(range(len(genes)))[::-1]
    ax.axvline(0, color="#888888", lw=0.6)
    for i, g in enumerate(genes):
        color = "#b2182b" if g == "CLDN4" else ("#e08214" if g == "CLDN7" else "#2166ac")
        ax.scatter([by[g]["attenuation"]], [ypos[i]], c=color, s=38, zorder=3)
        xs = [r["attenuation"] for r in cohort_rows if r["conditioner"] == g]
        ax.scatter(xs, np.full(len(xs), ypos[i]), s=12, c="#999999", zorder=2)
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes)
    ax.set_xlabel("Change in ρ  (partial − unadjusted)")
    ax.set_title(title)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    c4 = pd.read_csv(TAB / "c4_meta_attenuation.tsv", sep="\t")
    c4c = pd.read_csv(TAB / "c4_cohort_attenuation.tsv", sep="\t")
    tcga = pd.read_csv(TAB / "tcga_meta_attenuation.tsv", sep="\t")
    tcgac = pd.read_csv(TAB / "tcga_cohort_attenuation.tsv", sep="\t")
    c4_rows = c4[(c4.subset == "all4") & (c4.score == "pct")].to_dict("records")
    tcga_rows = tcga[(tcga.block == "funnel7") & (tcga.outcome == "tnk")].to_dict("records")
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6))
    _forest(
        axes[0],
        c4_rows,
        "Concordant-4  (n = 65, % positive)",
        "DL Spearman ρ   TACSTD2 vs T/NK fraction",
    )
    _forest(
        axes[1],
        tcga_rows,
        "TCGA funnel-7  (n = 3,444)",
        "DL Spearman ρ   TACSTD2 vs T/NK score",
    )
    fig.tight_layout()
    fig.savefig(FIG / "partial_rho_forest.png", dpi=160)
    fig.savefig(FIG / "partial_rho_forest.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    _forest(axes[0], c4_rows, "Concordant-4 partial ρ", "DL Spearman ρ  TACSTD2 %pos vs T/NK")
    c4_cohort = c4c[(c4c.subset == "all4") & (c4c.score == "pct") & (c4c.conditioner != "none")].to_dict("records")
    _shift(axes[1], c4_rows, c4_cohort, "Concordant-4 shift")
    fig.tight_layout()
    fig.savefig(FIG / "c4_attenuation.png", dpi=160)
    fig.savefig(FIG / "c4_attenuation.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    _forest(axes[0], tcga_rows, "TCGA funnel-7 partial ρ", "DL Spearman ρ  TACSTD2 vs T/NK score")
    tcga_cohort = tcgac[(tcgac.block == "funnel7") & (tcgac.outcome == "tnk") & (tcgac.conditioner != "none")].to_dict("records")
    _shift(axes[1], tcga_rows, tcga_cohort, "TCGA funnel-7 shift")
    fig.tight_layout()
    fig.savefig(FIG / "tcga_attenuation.png", dpi=160)
    fig.savefig(FIG / "tcga_attenuation.pdf")
    plt.close(fig)
    print("wrote figures")


if __name__ == "__main__":
    main()
