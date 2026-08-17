#!/usr/bin/env python3
"""Figures for CLDN4-high vs low T/NK and B composition."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "results"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PRIMARY = [
    "GSE207422_post_A3mal",
    "GSE241934_IIT",
    "GSE241934_RWC",
    "GSE291670",
    "GSE205335",
    "GSE253013_tumor",
]


def fig_fractions() -> None:
    table = pd.read_csv(OUT / "composition_table.tsv", sep="\t")
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4), sharey=False)
    axes = axes.ravel()
    for ax, sl in zip(axes, PRIMARY):
        sub = table[(table["slice"] == sl) & table["eligible"].astype(bool)].copy()
        for i, (comp, col) in enumerate([("T/NK", "frac_TNK"), ("B", "frac_B")]):
            if sub[col].isna().all():
                col = "frac_B_plasma"
                if comp == "B":
                    comp = "B+plasma"
            y0 = sub.loc[sub["cldn4_high"] == 0, col].to_numpy(dtype=float)
            y1 = sub.loc[sub["cldn4_high"] == 1, col].to_numpy(dtype=float)
            x0 = np.full(len(y0), i * 2 + 0.7)
            x1 = np.full(len(y1), i * 2 + 1.3)
            ax.scatter(x0, y0, s=18, c="#4C78A8", alpha=0.85, zorder=3)
            ax.scatter(x1, y1, s=18, c="#E45756", alpha=0.85, zorder=3)
            if len(y0):
                ax.hlines(np.median(y0), i * 2 + 0.45, i * 2 + 0.95, color="#4C78A8", lw=2)
            if len(y1):
                ax.hlines(np.median(y1), i * 2 + 1.05, i * 2 + 1.55, color="#E45756", lw=2)
        ax.set_xticks([1, 3])
        ax.set_xticklabels(["T/NK\nlow|high", "B\nlow|high"], fontsize=8)
        ax.set_title(f"{sl}\nn={int(sub.shape[0])}", fontsize=9)
        ax.set_ylabel("fraction" if ax in (axes[0], axes[3]) else "")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.legend(
        [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#4C78A8", markersize=7),
         plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#E45756", markersize=7)],
        ["CLDN4-low", "CLDN4-high"],
        loc="lower center",
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Patient T/NK and B fractions vs malignant CLDN4 median split", fontsize=11)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig(FIG / "fig1_fractions_high_vs_low.png", dpi=160)
    fig.savefig(FIG / "fig1_fractions_high_vs_low.pdf")
    plt.close(fig)


def fig_honest_n() -> None:
    hon = pd.read_csv(OUT / "honest_n.tsv", sep="\t")
    hon = hon[hon["slice"].isin(PRIMARY + ["GSE207422_drmref"])].copy()
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    y = np.arange(len(hon))
    ax.barh(y, hon["n_extracted"], color="#D8D8D8", label="extracted")
    ax.barh(y, hon["n_cldn4_eligible"], color="#4C78A8", label="CLDN4-eligible")
    ax.set_yticks(y)
    ax.set_yticklabels(hon["slice"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("patients (not cells)")
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Honest n: CLDN4-eligible patients per slice")
    for i, row in enumerate(hon.itertuples()):
        ax.text(row.n_cldn4_eligible + 0.3, i, f"{row.n_high} vs {row.n_low}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_honest_n.png", dpi=160)
    fig.savefig(FIG / "fig2_honest_n.pdf")
    plt.close(fig)


def fig_alr_forest() -> None:
    t = pd.read_csv(OUT / "tests_full.tsv", sep="\t")
    t = t[t["slice"].isin(PRIMARY + ["ICI_pool_cohort_centered"])].copy()
    t = t.sort_values(["compartment", "slice"])
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    y = np.arange(len(t))
    colors = ["#4C78A8" if c == "TNK" else "#F58518" for c in t["compartment"]]
    ax.axvline(0, color="#888", lw=0.8)
    ax.errorbar(t["alr_effect"], y, xerr=t["alr_se"].fillna(0), fmt="none", ecolor="#999", elinewidth=1)
    ax.scatter(t["alr_effect"], y, c=colors, s=28, zorder=3)
    labels = [f"{s} · {c}  n={n}" for s, c, n in zip(t["slice"], t["compartment"], t["n"])]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("ALR slope (CLDN4-high vs low)")
    ax.set_title("ALR effects; primary p is patient permutation (not OLS)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_alr_forest.png", dpi=160)
    fig.savefig(FIG / "fig3_alr_forest.pdf")
    plt.close(fig)


def main() -> None:
    fig_fractions()
    fig_honest_n()
    fig_alr_forest()
    print(f"wrote figures in {FIG}")


if __name__ == "__main__":
    main()
