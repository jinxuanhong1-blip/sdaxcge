#!/usr/bin/env python3
"""Figures for the TISMO lung-only Tacstd2/Cldn4 hunt."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "hunt_tismo_lung"
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def waterfall() -> None:
    g = pd.read_csv(RES / "icb_paired_group_stats.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=False)
    for ax, gene in zip(axes, ["Tacstd2", "Cldn4"]):
        t = g[g["gene"] == gene].sort_values("delta_mean")
        colors = []
        for _, r in t.iterrows():
            if r["is_lung"]:
                colors.append("#d62728")
            elif r["delta_mean"] > 0:
                colors.append("#4c78a8")
            else:
                colors.append("#9ecae1")
        ax.bar(np.arange(len(t)), t["delta_mean"], color=colors, width=1.0)
        ax.axhline(0, color="k", lw=0.6)
        n_up = int((t["direction_mean"] == "up").sum())
        ax.set_title(f"{gene}  {n_up}/{len(t)} groups up  (lung in red)")
        ax.set_xlabel("ICB comparison groups (sorted Δ)")
        ax.set_ylabel("mean(ICB) − mean(baseline)\nlog2(TPM+1)")
        ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(FIG / "all_cancer_delta_waterfall.png", dpi=160)
    plt.close(fig)


def lung_strips() -> None:
    s = pd.read_csv(RES / "lung_icb_per_sample.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    order = [
        "LLC_GSE155972_antiCTLA4&antiPD1",
        "LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1",
    ]
    labels = ["LLC WT", "LLC Setdb1-KO"]
    rng = np.random.default_rng(0)
    for ax, gene in zip(axes, ["Tacstd2", "Cldn4"]):
        sub = s[s["geneID"] == gene]
        xs, ys, cs = [], [], []
        for i, grp in enumerate(order):
            g = sub[sub["group"] == grp]
            for base, x0, col in [(1, i - 0.18, "#4c78a8"), (0, i + 0.18, "#d62728")]:
                vals = g.loc[g["Baseline"] == base, "value"].to_numpy()
                jitter = rng.normal(0, 0.04, size=len(vals))
                xs.extend(x0 + jitter)
                ys.extend(vals)
                cs.extend([col] * len(vals))
                ax.hlines(np.mean(vals), x0 - 0.12, x0 + 0.12, color=col, lw=2)
        ax.scatter(xs, ys, c=cs, s=28, alpha=0.85, edgecolors="k", linewidths=0.3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(labels)
        ax.set_title(gene + "  GSE155972")
        ax.set_ylabel("TISMO log2(TPM+1)")
        ax.scatter([], [], c="#4c78a8", label="baseline")
        ax.scatter([], [], c="#d62728", label="ICB")
        ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "lung_gse155972_strips.png", dpi=160)
    plt.close(fig)


def sign_bars() -> None:
    sm = pd.read_csv(RES / "icb_paired_sign_summary.csv")
    want = [
        "Tacstd2_all_cancer",
        "Cldn4_all_cancer",
        "Actb_all_cancer",
        "Cd8a_all_cancer",
        "Tacstd2_lung_only",
        "Cldn4_lung_only",
    ]
    sm = sm[sm["label"].isin(want)].set_index("label").loc[want]
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    x = np.arange(len(sm))
    ax.bar(x, sm["n_up_mean"], color="#d62728", label="up")
    ax.bar(x, sm["n_down_mean"], bottom=sm["n_up_mean"], color="#4c78a8", label="down")
    ax.bar(
        x,
        sm["n_tie_mean"],
        bottom=sm["n_up_mean"] + sm["n_down_mean"],
        color="#bbbbbb",
        label="tie",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace("_", "\n") for l in sm.index], fontsize=8)
    ax.set_ylabel("ICB comparison groups")
    ax.set_title("Direction of mean ICB vs baseline (TISMO official groups)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "sign_counts.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    waterfall()
    lung_strips()
    sign_bars()
    print("figures in", FIG)
