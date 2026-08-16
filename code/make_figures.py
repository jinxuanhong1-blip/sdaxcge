#!/usr/bin/env python3
"""Figures for Claim A1: forest plot of purity-adjusted partial Spearman and
pooled rank scatter plots."""
import os
from math import sqrt, exp, log

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "claim_A1")

SIGS = ["immune_tcell_effector", "cytotoxic", "exhaustion"]
SIG_LABEL = {"immune_tcell_effector": "Immune (T-cell effector)",
             "cytotoxic": "Cytotoxic", "exhaustion": "Exhaustion / checkpoint"}


def fisher_ci(r, n):
    if n is None or n <= 3 or r is None or np.isnan(r):
        return np.nan, np.nan
    z = 0.5 * log((1 + r) / (1 - r))
    se = 1.0 / sqrt(n - 3)
    lo, hi = z - 1.96 * se, z + 1.96 * se
    b = lambda x: (exp(2 * x) - 1) / (exp(2 * x) + 1)
    return b(lo), b(hi)


def forest():
    per = pd.read_csv(os.path.join(OUT, "per_cohort_correlations.csv"))
    pool = pd.read_csv(os.path.join(OUT, "pooled_correlations.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    for ax, sig in zip(axes, SIGS):
        entries = []
        for _, r in per[per.signature == sig].iterrows():
            entries.append((r["cohort"], r["partial_rho_purityadj"], r["n_partial"], "cohort"))
        for _, r in pool[pool.signature == sig].iterrows():
            entries.append((r["pool"], r["partial_rho_purityadj"], r["n_partial"], "pool"))
        entries = entries[::-1]
        ys = range(len(entries))
        for y, (name, rho, n, kind) in zip(ys, entries):
            lo, hi = fisher_ci(rho, n)
            color = "#b2182b" if kind == "pool" else "#2166ac"
            ax.errorbar(rho, y, xerr=[[rho - lo], [hi - rho]], fmt="o",
                        color=color, capsize=3, ms=6 if kind == "pool" else 5)
            ax.text(0.02, y, f"n={int(n)}", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=7, color="#555")
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set_yticks(list(ys))
        ax.set_yticklabels([e[0] for e in entries], fontsize=8)
        ax.set_title(SIG_LABEL[sig], fontsize=10)
        ax.set_xlabel("Partial Spearman rho (purity-adjusted)")
    fig.suptitle("TACSTD2 vs immune signatures, purity-adjusted partial Spearman "
                 "(negative = replicates claim A1)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "forest_partial_spearman.png"), dpi=150)
    plt.close(fig)


def scatter():
    frames = []
    for cid in ["TCGA_LUAD", "TCGA_LUSC", "OncoSG_LUAD"]:
        df = pd.read_csv(os.path.join(OUT, f"sample_data_{cid}.csv"))
        df["cohort"] = cid
        frames.append(df)
    d = pd.concat(frames, ignore_index=True)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, sig in zip(axes, SIGS):
        sub = d.dropna(subset=["TACSTD2", sig])
        xr = sub["TACSTD2"].rank()
        yr = sub[sig].rank()
        ax.scatter(xr, yr, s=6, alpha=0.35, color="#2166ac", edgecolors="none")
        z = np.polyfit(xr, yr, 1)
        xs = np.linspace(xr.min(), xr.max(), 50)
        ax.plot(xs, np.polyval(z, xs), color="#b2182b", lw=2)
        ax.set_title(SIG_LABEL[sig], fontsize=10)
        ax.set_xlabel("TACSTD2 rank")
        ax.set_ylabel(f"{SIG_LABEL[sig]} rank")
    fig.suptitle("Pooled TCGA+OncoSG (n=1145): rank-rank association "
                 "(downward slope = negative correlation)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "scatter_pooled_rank.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    forest()
    scatter()
    print("figures written to", OUT)
