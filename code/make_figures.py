#!/usr/bin/env python3
"""Figures for Claim A1: ESTIMATE vs ABSOLUTE purity-adjusted partial Spearman."""
import os
from math import exp, log, sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "claim_A1")

SIGS = ["immune_tcell_effector", "cytotoxic", "exhaustion"]
SIG_LABEL = {
    "immune_tcell_effector": "Immune (T-cell effector)",
    "cytotoxic": "Cytotoxic",
    "exhaustion": "Exhaustion / checkpoint",
}


def fisher_ci(r, n):
    if n is None or n <= 3 or r is None or np.isnan(r):
        return np.nan, np.nan
    z = 0.5 * log((1 + r) / (1 - r))
    se = 1.0 / sqrt(n - 3)
    lo, hi = z - 1.96 * se, z + 1.96 * se
    back = lambda x: (exp(2 * x) - 1) / (exp(2 * x) + 1)
    return back(lo), back(hi)


def forest_methods():
    pool = pd.read_csv(os.path.join(OUT, "pooled_correlations.csv"))
    keep_pools = ["PanCan_LUAD_LUSC", "PanCan_plus_OncoSG",
                  "Firehose_LUAD_LUSC", "Firehose_plus_OncoSG"]
    keep_methods = ["ESTIMATE", "ABSOLUTE", "CPE", "ESTIMATE_or_ABSOLUTE"]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 6.2), sharex=True)
    for ax, sig in zip(axes, SIGS):
        sub = pool[(pool.signature == sig) & (pool.pool.isin(keep_pools))
                   & (pool.purity_method.isin(keep_methods))].copy()
        # drop ESTIMATE+OncoSG rows that silently drop OncoSG (same n as TCGA-only)
        sub = sub[~((sub.purity_method == "ESTIMATE") & sub.pool.str.contains("OncoSG"))]
        sub["label"] = sub["pool"] + " / " + sub["purity_method"]
        sub = sub.iloc[::-1]
        ys = range(len(sub))
        for y, (_, r) in zip(ys, sub.iterrows()):
            lo, hi = fisher_ci(r["partial_rho"], r["n_partial"])
            color = "#b2182b" if "OncoSG" in r["pool"] else "#2166ac"
            ax.errorbar(r["partial_rho"], y,
                        xerr=[[r["partial_rho"] - lo], [hi - r["partial_rho"]]],
                        fmt="o", color=color, capsize=2.5, ms=5)
            ax.text(0.01, y, f"n={int(r['n_partial'])}",
                    transform=ax.get_yaxis_transform(),
                    va="center", fontsize=6.5, color="#555")
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set_yticks(list(ys))
        ax.set_yticklabels(sub["label"].tolist(), fontsize=7)
        ax.set_title(SIG_LABEL[sig], fontsize=10)
        ax.set_xlabel("Partial Spearman rho (purity-adj)")
    fig.suptitle("TACSTD2 vs immune signatures after ESTIMATE / ABSOLUTE / CPE "
                 "purity adjustment (negative = claim direction)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "forest_partial_spearman.png"), dpi=150)
    plt.close(fig)


def scatter():
    frames = []
    for cid in ["TCGA_LUAD_PanCan", "TCGA_LUSC_PanCan", "OncoSG_LUAD"]:
        path = os.path.join(OUT, f"sample_data_{cid}.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
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
    fig.suptitle("PanCan LUAD+LUSC+OncoSG rank-rank (unadjusted; downward = negative)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "scatter_pooled_rank.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    forest_methods()
    scatter()
    print("figures written to", OUT)
