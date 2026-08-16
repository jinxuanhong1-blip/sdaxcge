#!/usr/bin/env python3
"""Figures for the TISMO Tacstd2 baseline-vs-ICB replication."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

GENE = "Tacstd2"
UP, DOWN = "#c0392b", "#2c6fbb"


def fig_paired_slopes(cohorts: pd.DataFrame, path: Path) -> None:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 6), gridspec_kw={"width_ratios": [2, 1]})

    for _, r in cohorts.iterrows():
        c = UP if r.delta > 0 else DOWN
        ax.plot([0, 1], [r.mean_baseline, r.mean_icb], color=c, alpha=0.45, lw=1.1,
                marker="o", ms=3)
    ax.plot([0, 1], [cohorts.mean_baseline.mean(), cohorts.mean_icb.mean()],
            color="black", lw=3, marker="o", ms=8, label="mean of cohorts", zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Baseline\n(control arm)", "ICB-treated"])
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylabel(f"{GENE} expression (TISMO normalised units)")
    n_up = int((cohorts.delta > 0).sum())
    ax.set_title(f"{GENE}: paired within-study contrast\n"
                 f"{n_up}/{len(cohorts)} cohorts up  "
                 f"({cohorts.mean_baseline.mean():.2f} \u2192 {cohorts.mean_icb.mean():.2f})")
    ax.legend(frameon=False, loc="upper left")

    order = cohorts.sort_values("delta")
    colors = [UP if d > 0 else DOWN for d in order.delta]
    ax2.barh(np.arange(len(order)), order.delta, color=colors)
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_yticks([])
    ax2.set_xlabel("ICB \u2212 baseline")
    ax2.set_title("Per-cohort difference")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_by_cancer_type(cohorts: pd.DataFrame, path: Path) -> None:
    grp = cohorts.groupby("cancer_type")
    order = grp["delta"].mean().sort_values()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    data = [grp.get_group(ct)["delta"].to_numpy() for ct in order.index]
    ax.axvline(0, color="black", lw=0.8, zorder=0)
    for i, (ct, vals) in enumerate(zip(order.index, data)):
        jitter = np.random.default_rng(0).normal(0, 0.06, len(vals))
        ax.scatter(vals, i + jitter, s=28,
                   color=[UP if v > 0 else DOWN for v in vals], alpha=0.8, zorder=3)
        ax.plot([vals.mean()], [i], marker="|", ms=26, color="black", zorder=4)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{ct}  (n={len(v)})" for ct, v in zip(order.index, data)])
    ax.set_xlabel(f"{GENE}: ICB \u2212 baseline, per cohort")
    ax.set_title(f"{GENE} change by cancer type\n"
                 "lung is represented by 2 LLC cohorts from a single study")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_null(null_df: pd.DataFrame, target: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    ax = axes[0]
    ax.hist(null_df["n_up"], bins=np.arange(0, 68, 2), color="#95a5a6", edgecolor="white")
    ax.axvline(target["n_up"], color=UP, lw=2.5,
               label=f"{GENE} = {target['n_up']}")
    ax.axvline(len(null_df) and null_df["n_up"].median(), color="black", ls="--", lw=1.5,
               label=f"null median = {null_df['n_up'].median():.0f}")
    ax.set_xlabel("cohorts with higher ICB value")
    ax.set_ylabel("random genes")
    ax.set_title("Is 49/64 unusual?")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.hist(null_df["mean_delta"], bins=40, color="#95a5a6", edgecolor="white")
    ax.axvline(target["mean_delta"], color=UP, lw=2.5, label=f"{GENE}")
    ax.axvline(0, color="black", ls="--", lw=1)
    ax.set_xlabel("mean (ICB \u2212 baseline)")
    ax.set_title("Effect size vs random genes")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    p = null_df["wilcoxon_p"].clip(lower=1e-12)
    ax.hist(np.log10(p), bins=40, color="#95a5a6", edgecolor="white")
    ax.axvline(np.log10(target["wilcoxon_p"]), color=UP, lw=2.5, label=f"{GENE}")
    ax.axvline(np.log10(0.05), color="black", ls="--", lw=1, label="p = 0.05")
    ax.set_xlabel("log10 Wilcoxon p")
    frac = float((null_df["wilcoxon_p"] < 0.05).mean())
    ax.set_title(f"{frac:.0%} of random genes reach p < 0.05")
    ax.legend(frameon=False, fontsize=9)

    for ax in axes:
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_reference(ref_df: pd.DataFrame, target: dict, path: Path) -> None:
    df = pd.concat([ref_df, pd.DataFrame([{"gene": GENE, **target}])], ignore_index=True)
    df = df.sort_values("mean_delta")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    colors = [UP if g == GENE else "#7f8c8d" for g in df["gene"]]
    ax.barh(df["gene"], df["mean_delta"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    for y, (d, p) in enumerate(zip(df["mean_delta"], df["wilcoxon_p"])):
        ax.text(d + (0.02 if d >= 0 else -0.02), y, f"p={p:.2g}",
                va="center", ha="left" if d >= 0 else "right", fontsize=8)
    ax.set_xlabel("mean (ICB \u2212 baseline) across cohorts")
    ax.set_title("Tacstd2 against ICB pharmacodynamic markers (Cd274, Cd8a, Ifng, Pdcd1),\n"
                 "epithelial-content markers (Epcam, Krt8) and a housekeeper (Actb)")
    ax.margins(x=0.18)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> int:
    tc.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    cohorts = pd.read_csv(tc.TABLES_DIR / f"cohort_level_{GENE}.csv")

    fig_paired_slopes(cohorts, tc.FIGURES_DIR / f"{GENE}_paired_slopes.png")
    fig_by_cancer_type(cohorts, tc.FIGURES_DIR / f"{GENE}_by_cancer_type.png")

    calib_path = tc.TABLES_DIR / "null_calibration.json"
    null_path = tc.TABLES_DIR / "null_panel_statistics.csv"
    if calib_path.exists() and null_path.exists():
        target = json.loads(calib_path.read_text())["target"]
        fig_null(pd.read_csv(null_path), target, tc.FIGURES_DIR / f"{GENE}_null_calibration.png")
        ref_path = tc.TABLES_DIR / "reference_gene_statistics.csv"
        if ref_path.exists():
            fig_reference(pd.read_csv(ref_path), target,
                          tc.FIGURES_DIR / f"{GENE}_reference_genes.png")

    for p in sorted(tc.FIGURES_DIR.glob("*.png")):
        print(f"  wrote {p.relative_to(tc.REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
