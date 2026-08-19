#!/usr/bin/env python3
"""Lead paper figure: barrier index and infiltration-depth AUC, not nearest-µm."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
            "savefig.bbox": "tight",
        }
    )


def _jitter(n, rng, scale=0.10):
    return rng.uniform(-scale, scale, n)


def plot_lead(sample_csv: Path, domain_csv: Path, curve_csv: Path, out: Path):
    _style()
    sam = pd.read_csv(sample_csv)
    dom = pd.read_csv(domain_csv)
    cur = pd.read_csv(curve_csv)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.4))

    # ---- A: CosMx barrier, 8 sections ----
    ax = axes[0, 0]
    cs = sam[sam.platform == "cosmx"]
    for j, col, color, lab in (
        (0, "barrier_high", "#b2182b", "CLDN4-high"),
        (1, "barrier_low", "#2166ac", "CLDN4-low"),
    ):
        v = pd.to_numeric(cs[col], errors="coerce").dropna().to_numpy()
        ax.scatter(np.full(len(v), j) + _jitter(len(v), rng), v, s=28, color=color, zorder=3, label=lab)
        ax.errorbar(j, np.mean(v), yerr=np.std(v, ddof=1) / np.sqrt(len(v)), fmt="o", color="k", ms=6, zorder=4)
    ax.axhline(0, color="0.55", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-high", "CLDN4-low"])
    ax.set_ylabel("Barrier index")
    ax.set_title("A  CosMx barrier (8 sections)")
    ax.legend(frameon=False, loc="lower left")

    # ---- B: CosMx AUC, 8 sections ----
    ax = axes[0, 1]
    for j, col, color, lab in (
        (0, "auc_high", "#b2182b", "CLDN4-high"),
        (1, "auc_low", "#2166ac", "CLDN4-low"),
    ):
        v = pd.to_numeric(cs[col], errors="coerce").dropna().to_numpy() / 1000.0
        ax.scatter(np.full(len(v), j) + _jitter(len(v), rng, 0.08), v, s=28, color=color, zorder=3, label=lab)
        ax.errorbar(j, np.mean(v), yerr=np.std(v, ddof=1) / np.sqrt(len(v)), fmt="o", color="k", ms=6, zorder=4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-high", "CLDN4-low"])
    ax.set_ylabel("CD8 AUC₀–₂₀₀ (10³ cells·µm / mm²)")
    ax.set_title("B  CosMx infiltration-depth AUC")

    # ---- C: Visium barrier ----
    ax = axes[1, 0]
    vs = sam[sam.platform == "visium"]
    for j, col, color, lab in (
        (0, "barrier_high", "#b2182b", "CLDN4-high"),
        (1, "barrier_low", "#2166ac", "CLDN4-low"),
    ):
        v = pd.to_numeric(vs[col], errors="coerce").dropna().to_numpy()
        ax.scatter(np.full(len(v), j) + _jitter(len(v), rng, 0.12), v, s=16, color=color, alpha=0.8, zorder=3)
        if len(v):
            ax.errorbar(j, np.mean(v), yerr=np.std(v, ddof=1) / np.sqrt(len(v)), fmt="o", color="k", ms=6, zorder=4)
    ax.axhline(0, color="0.55", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4-high", "CLDN4-low"])
    ax.set_ylabel("Barrier index (100 µm rim)")
    ax.set_title("C  Visium LUAD barrier (26 slides)")

    # ---- D: infiltration-depth curves CosMx ----
    ax = axes[1, 1]
    for cls, color, lab in (("high", "#b2182b", "CLDN4-high"), ("low", "#2166ac", "CLDN4-low")):
        sub = cur[(cur.platform == "cosmx") & (cur.cldn4_class == cls)]
        if sub.empty:
            continue
        g = sub.groupby("dist_um", as_index=False).cd8.agg(["mean", "sem", "count"])
        ax.plot(g["dist_um"], g["mean"], color=color, lw=2.0, label=lab)
        ax.fill_between(g["dist_um"], g["mean"] - g["sem"], g["mean"] + g["sem"], color=color, alpha=0.18)
    ax.axvline(0, color="0.35", lw=0.8)
    ax.axvline(50, color="0.6", ls="--", lw=0.8)
    ax.set_xlim(0, 200)
    ax.set_xlabel("Inward distance from margin (µm)")
    ax.set_ylabel("CD8+ cells / mm²")
    ax.set_title("D  CosMx depth curve (0–200 µm)")
    ax.legend(frameon=False)

    fig.suptitle("Barrier index and CD8 infiltration-depth AUC: CLDN4-high vs CLDN4-low tumor domains", y=1.02)
    fig.tight_layout()
    fig.savefig(out / "fig1_barrier_and_auc.png")
    fig.savefig(out / "fig1_barrier_and_auc.pdf")
    plt.close(fig)

    # keep a standalone depth-curve figure as fig2
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    for ax, plat, title, ylab in (
        (axes[0], "cosmx", "CosMx NSCLC (8 sections)", "CD8+ cells / mm²"),
        (axes[1], "visium", "Visium LUAD GSE307534 (invasive)", "CD8A/B module (log1p)"),
    ):
        for cls, color, lab in (("high", "#b2182b", "CLDN4-high domains"), ("low", "#2166ac", "CLDN4-low domains")):
            sub = cur[(cur.platform == plat) & (cur.cldn4_class == cls)]
            if sub.empty:
                continue
            g = sub.groupby("dist_um", as_index=False).cd8.agg(["mean", "sem"])
            ax.plot(g["dist_um"], g["mean"], color=color, lw=2.0, label=lab)
            ax.fill_between(g["dist_um"], g["mean"] - g["sem"], g["mean"] + g["sem"], color=color, alpha=0.18)
        ax.axvline(50, color="0.6", ls="--", lw=0.8)
        ax.set_xlim(0, 200)
        ax.set_xlabel("Inward distance from tumor margin (µm)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(frameon=False)
    fig.suptitle("CD8 infiltration-depth AUC support: curves by CLDN4 domain class", y=1.03)
    fig.tight_layout()
    fig.savefig(out / "fig2_infiltration_depth_curves.png")
    fig.savefig(out / "fig2_infiltration_depth_curves.pdf")
    plt.close(fig)
    print("wrote lead figures to", out)


if __name__ == "__main__":
    root = Path("/workspace/results")
    plot_lead(root / "tables" / "sample_metrics.csv", root / "tables" / "domain_metrics.csv", root / "tables" / "infiltration_curves.csv", root / "figures")
