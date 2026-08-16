#!/usr/bin/env python3
"""Extra figures: leftover CLDN4/TACSTD2 vs immune neighborhood."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "leftover_spatial" / "tables"
FIG = ROOT / "results" / "leftover_spatial" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png")
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)


def visium_samespot():
    df = pd.read_csv(TAB / "visium_samespot.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharey=True)
    for ax, series, title in [
        (axes[0], "GSE263196", "GSE263196 Visium SCLC"),
        (axes[1], "GSE273378", "GSE273378 Visium stage I LUAD"),
    ]:
        sub = df[df.series == series]
        xs = np.arange(len(sub))
        ax.axhline(0, color="0.6", lw=0.8)
        ax.plot(xs, sub["same_CLDN4_TB_rho"], "o", label="CLDN4 vs T+B", color="#1f4e79")
        ax.plot(xs, sub["same_TACSTD2_TB_rho"], "s", label="TACSTD2 vs T+B", color="#c45c26")
        ax.set_xticks(xs)
        ax.set_xticklabels(sub["section"], rotation=60, ha="right", fontsize=7)
        ax.set_title(title)
        ax.set_ylabel("same-spot Spearman ρ")
        ax.legend(frameon=False, fontsize=7)
    _save(fig, "visium_samespot_rho")


def visium_rings():
    df = pd.read_csv(TAB / "visium_rings.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharey=True)
    for ax, series, title in [
        (axes[0], "GSE263196", "GSE263196 hex rings"),
        (axes[1], "GSE273378", "GSE273378 hex rings"),
    ]:
        sub = df[df.series == series]
        for gene, col, color in [
            ("CLDN4", "CLDN4_vs_neiTB_rho", "#1f4e79"),
            ("TACSTD2", "TACSTD2_vs_neiTB_rho", "#c45c26"),
        ]:
            med = sub.groupby("ring")[col].median()
            ax.plot(med.index, med.values, "-o", label=f"{gene} vs nei T+B", color=color)
        for gene, col, color in [
            ("CLDN4 partial", "partial_CLDN4_vs_neiTB_ctrl_broad_rho", "#7aa6c2"),
            ("TACSTD2 partial", "partial_TACSTD2_vs_neiTB_ctrl_broad_rho", "#e0a07a"),
        ]:
            med = sub.groupby("ring")[col].median()
            ax.plot(med.index, med.values, "--s", label=gene, color=color, ms=4)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xticks([1, 2, 3])
        ax.set_xlabel("hex ring")
        ax.set_title(title)
        ax.set_ylabel("median section Spearman ρ")
        ax.legend(frameon=False, fontsize=6.5)
    _save(fig, "visium_neighborhood_rings")


def geomx_bars():
    df = pd.read_csv(TAB / "geomx_correlations.csv")
    # pick primary subsets
    keep = [
        ("GSE265899", "all_AOI", "GSE265899 all AOI"),
        ("GSE265899", "tumor", "GSE265899 tumor"),
        ("GSE289483", "all_AOI", "GSE289483 all AOI"),
        ("GSE289483", "segmentation=CD45-", "GSE289483 CD45−"),
        ("GSE334014", "all_AOI", "GSE334014 all AOI"),
        ("GSE334014", "PanCK", "GSE334014 PanCK"),
        ("GSE326968", "all_AOI", "GSE326968 CLAD"),
    ]
    rows = []
    labels = []
    for series, subset, lab in keep:
        hit = df[(df.series == series) & (df.subset == subset)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        rows.append([r.get("CLDN4_vs_TB_rho", np.nan), r.get("TACSTD2_vs_TB_rho", np.nan)])
        labels.append(f"{lab}\nn={int(r['n'])}")
    arr = np.array(rows, dtype=float)
    fig, ax = plt.subplots(figsize=(8.6, 3.8))
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, arr[:, 0], w, label="CLDN4 vs T+B", color="#1f4e79")
    ax.bar(x + w / 2, arr[:, 1], w, label="TACSTD2 vs T+B", color="#c45c26")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Spearman ρ")
    ax.set_title("Leftover GeoMx: target vs same-AOI T+B")
    ax.legend(frameon=False)
    _save(fig, "geomx_target_vs_TB")


def cosmx_fov():
    fov = pd.read_csv(TAB / "gse276083_fov_means.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))
    for ax, xcol, title in [
        (axes[0], "CLDN4", "GSE276083 CosMx FOV: CLDN4 vs T+B"),
        (axes[1], "TACSTD2", "GSE276083 CosMx FOV: TACSTD2 vs T+B"),
    ]:
        ax.scatter(fov[xcol], fov["TB"], s=18, alpha=0.75, c="#1f4e79")
        ax.set_xlabel(f"FOV mean {xcol}")
        ax.set_ylabel("FOV mean T+B")
        ax.set_title(title)
    _save(fig, "cosmx_fov_scatter")


def main():
    visium_samespot()
    visium_rings()
    geomx_bars()
    cosmx_fov()
    print("figures written", FIG)


if __name__ == "__main__":
    main()
