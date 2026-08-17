#!/usr/bin/env python3
"""Extra figures: Cldn4 vs T/exclusion and Cldn4-high + immune-low cuts."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import OUT, high_low_cut

FIG = OUT / "figures"
plt.rcParams.update({"figure.dpi": 140, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def _save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png")
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)


def scatter_cut(df, xcol, ycol, title, name, hue=None):
    if df.empty or xcol not in df or ycol not in df:
        return
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 8:
        return
    mask, meta = high_low_cut(x, y, how="median")
    qmask, qmeta = high_low_cut(x, y, how="quartile")
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    ax.scatter(x[ok & ~mask], y[ok & ~mask], s=6, c="#9aa3ad", alpha=0.35, linewidths=0, label="other")
    ax.scatter(x[mask], y[mask], s=8, c="#c0392b", alpha=0.55, linewidths=0, label=f"Cldn4-high + T-low (median) n={meta['n_hi_lo']}")
    if meta.get("cldn4_cut") is not None:
        ax.axvline(meta["cldn4_cut"], color="#c0392b", ls="--", lw=0.8)
        ax.axhline(meta["immune_cut"], color="#2980b9", ls="--", lw=0.8)
    ax.set_xlabel("Cldn4 log1p(UMI)")
    ax.set_ylabel("T-cell score (mean log1p)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="best", fontsize=7)
    ax.text(
        0.02,
        0.02,
        f"median cut {meta['frac_hi_lo']:.1%}  |  Q3/Q1 cut {qmeta['frac_hi_lo']:.1%} (n={qmeta['n_hi_lo']})",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
    )
    _save(fig, name)


def scrna_figures():
    p = OUT / "scrna_cell_scores.tsv.gz"
    if not p.exists():
        return
    df = pd.read_csv(p, sep="\t")
    for sample, sub in df.groupby("sample"):
        scatter_cut(
            sub,
            "Cldn4",
            "Tscore",
            f"{sample}: all cells",
            f"cut_{sample}_all",
        )
        epi = sub[sub["label"] == "epithelial"]
        scatter_cut(
            epi,
            "Cldn4",
            "Tscore",
            f"{sample}: epithelial cells",
            f"cut_{sample}_epithelial",
        )
    # compartment Cldn4 per sample
    inv = OUT / "scrna_sample_inventory.tsv"
    if inv.exists():
        s = pd.read_csv(inv, sep="\t")
        ok = s[s.get("status", "ok") == "ok"] if "status" in s else s
        if not ok.empty and "epithelial_Cldn4_mean" in ok:
            fig, ax = plt.subplots(figsize=(6.2, 3.6))
            x = np.arange(len(ok))
            ax.bar(x - 0.18, ok["epithelial_Cldn4_mean"], 0.36, label="epithelial Cldn4", color="#2c3e50")
            ax.bar(x + 0.18, ok["tnk_Cldn4_mean"], 0.36, label="T/NK Cldn4", color="#7f8c8d")
            ax.set_xticks(x)
            ax.set_xticklabels(ok["sample"], rotation=30, ha="right")
            ax.set_ylabel("mean log1p Cldn4")
            ax.legend(frameon=False)
            ax.set_title("Cldn4 in epithelium vs T/NK (honest library n)")
            _save(fig, "scrna_cldn4_epi_vs_tnk")
            fig, ax = plt.subplots(figsize=(6.2, 3.6))
            ax.bar(x - 0.18, ok["epithelial_Cldn4_mean"], 0.36, label="epithelial Cldn4", color="#2c3e50")
            ax.bar(x + 0.18, ok["frac_tnk"], 0.36, label="T/NK fraction", color="#2980b9")
            ax.set_xticks(x)
            ax.set_xticklabels(ok["sample"], rotation=30, ha="right")
            ax.set_ylabel("Cldn4 mean  |  T/NK fraction")
            ax.legend(frameon=False)
            ax.set_title("Exclusion-style: epithelial Cldn4 vs T/NK fraction")
            _save(fig, "scrna_cldn4_vs_tnk_frac")


def spatial_figures():
    p = OUT / "spatial_spot_scores.tsv.gz"
    if not p.exists():
        return
    df = pd.read_csv(p, sep="\t")
    for sample, sub in df.groupby("sample"):
        scatter_cut(sub, "Cldn4", "Tscore", f"GSE261890 {sample} spots/bins", f"spatial_cut_{sample}")
        if {"x", "y"}.issubset(sub.columns) and sub["x"].notna().any():
            mask, meta = high_low_cut(sub["Cldn4"], sub["Tscore"], how="median")
            fig, ax = plt.subplots(figsize=(5.0, 4.4))
            ax.scatter(sub.loc[~mask, "x"], sub.loc[~mask, "y"], s=2, c="#bdc3c7", linewidths=0)
            ax.scatter(sub.loc[mask, "x"], sub.loc[mask, "y"], s=3, c="#c0392b", linewidths=0, label=f"Cldn4-high + T-low n={meta['n_hi_lo']}")
            ax.set_aspect("equal")
            ax.set_title(f"{sample}: Cldn4-high + immune-low bins")
            ax.legend(frameon=False, fontsize=7)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            _save(fig, f"spatial_map_cut_{sample}")


def main():
    scrna_figures()
    spatial_figures()
    print("figures in", FIG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
