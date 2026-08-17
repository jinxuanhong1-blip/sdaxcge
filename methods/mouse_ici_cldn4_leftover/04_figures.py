#!/usr/bin/env python3
"""Extra figures: Cldn4 vs T/exclusion and Cldn4-detected + T-low cuts."""
from __future__ import annotations

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


def scatter_cut(df, xcol, ycol, title, name):
    if df.empty or xcol not in df or ycol not in df:
        return None
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 8:
        return None
    mask, meta = high_low_cut(x, y, how="detected")
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    ax.scatter(x[ok & ~mask], y[ok & ~mask], s=6, c="#9aa3ad", alpha=0.35, linewidths=0, label="other")
    ax.scatter(
        x[mask],
        y[mask],
        s=10,
        c="#c0392b",
        alpha=0.7,
        linewidths=0,
        label=f"Cldn4+ & T-low n={meta['n_hi_lo']} / {meta['n_cldn4_pos']} Cldn4+",
    )
    if meta.get("immune_cut") is not None:
        ax.axhline(meta["immune_cut"], color="#2980b9", ls="--", lw=0.8)
    ax.axvline(0, color="#888", ls=":", lw=0.6)
    ax.set_xlabel("Cldn4 log1p(UMI)")
    ax.set_ylabel("T-cell score (mean log1p)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="best", fontsize=7)
    ax.text(
        0.02,
        0.02,
        f"{meta['frac_of_cldn4_pos']:.0%} of Cldn4+ are T-low (cut T≤{meta['immune_cut']:.3g})",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
    )
    _save(fig, name)
    return meta


def scrna_figures():
    p = OUT / "scrna_cell_scores.tsv.gz"
    if not p.exists():
        return
    df = pd.read_csv(p, sep="\t")
    for sample, sub in df.groupby("sample"):
        epi = sub[sub["label"] == "epithelial"]
        scatter_cut(epi, "Cldn4", "Tscore", f"{sample}: epithelial cells", f"cut_{sample}_epithelial")
    inv = OUT / "scrna_sample_inventory.tsv"
    if not inv.exists():
        return
    s = pd.read_csv(inv, sep="\t")
    ok = s[s["status"] == "ok"] if "status" in s else s
    if ok.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    x = np.arange(len(ok))
    ax.bar(x - 0.18, ok["epithelial_Cldn4_mean"], 0.36, label="epithelial Cldn4", color="#2c3e50")
    ax.bar(x + 0.18, ok["tnk_Cldn4_mean"], 0.36, label="T/NK Cldn4", color="#7f8c8d")
    ax.set_xticks(x)
    ax.set_xticklabels(ok["sample"], rotation=30, ha="right")
    ax.set_ylabel("mean log1p Cldn4")
    ax.legend(frameon=False)
    ax.set_title("Cldn4 in epithelium vs T/NK (library n as labeled)")
    _save(fig, "scrna_cldn4_epi_vs_tnk")
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for gse, sub in ok.groupby("gse"):
        ax.scatter(sub["epithelial_Cldn4_mean"], sub["frac_tnk"], s=50, label=gse)
        for _, r in sub.iterrows():
            ax.annotate(r["sample"], (r["epithelial_Cldn4_mean"], r["frac_tnk"]), fontsize=6, xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("epithelial Cldn4 mean")
    ax.set_ylabel("T/NK fraction")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("Exclusion-style: epi Cldn4 vs T/NK fraction")
    _save(fig, "scrna_cldn4_vs_tnk_frac")


def spatial_figures():
    p = OUT / "spatial_spot_scores.tsv.gz"
    if not p.exists():
        return
    df = pd.read_csv(p, sep="\t")
    for sample, sub in df.groupby("sample"):
        scatter_cut(sub, "Cldn4", "Tscore", f"GSE261890 {sample} L7 bins", f"spatial_cut_{sample}")
        if {"x", "y"}.issubset(sub.columns) and sub["x"].notna().any():
            mask, meta = high_low_cut(sub["Cldn4"], sub["Tscore"], how="detected")
            fig, ax = plt.subplots(figsize=(5.0, 4.4))
            ax.scatter(sub.loc[~mask, "x"], sub.loc[~mask, "y"], s=2, c="#bdc3c7", linewidths=0)
            ax.scatter(
                sub.loc[mask, "x"],
                sub.loc[mask, "y"],
                s=6,
                c="#c0392b",
                linewidths=0,
                label=f"Cldn4+ & T-low n={meta['n_hi_lo']}",
            )
            ax.set_aspect("equal")
            ax.set_title(f"{sample}: Cldn4-detected + T-low bins")
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
