#!/usr/bin/env python3
"""Redraw PPT figures from copied CosMx paper-funnel numbers.

Does not load the CosMx h5ad. Does not recompute 0.36/0.52.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_funnel_cosmx_spatial" / "figures"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    radii = ["50 µm", "100 µm"]
    ratios = [0.36, 0.52]
    bars = ax.bar(radii, ratios, color=["#2c7bb6", "#abd9e9"], edgecolor="black", width=0.55)
    ax.axhline(1.0, color="#666666", ls="--", lw=1, label="parity (hi=lo)")
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("CLDN4-high / CLDN4-low\ncytotoxic neighbor ratio")
    ax.set_title("LOCKED CosMx CLDN4 exclusion\n(not recomputed on this PR)")
    for b, r in zip(bars, ratios):
        ax.text(
            b.get_x() + b.get_width() / 2,
            r + 0.03,
            f"{r:.2f}",
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )
    ax.text(
        0.5,
        0.08,
        "8/8 sections · 5/5 patients · sign P=0.031\n"
        "exclusion, not muzzling (GZMB/PRF1/NKG7/IFNG hi/lo 1.11–1.22, 0/8↓)",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="#333333",
    )
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_locked_cldn4_cytotoxic_ratios.png", dpi=200)
    fig.savefig(OUT / "fig1_locked_cldn4_cytotoxic_ratios.pdf")
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    ax = axes[0]
    rads = np.array([10, 20, 50, 100])
    cd8nk = [0.455, 0.670, 0.854, 0.947]
    imm = [0.519, 0.643, 0.848, 0.957]
    ax.plot(rads, cd8nk, "o-", color="#d7191c", lw=2, label="CD8+NK count")
    ax.plot(rads, imm, "s--", color="#fdae61", lw=2, label="immune fraction")
    ax.axhline(1.0, color="#666", ls=":", lw=1)
    ax.set_xlabel("radius (µm)")
    ax.set_ylabel("TACSTD2-high / low ratio")
    ax.set_title("TACSTD2-high cold neighbors\n(PR #726; median split)")
    ax.set_xticks(rads)
    ax.set_ylim(0.3, 1.15)
    for x, y, nsec in zip(rads, cd8nk, ["8/8", "8/8", "7/8", "4/8"]):
        ax.annotate(nsec, (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=7, color="#d7191c")
    ax.legend(frameon=False, fontsize=8)
    ax.text(0.02, 0.02, "8/8 & 5/5 at 10–20 µm; sign P=0.031", transform=ax.transAxes, fontsize=8)

    ax = axes[1]
    labels = ["CD8+NK\n10 µm", "CD8+NK\n20 µm", "Immune\n10 µm", "Immune\n20 µm"]
    vals = [0.455, 0.670, 0.519, 0.643]
    colors = ["#d7191c", "#d7191c", "#fdae61", "#fdae61"]
    bars = ax.bar(labels, vals, color=colors, edgecolor="black")
    ax.axhline(1.0, color="#666", ls="--", lw=1)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("high / low ratio")
    ax.set_title("Concordant short-range\n(8/8 sections, 5/5 patients)")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_tacstd2_cold_neighbors.png", dpi=200)
    fig.savefig(OUT / "fig2_tacstd2_cold_neighbors.pdf")
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7))
    ax = axes[0]
    patients = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
    rhos = [-0.453, -0.057, -0.201, -0.361, -0.067]
    colors = ["#1b9e77" if r < 0 else "#d95f02" for r in rhos]
    ax.barh(patients[::-1], rhos[::-1], color=colors[::-1], edgecolor="black")
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(-0.228, color="#7570b3", ls="--", lw=1.5, label="patient-mean ρ=−0.228")
    ax.set_xlabel("Spearman ρ (program vs CD8 count @ 50 µm)")
    ax.set_title("CLDN4 epithelial program cold\n(PR #643; 5/5 patients)")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.text(
        0.98,
        0.05,
        "toroidal-shift p=0.005\nhigh/low CD8 count ≈0.61×",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
    )

    ax = axes[1]
    labels2 = ["Program\nCD8 hi/lo", "Contact OR\nCLDN4-hi", "Locked\n50 µm", "Locked\n100 µm"]
    vals2 = [0.61, 0.584, 0.36, 0.52]
    cols2 = ["#7570b3", "#1b9e77", "#2c7bb6", "#abd9e9"]
    bars = ax.bar(labels2, vals2, color=cols2, edgecolor="black")
    ax.axhline(1.0, color="#666", ls="--", lw=1)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("ratio or odds ratio (<1 = colder)")
    ax.set_title("CLDN4 / program cold metrics\n(CosMx He2022)")
    for b, v in zip(bars, vals2):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.2f}", ha="center", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_cldn4_tj_program_cold.png", dpi=200)
    fig.savefig(OUT / "fig3_cldn4_tj_program_cold.pdf")
    plt.close()

    fig = plt.figure(figsize=(10.5, 4.2))
    gs = fig.add_gridspec(1, 3, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    bars = ax1.bar(["50 µm", "100 µm"], [0.36, 0.52], color=["#2c7bb6", "#abd9e9"], edgecolor="k")
    ax1.axhline(1, color="#666", ls="--", lw=1)
    ax1.set_ylim(0, 1.2)
    ax1.set_title("A. LOCKED CLDN4\ncytotoxic neighbors", fontsize=10)
    ax1.set_ylabel("high/low ratio")
    for b, v in zip(bars, [0.36, 0.52]):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.04, f"{v:.2f}", ha="center", fontweight="bold")
    ax1.text(0.5, 0.02, "8/8 · 5/5 · P=0.031", transform=ax1.transAxes, ha="center", fontsize=8)

    bars = ax2.bar(
        ["CD8+NK\n10 µm", "CD8+NK\n20 µm", "Imm\n10 µm", "Imm\n20 µm"],
        [0.455, 0.670, 0.519, 0.643],
        color=["#d7191c", "#d7191c", "#fdae61", "#fdae61"],
        edgecolor="k",
    )
    ax2.axhline(1, color="#666", ls="--", lw=1)
    ax2.set_ylim(0, 1.2)
    ax2.set_title("B. TACSTD2-high\ncold neighbors", fontsize=10)
    ax2.set_ylabel("high/low ratio")
    for b, v in zip(bars, [0.455, 0.670, 0.519, 0.643]):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.04, f"{v:.2f}", ha="center", fontsize=8, fontweight="bold")
    ax2.text(0.5, 0.02, "8/8 · 5/5 · P=0.031 (#726)", transform=ax2.transAxes, ha="center", fontsize=8)

    bars = ax3.bar(
        ["Program\nρ (abs)", "Program\nCD8 ratio", "Contact\nOR"],
        [0.228, 0.61, 0.584],
        color=["#7570b3", "#7570b3", "#1b9e77"],
        edgecolor="k",
    )
    ax3.set_ylim(0, 1.0)
    ax3.set_title("C. CLDN4/TJ program cold\n(ρ, ratio, OR)", fontsize=10)
    ax3.set_ylabel("effect size")
    for b, v, lab in zip(bars, [0.228, 0.61, 0.584], ["|ρ|=0.23", "0.61×", "OR 0.58"]):
        ax3.text(b.get_x() + b.get_width() / 2, v + 0.03, lab, ha="center", fontsize=8, fontweight="bold")
    ax3.text(0.5, 0.02, "5/5 · shift p=0.005 (#643/#567)", transform=ax3.transAxes, ha="center", fontsize=8)

    fig.suptitle(
        "CosMx He2022 spatial corroboration — PPT funnel (numbers copied; 0.36/0.52 locked)",
        fontsize=11,
        y=1.02,
    )
    fig.savefig(OUT / "fig0_ppt_combined_strip.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / "fig0_ppt_combined_strip.pdf", bbox_inches="tight")
    plt.close()
    print(f"wrote figures under {OUT}")


if __name__ == "__main__":
    main()
