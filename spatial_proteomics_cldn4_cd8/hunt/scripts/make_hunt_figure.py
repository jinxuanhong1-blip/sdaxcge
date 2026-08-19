#!/usr/bin/env python3
"""Hunt-funnel figure. No spatial metrics — there is no qualified dataset."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = Path(__file__).resolve().parents[2] / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Counts are hunt tallies, not biology.
STEPS = [
    ("Public series / kits / papers screened", 41, "#4C6A92"),
    ("CD8 protein + public cell XY", 14, "#6B8F71"),
    ("CLDN4 protein on the same panel", 0, "#C44E52"),
    ("Lung among those 0", 0, "#C44E52"),
    ("Non-lung epithelial fallback with both", 0, "#C44E52"),
]


def main() -> None:
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    labels = [s[0] for s in STEPS]
    values = [s[1] for s in STEPS]
    colors = [s[2] for s in STEPS]
    y = list(range(len(STEPS) - 1, -1, -1))
    bars = ax.barh(y, values, color=colors, edgecolor="black", linewidth=0.6, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Count")
    ax.set_xlim(0, 48)
    ax.set_title(
        "Public spatial proteomics hunt: CLDN4 protein + CD8 protein + XY\n"
        "0 qualified datasets (0 lung)",
        loc="left",
        fontsize=12,
    )
    for bar, val in zip(bars, values):
        ax.text(val + 0.6, bar.get_y() + bar.get_height() / 2, str(val), va="center", fontsize=11)
    ax.axvline(0, color="0.3", linewidth=0.5)
    note = (
        "CD8+XY without CLDN4 is the typical IMC/CODEX/CosMx-protein/Xenium-protein panel.\n"
        "No nearest-CD8 / radius / mixing metrics were computed."
    )
    ax.text(0.0, -1.35, note, transform=ax.get_xaxis_transform(), fontsize=8.5, color="0.25")
    legend = [
        mpatches.Patch(color="#4C6A92", label="Screened"),
        mpatches.Patch(color="#6B8F71", label="CD8 protein present"),
        mpatches.Patch(color="#C44E52", label="CLDN4+CD8 same panel = 0"),
    ]
    ax.legend(handles=legend, loc="lower right", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.42, right=0.98, top=0.86, bottom=0.22)
    png = OUT / "fig1_hunt_funnel.png"
    pdf = OUT / "fig1_hunt_funnel.pdf"
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
