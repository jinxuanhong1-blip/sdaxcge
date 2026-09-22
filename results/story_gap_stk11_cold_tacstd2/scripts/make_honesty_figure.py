#!/usr/bin/env python3
"""Display-only honesty matrix for the disease-context slide. No statistics recomputed."""
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = [
        ("Skoulidis 2018\n(literature)", "yes", "KL ORR 7.4%\nCD8↓ PROSPECT", "na", "not about\nantigen", "Human lit"),
        ("MSK KRAS±STK11\nDCB (PR #49)", "yes", "OR 0.24 / 0.38\nMH≈0.36", "na", "genotype→ICI\nonly", "Human bulk"),
        ("GSE137244\nKL vs KP lines", "soft", "IFN soft only\n(not the lock)", "yes", "Tacstd2 +3.24\nCldn4 +5.57\nMW p=0.00794", "Mouse\ncell-line bulk"),
        ("GSE179500\nLkb1XTR", "na", "sorted tumor\nbulk", "yes", "LKB1-off↑\npadj 0.032", "Mouse\ntumor bulk"),
        ("GSE179502\nscRNA", "no", "immune FACS\ndepleted", "partial", "Cldn4↓ restore\nTacstd2 p=0.10", "Mouse\nscRNA"),
        ("GSE180963\nscRNA", "no", "immune ~82%\nboth; n=1", "partial", "epi Tacstd2\nselectivity", "Mouse\nscRNA"),
        ("TCGA-LUAD\nSTK11", "lit", "cold via lit;\nantigen≠cold", "no", "TACSTD2 −0.51\nCLDN4 −0.47\n(REVERSE)", "Human bulk"),
        ("GSE280232\nscRNA", "na", "post-ICB\nT-cell rich", "no", "epi yes;\ngenotype unpowered", "Human\nscRNA"),
    ]
    color = {
        "yes": "#2E7D32",
        "partial": "#F9A825",
        "soft": "#C0CA33",
        "no": "#C62828",
        "na": "#9E9E9E",
        "lit": "#546E7A",
    }
    fig, ax = plt.subplots(figsize=(11.5, 8.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.text(5.0, 9.55, "Disease context honesty matrix — bulk vs scRNA",
            ha="center", va="top", fontsize=14, fontweight="bold", color="#1A237E")
    ax.text(5.0, 9.15, "Green = supports slide limb · Yellow = partial · Red = reverse / no · Gray = N/A",
            ha="center", va="top", fontsize=9, color="#455A64")
    for x, t in [(2.8, "COLD / ICI-resistant"), (6.6, "Tacstd2 / Cldn4 HIGH"), (9.0, "Assay")]:
        ax.text(x, 8.7, t, ha="center", va="center", fontsize=10, fontweight="bold")
    row_h = 0.95
    y0 = 8.2
    for i, (lab, cs, ct, as_, at, assay) in enumerate(cells):
        y = y0 - i * row_h
        ax.text(0.15, y - 0.35, lab, ha="left", va="center", fontsize=8, fontweight="bold")
        for x, status, text in [(2.8, cs, ct), (6.6, as_, at)]:
            box = FancyBboxPatch(
                (x - 1.55, y - 0.85), 3.1, 0.8,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                facecolor=color[status], edgecolor="white", alpha=0.85, linewidth=1.5,
            )
            ax.add_patch(box)
            ax.text(x, y - 0.45, text, ha="center", va="center", fontsize=7.2,
                    color="white", fontweight="bold")
        ax.text(9.0, y - 0.45, assay, ha="center", va="center", fontsize=7.5, color="#263238")
    ax.text(
        5.0, 0.35,
        "LOCK antigen on GSE137244 (mouse cell-line bulk). Human TCGA bulk is REVERSE for antigen. Keep bulk ≠ scRNA.",
        ha="center", va="center", fontsize=8.5, style="italic", color="#B71C1C",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFEBEE", edgecolor="#EF9A9A"),
    )
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"honesty_matrix.{ext}", dpi=160, bbox_inches="tight", facecolor="white")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
