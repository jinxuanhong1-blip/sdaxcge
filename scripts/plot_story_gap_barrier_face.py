#!/usr/bin/env python3
"""Plot the crude Q4 vs Q1 barrier-ligand deltas copied from PR #716.

Does not refit cohorts. PR #711 is a different score and is not drawn on this axis.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results" / "story_gap_barrier_face" / "evidence_table.tsv"
OUT_DIR = ROOT / "results" / "story_gap_barrier_face" / "figures"

LIGANDS = ["CDH1", "NECTIN2", "F11R", "LGALS9"]
GATES = ["TACSTD2", "CLDN4"]
COLORS = {"TACSTD2": "#0072B2", "CLDN4": "#E69F00"}


def load_bars():
    rows = {}
    with TABLE.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            rec = dict(zip(header, line.rstrip("\n").split("\t")))
            if (
                rec["block"] == "same_scale"
                and rec["split"] == "q4q1"
                and rec["contrast"] == "crude"
                and rec["score"] == "expr_prop_pp"
                and rec["item"] in LIGANDS
            ):
                rows[(rec["gate"], rec["item"])] = float(rec["delta"])
    missing = [(g, lig) for g in GATES for lig in LIGANDS if (g, lig) not in rows]
    if missing:
        raise SystemExit(f"missing rows: {missing}")
    return rows


def main():
    rows = load_bars()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(LIGANDS))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=160)

    for i, gate in enumerate(GATES):
        vals = [rows[(gate, lig)] for lig in LIGANDS]
        offset = -width / 2 if i == 0 else width / 2
        bars = ax.bar(
            x + offset,
            vals,
            width,
            label=f"{gate} Q4 vs Q1",
            color=COLORS[gate],
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
        for bar, val in zip(bars, vals):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val + 0.6,
                f"{val:.1f}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#222222",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(LIGANDS, fontsize=11)
    ax.set_ylabel("Expression-proportion Δ (percentage points)", fontsize=10)
    ax.set_ylim(0, 42)
    ax.set_title(
        "Same barrier face, crude Δ\n"
        "Family mean  +25.19 pp (TACSTD2, n=63)    +25.60 pp (CLDN4, n=64)",
        fontsize=12,
        loc="left",
        pad=10,
    )
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle=":", color="#cccccc", zorder=0)
    ax.set_axisbelow(True)
    fig.text(
        0.01,
        0.01,
        "Concordant-4 malignant → T/NK. PR #716 crude Q4 vs Q1 only. "
        "PR #711 CellPhoneDB +0.25 is a different gate and is not on this axis.",
        fontsize=7.5,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    png = OUT_DIR / "barrier_face_crude_delta.png"
    pdf = OUT_DIR / "barrier_face_crude_delta.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    print(png)


if __name__ == "__main__":
    main()
