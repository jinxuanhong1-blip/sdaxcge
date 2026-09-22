#!/usr/bin/env python3
"""PPT figures for the TJ gene screen story-gap table. Compile-only numbers."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT = Path(__file__).resolve().parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Quoted from PR #748 / #644 — do not invent
C4 = [
    ("CLDN4", -0.478, 1, "lead"),
    ("CLDN7", -0.352, 2, "ns_sep"),
    ("CLDN3", -0.346, 3, "ok"),
    ("CDH1", -0.277, 4, "ok"),
    ("OCLN", -0.157, 5, "ok"),
    ("F11R", -0.070, 6, "ok"),
    ("CLDN1", +0.225, None, "opp"),
]
TCGA = [
    ("F11R", -0.131, 1),
    ("CDH1", -0.126, 2),
    ("CLDN7", -0.109, 3),
    ("OCLN", -0.101, 4),
    ("CLDN4", -0.081, 5),
    ("CLDN3", -0.060, 6),
]
COSMX = [
    ("CDH1", 0.786, True),
    ("CLDN4", 0.800, False),
]


def fig_why_wins():
    fig = plt.figure(figsize=(11.5, 7.2), facecolor="#f7f4ef")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.28,
                          left=0.07, right=0.98, top=0.88, bottom=0.08)

    fig.suptitle(
        "Why CLDN4 wins the screen — after Tacstd2 → TJ (public)",
        fontsize=15, fontweight="bold", color="#1a1a1a", y=0.97,
    )
    fig.text(
        0.5, 0.915,
        "Concordant-4 membership screen leads on CLDN4. TCGA/CosMx do not uniquely pin CLDN4. Never fabricate off-panel CosMx.",
        ha="center", fontsize=9, color="#444",
    )

    # --- panel A: C4 ranks ---
    ax = fig.add_subplot(gs[0, 0])
    genes = [g for g, *_ in C4]
    rhos = [r for _, r, *_ in C4]
    colors = []
    for g, r, rank, tag in C4:
        if g == "CLDN4":
            colors.append("#b45309")
        elif tag == "opp":
            colors.append("#6b7280")
        else:
            colors.append("#78716c")
    y = np.arange(len(genes))[::-1]
    ax.barh(y, rhos, color=colors, edgecolor="none", height=0.72)
    ax.axvline(0, color="#1a1a1a", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel("KRT-partial Spearman ρ vs patient T/NK  (n=65)", fontsize=9)
    ax.set_title("A. Concordant-4 — CLDN4 #1", fontsize=11, fontweight="bold", loc="left")
    for yi, (g, r, rank, tag) in zip(y, C4):
        label = f"{r:+.3f}" if rank is None else f"#{rank}  {r:+.3f}"
        ax.text(r + (0.02 if r >= 0 else -0.02), yi, label,
                va="center", ha="left" if r >= 0 else "right", fontsize=8, color="#1a1a1a")
    ax.set_xlim(-0.65, 0.40)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f7f4ef")

    # --- panel B: TCGA ---
    ax = fig.add_subplot(gs[0, 1])
    genes = [g for g, *_ in TCGA]
    rhos = [r for _, r, *_ in TCGA]
    colors = ["#b45309" if g == "CLDN4" else "#78716c" for g in genes]
    y = np.arange(len(genes))[::-1]
    ax.barh(y, rhos, color=colors, edgecolor="none", height=0.72)
    ax.axvline(0, color="#1a1a1a", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel("KRT-partial Spearman ρ vs CD8  (8 cohorts, N=3945)", fontsize=9)
    ax.set_title("B. TCGA — CLDN4 #5 / 6 (F11R leads)", fontsize=11, fontweight="bold", loc="left")
    for yi, (g, r, rank) in zip(y, TCGA):
        ax.text(r - 0.004, yi, f"#{rank}  {r:+.3f}", va="center", ha="right", fontsize=8)
    ax.set_xlim(-0.18, 0.02)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f7f4ef")

    # --- panel C: CosMx ---
    ax = fig.add_subplot(gs[1, 0])
    genes = [g for g, *_ in COSMX]
    ratios = [r for _, r, *_ in COSMX]
    colors = ["#0f766e" if ok else "#b45309" for _, _, ok in COSMX]
    x = np.arange(len(genes))
    bars = ax.bar(x, ratios, color=colors, edgecolor="none", width=0.55)
    ax.axhline(1.0, color="#1a1a1a", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=10)
    ax.set_ylabel("Q4/Q1 CD8+NK neighbor ratio @50 µm", fontsize=9)
    ax.set_title("C. CosMx — only CLDN4 & CDH1 on 960; CDH1 colder", fontsize=11, fontweight="bold", loc="left")
    ax.set_ylim(0.5, 1.15)
    for xi, (g, r, ok) in enumerate(COSMX):
        note = "8/8 · 5/5" if ok else "4/8 · 3/5"
        ax.text(xi, r + 0.025, f"{r:.3f}\n{note}", ha="center", va="bottom", fontsize=8)
    ax.text(0.5, 0.08, "OFF panel (not invented): CLDN1/3/5/7 · OCLN · F11R · TJP1",
            transform=ax.transAxes, ha="center", fontsize=8, color="#6b7280")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f7f4ef")

    # --- panel D: verdict table ---
    ax = fig.add_subplot(gs[1, 1])
    ax.set_axis_off()
    ax.set_title("D. Verdict for the talk track", fontsize=11, fontweight="bold", loc="left")
    rows = [
        ["Layer", "Lead", "CLDN4", "Say on slide?"],
        ["Concordant-4", "CLDN4 ρ=−0.478", "#1 / scored", "YES — screen win"],
        ["TCGA CD8", "F11R ρ=−0.131", "#5 / 6", "NO unique pin"],
        ["CosMx same-cut", "CDH1 0.786", "#2 / 2 on-panel", "NO unique pin"],
        ["Locked CosMx", "CLDN4 0.36/0.52", "prior lock", "keep; do not swap"],
        ["CLDN1 / 5 / TJP1", "—", "opp / GAP / GAP", "do not invent"],
    ]
    table = ax.table(
        cellText=rows[1:],
        colLabels=rows[0],
        loc="center",
        cellLoc="left",
        colLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.05, 1.55)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#d6d3d1")
        if r == 0:
            cell.set_facecolor("#1c1917")
            cell.set_text_props(color="white", fontweight="bold")
        elif r == 1:
            cell.set_facecolor("#ffedd5")
        else:
            cell.set_facecolor("#fafaf9")

    fig.text(
        0.5, 0.015,
        "Quoted: PR #748 (C4/TCGA/CosMx 6-gene), #644 (CLDN1), #741 (DEG), #698 (locked CosMx 0.36/0.52), #643 (panel). Compile-only.",
        ha="center", fontsize=7.5, color="#78716c",
    )

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"why_cldn4_wins_screen.{ext}", dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def fig_compact_table():
    fig, ax = plt.subplots(figsize=(12.2, 5.4), facecolor="#f7f4ef")
    ax.set_axis_off()
    ax.set_title(
        "TJ gene screen after Tacstd2-high DEG — why CLDN4 wins (concordant-4)",
        fontsize=13, fontweight="bold", pad=18,
    )
    rows = [
        ["Gene", "DEG logFC", "C4 ρ (KRT)", "C4 rank", "TCGA ρ (KRT·CD8)", "TCGA rank", "CosMx Q4/Q1@50µm"],
        ["CLDN4", "+1.632", "−0.478", "#1 WINS", "−0.081", "#5", "0.800 (#2 of 2)"],
        ["CLDN7", "+0.713", "−0.352", "#2", "−0.109", "#3", "OFF 960"],
        ["CLDN3", "+0.652", "−0.346", "#3", "−0.060", "#6", "OFF 960"],
        ["CDH1", "+1.032", "−0.277", "#4", "−0.126", "#2", "0.786 (#1; 8/8)"],
        ["OCLN", "+0.989", "−0.157", "#5", "−0.101", "#4", "OFF 960"],
        ["F11R", "+0.447", "−0.070", "#6", "−0.131", "#1", "OFF 960"],
        ["CLDN1", "+1.603", "+0.225", "loses (opp.)", "not scored", "GAP", "OFF 960"],
        ["CLDN5", "—", "not scored", "GAP", "not scored", "GAP", "OFF 960"],
        ["TJP1", "+0.593", "not scored solo", "GAP", "module only", "module", "OFF 960"],
    ]
    table = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.15, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#d6d3d1")
        if r == 0:
            cell.set_facecolor("#1c1917")
            cell.set_text_props(color="white", fontweight="bold")
        elif r == 1:
            cell.set_facecolor("#ffedd5")
            cell.set_text_props(fontweight="bold")
        elif r in (8, 9, 10):
            cell.set_facecolor("#f5f5f4")
            cell.set_text_props(color="#57534e")
        else:
            cell.set_facecolor("#fafaf9")
    fig.text(
        0.5, 0.04,
        "Win criterion = most negative concordant-4 KRT-partial %pos vs T/NK among scored genes. "
        "Locked unadj CLDN4 ρ=−0.531. Locked CosMx CLDN4 0.36/0.52 not replaced. Sources: #748 #644 #741 #698 #643.",
        ha="center", fontsize=7.5, color="#57534e",
    )
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"ppt_table_compact.{ext}", dpi=180, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    fig_why_wins()
    fig_compact_table()
    print("wrote", OUT)
