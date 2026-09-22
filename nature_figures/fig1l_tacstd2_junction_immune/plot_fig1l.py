#!/usr/bin/env python3
"""Fig. 1l. Three separate concordant-4 relations. Data visualization only.

Every drawn estimate is recomputed from source_data/units.tsv and checked
against source_data/effects.tsv. Cohort intervals use the same Fisher-z
formula as the pool. The Part 2 CLDN4 percent-positive analysis is not drawn.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from compute_scores import cohort_ci, fisher_z_dl, spearman
from pubstyle import (
    COHORT_COLORS,
    COHORT_ORDER,
    PALETTE,
    apply_publication_style,
    fmt_p,
    fmt_rho,
    mm_to_in,
    style_ax,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"

for face in (
    "/usr/share/fonts/truetype/croscore/Arimo-Regular.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Bold.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Italic.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-BoldItalic.ttf",
):
    font_manager.fontManager.addfont(face)

PRIMARY = [
    ("tacstd2_vs_junction", "junction_hallmark", "tacstd2", "junction_hallmark", "Tumour TACSTD2\nvs junction score"),
    ("tacstd2_vs_tnk", "junction_hallmark", "tacstd2", "frac_tnk", "Tumour TACSTD2\nvs T/NK fraction"),
    ("junction_vs_tnk", "junction_hallmark", "junction_hallmark", "frac_tnk", "Junction score\nvs T/NK fraction"),
]

Y_COHORT = np.array([4.0, 3.0, 2.0, 1.0])
Y_POOL = 0.0
UNIT_WORD = {
    "GSE123902": "donors",
    "GSE131907": "samples",
    "GSE205335": "patients",
    "GSE189357": "patients",
}


def load_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / name, sep="\t")


def require_close(got: float, expected: float, label: str) -> None:
    if not np.isfinite(got) or abs(got - expected) > 1e-8:
        raise SystemExit(f"{label}: recomputed {got} != table {expected}")


def relation_table(units: pd.DataFrame, xcol: str, ycol: str) -> pd.DataFrame:
    rows = []
    rhos, ns = [], []
    for cohort in COHORT_ORDER:
        m = units.loc[units["cohort"] == cohort]
        rho, p, n = spearman(m[xcol], m[ycol])
        lo, hi = cohort_ci(rho, n)
        rows.append({"cohort": cohort, "n": n, "rho": rho, "p": p, "ci95_lo": lo, "ci95_hi": hi, "kind": "cohort"})
        rhos.append(rho)
        ns.append(n)
    pool = fisher_z_dl(rhos, ns)
    rows.append(
        {
            "cohort": "Pooled",
            "n": pool["n"],
            "rho": pool["rho"],
            "p": pool["p"],
            "ci95_lo": pool["ci95_lo"],
            "ci95_hi": pool["ci95_hi"],
            "I2": pool["I2"],
            "kind": "pooled_dl_fisher_z",
        }
    )
    return pd.DataFrame(rows)


def check_against_effects(fresh: pd.DataFrame, effects: pd.DataFrame, relation: str, score: str) -> None:
    stored = effects.loc[(effects["relation"] == relation) & (effects["score"] == score)]
    for _, row in fresh.iterrows():
        if row["kind"] == "cohort":
            hit = stored.loc[(stored["kind"] == "cohort") & (stored["cohort"] == row["cohort"])].iloc[0]
        else:
            hit = stored.loc[stored["kind"] == "pooled_dl_fisher_z"].iloc[0]
        for key in ("n", "rho", "p", "ci95_lo", "ci95_hi"):
            require_close(float(row[key]), float(hit[key]), f"{relation} {score} {row['cohort']} {key}")
        if row["kind"] == "pooled_dl_fisher_z":
            require_close(float(row["I2"]), float(hit["I2"]), f"{relation} {score} I2")


def forest(ax, table: pd.DataFrame, title: str, show_labels: bool) -> None:
    ax.axvline(0, color=PALETTE["neutral_light"], lw=0.6, zorder=0)
    for y, cohort in zip(Y_COHORT, COHORT_ORDER):
        row = table.loc[table["cohort"] == cohort].iloc[0]
        color = COHORT_COLORS[cohort]
        ax.plot(
            [row["ci95_lo"], row["ci95_hi"]],
            [y, y],
            color=color,
            lw=1.05,
            solid_capstyle="round",
            zorder=2,
        )
        ax.plot(
            row["rho"],
            y,
            "o",
            color=color,
            markersize=4.4,
            markeredgecolor="white",
            markeredgewidth=0.4,
            zorder=3,
        )
    pool = table.loc[table["kind"] == "pooled_dl_fisher_z"].iloc[0]
    diamond_x = [pool["ci95_lo"], pool["rho"], pool["ci95_hi"], pool["rho"]]
    diamond_y = [Y_POOL, Y_POOL + 0.16, Y_POOL, Y_POOL - 0.16]
    ax.fill(
        diamond_x,
        diamond_y,
        facecolor=PALETTE["neutral_black"],
        edgecolor=PALETTE["neutral_black"],
        lw=0.3,
        zorder=3,
        closed=True,
    )
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylim(-0.72, 4.55)
    ax.set_yticks(list(Y_COHORT) + [Y_POOL])
    if show_labels:
        labels = []
        for cohort in COHORT_ORDER:
            n = int(table.loc[table["cohort"] == cohort, "n"].iloc[0])
            labels.append(f"{cohort}\n{n} {UNIT_WORD[cohort]}")
        labels.append(f"Pooled\n{int(pool['n'])} units")
        ax.set_yticklabels(labels)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Spearman ρ")
    ax.set_title(title, loc="left", pad=6, color=PALETTE["neutral_black"], fontsize=7)
    style_ax(ax)
    # Numeric ρ sits in the right margin so it does not cross the interval.
    for y, cohort in list(zip(Y_COHORT, COHORT_ORDER)) + [(Y_POOL, "Pooled")]:
        row = table.loc[table["cohort"] == cohort].iloc[0]
        weight = "bold" if cohort == "Pooled" else "regular"
        extra = ""
        if cohort == "Pooled":
            i2 = float(row["I2"])
            i2_txt = "0" if i2 < 0.05 else f"{i2:.0f}"
            extra = f"\nI² = {i2_txt}%"
        ax.text(
            1.06,
            y,
            fmt_rho(float(row["rho"])) + extra,
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=6,
            fontweight=weight,
            color=PALETTE["neutral_black"],
            linespacing=1.05,
            clip_on=False,
        )


def panel_letter(fig, ax, letter: str, xshift: float = -0.012, yshift: float = 0.012) -> None:
    pos = ax.get_position()
    fig.text(
        pos.x0 + xshift,
        pos.y1 + yshift,
        letter,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=PALETTE["neutral_black"],
    )


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf", "svg"):
        kw = {"dpi": 600} if ext == "png" else {}
        fig.savefig(ROOT / f"{stem}.{ext}", bbox_inches="tight", pad_inches=0.04, facecolor="white", **kw)


def plot_main(tables: list[pd.DataFrame]) -> None:
    apply_publication_style()
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(mm_to_in(183), mm_to_in(78)),
        gridspec_kw={"wspace": 0.55, "left": 0.14, "right": 0.90, "top": 0.78, "bottom": 0.16},
    )
    letters = ["a", "b", "c"]
    for i, (ax, table, spec) in enumerate(zip(axes, tables, PRIMARY)):
        forest(ax, table, spec[4], show_labels=(i == 0))
        panel_letter(fig, ax, letters[i])
    save(fig, "Fig1l")
    plt.close(fig)


def plot_sensitivity(tables: dict[tuple[str, str], pd.DataFrame]) -> None:
    apply_publication_style()
    specs = [
        ("a", "junction_hallmark_no_cldn4", "tacstd2_vs_junction", "Hallmark, CLDN4 removed\nTACSTD2 vs junction score"),
        ("b", "junction_hallmark_no_cldn4", "junction_vs_tnk", "Hallmark, CLDN4 removed\nJunction score vs T/NK"),
        ("c", "junction_kegg", "tacstd2_vs_junction", "KEGG tight junction\nTACSTD2 vs junction score"),
        ("d", "junction_kegg", "junction_vs_tnk", "KEGG tight junction\nJunction score vs T/NK"),
    ]
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(mm_to_in(183), mm_to_in(132)),
        gridspec_kw={"wspace": 0.62, "hspace": 0.55, "left": 0.14, "right": 0.90, "top": 0.90, "bottom": 0.08},
    )
    for ax, (letter, score, relation, title) in zip(axes.ravel(), specs):
        forest(ax, tables[(relation, score)], title, show_labels=(letter in {"a", "c"}))
        panel_letter(fig, ax, letter)
    save(fig, "Fig1l_ED_sensitivity")
    plt.close(fig)


def plot_scatters(units: pd.DataFrame, tables: list[pd.DataFrame]) -> None:
    apply_publication_style()
    rows = [
        ("tacstd2", "junction_hallmark", "Tumour TACSTD2", "Junction score", "log2 (TMM-CPM+1)", "mean log2 (TMM-CPM+1)"),
        ("tacstd2", "frac_tnk", "Tumour TACSTD2", "T/NK fraction", "log2 (TMM-CPM+1)", "fraction of cells"),
        ("junction_hallmark", "frac_tnk", "Junction score", "T/NK fraction", "mean log2 (TMM-CPM+1)", "fraction of cells"),
    ]
    fig, axes = plt.subplots(
        3,
        4,
        figsize=(mm_to_in(183), mm_to_in(142)),
        gridspec_kw={"wspace": 0.42, "hspace": 0.72, "left": 0.10, "right": 0.99, "top": 0.90, "bottom": 0.07},
    )
    for r, (xcol, ycol, xlab, ylab, xunit, yunit) in enumerate(rows):
        xmin = float(units[xcol].min())
        xmax = float(units[xcol].max())
        ymin = float(units[ycol].min())
        ymax = float(units[ycol].max())
        xpad = 0.06 * (xmax - xmin)
        ypad = 0.08 * (ymax - ymin)
        for c, cohort in enumerate(COHORT_ORDER):
            ax = axes[r, c]
            m = units.loc[units["cohort"] == cohort]
            ax.scatter(
                m[xcol],
                m[ycol],
                s=16,
                c=COHORT_COLORS[cohort],
                edgecolors="white",
                linewidths=0.3,
                zorder=3,
            )
            ax.set_xlim(xmin - xpad, xmax + xpad)
            ax.set_ylim(ymin - ypad, ymax + ypad)
            style_ax(ax)
            row = tables[r].loc[tables[r]["cohort"] == cohort].iloc[0]
            stats_line = f"{fmt_rho(float(row['rho']))}   {fmt_p(float(row['p']))}   n = {int(row['n'])}"
            if r == 0:
                ax.set_title(
                    f"{cohort}\n{stats_line}",
                    fontsize=6.5,
                    color=PALETTE["neutral_black"],
                    pad=2,
                    loc="left",
                )
            else:
                ax.set_title(stats_line, fontsize=6, color=PALETTE["neutral_dark"], pad=2, loc="left")
            if c == 0:
                ax.set_ylabel(f"{ylab}\n{yunit}")
            else:
                ax.set_ylabel("")
            if r == 2:
                ax.set_xlabel(f"{xlab}\n{xunit}")
            else:
                ax.set_xlabel("")
            ax.tick_params(axis="x", labelsize=5.5)
            ax.tick_params(axis="y", labelsize=5.5)
        # Row letter on the first panel of each relation.
        panel_letter(fig, axes[r, 0], "abc"[r], xshift=-0.055, yshift=0.028)
    save(fig, "Fig1l_ED_scatters")
    plt.close(fig)


def main() -> None:
    units = load_tsv("units.tsv")
    units = units.loc[units["included"] == True].copy()  # noqa: E712
    if len(units) != 64:
        raise SystemExit(f"expected 64 included units, found {len(units)}")
    effects = load_tsv("effects.tsv")
    tables = []
    cache: dict[tuple[str, str], pd.DataFrame] = {}
    needed = set()
    for relation, score, xcol, ycol, _title in PRIMARY:
        needed.add((relation, score, xcol, ycol))
    needed.update(
        [
            ("tacstd2_vs_junction", "junction_hallmark_no_cldn4", "tacstd2", "junction_hallmark_no_cldn4"),
            ("junction_vs_tnk", "junction_hallmark_no_cldn4", "junction_hallmark_no_cldn4", "frac_tnk"),
            ("tacstd2_vs_junction", "junction_kegg", "tacstd2", "junction_kegg"),
            ("junction_vs_tnk", "junction_kegg", "junction_kegg", "frac_tnk"),
        ]
    )
    for relation, score, xcol, ycol in needed:
        fresh = relation_table(units, xcol, ycol)
        check_against_effects(fresh, effects, relation, score)
        cache[(relation, score)] = fresh
    for relation, score, _x, _y, _title in PRIMARY:
        tables.append(cache[(relation, score)])
    plot_main(tables)
    plot_sensitivity(cache)
    plot_scatters(units, tables)


if __name__ == "__main__":
    main()
