#!/usr/bin/env python3
"""Figure 5. Concordant-4 malignant CLDN4 percent-positive versus T/NK fraction.

Data visualization only. Every annotated estimate is read from the source TSVs.
Cohort confidence intervals are not in the tables and are not drawn.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

import pubstyle
from pubstyle import (
    COHORT_COLORS,
    COHORT_ORDER,
    PALETTE,
    apply_publication_style,
    fmt_p,
    fmt_rho,
    forest_diamond,
    forest_points,
    mm_to_in,
    style_ax,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
POOLED_COHORT = "GSE123902+GSE131907+GSE205335+GSE189357"

# Display order for the subordinate gene-set panel. Immune programs, then
# tight junction. Not sorted by the observed effect.
FAMILY_ORDER = ["IFN", "chemokine", "MHC-I/APM", "TJ"]
FAMILY_LABEL = {
    "IFN": "IFN",
    "chemokine": "Chemokine",
    "MHC-I/APM": "MHC-I/APM",
    "TJ": "Tight junction",
}

ARIMO = [
    "/usr/share/fonts/truetype/croscore/Arimo-Regular.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Bold.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Italic.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-BoldItalic.ttf",
]


def load_tsv(name: str) -> list[dict[str, str]]:
    path = DATA / name
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fnum(value: str | None) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    rx, ry = rankdata(x), rankdata(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def fmt_signed(value: float, digits: int = 2) -> str:
    sign = "−" if value < 0 else ""
    return f"{sign}{abs(value):.{digits}f}"


def use_arimo() -> None:
    """Arimo is the Arial-metric face installed here and contains the superscripts."""
    missing = [path for path in ARIMO if not Path(path).is_file()]
    require(not missing, "Arimo font files are missing: " + ", ".join(missing))
    for path in ARIMO:
        font_manager.fontManager.addfont(path)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arimo"]


def load_locked():
    units = load_tsv("tnk_units.tsv")
    pooled_rows = load_tsv("tnk_pooled.tsv")
    single_rows = load_tsv("tnk_singles.tsv")
    family_rows = load_tsv("family_de.tsv")

    require(len(units) == 65, f"expected 65 units, found {len(units)}")
    for row in units:
        require(fnum(row["cldn4_pct"]) is not None, "missing cldn4_pct")
        require(fnum(row["frac_tnk"]) is not None, "missing frac_tnk")
        require(row["cohort"] in COHORT_ORDER, f"unexpected cohort {row['cohort']}")

    pooled_hits = [
        row
        for row in pooled_rows
        if row["kind"] == "concordant4" and row["score"] == "pct" and row["k"] == "4"
    ]
    require(len(pooled_hits) == 1, "expected one concordant-4 pct pooled row")
    pooled = pooled_hits[0]
    require(int(pooled["n"]) == 65, "pooled n is not 65")
    require(float(pooled["I2"]) == 0.0, "pooled I² is not 0")
    require(pooled["combo"] == POOLED_COHORT, "pooled combo does not match the four cohorts")

    singles = {
        row["combo"]: row
        for row in single_rows
        if row["kind"] == "single" and row["score"] == "pct"
    }
    require(set(singles) == set(COHORT_ORDER), "pct single-cohort rows do not match the four cohorts")
    require(sum(int(singles[c]["n"]) for c in COHORT_ORDER) == 65, "cohort n does not sum to 65")

    for cohort in COHORT_ORDER:
        rows = [row for row in units if row["cohort"] == cohort]
        require(len(rows) == int(singles[cohort]["n"]), f"{cohort} unit count != singles n")
        rho = spearman(
            [float(row["cldn4_pct"]) for row in rows],
            [float(row["frac_tnk"]) for row in rows],
        )
        require(
            abs(rho - float(singles[cohort]["rho"])) < 1e-9,
            f"{cohort} Spearman in tnk_singles does not match tnk_units",
        )

    q1 = [float(row["frac_tnk"]) for row in units if row["quartile"] == "Q1"]
    q4 = [float(row["frac_tnk"]) for row in units if row["quartile"] == "Q4"]
    require(len(q1) == int(pooled["n_q1"]) == 19, "Q1 n is not 19")
    require(len(q4) == int(pooled["n_q4"]) == 16, "Q4 n is not 16")
    require(abs(median(q1) - float(pooled["median_q1"])) < 1e-12, "Q1 median mismatch")
    require(abs(median(q4) - float(pooled["median_q4"])) < 1e-12, "Q4 median mismatch")
    require(abs((median(q4) - median(q1)) - float(pooled["delta_median"])) < 1e-12, "delta median mismatch")

    families = []
    for name in FAMILY_ORDER:
        hits = [
            row
            for row in family_rows
            if row["family"] == name
            and row["split"] == "q4q1"
            and row["cohort"] == POOLED_COHORT
        ]
        require(len(hits) == 1, f"missing combined q4q1 row for {name}")
        families.append(hits[0])
    require(int(float(families[0]["n"])) == 34, "family-panel n is not 34")
    require(int(float(families[0]["n_q1"])) == 18, "family-panel Q1 n is not 18")
    require(int(float(families[0]["n_q4"])) == 16, "family-panel Q4 n is not 16")

    # The count-matrix filter drops one Q1 unit. Do not relabel that n as 65.
    missing_counts = [row for row in units if row["in_count_matrix"] != "True"]
    require(len(missing_counts) == 1, "expected exactly one unit without a count matrix")
    require(missing_counts[0]["quartile"] == "Q1", "the unit without counts is not Q1")

    rho = float(pooled["rho"])
    p_value = float(pooled["p"])
    lo = float(pooled["ci95_lo"])
    hi = float(pooled["ci95_hi"])
    require(np.isfinite(lo) and np.isfinite(hi) and lo < rho < hi, "pooled 95% CI is unusable")
    require(fmt_rho(rho) == "−0.53", f"unexpected pooled rho display: {fmt_rho(rho)}")
    require(fmt_p(p_value) == "P = 1.6 × 10⁻⁵", f"unexpected pooled P display: {fmt_p(p_value)}")

    return {
        "units": units,
        "pooled": pooled,
        "singles": singles,
        "families": families,
        "q1": q1,
        "q4": q4,
    }


def cohort_n(singles: dict[str, dict[str, str]], cohort: str) -> int:
    return int(singles[cohort]["n"])


def marker_area_size(n: int) -> float:
    """Marker diameter in points. Area scales with unit n; this is not a meta-analytic weight."""
    return 1.48 * float(np.sqrt(n))


def style_panel(ax) -> None:
    style_ax(ax)
    ax.set_axisbelow(True)


def draw_scatter(ax, units, pooled) -> None:
    for cohort in COHORT_ORDER:
        rows = [row for row in units if row["cohort"] == cohort]
        ax.scatter(
            [float(row["cldn4_pct"]) for row in rows],
            [float(row["frac_tnk"]) for row in rows],
            s=22,
            c=COHORT_COLORS[cohort],
            edgecolors=PALETTE["neutral_black"],
            linewidths=0.3,
            zorder=3,
            clip_on=True,
        )

    x = np.array([float(row["cldn4_pct"]) for row in units], dtype=float)
    y = np.array([float(row["frac_tnk"]) for row in units], dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(float(x.min()), float(x.max()), 50)
    ax.plot(
        x_line,
        slope * x_line + intercept,
        color=PALETTE["neutral_mid"],
        lw=0.9,
        zorder=2,
        solid_capstyle="round",
    )

    ax.set_title(
        f"Pooled ρ = {fmt_rho(float(pooled['rho']))}      {fmt_p(float(pooled['p']))}\n"
        f"n = {int(pooled['n'])} units      I² = 0%",
        loc="left",
        fontsize=7,
        color=PALETTE["neutral_black"],
        pad=8,
    )

    ax.set_xlim(-4, 104)
    ax.set_ylim(-0.03, 0.98)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel("Malignant CLDN4⁺ cells (%)")
    ax.set_ylabel("T/NK fraction")
    style_panel(ax)


def draw_quartile(ax, units, pooled) -> None:
    q1 = [float(row["frac_tnk"]) for row in units if row["quartile"] == "Q1"]
    q4 = [float(row["frac_tnk"]) for row in units if row["quartile"] == "Q4"]
    box = ax.boxplot(
        [q1, q4],
        positions=[0, 1],
        widths=0.52,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": PALETTE["neutral_black"], "linewidth": 1.15},
        whiskerprops={"color": PALETTE["neutral_dark"], "linewidth": 0.7},
        capprops={"color": PALETTE["neutral_dark"], "linewidth": 0.7},
        boxprops={
            "facecolor": "white",
            "edgecolor": PALETTE["neutral_black"],
            "linewidth": 0.7,
        },
        zorder=1,
    )
    for patch in box["boxes"]:
        patch.set_zorder(1)

    rng = np.random.default_rng(21)
    for cohort in COHORT_ORDER:
        for xpos, quartile in ((0, "Q1"), (1, "Q4")):
            values = [
                float(row["frac_tnk"])
                for row in units
                if row["cohort"] == cohort and row["quartile"] == quartile
            ]
            if not values:
                continue
            jitter = rng.uniform(-0.13, 0.13, size=len(values))
            ax.scatter(
                np.full(len(values), xpos) + jitter,
                values,
                s=16,
                c=COHORT_COLORS[cohort],
                edgecolors=PALETTE["neutral_black"],
                linewidths=0.3,
                zorder=3,
            )

    med_q1 = float(pooled["median_q1"])
    med_q4 = float(pooled["median_q4"])
    ax.text(0.34, med_q1, fmt_signed(med_q1), va="center", ha="left", fontsize=6, color=PALETTE["neutral_dark"])
    ax.text(1.34, med_q4, fmt_signed(med_q4), va="center", ha="left", fontsize=6, color=PALETTE["neutral_dark"])

    bracket_y, bracket_h = 0.955, 0.02
    ax.plot(
        [0, 0, 1, 1],
        [bracket_y, bracket_y + bracket_h, bracket_y + bracket_h, bracket_y],
        color=PALETTE["neutral_black"],
        lw=0.6,
        clip_on=False,
        zorder=2,
    )
    annotation = (
        f"Δ median = {fmt_signed(float(pooled['delta_median']))}\n"
        f"rank-biserial r = {fmt_signed(float(pooled['r_rb']))}\n"
        f"{fmt_p(float(pooled['p_q4q1']))}"
    )
    ax.text(
        0.5,
        bracket_y + bracket_h + 0.012,
        annotation,
        ha="center",
        va="bottom",
        fontsize=6,
        color=PALETTE["neutral_black"],
        linespacing=1.2,
        clip_on=False,
    )

    ax.set_xlim(-0.55, 1.62)
    ax.set_ylim(0, 1.28)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Q1\nlow CLDN4\nn = {int(pooled['n_q1'])}", f"Q4\nhigh CLDN4\nn = {int(pooled['n_q4'])}"])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylabel("T/NK fraction")
    ax.set_xlabel("Within-cohort CLDN4 quartile")
    style_panel(ax)


def draw_forest(ax, ax_table, singles, pooled) -> None:
    cohort_y = {cohort: float(len(COHORT_ORDER) - 1 - i) for i, cohort in enumerate(COHORT_ORDER)}
    pool_y = -1.15

    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.6, zorder=0)
    ax.axhline(-0.45, color=PALETTE["neutral_light"], lw=0.6, zorder=0)
    ax_table.axhline(-0.45, color=PALETTE["neutral_light"], lw=0.6, zorder=0)

    estimates, ys, colors, sizes, ns = [], [], [], [], []
    for cohort in COHORT_ORDER:
        row = singles[cohort]
        y = cohort_y[cohort]
        n = int(row["n"])
        estimates.append(float(row["rho"]))
        ys.append(y)
        colors.append(COHORT_COLORS[cohort])
        sizes.append(marker_area_size(n))
        ns.append(n)
        ax_table.text(0.02, y, fmt_rho(float(row["rho"])), ha="left", va="center", fontsize=6, color=PALETTE["neutral_black"])
        ax_table.text(0.40, y, fmt_p(float(row["p"])), ha="left", va="center", fontsize=6, color=PALETTE["neutral_black"])

    for y, est, n, color, size in zip(ys, estimates, ns, colors, sizes):
        forest_points(ax, [y], [est], [n], [color], markersize=size)

    rho = float(pooled["rho"])
    lo = float(pooled["ci95_lo"])
    hi = float(pooled["ci95_hi"])
    forest_diamond(ax, pool_y, rho, lo, hi, PALETTE["neutral_black"], height=0.84)
    ax.text(
        rho,
        pool_y - 1.22,
        f"I² = 0%\n95% CI {fmt_signed(lo)} to {fmt_signed(hi)}",
        ha="center",
        va="top",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        linespacing=1.3,
    )

    ax_table.text(0.02, pool_y, fmt_rho(rho), ha="left", va="center", fontsize=6, fontweight="bold", color=PALETTE["neutral_black"])
    ax_table.text(
        0.40,
        pool_y,
        fmt_p(float(pooled["p"])),
        ha="left",
        va="center",
        fontsize=6,
        fontweight="bold",
        color=PALETTE["neutral_black"],
    )
    header_y = max(cohort_y.values()) + 0.72
    ax_table.text(0.02, header_y, "ρ", ha="left", va="center", fontsize=6, fontweight="bold", color=PALETTE["neutral_dark"])
    ax_table.text(0.40, header_y, "P", ha="left", va="center", fontsize=6, fontweight="bold", color=PALETTE["neutral_dark"])

    yticks = [cohort_y[c] for c in COHORT_ORDER] + [pool_y]
    ylabels = [f"{c}  ({cohort_n(singles, c)})" for c in COHORT_ORDER] + [f"Pooled  ({int(pooled['n'])})"]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    for label in ax.get_yticklabels():
        label.set_color(PALETTE["neutral_black"])
    ax.get_yticklabels()[-1].set_fontweight("bold")

    ax.set_xlim(-1.12, 0.18)
    ax.set_ylim(pool_y - 2.35, header_y + 0.42)
    ax.xaxis.labelpad = 3
    ax.set_xticks([-1, -0.5, 0])
    ax.set_xticklabels(["−1", "−0.5", "0"])
    ax.set_xlabel("Spearman ρ")
    style_panel(ax)

    ax_table.set_ylim(ax.get_ylim())
    ax_table.set_xlim(0, 1)
    ax_table.axis("off")


def draw_families(ax, families) -> None:
    ys = np.arange(len(families))[::-1]
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.6, zorder=0)
    for y, row in zip(ys, families):
        estimate = float(row["logFC"])
        se = float(row["se"])
        require(np.isfinite(estimate) and np.isfinite(se) and se > 0, f"bad SE for {row['family']}")
        ax.errorbar(
            estimate,
            y,
            xerr=se,
            fmt="o",
            color=PALETTE["blue_main"],
            markersize=4.4,
            markeredgecolor="white",
            markeredgewidth=0.35,
            elinewidth=0.85,
            capsize=1.8,
            capthick=0.7,
            zorder=3,
        )
        ax.text(
            0.78,
            y,
            fmt_p(float(row["p"])),
            ha="left",
            va="center",
            fontsize=6,
            color=PALETTE["neutral_black"],
        )

    n_units = int(float(families[0]["n"]))
    n_q1 = int(float(families[0]["n_q1"]))
    n_q4 = int(float(families[0]["n_q4"]))
    ax.text(
        0.0,
        1.08,
        f"n = {n_units} units ({n_q1} Q1, {n_q4} Q4)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        clip_on=False,
    )

    ax.set_yticks(ys)
    ax.set_yticklabels([FAMILY_LABEL[row["family"]] for row in families])
    ax.set_ylim(-0.65, len(families) - 0.25)
    ax.set_xlim(-1.55, 2.15)
    ax.set_xticks([-1.5, -1.0, -0.5, 0, 0.5])
    ax.set_xticklabels(["−1.5", "−1", "−0.5", "0", "0.5"])
    ax.set_xlabel("Mean logFC ± SE (Q4 versus Q1)")
    style_panel(ax)


def cohort_handles(singles) -> list[Line2D]:
    handles = []
    for cohort in COHORT_ORDER:
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="None",
                markersize=4.6,
                markerfacecolor=COHORT_COLORS[cohort],
                markeredgecolor=PALETTE["neutral_black"],
                markeredgewidth=0.3,
                label=f"{cohort} (n = {cohort_n(singles, cohort)})",
            )
        )
    handles.append(
        Line2D(
            [],
            [],
            color=PALETTE["neutral_mid"],
            lw=0.9,
            label="Linear guide",
        )
    )
    return handles


def place_letters(fig, axes) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, letter in zip(axes, "abcd"):
        boxes = [label.get_window_extent(renderer) for label in ax.get_yticklabels() if label.get_text()]
        ylabel = ax.yaxis.label.get_window_extent(renderer)
        left = min([box.x0 for box in boxes] + [ylabel.x0, ax.get_window_extent(renderer).x0])
        left_fig = left / fig.bbox.width
        top = ax.get_position().y1
        fig.text(
            max(0.004, left_fig - 0.028),
            top + 0.012,
            letter,
            fontsize=8,
            fontweight="bold",
            ha="left",
            va="bottom",
            color=PALETTE["neutral_black"],
        )


def build() -> plt.Figure:
    locked = load_locked()
    fig = plt.figure(figsize=(mm_to_in(183), mm_to_in(118)))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.48, 1.0],
        height_ratios=[1.22, 1.0],
        left=0.168,
        right=0.988,
        top=0.90,
        bottom=0.145,
        wspace=0.50,
        hspace=0.52,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    forest_grid = grid[1, 0].subgridspec(1, 2, width_ratios=[1.55, 0.78], wspace=0.06)
    ax_c = fig.add_subplot(forest_grid[0, 0])
    ax_c_table = fig.add_subplot(forest_grid[0, 1], sharey=ax_c)
    ax_d = fig.add_subplot(grid[1, 1])

    draw_scatter(ax_a, locked["units"], locked["pooled"])
    draw_quartile(ax_b, locked["units"], locked["pooled"])
    draw_forest(ax_c, ax_c_table, locked["singles"], locked["pooled"])
    draw_families(ax_d, locked["families"])

    fig.legend(
        handles=cohort_handles(locked["singles"]),
        loc="lower center",
        bbox_to_anchor=(0.54, 0.012),
        ncol=5,
        frameon=False,
        fontsize=6,
        handletextpad=0.35,
        columnspacing=0.9,
        borderaxespad=0.0,
    )
    place_letters(fig, [ax_a, ax_b, ax_c, ax_d])
    return fig


def save(fig: plt.Figure) -> None:
    for ext in ("pdf", "svg"):
        fig.savefig(ROOT / f"Fig5.{ext}", bbox_inches="tight", pad_inches=0.02, facecolor="white")
    fig.savefig(ROOT / "Fig5.png", dpi=600, bbox_inches="tight", pad_inches=0.02, facecolor="white")


def main() -> None:
    apply_publication_style()
    use_arimo()
    fig = build()
    save(fig)
    print("Wrote Fig5.pdf, Fig5.svg, Fig5.png")


if __name__ == "__main__":
    main()
