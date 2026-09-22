#!/usr/bin/env python3
"""Nature-style figure: barrier/inhibitory ligand mean Δ and family DE.

Plots tabulated values only (family_barrier_inhibitory.tsv, family_de.tsv).
Error bars are ± 1 SE from the source `se` column. No intervals are invented.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

for _font in (
    "/usr/share/fonts/truetype/croscore/Arimo-Regular.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Bold.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Italic.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-BoldItalic.ttf",
):
    font_manager.fontManager.addfont(_font)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_3": "#8BCF8B",
    "red_strong": "#B64342",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "band": "#F4F5F6",
}

LIGAND_COLORS = {
    "F11R": PALETTE["blue_main"],
    "NECTIN2-TIGIT": PALETTE["teal"],
    "CDH1": PALETTE["violet"],
    "LGALS9": PALETTE["red_strong"],
    "family": PALETTE["neutral_black"],
}

COHORT_COLORS = {
    "GSE123902": PALETTE["blue_main"],
    "GSE131907": PALETTE["teal"],
    "GSE205335": PALETTE["violet"],
    "GSE189357": PALETTE["green_3"],
}
COHORT_ORDER = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
POOLED_KEY = "GSE123902+GSE131907+GSE205335+GSE189357"
POOLED_COLOR = PALETTE["neutral_black"]

# Display order: strongest negative pooled Q4-vs-Q1 logFC at the top.
FAMILY_ORDER = ["chemokine", "MHC-I/APM", "IFN", "TJ"]

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arimo", "Arial", "Liberation Sans", "Noto Sans", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
            "pdf.compression": 9,
            "axes.unicode_minus": False,
        }
    )


def mm(v: float) -> float:
    return v / 25.4


def unicode_sup(n: int) -> str:
    return str(int(n)).translate(_SUPERSCRIPT)


def fmt_p(p: float, digits: int = 1) -> str:
    """Numeric P/FDR text. Does not prepend a label."""
    p = float(p)
    if p <= 0:
        return "< 10⁻³⁰⁰"
    if p >= 0.001:
        if p >= 0.1:
            return f"{p:.2f}"
        return f"{p:.3f}"
    exp = int(np.floor(np.log10(p)))
    mant = p / (10**exp)
    return f"{mant:.{digits}f} × 10{unicode_sup(exp)}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fnum(row: dict[str, str], key: str) -> float:
    raw = row[key].strip()
    if raw == "":
        raise ValueError(f"Missing {key} in row {row}")
    return float(raw)


def style_ax(ax) -> None:
    ax.tick_params(axis="both", which="major", labelsize=6, length=2.4, width=0.8, pad=1.6)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
        ax.spines[spine].set_color(PALETTE["neutral_black"])
    ax.xaxis.label.set_color(PALETTE["neutral_black"])
    ax.yaxis.label.set_color(PALETTE["neutral_black"])
    ax.tick_params(colors=PALETTE["neutral_black"])
    ax.title.set_color(PALETTE["neutral_black"])


def pair_label(row: dict[str, str]) -> str:
    ligand = row["ligand"]
    receptor = row["receptor"].replace("_", "/")
    return f"{ligand}–{receptor}"


def load_barrier() -> tuple[list[dict[str, str]], dict[str, str]]:
    rows = read_tsv(DATA / "family_barrier_inhibitory.tsv")
    pairs = [r for r in rows if r["row"] == "pair"]
    aggs = [r for r in rows if r["row"] == "FAMILY_AGGREGATE"]
    if len(pairs) != 7 or len(aggs) != 1:
        raise SystemExit(f"Unexpected barrier rows: {len(pairs)} pairs, {len(aggs)} aggregates")
    for r in pairs + aggs:
        if r["split"] != "q4q1":
            raise SystemExit(f"Unexpected split {r['split']}")
        n_sum = sum(int(r[k]) for k in ("n_gse123902", "n_gse131907", "n_gse205335", "n_gse189357"))
        if n_sum != int(r["n_units"]):
            raise SystemExit(f"Cohort n does not sum to n_units for {r['interaction_name']}")
    pairs.sort(key=lambda r: fnum(r, "mean_delta"), reverse=True)
    return pairs, aggs[0]


def load_de() -> dict[tuple[str, str, str], dict[str, str]]:
    rows = read_tsv(DATA / "family_de.tsv")
    out: dict[tuple[str, str, str], dict[str, str]] = {}
    for r in rows:
        key = (r["family"], r["split"], r["cohort"])
        if key in out:
            raise SystemExit(f"Duplicate DE row {key}")
        out[key] = r
    return out


def draw_mean_se(ax, x: float, y: float, se: float, color: str, marker: str, filled: bool, ms: float) -> None:
    ax.plot([x - se, x + se], [y, y], color=color, lw=0.95, solid_capstyle="butt", zorder=2, clip_on=False)
    cap = 0.11
    ax.plot([x - se, x - se], [y - cap, y + cap], color=color, lw=0.75, zorder=2, clip_on=False)
    ax.plot([x + se, x + se], [y - cap, y + cap], color=color, lw=0.75, zorder=2, clip_on=False)
    ax.plot(
        x,
        y,
        marker,
        markersize=ms,
        color=color,
        markerfacecolor=color if filled else "white",
        markeredgecolor=color,
        markeredgewidth=0.7,
        zorder=4,
        clip_on=False,
    )


def draw_median_tick(ax, x: float, y: float, color: str) -> None:
    ax.plot([x, x], [y - 0.16, y + 0.16], color=color, lw=1.15, solid_capstyle="butt", zorder=3, clip_on=False)


def panel_letter(ax, letter: str, x: float, y: float) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=PALETTE["neutral_black"],
        clip_on=False,
    )


def plot_barrier(ax, pairs: list[dict[str, str]], agg: dict[str, str]) -> None:
    # Top = largest mean Δ. Aggregate is the summary row, separated by a short gap.
    y_pairs = [1.25 + i for i in range(len(pairs) - 1, -1, -1)]
    y_agg = 0.0

    for y in y_pairs:
        ax.axhspan(y - 0.44, y + 0.44, color=PALETTE["band"], zorder=0, linewidth=0)
    ax.axhspan(y_agg - 0.44, y_agg + 0.44, color=PALETTE["band"], zorder=0, linewidth=0)
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.6, zorder=1, linestyle=(0, (2.2, 1.6)))
    ax.axhline(0.62, color=PALETTE["neutral_light"], lw=0.6, zorder=1)

    for y, row in zip(y_pairs, pairs):
        color = LIGAND_COLORS[row["axis"]]
        mean = fnum(row, "mean_delta")
        se = fnum(row, "se")
        med = fnum(row, "median_delta")
        draw_mean_se(ax, mean, y, se, color, "o", True, 4.6)
        draw_median_tick(ax, med, y, color)
        _label_barrier_row(ax, y, pair_label(row), row, color, bold=False)

    mean = fnum(agg, "mean_delta")
    se = fnum(agg, "se")
    med = fnum(agg, "median_delta")
    draw_mean_se(ax, mean, y_agg, se, LIGAND_COLORS["family"], "D", True, 5.0)
    draw_median_tick(ax, med, y_agg, LIGAND_COLORS["family"])
    _label_barrier_row(ax, y_agg, "Family aggregate", agg, LIGAND_COLORS["family"], bold=True)

    # Meta P sits to the right of the longest SE bar (F11R mean+SE = 0.211).
    p_x = 0.238
    for y, row in zip(y_pairs, pairs):
        ax.text(
            p_x,
            y,
            fmt_p(fnum(row, "meta_p")),
            ha="left",
            va="center",
            fontsize=6,
            color=PALETTE["neutral_black"],
            clip_on=False,
        )
    ax.text(
        p_x,
        y_agg,
        fmt_p(fnum(agg, "meta_p")),
        ha="left",
        va="center",
        fontsize=6,
        fontweight="bold",
        color=PALETTE["neutral_black"],
        clip_on=False,
    )
    ax.text(
        p_x,
        y_pairs[0] + 0.50,
        "Meta P",
        ha="left",
        va="center",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        clip_on=False,
    )

    ax.set_xlim(-0.012, 0.385)
    ax.set_ylim(-0.72, y_pairs[0] + 1.05)
    ax.set_yticks([])
    ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    ax.set_xticklabels(["0", "0.05", "0.10", "0.15", "0.20"])
    ax.set_xlabel("Mean Δ (Q4 − Q1)")
    ax.set_title("Barrier and inhibitory ligand–receptor pairs", loc="left", pad=10, fontsize=7)
    style_ax(ax)

    mean_handle = mpl.lines.Line2D(
        [0],
        [0],
        marker="o",
        color=PALETTE["neutral_dark"],
        markerfacecolor=PALETTE["neutral_dark"],
        markeredgecolor=PALETTE["neutral_dark"],
        markersize=4.2,
        lw=0.95,
        label="Mean ± 1 SE",
    )
    med_handle = mpl.lines.Line2D(
        [0],
        [0],
        marker="|",
        color=PALETTE["neutral_dark"],
        markeredgewidth=1.15,
        markersize=8,
        lw=0,
        label="Median",
    )
    ax.legend(
        handles=[mean_handle, med_handle],
        loc="lower right",
        bbox_to_anchor=(1.0, 1.005),
        ncol=2,
        frameon=False,
        fontsize=6,
        handlelength=1.5,
        columnspacing=1.0,
        borderaxespad=0.0,
    )


def _label_barrier_row(ax, y: float, name: str, row: dict[str, str], color: str, bold: bool) -> None:
    n = int(row["n_units"])
    n_pos = int(row["n_pos"])
    n_neg = int(row["n_neg"])
    weight = "bold" if bold else "regular"
    ax.text(
        -0.018,
        y + 0.13,
        name,
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=6.5,
        fontweight=weight,
        color=color,
        clip_on=False,
    )
    ax.text(
        -0.018,
        y - 0.13,
        f"n = {n}    {n_pos} pos / {n_neg} neg",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="top",
        fontsize=5.5,
        color=PALETTE["neutral_dark"],
        clip_on=False,
    )


def _de_offsets(split: str) -> list[tuple[str, float]]:
    if split == "q4q1":
        cohorts = ["GSE123902", "GSE131907", "GSE205335", "Pooled"]
    elif split == "continuous":
        cohorts = ["GSE123902", "GSE131907", "GSE205335", "GSE189357", "Pooled"]
    else:
        raise ValueError(split)
    k = len(cohorts)
    # Keep the cluster inside ±0.38 so adjacent families (pitch 1) do not touch.
    offsets = np.linspace(-0.34, 0.34, k)
    return list(zip(cohorts, offsets))


def plot_de(ax, de: dict[tuple[str, str, str], dict[str, str]], split: str, title: str) -> None:
    y_of = {fam: len(FAMILY_ORDER) - 1 - i for i, fam in enumerate(FAMILY_ORDER)}
    for fam, y in y_of.items():
        ax.axhspan(y - 0.46, y + 0.46, color=PALETTE["band"], zorder=0, linewidth=0)
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.6, zorder=1, linestyle=(0, (2.2, 1.6)))

    ns: dict[str, int] = {}
    for fam in FAMILY_ORDER:
        y0 = y_of[fam]
        for cohort, off in _de_offsets(split):
            if cohort == "Pooled":
                row = de[(fam, split, POOLED_KEY)]
                color = POOLED_COLOR
                marker = "D"
                ms = 4.3
            else:
                row = de[(fam, split, cohort)]
                color = COHORT_COLORS[cohort]
                marker = "o"
                ms = 3.8
            logfc = fnum(row, "logFC")
            se = fnum(row, "se")
            fdr = fnum(row, "fdr")
            n = int(row["n"])
            if cohort not in ns:
                ns[cohort] = n
            elif ns[cohort] != n:
                raise SystemExit(f"n differs across families for {split} {cohort}")
            # Slightly shorter caps than panel a; five rows share one family band.
            ax.plot(
                [logfc - se, logfc + se],
                [y0 + off, y0 + off],
                color=color,
                lw=0.85,
                zorder=2,
                solid_capstyle="butt",
            )
            cap = 0.045
            yy = y0 + off
            ax.plot([logfc - se, logfc - se], [yy - cap, yy + cap], color=color, lw=0.65, zorder=2)
            ax.plot([logfc + se, logfc + se], [yy - cap, yy + cap], color=color, lw=0.65, zorder=2)
            ax.plot(
                logfc,
                yy,
                marker,
                markersize=ms,
                markerfacecolor=color if fdr < 0.05 else "white",
                markeredgecolor=color,
                markeredgewidth=0.95,
                zorder=4,
            )

    ax.set_yticks([y_of[f] for f in FAMILY_ORDER])
    ax.set_yticklabels(FAMILY_ORDER)
    ax.set_xlim(-2.55, 1.05)
    ax.set_ylim(-0.62, 3.72)
    ax.set_xticks([-2, -1, 0, 1])
    ax.set_xlabel("logFC")
    ax.set_title(title, loc="left", pad=6, fontsize=7)
    style_ax(ax)
    ax.tick_params(axis="y", length=0, pad=2)

    # Sample size is constant across families within a split; keep it in the legend.
    handles = []
    labels = []
    cohorts = [c for c, _ in _de_offsets(split)]
    for cohort in cohorts:
        color = POOLED_COLOR if cohort == "Pooled" else COHORT_COLORS[cohort]
        marker = "D" if cohort == "Pooled" else "o"
        handles.append(
            mpl.lines.Line2D(
                [0],
                [0],
                marker=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=3.8,
                lw=0,
            )
        )
        labels.append(f"{cohort if cohort != 'Pooled' else 'Pooled'} (n = {ns[cohort]})")
    ax.legend(
        handles,
        labels,
        loc="upper right",
        frameon=False,
        fontsize=5.5,
        handlelength=0.8,
        handletextpad=0.35,
        borderaxespad=0.15,
        labelspacing=0.28,
    )


def build() -> plt.Figure:
    apply_style()
    pairs, agg = load_barrier()
    de = load_de()

    fig = plt.figure(figsize=(mm(183), mm(158)))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.18, 1.0],
        width_ratios=[1.0, 1.0],
        left=0.210,
        right=0.978,
        top=0.925,
        bottom=0.112,
        hspace=0.48,
        wspace=0.32,
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1], sharey=ax_b)

    plot_barrier(ax_a, pairs, agg)
    plot_de(ax_b, de, "q4q1", "Q4 versus Q1")
    plot_de(ax_c, de, "continuous", "Continuous")
    ax_c.tick_params(labelleft=False)

    panel_letter(ax_a, "a", -0.175, 1.01)
    panel_letter(ax_b, "b", -0.28, 1.02)
    panel_letter(ax_c, "c", -0.08, 1.02)

    fig.text(
        0.205,
        0.028,
        "Filled marker, tabulated FDR < 0.05. Open marker, FDR ≥ 0.05. Bars, ± 1 SE. "
        "GSE189357 has no separate Q4-versus-Q1 rows in this table.",
        ha="left",
        va="bottom",
        fontsize=6,
        color=PALETTE["neutral_dark"],
    )
    return fig


def save(fig: plt.Figure) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / "fig_barrier_ligands"
    for ext in ("svg", "pdf"):
        fig.savefig(stem.with_suffix(f".{ext}"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    # Small preview for visual QA; same pixels as the publication png, just resampled on write.
    preview = OUT / "preview"
    preview.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview / "fig_barrier_ligands.png", dpi=160, bbox_inches="tight", pad_inches=0.04, facecolor="white")


def write_caption() -> None:
    """Caption states only values present in the two source tables."""
    pairs, agg = load_barrier()
    de = load_de()
    lines = [
        "Figure. Barrier and inhibitory ligand–receptor shifts and family differential expression.",
        "",
        "a, Mean Δ for seven barrier/inhibitory ligand–receptor pairs and the family aggregate "
        "(family_barrier_inhibitory.tsv; split q4q1). Circles (pairs) and the diamond (family aggregate) "
        "are mean_delta; horizontal bars are ± 1 SE (column se); vertical ticks are median_delta. "
        "Positive mean_delta matches the tabulated direction high > low (observed_direction; thesis_expect high>low). "
        "The axis is labeled Q4 − Q1 because the split is q4q1. Meta P is the tabulated meta_p "
        "(not recomputed). Row color is the ligand axis (F11R, NECTIN2–TIGIT, CDH1, LGALS9). "
        "n is n_units; pos/neg are n_pos/n_neg. Cohort unit counts sum to n_units but are not plotted, "
        "because the table has no per-cohort mean_delta. For LGALS9–HAVCR2, n_pos + n_neg = 26 and n_units = 27. "
        "Tabulated I² is 0 or < 10⁻¹³ on every row.",
        "",
        "b, c, Family logFC (family_de.tsv) for the Q4-versus-Q1 split and the continuous split. "
        "Point is logFC; bars are ± 1 SE (column se). Fill denotes the tabulated FDR: filled, FDR < 0.05; "
        "open, FDR ≥ 0.05. Cohort colors: GSE123902, GSE131907, GSE205335, GSE189357. "
        "Diamonds are the pooled row whose cohort field is GSE123902+GSE131907+GSE205335+GSE189357. "
        "GSE189357 has continuous rows and no separate Q4-versus-Q1 rows. "
        "The three Q4-versus-Q1 cohort rows have n = 7, 11 and 11 (sum 29); the pooled Q4-versus-Q1 row has n = 34. "
        "Continuous cohort n sums to the pooled n (13 + 21 + 21 + 9 = 64). "
        "n_genes is constant within a family: chemokine 25, MHC-I/APM 21, IFN 221, TJ 194. "
        "No confidence intervals, effect sizes, or cohort-level ligand Δ values were imputed.",
        "",
        "Source files: data/family_barrier_inhibitory.tsv, data/family_de.tsv.",
        "",
        "Plotted barrier rows (mean_delta, se, median_delta, meta_p, n_units, n_pos, n_neg):",
    ]
    for row in pairs + [agg]:
        name = "Family aggregate" if row["row"] == "FAMILY_AGGREGATE" else pair_label(row)
        lines.append(
            f"  {name}: {fnum(row, 'mean_delta'):.6g}, {fnum(row, 'se'):.6g}, "
            f"{fnum(row, 'median_delta'):.6g}, {fnum(row, 'meta_p'):.6g}, "
            f"{int(row['n_units'])}, {int(row['n_pos'])}, {int(row['n_neg'])}"
        )
    lines.append("")
    lines.append("Plotted family logFC (split, cohort, logFC, se, fdr, n):")
    for split in ("q4q1", "continuous"):
        for fam in FAMILY_ORDER:
            cohorts = [c for c, _ in _de_offsets(split)]
            for cohort in cohorts:
                key_cohort = POOLED_KEY if cohort == "Pooled" else cohort
                row = de[(fam, split, key_cohort)]
                lines.append(
                    f"  {fam} | {split} | {key_cohort}: {fnum(row, 'logFC'):.6g}, "
                    f"{fnum(row, 'se'):.6g}, {fnum(row, 'fdr'):.6g}, {int(row['n'])}"
                )
    (OUT / "figure_caption.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    fig = build()
    save(fig)
    plt.close(fig)
    write_caption()


if __name__ == "__main__":
    main()
