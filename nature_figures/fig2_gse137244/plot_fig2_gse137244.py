#!/usr/bin/env python3
"""Fig. 2 — GSE137244 bulk KL vs KP: Tacstd2, Cldn4, epithelial tight junction.

Bulk cell-line RNA-seq only. Not single-cell. Not the private 8KL atlas.
Every plotted point is read from source_data/. The script stops if those
tables do not recover Tacstd2 Δ = +3.24, Cldn4 Δ = +5.57, and two-sided
Mann–Whitney P = 0.00794 at n = 5 vs 5.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
SRC = HERE / "source_data"

KP = "#3775BA"
KL = "#B64342"
INK = "#272727"
MEAN = "#272727"

# Genes drawn in the heatmap. Membership is fixed; order is by Δ from the table.
# Tacstd2 is the surface gene in panels a–b, not a core TJ gene.
# Cldn6 is a claudin callout and is not one of the 18 core genes.
KL_BLOCK = ("Tacstd2", "Cldn6", "Cldn4", "Cldn7", "Cldn3", "Cldn1", "Ocln", "Tjp1")
KP_BLOCK = ("Cldn2", "Cldn5", "Cldn18")


def register_arial() -> None:
    fonts = [
        Path("/home/ubuntu/.local/share/fonts/Arial.TTF"),
        Path("/home/ubuntu/.local/share/fonts/Arialbd.TTF"),
        Path("/home/ubuntu/.local/share/fonts/Ariali.TTF"),
        Path("/home/ubuntu/.local/share/fonts/Arialbi.TTF"),
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf"),
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf"),
        Path("/Library/Fonts/Arial.ttf"),
    ]
    found = False
    for path in fonts:
        if path.is_file():
            mpl.font_manager.fontManager.addfont(str(path))
            found = True
    if not found:
        raise SystemExit("Arial is not installed; refusing to substitute another family.")
    chosen = mpl.font_manager.findfont("Arial", fallback_to_default=False)
    if "arial" not in Path(chosen).name.lower():
        raise SystemExit(f"Matplotlib did not select Arial (got {chosen}).")


def apply_style() -> None:
    register_arial()
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "pdf.compression": 9,
        }
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fmt_delta(delta: float) -> str:
    text = f"{delta:+.2f}"
    return text.replace("-", "−")


def fmt_p(p: float) -> str:
    return f"P = {p:.5f}"


def exact_mw(kl: np.ndarray, kp: np.ndarray) -> tuple[float, float]:
    result = mannwhitneyu(kl, kp, alternative="two-sided", method="exact")
    return float(result.statistic), float(result.pvalue)


def load() -> dict:
    samples = read_tsv(SRC / "samples.tsv")
    if any(row["arm"] not in {"KP", "KL"} for row in samples):
        raise SystemExit("samples.tsv contains a library that is not KP or KL.")
    if len(samples) != 10:
        raise SystemExit(f"Expected 10 tumor libraries, found {len(samples)}.")

    watch = read_tsv(SRC / "watch_genes_log2.tsv")
    by_short = {row[""]: row for row in watch}
    if set(by_short) != {row["short"] for row in samples}:
        raise SystemExit("watch_genes_log2.tsv libraries do not match samples.tsv.")

    core_genes = [
        line.strip()
        for line in (SRC / "epithelial_tj_core_genes.txt").read_text().splitlines()
        if line.strip()
    ]
    if len(core_genes) != 18 or len(set(core_genes)) != 18:
        raise SystemExit("Epithelial TJ core must be 18 unique symbols.")

    scores = read_tsv(SRC / "mean_signature_scores.tsv")
    core_row = next(row for row in scores if row["set"] == "EPITHELIAL_TJ_CORE")
    if int(core_row["n_genes_in_matrix"]) != 18:
        raise SystemExit("EPITHELIAL_TJ_CORE is not the 18-gene score.")

    libraries = []
    for row in samples:
        short = row["short"]
        libraries.append(
            {
                "library": row["library"],
                "short": short,
                "label": short.split(" ", 1)[1],
                "gsm": row["gsm"],
                "arm": row["arm"],
                "genotype": row["genotype"],
                "watch": by_short[short],
                "tj_core": float(core_row[short]),
            }
        )

    callouts = {row["gene"]: row for row in read_tsv(SRC / "tj_gene_callouts.tsv")}
    locked = {row["gene"]: row for row in read_tsv(SRC / "locked_genes.tsv")}
    return {
        "libraries": libraries,
        "core_genes": core_genes,
        "core_row": core_row,
        "callouts": callouts,
        "locked": locked,
    }


def values_for(data: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    kp, kl = [], []
    for lib in data["libraries"]:
        if key == "TJ_CORE":
            value = lib["tj_core"]
        else:
            value = float(lib["watch"][key])
        (kl if lib["arm"] == "KL" else kp).append(value)
    return np.asarray(kl, dtype=float), np.asarray(kp, dtype=float)


def summarise(data: dict, key: str, kind: str) -> dict:
    kl, kp = values_for(data, key)
    if len(kl) != 5 or len(kp) != 5:
        raise SystemExit(f"{key} does not have 5 vs 5 values.")
    if np.any(kl < 0) or np.any(kp < 0):
        raise SystemExit(f"{key} has a negative log2(FPKM+1) value.")
    delta = float(kl.mean() - kp.mean())
    statistic, pvalue = exact_mw(kl, kp)
    separated = bool(kp.max() < kl.min() or kl.max() < kp.min())
    record = {
        "feature": "Epithelial TJ core" if key == "TJ_CORE" else key,
        "key": key,
        "kind": kind,
        "n_kp": 5,
        "n_kl": 5,
        "mean_kp": kp.mean(),
        "mean_kl": kl.mean(),
        "delta": delta,
        "U": statistic,
        "p_mw_two_sided": pvalue,
        "complete_separation": separated,
        "min_kp": float(kp.min()),
        "max_kp": float(kp.max()),
        "min_kl": float(kl.min()),
        "max_kl": float(kl.max()),
    }
    if key in data["callouts"]:
        published = float(data["callouts"][key]["delta_log2"])
        if abs(published - delta) > 1e-6:
            raise SystemExit(f"{key} Δ {delta} disagrees with tj_gene_callouts {published}.")
        published_p = float(data["callouts"][key]["mw_p"])
        if abs(published_p - pvalue) > 1e-9:
            raise SystemExit(f"{key} P disagrees with tj_gene_callouts.")
    if key in data["locked"]:
        published = float(data["locked"][key]["delta_log2"])
        if abs(published - delta) > 1e-6:
            raise SystemExit(f"{key} Δ disagrees with locked_genes.")
    if key == "TJ_CORE":
        published = float(data["core_row"]["library_delta"])
        if abs(published - delta) > 1e-6:
            raise SystemExit("TJ core Δ disagrees with mean_signature_scores.")
        if abs(float(data["core_row"]["library_p"]) - pvalue) > 1e-9:
            raise SystemExit("TJ core P disagrees with mean_signature_scores.")
    return record


def check_locked(summaries: dict[str, dict]) -> None:
    expected = {
        "Tacstd2": "+3.24",
        "Cldn4": "+5.57",
        "TJ_CORE": "+1.63",
    }
    for key, text in expected.items():
        got = fmt_delta(summaries[key]["delta"])
        if got != text:
            raise SystemExit(f"{key} displays as {got}, expected {text}.")
        if fmt_p(summaries[key]["p_mw_two_sided"]) != "P = 0.00794":
            raise SystemExit(f"{key} P is not the locked 0.00794.")
        if not summaries[key]["complete_separation"]:
            raise SystemExit(f"{key} does not completely separate; refusing a floor P.")
        if summaries[key]["U"] != 25:
            raise SystemExit(f"{key} U is {summaries[key]['U']}, expected 25.")


def style_strip(ax) -> None:
    ax.tick_params(axis="both", labelsize=6, length=2.4, width=0.7, pad=1.5, colors=INK)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    ax.yaxis.label.set_color(INK)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.7)
        ax.spines[spine].set_color(INK)
    ax.set_xlim(-0.62, 1.62)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["KP\nn = 5", "KL\nn = 5"])
    ax.get_xticklabels()[0].set_color(KP)
    ax.get_xticklabels()[1].set_color(KL)
    ax.set_ylabel("log2(FPKM + 1)")


def draw_strip(ax, data: dict, key: str, summary: dict, title: str, italic: bool) -> None:
    offsets = np.array([-0.16, -0.08, 0.0, 0.08, 0.16])
    kp_libs = [lib for lib in data["libraries"] if lib["arm"] == "KP"]
    kl_libs = [lib for lib in data["libraries"] if lib["arm"] == "KL"]
    for center, libs, color in ((0, kp_libs, KP), (1, kl_libs, KL)):
        ys = []
        for offset, lib in zip(offsets, libs):
            value = lib["tj_core"] if key == "TJ_CORE" else float(lib["watch"][key])
            ys.append(value)
            ax.plot(
                center + offset,
                value,
                "o",
                markersize=5.4,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=0.45,
                clip_on=False,
                zorder=3,
            )
        ax.plot(
            [center - 0.30, center + 0.30],
            [float(np.mean(ys)), float(np.mean(ys))],
            color=MEAN,
            lw=1.15,
            solid_capstyle="butt",
            zorder=2,
            clip_on=False,
        )
    ymax = max(summary["max_kp"], summary["max_kl"])
    floor = -0.55 if min(summary["min_kp"], summary["min_kl"]) < 0.05 else 0.0
    ax.set_ylim(floor, ymax * 1.14)
    span = ax.get_ylim()[1] - floor
    step = 2.0
    ticks = np.arange(0, ax.get_ylim()[1] + 1e-6, step)
    ticks = ticks[(ticks >= 0) & (ticks <= ax.get_ylim()[1] - 0.08 * span)]
    ax.set_yticks(ticks)
    style_strip(ax)
    ax.set_title(
        title,
        fontsize=7,
        pad=20,
        color=INK,
        fontstyle="italic" if italic else "normal",
    )
    ax.text(
        0.5,
        1.02,
        f"Δ = {fmt_delta(summary['delta'])}\n{fmt_p(summary['p_mw_two_sided'])}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6,
        color=INK,
        linespacing=1.15,
        clip_on=False,
        zorder=4,
    )


def panel_letter(ax, letter: str, x: float = -0.30, y: float = 1.36) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=INK,
        clip_on=False,
    )


def draw_heatmap(ax, data: dict, summaries: dict[str, dict]) -> list[str]:
    kl_genes = sorted(KL_BLOCK, key=lambda gene: summaries[gene]["delta"], reverse=True)
    kp_genes = sorted(KP_BLOCK, key=lambda gene: summaries[gene]["delta"], reverse=True)
    # Top of the axes is the KL-higher block. A gap separates the KP-higher block.
    gap = 0.62
    y_of: dict[str, float] = {}
    y = 0.0
    for gene in reversed(kp_genes):
        y_of[gene] = y
        y += 1.0
    y += gap
    for gene in reversed(kl_genes):
        y_of[gene] = y
        y += 1.0
    top = y

    cmap = LinearSegmentedColormap.from_list(
        "gse_blue",
        ["#F4F6F8", "#D5E3F2", "#7FA6D0", "#2E6399", "#0F4D92"],
    )
    norm = mpl.colors.Normalize(vmin=0, vmax=10)
    libs = data["libraries"]
    for gene, gy in y_of.items():
        for j, lib in enumerate(libs):
            value = float(lib["watch"][gene])
            ax.add_patch(
                Rectangle(
                    (j, gy),
                    1,
                    0.92,
                    facecolor=cmap(norm(value)),
                    edgecolor="white",
                    linewidth=0.7,
                    zorder=2,
                )
            )
        ax.text(
            len(libs) + 0.18,
            gy + 0.46,
            fmt_delta(summaries[gene]["delta"]),
            ha="left",
            va="center",
            fontsize=6,
            color=INK,
            clip_on=False,
        )

    # Genotype strip, tall enough for 6 pt type.
    strip_h = 1.15
    strip_y = top + 0.18
    for j, lib in enumerate(libs):
        ax.add_patch(
            Rectangle(
                (j, strip_y),
                1,
                strip_h,
                facecolor=KP if lib["arm"] == "KP" else KL,
                edgecolor="white",
                linewidth=0.7,
                zorder=2,
            )
        )
    ax.text(2.5, strip_y + strip_h / 2, "KP   n = 5", ha="center", va="center", fontsize=6, color="white", zorder=3)
    ax.text(7.5, strip_y + strip_h / 2, "KL   n = 5", ha="center", va="center", fontsize=6, color="white", zorder=3)
    ax.plot([5, 5], [0, top], color="white", lw=1.1, zorder=4)
    for genes in (kl_genes, kp_genes):
        y0 = min(y_of[gene] for gene in genes)
        y1 = max(y_of[gene] for gene in genes) + 0.92
        ax.add_patch(
            Rectangle((0, y0), len(libs), y1 - y0, fill=False, edgecolor=INK, linewidth=0.6, zorder=5)
        )

    kl_ys = [y_of[gene] + 0.46 for gene in kl_genes]
    kp_ys = [y_of[gene] + 0.46 for gene in kp_genes]
    ax.text(
        -1.05,
        float(np.mean(kl_ys)),
        "Higher in KL",
        rotation=90,
        ha="center",
        va="center",
        fontsize=6,
        color=INK,
        clip_on=False,
    )
    ax.text(
        -1.05,
        float(np.mean(kp_ys)),
        "Higher in KP",
        rotation=90,
        ha="center",
        va="center",
        fontsize=6,
        color=INK,
        clip_on=False,
    )

    ax.set_xlim(0, len(libs))
    ax.set_ylim(-0.02, strip_y + strip_h)
    ax.set_xticks([j + 0.5 for j in range(len(libs))])
    ax.set_xticklabels([lib["label"] for lib in libs], rotation=45, ha="right", rotation_mode="anchor")
    label_genes = list(kp_genes) + list(kl_genes)
    ax.set_yticks([y_of[gene] + 0.46 for gene in label_genes])
    ax.set_yticklabels(label_genes, fontstyle="italic")
    ax.tick_params(axis="both", length=0, pad=2, labelsize=6, colors=INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        len(libs) + 0.18,
        strip_y + strip_h / 2,
        "Δ",
        ha="left",
        va="center",
        fontsize=6,
        color=INK,
        clip_on=False,
    )
    ax.set_xlabel("Library", fontsize=7, labelpad=2, color=INK)
    return [cmap, norm]


def export_tables(data: dict, summaries: dict[str, dict]) -> None:
    plotted = []
    for lib in data["libraries"]:
        for key, panel in (
            ("Tacstd2", "a"),
            ("Cldn4", "b"),
            ("TJ_CORE", "c"),
        ):
            value = lib["tj_core"] if key == "TJ_CORE" else float(lib["watch"][key])
            plotted.append(
                {
                    "panel": panel,
                    "feature": "Epithelial_TJ_core" if key == "TJ_CORE" else key,
                    "library": lib["library"],
                    "label": lib["label"],
                    "gsm": lib["gsm"],
                    "arm": lib["arm"],
                    "genotype": lib["genotype"],
                    "log2_fpkm_plus1": f"{value:.10f}",
                }
            )
        for gene in KL_BLOCK + KP_BLOCK:
            plotted.append(
                {
                    "panel": "d",
                    "feature": gene,
                    "library": lib["library"],
                    "label": lib["label"],
                    "gsm": lib["gsm"],
                    "arm": lib["arm"],
                    "genotype": lib["genotype"],
                    "log2_fpkm_plus1": f"{float(lib['watch'][gene]):.10f}",
                }
            )
    write_tsv(
        SRC / "fig2_plotted_values.tsv",
        plotted,
        ["panel", "feature", "library", "label", "gsm", "arm", "genotype", "log2_fpkm_plus1"],
    )
    summary_rows = []
    for key in ("Tacstd2", "Cldn4", "TJ_CORE", *sorted(KL_BLOCK, key=lambda g: -summaries[g]["delta"]), *sorted(KP_BLOCK, key=lambda g: -summaries[g]["delta"])):
        row = summaries[key]
        summary_rows.append(
            {
                "feature": row["feature"],
                "kind": row["kind"],
                "n_kp": row["n_kp"],
                "n_kl": row["n_kl"],
                "mean_kp": f"{row['mean_kp']:.6f}",
                "mean_kl": f"{row['mean_kl']:.6f}",
                "delta_kl_minus_kp": f"{row['delta']:.6f}",
                "delta_display": fmt_delta(row["delta"]),
                "U": f"{row['U']:.0f}",
                "p_mw_two_sided": f"{row['p_mw_two_sided']:.10f}",
                "p_display": fmt_p(row["p_mw_two_sided"]),
                "complete_separation": row["complete_separation"],
            }
        )
    # de-duplicate TJ genes that are also headline features
    seen = set()
    unique = []
    for row in summary_rows:
        if row["feature"] in seen and row["kind"] != "headline":
            continue
        # Tacstd2/Cldn4 appear as headline and again in the heatmap block.
        tag = (row["feature"], row["kind"])
        if tag in seen:
            continue
        seen.add(tag)
        unique.append(row)
    write_tsv(
        SRC / "fig2_summary.tsv",
        unique,
        [
            "feature",
            "kind",
            "n_kp",
            "n_kl",
            "mean_kp",
            "mean_kl",
            "delta_kl_minus_kp",
            "delta_display",
            "U",
            "p_mw_two_sided",
            "p_display",
            "complete_separation",
        ],
    )

    watch_genes = set(data["libraries"][0]["watch"]) - {""}
    coverage = []
    for gene in data["core_genes"]:
        coverage.append(
            {
                "gene": gene,
                "in_epithelial_tj_core": True,
                "per_library_in_watch_genes": gene in watch_genes,
                "drawn_as_heatmap_cell": gene in KL_BLOCK or gene in KP_BLOCK,
            }
        )
    for gene in ("Tacstd2", "Cldn6", "Cldn2", "Cldn5", "Cldn18"):
        coverage.append(
            {
                "gene": gene,
                "in_epithelial_tj_core": False,
                "per_library_in_watch_genes": True,
                "drawn_as_heatmap_cell": True,
            }
        )
    write_tsv(
        SRC / "heatmap_gene_coverage.tsv",
        coverage,
        ["gene", "in_epithelial_tj_core", "per_library_in_watch_genes", "drawn_as_heatmap_cell"],
    )


def audit_layout(fig) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = []
    for text in fig.findobj(mpl.text.Text):
        label = text.get_text().replace("\n", " | ")
        if not label.strip():
            continue
        bbox = text.get_window_extent(renderer)
        texts.append((label, bbox, text.get_fontsize()))
    problems = []
    for i, (label, bbox, size) in enumerate(texts):
        if size < 6:
            problems.append(f"fontsize {size}: {label}")
        for other, obox, _ in texts[i + 1 :]:
            overlap = bbox.overlaps(obox)
            if not overlap:
                continue
            # Identical strings stacked on purpose (none expected).
            ix0 = max(bbox.x0, obox.x0)
            iy0 = max(bbox.y0, obox.y0)
            ix1 = min(bbox.x1, obox.x1)
            iy1 = min(bbox.y1, obox.y1)
            area = max(0, ix1 - ix0) * max(0, iy1 - iy0)
            if area > 8:
                problems.append(f"overlap {area:.0f}px '{label}' × '{other}'")
    if problems:
        details = "\n".join(f" - {item}" for item in problems)
        raise SystemExit(f"Figure layout problem:\n{details}")


def build() -> None:
    apply_style()
    data = load()
    summaries = {
        "Tacstd2": summarise(data, "Tacstd2", "headline"),
        "Cldn4": summarise(data, "Cldn4", "headline"),
        "TJ_CORE": summarise(data, "TJ_CORE", "headline"),
    }
    for gene in KL_BLOCK + KP_BLOCK:
        if gene not in summaries:
            summaries[gene] = summarise(data, gene, "heatmap")
    check_locked(summaries)
    for gene in KP_BLOCK:
        if summaries[gene]["U"] != 0:
            raise SystemExit(f"{gene} was expected to be higher in every KP library.")
        if fmt_p(summaries[gene]["p_mw_two_sided"]) != "P = 0.00794":
            raise SystemExit(f"{gene} P is not 0.00794.")

    fig = plt.figure(figsize=(183 / 25.4, 128 / 25.4))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.05, 1.18],
        hspace=0.52,
        left=0.10,
        right=0.90,
        top=0.82,
        bottom=0.09,
    )
    top = outer[0].subgridspec(1, 3, wspace=0.55)
    bottom = outer[1].subgridspec(1, 2, width_ratios=[1, 0.035], wspace=0.05)
    axes = [fig.add_subplot(top[0, i]) for i in range(3)]
    heat = fig.add_subplot(bottom[0, 0])
    cax = fig.add_subplot(bottom[0, 1])

    draw_strip(axes[0], data, "Tacstd2", summaries["Tacstd2"], "Tacstd2", True)
    draw_strip(axes[1], data, "Cldn4", summaries["Cldn4"], "Cldn4", True)
    draw_strip(axes[2], data, "TJ_CORE", summaries["TJ_CORE"], "Epithelial TJ core", False)
    axes[2].set_ylabel("Mean log2(FPKM + 1)")
    for ax, letter in zip(axes, "abc"):
        panel_letter(ax, letter)
    panel_letter(heat, "d", x=-0.055, y=1.08)

    cmap, norm = draw_heatmap(heat, data, summaries)
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, cax=cax)
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor(INK)
    cbar.set_ticks([0, 2, 4, 6, 8, 10])
    cbar.ax.tick_params(labelsize=6, width=0.6, length=2.2, colors=INK)
    cbar.set_label("log2(FPKM + 1)", fontsize=7, color=INK)

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=KP, markeredgecolor="white", markeredgewidth=0.4, markersize=5.2, label="KP  KrasG12D;Trp53"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=KL, markeredgecolor="white", markeredgewidth=0.4, markersize=5.2, label="KL  KrasG12D;Lkb1"),
        Line2D([0], [0], color=MEAN, lw=1.15, label="Mean"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=3,
        frameon=False,
        fontsize=6,
        handletextpad=0.4,
        columnspacing=1.2,
        borderaxespad=0,
    )

    export_tables(data, summaries)
    audit_layout(fig)
    written = []
    for ext in ("pdf", "svg", "png"):
        path = HERE / f"fig2_gse137244.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else None, bbox_inches="tight", pad_inches=0.04)
        written.append(path)
    plt.close(fig)

    panel_map = {
        "figure_id": "fig2_gse137244",
        "title": "GSE137244 bulk KL vs KP: Tacstd2, Cldn4, and epithelial tight junction",
        "assay": "bulk RNA-seq, mouse lung-cancer cell lines",
        "not": ["single-cell RNA-seq", "private 8KL", "methods schematic", "normal lung"],
        "accession": "GSE137244",
        "pmid": "34142094",
        "transform": "log2(FPKM+1)",
        "contrast": "mean(KL) - mean(KP)",
        "test": "two-sided exact Mann-Whitney U",
        "n_libraries": {"KP": 5, "KL": 5},
        "kp_note": "Five libraries of B6AL10 (KrasG12D;Trp53), not five independent lines.",
        "held_out": "GSM4073826 normal lung",
        "panels": [
            {
                "letter": "a",
                "title": "Tacstd2",
                "geom": "library dots and mean line",
                "n": "5 vs 5",
                "delta_display": fmt_delta(summaries["Tacstd2"]["delta"]),
                "p_display": fmt_p(summaries["Tacstd2"]["p_mw_two_sided"]),
                "source": "source_data/watch_genes_log2.tsv",
            },
            {
                "letter": "b",
                "title": "Cldn4",
                "geom": "library dots and mean line",
                "n": "5 vs 5",
                "delta_display": fmt_delta(summaries["Cldn4"]["delta"]),
                "p_display": fmt_p(summaries["Cldn4"]["p_mw_two_sided"]),
                "source": "source_data/watch_genes_log2.tsv",
            },
            {
                "letter": "c",
                "title": "Epithelial TJ core (18-gene mean)",
                "geom": "library dots and mean line",
                "n": "5 vs 5",
                "genes": data["core_genes"],
                "delta_display": fmt_delta(summaries["TJ_CORE"]["delta"]),
                "p_display": fmt_p(summaries["TJ_CORE"]["p_mw_two_sided"]),
                "source": "source_data/mean_signature_scores.tsv row EPITHELIAL_TJ_CORE",
                "not_recomputed": "Score is the published per-library mean. Genes without per-library rows are not imputed.",
            },
            {
                "letter": "d",
                "title": "Per-library log2(FPKM+1)",
                "geom": "heatmap of published per-library values",
                "higher_in_kl": sorted(KL_BLOCK, key=lambda g: summaries[g]["delta"], reverse=True),
                "higher_in_kp": sorted(KP_BLOCK, key=lambda g: summaries[g]["delta"], reverse=True),
                "omitted_core_genes_without_per_library_values": [
                    gene for gene in data["core_genes"] if gene not in KL_BLOCK and gene not in KP_BLOCK
                ],
                "source": "source_data/watch_genes_log2.tsv",
            },
        ],
        "outputs": [path.name for path in written],
    }
    (HERE / "panel_map.json").write_text(json.dumps(panel_map, indent=2) + "\n")
    for path in written:
        print(path)


if __name__ == "__main__":
    build()
