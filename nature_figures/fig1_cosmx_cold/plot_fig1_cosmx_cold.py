#!/usr/bin/env python3
"""Fig. 1 — CosMx NSCLC: TROP2/CLDN4-high neighborhoods and CD8+NK neighbors.

Drawn only from section-level tables already in the repository.
Does not recompute neighborhoods and does not draw the cited 0.36/0.52
handoff, which has no section rows here.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source"
OUT = ROOT

PALETTE = {
    "blue_main": "#0F4D92",
    "red_strong": "#B64342",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "teal": "#42949E",
    "violet": "#9A4D8E",
}

PATIENT_COLOR = {
    "Lung5": PALETTE["blue_main"],
    "Lung6": PALETTE["red_strong"],
    "Lung9": PALETTE["teal"],
    "Lung12": PALETTE["violet"],
    "Lung13": PALETTE["neutral_dark"],
}


def canonical_patient(value: str) -> str:
    text = value.strip()
    if text.startswith("P") and text[1:].isdigit():
        return "Lung" + text[1:]
    return text

PATIENT_ORDER = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
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
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
            "pdf.compression": 9,
        }
    )


def read_csv(name: str) -> list[dict]:
    with (SRC / name).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if "patient" in row:
            row["patient"] = canonical_patient(row["patient"])
    return rows


def rep_of(sample: str) -> str:
    text = sample.replace("Rep", "R").replace("_", " ")
    if "R1" in text:
        return "R1"
    if "R2" in text:
        return "R2"
    if "R3" in text:
        return "R3"
    return ""


def marker_of(rep: str) -> str:
    return {"R1": "o", "R2": "s", "R3": "D"}.get(rep, "o")


def sign_p_one_sided(n_down: int, n: int) -> float:
    """Exact one-sided sign probability under a fair coin, no ties."""
    if n_down != n or n <= 0:
        # binomial tail P(X >= n_down)
        from math import comb

        return sum(comb(n, k) for k in range(n_down, n + 1)) / (2**n)
    return 1.0 / (2**n)


def equal_weight_ratio(pairs: list[tuple[float, float]]) -> float:
    high = float(np.mean([h for h, _ in pairs]))
    low = float(np.mean([l for _, l in pairs]))
    return high / low


def patient_means(rows: list[dict]) -> dict[str, tuple[float, float]]:
    buckets: dict[str, list[tuple[float, float]]] = {p: [] for p in PATIENT_ORDER}
    for row in rows:
        buckets[row["patient"]].append((float(row["high_mean"]), float(row["low_mean"])))
    out = {}
    for patient, pairs in buckets.items():
        if not pairs:
            continue
        out[patient] = (
            float(np.mean([h for h, _ in pairs])),
            float(np.mean([l for _, l in pairs])),
        )
    return out


def concordance(rows: list[dict]) -> dict:
    n_sec = sum(float(r["high_mean"]) < float(r["low_mean"]) for r in rows)
    patients = patient_means(rows)
    n_pat = sum(h < l for h, l in patients.values())
    return {
        "n_sections": len(rows),
        "sections_high_lt_low": n_sec,
        "n_patients": len(patients),
        "patients_high_lt_low": n_pat,
        "ratio_of_means": equal_weight_ratio(
            [(float(r["high_mean"]), float(r["low_mean"])) for r in rows]
        ),
        "patient_ratio_of_means": equal_weight_ratio(list(patients.values())),
        "sign_p_patients": sign_p_one_sided(n_pat, len(patients)),
        "sign_p_sections": sign_p_one_sided(n_sec, len(rows)),
        "median_section_ratio": float(np.median([float(r["ratio"]) for r in rows])),
    }


def style_ax(ax) -> None:
    ax.tick_params(axis="both", which="major", labelsize=6, length=2.2, width=0.6, pad=1.4)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.6)
        ax.spines[spine].set_color(PALETTE["neutral_black"])
    ax.xaxis.label.set_color(PALETTE["neutral_black"])
    ax.yaxis.label.set_color(PALETTE["neutral_black"])
    ax.tick_params(colors=PALETTE["neutral_black"])
    ax.title.set_color(PALETTE["neutral_black"])


def panel_letter(fig, ax, letter: str) -> None:
    pos = ax.get_position()
    fig.text(
        pos.x0 - 0.018,
        pos.y1 + 0.012,
        letter,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=PALETTE["neutral_black"],
    )


def select(rows: list[dict], **want) -> list[dict]:
    out = []
    for row in rows:
        if all(str(row[k]) == str(v) for k, v in want.items()):
            out.append(row)
    return out


def draw_paired(ax, rows: list[dict], radii: list[str], radius_labels: list[str]) -> dict:
    grouped = {rad: select(rows, radius_um=rad) for rad in radii}
    positions = {}
    x = 0.0
    xticks = []
    xticklabels = []
    for rad in radii:
        positions[rad] = (x, x + 0.78)
        xticks.extend(positions[rad])
        xticklabels.extend(["Low", "High"])
        x += 2.05

    for rad in radii:
        by_sample = {r["sample"]: r for r in grouped[rad]}
        for sample, row in by_sample.items():
            color = PATIENT_COLOR[row["patient"]]
            marker = marker_of(rep_of(sample))
            xl, xh = positions[rad]
            low = float(row["low_mean"])
            high = float(row["high_mean"])
            ax.plot([xl, xh], [low, high], color=color, lw=0.7, zorder=2, alpha=0.95)
            ax.plot(
                [xl, xh],
                [low, high],
                linestyle="none",
                marker=marker,
                color=color,
                markersize=3.8,
                markeredgecolor="white",
                markeredgewidth=0.35,
                zorder=3,
            )
        lows = [float(r["low_mean"]) for r in grouped[rad]]
        highs = [float(r["high_mean"]) for r in grouped[rad]]
        xl, xh = positions[rad]
        mean_low = float(np.mean(lows))
        mean_high = float(np.mean(highs))
        for xpos, mean in ((xl, mean_low), (xh, mean_high)):
            ax.plot([xpos - 0.16, xpos + 0.16], [mean, mean], color="black", lw=1.05, zorder=4, solid_capstyle="butt")

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yscale("log")
    ax.set_xlim(-0.45, positions[radii[-1]][1] + 0.45)
    style_ax(ax)

    # Radius labels sit under the Low/High pair.
    for rad, label in zip(radii, radius_labels):
        xl, xh = positions[rad]
        ax.text(
            (xl + xh) / 2,
            -0.22,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=6,
            color=PALETTE["neutral_black"],
        )
    stats = {rad: concordance(grouped[rad]) for rad in radii}
    return stats


def annotate_ratios(ax, rows: list[dict], radii: list[str], y_mult: float) -> None:
    """Place the ratio of equal-weight section means above each pair."""
    # Recompute x positions the same way as draw_paired.
    x = 0.0
    for rad in radii:
        local = [r for r in rows if r["radius_um"] == rad]
        local_max = max(max(float(r["high_mean"]), float(r["low_mean"])) for r in local)
        xl, xh = x, x + 0.78
        sub = select(rows, radius_um=rad)
        ratio = concordance(sub)["ratio_of_means"]
        ax.text(
            (xl + xh) / 2,
            local_max * y_mult,
            f"{ratio:.2f}×",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color=PALETTE["neutral_black"],
            clip_on=False,
        )
        x += 2.05


def draw_ratio_strip(ax, blocks: list[tuple[str, list[dict], list[str]]]) -> None:
    x = 0.0
    xticks = []
    xticklabels = []
    group_spans = []
    for title, rows, radii in blocks:
        start = x
        for rad in radii:
            sub = select(rows, radius_um=rad)
            ordered = sorted(sub, key=lambda r: (PATIENT_ORDER.index(r["patient"]), rep_of(r["sample"])))
            n = len(ordered)
            for i, row in enumerate(ordered):
                jitter = (i - (n - 1) / 2) * 0.055
                ax.plot(
                    x + jitter,
                    float(row["ratio"]),
                    marker=marker_of(rep_of(row["sample"])),
                    color=PATIENT_COLOR[row["patient"]],
                    markersize=3.5,
                    markeredgecolor="white",
                    markeredgewidth=0.3,
                    linestyle="none",
                    zorder=3,
                )
            ratio_means = concordance(sub)["ratio_of_means"]
            ax.plot(
                x,
                ratio_means,
                marker="D",
                color="black",
                markersize=3.2,
                linestyle="none",
                zorder=4,
            )
            n_down = concordance(sub)["sections_high_lt_low"]
            ax.text(
                x,
                1.32,
                f"{n_down}/8",
                ha="center",
                va="bottom",
                fontsize=5.5,
                color=PALETTE["neutral_dark"],
                clip_on=False,
            )
            xticks.append(x)
            xticklabels.append(str(rad))
            x += 0.72
        group_spans.append((title, start, x - 0.72))
        x += 0.55

    ax.axhline(1.0, color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), zorder=1)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlim(-0.45, x - 0.55 + 0.35)
    ax.set_ylim(0, 1.58)
    ax.set_ylabel("Section ratio, high / low")
    ax.set_xlabel("Radius (µm)")
    for title, left, right in group_spans:
        ax.text(
            (left + right) / 2,
            -0.42,
            title,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=6,
            color=PALETTE["neutral_black"],
            clip_on=False,
        )
    style_ax(ax)


def draw_muzzling(ax, rows: list[dict]) -> None:
    genes = ["GZMB", "PRF1", "NKG7", "IFNG"]
    radii = ["50", "100"]
    offsets = {"50": -0.16, "100": 0.16}
    for gi, gene in enumerate(genes):
        for rad in radii:
            sub = [r for r in rows if r["gene"] == gene and r["radius_um"] == rad]
            ordered = sorted(sub, key=lambda r: (PATIENT_ORDER.index(r["patient"]), rep_of(r["sample"])))
            n = len(ordered)
            for i, row in enumerate(ordered):
                jitter = (i - (n - 1) / 2) * 0.028
                face = PATIENT_COLOR[row["patient"]] if rad == "100" else "white"
                ax.plot(
                    gi + offsets[rad] + jitter,
                    float(row["ratio_cpm"]),
                    marker=marker_of(rep_of(row["sample"])),
                    markerfacecolor=face,
                    markeredgecolor=PATIENT_COLOR[row["patient"]],
                    markersize=3.4,
                    markeredgewidth=0.7,
                    linestyle="none",
                    zorder=3,
                )
            med = float(np.median([float(r["ratio_cpm"]) for r in sub]))
            ax.plot(
                [gi + offsets[rad] - 0.08, gi + offsets[rad] + 0.08],
                [med, med],
                color="black",
                lw=1.0,
                zorder=4,
                solid_capstyle="butt",
            )
    ax.axhline(1.0, color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), zorder=1)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes)
    ax.set_xlim(-0.55, len(genes) - 0.45)
    ax.set_ylabel("CPM ratio, high / low")
    ax.set_xlabel("Transcript in nearby CD8 or NK cells")
    style_ax(ax)
    # Radius legend proxies.
    ax.plot([], [], marker="o", markerfacecolor="white", markeredgecolor=PALETTE["neutral_dark"],
            linestyle="none", markersize=3.6, markeredgewidth=0.7, label="50 µm")
    ax.plot([], [], marker="o", markerfacecolor=PALETTE["neutral_dark"], markeredgecolor=PALETTE["neutral_dark"],
            linestyle="none", markersize=3.6, label="100 µm")
    ax.legend(loc="upper right", fontsize=5.5, handletextpad=0.3, borderaxespad=0.2)


def patient_legend(fig) -> None:
    handles = []
    for patient in PATIENT_ORDER:
        handles.append(
            mpl.lines.Line2D(
                [],
                [],
                marker="o",
                color=PATIENT_COLOR[patient],
                linestyle="none",
                markersize=3.8,
                markeredgecolor="white",
                markeredgewidth=0.3,
                label=patient.replace("Lung", "Lung "),
            )
        )
    handles.append(
        mpl.lines.Line2D(
            [],
            [],
            marker="s",
            color=PALETTE["neutral_dark"],
            linestyle="none",
            markersize=3.6,
            label="Replicate 2",
        )
    )
    handles.append(
        mpl.lines.Line2D(
            [],
            [],
            marker="D",
            color=PALETTE["neutral_dark"],
            linestyle="none",
            markersize=3.4,
            label="Replicate 3",
        )
    )
    handles.append(
        mpl.lines.Line2D([], [], color="black", lw=1.05, label="Equal-weight mean")
    )
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=8,
        frameon=False,
        fontsize=6,
        bbox_to_anchor=(0.5, 0.005),
        handletextpad=0.3,
        columnspacing=1.0,
    )


def main() -> None:
    apply_style()
    marginal = read_csv("marginal_cd8nk_by_section.csv")
    q4 = read_csv("cldn4_q4q1_cd8nk_by_section.csv")
    muzzling = read_csv("cldn4_nearby_effector_cpm.csv")

    tac = select(marginal, gene="TACSTD2")
    cld = select(marginal, gene="CLDN4")

    fig = plt.figure(figsize=(180 / 25.4, 128 / 25.4))
    gs = fig.add_gridspec(
        2,
        2,
        left=0.07,
        right=0.985,
        top=0.90,
        bottom=0.14,
        wspace=0.34,
        hspace=0.62,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    tac_radii = ["10", "20"]
    stats_a = draw_paired(ax_a, tac, tac_radii, ["10 µm", "20 µm"])
    annotate_ratios(ax_a, tac, tac_radii, y_mult=1.45)
    ax_a.set_ylim(4e-4, 0.45)
    ax_a.set_ylabel("CD8+NK neighbors per tumor cell")
    ax_a.set_title(
        "TACSTD2, section median\n8/8 sections, 5/5 patients, sign P = 0.031",
        loc="left",
        pad=4,
        fontsize=7,
    )

    q_radii = ["20", "40"]
    stats_b = draw_paired(ax_b, q4, q_radii, ["20 µm", "40 µm"])
    annotate_ratios(ax_b, q4, q_radii, y_mult=1.55)
    ax_b.set_ylim(0.03, 12)
    ax_b.set_ylabel("CD8+NK neighbors per tumor cell")
    ax_b.set_title(
        "CLDN4, top vs bottom quartile\n8/8 sections, 5/5 patients, sign P = 0.031",
        loc="left",
        pad=4,
        fontsize=7,
    )

    draw_ratio_strip(
        ax_c,
        [
            ("TACSTD2 median", tac, ["10", "20", "50", "100"]),
            ("CLDN4 median", cld, ["10", "20", "50", "100"]),
        ],
    )
    ax_c.set_title("Same median split across radius", loc="left", pad=2, fontsize=7)
    ax_c.text(
        0.98,
        0.04,
        "Diamond: ratio of means",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.5,
        color=PALETTE["neutral_dark"],
    )

    draw_muzzling(ax_d, muzzling)
    ax_d.set_title("Effector transcripts near CLDN4-high tumor", loc="left", pad=2, fontsize=7)
    ax_d.set_ylim(0.2, 1.85)

    panel_letter(fig, ax_a, "a")
    panel_letter(fig, ax_b, "b")
    panel_letter(fig, ax_c, "c")
    panel_letter(fig, ax_d, "d")
    patient_legend(fig)

    summary = {
        "tacstd2_median_cd8nk": {rad: stats_a[rad] for rad in tac_radii},
        "cldn4_median_cd8nk": {rad: concordance(select(cld, radius_um=rad)) for rad in ["10", "20", "50", "100"]},
        "cldn4_q4q1_cd8nk": {rad: stats_b[rad] for rad in q_radii},
        "cldn4_q4q1_cd8nk_60_80": {
            rad: concordance(select(q4, radius_um=rad)) for rad in ["60", "80"]
        },
        "cited_handoff_not_plotted": {
            "ratio_50um": 0.36,
            "ratio_100um": 0.52,
            "reason": "No section-level rows for this summary are in the repository tables.",
        },
    }
    # JSON-friendly floats
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, float):
            return round(obj, 6)
        return obj

    expected = {
        ("tacstd2_median_cd8nk", "10"): 0.455420,
        ("tacstd2_median_cd8nk", "20"): 0.670165,
        ("cldn4_q4q1_cd8nk", "20"): 0.459808,
        ("cldn4_q4q1_cd8nk", "40"): 0.574290,
        ("cldn4_median_cd8nk", "50"): 0.810267,
        ("cldn4_median_cd8nk", "100"): 0.911636,
    }
    for (block, rad), target in expected.items():
        got = summary[block][rad]["ratio_of_means"]
        if abs(got - target) > 1e-4:
            raise SystemExit(f"ratio drift {block} {rad}: {got} != {target}")
        if summary[block][rad]["patients_high_lt_low"] == summary[block][rad]["n_patients"]:
            if abs(summary[block][rad]["sign_p_patients"] - 0.03125) > 1e-9:
                raise SystemExit(f"sign P drift {block} {rad}")

    (OUT / "source_stats.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")

    for ext in ("pdf", "svg", "png"):
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.04, "facecolor": "white"}
        if ext == "png":
            kwargs["dpi"] = 600
        fig.savefig(OUT / f"fig1_cosmx_cold.{ext}", **kwargs)
    plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
