"""Nature-figure style used by the TACSTD2–junction–T/NK forests."""

from __future__ import annotations

import matplotlib as mpl
import numpy as np

PALETTE = {
    "blue_main": "#0F4D92",
    "green_3": "#8BCF8B",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "teal": "#42949E",
    "violet": "#9A4D8E",
}

COHORT_COLORS = {
    "GSE123902": PALETTE["blue_main"],
    "GSE131907": PALETTE["teal"],
    "GSE205335": PALETTE["violet"],
    "GSE189357": PALETTE["green_3"],
}

COHORT_ORDER = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def apply_publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arimo",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
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
        }
    )


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def unicode_sup(n: int) -> str:
    return str(int(n)).translate(_SUPERSCRIPT)


def fmt_p(p: float) -> str:
    p = float(p)
    if p < 0.001:
        exp = int(np.floor(np.log10(p)))
        mant = p / (10**exp)
        return f"P = {mant:.1f} × 10{unicode_sup(exp)}"
    if p >= 0.1:
        return f"P = {p:.2f}"
    return f"P = {p:.3f}"


def fmt_rho(rho: float) -> str:
    sign = "−" if rho < 0 else ""
    return f"{sign}{abs(rho):.2f}"


def style_ax(ax) -> None:
    ax.tick_params(axis="both", which="major", labelsize=6, length=2.2, width=0.6, pad=1.5)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.6)
        ax.spines[spine].set_color(PALETTE["neutral_black"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.label.set_color(PALETTE["neutral_black"])
    ax.yaxis.label.set_color(PALETTE["neutral_black"])
    ax.tick_params(colors=PALETTE["neutral_black"])
