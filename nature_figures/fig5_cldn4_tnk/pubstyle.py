"""Shared Nature-figure style, palette, and export helpers (Python only)."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_3": "#8BCF8B",
    "red_strong": "#B64342",
    "red_1": "#F6CFCB",
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

SOURCE_DIR = Path("/workspace/nature-figures/source")
OUT_DIR = Path("/workspace/nature-figures/out")
PREVIEW_DIR = OUT_DIR / "preview"

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def apply_publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
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
        }
    )


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def unicode_sup(n: int) -> str:
    return str(int(n)).translate(_SUPERSCRIPT)


def fmt_p(p: float, digits: int = 1) -> str:
    """P-value as unicode scientific notation (no mathtext, keeps glyphs >= 5 pt)."""
    p = float(p)
    if p <= 0:
        return "P < 10⁻³⁰⁰"
    if p >= 0.001:
        if p >= 0.1:
            return f"P = {p:.2f}"
        return f"P = {p:.3f}"
    exp = int(np.floor(np.log10(p)))
    mant = p / (10 ** exp)
    return f"P = {mant:.{digits}f} × 10{unicode_sup(exp)}"


def fmt_rho(rho: float, digits: int = 2) -> str:
    sign = "−" if rho < 0 else ""
    return f"{sign}{abs(rho):.{digits}f}"


def add_panel_letter(fig, ax, letter: str, xshift: float = -0.055, yshift: float = 0.018) -> None:
    """Bold 8 pt lowercase panel letter placed outside the axes."""
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


def save_pub(fig, stem: str, dpi: int = 600) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for ext in ("svg", "pdf"):
        path = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.03, facecolor="white")
        written.append(path)
    for ext in ("png", "tiff"):
        path = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.03, facecolor="white")
        written.append(path)
    preview = PREVIEW_DIR / f"{stem}.png"
    fig.savefig(preview, dpi=dpi, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    written.append(preview)
    return written


def forest_points(ax, y, estimates, n_list, colors, ci_lo=None, ci_hi=None, markersize=4.2):
    """Point estimates; draw CI only when both bounds are finite (never invent)."""
    for i, (yy, est, n, col) in enumerate(zip(y, estimates, n_list, colors)):
        if ci_lo is not None and ci_hi is not None:
            lo = ci_lo[i]
            hi = ci_hi[i]
            if lo is not None and hi is not None and np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [yy, yy], color=col, lw=1.1, solid_capstyle="round", zorder=2)
        ax.plot(
            est,
            yy,
            "o",
            color=col,
            markersize=markersize,
            markeredgecolor="white",
            markeredgewidth=0.35,
            zorder=3,
        )


def forest_diamond(ax, y, x, lo, hi, color, height=0.32):
    xs = [lo, x, hi, x]
    ys = [y, y + height / 2.0, y, y - height / 2.0]
    ax.fill(xs, ys, facecolor=color, edgecolor=color, lw=0.4, zorder=3, closed=True)
