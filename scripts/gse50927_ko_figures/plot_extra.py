#!/usr/bin/env python3
"""ADDITIVE high-contrast extra pack for GSE50927 Cldn4 KO.

Redraws locked numbers only. Does not compute new p / NES / FDR.
GSE50927 is naive mouse whole lung (not lung cancer).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "methods" / "gse50927_ko_figures"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)

# High-contrast dark pack
BG = "#070B12"
INK = "#F4F7FB"
MUTED = "#A8B3C7"
GRID = "#243044"
UP = "#FFC857"
DOWN = "#3D8BFF"
TARGET = "#FF3B6B"
IFN = "#FFB703"
MHC = "#2EE6D6"
TJ = "#8B93A7"
POS_NES = "#FFC857"
NEG_NES = "#4C7CFF"


def style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": BG,
            "savefig.facecolor": BG,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": MUTED,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, stem: str) -> Path:
    png = FIG / f"{stem}.png"
    pdf = FIG / f"{stem}.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png.relative_to(ROOT)}")
    return png


def banner(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.0,
        1.02,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=MUTED,
    )


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    ext = pd.read_csv(TAB / "ifn_mhc_extended.tsv", sep="\t")
    core = pd.read_csv(TAB / "ifn_mhc_core8.tsv", sep="\t")
    gsea = pd.read_csv(TAB / "gsea_headline.tsv", sep="\t")
    qc = pd.read_csv(TAB / "cldn4_qc.tsv", sep="\t")
    locked = json.loads((TAB / "locked_stats.json").read_text())
    return ext, core, gsea, qc, locked


def fig1_heatmap(ext: pd.DataFrame, locked: dict) -> None:
    """IFN/MHC heatmap of the 39 detected genes (34 up)."""
    det = ext[ext["direction"] != "not detected"].copy()
    det["logFC"] = pd.to_numeric(det["logFC"])
    det = det.sort_values(["arm", "logFC"], ascending=[True, False])
    n = len(det)
    n_up = int((det["logFC"] > 0).sum())
    assert n == 39 and n_up == 34, (n, n_up)

    # Two-column layout: ISG then MHC, each a 1-wide heatmap
    isg = det[det["arm"] == "ISG"].sort_values("logFC", ascending=False)
    mhc = det[det["arm"] == "MHC-I/APM"].sort_values("logFC", ascending=False)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "hc", ["#1A4CFF", "#0C1220", "#FFC857"], N=256
    )
    vmax = 1.9
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.2, 11.4),
        gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.62},
    )
    for ax, block, title, accent in (
        (axes[0], isg, "ISG", IFN),
        (axes[1], mhc, "MHC-I / APM", MHC),
    ):
        mat = block["logFC"].to_numpy().reshape(-1, 1)
        ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks([])
        ax.set_yticks(range(len(block)))
        labels = [
            f"{h}  {m}" if h != m else h
            for h, m in zip(block["human"], block["matched"].fillna(""))
        ]
        ax.set_yticklabels(labels, fontsize=8, fontfamily="monospace")
        ax.set_title(f"{title}  {int((block['logFC']>0).sum())}/{len(block)} up", color=accent, pad=8)
        for i, (fc, core) in enumerate(zip(block["logFC"], block["core"])):
            col = INK if abs(fc) < 0.55 else BG
            star = " ▸" if core == "yes" else ""
            ax.text(0, i, f"{fc:+.2f}{star}", ha="center", va="center", fontsize=7.5, color=col, fontweight="bold")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)

    fig.suptitle(
        f"GSE50927 Cldn4 KO  ·  IFN/MHC panel  ·  {n_up}/{n} genes UP",
        color=INK,
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.012,
        "Naive mouse whole lung (not lung cancer)  ·  n=1 vs 1  ·  author edgeR logFC  ·  ▸ = core 8  ·  no new p",
        ha="center",
        color=MUTED,
        fontsize=8.5,
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.04, shrink=0.55)
    cbar.set_label("logFC (KO − WT)", color=INK)
    cbar.ax.yaxis.set_tick_params(color=INK)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=INK)
    save(fig, "fig_extra1_ifn_mhc_heatmap")


def fig2_nes_bar(gsea: pd.DataFrame, locked: dict) -> None:
    gsea = gsea.copy()
    colors = [POS_NES if n > 0 else NEG_NES for n in gsea["nes"]]
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    y = np.arange(len(gsea))
    ax.barh(y, gsea["nes"], color=colors, edgecolor=INK, linewidth=0.6, height=0.72, zorder=3)
    ax.axvline(0, color=INK, lw=1.0, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(gsea["label"], fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel("NES  (positive = up after Cldn4 loss)")
    ax.set_xlim(-2.05, 2.35)
    ax.xaxis.grid(True, ls=":", color=GRID, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("GSE50927 prerank GSEA  ·  locked NES (not recomputed)")
    banner(ax, "Naive whole lung, not lung cancer  ·  n=1 vs 1  ·  BH-FDR inside 5 headline sets  ·  seed=42")
    for i, row in gsea.iterrows():
        nes = float(row["nes"])
        fdr = float(row["fdr"])
        star = "  ★ FDR<0.05" if fdr < 0.05 else ""
        label = f"{nes:+.2f}   FDR {fdr:.3f}{star}"
        ha = "left" if nes >= 0 else "right"
        x = nes + (0.06 if nes >= 0 else -0.06)
        ax.text(x, i, label, va="center", ha=ha, fontsize=9, color=INK, fontweight="bold")
    ax.text(
        0.99,
        0.02,
        "IFN-γ NES +1.51  (locked)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=IFN,
        fontweight="bold",
    )
    save(fig, "fig_extra2_nes_bar")


def fig3_cldn4_qc(qc: pd.DataFrame, locked: dict) -> None:
    # Expressed claudins + epithelial/junction context
    show = qc[qc["class"].isin(["target", "claudin", "epithelial"])].copy()
    show = show.sort_values("logFC")
    colors = [TARGET if s == "Cldn4" else (UP if fc > 0 else DOWN) for s, fc in zip(show["symbol"], show["logFC"])]
    fig, ax = plt.subplots(figsize=(10.4, 6.4))
    y = np.arange(len(show))
    ax.barh(y, show["logFC"], color=colors, edgecolor=INK, linewidth=0.4, height=0.7, zorder=3)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(show["symbol"], fontsize=11)
    ax.set_xlabel("author edgeR logFC (Cldn4 KO − WT)")
    ax.set_title("Cldn4 QC  ·  target collapses; other claudins do not")
    banner(ax, "GSE50927 naive whole lung (not a tumour KO)  ·  n=1 vs 1  ·  author FDR shown as deposited  ·  no new p")
    ax.xaxis.grid(True, ls=":", color=GRID)
    ax.set_axisbelow(True)
    for i, row in show.reset_index(drop=True).iterrows():
        fc = float(row["logFC"])
        fdr = float(row["author_FDR"])
        if row["symbol"] == "Cldn4":
            txt = f"{fc:+.2f}   author FDR {fdr:.2e}"
            ax.text(fc + 0.12, i, txt, va="center", ha="left", fontsize=9, color=TARGET, fontweight="bold")
        elif fdr < 0.05:
            ax.text(fc + (0.08 if fc >= 0 else -0.08), i, f"{fc:+.2f}", va="center",
                    ha="left" if fc >= 0 else "right", fontsize=8, color=INK)
    ax.set_xlim(-7.4, 1.6)
    ax.text(
        0.98,
        0.06,
        "Cldn4 logFC = −6.06\nlogCPM = 2.25\nauthor FDR = 4.07×10⁻²⁶",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color=TARGET,
        fontweight="bold",
        linespacing=1.45,
    )
    save(fig, "fig_extra3_cldn4_qc")


def fig4_lollipop(ext: pd.DataFrame) -> None:
    det = ext[ext["direction"] != "not detected"].copy()
    det["logFC"] = pd.to_numeric(det["logFC"])
    det = det.sort_values("logFC")
    fig, ax = plt.subplots(figsize=(8.8, 11.0))
    y = np.arange(len(det))
    cols = [UP if v > 0 else DOWN for v in det["logFC"]]
    ax.hlines(y, 0, det["logFC"], color=cols, lw=1.6, zorder=2)
    ax.scatter(det["logFC"], y, c=cols, s=36, zorder=3, edgecolors=INK, linewidths=0.4)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(y)
    labels = [f"{h} ({m})" if m and m != h else h for h, m in zip(det["human"], det["matched"].fillna(""))]
    ax.set_yticklabels(labels, fontsize=8, fontfamily="monospace")
    ax.set_xlabel("author edgeR logFC (KO − WT)")
    ax.set_title("34 / 39 IFN+MHC genes UP after Cldn4 loss")
    banner(ax, "Extended ISG + MHC-I/APM panel  ·  4 human symbols not on the mouse table  ·  no new p")
    ax.xaxis.grid(True, ls=":", color=GRID)
    n_up = int((det["logFC"] > 0).sum())
    n_dn = int((det["logFC"] < 0).sum())
    ax.text(
        0.98,
        0.02,
        f"{n_up} up   {n_dn} down   median +0.43 (locked)",
        transform=ax.transAxes,
        ha="right",
        fontsize=9,
        color=UP,
        fontweight="bold",
    )
    save(fig, "fig_extra4_panel_lollipop")


def fig5_core8(core: pd.DataFrame) -> None:
    core = core.copy()
    core["logFC"] = pd.to_numeric(core["logFC"])
    core = core.sort_values("logFC")
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    y = np.arange(len(core))
    cols = [UP if v > 0 else DOWN for v in core["logFC"]]
    ax.barh(y, core["logFC"], color=cols, edgecolor=INK, linewidth=0.5, height=0.68, zorder=3)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{h}  →  {m}" for h, m in zip(core["human"], core["matched"])], fontsize=11)
    ax.set_xlabel("author edgeR logFC (KO − WT)")
    ax.set_title("Core 8 (user panel)  ·  6 / 8 UP  ·  author FDR only")
    banner(ax, "IFI27 / OAS2 / IFIT1 / MX1 / ISG15 / HLA-A / TAP1 / TAP2  ·  GSE50927 whole lung  ·  no new p")
    ax.xaxis.grid(True, ls=":", color=GRID)
    for i, row in core.reset_index(drop=True).iterrows():
        fc = float(row["logFC"])
        fdr = float(row["author_FDR"])
        mark = "  ★ author FDR<0.05" if fdr < 0.05 else ""
        ax.text(
            fc + (0.05 if fc >= 0 else -0.05),
            i,
            f"{fc:+.2f}{mark}",
            va="center",
            ha="left" if fc >= 0 else "right",
            fontsize=9,
            color=INK,
            fontweight="bold",
        )
    save(fig, "fig_extra5_core8")


def fig6_qc_cpm(qc: pd.DataFrame) -> None:
    """Second QC: Cldn4 vs housekeeping / epithelial abundance (logCPM) + logFC."""
    pick = qc[qc["symbol"].isin(["Cldn4", "Tacstd2", "Epcam", "Cldn3", "Cldn7", "Actb", "Tjp1", "Isg15", "B2m", "Nlrc5"])]
    # Isg15/B2m/Nlrc5 are in the deposited QC extract
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    x = pick["logCPM"].to_numpy()
    y = pick["logFC"].to_numpy()
    names = pick["symbol"].to_numpy()
    for xi, yi, name in zip(x, y, names):
        col = TARGET if name == "Cldn4" else (UP if yi > 0.3 else (DOWN if yi < 0 else MUTED))
        ax.scatter([xi], [yi], s=160 if name == "Cldn4" else 90, c=col, edgecolors=INK, linewidths=0.6, zorder=3)
        ax.annotate(name, (xi, yi), textcoords="offset points", xytext=(7, 5), fontsize=9,
                    color=TARGET if name == "Cldn4" else INK, fontweight="bold" if name == "Cldn4" else "regular")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xlabel("author logCPM (abundance)")
    ax.set_ylabel("author logFC (KO − WT)")
    ax.set_title("Cldn4 QC  ·  detected, collapsed, not a low-count artefact")
    banner(ax, "Same deposited edgeR table  ·  Isg15 / B2m / Nlrc5 shown as IFN/MHC context  ·  no new p")
    ax.grid(True, ls=":", color=GRID)
    save(fig, "fig_extra6_cldn4_qc_cpm")


def main() -> None:
    style()
    ext, core, gsea, qc, locked = load()
    fig1_heatmap(ext, locked)
    fig2_nes_bar(gsea, locked)
    fig3_cldn4_qc(qc, locked)
    fig4_lollipop(ext)
    fig5_core8(core)
    fig6_qc_cpm(qc)
    n = len(list(FIG.glob("fig_extra*.png")))
    print(f"extra PNG count: {n}")
    if n < 4:
        raise SystemExit(f"need ≥4 extra figures, got {n}")


if __name__ == "__main__":
    main()
