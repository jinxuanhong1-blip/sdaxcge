#!/usr/bin/env python3
"""Redraw already-significant public CLDN4 extras as a high-contrast slide pack.

No re-audit. No new tests. Numbers are locked from prior public extras
(PRs 239, 251, 254, 287, 289, 290, 292). Fisher-z intervals used only as
display whiskers from published n and ρ — they are not new p-values.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
TABLE = HERE / "locked_numbers.tsv"
FIGDIR = HERE / "figures"

BG = "#0B1220"
PANEL = "#121A2B"
INK = "#F4F7FB"
MUTED = "#A8B3C7"
GRID = "#2A3550"
GOLD = "#F5C542"
CYAN = "#4CC9F0"
MINT = "#3DDC97"
CORAL = "#FF4D6D"
AMBER = "#FFB703"
SLATE = "#7D8CA8"

W, H = 13.333, 7.5


def fmt_p(p: float) -> str:
    if p < 1e-4:
        return f"{p:.1e}"
    if p < 0.01:
        return f"{p:.3g}"
    return f"{p:.2f}"


def fmt_rho(r: float) -> str:
    sign = "+" if r > 0 else "−"
    return f"{sign}{abs(r):.3f}"


def apply_slide(fig) -> None:
    fig.patch.set_facecolor(BG)


def panel_ax(ax, title: str | None = None) -> None:
    ax.set_facecolor(PANEL)
    for sp in ax.spines.values():
        sp.set_color(GRID)
        sp.set_linewidth(1.0)
    ax.tick_params(colors=INK, labelsize=11)
    ax.yaxis.label.set_color(INK)
    ax.xaxis.label.set_color(INK)
    if title:
        ax.set_title(title, color=INK, fontsize=14, pad=10, loc="left", fontweight="semibold")


def footer(fig, text: str) -> None:
    fig.text(0.018, 0.018, text, color=SLATE, fontsize=8.2, va="bottom", ha="left")


def header(fig, kicker: str, title: str) -> None:
    fig.text(0.03, 0.955, kicker.upper(), color=GOLD, fontsize=10.5, fontweight="bold", va="top")
    fig.text(0.03, 0.905, title, color=INK, fontsize=20, fontweight="bold", va="top")


def load_locked() -> pd.DataFrame:
    df = pd.read_csv(TABLE, sep="\t")
    for col in ("n", "effect", "effect_lo", "effect_hi", "p"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def save(fig, stem: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{stem}.png", dpi=220, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.pdf", facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def fig01_forest(df: pd.DataFrame) -> None:
    rows = df[df["figure_id"] == "forest"].copy().reset_index(drop=True)
    # Display order: immune-low (negative) first, then IFN-positive extras
    order = [
        "CLDN4 vs T+B all AOI",
        "author %pos CLDN4 vs T/NK",
        "CLDN4 protein vs GEP18 RNA",
        "CLDN4 protein vs CD8A RNA",
        "CLDN4 protein vs ImmuneScore",
        "CLDN4 vs CXCL13+ T fraction",
        "CLDN4 vs T+B tumor AOI",
        "CLDN4 vs CD8A",
        "CLDN4 vs IFN-compact",
        "CLDN4 vs MHC-I",
    ]
    rows["_ord"] = rows["contrast"].map({k: i for i, k in enumerate(order)})
    rows = rows.sort_values("_ord", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.42, right=0.78, top=0.78, bottom=0.12)
    panel_ax(ax)
    header(
        fig,
        "Extra figure 1  ·  one forest",
        "Public CLDN4 extras — Spearman ρ with published n / p",
    )

    y = np.arange(len(rows))
    for i, r in rows.iterrows():
        color = MINT if r["effect"] > 0 else CORAL
        ax.plot(
            [r["effect_lo"], r["effect_hi"]],
            [i, i],
            color=color,
            lw=3.2,
            solid_capstyle="round",
            zorder=2,
        )
        ax.scatter(
            [r["effect"]],
            [i],
            s=70,
            color=color,
            edgecolor=INK,
            linewidth=0.8,
            zorder=3,
        )
        left = f"{r['accession']}\n{r['contrast']}"
        ax.text(-1.18, i, left, color=INK, fontsize=10.2, va="center", ha="left", linespacing=1.25)
        right = f"n={int(r['n'])}   ρ={fmt_rho(r['effect'])}   p={fmt_p(r['p'])}"
        ax.text(1.12, i, right, color=INK, fontsize=10.5, va="center", ha="left", fontweight="medium")

    ax.axvline(0, color=GOLD, lw=1.2, ls="--", zorder=1)
    ax.set_xlim(-1.20, 1.10)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("Spearman ρ  (whiskers = published CI, or Fisher-z display from n, ρ)", color=MUTED, fontsize=10)
    ax.set_xticks([-1.0, -0.5, 0, 0.5, 1.0])
    footer(
        fig,
        "Locked from PRs 239 / 251 / 289 / 290 / 292. CPTAC and GSE285029 CIs are published. "
        "GSE218989, PR290, GSE265899 whiskers are Fisher-z 95% display intervals from the published n and ρ — not a new test. "
        "No p-value was invented.",
    )
    save(fig, "fig01_forest_spearman")


def fig02_gse50927(df: pd.DataFrame) -> None:
    rows = df[(df["figure_id"] == "gse50927") & (df["metric"] == "NES")].copy()
    log2 = df[(df["figure_id"] == "gse50927") & (df["metric"] == "cldn4_log2FC")].iloc[0]
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.22, right=0.92, top=0.72, bottom=0.16)
    panel_ax(ax)
    header(fig, "Extra figure 2  ·  GSE50927 IFN NES", "Cldn4 loss in naive whole lung lifts IFN / MHC-I NES")

    labels = ["Hallmark IFN-γ", "Hallmark IFN-α", "MHC-I / APM"]
    nes = rows["effect"].to_numpy()
    fdr = ["FDR 0.005", "FDR 0.010", "FDR 0.050"]
    colors = [MINT, CYAN, GOLD]
    bars = ax.barh(labels[::-1], nes[::-1], color=colors[::-1], height=0.55, edgecolor=INK, linewidth=0.6)
    for bar, val, lab in zip(bars, nes[::-1], fdr[::-1]):
        ax.text(
            val + 0.04,
            bar.get_y() + bar.get_height() / 2,
            f"NES {val:+.3f}   {lab}",
            color=INK,
            va="center",
            fontsize=13,
            fontweight="medium",
        )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(0, 2.15)
    ax.set_xlabel("NES after Cldn4 loss  (positive = up on the loss rank)", color=MUTED)
    ax.text(
        0.0,
        -0.95,
        f"GSE50927 naive whole lung  ·  n=1 vs 1  ·  Cldn4 log2FC {log2['display_effect']}  ·  author FDR {fmt_p(log2['p'])}\n"
        "Not a tumour KO. Mouse mixed-cell bulk. Prerank 1000 permutations, seed=42. BH-FDR inside five headline sets.",
        color=MUTED,
        fontsize=11,
        transform=ax.get_xaxis_transform(),
        va="top",
    )
    footer(fig, "Source PR 287  ·  methods/cldn4_ko_gsea/tables/gsea_headline.tsv  ·  honest NES / FDR only.")
    save(fig, "fig02_gse50927_ifn_nes")


def fig03_cptac(df: pd.DataFrame) -> None:
    rows = df[(df["figure_id"] == "forest") & (df["accession"] == "CPTAC-LSCC")].copy()
    part = df[df["metric"] == "partial_spearman"].iloc[0]
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.28, right=0.88, top=0.74, bottom=0.16)
    panel_ax(ax)
    header(fig, "Extra figure 3  ·  CPTAC LSCC protein −0.43", "CLDN4 protein is inverse with ImmuneScore / GEP18 / CD8A")

    rows = rows.iloc[::-1]
    y = np.arange(len(rows))
    ax.barh(
        y,
        rows["effect"],
        color=CORAL,
        height=0.52,
        edgecolor=INK,
        linewidth=0.5,
        xerr=[rows["effect"] - rows["effect_lo"], rows["effect_hi"] - rows["effect"]],
        error_kw={"ecolor": INK, "elinewidth": 1.4, "capsize": 4},
    )
    ax.set_yticks(y)
    ax.set_yticklabels(rows["contrast"], color=INK, fontsize=12)
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(-0.85, 0.15)
    for i, r in enumerate(rows.itertuples()):
        ax.text(
            0.02,
            i,
            f"n={int(r.n)}   ρ={fmt_rho(r.effect)}   p={fmt_p(r.p)}",
            color=INK,
            va="center",
            fontsize=12,
        )
    ax.set_xlabel("Spearman ρ  (published bootstrap 95% CI)", color=MUTED)
    ax.text(
        0.0,
        -0.85,
        f"n=78 / 108 LSCC tumors with CLDN4 protein. Partial ρ | WES purity vs ImmuneScore = {fmt_rho(part['effect'])} "
        f"(n=77, p={fmt_p(part['p'])}). TJ-15 protein vs ImmuneScore is null (ρ=−0.082, n=108, p=0.40) and is not drawn.",
        color=MUTED,
        fontsize=11,
        transform=ax.get_xaxis_transform(),
        va="top",
    )
    footer(fig, "Source PR 289  ·  methods/cptac_lusc_cldn4_protein/tables/associations.tsv  ·  no invented p.")
    save(fig, "fig03_cptac_lscc_protein")


def fig04_gse218989(df: pd.DataFrame) -> None:
    rho = df[(df["accession"] == "GSE218989") & (df["metric"] == "spearman_rho")].iloc[0]
    mwu = df[(df["accession"] == "GSE218989") & (df["metric"] == "mwu_auc")].iloc[0]
    cox = df[(df["accession"] == "GSE218989") & (df["metric"] == "hr")].iloc[0]
    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(fig, "Extra figure 4  ·  GSE218989 CLDN4 vs CD8", "CLDN4 tracks CD8-low, not ICI objective response")

    cards = [
        (f"ρ = {fmt_rho(rho['effect'])}", f"CLDN4 vs CD8A\nn={int(rho['n'])} patients\np={fmt_p(rho['p'])}", CORAL),
        (f"AUC = {mwu['effect']:.2f}", f"CLDN4 vs GEO response\n168 R / 187 NR\nMWU p={fmt_p(mwu['p'])}  (null)", SLATE),
        (f"HR = {cox['effect']:.2f}", f"CLDN4 OS per SD\n265 deaths / n=355\n95% CI {cox['effect_lo']:.2f}–{cox['effect_hi']:.2f}  p={fmt_p(cox['p'])}", AMBER),
    ]
    for i, (big, sub, color) in enumerate(cards):
        x0 = 0.07 + i * 0.31
        box = FancyBboxPatch(
            (x0, 0.28),
            0.28,
            0.48,
            boxstyle="round,pad=0.018,rounding_size=0.02",
            transform=fig.transFigure,
            facecolor=PANEL,
            edgecolor=GRID,
            linewidth=1.4,
        )
        fig.patches.append(box)
        fig.text(x0 + 0.14, 0.62, big, color=color, fontsize=28, ha="center", va="center", fontweight="bold")
        fig.text(x0 + 0.14, 0.42, sub, color=INK, fontsize=13, ha="center", va="center", linespacing=1.45)
    fig.text(
        0.03,
        0.14,
        "Patient is the unit. GEO GSE218989 SMC/KAIST TPM (19,916 genes × 355). Histology is not deposited. "
        "Response is null; the CD8 inverse is the extra that is redrawn here.",
        color=MUTED,
        fontsize=11,
    )
    footer(fig, "Source PR 292  ·  methods/gse218989_cldn4_ici/tables/stats.tsv  ·  honest n / ρ / p only.")
    save(fig, "fig04_gse218989_cd8")


def fig05_pr290(df: pd.DataFrame) -> None:
    rows = df[(df["source_pr"] == 290) & (df["figure_id"] == "forest")].copy()
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.12, right=0.92, top=0.72, bottom=0.18)
    panel_ax(ax)
    header(fig, "Extra figure 5  ·  PR290 combos n=43 / −0.48 and CXCL13 n=60 / −0.43", "CLDN4-first scRNA combinations that recover ρ < 0")

    labels = [
        "Author %pos vs T/NK\nGSE131907 + GSE205335",
        "CXCL13+ T fraction\nGSE148071 + GSE207422 + GSE253013",
    ]
    vals = rows["effect"].to_numpy()
    ns = rows["n"].to_numpy()
    ps = rows["p"].to_numpy()
    x = np.arange(2)
    ax.bar(x, vals, color=[CORAL, AMBER], width=0.55, edgecolor=INK, linewidth=0.6)
    ax.axhline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK, fontsize=12)
    ax.set_ylim(-0.75, 0.12)
    ax.set_ylabel("Spearman ρ (CLDN4 vs immune)", color=INK)
    for i, (v, n, p) in enumerate(zip(vals, ns, ps)):
        ax.text(i, v - 0.06, f"n={int(n)}\nρ={fmt_rho(v)}\np={fmt_p(p)}\nI²=0%", ha="center", va="top", color=INK, fontsize=13)
    ax.text(
        0.5,
        -0.92,
        "Full-pool primary mixed grid is weaker (N=145, ρ=−0.137, p=0.129) and is not claimed as the answer. "
        "These two highlighted combos are the extras named in the request.",
        color=MUTED,
        fontsize=11,
        ha="center",
        transform=ax.get_xaxis_transform(),
    )
    footer(fig, "Source PR 290  ·  methods/scrna_cldn4_combo/tables/highlighted_combos.tsv  ·  no invented p.")
    save(fig, "fig05_pr290_combos")


def fig06_gse285029(df: pd.DataFrame) -> None:
    rows = df[(df["accession"] == "GSE285029") & (df["figure_id"] == "forest")].copy()
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.16, right=0.90, top=0.72, bottom=0.16)
    panel_ax(ax)
    header(fig, "Extra figure 6  ·  GSE285029 CLDN4 vs IFN", "CLDN4-high sits on IFN / MHC-I, not on CD8")

    labels = ["CLDN4 vs IFN-compact", "CLDN4 vs MHC-I"]
    vals = rows["effect"].to_numpy()
    lo = rows["effect_lo"].to_numpy()
    hi = rows["effect_hi"].to_numpy()
    ns = rows["n"].to_numpy()
    ps = rows["p"].to_numpy()
    x = np.arange(2)
    ax.bar(x, vals, color=[CYAN, MINT], width=0.5, edgecolor=INK, linewidth=0.6)
    ax.errorbar(x, vals, yerr=[vals - lo, hi - vals], fmt="none", ecolor=INK, elinewidth=1.6, capsize=5)
    ax.axhline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK, fontsize=13)
    ax.set_ylim(-0.05, 0.55)
    ax.set_ylabel("Spearman ρ  (published table CI)", color=INK)
    for i, (v, n, p) in enumerate(zip(vals, ns, ps)):
        ax.text(i, v + 0.045, f"n={int(n)}   ρ={fmt_rho(v)}   p={fmt_p(p)}", ha="center", color=INK, fontsize=13)
    ax.text(
        0.5,
        -0.18,
        "n=234 pre-ICI NSCLC WTS (Koh 2025). CLDN4 vs CD8A is null (ρ=+0.067, p=0.31) and is not a claimed extra. "
        "Epithelial residual drops IFN-compact to partial p≈0.12 — shown only as the published caveat.",
        color=MUTED,
        fontsize=11,
        ha="center",
        transform=ax.get_xaxis_transform(),
    )
    footer(fig, "Source PR 239  ·  results/w200/A11_GSE285029/correlations.csv  ·  honest n / ρ / p only.")
    save(fig, "fig06_gse285029_ifn")


def fig07_spatial(df: pd.DataFrame) -> None:
    rows = df[(df["accession"] == "GSE265899") & (df["figure_id"] == "forest")].copy()
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.14, right=0.90, top=0.72, bottom=0.16)
    panel_ax(ax)
    header(fig, "Extra figure 7  ·  spatial GSE265899", "GeoMx CLDN4 vs T+B neighborhood")

    labels = ["All AOI\nn=95", "Tumor AOI\nn=48"]
    vals = rows["effect"].to_numpy()
    ps = rows["p"].to_numpy()
    x = np.arange(2)
    ax.bar(x, vals, color=[CORAL, AMBER], width=0.5, edgecolor=INK, linewidth=0.6)
    ax.axhline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK, fontsize=14)
    ax.set_ylim(-1.05, 0.12)
    ax.set_ylabel("Spearman ρ (CLDN4 vs T+B)", color=INK)
    for i, (v, p) in enumerate(zip(vals, ps)):
        ax.text(i, v - 0.08, f"ρ={fmt_rho(v)}\np={fmt_p(p)}", ha="center", va="top", color=INK, fontsize=14)
    ax.text(
        0.5,
        -0.95,
        "Immune AOIs in the same series are NS for CLDN4 vs T+B (ρ=−0.180, n=47, p=0.23) and are not drawn as a hit. "
        "All-AOI ρ mixes compartment; tumor-AOI is the cleaner extra.",
        color=MUTED,
        fontsize=11,
        ha="center",
        transform=ax.get_xaxis_transform(),
    )
    footer(fig, "Source PR 251  ·  results/leftover_spatial/tables/geomx_correlations.csv  ·  no invented p.")
    save(fig, "fig07_gse265899_spatial")


def fig08_gse312098(df: pd.DataFrame) -> None:
    rows = df[df["figure_id"] == "gse312098"].copy()
    cldn = rows[rows["metric"] == "log2FC"].iloc[0]
    sets = rows[rows["metric"] == "mwu_direction"]
    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(fig, "Extra figure 8  ·  GSE312098 CLDN4 down", "CX-1 IMMU132 2 d: CLDN4 down, IFN / APM up")

    box = FancyBboxPatch(
        (0.07, 0.30),
        0.36,
        0.48,
        boxstyle="round,pad=0.018,rounding_size=0.02",
        transform=fig.transFigure,
        facecolor=PANEL,
        edgecolor=GRID,
        linewidth=1.4,
    )
    fig.patches.append(box)
    fig.text(0.25, 0.64, "CLDN4 log2FC", color=MUTED, ha="center", fontsize=13)
    fig.text(0.25, 0.54, "−0.86", color=CORAL, ha="center", fontsize=42, fontweight="bold")
    fig.text(
        0.25,
        0.40,
        f"n=3 vs 3 FPKM\np={fmt_p(cldn['p'])}   q=0.0043\nCX-1 CRC, IMMU132 2 d",
        color=INK,
        ha="center",
        fontsize=13,
        linespacing=1.45,
    )

    ax = fig.add_axes([0.50, 0.28, 0.44, 0.50])
    panel_ax(ax, "Gene-set direction after IMMU132 (published MW)")
    labels = ["IFN set\n55 genes", "APM\n19 genes", "C4 IFN/MHC-I\n6 genes"]
    ps = sets["p"].to_numpy()
    y = np.arange(3)
    ax.barh(y, -np.log10(ps), color=[MINT, CYAN, GOLD], height=0.55, edgecolor=INK, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK, fontsize=11)
    ax.set_xlabel("−log10 p  (published MW on log2FC)", color=MUTED)
    for i, (p, lab) in enumerate(zip(ps, sets["display_effect"])):
        ax.text(-np.log10(p) + 0.08, i, f"{lab}   p={fmt_p(p)}", color=INK, va="center", fontsize=11)
    ax.set_xlim(0, 9.2)
    fig.text(
        0.03,
        0.14,
        "Taken as given from PR 254 (not re-audited). Not SKB264. Not lung. Junction set is null (MW p=0.99) and is not drawn as a hit.",
        color=MUTED,
        fontsize=11,
    )
    footer(fig, "Source PR 254  ·  results/c_public_trop2_adc_analogs/tables  ·  honest log2FC / p / q only.")
    save(fig, "fig08_gse312098_cldn4_down")


def _mini_bar(ax, labels, values, colors, annots, xlim, xlabel) -> None:
    panel_ax(ax)
    y = np.arange(len(values))
    ax.barh(y, values, color=colors, height=0.62, edgecolor=INK, linewidth=0.4)
    ax.axvline(0, color=GOLD, lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK, fontsize=8.6)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, color=MUTED, fontsize=8)
    ax.tick_params(labelsize=8)
    for i, (v, t) in enumerate(zip(values, annots)):
        ha = "left" if v >= 0 else "right"
        off = 0.03 if v >= 0 else -0.03
        ax.text(v + off, i, t, color=INK, va="center", ha=ha, fontsize=8.2)


def fig09_multipanel(df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(fig, "Extra figure set  ·  one multi-panel", "Seven already-significant public CLDN4 extras")
    gs = fig.add_gridspec(2, 4, left=0.05, right=0.98, top=0.80, bottom=0.10, wspace=0.55, hspace=0.55)

    # A GSE50927 NES
    ax = fig.add_subplot(gs[0, 0])
    nes = df[(df["figure_id"] == "gse50927") & (df["metric"] == "NES")]
    _mini_bar(
        ax,
        ["IFN-γ", "IFN-α", "MHC-I"],
        nes["effect"].to_list(),
        [MINT, CYAN, GOLD],
        [f"NES {v:+.2f}\n{d}" for v, d in zip(nes["effect"], nes["display_p"])],
        (0, 2.4),
        "GSE50927 NES  n=1 vs 1",
    )
    ax.set_title("A  GSE50927 IFN NES", color=GOLD, fontsize=10, loc="left")

    # B CPTAC
    ax = fig.add_subplot(gs[0, 1])
    cpt = df[(df["figure_id"] == "forest") & (df["accession"] == "CPTAC-LSCC")]
    _mini_bar(
        ax,
        ["ImmuneScore", "GEP18", "CD8A"],
        cpt["effect"].to_list(),
        [CORAL] * 3,
        [f"ρ={fmt_rho(v)}  n=78\np={fmt_p(p)}" for v, p in zip(cpt["effect"], cpt["p"])],
        (-0.75, 0.15),
        "CPTAC LSCC protein",
    )
    ax.set_title("B  CPTAC protein −0.43", color=GOLD, fontsize=10, loc="left")

    # C GSE218989
    ax = fig.add_subplot(gs[0, 2])
    r = df[(df["accession"] == "GSE218989") & (df["metric"] == "spearman_rho")].iloc[0]
    _mini_bar(
        ax,
        ["CLDN4 vs CD8A"],
        [r["effect"]],
        [CORAL],
        [f"n=355  ρ={fmt_rho(r['effect'])}\np={fmt_p(r['p'])}"],
        (-0.45, 0.12),
        "GSE218989 patient TPM",
    )
    ax.set_title("C  GSE218989 vs CD8", color=GOLD, fontsize=10, loc="left")

    # D PR290
    ax = fig.add_subplot(gs[0, 3])
    pr = df[(df["source_pr"] == 290) & (df["figure_id"] == "forest")]
    _mini_bar(
        ax,
        ["T/NK combo", "CXCL13+"],
        pr["effect"].to_list(),
        [CORAL, AMBER],
        [f"n={int(n)}  ρ={fmt_rho(v)}\np={fmt_p(p)}" for n, v, p in zip(pr["n"], pr["effect"], pr["p"])],
        (-0.75, 0.12),
        "PR290 highlighted combos",
    )
    ax.set_title("D  PR290 n=43 / n=60", color=GOLD, fontsize=10, loc="left")

    # E GSE285029
    ax = fig.add_subplot(gs[1, 0])
    g = df[(df["accession"] == "GSE285029") & (df["figure_id"] == "forest")]
    _mini_bar(
        ax,
        ["IFN-compact", "MHC-I"],
        g["effect"].to_list(),
        [CYAN, MINT],
        [f"n=234  ρ={fmt_rho(v)}\np={fmt_p(p)}" for v, p in zip(g["effect"], g["p"])],
        (-0.05, 0.45),
        "GSE285029 WTS",
    )
    ax.set_title("E  GSE285029 vs IFN", color=GOLD, fontsize=10, loc="left")

    # F spatial
    ax = fig.add_subplot(gs[1, 1])
    s = df[(df["accession"] == "GSE265899") & (df["figure_id"] == "forest")]
    _mini_bar(
        ax,
        ["All AOI", "Tumor AOI"],
        s["effect"].to_list(),
        [CORAL, AMBER],
        [f"n={int(n)}  ρ={fmt_rho(v)}\np={fmt_p(p)}" for n, v, p in zip(s["n"], s["effect"], s["p"])],
        (-1.05, 0.15),
        "GSE265899 GeoMx",
    )
    ax.set_title("F  spatial GSE265899", color=GOLD, fontsize=10, loc="left")

    # G GSE312098
    ax = fig.add_subplot(gs[1, 2])
    cldn = df[(df["figure_id"] == "gse312098") & (df["metric"] == "log2FC")].iloc[0]
    _mini_bar(
        ax,
        ["CLDN4"],
        [cldn["effect"]],
        [CORAL],
        [f"n=3 vs 3  log2FC −0.86\np={fmt_p(cldn['p'])}"],
        (-1.2, 0.15),
        "GSE312098 CX-1 2 d",
    )
    ax.set_title("G  GSE312098 CLDN4 down", color=GOLD, fontsize=10, loc="left")

    # H note card
    ax = fig.add_subplot(gs[1, 3])
    ax.set_facecolor(PANEL)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_title("H  What is not claimed", color=GOLD, fontsize=10, loc="left")
    ax.text(
        0.06,
        0.88,
        "GSE50927 is mouse whole lung,\nnot a tumour KO.\n\n"
        "GSE218989 response is null\n(AUC 0.47, p=0.40).\n\n"
        "GSE312098 is IMMU132 CRC,\nnot SKB264 / not lung.\n\n"
        "No p-value was invented.",
        color=INK,
        fontsize=9.2,
        va="top",
        ha="left",
        transform=ax.transAxes,
        linespacing=1.35,
    )

    footer(fig, "Additive redraw only. Locked numbers in locked_numbers.tsv. Sources: PRs 239, 251, 254, 287, 289, 290, 292.")
    save(fig, "fig09_multipanel_extras")


def main() -> None:
    df = load_locked()
    fig01_forest(df)
    fig02_gse50927(df)
    fig03_cptac(df)
    fig04_gse218989(df)
    fig05_pr290(df)
    fig06_gse285029(df)
    fig07_spatial(df)
    fig08_gse312098(df)
    fig09_multipanel(df)
    print(f"wrote {len(list(FIGDIR.glob('*.png')))} png under {FIGDIR}")


if __name__ == "__main__":
    main()
