#!/usr/bin/env python3
"""PAPER FUNNEL: PPT-ready corroboration pack for TJ/CLDN4 after Tacstd2 DEG.

Redraw only. No re-download. No new tests. Every number is locked from
PRs 539 / 579 / 289 / 693 / 685 / 698 and the public CosMx handoff.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

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
CLDN = "#FF6B4A"
TJ = "#7B8CDE"

W, H = 13.333, 7.5


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "—"
    if p < 1e-4:
        return f"{p:.1e}"
    if p < 0.01:
        return f"{p:.3g}"
    if p < 0.05:
        return f"{p:.3f}".rstrip("0").rstrip(".")
    return f"{p:.2f}"


def fmt_rho(r: float) -> str:
    if pd.isna(r):
        return "—"
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
    fig.text(0.018, 0.018, text, color=SLATE, fontsize=8.0, va="bottom", ha="left")


def header(fig, kicker: str, title: str) -> None:
    fig.text(0.03, 0.955, kicker.upper(), color=GOLD, fontsize=10.5, fontweight="bold", va="top")
    fig.text(0.03, 0.905, title, color=INK, fontsize=19.5, fontweight="bold", va="top")


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


def card(fig, x, y, w, h, big, sub, color=GOLD):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.016",
        transform=fig.transFigure,
        facecolor=PANEL,
        edgecolor=GRID,
        linewidth=1.3,
    )
    fig.patches.append(box)
    fig.text(x + w / 2, y + h * 0.62, big, color=color, fontsize=22, ha="center", va="center", fontweight="bold")
    fig.text(x + w / 2, y + h * 0.28, sub, color=INK, fontsize=11.5, ha="center", va="center", linespacing=1.35)


def fig01_funnel(df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(
        fig,
        "Paper funnel  ·  slide 1",
        "After Tacstd2 DEG: corroborate TJ / CLDN4 as the immune-inverse feature",
    )

    stages = [
        ("1  ENTRY", "GSE137244 KL>KP\nTacstd2 +3.24\nCldn4 +5.57", GOLD, "#685"),
        ("2  PATIENT", "Concordant-4\nCLDN4–T/NK\nρ = −0.531", CORAL, "#539"),
        ("3  SPACE", "CosMx exclusion\n0.36 / 0.52\nnot muzzling", CYAN, "#698"),
        ("4  BARRIER", "LIANA barrier\nfamily Δ > 0\nq ≤ 4.9e−8", AMBER, "#579"),
        ("5  PROTEIN", "CPTAC LSCC\nCLDN4 −0.43\nTJ-15 null", MINT, "#289"),
        ("6  SIGNATURE", "Sig #693\nssGSEA −0.533\nI² = 0%", CLDN, "#693"),
    ]
    for i, (kicker, body, color, pr) in enumerate(stages):
        x = 0.035 + i * 0.162
        box = FancyBboxPatch(
            (x, 0.38),
            0.148,
            0.42,
            boxstyle="round,pad=0.01,rounding_size=0.014",
            transform=fig.transFigure,
            facecolor=PANEL,
            edgecolor=color,
            linewidth=1.8,
        )
        fig.patches.append(box)
        fig.text(x + 0.074, 0.74, kicker, color=color, fontsize=11, ha="center", fontweight="bold")
        fig.text(x + 0.074, 0.58, body, color=INK, fontsize=12.2, ha="center", va="center", linespacing=1.4)
        fig.text(x + 0.074, 0.42, f"PR {pr}", color=SLATE, fontsize=10, ha="center")
        if i < len(stages) - 1:
            fig.patches.append(
                FancyArrowPatch(
                    (x + 0.148, 0.59),
                    (x + 0.162, 0.59),
                    transform=fig.transFigure,
                    arrowstyle="-|>",
                    mutation_scale=12,
                    color=GOLD,
                    lw=1.4,
                )
            )

    fig.text(
        0.03,
        0.28,
        "Thesis for the pack: prioritize CLDN4 within TJ.",
        color=GOLD,
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.03,
        0.14,
        "CPTAC LSCC: CLDN4 protein vs ImmuneScore ρ=−0.432 (n=78); TJ-15 protein vs ImmuneScore ρ=−0.082 (n=108, p=0.40).\n"
        "No new analysis. Numbers locked from the cited PRs. Private 8-KL / KD co-culture are out of this public pack.",
        color=MUTED,
        fontsize=11.5,
        linespacing=1.45,
    )
    footer(fig, "Redraw pack  ·  methods/paper_funnel_tj_cldn4_corroboration/locked_numbers.tsv  ·  never fabricate.")
    save(fig, "fig01_funnel_overview")


def fig02_entry_deg(df: pd.DataFrame) -> None:
    rows = df[df["figure_id"] == "entry_deg"].copy()
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.22, right=0.90, top=0.72, bottom=0.18)
    panel_ax(ax)
    header(fig, "Paper funnel  ·  slide 2  ·  entry DEG", "GSE137244 KL vs KP: Tacstd2 and Cldn4 rise together")

    order = [
        "Tacstd2 KL vs KP",
        "Cldn4 KL vs KP",
        "TJ_TISMO 7-gene KL vs KP",
        "TJ_EPITHELIAL 18-gene KL vs KP",
    ]
    rows["_ord"] = rows["contrast"].map({k: i for i, k in enumerate(order)})
    rows = rows.sort_values("_ord").reset_index(drop=True)
    colors = [GOLD, CLDN, TJ, CYAN]
    y = np.arange(len(rows))
    ax.barh(y, rows["effect"], color=colors, height=0.55, edgecolor=INK, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(
        ["Tacstd2", "Cldn4", "TJ TISMO (7)", "TJ epithelial (18)"],
        color=INK,
        fontsize=13,
    )
    for i, r in enumerate(rows.itertuples()):
        ax.text(
            r.effect + 0.08,
            i,
            f"Δ = {fmt_rho(r.effect).replace('+','+')}   MW p={fmt_p(r.p)}   complete KL>KP",
            color=INK,
            va="center",
            fontsize=12,
        )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(0, 7.2)
    ax.set_xlabel("Δ mean log2(FPKM+1), KL − KP", color=MUTED)
    ax.text(
        0.0,
        -0.9,
        "Honest n = 5 libraries vs 5 libraries. Exact MW p = 0.00794 is the 5-vs-5 floor.\n"
        "Handoff TJ +3.03 is a different average; the TISMO 7-gene module here is +3.269 (PR 685).",
        color=MUTED,
        fontsize=11,
        transform=ax.get_xaxis_transform(),
        va="top",
    )
    footer(fig, "Source PR 685  ·  methods/gse137244_transfer_signatures/tables/module_contrasts.tsv  ·  locked Tacstd2/Cldn4.")
    save(fig, "fig02_gse137244_entry")


def fig03_concordant4(df: pd.DataFrame) -> None:
    meta = df[(df["figure_id"] == "concordant4") & (df["metric"] == "spearman_meta_rho")].iloc[0]
    singles = df[(df["figure_id"] == "concordant4") & (df["metric"] == "spearman_rho")].copy()
    q = df[(df["figure_id"] == "concordant4") & (df["metric"] == "rank_biserial")].iloc[0]

    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.28, right=0.72, top=0.74, bottom=0.14)
    panel_ax(ax)
    header(
        fig,
        "Paper funnel  ·  slide 3  ·  locked ρ = −0.531",
        "Concordant-4: malignant CLDN4 %pos vs patient T/NK fraction",
    )

    order = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
    singles["_ord"] = singles["accession"].map({k: i for i, k in enumerate(order)})
    singles = singles.sort_values("_ord", ascending=False).reset_index(drop=True)
    y = np.arange(len(singles))
    for i, r in singles.iterrows():
        ax.scatter([r["effect"]], [i], s=90, color=CORAL, edgecolor=INK, zorder=3)
        ax.text(
            0.05,
            i,
            f"n={int(r['n'])}   ρ={fmt_rho(r['effect'])}   p={fmt_p(r['p'])}",
            color=INK,
            va="center",
            fontsize=12,
        )
    # meta diamond
    ax.axhspan(-1.15, -0.45, color="#1A2438", zorder=0)
    ax.scatter([meta["effect"]], [-0.8], s=160, marker="D", color=GOLD, edgecolor=INK, zorder=4)
    ax.plot([meta["effect_lo"], meta["effect_hi"]], [-0.8, -0.8], color=GOLD, lw=3.2, zorder=3)
    ax.text(
        0.05,
        -0.8,
        f"DL meta n=65   ρ={fmt_rho(meta['effect'])}   p={fmt_p(meta['p'])}   I²=0%",
        color=GOLD,
        va="center",
        fontsize=12.5,
        fontweight="bold",
    )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--", zorder=1)
    ax.set_xlim(-1.05, 0.35)
    ax.set_ylim(-1.5, len(singles) - 0.3)
    ax.set_yticks(list(y) + [-0.8])
    ax.set_yticklabels(list(singles["accession"]) + ["META"], color=INK, fontsize=12)
    ax.set_xlabel("Spearman ρ  (meta whiskers = published 95% CI from PR 539)", color=MUTED)
    fig.text(
        0.03,
        0.065,
        f"Stacked Q4 vs Q1 rank-biserial r = {fmt_rho(q['effect'])} (19/16, p={fmt_p(q['p'])}). "
        "Four cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.",
        color=MUTED,
        fontsize=10.5,
    )
    footer(fig, "Source PR 539  ·  methods/seurat_concordant4_cldn4/FINDING.md  ·  patient/donor/sample is the unit.")
    save(fig, "fig03_concordant4_rho")


def fig04_cosmx(df: pd.DataFrame) -> None:
    r50 = df[df["contrast"] == "cytotoxic neighbor ratio 50 um"].iloc[0]
    r100 = df[df["contrast"] == "cytotoxic neighbor ratio 100 um"].iloc[0]
    muz = df[df["contrast"].str.contains("GZMB", regex=False)].iloc[0]

    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(
        fig,
        "Paper funnel  ·  slide 4  ·  CosMx exclusion",
        "CLDN4-high tumor niches have fewer cytotoxic neighbors — not muzzled ones",
    )

    card(
        fig,
        0.06,
        0.42,
        0.27,
        0.36,
        f"{r50['display_effect']}",
        f"cytotoxic neighbor ratio\nat 50 um\n8/8 · 5/5 · sign P={fmt_p(r50['p'])}",
        CORAL,
    )
    card(
        fig,
        0.365,
        0.42,
        0.27,
        0.36,
        f"{r100['display_effect']}",
        f"cytotoxic neighbor ratio\nat 100 um\n8/8 · 5/5 · sign P={fmt_p(r100['p'])}",
        CORAL,
    )
    card(
        fig,
        0.67,
        0.42,
        0.27,
        0.36,
        muz["display_effect"],
        "nearby effector\nGZMB/PRF1/NKG7/IFNG\nhi/lo  ·  0/8 decline",
        MINT,
    )

    fig.text(
        0.03,
        0.28,
        "Write: exclusion, not muzzling.",
        color=GOLD,
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.03,
        0.14,
        "He et al. 2022 CosMx NSCLC (figshare 25976224; 8 sections / 5 patients).\n"
        "Locked ratios reaffirmed in PR 698 (sensitivity does not replace 0.36 / 0.52). "
        "Effector hi/lo range is the public handoff lock.",
        color=MUTED,
        fontsize=12,
        linespacing=1.45,
    )
    footer(fig, "Sources PR 698 + PUBLIC_HANDOFF CosMx lock  ·  Visium same-spot correlation is not spatial exclusion.")
    save(fig, "fig04_cosmx_exclusion")


def fig05_liana(df: pd.DataFrame) -> None:
    rows = df[df["figure_id"] == "liana"].copy()
    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.34, right=0.88, top=0.72, bottom=0.16)
    panel_ax(ax)
    header(
        fig,
        "Paper funnel  ·  slide 5  ·  LIANA barrier",
        "CLDN4-high malignant → barrier/exclusion family is higher on both receivers",
    )

    labels = []
    vals = []
    ps = []
    for _, r in rows.iterrows():
        short = r["contrast"].replace(" barrier_exclusion ", " · ")
        labels.append(short)
        vals.append(r["effect"])
        ps.append(r["p"])
    y = np.arange(len(labels))
    ax.barh(y, vals, color=AMBER, height=0.55, edgecolor=INK, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK, fontsize=11)
    for i, (v, p) in enumerate(zip(vals, ps)):
        ax.text(v + 0.002, i, f"Δ={v:+.3f}   p={fmt_p(p)}", color=INK, va="center", fontsize=11)
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(0, 0.11)
    ax.set_xlabel("Family mean Δ (CLDN4-high − low)  ·  positive = barrier higher from CLDN4-high", color=MUTED)
    ax.text(
        0.0,
        -0.9,
        "Consensus barrier hits include NECTIN2–CD96/TIGIT, CDH1–integrin, F11R–LFA1, LGALS9–CD44/CD45, HLA-E–NKG2A.\n"
        "CXCL9/10/11–CXCR3 are not barrier hits. Expression LR contrast — not spatial exclusion.",
        color=MUTED,
        fontsize=11,
        transform=ax.get_xaxis_transform(),
        va="top",
    )
    footer(fig, "Source PR 579  ·  methods/liana_consensus_concordant4_cldn4/results/tables/family_tests.tsv  ·  n=64 units.")
    save(fig, "fig05_liana_barrier")


def fig06_cptac_priority(df: pd.DataFrame) -> None:
    """Key slide: prioritize CLDN4 within TJ using CPTAC LSCC protein."""
    rows = df[(df["figure_id"] == "cptac") & (df["metric"] == "spearman_rho")].copy()
    part = df[(df["figure_id"] == "cptac") & (df["metric"] == "partial_spearman")].iloc[0]

    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.36, right=0.78, top=0.72, bottom=0.16)
    panel_ax(ax)
    header(
        fig,
        "Paper funnel  ·  slide 6  ·  prioritize CLDN4 within TJ",
        "CPTAC LSCC protein: CLDN4 is immune-inverse; TJ-15 is not",
    )

    order = [
        "CLDN4 protein vs ImmuneScore",
        "CLDN4 protein vs GEP18 RNA",
        "CLDN4 protein vs CD8A RNA",
        "TJ-15 protein vs ImmuneScore",
        "TJ-15 protein vs GEP18 RNA",
        "TJ-15 protein vs CD8A RNA",
    ]
    rows["_ord"] = rows["contrast"].map({k: i for i, k in enumerate(order)})
    rows = rows.sort_values("_ord", ascending=False).reset_index(drop=True)
    y = np.arange(len(rows))
    for i, r in rows.iterrows():
        is_cldn = "CLDN4" in r["contrast"]
        color = CLDN if is_cldn else SLATE
        ax.plot([r["effect_lo"], r["effect_hi"]], [i, i], color=color, lw=3.0, solid_capstyle="round", zorder=2)
        ax.scatter([r["effect"]], [i], s=75, color=color, edgecolor=INK, zorder=3)
        ax.text(
            0.18,
            i,
            f"n={int(r['n'])}   ρ={fmt_rho(r['effect'])}   p={fmt_p(r['p'])}",
            color=INK if is_cldn else MUTED,
            va="center",
            fontsize=11.5,
            fontweight="bold" if is_cldn else "normal",
        )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--", zorder=1)
    ax.set_xlim(-0.75, 0.35)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(rows["contrast"], color=INK, fontsize=10.5)
    ax.set_xlabel("Spearman ρ  (published bootstrap 95% CI)", color=MUTED)
    ax.text(
        0.0,
        -0.85,
        f"CLDN4 protein partial | WES vs ImmuneScore = {fmt_rho(part['effect'])} (n=77, p={fmt_p(part['p'])}). "
        "Treatment-naive surgical LSCC. No ICI labels.\n"
        "This is the public protein argument to prioritize CLDN4 inside the TJ story — not the private KD co-culture.",
        color=MUTED,
        fontsize=11,
        transform=ax.get_xaxis_transform(),
        va="top",
    )
    footer(fig, "Source PR 289  ·  methods/cptac_lusc_cldn4_protein/tables/associations.tsv  ·  CLDN4 vs TJ-15 contrast.")
    save(fig, "fig06_cptac_cldn4_vs_tj15")


def fig07_sig693(df: pd.DataFrame) -> None:
    meta = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("vs CD8A")) & (df["metric"] == "spearman_meta_rho")].iloc[0]
    imm = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("ImmuneScore")) & (df["metric"] == "spearman_meta_rho")].iloc[0]
    base = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("baseline"))].iloc[0]
    cohorts = df[(df["figure_id"] == "sig693") & (df["metric"] == "spearman_rho")].copy()

    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.30, right=0.72, top=0.72, bottom=0.14)
    panel_ax(ax)
    header(
        fig,
        "Paper funnel  ·  slide 7  ·  signature #693",
        "CLDN4-high 221-gene signature: max |ρ| vs CD8A = −0.533 (I²=0%)",
    )

    order = ["OncoSG", "GSE273377 discovery", "GSE273377 validation", "GSE282774", "GSE233774"]
    cohorts["_ord"] = cohorts["accession"].map({k: i for i, k in enumerate(order)})
    cohorts = cohorts.sort_values("_ord", ascending=False).reset_index(drop=True)
    y = np.arange(len(cohorts))
    for i, r in cohorts.iterrows():
        ax.scatter([r["effect"]], [i], s=80, color=CLDN, edgecolor=INK, zorder=3)
        ax.text(
            0.02,
            i,
            f"n={int(r['n'])}   ρ={fmt_rho(r['effect'])}   p={fmt_p(r['p'])}",
            color=INK,
            va="center",
            fontsize=11.5,
        )
    ax.axhspan(-1.2, -0.45, color="#1A2438", zorder=0)
    ax.scatter([meta["effect"]], [-0.82], s=150, marker="D", color=GOLD, edgecolor=INK, zorder=4)
    ax.plot([meta["effect_lo"], meta["effect_hi"]], [-0.82, -0.82], color=GOLD, lw=3.0, zorder=3)
    ax.text(
        0.02,
        -0.82,
        f"meta n={int(meta['n'])}   ρ={fmt_rho(meta['effect'])}   I²=0%",
        color=GOLD,
        va="center",
        fontsize=12,
        fontweight="bold",
    )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(-0.85, 0.25)
    ax.set_ylim(-1.55, len(cohorts) - 0.3)
    ax.set_yticks(list(y) + [-0.82])
    ax.set_yticklabels(list(cohorts["accession"]) + ["META CD8A"], color=INK, fontsize=11)
    ax.set_xlabel("Spearman ρ  (winning spec: ssGSEA α=0.75, first 163 genes, unadjusted)", color=MUTED)
    fig.text(
        0.03,
        0.065,
        f"ImmuneScore at the same spec: meta ρ={fmt_rho(imm['effect'])}. "
        f"Locked 221-gene z-mean baseline: ρ={fmt_rho(base['effect'])} (I²=71%). "
        "Meta p belongs to a maximized |ρ|.",
        color=MUTED,
        fontsize=10.5,
    )
    footer(fig, "Source PR 693  ·  methods/cldn4_sig_max_rho/  ·  gene order fixed from PR 590; size is the only knob.")
    save(fig, "fig07_signature_693")


def fig08_master_forest(df: pd.DataFrame) -> None:
    """One-slide corroboration forest across layers."""
    rows = []

    def add(label, effect, lo, hi, n, p, color):
        rows.append(dict(label=label, effect=effect, lo=lo, hi=hi, n=n, p=p, color=color))

    c4 = df[(df["figure_id"] == "concordant4") & (df["metric"] == "spearman_meta_rho")].iloc[0]
    add("Concordant-4 CLDN4–T/NK", c4["effect"], c4["effect_lo"], c4["effect_hi"], int(c4["n"]), c4["p"], CORAL)

    s = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("vs CD8A")) & (df["metric"] == "spearman_meta_rho")].iloc[0]
    add("Sig #693 vs CD8A", s["effect"], s["effect_lo"], s["effect_hi"], int(s["n"]), s["p"], CLDN)

    si = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("ImmuneScore")) & (df["metric"] == "spearman_meta_rho")].iloc[0]
    add("Sig #693 vs ImmuneScore", si["effect"], si["effect_lo"], si["effect_hi"], int(si["n"]), si["p"], CLDN)

    for contrast, color in [
        ("CLDN4 protein vs ImmuneScore", CLDN),
        ("CLDN4 protein vs GEP18 RNA", CLDN),
        ("CLDN4 protein vs CD8A RNA", CLDN),
        ("TJ-15 protein vs ImmuneScore", SLATE),
    ]:
        r = df[(df["figure_id"] == "cptac") & (df["contrast"] == contrast)].iloc[0]
        add(contrast.replace(" protein", ""), r["effect"], r["effect_lo"], r["effect_hi"], int(r["n"]), r["p"], color)

    fig, ax = plt.subplots(figsize=(W, H))
    apply_slide(fig)
    fig.subplots_adjust(left=0.34, right=0.74, top=0.76, bottom=0.12)
    panel_ax(ax)
    header(
        fig,
        "Paper funnel  ·  slide 8  ·  master corroboration",
        "CLDN4 immune-inverse across patient / protein / signature layers",
    )

    y = np.arange(len(rows))
    for i, r in enumerate(reversed(rows)):
        ax.plot([r["lo"], r["hi"]], [i, i], color=r["color"], lw=3.0, solid_capstyle="round", zorder=2)
        ax.scatter([r["effect"]], [i], s=70, color=r["color"], edgecolor=INK, zorder=3)
        ax.text(
            0.22,
            i,
            f"n={r['n']}   ρ={fmt_rho(r['effect'])}   p={fmt_p(r['p'])}",
            color=INK,
            va="center",
            fontsize=11,
        )
    ax.axvline(0, color=GOLD, lw=1.2, ls="--")
    ax.set_xlim(-0.85, 0.35)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in reversed(rows)], color=INK, fontsize=10.8)
    ax.set_xlabel("Spearman ρ  (published / locked CIs only)", color=MUTED)
    footer(
        fig,
        "Locked from PRs 539 / 289 / 693. CosMx ratios and LIANA barrier Δ are on slides 4–5 (not Spearman). "
        "TJ-15 null is the within-TJ priority argument.",
    )
    save(fig, "fig08_master_forest")


def fig09_multipanel(df: pd.DataFrame) -> None:
    """Compact multipanel for one PPT appendix slide."""
    fig = plt.figure(figsize=(W, H))
    apply_slide(fig)
    header(fig, "Paper funnel  ·  slide 9  ·  multipanel", "Five locked corroboration limbs — CLDN4 first inside TJ")

    # A entry
    ax1 = fig.add_axes([0.05, 0.48, 0.28, 0.34])
    panel_ax(ax1, "A  GSE137244 entry")
    entry = df[df["figure_id"] == "entry_deg"].copy()
    keep = entry[entry["contrast"].isin(["Tacstd2 KL vs KP", "Cldn4 KL vs KP"])]
    ax1.barh(["Tacstd2", "Cldn4"], keep["effect"].values[::-1], color=[GOLD, CLDN], height=0.5)
    ax1.set_xlim(0, 7)
    ax1.tick_params(colors=INK, labelsize=9)
    ax1.set_xlabel("Δ log2", color=MUTED, fontsize=9)

    # B concordant4
    ax2 = fig.add_axes([0.38, 0.48, 0.28, 0.34])
    panel_ax(ax2, "B  Concordant-4")
    meta = df[(df["figure_id"] == "concordant4") & (df["metric"] == "spearman_meta_rho")].iloc[0]
    ax2.barh(["meta ρ"], [meta["effect"]], color=CORAL, height=0.45, xerr=[[meta["effect"] - meta["effect_lo"]], [meta["effect_hi"] - meta["effect"]]], error_kw={"ecolor": INK, "capsize": 4})
    ax2.axvline(0, color=GOLD, ls="--", lw=1)
    ax2.set_xlim(-0.85, 0.15)
    ax2.tick_params(colors=INK, labelsize=9)
    ax2.text(0.02, 0, f"{fmt_rho(meta['effect'])}\nn=65", color=INK, va="center", fontsize=11)

    # C CosMx
    ax3 = fig.add_axes([0.71, 0.48, 0.26, 0.34])
    panel_ax(ax3, "C  CosMx ratios")
    r50 = df[df["contrast"] == "cytotoxic neighbor ratio 50 um"].iloc[0]
    r100 = df[df["contrast"] == "cytotoxic neighbor ratio 100 um"].iloc[0]
    ax3.bar(["50 um", "100 um"], [r50["effect"], r100["effect"]], color=CYAN, edgecolor=INK)
    ax3.axhline(1.0, color=GOLD, ls="--", lw=1)
    ax3.set_ylim(0, 1.15)
    ax3.tick_params(colors=INK, labelsize=9)
    ax3.set_ylabel("hi/lo ratio", color=MUTED, fontsize=9)

    # D CPTAC priority
    ax4 = fig.add_axes([0.05, 0.08, 0.42, 0.32])
    panel_ax(ax4, "D  CPTAC LSCC: CLDN4 vs TJ-15")
    cldn = df[(df["figure_id"] == "cptac") & (df["contrast"] == "CLDN4 protein vs ImmuneScore")].iloc[0]
    tj = df[(df["figure_id"] == "cptac") & (df["contrast"] == "TJ-15 protein vs ImmuneScore")].iloc[0]
    ax4.barh(["CLDN4", "TJ-15"], [cldn["effect"], tj["effect"]], color=[CLDN, SLATE], height=0.45)
    ax4.axvline(0, color=GOLD, ls="--", lw=1)
    ax4.set_xlim(-0.6, 0.2)
    ax4.tick_params(colors=INK, labelsize=9)
    ax4.set_xlabel("ρ vs ImmuneScore", color=MUTED, fontsize=9)

    # E signature
    ax5 = fig.add_axes([0.55, 0.08, 0.42, 0.32])
    panel_ax(ax5, "E  Signature #693 vs CD8A")
    s = df[(df["figure_id"] == "sig693") & (df["contrast"].str.contains("vs CD8A")) & (df["metric"] == "spearman_meta_rho")].iloc[0]
    ax5.barh(["meta ρ"], [s["effect"]], color=CLDN, height=0.45, xerr=[[s["effect"] - s["effect_lo"]], [s["effect_hi"] - s["effect"]]], error_kw={"ecolor": INK, "capsize": 4})
    ax5.axvline(0, color=GOLD, ls="--", lw=1)
    ax5.set_xlim(-0.7, 0.1)
    ax5.tick_params(colors=INK, labelsize=9)
    ax5.text(0.02, 0, f"{fmt_rho(s['effect'])}  I²=0%\nn sum=420", color=INK, va="center", fontsize=11)

    footer(fig, "Locked redraw  ·  PRs 685 / 539 / 698 / 289 / 693  ·  LIANA barrier on slide 5.")
    save(fig, "fig09_multipanel")


def main() -> None:
    df = load_locked()
    fig01_funnel(df)
    fig02_entry_deg(df)
    fig03_concordant4(df)
    fig04_cosmx(df)
    fig05_liana(df)
    fig06_cptac_priority(df)
    fig07_sig693(df)
    fig08_master_forest(df)
    fig09_multipanel(df)
    print(f"Wrote figures to {FIGDIR}")


if __name__ == "__main__":
    main()
