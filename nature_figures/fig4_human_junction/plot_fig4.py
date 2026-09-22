#!/usr/bin/env python3
"""Fig. 4 — public concordant-4 TACSTD2-high malignant DEG to apical junction.

Reads only tables in source_data/. Does not re-fit ORA, GSEA, or OLS.
The expression heatmap is a display of log2(TMM-CPM+1) values copied from
the PR #741 malignant pseudobulk (TACSTD2 matched to units_tacstd2.tsv).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source_data"
OUT = ROOT

for _ttf in (
    "/usr/share/fonts/truetype/croscore/Arimo-Regular.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Bold.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Italic.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-BoldItalic.ttf",
):
    if Path(_ttf).exists():
        font_manager.fontManager.addfont(_ttf)

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

# Published focal junction genes (FINDING core panel plus CDH1). Not an expanded search.
FOCAL_JUNCTION = ["CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "TJP1", "CDH1"]

TERM_LABEL = {
    "HALLMARK_APICAL_JUNCTION": "Hallmark apical junction",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Hallmark EMT",
    "GOBP_KERATINIZATION": "GO keratinization",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Hallmark inflammatory response",
    "HALLMARK_KRAS_SIGNALING_UP": "Hallmark KRAS signaling up",
    "HALLMARK_P53_PATHWAY": "Hallmark p53 pathway",
    "HALLMARK_GLYCOLYSIS": "Hallmark glycolysis",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "Hallmark TNFα via NF-κB",
    "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER": "GO skin barrier",
    "KRT_EPITHELIAL": "Epithelial keratins",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION": "GO tight junction organization",
    "HALLMARK_HYPOXIA": "Hallmark hypoxia",
    "KEGG_TIGHT_JUNCTION": "KEGG tight junction",
}

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arimo", "Arial", "Liberation Sans", "DejaVu Sans"],
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


def fmt_p(p: float) -> str:
    p = float(p)
    if p >= 0.001:
        return f"P = {p:.3f}" if p < 0.1 else f"P = {p:.2f}"
    exp = int(np.floor(np.log10(p)))
    mant = p / (10**exp)
    return f"P = {mant:.1f} × 10{unicode_sup(exp)}"


def fmt_fdr(p: float) -> str:
    p = float(p)
    if p >= 0.01:
        return f"{p:.2f}"
    exp = int(np.floor(np.log10(p)))
    mant = p / (10**exp)
    return f"{mant:.2f} × 10{unicode_sup(exp)}"


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


def add_panel_letter(fig, ax, letter: str, xshift: float = -0.02, yshift: float = 0.012) -> None:
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


def load_tables() -> dict:
    de = pd.read_csv(SRC / "de_q4q1_stacked.tsv.gz", sep="\t")
    ora = pd.read_csv(SRC / "ora_up.tsv", sep="\t")
    barrier = pd.read_csv(SRC / "ora_barrier_subset.tsv", sep="\t")
    focal = pd.read_csv(SRC / "de_focal_genes.tsv", sep="\t")
    n_honest = pd.read_csv(SRC / "n_honest.tsv", sep="\t")
    gsea = pd.read_csv(SRC / "gsea_prerank_q4q1.tsv", sep="\t")
    units = pd.read_csv(SRC / "units_tacstd2.tsv", sep="\t")
    expr = pd.read_csv(SRC / "tj_gene_log2_tmm_cpm_q1q4.tsv", sep="\t")
    return {
        "de": de,
        "ora": ora,
        "barrier": barrier,
        "focal": focal,
        "n_honest": n_honest,
        "gsea": gsea,
        "units": units,
        "expr": expr,
    }


def assert_published(tables: dict) -> dict:
    """Stop if a plotted claim would not match the copied PR #741 tables."""
    de = tables["de"]
    ora = tables["ora"].reset_index(drop=True)
    barrier = tables["barrier"].set_index("term")
    n_honest = tables["n_honest"]
    gsea = tables["gsea"]
    units = tables["units"]
    expr = tables["expr"]
    focal = tables["focal"].set_index("gene")

    if list(ora["term"].head(1)) != ["HALLMARK_APICAL_JUNCTION"]:
        raise SystemExit("ORA rank 1 is not HALLMARK_APICAL_JUNCTION; refusing to plot.")
    apical = ora.iloc[0]
    if int(barrier.loc["HALLMARK_APICAL_JUNCTION", "rank_among_all_ora"]) != 1:
        raise SystemExit("Barrier table does not rank apical junction as 1.")
    if abs(float(apical["enrichment"]) - 5.889814884490932) > 1e-9:
        raise SystemExit("Apical-junction enrichment does not match ora_up.tsv.")
    if int(barrier.loc["KEGG_TIGHT_JUNCTION", "rank_among_all_ora"]) != 31:
        raise SystemExit("KEGG tight junction rank is not 31; refusing to relabel it.")
    if int(barrier.loc["GOBP_TIGHT_JUNCTION_ORGANIZATION", "rank_among_all_ora"]) != 11:
        raise SystemExit("GO tight junction organization rank is not 11.")

    tac = de.loc[de["gene"] == "TACSTD2"].iloc[0]
    if abs(float(tac["logFC"]) - 3.9040230854753495) > 1e-9:
        raise SystemExit("TACSTD2 logFC does not match the DEG table.")
    n_fdr = int((de["fdr"] < 0.05).sum())
    if n_fdr != 1 or de.loc[de["fdr"] < 0.05, "gene"].iloc[0] != "TACSTD2":
        raise SystemExit("FDR < 0.05 set is not TACSTD2 alone.")
    if len(de) != 24083:
        raise SystemExit(f"Unexpected DEG gene count: {len(de)}")
    if int(apical["n_query"]) != 274 or int(apical["n_bg"]) != 24083:
        raise SystemExit("ORA query or background n does not match.")

    n_q1 = int(n_honest["n_q1"].sum())
    n_q4 = int(n_honest["n_q4"].sum())
    if (n_q1, n_q4) != (19, 15):
        raise SystemExit(f"Quartile n is {n_q4} vs {n_q1}, expected 15 vs 19.")
    if int(n_honest["n_vector"].sum()) != 64:
        raise SystemExit("Expression n is not 64.")

    pos = gsea.loc[gsea["nes"] > 0].sort_values("nes", ascending=False).reset_index(drop=True)
    pos.index = pos.index + 1
    aj_rank = int(pos.index[pos["term"] == "HALLMARK_APICAL_JUNCTION"][0])
    kegg_rank = int(pos.index[pos["term"] == "KEGG_TIGHT_JUNCTION"][0])
    if aj_rank != 8 or kegg_rank != 22:
        raise SystemExit(f"GSEA ranks drifted: apical {aj_rank}, KEGG {kegg_rank}.")

    for gene in FOCAL_JUNCTION:
        if gene not in focal.index:
            raise SystemExit(f"{gene} missing from de_focal_genes.tsv.")

    # Expression matrix must be the published TACSTD2 values, not a new quantification.
    tac_expr = expr.loc[expr["gene"] == "TACSTD2", ["patient", "log2_tmm_cpm", "tacstd2_published"]]
    if tac_expr.empty:
        raise SystemExit("Heatmap source is missing TACSTD2 rows.")
    delta = (tac_expr["log2_tmm_cpm"] - tac_expr["tacstd2_published"]).abs().max()
    if delta > 1e-6:
        raise SystemExit(f"Heatmap TACSTD2 disagrees with units_tacstd2.tsv (max {delta}).")
    q_patients = set(units.loc[units["quartile"].isin(["Q1", "Q4"]), "patient"])
    if set(tac_expr["patient"]) != q_patients:
        raise SystemExit("Heatmap patients are not the published Q1/Q4 units.")

    apical_genes = set(str(apical["overlap_genes"]).split(","))
    if "CLDN4" in apical_genes or "TACSTD2" in apical_genes:
        raise SystemExit("CLDN4 or TACSTD2 is inside the apical-junction ORA overlap; check the table.")

    return {
        "n_q1": n_q1,
        "n_q4": n_q4,
        "n_expr": 64,
        "n_genes": 24083,
        "n_up": 274,
        "apical": apical,
        "apical_genes": apical_genes,
        "kegg": barrier.loc["KEGG_TIGHT_JUNCTION"],
        "go_tj": barrier.loc["GOBP_TIGHT_JUNCTION_ORGANIZATION"],
        "gsea_apical": pos.loc[aj_rank],
        "gsea_kegg": pos.loc[kegg_rank],
        "tac": tac,
        "aj_rank": aj_rank,
        "kegg_gsea_rank": kegg_rank,
    }


def gene_blocks(tables: dict, facts: dict) -> list[dict]:
    de = tables["de"].set_index("gene")
    barrier = tables["barrier"].set_index("term")
    apical = list(facts["apical_genes"])
    go_genes = str(barrier.loc["GOBP_TIGHT_JUNCTION_ORGANIZATION", "overlap_genes"]).split(",")
    kegg_genes = str(barrier.loc["KEGG_TIGHT_JUNCTION", "overlap_genes"]).split(",")
    tj_ora = sorted(set(go_genes) | set(kegg_genes))
    focal_only = [g for g in FOCAL_JUNCTION if g not in set(tj_ora) and g not in facts["apical_genes"]]

    def order(genes: list[str]) -> list[str]:
        return sorted(genes, key=lambda g: float(de.loc[g, "logFC"]), reverse=True)

    blocks = [
        {
            "key": "apical",
            "header": "Hallmark apical junction overlap · public ORA rank 1",
            "genes": order(apical),
            "color": PALETTE["blue_main"],
        },
        {
            "key": "tj_ora",
            "header": "Tight-junction ORA overlap (GO organization and KEGG)",
            "genes": order(tj_ora),
            "color": PALETTE["teal"],
        },
        {
            "key": "focal",
            "header": "Focal junction genes outside this ORA query",
            "genes": order(focal_only),
            "color": PALETTE["neutral_dark"],
        },
    ]
    for block in blocks:
        missing = [g for g in block["genes"] if g not in de.index]
        if missing:
            raise SystemExit(f"Genes missing from DEG table: {missing}")
    return blocks


def cohort_quartile_columns(n_honest: pd.DataFrame) -> list[dict]:
    cols = []
    for cohort in COHORT_ORDER:
        row = n_honest.loc[n_honest["cohort"] == cohort].iloc[0]
        for quartile, n_key in (("Q1", "n_q1"), ("Q4", "n_q4")):
            cols.append(
                {
                    "cohort": cohort,
                    "quartile": quartile,
                    "n": int(row[n_key]),
                }
            )
    return cols


def median_matrix(expr: pd.DataFrame, genes: list[str], columns: list[dict]) -> np.ndarray:
    mat = np.full((len(genes), len(columns)), np.nan)
    for i, gene in enumerate(genes):
        sub = expr.loc[expr["gene"] == gene]
        for j, col in enumerate(columns):
            vals = sub.loc[
                (sub["cohort"] == col["cohort"]) & (sub["quartile"] == col["quartile"]),
                "log2_tmm_cpm",
            ]
            if len(vals) != col["n"]:
                raise SystemExit(
                    f"{gene} {col['cohort']} {col['quartile']}: n={len(vals)} expected {col['n']}"
                )
            mat[i, j] = float(vals.median())
    return mat


def row_center(mat: np.ndarray) -> np.ndarray:
    out = np.full_like(mat, np.nan)
    for i in range(mat.shape[0]):
        row = mat[i]
        if np.isnan(row).all():
            continue
        out[i] = row - np.nanmean(row)
    return out


def draw_volcano(ax, de: pd.DataFrame, facts: dict, blocks: list[dict]) -> None:
    de = de.copy()
    de["nlp"] = -np.log10(de["p"].clip(lower=1e-300))
    apical = facts["apical_genes"]
    tj = set(blocks[1]["genes"]) | set(blocks[2]["genes"])
    is_tac = de["gene"] == "TACSTD2"
    is_ap = de["gene"].isin(apical)
    is_tj = de["gene"].isin(tj)
    bg = ~(is_tac | is_ap | is_tj)
    ax.scatter(
        de.loc[bg, "logFC"],
        de.loc[bg, "nlp"],
        s=3.2,
        c=PALETTE["neutral_light"],
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        de.loc[is_tj, "logFC"],
        de.loc[is_tj, "nlp"],
        s=16,
        c=PALETTE["teal"],
        linewidths=0.3,
        edgecolors="white",
        zorder=3,
    )
    ax.scatter(
        de.loc[is_ap, "logFC"],
        de.loc[is_ap, "nlp"],
        s=18,
        c=PALETTE["blue_main"],
        linewidths=0.3,
        edgecolors="white",
        zorder=4,
    )
    ax.scatter(
        de.loc[is_tac, "logFC"],
        de.loc[is_tac, "nlp"],
        s=28,
        c=PALETTE["red_strong"],
        linewidths=0.35,
        edgecolors="white",
        zorder=5,
    )
    ax.axvline(0.25, color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), zorder=2)
    ax.axvline(-0.25, color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), zorder=2)
    ax.axhline(-np.log10(0.01), color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), zorder=2)

    # Hand-tuned offsets (points) so labels sit in open space.
    offsets = {
        "TACSTD2": (-28, 2),
        "ICAM1": (-20, 5),
        "TGFBI": (5, 3),
        "MPZL2": (5, -7),
        "CDH3": (6, 2),
        "NECTIN4": (6, -8),
        "CLDN4": (5, -2),
        "OCLN": (-18, -2),
        "TJP1": (5, -7),
    }
    colors = {
        "TACSTD2": PALETTE["red_strong"],
        "CLDN4": PALETTE["teal"],
        "OCLN": PALETTE["teal"],
        "TJP1": PALETTE["teal"],
    }
    by = de.set_index("gene")
    for gene, (dx, dy) in offsets.items():
        row = by.loc[gene]
        ax.annotate(
            gene,
            xy=(row["logFC"], row["nlp"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=6,
            color=colors.get(gene, PALETTE["blue_main"]),
            fontstyle="italic",
            ha="left" if dx >= 0 else "right",
            va="center",
            arrowprops={
                "arrowstyle": "-",
                "color": PALETTE["neutral_mid"],
                "lw": 0.4,
            },
            zorder=6,
        )

    ax.set_xlim(-2.7, 4.7)
    ax.set_ylim(-0.15, 7.05)
    ax.set_xlabel("log2 fold change (TACSTD2 Q4 − Q1)")
    ax.set_ylabel("−log10 P")
    ax.set_title("Malignant pseudobulk DEG", loc="left", pad=6, fontweight="regular", color=PALETTE["neutral_black"])
    style_ax(ax)
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["red_strong"], markeredgecolor="white", markersize=5, label="TACSTD2"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["blue_main"], markeredgecolor="white", markersize=4.5, label="Apical-junction overlap"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["teal"], markeredgecolor="white", markersize=4.5, label="Tight-junction genes"),
        Line2D([0], [0], color=PALETTE["neutral_mid"], lw=0.6, ls=(0, (2, 1.4)), label="ORA query"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2, handlelength=1.4, labelspacing=0.25)
    ax.text(
        0.98,
        0.98,
        f"15 Q4 vs 19 Q1\nFDR < 0.05: TACSTD2 only\n{fmt_p(float(facts['tac']['p']))}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        linespacing=1.25,
    )


def draw_ora(ax, ora: pd.DataFrame) -> None:
    top = ora.head(12).copy().reset_index(drop=True)
    top["rank"] = np.arange(1, len(top) + 1)
    top["label"] = top["term"].map(lambda t: TERM_LABEL.get(t, t))
    top.loc[top["term"] == "HALLMARK_APICAL_JUNCTION", "label"] = "Hallmark apical junction (rank 1)"
    top["nlp"] = -np.log10(top["fdr"].clip(lower=1e-300))
    # Rank 1 at the top.
    y = (len(top) - top["rank"]).to_numpy()
    cmap = LinearSegmentedColormap.from_list("fdr_blue", ["#E6EEF6", PALETTE["blue_main"]])
    nlp = top["nlp"].to_numpy()
    norm = plt.Normalize(vmin=float(nlp.min()) - 0.15, vmax=float(nlp.max()))
    colors = cmap(norm(nlp))
    sizes = (top["n_overlap"].to_numpy().astype(float) ** 1.55) * 4.2
    edge = []
    lw = []
    for term in top["term"]:
        if term == "HALLMARK_APICAL_JUNCTION":
            edge.append(PALETTE["red_strong"])
            lw.append(1.15)
        elif term == "GOBP_TIGHT_JUNCTION_ORGANIZATION":
            edge.append(PALETTE["teal"])
            lw.append(1.05)
        else:
            edge.append("white")
            lw.append(0.35)
    ax.scatter(
        top["enrichment"],
        y,
        s=sizes,
        c=colors,
        edgecolors=edge,
        linewidths=lw,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["label"])
    for tick, term in zip(ax.get_yticklabels(), top["term"]):
        if term == "HALLMARK_APICAL_JUNCTION":
            tick.set_color(PALETTE["blue_main"])
            tick.set_fontweight("bold")
        elif term == "GOBP_TIGHT_JUNCTION_ORGANIZATION":
            tick.set_color(PALETTE["teal"])
            tick.set_fontweight("bold")
        else:
            tick.set_color(PALETTE["neutral_black"])
    ax.set_xlabel("Enrichment (overlap / expected)")
    ax.set_ylim(-0.65, len(top) - 0.35)
    ax.set_xlim(0, 17.2)
    ax.set_title(
        "Public ORA, TACSTD2-high up genes\nP < 0.01, |logFC| > 0.25, n = 274 (TACSTD2 held out)",
        loc="left",
        pad=4,
        color=PALETTE["neutral_black"],
        linespacing=1.25,
    )
    style_ax(ax)
    ax.tick_params(axis="y", length=0, pad=2)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = ax.figure.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.055, pad=0.16, aspect=28)
    cbar.ax.tick_params(labelsize=5.5, length=2, width=0.5, pad=1)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label("−log10 FDR", fontsize=6, color=PALETTE["neutral_dark"])
    # Size key sits in the open upper-right of the panel (enrichment > 8, ranks 1–2).
    ax.scatter(
        [13.2, 14.6, 16.1],
        [10.15, 10.15, 10.15],
        s=(np.array([5, 8, 13]) ** 1.55) * 4.2,
        c="white",
        edgecolors=PALETTE["neutral_dark"],
        linewidths=0.5,
        zorder=2,
    )
    for x, n in ((13.2, "5"), (14.6, "8"), (16.1, "13")):
        ax.text(x, 9.55, n, ha="center", va="top", fontsize=5.5, color=PALETTE["neutral_dark"])
    ax.text(14.65, 10.72, "Overlap, genes", ha="center", va="bottom", fontsize=5.5, color=PALETTE["neutral_dark"])


def draw_heatmap(ax, tables: dict, blocks: list[dict]) -> np.ndarray:
    columns = cohort_quartile_columns(tables["n_honest"])
    # Build display rows: header (nan) then genes.
    row_meta = []
    for block in blocks:
        row_meta.append({"kind": "header", "label": block["header"], "color": block["color"], "gene": None})
        for gene in block["genes"]:
            row_meta.append({"kind": "gene", "label": gene, "color": block["color"], "gene": gene})
    genes = [r["gene"] for r in row_meta if r["kind"] == "gene"]
    med = median_matrix(tables["expr"], genes, columns)
    centered = row_center(med)
    # Map gene rows into the full grid.
    grid = np.full((len(row_meta), len(columns)), np.nan)
    logfc = []
    pvals = []
    de = tables["de"].set_index("gene")
    gi = 0
    for i, row in enumerate(row_meta):
        if row["kind"] == "gene":
            grid[i] = centered[gi]
            logfc.append(float(de.loc[row["gene"], "logFC"]))
            pvals.append(float(de.loc[row["gene"], "p"]))
            gi += 1
        else:
            logfc.append(np.nan)
            pvals.append(np.nan)

    cmap = LinearSegmentedColormap.from_list(
        "delta",
        [PALETTE["blue_main"], "#F7F7F7", PALETTE["red_strong"]],
    )
    cmap.set_bad("#FFFFFF")
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5)
    ax.imshow(grid, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_xticks(range(len(columns)))
    xticklabels = []
    for col in columns:
        xticklabels.append(f"{col['quartile']}\n{col['n']}")
    ax.set_xticklabels(xticklabels, fontsize=6)
    for tick, col in zip(ax.get_xticklabels(), columns):
        tick.set_color(COHORT_COLORS[col["cohort"]])
    # Cohort names centered over each pair.
    for i, cohort in enumerate(COHORT_ORDER):
        ax.text(
            i * 2 + 0.5,
            -0.85,
            cohort,
            ha="center",
            va="bottom",
            fontsize=6,
            color=COHORT_COLORS[cohort],
            clip_on=False,
        )
    ax.set_yticks(range(len(row_meta)))
    ax.set_yticklabels([r["label"] if r["kind"] == "gene" else "" for r in row_meta], fontsize=6)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")
        tick.set_color(PALETTE["neutral_black"])
    for i, row in enumerate(row_meta):
        if row["kind"] == "header":
            ax.text(
                0.02,
                i,
                row["label"],
                ha="left",
                va="center",
                fontsize=6,
                fontweight="bold",
                color=row["color"],
                clip_on=False,
            )
    for boundary in (1.5, 3.5, 5.5):
        ax.axvline(boundary, color="white", lw=1.4, zorder=3)
    # OLS logFC from the published DEG table, not from the medians.
    for i, (lf, pv) in enumerate(zip(logfc, pvals)):
        if not np.isfinite(lf):
            continue
        ax.text(
            len(columns) - 0.5 + 0.42,
            i,
            f"{lf:+.2f}",
            va="center",
            ha="left",
            fontsize=5.5,
            color=PALETTE["neutral_black"],
            fontweight="bold" if pv < 0.01 else "regular",
            clip_on=False,
        )
    ax.text(
        len(columns) - 0.5 + 0.42,
        -0.85,
        "OLS logFC",
        ha="left",
        va="bottom",
        fontsize=6,
        color=PALETTE["neutral_black"],
        clip_on=False,
    )
    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(len(row_meta) - 0.5, -1.35)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "Junction genes, TACSTD2 Q4 versus Q1",
        loc="left",
        pad=2,
        color=PALETTE["neutral_black"],
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = ax.figure.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.035, pad=0.12, aspect=40)
    cbar.ax.tick_params(labelsize=5.5, length=2, width=0.5, pad=1)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label(
        "Row-centered median log2(TMM-CPM+1), scale ±2.5. Bold OLS logFC: nominal P < 0.01.",
        fontsize=6,
        color=PALETTE["neutral_dark"],
    )
    return grid


def write_caption(facts: dict, focal: pd.DataFrame) -> None:
    apical = facts["apical"]
    kegg = facts["kegg"]
    go_tj = facts["go_tj"]
    gsea_aj = facts["gsea_apical"]
    gsea_kegg = facts["gsea_kegg"]
    tac = facts["tac"]
    cldn4 = focal.set_index("gene").loc["CLDN4"]
    text = f"""# Fig. 4 — Human public TACSTD2-high malignant cells and apical junction

**Fig. 4 | Hallmark apical junction is the leading public ORA term in TACSTD2-high malignant cells.**

Concordant-4 malignant pseudobulk only (GSE123902, GSE131907, GSE205335, GSE189357). The split is within-cohort TACSTD2 expression quartile, not CLDN4. Differential expression is cohort-adjusted OLS on log2(TMM-CPM+1). Positive log2 fold change is higher in TACSTD2 Q4. Expression n = 64 units; the Q4 versus Q1 contrast is 15 versus 19. These panels are the public result in PR #741. They are not a private 8KL KEGG tight-junction rank 1.

**a**, Volcano of 24,083 genes. Dashed lines are the ORA query rule (nominal P < 0.01 and |logFC| > 0.25). TACSTD2 logFC = {float(tac['logFC']):+.2f}, {fmt_p(float(tac['p']))}, FDR = {float(tac['fdr']):.3f}. It is the only gene with FDR < 0.05. Blue points are the 13 Hallmark apical-junction genes that overlap the up-gene query. Teal points are tight-junction ORA overlaps and the focal junction genes. CLDN4 (logFC {float(cldn4['logFC']):+.2f}, nominal P = {float(cldn4['p']):.3f}) is labeled and is not in the apical-junction overlap.

**b**, Over-representation of 274 up genes (TACSTD2 held out) on a background of 24,083 genes. Hallmark apical junction is public ORA rank 1 ({int(apical['n_overlap'])}/{int(apical['n_set_in_bg'])} genes, enrichment {float(apical['enrichment']):.2f}, FDR {fmt_fdr(float(apical['fdr']))}). GO tight-junction organization is rank {int(go_tj['rank_among_all_ora'])} (enrichment {float(go_tj['enrichment']):.2f}, FDR {fmt_fdr(float(go_tj['fdr']))}). KEGG tight junction is rank {int(kegg['rank_among_all_ora'])} (enrichment {float(kegg['enrichment']):.2f}, FDR {fmt_fdr(float(kegg['fdr']))}) and is not rank 1. Bubble color is −log10 FDR; size is overlap count. The same contrast, preranked GSEA (OLS t, TACSTD2 removed from the ranking, 1,000 permutations, seed 42), places Hallmark apical junction at NES {float(gsea_aj['nes']):+.2f}, rank {facts['aj_rank']} of positive NES (FDR {fmt_fdr(float(gsea_aj['fdr']))}), and KEGG tight junction at NES {float(gsea_kegg['nes']):+.2f}, rank {facts['kegg_gsea_rank']} (FDR {fmt_fdr(float(gsea_kegg['fdr']))}). GSEA is not drawn as a separate panel.

**c**, Within-cohort median log2(TMM-CPM+1) for the same gene groups, row-centered within each gene for display (color truncated at ±2.5). Column n is the published Q1 or Q4 count. Numbers at the right are the published cohort-adjusted OLS logFC; bold marks nominal P < 0.01. Medians are descriptive. They are not a second differential-expression test. GSE189357 Q4 is two patients.

Source tables are copies from PR #741 (`methods/concordant4_tacstd2_malignant_deg/tables/`). The heatmap matrix is log2(TMM-CPM+1) from that page’s malignant UMI sums; TACSTD2 matches `units_tacstd2.tsv`.
"""
    (OUT / "captions.md").write_text(text)


def write_panel_map(facts: dict) -> None:
    apical = facts["apical"]
    kegg = facts["kegg"]
    go_tj = facts["go_tj"]
    rows = [
        {
            "panel": "a",
            "content": "Volcano, cohort-adjusted OLS logFC vs -log10 P, TACSTD2 Q4 vs Q1",
            "n": "15 Q4 vs 19 Q1; 24083 genes; expression units 64",
            "highlight": "TACSTD2; Hallmark apical-junction ORA overlap; tight-junction genes",
            "key_number": f"TACSTD2 logFC {float(facts['tac']['logFC']):.3f}; FDR<0.05 only TACSTD2",
            "source_file": "source_data/de_q4q1_stacked.tsv.gz",
            "source_pr": "741",
            "label_constraint": "Public concordant-4 malignant pseudobulk. Not private 8KL.",
        },
        {
            "panel": "b",
            "content": "ORA bubble, top 12 up-gene terms; Hallmark apical junction called as rank 1",
            "n": "274 up genes; background 24083; TACSTD2 held out",
            "highlight": "HALLMARK_APICAL_JUNCTION rank 1; GO TJ organization rank 11; KEGG TJ annotated as rank 31",
            "key_number": (
                f"apical enrichment {float(apical['enrichment']):.2f} FDR {float(apical['fdr']):.3e}; "
                f"GO TJ rank {int(go_tj['rank_among_all_ora'])}; "
                f"KEGG TJ rank {int(kegg['rank_among_all_ora'])} FDR {float(kegg['fdr']):.3f}"
            ),
            "source_file": "source_data/ora_up.tsv; source_data/ora_barrier_subset.tsv",
            "source_pr": "741",
            "label_constraint": "Public ORA rank 1 is Hallmark apical junction. Do not call this KEGG TJ #1 or private 8KL KEGG TJ #1.",
        },
        {
            "panel": "c",
            "content": "Heatmap of within-cohort median log2(TMM-CPM+1), row-centered; OLS logFC printed from the DEG table",
            "n": "Q1/Q4 per cohort as in n_honest.tsv (4/3, 6/5, 6/5, 3/2)",
            "highlight": "Apical-junction overlap block; tight-junction ORA overlap; focal junction genes outside that ORA",
            "key_number": "Display only. Inferential logFC is the printed OLS column.",
            "source_file": "source_data/tj_gene_log2_tmm_cpm_q1q4.tsv; source_data/de_q4q1_stacked.tsv.gz; source_data/n_honest.tsv",
            "source_pr": "741",
            "label_constraint": "Medians are descriptive. CLDN4 is not part of the apical-junction ORA overlap.",
        },
        {
            "panel": "caption_only",
            "content": "Preranked GSEA NES ranks on the same Q4 vs Q1 t statistic",
            "n": "1000 permutations; TACSTD2 dropped from the ranking; 36 positive-NES sets",
            "highlight": "Hallmark apical junction NES rank 8; KEGG tight junction NES rank 22",
            "key_number": (
                f"apical NES {float(facts['gsea_apical']['nes']):.3f} rank {facts['aj_rank']}; "
                f"KEGG NES {float(facts['gsea_kegg']['nes']):.3f} rank {facts['kegg_gsea_rank']}"
            ),
            "source_file": "source_data/gsea_prerank_q4q1.tsv",
            "source_pr": "741",
            "label_constraint": "GSEA does not put junction first. Not plotted as its own panel.",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "panel_map.tsv", sep="\t", index=False)


def save(fig) -> None:
    stem = OUT / "fig4_human_junction"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.04, facecolor="white")


def main() -> None:
    apply_style()
    tables = load_tables()
    facts = assert_published(tables)
    blocks = gene_blocks(tables, facts)
    fig = plt.figure(figsize=(mm_to_in(183), mm_to_in(202)))
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.0, 1.46, 0.10],
        width_ratios=[1.0, 1.22],
        left=0.08,
        right=0.90,
        top=0.965,
        bottom=0.02,
        wspace=0.62,
        hspace=0.46,
    )
    ax_vol = fig.add_subplot(gs[0, 0])
    ax_ora = fig.add_subplot(gs[0, 1])
    ax_hm = fig.add_subplot(gs[1, :])
    ax_note = fig.add_subplot(gs[2, :])
    ax_note.axis("off")
    draw_volcano(ax_vol, tables["de"], facts, blocks)
    draw_ora(ax_ora, tables["ora"])
    draw_heatmap(ax_hm, tables, blocks)
    kegg = facts["kegg"]
    ax_note.text(
        0.5,
        0.55,
        "Public ORA rank 1 is Hallmark apical junction, not KEGG tight junction.\n"
        f"KEGG tight junction is rank {int(kegg['rank_among_all_ora'])} "
        f"(enrichment {float(kegg['enrichment']):.2f}, FDR {fmt_fdr(float(kegg['fdr']))}). "
        f"GSEA on this contrast: apical junction NES rank {facts['aj_rank']}; "
        f"KEGG tight junction NES rank {facts['kegg_gsea_rank']}. Not a private 8KL result.",
        ha="center",
        va="center",
        fontsize=6,
        color=PALETTE["neutral_dark"],
        linespacing=1.35,
    )
    fig.canvas.draw()
    add_panel_letter(fig, ax_vol, "a", xshift=-0.035, yshift=0.006)
    add_panel_letter(fig, ax_ora, "b", xshift=-0.05, yshift=0.006)
    add_panel_letter(fig, ax_hm, "c", xshift=-0.02, yshift=0.012)
    save(fig)
    write_caption(facts, tables["focal"])
    write_panel_map(facts)
    plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
