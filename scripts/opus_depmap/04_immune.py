#!/usr/bin/env python3
"""TACSTD2 / CLDN4 versus IFN, MHC-I and tumour-cell immune genes in lung lines.

Three complementary tests, all on open DepMap Public 24Q4 files:

A. Expression vs expression (CCLE RNA-seq, log2(TPM+1)).
   Spearman correlation of TACSTD2 or CLDN4 with each panel gene, plus two
   composite scores (MHC-I antigen-presentation mean z-score; IFN-response
   mean z-score).

B. Dependency vs immune expression.
   Spearman correlation of TACSTD2 or CLDN4 Chronos gene effect with each
   panel gene's expression and with the two composite scores. This is the
   test implied by "dependency vs IFN/MHC-I/immune genes".

C. Dependency vs immune-gene dependency.
   Pearson correlation of TACSTD2 or CLDN4 Chronos profiles with the Chronos
   profiles of the same panel genes (targeted co-dependency, not a genome-wide
   scan).

Every row reports n, effect size, a two-sided p-value and a BH q-value
computed *within that test block*. NSCLC-only and SCLC-adjusted (partial
Spearman controlling for a binary NSCLC indicator) versions of A and B are
included because neuroendocrine vs NSCLC identity is a known confounder of
both epithelial and MHC-I expression.

Panel membership is pre-specified from textbook antigen-presentation / IFN
signalling / tumour-cell ligand lists. T-cell lineage markers are included
only as negative-control genes that should be near-absent in cancer lines.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    FIGURES,
    GENES_OF_INTEREST,
    NOTES,
    RELEASE,
    TABLES,
    bh_fdr,
    ensure_dirs,
    fisher_ci,
    load_expression,
    load_gene_effect,
    load_models,
    lung_cohort,
    partial_spearman,
)

# Pre-specified panels. Symbols must match DepMap gene-effect / expression
# headers after Entrez stripping. Genes missing from a matrix are dropped
# and recorded, not silently invented.
PANELS: dict[str, list[str]] = {
    "MHC-I antigen presentation": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "TAPBP",
        "NLRC5",
        "PSMB8",
        "PSMB9",
        "PSMB10",
        "ERAP1",
        "ERAP2",
        "CALR",
        "CANX",
        "PDIA3",
    ],
    "IFN signalling": [
        "IFNGR1",
        "IFNGR2",
        "IFNAR1",
        "IFNAR2",
        "JAK1",
        "JAK2",
        "TYK2",
        "STAT1",
        "STAT2",
        "IRF1",
        "IRF9",
        "SOCS1",
    ],
    "IFN response / ISG": [
        "ISG15",
        "MX1",
        "OAS1",
        "IFI27",
        "IFI44",
        "IFIT1",
        "IFIT3",
        "CXCL9",
        "CXCL10",
        "IDO1",
        "GBP1",
        "IFI6",
        "RSAD2",
        "IFITM1",
    ],
    "tumour-cell immune ligand": [
        "CD274",
        "PDCD1LG2",
        "CD47",
        "PVR",
        "NECTIN2",
        "LGALS9",
        "TNFRSF14",
        "HLA-E",
        "HLA-G",
    ],
    "epithelial positive control": [
        "EPCAM",
        "CDH1",
        "KRT8",
        "KRT18",
        "KRT19",
        "CLDN3",
        "CLDN7",
    ],
    "lymphoid negative control": ["CD3E", "CD8A", "CD4", "CD19"],
}

COMPOSITE_MEMBERS = {
    "MHC-I score": PANELS["MHC-I antigen presentation"],
    "IFN-response score": PANELS["IFN response / ISG"],
}


def present(symbols: list[str], columns) -> list[str]:
    return [g for g in symbols if g in columns]


def zscore_mean(frame: pd.DataFrame) -> pd.Series:
    """Mean of column-wise z-scores; requires at least 80% non-missing genes."""
    z = (frame - frame.mean(axis=0)) / frame.std(axis=0, ddof=1)
    score = z.mean(axis=1)
    frac = frame.notna().mean(axis=1)
    score = score.where(frac >= 0.8)
    return score


def spearman_row(x: pd.Series, y: pd.Series) -> dict:
    pair = pd.concat([x, y], axis=1).dropna()
    n = int(len(pair))
    if n < 8:
        return {"n": n, "r": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    r, p = stats.spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
    lo, hi = fisher_ci(float(r), n)
    return {"n": n, "r": float(r), "ci_low": lo, "ci_high": hi, "p": float(p)}


def pearson_row(x: pd.Series, y: pd.Series) -> dict:
    pair = pd.concat([x, y], axis=1).dropna()
    n = int(len(pair))
    if n < 8:
        return {"n": n, "r": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    r, p = stats.pearsonr(pair.iloc[:, 0], pair.iloc[:, 1])
    lo, hi = fisher_ci(float(r), n)
    return {"n": n, "r": float(r), "ci_low": lo, "ci_high": hi, "p": float(p)}


def annotate_fdr(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    parts = []
    for _, sub in frame.groupby(group_cols, dropna=False):
        sub = sub.copy()
        sub["q"] = bh_fdr(sub["p"].to_numpy())
        parts.append(sub)
    return pd.concat(parts).sort_index()


def main() -> int:
    ensure_dirs()
    models = load_models()
    lung = lung_cohort(models)
    expr = load_expression()
    gene_effect = load_gene_effect()

    lung_ex = lung.index.intersection(expr.index)
    lung_ge = lung.index.intersection(gene_effect.index)
    lung_both = lung_ex.intersection(lung_ge)
    nsclc_ex = lung.index[lung.LungGroup == "NSCLC"].intersection(expr.index)
    nsclc_both = lung.index[lung.LungGroup == "NSCLC"].intersection(lung_both)
    sclc_ex = lung.index[lung.LungGroup == "SCLC"].intersection(expr.index)

    # Inventory of panel genes actually present
    inventory = []
    for panel, genes in PANELS.items():
        for g in genes:
            inventory.append(
                {
                    "panel": panel,
                    "gene": g,
                    "in_expression": g in expr.columns,
                    "in_gene_effect": g in gene_effect.columns,
                    "lung_expr_n": int(expr.loc[lung_ex, g].notna().sum()) if g in expr.columns else 0,
                    "lung_expr_mean": (
                        float(expr.loc[lung_ex, g].mean()) if g in expr.columns else np.nan
                    ),
                    "lung_expr_sd": (
                        float(expr.loc[lung_ex, g].std(ddof=1)) if g in expr.columns else np.nan
                    ),
                }
            )
    inv = pd.DataFrame(inventory)
    inv.to_csv(TABLES / "immune_panel_inventory.csv", index=False)
    missing_expr = inv.loc[~inv["in_expression"], "gene"].tolist()
    missing_ge = inv.loc[~inv["in_gene_effect"], "gene"].tolist()
    print(f"panel genes missing from expression: {missing_expr or 'none'}")
    print(f"panel genes missing from gene effect: {missing_ge or 'none'}")

    # Composite scores on lung expression
    mhc_genes = present(COMPOSITE_MEMBERS["MHC-I score"], expr.columns)
    ifn_genes = present(COMPOSITE_MEMBERS["IFN-response score"], expr.columns)
    scores = pd.DataFrame(
        {
            "MHC-I score": zscore_mean(expr.loc[lung_ex, mhc_genes]),
            "IFN-response score": zscore_mean(expr.loc[lung_ex, ifn_genes]),
        }
    )
    scores["is_NSCLC"] = (lung.reindex(scores.index)["LungGroup"] == "NSCLC").astype(float)
    scores.to_csv(TABLES / "immune_composite_scores_lung.csv")

    # ------------------------------------------------------------------
    # A. expression vs expression
    # ------------------------------------------------------------------
    expr_rows = []
    partners = []
    for panel, genes in PANELS.items():
        for g in present(genes, expr.columns):
            partners.append((panel, g))
    partners += [("composite", "MHC-I score"), ("composite", "IFN-response score")]

    cohorts = {
        "lung": lung_ex,
        "NSCLC": nsclc_ex,
        "SCLC": sclc_ex,
    }
    for query in GENES_OF_INTEREST:
        for cohort_name, ids in cohorts.items():
            qx = expr.loc[ids, query]
            for panel, partner in partners:
                if partner in scores.columns:
                    y = scores.loc[ids.intersection(scores.index), partner]
                else:
                    y = expr.loc[ids, partner]
                row = spearman_row(qx, y)
                row.update(
                    {
                        "test": "expr_vs_expr",
                        "query": query,
                        "partner": partner,
                        "panel": panel,
                        "cohort": cohort_name,
                    }
                )
                expr_rows.append(row)
            # partial Spearman controlling for NSCLC vs not, lung only
            if cohort_name == "lung":
                for panel, partner in partners:
                    if partner in scores.columns:
                        y = scores.loc[ids.intersection(scores.index), partner]
                    else:
                        y = expr.loc[ids, partner]
                    r, p, n = partial_spearman(qx, y, scores.reindex(ids)["is_NSCLC"])
                    expr_rows.append(
                        {
                            "test": "expr_vs_expr_partial_NSCLC",
                            "query": query,
                            "partner": partner,
                            "panel": panel,
                            "cohort": "lung_partial_NSCLC",
                            "n": n,
                            "r": r,
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "p": p,
                        }
                    )

    expr_frame = pd.DataFrame(expr_rows)
    expr_frame = annotate_fdr(expr_frame, ["test", "query", "cohort"])
    expr_frame = expr_frame[
        ["test", "query", "partner", "panel", "cohort", "n", "r", "ci_low", "ci_high", "p", "q"]
    ]
    expr_frame.to_csv(TABLES / "immune_expr_vs_expr.csv", index=False)

    # ------------------------------------------------------------------
    # B. dependency (gene effect) vs immune expression
    # ------------------------------------------------------------------
    dep_rows = []
    dep_cohorts = {
        "lung": lung_both,
        "NSCLC": nsclc_both,
    }
    for query in GENES_OF_INTEREST:
        for cohort_name, ids in dep_cohorts.items():
            qg = gene_effect.loc[ids, query]
            for panel, partner in partners:
                if partner in scores.columns:
                    y = scores.loc[ids.intersection(scores.index), partner]
                else:
                    y = expr.loc[ids, partner]
                row = spearman_row(qg, y)
                row.update(
                    {
                        "test": "geneEffect_vs_expr",
                        "query": query,
                        "partner": partner,
                        "panel": panel,
                        "cohort": cohort_name,
                    }
                )
                dep_rows.append(row)
            if cohort_name == "lung":
                for panel, partner in partners:
                    if partner in scores.columns:
                        y = scores.loc[ids.intersection(scores.index), partner]
                    else:
                        y = expr.loc[ids, partner]
                    r, p, n = partial_spearman(qg, y, scores.reindex(ids)["is_NSCLC"])
                    dep_rows.append(
                        {
                            "test": "geneEffect_vs_expr_partial_NSCLC",
                            "query": query,
                            "partner": partner,
                            "panel": panel,
                            "cohort": "lung_partial_NSCLC",
                            "n": n,
                            "r": r,
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "p": p,
                        }
                    )
    dep_frame = pd.DataFrame(dep_rows)
    dep_frame = annotate_fdr(dep_frame, ["test", "query", "cohort"])
    dep_frame = dep_frame[
        ["test", "query", "partner", "panel", "cohort", "n", "r", "ci_low", "ci_high", "p", "q"]
    ]
    dep_frame.to_csv(TABLES / "immune_geneEffect_vs_expr.csv", index=False)

    # ------------------------------------------------------------------
    # C. dependency vs immune-gene dependency
    # ------------------------------------------------------------------
    ge_rows = []
    ge_partners = []
    for panel, genes in PANELS.items():
        for g in present(genes, gene_effect.columns):
            ge_partners.append((panel, g))
    for query in GENES_OF_INTEREST:
        for cohort_name, ids in (
            ("lung", lung_ge),
            ("NSCLC", lung.index[lung.LungGroup == "NSCLC"].intersection(gene_effect.index)),
            ("pan-cancer", gene_effect.index),
        ):
            qg = gene_effect.loc[ids, query]
            for panel, partner in ge_partners:
                if partner == query:
                    continue
                row = pearson_row(qg, gene_effect.loc[ids, partner])
                row.update(
                    {
                        "test": "geneEffect_vs_geneEffect",
                        "query": query,
                        "partner": partner,
                        "panel": panel,
                        "cohort": cohort_name,
                    }
                )
                ge_rows.append(row)
    ge_frame = pd.DataFrame(ge_rows)
    ge_frame = annotate_fdr(ge_frame, ["test", "query", "cohort"])
    ge_frame = ge_frame[
        ["test", "query", "partner", "panel", "cohort", "n", "r", "ci_low", "ci_high", "p", "q"]
    ]
    ge_frame.to_csv(TABLES / "immune_geneEffect_vs_geneEffect.csv", index=False)

    # ------------------------------------------------------------------
    # Compact headline table: composites + MHC-I core + IFN core
    # ------------------------------------------------------------------
    headline_partners = [
        "MHC-I score",
        "IFN-response score",
        "B2M",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "NLRC5",
        "TAP1",
        "STAT1",
        "IRF1",
        "IFNGR1",
        "CD274",
        "EPCAM",
        "CDH1",
    ]
    headline = pd.concat(
        [
            expr_frame[expr_frame.partner.isin(headline_partners)],
            dep_frame[dep_frame.partner.isin(headline_partners)],
        ],
        ignore_index=True,
    )
    headline.to_csv(TABLES / "immune_headline.csv", index=False)

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_genes = (
        present(PANELS["MHC-I antigen presentation"], expr.columns)
        + present(PANELS["IFN signalling"], expr.columns)
        + present(PANELS["IFN response / ISG"], expr.columns)
        + present(PANELS["tumour-cell immune ligand"], expr.columns)
        + present(PANELS["epithelial positive control"], expr.columns)
    )
    heat = (
        expr_frame[
            (expr_frame.test == "expr_vs_expr")
            & (expr_frame.cohort == "lung")
            & (expr_frame.partner.isin(plot_genes))
        ]
        .pivot(index="partner", columns="query", values="r")
        .reindex(plot_genes)
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 10), gridspec_kw={"width_ratios": [1.1, 1.1, 1.4]})

    ax = axes[0]
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(heat.columns.tolist())
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index.tolist(), fontsize=7)
    ax.set_title("Spearman r\nexpression vs expression\nlung lines")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[1]
    heat2 = (
        dep_frame[
            (dep_frame.test == "geneEffect_vs_expr")
            & (dep_frame.cohort == "lung")
            & (dep_frame.partner.isin(plot_genes))
        ]
        .pivot(index="partner", columns="query", values="r")
        .reindex(plot_genes)
    )
    im2 = ax.imshow(heat2.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.5, vmax=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(heat2.columns.tolist())
    ax.set_yticks(range(len(heat2.index)))
    ax.set_yticklabels(heat2.index.tolist(), fontsize=7)
    ax.set_title("Spearman r\ngene effect vs expression\nlung lines with CRISPR")
    fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[2]
    for query, colour, marker in (("TACSTD2", "#1a9850", "o"), ("CLDN4", "#762a83", "s")):
        x = expr.loc[lung_ex, query]
        y = scores.loc[lung_ex, "MHC-I score"]
        grp = lung.reindex(lung_ex)["LungGroup"]
        for gname, face in (("NSCLC", colour), ("SCLC", "white"), ("Other lung", "grey")):
            mask = grp == gname
            ax.scatter(
                x[mask],
                y[mask],
                s=22,
                marker=marker,
                facecolors=face if gname != "SCLC" else "white",
                edgecolors=colour,
                alpha=0.75,
                label=f"{query} {gname}",
            )
    ax.set_xlabel("query expression, log2(TPM+1)")
    ax.set_ylabel("MHC-I composite score (mean z)")
    ax.set_title("TACSTD2 / CLDN4 expression vs MHC-I score")
    ax.legend(fontsize=7, loc="best")
    fig.suptitle(f"{RELEASE} · lung lines · TACSTD2 / CLDN4 vs IFN and MHC-I", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_immune.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Note — numbers come from the frames just written, not from memory
    # ------------------------------------------------------------------
    def pick(frame, test, query, partner, cohort):
        hit = frame[
            (frame.test == test)
            & (frame.query == query)
            & (frame.partner == partner)
            & (frame.cohort == cohort)
        ]
        return None if hit.empty else hit.iloc[0]

    lines = [
        "# 04 · TACSTD2 / CLDN4 versus IFN, MHC-I and immune genes",
        "",
        f"Release {RELEASE}. Lung lines with expression: **{len(lung_ex)}** "
        f"(NSCLC {len(nsclc_ex)}, SCLC {len(sclc_ex)}). Lung lines with both "
        f"CRISPR and expression: **{len(lung_both)}** (NSCLC {len(nsclc_both)}).",
        "",
        f"MHC-I composite uses {len(mhc_genes)} genes; IFN-response composite uses "
        f"{len(ifn_genes)} genes. Missing panel genes are listed in "
        "`immune_panel_inventory.csv`.",
        "",
        "## Composite scores (headline)",
        "",
        "| test | query | partner | cohort | n | r | 95% CI | p | q |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for test, frame, partner, cohort in (
        ("expr_vs_expr", expr_frame, "MHC-I score", "lung"),
        ("expr_vs_expr", expr_frame, "IFN-response score", "lung"),
        ("expr_vs_expr", expr_frame, "MHC-I score", "NSCLC"),
        ("expr_vs_expr_partial_NSCLC", expr_frame, "MHC-I score", "lung_partial_NSCLC"),
        ("geneEffect_vs_expr", dep_frame, "MHC-I score", "lung"),
        ("geneEffect_vs_expr", dep_frame, "IFN-response score", "lung"),
        ("geneEffect_vs_expr", dep_frame, "MHC-I score", "NSCLC"),
        ("geneEffect_vs_expr_partial_NSCLC", dep_frame, "MHC-I score", "lung_partial_NSCLC"),
    ):
        for query in GENES_OF_INTEREST:
            r = pick(frame, test, query, partner, cohort)
            if r is None:
                continue
            ci = (
                f"[{r.ci_low:.3f}, {r.ci_high:.3f}]"
                if np.isfinite(r.ci_low)
                else "n/a (partial)"
            )
            lines.append(
                f"| {test} | {query} | {partner} | {cohort} | {int(r.n)} | {r.r:.3f} | "
                f"{ci} | {r.p:.3g} | {r.q:.3g} |"
            )

    lines += [
        "",
        "Full per-gene tables: `immune_expr_vs_expr.csv`, `immune_geneEffect_vs_expr.csv`, "
        "`immune_geneEffect_vs_geneEffect.csv`, `immune_headline.csv`. "
        "Figure: `fig3_immune.png`.",
        "",
        "A negative r in test B means higher immune-gene expression accompanies a more "
        "negative (more dependent) Chronos score. Given that neither TACSTD2 nor CLDN4 "
        "is a CRISPR dependency in these lines, large effects in B are not expected.",
    ]
    (NOTES / "04_immune.md").write_text("\n".join(lines) + "\n")

    # print headline composites for the operator
    show = headline[
        (headline.partner.isin(["MHC-I score", "IFN-response score"]))
        & (headline.cohort.isin(["lung", "NSCLC", "lung_partial_NSCLC"]))
    ]
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
