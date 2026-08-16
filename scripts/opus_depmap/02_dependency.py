#!/usr/bin/env python3
"""Is TACSTD2 or CLDN4 a CRISPR dependency in lung cell lines?

Tests performed
---------------
1. Location of the gene-effect distribution relative to the two DepMap control
   gene sets (Hart non-essentials and Hart/Blomen common essentials).
2. One-sample Wilcoxon test of lung gene effect against 0 (the non-essential
   anchor of the Chronos scale), with the effect size expressed in units of the
   non-essential control spread.
3. Lung vs non-lung Mann-Whitney with Cliff's delta, and the genome-wide
   percentile of that effect size so the lung comparison is calibrated against
   all ~18k genes rather than read in isolation.
4. Counts of dependent lines using DepMap's probability > 0.5 convention.
5. Whether expression of the gene predicts its own gene effect (i.e. is there a
   high-expressing subset that does depend on it?).

Writes tables to results/opus_depmap/tables/ and a note to notes/opus_depmap/.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DATA,
    FIGURES,
    GENES_OF_INTEREST,
    NOTES,
    RELEASE,
    TABLES,
    bh_fdr,
    ensure_dirs,
    fisher_ci,
    load_expression,
    load_gene_dependency,
    load_gene_effect,
    load_models,
    lung_cohort,
    mannwhitney,
    strip_entrez,
)

# Reference genes used for orientation: strong, well known lung-lineage
# dependencies plus a pan-essential and a neutral control.
BENCHMARKS = ["RPL7", "EEF2", "KRAS", "EGFR", "SOX2", "NKX2-1", "TP53", "OR2T1"]


def control_sets() -> tuple[list[str], list[str]]:
    ess = pd.read_csv(DATA / "AchillesCommonEssentialControls.csv")
    non = pd.read_csv(DATA / "AchillesNonessentialControls.csv")
    ess_genes = strip_entrez(ess.iloc[:, 0].tolist())
    non_genes = strip_entrez(non.iloc[:, 0].tolist())
    return ess_genes, non_genes


def describe(values: pd.Series) -> dict:
    v = values.dropna().astype(float)
    return {
        "n": int(v.size),
        "mean": float(v.mean()),
        "sd": float(v.std(ddof=1)) if v.size > 1 else np.nan,
        "min": float(v.min()) if v.size else np.nan,
        "q05": float(v.quantile(0.05)) if v.size else np.nan,
        "median": float(v.median()) if v.size else np.nan,
        "q95": float(v.quantile(0.95)) if v.size else np.nan,
        "max": float(v.max()) if v.size else np.nan,
    }


def main() -> int:
    ensure_dirs()
    models = load_models()
    lung = lung_cohort(models)
    gene_effect = load_gene_effect()
    gene_dep = load_gene_dependency()
    expr = load_expression()

    lung_ids = lung.index.intersection(gene_effect.index)
    other_ids = gene_effect.index.difference(lung.index)
    nsclc_ids = lung.index[lung.LungGroup == "NSCLC"].intersection(gene_effect.index)
    sclc_ids = lung.index[lung.LungGroup == "SCLC"].intersection(gene_effect.index)

    ess_genes, non_genes = control_sets()
    ess_present = [g for g in ess_genes if g in gene_effect.columns]
    non_present = [g for g in non_genes if g in gene_effect.columns]

    # Non-essential control spread in lung lines: the technical scale against
    # which a "small" gene effect should be judged.
    non_lung = gene_effect.loc[lung_ids, non_present]
    non_gene_means = non_lung.mean(axis=0)
    non_null_sd = float(non_gene_means.std(ddof=1))
    non_null_mean = float(non_gene_means.mean())
    ess_gene_means = gene_effect.loc[lung_ids, ess_present].mean(axis=0)

    # ------------------------------------------------------------------
    # 1-4: per-gene summary
    # ------------------------------------------------------------------
    rows = []
    for gene in GENES_OF_INTEREST + BENCHMARKS:
        if gene not in gene_effect.columns:
            continue
        vals = gene_effect[gene]
        lung_v = vals.loc[lung_ids].dropna()
        other_v = vals.loc[other_ids].dropna()
        row = {"gene": gene}
        row.update({f"lung_{k}": v for k, v in describe(lung_v).items()})
        row.update({f"pan_{k}": v for k, v in describe(vals).items()})
        row["nsclc_mean"] = float(vals.loc[nsclc_ids].mean())
        row["sclc_mean"] = float(vals.loc[sclc_ids].mean())

        # Wilcoxon signed rank of lung values against 0
        if lung_v.size > 10:
            row["lung_wilcoxon_vs0_p"] = float(
                stats.wilcoxon(lung_v.to_numpy(), alternative="two-sided").pvalue
            )
        else:
            row["lung_wilcoxon_vs0_p"] = np.nan
        # Effect size in units of the non-essential control gene spread
        row["lung_mean_z_vs_nonessential_controls"] = (
            float(lung_v.mean()) - non_null_mean
        ) / non_null_sd

        p, delta, shift = mannwhitney(lung_v.to_numpy(), other_v.to_numpy())
        row["lung_vs_other_p"] = p
        row["lung_vs_other_cliffs_delta"] = delta
        row["lung_vs_other_median_shift"] = shift

        p2, delta2, shift2 = mannwhitney(
            vals.loc[nsclc_ids].dropna().to_numpy(), vals.loc[sclc_ids].dropna().to_numpy()
        )
        row["nsclc_vs_sclc_p"] = p2
        row["nsclc_vs_sclc_cliffs_delta"] = delta2

        dep = gene_dep[gene] if gene in gene_dep.columns else pd.Series(dtype=float)
        row["lung_n_dependent_prob_gt_0.5"] = (
            int((dep.loc[lung_ids] > 0.5).sum()) if not dep.empty else np.nan
        )
        row["pan_n_dependent_prob_gt_0.5"] = int((dep > 0.5).sum()) if not dep.empty else np.nan
        row["lung_max_dep_prob"] = float(dep.loc[lung_ids].max()) if not dep.empty else np.nan
        row["lung_n_strongly_negative_lt_-0.5"] = int((lung_v < -0.5).sum())
        rows.append(row)

    summary = pd.DataFrame(rows).set_index("gene")

    # BH correction across the tests that are actually part of the primary
    # question (the two target genes), keeping benchmarks as context only.
    focus = summary.loc[GENES_OF_INTEREST]
    for col in ["lung_wilcoxon_vs0_p", "lung_vs_other_p", "nsclc_vs_sclc_p"]:
        summary.loc[GENES_OF_INTEREST, col.replace("_p", "_q")] = bh_fdr(
            focus[col].to_numpy()
        )
    summary.to_csv(TABLES / "dependency_summary.csv")

    # ------------------------------------------------------------------
    # Genome-wide calibration
    # ------------------------------------------------------------------
    lung_means = gene_effect.loc[lung_ids].mean(axis=0)
    lung_mins = gene_effect.loc[lung_ids].min(axis=0)
    dep_counts = (gene_dep.loc[lung_ids] > 0.5).sum(axis=0)
    calib_rows = []
    for gene in GENES_OF_INTEREST:
        calib_rows.append(
            {
                "gene": gene,
                "lung_mean_gene_effect": float(lung_means[gene]),
                "percentile_of_mean_among_all_genes": float(
                    stats.percentileofscore(lung_means.dropna(), lung_means[gene])
                ),
                "n_genes_more_negative_mean": int((lung_means < lung_means[gene]).sum()),
                "lung_min_gene_effect": float(lung_mins[gene]),
                "percentile_of_min_among_all_genes": float(
                    stats.percentileofscore(lung_mins.dropna(), lung_mins[gene])
                ),
                "n_dependent_lung_lines": int(dep_counts[gene]),
                "percentile_of_dependent_count": float(
                    stats.percentileofscore(dep_counts, dep_counts[gene])
                ),
                "genes_tested": int(lung_means.notna().sum()),
            }
        )
    calibration = pd.DataFrame(calib_rows)
    calibration.to_csv(TABLES / "dependency_genomewide_calibration.csv", index=False)

    # ------------------------------------------------------------------
    # 5: does expression of the gene predict its own gene effect?
    # ------------------------------------------------------------------
    expr_rows = []
    shared = lung_ids.intersection(expr.index)
    for gene in GENES_OF_INTEREST:
        e = expr.loc[shared, gene].astype(float)
        g = gene_effect.loc[shared, gene].astype(float)
        ok = e.notna() & g.notna()
        r, p = stats.spearmanr(e[ok], g[ok])
        lo, hi = fisher_ci(float(r), int(ok.sum()))
        hi_q = e[ok] >= e[ok].quantile(0.75)
        lo_q = e[ok] <= e[ok].quantile(0.25)
        pmw, delta, shift = mannwhitney(g[ok][hi_q].to_numpy(), g[ok][lo_q].to_numpy())
        expr_rows.append(
            {
                "gene": gene,
                "n": int(ok.sum()),
                "spearman_r_expr_vs_geneEffect": float(r),
                "ci_low": lo,
                "ci_high": hi,
                "p": float(p),
                "mean_geneEffect_top_expr_quartile": float(g[ok][hi_q].mean()),
                "mean_geneEffect_bottom_expr_quartile": float(g[ok][lo_q].mean()),
                "top_vs_bottom_quartile_p": pmw,
                "top_vs_bottom_cliffs_delta": delta,
            }
        )
    expr_dep = pd.DataFrame(expr_rows)
    expr_dep["q"] = bh_fdr(expr_dep["p"].to_numpy())
    expr_dep.to_csv(TABLES / "dependency_vs_expression.csv", index=False)

    # ------------------------------------------------------------------
    # Figure: gene effect distributions against controls
    # ------------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    ax.hist(ess_gene_means, bins=60, alpha=0.55, label=f"common essentials (n={len(ess_present)})", color="#b2182b")
    ax.hist(non_gene_means, bins=60, alpha=0.55, label=f"non-essential controls (n={len(non_present)})", color="#4393c3")
    for gene, colour in zip(GENES_OF_INTEREST, ["#1a9850", "#762a83"]):
        ax.axvline(lung_means[gene], color=colour, lw=2, label=f"{gene} (mean {lung_means[gene]:.3f})")
    ax.set_xlabel("mean Chronos gene effect across lung lines")
    ax.set_ylabel("number of genes")
    ax.set_title(f"Lung lines (n={len(lung_ids)}): control calibration")
    ax.legend(fontsize=8)

    ax = axes[1]
    data, labels = [], []
    for gene in GENES_OF_INTEREST + ["EGFR", "KRAS", "EEF2"]:
        data.append(gene_effect.loc[lung_ids, gene].dropna().to_numpy())
        labels.append(gene)
    parts = ax.violinplot(data, showmedians=True, widths=0.8)
    for pc in parts["bodies"]:
        pc.set_alpha(0.6)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=20)
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.axhline(-1, color="#b2182b", lw=0.8, ls=":")
    ax.set_ylabel("Chronos gene effect")
    ax.set_title("Per-line gene effect in lung lines")
    fig.suptitle(f"{RELEASE} · CRISPR dependency of TACSTD2 and CLDN4", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_dependency_distributions.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Note
    # ------------------------------------------------------------------
    lines = [
        "# 02 · Dependency of TACSTD2 and CLDN4 in lung lines",
        "",
        f"Release {RELEASE}. Lung models with CRISPR screens: **{len(lung_ids)}** "
        f"(NSCLC {len(nsclc_ids)}, SCLC {len(sclc_ids)}); non-lung comparison set: "
        f"{len(other_ids)} models.",
        "",
        "Chronos scale reminder: 0 = median non-essential gene, -1 = median common "
        "essential gene. In these lung lines the mean gene effect of the "
        f"{len(non_present)} non-essential control genes is {non_null_mean:.4f} with a "
        f"between-gene SD of {non_null_sd:.4f}; the {len(ess_present)} common essential "
        f"controls average {ess_gene_means.mean():.3f}. Those two numbers set the scale "
        "for everything below.",
        "",
        "## Headline numbers",
        "",
        "| gene | lung mean (SD) | lung min | dependent lung lines (p>0.5) | dependent lines pan-cancer | mean in non-essential SD units |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for gene in GENES_OF_INTEREST:
        r = summary.loc[gene]
        lines.append(
            f"| {gene} | {r['lung_mean']:.4f} ({r['lung_sd']:.3f}) | {r['lung_min']:.3f} | "
            f"{int(r['lung_n_dependent_prob_gt_0.5'])} / {len(lung_ids)} | "
            f"{int(r['pan_n_dependent_prob_gt_0.5'])} / {gene_effect.shape[0]} | "
            f"{r['lung_mean_z_vs_nonessential_controls']:+.2f} |"
        )
    lines += [
        "",
        "Full statistics, including the benchmark genes, are in "
        "`results/opus_depmap/tables/dependency_summary.csv`; genome-wide calibration is "
        "in `dependency_genomewide_calibration.csv`; the expression-versus-own-dependency "
        "test is in `dependency_vs_expression.csv`.",
        "",
        "Figure: `results/opus_depmap/figures/fig1_dependency_distributions.png`.",
    ]
    (NOTES / "02_dependency.md").write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 200)
    print(summary[
        [
            "lung_n",
            "lung_mean",
            "lung_sd",
            "lung_min",
            "lung_n_dependent_prob_gt_0.5",
            "pan_n_dependent_prob_gt_0.5",
            "lung_mean_z_vs_nonessential_controls",
            "lung_vs_other_p",
            "lung_vs_other_cliffs_delta",
        ]
    ].to_string())
    print()
    print(calibration.to_string(index=False))
    print()
    print(expr_dep.to_string(index=False))
    print()
    print(f"non-essential control between-gene SD in lung: {non_null_sd:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
