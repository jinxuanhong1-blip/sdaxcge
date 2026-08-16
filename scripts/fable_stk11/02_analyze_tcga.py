#!/usr/bin/env python3
"""Genotype <-> expression axis (TCGA).

Question: do TACSTD2 (TROP2) and CLDN4 (Claudin-4) mRNA levels differ by
STK11 / KEAP1 / KRAS mutation status in lung cancer?

Primary cohort: TCGA-LUAD (adenocarcinoma). TCGA-LUSC reported as a secondary
check (STK11/KRAS are rare in LUSC so only KEAP1 is informative there).
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fable_common as fc
import fable_stats as fs

EXPR_GENES = fc.TARGET_GENES + ["CD274", "CD8A"]


def load(label):
    expr = pd.read_csv(os.path.join(fc.PROC, f"{label}_expression.csv"))
    geno = pd.read_csv(os.path.join(fc.PROC, f"{label}_genotype.csv"))
    df = expr.merge(geno, on="sampleId", how="inner")
    # log2(RSEM + 1) for expression genes
    for g in EXPR_GENES:
        if g in df.columns:
            df[f"log2_{g}"] = np.log2(df[g].clip(lower=0) + 1)
    return df


def analyze_cohort(label):
    df = load(label)
    rows = []
    for gt in fc.GENOTYPE_GENES:
        flag = f"{gt}_mut"
        for tgt in fc.TARGET_GENES:
            col = f"log2_{tgt}"
            mut = df.loc[df[flag] == 1, col]
            wt = df.loc[df[flag] == 0, col]
            r = fs.group_compare(mut, wt)
            r.update({"cohort": label, "genotype": gt, "target": tgt})
            rows.append(r)
    res = pd.DataFrame(rows)
    res["fdr_bh"] = fs.bh_fdr(res["mwu_p"].fillna(1.0))
    return df, res


def comutation_analysis(df, label):
    """Within KRAS-mutant tumors, effect of STK11 (KL) or KEAP1 co-mutation."""
    rows = []
    kras = df[df.KRAS_mut == 1]
    for co in ["STK11", "KEAP1"]:
        for tgt in fc.TARGET_GENES:
            col = f"log2_{tgt}"
            mut = kras.loc[kras[f"{co}_mut"] == 1, col]
            wt = kras.loc[kras[f"{co}_mut"] == 0, col]
            r = fs.group_compare(mut, wt)
            r.update({"cohort": label, "background": "KRAS-mut",
                      "co_mutation": co, "target": tgt})
            rows.append(r)
    return pd.DataFrame(rows)


def immune_correlations(df, label):
    rows = []
    for tgt in fc.TARGET_GENES:
        for ctx in ["CD8A", "CD274"]:
            x = df[f"log2_{tgt}"]
            y = df[f"log2_{ctx}"]
            m = x.notna() & y.notna()
            rho, p = stats.spearmanr(x[m], y[m])
            rows.append({"cohort": label, "target": tgt, "context": ctx,
                         "spearman_rho": float(rho), "p": float(p), "n": int(m.sum())})
    return pd.DataFrame(rows)


def boxplot(df, label):
    fig, axes = plt.subplots(len(fc.TARGET_GENES), len(fc.GENOTYPE_GENES),
                             figsize=(11, 7), squeeze=False)
    for i, tgt in enumerate(fc.TARGET_GENES):
        for j, gt in enumerate(fc.GENOTYPE_GENES):
            ax = axes[i][j]
            col = f"log2_{tgt}"
            flag = f"{gt}_mut"
            wt = df.loc[df[flag] == 0, col].dropna()
            mut = df.loc[df[flag] == 1, col].dropna()
            ax.boxplot([wt, mut], showfliers=False,
                       tick_labels=[f"WT\n(n={wt.size})", f"MUT\n(n={mut.size})"])
            for k, vals in enumerate([wt, mut], start=1):
                x = np.random.normal(k, 0.05, size=vals.size)
                ax.scatter(x, vals, s=6, alpha=0.35,
                           color="#377eb8" if k == 1 else "#e41a1c")
            try:
                _, p = stats.mannwhitneyu(mut, wt, alternative="two-sided")
                ax.set_title(f"{tgt} by {gt}  (p={p:.1e})", fontsize=9)
            except ValueError:
                ax.set_title(f"{tgt} by {gt}", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"log2 {tgt}")
    fig.suptitle(f"{label}: ADC-target expression by genotype", fontsize=12)
    fig.tight_layout()
    out = os.path.join(fc.FIG, f"tcga_{label}_expr_by_genotype.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    all_res, all_co, all_corr = [], [], []
    for label in ["TCGA-LUAD", "TCGA-LUSC"]:
        df, res = analyze_cohort(label)
        all_res.append(res)
        all_co.append(comutation_analysis(df, label))
        all_corr.append(immune_correlations(df, label))
        fig = boxplot(df, label)
        print(f"[{label}] n={df.shape[0]} figure -> {os.path.relpath(fig, fc.REPO)}")

    res = pd.concat(all_res, ignore_index=True)
    co = pd.concat(all_co, ignore_index=True)
    corr = pd.concat(all_corr, ignore_index=True)
    res.to_csv(os.path.join(fc.TAB, "tcga_genotype_vs_expression.csv"), index=False)
    co.to_csv(os.path.join(fc.TAB, "tcga_kras_comutation_vs_expression.csv"), index=False)
    corr.to_csv(os.path.join(fc.TAB, "tcga_target_immune_correlation.csv"), index=False)

    cols = ["cohort", "genotype", "target", "n_mut", "n_wt", "median_mut",
            "median_wt", "log2fc_median", "cliffs_delta", "mwu_p", "fdr_bh"]
    print("\n=== TCGA genotype vs ADC-target expression (log2 RSEM) ===")
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(res[cols].round(4).to_string(index=False))
    print("\n=== KRAS-mutant co-mutation (STK11/KEAP1) vs expression ===")
    print(co[["cohort", "co_mutation", "target", "n_mut", "n_wt",
              "log2fc_median", "cliffs_delta", "mwu_p"]].round(4).to_string(index=False))
    print("\n=== Target vs immune-context correlation (Spearman) ===")
    print(corr.round(4).to_string(index=False))


if __name__ == "__main__":
    np.random.seed(0)
    main()
