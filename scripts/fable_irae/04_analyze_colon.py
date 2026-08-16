#!/usr/bin/env python3
"""GSE206300 - ICI colitis colon *epithelial* single-nucleus RNA-seq.

Primary analysis for the epithelial markers TACSTD2 (TROP2) and CLDN4 (claudin-4).
Design: 12 irColitis "Case" donors vs 14 "Control" donors (single-nucleus, epithelial
compartment). X = raw UMI counts.

Approach (standard donor-level pseudobulk):
  1. For each gene of interest + housekeeping/epithelial controls, pull the raw count
     column (memory-safe backed slice).
  2. Per donor: pseudobulk CP10K = 1e4 * sum(gene counts) / sum(all counts), then log1p.
     Also per-cell detection rate (fraction of nuclei with >=1 UMI).
  3. Compare Case vs Control across donors with Mann-Whitney U (two-sided) and report
     log2 fold-change of mean CP10K, Cliff's delta effect size.

Outputs (small, committed):
  results/fable_irae/tables/colon_GSE206300_pseudobulk.tsv
  results/fable_irae/tables/colon_GSE206300_de_summary.tsv
  results/fable_irae/figures/colon_GSE206300_*.png
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "results" / "fable_irae"
RAW = ROOT / "data" / "raw" / "GSE206300_ircolitis-tissue-epithelial.h5ad"
TAB = ROOT / "tables"; TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"; FIG.mkdir(parents=True, exist_ok=True)

GOI = ["TACSTD2", "CLDN4"]
CONTROLS = ["EPCAM", "CLDN3", "CLDN7", "ACTB", "PTPRC"]  # epithelial + housekeeping + immune-neg
ALL_GENES = GOI + CONTROLS


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def main():
    a = ad.read_h5ad(RAW, backed="r")
    syms = np.array([str(v).split("|")[-1] for v in a.var_names])
    gene_idx = {}
    for g in ALL_GENES:
        hit = np.where(syms == g)[0]
        if len(hit):
            gene_idx[g] = int(hit[0])
        else:
            print(f"  WARNING gene not found: {g}")

    obs = a.obs[["donor", "case", "class2", "drug", "n_counts"]].copy()
    obs["n_counts"] = obs["n_counts"].astype(float)

    # pull gene count columns (dense small matrix cells x len(genes))
    idx_list = [gene_idx[g] for g in gene_idx]
    sub = a[:, idx_list].to_memory()
    X = sub.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    counts = pd.DataFrame(X, columns=list(gene_idx.keys()), index=obs.index)
    a.file.close()

    counts = counts.join(obs)

    # donor-level pseudobulk
    donor_meta = counts.groupby("donor", observed=True).agg(
        case=("case", "first"), class2=("class2", "first"), drug=("drug", "first"),
        n_cells=("n_counts", "size"), total_counts=("n_counts", "sum"))
    rows = []
    for donor, g in counts.groupby("donor", observed=True):
        tot = g["n_counts"].sum()
        rec = {"donor": donor}
        for gene in gene_idx:
            gsum = g[gene].sum()
            rec[f"{gene}_cp10k"] = 1e4 * gsum / tot
            rec[f"{gene}_log1p_cp10k"] = np.log1p(1e4 * gsum / tot)
            rec[f"{gene}_detrate"] = float((g[gene] > 0).mean())
        rows.append(rec)
    pdb = pd.DataFrame(rows).set_index("donor").join(donor_meta)
    pdb.to_csv(TAB / "colon_GSE206300_pseudobulk.tsv", sep="\t")
    print(f"[write] colon pseudobulk ({pdb.shape[0]} donors)")

    # DE Case vs Control across donors
    de_rows = []
    case = pdb[pdb["case"] == "Case"]
    ctrl = pdb[pdb["case"] == "Control"]
    for gene in gene_idx:
        col = f"{gene}_log1p_cp10k"
        cv, kv = case[col].values, ctrl[col].values
        u, p = stats.mannwhitneyu(cv, kv, alternative="two-sided")
        # log2FC on mean CP10K (linear)
        mean_case = case[f"{gene}_cp10k"].mean()
        mean_ctrl = ctrl[f"{gene}_cp10k"].mean()
        l2fc = np.log2((mean_case + 1e-6) / (mean_ctrl + 1e-6))
        de_rows.append(dict(
            gene=gene, group="Case_vs_Control",
            n_case=len(cv), n_ctrl=len(kv),
            mean_cp10k_case=mean_case, mean_cp10k_ctrl=mean_ctrl,
            log2FC_case_over_ctrl=l2fc,
            detrate_case=case[f"{gene}_detrate"].mean(),
            detrate_ctrl=ctrl[f"{gene}_detrate"].mean(),
            mwu_U=u, p_value=p, cliffs_delta=cliffs_delta(cv, kv)))
    de = pd.DataFrame(de_rows)
    # BH FDR across the two genes of interest only
    from statsmodels.stats.multitest import multipletests
    goi_mask = de["gene"].isin(GOI)
    de.loc[goi_mask, "fdr_bh"] = multipletests(de.loc[goi_mask, "p_value"], method="fdr_bh")[1]
    de.to_csv(TAB / "colon_GSE206300_de_summary.tsv", sep="\t", index=False)
    print("[write] colon DE summary")
    print(de[["gene", "mean_cp10k_case", "mean_cp10k_ctrl", "log2FC_case_over_ctrl",
              "p_value", "cliffs_delta"]].to_string(index=False))

    # Figures: per-donor boxplot for GOI
    fig, axes = plt.subplots(1, len(GOI), figsize=(4 * len(GOI), 4))
    for ax, gene in zip(np.atleast_1d(axes), GOI):
        col = f"{gene}_log1p_cp10k"
        data = [ctrl[col].values, case[col].values]
        ax.boxplot(data, tick_labels=["Control", "irColitis"], showfliers=False)
        for i, d in enumerate(data, 1):
            ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=18, alpha=0.7, color="#c0392b" if i == 2 else "#2c7fb8")
        p = de.loc[de.gene == gene, "p_value"].values[0]
        ax.set_title(f"{gene} (colon epithelium)\nMWU p={p:.3g}")
        ax.set_ylabel("pseudobulk log1p(CP10K)")
    fig.suptitle("GSE206300 ICI colitis — epithelial TACSTD2 / CLDN4 (donor-level)")
    fig.tight_layout()
    fig.savefig(FIG / "colon_GSE206300_boxplots.png", dpi=140)
    print("[write] colon figure")


if __name__ == "__main__":
    main()
