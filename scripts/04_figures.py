"""Figures for Claim A10. Reads only tables already written under results/claim_A10/."""

from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _load(name):
    p = os.path.join(C.TABLES, name)
    if not os.path.exists(p):
        return None
    return pd.read_csv(p)


def fig_pairwise_heatmaps():
    df = _load("human_bulk_pairwise_spearman.csv")
    if df is None:
        return
    genes = C.CLAIM_GENES_H
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, ds in zip(axes, ["GTEx-lung", "TCGA-LUAD"]):
        sub = df[df.dataset == ds]
        M = pd.DataFrame(np.nan, index=genes, columns=genes)
        for _, r in sub.iterrows():
            if r.gene_a in genes and r.gene_b in genes:
                M.loc[r.gene_a, r.gene_b] = r.spearman_rho
                M.loc[r.gene_b, r.gene_a] = r.spearman_rho
        np.fill_diagonal(M.values, 1.0)
        im = ax.imshow(M.to_numpy(), vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_xticks(range(len(genes)))
        ax.set_yticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=45, ha="right")
        ax.set_yticklabels(genes)
        ax.set_title(ds)
        for i in range(len(genes)):
            for j in range(len(genes)):
                v = M.iloc[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                            color="white" if abs(v) > 0.55 else "black")
    fig.colorbar(im, ax=axes, shrink=0.8, label="Spearman ρ")
    fig.suptitle("Claim-gene pairwise Spearman (human bulk lung)", y=1.02)
    fig.savefig(os.path.join(C.FIGURES, "fig1_human_bulk_pairwise.png"))
    plt.close(fig)


def fig_module_forest():
    df = _load("human_bulk_module_score_assoc.csv")
    if df is None:
        return
    keep = C.TARGET_H + [C.NKX_H] + ["EPCAM", "SFTPC", "SCGB1A1", "KRT5", "PTPRC"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = 0
    yticks, ylabels = [], []
    for g in keep:
        for ds, color in [("GTEx-lung", "#1f4e79"), ("TCGA-LUAD", "#b85c38")]:
            r = df[(df.y == g) & (df.dataset == ds)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.plot([r.ci_low, r.ci_high], [y, y], color=color, lw=1.6)
            ax.plot(r.spearman_rho, y, "o", color=color, ms=5)
            yticks.append(y)
            ylabels.append(f"{g}  {ds}")
            y -= 1
        y -= 0.35
    ax.axvline(0, color="0.5", lw=0.8)
    ax.axvline(0.3, color="0.8", ls="--", lw=0.7)
    ax.axvline(-0.3, color="0.8", ls="--", lw=0.7)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Spearman ρ of MODULE score vs gene (95% CI)")
    ax.set_title("P2/P3: module vs targets, NKX2-1, and composition markers")
    fig.savefig(os.path.join(C.FIGURES, "fig2_module_score_forest.png"))
    plt.close(fig)


def fig_partial():
    df = _load("human_bulk_partial_correlation.csv")
    if df is None:
        return
    pairs = [(a, b) for a in C.MODULE_H for b in C.TARGET_H]
    pairs += [(a, C.NKX_H) for a in C.MODULE_H + C.TARGET_H]
    cov_order = ["none", "EPCAM", "EPCAM+SFTPC+SCGB1A1", "all_lineage_markers"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), sharey=True)
    for ax, ds in zip(axes, ["GTEx-lung", "TCGA-LUAD"]):
        sub = df[df.dataset == ds]
        for i, (a, b) in enumerate(pairs):
            ys = []
            for cov in cov_order:
                r = sub[(sub.gene_a == a) & (sub.gene_b == b) & (sub.covariates == cov)]
                if r.empty:
                    r = sub[(sub.gene_a == b) & (sub.gene_b == a) & (sub.covariates == cov)]
                ys.append(float(r.iloc[0].partial_rho) if len(r) else np.nan)
            ax.plot(range(len(cov_order)), ys, "-o", ms=3, lw=0.9, alpha=0.85,
                    label=f"{a}–{b}")
        ax.axhline(0, color="0.5", lw=0.7)
        ax.set_xticks(range(len(cov_order)))
        ax.set_xticklabels(["none", "EPCAM", "+SFTPC\n+SCGB1A1", "all lineage"], fontsize=7)
        ax.set_title(ds)
        ax.set_ylabel("partial Spearman ρ")
    axes[1].legend(fontsize=6, ncol=2, loc="upper right", frameon=False)
    fig.suptitle("Do claim associations survive composition adjustment?")
    fig.savefig(os.path.join(C.FIGURES, "fig3_partial_correlation.png"))
    plt.close(fig)


def fig_celltype_means():
    df = _load("census_human_celltype_means.csv")
    if df is None or df.empty:
        return
    # Prefer normal tissue; fall back to whatever is there
    sub = df[df.disease == "normal"].copy()
    if sub.empty:
        sub = df.copy()
    genes = C.CLAIM_GENES_H
    # rank cell types by TACSTD2
    if "TACSTD2" in sub.columns:
        sub = sub.sort_values("TACSTD2", ascending=False)
    sub = sub.head(16)
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    x = np.arange(len(sub))
    width = 0.12
    colors = ["#4c78a8", "#f58518", "#54a24b", "#e45756", "#b279a2", "#72b7b2", "#ff9da6"]
    for i, g in enumerate(genes):
        if g not in sub.columns:
            continue
        ax.bar(x + (i - 3) * width, sub[g].to_numpy(), width, label=g, color=colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c} (n={n})" for c, n in zip(sub.cell_type, sub.n_cells)],
                       rotation=40, ha="right", fontsize=7)
    ax.set_ylabel("mean log1p(CP10k)")
    ax.legend(ncol=4, fontsize=7, frameon=False)
    ax.set_title("Human lung scRNA (Census): claim-gene means by epithelial type (normal)")
    fig.savefig(os.path.join(C.FIGURES, "fig4_census_human_celltype_means.png"))
    plt.close(fig)


def fig_within_type():
    df = _load("census_human_within_type_spearman.csv")
    if df is None or df.empty:
        return
    focus = [("ELF3", "TACSTD2"), ("ELF3", "CLDN4"), ("GRHL1", "TACSTD2"),
             ("KLF4", "TACSTD2"), ("TFAP2A", "TACSTD2"),
             ("NKX2-1", "TACSTD2"), ("NKX2-1", "CLDN4"), ("NKX2-1", "ELF3"),
             ("TACSTD2", "CLDN4")]
    types = ["basal cell", "club cell", "ciliated cell", "lung secretory cell",
             "pulmonary alveolar type 2 cell", "pulmonary alveolar type 1 cell",
             "malignant cell"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.0), sharey=True)
    for ax, dis in zip(axes, ["normal", "lung adenocarcinoma"]):
        sub = df[df.disease == dis]
        mat = pd.DataFrame(index=types, columns=[f"{a}–{b}" for a, b in focus], dtype=float)
        for a, b in focus:
            for ct in types:
                r = sub[(((sub.gene_a == a) & (sub.gene_b == b)) |
                         ((sub.gene_a == b) & (sub.gene_b == a))) & (sub.cell_type == ct)]
                if len(r):
                    mat.loc[ct, f"{a}–{b}"] = r.iloc[0].spearman_rho
        im = ax.imshow(mat.to_numpy(dtype=float), vmin=-0.6, vmax=0.6, cmap="RdBu_r",
                       aspect="auto")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(types)))
        ax.set_yticklabels(types, fontsize=8)
        ax.set_title(dis)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.iloc[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=axes, shrink=0.7, label="within-type Spearman ρ")
    fig.suptitle("P6: within-cell-type correlations (Census human lung epithelium)")
    fig.savefig(os.path.join(C.FIGURES, "fig5_census_within_type.png"))
    plt.close(fig)


def fig_perturbation():
    df = _load("perturbation_contrasts.csv")
    if df is None or df.empty:
        return
    gene_order = ["ELF3", "Elf3", "GRHL1", "Grhl1", "KLF4", "Klf4",
                  "TFAP2A", "Tfap2a", "TACSTD2", "Tacstd2", "CLDN4", "Cldn4",
                  "NKX2-1", "Nkx2-1"]
    genes = [g for g in gene_order if g in set(df.gene)]
    ds = list(df.dataset.unique())
    M = pd.DataFrame(index=ds, columns=genes, dtype=float)
    for _, r in df.iterrows():
        M.loc[r.dataset, r.gene] = r.log2fc_low_minus_high
    fig, ax = plt.subplots(figsize=(8.2, max(4.0, 0.32 * len(ds) + 1.5)))
    im = ax.imshow(M.to_numpy(dtype=float), vmin=-2, vmax=2, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(ds)))
    ax.set_yticklabels(ds, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iloc[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 1.1 else "black")
    fig.colorbar(im, ax=ax, shrink=0.7, label="log2FC (NKX-low − control)")
    ax.set_title("P4: NKX2-1 loss/low vs control (red = claim direction)")
    fig.savefig(os.path.join(C.FIGURES, "fig6_perturbation_log2fc.png"))
    plt.close(fig)


def fig_tcga_strata():
    df = _load("human_tcga_nkx2-1_low_vs_high.csv")
    if df is None or df.empty:
        return
    genes = C.CLAIM_GENES_H + ["EPCAM", "SFTPC", "SCGB1A1", "KRT5"]
    sub = df[df.gene.isin(genes)].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    colors = ["#4c78a8" if g in C.MODULE_H else "#b279a2" if g in C.TARGET_H
              else "#333333" if g == C.NKX_H else "0.55" for g in sub.gene]
    ax.barh(sub.gene[::-1], sub.log2fc_low_minus_high[::-1], color=colors[::-1])
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlabel("log2FC (NKX2-1-low tumours − NKX2-1-high tumours)")
    ax.set_title("TCGA-LUAD: claim genes in NKX2-1 quartile split")
    fig.savefig(os.path.join(C.FIGURES, "fig7_tcga_nkx_strata.png"))
    plt.close(fig)


def main():
    fig_pairwise_heatmaps()
    fig_module_forest()
    fig_partial()
    fig_celltype_means()
    fig_within_type()
    fig_perturbation()
    fig_tcga_strata()
    print("figures written to", C.FIGURES, flush=True)


if __name__ == "__main__":
    main()
