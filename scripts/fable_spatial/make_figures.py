"""Figures for the fable_spatial slice (saved to results/fable_spatial/)."""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import ensure_results

OUT = ensure_results()
DATA = "/workspace/data/fable_spatial/gse271689"


def fig_visium_rhos():
    em = pd.read_csv(os.path.join(OUT, "emtab13530_per_section.csv"))
    gs = pd.read_csv(os.path.join(OUT, "gse189487_per_section.csv"))
    gs["tissue"] = "LUAD (" + gs["stage"] + ")"
    em["dataset"] = "E-MTAB-13530"
    gs["dataset"] = "GSE189487"
    df = pd.concat([em, gs])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, nb, title in zip(axes, ["immune_broad", "t_effector"],
                             ["broad immune neighborhood", "T-effector neighborhood"]):
        sub = df[df.neighborhood == nb]
        groups, labels = [], []
        for (ds, tis, gene), g in sub.groupby(["dataset", "tissue", "gene"]):
            groups.append(g.partial_rho.values)
            labels.append(f"{ds[:12]} {tis} {gene} (n={len(g)})")
        pos = np.arange(len(groups))
        for i, vals in enumerate(groups):
            ax.scatter(vals, np.full(len(vals), pos[i]), s=14, alpha=0.65,
                       c=["#c0392b" if v > 0 else "#2980b9" for v in vals])
            ax.plot([np.median(vals)] * 2, [pos[i] - 0.3, pos[i] + 0.3], c="k", lw=2)
        ax.axvline(0, c="grey", ls="--", lw=1)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("partial Spearman rho (per section)")
        ax.set_title(title, fontsize=10)
    fig.suptitle("Epithelial-spot TACSTD2/CLDN4 vs neighborhood immune score (Visium)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_visium_neighborhood_rhos.png"), dpi=200)
    plt.close(fig)


def fig_yale_km():
    from lifelines import KaplanMeierFitter
    xl = pd.ExcelFile(os.path.join(DATA, "moesm10.xlsx"))
    yt = pd.read_excel(xl, "source_data_Figure_6b")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
        d = yt[["OS_Days", "OS_Index", gene]].dropna()
        d = d[d.OS_Days > 0]
        expr = np.log2(d[gene].astype(float) + 1)
        hi = expr > expr.median()
        for mask, lab, col in [(hi, "high", "#c0392b"), (~hi, "low", "#2980b9")]:
            km = KaplanMeierFitter().fit(d.OS_Days[mask], d.OS_Index[mask],
                                         label=f"{gene} {lab} (n={mask.sum()})")
            km.plot_survival_function(ax=ax, ci_show=True, color=col)
        ax.set_xlabel("days")
        ax.set_ylabel("OS probability")
        ax.set_title(f"Yale WTA tumor compartment: {gene}", fontsize=10)
    fig.suptitle("GSE271689 Yale first-line ICI cohort, median split", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_gse271689_yale_km.png"), dpi=200)
    plt.close(fig)


def fig_tumour_vs_normal():
    ex = pd.read_csv(os.path.join(OUT, "emtab13530_epi_expression.csv"))
    pt = ex.groupby(["patient", "tissue"])[["TACSTD2_mean_epi", "CLDN4_mean_epi"]].mean().reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    order = ["donor", "non_involved", "tumour"]
    for ax, g in zip(axes, ["TACSTD2", "CLDN4"]):
        vals = [pt[pt.tissue == t][f"{g}_mean_epi"].values for t in order]
        ax.boxplot(vals, tick_labels=order, widths=0.5)
        for i, v in enumerate(vals):
            ax.scatter(np.random.default_rng(1).normal(i + 1, 0.05, len(v)), v, s=18, alpha=0.7)
        ax.set_ylabel(f"{g} (mean log1p CP10K, epithelial spots)")
        ax.set_title(g)
    fig.suptitle("E-MTAB-13530: epithelial-spot expression by tissue", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_emtab13530_expression_by_tissue.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    fig_visium_rhos()
    fig_yale_km()
    fig_tumour_vs_normal()
    print("figures written")
