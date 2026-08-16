"""Figures for the George bulk analysis."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from . import bulk_george as bg
from . import signatures as sig

FIGS = bg.FIGS
SUBTYPE_ORDER = ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]
SUBTYPE_COLORS = {"SCLC-A": "#4C72B0", "SCLC-N": "#DD8452",
                  "SCLC-P": "#55A868", "SCLC-I": "#C44E52"}


def _annot_p(ax, p, x=0.98, y=0.02):
    ax.text(x, y, f"p={p:.3g}", ha="right", va="bottom",
            transform=ax.transAxes, fontsize=8, color="0.3")


def fig_box_by_subtype(expr, subtype):
    genes = ["TACSTD2", "CLDN4"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for ax, gene in zip(axes, genes):
        data, labels = [], []
        for s in SUBTYPE_ORDER:
            v = expr.loc[gene, subtype[subtype == s].index].astype(float).values
            if v.size:
                data.append(v)
                labels.append(f"{s}\n(n={v.size})")
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.6,
                        showfliers=False)
        for patch, s in zip(bp["boxes"], SUBTYPE_ORDER):
            patch.set_facecolor(SUBTYPE_COLORS[s])
            patch.set_alpha(0.55)
        for i, v in enumerate(data, 1):
            ax.scatter(np.random.normal(i, 0.06, size=len(v)), v, s=12,
                       color="0.2", alpha=0.6, zorder=3)
        # SCLC-I vs rest MWU
        iv = expr.loc[gene, subtype[subtype == "SCLC-I"].index].astype(float).values
        rv = expr.loc[gene, subtype[subtype != "SCLC-I"].index].astype(float).values
        if iv.size >= 3 and rv.size >= 3:
            _, p = stats.mannwhitneyu(iv, rv, alternative="two-sided")
            _annot_p(ax, p)
        ax.set_title(f"{gene} expression by SCLC subtype")
        ax.set_ylabel("log2(FPKM+1)")
    fig.suptitle("George et al. 2015 (n=81) — ADC targets across subtypes\n"
                 "p = SCLC-I vs all-other (Mann-Whitney)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIGS / "george_box_TACSTD2_CLDN4_by_subtype.png", dpi=150)
    plt.close(fig)


def fig_scatter_vs_immune(expr, scores):
    pairs = [("TACSTD2", "GEP_Ayers_Tcell_inflamed"),
             ("TACSTD2", "CD8A_log2"),
             ("CLDN4", "GEP_Ayers_Tcell_inflamed"),
             ("CLDN4", "CD8A_log2")]
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for ax, (gene, axis) in zip(axes.ravel(), pairs):
        x = scores[axis].astype(float)
        y = expr.loc[gene].reindex(x.index).astype(float)
        c = scores["subtype"].map(SUBTYPE_COLORS)
        ax.scatter(x, y, c=c, s=28, alpha=0.8, edgecolor="0.3", linewidth=0.4)
        rho, p = stats.spearmanr(x, y)
        # trend line
        b, a = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, a + b * xs, color="0.3", lw=1, ls="--")
        ax.set_xlabel(axis.replace("_", " "))
        ax.set_ylabel(f"{gene} log2(FPKM+1)")
        ax.set_title(f"{gene} vs {axis.replace('_',' ')}\nSpearman rho={rho:.2f}, p={p:.3g}",
                     fontsize=9)
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=SUBTYPE_COLORS[s],
                          label=s) for s in SUBTYPE_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("George et al. 2015 — ADC targets vs immune axis", fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(FIGS / "george_scatter_ADC_vs_immune.png", dpi=150)
    plt.close(fig)


def fig_corr_heatmap(corr):
    axes_keep = ["GEP_Ayers_Tcell_inflamed", "CYT_log2", "CD8A_log2",
                 "HLA_antigen_presentation", "Immune_checkpoints",
                 "STING_chemokines", "EMT_mesenchymal", "Epithelial",
                 "NE_score", "axis_A", "axis_N", "axis_P", "axis_I"]
    genes = ["TACSTD2", "CLDN4"]
    mat = pd.DataFrame(index=genes, columns=axes_keep, dtype=float)
    qmat = pd.DataFrame(index=genes, columns=axes_keep, dtype=float)
    for _, r in corr.iterrows():
        if r["gene"] in genes and r["axis"] in axes_keep:
            mat.loc[r["gene"], r["axis"]] = r["spearman_rho"]
            qmat.loc[r["gene"], r["axis"]] = r["q_value_BH"]
    fig, ax = plt.subplots(figsize=(11, 3.2))
    im = ax.imshow(mat.values.astype(float), cmap="RdBu_r", vmin=-0.6, vmax=0.6,
                   aspect="auto")
    ax.set_xticks(range(len(axes_keep)))
    ax.set_xticklabels([a.replace("_", " ") for a in axes_keep], rotation=45,
                       ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    for i, g in enumerate(genes):
        for j, a in enumerate(axes_keep):
            rho = mat.loc[g, a]
            q = qmat.loc[g, a]
            star = "*" if (pd.notna(q) and q < 0.05) else ""
            ax.text(j, i, f"{rho:.2f}{star}", ha="center", va="center",
                    fontsize=7, color="black")
    fig.colorbar(im, ax=ax, label="Spearman rho", fraction=0.02, pad=0.01)
    ax.set_title("ADC target vs immune / lineage axes (George 2015). * = BH-FDR q<0.05",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "george_corr_heatmap.png", dpi=150)
    plt.close(fig)


def run():
    res = bg.run()
    expr, scores, subtype, corr = res["expr"], res["scores"], res["subtype"], res["corr"]
    fig_box_by_subtype(expr, subtype)
    fig_scatter_vs_immune(expr, scores)
    fig_corr_heatmap(corr)
    print("Figures written to", FIGS)


if __name__ == "__main__":
    run()
