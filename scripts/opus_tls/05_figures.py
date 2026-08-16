"""Figures for the TACSTD2 vs TLS / B-cell / T-cell slice.

All numbers are read from the tables written by 03/04. No recomputation of
p-values here except for the scatter-plot annotations (which re-read the
cached scored cohorts).
"""

from __future__ import annotations

import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from opus_tls_lib import PROC_DIR, RESULTS_DIR, fmt_p

FIG = os.path.join(RESULTS_DIR, "figures")
TAB = os.path.join(RESULTS_DIR, "tables")
os.makedirs(FIG, exist_ok=True)

SIG_ORDER = ["TLS_Cabrita", "TLS_12chemokine", "TLS_imprint", "B_cell",
             "Plasma_cell", "Tfh", "T_cell_CD8", "IFNg_Ayers"]
COH_ORDER = ["TCGA-LUSC", "TCGA-LUAD", "GSE72094", "GSE81089",
             "GSE207422", "GSE135222", "GSE126044", "GSE190265"]


def _style() -> None:
    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 140, "savefig.bbox": "tight", "savefig.facecolor": "white",
    })


def fig_heatmap_adj() -> None:
    dec = pd.read_csv(os.path.join(TAB, "where_it_holds.csv"))
    sub = dec[(dec.gene == "TACSTD2") & (dec.signature.isin(SIG_ORDER))].copy()
    mat = sub.pivot(index="signature", columns="cohort", values="rho_adj")
    mat = mat.reindex(index=SIG_ORDER, columns=[c for c in COH_ORDER if c in mat.columns])
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    vmax = 0.45
    im = ax.imshow(mat.values.astype(float), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if not np.isfinite(v):
                continue
            row = sub[(sub.signature == mat.index[i]) & (sub.cohort == mat.columns[j])].iloc[0]
            star = "*" if row["holds_after_purity"] else ("+" if row["opposite_after_purity"] else "")
            ax.text(j, i, f"{v:.2f}{star}", ha="center", va="center", fontsize=7,
                    color="white" if abs(v) > 0.28 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="partial Spearman ρ")
    ax.set_title("TACSTD2 vs immune signatures after purity / epithelial correction\n"
                 "* holds (ρ<0, p<0.05, n≥40)   + opposite (ρ>0, p<0.05, n≥40)")
    fig.savefig(os.path.join(FIG, "fig1_heatmap_tacstd2_adj.png"))
    plt.close()


def fig_forest_lusc_vs_luad() -> None:
    corr = pd.read_csv(os.path.join(TAB, "correlations_by_cohort.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), sharey=True)
    for ax, coh, title in [
        (axes[0], "TCGA-LUSC", "TCGA-LUSC  (ABSOLUTE purity)"),
        (axes[1], "TCGA-LUAD", "TCGA-LUAD  (ABSOLUTE purity)"),
    ]:
        sub = corr[(corr.cohort == coh) & (corr.gene == "TACSTD2")
                   & (corr.signature.isin(SIG_ORDER))].set_index("signature").loc[SIG_ORDER]
        y = np.arange(len(SIG_ORDER))
        ax.errorbar(sub["rho"], y - 0.15, xerr=[sub["rho"] - sub["ci_low"], sub["ci_high"] - sub["rho"]],
                    fmt="o", color="#888888", ms=4, label="crude", capsize=2)
        ax.scatter(sub["rho_adj_purity"], y + 0.15, color="#b2182b", s=28, zorder=3, label="adj. ABSOLUTE")
        ax.axvline(0, color="k", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(SIG_ORDER)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title)
        ax.set_xlim(-0.45, 0.25)
        ax.legend(frameon=False, loc="lower left")
        for i, sig in enumerate(SIG_ORDER):
            p = sub.loc[sig, "p_adj_purity"]
            ax.text(0.22, i + 0.15, fmt_p(p), va="center", ha="right", fontsize=7, color="#b2182b")
    fig.suptitle("TACSTD2–immune association: crude vs ABSOLUTE-purity partial Spearman", y=1.02)
    fig.savefig(os.path.join(FIG, "fig2_forest_tcga.png"))
    plt.close()


def fig_scatter_lusc() -> None:
    with open(os.path.join(PROC_DIR, "scored_cohorts.pkl"), "rb") as fh:
        slim = pickle.load(fh)
    c = slim["TCGA-LUSC"]
    x = c["focus_expr"].loc["TACSTD2"]
    pur = c["pheno"]["purity_absolute"]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3))
    for ax, sig, ylab in [
        (axes[0], "TLS_Cabrita", "TLS (Cabrita)"),
        (axes[1], "B_cell", "B-cell signature"),
        (axes[2], "T_cell_CD8", "CD8 T-cell signature"),
    ]:
        y = c["scores"][sig]
        ok = x.notna() & y.notna() & pur.notna()
        sc = ax.scatter(x[ok], y[ok], c=pur[ok], cmap="viridis", s=10, alpha=0.75, vmin=0.2, vmax=0.95)
        rho, p = stats.spearmanr(x[ok], y[ok])
        ax.set_xlabel("TACSTD2  log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(f"crude ρ={rho:.2f}  p={fmt_p(p)}")
    fig.colorbar(sc, ax=axes, fraction=0.02, pad=0.02, label="ABSOLUTE purity")
    fig.suptitle("TCGA-LUSC  n=501  (colour = tumour purity)", y=1.04)
    fig.savefig(os.path.join(FIG, "fig3_scatter_lusc.png"))
    plt.close()


def fig_null_hist() -> None:
    tn = pd.read_csv(os.path.join(TAB, "transcriptome_wide_null.csv"))
    with open(os.path.join(PROC_DIR, "null_distributions.pkl"), "rb") as fh:
        dists = pickle.load(fh)
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.6), sharex=False)
    pairs = [
        ("TCGA-LUSC", "T_cell_CD8"), ("TCGA-LUSC", "B_cell"), ("TCGA-LUSC", "TLS_Cabrita"),
        ("TCGA-LUAD", "T_cell_CD8"), ("TCGA-LUAD", "B_cell"), ("TCGA-LUAD", "Plasma_cell"),
    ]
    for ax, (coh, sig) in zip(axes.ravel(), pairs):
        key = (coh, sig, "adjusted")
        if key not in dists:
            ax.set_visible(False)
            continue
        s = dists[key].dropna()
        ax.hist(s, bins=80, color="#cccccc", edgecolor="none")
        row = tn[(tn.cohort == coh) & (tn.signature == sig) & (tn.gene == "TACSTD2")]
        if len(row):
            v = float(row.iloc[0]["rho_adjusted"])
            pct = float(row.iloc[0]["adj_percentile"])
            ax.axvline(v, color="#b2182b", lw=1.6)
            ax.text(0.02, 0.95, f"TACSTD2 ρ={v:.2f}\npercentile {pct:.1f}",
                    transform=ax.transAxes, va="top", fontsize=7, color="#b2182b")
        ax.set_title(f"{coh}  {sig}\n(purity-adjusted, all genes)", fontsize=8)
        ax.set_xlabel("partial Spearman ρ")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_transcriptome_null.png"))
    plt.close()


def fig_gse72094() -> None:
    corr = pd.read_csv(os.path.join(TAB, "correlations_by_cohort.csv"))
    sub = corr[(corr.cohort == "GSE72094") & (corr.gene == "TACSTD2")
               & (corr.signature.isin(SIG_ORDER))].set_index("signature").loc[SIG_ORDER]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    y = np.arange(len(SIG_ORDER))
    ax.errorbar(sub["rho"], y - 0.15,
                xerr=[sub["rho"] - sub["ci_low"], sub["ci_high"] - sub["rho"]],
                fmt="o", color="#888888", ms=4, label="crude", capsize=2)
    ax.scatter(sub["rho_adj_epithelial"], y + 0.15, color="#b2182b", s=28,
               zorder=3, label="adj. epithelial score")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(SIG_ORDER)
    ax.set_xlabel("Spearman ρ")
    ax.set_title("GSE72094 LUAD microarray  n=442\nTACSTD2 vs immune signatures")
    ax.legend(frameon=False)
    for i, sig in enumerate(SIG_ORDER):
        ax.text(0.02, i + 0.15, fmt_p(sub.loc[sig, "p_adj_epithelial"]),
                va="center", fontsize=7, color="#b2182b", transform=ax.get_yaxis_transform())
    fig.savefig(os.path.join(FIG, "fig5_forest_gse72094.png"))
    plt.close()


def fig_verdict_table() -> None:
    """Render a compact text figure of the hold/null/opposite calls."""
    dec = pd.read_csv(os.path.join(TAB, "where_it_holds.csv"))
    sigs = ["TLS_Cabrita", "B_cell", "Plasma_cell", "T_cell_CD8"]
    cohs = [c for c in COH_ORDER]
    sub = dec[(dec.gene == "TACSTD2") & (dec.signature.isin(sigs))]
    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    ax.set_xlim(-0.5, len(cohs) - 0.5)
    ax.set_ylim(-0.5, len(sigs) - 0.5)
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels(cohs, rotation=30, ha="right")
    ax.set_yticks(range(len(sigs)))
    ax.set_yticklabels(sigs)
    color = {"HOLDS": "#2166ac", "NO_EVIDENCE": "#f0f0f0", "OPPOSITE": "#b2182b",
             "UNDERPOWERED": "#dddddd"}
    for i, sig in enumerate(sigs):
        for j, coh in enumerate(cohs):
            row = sub[(sub.signature == sig) & (sub.cohort == coh)]
            if row.empty:
                continue
            v = row.iloc[0]["verdict"]
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                       facecolor=color[v], edgecolor="white", lw=1.2))
            txt = f"{row.iloc[0]['rho_adj']:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if v in ("HOLDS", "OPPOSITE") else "black")
    ax.set_title("TACSTD2 after purity correction   blue=HOLDS (ρ<0 p<0.05 n≥40)  "
                 "red=OPPOSITE  grey=no evidence / n<40")
    ax.invert_yaxis()
    fig.savefig(os.path.join(FIG, "fig6_verdict_grid.png"))
    plt.close()


def main() -> None:
    _style()
    fig_heatmap_adj()
    fig_forest_lusc_vs_luad()
    fig_scatter_lusc()
    fig_null_hist()
    fig_gse72094()
    fig_verdict_table()
    print("[done] figures in", FIG)


if __name__ == "__main__":
    main()
