"""Figures for the Chan atlas analysis."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from . import chan_atlas as ca
from . import bulk_george as bg

FIGS = bg.FIGS
SUBTYPE_COLORS = {"SCLC-A": "#4C72B0", "SCLC-N": "#DD8452",
                  "SCLC-P": "#55A868", "SCLC": "#937860"}
COMP_COLORS = {"Epithelial": "#4C72B0", "Lymphoid": "#C44E52",
               "Myeloid": "#DD8452", "Mesenchymal": "#55A868"}


def fig_detect_by_compartment(comp):
    genes = ["TACSTD2", "CLDN4"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    comps = list(comp["compartment"])
    x = np.arange(len(comps))
    for ax, g in zip(axes, genes):
        y = comp[f"{g}_detect_frac"].to_numpy() * 100
        colors = [COMP_COLORS[c] for c in comps]
        ax.bar(x, y, color=colors, edgecolor="0.2", width=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(comps, rotation=20, ha="right")
        ax.set_ylabel("% cells detected (>0)")
        ax.set_title(f"{g} detection in SCLC samples")
        ax.set_ylim(0, max(15, y.max() * 1.25) if np.isfinite(y).any() else 15)
        for i, (yi, n) in enumerate(zip(y, comp["n_cells"])):
            ax.text(i, yi + 0.3, f"n={n:,}", ha="center", va="bottom", fontsize=7)
    fig.suptitle("Chan atlas — ADC-target detection by cell compartment (SCLC samples only)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIGS / "chan_detect_by_compartment.png", dpi=150)
    plt.close(fig)


def fig_detect_by_tumor_label(subt):
    genes = ["TACSTD2", "CLDN4"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    labs = list(subt["chan_tumor_label"])
    x = np.arange(len(labs))
    for ax, g in zip(axes, genes):
        y = subt[f"{g}_detect_frac"].to_numpy() * 100
        colors = [SUBTYPE_COLORS.get(l, "0.5") for l in labs]
        ax.bar(x, y, color=colors, edgecolor="0.2", width=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{l}\n(n={n:,})" for l, n in zip(labs, subt["n_cells"])],
                           fontsize=8)
        ax.set_ylabel("% tumor cells detected (>0)")
        ax.set_title(f"{g} in Chan tumor-cell labels")
        ymax = max(10, np.nanmax(y) * 1.3) if np.isfinite(y).any() else 10
        ax.set_ylim(0, ymax)
    fig.suptitle("Chan atlas — no SCLC-I tumor-cell label exists; A/N/P only",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIGS / "chan_detect_by_tumor_label.png", dpi=150)
    plt.close(fig)


def fig_patient_scatter(pat):
    use = pat[pat["keep_for_stats"]].copy()
    pairs = [("tumor_TACSTD2_mean", "frac_immune", "TACSTD2 (tumor-cell mean)"),
             ("tumor_CLDN4_mean", "frac_immune", "CLDN4 (tumor-cell mean)"),
             ("tumor_TACSTD2_mean", "frac_tcell", "TACSTD2 (tumor-cell mean)"),
             ("tumor_CLDN4_mean", "frac_tcell", "CLDN4 (tumor-cell mean)")]
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 8.0))
    for ax, (gcol, axcol, ylab) in zip(axes.ravel(), pairs):
        x = use[axcol].astype(float)
        y = use[gcol].astype(float)
        c = use["dominant_tumor_label"].map(SUBTYPE_COLORS).fillna("0.5")
        ax.scatter(x, y, c=c, s=36, alpha=0.85, edgecolor="0.25", linewidth=0.4)
        m = x.notna() & y.notna()
        if m.sum() >= 6:
            rho, p = stats.spearmanr(x[m], y[m])
            b, a = np.polyfit(x[m], y[m], 1)
            xs = np.linspace(x[m].min(), x[m].max(), 40)
            ax.plot(xs, a + b * xs, color="0.3", lw=1, ls="--")
            ax.set_title(f"rho={rho:.2f}, p={p:.3g}, n={int(m.sum())}", fontsize=9)
        ax.set_xlabel(axcol.replace("_", " "))
        ax.set_ylabel(ylab)
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=SUBTYPE_COLORS[s],
                          label=s) for s in ["SCLC-A", "SCLC-N", "SCLC-P"]]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Chan atlas — patient-level tumor ADC-target vs immune fraction\n"
                 "(SCLC donors with ≥30 tumor cells; cells are not the unit of inference)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(FIGS / "chan_patient_ADC_vs_immune.png", dpi=150)
    plt.close(fig)


def run():
    res = ca.run()
    fig_detect_by_compartment(res["comp"])
    fig_detect_by_tumor_label(res["subt"])
    fig_patient_scatter(res["pat"])
    print("Chan figures written to", FIGS)


if __name__ == "__main__":
    run()
