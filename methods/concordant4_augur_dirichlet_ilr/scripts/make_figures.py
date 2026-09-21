#!/usr/bin/env python3
"""Figures for the concordant-4 compositional and Augur results."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

COHORT_COLOR = {
    "GSE123902": "#0072B2",
    "GSE131907": "#E69F00",
    "GSE205335": "#009E73",
    "GSE189357": "#CC79A7",
}


def _style():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
            "savefig.dpi": 200,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _rho_label(fraction: str) -> str:
    sp = pd.read_csv(TAB / "fraction_spearman.tsv", sep="\t")
    hit = sp[(sp.fraction == fraction) & (sp.dataset == "DL_meta")].iloc[0]
    return f"DL ρ = {hit.rho:+.2f}"


def fig_fractions(tab: pd.DataFrame) -> None:
    specs = [
        ("frac_T", "T fraction", _rho_label("T")),
        ("frac_NK", "NK fraction", _rho_label("NK")),
        ("frac_myeloid", "Myeloid fraction", _rho_label("myeloid")),
        ("frac_malignant", "Malignant fraction", _rho_label("malignant")),
    ]
    tab = tab.copy()
    tab["frac_malignant"] = tab.n_malignant / tab.n_cells
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), sharex=True)
    for ax, (col, title, rho) in zip(axes.ravel(), specs):
        for ds, sub in tab.groupby("dataset"):
            ax.scatter(
                sub.mal_CLDN4_pct,
                sub[col],
                s=28,
                c=COHORT_COLOR[ds],
                label=ds,
                alpha=0.9,
                linewidths=0,
            )
        ax.set_title(f"{title}  ({rho})")
        ax.set_ylabel("fraction of cells")
        ax.set_xlabel("Malignant CLDN4 % positive")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Concordant-4 patient units (n = 65)", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "fig_fractions_vs_cldn4.png", bbox_inches="tight")
    fig.savefig(FIG / "fig_fractions_vs_cldn4.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_ilr(ols: pd.DataFrame) -> None:
    want = [
        ("T_NK_myeloid_Rest", "immune_focus_vs_rest", "T/NK/myeloid vs rest"),
        ("T_NK_myeloid_Rest", "lymphoid_vs_myeloid", "T/NK vs myeloid"),
        ("T_NK_myeloid_Rest", "T_vs_NK", "T vs NK"),
    ]
    rows = []
    for simplex, balance, label in want:
        hit = ols[(ols.simplex == simplex) & (ols.balance == balance)].iloc[0]
        rows.append((label, hit.beta_per_sd, hit.se_hc1, hit.p_perm))
    labels = [r[0] for r in rows][::-1]
    beta = np.array([r[1] for r in rows][::-1])
    se = np.array([r[2] for r in rows][::-1])
    p = [r[3] for r in rows][::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    y = np.arange(len(labels))
    ax.axvline(0, color="#666666", lw=0.8)
    ax.errorbar(beta, y, xerr=1.96 * se, fmt="o", color="#0072B2", ecolor="#0072B2", capsize=3, ms=7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("ILR change per within-cohort SD of malignant CLDN4 %")
    for yi, (b, pv) in enumerate(zip(beta, p)):
        ax.text(b, yi + 0.18, f"perm p = {pv:.3f}", ha="center", va="bottom", fontsize=8, color="#333333")
    ax.set_title("Pre-specified balances, cohort-adjusted (n = 65)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_ilr_forest.png", bbox_inches="tight")
    fig.savefig(FIG / "fig_ilr_forest.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_dirichlet(dm: pd.DataFrame) -> None:
    sub = dm[dm.simplex == "T_NK_myeloid_Rest"].copy()
    order = ["T", "NK", "myeloid", "Rest"]
    sub["part"] = pd.Categorical(sub.part, order, ordered=True)
    sub = sub.sort_values("part")
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    colors = ["#0072B2", "#56B4E9", "#E69F00", "#999999"]
    x = np.arange(len(sub))
    heights = sub["delta_pi_plus_minus_0.5sd"].to_numpy()
    ax.bar(x, heights, color=colors, width=0.72)
    ax.axhline(0, color="#333333", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["T", "NK", "Myeloid", "Rest"])
    for xi, (h, pv) in enumerate(zip(heights, sub.p_perm)):
        if not np.isfinite(pv):
            continue
        label = f"p = {pv:.3f}"
        y = h + (0.006 if h >= 0 else -0.006)
        va = "bottom" if h >= 0 else "top"
        ax.text(xi, y, label, ha="center", va=va, fontsize=8, color="#333333")
    ax.set_ylabel("Predicted fraction change\n(+0.5 SD minus −0.5 SD CLDN4)")
    ax.set_title("Dirichlet-multinomial mean (n = 65)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_dirichlet_delta.png", bbox_inches="tight")
    fig.savefig(FIG / "fig_dirichlet_delta.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_augur(path: Path) -> None:
    if not path.exists():
        return
    tab = pd.read_csv(path, sep="\t")
    tab = tab.dropna(subset=["unit_auc"])
    if tab.empty:
        return
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    y = np.arange(len(tab))[::-1]
    ax.axvline(0.5, color="#666666", lw=0.8)
    colors = ["#009E73" if r.cell_type == "malignant" else "#0072B2" for r in tab.itertuples()]
    ax.barh(y, tab.unit_auc, color=colors, height=0.65, xerr=tab.unit_auc_sd, capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.cell_type}  (p={r.p_perm:.2f})" for r in tab.itertuples()]
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("AUC of unit-mean out-of-fold probability")
    ax.set_title("Augur-style priority, patient-blocked CV (n = 35 Q1/Q4 units)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_augur_priority.png", bbox_inches="tight")
    fig.savefig(FIG / "fig_augur_priority.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    _style()
    tab = pd.read_csv(TAB / "composition_counts.tsv", sep="\t")
    fig_fractions(tab)
    fig_ilr(pd.read_csv(TAB / "ilr_ols.tsv", sep="\t"))
    fig_dirichlet(pd.read_csv(TAB / "dirichlet_multinomial.tsv", sep="\t"))
    fig_augur(TAB / "augur_priority.tsv")
    print("figures written", flush=True)


if __name__ == "__main__":
    main()
