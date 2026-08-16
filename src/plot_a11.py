"""Figures for A11. All numbers are read from the tables; nothing is recomputed."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C
from stats_utils import spearman, partial_spearman

TAB = C.OUT / "tables"
FIG = C.OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 160,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

AXIS_COLOR = {
    "CD47": "#c0392b",
    "Galectin": "#8e44ad",
    "Nectin": "#2980b9",
    "TGFB": "#16a085",
    "negative_control": "#7f8c8d",
}


def _load_primary():
    g = pd.read_csv(TAB / "10_bulk_gene_level.csv")
    return g[~g.is_control].copy(), g[g.is_control].copy()


def fig01_cd47_scatter():
    """The money plot: CD47 vs TACSTD2, raw and purity-residual, both histologies."""
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2))
    for i, cohort in enumerate(["LUAD", "LUSC"]):
        mat = pd.read_parquet(C.DATA / "derived" / f"bulk_{cohort}_panel.parquet")
        samp = pd.read_parquet(C.DATA / "derived" / f"bulk_{cohort}_samples.parquet")
        keep = samp.sample_type == "PrimaryTumor"
        mat, samp = mat.loc[keep], samp.loc[keep]
        ok = samp["purity"].notna()
        x = mat.loc[ok, "TACSTD2"].values
        y = mat.loc[ok, "CD47"].values
        p = samp.loc[ok, "purity"].values

        r = spearman(x, y)
        ax = axes[i, 0]
        ax.scatter(x, y, s=10, c="#2c3e50", alpha=0.35, linewidths=0)
        ax.set_xlabel("TACSTD2  log2(TPM+1)")
        ax.set_ylabel("CD47  log2(TPM+1)")
        ax.set_title(f"{cohort}  n={r['n']}   ρ={r['rho']:.2f} "
                     f"[{r['lo']:.2f}, {r['hi']:.2f}]")

        # Rank-residuals after purity (same transform as the partial Spearman).
        from scipy.stats import rankdata
        rx, ry, rp = rankdata(x), rankdata(y), rankdata(p)
        design = np.column_stack([np.ones(len(x)), rp])
        ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
        ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
        adj = partial_spearman(x, y, [p])
        ax = axes[i, 1]
        ax.scatter(ex, ey, s=10, c="#c0392b", alpha=0.35, linewidths=0)
        ax.axhline(0, color="0.7", lw=0.6)
        ax.axvline(0, color="0.7", lw=0.6)
        ax.set_xlabel("TACSTD2 rank residual | purity")
        ax.set_ylabel("CD47 rank residual | purity")
        ax.set_title(f"{cohort}  partial ρ={adj['rho']:.2f} "
                     f"[{adj['lo']:.2f}, {adj['hi']:.2f}]")

    fig.suptitle("CD47 vs TACSTD2 in TCGA primary tumours, before and after ABSOLUTE purity",
                 y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "01_cd47_vs_tacstd2_scatter.png", bbox_inches="tight")
    fig.savefig(FIG / "01_cd47_vs_tacstd2_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def fig02_forest():
    """Naive vs purity-adjusted rho for the primary panel, LUAD and LUSC."""
    panel, ctrl = _load_primary()
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 8.6), sharex=True)
    for ax, cohort in zip(axes, ["LUAD", "LUSC"]):
        d = panel[panel.cohort == cohort].copy()
        d = d.sort_values(["axis", "rho_adj_purity"])
        y = np.arange(len(d))
        colors = d.axis.map(AXIS_COLOR)
        ax.axvline(0, color="0.5", lw=0.7)
        ax.axvline(0.3, color="0.75", ls="--", lw=0.7)
        ax.axvline(-0.3, color="0.75", ls="--", lw=0.7)
        # naive as thin grey, adjusted as coloured
        ax.hlines(y, d.rho_naive_lo, d.rho_naive_hi, color="0.75", lw=1.2)
        ax.plot(d.rho_naive, y, "o", color="0.55", ms=3.5, zorder=3)
        ax.hlines(y, d.rho_adj_lo, d.rho_adj_hi, color=colors, lw=2.0)
        ax.scatter(d.rho_adj_purity, y, c=colors, s=18, zorder=4, edgecolors="none")
        labels = [f"{g}  {a}" for g, a in zip(d.gene, d.axis)]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Spearman ρ with TACSTD2")
        ax.set_title(cohort)
        ax.set_xlim(-0.7, 0.7)
    handles = [Line2D([0], [0], marker="o", color="0.55", ls="", label="naive"),
               Line2D([0], [0], marker="o", color="#c0392b", ls="", label="purity-adjusted")]
    axes[1].legend(handles=handles, loc="lower right")
    fig.suptitle("Primary panel vs TACSTD2. Grey = naive; colour = partial ρ | ABSOLUTE purity. "
                 "Dashed lines = pre-specified |ρ|=0.3", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "02_forest_naive_vs_purity.png", bbox_inches="tight")
    fig.savefig(FIG / "02_forest_naive_vs_purity.pdf", bbox_inches="tight")
    plt.close(fig)


def fig03_attenuation_and_controls():
    """Does purity actually move anything? And where do housekeepers sit?"""
    panel, ctrl = _load_primary()
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2))

    ax = axes[0]
    for cohort, mk in [("LUAD", "o"), ("LUSC", "s")]:
        d = panel[panel.cohort == cohort]
        ax.scatter(d.rho_naive, d.rho_adj_purity, c=d.axis.map(AXIS_COLOR),
                   marker=mk, s=28, alpha=0.85, label=cohort, edgecolors="none")
        # highlight CD47
        cd = d[d.gene == "CD47"]
        ax.scatter(cd.rho_naive, cd.rho_adj_purity, facecolors="none",
                   edgecolors="k", s=90, lw=1.2, zorder=5)
    ax.axline((0, 0), slope=1, color="0.6", lw=0.7)
    ax.axhline(0, color="0.8", lw=0.5)
    ax.axvline(0, color="0.8", lw=0.5)
    ax.set_xlabel("naive ρ")
    ax.set_ylabel("purity-adjusted ρ")
    ax.set_title("Purity barely moves any estimate\n(open circle = CD47)")
    ax.legend(loc="upper left")
    ax.set_xlim(-0.4, 0.6)
    ax.set_ylim(-0.4, 0.6)

    ax = axes[1]
    for cohort, mk in [("LUAD", "o"), ("LUSC", "s")]:
        c = ctrl[ctrl.cohort == cohort]
        ax.scatter(np.zeros(len(c)) + (0 if cohort == "LUAD" else 1),
                   c.rho_adj_purity, marker=mk, c="#7f8c8d", s=22, alpha=0.8)
        cd = panel[(panel.cohort == cohort) & (panel.gene == "CD47")]
        ax.scatter([0 if cohort == "LUAD" else 1], cd.rho_adj_purity,
                   marker=mk, c="#c0392b", s=50, zorder=4, label=f"CD47 {cohort}")
    ax.axhline(0, color="0.6", lw=0.7)
    ax.axhline(0.3, color="0.75", ls="--", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LUAD", "LUSC"])
    ax.set_ylabel("purity-adjusted ρ with TACSTD2")
    ax.set_title("CD47 (red) vs 10 housekeeping genes (grey)")
    ax.set_xlim(-0.5, 1.5)
    fig.tight_layout()
    fig.savefig(FIG / "03_purity_shift_and_housekeepers.png", bbox_inches="tight")
    fig.savefig(FIG / "03_purity_shift_and_housekeepers.pdf", bbox_inches="tight")
    plt.close(fig)


def fig04_purity_strata():
    s = pd.read_csv(TAB / "12_bulk_purity_strata.csv")
    s = s[s.feature == "module_CD47"].copy()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    order = ["low_purity", "mid_purity", "high_purity"]
    x = np.arange(len(order))
    width = 0.35
    for i, cohort in enumerate(["LUAD", "LUSC"]):
        d = s[s.cohort == cohort].set_index("stratum").loc[order]
        xx = x + (i - 0.5) * width
        ax.errorbar(xx, d.rho, yerr=[d.rho - d.lo, d.hi - d.rho],
                    fmt="o", ms=6, capsize=3, label=cohort,
                    color="#c0392b" if cohort == "LUAD" else "#2980b9")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.axhline(0.3, color="0.75", ls="--", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(["low purity", "mid", "high purity"])
    ax.set_ylabel("Spearman ρ  CD47 vs TACSTD2")
    ax.set_title("CD47–TACSTD2 inside ABSOLUTE-purity tertiles (pre-specified A3)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "04_cd47_purity_strata.png", bbox_inches="tight")
    fig.savefig(FIG / "04_cd47_purity_strata.pdf", bbox_inches="tight")
    plt.close(fig)


def fig05_tumor_normal():
    tn = pd.read_csv(TAB / "13_bulk_tumor_vs_normal.csv")
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    genes = ["TACSTD2", "CD47", "EPCAM"]
    x = np.arange(len(genes))
    width = 0.35
    for i, cohort in enumerate(["LUAD", "LUSC"]):
        d = tn[(tn.cohort == cohort) & (tn.gene.isin(genes))].set_index("gene").loc[genes]
        ax.bar(x + (i - 0.5) * width, d.log2fc_tumor_vs_normal, width=width,
               label=cohort, color="#c0392b" if cohort == "LUAD" else "#2980b9",
               alpha=0.85)
    ax.axhline(0, color="0.4", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.set_ylabel("mean log2(TPM+1)  tumour − adjacent normal")
    ax.set_title("Context: CD47 is lower in tumour than adjacent lung, in both histologies")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "05_tumor_vs_normal.png", bbox_inches="tight")
    fig.savefig(FIG / "05_tumor_vs_normal.pdf", bbox_inches="tight")
    plt.close(fig)


def fig06_sc_pseudobulk():
    """Patient-level malignant-cell CD47 vs TACSTD2, coloured by site.

    The pooled cloud looks tight; colouring by site is the honest view.
    """
    pb = pd.read_csv(TAB / "21_sc_malignant_pseudobulk_samples.csv")
    # Use the mixed-no-brain slice plus brain, plotted together, one point
    # per sample (samples can appear in only one slice).
    keep = pb.slice.isin(["tumor_sites_mixed_no_brain", "brain_mets",
                          "primary_tLung_tS", "mets_mLN_tLB"])
    # Prefer the site-specific rows so each patient appears once.
    one = pb[pb.slice.isin(["primary_tLung_tS", "mets_mLN_tLB", "brain_mets"])]
    one = one.drop_duplicates("Sample")

    color = {"primary_tLung_tS": "#c0392b",
             "mets_mLN_tLB": "#2980b9",
             "brain_mets": "#7f8c8d"}
    label = {"primary_tLung_tS": "primary tLung (tS1/2/3)",
             "mets_mLN_tLB": "LN / bronchus mets",
             "brain_mets": "brain mets"}

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    for sl, g in one.groupby("slice"):
        ax.scatter(g.TACSTD2, g.CD47, s=np.clip(g.n_cells / 8, 18, 90),
                   c=color[sl], label=label[sl], alpha=0.85, edgecolors="k",
                   linewidths=0.4)
    ax.set_xlabel("TACSTD2  patient mean in malignant cells")
    ax.set_ylabel("CD47  patient mean in malignant cells")
    ax.set_title("GSE131907  malignant-cell pseudobulk\n"
                 "point size ~ n cells; site is not a nuisance, it is structure")
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "06_sc_malignant_pseudobulk.png", bbox_inches="tight")
    fig.savefig(FIG / "06_sc_malignant_pseudobulk.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    fig01_cd47_scatter()
    fig02_forest()
    fig03_attenuation_and_controls()
    fig04_purity_strata()
    fig05_tumor_normal()
    fig06_sc_pseudobulk()
    print("wrote figures to", FIG)


if __name__ == "__main__":
    main()
