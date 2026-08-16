"""Figures for results/claim_B6/."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

OUT = "results/claim_B6/figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _strip(ax, df, xcol, ycol, order, colors):
    rng = np.random.default_rng(0)
    for i, lab in enumerate(order):
        y = df.loc[df[xcol] == lab, ycol].dropna().values
        x = np.full(len(y), i) + rng.uniform(-0.12, 0.12, len(y))
        ax.scatter(x, y, s=22, c=colors[i], alpha=0.85, edgecolors="none", zorder=3)
        if len(y):
            ax.hlines(np.median(y), i - 0.22, i + 0.22, colors="black", lw=1.6, zorder=4)
    ax.axhline(0, color="0.6", lw=0.8, ls="--")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=20, ha="right")


def fig_samespot():
    em = pd.read_csv("results/claim_B6/tables/emtab13530_per_section.csv")
    g = pd.read_csv("results/claim_B6/tables/gse189487_per_section.csv")
    em["cohort"] = em["group"].map({
        "tumor": "E-MTAB tumor",
        "adjacent_normal": "E-MTAB adjacent",
        "healthy": "E-MTAB healthy",
    })
    g["cohort"] = "GSE189487 LUAD"
    df = pd.concat([em, g], ignore_index=True)
    order = ["E-MTAB tumor", "GSE189487 LUAD", "E-MTAB adjacent", "E-MTAB healthy"]
    colors = ["#b2182b", "#d6604d", "#4393c3", "#92c5de"]
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), sharey=True)
    for ax, col, title in zip(
        axes,
        ["rho_epi_vs_TB", "rho_epi_vs_T", "rho_epi_vs_B"],
        ["CLDN4/TACSTD2 vs T+B", "vs T-cell score", "vs B-cell score"],
    ):
        _strip(ax, df, "cohort", col, order, colors)
        ax.set_title(title)
        ax.set_ylabel("Spearman ρ (same spot)" if ax is axes[0] else "")
    fig.suptitle("Same-spot anti-colocalization (Visium)", y=1.03)
    fig.savefig(f"{OUT}/visium_samespot_rho.png")
    fig.savefig(f"{OUT}/visium_samespot_rho.pdf")
    plt.close()


def fig_neighborhood():
    nb = pd.read_csv("results/claim_B6/tables/visium_neighborhood_rings.csv")
    nb["cohort"] = np.where(
        nb["dataset"] == "GSE189487", "GSE189487 LUAD",
        nb["group"].map({
            "tumor": "E-MTAB tumor",
            "adjacent_normal": "E-MTAB adjacent",
            "healthy": "E-MTAB healthy",
        }),
    )
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), sharey=True)
    colors = {1: "#b2182b", 2: "#ef8a62", 3: "#fddbc7"}
    for ax, cohort in zip(axes, ["E-MTAB tumor", "GSE189487 LUAD", "E-MTAB adjacent"]):
        sub = nb[nb["cohort"] == cohort]
        for k in (1, 2, 3):
            y = sub.loc[sub["ring"] == k, "rho_CLDN4_vs_TBnb"].dropna().values
            rng = np.random.default_rng(k)
            x = np.full(len(y), k) + rng.uniform(-0.12, 0.12, len(y))
            ax.scatter(x, y, s=20, c=colors[k], alpha=0.9, edgecolors="none")
            if len(y):
                ax.hlines(np.median(y), k - 0.2, k + 0.2, colors="black", lw=1.5)
        ax.axhline(0, color="0.6", lw=0.8, ls="--")
        ax.set_xticks([1, 2, 3])
        ax.set_xlabel("hex ring")
        ax.set_title(cohort)
        if ax is axes[0]:
            ax.set_ylabel("Spearman ρ  CLDN4 vs neighbor T+B")
    fig.suptitle("Visium hex-ring neighborhood (raw, all spots)", y=1.03)
    fig.savefig(f"{OUT}/visium_neighborhood_rings.png")
    fig.savefig(f"{OUT}/visium_neighborhood_rings.pdf")
    plt.close()


def fig_pr110like():
    df = pd.read_csv("results/claim_B6/tables/visium_pr110like_rings.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=True)
    for ax, ds, title in zip(
        axes,
        ["E-MTAB-13530", "GSE189487"],
        ["E-MTAB tumor (n=20)", "GSE189487 (n=6)"],
    ):
        sub = df[df["dataset"] == ds]
        for k in (1, 2, 3):
            y = sub.loc[sub["ring"] == k, "rho_partialEpi"].dropna().values
            rng = np.random.default_rng(k)
            x = np.full(len(y), k) + rng.uniform(-0.12, 0.12, len(y))
            ax.scatter(x, y, s=22, c="#b2182b", alpha=0.85, edgecolors="none")
            if len(y):
                ax.hlines(np.median(y), k - 0.2, k + 0.2, colors="black", lw=1.5)
        ax.axhline(0, color="0.6", lw=0.8, ls="--")
        ax.set_xticks([1, 2, 3])
        ax.set_xlabel("hex ring")
        ax.set_title(title)
        if ax is axes[0]:
            ax.set_ylabel("partial ρ  CLDN4 vs neighbor T+B\n(index = epi-high ∩ not TB-rich)")
    fig.suptitle("PR#110-like filter, T/B neighborhood (this analysis — not PR #110)", y=1.04)
    fig.savefig(f"{OUT}/visium_pr110like_partial.png")
    fig.savefig(f"{OUT}/visium_pr110like_partial.pdf")
    plt.close()


def fig_geomx():
    aoi = pd.read_csv("results/claim_B6/tables/gse271689_geomx_per_aoi.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))
    order = ["CK", "CD68", "CD45"]
    cols = {"CK": "#b2182b", "CD68": "#f4a582", "CD45": "#2166ac"}
    for ax, col, title in zip(axes, ["epi_score", "tb_score"],
                              ["CLDN4/TACSTD2 score", "T+B score"]):
        data = [aoi.loc[aoi["cell_type"] == c, col].dropna().values for c in order]
        bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55,
                        medianprops=dict(color="black", lw=1.4),
                        flierprops=dict(marker=".", ms=3, alpha=0.4))
        for patch, c in zip(bp["boxes"], order):
            patch.set_facecolor(cols[c])
            patch.set_alpha(0.7)
        ax.set_title(title)
        ax.set_xlabel("GeoMx compartment (antibody segment)")
    fig.suptitle("GSE271689 GeoMx: scores by pre-defined compartment", y=1.03)
    fig.savefig(f"{OUT}/geomx_compartment_scores.png")
    fig.savefig(f"{OUT}/geomx_compartment_scores.pdf")
    plt.close()


def fig_cosmx():
    nb = pd.read_csv("results/claim_B6/tables/gse287472_cosmx_neighborhood.csv")
    same = pd.read_csv("results/claim_B6/tables/gse287472_cosmx_samecell.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    # same-cell as radius 0
    xs, ys, labs = [], [], []
    for i, r in same.iterrows():
        xs.append(-0.3 if r["scope"] == "all_QC_cells" else 0.3)
        ys.append(r["rho_CLDN4_vs_TB"])
        labs.append("same-cell\n" + r["scope"].replace("_", " "))
    ax.scatter([-0.15], [same.loc[same.scope == "all_QC_cells", "rho_CLDN4_vs_TB"].iloc[0]],
               s=50, c="#4d4d4d", label="same-cell, all QC", zorder=3)
    ax.scatter([0.15], [same.loc[same.scope == "PanCK_high", "rho_CLDN4_vs_TB"].iloc[0]],
               s=50, c="#b2182b", label="same-cell, PanCK-high", zorder=3)
    for r, c, lab in [(25, "#80cdc1", "all QC"), (50, "#35978f", None), (100, "#01665e", None)]:
        y = nb.loc[(nb.radius_um == r) & (nb.scope == "all_QC_cells"), "rho_CLDN4_vs_TBnb"].iloc[0]
        ax.scatter([r], [y], s=50, c="#35978f", zorder=3)
    for r in (25, 50, 100):
        y = nb.loc[(nb.radius_um == r) & (nb.scope == "PanCK_high"), "rho_CLDN4_vs_TBnb"].iloc[0]
        ax.scatter([r], [y], s=50, c="#b2182b", zorder=3)
    ax.axhline(0, color="0.6", lw=0.8, ls="--")
    ax.set_xticks([0, 25, 50, 100])
    ax.set_xticklabels(["same cell", "25 µm", "50 µm", "100 µm"])
    ax.set_ylabel("Spearman ρ  CLDN4 vs T+B (or neighbor T+B)")
    ax.set_title("GSE287472 CosMx LUAD (1 patient) — no neighborhood exclusion")
    ax.legend(frameon=False, loc="upper right")
    ax.set_ylim(-0.15, 0.15)
    fig.savefig(f"{OUT}/cosmx_neighborhood_rho.png")
    fig.savefig(f"{OUT}/cosmx_neighborhood_rho.pdf")
    plt.close()


if __name__ == "__main__":
    fig_samespot()
    fig_neighborhood()
    fig_pr110like()
    fig_geomx()
    fig_cosmx()
    print("wrote figures to", OUT)
