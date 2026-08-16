#!/usr/bin/env python3
"""Step 6 - figures for the LUSC / ICI TACSTD2-CLDN4 slice."""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import FIGURES_DIR, LOGS_DIR, TABLES_DIR

plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 140,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})


def save(fig, name: str) -> None:
    path = FIGURES_DIR / name
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path.name}")


def fig_histology():
    clin = pd.read_csv(TABLES_DIR / "cohort_clinical_with_histology.csv")
    clin["histology"] = clin["histology"].fillna("NA")
    path = clin[clin["histology"].isin(["LUSC", "non-LUSC"])].copy()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    cohorts = list(path["cohort"].unique())
    rng = np.random.default_rng(0)
    for i, cohort in enumerate(cohorts):
        sub = path[path["cohort"] == cohort]
        for label, color in (("LUSC", "#b85c38"), ("non-LUSC", "#2c6e8a")):
            y = sub.loc[sub["histology"] == label, "histology_delta"].to_numpy()
            if not len(y):
                continue
            x = i + rng.uniform(-0.12, 0.12, size=len(y))
            ax.scatter(x, y, s=18, alpha=0.75, c=color, edgecolors="none",
                       label=label if i == 0 else None)
    ax.axhline(0.20, color="#b85c38", ls="--", lw=0.8, label="inferred LUSC gate")
    ax.axhline(-0.20, color="#2c6e8a", ls="--", lw=0.8, label="inferred non-LUSC gate")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts, rotation=20, ha="right")
    ax.set_ylabel("squamous − adenocarcinoma marker score")
    ax.set_title("Pathology-annotated ICI samples (classifier is not used for these labels)")
    ax.legend(loc="lower left", fontsize=8)
    save(fig, "fig01_histology_scores.png")


def fig_inventory():
    inv = pd.read_csv(TABLES_DIR / "analysis_set_inventory.csv")
    primary = inv[inv["analysis_set"] == "primary_pathology_LUSC"].copy()
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    x = np.arange(len(primary))
    ax.bar(x - 0.18, primary["n_benefit"], 0.36, color="#3d8b6e", label="benefit")
    ax.bar(x + 0.18, primary["n_nobenefit"], 0.36, color="#a65d4e", label="no benefit")
    ax.set_xticks(x)
    ax.set_xticklabels(primary["cohort"])
    ax.set_ylabel("pre-treatment LUSC samples")
    ax.set_title("Primary set: pathology-confirmed LUSC with a binary ICI endpoint")
    ax.legend()
    save(fig, "fig02_primary_inventory.png")


def fig_response_boxes():
    samples = pd.read_csv(TABLES_DIR / "sample_level_scores.csv")
    samples["histology_final"] = samples["histology_final"].fillna("NA")
    use = samples[(samples["timepoint"] == "pre-treatment")
                  & (samples["histology_tier"] == "pathology")
                  & (samples["histology_final"] == "LUSC")
                  & samples["benefit"].notna()].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), sharey=False)
    for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
        parts = []
        labels = []
        for cohort, sub in use.groupby("cohort"):
            if sub[gene].notna().sum() < 4:
                continue
            for ben, tag in ((1, "B"), (0, "NB")):
                vals = sub.loc[sub["benefit"] == ben, gene].dropna().to_numpy()
                if len(vals) < 2:
                    continue
                parts.append(vals)
                labels.append(f"{cohort}\n{tag} n={len(vals)}")
        if not parts:
            ax.set_title(f"{gene}: no eligible cohort")
            continue
        bp = ax.boxplot(parts, tick_labels=labels, patch_artist=True, widths=0.6)
        for i, box in enumerate(bp["boxes"]):
            is_benefit = "\nB " in labels[i]
            box.set_facecolor("#3d8b6e" if is_benefit else "#a65d4e")
            box.set_alpha(0.7)
        ax.set_title(gene)
        ax.set_ylabel("log2 expression (cohort native scale)")
        ax.tick_params(axis="x", labelsize=7)
    fig.suptitle("Pathology LUSC: expression by ICI benefit (B) vs no benefit (NB)", y=1.02)
    save(fig, "fig03_response_boxplots.png")


def fig_forest():
    resp = pd.read_csv(TABLES_DIR / "response_by_cohort.csv")
    meta = pd.read_csv(TABLES_DIR / "response_meta.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharex=True)
    for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
        rows = resp[(resp["analysis_set"] == "primary_pathology_LUSC") & (resp["gene"] == gene)]
        m = meta[(meta["analysis_set"] == "primary_pathology_LUSC") & (meta["gene"] == gene)]
        names, est, lo, hi = [], [], [], []
        for _, r in rows.iterrows():
            names.append(f"{r['cohort']} ({int(r['n_benefit'])}/{int(r['n_nobenefit'])})")
            est.append(r["hedges_g"])
            lo.append(r["hedges_g"] - 1.96 * r["hedges_g_se"])
            hi.append(r["hedges_g"] + 1.96 * r["hedges_g_se"])
        if len(m):
            names.append(f"RE meta (k={int(m.iloc[0]['k'])})")
            est.append(m.iloc[0]["estimate"])
            lo.append(m.iloc[0]["ci_low"])
            hi.append(m.iloc[0]["ci_high"])
        y = np.arange(len(names))[::-1]
        ax.axvline(0, color="0.5", lw=0.8)
        for i, (e, a, b) in enumerate(zip(est, lo, hi)):
            color = "#222" if i < len(names) - 1 else "#b85c38"
            ax.plot([a, b], [y[i], y[i]], color=color, lw=1.4)
            ax.plot(e, y[i], "o", color=color, ms=5)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("Hedges' g (benefit − no benefit)")
        ax.set_title(gene)
    fig.suptitle("Primary pathology-LUSC response meta-analysis", y=1.02)
    save(fig, "fig04_response_forest.png")


def fig_immune():
    ici = pd.read_csv(TABLES_DIR / "immune_meta.csv")
    tcga = pd.read_csv(TABLES_DIR / "tcga_lusc_immune.csv")
    sigs = ["CYT", "IFNG_6", "TIS_18", "CD8_T_CELL", "CHECKPOINT",
            "ANTIGEN_PRESENTATION", "TGFB_FTBRS", "EPITHELIAL", "PROLIFERATION"]
    genes = ["TACSTD2", "CLDN4"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    for ax, (title, frame, rho_col) in zip(axes, (
        ("ICI pathology LUSC (RE meta ρ)",
         ici[ici["analysis_set"] == "primary_pathology_LUSC"], "rho"),
        ("TCGA-LUSC (Spearman ρ)", tcga, "rho"),
    )):
        mat = np.full((len(genes), len(sigs)), np.nan)
        for i, gene in enumerate(genes):
            for j, sig in enumerate(sigs):
                hit = frame[(frame["gene"] == gene) & (frame["signature"] == sig)]
                if len(hit) and pd.notna(hit.iloc[0][rho_col]):
                    mat[i, j] = hit.iloc[0][rho_col]
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_xticks(range(len(sigs)))
        ax.set_xticklabels(sigs, rotation=40, ha="right", fontsize=7)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.set_title(title)
        for i in range(len(genes)):
            for j in range(len(sigs)):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("TACSTD2 / CLDN4 versus immune and epithelial programmes", y=1.02)
    save(fig, "fig05_immune_correlations.png")


def fig_loo():
    loo = pd.read_csv(TABLES_DIR / "response_leave_one_cohort_out.csv")
    use = loo[loo["analysis_set"] == "primary_pathology_LUSC"]
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    for gene, color in (("TACSTD2", "#b85c38"), ("CLDN4", "#2c6e8a")):
        sub = use[use["gene"] == gene]
        if not len(sub):
            continue
        y = np.arange(len(sub))
        offset = -0.15 if gene == "TACSTD2" else 0.15
        ax.errorbar(sub["estimate"], y + offset,
                    xerr=[sub["estimate"] - sub["ci_low"], sub["ci_high"] - sub["estimate"]],
                    fmt="o", color=color, label=gene, capsize=2)
        ax.set_yticks(y)
        ax.set_yticklabels([f"drop {d}" for d in sub["dropped"]])
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("RE meta Hedges' g after dropping one cohort")
    ax.set_title("Leave-one-cohort-out (primary pathology LUSC)")
    ax.legend()
    save(fig, "fig06_leave_one_cohort_out.png")


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_histology()
    fig_inventory()
    fig_response_boxes()
    fig_forest()
    fig_immune()
    fig_loo()
    (LOGS_DIR / "06_figures.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
