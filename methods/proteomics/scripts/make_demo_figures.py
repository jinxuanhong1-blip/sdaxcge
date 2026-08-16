#!/usr/bin/env python3
"""Publication-style demo figures for the CPTAC/TMT proteomics playbook."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib_io import (  # noqa: E402
    gene_presence,
    lookup_gene_rows,
    normalize_sample_id,
    read_expression_matrix,
    read_linkedomics_phenotype,
)

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
COLORS = {
    "LUAD": "#1f4e79",
    "LSCC": "#b85c38",
    "TACSTD2": "#2a6f97",
    "CLDN4": "#c1121f",
    "present": "#2a9d8f",
    "absent": "#c1121f",
    "Tumor": "#264653",
    "NAT": "#e9c46a",
}


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_missingness(protein: dict[str, pd.DataFrame], rna: dict[str, pd.DataFrame], out: Path) -> None:
    panel = ["TACSTD2", "CLDN1", "CLDN3", "CLDN4", "CLDN5", "CLDN7", "CLDN18", "EPCAM"]
    cohorts = list(protein)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), gridspec_kw={"width_ratios": [1.15, 1]})

    present = np.zeros((len(panel), len(cohorts) * 2))
    labels = []
    col = 0
    for cohort in cohorts:
        labels.append(f"{cohort}\nprotein")
        for i, g in enumerate(panel):
            present[i, col] = 1.0 if gene_presence(protein[cohort], g)["present"] else 0.0
        col += 1
        labels.append(f"{cohort}\nRNA")
        for i, g in enumerate(panel):
            present[i, col] = 1.0 if gene_presence(rna[cohort], g)["present"] else 0.0
        col += 1

    ax = axes[0]
    cmap = plt.matplotlib.colors.ListedColormap(["#f4d6d6", "#9ad0c2"])
    ax.imshow(present, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(len(panel)))
    ax.set_yticklabels(panel)
    ax.set_title("Row present in table (not NA rate)")
    for i, g in enumerate(panel):
        for j in range(present.shape[1]):
            ax.text(j, i, "yes" if present[i, j] else "no", ha="center", va="center", fontsize=7,
                    color="#1b1b1b")

    ax = axes[1]
    for cohort, mat in protein.items():
        na = mat.isna().mean(axis=1)
        ax.hist(na, bins=30, range=(0, 1), histtype="step", linewidth=1.6,
                label=f"{cohort} protein", color=COLORS[cohort])
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set_xlabel("Gene-wise NA fraction")
    ax.set_ylabel("Genes")
    ax.set_title("NArm already dropped high-NA genes")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("CLDN4 is absent as a row in TMT gene-abundance (RNA is not)", fontsize=12, y=1.02)
    _save(fig, out)


def fig_log2(protein: dict[str, pd.DataFrame], rna: dict[str, pd.DataFrame], out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=False)
    ax = axes[0]
    for cohort, mat in protein.items():
        v = pd.Series(mat.to_numpy().ravel()).dropna().sample(min(20000, mat.size), random_state=0)
        ax.hist(v, bins=50, density=True, histtype="step", lw=1.5, label=cohort, color=COLORS[cohort])
    ax.set_title("TMT protein (already log2-ratio)")
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")
    ax.legend(frameon=False)
    ax.axvline(0, color="0.4", lw=0.7, ls="--")

    ax = axes[1]
    for cohort, mat in rna.items():
        v = pd.Series(mat.to_numpy().ravel()).dropna().sample(min(20000, mat.size), random_state=0)
        ax.hist(v, bins=50, density=True, histtype="step", lw=1.5, label=cohort, color=COLORS[cohort])
    ax.set_title("RNA tables (already log2; scales differ)")
    ax.set_xlabel("Value")
    ax.legend(frameon=False)
    fig.suptitle("Do not double-log. LUAD vs LSCC RNA are not on one scale.", fontsize=12, y=1.02)
    _save(fig, out)


def fig_concordance(pairs: pd.DataFrame, stats_df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0))
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
        sub = pairs[pairs["gene"] == gene] if not pairs.empty else pairs
        if sub.empty:
            ax.set_title(f"{gene}: no paired protein–RNA")
            ax.text(0.5, 0.5, "protein row absent\n(CLDN4 TMT/NArm)", ha="center", va="center",
                    transform=ax.transAxes, color=COLORS["CLDN4"])
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        for cohort, g in sub.groupby("cohort"):
            ax.scatter(g["rna"], g["protein"], s=18, alpha=0.75, label=cohort, color=COLORS.get(cohort, "0.3"))
            if g.shape[0] >= 5:
                lr = stats.linregress(g["rna"], g["protein"])
                xs = np.linspace(g["rna"].min(), g["rna"].max(), 50)
                ax.plot(xs, lr.intercept + lr.slope * xs, color=COLORS.get(cohort, "0.3"), lw=1)
        st = stats_df[stats_df["gene"] == gene]
        bits = [f"{r.cohort} ρ={r.spearman_rho:.2f} n={int(r.n_paired)}" for r in st.itertuples() if pd.notna(r.spearman_rho)]
        ax.set_title(gene + (("\n" + "; ".join(bits)) if bits else ""))
        ax.set_xlabel("RNA (log2 table)")
        ax.set_ylabel("Protein TMT log2-ratio")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Protein–RNA concordance is only defined when both rows exist", fontsize=12, y=1.03)
    _save(fig, out)


def fig_limma(fit: pd.DataFrame, out: Path, highlight=("TACSTD2", "CLDN4", "EPCAM", "KRT7")) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    x = fit["logFC"]
    y = -np.log10(fit["P.Value"].clip(lower=1e-300))
    sig = fit["adj.P.Val"] < 0.05
    ax.scatter(x[~sig], y[~sig], s=8, c="#cfcfcf", linewidths=0)
    ax.scatter(x[sig], y[sig], s=10, c="#457b9d", linewidths=0, label="BH < 0.05")
    ax.axvline(0, color="0.5", lw=0.6)
    ax.axhline(-np.log10(0.05), color="0.5", lw=0.6, ls="--")
    for gene in highlight:
        hit = fit[fit["gene"].astype(str).str.upper() == gene.upper()]
        if hit.empty:
            continue
        r = hit.iloc[0]
        ax.scatter([r.logFC], [-np.log10(max(r["P.Value"], 1e-300))], s=40, c=COLORS.get(gene, "#e63939"), zorder=3)
        ax.annotate(gene, (r.logFC, -np.log10(max(r["P.Value"], 1e-300))),
                    textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax.set_xlabel("logFC Tumor − NAT (TMT log2-ratio)")
    ax.set_ylabel("−log10 P (limma eBayes)")
    ax.set_title("limma on already-logged TMT (not voom / not limma-trend-on-counts)")
    ax.legend(frameon=False, loc="upper left")
    _save(fig, out)


def fig_immune_join(joined: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    prot_col = next((c for c in joined.columns if c.startswith("protein_") and c.endswith("TACSTD2")), None)
    cluster_col = next((c for c in joined.columns if "Immune_Cluster" in c or "Immune.Cluster" in c), None)
    if cluster_col is None:
        cluster_col = next((c for c in joined.columns if "Cluster" in c and "Immune" in c), None)
    ax = axes[0]
    if prot_col and cluster_col:
        d = joined[[prot_col, cluster_col]].dropna()
        order = [x for x in ["Cold Tumor", "Warm Tumor", "Hot Tumor"] if x in set(d[cluster_col])]
        if not order:
            order = list(d[cluster_col].astype(str).value_counts().index)
        data = [d.loc[d[cluster_col] == k, prot_col].astype(float) for k in order]
        bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55)
        for patch, color in zip(bp["boxes"], ["#8d99ae", "#e9c46a", "#e76f51"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)
        ax.set_ylabel("TACSTD2 protein (TMT log2-ratio)")
        ax.set_title("LSCC Immune.Cluster.rna (treatment-naive)")
        if len(data) >= 2:
            try:
                stat, p = stats.kruskal(*[x for x in data if len(x)])
                ax.text(0.02, 0.98, f"Kruskal–Wallis p={p:.3g}", transform=ax.transAxes, va="top", fontsize=8)
            except Exception:
                pass
    else:
        ax.text(0.5, 0.5, "no Immune.Cluster column", ha="center")
    ax = axes[1]
    xcell = next((c for c in joined.columns if "xCell_ImmuneScore" in c), None)
    cib = next((c for c in joined.columns if "CIBERSORT_Absolute" in c), None)
    score_col = xcell or cib
    if prot_col and score_col:
        d = joined[[prot_col, score_col]].dropna()
        ax.scatter(d[score_col], d[prot_col], s=16, alpha=0.75, c=COLORS["LSCC"])
        if d.shape[0] >= 5:
            rho, p = stats.spearmanr(d[score_col], d[prot_col])
            label = "xCell ImmuneScore" if xcell else "CIBERSORT Absolute"
            ax.set_title(f"TACSTD2 protein vs {label}\nSpearman ρ={rho:.2f} p={p:.3g}")
        ax.set_xlabel(score_col.replace("pheno_", ""))
        ax.set_ylabel("TACSTD2 protein")
    else:
        ax.text(0.5, 0.5, "no xCell/CIBERSORT column", ha="center")
    fig.suptitle("Join is valid. Interpreting this as ICI response is not.", fontsize=12, y=1.03)
    _save(fig, out)


def fig_naive_not_ici(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.8, 4.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    boxes = [
        (0.3, 3.3, 2.6, 2.0, "#dbe7f3", "CPTAC LUAD / LSCC\nTMT protein + RNA\nsurgical, mostly\ntreatment-naive"),
        (3.7, 3.3, 2.6, 2.0, "#f4e1d2", "What you may test\nprotein↔RNA\nprotein↔xCell/CIBERSORT\ntumor vs NAT (limma)"),
        (7.1, 3.3, 2.6, 2.0, "#f6d5d5", "What you may NOT claim\nICI ORR / DCR / RECIST\nPFS/OS on PD-1/PD-L1\nBessede 2024 replication"),
        (0.3, 0.4, 4.4, 2.2, "#e8f0e8", "Need ICI-labeled cohorts\nOAK/POPLAR (often EGA)\nGEO ICI RNA (e.g. GSE135222)\nPRIDE only if ICI + protein IDs"),
        (5.3, 0.4, 4.4, 2.2, "#fff3d6", "CLDN4 special case\nRNA yes, TMT row often no\nDo not impute protein\nDo not call ICI from NAT protein"),
    ]
    for x, y, w, h, c, t in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=c, edgecolor="#333", lw=0.8, transform=ax.transData))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8.5)
    ax.annotate("", xy=(3.6, 4.3), xytext=(3.0, 4.3), arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.annotate("", xy=(7.0, 4.3), xytext=(6.4, 4.3), arrowprops=dict(arrowstyle="->", lw=1.2, color="#c1121f"))
    ax.set_title("Treatment-naive CPTAC protein ≠ ICI response", fontsize=13, pad=8)
    _save(fig, out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protein-luad", type=Path, required=True)
    ap.add_argument("--protein-lscc", type=Path, required=True)
    ap.add_argument("--rna-luad", type=Path, required=True)
    ap.add_argument("--rna-lscc", type=Path, required=True)
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--concord-stats", type=Path, required=True)
    ap.add_argument("--limma", type=Path, required=True)
    ap.add_argument("--joined", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    protein = {
        "LUAD": read_expression_matrix(args.protein_luad),
        "LSCC": read_expression_matrix(args.protein_lscc),
    }
    rna = {
        "LUAD": read_expression_matrix(args.rna_luad),
        "LSCC": read_expression_matrix(args.rna_lscc),
    }
    pairs = pd.read_csv(args.pairs, sep="\t")
    cstats = pd.read_csv(args.concord_stats, sep="\t")
    fit = pd.read_csv(args.limma, sep="\t")
    joined = pd.read_csv(args.joined, sep="\t")

    fig_missingness(protein, rna, args.outdir / "01_cldn4_missingness.png")
    fig_log2(protein, rna, args.outdir / "02_log2_scale_check.png")
    fig_concordance(pairs, cstats, args.outdir / "03_protein_rna_concordance.png")
    fig_limma(fit, args.outdir / "04_limma_tumor_vs_nat.png")
    fig_immune_join(joined, args.outdir / "05_join_immune_deconv.png")
    fig_naive_not_ici(args.outdir / "06_naive_not_ici.png")
    print(f"figures -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
