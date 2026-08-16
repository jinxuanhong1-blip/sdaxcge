#!/usr/bin/env python3
"""GSE180963: extract Tacstd2 (TROP2) and Cldn4 in K vs KL lung GEMM scRNA.

Public processed 10x count matrices only (GEO supplementary). n = 2 samples
(1 KrasG12D/+ and 1 KrasG12D/+;Lkb1fl/fl). Every cross-genotype comparison is
confounded with mouse/library; p-values are descriptive (cell-level
pseudoreplication), not inferential.

Cldn4 is scored as a *target gene*, not used to define epithelium (avoids
circularity). Epithelial calling uses Epcam / keratins / Cldn18 / Cdh1.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr

warnings.simplefilter("ignore", category=FutureWarning)

REPO = Path(__file__).resolve().parents[2]
CANDIDATE_DATA = [
    Path("/tmp/gse180963"),
    REPO / "data" / "w200" / "GSE180963",
]
OUT = REPO / "results" / "w200" / "GSE180963"
FIG = OUT / "figures"
TAB = OUT / "tables"
TARGETS = ("Tacstd2", "Cldn4")
RNG = 0

# Lineage markers. Cldn4 is intentionally omitted from epithelium (it is a
# target). Sftpc is listed only as an ambient-RNA flag, not as a caller.
MARKERS = {
    "Epithelial/Tumor": ["Epcam", "Krt8", "Krt18", "Krt19", "Cldn18", "Cdh1", "Nkx2-1"],
    "T/NK": ["Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Trbc2", "Nkg7", "Ncr1", "Gzmb"],
    "B": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ighm"],
    "Myeloid": ["Lyz2", "Itgam", "Cd68", "C1qa", "Adgre1", "Csf1r", "S100a8", "S100a9"],
    "Endothelial": ["Pecam1", "Cldn5", "Cdh5", "Flt1"],
    "Fibroblast": ["Col1a1", "Col1a2", "Dcn", "Pdgfra"],
}
IMMUNE = {"T/NK", "B", "Myeloid"}
# Tight-junction / barrier genes reported alongside Cldn4 (not used for calling).
TJ_GENES = ["Cldn3", "Cldn7", "Ocln", "Tjp1", "Epcam"]


def find_data() -> Path:
    for p in CANDIDATE_DATA:
        if (p / "K" / "matrix.mtx").exists() and (p / "KL" / "matrix.mtx").exists():
            return p
    raise FileNotFoundError(
        "GSE180963 matrices not found. Run scripts/w200_gse180963/01_download.py"
    )


def load_sample(root: Path, name: str) -> sc.AnnData:
    d = root / name
    mtx = scipy.io.mmread(d / "matrix.mtx").tocsr()  # genes x cells
    genes = pd.read_csv(d / "genes.tsv", sep="\t", header=None)[0].astype(str).to_numpy()
    barcodes = pd.read_csv(d / "barcodes.tsv", sep="\t", header=None)[0].astype(str).to_numpy()
    ad = sc.AnnData(X=sp.csr_matrix(mtx.T))
    ad.var_names = genes
    ad.obs_names = [f"{name}_{b}" for b in barcodes]
    ad.obs["sample"] = name
    ad.obs["genotype"] = name
    ad.var_names_make_unique()
    return ad


def expr_vec(adata: sc.AnnData, gene: str, layer: str | None = None) -> np.ndarray:
    if layer is None:
        x = adata.raw[:, gene].X
    else:
        x = adata[:, gene].layers[layer]
    if sp.issparse(x):
        return np.asarray(x.todense()).ravel()
    return np.asarray(x).ravel()


def summarize_gene(obs: pd.DataFrame, gene: str) -> pd.Series:
    pos = obs[f"{gene}_pos"]
    return pd.Series(
        {
            "n_cells": int(len(obs)),
            "n_pos": int(pos.sum()),
            "pct_pos": float(100.0 * pos.mean()),
            "mean_lognorm": float(obs[f"{gene}_lognorm"].mean()),
            "median_lognorm": float(obs[f"{gene}_lognorm"].median()),
        }
    )


def contrast(obs: pd.DataFrame, gene: str, mask_a, mask_b, label: str) -> dict:
    a = obs.loc[mask_a, f"{gene}_lognorm"].to_numpy()
    b = obs.loc[mask_b, f"{gene}_lognorm"].to_numpy()
    pos_a = obs.loc[mask_a, f"{gene}_pos"].to_numpy()
    pos_b = obs.loc[mask_b, f"{gene}_pos"].to_numpy()
    mean_a = float(a.mean()) if len(a) else float("nan")
    mean_b = float(b.mean()) if len(b) else float("nan")
    # +eps so zeros do not explode; this is a descriptive fold, not a model.
    log2fc = float(np.log2((mean_a + 1e-3) / (mean_b + 1e-3)))
    if len(a) and len(b):
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial: positive => a > b
        rbc = float(1.0 - (2.0 * u) / (len(a) * len(b)))
    else:
        u, p, rbc = float("nan"), float("nan"), float("nan")
    return {
        "gene": gene,
        "contrast": label,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": mean_a,
        "mean_b": mean_b,
        "pct_pos_a": float(100.0 * pos_a.mean()) if len(pos_a) else float("nan"),
        "pct_pos_b": float(100.0 * pos_b.mean()) if len(pos_b) else float("nan"),
        "log2FC_mean(a/b)": log2fc,
        "mannwhitney_U": float(u),
        "mannwhitney_p": float(p),
        "rank_biserial_a_gt_b": rbc,
        "note": "cell-level p; n=1 mouse/genotype; descriptive only",
    }


def savefig(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIG / name, dpi=150)
    plt.close()


def main() -> None:
    np.random.seed(RNG)
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 1
    sc.settings.figdir = FIG
    sc.settings.set_figure_params(dpi=120, dpi_save=150, frameon=False)

    data_root = find_data()
    log: dict = {
        "series": "GSE180963",
        "data_root": str(data_root),
        "n_samples": 2,
        "n_mice_per_genotype": 1,
        "targets": list(TARGETS),
        "honesty": [
            "n=2 (1 K + 1 KL). Genotype is fully confounded with mouse and the single pooled 10x library.",
            "Cell-level p-values are pseudoreplication and are reported as descriptive only.",
            "Matrices are author-prefiltered (500-6000 genes, mito<20%). No SoupX/CellBender possible.",
            "Cldn4 is not used to define epithelium.",
        ],
    }

    k = load_sample(data_root, "K")
    kl = load_sample(data_root, "KL")
    adata = sc.concat([k, kl], join="inner")
    adata.obs_names_make_unique()
    log["cells_loaded"] = {"K": int(k.n_obs), "KL": int(kl.n_obs), "total": int(adata.n_obs)}
    log["genes_loaded"] = int(adata.n_vars)
    log["targets_present"] = {g: bool(g in adata.var_names) for g in TARGETS}
    missing = [g for g in TARGETS if g not in adata.var_names]
    if missing:
        raise SystemExit(f"target genes missing from matrix: {missing}")

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    (
        adata.obs.groupby("genotype")[["n_genes_by_counts", "total_counts", "pct_counts_mt"]]
        .agg(["median", "min", "max"])
        .to_csv(TAB / "qc_summary.csv")
    )

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=500)
    adata = adata[(adata.obs["n_genes_by_counts"] <= 6000) & (adata.obs["pct_counts_mt"] < 20)].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    log["cells_after_refilter"] = {"before": int(n_before), "after": int(adata.n_obs)}

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # No batch integration: genotype == sample. Integrating would erase the contrast.
    proc = adata.copy()
    sc.pp.highly_variable_genes(proc, n_top_genes=2000, flavor="seurat")
    proc = proc[:, proc.var["highly_variable"]].copy()
    sc.pp.scale(proc, max_value=10)
    sc.tl.pca(proc, n_comps=50, svd_solver="arpack", random_state=RNG)
    sc.pp.neighbors(proc, n_neighbors=15, n_pcs=30, random_state=RNG)
    sc.tl.umap(proc, random_state=RNG)
    sc.tl.leiden(proc, resolution=1.0, random_state=RNG, flavor="igraph", n_iterations=2, directed=False)
    adata.obs["leiden"] = proc.obs["leiden"].values
    adata.obsm["X_umap"] = proc.obsm["X_umap"]
    log["n_leiden_clusters"] = int(adata.obs["leiden"].nunique())

    present = {ct: [g for g in gs if g in adata.raw.var_names] for ct, gs in MARKERS.items()}
    log["markers_present"] = present
    for ct, gs in present.items():
        sc.tl.score_genes(adata, gs, score_name="score_" + ct.replace("/", "_"), use_raw=True)
    score_cols = ["score_" + ct.replace("/", "_") for ct in MARKERS]
    adata.obs["lineage_percell"] = [list(MARKERS)[i] for i in adata.obs[score_cols].to_numpy().argmax(1)]
    cluster_label = adata.obs.groupby("leiden")["lineage_percell"].agg(lambda s: s.value_counts().idxmax())
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype(str)
    adata.obs["compartment"] = adata.obs["cell_type"].map(
        lambda ct: "Epithelial/Tumor" if ct == "Epithelial/Tumor" else ("Immune" if ct in IMMUNE else "Stromal")
    )
    (
        adata.obs.groupby("leiden")[score_cols]
        .mean()
        .join(cluster_label.rename("assigned"))
        .to_csv(TAB / "cluster_lineage_scores.csv")
    )

    # Clustering-free epithelial rule (no Cldn4, no Sftpc).
    epcam = expr_vec(adata, "Epcam") if "Epcam" in adata.raw.var_names else np.zeros(adata.n_obs)
    ptprc = expr_vec(adata, "Ptprc") if "Ptprc" in adata.raw.var_names else np.zeros(adata.n_obs)
    krt = np.zeros(adata.n_obs)
    for g in ("Krt8", "Krt18", "Krt19", "Cldn18"):
        if g in adata.raw.var_names:
            krt = np.maximum(krt, expr_vec(adata, g))
    adata.obs["Epcam_pos"] = (expr_vec(adata, "Epcam", layer="counts") > 0).astype(int) if "Epcam" in adata.var_names else 0
    adata.obs["Ptprc_pos"] = (expr_vec(adata, "Ptprc", layer="counts") > 0).astype(int) if "Ptprc" in adata.var_names else 0
    adata.obs["Sftpc_pos"] = (expr_vec(adata, "Sftpc", layer="counts") > 0).astype(int) if "Sftpc" in adata.var_names else 0
    adata.obs["epi_marker"] = ((epcam > 0) & (krt > 0) & (ptprc == 0)).astype(int)

    counts_ct = pd.crosstab(adata.obs["genotype"], adata.obs["cell_type"])
    frac_ct = counts_ct.div(counts_ct.sum(1), axis=0)
    counts_ct.to_csv(TAB / "celltype_counts_by_genotype.csv")
    frac_ct.round(4).to_csv(TAB / "celltype_fraction_by_genotype.csv")
    counts_comp = pd.crosstab(adata.obs["genotype"], adata.obs["compartment"])
    frac_comp = counts_comp.div(counts_comp.sum(1), axis=0)
    counts_comp.to_csv(TAB / "compartment_counts_by_genotype.csv")
    frac_comp.round(4).to_csv(TAB / "compartment_fraction_by_genotype.csv")
    log["immune_fraction_cluster"] = {g: float(frac_comp.loc[g, "Immune"]) for g in frac_comp.index if "Immune" in frac_comp.columns}

    ambient_rows = []
    for g, sub in adata.obs.groupby("genotype"):
        ambient_rows.append(
            {
                "genotype": g,
                "n_cells": int(len(sub)),
                "pct_Ptprc_pos": float(100.0 * sub["Ptprc_pos"].mean()),
                "pct_Sftpc_pos_AMBIENT_FLAG": float(100.0 * sub["Sftpc_pos"].mean()),
                "n_epi_cluster": int((sub["compartment"] == "Epithelial/Tumor").sum()),
                "n_epi_marker": int(sub["epi_marker"].sum()),
                "n_immune_cluster": int((sub["compartment"] == "Immune").sum()),
            }
        )
    pd.DataFrame(ambient_rows).to_csv(TAB / "immune_dominance_ambient_qc.csv", index=False)
    log["ambient_qc"] = ambient_rows

    # Target gene scores
    for gene in TARGETS:
        adata.obs[f"{gene}_lognorm"] = expr_vec(adata, gene)
        adata.obs[f"{gene}_pos"] = (expr_vec(adata, gene, layer="counts") > 0).astype(int)

    by_ct = []
    by_comp = []
    for gene in TARGETS:
        by_ct.append(
            adata.obs.groupby(["genotype", "cell_type"], observed=True)
            .apply(lambda d: summarize_gene(d, gene), include_groups=False)
            .reset_index()
            .assign(gene=gene)
        )
        by_comp.append(
            adata.obs.groupby(["genotype", "compartment"], observed=True)
            .apply(lambda d: summarize_gene(d, gene), include_groups=False)
            .reset_index()
            .assign(gene=gene)
        )
    pd.concat(by_ct, ignore_index=True).to_csv(TAB / "targets_by_celltype_genotype.csv", index=False)
    pd.concat(by_comp, ignore_index=True).to_csv(TAB / "targets_by_compartment_genotype.csv", index=False)

    contrasts = []
    for gene in TARGETS:
        for geno in ("K", "KL"):
            gmask = adata.obs["genotype"] == geno
            contrasts.append(
                contrast(
                    adata.obs,
                    gene,
                    gmask & (adata.obs["compartment"] == "Epithelial/Tumor"),
                    gmask & (adata.obs["compartment"] == "Immune"),
                    f"{geno}: Epithelial/Tumor vs Immune (cluster)",
                )
            )
            contrasts.append(
                contrast(
                    adata.obs,
                    gene,
                    gmask & (adata.obs["epi_marker"] == 1),
                    gmask & (adata.obs["Ptprc_pos"] == 1),
                    f"{geno}: epi(marker) vs Ptprc+ (annotation-independent)",
                )
            )
        contrasts.append(
            contrast(
                adata.obs,
                gene,
                (adata.obs["genotype"] == "KL") & (adata.obs["compartment"] == "Epithelial/Tumor"),
                (adata.obs["genotype"] == "K") & (adata.obs["compartment"] == "Epithelial/Tumor"),
                "Epithelial/Tumor: KL vs K (CONFOUNDED with sample)",
            )
        )
    pd.DataFrame(contrasts).to_csv(TAB / "target_contrasts.csv", index=False)
    log["contrasts"] = contrasts

    # Tacstd2–Cldn4 co-expression inside epithelium
    co_rows = []
    for geno in ("K", "KL"):
        for defn, mask in (
            ("cluster", (adata.obs["genotype"] == geno) & (adata.obs["compartment"] == "Epithelial/Tumor")),
            ("marker", (adata.obs["genotype"] == geno) & (adata.obs["epi_marker"] == 1)),
        ):
            sub = adata.obs.loc[mask]
            n = int(len(sub))
            if n == 0:
                continue
            both = int(((sub["Tacstd2_pos"] == 1) & (sub["Cldn4_pos"] == 1)).sum())
            t_only = int(((sub["Tacstd2_pos"] == 1) & (sub["Cldn4_pos"] == 0)).sum())
            c_only = int(((sub["Tacstd2_pos"] == 0) & (sub["Cldn4_pos"] == 1)).sum())
            neither = n - both - t_only - c_only
            if n >= 5:
                rho, rp = spearmanr(sub["Tacstd2_lognorm"], sub["Cldn4_lognorm"])
            else:
                rho, rp = float("nan"), float("nan")
            co_rows.append(
                {
                    "genotype": geno,
                    "epithelial_definition": defn,
                    "n_epi": n,
                    "n_Tacstd2_pos": int(sub["Tacstd2_pos"].sum()),
                    "n_Cldn4_pos": int(sub["Cldn4_pos"].sum()),
                    "n_both_pos": both,
                    "n_Tacstd2_only": t_only,
                    "n_Cldn4_only": c_only,
                    "n_neither": neither,
                    "pct_both_of_epi": float(100.0 * both / n),
                    "pct_Cldn4_among_Tacstd2pos": float(100.0 * both / sub["Tacstd2_pos"].sum())
                    if sub["Tacstd2_pos"].sum()
                    else float("nan"),
                    "spearman_rho": float(rho),
                    "spearman_p": float(rp),
                }
            )
    pd.DataFrame(co_rows).to_csv(TAB / "tacstd2_cldn4_coexpression_epithelial.csv", index=False)
    log["coexpression"] = co_rows

    # Tight-junction companion genes in KL epithelium (descriptive)
    tj_rows = []
    for gene in TJ_GENES:
        if gene not in adata.raw.var_names:
            continue
        logn = expr_vec(adata, gene)
        pos = expr_vec(adata, gene, layer="counts") > 0
        for geno in ("K", "KL"):
            mask = (adata.obs["genotype"] == geno) & (adata.obs["compartment"] == "Epithelial/Tumor")
            tj_rows.append(
                {
                    "gene": gene,
                    "genotype": geno,
                    "n_epi_cluster": int(mask.sum()),
                    "pct_pos": float(100.0 * pos[mask.to_numpy()].mean()) if mask.sum() else float("nan"),
                    "mean_lognorm": float(logn[mask.to_numpy()].mean()) if mask.sum() else float("nan"),
                }
            )
    pd.DataFrame(tj_rows).to_csv(TAB / "tj_companion_genes_epithelial.csv", index=False)

    # ---- figures -------------------------------------------------------------
    umap = adata.obsm["X_umap"]
    pal = {"K": "#4C78A8", "KL": "#F58518", "Epithelial/Tumor": "#E45756", "Immune": "#54A24B", "Stromal": "#B279A2"}

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for g, c in (("K", pal["K"]), ("KL", pal["KL"])):
        m = adata.obs["genotype"] == g
        ax.scatter(umap[m, 0], umap[m, 1], s=2, c=c, label=g, linewidths=0, alpha=0.6)
    ax.legend(markerscale=4, frameon=False)
    ax.set_title("GSE180963 UMAP by genotype (n=1 each)")
    ax.set_xticks([])
    ax.set_yticks([])
    savefig("umap_genotype.png")

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for g, c in pal.items():
        if g in ("K", "KL"):
            continue
        m = adata.obs["compartment"] == g
        ax.scatter(umap[m, 0], umap[m, 1], s=2, c=c, label=g, linewidths=0, alpha=0.6)
    ax.legend(markerscale=4, frameon=False)
    ax.set_title("UMAP by compartment")
    ax.set_xticks([])
    ax.set_yticks([])
    savefig("umap_compartment.png")

    for gene in TARGETS:
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        v = adata.obs[f"{gene}_lognorm"].to_numpy()
        order = np.argsort(v)
        sca = ax.scatter(umap[order, 0], umap[order, 1], c=v[order], s=2, cmap="magma", linewidths=0)
        fig.colorbar(sca, ax=ax, label=f"{gene} log1p CP10k")
        ax.set_title(f"UMAP {gene}")
        ax.set_xticks([])
        ax.set_yticks([])
        savefig(f"umap_{gene}.png")

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0), sharey=True)
    for ax, gene in zip(axes, TARGETS):
        data, labels, colors = [], [], []
        for geno, hatch in (("K", False), ("KL", True)):
            for comp, c in (("Epithelial/Tumor", pal["Epithelial/Tumor"]), ("Immune", pal["Immune"])):
                vals = adata.obs.loc[
                    (adata.obs["genotype"] == geno) & (adata.obs["compartment"] == comp),
                    f"{gene}_lognorm",
                ]
                data.append(vals)
                labels.append(f"{geno}\n{comp.split('/')[0]}")
                colors.append(c)
        parts = ax.violinplot(data, showmeans=True, showextrema=False)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(colors[i])
            body.set_alpha(0.7)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(gene)
        ax.set_ylabel("log1p CP10k" if gene == "Tacstd2" else "")
    fig.suptitle("Target expression by compartment × genotype (n=1 mouse each)")
    savefig("violin_targets_compartment_genotype.png")

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0), sharex=True, sharey=True)
    for ax, geno in zip(axes, ("K", "KL")):
        m = (adata.obs["genotype"] == geno) & (adata.obs["compartment"] == "Epithelial/Tumor")
        ax.scatter(
            adata.obs.loc[m, "Tacstd2_lognorm"],
            adata.obs.loc[m, "Cldn4_lognorm"],
            s=12,
            c=pal[geno],
            alpha=0.7,
            linewidths=0,
        )
        ax.set_title(f"{geno} epithelium (cluster, n={int(m.sum())})")
        ax.set_xlabel("Tacstd2 log1p CP10k")
        if geno == "K":
            ax.set_ylabel("Cldn4 log1p CP10k")
    fig.suptitle("Tacstd2 vs Cldn4 in cluster-called epithelium")
    savefig("scatter_tacstd2_vs_cldn4_epithelial.png")

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    frac_comp.plot(kind="bar", stacked=True, ax=ax, color=[pal[c] for c in frac_comp.columns], rot=0)
    ax.set_ylabel("fraction of cells")
    ax.set_title("Compartment mix (immune-dominated; n=1 each)")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1))
    savefig("bar_compartment_fraction.png")

    # compact numeric summary figure
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.axis("off")
    lines = ["GSE180963 Tacstd2 / Cldn4 — key numbers (descriptive; n=1 mouse/genotype)"]
    for c in contrasts:
        if "Epithelial/Tumor vs Immune (cluster)" in c["contrast"]:
            lines.append(
                f"{c['gene']} {c['contrast'].split(':')[0]}: "
                f"epi {c['pct_pos_a']:.1f}% pos (n={c['n_a']}) vs immune {c['pct_pos_b']:.1f}% "
                f"(n={c['n_b']}); log2FC={c['log2FC_mean(a/b)']:.2f}"
            )
    ax.text(0.01, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=8)
    savefig("key_numbers.png")

    (OUT / "key_stats.json").write_text(json.dumps(log, indent=2) + "\n")
    print(json.dumps({k: log[k] for k in ("cells_loaded", "targets_present", "immune_fraction_cluster", "ambient_qc")}, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
