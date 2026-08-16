#!/usr/bin/env python3
"""
Target hunt in GSE180963: Tacstd2 (Trop2) in "cold" KL lung tumors vs immune cells.

Dataset
-------
GSE180963 - "Single cell RNA sequencing of tumor sections from GEMM harboring
KrasG12D/+ or KrasG12D/+Lkb1fl/fl (KL) mutation." (Southern Medical University).

Two samples, one mouse genotype each:
  - GSM5481386 = K  = KrasG12D/+                 (immune-normal / control)
  - GSM5481387 = KL = KrasG12D/+ ; Lkb1fl/fl     (LKB1-loss, "immune-desert"/cold)

LKB1 (Stk11) loss is reported to produce "cold", immune-desert lung tumors. The
question here: in the cold KL tumor microenvironment, is Tacstd2/Trop2 expressed
selectively by the tumor (epithelial) compartment relative to immune cells? That
selectivity is what would make Trop2 attractive as a tumor-directed target (e.g.
the antibody-drug conjugate sacituzumab govitecan).

This script is intentionally conservative and prints/records every honesty caveat
(n=1 per genotype, genotype fully confounded with sample/batch, author-prefiltered
matrices, marker-based automated annotation).
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

warnings.simplefilter("ignore", category=FutureWarning)

DATA = Path("/workspace/data/gse180963")
OUT = Path("/workspace/results/hunt_gse180963")
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
sc.settings.figdir = FIG
sc.settings.set_figure_params(dpi=120, dpi_save=150, frameon=False)
RNG = 0
np.random.seed(RNG)

GENE = "Tacstd2"

# ---------------------------------------------------------------------------
# Marker sets for lineage annotation (mouse symbols). Kept small + canonical.
# ---------------------------------------------------------------------------
MARKERS = {
    "Epithelial/Tumor": ["Epcam", "Sftpc", "Sftpb", "Nkx2-1", "Krt8", "Krt18",
                          "Cldn18", "Sfta2", "Scgb1a1", "Krt19", "Cdh1"],
    "T/NK":             ["Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Trbc2",
                          "Nkg7", "Klrb1c", "Ncr1", "Gzmb"],
    "B":                ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ighm", "Igkc"],
    "Myeloid":          ["Lyz2", "Itgam", "Cd68", "C1qa", "C1qb", "Adgre1",
                          "Csf1r", "S100a8", "S100a9", "Itgax"],
    "Endothelial":      ["Pecam1", "Cldn5", "Cdh5", "Flt1", "Egfl7"],
    "Fibroblast":       ["Col1a1", "Col1a2", "Dcn", "Pdgfra", "Lum"],
}
IMMUNE_LINEAGES = {"T/NK", "B", "Myeloid"}


def load_sample(name, folder):
    d = DATA / folder / folder
    mtx = scipy.io.mmread(d / "matrix.mtx").tocsr()          # genes x cells
    genes = pd.read_csv(d / "genes.tsv", sep="\t", header=None)[0].astype(str).values
    barcodes = pd.read_csv(d / "barcodes.tsv", sep="\t", header=None)[0].astype(str).values
    ad = sc.AnnData(X=sp.csr_matrix(mtx.T))                  # cells x genes
    ad.var_names = genes
    ad.obs_names = [f"{name}_{b}" for b in barcodes]
    ad.obs["sample"] = name
    ad.obs["genotype"] = name
    ad.var_names_make_unique()
    return ad


def main():
    log = {}

    # ---- load + concatenate --------------------------------------------------
    k = load_sample("K", "K")
    kl = load_sample("KL", "KL")
    adata = sc.concat([k, kl], join="inner", label=None)
    adata.obs_names_make_unique()
    log["cells_loaded"] = {"K": int(k.n_obs), "KL": int(kl.n_obs), "total": int(adata.n_obs)}
    log["genes_loaded"] = int(adata.n_vars)

    # ---- QC metrics ----------------------------------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    qc = (adata.obs.groupby("genotype")[["n_genes_by_counts", "total_counts", "pct_counts_mt"]]
          .agg(["median", "min", "max"]))
    qc.to_csv(TAB / "qc_summary.csv")

    # Defensive re-filtering to the authors' stated thresholds (matrices are
    # already filtered, so this should remove very few / no cells).
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=500)
    adata = adata[(adata.obs["n_genes_by_counts"] <= 6000) &
                  (adata.obs["pct_counts_mt"] < 20)].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    log["cells_after_refilter"] = {"before": int(n_before), "after": int(adata.n_obs)}

    # ---- keep raw counts, normalize -----------------------------------------
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # log-normalized full matrix for plotting / DE

    # ---- HVG / PCA / neighbors / UMAP / clustering ---------------------------
    # NOTE: genotype is fully confounded with sample. We deliberately DO NOT
    # batch-integrate, because "correcting" the single K vs single KL sample
    # would erase the very biology we want to describe. Clusters may therefore
    # be partly sample-driven; annotation is done at the lineage level, which is
    # robust to this.
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

    # ---- lineage scoring + per-cluster majority-vote annotation -------------
    present = {ct: [g for g in gs if g in adata.raw.var_names] for ct, gs in MARKERS.items()}
    log["markers_present"] = {ct: gs for ct, gs in present.items()}

    def score_key(ct):  # avoid '/' which h5py disallows in obs keys
        return "score_" + ct.replace("/", "_")

    for ct, gs in present.items():
        sc.tl.score_genes(adata, gs, score_name=score_key(ct), use_raw=True)

    score_cols = [score_key(ct) for ct in MARKERS]
    per_cell = adata.obs[score_cols].values
    adata.obs["cell_lineage_percell"] = [list(MARKERS)[i] for i in per_cell.argmax(1)]

    # cluster label = most common per-cell lineage in the cluster
    cluster_label = (adata.obs.groupby("leiden")["cell_lineage_percell"]
                     .agg(lambda s: s.value_counts().idxmax()))
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype(str)

    # compartment grouping
    def compartment(ct):
        if ct == "Epithelial/Tumor":
            return "Epithelial/Tumor"
        if ct in IMMUNE_LINEAGES:
            return "Immune"
        return "Stromal"
    adata.obs["compartment"] = adata.obs["cell_type"].map(compartment)

    # audit: mean marker-score per cluster
    (adata.obs.groupby("leiden")[score_cols].mean()
     .join(cluster_label.rename("assigned"))
     .to_csv(TAB / "cluster_lineage_scores.csv"))

    # ---- composition tables --------------------------------------------------
    comp = (pd.crosstab(adata.obs["genotype"], adata.obs["cell_type"]))
    comp_frac = comp.div(comp.sum(1), axis=0)
    comp.to_csv(TAB / "celltype_counts_by_genotype.csv")
    comp_frac.round(4).to_csv(TAB / "celltype_fraction_by_genotype.csv")

    compart = pd.crosstab(adata.obs["genotype"], adata.obs["compartment"])
    compart_frac = compart.div(compart.sum(1), axis=0)
    compart_frac.round(4).to_csv(TAB / "compartment_fraction_by_genotype.csv")
    log["immune_fraction"] = {g: float(compart_frac.loc[g, "Immune"])
                              for g in compart_frac.index if "Immune" in compart_frac.columns}

    # =====================================================================
    # Tacstd2 / Trop2 analysis
    # =====================================================================
    lognorm = np.asarray(adata.raw[:, GENE].X.todense()).ravel()
    counts = np.asarray(adata[:, GENE].layers["counts"].todense()).ravel()
    adata.obs[f"{GENE}_lognorm"] = lognorm
    adata.obs[f"{GENE}_pos"] = (counts > 0).astype(int)

    def summarize(df):
        return pd.Series({
            "n_cells": len(df),
            "pct_expressing": 100.0 * (df[f"{GENE}_pos"].mean()),
            "mean_lognorm": df[f"{GENE}_lognorm"].mean(),
            "mean_lognorm_in_pos": df.loc[df[f"{GENE}_pos"] == 1, f"{GENE}_lognorm"].mean(),
        })

    by_ct_geno = (adata.obs.groupby(["genotype", "cell_type"], observed=True)
                  .apply(summarize).reset_index())
    by_ct_geno.to_csv(TAB / f"{GENE}_by_celltype_genotype.csv", index=False)

    by_compart_geno = (adata.obs.groupby(["genotype", "compartment"], observed=True)
                       .apply(summarize).reset_index())
    by_compart_geno.to_csv(TAB / f"{GENE}_by_compartment_genotype.csv", index=False)

    # ---- Core contrast: tumor(epi) vs immune, within KL (the "cold" tumor) ---
    kl_obs = adata.obs[adata.obs["genotype"] == "KL"]
    epi_kl = kl_obs.loc[kl_obs["compartment"] == "Epithelial/Tumor", f"{GENE}_lognorm"].values
    imm_kl = kl_obs.loc[kl_obs["compartment"] == "Immune", f"{GENE}_lognorm"].values

    def contrast(a, b, name_a, name_b):
        res = {
            "group_a": name_a, "n_a": int(len(a)),
            "group_b": name_b, "n_b": int(len(b)),
            "mean_a": float(np.mean(a)) if len(a) else float("nan"),
            "mean_b": float(np.mean(b)) if len(b) else float("nan"),
            "pct_pos_a": float(100 * np.mean(a > 0)) if len(a) else float("nan"),
            "pct_pos_b": float(100 * np.mean(b > 0)) if len(b) else float("nan"),
        }
        res["log2FC_mean(a/b)"] = float(np.log2((res["mean_a"] + 1e-9) / (res["mean_b"] + 1e-9)))
        if len(a) > 0 and len(b) > 0:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            # rank-biserial effect size
            res["mannwhitney_U"] = float(u)
            res["mannwhitney_p"] = float(p)
            res["rank_biserial"] = float(1 - 2 * u / (len(a) * len(b)))
        return res

    contrasts = []
    contrasts.append({"contrast": "KL: Epithelial/Tumor vs Immune",
                      **contrast(epi_kl, imm_kl, "KL Epithelial/Tumor", "KL Immune")})

    # secondary: same contrast in K, and KL vs K epithelial
    k_obs = adata.obs[adata.obs["genotype"] == "K"]
    epi_k = k_obs.loc[k_obs["compartment"] == "Epithelial/Tumor", f"{GENE}_lognorm"].values
    imm_k = k_obs.loc[k_obs["compartment"] == "Immune", f"{GENE}_lognorm"].values
    contrasts.append({"contrast": "K: Epithelial/Tumor vs Immune",
                      **contrast(epi_k, imm_k, "K Epithelial/Tumor", "K Immune")})
    contrasts.append({"contrast": "Epithelial/Tumor: KL vs K (CONFOUNDED w/ sample)",
                      **contrast(epi_kl, epi_k, "KL Epithelial/Tumor", "K Epithelial/Tumor")})
    pd.DataFrame(contrasts).to_csv(TAB / f"{GENE}_contrasts.csv", index=False)
    log["contrasts"] = contrasts

    # =====================================================================
    # Annotation-INDEPENDENT robustness checks (do not rely on clustering)
    # =====================================================================
    def raw_counts(gene):
        return np.asarray(adata[:, gene].layers["counts"].todense()).ravel()

    ptprc = raw_counts("Ptprc")
    adata.obs["Ptprc_pos"] = ptprc > 0
    epcam = raw_counts("Epcam")
    krt_any = np.zeros(adata.n_obs, dtype=bool)
    for g in ["Krt8", "Krt18", "Cldn18", "Krt19"]:
        if g in adata.var_names:
            krt_any |= raw_counts(g) > 0
    # marker-threshold epithelial: epithelial-marked AND not Ptprc(CD45)+
    adata.obs["epi_marker_def"] = (epcam > 0) & krt_any & (~adata.obs["Ptprc_pos"].values)

    # immune-dominance / ambient QC (honesty)
    dom = {}
    for g in adata.obs["genotype"].unique():
        m = adata.obs["genotype"] == g
        dom[g] = {
            "n_cells": int(m.sum()),
            "pct_Ptprc_pos": float(100 * adata.obs.loc[m, "Ptprc_pos"].mean()),
            "pct_Sftpc_pos_AMBIENT_FLAG": float(100 * (raw_counts("Sftpc")[m.values] > 0).mean()),
            "n_epi_marker_def": int(adata.obs.loc[m, "epi_marker_def"].sum()),
        }
    log["immune_dominance_and_ambient"] = dom
    pd.DataFrame(dom).T.to_csv(TAB / "immune_dominance_ambient_qc.csv")

    ind_contrasts = []
    for g in ["KL", "K"]:
        m = (adata.obs["genotype"] == g).values
        lg = adata.obs[f"{GENE}_lognorm"].values
        # non-immune vs immune (Ptprc)
        a = lg[m & (~adata.obs["Ptprc_pos"].values)]
        b = lg[m & (adata.obs["Ptprc_pos"].values)]
        ind_contrasts.append({"contrast": f"{g}: Ptprc- (non-immune) vs Ptprc+ (immune)",
                              **contrast(a, b, f"{g} Ptprc-neg", f"{g} Ptprc-pos")})
        # marker-threshold epithelial vs immune
        a2 = lg[m & adata.obs["epi_marker_def"].values]
        b2 = lg[m & adata.obs["Ptprc_pos"].values]
        ind_contrasts.append({"contrast": f"{g}: epithelial(marker-def) vs Ptprc+ immune",
                              **contrast(a2, b2, f"{g} epi(marker)", f"{g} immune")})
    pd.DataFrame(ind_contrasts).to_csv(TAB / f"{GENE}_contrasts_annotation_independent.csv", index=False)
    log["contrasts_annotation_independent"] = ind_contrasts

    # ---- figures -------------------------------------------------------------
    sc.pl.umap(adata, color=["genotype"], save="_genotype.png", show=False, title="Genotype")
    sc.pl.umap(adata, color=["cell_type"], save="_celltype.png", show=False, title="Cell type (marker-based)")
    sc.pl.umap(adata, color=[GENE], save=f"_{GENE}.png", show=False,
               use_raw=True, color_map="viridis", title=f"{GENE} (log-norm)")
    sc.pl.umap(adata, color=["compartment"], save="_compartment.png", show=False, title="Compartment")

    marker_flat = []
    for gs in present.values():
        marker_flat += gs[:4]
    sc.pl.dotplot(adata, var_names=[g for g in marker_flat if g in adata.raw.var_names],
                  groupby="cell_type", use_raw=True, standard_scale="var",
                  save="_markers.png", show=False)

    sc.pl.violin(adata, keys=[GENE], groupby="cell_type", use_raw=True, rotation=90,
                 save=f"_{GENE}_by_celltype.png", show=False)

    # Tacstd2 by compartment split by genotype
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    order = ["Epithelial/Tumor", "Immune", "Stromal"]
    plotdf = adata.obs[[f"{GENE}_lognorm", "compartment", "genotype"]].copy()
    import seaborn as sns
    sns.violinplot(data=plotdf, x="compartment", y=f"{GENE}_lognorm", hue="genotype",
                   order=order, cut=0, inner="quartile", ax=ax)
    ax.set_title(f"{GENE} (Trop2) log-norm by compartment and genotype")
    ax.set_ylabel(f"{GENE} log-norm expr")
    fig.tight_layout()
    fig.savefig(FIG / f"violin_{GENE}_compartment_genotype.png", dpi=150)
    plt.close(fig)

    # ---- persist -------------------------------------------------------------
    adata.write(OUT / "adata_gse180963_processed.h5ad")
    with open(OUT / "metrics.json", "w") as f:
        json.dump(log, f, indent=2)

    print("\n==== KEY RESULT (KL, cold tumor) ====")
    print(json.dumps(contrasts[0], indent=2))
    print("\nImmune fraction by genotype:", log.get("immune_fraction"))
    print("Done. Outputs in", OUT)


if __name__ == "__main__":
    main()
