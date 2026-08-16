#!/usr/bin/env python
"""B6 (extra): CLDN4-high epithelial/tumor regions vs immune niches in public Visium lung.

Dataset: 10x Genomics "Human Lung Cancer (FFPE)" — CytAssist Visium, squamous cell
carcinoma (Space Ranger 2.0.0 outputs, CC BY 4.0). NOT E-MTAB-13530.
Only small processed outputs are used (filtered matrix + spatial folder);
huge raw files (3 GB tissue_image.tif, 1.6 GB cloupe, 377 MB molecule_info.h5,
FASTQs) are deliberately skipped — see SKIPPED_FILES.md in the results folder.

Outputs -> results/w200/B6_visium_extra/
"""

import json
import os
import shutil
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import stats

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
RNG = 0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "CytAssist_FFPE_Human_Lung_SCC")
OUT = os.path.join(ROOT, "results", "w200", "B6_visium_extra")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
for d in (OUT, FIG, TAB):
    os.makedirs(d, exist_ok=True)

PREFIX = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"

# ---------------------------------------------------------------- load
# Arrange files into the standard Space Ranger layout expected by sc.read_visium.
sr_dir = os.path.join(DATA, "spaceranger_layout")
os.makedirs(os.path.join(sr_dir, "spatial"), exist_ok=True)
if not os.path.exists(os.path.join(sr_dir, "filtered_feature_bc_matrix.h5")):
    shutil.copy(
        os.path.join(DATA, f"{PREFIX}_filtered_feature_bc_matrix.h5"),
        os.path.join(sr_dir, "filtered_feature_bc_matrix.h5"),
    )
for f in os.listdir(os.path.join(DATA, "spatial")):
    src = os.path.join(DATA, "spatial", f)
    dst = os.path.join(sr_dir, "spatial", f)
    if not os.path.exists(dst) and not f.endswith(".tiff"):
        shutil.copy(src, dst)

adata = sc.read_visium(sr_dir, count_file="filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
print(f"Loaded: {adata.shape[0]} spots x {adata.shape[1]} genes")

# ---------------------------------------------------------------- QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
n0 = adata.n_obs
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=250)
sc.pp.filter_genes(adata, min_cells=5)
qc_note = f"QC: {n0} -> {adata.n_obs} spots (min_counts=500, min_genes=250); {adata.n_vars} genes kept (min_cells=5)"
print(qc_note)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ---------------------------------------------------------------- clustering
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.pca(adata, n_comps=30, random_state=RNG)
sc.pp.neighbors(adata, n_neighbors=15, random_state=RNG)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden", random_state=RNG,
             flavor="igraph", n_iterations=2, directed=False)
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

# ---------------------------------------------------------------- signatures
immune_sets = {
    "Tcell": ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "IL7R"],
    "Bcell_plasma": ["MS4A1", "CD79A", "CD79B", "IGHM", "MZB1", "JCHAIN", "IGKC"],
    "Myeloid": ["CD68", "LYZ", "CD14", "FCGR3A", "ITGAX", "C1QA", "C1QB", "AIF1"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"],
}
immune_all = sorted({g for gs in immune_sets.values() for g in gs if g in adata.var_names})
panimmune = [g for g in ["PTPRC"] + immune_all if g in adata.var_names]
epi_markers = [g for g in ["EPCAM", "KRT5", "KRT6A", "KRT17", "SFTPC", "SFTPB", "NAPSA"]
               if g in adata.var_names]

sc.tl.score_genes(adata, panimmune, score_name="immune_score", random_state=RNG)
for name, genes in immune_sets.items():
    gs = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, gs, score_name=f"{name}_score", random_state=RNG)
sc.tl.score_genes(adata, epi_markers, score_name="epithelial_score", random_state=RNG)

assert "CLDN4" in adata.var_names, "CLDN4 not in filtered gene set"
adata.obs["CLDN4"] = np.asarray(adata[:, "CLDN4"].X.todense()).ravel()

# ---------------------------------------------------------------- niche calls
# Immune niche: spots in the top quartile of the pan-immune score AND belonging to
# Leiden clusters whose mean immune score is above the 75th percentile of cluster
# means (spot-level signal backed by a coherent cluster — more robust than either alone).
cl_imm = adata.obs.groupby("leiden")["immune_score"].mean()
immune_clusters = cl_imm[cl_imm > cl_imm.quantile(0.75)].index.tolist()
imm_thresh = adata.obs["immune_score"].quantile(0.75)
adata.obs["immune_niche"] = (
    (adata.obs["immune_score"] > imm_thresh)
    & adata.obs["leiden"].isin(immune_clusters)
).values

# CLDN4-high: top quartile of CLDN4 (log-normalized), excluding immune-niche spots
cldn4_thresh = adata.obs["CLDN4"].quantile(0.75)
adata.obs["cldn4_high"] = (adata.obs["CLDN4"] > cldn4_thresh).values

cat = np.where(
    adata.obs["immune_niche"] & adata.obs["cldn4_high"], "CLDN4-high & immune",
    np.where(adata.obs["cldn4_high"], "CLDN4-high",
             np.where(adata.obs["immune_niche"], "Immune niche", "Other")),
)
adata.obs["niche"] = pd.Categorical(
    cat, categories=["CLDN4-high", "Immune niche", "CLDN4-high & immune", "Other"]
)
niche_counts = adata.obs["niche"].value_counts()
print(niche_counts)

# ---------------------------------------------------------------- statistics
res = {}

# Spot-level association CLDN4 vs immune score
rho, p_rho = stats.spearmanr(adata.obs["CLDN4"], adata.obs["immune_score"])
res["spearman_CLDN4_vs_immune_score"] = {"rho": float(rho), "p": float(p_rho)}

# Fisher exact: CLDN4-high x immune-niche co-membership
ct = pd.crosstab(adata.obs["cldn4_high"], adata.obs["immune_niche"])
odds, p_fish = stats.fisher_exact(ct.values)
res["fisher_cldn4high_x_immuneniche"] = {"odds_ratio": float(odds), "p": float(p_fish)}
ct.to_csv(os.path.join(TAB, "contingency_cldn4high_x_immuneniche.csv"))

# CLDN4 expression by niche (Mann-Whitney, CLDN4-high excluded by construction;
# compare immune niche vs other)
mw = stats.mannwhitneyu(
    adata.obs.loc[adata.obs["niche"] == "Immune niche", "CLDN4"],
    adata.obs.loc[adata.obs["niche"] == "Other", "CLDN4"],
)
res["mannwhitney_CLDN4_immune_vs_other"] = {"U": float(mw.statistic), "p": float(mw.pvalue)}

# Immune sub-lineage scores in CLDN4-high vs immune niche vs other
rows = []
for sig in ["immune_score", "Tcell_score", "Bcell_plasma_score", "Myeloid_score",
            "NK_score", "epithelial_score", "CLDN4"]:
    for grp in adata.obs["niche"].cat.categories:
        v = adata.obs.loc[adata.obs["niche"] == grp, sig]
        if len(v):
            rows.append({"signature": sig, "niche": grp, "n": len(v),
                         "mean": v.mean(), "median": v.median()})
pd.DataFrame(rows).to_csv(os.path.join(TAB, "signature_means_by_niche.csv"), index=False)

# ---------------------------------------------------------------- spatial stats
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
sq.gr.nhood_enrichment(adata, cluster_key="niche", seed=RNG)
z = pd.DataFrame(
    adata.uns["niche_nhood_enrichment"]["zscore"],
    index=adata.obs["niche"].cat.categories,
    columns=adata.obs["niche"].cat.categories,
)
z.to_csv(os.path.join(TAB, "nhood_enrichment_zscores.csv"))
res["nhood_zscore_CLDN4high_vs_ImmuneNiche"] = float(z.loc["CLDN4-high", "Immune niche"])
res["nhood_zscore_CLDN4high_self"] = float(z.loc["CLDN4-high", "CLDN4-high"])
res["nhood_zscore_ImmuneNiche_self"] = float(z.loc["Immune niche", "Immune niche"])

# Distance from CLDN4-high spots to nearest immune-niche spot vs random spots
from scipy.spatial import cKDTree

coords = adata.obsm["spatial"].astype(float)
imm_xy = coords[adata.obs["immune_niche"].values]
tree = cKDTree(imm_xy)
d_cldn4, _ = tree.query(coords[(adata.obs["niche"] == "CLDN4-high").values])
d_other, _ = tree.query(coords[(adata.obs["niche"] == "Other").values])
mw_d = stats.mannwhitneyu(d_cldn4, d_other)
# convert px to um via spot spacing: use scalefactor
sf = json.load(open(os.path.join(sr_dir, "spatial", "scalefactors_json.json")))
px_per_um = sf["spot_diameter_fullres"] / 55.0  # Visium spot = 55 um
res["dist_to_immune_um"] = {
    "cldn4_high_median": float(np.median(d_cldn4) / px_per_um),
    "other_median": float(np.median(d_other) / px_per_um),
    "mannwhitney_p": float(mw_d.pvalue),
}

# ---------------------------------------------------------------- DE
sc.tl.rank_genes_groups(
    adata, "niche", groups=["CLDN4-high"], reference="Immune niche", method="wilcoxon"
)
de = sc.get.rank_genes_groups_df(adata, group="CLDN4-high")
de.to_csv(os.path.join(TAB, "DE_CLDN4high_vs_immune_niche.csv"), index=False)
res["top10_up_in_CLDN4high"] = de.head(10)["names"].tolist()
res["top10_up_in_immune_niche"] = de.tail(10)["names"].tolist()[::-1]

# ---------------------------------------------------------------- figures
niche_colors = {
    "CLDN4-high": "#d62728", "Immune niche": "#1f77b4",
    "CLDN4-high & immune": "#9467bd", "Other": "#d3d3d3",
}
adata.uns["niche_colors"] = [niche_colors[c] for c in adata.obs["niche"].cat.categories]

fig, axs = plt.subplots(2, 2, figsize=(16, 14))
sc.pl.spatial(adata, color="CLDN4", ax=axs[0, 0], show=False, size=1.4, cmap="Reds",
              title="CLDN4 (log-norm)")
sc.pl.spatial(adata, color="immune_score", ax=axs[0, 1], show=False, size=1.4,
              cmap="Blues", title="Pan-immune score")
sc.pl.spatial(adata, color="niche", ax=axs[1, 0], show=False, size=1.4,
              title="Niche assignment")
sc.pl.spatial(adata, color="leiden", ax=axs[1, 1], show=False, size=1.4,
              title="Leiden clusters")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "01_spatial_overview.png"), dpi=150, bbox_inches="tight")
plt.close()

fig, axs = plt.subplots(1, 3, figsize=(18, 5))
axs[0].scatter(adata.obs["CLDN4"], adata.obs["immune_score"], s=4, alpha=0.4,
               c=[niche_colors[c] for c in adata.obs["niche"]])
axs[0].set_xlabel("CLDN4 (log-norm)")
axs[0].set_ylabel("Pan-immune score")
axs[0].set_title(f"Spearman rho = {rho:.2f} (p = {p_rho:.1e})")

order = ["CLDN4-high", "CLDN4-high & immune", "Immune niche", "Other"]
data_v = [adata.obs.loc[adata.obs["niche"] == g, "immune_score"] for g in order]
bp = axs[1].violinplot(data_v, showmedians=True)
axs[1].set_xticks(range(1, 5), order, rotation=20)
axs[1].set_ylabel("Pan-immune score")
axs[1].set_title("Immune score by niche")

sub = pd.read_csv(os.path.join(TAB, "signature_means_by_niche.csv"))
piv = sub.pivot(index="signature", columns="niche", values="mean").loc[
    ["Tcell_score", "Bcell_plasma_score", "Myeloid_score", "NK_score",
     "epithelial_score", "CLDN4"]]
im = axs[2].imshow(piv.values, cmap="RdBu_r", aspect="auto",
                   vmin=-np.nanmax(np.abs(piv.values)),
                   vmax=np.nanmax(np.abs(piv.values)))
axs[2].set_xticks(range(piv.shape[1]), piv.columns, rotation=20)
axs[2].set_yticks(range(piv.shape[0]), piv.index)
plt.colorbar(im, ax=axs[2], label="mean score")
axs[2].set_title("Signature means by niche")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "02_cldn4_vs_immune_stats.png"), dpi=150,
            bbox_inches="tight")
plt.close()

fig, axs = plt.subplots(1, 2, figsize=(13, 5))
im = axs[0].imshow(z.values, cmap="coolwarm",
                   vmin=-np.abs(z.values).max(), vmax=np.abs(z.values).max())
axs[0].set_xticks(range(len(z)), z.columns, rotation=20)
axs[0].set_yticks(range(len(z)), z.index)
for i in range(len(z)):
    for j in range(len(z)):
        axs[0].text(j, i, f"{z.values[i, j]:.0f}", ha="center", va="center", fontsize=9)
plt.colorbar(im, ax=axs[0], label="z-score")
axs[0].set_title("Neighborhood enrichment (z)")

axs[1].hist(d_cldn4 / px_per_um, bins=40, alpha=0.6, density=True,
            label=f"CLDN4-high (median {np.median(d_cldn4)/px_per_um:.0f} um)",
            color="#d62728")
axs[1].hist(d_other / px_per_um, bins=40, alpha=0.6, density=True,
            label=f"Other (median {np.median(d_other)/px_per_um:.0f} um)",
            color="#7f7f7f")
axs[1].set_xlabel("Distance to nearest immune-niche spot (um)")
axs[1].set_ylabel("density")
axs[1].legend()
axs[1].set_title(f"Proximity to immune niches (MWU p = {mw_d.pvalue:.1e})")
plt.tight_layout()
plt.savefig(os.path.join(FIG, "03_spatial_relationship.png"), dpi=150,
            bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- save
adata.obs[["leiden", "CLDN4", "immune_score", "Tcell_score", "Bcell_plasma_score",
           "Myeloid_score", "NK_score", "epithelial_score", "cldn4_high",
           "immune_niche", "niche", "total_counts", "n_genes_by_counts"]].to_csv(
    os.path.join(TAB, "spot_annotations.csv"))

res["dataset"] = {
    "name": "10x Genomics Human Lung Cancer (FFPE), CytAssist Visium, Squamous Cell Carcinoma",
    "source": "https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard",
    "spaceranger": "2.0.0",
    "license": "CC BY 4.0",
    "not_emtab13530": True,
}
res["qc"] = qc_note
res["n_spots_by_niche"] = {str(k): int(v) for k, v in niche_counts.items()}
res["immune_clusters_leiden"] = immune_clusters
res["thresholds"] = {"cldn4_lognorm_q75": float(cldn4_thresh),
                     "immune_score_q75": float(imm_thresh)}
with open(os.path.join(OUT, "stats_summary.json"), "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
print("DONE")
