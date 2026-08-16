#!/usr/bin/env python3
"""Main analysis: tumor-cell TACSTD2/CLDN4 vs CD8 infiltration, TLS
(12-chemokine) score, and ICI response in GSE207422 (MPR) and GSE205335
(RECIST).

Inputs (produced by 00/01/02 scripts):
  /tmp/data_scrna_ici/gse207422_panel_cells.tsv.gz
  /tmp/data_scrna_ici/gse205335_panel_cells.tsv.gz
  /tmp/data_scrna_ici/GSE207422_NSCLC_scRNAseq_metadata.xlsx
  /tmp/data_scrna_ici/GSE205335_Lung_IO_CellIdentity.txt.gz
  results/fable_scrna_ici/gse205335_sample_table.tsv

Outputs under results/fable_scrna_ici/.
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "/tmp/data_scrna_ici"
RES = "/workspace/results/fable_scrna_ici"
PANEL = pd.read_csv("/workspace/scripts/fable_scrna_ici/gene_panel.tsv", sep="\t")

TLS12 = PANEL.loc[PANEL.category == "tls12", "gene"].tolist()
MIN_TUMOR_CELLS = 30
MIN_TOTAL_UMI = 200

LINEAGE_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17",
                    "ELF3", "CDH1", "MUC1"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1", "IGHM"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1", "ITGAX"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}

stats_rows = []


def add_stat(dataset, analysis, group_or_x, metric, n, stat_name, stat, p, note=""):
    stats_rows.append(dict(dataset=dataset, analysis=analysis,
                           comparison=group_or_x, metric=metric, n=n,
                           stat=stat_name, value=stat, p_value=p, note=note))


def lognorm(counts, total):
    """log1p CP10K."""
    return np.log1p(counts / total.values[:, None] * 1e4)


def pseudobulk_tls(df_counts, sample_col, genes):
    """Per-sample pseudobulk log2-CPM mean over TLS genes + CXCL13 alone."""
    agg = df_counts.groupby(sample_col)[genes + ["total_umi"]].sum()
    cpm = agg[genes].div(agg["total_umi"], axis=0) * 1e6
    out = pd.DataFrame(index=agg.index)
    out["tls12_score"] = np.log2(cpm[genes] + 1).mean(axis=1)
    out["cxcl13_log2cpm"] = np.log2(cpm["CXCL13"] + 1)
    return out


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return u, p


def spear(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return r, p, int(m.sum())


# ----------------------------------------------------------------------
# GSE207422 (neoadjuvant anti-PD-1 + chemo NSCLC, MPR endpoint)
# ----------------------------------------------------------------------
print("=== GSE207422 ===", flush=True)
cells7 = pd.read_csv(f"{DATA}/gse207422_panel_cells.tsv.gz", sep="\t")
cells7["sample"] = cells7["barcode"].str.rsplit("_", n=1).str[0]
cells7 = cells7[cells7.total_umi >= MIN_TOTAL_UMI].copy()
genes7 = [g for g in PANEL.gene if g in cells7.columns]

meta7 = pd.read_excel(f"{DATA}/GSE207422_NSCLC_scRNAseq_metadata.xlsx").iloc[:15]
meta7 = meta7.rename(columns={"Sample": "sample"})
meta7["mpr_group"] = meta7["Pathologic Response"].map(
    {"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR", "NE": np.nan})

# marker-score lineage assignment
ln = lognorm(cells7[genes7], cells7["total_umi"])
ln = pd.DataFrame(ln, columns=genes7, index=cells7.index)
scores = pd.DataFrame({lin: ln[[g for g in gs if g in ln.columns]].mean(axis=1)
                       for lin, gs in LINEAGE_MARKERS.items()})
cells7["lineage"] = scores.idxmax(axis=1)
cells7["max_score"] = scores.max(axis=1)
cells7.loc[cells7["max_score"] <= 0, "lineage"] = "Unassigned"

tscore = ln[["CD3D", "CD3E", "TRAC", "CD2"]].mean(axis=1)
cd8score = ln[["CD8A", "CD8B"]].mean(axis=1)
is_t = (cells7["lineage"] == "T/NK") & (tscore > 0)
cells7["is_cd8t"] = is_t & (cd8score > 0)

# classification QC: mean marker expression per lineage
qc = ln[["PTPRC", "EPCAM", "CD3E", "CD8A", "LYZ", "COL1A1", "PECAM1",
         "MS4A1", "TACSTD2", "CLDN4"]].groupby(cells7["lineage"]).mean().round(3)
qc.insert(0, "n_cells", cells7.groupby("lineage").size())
qc.to_csv(f"{RES}/gse207422_lineage_marker_qc.tsv", sep="\t")
print(qc.to_string())

epi7 = cells7[cells7.lineage == "Epithelial"].copy()
epi_ln7 = ln.loc[epi7.index]

samp7 = pd.DataFrame(index=sorted(cells7["sample"].unique()))
samp7["n_cells"] = cells7.groupby("sample").size()
samp7["n_epithelial"] = epi7.groupby("sample").size()
samp7["n_cd8t"] = cells7.groupby("sample")["is_cd8t"].sum()
samp7["cd8_frac_all"] = samp7["n_cd8t"] / samp7["n_cells"]
for g in ["TACSTD2", "CLDN4"]:
    samp7[f"tumor_{g}_mean"] = epi_ln7[g].groupby(epi7["sample"]).mean()
    samp7[f"tumor_{g}_pct_pos"] = (epi7[g] > 0).groupby(epi7["sample"]).mean()
samp7 = samp7.join(pseudobulk_tls(cells7, "sample", TLS12))
samp7 = samp7.merge(meta7[["sample", "Patient", "Resource",
                           "Pathologic Response", "mpr_group", "RECIST"]],
                    left_index=True, right_on="sample").set_index("sample")
enough = samp7["n_epithelial"].fillna(0) >= MIN_TUMOR_CELLS
samp7.loc[~enough, [c for c in samp7.columns if c.startswith("tumor_")]] = np.nan
samp7.round(4).to_csv(f"{RES}/gse207422_sample_metrics.tsv", sep="\t")
print(samp7.round(3).to_string())

post7 = samp7[samp7.Resource == "Post-treatment surgery"]
for metric in ["tumor_TACSTD2_mean", "tumor_TACSTD2_pct_pos",
               "tumor_CLDN4_mean", "tumor_CLDN4_pct_pos",
               "cd8_frac_all", "tls12_score", "cxcl13_log2cpm"]:
    a = post7.loc[post7.mpr_group == "MPR", metric].dropna()
    b = post7.loc[post7.mpr_group == "NMPR", metric].dropna()
    u, p = mwu(a, b)
    add_stat("GSE207422", "MPR_vs_NMPR_posttreat", "MPR vs NMPR", metric,
             f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             f"median MPR={a.median():.3f} NMPR={b.median():.3f}"
             if len(a) and len(b) else "")

for g in ["TACSTD2", "CLDN4"]:
    for target in ["cd8_frac_all", "tls12_score", "cxcl13_log2cpm"]:
        r, p, n = spear(samp7[f"tumor_{g}_mean"].values, samp7[target].values)
        add_stat("GSE207422", "spearman_sample", f"tumor_{g}_mean~{target}",
                 target, n, "spearman_rho", r, p)

# cell-level (exploratory; pseudoreplicated)
epi_post = epi7.merge(samp7[["mpr_group", "Resource"]], left_on="sample",
                      right_index=True)
epi_post = epi_post[epi_post.Resource == "Post-treatment surgery"]
for g in ["TACSTD2", "CLDN4"]:
    vals = np.log1p(epi_post[g] / epi_post["total_umi"] * 1e4)
    a = vals[epi_post.mpr_group == "MPR"]
    b = vals[epi_post.mpr_group == "NMPR"]
    u, p = mwu(a, b)
    add_stat("GSE207422", "cell_level_MPR_vs_NMPR", "MPR vs NMPR",
             f"epi_{g}_lognorm", f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             "pseudoreplicated, exploratory; "
             f"median MPR={a.median():.3f} NMPR={b.median():.3f}")

# ----------------------------------------------------------------------
# GSE205335 (lung cancer ICI, RECIST endpoint, authors' annotations)
# ----------------------------------------------------------------------
print("=== GSE205335 ===", flush=True)
cells5 = pd.read_csv(f"{DATA}/gse205335_panel_cells.tsv.gz", sep="\t")
ident5 = pd.read_csv(f"{DATA}/GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
cells5 = cells5.merge(ident5, on="barcode", how="inner")
samp_tab5 = pd.read_csv(f"{RES}/gse205335_sample_table.tsv", sep="\t")
cells5 = cells5.merge(samp_tab5[["orig.ident", "patient", "tissue", "recist",
                                 "cancer_subtype"]], on="orig.ident")
cells5 = cells5[cells5.total_umi >= MIN_TOTAL_UMI].copy()
genes5 = [g for g in PANEL.gene if g in cells5.columns]

cells5["is_tumor"] = cells5["lineage.sub"] == "Malignant cells"
cells5["is_cd8t"] = cells5["lineage.sub"] == "CD8+ T cells"
cells5["is_cd8tex"] = cells5["celltype"] == "CD8+ TEX"

tum5 = cells5[cells5.is_tumor].copy()
tum_ln5 = pd.DataFrame(lognorm(tum5[["TACSTD2", "CLDN4"]], tum5["total_umi"]),
                       columns=["TACSTD2", "CLDN4"], index=tum5.index)

samp5 = pd.DataFrame(index=sorted(cells5["orig.ident"].unique()))
samp5["n_cells"] = cells5.groupby("orig.ident").size()
samp5["n_tumor"] = tum5.groupby("orig.ident").size()
samp5["n_cd8t"] = cells5.groupby("orig.ident")["is_cd8t"].sum()
samp5["n_cd8tex"] = cells5.groupby("orig.ident")["is_cd8tex"].sum()
samp5["cd8_frac_all"] = samp5["n_cd8t"] / samp5["n_cells"]
samp5["cd8tex_frac_all"] = samp5["n_cd8tex"] / samp5["n_cells"]
for g in ["TACSTD2", "CLDN4"]:
    samp5[f"tumor_{g}_mean"] = tum_ln5[g].groupby(tum5["orig.ident"]).mean()
    samp5[f"tumor_{g}_pct_pos"] = (tum5[g] > 0).groupby(tum5["orig.ident"]).mean()
samp5 = samp5.join(pseudobulk_tls(cells5, "orig.ident", TLS12))
samp5 = samp5.merge(samp_tab5.set_index("orig.ident")[["patient", "tissue",
                    "recist", "cancer_subtype"]], left_index=True,
                    right_index=True)
enough5 = samp5["n_tumor"].fillna(0) >= MIN_TUMOR_CELLS
samp5.loc[~enough5, [c for c in samp5.columns if c.startswith("tumor_")]] = np.nan
samp5["is_tumor_tissue"] = ~samp5["tissue"].str.startswith("Normal")
samp5["response"] = samp5["recist"].map({"PR": "R", "SD": "NR", "PD": "NR"})
samp5.round(4).to_csv(f"{RES}/gse205335_sample_metrics.tsv", sep="\t")
print(samp5.round(3).to_string())

ana5 = samp5[samp5.is_tumor_tissue & samp5.response.notna()]
# patient-level aggregation (samples per patient are not independent)
pat5 = ana5.groupby("patient").agg(
    response=("response", "first"), cancer_subtype=("cancer_subtype", "first"),
    **{c: (c, "mean") for c in
       ["tumor_TACSTD2_mean", "tumor_TACSTD2_pct_pos", "tumor_CLDN4_mean",
        "tumor_CLDN4_pct_pos", "cd8_frac_all", "cd8tex_frac_all",
        "tls12_score", "cxcl13_log2cpm"]})
pat5.round(4).to_csv(f"{RES}/gse205335_patient_metrics.tsv", sep="\t")

for level, dfx in [("patient", pat5), ("sample", ana5)]:
    for metric in ["tumor_TACSTD2_mean", "tumor_TACSTD2_pct_pos",
                   "tumor_CLDN4_mean", "tumor_CLDN4_pct_pos", "cd8_frac_all",
                   "cd8tex_frac_all", "tls12_score", "cxcl13_log2cpm"]:
        a = dfx.loc[dfx.response == "R", metric].dropna()
        b = dfx.loc[dfx.response == "NR", metric].dropna()
        u, p = mwu(a, b)
        add_stat("GSE205335", f"R_vs_NR_{level}", "PR vs SD/PD", metric,
                 f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
                 f"median R={a.median():.3f} NR={b.median():.3f}"
                 if len(a) and len(b) else "")

# NSCLC-only sensitivity (exclude SCLC/NUT)
nsclc = pat5[pat5.cancer_subtype.isin(["ADC", "SQ"])]
for metric in ["tumor_TACSTD2_mean", "tumor_CLDN4_mean", "cd8_frac_all",
               "tls12_score"]:
    a = nsclc.loc[nsclc.response == "R", metric].dropna()
    b = nsclc.loc[nsclc.response == "NR", metric].dropna()
    u, p = mwu(a, b)
    add_stat("GSE205335", "R_vs_NR_patient_NSCLConly", "PR vs SD/PD", metric,
             f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             f"median R={a.median():.3f} NR={b.median():.3f}"
             if len(a) and len(b) else "")

for g in ["TACSTD2", "CLDN4"]:
    for target in ["cd8_frac_all", "cd8tex_frac_all", "tls12_score",
                   "cxcl13_log2cpm"]:
        sub = samp5[samp5.is_tumor_tissue]
        r, p, n = spear(sub[f"tumor_{g}_mean"].values, sub[target].values)
        add_stat("GSE205335", "spearman_sample_tumor_tissue",
                 f"tumor_{g}_mean~{target}", target, n, "spearman_rho", r, p)

# cell-level exploratory
tumx = tum5.merge(samp5[["response", "is_tumor_tissue"]],
                  left_on="orig.ident", right_index=True)
tumx = tumx[tumx.is_tumor_tissue & tumx.response.notna()]
for g in ["TACSTD2", "CLDN4"]:
    vals = np.log1p(tumx[g] / tumx["total_umi"] * 1e4)
    a = vals[tumx.response == "R"]
    b = vals[tumx.response == "NR"]
    u, p = mwu(a, b)
    add_stat("GSE205335", "cell_level_R_vs_NR", "PR vs SD/PD",
             f"tumor_{g}_lognorm", f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             "pseudoreplicated, exploratory; "
             f"median R={a.median():.3f} NR={b.median():.3f}")

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(f"{RES}/stats_summary.tsv", sep="\t", index=False)
print(stats_df.to_string())

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def box_strip(ax, df, group_col, val_col, groups, title):
    data = [df.loc[df[group_col] == g, val_col].dropna() for g in groups]
    ax.boxplot(data, labels=[f"{g}\n(n={len(d)})" for g, d in zip(groups, data)],
               showfliers=False)
    for i, d in enumerate(data):
        ax.scatter(np.random.normal(i + 1, 0.06, len(d)), d, s=25, alpha=0.8,
                   zorder=3)
    ax.set_title(title, fontsize=9)


fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for j, metric in enumerate(["tumor_TACSTD2_mean", "tumor_CLDN4_mean",
                            "cd8_frac_all", "tls12_score"]):
    box_strip(axes[0, j], post7.reset_index(), "mpr_group", metric,
              ["MPR", "NMPR"], f"GSE207422 post-tx\n{metric}")
    box_strip(axes[1, j], pat5.reset_index(), "response", metric,
              ["R", "NR"], f"GSE205335 patient-level\n{metric}")
fig.suptitle("Tumor TACSTD2/CLDN4, CD8 fraction, TLS score vs ICI response")
fig.tight_layout()
fig.savefig(f"{RES}/fig1_response_boxplots.png", dpi=150)

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
pairs = [("tumor_TACSTD2_mean", "cd8_frac_all"),
         ("tumor_TACSTD2_mean", "tls12_score"),
         ("tumor_CLDN4_mean", "cd8_frac_all"),
         ("tumor_CLDN4_mean", "tls12_score")]
for j, (xv, yv) in enumerate(pairs):
    for i, (name, dfx) in enumerate([("GSE207422", samp7),
                                     ("GSE205335 tumor tissues",
                                      samp5[samp5.is_tumor_tissue])]):
        ax = axes[i, j]
        m = dfx[[xv, yv]].dropna()
        ax.scatter(m[xv], m[yv], s=30)
        r, p, n = spear(dfx[xv].values, dfx[yv].values)
        ax.set_xlabel(xv, fontsize=8)
        ax.set_ylabel(yv, fontsize=8)
        ax.set_title(f"{name}\nrho={r:.2f} p={p:.3g} n={n}", fontsize=9)
fig.suptitle("Sample-level Spearman: tumor TACSTD2/CLDN4 vs CD8 / TLS")
fig.tight_layout()
fig.savefig(f"{RES}/fig2_correlations.png", dpi=150)
print("figures written")
