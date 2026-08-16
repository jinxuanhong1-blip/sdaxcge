#!/usr/bin/env python3
"""Malignant-restricted TACSTD2 vs T/NK fraction and ICI response.

GSE207422 (Hu et al. Genome Med 2023, PMID 36869384):
  Authors identified malignant epithelium with CopyKAT (stromal reference).
  GEO does not provide per-cell CopyKAT labels, so malignant-like cells are
  approximated as epithelial-lineage cells that lack normal-lung markers
  (alveolar SFTPA1/SFTPA2/SFTPB/AGER, club SCGB1A1, ciliated TPPP3/FOXJ1)
  used by the paper to annotate non-malignant clusters. TACSTD2 is scored
  only in those cells. pCR is grouped with MPR (authors' grouping). Samples
  with <10 malignant-like cells are dropped from TACSTD2 means (authors
  dropped one NMPR sample with <10 malignant cells).

GSE205335 (Park/Ahn/Lee lung ICI atlas):
  Uses the authors' published 'Malignant cells' labels. RECIST PR = R,
  SD/PD = NR. Normal tissues excluded from response tests.

No numbers are invented; all statistics are computed from the processed
UMI matrices downloaded from GEO.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "/tmp/data_scrna_ici"
RES = "/workspace/results/fable_scrna_ici"
os.makedirs(RES, exist_ok=True)
PANEL = pd.read_csv("/workspace/scripts/fable_scrna_ici/gene_panel.tsv", sep="\t")
TLS12 = PANEL.loc[PANEL.category == "tls12", "gene"].tolist()

MIN_UMI = 200
MIN_MALIG = 10          # Hu et al. dropped samples with <10 malignant cells
MIN_MALIG_STRICT = 30

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
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPD", "AGER", "NAPSA",
               "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
TUMOR_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CEACAM5",
             "CEACAM6", "MUC1", "ELF3"]

stats_rows = []


def add_stat(dataset, analysis, comparison, metric, n, stat_name, value, p, note=""):
    stats_rows.append(dict(dataset=dataset, analysis=analysis,
                           comparison=comparison, metric=metric, n=n,
                           stat=stat_name, value=value, p_value=p, note=note))


def log1p_cp10k(counts, total):
    return np.log1p(np.asarray(counts, float) / np.asarray(total, float)[:, None] * 1e4)


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    return stats.mannwhitneyu(a, b, alternative="two-sided")


def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def present(genes, columns):
    return [g for g in genes if g in columns]


# ======================================================================
# GSE207422
# ======================================================================
print("=== GSE207422 ===", flush=True)
cells7 = pd.read_csv(f"{DATA}/gse207422_panel_cells.tsv.gz", sep="\t")
extra_path = f"{DATA}/gse207422_extra_genes.tsv.gz"
if os.path.exists(extra_path):
    extra = pd.read_csv(extra_path, sep="\t")
    newcols = [c for c in extra.columns if c != "barcode" and c not in cells7.columns]
    cells7 = cells7.merge(extra[["barcode"] + newcols], on="barcode", how="left")
    print("merged extra genes:", newcols)
else:
    print("WARNING: extra-gene file missing; normal-lung filter uses panel genes only")

cells7["sample"] = cells7["barcode"].str.rsplit("_", n=1).str[0]
cells7 = cells7.loc[cells7.total_umi >= MIN_UMI].copy()
count_genes = [c for c in cells7.columns
               if c not in ("barcode", "total_umi", "sample")]

meta7 = pd.read_excel(f"{DATA}/GSE207422_NSCLC_scRNAseq_metadata.xlsx").iloc[:15]
meta7 = meta7.rename(columns={"Sample": "sample"})
# authors: pCR counted as MPR; P01 pre-tx is NE
meta7["mpr_group"] = meta7["Pathologic Response"].map(
    {"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR"})
meta7["timing"] = meta7["Resource"].map({
    "Pre-treatment biopsy": "pre", "Post-treatment surgery": "post"})

ln7 = pd.DataFrame(log1p_cp10k(cells7[count_genes], cells7["total_umi"]),
                   columns=count_genes, index=cells7.index)
scores = pd.DataFrame({lin: ln7[present(gs, ln7.columns)].mean(axis=1)
                       for lin, gs in LINEAGE_MARKERS.items()})
cells7["lineage"] = scores.idxmax(axis=1)
cells7["max_score"] = scores.max(axis=1)
cells7.loc[cells7["max_score"] <= 0, "lineage"] = "Unassigned"
cells7["is_tnk"] = cells7["lineage"] == "T/NK"
cd8s = ln7[present(["CD8A", "CD8B"], ln7.columns)].mean(axis=1)
ts = ln7[present(["CD3D", "CD3E", "TRAC", "CD2"], ln7.columns)].mean(axis=1)
cells7["is_cd8t"] = cells7["is_tnk"] & (ts > 0) & (cd8s > 0)

norm_g = present(NORMAL_LUNG, ln7.columns)
tum_g = present(TUMOR_EPI, ln7.columns)
print("normal-lung genes used:", norm_g)
print("tumor-epi genes used:", tum_g)
cells7["normal_lung_score"] = ln7[norm_g].mean(axis=1) if norm_g else 0
cells7["tumor_epi_score"] = ln7[tum_g].mean(axis=1)
# malignant-like: epithelial + tumor-epi program exceeds normal-lung program
cells7["is_epithelial"] = cells7["lineage"] == "Epithelial"
cells7["is_malignant"] = (cells7["is_epithelial"]
                          & (cells7["tumor_epi_score"] > cells7["normal_lung_score"])
                          & (cells7["normal_lung_score"] < 0.4))
cells7["is_normal_epi"] = cells7["is_epithelial"] & ~cells7["is_malignant"]

qc_genes = present(["PTPRC", "EPCAM", "CD3E", "CD8A", "LYZ", "COL1A1",
                    "PECAM1", "MS4A1", "TACSTD2", "CLDN4", "SFTPA1",
                    "SFTPB", "NAPSA"] + norm_g, ln7.columns)
qc = ln7[qc_genes].groupby(cells7["lineage"]).mean().round(3)
qc.insert(0, "n_cells", cells7.groupby("lineage").size())
qc.to_csv(f"{RES}/gse207422_lineage_marker_qc.tsv", sep="\t")
print(qc.to_string())

epi_qc = (ln7.loc[cells7.is_epithelial, present(
    ["TACSTD2", "CLDN4", "EPCAM"] + norm_g, ln7.columns)]
          .groupby(cells7.loc[cells7.is_epithelial, "is_malignant"])
          .mean().round(3))
epi_qc.insert(0, "n_cells",
              cells7.loc[cells7.is_epithelial].groupby("is_malignant").size())
epi_qc.to_csv(f"{RES}/gse207422_malignant_vs_normal_epi_qc.tsv", sep="\t")
print("malignant vs normal epi:\n", epi_qc.to_string())

mal7 = cells7[cells7.is_malignant]
epi7 = cells7[cells7.is_epithelial]
mal_ln = ln7.loc[mal7.index]
epi_ln = ln7.loc[epi7.index]

samp7 = pd.DataFrame(index=sorted(cells7["sample"].unique()))
samp7["n_cells"] = cells7.groupby("sample").size()
samp7["n_epithelial"] = cells7.groupby("sample")["is_epithelial"].sum()
samp7["n_malignant"] = cells7.groupby("sample")["is_malignant"].sum()
samp7["n_normal_epi"] = cells7.groupby("sample")["is_normal_epi"].sum()
samp7["n_tnk"] = cells7.groupby("sample")["is_tnk"].sum()
samp7["n_cd8t"] = cells7.groupby("sample")["is_cd8t"].sum()
samp7["tnk_frac"] = samp7["n_tnk"] / samp7["n_cells"]
samp7["cd8_frac"] = samp7["n_cd8t"] / samp7["n_cells"]
for g in ["TACSTD2", "CLDN4"]:
    samp7[f"malig_{g}_mean"] = mal_ln[g].groupby(mal7["sample"]).mean()
    samp7[f"malig_{g}_pct"] = (mal7[g] > 0).groupby(mal7["sample"]).mean()
    samp7[f"epi_{g}_mean"] = epi_ln[g].groupby(epi7["sample"]).mean()
# TLS 12-chemokine pseudobulk
agg = cells7.groupby("sample")[TLS12 + ["total_umi"]].sum()
cpm = agg[TLS12].div(agg["total_umi"], axis=0) * 1e6
samp7["tls12_score"] = np.log2(cpm[TLS12] + 1).mean(axis=1)
samp7["cxcl13_log2cpm"] = np.log2(cpm["CXCL13"] + 1)
samp7 = samp7.merge(meta7[["sample", "Patient", "Resource", "timing",
                           "Pathologic Response", "mpr_group", "RECIST",
                           "Pathology"]],
                    left_index=True, right_on="sample").set_index("sample")

# mask TACSTD2 means when too few malignant cells (authors used <10)
low = samp7["n_malignant"].fillna(0) < MIN_MALIG
for c in [c for c in samp7.columns if c.startswith("malig_")]:
    samp7.loc[low, c] = np.nan
samp7.round(4).to_csv(f"{RES}/gse207422_sample_metrics.tsv", sep="\t")
print(samp7.round(3).to_string())

# --- tests ---
def report_mwu(df, group_col, g1, g2, metric, analysis, note=""):
    a = df.loc[df[group_col] == g1, metric].dropna()
    b = df.loc[df[group_col] == g2, metric].dropna()
    u, p = mwu(a, b)
    add_stat("GSE207422", analysis, f"{g1} vs {g2}", metric,
             f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             (f"median {g1}={a.median():.4f} {g2}={b.median():.4f}; "
              f"mean {g1}={a.mean():.4f} {g2}={b.mean():.4f}; {note}").strip("; "))
    return a, b, u, p


post = samp7[samp7.timing == "post"]
for metric in ["malig_TACSTD2_mean", "malig_TACSTD2_pct", "malig_CLDN4_mean",
               "epi_TACSTD2_mean", "tnk_frac", "cd8_frac", "tls12_score"]:
    report_mwu(post, "mpr_group", "NMPR", "MPR", metric,
               "post_NMPR_vs_MPR",
               "post-treatment surgery only; pCR counted as MPR")

# direction check for the claimed NMPR > MPR TACSTD2
a, b, u, p = report_mwu(post, "mpr_group", "NMPR", "MPR", "malig_TACSTD2_mean",
                        "CLAIM_NMPR_gt_MPR_malig_TACSTD2",
                        "one-sided claim NMPR>MPR is supported only if median/mean NMPR>MPR")

# Spearman: per-patient malignant TACSTD2 vs T/NK (1 sample / patient)
for subset_name, dfx in [
    ("all_with_malig", samp7),
    ("post_with_malig", post),
    ("post_n_malig_ge30", post[post.n_malignant.fillna(0) >= MIN_MALIG_STRICT]),
]:
    for xv, yv in [("malig_TACSTD2_mean", "tnk_frac"),
                   ("malig_TACSTD2_mean", "cd8_frac"),
                   ("malig_TACSTD2_mean", "tls12_score"),
                   ("malig_CLDN4_mean", "tnk_frac"),
                   ("epi_TACSTD2_mean", "tnk_frac")]:
        r, p, n = spear(dfx[xv], dfx[yv])
        add_stat("GSE207422", f"spearman_{subset_name}", f"{xv}~{yv}",
                 yv, n, "spearman_rho", r, p)

# cell-level exploratory (pseudoreplicated)
mal_post = mal7.merge(samp7[["mpr_group", "timing"]], left_on="sample",
                      right_index=True)
mal_post = mal_post[mal_post.timing == "post"]
vals = np.log1p(mal_post["TACSTD2"] / mal_post["total_umi"] * 1e4)
aa = vals[mal_post.mpr_group == "NMPR"]
bb = vals[mal_post.mpr_group == "MPR"]
u, p = mwu(aa, bb)
add_stat("GSE207422", "cell_level_exploratory", "NMPR vs MPR",
         "malig_TACSTD2_lognorm", f"{len(aa)}v{len(bb)}", "MannWhitneyU", u, p,
         f"PSEUDOREPLICATED; median NMPR={aa.median():.4f} MPR={bb.median():.4f}")

# ======================================================================
# GSE205335 — author malignant labels
# ======================================================================
print("=== GSE205335 ===", flush=True)
cells5 = pd.read_csv(f"{DATA}/gse205335_panel_cells.tsv.gz", sep="\t")
ident5 = pd.read_csv(f"{DATA}/GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
cells5 = cells5.merge(ident5, on="barcode", how="inner")
samp_tab5 = pd.read_csv(f"{RES}/gse205335_sample_table.tsv", sep="\t")
cells5 = cells5.merge(samp_tab5[["orig.ident", "patient", "tissue", "recist",
                                 "cancer_subtype"]], on="orig.ident")
cells5 = cells5.loc[cells5.total_umi >= MIN_UMI].copy()
cells5["is_malignant"] = cells5["lineage.sub"] == "Malignant cells"
cells5["is_tnk"] = cells5["lineage.total"] == "T/NK cells"
cells5["is_cd8t"] = cells5["lineage.sub"] == "CD8+ T cells"

mal5 = cells5[cells5.is_malignant]
mal_ln5 = pd.DataFrame(log1p_cp10k(mal5[["TACSTD2", "CLDN4"]], mal5["total_umi"]),
                       columns=["TACSTD2", "CLDN4"], index=mal5.index)

samp5 = pd.DataFrame(index=sorted(cells5["orig.ident"].unique()))
samp5["n_cells"] = cells5.groupby("orig.ident").size()
samp5["n_malignant"] = cells5.groupby("orig.ident")["is_malignant"].sum()
samp5["n_tnk"] = cells5.groupby("orig.ident")["is_tnk"].sum()
samp5["n_cd8t"] = cells5.groupby("orig.ident")["is_cd8t"].sum()
samp5["tnk_frac"] = samp5["n_tnk"] / samp5["n_cells"]
samp5["cd8_frac"] = samp5["n_cd8t"] / samp5["n_cells"]
for g in ["TACSTD2", "CLDN4"]:
    samp5[f"malig_{g}_mean"] = mal_ln5[g].groupby(mal5["orig.ident"]).mean()
    samp5[f"malig_{g}_pct"] = (mal5[g] > 0).groupby(mal5["orig.ident"]).mean()
agg5 = cells5.groupby("orig.ident")[TLS12 + ["total_umi"]].sum()
cpm5 = agg5[TLS12].div(agg5["total_umi"], axis=0) * 1e6
samp5["tls12_score"] = np.log2(cpm5[TLS12] + 1).mean(axis=1)
samp5["cxcl13_log2cpm"] = np.log2(cpm5["CXCL13"] + 1)
samp5 = samp5.merge(samp_tab5.set_index("orig.ident")[["patient", "tissue",
                    "recist", "cancer_subtype"]], left_index=True,
                    right_index=True)
low5 = samp5["n_malignant"].fillna(0) < MIN_MALIG
for c in [c for c in samp5.columns if c.startswith("malig_")]:
    samp5.loc[low5, c] = np.nan
samp5["is_tumor_tissue"] = ~samp5["tissue"].astype(str).str.startswith("Normal")
samp5["response"] = samp5["recist"].map({"PR": "R", "SD": "NR", "PD": "NR"})
samp5.round(4).to_csv(f"{RES}/gse205335_sample_metrics.tsv", sep="\t")
print(samp5.round(3).to_string())

ana5 = samp5[samp5.is_tumor_tissue & samp5.response.notna()].copy()
pat5 = ana5.groupby("patient").agg(
    response=("response", "first"),
    cancer_subtype=("cancer_subtype", "first"),
    n_samples=("n_cells", "size"),
    n_malignant=("n_malignant", "sum"),
    **{c: (c, "mean") for c in
       ["malig_TACSTD2_mean", "malig_TACSTD2_pct", "malig_CLDN4_mean",
        "malig_CLDN4_pct", "tnk_frac", "cd8_frac", "tls12_score",
        "cxcl13_log2cpm"]})
pat5.round(4).to_csv(f"{RES}/gse205335_patient_metrics.tsv", sep="\t")
print(pat5.round(3).to_string())

for level, dfx, ds_note in [
    ("patient", pat5, "mean of tumor-tissue samples per patient"),
    ("sample", ana5, "tumor tissues only; samples not independent"),
]:
    for metric in ["malig_TACSTD2_mean", "malig_TACSTD2_pct",
                   "malig_CLDN4_mean", "tnk_frac", "cd8_frac", "tls12_score"]:
        a = dfx.loc[dfx.response == "NR", metric].dropna()
        b = dfx.loc[dfx.response == "R", metric].dropna()
        u, p = mwu(a, b)
        add_stat("GSE205335", f"NR_vs_R_{level}", "SD/PD vs PR", metric,
                 f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
                 f"median NR={a.median():.4f} R={b.median():.4f}; {ds_note}"
                 if len(a) and len(b) else ds_note)

nsclc = pat5[pat5.cancer_subtype.isin(["ADC", "SQ"])]
for metric in ["malig_TACSTD2_mean", "tnk_frac", "tls12_score"]:
    a = nsclc.loc[nsclc.response == "NR", metric].dropna()
    b = nsclc.loc[nsclc.response == "R", metric].dropna()
    u, p = mwu(a, b)
    add_stat("GSE205335", "NR_vs_R_patient_NSCLC", "SD/PD vs PR", metric,
             f"{len(a)}v{len(b)}", "MannWhitneyU", u, p,
             f"median NR={a.median():.4f} R={b.median():.4f}; ADC+SQ only"
             if len(a) and len(b) else "ADC+SQ only")

for subset_name, dfx in [
    ("patient_tumor_tissue", pat5),
    ("patient_NSCLC", nsclc),
    ("sample_tumor_tissue", ana5),
]:
    for xv, yv in [("malig_TACSTD2_mean", "tnk_frac"),
                   ("malig_TACSTD2_mean", "cd8_frac"),
                   ("malig_TACSTD2_mean", "tls12_score"),
                   ("malig_CLDN4_mean", "tnk_frac")]:
        r, p, n = spear(dfx[xv], dfx[yv])
        add_stat("GSE205335", f"spearman_{subset_name}", f"{xv}~{yv}",
                 yv, n, "spearman_rho", r, p)

malx = mal5.merge(samp5[["response", "is_tumor_tissue"]],
                  left_on="orig.ident", right_index=True)
malx = malx[malx.is_tumor_tissue & malx.response.notna()]
vals = np.log1p(malx["TACSTD2"] / malx["total_umi"] * 1e4)
aa = vals[malx.response == "NR"]
bb = vals[malx.response == "R"]
u, p = mwu(aa, bb)
add_stat("GSE205335", "cell_level_exploratory", "SD/PD vs PR",
         "malig_TACSTD2_lognorm", f"{len(aa)}v{len(bb)}", "MannWhitneyU", u, p,
         f"PSEUDOREPLICATED; median NR={aa.median():.4f} R={bb.median():.4f}")

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(f"{RES}/stats_summary.tsv", sep="\t", index=False)
print("\n=== STATS ===")
print(stats_df.to_string())

# ======================================================================
# Figures
# ======================================================================
rng = np.random.default_rng(0)


def box_strip(ax, df, group_col, val_col, groups, title):
    data = [df.loc[df[group_col] == g, val_col].dropna().to_numpy()
            for g in groups]
    labels = [f"{g}\n(n={len(d)})" for g, d in zip(groups, data)]
    ax.boxplot([d if len(d) else [np.nan] for d in data], showfliers=False)
    ax.set_xticklabels(labels)
    for i, d in enumerate(data):
        if len(d):
            ax.scatter(rng.normal(i + 1, 0.06, len(d)), d, s=28, alpha=0.85,
                       zorder=3)
    ax.set_title(title, fontsize=9)


fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
box_strip(axes[0, 0], post.reset_index(), "mpr_group", "malig_TACSTD2_mean",
          ["MPR", "NMPR"], "GSE207422 post-tx\nmalignant TACSTD2")
box_strip(axes[0, 1], post.reset_index(), "mpr_group", "tnk_frac",
          ["MPR", "NMPR"], "GSE207422 post-tx\nT/NK fraction")
box_strip(axes[0, 2], post.reset_index(), "mpr_group", "tls12_score",
          ["MPR", "NMPR"], "GSE207422 post-tx\nTLS12 score")
box_strip(axes[1, 0], pat5.reset_index(), "response", "malig_TACSTD2_mean",
          ["R", "NR"], "GSE205335 patient\nmalignant TACSTD2")
box_strip(axes[1, 1], pat5.reset_index(), "response", "tnk_frac",
          ["R", "NR"], "GSE205335 patient\nT/NK fraction")
box_strip(axes[1, 2], pat5.reset_index(), "response", "tls12_score",
          ["R", "NR"], "GSE205335 patient\nTLS12 score")
fig.suptitle("Malignant-restricted TACSTD2 and T/NK vs ICI response")
fig.tight_layout()
fig.savefig(f"{RES}/fig1_response_boxplots.png", dpi=150)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax, (name, dfx, xv, yv) in zip(axes, [
    ("GSE207422 patients (1 sample each)", samp7,
     "malig_TACSTD2_mean", "tnk_frac"),
    ("GSE205335 patients (tumor tissues)", pat5,
     "malig_TACSTD2_mean", "tnk_frac"),
]):
    m = dfx[[xv, yv]].dropna()
    ax.scatter(m[xv], m[yv], s=36)
    r, p, n = spear(dfx[xv], dfx[yv])
    ax.set_xlabel("malignant TACSTD2 mean (log1p CP10K)")
    ax.set_ylabel("T/NK fraction")
    ax.set_title(f"{name}\nSpearman ρ={r:.3f}  p={p:.3g}  n={n}")
fig.suptitle("Per-patient malignant TACSTD2 vs T/NK fraction")
fig.tight_layout()
fig.savefig(f"{RES}/fig2_tacstd2_vs_tnk.png", dpi=150)
print("figures written")
