#!/usr/bin/env python3
"""Malignant/epithelial-restricted TACSTD2 and CLDN4 vs T/NK (CD8, TLS) and ICI response.

GSE207422 (Hu et al. Genome Med 2023, PMID 36869384)
  Neoadjuvant PD-1 + chemo, pathologic response (MPR/pCR vs NMPR).
  GEO has sample-level clinical metadata only. Authors used CopyKAT for
  malignant epithelium; those labels are not on GEO. Malignant-like cells
  are a marker proxy: epithelial-lineage AND low normal-lung (alveolar /
  club / ciliated) score. pCR is grouped with MPR (authors' grouping).
  Primary test = 12 post-treatment surgical samples. Samples with <10
  malignant-like cells are dropped from TACSTD2/CLDN4 means.

GSE205335 (Park/Ahn/Lee lung ICI atlas)
  Uses the authors' published `Malignant cells` labels (lineage.sub).
  RECIST PR = R, SD/PD = NR. No MPR field in SOFT or the cell table.
  Normal tissues excluded from response tests. Patient is the unit.

All statistics are computed from the downloaded processed UMI matrices.
Cell-level tests are labeled exploratory (pseudoreplication).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/scrna_data")
RES = Path("/workspace/results/scrna")
RES.mkdir(parents=True, exist_ok=True)

MIN_UMI = 200
MIN_MALIG = 10
MIN_MALIG_PT = 20  # GSE205335 patient-level gate (matches published-label analyses)

LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17",
                   "ELF3", "CDH1", "MUC1"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1", "IGHM"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1", "ITGAX"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "AGER", "NAPSA",
               "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
CD8_MARKERS = ["CD8A", "CD8B"]
TLS_CHEMO = ["CXCL13", "CCL19", "CCL21", "CXCL9", "CXCL10"]

stats_rows: list[dict] = []


def add_stat(dataset, analysis, comparison, metric, n, stat_name, value, p, note=""):
    stats_rows.append(dict(
        dataset=dataset, analysis=analysis, comparison=comparison, metric=metric,
        n=n, stat=stat_name, value=value, p_value=p, note=note,
    ))


def present(genes, columns):
    return [g for g in genes if g in columns]


def log1p_cp10k(counts: pd.DataFrame, total: pd.Series) -> pd.DataFrame:
    tot = np.asarray(total, float)
    tot[tot <= 0] = np.nan
    return pd.DataFrame(
        np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4),
        index=counts.index, columns=counts.columns,
    )


def mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, len(a), len(b)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u), float(p), len(a), len(b)


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_r(r):
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def annotate_box(ax, a, b, p):
    ymax = max(np.nanmax(a) if len(a) else 0, np.nanmax(b) if len(b) else 0)
    ax.text(0.5, 0.98, f"MWU p={fmt_p(p)}", transform=ax.transAxes,
            ha="center", va="top", fontsize=9)


# ======================================================================
# GSE207422
# ======================================================================
print("=== GSE207422 ===", flush=True)
cells = pd.read_csv(DATA / "gse207422_panel_cells.tsv.gz", sep="\t")
cells["sample"] = cells["barcode"].str.rsplit("_", n=1).str[0]
n_raw = len(cells)
cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
print(f"cells raw={n_raw} after UMI>={MIN_UMI}: {len(cells)}", flush=True)

count_genes = [c for c in cells.columns if c not in ("barcode", "total_umi", "sample")]
ln = log1p_cp10k(cells[count_genes], cells["total_umi"])

scores = pd.DataFrame(
    {lin: ln[present(gs, ln.columns)].mean(axis=1) for lin, gs in LINEAGE.items()},
    index=cells.index,
)
cells["lineage"] = scores.idxmax(axis=1)
cells["max_score"] = scores.max(axis=1)
nl_cols = present(NORMAL_LUNG, ln.columns)
cells["normal_lung"] = ln[nl_cols].mean(axis=1) if nl_cols else 0.0
# malignant-like: epithelial lineage, not dominated by normal-lung program
# threshold: normal_lung score below the 75th percentile of epithelial cells
epi_mask = cells["lineage"] == "Epithelial"
nl_cut = float(cells.loc[epi_mask, "normal_lung"].quantile(0.75)) if epi_mask.any() else 0.0
cells["malignant_like"] = epi_mask & (cells["normal_lung"] <= nl_cut)
cells["is_tnk"] = cells["lineage"] == "T/NK"
cells["is_b"] = cells["lineage"] == "B/Plasma"
cd8_cols = present(CD8_MARKERS, ln.columns)
cells["cd8_score"] = ln[cd8_cols].mean(axis=1) if cd8_cols else 0.0
# CD8-like among T/NK: CD8 score > CD4 if present
if "CD4" in ln.columns and cd8_cols:
    cells["is_cd8"] = cells["is_tnk"] & (cells["cd8_score"] > ln["CD4"])
else:
    cells["is_cd8"] = cells["is_tnk"] & (cells["cd8_score"] > cells.loc[cells["is_tnk"], "cd8_score"].median())
tls_cols = present(TLS_CHEMO, ln.columns)
cells["tls_chemo"] = ln[tls_cols].mean(axis=1) if tls_cols else 0.0

print("lineage counts:\n", cells["lineage"].value_counts().to_string())
print(f"malignant_like n={int(cells.malignant_like.sum())}  nl_cut={nl_cut:.4f}")
print(f"T/NK n={int(cells.is_tnk.sum())}  CD8-like n={int(cells.is_cd8.sum())}  B/Plasma n={int(cells.is_b.sum())}")

# sanity: TACSTD2 should be high in epithelial / malignant-like, low in T/NK
for gene in ("TACSTD2", "CLDN4"):
    if gene not in ln.columns:
        print(f"WARNING: {gene} missing from matrix")
        continue
    cells[f"{gene}_log1p"] = ln[gene]
    cells[f"{gene}_pos"] = cells[gene] > 0

meta = pd.read_excel(DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
meta = meta.rename(columns={"Sample": "sample"})
meta["mpr_group"] = meta["Pathologic Response"].map({"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR"})
meta["timing"] = meta["Resource"].map({
    "Pre-treatment biopsy": "pre",
    "Post-treatment surgery": "post",
})
meta["recist_bin"] = meta["RECIST"].map({"PR": "PR_CR", "CR": "PR_CR", "SD": "SD", "PD": "SD"})

rows = []
for sample, g in cells.groupby("sample"):
    n = len(g)
    n_mal = int(g.malignant_like.sum())
    n_epi = int((g.lineage == "Epithelial").sum())
    n_tnk = int(g.is_tnk.sum())
    n_cd8 = int(g.is_cd8.sum())
    n_b = int(g.is_b.sum())
    rec = {
        "sample": sample,
        "n_cells": n,
        "n_malignant_like": n_mal,
        "n_epithelial": n_epi,
        "n_tnk": n_tnk,
        "n_cd8": n_cd8,
        "n_b_plasma": n_b,
        "frac_tnk": n_tnk / n,
        "frac_cd8": n_cd8 / n,
        "frac_b_plasma": n_b / n,
        "frac_tls_proxy": (n_b + n_cd8) / n,
        "mean_tls_chemo": float(g.tls_chemo.mean()),
    }
    mal = g.loc[g.malignant_like]
    for gene in ("TACSTD2", "CLDN4"):
        col = f"{gene}_log1p"
        if col not in g.columns:
            rec[f"mal_{gene}_mean"] = np.nan
            rec[f"mal_{gene}_pct_pos"] = np.nan
            continue
        if n_mal >= MIN_MALIG:
            rec[f"mal_{gene}_mean"] = float(mal[col].mean())
            rec[f"mal_{gene}_pct_pos"] = float((mal[gene] > 0).mean() * 100)
        else:
            rec[f"mal_{gene}_mean"] = np.nan
            rec[f"mal_{gene}_pct_pos"] = np.nan
        rec[f"tnk_{gene}_mean"] = float(g.loc[g.is_tnk, col].mean()) if n_tnk else np.nan
    rows.append(rec)

samp = pd.DataFrame(rows).merge(meta, on="sample", how="left")
samp.to_csv(RES / "gse207422_sample_table.tsv", sep="\t", index=False)
print(samp[["sample", "Patient", "timing", "mpr_group", "RECIST", "n_cells",
            "n_malignant_like", "n_tnk", "mal_TACSTD2_mean", "mal_CLDN4_mean",
            "frac_tnk", "frac_cd8", "frac_b_plasma"]].to_string(index=False))

# --- primary: post-treatment, MPR vs NMPR, sample unit ---
post = samp.loc[samp.timing == "post"].copy()
post_ok = post.loc[post.mal_TACSTD2_mean.notna()].copy()
print(f"\npost-tx samples={len(post)} with mal>={MIN_MALIG}: {len(post_ok)}")
print(post_ok.mpr_group.value_counts(dropna=False).to_string())

for gene in ("TACSTD2", "CLDN4"):
    metric = f"mal_{gene}_mean"
    a = post_ok.loc[post_ok.mpr_group == "MPR", metric]
    b = post_ok.loc[post_ok.mpr_group == "NMPR", metric]
    u, p, na, nb = mwu(a, b)
    add_stat("GSE207422", "response_mpr_post", f"NMPR vs MPR (post, n_mal>={MIN_MALIG})",
             metric, f"{na}+{nb}", "mannwhitney_u", u, p,
             f"MPR median={a.median():.4f} (n={na}); NMPR median={b.median():.4f} (n={nb}); "
             f"primary unit=sample; pCR counted as MPR")
    print(f"  {gene} MPR n={na} med={a.median():.4f} | NMPR n={nb} med={b.median():.4f} | U={u} p={fmt_p(p)}")

    for immune, ilabel in (("frac_tnk", "T/NK fraction"),
                           ("frac_cd8", "CD8-like fraction"),
                           ("frac_b_plasma", "B/Plasma fraction"),
                           ("mean_tls_chemo", "TLS chemokine score")):
        r, p_s, n_s = spear(post_ok[metric], post_ok[immune])
        add_stat("GSE207422", "spearman_post", f"malignant {gene} vs {ilabel}",
                 metric, n_s, "spearman_rho", r, p_s,
                 "post-treatment samples with >=10 malignant-like cells")
        print(f"    vs {ilabel}: n={n_s} rho={fmt_r(r)} p={fmt_p(p_s)}")

# residual tumor (continuous pathologic burden) vs malignant TACSTD2
rt = pd.to_numeric(post_ok["Residual Tumor"], errors="coerce")
for gene in ("TACSTD2", "CLDN4"):
    r, p_s, n_s = spear(post_ok[f"mal_{gene}_mean"], rt)
    add_stat("GSE207422", "spearman_post", f"malignant {gene} vs residual tumor fraction",
             f"mal_{gene}_mean", n_s, "spearman_rho", r, p_s,
             "post-tx; residual tumor from GEO metadata")
    print(f"  {gene} vs residual tumor: n={n_s} rho={fmt_r(r)} p={fmt_p(p_s)}")

# immune fractions vs MPR (sanity / companion)
for immune, ilabel in (("frac_tnk", "T/NK fraction"),
                       ("frac_cd8", "CD8-like fraction"),
                       ("frac_b_plasma", "B/Plasma fraction")):
    a = post.loc[post.mpr_group == "MPR", immune]
    b = post.loc[post.mpr_group == "NMPR", immune]
    u, p, na, nb = mwu(a, b)
    add_stat("GSE207422", "immune_vs_mpr_post", f"NMPR vs MPR {ilabel}",
             immune, f"{na}+{nb}", "mannwhitney_u", u, p,
             f"MPR med={a.median():.4f}; NMPR med={b.median():.4f}")
    print(f"  {ilabel} MPR vs NMPR: U={u} p={fmt_p(p)}")

# exploratory: all samples including pre-tx
all_ok = samp.loc[samp.mal_TACSTD2_mean.notna() & samp.mpr_group.isin(["MPR", "NMPR"])]
for gene in ("TACSTD2", "CLDN4"):
    a = all_ok.loc[all_ok.mpr_group == "MPR", f"mal_{gene}_mean"]
    b = all_ok.loc[all_ok.mpr_group == "NMPR", f"mal_{gene}_mean"]
    u, p, na, nb = mwu(a, b)
    add_stat("GSE207422", "response_mpr_all_timing", "NMPR vs MPR (pre+post)",
             f"mal_{gene}_mean", f"{na}+{nb}", "mannwhitney_u", u, p,
             "exploratory: mixes pre-biopsy and post-surgery")

# paired malignant vs T/NK TACSTD2 (compartment specificity)
pair = post_ok.dropna(subset=["mal_TACSTD2_mean", "tnk_TACSTD2_mean"])
if len(pair) >= 4:
    w, p_w = stats.wilcoxon(pair["mal_TACSTD2_mean"], pair["tnk_TACSTD2_mean"])
    add_stat("GSE207422", "compartment", "malignant vs T/NK TACSTD2 (paired sample)",
             "mal_TACSTD2_mean", len(pair), "wilcoxon_w", float(w), float(p_w),
             f"mal med={pair.mal_TACSTD2_mean.median():.4f}; T/NK med={pair.tnk_TACSTD2_mean.median():.4f}")
    print(f"  paired mal vs T/NK TACSTD2: n={len(pair)} W={w} p={fmt_p(p_w)}")

# figures — GSE207422
fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0), sharey=False)
for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
    metric = f"mal_{gene}_mean"
    plot_df = post_ok.dropna(subset=[metric, "mpr_group"])
    groups = ["MPR", "NMPR"]
    data = [plot_df.loc[plot_df.mpr_group == g, metric].to_numpy() for g in groups]
    bp = ax.boxplot(data, labels=groups, patch_artist=True, widths=0.55)
    colors = ["#4C9F70", "#C44E52"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    for i, g in enumerate(groups, start=1):
        y = plot_df.loc[plot_df.mpr_group == g, metric]
        ax.scatter(np.random.default_rng(1).normal(i, 0.04, size=len(y)), y,
                   c="black", s=22, zorder=3)
    _, p, _, _ = mwu(data[0], data[1])
    annotate_box(ax, data[0], data[1], p)
    ax.set_title(f"malignant-like {gene}")
    ax.set_ylabel("mean log1p CP10K")
    ax.set_xlabel("pathologic response (post-tx)")
fig.suptitle("GSE207422 post-tx: malignant-like TACSTD2 / CLDN4 vs MPR", fontsize=11)
fig.tight_layout()
fig.savefig(RES / "gse207422_mal_tacstd2_cldn4_vs_mpr.png", dpi=160)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
for ax, immune, ilabel in zip(
    axes,
    ("frac_tnk", "frac_cd8", "frac_b_plasma"),
    ("T/NK fraction", "CD8-like fraction", "B/Plasma (TLS proxy)"),
):
    x = post_ok["mal_TACSTD2_mean"]
    y = post_ok[immune]
    color = post_ok["mpr_group"].map({"MPR": "#4C9F70", "NMPR": "#C44E52"}).fillna("#888")
    ax.scatter(x, y, c=color, s=36, edgecolor="k", linewidth=0.4)
    r, p_s, n_s = spear(x, y)
    ax.set_title(f"n={n_s}  ρ={fmt_r(r)}  p={fmt_p(p_s)}")
    ax.set_xlabel("malignant-like TACSTD2 (log1p CP10K)")
    ax.set_ylabel(ilabel)
fig.suptitle("GSE207422 post-tx: malignant TACSTD2 vs immune fractions (green=MPR, red=NMPR)", fontsize=10)
fig.tight_layout()
fig.savefig(RES / "gse207422_mal_tacstd2_vs_immune.png", dpi=160)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
for ax, immune, ilabel in zip(
    axes,
    ("frac_tnk", "frac_cd8", "frac_b_plasma"),
    ("T/NK fraction", "CD8-like fraction", "B/Plasma (TLS proxy)"),
):
    x = post_ok["mal_CLDN4_mean"]
    y = post_ok[immune]
    color = post_ok["mpr_group"].map({"MPR": "#4C9F70", "NMPR": "#C44E52"}).fillna("#888")
    ax.scatter(x, y, c=color, s=36, edgecolor="k", linewidth=0.4)
    r, p_s, n_s = spear(x, y)
    ax.set_title(f"n={n_s}  ρ={fmt_r(r)}  p={fmt_p(p_s)}")
    ax.set_xlabel("malignant-like CLDN4 (log1p CP10K)")
    ax.set_ylabel(ilabel)
fig.suptitle("GSE207422 post-tx: malignant CLDN4 vs immune fractions (green=MPR, red=NMPR)", fontsize=10)
fig.tight_layout()
fig.savefig(RES / "gse207422_mal_cldn4_vs_immune.png", dpi=160)
plt.close(fig)

# lineage bar
fig, ax = plt.subplots(figsize=(7.2, 3.6))
order = cells["lineage"].value_counts()
ax.bar(order.index.astype(str), order.values, color="#4C78A8")
ax.set_ylabel("cells (UMI≥200)")
ax.set_title("GSE207422 marker-based lineage (argmax score)")
ax.tick_params(axis="x", rotation=30)
fig.tight_layout()
fig.savefig(RES / "gse207422_lineage_counts.png", dpi=160)
plt.close(fig)

# compartment TACSTD2
fig, ax = plt.subplots(figsize=(5.2, 3.8))
comp = [
    cells.loc[cells.malignant_like, "TACSTD2_log1p"].dropna(),
    cells.loc[cells.is_tnk, "TACSTD2_log1p"].dropna(),
]
ax.boxplot(comp, labels=["malignant-like", "T/NK"], showfliers=False)
ax.set_ylabel("TACSTD2 log1p CP10K (cells)")
ax.set_title("GSE207422: TACSTD2 is epithelial-restricted\n(cell-level; exploratory)")
fig.tight_layout()
fig.savefig(RES / "gse207422_tacstd2_compartment.png", dpi=160)
plt.close(fig)

print("GSE207422 figures written", flush=True)


# ======================================================================
# GSE205335
# ======================================================================
print("\n=== GSE205335 ===", flush=True)
ident = pd.read_csv(DATA / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
panel335 = pd.read_csv(DATA / "gse205335_panel_cells.tsv.gz", sep="\t")
print(f"identity {ident.shape}  panel {panel335.shape}", flush=True)

# SOFT sample sheet
soft_recs = []
cur = None
with gzip.open(DATA / "GSE205335_family.soft.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("^SAMPLE = "):
            if cur:
                soft_recs.append(cur)
            cur = {"gsm": line.split(" = ", 1)[1].strip()}
        elif cur is not None:
            if line.startswith("!Sample_title = "):
                cur["title"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_description = ") and "description" not in cur:
                cur["description"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                v = line.split(" = ", 1)[1].strip()
                if ": " in v:
                    k, x = v.split(": ", 1)
                    cur[k] = x
    if cur:
        soft_recs.append(cur)
soft = pd.DataFrame(soft_recs)
soft["orig_key"] = soft["description"].str.replace("_", "-", regex=False)
# CellIdentity orig.ident = DESCRIPTION-with-dashes + -3P/-5P
ident["orig_key"] = ident["orig.ident"].str.replace(r"-[35]P$", "", regex=True)
merged_map = ident[["orig.ident", "orig_key"]].drop_duplicates().merge(
    soft, on="orig_key", how="left",
)
print("SOFT map unmatched orig.ident:",
      merged_map.loc[merged_map.gsm.isna(), "orig.ident"].tolist())
merged_map.to_csv(RES / "gse205335_sample_map.tsv", sep="\t", index=False)

cells335 = panel335.merge(ident, on="barcode", how="inner")
cells335 = cells335.merge(
    merged_map[["orig.ident", "gsm", "patient", "tissue", "cancer subtype",
                "tumor stage", "recist", "title"]],
    on="orig.ident", how="left",
)
cells335 = cells335.loc[cells335.total_umi >= MIN_UMI].copy()
print(f"merged cells {len(cells335)}  recist missing {cells335.recist.isna().sum()}")

count_genes335 = [c for c in panel335.columns if c not in ("barcode", "total_umi")]
ln335 = log1p_cp10k(cells335[present(count_genes335, cells335.columns)], cells335["total_umi"])
for gene in ("TACSTD2", "CLDN4"):
    if gene in ln335.columns:
        cells335[f"{gene}_log1p"] = ln335[gene]

cells335["is_mal"] = cells335["lineage.sub"] == "Malignant cells"
cells335["is_tnk"] = cells335["lineage.total"] == "T/NK cells"
cells335["is_cd8"] = cells335["lineage.sub"] == "CD8+ T cells"
cells335["is_b"] = cells335["lineage.total"] == "B/Plasma cells"
cells335["is_nk"] = cells335["lineage.sub"] == "NK cells"

# drop normal tissues from response / infiltration tests
NORMAL_TISSUE = {"Normal Brain", "Normal LN", "Normal Lung"}
cells335["is_normal_tissue"] = cells335["tissue"].isin(NORMAL_TISSUE)
tumor_cells = cells335.loc[~cells335.is_normal_tissue].copy()

print("lineage.sub (tumor tissues):\n", tumor_cells["lineage.sub"].value_counts(dropna=False).head(12).to_string())
print("recist patients:", tumor_cells.groupby("patient")["recist"].first().value_counts(dropna=False).to_dict())

pt_rows = []
for patient, g in tumor_cells.groupby("patient"):
    n = len(g)
    n_mal = int(g.is_mal.sum())
    n_tnk = int(g.is_tnk.sum())
    n_cd8 = int(g.is_cd8.sum())
    n_b = int(g.is_b.sum())
    rec = {
        "patient": patient,
        "recist": g["recist"].iloc[0],
        "cancer_subtype": g["cancer subtype"].iloc[0],
        "tissue": ",".join(sorted(g["tissue"].dropna().unique())),
        "n_cells": n,
        "n_malignant": n_mal,
        "n_tnk": n_tnk,
        "n_cd8": n_cd8,
        "n_b_plasma": n_b,
        "frac_tnk": n_tnk / n,
        "frac_cd8": n_cd8 / n,
        "frac_b_plasma": n_b / n,
    }
    rec["response"] = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR"}.get(rec["recist"], "NE")
    mal = g.loc[g.is_mal]
    tnk = g.loc[g.is_tnk]
    for gene in ("TACSTD2", "CLDN4"):
        col = f"{gene}_log1p"
        if col not in g.columns:
            continue
        if n_mal >= MIN_MALIG_PT:
            rec[f"mal_{gene}_mean"] = float(mal[col].mean())
            rec[f"mal_{gene}_pct_pos"] = float((mal[gene] > 0).mean() * 100)
        else:
            rec[f"mal_{gene}_mean"] = np.nan
            rec[f"mal_{gene}_pct_pos"] = np.nan
        rec[f"tnk_{gene}_mean"] = float(tnk[col].mean()) if n_tnk >= 10 else np.nan
    pt_rows.append(rec)

pt = pd.DataFrame(pt_rows)
pt.to_csv(RES / "gse205335_patient_table.tsv", sep="\t", index=False)
print(pt.to_string(index=False))

pt_ok = pt.loc[pt.mal_TACSTD2_mean.notna()].copy()
pt_rn = pt_ok.loc[pt_ok.response.isin(["R", "NR"])].copy()
print(f"\npatients with mal>={MIN_MALIG_PT}: {len(pt_ok)}  R/NR: {len(pt_rn)}")
print(pt_rn.response.value_counts().to_string())

for gene in ("TACSTD2", "CLDN4"):
    metric = f"mal_{gene}_mean"
    a = pt_rn.loc[pt_rn.response == "R", metric]
    b = pt_rn.loc[pt_rn.response == "NR", metric]
    u, p, na, nb = mwu(a, b)
    add_stat("GSE205335", "response_recist", "NR (SD/PD) vs R (PR)",
             metric, f"{na}+{nb}", "mannwhitney_u", u, p,
             f"R median={a.median():.4f} (n={na}); NR median={b.median():.4f} (n={nb}); "
             f"RECIST not MPR; unit=patient; min {MIN_MALIG_PT} malignant cells")
    print(f"  {gene} R n={na} med={a.median():.4f} | NR n={nb} med={b.median():.4f} | U={u} p={fmt_p(p)}")

    for immune, ilabel in (("frac_tnk", "T/NK fraction"),
                           ("frac_cd8", "CD8 fraction"),
                           ("frac_b_plasma", "B/Plasma fraction")):
        r, p_s, n_s = spear(pt_ok[metric], pt_ok[immune])
        add_stat("GSE205335", "spearman_all_evaluable", f"malignant {gene} vs {ilabel}",
                 metric, n_s, "spearman_rho", r, p_s,
                 "tumor tissues; patients with >=20 malignant cells; includes NE")
        print(f"    vs {ilabel} (all eval): n={n_s} rho={fmt_r(r)} p={fmt_p(p_s)}")
        r2, p2, n2 = spear(pt_rn[metric], pt_rn[immune])
        add_stat("GSE205335", "spearman_recist_r_nr", f"malignant {gene} vs {ilabel}",
                 metric, n2, "spearman_rho", r2, p2,
                 "R+NR only")
        print(f"    vs {ilabel} (R/NR): n={n2} rho={fmt_r(r2)} p={fmt_p(p2)}")

# NSCLC-only sensitivity (ADC+SQ)
nsclc = pt_rn.loc[pt_rn.cancer_subtype.isin(["ADC", "SQ"])]
print(f"NSCLC ADC+SQ R/NR n={len(nsclc)}")
for gene in ("TACSTD2", "CLDN4"):
    metric = f"mal_{gene}_mean"
    a = nsclc.loc[nsclc.response == "R", metric]
    b = nsclc.loc[nsclc.response == "NR", metric]
    u, p, na, nb = mwu(a, b)
    add_stat("GSE205335", "response_recist_nsclc", "NR vs R (ADC+SQ only)",
             metric, f"{na}+{nb}", "mannwhitney_u", u, p,
             "sensitivity; SCLC/NUT excluded")
    r, p_s, n_s = spear(nsclc[metric], nsclc["frac_tnk"])
    add_stat("GSE205335", "spearman_nsclc", f"malignant {gene} vs T/NK fraction",
             metric, n_s, "spearman_rho", r, p_s, "ADC+SQ R+NR")
    print(f"  NSCLC {gene} R vs NR p={fmt_p(p)}; vs T/NK n={n_s} rho={fmt_r(r)} p={fmt_p(p_s)}")

# paired compartment
pair335 = pt_ok.dropna(subset=["mal_TACSTD2_mean", "tnk_TACSTD2_mean"])
if len(pair335) >= 4:
    w, p_w = stats.wilcoxon(pair335["mal_TACSTD2_mean"], pair335["tnk_TACSTD2_mean"])
    add_stat("GSE205335", "compartment", "malignant vs T/NK TACSTD2 (paired patient)",
             "mal_TACSTD2_mean", len(pair335), "wilcoxon_w", float(w), float(p_w),
             f"mal med={pair335.mal_TACSTD2_mean.median():.4f}; T/NK med={pair335.tnk_TACSTD2_mean.median():.4f}")
    print(f"  paired mal vs T/NK TACSTD2: n={len(pair335)} W={w} p={fmt_p(p_w)}")

# figures GSE205335
fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0))
for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
    metric = f"mal_{gene}_mean"
    plot_df = pt_rn.dropna(subset=[metric])
    groups = ["R", "NR"]
    data = [plot_df.loc[plot_df.response == g, metric].to_numpy() for g in groups]
    bp = ax.boxplot(data, labels=groups, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], ["#4C9F70", "#C44E52"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(2)
    for i, g in enumerate(groups, start=1):
        y = plot_df.loc[plot_df.response == g, metric]
        ax.scatter(rng.normal(i, 0.04, size=len(y)), y, c="black", s=22, zorder=3)
    _, p, _, _ = mwu(data[0], data[1])
    annotate_box(ax, data[0], data[1], p)
    ax.set_title(f"malignant {gene}")
    ax.set_ylabel("mean log1p CP10K")
    ax.set_xlabel("RECIST (PR=R, SD/PD=NR)")
fig.suptitle("GSE205335: malignant TACSTD2 / CLDN4 vs RECIST (not MPR)", fontsize=11)
fig.tight_layout()
fig.savefig(RES / "gse205335_mal_tacstd2_cldn4_vs_recist.png", dpi=160)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
for ax, immune, ilabel in zip(
    axes,
    ("frac_tnk", "frac_cd8", "frac_b_plasma"),
    ("T/NK fraction", "CD8 fraction", "B/Plasma fraction"),
):
    x = pt_ok["mal_TACSTD2_mean"]
    y = pt_ok[immune]
    color = pt_ok["response"].map({"R": "#4C9F70", "NR": "#C44E52", "NE": "#888888"})
    ax.scatter(x, y, c=color, s=36, edgecolor="k", linewidth=0.4)
    r, p_s, n_s = spear(x, y)
    ax.set_title(f"n={n_s}  ρ={fmt_r(r)}  p={fmt_p(p_s)}")
    ax.set_xlabel("malignant TACSTD2 (log1p CP10K)")
    ax.set_ylabel(ilabel)
fig.suptitle("GSE205335: malignant TACSTD2 vs immune fractions (green=PR, red=SD/PD, grey=NE)", fontsize=10)
fig.tight_layout()
fig.savefig(RES / "gse205335_mal_tacstd2_vs_immune.png", dpi=160)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
for ax, immune, ilabel in zip(
    axes,
    ("frac_tnk", "frac_cd8", "frac_b_plasma"),
    ("T/NK fraction", "CD8 fraction", "B/Plasma fraction"),
):
    x = pt_ok["mal_CLDN4_mean"]
    y = pt_ok[immune]
    color = pt_ok["response"].map({"R": "#4C9F70", "NR": "#C44E52", "NE": "#888888"})
    ax.scatter(x, y, c=color, s=36, edgecolor="k", linewidth=0.4)
    r, p_s, n_s = spear(x, y)
    ax.set_title(f"n={n_s}  ρ={fmt_r(r)}  p={fmt_p(p_s)}")
    ax.set_xlabel("malignant CLDN4 (log1p CP10K)")
    ax.set_ylabel(ilabel)
fig.suptitle("GSE205335: malignant CLDN4 vs immune fractions", fontsize=10)
fig.tight_layout()
fig.savefig(RES / "gse205335_mal_cldn4_vs_immune.png", dpi=160)
plt.close(fig)

print("GSE205335 figures written", flush=True)

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(RES / "stats.tsv", sep="\t", index=False)
print("\n=== STATS ===")
print(stats_df.to_string(index=False))
print(f"\nwrote {RES / 'stats.tsv'}")
