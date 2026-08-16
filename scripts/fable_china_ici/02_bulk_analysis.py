#!/usr/bin/env python3
"""TACSTD2 (TROP2) / CLDN4 vs ICI response in open East-Asian NSCLC bulk expression sets.

Datasets (downloaded by 01_download_data.sh):
  - GSE207422  Shanghai Pulmonary Hospital (CN), neoadjuvant anti-PD-1 + chemo,
               tumor bulk RNA-seq (log2TPM), pathologic response MPR/NMPR + RECIST.
               Primary analysis restricted to pre-treatment biopsies.
  - GSE126044  Yonsei (KR), anti-PD-1 monotherapy, tumor RNA-seq counts, R vs NR.
  - GSE135222  KR, anti-PD-1/PD-L1, tumor RNA-seq (normalized), PFS event/time.
               DCB = PFS >= 180 days; also KM/Cox on median split.
  - GSE260770  Guangzhou Med Univ (CN), sintilimab, plasma exosomal mRNA FPKM,
               R vs NR (EXPLORATORY: blood exosome, not tumor tissue).

Outputs -> results/fable_china_ici/{tables,figures}
"""

import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.environ.get("DATA_DIR", "/tmp/fable_china_ici_data")
ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
TABLES = os.path.join(ROOT, "results", "fable_china_ici", "tables")
FIGS = os.path.join(ROOT, "results", "fable_china_ici", "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
ENSEMBL = {"ENSG00000184292": "TACSTD2", "ENSG00000189143": "CLDN4"}


def parse_series_matrix(path):
    """Return DataFrame: rows = samples, columns = title/geo + characteristics."""
    rows = {}
    chars = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                rows["title"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                rows["gsm"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":", 1)[0].strip()
                chars.append((key, [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]))
    df = pd.DataFrame(rows)
    for key, vals in chars:
        df[key] = vals
    return df


def mw_stats(dataset, gene, expr, group, pos_label, neg_label, unit, note=""):
    """Mann-Whitney U + AUC for expr in pos vs neg group."""
    x = np.asarray(expr, dtype=float)
    g = np.asarray(group)
    pos = x[g == pos_label]
    neg = x[g == neg_label]
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = u / (len(pos) * len(neg))  # AUC for "higher expr predicts pos_label"
    return dict(
        dataset=dataset, gene=gene, unit=unit,
        group_pos=pos_label, n_pos=len(pos), median_pos=round(float(np.median(pos)), 3),
        group_neg=neg_label, n_neg=len(neg), median_neg=round(float(np.median(neg)), 3),
        mannwhitney_p=round(float(p), 4), auc_pos_high=round(float(auc), 3), note=note,
    )


stats_rows = []
long_rows = []


def add_long(dataset, df, group_col):
    for _, r in df.iterrows():
        long_rows.append(dict(dataset=dataset, sample=r["sample"], group=r[group_col],
                              TACSTD2=r["TACSTD2"], CLDN4=r["CLDN4"]))


# ------------------------------------------------------------------ GSE207422
expr = pd.read_csv(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
meta = pd.read_excel(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
meta = meta.dropna(subset=["Patient"])
meta["mpr"] = np.where(meta["Pathologic Response"].str.startswith(("MPR", "pCR")), "MPR", "NMPR")
sub = expr.loc[GENES].T
sub.index.name = "sample"
sub = sub.reset_index()
m = sub.merge(meta, left_on="sample", right_on="Sample", how="inner")
pre = m[m["Resource"].str.contains("Pre", case=False)].copy()
print(f"GSE207422: matrix n={expr.shape[1]}, with metadata n={len(m)}, pre-treatment n={len(pre)} "
      f"(MPR {sum(pre['mpr'] == 'MPR')} / NMPR {sum(pre['mpr'] == 'NMPR')})")
for gene in GENES:
    stats_rows.append(mw_stats("GSE207422 (CN, neoadj PD-1+chemo, pre-tx tumor)", gene,
                               pre[gene], pre["mpr"], "MPR", "NMPR", "log2TPM",
                               "pathologic response, pre-treatment biopsies only"))
add_long("GSE207422_pre", pre.assign(group=pre["mpr"]), "group")
gse207422_pre = pre

# ------------------------------------------------------------------ GSE126044
cnt = pd.read_csv(f"{DATA}/GSE126044_counts.txt.gz", sep="\t", index_col=0)
cpm = np.log2(cnt / cnt.sum(axis=0) * 1e6 + 1)
sm = parse_series_matrix(f"{DATA}/GSE126044_series_matrix.txt.gz")
sm["sample"] = sm["title"].str.replace("RNA-seq_", "", regex=False)
sub = cpm.loc[GENES].T
sub.index.name = "sample"
sub = sub.reset_index().merge(sm[["sample", "patient response"]], on="sample")
sub["resp"] = sub["patient response"].map({"responder": "R", "non-responder": "NR"})
print(f"GSE126044: n={len(sub)} (R {sum(sub['resp'] == 'R')} / NR {sum(sub['resp'] == 'NR')})")
for gene in GENES:
    stats_rows.append(mw_stats("GSE126044 (KR, anti-PD-1, tumor)", gene,
                               sub[gene], sub["resp"], "R", "NR", "log2CPM", "RECIST responder vs non-responder"))
add_long("GSE126044", sub.assign(group=sub["resp"]), "group")
gse126044 = sub

# ------------------------------------------------------------------ GSE135222
expr = pd.read_csv(f"{DATA}/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0)
expr.index = expr.index.str.split(".").str[0]
sub = np.log2(expr.loc[list(ENSEMBL)].rename(index=ENSEMBL).T + 1)
sub.index = sub.index.str.replace("NSCLC", "NSCLC ", regex=False)
sub.index.name = "title"
sub = sub.reset_index()
sm = parse_series_matrix(f"{DATA}/GSE135222_series_matrix.txt.gz")
sm = sm.rename(columns={"progression-free survival (pfs)": "pfs_event"})
sub = sub.merge(sm[["title", "pfs_event", "pfs.time"]], on="title")
sub["pfs_event"] = sub["pfs_event"].astype(int)
sub["pfs_days"] = sub["pfs.time"].astype(float)
sub["dcb"] = np.where(sub["pfs_days"] >= 180, "DCB", "NDB")  # all PFS<180 here are events
sub["sample"] = sub["title"]
print(f"GSE135222: n={len(sub)} (DCB {sum(sub['dcb'] == 'DCB')} / NDB {sum(sub['dcb'] == 'NDB')})")
for gene in GENES:
    stats_rows.append(mw_stats("GSE135222 (KR, anti-PD-1/PD-L1, tumor)", gene,
                               sub[gene], sub["dcb"], "DCB", "NDB", "log2(norm+1)",
                               "DCB = PFS >= 180 days"))
add_long("GSE135222", sub.assign(group=sub["dcb"]), "group")

# KM / Cox on median split
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
for ax, gene in zip(axes, GENES):
    hi = sub[gene] >= sub[gene].median()
    lr = logrank_test(sub.loc[hi, "pfs_days"], sub.loc[~hi, "pfs_days"],
                      sub.loc[hi, "pfs_event"], sub.loc[~hi, "pfs_event"])
    cph = CoxPHFitter().fit(sub[[gene, "pfs_days", "pfs_event"]], "pfs_days", "pfs_event")
    hr = float(np.exp(cph.params_[gene]))
    cox_p = float(cph.summary.loc[gene, "p"])
    stats_rows.append(dict(
        dataset="GSE135222 (KR, anti-PD-1/PD-L1, tumor)", gene=gene, unit="log2(norm+1)",
        group_pos=f"{gene} high (>=median)", n_pos=int(hi.sum()),
        median_pos=np.nan, group_neg=f"{gene} low", n_neg=int((~hi).sum()), median_neg=np.nan,
        mannwhitney_p=np.nan, auc_pos_high=np.nan,
        note=f"PFS logrank p={lr.p_value:.3f}; Cox per-unit HR={hr:.2f} (p={cox_p:.3f})"))
    for mask, lab, color in [(hi, "high", "#c0392b"), (~hi, "low", "#2980b9")]:
        km = KaplanMeierFitter().fit(sub.loc[mask, "pfs_days"], sub.loc[mask, "pfs_event"],
                                     label=f"{gene} {lab} (n={mask.sum()})")
        km.plot_survival_function(ax=ax, color=color, ci_show=False)
    ax.set_title(f"GSE135222 PFS by {gene} (logrank p={lr.p_value:.3f})", fontsize=10)
    ax.set_xlabel("Days")
    ax.set_ylabel("PFS probability")
fig.tight_layout()
fig.savefig(f"{FIGS}/gse135222_km_pfs.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ GSE260770 (exploratory)
expr = pd.read_csv(f"{DATA}/GSE260770_mRNA_FPKM.txt.gz", sep="\t")
expr = expr.set_index("Symbol").drop(columns=["#ID"])
sub = np.log2(expr.loc[GENES].T + 1)
sub.index.name = "sample"
sub = sub.reset_index()
sm = parse_series_matrix(f"{DATA}/GSE260770_series_matrix.txt.gz")
sm["resp"] = np.where(sm["group"].str.startswith("Responsed"), "R", "NR")
sub = sub.merge(sm[["title", "resp"]], left_on="sample", right_on="title")
print(f"GSE260770: n={len(sub)} (R {sum(sub['resp'] == 'R')} / NR {sum(sub['resp'] == 'NR')})")
for gene in GENES:
    stats_rows.append(mw_stats("GSE260770 (CN, sintilimab GGO, plasma exosome; EXPLORATORY)", gene,
                               sub[gene], sub["resp"], "R", "NR", "log2(FPKM+1)",
                               "plasma exosomal mRNA, not tumor tissue"))
add_long("GSE260770_exosome", sub.assign(group=sub["resp"]), "group")

# ------------------------------------------------------------------ outputs
stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(f"{TABLES}/bulk_stats.csv", index=False)
pd.DataFrame(long_rows).to_csv(f"{TABLES}/bulk_expression_by_sample.csv", index=False)
print("\n== bulk_stats.csv ==")
print(stats_df.to_string(index=False))

# boxplot grid
long_df = pd.DataFrame(long_rows)
panel_defs = [
    ("GSE207422_pre", "GSE207422 (CN)\npre-tx tumor, MPR vs NMPR", ["MPR", "NMPR"]),
    ("GSE126044", "GSE126044 (KR)\ntumor, R vs NR", ["R", "NR"]),
    ("GSE135222", "GSE135222 (KR)\ntumor, DCB vs NDB", ["DCB", "NDB"]),
    ("GSE260770_exosome", "GSE260770 (CN)\nplasma exosome, R vs NR", ["R", "NR"]),
]
fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharey=False)
for j, (ds, label, order) in enumerate(panel_defs):
    d = long_df[long_df.dataset == ds]
    for i, gene in enumerate(GENES):
        ax = axes[i, j]
        vals = [d.loc[d.group == g, gene].astype(float).values for g in order]
        bp = ax.boxplot(vals, tick_labels=[f"{g}\n(n={len(v)})" for g, v in zip(order, vals)],
                        widths=0.55, showfliers=False, patch_artist=True)
        for patch, color in zip(bp["boxes"], ["#f6b0a0", "#a8c8e8"]):
            patch.set_facecolor(color)
        for k, v in enumerate(vals):
            ax.scatter(np.random.default_rng(1).normal(k + 1, 0.06, len(v)), v,
                       s=18, color="#333333", alpha=0.75, zorder=3)
        row = stats_df[(stats_df.dataset.str.startswith(ds.split("_")[0])) & (stats_df.gene == gene)
                       & stats_df.mannwhitney_p.notna()]
        p = row.mannwhitney_p.iloc[0] if len(row) else np.nan
        ax.set_title(f"{gene}  (MW p={p:.3f})", fontsize=10)
        if i == 0:
            ax.text(0.5, 1.22, label, ha="center", va="bottom", transform=ax.transAxes, fontsize=9)
        ax.set_ylabel("expression" if j == 0 else "")
fig.suptitle("TACSTD2 / CLDN4 vs ICI response — open East-Asian NSCLC cohorts", y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{FIGS}/bulk_boxplots.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nFigures -> {FIGS}")
