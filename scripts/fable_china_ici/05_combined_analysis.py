#!/usr/bin/env python3
"""Combined TACSTD2/CLDN4 analysis on open processed matrices.

Adds:
  - GSE329813 (PUMC, CN): GeoMx DSP, neoadjuvant pembrolizumab + chemo, 22 pts,
    TACSTD2 on panel (CLDN4 not on GeoMx IO panel). Patient-mean of primary
    tumor-bed ROIs vs MPR/NMPR; ROI-level as sensitivity.
  - GSE207422 RECIST (CR/PR vs SD) on the same pre-treatment bulk biopsies.
  - Stouffer meta-analysis + forest plot of responder-high AUC across
    pre-treatment tumor cohorts (GSE207422 MPR, GSE126044, GSE135222, GSE329813).

Reads existing tables from 02_bulk_analysis.py when present; recomputes
GSE207422 RECIST and GSE329813 from cached downloads.

Outputs -> results/fable_china_ici/{tables,figures}
"""

import gzip
import os
import re

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


def parse_series_matrix(path):
    rows, chars = {}, []
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
    x = np.asarray(expr, dtype=float)
    g = np.asarray(group)
    pos, neg = x[g == pos_label], x[g == neg_label]
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = u / (len(pos) * len(neg))
    return dict(
        dataset=dataset, gene=gene, unit=unit,
        group_pos=pos_label, n_pos=len(pos), median_pos=round(float(np.median(pos)), 3),
        group_neg=neg_label, n_neg=len(neg), median_neg=round(float(np.median(neg)), 3),
        mannwhitney_p=round(float(p), 4), auc_pos_high=round(float(auc), 3),
        auc_raw=float(auc), p_raw=float(p), n_total=len(pos) + len(neg), note=note,
    )


new_rows = []

# ------------------------------------------------------------------ GSE207422 RECIST
expr = pd.read_csv(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
meta = pd.read_excel(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx").dropna(subset=["Patient"])
sub = expr.loc[["TACSTD2", "CLDN4"]].T.reset_index().rename(columns={"index": "sample"})
m = sub.merge(meta, left_on="sample", right_on="Sample")
m["recist_bin"] = np.where(m["RECIST"].isin(["CR", "PR"]), "CR/PR", "SD")
print(f"GSE207422 RECIST: CR/PR {sum(m.recist_bin == 'CR/PR')} / SD {sum(m.recist_bin == 'SD')}")
for gene in ["TACSTD2", "CLDN4"]:
    new_rows.append(mw_stats(
        "GSE207422 (CN, neoadj PD-1+chemo, pre-tx tumor, RECIST)", gene,
        m[gene], m["recist_bin"], "CR/PR", "SD", "log2TPM",
        "RECIST CR/PR vs SD on same 24 pre-treatment biopsies"))

# ------------------------------------------------------------------ GSE329813 GeoMx
expr = pd.read_csv(f"{DATA}/GSE329813_processed_data_file_normalized_data.csv.gz")
expr = expr.set_index("Row")
sm = parse_series_matrix(f"{DATA}/GSE329813_series_matrix.txt.gz")
parsed = []
for t in sm["title"]:
    mm = re.match(r"ROI (\d+), Patient (\d+), (.+), (MPR|NMPR)", t)
    if not mm:
        raise ValueError(f"unparsed title: {t}")
    parsed.append(dict(roi=f"ROI {mm.group(1)}", patient=int(mm.group(2)),
                       tissue=mm.group(3), resp=mm.group(4)))
roi_meta = pd.DataFrame(parsed)
# matrix columns are 'ROI 31' etc.
common = [c for c in expr.columns if c in set(roi_meta.roi)]
print(f"GSE329813: {len(common)} ROIs matched, genes={expr.shape[0]}, TACSTD2={'TACSTD2' in expr.index}, CLDN4={'CLDN4' in expr.index}")
roi = roi_meta[roi_meta.roi.isin(common)].copy()
roi["TACSTD2"] = expr.loc["TACSTD2", roi.roi].astype(float).values
tumor = roi[roi.tissue == "Primary tumor bed"]
pt = tumor.groupby("patient").agg(TACSTD2=("TACSTD2", "mean"), resp=("resp", "first"), n_roi=("TACSTD2", "size")).reset_index()
print(f"  tumor-bed patients: {len(pt)} (MPR {sum(pt.resp == 'MPR')} / NMPR {sum(pt.resp == 'NMPR')})")
new_rows.append(mw_stats(
    "GSE329813 (CN, PUMC, neoadj pembro+chemo, GeoMx tumor-bed, patient-mean)",
    "TACSTD2", pt["TACSTD2"], pt["resp"], "MPR", "NMPR", "GeoMx Q3-norm",
    "patient mean of primary tumor-bed ROIs; CLDN4 not on GeoMx panel"))
new_rows.append(mw_stats(
    "GSE329813 (CN, PUMC, GeoMx tumor-bed ROI-level; sensitivity)",
    "TACSTD2", tumor["TACSTD2"], tumor["resp"], "MPR", "NMPR", "GeoMx Q3-norm",
    "ROI-level (not independent); CLDN4 not on panel"))
pt.to_csv(f"{TABLES}/gse329813_tacstd2_by_patient.csv", index=False)
roi.to_csv(f"{TABLES}/gse329813_tacstd2_by_roi.csv", index=False)

# ------------------------------------------------------------------ merge with existing bulk_stats
exist = pd.read_csv(f"{TABLES}/bulk_stats.csv")
# drop leftover raw cols if re-run
for c in ["auc_raw", "p_raw", "n_total"]:
    if c in exist.columns:
        exist = exist.drop(columns=c)
add = pd.DataFrame(new_rows)
# persist new rows without raw helper cols into an extended stats file
keep_cols = ["dataset", "gene", "unit", "group_pos", "n_pos", "median_pos",
             "group_neg", "n_neg", "median_neg", "mannwhitney_p", "auc_pos_high", "note"]
ext = pd.concat([exist[keep_cols], add[keep_cols]], ignore_index=True)
ext.to_csv(f"{TABLES}/bulk_stats_extended.csv", index=False)

# ------------------------------------------------------------------ Stouffer meta
# Pre-tx tumor (independent patients): GSE207422 MPR, GSE126044, GSE135222 DCB
# GSE329813 is POST-resection GeoMx — shown on the forest plot but excluded from
# the pre-tx Stouffer (residual-tumor content confounds TACSTD2 after MPR).
meta_keys = [
    ("GSE207422 (CN, neoadj PD-1+chemo, pre-tx tumor)", "TACSTD2"),
    ("GSE207422 (CN, neoadj PD-1+chemo, pre-tx tumor)", "CLDN4"),
    ("GSE126044 (KR, anti-PD-1, tumor)", "TACSTD2"),
    ("GSE126044 (KR, anti-PD-1, tumor)", "CLDN4"),
    ("GSE135222 (KR, anti-PD-1/PD-L1, tumor)", "TACSTD2"),
    ("GSE135222 (KR, anti-PD-1/PD-L1, tumor)", "CLDN4"),
]
post_keys = [
    ("GSE329813 (CN, PUMC, neoadj pembro+chemo, GeoMx tumor-bed, patient-mean)", "TACSTD2"),
]
# attach raw p/auc from exist + add
all_for_meta = []
for _, r in exist.iterrows():
    if pd.isna(r.get("mannwhitney_p")):
        continue
    all_for_meta.append(r)
for r in new_rows:
    all_for_meta.append(pd.Series(r))
lookup = {(r["dataset"], r["gene"]): r for r in all_for_meta}

def _signed(r):
    p = float(r["p_raw"] if "p_raw" in r and pd.notna(r.get("p_raw")) else r["mannwhitney_p"])
    auc = float(r["auc_raw"] if "auc_raw" in r and pd.notna(r.get("auc_raw")) else r["auc_pos_high"])
    z = stats.norm.ppf(1 - max(p, 1e-15) / 2.0) * np.sign(auc - 0.5)
    n = int(r["n_pos"] + r["n_neg"]) if "n_total" not in r or pd.isna(r.get("n_total")) else int(r["n_total"])
    return p, auc, float(z), n


meta_out = []
for gene in ["TACSTD2", "CLDN4"]:
    zs, ns = [], []
    for ds, g in meta_keys:
        if g != gene:
            continue
        p, auc, z, n = _signed(lookup[(ds, g)])
        zs.append(z)
        ns.append(n)
        meta_out.append(dict(gene=gene, dataset=ds, n=n, auc=round(auc, 3), p=round(p, 4),
                             signed_z=round(z, 3), stratum="pre-tx tumor"))
    Z = float(np.sum(zs) / np.sqrt(len(zs)))
    p_meta = float(2 * (1 - stats.norm.cdf(abs(Z))))
    meta_out.append(dict(gene=gene, dataset="STOUFFER_META (pre-tx tumor cohorts)",
                         n=int(np.sum(ns)), auc=np.nan, p=round(p_meta, 4),
                         signed_z=round(Z, 3),
                         stratum=f"k={len(zs)} independent pre-tx cohorts; +z = higher in responders"))
    print(f"Stouffer pre-tx {gene}: Z={Z:.3f} p={p_meta:.4f} k={len(zs)} N={sum(ns)}")
for ds, g in post_keys:
    p, auc, z, n = _signed(lookup[(ds, g)])
    meta_out.append(dict(gene=g, dataset=ds, n=n, auc=round(auc, 3), p=round(p, 4),
                         signed_z=round(z, 3), stratum="post-tx GeoMx (residual-tumor confound)"))

pd.DataFrame(meta_out).to_csv(f"{TABLES}/stouffer_meta.csv", index=False)

# ------------------------------------------------------------------ forest plot
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharex=True)
labels_short = {
    "GSE207422 (CN, neoadj PD-1+chemo, pre-tx tumor)": "GSE207422 CN pre-tx bulk",
    "GSE126044 (KR, anti-PD-1, tumor)": "GSE126044 KR pre-tx bulk",
    "GSE135222 (KR, anti-PD-1/PD-L1, tumor)": "GSE135222 KR pre-tx bulk",
    "GSE329813 (CN, PUMC, neoadj pembro+chemo, GeoMx tumor-bed, patient-mean)": "GSE329813 CN post-tx GeoMx*",
}
for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
    rows = [r for r in meta_out if r["gene"] == gene and not r["dataset"].startswith("STOUFFER")]
    y = np.arange(len(rows))[::-1]
    aucs = [r["auc"] for r in rows]
    names = [f"{labels_short.get(r['dataset'], r['dataset'][:28])}  n={r['n']}" for r in rows]
    ax.axvline(0.5, color="#888", ls="--", lw=1)
    for yi, a, r in zip(y, aucs, rows):
        col = "#8e44ad" if "post-tx" in r.get("stratum", "") else "#c0392b"
        ax.plot([0.5, a], [yi, yi], color=col, lw=2)
        ax.scatter([a], [yi], s=70, color=col, zorder=3)
        ax.text(min(a + 0.015, 0.72), yi, f"AUC={a:.2f} p={r['p']:.3f}", va="center", fontsize=8)
    meta_r = [r for r in meta_out if r["gene"] == gene and r["dataset"].startswith("STOUFFER")]
    title = f"{gene}"
    if meta_r:
        title += f"   pre-tx Stouffer Z={meta_r[0]['signed_z']:.2f} p={meta_r[0]['p']:.3f}"
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlim(0.05, 0.90)
    ax.set_xlabel("AUC (higher expression → responder)")
    ax.set_title(title, fontsize=10)
fig.suptitle("TACSTD2/CLDN4 vs ICI response — open East-Asian tumor cohorts (*post-tx residual-tumor confound)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIGS}/forest_auc_pretreatment.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ GSE329813 boxplot
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8))
for ax, (lab, d, gcol) in zip(axes, [
    ("patient-mean tumor-bed", pt, "resp"),
    ("ROI-level tumor-bed", tumor, "resp"),
]):
    order = ["MPR", "NMPR"]
    vals = [d.loc[d[gcol] == g, "TACSTD2"].astype(float).values for g in order]
    bp = ax.boxplot(vals, tick_labels=[f"{g}\n(n={len(v)})" for g, v in zip(order, vals)],
                    widths=0.55, showfliers=False, patch_artist=True)
    for patch, color in zip(bp["boxes"], ["#f6b0a0", "#a8c8e8"]):
        patch.set_facecolor(color)
    rng = np.random.default_rng(4)
    for k, v in enumerate(vals):
        ax.scatter(rng.normal(k + 1, 0.06, len(v)), v, s=18, color="#333", alpha=0.75, zorder=3)
    u, p = stats.mannwhitneyu(vals[0], vals[1], alternative="two-sided")
    ax.set_title(f"TACSTD2 {lab}\nMW p={p:.3f}", fontsize=10)
    ax.set_ylabel("GeoMx Q3-normalized")
fig.suptitle("GSE329813 (PUMC): neoadjuvant pembrolizumab + chemo, GeoMx DSP", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIGS}/gse329813_tacstd2.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ GSE207422 RECIST box
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6))
for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
    order = ["CR/PR", "SD"]
    vals = [m.loc[m.recist_bin == g, gene].astype(float).values for g in order]
    bp = ax.boxplot(vals, tick_labels=[f"{g}\n(n={len(v)})" for g, v in zip(order, vals)],
                    widths=0.55, showfliers=False, patch_artist=True)
    for patch, color in zip(bp["boxes"], ["#f6b0a0", "#a8c8e8"]):
        patch.set_facecolor(color)
    rng = np.random.default_rng(5)
    for k, v in enumerate(vals):
        ax.scatter(rng.normal(k + 1, 0.06, len(v)), v, s=18, color="#333", alpha=0.75, zorder=3)
    u, p = stats.mannwhitneyu(vals[0], vals[1], alternative="two-sided")
    ax.set_title(f"{gene} (MW p={p:.3f})", fontsize=10)
    ax.set_ylabel("log2TPM")
fig.suptitle("GSE207422 pre-tx bulk: TACSTD2/CLDN4 vs RECIST (CR/PR vs SD)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIGS}/gse207422_recist.png", dpi=200)
plt.close(fig)

print("\n== new/extended stats ==")
print(add[keep_cols].to_string(index=False))
print(f"\nFigures -> {FIGS}")
