#!/usr/bin/env python3
"""fable_adc slice — expression vs outcome analysis for TROP2/CLDN ADC targets in lung.

Uses three OPEN GEO NSCLC cohorts (ICI or chemo-ICI) as public surrogates to relate the
ADC-relevant genes to immunotherapy outcome:

    TACSTD2 (TROP2)  -> SG / Dato-DXd target
    CLDN18           -> Claudin-18 (18.2 is therapeutic epitope; bulk cannot resolve isoform)
    TOP1             -> SG payload (SN-38) / Dato-DXd payload (DXd) target
    CD274 (PD-L1)    -> ICI context / SG+pembro combination rationale

Cohorts:
    GSE126044  anti-PD-1 mono, responder vs non-responder  (raw counts, symbols)
    GSE135222  anti-PD-1/PD-L1, PFS event + time           (TPM, Ensembl)
    GSE207422  neoadjuvant chemo+anti-PD-1, MPR vs NMPR    (log2TPM, symbols)

No TROP2-ADC patient-level omics are public; findings here describe target biology & ICI outcome,
not ADC response. Outputs -> results/fable_adc/.
"""
from __future__ import annotations
import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results" / "fable_adc" / "raw"
OUT = ROOT / "results" / "fable_adc"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

GENES = {
    "TACSTD2": "ENSG00000184292",
    "CLDN18": "ENSG00000066405",
    "TOP1": "ENSG00000198900",
    "CD274": "ENSG00000120217",
}
GENE_LABEL = {
    "TACSTD2": "TACSTD2 (TROP2)",
    "CLDN18": "CLDN18 (Claudin-18)",
    "TOP1": "TOP1 (SN-38/DXd target)",
    "CD274": "CD274 (PD-L1)",
}

rng = np.random.default_rng(0)


# --------------------------------------------------------------------------- helpers
def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Return a dataframe: rows=samples, columns=characteristic key/value fields."""
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    rec: dict[str, list[str]] = {}
    geo = titles = None
    char_rows = []
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            geo = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            char_rows.append([x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
    n = len(geo)
    df = pd.DataFrame({"gsm": geo, "title": titles})
    for row in char_rows:
        # each entry looks like "key: value"; key can repeat -> disambiguate
        keys = [c.split(":", 1)[0].strip() if ":" in c else "field" for c in row]
        key = keys[0] if keys else "field"
        base = key
        i = 2
        while base in df.columns:
            base = f"{key}_{i}"
            i += 1
        df[base] = [c.split(":", 1)[1].strip() if ":" in c else c for c in row]
    return df


def mannwhitney(group_a, group_b):
    a = np.asarray(group_a, float)
    b = np.asarray(group_b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # rank-biserial effect size
    rbc = 1 - (2 * u) / (len(a) * len(b))
    return p, rbc


def save_box(data_by_group, title, ylabel, fname, annotate=None):
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    labels = list(data_by_group.keys())
    vals = [np.asarray(data_by_group[k], float) for k in labels]
    bp = ax.boxplot(vals, tick_labels=labels, showfliers=False, widths=0.55, patch_artist=True)
    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    for i, v in enumerate(vals, start=1):
        jitter = rng.normal(0, 0.05, size=len(v))
        ax.scatter(np.full(len(v), i) + jitter, v, s=18, color=colors[(i - 1) % len(colors)],
                   edgecolor="black", linewidth=0.3, zorder=3)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)
    if annotate:
        ax.text(0.5, 0.98, annotate, transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- cohorts
def load_gse126044():
    counts = pd.read_csv(RAW / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    # CPM log2 normalization
    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    meta = parse_series_matrix(RAW / "GSE126044_series_matrix.txt.gz")
    # map title RNA-seq_Dis_XX -> counts column Dis_XX
    meta["col"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    resp_col = [c for c in meta.columns if c.lower().startswith("patient response")]
    meta["response"] = meta[resp_col[0]] if resp_col else meta.filter(like="response").iloc[:, 0]
    meta["response"] = meta["response"].str.strip()
    meta = meta[meta["col"].isin(logcpm.columns)].copy()
    return logcpm, meta


def load_gse135222():
    exp = pd.read_csv(RAW / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0)
    exp.index = [i.split(".")[0] for i in exp.index]  # strip Ensembl version
    log = np.log2(exp + 1)
    meta = parse_series_matrix(RAW / "GSE135222_series_matrix.txt.gz")
    meta["key"] = meta["title"].str.extract(r"(\d+)")[0]
    # expression columns like NSCLC378
    col_key = {re.sub(r"\D", "", c): c for c in exp.columns}
    meta["col"] = meta["key"].map(col_key)
    pfs_evt = [c for c in meta.columns if c.lower().startswith("progression-free survival")]
    pfs_t = [c for c in meta.columns if c.lower().startswith("pfs.time")]
    meta["pfs_event"] = pd.to_numeric(meta[pfs_evt[0]], errors="coerce")
    meta["pfs_time"] = pd.to_numeric(meta[pfs_t[0]], errors="coerce")
    meta = meta[meta["col"].notna()].copy()
    return log, meta


def load_gse207422():
    exp = pd.read_csv(RAW / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
    meta = pd.read_excel(RAW / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    meta.columns = [str(c).strip() for c in meta.columns]
    meta = meta[meta["Sample"].isin(exp.columns)].copy()
    meta["mpr"] = np.where(meta["Pathologic Response"].str.contains("MPR", case=False, na=False)
                           & ~meta["Pathologic Response"].str.contains("NMPR", case=False, na=False),
                           "MPR", "NMPR")
    meta["residual"] = pd.to_numeric(meta["Residual Tumor"], errors="coerce")
    return exp, meta


# --------------------------------------------------------------------------- analyses
results = {"genes": GENE_LABEL, "cohorts": {}}
diff_rows = []
corr_rows = []
surv_rows = []


def gene_vec_symbol(mat, gene, cols):
    if gene not in mat.index:
        return None
    return mat.loc[gene, cols].astype(float)


def gene_vec_ensembl(mat, gene, cols):
    ens = GENES[gene]
    if ens not in mat.index:
        return None
    return mat.loc[ens, cols].astype(float)


# ---- GSE126044: responder vs non-responder
log126, meta126 = load_gse126044()
c126 = {"n": int(meta126.shape[0]),
        "responders": int((meta126["response"] == "responder").sum()),
        "non_responders": int((meta126["response"] == "non-responder").sum())}
resp_cols = meta126.loc[meta126["response"] == "responder", "col"].tolist()
nonr_cols = meta126.loc[meta126["response"] == "non-responder", "col"].tolist()
for g in GENES:
    v = gene_vec_symbol(log126, g, meta126["col"].tolist())
    if v is None:
        continue
    a = log126.loc[g, resp_cols].astype(float)
    b = log126.loc[g, nonr_cols].astype(float)
    p, rbc = mannwhitney(a, b)
    diff_rows.append(dict(cohort="GSE126044", contrast="responder_vs_nonresponder", gene=g,
                          group1="responder", n1=len(a), median1=float(np.median(a)),
                          group2="non-responder", n2=len(b), median2=float(np.median(b)),
                          p_mannwhitney=p, rank_biserial=rbc))
    save_box({"Responder": a.values, "Non-resp.": b.values},
             f"GSE126044 anti-PD-1\n{GENE_LABEL[g]}", "log2(CPM+1)",
             f"gse126044_{g}_response.png",
             annotate=f"Mann-Whitney p={p:.3g}")
# TROP2 vs PD-L1 correlation
for g1, g2 in [("TACSTD2", "CD274")]:
    x = log126.loc[g1, meta126["col"]].astype(float)
    y = log126.loc[g2, meta126["col"]].astype(float)
    rho, pr = stats.spearmanr(x, y)
    corr_rows.append(dict(cohort="GSE126044", gene_x=g1, gene_y=g2, spearman_rho=rho, p=pr, n=len(x)))
results["cohorts"]["GSE126044"] = c126

# ---- GSE135222: PFS survival by gene median split + Cox
log135, meta135 = load_gse135222()
meta135 = meta135.dropna(subset=["pfs_event", "pfs_time"]).copy()
c135 = {"n": int(meta135.shape[0]),
        "events": int(meta135["pfs_event"].sum())}
cols135 = meta135["col"].tolist()
kmf = KaplanMeierFitter()
for g in GENES:
    v = gene_vec_ensembl(log135, g, cols135)
    if v is None:
        continue
    vals = v.values.astype(float)
    med = np.median(vals)
    high = vals >= med
    T = meta135["pfs_time"].values.astype(float)
    E = meta135["pfs_event"].values.astype(float)
    lr = logrank_test(T[high], T[~high], E[high], E[~high])
    # Cox with continuous z-scored expression
    z = (vals - vals.mean()) / (vals.std(ddof=0) + 1e-9)
    cph_df = pd.DataFrame({"T": T, "E": E, "expr": z})
    try:
        cph = CoxPHFitter().fit(cph_df, "T", "E")
        hr = float(np.exp(cph.params_["expr"]))
        hr_p = float(cph.summary.loc["expr", "p"])
        hr_lo = float(np.exp(cph.confidence_intervals_.loc["expr"].iloc[0]))
        hr_hi = float(np.exp(cph.confidence_intervals_.loc["expr"].iloc[1]))
    except Exception:
        hr = hr_p = hr_lo = hr_hi = np.nan
    surv_rows.append(dict(cohort="GSE135222", gene=g, outcome="PFS",
                          logrank_p=float(lr.p_value),
                          n_high=int(high.sum()), n_low=int((~high).sum()),
                          cox_HR_per_SD=hr, cox_HR_p=hr_p, cox_HR_low=hr_lo, cox_HR_high=hr_hi))
    # KM plot
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    kmf.fit(T[high], E[high], label=f"{g} high (n={high.sum()})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="#C44E52")
    kmf.fit(T[~high], E[~high], label=f"{g} low (n={(~high).sum()})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="#4C72B0")
    ax.set_title(f"GSE135222 anti-PD-(L)1 PFS\n{GENE_LABEL[g]} (median split)", fontsize=9)
    ax.set_xlabel("Days")
    ax.set_ylabel("PFS probability")
    ax.text(0.98, 0.95, f"log-rank p={lr.p_value:.3g}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / f"gse135222_{g}_KM_PFS.png", dpi=150)
    plt.close(fig)
for g1, g2 in [("TACSTD2", "CD274")]:
    x = log135.loc[GENES[g1], cols135].astype(float)
    y = log135.loc[GENES[g2], cols135].astype(float)
    rho, pr = stats.spearmanr(x, y)
    corr_rows.append(dict(cohort="GSE135222", gene_x=g1, gene_y=g2, spearman_rho=rho, p=pr, n=len(x)))
results["cohorts"]["GSE135222"] = c135

# ---- GSE207422: MPR vs NMPR + residual-tumor correlation
exp207, meta207 = load_gse207422()
c207 = {"n": int(meta207.shape[0]),
        "MPR": int((meta207["mpr"] == "MPR").sum()),
        "NMPR": int((meta207["mpr"] == "NMPR").sum())}
mpr_cols = meta207.loc[meta207["mpr"] == "MPR", "Sample"].tolist()
nmpr_cols = meta207.loc[meta207["mpr"] == "NMPR", "Sample"].tolist()
for g in GENES:
    if g not in exp207.index:
        continue
    a = exp207.loc[g, mpr_cols].astype(float)
    b = exp207.loc[g, nmpr_cols].astype(float)
    p, rbc = mannwhitney(a, b)
    diff_rows.append(dict(cohort="GSE207422", contrast="MPR_vs_NMPR", gene=g,
                          group1="MPR", n1=len(a), median1=float(np.median(a)),
                          group2="NMPR", n2=len(b), median2=float(np.median(b)),
                          p_mannwhitney=p, rank_biserial=rbc))
    save_box({"MPR": a.values, "NMPR": b.values},
             f"GSE207422 neoadj chemo-IO\n{GENE_LABEL[g]}", "log2 TPM",
             f"gse207422_{g}_MPR.png", annotate=f"Mann-Whitney p={p:.3g}")
    # residual tumor correlation
    resid = meta207.set_index("Sample").loc[list(a.index) + list(b.index), "residual"].astype(float)
    expr_all = exp207.loc[g, list(a.index) + list(b.index)].astype(float)
    ok = resid.notna().values
    if ok.sum() >= 4:
        rho, pr = stats.spearmanr(expr_all.values[ok], resid.values[ok])
        corr_rows.append(dict(cohort="GSE207422", gene_x=g, gene_y="ResidualTumorFraction",
                              spearman_rho=rho, p=pr, n=int(ok.sum())))
for g1, g2 in [("TACSTD2", "CD274")]:
    x = exp207.loc[g1, meta207["Sample"]].astype(float)
    y = exp207.loc[g2, meta207["Sample"]].astype(float)
    rho, pr = stats.spearmanr(x, y)
    corr_rows.append(dict(cohort="GSE207422", gene_x=g1, gene_y=g2, spearman_rho=rho, p=pr, n=len(x)))
results["cohorts"]["GSE207422"] = c207

# --------------------------------------------------------------------------- write tables
diff_df = pd.DataFrame(diff_rows)
surv_df = pd.DataFrame(surv_rows)
corr_df = pd.DataFrame(corr_rows)
diff_df.to_csv(TAB / "differential_expression_by_outcome.csv", index=False)
surv_df.to_csv(TAB / "survival_association_GSE135222.csv", index=False)
corr_df.to_csv(TAB / "correlations.csv", index=False)

# cross-cohort TROP2 expression summary panel
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
panels = [
    ("GSE126044\n(anti-PD-1, log2 CPM)",
     {"Responder": log126.loc["TACSTD2", resp_cols].astype(float).values,
      "Non-resp.": log126.loc["TACSTD2", nonr_cols].astype(float).values}),
    ("GSE135222\n(anti-PD-(L)1, log2 TPM)",
     {"All": log135.loc[GENES["TACSTD2"], cols135].astype(float).values}),
    ("GSE207422\n(neoadj chemo-IO, log2 TPM)",
     {"MPR": exp207.loc["TACSTD2", mpr_cols].astype(float).values,
      "NMPR": exp207.loc["TACSTD2", nmpr_cols].astype(float).values}),
]
for ax, (title, d) in zip(axes, panels):
    labels = list(d.keys())
    vals = [d[k] for k in labels]
    ax.boxplot(vals, tick_labels=labels, showfliers=False, patch_artist=True)
    for i, v in enumerate(vals, start=1):
        ax.scatter(np.full(len(v), i) + rng.normal(0, 0.05, len(v)), v, s=16,
                   color="#4C72B0", edgecolor="black", linewidth=0.3, zorder=3)
    ax.set_title("TROP2 (TACSTD2)\n" + title, fontsize=9)
    ax.set_ylabel("expression")
fig.tight_layout()
fig.savefig(FIG / "trop2_expression_across_cohorts.png", dpi=150)
plt.close(fig)

results["summary_tables"] = {
    "differential_expression_by_outcome": "tables/differential_expression_by_outcome.csv",
    "survival_association_GSE135222": "tables/survival_association_GSE135222.csv",
    "correlations": "tables/correlations.csv",
}
with open(OUT / "analysis_summary.json", "w") as fh:
    json.dump(results, fh, indent=2)

# console recap
pd.set_option("display.width", 160, "display.max_columns", 20)
print("\n=== cohort sizes ===")
print(json.dumps(results["cohorts"], indent=2))
print("\n=== differential expression by outcome ===")
print(diff_df.to_string(index=False))
print("\n=== survival association (GSE135222 PFS) ===")
print(surv_df.to_string(index=False))
print("\n=== correlations ===")
print(corr_df.to_string(index=False))
print(f"\nOutputs written under {OUT}")
