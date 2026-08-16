#!/usr/bin/env python3
"""Analyze TACSTD2 / CLDN4 expression vs ICI outcomes across verified GEO series.

Datasets and outcome models
---------------------------
* GSE126044  tumor RNA-seq (raw counts)  -> responder vs non-responder (anti-PD-1)
* GSE135222  tumor RNA-seq (TPM)         -> progression-free survival (anti-PD-1/PD-L1)
* GSE182328  lung-tumor RNA-seq (counts) -> Akkermansia group (ICI-prognostic surrogate)
* GSE111414  PBMC CD8+ T cells (counts)  -> responder vs non-responder (QC / epithelial-null control)
* GSE136961  Oncomine 395-gene panel     -> EXCLUDED (TACSTD2/CLDN4 not on panel)

Normalisation
-------------
* Counts matrices -> log2 CPM = log2(1e6 * count / library_size + 1)
* TPM matrix      -> log2(TPM + 1)

Outputs
-------
Per-sample expression+outcome tables, a combined statistics table, and figures,
all under results/fable_geo_2019_2021/.
"""
import gzip
import json
import math
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
RES = ROOT / "results" / "fable_geo_2019_2021"
DL = RES / "downloads"
CLIN = RES / "clinical"
FIG = RES / "figures"
TAB = RES / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
ENS = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}

results_rows = []
persample_frames = {}


def log2cpm(counts_df):
    lib = counts_df.sum(axis=0)
    cpm = counts_df.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def auc_mannwhitney(pos, neg):
    """AUC via Mann-Whitney U statistic (prob a positive > a negative)."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan")
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = u / (len(pos) * len(neg))
    return auc, p


# --------------------------------------------------------------------------
# GSE126044 : responder vs non-responder (tumor RNA-seq counts)
# --------------------------------------------------------------------------
def analyze_gse126044():
    gse = "GSE126044"
    counts = pd.read_csv(DL / gse / f"{gse}_counts.txt.gz", sep="\t", index_col=0)
    clin = pd.read_csv(CLIN / f"{gse}_clinical.csv")
    # map column "Dis_XX" -> response using title "RNA-seq_Dis_XX"
    resp = {}
    for _, r in clin.iterrows():
        key = r["title"].replace("RNA-seq_", "")
        resp[key] = r["patient_response"]
    logcpm = log2cpm(counts)
    df = pd.DataFrame({"sample": counts.columns})
    df["response"] = df["sample"].map(resp)
    for g in GENES:
        df[g] = logcpm.loc[g, df["sample"]].values
    df.to_csv(TAB / f"{gse}_expr_outcome.csv", index=False)
    persample_frames[gse] = df
    for g in GENES:
        pos = df.loc[df.response == "responder", g]
        neg = df.loc[df.response == "non-responder", g]
        auc, p = auc_mannwhitney(pos, neg)
        results_rows.append({
            "dataset": gse, "compartment": "tumor", "gene": g,
            "outcome": "response (R vs NR)", "n": len(df),
            "n_group1": len(pos), "n_group2": len(neg),
            "median_responder": round(float(pos.median()), 3),
            "median_nonresponder": round(float(neg.median()), 3),
            "effect": f"AUC={auc:.3f}", "auc": round(auc, 3),
            "p_value": round(p, 4),
            "direction": "higher in responders" if pos.median() > neg.median() else "higher in non-responders",
        })
    _box(df, "response", GENES, gse, "log2 CPM", order=["responder", "non-responder"])
    return df


# --------------------------------------------------------------------------
# GSE135222 : PFS survival (tumor RNA-seq TPM, Ensembl ids)
# --------------------------------------------------------------------------
def analyze_gse135222():
    gse = "GSE135222"
    tpm = pd.read_csv(DL / gse / f"{gse}_GEO_RNA-seq_omicslab_exp.tsv.gz",
                      sep="\t", index_col=0)
    # strip Ensembl version
    tpm.index = [i.split(".")[0] for i in tpm.index]
    clin = pd.read_csv(CLIN / f"{gse}_clinical.csv")
    # map title "NSCLC 990" -> column "NSCLC990"
    clin["col"] = clin["title"].str.replace(" ", "", regex=False)
    clin = clin[clin["col"].isin(tpm.columns)].copy()
    df = clin[["col", "progression_free_survival_pfs", "pfs_time"]].rename(
        columns={"col": "sample", "progression_free_survival_pfs": "pfs_event",
                 "pfs_time": "pfs_time"})
    df["pfs_event"] = df["pfs_event"].astype(int)
    df["pfs_time"] = df["pfs_time"].astype(float)
    for g in GENES:
        df[g] = np.log2(tpm.loc[ENS[g], df["sample"]].values.astype(float) + 1)
    df.to_csv(TAB / f"{gse}_expr_outcome.csv", index=False)
    persample_frames[gse] = df

    for g in GENES:
        # continuous Cox
        cph = CoxPHFitter()
        cph.fit(df[[g, "pfs_time", "pfs_event"]], duration_col="pfs_time",
                event_col="pfs_event")
        hr = math.exp(cph.params_[g])
        p_cox = cph.summary.loc[g, "p"]
        # median split KM + logrank
        med = df[g].median()
        high = df[df[g] > med]
        low = df[df[g] <= med]
        lr = logrank_test(high["pfs_time"], low["pfs_time"],
                          high["pfs_event"], low["pfs_event"])
        results_rows.append({
            "dataset": gse, "compartment": "tumor", "gene": g,
            "outcome": "PFS (Cox continuous / KM median-split)",
            "n": len(df), "n_group1": len(high), "n_group2": len(low),
            "median_responder": "", "median_nonresponder": "",
            "effect": f"HR={hr:.3f} (per log2TPM); logrank p={lr.p_value:.4f}",
            "auc": "", "p_value": round(float(p_cox), 4),
            "direction": "higher expr -> shorter PFS" if hr > 1 else "higher expr -> longer PFS",
        })
        _km(high, low, g, gse, med)
    return df


# --------------------------------------------------------------------------
# GSE182328 : Akkermansia group (ICI-prognostic surrogate) - lung tumor counts
# --------------------------------------------------------------------------
def analyze_gse182328():
    gse = "GSE182328"
    counts = pd.read_csv(DL / gse / f"{gse}_Gene_counts_matrix.txt.gz",
                         sep="\t", index_col=0)
    clin = pd.read_csv(CLIN / f"{gse}_clinical.csv")
    # columns S01008.. align to clinical rows in same order
    order = {t: a for t, a in zip(clin["title"], clin["akk_metaominer"])}
    stage = {t: s for t, s in zip(clin["title"], clin["disease_state"])}
    logcpm = log2cpm(counts)
    df = pd.DataFrame({"sample": counts.columns})
    df["akk"] = df["sample"].map(order)
    df["stage"] = df["sample"].map(stage)
    for g in GENES:
        df[g] = logcpm.loc[g, df["sample"]].values
    df.to_csv(TAB / f"{gse}_expr_outcome.csv", index=False)
    persample_frames[gse] = df
    for g in GENES:
        pos = df.loc[df.akk == "detectable", g]
        neg = df.loc[df.akk == "not_detectable", g]
        auc, p = auc_mannwhitney(pos, neg)
        results_rows.append({
            "dataset": gse, "compartment": "tumor", "gene": g,
            "outcome": "Akkermansia detectable vs not (ICI-prognostic surrogate)",
            "n": len(df), "n_group1": len(pos), "n_group2": len(neg),
            "median_responder": round(float(pos.median()), 3),
            "median_nonresponder": round(float(neg.median()), 3),
            "effect": f"AUC={auc:.3f}", "auc": round(auc, 3),
            "p_value": round(p, 4),
            "direction": "higher in Akk+" if pos.median() > neg.median() else "higher in Akk-",
        })
    _box(df, "akk", GENES, gse, "log2 CPM",
         order=["detectable", "not_detectable"])
    return df


# --------------------------------------------------------------------------
# GSE111414 : peripheral CD8+ T cells - epithelial-null control
# --------------------------------------------------------------------------
def analyze_gse111414():
    gse = "GSE111414"
    counts = pd.read_csv(DL / gse / f"{gse}_gene_counts.csv.gz", index_col=1)
    meta_cols = ["ENTREZID", "GENENAME", "LENGTH"]
    expr = counts.drop(columns=[c for c in meta_cols if c in counts.columns])
    clin = pd.read_csv(CLIN / f"{gse}_clinical.csv")
    # columns look like 373_N1 ; map via subject id digits + timepoint.
    # response taken from source_name ("responder_*" / "nonresponder_*").
    resp = {}
    for _, r in clin.iterrows():
        sid = r["subject_id"].split()[-1]  # e.g. BS373
        num = "".join(ch for ch in sid if ch.isdigit())
        tp = r["timepoint"]
        label = "responder" if str(r["source_name"]).startswith("responder") else "non-responder"
        resp[f"{num}_{tp}"] = label
    logcpm = log2cpm(expr)
    df = pd.DataFrame({"sample": expr.columns})
    df["response"] = df["sample"].map(resp)
    for g in GENES:
        df[g] = logcpm.loc[g, df["sample"]].values if g in logcpm.index else np.nan
    df.to_csv(TAB / f"{gse}_expr_outcome.csv", index=False)
    persample_frames[gse] = df
    for g in GENES:
        vals = df[g].dropna()
        pos = df.loc[df.response == "responder", g].dropna()
        neg = df.loc[df.response == "non-responder", g].dropna()
        auc, p = auc_mannwhitney(pos, neg)
        results_rows.append({
            "dataset": gse, "compartment": "PBMC CD8+ T cell", "gene": g,
            "outcome": "response (R vs NR) [epithelial-null control]",
            "n": len(df), "n_group1": len(pos), "n_group2": len(neg),
            "median_responder": round(float(pos.median()), 3) if len(pos) else "",
            "median_nonresponder": round(float(neg.median()), 3) if len(neg) else "",
            "effect": f"mean log2CPM={vals.mean():.3f} (near-zero: epithelial gene not expressed in CD8 T)",
            "auc": round(auc, 3) if not math.isnan(auc) else "",
            "p_value": round(p, 4) if not math.isnan(p) else "",
            "direction": "not interpretable (compartment lacks epithelium)",
        })
    return df


# --------------------------------------------------------------------------
# plotting helpers
# --------------------------------------------------------------------------
def _box(df, group, genes, gse, ylab, order):
    fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4))
    if len(genes) == 1:
        axes = [axes]
    for ax, g in zip(axes, genes):
        data = [df.loc[df[group] == lvl, g].values for lvl in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        for i, lvl in enumerate(order, start=1):
            y = df.loc[df[group] == lvl, g].values
            x = np.random.normal(i, 0.05, size=len(y))
            ax.scatter(x, y, alpha=0.7, s=18, color="#377eb8")
        ax.set_title(f"{gse}: {g}")
        ax.set_ylabel(ylab)
        ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(FIG / f"{gse}_{group}_boxplot.png", dpi=130)
    plt.close(fig)


def _km(high, low, g, gse, med):
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(5, 4))
    kmf.fit(high["pfs_time"], high["pfs_event"], label=f"{g} high (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=False)
    kmf.fit(low["pfs_time"], low["pfs_event"], label=f"{g} low (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=False)
    ax.set_title(f"{gse}: PFS by {g} (median split)")
    ax.set_xlabel("PFS time (days)")
    ax.set_ylabel("PFS probability")
    fig.tight_layout()
    fig.savefig(FIG / f"{gse}_{g}_KM.png", dpi=130)
    plt.close(fig)


def main():
    analyze_gse126044()
    analyze_gse135222()
    analyze_gse182328()
    analyze_gse111414()
    res = pd.DataFrame(results_rows)
    res.to_csv(TAB / "combined_TACSTD2_CLDN4_vs_ICI_outcomes.csv", index=False)
    print(res.to_string(index=False))
    print("\nWrote", TAB / "combined_TACSTD2_CLDN4_vs_ICI_outcomes.csv")
    print("Figures in", FIG)


if __name__ == "__main__":
    main()
