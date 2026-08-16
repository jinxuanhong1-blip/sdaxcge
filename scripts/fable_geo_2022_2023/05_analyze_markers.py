#!/usr/bin/env python3
"""
Association of TACSTD2 (Trop-2) and CLDN4 (Claudin-4) tumour expression with
immune-checkpoint-inhibitor (ICI) outcomes.

Cohorts actually usable for these two genes (see notes/verification):
  * PRIMARY (lung)      GSE207422 - NSCLC, neoadjuvant anti-PD-1 + chemo,
                        n=24 baseline (pre-treatment) tumours, bulk RNA-seq
                        log2 TPM.  Outcome = major pathologic response
                        (MPR incl. pCR) vs non-MPR (NMPR); + residual tumour %,
                        + RECIST.
  * CROSS-CHECK (non-lung, methodological only) GSE243238 - acral melanoma,
                        ICI-treated, n=20 baseline tumours, bulk RNA-seq raw
                        counts.  Outcome = clinical benefit (CB) yes/no.

The two big NSCLC ICI survival cohorts (GSE161537, GSE162520) use a ~2.5k-gene
targeted immune panel that does NOT contain TACSTD2 or CLDN4, so they cannot be
used for these markers (documented, not analysed here).

Outputs: results/fable_geo_2022_2023/analysis/*
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
dl = import_module("04_download_and_parse")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "fable_geo_2022_2023"
DATA = OUT / "data"
ANA = OUT / "analysis"
ANA.mkdir(parents=True, exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
RESULTS = []


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def group_test(gene, vals_by_group, g1, g2, cohort, outcome):
    a = np.asarray(vals_by_group[g1], float)
    b = np.asarray(vals_by_group[g2], float)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cliffs_delta(a, b)
    row = {
        "cohort": cohort, "gene": gene, "outcome": outcome,
        "group1": g1, "n1": len(a), "median1": round(float(np.median(a)), 3),
        "group2": g2, "n2": len(b), "median2": round(float(np.median(b)), 3),
        "mannwhitney_U": float(u), "p_value": round(float(p), 4),
        "cliffs_delta": round(float(d), 3),
    }
    RESULTS.append(row)
    return row


def boxplot(gene, groups, data, title, fname, ylabel):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    vals = [data[g] for g in groups]
    bp = ax.boxplot(vals, tick_labels=groups, patch_artist=True, widths=0.55)
    colors = ["#4C72B0", "#C44E52", "#55A868"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    for i, g in enumerate(groups, 1):
        y = data[g]
        x = np.random.normal(i, 0.05, len(y))
        ax.scatter(x, y, color="black", s=16, zorder=3, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)


def scatter(gene, x, y, title, fname, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    ax.scatter(x, y, color="#4C72B0", s=28, alpha=0.8)
    if len(x) > 2:
        rho, p = stats.spearmanr(x, y)
        b, a = np.polyfit(x, y, 1)
        xs = np.linspace(min(x), max(x), 50)
        ax.plot(xs, b * xs + a, color="#C44E52", lw=1.5)
        ax.set_title(f"{title}\nSpearman rho={rho:.2f}, p={p:.3f}", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------- GSE207422
def analyze_gse207422():
    expr = pd.read_csv(
        DATA / "GSE207422__GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        sep="\t", index_col=0)
    meta = pd.read_excel(DATA / "GSE207422_bulk_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"])
    meta = meta[meta["Sample"].astype(str).str.startswith("R")].copy()
    meta = meta[meta["Sample"].isin(expr.columns)]
    meta["responder"] = np.where(
        meta["Pathologic Response"].str.contains("MPR", na=False) &
        ~meta["Pathologic Response"].str.contains("NMPR", na=False),
        "MPR", "NMPR")
    meta["recist_resp"] = np.where(
        meta["RECIST"].isin(["CR", "PR"]), "CR/PR", "SD/PD")
    per = meta[["Sample", "Patient", "Pathologic Response", "responder",
                "RECIST", "recist_resp", "Residual Tumor"]].copy()
    for g in GENES:
        per[g + "_log2TPM"] = per["Sample"].map(expr.loc[g])
    per.to_csv(ANA / "GSE207422_per_sample.csv", index=False)

    for g in GENES:
        col = g + "_log2TPM"
        by = {grp: per.loc[per.responder == grp, col].tolist()
              for grp in ["MPR", "NMPR"]}
        group_test(g, by, "MPR", "NMPR", "GSE207422", "pathologic_response")
        boxplot(g, ["MPR", "NMPR"], by,
                f"GSE207422 (NSCLC, n={len(per)})\n{g} vs pathologic response",
                f"GSE207422_{g}_MPRvsNMPR.png", f"{g} (log2 TPM)")
        # RECIST
        by2 = {grp: per.loc[per.recist_resp == grp, col].tolist()
               for grp in ["CR/PR", "SD/PD"]}
        if all(len(v) > 0 for v in by2.values()):
            group_test(g, by2, "CR/PR", "SD/PD", "GSE207422", "RECIST")
        # residual tumour % correlation (continuous outcome; higher = worse)
        sub = per.dropna(subset=["Residual Tumor", col])
        xv = sub["Residual Tumor"].astype(float).values
        yv = sub[col].astype(float).values
        rho, prho = stats.spearmanr(xv, yv)
        RESULTS.append({
            "cohort": "GSE207422", "gene": g, "outcome": "residual_tumor_frac",
            "group1": "spearman", "n1": len(xv), "median1": "",
            "group2": "", "n2": "", "median2": "",
            "mannwhitney_U": "", "p_value": round(float(prho), 4),
            "cliffs_delta": round(float(rho), 3),
        })
        scatter(g, xv, yv,
                f"GSE207422 {g} vs residual tumour",
                f"GSE207422_{g}_residual.png",
                "Residual tumour fraction", f"{g} (log2 TPM)")
    return per


# ---------------------------------------------------------------- GSE243238
def analyze_gse243238():
    counts = pd.read_csv(
        DATA / "GSE243238__acral.20.counts.raw.mat.csv.gz", index_col=0)
    # CPM -> log2
    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    meta = dl.parse_series_matrix_meta("GSE243238")
    # map cbl id -> clinical benefit
    cbl = meta["cbl id"] if "cbl id" in meta else None
    cb = meta["clinical benefit_(cb)"]
    m = pd.DataFrame({"gsm": meta.index, "cbl_id": cbl.values,
                      "cb": cb.values})
    m = m[m["cbl_id"].isin(logcpm.columns)].copy()
    m = m[m["cb"].isin(["Yes", "No"])].copy()
    for g in GENES:
        if g in logcpm.index:
            m[g + "_log2CPM"] = m["cbl_id"].map(logcpm.loc[g])
    m.to_csv(ANA / "GSE243238_per_sample.csv", index=False)
    for g in GENES:
        col = g + "_log2CPM"
        if col not in m:
            continue
        by = {grp: m.loc[m.cb == grp, col].tolist() for grp in ["Yes", "No"]}
        if all(len(v) >= 2 for v in by.values()):
            group_test(g, by, "Yes", "No", "GSE243238", "clinical_benefit")
            boxplot(g, ["Yes", "No"], by,
                    f"GSE243238 (acral melanoma, n={len(m)})\n{g} vs clinical benefit",
                    f"GSE243238_{g}_CB.png", f"{g} (log2 CPM)")
    return m


def main():
    p1 = analyze_gse207422()
    print("GSE207422 per-sample:\n",
          p1[["Sample", "responder", "RECIST", "Residual Tumor",
              "TACSTD2_log2TPM", "CLDN4_log2TPM"]].to_string(index=False),
          file=sys.stderr)
    m2 = analyze_gse243238()
    print("\nGSE243238 per-sample:\n", m2.to_string(index=False),
          file=sys.stderr)

    res = pd.DataFrame(RESULTS)
    res.to_csv(ANA / "marker_outcome_tests.csv", index=False)
    print("\n=== Association tests ===", file=sys.stderr)
    print(res.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
