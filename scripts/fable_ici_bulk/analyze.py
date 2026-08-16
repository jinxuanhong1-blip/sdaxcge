#!/usr/bin/env python3
"""TACSTD2 and CLDN4 expression vs ICI response in open human bulk lung cohorts.

Usable cohorts (see notes/fable_ici_bulk/geo_verification.md):
  GSE126044  -- NSCLC anti-PD-1, raw counts, binary responder/non-responder
  GSE135222  -- NSCLC anti-PD-(L)1, TPM, PFS survival
  GSE166449  -- lung cancer immunotherapy, TPM, binary responder/non-responder
  GSE207422  -- NSCLC neoadjuvant anti-PD-1+chemo, bulk log2TPM, MPR/NMPR

Unusable (target genes not on targeted immune panels):
  GSE136961 (Oncomine Immune Response, 395 genes), GSE93157 (NanoString 730).
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(__file__))
from geo_utils import parse_series_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA = "/tmp/ici_bulk_data"
RES = "/workspace/results/fable_ici_bulk"
os.makedirs(RES, exist_ok=True)
os.makedirs(os.path.join(RES, "figures"), exist_ok=True)
os.makedirs(os.path.join(RES, "tables"), exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
ENSG = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}

binary_rows = []   # collected stats for binary comparisons
survival_rows = []  # collected stats for survival
long_records = []  # per-sample expression + label, for reproducibility


def mannwhitney_block(gse, gene, values_resp, values_non, unit):
    x = np.asarray(values_resp, float)
    y = np.asarray(values_non, float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    n1, n2 = len(x), len(y)
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    auc = U / (n1 * n2)  # P(responder value > non-responder value)
    rbc = 2 * auc - 1     # rank-biserial correlation
    return {
        "dataset": gse,
        "gene": gene,
        "unit": unit,
        "n_responder": n1,
        "n_nonresponder": n2,
        "median_responder": float(np.median(x)),
        "median_nonresponder": float(np.median(y)),
        "mean_responder": float(np.mean(x)),
        "mean_nonresponder": float(np.mean(y)),
        "mannwhitney_U": float(U),
        "auc_resp_gt_nonresp": float(auc),
        "rank_biserial": float(rbc),
        "p_value": float(p),
    }


def boxplot_binary(gse, gene, resp, non, unit, fname):
    fig, ax = plt.subplots(figsize=(3.4, 4))
    data = [np.asarray(non, float), np.asarray(resp, float)]
    labels = [f"Non-resp\n(n={len(non)})", f"Resp\n(n={len(resp)})"]
    bp = ax.boxplot(data, tick_labels=labels, widths=0.6, showfliers=False,
                    patch_artist=True)
    for patch, c in zip(bp["boxes"], ["#c9d6df", "#f6b26b"]):
        patch.set_facecolor(c)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(rng.normal(i, 0.06, len(d)), d, s=18, color="#333",
                   alpha=0.7, zorder=3)
    ax.set_ylabel(f"{gene} ({unit})")
    ax.set_title(f"{gse}: {gene}")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "figures", fname), dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# GSE126044 -- raw counts -> log2 CPM ; binary responder/non-responder
# ---------------------------------------------------------------------------
def run_gse126044():
    gse = "GSE126044"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    resp_map = {}
    for s in meta["samples"]:
        # title like "RNA-seq_Dis_01"; counts columns are "Dis_01"
        col = s["title"].replace("RNA-seq_", "")
        resp_map[col] = s.get("patient response")
    counts = pd.read_csv(f"{DATA}/{gse}_counts.txt.gz", sep="\t", index_col=0)
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1)
    for gene in GENES:
        vals = logcpm.loc[gene]
        resp = [vals[c] for c in vals.index if resp_map.get(c) == "responder"]
        non = [vals[c] for c in vals.index if resp_map.get(c) == "non-responder"]
        for c in vals.index:
            long_records.append({"dataset": gse, "sample": c, "gene": gene,
                                 "value": float(vals[c]), "unit": "log2CPM",
                                 "group": resp_map.get(c)})
        binary_rows.append(mannwhitney_block(gse, gene, resp, non, "log2CPM"))
        boxplot_binary(gse, gene, resp, non, "log2CPM", f"{gse}_{gene}.png")


# ---------------------------------------------------------------------------
# GSE166449 -- TPM -> log2(TPM+1) ; binary responder/non-responder
# ---------------------------------------------------------------------------
def run_gse166449():
    gse = "GSE166449"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    resp_map = {}
    for s in meta["samples"]:
        col = s.get("description")  # SMC_IO_C0xx == expression column
        title = s["title"].lower()
        if "nonresponder" in title:
            resp_map[col] = "non-responder"
        elif "responder" in title:
            resp_map[col] = "responder"
    tpm = pd.read_csv(f"{DATA}/{gse}_TPM.txt.gz", sep="\t", index_col=0)
    logtpm = np.log2(tpm + 1)
    for gene in GENES:
        vals = logtpm.loc[gene]
        resp = [vals[c] for c in vals.index if resp_map.get(c) == "responder"]
        non = [vals[c] for c in vals.index if resp_map.get(c) == "non-responder"]
        for c in vals.index:
            long_records.append({"dataset": gse, "sample": c, "gene": gene,
                                 "value": float(vals[c]), "unit": "log2(TPM+1)",
                                 "group": resp_map.get(c)})
        binary_rows.append(mannwhitney_block(gse, gene, resp, non, "log2(TPM+1)"))
        boxplot_binary(gse, gene, resp, non, "log2(TPM+1)", f"{gse}_{gene}.png")


# ---------------------------------------------------------------------------
# GSE207422 -- bulk log2TPM (already log2) ; MPR vs NMPR (from metadata xlsx)
# ---------------------------------------------------------------------------
def run_gse207422():
    gse = "GSE207422"
    md = pd.read_excel(f"{DATA}/{gse}_metadata.xlsx")
    md = md.dropna(subset=["Sample", "Pathologic Response"])
    md = md[md["Pathologic Response"].astype(str).str.contains("MPR")]
    def grp(x):
        x = str(x)
        if x.startswith("NMPR"):
            return "non-responder"
        if "MPR" in x:  # MPR or MPR (pCR)
            return "responder"
        return None
    resp_map = {row["Sample"]: grp(row["Pathologic Response"])
                for _, row in md.iterrows()}
    expr = pd.read_csv(f"{DATA}/{gse}_log2TPM.txt.gz", sep="\t", index_col=0)
    cols = [c for c in expr.columns if c in resp_map]
    for gene in GENES:
        vals = expr.loc[gene, cols]
        resp = [vals[c] for c in cols if resp_map.get(c) == "responder"]
        non = [vals[c] for c in cols if resp_map.get(c) == "non-responder"]
        for c in cols:
            long_records.append({"dataset": gse, "sample": c, "gene": gene,
                                 "value": float(vals[c]), "unit": "log2TPM",
                                 "group": resp_map.get(c)})
        binary_rows.append(mannwhitney_block(gse, gene, resp, non, "log2TPM"))
        boxplot_binary(gse, gene, resp, non, "log2TPM", f"{gse}_{gene}.png")


# ---------------------------------------------------------------------------
# GSE135222 -- TPM -> log2(TPM+1) ; PFS survival (Cox + median-split logrank)
# ---------------------------------------------------------------------------
def run_gse135222():
    gse = "GSE135222"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    surv = {}
    for s in meta["samples"]:
        col = s["title"].replace(" ", "")  # "NSCLC 990" -> "NSCLC990"
        pfs_event = s.get("progression-free survival (pfs)")
        pfs_time = s.get("pfs.time")
        if pfs_event is not None and pfs_time is not None:
            surv[col] = (int(pfs_event), float(pfs_time))
    tpm = pd.read_csv(f"{DATA}/{gse}_exp.tsv.gz", sep="\t", index_col=0)
    tpm.index = [i.split(".")[0] for i in tpm.index]  # strip ENSG version
    logtpm = np.log2(tpm + 1)
    for gene in GENES:
        row = logtpm.loc[ENSG[gene]]
        recs = []
        for c in row.index:
            if c in surv:
                ev, t = surv[c]
                recs.append({"sample": c, "expr": float(row[c]),
                             "event": ev, "time": t})
                long_records.append({"dataset": gse, "sample": c, "gene": gene,
                                     "value": float(row[c]), "unit": "log2(TPM+1)",
                                     "group": f"PFSevent={ev}"})
        df = pd.DataFrame(recs)
        # Cox on continuous expression
        cph = CoxPHFitter()
        cph.fit(df[["time", "event", "expr"]], duration_col="time",
                event_col="event")
        hr = float(np.exp(cph.params_["expr"]))
        ci_low = float(np.exp(cph.confidence_intervals_.loc["expr"].iloc[0]))
        ci_high = float(np.exp(cph.confidence_intervals_.loc["expr"].iloc[1]))
        cox_p = float(cph.summary.loc["expr", "p"])
        # Median split logrank + KM
        med = df["expr"].median()
        high = df[df["expr"] >= med]
        low = df[df["expr"] < med]
        lr = logrank_test(high["time"], low["time"],
                          event_observed_A=high["event"],
                          event_observed_B=low["event"])
        survival_rows.append({
            "dataset": gse, "gene": gene, "unit": "log2(TPM+1)",
            "n": len(df), "n_events": int(df["event"].sum()),
            "cox_HR_per_log2unit": hr, "cox_CI_low": ci_low,
            "cox_CI_high": ci_high, "cox_p": cox_p,
            "logrank_median_split_p": float(lr.p_value),
            "n_high": len(high), "n_low": len(low),
        })
        # KM plot
        fig, ax = plt.subplots(figsize=(4.2, 4))
        kmf = KaplanMeierFitter()
        kmf.fit(high["time"], high["event"], label=f"High {gene} (n={len(high)})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="#f6b26b")
        kmf.fit(low["time"], low["event"], label=f"Low {gene} (n={len(low)})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="#6fa8dc")
        ax.set_title(f"{gse}: {gene} (logrank p={lr.p_value:.3f})")
        ax.set_xlabel("PFS time (days)")
        ax.set_ylabel("PFS probability")
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "figures", f"{gse}_{gene}_KM.png"), dpi=140)
        plt.close(fig)


def main():
    run_gse126044()
    run_gse166449()
    run_gse207422()
    run_gse135222()

    bdf = pd.DataFrame(binary_rows)
    sdf = pd.DataFrame(survival_rows)

    # BH-FDR across all primary p-values (binary MWU + survival Cox)
    pvals = list(bdf["p_value"]) + list(sdf["cox_p"])
    rej, q, _, _ = multipletests(pvals, method="fdr_bh")
    bdf["p_fdr_bh"] = q[:len(bdf)]
    sdf["cox_p_fdr_bh"] = q[len(bdf):]

    bdf.to_csv(os.path.join(RES, "tables", "binary_response_stats.csv"), index=False)
    sdf.to_csv(os.path.join(RES, "tables", "survival_stats.csv"), index=False)
    pd.DataFrame(long_records).to_csv(
        os.path.join(RES, "tables", "per_sample_expression.csv"), index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print("=== Binary response (Mann-Whitney U) ===")
    print(bdf.to_string(index=False))
    print("\n=== Survival (GSE135222, PFS) ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
