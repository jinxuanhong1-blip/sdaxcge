#!/usr/bin/env python3
"""Malignant TACSTD2 / CLDN4 vs ICI response and vs T/NK.

Cohorts (accessions taken from live GEO records, never invented):

  GSE291670 (2025-03-31) — 6 NSCLC tumors after neoadjuvant anlotinib +
      camrelizumab (PD-1). Sample titles encode MPR vs Non-MPR. Appeared in
      a prior 2024–2026 candidate TSV but was never analyzed for this endpoint.

  GSE325414 (2026-06-10) — NEW 2026 NSCLC BD-Rhapsody scRNA (not ICI-treated;
      PEF ablation treat-and-resect). Author cell-type labels. Used for
      TACSTD2/CLDN4 vs T/NK (and vs TLS) because no NEW 2026 lung-primary
      ICI scRNA series with response labels exists beyond already-listed GSEs.

Outputs: results/hunt_geo2026_sc/
"""
from __future__ import annotations

import gzip
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr

sc.settings.verbosity = 1
DATA = "/workspace/data"
OUT = "/workspace/results/hunt_geo2026_sc"
os.makedirs(OUT, exist_ok=True)

# Verbatim from GEO sample titles / characteristics (GSE291670).
GSE291670_SAMPLES = {
    "GSM8839599_MPR-1": "MPR",
    "GSM8839600_MPR-2": "MPR",
    "GSM8839601_MPR-3": "MPR",
    "GSM8839602_Non-MPR-1": "non-MPR",
    "GSM8839603_Non-MPR-2": "non-MPR",
    "GSM8839604_Non-MPR-3": "non-MPR",
}

MARKERS = {
    "Malignant/Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "ELF3"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4", "NKG7", "GNLY", "KLRD1", "TRAC"],
    "B/Plasma": ["CD79A", "MS4A1", "MZB1", "IGHG1", "CD19"],
    "Myeloid": ["LYZ", "CD68", "CD14", "FCGR3A", "C1QA", "AIF1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3"],
}

TNK_MAJOR = {"Tcell_CD4", "Tcell_CD8", "NKcell", "Tcell_CD4-_CD8-", "Tcell_CD4+_CD8+"}
TNK_SUB = {"Lymphoid"}
MAL_SUB = {"EpithelialCells_TumorCells"}


def mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return {"n_a": int(len(a)), "n_b": int(len(b)), "p": None, "U": None}
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    return {"n_a": int(len(a)), "n_b": int(len(b)), "p": float(p), "U": float(U)}


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return {"n": int(m.sum()), "rho": None, "p": None}
    rho, p = spearmanr(x[m], y[m])
    return {"n": int(m.sum()), "rho": float(rho), "p": float(p)}


# ---------------------------------------------------------------------------
# GSE291670
# ---------------------------------------------------------------------------
def read_10x_sample(prefix):
    d = os.path.join(DATA, "GSE291670")
    M = scipy.io.mmread(os.path.join(d, f"{prefix}_matrix.mtx.gz")).tocsr().astype(np.float32)
    with gzip.open(os.path.join(d, f"{prefix}_barcodes.tsv.gz"), "rt") as f:
        barcodes = [l.strip().split("\t")[0] for l in f]
    genes = []
    with gzip.open(os.path.join(d, f"{prefix}_features.tsv.gz"), "rt") as f:
        for l in f:
            parts = l.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    ad = sc.AnnData(M.T.tocsr())
    ad.var_names = pd.Index(genes)
    ad.var_names_make_unique()
    ad.obs_names = [f"{prefix}|{b}" for b in barcodes]
    ad.obs["sample"] = prefix
    return ad


def analyze_gse291670():
    print("==== GSE291670 ====")
    ads = []
    for prefix, resp in GSE291670_SAMPLES.items():
        a = read_10x_sample(prefix)
        a.obs["response"] = resp
        a.obs["responder"] = "Responder" if resp == "MPR" else "Non-responder"
        print(f"  {prefix}: {a.n_obs} cells ({resp})")
        ads.append(a)
    ad = sc.concat(ads, join="outer")
    ad.obs_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    ad = ad[ad.obs["pct_counts_mt"] < 20].copy()
    ad = ad[ad.obs["n_genes_by_counts"] < 8000].copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw = ad
    hv = ad.copy()
    sc.pp.highly_variable_genes(hv, n_top_genes=2000)
    hv = hv[:, hv.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=30)
    sc.pp.neighbors(hv, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(hv, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    ad.obs["leiden"] = hv.obs["leiden"].values
    for lin, gl in MARKERS.items():
        present = [g for g in gl if g in ad.var_names]
        sc.tl.score_genes(ad, present, score_name=f"score_{lin}")
    score_cols = [f"score_{l}" for l in MARKERS]
    cl_scores = ad.obs.groupby("leiden", observed=True)[score_cols].mean()
    cl_assign = cl_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
    ad.obs["lineage"] = ad.obs["leiden"].map(cl_assign).astype(str)

    for g in ("TACSTD2", "CLDN4"):
        ad.obs[g] = np.asarray(ad[:, g].X.todense()).ravel() if g in ad.var_names else np.nan

    rows = []
    for s, sub in ad.obs.groupby("sample", observed=True):
        mal = sub[sub["lineage"] == "Malignant/Epithelial"]
        rows.append({
            "sample": s,
            "response": sub["response"].iloc[0],
            "responder": sub["responder"].iloc[0],
            "n_cells": int(len(sub)),
            "n_malignant": int(len(mal)),
            "frac_malignant": float(len(mal) / len(sub)),
            "n_TNK": int((sub["lineage"] == "T/NK").sum()),
            "frac_TNK": float((sub["lineage"] == "T/NK").mean()),
            "malignant_TACSTD2_mean": float(mal["TACSTD2"].mean()) if len(mal) else np.nan,
            "malignant_CLDN4_mean": float(mal["CLDN4"].mean()) if len(mal) else np.nan,
            "malignant_TACSTD2_pct_pos": float((mal["TACSTD2"] > 0).mean() * 100) if len(mal) else np.nan,
            "malignant_CLDN4_pct_pos": float((mal["CLDN4"] > 0).mean() * 100) if len(mal) else np.nan,
        })
    ps = pd.DataFrame(rows).sort_values(["responder", "sample"])
    ps.to_csv(f"{OUT}/GSE291670_per_sample_summary.csv", index=False)
    (ad.obs.groupby(["sample", "lineage"], observed=True).size()
     .unstack(fill_value=0).to_csv(f"{OUT}/GSE291670_lineage_counts.csv"))

    r = ps[ps["responder"] == "Responder"]
    nr = ps[ps["responder"] == "Non-responder"]
    res = {
        "cohort": "GSE291670",
        "pdat": "2025/03/31",
        "already_listed_as_candidate": True,
        "analyzed_for_tacstd2_cldn4_vs_tnk_before": False,
        "regimen": "neoadjuvant anlotinib + camrelizumab (PD-1), from GEO sample characteristics",
        "n_samples": int(ad.obs["sample"].nunique()),
        "n_cells_qc": int(ad.n_obs),
        "n_malignant": int((ad.obs["lineage"] == "Malignant/Epithelial").sum()),
        "n_TNK": int((ad.obs["lineage"] == "T/NK").sum()),
        "lineage_counts": ad.obs["lineage"].value_counts().to_dict(),
        "per_sample_TACSTD2_vs_MPR": mwu(r["malignant_TACSTD2_mean"], nr["malignant_TACSTD2_mean"]),
        "per_sample_CLDN4_vs_MPR": mwu(r["malignant_CLDN4_mean"], nr["malignant_CLDN4_mean"]),
        "mean_TACSTD2_MPR": float(r["malignant_TACSTD2_mean"].mean()),
        "mean_TACSTD2_nonMPR": float(nr["malignant_TACSTD2_mean"].mean()),
        "mean_CLDN4_MPR": float(r["malignant_CLDN4_mean"].mean()),
        "mean_CLDN4_nonMPR": float(nr["malignant_CLDN4_mean"].mean()),
        "spearman_TACSTD2_vs_fracTNK": spear(ps["malignant_TACSTD2_mean"], ps["frac_TNK"]),
        "spearman_CLDN4_vs_fracTNK": spear(ps["malignant_CLDN4_mean"], ps["frac_TNK"]),
    }
    mal = ad.obs[ad.obs["lineage"] == "Malignant/Epithelial"]
    res["cell_TACSTD2_vs_MPR"] = mwu(
        mal.loc[mal["responder"] == "Responder", "TACSTD2"],
        mal.loc[mal["responder"] == "Non-responder", "TACSTD2"])
    res["cell_CLDN4_vs_MPR"] = mwu(
        mal.loc[mal["responder"] == "Responder", "CLDN4"],
        mal.loc[mal["responder"] == "Non-responder", "CLDN4"])
    json.dump(res, open(f"{OUT}/GSE291670_results.json", "w"), indent=2)
    ad.obs[["sample", "response", "responder", "leiden", "lineage", "TACSTD2", "CLDN4",
            "n_genes_by_counts", "pct_counts_mt"]].to_csv(
        f"{OUT}/GSE291670_cell_annotations.csv.gz", compression="gzip")
    print(json.dumps(res, indent=2))
    return ps, res


# ---------------------------------------------------------------------------
# GSE325414 — stream only the genes we need (264M nnz full matrix)
# ---------------------------------------------------------------------------
def stream_gene_subset(mtx_gz, n_genes, n_cells, keep_idx):
    """Return cells x keep_genes CSR from a genes x cells MatrixMarket file."""
    keep_set = set(keep_idx)
    pos = {g: i for i, g in enumerate(keep_idx)}
    buckets = [[] for _ in keep_idx]  # list of (col, val) per kept gene
    with gzip.open(mtx_gz, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            break  # header dims already known
        for line in f:
            r, c, v = line.split()
            gi = int(r) - 1
            if gi in keep_set:
                buckets[pos[gi]].append((int(c) - 1, float(v)))
    data, indices, indptr = [], [], [0]
    # build genes x cells then transpose — or build cells x genes via COO
    rows, cols, vals = [], [], []
    for new_j, pairs in enumerate(buckets):
        for ci, v in pairs:
            rows.append(ci)
            cols.append(new_j)
            vals.append(v)
    X = sp.coo_matrix((vals, (rows, cols)), shape=(n_cells, len(keep_idx))).tocsr()
    X.data = X.data.astype(np.float32)
    return X


def analyze_gse325414():
    print("==== GSE325414 ====")
    feat_path = f"{DATA}/GSE325414/GSE325414_features.tsv.gz"
    bc_path = f"{DATA}/GSE325414/GSE325414_barcodes.tsv.gz"
    meta_path = f"{DATA}/GSE325414/GSE325414_metadata_individual_cells.txt.gz"
    mtx_path = f"{DATA}/GSE325414/GSE325414_matrix.mtx.gz"

    genes = [l.strip().split("\t")[0] for l in gzip.open(feat_path, "rt")]
    barcodes = [l.strip().split("\t")[0] for l in gzip.open(bc_path, "rt")]
    want = ["TACSTD2", "CLDN4", "EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC",
            "CD3D", "CD3E", "CD8A", "NKG7", "GNLY"]
    keep_names = [g for g in want if g in genes]
    keep_idx = [genes.index(g) for g in keep_names]
    print(f"  streaming {len(keep_names)} genes from {len(genes)} x {len(barcodes)} mtx")
    X = stream_gene_subset(mtx_path, len(genes), len(barcodes), keep_idx)
    ad = sc.AnnData(X)
    ad.var_names = keep_names
    ad.obs_names = barcodes

    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("cell")
    # align
    common = ad.obs_names.intersection(meta.index)
    ad = ad[common].copy()
    meta = meta.loc[ad.obs_names]
    ad.obs["donor"] = meta["donor"].values
    ad.obs["sample"] = meta["sample"].values
    ad.obs["major"] = meta["major.populations"].astype(str).values
    ad.obs["sub1"] = meta["sub.pop.level1"].astype(str).values
    ad.obs["sample_type"] = meta["sample.type.3"].astype(str).values
    ad.obs["mTLS"] = meta["mTLS"].astype(str).values

    # lineage from author labels (prefer major, fall back to sub1)
    lin = np.array(["Other"] * ad.n_obs, dtype=object)
    lin[np.isin(ad.obs["sub1"].values, list(MAL_SUB))] = "Malignant/Epithelial"
    lin[np.isin(ad.obs["major"].values, list(TNK_MAJOR))] = "T/NK"
    lin[(ad.obs["major"].values == "NA") & np.isin(ad.obs["sub1"].values, list(TNK_SUB))] = "T/NK"
    ad.obs["lineage"] = lin

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    for g in ("TACSTD2", "CLDN4"):
        ad.obs[g] = np.asarray(ad[:, g].X.todense()).ravel()

    # per-sample (author sample IDs)
    rows = []
    for s, sub in ad.obs.groupby("sample", observed=True):
        mal = sub[sub["lineage"] == "Malignant/Epithelial"]
        rows.append({
            "sample": s,
            "donor": sub["donor"].iloc[0],
            "sample_type": sub["sample_type"].iloc[0],
            "mTLS": sub["mTLS"].iloc[0],
            "n_cells": int(len(sub)),
            "n_malignant": int(len(mal)),
            "n_TNK": int((sub["lineage"] == "T/NK").sum()),
            "frac_TNK": float((sub["lineage"] == "T/NK").mean()),
            "malignant_TACSTD2_mean": float(mal["TACSTD2"].mean()) if len(mal) else np.nan,
            "malignant_CLDN4_mean": float(mal["CLDN4"].mean()) if len(mal) else np.nan,
            "malignant_TACSTD2_pct_pos": float((mal["TACSTD2"] > 0).mean() * 100) if len(mal) else np.nan,
            "malignant_CLDN4_pct_pos": float((mal["CLDN4"] > 0).mean() * 100) if len(mal) else np.nan,
        })
    ps = pd.DataFrame(rows).sort_values(["donor", "sample"])
    ps.to_csv(f"{OUT}/GSE325414_per_sample_summary.csv", index=False)

    # donor-level (collapse samples) for a cleaner T/NK correlation
    don = (ps.groupby("donor")
           .apply(lambda g: pd.Series({
               "n_cells": g["n_cells"].sum(),
               "n_malignant": g["n_malignant"].sum(),
               "n_TNK": g["n_TNK"].sum(),
               "frac_TNK": g["n_TNK"].sum() / g["n_cells"].sum(),
               "malignant_TACSTD2_mean": np.nanmean(g["malignant_TACSTD2_mean"]),
               "malignant_CLDN4_mean": np.nanmean(g["malignant_CLDN4_mean"]),
           }), include_groups=False)
           .reset_index())
    don.to_csv(f"{OUT}/GSE325414_per_donor_summary.csv", index=False)

    ok = ps.dropna(subset=["malignant_TACSTD2_mean"])
    res = {
        "cohort": "GSE325414",
        "pdat": "2026/06/10",
        "already_listed_as_candidate": False,
        "note": "NEW 2026 NSCLC scRNA; PEF ablation, not ICI. No response labels. TACSTD2/CLDN4 vs T/NK only.",
        "n_cells": int(ad.n_obs),
        "n_samples": int(ad.obs["sample"].nunique()),
        "n_donors": int(ad.obs["donor"].nunique()),
        "n_malignant": int((ad.obs["lineage"] == "Malignant/Epithelial").sum()),
        "n_TNK": int((ad.obs["lineage"] == "T/NK").sum()),
        "lineage_counts": ad.obs["lineage"].value_counts().to_dict(),
        "spearman_sample_TACSTD2_vs_fracTNK": spear(ok["malignant_TACSTD2_mean"], ok["frac_TNK"]),
        "spearman_sample_CLDN4_vs_fracTNK": spear(ok["malignant_CLDN4_mean"], ok["frac_TNK"]),
        "spearman_donor_TACSTD2_vs_fracTNK": spear(don["malignant_TACSTD2_mean"], don["frac_TNK"]),
        "spearman_donor_CLDN4_vs_fracTNK": spear(don["malignant_CLDN4_mean"], don["frac_TNK"]),
    }
    # TLS High vs Low among samples with malignant cells
    hi = ok[ok["mTLS"] == "High"]["malignant_TACSTD2_mean"]
    lo = ok[ok["mTLS"] == "Low"]["malignant_TACSTD2_mean"]
    res["TACSTD2_mTLS_High_vs_Low"] = mwu(hi, lo)
    res["CLDN4_mTLS_High_vs_Low"] = mwu(
        ok[ok["mTLS"] == "High"]["malignant_CLDN4_mean"],
        ok[ok["mTLS"] == "Low"]["malignant_CLDN4_mean"])
    json.dump(res, open(f"{OUT}/GSE325414_results.json", "w"), indent=2)
    print(json.dumps(res, indent=2))
    return ps, don, res


def make_plots(ps291, ps325):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    ax = axes[0, 0]
    for status, col, x in [("Responder", "#2c7fb8", 0), ("Non-responder", "#de2d26", 1)]:
        vals = ps291[ps291["responder"] == status]["malignant_TACSTD2_mean"].dropna().values
        ax.scatter(np.full(len(vals), x) + np.random.normal(0, 0.04, len(vals)),
                   vals, color=col, s=70, edgecolor="k", zorder=3, label=status)
        if len(vals):
            ax.hlines(vals.mean(), x - 0.2, x + 0.2, color=col, lw=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MPR", "non-MPR"])
    ax.set_ylabel("Mean malignant TACSTD2 (log-norm)")
    ax.set_title("GSE291670: malignant TACSTD2 vs MPR")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for status, col, x in [("Responder", "#2c7fb8", 0), ("Non-responder", "#de2d26", 1)]:
        vals = ps291[ps291["responder"] == status]["malignant_CLDN4_mean"].dropna().values
        ax.scatter(np.full(len(vals), x) + np.random.normal(0, 0.04, len(vals)),
                   vals, color=col, s=70, edgecolor="k", zorder=3)
        if len(vals):
            ax.hlines(vals.mean(), x - 0.2, x + 0.2, color=col, lw=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MPR", "non-MPR"])
    ax.set_ylabel("Mean malignant CLDN4 (log-norm)")
    ax.set_title("GSE291670: malignant CLDN4 vs MPR")

    ax = axes[1, 0]
    ax.scatter(ps291["frac_TNK"], ps291["malignant_TACSTD2_mean"], s=70, edgecolor="k",
               c=ps291["responder"].map({"Responder": "#2c7fb8", "Non-responder": "#de2d26"}))
    ax.set_xlabel("T/NK fraction")
    ax.set_ylabel("Mean malignant TACSTD2")
    ax.set_title("GSE291670: TACSTD2 vs T/NK")

    ax = axes[1, 1]
    ok = ps325.dropna(subset=["malignant_TACSTD2_mean"])
    ax.scatter(ok["frac_TNK"], ok["malignant_TACSTD2_mean"], s=28, alpha=0.75, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("T/NK fraction (per sample)")
    ax.set_ylabel("Mean malignant TACSTD2")
    ax.set_title("GSE325414 (2026 NSCLC, not ICI): TACSTD2 vs T/NK")

    plt.tight_layout()
    plt.savefig(f"{OUT}/tacstd2_cldn4_summary.png", dpi=150)
    print(f"wrote {OUT}/tacstd2_cldn4_summary.png")


def main():
    ps291, res291 = analyze_gse291670()
    ps325, don325, res325 = analyze_gse325414()
    make_plots(ps291, ps325)
    json.dump([res291, res325], open(f"{OUT}/combined_results.json", "w"), indent=2)


if __name__ == "__main__":
    main()
