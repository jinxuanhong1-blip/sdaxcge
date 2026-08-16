#!/usr/bin/env python3
"""
GSE207422 bulk RNA-seq (pre-treatment NSCLC neoadjuvant PD-1 + chemo).
Tests two target findings:
  A) TACSTD2 higher in NMPR than MPR (Mann-Whitney U, one-sided NMPR>MPR reported alongside two-sided).
  B) Per-patient TACSTD2 vs T/NK signature Spearman rho is negative.

Bulk here is a proxy for the malignant compartment: TACSTD2 (TROP2) is an
epithelial/tumor marker, so in bulk tumor it primarily reflects malignant cells.
This is an honest, conservative check; the scRNA analysis isolates malignant cells directly.
"""
import gzip
import json
import numpy as np
import pandas as pd
from scipy import stats

DATA = "results/hunt_neoadj/data"
OUT = "results/hunt_neoadj/tables"

# Canonical T and NK lineage markers (present-in-matrix subset used).
TNK_MARKERS = [
    "CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "TRBC2",
    "CD8A", "CD8B", "CD4", "IL7R",
    "NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "GZMB", "GZMK", "PRF1",
]

def load_matrix():
    df = pd.read_csv(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
                     sep="\t", index_col=0)
    df = df[~df.index.duplicated(keep="first")]
    return df

def main():
    expr = load_matrix()  # genes x samples, log2 TPM
    meta = pd.read_excel(f"{DATA}/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"])
    meta = meta[meta["Sample"].isin(expr.columns)].copy()

    # Map pathologic response to MPR vs NMPR (MPR group includes pCR / MPR(pCR)).
    def grp(x):
        x = str(x)
        if x.startswith("NMPR"):
            return "NMPR"
        if x.startswith("MPR") or "pCR" in x:
            return "MPR"
        return np.nan
    meta["group"] = meta["Pathologic Response"].map(grp)
    meta = meta.dropna(subset=["group"])

    samples = meta["Sample"].tolist()
    expr = expr[samples]

    # ---- Finding A: TACSTD2 NMPR vs MPR ----
    tac = expr.loc["TACSTD2"]
    d = meta.set_index("Sample").copy()
    d["TACSTD2"] = tac.reindex(d.index).values
    nmpr = d.loc[d.group == "NMPR", "TACSTD2"].astype(float).values
    mpr = d.loc[d.group == "MPR", "TACSTD2"].astype(float).values

    u2, p2 = stats.mannwhitneyu(nmpr, mpr, alternative="two-sided")
    u1, p1 = stats.mannwhitneyu(nmpr, mpr, alternative="greater")  # NMPR > MPR
    # rank-biserial effect size
    n1, n2 = len(nmpr), len(mpr)
    rbc = 1 - 2 * u2 / (n1 * n2) if (n1 and n2) else np.nan
    rbc = -rbc  # orient so positive = NMPR>MPR (since U computed with nmpr first)
    findingA = {
        "n_NMPR": int(n1), "n_MPR": int(n2),
        "median_TACSTD2_NMPR": float(np.median(nmpr)),
        "median_TACSTD2_MPR": float(np.median(mpr)),
        "mean_TACSTD2_NMPR": float(np.mean(nmpr)),
        "mean_TACSTD2_MPR": float(np.mean(mpr)),
        "mannwhitney_U": float(u2),
        "p_two_sided": float(p2),
        "p_one_sided_NMPR_gt_MPR": float(p1),
        "rank_biserial_NMPR_vs_MPR": float(rbc),
        "direction_matches_expected(NMPR>MPR)": bool(np.median(nmpr) > np.median(mpr)),
    }

    # ---- Finding B: per-patient TACSTD2 vs T/NK signature (Spearman) ----
    present = [g for g in TNK_MARKERS if g in expr.index]
    missing = [g for g in TNK_MARKERS if g not in expr.index]
    tnk_score = expr.loc[present].mean(axis=0)  # per-sample mean log2TPM
    tac_all = expr.loc["TACSTD2"]
    rho, pB = stats.spearmanr(tac_all.values, tnk_score.reindex(tac_all.index).values)
    r_pear, pB_pear = stats.pearsonr(tac_all.values, tnk_score.reindex(tac_all.index).values)
    findingB = {
        "n_samples": int(expr.shape[1]),
        "n_patients": int(meta["Patient"].nunique()),
        "tnk_markers_used": present,
        "tnk_markers_missing": missing,
        "spearman_rho": float(rho),
        "spearman_p": float(pB),
        "pearson_r": float(r_pear),
        "pearson_p": float(pB_pear),
        "direction_matches_expected(negative)": bool(rho < 0),
    }

    # Per-sample export
    tab = d[["Patient", "group", "TACSTD2"]].copy()
    tab["TNK_score"] = tnk_score.reindex(tab.index).values
    tab.to_csv(f"{OUT}/GSE207422_bulk_per_sample.csv")

    result = {
        "dataset": "GSE207422",
        "assay": "bulk RNA-seq (log2 TPM), all pre-treatment biopsies, 1 sample/patient",
        "findingA_TACSTD2_NMPR_vs_MPR": findingA,
        "findingB_TACSTD2_vs_TNK_spearman": findingB,
    }
    with open(f"{OUT}/GSE207422_bulk_results.json", "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
