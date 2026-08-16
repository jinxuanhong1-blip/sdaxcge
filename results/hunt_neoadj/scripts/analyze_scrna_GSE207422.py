#!/usr/bin/env python3
"""
GSE207422 scRNA-seq: marker-based malignant vs T/NK annotation on the public UMI matrix.

Does NOT download FASTQ. Streams the GEO supplementary dense UMI TSV and keeps
only a small lineage-marker gene set (plus TACSTD2).

Finding A: among epithelial/malignant-like cells, per-sample mean TACSTD2
           (log1p UMI) is higher in NMPR than MPR.
Finding B: per-patient (one post-tx surgery sample preferred; else the only
           available sample) malignant TACSTD2 vs T/NK cell fraction, Spearman rho.

Annotation (deliberately simple, fully specified):
  epithelial/malignant-like: (EPCAM>0 or KRT19>0 or KRT18>0 or KRT8>0)
                             AND PTPRC==0
  T/NK: (CD3D>0 or CD3E>0 or NKG7>0 or GNLY>0) AND PTPRC>0
  discarded from both if they meet both rules (should be rare).

pCR is grouped with MPR. NE is excluded from Finding A.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("results/hunt_neoadj/data")
OUT = Path("results/hunt_neoadj/tables")
SRC = DATA / "GSE207422_scrna_markers.tsv.gz"
META = DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"

MARKER_GENES = [
    "TACSTD2",
    "EPCAM",
    "KRT19",
    "KRT18",
    "KRT8",
    "KRT7",
    "KRT5",
    "KRT17",
    "NAPSA",
    "NKX2-1",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD3G",
    "CD2",
    "CD8A",
    "CD4",
    "IL7R",
    "NKG7",
    "GNLY",
    "KLRD1",
    "NCAM1",
    "GZMB",
    "PRF1",
    "MS4A1",
    "CD79A",
    "LYZ",
    "CD68",
    "PECAM1",
    "COL1A1",
]


def load_marker_matrix() -> pd.DataFrame:
    """Return cells x genes UMI counts from the pre-extracted marker TSV."""
    raw = pd.read_csv(SRC, sep="\t", index_col=0)
    raw = raw[~raw.index.duplicated(keep="first")]
    missing = [g for g in MARKER_GENES if g not in raw.index]
    print("loaded", list(raw.index), "missing", missing, "n_cells", raw.shape[1], flush=True)
    df = raw.T.astype(np.float32)
    df.index.name = "cell"
    return df


def response_group(x) -> str | float:
    x = str(x)
    if x.startswith("NMPR"):
        return "NMPR"
    if x.startswith("MPR") or x == "pCR" or "pCR" in x:
        return "MPR"
    return np.nan


def main():
    expr = load_marker_matrix()
    expr["sample"] = [c.rsplit("_", 1)[0] for c in expr.index]

    epi = (
        ((expr.get("EPCAM", 0) > 0) | (expr.get("KRT19", 0) > 0) | (expr.get("KRT18", 0) > 0) | (expr.get("KRT8", 0) > 0))
        & (expr.get("PTPRC", 0) == 0)
    )
    tnk = (
        ((expr.get("CD3D", 0) > 0) | (expr.get("CD3E", 0) > 0) | (expr.get("NKG7", 0) > 0) | (expr.get("GNLY", 0) > 0))
        & (expr.get("PTPRC", 0) > 0)
    )
    both = epi & tnk
    epi = epi & ~both
    tnk = tnk & ~both
    expr["is_epi"] = epi.astype(int)
    expr["is_tnk"] = tnk.astype(int)
    expr["tac_log1p"] = np.log1p(expr["TACSTD2"].astype(float))

    meta = pd.read_excel(META)
    meta = meta.dropna(subset=["Sample"]).copy()
    meta["group"] = meta["Pathologic Response"].map(response_group)
    meta["Resource"] = meta["Resource"].astype(str)
    meta = meta.set_index("Sample")

    # per-sample summaries
    rows = []
    for sample, sub in expr.groupby("sample"):
        n = len(sub)
        n_epi = int(sub["is_epi"].sum())
        n_tnk = int(sub["is_tnk"].sum())
        if n_epi > 0:
            tac_mean = float(sub.loc[sub["is_epi"] == 1, "tac_log1p"].mean())
            tac_mean_raw = float(sub.loc[sub["is_epi"] == 1, "TACSTD2"].mean())
            frac_tac_pos = float((sub.loc[sub["is_epi"] == 1, "TACSTD2"] > 0).mean())
        else:
            tac_mean = np.nan
            tac_mean_raw = np.nan
            frac_tac_pos = np.nan
        rec = {
            "sample": sample,
            "n_cells": n,
            "n_epi": n_epi,
            "n_tnk": n_tnk,
            "frac_epi": n_epi / n if n else np.nan,
            "frac_tnk": n_tnk / n if n else np.nan,
            "malignant_TACSTD2_mean_log1p": tac_mean,
            "malignant_TACSTD2_mean_UMI": tac_mean_raw,
            "malignant_TACSTD2_frac_pos": frac_tac_pos,
        }
        if sample in meta.index:
            rec.update(
                {
                    "patient": meta.loc[sample, "Patient"],
                    "resource": meta.loc[sample, "Resource"],
                    "pathologic_response": meta.loc[sample, "Pathologic Response"],
                    "group": meta.loc[sample, "group"],
                    "residual_tumor": meta.loc[sample, "Residual Tumor"],
                    "pathology": meta.loc[sample, "Pathology"],
                }
            )
        rows.append(rec)
    per = pd.DataFrame(rows).sort_values("sample")
    per.to_csv(OUT / "GSE207422_scrna_per_sample.csv", index=False)

    # Finding A: samples with >=20 epithelial cells and a defined MPR/NMPR group
    a = per[(per["n_epi"] >= 20) & (per["group"].isin(["MPR", "NMPR"]))].copy()
    # Prefer post-treatment surgery when both exist; here samples are 1:1 with patients
    nmpr = a.loc[a.group == "NMPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    mpr = a.loc[a.group == "MPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    findingA = {"n_samples_used": int(len(a)), "n_NMPR": int(len(nmpr)), "n_MPR": int(len(mpr))}
    if len(nmpr) >= 2 and len(mpr) >= 2:
        u2, p2 = stats.mannwhitneyu(nmpr, mpr, alternative="two-sided")
        u1, p1 = stats.mannwhitneyu(nmpr, mpr, alternative="greater")
        findingA.update(
            {
                "median_NMPR": float(np.median(nmpr)),
                "median_MPR": float(np.median(mpr)),
                "mean_NMPR": float(np.mean(nmpr)),
                "mean_MPR": float(np.mean(mpr)),
                "mannwhitney_U": float(u2),
                "p_two_sided": float(p2),
                "p_one_sided_NMPR_gt_MPR": float(p1),
                "direction_matches_expected(NMPR>MPR)": bool(np.median(nmpr) > np.median(mpr)),
                "samples_used": a["sample"].tolist(),
            }
        )
    else:
        findingA["note"] = "insufficient groups after epithelial-cell filter"

    # Finding A sensitivity: post-treatment surgery only
    a_post = a[a["resource"].str.contains("Post", case=False, na=False)]
    nmpr_p = a_post.loc[a_post.group == "NMPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    mpr_p = a_post.loc[a_post.group == "MPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    findingA_post = {"n_NMPR": int(len(nmpr_p)), "n_MPR": int(len(mpr_p))}
    if len(nmpr_p) >= 2 and len(mpr_p) >= 2:
        u2, p2 = stats.mannwhitneyu(nmpr_p, mpr_p, alternative="two-sided")
        u1, p1 = stats.mannwhitneyu(nmpr_p, mpr_p, alternative="greater")
        findingA_post.update(
            {
                "median_NMPR": float(np.median(nmpr_p)),
                "median_MPR": float(np.median(mpr_p)),
                "p_two_sided": float(p2),
                "p_one_sided_NMPR_gt_MPR": float(p1),
                "direction_matches_expected(NMPR>MPR)": bool(np.median(nmpr_p) > np.median(mpr_p)),
            }
        )

    # Finding B: one row per patient. Prefer post-treatment surgery sample.
    bsrc = per[(per["n_epi"] >= 20) & (per["n_tnk"] >= 20) & per["patient"].notna()].copy()
    bsrc["pref"] = bsrc["resource"].str.contains("Post", case=False, na=False).astype(int)
    bsrc = bsrc.sort_values(["patient", "pref"], ascending=[True, False])
    bsrc = bsrc.drop_duplicates("patient", keep="first")
    if len(bsrc) >= 4:
        rho, pB = stats.spearmanr(
            bsrc["malignant_TACSTD2_mean_log1p"].astype(float),
            bsrc["frac_tnk"].astype(float),
        )
        r_pear, p_pear = stats.pearsonr(
            bsrc["malignant_TACSTD2_mean_log1p"].astype(float),
            bsrc["frac_tnk"].astype(float),
        )
        findingB = {
            "n_patients": int(len(bsrc)),
            "spearman_rho": float(rho),
            "spearman_p": float(pB),
            "pearson_r": float(r_pear),
            "pearson_p": float(p_pear),
            "direction_matches_expected(negative)": bool(rho < 0),
            "patients": bsrc["patient"].tolist(),
        }
    else:
        findingB = {"n_patients": int(len(bsrc)), "note": "too few patients"}

    result = {
        "dataset": "GSE207422",
        "assay": "scRNA-seq public UMI matrix; marker-based epithelial (PTPRC-) vs T/NK (PTPRC+)",
        "n_cells_total": int(len(expr)),
        "n_epi": int(expr["is_epi"].sum()),
        "n_tnk": int(expr["is_tnk"].sum()),
        "n_both_excluded": int(both.sum()),
        "min_epi_cells_per_sample": 20,
        "findingA_all_evaluable_samples": findingA,
        "findingA_post_treatment_only": findingA_post,
        "findingB_per_patient_TACSTD2_vs_TNK_frac": findingB,
    }
    (OUT / "GSE207422_scrna_results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
