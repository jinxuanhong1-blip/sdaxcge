#!/usr/bin/env python3
"""
GSE329813 GeoMx DSP (not scRNA): neoadjuvant pembrolizumab + chemo, MPR primary endpoint.

Public processed normalized matrix (genes x ROIs). ROI-to-patient/MPR map is parsed
from GEO GSM sample titles ('ROI N, Patient P, <site>, MPR|non-MPR|...').

Finding A: per-patient mean TACSTD2 in primary-tumor-bed ROIs, NMPR vs MPR.
Finding B: per-patient TACSTD2 vs T/NK signature (same marker mean) on those ROIs.

This is a spatial-bulk analog, labeled as such. It is not single-cell malignant isolation.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

META = Path("results/hunt_neoadj/metadata/GSE329813.gsm.soft.txt")
MAT = Path("results/hunt_neoadj/data/GSE329813_processed_data_file_normalized_data.csv.gz")
OUT = Path("results/hunt_neoadj/tables")

TNK = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1", "GZMB", "PRF1"]


def parse_gsm_titles() -> pd.DataFrame:
    titles = []
    cur = None
    for line in META.read_text().splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split(" = ", 1)[1] if " = " in line else line}
        elif line.startswith("!Sample_title = ") and cur is not None:
            cur["title"] = line.split(" = ", 1)[1]
            titles.append(cur)
            cur = None
    rows = []
    for t in titles:
        title = t["title"]
        # e.g. ROI 1, Patient 1, Primary tumor bed, MPR
        m = re.match(r"ROI\s+(\d+),\s+Patient\s+(\d+),\s+(.+),\s+(\S+)\s*$", title)
        if not m:
            rows.append({"gsm": t["gsm"], "title": title, "parse_ok": False})
            continue
        roi, patient, site, resp = m.groups()
        grp = np.nan
        if resp.upper() in {"NMPR", "NON-MPR", "NONMPR"} or resp.lower().startswith("non"):
            grp = "NMPR"
        elif resp.upper() in {"MPR", "PCR"} or resp.lower() == "pcr":
            grp = "MPR"
        rows.append(
            {
                "gsm": t["gsm"],
                "title": title,
                "roi": f"ROI {roi}",
                "patient": f"P{patient}",
                "site": site,
                "response_raw": resp,
                "group": grp,
                "parse_ok": True,
            }
        )
    return pd.DataFrame(rows)


def main():
    roi_meta = parse_gsm_titles()
    roi_meta.to_csv(OUT / "GSE329813_roi_metadata.csv", index=False)
    print("parsed", roi_meta.parse_ok.sum(), "/", len(roi_meta), "titles")
    print(roi_meta.group.value_counts(dropna=False).to_string())
    print(roi_meta.site.value_counts(dropna=False).head(10).to_string())

    expr = pd.read_csv(MAT, index_col=0)
    expr.index = expr.index.astype(str)
    print("matrix", expr.shape, "TACSTD2" in expr.index)

    genes_present = [g for g in ["TACSTD2"] + TNK if g in expr.index]
    missing = [g for g in ["TACSTD2"] + TNK if g not in expr.index]
    print("present", genes_present, "missing", missing)

    # restrict to ROIs present in both
    common = [c for c in expr.columns if c in set(roi_meta.roi)]
    print("ROIs in both", len(common), "matrix cols", expr.shape[1], "meta rois", roi_meta.roi.nunique())

    md = roi_meta.set_index("roi")
    tumor = md[md["site"].str.contains("Primary tumor", case=False, na=False)].copy()
    tumor = tumor[tumor.index.isin(expr.columns)]

    per_roi = tumor.copy()
    per_roi["TACSTD2"] = expr.loc["TACSTD2", per_roi.index].astype(float).values if "TACSTD2" in expr.index else np.nan
    tnk_use = [g for g in TNK if g in expr.index]
    per_roi["TNK"] = expr.loc[tnk_use, per_roi.index].astype(float).mean(axis=0).values
    per_roi.to_csv(OUT / "GSE329813_tumorbed_per_roi.csv")

    # patient-level mean over tumor-bed ROIs
    rows = []
    for pid, sub in per_roi.groupby("patient"):
        rows.append(
            {
                "patient": pid,
                "n_rois": int(len(sub)),
                "group": sub["group"].iloc[0],
                "TACSTD2": float(sub["TACSTD2"].mean()),
                "TNK": float(sub["TNK"].mean()),
            }
        )
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "GSE329813_tumorbed_per_patient.csv", index=False)

    a = per[per.group.isin(["MPR", "NMPR"])]
    nmpr = a.loc[a.group == "NMPR", "TACSTD2"].values
    mpr = a.loc[a.group == "MPR", "TACSTD2"].values
    findingA = {"n_NMPR": int(len(nmpr)), "n_MPR": int(len(mpr)), "n_patients": int(len(a))}
    if len(nmpr) >= 2 and len(mpr) >= 2:
        u2, p2 = stats.mannwhitneyu(nmpr, mpr, alternative="two-sided")
        _, p1 = stats.mannwhitneyu(nmpr, mpr, alternative="greater")
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
            }
        )
    if len(a) >= 4:
        rho, pB = stats.spearmanr(a["TACSTD2"], a["TNK"])
        r, pP = stats.pearsonr(a["TACSTD2"], a["TNK"])
        findingB = {
            "n_patients": int(len(a)),
            "tnk_markers_used": tnk_use,
            "spearman_rho": float(rho),
            "spearman_p": float(pB),
            "pearson_r": float(r),
            "pearson_p": float(pP),
            "direction_matches_expected(negative)": bool(rho < 0),
        }
    else:
        findingB = {"n_patients": int(len(a))}

    result = {
        "dataset": "GSE329813",
        "assay": "GeoMx DSP normalized counts (NOT scRNA); primary tumor-bed ROIs aggregated per patient",
        "genes_missing": missing,
        "findingA": findingA,
        "findingB": findingB,
    }
    (OUT / "GSE329813_results.json").write_text(json.dumps(result, indent=2))
    print(per.to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
