#!/usr/bin/env python3
"""
GSE248378 bulk FPKM: post-neoadjuvant durvalumab +/- SBRT resected NSCLC.

GEO sample characteristics have treatment arm (Arm1/Arm2) and histology, NOT MPR/NMPR.
Finding A cannot be tested from public metadata.

Finding B: per-sample TACSTD2 vs T/NK signature Spearman (bulk proxy).
Sample titles with a trailing R are NOT treated as MPR labels.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("results/hunt_neoadj/data")
META = Path("results/hunt_neoadj/metadata/GSE248378.gsm.soft.txt")
OUT = Path("results/hunt_neoadj/tables")
TNK = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4", "IL7R",
       "NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "GZMB", "GZMK", "PRF1"]


def parse_gsm() -> pd.DataFrame:
    recs = []
    cur = {}
    for line in META.read_text().splitlines():
        if line.startswith("^SAMPLE"):
            if cur.get("title"):
                recs.append(cur)
            cur = {"gsm": line.split(" = ", 1)[1] if " = " in line else ""}
        elif line.startswith("!Sample_title = "):
            cur["title"] = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = "):
            kv = line.split(" = ", 1)[1]
            if ": " in kv:
                k, v = kv.split(": ", 1)
                cur[k] = v
    if cur.get("title"):
        recs.append(cur)
    return pd.DataFrame(recs)


def main():
    expr = pd.read_csv(DATA / "GSE248378_Durva_Post_FPKMs.txt.gz", sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    meta = parse_gsm()
    meta.to_csv(OUT / "GSE248378_sample_metadata.csv", index=False)

    present = [g for g in TNK if g in expr.index]
    missing = [g for g in TNK if g not in expr.index]
    if "TACSTD2" not in expr.index:
        result = {"dataset": "GSE248378", "error": "TACSTD2 absent from FPKM matrix"}
        (OUT / "GSE248378_bulk_results.json").write_text(json.dumps(result, indent=2))
        print(result)
        return

    samples = [c for c in expr.columns]
    tac = expr.loc["TACSTD2", samples].astype(float)
    tnk = expr.loc[present, samples].astype(float).mean(axis=0)
    rho, pB = stats.spearmanr(tac, tnk)
    r, pP = stats.pearsonr(tac, tnk)
    tab = pd.DataFrame({"sample": samples, "TACSTD2": tac.values, "TNK": tnk.values})
    tab = tab.merge(meta, left_on="sample", right_on="title", how="left")
    tab.to_csv(OUT / "GSE248378_bulk_per_sample.csv", index=False)

    result = {
        "dataset": "GSE248378",
        "assay": "bulk FPKM, post-neoadjuvant durvalumab +/- SBRT resected tumors",
        "findingA": {
            "tested": False,
            "reason": "GEO GSM characteristics have treatment arm and histology only; no MPR/NMPR/pCR labels",
        },
        "findingB": {
            "n_samples": int(len(samples)),
            "tnk_markers_used": present,
            "tnk_markers_missing": missing,
            "spearman_rho": float(rho),
            "spearman_p": float(pB),
            "pearson_r": float(r),
            "pearson_p": float(pP),
            "direction_matches_expected(negative)": bool(rho < 0),
        },
    }
    (OUT / "GSE248378_bulk_results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
