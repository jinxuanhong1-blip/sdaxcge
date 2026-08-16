#!/usr/bin/env python3
"""
GSE291670: 6 NSCLC tumors after neoadjuvant camrelizumab + anlotinib.
Public 10x MTX in GSE291670_RAW.tar. Sample titles are MPR-1/2/3 and Non-MPR-1/2/3.

Marker-based epithelial (PTPRC-) vs T/NK (PTPRC+). Same rules as GSE207422.
n=3 vs 3: stats are reported but underpowered; direction is the honest takeaway.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("results/hunt_neoadj/data/GSE291670")
OUT = Path("results/hunt_neoadj/tables")

MARKERS = [
    "TACSTD2", "EPCAM", "KRT19", "KRT18", "KRT8", "PTPRC",
    "CD3D", "CD3E", "NKG7", "GNLY",
]


def gene_rows(features_path: Path) -> dict[str, int]:
    out = {}
    with gzip.open(features_path, "rt") as f:
        for i, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            for g in MARKERS:
                if g in parts[:2]:
                    out[g] = i
    return out


def extract_genes(mtx_path: Path, rows: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
    want = {v: k for k, v in rows.items()}
    vals = {g: np.zeros(n_cells, dtype=np.float32) for g in rows}
    with gzip.open(mtx_path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            break
        for line in f:
            r, c, v = line.split()
            ri = int(r)
            if ri in want:
                vals[want[ri]][int(c) - 1] = float(v)
    return vals


def main():
    samples = []
    for mtx in sorted(DATA.glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        # GSM8839599_MPR-1
        label = stem.split("_", 1)[1]  # MPR-1
        group = "MPR" if label.startswith("MPR") else "NMPR"
        feat = DATA / f"{stem}_features.tsv.gz"
        bc = DATA / f"{stem}_barcodes.tsv.gz"
        n_cells = sum(1 for _ in gzip.open(bc, "rt"))
        rows = gene_rows(feat)
        print(label, "n_cells", n_cells, "genes_found", rows, flush=True)
        expr = extract_genes(mtx, rows, n_cells)
        df = pd.DataFrame(expr)
        epi = (
            ((df.get("EPCAM", 0) > 0) | (df.get("KRT19", 0) > 0) | (df.get("KRT18", 0) > 0) | (df.get("KRT8", 0) > 0))
            & (df.get("PTPRC", 0) == 0)
        )
        tnk = (
            ((df.get("CD3D", 0) > 0) | (df.get("CD3E", 0) > 0) | (df.get("NKG7", 0) > 0) | (df.get("GNLY", 0) > 0))
            & (df.get("PTPRC", 0) > 0)
        )
        both = epi & tnk
        epi = epi & ~both
        tnk = tnk & ~both
        n_epi = int(epi.sum())
        n_tnk = int(tnk.sum())
        if n_epi > 0:
            tac = float(np.log1p(df.loc[epi, "TACSTD2"]).mean())
            tac_frac = float((df.loc[epi, "TACSTD2"] > 0).mean())
        else:
            tac = np.nan
            tac_frac = np.nan
        samples.append(
            {
                "sample": label,
                "group": group,
                "n_cells": n_cells,
                "n_epi": n_epi,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n_cells if n_cells else np.nan,
                "malignant_TACSTD2_mean_log1p": tac,
                "malignant_TACSTD2_frac_pos": tac_frac,
                "genes_found": ",".join(sorted(rows)),
            }
        )

    per = pd.DataFrame(samples)
    per.to_csv(OUT / "GSE291670_scrna_per_sample.csv", index=False)

    a = per[per["n_epi"] >= 20]
    nmpr = a.loc[a.group == "NMPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    mpr = a.loc[a.group == "MPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    findingA = {"n_NMPR": int(len(nmpr)), "n_MPR": int(len(mpr))}
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
    b = a[a["n_tnk"] >= 20]
    if len(b) >= 4:
        rho, pB = stats.spearmanr(b["malignant_TACSTD2_mean_log1p"], b["frac_tnk"])
        r, pP = stats.pearsonr(b["malignant_TACSTD2_mean_log1p"], b["frac_tnk"])
        findingB = {
            "n_patients": int(len(b)),
            "spearman_rho": float(rho),
            "spearman_p": float(pB),
            "pearson_r": float(r),
            "pearson_p": float(pP),
            "direction_matches_expected(negative)": bool(rho < 0),
        }
    else:
        findingB = {"n_patients": int(len(b)), "note": "too few"}

    result = {
        "dataset": "GSE291670",
        "assay": "scRNA 10x MTX from GSE291670_RAW.tar; marker-based Epi vs T/NK",
        "n_samples": int(len(per)),
        "findingA": findingA,
        "findingB": findingB,
        "note": "n=3 vs 3. Mann-Whitney p-values are reported but have almost no power.",
    }
    (OUT / "GSE291670_scrna_results.json").write_text(json.dumps(result, indent=2))
    print(per.to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
