#!/usr/bin/env python3
"""
GSE241934 scRNA (NEOTIDE/CTONG2104 + real-world neoadjuvant IO+chemo).

Public processed MTX + author cell-type and MPR annotations.
Finding A: epithelial TACSTD2 (mean log1p UMI) NMPR > MPR (pCR counted as MPR).
Finding B: per-patient epithelial TACSTD2 vs T/NK fraction, Spearman rho.

IIT = EGFR-mutant trial cohort (11 pts). REAL = wild-type real-world (34 pts).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("results/hunt_neoadj/data")
OUT = Path("results/hunt_neoadj/tables")

GENES = ["TACSTD2"]  # only need TACSTD2; T/NK fraction comes from author labels


def gene_index(features_path: Path, gene: str) -> int:
    """1-based MTX row index of gene symbol (column 1 or 0)."""
    with gzip.open(features_path, "rt") as f:
        for i, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            if gene in parts[:2]:
                return i
    raise KeyError(gene)


def extract_gene_from_mtx(mtx_path: Path, gene_row: int, n_cells: int) -> np.ndarray:
    """Stream MTX (genes x cells, 1-based) and return dense UMI vector for one gene."""
    vals = np.zeros(n_cells, dtype=np.float32)
    with gzip.open(mtx_path, "rt") as f:
        # skip comments / header
        for line in f:
            if line.startswith("%"):
                continue
            # first non-comment is dimensions
            _nr, _nc, _nz = line.split()
            break
        for line in f:
            r, c, v = line.split()
            if int(r) == gene_row:
                vals[int(c) - 1] = float(v)
    return vals


def load_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    return df


def group_of(x) -> str | float:
    x = str(x)
    if x in {"non-MPR", "NMPR", "nonMPR"}:
        return "NMPR"
    if x in {"MPR", "pCR", "MPR (pCR)"}:
        return "MPR"
    return np.nan


def analyze_cohort(name: str, mtx: Path, features: Path, barcodes: Path, meta: Path) -> dict:
    print(f"=== {name} ===", flush=True)
    md = load_meta(meta)
    bc = pd.read_csv(barcodes, sep="\t", header=None)
    print("meta", md.shape, "barcodes", bc.shape, "meta cols", list(md.columns)[:8], flush=True)

    # Align cells: prefer cellID in meta vs barcode strings
    if "cellID" in md.columns:
        md = md.set_index("cellID", drop=False)
    else:
        md = md.set_index(md.columns[0], drop=False)

    n_cells = len(bc)
    gi = gene_index(features, "TACSTD2")
    print("TACSTD2 mtx row", gi, "n_cells", n_cells, flush=True)
    tac = extract_gene_from_mtx(mtx, gi, n_cells)
    print("TACSTD2 nonzero", int((tac > 0).sum()), "max", float(tac.max()), flush=True)

    # barcode alignment
    barcodes_list = bc.iloc[:, 0].astype(str).tolist()
    # IIT barcodes look like AAACCTGAGCTATGCT-1; meta cellID like P343_AAACCTGAGCTATGCT-1
    if barcodes_list[0] not in md.index and "cellID" in md.columns:
        # try sample-prefixed
        if "sampleID" in md.columns:
            # build from mtx order: if meta is in same order as barcodes
            if len(md) == n_cells:
                print("aligning by row order (equal lengths)", flush=True)
                md = md.reset_index(drop=True)
                md["TACSTD2"] = tac
            else:
                raise RuntimeError(f"cannot align {name}: {len(md)} vs {n_cells}")
        else:
            raise RuntimeError(f"cannot align {name}")
    else:
        md = md.reindex(barcodes_list)
        md["TACSTD2"] = tac

    md["tac_log1p"] = np.log1p(md["TACSTD2"].astype(float))
    major = md["major.cell.type"] if "major.cell.type" in md.columns else md["major_cell_type"]
    md["is_epi"] = major.eq("Epi")
    md["is_tnk"] = major.isin(["T", "NK"])
    md["group"] = md["Pathological Response"].map(group_of)
    md["sampleID"] = md["sampleID"].astype(str)

    rows = []
    for sid, sub in md.groupby("sampleID"):
        n = len(sub)
        n_epi = int(sub["is_epi"].sum())
        n_tnk = int(sub["is_tnk"].sum())
        if n_epi > 0:
            tac_mean = float(sub.loc[sub["is_epi"], "tac_log1p"].mean())
            tac_frac = float((sub.loc[sub["is_epi"], "TACSTD2"] > 0).mean())
        else:
            tac_mean = np.nan
            tac_frac = np.nan
        rows.append(
            {
                "cohort": name,
                "sampleID": sid,
                "n_cells": n,
                "n_epi": n_epi,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n if n else np.nan,
                "malignant_TACSTD2_mean_log1p": tac_mean,
                "malignant_TACSTD2_frac_pos": tac_frac,
                "group": sub["group"].iloc[0],
                "pathologic_response": sub["Pathological Response"].iloc[0],
                "PD1": sub["PD1"].iloc[0] if "PD1" in sub.columns else None,
                "EGFR": sub["EGFR"].iloc[0] if "EGFR" in sub.columns else None,
                "Histology": sub["Histology"].iloc[0] if "Histology" in sub.columns else None,
            }
        )
    per = pd.DataFrame(rows)
    return per


def stats_from_per(per: pd.DataFrame, label: str, min_epi: int = 20) -> dict:
    a = per[(per["n_epi"] >= min_epi) & (per["group"].isin(["MPR", "NMPR"]))].copy()
    nmpr = a.loc[a.group == "NMPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    mpr = a.loc[a.group == "MPR", "malignant_TACSTD2_mean_log1p"].astype(float).values
    out = {"label": label, "min_epi": min_epi, "n_NMPR": int(len(nmpr)), "n_MPR": int(len(mpr))}
    if len(nmpr) >= 2 and len(mpr) >= 2:
        u2, p2 = stats.mannwhitneyu(nmpr, mpr, alternative="two-sided")
        _, p1 = stats.mannwhitneyu(nmpr, mpr, alternative="greater")
        out.update(
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
    b = a[a["n_tnk"] >= 20].copy()
    if len(b) >= 4:
        rho, pB = stats.spearmanr(b["malignant_TACSTD2_mean_log1p"], b["frac_tnk"])
        r, pP = stats.pearsonr(b["malignant_TACSTD2_mean_log1p"], b["frac_tnk"])
        out["findingB"] = {
            "n_patients": int(len(b)),
            "spearman_rho": float(rho),
            "spearman_p": float(pB),
            "pearson_r": float(r),
            "pearson_p": float(pP),
            "direction_matches_expected(negative)": bool(rho < 0),
        }
    else:
        out["findingB"] = {"n_patients": int(len(b)), "note": "too few"}
    return out


def main():
    iit = analyze_cohort(
        "IIT_EGFRmut",
        DATA / "GSE241934_IIT_Matrix.mtx.gz",
        DATA / "GSE241934_IIT_features.tsv.gz",
        DATA / "GSE241934_IIT_barcodes.tsv.gz",
        DATA / "GSE241934_IIT_Meta.txt.gz",
    )
    real = analyze_cohort(
        "REAL_WT",
        DATA / "GSE241934_Real_Matrix.mtx.gz",
        DATA / "GSE241934_RWC_features.tsv.gz",
        DATA / "GSE241934_RWC_barcodes.tsv.gz",
        DATA / "GSE241934_Real_Meta.txt.gz",
    )
    per = pd.concat([iit, real], ignore_index=True)
    per.to_csv(OUT / "GSE241934_scrna_per_sample.csv", index=False)

    result = {
        "dataset": "GSE241934",
        "assay": "scRNA public MTX; author Epi / T / NK labels; pCR grouped with MPR",
        "IIT": stats_from_per(iit, "IIT_EGFRmut"),
        "REAL": stats_from_per(real, "REAL_WT"),
        "COMBINED": stats_from_per(per, "IIT+REAL"),
        "IIT_minEpi5": stats_from_per(iit, "IIT_minEpi5", min_epi=5),
        "COMBINED_minEpi5": stats_from_per(per, "IIT+REAL_minEpi5", min_epi=5),
    }
    (OUT / "GSE241934_scrna_results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
