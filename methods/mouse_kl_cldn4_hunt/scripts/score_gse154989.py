#!/usr/bin/env python3
"""Mouse-level Cldn4 vs T/NK and epithelial IFN/MHC in GSE154989 (KP/K plate scRNA)."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import CLDN4, IFN_MHC, T_NK  # noqa: E402
from score_utils import association_table, mean_score, present_genes, zscore_rows  # noqa: E402

RAW = Path("/tmp/geo_dl")
OUT = Path("/workspace/analysis/GSE154989")


def biological_mouse(mouse_id: str) -> str:
    """Collapse tumor suffixes (…_m2_T4) to one mouse (…_m2)."""
    return re.sub(r"_T\d+$", "", str(mouse_id))


def genotype_from_id(mouse_id: str) -> str:
    s = str(mouse_id)
    if s.startswith("KP"):
        return "KP"
    if s.startswith("K_"):
        return "K"
    if s.startswith("T"):
        return "T_early"
    return "other"


def load_sparse_gene_matrix(wanted: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = pd.read_csv(RAW / "GSE154989_mmLungPlate_fQC_geneTable.csv.gz")
    smp = pd.read_csv(RAW / "GSE154989_mmLungPlate_fQC_smpTable.csv.gz")
    genes["geneSymbol"] = genes["geneSymbol"].astype(str)
    # first occurrence of each symbol
    symbol_to_row = {}
    for i, sym in enumerate(genes["geneSymbol"]):
        if sym not in symbol_to_row:
            symbol_to_row[sym] = i
    use_syms = [g for g in wanted if g in symbol_to_row]
    row_map = {symbol_to_row[s] + 1: s for s in use_syms}  # 1-based i
    want_i = np.array(list(row_map.keys()), dtype=float)

    with h5py.File(RAW / "GSE154989_mmLungPlate_fQC_dSp_normTPM.h5", "r") as f:
        i = f["i"][0]
        j = f["j"][0]
        v = f["v"][0]
    mask = np.isin(i, want_i)
    i, j, v = i[mask], j[mask], v[mask]
    mat = pd.DataFrame(0.0, index=use_syms, columns=smp["sampleID"].astype(str))
    for gi, cj, val in zip(i.astype(int), j.astype(int), v):
        mat.loc[row_map[gi], smp["sampleID"].iloc[cj - 1]] = float(val)
    return mat, smp, genes


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wanted = list(dict.fromkeys(CLDN4 + T_NK + IFN_MHC))
    expr, smp, genes = load_sparse_gene_matrix(wanted)
    present = {
        "Cldn4": present_genes(expr.index, CLDN4),
        "T_NK": present_genes(expr.index, T_NK),
        "IFN_MHC": present_genes(expr.index, IFN_MHC),
    }
    z = zscore_rows(expr)
    cell = smp.copy()
    cell["sampleID"] = cell["sampleID"].astype(str)
    cell = cell.set_index("sampleID")
    cell["cldn4"] = expr.loc["Cldn4"] if "Cldn4" in expr.index else np.nan
    cell["tnk_score"] = mean_score(z, present["T_NK"])
    cell["ifn_mhc_score"] = mean_score(z, present["IFN_MHC"])
    cell["bio_mouse"] = cell["mouseID"].map(biological_mouse)
    cell["genotype"] = cell["mouseID"].map(genotype_from_id)
    cell["n_genes_tnk"] = len(present["T_NK"])
    cell["n_genes_ifn_mhc"] = len(present["IFN_MHC"])

    # Mouse-level: mean over cells. Honest: one row per biological mouse.
    g = cell.groupby("bio_mouse", sort=False)
    mouse = pd.DataFrame(
        {
            "bio_mouse": g.size().index,
            "genotype": g["genotype"].first(),
            "n_cells": g.size().values,
            "n_tumor_ids": g["mouseID"].nunique().values,
            "timesimple": g["timesimple"].first(),
            "cldn4_mean": g["cldn4"].mean().values,
            "cldn4_pct_pos": g["cldn4"].apply(lambda s: float((s > 0).mean())).values,
            "tnk_score": g["tnk_score"].mean().values,
            "ifn_mhc_score": g["ifn_mhc_score"].mean().values,
            "tnk_genes_used": ",".join(present["T_NK"]),
            "ifn_mhc_genes_used": ",".join(present["IFN_MHC"]),
        }
    )
    # Require enough cells for a stable mouse mean.
    mouse["pass_min_cells"] = mouse["n_cells"] >= 10
    scored = mouse.loc[mouse["pass_min_cells"]].copy()

    assoc_all = association_table(scored, "cldn4_mean", ["tnk_score", "ifn_mhc_score"])
    assoc_all.insert(0, "subset", "all_pass_min_cells")
    parts = [assoc_all]
    for geno, sub in scored.groupby("genotype"):
        a = association_table(sub, "cldn4_mean", ["tnk_score", "ifn_mhc_score"])
        a.insert(0, "subset", f"genotype={geno}")
        parts.append(a)
    assoc = pd.concat(parts, ignore_index=True)

    genes_used = pd.DataFrame(
        {
            "set": ["Cldn4", "T_NK", "IFN_MHC"],
            "requested": [",".join(CLDN4), ",".join(T_NK), ",".join(IFN_MHC)],
            "present": [
                ",".join(present["Cldn4"]),
                ",".join(present["T_NK"]),
                ",".join(present["IFN_MHC"]),
            ],
            "n_present": [len(present["Cldn4"]), len(present["T_NK"]), len(present["IFN_MHC"])],
        }
    )

    mouse.to_csv(OUT / "mouse_level_scores.tsv", sep="\t", index=False)
    assoc.to_csv(OUT / "high_vs_low_summary.tsv", sep="\t", index=False)
    genes_used.to_csv(OUT / "genes_used.tsv", sep="\t", index=False)
    cell.reset_index().to_csv(OUT / "cell_level_scores.tsv", sep="\t", index=False)

    n_all = int(mouse["bio_mouse"].nunique())
    n_scored = int(scored["bio_mouse"].nunique())
    n_tumors = int(cell["mouseID"].nunique())
    print(
        f"GSE154989: {len(cell)} cells, {n_tumors} tumor IDs, {n_all} bio mice, "
        f"{n_scored} mice with >=10 cells. Cldn4 present={bool(present['Cldn4'])}"
    )
    print(scored.groupby("genotype").size().to_string())
    print(assoc.to_string(index=False))


if __name__ == "__main__":
    main()
