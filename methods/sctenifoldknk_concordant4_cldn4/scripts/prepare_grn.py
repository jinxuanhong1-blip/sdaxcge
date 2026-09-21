#!/usr/bin/env python3
"""Choose the scTenifoldKnk gene universe on the QC malignant subsample.

Background genes are the most variable non-family genes with detection >= 5%.
IFN / MHC-I / TJ genes that pass the same detection filter are added.
CLDN4 is required. The KO gene itself is not treated as a family member.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contrastlib import DATASETS, ROOT

MIN_PCT = 0.05
N_BACKGROUND = 400
MIN_BACKGROUND = 250
MAX_GENES = 750
SETS = json.loads((ROOT / "data" / "gene_sets.json").read_text())
IFN = set(SETS["IFN"])
MHC = set(SETS["MHC_I_APM"])
TJ = set(SETS["TJ"])
FAMILY = IFN | MHC | TJ


def load_grn(folder: Path):
    mat = spio.mmread(folder / "grn_matrix.mtx").tocsr()
    genes = pd.read_csv(folder / "grn_genes.tsv", header=None)[0].astype(str).str.upper().to_numpy()
    cells = pd.read_csv(folder / "grn_cells.tsv", sep="\t")
    if mat.shape[0] != len(genes) or mat.shape[1] != len(cells):
        raise SystemExit(f"shape mismatch in {folder}: {mat.shape} genes {len(genes)} cells {len(cells)}")
    if not pd.Index(genes).is_unique:
        uniq, inv = np.unique(genes, return_inverse=True)
        coo = mat.tocoo()
        mat = sparse.coo_matrix((coo.data, (inv[coo.row], coo.col)), shape=(len(uniq), mat.shape[1])).tocsr()
        genes = uniq
    return mat, genes, cells


def select_genes(mat: sparse.csr_matrix, genes: np.ndarray) -> pd.DataFrame:
    n_cells = mat.shape[1]
    detected = np.asarray((mat > 0).sum(axis=1)).ravel()
    pct = detected / n_cells
    lib = np.asarray(mat.sum(axis=0)).ravel()
    lib = np.where(lib == 0, 1.0, lib)
    # variance of log1p(CPM) on genes that pass detection, to limit the dense block
    pass_det = pct >= MIN_PCT
    if int(pass_det.sum()) < 50:
        raise SystemExit("too few genes pass detection")
    sub = mat[pass_det]
    cpm = sub.multiply(1e6 / lib)
    logc = np.log1p(np.asarray(cpm.todense()))
    var = logc.var(axis=1)
    det_genes = genes[pass_det]
    info = pd.DataFrame(
        {
            "gene": det_genes,
            "pct": pct[pass_det],
            "variance": var,
            "detected_cells": detected[pass_det],
        }
    )
    info["in_IFN"] = info["gene"].isin(IFN)
    info["in_MHC"] = info["gene"].isin(MHC)
    info["in_TJ"] = info["gene"].isin(TJ)
    info["in_family"] = info["in_IFN"] | info["in_MHC"] | info["in_TJ"]
    if "CLDN4" not in set(info["gene"]):
        raise SystemExit("CLDN4 is below the 5% detection filter or absent")
    family = info.loc[info["in_family"]].copy()
    background = info.loc[~info["in_family"]].sort_values(["variance", "gene"], ascending=[False, True])
    n_bg = N_BACKGROUND
    chosen_bg = background.head(n_bg)
    chosen_fam = family
    total = len(chosen_bg) + len(chosen_fam) + (0 if "CLDN4" in set(chosen_fam["gene"]) or "CLDN4" in set(chosen_bg["gene"]) else 1)
    # CLDN4 is not in the families. It may already be in the background.
    while total > MAX_GENES and len(chosen_bg) > MIN_BACKGROUND:
        chosen_bg = chosen_bg.iloc[:-1]
        total = len(chosen_bg) + len(chosen_fam) + (0 if (chosen_bg["gene"] == "CLDN4").any() else 1)
    parts = [chosen_bg, chosen_fam]
    if not ((chosen_bg["gene"] == "CLDN4").any() or (chosen_fam["gene"] == "CLDN4").any()):
        parts.append(info.loc[info["gene"] == "CLDN4"])
    out = pd.concat(parts, ignore_index=True).drop_duplicates("gene")
    out["role"] = np.where(out["gene"] == "CLDN4", "ko_gene", np.where(out["in_family"], "family", "background"))
    return out.sort_values("gene")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, default=Path("/tmp/c4_work"))
    ap.add_argument("--dataset", required=True, choices=DATASETS)
    args = ap.parse_args()
    folder = args.work / args.dataset
    mat, genes, cells = load_grn(folder)
    print(f"{args.dataset} raw grn {mat.shape} nnz {mat.nnz}", flush=True)
    info = select_genes(mat, genes)
    keep = info["gene"].tolist()
    gene_to_i = {g: i for i, g in enumerate(genes)}
    ix = [gene_to_i[g] for g in keep]
    sub = mat[ix]
    # drop cells that became empty
    cell_sum = np.asarray(sub.sum(axis=0)).ravel()
    sub = sub[:, cell_sum > 0]
    cells = cells.loc[cell_sum > 0].reset_index(drop=True)
    out = folder / "knk_input"
    out.mkdir(exist_ok=True)
    spio.mmwrite(str(out / "counts.mtx"), sub, field="real")
    pd.Series(keep).to_csv(out / "genes.tsv", index=False, header=False)
    cells.to_csv(out / "cells.tsv", sep="\t", index=False)
    info.to_csv(folder / "gene_universe.tsv", sep="\t", index=False)
    print(
        f"{args.dataset} universe {len(keep)} "
        f"family {(info.role=='family').sum()} background {(info.role=='background').sum()} "
        f"cells {sub.shape[1]}",
        flush=True,
    )


if __name__ == "__main__":
    main()
