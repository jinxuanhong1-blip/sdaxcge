#!/usr/bin/env python3
"""Two-pass extract of GSE207422 public UMI → epithelial h5ad.

Pass 1: lineage + A3-malignant-like from locked markers (same rule as the
given A3 / dual-high slices). Pass 2: keep the full gene axis for epithelial
barcodes only. Gzip is ~175 MB; uncompressed text is never materialized.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import A3_NORMAL_LUNG, LINEAGES, MARKER_GENES, PAPER_GROUP  # noqa: E402


def _score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([_score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def pass1_markers(matrix: Path, wanted: set[str]) -> dict:
    t0 = time.time()
    rows: dict[str, np.ndarray] = {}
    ngenes = 0
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        ncells = len(cells)
        total = np.zeros(ncells, dtype=np.int64)
        n_nonzero_genes = np.zeros(ncells, dtype=np.int32)
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != ncells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {ncells}")
            total += arr
            n_nonzero_genes += arr > 0
            if gene in wanted:
                rows[gene] = arr
            if ngenes % 4000 == 0:
                print(
                    f"  pass1 {ngenes} genes, {len(rows)}/{len(wanted)} kept, "
                    f"{time.time() - t0:.0f}s",
                    flush=True,
                )
    missing = sorted(wanted - set(rows))
    print(
        f"pass1 done: {ngenes} genes × {ncells} cells in {time.time() - t0:.0f}s; "
        f"missing {missing}",
        flush=True,
    )
    return {
        "cells": cells,
        "total": total,
        "n_nonzero_genes": n_nonzero_genes,
        "expr": rows,
        "n_genes_in_matrix": ngenes,
        "missing_markers": missing,
    }


def pass2_epithelial(matrix: Path, epi_idx: np.ndarray) -> tuple[csr_matrix, np.ndarray]:
    t0 = time.time()
    data_chunks: list[np.ndarray] = []
    indices_chunks: list[np.ndarray] = []
    indptr = [0]
    nnz = 0
    var_names: list[str] = []
    ngenes = 0
    n_epi = int(epi_idx.size)
    with gzip.open(matrix, "rt") as fh:
        fh.readline()
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)[epi_idx]
            nz = np.flatnonzero(arr)
            if nz.size:
                data_chunks.append(arr[nz].astype(np.int32, copy=False))
                indices_chunks.append(nz.astype(np.int32, copy=False))
                nnz += int(nz.size)
            indptr.append(nnz)
            var_names.append(gene)
            if ngenes % 4000 == 0:
                print(
                    f"  pass2 {ngenes} genes, nnz={nnz}, {time.time() - t0:.0f}s",
                    flush=True,
                )
    data = (
        np.concatenate(data_chunks)
        if data_chunks
        else np.array([], dtype=np.int32)
    )
    indices = (
        np.concatenate(indices_chunks)
        if indices_chunks
        else np.array([], dtype=np.int32)
    )
    X_genes = csr_matrix(
        (data, indices, np.asarray(indptr, dtype=np.int64)),
        shape=(ngenes, n_epi),
        dtype=np.int32,
    )
    X = X_genes.T.tocsr()
    print(
        f"pass2 done: {X.shape[0]} cells × {X.shape[1]} genes, nnz={X.nnz}, "
        f"{time.time() - t0:.0f}s",
        flush=True,
    )
    return X, np.array(var_names, dtype=object)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp/gse207422_traj_data"))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_xlsx = args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    out = args.out or (args.workdir / "epithelium.h5ad")

    p1 = pass1_markers(matrix, set(MARKER_GENES))
    cells = p1["cells"]
    n = len(cells)
    expr = p1["expr"]
    lineage = assign_lineage(expr, n)
    is_epi = lineage == "epithelial"
    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_malig = is_epi & (a3_normal == 0)
    is_leftover = is_epi & (a3_normal > 0)

    epi_idx = np.flatnonzero(is_epi)
    print(
        f"lineage epithelial={int(is_epi.sum())} A3-malignant={int(is_malig.sum())} "
        f"leftover={int(is_leftover.sum())} / {n}",
        flush=True,
    )
    if epi_idx.size < 50:
        raise SystemExit(f"too few epithelial cells: {epi_idx.size}")

    X, var_names = pass2_epithelial(matrix, epi_idx)

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells[epi_idx]])
    obs = pd.DataFrame(
        {
            "cell": cells[epi_idx],
            "Sample": sample_ids,
            "patient": np.array(
                [s.replace("BD_immune", "P") for s in sample_ids], dtype=object
            ),
            "paper_group": pd.Series(sample_ids).map(PAPER_GROUP).to_numpy(),
            "lineage": lineage[epi_idx],
            "is_malig_a3": is_malig[epi_idx],
            "is_leftover_epi": is_leftover[epi_idx],
            "n_umi": p1["total"][epi_idx],
            "n_genes": p1["n_nonzero_genes"][epi_idx],
        }
    )
    obs["compartment"] = np.where(obs["is_malig_a3"], "A3_malignant", "leftover_epi")
    obs["timing"] = np.where(obs["paper_group"].astype(str).eq("TN"), "pre", "post")
    obs.index = obs["cell"].astype(str)

    if meta_xlsx.exists():
        meta = load_sample_meta(meta_xlsx)
        keep = [
            c
            for c in [
                "Sample",
                "Patient",
                "Resource",
                "Sex",
                "Age",
                "Clinical Stage",
                "Pathology",
                "PD1 Antibody",
                "Chemotherapy",
                "Pathologic Response",
                "Residual Tumor",
                "RECIST",
                "path_response",
            ]
            if c in meta.columns
        ]
        obs = obs.merge(meta[keep], on="Sample", how="left")
        for col in list(obs.columns):
            if col in {"n_umi", "n_genes"} or pd.api.types.is_bool_dtype(obs[col]):
                continue
            if pd.api.types.is_numeric_dtype(obs[col]) and col not in {
                "Residual Tumor",
                "Age",
            }:
                continue
            obs[col] = obs[col].map(lambda x: "" if pd.isna(x) else str(x))
        obs.index = obs["cell"].astype(str)

    import anndata as ad

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(var_names, name="gene")),
    )
    adata.uns["extract"] = {
        "accession": "GSE207422",
        "n_cells_matrix": int(n),
        "n_genes_matrix": int(p1["n_genes_in_matrix"]),
        "n_epithelial": int(is_epi.sum()),
        "n_malig_a3": int(is_malig.sum()),
        "n_leftover_epi": int(is_leftover.sum()),
        "missing_markers": p1["missing_markers"],
        "a3_normal_lung": list(A3_NORMAL_LUNG),
        "author_cell_labels_public": False,
        "public_umi_only": True,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    (out.with_suffix(".extract.json")).write_text(
        json.dumps(adata.uns["extract"], indent=2)
    )
    print(f"wrote {out}  {adata.n_obs} × {adata.n_vars}", flush=True)


if __name__ == "__main__":
    main()
