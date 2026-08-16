#!/usr/bin/env python3
"""Stream the GSE207422 genes×cells UMI TSV into a compact HVG+marker matrix.

The gzip is ~175 MB; uncompressed text is a few GB and is never materialized.
Pass 1: library sizes + per-gene mean/variance + marker rows.
Pass 2: write CSR counts for HVGs ∪ markers.
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse

MARKERS = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT17",
    "CDH1",
    "DST",
    "SFTPA2",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "TPPP3",
    "CD3D",
    "CD3E",
    "CD2",
    "TRAC",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "FGFBP2",
    "PTPRC",
    "LYZ",
    "CD68",
    "CD14",
    "CSF3R",
    "CD79A",
    "MS4A1",
    "MZB1",
    "IGHG1",
    "COL1A1",
    "DCN",
    "PECAM1",
    "VWF",
    "KIT",
    "LILRA4",
    "MKI67",
    "CX3CL1",
    "CD74",
    "HLA-DRA",
]


def _open_matrix(path: Path):
    fh = gzip.open(path, "rt")
    header = fh.readline().rstrip("\n").split("\t")
    if header[0] != "Gene":
        raise SystemExit(f"unexpected header[0]={header[0]!r}")
    cells = np.array(header[1:], dtype=object)
    return fh, cells


def pass1(matrix: Path) -> dict:
    wanted = set(MARKERS)
    t0 = time.time()
    fh, cells = _open_matrix(matrix)
    n_cells = len(cells)
    total = np.zeros(n_cells, dtype=np.int64)
    n_genes = 0
    gene_names: list[str] = []
    gene_sum: list[float] = []
    gene_sumsq: list[float] = []
    gene_nnz: list[int] = []
    marker_rows: dict[str, np.ndarray] = {}
    try:
        for line in fh:
            n_genes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != n_cells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
            total += arr
            s = float(arr.sum())
            gene_names.append(gene)
            gene_sum.append(s)
            gene_sumsq.append(float(np.dot(arr.astype(np.float64), arr.astype(np.float64))))
            gene_nnz.append(int((arr > 0).sum()))
            if gene in wanted:
                marker_rows[gene] = arr
            if n_genes % 4000 == 0:
                print(
                    f"  pass1 {n_genes} genes, {len(marker_rows)} markers, {time.time() - t0:.0f}s",
                    flush=True,
                )
    finally:
        fh.close()
    print(f"pass1 done: {n_genes} genes × {n_cells} cells in {time.time() - t0:.0f}s", flush=True)
    return {
        "cells": cells,
        "total": total,
        "gene_names": np.array(gene_names, dtype=object),
        "gene_sum": np.array(gene_sum, dtype=np.float64),
        "gene_sumsq": np.array(gene_sumsq, dtype=np.float64),
        "gene_nnz": np.array(gene_nnz, dtype=np.int32),
        "marker_rows": marker_rows,
        "n_genes": n_genes,
        "n_cells": n_cells,
    }


def select_hvgs(stats: dict, n_hvg: int) -> list[str]:
    n = stats["n_cells"]
    mean = stats["gene_sum"] / n
    var = stats["gene_sumsq"] / n - mean**2
    var = np.maximum(var, 0)
    disp = var / np.maximum(mean, 1e-12)
    ok = (stats["gene_nnz"] >= 20) & (mean >= 0.01)
    names = stats["gene_names"]
    order = np.argsort(-disp)
    hvgs: list[str] = []
    for i in order:
        if not ok[i]:
            continue
        hvgs.append(str(names[i]))
        if len(hvgs) >= n_hvg:
            break
    keep = []
    seen = set()
    for g in MARKERS + hvgs:
        if g in seen:
            continue
        if g in set(names.tolist()) or g in stats["marker_rows"]:
            keep.append(g)
            seen.add(g)
    return keep


def pass2(matrix: Path, keep_genes: list[str], n_cells: int) -> tuple[list[str], sparse.csr_matrix]:
    wanted = set(keep_genes)
    t0 = time.time()
    fh, _cells = _open_matrix(matrix)
    rows = []
    genes = []
    n = 0
    try:
        for line in fh:
            n += 1
            i = line.find("\t")
            gene = line[:i]
            if gene not in wanted:
                continue
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            rows.append(sparse.csr_matrix(arr.reshape(1, -1)))
            genes.append(gene)
            if len(genes) % 200 == 0:
                print(f"  pass2 kept {len(genes)}/{len(keep_genes)}, {time.time() - t0:.0f}s", flush=True)
    finally:
        fh.close()
    if not rows:
        raise SystemExit("pass2 kept 0 genes")
    mat = sparse.vstack(rows, format="csr")  # genes × cells
    print(f"pass2 done: {mat.shape} in {time.time() - t0:.0f}s", flush=True)
    # preserve requested order
    pos = {g: i for i, g in enumerate(genes)}
    ordered = [g for g in keep_genes if g in pos]
    idx = [pos[g] for g in ordered]
    return ordered, mat[idx]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--n-hvg", type=int, default=2500)
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    stats = pass1(matrix)
    keep = select_hvgs(stats, args.n_hvg)
    missing_markers = [g for g in MARKERS if g not in set(stats["gene_names"])]
    genes, gene_x_cell = pass2(matrix, keep, stats["n_cells"])
    cell_x_gene = gene_x_cell.T.tocsr()
    out = args.workdir / "hvg_marker_counts.npz"
    sparse.save_npz(args.workdir / "counts_csr.npz", cell_x_gene)
    np.savez_compressed(
        out,
        cells=stats["cells"],
        genes=np.array(genes, dtype=object),
        total=stats["total"],
        n_genes_in_matrix=np.array([stats["n_genes"]]),
        missing_markers=np.array(missing_markers, dtype=object),
    )
    # also stash marker-only dense for lineage (small)
    marker_genes = [g for g in MARKERS if g in stats["marker_rows"]]
    marker_mat = np.vstack([stats["marker_rows"][g] for g in marker_genes]).astype(np.int32)
    np.savez_compressed(
        args.workdir / "markers.npz",
        cells=stats["cells"],
        genes=np.array(marker_genes, dtype=object),
        mat=marker_mat,
        total=stats["total"],
    )
    meta = {
        "n_cells": int(stats["n_cells"]),
        "n_genes_in_matrix": int(stats["n_genes"]),
        "n_genes_kept": len(genes),
        "n_hvg_requested": args.n_hvg,
        "missing_markers": missing_markers,
        "counts_shape": list(cell_x_gene.shape),
        "counts_nnz": int(cell_x_gene.nnz),
    }
    (args.workdir / "stream_summary.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
