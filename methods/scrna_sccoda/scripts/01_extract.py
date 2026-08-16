#!/usr/bin/env python3
"""Extract a small gene panel from public GEO matrices. Never materializes the full matrix."""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.annotations import MARKER_PANEL  # noqa: E402


def extract_gse207422(matrix: Path, out: Path) -> None:
    wanted = set(MARKER_PANEL)
    t0 = time.time()
    rows: dict[str, np.ndarray] = {}
    ngenes = 0
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        ncells = len(cells)
        total = np.zeros(ncells, dtype=np.int64)
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != ncells:
                raise SystemExit(f"{gene}: {arr.size} != {ncells}")
            total += arr
            if gene in wanted:
                rows[gene] = arr
            if ngenes % 5000 == 0:
                print(f"  GSE207422 {ngenes} genes, kept {len(rows)}, {time.time()-t0:.0f}s", flush=True)
    genes = np.array([g for g in MARKER_PANEL if g in rows], dtype=object)
    mat = np.vstack([rows[g] for g in genes]).astype(np.int32)
    np.savez_compressed(out, cells=cells, total=total, genes=genes, mat=mat, n_genes=np.array([ngenes]))
    print(f"GSE207422: {ngenes} x {ncells}; kept {list(genes)}; {time.time()-t0:.0f}s -> {out}", flush=True)


def _feature_index(path: Path) -> dict[str, int]:
    idx = {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            gene = line.split("\t", 1)[0]
            idx[gene] = i
    return idx


def _barcodes(path: Path) -> np.ndarray:
    with gzip.open(path, "rt") as fh:
        return np.array([ln.rstrip("\n") for ln in fh], dtype=object)


def extract_mtx(mtx: Path, features: Path, barcodes: Path, out: Path) -> None:
    feat = _feature_index(features)
    wanted_row = {feat[g]: g for g in MARKER_PANEL if g in feat}
    print(f"  MTX wanted rows: {wanted_row}", flush=True)
    cells = _barcodes(barcodes)
    ncells = len(cells)
    store = {g: np.zeros(ncells, dtype=np.int32) for g in wanted_row.values()}
    t0 = time.time()
    with gzip.open(mtx, "rt") as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise SystemExit(f"not MTX: {mtx}")
        dims = fh.readline()
        while dims.startswith("%"):
            dims = fh.readline()
        nrows, ncols, nnz = map(int, dims.split())
        if ncols != ncells:
            raise SystemExit(f"ncols {ncols} != barcodes {ncells}")
        seen = 0
        for line in fh:
            seen += 1
            i_s, j_s, v_s = line.split()
            i = int(i_s)
            if i in wanted_row:
                store[wanted_row[i]][int(j_s) - 1] = int(v_s)
            if seen % 20_000_000 == 0:
                print(f"  {mtx.name} {seen}/{nnz} nnz {time.time()-t0:.0f}s", flush=True)
    genes = np.array([g for g in MARKER_PANEL if g in store], dtype=object)
    mat = np.vstack([store[g] for g in genes]).astype(np.int32)
    np.savez_compressed(
        out,
        cells=cells,
        genes=genes,
        mat=mat,
        n_genes=np.array([nrows]),
        n_cells=np.array([ncols]),
    )
    print(f"{mtx.name}: kept {list(genes)}; {time.time()-t0:.0f}s -> {out}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--which", choices=["gse207422", "gse241934_iit", "gse241934_rwc", "all"], default="all")
    args = ap.parse_args()
    d207 = args.root / "data" / "GSE207422"
    d241 = args.root / "data" / "GSE241934"
    if args.which in {"gse207422", "all"}:
        extract_gse207422(
            d207 / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
            d207 / "extracted_markers.npz",
        )
    if args.which in {"gse241934_iit", "all"}:
        extract_mtx(
            d241 / "GSE241934_IIT_Matrix.mtx.gz",
            d241 / "GSE241934_IIT_features.tsv.gz",
            d241 / "GSE241934_IIT_barcodes.tsv.gz",
            d241 / "extracted_IIT_markers.npz",
        )
    if args.which in {"gse241934_rwc", "all"}:
        extract_mtx(
            d241 / "GSE241934_Real_Matrix.mtx.gz",
            d241 / "GSE241934_RWC_features.tsv.gz",
            d241 / "GSE241934_RWC_barcodes.tsv.gz",
            d241 / "extracted_RWC_markers.npz",
        )


if __name__ == "__main__":
    main()
