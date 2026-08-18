#!/usr/bin/env python3
"""Stream public GSE207422 genes×cells UMI TSV into 10x MTX for Seurat ReadMtx.

Does not load a dense genes×cells array. One gzip pass; sparse COO in memory.
"""
from __future__ import annotations

import gzip
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

REPO = Path(__file__).resolve().parents[3]
DEFAULT_SRC = REPO / "data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
DEFAULT_OUT = REPO / "data/GSE207422/tenx"


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    SRC = args.src
    OUT = args.out
    if not SRC.exists() or SRC.stat().st_size == 0:
        raise SystemExit(f"UMI matrix missing: {SRC} — stop.")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data_chunks: list[np.ndarray] = []
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    genes: list[str] = []

    with gzip.open(SRC, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        n_cells = len(barcodes)
        print(f"cells={n_cells} header0={header[0]!r}", flush=True)
        for gi, line in enumerate(fh):
            tab = line.find("\t")
            gene = line[:tab]
            arr = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.int32)
            if arr.size != n_cells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
            nz = np.flatnonzero(arr)
            if nz.size:
                data_chunks.append(arr[nz].astype(np.int32, copy=False))
                row_chunks.append(np.full(nz.size, gi, dtype=np.int32))
                col_chunks.append(nz.astype(np.int32, copy=False))
            genes.append(gene)
            if (gi + 1) % 4000 == 0:
                nnz = sum(int(x.size) for x in data_chunks)
                print(f"  {gi+1} genes nnz={nnz} {time.time()-t0:.0f}s", flush=True)

    n_genes = len(genes)
    if not data_chunks:
        raise SystemExit("no nonzeros")
    data = np.concatenate(data_chunks)
    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    del data_chunks, row_chunks, col_chunks
    mat = sparse.csc_matrix((data, (rows, cols)), shape=(n_genes, n_cells), dtype=np.int32)
    del data, rows, cols
    print(
        f"sparse {n_genes} x {n_cells} nnz={mat.nnz} "
        f"after {time.time()-t0:.0f}s",
        flush=True,
    )

    feat = OUT / "features.tsv"
    bc = OUT / "barcodes.tsv"
    mtx = OUT / "matrix.mtx"
    with feat.open("w") as fh:
        for g in genes:
            fh.write(f"{g}\t{g}\tGene Expression\n")
    with bc.open("w") as fh:
        for b in barcodes:
            fh.write(f"{b}\n")
    mmwrite(str(mtx), mat)
    # gzip mtx (Seurat ReadMtx accepts .mtx.gz)
    import subprocess

    subprocess.check_call(["gzip", "-f", str(mtx)])
    subprocess.check_call(["gzip", "-f", str(feat)])
    subprocess.check_call(["gzip", "-f", str(bc)])
    print(f"wrote {OUT} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
