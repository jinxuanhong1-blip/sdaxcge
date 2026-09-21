#!/usr/bin/env python3
"""Stream GSE207422 genes x cells UMI. Keep the NHEJ/STING/IFN panel.

Gzip is ~175 MB. The uncompressed text is never written to disk.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import requested_symbols


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    out = args.workdir / "extracted_nhej_sting_ifn.npz"
    wanted = set(requested_symbols())
    if not matrix.exists():
        raise SystemExit(f"missing {matrix}; run download.py")

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
            total += arr.astype(np.int64)
            n_nonzero_genes += arr > 0
            if gene in wanted:
                rows[gene] = arr
            if ngenes % 4000 == 0:
                print(
                    f"  {ngenes} genes, {len(rows)}/{len(wanted)} kept, {time.time() - t0:.0f}s",
                    flush=True,
                )

    order = [g for g in requested_symbols() if g in rows]
    missing = [g for g in requested_symbols() if g not in rows]
    genes = np.array(order, dtype=object)
    mat = np.vstack([rows[g] for g in order]).astype(np.int32) if order else np.zeros((0, ncells), dtype=np.int32)
    np.savez_compressed(
        out,
        cells=cells,
        total=total,
        n_nonzero_genes=n_nonzero_genes,
        genes=genes,
        mat=mat,
        n_genes_in_matrix=np.array([ngenes]),
        missing=np.array(missing, dtype=object),
    )
    gate = {g: g in rows for g in ("TACSTD2", "CLDN4")}
    print(
        f"parsed {ngenes} genes x {ncells} cells in {time.time() - t0:.0f}s; "
        f"kept {len(order)}; missing {missing}",
        flush=True,
    )
    print(f"presence gate TACSTD2={gate['TACSTD2']} CLDN4={gate['CLDN4']}", flush=True)
    print(f"saved {out}", flush=True)


if __name__ == "__main__":
    main()
