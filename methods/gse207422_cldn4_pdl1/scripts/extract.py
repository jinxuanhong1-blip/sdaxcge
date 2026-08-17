#!/usr/bin/env python3
"""Stream GSE207422 genes×cells UMI; keep CLDN4, CD274, HLA, lineage genes.

Gzip is ~175 MB. Uncompressed text is never materialized.
"""

from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path

import numpy as np

WANTED = [
    "CLDN4",
    "CD274",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "HLA-DRA",
    "HLA-DRB1",
    "CD74",
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "SFTPA2",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "TPPP3",
    "CD3D",
    "CD3E",
    "CD2",
    "NKG7",
    "GNLY",
    "FGFBP2",
    "CD79A",
    "MS4A1",
    "IGHG1",
    "MZB1",
    "LYZ",
    "CD68",
    "CD14",
    "CSF3R",
    "COL1A1",
    "DCN",
    "VWF",
    "PECAM1",
    "KIT",
    "PTPRC",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    out = args.workdir / "extracted_cldn4_pdl1_markers.npz"
    wanted = set(WANTED)

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
                    f"  {ngenes} genes, {len(rows)}/{len(wanted)} kept, {time.time() - t0:.0f}s",
                    flush=True,
                )

    missing = sorted(wanted - set(rows))
    genes = np.array([g for g in WANTED if g in rows], dtype=object)
    mat = np.vstack([rows[g] for g in genes]).astype(np.int32)
    np.savez_compressed(
        out,
        cells=cells,
        total=total,
        n_nonzero_genes=n_nonzero_genes,
        genes=genes,
        mat=mat,
        n_genes_in_matrix=np.array([ngenes]),
    )
    print(
        f"parsed {ngenes} genes × {ncells} cells in {time.time() - t0:.0f}s; "
        f"kept {len(genes)}; missing {missing}; saved {out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
