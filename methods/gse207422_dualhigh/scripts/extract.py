#!/usr/bin/env python3
"""Stream GSE207422 genes×cells UMI; keep CLDN4, companion TACSTD2, lineage, T-state genes.

Gzip is ~175 MB. Uncompressed text is never materialized.
"""

from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path

import numpy as np

WANTED = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "CDH1",
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
    "CAPS",
    "CD3D",
    "CD3E",
    "CD3G",
    "CD2",
    "TRAC",
    "CD4",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "KLRF1",
    "NCR1",
    "FGFBP2",
    "NCAM1",
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
    "CXCL13",
    "GZMB",
    "GZMA",
    "GZMK",
    "PRF1",
    "IFNG",
    "MKI67",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    out = args.workdir / "extracted_cldn4_markers.npz"
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
