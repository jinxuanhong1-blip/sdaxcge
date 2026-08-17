#!/usr/bin/env python3
"""Stream GSE207422 genes×cells UMI and keep the CLDN4 / T/NK score panel.

The gzip is ~176 MB; uncompressed text is never written to disk.
"""
from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path

import numpy as np

PANEL = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT5",
    "KRT7",
    "KRT8",
    "KRT17",
    "KRT18",
    "KRT19",
    "ELF3",
    "CDH1",
    "MUC1",
    "CEACAM5",
    "CEACAM6",
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPD",
    "AGER",
    "NAPSA",
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
    "CAPS",
    "CD3D",
    "CD3E",
    "CD3G",
    "TRAC",
    "CD2",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "KLRF1",
    "FGFBP2",
    "GZMB",
    "PRF1",
    "GZMA",
    "IFNG",
    "PDCD1",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "TOX",
    "CTLA4",
    "PTPRC",
    "CD79A",
    "CD79B",
    "MS4A1",
    "JCHAIN",
    "MZB1",
    "LYZ",
    "CD68",
    "CD14",
    "C1QA",
    "C1QB",
    "FCN1",
    "TPSAB1",
    "TPSB2",
    "CPA3",
    "MS4A2",
    "PECAM1",
    "VWF",
    "CLDN5",
    "CDH5",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "TAGLN",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse207422_cldn4_exh"))
    args = ap.parse_args()
    matrix = args.datadir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    out = args.datadir / "extracted_panel.npz"
    wanted = set(PANEL)
    t0 = time.time()
    rows: dict[str, np.ndarray] = {}
    ngenes = 0
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        ncells = len(cells)
        total = np.zeros(ncells, dtype=np.int64)
        print(f"cells {ncells}", flush=True)
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != ncells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {ncells}")
            total += arr
            if gene in wanted:
                rows[gene] = arr
            if ngenes % 4000 == 0:
                print(
                    f"  {ngenes} genes, {len(rows)}/{len(wanted)} kept, "
                    f"{time.time() - t0:.0f}s",
                    flush=True,
                )
    missing = sorted(wanted - set(rows))
    genes = np.array([g for g in PANEL if g in rows], dtype=object)
    mat = np.vstack([rows[g] for g in genes]).astype(np.int32)
    np.savez_compressed(
        out,
        cells=cells,
        total=total,
        genes=genes,
        mat=mat,
        n_genes_in_matrix=np.array([ngenes]),
        missing=np.array(missing, dtype=object),
    )
    print(
        f"parsed {ngenes} genes × {ncells} cells in {time.time() - t0:.0f}s; "
        f"kept {len(genes)}; missing {missing}; saved {out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
