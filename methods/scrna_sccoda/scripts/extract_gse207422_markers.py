#!/usr/bin/env python3
"""Keep lineage-marker + TACSTD2 rows from the GSE207422 dense UMI TSV."""
from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
OUT = ROOT / "data" / "GSE207422" / "GSE207422_markers.tsv.gz"

KEEP = {
    "TACSTD2", "EPCAM", "KRT19", "KRT18", "KRT8", "KRT7", "KRT5", "KRT17",
    "NAPSA", "NKX2-1", "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD4",
    "IL7R", "NKG7", "GNLY", "KLRD1", "NCAM1", "GZMB", "PRF1", "MS4A1", "CD79A",
    "LYZ", "CD68", "FCGR3A", "CD14", "PECAM1", "VWF", "COL1A1", "DCN",
    "PDGFRA", "TPSAB1", "KIT",
}


def main() -> None:
    n = 0
    with gzip.open(SRC, "rt") as fin, gzip.open(OUT, "wt") as fout:
        header = fin.readline()
        fout.write(header)
        n += 1
        for line in fin:
            gene = line.split("\t", 1)[0]
            if gene in KEEP:
                fout.write(line)
                n += 1
    print("wrote", OUT, "rows_including_header", n)


if __name__ == "__main__":
    main()
