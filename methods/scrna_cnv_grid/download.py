#!/usr/bin/env python3
"""Download public GSE207422 and GSE241934 processed files (GEO only)."""
from __future__ import annotations

import urllib.request
from pathlib import Path

GSE207422 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
GSE241934 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl"

FILES = {
    "data/GSE207422": [
        (f"{GSE207422}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"),
        (f"{GSE207422}/GSE207422_NSCLC_scRNAseq_metadata.xlsx", "GSE207422_NSCLC_scRNAseq_metadata.xlsx"),
    ],
    "data/GSE241934": [
        (f"{GSE241934}/GSE241934_IIT_Matrix.mtx.gz", "GSE241934_IIT_Matrix.mtx.gz"),
        (f"{GSE241934}/GSE241934_IIT_Meta.txt.gz", "GSE241934_IIT_Meta.txt.gz"),
        (f"{GSE241934}/GSE241934_IIT_barcodes.tsv.gz", "GSE241934_IIT_barcodes.tsv.gz"),
        (f"{GSE241934}/GSE241934_IIT_features.tsv.gz", "GSE241934_IIT_features.tsv.gz"),
        (f"{GSE241934}/GSE241934_Real_Matrix.mtx.gz", "GSE241934_Real_Matrix.mtx.gz"),
        (f"{GSE241934}/GSE241934_Real_Meta.txt.gz", "GSE241934_Real_Meta.txt.gz"),
        (f"{GSE241934}/GSE241934_RWC_barcodes.tsv.gz", "GSE241934_RWC_barcodes.tsv.gz"),
        (f"{GSE241934}/GSE241934_RWC_features.tsv.gz", "GSE241934_RWC_features.tsv.gz"),
    ],
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"downloading {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    print(f"  -> {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    for folder, items in FILES.items():
        for url, name in items:
            fetch(url, Path(folder) / name)


if __name__ == "__main__":
    main()
