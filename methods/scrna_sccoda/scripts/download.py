#!/usr/bin/env python3
"""Download public GEO supplementary files used by this folder. No FASTQ."""
from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

FILES = {
    "GSE207422": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    ],
    "GSE241934": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
    ],
    "GSE291670": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291670/suppl/GSE291670_RAW.tar",
    ],
}


def main() -> None:
    for cohort, urls in FILES.items():
        dest = DATA / cohort
        dest.mkdir(parents=True, exist_ok=True)
        for url in urls:
            name = url.rsplit("/", 1)[-1]
            path = dest / name
            if path.exists() and path.stat().st_size > 0:
                print("have", path)
                continue
            print("get", url, flush=True)
            urllib.request.urlretrieve(url, path)
            print(" wrote", path, path.stat().st_size)


if __name__ == "__main__":
    main()
