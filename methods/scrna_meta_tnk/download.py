#!/usr/bin/env python3
"""Download processed GEO files used by this meta (public FTP only; skip EGA/dbGaP)."""
from __future__ import annotations

import urllib.request
from pathlib import Path

DATA = Path("/tmp/scrna_meta_tnk")
FILES = {
    "GSE241934": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
    ],
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print("have", dest)
        return
    print("get", url)
    urllib.request.urlretrieve(url, dest)


def main() -> None:
    for cohort, urls in FILES.items():
        for url in urls:
            fetch(url, DATA / cohort / Path(url).name)
    print("GSE207422 is taken as given (no download).")
    print("Other must-try series use deposited processed-GEO patient tables under data/existing/.")


if __name__ == "__main__":
    main()
