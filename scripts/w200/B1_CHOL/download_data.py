#!/usr/bin/env python3
"""Download open-access inputs for the B1_CHOL TACSTD2–CLDN4 analysis.

All files are fully open (no dbGaP / controlled access):

  1. UCSC Xena TCGA hub HiSeqV2 expression for TCGA-CHOL
     (log2(norm_count+1), HGNC symbols). Same matrix family as B1_BRCA.
  2. Bausch-Fluck et al. 2018 in-silico surfaceome table S3, vendored copy
     from steveneschrich/surfaceome (the ETH host now serves HTML/LFS).
  3. ABSOLUTE tumor purity/ploidy, PanCanAtlas open supplement (GDC API).
  4. Xena GDC hub clinical phenotype for TCGA-CHOL (histology check).

Raw matrices land in DATA_DIR (default /tmp/b1_chol_data; override with
B1_CHOL_DATA) and are intentionally not committed.

Usage:  python3 scripts/w200/B1_CHOL/download_data.py
"""
from __future__ import annotations

import hashlib
import os
import sys
import urllib.request

DATA_DIR = os.environ.get("B1_CHOL_DATA", "/tmp/b1_chol_data")

FILES = {
    "CHOL.HiSeqV2.gz": (
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.CHOL.sampleMap%2FHiSeqV2.gz"
    ),
    "table_S3_surfaceome.xlsx": (
        "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
        "main/data-raw/surfy/table_S3_surfaceome.xlsx"
    ),
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
    "TCGA-CHOL.clinical.tsv.gz": (
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
        "TCGA-CHOL.clinical.tsv.gz"
    ),
}


def fetch(name: str, url: str) -> None:
    dest = os.path.join(DATA_DIR, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"[skip] {name} already present ({os.path.getsize(dest):,} bytes)")
        return
    print(f"[get ] {url}")
    os.makedirs(DATA_DIR, exist_ok=True)
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    md5 = hashlib.md5(open(dest, "rb").read()).hexdigest()
    print(f"[done] {name}  size={os.path.getsize(dest):,}  md5={md5}")


def main() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)
    for name, url in FILES.items():
        fetch(name, url)
    print(f"\nAll files in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
