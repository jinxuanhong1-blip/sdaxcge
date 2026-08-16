#!/usr/bin/env python3
"""Download open-access inputs for the B1_THCA TACSTD2–CLDN4 analysis.

All files are fully open (no dbGaP / controlled access):

  * UCSC Xena GDC hub STAR TPM for TCGA-THCA. Values are log2(TPM+1),
    GENCODE v36 Ensembl gene IDs. (A previous draft treated this matrix as
    log2(TPM+0.001); the observed floor is 0.0, which matches log2(TPM+1).)
  * GENCODE v36 probemap (Ensembl id → HGNC symbol).
  * GDC clinical phenotype (histology / morphology).
  * Bausch-Fluck 2018 SURFY in-silico surfaceome workbook. The canonical
    host (wlab.ethz.ch/surfaceome) now serves HTML for the file path, so
    we pull the identical workbook vendored in steveneschrich/surfaceome.
  * ABSOLUTE tumor purity from the open PanCanAtlas supplement (GDC API).

Raw matrices are intentionally NOT committed. Default destination is
data/raw/ (override with B1_THCA_DATA).
"""
from __future__ import annotations

import hashlib
import os
import sys
import urllib.request

DATA_DIR = os.environ.get(
    "B1_THCA_DATA",
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "raw"),
)
DATA_DIR = os.path.abspath(DATA_DIR)

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"

FILES = {
    "TCGA-THCA.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-THCA.star_tpm.tsv.gz",
    "TCGA-THCA.clinical.tsv.gz": f"{GDC_HUB}/TCGA-THCA.clinical.tsv.gz",
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    "table_S3_surfaceome.xlsx": (
        "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
        "main/data-raw/surfy/table_S3_surfaceome.xlsx"
    ),
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
}


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str, url: str) -> None:
    dest = os.path.join(DATA_DIR, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"[skip] {name} already present  size={os.path.getsize(dest):,}  md5={md5(dest)}")
        return
    print(f"[get ] {url}")
    os.makedirs(DATA_DIR, exist_ok=True)
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    print(f"[done] {name}  size={os.path.getsize(dest):,}  md5={md5(dest)}")


def main() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)
    for name, url in FILES.items():
        fetch(name, url)
    print(f"\nAll files in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
