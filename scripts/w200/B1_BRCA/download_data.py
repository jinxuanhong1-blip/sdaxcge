"""Download public inputs for the TCGA-BRCA B1 analog.

Four files land in ``data/`` (or ``$B1_BRCA_DATA``). The ~64 MB expression
matrix is git-ignored; the other three are small enough to re-fetch cheaply.

1. UCSC Xena TCGA hub ``TCGA.BRCA.sampleMap/HiSeqV2`` — RNA-seq,
   log2(norm_count + 1), HGNC-symbol rows. Same matrix as the first-pass
   ranking, so the headline ranks stay comparable.
2. Xena ``BRCA_clinicalMatrix`` — PAM50Call_RNAseq, histology, ER/PR/HER2.
3. PanCanAtlas ABSOLUTE purity (GDC open file
   ``4f277128-f793-4354-a13d-30cc7fe9f6b5``).
4. Bausch-Fluck 2018 in-silico surfaceome table S3, pulled from the
   steveneschrich/surfaceome GitHub mirror (the ETH host now serves an SPA
   HTML page at the original file path).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("B1_BRCA_DATA", ROOT / "data"))

FILES = {
    "TCGA-BRCA.HiSeqV2.gz": (
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz"
    ),
    "BRCA_clinicalMatrix": (
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"
    ),
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
    "table_S3_surfaceome.xlsx": (
        "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
        "main/data-raw/surfy/table_S3_surfaceome.xlsx"
    ),
}


def _download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        tmp.replace(dest)
    print(f"[ok  ] {dest.name} ({dest.stat().st_size:,} bytes)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        _download(url, DATA_DIR / name, args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
