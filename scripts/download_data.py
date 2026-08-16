"""Download the raw inputs needed for the TACSTD2 surfaceome co-expression analysis.

Two inputs are fetched into ``data/`` (git-ignored, because the expression
matrix is ~64 MB):

1. TCGA-BRCA RNA-seq gene-expression matrix from the UCSC Xena TCGA hub
   (dataset ``TCGA.BRCA.sampleMap/HiSeqV2``). Values are log2(norm_count + 1),
   rows are HGNC gene symbols, columns are TCGA sample barcodes. Using the
   symbol-indexed matrix avoids an Ensembl->symbol mapping step and makes the
   focus genes (TACSTD2, CLDN4) directly addressable.

2. The "in silico human surfaceome" master table (Bausch-Fluck et al.,
   PNAS 2018). The canonical host (wlab.ethz.ch) now serves a single-page app
   that returns HTML for the file path, so we pull the identical workbook that
   is vendored in the steveneschrich/surfaceome R package on GitHub.

Re-running is cheap: existing files are left in place unless --force is given.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

EXPRESSION_URL = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz"
)
EXPRESSION_FILE = DATA_DIR / "TCGA-BRCA.HiSeqV2.gz"

SURFACEOME_URL = (
    "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
    "main/data-raw/surfy/table_S3_surfaceome.xlsx"
)
SURFACEOME_FILE = DATA_DIR / "table_S3_surfaceome.xlsx"


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
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    args = ap.parse_args(argv)
    _download(EXPRESSION_URL, EXPRESSION_FILE, args.force)
    _download(SURFACEOME_URL, SURFACEOME_FILE, args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
