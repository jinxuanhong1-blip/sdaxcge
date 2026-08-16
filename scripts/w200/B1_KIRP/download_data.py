#!/usr/bin/env python3
"""Download the two public inputs for the B1 analog in TCGA-KIRP.

This slice is the KIRP counterpart of the BRCA B1 analog
(``results/w200/B1_BRCA``): rank every surfaceome gene by co-expression
with TACSTD2 (TROP2) and report where CLDN4 actually lands.

Inputs (written to ``data/``, git-ignored):

1. UCSC Xena TCGA hub ``TCGA.KIRP.sampleMap/HiSeqV2``.
   Values are log2(norm_count + 1). Rows are HGNC symbols, columns are
   TCGA barcodes. Same matrix family as the BRCA B1 analog, so the KIRP
   rank is comparable without an Ensembl-to-symbol remap.

2. Bausch-Fluck et al. 2018 in-silico human surfaceome, table S3.
   The ETH / Wollscheid host now serves a Git-LFS pointer for the xlsx,
   so we use the identical workbook vendored in
   steveneschrich/surfaceome (GitHub).

Re-running is cheap: existing files are left in place unless --force.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[3] / "data"

EXPRESSION_URL = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.KIRP.sampleMap%2FHiSeqV2.gz"
)
EXPRESSION_FILE = DATA_DIR / "TCGA-KIRP.HiSeqV2.gz"

SURFACEOME_URL = (
    "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
    "main/data-raw/surfy/table_S3_surfaceome.xlsx"
)
SURFACEOME_FILE = DATA_DIR / "table_S3_surfaceome.xlsx"

UA = "w200-B1-KIRP/1.0 (public TCGA + SURFY download)"


def _download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
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
