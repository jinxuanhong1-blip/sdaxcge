#!/usr/bin/env python3
"""Download GSE207422 processed GEO supplementary files (author UMI matrix + sample metadata).

No barcode-level author cell-type file exists on GEO or in the Genome Medicine
supplements (checked Additional files 1–4). Cell labels are reconstructed in
analyze.py from the paper's published canonical-marker scheme.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
FILES = [
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        dest = args.outdir / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"exists {dest} ({dest.stat().st_size} bytes)")
            continue
        url = f"{BASE}/{name}"
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, dest)
        print(f"wrote {dest} ({dest.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
