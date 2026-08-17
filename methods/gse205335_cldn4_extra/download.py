#!/usr/bin/env python3
"""Download open GEO processed files for GSE205335. No EGA raw."""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FILES = {
    "GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
        "suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
        "suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
        "soft/GSE205335_family.soft.gz"
    ),
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/gse205335"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
            continue
        print(f"get {name}", flush=True)
        urllib.request.urlretrieve(url, dest)
        print(f"  wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
