#!/usr/bin/env python3
"""Download public GEO processed files. No EGA / controlled raw."""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
        "suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    ),
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
        "suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    ),
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/geo"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
            continue
        print(f"GET {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
        print(f"  wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
