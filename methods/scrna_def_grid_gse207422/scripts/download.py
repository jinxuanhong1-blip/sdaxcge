#!/usr/bin/env python3
"""Download public GSE207422 UMI matrix + sample metadata (GEO only)."""
from __future__ import annotations

import urllib.request
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
FILES = [
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
]


def main() -> None:
    out = Path("data/GSE207422")
    out.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        dest = out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
            continue
        url = f"{GEO}/{name}"
        print(f"downloading {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
        print(f"  -> {dest} ({dest.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
