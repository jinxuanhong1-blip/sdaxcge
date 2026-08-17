#!/usr/bin/env python3
"""Download GSE293914 public series matrix + GSM8893263 MTX (no SRA)."""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FILES = {
    "GSE293914_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293914/matrix/"
        "GSE293914_series_matrix.txt.gz"
    ),
    "GSM8893263_barcodes.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8893nnn/GSM8893263/suppl/"
        "GSM8893263_barcodes.tsv.gz"
    ),
    "GSM8893263_features.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8893nnn/GSM8893263/suppl/"
        "GSM8893263_features.tsv.gz"
    ),
    "GSM8893263_matrix.mtx.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8893nnn/GSM8893263/suppl/"
        "GSM8893263_matrix.mtx.gz"
    ),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE293914"))
    args = ap.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.datadir / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"have {dest} ({dest.stat().st_size:,} bytes)")
            continue
        print(f"get {name}", flush=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)
        print(f"  -> {dest} ({dest.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
