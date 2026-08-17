#!/usr/bin/env python3
"""Download public processed UMIs for GSE131907 + GSE205335.

Skip the 2.86 GB GSE131907 log2TPM text matrix and all EGA/FASTQ.
Do not download GSE207422 or GSE148071.
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FILES = {
    "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "gse131907/GSE131907_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
        "GSE131907_series_matrix.txt.gz"
    ),
    "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "gse205335/GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
        "GSE205335_family.soft.gz"
    ),
}

SKIP = {
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": "2.86 GB; UMI is sufficient",
    "EGAD00001005054": "controlled-access FASTQ",
    "EGAD00001008703": "controlled-access FASTQ",
    "GSE207422": "not the winning pair; do not merge",
    "GSE148071": "not the winning pair; do not merge",
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("/tmp/geo_winpair"))
    args = parser.parse_args()
    for name, url in FILES.items():
        fetch(url, args.out / name)
    print("skipped:", SKIP)


if __name__ == "__main__":
    main()
