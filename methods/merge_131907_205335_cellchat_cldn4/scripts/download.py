#!/usr/bin/env python3
"""Download public processed UMIs for GSE131907 + GSE205335.

Skip the 2.86 GB GSE131907 log2TPM text matrix and all EGA/FASTQ.
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
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_merge_131907_205335"))
    args = p.parse_args()
    for name, url in FILES.items():
        fetch(url, args.out / name)


if __name__ == "__main__":
    main()
