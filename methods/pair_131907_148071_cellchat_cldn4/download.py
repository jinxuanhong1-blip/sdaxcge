#!/usr/bin/env python3
"""Download public GEO processed files for GSE131907 + GSE148071.

No EGA raw. Skip the 2.86 GB GSE131907 log2TPM text matrix.
"""
from __future__ import annotations

import argparse
import tarfile
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
    "gse148071/GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "gse148071/GSE148071_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/"
        "GSE148071_series_matrix.txt.gz"
    ),
}


def extract_tar(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.rglob("*_exp.txt.gz")):
        return
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_131907_148071"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for rel, url in FILES.items():
        dest = args.out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
        else:
            print(f"GET {url}", flush=True)
            tmp = dest.with_suffix(dest.suffix + ".partial")
            urllib.request.urlretrieve(url, tmp)
            tmp.replace(dest)
            print(f"  wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    extract_tar(args.out / "gse148071/GSE148071_RAW.tar", args.out / "gse148071/GSE148071_files")
    n = len(list((args.out / "gse148071/GSE148071_files").rglob("*_exp.txt.gz")))
    print(f"extracted GSE148071 exp matrices: {n}", flush=True)


if __name__ == "__main__":
    main()
