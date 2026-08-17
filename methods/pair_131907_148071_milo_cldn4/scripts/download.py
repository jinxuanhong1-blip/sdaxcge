#!/usr/bin/env python3
"""Download public GEO files for the GSE131907 + GSE148071 pair.

Not GSE205335. Not the triple. UMI / raw counts only.
"""

from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

GSE131907 = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE131907_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
        "GSE131907_series_matrix.txt.gz"
    ),
}

GSE148071 = {
    "GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "GSE148071_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/"
        "GSE148071_series_matrix.txt.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def extract_tar(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    existing = list(dest.glob("*_exp.txt.gz")) or list(dest.rglob("*_exp.txt.gz"))
    if existing:
        print(f"already extracted {len(existing)} exp files in {dest}", flush=True)
        return
    print(f"extract {tar_path} -> {dest}", flush=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    n = len(list(dest.glob("*_exp.txt.gz")) or list(dest.rglob("*_exp.txt.gz")))
    print(f"extracted {n} exp files", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=Path("data"))
    args = ap.parse_args()
    d131 = args.data_root / "GSE131907"
    d148 = args.data_root / "GSE148071"
    for name, url in GSE131907.items():
        fetch(url, d131 / name)
    for name, url in GSE148071.items():
        fetch(url, d148 / name)
    extract_tar(d148 / "GSE148071_RAW.tar", d148 / "files")


if __name__ == "__main__":
    main()
