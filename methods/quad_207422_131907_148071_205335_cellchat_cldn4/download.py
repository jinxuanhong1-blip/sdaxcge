#!/usr/bin/env python3
"""Download public processed UMIs for the CLDN4 quad merge.

Cohorts: GSE207422 + GSE131907 + GSE148071 + GSE205335.
Skip any file >2 GB. Skip GSE131907 log2TPM (2.86 GB) and all EGA raw.
"""
from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

MAX_BYTES = 2 * 1024 * 1024 * 1024

FILES = {
    "GSE207422": {
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
        ),
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
        ),
    },
    "GSE131907": {
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        ),
        "GSE131907_series_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "matrix/GSE131907_series_matrix.txt.gz"
        ),
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
        ),
    },
    "GSE148071": {
        "GSE148071_RAW.tar": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/"
            "suppl/GSE148071_RAW.tar"
        ),
        "GSE148071_series_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/"
            "matrix/GSE148071_series_matrix.txt.gz"
        ),
    },
    "GSE205335": {
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
    },
}


def content_length(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        value = resp.headers.get("Content-Length")
    return int(value) if value else None


def fetch(url: str, dest: Path) -> None:
    size = content_length(url)
    if size is not None and size > MAX_BYTES:
        print(f"SKIP >2GB {dest.name} ({size} bytes)", flush=True)
        return
    if dest.exists() and dest.stat().st_size > 1000:
        if size is None or dest.stat().st_size >= size * 0.98:
            print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
            return
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    if tmp.stat().st_size > MAX_BYTES:
        tmp.unlink()
        print(f"SKIP downloaded >2GB {dest.name}", flush=True)
        return
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def extract_gse148071(data_dir: Path) -> None:
    tar_path = data_dir / "GSE148071_RAW.tar"
    dest = data_dir / "GSE148071_files"
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.rglob("*_exp.txt.gz")):
        n = len(list(dest.rglob("*_exp.txt.gz")))
        print(f"GSE148071 already extracted ({n} exp matrices)", flush=True)
        return
    if not tar_path.exists():
        print("no GSE148071_RAW.tar to extract", flush=True)
        return
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    n = len(list(dest.rglob("*_exp.txt.gz")))
    print(f"extracted GSE148071 exp matrices: {n}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("/tmp/quad_cldn4_geo"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for cohort, files in FILES.items():
        print(f"== {cohort} ==", flush=True)
        for name, url in files.items():
            fetch(url, args.out / name)
    extract_gse148071(args.out)


if __name__ == "__main__":
    main()
