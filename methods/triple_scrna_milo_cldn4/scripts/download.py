#!/usr/bin/env python3
"""Download public processed GEO files for the triple Milo CLDN4 merge.

GSE131907 + GSE148071 + GSE205335 only. GSE207422 is a different agent
and is not fetched. EGA FASTQ is controlled-access and is not used.
"""

from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

FILES = {
    "GSE131907": [
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
            "GSE131907_series_matrix.txt.gz",
        ),
    ],
    "GSE148071": [
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/GSE148071_RAW.tar",
            "GSE148071_RAW.tar",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/GSE148071_series_matrix.txt.gz",
            "GSE148071_series_matrix.txt.gz",
        ),
    ],
    "GSE205335": [
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
            "GSE205335_Lung_IO_UMI_matrix.rds.gz",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
            "GSE205335_Lung_IO_CellIdentity.txt.gz",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
            "GSE205335_family.soft.gz",
        ),
    ],
}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    return dest


def extract_gse148071(tar_path: Path, dest: Path) -> None:
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
    ap.add_argument("--datadir", type=Path, default=Path("data"))
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["GSE131907", "GSE148071", "GSE205335"],
        choices=list(FILES),
    )
    args = ap.parse_args()
    for ds in args.datasets:
        ddir = args.datadir / ds
        for url, name in FILES[ds]:
            fetch(url, ddir / name)
        if ds == "GSE148071":
            extract_gse148071(ddir / "GSE148071_RAW.tar", ddir / "files")


if __name__ == "__main__":
    main()
