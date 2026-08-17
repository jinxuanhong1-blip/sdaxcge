#!/usr/bin/env python3
"""Download open processed files for the triple-merge CLDN4 AUCell slice.

GSE131907 / GSE205335: GEO processed UMI + author identity.
GSE148071: TISCH2 expression.h5 + CellMetainfo (GEO 10x is per-sample RAW).
Does not fetch EGA raw (EGAD00001005054 / EGAD00001008703).
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "gse131907": {
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
    "gse205335": {
        "GSE205335_Lung_IO_CellIdentity.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
        ),
        "GSE205335_family.soft.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "soft/GSE205335_family.soft.gz"
        ),
        "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
        ),
    },
    "gse148071": {
        "NSCLC_GSE148071_CellMetainfo_table.tsv": (
            "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/"
            "NSCLC_GSE148071_CellMetainfo_table.tsv"
        ),
        "NSCLC_GSE148071_expression.h5": (
            "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/"
            "NSCLC_GSE148071_expression.h5"
        ),
    },
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"get {dest.name}", flush=True)
    subprocess.run(
        ["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(dest), url],
        check=True,
    )
    print(f"  wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("/tmp"))
    args = p.parse_args()
    for cohort, files in FILES.items():
        out = args.root / cohort
        for name, url in files.items():
            fetch(url, out / name)


if __name__ == "__main__":
    main()
