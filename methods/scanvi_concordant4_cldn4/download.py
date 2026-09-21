#!/usr/bin/env python3
"""Download the four public concordant-4 matrices. Nothing else."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "GSE123902_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar",
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "GSE205335_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
    "GSE189357_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar",
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_c4"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest}", flush=True)
            continue
        subprocess.run(
            ["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-C", "-", "-o", str(dest), url],
            check=True,
        )


if __name__ == "__main__":
    main()
