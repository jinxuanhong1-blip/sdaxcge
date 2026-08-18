#!/usr/bin/env python3
"""Download public processed matrices for the concordant-4 atlas.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
No GSE148071 / GSE127465 / GSE154826 / GSE207422. No FASTQ / EGA.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = Path("/tmp/geo_atlas")

FILES = {
    "GSE123902/GSE123902_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
    ),
    "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "GSE205335/GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
        "GSE205335_family.soft.gz"
    ),
    "GSE189357/GSE189357_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/"
        "GSE189357_RAW.tar"
    ),
}

SKIPPED = [
    "GSE148071",
    "GSE127465",
    "GSE154826",
    "GSE207422",
    "author 36.5 GB GSE123902 H5",
    "GSE131907 2.9 GB log2TPM text",
    "EGA FASTQ EGAD00001008703",
    "GSE189357 SRA FASTQ / GSE189487 spatial",
]


def curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = [
        "curl",
        "-fL",
        "--retry",
        "5",
        "--retry-delay",
        "8",
        "--retry-all-errors",
        "-C",
        "-",
        "-o",
        str(tmp),
        url,
    ]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()
    for rel, url in FILES.items():
        curl_download(url, args.out / rel)
    print("skipped:", ", ".join(SKIPPED))


if __name__ == "__main__":
    main()
