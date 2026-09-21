#!/usr/bin/env python3
"""Download public processed matrices for concordant-4 CytoTRACE2.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
No GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.
No FASTQ. No GSE131907 log2TPM text. No author 36.5 GB GSE123902 H5.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

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


def curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = [
        "curl", "-fL", "--retry", "5", "--retry-delay", "8", "--retry-all-errors",
        "-C", "-", "-o", str(tmp), url,
    ]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/concordant4_ct2"))
    args = p.parse_args()
    for rel, url in FILES.items():
        curl_download(url, args.out / rel)


if __name__ == "__main__":
    main()
