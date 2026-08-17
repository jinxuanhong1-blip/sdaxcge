#!/usr/bin/env python3
"""Download public processed files for the winning-pair propeller + CellPhoneDB run.

GSE131907: author annotation + raw UMI text (skip 2.86 GB log2TPM).
GSE205335: author identity + processed UMI RDS (skip EGA raw).
CellPhoneDB v5: interaction / gene / complex tables.
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FILES = [
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    ),
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    ),
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz",
    ),
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz",
    ),
    (
        "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv",
        "cellphonedb/interaction_input.csv",
    ),
    (
        "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv",
        "cellphonedb/gene_input.csv",
    ),
    (
        "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv",
        "cellphonedb/complex_input.csv",
    ),
]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"OK {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp"))
    args = ap.parse_args()
    for url, rel in FILES:
        fetch(url, args.workdir / rel)


if __name__ == "__main__":
    main()
