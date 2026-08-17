#!/usr/bin/env python3
"""Download public GEO processed files and the NicheNet-v2 ligand–target RDS.

Skipped: EGA raw FASTQ, GSE131907 2.86 GB log2TPM, GSE253013 9.3 GB RDS.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_io import CACHE, DATA, download, log  # noqa: E402

FILES = [
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        ),
        "dest": CACHE / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "min_bytes": 100_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "matrix/GSE131907_series_matrix.txt.gz"
        ),
        "dest": CACHE / "GSE131907_series_matrix.txt.gz",
        "min_bytes": 2_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
        ),
        "dest": CACHE / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "min_bytes": 100_000_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
        ),
        "dest": CACHE / "GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "min_bytes": 100_000_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
        ),
        "dest": CACHE / "GSE205335_Lung_IO_CellIdentity.txt.gz",
        "min_bytes": 10_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "soft/GSE205335_family.soft.gz"
        ),
        "dest": CACHE / "GSE205335_family.soft.gz",
        "min_bytes": 10_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
        ),
        "dest": CACHE / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "min_bytes": 50_000_000,
    },
    {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
        ),
        "dest": DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "min_bytes": 5_000,
    },
    {
        "url": (
            "https://zenodo.org/records/7074291/files/"
            "ligand_target_matrix_nsga2r_final.rds?download=1"
        ),
        "dest": CACHE / "ligand_target_matrix_nsga2r_final.rds",
        "min_bytes": 200_000_000,
    },
]


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    for spec in FILES:
        download(spec["url"], spec["dest"], spec["min_bytes"])
    log("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
