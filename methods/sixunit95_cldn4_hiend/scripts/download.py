#!/usr/bin/env python3
"""Download public processed files for the locked 6-unit combo.

Skip the 9.3 GB GSE253013 RDS. Skip GSE131907 2.86 GB log2TPM.
No FASTQ / EGA / SRA.
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

# Hard skip: GEO lists only this processed matrix for GSE253013.
SKIP_9GB = {
    "GSE253013_all_luad_garnett_temp.rds.gz": 9_966_268_875,
}

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
    "GSE205335": {
        "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
        ),
        "GSE205335_Lung_IO_CellIdentity.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
        ),
        "GSE205335_family.soft.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
            "soft/GSE205335_family.soft.gz"
        ),
    },
    "GSE291670": {
        "GSE291670_RAW.tar": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291670/"
            "suppl/GSE291670_RAW.tar"
        ),
    },
    "GSE131907": {
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
        ),
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        ),
        "GSE131907_series_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/"
            "matrix/GSE131907_series_matrix.txt.gz"
        ),
    },
    "GSE325414": {
        "GSE325414_barcodes.tsv.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325414/"
            "suppl/GSE325414_barcodes.tsv.gz"
        ),
        "GSE325414_features.tsv.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325414/"
            "suppl/GSE325414_features.tsv.gz"
        ),
        "GSE325414_matrix.mtx.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325414/"
            "suppl/GSE325414_matrix.mtx.gz"
        ),
        "GSE325414_metadata_individual_cells.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325414/"
            "suppl/GSE325414_metadata_individual_cells.txt.gz"
        ),
    },
}

MAX_BYTES = 2_000_000_000  # refuse anything ≥ 2 GB


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"have {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    if dest.name in SKIP_9GB:
        raise SystemExit(f"refusing 9GB file {dest.name}")
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    size = tmp.stat().st_size
    if size >= MAX_BYTES:
        tmp.unlink()
        raise SystemExit(f"refusing {dest.name}: {size} bytes ≥ 2 GB")
    tmp.replace(dest)
    print(f"wrote {dest} ({size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--out",
        type=Path,
        default=Path("/tmp/sixunit95_raw"),
    )
    args = p.parse_args()
    for unit, files in FILES.items():
        for name, url in files.items():
            fetch(url, args.out / unit / name)
    print("skipped GSE253013 9.3 GB RDS (only public processed matrix)", flush=True)


if __name__ == "__main__":
    main()
