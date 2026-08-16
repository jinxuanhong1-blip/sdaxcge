#!/usr/bin/env python3
"""Download public GEO supplementary files for GSE207422 and GSE241934."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

DATA = Path("/tmp/scrna_cytotrace")
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
    "GSE241934_IIT_Meta.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
    "GSE241934_IIT_barcodes.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
    "GSE241934_IIT_features.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
    "GSE241934_IIT_Matrix.mtx.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
    "GSE241934_Real_Meta.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
    "GSE241934_RWC_barcodes.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
    "GSE241934_RWC_features.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
    "GSE241934_Real_Matrix.mtx.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
}


def fetch(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = [
        "curl", "-L", "--retry", "8", "--retry-delay", "4",
        "--retry-all-errors", "-C", "-", "--fail", "-o", str(tmp), url,
    ]
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd)
    tmp.rename(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    for name, url in FILES.items():
        fetch(url, DATA / name)
    print("download complete", flush=True)


if __name__ == "__main__":
    sys.exit(main())
