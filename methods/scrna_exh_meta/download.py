#!/usr/bin/env python3
"""Download public GEO files used by methods/scrna_exh_meta.

Does not invent accessions. Large matrices stay under --datadir (default /tmp).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
    "GSE205335_Lung_IO_UMI_matrix.rds.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
    "GSE205335_Lung_IO_CellIdentity.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
    "GSE241934_IIT_Matrix.mtx.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
    "GSE241934_IIT_Meta.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
    "GSE241934_IIT_barcodes.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
    "GSE241934_IIT_features.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
    "GSE241934_Real_Matrix.mtx.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
    "GSE241934_Real_Meta.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
    "GSE241934_RWC_barcodes.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
    "GSE241934_RWC_features.tsv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
    "GSE291670_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291670/suppl/GSE291670_RAW.tar",
    "GSE233203_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233203/suppl/GSE233203_RAW.tar",
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    part = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "curl", "-fL", "--retry", "5", "--retry-delay", "4",
        "-A", "Mozilla/5.0", "-o", str(part), url,
    ]
    print("GET", dest.name, flush=True)
    subprocess.check_call(cmd)
    part.rename(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/scrna_exh_meta"))
    args = ap.parse_args()
    for name, url in FILES.items():
        try:
            fetch(url, args.datadir / name)
        except subprocess.CalledProcessError as exc:
            print(f"FAILED {name}: {exc}", file=sys.stderr)
            raise
    print("download complete", args.datadir)


if __name__ == "__main__":
    main()
