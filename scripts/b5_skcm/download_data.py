#!/usr/bin/env python3
"""Download the public files used by B5_SKCM. Re-runnable; skips files that already exist."""

from __future__ import annotations

import urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parents[2] / "data" / "raw"
URLS = {
    "GSE78220_PatientFPKM.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx",
    "GSE78220_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/GSE78220_series_matrix.txt.gz",
    "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
    "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
    "GSE91061_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz",
    "GSE115821_MGH_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115821/suppl/GSE115821_MGH_counts.csv.gz",
    "GSE115821_family.soft.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115821/soft/GSE115821_family.soft.gz",
    "Homo_sapiens.gene_info.gz":
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "TCGA-SKCM.star_tpm.tsv.gz":
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.star_tpm.tsv.gz",
    "gencode.v36.annotation.gtf.gene.probemap":
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
}


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        dest = RAW / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"skip {name} ({dest.stat().st_size} bytes)")
            continue
        print(f"get  {name}")
        urllib.request.urlretrieve(url, dest)
        print(f"     {dest.stat().st_size} bytes")


if __name__ == "__main__":
    main()
