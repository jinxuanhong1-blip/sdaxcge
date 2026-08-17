#!/usr/bin/env python3
"""Download processed matrices used by analyze.py (not committed)."""

from __future__ import annotations

import urllib.request
from pathlib import Path

OUT = Path("/tmp/trop2_extra")
OUT.mkdir(parents=True, exist_ok=True)
AE = OUT / "emtab16433"
AE.mkdir(exist_ok=True)

FILES = {
    OUT / "GSE311016_gene_fpkm.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311016/suppl/GSE311016_gene_fpkm.txt.gz",
    OUT / "GSE311016_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311016/matrix/GSE311016_series_matrix.txt.gz",
    OUT / "GSE304294_gene_fpkm.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304294/suppl/GSE304294_gene_fpkm.txt.gz",
    OUT / "GSE304294_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304294/matrix/GSE304294_series_matrix.txt.gz",
    AE / "treatment_HD4246_SG_metadata.tsv":
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-16433/treatment_HD4246_SG_metadata.tsv",
    AE / "treatment_HD4246_SG_features.tsv":
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-16433/treatment_HD4246_SG_features.tsv",
    AE / "treatment_HD4246_SG_counts.txt":
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-16433/treatment_HD4246_SG_counts.txt",
}


def main() -> None:
    for dest, url in FILES.items():
        if dest.exists() and dest.stat().st_size > 1000:
            print("have", dest, dest.stat().st_size)
            continue
        print("get", url)
        urllib.request.urlretrieve(url, dest)
        print(" wrote", dest, dest.stat().st_size)


if __name__ == "__main__":
    main()
