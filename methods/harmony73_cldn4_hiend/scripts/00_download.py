#!/usr/bin/env python3
"""Download public processed matrices for the Harmony n=73 cohorts.

Does **not** download GSE253013_all_luad_garnett_temp.rds.gz (~9 GB).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

DATA = Path("/tmp/harmony73_data")
URLS = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "GSE127465_human_cell_metadata_54773x25.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_human_cell_metadata_54773x25.tsv.gz"
    ),
    "GSE127465_gene_names_human_41861.tsv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_gene_names_human_41861.tsv.gz"
    ),
    "GSE127465_human_counts_normalized_54773x41861.mtx.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_human_counts_normalized_54773x41861.mtx.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"have {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return
    cmd = ["curl", "-L", "--retry", "5", "--retry-delay", "8", "-C", "-", "-o", str(dest), url]
    print("GET", url, "->", dest, flush=True)
    subprocess.check_call(cmd)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        fetch(url, DATA / name)
    tar = DATA / "GSE148071_RAW.tar"
    dest = DATA / "GSE148071_files"
    dest.mkdir(exist_ok=True)
    if not any(dest.glob("*_exp.txt.gz")):
        subprocess.check_call(["tar", "-xf", str(tar), "-C", str(dest)])
    print("downloads complete; GSE253013 RDS was not requested", flush=True)


if __name__ == "__main__":
    main()
