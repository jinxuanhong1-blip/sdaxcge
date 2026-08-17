#!/usr/bin/env python3
"""Download public processed matrices for the GSE123902 + GSE131907 pair.

Gates:
- GSE123902 GEO RAW.tar (dense unnormalized CSVs) is 90.4 MB.
- Author annotated H5 (36.5 GB) is skipped (>2 GB).
- GSE131907 raw UMI txt.gz is 389.8 MB. The 2.9 GB log2TPM text is skipped.
- GSE123904 SuperSeries RAW.tar (211.6 MB) includes mouse GSE123903; human
  matrices are taken from GSE123902.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "gse123902/GSE123902_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
    ),
    "gse123902/GSE123902_GEO_README.rtf": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_GEO_README.rtf"
    ),
    "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "gse131907/GSE131907_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
        "GSE131907_series_matrix.txt.gz"
    ),
}

SKIPPED = {
    "PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5": {
        "url": "https://s3.amazonaws.com/dp-lab-data-public/lung-development-cancer-progression/PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5",
        "bytes": 36489164157,
        "reason": "author annotated H5 is 36.5 GB (>2 GB public-processed gate)",
    },
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
        "bytes": 3110000000,
        "reason": "2.9 GB log2TPM text; UMI matrix exists and is <2 GB",
    },
}


def curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = ["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(tmp), url]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_123902_131907"))
    args = p.parse_args()
    for rel, url in FILES.items():
        curl_download(url, args.out / rel)
    print("skipped:")
    for name, rec in SKIPPED.items():
        print(f"  {name}: {rec['reason']}")


if __name__ == "__main__":
    main()
