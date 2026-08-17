#!/usr/bin/env python3
"""Download public GSE131907 processed UMI + author annotation (GEO).

Skip the 2.86 GB log2TPM text matrix and EGA raw FASTQ (EGAD00001005054).
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"

FILES = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": SUPP + "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    "GSE131907_Lung_Cancer_Feature_Summary.xlsx": SUPP + "GSE131907_Lung_Cancer_Feature_Summary.xlsx",
    "GSE131907_series_matrix.txt.gz": MATRIX + "GSE131907_series_matrix.txt.gz",
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": SUPP + "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/gse131907"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
            continue
        print(f"GET {url}", flush=True)
        tmp = dest.with_suffix(dest.suffix + ".partial")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(dest)
        print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
