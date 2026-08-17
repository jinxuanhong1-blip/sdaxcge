#!/usr/bin/env python3
"""Download public processed UMI for the triple merge. Skip files >2 GB."""
from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

MAX_BYTES = 2 * 1024**3

FILES = {
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
    "GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
}

SKIP = {
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": 3.066e9,
}


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        cl = resp.headers.get("Content-Length")
    return int(cl) if cl else None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/geo"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
            continue
        size = head_size(url)
        if size is not None and size > MAX_BYTES:
            print(f"SKIP {name} ({size/1e9:.2f} GB > 2 GB)", flush=True)
            continue
        print(f"GET {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
        print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    for name, size in SKIP.items():
        print(f"SKIP listed {name} ({size/1e9:.2f} GB)", flush=True)
    tar = args.out / "GSE148071_RAW.tar"
    dest = args.out / "GSE148071_files"
    if tar.exists():
        dest.mkdir(parents=True, exist_ok=True)
        if not any(dest.glob("*_exp.txt.gz")):
            with tarfile.open(tar, "r") as tf:
                tf.extractall(dest)
        print(f"extracted {len(list(dest.glob('*exp.txt.gz')))} GSE148071 matrices", flush=True)


if __name__ == "__main__":
    main()
