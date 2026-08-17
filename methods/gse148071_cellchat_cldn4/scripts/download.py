#!/usr/bin/env python3
"""Download public GSE148071 per-sample UMI matrices + series matrix (GEO)."""
from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

FILES = {
    "GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "GSE148071_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/"
        "GSE148071_series_matrix.txt.gz"
    ),
}


def extract_tar(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.glob("*_exp.txt.gz")) or any(dest.rglob("*_exp.txt.gz")):
        return
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data",
    )
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        else:
            print(f"GET {url}", flush=True)
            urllib.request.urlretrieve(url, dest)
            print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    extract_tar(args.out / "GSE148071_RAW.tar", args.out / "GSE148071_files")
    n = len(list((args.out / "GSE148071_files").rglob("*_exp.txt.gz")))
    print(f"extracted exp matrices: {n}", flush=True)


if __name__ == "__main__":
    main()
