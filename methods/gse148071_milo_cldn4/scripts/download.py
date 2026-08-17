#!/usr/bin/env python3
"""Download public processed GEO files for GSE148071 (Wu et al. 2021)."""

from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071"
FILES = {
    "GSE148071_RAW.tar": f"{BASE}/suppl/GSE148071_RAW.tar",
    "GSE148071_series_matrix.txt.gz": f"{BASE}/matrix/GSE148071_series_matrix.txt.gz",
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)")


def extract_tar(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    existing = list(dest.glob("*_exp.txt.gz")) or list(dest.rglob("*_exp.txt.gz"))
    if existing:
        print(f"already extracted {len(existing)} exp files in {dest}")
        return
    print(f"extract {tar_path} -> {dest}")
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    n = len(list(dest.glob("*_exp.txt.gz")) or list(dest.rglob("*_exp.txt.gz")))
    print(f"extracted {n} exp files")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE148071"))
    ap.add_argument("--skip-extract", action="store_true")
    args = ap.parse_args()
    for name, url in FILES.items():
        fetch(url, args.datadir / name)
    if not args.skip_extract:
        extract_tar(args.datadir / "GSE148071_RAW.tar", args.datadir / "files")


if __name__ == "__main__":
    main()
