#!/usr/bin/env python3
"""Download public GSE131907 processed files needed for CLDN4 vs potency."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

DATA = Path("/tmp/gse131907")
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    "GSE131907_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
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
