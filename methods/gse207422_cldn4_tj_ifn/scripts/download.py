#!/usr/bin/env python3
"""Download public GEO supplementary files for GSE207422."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
}


def fetch(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(url, args.workdir / name)
    print("download complete", flush=True)


if __name__ == "__main__":
    main()
