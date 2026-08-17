#!/usr/bin/env python3
"""Download public GSE207422 processed UMI + sample metadata (GEO only)."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    ),
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    part = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "curl", "-fL", "--retry", "5", "--retry-delay", "4",
        "-A", "Mozilla/5.0", "-o", str(part), url,
    ]
    print("GET", dest.name, flush=True)
    subprocess.check_call(cmd)
    part.rename(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse207422_cldn4_exh"))
    args = ap.parse_args()
    for name, url in FILES.items():
        try:
            fetch(url, args.datadir / name)
        except subprocess.CalledProcessError as exc:
            print(f"FAILED {name}: {exc}", file=sys.stderr)
            raise
    print("download complete", args.datadir)


if __name__ == "__main__":
    main()
