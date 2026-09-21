#!/usr/bin/env python3
"""Download public GSE316655 10x MTX/TSV (no SRA)."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

FTP = "https://ftp.ncbi.nlm.nih.gov/geo"

FILES = [
    (
        f"{FTP}/series/GSE316nnn/GSE316655/matrix/GSE316655_series_matrix.txt.gz",
        "GSE316655_series_matrix.txt.gz",
    ),
    (
        f"{FTP}/series/GSE316nnn/GSE316655/suppl/filelist.txt",
        "filelist.txt",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457798/suppl/GSM9457798_ND_anti_B2_barcodes.tsv.gz",
        "GSM9457798_ND_anti_B2_barcodes.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457798/suppl/GSM9457798_ND_anti_B2_features.tsv.gz",
        "GSM9457798_ND_anti_B2_features.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457798/suppl/GSM9457798_ND_anti_B2_matrix.mtx.gz",
        "GSM9457798_ND_anti_B2_matrix.mtx.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457799/suppl/GSM9457799_ND_CTR_barcodes.tsv.gz",
        "GSM9457799_ND_CTR_barcodes.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457799/suppl/GSM9457799_ND_CTR_features.tsv.gz",
        "GSM9457799_ND_CTR_features.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457799/suppl/GSM9457799_ND_CTR_matrix.mtx.gz",
        "GSM9457799_ND_CTR_matrix.mtx.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457800/suppl/GSM9457800_PB_anti_B2_barcodes.tsv.gz",
        "GSM9457800_PB_anti_B2_barcodes.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457800/suppl/GSM9457800_PB_anti_B2_features.tsv.gz",
        "GSM9457800_PB_anti_B2_features.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457800/suppl/GSM9457800_PB_anti_B2_matrix.mtx.gz",
        "GSM9457800_PB_anti_B2_matrix.mtx.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457801/suppl/GSM9457801_PB_CTR_barcodes.tsv.gz",
        "GSM9457801_PB_CTR_barcodes.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457801/suppl/GSM9457801_PB_CTR_features.tsv.gz",
        "GSM9457801_PB_CTR_features.tsv.gz",
    ),
    (
        f"{FTP}/samples/GSM9457nnn/GSM9457801/suppl/GSM9457801_PB_CTR_matrix.mtx.gz",
        "GSM9457801_PB_CTR_matrix.mtx.gz",
    ),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("data/raw"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for url, name in FILES:
        dest = args.outdir / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"exists {dest} ({dest.stat().st_size} bytes)")
            continue
        print(f"GET {url}")
        urllib.request.urlretrieve(url, dest)
        print(f"  -> {dest} ({dest.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
