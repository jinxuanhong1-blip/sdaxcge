#!/usr/bin/env python3
"""Download processed public ICI-irAE datasets used in this slice.

Datasets (all human, immune checkpoint inhibitor context):
  - GSE319496 : whole-blood bulk RNA-seq, mRCC nivo+ipi, irAE Yes/No labels
                (processed raw counts CSV, ~0.7 MB). BULK irAE case/control.
  - GSE206300 : ircolitis colon *epithelial* single-nucleus RNA-seq h5ad
                (ICI colitis vs control). Epithelial -> TACSTD2/CLDN4 measurable.
  - GSE277136 : bronchoalveolar lavage fluid scRNA-seq of ICI-related
                *pneumonitis* (lung). Measurability of epithelial markers in BAL.

Everything lands under results/fable_irae/data/raw/ (gitignored). We keep the
on-disk processed footprint under 2 GB and never commit the large matrices.
"""
import hashlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_geo import fetch_url, series_suppl_url  # noqa: E402

RAW = Path(__file__).resolve().parents[2] / "results" / "fable_irae" / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

FILES = [
    ("GSE319496", "GSE319496_GEO_raw_counts_SampleID.csv.gz"),
    ("GSE206300", "GSE206300_ircolitis-tissue-epithelial.h5ad.gz"),
    ("GSE277136", "GSE277136_BLFplusNM.h5ad"),
]

# series matrices (small) for phenotype extraction
SERIES_MATRIX = ["GSE319496", "GSE277136"]


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  exists: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    print(f"  GET {url}")
    data = fetch_url(url, timeout=600, binary=True)
    dest.write_bytes(data)
    print(f"  saved {dest.name} ({len(data)/1e6:.1f} MB)")
    return dest


def main():
    for acc, fname in FILES:
        download(series_suppl_url(acc, fname), RAW / fname)
    for acc in SERIES_MATRIX:
        stub = acc[:-3] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/matrix/{acc}_series_matrix.txt.gz"
        download(url, RAW / f"{acc}_series_matrix.txt.gz")

    total = sum(p.stat().st_size for p in RAW.glob("*"))
    print(f"\nTotal raw footprint: {total/1e9:.2f} GB")
    manifest = RAW.parent / "download_manifest.tsv"
    with manifest.open("w") as f:
        f.write("file\tbytes\tmd5\n")
        for p in sorted(RAW.glob("*")):
            md5 = hashlib.md5(p.read_bytes()).hexdigest() if p.stat().st_size < 3e8 else "skipped_large"
            f.write(f"{p.name}\t{p.stat().st_size}\t{md5}\n")
    print(f"manifest -> {manifest}")


if __name__ == "__main__":
    main()
