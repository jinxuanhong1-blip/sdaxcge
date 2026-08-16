#!/usr/bin/env python3
"""Download public processed GEO files used by methods/scrna_milo."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

GSE207422 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
GSE241934 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl"

FILES = {
    "GSE207422": [
        f"{GSE207422}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        f"{GSE207422}/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
    ],
    "GSE241934": [
        f"{GSE241934}/GSE241934_IIT_Matrix.mtx.gz",
        f"{GSE241934}/GSE241934_IIT_barcodes.tsv.gz",
        f"{GSE241934}/GSE241934_IIT_features.tsv.gz",
        f"{GSE241934}/GSE241934_IIT_Meta.txt.gz",
        f"{GSE241934}/GSE241934_Real_Matrix.mtx.gz",
        f"{GSE241934}/GSE241934_RWC_barcodes.tsv.gz",
        f"{GSE241934}/GSE241934_RWC_features.tsv.gz",
        f"{GSE241934}/GSE241934_Real_Meta.txt.gz",
    ],
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data"))
    ap.add_argument("--dataset", choices=["GSE207422", "GSE241934", "all"], default="all")
    args = ap.parse_args()
    keys = ["GSE207422", "GSE241934"] if args.dataset == "all" else [args.dataset]
    for key in keys:
        for url in FILES[key]:
            fetch(url, args.datadir / key / Path(url).name)


if __name__ == "__main__":
    main()
