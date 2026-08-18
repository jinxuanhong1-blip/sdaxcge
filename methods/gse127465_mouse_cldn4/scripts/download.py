#!/usr/bin/env python3
"""Download GSE127465 *mouse* processed files only. Human MTX is not fetched."""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
FILES = [
    "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz",
    "GSE127465_mouse_cell_metadata_15939x12.tsv.gz",
    "GSE127465_gene_names_mouse_28205.tsv.gz",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("/tmp/gse127465_mouse"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"exists {dest} ({dest.stat().st_size} B)")
            continue
        url = BASE + name
        print(f"GET {url}")
        urllib.request.urlretrieve(url, dest)
        print(f"wrote {dest} ({dest.stat().st_size} B)")


if __name__ == "__main__":
    main()
