#!/usr/bin/env python3
"""Download public GSE127465 human inDrops processed files (Zilionis 2019)."""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
FILES = {
    "GSE127465_human_cell_metadata_54773x25.tsv.gz": 4092351,
    "GSE127465_gene_names_human_41861.tsv.gz": 104690,
    "GSE127465_human_counts_normalized_54773x41861.mtx.gz": 528303938,
}


def fetch(url: str, dest: Path, expect: int | None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        if expect is None or dest.stat().st_size == expect:
            print(f"exists {dest} ({dest.stat().st_size} B)", flush=True)
            return dest
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    print(f"wrote {dest} ({dest.stat().st_size} B)", flush=True)
    return dest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/gse127465_slingshot")
    args = ap.parse_args()
    out = Path(args.out)
    for name, nbytes in FILES.items():
        fetch(BASE + name, out / name, nbytes)


if __name__ == "__main__":
    main()
