#!/usr/bin/env python3
"""Download public GSE131907 annotation + raw UMI (Kim 2020).

Skip the 2.86 GB log2TPM text matrix and EGA/FASTQ.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
FILES = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        f"{BASE}/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907_series_matrix.txt.gz": (
        f"{BASE}/matrix/GSE131907_series_matrix.txt.gz"
    ),
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        f"{BASE}/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "sdaxcge-gse131907-slingshot-real-cldn4/1.0 "
                "(public GEO extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
            )
        },
    )
    with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.suffix == ".partial":
            continue
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                h.update(chunk)
        rows.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": h.hexdigest(),
            }
        )
    out = root / "DOWNLOAD_MANIFEST.json"
    out.write_text(
        json.dumps(
            {
                "accession": "GSE131907",
                "paper": "Kim et al. Nat Commun 2020 PMID 32385277",
                "skipped": [
                    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
                    "EGA FASTQ",
                    "GSE207422",
                    "GSE205335",
                ],
                "files": rows,
            },
            indent=2,
        )
    )
    print(f"manifest {out}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/gse131907_slingshot_data"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(url, args.out / name)
    manifest(args.out)


if __name__ == "__main__":
    main()
