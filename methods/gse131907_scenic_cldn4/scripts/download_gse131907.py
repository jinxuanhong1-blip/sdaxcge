#!/usr/bin/env python3
"""Download in-budget GSE131907 processed files.

Uses the public UMI matrix (0.38 GB) + author cell annotation.
Skips the 2.86 GB log2TPM text and EGA raw FASTQ.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
SIZE_BUDGET = 2 * 1024**3

FILES = {
    "GSE131907_Lung_Cancer_cell_annotation.txt.gz": SUPP + "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    "GSE131907_series_matrix.txt.gz": MATRIX + "GSE131907_series_matrix.txt.gz",
    "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": SUPP + "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
}
SKIP = {
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": "2.86 GB; UMI is sufficient",
    "EGAD00001005054": "controlled-access FASTQ; not a processed GEO supplement",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[cache] {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"[get] {url} -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/gse131907_scenic_cldn4/geo"))
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "accession": "GSE131907",
        "pmid": "32385277",
        "size_budget_gb": 2.0,
        "skipped": SKIP,
        "files": {},
    }
    for name, url in FILES.items():
        dest = out / name
        fetch(url, dest)
        size = dest.stat().st_size
        if size >= SIZE_BUDGET:
            raise SystemExit(f"{name} is {size} bytes; over 2 GB budget")
        manifest["files"][name] = {
            "url": url,
            "bytes": size,
            "sha256": sha256(dest),
        }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
