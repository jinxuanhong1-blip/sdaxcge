#!/usr/bin/env python3
"""Download in-budget processed GEO files for GSE131907 + GSE205335.

Skips GSE131907 2.86 GB log2TPM text, EGA FASTQ, and GSE207422.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

SIZE_BUDGET = 2 * 1024**3

FILES = {
    "GSE131907": {
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
            "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        ),
        "GSE131907_series_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
            "GSE131907_series_matrix.txt.gz"
        ),
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
            "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
        ),
    },
    "GSE205335": {
        "GSE205335_Lung_IO_CellIdentity.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
            "GSE205335_Lung_IO_CellIdentity.txt.gz"
        ),
        "GSE205335_family.soft.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
            "GSE205335_family.soft.gz"
        ),
        "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
            "GSE205335_Lung_IO_UMI_matrix.rds.gz"
        ),
    },
}

SKIPPED = {
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": "2.86 GB; UMI is sufficient",
    "EGAD00001005054": "controlled-access FASTQ; not used",
    "EGAD00001008703": "controlled-access FASTQ; not used",
    "GSE207422": "out of scope (no GSE207422-only)",
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
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/merge_131907_205335_scenic/geo"))
    args = ap.parse_args()
    manifest = {"size_budget_gb": 2.0, "skipped": SKIPPED, "files": {}}
    for acc, files in FILES.items():
        acc_dir = args.outdir / acc
        acc_dir.mkdir(parents=True, exist_ok=True)
        for name, url in files.items():
            dest = acc_dir / name
            fetch(url, dest)
            size = dest.stat().st_size
            if size >= SIZE_BUDGET:
                raise SystemExit(f"{name} is {size} bytes; over 2 GB budget")
            manifest["files"][f"{acc}/{name}"] = {
                "url": url,
                "bytes": size,
                "sha256": sha256(dest),
            }
    (args.outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v["bytes"] for k, v in manifest["files"].items()}, indent=2))


if __name__ == "__main__":
    main()
