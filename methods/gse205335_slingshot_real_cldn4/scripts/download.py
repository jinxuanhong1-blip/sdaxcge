#!/usr/bin/env python3
"""Download public processed UMIs for GSE205335 (Ahn/Lee ICI).

Skip EGA FASTQ (EGAD00001008703). No other accessions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

FILES = {
    "GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
        "GSE205335_family.soft.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.suffix == ".partial":
            continue
        if path.name == "DOWNLOAD_MANIFEST.json":
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
                "accessions": ["GSE205335"],
                "citation": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                "skipped": ["EGA FASTQ EGAD00001008703", "GSE131907", "GSE207422"],
                "files": rows,
            },
            indent=2,
        )
    )
    print(f"manifest {out}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/gse205335_slingshot_real"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(url, args.out / name)
    manifest(args.out)


if __name__ == "__main__":
    main()
