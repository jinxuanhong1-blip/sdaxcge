#!/usr/bin/env python3
"""Download public processed UMIs for GSE205335 + GSE207422.

Skip GSE148071, GSE179994, EGA/FASTQ, GSA-Human HRA001033, and any file ≥2 GB.
Not a 7-pool.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

SIZE_BUDGET = 2 * 1024**3

FILES = {
    "gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "gse205335/GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
        "GSE205335_family.soft.gz"
    ),
    "gse207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    ),
    "gse207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
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
    for path in sorted(root.rglob("*")):
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
                "file": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": h.hexdigest(),
            }
        )
    out = root / "DOWNLOAD_MANIFEST.json"
    out.write_text(
        json.dumps(
            {
                "accessions": ["GSE205335", "GSE207422"],
                "skipped": [
                    "GSE148071",
                    "GSE179994",
                    "7-pool / Harmony joint object",
                    "EGA FASTQ",
                    "HRA001033",
                ],
                "size_budget_bytes": SIZE_BUDGET,
                "files": rows,
            },
            indent=2,
        )
    )
    print(f"manifest {out}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/ici_pair_205335_207422"))
    args = p.parse_args()
    for name, url in FILES.items():
        fetch(url, args.out / name)
        dest = args.out / name
        if dest.exists() and dest.stat().st_size >= SIZE_BUDGET:
            raise SystemExit(f"{dest} exceeds 2 GB budget")
    manifest(args.out)


if __name__ == "__main__":
    main()
