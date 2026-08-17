#!/usr/bin/env python3
"""Download the public GSE207422 processed UMI matrix + sample metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    ),
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    ),
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
    urllib.request.urlretrieve(url, dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("data/scrna_scenic_cldn4/GSE207422"))
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"accession": "GSE207422", "files": {}}
    for name, url in FILES.items():
        dest = out / name
        fetch(url, dest)
        manifest["files"][name] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
        }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
