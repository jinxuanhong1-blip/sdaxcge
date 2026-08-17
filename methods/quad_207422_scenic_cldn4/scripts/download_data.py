#!/usr/bin/env python3
"""Download the four public QUAD matrices used by this folder. Not committed."""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

FILES = {
    "GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    ),
    "GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    ),
    "GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ),
    "GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
        "GSE205335_Lung_IO_CellIdentity.txt.gz"
    ),
    "GSE205335/GSE205335_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/"
        "GSE205335_family.soft.gz"
    ),
    "GSE148071/NSCLC_GSE148071_expression.h5": (
        "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/"
        "NSCLC_GSE148071_expression.h5"
    ),
    "GSE148071/NSCLC_GSE148071_CellMetainfo_table.tsv": (
        "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/"
        "NSCLC_GSE148071_CellMetainfo_table.tsv"
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
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/quad_scenic"))
    args = ap.parse_args()
    manifest = {"files": {}}
    for rel, url in FILES.items():
        dest = args.outdir / rel
        fetch(url, dest)
        manifest["files"][rel] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
        }
    (args.outdir / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps({k: v["bytes"] for k, v in manifest["files"].items()}, indent=2))


if __name__ == "__main__":
    main()
