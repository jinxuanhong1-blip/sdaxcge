#!/usr/bin/env python3
"""Download public processed files for the triple merge + CellPhoneDB v5."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

FILES = [
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "dest": "GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "role": "gse131907_umi",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "dest": "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "role": "gse131907_annot",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
        "dest": "GSE131907/GSE131907_series_matrix.txt.gz",
        "role": "gse131907_series",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/GSE148071_RAW.tar",
        "dest": "GSE148071/GSE148071_RAW.tar",
        "role": "gse148071_tar",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/GSE148071_series_matrix.txt.gz",
        "dest": "GSE148071/GSE148071_series_matrix.txt.gz",
        "role": "gse148071_series",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "dest": "GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "role": "gse205335_umi",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "dest": "GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "role": "gse205335_ident",
    },
    {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
        "dest": "GSE205335/GSE205335_family.soft.gz",
        "role": "gse205335_soft",
    },
    {
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv",
        "dest": "cellphonedb/interaction_input.csv",
        "role": "cpdb_interactions",
    },
    {
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv",
        "dest": "cellphonedb/gene_input.csv",
        "role": "cpdb_genes",
    },
    {
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv",
        "dest": "cellphonedb/complex_input.csv",
        "role": "cpdb_complexes",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    print(f"GET  {url}", flush=True)
    wget = shutil.which("wget")
    if wget:
        subprocess.check_call(
            [
                wget,
                "-c",
                "--tries=8",
                "--waitretry=15",
                "--timeout=120",
                "--progress=dot:giga",
                "-O",
                str(dest),
                url,
            ]
        )
    else:
        tmp = dest.with_suffix(dest.suffix + ".partial")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(dest)
    print(f"OK   {dest} ({dest.stat().st_size} bytes)", flush=True)


def extract_tar(tar_path: Path, dest: Path) -> int:
    dest.mkdir(parents=True, exist_ok=True)
    existing = list(dest.rglob("*_exp.txt.gz"))
    if existing:
        return len(existing)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    return len(list(dest.rglob("*_exp.txt.gz")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    args = ap.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in FILES:
        dest = args.datadir / spec["dest"]
        fetch(spec["url"], dest)
        rows.append(
            {
                "file": spec["dest"],
                "role": spec["role"],
                "bytes": dest.stat().st_size,
                "sha256": sha256(dest),
                "url": spec["url"],
            }
        )
    n_exp = extract_tar(
        args.datadir / "GSE148071/GSE148071_RAW.tar",
        args.datadir / "GSE148071/GSE148071_files",
    )
    print(f"GSE148071 exp matrices: {n_exp}", flush=True)
    manifest = args.datadir / "file_manifest.json"
    manifest.write_text(json.dumps(rows, indent=2))
    print(f"wrote {manifest}", flush=True)


if __name__ == "__main__":
    main()
