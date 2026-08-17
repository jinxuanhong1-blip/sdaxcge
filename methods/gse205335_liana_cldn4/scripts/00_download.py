#!/usr/bin/env python3
"""Download public GSE205335 processed UMI + identity + SOFT, and CellPhoneDB v5.

Public processed files only. EGA raw (EGAD00001008703) is not accessed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/"
CPDB = "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/"

FILES = [
    {
        "url": GEO + "suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "dest": "GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "role": "processed_umi",
    },
    {
        "url": GEO + "suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "dest": "GSE205335_Lung_IO_CellIdentity.txt.gz",
        "role": "author_identity",
    },
    {
        "url": GEO + "soft/GSE205335_family.soft.gz",
        "dest": "GSE205335_family.soft.gz",
        "role": "geo_soft",
    },
    {
        "url": CPDB + "interaction_input.csv",
        "dest": "cellphonedb/interaction_input.csv",
        "role": "cpdb_interactions",
    },
    {
        "url": CPDB + "gene_input.csv",
        "dest": "cellphonedb/gene_input.csv",
        "role": "cpdb_genes",
    },
    {
        "url": CPDB + "complex_input.csv",
        "dest": "cellphonedb/complex_input.csv",
        "role": "cpdb_complexes",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET  {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"OK   {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=Path, default=Path("/tmp/gse205335"))
    args = parser.parse_args()
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
    manifest = args.datadir / "file_manifest.json"
    manifest.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"wrote {manifest}", flush=True)
    print("skipped_raw_ega=EGAD00001008703", flush=True)


if __name__ == "__main__":
    main()
