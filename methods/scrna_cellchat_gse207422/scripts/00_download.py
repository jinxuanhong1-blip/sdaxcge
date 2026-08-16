#!/usr/bin/env python3
"""Download public GSE207422 processed UMI + sample sheet, and CellPhoneDB v5 tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
CPDB = "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/"

FILES = [
    {
        "url": GEO + "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "dest": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "role": "processed_umi",
    },
    {
        "url": GEO + "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "dest": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "role": "sample_metadata",
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
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET  {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"OK   {dest} ({dest.stat().st_size} bytes)")


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
    manifest = args.datadir / "file_manifest.json"
    manifest.write_text(json.dumps(rows, indent=2))
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
