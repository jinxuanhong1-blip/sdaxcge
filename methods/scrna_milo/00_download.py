#!/usr/bin/env python3
"""Download public GSE207422 processed scRNA files (GEO + paper Table S1)."""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
ESM_BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs13073-023-01164-9/MediaObjects/"
)

FILES = [
    {
        "file": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": GEO_SUPPL + "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "role": "processed_umi",
    },
    {
        "file": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": GEO_SUPPL + "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "role": "sample_metadata",
    },
    {
        "file": "13073_2023_1164_MOESM1_ESM.xlsx",
        "url": ESM_BASE + "13073_2023_1164_MOESM1_ESM.xlsx",
        "role": "paper_table_s1_clinical",
    },
]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"PRESENT {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {dest.name}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"SAVED {dest.name} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for spec in FILES:
        dest = args.workdir / spec["file"]
        fetch(spec["url"], dest)
        manifest.append(
            {
                "file": spec["file"],
                "role": spec["role"],
                "bytes": dest.stat().st_size,
                "url": spec["url"],
            }
        )
    (args.workdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
