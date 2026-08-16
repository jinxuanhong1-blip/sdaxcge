#!/usr/bin/env python3
"""Download public GEO files that fit the 2 GB processed-file budget."""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIZE_BUDGET = 2 * 1024**3

FILES = [
    {
        "cohort": "GSE207422",
        "file": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE207422",
        "file": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_IIT_Meta.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_Real_Meta.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_IIT_features.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_IIT_barcodes.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_RWC_features.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_RWC_barcodes.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_IIT_Matrix.mtx.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE241934",
        "file": "GSE241934_Real_Matrix.mtx.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
        "wanted": True,
    },
    {
        "cohort": "GSE253013",
        "file": "GSE253013_all_luad_garnett_temp.rds.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/suppl/GSE253013_all_luad_garnett_temp.rds.gz",
        "wanted": False,
        "reason": "9.3 GB > 2 GB budget; use public-derived patient tables",
    },
]


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    rows = []
    for spec in FILES:
        dest = args.root / "data" / spec["cohort"] / spec["file"]
        size = dest.stat().st_size if dest.exists() else head_size(spec["url"])
        within = True if size is None else size < SIZE_BUDGET
        if spec["wanted"] and within:
            action = "download"
        elif not spec["wanted"]:
            action = "skip_catalog"
        else:
            action = "skip_over_2gb"
        if dest.exists() and dest.stat().st_size > 0 and action == "download":
            action = "already_present"
            size = dest.stat().st_size
        rows.append({**spec, "size_bytes": size, "action": action, "within_budget": within})
        print(f"{action.upper():16} {spec['file']} ({size})")
        if action == "download":
            fetch(spec["url"], dest)
    (args.root / "results").mkdir(parents=True, exist_ok=True)
    (args.root / "results" / "download_manifest.json").write_text(json.dumps(rows, indent=2, default=str))


if __name__ == "__main__":
    main()
