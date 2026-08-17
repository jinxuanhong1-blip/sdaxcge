#!/usr/bin/env python3
"""Download public processed files for the ICI trio (skip anything ≥2 GB)."""
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIZE_BUDGET = 2 * 1024**3

FILES = [
    {
        "cohort": "GSE179994",
        "file": "GSE179994_Tcell.metadata.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_Tcell.metadata.tsv.gz",
        "role": "tcell_metadata_audit",
        "wanted": True,
    },
    {
        "cohort": "GSE179994",
        "file": "GSE179994_all.Tcell.rawCounts.rds.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz",
        "role": "tcell_only_umi_not_downloaded",
        "wanted": False,
    },
    {
        "cohort": "GSE207422",
        "file": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "role": "processed_umi",
        "wanted": True,
    },
    {
        "cohort": "GSE207422",
        "file": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "role": "sample_metadata",
        "wanted": True,
    },
    {
        "cohort": "GSE205335",
        "file": "GSE205335_Lung_IO_CellIdentity.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "role": "author_identity",
        "wanted": True,
    },
    {
        "cohort": "GSE205335",
        "file": "GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "role": "processed_umi",
        "wanted": True,
    },
    {
        "cohort": "GSE205335",
        "file": "GSE205335_family.soft.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
        "role": "sample_metadata",
        "wanted": True,
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
    print(f"GET {url} -> {dest}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    args = ap.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)

    rows = []
    for spec in FILES:
        dest = args.datadir / spec["cohort"] / spec["file"]
        size = head_size(spec["url"])
        within = True if size is None else size < SIZE_BUDGET
        if spec["wanted"] and within:
            action = "download"
        elif not within:
            action = "skip_over_2gb"
        else:
            action = "skip_not_usable"
        if dest.exists() and dest.stat().st_size > 0 and action == "download":
            action = "already_present"
            size = size or dest.stat().st_size
        rows.append(
            {
                "cohort": spec["cohort"],
                "file": spec["file"],
                "role": spec["role"],
                "size_bytes": size,
                "size_mb": None if size is None else round(size / 1024**2, 1),
                "within_2gb": within,
                "action": action,
                "url": spec["url"],
            }
        )
        if action == "download":
            fetch(spec["url"], dest)
        else:
            print(f"{action.upper()} {spec['file']} ({size} bytes)", flush=True)

    out = args.datadir / "download_manifest.tsv"
    with out.open("w") as fh:
        fh.write("cohort\tfile\trole\tsize_bytes\tsize_mb\twithin_2gb\taction\turl\n")
        for row in rows:
            fh.write(
                f"{row['cohort']}\t{row['file']}\t{row['role']}\t{row['size_bytes']}\t"
                f"{row['size_mb']}\t{row['within_2gb']}\t{row['action']}\t{row['url']}\n"
            )
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
