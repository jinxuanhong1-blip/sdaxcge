#!/usr/bin/env python3
"""Download the three public CXCL13+ trio scRNA matrices (GEO)."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "GSE148071": [
        (
            "GSE148071_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/GSE148071_RAW.tar",
        ),
        (
            "GSE148071_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/GSE148071_series_matrix.txt.gz",
        ),
    ],
    "GSE207422": [
        (
            "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        ),
        (
            "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        ),
    ],
    "GSE253013": [
        (
            "GSE253013_all_luad_garnett_temp.rds.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/suppl/GSE253013_all_luad_garnett_temp.rds.gz",
        ),
        (
            "GSE253013_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/matrix/GSE253013_series_matrix.txt.gz",
        ),
    ],
}


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "wget",
        "-c",
        "--tries=8",
        "--waitretry=8",
        "--timeout=60",
        "--progress=dot:giga",
        "-O",
        str(tmp),
        url,
    ]
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd)
    tmp.rename(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/cxcl13_trio_data"))
    p.add_argument("--cohorts", default="GSE148071,GSE207422,GSE253013")
    args = p.parse_args()
    want = {c.strip() for c in args.cohorts.split(",") if c.strip()}
    for cohort, items in FILES.items():
        if cohort not in want:
            continue
        for name, url in items:
            wget(url, args.out / cohort / name)


if __name__ == "__main__":
    main()
