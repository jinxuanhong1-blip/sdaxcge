#!/usr/bin/env python3
"""Download public processed matrices for GSE123902 + GSE131907.

Skip the 36.5 GB Laughney annotated H5, the 2.9 GB GSE131907 log2TPM text,
EGA/FASTQ, GSE148071, and GSE207422.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

FILES = {
    "gse123902/GSE123902_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
    ),
    "gse123902/GSE123902_GEO_README.rtf": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_GEO_README.rtf"
    ),
    "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "gse131907/GSE131907_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
        "GSE131907_series_matrix.txt.gz"
    ),
    "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
}

SKIPPED = {
    "PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5": {
        "bytes": 36489164157,
        "reason": "author annotated H5 is 36.5 GB (>2 GB public-processed gate)",
    },
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz": {
        "reason": "2.9 GB log2TPM text; UMI matrix exists and is <2 GB",
    },
    "GSE148071": {"reason": "explicitly not added (dilutes the PR #459 pair)"},
    "GSE207422": {"reason": "explicitly not added"},
}


def curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = ["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(tmp), url]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix == ".partial":
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
                "accessions": ["GSE123902", "GSE131907"],
                "skipped": SKIPPED,
                "files": rows,
            },
            indent=2,
        )
    )
    print(f"manifest {out}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_123902_131907"))
    args = p.parse_args()
    for rel, url in FILES.items():
        curl_download(url, args.out / rel)
    print("skipped:")
    for name, rec in SKIPPED.items():
        print(f"  {name}: {rec['reason']}")
    manifest(args.out)


if __name__ == "__main__":
    main()
