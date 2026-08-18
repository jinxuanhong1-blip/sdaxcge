#!/usr/bin/env python3
"""Download public processed UMIs for GSE131907 + GSE189357.

Skip the 2.86 GB GSE131907 log2TPM text, EGA/FASTQ, and GSE148071.
Refuse anything ≥ 2 GB.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

FILES = {
    "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "gse189357/GSE189357_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/"
        "GSE189357_RAW.tar"
    ),
    "gse189357/GSE189357_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/matrix/"
        "GSE189357_series_matrix.txt.gz"
    ),
}

MAX_BYTES = 2_000_000_000


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
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
    subprocess.check_call(cmd)
    size = tmp.stat().st_size
    if size >= MAX_BYTES:
        tmp.unlink()
        raise SystemExit(f"refusing {dest.name}: {size} bytes ≥ 2 GB")
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
                "accessions": ["GSE131907", "GSE189357"],
                "skipped": [
                    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
                    "GSE148071",
                    "EGA FASTQ",
                    "SRA FASTQ",
                ],
                "files": rows,
            },
            indent=2,
        )
    )
    print(f"manifest {out}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_131907_189357"))
    args = p.parse_args()
    for name, url in FILES.items():
        fetch(url, args.out / name)
    manifest(args.out)
    print("skipped GSE148071 (do not add).", flush=True)


if __name__ == "__main__":
    main()
