#!/usr/bin/env python3
"""Download public processed matrices for the concordant four.

GSE123902 + GSE131907 + GSE205335 + GSE189357 only.
Do NOT add GSE148071 / GSE127465 / CD45-only / 7-pool.
No FASTQ / EGA / SRA. Refuse anything ≥ 2 GB.
Skip the 2.86 GB GSE131907 log2TPM text and the 36.5 GB Laughney H5.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

FILES = {
    "GSE123902": [
        (
            "GSE123902_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar",
        ),
    ],
    "GSE131907": [
        (
            "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        ),
        (
            "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        ),
    ],
    "GSE205335": [
        (
            "GSE205335_Lung_IO_UMI_matrix.rds.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        ),
        (
            "GSE205335_Lung_IO_CellIdentity.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        ),
        (
            "GSE205335_family.soft.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
        ),
    ],
    "GSE189357": [
        (
            "GSE189357_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar",
        ),
    ],
}

MAX_BYTES = 2_000_000_000
REFUSED = (
    "GSE148071",
    "GSE127465",
    "GSE229353",
    "CD45-only",
    "7-pool",
)


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
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
    print(f"OK {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/concordant4_raw"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    prov = []
    for cohort, items in FILES.items():
        d = args.outdir / cohort
        d.mkdir(parents=True, exist_ok=True)
        for name, url in items:
            dest = d / name
            print(f"GET {cohort} {name}", flush=True)
            wget(url, dest)
            prov.append(
                {
                    "cohort": cohort,
                    "filename": name,
                    "url": url,
                    "bytes": dest.stat().st_size,
                }
            )
    (args.outdir / "provenance_downloads.json").write_text(json.dumps(prov, indent=2) + "\n")
    print("wrote", args.outdir / "provenance_downloads.json", flush=True)
    print("refused:", ", ".join(REFUSED), flush=True)


if __name__ == "__main__":
    main()
