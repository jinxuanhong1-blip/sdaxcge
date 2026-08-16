#!/usr/bin/env python3
"""Download public GSE207422 processed files and NicheNet-v2 priors.

GSE253013 is not fetched: 9.3 GB RDS, no MPR/NMPR labels (see playbook).
"""
from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/scrna_nichenet")
DATA.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

FILES = [
    {
        "name": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
        ),
        "dest": CACHE / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "min_bytes": 100_000_000,
    },
    {
        "name": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
        ),
        "dest": DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "min_bytes": 5_000,
    },
    {
        "name": "lr_network_human_21122021.rds",
        "url": (
            "https://zenodo.org/records/7074291/files/"
            "lr_network_human_21122021.rds?download=1"
        ),
        "dest": CACHE / "lr_network_human_21122021.rds",
        "min_bytes": 10_000,
    },
    {
        "name": "ligand_target_matrix_nsga2r_final.rds",
        "url": (
            "https://zenodo.org/records/7074291/files/"
            "ligand_target_matrix_nsga2r_final.rds?download=1"
        ),
        "dest": CACHE / "ligand_target_matrix_nsga2r_final.rds",
        "min_bytes": 200_000_000,
    },
]


def sha256_file(path: Path, nbytes: int = 1_000_000) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()[:16]


def download(url: str, dest: Path, min_bytes: int, attempts: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    wait = 4
    last_err = None
    for i in range(1, attempts + 1):
        try:
            print(f"GET {url} -> {dest} (try {i}/{attempts})", flush=True)
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                written = 0
                with tmp.open("wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        if chunk:
                            fh.write(chunk)
                            written += len(chunk)
            if written < min_bytes:
                raise RuntimeError(f"too small: {written} < {min_bytes}")
            tmp.replace(dest)
            print(
                f"wrote {dest} {dest.stat().st_size} bytes sha256_1MB={sha256_file(dest)}",
                flush=True,
            )
            return
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            print(f"fail try {i}: {exc}", flush=True)
            if tmp.exists():
                tmp.unlink()
            if i < attempts:
                time.sleep(wait)
                wait *= 2
    raise RuntimeError(f"download failed {url}: {last_err}")


def main() -> int:
    for spec in FILES:
        download(spec["url"], spec["dest"], spec["min_bytes"])
    print("ok", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
