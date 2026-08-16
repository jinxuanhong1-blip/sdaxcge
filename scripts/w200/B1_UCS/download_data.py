#!/usr/bin/env python3
"""Download open-access inputs for the TCGA-UCS B1 surface-rank analog.

All files are fully open (no dbGaP / controlled access):

  * UCSC Xena GDC hub: STAR TPM (log2(TPM+1), GENCODE v36), clinical,
    and the GENCODE v36 gene probemap.
  * UCSC Xena TCGA hub: legacy HiSeqV2 symbol-indexed matrix (robustness).
  * Bausch-Fluck 2018 in-silico surfaceome table S3, from the
    steveneschrich/surfaceome GitHub mirror (the ETH host now serves a
    single-page app instead of the xlsx).
  * ABSOLUTE tumor purity from the open PanCanAtlas supplement via GDC.

Raw matrices land in DATA_DIR (default ``data/``, override with B1_UCS_DATA)
and are git-ignored. Re-running is cheap: existing non-empty files are kept
unless --force is given.

Usage:  python3 scripts/w200/B1_UCS/download_data.py
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("B1_UCS_DATA", ROOT / "data"))

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
TCGA_HUB = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"

FILES = {
    "TCGA-UCS.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-UCS.star_tpm.tsv.gz",
    "TCGA-UCS.clinical.tsv.gz": f"{GDC_HUB}/TCGA-UCS.clinical.tsv.gz",
    "gencode.v36.annotation.gtf.gene.probemap": (
        f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap"
    ),
    "TCGA-UCS.HiSeqV2.gz": f"{TCGA_HUB}/TCGA.UCS.sampleMap%2FHiSeqV2.gz",
    "table_S3_surfaceome.xlsx": (
        "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
        "main/data-raw/surfy/table_S3_surfaceome.xlsx"
    ),
    "tcga_absolute_purity.txt": (
        "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
    ),
}


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str, url: str, force: bool) -> None:
    dest = DATA_DIR / name
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {name} already present ({dest.stat().st_size:,} bytes)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    tmp.replace(dest)
    print(f"[ok  ] {name}  size={dest.stat().st_size:,}  md5={_md5(dest)}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(name, url, args.force)
    print(f"\nAll files in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
