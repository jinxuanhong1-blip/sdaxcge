#!/usr/bin/env python3
"""Download public TCGA-LUAD HiSeqV2 + MDACC ESTIMATE for the CLDN4–CXCL leftover.

Does not re-audit Q4 vs CD8 (PR #352) or TJ7 vs CD8 (PR #69 / #90).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "data" / "tcga_luad_cldn4_cxcl"

HISEQ_URLS = [
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
]
EST_URLS = [
    "https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
    "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
]


def download(urls: list[str], dest: Path, min_bytes: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"cached {dest.name} ({dest.stat().st_size} bytes)")
        return
    last = None
    for url in urls:
        for attempt in range(5):
            try:
                req = Request(url, headers={"User-Agent": "tcga-luad-cldn4-cxcl/1.0"})
                with urlopen(req, timeout=300) as r:
                    data = r.read()
                if len(data) < min_bytes:
                    raise RuntimeError(f"too small: {len(data)} bytes from {url}")
                tmp = dest.with_suffix(dest.suffix + ".tmp")
                tmp.write_bytes(data)
                tmp.replace(dest)
                print(f"downloaded {dest.name} ({len(data)} bytes) <- {url}")
                return
            except Exception as e:  # noqa: BLE001
                last = e
                print(f"retry {attempt + 1} {url}: {e}")
                time.sleep(2**attempt)
    raise SystemExit(f"failed to download {dest}: {last}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=DATA)
    args = ap.parse_args()
    download(HISEQ_URLS, args.datadir / "TCGA.LUAD.HiSeqV2.gz", min_bytes=1_000_000)
    download(EST_URLS, args.datadir / "MDACC_estimate_LUAD_RNAseqV2.txt", min_bytes=200)


if __name__ == "__main__":
    main()
