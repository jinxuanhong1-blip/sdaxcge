#!/usr/bin/env python3
"""Download Mirhadi et al. Nat Commun 2022 Supplementary Data 1.

Public NSCLC PDX TMT proteome (PXD016579). The xlsx is the normalized
log2 protein matrix plus clinical histology. Raw PRIDE files are not used.
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-022-29444-9/MediaObjects/"
    "41467_2022_29444_MOESM4_ESM.xlsx"
)
FILENAME = "41467_2022_29444_MOESM4_ESM.xlsx"


def download(outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / FILENAME
    if dest.exists() and dest.stat().st_size > 1_000_000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)")
        return dest
    print(f"GET {URL}")
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        data = resp.read()
    dest.write_bytes(data)
    print(f"wrote {dest} ({len(data)} bytes)")
    return dest


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type=Path, default=Path("data/pdx_cldn4_dnapk_sting"))
    args = p.parse_args()
    download(args.outdir)


if __name__ == "__main__":
    main()
