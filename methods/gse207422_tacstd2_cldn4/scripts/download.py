#!/usr/bin/env python3
"""Re-download the public GSE207422 processed UMI matrix and sample metadata.

Hu et al., Genome Medicine 2023 (PMID 36869384). GEO deposits the author UMI
matrix and a sample-level clinical table. No barcode-level cell-type file.
"""

from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
FILES = {
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz": 184_001_817,
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx": None,
}
UA = "Mozilla/5.0 (research; GSE207422 TACSTD2/CLDN4 reanalysis)"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("data/GSE207422"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, expect in FILES.items():
        dest = args.outdir / name
        if dest.exists() and dest.stat().st_size > 0 and (expect is None or dest.stat().st_size == expect):
            print(f"exists {dest} ({dest.stat().st_size} bytes) sha256={sha256(dest)}")
            continue
        url = f"{BASE}/{name}"
        print(f"downloading {url}")
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as response, dest.open("wb") as handle:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                handle.write(chunk)
        print(f"wrote {dest} ({dest.stat().st_size} bytes) sha256={sha256(dest)}")


if __name__ == "__main__":
    main()
