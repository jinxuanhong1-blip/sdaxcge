"""Download public inputs for the TCGA-UCEC B1 analog (TACSTD2 surfaceome ranking).

Files land in ``data/`` (git-ignored). Re-running is cheap: existing files are
left in place unless ``--force`` is given.

1. TCGA-UCEC STAR FPKM-UQ matrix from the UCSC Xena GDC hub
   (``TCGA-UCEC.star_fpkm-uq.tsv.gz``). Values are log2(fpkm-uq + 1),
   rows are versioned Ensembl gene IDs.
2. GENCODE v36 gene probeMap (Ensembl id → HGNC symbol) from the same hub.
3. Bausch-Fluck et al. 2018 in-silico surfaceome table S3. The ETH host now
   serves a single-page app for that path, so we pull the identical workbook
   vendored in the steveneschrich/surfaceome R package (same source as the
   BRCA B1 analog).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

DATA_DIR = Path(__file__).resolve().parents[3] / "data"

FILES = {
    "TCGA-UCEC.star_fpkm-uq.tsv.gz": (
        "https://gdc.xenahubs.net/download/TCGA-UCEC.star_fpkm-uq.tsv.gz"
    ),
    "gencode.v36.annotation.gtf.gene.probemap": (
        "https://gdc.xenahubs.net/download/gencode.v36.annotation.gtf.gene.probemap"
    ),
    "table_S3_surfaceome.xlsx": (
        "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
        "main/data-raw/surfy/table_S3_surfaceome.xlsx"
    ),
}


def _download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and dest.stat().st_size > 1000 and not force:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        tmp.replace(dest)
    print(f"[ok  ] {dest.name} ({dest.stat().st_size:,} bytes)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    args = ap.parse_args(argv)
    for name, url in FILES.items():
        _download(url, DATA_DIR / name, args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
