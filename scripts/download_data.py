#!/usr/bin/env python3
"""Download inputs for the TCGA-COAD TACSTD2 surfaceome co-expression analysis.

Two inputs land in ``data/`` (git-ignored; the expression matrix is ~110 MB):

1. TCGA-COAD RNA-seq from the NCI GDC open-access API
   (workflow ``STAR - Counts``, GENCODE v36, ``tpm_unstranded``).
   Built by ``src/fetch_gdc_coad.py`` into ``data/coad_star_tpm.parquet``
   plus a sample manifest and gene annotation.

2. The "in silico human surfaceome" master table (Bausch-Fluck et al.,
   PNAS 2018, table S3). The Wollscheid-lab host now serves a single-page
   app for the original path, so we pull the identical workbook vendored
   in the steveneschrich/surfaceome R package on GitHub — the same file
   used by the BRCA B1 analog.

Re-running is cheap: existing files are left in place unless --force is given.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"

SURFACEOME_URL = (
    "https://raw.githubusercontent.com/steveneschrich/surfaceome/"
    "main/data-raw/surfy/table_S3_surfaceome.xlsx"
)
SURFACEOME_FILE = DATA_DIR / "table_S3_surfaceome.xlsx"

# Independent robustness matrix (same Xena HiSeqV2 build used by B1_BRCA).
XENA_URL = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.COAD.sampleMap%2FHiSeqV2.gz"
)
XENA_FILE = DATA_DIR / "TCGA-COAD.HiSeqV2.gz"

EXPR_FILES = (
    DATA_DIR / "coad_star_tpm.parquet",
    DATA_DIR / "coad_sample_manifest.tsv",
    DATA_DIR / "coad_gene_annotation.tsv",
)


def _download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    with requests.get(url, stream=True, timeout=180) as r:
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
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args(argv)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    missing = [p for p in EXPR_FILES if not p.exists()]
    if missing or args.force:
        cmd = [
            sys.executable,
            str(ROOT / "src" / "fetch_gdc_coad.py"),
            "--outdir",
            str(DATA_DIR),
            "--workers",
            str(args.workers),
        ]
        print(f"[run ] {' '.join(cmd)}")
        subprocess.check_call(cmd)
    else:
        print("[skip] GDC COAD STAR-Counts matrix already present")

    _download(SURFACEOME_URL, SURFACEOME_FILE, args.force)
    _download(XENA_URL, XENA_FILE, args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
