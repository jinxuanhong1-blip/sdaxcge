#!/usr/bin/env python3
"""Download the public GSE253013 processed object (GEO supplementary RDS).

GEO series matrix is metadata-only (6.7 KB; no expression). The only public
processed matrix is GSE253013_all_luad_garnett_temp.rds.gz (9.3 GB). That file
is double-gzipped XDR RDS (R 3.5.2, format v2).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/suppl/"
    "GSE253013_all_luad_garnett_temp.rds.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/matrix/"
    "GSE253013_series_matrix.txt.gz"
)
EXPECTED_RDS_BYTES = 9966268875


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["wget", "-c", "-O", str(dest), url]
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("data/gse253013"),
        help="Directory for GEO downloads (not committed).",
    )
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    rds = out / "GSE253013_all_luad_garnett_temp.rds.gz"
    matrix = out / "GSE253013_series_matrix.txt.gz"
    wget(MATRIX_URL, matrix)
    wget(GEO_URL, rds)

    size = rds.stat().st_size
    sha = hashlib.sha256()
    with rds.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            sha.update(chunk)

    manifest = {
        "dataset": "GSE253013",
        "pmid": "38335304",
        "title": (
            "Single-cell RNA sequencing identifies dysregulated NOTCH3 "
            "signaling in tumor stroma of lung adenocarcinoma"
        ),
        "rds_url": GEO_URL,
        "rds_path": str(rds),
        "rds_bytes": size,
        "rds_bytes_expected": EXPECTED_RDS_BYTES,
        "rds_sha256": sha.hexdigest(),
        "series_matrix_url": MATRIX_URL,
        "series_matrix_bytes": matrix.stat().st_size,
        "series_matrix_has_expression": False,
        "note": (
            "89 GSM 10x lanes from 9 treatment-naive LUAD patients "
            "(tumor + adjacent non-tumor lung). Not an ICI/neoadjuvant "
            "response cohort. No public MPR/R labels in the series matrix."
        ),
    }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    if size != EXPECTED_RDS_BYTES:
        print(f"WARNING: size {size} != expected {EXPECTED_RDS_BYTES}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
