#!/usr/bin/env python3
"""Download public processed UMIs for GSE131907 + GSE205335.

Both gzipped matrices are <2 GB. The 3 GB GSE131907 log2TPM text is the
same cells in another normalization and is not downloaded. EGA FASTQ
(EGAD00001005054, EGAD00001008703) is controlled-access and is not used.
GSE207422 is out of scope.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

GSE131907_SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl"
GSE131907_MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix"
GSE205335_SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl"
GSE205335_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"

FILES = {
    "GSE131907": [
        f"{GSE131907_SUPP}/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        f"{GSE131907_SUPP}/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        f"{GSE131907_MATRIX}/GSE131907_series_matrix.txt.gz",
    ],
    "GSE205335": [
        f"{GSE205335_SUPP}/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        f"{GSE205335_SUPP}/GSE205335_Lung_IO_CellIdentity.txt.gz",
        GSE205335_SOFT,
    ],
}

SKIPPED = [
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz (>2 GB extra; same cells)",
    "EGAD00001005054",
    "EGAD00001008703",
    "GSE207422",
]


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    return dest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp"))
    args = ap.parse_args()
    for gse, urls in FILES.items():
        dest_dir = args.datadir / gse
        for url in urls:
            fetch(url, dest_dir / Path(url).name)
    print("skipped:", "; ".join(SKIPPED), flush=True)


if __name__ == "__main__":
    main()
