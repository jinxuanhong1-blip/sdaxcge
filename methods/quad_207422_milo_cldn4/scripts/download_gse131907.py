#!/usr/bin/env python3
"""Download public GSE131907 processed files used by this folder.

UMI text + annotation + series matrix only. The 3 GB log2TPM text is the
same cells in another normalization and is not needed for Milo. EGA FASTQ
is controlled-access and is not used. GSE207422 is a different agent.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl"
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix"

FILES = [
    f"{SUPP}/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    f"{SUPP}/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    f"{MATRIX}/GSE131907_series_matrix.txt.gz",
]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE131907"))
    args = ap.parse_args()
    for url in FILES:
        fetch(url, args.datadir / Path(url).name)


if __name__ == "__main__":
    main()
