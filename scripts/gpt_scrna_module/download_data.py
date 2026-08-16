#!/usr/bin/env python3
"""Download the open processed inputs used by analyze_atlases.py."""

from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path


FILES = {
    "GSE131907_cell_annotation.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        1_886_187,
        "9d017652a31a489c1f62cacc436936d48a5b07bb554fc0dbe602da7474532205",
    ),
    "GSE131907_raw_UMI_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        408_736_818,
        "b44a9c6635d0b59ec66b5eb6b40b40e7292c3a75aaf8afc4f933eb0550a26372",
    ),
    "GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar",
        180_193_280,
        "6bf70027c88bc0659c25dc8e5251f70351fce6b90913bf4ff33ca821df6e906f",
    ),
    "GSE154826_cell_metadata.csv": (
        "https://raw.githubusercontent.com/effiken/Leader_et_al/master/"
        "input_tables/cell_metadata.csv",
        13_037_071,
        "e7000734f1da8c8e7c5027c54ac1c7d8222179f84b2f575014af2a4cc6629add",
    ),
    "GSE154826_annots_list.csv": (
        "https://raw.githubusercontent.com/effiken/Leader_et_al/master/"
        "input_tables/annots_list.csv",
        1_818,
        "e8f28b7deba93d1a791d768a8e11ba7c4ed36c71eb960a51718f2656cffc99bd",
    ),
    "GSE154826_sample_metadata.csv": (
        "https://raw.githubusercontent.com/effiken/Leader_et_al/master/"
        "input_tables/table_s1_sample_table.csv",
        8_893,
        "e4c78eab1a986bd966972ab92dffbc66e1ddb4409e241f8212fe715861ddf914",
    ),
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    for name, (url, expected_bytes, expected_sha) in FILES.items():
        destination = args.output / name
        if not destination.exists():
            print(f"Downloading {name}", flush=True)
            urllib.request.urlretrieve(url, destination)
        observed = (destination.stat().st_size, digest(destination))
        expected = (expected_bytes, expected_sha)
        if observed != expected:
            raise RuntimeError(f"Integrity failure for {name}: {observed} != {expected}")
        print(f"Verified {name}: {observed[0]} bytes")

    annotation_gz = args.output / "GSE131907_cell_annotation.txt.gz"
    annotation_txt = args.output / "GSE131907_cell_annotation.txt"
    if not annotation_txt.exists():
        import gzip
        import shutil

        with gzip.open(annotation_gz, "rb") as source, annotation_txt.open("wb") as target:
            shutil.copyfileobj(source, target)


if __name__ == "__main__":
    main()
