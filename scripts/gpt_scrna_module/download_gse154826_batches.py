#!/usr/bin/env python3
"""Download GSE154826 per-sample GEO MTX archives used for tumor samples."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import pandas as pd

BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/"
    "GSE154826_amp_batch_ID_{batch}.tar.gz"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.source / "GSE154826_sample_metadata.csv")
    tumor = meta[(meta["tissue"] == "Tumor") & (meta["Use.in.Clustering.Model."] == "Yes")]
    batches = sorted({int(x) for x in tumor["amp_batch_ID"]})
    for batch in batches:
        name = f"GSE154826_amp_batch_ID_{batch}.tar.gz"
        destination = args.output / name
        if destination.exists() and destination.stat().st_size > 0:
            print(f"exists {name} {destination.stat().st_size}", flush=True)
            continue
        url = BASE.format(batch=batch)
        print(f"Downloading {name}", flush=True)
        urllib.request.urlretrieve(url, destination)
        print(f"saved {name} {destination.stat().st_size}", flush=True)


if __name__ == "__main__":
    main()
