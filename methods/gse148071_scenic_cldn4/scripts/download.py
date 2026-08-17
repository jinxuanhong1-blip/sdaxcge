#!/usr/bin/env python3
"""Download public GSE148071 processed UMI matrices + series matrix (GEO).

Public processed only. Raw FASTQ is not fetched. GSE207422 / GSE131907
are not used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

FILES = {
    "GSE148071_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "GSE148071_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/"
        "GSE148071_series_matrix.txt.gz"
    ),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_tar(tar_path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.rglob("*_exp.txt.gz")):
        return
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data",
    )
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {"accession": "GSE148071", "files": {}}
    for name, url in FILES.items():
        dest = args.out / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        else:
            print(f"GET {url}", flush=True)
            urllib.request.urlretrieve(url, dest)
            print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
        manifest["files"][name] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
        }
    extract_tar(args.out / "GSE148071_RAW.tar", args.out / "GSE148071_files")
    n = len(list((args.out / "GSE148071_files").rglob("*_exp.txt.gz")))
    manifest["n_exp_matrices"] = n
    (args.out / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(f"extracted exp matrices: {n}", flush=True)


if __name__ == "__main__":
    main()
