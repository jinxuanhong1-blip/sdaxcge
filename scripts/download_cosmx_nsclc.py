#!/usr/bin/env python3
"""Download official CosMx NSCLC 960-plex flat files (He et al. 2022).

8 samples / 5 patients from NanoString public S3. Extracts only
exprMat, metadata, and fov_positions CSVs. Does not use private 8-KL data.
"""

from __future__ import annotations

import argparse
import os
import sys
import tarfile
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "cosmx_nsclc")

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]

BASE = "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed"
KEEP_SUFFIXES = (
    "exprMat_file.csv",
    "metadata_file.csv",
    "fov_positions_file.csv",
)


def sample_url(sample: str) -> str:
    return f"{BASE}/{sample}/{sample}+SMI+Flat+data.tar.gz"


def needed_csvs(sample_dir: str, sample: str) -> bool:
    if not os.path.isdir(sample_dir):
        return False
    found = {s: False for s in KEEP_SUFFIXES}
    for root, _, files in os.walk(sample_dir):
        for fn in files:
            for s in KEEP_SUFFIXES:
                if fn.endswith(s):
                    found[s] = True
    return all(found.values())


def download(url: str, dest: str, retries: int = 4) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    delay = 4
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            print(f"GET {url} -> {dest} (attempt {attempt})", flush=True)
            urllib.request.urlretrieve(url, dest)
            if os.path.getsize(dest) < 1_000_000:
                raise RuntimeError(f"download too small: {os.path.getsize(dest)} bytes")
            return
        except Exception as err:
            last_err = err
            print(f"download failed: {err}", flush=True)
            if attempt < retries:
                time.sleep(delay)
                delay *= 2
    raise RuntimeError(f"failed to download {url}: {last_err}")


def extract_flat(tar_path: str, dest: str) -> list[str]:
    os.makedirs(dest, exist_ok=True)
    kept = []
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar:
            if not member.isfile():
                continue
            base = os.path.basename(member.name)
            if not any(base.endswith(s) for s in KEEP_SUFFIXES):
                continue
            member.name = base
            tar.extract(member, dest, filter="data")
            kept.append(os.path.join(dest, base))
            print(f"  extracted {base}", flush=True)
    if len(kept) < 2:
        raise RuntimeError(f"expected exprMat+metadata in {tar_path}, got {kept}")
    return kept


def ensure_sample(sample: str, keep_tarball: bool = False) -> str:
    sample_dir = os.path.join(DATA, sample)
    if needed_csvs(sample_dir, sample):
        print(f"{sample}: already present", flush=True)
        return sample_dir
    os.makedirs(sample_dir, exist_ok=True)
    tar_path = os.path.join(DATA, f"{sample}.tar.gz")
    download(sample_url(sample), tar_path)
    extract_flat(tar_path, sample_dir)
    if not keep_tarball:
        os.remove(tar_path)
    if not needed_csvs(sample_dir, sample):
        raise RuntimeError(f"{sample}: missing required CSVs after extract")
    return sample_dir


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="*", default=SAMPLES)
    parser.add_argument("--keep-tarball", action="store_true")
    args = parser.parse_args()
    os.makedirs(DATA, exist_ok=True)
    for sample in args.samples:
        if sample not in SAMPLES:
            raise SystemExit(f"unknown sample {sample}; official set is {SAMPLES}")
        ensure_sample(sample, keep_tarball=args.keep_tarball)
    print("all requested official CosMx NSCLC flat files ready", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
