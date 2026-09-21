#!/usr/bin/env python3
"""Download official CosMx NSCLC 960-plex flat files (He et al. 2022).

Eight sections, five tissues, from the NanoString public S3 release.
Keeps exprMat, metadata, and fov position tables. Does not use private data.
"""

from __future__ import annotations

import argparse
import os
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
KEEP = ("exprMat_file.csv", "metadata_file.csv", "fov_positions_file.csv")


def ready(sample_dir: str) -> bool:
    if not os.path.isdir(sample_dir):
        return False
    found = {s: False for s in KEEP}
    for fn in os.listdir(sample_dir):
        for s in KEEP:
            if fn.endswith(s):
                found[s] = True
    return all(found.values())


def download(url: str, dest: str) -> None:
    delay = 4
    last = None
    for attempt in range(1, 5):
        try:
            print(f"GET {url}", flush=True)
            urllib.request.urlretrieve(url, dest)
            if os.path.getsize(dest) < 1_000_000:
                raise RuntimeError(f"download too small: {os.path.getsize(dest)}")
            return
        except Exception as err:
            last = err
            print(f"download failed: {err}", flush=True)
            if attempt < 4:
                time.sleep(delay)
                delay *= 2
    raise RuntimeError(f"failed to download {url}: {last}")


def extract(tar_path: str, dest: str) -> None:
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar:
            if not member.isfile():
                continue
            base = os.path.basename(member.name)
            if not any(base.endswith(s) for s in KEEP):
                continue
            member.name = base
            tar.extract(member, dest, filter="data")
            print(f"  extracted {base}", flush=True)


def ensure(sample: str) -> None:
    sample_dir = os.path.join(DATA, sample)
    if ready(sample_dir):
        print(f"{sample}: already present", flush=True)
        return
    os.makedirs(sample_dir, exist_ok=True)
    tar_path = os.path.join(DATA, f"{sample}.tar.gz")
    url = f"{BASE}/{sample}/{sample}+SMI+Flat+data.tar.gz"
    download(url, tar_path)
    extract(tar_path, sample_dir)
    os.remove(tar_path)
    if not ready(sample_dir):
        raise RuntimeError(f"{sample}: missing csv after extract")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="*", default=SAMPLES)
    args = parser.parse_args()
    for sample in args.samples:
        if sample not in SAMPLES:
            raise SystemExit(f"unknown sample {sample}")
        ensure(sample)
    print("official CosMx NSCLC flat files ready", flush=True)


if __name__ == "__main__":
    main()
