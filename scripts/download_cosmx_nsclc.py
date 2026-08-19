#!/usr/bin/env python3
"""Download official CosMx NSCLC 960-plex flat files (He et al. 2022).

8 samples / 5 patients from NanoString public S3. Streams each tarball and
extracts only exprMat, metadata, and fov_positions CSVs (no images).
Does not use private 8-KL data.
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


def needed_csvs(sample_dir: str) -> bool:
    if not os.path.isdir(sample_dir):
        return False
    found = {s: False for s in KEEP_SUFFIXES}
    for _root, _dirs, files in os.walk(sample_dir):
        for fn in files:
            for s in KEEP_SUFFIXES:
                if fn.endswith(s):
                    found[s] = True
    return all(found.values())


def stream_extract(url: str, dest: str, retries: int = 4) -> list[str]:
    os.makedirs(dest, exist_ok=True)
    delay = 4
    last_err = None
    for attempt in range(1, retries + 1):
        kept: list[str] = []
        try:
            print(f"STREAM {url} (attempt {attempt})", flush=True)
            req = urllib.request.Request(url, headers={"User-Agent": "cosmx-nsclc-download/1.0"})
            with urllib.request.urlopen(req, timeout=600) as resp:
                with tarfile.open(fileobj=resp, mode="r|gz") as tar:
                    for member in tar:
                        if not member.isfile():
                            continue
                        base = os.path.basename(member.name)
                        if not any(base.endswith(s) for s in KEEP_SUFFIXES):
                            continue
                        dest_path = os.path.join(dest, base)
                        src = tar.extractfile(member)
                        if src is None:
                            continue
                        with open(dest_path, "wb") as fh:
                            while True:
                                chunk = src.read(1024 * 1024)
                                if not chunk:
                                    break
                                fh.write(chunk)
                        kept.append(dest_path)
                        print(f"  extracted {base} ({os.path.getsize(dest_path)} bytes)", flush=True)
            if len(kept) < 2:
                raise RuntimeError(f"expected exprMat+metadata, got {kept}")
            return kept
        except Exception as err:
            last_err = err
            print(f"download failed: {err}", flush=True)
            if attempt < retries:
                time.sleep(delay)
                delay *= 2
    raise RuntimeError(f"failed to download {url}: {last_err}")


def ensure_sample(sample: str) -> str:
    sample_dir = os.path.join(DATA, sample)
    if needed_csvs(sample_dir):
        print(f"{sample}: already present", flush=True)
        return sample_dir
    os.makedirs(sample_dir, exist_ok=True)
    stream_extract(sample_url(sample), sample_dir)
    if not needed_csvs(sample_dir):
        raise RuntimeError(f"{sample}: missing required CSVs after extract")
    return sample_dir


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="*", default=SAMPLES)
    args = parser.parse_args()
    os.makedirs(DATA, exist_ok=True)
    for sample in args.samples:
        if sample not in SAMPLES:
            raise SystemExit(f"unknown sample {sample}; official set is {SAMPLES}")
        ensure_sample(sample)
    print("all requested official CosMx NSCLC flat files ready", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
