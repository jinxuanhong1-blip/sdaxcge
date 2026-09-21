#!/usr/bin/env python3
"""Download the public He et al. 2022 CosMx NSCLC clustered object.

figshare 25976224, file cosmx_human_nsclc_clustered.h5ad (765,771 cells).
Does not fetch or write any private 8-KL matrices.
"""

from __future__ import annotations

import os
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
URL = "https://ndownloader.figshare.com/files/46841842"
MIN_BYTES = 2_000_000_000


def download(url: str, dest: str, retries: int = 4) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    delay = 4
    last = None
    for attempt in range(1, retries + 1):
        try:
            print(f"GET {url} -> {dest} (attempt {attempt})", flush=True)
            urllib.request.urlretrieve(url, dest)
            size = os.path.getsize(dest)
            if size < MIN_BYTES:
                raise RuntimeError(f"download too small: {size} bytes")
            print(f"wrote {size} bytes", flush=True)
            return
        except Exception as err:
            last = err
            print(f"download failed: {err}", flush=True)
            if attempt < retries:
                time.sleep(delay)
                delay *= 2
    raise RuntimeError(f"failed to download {url}: {last}")


def main() -> int:
    if os.path.exists(DEST) and os.path.getsize(DEST) >= MIN_BYTES:
        print(f"already present: {DEST}", flush=True)
        return 0
    download(URL, DEST)
    return 0


if __name__ == "__main__":
    sys.exit(main())
