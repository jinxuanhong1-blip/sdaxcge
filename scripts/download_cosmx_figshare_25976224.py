#!/usr/bin/env python3
"""Download the public CosMx human NSCLC clustered object (figshare 25976224).

He et al. 2022 / the public reprocess deposited as
`cosmx_human_nsclc_clustered.h5ad` (765,771 cells, 960 genes, 8 sections, 5 donors).
Does not download or write any private cohort.
"""

from __future__ import annotations

import os
import sys
import urllib.request

URL = "https://ndownloader.figshare.com/files/46841842"
EXPECTED_MIN_BYTES = 2_000_000_000


def main() -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dest_dir = os.path.join(root, "data", "cosmx")
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, "cosmx_human_nsclc_clustered.h5ad")
    if os.path.exists(dest) and os.path.getsize(dest) >= EXPECTED_MIN_BYTES:
        print(f"already present: {dest} ({os.path.getsize(dest)} bytes)")
        return 0
    print(f"GET {URL}", flush=True)
    urllib.request.urlretrieve(URL, dest)
    size = os.path.getsize(dest)
    if size < EXPECTED_MIN_BYTES:
        raise SystemExit(f"download too small: {size} bytes at {dest}")
    print(f"wrote {dest} ({size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
