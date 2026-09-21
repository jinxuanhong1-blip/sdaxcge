#!/usr/bin/env python3
"""Download the clustered He et al. 2022 CosMx NSCLC object.

figshare 25976224. Eight sections, five patients, 765,771 cells.
The h5ad is gitignored (about 2.8 GB).
"""

from __future__ import annotations

import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
URL = "https://ndownloader.figshare.com/files/46841842"
EXPECTED_BYTES = 2_755_776_882


def main() -> int:
    if os.path.isfile(DEST) and os.path.getsize(DEST) == EXPECTED_BYTES:
        print(f"already present: {DEST}")
        return 0
    os.makedirs(os.path.dirname(DEST), exist_ok=True)
    print(f"GET {URL}", flush=True)
    urllib.request.urlretrieve(URL, DEST)
    got = os.path.getsize(DEST)
    if got != EXPECTED_BYTES:
        raise SystemExit(f"size {got} != expected {EXPECTED_BYTES}")
    print(f"wrote {DEST} ({got} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
