#!/usr/bin/env python3
"""Download the public He 2022 CosMx NSCLC object used by the IFN/STAT1 script.

figshare 25976224, file cosmx_human_nsclc_clustered.h5ad (CellCharter release).
Does not download private 8-KL data.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "cosmx" / "cosmx_human_nsclc_clustered.h5ad"
URL = "https://ndownloader.figshare.com/files/46841842"


def main() -> None:
    DEST.parent.mkdir(parents=True, exist_ok=True)
    if DEST.exists() and DEST.stat().st_size > 1_000_000_000:
        print(f"already present: {DEST}")
        return
    print(f"GET {URL}")
    urllib.request.urlretrieve(URL, DEST)
    print(f"wrote {DEST} ({DEST.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
