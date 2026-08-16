#!/usr/bin/env python3
"""Download GSE245459 processed FPKM (GEO supplementary, ~14 MB)."""

from __future__ import annotations

import urllib.request
from pathlib import Path

URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245459/"
    "suppl/GSE245459_fpkm.anno.txt.gz"
)
OUT = (
    Path(__file__).resolve().parents[3]
    / "results"
    / "w200"
    / "C4_GSE245459"
    / "raw"
    / "GSE245459_fpkm.anno.txt.gz"
)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists() and OUT.stat().st_size > 1_000_000:
        print(f"already present: {OUT} ({OUT.stat().st_size} bytes)")
        return
    print(f"GET {URL}")
    urllib.request.urlretrieve(URL, OUT)
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
