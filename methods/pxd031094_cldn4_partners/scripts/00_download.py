#!/usr/bin/env python3
"""Download the public PXD031094 MaxQuant SEARCH table for CLDN4 CoIP.

Skips RAW spectra. Source: PRIDE FTP HTTPS mirror.
"""
from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

URL = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD031094/"
    "proteinGroups_Cldn4.txt"
)
# PRIDE Archive SHA-1 for proteinGroups_Cldn4.txt (SEARCH)
PRIDE_SHA1 = "4dbde6c0d827105412ef1f34cbb969b8dd8acb36"
EXPECTED_BYTES = 3950383

OUT = Path(__file__).resolve().parents[1] / "data" / "proteinGroups_Cldn4.txt"


def sha1(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists() and OUT.stat().st_size == EXPECTED_BYTES and sha1(OUT) == PRIDE_SHA1:
        print(f"already present {OUT} ({OUT.stat().st_size} bytes) sha1={PRIDE_SHA1}")
        return 0
    print(f"GET {URL}")
    urllib.request.urlretrieve(URL, OUT)
    digest = sha1(OUT)
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes) sha1={digest}")
    if digest != PRIDE_SHA1:
        print(f"WARNING: SHA-1 {digest} != PRIDE {PRIDE_SHA1}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
