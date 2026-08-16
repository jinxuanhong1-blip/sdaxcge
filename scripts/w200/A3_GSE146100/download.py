#!/usr/bin/env python3
"""Download the public GSE146100 normalized matrix (not committed)."""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146100/"
    "suppl/GSE146100_NormData.txt.gz"
)
OUT = Path("data/GSE146100_NormData.txt.gz")
# Size-only check; GEO does not publish a stable checksum for this file.
EXPECTED_BYTES = 229_465_901


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists() and OUT.stat().st_size == EXPECTED_BYTES:
        print(f"already present: {OUT} ({OUT.stat().st_size} bytes)")
        print(f"sha256={sha256(OUT)}")
        return
    print(f"downloading {URL}")
    urllib.request.urlretrieve(URL, OUT)
    size = OUT.stat().st_size
    print(f"wrote {OUT} ({size} bytes) sha256={sha256(OUT)}")
    if size != EXPECTED_BYTES:
        raise SystemExit(f"unexpected size {size} != {EXPECTED_BYTES}")


if __name__ == "__main__":
    main()
