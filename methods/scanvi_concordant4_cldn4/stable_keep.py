#!/usr/bin/env python3
"""Deterministic 10% keep. Same rule in Python and the GSE205335 R extract.

A QC-pass cell is kept when crc32(key) % 1_000_000 < 100_000.
key is dataset-specific but always "{unit_id}|{barcode}".
"""
from __future__ import annotations

import sys
import zlib

THRESHOLD = 100_000
MODULUS = 1_000_000


def stable_keep(key: str) -> bool:
    h = zlib.crc32(key.encode()) & 0xFFFFFFFF
    return (h % MODULUS) < THRESHOLD


def main() -> None:
    for line in sys.stdin:
        key = line.rstrip("\n")
        if key:
            print(1 if stable_keep(key) else 0)


if __name__ == "__main__":
    main()
