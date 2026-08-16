#!/usr/bin/env python3
"""Download GSE131907 processed GEO supplements used by this hunt."""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size} bytes)")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[get] {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    print(f"[ok]  {dest.name} ({dest.stat().st_size} bytes)")


def main() -> None:
    C.DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in C.DOWNLOADS.items():
        fetch(url, C.DATA_DIR / name)


if __name__ == "__main__":
    main()
