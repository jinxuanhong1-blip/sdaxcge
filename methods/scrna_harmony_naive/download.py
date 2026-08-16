#!/usr/bin/env python3
"""Download public GEO supplements used by the naive Harmony joint object."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from config import DATA, FILES, URLS


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"have {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return
    cmd = [
        "curl",
        "-L",
        "--retry",
        "5",
        "--retry-delay",
        "8",
        "--retry-all-errors",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    print("GET", url, "->", dest, flush=True)
    subprocess.check_call(cmd)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    for key, url in URLS.items():
        fetch(url, FILES[key])
    print("downloads complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
