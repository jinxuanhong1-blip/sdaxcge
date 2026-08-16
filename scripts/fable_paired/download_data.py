"""Reproducibly download every raw GEO file used by this analysis.

Files land in RAW_DIR (default /tmp/fable_raw), OUTSIDE the repo, so the
committed processed artifacts remain small. Re-running is idempotent: existing,
non-empty files are skipped unless FABLE_FORCE_DOWNLOAD=1.
"""
from __future__ import annotations

import os
import sys
import time
import urllib.request

from config import GEO_BASE, GEO_FILES, raw_path


def _download(url: str, dest, retries: int = 4) -> None:
    delay = 4
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=600) as resp, open(dest, "wb") as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            return
        except Exception as exc:  # network hiccups -> exponential backoff
            if attempt == retries:
                raise
            print(f"  retry {attempt} after error: {exc}", file=sys.stderr)
            time.sleep(delay)
            delay *= 2


def main() -> None:
    force = os.environ.get("FABLE_FORCE_DOWNLOAD") == "1"
    for key, (rel, fname) in GEO_FILES.items():
        dest = raw_path(key)
        if dest.exists() and dest.stat().st_size > 0 and not force:
            print(f"[skip] {key}: {dest} ({dest.stat().st_size/1e6:.1f} MB)")
            continue
        url = f"{GEO_BASE}/{rel}/{fname}"
        print(f"[get ] {key}: {url}")
        _download(url, dest)
        print(f"       -> {dest} ({dest.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
