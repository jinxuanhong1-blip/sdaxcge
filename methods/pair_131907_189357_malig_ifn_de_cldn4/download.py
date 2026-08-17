#!/usr/bin/env python3
"""Download public processed GSE189357 MTX (Zhu/Wang AIS–IAC).

GSE189357_RAW.tar ~624 MB. No FASTQ. Public GEO only. GSE131907 counts are
the committed author-malignant UMI-sum from PR #456 / #472 (not re-downloaded).
"""
from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path

FILES = {
    "GSE189357_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar",
        500_000_000,
    ),
}


def fetch(url: str, dest: Path, min_bytes: int, retries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    delay = 4
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        tmp = dest.with_suffix(dest.suffix + ".partial")
        try:
            print(f"GET {url} -> {dest} (try {attempt}/{retries})", flush=True)
            if tmp.exists():
                tmp.unlink()
            urllib.request.urlretrieve(url, tmp)
            size = tmp.stat().st_size
            if size < min_bytes:
                raise RuntimeError(f"too small: {size} < {min_bytes}")
            tmp.replace(dest)
            print(f"wrote {dest} ({size} bytes)", flush=True)
            return
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            print(f"FAIL {url}: {exc}", flush=True)
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                time.sleep(delay)
                delay *= 2
    raise SystemExit(f"could not download {url}: {last_err}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_pair_131907_189357"))
    args = p.parse_args()
    for name, (url, min_bytes) in FILES.items():
        fetch(url, args.out / name, min_bytes)
    print("OK", args.out)


if __name__ == "__main__":
    main()
