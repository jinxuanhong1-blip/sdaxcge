#!/usr/bin/env python3
"""Download public processed matrices for GSE123902 + GSE189357.

GSE123902_RAW.tar ~90 MB (Laughney 2020 dense CSV).
GSE189357_RAW.tar ~624 MB (Zhu/Wang AIS–IAC 10x MTX).
No FASTQ. Public GEO only.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
import urllib.request
from pathlib import Path

FILES = {
    "GSE123902_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar",
        80_000_000,
    ),
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
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_pair_123902_189357"))
    args = p.parse_args()
    for name, (url, min_bytes) in FILES.items():
        fetch(url, args.out / name, min_bytes)
    print("OK", args.out, file=sys.stderr)
    # keep a copy listing
    listing = args.out / "DOWNLOAD.txt"
    listing.write_text(
        "\n".join(f"{p.name}\t{p.stat().st_size}" for p in sorted(args.out.glob('*.tar'))) + "\n"
    )
    if shutil.which("tar"):
        print("tars ready; analyze.py reads them in place", flush=True)


if __name__ == "__main__":
    main()
