#!/usr/bin/env python3
"""Download GSE123902 processed dense CSVs only (<2 GB). No GSE148071."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
MAX_BYTES = 2_000_000_000


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = [
        "wget",
        "-c",
        "--tries=8",
        "--waitretry=8",
        "--timeout=60",
        "--progress=dot:giga",
        "-O",
        str(tmp),
        url,
    ]
    subprocess.check_call(cmd)
    size = tmp.stat().st_size
    if size >= MAX_BYTES:
        tmp.unlink()
        raise SystemExit(f"refusing {dest.name}: {size} bytes ≥ 2 GB")
    tmp.replace(dest)
    print(f"OK {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/pair_123902_205335_malig"))
    args = ap.parse_args()
    dest = args.outdir / "GSE123902" / "GSE123902_RAW.tar"
    print("GET GSE123902_RAW.tar", flush=True)
    wget(URL, dest)
    print("skipped GSE148071. skipped T/NK rebuild.", flush=True)


if __name__ == "__main__":
    main()
