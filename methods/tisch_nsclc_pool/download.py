#!/usr/bin/env python3
"""Download public TISCH2 NSCLC expression.h5 + CellMetainfo_table.tsv.

Does not clone TISCH or any other git repo. Writes provenance JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATASETS, TISCH_BASE


def sha256_file(path: Path, max_bytes: int | None = None) -> str:
    h = hashlib.sha256()
    n = 0
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
            n += len(chunk)
            if max_bytes is not None and n >= max_bytes:
                break
    return h.hexdigest()


def curl_download(url: str, dest: Path, retries: int = 4) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return {
            "url": url,
            "path": str(dest),
            "bytes": dest.stat().st_size,
            "skipped": True,
            "ok": True,
        }
    delay = 4
    last_err = None
    for attempt in range(1, retries + 1):
        cmd = [
            "curl",
            "-fL",
            "--retry",
            "3",
            "--retry-delay",
            "4",
            "-o",
            str(dest),
            url,
        ]
        try:
            subprocess.run(cmd, check=True)
            if dest.exists() and dest.stat().st_size > 0:
                return {
                    "url": url,
                    "path": str(dest),
                    "bytes": dest.stat().st_size,
                    "skipped": False,
                    "ok": True,
                    "attempt": attempt,
                }
            last_err = "empty file"
        except subprocess.CalledProcessError as e:
            last_err = str(e)
            if dest.exists():
                dest.unlink()
        time.sleep(delay)
        delay *= 2
    return {
        "url": url,
        "path": str(dest),
        "ok": False,
        "error": last_err,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/tisch_nsclc_pool")
    ap.add_argument("--out-dir", default="methods/tisch_nsclc_pool")
    ap.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="Optional dataset ids to download (default: all catalog download_h5=True plus all CellMetainfo)",
    )
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    wanted = set(args.only) if args.only else set(DATASETS)
    provenance: dict = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "base": TISCH_BASE,
        "files": {},
    }

    for ds, spec in DATASETS.items():
        if ds not in wanted:
            continue
        files = ["CellMetainfo_table.tsv"]
        if spec.get("download_h5", False):
            files.append("expression.h5")
        for fn in files:
            url = f"{TISCH_BASE}/{ds}/{ds}_{fn}"
            dest = data_dir / f"{ds}_{fn}"
            print(f">> {ds} {fn}", flush=True)
            rec = curl_download(url, dest)
            if rec.get("ok") and dest.exists():
                rec["sha256"] = sha256_file(dest)
            provenance["files"][f"{ds}_{fn}"] = rec
            print(
                f"   ok={rec.get('ok')} bytes={rec.get('bytes')} err={rec.get('error')}",
                flush=True,
            )

    dest = out_dir / "download_provenance.json"
    dest.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {dest}")
    n_fail = sum(1 for r in provenance["files"].values() if not r.get("ok"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
