#!/usr/bin/env python3
"""Step 1 - fetch every public input this slice depends on.

Downloads go to an out-of-tree cache (``OPUS_LUSC_CACHE``); only a manifest with
sizes and SHA-256 digests is written into the repository so that the exact
inputs can be re-identified later.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

from config import COHORTS, MANIFEST_DIR, RAW_DIR, TCGA_FILES

USER_AGENT = "opus-lusc-ici-slice/1.0 (public data reanalysis)"
CHUNK = 1 << 20


def fetch(url: str, dest_name: str, retries: int = 4) -> dict:
    dest = RAW_DIR / dest_name
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached  {dest_name} ({dest.stat().st_size/1e6:.1f} MB)")
        return describe(url, dest)

    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as resp, open(dest, "wb") as fh:
                while True:
                    block = resp.read(CHUNK)
                    if not block:
                        break
                    fh.write(block)
            print(f"  fetched {dest_name} ({dest.stat().st_size/1e6:.1f} MB)")
            return describe(url, dest)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_err = exc
            wait = 4 * (2**attempt)
            print(f"  retry {attempt+1}/{retries} for {dest_name} after {wait}s ({exc})")
            time.sleep(wait)
    raise RuntimeError(f"could not download {url}: {last_err}")


def describe(url: str, path) -> dict:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            digest.update(block)
    return {
        "url": url,
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def main() -> int:
    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cache_dir": str(RAW_DIR),
        "files": [],
    }

    for cohort in COHORTS:
        print(f"[{cohort.key}] {cohort.gse} - {cohort.label}")
        matrix = f"{cohort.gse}_series_matrix.txt.gz"
        manifest["files"].append(
            fetch(f"{cohort.ftp_dir}/matrix/{matrix}", matrix) | {"cohort": cohort.key, "kind": "metadata"}
        )
        for supp in cohort.supplementary:
            manifest["files"].append(
                fetch(f"{cohort.ftp_dir}/suppl/{supp}", supp) | {"cohort": cohort.key, "kind": "expression"}
            )

    print("[TCGA] UCSC Xena GDC hub")
    for key, url in TCGA_FILES.items():
        name = url.rsplit("/", 1)[-1]
        manifest["files"].append(fetch(url, name) | {"cohort": "TCGA", "kind": key})

    total = sum(f["bytes"] for f in manifest["files"])
    manifest["total_bytes"] = total
    print(f"\nTotal raw download: {total/1e9:.2f} GB across {len(manifest['files'])} files")

    out = MANIFEST_DIR / "raw_downloads.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Manifest written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
