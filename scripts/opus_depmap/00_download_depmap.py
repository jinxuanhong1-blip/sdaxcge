#!/usr/bin/env python3
"""Download the open DepMap 24Q4 Public files used by this slice.

Source of truth is the figshare record for the DepMap 24Q4 Public release
(article 27993248, DOI 10.6084/m9.figshare.27993248). Only openly licensed
release files are fetched. Each file is verified against the md5 supplied by
the figshare API and a manifest with sizes/checksums is written next to the
data so the analysis is reproducible from a clean checkout.

Usage:
    python scripts/opus_depmap/00_download_depmap.py [--outdir data/opus_depmap]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

FIGSHARE_ARTICLE = 27993248
FIGSHARE_DOI = "10.6084/m9.figshare.27993248"
RELEASE = "DepMap Public 24Q4"
FILES_API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}/files?page_size=200"
DOWNLOAD_URL = "https://ndownloader.figshare.com/files/{file_id}"

WANTED = [
    "Model.csv",
    "OmicsProfiles.csv",
    "CRISPRGeneEffect.csv",
    "CRISPRGeneDependency.csv",
    "CRISPRInferredCommonEssentials.csv",
    "AchillesCommonEssentialControls.csv",
    "AchillesNonessentialControls.csv",
    "CRISPRScreenMap.csv",
    "ScreenGeneEffect.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "README.txt",
]


def fetch_file_index() -> dict[str, dict]:
    with urllib.request.urlopen(FILES_API, timeout=120) as resp:
        entries = json.load(resp)
    return {e["name"]: e for e in entries}


def md5sum(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path, attempts: int = 4) -> None:
    for attempt in range(1, attempts + 1):
        try:
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url, timeout=300) as resp, tmp.open("wb") as out:
                while True:
                    block = resp.read(1 << 22)
                    if not block:
                        break
                    out.write(block)
            tmp.replace(dest)
            return
        except Exception as exc:  # network flakiness on large files
            wait = 4 * 2 ** (attempt - 1)
            print(f"    attempt {attempt} failed ({exc}); retrying in {wait}s", flush=True)
            if attempt == attempts:
                raise
            time.sleep(wait)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/opus_depmap")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    index = fetch_file_index()
    manifest = {
        "release": RELEASE,
        "figshare_article": FIGSHARE_ARTICLE,
        "figshare_doi": FIGSHARE_DOI,
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "files": [],
    }

    for name in WANTED:
        if name not in index:
            print(f"!! {name} not present in figshare record", flush=True)
            return 2
        entry = index[name]
        dest = outdir / name
        expected = entry["supplied_md5"]
        if dest.exists() and md5sum(dest) == expected:
            print(f"== {name} already present and verified", flush=True)
        else:
            size_mb = entry["size"] / 1e6
            print(f">> downloading {name} ({size_mb:.1f} MB)", flush=True)
            download(DOWNLOAD_URL.format(file_id=entry["id"]), dest)
            got = md5sum(dest)
            if got != expected:
                print(f"!! md5 mismatch for {name}: {got} != {expected}", flush=True)
                return 3
            print(f"   ok {name}", flush=True)
        manifest["files"].append(
            {
                "name": name,
                "figshare_file_id": entry["id"],
                "size_bytes": entry["size"],
                "md5": expected,
                "url": DOWNLOAD_URL.format(file_id=entry["id"]),
            }
        )

    (outdir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("manifest written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
