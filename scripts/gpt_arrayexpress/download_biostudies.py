#!/usr/bin/env python3
"""Download eligible open processed files from a generated search manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path

LIMIT = 2_000_000_000
DOWNLOAD_DECISIONS = {"include", "include_context", "include_target"}
METADATA_DECISIONS = DOWNLOAD_DECISIONS | {"metadata_only", "metadata_only_target"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, expected_size: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size == expected_size:
        return
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "gpt-arrayexpress/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as output:
                while chunk := response.read(4 * 1024 * 1024):
                    output.write(chunk)
            if partial.stat().st_size != expected_size:
                raise IOError(f"size mismatch for {destination}: {partial.stat().st_size} != {expected_size}")
            partial.replace(destination)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("results/gpt_arrayexpress/search_manifest.json"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/gpt_arrayexpress/downloads"))
    parser.add_argument("--download-manifest", type=Path, default=Path("results/gpt_arrayexpress/download_manifest.json"))
    args = parser.parse_args()

    search = json.loads(args.manifest.read_text())
    records = []
    for study in search["studies"]:
        if study["decision"] not in METADATA_DECISIONS:
            continue
        accession = study["accession"]
        eligible = []
        for item in study["all_files"]:
            is_metadata = item["path"].endswith((".idf.txt", ".sdrf.txt"))
            if item["processed"] and study["decision"] in DOWNLOAD_DECISIONS or is_metadata:
                eligible.append({**item, "role": "processed" if item["processed"] else "metadata"})
        for item in eligible:
            if not (0 < item["size"] < LIMIT):
                records.append({**item, "accession": accession, "status": "skipped_size"})
                continue
            destination = args.out_dir / accession / item["path"]
            download(item["download_url"], destination, item["size"])
            records.append({
                **item,
                "accession": accession,
                "status": "downloaded",
                "local_path": str(destination),
                "sha256": sha256(destination),
            })
            print(f"{accession}: {item['path']} ({item['size']} bytes)")

    payload = {
        "limit_bytes_strictly_less_than": LIMIT,
        "file_count": sum(r["status"] == "downloaded" for r in records),
        "total_bytes": sum(r["size"] for r in records if r["status"] == "downloaded"),
        "files": records,
    }
    args.download_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.download_manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {args.download_manifest}")


if __name__ == "__main__":
    main()
