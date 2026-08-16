#!/usr/bin/env python3
"""Download only allowlisted, public, processed files smaller than 2 GiB."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import urllib.request
from pathlib import Path

LIMIT = 2 * 1024**3
FILES = (
    {
        "id": "locuszoom-74850",
        "url": "https://my.locuszoom.org/gwas/74850/data/",
        "landing_page": "https://my.locuszoom.org/gwas/74850",
        "filename": "locuszoom_74850_summary_stats.gz",
        "access": "open",
        "processed": True,
        "scope": "Pooled GWAS including IMpower110, IMpower130, and IMpower150",
    },
    {
        "id": "locuszoom-743668",
        "url": "https://my.locuszoom.org/gwas/743668/data/",
        "landing_page": "https://my.locuszoom.org/gwas/743668",
        "filename": "locuszoom_743668_taxane_summary_stats.gz",
        "access": "open",
        "processed": True,
        "scope": "Pooled taxane-subcohort GWAS including relevant IMpower trials",
    },
)


def inspect(item: dict[str, object]) -> tuple[int, str]:
    request = urllib.request.Request(
        str(item["url"]),
        method="HEAD",
        headers={"User-Agent": "gpt-impower-open-processed-downloader/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        size = int(response.headers.get("Content-Length", "-1"))
        media_type = response.headers.get_content_type()
    if item["access"] != "open" or item["processed"] is not True:
        raise RuntimeError(f"{item['id']}: not explicitly open and processed")
    if size < 0:
        raise RuntimeError(f"{item['id']}: server did not provide Content-Length")
    if size >= LIMIT:
        raise RuntimeError(f"{item['id']}: {size} bytes exceeds the <2 GiB rule")
    if media_type not in {"application/gzip", "application/octet-stream"}:
        raise RuntimeError(f"{item['id']}: unexpected media type {media_type}")
    return size, media_type


def download(item: dict[str, object], destination: Path, expected: int) -> str:
    partial = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    request = urllib.request.Request(
        str(item["url"]),
        headers={"User-Agent": "gpt-impower-open-processed-downloader/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response, partial.open(
            "wb"
        ) as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
                digest.update(chunk)
        actual = partial.stat().st_size
        if actual != expected:
            raise RuntimeError(
                f"{item['id']}: expected {expected} bytes, downloaded {actual}"
            )
        os.replace(partial, destination)
        return digest.hexdigest()
    finally:
        partial.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/gpt_impower/open_processed"),
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Validate public headers without downloading payloads.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for item in FILES:
        size, media_type = inspect(item)
        destination = args.output_dir / str(item["filename"])
        sha256 = ""
        status = "validated_not_downloaded"
        if not args.metadata_only:
            sha256 = download(item, destination, size)
            status = "downloaded"
        rows.append(
            {
                **item,
                "size_bytes": size,
                "media_type": media_type,
                "status": status,
                "sha256": sha256,
                "local_path": str(destination),
            }
        )

    fields = tuple(rows[0])
    with (args.output_dir / "manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"files": len(rows), "status": rows[0]["status"]}))


if __name__ == "__main__":
    main()
