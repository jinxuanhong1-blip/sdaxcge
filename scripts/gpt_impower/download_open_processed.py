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
        "id": "S-EPMC10115641-MOESM2",
        "url": "https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/641/S-EPMC10115641/Files/41591_2023_2226_MOESM2_ESM.xlsx",
        "landing_page": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC10115641",
        "filename": "impower150_ctdna_supplementary_tables.xlsx",
        "access": "open",
        "processed": True,
        "scope": "Direct IMpower150 ctDNA cohort and analysis summary tables",
    },
    {
        "id": "S-EPMC11316765-MOESM4",
        "url": "https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/765/S-EPMC11316765/Files/41467_2024_51316_MOESM4_ESM.xlsx",
        "landing_page": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11316765",
        "filename": "impower150_ctdna_source_data.xlsx",
        "access": "open",
        "processed": True,
        "scope": "Direct IMpower150 figure-level ctDNA source data",
    },
    {
        "id": "S-EPMC12775477-MOESM4",
        "url": "https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/477/S-EPMC12775477/Files/41467_2025_66803_MOESM4_ESM.xlsx",
        "landing_page": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12775477",
        "filename": "impower150_transcriptomic_source_data.xlsx",
        "access": "open",
        "processed": True,
        "scope": "Direct IMpower150 figure-level transcriptomic source data",
    },
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
    {
        "id": "PGS000759",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000759/ScoringFiles/PGS000759.txt.gz",
        "landing_page": "https://www.pgscatalog.org/score/PGS000759/",
        "filename": "PGS000759.txt.gz",
        "access": "open",
        "processed": True,
        "scope": "Pooled hypothyroidism PRS evaluated across trials including IMpower130",
    },
    {
        "id": "PGS000760",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000760/ScoringFiles/PGS000760.txt.gz",
        "landing_page": "https://www.pgscatalog.org/score/PGS000760/",
        "filename": "PGS000760.txt.gz",
        "access": "open",
        "processed": True,
        "scope": "Pooled hypothyroidism PRS evaluated across trials including IMpower130",
    },
    {
        "id": "PGS000761",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000761/ScoringFiles/PGS000761.txt.gz",
        "landing_page": "https://www.pgscatalog.org/score/PGS000761/",
        "filename": "PGS000761.txt.gz",
        "access": "open",
        "processed": True,
        "scope": "Pooled hypothyroidism PRS evaluated across trials including IMpower130",
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
    allowed_media_types = {
        "application/gzip",
        "application/octet-stream",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        "application/x-gzip",
    }
    if media_type not in allowed_media_types:
        raise RuntimeError(f"{item['id']}: unexpected media type {media_type}")
    return size, media_type


def download(item: dict[str, object], destination: Path, expected: int) -> str:
    if destination.exists() and destination.stat().st_size == expected:
        digest = hashlib.sha256()
        with destination.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()

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
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"files": len(rows), "status": rows[0]["status"]}))


if __name__ == "__main__":
    main()
