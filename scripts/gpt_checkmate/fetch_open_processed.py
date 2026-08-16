#!/usr/bin/env python3
"""Download only the verified open processed CheckMate files below 2 GB."""

from __future__ import annotations

import csv
import hashlib
import sys
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

MAX_BYTES = 2_000_000_000
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "gpt_checkmate" / "open_processed"
MANIFEST = ROOT / "results" / "gpt_checkmate" / "open_processed_manifest.tsv"
SOURCES = [
    {
        "trial": "153",
        "name": "checkmate153_supplementary_table_1.xlsx",
        "description": "Differential-expression analysis summarized results",
        "url": (
            "https://media.springernature.com/original/springer-static/esm/"
            "art%3A10.1038%2Fs41591-024-03240-y/MediaObjects/"
            "41591_2024_3240_MOESM3_ESM.xlsx"
        ),
        "doi": "10.1038/s41591-024-03240-y",
    },
    {
        "trial": "153",
        "name": "checkmate153_supplementary_table_2.xlsx",
        "description": "Full immunogenicity-screen table used for figure generation",
        "url": (
            "https://media.springernature.com/original/springer-static/esm/"
            "art%3A10.1038%2Fs41591-024-03240-y/MediaObjects/"
            "41591_2024_3240_MOESM4_ESM.xlsx"
        ),
        "doi": "10.1038/s41591-024-03240-y",
    },
]


def download(source: dict[str, str]) -> tuple[Path, int, str]:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    destination = OUTPUT / source["name"]
    temporary = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    size = 0
    request = urllib.request.Request(
        source["url"], headers={"User-Agent": "gpt-checkmate-downloader/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) >= MAX_BYTES:
                raise ValueError(f"{source['name']} declares {declared} bytes (limit {MAX_BYTES})")
            with temporary.open("wb") as handle:
                while chunk := response.read(1024 * 1024):
                    size += len(chunk)
                    if size >= MAX_BYTES:
                        raise ValueError(f"{source['name']} reached the {MAX_BYTES}-byte limit")
                    digest.update(chunk)
                    handle.write(chunk)
        if not zipfile.is_zipfile(temporary):
            raise ValueError(f"{source['name']} is not a valid XLSX/ZIP file")
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return destination, size, digest.hexdigest()


def workbook_summary(path: Path) -> tuple[int, str]:
    namespace = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    names = [
        sheet.attrib["name"]
        for sheet in workbook.findall("m:sheets/m:sheet", namespaces=namespace)
    ]
    return len(names), ";".join(names)


def main() -> int:
    rows: list[dict[str, object]] = []
    for source in SOURCES:
        path, size, sha256 = download(source)
        sheet_count, sheet_names = workbook_summary(path)
        rows.append(
            {
                "trial": source["trial"],
                "description": source["description"],
                "url": source["url"],
                "doi": source["doi"],
                "path": path.relative_to(ROOT).as_posix(),
                "bytes": size,
                "sha256": sha256,
                "sheet_count": sheet_count,
                "sheet_names": sheet_names,
                "downloaded_at": datetime.now(timezone.utc)
                .replace(microsecond=0)
                .isoformat(),
                "access": "open",
                "data_level": "processed",
            }
        )

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "trial",
        "path",
        "bytes",
        "sha256",
        "access",
        "data_level",
        "description",
        "sheet_count",
        "sheet_names",
        "doi",
        "url",
        "downloaded_at",
    ]
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Downloaded {len(rows)} open processed files; manifest: {MANIFEST}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
