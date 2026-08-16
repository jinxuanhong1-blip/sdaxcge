#!/usr/bin/env python3
"""Audit public EGA metadata for the requested lung-cancer trials.

This reads public metadata only; it never authenticates or downloads controlled data.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

API = "https://metadata.ega-archive.org"
TRIAL_TERMS = {
    "IMpower150": ("impower150", "impower 150"),
    "IMpower110": ("impower110", "impower 110"),
    "IMpower130": ("impower130", "impower 130"),
    "PACIFIC": ("pacific trial", "pacific study", "nct02125461"),
}
KNOWN_DATASETS = (
    "EGAD00001009725",
    "EGAD00001009726",
    "EGAD00001009764",
    "EGAD50000000273",
    "EGAD50000001814",
)


def get_json(path: str, attempts: int = 4) -> Any:
    url = f"{API}/{path.lstrip('/')}"
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "gpt-impower-public-metadata-audit/1.0"}
            )
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except (OSError, json.JSONDecodeError):
            if attempt == attempts - 1:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def object_text(item: dict[str, Any]) -> str:
    return " ".join(
        str(item.get(key) or "") for key in ("accession_id", "title", "description")
    ).lower()


def matching_trials(item: dict[str, Any]) -> list[str]:
    text = object_text(item)
    return [
        trial
        for trial, terms in TRIAL_TERMS.items()
        if any(term in text for term in terms)
    ]


def scan_collection(collection: str, page_size: int) -> tuple[list[dict[str, Any]], int]:
    matches: list[dict[str, Any]] = []
    offset = 0
    while True:
        query = urllib.parse.urlencode({"limit": page_size, "offset": offset})
        page = get_json(f"{collection}?{query}")
        for item in page:
            trials = matching_trials(item)
            if trials:
                matches.append(
                    {"collection": collection, "matched_trials": trials, **item}
                )
        offset += len(page)
        if len(page) < page_size:
            return matches, offset


def inspect_known_datasets() -> list[dict[str, Any]]:
    records = []
    for accession in KNOWN_DATASETS:
        metadata = get_json(f"datasets/{accession}")
        records.append(
            {
                "dataset": metadata,
                "studies": get_json(f"datasets/{accession}/studies"),
                "files": get_json(f"datasets/{accession}/files"),
            }
        )
    return records


def write_match_tsv(path: Path, matches: list[dict[str, Any]]) -> None:
    fields = (
        "collection",
        "matched_trials",
        "accession_id",
        "title",
        "access_type",
        "is_released",
        "description",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in matches:
            writer.writerow(
                {
                    field: (
                        ",".join(item[field])
                        if field == "matched_trials"
                        else item.get(field, "")
                    )
                    for field in fields
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/gpt_impower/ega_audit"),
    )
    parser.add_argument("--page-size", type=int, default=1000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_matches: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for collection in ("datasets", "studies"):
        matches, count = scan_collection(collection, args.page_size)
        all_matches.extend(matches)
        counts[collection] = count

    known = inspect_known_datasets()
    snapshot = {
        "api": API,
        "query_terms": TRIAL_TERMS,
        "objects_scanned": counts,
        "matches": all_matches,
        "known_impower150_datasets": known,
    }
    (args.output_dir / "ega_audit.json").write_text(
        json.dumps(snapshot, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_match_tsv(args.output_dir / "ega_name_matches.tsv", all_matches)
    print(json.dumps({"objects_scanned": counts, "matches": len(all_matches)}))


if __name__ == "__main__":
    main()
