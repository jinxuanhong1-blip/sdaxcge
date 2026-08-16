#!/usr/bin/env python3
"""Scan EGA public metadata for the requested CheckMate lung trials."""

from __future__ import annotations

import csv
import json
import sys
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

API = "https://metadata.ega-archive.org"
OUTPUT = Path(__file__).resolve().parents[2] / "results" / "gpt_checkmate"
LIMIT = 100_000
TRIAL_ALIASES = {
    "017": ["checkmate 017", "checkmate-017", "checkmate017", "ca209-017", "ca209017"],
    "057": ["checkmate 057", "checkmate-057", "checkmate057", "ca209-057", "ca209057"],
    "227": ["checkmate 227", "checkmate-227", "checkmate227", "ca209-227", "ca209227"],
    "9LA": ["checkmate 9la", "checkmate-9la", "checkmate9la", "ca209-9la", "ca2099la"],
    "816": ["checkmate 816", "checkmate-816", "checkmate816", "ca209-816", "ca209816"],
    "153": ["checkmate 153", "checkmate-153", "checkmate153", "ca209-153", "ca209153"],
}


def get_json(path: str) -> object:
    request = urllib.request.Request(
        f"{API}{path}", headers={"User-Agent": "gpt-checkmate-catalog/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def searchable_text(record: dict[str, object]) -> str:
    return " ".join(
        str(record.get(field) or "") for field in ("title", "description", "abstract")
    ).lower()


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    collections: dict[str, list[dict[str, object]]] = {}
    matches: dict[str, dict[str, list[dict[str, object]]]] = {
        trial: {"studies": [], "datasets": []} for trial in TRIAL_ALIASES
    }

    for kind in ("studies", "datasets"):
        records = get_json(f"/{kind}?limit={LIMIT}&offset=0")
        if not isinstance(records, list):
            raise TypeError(f"Unexpected EGA response for {kind}")
        if len(records) >= LIMIT:
            raise RuntimeError(f"EGA {kind} response reached scan limit {LIMIT}")
        collections[kind] = records
        for record in records:
            text = searchable_text(record)
            for trial, aliases in TRIAL_ALIASES.items():
                if any(alias in text for alias in aliases):
                    matches[trial][kind].append(record)

    dataset_ids = sorted(
        {
            record["accession_id"]
            for trial in matches.values()
            for record in trial["datasets"]
        }
    )
    details: dict[str, object] = {}
    for accession in dataset_ids:
        dataset = get_json(f"/datasets/{accession}")
        studies = get_json(f"/datasets/{accession}/studies")
        dacs = get_json(f"/datasets/{accession}/dacs")
        files = get_json(f"/datasets/{accession}/files")
        if not isinstance(files, list):
            raise TypeError(f"Unexpected EGA file response for {accession}")
        sizes = [int(row["filesize"]) for row in files]
        extensions = Counter(str(row["extension"]) for row in files)
        details[accession] = {
            "dataset": dataset,
            "studies": studies,
            "dacs": dacs,
            "file_summary": {
                "count": len(files),
                "total_bytes": sum(sizes),
                "min_bytes": min(sizes) if sizes else None,
                "max_bytes": max(sizes) if sizes else None,
                "extensions": dict(sorted(extensions.items())),
                "files_under_2GB_decimal": sum(size < 2_000_000_000 for size in sizes),
                "note": "File-level records are metadata only; controlled files were not downloaded.",
            },
        }

    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    catalog = {
        "generated_at": generated,
        "api": API,
        "scan": {
            "method": "case-insensitive literal alias matching over title, description, and abstract",
            "study_records_scanned": len(collections["studies"]),
            "dataset_records_scanned": len(collections["datasets"]),
            "aliases": TRIAL_ALIASES,
        },
        "matches": matches,
        "dataset_details": details,
    }
    (OUTPUT / "ega_catalog.json").write_text(
        json.dumps(catalog, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    with (OUTPUT / "ega_catalog.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "trial",
                "study_accessions",
                "dataset_accessions",
                "access_type",
                "sample_count",
                "file_count",
                "total_bytes",
                "file_types",
            ]
        )
        for trial in TRIAL_ALIASES:
            studies = matches[trial]["studies"]
            datasets = matches[trial]["datasets"]
            dataset_accessions = [str(row["accession_id"]) for row in datasets]
            summaries = [details[accession] for accession in dataset_accessions]
            writer.writerow(
                [
                    trial,
                    ";".join(str(row["accession_id"]) for row in studies) or "-",
                    ";".join(dataset_accessions) or "-",
                    ";".join(
                        str(item["dataset"].get("access_type") or "") for item in summaries
                    )
                    or "-",
                    sum(int(item["dataset"].get("num_samples") or 0) for item in summaries),
                    sum(int(item["file_summary"]["count"]) for item in summaries),
                    sum(int(item["file_summary"]["total_bytes"]) for item in summaries),
                    ";".join(
                        ",".join(item["file_summary"]["extensions"]) for item in summaries
                    )
                    or "-",
                ]
            )

    print(
        f"Wrote EGA catalog: {len(collections['studies'])} studies, "
        f"{len(collections['datasets'])} datasets scanned"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
