#!/usr/bin/env python3
"""Write a compact leftover catalog TSV from the search manifest."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("results/w200/AE_leftover/search_manifest.json"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/w200/AE_leftover/catalog.tsv"),
    )
    args = parser.parse_args()
    payload = json.loads(args.manifest.read_text())
    rows = []
    for s in payload["studies"]:
        rows.append(
            {
                "accession": s["accession"],
                "decision": s["decision"],
                "title": s.get("title") or "",
                "n_processed": s.get("n_processed"),
                "processed_bytes": s.get("processed_bytes"),
                "geo_matches": ";".join(s.get("geo_matches") or []),
                "rationale": s.get("rationale") or "",
                "study_url": s.get("study_url") or "",
            }
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out}: {len(rows)} rows")


if __name__ == "__main__":
    main()
