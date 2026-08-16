#!/usr/bin/env python3
"""Replace the multi-megabyte probe dumps with a slim, commitable leftover audit table.

The full `geo_probe.json` / search dumps stay on the local cache and can be regenerated
by scripts 01–03. What we keep in git is one row per leftover-relevant accession.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"

DROP = (
    "geo_probe.json",
    "geo_search_candidates.json",
    "geo_triage_scored.json",
    "geo_triage_shortlist.json",
    "dataset_excluded.json",
    "dataset_shortlist.json",
    "supp_clinical_scan.json",
)


def main() -> int:
    probe_path = OUT / "geo_probe.json"
    if not probe_path.exists():
        print("[slim] geo_probe.json already gone")
        return 0
    probe = json.loads(probe_path.read_text())
    slim = []
    for r in probe:
        slim.append(
            {
                "accession": r["accession"],
                "title": r["title"],
                "n_characteristics": max(
                    (len(c.split("\t")) - 1 for c in r["characteristics"]), default=0
                ),
                "has_response_annotation": r["has_response_annotation"],
                "mentions_ici": r["mentions_ici"],
                "matrix_data_rows": r["matrix_data_rows"],
                "n_supp_files": r["n_supp_files"],
                "processed_candidates": [f["name"] for f in r["processed_candidates"]],
                "oversized_files": [f["name"] for f in r["oversized_files"]],
                "pubmed": r["pubmed"],
                "platform_ids": r["platform_ids"],
            }
        )
    (OUT / "geo_probe_slim.json").write_text(json.dumps(slim, indent=2))
    for name in DROP:
        p = OUT / name
        if p.exists():
            p.unlink()
            print(f"[slim] dropped {name}")
    print(f"[slim] kept geo_probe_slim.json ({len(slim)} series)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
