#!/usr/bin/env python3
"""Fetch GEO series SOFT (full) for candidate accessions and summarize supplementary files."""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

OUT = Path("results/hunt_neoadj/metadata")
CANDIDATES = [
    "GSE207422",
    "GSE243013",
    "GSE146100",
    "GSE229353",
    "GSE248378",
    "GSE241934",
    "GSE225620",
    "GSE280232",
    "GSE274934",
    "GSE291670",
    "GSE299111",
    "GSE300685",
    "GSE337519",
    "GSE249568",
    "GSE250509",
    "GSE240033",
]


def fetch(acc: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full"
    req = urllib.request.Request(url, headers={"User-Agent": "hunt-neoadj/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", errors="replace")


def parse(text: str) -> dict:
    fields = {}
    samples = []
    suppl = []
    pubmed = []
    for line in text.splitlines():
        if not line.startswith("!"):
            continue
        if " = " not in line:
            continue
        k, v = line[1:].split(" = ", 1)
        if k == "Series_sample_id":
            samples.append(v)
        elif k == "Series_supplementary_file":
            suppl.append(v)
        elif k == "Series_pubmed_id":
            pubmed.append(v)
        else:
            fields.setdefault(k, []).append(v)
    return {
        "title": " | ".join(fields.get("Series_title", [])),
        "type": " | ".join(fields.get("Series_type", [])),
        "status": " | ".join(fields.get("Series_status", [])),
        "submission_date": " | ".join(fields.get("Series_submission_date", [])),
        "last_update": " | ".join(fields.get("Series_last_update_date", [])),
        "pubmed": pubmed,
        "n_samples": len(samples),
        "platform": fields.get("Series_platform_id", []),
        "overall_design": " ".join(fields.get("Series_overall_design", [])),
        "summary": " ".join(fields.get("Series_summary", [])),
        "supplementary": suppl,
        "organism": fields.get("Series_sample_organism", fields.get("Series_sample_taxid", [])),
    }


def main():
    summaries = {}
    for acc in CANDIDATES:
        print("fetch", acc, flush=True)
        try:
            text = fetch(acc)
        except Exception as e:
            summaries[acc] = {"error": str(e)}
            print("  ERROR", e, flush=True)
            time.sleep(0.4)
            continue
        (OUT / f"{acc}.full.soft.txt").write_text(text)
        summaries[acc] = parse(text)
        print(
            f"  n={summaries[acc]['n_samples']} type={summaries[acc]['type'][:80]} "
            f"suppl={len(summaries[acc]['supplementary'])}",
            flush=True,
        )
        time.sleep(0.35)
    (OUT / "candidate_series_summary.json").write_text(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
