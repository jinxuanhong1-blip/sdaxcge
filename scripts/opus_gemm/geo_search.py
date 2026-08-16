#!/usr/bin/env python3
"""Query GEO (DataSets db) for candidate mouse lung-cancer immunotherapy studies.

Writes one TSV row per GEO series with the metadata needed to triage
candidates by hand: accession, title, summary, platform, sample count and
supplementary-file types.

Usage:
    python3 scripts/opus_gemm/geo_search.py --out results/opus_gemm/geo_candidates.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = [
    # Autochthonous GEMMs
    'KrasLSL-G12D AND p53 AND lung AND (PD-1 OR PD-L1 OR CTLA-4 OR immunotherapy)',
    '(KP OR "Kras;p53") AND lung adenocarcinoma AND mouse AND checkpoint blockade',
    '(RPM OR RPP OR "Rb1 Trp53") AND small cell lung',
    'small cell lung cancer mouse AND (PD-1 OR PD-L1 OR immunotherapy)',
    'EGFR AND lung AND mouse AND (PD-1 OR PD-L1 OR immunotherapy)',
    # Transplantable syngeneic lines derived from lung GEMMs
    '344SQ',
    '344P OR 393P OR 531LN2',
    'CMT167',
    'CMT167 AND (PD-1 OR PD-L1)',
    'KPB25L OR KP7B OR "KP cell line" AND lung',
    'orthotopic lung mouse AND anti-PD-1 AND RNA-seq',
    'lung adenocarcinoma mouse AND immune checkpoint AND single cell',
    # Marker-driven
    'Tacstd2 AND lung',
    'Cldn4 AND lung',
    'Trop2 AND lung AND mouse',
]


def eutils(endpoint: str, params: dict) -> dict | str:
    params = dict(params)
    params.setdefault("retmode", "json")
    url = f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                raw = fh.read().decode("utf-8", "replace")
            return json.loads(raw) if params.get("retmode") == "json" else raw
        except Exception as exc:  # transient NCBI throttling / network blips
            if attempt == 4:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def esearch(term: str, retmax: int = 60) -> list[str]:
    res = eutils("esearch", {"db": "gds", "term": term, "retmax": retmax})
    return res["esearchresult"].get("idlist", [])


def esummary(uids: list[str]) -> list[dict]:
    if not uids:
        return []
    res = eutils("esummary", {"db": "gds", "id": ",".join(uids)})
    result = res.get("result", {})
    return [result[u] for u in result.get("uids", [])]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seen: dict[str, dict] = {}
    for term in QUERIES:
        print(f"[search] {term}", file=sys.stderr)
        uids = esearch(term)
        time.sleep(0.4)
        for rec in esummary(uids):
            acc = rec.get("accession", "")
            if not acc.startswith("GSE"):
                continue  # skip GSM/GPL hits
            hit = seen.setdefault(
                acc,
                {
                    "accession": acc,
                    "title": rec.get("title", "").replace("\t", " "),
                    "taxon": rec.get("taxon", ""),
                    "gdstype": rec.get("gdsType", ""),
                    "n_samples": rec.get("n_samples", ""),
                    "pdat": rec.get("PDAT", ""),
                    "suppfile": rec.get("suppFile", ""),
                    "summary": rec.get("summary", "").replace("\t", " ").replace("\n", " "),
                    "queries": [],
                },
            )
            hit["queries"].append(term)
        time.sleep(0.4)

    rows = sorted(seen.values(), key=lambda r: r["accession"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "accession", "title", "taxon", "gdstype", "n_samples",
                "pdat", "suppfile", "queries", "summary",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in rows:
            r = dict(r)
            r["queries"] = " | ".join(r["queries"])
            w.writerow(r)
    print(f"wrote {len(rows)} series -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
