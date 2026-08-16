#!/usr/bin/env python3
"""Search NCBI GEO (and log ArrayExpress queries) for CLDN4 / TACSTD2 perturbation series.

Writes:
  results/kdko/geo_search_raw.tsv   -- every GEO DataSet/Series hit with esummary metadata
  results/kdko/geo_search_queries.tsv -- query -> hit count audit trail

No accession is ever synthesised: every row here comes straight from E-utilities.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results", "kdko")
OUT = os.path.abspath(OUT)
TOOL = {"tool": "cldn4_trop2_kdko", "email": "cloud-agent@example.org"}

# Gene-level free-text queries. Kept deliberately broad; filtering happens later.
QUERIES = [
    # --- CLDN4 ---
    "CLDN4 AND knockdown", "CLDN4 AND knockout", "CLDN4 AND shRNA", "CLDN4 AND siRNA",
    "CLDN4 AND CRISPR", "CLDN4 AND silencing", "CLDN4 AND overexpression",
    "claudin-4 AND knockdown", "claudin 4 AND knockout", "claudin-4 AND siRNA",
    "Cldn4 AND mouse AND knockout", "Cldn4 AND deficient",
    # --- TACSTD2 / TROP2 ---
    "TACSTD2 AND knockdown", "TACSTD2 AND knockout", "TACSTD2 AND shRNA",
    "TACSTD2 AND siRNA", "TACSTD2 AND CRISPR", "TACSTD2 AND silencing",
    "TACSTD2 AND overexpression",
    "TROP2 AND knockdown", "TROP2 AND knockout", "TROP2 AND shRNA", "TROP2 AND siRNA",
    "TROP2 AND CRISPR", "TROP-2 AND knockdown", "Trop2 AND depletion",
    "Tacstd2 AND mouse AND knockout", "Tacstd2 AND deficient",
    "TROP2 AND sacituzumab",  # resistance models often carry TACSTD2 KO arms
    # --- broader tight junction / claudin context (for the "controlled sets" list) ---
    "claudin AND knockout AND lung",
    "tight junction AND knockdown AND interferon",
]


def eutils(endpoint: str, **params) -> dict | str:
    params = {**TOOL, **params}
    url = f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                raw = fh.read().decode("utf-8", "replace")
            if params.get("retmode") == "json":
                return json.loads(raw)
            return raw
        except Exception as exc:  # noqa: BLE001
            if attempt == 4:
                raise
            sys.stderr.write(f"retry {attempt} {endpoint}: {exc}\n")
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def search(term: str) -> list[str]:
    res = eutils("esearch", db="gds", term=term, retmax=200, retmode="json")
    return res["esearchresult"].get("idlist", [])


def summary(uids: list[str]) -> list[dict]:
    if not uids:
        return []
    out = []
    for i in range(0, len(uids), 100):
        chunk = uids[i:i + 100]
        res = eutils("esummary", db="gds", id=",".join(chunk), retmode="json")
        result = res.get("result", {})
        for uid in result.get("uids", []):
            out.append(result[uid])
        time.sleep(0.4)
    return out


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    seen: dict[str, dict] = {}
    query_log = []

    for q in QUERIES:
        uids = search(q)
        query_log.append({"query": q, "n_hits": len(uids)})
        print(f"{len(uids):4d}  {q}", flush=True)
        for rec in summary(uids):
            acc = rec.get("accession", "")
            if not acc:
                continue
            if acc in seen:
                seen[acc]["matched_queries"].add(q)
                continue
            rec["matched_queries"] = {q}
            seen[acc] = rec
        time.sleep(0.4)

    fields = [
        "accession", "gds_type", "entrytype", "taxon", "n_samples", "gpl", "pdat",
        "title", "summary", "suppfile", "ftplink", "matched_queries",
    ]
    with open(os.path.join(OUT, "geo_search_raw.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for acc, rec in sorted(seen.items()):
            row = {k: rec.get(k, "") for k in fields}
            row["accession"] = acc
            row["n_samples"] = rec.get("n_samples", "")
            row["matched_queries"] = ";".join(sorted(rec["matched_queries"]))
            for k, v in row.items():
                if isinstance(v, str):
                    row[k] = " ".join(v.split())
            w.writerow(row)

    with open(os.path.join(OUT, "geo_search_queries.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["query", "n_hits"], delimiter="\t")
        w.writeheader()
        w.writerows(query_log)

    print(f"\n{len(seen)} unique GEO entries -> {OUT}/geo_search_raw.tsv")


if __name__ == "__main__":
    main()
