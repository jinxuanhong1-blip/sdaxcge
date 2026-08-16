#!/usr/bin/env python3
"""Round-2 GEO search: sample-level and lung-focused queries + ArrayExpress/BioStudies.

Round 1 (01_search_geo.py) queried gene+method free text on the `gds` database. Many
TACSTD2 hits there were false positives (gene mentioned only in a supplementary gene
list). Round 2 does two extra things:
  1. queries the `gds` DataSets index with sample-title-style and tissue-specific terms;
  2. queries EBI BioStudies (ArrayExpress collection) for the same gene symbols.

Writes results/kdko/geo_search_round2.tsv and results/kdko/arrayexpress_search.tsv
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "results", "kdko"))
TOOL = {"tool": "cldn4_trop2_kdko", "email": "cloud-agent@example.org"}

GEO_QUERIES = [
    # sample-title style (GEO indexes GSM titles in the gds db)
    'CLDN4[Description] AND ("knock down"[All Fields] OR knockdown[All Fields])',
    'TACSTD2[Description]', 'CLDN4[Description]',
    '"claudin-4"[Title]', '"claudin 4"[Title]', 'TROP2[Title]', 'TROP-2[Title]',
    'TACSTD2[Title]', 'CLDN4[Title]',
    # lung focus
    'CLDN4 AND lung AND (knockdown OR knockout OR shRNA OR siRNA OR CRISPR)',
    'TACSTD2 AND lung AND (knockdown OR knockout OR shRNA OR siRNA OR CRISPR)',
    'TROP2 AND lung AND (knockdown OR knockout OR CRISPR)',
    'Cldn4 AND lung AND mouse',
    'Tacstd2 AND lung AND mouse',
    # other perturbation vocabularies
    'CLDN4 AND (depletion OR ablation OR "loss of function")',
    'TACSTD2 AND (depletion OR ablation OR "loss of function")',
    'TROP2 AND immune AND exclusion',
    'claudin AND immune exclusion',
    'TACSTD2 AND interferon',
    'CLDN4 AND interferon',
    'TROP2 AND checkpoint blockade',
    # degrader / inducible systems
    'CLDN4 AND doxycycline', 'TACSTD2 AND doxycycline',
    'TROP2 AND organoid AND knockout',
]

AE_QUERIES = ["CLDN4 knockdown", "CLDN4 knockout", "claudin-4 siRNA",
              "TACSTD2 knockdown", "TACSTD2 knockout", "TROP2 knockdown",
              "TROP2 knockout", "Trop2 shRNA"]

GENE_PAT = re.compile(
    r"\b(CLDN4|claudin[- ]?4|TACSTD2|TROP[- ]?2|GA733|EGP-?1|M1S1)\b", re.I)


def fetch(url: str, as_json: bool = True):
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cldn4-trop2-kdko/1.0"})
            with urllib.request.urlopen(req, timeout=90) as fh:
                raw = fh.read().decode("utf-8", "replace")
            return json.loads(raw) if as_json else raw
        except Exception as exc:  # noqa: BLE001
            if attempt == 4:
                sys.stderr.write(f"FAILED {url}: {exc}\n")
                return None
            time.sleep(2 ** attempt)
    return None


def eutils(endpoint: str, **params):
    return fetch(f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode({**TOOL, **params}))


def geo_round2() -> None:
    seen: dict[str, dict] = {}
    for q in GEO_QUERIES:
        res = eutils("esearch", db="gds", term=q, retmax=200, retmode="json")
        uids = res["esearchresult"].get("idlist", []) if res else []
        print(f"{len(uids):4d}  {q}", flush=True)
        for i in range(0, len(uids), 100):
            s = eutils("esummary", db="gds", id=",".join(uids[i:i + 100]), retmode="json")
            if not s:
                continue
            r = s.get("result", {})
            for uid in r.get("uids", []):
                rec = r[uid]
                acc = rec.get("accession", "")
                if not acc:
                    continue
                if acc in seen:
                    seen[acc]["matched_queries"].add(q)
                else:
                    rec["matched_queries"] = {q}
                    seen[acc] = rec
            time.sleep(0.4)
        time.sleep(0.4)

    fields = ["accession", "entrytype", "taxon", "n_samples", "gds_type", "gpl",
              "pdat", "title", "summary", "suppfile", "gene_mentioned", "matched_queries"]
    with open(os.path.join(OUT, "geo_search_round2.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for acc, rec in sorted(seen.items()):
            row = {k: rec.get(k, "") for k in fields}
            row["accession"] = acc
            row["matched_queries"] = ";".join(sorted(rec["matched_queries"]))
            row["gene_mentioned"] = bool(
                GENE_PAT.search(f"{rec.get('title','')} {rec.get('summary','')}"))
            for k, v in row.items():
                if isinstance(v, str):
                    row[k] = " ".join(v.split())
            w.writerow(row)
    print(f"round2: {len(seen)} unique entries")


def arrayexpress() -> None:
    rows = []
    for q in AE_QUERIES:
        url = ("https://www.ebi.ac.uk/biostudies/api/v1/search?query="
               + urllib.parse.quote(q) + "&collection=arrayexpress&pageSize=50")
        res = fetch(url)
        hits = res.get("hits", []) if res else []
        print(f"{len(hits):4d}  AE: {q}", flush=True)
        for h in hits:
            rows.append({
                "query": q,
                "accession": h.get("accession", ""),
                "title": " ".join(str(h.get("title", "")).split()),
                "organism": " ".join(str(h.get("organism", "")).split()),
                "release_date": h.get("release_date", ""),
                "gene_mentioned": bool(GENE_PAT.search(str(h.get("title", "")))),
            })
        time.sleep(0.5)
    with open(os.path.join(OUT, "arrayexpress_search.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else
                           ["query", "accession", "title", "organism",
                            "release_date", "gene_mentioned"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"arrayexpress: {len(rows)} rows")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    geo_round2()
    arrayexpress()
