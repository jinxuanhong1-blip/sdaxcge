#!/usr/bin/env python3
"""Re-run the 2026-09-21 CLDN4 KD/KO expression-dataset queries.

Writes a fresh query log under results/cldn4_kd_match_wave_2026/rerun/.
Does not download matrices and does not score IFN/APM. The frozen inventory
from the 2026-09-21 run is queries.tsv and inventory.tsv beside RESULTS.md.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "cldn4_kd_match_wave_2026" / "rerun"

NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1/search"
NGDC = "https://ngdc.cncb.ac.cn/search/api/specific"

GEO_QUERIES = [
    '(CLDN4 OR Cldn4 OR "claudin 4" OR "claudin-4") AND (knockdown OR knockout OR CRISPR OR shRNA OR siRNA)',
    '(CLDN4 OR Cldn4 OR "claudin-4") AND (shRNA OR siRNA OR CRISPR OR knockout OR knockdown) AND ("2024"[PDAT] : "2027"[PDAT])',
    '(CLDN4[All Fields] OR Cldn4[All Fields]) AND ("2024/01/01"[PDAT] : "2026/12/31"[PDAT])',
    '"CLDN4-/-"[All Fields] OR "CLDN4 knockout"[All Fields] OR "shCLDN4"[All Fields] OR "siCLDN4"[All Fields]',
]
SRA_QUERIES = [
    'CLDN4[All Fields] AND (knockout[All Fields] OR knockdown[All Fields] OR CRISPR[All Fields] OR shRNA[All Fields])',
    '(CLDN4 OR Cldn4) AND (knockdown OR knockout OR CRISPR OR shRNA OR siRNA) AND ("2024"[PDAT] : "2026"[PDAT])',
]
NGDC_QUERIES = [
    ("gsa", "CLDN4"),
    ("gsa", "Cldn4"),
    ("hra", "CLDN4"),
    ("hra", "Cldn4"),
    ("cra", "CLDN4"),
    ("omix", "CLDN4"),
    ("bioproject", "CLDN4"),
]


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "cldn4-kd-match-wave/1.0"})
    with urllib.request.urlopen(req, timeout=90) as response:
        return response.read()


def esearch(db: str, term: str) -> tuple[str, str]:
    url = NCBI + "?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": 20, "retmode": "json"}
    )
    result = json.loads(fetch(url))["esearchresult"]
    return str(result.get("count")), ",".join(result.get("idlist", [])[:20])


def ngdc(db: str, query: str) -> tuple[str, str]:
    url = NGDC + "?" + urllib.parse.urlencode(
        {"db": db, "q": query, "sort": "desc", "length": 15}
    )
    block = json.loads(fetch(url))["result"]["data"]
    rows = block.get("data") or []
    ids = []
    for row in rows:
        attrs = row.get("attrs") or {}
        ids.append(str(row.get("accession") or attrs.get("Accession") or row.get("id") or ""))
    return str(block.get("recordsTotal")), ",".join(ids)


def biostudies(query: str) -> tuple[str, str]:
    url = BIOSTUDIES + "?" + urllib.parse.urlencode({"query": query, "pageSize": 100, "page": 1})
    payload = json.loads(fetch(url))
    accessions = []
    for hit in payload.get("hits") or []:
        acc = hit.get("accession") or ""
        if acc.startswith("E-"):
            accessions.append(acc)
    return str(payload.get("totalHits")), ",".join(accessions)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["source\tquery\tcount\tids"]
    for term in GEO_QUERIES:
        count, ids = esearch("gds", term)
        lines.append(f"GEO/GDS\t{term}\t{count}\t{ids}")
    for term in SRA_QUERIES:
        count, ids = esearch("sra", term)
        lines.append(f"SRA\t{term}\t{count}\t{ids}")
    query = 'CLDN4 AND (knockdown OR knockout OR CRISPR OR shRNA OR siRNA)'
    count, ids = biostudies(query)
    lines.append(f"BioStudies\t{query}\t{count}\t{ids}")
    for db, query in NGDC_QUERIES:
        count, ids = ngdc(db, query)
        lines.append(f"NGDC/{db}\t{query}\t{count}\t{ids}")
    path = OUT / "query_counts.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
