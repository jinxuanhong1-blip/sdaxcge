#!/usr/bin/env python3
"""Harvest candidate literature for claim C3 from Europe PMC.

C3: "PD-L1/CD274 rises after CLDN4 loss, or after TROP2 ADC exposure."

The harvest is deliberately over-inclusive. Screening/scoring happens in
02_screen_literature.py; nothing here is treated as evidence for the claim.

Outputs
-------
logs/europepmc_raw/<qid>.json   raw API payloads, one file per query
catalog/records.tsv             deduplicated union of all hits
catalog/query_log.tsv           per-query hit counts and retrieval provenance
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "logs", "europepmc_raw")
CATALOG = os.path.join(ROOT, "catalog")

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = "C3-PDL1-catalog/1.0 (research literature harvest; contact: repo maintainer)"
PAGE = 100
MAX_PER_QUERY = 400

# arm A  = CLDN4 loss -> CD274
# arm B  = TROP2 ADC -> CD274
# arm M  = mechanistic bridge / prior that either arm would have to run through
# arm X  = adjacent claims that are easy to mistake for C3 (guards against
#          confirmation by lookalike)
QUERIES = [
    # ---------------- arm A: CLDN4 loss ----------------
    ("A01", "A", 'TITLE_ABS:"CLDN4" AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274")'),
    ("A02", "A", 'TITLE_ABS:"claudin-4" AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274")'),
    ("A03", "A", '"CLDN4" AND ("PD-L1" OR "CD274")'),
    ("A04", "A", 'TITLE_ABS:("CLDN4" OR "claudin-4" OR "claudin 4") AND TITLE_ABS:("knockdown" OR "knockout" OR "silencing" OR "CRISPR" OR "siRNA" OR "shRNA" OR "depletion" OR "loss")'),
    ("A05", "A", 'TITLE_ABS:("CLDN4" OR "claudin-4") AND TITLE_ABS:("immune evasion" OR "immune escape" OR "immune checkpoint" OR "tumor microenvironment" OR "T cell" OR "T-cell")'),
    ("A06", "A", 'TITLE_ABS:("claudin") AND TITLE_ABS:("PD-L1" OR "CD274") AND TITLE_ABS:("knockdown" OR "knockout" OR "silencing" OR "CRISPR" OR "siRNA" OR "shRNA" OR "overexpression")'),
    ("A07", "A", 'TITLE_ABS:("tight junction") AND TITLE_ABS:("PD-L1" OR "CD274")'),
    ("A08", "A", 'TITLE_ABS:("CLDN4" OR "claudin-4") AND TITLE_ABS:("interferon" OR "IFN" OR "STAT1" OR "NF-kB" OR "NF-kappaB")'),

    # ---------------- arm B: TROP2 ADC ----------------
    ("B01", "B", 'TITLE_ABS:("TROP2" OR "TROP-2" OR "TACSTD2") AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274")'),
    ("B02", "B", 'TITLE_ABS:("sacituzumab") AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274" OR "PD-1" OR "immune")'),
    ("B03", "B", 'TITLE_ABS:("datopotamab" OR "Dato-DXd" OR "DS-1062") AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274" OR "PD-1" OR "immune")'),
    ("B04", "B", 'TITLE_ABS:("SN-38" OR "SN38") AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274")'),
    ("B05", "B", 'TITLE_ABS:("deruxtecan" OR "DXd") AND TITLE_ABS:("PD-L1" OR "PDL1" OR "CD274")'),
    ("B06", "B", 'TITLE_ABS:("antibody-drug conjugate" OR "antibody drug conjugate" OR "ADC") AND TITLE_ABS:("PD-L1" OR "CD274") AND TITLE_ABS:("upregulat*" OR "up-regulat*" OR "induc*" OR "increase*")'),
    ("B07", "B", 'TITLE_ABS:("topoisomerase I" OR "topoisomerase 1" OR "irinotecan" OR "topotecan" OR "camptothecin") AND TITLE_ABS:("PD-L1" OR "CD274")'),
    ("B08", "B", 'TITLE_ABS:("sacituzumab tirumotecan" OR "SKB264" OR "MK-2870") AND TITLE_ABS:("PD-L1" OR "immune" OR "pembrolizumab")'),
    ("B09", "B", '("sacituzumab govitecan" OR "datopotamab deruxtecan") AND ("PD-L1" OR "CD274")'),
    ("B10", "B", 'TITLE_ABS:("TROP2" OR "TACSTD2") AND TITLE_ABS:("immunogenic cell death" OR "STING" OR "cGAS" OR "interferon")'),
    ("B11", "B", 'TITLE_ABS:("antibody-drug conjugate" OR "ADC") AND TITLE_ABS:("immunogenic cell death") AND TITLE_ABS:("PD-1" OR "PD-L1")'),

    # ---------------- arm M: mechanism the claim must route through ----------
    ("M01", "M", 'TITLE_ABS:("DNA damage") AND TITLE_ABS:("PD-L1" OR "CD274") AND TITLE_ABS:("STING" OR "cGAS" OR "interferon")'),
    ("M02", "M", 'TITLE_ABS:("CD274") AND TITLE_ABS:("transcriptional regulation" OR "promoter" OR "3\'UTR" OR "3\' UTR" OR "mRNA stability")'),
    ("M03", "M", 'TITLE_ABS:("EMT" OR "epithelial-mesenchymal transition") AND TITLE_ABS:("PD-L1" OR "CD274") AND TITLE_ABS:("claudin" OR "E-cadherin" OR "tight junction")'),

    # ---------------- arm X: lookalike claims to keep separate ---------------
    ("X01", "X", 'TITLE_ABS:("claudin-low") AND TITLE_ABS:("PD-L1" OR "CD274" OR "immune")'),
    ("X02", "X", 'TITLE_ABS:("CLDN18.2" OR "claudin-18.2" OR "CLDN18") AND TITLE_ABS:("PD-L1" OR "CD274")'),
    ("X03", "X", 'TITLE_ABS:("CLDN4" OR "claudin-4") AND TITLE_ABS:("prognosis" OR "prognostic" OR "survival")'),
]

FIELDS = [
    "id", "source", "pmid", "pmcid", "doi", "title", "authorString",
    "pubYear", "pubType", "isOpenAccess", "citedByCount", "abstractText",
]


def fetch(query: str, max_results: int = MAX_PER_QUERY):
    """Page through Europe PMC with cursorMark. Returns (hit_count, results)."""
    out, cursor, hit_count = [], "*", None
    while True:
        params = {
            "query": query,
            "format": "json",
            "pageSize": str(PAGE),
            "resultType": "core",
            "cursorMark": cursor,
        }
        url = f"{EPMC}?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=90) as fh:
                    payload = json.load(fh)
                break
            except Exception as exc:  # transient API/network failure
                if attempt == 3:
                    raise
                sys.stderr.write(f"  retry {attempt + 1} after {exc}\n")
                time.sleep(2 ** attempt * 2)
        if hit_count is None:
            hit_count = payload.get("hitCount", 0)
        batch = payload.get("resultList", {}).get("result", [])
        out.extend(batch)
        next_cursor = payload.get("nextCursorMark")
        if not batch or not next_cursor or next_cursor == cursor or len(out) >= max_results:
            break
        cursor = next_cursor
        time.sleep(0.34)  # stay well under the EPMC courtesy rate
    return hit_count, out[:max_results]


def norm_key(rec: dict) -> str:
    """Stable dedup key: prefer PMID, then DOI, then EPMC id."""
    if rec.get("pmid"):
        return "PMID:" + str(rec["pmid"])
    if rec.get("doi"):
        return "DOI:" + str(rec["doi"]).lower()
    return f"{rec.get('source', '?')}:{rec.get('id', '?')}"


def clean(value) -> str:
    if value is None:
        return ""
    return " ".join(str(value).replace("\t", " ").split())


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(CATALOG, exist_ok=True)

    merged: "OrderedDict[str, dict]" = OrderedDict()
    query_rows = []

    for qid, arm, query in QUERIES:
        sys.stderr.write(f"[{qid}] {query}\n")
        hit_count, results = fetch(query)
        sys.stderr.write(f"  hitCount={hit_count} retrieved={len(results)}\n")
        with open(os.path.join(RAW, f"{qid}.json"), "w") as fh:
            json.dump(
                {"qid": qid, "arm": arm, "query": query,
                 "hit_count": hit_count, "retrieved": len(results),
                 "results": results},
                fh, indent=1,
            )
        query_rows.append((qid, arm, query, hit_count, len(results)))

        for rec in results:
            key = norm_key(rec)
            if key not in merged:
                slim = {f: rec.get(f) for f in FIELDS}
                slim["_key"] = key
                slim["_queries"] = []
                slim["_arms"] = set()
                merged[key] = slim
            merged[key]["_queries"].append(qid)
            merged[key]["_arms"].add(arm)

    cols = ["rec_key", "pmid", "pmcid", "doi", "year", "open_access",
            "cited_by", "pub_type", "arms", "queries", "n_queries",
            "title", "journal_authors", "abstract"]
    path = os.path.join(CATALOG, "records.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for rec in merged.values():
            fh.write("\t".join([
                rec["_key"],
                clean(rec.get("pmid")),
                clean(rec.get("pmcid")),
                clean(rec.get("doi")),
                clean(rec.get("pubYear")),
                clean(rec.get("isOpenAccess")),
                clean(rec.get("citedByCount")),
                clean(rec.get("pubType")),
                "|".join(sorted(rec["_arms"])),
                "|".join(rec["_queries"]),
                str(len(rec["_queries"])),
                clean(rec.get("title")),
                clean(rec.get("authorString"))[:200],
                clean(rec.get("abstractText")),
            ]) + "\n")

    with open(os.path.join(CATALOG, "query_log.tsv"), "w") as fh:
        fh.write("qid\tarm\tquery\thit_count\tretrieved\n")
        for row in query_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    sys.stderr.write(f"\nunique records: {len(merged)} -> {path}\n")


if __name__ == "__main__":
    main()
