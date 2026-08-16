#!/usr/bin/env python3
"""Broader GEO GDS search (not title-restricted) for missed 2023-2026 series."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/hunt_neoadj/metadata")
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

QUERIES = {
    "gds_neoadj_lung_io_any_2023_2026": (
        '(lung OR NSCLC OR "non-small cell lung" OR LUAD OR LUSC) AND '
        '(neoadjuvant OR "neo-adjuvant" OR "induction immunotherapy") AND '
        '(immunotherapy OR "PD-1" OR "PD-L1" OR pembrolizumab OR nivolumab OR '
        'sintilimab OR camrelizumab OR tislelizumab OR toripalimab OR durvalumab OR '
        '"checkpoint") AND 2023:2026[PDAT] AND "gse"[ETYP]'
    ),
    "gds_scrna_lung_mpr_2023_2026": (
        '(lung OR NSCLC) AND (scRNA OR "single-cell" OR "single cell") AND '
        '(MPR OR "major pathologic" OR "pathological response" OR pCR) AND '
        '2023:2026[PDAT] AND "gse"[ETYP]'
    ),
    "gds_scrna_lung_neoadj_2023_2026": (
        '(lung OR NSCLC) AND (scRNA OR "single-cell RNA" OR "single cell RNA") AND '
        '(neoadjuvant OR "neo-adjuvant") AND 2023:2026[PDAT] AND "gse"[ETYP]'
    ),
}


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "hunt-neoadj/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def esearch(term: str) -> list[str]:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": 300, "retmode": "json"})
    return json.loads(get(f"{ESEARCH}?{q}"))["esearchresult"].get("idlist", [])


def esummary(ids: list[str]) -> dict:
    out = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i : i + 40]
        q = urllib.parse.urlencode({"db": "gds", "id": ",".join(chunk), "retmode": "json"})
        out.update(json.loads(get(f"{ESUMMARY}?{q}")).get("result", {}))
        time.sleep(0.35)
    return out


def main():
    all_ids = set()
    hits = {}
    for name, term in QUERIES.items():
        ids = esearch(term)
        hits[name] = {"n": len(ids), "ids": ids, "term": term}
        all_ids.update(ids)
        print(name, len(ids), flush=True)
        time.sleep(0.35)
    summaries = esummary(sorted(all_ids))
    records = []
    for uid, rec in summaries.items():
        if uid == "uids" or not isinstance(rec, dict):
            continue
        acc = rec.get("accession") or ""
        if not str(acc).startswith("GSE"):
            continue
        records.append(
            {
                "accession": acc,
                "pdat": rec.get("pdat"),
                "n_samples": rec.get("n_samples"),
                "gdstype": rec.get("gdstype"),
                "title": rec.get("title"),
                "summary": rec.get("summary"),
            }
        )
    records.sort(key=lambda r: (r.get("pdat") or "", r["accession"]))
    payload = {"queries": hits, "records": records}
    (OUT / "geo_search_broader.json").write_text(json.dumps(payload, indent=2))
    lines = ["accession\tpdat\tn_samples\tgdstype\ttitle"]
    for r in records:
        title = (r.get("title") or "").replace("\t", " ").replace("\n", " ")
        lines.append(f"{r['accession']}\t{r.get('pdat')}\t{r.get('n_samples')}\t{r.get('gdstype')}\t{title}")
    (OUT / "geo_search_broader.tsv").write_text("\n".join(lines) + "\n")
    print("GSE records", len(records))


if __name__ == "__main__":
    main()
