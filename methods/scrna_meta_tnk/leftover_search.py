#!/usr/bin/env python3
"""GEO E-utilities search: 2023–2026 human lung scRNA series (public catalog only)."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / "results"
TERM = (
    '(lung[Title] OR NSCLC[Title] OR LUAD[Title] OR LUSC[Title] OR "non-small cell"[Title]) '
    'AND (scRNA-seq[All Fields] OR "single-cell"[Title] OR "single cell"[Title]) '
    'AND Homo sapiens[Organism] AND ("2023/01/01"[PDAT] : "2026/12/31"[PDAT]) AND gse[ETYP]'
)


def esearch(term: str, retmax: int = 300) -> list[str]:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q, timeout=60
    ) as r:
        js = json.load(r)
    return js.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str]) -> list[dict]:
    rows = []
    for i in range(0, len(ids), 20):
        batch = ids[i : i + 20]
        q = urllib.parse.urlencode({"db": "gds", "id": ",".join(batch), "retmode": "json"})
        with urllib.request.urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + q, timeout=90
        ) as r:
            sm = json.load(r)
        result = sm.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            rows.append(
                {
                    "accession": rec.get("accession"),
                    "pdat": rec.get("pdat"),
                    "n_samples": rec.get("n_samples"),
                    "gdstype": rec.get("gdstype"),
                    "title": rec.get("title"),
                    "summary": (rec.get("summary") or "")[:500],
                }
            )
        time.sleep(0.35)
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ids = esearch(TERM)
    rows = esummary(ids)
    (OUT / "geo_leftover_2023_2026.json").write_text(json.dumps(rows, indent=2))
    print("wrote", len(rows), "series")


if __name__ == "__main__":
    main()
