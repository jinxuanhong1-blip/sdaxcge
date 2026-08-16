#!/usr/bin/env python3
"""GEO leftover search: 2023–2026 public lung neoadjuvant/ICI scRNA with MPR in title or metadata."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

MUST_TRY = [
    "GSE207422",
    "GSE241934",
    "GSE291670",
    "GSE205335",
    "GSE243013",
]

QUERIES = [
    (
        "mpr_scrna_lung_2023_2026",
        '(((lung OR NSCLC OR "non-small cell") AND (neoadjuvant OR "PD-1" OR "PD-L1" OR immunotherapy OR ICI) AND (MPR OR "major pathologic" OR "major pathological") AND (scRNA-seq OR "single-cell" OR "single cell")) AND ("2023/01/01"[PDAT] : "2026/12/31"[PDAT])) AND "gse"[Entry Type]',
    ),
    (
        "recist_scrna_lung_2023_2026",
        '(((lung OR NSCLC) AND (immunotherapy OR "PD-1" OR ICI) AND (RECIST OR responder OR "non-responder") AND (scRNA-seq OR "single-cell")) AND ("2023/01/01"[PDAT] : "2026/12/31"[PDAT])) AND "gse"[Entry Type]',
    ),
    (
        "broader_neoadj_scrna_2023_2026",
        '((NSCLC OR "lung cancer") AND neoadjuvant AND (scRNA-seq OR "single-cell RNA") AND ("2023/01/01"[PDAT] : "2026/12/31"[PDAT])) AND "gse"[Entry Type]',
    ),
]


def esearch(term: str, retmax: int = 80) -> list[str]:
    q = urllib.parse.urlencode(
        {
            "db": "gds",
            "term": term,
            "retmax": str(retmax),
            "retmode": "json",
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}"
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.loads(r.read().decode())
    return data.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    out = []
    for i in range(0, len(ids), 20):
        chunk = ids[i : i + 20]
        q = urllib.parse.urlencode(
            {
                "db": "gds",
                "id": ",".join(chunk),
                "retmode": "json",
            }
        )
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{q}"
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.loads(r.read().decode())
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            out.append(
                {
                    "uid": uid,
                    "accession": rec.get("accession"),
                    "title": rec.get("title"),
                    "summary": rec.get("summary"),
                    "n_samples": rec.get("n_samples"),
                    "taxon": rec.get("taxon"),
                    "gdstype": rec.get("gdstype"),
                    "pdat": rec.get("pdat"),
                    "ftp": rec.get("ftplink"),
                }
            )
        time.sleep(0.35)
    return out


def main() -> None:
    all_hits: dict[str, dict] = {}
    by_query = {}
    for name, term in QUERIES:
        print(f"search {name}", flush=True)
        try:
            ids = esearch(term)
        except Exception as exc:
            print(f"  FAIL {exc}", flush=True)
            by_query[name] = {"error": str(exc), "ids": []}
            continue
        print(f"  {len(ids)} ids", flush=True)
        recs = esummary(ids)
        by_query[name] = {"n": len(recs), "accessions": [r.get("accession") for r in recs]}
        for rec in recs:
            acc = rec.get("accession") or rec.get("uid")
            all_hits[acc] = rec
        time.sleep(0.35)

    rows = []
    for acc, rec in sorted(all_hits.items(), key=lambda x: str(x[0])):
        title = (rec.get("title") or "") + " " + (rec.get("summary") or "")
        low = title.lower()
        rows.append(
            {
                **rec,
                "must_try": acc in MUST_TRY,
                "mentions_mpr": ("mpr" in low) or ("major pathologic" in low) or ("major pathological" in low),
                "mentions_recist": "recist" in low,
                "mentions_scrna": ("single-cell" in low) or ("single cell" in low) or ("scrna" in low),
                "mentions_lung": ("lung" in low) or ("nsclc" in low),
            }
        )

    (OUT / "geo_leftover_search.json").write_text(json.dumps({"queries": by_query, "hits": rows}, indent=2))
    # TSV
    keys = [
        "accession",
        "pdat",
        "n_samples",
        "gdstype",
        "must_try",
        "mentions_mpr",
        "mentions_recist",
        "mentions_scrna",
        "mentions_lung",
        "title",
    ]
    lines = ["\t".join(keys)]
    for r in rows:
        lines.append("\t".join(str(r.get(k, "")).replace("\t", " ")[:300] for k in keys))
    (OUT / "geo_leftover_search.tsv").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(rows)} unique hits", flush=True)


if __name__ == "__main__":
    main()
