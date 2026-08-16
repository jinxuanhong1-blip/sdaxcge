#!/usr/bin/env python3
"""Search GEO (db=gds) for human lung-cancer ICI series published 2024-2026.

Uses NCBI E-utilities (esearch + esummary). Excludes a fixed list of GSEs
already analyzed in prior slices. Writes candidate table to
results/fable_geo_2026/geo_search_candidates.tsv and the raw esummary JSON
to results/fable_geo_2026/geo_search_raw.json.

No accessions are invented: everything comes directly from the live
E-utilities responses.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2026"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE = {
    "GSE126044", "GSE135222", "GSE136961", "GSE166449",
    "GSE93157", "GSE207422", "GSE205335",
}

LUNG_TERMS = [
    '"lung cancer"', '"lung adenocarcinoma"', '"lung squamous"', "NSCLC",
    '"non-small cell lung"', '"non small cell lung"', "SCLC",
    '"small cell lung"', "LUAD", "LUSC",
]
ICI_TERMS = [
    '"immune checkpoint"', '"checkpoint inhibitor"', '"checkpoint blockade"',
    "immunotherapy", '"anti-PD-1"', '"anti-PD-L1"', '"anti-PD1"', '"anti-PDL1"',
    "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "cemiplimab",
    "ipilimumab", "tislelizumab", "sintilimab", "camrelizumab", "toripalimab",
    "sugemalimab", "serplulimab", "PD-1[Title]", "PD-L1[Title]",
]

query = (
    "(" + " OR ".join(LUNG_TERMS) + ") AND (" + " OR ".join(ICI_TERMS) + ")"
    ' AND "Homo sapiens"[Organism] AND gse[Entry Type]'
    ' AND ("2024/01/01"[PDAT] : "2026/12/31"[PDAT])'
)


def eutils(endpoint: str, params: dict) -> dict:
    url = f"{BASE}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa: BLE001 - retry on transient network/API errors
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
            print(f"  retry {attempt + 1} after {e}")
    raise RuntimeError("unreachable")


def main() -> None:
    print("Query:", query)
    es = eutils("esearch", {
        "db": "gds", "term": query, "retmax": 5000, "retmode": "json",
        "usehistory": "y",
    })
    ids = es["esearchresult"]["idlist"]
    print(f"esearch hits: {es['esearchresult']['count']} (fetched {len(ids)})")

    docs = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        summ = eutils("esummary", {
            "db": "gds", "id": ",".join(chunk), "retmode": "json",
        })
        for uid in summ["result"]["uids"]:
            docs[uid] = summ["result"][uid]
        print(f"  esummary {i + len(chunk)}/{len(ids)}")
        time.sleep(0.4)

    (OUT_DIR / "geo_search_raw.json").write_text(json.dumps(docs, indent=1))

    rows = []
    for uid, d in docs.items():
        acc = d.get("accession", "")
        if not acc.startswith("GSE") or acc in EXCLUDE:
            continue
        rows.append({
            "accession": acc,
            "title": d.get("title", "").replace("\t", " "),
            "gdstype": d.get("gdstype", ""),
            "n_samples": d.get("n_samples", ""),
            "pdat": d.get("pdat", ""),
            "suppfile": d.get("suppfile", ""),
            "ftplink": d.get("ftplink", ""),
            "summary": d.get("summary", "").replace("\t", " ").replace("\n", " "),
        })
    rows.sort(key=lambda r: r["accession"])

    out = OUT_DIR / "geo_search_candidates.tsv"
    cols = ["accession", "title", "gdstype", "n_samples", "pdat", "suppfile",
            "ftplink", "summary"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {len(rows)} candidate GSEs to {out}")


if __name__ == "__main__":
    main()
