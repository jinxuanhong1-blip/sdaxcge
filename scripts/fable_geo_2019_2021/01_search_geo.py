#!/usr/bin/env python3
"""Exhaustively search GEO (NCBI gds database) for human lung ICI series (2019-2021).

Strategy
--------
We issue several complementary Entrez esearch queries against the `gds`
database, each restricted to:
  * Entry type GSE (series, not samples/platforms/datasets)
  * Organism Homo sapiens
  * Publication date window 2019/01/01 - 2021/12/31 (PDAT)
and combining lung-cancer terms with immune-checkpoint-inhibitor terms.

All returned GEO accessions are REAL accessions reported by NCBI. We never
invent IDs. The union of hits is written to a candidates table for manual
verification in the next stage.
"""
import time
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2019_2021"
NOTES_DIR = Path(__file__).resolve().parents[2] / "notes" / "fable_geo_2019_2021"
OUT_DIR.mkdir(parents=True, exist_ok=True)
NOTES_DIR.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

LUNG_TERMS = [
    "lung", "NSCLC", "non-small cell lung", "non small cell lung",
    "lung adenocarcinoma", "LUAD", "lung squamous", "LUSC",
    "small cell lung", "SCLC", "pulmonary carcinoma",
]

ICI_TERMS = [
    "immune checkpoint", "checkpoint inhibitor", "checkpoint blockade",
    "immunotherapy", "PD-1", "PD1", "PD-L1", "PDL1", "CTLA-4", "CTLA4",
    "anti-PD-1", "anti-PD-L1", "nivolumab", "pembrolizumab", "atezolizumab",
    "durvalumab", "avelumab", "ipilimumab", "cemiplimab",
]

DATE_FILTER = '("2019/01/01"[PDAT] : "2021/12/31"[PDAT])'
COMMON = '"Homo sapiens"[Organism] AND "gse"[Entry Type]'


def _or(terms):
    return "(" + " OR ".join(f'"{t}"[All Fields]' for t in terms) + ")"


def esearch(term, retmax=1000):
    params = {
        "db": "gds",
        "term": term,
        "retmax": str(retmax),
        "retmode": "json",
    }
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = json.load(r)
            return data.get("esearchresult", {})
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  esearch retry {attempt} after {wait}s ({e})", file=sys.stderr)
            time.sleep(wait)
    return {}


def main():
    lung = _or(LUNG_TERMS)
    ici = _or(ICI_TERMS)
    term = f"{lung} AND {ici} AND {COMMON} AND {DATE_FILTER}"
    print("Primary combined query:\n", term)
    res = esearch(term, retmax=2000)
    ids = res.get("idlist", [])
    count = res.get("count")
    print(f"Combined query returned count={count}, ids fetched={len(ids)}")

    all_ids = set(ids)

    # Also run per-drug queries to catch series indexed only under a drug name,
    # to be as exhaustive as possible.
    for drug in ["nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
                 "ipilimumab", "avelumab", "cemiplimab"]:
        t = f'{lung} AND "{drug}"[All Fields] AND {COMMON} AND {DATE_FILTER}'
        r = esearch(t, retmax=1000)
        new = set(r.get("idlist", []))
        print(f"  drug={drug}: count={r.get('count')} ids={len(new)} new={len(new - all_ids)}")
        all_ids |= new
        time.sleep(0.4)

    print(f"\nTotal unique gds UIDs: {len(all_ids)}")
    (OUT_DIR / "search_uids.json").write_text(
        json.dumps({"query": term, "uids": sorted(all_ids)}, indent=2))
    print("Wrote", OUT_DIR / "search_uids.json")


if __name__ == "__main__":
    main()
