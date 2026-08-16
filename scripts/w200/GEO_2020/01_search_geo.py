#!/usr/bin/env python3
"""Exhaustive GEO search for human lung ICI series with PDAT in calendar 2020.

Same query family as the 2019–2021 / 2022–2023 waves, but the date window is
strictly 2020/01/01–2020/12/31. Every returned UID is a real NCBI gds id.
No accessions are invented.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[3] / "results" / "w200" / "GEO_2020"
OUT_DIR.mkdir(parents=True, exist_ok=True)

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
DATE_FILTER = '("2020/01/01"[PDAT] : "2020/12/31"[PDAT])'
COMMON = '"Homo sapiens"[Organism] AND "gse"[Entry Type]'


def _or(terms):
    return "(" + " OR ".join(f'"{t}"[All Fields]' for t in terms) + ")"


def esearch(term, retmax=2000):
    params = {"db": "gds", "term": term, "retmax": str(retmax), "retmode": "json"}
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r).get("esearchresult", {})
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
    res = esearch(term)
    ids = set(res.get("idlist", []))
    print(f"Combined query count={res.get('count')} fetched={len(ids)}")

    per_drug = {}
    for drug in ["nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
                 "ipilimumab", "avelumab", "cemiplimab"]:
        t = f'{lung} AND "{drug}"[All Fields] AND {COMMON} AND {DATE_FILTER}'
        r = esearch(t, retmax=1000)
        new = set(r.get("idlist", []))
        per_drug[drug] = {"count": r.get("count"), "n": len(new), "new": len(new - ids)}
        print(f"  drug={drug}: count={r.get('count')} ids={len(new)} new={len(new - ids)}")
        ids |= new
        time.sleep(0.35)

    payload = {
        "query": term,
        "date_window": "2020/01/01-2020/12/31",
        "organism": "Homo sapiens",
        "entry_type": "gse",
        "n_uids": len(ids),
        "per_drug": per_drug,
        "uids": sorted(ids),
    }
    out = OUT_DIR / "search_uids.json"
    out.write_text(json.dumps(payload, indent=2))
    print(f"Total unique gds UIDs: {len(ids)}")
    print("Wrote", out)


if __name__ == "__main__":
    main()
