#!/usr/bin/env python3
"""Search GEO (NCBI gds) for human lung ICI series published 2017-2018.

This is the 2017-2018 slice of the ongoing TACSTD2/CLDN4 x lung-ICI GEO
sweep (earlier slices covered 2019-2021, 2022-2023, and a remaining sweep).

Strategy: Entrez esearch on `gds`, restricted to GSE entries, Homo sapiens,
PDAT 2017/01/01 - 2018/12/31, combining lung-cancer terms with
immune-checkpoint-inhibitor terms; per-drug queries added for exhaustiveness.
All accessions come straight from NCBI; none are invented. esummary metadata
for every hit is stored for triage.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "results" / "w200" / "GEO_2017_2018"
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

DATE_FILTER = '("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
COMMON = '"Homo sapiens"[Organism] AND "gse"[Entry Type]'


def _or(terms):
    return "(" + " OR ".join(f'"{t}"[All Fields]' for t in terms) + ")"


def _get_json(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  retry {attempt} after {wait}s ({e})", file=sys.stderr)
            time.sleep(wait)
    return {}


def esearch(term, retmax=2000):
    params = {"db": "gds", "term": term, "retmax": str(retmax), "retmode": "json"}
    return _get_json(f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)).get(
        "esearchresult", {})


def esummary(uids):
    params = {"db": "gds", "id": ",".join(uids), "retmode": "json"}
    return _get_json(f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params))


def main():
    lung = _or(LUNG_TERMS)
    ici = _or(ICI_TERMS)
    term = f"{lung} AND {ici} AND {COMMON} AND {DATE_FILTER}"
    print("Primary combined query:\n", term)
    res = esearch(term)
    all_ids = set(res.get("idlist", []))
    print(f"Combined query: count={res.get('count')} ids={len(all_ids)}")

    for drug in ["nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
                 "ipilimumab", "avelumab", "cemiplimab"]:
        t = f'{lung} AND "{drug}"[All Fields] AND {COMMON} AND {DATE_FILTER}'
        r = esearch(t)
        new = set(r.get("idlist", []))
        print(f"  drug={drug}: count={r.get('count')} new={len(new - all_ids)}")
        all_ids |= new
        time.sleep(0.4)

    print(f"Total unique gds UIDs: {len(all_ids)}")
    (OUT_DIR / "search_uids.json").write_text(
        json.dumps({"query": term, "uids": sorted(all_ids)}, indent=2))

    uids = sorted(all_ids)
    records = []
    for i in range(0, len(uids), 20):
        result = esummary(uids[i:i + 20]).get("result", {})
        for uid in result.get("uids", []):
            r = result[uid]
            records.append({
                "uid": uid,
                "accession": r.get("accession"),
                "title": r.get("title"),
                "summary": r.get("summary"),
                "taxon": r.get("taxon"),
                "n_samples": r.get("n_samples"),
                "gdsType": r.get("gdsType"),
                "gpl": r.get("gpl"),
                "pdat": r.get("pdat"),
                "suppFile": r.get("suppFile"),
                "ftplink": r.get("ftplink"),
                "pubmedids": r.get("pubmedids"),
            })
        print(f"esummary {min(i + 20, len(uids))}/{len(uids)}")
        time.sleep(0.4)

    records.sort(key=lambda x: x.get("accession") or "")
    (OUT_DIR / "candidates_metadata.json").write_text(json.dumps(records, indent=2))
    print(f"Wrote {len(records)} records to candidates_metadata.json")


if __name__ == "__main__":
    main()
