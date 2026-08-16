#!/usr/bin/env python3
"""Search GEO (NCBI `gds`) for human lung immunotherapy series in the 2019 slice.

Scope
-----
"2019 slice" = GEO series whose GEO release date (PDAT) falls in calendar 2019.
Submission date is not an indexed Entrez field, so it is recovered later from
the SOFT record (02_fetch_metadata.py) and reported alongside PDAT.

This is deliberately a *wider* net than the earlier 2019-2021 sweep
(`results/fable_geo_2019_2021/`): more drug names (including the China-approved
PD-1 antibodies), more treatment phrasings, and target-gene-driven queries.
The point is to surface 2019 series that a narrower lung x checkpoint query
misses, not to re-derive the same list.

Every accession written out is returned by NCBI; none are invented.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "w200" / "GEO_2019" / "search"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATE_FILTER = '("2019/01/01"[PDAT] : "2019/12/31"[PDAT])'
COMMON = '"Homo sapiens"[Organism] AND "gse"[Entry Type]'

LUNG_TERMS = [
    "lung", "NSCLC", "non-small cell lung", "non small cell lung",
    "lung adenocarcinoma", "LUAD", "lung squamous", "LUSC",
    "small cell lung", "SCLC", "pulmonary carcinoma", "pleural",
    "bronchial carcinoma", "lung carcinoma", "lung cancer",
]

ICI_TERMS = [
    "immune checkpoint", "checkpoint inhibitor", "checkpoint blockade",
    "immunotherapy", "immuno-oncology", "PD-1", "PD1", "PD-L1", "PDL1",
    "PDCD1", "CD274", "CTLA-4", "CTLA4", "anti-PD-1", "anti-PD-L1",
    "nivolumab", "pembrolizumab", "atezolizumab", "durvalumab", "avelumab",
    "ipilimumab", "cemiplimab", "tremelimumab", "camrelizumab", "sintilimab",
    "tislelizumab", "toripalimab", "immune-related adverse", "ICB", "ICI",
    "neoadjuvant immunotherapy", "chemoimmunotherapy", "adjuvant immunotherapy",
]

# Target-gene driven queries: any 2019 human lung series that explicitly
# mentions TROP2/TACSTD2 or CLDN4, even without an ICI keyword, is worth a look
# because the treatment context is often only in the sample metadata.
TARGET_TERMS = [
    "TACSTD2", "TROP2", "TROP-2", "CLDN4", "claudin-4", "claudin 4",
    "sacituzumab", "datopotamab", "Dato-DXd", "SKB264", "sacituzumab govitecan",
]


def _or(terms):
    return "(" + " OR ".join(f'"{t}"[All Fields]' for t in terms) + ")"


def esearch(term, retmax=2000):
    params = {"db": "gds", "term": term, "retmax": str(retmax), "retmode": "json"}
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r).get("esearchresult", {})
        except Exception as exc:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  retry {attempt} in {wait}s ({exc})", file=sys.stderr)
            time.sleep(wait)
    return {}


def main():
    lung = _or(LUNG_TERMS)
    ici = _or(ICI_TERMS)
    queries = {}

    queries["lung_x_ici"] = f"{lung} AND {ici} AND {COMMON} AND {DATE_FILTER}"
    queries["lung_x_target"] = f"{lung} AND {_or(TARGET_TERMS)} AND {COMMON} AND {DATE_FILTER}"
    # Treatment-context series are sometimes titled by trial/《response》wording only.
    queries["lung_x_response"] = (
        f"{lung} AND "
        '("responder"[All Fields] OR "non-responder"[All Fields] OR '
        '"durable clinical benefit"[All Fields] OR "progression-free survival"[All Fields] OR '
        '"treatment response"[All Fields] OR "resistance to therapy"[All Fields]) '
        f"AND {COMMON} AND {DATE_FILTER}"
    )
    for drug in ["nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
                 "ipilimumab", "avelumab", "cemiplimab", "tremelimumab",
                 "camrelizumab", "sintilimab", "tislelizumab", "toripalimab"]:
        queries[f"drug_{drug}"] = (
            f'{lung} AND "{drug}"[All Fields] AND {COMMON} AND {DATE_FILTER}'
        )

    all_ids, per_query = set(), {}
    for name, term in queries.items():
        res = esearch(term)
        ids = set(res.get("idlist", []))
        new = ids - all_ids
        per_query[name] = {"term": term, "count": res.get("count"), "n_ids": len(ids),
                           "n_new": len(new), "uids": sorted(ids)}
        all_ids |= ids
        print(f"{name:26s} count={res.get('count'):>5} ids={len(ids):>4} new={len(new):>4}")
        time.sleep(0.4)

    out = {"date_filter": DATE_FILTER, "queries": per_query, "uids": sorted(all_ids)}
    (OUT_DIR / "search_uids.json").write_text(json.dumps(out, indent=2))
    print(f"\nTotal unique gds UIDs in 2019 slice: {len(all_ids)}")
    print("Wrote", OUT_DIR / "search_uids.json")


if __name__ == "__main__":
    main()
