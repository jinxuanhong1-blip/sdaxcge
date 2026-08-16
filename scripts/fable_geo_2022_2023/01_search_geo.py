#!/usr/bin/env python3
"""
Exhaustive GEO search for human lung-cancer ICI (immune checkpoint inhibitor)
series published in 2022-2023.

Strategy
--------
GEO is queried through NCBI E-utilities (db=gds). Because a single boolean
string is easy to make either too narrow or too broad, we run a *union* of
several complementary queries (different ICI vocabulary, different lung-cancer
vocabulary) and de-duplicate by GSE accession. Every query is restricted to:

  * Entry type    : GSE series           ("gse"[ETYP])
  * Organism      : Homo sapiens
  * Publication   : 2022/01/01 - 2023/12/31   ([PDAT])

The raw union of series accessions + their esummary metadata is written to
results/fable_geo_2022_2023/geo_search_candidates.csv so the (heavier)
verification stage can run offline against a stable snapshot.
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUTDIR = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2022_2023"
NOTEDIR = Path(__file__).resolve().parents[2] / "notes" / "fable_geo_2022_2023"
OUTDIR.mkdir(parents=True, exist_ok=True)
NOTEDIR.mkdir(parents=True, exist_ok=True)

TOOL = "fable_geo_2022_2023"
EMAIL = "cursoragent@cursor.com"

DATE_FILTER = '("2022/01/01"[PDAT] : "2023/12/31"[PDAT])'
ORG = 'Homo sapiens[Organism]'
ETYP = 'gse[ETYP]'

# ICI vocabulary (immune checkpoint inhibitors / immunotherapy)
ICI_TERMS = [
    "immune checkpoint",
    "checkpoint inhibitor",
    "checkpoint blockade",
    "immunotherapy",
    "anti-PD-1", "anti-PD1", "PD-1", "PD1",
    "anti-PD-L1", "anti-PDL1", "PD-L1", "PDL1",
    "CTLA-4", "CTLA4",
    "nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
    "cemiplimab", "avelumab", "ipilimumab", "tislelizumab",
    "sintilimab", "camrelizumab", "toripalimab",
]

# Lung-cancer vocabulary
LUNG_TERMS = [
    "lung",
    "NSCLC", "non-small cell lung", "non small cell lung",
    "LUAD", "lung adenocarcinoma",
    "LUSC", "lung squamous",
    "SCLC", "small cell lung",
    "pulmonary",
]


def eutils_get(endpoint, params, retries=4):
    params = {**params, "tool": TOOL, "email": EMAIL}
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"eutils failed: {endpoint} {params}: {last}")


def esearch_gds(term, retmax=2000):
    txt = eutils_get("esearch.fcgi", {
        "db": "gds", "term": term, "retmax": retmax, "retmode": "json",
    })
    data = json.loads(txt)
    return data["esearchresult"].get("idlist", [])


def build_terms():
    lung_block = "(" + " OR ".join(f'"{t}"' for t in LUNG_TERMS) + ")"
    ici_block = "(" + " OR ".join(f'"{t}"' for t in ICI_TERMS) + ")"
    base = f"{ETYP} AND {ORG} AND {DATE_FILTER}"
    # Query 1: broad union (most recall)
    q_union = f"{base} AND {lung_block} AND {ici_block}"
    queries = {"union_broad": q_union}
    # Per-ICI-drug queries to catch series where the vocabulary block
    # phrase parsing behaves unexpectedly (belt-and-suspenders for recall).
    for drug in ["nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
                 "immune checkpoint", "immunotherapy", "PD-1", "PD-L1"]:
        queries[f"ici_{drug.replace(' ', '_')}"] = (
            f'{base} AND {lung_block} AND "{drug}"')
    return queries


def esummary_gds(uids):
    """Fetch esummary metadata for a list of UIDs (batched)."""
    out = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        txt = eutils_get("esummary.fcgi", {
            "db": "gds", "id": ",".join(chunk), "retmode": "json",
        })
        data = json.loads(txt)
        res = data.get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
        time.sleep(0.34)
    return out


def main():
    queries = build_terms()
    uid_to_queries = {}
    for name, term in queries.items():
        ids = esearch_gds(term)
        print(f"[{name}] {len(ids)} hits", file=sys.stderr)
        for uid in ids:
            uid_to_queries.setdefault(uid, set()).add(name)
        time.sleep(0.34)

    all_uids = sorted(uid_to_queries)
    print(f"Total unique UIDs: {len(all_uids)}", file=sys.stderr)

    meta = esummary_gds(all_uids)

    rows = []
    for uid in all_uids:
        m = meta.get(uid, {})
        acc = m.get("accession", "")
        rows.append({
            "uid": uid,
            "accession": acc,
            "title": m.get("title", ""),
            "gpl": m.get("gpl", ""),
            "taxon": m.get("taxon", ""),
            "entrytype": m.get("entrytype", ""),
            "gdstype": m.get("gdstype", ""),
            "n_samples": m.get("n_samples", ""),
            "pdat": m.get("pdat", ""),
            "suppfile": m.get("suppfile", ""),
            "ftplink": m.get("ftplink", ""),
            "summary": (m.get("summary", "") or "").replace("\n", " "),
            "matched_queries": ";".join(sorted(uid_to_queries[uid])),
        })

    # Keep only GSE series entries
    rows = [r for r in rows if r["entrytype"] == "GSE" or r["accession"].startswith("GSE")]
    rows.sort(key=lambda r: r["accession"])

    out_csv = OUTDIR / "geo_search_candidates.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Also record the exact queries used for reproducibility
    (NOTEDIR / "search_queries.json").write_text(
        json.dumps(queries, indent=2, ensure_ascii=False))

    print(f"Wrote {len(rows)} candidate series -> {out_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
