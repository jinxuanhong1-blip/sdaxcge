#!/usr/bin/env python3
"""Discover public ICI irAE (immune-related adverse event) transcriptomic datasets on GEO.

Queries NCBI E-utilities (esearch/esummary) on the GEO DataSets (gds) database for
immune checkpoint inhibitor (ICI) associated adverse events, with emphasis on lung
cancer contexts and irAEs such as pneumonitis and colitis.

Outputs a candidate table under results/fable_irae/tables/geo_candidates.tsv
and raw esummary json under results/fable_irae/data/.
No large downloads here; this only fetches metadata.
"""
import json
import time
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = Path(__file__).resolve().parents[2] / "results" / "fable_irae"
DATA_DIR = OUT_DIR / "data"
TAB_DIR = OUT_DIR / "tables"
for d in (DATA_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Query terms crafted around ICI + irAE + (lung / pneumonitis / colitis).
QUERIES = {
    "ici_irae": '("immune-related adverse event"[All Fields] OR "immune related adverse"[All Fields] OR irAE[All Fields]) AND "Homo sapiens"[Organism]',
    "ici_pneumonitis": '(pneumonitis[All Fields] AND (immunotherapy[All Fields] OR "checkpoint"[All Fields] OR PD-1[All Fields] OR PD-L1[All Fields] OR nivolumab[All Fields] OR pembrolizumab[All Fields])) AND "Homo sapiens"[Organism]',
    "ici_colitis": '(colitis[All Fields] AND (immunotherapy[All Fields] OR "checkpoint"[All Fields] OR PD-1[All Fields] OR CTLA-4[All Fields] OR ipilimumab[All Fields] OR nivolumab[All Fields])) AND "Homo sapiens"[Organism]',
    "lung_ici": '((NSCLC[All Fields] OR "non-small cell lung"[All Fields] OR "lung cancer"[All Fields]) AND ("checkpoint"[All Fields] OR PD-1[All Fields] OR PD-L1[All Fields] OR nivolumab[All Fields] OR pembrolizumab[All Fields] OR atezolizumab[All Fields]) AND (toxicity[All Fields] OR "adverse event"[All Fields] OR irAE[All Fields] OR pneumonitis[All Fields])) AND "Homo sapiens"[Organism]',
}


def eutils_get(endpoint, params):
    params = dict(params)
    params.setdefault("retmode", "json")
    url = f"{BASE}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa
            wait = 4 * (2 ** attempt)
            print(f"  retry {attempt+1} after {wait}s ({e})", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"failed: {url}")


def main():
    all_ids = {}
    for name, term in QUERIES.items():
        print(f"[esearch] {name}")
        res = eutils_get("esearch", {"db": "gds", "term": term, "retmax": 200})
        idlist = res.get("esearchresult", {}).get("idlist", [])
        print(f"  {len(idlist)} hits")
        all_ids[name] = idlist
        time.sleep(0.4)

    (DATA_DIR / "esearch_idlists.json").write_text(json.dumps(all_ids, indent=2))

    # Union of ids, fetch summaries.
    uid_set = sorted({uid for ids in all_ids.values() for uid in ids})
    print(f"[esummary] {len(uid_set)} unique UIDs")
    summaries = {}
    for i in range(0, len(uid_set), 100):
        chunk = uid_set[i:i + 100]
        res = eutils_get("esummary", {"db": "gds", "id": ",".join(chunk)})
        result = res.get("result", {})
        for uid in result.get("uids", []):
            summaries[uid] = result[uid]
        time.sleep(0.4)

    (DATA_DIR / "esummary_gds.json").write_text(json.dumps(summaries, indent=2))

    # Build candidate table.
    rows = []
    which_query = {}
    for name, ids in all_ids.items():
        for uid in ids:
            which_query.setdefault(uid, []).append(name)

    for uid, s in summaries.items():
        rows.append({
            "uid": uid,
            "accession": s.get("accession", ""),
            "gdstype": s.get("gdsType", ""),
            "entrytype": s.get("entryType", ""),
            "n_samples": s.get("n_samples", ""),
            "platform": s.get("gpl", ""),
            "taxon": s.get("taxon", ""),
            "pdat": s.get("PDAT", ""),
            "matched_queries": ";".join(which_query.get(uid, [])),
            "title": (s.get("title", "") or "").replace("\t", " ").replace("\n", " "),
            "summary": (s.get("summary", "") or "").replace("\t", " ").replace("\n", " ")[:400],
        })

    # Sort: prefer Series (GSE), with samples, recent.
    def sortkey(r):
        is_series = 1 if str(r["entrytype"]).upper() == "GSE" else 0
        try:
            ns = int(r["n_samples"])
        except Exception:
            ns = 0
        return (is_series, ns)

    rows.sort(key=sortkey, reverse=True)

    import csv
    out = TAB_DIR / "geo_candidates.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[write] {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
