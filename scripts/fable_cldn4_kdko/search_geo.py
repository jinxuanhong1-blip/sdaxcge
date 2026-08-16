"""Exhaustively search GEO DataSets (gds) for CLDN4 perturbation studies.

Strategy: multiple broad queries combining CLDN4/Claudin-4 with
knockdown/knockout/shRNA/siRNA/CRISPR/silencing terms, restricted to
GSE series (Entry Type = GSE). Human + mouse (and any organism, filtered
later). Aggregate unique GEO accessions and dump summaries for manual
verification.
"""
import json
import re
from collections import OrderedDict

from eutils import esearch, esummary

# Perturbation vocabulary.
PERTURB = [
    "knockdown", "knock-down", "knockout", "knock-out",
    "shRNA", "siRNA", "CRISPR", "sgRNA", "silencing", "silenced",
    "depletion", "depleted", "deficient", "loss of function",
    "gene editing", "ablation",
]

GENE = ['CLDN4', 'claudin-4', 'claudin 4', '"claudin 4"']

QUERIES = []
# Broad: CLDN4 + any perturbation term, GSE only
pert_or = " OR ".join(f'"{p}"' for p in PERTURB)
for g in GENE:
    QUERIES.append(f'({g}[All Fields]) AND ({pert_or}) AND "gse"[Entry Type]')
# Also a very broad net: CLDN4 anywhere in title/abstract as expression profiling
QUERIES.append('(CLDN4[Title] OR "claudin-4"[Title] OR "claudin 4"[Title]) AND "gse"[Entry Type]')
QUERIES.append('(CLDN4[Description] OR "claudin-4"[Description]) AND ("knockdown" OR "knockout" OR "shRNA" OR "siRNA" OR "CRISPR" OR "overexpression") AND "gse"[Entry Type]')

def main():
    all_ids = OrderedDict()
    for q in QUERIES:
        r = esearch("gds", q, retmax=1000)
        ids = r.get("idlist", [])
        print(f"[{r.get('count')}] {q}")
        for i in ids:
            all_ids.setdefault(i, q)
    print(f"\nTotal unique gds UIDs: {len(all_ids)}")

    summ = esummary("gds", list(all_ids.keys()))
    rows = []
    for uid, meta in summ.items():
        acc = meta.get("accession", "")
        rows.append({
            "uid": uid,
            "accession": acc,
            "title": meta.get("title", ""),
            "summary": meta.get("summary", ""),
            "taxon": meta.get("taxon", ""),
            "gdsType": meta.get("gdsType", ""),
            "n_samples": meta.get("n_samples", ""),
            "entryType": meta.get("entryType", ""),
            "pdat": meta.get("PDAT", ""),
            "found_by": all_ids.get(uid, ""),
        })
    rows.sort(key=lambda x: x["accession"])
    with open("../../results/fable_cldn4_kdko/geo_candidates_raw.json", "w") as f:
        json.dump(rows, f, indent=1)

    # Print a compact table for manual triage.
    print(f"\n{'ACC':14} {'TAXON':22} {'N':4} TITLE")
    for r in rows:
        if r["entryType"] != "GSE":
            continue
        print(f"{r['accession']:14} {r['taxon'][:22]:22} {str(r['n_samples']):4} {r['title'][:80]}")

if __name__ == "__main__":
    main()
