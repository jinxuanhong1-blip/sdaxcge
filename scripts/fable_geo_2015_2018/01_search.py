"""Exhaustive GEO search for 2015-2018 human lung-cancer immune-checkpoint
(ICI / PD-1 / PD-L1 / CTLA-4) SERIES (GSE).

Strategy: run several complementary queries against the GEO DataSets (gds)
database, union the resulting series, then fetch esummary metadata for each.
Filtering / verification is done in 02_verify.py so this script keeps the raw,
un-invented union of real accessions returned by NCBI.

Outputs:
  results/fable_geo_2015_2018/search_raw_summaries.json
  results/fable_geo_2015_2018/search_candidates.tsv
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eutils import esearch, esummary  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")
os.makedirs(OUT, exist_ok=True)

DATE = '("2015/01/01"[PDAT] : "2018/12/31"[PDAT])'
HUMAN = '"Homo sapiens"[Organism]'
GSE = 'gse[Entry Type]'

LUNG = ('(lung[Title] OR NSCLC[Title] OR "non-small cell lung"[All Fields] OR '
        '"lung cancer"[All Fields] OR "lung carcinoma"[All Fields] OR '
        '"lung adenocarcinoma"[All Fields] OR "lung neoplasms"[All Fields] OR '
        '"lung"[All Fields])')

ICI = ('("PD-1"[All Fields] OR "PD1"[All Fields] OR "PDCD1"[All Fields] OR '
       '"PD-L1"[All Fields] OR "PDL1"[All Fields] OR "CD274"[All Fields] OR '
       '"CTLA-4"[All Fields] OR "CTLA4"[All Fields] OR '
       '"immune checkpoint"[All Fields] OR "checkpoint inhibitor"[All Fields] OR '
       '"checkpoint blockade"[All Fields] OR "immunotherapy"[All Fields] OR '
       '"anti-PD-1"[All Fields] OR "anti-PD-L1"[All Fields] OR '
       '"nivolumab"[All Fields] OR "pembrolizumab"[All Fields] OR '
       '"atezolizumab"[All Fields] OR "durvalumab"[All Fields] OR '
       '"avelumab"[All Fields] OR "ipilimumab"[All Fields] OR '
       '"tremelimumab"[All Fields])')

QUERIES = {
    "lung_ici_human_gse_2015_2018": f"{LUNG} AND {ICI} AND {HUMAN} AND {GSE} AND {DATE}",
    # Broader: ICI drugs + lung, no strict organism token (catch mislabeled organism)
    "lung_ici_gse_2015_2018": f"{LUNG} AND {ICI} AND {GSE} AND {DATE}",
}


def main():
    all_uids = {}
    for name, term in QUERIES.items():
        uids = esearch("gds", term)
        print(f"[{name}] {len(uids)} UIDs", flush=True)
        all_uids[name] = uids

    union = sorted(set(u for lst in all_uids.values() for u in lst), key=int)
    print(f"union UIDs: {len(union)}", flush=True)

    summaries = esummary("gds", union)

    with open(os.path.join(OUT, "search_raw_summaries.json"), "w") as f:
        json.dump({"queries": QUERIES, "uids_by_query": all_uids,
                   "summaries": summaries}, f, indent=1)

    # Flatten to a candidate TSV (series only).
    rows = []
    for uid, s in summaries.items():
        if s.get("entrytype") != "GSE":
            continue
        rows.append({
            "uid": uid,
            "accession": s.get("accession", ""),
            "gpl": s.get("gpl", ""),
            "n_samples": s.get("n_samples", ""),
            "taxon": s.get("taxon", ""),
            "pdat": s.get("pdat", ""),
            "suppfile": s.get("suppfile", ""),
            "ftplink": s.get("ftplink", ""),
            "title": (s.get("title", "") or "").replace("\t", " ").replace("\n", " "),
            "summary": (s.get("summary", "") or "").replace("\t", " ").replace("\n", " "),
        })
    rows.sort(key=lambda r: r["accession"])
    cols = ["accession", "uid", "taxon", "n_samples", "pdat", "gpl",
            "suppfile", "ftplink", "title", "summary"]
    with open(os.path.join(OUT, "search_candidates.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"GSE candidates: {len(rows)} -> results/fable_geo_2015_2018/search_candidates.tsv")


if __name__ == "__main__":
    main()
