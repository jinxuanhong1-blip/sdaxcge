"""Exhaustive GEO search for leftover 2010-2014 lung ICI / PD-1 / PD-L1 SERIES.

This window is leftover relative to the already-mined GEO slices
(2015-2018, 2019-2021, 2022-2023, 2026). No prior 2010-2014 GEO ICI PR
exists in this repo, so leftover = the real NCBI union for this date range.

Strategy: complementary queries against GEO DataSets (gds), union the
series, fetch esummary. Filtering is in 02_verify.py so this script keeps
the raw, un-invented union of real accessions.

Outputs (under results/w200/GEO_2010_2014/):
  search_raw_summaries.json
  search_candidates.tsv
  search_manifest.json
"""
import json
import os
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eutils import esearch, esummary  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")
os.makedirs(OUT, exist_ok=True)

# Same date syntax as the 2015-2018 sibling. Year-only [PDAT] is unreliable
# on gds; the quoted YYYY/MM/DD range is what NCBI actually indexes.
DATE = '("2010/01/01"[PDAT] : "2014/12/31"[PDAT])'
HUMAN = '"Homo sapiens"[Organism]'
GSE = "gse[Entry Type]"

LUNG = (
    "(lung[Title] OR NSCLC[Title] OR \"non-small cell lung\"[All Fields] OR "
    "\"lung cancer\"[All Fields] OR \"lung carcinoma\"[All Fields] OR "
    "\"lung adenocarcinoma\"[All Fields] OR \"lung neoplasms\"[All Fields] OR "
    "\"lung\"[All Fields])"
)

ICI = (
    '("PD-1"[All Fields] OR "PD1"[All Fields] OR "PDCD1"[All Fields] OR '
    '"PD-L1"[All Fields] OR "PDL1"[All Fields] OR "CD274"[All Fields] OR '
    '"CTLA-4"[All Fields] OR "CTLA4"[All Fields] OR '
    '"immune checkpoint"[All Fields] OR "checkpoint inhibitor"[All Fields] OR '
    '"checkpoint blockade"[All Fields] OR "immunotherapy"[All Fields] OR '
    '"anti-PD-1"[All Fields] OR "anti-PD-L1"[All Fields] OR '
    '"nivolumab"[All Fields] OR "pembrolizumab"[All Fields] OR '
    '"atezolizumab"[All Fields] OR "durvalumab"[All Fields] OR '
    '"avelumab"[All Fields] OR "ipilimumab"[All Fields] OR '
    '"tremelimumab"[All Fields])'
)

# Literal PD-1/PD-L1 text probe (recorded because it is expected to be empty
# in this window — the phrasing is largely post-2014).
STRICT_PD = "(lung) AND (PD-L1 OR PD-1 OR CD274 OR PDCD1) AND gse[Entry Type] AND 2010:2014[PDAT]"

QUERIES = {
    "lung_ici_human_gse_2010_2014": f"{LUNG} AND {ICI} AND {HUMAN} AND {GSE} AND {DATE}",
    "lung_ici_gse_2010_2014": f"{LUNG} AND {ICI} AND {GSE} AND {DATE}",
    "strict_pdl1_text_2010_2014": STRICT_PD,
}


def main():
    all_uids = {}
    for name, term in QUERIES.items():
        uids = esearch("gds", term)
        print(f"[{name}] {len(uids)} UIDs", flush=True)
        all_uids[name] = uids

    # Union of the two ICI queries; the strict-text probe is recorded separately.
    union = sorted(
        set(all_uids["lung_ici_human_gse_2010_2014"] + all_uids["lung_ici_gse_2010_2014"]),
        key=int,
    )
    print(f"union UIDs (ICI queries): {len(union)}", flush=True)

    summaries = esummary("gds", union) if union else {}

    with open(os.path.join(OUT, "search_raw_summaries.json"), "w") as f:
        json.dump(
            {"queries": QUERIES, "uids_by_query": all_uids, "summaries": summaries},
            f,
            indent=1,
        )

    rows = []
    for uid, s in summaries.items():
        if s.get("entrytype") != "GSE":
            continue
        rows.append(
            {
                "uid": uid,
                "accession": s.get("accession", ""),
                "gpl": s.get("gpl", ""),
                "n_samples": s.get("n_samples", ""),
                "taxon": s.get("taxon", ""),
                "pdat": s.get("pdat", ""),
                "gdstype": s.get("gdstype", ""),
                "suppfile": s.get("suppfile", ""),
                "ftplink": s.get("ftplink", ""),
                "pubmedids": ";".join(str(p) for p in s.get("pubmedids", [])),
                "title": (s.get("title", "") or "").replace("\t", " ").replace("\n", " "),
                "summary": (s.get("summary", "") or "").replace("\t", " ").replace("\n", " "),
            }
        )
    rows.sort(key=lambda r: r["accession"])
    cols = [
        "accession",
        "uid",
        "taxon",
        "n_samples",
        "pdat",
        "gpl",
        "gdstype",
        "pubmedids",
        "suppfile",
        "ftplink",
        "title",
        "summary",
    ]
    with open(os.path.join(OUT, "search_candidates.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "queries": QUERIES,
        "n_uids_by_query": {k: len(v) for k, v in all_uids.items()},
        "n_union_ici": len(union),
        "n_gse_candidates": len(rows),
        "accessions": [r["accession"] for r in rows],
        "note": (
            "Every accession is a live NCBI gds hit. The strict PD-1/PD-L1 text "
            "probe is recorded separately and is not required for the union."
        ),
    }
    with open(os.path.join(OUT, "search_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"GSE candidates: {len(rows)} -> results/w200/GEO_2010_2014/search_candidates.tsv")


if __name__ == "__main__":
    main()
