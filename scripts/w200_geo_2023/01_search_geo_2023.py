#!/usr/bin/env python3
"""
Stage 01 - exhaustive GEO search for human lung-cancer ICI series with a 2023
publication date.

This re-runs the search independently of the earlier `fable_geo_2022_2023`
slice rather than trusting its candidate list, so that series the earlier
sweep missed can still surface. The union of a broad lung x ICI query plus
per-drug and per-histology queries is de-duplicated by GSE accession.

Output: results/w200/GEO_2023/search_candidates_2023.csv
        notes/w200_geo_2023/search_queries.json
"""
import csv
import json
import re

from geo_common import OUT, NOTES, esearch, esummary

DATE = '("2023/01/01"[PDAT] : "2023/12/31"[PDAT])'
BASE = f'gse[ETYP] AND "Homo sapiens"[Organism] AND {DATE}'

ICI_TERMS = [
    "immune checkpoint", "checkpoint inhibitor", "checkpoint blockade",
    "immunotherapy", "anti-PD-1", "anti-PD1", "PD-1", "PD1",
    "anti-PD-L1", "anti-PDL1", "PD-L1", "PDL1", "CTLA-4", "CTLA4",
    "nivolumab", "pembrolizumab", "atezolizumab", "durvalumab",
    "cemiplimab", "avelumab", "ipilimumab", "tislelizumab",
    "sintilimab", "camrelizumab", "toripalimab", "serplulimab",
    "penpulimab", "sugemalimab", "adebrelimab", "neoadjuvant immunotherapy",
]

LUNG_TERMS = [
    "lung", "NSCLC", "non-small cell lung", "non small cell lung",
    "LUAD", "lung adenocarcinoma", "LUSC", "lung squamous",
    "SCLC", "small cell lung", "pulmonary carcinoma",
]


def main():
    queries = []
    uids = set()

    lung_or = " OR ".join(f'"{t}"' for t in LUNG_TERMS)
    ici_or = " OR ".join(f'"{t}"' for t in ICI_TERMS)

    broad = f'{BASE} AND ({lung_or}) AND ({ici_or})'
    queries.append(broad)
    uids.update(esearch(broad))

    # Per-drug sweeps: a drug name alone can retrieve series whose free text
    # never says "lung" in the fields the broad query matched.
    for drug in ICI_TERMS:
        q = f'{BASE} AND ({lung_or}) AND "{drug}"'
        queries.append(q)
        uids.update(esearch(q))

    # Per-histology sweeps against the generic ICI vocabulary.
    for lt in LUNG_TERMS:
        q = f'{BASE} AND "{lt}" AND ({ici_or})'
        queries.append(q)
        uids.update(esearch(q))

    print(f"unique gds UIDs: {len(uids)}")
    summ = esummary(sorted(uids))

    rows = []
    for uid, rec in summ.items():
        acc = rec.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        rows.append({
            "accession": acc,
            "uid": uid,
            "title": (rec.get("title") or "").replace("\n", " "),
            "summary": (rec.get("summary") or "").replace("\n", " "),
            "pdat": rec.get("pdat") or rec.get("PDAT") or "",
            "taxon": rec.get("taxon", ""),
            "n_samples": rec.get("n_samples", ""),
            "gdsType": rec.get("gdstype") or rec.get("gdsType") or "",
            "gpl": rec.get("gpl") or rec.get("GPL") or "",
        })

    # Keep only true 2023 publication dates (esearch PDAT ranges can bleed).
    rows = [r for r in rows if str(r["pdat"]).startswith("2023")]
    rows.sort(key=lambda r: int(re.sub(r"\D", "", r["accession"]) or 0))

    fields = ["accession", "uid", "title", "summary", "pdat", "taxon",
              "n_samples", "gdsType", "gpl"]
    path = OUT / "search_candidates_2023.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (NOTES / "search_queries.json").write_text(
        json.dumps({"queries": queries, "n_uids": len(uids),
                    "n_gse_2023": len(rows)}, indent=2))
    print(f"wrote {path} ({len(rows)} GSE series dated 2023)")


if __name__ == "__main__":
    main()
