#!/usr/bin/env python3
"""Triage GEO search candidates: score by keywords suggesting a patient ICI
cohort with response annotation vs cell-line/mechanistic studies.

Writes results/w200/GEO_2026/geo_triage.tsv sorted by score desc.
"""
import csv
import re
from pathlib import Path

RES_DIR = Path(__file__).resolve().parents[2] / "results" / "w200" / "GEO_2026"

POS = {
    r"\bresponder": 4, r"non-?responder": 5, r"\bresponse\b": 2,
    r"pre-?treatment": 2, r"post-?treatment": 1, r"baseline": 1,
    r"biops": 2, r"\bpatients?\b": 2, r"\bcohort\b": 2, r"clinical": 1,
    r"progression-?free|PFS": 2, r"overall survival|\bOS\b": 1,
    r"durable clinical benefit|\bDCB\b": 4, r"\bMPR\b|major pathologic": 3,
    r"pathologic(al)? (complete )?response|\bpCR\b": 3,
    r"RECIST": 4, r"neoadjuvant": 2, r"peripheral blood|PBMC|plasma": 1,
    r"resistan": 1, r"\bORR\b": 3,
}
NEG = {
    r"cell lines?": -3, r"\bmice\b|\bmouse\b|murine|xenograft": -3,
    r"knock-?down|knock-?out|overexpress|shRNA|siRNA|CRISPR": -3,
    r"co-?culture": -2, r"organoid": -1, r"\bin vitro\b": -2,
}


def main() -> None:
    rows = []
    with (RES_DIR / "geo_search_candidates.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            text = (r["title"] + " " + r["summary"]).lower()
            score = 0
            hits = []
            for pat, w in {**POS, **NEG}.items():
                if re.search(pat, text, re.I):
                    score += w
                    hits.append(pat)
            rows.append({
                "accession": r["accession"], "score": score,
                "n_samples": r["n_samples"], "gdstype": r["gdstype"],
                "pdat": r["pdat"], "suppfile": r["suppfile"],
                "title": r["title"], "hits": ";".join(hits),
            })
    rows.sort(key=lambda x: -x["score"])
    cols = ["accession", "score", "n_samples", "pdat", "gdstype", "suppfile",
            "title", "hits"]
    with (RES_DIR / "geo_triage.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote", RES_DIR / "geo_triage.tsv")
    for r in rows[:40]:
        print(f"{r['score']:>3} {r['accession']:<11} n={r['n_samples']:<4} "
              f"{r['pdat']} {r['title'][:90]}")


if __name__ == "__main__":
    main()
