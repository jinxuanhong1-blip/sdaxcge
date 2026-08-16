#!/usr/bin/env python3
"""Triage the GEO search hits down to plausible bulk human lung ICI cohorts.

Scores each candidate on lung specificity, ICI treatment, response/outcome annotation, and
bulk-vs-single-cell platform, then emits a ranked shortlist for manual accession verification.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "results" / "opus_geo_leftover"

SINGLE_CELL = re.compile(
    r"single[- ]cell|single cell|scRNA|snRNA|10x genomics|10X Genomics|smart-seq|CITE-seq|"
    r"single-nucleus|spatial transcriptom|Visium|GeoMx|CosMx|MERFISH|scATAC",
    re.I,
)
LUNG = re.compile(r"\blung\b|NSCLC|non-small[- ]cell|SCLC|small[- ]cell lung|pulmonary|LUAD|LUSC", re.I)
ICI = re.compile(
    r"immunotherap|immune[- ]checkpoint|checkpoint inhibit|checkpoint block|anti-?PD-?L?1|"
    r"PD-?1\b|PD-?L1\b|CTLA-?4|nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|"
    r"ipilimumab|cemiplimab|camrelizumab|sintilimab|tislelizumab|toripalimab|\bICI\b|\bICB\b|"
    r"serplulimab|penpulimab|adebrelimab|sugemalimab",
    re.I,
)
RESPONSE = re.compile(
    r"respond|response|responder|non-?responder|resistan|refractor|efficacy|"
    r"progression[- ]free|\bPFS\b|\bOS\b|overall survival|durable|benefit|RECIST|"
    r"\bpCR\b|major pathologic|pathological response|\bMPR\b|outcome|sensitiv",
    re.I,
)
NEGATIVE_ORG = re.compile(
    r"melanoma|bladder|urothelial|renal cell|kidney|gastric|esophage|colorect|"
    r"hepatocellular|breast|ovarian|head and neck|glioma|lymphoma|myeloma|leukemi|"
    r"nasopharyn|cervical|prostate|pancrea|mesotheliom|thymom|sarcoma",
    re.I,
)
EXPRESSION = re.compile(r"Expression profiling", re.I)
MOUSE_ONLY = re.compile(r"^Mus musculus$|^Rattus", re.I)


def main() -> int:
    rows = json.loads((OUT / "geo_search_candidates.json").read_text())
    scored = []
    for r in rows:
        text = f"{r['title']} {r['summary']}"
        gds = r["gdsType"] or ""
        if MOUSE_ONLY.search(r.get("taxon", "") or ""):
            continue
        if not EXPRESSION.search(gds):
            continue
        sc = bool(SINGLE_CELL.search(text)) or bool(SINGLE_CELL.search(gds))
        score = 0
        score += 3 if LUNG.search(text) else 0
        score += 3 if ICI.search(text) else 0
        score += 2 if RESPONSE.search(text) else 0
        score -= 4 if sc else 0
        # Pan-cancer cohorts are still usable if lung is present, but deprioritise pure other-tumour studies.
        if NEGATIVE_ORG.search(text) and not LUNG.search(text):
            score -= 5
        n = r.get("n_samples") or 0
        try:
            n = int(n)
        except (TypeError, ValueError):
            n = 0
        score += 1 if n >= 10 else 0
        r = {**r, "score": score, "single_cell_flag": sc, "n": n}
        scored.append(r)

    scored.sort(key=lambda x: (-x["score"], x["accession"]))
    shortlist = [r for r in scored if r["score"] >= 7]
    (OUT / "geo_triage_scored.json").write_text(json.dumps(scored, indent=2))
    (OUT / "geo_triage_shortlist.json").write_text(json.dumps(shortlist, indent=2))

    for r in shortlist:
        print(
            f"{r['accession']}\tscore={r['score']}\tn={r['n']}\tsc={int(r['single_cell_flag'])}\t"
            f"{r['gdsType'][:60]}\t{r['title'][:110]}"
        )
    print(f"\n[triage] {len(scored)} expression series, {len(shortlist)} shortlisted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
