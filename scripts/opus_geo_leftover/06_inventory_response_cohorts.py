#!/usr/bin/env python3
"""Full inventory of probed series that are human lung, ICI-related, and carry an outcome field.

This is deliberately looser than 04_select_datasets.py: it lists every lung/ICI series whose GEO
sample characteristics contain any response or survival field, with the compartment, platform and
open-file situation spelled out. It is the audit table behind the final dataset choice.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
MAX_BYTES = 2 * 1024**3

LUNG = re.compile(r"\blung\b|NSCLC|non[- ]small[- ]cell|\bSCLC\b|LUAD|LUSC|pulmonary", re.I)
ICI = re.compile(
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|cemiplimab|"
    r"camrelizumab|sintilimab|tislelizumab|toripalimab|anti-?PD-?1|anti-?PD-?L1|anti-?PD-?\(L\)1|"
    r"immunotherap|immune checkpoint|checkpoint (inhibit|block)|\bICI\b|\bICB\b|PD-?1/PD-?L1",
    re.I,
)
OUTCOME_FIELDS = {
    "recist": re.compile(r"recist|best response|\bBOR\b", re.I),
    "responder_label": re.compile(r"response:|responder|non-?responder", re.I),
    "survival": re.compile(r"os \(|pfs \(|overall survival|progression.free|survival", re.I),
    "pathologic": re.compile(r"\bMPR\b|pathologic|\bpCR\b", re.I),
    "recurrence": re.compile(r"recurrence|relapse", re.I),
}
COMPARTMENT = [
    ("platelet", re.compile(r"thrombocyte|platelet|\bTEP\b", re.I)),
    ("blood", re.compile(r"\bPBMC\b|peripheral blood|whole blood|buffy|\bserum\b|plasma", re.I)),
    ("stool", re.compile(r"\bstool\b|fecal|faecal|gut microbio", re.I)),
    ("airway", re.compile(r"\bBALF\b|bronchoalveolar|sputum", re.I)),
    ("tumour", re.compile(r"tumor|tumour|biops|resect|FFPE|lung tissue|NSCLC", re.I)),
]
SINGLE_CELL = re.compile(
    r"single[- ]cell|scRNA|snRNA|\b10x\b|CITE-seq|spatial|Visium|GeoMx|digital spatial|CosMx|\bDSP\b",
    re.I,
)


def compartment_of(text: str) -> str:
    for name, pat in COMPARTMENT:
        if pat.search(text):
            return name
    return "unknown"


def main() -> int:
    probe = json.loads((OUT / "geo_probe.json").read_text())
    cands = {c["accession"]: c for c in json.loads((OUT / "geo_triage_scored.json").read_text())}

    rows = []
    for r in probe:
        acc = r["accession"]
        cand = cands.get(acc, {})
        chars = " ".join(r["characteristics"])
        blob = " ".join([r["title"], cand.get("summary", ""), chars, cand.get("sample_titles", "")])
        if not (LUNG.search(blob) and ICI.search(blob)):
            continue
        found = [k for k, pat in OUTCOME_FIELDS.items() if pat.search(chars)]
        if not found:
            continue
        n = max((len(c.split("\t")) - 1 for c in r["characteristics"]), default=0)
        usable = [f for f in r["processed_candidates"] if (f["bytes"] or 0) <= MAX_BYTES]
        rows.append(
            {
                "accession": acc,
                "title": r["title"],
                "n_samples": n,
                "compartment": compartment_of(chars) if chars else compartment_of(blob),
                "assay": "single-cell/spatial" if SINGLE_CELL.search(blob) else "bulk",
                "platform": ";".join(r["platform_ids"]),
                "outcome_fields": ";".join(found),
                "series_matrix_rows": r["matrix_data_rows"],
                "open_processed_files": ";".join(f["name"] for f in usable),
                "oversized_files": ";".join(f["name"] for f in r["oversized_files"]),
                "pubmed": ";".join(r["pubmed"]),
            }
        )

    rows.sort(key=lambda x: (x["assay"] != "bulk", x["compartment"] != "tumour", -x["n_samples"]))
    path = OUT / "lung_ici_response_inventory.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for x in rows:
        print(
            f"{x['accession']:<11} n={x['n_samples']:>3} {x['compartment']:<9} {x['assay']:<19} "
            f"{x['outcome_fields']:<34} {x['title'][:58]}"
        )
    print(f"\n[inventory] {len(rows)} lung/ICI series with an outcome field -> {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
