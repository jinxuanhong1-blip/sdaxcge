#!/usr/bin/env python3
"""Rule-based triage of all 90 candidate GEO series into auditable categories.

Every candidate accession is REAL (from the esearch/esummary output). This step
classifies each series by relevance for a tumor-expression vs ICI-outcome
analysis, so the "exhaustive + verify" requirement is fully documented.
"""
import json
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "fable_geo_2019_2021"

recs = json.loads((RES / "candidates_metadata.json").read_text())

# Series that passed manual verification and were downloaded.
ANALYZED = {
    "GSE126044": "tumor RNA-seq counts; responder vs non-responder (anti-PD-1); TACSTD2/CLDN4 present",
    "GSE135222": "tumor RNA-seq TPM; PFS (anti-PD-1/PD-L1); TACSTD2/CLDN4 present",
    "GSE182328": "lung-tumor RNA-seq counts; Akkermansia ICI-prognostic surrogate; TACSTD2/CLDN4 present",
    "GSE111414": "PBMC CD8+ T-cell counts; responder vs non-responder; epithelial-null control",
}
VERIFIED_NO_TARGET = {
    "GSE136961": "NSCLC anti-PD-1 DCB/NDB, but Oncomine 395-gene immune panel lacks TACSTD2/CLDN4",
}

LUNG = re.compile(r"lung|nsclc|non-small|luad|lusc|sclc|pulmonary|adenocarcinoma|squamous", re.I)
ICI = re.compile(r"immune checkpoint|checkpoint|immunotherap|pd-?1|pd-?l1|ctla|nivolumab|pembrolizumab|atezolizumab|durvalumab|ipilimumab|avelumab|cemiplimab|anti-pd", re.I)
SC = re.compile(r"single[- ]cell|scrna|single cell|cite-?seq|10x|citeseq", re.I)
METH = re.compile(r"methylat|epigenetic|methylom", re.I)
CELLLINE = re.compile(r"cell line|A549|H1975|HCC515|knockdown|CRISPR|dox-inducible|in vitro", re.I)


def classify(r):
    acc = r["accession"]
    if acc in ANALYZED:
        return "ANALYZED", ANALYZED[acc]
    if acc in VERIFIED_NO_TARGET:
        return "VERIFIED_NO_TARGET_GENES", VERIFIED_NO_TARGET[acc]
    text = f"{r.get('title','')} {r.get('summary','')}"
    taxon = (r.get("taxon") or "")
    if "Homo sapiens" not in taxon:
        return "EXCLUDED_non_human", f"taxon={taxon}"
    lung = bool(LUNG.search(text))
    ici = bool(ICI.search(text))
    if not lung:
        return "EXCLUDED_not_lung", "no lung/NSCLC term in title+summary"
    if not ici:
        return "EXCLUDED_not_ICI", "lung but no ICI/immunotherapy context"
    # lung + ICI but not suitable bulk-tumor-expression + outcome
    if SC.search(text):
        return "LUNG_ICI_single_cell", "lung+ICI but single-cell (not patient-level bulk outcome)"
    if METH.search(text):
        return "LUNG_ICI_methylation", "lung+ICI but methylation/epigenetic assay"
    if CELLLINE.search(text):
        return "LUNG_ICI_cellline_invitro", "lung+ICI but cell-line/in-vitro"
    return "LUNG_ICI_other", "lung+ICI; needs case-by-case review (small n, targeted, or no per-sample outcome)"


rows = []
for r in recs:
    cat, reason = classify(r)
    rows.append({
        "accession": r["accession"], "category": cat, "reason": reason,
        "taxon": r.get("taxon"), "n_samples": r.get("n_samples"),
        "pubmed": ";".join(r.get("pubmedids") or []),
        "title": r.get("title"),
    })

order = {"ANALYZED": 0, "VERIFIED_NO_TARGET_GENES": 1, "LUNG_ICI_other": 2,
         "LUNG_ICI_single_cell": 3, "LUNG_ICI_methylation": 4,
         "LUNG_ICI_cellline_invitro": 5, "EXCLUDED_not_ICI": 6,
         "EXCLUDED_not_lung": 7, "EXCLUDED_non_human": 8}
rows.sort(key=lambda x: (order.get(x["category"], 99), x["accession"]))

out = RES / "tables" / "triage_all_candidates.csv"
out.parent.mkdir(parents=True, exist_ok=True)
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "category", "reason",
                                      "taxon", "n_samples", "pubmed", "title"])
    w.writeheader()
    w.writerows(rows)

from collections import Counter
counts = Counter(r["category"] for r in rows)
print(f"Total candidates: {len(rows)}")
for cat, n in sorted(counts.items(), key=lambda kv: order.get(kv[0], 99)):
    print(f"  {cat}: {n}")
print("Wrote", out)
