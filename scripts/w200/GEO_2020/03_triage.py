#!/usr/bin/env python3
"""Triage 2020 GEO hits into leftover vs already-done, with honest reasons.

A series is a leftover only if:
  * PDAT is in calendar 2020
  * it was not already analyzed for TACSTD2/CLDN4 vs ICI outcome
    (prior waves may have listed it without measuring the genes)

Rule-based flags are a first pass. Case-by-case leftover review happens in
04_probe_leftovers.py against series-matrix headers and supplementary files.
"""
import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "w200" / "GEO_2020"

# Already measured TACSTD2/CLDN4 vs ICI outcome in a prior wave.
ALREADY_ANALYZED = {
    "GSE126044": "2019-2021 + ici-bulk waves: tumor RNA-seq, R vs NR, both genes present",
    "GSE135222": "2019-2021 + ici-bulk waves: tumor RNA-seq, PFS (PDAT 2019, not 2020 leftover)",
    "GSE136961": "2019-2021 + ici-bulk: NSCLC anti-PD-1 DCB/NDB, Oncomine 395-gene panel lacks TACSTD2/CLDN4",
    "GSE111414": "2019-2021: PBMC CD8+ T cells, epithelial-null control",
    "GSE182328": "2019-2021: Akkermansia surrogate (PDAT later than 2020)",
    "GSE166449": "ici-bulk: 22 lung, R vs NR, both genes present (PDAT 2021)",
    "GSE207422": "2022-2023 + ici-bulk: neoadjuvant anti-PD-1 + chemo, MPR vs NMPR",
}

LUNG = re.compile(
    r"lung|nsclc|non-small|luad|lusc|sclc|pulmonary|adenocarcinoma|squamous",
    re.I,
)
ICI = re.compile(
    r"immune checkpoint|checkpoint|immunotherap|pd-?1|pd-?l1|ctla|"
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|ipilimumab|"
    r"avelumab|cemiplimab|anti-pd|toripalimab|sintilimab|camrelizumab",
    re.I,
)
SC = re.compile(r"single[- ]cell|scrna|single cell|cite-?seq|10x|citeseq", re.I)
METH = re.compile(r"methylat|epigenetic|methylom|chip-?seq", re.I)
CELLLINE = re.compile(
    r"cell line|A549|H1975|HCC515|knockdown|CRISPR|dox-inducible|in vitro",
    re.I,
)
NOT_TUMOR_LUNG = re.compile(
    r"oral tongue|glioblastoma|melanoma|breast|prostate|pancreatic|"
    r"urothelial|sarcoidosis|HIV infection|SARS-CoV-2|COVID",
    re.I,
)


def load_prior():
    path = RES / "prior_gse_coverage.tsv"
    prior = {}
    if not path.exists():
        return prior
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            prior[row["accession"]] = row
    return prior


def classify(r, prior_row):
    acc = r["accession"]
    text = f"{r.get('title', '')} {r.get('summary', '')}"
    taxon = r.get("taxon") or ""
    pdat = r.get("pdat") or ""
    in_2020 = pdat.startswith("2020")

    if acc in ALREADY_ANALYZED:
        return "ALREADY_ANALYZED", ALREADY_ANALYZED[acc]
    if "Homo sapiens" not in taxon:
        return "EXCLUDED_non_human", f"taxon={taxon}"
    lung = bool(LUNG.search(text))
    ici = bool(ICI.search(text))
    if not lung:
        return "EXCLUDED_not_lung", "no lung/NSCLC term in title+summary"
    if NOT_TUMOR_LUNG.search(r.get("title") or "") and not re.search(
        r"lung cancer|NSCLC|LUAD|LUSC|SCLC", text, re.I
    ):
        return "EXCLUDED_not_lung_cancer", "lung word present but disease is not lung cancer"
    if not ici:
        return "EXCLUDED_not_ICI", "lung but no ICI/immunotherapy context"
    if SC.search(text):
        return "LUNG_ICI_single_cell", "lung+ICI but single-cell (not patient-level bulk outcome)"
    if METH.search(text):
        return "LUNG_ICI_methylation", "lung+ICI but methylation/epigenetic/ChIP assay"
    if CELLLINE.search(text):
        return "LUNG_ICI_cellline_invitro", "lung+ICI but cell-line/in-vitro"
    leftover_flag = "leftover_2020" if in_2020 else "pdat_not_2020"
    return (
        "LUNG_ICI_LEFTOVER_CANDIDATE",
        f"{leftover_flag}; needs series-matrix + suppl probe "
        f"(prior_waves={prior_row.get('prior_geo_waves') or 'none'})",
    )


def main():
    recs = json.loads((RES / "candidates_metadata.json").read_text())
    prior = load_prior()
    rows = []
    for r in recs:
        acc = r["accession"]
        prior_row = prior.get(acc, {})
        cat, reason = classify(r, prior_row)
        rows.append({
            "accession": acc,
            "pdat": r.get("pdat"),
            "category": cat,
            "reason": reason,
            "taxon": r.get("taxon"),
            "n_samples": r.get("n_samples"),
            "gdsType": r.get("gdsType"),
            "gpl": r.get("gpl"),
            "pubmed": ";".join(r.get("pubmedids") or []),
            "suppFile": r.get("suppFile"),
            "prior_geo_waves": prior_row.get("prior_geo_waves", ""),
            "in_dedicated_path": prior_row.get("in_dedicated_path", "N"),
            "title": r.get("title"),
        })

    order = {
        "ALREADY_ANALYZED": 0,
        "LUNG_ICI_LEFTOVER_CANDIDATE": 1,
        "LUNG_ICI_single_cell": 2,
        "LUNG_ICI_methylation": 3,
        "LUNG_ICI_cellline_invitro": 4,
        "EXCLUDED_not_ICI": 5,
        "EXCLUDED_not_lung_cancer": 6,
        "EXCLUDED_not_lung": 7,
        "EXCLUDED_non_human": 8,
    }
    rows.sort(key=lambda x: (order.get(x["category"], 99), x["accession"]))
    out = RES / "tables" / "triage_2020.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys()) if rows else []
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    counts = Counter(r["category"] for r in rows)
    print(f"Total 2020-search candidates: {len(rows)}")
    for cat, n in sorted(counts.items(), key=lambda kv: order.get(kv[0], 99)):
        print(f"  {cat}: {n}")
    leftover = [r for r in rows if r["category"] == "LUNG_ICI_LEFTOVER_CANDIDATE"]
    print(f"Leftover candidates: {len(leftover)}")
    for r in leftover:
        print(f"  {r['accession']} n={r['n_samples']} pdat={r['pdat']}  {r['title'][:80]}")
    print("Wrote", out)


if __name__ == "__main__":
    main()
