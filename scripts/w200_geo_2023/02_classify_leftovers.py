#!/usr/bin/env python3
"""
Stage 02 - classify every 2023 candidate and mark leftovers.

"Leftover" is defined against the earlier fable_geo_2022_2023 slice:
  * GSE207422 is the only 2023 *lung* series that was actually tested for
    TACSTD2/CLDN4 vs ICI outcome. It is *already analysed*, not leftover.
  * GSE243238 is a 2023 melanoma cross-check (not lung). Not leftover.
  * Every other 2023 human lung ICI series is leftover, including series
    the earlier sweep listed but never probed for gene presence / outcome.

This stage is text-only (title + summary + gdsType). Stage 03 fetches
GEO metadata and supplementary files for leftover series.

Output: results/w200/GEO_2023/leftover_classification.csv
"""
import csv
import re
from pathlib import Path

from geo_common import OUT

ALREADY_ANALYSED = {
    "GSE207422": "already analysed in fable_geo_2022_2023 (NSCLC neoadjuvant anti-PD-1+chemo, TACSTD2/CLDN4 vs MPR)",
    "GSE243238": "already analysed as non-lung melanoma cross-check in fable_geo_2022_2023",
}

LUNG_RE = re.compile(
    r"\b(lung|nsclc|non[- ]?small[- ]?cell|luad|lusc|sclc|"
    r"small[- ]cell lung|pulmonary|lung adenocarcinoma|lung squamous|"
    r"adenocarcinoma of the lung)\b",
    re.I)
# False-positive ICI hits: PD-1/PD-L1 as a gene/stain without treatment.
ICI_RE = re.compile(
    r"\b(immune checkpoint|checkpoint inhibitor|checkpoint blockade|"
    r"immunotherap|anti[- ]?pd[- ]?l?1|anti[- ]?ctla|"
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|cemiplimab|"
    r"avelumab|ipilimumab|tislelizumab|sintilimab|camrelizumab|"
    r"toripalimab|serplulimab|penpulimab|sugemalimab|adebrelimab|"
    r"\bici\b|pd[- ]?1 blockade|pd[- ]?l1 blockade)\b",
    re.I)
# Weaker ICI cue: PD-1/PD-L1 mentioned as a treatment context, not just a stain.
PD_TREAT_RE = re.compile(
    r"\b(anti[- ]?pd|pd[- ]?[l1]{0,2} (blockade|inhibitor|therapy|mab)|"
    r"checkpoint)\b",
    re.I)
EXPR_RE = re.compile(
    r"(expression profiling|rna-?seq|microarray|transcriptom|scrna|"
    r"single[- ]cell rna|counts|tpm|fpkm)",
    re.I)
SC_RE = re.compile(r"\b(single[- ]cell|scRNA|snRNA|CITE-seq|scTCR|10[xX])\b", re.I)
CELL_RE = re.compile(r"\b(cell[- ]line|A549|H1975|H1299|PC9|in vitro)\b", re.I)
MOUSE_RE = re.compile(r"\b(mouse|mice|murine|xenograft|organoid)\b", re.I)
PLATELET_RE = re.compile(r"\b(platelet|TEP|thrombocyte)\b", re.I)
PBMC_RE = re.compile(r"\b(PBMC|peripheral blood|whole blood|circulating)\b", re.I)
ATAC_RE = re.compile(r"\b(ATAC-?seq|ChIP-?seq|cut&tag|Hi-C)\b", re.I)
FLOW_RE = re.compile(r"\b(flow cytometr|CyTOF|spectral flow)\b", re.I)
PANEL_RE = re.compile(r"\b(NanoString|nCounter|targeted panel|IO360|PanCancer)\b", re.I)
OUTCOME_RE = re.compile(
    r"\b(response|responder|non[- ]?responder|resistan|sensitiv|"
    r"progression[- ]free|overall survival|\bpfs\b|\bos\b|\borr\b|"
    r"clinical benefit|outcome|efficacy|recist|MPR|pCR|pathologic response|"
    r"major pathologic)\b",
    re.I)


def classify_text(title, summary, gdstype, taxon):
    text = f"{title} {summary} {gdstype}"
    is_human = "Homo sapiens" in str(taxon)
    is_lung = bool(LUNG_RE.search(text))
    is_ici = bool(ICI_RE.search(text) or PD_TREAT_RE.search(text))
    is_expr = bool(EXPR_RE.search(text) or "Expression profiling" in str(gdstype))
    has_outcome_text = bool(OUTCOME_RE.search(text))

    cats = []
    if SC_RE.search(text) or "single cell" in str(gdstype).lower():
        cats.append("single_cell")
    if CELL_RE.search(text):
        cats.append("cell_line")
    if MOUSE_RE.search(text) and not is_human:
        cats.append("mouse")
    elif MOUSE_RE.search(text):
        cats.append("has_mouse_or_xeno")
    if PLATELET_RE.search(text):
        cats.append("platelet")
    if PBMC_RE.search(text):
        cats.append("blood_or_pbmc")
    if ATAC_RE.search(text):
        cats.append("chromatin")
    if FLOW_RE.search(text):
        cats.append("flow_or_cytof")
    if PANEL_RE.search(text):
        cats.append("targeted_panel")
    if not cats and is_expr:
        cats.append("bulk_or_other")
    if not cats:
        cats.append("other")

    relevant = is_human and is_lung and is_ici and is_expr
    return {
        "is_human": is_human,
        "is_lung": is_lung,
        "is_ici": is_ici,
        "is_expression": is_expr,
        "mentions_outcome": has_outcome_text,
        "relevant_text": relevant,
        "study_category": ";".join(cats),
    }


def main():
    cand = OUT / "search_candidates_2023.csv"
    rows = list(csv.DictReader(cand.open()))
    out_rows = []
    for r in rows:
        flags = classify_text(r["title"], r["summary"], r["gdsType"], r["taxon"])
        acc = r["accession"]
        leftover = acc not in ALREADY_ANALYSED
        already = ALREADY_ANALYSED.get(acc, "")
        if leftover and flags["relevant_text"]:
            leftover_status = "leftover_relevant"
        elif leftover and flags["is_lung"] and flags["is_ici"]:
            leftover_status = "leftover_lung_ici_not_expr_or_not_human"
        elif leftover:
            leftover_status = "leftover_not_relevant"
        else:
            leftover_status = "already_analysed"
        out_rows.append({
            **{k: r[k] for k in ("accession", "title", "summary", "pdat",
                                 "taxon", "n_samples", "gdsType", "gpl")},
            **flags,
            "leftover": leftover,
            "leftover_status": leftover_status,
            "already_analysed_note": already,
        })

    path = OUT / "leftover_classification.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    n = len(out_rows)
    n_left = sum(1 for x in out_rows if x["leftover"])
    n_rel = sum(1 for x in out_rows if x["leftover_status"] == "leftover_relevant")
    print(f"candidates={n} leftover={n_left} leftover_relevant={n_rel}")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
