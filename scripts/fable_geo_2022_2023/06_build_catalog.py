#!/usr/bin/env python3
"""
Assemble the final verified catalog of GEO 2022-2023 human lung ICI series.

Merges the search table, the supplementary-file verification, and the
series-matrix characteristics, then adds a study-category heuristic and a
verdict on usability for a TACSTD2/CLDN4 vs ICI-outcome analysis.

Outputs:
  results/fable_geo_2022_2023/geo_lung_ici_catalog.csv   (master table)
  notes/fable_geo_2022_2023/catalog_summary.md           (human summary)
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "fable_geo_2022_2023"
NOTES = ROOT / "notes" / "fable_geo_2022_2023"
NOTES.mkdir(parents=True, exist_ok=True)

# Manual verdicts for the series we inspected in depth.
MANUAL = {
    "GSE207422": ("bulk_tumor_rnaseq", True, True,
                  "PRIMARY: NSCLC neoadjuvant anti-PD-1+chemo, 24 baseline bulk "
                  "tumours, whole-transcriptome log2TPM, MPR/NMPR + RECIST + "
                  "residual tumour %. Both genes present."),
    "GSE162520": ("targeted_panel", True, False,
                  "NSCLC PD-1/PD-L1 cohort with OS/PFS (n=92) but expression is "
                  "a ~2.5k-gene targeted panel WITHOUT TACSTD2/CLDN4."),
    "GSE161537": ("targeted_panel", True, False,
                  "NSCLC immunotherapy cohort with RECIST/OS/PFS (n=82) but "
                  "~2.5k-gene targeted panel WITHOUT TACSTD2/CLDN4."),
    "GSE216297": ("platelet_rnaseq", True, False,
                  "NSCLC baseline platelet (TEP) RNA for nivolumab response "
                  "(n=286); platelet not tumour, matrix in .RData, no per-sample "
                  "response label in GEO metadata."),
    "GSE243238": ("bulk_tumor_rnaseq", False, True,
                  "CROSS-CHECK (NOT lung): acral melanoma, ICI-treated, bulk raw "
                  "counts, clinical-benefit label. Both genes present."),
    "GSE235500": ("single_cell", True, False,
                  "Tumour-infiltrating Treg scRNA-seq (checkpoint blockade); "
                  "sorted immune cells, not bulk tumour epithelium."),
    "GSE235603": ("single_cell", True, False,
                  "SuperSeries of Treg scRNA-seq; sorted immune cells."),
    "GSE185206": ("single_cell", True, False,
                  "scRNA/TCR-seq of T cells during checkpoint blockade."),
    "GSE160903": ("single_cell", True, False,
                  "scRNA of T cells (mouse+human); mechanistic."),
    "GSE214992": ("cell_line", True, False,
                  "NSCLC cell lines +/- CD8 T-cell killing; not patient outcome."),
    "GSE190731": ("xenograft", True, False,
                  "EGFR-mutant NSCLC xenografts, durvalumab/oleclumab; not "
                  "patient outcome."),
    "GSE248830": ("targeted_panel", True, False,
                  "Brain-metastasis (breast+lung) NanoString-scale panel; not an "
                  "ICI response cohort."),
}

CAT_RE = [
    ("single_cell", re.compile(r"single[- ]cell|scRNA|10x|10X|CITE-seq|TCR-seq", re.I)),
    ("cell_line", re.compile(r"cell line|cell-line", re.I)),
    ("xenograft", re.compile(r"xenograft|PDX|mouse model", re.I)),
]


def categorize(title, summary, gdstype):
    text = f"{title} {summary} {gdstype}"
    for name, rx in CAT_RE:
        if rx.search(text):
            return name
    return "bulk_or_other"


def main():
    ver = pd.read_csv(OUT / "geo_verified.csv")
    chars = pd.read_csv(OUT / "series_characteristics_summary.csv")
    m = ver.merge(chars[["accession", "char_keys", "outcome_keys",
                         "has_outcome_annotation", "n_table_rows"]],
                  on="accession", how="left")

    cats, is_lung, usable, verdicts = [], [], [], []
    for _, r in m.iterrows():
        acc = r["accession"]
        if acc in MANUAL:
            cat, lung, use, note = MANUAL[acc]
        else:
            cat = categorize(r["title"], r.get("summary", ""), r["gdstype"])
            lung = bool(r["is_lung"])
            use = ""  # not deeply verified
            note = ""
        cats.append(cat)
        is_lung.append(lung)
        usable.append(use)
        verdicts.append(note)
    m["study_category"] = cats
    m["is_lung_confirmed"] = is_lung
    m["usable_for_TACSTD2_CLDN4"] = usable
    m["verdict"] = verdicts

    cols = ["accession", "title", "pdat", "n_samples", "taxon",
            "is_human", "is_lung_confirmed", "is_ici", "is_expression",
            "relevant", "study_category", "outcome_keys",
            "has_processed_suppl", "suppl_total_mb", "under_2gb",
            "usable_for_TACSTD2_CLDN4", "verdict", "suppl_files"]
    cat = m[cols].sort_values(["relevant", "n_samples"], ascending=[False, False])
    cat.to_csv(OUT / "geo_lung_ici_catalog.csv", index=False)

    rel = cat[cat.relevant]
    lines = []
    lines.append("# GEO 2022-2023 human lung ICI series - verified catalog\n")
    lines.append(f"- Candidate GSE series returned by the exhaustive search: "
                 f"**{len(cat)}**\n")
    lines.append(f"- Passing text relevance (human + lung + ICI + expression): "
                 f"**{len(rel)}**\n")
    lines.append(f"- With per-sample outcome annotation in GEO metadata: "
                 f"**{int((rel.outcome_keys.fillna('')!='').sum())}**\n")
    lines.append(f"- Open, processed supplementary download <2 GB: "
                 f"**{int((rel.has_processed_suppl & rel.under_2gb).sum())}**\n")
    lines.append("\n## Deeply verified series\n")
    lines.append("| GSE | n | category | lung | usable for TACSTD2/CLDN4 | note |\n")
    lines.append("|---|---|---|---|---|---|\n")
    for acc in MANUAL:
        row = cat[cat.accession == acc]
        if row.empty:
            continue
        r = row.iloc[0]
        lines.append(f"| {acc} | {r['n_samples']} | {r['study_category']} | "
                     f"{r['is_lung_confirmed']} | {r['usable_for_TACSTD2_CLDN4']} | "
                     f"{r['verdict']} |\n")
    lines.append("\n## All relevant series (sorted by n)\n")
    lines.append("| GSE | n | date | category | outcome keys |\n")
    lines.append("|---|---|---|---|---|\n")
    for _, r in rel.iterrows():
        lines.append(f"| {r['accession']} | {r['n_samples']} | {r['pdat']} | "
                     f"{r['study_category']} | {str(r['outcome_keys'])[:80]} |\n")
    (NOTES / "catalog_summary.md").write_text("".join(lines))
    print("Wrote catalog + summary.")
    print(cat[["accession", "n_samples", "study_category",
               "usable_for_TACSTD2_CLDN4"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
