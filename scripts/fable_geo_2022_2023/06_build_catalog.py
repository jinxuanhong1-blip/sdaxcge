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
    # ---- leftover pass (every previously unverified relevant / recovered series)
    "GSE221733": ("spatial_dsp", True, True,
                  "RECOVERED leftover: NSCLC GeoMx DSP CTA, immunotherapy-treated, "
                  "41 pts with Responder/Non-responder. Missed by first-pass ICI "
                  "regex (immunotherapy != immunotherap\\b). TACSTD2 on panel; "
                  "CLDN4 absent. Analysed patient-level PanCK+ TACSTD2."),
    "GSE221322": ("spatial_protein", True, False,
                  "Sister DSP protein panel of GSE221733; no TACSTD2/CLDN4 proteins."),
    "GSE248378": ("bulk_tumor_rnaseq", True, True,
                  "LEFTOVER PRIMARY: neoadjuvant durvalumab +/- SBRT, 29 post-Rx "
                  "non-MPR tumours, FPKM with TACSTD2+CLDN4. Recurrence joined from "
                  "Nat Commun source data (DurvaNNN). 45-M-PO arm-discordant."),
    "GSE193049": ("balf_rnaseq", True, True,
                  "LEFTOVER: BALF (not tumour) RNA, PD-1 responders vs non-responders "
                  "n=7; both genes present. Analysed with that caveat."),
    "GSE189045": ("exosomal_mirna", True, False,
                  "Serum exosomal miRNA with anti-PD-1/PD-L1 response labels; no mRNA "
                  "so TACSTD2/CLDN4 cannot be measured."),
    "GSE250262": ("nanostring_io360", True, False,
                  "NSCLC tumour NanoString IO360 RCC; panel has EPCAM not TACSTD2/CLDN4; "
                  "GEO metadata has no ICI response/survival label."),
    "GSE185204": ("single_cell", True, False,
                  "scRNA-seq of CD3+ / multimer+ T cells during ICB; not bulk tumour."),
    "GSE186446": ("bulk_or_other", True, False,
                  "Regional featureCounts from 3 ICB patients (T-cell lineage paper); "
                  "both genes present but no per-sample ICI outcome."),
    "GSE235048": ("pbmc_rnaseq", True, False,
                  "PBMC TPM from NSCLC on various ICI regimens; TACSTD2 present, no "
                  "response label (treatment only)."),
    "GSE228419": ("sorted_tcell", True, False,
                  "Sorted CD8 T cells (PBMC/TIL); both genes present; no ICI outcome."),
    "GSE164146": ("single_cell", True, False,
                  "Treg/CD4 scRNA-seq in lung cancer; not bulk tumour epithelium."),
    "GSE212622": ("infection", False, False,
                  "False positive: murine trypanosome lung infection, not ICI."),
    "GSE224099": ("sorted_tcell", False, False,
                  "Melanoma (not lung) pre-ICI T cells; not a lung tumour cohort."),
    "GSE193719": ("cell_line", True, False,
                  "NSCLC cell lines +/- miR-455-5p; no patient ICI outcome."),
    "GSE189804": ("cell_line", True, False,
                  "A549 NRF2 KO RNA-seq; no patient ICI outcome."),
    "GSE217451": ("cell_line", True, False,
                  "H1650 hMENA siRNA; no patient ICI outcome."),
    "GSE224216": ("cell_line", True, False,
                  "H2030 hMENA siRNA RNA-seq; no patient ICI outcome."),
    "GSE150255": ("cell_line", True, False,
                  "A549/HCC827 IFN-γ time course; no patient ICI outcome."),
    "GSE194350": ("cell_line", True, False,
                  "A549 CRISPR/in-vivo MEN1 screen RNA-seq; no patient ICI outcome."),
    "GSE195770": ("ipf", False, False,
                  "False positive: idiopathic pulmonary fibrosis, not lung cancer ICI."),
    "GSE218402": ("cell_line", True, False,
                  "A549 + PBMC norepinephrine/adenosine; no patient ICI outcome."),
    "GSE224246": ("cell_line", True, False,
                  "H1650 ATAC-seq (hMENA); chromatin, not mRNA outcome cohort."),
    "GSE238006": ("cell_line", True, False,
                  "SCLC cell-line decitabine RNA-seq; no patient ICI outcome."),
    "GSE250254": ("sorted_tcell", True, False,
                  "Sorted CD4 memory T cells (COPD/NSCLC); no ICI outcome."),
    "GSE229353": ("single_cell", True, False,
                  "CD45+ scRNA-seq, chemo vs anti-PD-1+chemo; sorted immune cells."),
    "GSE178521": ("cell_line", True, False,
                  "NCI-H460 PRMT5 knockdown; no patient ICI outcome."),
    "GSE192591": ("sorted_tcell", True, False,
                  "Blood CD8 T cells +/- IL-27; no tumour ICI outcome."),
    "GSE192790": ("cell_line", True, False,
                  "A549 tumoursphere RNA-seq; no patient ICI outcome."),
    "GSE197236": ("cell_line", True, False,
                  "A549 radioresistance lines; not ICI."),
    "GSE198099": ("single_cell", True, False,
                  "NSCLC TME scRNA-seq, n=2 patients; no ICI outcome label."),
    "GSE213590": ("cell_line", True, False,
                  "PC-9 c-Jun overexpression; no patient ICI outcome."),
    "GSE213902": ("single_cell", True, False,
                  "PBMC CD3+ scRNA/TCR on chemo-IO; circulating T cells, not tumour "
                  "TACSTD2/CLDN4."),
    "GSE223779": ("organoid", True, False,
                  "ALK+ tumour organoids; treatment field is sample origin, not ICI "
                  "outcome."),
    "GSE193707": ("bulk_tumor_rnaseq", True, False,
                  "Mediastinal LN after neoadjuvant chemo (not ICI); n=4."),
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
    recov = cat[cat.verdict.astype(str).str.contains("RECOVERED|LEFTOVER", na=False)]
    lines.append(f"- Leftover / recovered series given a deep verdict this pass: "
                 f"**{len(recov)}**\n")
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
