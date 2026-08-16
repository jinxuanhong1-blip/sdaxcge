#!/usr/bin/env python3
"""Build the leftover-2019 triage table.

Every row is a real GEO series returned by 01_search_geo.py (PDAT 2019).
GSE135222 is the only 2019 series already analyzed for TACSTD2/CLDN4 vs ICI
outcomes in cursor/fable-geo-2019-2021-c71e; it is listed as SKIP_FABLE_COVERED
and is not re-analyzed here.

Verdicts are from deposited series/sample metadata (esummary + SOFT + series
matrix / GSM brief), not from invented clinical labels.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEARCH = ROOT / "results" / "w200" / "GEO_2019" / "search"
TRIAGE = ROOT / "results" / "w200" / "GEO_2019" / "triage"
TABLES = ROOT / "results" / "w200" / "GEO_2019" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)
TRIAGE.mkdir(parents=True, exist_ok=True)

# Case-by-case review. Keys must match real accessions in candidates_metadata.json.
REVIEW = {
    "GSE135222": {
        "verdict": "SKIP_FABLE_COVERED",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "NSCLC tumor RNA-seq with anti-PD-1/PD-L1 PFS; TACSTD2/CLDN4 already "
            "tested in cursor/fable-geo-2019-2021-c71e (Cox HR ~1.04/1.05, n=27). "
            "Not re-analyzed."
        ),
    },
    "GSE119144": {
        "verdict": "LUNG_METHYLATION_NO_ICI_LABELS",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Illumina methylation companion of the same Jung et al. paper as "
            "GSE135222 (n=60 NSCLC). Sample characteristics are only "
            "'disease state: tumor'. No ICI drug, response, or PFS is deposited. "
            "Cannot test TACSTD2/CLDN4 promoter methylation vs ICI outcome from GEO."
        ),
    },
    "GSE111894": {
        "verdict": "LUNG_SC_T_CELL_NO_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "scRNA-seq of CD3+ tissue-resident memory T cells from lung tumor / "
            "non-involved lung (Clarke et al.). No ICI treatment or outcome. "
            "Epithelial TACSTD2/CLDN4 are not a T-cell phenotype."
        ),
    },
    "GSE111898": {
        "verdict": "LUNG_SC_T_CELL_NO_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "SuperSeries of GSE111894 (sorted CD8 TRM/non-TRM). Same limitation: "
            "no ICI, T-cell compartment only."
        ),
    },
    "GSE117049": {
        "verdict": "LUNG_NOT_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Early-stage LUAD hypoxia lncRNA (NLUCAT1) SuperSeries. Sample "
            "characteristics include histology/stage/recurrence/vital status but "
            "no immunotherapy. Hit the search only via the substring 'ICI' inside "
            "unrelated words (e.g. cisplatin-response). Not an ICI cohort."
        ),
    },
    "GSE117570": {
        "verdict": "LUNG_SC_NO_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "10x scRNA-seq, 4 treatment-naive early NSCLC tumor/normal pairs. "
            "Mentions immunotherapy only as future design insight. Series matrix "
            "has no gene table. No ICI exposure or outcome."
        ),
    },
    "GSE118933": {
        "verdict": "NOT_CANCER_IPF",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Bulk RNA-seq of invasive fibroblasts from idiopathic pulmonary "
            "fibrosis. PD-L1 is discussed as a fibrosis checkpoint, not as "
            "lung-cancer ICI therapy. No tumors, no ICI outcomes."
        ),
    },
    "GSE120028": {
        "verdict": "NOT_LUNG_PRIMARY_NOT_PD1",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "IT1208 anti-CD4 depleting antibody phase 1, 3'-SAGE-seq. Sources: "
            "gastric, gastro-esophageal, CRC (including one bronchus metastasis), "
            "pancreatic liver met, non-cancer esophagus/stomach. No primary lung "
            "cancer. Agent is anti-CD4, not PD-1/PD-L1/CTLA-4. Series matrix has "
            "no gene table. No deposited response labels."
        ),
    },
    "GSE120101": {
        "verdict": "NOT_LUNG_PRIMARY_NOT_PD1",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "TCR-seq companion of GSE120028 (same IT1208 trial). Blood/TIL CD4/CD8 "
            "repertoires, not tumor epithelial expression. Same non-lung, non-PD-1 "
            "limitations."
        ),
    },
    "GSE124574": {
        "verdict": "NOT_LUNG",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": "NanoString of melanoma NK-cell signatures with PFS. Not lung.",
    },
    "GSE127462": {
        "verdict": "NOT_LUNG_ICI_COHORT",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "CIBERSORTx methods paper, Affymetrix bulk of follicular lymphoma and "
            "other dissection samples. Mentions immunotherapy only in the methods "
            "outlook. Not a lung ICI outcome cohort."
        ),
    },
    "GSE127465": {
        "verdict": "LUNG_SC_MYELOID_NO_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "scRNA-seq of tumor-infiltrating myeloid cells in human/mouse lung "
            "cancer. Human samples are surgical NSCLC (treatment-prior flag = 0). "
            "Immunotherapy is mentioned as a future TIM target, not as the "
            "treatment under study. No ICI outcome; not bulk epithelium."
        ),
    },
    "GSE127471": {
        "verdict": "LUNG_SC_N1_NO_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Single scRNA-seq PBMC sample from one NSCLC patient, deposited as a "
            "CIBERSORTx reference. No ICI, n=1."
        ),
    },
    "GSE127825": {
        "verdict": "NOT_LUNG",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "RNA-seq of human thymic / medullary thymic epithelial cells for "
            "tumor-specific antigen discovery. Not lung cancer, not ICI-treated."
        ),
    },
    "GSE129380": {
        "verdict": "CELLLINE_CHIP",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "ChIP-seq (H3K27me3/H3K4me3/H3K27ac) including mouse SCLC line RP-116 "
            "plus non-lung lines. No patient ICI outcomes; not expression of "
            "TACSTD2/CLDN4."
        ),
    },
    "GSE129381": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "RNA-seq of K562 and related lines under polycomb / IFN perturbation. "
            "SCLC mentioned in the SuperSeries; no human lung ICI cohort."
        ),
    },
    "GSE133605": {
        "verdict": "LUNG_N1_NOT_ICI",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "n=1 third-party reanalysis of an ALK+/PD-L1+ LUAD biopsy. PD-L1 is a "
            "tumor marker, not ICI treatment. Series matrix has no gene table. "
            "Cannot test a gene–outcome association."
        ),
    },
    "GSE135164": {
        "verdict": "LUNG_PDC_NO_ICI_OUTCOME",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "15 NSCLC patient-derived cell cultures. Treatment history includes "
            "naive / chemo / TKI; exactly one line (PDC6) is annotated "
            "'Chemotherapy, immunotherapy' with no drug name, no response, and "
            "no survival. Series matrix has no gene table. n=1 immunotherapy "
            "exposure is not an ICI outcome analysis."
        ),
    },
    "GSE135976": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Agilent array, n=4 Hs746T (gastric MET-amplified) MET siRNA/inhibitor. "
            "Not lung; not ICI. Immunotherapy is mentioned only as a hypothetical "
            "combination."
        ),
    },
    "GSE136932": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": "A549/H1975 pirfenidone time-course arrays. Checkpoint therapy is mentioned only as clinical background. No ICI treatment arm.",
    },
    "GSE136933": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": "A549 pirfenidone 6 h arrays. Same in-vitro, non-ICI limitation as GSE136932.",
    },
    "GSE136934": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": "A549 pirfenidone 24 h arrays. Same in-vitro, non-ICI limitation as GSE136932.",
    },
    "GSE138571": {
        "verdict": "MOUSE_ATAC_SCREEN",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "ATAC-seq of mouse KP cells ± Asf1a KO and human H2009. The paper "
            "studies anti-PD-1 sensitization in mice. No human lung ICI expression "
            "cohort and no TACSTD2/CLDN4 outcome labels."
        ),
    },
    "GSE139555": {
        "verdict": "MIXED_SC_T_CELL_NO_ICI_LABELS",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "scRNA/TCR of CD3+/CD45+ immune cells. 13/32 GSMs are lung "
            "(LUAD/LUSC/LCNEC tumor/normal); the rest are endometrial, CRC, RCC. "
            "No ICI drug or response field is deposited. TACSTD2/CLDN4 are "
            "epithelial and are not informative in sorted T/immune cells."
        ),
    },
    "GSE142620": {
        "verdict": "CELLLINE_INVITRO",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": "A549 chronic IL-1β EMT-memory RNA-seq. PD-L1 appears as a gene of interest, not as therapy. No ICI outcomes.",
    },
    "GSE71799": {
        "verdict": "NOT_CANCER",
        "analyzable_tacstd2_cldn4_ici": "no",
        "reason": (
            "Cystic-fibrosis plasma functional genomics (2015 submission, PDAT "
            "2019/01/01). Hit via 'responder cells' wording. Not lung cancer, not ICI."
        ),
    },
}


def main():
    recs = json.loads((SEARCH / "candidates_metadata.json").read_text())
    rows = []
    missing = []
    for r in recs:
        acc = r["accession"]
        rev = REVIEW.get(acc)
        if rev is None:
            missing.append(acc)
            rev = {
                "verdict": "UNREVIEWED",
                "analyzable_tacstd2_cldn4_ici": "unknown",
                "reason": "No case-by-case review recorded.",
            }
        rows.append({
            "accession": acc,
            "pdat": r.get("pdat", ""),
            "submission_date": r.get("soft_submission_date", ""),
            "n_samples": r.get("n_samples", ""),
            "taxon": r.get("taxon", ""),
            "assay": r.get("soft_type") or r.get("gdsType", ""),
            "pubmed": r.get("pubmed", ""),
            "title": r.get("title", ""),
            "verdict": rev["verdict"],
            "analyzable_tacstd2_cldn4_ici": rev["analyzable_tacstd2_cldn4_ici"],
            "reason": rev["reason"],
        })
    if missing:
        raise SystemExit(f"REVIEW missing accessions: {missing}")

    out = TABLES / "triage_leftover_2019.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    n_yes = sum(1 for row in rows if row["analyzable_tacstd2_cldn4_ici"] == "yes")
    summary = {
        "n_series_2019_search": len(rows),
        "n_analyzable_leftover": n_yes,
        "verdict_counts": counts,
        "fable_skipped": ["GSE135222"],
        "fable_2019_2021_analyzed_but_not_2019_pdat": [
            "GSE126044 (PDAT 2020/02/03)",
            "GSE136961 (PDAT 2020/02/03; Oncomine panel lacks TACSTD2/CLDN4)",
            "GSE111414 (PDAT 2021/03/01; CD8 PBMC control)",
            "GSE182328 (PDAT 2021/08/20)",
        ],
        "honest_result": (
            "After skipping GSE135222, zero leftover 2019 GEO series have "
            "human lung-cancer bulk expression of TACSTD2/CLDN4 plus a deposited "
            "per-sample ICI outcome."
        ),
    }
    (TRIAGE / "triage_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("Wrote", out)


if __name__ == "__main__":
    main()
