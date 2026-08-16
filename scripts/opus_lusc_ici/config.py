"""Shared configuration for the LUSC / ICI TACSTD2-CLDN4 slice.

All paths, cohort definitions and analysis constants live here so that every
step of the pipeline reads from a single source of truth.

Scope of this slice (do not write outside these trees):
    notes/opus_lusc_ici/
    scripts/opus_lusc_ici/
    results/opus_lusc_ici/

Bulk downloads are kept in an out-of-tree cache so the repository only ever
carries small, processed artefacts.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SLICE = "opus_lusc_ici"

SCRIPTS_DIR = REPO_ROOT / "scripts" / SLICE
NOTES_DIR = REPO_ROOT / "notes" / SLICE
RESULTS_DIR = REPO_ROOT / "results" / SLICE

DATA_DIR = RESULTS_DIR / "data"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
LOGS_DIR = RESULTS_DIR / "logs"
MANIFEST_DIR = RESULTS_DIR / "manifest"

# Raw downloads never enter the repository.
CACHE_DIR = Path(os.environ.get("OPUS_LUSC_CACHE", "/tmp/opus_lusc_ici_cache"))
RAW_DIR = CACHE_DIR / "raw"

for _d in (DATA_DIR, TABLES_DIR, FIGURES_DIR, LOGS_DIR, MANIFEST_DIR, RAW_DIR):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20260816

# Genes of interest. TROP2 (TACSTD2) and CLDN4 are the ADC-relevant surface
# targets under study; the remaining entries are used as interpretive anchors.
TARGET_GENES = ["TACSTD2", "CLDN4"]
COMPANION_GENES = ["CLDN3", "CLDN7", "EPCAM", "CD274", "PDCD1", "CTLA4"]

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
XENA_GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"


@dataclass(frozen=True)
class Cohort:
    """One public ICI cohort contributing lung-squamous samples."""

    key: str                  # short internal identifier
    gse: str                  # GEO series accession
    label: str                # human readable label used in figures/tables
    setting: str              # advanced | neoadjuvant
    regimen: str
    histology_source: str     # 'pathology' or 'inferred'
    outcome: str              # what the cohort contributes to the response meta-analysis
    supplementary: tuple = field(default=())

    @property
    def ftp_dir(self) -> str:
        stem = self.gse[:-3] + "nnn"
        return f"{GEO_FTP}/{stem}/{self.gse}"


COHORTS = [
    Cohort(
        key="FRANCE3",
        gse="GSE190265",
        label="Dijon France-3 (GSE190265)",
        setting="advanced",
        regimen="anti-PD-1 monotherapy (1st/2nd line)",
        histology_source="pathology",
        outcome="PFS",
        supplementary=(
            "GSE190265_TPM_France3.csv.gz",
            "GSE190265_samples_info_France3.csv.gz",
        ),
    ),
    Cohort(
        key="FRANCE4",
        gse="GSE190266",
        label="Dijon France-4 (GSE190266)",
        setting="advanced",
        regimen="anti-PD-1 monotherapy (1st/2nd line)",
        histology_source="pathology",
        outcome="PFS (6-month landmark)",
        supplementary=("GSE190266_TPM_France4.csv.gz",),
    ),
    Cohort(
        key="NEOCHEMO",
        gse="GSE207422",
        label="Neoadjuvant chemo-ICI (GSE207422)",
        setting="neoadjuvant",
        regimen="anti-PD-1 + platinum chemotherapy",
        histology_source="pathology",
        outcome="major pathological response",
        supplementary=(
            "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
            "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
        ),
    ),
    Cohort(
        key="PLA",
        gse="GSE283829",
        label="PD1-PD-L1 PLA cohort (GSE283829)",
        setting="advanced",
        regimen="anti-PD-1/PD-L1",
        histology_source="pathology",
        outcome="RECIST best response",
        supplementary=("GSE283829_raw_express_matrix_all_samples.txt.gz",),
    ),
    Cohort(
        key="DURVART",
        gse="GSE253564",
        label="Neoadjuvant durvalumab +/- SBRT (GSE253564)",
        setting="neoadjuvant",
        regimen="durvalumab +/- stereotactic radiation",
        histology_source="pathology",
        outcome="none (immune-context only)",
        supplementary=("GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",),
    ),
    Cohort(
        key="SMC",
        gse="GSE135222",
        label="SMC Korea (GSE135222)",
        setting="advanced",
        regimen="anti-PD-1/PD-L1",
        histology_source="inferred",
        outcome="PFS",
        supplementary=("GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",),
    ),
    Cohort(
        key="YUHS",
        gse="GSE126044",
        label="Yonsei Korea (GSE126044)",
        setting="advanced",
        regimen="anti-PD-1",
        histology_source="inferred",
        outcome="RECIST responder / non-responder",
        supplementary=("GSE126044_counts.txt.gz",),
    ),
]

COHORTS_BY_KEY = {c.key: c for c in COHORTS}

# Cohorts whose squamous label comes from the pathology report recorded in GEO.
PATHOLOGY_COHORTS = [c.key for c in COHORTS if c.histology_source == "pathology"]
INFERRED_COHORTS = [c.key for c in COHORTS if c.histology_source == "inferred"]

# Cohorts that contribute a response/benefit endpoint to the meta-analysis.
RESPONSE_COHORTS = [c.key for c in COHORTS if c.outcome != "none (immune-context only)"]

TCGA_FILES = {
    "lusc_tpm": f"{XENA_GDC}/TCGA-LUSC.star_tpm.tsv.gz",
    "luad_tpm": f"{XENA_GDC}/TCGA-LUAD.star_tpm.tsv.gz",
    "lusc_clinical": f"{XENA_GDC}/TCGA-LUSC.clinical.tsv.gz",
    "probemap": f"{XENA_GDC}/gencode.v36.annotation.gtf.gene.probemap",
}

# Guard rail from the task brief: processed artefacts committed to the repo
# must stay well under 2 GB.
MAX_PROCESSED_BYTES = 2 * 1024**3
# GitHub refuses single files above 100 MB; stay comfortably below.
MAX_SINGLE_FILE_BYTES = 90 * 1024**2
