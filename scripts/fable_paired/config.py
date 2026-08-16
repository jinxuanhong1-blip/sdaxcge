"""Shared configuration for the fable_paired ICI lung TACSTD2/CLDN4 analysis.

All heavy raw downloads live outside the repo (default: /tmp/fable_raw) so that
committed, "processed" artifacts stay comfortably under the 2 GB budget. Only the
compact tables and figures written under results/fable_paired/ are version
controlled.
"""
from __future__ import annotations

import os
from pathlib import Path

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "fable_paired"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
NOTES_DIR = REPO_ROOT / "notes" / "fable_paired"

# Raw downloads (NOT committed). Override with FABLE_RAW_DIR if desired.
RAW_DIR = Path(os.environ.get("FABLE_RAW_DIR", "/tmp/fable_raw"))

for _d in (RESULTS_DIR, TABLES_DIR, FIGURES_DIR, NOTES_DIR, RAW_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------
# Genes of interest
# ----------------------------------------------------------------------------
TARGET_GENES = ["TACSTD2", "CLDN4"]  # TROP2 and Claudin-4

# Markers used to classify cells in scRNA data.
EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
IMMUNE_MARKERS = ["PTPRC", "CD3D", "CD68", "LYZ"]
CONTEXT_GENES = ["CLDN3", "CLDN7", "KRT17", "SFTPC"]

# Full symbol panel extracted from the scRNA matrix in a single streaming pass.
SCRNA_PANEL = sorted(set(TARGET_GENES + EPI_MARKERS + IMMUNE_MARKERS + CONTEXT_GENES))

# Cross-namespace identifiers for the target genes.
ENSEMBL_IDS = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}
ENTREZ_IDS = {"TACSTD2": "4070", "CLDN4": "1364"}

# ----------------------------------------------------------------------------
# Dataset registry: accession -> (relative GEO path, filename)
# ----------------------------------------------------------------------------
GEO_FILES = {
    # GSE207422 - NSCLC neoadjuvant anti-PD-1 + chemo (whole transcriptome)
    "GSE207422_bulk_expr": (
        "series/GSE207nnn/GSE207422/suppl",
        "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    ),
    "GSE207422_bulk_meta": (
        "series/GSE207nnn/GSE207422/suppl",
        "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
    ),
    "GSE207422_sc_expr": (
        "series/GSE207nnn/GSE207422/suppl",
        "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
    ),
    "GSE207422_sc_meta": (
        "series/GSE207nnn/GSE207422/suppl",
        "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
    ),
    # GSE126044 - NSCLC anti-PD-1 baseline (counts, responder/non-responder)
    "GSE126044_counts": (
        "series/GSE126nnn/GSE126044/suppl",
        "GSE126044_counts.txt.gz",
    ),
    "GSE126044_matrix": (
        "series/GSE126nnn/GSE126044/matrix",
        "GSE126044_series_matrix.txt.gz",
    ),
    # GSE135222 - NSCLC anti-PD-(L)1 baseline (TPM, PFS)
    "GSE135222_tpm": (
        "series/GSE135nnn/GSE135222/suppl",
        "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    ),
    "GSE135222_matrix": (
        "series/GSE135nnn/GSE135222/matrix",
        "GSE135222_series_matrix.txt.gz",
    ),
    # GSE91061 - Riaz melanoma, paired Pre/On nivolumab (FPKM, Entrez)
    "GSE91061_fpkm": (
        "series/GSE91nnn/GSE91061/suppl",
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
    ),
    "GSE91061_matrix": (
        "series/GSE91nnn/GSE91061/matrix",
        "GSE91061_series_matrix.txt.gz",
    ),
    # GSE248249 - Memon NSCLC paired pre / acquired-resistance (Clariom D)
    "GSE248249_matrix": (
        "series/GSE248nnn/GSE248249/matrix",
        "GSE248249_series_matrix.txt.gz",
    ),
    # GSE246922 - mouse KP / LLC1 ICB-resistant RNA (VST)
    "GSE246922_KP": (
        "series/GSE246nnn/GSE246922/suppl",
        "GSE246922_KP_RNA_counts_vst.csv.gz",
    ),
    "GSE246922_LLC1": (
        "series/GSE246nnn/GSE246922/suppl",
        "GSE246922_LLC1_RNA_counts_vst.csv.gz",
    ),
}

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"


def raw_path(key: str) -> Path:
    """Return the local path to a registered raw file."""
    _, fname = GEO_FILES[key]
    return RAW_DIR / fname
