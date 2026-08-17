"""Pre-registered gene sets and cohort cell-type rules (CLDN4-only).

Cytotoxicity / IFN / exhaustion lists are a priori (not mined from these
matrices). NicheNet-v2 prior scores are unsigned regulatory potential —
direction of a T/NK program is tested at the patient level.

TACSTD2 is extracted as a companion column and is never used to define
senders or patient groups.
"""

from __future__ import annotations

# --- a priori T/NK programs ---

CYTOTOXICITY = [
    "GZMB",
    "GZMA",
    "GZMH",
    "GZMK",
    "PRF1",
    "GNLY",
    "NKG7",
    "IFNG",
    "FASLG",
    "TNF",
    "CST7",
    "FGFBP2",
    "KLRK1",
    "GZMM",
]

# Compact type-II / shared ISG program. Not Hallmark IFNG-response
# (too large for a panel-restricted NicheNet background).
IFN = [
    "IFNG",
    "STAT1",
    "IRF1",
    "IRF7",
    "ISG15",
    "MX1",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "OAS1",
    "IFI6",
    "RSAD2",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "GBP1",
]

EXHAUSTION = [
    "PDCD1",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "TOX",
    "CTLA4",
    "ENTPD1",
    "LAYN",
    "CXCL13",
    "CD38",
    "CD244",
    "TOX2",
    "CD160",
    "BTLA",
]

LINEAGE = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B",
    "CD4",
    "FOXP3",
    "NKG7",
    "GNLY",
    "KLRD1",
    "NCR1",
]

# Extra T/NK genes for empirical DE background (patient-level Wilcoxon).
EXTRA_TNK = [
    "IL2",
    "IL2RA",
    "IL7R",
    "TCF7",
    "LEF1",
    "CCR7",
    "SELL",
    "CXCR5",
    "CXCR6",
    "CX3CR1",
    "ITGAE",
    "ZNF683",
    "EOMES",
    "TBX21",
    "PRDM1",
    "ID2",
    "ID3",
    "BATF",
    "IRF4",
    "IRF9",
    "STAT2",
    "NR4A1",
    "NR4A2",
    "NR4A3",
    "TOX",
    "MKI67",
    "TOP2A",
    "CD27",
    "CD28",
    "ICOS",
    "TNFRSF9",
    "TNFRSF18",
    "CD69",
    "HLA-DRA",
    "CD74",
    "OAS2",
    "OAS3",
    "MX2",
    "IFI27",
    "IFI44",
    "IFI44L",
    "ISG20",
    "GBP2",
    "GBP5",
    "CCL4",
    "CCL5",
    "XCL1",
    "XCL2",
    "IL10",
    "TGFB1",
    "IFNG",
    "TNF",
]

# --- GSE207422 (Hu / DRMref) ---
GSE207422_MALIGNANT = "Malignant cells"
GSE207422_TNK = ("CD8+ T cells", "CD4+ T cells", "NK cells")
MPR_LABELS = {"MPR", "pCR", "MPR (pCR)"}
NMPR_LABELS = {"NMPR"}

# --- GSE131907 (Kim et al.) ---
GSE131907_TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
GSE131907_MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
GSE131907_TNK_TYPES = {"T lymphocytes", "NK cells"}

# --- GSE205335 (Hu / Ahn / Lee) ---
GSE205335_T_SUBS = {"CD4+ T cells", "CD8+ T cells"}
GSE205335_NK_SUBS = {"NK cells"}

# Inclusion gates (cells are counts; patient is the unit)
MIN_MALIG_PER_STATE = 10
MIN_TNK = 20
MIN_MALIG_PROGRAM = 20  # between-patient program / Q4Q1 row
DETECT_FRAC = 0.10
THIN_Q_N = 8
THIN_TAIL = 3
