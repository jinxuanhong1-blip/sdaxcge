"""Pre-registered gene sets and cell-type rules (CLDN4-only senders).

Cytotoxicity / IFN lists are a priori (not mined from this matrix).
NicheNet prior scores are unsigned regulatory potential — direction of a
T/NK program is tested at the patient level, not assumed from the prior.

This folder is CLDN4-only. TACSTD2 is extracted for a companion column
and is never used to define the sender set.
"""

# DRMref labels used as senders / receivers
MALIGNANT = "Malignant cells"
TNK_TYPES = ("CD8+ T cells", "CD4+ T cells", "NK cells")
CD8_TYPE = "CD8+ T cells"
NK_TYPE = "NK cells"

# Hu et al. Genome Med 2023: pCR is a major pathologic response
MPR_LABELS = {"MPR", "pCR", "MPR (pCR)"}
NMPR_LABELS = {"NMPR"}

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

# Compact a priori T/NK IFN / ISG program (type-II + shared ISGs).
# Not Hallmark IFNG-response (too large for a panel-restricted background).
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

# Sensitivity only — not a primary estimand in this folder
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

# Extra T/NK genes often in ICI / IFN programs (empirical DE background)
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
