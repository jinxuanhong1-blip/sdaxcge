"""Pre-registered gene sets for CLDN4-only MultiNicheNet (winning pair).

Cytotoxicity / IFN lists are a priori (not mined from these matrices).
NicheNet-v2 prior scores are unsigned regulatory potential.

This folder is CLDN4-only. TACSTD2 is never used to define senders.
GSE207422 is not used (PR #334 NS).
"""

# GSE131907 (Kim et al. Nat Commun 2020) author labels
KIM_TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain", "PE")
KIM_MALIGNANT_SUBTYPES = ("Malignant cells", "tS1", "tS2", "tS3")
KIM_TNK_TYPES = ("T lymphocytes", "NK cells")

# GSE205335 (Hu / Ahn / Lee) author labels
HU_MALIGNANT_SUB = "Malignant cells"
HU_TNK_TOTAL = "T/NK cells"

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

# Sensitivity only
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
