"""Pre-registered gene sets and cell-type rules.

Cytotoxicity / exhaustion lists are a priori (not mined from this matrix).
NicheNet prior scores are unsigned regulatory potential — direction of a
T/NK program is tested at the patient level, not assumed from the prior.
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
    "TACSTD2",
    "CLDN4",
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

# Extra T/NK genes often in ICI programs (for empirical DE background)
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
    "STAT1",
    "STAT3",
    "IFIT1",
    "ISG15",
    "MX1",
    "CCL4",
    "CCL5",
    "XCL1",
    "XCL2",
    "IL10",
    "TGFB1",
    "IFNG",
    "TNF",
]
