#!/usr/bin/env python3
"""Fixed gene lists for the C-CPE / CLDN4 pharmacologic analog.

Lists are human HGNC symbols, frozen before looking at fold-changes.
They match the C4 GSE22493 IFN/MHC-I/APM slice plus a pre-specified
tight-junction panel used in public TROP2-ADC analog work.

TACSTD2 is requested. GPL10555 does not annotate it (only TACSTD1 /
EPCAM). The analysis must report that absence rather than invent a
value.
"""

from __future__ import annotations

# User-named priority IFN / MHC-I panel (same six genes as C4 GSE22493).
PRIORITY = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# MHC-I / antigen-presentation machinery.
APM = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "NLRC5",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "ERAP1",
    "ERAP2",
    "CALR",
    "CANX",
    "PDIA3",
]

# Broader type-I ISG + IFN signaling + chemokine set (same as C4).
IFN_IMMUNE = [
    "ISG15",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFIT5",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "IFI6",
    "IFI27",
    "IFI35",
    "IFI44",
    "IFI44L",
    "IFI16",
    "DDX58",
    "IFIH1",
    "DDX60",
    "XAF1",
    "HERC5",
    "USP18",
    "CMPK2",
    "BST2",
    "GBP1",
    "GBP2",
    "GBP4",
    "GBP5",
    "EIF2AK2",
    "ZBP1",
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF7",
    "IRF9",
    "JAK2",
    "SOCS1",
    "SOCS3",
    "B2M",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "NLRC5",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "TAPBP",
    "CXCL10",
    "CXCL11",
    "CXCL9",
    "CCL5",
    "CCL2",
    "IL6",
    "TNF",
    "NFKB1",
    "RELA",
    "IL15",
    "IL32",
    "IFNB1",
    "IFNG",
    "IFNAR1",
    "IFNAR2",
    "IFNGR1",
    "IFNGR2",
    "IFNL1",
]

# Tight-junction / epithelial barrier. C-CPE binds CLDN3/4.
TJ = [
    "CLDN1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN5",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CLDN10",
    "CLDN12",
    "CLDN18",
    "TJP1",
    "TJP2",
    "TJP3",
    "OCLN",
    "MARVELD2",
    "MARVELD3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CDH1",
    "CTNNB1",
    "CTNNA1",
    "EPCAM",
    "CRB3",
]

# Axis genes requested in the analog brief.
AXIS = ["CLDN4", "CLDN3", "TACSTD2"]

# Historical / alias remaps used only when ORF or DESCRIPTION prefix
# is the old symbol. Labelled in tables; not a free-text salvage.
ALIASES = {
    "G1P2": "ISG15",  # IFI-15K
    "IFI15": "ISG15",
    "G1P3": "IFI6",
    "TROP2": "TACSTD2",
    "GA733-1": "TACSTD2",
    "M1S1": "TACSTD2",
    "EGP1": "TACSTD2",
}

# Hallmark sets used for optional preranked GSEA (descriptive).
HALLMARK_FOCUS = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
]
