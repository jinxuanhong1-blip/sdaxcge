"""Pre-specified gene sets for extra public TROP2-ADC RNA (not GSE312098).

Frozen before scoring extra matrices. IMMU132 = sacituzumab govitecan, not SKB264.
"""

# Compact IFN / MHC-I panel used on the given CX-1 analog.
C4_IFN_MHCI = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

MHC_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5"]

APM = [
    "B2M",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "NLRC5",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "HLA-G",
    "CALR",
    "CANX",
    "ERAP1",
    "ERAP2",
    "PDIA3",
]

IFN = [
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
    "JAK1",
    "JAK2",
    "SOCS1",
    "SOCS3",
    "CXCL10",
    "CXCL11",
    "CXCL9",
    "CCL5",
    "IFNB1",
    "IFNG",
    "IFNAR1",
    "IFNAR2",
    "IFNGR1",
    "IFNGR2",
    "IFNL1",
]

# KEGG hsa04530 tight-junction core + common epithelial barrier genes.
# CLDN4 is the primary gene; the set is not a substitute for it.
KEGG_TJ = [
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
    "CRB3",
]

TARGETS = ["TACSTD2", "CLDN4"]


def all_sets():
    return {
        "C4_IFN_MHCI": C4_IFN_MHCI,
        "MHC_I": MHC_I,
        "APM": APM,
        "IFN": IFN,
        "KEGG_TJ": KEGG_TJ,
    }
