"""Pre-specified gene sets for the C4 analog (GSE245459 SKOV3 shTACSTD2).

C4 private claim: CLDN4 KD opens IFN/MHC-I (IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A).
This analog asks whether TACSTD2 KD in SKOV3 moves CLDN4 and the same IFN/APM
programs. Sets are frozen before looking at fold-changes.
"""

# Perturbation target and the C4 junction gene.
TARGETS = ["TACSTD2", "CLDN4"]

# Exact C4 claim panel (IFN/MHC-I "opening").
C4_IFN_MHCI = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# Antigen-presentation machinery (APM / MHC-I pathway).
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

# Broader type-I/II interferon / ISG core (excludes APM overlap except STAT/IRF).
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

# Tight-junction / epithelial barrier (CLDN4 is the C4 gene; kept as internal).
JUNCTION = [
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

# Paper-claimed TACSTD2 → Rap1/PI3K/AKT axis (negative-control / context, not C4).
RAP1_PI3K_AKT = [
    "RAP1A",
    "RAP1B",
    "RAPGEF3",
    "PIK3CA",
    "PIK3CB",
    "PIK3R1",
    "AKT1",
    "AKT2",
    "AKT3",
    "MTOR",
    "PTEN",
    "GSK3B",
    "FOXO1",
    "FOXO3",
]


def all_sets():
    return {
        "C4_IFN_MHCI": C4_IFN_MHCI,
        "APM": APM,
        "IFN": IFN,
        "JUNCTION": JUNCTION,
        "RAP1_PI3K_AKT": RAP1_PI3K_AKT,
    }
