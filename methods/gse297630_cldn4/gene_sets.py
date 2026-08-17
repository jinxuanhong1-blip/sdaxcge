"""Pre-specified IFN / MHC panels for the GSE297630 public matrix.

The deposited Clariom S processed table uses human gene symbols (HLA-A, not
H2-K1). Lists match that column so scoring stays on the same matrix that
gives Tacstd2 +2.05.
"""

# Given / additive targets
TARGETS = ["TACSTD2", "CLDN4"]

# Six-gene IFN/MHC-I panel used across the CLDN4-loss public slices
CORE6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# Broader type-I ISG + IFN signaling + chemokine set (human symbols)
IFN_ISG = [
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
    "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3",
    "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "IFI16",
    "DDX58", "IFIH1", "DDX60", "XAF1", "HERC5", "USP18", "CMPK2", "BST2",
    "GBP1", "GBP2", "GBP4", "GBP5", "EIF2AK2", "ZBP1",
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "JAK2", "SOCS1", "SOCS3",
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "CCL2",
    "IFNB1", "IFNG", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
]

# Classical MHC-I + B2M (the clean MHC call on this matrix)
MHC_I = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M"]

# MHC-II + master regulator
MHC_II = [
    "CIITA", "CD74",
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
]

# Antigen-presentation machinery / immunoproteasome (not classical MHC-I)
APM = [
    "NLRC5", "TAP1", "TAP2", "TAPBP",
    "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3",
]

# Union used as "IFN/MHC" in FINDING.md
MHC_APM = MHC_I + MHC_II + APM
IFN_MHC = IFN_ISG + MHC_APM

SETS = {
    "CORE6": CORE6,
    "IFN_ISG": IFN_ISG,
    "MHC_I": MHC_I,
    "MHC_II": MHC_II,
    "APM": APM,
    "MHC_APM": MHC_APM,
    "IFN_MHC": IFN_MHC,
}
