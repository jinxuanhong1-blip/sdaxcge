"""Gene panels for the NHEJ-loss → IFN/STING/APM bridge.

Symbols are current HGNC names. Older annotations (TMEM173, MB21D1) are
resolved to these names before testing.
"""

# On-target NHEJ components. The perturbed gene must fall.
NHEJ = ["XRCC6", "XRCC5", "PRKDC", "XRCC4", "LIG4", "NHEJ1", "DCLRE1C", "PAXX", "APLF"]

# Cytosolic DNA-sensing / STING induction machinery (the pathway, not the ISG output).
STING = [
    "CGAS", "STING1", "TBK1", "IKBKE", "IRF3", "IRF7", "IFI16", "ZBP1",
    "DDX58", "IFIH1", "MAVS", "IFNB1", "CXCL10",
]

# Type-I ISG core used in the earlier CLDN4-loss panel, plus STAT1.
IFN_CORE = [
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
    "IFIT1", "IFIT2", "IFIT3", "IFI27", "IFI44", "IFI44L", "IFI6",
    "USP18", "HERC5", "DDX58", "IFIH1", "XAF1", "BST2", "ISG20",
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "GBP1", "GBP2",
]

# The eight genes previously treated as the directional IFN/APM core.
CORE8 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A", "TAP1", "TAP2"]

# MHC-I antigen-presentation machinery.
APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M",
    "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10",
    "NLRC5", "ERAP1", "CALR", "PDIA3", "B2M",
]

# Display order for the focus heatmap (unique, current symbols).
FOCUS_DISPLAY = [
    "XRCC6", "XRCC5", "PRKDC", "XRCC4", "LIG4",
    "CGAS", "STING1", "TBK1", "IRF3", "IRF7", "IFNB1",
    "ISG15", "MX1", "IFI27", "IFIT1", "OAS1", "OAS2", "RSAD2", "STAT1", "CXCL10",
    "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5",
]

PANELS = {
    "STING": STING,
    "IFN_CORE": IFN_CORE,
    "CORE8": CORE8,
    "APM": list(dict.fromkeys(APM)),
}
