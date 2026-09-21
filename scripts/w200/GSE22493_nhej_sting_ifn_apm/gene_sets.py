#!/usr/bin/env python3
"""Pre-specified panels for the GSE22493 CLDN4 ovarian array.

Lists were fixed before reading log-ratios. NHEJ is classical
non-homologous end joining, not homologous recombination. STING is the
cGAS–STING axis, not the interferon program. IFN_IMMUNE and APM are the
same human lists used in the earlier C4 slice of this accession, so this
wave does not swap in a friendlier interferon set.

Symbols are current HGNC. Historical ORF symbols on GPL10555 are mapped
in ALIASES inside the analysis script, not by dropping genes from these lists.
"""

from __future__ import annotations

# Classical NHEJ core machinery.
NHEJ = [
    "XRCC6",  # Ku70
    "XRCC5",  # Ku80
    "PRKDC",  # DNA-PKcs
    "LIG4",
    "XRCC4",
    "NHEJ1",  # XLF
    "DCLRE1C",  # Artemis
    "PAXX",
]

# End-processing / accessory factors often listed beside NHEJ.
# Secondary: not mixed into the primary NHEJ set test.
NHEJ_ACCESSORY = [
    "PNKP",
    "APLF",
    "APTX",
    "POLL",
    "POLM",
    "DNTT",
    "WRN",
]

# cGAS–STING axis. STING1 and TBK1 stay on the list even if the 2010
# platform does not annotate them.
STING = [
    "CGAS",
    "STING1",
    "TBK1",
    "IRF3",
    "IKBKE",
]

# Cytosolic-DNA regulators discussed next to STING. Secondary.
STING_REG = [
    "TREX1",
    "ENPP1",
    "IFI16",
    "DDX41",
    "SAMHD1",
    "RNASEH2A",
    "RNASEH2B",
    "RNASEH2C",
]

# MHC-I / antigen-presentation machinery. Same 16-gene list as C4.
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

# Broader IFN / ISG / IFN-signaling set. Same 73-gene list as C4.
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

PERTURBATION = ["CLDN4"]

PRIMARY = {
    "NHEJ": NHEJ,
    "STING": STING,
    "IFN": IFN,
    "APM": APM,
}

SECONDARY = {
    "NHEJ_ACCESSORY": NHEJ_ACCESSORY,
    "STING_REG": STING_REG,
}
