"""A priori gene sets for the GSE316655 CLDN4-related extract.

Symbols are HGNC. Cell Ranger GRCh38-2020-A may still use older names
(TMEM173, MB21D1). Aliases are resolved at load time; the table records
which symbol was present. Missing genes are dropped from scores, not imputed.

IFN and APM lists match the panel already used in the public CLDN4 IFN/APM
analyses in this repo (ISG_CORE, MHC1_APM). NHEJ is canonical classical NHEJ.
STING is the cGAS–STING DNA-sensing core plus proximal regulators, without
the ISG outputs (those sit in IFN).
"""

from __future__ import annotations

# Classical NHEJ core (Ku–DNA-PKcs–XRCC4–LIG4–XLF–PAXX–Artemis).
NHEJ = [
    "XRCC6",  # Ku70
    "XRCC5",  # Ku80
    "PRKDC",  # DNA-PKcs
    "XRCC4",
    "LIG4",
    "NHEJ1",  # XLF
    "PAXX",
    "DCLRE1C",  # Artemis
]

# cGAS–STING. ISG outputs are intentionally not here.
STING = [
    "CGAS",
    "STING1",
    "TBK1",
    "IKBKE",
    "IRF3",
    "TREX1",
    "ENPP1",
    "IFI16",
    "DDX41",
    "ZBP1",
]

# Type-I ISG panel used by the repo's CLDN4 IFN analyses.
IFN = [
    "ISG15",
    "IFI6",
    "IFI27",
    "IFI44",
    "IFI44L",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFIT5",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "USP18",
    "BST2",
    "XAF1",
    "STAT1",
    "STAT2",
    "IRF7",
    "IRF9",
    "DDX58",
    "IFIH1",
    "SAMD9",
    "SAMD9L",
    "HERC5",
    "HERC6",
    "EPSTI1",
    "CMPK2",
    "PARP9",
    "DTX3L",
    "LY6E",
    "SP100",
    "SP110",
    "PLSCR1",
]

# MHC-I antigen processing and presentation.
APM = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "B2M",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "PSME1",
    "PSME2",
    "NLRC5",
    "ERAP1",
    "ERAP2",
    "CALR",
    "PDIA3",
    "CANX",
    "IRF1",
]

PROGRAMS = {
    "NHEJ": NHEJ,
    "STING": STING,
    "IFN": IFN,
    "APM": APM,
}

# Older symbols that GRCh38-2020-A may still carry.
ALIASES = {
    "CGAS": ["CGAS", "MB21D1", "C6orf150"],
    "STING1": ["STING1", "TMEM173"],
    "NHEJ1": ["NHEJ1", "XLF"],
    "DCLRE1C": ["DCLRE1C", "ARTEMIS"],
    "PAXX": ["PAXX", "C9orf142"],
    "DDX58": ["DDX58", "RIGI"],
}

LINEAGE = {
    "myeloid": ["LYZ", "CD14", "CD68", "CSF1R", "FCGR3A", "S100A8", "S100A9", "CST3"],
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "FGFBP2"],
    "B": ["CD79A", "MS4A1", "CD19", "IGHM"],
    "pDC": ["IL3RA", "CLEC4C", "LILRA4"],
    "melanoma": ["MLANA", "PMEL", "TYR", "MITF", "DCT"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
}

ANCHORS = ["CLDN4", "CLDN18", "LILRB2", "PTPRC"]


def lookup_names(symbol: str) -> list[str]:
    return ALIASES.get(symbol, [symbol])


def panel_symbols() -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for block in (NHEJ, STING, IFN, APM, ANCHORS):
        for g in block:
            if g not in seen:
                seen.add(g)
                out.append(g)
    for genes in LINEAGE.values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out
