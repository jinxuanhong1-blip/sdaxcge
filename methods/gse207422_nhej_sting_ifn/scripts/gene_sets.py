"""Gene sets for the GSE207422 NHEJ / STING / IFN split.

The three primary modules do not share genes.

NHEJ is classical non-homologous end joining (Ku, DNA-PKcs, ligation).
MRN and HR genes are excluded so the score is not generic DNA repair
and does not overlap cytosolic-DNA sensing.

STING is the proximal cGAS–STING machinery (sensor, adaptor, kinase, IRF3).
ISGs are excluded. Gene symbols that have changed (CGAS/MB21D1,
STING1/TMEM173) are both requested; scoring keeps the symbol present
in the matrix.

IFN effectors are the prior 40-gene type-I ISG core with the sensor and
signaling genes removed (DDX58, IFIH1, IRF7, IRF9, STAT1, STAT2). The
unsplit 40-gene list is retained only so this slice can be compared
with the earlier GSE207422 IFN result.

TACSTD2 and CLDN4 are the presence gate and the splitters. Neither
gene is a member of any scored module.
"""

from __future__ import annotations

# Classical NHEJ. No MRN, no HR, no ISGs.
NHEJ = [
    "XRCC6",
    "XRCC5",
    "PRKDC",
    "LIG4",
    "XRCC4",
    "NHEJ1",
    "DCLRE1C",
    "PAXX",
    "PNKP",
    "APLF",
    "POLL",
    "POLM",
]

# Proximal cGAS–STING. Aliases are resolved at score time.
STING_ALIASES = {
    "CGAS": ("CGAS", "MB21D1"),
    "STING1": ("STING1", "TMEM173"),
}
STING = ["CGAS", "STING1", "TBK1", "IKBKE", "IRF3"]

# Negative regulators of cytosolic DNA / cGAMP. Scored separately.
# Not part of the STING activation module (opposite direction).
STING_NEGATIVE = ["TREX1", "ENPP1"]

# RIG-I-like RNA sensing. Not STING. Kept so a STING result is not
# read as generic innate sensing.
RNA_SENSING = ["DDX58", "IFIH1", "MAVS", "DHX58"]

# Prior 40-gene type-I ISG core (align_cldn4_ifn / GSE207422 TJ-IFN slice).
ISG_CORE40 = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]

# Signaling / sensor genes removed from the effector module.
IFN_SIGNAL_EXCLUDED = ["DDX58", "IFIH1", "IRF7", "IRF9", "STAT1", "STAT2"]

IFN_EFFECTOR = [g for g in ISG_CORE40 if g not in set(IFN_SIGNAL_EXCLUDED)]

# Specificity controls. Not primary endpoints.
HR = ["RAD51", "BRCA1", "BRCA2", "PALB2", "RAD51C", "XRCC2", "XRCC3", "RPA1"]
OXPHOS = [
    "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
    "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1",
]
PROLIF = ["MKI67", "TOP2A", "PCNA", "MCM2", "TYMS"]

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
    "mast": ["KIT"],
}

A3_NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]

FOCAL = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC"]

PRIMARY_MODULES = ("nhej", "sting", "ifn_effector")

# STING machinery without IRF3. IRF3 is widely detected and is not
# STING-restricted, so this is a sensitivity, not the primary module.
STING_NO_IRF3 = ["CGAS", "STING1", "TBK1", "IKBKE"]

MODULES = {
    "nhej": NHEJ,
    "sting": STING,
    "sting_no_irf3": STING_NO_IRF3,
    "ifn_effector": IFN_EFFECTOR,
    "ifn_isg40": ISG_CORE40,
    "sting_negative": STING_NEGATIVE,
    "rna_sensing": RNA_SENSING,
    "hr": HR,
    "oxphos": OXPHOS,
    "prolif": PROLIF,
}


def requested_symbols() -> list[str]:
    """Every symbol to pull from the UMI, including STING aliases."""
    genes: list[str] = []
    seen: set[str] = set()

    def add(g: str) -> None:
        if g not in seen:
            seen.add(g)
            genes.append(g)

    for g in FOCAL + A3_NORMAL_LUNG:
        add(g)
    for block in MODULES.values():
        for g in block:
            if g in STING_ALIASES:
                for alias in STING_ALIASES[g]:
                    add(alias)
            else:
                add(g)
    for geneset in LINEAGES.values():
        for g in geneset:
            add(g)
    return genes
