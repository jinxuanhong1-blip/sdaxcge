"""Gene sets for GSE207422 malignant CLDN4-high vs low TJ/keratin and IFN.

CLDN4 is the splitter and is never a member of the scored TJ module.
Keratin lists match the public CytoTRACE-like slice (simple + basal).
IFN lists match the public CLDN4/IFN alignment core (ISG_CORE, USER_CORE6, MHC1_APM).
"""

from __future__ import annotations

# Tight-junction strand / polarity genes used as a compact barrier module.
# CLDN4 is intentionally absent (circular if scored on a CLDN4 split).
TJ_NO_CLDN4 = [
    "CLDN1",
    "CLDN3",
    "CLDN7",
    "OCLN",
    "TJP1",
    "TJP2",
    "F11R",
    "PARD3",
    "MARVELD2",
    "CGN",
    "CRB3",
    "JAM3",
]

# Simple-epithelium keratins (LUAD-typical)
KERATIN_SIMPLE = ["KRT7", "KRT8", "KRT18", "KRT19"]

# Basal / squamous keratins
KERATIN_BASAL = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"]

KERATIN_ALL = KERATIN_SIMPLE + KERATIN_BASAL

# Combined TJ/keratin module (primary). CLDN4 excluded.
TJ_KERATIN = TJ_NO_CLDN4 + KERATIN_ALL

# User six-gene IFN/MHC-I panel from the public CLDN4/IFN alignment.
USER_CORE6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# Canonical type-I ISG core (same list as align_cldn4_ifn).
ISG_CORE = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]

# MHC class I antigen processing (companion, not the primary IFN module).
MHC1_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
    "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1",
]

# Compact OXPHOS-ish control (specificity: not everything should move).
CTRL_OXPHOS = [
    "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
    "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1",
]

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

# A3 given malignant-like: epithelial AND zero UMI on this panel.
A3_NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]

BROAD_NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER",
    "SCGB1A1", "SCGB3A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS",
]

FOCAL = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC", "MKI67"]


def all_panel_genes() -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    for block in (
        FOCAL,
        TJ_NO_CLDN4,
        KERATIN_ALL,
        USER_CORE6,
        ISG_CORE,
        MHC1_APM,
        CTRL_OXPHOS,
        A3_NORMAL_LUNG,
        BROAD_NORMAL_LUNG,
        *[LINEAGES[k] for k in LINEAGES],
    ):
        for g in block:
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


MODULES = {
    "tj_keratin": TJ_KERATIN,
    "tj_no_cldn4": TJ_NO_CLDN4,
    "keratin": KERATIN_ALL,
    "keratin_simple": KERATIN_SIMPLE,
    "keratin_basal": KERATIN_BASAL,
    "ifn_isg": ISG_CORE,
    "ifn_core6": USER_CORE6,
    "mhc1_apm": MHC1_APM,
    "ctrl_oxphos": CTRL_OXPHOS,
}
