"""Gene sets for GSE207422 Palantir destinies (CLDN4-only).

CLDN4 is the readout, never a terminal definition and never a member of the
barrier module. TACSTD2 is a companion only — never a dual-high gate.
"""

from __future__ import annotations

# Tight-junction / polarity. CLDN4 intentionally absent (circular otherwise).
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

KERATIN_SIMPLE = ["KRT7", "KRT8", "KRT18", "KRT19"]
KERATIN_BASAL = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"]
KERATIN_ALL = KERATIN_SIMPLE + KERATIN_BASAL

# Combined barrier/keratin (primary barrier score). No CLDN4.
BARRIER_NO_CLDN4 = TJ_NO_CLDN4 + KERATIN_ALL

# Compact transitional / KAC-like waypoint (CLDN4 still excluded).
TRANSITIONAL_NO_CLDN4 = ["KRT8", "KRT19", "CDKN1A", "PLAUR"]

USER_CORE6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

ISG_CORE = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]

MHC1_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
    "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1",
]

CTRL_OXPHOS = [
    "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
    "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1",
]

AT2_RESIDUAL = ["SFTPB", "NAPSA", "LAMP3", "ABCA3", "SFTPA1", "SFTPC"]

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

FOCAL = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC", "MKI67", "CDKN1A", "PLAUR"]


def all_panel_genes() -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    for block in (
        FOCAL,
        TJ_NO_CLDN4,
        KERATIN_ALL,
        TRANSITIONAL_NO_CLDN4,
        USER_CORE6,
        ISG_CORE,
        MHC1_APM,
        CTRL_OXPHOS,
        AT2_RESIDUAL,
        A3_NORMAL_LUNG,
        *[LINEAGES[k] for k in LINEAGES],
    ):
        for g in block:
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


MODULES = {
    "barrier_no_cldn4": BARRIER_NO_CLDN4,
    "tj_no_cldn4": TJ_NO_CLDN4,
    "keratin": KERATIN_ALL,
    "transitional_no_cldn4": TRANSITIONAL_NO_CLDN4,
    "ifn_isg": ISG_CORE,
    "ifn_core6": USER_CORE6,
    "mhc1_apm": MHC1_APM,
    "ctrl_oxphos": CTRL_OXPHOS,
    "at2_residual": AT2_RESIDUAL,
}
