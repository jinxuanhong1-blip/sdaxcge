"""Locked IFN / MHC / TJ TF lists and companion gene-set modules.

CLDN4 is the splitter and is never a member of a scored regulon or TJ module.
TACSTD2 is a companion gene only — never a gate (no dual-high).
"""
from __future__ import annotations

# --- TF programs (public prior AUCell) ---
IFN_TFS = ["STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IRF9"]
MHC_TFS = ["NLRC5", "CIITA", "RFX5", "RFXANK", "RFXAP"]
TJ_TFS = ["GRHL1", "GRHL2", "ELF3", "KLF4", "OVOL1", "OVOL2"]
CTRL_TFS = ["HIF1A"]
ALL_TFS = sorted(set(IFN_TFS + MHC_TFS + TJ_TFS + CTRL_TFS))

TF_PROGRAM = {}
for g in IFN_TFS:
    TF_PROGRAM[g] = "IFN"
for g in MHC_TFS:
    TF_PROGRAM[g] = "MHC"
for g in TJ_TFS:
    TF_PROGRAM[g] = "TJ"
for g in CTRL_TFS:
    TF_PROGRAM[g] = "CTRL"
# IRF1 sits on both IFN and MHC; keep IFN as the primary bucket, MHC as extra.
TF_PROGRAM["IRF1"] = "IFN"

# --- Companion gene-set modules (AUCell on curated lists, not a TF regulon) ---
# Tight-junction / polarity. CLDN4 intentionally absent.
TJ_NO_CLDN4 = [
    "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1", "TJP2", "F11R",
    "PARD3", "MARVELD2", "CGN", "CRB3", "JAM3",
]
KERATIN_SIMPLE = ["KRT7", "KRT8", "KRT18", "KRT19"]
KERATIN_BASAL = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"]

# Type-I ISG core (same list as the public CLDN4/IFN alignment).
ISG_CORE = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]
USER_CORE6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

MHC1_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
    "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1",
]

CTRL_OXPHOS = [
    "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
    "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1",
]

GENESET_MODULES = {
    "gs_tj_no_cldn4": TJ_NO_CLDN4,
    "gs_ifn_isg": ISG_CORE,
    "gs_ifn_core6": USER_CORE6,
    "gs_mhc1_apm": MHC1_APM,
    "gs_oxphos": CTRL_OXPHOS,
    "gs_keratin_simple": KERATIN_SIMPLE,
}

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}

NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "AGER", "NAPSA",
    "SCGB1A1", "SCGB3A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS",
]

FOCAL = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "MKI67"]
HOLD_OUT = {"CLDN4"}  # never inside a scored regulon / TJ set


def all_always_genes() -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    blocks = [
        FOCAL,
        ALL_TFS,
        TJ_NO_CLDN4,
        KERATIN_SIMPLE,
        KERATIN_BASAL,
        ISG_CORE,
        USER_CORE6,
        MHC1_APM,
        CTRL_OXPHOS,
        NORMAL_LUNG,
        *[LINEAGE_MARKERS[k] for k in LINEAGE_MARKERS],
    ]
    for block in blocks:
        for g in block:
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out
