"""Shared gene panel for neoadjuvant lung scRNA Harmony.

Lineage + target genes that must be kept if present. Extra HVGs are
added at extract time from the intersection of deposited symbols.
"""

from __future__ import annotations

# Dual targets (User A3 taken as given; scored in malignant-like cells only)
TARGETS = ["TACSTD2", "CLDN4"]

# Epithelial / keratin / junction (malignant-like scoring + restriction)
EPITHELIAL = [
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT5",
    "KRT6A",
    "KRT14",
    "KRT17",
    "MUC1",
    "CDH1",
    "ELF3",
    "NKX2-1",
    "NAPSA",
    "SFTPB",
    "TP63",
    "SOX2",
    "CLDN3",
    "CLDN7",
    "CLDN18",
    "TJP1",
    "OCLN",
    "F11R",
]

# Normal-lung programs used to drop AT2 / AT1 / club / ciliated from "malignant-like"
NORMAL_LUNG = [
    "SFTPA1",
    "SFTPA2",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "TPPP3",
    "FOXJ1",
    "CAPS",
    "PIFO",
]

# T / NK
T_MARKERS = ["CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B", "TRAC", "IL7R", "FOXP3", "IL2RA"]
NK_MARKERS = ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "FGFBP2", "KLRF1", "NCR1"]

# Other lineages (Harmony conservation)
B_MARKERS = ["CD79A", "MS4A1", "CD19", "MZB1", "JCHAIN", "IGHG1", "IGHA1", "CD79B"]
MYELOID = [
    "LYZ",
    "CD14",
    "FCGR3A",
    "CD68",
    "C1QA",
    "C1QB",
    "S100A8",
    "S100A9",
    "MARCO",
    "APOE",
    "CSF1R",
    "FCN1",
]
DC_MARKERS = ["CLEC9A", "CD1C", "CLEC10A", "LILRA4", "IRF8", "ITGAX"]
ENDO = ["PECAM1", "VWF", "CLDN5", "ENG", "CDH5", "FLT1"]
FIBRO = ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "PDGFRA", "ACTA2"]
MAST = ["TPSAB1", "TPSB2", "CPA3", "KIT"]
CYCLING = ["MKI67", "TOP2A", "STMN1", "PCNA"]
IMMUNE_PAN = ["PTPRC"]
HOUSEKEEPING = ["ACTB", "GAPDH", "MALAT1", "B2M"]

# Extra commonly variable lung TME genes to thicken the shared space
EXTRA_LUNG = [
    "CXCL9",
    "CXCL10",
    "CXCL13",
    "CCL5",
    "CCL19",
    "GZMA",
    "GZMB",
    "GZMK",
    "PRF1",
    "IFNG",
    "PDCD1",
    "CTLA4",
    "LAG3",
    "HAVCR2",
    "TIGIT",
    "TOX",
    "CD274",
    "PDCD1LG2",
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "STAT1",
    "IRF1",
    "ISG15",
    "MX1",
    "IFI6",
    "IFI27",
    "CD47",
    "EGFR",
    "KRAS",
    "ALK",
    "MET",
    "ERBB2",
    "FGFR1",
    "PIK3CA",
    "TP53",
    "CDKN2A",
    "RB1",
    "MYC",
    "VIM",
    "FN1",
    "ZEB1",
    "SNAI1",
    "TWIST1",
    "CD44",
    "ALDH1A1",
    "SOX9",
    "FOXA2",
    "GATA6",
    "ASCL1",
    "NEUROD1",
    "INSM1",
    "CHGA",
    "SYP",
    "NCAM1",
]


def core_panel() -> list[str]:
    seen: list[str] = []
    for block in (
        TARGETS,
        EPITHELIAL,
        NORMAL_LUNG,
        T_MARKERS,
        NK_MARKERS,
        B_MARKERS,
        MYELOID,
        DC_MARKERS,
        ENDO,
        FIBRO,
        MAST,
        CYCLING,
        IMMUNE_PAN,
        HOUSEKEEPING,
        EXTRA_LUNG,
    ):
        for g in block:
            if g not in seen:
                seen.append(g)
    return seen


LINEAGE_SCORES = {
    "epithelial": EPITHELIAL[:8] + TARGETS,
    "t": T_MARKERS[:6],
    "nk": NK_MARKERS[:5],
    "b": B_MARKERS[:5],
    "myeloid": MYELOID[:8],
    "endothelial": ENDO[:4],
    "fibroblast": FIBRO[:5],
    "mast": MAST[:3],
}
