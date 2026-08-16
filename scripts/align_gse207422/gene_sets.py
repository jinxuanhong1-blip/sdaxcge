"""Literature gene sets used for the GSE207422 USER-ALIGN analysis.

These are starting panels, not a claim that the user's unpublished slide used
exactly these symbols. Each set is cited in WRITEUP.md.
"""

# Rooney et al. Cell 2015 cytolytic activity (CYT).
CYT = ["GZMA", "PRF1"]

# Broader cytotoxic effector program (CD8 / NK shared).
CYTOTOXIC = [
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7",
    "KLRD1", "KLRK1", "IFNG", "FASLG",
]

# T-cell exhaustion / co-inhibitory checkpoints.
EXHAUSTION = [
    "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "BTLA",
]

# Tight-junction / apical-junction program (user TJ intersection + core TJ).
TIGHT_JUNCTION = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN8",
    "TJP1", "TJP2", "OCLN", "F11R", "PARD3",
    "CGN", "MARVELD2", "CRB3", "CDH1", "EPCAM",
]

# Keratinization / squamous differentiation (GO keratinization core + lung SQ markers).
KERATINIZATION = [
    "KRT1", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17",
    "SPRR1A", "SPRR1B", "SPRR2A", "SPRR2D", "SPRR3",
    "IVL", "LORICRIN", "TGM1", "TGM3", "DSG1", "DSC1", "PKP1",
    "EVPL", "CSTA", "CNFN", "S100A7", "S100A8", "S100A9",
]

# Lineage markers used only for annotation sanity checks.
LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "SFN"],
    "T": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA"],
    "B_Plasma": ["MS4A1", "CD79A", "MZB1"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
}
