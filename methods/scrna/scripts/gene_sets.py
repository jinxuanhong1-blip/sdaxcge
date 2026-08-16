"""Curated gene sets for the lung-cancer ICI scRNA playbook (Python).

Edit / extend to taste and CITE the source when you use one in a manuscript.
These are starting points, not gospel. See playbook.md sections 6-7.
"""

# --- Tumor epithelial program of interest ---------------------------------
TACSTD2_CLDN4_JUNCTION = [
    "TACSTD2", "CLDN4", "CLDN3", "CLDN7", "CLDN18",
    "CDH1", "TJP1", "OCLN", "F11R", "EPCAM", "ELF3",
]

EPITHELIAL_LINEAGE = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "SFTPC", "SCGB1A1"]

# --- CD8 T-cell states ------------------------------------------------------
CD8_CYTOTOXIC = ["CD8A", "CD8B", "GZMB", "GZMK", "GZMH", "PRF1", "IFNG", "NKG7", "KLRG1"]
CD8_EXHAUSTION = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1"]
TRM = ["ITGAE", "ZNF683", "CXCR6", "ITGA1"]

# --- TLS / follicular -------------------------------------------------------
# 12-chemokine TLS signature (adapt / cite Coppola/Prabhakaran-style panels).
TLS_CHEMOKINE = ["CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19",
                 "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"]
TLS_CORE = ["CXCL13", "CCL19", "CCL21", "CR2", "CXCR5", "LTB", "SELL", "MS4A1"]

CXCL13 = ["CXCL13"]

# --- Lineage markers for coarse annotation ---------------------------------
LINEAGE_MARKERS = {
    "T/NK":        ["CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma":    ["MS4A1", "CD79A", "CD79B", "BANK1", "MZB1", "IGHG1", "XBP1"],
    "Myeloid":     ["LYZ", "CD68", "CD14", "FCGR3A", "C1QA", "C1QB", "SPP1"],
    "DC":          ["CLEC9A", "CD1C", "LAMP3", "LILRA4"],
    "Mast":        ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Epithelial":  EPITHELIAL_LINEAGE,
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CLEC14A"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "ACTA2"],
}

ALL_MODULES = {
    "TACSTD2_CLDN4_junction": TACSTD2_CLDN4_JUNCTION,
    "CD8_cytotoxic": CD8_CYTOTOXIC,
    "CD8_exhaustion": CD8_EXHAUSTION,
    "TRM": TRM,
    "TLS_core": TLS_CORE,
    "TLS_chemokine": TLS_CHEMOKINE,
    "CXCL13": CXCL13,
}
