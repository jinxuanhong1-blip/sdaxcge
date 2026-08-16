#!/usr/bin/env python3
"""Pre-specified human gene sets for GSE302284 (and reusable TROP2 / TJ / IFN scoring).

SKB264 / CLDN4 thesis is taken as given. These lists are fixed before looking at
GSE302284 values. Symbols are HGNC.

Tight-junction (TJ) signature: CLDN1/4/7, F11R, PARD3, OCLN, TJP1
  (user-specified composite; mean of per-cell log1p(CP10k) of genes present).

IFN / MHC-I / APM: same human panels used in the C4 public-analog slices
  (PRIORITY6, MHC1_APM, IFN_ISG).
"""

from __future__ import annotations

# Targets
TARGETS = ["TACSTD2", "CLDN4"]

# User-specified tight-junction signature
TJ_SIG = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3", "OCLN", "TJP1"]

# Priority IFN / MHC-I genes (C4 CORE6)
PRIORITY6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# MHC-I antigen-processing / presentation machinery
MHC1_APM = [
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

# Broader type-I ISG + IFN signaling + IFN-inducible chemokines
IFN_ISG = [
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
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CCL5",
    "IFNB1",
    "IFNG",
    "IFNAR1",
    "IFNAR2",
    "IFNGR1",
    "IFNGR2",
]

# Epithelial / immune markers used only for patient-sample compartmenting
EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"]
IMMUNE_MARKERS = ["PTPRC"]

# Extra genes reported individually (not part of set tests)
REPORT_EXTRA = ["CXCL9", "CXCL10", "EPCAM", "CD274"]

SETS = {
    "TJ_SIG": TJ_SIG,
    "PRIORITY6": PRIORITY6,
    "MHC1_APM": MHC1_APM,
    "IFN_ISG": IFN_ISG,
}

FOCAL_GENES = list(
    dict.fromkeys(TARGETS + TJ_SIG + PRIORITY6 + REPORT_EXTRA)
)
