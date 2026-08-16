#!/usr/bin/env python3
"""Fixed gene lists for C4 GSE22493 IFN / MHC-I / APM.

Priority six genes are the user-named panel. APM is a predefined MHC-I /
antigen-presentation machinery panel. IFN_IMMUNE is the broader type-I
ISG + IFN signaling + chemokine set used in the parallel CLDN4-loss
slices. Lists are human HGNC symbols.
"""

from __future__ import annotations

# User-named priority panel (C4 request: IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A,
# with TAP1/TAP2 called out alongside MHC-I/APM).
PRIORITY = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# MHC-I / antigen-presentation machinery (classic APM + editors).
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

# Broader IFN / ISG / IFN-signaling / chemokine set (human).
IFN_IMMUNE = [
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

# Knockdown sanity-check gene (not part of the IFN claim).
PERTURBATION = ["CLDN4"]
