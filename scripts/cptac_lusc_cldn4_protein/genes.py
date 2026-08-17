"""Locked gene sets for CPTAC-LSCC (LUSC) CLDN4 / TJ protein vs immune.

TJ-15 matches the structural 15-gene set used for CPTAC LUAD protein
(PR #245). That LUAD TJ-15 vs ImmuneScore result (ρ=−0.30, n=110) is
taken as given and is not re-estimated here.
"""

from __future__ import annotations

# symbol -> preferred Ensembl gene id, then alternates seen in GENCODE.
ENSEMBL = {
    "CLDN1": ["ENSG00000163347"],
    "CLDN3": ["ENSG00000165215"],
    "CLDN4": ["ENSG00000189143"],
    "CLDN7": ["ENSG00000181885", "ENSG00000288292"],
    "OCLN": ["ENSG00000197822", "ENSG00000273814"],
    "TJP1": ["ENSG00000104067", "ENSG00000277401"],
    "TJP2": ["ENSG00000119139"],
    "TJP3": ["ENSG00000105289"],
    "F11R": ["ENSG00000158769"],
    "JAM2": ["ENSG00000154721"],
    "JAM3": ["ENSG00000166086"],
    "MARVELD2": ["ENSG00000152939", "ENSG00000274671"],
    "MARVELD3": ["ENSG00000140832"],
    "CGN": ["ENSG00000143375"],
    "CGNL1": ["ENSG00000128849"],
    "CD8A": ["ENSG00000153563"],
    "CD8B": ["ENSG00000172116"],
    "GZMA": ["ENSG00000145649"],
    "PRF1": ["ENSG00000180644"],
    "CCL5": ["ENSG00000271503", "ENSG00000274233"],
    "CD27": ["ENSG00000139193"],
    "CD274": ["ENSG00000120217"],
    "CD276": ["ENSG00000103855"],
    "CMKLR1": ["ENSG00000174600"],
    "CXCL9": ["ENSG00000138755"],
    "CXCR6": ["ENSG00000172215"],
    "HLA-DQA1": ["ENSG00000196735", "ENSG00000237541"],
    "HLA-DRB1": ["ENSG00000196126"],
    "HLA-E": ["ENSG00000204592"],
    "IDO1": ["ENSG00000131203"],
    "LAG3": ["ENSG00000089692"],
    "NKG7": ["ENSG00000105374"],
    "PDCD1LG2": ["ENSG00000197646"],
    "PSMB10": ["ENSG00000205220"],
    "STAT1": ["ENSG00000115415"],
    "TIGIT": ["ENSG00000181847"],
    "TACSTD2": ["ENSG00000184292"],
}

# Structural tight-junction set (same 15 genes as PR #245 LUAD protein).
TJ_15 = [
    "CLDN1",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "OCLN",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "MARVELD2",
    "MARVELD3",
    "CGN",
    "CGNL1",
]

# Claim-page 7-gene module. Sensitivity, not a search.
TJ_7 = ["CLDN1", "CLDN4", "CLDN7", "F11R", "TJP1", "TJP2", "OCLN"]

CD8 = ["CD8A", "CD8B"]
CYT = ["GZMA", "PRF1"]
GEP18 = [
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
]

SYNONYMS = {
    "F11R": ["JAM1", "JAMA"],
    "PDCD1LG2": ["PDL2"],
    "HLA-DQA1": ["HLA.DQA1"],
    "HLA-DRB1": ["HLA.DRB1"],
    "HLA-E": ["HLA.E"],
}

PHENO_ALIASES = {
    "ImmuneScore": [
        "ESTIMATE_ImmuneScore",
        "ImmuneScore",
        "ESTIMATE_immune_score",
    ],
    "StromalScore": [
        "ESTIMATE_StromalScore",
        "StromalScore",
        "ESTIMATE_stromal_score",
    ],
    "ESTIMATEScore": [
        "ESTIMATE_ESTIMATEScore",
        "ESTIMATEScore",
        "ESTIMATE_score",
    ],
    "CD8_CIBERSORT": [
        "CIBERSORT_T_cell_CD8+",
        "CIBERSORT_T.cells.CD8",
        "CIBERSORT_T_cells_CD8",
        "CIBERSORT_CD8",
    ],
    "CD8_xCell": [
        "xCell_T_cell_CD8+",
        "xCell_CD8+_T-cells",
        "xCell_CD8_Tcells",
        "xCell_T_cell_CD8",
    ],
    "Immune_xCell": [
        "xCell_immune_score",
        "xCell_ImmuneScore",
    ],
    "WES_purity": [
        "WES_purity",
        "purity_WES",
        "wES_purity",
    ],
    "WGS_purity": [
        "WGS_purity",
        "purity_WGS",
    ],
}
