"""Canonical GRCh38 Ensembl gene IDs used in this slice.

Targets are the master-task IDs. Signature IDs were resolved against
Ensembl REST /lookup/symbol/homo_sapiens on 2026-08-16.
"""

TARGETS = {
    "TACSTD2": "ENSG00000184292",  # TROP2 / UniProt P09758
    "CLDN4": "ENSG00000189143",
}

# Ayers 2017 JCI 18-gene T-cell inflamed GEP
GEP_TCELL_INFLAMED = {
    "CCL5": "ENSG00000271503",
    "CD27": "ENSG00000139193",
    "CD274": "ENSG00000120217",
    "CD276": "ENSG00000103855",
    "CD8A": "ENSG00000153563",
    "CMKLR1": "ENSG00000174600",
    "CXCL9": "ENSG00000138755",
    "CXCR6": "ENSG00000172215",
    "HLA-DQA1": "ENSG00000196735",
    "HLA-DRB1": "ENSG00000196126",
    "HLA-E": "ENSG00000204592",
    "IDO1": "ENSG00000131203",
    "LAG3": "ENSG00000089692",
    "NKG7": "ENSG00000105374",
    "PDCD1LG2": "ENSG00000197646",
    "PSMB10": "ENSG00000205220",
    "STAT1": "ENSG00000115415",
    "TIGIT": "ENSG00000181847",
}

IFNG_6GENE = {
    "IFNG": "ENSG00000111537",
    "STAT1": "ENSG00000115415",
    "IDO1": "ENSG00000131203",
    "CXCL9": "ENSG00000138755",
    "CXCL10": "ENSG00000169245",
    "HLA-DRA": "ENSG00000204287",
}

CD8_TCELL = {
    "CD8A": "ENSG00000153563",
    "CD8B": "ENSG00000172116",
    "GZMK": "ENSG00000113088",
}

CYTOTOXIC_EFFECTOR = {
    "GZMA": "ENSG00000145649",
    "GZMB": "ENSG00000100453",
    "PRF1": "ENSG00000180644",
    "NKG7": "ENSG00000105374",
    "GNLY": "ENSG00000115523",
    "KLRD1": "ENSG00000134539",
}

CYT_ROONEY = {
    "GZMA": "ENSG00000145649",
    "PRF1": "ENSG00000180644",
}

TLS_12CHEMOKINE = {
    "CCL2": "ENSG00000108691",
    "CCL3": "ENSG00000277632",
    "CCL4": "ENSG00000275302",
    "CCL5": "ENSG00000271503",
    "CCL8": "ENSG00000108700",
    "CCL18": "ENSG00000275385",
    "CCL19": "ENSG00000172724",
    "CCL21": "ENSG00000137077",
    "CXCL9": "ENSG00000138755",
    "CXCL10": "ENSG00000169245",
    "CXCL11": "ENSG00000169248",
    "CXCL13": "ENSG00000156234",
}

TGFB_EXCLUSION = {
    "TGFB1": "ENSG00000105329",
    "TGFB2": "ENSG00000092969",
    "TGFB3": "ENSG00000119699",
    "TGFBR1": "ENSG00000106799",
    "TGFBR2": "ENSG00000163513",
    "LTBP1": "ENSG00000049323",
}

IMMUNE_GENERAL = {
    "PTPRC": "ENSG00000081237",
    "CD3D": "ENSG00000167286",
    "CD3E": "ENSG00000198851",
    "CD2": "ENSG00000116824",
    "CD4": "ENSG00000010610",
    "CD8A": "ENSG00000153563",
    "CD19": "ENSG00000177455",
    "MS4A1": "ENSG00000156738",
}

SIGNATURES = {
    "GEP_Tcell_inflamed": GEP_TCELL_INFLAMED,
    "IFNG_6gene": IFNG_6GENE,
    "CD8_Tcell": CD8_TCELL,
    "Cytotoxic_effector": CYTOTOXIC_EFFECTOR,
    "CYT_cytolytic": CYT_ROONEY,
    "TLS_12chemokine": TLS_12CHEMOKINE,
    "TGFB_exclusion": TGFB_EXCLUSION,
    "Immune_general": IMMUNE_GENERAL,
}

# Immune / IO-relevant columns from LSCC_phenotype.txt (freeze-derived).
PHENOTYPE_IMMUNE = [
    "CIBERSORT_B_cell_naive",
    "CIBERSORT_B_cell_memory",
    "CIBERSORT_B_cell_plasma",
    "CIBERSORT_T_cell_CD8+",
    "CIBERSORT_T_cell_CD4+_naive",
    "CIBERSORT_T_cell_CD4+_memory_resting",
    "CIBERSORT_T_cell_CD4+_memory_activated",
    "CIBERSORT_T_cell_follicular_helper",
    "CIBERSORT_T_cell_regulatory_(Tregs)",
    "CIBERSORT_T_cell_gamma_delta",
    "CIBERSORT_NK_cell_resting",
    "CIBERSORT_NK_cell_activated",
    "CIBERSORT_Monocyte",
    "CIBERSORT_Macrophage_M0",
    "CIBERSORT_Macrophage_M1",
    "CIBERSORT_Macrophage_M2",
    "CIBERSORT_Myeloid_dendritic_cell_resting",
    "CIBERSORT_Myeloid_dendritic_cell_activated",
    "CIBERSORT_Mast_cell_activated",
    "CIBERSORT_Mast_cell_resting",
    "CIBERSORT_Eosinophil",
    "CIBERSORT_Neutrophil",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore",
    "xCell_T_cell_CD8+",
    "xCell_T_cell_CD8+_effector_memory",
    "xCell_T_cell_CD4+_(non-regulatory)",
    "xCell_T_cell_regulatory_(Tregs)",
    "xCell_Macrophage_M1",
    "xCell_Macrophage_M2",
    "xCell_NK_cell",
    "xCell_B_cell",
    "xCell_immune_score",
    "xCell_stroma_score",
    "xCell_microenvironment_score",
    "xCell_Cancer_associated_fibroblast",
    "xCell_Myeloid_dendritic_cell",
    "PROGENy_JAK-STAT",
    "PROGENy_NFkB",
    "PROGENy_TGFb",
    "PROGENy_TNFa",
    "PROGENy_Trail",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "TMB",
]
