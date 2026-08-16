"""Gene signatures used in the SCLC ICI hunt.

All lists are transcribed from the primary literature. Where a published score
uses proprietary / weighted coefficients (e.g. the Ayers T-cell-inflamed GEP or
the Zhang neuroendocrine score), we DO NOT claim to reproduce the exact weighted
score. Instead we compute a transparent, unweighted mean z-score over the gene
set and label it as such. This is stated openly in the report so results are not
over-sold as identical to the original publications.
"""

# --- SCLC subtype-defining transcription factors (Rudin et al. 2019) ---
SUBTYPE_TFS = {
    "ASCL1": "SCLC-A",
    "NEUROD1": "SCLC-N",
    "POU2F3": "SCLC-P",
}
YAP1 = "YAP1"  # historically proposed SCLC-Y; not an independent subtype (Baine 2020)

# --- Genes of interest: ADC / surface antigen targets ---
# TACSTD2 = TROP2 (target of sacituzumab govitecan, datopotamab deruxtecan)
# CLDN4   = Claudin-4 (claudin-targeting ADC / bispecific interest)
GENES_OF_INTEREST = ["TACSTD2", "CLDN4"]

# Context ADC / surface targets discussed by Gay et al. 2021 for SCLC subtypes
CONTEXT_SURFACE_TARGETS = ["DLL3", "CEACAM5", "SSTR2", "MICA", "EPCAM"]

# --- Ayers et al. 2017 (JCI) 18-gene IFN-gamma T-cell-inflamed GEP ---
# Used by Gay et al. 2021 (their Figure 3F/G) to define the SCLC-I immune axis.
# We compute an UNWEIGHTED mean z-score over these genes (see module docstring).
AYERS_TCELL_INFLAMED_GEP = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]

# --- Rooney et al. 2015 cytolytic activity (CYT) score ---
# Defined exactly as geometric mean of GZMA and PRF1 expression (TPM/FPKM).
CYT_GENES = ["GZMA", "PRF1"]

# --- Immune checkpoints highlighted as higher in SCLC-I (Gay 2021) ---
IMMUNE_CHECKPOINTS = [
    "CD274", "PDCD1", "CTLA4", "CD80", "CD86", "CD38", "IDO1", "TIGIT",
    "VSIR", "ICOS", "LAG3", "HAVCR2",  # VSIR = C10orf54/VISTA
]

# --- Antigen presentation / HLA (higher in SCLC-I) ---
HLA_ANTIGEN_PRESENTATION = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
    "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5",
]

# STING-induced T-cell attractant chemokines (Gay 2021, Sup Fig 3O-P)
STING_CHEMOKINES = ["CCL5", "CXCL10", "CXCL9", "CXCL11"]

# --- Neuroendocrine (NE) axis, transparent proxy ---
# Canonical NE and non-NE markers. We deliberately EXCLUDE the genes of interest
# (TACSTD2, CLDN4) from the non-NE side to avoid circularity, since TACSTD2 in
# particular appears in some published non-NE signatures.
NE_UP = [
    "ASCL1", "INSM1", "CHGA", "CHGB", "SYP", "NCAM1", "GRP", "CALCA",
    "SCG2", "SCG3", "TAGLN3", "SEZ6", "KIF1A", "MYT1", "CELF3", "BEX1",
]
NON_NE_UP = [
    "YAP1", "VIM", "CAV1", "CAV2", "ANXA1", "S100A10", "EPHA2", "TGFBI",
    "AHNAK", "MYOF", "PLAU", "ITGB4", "CCND1", "LGALS3",
]  # note: TACSTD2 / CLDN4 intentionally omitted

# --- EMT axis (SCLC-I is described as the most mesenchymal subtype) ---
EMT_MESENCHYMAL = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2"]
EPITHELIAL = ["CDH1", "EPCAM", "CLDN3", "CLDN7", "KRT8", "KRT18", "KRT19"]

# Convenience registry of "mean z-score" signatures for scoring
MEANZ_SIGNATURES = {
    "GEP_Ayers_Tcell_inflamed": AYERS_TCELL_INFLAMED_GEP,
    "HLA_antigen_presentation": HLA_ANTIGEN_PRESENTATION,
    "Immune_checkpoints": IMMUNE_CHECKPOINTS,
    "STING_chemokines": STING_CHEMOKINES,
    "EMT_mesenchymal": EMT_MESENCHYMAL,
    "Epithelial": EPITHELIAL,
}
