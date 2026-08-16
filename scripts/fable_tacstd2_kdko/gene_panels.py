"""Curated gene panels for the TROP2/TACSTD2 loss-of-function analysis.

Symbols are given in HUMAN (upper-case) form. For mouse datasets we match
case-insensitively and also map via mygene, so Cldn4/CLDN4 etc. both resolve.

Three biological axes are interrogated:
  1. CLDN4 (headline gene) + broader claudin / tight-junction module
  2. junctions (tight junction, adherens, desmosome)
  3. IFN / immune (interferon-stimulated genes, antigen presentation, chemokines,
     T-cell / cytotoxicity markers relevant in bulk-tumour datasets)
"""

# ---- Headline ----
HEADLINE = ["CLDN4"]

# ---- Tight junction (claudins + scaffold/occludin) ----
TIGHT_JUNCTION = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN6", "CLDN7", "CLDN9", "CLDN10",
    "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN18", "CLDN23",
    "OCLN", "TJP1", "TJP2", "TJP3", "MARVELD2", "MARVELD3",
    "F11R", "JAM2", "JAM3", "CGN", "CGNL1", "CRB3", "PARD3", "PATJ",
]

# ---- Adherens junction ----
ADHERENS = [
    "CDH1", "CDH2", "CTNNB1", "CTNNA1", "CTNND1",
    "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4", "AFDN", "EPCAM",
]

# ---- Desmosome ----
DESMOSOME = [
    "DSP", "DSG1", "DSG2", "DSG3", "DSG4", "DSC1", "DSC2", "DSC3",
    "PKP1", "PKP2", "PKP3", "JUP", "PERP", "PPL", "EVPL",
]

# ---- IFN / interferon-stimulated genes ----
IFN_ISG = [
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20", "IFI6",
    "IFI27", "IFI35", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFITM1",
    "IFITM3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "USP18",
    "DDX58", "IFIH1", "BST2", "XAF1", "HERC5", "HERC6", "SAMD9", "SAMD9L",
    "PARP9", "PARP14", "DTX3L", "EIF2AK2", "IRF8", "STAT2",
]

# ---- IFN-driven chemokines / GBPs ----
IFN_CHEMOKINE = [
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "CXCL13",
    "GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6",
]

# ---- Antigen presentation (human + mouse H2 aliases) ----
ANTIGEN_PRESENTATION = [
    # human
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5", "CIITA",
    "HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DPA1",
    # mouse H2 orthologs
    "H2-K1", "H2-D1", "H2-Q7", "H2-T23", "H2-Aa", "H2-Ab1", "H2-Eb1",
]

# ---- T-cell / cytotoxicity (meaningful in bulk tumour datasets) ----
TCELL_CYTOTOX = [
    "CD8A", "CD8B", "CD3E", "CD3D", "CD3G", "GZMA", "GZMB", "GZMK", "PRF1",
    "IFNG", "NKG7", "KLRG1", "KLRB1", "PDCD1", "CTLA4", "CD274", "PDCD1LG2",
    "FOXP3", "LAG3", "HAVCR2", "TIGIT", "CXCR3",
]

PANELS = {
    "tight_junction": TIGHT_JUNCTION,
    "adherens": ADHERENS,
    "desmosome": DESMOSOME,
    "ifn_isg": IFN_ISG,
    "ifn_chemokine": IFN_CHEMOKINE,
    "antigen_presentation": ANTIGEN_PRESENTATION,
    "tcell_cytotox": TCELL_CYTOTOX,
}

# Higher-level groupings used for reporting the three headline axes.
AXES = {
    "CLDN4_headline": HEADLINE,
    "junctions": sorted(set(TIGHT_JUNCTION + ADHERENS + DESMOSOME)),
    "ifn_immune": sorted(set(IFN_ISG + IFN_CHEMOKINE + ANTIGEN_PRESENTATION + TCELL_CYTOTOX)),
}

# Perturbation gene (for QC that the KD/KO actually worked).
PERTURBATION = ["TACSTD2", "TROP2"]


def all_panel_symbols():
    s = set(HEADLINE)
    for v in PANELS.values():
        s.update(v)
    s.update(PERTURBATION)
    return sorted(s)


def upper(sym):
    return sym.upper().replace("_", "-")
