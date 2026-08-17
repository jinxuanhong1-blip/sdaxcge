"""Locked gene sets for GSE205335 malignant CLDN4-only Slingshot/PAGA.

CLDN4 is the readout. It is NOT in the barrier/keratin score.
IFN is a compact Hallmark-like ISG panel (no MSigDB download).
TACSTD2 is a comparator only — no dual-high gate.
"""

STATES = {
    "AT2": ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"),
    "AT1": ("AGER", "PDPN", "CAV1"),
    "club": ("SCGB1A1", "SCGB3A2", "SCGB3A1"),
    "basal": ("KRT5", "KRT15", "TP63", "NGFR"),
    "ciliated": ("FOXJ1", "TPPP3", "PIFO"),
    # CLDN4 excluded on purpose.
    "barrier_keratin": (
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT7",
        "CDKN1A",
        "PLAUR",
    ),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
    # Compact Hallmark IFNα/γ-style ISGs used in the pair IFN DE folders.
    "IFN": (
        "STAT1",
        "IRF1",
        "ISG15",
        "MX1",
        "OAS1",
        "IFIT1",
        "IFIT3",
        "IFI6",
        "BST2",
        "GBP1",
        "CXCL9",
        "CXCL10",
        "IFI27",
        "RSAD2",
        "IFITM1",
        "OAS2",
        "IRF7",
        "IFI44L",
        "MX2",
        "ISG20",
    ),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")
