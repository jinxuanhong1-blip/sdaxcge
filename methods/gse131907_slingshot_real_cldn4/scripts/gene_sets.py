"""Locked gene sets for GSE131907 REAL Slingshot, CLDN4-only.

CLDN4 is the readout. It is NOT in the barrier/keratin score and NOT in
the IFN score (no circularity). TACSTD2 is a comparator only — no dual-high
gate. This is not a redo of PR #325 (DPT/PAGA).
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
    # Compact IFN / ISG panel. CLDN4 is not a member.
    "IFN": (
        "STAT1",
        "STAT2",
        "IRF1",
        "IRF7",
        "IRF9",
        "ISG15",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "IFI6",
        "IFI27",
        "IFI44",
        "IFI44L",
        "BST2",
        "GBP1",
        "GBP2",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "IDO1",
        "TAP1",
        "PSMB8",
        "PSMB9",
    ),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC", "STAT1")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

AUTHOR_AT2 = {"AT2"}
AUTHOR_CLUB = {"Club"}
AUTHOR_AT1 = {"AT1"}
AUTHOR_CILIATED = {"Ciliated"}
AUTHOR_TUMOR_STATE = {"tS1", "tS2", "tS3", "Malignant cells"}
AUTHOR_MALIGNANT = {"tS1", "tS2", "tS3", "Malignant cells"}
