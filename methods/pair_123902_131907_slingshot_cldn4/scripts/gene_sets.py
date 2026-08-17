"""Locked gene sets for the GSE123902+GSE131907 CLDN4-only trajectory.

CLDN4 is the readout. It is NOT in the barrier/keratin score.
TACSTD2 is a comparator only — no dual-high gate.
IFN is a compact ISG / IFNG-response core (not the 200-gene Hallmark dump).
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
    # Compact IFN / ISG core. CLDN4 is not in this list.
    "IFN": (
        "STAT1",
        "STAT2",
        "IRF1",
        "IRF7",
        "IRF9",
        "ISG15",
        "ISG20",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "IFITM1",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "OAS3",
        "OASL",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "IDO1",
        "GBP1",
        "GBP2",
        "IFI27",
        "IFI44",
        "IFI44L",
        "IFI6",
        "IFI16",
        "RSAD2",
        "SAMD9",
        "SAMD9L",
        "BST2",
        "B2M",
        "TAP1",
        "PSMB8",
        "PSMB9",
    ),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

AUTHOR_AT2 = {"AT2"}
AUTHOR_CLUB = {"Club"}
AUTHOR_TUMOR_STATE = {
    "tS1",
    "tS2",
    "tS3",
    "Malignant cells",
}

EPI_MARKERS_123902 = ("EPCAM", "KRT8", "KRT18", "KRT19")
