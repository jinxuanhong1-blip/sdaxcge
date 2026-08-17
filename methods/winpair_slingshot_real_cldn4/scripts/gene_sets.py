"""Locked gene sets for the winning-pair CLDN4 real Slingshot trajectory.

CLDN4 is the readout. It is NOT in the barrier/keratin or IFN scores.
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
    # Compact Hallmark IFN-γ / ISG set. No CLDN4. No TACSTD2.
    "IFN": (
        "STAT1",
        "IRF1",
        "IRF7",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "OAS1",
        "OAS2",
        "MX1",
        "IFI27",
        "IFI44",
        "IFI44L",
        "RSAD2",
        "IFI6",
        "GBP1",
        "IDO1",
        "CXCL9",
        "CXCL10",
        "CXCL11",
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
