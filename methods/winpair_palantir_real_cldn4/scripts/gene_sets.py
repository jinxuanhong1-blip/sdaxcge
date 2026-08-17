"""Locked gene sets for winning-pair REAL Palantir, CLDN4 only.

CLDN4 is the readout. It is NOT in the barrier score and NOT in the IFN score.
TACSTD2 is a comparator only — no dual-high gate.
IFN core excludes LAMP3 (AT2 overlap) and CDKN1A (barrier overlap).
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
    # Compact IFN / ISG core. No CLDN4. No LAMP3. No CDKN1A.
    "IFN": (
        "STAT1",
        "IRF1",
        "IRF7",
        "IRF9",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "OAS1",
        "OAS2",
        "MX1",
        "MX2",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "IDO1",
        "TAP1",
        "GBP1",
        "IFI44L",
        "IFI27",
        "IFI44",
        "RSAD2",
        "USP18",
        "EPSTI1",
        "SAMD9",
    ),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

AUTHOR_AT2 = {"AT2"}
AUTHOR_CLUB = {"Club"}
AUTHOR_AIRWAY = {"Club", "Ciliated", "Basal"}
AUTHOR_TUMOR_STATE = {
    "tS1",
    "tS2",
    "tS3",
    "Malignant cells",
    "Malignant",
    "Cancer cells",
    "Tumor",
}
