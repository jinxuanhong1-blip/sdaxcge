"""Locked gene sets for the GSE148071 CLDN4-primary PAGA slice.

CLDN4 is the readout. It is NOT in the barrier/keratin score (avoids
circularity). TACSTD2 is a comparator only — this is not a TACSTD2 redo.
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
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

# TISCH2 major-lineage tokens (strip trailing spaces at load).
TISCH_MALIGNANT = {"Malignant"}
TISCH_LEFTOVER_EPI = {"Alveolar", "Basal", "Epithelial", "AT1", "AT2", "Ciliated", "Club"}
TISCH_EPITHELIAL = TISCH_MALIGNANT | TISCH_LEFTOVER_EPI
