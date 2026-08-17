"""Locked gene sets. CLDN4 is the readout and is excluded from barrier/keratin.

No TACSTD2∩CLDN4 dual-high gate. TACSTD2 is a comparator only.
"""

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "EPCAM", "PTPRC", "MKI67")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

# CLDN4 is intentionally absent.
BARRIER_KERATIN = ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR")
MALIGNANT_LIKE = ("CEACAM5", "CEACAM6", "MKI67")
AT2 = ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3")
