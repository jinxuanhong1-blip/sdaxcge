"""Pre-specified epithelial state scores for the GSE131907 PAGA slice.

Lists are locked before DPT. Genes missing from the matrix are recorded as
absent and dropped from that score — they are not replaced.
TACSTD2 is intentionally NOT inside the barrier/keratin set.
"""

from __future__ import annotations

STATES: dict[str, tuple[str, ...]] = {
    "AT2": ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"),
    "AT1": ("AGER", "PDPN", "CAV1"),
    "club": ("SCGB1A1", "SCGB3A2", "SCGB3A1"),
    "basal": ("KRT5", "KRT15", "TP63", "NGFR"),
    "ciliated": ("FOXJ1", "TPPP3", "PIFO"),
    "barrier_keratin": (
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT7",
        "CLDN4",
        "CDKN1A",
        "PLAUR",
    ),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
}

FOCAL = ("TACSTD2", "CLDN4")
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

# Author labels used as-is (Kim et al. 2020).
AUTHOR_AT2 = {"AT2"}
AUTHOR_CLUB = {"Club"}
AUTHOR_AT1 = {"AT1"}
AUTHOR_CILIATED = {"Ciliated"}
AUTHOR_TUMOR_STATE = {"tS1", "tS2", "tS3", "Malignant cells"}

AXIS_LABELS = {
    "AT2": "AT2",
    "AT1": "AT1",
    "Club": "club",
    "Ciliated": "ciliated",
    "tS1": "malignant-like (tS1)",
    "tS2": "malignant-like (tS2)",
    "tS3": "malignant-like (tS3)",
    "Malignant cells": "malignant-like",
    "Undetermined": "undetermined",
}
