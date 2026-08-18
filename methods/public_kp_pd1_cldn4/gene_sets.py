"""Cldn4-only mouse gene sets. T/NK and IFN/MHC are correlates, not the claim."""

CLDN4_SYMBOL = "Cldn4"
CLDN4_ENS = "ENSMUSG00000047501"

# Same T/NK list used on the public mouse KL hunt / leftover ICI pages.
TNK_SYMBOLS = [
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Cd2",
    "Cd8a",
    "Cd8b1",
    "Cd4",
    "Nkg7",
    "Gzma",
    "Gzmb",
    "Prf1",
    "Klrb1c",
    "Ncr1",
    "Klrd1",
    "Klrc1",
    "Ifng",
]

# Compact scRNA gate (leftover mouse ICI page).
SCRNA_TNK_GATE = ["Cd3d", "Cd3e", "Cd8a", "Nkg7", "Ncr1"]
SCRNA_TSCORE = ["Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Nkg7", "Gzmb", "Prf1", "Ifng"]
