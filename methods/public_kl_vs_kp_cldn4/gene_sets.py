"""Cldn4-only mouse gene sets for public KL vs KP/K contrasts.

Do not expand Cldn4 to other claudins. Do not define Tacstd2∩Cldn4 dual-high.
"""

CLDN4 = ["Cldn4"]

# Cytotoxic / T and NK lineage transcripts (mouse symbols).
T_NK = [
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

# Compact epithelial IFN response.
IFN = [
    "Stat1",
    "Stat2",
    "Irf1",
    "Irf7",
    "Irf9",
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Mx1",
    "Oasl2",
    "Rsad2",
    "Ifih1",
    "Ddx58",
    "Ifnb1",
]

# Compact MHC / APM (mouse).
MHC = [
    "B2m",
    "H2-K1",
    "H2-D1",
    "H2-Q4",
    "H2-Q6",
    "H2-Q7",
    "H2-Aa",
    "H2-Ab1",
    "H2-Eb1",
    "Tap1",
    "Tap2",
    "Psmb8",
    "Psmb9",
    "Nlrc5",
    "Ciita",
]

IFN_MHC = IFN + MHC

# Marker genes used only to classify leftover compartments in GSE267321
# (non-malignant digest). Not a Cldn4 gate.
T_MARKERS = ["Cd3d", "Cd3e", "Cd3g", "Cd8a"]
NK_MARKERS = ["Nkg7", "Ncr1", "Klrb1c"]
EPI_MARKERS = ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"]
HOST_LUNG = ["Sftpc", "Scgb1a1", "Ager"]
