"""CLDN4-only mouse gene sets for T/NK and epithelial IFN/MHC scores."""

# Claudin-4 only (do not expand to other claudins).
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

# Epithelial / tumor IFN response + MHC antigen presentation.
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
