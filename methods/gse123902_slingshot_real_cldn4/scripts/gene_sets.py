"""Locked gene sets for GSE123902 CLDN4-only Slingshot/PAGA.

CLDN4 is the readout. It is NOT in the barrier/keratin score.
TACSTD2 is a comparator only — no dual-high gate.
IFN is Hallmark IFNα ∩ IFNγ core (CLDN4 is not in this set).
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
    "IFN": (
        "STAT1",
        "STAT2",
        "IRF1",
        "IRF7",
        "IRF9",
        "ISG15",
        "ISG20",
        "MX1",
        "MX2",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "OAS1",
        "OAS2",
        "OASL",
        "IFI27",
        "IFI44",
        "IFI44L",
        "IFI6",
        "RSAD2",
        "CXCL9",
        "CXCL10",
        "CXCL11",
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
EPI_GATE = ("EPCAM", "KRT8", "KRT18", "KRT19")
