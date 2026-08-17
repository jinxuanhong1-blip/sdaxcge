"""Locked gene sets for REAL Slingshot/PAGA on GSE148071, CLDN4 only.

CLDN4 is the readout. It is NOT in the barrier/keratin score and NOT in
the IFN score. TACSTD2 is a comparator only — no dual-high gate.

IFN = Hallmark IFNα ∩ IFNγ intersection (Liberzon Cell Syst 2015;
MSigDB Hallmark 2020 via the locked a8 freeze). Compact enough for
single-cell mean scoring; not a 200-gene IFNγ dump.
"""

STATES = {
    "AT2": ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"),
    "AT1": ("AGER", "PDPN", "CAV1"),
    "club": ("SCGB1A1", "SCGB3A2", "SCGB3A1"),
    "basal": ("KRT5", "KRT15", "TP63", "NGFR"),
    "ciliated": ("FOXJ1", "TPPP3", "PIFO"),
    # CLDN4 excluded on purpose (same lock as PR #394).
    "barrier_keratin": (
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT7",
        "CDKN1A",
        "PLAUR",
    ),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
    # Hallmark IFNA ∩ IFNG (73). WARS1/WARS both accepted at score time.
    "IFN": (
        "ADAR",
        "B2M",
        "BATF2",
        "BST2",
        "C1S",
        "CASP1",
        "CASP8",
        "CD74",
        "CMPK2",
        "CMTR1",
        "CXCL10",
        "CXCL11",
        "DDX60",
        "DHX58",
        "EIF2AK2",
        "EPSTI1",
        "GBP4",
        "HELZ2",
        "HERC6",
        "IFI27",
        "IFI30",
        "IFI35",
        "IFI44",
        "IFI44L",
        "IFIH1",
        "IFIT2",
        "IFIT3",
        "IFITM2",
        "IFITM3",
        "IL15",
        "IL4R",
        "IL7",
        "IRF1",
        "IRF2",
        "IRF7",
        "IRF9",
        "ISG15",
        "ISG20",
        "LAP3",
        "LGALS3BP",
        "LY6E",
        "MX1",
        "NMI",
        "OASL",
        "OGFR",
        "PARP12",
        "PARP14",
        "PLSCR1",
        "PNPT1",
        "PSMA3",
        "PSMB8",
        "PSMB9",
        "PSME1",
        "PSME2",
        "RIPK2",
        "RNF31",
        "RSAD2",
        "RTP4",
        "SAMD9L",
        "SLC25A28",
        "SP110",
        "STAT2",
        "TAP1",
        "TDRD7",
        "TRAFD1",
        "TRIM14",
        "TRIM21",
        "TRIM25",
        "TRIM26",
        "TXNIP",
        "UBE2L6",
        "USP18",
        "WARS1",
    ),
}

IFN_ALIASES = {"WARS1": ("WARS",)}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

TISCH_MALIGNANT = {"Malignant"}
TISCH_LEFTOVER_EPI = {"Alveolar", "Basal", "Epithelial", "AT1", "AT2", "Ciliated", "Club"}
TISCH_EPITHELIAL = TISCH_MALIGNANT | TISCH_LEFTOVER_EPI

# PR #394 given result — do not re-audit.
PR394_BARRIER_RHO = 0.590
PR394_BARRIER_N = 42
PR394_BARRIER_P = 3.87e-05
