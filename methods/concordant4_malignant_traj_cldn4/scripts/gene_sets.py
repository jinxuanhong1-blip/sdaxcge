"""Gene sets for concordant-4 malignant AT2 → barrier trajectories.

CLDN4 is the readout along the path. It is not in the barrier score and not
in the IFN score. TACSTD2 is not a gate.
"""

AT2 = ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3")
BARRIER = ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR")  # CLDN4 out
IFN = (
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
)
MARKER_EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")
FOCAL = ("CLDN4",)
CONTROL = ("SFTPC",)
