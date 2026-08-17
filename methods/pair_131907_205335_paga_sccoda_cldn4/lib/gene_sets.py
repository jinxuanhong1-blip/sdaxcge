"""Locked gene sets. CLDN4 is the readout; it is not in barrier/keratin."""

STATES = {
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
        "CDKN1A",
        "PLAUR",
    ),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
    "tnk_cytotoxic": ("GZMB", "PRF1", "NKG7", "GNLY"),
    "tnk_exh": ("HAVCR2", "LAG3", "TIGIT", "PDCD1", "TOX"),
    "tnk_naive": ("CCR7", "SELL", "TCF7", "IL7R"),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
