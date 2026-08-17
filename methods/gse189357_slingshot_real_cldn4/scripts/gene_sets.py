"""Locked gene sets for GSE189357 CLDN4-only REAL Slingshot/PAGA.

CLDN4 is the readout. It is NOT in the barrier/keratin score.
TACSTD2 is a comparator only — no dual-high gate.
IFN is Hallmark IFNα ∪ IFNγ (CLDN4 is not an IFN gene).
"""

from __future__ import annotations

import json
from pathlib import Path

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
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67", "UBE2C"),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC", "UBE2C")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")
EPI_POS = ("EPCAM", "KRT8", "KRT18", "KRT19")

STAGE_ORDER = ("AIS", "MIA", "IAC")


def load_ifn() -> dict[str, tuple[str, ...]]:
    path = Path(__file__).with_name("ifn_genes.json")
    blob = json.loads(path.read_text())
    return {
        "IFN": tuple(blob["hallmark_ifn_union"]),
        "IFN_core": tuple(blob["ifn_core"]),
    }
