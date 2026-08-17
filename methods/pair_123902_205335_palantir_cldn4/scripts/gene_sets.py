"""Locked gene sets. CLDN4 is the readout and is held out of barrier.

IFN = Hallmark IFNα ∪ IFNγ from the pair IFN DE a8_sets.json (already −1.05).
No TACSTD2∩CLDN4 dual-high gate.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
_IFN = json.loads((HERE / "data" / "ifn_sets.json").read_text())

STATES = {
    "AT2": ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"),
    "AT1": ("AGER", "PDPN", "CAV1"),
    "club": ("SCGB1A1", "SCGB3A2", "SCGB3A1"),
    "basal": ("KRT5", "KRT15", "TP63", "NGFR"),
    "ciliated": ("FOXJ1", "TPPP3", "PIFO"),
    # CLDN4 excluded on purpose.
    "barrier": ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR"),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")
EPI_POS = ("EPCAM", "KRT8", "KRT18", "KRT19")

IFN = tuple(_IFN["IFN"])
MHC_I_APM = tuple(_IFN["MHC_I_APM"])
