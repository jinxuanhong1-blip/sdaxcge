"""Locked gene sets for the GSE131907+GSE189357 CLDN4 trajectory.

CLDN4 is the readout. It is NOT in the barrier/keratin score and not in IFN.
TACSTD2 is a comparator only — no dual-high gate.
IFN primary = Hallmark IFNα ∪ IFNγ mean (CLDN4 not a member).
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
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
}

# Compact IFN core used as a secondary score (also no CLDN4).
IFN_CORE = (
    "STAT1",
    "IRF1",
    "IRF9",
    "ISG15",
    "MX1",
    "OAS1",
    "IFI27",
    "IFI44",
    "IFI44L",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "B2M",
    "TAP1",
    "PSMB8",
    "PSMB9",
    "HLA-A",
    "HLA-B",
    "HLA-C",
)

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

AUTHOR_AT2 = {"AT2"}
AUTHOR_CLUB = {"Club"}
AUTHOR_TUMOR_STATE = {
    "tS1",
    "tS2",
    "tS3",
    "Malignant cells",
}


def hallmark_ifn_union(path: Path | None = None) -> tuple[str, ...]:
    if path is None:
        path = Path(__file__).resolve().parents[1] / "data" / "ifn_sets.json"
    raw = json.loads(path.read_text())
    genes = set(raw["HALLMARK_INTERFERON_ALPHA_RESPONSE"]) | set(
        raw["HALLMARK_INTERFERON_GAMMA_RESPONSE"]
    )
    genes.discard("CLDN4")
    return tuple(sorted(genes))
