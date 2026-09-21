"""Locked protein list. Not a search.

TROP2 is TACSTD2. The mediator under test is CLDN4. EPCAM is the epithelial
surface control, run through the same equations. Outcomes are CD8A protein and
the MHC-I score already used on this freeze (mean of per-gene z-scores of
HLA-A, HLA-B, and HLA-C). B2M is recorded for coverage and is not added to
that score.
"""

from __future__ import annotations

ENSEMBL: dict[str, list[str]] = {
    "TACSTD2": ["ENSG00000184292"],
    "CLDN4": ["ENSG00000189143"],
    "EPCAM": ["ENSG00000119888"],
    "CD8A": ["ENSG00000153563"],
    "HLA-A": ["ENSG00000206503"],
    "HLA-B": ["ENSG00000234745"],
    "HLA-C": ["ENSG00000204525"],
    "B2M": ["ENSG00000166710"],
}

MHC1_GENES = ["HLA-A", "HLA-B", "HLA-C"]
REQUIRED = ["TACSTD2", "CLDN4", "EPCAM", "CD8A", *MHC1_GENES]
OUTCOMES = ["CD8A", "MHC1"]
EXPOSURE = "TACSTD2"
MEDIATOR = "CLDN4"
CONTROL = "EPCAM"
