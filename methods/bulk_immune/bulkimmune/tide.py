"""TIDE (Jiang et al., Nature Medicine 2018) via the official ``tidepy`` package.

TIDE is a *trained* model, not a gene-set score. The official implementation
(``tidepy``, MIT licence) ships its own weights (``model.pkl``) and gene-ID
table, so the scores computed here are the published TIDE scores, not a
reimplementation. What this wrapper adds is the input-scale contract and the
NSCLC-specific caveats that the playbook relies on.

Input contract (from ``tidepy.pred.TIDE``):
    genes x samples, HGNC / Ensembl / Entrez. If the matrix is not already
    row-centred log2 expression, tidepy applies ``log2(x+1)`` and subtracts
    each gene's mean across the samples in hand. That last step makes TIDE
    **cohort-relative**: a sample's Exclusion score depends on who else is in
    the matrix. Never concatenate two cohorts and then split the scores.

Cancer argument:
    ``"NSCLC"`` averages the LUAD and LUSC signature SDs. ``"Melanoma"`` uses
    the SKCM model. ``"Other"`` falls back to melanoma. Use ``NSCLC`` here.

``pretreat``:
    If the patient has already received ICI, TIDE uses Exclusion for every
    sample. For treatment-naive biopsies (GSE126044, GSE135222, TCGA) leave
    this False so that CTL-high samples are scored on the Dysfunction axis.

The web server at tide.dfci.harvard.edu is the same model; this wrapper is
offline-equivalent and should be preferred for a reproducible methods section.
"""

from __future__ import annotations

from typing import Literal

import pandas as pd

__all__ = ["tide_score", "CANCER_CHOICES"]

CANCER_CHOICES = ("NSCLC", "Melanoma", "Other")


def tide_score(
    expr: pd.DataFrame,
    cancer: Literal["NSCLC", "Melanoma", "Other"] = "NSCLC",
    pretreat: bool = False,
    ignore_norm: bool = False,
    force_normalize: bool = False,
) -> pd.DataFrame:
    """Run official TIDE. Returns samples x scores.

    Columns: No benefits, Responder, TIDE, IFNG, MSI Score, CD274, CD8,
    CTL.flag, Dysfunction, Exclusion, MDSC, CAF, TAM M2, CTL.
    """
    try:
        from tidepy.pred import TIDE
    except ImportError as exc:
        raise ImportError(
            "tidepy is required for TIDE scores. Install with: pip install tidepy"
        ) from exc

    if cancer not in CANCER_CHOICES:
        raise ValueError(f"cancer must be one of {CANCER_CHOICES}")

    result = TIDE(
        expression=expr.copy(),
        cancer=cancer,
        pretreat=pretreat,
        ignore_norm=ignore_norm,
        force_normalize=force_normalize,
    )
    # tidepy sorts by TIDE descending; restore input sample order.
    result = result.reindex(expr.columns)
    result.attrs["cancer"] = cancer
    result.attrs["pretreat"] = pretreat
    return result
