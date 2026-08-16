"""ESTIMATE (Yoshihara et al., Nature Communications 2013).

Faithful reimplementation of ``estimate::estimateScore`` (R-Forge estimate
1.0.13, GPL-2). Two properties of the original that are easy to get wrong and
are preserved here:

1. **The background matters.** The original first calls ``filterCommonGenes``,
   restricting the matrix to the 10,412 genes common to the platforms the model
   was built on. The enrichment score is a rank walk over whatever background it
   is given, so scores computed on all ~20,000 genes are *not* on the same scale
   as published ESTIMATE scores.
2. **The score is an unnormalised ssGSEA integral** (sum of the running
   difference of the two CDFs), with alpha = 0.25 on rank-normalised values.
   Because ``Pn`` is divided by the sum of the weights, ESTIMATE's
   ``10000*rank/n_genes`` rescaling cancels exactly, so the statistic is
   identical to ``GSVA::gsva(method="ssgsea", ssgsea.norm=FALSE)`` up to tie
   handling. That equivalence is exercised in ``tests/test_ssgsea.py``.

Tumour purity is the Affymetrix-calibrated transform
``cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)``. The original package only
emits it for ``platform="affymetrix"``; it is *not* calibrated for RNA-seq, so
``tumor_purity`` here must be requested explicitly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .ssgsea import ssgsea

__all__ = ["load_signatures", "load_common_genes", "filter_common_genes", "estimate_score",
           "PURITY_INTERCEPT", "PURITY_SLOPE"]

PURITY_INTERCEPT = 0.6049872018
PURITY_SLOPE = 0.0001467884


def load_signatures(path: str | Path) -> dict[str, list[str]]:
    """Load the stromal/immune signatures (GMT written by 00_fetch_resources)."""
    sets: dict[str, list[str]] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def load_common_genes(path: str | Path) -> list[str]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    column = "GeneSymbol" if "GeneSymbol" in table.columns else table.columns[1]
    return [g for g in table[column].dropna().unique() if g]


def filter_common_genes(
    expr: pd.DataFrame, common_genes: list[str], verbose: bool = True
) -> pd.DataFrame:
    """Restrict to the ESTIMATE common-gene background (10,412 genes)."""
    keep = [g for g in common_genes if g in expr.index]
    if verbose:
        print(
            f"[estimate] merged dataset includes {len(keep)} genes "
            f"({len(common_genes) - len(keep)} common genes missing)"
        )
    if len(keep) < 0.5 * len(common_genes):
        raise ValueError(
            f"only {len(keep)}/{len(common_genes)} ESTIMATE common genes found; "
            "check the identifier namespace before trusting these scores"
        )
    return expr.loc[keep]


def estimate_score(
    expr: pd.DataFrame,
    signatures: dict[str, list[str]],
    common_genes: list[str] | None = None,
    tumor_purity: bool = False,
    platform: str = "illumina",
) -> pd.DataFrame:
    """Stromal / Immune / ESTIMATE scores (and optionally purity).

    Parameters
    ----------
    expr
        genes x samples. Any within-sample monotone scale (the statistic ranks
        each column), but between-sample normalisation must already be done.
    common_genes
        if given, the matrix is filtered to this background first (recommended,
        and required to reproduce published values).
    tumor_purity
        emit the Affymetrix-calibrated purity transform. Only meaningful for
        Affymetrix-like data; on RNA-seq treat it as an ordinal proxy at best.
    """
    work = expr if common_genes is None else filter_common_genes(expr, common_genes)

    # ESTIMATE ranks with ties.method="average" and keeps the fractional ranks.
    scores = ssgsea(
        work,
        signatures,
        alpha=0.25,
        normalize="none",
        tie_method="average",
    )

    expected = {"StromalSignature", "ImmuneSignature"}
    if not expected.issubset(scores.index):
        raise ValueError(f"expected signatures {expected}, got {list(scores.index)}")

    out = pd.DataFrame(
        {
            "StromalScore": scores.loc["StromalSignature"],
            "ImmuneScore": scores.loc["ImmuneSignature"],
        }
    )
    out["ESTIMATEScore"] = out["StromalScore"] + out["ImmuneScore"]

    if tumor_purity:
        if platform != "affymetrix":
            print(
                "[estimate] WARNING: the purity transform was calibrated on "
                "Affymetrix U133A; on RNA-seq it is an uncalibrated proxy."
            )
        purity = np.cos(PURITY_INTERCEPT + PURITY_SLOPE * out["ESTIMATEScore"])
        out["TumorPurity"] = purity.where(purity >= 0)

    out.attrs["overlap"] = scores.attrs.get("overlap", {})
    out.attrs["n_genes_background"] = scores.attrs.get("n_genes_background")
    return out
