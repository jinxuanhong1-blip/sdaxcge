"""Single-sample gene set enrichment analysis (ssGSEA).

Implements the Barbie et al. (2009) statistic exactly as it is computed by
``GSVA::gsva(..., method="ssgsea")`` (Hanzelmann et al. 2013), which is also the
statistic used inside ESTIMATE, xCell and TIP:

    1. Rank the genes within each sample.
    2. Weight the ranks by ``|rank| ** alpha`` (alpha = 0.25 by default).
    3. Walk down the ranking and accumulate two step CDFs, one over the genes in
       the set (weighted) and one over the genes outside the set (unweighted).
    4. The enrichment score is the *sum* of the difference of the two CDFs
       (an integral, not the maximum deviation used by classical GSEA).

Two knobs matter for reproducing published numbers and are made explicit here
rather than hidden:

``tie_method``
    ``"average_int"`` reproduces GSVA (which computes average ranks and then
    truncates them to integers) and ``"average"`` reproduces ESTIMATE.
    The difference only appears when a sample has tied expression values --
    which, for RNA-seq with many zeros, is the rule rather than the exception.

``normalize``
    ``"none"`` returns the raw enrichment score. ``"global_range"`` reproduces
    the GSVA default (``ssgsea.norm=TRUE``), which divides every score by the
    range of the *whole* score matrix. That normalisation makes scores depend on
    which samples were run together, so it must not be used when scores from
    different batches or cohorts will be compared. See playbook section 8.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = ["rank_matrix", "ssgsea", "ssgsea_scores"]


def rank_matrix(expr: pd.DataFrame, tie_method: str = "average_int") -> np.ndarray:
    """Rank genes within each sample (column), largest expression = largest rank.

    Parameters
    ----------
    expr
        genes x samples expression matrix.
    tie_method
        ``"average_int"``  average ranks truncated toward zero (GSVA behaviour),
        ``"average"``      average ranks kept as floats (ESTIMATE behaviour),
        ``"min"``/``"max"``/``"first"`` pass straight through to scipy.
    """
    from scipy.stats import rankdata

    values = np.asarray(expr, dtype=float)
    if tie_method in ("average_int", "average"):
        ranks = np.apply_along_axis(lambda col: rankdata(col, method="average"), 0, values)
        if tie_method == "average_int":
            # R's as.integer() truncates toward zero, it does not round.
            ranks = np.trunc(ranks)
    else:
        ranks = np.apply_along_axis(lambda col: rankdata(col, method=tie_method), 0, values)
    return ranks


def _es_for_set(
    order: np.ndarray,
    weights_sorted: np.ndarray,
    member_sorted: np.ndarray,
    n_out: int,
) -> np.ndarray:
    """Enrichment score for one gene set across all samples.

    ``order`` is unused directly; the caller pre-sorts ``weights_sorted`` and
    ``member_sorted`` (both genes x samples, ordered by decreasing rank).
    """
    weighted_in = weights_sorted * member_sorted
    denom_in = weighted_in.sum(axis=0)
    # A gene set with zero total weight cannot produce a meaningful walk.
    denom_in = np.where(denom_in == 0, np.nan, denom_in)
    cdf_in = np.cumsum(weighted_in, axis=0) / denom_in
    cdf_out = np.cumsum(1.0 - member_sorted, axis=0) / n_out
    return (cdf_in - cdf_out).sum(axis=0)


def ssgsea(
    expr: pd.DataFrame,
    gene_sets: Mapping[str, Iterable[str]],
    alpha: float = 0.25,
    normalize: str = "none",
    tie_method: str = "average_int",
    min_size: int = 1,
    max_size: int | None = None,
    drop_constant: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """Compute ssGSEA scores.

    Parameters
    ----------
    expr
        genes x samples matrix. Should already be on a within-sample comparable
        scale (log2 TPM/CPM or normalised intensities). ssGSEA is rank based, so
        any strictly monotone transform of a sample leaves the score unchanged;
        it is *not* invariant to between-sample normalisation.
    gene_sets
        mapping of set name -> gene identifiers, matching ``expr`` row names.
    alpha
        rank exponent (GSVA's ``tau``). 0.25 is the published default.
    normalize
        ``"none"`` or ``"global_range"``.
    min_size, max_size
        sets whose *observed* overlap with ``expr`` falls outside these bounds
        are dropped (and reported in the returned frame's ``attrs``).

    Returns
    -------
    DataFrame
        gene sets x samples. ``df.attrs["overlap"]`` holds the number of genes of
        each set that were found, and ``df.attrs["dropped"]`` the discarded sets.
        Always report the overlap: an IFN-gamma score computed on 12/18 genes is
        not the published signature.
    """
    if normalize not in ("none", "global_range"):
        raise ValueError(f"unknown normalize={normalize!r}")

    work = expr
    if drop_constant:
        keep = work.std(axis=1, ddof=1) > 0
        if verbose and (~keep).any():
            print(f"[ssgsea] dropping {int((~keep).sum())} constant genes")
        work = work.loc[keep]

    if work.index.has_duplicates:
        raise ValueError(
            "expression matrix has duplicate gene identifiers; collapse them first "
            "(bulkimmune.genes.collapse_duplicates)"
        )

    n_genes, n_samples = work.shape
    ranks = rank_matrix(work, tie_method=tie_method)
    weights = np.abs(ranks) ** alpha

    # Decreasing rank order per sample; ties keep their matrix order, as in R's
    # order(..., decreasing=TRUE) which is stable.
    order = np.argsort(-ranks, axis=0, kind="stable")
    cols = np.arange(n_samples)
    weights_sorted = weights[order, cols]

    index = {g: i for i, g in enumerate(work.index)}
    scores: dict[str, np.ndarray] = {}
    overlap: dict[str, int] = {}
    dropped: dict[str, str] = {}

    for name, members in gene_sets.items():
        idx = np.fromiter(
            (index[g] for g in dict.fromkeys(members) if g in index), dtype=int
        )
        overlap[name] = int(idx.size)
        if idx.size < max(1, min_size):
            dropped[name] = f"overlap {idx.size} < min_size {min_size}"
            continue
        if max_size is not None and idx.size > max_size:
            dropped[name] = f"overlap {idx.size} > max_size {max_size}"
            continue
        n_out = n_genes - idx.size
        if n_out <= 0:
            dropped[name] = "set covers the whole expression matrix"
            continue
        member = np.zeros(n_genes, dtype=float)
        member[idx] = 1.0
        member_sorted = member[order]
        scores[name] = _es_for_set(order, weights_sorted, member_sorted, n_out)

    if not scores:
        raise ValueError("no gene set survived filtering; check identifier namespace")

    out = pd.DataFrame(scores, index=work.columns).T

    if normalize == "global_range":
        span = np.nanmax(out.to_numpy()) - np.nanmin(out.to_numpy())
        if span > 0:
            out = out / span

    out.attrs["overlap"] = overlap
    out.attrs["dropped"] = dropped
    out.attrs["alpha"] = alpha
    out.attrs["normalize"] = normalize
    out.attrs["tie_method"] = tie_method
    out.attrs["n_genes_background"] = n_genes
    return out


def ssgsea_scores(
    expr: pd.DataFrame,
    gene_sets: Mapping[str, Sequence[str]],
    **kwargs,
) -> pd.DataFrame:
    """Convenience wrapper returning samples x gene sets (tidy orientation)."""
    res = ssgsea(expr, gene_sets, **kwargs)
    out = res.T
    out.attrs.update(res.attrs)
    return out
