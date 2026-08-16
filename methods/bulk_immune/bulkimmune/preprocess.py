"""Expression normalisation, on the scale each downstream method expects.

The single most common error in this analysis is feeding a method a matrix on
the wrong scale. The rules implemented here:

===================  ==========================================  ==============
method               required input                               linear/log
===================  ==========================================  ==============
CIBERSORT(x)         non-log, per-sample comparable (TPM/RMA)     linear
quanTIseq-style      TPM                                          linear
ESTIMATE             any within-sample monotone scale (rank based) either
MCP-counter          log2 expression (scores are means of logs)   log2
xCell                any within-sample monotone scale (rank based) either
ssGSEA               any within-sample monotone scale (rank based) either
TIDE                 log2(TPM+1), row-centred across samples      log2
===================  ==========================================  ==============

Rank-based methods are invariant to *within-sample* monotone transforms but not
to *between-sample* normalisation, so library-size normalisation still matters
for all of them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "cpm",
    "tpm_from_counts",
    "log2p1",
    "filter_low_expression",
    "quantile_normalise",
    "row_center",
    "looks_log_scale",
    "detect_scale",
]


def cpm(counts: pd.DataFrame, log: bool = False, prior_count: float = 1.0) -> pd.DataFrame:
    """Counts per million. ``log=True`` gives log2(CPM + prior_count)."""
    lib = counts.sum(axis=0)
    if (lib <= 0).any():
        raise ValueError(f"samples with zero library size: {list(lib[lib <= 0].index)}")
    out = counts.div(lib, axis=1) * 1e6
    return np.log2(out + prior_count) if log else out


def tpm_from_counts(
    counts: pd.DataFrame,
    gene_lengths: pd.Series,
    log: bool = False,
    prior_count: float = 1.0,
) -> pd.DataFrame:
    """TPM from counts and effective gene lengths (kb-independent scaling)."""
    common = counts.index.intersection(gene_lengths.dropna().index)
    if len(common) == 0:
        raise ValueError("no overlap between counts and gene_lengths")
    sub = counts.loc[common]
    rate = sub.div(gene_lengths.loc[common].astype(float), axis=0)
    out = rate.div(rate.sum(axis=0), axis=1) * 1e6
    return np.log2(out + prior_count) if log else out


def log2p1(expr: pd.DataFrame) -> pd.DataFrame:
    return np.log2(expr + 1.0)


def filter_low_expression(
    expr: pd.DataFrame,
    min_value: float = 1.0,
    min_samples: int | float = 0.2,
    keep: list[str] | None = None,
) -> pd.DataFrame:
    """Drop genes below ``min_value`` in all but a few samples.

    ``min_samples`` is a count when >= 1 and a fraction of samples when < 1.
    ``keep`` protects genes you must not lose (e.g. the target genes and the
    signature genes) -- filtering a signature gene out is indistinguishable from
    it being absent, and quietly shrinks the signature.
    """
    n = expr.shape[1]
    thresh = int(np.ceil(min_samples * n)) if min_samples < 1 else int(min_samples)
    mask = (expr >= min_value).sum(axis=1) >= thresh
    if keep:
        mask = mask | expr.index.isin(keep)
    return expr.loc[mask]


def quantile_normalise(expr: pd.DataFrame) -> pd.DataFrame:
    """Column-wise quantile normalisation (CIBERSORT's default for arrays)."""
    ranks = expr.rank(axis=0, method="average")
    mean_sorted = pd.Series(
        np.sort(expr.to_numpy(), axis=0).mean(axis=1), index=np.arange(1, expr.shape[0] + 1)
    )
    lo = np.floor(ranks.to_numpy()).astype(int)
    hi = np.ceil(ranks.to_numpy()).astype(int)
    frac = ranks.to_numpy() - lo
    values = mean_sorted.reindex(lo.ravel()).to_numpy().reshape(lo.shape) * (1 - frac) + (
        mean_sorted.reindex(hi.ravel()).to_numpy().reshape(hi.shape) * frac
    )
    return pd.DataFrame(values, index=expr.index, columns=expr.columns)


def row_center(expr: pd.DataFrame) -> pd.DataFrame:
    """Subtract each gene's mean across samples (TIDE's expected input)."""
    return expr.sub(expr.mean(axis=1), axis=0)


def looks_log_scale(expr: pd.DataFrame, threshold: float = 50.0) -> bool:
    """Heuristic: log-scale matrices have small maxima."""
    return float(np.nanmax(expr.to_numpy())) < threshold


def detect_scale(expr: pd.DataFrame) -> dict[str, float | bool | str]:
    """Describe the scale of a matrix so the caller can assert what it expects.

    Returns the column sums (1e6 => TPM/CPM), the maximum, whether values look
    log-transformed, and whether the values are integers (=> raw counts).
    """
    arr = expr.to_numpy(dtype=float)
    col_sums = np.nansum(arr, axis=0)
    finite = arr[np.isfinite(arr)]
    is_integer = bool(np.allclose(finite, np.round(finite)))
    med_sum = float(np.median(col_sums))
    if abs(med_sum - 1e6) / 1e6 < 0.02:
        guess = "TPM/CPM (columns sum to 1e6)"
    elif is_integer and med_sum > 1e5:
        guess = "counts"
    elif looks_log_scale(expr):
        guess = "log-transformed"
    else:
        guess = "linear, not sum-normalised (FPKM or similar)"
    return {
        "median_column_sum": med_sum,
        "min_column_sum": float(np.min(col_sums)),
        "max_column_sum": float(np.max(col_sums)),
        "max_value": float(np.nanmax(arr)),
        "min_value": float(np.nanmin(arr)),
        "is_integer": is_integer,
        "looks_log": looks_log_scale(expr),
        "guess": guess,
    }
