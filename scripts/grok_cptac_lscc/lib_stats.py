"""Shared statistics for the CPTAC LSCC TACSTD2/CLDN4 slice."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def strip_ensembl_version(ids: Iterable[str]) -> pd.Index:
    return pd.Index([str(x).split(".")[0] for x in ids])


def match_ensembl_row(index: Iterable[str], ensembl_id: str) -> str | None:
    """Return the first versioned index label whose prefix matches ensembl_id."""
    prefix = ensembl_id.split(".")[0]
    for raw in index:
        if str(raw).split(".")[0] == prefix:
            return str(raw)
    return None


def spearman_pair(x: pd.Series, y: pd.Series) -> dict:
    """Pairwise-complete Spearman. n < 4 -> NaN (parent convention)."""
    d = pd.concat([x, y], axis=1).dropna()
    n = int(d.shape[0])
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    a = d.iloc[:, 0].to_numpy()
    b = d.iloc[:, 1].to_numpy()
    if np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(a, b)
    return {"n": n, "rho": float(rho), "p": float(p)}


def mannwhitney_two_sided(a: pd.Series, b: pd.Series) -> dict:
    """Two-sided MWU + rank-biserial r = 1 - 2U/(n_a*n_b). n < 2 per group -> NaN."""
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    n_a, n_b = int(a.size), int(b.size)
    if n_a < 2 or n_b < 2:
        return {"n_a": n_a, "n_b": n_b, "U": np.nan, "p": np.nan, "rank_biserial": np.nan}
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 1.0 - (2.0 * float(res.statistic)) / (n_a * n_b)
    return {
        "n_a": n_a,
        "n_b": n_b,
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "rank_biserial": float(r),
    }


def wilcoxon_signed_rank(a: pd.Series, b: pd.Series) -> dict:
    """Paired Wilcoxon on aligned case IDs. n < 6 -> NaN."""
    d = pd.concat([a, b], axis=1).dropna()
    n = int(d.shape[0])
    if n < 6:
        return {"n": n, "W": np.nan, "p": np.nan, "median_delta": np.nan}
    delta = d.iloc[:, 0] - d.iloc[:, 1]
    if np.allclose(delta.values, 0):
        return {"n": n, "W": np.nan, "p": np.nan, "median_delta": 0.0}
    res = stats.wilcoxon(d.iloc[:, 0], d.iloc[:, 1], alternative="two-sided", zero_method="wilcox")
    return {
        "n": n,
        "W": float(res.statistic),
        "p": float(res.pvalue),
        "median_delta": float(delta.median()),
    }


def bh_fdr(pvalues: Iterable[float]) -> list[float]:
    """Benjamini-Hochberg FDR. NaNs preserved; ranking uses finite p only."""
    p = np.asarray(list(pvalues), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty(m, dtype=float)
    restored[order] = q
    out[ok] = restored
    return out.tolist()


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0)
    sd = sd.replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0)


def signature_score(expr_genes_x_samples: pd.DataFrame, gene_ids: list[str]) -> tuple[pd.Series, str]:
    """Mean of per-gene z-scores. Returns (score, coverage used/total)."""
    present = [g for g in gene_ids if g in expr_genes_x_samples.index]
    if not present:
        return pd.Series(np.nan, index=expr_genes_x_samples.columns), f"0/{len(gene_ids)}"
    z = zscore_rows(expr_genes_x_samples.loc[present])
    score = z.mean(axis=0, skipna=True)
    return score, f"{len(present)}/{len(gene_ids)}"


def median_split(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    med = x.median()
    out = pd.Series(np.nan, index=x.index, dtype=object)
    out[x >= med] = "high"
    out[x < med] = "low"
    out[x.isna()] = np.nan
    return out


def safe_log10_p(p: float) -> float:
    if p is None or not np.isfinite(p) or p <= 0:
        return np.nan
    return -math.log10(p)
