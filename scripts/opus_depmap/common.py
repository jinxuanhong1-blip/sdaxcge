"""Shared paths, loaders and statistics helpers for the TACSTD2/CLDN4 lung slice.

Everything here operates on the open DepMap Public 24Q4 files fetched by
``00_download_depmap.py``. Nothing in this slice writes outside
``notes/opus_depmap``, ``scripts/opus_depmap`` and ``results/opus_depmap``
(plus the git-ignored local data cache).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "opus_depmap"
CACHE = DATA / "cache"
RESULTS = REPO / "results" / "opus_depmap"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
NOTES = REPO / "notes" / "opus_depmap"

RELEASE = "DepMap Public 24Q4"

GENES_OF_INTEREST = ["TACSTD2", "CLDN4"]

_SYMBOL_RE = re.compile(r"^(?P<symbol>.+?)\s+\((?P<entrez>\d+)\)$")


def ensure_dirs() -> None:
    for path in (CACHE, TABLES, FIGURES, NOTES):
        path.mkdir(parents=True, exist_ok=True)


def strip_entrez(columns) -> list[str]:
    """'TACSTD2 (4070)' -> 'TACSTD2'; leave already-plain symbols untouched."""
    out = []
    for col in columns:
        match = _SYMBOL_RE.match(str(col))
        out.append(match.group("symbol") if match else str(col))
    return out


def _load_matrix(csv_name: str, cache_name: str) -> pd.DataFrame:
    cache_path = CACHE / cache_name
    if cache_path.exists():
        return pd.read_parquet(cache_path)
    CACHE.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(DATA / csv_name, index_col=0)
    frame.index.name = "ModelID"
    frame.columns = strip_entrez(frame.columns)
    # A handful of symbols repeat across Entrez ids; keep the first occurrence
    # so downstream lookups by symbol stay unambiguous.
    frame = frame.loc[:, ~frame.columns.duplicated()]
    frame = frame.astype("float32")
    frame.to_parquet(cache_path)
    return frame


def load_gene_effect() -> pd.DataFrame:
    """Chronos gene effect, models x genes. 0 = neutral, -1 = median common essential."""
    return _load_matrix("CRISPRGeneEffect.csv", "gene_effect.parquet")


def load_gene_dependency() -> pd.DataFrame:
    """Posterior probability that a gene is a dependency in a given model."""
    return _load_matrix("CRISPRGeneDependency.csv", "gene_dependency.parquet")


def load_expression() -> pd.DataFrame:
    """log2(TPM+1) protein-coding expression, models x genes."""
    return _load_matrix(
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "expression.parquet"
    )


def load_models() -> pd.DataFrame:
    models = pd.read_csv(DATA / "Model.csv")
    return models.set_index("ModelID")


def lung_cohort(models: pd.DataFrame) -> pd.DataFrame:
    """Lung-lineage cancer models with a coarse NSCLC / SCLC / other label."""
    lung = models[models["OncotreeLineage"] == "Lung"].copy()
    lung = lung[lung["OncotreePrimaryDisease"] != "Non-Cancerous"]
    lung["LungGroup"] = np.where(
        lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer",
        "NSCLC",
        np.where(
            lung["OncotreeSubtype"].astype(str).str.contains("Small Cell Lung", case=False),
            "SCLC",
            "Other lung",
        ),
    )
    return lung


# --------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs propagate."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape, np.nan)
    ok = np.isfinite(pvals)
    p = pvals[ok]
    if p.size == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    n = p.size
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Rank-biserial / Cliff's delta effect size derived from Mann-Whitney U."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (a.size * b.size) - 1.0


def mannwhitney(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Return (U p-value, Cliff's delta, Hodges-Lehmann location shift a-b)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return np.nan, np.nan, np.nan
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    delta = cliffs_delta(a, b)
    # Hodges-Lehmann estimator, subsampled for large inputs to bound memory.
    rng = np.random.default_rng(0)
    aa = a if a.size <= 400 else rng.choice(a, 400, replace=False)
    bb = b if b.size <= 400 else rng.choice(b, 400, replace=False)
    shift = float(np.median(aa[:, None] - bb[None, :]))
    return float(p), float(delta), shift


def fisher_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Fisher z confidence interval for a correlation coefficient."""
    if not np.isfinite(r) or n < 4 or abs(r) >= 1:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    crit = stats.norm.ppf(1 - alpha / 2)
    return (float(np.tanh(z - crit * se)), float(np.tanh(z + crit * se)))


def corr_vector_vs_matrix(
    y: pd.Series, mat: pd.DataFrame, method: str = "pearson"
) -> pd.DataFrame:
    """Correlate one vector against every column of a matrix on shared rows.

    Uses complete observations per column, so columns with different missingness
    keep their own n. Returns r, n, two-sided p (t approximation) and BH q.
    """
    shared = y.index.intersection(mat.index)
    y = y.loc[shared].astype(float)
    x = mat.loc[shared].astype(float)
    if method == "spearman":
        y = y.rank()
        x = x.rank()
    mask = np.isfinite(x.to_numpy()) & np.isfinite(y.to_numpy())[:, None]
    xv = np.where(mask, x.to_numpy(), np.nan)
    yv = np.where(mask, y.to_numpy()[:, None], np.nan)
    n = mask.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        xm = np.nanmean(xv, axis=0)
        ym = np.nanmean(yv, axis=0)
        xc = xv - xm
        yc = yv - ym
        cov = np.nansum(xc * yc, axis=0)
        denom = np.sqrt(np.nansum(xc**2, axis=0) * np.nansum(yc**2, axis=0))
        r = np.where(denom > 0, cov / denom, np.nan)
        r = np.clip(r, -1.0, 1.0)
        t = r * np.sqrt(np.maximum(n - 2, 0) / np.maximum(1 - r**2, 1e-300))
        p = 2 * stats.t.sf(np.abs(t), df=np.maximum(n - 2, 1))
    p = np.where(n >= 4, p, np.nan)
    out = pd.DataFrame(
        {"gene": mat.columns, "r": r, "n": n, "p": p},
        index=range(len(mat.columns)),
    )
    out["q"] = bh_fdr(out["p"].to_numpy())
    return out


def partial_spearman(
    x: pd.Series, y: pd.Series, z: pd.Series
) -> tuple[float, float, int]:
    """Spearman partial correlation of x and y controlling for z.

    Computed as the Pearson correlation of the residuals of rank(x) and rank(y)
    after regressing each on rank(z); p-value uses df = n - 3.
    """
    frame = pd.concat([x, y, z], axis=1).dropna()
    n = len(frame)
    if n < 6:
        return np.nan, np.nan, n
    ranks = frame.rank()
    zr = ranks.iloc[:, 2].to_numpy()
    design = np.column_stack([np.ones(n), zr])

    def resid(v: np.ndarray) -> np.ndarray:
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    rx = resid(ranks.iloc[:, 0].to_numpy())
    ry = resid(ranks.iloc[:, 1].to_numpy())
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    if denom == 0:
        return np.nan, np.nan, n
    r = float((rx * ry).sum() / denom)
    r = float(np.clip(r, -1.0, 1.0))
    df = n - 3
    t = r * np.sqrt(df / max(1 - r**2, 1e-300))
    p = float(2 * stats.t.sf(abs(t), df=df))
    return r, p, n
