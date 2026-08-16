"""Batch correction for bulk immune scores -- what to correct, and what not to.

Immune deconvolution scores are *derived* features. Correcting the expression
matrix with ComBat and then scoring is not the same as scoring and then
adjusting the scores, and both are different from including batch as a
covariate in the TACSTD2/CLDN4 correlation. The playbook's rule, implemented
here:

1. **Never ComBat-correct a CIBERSORT mixture.** Relative-mode SVR assumes a
   linear mixture on the original scale. ComBat's location/scale shift breaks
   that assumption and the fractions become uninterpretable.
2. **Rank-based scores (ssGSEA, ESTIMATE, xCell) are cohort-relative.**
   Running them on a ComBat-merged matrix is legitimate *if and only if* the
   samples will be analysed together afterwards. Do not ComBat two ICI
   cohorts, score, and then compare the score distributions across cohorts --
   the ranks have already been mixed.
3. **Prefer residualising the scores.** For a TACSTD2 ~ ImmuneScore question,
   fit ``ImmuneScore ~ batch + purity`` and correlate the residual with
   TACSTD2. That answers "is the association independent of batch" without
   rewriting the expression matrix.
4. **RNA-seq and microarray are not a ComBat problem, they are a method
   problem.** LM22 was built on Affymetrix; xCell ships separate spillover
   matrices; ESTIMATE purity is Affymetrix-only. Run the platform-appropriate
   configuration, do not force the matrices onto a common scale and pretend
   the scores are commensurate.

The ComBat implementation below is the parametric empirical-Bayes procedure
of Johnson, Li & Rabinovic (Biostatistics 2007), sufficient for the
location/scale adjustment used in the demo. It is not sva::ComBat
feature-for-feature (no nonparametric option, no reference-batch mode).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["combat_parametric", "residualise", "batch_association"]


def _design(batch: pd.Series, covariates: pd.DataFrame | None) -> tuple[np.ndarray, int]:
    batch = batch.astype("category")
    dummies = pd.get_dummies(batch, drop_first=False).to_numpy(dtype=float)
    if covariates is None or covariates.shape[1] == 0:
        return dummies, dummies.shape[1]
    cov = covariates.reindex(batch.index).to_numpy(dtype=float)
    if np.isnan(cov).any():
        raise ValueError("covariates contain NaN; impute or drop those samples first")
    return np.column_stack([dummies, cov]), dummies.shape[1]


def combat_parametric(
    expr: pd.DataFrame,
    batch: pd.Series,
    covariates: pd.DataFrame | None = None,
    mean_only: bool = False,
    prior_plots: bool = False,
) -> pd.DataFrame:
    """Parametric ComBat. ``expr`` is genes x samples; ``batch`` is a sample-aligned series.

    Samples in a singleton batch are returned unchanged (there is nothing to
    estimate). Genes with zero variance inside a batch are left unadjusted
    for that batch.
    """
    batch = batch.loc[expr.columns]
    counts = batch.value_counts()
    if (counts < 2).all():
        raise ValueError("every batch has n<2; ComBat cannot estimate batch effects")

    Y = expr.to_numpy(dtype=float)
    n_genes, n_samples = Y.shape
    design, n_batch = _design(batch, covariates)
    # OLS fit of gene ~ batch + covariates
    beta, *_ = np.linalg.lstsq(design, Y.T, rcond=None)
    beta = beta.T  # genes x design
    grand = Y.mean(axis=1, keepdims=True)
    # Standardise as in Johnson 2007
    batch_levels = list(pd.Categorical(batch).categories)
    batch_idx = {b: np.where(batch.to_numpy() == b)[0] for b in batch_levels}

    # Residual variance after the design fit
    fitted = design @ beta.T
    resid = Y - fitted.T
    var = resid.var(axis=1, ddof=1)
    var[var == 0] = 1.0
    Z = (Y - grand) / np.sqrt(var)[:, None]

    # Per-batch location / scale on the standardised matrix
    gamma_hat = np.zeros((n_genes, len(batch_levels)))
    delta_hat = np.ones((n_genes, len(batch_levels)))
    for j, b in enumerate(batch_levels):
        idx = batch_idx[b]
        gamma_hat[:, j] = Z[:, idx].mean(axis=1)
        if not mean_only and idx.size > 1:
            d = Z[:, idx].var(axis=1, ddof=1)
            d[d == 0] = 1.0
            delta_hat[:, j] = d

    # Empirical Bayes shrinkage of gamma (location)
    gamma_bar = gamma_hat.mean(axis=0)
    t2 = gamma_hat.var(axis=0, ddof=1)
    t2[t2 == 0] = 1e-6
    # Inverse-gamma prior on delta
    # method of moments as in the original ComBat
    def _ig_moments(x):
        m = x.mean()
        v = x.var(ddof=1)
        if v <= 0 or m <= 0:
            return 1.0, 1.0
        a = (m ** 2 / v) + 2.0
        b = m * (a - 1.0)
        return a, b

    gamma_star = np.zeros_like(gamma_hat)
    delta_star = np.ones_like(delta_hat)
    for j, b in enumerate(batch_levels):
        idx = batch_idx[b]
        n_j = idx.size
        if n_j < 2:
            gamma_star[:, j] = 0.0
            delta_star[:, j] = 1.0
            continue
        denom = n_j * t2[j] + delta_hat[:, j]
        gamma_star[:, j] = (n_j * t2[j] * gamma_hat[:, j] + delta_hat[:, j] * gamma_bar[j]) / denom
        if mean_only:
            delta_star[:, j] = 1.0
        else:
            a, bpar = _ig_moments(delta_hat[:, j])
            # posterior scale: (prior + SSE) / (a - 1 + n/2)  -- Johnson eq.
            sse = ((Z[:, idx] - gamma_star[:, j][:, None]) ** 2).sum(axis=1)
            delta_star[:, j] = (bpar + 0.5 * sse) / (a - 1.0 + 0.5 * n_j)
            delta_star[delta_star[:, j] <= 0, j] = 1.0

    adjusted = Y.copy()
    for j, b in enumerate(batch_levels):
        idx = batch_idx[b]
        if idx.size < 2:
            continue
        Z_adj = (Z[:, idx] - gamma_star[:, j][:, None]) / np.sqrt(delta_star[:, j])[:, None]
        adjusted[:, idx] = Z_adj * np.sqrt(var)[:, None] + grand

    out = pd.DataFrame(adjusted, index=expr.index, columns=expr.columns)
    if prior_plots:
        out.attrs["gamma_hat"] = gamma_hat
        out.attrs["gamma_star"] = gamma_star
        out.attrs["delta_hat"] = delta_hat
        out.attrs["delta_star"] = delta_star
    return out


def residualise(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    """OLS residuals of each score column on the covariate matrix.

    Use this to ask "is TACSTD2 still associated with ImmuneScore after
    removing batch and purity" without rewriting the expression matrix.
    Categorical covariates should already be dummy-coded, or passed as
    object/category columns (they are dummy-coded here, dropping the first
    level).
    """
    cov = covariates.reindex(scores.index)
    frames = []
    for col in cov.columns:
        s = cov[col]
        if pd.api.types.is_numeric_dtype(s):
            frames.append(s.astype(float).rename(col))
        else:
            frames.append(pd.get_dummies(s, prefix=col, drop_first=True).astype(float))
    X = pd.concat(frames, axis=1) if frames else pd.DataFrame(index=scores.index)
    X = X.dropna(axis=0, how="any")
    X.insert(0, "_intercept", 1.0)
    common = scores.index.intersection(X.index)
    Y = scores.loc[common].to_numpy(dtype=float)
    A = X.loc[common].to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(A, Y, rcond=None)
    resid = Y - A @ beta
    return pd.DataFrame(resid, index=common, columns=scores.columns)


def batch_association(
    scores: pd.DataFrame,
    batch: pd.Series,
) -> pd.DataFrame:
    """Kruskal-Wallis (or Mann-Whitney if 2 levels) of each score vs batch.

    A score that is more associated with batch than with TACSTD2 is not a
    usable immune readout in that cohort -- report the test, do not "fix"
    it by ComBat and then claim a TACSTD2 association.
    """
    from scipy.stats import kruskal, mannwhitneyu

    batch = batch.loc[scores.index]
    rows = []
    levels = list(pd.unique(batch.dropna()))
    for col in scores.columns:
        groups = [scores.loc[batch == lv, col].dropna().to_numpy() for lv in levels]
        groups = [g for g in groups if g.size]
        if len(groups) < 2:
            stat, p = np.nan, np.nan
            test = "skipped"
        elif len(groups) == 2:
            stat, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
            test = "mannwhitneyu"
        else:
            stat, p = kruskal(*groups)
            test = "kruskal"
        rows.append({"score": col, "test": test, "n_levels": len(groups), "statistic": stat, "p": p})
    return pd.DataFrame(rows)
