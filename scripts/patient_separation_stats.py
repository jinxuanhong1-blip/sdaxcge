#!/usr/bin/env python3
"""Patient-unit median contrasts, bootstrap intervals, and permutation p-values.

CosMx donors and concordant-4 donors/samples are the sampling units.
A confidence interval or a p-value is never computed with one cell as one draw.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np
from scipy import stats


def labels_within_groups(
    values: np.ndarray,
    groups: np.ndarray,
    mode: str,
    tie_ids: np.ndarray,
) -> np.ndarray:
    """Label rows inside each group.

    1 = high CLDN4, 0 = low CLDN4, -1 = unused middle.
    Ties break toward the high arm when tie_ids is larger.
    ``detected`` is count > 0 versus count == 0 (no rank cut).
    ``median`` / ``tertile`` / ``quartile`` keep equal-sized extremes;
    an odd middle row is unused.
    """
    values = np.asarray(values)
    groups = np.asarray(groups)
    tie_ids = np.asarray(tie_ids)
    labels = np.full(values.shape[0], -1, dtype=np.int8)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        v = values[idx]
        if mode == "detected":
            labels[idx] = (v > 0).astype(np.int8)
            continue
        order = np.lexsort((tie_ids[idx], v))
        ranked = idx[order]
        n = int(ranked.size)
        if mode == "median":
            k = n // 2
        elif mode == "tertile":
            k = n // 3
        elif mode == "quartile":
            k = n // 4
        else:
            raise ValueError(f"unknown label mode {mode}")
        if k < 1:
            continue
        labels[ranked[:k]] = 0
        labels[ranked[-k:]] = 1
    return labels


def permute_within_groups(values: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Permute values inside each group. Group sizes stay fixed."""
    out = np.array(values, copy=True)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        out[idx] = rng.permutation(out[idx])
    return out


def bootstrap_median_ci(
    values: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 4999,
) -> tuple[float, float, float]:
    """Point median and 95% percentile bootstrap interval.

    With five patient values the bootstrap median is one of those five numbers,
    and P(bootstrap median equals the sample max) is about 0.059, which is
    above 0.025. The 95% interval is then exactly the sample range. That
    identity is used directly so a Monte Carlo draw cannot jitter the bound.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.any(~np.isfinite(values)):
        return math.nan, math.nan, math.nan
    med = float(np.median(values))
    if values.size == 5:
        return med, float(np.min(values)), float(np.max(values))
    draws = rng.integers(0, values.size, size=(n_boot, values.size))
    stats_b = np.median(values[draws], axis=1)
    lo, hi = np.quantile(stats_b, [0.025, 0.975], method="linear")
    return med, float(lo), float(hi)


def percentile_interval(replicates: np.ndarray) -> tuple[float, float, int]:
    """95% percentile interval. Non-finite replicates are dropped."""
    reps = np.asarray(replicates, dtype=float)
    finite = reps[np.isfinite(reps)]
    if finite.size < 2:
        return math.nan, math.nan, int(finite.size)
    lo, hi = np.quantile(finite, [0.025, 0.975], method="linear")
    return float(lo), float(hi), int(finite.size)


def exclusion_magnitude(median: float, hi: float) -> float:
    """|median| when the interval lies entirely below 0, otherwise 0.

    Positive CLDN4-minus-low contrasts are not exclusion and do not compete.
    """
    if not np.isfinite(median) or not np.isfinite(hi):
        return 0.0
    if median < 0.0 and hi < 0.0:
        return float(-median)
    return 0.0


def pick_winner(rows: Sequence[dict]) -> dict | None:
    """Largest |median| among exclusion intervals. Narrower interval breaks ties."""
    ok = [r for r in rows if exclusion_magnitude(r["median"], r["hi"]) > 0.0 and r.get("eligible", True)]
    if not ok:
        return None

    def key(r: dict) -> tuple[float, float, float]:
        width = float(r["hi"] - r["lo"])
        return (abs(float(r["median"])), -width, -float(r["hi"]))

    ok.sort(key=key)
    return ok[-1]


def perm_p_ge(null_stats: np.ndarray, observed: float) -> float:
    """One-sided permutation p for a larger-is-more-extreme statistic. Includes the observed draw."""
    null_stats = np.asarray(null_stats, dtype=float)
    return float((1 + np.sum(null_stats >= observed)) / (null_stats.size + 1))


def corrected_log2_ratio(high: np.ndarray, low: np.ndarray, eps: float) -> np.ndarray:
    """log2((high+eps)/(low+eps)). eps keeps a zero high-arm finite and shrinks tiny denominators."""
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    return np.log2((high + eps) / (low + eps))


def rollup_section_means(
    section_high: np.ndarray,
    section_low: np.ndarray,
    sample_patient: np.ndarray,
    n_patients: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Unweighted mean of section means inside each patient. Shape (n_patients, K)."""
    section_high = np.asarray(section_high, dtype=float)
    section_low = np.asarray(section_low, dtype=float)
    k = section_high.shape[1]
    high = np.full((n_patients, k), np.nan)
    low = np.full((n_patients, k), np.nan)
    for p in range(n_patients):
        m = np.asarray(sample_patient) == p
        if not np.any(m):
            continue
        high[p] = _nanmean_stack(section_high[m])
        low[p] = _nanmean_stack(section_low[m])
    return high, low


def _nanmean_stack(block: np.ndarray) -> np.ndarray:
    cnt = np.sum(np.isfinite(block), axis=0)
    total = np.nansum(block, axis=0)
    out = np.full(block.shape[1], np.nan, dtype=float)
    ok = cnt > 0
    out[ok] = total[ok] / cnt[ok]
    return out


def pairwise_median_diff(y: np.ndarray, lab: np.ndarray, cohort: np.ndarray) -> float:
    """Median of within-cohort pairwise differences, high minus low.

    Every cohort must contribute at least two patients on each arm.
    Pairs are not formed across cohorts.
    """
    chunks: list[np.ndarray] = []
    y = np.asarray(y, dtype=float)
    lab = np.asarray(lab)
    cohort = np.asarray(cohort)
    for c in np.unique(cohort):
        h = y[(cohort == c) & (lab == 1)]
        low = y[(cohort == c) & (lab == 0)]
        if h.size < 2 or low.size < 2:
            return math.nan
        chunks.append((h[:, None] - low[None, :]).ravel())
    return float(np.median(np.concatenate(chunks)))


def pooled_median_diff(y: np.ndarray, lab: np.ndarray, cohort: np.ndarray | None = None) -> float:
    """Difference of medians after within-cohort labels are pooled. ``cohort`` is unused."""
    del cohort
    y = np.asarray(y, dtype=float)
    lab = np.asarray(lab)
    h = y[lab == 1]
    low = y[lab == 0]
    if h.size < 2 or low.size < 2:
        return math.nan
    return float(np.median(h) - np.median(low))


def median_of_cohort_diffs(y: np.ndarray, lab: np.ndarray, cohort: np.ndarray) -> float:
    """Median across cohorts of (median high − median low). Each cohort is one number."""
    y = np.asarray(y, dtype=float)
    lab = np.asarray(lab)
    cohort = np.asarray(cohort)
    diffs = []
    for c in np.unique(cohort):
        h = y[(cohort == c) & (lab == 1)]
        low = y[(cohort == c) & (lab == 0)]
        if h.size < 2 or low.size < 2:
            return math.nan
        diffs.append(float(np.median(h) - np.median(low)))
    return float(np.median(diffs))


ESTIMANDS: dict[str, Callable[[np.ndarray, np.ndarray, np.ndarray], float]] = {
    "hl_within_cohort": pairwise_median_diff,
    "pooled_diff_of_medians": pooled_median_diff,
    "median_of_cohort_diffs": median_of_cohort_diffs,
}


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    rho, _ = stats.spearmanr(x, y)
    return float(rho)


def der_simonian_laird(rhos: Sequence[float], ns: Sequence[int]) -> dict:
    """Random-effects pool of Spearman rhos on the Fisher-z scale."""
    # A bootstrap draw can hit |rho|=1 when a small cohort resamples to ties.
    # Clip before arctanh so one draw does not make the pooled z infinite.
    z = np.array([np.arctanh(np.clip(r, -0.999999, 0.999999)) for r in rhos], dtype=float)
    n = np.array(list(ns), dtype=float)
    if np.any(n <= 3):
        raise ValueError("Spearman variance 1/(n-3) needs n > 3")
    v = 1.0 / (n - 3.0)
    w = 1.0 / v
    zbar = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - zbar) ** 2))
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    k = z.size
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (v + tau2)
    z_hat = float(np.sum(w_star * z) / np.sum(w_star))
    se = float(np.sqrt(1.0 / np.sum(w_star)))
    rho = float(np.tanh(z_hat))
    zstat = z_hat / se if se > 0 else math.nan
    p = float(math.erfc(abs(zstat) / math.sqrt(2.0))) if np.isfinite(zstat) else math.nan
    ci_lo = float(np.tanh(z_hat - 1.96 * se)) if se > 0 else math.nan
    ci_hi = float(np.tanh(z_hat + 1.96 * se)) if se > 0 else math.nan
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {
        "rho": rho,
        "se_z": se,
        "z": z_hat,
        "p": p,
        "tau2": tau2,
        "i2": i2,
        "q": q,
        "k": int(k),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
    }


def cohort_spearman_dl(score: np.ndarray, outcome: np.ndarray, cohort: np.ndarray) -> dict:
    rhos = []
    ns = []
    names = []
    for c in np.unique(cohort):
        m = cohort == c
        rhos.append(spearman_rho(score[m], outcome[m]))
        ns.append(int(m.sum()))
        names.append(c)
    out = der_simonian_laird(rhos, ns)
    out["cohort_rho"] = rhos
    out["cohort_n"] = ns
    out["cohort"] = names
    return out


def resample_within_cohort(
    cohort: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stratified patient bootstrap indices. Each cohort keeps its sample size."""
    parts = []
    for c in np.unique(cohort):
        idx = np.flatnonzero(cohort == c)
        parts.append(idx[rng.integers(0, idx.size, size=idx.size)])
    return np.concatenate(parts)
