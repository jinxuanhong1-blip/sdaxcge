"""Shared helpers for the fable_neoadj TACSTD2/CLDN4-vs-MPR analysis.

All datasets are open neoadjuvant lung (NSCLC) PD-1/PD-L1 +/- chemo cohorts
with pathologic-response (MPR / pCR) annotation. Helpers here keep statistics
and plotting consistent across the per-dataset scripts.
"""
from __future__ import annotations

import os
import numpy as np
from scipy import stats

GENES = ["TACSTD2", "CLDN4"]

# Repo-relative output roots (all outputs must live under results/fable_neoadj)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, "results", "fable_neoadj")
DATA = os.path.join(RES, "data")
FIG = os.path.join(RES, "figures")
TAB = os.path.join(RES, "tables")
for _d in (FIG, TAB):
    os.makedirs(_d, exist_ok=True)


def cliffs_delta(a, b):
    """Cliff's delta effect size for group a vs b (positive => a tends higher)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def mwu_stats(pos, neg):
    """Mann-Whitney U comparing responders (pos) vs non-responders (neg).

    Returns dict with U, p (two-sided), rank-biserial AUC = P(pos>neg), and
    Cliff's delta. AUC>0.5 => marker higher in responders.
    """
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos = pos[~np.isnan(pos)]
    neg = neg[~np.isnan(neg)]
    out = {
        "n_pos": int(len(pos)),
        "n_neg": int(len(neg)),
        "median_pos": float(np.median(pos)) if len(pos) else np.nan,
        "median_neg": float(np.median(neg)) if len(neg) else np.nan,
    }
    if len(pos) < 2 or len(neg) < 2:
        out.update({"U": np.nan, "p": np.nan, "auc": np.nan, "cliffs_delta": np.nan})
        return out
    U, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = U / (len(pos) * len(neg))  # ROC-AUC of marker ranking responders above non-responders
    out.update({"U": float(U), "p": float(p), "auc": float(auc),
                "cliffs_delta": float(cliffs_delta(pos, neg))})
    return out


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return {"rho": np.nan, "p": np.nan, "n": int(m.sum())}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"rho": float(rho), "p": float(p), "n": int(m.sum())}


def stouffer(pvals, signs, weights=None):
    """Combine two-sided p-values with directions into one signed Z / p.

    signs: +1 if marker higher in responders in that study, -1 otherwise.
    """
    pvals = np.asarray(pvals, float)
    signs = np.asarray(signs, float)
    m = ~np.isnan(pvals) & ~np.isnan(signs)
    pvals, signs = pvals[m], signs[m]
    if len(pvals) == 0:
        return {"z": np.nan, "p": np.nan, "k": 0}
    if weights is None:
        weights = np.ones(len(pvals))
    else:
        weights = np.asarray(weights, float)[m]
    # two-sided p -> one-sided z in the observed direction
    z_one = stats.norm.isf(pvals / 2.0) * signs
    z = np.sum(weights * z_one) / np.sqrt(np.sum(weights ** 2))
    p = 2 * stats.norm.sf(abs(z))
    return {"z": float(z), "p": float(p), "k": int(len(pvals))}


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR; NaNs preserved."""
    p = np.asarray(pvals, float)
    out = np.full_like(p, np.nan, dtype=float)
    idx = np.where(~np.isnan(p))[0]
    if len(idx) == 0:
        return out
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(pv)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    res = np.empty(n)
    res[order] = q
    out[idx] = res
    return out
