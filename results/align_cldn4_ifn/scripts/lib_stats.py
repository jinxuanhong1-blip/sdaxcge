#!/usr/bin/env python3
"""Shared statistics for the CLDN4 / TACSTD2 IFN-direction replication."""
import itertools

import numpy as np
import pandas as pd
from scipy import stats


# ------------------------------------------------------------ normalisation
def size_factors(counts):
    """DESeq2 median-of-ratios size factors (counts: genes x samples)."""
    x = counts.to_numpy(float)
    keep = (x > 0).all(axis=1)
    if keep.sum() < 50:  # too sparse; fall back to library size
        lib = x.sum(axis=0)
        return pd.Series(lib / np.exp(np.mean(np.log(lib))), index=counts.columns)
    logx = np.log(x[keep])
    ref = logx.mean(axis=1, keepdims=True)
    sf = np.exp(np.median(logx - ref, axis=0))
    return pd.Series(sf, index=counts.columns)


def cpm_log(counts, prior=1.0):
    sf = size_factors(counts)
    norm = counts.to_numpy(float) / sf.to_numpy()[None, :]
    norm = norm / norm.sum(axis=0, keepdims=True) * 1e6
    return pd.DataFrame(np.log2(norm + prior), index=counts.index, columns=counts.columns)


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


# ------------------------------------------------------- per-gene contrasts
def per_gene(logmat, grp_test, grp_ref):
    """log2FC (test - ref) and Welch t-test per gene on a log2 matrix."""
    a = logmat[grp_test].to_numpy(float)
    b = logmat[grp_ref].to_numpy(float)
    lfc = a.mean(axis=1) - b.mean(axis=1)
    if a.shape[1] > 1 and b.shape[1] > 1:
        t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    else:
        t = np.full(lfc.shape, np.nan)
        p = np.full(lfc.shape, np.nan)
    return pd.DataFrame({"log2FC": lfc, "t": t, "p": p, "fdr": bh(p)},
                        index=logmat.index)


# -------------------------------------------------- competitive set testing
def competitive(lfc, members):
    """Are the set's log2FC shifted up relative to all other genes?

    Mann-Whitney U on log2FC (set vs background). This is a *competitive*
    test: it asks about the set relative to the rest of the transcriptome and
    is therefore insensitive to global shifts, but inter-gene correlation
    within IFN modules makes its nominal p optimistic. Reported alongside a
    plain effect size (median log2FC difference) and the AUC.
    """
    lfc = lfc.dropna()
    inset = lfc.index.isin(members)
    n_in = int(inset.sum())
    if n_in < 3:
        return dict(n_set=n_in, n_bg=int((~inset).sum()), auc=np.nan,
                    p_up=np.nan, p_two=np.nan, med_set=np.nan, med_bg=np.nan,
                    mean_set=np.nan, mean_bg=np.nan, delta_med=np.nan)
    a, b = lfc[inset].to_numpy(), lfc[~inset].to_numpy()
    u_up, p_up = stats.mannwhitneyu(a, b, alternative="greater")
    _, p_two = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_set=n_in, n_bg=b.size, auc=u_up / (a.size * b.size),
                p_up=p_up, p_two=p_two,
                med_set=float(np.median(a)), med_bg=float(np.median(b)),
                mean_set=float(a.mean()), mean_bg=float(b.mean()),
                delta_med=float(np.median(a) - np.median(b)))


# ------------------------------------------- sample-level score + permutation
def sample_scores(logmat, members):
    """Mean z-score (across samples) of set genes -> one score per sample."""
    sub = logmat.loc[logmat.index.isin(members)]
    sub = sub.loc[sub.std(axis=1) > 0]
    if sub.shape[0] < 3:
        return None
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
    return z.mean(axis=0)


def perm_test(scores, grp_test, grp_ref, max_exact=20000):
    """Exact (or sampled) label-permutation test on the group score difference."""
    a = scores[grp_test].to_numpy(float)
    b = scores[grp_ref].to_numpy(float)
    obs = a.mean() - b.mean()
    allv = np.concatenate([a, b])
    n = allv.size
    k = a.size
    combos = list(itertools.combinations(range(n), k))
    if len(combos) > max_exact:
        rng = np.random.default_rng(0)
        combos = [tuple(rng.choice(n, k, replace=False)) for _ in range(max_exact)]
        exact = False
    else:
        exact = True
    diffs = np.array([allv[list(c)].mean() - allv[[i for i in range(n) if i not in c]].mean()
                      for c in combos])
    p_up = (np.sum(diffs >= obs - 1e-12)) / diffs.size
    p_two = (np.sum(np.abs(diffs) >= abs(obs) - 1e-12)) / diffs.size
    return dict(delta_score=float(obs), perm_p_up=float(p_up),
                perm_p_two=float(p_two), n_perm=int(diffs.size), exact=exact)
