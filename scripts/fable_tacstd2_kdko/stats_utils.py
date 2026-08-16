"""Statistics helpers: per-gene two-group tests and competitive gene-set tests.

Kept deliberately transparent so the numbers are auditable:
  * per-gene: Welch's t-test + Mann-Whitney U on log-scale expression, log2FC as
    difference of group means in log2 space, Cohen's d effect size.
  * multiple testing: Benjamini-Hochberg FDR across all tested genes.
  * gene-set: a simplified CAMERA-style *competitive* test (Mann-Whitney of the
    per-gene ranking statistic for set members vs all other genes) plus a
    *self-contained* one-sample Wilcoxon that the set's log2FC distribution is
    shifted from zero. Directionality is reported explicitly.
"""
import numpy as np
from scipy import stats


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def cohens_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sp2 = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if sp2 <= 0:
        return np.nan
    return (a.mean() - b.mean()) / np.sqrt(sp2)


def per_gene_two_group(logmat, case_idx, ctrl_idx):
    """logmat: (genes x samples) log-scale array. Returns dict of arrays.

    log2FC = mean(case) - mean(ctrl)  (so >0 = up upon perturbation/KO/KD).
    """
    case = logmat[:, case_idx]
    ctrl = logmat[:, ctrl_idx]
    mean_case = case.mean(axis=1)
    mean_ctrl = ctrl.mean(axis=1)
    log2fc = mean_case - mean_ctrl

    n = logmat.shape[0]
    tp = np.full(n, np.nan)
    up = np.full(n, np.nan)
    d = np.full(n, np.nan)
    for i in range(n):
        a = case[i]; b = ctrl[i]
        # Welch t-test
        if np.std(a) + np.std(b) == 0:
            tp[i] = 1.0
        else:
            tp[i] = stats.ttest_ind(a, b, equal_var=False).pvalue
        # Mann-Whitney (rank based)
        try:
            up[i] = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            up[i] = 1.0
        d[i] = cohens_d(a, b)
    tp = np.nan_to_num(tp, nan=1.0)
    up = np.nan_to_num(up, nan=1.0)
    return {
        "log2FC": log2fc,
        "mean_case": mean_case,
        "mean_ctrl": mean_ctrl,
        "t_pvalue": tp,
        "t_fdr": bh_fdr(tp),
        "mwu_pvalue": up,
        "cohens_d": d,
    }


def competitive_set_test(stat_all, in_set_mask):
    """CAMERA-like competitive test: are set members' ranking stats shifted vs
    the rest of the genome? Uses Mann-Whitney U. Returns two-sided p, a directional
    call, and a rank-biserial effect size (AUC-like, in [-1,1])."""
    x = np.asarray(stat_all, float)
    m = np.asarray(in_set_mask, bool)
    inset = x[m]
    rest = x[~m]
    inset = inset[np.isfinite(inset)]
    rest = rest[np.isfinite(rest)]
    if len(inset) < 3 or len(rest) < 3:
        return {"n_in_set": len(inset), "p_two_sided": np.nan,
                "direction": "NA", "rank_biserial": np.nan,
                "median_set": np.nan, "median_rest": np.nan}
    U, p = stats.mannwhitneyu(inset, rest, alternative="two-sided")
    # rank-biserial correlation = 2*U/(n1*n2) - 1  (prob set > rest, rescaled)
    rb = 2 * U / (len(inset) * len(rest)) - 1
    direction = "up" if np.median(inset) > np.median(rest) else "down"
    return {
        "n_in_set": int(len(inset)),
        "p_two_sided": float(p),
        "direction": direction,
        "rank_biserial": float(rb),
        "median_set": float(np.median(inset)),
        "median_rest": float(np.median(rest)),
    }


def selfcontained_set_test(log2fc_set):
    """Self-contained: is the set's log2FC distribution shifted from 0?
    Wilcoxon signed-rank (two-sided) + mean log2FC + sign counts."""
    v = np.asarray(log2fc_set, float)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return {"n": len(v), "wilcoxon_p": np.nan, "mean_log2FC": np.nan,
                "n_up": int(np.sum(v > 0)), "n_down": int(np.sum(v < 0))}
    try:
        p = stats.wilcoxon(v).pvalue
    except ValueError:
        p = np.nan
    return {
        "n": int(len(v)),
        "wilcoxon_p": float(p),
        "mean_log2FC": float(np.mean(v)),
        "median_log2FC": float(np.median(v)),
        "n_up": int(np.sum(v > 0)),
        "n_down": int(np.sum(v < 0)),
    }
