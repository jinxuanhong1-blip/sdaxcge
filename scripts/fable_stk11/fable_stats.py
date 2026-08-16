"""Small, dependency-light statistics helpers used across the slice."""
from __future__ import annotations

import numpy as np
from scipy import stats


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR. Returns array aligned with input order."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def cliffs_delta(a, b):
    """Cliff's delta effect size for a vs b (non-parametric)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return np.nan
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (a.size * b.size)


def cohens_d(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size < 2 or b.size < 2:
        return np.nan
    na, nb = a.size, b.size
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return np.nan
    return (a.mean() - b.mean()) / sp


def group_compare(vals_mut, vals_wt):
    """Compare a continuous variable between mutant vs wild-type groups."""
    a = np.asarray(vals_mut, float)
    b = np.asarray(vals_wt, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    res = {
        "n_mut": int(a.size), "n_wt": int(b.size),
        "median_mut": float(np.median(a)) if a.size else np.nan,
        "median_wt": float(np.median(b)) if b.size else np.nan,
        "mean_mut": float(a.mean()) if a.size else np.nan,
        "mean_wt": float(b.mean()) if b.size else np.nan,
    }
    if a.size >= 3 and b.size >= 3:
        u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
        t, pt = stats.ttest_ind(a, b, equal_var=False)
        res.update({
            "mwu_p": float(pu), "ttest_p": float(pt),
            "log2fc_median": res["median_mut"] - res["median_wt"],
            "cliffs_delta": float(cliffs_delta(a, b)),
            "cohens_d": float(cohens_d(a, b)),
        })
    else:
        res.update({"mwu_p": np.nan, "ttest_p": np.nan,
                    "log2fc_median": np.nan, "cliffs_delta": np.nan,
                    "cohens_d": np.nan})
    return res


def fisher_dcb(mut_flags, dcb_yes):
    """2x2 Fisher exact: rows=mutant/wt, cols=DCB yes/no."""
    mut = np.asarray(mut_flags, bool)
    dcb = np.asarray(dcb_yes, bool)
    a = int((mut & dcb).sum())     # mutant & DCB
    b = int((mut & ~dcb).sum())    # mutant & noDCB
    c = int((~mut & dcb).sum())    # wt & DCB
    d = int((~mut & ~dcb).sum())   # wt & noDCB
    table = [[a, b], [c, d]]
    orr, p = stats.fisher_exact(table, alternative="two-sided")
    return {
        "table": table,
        "dcb_rate_mut": a / (a + b) if (a + b) else np.nan,
        "dcb_rate_wt": c / (c + d) if (c + d) else np.nan,
        "odds_ratio": float(orr), "fisher_p": float(p),
        "n_mut": a + b, "n_wt": c + d,
    }


def mantel_haenszel_or(tables):
    """Mantel-Haenszel pooled odds ratio across a list of 2x2 [[a,b],[c,d]]."""
    num = den = 0.0
    for t in tables:
        (a, b), (c, d) = t
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
    return float(num / den) if den else np.nan
