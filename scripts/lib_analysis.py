#!/usr/bin/env python3
"""Shared helpers for the claim A8/A9 TROP2 keratin / tight-junction analysis.

Design (identical logic applied to every dataset so results are comparable):

  * TROP2 = TACSTD2 expression.
  * Ranking metric for pre-ranked GSEA = Spearman correlation of each gene with
    TACSTD2 across samples/cells. This is the continuous, cut-off-free version of
    "up in TROP2-high" and is what we feed to GSEA prerank.
  * TROP2-high vs TROP2-low contrast = top vs bottom tertile of TACSTD2. We report
    log2 fold-change, Welch t-test, Mann-Whitney U and BH-FDR for every gene. This
    gives an explicit "high vs low" call for the CLDN1/CLDN4/CLDN7/F11R/PARD3 panel.

A gene is called "up in TROP2-high" (per dataset) when it is concordant
(rho > 0 AND log2FC_high_vs_low > 0) and significant (tertile Welch FDR < 0.05).
Everything is written out so failures are visible -- no cherry-picking.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

PANEL = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]
TROP2 = "TACSTD2"
GMT = "data/genesets/claim_A8A9_genesets.gmt"
FDR_THRESHOLD = 0.05


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pp = p[ok]
    n = pp.size
    if n == 0:
        return q
    order = np.argsort(pp)
    ranked = pp[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return q


def per_gene_stats(expr, trop2, min_frac_expressed=0.0, tertile_q=(1/3, 2/3)):
    """expr: genes x samples DataFrame (log-space). trop2: Series over samples.

    Returns a per-gene DataFrame with Spearman vs TACSTD2 and top/bottom-tertile
    TROP2-high vs TROP2-low differential statistics.
    """
    expr = expr.loc[:, trop2.index]
    trop2 = trop2.astype(float)
    X = expr.to_numpy(dtype=float)
    genes = expr.index.to_numpy()
    t = trop2.to_numpy(dtype=float)

    # optional expression filter (fraction of samples/cells with value > 0)
    if min_frac_expressed > 0:
        frac = (X > 0).mean(axis=1)
        keep = frac >= min_frac_expressed
        X, genes = X[keep], genes[keep]

    # ---- Spearman correlation of each gene with TACSTD2 ----
    tr = stats.rankdata(t)
    tr = (tr - tr.mean())
    tr_norm = tr / np.sqrt((tr ** 2).sum())
    Xr = np.apply_along_axis(stats.rankdata, 1, X)
    Xr = Xr - Xr.mean(axis=1, keepdims=True)
    denom = np.sqrt((Xr ** 2).sum(axis=1))
    denom[denom == 0] = np.nan
    rho = (Xr @ tr_norm) / denom
    n = t.size
    with np.errstate(invalid="ignore", divide="ignore"):
        tval = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    rho_p = 2 * stats.t.sf(np.abs(tval), df=n - 2)

    # ---- tertile high vs low ----
    lo_thr, hi_thr = np.quantile(t, tertile_q)
    low_mask = t <= lo_thr
    high_mask = t >= hi_thr
    Xl, Xh = X[:, low_mask], X[:, high_mask]
    mean_low = Xl.mean(axis=1)
    mean_high = Xh.mean(axis=1)
    log2fc = mean_high - mean_low  # expr already log2 space
    tt, tp = stats.ttest_ind(Xh, Xl, axis=1, equal_var=False)
    # Mann-Whitney per gene (vectorised via ranking is complex; loop is fine here)
    mwu_p = np.empty(X.shape[0])
    for i in range(X.shape[0]):
        try:
            mwu_p[i] = stats.mannwhitneyu(Xh[i], Xl[i], alternative="two-sided").pvalue
        except ValueError:
            mwu_p[i] = np.nan

    df = pd.DataFrame({
        "gene": genes,
        "mean_high": mean_high, "mean_low": mean_low,
        "log2FC_high_vs_low": log2fc,
        "spearman_rho": rho, "spearman_p": rho_p,
        "welch_t": tt, "welch_p": tp,
        "mwu_p": mwu_p,
    }).set_index("gene")
    df["spearman_fdr"] = bh_fdr(df["spearman_p"].to_numpy())
    df["welch_fdr"] = bh_fdr(df["welch_p"].to_numpy())
    df["up_in_trop2_high"] = (
        (df["spearman_rho"] > 0) & (df["log2FC_high_vs_low"] > 0)
        & (df["welch_fdr"] < FDR_THRESHOLD)
    )
    df["n_high"] = int(high_mask.sum())
    df["n_low"] = int(low_mask.sum())
    df["n_total"] = int(n)
    return df.sort_values("spearman_rho", ascending=False)


def run_prerank_gsea(stat_df, outdir, seed=7):
    """Pre-ranked GSEA on the Spearman-rho ranking against the local GMT."""
    import gseapy
    rnk = stat_df["spearman_rho"].dropna()
    rnk = rnk[rnk.index != TROP2]              # exclude the ranking gene itself
    rnk = rnk[~rnk.index.duplicated()].sort_values(ascending=False)
    res = gseapy.prerank(
        rnk=rnk.reset_index().rename(columns={"gene": 0, "spearman_rho": 1}),
        gene_sets=GMT, min_size=3, max_size=1000, permutation_num=1000,
        seed=seed, outdir=None, no_plot=True, threads=4,
    )
    tbl = res.res2d.copy()
    tbl.to_csv(f"{outdir}/gsea_prerank.csv", index=False)
    return tbl


def load_gmt(path=GMT):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = parts[2:]
    return sets
