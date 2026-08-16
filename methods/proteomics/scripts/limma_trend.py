#!/usr/bin/env python3
"""limma / limma-trend (Smyth eBayes) for already-logged TMT or log-CPM RNA.

Use
---
* TMT log2-ratio (CPTAC/LinkedOmics): ``--method limma`` (optionally ``--trend``).
  Do **not** run voom on these values.
* RNA raw counts: convert to log2-CPM (+ prior) and use ``--method limma-trend``.
* RNA already log2(FPKM/UQ): ``--method limma``; still not voom.

This is a compact Python port of the Smyth (2004) moderated t + optional
mean–variance trend. For production, prefer Bioconductor ``limma``
(``limma_trend.R`` in this folder).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.special import digamma, polygamma

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib_io import apply_log2_if_needed, normalize_sample_id, read_expression_matrix  # noqa: E402


def counts_to_logcpm(counts: pd.DataFrame, prior: float = 2.0) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + prior)


def _design(sample_ids: list[str], group: pd.Series) -> tuple[np.ndarray, list[str]]:
    g = group.reindex(sample_ids)
    if g.isna().any():
        missing = g[g.isna()].index.tolist()
        raise SystemExit(f"group labels missing for: {missing[:8]}")
    levels = list(pd.Index(g.astype(str).unique()))
    X = np.column_stack([(g.astype(str) == lev).astype(float) for lev in levels])
    return X, levels


def _fit_ols(Y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Y: genes x samples. Returns beta (genes x p), resid_ss, df, unscaled_se factor."""
    n, p = X.shape
    xtx = X.T @ X
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    hat = xtx_inv @ X.T
    beta = (hat @ Y.T).T  # genes x p
    fitted = beta @ X.T
    resid = Y - fitted
    df = float(n - np.linalg.matrix_rank(X))
    rss = np.sum(resid**2, axis=1)
    return beta, rss, df, xtx_inv


def _fit_f_dist(s2: np.ndarray, df: float) -> tuple[float, float]:
    """Moment estimators for scaled-F as in limma::fitFDist (simplified)."""
    s2 = np.asarray(s2, dtype=float)
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    if s2.size < 10:
        return 0.0, float(np.median(s2) if s2.size else 1.0)
    z = np.log(s2)
    e = float(z.mean())
    v = float(z.var(ddof=1))
    # trigamma(df/2) correction roughly as limma
    trigamma_df = float(polygamma(1, df / 2.0)) if df > 2 else 1.0
    evar = v - trigamma_df
    if evar <= 0:
        d0 = np.inf
    else:
        # solve trigamma(d0/2) ≈ evar
        def obj(d0: float) -> float:
            return (float(polygamma(1, d0 / 2.0)) - evar) ** 2

        res = minimize_scalar(obj, bounds=(1e-2, 1e4), method="bounded")
        d0 = float(res.x)
    digamma_df = float(digamma(df / 2.0)) if df > 0 else 0.0
    if np.isfinite(d0):
        s0 = float(np.exp(e - digamma(d0 / 2.0) + digamma_df))
    else:
        s0 = float(np.exp(e + digamma_df))
    return d0, s0


def _loess_trend(mean_expr: np.ndarray, s2: np.ndarray, frac: float = 0.3) -> np.ndarray:
    """Lowess mean–variance trend; fallback to running median if needed."""
    order = np.argsort(mean_expr)
    x = mean_expr[order]
    y = np.log(np.clip(s2[order], 1e-12, None))
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        sm = lowess(y, x, frac=frac, return_sorted=True)
        pred = np.interp(mean_expr, sm[:, 0], sm[:, 1])
    except Exception:
        win = max(11, int(len(y) * frac) | 1)
        ser = pd.Series(y).rolling(win, center=True, min_periods=5).median()
        ser = ser.bfill().ffill()
        pred = np.empty_like(mean_expr)
        pred[order] = ser.to_numpy()
        return np.exp(pred)
    return np.exp(pred)


def ebayes(
    Y: pd.DataFrame,
    X: np.ndarray,
    coef: int,
    trend: bool = False,
) -> pd.DataFrame:
    beta, rss, df, xtx_inv = _fit_ols(Y.to_numpy(dtype=float), X)
    s2 = rss / df if df > 0 else np.full(Y.shape[0], np.nan)
    mean_expr = np.nanmean(Y.to_numpy(dtype=float), axis=1)
    if trend:
        s2_trend = _loess_trend(mean_expr, s2)
        # residual about the trend, then squeeze
        ratio = s2 / np.clip(s2_trend, 1e-12, None)
        d0, s0 = _fit_f_dist(ratio, df)
        if np.isfinite(d0):
            s2_post = ((d0 * s0 * s2_trend) + (df * s2)) / (d0 + df)
        else:
            s2_post = s2_trend
    else:
        d0, s0 = _fit_f_dist(s2, df)
        if np.isfinite(d0):
            s2_post = (d0 * s0 + df * s2) / (d0 + df)
        else:
            s2_post = np.full_like(s2, s0)

    var_unscaled = float(xtx_inv[coef, coef])
    se = np.sqrt(np.clip(s2_post, 1e-18, None) * var_unscaled)
    tstat = beta[:, coef] / se
    df_total = df + (0.0 if not np.isfinite(d0) else d0)
    pval = 2 * stats.t.sf(np.abs(tstat), df_total)
    out = pd.DataFrame(
        {
            "gene": Y.index,
            "logFC": beta[:, coef],
            "AveExpr": mean_expr,
            "t": tstat,
            "P.Value": pval,
            "s2": s2,
            "s2_post": s2_post,
            "df_total": df_total,
        }
    )
    out["adj.P.Val"] = _bh(out["P.Value"])
    return out.sort_values("P.Value")


def _bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    n = p.notna().sum()
    order = p.argsort()
    ranked = pd.Series(np.nan, index=p.index)
    pv = p.iloc[order].to_numpy()
    adj = np.empty_like(pv)
    prev = 1.0
    for i in range(len(pv) - 1, -1, -1):
        if np.isnan(pv[i]):
            adj[i] = np.nan
            continue
        val = min(prev, pv[i] * n / (i + 1))
        adj[i] = val
        prev = val
    ranked.iloc[order] = adj
    return ranked.clip(upper=1.0)


def align_group(matrix: pd.DataFrame, group_tsv: Path, sample_col: str, group_col: str) -> pd.Series:
    g = pd.read_csv(group_tsv, sep="\t")
    g[sample_col] = g[sample_col].map(normalize_sample_id)
    series = g.set_index(sample_col)[group_col]
    cols = [normalize_sample_id(c) for c in matrix.columns]
    matrix.columns = cols
    keep = [c for c in cols if c in series.index]
    return series.loc[keep]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--group-tsv", type=Path, required=True, help="sample_id + group columns")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--group-col", default="group")
    ap.add_argument("--contrast", nargs=2, default=["Tumor", "NAT"], metavar=("NUM", "DEN"))
    ap.add_argument("--method", choices=["limma", "limma-trend"], default="limma")
    ap.add_argument("--counts", action="store_true", help="Treat matrix as raw counts; convert to log-CPM")
    ap.add_argument("--min-obs", type=float, default=0.7, help="Keep genes observed in this fraction of samples")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    mat = read_expression_matrix(args.matrix)
    mat.columns = [normalize_sample_id(c) for c in mat.columns]
    if args.counts:
        mat = counts_to_logcpm(mat.fillna(0))
        log_info = {"already_log": True, "reason": "converted counts -> log2-CPM"}
    else:
        mat, log_info = apply_log2_if_needed(mat, name=args.matrix.name)
        if not log_info.get("already_log", True):
            print("applied log2(x+1); for counts prefer --counts (limma-trend on log-CPM)", file=sys.stderr)

    group = align_group(mat, args.group_tsv, args.sample_col, args.group_col)
    mat = mat.loc[:, group.index]
    keep = mat.notna().mean(axis=1) >= args.min_obs
    mat = mat.loc[keep].copy()
    # gene-median impute for limma (document; better than dropping every incomplete gene)
    mat = mat.apply(lambda s: s.fillna(s.median()), axis=1)

    X, levels = _design(list(mat.columns), group)
    num, den = args.contrast
    if num not in levels or den not in levels:
        raise SystemExit(f"contrast levels {args.contrast} not in {levels}")
    # coef of num - den via contrast on treatment-coded intercept-free design
    # beta_num - beta_den
    c = np.zeros(len(levels))
    c[levels.index(num)] = 1
    c[levels.index(den)] = -1
    # rewrite as single coefficient by rotating design
    # Use OLS on contrast: fit full, then contrast
    beta, rss, df, xtx_inv = _fit_ols(mat.to_numpy(dtype=float), X)
    # attach contrast as extra column by evaluating c
    # Reuse ebayes on a 2-level subset if both groups present
    mask = group.astype(str).isin(args.contrast)
    Y = mat.loc[:, mask.index[mask]]
    g2 = group.loc[Y.columns]
    X2, lev2 = _design(list(Y.columns), g2)
    # put NUM last so coef=-1 is DEN, coef=1 is NUM if we use NUM as second dummy
    # Force order [DEN, NUM] so coef 1 is NUM vs DEN when using treatment contrast.
    # Our design is one-hot without intercept, so coef_NUM - coef_DEN.
    trend = args.method == "limma-trend"
    # Compute full one-hot eBayes then contrast
    fit_num = ebayes(Y, X2, coef=lev2.index(num), trend=trend)
    fit_den = ebayes(Y, X2, coef=lev2.index(den), trend=trend)
    # logFC = Ave(num) - Ave(den) ≈ beta_num - beta_den; rebuild from stored logFC
    merged = fit_num.set_index("gene")[["logFC", "AveExpr", "s2", "s2_post", "df_total"]].rename(
        columns={"logFC": "beta_num"}
    )
    merged["beta_den"] = fit_den.set_index("gene")["logFC"]
    merged["logFC"] = merged["beta_num"] - merged["beta_den"]
    # SE of difference from posterior variances (conservative, independent-ish)
    # Better: use contrast variance c' (X'X)^-1 c * s2_post
    _, _, _, xtx2 = _fit_ols(Y.to_numpy(dtype=float), X2)
    cvec = np.zeros(len(lev2))
    cvec[lev2.index(num)] = 1
    cvec[lev2.index(den)] = -1
    var_unscaled = float(cvec @ xtx2 @ cvec)
    se = np.sqrt(np.clip(merged["s2_post"].to_numpy(), 1e-18, None) * var_unscaled)
    tstat = merged["logFC"].to_numpy() / se
    pval = 2 * stats.t.sf(np.abs(tstat), merged["df_total"].to_numpy())
    out = pd.DataFrame(
        {
            "gene": merged.index,
            "logFC": merged["logFC"].to_numpy(),
            "AveExpr": merged["AveExpr"].to_numpy(),
            "t": tstat,
            "P.Value": pval,
            "method": args.method,
            "contrast": f"{num}-{den}",
        }
    )
    out["adj.P.Val"] = _bh(out["P.Value"])
    out = out.sort_values("P.Value")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}  n_genes={out.shape[0]}  method={args.method}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
