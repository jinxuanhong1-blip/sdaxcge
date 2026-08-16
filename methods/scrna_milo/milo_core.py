"""Milo-equivalent neighborhood DA (Dann et al. 2022), without miloR/edgeR.

Neighborhood construction follows milopy/miloR:
  1. sample a fraction of graph vertices
  2. refine each sampled neighborhood to the cell nearest the median embedding
  3. neighborhood = KNN of the refined index (binary connectivities)

DA uses a shared-dispersion negative-binomial GLM + quasi-likelihood F test
on the nhood × sample count matrix, with log(sample library size) as offset.
SpatialFDR is the density-weighted BH procedure from cydar/miloR
(weight = 1 / k-th neighbor distance of the index cell).

This is not a drop-in of the Bioconductor binaries. The SpatialFDR formula
and the sampling scheme are the published ones; the GLM is a Python QLF
stand-in for edgeR::glmQLFTest.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse, stats


def make_nhoods(
    connectivities: sparse.spmatrix,
    distances: sparse.spmatrix,
    embedding: np.ndarray,
    n_neighbors: int,
    prop: float = 0.1,
    seed: int = 42,
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Return (cells × nhoods binary matrix, index cell ixs, kth distances)."""
    n_cells = connectivities.shape[0]
    knn = connectivities.tocsr().copy()
    knn.data = np.ones_like(knn.data, dtype=np.float32)

    rng = np.random.default_rng(seed)
    n_ixs = max(int(np.round(n_cells * prop)), 1)
    random_vertices = np.sort(rng.choice(n_cells, size=n_ixs, replace=False))

    refined = np.empty(len(random_vertices), dtype=np.int64)
    for i, v in enumerate(random_vertices):
        start, end = knn.indptr[v], knn.indptr[v + 1]
        cols = knn.indices[start:end]
        if cols.size == 0:
            refined[i] = v
            continue
        coords = embedding[cols]
        med = np.median(coords, axis=0)
        d2 = np.sum((coords - med) ** 2, axis=1)
        refined[i] = int(cols[int(np.argmin(d2))])
    index_ixs = np.unique(refined)

    nhoods = knn[:, index_ixs].tocsr()
    # include the index cell itself
    nhoods = nhoods.tolil()
    for j, ix in enumerate(index_ixs):
        nhoods[ix, j] = 1
    nhoods = nhoods.tocsr()
    nhoods.data = np.ones_like(nhoods.data, dtype=np.float32)

    dist = distances.tocsr()
    kth = np.zeros(len(index_ixs), dtype=np.float64)
    for j, ix in enumerate(index_ixs):
        start, end = dist.indptr[ix], dist.indptr[ix + 1]
        if end > start:
            kth[j] = float(dist.data[start:end].max())
        else:
            kth[j] = np.nan
    return nhoods, index_ixs, kth


def count_nhoods(
    nhoods: sparse.spmatrix,
    sample_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (nhoods × samples counts, sample names)."""
    samples = np.array(sorted(pd.unique(sample_ids)))
    sample_index = {s: i for i, s in enumerate(samples)}
    dummy = sparse.lil_matrix((len(sample_ids), len(samples)), dtype=np.float32)
    for i, s in enumerate(sample_ids):
        dummy[i, sample_index[s]] = 1
    dummy = dummy.tocsr()
    counts = nhoods.T.dot(dummy).toarray().astype(np.float64)
    return counts, samples


def graph_spatial_fdr(pvalues: np.ndarray, kth_distance: np.ndarray) -> np.ndarray:
    """Weighted BH SpatialFDR (cydar / miloR graphSpatialFDR)."""
    p = np.asarray(pvalues, dtype=np.float64)
    dist = np.asarray(kth_distance, dtype=np.float64)
    w = np.divide(1.0, dist, out=np.zeros_like(dist), where=dist > 0)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    keep = np.isfinite(p) & np.isfinite(w) & (w > 0)
    if keep.sum() == 0:
        return out
    p_k = p[keep]
    w_k = w[keep]
    order = np.argsort(p_k, kind="mergesort")
    p_s = p_k[order]
    w_s = w_k[order]
    q = (w_s.sum() * p_s) / np.maximum(np.cumsum(w_s), np.finfo(float).tiny)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    adj_pk = np.empty_like(p_s)
    adj_pk[order] = q
    out[keep] = adj_pk
    return out


def _design_matrix(condition: np.ndarray, alt_level: str) -> tuple[np.ndarray, list[str]]:
    cond = np.asarray(condition).astype(str)
    levels = [str(x) for x in pd.unique(cond)]
    if len(levels) != 2:
        raise ValueError(f"DA requires exactly 2 condition levels, got {levels}")
    if alt_level not in levels:
        raise ValueError(f"alt_level {alt_level!r} not in {levels}")
    ref = [x for x in levels if x != alt_level][0]
    x = np.column_stack(
        [
            np.ones(len(cond), dtype=np.float64),
            (cond == alt_level).astype(np.float64),
        ]
    )
    return x, [ref, alt_level]


def _irls_nb(
    y: np.ndarray,
    X: np.ndarray,
    offset: np.ndarray,
    dispersion: np.ndarray,
    max_iter: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized IRLS for log-link NB/Poisson across neighborhoods.

    y: (n_nhoods, n_samples)
    X: (n_samples, p)
    offset: (n_samples,)
    dispersion: (n_nhoods,)  NB variance = mu + disp * mu^2; 0 = Poisson
    """
    n_h, n_s = y.shape
    p = X.shape[1]
    log_y = np.log(y + 0.5) - offset
    xtx = X.T @ X
    try:
        beta = np.linalg.solve(xtx, X.T @ log_y.T).T
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(X, log_y.T, rcond=None)[0].T

    disp = np.asarray(dispersion, dtype=np.float64).reshape(-1, 1)
    for _ in range(max_iter):
        eta = beta @ X.T + offset
        mu = np.exp(np.clip(eta, -20.0, 20.0))
        var = np.maximum(mu + disp * mu**2, 1e-8)
        w = (mu**2) / var
        z = (eta - offset) + (y - mu) / np.maximum(mu, 1e-8)
        xtwx = np.einsum("sp,hs,sq->hpq", X, w, X)
        xtwz = np.einsum("sp,hs,hs->hp", X, w, z)
        # ridge for singular nhoods
        xtwx = xtwx + 1e-6 * np.eye(p)[None, :, :]
        beta = np.linalg.solve(xtwx, xtwz[..., None]).squeeze(-1)
    eta = beta @ X.T + offset
    mu = np.exp(np.clip(eta, -20.0, 20.0))
    return beta, mu


def _nb_loglik(y: np.ndarray, mu: np.ndarray, disp: np.ndarray) -> np.ndarray:
    """Per-nhood NB log-likelihood (up to a constant in y). disp=0 → Poisson."""
    mu = np.maximum(mu, 1e-8)
    d = np.asarray(disp, dtype=np.float64).reshape(-1, 1)
    poisson = np.sum(y * np.log(mu) - mu, axis=1)
    # NB: r = 1/disp, p = r / (r+mu)
    r = np.where(d > 1e-8, 1.0 / np.maximum(d, 1e-8), np.inf)
    finite = np.isfinite(r[:, 0])
    ll = poisson.copy()
    if finite.any():
        rr = r[finite]
        mm = mu[finite]
        yy = y[finite]
        ll[finite] = np.sum(
            stats.nbinom.logpmf(np.clip(yy, 0, None).astype(int), rr, rr / (rr + mm)),
            axis=1,
        )
    return ll


def da_nhoods_qlf(
    counts: np.ndarray,
    condition: np.ndarray,
    lib_size: np.ndarray,
    kth_distance: np.ndarray,
    alt_level: str,
    min_total: int = 20,
    min_samples: int = 3,
) -> pd.DataFrame:
    """NB-QLF DA for a two-level condition. Returns one row per neighborhood."""
    y = np.asarray(counts, dtype=np.float64)
    n_h, n_s = y.shape
    cond = np.asarray(condition)
    lib = np.asarray(lib_size, dtype=np.float64)
    keep_s = lib > 0
    y = y[:, keep_s]
    cond = cond[keep_s]
    lib = lib[keep_s]
    offset = np.log(np.maximum(lib, 1.0))

    n_per = (y > 0).sum(axis=1)
    tot = y.sum(axis=1)
    testable = (tot >= min_total) & (n_per >= min_samples)
    # both condition levels must appear
    levels = pd.unique(cond.astype(str))
    if len(levels) != 2:
        raise ValueError(f"need 2 condition levels after lib-size filter, got {list(levels)}")

    X, level_names = _design_matrix(cond.astype(str), alt_level=alt_level)
    X_red = X[:, :1]

    y_fit = y.copy()
    y_fit[~testable] = 0.0

    # Poisson IRLS → method-of-moments common + tagwise dispersion
    beta0, mu0 = _irls_nb(y, X, offset, dispersion=np.zeros(n_h), max_iter=12)
    resid2 = (y - mu0) ** 2
    mom = (resid2 - mu0) / np.maximum(mu0**2, 1e-8)
    tagwise = np.clip(np.nanmedian(mom, axis=1), 1e-6, 20.0)
    tagwise[~testable] = np.nan
    common = float(np.nanmedian(tagwise)) if np.isfinite(tagwise).any() else 0.1
    # edgeR-like squeeze toward common (moderate prior)
    disp = np.where(testable, (tagwise + common) / 2.0, common)

    beta, mu = _irls_nb(y, X, offset, dispersion=disp, max_iter=15)
    beta_r, mu_r = _irls_nb(y, X_red, offset, dispersion=disp, max_iter=15)

    ll_f = _nb_loglik(y, mu, disp)
    ll_r = _nb_loglik(y, mu_r, disp)
    lrt = np.clip(2.0 * (ll_f - ll_r), 0.0, None)

    # QL dispersion from Pearson residuals of the full NB mean
    var = np.maximum(mu + disp.reshape(-1, 1) * mu**2, 1e-8)
    pearson = np.sum((y - mu) ** 2 / var, axis=1)
    df_res = max(y.shape[1] - X.shape[1], 1)
    ql_disp = np.maximum(pearson / df_res, 1e-4)
    f_stat = (lrt / 1.0) / ql_disp
    pval = stats.f.sf(f_stat, 1, df_res)

    logfc = beta[:, 1] / np.log(2.0)  # log2 fold change of the alt level
    # logCPM-like: mean log2 (count/lib * 1e6)
    cpm = np.log2(y / lib[None, :] * 1e6 + 1.0).mean(axis=1)

    pval = np.where(testable, pval, np.nan)
    f_stat = np.where(testable, f_stat, np.nan)
    lrt = np.where(testable, lrt, np.nan)
    logfc = np.where(testable, logfc, np.nan)

    spatial = graph_spatial_fdr(pval, kth_distance)
    # ordinary BH for comparison (not the Milo quantity)
    bh = np.full(n_h, np.nan)
    if testable.any():
        from statsmodels.stats.multitest import multipletests

        _, bh_t, _, _ = multipletests(pval[testable], method="fdr_bh")
        bh[testable] = bh_t

    return pd.DataFrame(
        {
            "nhood": np.arange(n_h),
            "n_cells": tot,
            "n_samples_nonzero": n_per,
            "testable": testable,
            "logFC": logfc,
            "logCPM": cpm,
            "F": f_stat,
            "PValue": pval,
            "FDR_BH": bh,
            "SpatialFDR": spatial,
            "dispersion": disp,
            "ql_dispersion": ql_disp,
            "ref_level": level_names[0],
            "alt_level": level_names[1],
            "n_samples": n_s,
            "n_ref": int((cond.astype(str) == level_names[0]).sum()),
            "n_alt": int((cond.astype(str) == level_names[1]).sum()),
        }
    )


@dataclass
class NhoodDA:
    embedding: str
    contrast: str
    results: pd.DataFrame
    counts: np.ndarray
    samples: np.ndarray
    index_ixs: np.ndarray
    kth_distance: np.ndarray
