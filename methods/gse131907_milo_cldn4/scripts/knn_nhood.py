"""Documented kNN-neighborhood abundance test (miloR fallback).

miloR / edgeR are not required. This module reimplements the *graph and
multiple-testing* pieces of Milo (Dann et al., Nat Biotechnol 2022) and
replaces the edgeR quasi-likelihood DA GLM with sample-level Welch t-tests
or Spearman tests. Sample is the independent unit.

SpatialFDR is miloR::graphSpatialFDR with weighting='k-distance'
(Dann et al. 2022; Lun et al. cydar):

    w = 1 / k_distance
    order p ascending; carry w in that order
    adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
    clip to [0, 1]

This is not a claim that the Python numbers match miloR bit-for-bit
(index sampling and the DA model differ). Report both BH-FDR and
SpatialFDR, and the neighbourhood / sample counts actually tested.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

EPS = 1e-8


def spatial_fdr_kdistance(pvalues: np.ndarray, k_distance: np.ndarray) -> np.ndarray:
    """miloR::graphSpatialFDR, weighting='k-distance' (Dann / cydar)."""
    p = np.asarray(pvalues, dtype=float)
    d = np.asarray(k_distance, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p) & np.isfinite(d) & (d > 0)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    w = 1.0 / d[ok]
    w[~np.isfinite(w)] = 1.0
    o = np.argsort(pv, kind="mergesort")
    pv_o = pv[o]
    w_o = w[o]
    raw = w_o.sum() * pv_o / np.cumsum(w_o)
    adj_o = np.minimum(1.0, np.minimum.accumulate(raw[::-1])[::-1])
    adj = np.empty_like(adj_o)
    adj[o] = adj_o
    out[ok] = adj
    return out


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = pv.size
    o = np.argsort(pv, kind="mergesort")
    ranked = pv[o]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum(1.0, np.minimum.accumulate(q[::-1])[::-1])
    adj = np.empty_like(q)
    adj[o] = q
    out[ok] = adj
    return out


def select_hvg(
    gene_mean: np.ndarray,
    gene_var: np.ndarray,
    n_hvg: int = 2000,
    min_mean: float = 0.01,
    max_mean: float = 50.0,
) -> np.ndarray:
    """Dispersion rank among genes with intermediate mean UMI."""
    mu = np.asarray(gene_mean, dtype=float)
    va = np.asarray(gene_var, dtype=float)
    ok = np.isfinite(mu) & np.isfinite(va) & (mu >= min_mean) & (mu <= max_mean) & (va > 0)
    disp = np.full(mu.shape, -np.inf, dtype=float)
    disp[ok] = va[ok] / mu[ok]
    n_keep = min(int(n_hvg), int(ok.sum()))
    if n_keep <= 0:
        raise RuntimeError("no genes passed HVG mean filter")
    return np.argpartition(-disp, n_keep - 1)[:n_keep]


def pca_knn(
    log_cp10k: np.ndarray,
    sample: np.ndarray,
    n_pcs: int = 30,
    k: int = 30,
    random_state: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (pca_centered, knn_idx [n,k] excluding self, knn_dist)."""
    z = StandardScaler(with_mean=True, with_std=True).fit_transform(log_cp10k)
    n_pcs = min(n_pcs, z.shape[0] - 1, z.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="randomized", random_state=random_state).fit_transform(z)
    for s in np.unique(sample):
        m = sample == s
        if m.sum() > 1:
            pcs[m] -= pcs[m].mean(axis=0)
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean")
    nn.fit(pcs)
    dist, idx = nn.kneighbors(pcs)
    return pcs.astype(np.float32), idx[:, 1:].astype(np.int32), dist[:, 1:].astype(np.float32)


def refine_indices(
    pcs: np.ndarray,
    knn_idx: np.ndarray,
    prop: float = 0.1,
    random_state: int = 1,
) -> np.ndarray:
    """Milo refined sampling: random seeds → kNN centroid → nearest member."""
    rng = np.random.default_rng(random_state)
    n = pcs.shape[0]
    n_seed = max(1, int(round(prop * n)))
    seeds = rng.choice(n, size=n_seed, replace=False)
    refined = np.empty(n_seed, dtype=np.int32)
    for i, s in enumerate(seeds):
        memb = knn_idx[s]
        centroid = pcs[memb].mean(axis=0)
        d = np.linalg.norm(pcs[memb] - centroid, axis=1)
        refined[i] = memb[int(np.argmin(d))]
    return np.unique(refined)


@dataclass
class Nhoods:
    index: np.ndarray
    members: np.ndarray
    k_distance: np.ndarray
    k: int


def make_nhoods(knn_idx: np.ndarray, knn_dist: np.ndarray, indices: np.ndarray) -> Nhoods:
    k = knn_idx.shape[1]
    return Nhoods(
        index=np.asarray(indices, dtype=np.int32),
        members=knn_idx[indices],
        k_distance=knn_dist[indices, k - 1],
        k=k,
    )


def nhood_membership(nhoods: Nhoods) -> list[np.ndarray]:
    out = []
    for i, idx in enumerate(nhoods.index):
        out.append(np.unique(np.concatenate(([idx], nhoods.members[i]))))
    return out


def count_matrix(members: list[np.ndarray], sample: np.ndarray, sample_levels: list[str]) -> np.ndarray:
    smap = {s: j for j, s in enumerate(sample_levels)}
    codes = np.array([smap[s] for s in sample], dtype=np.int32)
    mat = np.zeros((len(members), len(sample_levels)), dtype=np.int32)
    for i, cells in enumerate(members):
        for j in codes[cells]:
            mat[i, j] += 1
    return mat


def sample_sizes(sample: np.ndarray, sample_levels: list[str]) -> np.ndarray:
    return np.array([(sample == s).sum() for s in sample_levels], dtype=np.float64)


def welch_da(
    counts: np.ndarray,
    sizes: np.ndarray,
    group: np.ndarray,
    group_a: str,
    group_b: str,
    min_per_group: int = 2,
    min_samples: int = 3,
) -> pd.DataFrame:
    """Per-nhood Welch t-test of sample proportions. logFC = B − A on log2 mean prop."""
    prop = counts / np.maximum(sizes, 1.0)
    a = group == group_a
    b = group == group_b
    rows = []
    for i in range(prop.shape[0]):
        pa = prop[i, a]
        pb = prop[i, b]
        n_present = int(((counts[i] > 0) & (a | b)).sum())
        n_a = int((counts[i, a] > 0).sum())
        n_b = int((counts[i, b] > 0).sum())
        rec = {
            "nhood": i,
            "n_samples_present": n_present,
            "n_samples_A": n_a,
            "n_samples_B": n_b,
            "n_cells": int(counts[i].sum()),
            "mean_prop_A": float(pa.mean()) if a.any() else np.nan,
            "mean_prop_B": float(pb.mean()) if b.any() else np.nan,
            "logFC_B_minus_A": np.nan,
            "t": np.nan,
            "p": np.nan,
            "testable": False,
        }
        if n_a >= min_per_group and n_b >= min_per_group and n_present >= min_samples:
            rec["testable"] = True
            rec["logFC_B_minus_A"] = float(np.log2(pb.mean() + EPS) - np.log2(pa.mean() + EPS))
            if pa.std(ddof=1) == 0 and pb.std(ddof=1) == 0:
                rec["t"] = 0.0
                rec["p"] = 1.0
            else:
                t, p = stats.ttest_ind(pb, pa, equal_var=False)
                rec["t"] = float(t)
                rec["p"] = float(p)
        rows.append(rec)
    return pd.DataFrame(rows)


def spearman_da(
    counts: np.ndarray,
    sizes: np.ndarray,
    covariate: np.ndarray,
    min_samples: int = 5,
) -> pd.DataFrame:
    """Per-nhood Spearman of sample proportion vs a sample-level covariate."""
    prop = counts / np.maximum(sizes, 1.0)
    rows = []
    for i in range(prop.shape[0]):
        y = prop[i]
        x = covariate
        m = np.isfinite(x) & np.isfinite(y) & (counts[i] > 0)
        rec = {
            "nhood": i,
            "n_samples_present": int(m.sum()),
            "n_cells": int(counts[i].sum()),
            "spearman_rho": np.nan,
            "p": np.nan,
            "testable": False,
        }
        if m.sum() >= min_samples and np.unique(x[m]).size >= 3:
            rec["testable"] = True
            rho, p = stats.spearmanr(x[m], y[m])
            rec["spearman_rho"] = float(rho)
            rec["p"] = float(p)
        rows.append(rec)
    return pd.DataFrame(rows)


def attach_fdr(df: pd.DataFrame, k_distance: np.ndarray, p_col: str = "p") -> pd.DataFrame:
    out = df.copy()
    p = out[p_col].to_numpy(dtype=float)
    out["BH_FDR"] = bh_fdr(p)
    out["SpatialFDR"] = spatial_fdr_kdistance(p, k_distance)
    return out


def nhood_composition(
    members: list[np.ndarray],
    is_tnk: np.ndarray,
    is_malig: np.ndarray,
    gene_log: np.ndarray,
    lineage: np.ndarray,
    gene_prefix: str = "cldn4",
) -> pd.DataFrame:
    rows = []
    for i, cells in enumerate(members):
        tnk = is_tnk[cells]
        mal = is_malig[cells]
        rows.append(
            {
                "nhood": i,
                "n_cells": int(cells.size),
                "n_tnk": int(tnk.sum()),
                "n_malig": int(mal.sum()),
                "frac_tnk": float(tnk.mean()),
                "frac_malig": float(mal.mean()),
                f"{gene_prefix}_all_mean": float(gene_log[cells].mean()),
                f"{gene_prefix}_malig_mean": float(gene_log[cells[mal]].mean()) if mal.any() else np.nan,
                "dominant_lineage": _mode(lineage[cells]),
            }
        )
    return pd.DataFrame(rows)


def _mode(x: np.ndarray) -> str:
    vals, counts = np.unique(x.astype(str), return_counts=True)
    return str(vals[counts.argmax()])


def greedy_independent(members: list[np.ndarray], max_shared: int = 0) -> np.ndarray:
    """Greedy low-overlap subset. max_shared=0 → disjoint neighbourhoods."""
    order = np.argsort([-m.size for m in members])
    chosen: list[int] = []
    covered: set[int] = set()
    for i in order:
        cells = members[i]
        shared = sum(1 for c in cells if c in covered)
        if shared <= max_shared:
            chosen.append(int(i))
            covered.update(int(c) for c in cells)
    return np.array(sorted(chosen), dtype=int)


def sample_paired_tnk_by_gene(
    members: list[np.ndarray],
    sample: np.ndarray,
    is_tnk: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
) -> pd.DataFrame:
    """Per-sample T/NK fraction among cells that sit in gene-high vs low nhoods."""
    sample_levels = list(pd.unique(sample))
    rows = []
    for s in sample_levels:
        high_cells: set[int] = set()
        low_cells: set[int] = set()
        for i, cells in enumerate(members):
            if high[i]:
                high_cells.update(int(c) for c in cells if sample[c] == s)
            if low[i]:
                low_cells.update(int(c) for c in cells if sample[c] == s)

        def frac(ids: set[int]) -> float:
            if not ids:
                return np.nan
            arr = np.fromiter(ids, dtype=np.int32)
            return float(is_tnk[arr].mean())

        rows.append(
            {
                "sample": s,
                "n_cells_high_nhood": len(high_cells),
                "n_cells_low_nhood": len(low_cells),
                "frac_tnk_high_nhood": frac(high_cells),
                "frac_tnk_low_nhood": frac(low_cells),
            }
        )
    return pd.DataFrame(rows)


def spearman_safe(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return {"n": int(m.sum()), "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": int(m.sum()), "rho": float(rho), "p": float(p)}


def wilcoxon_paired(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    rec = {"n": int(m.sum()), "median_a": None, "median_b": None, "p": None}
    if m.sum() < 3:
        return rec
    rec["median_a"] = float(np.median(a[m]))
    rec["median_b"] = float(np.median(b[m]))
    try:
        rec["p"] = float(stats.wilcoxon(a[m], b[m], alternative="two-sided").pvalue)
    except ValueError:
        rec["p"] = None
    return rec


def write_json(path, obj) -> None:
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_json_default)


def _json_default(x):
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    raise TypeError(type(x))


def self_test() -> None:
    """Tiny numerical check of SpatialFDR monotonicity and BH on independent p."""
    rng = np.random.default_rng(0)
    p = rng.uniform(0, 1, 200)
    d = rng.uniform(0.2, 2.0, 200)
    q = spatial_fdr_kdistance(p, d)
    assert np.all((q >= 0) | ~np.isfinite(q))
    assert np.all((q <= 1) | ~np.isfinite(q))
    q2 = spatial_fdr_kdistance(p, np.ones_like(p))
    bh = bh_fdr(p)
    if not np.allclose(q2, bh, atol=1e-8, equal_nan=True):
        raise AssertionError("equal-weight SpatialFDR should match BH")
    print("knn_nhood.self_test OK")


if __name__ == "__main__":
    self_test()
