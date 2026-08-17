"""Street et al. 2018 Slingshot — Python equivalent.

Two official steps, no diffusion-pseudotime fallback:

1. getLineages: MST on cluster centroids in the reduced space.
2. getCurves: Hastie–Stuetzle principal curves along each lineage,
   initialized from the centroid polyline.

This is the clock used when R/Bioconductor slingshot is unavailable.
It is not scanpy DPT.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist


def cluster_centroids(embedding: np.ndarray, labels: np.ndarray) -> dict[str, np.ndarray]:
    cents: dict[str, np.ndarray] = {}
    for lab in pd.unique(labels):
        mask = labels == lab
        if int(mask.sum()) == 0:
            continue
        cents[str(lab)] = np.asarray(embedding[mask], dtype=float).mean(axis=0)
    return cents


def _mst_edges(cents: dict[str, np.ndarray]) -> list[tuple[str, str, float]]:
    labs = list(cents.keys())
    if len(labs) < 2:
        return []
    mat = np.vstack([cents[k] for k in labs])
    dist = cdist(mat, mat, metric="euclidean")
    np.fill_diagonal(dist, 0.0)
    mst = minimum_spanning_tree(dist).toarray()
    edges: list[tuple[str, str, float]] = []
    n = len(labs)
    for i in range(n):
        for j in range(n):
            if mst[i, j] > 0:
                edges.append((labs[i], labs[j], float(mst[i, j])))
            elif mst[j, i] > 0 and i < j:
                edges.append((labs[i], labs[j], float(mst[j, i])))
    return edges


def _adjacency(edges: list[tuple[str, str, float]]) -> dict[str, list[str]]:
    adj: dict[str, list[str]] = {}
    for a, b, _ in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    return adj


def _path(adj: dict[str, list[str]], start: str, end: str) -> list[str] | None:
    if start == end:
        return [start]
    seen = {start}
    stack = [(start, [start])]
    while stack:
        node, path = stack.pop()
        for nxt in adj.get(node, []):
            if nxt in seen:
                continue
            if nxt == end:
                return path + [nxt]
            seen.add(nxt)
            stack.append((nxt, path + [nxt]))
    return None


def get_lineages(
    embedding: np.ndarray,
    labels: np.ndarray,
    start_cluster: str,
) -> dict[str, Any]:
    """MST lineages from the start cluster to every leaf."""
    labels = np.asarray(labels).astype(str)
    cents = cluster_centroids(embedding, labels)
    if start_cluster not in cents:
        raise ValueError(f"start cluster {start_cluster!r} missing from centroids")
    edges = _mst_edges(cents)
    adj = _adjacency(edges)
    for lab in cents:
        adj.setdefault(lab, [])
    degrees = {k: len(v) for k, v in adj.items()}
    leaves = [k for k, d in degrees.items() if d <= 1 and k != start_cluster]
    if not leaves:
        others = [k for k in cents if k != start_cluster]
        leaves = others
    lineages: list[list[str]] = []
    for leaf in sorted(leaves, key=lambda x: (len(_path(adj, start_cluster, x) or []), x)):
        path = _path(adj, start_cluster, leaf)
        if path and len(path) >= 2:
            lineages.append(path)
    if not lineages and len(cents) >= 2:
        # degenerate: one other cluster
        other = next(k for k in cents if k != start_cluster)
        lineages.append([start_cluster, other])
    return {
        "centroids": {k: v.tolist() for k, v in cents.items()},
        "edges": [{"a": a, "b": b, "dist": d} for a, b, d in edges],
        "start_cluster": start_cluster,
        "leaves": leaves,
        "lineages": lineages,
        "n_clusters": int(len(cents)),
        "n_lineages": int(len(lineages)),
        "method": "slingshot.getLineages MST (Street 2018 Python equivalent)",
    }


def _polyline_project(X: np.ndarray, curve: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project points onto a piecewise-linear curve.

    Returns (projected_points, arc_length, segment_index).
    """
    segs = curve[1:] - curve[:-1]
    seg_len = np.linalg.norm(segs, axis=1)
    seg_len = np.maximum(seg_len, 1e-12)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    n = X.shape[0]
    best_t = np.empty(n, dtype=float)
    best_proj = np.empty_like(X, dtype=float)
    best_seg = np.empty(n, dtype=int)
    best_d2 = np.full(n, np.inf)
    for i, (p0, vec, length) in enumerate(zip(curve[:-1], segs, seg_len)):
        u = vec / length
        rel = X - p0
        alpha = np.clip(rel @ u, 0.0, length)
        proj = p0 + np.outer(alpha, u)
        d2 = np.sum((X - proj) ** 2, axis=1)
        better = d2 < best_d2
        best_d2[better] = d2[better]
        best_proj[better] = proj[better]
        best_t[better] = cum[i] + alpha[better]
        best_seg[better] = i
    return best_proj, best_t, best_seg


def _resample_curve(proj: np.ndarray, t: np.ndarray, n_knots: int) -> np.ndarray:
    order = np.argsort(t)
    t_s = t[order]
    p_s = proj[order]
    # unique-ish quantiles along occupied arc
    lo, hi = float(t_s[0]), float(t_s[-1])
    if hi <= lo:
        return np.vstack([p_s[0], p_s[-1]])
    qs = np.linspace(lo, hi, n_knots)
    knots = []
    for q in qs:
        w = np.exp(-((t_s - q) ** 2) / (2.0 * ((hi - lo) / max(n_knots, 2)) ** 2 + 1e-12))
        w = w / w.sum()
        knots.append(w @ p_s)
    curve = np.asarray(knots, dtype=float)
    # drop near-duplicates
    keep = [0]
    for i in range(1, len(curve)):
        if np.linalg.norm(curve[i] - curve[keep[-1]]) > 1e-8:
            keep.append(i)
    return curve[keep]


def fit_principal_curve(
    X: np.ndarray,
    init: np.ndarray,
    n_knots: int = 24,
    n_iter: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """Hastie–Stuetzle principal curve, initialized from a centroid polyline."""
    curve = np.asarray(init, dtype=float)
    if curve.shape[0] < 2:
        curve = np.vstack([curve, curve + 1e-6])
    t = np.zeros(X.shape[0], dtype=float)
    for _ in range(n_iter):
        proj, t, _ = _polyline_project(X, curve)
        curve = _resample_curve(proj, t, n_knots=n_knots)
        if curve.shape[0] < 2:
            break
    # final projection
    _, t, _ = _polyline_project(X, curve)
    # normalize to [0, 1] on this lineage
    tmin, tmax = float(np.min(t)), float(np.max(t))
    if tmax > tmin:
        t = (t - tmin) / (tmax - tmin)
    else:
        t = np.zeros_like(t)
    return curve, t


def get_curves(
    embedding: np.ndarray,
    labels: np.ndarray,
    lineage_info: dict[str, Any],
    n_knots: int = 24,
    n_iter: int = 12,
) -> dict[str, Any]:
    """Fit one principal curve per MST lineage. Cells outside the lineage are NA."""
    labels = np.asarray(labels).astype(str)
    cents = {k: np.asarray(v, dtype=float) for k, v in lineage_info["centroids"].items()}
    lineages = lineage_info["lineages"]
    n = embedding.shape[0]
    times = np.full((n, len(lineages)), np.nan, dtype=float)
    weights = np.zeros((n, len(lineages)), dtype=float)
    curves: list[list[list[float]]] = []
    lineage_n: list[int] = []
    for j, path in enumerate(lineages):
        member = np.isin(labels, path)
        n_mem = int(member.sum())
        lineage_n.append(n_mem)
        if n_mem < 8:
            curves.append([cents[c].tolist() for c in path if c in cents])
            continue
        init = np.vstack([cents[c] for c in path if c in cents])
        curve, t = fit_principal_curve(
            embedding[member], init, n_knots=n_knots, n_iter=n_iter
        )
        times[member, j] = t
        weights[member, j] = 1.0
        curves.append(curve.tolist())
    # shared clock: mean of assigned lineage times (Street: cells can sit on several)
    n_assigned = np.sum(np.isfinite(times), axis=1)
    shared = np.nanmean(times, axis=1)
    shared[n_assigned == 0] = np.nan
    return {
        "pseudotime": times,
        "shared_pseudotime": shared,
        "weights": weights,
        "curves": curves,
        "lineage_n_cells": lineage_n,
        "n_lineages": int(len(lineages)),
        "method": "slingshot.getCurves principal curves (Street 2018 Python equivalent)",
    }


def run_slingshot(
    embedding: np.ndarray,
    labels: np.ndarray,
    start_cluster: str,
) -> dict[str, Any]:
    lin = get_lineages(embedding, labels, start_cluster)
    cur = get_curves(embedding, labels, lin)
    return {**lin, **cur}
