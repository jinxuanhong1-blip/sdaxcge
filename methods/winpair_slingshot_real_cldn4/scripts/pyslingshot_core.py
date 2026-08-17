"""Vendored core of pyslingshot-bio v0.1.0 (MIT; omicverse/py-Slingshot).

Python port of Slingshot (Street et al. BMC Genomics 2018). Vendored so the
analysis does not depend on the optional ggplot2_py plotting extra that the
upstream package __init__ imports.

References:
- Street, K. et al. *Slingshot: cell lineage and pseudotime inference for
  single-cell transcriptomics.* BMC Genomics 19, 477 (2018).
- https://github.com/omicverse/py-Slingshot
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist


# ----- containers ---------------------------------------------------------- #


@dataclass
class SlingshotResult:
    """All Slingshot outputs in one object.

    Attributes:
        reduced_dim: (n × d) input coordinates.
        cluster_labels: (n,) cluster assignment for each cell.
        cluster_centers: dict cluster_label → (d,) centroid.
        mst_adj: (k × k) adjacency matrix of the MST between clusters.
        cluster_order: (k,) ordering used for the MST adjacency.
        lineages: list of cluster-label sequences (each a path from root to leaf).
        weights: (n × n_lineages) per-cell weight per lineage.
        pseudotime: (n × n_lineages) per-cell pseudotime per lineage.
        curves: list of (m × d) principal-curve point matrices.
    """

    reduced_dim: np.ndarray
    cluster_labels: np.ndarray
    cluster_centers: dict
    mst_adj: np.ndarray
    cluster_order: list
    lineages: list[list]
    weights: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    pseudotime: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    curves: list[np.ndarray] = field(default_factory=list)


# ----- getLineages --------------------------------------------------------- #


def _compute_centers(X: np.ndarray, labels: np.ndarray) -> dict:
    return {c: X[labels == c].mean(axis=0) for c in np.unique(labels)}


def _build_mst(centers: dict) -> tuple[np.ndarray, list]:
    """MST on centroid pairwise distances. Returns symmetric (k×k) adjacency."""
    order = list(centers.keys())
    C = np.stack([centers[c] for c in order])
    D = cdist(C, C, "euclidean")
    mst = minimum_spanning_tree(D).toarray()
    mst = np.maximum(mst, mst.T)  # symmetrize
    return mst, order


def _enumerate_root_to_leaf_paths(mst_adj: np.ndarray, root_idx: int) -> list[list[int]]:
    """DFS from root: every distinct leaf path becomes a lineage."""
    n = mst_adj.shape[0]
    paths: list[list[int]] = []

    def dfs(node: int, parent: int, path: list[int]):
        # Find children
        children = [j for j in range(n)
                    if mst_adj[node, j] > 0 and j != parent]
        if not children:
            paths.append(path + [node])
            return
        for ch in children:
            dfs(ch, node, path + [node])

    dfs(root_idx, -1, [])
    return paths


def _pick_root(mst_adj: np.ndarray, order: list, start_cluster=None) -> int:
    """Pick root cluster — user-given if possible, else max-average-path-length leaf."""
    if start_cluster is not None and start_cluster in order:
        return order.index(start_cluster)

    n = mst_adj.shape[0]
    if n == 1:
        return 0

    # Leaves only (degree 1 nodes)
    degrees = (mst_adj > 0).sum(axis=1)
    leaves = [i for i in range(n) if degrees[i] == 1]
    if not leaves:
        return 0

    # For each leaf as candidate root, compute mean shortest-path length to other leaves
    from scipy.sparse.csgraph import shortest_path
    sp = shortest_path((mst_adj > 0).astype(float))
    best_leaf = leaves[0]
    best_avg = -1.0
    for leaf in leaves:
        avg = sp[leaf, [l for l in leaves if l != leaf]].mean()
        if avg > best_avg:
            best_avg = avg
            best_leaf = leaf
    return best_leaf


def getLineages(
    data: np.ndarray | pd.DataFrame,
    cluster_labels: np.ndarray | list,
    *,
    start_cluster=None,
) -> SlingshotResult:
    """1:1 port of slingshot::getLineages.

    Args:
        data: (n × d) reduced-dim coords.
        cluster_labels: (n,) cluster IDs.
        start_cluster: optional root cluster ID.

    Returns:
        SlingshotResult with .lineages populated (no curves yet).
    """
    X = np.asarray(data, dtype=np.float64)
    labels = np.asarray(cluster_labels)
    if labels.dtype.kind in {"i", "u"}:
        labels = labels.astype(int)
    centers = _compute_centers(X, labels)
    mst_adj, order = _build_mst(centers)
    root_idx = _pick_root(mst_adj, order, start_cluster=start_cluster)
    leaf_paths_idx = _enumerate_root_to_leaf_paths(mst_adj, root_idx)
    lineages = [[order[i] for i in path] for path in leaf_paths_idx]
    # Sort by length descending
    lineages.sort(key=len, reverse=True)
    return SlingshotResult(
        reduced_dim=X,
        cluster_labels=labels,
        cluster_centers=centers,
        mst_adj=mst_adj,
        cluster_order=order,
        lineages=lineages,
    )


# ----- principal curve (Hastie-Stuetzle, per lineage) --------------------- #


def _principal_curve(
    X: np.ndarray,
    init: np.ndarray,
    *,
    weights: np.ndarray | None = None,
    max_iter: int = 10,
    smoother_span: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Hastie-Stuetzle principal curve via lowess smoothing.

    Args:
        X: (n × d) cells.
        init: (m × d) initial curve points.
        weights: (n,) per-cell weight (1 = full, 0 = ignored).

    Returns:
        s: (m × d) refined curve.
        lambda_proj: (n,) per-cell arc-length along the curve.
        ord_idx: (m,) reordering of init by arc-length (always sorted ascending).
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    s = init.copy()
    if weights is None:
        weights = np.ones(X.shape[0])
    fit_idx = np.flatnonzero(weights > 0)
    if fit_idx.size < 10:
        fit_idx = np.arange(X.shape[0])
    X_fit = X[fit_idx]

    for _ in range(max_iter):
        # Project each contributing cell onto current curve segment
        D2 = cdist(X_fit, s, "sqeuclidean")
        nearest = np.argmin(D2, axis=1)
        # Approximate arc-length along the curve up to each segment endpoint
        seg_lens = np.sqrt(((np.diff(s, axis=0)) ** 2).sum(axis=1))
        cum_arc = np.concatenate([[0.0], np.cumsum(seg_lens)])

        # Per-cell projection parameter — use segment indices' arc length
        lam = cum_arc[nearest]
        # Refine via a per-dim lowess smoother of X vs lam
        order = np.argsort(lam)
        lam_sorted = lam[order]
        X_sorted = X_fit[order]
        w_sorted = weights[fit_idx][order]
        # Smooth each dim. Cap the curve grid — n_cells points is O(n²) and hangs
        # on 10k+ cell objects (Street 2018 uses a much smaller smoother grid).
        n_grid = int(min(80, max(20, X_fit.shape[0] // 50)))
        new_s = np.zeros((n_grid, X.shape[1]))
        new_lam_grid = np.linspace(lam_sorted.min(), lam_sorted.max(), new_s.shape[0])
        for d in range(X.shape[1]):
            sm = lowess(
                X_sorted[:, d],
                lam_sorted,
                frac=smoother_span,
                it=0,
                return_sorted=True,
            )
            # Interpolate onto new_lam_grid
            new_s[:, d] = np.interp(new_lam_grid, sm[:, 0], sm[:, 1])
        # Reorder curve to ensure monotone arc-length
        new_seg = np.sqrt(((np.diff(new_s, axis=0)) ** 2).sum(axis=1))
        if (new_seg < 1e-12).all():
            break
        s = new_s

    # Final projection
    D2 = cdist(X, s, "sqeuclidean")
    nearest = np.argmin(D2, axis=1)
    seg_lens = np.sqrt(((np.diff(s, axis=0)) ** 2).sum(axis=1))
    cum_arc = np.concatenate([[0.0], np.cumsum(seg_lens)])
    lam = cum_arc[nearest]
    return s, lam, np.arange(s.shape[0])


# ----- getCurves ----------------------------------------------------------- #


def getCurves(
    lineages_result: SlingshotResult,
    *,
    max_iter: int = 10,
    smoother_span: float = 0.5,
) -> SlingshotResult:
    """1:1 port of slingshot::getCurves (simplified — fits each lineage independently).

    Args:
        lineages_result: output of ``getLineages``.
        max_iter: principal-curve refinement iterations.
        smoother_span: LOWESS span (0–1).

    Returns:
        the same result with ``weights``, ``pseudotime``, ``curves`` populated.
    """
    sr = lineages_result
    X = sr.reduced_dim
    centers = sr.cluster_centers
    lineages = sr.lineages
    n = X.shape[0]
    L = len(lineages)

    weights = np.zeros((n, L))
    pseudotime = np.full((n, L), np.nan)
    curves: list[np.ndarray] = []

    for li, lineage in enumerate(lineages):
        # Cells in this lineage = those whose cluster is in the path
        in_lineage = np.isin(sr.cluster_labels, lineage)
        weights[:, li] = in_lineage.astype(float)
        # Initial curve = sequence of cluster centroids
        init = np.stack([centers[c] for c in lineage])
        # Cells contributing weight
        w = weights[:, li]
        curve, lam, _ = _principal_curve(
            X, init, weights=w, max_iter=max_iter,
            smoother_span=smoother_span,
        )
        curves.append(curve)
        pseudotime[in_lineage, li] = lam[in_lineage]

    sr.weights = weights
    sr.pseudotime = pseudotime
    sr.curves = curves
    return sr


# ----- one-shot ------------------------------------------------------------ #


def slingshot(
    data: np.ndarray | pd.DataFrame,
    cluster_labels: np.ndarray | list,
    *,
    start_cluster=None,
    max_iter: int = 10,
    smoother_span: float = 0.5,
) -> SlingshotResult:
    """End-to-end: getLineages → getCurves."""
    sr = getLineages(data, cluster_labels, start_cluster=start_cluster)
    return getCurves(sr, max_iter=max_iter, smoother_span=smoother_span)
