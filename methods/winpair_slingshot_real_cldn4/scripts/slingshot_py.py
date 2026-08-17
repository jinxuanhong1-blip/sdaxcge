"""Documented Python equivalent of Slingshot (Street et al. BMC Genomics 2018).

Primary path: ``pyslingshot-bio`` (module ``pyslingshot``), a NumPy/SciPy port of
R Slingshot with reported Spearman ≥0.99 vs R on a Y-shaped fixture
(https://github.com/omicverse/py-Slingshot; Street et al. 2018).

Fallback (same paper, same steps): MST over cluster centroids → lineages from a
named start cluster to MST leaves → per-lineage principal curve (Hastie–Stuetzle
via LOWESS) → arc-length pseudotime. This is used only if ``pyslingshot-bio``
cannot be imported or returns no lineages. It is **not** diffusion pseudotime.

R Slingshot / tradeSeq are preferred when ``Rscript`` + Bioconductor ``slingshot``
are installed; this module does not pretend they ran.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from statsmodels.nonparametric.smoothers_lowess import lowess


@dataclass
class SlingshotResult:
    engine: str
    start_cluster: str
    lineages: list[list[str]]
    pseudotime: np.ndarray  # (n_cells, n_lineages) NaN if unassigned
    lineage_names: list[str]
    cluster_centroids: dict[str, np.ndarray]
    curves_embed: dict[str, np.ndarray]  # lineage -> (n_points, 2) for plotting
    notes: list[str] = field(default_factory=list)

    @property
    def shared_pseudotime(self) -> np.ndarray:
        """Mean of finite lineage times (Street shared-pseudotime analogue)."""
        pt = self.pseudotime.astype(float)
        n = np.isfinite(pt).sum(axis=1)
        s = np.nansum(pt, axis=1)
        out = np.full(pt.shape[0], np.nan)
        ok = n > 0
        out[ok] = s[ok] / n[ok]
        return out


def r_slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "available": False,
            "reason": "Rscript not on PATH",
            "used": "pyslingshot-bio or local Street-2018 MST+principal-curve port",
        }
    try:
        proc = subprocess.run(
            [rscript, "-e", 'packageVersion("slingshot")'],
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "available": False,
            "reason": f"Rscript slingshot probe failed: {exc}",
            "used": "pyslingshot-bio or local Street-2018 port",
        }
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:300],
            "used": "pyslingshot-bio or local Street-2018 port",
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


def _centroids(X: np.ndarray, labels: np.ndarray) -> dict[str, np.ndarray]:
    out = {}
    for lab in np.unique(labels):
        mask = labels == lab
        if mask.sum() == 0:
            continue
        out[str(lab)] = X[mask].mean(axis=0)
    return out


def _mst_lineages(centroids: dict[str, np.ndarray], start: str) -> list[list[str]]:
    labs = list(centroids.keys())
    if start not in labs:
        raise ValueError(f"start cluster {start} not in {labs}")
    coords = np.vstack([centroids[k] for k in labs])
    d = cdist(coords, coords, metric="euclidean")
    mst = minimum_spanning_tree(d).toarray()
    adj: dict[str, list[str]] = {k: [] for k in labs}
    for i, a in enumerate(labs):
        for j, b in enumerate(labs):
            if mst[i, j] > 0 or mst[j, i] > 0:
                adj[a].append(b)
    # DFS paths from start to leaves
    lineages: list[list[str]] = []

    def walk(node: str, parent: str | None, path: list[str]) -> None:
        kids = [c for c in adj[node] if c != parent]
        if not kids:
            lineages.append(path[:])
            return
        for c in kids:
            walk(c, node, path + [c])

    walk(start, None, [start])
    if not lineages:
        lineages = [[start]]
    return lineages


def _project_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, np.ndarray, float]:
    ab = b - a
    denom = float(np.dot(ab, ab)) + 1e-12
    t = float(np.clip(np.dot(p - a, ab) / denom, 0.0, 1.0))
    proj = a + t * ab
    dist = float(np.linalg.norm(p - proj))
    return t, proj, dist


def _polyline_project(X: np.ndarray, vertices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Project points onto a polyline; return arc-length and squared distance."""
    segs = vertices[1:] - vertices[:-1]
    seglen = np.linalg.norm(segs, axis=1)
    seglen = np.maximum(seglen, 1e-12)
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    n = X.shape[0]
    best_s = np.zeros(n)
    best_d = np.full(n, np.inf)
    for i in range(len(seglen)):
        a, b = vertices[i], vertices[i + 1]
        t, proj, dist = (
            np.empty(n),
            np.empty_like(X),
            np.empty(n),
        )
        ab = b - a
        denom = float(np.dot(ab, ab)) + 1e-12
        tt = ((X - a) @ ab) / denom
        tt = np.clip(tt, 0.0, 1.0)
        proj = a + tt[:, None] * ab
        dist = np.linalg.norm(X - proj, axis=1)
        s = cum[i] + tt * seglen[i]
        better = dist < best_d
        best_d[better] = dist[better]
        best_s[better] = s[better]
    return best_s, best_d


def _principal_curve(X: np.ndarray, init_vertices: np.ndarray, n_iter: int = 8, frac: float = 0.25) -> np.ndarray:
    """Hastie–Stuetzle principal curve, initialized from the centroid polyline."""
    verts = init_vertices.copy()
    if X.shape[0] < 20:
        return verts
    for _ in range(n_iter):
        s, _ = _polyline_project(X, verts)
        order = np.argsort(s)
        s_ord = s[order]
        # unique-ish grid
        n_grid = min(60, max(8, X.shape[0] // 40))
        grid = np.linspace(s_ord[0], s_ord[-1], n_grid)
        new_verts = []
        for d in range(X.shape[1]):
            y = lowess(X[order, d], s_ord, frac=frac, return_sorted=True, it=1)
            # y is (n, 2) sorted by x
            new_verts.append(np.interp(grid, y[:, 0], y[:, 1]))
        verts = np.column_stack(new_verts)
        # drop near-duplicates
        keep = [0]
        for i in range(1, len(verts)):
            if np.linalg.norm(verts[i] - verts[keep[-1]]) > 1e-6:
                keep.append(i)
        verts = verts[keep]
        if len(verts) < 2:
            return init_vertices
    return verts


def _embed_curve(curve: np.ndarray, embed: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Map a high-D curve onto a 2D embedding by nearest-cell smoothing."""
    # nearest cells to sampled curve points
    n_pts = min(80, max(10, curve.shape[0]))
    idx = np.linspace(0, curve.shape[0] - 1, n_pts).astype(int)
    pts = curve[idx]
    d = cdist(pts, X)
    nn = d.argmin(axis=1)
    return embed[nn]


def run_local_street(
    X: np.ndarray,
    labels: np.ndarray,
    start_cluster: str,
    embed2d: np.ndarray | None = None,
) -> SlingshotResult:
    labels = np.asarray(labels).astype(str)
    cents = _centroids(X, labels)
    lineages = _mst_lineages(cents, start_cluster)
    n = X.shape[0]
    pt = np.full((n, len(lineages)), np.nan)
    curves_embed: dict[str, np.ndarray] = {}
    names = []
    for li, lin in enumerate(lineages):
        name = "Lineage" + str(li + 1)
        names.append(name)
        verts = np.vstack([cents[c] for c in lin])
        # cells in these clusters
        mask = np.isin(labels, lin)
        if mask.sum() < 5:
            continue
        curve = _principal_curve(X[mask], verts)
        s, dist = _polyline_project(X, curve)
        # assign cells closer to this lineage than a loose threshold:
        # keep cells whose nearest lineage cluster is on this path
        pt[mask, li] = s[mask]
        # normalize 0-1 within lineage
        finite = np.isfinite(pt[:, li])
        if finite.any():
            lo, hi = np.nanmin(pt[finite, li]), np.nanmax(pt[finite, li])
            if hi > lo:
                pt[finite, li] = (pt[finite, li] - lo) / (hi - lo)
        if embed2d is not None:
            curves_embed[name] = _embed_curve(curve, embed2d, X)
        else:
            curves_embed[name] = verts[:, :2] if verts.shape[1] >= 2 else verts
    return SlingshotResult(
        engine="local_street2018_mst_principal_curve",
        start_cluster=start_cluster,
        lineages=lineages,
        pseudotime=pt,
        lineage_names=names,
        cluster_centroids=cents,
        curves_embed=curves_embed,
        notes=[
            "Street et al. 2018 steps: MST centroids, lineages from start to leaves, LOWESS principal curve.",
            "Not diffusion pseudotime. Independent per-lineage curves (no shared-stem shrinkage).",
        ],
    )


def run_pyslingshot_bio(
    X: np.ndarray,
    labels: np.ndarray,
    start_cluster: str,
    embed2d: np.ndarray | None = None,
) -> SlingshotResult:
    # Vendored pyslingshot-bio core (Street 2018 port). Package __init__ needs ggplot2_py.
    from pyslingshot_core import slingshot as _slingshot

    labels = np.asarray(labels).astype(str)
    sr = _slingshot(X, labels, start_cluster=start_cluster, max_iter=5, smoother_span=0.4)
    pt = np.asarray(sr.pseudotime, dtype=float)
    if pt.ndim == 1:
        pt = pt[:, None]
    # lineages
    lineages = []
    if hasattr(sr, "lineages"):
        raw = sr.lineages
        if isinstance(raw, dict):
            for v in raw.values():
                lineages.append([str(x) for x in v])
        else:
            for v in raw:
                lineages.append([str(x) for x in v])
    names = [f"Lineage{i+1}" for i in range(pt.shape[1])]
    cents = _centroids(X, labels)
    curves_embed: dict[str, np.ndarray] = {}
    if embed2d is not None:
        for i, name in enumerate(names):
            col = pt[:, i]
            finite = np.isfinite(col)
            if finite.sum() < 10:
                continue
            order = np.argsort(col[finite])
            coords = embed2d[finite][order]
            # downsample + smooth
            n_pts = min(80, max(15, finite.sum() // 40))
            idx = np.linspace(0, len(coords) - 1, n_pts).astype(int)
            curves_embed[name] = coords[idx]
    return SlingshotResult(
        engine="pyslingshot-bio (Street et al. 2018 Python port)",
        start_cluster=start_cluster,
        lineages=lineages or [["unknown"]],
        pseudotime=pt,
        lineage_names=names,
        cluster_centroids=cents,
        curves_embed=curves_embed,
        notes=[
            "Installed pyslingshot-bio; module pyslingshot. Documented R-parity port of Slingshot.",
            "https://github.com/omicverse/py-Slingshot",
        ],
    )


def run_slingshot(
    X: np.ndarray,
    labels: np.ndarray,
    start_cluster: str,
    embed2d: np.ndarray | None = None,
) -> SlingshotResult:
    """Run a real Slingshot lineage fit. Never silently falls back to DPT."""
    notes = []
    try:
        res = run_pyslingshot_bio(X, labels, start_cluster, embed2d=embed2d)
        if res.pseudotime.size and np.isfinite(res.pseudotime).any():
            return res
        notes.append("pyslingshot-bio returned empty pseudotime; using local Street port")
    except Exception as exc:  # noqa: BLE001 — documented fallback
        notes.append(f"pyslingshot-bio failed ({type(exc).__name__}: {exc}); using local Street port")
    res = run_local_street(X, labels, start_cluster, embed2d=embed2d)
    res.notes = notes + res.notes
    return res
