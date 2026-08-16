"""Visium hex-ring neighborhood scores.

Ring 1 neighbors of (row, col) on the Visium honeycomb:
  (row, col±2), (row±1, col±1)
Ring k = graph distance k (self excluded).
"""
from __future__ import annotations

from collections import defaultdict, deque
import numpy as np
import pandas as pd
from scipy import stats


def hex_adj(rows, cols):
    """Map (row,col) -> list of neighbor (row,col) that exist in the section."""
    present = set(zip(rows.astype(int), cols.astype(int)))
    adj = {}
    for r, c in present:
        cand = [(r, c - 2), (r, c + 2),
                (r - 1, c - 1), (r - 1, c + 1),
                (r + 1, c - 1), (r + 1, c + 1)]
        adj[(r, c)] = [p for p in cand if p in present]
    return adj


def rings_from_adj(adj, max_k=3):
    """For each node, return {k: list of nodes at graph distance k}."""
    out = {}
    for src in adj:
        dist = {src: 0}
        q = deque([src])
        buckets = {k: [] for k in range(1, max_k + 1)}
        while q:
            u = q.popleft()
            du = dist[u]
            if du == max_k:
                continue
            for v in adj[u]:
                if v not in dist:
                    dist[v] = du + 1
                    if dist[v] <= max_k:
                        buckets[dist[v]].append(v)
                        q.append(v)
        out[src] = buckets
    return out


def neighbor_means(coord_keys, values, ring_map, k, min_n=3):
    """Mean of `values` over ring-k neighbors. coord_keys aligned with values."""
    key2i = {k_: i for i, k_ in enumerate(coord_keys)}
    out = np.full(len(values), np.nan)
    for i, src in enumerate(coord_keys):
        nbrs = ring_map[src][k]
        if len(nbrs) < min_n:
            continue
        idx = [key2i[n] for n in nbrs]
        out[i] = np.nanmean(values[idx])
    return out


def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 20 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(a[m], b[m])
    return float(rho), float(p), int(m.sum())


def partial_spearman(x, y, z):
    """Spearman of residuals after ranking and regressing out z."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if m.sum() < 20:
        return np.nan, np.nan, int(m.sum())
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])
    # residualize
    def resid(a, b):
        b1 = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ coef
    rxr = resid(rx, rz)
    ryr = resid(ry, rz)
    if np.std(rxr) == 0 or np.std(ryr) == 0:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(rxr, ryr)
    return float(rho), float(p), int(m.sum())


def neighborhood_rows(sample, group, dataset, rows, cols, cldn4, epi, t, b, tb,
                      max_k=3):
    """Return list of dicts, one per ring, for this section."""
    keys = list(zip(rows.astype(int), cols.astype(int)))
    adj = hex_adj(rows, cols)
    rmap = rings_from_adj(adj, max_k=max_k)
    out = []
    for k in range(1, max_k + 1):
        t_nb = neighbor_means(keys, t, rmap, k)
        b_nb = neighbor_means(keys, b, rmap, k)
        tb_nb = neighbor_means(keys, tb, rmap, k)
        row = dict(sample=sample, group=group, dataset=dataset, ring=k,
                   n_spots=len(cldn4))
        for lab, idx, nb in [
            ("CLDN4_vs_Tnb", cldn4, t_nb),
            ("CLDN4_vs_Bnb", cldn4, b_nb),
            ("CLDN4_vs_TBnb", cldn4, tb_nb),
            ("epi_vs_TBnb", epi, tb_nb),
        ]:
            rho, p, n = spearman(idx, nb)
            row[f"rho_{lab}"] = rho
            row[f"p_{lab}"] = p
            row[f"n_{lab}"] = n
        # partial: CLDN4 vs neighbor TB, controlling for own epithelial score
        # (removes "this is just a tumor-core vs stroma" composition leak)
        pr, pp, pn = partial_spearman(cldn4, tb_nb, epi)
        row["rho_CLDN4_vs_TBnb_partialEpi"] = pr
        row["p_CLDN4_vs_TBnb_partialEpi"] = pp
        row["n_partial"] = pn
        out.append(row)
    return out
