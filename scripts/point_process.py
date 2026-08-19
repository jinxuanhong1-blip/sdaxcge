"""Cross-type Ripley's K, L, and pair-correlation g(r).

Homogeneous (Lotwick–Silverman / Ripley cross-K):
    K_AB(r) = |W| / (nA nB) * sum_i sum_j 1{d_ij <= r}
    L(r)    = sqrt(K(r) / pi)
    g(r)    = (1 / (2 pi r)) dK/dr
            ≈ |W| / (nA nB 2 pi r dr) * N_pairs in [r, r+dr)

Inhomogeneous (Baddeley–Møller–Waagepetersen):
    K_inhom(r) = 1/|W| * sum_i sum_j 1{d_ij <= r} / (lambda_A(x_i) lambda_B(x_j))
    g_inhom from ring sums of the same weights.

lambda for CLDN4-high tumor cells is the kernel intensity of ALL
epithelial/tumor cells, not the CLDN4-high subset, so exclusion is not
explained merely by “tumor is dense”.

Border (reduced-sample) edge correction: only origins with distance-to-edge
>= r contribute at that r; weights are renormalised by n_keep / n_A.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def rectangle_window(xy: np.ndarray, pad: float = 0.0):
    xmin, ymin = xy.min(axis=0) - pad
    xmax, ymax = xy.max(axis=0) + pad
    area = float((xmax - xmin) * (ymax - ymin))
    return np.array([xmin, xmax, ymin, ymax], dtype=float), area


def dist_to_border(xy: np.ndarray, win: np.ndarray) -> np.ndarray:
    xmin, xmax, ymin, ymax = win
    return np.minimum.reduce(
        [xy[:, 0] - xmin, xmax - xy[:, 0], xy[:, 1] - ymin, ymax - xy[:, 1]]
    )


def intensity_at_points(
    source_xy: np.ndarray,
    query_xy: np.ndarray,
    win: np.ndarray,
    bandwidth: float,
    pix: float = 10.0,
    floor_q: float = 0.05,
) -> np.ndarray:
    """Gaussian-smoothed 2D histogram intensity (points per µm²) at query xy."""
    xmin, xmax, ymin, ymax = win
    area = max((xmax - xmin) * (ymax - ymin), 1.0)
    n_src = source_xy.shape[0]
    if n_src < 3 or query_xy.shape[0] == 0:
        return np.full(query_xy.shape[0], max(n_src, 1) / area)
    nx = max(8, int(np.ceil((xmax - xmin) / pix)))
    ny = max(8, int(np.ceil((ymax - ymin) / pix)))
    H, xe, ye = np.histogram2d(
        source_xy[:, 0],
        source_xy[:, 1],
        bins=[nx, ny],
        range=[[xmin, xmax], [ymin, ymax]],
    )
    sigma = max(bandwidth / pix, 0.5)
    dens = gaussian_filter(H.astype(float), sigma=sigma, mode="constant") / (pix * pix)
    ix = np.clip(np.searchsorted(xe, query_xy[:, 0], side="right") - 1, 0, dens.shape[0] - 1)
    iy = np.clip(np.searchsorted(ye, query_xy[:, 1], side="right") - 1, 0, dens.shape[1] - 1)
    lam = dens[ix, iy]
    pos = lam[lam > 0]
    floor = np.quantile(pos, floor_q) if pos.size else (n_src / area)
    return np.maximum(lam, floor)


def pairwise_dist(xy_a: np.ndarray, xy_b: np.ndarray) -> np.ndarray:
    """(nA x nB) Euclidean distances. Empty-safe."""
    if xy_a.size == 0 or xy_b.size == 0:
        return np.zeros((xy_a.shape[0], xy_b.shape[0]), dtype=float)
    d = xy_a[:, None, :] - xy_b[None, :, :]
    return np.sqrt((d * d).sum(axis=2))


def _pair_contrib(dist: np.ndarray, r: float, w_row: np.ndarray, w_col: np.ndarray, keep: np.ndarray) -> float:
    if not np.any(keep):
        return np.nan
    mask = dist[keep] <= r  # (n_keep, nB)
    # contrib per kept origin: w_row[i] * sum_j mask_ij * w_col[j]
    row = (mask * w_col[None, :]).sum(axis=1)
    n_a = dist.shape[0]
    n_keep = int(keep.sum())
    return float(np.dot(w_row[keep], row) * (n_a / n_keep))


def cross_k_curve(
    dist: np.ndarray,
    r: np.ndarray,
    origin_border: np.ndarray,
    w_row: np.ndarray,
    w_col: np.ndarray,
) -> np.ndarray:
    k = np.empty(r.size, dtype=float)
    for i, ri in enumerate(r):
        keep = origin_border >= ri
        k[i] = _pair_contrib(dist, float(ri), w_row, w_col, keep)
    return k


def cross_g_rings(
    dist: np.ndarray,
    r_lo: np.ndarray,
    r_hi: np.ndarray,
    origin_border: np.ndarray,
    w_row: np.ndarray,
    w_col: np.ndarray,
) -> np.ndarray:
    """Ring estimator: g(r) = [K(r_hi) - K(r_lo)] / (2 pi r_mid dr) using same weights as K."""
    r_mid = 0.5 * (r_lo + r_hi)
    dr = r_hi - r_lo
    g = np.empty(r_lo.size, dtype=float)
    for i in range(r_lo.size):
        keep = origin_border >= r_hi[i]
        k_hi = _pair_contrib(dist, float(r_hi[i]), w_row, w_col, keep)
        k_lo = _pair_contrib(dist, float(r_lo[i]), w_row, w_col, keep)
        denom = 2.0 * np.pi * r_mid[i] * dr[i]
        g[i] = (k_hi - k_lo) / denom if np.isfinite(k_hi) and np.isfinite(k_lo) else np.nan
    return g


def homogeneous_weights(n_a: int, n_b: int, area: float):
    if n_a == 0 or n_b == 0:
        return np.zeros(n_a), np.ones(max(n_b, 1))
    w_row = np.full(n_a, area / (n_a * n_b), dtype=float)
    w_col = np.ones(n_b, dtype=float)
    return w_row, w_col


def inhomogeneous_weights(lam_a: np.ndarray, lam_b: np.ndarray, area: float):
    # K = 1/|W| sum 1/(la lb)  => w_row_i * w_col_j with w_row = 1/(|W| la), w_col = 1/lb
    w_row = 1.0 / (area * np.asarray(lam_a, dtype=float))
    w_col = 1.0 / np.asarray(lam_b, dtype=float)
    return w_row, w_col


def k_to_l(k: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(k, 0.0) / np.pi)


def empirical_p(obs: float, null: np.ndarray, alternative: str = "less") -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return np.nan
    ext = np.sum(null <= obs) if alternative == "less" else np.sum(null >= obs)
    return float(1.0 + ext) / float(1.0 + null.size)
