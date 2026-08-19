#!/usr/bin/env python3
"""CSR sanity check: two independent uniform patterns should have K ≈ πr² and g ≈ 1."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from point_process import (
    cross_g_rings,
    cross_k_curve,
    dist_to_border,
    homogeneous_weights,
    k_to_l,
    pairwise_dist,
    rectangle_window,
)


def main():
    rng = np.random.default_rng(1)
    win_xy = np.array([[0.0, 0.0], [800.0, 600.0]])
    win, area = rectangle_window(win_xy)
    r = np.arange(20.0, 121.0, 10.0)
    r_lo, r_hi = r[:-1], r[1:]
    r_mid = 0.5 * (r_lo + r_hi)
    ratios = []
    gs = []
    for _ in range(8):
        a = rng.uniform([0, 0], [800, 600], size=(250, 2))
        b = rng.uniform([0, 0], [800, 600], size=(180, 2))
        dist = pairwise_dist(a, b)
        border = dist_to_border(a, win)
        w_row, w_col = homogeneous_weights(len(a), len(b), area)
        k = cross_k_curve(dist, r, border, w_row, w_col)
        g = cross_g_rings(dist, r_lo, r_hi, border, w_row, w_col)
        theo = np.pi * r**2
        ratios.append(np.nanmean(k / theo))
        gs.append(np.nanmean(g))
        _ = k_to_l(k)
    print("mean K / (pi r^2) over sims:", float(np.mean(ratios)), "sd", float(np.std(ratios)))
    print("mean g(r) over sims:", float(np.mean(gs)), "sd", float(np.std(gs)))
    if not (0.85 < np.mean(ratios) < 1.15):
        raise SystemExit("CSR K calibration failed")
    if not (0.75 < np.mean(gs) < 1.25):
        raise SystemExit("CSR g calibration failed")
    print("CSR calibration OK")


if __name__ == "__main__":
    main()
