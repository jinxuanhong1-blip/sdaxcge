#!/usr/bin/env python3
"""Sanity checks for SpatialFDR (no scRNA required)."""

from __future__ import annotations

import numpy as np
from milo_core import graph_spatial_fdr


def test_uniform_equal_weights() -> None:
    rng = np.random.default_rng(0)
    p = rng.uniform(0, 1, 200)
    w_dist = np.ones(200)
    q = graph_spatial_fdr(p, w_dist)
    assert np.all(q >= p - 1e-12)
    assert np.all(q <= 1 + 1e-12)
    assert np.isfinite(q).all()


def test_nan_propagation() -> None:
    p = np.array([0.01, np.nan, 0.2])
    d = np.array([1.0, 1.0, 0.0])  # 0 distance → weight inf → dropped
    q = graph_spatial_fdr(p, d)
    assert np.isnan(q[1])
    assert np.isnan(q[2])
    assert q[0] <= 1


if __name__ == "__main__":
    test_uniform_equal_weights()
    test_nan_propagation()
    print("ok")
