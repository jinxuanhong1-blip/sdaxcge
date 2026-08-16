#!/usr/bin/env python3
"""Synthetic DM self-test: recover a planted TNK-down effect."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.dirichlet_multinomial import fit_dm_glm  # noqa: E402
from lib.fdr import bh_fdr  # noqa: E402


def main() -> None:
    rng = np.random.default_rng(7)
    celltypes = ["TNK", "TLS", "myeloid", "epithelial", "stromal"]
    n = 40
    x = np.array([0] * 20 + [1] * 20, dtype=float)
    # baseline pi; TACSTD2-high lowers TNK, raises epithelial
    pi0 = np.array([0.35, 0.15, 0.20, 0.15, 0.15])
    pi1 = np.array([0.18, 0.12, 0.20, 0.35, 0.15])
    counts = []
    for i in range(n):
        pi = pi1 if x[i] else pi0
        alpha = 80 * pi
        p = rng.dirichlet(alpha)
        counts.append(rng.multinomial(2000, p))
    counts = np.asarray(counts, dtype=float)
    res = fit_dm_glm(counts, x, celltypes, reference="stromal", ridge=0.2)
    q = bh_fdr(res.p_slope)
    tnk = celltypes.index("TNK")
    print("converged", res.converged, "precision", round(res.precision, 2))
    print("slopes", dict(zip(celltypes, np.round(res.slope, 3))))
    print("p", dict(zip(celltypes, np.round(res.p_slope, 4))))
    print("q", dict(zip(celltypes, np.round(q, 4))))
    if not (res.converged and res.slope[tnk] < 0 and q[tnk] < 0.10):
        raise SystemExit("self-test failed: planted TNK-down not recovered")
    print("SELFTEST_OK")


if __name__ == "__main__":
    main()
