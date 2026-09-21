"""Numeric checks for the ILR basis and Dirichlet-multinomial recovery."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compositional import (  # noqa: E402
    SBP_TNK_MYE_REST,
    SBP_WITHIN_IMMUNE,
    design_matrix,
    fit_dirichlet_multinomial,
    ilr_from_sbp,
    proportions_with_pseudocount,
    sbp_orthonorm_check,
)


def test_sbp_is_orthonormal() -> None:
    assert sbp_orthonorm_check(SBP_TNK_MYE_REST)
    assert sbp_orthonorm_check(SBP_WITHIN_IMMUNE)


def test_ilr_sign() -> None:
    # More T/NK/myeloid than rest -> positive first balance.
    high = np.array([[40, 20, 20, 20]], dtype=float)
    low = np.array([[5, 5, 5, 85]], dtype=float)
    p = proportions_with_pseudocount(np.vstack([high, low]), pseudo=0.5)
    z = ilr_from_sbp(p, SBP_TNK_MYE_REST)
    assert z[0, 0] > z[1, 0]
    # More T than NK -> positive third balance.
    t_hi = proportions_with_pseudocount(np.array([[30, 5, 10, 55]], dtype=float))
    nk_hi = proportions_with_pseudocount(np.array([[5, 30, 10, 55]], dtype=float))
    assert ilr_from_sbp(t_hi, SBP_TNK_MYE_REST)[0, 2] > ilr_from_sbp(nk_hi, SBP_TNK_MYE_REST)[0, 2]


def test_dm_recovers_planted_drop() -> None:
    rng = np.random.default_rng(0)
    n = 80
    cohorts = np.array(["A"] * 40 + ["B"] * 40)
    z = np.concatenate([np.linspace(-1.5, 1.5, 40), np.linspace(-1.5, 1.5, 40)])
    x, names = design_matrix(cohorts, z, ref="A")
    # Reference = Rest. Negative coef on T and NK, near-zero on Myeloid.
    beta = np.zeros((3, x.shape[1]))
    beta[0, -1] = -0.8  # T
    beta[1, -1] = -0.7  # NK
    beta[2, -1] = -0.05  # Myeloid
    eta = x @ beta.T
    eta_full = np.concatenate([eta, np.zeros((n, 1))], axis=1)
    mu = np.exp(eta_full - eta_full.max(axis=1, keepdims=True))
    mu = mu / mu.sum(axis=1, keepdims=True)
    phi = 30.0
    counts = np.vstack([rng.multinomial(400, m) for m in mu])
    fit = fit_dirichlet_multinomial(
        counts, x, ["T", "NK", "Myeloid", "Rest"], names, n_starts=2, seed=1
    )
    assert fit.success
    assert fit.beta[0, -1] < -0.3
    assert fit.beta[1, -1] < -0.3
    assert abs(fit.beta[2, -1]) < abs(fit.beta[0, -1])


if __name__ == "__main__":
    test_sbp_is_orthonormal()
    test_ilr_sign()
    test_dm_recovers_planted_drop()
    print("ok")
