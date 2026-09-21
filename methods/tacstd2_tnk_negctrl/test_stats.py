"""Checks for partial Spearman and the attenuation contrast. No cohort data."""

import numpy as np

from stats import dl_meta, partial_spearman, spearman


def test_partial_removes_a_shared_covariate():
    rng = np.random.default_rng(0)
    n = 400
    z = rng.normal(size=n)
    x = z + 0.15 * rng.normal(size=n)
    y = -z + 0.15 * rng.normal(size=n)
    raw = spearman(x, y)
    adj = partial_spearman(x, y, [z])
    assert raw < -0.7
    assert abs(adj) < 0.2


def test_irrelevant_covariate_does_not_attenuate():
    rng = np.random.default_rng(1)
    n = 400
    x = rng.normal(size=n)
    y = -0.8 * x + rng.normal(size=n)
    z = rng.normal(size=n)
    raw = spearman(x, y)
    adj = partial_spearman(x, y, [z])
    assert raw < -0.5
    assert abs((adj - raw) - 0) < 0.08


def test_dl_recovers_a_common_rho():
    meta = dl_meta([-0.4, -0.42, -0.38, -0.41], [80, 90, 70, 100], k_cov=0)
    assert abs(meta["rho"] + 0.40) < 0.03
    assert meta["k"] == 4
    assert meta["I2"] < 0.05
    assert meta["p"] < 1e-6


if __name__ == "__main__":
    test_partial_removes_a_shared_covariate()
    test_irrelevant_covariate_does_not_attenuate()
    test_dl_recovers_a_common_rho()
    print("ok")
