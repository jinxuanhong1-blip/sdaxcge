"""Unit checks for partial Spearman and BH-FDR. No TCGA download."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats import bh_fdr, partial_pearson, partial_spearman, random_effects_meta


class StatsTests(unittest.TestCase):
    def test_marginal_matches_spearman(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=200)
        y = 0.5 * x + rng.normal(size=200)
        rho, p, n = partial_spearman(x, y)
        from scipy import stats

        r, sp = stats.spearmanr(x, y)
        self.assertEqual(n, 200)
        self.assertAlmostEqual(rho, r, places=10)
        self.assertAlmostEqual(p, sp, places=6)

    def test_shared_covariate_is_removed(self):
        rng = np.random.default_rng(1)
        z = rng.normal(size=400)
        x = z + 0.01 * rng.normal(size=400)
        y = z + 0.01 * rng.normal(size=400)
        rho_m, _, _ = partial_spearman(x, y)
        rho_p, p, n = partial_spearman(x, y, [z])
        self.assertGreater(rho_m, 0.9)
        self.assertLess(abs(rho_p), 0.2)
        self.assertEqual(n, 400)
        self.assertGreater(p, 0.01)

    def test_residual_association_survives_nuisance(self):
        rng = np.random.default_rng(2)
        z = rng.normal(size=500)
        x = rng.normal(size=500)
        y = 0.8 * x + 0.5 * z + rng.normal(size=500)
        rho, p, _ = partial_spearman(x, y, [z])
        self.assertGreater(rho, 0.5)
        self.assertLess(p, 1e-6)

    def test_bh_monotonic_and_bounds(self):
        q = bh_fdr([0.001, 0.04, 0.2, np.nan])
        self.assertTrue(np.isnan(q[3]))
        self.assertLessEqual(q[0], q[1] + 1e-15)
        self.assertTrue(np.all(q[:3] <= 1))
        self.assertGreater(q[0], 0)

    def test_partial_pearson_removes_linear_covariate(self):
        rng = np.random.default_rng(3)
        z = rng.normal(size=500)
        x = rng.normal(size=500)
        y = 0.7 * x + 1.2 * z + rng.normal(scale=0.4, size=500)
        rho, p, n = partial_pearson(x, y, [z])
        self.assertGreater(rho, 0.6)
        self.assertLess(p, 1e-8)
        self.assertEqual(n, 500)
        rho_s, _, _ = partial_spearman(x, y, [z])
        self.assertGreater(rho_s, 0.5)

    def test_meta_positive_when_all_positive(self):
        meta = random_effects_meta([0.4, 0.35, 0.5, 0.3], [100, 120, 80, 90], k=0)
        self.assertGreater(meta["rho"], 0.3)
        self.assertLess(meta["p"], 1e-6)
        self.assertEqual(meta["n_cohorts"], 4)


if __name__ == "__main__":
    unittest.main()
