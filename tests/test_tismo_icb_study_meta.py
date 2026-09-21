"""Locks for the TISMO ICB study-level REML analysis."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "methods" / "tismo_icb_study_meta"))

import analyze  # noqa: E402


class RemlMathTests(unittest.TestCase):
    def test_homogeneous_tau_is_zero(self):
        y = np.array([0.2, 0.2, 0.2])
        v = np.array([0.05, 0.05, 0.05])
        fit = analyze.fit_meta(y, v, method="REML")
        self.assertEqual(fit["tau2"], 0.0)
        self.assertAlmostEqual(fit["mu"], 0.2, places=10)
        self.assertAlmostEqual(fit["se"], math.sqrt(0.05 / 3), places=10)

    def test_q_and_i2_when_underdispersed(self):
        y = np.array([0.2, 0.3, 0.1])
        v = np.array([0.05, 0.05, 0.05])
        fit = analyze.fit_meta(y, v, method="REML")
        # FE mean is 0.2; Q = (0 + 0.1^2 + 0.1^2) / 0.05 = 0.4; df = 2; I^2 = 0.
        self.assertAlmostEqual(fit["mu_fe"], 0.2, places=10)
        self.assertAlmostEqual(fit["Q"], 0.4, places=10)
        self.assertEqual(fit["i2"], 0.0)
        self.assertEqual(fit["tau2"], 0.0)

    def test_modified_kh_is_at_least_as_wide_as_wald(self):
        y = np.array([0.1, 0.2, 0.0, 0.4, -0.1])
        v = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        fit = analyze.fit_meta(y, v, method="REML")
        wald = fit["ci_high"] - fit["ci_low"]
        mkh = fit["mkh_high"] - fit["mkh_low"]
        self.assertGreaterEqual(mkh + 1e-9, wald)

    def test_heterogeneous_tau_positive(self):
        y = np.array([0.0, 2.0, 0.0, 2.0])
        v = np.array([0.01, 0.01, 0.01, 0.01])
        fit = analyze.fit_meta(y, v, method="REML")
        self.assertGreater(fit["tau2"], 0.5)
        self.assertAlmostEqual(fit["mu"], 1.0, places=6)
        self.assertGreater(fit["i2"], 0.9)

    def test_shared_baseline_variance_matches_analytic(self):
        # Two slices, identical 10 baseline mice, disjoint treatment arms of 10.
        # sigma^2 = 1. Var(mean delta) = 0.15, independence formula = 0.10.
        base = frozenset(f"b{i}" for i in range(10))
        a = {"base_ids": base, "icb_ids": frozenset(f"t{i}" for i in range(10))}
        b = {"base_ids": base, "icb_ids": frozenset(f"u{i}" for i in range(10))}
        factor = analyze.cov_factor(a, b)
        self.assertAlmostEqual(factor, 0.1, places=10)
        sig = np.array([[0.2, 0.1], [0.1, 0.2]])
        ones = np.ones(2)
        v = float(ones @ sig @ ones) / 4
        self.assertAlmostEqual(v, 0.15, places=10)


class LockedReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = analyze.run()

    def test_tacstd2_49_of_64(self):
        ref = self.payload["blocks"]["Tacstd2"]["reference"]
        self.assertTrue(ref["matches_lock"])
        self.assertEqual(ref["n_up"], 49)
        self.assertEqual(ref["n_down"], 15)
        self.assertEqual(ref["n_tie"], 0)
        self.assertEqual(ref["n_slices"], 64)
        self.assertAlmostEqual(ref["wilcoxon_stat"], 439.0, places=6)
        self.assertAlmostEqual(ref["wilcoxon_p_two_sided"], 5.8398482923347564e-05, places=12)

    def test_cldn4_weaker_than_tacstd2_on_the_lock(self):
        cld = self.payload["blocks"]["Cldn4"]["reference"]
        tac = self.payload["blocks"]["Tacstd2"]["reference"]
        self.assertTrue(cld["matches_lock"])
        self.assertEqual(cld["n_up"], 34)
        self.assertEqual(cld["n_down"], 28)
        self.assertEqual(cld["n_tie"], 2)
        self.assertGreater(cld["wilcoxon_p_two_sided"], 0.05)
        self.assertLess(tac["wilcoxon_p_two_sided"], 1e-4)
        self.assertLess(cld["n_up"], tac["n_up"])

    def test_tj_lock(self):
        ref = self.payload["blocks"]["TJ"]["reference"]
        self.assertTrue(ref["matches_lock"])
        self.assertEqual(ref["n_up"], 37)
        self.assertEqual(ref["n_down"], 27)

    def test_twenty_two_studies_and_no_fake_kl(self):
        studies = self.payload["blocks"]["Tacstd2"]["studies"]
        self.assertEqual(len(studies), 22)
        lung = studies[studies["is_lung"]]
        self.assertEqual(list(lung["gse_id"]), ["GSE155972"])
        self.assertEqual(list(lung["lines"]), ["LLC"])
        audit = self.payload["audit"]
        self.assertFalse(audit["cellLineMeta_mentions_stk11"])
        self.assertFalse(audit["cellLineMeta_mentions_lkb1"])
        self.assertFalse(audit["vivoMeta_mentions_stk11"])
        self.assertFalse(audit["vivoMeta_mentions_lkb1"])
        self.assertEqual(audit["llc_cancer"], "Lung carcinoma")
        self.assertIn("Mammary", audit["kpb25l_cancer"])

    def test_loso_has_one_row_per_study(self):
        for marker in ("Tacstd2", "Cldn4", "TJ"):
            loso = self.payload["blocks"][marker]["loso"]
            self.assertEqual(len(loso), 22)
            self.assertEqual(len({row["left_out"] for row in loso}), 22)

    def test_study_sign_count_keeps_cldn4_weaker(self):
        tac = self.payload["blocks"]["Tacstd2"]["study_signs"]
        cld = self.payload["blocks"]["Cldn4"]["study_signs"]
        self.assertEqual(tac["n_up"], 18)
        self.assertEqual(cld["n_up"], 12)
        self.assertLess(tac["binomial_p_two_sided"], 0.05)
        self.assertGreater(cld["binomial_p_two_sided"], 0.05)
        self.assertLess(
            self.payload["blocks"]["Tacstd2"]["study_wilcoxon"]["p_two_sided"],
            self.payload["blocks"]["Cldn4"]["study_wilcoxon"]["p_two_sided"],
        )

    def test_primary_reml_does_not_let_a_floor_study_dominate(self):
        for marker in ("Tacstd2", "Cldn4", "TJ"):
            block = self.payload["blocks"][marker]
            self.assertLess(block["max_weight_pct"], 40.0)
            self.assertGreater(block["fit"]["mkh_high"] - block["fit"]["mkh_low"], block["fit"]["ci_high"] - block["fit"]["ci_low"] - 1e-9)
        cld = self.payload["blocks"]["Cldn4"]
        tac = self.payload["blocks"]["Tacstd2"]
        self.assertGreater(cld["technical_max_weight_pct"], 90.0)
        self.assertLess(tac["fit"]["mkh_low"], tac["fit"]["mu"])
        self.assertGreater(tac["fit"]["mkh_low"], 0.0)
        self.assertLess(cld["fit"]["mkh_low"], 0.0)
        self.assertGreater(cld["fit"]["mu"], 0.0)
        self.assertLess(cld["fit"]["mu"], tac["fit"]["mu"])

    def test_tacstd2_loso_is_sensitive_to_one_large_study(self):
        lows = [row["ci_low"] for row in self.payload["blocks"]["Tacstd2"]["loso"]]
        self.assertTrue(any(lo <= 0 for lo in lows))
        self.assertGreater(sum(lo > 0 for lo in lows), 18)
        cld_lows = [row["ci_low"] for row in self.payload["blocks"]["Cldn4"]["loso"]]
        self.assertTrue(all(lo < 0 for lo in cld_lows))


if __name__ == "__main__":
    unittest.main()
