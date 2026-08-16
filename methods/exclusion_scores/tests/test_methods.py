from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


score = load_script("score_signatures")
tacstd2_test = load_script("test_tacstd2_exclusion")
thompson = load_script("score_thompson_emt")
spatial = load_script("score_spatial_ies")


class ScoreTests(unittest.TestCase):
    def test_weighted_mean_and_coverage(self):
        expression = {"A": [2.0], "B": [8.0]}
        value, present, expected = score.mean_score(
            [("A", 1.0), ("B", -1.0), ("C", 1.0)], expression, 0
        )
        self.assertEqual(value, -3.0)
        self.assertEqual((present, expected), (2, 3))

    def test_ips_mapping_boundaries(self):
        rows = score.read_ips(ROOT / "gene_sets" / "ips_genes.tsv")
        genes = sorted({row["gene"] for row in rows})
        expression = {gene: [float(index)] for index, gene in enumerate(genes)}
        result, present, expected = score.ips_score(expression, 0, rows)
        self.assertEqual(present, expected)
        self.assertGreaterEqual(result["IPS"], 0)
        self.assertLessEqual(result["IPS"], 10)

    def test_end_to_end_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            matrix = directory / "expression.tsv"
            output = directory / "scores.tsv"
            matrix.write_text(
                "gene\tS1\tS2\n"
                "IFNG\t1\t2\nSTAT1\t2\t3\nIDO1\t3\t4\n"
                "CXCL9\t4\t5\nCXCL10\t5\t6\nHLA-DRA\t6\t7\n"
                "B2M\t8\t9\nTAP1\t9\t10\nTAP2\t10\t11\n"
                "PDCD1\t1\t2\nAIM2\t2\t3\nCCR2\t3\t4\n",
                encoding="utf-8",
            )
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "score_signatures.py"),
                    str(matrix),
                    str(output),
                    "--signature",
                    "IFNG6_mean",
                ],
                check=True,
            )
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(float(rows[0]["IFNG6_mean"]), 3.5)
            self.assertEqual(rows[0]["IFNG6_mean_coverage"], "6/6")


class Tacstd2Tests(unittest.TestCase):
    def test_monotonic_association(self):
        marker = {f"S{i}": float(i) for i in range(1, 9)}
        scores = {f"S{i}": {"Exclusion": float(i * 2)} for i in range(1, 9)}
        result = tacstd2_test.analyze(marker, scores, "Exclusion", 999, 42)
        self.assertEqual(result["spearman_rho"], 1.0)
        self.assertGreater(result["high_minus_low"], 0)
        self.assertLess(result["permutation_p_two_sided"], 0.05)

    def test_cli_json(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            matrix = directory / "expression.tsv"
            scores = directory / "scores.tsv"
            output = directory / "result.json"
            matrix.write_text(
                "gene\tS1\tS2\tS3\tS4\tS5\nTACSTD2\t1\t2\t3\t4\t5\n",
                encoding="utf-8",
            )
            scores.write_text(
                "sample\tExclusion\nS1\t1\nS2\t2\nS3\t3\nS4\t4\nS5\t5\n",
                encoding="utf-8",
            )
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "test_tacstd2_exclusion.py"),
                    str(matrix),
                    str(scores),
                    str(output),
                    "--score",
                    "Exclusion",
                    "--permutations",
                    "99",
                ],
                check=True,
            )
            parsed = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(parsed["results"][0]["n"], 5)


class ThompsonTests(unittest.TestCase):
    def test_signed_emt_direction_and_formulas(self):
        expression = {
            gene: [1.0, 2.0, 3.0]
            for gene in set(thompson.MESENCHYMAL + thompson.INFLAMMATION)
        }
        expression.update(
            {gene: [3.0, 2.0, 1.0] for gene in thompson.EPITHELIAL}
        )
        rows, coverage = thompson.score(expression, 3)
        self.assertLess(rows[0]["Thompson_EMT"], 0)
        self.assertGreater(rows[2]["Thompson_EMT"], 0)
        self.assertAlmostEqual(
            rows[2]["Thompson_unweighted"],
            rows[2]["Thompson_inflammation"] - rows[2]["Thompson_EMT"],
        )
        self.assertEqual(coverage["Thompson_EMT_coverage"], "12/12")


class SpatialIesTests(unittest.TestCase):
    def test_signed_distance_direction(self):
        excluded, total = spatial.calculate(10.0, 40.0, 1.0)
        infiltrated, _ = spatial.calculate(40.0, 10.0, 1.0)
        balanced, _ = spatial.calculate(20.0, 20.0, 1.0)
        self.assertGreater(excluded, 0)
        self.assertLess(infiltrated, 0)
        self.assertEqual(balanced, 0)
        self.assertEqual(total, 50.0)
        self.assertAlmostEqual(excluded, -infiltrated)


if __name__ == "__main__":
    unittest.main()
