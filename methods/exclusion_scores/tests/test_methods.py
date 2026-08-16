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


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


score = load_script("score_signatures")
tacstd2_test = load_script("test_tacstd2_exclusion")


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


if __name__ == "__main__":
    unittest.main()
