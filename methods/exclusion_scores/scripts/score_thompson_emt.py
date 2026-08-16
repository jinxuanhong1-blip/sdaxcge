#!/usr/bin/env python3
"""Reproduce the Thompson NSCLC EMT/inflammation scores.

Genes are standardized across samples within the analyzed cohort, matching the
published model. Scores are therefore cohort-dependent.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

from score_signatures import read_expression


MESENCHYMAL = ("AGER", "FN1", "MMP2", "SNAI2", "VIM", "ZEB2")
EPITHELIAL = ("CDH1", "CDH3", "CLDN4", "EPCAM", "MAL2", "ST14")
INFLAMMATION = (
    "CCL5", "CCR5", "CD274", "CD3D", "CD3E", "CD8A", "CIITA", "CTLA4",
    "CXCL10", "CXCL11", "CXCL13", "CXCL9", "GZMA", "GZMB", "HLA-DRA",
    "HLA-DRB1", "HLA-E", "IDO1", "IL2RG", "ITGAL", "LAG3", "NKG7",
    "PDCD1", "PRF1", "PTPRC", "STAT1", "TAGAP",
)


def gene_zscores(values: list[float]) -> list[float]:
    if len(values) < 2:
        raise ValueError("Thompson scores require at least two samples")
    scale = statistics.stdev(values)
    if scale == 0:
        return [0.0] * len(values)
    center = statistics.fmean(values)
    return [(value - center) / scale for value in values]


def score(
    expression: dict[str, list[float]],
    sample_count: int,
) -> tuple[list[dict[str, float]], dict[str, str]]:
    required = set(MESENCHYMAL + EPITHELIAL + INFLAMMATION)
    missing = sorted(required - set(expression))
    if missing:
        raise ValueError(
            "exact Thompson scoring requires all genes; missing: " + ", ".join(missing)
        )
    z = {gene: gene_zscores(expression[gene]) for gene in required}
    rows = []
    for index in range(sample_count):
        emt = sum(z[gene][index] for gene in MESENCHYMAL) - sum(
            z[gene][index] for gene in EPITHELIAL
        )
        inflammation = sum(z[gene][index] for gene in INFLAMMATION)
        rows.append(
            {
                "Thompson_EMT": emt,
                "Thompson_inflammation": inflammation,
                "Thompson_unweighted": inflammation - emt,
                "Thompson_weighted": -0.60 * emt + 0.19 * inflammation,
            }
        )
    return rows, {
        "Thompson_EMT_coverage": "12/12",
        "Thompson_inflammation_coverage": "27/27",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="genes x samples log-expression TSV")
    parser.add_argument("output", type=Path, help="sample-level score TSV")
    args = parser.parse_args()
    samples, expression = read_expression(args.input)
    results, coverage = score(expression, len(samples))
    fields = [
        "sample",
        "Thompson_EMT",
        "Thompson_inflammation",
        "Thompson_unweighted",
        "Thompson_weighted",
        "Thompson_EMT_coverage",
        "Thompson_inflammation_coverage",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sample, values in zip(samples, results):
            writer.writerow(
                {
                    "sample": sample,
                    **{key: f"{value:.10g}" for key, value in values.items()},
                    **coverage,
                }
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
