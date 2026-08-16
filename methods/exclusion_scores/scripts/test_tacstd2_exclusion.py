#!/usr/bin/env python3
"""Test TACSTD2 association with exclusion-related scores.

Reports Spearman rho with a deterministic two-sided permutation P value and a
median-split group contrast. This is an association analysis, not a survival or
treatment-predictive model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float:
    left_center = statistics.fmean(left)
    right_center = statistics.fmean(right)
    numerator = sum(
        (x - left_center) * (y - right_center) for x, y in zip(left, right)
    )
    denominator = math.sqrt(
        sum((x - left_center) ** 2 for x in left)
        * sum((y - right_center) ** 2 for y in right)
    )
    return numerator / denominator if denominator else math.nan


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(average_ranks(left), average_ranks(right))


def permutation_p(
    left: list[float],
    right: list[float],
    observed: float,
    permutations: int,
    rng: random.Random,
) -> float:
    extreme = 0
    shuffled = right.copy()
    for _ in range(permutations):
        rng.shuffle(shuffled)
        candidate = spearman(left, shuffled)
        extreme += abs(candidate) >= abs(observed)
    return (extreme + 1) / (permutations + 1)


def read_gene(path: Path, gene: str) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter="\t")
        try:
            header = next(rows)
        except StopIteration:
            raise ValueError("expression matrix is empty")
        for row in rows:
            if row and row[0].strip().upper() == gene:
                return {sample: float(value) for sample, value in zip(header[1:], row[1:])}
    raise ValueError(f"{gene} was not found in expression matrix")


def read_scores(path: Path, names: list[str]) -> dict[str, dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not rows.fieldnames or "sample" not in rows.fieldnames:
            raise ValueError("score table requires a 'sample' column")
        missing = set(names) - set(rows.fieldnames)
        if missing:
            raise ValueError(f"score columns not found: {', '.join(sorted(missing))}")
        result = {}
        for row in rows:
            parsed = {}
            for name in names:
                if row[name] not in ("", "NA", "NaN"):
                    parsed[name] = float(row[name])
            result[row["sample"]] = parsed
        return result


def analyze(
    tacstd2: dict[str, float],
    scores: dict[str, dict[str, float]],
    score_name: str,
    permutations: int,
    seed: int,
) -> dict[str, float | int | str]:
    samples = [
        sample
        for sample in tacstd2
        if sample in scores and score_name in scores[sample]
    ]
    if len(samples) < 5:
        raise ValueError(f"{score_name}: fewer than five matched samples")
    x = [tacstd2[sample] for sample in samples]
    y = [scores[sample][score_name] for sample in samples]
    rho = spearman(x, y)
    if math.isnan(rho):
        raise ValueError(f"{score_name}: correlation undefined (constant values)")
    p_value = permutation_p(x, y, rho, permutations, random.Random(seed))
    median = statistics.median(x)
    low = [value for value, marker in zip(y, x) if marker <= median]
    high = [value for value, marker in zip(y, x) if marker > median]
    if not high:
        raise ValueError(f"{score_name}: median split produced an empty high group")
    return {
        "score": score_name,
        "n": len(samples),
        "spearman_rho": rho,
        "permutation_p_two_sided": p_value,
        "permutations": permutations,
        "tacstd2_median": median,
        "n_low": len(low),
        "n_high": len(high),
        "score_mean_low": statistics.fmean(low),
        "score_mean_high": statistics.fmean(high),
        "high_minus_low": statistics.fmean(high) - statistics.fmean(low),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expression", type=Path, help="genes x samples TSV")
    parser.add_argument("scores", type=Path, help="sample x score TSV")
    parser.add_argument("output", type=Path, help="JSON output")
    parser.add_argument(
        "--score",
        action="append",
        required=True,
        help="score column to test (repeatable)",
    )
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=613)
    args = parser.parse_args()
    if args.permutations < 99:
        raise ValueError("--permutations must be >=99")
    tacstd2 = read_gene(args.expression, "TACSTD2")
    scores = read_scores(args.scores, args.score)
    results = [
        analyze(tacstd2, scores, name, args.permutations, args.seed + index)
        for index, name in enumerate(args.score)
    ]
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "marker": "TACSTD2",
                "multiplicity_note": "Adjust across prespecified score tests (e.g. BH FDR).",
                "results": results,
            },
            handle,
            indent=2,
        )
        handle.write("\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
