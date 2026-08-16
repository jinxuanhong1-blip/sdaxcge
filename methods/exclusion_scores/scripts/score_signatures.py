#!/usr/bin/env python3
"""Score transparent immune/exclusion signatures and canonical IPS.

Input is a tab-delimited, genes-by-samples matrix. The first column contains
HGNC symbols; remaining columns are log-scale expression values.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SETS = ROOT / "gene_sets" / "signatures.tsv"
DEFAULT_IPS = ROOT / "gene_sets" / "ips_genes.tsv"


def read_expression(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter="\t")
        try:
            header = next(rows)
        except StopIteration:
            raise ValueError("expression matrix is empty")
        if len(header) < 2:
            raise ValueError("expression matrix needs a gene column and >=1 sample")
        samples = header[1:]
        if len(samples) != len(set(samples)):
            raise ValueError("sample names must be unique")
        expression: dict[str, list[float]] = {}
        for line_no, row in enumerate(rows, 2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"line {line_no}: expected {len(header)} columns")
            gene = row[0].strip().upper()
            if not gene:
                raise ValueError(f"line {line_no}: blank gene symbol")
            if gene in expression:
                raise ValueError(f"line {line_no}: duplicate gene {gene}")
            try:
                values = [float(value) for value in row[1:]]
            except ValueError as exc:
                raise ValueError(f"line {line_no}: non-numeric expression") from exc
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"line {line_no}: non-finite expression")
            expression[gene] = values
    if not expression:
        raise ValueError("expression matrix has no gene rows")
    return samples, expression


def read_signatures(path: Path) -> dict[str, list[tuple[str, float]]]:
    signatures: dict[str, list[tuple[str, float]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"signature", "gene", "weight"}
        if not rows.fieldnames or not required.issubset(rows.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in rows:
            signatures[row["signature"]].append(
                (row["gene"].strip().upper(), float(row["weight"]))
            )
    return dict(signatures)


def mean_score(
    members: list[tuple[str, float]],
    expression: dict[str, list[float]],
    sample_index: int,
) -> tuple[float, int, int]:
    present = [(expression[gene][sample_index], weight) for gene, weight in members if gene in expression]
    if not present:
        return math.nan, 0, len(members)
    denominator = sum(abs(weight) for _, weight in present)
    score = sum(value * weight for value, weight in present) / denominator
    return score, len(present), len(members)


def sample_zscores(expression: dict[str, list[float]], sample_index: int) -> dict[str, float]:
    values = [row[sample_index] for row in expression.values()]
    center = statistics.fmean(values)
    scale = statistics.stdev(values) if len(values) > 1 else 0.0
    if scale == 0:
        raise ValueError("IPS is undefined for a sample with zero expression variance")
    return {gene: (row[sample_index] - center) / scale for gene, row in expression.items()}


def read_ips(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"gene", "factor", "class", "weight"}.issubset(rows[0]):
        raise ValueError(f"{path}: malformed IPS gene table")
    return rows


def ips_score(
    expression: dict[str, list[float]],
    sample_index: int,
    rows: list[dict[str, str]],
) -> tuple[dict[str, float], int, int]:
    zscores = sample_zscores(expression, sample_index)
    factors: dict[str, list[float]] = defaultdict(list)
    factor_class: dict[str, str] = {}
    factor_weight: dict[str, float] = {}
    present_genes: set[str] = set()
    expected_genes = {row["gene"] for row in rows}
    for row in rows:
        gene = row["gene"]
        if gene not in zscores:
            continue
        factor = row["factor"]
        factors[factor].append(zscores[gene])
        factor_class[factor] = row["class"]
        factor_weight[factor] = float(row["weight"])
        present_genes.add(gene)

    class_values: dict[str, list[float]] = defaultdict(list)
    for factor, values in factors.items():
        class_values[factor_class[factor]].append(
            statistics.fmean(values) * factor_weight[factor]
        )
    scores = {
        class_name: statistics.fmean(class_values[class_name])
        if class_values[class_name]
        else math.nan
        for class_name in ("MHC", "CP", "EC", "SC")
    }
    if any(math.isnan(value) for value in scores.values()):
        aggregate = ips = math.nan
    else:
        aggregate = sum(scores.values())
        ips = 0.0 if aggregate <= 0 else 10.0 if aggregate >= 3 else round(aggregate * 10 / 3)
    scores["IPS_raw"] = aggregate
    scores["IPS"] = ips
    return scores, len(present_genes), len(expected_genes)


def format_number(value: float) -> str:
    return "NA" if math.isnan(value) else f"{value:.10g}"


def run(args: argparse.Namespace) -> None:
    samples, expression = read_expression(args.input)
    signatures = read_signatures(args.signatures)
    ips_rows = read_ips(args.ips_genes)
    selected = args.signature or list(signatures)
    unknown = sorted(set(selected) - set(signatures))
    if unknown:
        raise ValueError(f"unknown signature(s): {', '.join(unknown)}")

    output_rows: list[dict[str, str]] = []
    for index, sample in enumerate(samples):
        row: dict[str, str] = {"sample": sample}
        for name in selected:
            score, present, expected = mean_score(signatures[name], expression, index)
            row[name] = format_number(score)
            row[f"{name}_coverage"] = f"{present}/{expected}"
        ips, present, expected = ips_score(expression, index, ips_rows)
        for name, value in ips.items():
            row[name] = format_number(value)
        row["IPS_coverage"] = f"{present}/{expected}"
        output_rows.append(row)

    columns = ["sample"]
    for name in selected:
        columns.extend((name, f"{name}_coverage"))
    columns.extend(("MHC", "CP", "EC", "SC", "IPS_raw", "IPS", "IPS_coverage"))
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("input", type=Path, help="genes x samples TSV, log-scale expression")
    result.add_argument("output", type=Path, help="sample-level output TSV")
    result.add_argument(
        "--signature",
        action="append",
        help="signature to score (repeatable; default: all)",
    )
    result.add_argument("--signatures", type=Path, default=DEFAULT_SETS)
    result.add_argument("--ips-genes", type=Path, default=DEFAULT_IPS)
    return result


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
