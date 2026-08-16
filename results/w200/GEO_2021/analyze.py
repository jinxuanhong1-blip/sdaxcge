#!/usr/bin/env python3
"""Audit TACSTD2/CLDN4-positive cells in the GSE146100 normalized matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import sys
from collections import Counter
from pathlib import Path


EXPECTED_SHA256 = "b937b797bd701da06e1b616f72364ee9fdaecb6d3cf74be2b7ae667614c378b6"
TARGETS = ("TACSTD2", "CLDN4")
EPITHELIAL_SUPPORT = ("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1")
CONTEXT = ("PTPRC",)
GENES = TARGETS + EPITHELIAL_SUPPORT + CONTEXT


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize TACSTD2/CLDN4 detection by GSE146100 lesion."
    )
    parser.add_argument("matrix", type=Path, help="GSE146100_NormData.txt.gz")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("marker_summary.csv"),
        help="Output CSV (default: beside this script)",
    )
    parser.add_argument(
        "--skip-checksum",
        action="store_true",
        help="Allow a matrix whose SHA-256 differs from the archived copy.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_selected_rows(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    selected: dict[str, list[float]] = {}
    with gzip.open(path, "rt", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "Gene":
            raise ValueError("Expected a tab-delimited gene-by-cell matrix")
        cells = header[1:]
        for row in reader:
            if row and row[0] in GENES:
                if row[0] in selected:
                    raise ValueError(f"Duplicate gene row: {row[0]}")
                if len(row) != len(header):
                    raise ValueError(f"Wrong field count for {row[0]}")
                selected[row[0]] = [float(value) for value in row[1:]]

    missing = sorted(set(GENES) - selected.keys())
    if missing:
        raise ValueError(f"Missing required genes: {', '.join(missing)}")
    return cells, selected


def percentage(count: int, total: int) -> float:
    return 100.0 * count / total


def summarize(cells: list[str], values: dict[str, list[float]]) -> list[dict[str, object]]:
    lesions = [cell.split("_", 1)[0] for cell in cells]
    unexpected = sorted(set(lesions) - {"W1", "W2", "W3"})
    if unexpected:
        raise ValueError(f"Unexpected lesion prefixes: {', '.join(unexpected)}")

    rows: list[dict[str, object]] = []
    for lesion in ("W1", "W2", "W3", "all"):
        indices = [
            index
            for index, observed_lesion in enumerate(lesions)
            if lesion == "all" or lesion == observed_lesion
        ]
        n = len(indices)
        positive = {
            gene: {index for index in indices if values[gene][index] > 0}
            for gene in GENES
        }
        pair = positive["TACSTD2"] & positive["CLDN4"]
        supported = {
            index
            for index in pair
            if sum(index in positive[gene] for gene in EPITHELIAL_SUPPORT) >= 2
        }
        ptprc_pair = pair & positive["PTPRC"]
        row: dict[str, object] = {"lesion": lesion, "cells": n}
        for gene in TARGETS:
            count = len(positive[gene])
            row[f"{gene}_positive"] = count
            row[f"{gene}_percent"] = f"{percentage(count, n):.2f}"
        row["pair_positive"] = len(pair)
        row["pair_percent"] = f"{percentage(len(pair), n):.2f}"
        row["pair_epithelial_supported"] = len(supported)
        row["pair_epithelial_supported_percent"] = f"{percentage(len(supported), n):.2f}"
        row["pair_PTPRC_positive"] = len(ptprc_pair)
        row["pair_PTPRC_positive_percent"] = f"{percentage(len(ptprc_pair), n):.2f}"
        rows.append(row)
    return rows


def main() -> int:
    args = parse_args()
    observed_sha256 = sha256(args.matrix)
    if not args.skip_checksum and observed_sha256 != EXPECTED_SHA256:
        print(
            f"Checksum mismatch: expected {EXPECTED_SHA256}, got {observed_sha256}",
            file=sys.stderr,
        )
        return 2

    cells, values = read_selected_rows(args.matrix)
    lesion_counts = Counter(cell.split("_", 1)[0] for cell in cells)
    if sum(lesion_counts.values()) != 11_612:
        raise ValueError(f"Expected 11,612 cells, found {len(cells):,}")

    rows = summarize(cells, values)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output} from {len(cells):,} cells ({observed_sha256})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
