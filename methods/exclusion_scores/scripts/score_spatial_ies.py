#!/usr/bin/env python3
"""Compute a transparent reconstruction of the Sher spatial CD8 IES.

This is not the authors' exact unpublished implementation. It computes signed
perpendicular distance from y=x in log2 stromal-versus-epithelial CD8-density
coordinates; positive values indicate stromal enrichment/epithelial exclusion.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


REQUIRED = ("sample", "epithelial_cd8_density", "stromal_cd8_density")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(REQUIRED).issubset(reader.fieldnames):
            raise ValueError("input requires columns: " + ", ".join(REQUIRED))
        rows = list(reader)
    if not rows:
        raise ValueError("input has no samples")
    return rows


def choose_epsilon(values: list[float], requested: float | None) -> float:
    if requested is not None:
        if requested <= 0:
            raise ValueError("--epsilon must be positive")
        return requested
    positives = [value for value in values if value > 0]
    if not positives:
        raise ValueError("cannot derive epsilon when all densities are zero")
    return 0.5 * min(positives)


def calculate(
    epithelial: float, stromal: float, epsilon: float
) -> tuple[float, float]:
    if epithelial < 0 or stromal < 0:
        raise ValueError("CD8 densities cannot be negative")
    ies = (
        math.log2(stromal + epsilon) - math.log2(epithelial + epsilon)
    ) / math.sqrt(2)
    return ies, epithelial + stromal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="sample-level density TSV")
    parser.add_argument("output", type=Path, help="sample-level IES TSV")
    parser.add_argument(
        "--epsilon",
        type=float,
        help="positive pseudocount (default: half smallest positive density)",
    )
    args = parser.parse_args()
    rows = read_rows(args.input)
    parsed = []
    for row in rows:
        try:
            epithelial = float(row["epithelial_cd8_density"])
            stromal = float(row["stromal_cd8_density"])
        except ValueError as exc:
            raise ValueError(f"{row['sample']}: non-numeric density") from exc
        if not math.isfinite(epithelial) or not math.isfinite(stromal):
            raise ValueError(f"{row['sample']}: non-finite density")
        parsed.append((row["sample"], epithelial, stromal))
    epsilon = choose_epsilon(
        [value for _, epithelial, stromal in parsed for value in (epithelial, stromal)],
        args.epsilon,
    )
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "sample",
            "Spatial_IES_reconstructed",
            "total_cd8_density",
            "epsilon",
            "implementation_status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sample, epithelial, stromal in parsed:
            ies, total = calculate(epithelial, stromal, epsilon)
            writer.writerow(
                {
                    "sample": sample,
                    "Spatial_IES_reconstructed": f"{ies:.10g}",
                    "total_cd8_density": f"{total:.10g}",
                    "epsilon": f"{epsilon:.10g}",
                    "implementation_status": "reconstruction_not_author_exact",
                }
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
