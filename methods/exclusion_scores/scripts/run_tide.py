#!/usr/bin/env python3
"""Validate an expression matrix and run the canonical TIDEpy NSCLC model."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


REQUIRED_OUTPUTS = {
    "TIDE",
    "IFNG",
    "Dysfunction",
    "Exclusion",
    "MDSC",
    "CAF",
    "TAM M2",
}


def validate_input(path: Path) -> None:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter="\t")
        try:
            header = next(rows)
        except StopIteration:
            raise ValueError("input is empty")
        if len(header) < 3:
            raise ValueError("TIDE requires a gene column and preferably >=2 samples")
        genes = 0
        for line_no, row in enumerate(rows, 2):
            if len(row) != len(header):
                raise ValueError(f"line {line_no}: inconsistent number of columns")
            try:
                [float(value) for value in row[1:]]
            except ValueError as exc:
                raise ValueError(f"line {line_no}: non-numeric expression") from exc
            genes += 1
    if genes < 1000:
        raise ValueError(
            f"TIDE needs whole-transcriptome input; found only {genes} gene rows"
        )


def validate_output(path: Path) -> None:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
    missing = REQUIRED_OUTPUTS - fields
    if missing:
        raise ValueError(f"TIDE output is missing: {', '.join(sorted(missing))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="genes x samples expression TSV")
    parser.add_argument("output", type=Path, help="canonical TIDE output TSV")
    parser.add_argument("--pretreat", action="store_true", help="prior immunotherapy")
    parser.add_argument(
        "--executable",
        default="tidepy",
        help="TIDEpy executable (default: tidepy)",
    )
    args = parser.parse_args()

    validate_input(args.input)
    executable = shutil.which(args.executable)
    if not executable:
        raise ValueError(
            "canonical TIDEpy is not installed; install the authors' package "
            "from https://github.com/liulab-dfci/TIDEpy (v1.3), then rerun"
        )
    command = [
        executable,
        str(args.input),
        "-o",
        str(args.output),
        "-c",
        "NSCLC",
    ]
    if args.pretreat:
        command.append("--pretreat")
    subprocess.run(command, check=True)
    validate_output(args.output)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
