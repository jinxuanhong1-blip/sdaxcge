#!/usr/bin/env python3
"""Audit labels and summarize TACSTD2/CLDN4 in downloaded processed data."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

MOUSE_TARGETS = {
    "TACSTD2": "ENSMUSG00000051397",
    "CLDN4": "ENSMUSG00000047501",
}


def read_sdrf(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unique_assays(rows: list[dict]) -> list[dict]:
    seen = {}
    for row in rows:
        key = row.get("Comment[ENA_SAMPLE]") or row.get("Source Name")
        seen.setdefault(key, row)
    return list(seen.values())


def label_audit(root: Path) -> list[dict]:
    specs = {
        "E-MTAB-10633": ["Factor Value[compound]", "Factor Value[irradiate]"],
        "E-MTAB-13704": ["Factor Value[stimulus]", "Characteristics[tumours_numbers]"],
        "E-MTAB-15883": ["treatment", "cluster"],
        "E-MTAB-8867": ["Characteristics[clinical history]", "Factor Value[compound]", "Factor Value[disease staging]"],
        "E-MTAB-9451": ["Characteristics[individual]", "Characteristics[organism part]", "Factor Value[disease]"],
    }
    records = []
    for accession, fields in specs.items():
        if accession == "E-MTAB-15883":
            path = root / accession / "cell_metadata.tsv"
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            level = "cell"
        else:
            path = root / accession / f"{accession}.sdrf.txt"
            rows = unique_assays(read_sdrf(path))
            level = "sample"
        for field in fields:
            values = sorted({row.get(field, "") for row in rows if row.get(field, "")})
            records.append({
                "accession": accession,
                "level": level,
                "n_records": len(rows),
                "label": field,
                "n_values": len(values),
                "values": " | ".join(values),
            })
    return records


def analyze_13704(root: Path) -> list[dict]:
    accession = "E-MTAB-13704"
    sdrf = read_sdrf(root / accession / f"{accession}.sdrf.txt")
    sample_labels = {}
    for row in sdrf:
        column = row["Scan Name"].split("_", 1)[0]
        sample_labels[column] = {
            "group": row["Factor Value[stimulus]"],
            "sample": row["Source Name"],
        }

    matrix_path = root / accession / "GEMMS_raw_counts.csv"
    target_rows = {}
    with matrix_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for symbol, ensembl in MOUSE_TARGETS.items():
                if row["rowname"].split(".", 1)[0] == ensembl:
                    target_rows[symbol] = row
    output = []
    for symbol in MOUSE_TARGETS:
        row = target_rows.get(symbol)
        if row is None:
            output.append({"accession": accession, "gene": symbol, "group": "not_found", "n": 0})
            continue
        grouped = defaultdict(list)
        for column, labels in sample_labels.items():
            grouped[labels["group"]].append(float(row[column]))
        vehicle_mean = statistics.mean(grouped["vehicle"])
        for group, values in sorted(grouped.items()):
            output.append({
                "accession": accession,
                "gene": symbol,
                "group": group,
                "n": len(values),
                "mean": statistics.mean(values),
                "median": statistics.median(values),
                "min": min(values),
                "max": max(values),
                "log2_mean_ratio_vs_vehicle": math.log2((statistics.mean(values) + 1) / (vehicle_mean + 1)),
            })
    return output


def analyze_15883(root: Path) -> tuple[list[dict], dict]:
    accession = "E-MTAB-15883"
    base = root / accession
    feature_rows = {}
    with gzip.open(base / "features.tsv.gz", "rt") as handle:
        for index, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            symbol = fields[1].upper()
            if symbol in MOUSE_TARGETS:
                feature_rows[index] = symbol

    metadata = {}
    with (base / "cell_metadata.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["cell_barcode"]] = row

    relevant_columns = {}
    barcode_count = 0
    with gzip.open(base / "barcodes.tsv.gz", "rt") as handle:
        for index, line in enumerate(handle, 1):
            barcode_count = index
            barcode = line.rstrip("\n")
            if barcode in metadata:
                relevant_columns[index] = metadata[barcode]

    sums = defaultdict(float)
    nonzero = defaultdict(int)
    with gzip.open(base / "matrix.mtx.gz", "rt") as handle:
        for line in handle:
            if not line or line.startswith("%"):
                continue
            row_s, col_s, value_s = line.split()
            row_index, col_index = int(row_s), int(col_s)
            symbol = feature_rows.get(row_index)
            cell = relevant_columns.get(col_index)
            if symbol is None or cell is None:
                continue
            key = (symbol, cell["treatment"])
            value = float(value_s)
            sums[key] += value
            if value:
                nonzero[key] += 1

    group_sizes = defaultdict(int)
    for row in metadata.values():
        group_sizes[row["treatment"]] += 1
    output = []
    for symbol in MOUSE_TARGETS:
        for group in sorted(group_sizes):
            key = (symbol, group)
            output.append({
                "accession": accession,
                "gene": symbol,
                "group": group,
                "n": group_sizes[group],
                "mean": sums[key] / group_sizes[group],
                "nonzero_cells": nonzero[key],
                "fraction_nonzero": nonzero[key] / group_sizes[group],
            })
    diagnostics = {
        "feature_rows": {symbol: index for index, symbol in feature_rows.items()},
        "matrix_barcode_count": barcode_count,
        "annotated_cell_count": len(metadata),
        "matched_annotated_cells": len(relevant_columns),
    }
    return output, diagnostics


def audit_9451_panel(root: Path) -> dict:
    accession = "E-MTAB-9451"
    files = sorted((root / accession).glob("Combined_*_DBEC_MolsPerCell.csv"))
    panels = []
    for path in files:
        with path.open(newline="") as handle:
            for line in handle:
                if line.startswith("Cell_Index,"):
                    genes = {field.split("|", 1)[0].upper() for field in next(csv.reader([line]))[1:]}
                    panels.append(genes)
                    break
    return {
        "file_count": len(files),
        "same_panel_in_all_files": all(panel == panels[0] for panel in panels),
        "panel_gene_count": len(panels[0]) if panels else 0,
        "targets_assayed": {gene: bool(panels and gene in panels[0]) for gene in MOUSE_TARGETS},
        "interpretation": "Targeted immune panel; TACSTD2 and CLDN4 cannot be analyzed when absent from the panel.",
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("results/gpt_arrayexpress/downloads"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/gpt_arrayexpress"))
    args = parser.parse_args()

    labels = label_audit(args.root)
    expression = analyze_13704(args.root)
    expr_15883, diagnostics = analyze_15883(args.root)
    expression.extend(expr_15883)
    panel = audit_9451_panel(args.root)

    write_csv(args.out_dir / "label_audit.csv", labels)
    write_csv(args.out_dir / "target_expression_by_group.csv", expression)
    summary = {
        "E-MTAB-13704": {
            "analysis": "TACSTD2/CLDN4 abundance by treatment from the deposited processed matrix",
            "warning": "Despite the filename GEMMS_raw_counts.csv, values are non-integers; results are descriptive and no raw-count test was run.",
        },
        "E-MTAB-15883": {
            "analysis": "TACSTD2/CLDN4 abundance by treatment among deposited annotated cells",
            "diagnostics": diagnostics,
            "warning": "This is an immune-cell-focused dataset; absence/low expression is not evidence of tumor-cell absence.",
        },
        "E-MTAB-9451": panel,
        "metadata_only": ["E-MTAB-10633", "E-MTAB-8867"],
    }
    (args.out_dir / "analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote analysis outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
