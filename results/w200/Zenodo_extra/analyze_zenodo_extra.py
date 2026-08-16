#!/usr/bin/env python3
"""Audit target coverage in the leftover processed Zenodo NSCLC ICI cohort."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
TARGETS = ("TACSTD2", "CLDN4")
CONTEXT_GENES = ("EPCAM", "CDH1", "MUC1", "CEACAM6")
EXPECTED_MD5 = {
    "learn-data.csv": "90dcff5eba025880a7b55a925165ff96",
    "learn-metadata.csv": "57f426b6e54f47827fd88b5481e761b9",
    "validate-data.csv": "69d0f9ca7472ff459eebdc467e166764",
    "validate-metadata.csv": "4a9016ad18ce0cda0cad352182ddc382",
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def count_text(rows: list[dict[str, str]], field: str) -> str:
    counts = Counter((row.get(field) or "missing").strip() or "missing" for row in rows)
    return json.dumps(dict(sorted(counts.items())), sort_keys=True)


def main() -> None:
    manifest_rows: list[dict[str, object]] = []
    for name, expected in EXPECTED_MD5.items():
        observed = md5(DATA / name)
        if observed != expected:
            raise RuntimeError(f"Checksum mismatch for {name}: {observed} != {expected}")
        manifest_rows.append(
            {
                "file": name,
                "md5": observed,
                "bytes": (DATA / name).stat().st_size,
                "verified": True,
            }
        )

    coverage_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for cohort, stem in (("discovery", "learn"), ("validation", "validate")):
        expression_fields, expression = read_csv(DATA / f"{stem}-data.csv")
        metadata_fields, metadata = read_csv(DATA / f"{stem}-metadata.csv")
        if len(expression) != len(metadata):
            raise RuntimeError(
                f"{cohort}: expression rows ({len(expression)}) != metadata rows ({len(metadata)})"
            )
        if [row[""] for row in expression] != [row[""] for row in metadata]:
            raise RuntimeError(f"{cohort}: expression and metadata row identifiers do not align")

        genes = set(expression_fields)
        for gene in TARGETS + CONTEXT_GENES:
            coverage_rows.append(
                {
                    "cohort": cohort,
                    "gene": gene,
                    "requested_target": gene in TARGETS,
                    "present_exact_symbol": gene in genes,
                    "n_samples": len(expression),
                }
            )

        summary_rows.append(
            {
                "cohort": cohort,
                "n_samples": len(metadata),
                "n_expression_features_including_controls": len(expression_fields) - 1,
                "best_response_binary_counts": count_text(metadata, "best response binary"),
                "best_response_raw_counts": count_text(metadata, "best response IO"),
                "histology_counts": count_text(metadata, "entity"),
                "metadata_columns": len(metadata_fields) - 1,
            }
        )

    write_csv(
        ROOT / "source_manifest.csv",
        manifest_rows,
        ["file", "md5", "bytes", "verified"],
    )
    write_csv(
        ROOT / "target_coverage.csv",
        coverage_rows,
        ["cohort", "gene", "requested_target", "present_exact_symbol", "n_samples"],
    )
    write_csv(
        ROOT / "cohort_summary.csv",
        summary_rows,
        [
            "cohort",
            "n_samples",
            "n_expression_features_including_controls",
            "best_response_binary_counts",
            "best_response_raw_counts",
            "histology_counts",
            "metadata_columns",
        ],
    )

    result = {
        "record": "10.5281/zenodo.2635194",
        "assay": "NanoString nCounter PanCancer Immune Profiling panel",
        "population": "advanced NSCLC treated with anti-PD-1 immunotherapy",
        "n_discovery": summary_rows[0]["n_samples"],
        "n_validation": summary_rows[1]["n_samples"],
        "targets": {
            target: {
                "present_discovery": any(
                    row["cohort"] == "discovery"
                    and row["gene"] == target
                    and row["present_exact_symbol"]
                    for row in coverage_rows
                ),
                "present_validation": any(
                    row["cohort"] == "validation"
                    and row["gene"] == target
                    and row["present_exact_symbol"]
                    for row in coverage_rows
                ),
            }
            for target in TARGETS
        },
        "valid_target_analysis": False,
        "reason": "Neither requested gene is measured by the deposited targeted panel.",
        "surrogate_analysis_performed": False,
    }
    (ROOT / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
