#!/usr/bin/env python3
"""Validate and combine TACSTD2/CLDN4 multi-omic result tables.

Validates plug-in summary tables, harmonizes inferential z statistics, and
meta-analyzes only explicitly declared, compatible groups. The script uses
only the Python standard library and never pools patient-level records.

校验多组学汇总结果表，统一推断 z 统计量，并且只对预先声明且完全兼容的组
进行合并；不读取或拼接患者层级数据。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


REQUIRED_COLUMNS = (
    "dataset_id",
    "analysis_id",
    "modality",
    "access_class",
    "evidence_tier",
    "tumor_type",
    "treatment",
    "endpoint",
    "biological_unit",
    "score_name",
    "score_version",
    "contrast",
    "effect_type",
    "estimate",
    "se",
    "n_patients",
    "direction",
    "adjustment_set",
    "analysis_tier",
    "meta_group",
    "qc_status",
    "junction_measure",
    "endpoint_definition",
    "time_origin",
    "population",
    "participant_set_id",
)

VALID_VALUES = {
    "modality": {"bulk_rna", "scrna", "protein", "spatial", "multiomic"},
    "access_class": {"open", "controlled", "mixed_summary"},
    "evidence_tier": {"A", "B", "C", "D"},
    "effect_type": {
        "beta",
        "log_or",
        "log_hr",
        "interaction_beta",
        "log_ratio",
    },
    "analysis_tier": {"primary", "sensitivity", "exploratory"},
    "qc_status": {"pass", "warn", "fail"},
}

EXPECTED_DIRECTION = "higher_is_more_exclusion_or_resistance"

# A meta_group is rejected if any of these differ. Empty optional values still
# compare literally: agents should populate them consistently.
COMPATIBILITY_FIELDS = (
    "effect_type",
    "endpoint",
    "score_name",
    "score_version",
    "contrast",
    "direction",
    "tumor_type",
    "treatment",
    "biological_unit",
    "adjustment_set",
    "analysis_tier",
    "junction_measure",
    "endpoint_definition",
    "time_origin",
    "population",
)

OUTPUT_FIELDS = list(REQUIRED_COLUMNS) + [
    "source_file",
    "source_row",
    "z_value",
    "p_value_two_sided",
    "ci_low_95",
    "ci_high_95",
    "validation_status",
    "validation_messages",
]

EXCLUSION_FIELDS = (
    "meta_group",
    "dataset_id",
    "analysis_id",
    "reason_code",
    "reason",
)

META_FIELDS = (
    "meta_group",
    "k",
    "n_patients_sum",
    "effect_type",
    "endpoint",
    "score_name",
    "score_version",
    "contrast",
    "tumor_type",
    "treatment",
    "junction_measure",
    "estimate_fixed",
    "se_fixed",
    "ci_low_95",
    "ci_high_95",
    "z_value",
    "p_value_two_sided",
    "q_heterogeneity",
    "q_df",
    "i2_percent",
    "access_classes",
    "evidence_tiers",
    "datasets",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate multi-omic result tables and combine only compatible "
            "inverse-variance estimates."
        )
    )
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        type=Path,
        help="Input CSV or TSV result tables.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory (created if absent).",
    )
    parser.add_argument(
        "--allow-warn-meta",
        action="store_true",
        help="Allow input rows whose qc_status is 'warn' into meta-analysis.",
    )
    parser.add_argument(
        "--min-studies",
        type=int,
        default=2,
        help="Minimum independent datasets per meta-analysis (default: 2).",
    )
    return parser.parse_args(argv)


def detect_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    if path.suffix.lower() == ".csv":
        return ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error as exc:
        raise ValueError(f"cannot detect CSV/TSV delimiter: {path}") from exc


def load_rows(paths: Iterable[Path]) -> Tuple[List[Dict[str, str]], List[str]]:
    rows: List[Dict[str, str]] = []
    file_errors: List[str] = []
    for path in paths:
        if not path.is_file():
            file_errors.append(f"{path}: file not found")
            continue
        try:
            delimiter = detect_delimiter(path)
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter=delimiter)
                headers = reader.fieldnames or []
                missing = [column for column in REQUIRED_COLUMNS if column not in headers]
                if missing:
                    file_errors.append(
                        f"{path}: missing required columns: {', '.join(missing)}"
                    )
                    continue
                for row_number, raw in enumerate(reader, start=2):
                    row = {
                        str(key).strip(): (value or "").strip()
                        for key, value in raw.items()
                        if key is not None
                    }
                    row["source_file"] = str(path)
                    row["source_row"] = str(row_number)
                    rows.append(row)
        except (OSError, UnicodeError, csv.Error, ValueError) as exc:
            file_errors.append(f"{path}: {exc}")
    return rows, file_errors


def finite_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("not finite")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise ValueError("not positive")
    return number


def add_message(messages: List[Tuple[str, str]], level: str, text: str) -> None:
    messages.append((level, text))


def validate_rows(
    rows: List[Dict[str, str]],
) -> Tuple[List[Dict[str, str]], Counter]:
    analysis_counts = Counter(row.get("analysis_id", "") for row in rows)
    status_counts: Counter = Counter()

    for row in rows:
        messages: List[Tuple[str, str]] = []
        for column in REQUIRED_COLUMNS:
            # meta_group may intentionally be blank to prohibit synthesis.
            if column != "meta_group" and not row.get(column, ""):
                add_message(messages, "fail", f"missing:{column}")

        for column, allowed in VALID_VALUES.items():
            value = row.get(column, "")
            if value and value not in allowed:
                add_message(messages, "fail", f"invalid_{column}:{value}")

        estimate = se = None
        try:
            estimate = finite_float(row.get("estimate", ""))
        except (TypeError, ValueError):
            add_message(messages, "fail", "estimate_not_finite")
        try:
            se = finite_float(row.get("se", ""))
            if se <= 0:
                raise ValueError("not positive")
        except (TypeError, ValueError):
            add_message(messages, "fail", "se_must_be_finite_and_positive")
        try:
            positive_int(row.get("n_patients", ""))
        except (TypeError, ValueError):
            add_message(messages, "fail", "n_patients_must_be_positive_integer")

        if row.get("direction") and row["direction"] != EXPECTED_DIRECTION:
            add_message(messages, "fail", "effect_direction_not_harmonized")
        if analysis_counts[row.get("analysis_id", "")] > 1:
            add_message(messages, "fail", "duplicate_analysis_id")
        if row.get("biological_unit") != "patient":
            add_message(messages, "warn", "non_patient_biological_unit")
        if not row.get("meta_group"):
            add_message(messages, "warn", "blank_meta_group_no_synthesis")
        if row.get("qc_status") == "fail":
            add_message(messages, "fail", "upstream_qc_status_fail")
        elif row.get("qc_status") == "warn":
            add_message(messages, "warn", "upstream_qc_status_warn")

        if estimate is not None and se is not None and se > 0:
            z_value = estimate / se
            row["z_value"] = format_number(z_value)
            row["p_value_two_sided"] = format_number(math.erfc(abs(z_value) / math.sqrt(2)))
            row["ci_low_95"] = format_number(estimate - 1.959963984540054 * se)
            row["ci_high_95"] = format_number(estimate + 1.959963984540054 * se)
        else:
            for field in ("z_value", "p_value_two_sided", "ci_low_95", "ci_high_95"):
                row[field] = ""

        levels = {level for level, _ in messages}
        status = "fail" if "fail" in levels else "warn" if "warn" in levels else "pass"
        row["validation_status"] = status
        row["validation_messages"] = ";".join(text for _, text in messages)
        status_counts[status] += 1

    return rows, status_counts


def format_number(value: float) -> str:
    return f"{value:.12g}"


def exclusion(
    row: Dict[str, str], reason_code: str, reason: str
) -> Dict[str, str]:
    return {
        "meta_group": row.get("meta_group", ""),
        "dataset_id": row.get("dataset_id", ""),
        "analysis_id": row.get("analysis_id", ""),
        "reason_code": reason_code,
        "reason": reason,
    }


def compatibility_differences(rows: List[Dict[str, str]]) -> List[str]:
    differences = []
    for field in COMPATIBILITY_FIELDS:
        values = {row.get(field, "") for row in rows}
        if len(values) > 1:
            differences.append(f"{field}={sorted(values)}")
    return differences


def synthesize_groups(
    rows: List[Dict[str, str]], allow_warn: bool, min_studies: int
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    groups: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    exclusions: List[Dict[str, str]] = []
    summaries: List[Dict[str, str]] = []

    for row in rows:
        group = row.get("meta_group", "")
        if not group:
            exclusions.append(exclusion(row, "blank_meta_group", "No synthesis requested"))
        else:
            groups[group].append(row)

    for group_name, candidates in sorted(groups.items()):
        eligible: List[Dict[str, str]] = []
        for row in candidates:
            status = row["validation_status"]
            if status == "fail":
                exclusions.append(exclusion(row, "validation_fail", row["validation_messages"]))
            elif status == "warn" and not allow_warn:
                exclusions.append(
                    exclusion(
                        row,
                        "validation_warn",
                        "Warn row excluded; use --allow-warn-meta after review",
                    )
                )
            elif row.get("biological_unit") != "patient":
                exclusions.append(
                    exclusion(
                        row,
                        "non_patient_unit",
                        "Meta-analysis requires patient-level inferential units",
                    )
                )
            else:
                eligible.append(row)

        if not eligible:
            continue

        differences = compatibility_differences(eligible)
        if differences:
            reason = "Incompatible rows in declared meta_group: " + "; ".join(differences)
            exclusions.extend(
                exclusion(row, "incompatible_estimands", reason) for row in eligible
            )
            continue

        dataset_counts = Counter(row["dataset_id"] for row in eligible)
        repeated_datasets = sorted(
            dataset for dataset, count in dataset_counts.items() if count > 1
        )
        participant_ids = [
            row.get("participant_set_id", "")
            for row in eligible
            if row.get("participant_set_id", "")
        ]
        repeated_participants = sorted(
            value for value, count in Counter(participant_ids).items() if count > 1
        )
        if repeated_datasets or repeated_participants:
            details = []
            if repeated_datasets:
                details.append("repeated dataset_id=" + ",".join(repeated_datasets))
            if repeated_participants:
                details.append(
                    "repeated participant_set_id=" + ",".join(repeated_participants)
                )
            reason = "Potentially non-independent inputs: " + "; ".join(details)
            exclusions.extend(
                exclusion(row, "non_independent_inputs", reason) for row in eligible
            )
            continue

        if len(eligible) < min_studies:
            reason = f"Only {len(eligible)} eligible dataset(s); minimum is {min_studies}"
            exclusions.extend(
                exclusion(row, "too_few_studies", reason) for row in eligible
            )
            continue

        estimates = [finite_float(row["estimate"]) for row in eligible]
        ses = [finite_float(row["se"]) for row in eligible]
        weights = [1.0 / (se * se) for se in ses]
        weight_sum = sum(weights)
        pooled = sum(weight * estimate for weight, estimate in zip(weights, estimates)) / weight_sum
        pooled_se = math.sqrt(1.0 / weight_sum)
        z_value = pooled / pooled_se
        q_value = sum(
            weight * (estimate - pooled) ** 2
            for weight, estimate in zip(weights, estimates)
        )
        q_df = len(eligible) - 1
        i2 = max(0.0, (q_value - q_df) / q_value * 100.0) if q_value > 0 else 0.0
        first = eligible[0]

        summaries.append(
            {
                "meta_group": group_name,
                "k": str(len(eligible)),
                "n_patients_sum": str(sum(positive_int(row["n_patients"]) for row in eligible)),
                "effect_type": first["effect_type"],
                "endpoint": first["endpoint"],
                "score_name": first["score_name"],
                "score_version": first["score_version"],
                "contrast": first["contrast"],
                "tumor_type": first["tumor_type"],
                "treatment": first["treatment"],
                "junction_measure": first.get("junction_measure", ""),
                "estimate_fixed": format_number(pooled),
                "se_fixed": format_number(pooled_se),
                "ci_low_95": format_number(pooled - 1.959963984540054 * pooled_se),
                "ci_high_95": format_number(pooled + 1.959963984540054 * pooled_se),
                "z_value": format_number(z_value),
                "p_value_two_sided": format_number(
                    math.erfc(abs(z_value) / math.sqrt(2))
                ),
                "q_heterogeneity": format_number(q_value),
                "q_df": str(q_df),
                "i2_percent": format_number(i2),
                "access_classes": "|".join(sorted({row["access_class"] for row in eligible})),
                "evidence_tiers": "|".join(sorted({row["evidence_tier"] for row in eligible})),
                "datasets": "|".join(sorted(row["dataset_id"] for row in eligible)),
            }
        )

    return summaries, exclusions


def write_csv(path: Path, rows: List[Dict[str, str]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.min_studies < 2:
        print("--min-studies must be at least 2", file=sys.stderr)
        return 2

    rows, file_errors = load_rows(args.input)
    if not rows:
        for error in file_errors:
            print(error, file=sys.stderr)
        print("No valid input rows loaded.", file=sys.stderr)
        return 2

    rows, status_counts = validate_rows(rows)
    summaries, exclusions = synthesize_groups(
        rows, allow_warn=args.allow_warn_meta, min_studies=args.min_studies
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "harmonized_results.csv", rows, OUTPUT_FIELDS)
    write_csv(args.outdir / "meta_summary.csv", summaries, META_FIELDS)
    write_csv(args.outdir / "meta_exclusions.csv", exclusions, EXCLUSION_FIELDS)

    report = {
        "input_files": [str(path) for path in args.input],
        "input_rows": len(rows),
        "validation_counts": dict(sorted(status_counts.items())),
        "meta_groups_synthesized": len(summaries),
        "meta_exclusion_rows": len(exclusions),
        "allow_warn_meta": args.allow_warn_meta,
        "min_studies": args.min_studies,
        "file_errors": file_errors,
        "notes": [
            "z_value is estimate/se, not a cross-dataset patient score.",
            "Fixed-effect synthesis is performed only inside explicit compatible meta_group values.",
            "Review meta_exclusions.csv and all warn rows before interpretation.",
        ],
    }
    with (args.outdir / "qc_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 1 if file_errors or status_counts.get("fail", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
