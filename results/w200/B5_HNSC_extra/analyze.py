#!/usr/bin/env python3
"""Reproduce the deliberately small GSE93157 HNSCC analysis without dependencies."""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import math
from pathlib import Path
from statistics import mean


HERE = Path(__file__).resolve().parent
MATRIX = HERE / "data" / "GSE93157_series_matrix.txt.gz"
OR_LABELS = {"CR", "PR"}
CYT_GENES = ("GZMA", "PRF1")


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_matrix(path: Path) -> tuple[list[dict[str, str]], dict[str, list[float | None]]]:
    metadata: dict[str, list[str]] = {}
    characteristics: dict[str, list[str]] = {}
    expression: dict[str, list[float | None]] = {}
    in_table = False
    sample_ids: list[str] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                if line.startswith("!Sample_"):
                    fields = split_geo_line(line)
                    key, values = fields[0], fields[1:]
                    if key == "!Sample_characteristics_ch1":
                        characteristic = values[0].split(":", 1)[0].lower()
                        characteristics[characteristic] = [
                            value.split(":", 1)[1].strip() for value in values
                        ]
                    else:
                        metadata[key] = values
                continue

            fields = split_geo_line(line)
            if fields[0] == "ID_REF":
                sample_ids = fields[1:]
                continue
            expression[fields[0]] = [
                float(value) if value else None for value in fields[1:]
            ]

    samples = []
    for index, sample_id in enumerate(sample_ids):
        sample = {
            "sample_id": sample_id,
            "title": metadata["!Sample_title"][index],
            "cancer_type": metadata["!Sample_source_name_ch1"][index],
        }
        for key, values in characteristics.items():
            sample[key] = values[index]
        samples.append(sample)
    return samples, expression


def mean_difference(values: list[float], responder_indices: tuple[int, ...]) -> float:
    responders = [values[index] for index in responder_indices]
    others = [value for index, value in enumerate(values) if index not in responder_indices]
    return mean(responders) - mean(others)


def exact_mean_difference_test(
    values: list[float], responder_indices: tuple[int, ...]
) -> tuple[float, float]:
    observed = mean_difference(values, responder_indices)
    null = [
        mean_difference(values, group)
        for group in itertools.combinations(range(len(values)), len(responder_indices))
    ]
    p_value = sum(abs(value) >= abs(observed) - 1e-12 for value in null) / len(null)
    return observed, p_value


def ranks(values: list[float]) -> list[float]:
    result = [0.0] * len(values)
    order = sorted(range(len(values)), key=values.__getitem__)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        average_rank = (position + 1 + end) / 2
        for index in order[position:end]:
            result[index] = average_rank
        position = end
    return result


def correlation(left: list[float], right: list[float]) -> float:
    left_mean, right_mean = mean(left), mean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    denominator = math.sqrt(
        sum((x - left_mean) ** 2 for x in left)
        * sum((y - right_mean) ** 2 for y in right)
    )
    return numerator / denominator if denominator else 0.0


def exact_spearman(values: list[float], outcomes: list[float]) -> tuple[float, float]:
    value_ranks = ranks(values)
    outcome_ranks = ranks(outcomes)
    observed = correlation(value_ranks, outcome_ranks)
    null = [
        correlation(value_ranks, list(permutation))
        for permutation in itertools.permutations(outcome_ranks)
    ]
    p_value = sum(abs(value) >= abs(observed) - 1e-12 for value in null) / len(null)
    return observed, p_value


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(order, 1))):
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    samples, expression = parse_matrix(MATRIX)
    hnsc_indices = [
        index for index, sample in enumerate(samples) if sample["cancer_type"] == "HEADNECK"
    ]
    hnsc = [samples[index] for index in hnsc_indices]
    responder_indices = tuple(
        index for index, sample in enumerate(hnsc) if sample["best.resp"] in OR_LABELS
    )
    pfs = [float(sample["pfs"]) for sample in hnsc]

    patient_rows: list[dict[str, object]] = []
    for sample in hnsc:
        patient_rows.append(
            {
                "sample_id": sample["sample_id"],
                "patient": sample["title"],
                "best_response": sample["best.resp"],
                "objective_response": sample["best.resp"] in OR_LABELS,
                "pfs_months": sample["pfs"],
                "pfs_event": sample["pfse"],
                "biopsy": sample["biopsy"],
                "drug": sample["drug"],
                "age": sample["age"],
                "sex": sample["sex"],
                "ecog": sample["ecog"],
            }
        )
    write_csv(HERE / "patients.csv", patient_rows, list(patient_rows[0]))

    gene_rows: list[dict[str, object]] = []
    for gene, all_values in expression.items():
        values = [all_values[index] for index in hnsc_indices]
        if any(value is None for value in values):
            continue
        numeric = [float(value) for value in values if value is not None]
        difference, response_p = exact_mean_difference_test(numeric, responder_indices)
        rho, pfs_p = exact_spearman(numeric, pfs)
        gene_rows.append(
            {
                "gene": gene,
                "or_minus_non_or_log2": difference,
                "or_exact_p": response_p,
                "pfs_spearman_rho": rho,
                "pfs_exact_p": pfs_p,
            }
        )

    response_q = benjamini_hochberg([float(row["or_exact_p"]) for row in gene_rows])
    pfs_q = benjamini_hochberg([float(row["pfs_exact_p"]) for row in gene_rows])
    for row, or_q, survival_q in zip(gene_rows, response_q, pfs_q):
        row["or_bh_q"] = or_q
        row["pfs_bh_q"] = survival_q
    gene_rows.sort(key=lambda row: (-abs(float(row["or_minus_non_or_log2"])), row["gene"]))
    write_csv(
        HERE / "gene_statistics.csv",
        gene_rows,
        [
            "gene",
            "or_minus_non_or_log2",
            "or_exact_p",
            "or_bh_q",
            "pfs_spearman_rho",
            "pfs_exact_p",
            "pfs_bh_q",
        ],
    )

    cyt_values = [
        mean([float(expression[gene][index]) for gene in CYT_GENES])
        for index in hnsc_indices
    ]
    cyt_difference, cyt_response_p = exact_mean_difference_test(cyt_values, responder_indices)
    cyt_rho, cyt_pfs_p = exact_spearman(cyt_values, pfs)
    cyt_rows = [
        {
            "sample_id": sample["sample_id"],
            "best_response": sample["best.resp"],
            "pfs_months": sample["pfs"],
            "cyt_log2_mean": value,
        }
        for sample, value in zip(hnsc, cyt_values)
    ]
    write_csv(
        HERE / "cytolytic_score.csv",
        cyt_rows,
        ["sample_id", "best_response", "pfs_months", "cyt_log2_mean"],
    )

    summary = {
        "accession": "GSE93157",
        "platform": "NanoString PanCancer Immune Profiling Panel",
        "hnsc_n": len(hnsc),
        "objective_responders_n": len(responder_indices),
        "objective_nonresponders_n": len(hnsc) - len(responder_indices),
        "complete_case_genes_n": len(gene_rows),
        "minimum_possible_two_sided_or_p": 0.1,
        "genes_or_bh_q_below_0_05_n": sum(
            float(row["or_bh_q"]) < 0.05 for row in gene_rows
        ),
        "genes_pfs_bh_q_below_0_05_n": sum(
            float(row["pfs_bh_q"]) < 0.05 for row in gene_rows
        ),
        "cytolytic_score": {
            "genes": list(CYT_GENES),
            "or_minus_non_or_log2": cyt_difference,
            "or_exact_p": cyt_response_p,
            "pfs_spearman_rho": cyt_rho,
            "pfs_exact_p": cyt_pfs_p,
        },
        "top_or_effects": gene_rows[:10],
    }
    (HERE / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
