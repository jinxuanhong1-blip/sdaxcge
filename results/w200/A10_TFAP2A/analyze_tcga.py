#!/usr/bin/env python3
"""Reproduce the public TCGA lung expression comparison used in this folder."""

from __future__ import annotations

import csv
import json
import math
import statistics
import urllib.request
from pathlib import Path


API = "https://www.cbioportal.org/api"
GENES = {"TFAP2A": 7020, "TACSTD2": 4070, "CLDN4": 1364}
COHORTS = {
    "LUAD": "luad_tcga_pan_can_atlas_2018",
    "LUSC": "lusc_tcga_pan_can_atlas_2018",
}
OUT = Path(__file__).resolve().parent


def fetch(cohort_id: str) -> list[dict]:
    profile = f"{cohort_id}_rna_seq_v2_mrna"
    payload = json.dumps(
        {
            "entrezGeneIds": list(GENES.values()),
            "sampleListId": f"{cohort_id}_all",
        }
    ).encode()
    request = urllib.request.Request(
        f"{API}/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY",
        data=payload,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def quantile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    index = (len(ordered) - 1) * q
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = (start + end - 1) / 2 + 1
        for index in order[start:end]:
            result[index] = average_rank
        start = end
    return result


def pearson(left: list[float], right: list[float]) -> float:
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right)
    )
    denominator = math.sqrt(
        sum((x - left_mean) ** 2 for x in left)
        * sum((y - right_mean) ** 2 for y in right)
    )
    return numerator / denominator


def write_tsv(name: str, fieldnames: list[str], rows: list[dict]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    expression_rows: list[dict] = []
    summary_rows: list[dict] = []
    correlation_rows: list[dict] = []
    overlap_rows: list[dict] = []

    id_to_symbol = {entrez: symbol for symbol, entrez in GENES.items()}
    for cohort, cohort_id in COHORTS.items():
        records = fetch(cohort_id)
        by_sample: dict[str, dict[str, float]] = {}
        for record in records:
            symbol = id_to_symbol[record["entrezGeneId"]]
            by_sample.setdefault(record["sampleId"], {})[symbol] = record["value"]
        complete = {
            sample: values
            for sample, values in by_sample.items()
            if set(values) == set(GENES)
        }

        for sample, values in sorted(complete.items()):
            expression_rows.append(
                {"cohort": cohort, "sample_id": sample, **values}
            )

        for symbol in GENES:
            values = [sample[symbol] for sample in complete.values()]
            summary_rows.append(
                {
                    "cohort": cohort,
                    "gene": symbol,
                    "n": len(values),
                    "median_rsem": f"{statistics.median(values):.4f}",
                    "q1_rsem": f"{quantile(values, 0.25):.4f}",
                    "q3_rsem": f"{quantile(values, 0.75):.4f}",
                    "zero_count": sum(value == 0 for value in values),
                }
            )

        pairs = [
            ("TFAP2A", "TACSTD2"),
            ("TFAP2A", "CLDN4"),
            ("TACSTD2", "CLDN4"),
        ]
        for left_symbol, right_symbol in pairs:
            left = [sample[left_symbol] for sample in complete.values()]
            right = [sample[right_symbol] for sample in complete.values()]
            correlation_rows.append(
                {
                    "cohort": cohort,
                    "gene_1": left_symbol,
                    "gene_2": right_symbol,
                    "n": len(left),
                    "spearman_rho": f"{pearson(ranks(left), ranks(right)):.4f}",
                    "pearson_log2_r_plus_1": f"{pearson([math.log2(x + 1) for x in left], [math.log2(x + 1) for x in right]):.4f}",
                }
            )

            left_cut = quantile(left, 0.75)
            right_cut = quantile(right, 0.75)
            left_high = {
                sample for sample, values in complete.items()
                if values[left_symbol] >= left_cut
            }
            right_high = {
                sample for sample, values in complete.items()
                if values[right_symbol] >= right_cut
            }
            intersection = left_high & right_high
            union = left_high | right_high
            overlap_rows.append(
                {
                    "cohort": cohort,
                    "gene_1": left_symbol,
                    "gene_2": right_symbol,
                    "gene_1_top_quartile_n": len(left_high),
                    "gene_2_top_quartile_n": len(right_high),
                    "intersection_n": len(intersection),
                    "jaccard": f"{len(intersection) / len(union):.4f}",
                }
            )

    write_tsv(
        "tcga_sample_expression.tsv",
        ["cohort", "sample_id", *GENES],
        expression_rows,
    )
    write_tsv(
        "tcga_gene_summary.tsv",
        ["cohort", "gene", "n", "median_rsem", "q1_rsem", "q3_rsem", "zero_count"],
        summary_rows,
    )
    write_tsv(
        "tcga_correlations.tsv",
        [
            "cohort",
            "gene_1",
            "gene_2",
            "n",
            "spearman_rho",
            "pearson_log2_r_plus_1",
        ],
        correlation_rows,
    )
    write_tsv(
        "tcga_top_quartile_overlap.tsv",
        [
            "cohort",
            "gene_1",
            "gene_2",
            "gene_1_top_quartile_n",
            "gene_2_top_quartile_n",
            "intersection_n",
            "jaccard",
        ],
        overlap_rows,
    )


if __name__ == "__main__":
    main()
