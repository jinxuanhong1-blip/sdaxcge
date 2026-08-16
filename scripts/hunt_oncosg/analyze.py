#!/usr/bin/env python3
"""Purity-adjusted TACSTD2 versus IMSIG analysis in OncoSG LUAD."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy import stats


DATA_DIR = Path(os.environ.get("ONCOSG_DATA_DIR", "/tmp/hunt_oncosg_data"))
RESULTS_DIR = Path(
    os.environ.get("ONCOSG_RESULTS_DIR", "results/hunt_oncosg")
)
EXPRESSION = DATA_DIR / "expression_zscores.tsv"
CLINICAL = DATA_DIR / "clinical_sample.tsv"

METRICS = [
    ("IMSIG_B_CELLS", "B cells", True),
    ("IMSIG_INTERFERON", "interferon", True),
    ("IMSIG_MACROPHAGES", "macrophages", True),
    ("IMSIG_MONOCYTES", "monocytes", True),
    ("IMSIG_NEUTROPHILS", "neutrophils", True),
    ("IMSIG_NK_CELLS", "NK cells", True),
    ("IMSIG_PLASMA_CELLS", "plasma cells", True),
    ("IMSIG_T_CELLS", "T cells", True),
    ("IMSIG_PROLIFERATION", "proliferation", False),
    ("IMSIG_TRANSLATION", "translation", False),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def number(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA", "NaN"}:
        return None
    return float(value)


def read_tacstd2() -> tuple[dict[str, float], int, int]:
    with EXPRESSION.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_ids = header[2:]
        matches = [row for row in reader if row and row[0] == "TACSTD2"]
    if len(matches) != 1:
        raise ValueError(f"Expected one TACSTD2 row, found {len(matches)}")
    values = matches[0][2:]
    if len(values) != len(sample_ids):
        raise ValueError("TACSTD2 row and expression header have different lengths")
    expression = {
        sample_id: value
        for sample_id, raw in zip(sample_ids, values, strict=True)
        if (value := number(raw)) is not None
    }
    return expression, len(sample_ids), len(header)


def read_clinical() -> dict[str, dict[str, str]]:
    with CLINICAL.open(encoding="utf-8", newline="") as handle:
        for _ in range(4):
            next(handle)
        rows = csv.DictReader(handle, delimiter="\t")
        clinical = {row["SAMPLE_ID"]: row for row in rows}
    if len(clinical) != len(set(clinical)):
        raise ValueError("Duplicate SAMPLE_ID values in clinical data")
    return clinical


def residualize(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(covariate)), covariate])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def partial_spearman(
    x: np.ndarray, y: np.ndarray, purity: np.ndarray
) -> tuple[float, float, float, float]:
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    purity_rank = stats.rankdata(purity)
    rho = float(
        stats.pearsonr(
            residualize(x_rank, purity_rank),
            residualize(y_rank, purity_rank),
        ).statistic
    )
    degrees_freedom = len(x) - 3
    if abs(rho) == 1:
        p_value = 0.0
    else:
        statistic = rho * math.sqrt(degrees_freedom / (1 - rho * rho))
        p_value = float(2 * stats.t.sf(abs(statistic), degrees_freedom))
    # Fisher-z approximation for one controlled covariate.
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    half_width = stats.norm.ppf(0.975) / math.sqrt(len(x) - 4)
    low, high = np.tanh([z - half_width, z + half_width])
    return rho, p_value, float(low), float(high)


def bh_adjust(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = np.empty(count)
    running = 1.0
    for position in range(count - 1, -1, -1):
        index = int(order[position])
        candidate = p_values[index] * count / (position + 1)
        running = min(running, candidate)
        adjusted[index] = running
    return adjusted.tolist()


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.8g}"


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    expression, expression_samples, expression_columns = read_tacstd2()
    clinical = read_clinical()

    sample_rows: list[dict[str, object]] = []
    for sample_id, tacstd2 in expression.items():
        metadata = clinical.get(sample_id)
        if metadata is None:
            continue
        row: dict[str, object] = {
            "sample_id": sample_id,
            "TACSTD2_expression_zscore": tacstd2,
            "purity": number(metadata.get("PURITY")),
        }
        for metric, _, _ in METRICS:
            row[metric] = number(metadata.get(metric))
        sample_rows.append(row)

    results: list[dict[str, object]] = []
    for metric, label, primary in METRICS:
        complete = [
            row
            for row in sample_rows
            if row["purity"] is not None and row[metric] is not None
        ]
        x = np.asarray(
            [row["TACSTD2_expression_zscore"] for row in complete], dtype=float
        )
        y = np.asarray([row[metric] for row in complete], dtype=float)
        purity = np.asarray([row["purity"] for row in complete], dtype=float)
        if len(x) < 5:
            raise ValueError(f"Insufficient complete cases for {metric}: {len(x)}")
        marginal = stats.spearmanr(x, y)
        partial, partial_p, ci_low, ci_high = partial_spearman(x, y, purity)
        results.append(
            {
                "metric": metric,
                "label": label,
                "metric_class": "immune_primary" if primary else "nonimmune_control",
                "n_complete": len(x),
                "marginal_spearman_rho": float(marginal.statistic),
                "marginal_p_value": float(marginal.pvalue),
                "purity_partial_spearman_rho": partial,
                "partial_95ci_low": ci_low,
                "partial_95ci_high": ci_high,
                "partial_p_value": partial_p,
                "partial_bh_q_value_immune_family": None,
            }
        )

    primary_indices = [
        index
        for index, result in enumerate(results)
        if result["metric_class"] == "immune_primary"
    ]
    adjusted = bh_adjust(
        [float(results[index]["partial_p_value"]) for index in primary_indices]
    )
    for index, q_value in zip(primary_indices, adjusted, strict=True):
        results[index]["partial_bh_q_value_immune_family"] = q_value

    result_fields = list(results[0])
    with (RESULTS_DIR / "tacstd2_immune_partial.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=result_fields)
        writer.writeheader()
        for result in results:
            writer.writerow({key: fmt(value) if isinstance(value, float) else value
                             for key, value in result.items()})

    sample_fields = [
        "sample_id",
        "TACSTD2_expression_zscore",
        "purity",
        *[metric for metric, _, _ in METRICS],
    ]
    with (RESULTS_DIR / "analysis_samples.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields)
        writer.writeheader()
        for row in sample_rows:
            writer.writerow(
                {
                    key: fmt(value) if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )

    audit = {
        "study_id": "luad_oncosg_2020",
        "expression_profile": (
            "luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores"
        ),
        "expression_header_sample_count": expression_samples,
        "expression_total_column_count": expression_columns,
        "TACSTD2_nonmissing_count": len(expression),
        "clinical_sample_count": len(clinical),
        "matched_expression_clinical_count": len(sample_rows),
        "matched_with_nonmissing_purity_count": sum(
            row["purity"] is not None for row in sample_rows
        ),
        "primary_immune_hypothesis_count": len(primary_indices),
        "nonimmune_control_count": len(results) - len(primary_indices),
        "expression_sha256": sha256(EXPRESSION),
        "clinical_sha256": sha256(CLINICAL),
        "software": {
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
        },
    }
    with (RESULTS_DIR / "data_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")

    with (RESULTS_DIR / "source_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["artifact", "source_url", "sha256", "bytes"])
        for artifact, source_url in [
            (
                EXPRESSION,
                "https://media.githubusercontent.com/media/cBioPortal/datahub/"
                "master/public/luad_oncosg_2020/"
                "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
            ),
            (
                CLINICAL,
                "https://media.githubusercontent.com/media/cBioPortal/datahub/"
                "master/public/luad_oncosg_2020/data_clinical_sample.txt",
            ),
        ]:
            writer.writerow(
                [artifact.name, source_url, sha256(artifact), artifact.stat().st_size]
            )

    immune_results = [
        result for result in results if result["metric_class"] == "immune_primary"
    ]
    significant = [
        result
        for result in immune_results
        if float(result["partial_bh_q_value_immune_family"]) < 0.05
    ]
    strongest = max(
        immune_results, key=lambda result: abs(float(result["purity_partial_spearman_rho"]))
    )
    summary = {
        "n_significant_immune_metrics_bh_0_05": len(significant),
        "strongest_absolute_partial_association": {
            "metric": strongest["metric"],
            "rho": strongest["purity_partial_spearman_rho"],
            "p_value": strongest["partial_p_value"],
            "q_value": strongest["partial_bh_q_value_immune_family"],
        },
    }
    with (RESULTS_DIR / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
