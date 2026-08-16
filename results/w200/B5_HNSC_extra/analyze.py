#!/usr/bin/env python3
"""B5: honest TACSTD2/CLDN4 audit in leftover public HNSCC ICI cohorts."""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import math
from pathlib import Path
from statistics import mean, median


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
MARKERS = ("TACSTD2", "CLDN4")

# Locked from Liu et al. 2021 Supplementary Table S2 (mmc1.pdf).
# The paper's primary binary definition combines Responder and Stable as
# clinical benefit and compares them with Progressor.
LIU_OUTCOME = {
    "HN001": "Responder",
    "HN002": "Progressor",
    "HN003": "Progressor",
    "HN004": "Progressor",
    "HN005": "Responder",
    "HN007": "Progressor",
    "HN008": "Progressor",
    "HN009": "Responder",
    "HN010": "Stable",
    "HN012": "Stable",
    "HN014": "Stable",
}


def read_gene_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes = {}
        for row in reader:
            genes[row[0].strip('"')] = [float(value) if value else math.nan for value in row[1:]]
    return [value.strip('"') for value in header[1:]], genes


def read_geo_series_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    samples: list[str] = []
    genes: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            row = next(csv.reader([line], delimiter="\t", quotechar='"'))
            if row[0] == "ID_REF":
                samples = row[1:]
            else:
                genes[row[0]] = [float(value) if value else math.nan for value in row[1:]]
    return samples, genes


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        rank = (start + 1 + stop) / 2
        for index in order[start:stop]:
            ranks[index] = rank
        start = stop
    return ranks


def exact_mann_whitney(
    values: list[float], positive_indices: tuple[int, ...]
) -> tuple[float, float, float]:
    ranks = average_ranks(values)
    n_positive = len(positive_indices)
    n_negative = len(values) - n_positive

    def u_stat(group: tuple[int, ...]) -> float:
        return sum(ranks[index] for index in group) - n_positive * (n_positive + 1) / 2

    observed = u_stat(positive_indices)
    center = n_positive * n_negative / 2
    null = [
        u_stat(group)
        for group in itertools.combinations(range(len(values)), n_positive)
    ]
    p_value = sum(abs(value - center) >= abs(observed - center) - 1e-12 for value in null) / len(null)
    return observed, observed / (n_positive * n_negative), p_value


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(order, 1))):
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def comparison(
    gene: str,
    values: list[float],
    labels: list[str],
    positive_labels: set[str],
    contrast: str,
) -> dict[str, object]:
    positive = tuple(index for index, label in enumerate(labels) if label in positive_labels)
    negative = tuple(index for index in range(len(labels)) if index not in positive)
    u_stat, auc, p_value = exact_mann_whitney(values, positive)
    return {
        "dataset": "GSE179730",
        "contrast": contrast,
        "gene": gene,
        "unit": "log2(CPM+1)",
        "n_positive": len(positive),
        "n_negative": len(negative),
        "positive_median": median(values[index] for index in positive),
        "negative_median": median(values[index] for index in negative),
        "positive_mean": mean(values[index] for index in positive),
        "negative_mean": mean(values[index] for index in negative),
        "mann_whitney_u": u_stat,
        "auc_positive_gt_negative": auc,
        "rank_biserial": 2 * auc - 1,
        "exact_permutation_p": p_value,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t" if path.suffix == ".tsv" else ",",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    _, prat_genes = read_geo_series_matrix(DATA / "GSE93157_series_matrix.txt.gz")
    foy_samples, foy_genes = read_gene_matrix(DATA / "GSE159067_IHN_log2cpm_data.txt.gz")
    liu_samples, liu_genes = read_gene_matrix(DATA / "GSE179730_RNAseq-combinedCPM.txt.gz")

    # The deposited Liu file is linear CPM despite GEO metadata saying log2 CPM:
    # every complete sample column sums to approximately one million.
    column_sums = [sum(values[index] for values in liu_genes.values()) for index in range(len(liu_samples))]
    if not all(990_000 < value < 1_010_000 for value in column_sums):
        raise ValueError("GSE179730 matrix does not have expected linear-CPM column sums")

    pre_indices = [index for index, sample in enumerate(liu_samples) if sample.endswith(".Pre")]
    pre_samples = [liu_samples[index] for index in pre_indices]
    patient_ids = [sample.split(".")[0] for sample in pre_samples]
    if set(patient_ids) != set(LIU_OUTCOME):
        raise ValueError("GSE179730 pretreatment samples do not match locked outcome table")
    labels = [LIU_OUTCOME[patient] for patient in patient_ids]

    per_sample: list[dict[str, object]] = []
    marker_values: dict[str, list[float]] = {}
    for gene in MARKERS:
        if gene not in liu_genes:
            raise ValueError(f"{gene} missing from GSE179730")
        linear_cpm = [liu_genes[gene][index] for index in pre_indices]
        marker_values[gene] = [math.log2(value + 1) for value in linear_cpm]
        for sample, patient, outcome, cpm, transformed in zip(
            pre_samples, patient_ids, labels, linear_cpm, marker_values[gene]
        ):
            per_sample.append(
                {
                    "dataset": "GSE179730",
                    "sample": sample,
                    "patient": patient,
                    "timepoint": "pretreatment",
                    "paper_outcome": outcome,
                    "clinical_benefit": outcome in {"Responder", "Stable"},
                    "gene": gene,
                    "cpm": cpm,
                    "log2_cpm_plus_1": transformed,
                }
            )
    write_csv(HERE / "per_sample_expression.csv", per_sample)

    primary = [
        comparison(
            gene,
            marker_values[gene],
            labels,
            {"Responder", "Stable"},
            "clinical_benefit_vs_progression",
        )
        for gene in MARKERS
    ]
    q_values = benjamini_hochberg(
        [float(row["exact_permutation_p"]) for row in primary]
    )
    for row, q_value in zip(primary, q_values):
        row["bh_q_across_primary_markers"] = q_value

    sensitivity = [
        comparison(
            gene,
            marker_values[gene],
            labels,
            {"Responder"},
            "pathologic_responder_vs_stable_or_progression",
        )
        for gene in MARKERS
    ]
    for row in sensitivity:
        row["bh_q_across_primary_markers"] = ""
    write_csv(HERE / "marker_statistics.csv", primary + sensitivity)

    catalog = [
        {
            "accession": "GSE93157",
            "citation": "Prat et al. 2017; PMID 28487385",
            "design": "mixed-cancer anti-PD-1; HNSCC n=5",
            "assay": "NanoString PanCancer Immune 730",
            "TACSTD2": "absent",
            "CLDN4": "absent",
            "verdict": "UNUSABLE",
            "reason": "Neither target is measured; previously audited in lung ICI work.",
        },
        {
            "accession": "GSE159067",
            "citation": "Foy et al. 2022; CLB-IHN",
            "design": f"R/M HNSCC PD-1/PD-L1; n={len(foy_samples)}",
            "assay": "HTG EdgeSeq targeted panel",
            "TACSTD2": "present" if "TACSTD2" in foy_genes else "absent",
            "CLDN4": "present" if "CLDN4" in foy_genes else "absent",
            "verdict": "UNUSABLE",
            "reason": "Neither target is measured on the deposited targeted panel.",
        },
        {
            "accession": "GSE179730",
            "citation": "Liu et al. 2021; PMID 34755131",
            "design": "neoadjuvant nivolumab OCSCC; 11 pretreatment RNA-seq tumors",
            "assay": "bulk RNA-seq, deposited linear CPM",
            "TACSTD2": "present",
            "CLDN4": "present",
            "verdict": "USABLE_UNDERPOWERED",
            "reason": "Full-transcriptome leftover; analyzed without model fitting.",
        },
        {
            "accession": "GSE190575",
            "citation": "ALPHA study",
            "design": "pembrolizumab plus afatinib R/M HNSCC",
            "assay": "NanoString targeted RCC",
            "TACSTD2": "not audited here",
            "CLDN4": "not audited here",
            "verdict": "EXCLUDED_COMBINATION",
            "reason": "EGFR inhibitor combination prevents an ICI-only interpretation.",
        },
        {
            "accession": "EGAD50000002506",
            "citation": "EGAS50000001746",
            "design": "R/M HNSCC PD-1/PD-L1",
            "assay": "bulk RNA-seq FASTQ",
            "TACSTD2": "expected",
            "CLDN4": "expected",
            "verdict": "CONTROLLED_ACCESS",
            "reason": "EGA data require DAC approval; not openly downloadable.",
        },
    ]
    write_csv(HERE / "cohort_catalog.tsv", catalog)

    summary = {
        "task": "B5_HNSC_extra",
        "question": "Pretreatment TACSTD2 and CLDN4 expression versus HNSCC ICI benefit in public leftover cohorts",
        "primary_dataset": "GSE179730",
        "primary_definition": "Liu paper definition: Responder + Stable (clinical benefit) versus Progressor",
        "counts": {
            "pretreatment_tumors": len(pre_samples),
            "clinical_benefit": sum(label in {"Responder", "Stable"} for label in labels),
            "progressor": sum(label == "Progressor" for label in labels),
            "pathologic_responder": sum(label == "Responder" for label in labels),
            "stable": sum(label == "Stable" for label in labels),
            "prat_hnsc_patients": 5,
            "foy_patients": len(foy_samples),
        },
        "primary_results": primary,
        "sensitivity_results": sensitivity,
        "marker_detection": {
            gene: {
                "nonzero_pretreatment_n": sum(value > 0 for value in marker_values[gene]),
                "total_n": len(pre_samples),
            }
            for gene in MARKERS
        },
        "data_checks": {
            "GSE179730_column_sum_min": min(column_sums),
            "GSE179730_column_sum_max": max(column_sums),
            "GSE93157_TACSTD2_present": "TACSTD2" in prat_genes,
            "GSE93157_CLDN4_present": "CLDN4" in prat_genes,
            "GSE159067_TACSTD2_present": "TACSTD2" in foy_genes,
            "GSE159067_CLDN4_present": "CLDN4" in foy_genes,
        },
        "honest_verdict": {
            "label": "NULL_AND_UNDERPOWERED",
            "did_we_tune_filters_or_thresholds": False,
            "statement": (
                "Prat and Foy cannot test these markers because their targeted panels omit both genes. "
                "In the only analyzed full-transcriptome leftover, GSE179730 (n=11), neither marker "
                "is associated with the authors' clinical-benefit definition after exact testing; "
                "expression is sparse and the cohort is too small for validation."
            ),
        },
    }
    (HERE / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
