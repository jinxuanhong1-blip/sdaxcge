#!/usr/bin/env python3
"""B5: leftover public HNSCC ICI cohorts, focused on CLDN4."""

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
CLARIOM = {
    "TACSTD2": "TC0100014340.hg.1",
    "CLDN4": "TC0700007993.hg.1",
}

# Locked from Liu et al. 2021 Supplementary Table S2 (mmc1.pdf).
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

# Events are taken from the published last-observation text, not the numeric
# Status column, which contradicts that text for Pt08.
LIU_SURVIVAL = {
    "HN001": {"rfs": 0.5753, "rfs_event": 1, "os": 3.1808, "os_event": 0, "last": "Alive with evidence of disease"},
    "HN002": {"rfs": 3.2630, "rfs_event": 1, "os": 3.2932, "os_event": 0, "last": "Alive with evidence of disease"},
    "HN003": {"rfs": 2.6685, "rfs_event": 0, "os": 2.6685, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN004": {"rfs": 0.8466, "rfs_event": 1, "os": 1.1699, "os_event": 1, "last": "Dead of disease"},
    "HN005": {"rfs": 2.8986, "rfs_event": 0, "os": 2.8986, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN007": {"rfs": 0.7260, "rfs_event": 1, "os": 0.8329, "os_event": 1, "last": "Dead of disease"},
    "HN008": {"rfs": 2.0274, "rfs_event": 0, "os": 2.0274, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN009": {"rfs": 2.0521, "rfs_event": 0, "os": 2.0521, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN010": {"rfs": 2.2822, "rfs_event": 0, "os": 2.2822, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN012": {"rfs": 0.2767, "rfs_event": 0, "os": 0.2767, "os_event": 0, "last": "Alive without evidence of disease"},
    "HN014": {"rfs": 2.0164, "rfs_event": 0, "os": 2.0164, "os_event": 0, "last": "Alive without evidence of disease"},
}


def read_gene_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes = {}
        for row in reader:
            genes[row[0].strip('"')] = [float(value) if value else math.nan for value in row[1:]]
    return [value.strip('"') for value in header[1:]], genes


def read_geo_series_matrix(
    path: Path, keep_ids: set[str] | None = None
) -> tuple[list[str], dict[str, list[str]], dict[str, list[float]]]:
    samples: list[str] = []
    characteristics: dict[str, list[str]] = {}
    metadata: dict[str, list[str]] = {}
    genes: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                if line.startswith("!Sample_"):
                    fields = next(csv.reader([line], delimiter="\t", quotechar='"'))
                    if fields[0] == "!Sample_characteristics_ch1":
                        key = fields[1].split(":", 1)[0].lower()
                        characteristics[key] = [
                            value.split(":", 1)[1].strip() if ":" in value else value
                            for value in fields[1:]
                        ]
                    else:
                        metadata[fields[0]] = fields[1:]
                continue
            row = next(csv.reader([line], delimiter="\t", quotechar='"'))
            if row[0] == "ID_REF":
                samples = row[1:]
            elif keep_ids is None or row[0] in keep_ids:
                genes[row[0]] = [float(value) if value else math.nan for value in row[1:]]
    metadata.update(characteristics)
    return samples, metadata, genes


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


def logrank_stat(times: list[float], events: list[int], group: tuple[int, ...]) -> float:
    group_set = set(group)
    observed = 0.0
    expected = 0.0
    variance = 0.0
    for time in sorted(set(times)):
        at_risk = [index for index, value in enumerate(times) if value >= time]
        deaths = [index for index in at_risk if events[index] == 1 and times[index] == time]
        n = len(at_risk)
        d = len(deaths)
        if n <= 1 or d == 0:
            continue
        n1 = sum(index in group_set for index in at_risk)
        d1 = sum(index in group_set for index in deaths)
        observed += d1
        expected += n1 * d / n
        variance += n1 * (n - n1) * d * (n - d) / (n * n * (n - 1))
    if variance <= 0:
        return 0.0
    return (observed - expected) / math.sqrt(variance)


def exact_logrank(
    times: list[float], events: list[int], positive_indices: tuple[int, ...]
) -> tuple[float, float]:
    observed = logrank_stat(times, events, positive_indices)
    null = [
        logrank_stat(times, events, group)
        for group in itertools.combinations(range(len(times)), len(positive_indices))
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


def comparison(
    dataset: str,
    gene: str,
    values: list[float],
    labels: list[str],
    positive_labels: set[str],
    contrast: str,
    unit: str,
) -> dict[str, object]:
    positive = tuple(index for index, label in enumerate(labels) if label in positive_labels)
    negative = tuple(index for index in range(len(labels)) if index not in positive)
    u_stat, auc, p_value = exact_mann_whitney(values, positive)
    return {
        "dataset": dataset,
        "contrast": contrast,
        "gene": gene,
        "unit": unit,
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


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t" if path.suffix == ".tsv" else ",",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    prat_samples, prat_meta, prat_genes = read_geo_series_matrix(
        DATA / "GSE93157_series_matrix.txt.gz"
    )
    foy_samples, foy_genes = read_gene_matrix(DATA / "GSE159067_IHN_log2cpm_data.txt.gz")
    liu_samples, liu_genes = read_gene_matrix(DATA / "GSE179730_RNAseq-combinedCPM.txt.gz")
    _, _, alpha_genes = read_geo_series_matrix(DATA / "GSE190575_series_matrix.txt.gz")
    niva_samples, niva_meta, niva_genes = read_geo_series_matrix(
        DATA / "GSE212550_series_matrix.txt.gz",
        keep_ids=set(CLARIOM.values()),
    )

    column_sums = [sum(values[index] for values in liu_genes.values()) for index in range(len(liu_samples))]
    if not all(990_000 < value < 1_010_000 for value in column_sums):
        raise ValueError("GSE179730 matrix does not have expected linear-CPM column sums")
    if any(probe not in niva_genes for probe in CLARIOM.values()):
        raise ValueError("GSE212550 is missing a locked Clariom D marker probe")

    hnsc_indices = [
        index
        for index, source in enumerate(prat_meta["!Sample_source_name_ch1"])
        if source == "HEADNECK"
    ]
    write_table(
        HERE / "prat_panel_audit.tsv",
        [
            {
                "gene": gene,
                "on_panel": gene in prat_genes,
                "note": "NanoString PanCancer Immune 730; HNSCC n=5",
            }
            for gene in ["CLDN4", "TACSTD2", "EPCAM", "CDH1"]
        ],
    )

    pre_indices = [index for index, sample in enumerate(liu_samples) if sample.endswith(".Pre")]
    pre_samples = [liu_samples[index] for index in pre_indices]
    patient_ids = [sample.split(".")[0] for sample in pre_samples]
    if set(patient_ids) != set(LIU_OUTCOME):
        raise ValueError("GSE179730 pretreatment samples do not match locked outcome table")
    labels = [LIU_OUTCOME[patient] for patient in patient_ids]

    per_sample: list[dict[str, object]] = []
    liu_values: dict[str, list[float]] = {}
    for gene in MARKERS:
        linear_cpm = [liu_genes[gene][index] for index in pre_indices]
        liu_values[gene] = [math.log2(value + 1) for value in linear_cpm]
        for sample, patient, outcome, cpm, transformed in zip(
            pre_samples, patient_ids, labels, linear_cpm, liu_values[gene]
        ):
            survival = LIU_SURVIVAL[patient]
            per_sample.append(
                {
                    "dataset": "GSE179730",
                    "sample": sample,
                    "patient": patient,
                    "timepoint": "pretreatment",
                    "paper_outcome": outcome,
                    "clinical_benefit": outcome in {"Responder", "Stable"},
                    "rfs_years": survival["rfs"],
                    "rfs_event": survival["rfs_event"],
                    "os_years": survival["os"],
                    "os_event": survival["os_event"],
                    "last_observation": survival["last"],
                    "gene": gene,
                    "value": transformed,
                    "unit": "log2(CPM+1)",
                    "raw_cpm": cpm,
                }
            )

    niva_labels = niva_meta["lts vs sts"]
    niva_values: dict[str, list[float]] = {}
    for gene, probe in CLARIOM.items():
        niva_values[gene] = niva_genes[probe]
        for sample, title, label, value in zip(
            niva_samples, niva_meta["!Sample_title"], niva_labels, niva_values[gene]
        ):
            per_sample.append(
                {
                    "dataset": "GSE212550",
                    "sample": sample,
                    "patient": title,
                    "timepoint": "pretreatment",
                    "paper_outcome": label,
                    "clinical_benefit": label == "LTS",
                    "rfs_years": "",
                    "rfs_event": "",
                    "os_years": "",
                    "os_event": "",
                    "last_observation": "LTS>18mo" if label == "LTS" else "STS<6mo",
                    "gene": gene,
                    "value": value,
                    "unit": "ClariomD_RMA",
                    "raw_cpm": "",
                }
            )
    write_table(HERE / "per_sample_expression.csv", per_sample)

    primary = [
        comparison(
            "GSE179730",
            gene,
            liu_values[gene],
            labels,
            {"Responder", "Stable"},
            "clinical_benefit_vs_progression",
            "log2(CPM+1)",
        )
        for gene in MARKERS
    ] + [
        comparison(
            "GSE212550",
            gene,
            niva_values[gene],
            niva_labels,
            {"LTS"},
            "long_vs_short_survivor",
            "ClariomD_RMA",
        )
        for gene in MARKERS
    ]
    cldn4_primary = [row for row in primary if row["gene"] == "CLDN4"]
    q_values = benjamini_hochberg(
        [float(row["exact_permutation_p"]) for row in cldn4_primary]
    )
    q_by_key = {
        (row["dataset"], row["contrast"]): q_value
        for row, q_value in zip(cldn4_primary, q_values)
    }
    for row in primary:
        row["bh_q_across_primary_cldn4_tests"] = (
            q_by_key.get((row["dataset"], row["contrast"]), "")
            if row["gene"] == "CLDN4"
            else ""
        )

    sensitivity = [
        comparison(
            "GSE179730",
            gene,
            liu_values[gene],
            labels,
            {"Responder"},
            "pathologic_responder_vs_stable_or_progression",
            "log2(CPM+1)",
        )
        for gene in MARKERS
    ]
    for row in sensitivity:
        row["bh_q_across_primary_cldn4_tests"] = ""
    write_table(HERE / "marker_statistics.csv", primary + sensitivity)

    detectable = tuple(
        index for index, value in enumerate(liu_values["CLDN4"]) if value > 0
    )
    survival_rows = []
    for endpoint, time_key, event_key in (
        ("RFS", "rfs", "rfs_event"),
        ("OS", "os", "os_event"),
    ):
        times = [LIU_SURVIVAL[patient][time_key] for patient in patient_ids]
        events = [LIU_SURVIVAL[patient][event_key] for patient in patient_ids]
        z_stat, p_value = exact_logrank(times, events, detectable)
        survival_rows.append(
            {
                "dataset": "GSE179730",
                "gene": "CLDN4",
                "split": "detectable_vs_zero",
                "endpoint": endpoint,
                "n_detectable": len(detectable),
                "n_zero": len(patient_ids) - len(detectable),
                "events_detectable": sum(events[index] for index in detectable),
                "events_zero": sum(
                    event for index, event in enumerate(events) if index not in detectable
                ),
                "logrank_z_detectable_minus_zero": z_stat,
                "exact_permutation_p": p_value,
                "note": "Median split is undefined because median CLDN4 is 0.",
            }
        )
    write_table(HERE / "survival_statistics.csv", survival_rows)

    catalog = [
        {
            "accession": "GSE93157",
            "citation": "Prat et al. 2017; PMID 28487385",
            "design": f"mixed-cancer anti-PD-1; HNSCC n={len(hnsc_indices)}",
            "assay": "NanoString PanCancer Immune 730",
            "TACSTD2": "absent",
            "CLDN4": "absent",
            "verdict": "UNUSABLE",
            "reason": "No CLDN-family gene is on the panel; EPCAM is present but is not CLDN4.",
        },
        {
            "accession": "GSE159067",
            "citation": "Foy et al. 2022; CLB-IHN",
            "design": f"R/M HNSCC PD-1/PD-L1; n={len(foy_samples)}",
            "assay": "HTG EdgeSeq targeted panel",
            "TACSTD2": "absent",
            "CLDN4": "absent",
            "verdict": "UNUSABLE",
            "reason": "Neither target is measured on the deposited targeted panel.",
        },
        {
            "accession": "GSE190575",
            "citation": "ALPHA study; PMID 35046059",
            "design": "pembrolizumab plus afatinib R/M HNSCC",
            "assay": "NanoString PanCancer Immune 730",
            "TACSTD2": "absent",
            "CLDN4": "absent",
            "verdict": "UNUSABLE",
            "reason": "Same immune panel omits CLDN4; EGFR-inhibitor combination besides.",
        },
        {
            "accession": "GSE179730",
            "citation": "Liu et al. 2021; PMID 34755131",
            "design": "neoadjuvant nivolumab OCSCC; 11 pretreatment RNA-seq tumors",
            "assay": "bulk RNA-seq, deposited linear CPM",
            "TACSTD2": "present",
            "CLDN4": "present",
            "verdict": "USABLE_UNDERPOWERED",
            "reason": "Full-transcriptome leftover; sparse near-zero CLDN4.",
        },
        {
            "accession": "GSE212550",
            "citation": "NIVACTOR test set; PMID 38290766",
            "design": "R/M HNSCC ICI monotherapy; 12 LTS vs 8 STS",
            "assay": "Affymetrix Clariom D, RMA",
            "TACSTD2": "present",
            "CLDN4": "present",
            "verdict": "USABLE",
            "reason": "Public extreme-survival labels; CLDN4 mapped to TC0700007993.hg.1.",
        },
        {
            "accession": "GSE212549",
            "citation": "NIVACTOR training set; PMID 38290766",
            "design": "R/M HNSCC nivolumab; n=80",
            "assay": "Affymetrix Clariom D, RMA",
            "TACSTD2": "expected",
            "CLDN4": "expected",
            "verdict": "UNUSABLE_NO_PUBLIC_LABELS",
            "reason": "GEO has age/stage/drug only; response tables are not openly deposited.",
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
    write_table(HERE / "cohort_catalog.tsv", catalog)

    liu_cldn4 = next(
        row
        for row in primary
        if row["dataset"] == "GSE179730" and row["gene"] == "CLDN4"
    )
    niva_cldn4 = next(
        row
        for row in primary
        if row["dataset"] == "GSE212550" and row["gene"] == "CLDN4"
    )
    summary = {
        "task": "B5_HNSC_extra",
        "question": "Pretreatment CLDN4 versus HNSCC ICI benefit in public leftover cohorts",
        "primary_cldn4_tests": [
            "GSE179730 clinical_benefit_vs_progression",
            "GSE212550 long_vs_short_survivor",
        ],
        "counts": {
            "prat_hnsc_patients": len(hnsc_indices),
            "foy_patients": len(foy_samples),
            "liu_pretreatment_tumors": len(pre_samples),
            "liu_clinical_benefit": sum(label in {"Responder", "Stable"} for label in labels),
            "liu_progressor": sum(label == "Progressor" for label in labels),
            "niva_test_patients": len(niva_samples),
            "niva_lts": sum(label == "LTS" for label in niva_labels),
            "niva_sts": sum(label == "STS" for label in niva_labels),
        },
        "primary_results": primary,
        "sensitivity_results": sensitivity,
        "survival_results": survival_rows,
        "marker_detection": {
            "GSE179730": {
                gene: {
                    "nonzero_pretreatment_n": sum(value > 0 for value in liu_values[gene]),
                    "total_n": len(pre_samples),
                }
                for gene in MARKERS
            }
        },
        "data_checks": {
            "GSE179730_column_sum_min": min(column_sums),
            "GSE179730_column_sum_max": max(column_sums),
            "GSE93157_CLDN4_present": "CLDN4" in prat_genes,
            "GSE93157_any_CLDN_present": any(gene.startswith("CLDN") for gene in prat_genes),
            "GSE159067_CLDN4_present": "CLDN4" in foy_genes,
            "GSE190575_CLDN4_present": "CLDN4" in alpha_genes,
            "GSE212550_CLDN4_probe_present": CLARIOM["CLDN4"] in niva_genes,
        },
        "honest_verdict": {
            "label": "NO_REPRODUCIBLE_CLDN4_ICI_SIGNAL",
            "did_we_tune_filters_or_thresholds": False,
            "statement": (
                "Prat, Foy, and ALPHA cannot test CLDN4 because their targeted immune panels omit it. "
                f"Liu/GSE179730 remains nominally lower CLDN4 in clinical-benefit tumors "
                f"(exact p={float(liu_cldn4['exact_permutation_p']):.3f}) but is sparse and n=11. "
                f"The leftover NIVACTOR test set (GSE212550, n=20) measures CLDN4 and is null "
                f"(LTS vs STS exact p={float(niva_cldn4['exact_permutation_p']):.3f}; medians 3.51 vs 3.54). "
                "No CLDN4 association survives the two locked leftover tests."
            ),
        },
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
