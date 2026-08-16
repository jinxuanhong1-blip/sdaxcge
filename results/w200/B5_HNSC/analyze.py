#!/usr/bin/env python3
"""Audit public HNSCC ICI cohorts for a valid CLDN4-versus-ORR analysis.

Primary question (frozen): pretreatment tumor CLDN4 versus RECIST ORR
(CR/PR vs SD/PD). This script runs association tests only when expression and
patient-level RECIST labels are both present and joinable. No RECIST cohort
meets that bar.

A leftover neoadjuvant cohort (GSE179730) has joinable CLDN4 and Table S2
pathologic/hybrid outcomes. Those leftover tests are labeled as not RECIST ORR
and do not replace the primary result.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
import sqlite3
from pathlib import Path
from statistics import mean, median


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
CLDN4_GENE = "CLDN4"
CLDN4_ENTREZ = "1364"
CLDN4_PROBE_CLARIOM = "TC0700007993.hg.1"
RECIST_PR_CUTOFF = -30.0
RECIST_PD_CUTOFF = 20.0


def geo_metadata(path: Path) -> list[tuple[str, list[str]]]:
    rows: list[tuple[str, list[str]]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                row = next(csv.reader([line], delimiter="\t"))
                rows.append((row[0], [value.strip('"') for value in row[1:]]))
            if line.startswith("!series_matrix_table_begin"):
                break
    return rows


def characteristic(metadata: list[tuple[str, list[str]]], prefix: str) -> list[str]:
    prefix_lower = prefix.lower() + ":"
    for key, values in metadata:
        if key != "!Sample_characteristics_ch1" or not values:
            continue
        if values[0].lower().startswith(prefix_lower):
            return [value.split(":", 1)[1].strip() for value in values]
    return []


def sample_field(metadata: list[tuple[str, list[str]]], key: str) -> list[str]:
    for row_key, values in metadata:
        if row_key == key:
            return values
    return []


def response_like_fields(metadata: list[tuple[str, list[str]]]) -> list[str]:
    hits: list[str] = []
    for key, values in metadata:
        if key != "!Sample_characteristics_ch1" or not values:
            continue
        if any(term in values[0].lower() for term in ("response", "recist", "best.resp", "orr")):
            hits.append(values[0])
    return hits


def has_first_field(path: Path, target: str) -> bool:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            first = line.split("\t", 1)[0].strip('"')
            if first == target:
                return True
    return False


def first_field_values(path: Path, target: str) -> list[float]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            row = next(csv.reader([line], delimiter="\t"))
            if row and row[0].strip('"') == target:
                return [float(value) for value in row[1:] if value != ""]
    return []


def read_gene_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes = {
            row[0].strip('"'): [float(value) if value else math.nan for value in row[1:]]
            for row in reader
            if row
        }
    return [value.strip('"') for value in header[1:]], genes


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_liu_table(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    table = {row["patient_geo"]: row for row in rows}
    if len(table) != 12:
        raise ValueError(f"Locked Table S2 must contain 12 patients, found {len(table)}")
    return table


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
    p_value = sum(
        abs(value - center) >= abs(observed - center) - 1e-12 for value in null
    ) / len(null)
    return observed, observed / (n_positive * n_negative), p_value


def leftover_comparison(
    values: list[float],
    labels: list[str],
    positive_labels: set[str],
    contrast: str,
    definition: str,
    is_leftover_primary: bool,
    is_recist_orr: bool,
) -> dict[str, object]:
    positive = tuple(index for index, label in enumerate(labels) if label in positive_labels)
    negative = tuple(index for index in range(len(labels)) if index not in set(positive))
    u_stat, auc, p_value = exact_mann_whitney(values, positive)
    return {
        "dataset": "GSE179730",
        "gene": CLDN4_GENE,
        "unit": "log2(CPM+1)",
        "contrast": contrast,
        "definition": definition,
        "is_b5_primary": "no",
        "is_leftover_primary": "yes" if is_leftover_primary else "no",
        "is_recist_orr": "yes" if is_recist_orr else "no",
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
        "n_permutations": math.comb(len(values), len(positive)),
    }


def write_delimited(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    gse159_matrix = DATA / "raw" / "GSE159067_series_matrix.txt.gz"
    gse159_expr = DATA / "raw" / "GSE159067_IHN_log2cpm_data.txt.gz"
    gse931_matrix = DATA / "candidate_GSE93157" / "GSE93157_series_matrix.txt.gz"
    gse931_expr = DATA / "candidate_GSE93157" / "GSE93157_raw_data_values.txt.gz"
    gse190_matrix = DATA / "candidate_GSE190575" / "GSE190575_series_matrix.txt.gz"
    gse212_matrix = DATA / "candidate_GSE212549" / "GSE212549_series_matrix.txt.gz"
    gse212_sqlite = (
        DATA
        / "candidate_GSE212549"
        / "clariomdhumantranscriptcluster.db"
        / "inst"
        / "extdata"
        / "clariomdhumantranscriptcluster.sqlite"
    )
    gse212550_matrix = DATA / "candidate_GSE212550" / "GSE212550_series_matrix.txt.gz"
    gse301_matrix = DATA / "candidate_GSE301741" / "GSE301741_series_matrix.txt.gz"
    gse179_matrix = DATA / "candidate_GSE179730" / "GSE179730_series_matrix.txt.gz"
    gse179_expr = DATA / "candidate_GSE179730" / "GSE179730_RNAseq-combinedCPM.txt.gz"
    liu_table_path = ROOT / "liu_table_s2.tsv"

    required = [
        gse159_matrix,
        gse159_expr,
        gse931_matrix,
        gse931_expr,
        gse190_matrix,
        gse212_matrix,
        gse212_sqlite,
        gse212550_matrix,
        gse301_matrix,
        gse179_matrix,
        gse179_expr,
        liu_table_path,
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required audit inputs: {missing}")

    meta159 = geo_metadata(gse159_matrix)
    response159 = characteristic(meta159, "best response on immunotherapy (recist)")
    counts159 = {label: response159.count(label) for label in ("CR", "PR", "SD", "PD")}

    meta931 = geo_metadata(gse931_matrix)
    source931 = sample_field(meta931, "!Sample_source_name_ch1")
    response931 = characteristic(meta931, "best.resp")
    hn_indices = [i for i, source in enumerate(source931) if source == "HEADNECK"]
    hn_response931 = [response931[i] for i in hn_indices]

    meta190 = geo_metadata(gse190_matrix)
    cldn4_present190 = has_first_field(gse190_matrix, CLDN4_GENE)
    response_fields190 = response_like_fields(meta190)

    meta212 = geo_metadata(gse212_matrix)
    response_fields212 = response_like_fields(meta212)
    with sqlite3.connect(gse212_sqlite) as connection:
        cldn4_probes = [
            row[0]
            for row in connection.execute(
                "SELECT probe_id FROM probes WHERE gene_id = ?", (CLDN4_ENTREZ,)
            )
        ]
    if cldn4_probes != [CLDN4_PROBE_CLARIOM]:
        raise ValueError(f"Unexpected Clariom D CLDN4 mapping: {cldn4_probes}")
    cldn4_values212 = first_field_values(gse212_matrix, CLDN4_PROBE_CLARIOM)

    meta212550 = geo_metadata(gse212550_matrix)
    response_fields212550 = response_like_fields(meta212550)
    lts_sts = characteristic(meta212550, "lts vs sts")
    lts_sts_counts = {label: lts_sts.count(label) for label in sorted(set(lts_sts))}
    cldn4_values212550 = first_field_values(gse212550_matrix, CLDN4_PROBE_CLARIOM)

    meta301 = geo_metadata(gse301_matrix)
    response_fields301 = response_like_fields(meta301)

    meta179 = geo_metadata(gse179_matrix)
    response_fields179 = response_like_fields(meta179)
    liu_table = load_liu_table(liu_table_path)
    liu_samples, liu_genes = read_gene_matrix(gse179_expr)
    if CLDN4_GENE not in liu_genes:
        raise ValueError("CLDN4 missing from GSE179730")

    column_sums = [
        sum(values[index] for values in liu_genes.values())
        for index in range(len(liu_samples))
    ]
    if not all(990_000 < value < 1_010_000 for value in column_sums):
        raise ValueError("GSE179730 matrix does not have expected linear-CPM column sums")

    pre_indices = [index for index, sample in enumerate(liu_samples) if sample.endswith(".Pre")]
    pre_samples = [liu_samples[index] for index in pre_indices]
    patient_ids = [sample.split(".")[0] for sample in pre_samples]
    expected_rna = {
        patient for patient, row in liu_table.items() if row["pretreatment_rnaseq"] == "yes"
    }
    if set(patient_ids) != expected_rna:
        raise ValueError("GSE179730 pretreatment samples do not match locked Table S2 RNA flag")

    leftover_rows: list[dict[str, object]] = []
    leftover_values: list[float] = []
    leftover_labels: list[str] = []
    radiographic_pr = 0
    radiographic_pd = 0
    for sample, patient, index in zip(pre_samples, patient_ids, pre_indices):
        locked = liu_table[patient]
        cpm = liu_genes[CLDN4_GENE][index]
        transformed = math.log2(cpm + 1)
        radiographic = float(locked["radiographic_pct_scan1_vs_scan2"])
        if radiographic <= RECIST_PR_CUTOFF:
            radiographic_pr += 1
        if radiographic >= RECIST_PD_CUTOFF:
            radiographic_pd += 1
        leftover_values.append(transformed)
        leftover_labels.append(locked["table_s2_outcome"])
        leftover_rows.append(
            {
                "dataset": "GSE179730",
                "sample": sample,
                "patient_geo": patient,
                "patient_paper": locked["patient_paper"],
                "timepoint": "pretreatment",
                "table_s2_outcome": locked["table_s2_outcome"],
                "radiographic_pct_scan1_vs_scan2": locked["radiographic_pct_scan1_vs_scan2"],
                "pathologic_pct_scan1_vs_path": locked["pathologic_pct_scan1_vs_path"],
                "radiographic_recist_pr": "yes" if radiographic <= RECIST_PR_CUTOFF else "no",
                "leftover_orr_analog": "yes" if locked["table_s2_outcome"] == "Responder" else "no",
                "leftover_clinical_benefit_not_orr": (
                    "yes" if locked["table_s2_outcome"] in {"Responder", "Stable"} else "no"
                ),
                "cldn4_cpm": cpm,
                "cldn4_log2_cpm_plus_1": transformed,
            }
        )
    write_delimited(ROOT / "leftover_per_sample.tsv", leftover_rows)

    leftover_stats = [
        leftover_comparison(
            leftover_values,
            leftover_labels,
            {"Responder"},
            "leftover_pathologic_orr_analog",
            "Table S2 Responder versus Stable+Progressor (pathologic/hybrid ORR analog; not radiographic RECIST)",
            True,
            False,
        ),
        leftover_comparison(
            leftover_values,
            leftover_labels,
            {"Responder", "Stable"},
            "leftover_clinical_benefit_not_orr",
            "Table S2 Responder+Stable versus Progressor (clinical benefit / DCR analog; not ORR)",
            False,
            False,
        ),
    ]
    write_delimited(ROOT / "leftover_statistics.tsv", leftover_stats)

    audit_rows = [
        {
            "accession": "GSE159067",
            "cohort": "CLB-IHN",
            "setting": "R/M HNSCC, PD-1/PD-L1 therapy",
            "patients_or_samples": str(len(response159)),
            "patient_level_orr": "yes",
            "cldn4_measured": "no",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"RECIST labels public (CR={counts159['CR']}, PR={counts159['PR']}, "
                f"SD={counts159['SD']}, PD={counts159['PD']}), but CLDN4 is absent "
                "from the 2,559-transcript HTG panel."
            ),
        },
        {
            "accession": "GSE93157",
            "cohort": "Prat et al.",
            "setting": "Pan-cancer anti-PD-1; HNSCC subset",
            "patients_or_samples": str(len(hn_indices)),
            "patient_level_orr": "yes",
            "cldn4_measured": "no",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"HNSCC subset has {len(hn_indices)} cases and public best responses "
                f"({','.join(hn_response931)}), but CLDN4 is absent from the "
                "NanoString PanCancer Immune panel."
            ),
        },
        {
            "accession": "GSE190575",
            "cohort": "ALPHA",
            "setting": "R/M HNSCC, afatinib plus pembrolizumab",
            "patients_or_samples": str(len(sample_field(meta190, "!Sample_geo_accession"))),
            "patient_level_orr": "no",
            "cldn4_measured": "no",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"No response field in GEO metadata ({len(response_fields190)} found); "
                f"CLDN4 row present={cldn4_present190} in the deposited NanoString "
                "matrix; combination "
                "therapy is not an ICI-monotherapy analog."
            ),
        },
        {
            "accession": "GSE212549",
            "cohort": "NIVACTOR training",
            "setting": "R/M HNSCC, nivolumab",
            "patients_or_samples": str(len(sample_field(meta212, "!Sample_geo_accession"))),
            "patient_level_orr": "no",
            "cldn4_measured": "yes",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"CLDN4 is measured by {CLDN4_PROBE_CLARIOM} in {len(cldn4_values212)} "
                "samples, but GEO has no patient-level response field. The paper "
                "reports aggregate GE-cohort counts (PR=12, SD=14, PD=53, one NA) "
                "without a sample-to-response key."
            ),
        },
        {
            "accession": "GSE212550",
            "cohort": "NIVACTOR testing",
            "setting": "R/M HNSCC, ICI monotherapy; LTS vs STS survival",
            "patients_or_samples": str(len(sample_field(meta212550, "!Sample_geo_accession"))),
            "patient_level_orr": "no",
            "cldn4_measured": "yes",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"CLDN4 is measured by {CLDN4_PROBE_CLARIOM} in {len(cldn4_values212550)} "
                f"samples (range {min(cldn4_values212550):.2f}-{max(cldn4_values212550):.2f}). "
                f"GEO exposes only LTS vs STS survival ({', '.join(f'{k}={v}' for k, v in lts_sts_counts.items())}), "
                f"not RECIST; response-like fields={len(response_fields212550)}."
            ),
        },
        {
            "accession": "GSE301741",
            "cohort": "Mints et al.",
            "setting": "Neoadjuvant pembrolizumab, single-cell RNA-seq",
            "patients_or_samples": "16 patients / 58 libraries",
            "patient_level_orr": "no",
            "cldn4_measured": "yes",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"Endpoint is pathologic tumor regression rather than RECIST ORR, "
                f"and GEO metadata contains {len(response_fields301)} response fields. "
                "The publication identifies CLDN4 as an epithelial-state marker."
            ),
        },
        {
            "accession": "GSE179730",
            "cohort": "Liu / Knochelmann et al.",
            "setting": "Neoadjuvant nivolumab, resectable oral-cavity SCC",
            "patients_or_samples": f"{len(pre_samples)} pretreatment RNA-seq / 12 trial patients",
            "patient_level_orr": "no",
            "cldn4_measured": "yes",
            "expression_outcome_joinable": "leftover_pathologic_only",
            "decision": "leftover_not_recist",
            "reason": (
                f"CLDN4 is present. GEO has {len(response_fields179)} response fields. "
                "Table S2 locks 4 Responders / 3 Stable / 5 Progressors; Pt06 Responder "
                "has no pretreatment RNA (n=3 vs 8 leftover). Radiographic scan1-vs-scan2 "
                f"RECIST PR count among RNA cases={radiographic_pr}. Outcome is "
                "pathologic/hybrid size change, not radiographic RECIST ORR."
            ),
        },
    ]
    write_delimited(ROOT / "candidate_audit.tsv", audit_rows)

    leftover_orr = leftover_stats[0]
    leftover_dcr = leftover_stats[1]
    status = {
        "question": "Pretreatment tumor CLDN4 expression versus RECIST ORR in public HNSCC ICI data",
        "estimable": False,
        "association_tests_run": 0,
        "primary_result": "NOT ESTIMABLE FROM ELIGIBLE PUBLIC DATA",
        "why": (
            "No screened public cohort simultaneously provides a CLDN4 measurement "
            "and a patient-level RECIST response label that can be joined without inference."
        ),
        "anti_p_hacking_guardrail": (
            "No alternate gene, disease-control endpoint, survival endpoint, pathologic "
            "response endpoint, inferred labels, or outcome-optimized expression cutoff was "
            "substituted for the primary RECIST-ORR claim."
        ),
        "verified_facts": {
            "GSE159067_response_counts": counts159,
            "GSE159067_CLDN4_row_present": has_first_field(gse159_expr, CLDN4_GENE),
            "GSE93157_HNSCC_n": len(hn_indices),
            "GSE93157_CLDN4_row_present": has_first_field(gse931_expr, CLDN4_GENE),
            "GSE190575_CLDN4_row_present": cldn4_present190,
            "GSE212549_CLDN4_probe": CLDN4_PROBE_CLARIOM,
            "GSE212549_CLDN4_n": len(cldn4_values212),
            "GSE212549_CLDN4_range": [min(cldn4_values212), max(cldn4_values212)],
            "GSE212549_response_like_fields": len(response_fields212),
            "GSE212550_CLDN4_n": len(cldn4_values212550),
            "GSE212550_CLDN4_range": [min(cldn4_values212550), max(cldn4_values212550)],
            "GSE212550_LTS_STS_counts": lts_sts_counts,
            "GSE212550_response_like_fields": len(response_fields212550),
            "GSE179730_CLDN4_row_present": True,
            "GSE179730_GEO_response_fields": len(response_fields179),
            "GSE179730_column_sum_min": min(column_sums),
            "GSE179730_column_sum_max": max(column_sums),
            "GSE179730_pretreatment_n": len(pre_samples),
            "GSE179730_table_s2_counts_all_12": {
                "Responder": sum(row["table_s2_outcome"] == "Responder" for row in liu_table.values()),
                "Stable": sum(row["table_s2_outcome"] == "Stable" for row in liu_table.values()),
                "Progressor": sum(row["table_s2_outcome"] == "Progressor" for row in liu_table.values()),
            },
            "GSE179730_RNA_outcome_counts": {
                label: leftover_labels.count(label)
                for label in ("Responder", "Stable", "Progressor")
            },
            "GSE179730_CLDN4_nonzero_pretreatment_n": sum(row["cldn4_cpm"] > 0 for row in leftover_rows),
            "GSE179730_radiographic_RECIST_PR_among_RNA": radiographic_pr,
            "GSE179730_radiographic_RECIST_PD_among_RNA": radiographic_pd,
            "GSE179730_Pt06_responder_lacks_pretreatment_RNA": True,
        },
        "leftover": {
            "accession": "GSE179730",
            "role": "neoadjuvant leftover; not the primary RECIST-ORR analysis",
            "association_tests_run": 2,
            "gene": CLDN4_GENE,
            "substitute_gene_used": False,
            "label_source": "liu_table_s2.tsv independently verified from PMC8561238 mmc1.pdf Table S2",
            "pathologic_orr_analog": leftover_orr,
            "clinical_benefit_not_orr": leftover_dcr,
            "honest_leftover_verdict": (
                "Leftover pathologic-ORR analog is inconclusive "
                f"(n={leftover_orr['n_positive']} vs {leftover_orr['n_negative']}, "
                f"AUC={leftover_orr['auc_positive_gt_negative']}, "
                f"exact p={leftover_orr['exact_permutation_p']}). "
                "Radiographic RECIST ORR is not testable here (0 RNA-profiled cases "
                "meet -30% PR). The clinical-benefit contrast is DCR, not ORR, and "
                "is not the B5 primary claim."
            ),
        },
    }
    (ROOT / "analysis_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )

    checksums = sorted(
        (sha256(path), str(path.relative_to(ROOT)))
        for path in required
    )
    with (ROOT / "checksums.sha256").open("w", encoding="utf-8") as handle:
        for digest, relative_path in checksums:
            handle.write(f"{digest}  {relative_path}\n")

    print(status["primary_result"])
    print(status["leftover"]["honest_leftover_verdict"])
    print("Wrote candidate_audit.tsv, leftover tables, and analysis_status.json")


if __name__ == "__main__":
    main()
