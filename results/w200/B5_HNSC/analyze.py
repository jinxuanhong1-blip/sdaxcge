#!/usr/bin/env python3
"""Audit public HNSCC ICI cohorts for a valid CLDN4-versus-ORR analysis.

This script deliberately stops before association testing unless expression and
patient-level RECIST response are both present and joinable in the same public
cohort. It uses only processed public files; no outcome labels are inferred.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    gse301_matrix = DATA / "candidate_GSE301741" / "GSE301741_series_matrix.txt.gz"

    required = [
        gse159_matrix,
        gse159_expr,
        gse931_matrix,
        gse931_expr,
        gse190_matrix,
        gse212_matrix,
        gse212_sqlite,
        gse301_matrix,
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
    cldn4_present190 = has_first_field(gse190_matrix, "CLDN4")
    response_fields190 = [
        values[0]
        for key, values in meta190
        if key == "!Sample_characteristics_ch1"
        and values
        and any(term in values[0].lower() for term in ("response", "recist", "best.resp"))
    ]

    meta212 = geo_metadata(gse212_matrix)
    response_fields212 = [
        values[0]
        for key, values in meta212
        if key == "!Sample_characteristics_ch1"
        and values
        and any(term in values[0].lower() for term in ("response", "recist", "best.resp"))
    ]
    with sqlite3.connect(gse212_sqlite) as connection:
        cldn4_probes = [
            row[0]
            for row in connection.execute(
                "SELECT probe_id FROM probes WHERE gene_id = ?", ("1364",)
            )
        ]
    cldn4_probe212 = cldn4_probes[0] if len(cldn4_probes) == 1 else ""
    cldn4_values212 = first_field_values(gse212_matrix, cldn4_probe212)

    meta301 = geo_metadata(gse301_matrix)
    response_fields301 = [
        values[0]
        for key, values in meta301
        if key == "!Sample_characteristics_ch1"
        and values
        and any(term in values[0].lower() for term in ("response", "recist", "best.resp"))
    ]

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
            "cohort": "NIVACTOR",
            "setting": "R/M HNSCC, nivolumab",
            "patients_or_samples": str(len(sample_field(meta212, "!Sample_geo_accession"))),
            "patient_level_orr": "no",
            "cldn4_measured": "yes",
            "expression_outcome_joinable": "no",
            "decision": "exclude",
            "reason": (
                f"CLDN4 is measured by {cldn4_probe212} in {len(cldn4_values212)} "
                "samples, but GEO has no patient-level response field. The paper "
                "reports aggregate GE-cohort counts (PR=12, SD=14, PD=53, one NA) "
                "without a sample-to-response key."
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
    ]

    output_table = ROOT / "candidate_audit.tsv"
    with output_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

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
            "response endpoint, inferred labels, or outcome-optimized expression cutoff was substituted."
        ),
        "verified_facts": {
            "GSE159067_response_counts": counts159,
            "GSE159067_CLDN4_row_present": has_first_field(gse159_expr, "CLDN4"),
            "GSE93157_HNSCC_n": len(hn_indices),
            "GSE93157_CLDN4_row_present": has_first_field(gse931_expr, "CLDN4"),
            "GSE190575_CLDN4_row_present": cldn4_present190,
            "GSE212549_CLDN4_probe": cldn4_probe212,
            "GSE212549_CLDN4_n": len(cldn4_values212),
            "GSE212549_CLDN4_range": [
                min(cldn4_values212),
                max(cldn4_values212),
            ],
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
    print(f"Wrote {output_table.relative_to(ROOT)} and analysis_status.json")


if __name__ == "__main__":
    main()
