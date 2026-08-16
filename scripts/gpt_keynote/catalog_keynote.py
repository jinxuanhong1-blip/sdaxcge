#!/usr/bin/env python3
"""Build an evidence-based catalog of processed KEYNOTE lung expression data.

The script uses only the Python standard library.  It does not download any
biological matrix unless a curated record satisfies all three conditions:
KEYNOTE trial specimen, lung cancer, and anonymously open processed data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "gpt_keynote"
MAX_BYTES = 2_000_000_000

GEO_SEARCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=gds&term=KEYNOTE%5BAll%20Fields%5D&retmode=json&retmax=100"
)
GEO_SUMMARY = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    "?db=gds&id={ids}&retmode=json"
)

FIELDS = [
    "record",
    "trial",
    "lung_scope",
    "assay_or_data",
    "processed_patient_level_matrix",
    "access",
    "downloaded",
    "tacstd2",
    "cldn4",
    "evidence",
    "honest_note",
]

CATALOG = [
    {
        "record": "KEYNOTE-001 NSCLC cohort",
        "trial": "KEYNOTE-001 / NCT01295827",
        "lung_scope": "yes",
        "assay_or_data": "PD-L1 IHC in the pivotal NSCLC report; genomic biomarker work exists",
        "processed_patient_level_matrix": "not found",
        "access": "controlled/on-request",
        "downloaded": "no",
        "tacstd2": "not assessable",
        "cldn4": "not assessable",
        "evidence": "https://doi.org/10.1056/NEJMoa1501824",
        "honest_note": "No anonymous processed RNA/count matrix located for the NSCLC cohort.",
    },
    {
        "record": "GSE216069 and GSE216055",
        "trial": "KEYNOTE-001 / NCT01295827",
        "lung_scope": "no (KEYNOTE specimen is melanoma)",
        "assay_or_data": "processed snRNA-seq and Slide-seq V2 counts",
        "processed_patient_level_matrix": "yes, but out of scope",
        "access": "open",
        "downloaded": "no",
        "tacstd2": "not tested; out of scope",
        "cldn4": "not tested; out of scope",
        "evidence": "https://doi.org/10.1038/s41588-022-01268-9",
        "honest_note": (
            "The deposit also contains separate NSCLC method-validation samples, "
            "but the three sequential KEYNOTE-001 biopsies are explicitly melanoma."
        ),
    },
    {
        "record": "KEYNOTE-010",
        "trial": "KEYNOTE-010 / NCT01905657",
        "lung_scope": "yes (NSCLC)",
        "assay_or_data": "clinical/IHC; samples used in sponsor biomarker analyses",
        "processed_patient_level_matrix": "not found",
        "access": "controlled/on-request",
        "downloaded": "no",
        "tacstd2": "not assessable",
        "cldn4": "not assessable",
        "evidence": "https://doi.org/10.1016/S0140-6736(15)01281-7",
        "honest_note": "No open RNA, NanoString, or count matrix located.",
    },
    {
        "record": "KEYNOTE-028 SCLC cohort",
        "trial": "KEYNOTE-028 / NCT02054806",
        "lung_scope": "yes (small-cell lung cohort)",
        "assay_or_data": "T-cell-inflamed GEP (NanoString-derived) and TMB analyses",
        "processed_patient_level_matrix": "not found",
        "access": "controlled/on-request; publication summaries open",
        "downloaded": "no",
        "tacstd2": "not in published 18-gene GEP",
        "cldn4": "not in published 18-gene GEP",
        "evidence": "https://doi.org/10.1200/JCO.2018.78.2276",
        "honest_note": "Patient-level normalized counts/GEP matrix was not deposited openly.",
    },
    {
        "record": "KEYNOTE-042",
        "trial": "KEYNOTE-042 / NCT02220894",
        "lung_scope": "yes (NSCLC)",
        "assay_or_data": "PD-L1 IHC and exploratory WES/TMB",
        "processed_patient_level_matrix": "not found (not an RNA study)",
        "access": "controlled/on-request",
        "downloaded": "no",
        "tacstd2": "not assessable",
        "cldn4": "not assessable",
        "evidence": "https://doi.org/10.1016/j.annonc.2023.01.011",
        "honest_note": "The reported molecular analysis is WES, not expression counts.",
    },
    {
        "record": "KEYNOTE-189 and KEYNOTE-407",
        "trial": "NCT02578680 and NCT02775435",
        "lung_scope": "yes (NSCLC)",
        "assay_or_data": "exploratory WES/TMB and selected mutation status",
        "processed_patient_level_matrix": "not found (reported analysis is not RNA)",
        "access": "controlled/on-request",
        "downloaded": "no",
        "tacstd2": "not assessable",
        "cldn4": "not assessable",
        "evidence": "https://doi.org/10.1016/j.jtocrr.2022.100431",
        "honest_note": "MSD states genetic/exploratory biomarker access requires an approved plan and agreement.",
    },
    {
        "record": "KEYNOTE-495 / KeyImPaCT",
        "trial": "KEYNOTE-495 / NCT03516981",
        "lung_scope": "yes (NSCLC)",
        "assay_or_data": "tumor RNA-derived 18-gene T-cell-inflamed GEP plus TMB",
        "processed_patient_level_matrix": "not found",
        "access": "controlled/on-request; aggregate tables open",
        "downloaded": "no",
        "tacstd2": "not in 18-gene GEP",
        "cldn4": "not in 18-gene GEP",
        "evidence": "https://doi.org/10.1038/s41591-023-02385-6",
        "honest_note": "The paper reports biomarker groups/scores, not a gene-by-sample count matrix.",
    },
    {
        "record": "KEYNOTE-782",
        "trial": "KEYNOTE-782 / NCT03664024",
        "lung_scope": "yes (nonsquamous NSCLC)",
        "assay_or_data": "bulk tumor RNA-seq; 11 reported expression signatures; n=69 evaluable",
        "processed_patient_level_matrix": "not found",
        "access": "controlled/on-request; abstract open",
        "downloaded": "no",
        "tacstd2": "potentially measured by whole-transcriptome RNA-seq, values unavailable",
        "cldn4": "potentially measured by whole-transcriptome RNA-seq, values unavailable",
        "evidence": "https://doi.org/10.1200/JCO.2024.42.16_suppl.8578",
        "honest_note": "No repository accession or patient-level processed matrix was reported.",
    },
    {
        "record": "MSD/Vivli data-sharing route",
        "trial": "multiple sponsor KEYNOTE trials",
        "lung_scope": "yes",
        "assay_or_data": "anonymized clinical data; genetic/exploratory biomarkers by special review",
        "processed_patient_level_matrix": "not anonymously downloadable",
        "access": "application, scientific review, data-sharing agreement, secure portal/analysis",
        "downloaded": "no",
        "tacstd2": "request-dependent",
        "cldn4": "request-dependent",
        "evidence": "https://engagezone.msd.com/ds_documentation.php",
        "honest_note": "This is controlled access, not open data.",
    },
    {
        "record": "SU2C-MARK NSCLC cohort",
        "trial": "not a KEYNOTE trial",
        "lung_scope": "yes (NSCLC)",
        "assay_or_data": "processed bulk RNA-seq TPM; 152 QC-passing samples",
        "processed_patient_level_matrix": "yes, openly available",
        "access": "open processed data; raw data controlled",
        "downloaded": "no",
        "tacstd2": "expected from whole-transcriptome assay; matrix not inspected",
        "cldn4": "expected from whole-transcriptome assay; matrix not inspected",
        "evidence": "https://doi.org/10.5281/zenodo.11179623",
        "honest_note": "Useful near-match, but excluded because it is not KEYNOTE; archive is 1.18 GB.",
    },
]


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "gpt-keynote-catalog/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def refresh_geo() -> None:
    search = fetch_json(GEO_SEARCH)
    ids = search["esearchresult"]["idlist"]
    summaries = fetch_json(GEO_SUMMARY.format(ids=",".join(ids))) if ids else {"result": {}}
    payload = {
        "query_url": GEO_SEARCH,
        "query": "KEYNOTE[All Fields]",
        "search": search,
        "summaries": summaries,
        "interpretation": (
            "Literal GEO KEYNOTE hits include GSE216069/GSE216055, whose KEYNOTE-001 "
            "patient is melanoma, not lung. Other literal hits are unrelated uses of "
            "the ordinary word 'keynote'."
        ),
    }
    (RESULTS / "geo_literal_keynote_query.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def validate_downloads() -> list[dict]:
    # No current catalog row meets the complete eligibility intersection.
    eligible = [
        row
        for row in CATALOG
        if row["lung_scope"].startswith("yes")
        and row["processed_patient_level_matrix"].startswith("yes")
        and row["access"].startswith("open")
        and not row["trial"].startswith("not a KEYNOTE")
    ]
    if eligible:
        raise RuntimeError("A newly eligible row requires a curated URL and size check before download.")
    return []


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offline", action="store_true", help="skip live GEO API refresh")
    args = parser.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    write_tsv(RESULTS / "catalog.tsv", CATALOG, FIELDS)

    gene_rows = [
        {
            "record": row["record"],
            "trial": row["trial"],
            "TACSTD2": row["tacstd2"],
            "CLDN4": row["cldn4"],
            "basis": row["honest_note"],
        }
        for row in CATALOG
    ]
    write_tsv(
        RESULTS / "tacstd2_cldn4.tsv",
        gene_rows,
        ["record", "trial", "TACSTD2", "CLDN4", "basis"],
    )

    downloads = validate_downloads()
    write_tsv(
        RESULTS / "download_manifest.tsv",
        downloads,
        ["url", "filename", "bytes", "sha256", "license", "eligibility"],
    )
    if not args.offline:
        refresh_geo()

    checksums = []
    for path in sorted(RESULTS.iterdir()):
        if path.name == "checksums.sha256" or not path.is_file():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        checksums.append(f"{digest}  {path.name}")
    (RESULTS / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")

    print(f"Wrote {len(CATALOG)} catalog rows; {len(downloads)} eligible downloads.")


if __name__ == "__main__":
    main()
