#!/usr/bin/env python3
"""Honest inventory of must-try + 2024–2026 CD45+ / immune-only lung IO leftovers."""

from __future__ import annotations

import json

import pandas as pd

from config import RESULTS


ROWS = [
    {
        "accession": "GSE243013",
        "year": 2025,
        "library": "CD45+ immune atlas (computationally sorted; deposited immune MTX)",
        "n_patients_public": 243,
        "mpr_public": "yes (GEO metadata pathological_response: pCR/MPR/non-MPR)",
        "recist_public": "yes (GEO metadata radiological_response)",
        "tacstd2_cldn4_detectable": "yes (sparse; ~1.2% / 0.63% cells)",
        "include_in_meta": "yes",
        "reason": "Must-try. Only large public CD45+ lung neoadjuvant IO matrix with machine-readable MPR.",
    },
    {
        "accession": "GSE229353",
        "year": 2023,
        "library": "CD45+ bead-sorted immune 10x (tumor)",
        "n_patients_public": 7,
        "mpr_public": "not on GEO; yes in public paper Table S1 (3 MPR / 4 non-MPR)",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "tested in this extra (immune library leak)",
        "include_in_meta": "yes",
        "reason": "Must-try. Immune-only neoadjuvant libraries. MPR taken from public supplement, not invented.",
    },
    {
        "accession": "GSE154826",
        "year": 2020,
        "library": "CD45+ / CITE-seq early-stage NSCLC (Leader et al.)",
        "n_patients_public": 35,
        "mpr_public": "no",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "not tested for ICI (no ICI labels on scRNA patients)",
        "include_in_meta": "no",
        "reason": "Must-try. Treatment-naive resection atlas. ICI association in the paper is a separate bulk trial, not these 35 lesions. GEO family SOFT has no MPR/RECIST/ICI labels.",
    },
    {
        "accession": "Hui 2022 Cell Death Dis (no GEO)",
        "year": 2022,
        "library": "CD45+ bead-sorted (12 IIIA patients; 8 neoadjuvant pembro+chemo)",
        "n_patients_public": 0,
        "mpr_public": "paper describes responders vs non-responders; no public matrix",
        "recist_public": "no public matrix",
        "tacstd2_cldn4_detectable": "not public",
        "include_in_meta": "no",
        "reason": "Same Tianjin group as GSE229353. Data-availability statement is correspondence-only. No GEO accession in the article.",
    },
    {
        "accession": "GSE280232",
        "year": 2024,
        "library": "T-cell GEX+TCR (neoadjuvant ICB; KRAS/STK11), not CD45-all",
        "n_patients_public": 13,
        "mpr_public": "paper-level neoadjuvant ICB; not a CD45-all leak library",
        "recist_public": "not used",
        "tacstd2_cldn4_detectable": "not scored (T-sorted leftover, 1.2 GB RAW)",
        "include_in_meta": "no",
        "reason": "2024 leftover opened. Sorted T cells, not CD45+ immune-all. Wrong compartment for epithelial-leak fraction.",
    },
    {
        "accession": "GSE176021",
        "year": 2021,
        "library": "neoantigen-focused T cells (neoadjuvant nivo)",
        "n_patients_public": "trial subset",
        "mpr_public": "yes in Caushi/Forde papers",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "not scored (T-sorted, not CD45-all)",
        "include_in_meta": "no",
        "reason": "Immune leftover but CD8/T-focused, pre-2024, not CD45-all leak.",
    },
    {
        "accession": "GSE303680",
        "year": 2025,
        "library": "PD1-IL2v T-cell response in human lung cancer",
        "n_patients_public": 10,
        "mpr_public": "no neoadjuvant MPR design",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "not a CD45-all neoadjuvant atlas",
        "include_in_meta": "no",
        "reason": "2025 leftover opened. Cytokine-delivery T-cell study, not CD45+ MPR leak.",
    },
    {
        "accession": "GSE207422",
        "year": 2023,
        "library": "unsorted TME scRNA (epithelial + immune)",
        "n_patients_public": 15,
        "mpr_public": "yes",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "malignant RNA is User A3; taken as given",
        "include_in_meta": "no",
        "reason": "Not an immune-only library. A3 epithelial claim is not re-analyzed here.",
    },
    {
        "accession": "GSE241934 / GSE291670",
        "year": "2023/2025",
        "library": "unsorted TME scRNA with author epithelium",
        "n_patients_public": "35 / 6",
        "mpr_public": "yes",
        "recist_public": "no",
        "tacstd2_cldn4_detectable": "malignant/epithelial compartment (other extras)",
        "include_in_meta": "no",
        "reason": "Not CD45+ / immune-only libraries.",
    },
]


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(ROWS)
    df.to_csv(RESULTS / "inventory.tsv", sep="\t", index=False)
    (RESULTS / "inventory.json").write_text(json.dumps(ROWS, indent=2))
    print(df[["accession", "include_in_meta", "reason"]].to_string(index=False))


if __name__ == "__main__":
    main()
