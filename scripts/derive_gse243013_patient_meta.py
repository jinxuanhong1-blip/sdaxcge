#!/usr/bin/env python3
"""Collapse cell-level GEO metadata to one row per sampleID (no expression)."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path


FIELDS = [
    "sampleID",
    "n_cells",
    "pathological_response",
    "pathological_response_rate",
    "radiological_response",
    "cancer_type",
    "gender",
    "age",
    "smoking_history",
    "pre_treatment_staging",
    "anti-PD1_therapy",
    "chemotherapy",
    "targeted_therapy",
    "cycles",
]


def main() -> None:
    src = Path("results/noskip/GSE243013/source/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
    dst = Path("results/noskip/GSE243013/source/patient_metadata_from_geo.tsv")
    patients = {}
    with gzip.open(src, "rt", newline="") as f:
        for row in csv.DictReader(f):
            sid = row["sampleID"]
            if sid not in patients:
                rec = {k: row.get(k, "") for k in FIELDS if k != "n_cells"}
                rec["n_cells"] = 0
                patients[sid] = rec
            patients[sid]["n_cells"] += 1
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for sid in sorted(patients):
            w.writerow(patients[sid])
    print(f"wrote {dst} n={len(patients)}")


if __name__ == "__main__":
    main()
