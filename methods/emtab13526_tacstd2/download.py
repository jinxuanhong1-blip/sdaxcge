#!/usr/bin/env python3
"""Download public E-MTAB-13526 processed 10x matrices (CD235a- tumor only).

The author-annotated h5ads are 45 GB + 58 GB and are not required. This
script pulls the deposited Cell Ranger mtx files for the CD235a- (RBC-depleted,
not CD45-sorted) tumor lanes plus the shared features table and SDRF.

Public FTP (BioStudies/ArrayExpress):
  https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

FTP = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files"

# Unsorted / RBC-depleted tumor lanes (FACS = CD235a-). These are the only
# deposited tumor matrices that can contain both malignant epithelium and T/NK.
CD235A_TUMOR = [
    "P4_T2",
    "P4_T3",
    "P8_T2",
    "P15_T2",
    "P16_T2",
    "P17_T2",
    "P17_T3",
    "P18_T2",
    "P19_T2",
    "P20_T2",
    "P21_T1",
    "P21_T2",
    "P22_T1",
    "P23_T1",
    "P24_T1",
]


def curl(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return
    part = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "curl",
        "-fL",
        "--retry",
        "5",
        "--retry-delay",
        "8",
        "--retry-all-errors",
        "-o",
        str(part),
        url,
    ]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    part.rename(dest)
    print(f"OK {dest.name} ({dest.stat().st_size} bytes)", flush=True)


def write_sample_table(sdrf_path: Path, out_path: Path) -> None:
    with sdrf_path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    seen: dict[str, dict] = {}
    for row in rows:
        sn = row["Source Name"]
        if sn not in seen:
            seen[sn] = {
                "sample": sn,
                "patient": row["Characteristics[individual]"],
                "disease": row["Characteristics[disease]"],
                "facs": row["Characteristics[FACS]"],
                "sampling_site": row["Characteristics[sampling site]"],
                "tumor_grading": row["Characteristics[tumor grading]"],
                "sex": row["Characteristics[sex]"],
                "age": row["Characteristics[age]"],
                "subset": (
                    "cd235a_tumor"
                    if sn in CD235A_TUMOR
                    else "other_deposited_lane"
                ),
            }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(next(iter(seen.values())).keys())
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for key in sorted(seen):
            w.writerow(seen[key])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", type=Path, default=Path("data/emtab13526/raw"))
    args = ap.parse_args()
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)

    curl(f"{FTP}/E-MTAB-13526.sdrf.txt", dest / "E-MTAB-13526.sdrf.txt")
    curl(f"{FTP}/E-MTAB-13526.idf.txt", dest / "E-MTAB-13526.idf.txt")
    # Features are identical across lanes (Cell Ranger GRCh38 3.1.0).
    curl(f"{FTP}/P4_T2-features.tsv.gz", dest / "features.tsv.gz")
    if not (dest / "P4_T2-features.tsv.gz").exists():
        (dest / "P4_T2-features.tsv.gz").write_bytes((dest / "features.tsv.gz").read_bytes())

    for sample in CD235A_TUMOR:
        curl(f"{FTP}/{sample}-matrix.mtx.gz", dest / f"{sample}-matrix.mtx.gz")

    write_sample_table(dest / "E-MTAB-13526.sdrf.txt", dest.parent / "sample_metadata.tsv")
    manifest = {
        "accession": "E-MTAB-13526",
        "ftp": FTP,
        "subset": "CD235a- tumor Cell Ranger mtx (not the 45/58 GB annotated h5ads)",
        "n_lanes": len(CD235A_TUMOR),
        "lanes": CD235A_TUMOR,
        "skipped": [
            "10X_Lung_Tumour_Annotated_v2.h5ad (58.7 GB)",
            "10X_Lung_Healthy_Background_Annotated_v2.h5ad (45.5 GB)",
            "CD45+/MDSC-sorted tumor and all background/donor matrices",
        ],
    }
    (dest.parent / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
