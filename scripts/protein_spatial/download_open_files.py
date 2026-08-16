#!/usr/bin/env python3
"""Download open processed files only (skip raw MS/FASTQ and files >2GB)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from urllib.request import Request, urlopen

from config import CPTAC_BASE, CPTAC_COHORTS, CPTAC_PROTEIN_SUFFIX, DATA, MAX_BYTES

FILES = [
    # CPTAC freeze v1.2 (treatment-naive LUAD/LSCC)
    *[
        (
            f"cptac/{c}{CPTAC_PROTEIN_SUFFIX}",
            f"{CPTAC_BASE}/{c}/{c}{CPTAC_PROTEIN_SUFFIX}",
        )
        for c in CPTAC_COHORTS
    ],
    *[
        (f"cptac/{c}_{kind}.txt", f"{CPTAC_BASE}/{c}/{c}_{kind}.txt")
        for c in CPTAC_COHORTS
        for kind in ("phenotype", "survival", "meta")
    ],
    # PXD042091 processed SEARCH tables + spectral library (skip .wiff raw)
    (
        "pxd042091/2_dat_cs_all.txt",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/2_dat_cs_all.txt",
    ),
    (
        "pxd042091/6_data_tf_response.txt",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/6_data_tf_response.txt",
    ),
    (
        "pxd042091/NSClibrary.txt",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/NSClibrary.txt",
    ),
    # PXD059688 processed mzTab only (mgf 8.5GB and msf 37GB skipped)
    (
        "pxd059688/20230904_Tumor_AA.mzTab.gz",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.mzTab.gz",
    ),
    # GSE271689 GeoMx DCC + metadata (skip SRA FASTQ)
    (
        "gse271689/GSE271689_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/matrix/GSE271689_series_matrix.txt.gz",
    ),
    (
        "gse271689/GSE271689_family.soft.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/soft/GSE271689_family.soft.gz",
    ),
    (
        "gse271689/filelist.txt",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/filelist.txt",
    ),
    (
        "gse271689/GSE271689_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/GSE271689_RAW.tar",
    ),
    # Official NanoString Hs WTA v1.0 PKC (probe→gene map for DCC RTS IDs)
    (
        "geomx_pkc/Hs_R_NGS_WTA_v1.0.pkc",
        "https://zenodo.org/records/12752405/files/Hs_R_NGS_WTA_v1.0.pkc?download=1",
    ),
]

EMTAB_FTP = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files"
EMTAB_SECTIONS = [
    "D1_1",
    "D1_2",
    "D2_1",
    "D2_2",
    "P10_B1",
    "P10_B2",
    "P10_T1",
    "P10_T2",
    "P10_T3",
    "P10_T4",
    "P11_B1",
    "P11_B2",
    "P11_T1",
    "P11_T2",
    "P11_T3",
    "P11_T4",
    "P15_B1",
    "P15_B2",
    "P15_T1",
    "P15_T2",
    "P16_B1",
    "P16_B2",
    "P16_T1",
    "P16_T2",
    "P17_B1",
    "P17_B2",
    "P17_T1",
    "P17_T2",
    "P19_B1",
    "P19_B2",
    "P19_T1",
    "P19_T2",
    "P24_B1",
    "P24_B2",
    "P24_T1",
    "P24_T2",
    "P25_B1",
    "P25_B2",
    "P25_T1",
    "P25_T2",
]


def head_length(url: str) -> int | None:
    req = Request(url, method="HEAD")
    try:
        with urlopen(req, timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def curl(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest}")
        return
    size = head_length(url)
    if size is not None and size > MAX_BYTES:
        print(f"SKIP >2GB {size} {url}")
        return
    cmd = [
        "curl",
        "-fL",
        "--retry",
        "4",
        "--retry-delay",
        "4",
        "-o",
        str(dest),
        url,
    ]
    print("GET", url)
    subprocess.check_call(cmd)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    for rel, url in FILES:
        curl(url, DATA / rel)
    curl(f"{EMTAB_FTP}/E-MTAB-13530.sdrf.txt", DATA / "emtab13530/E-MTAB-13530.sdrf.txt")
    curl(f"{EMTAB_FTP}/E-MTAB-13530.idf.txt", DATA / "emtab13530/E-MTAB-13530.idf.txt")
    for sec in EMTAB_SECTIONS:
        name = f"{sec}-filtered_feature_bc_matrix.h5"
        curl(f"{EMTAB_FTP}/{name}", DATA / "emtab13530" / name)
    manifest = {
        "note": "Local cache only; not committed. Raw MS/FASTQ and files >2GB were skipped.",
        "skipped": [
            "PXD042091 *.wiff / *.wiff.scan (raw SWATH)",
            "PXD059688 20230904_Tumor_AA.mgf (8.49 GB)",
            "PXD059688 20230904_Tumor_AA.msf (37.9 GB)",
            "GSE271689 SRA FASTQ",
            "E-MTAB-13530 ENA FASTQ / spatial.tar / web_summary.html",
        ],
    }
    (DATA / "DOWNLOAD_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
