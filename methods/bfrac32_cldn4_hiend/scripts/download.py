#!/usr/bin/env python3
"""Download the two public scRNA matrices named in the given B-frac combo.

Cohorts (from methods/scrna_cldn4_combo and methods/scrna_tls_meta; not re-audited):
  GSE131907  Kim et al. Nat Commun 2020  — author-malignant n=21
  GSE241934 IIT  NEOTIDE neoadjuvant IO — author-epi n=11

GSE207422 is not downloaded (PR #400 B-frac is NS; this slice is the n=32 combo).
GSE241934 Real/RWC is not downloaded (not in the given combo).
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

FILES = {
    "GSE131907": [
        (
            "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        ),
        (
            "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        ),
        (
            "GSE131907_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
        ),
    ],
    "GSE241934_IIT": [
        (
            "GSE241934_IIT_Matrix.mtx.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        ),
        (
            "GSE241934_IIT_Meta.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        ),
        (
            "GSE241934_IIT_barcodes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        ),
        (
            "GSE241934_IIT_features.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        ),
    ],
}


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = [
        "wget",
        "-c",
        "--tries=8",
        "--waitretry=8",
        "--timeout=60",
        "--progress=dot:giga",
        "-O",
        str(tmp),
        url,
    ]
    subprocess.check_call(cmd)
    tmp.replace(dest)
    print(f"OK {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/bfrac32_cldn4_hiend"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    prov = []
    for cohort, items in FILES.items():
        d = args.outdir / cohort
        d.mkdir(parents=True, exist_ok=True)
        for name, url in items:
            dest = d / name
            print(f"GET {cohort} {name}", flush=True)
            wget(url, dest)
            prov.append(
                {
                    "cohort": cohort,
                    "filename": name,
                    "url": url,
                    "bytes": dest.stat().st_size,
                }
            )
    (args.outdir / "provenance_downloads.json").write_text(json.dumps(prov, indent=2) + "\n")
    print("wrote", args.outdir / "provenance_downloads.json", flush=True)


if __name__ == "__main__":
    main()
