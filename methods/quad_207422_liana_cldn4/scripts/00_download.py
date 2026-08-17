#!/usr/bin/env python3
"""Download public files for the QUAD CLDN4 LIANA merge.

GSE207422 + GSE131907 + GSE148071 + GSE205335 + CellPhoneDB v5.
No EGA / GSA raw FASTQ.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
DEFAULT_DATA = Path("/tmp/quad_207422_liana")

FILES = [
    {
        "cohort": "GSE207422",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "dest": "gse207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "role": "processed_umi",
    },
    {
        "cohort": "GSE207422",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "dest": "gse207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "role": "sample_metadata",
    },
    {
        "cohort": "GSE131907",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "dest": "gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "role": "author_annotation",
    },
    {
        "cohort": "GSE131907",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
        "dest": "gse131907/GSE131907_series_matrix.txt.gz",
        "role": "series_matrix",
    },
    {
        "cohort": "GSE131907",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "dest": "gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "role": "processed_umi",
    },
    {
        "cohort": "GSE148071",
        "url": "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/NSCLC_GSE148071_CellMetainfo_table.tsv",
        "dest": "gse148071/NSCLC_GSE148071_CellMetainfo_table.tsv",
        "role": "tisch2_metainfo",
    },
    {
        "cohort": "GSE148071",
        "url": "https://tisch.compbio.cn/static/data/NSCLC_GSE148071/NSCLC_GSE148071_expression.h5",
        "dest": "gse148071/NSCLC_GSE148071_expression.h5",
        "role": "tisch2_expression_h5",
    },
    {
        "cohort": "GSE205335",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "dest": "gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "role": "processed_umi",
    },
    {
        "cohort": "GSE205335",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "dest": "gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "role": "author_identity",
    },
    {
        "cohort": "GSE205335",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
        "dest": "gse205335/GSE205335_family.soft.gz",
        "role": "geo_soft",
    },
    {
        "cohort": "CellPhoneDB",
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv",
        "dest": "cellphonedb/interaction_input.csv",
        "role": "cpdb_interactions",
    },
    {
        "cohort": "CellPhoneDB",
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv",
        "dest": "cellphonedb/gene_input.csv",
        "role": "cpdb_genes",
    },
    {
        "cohort": "CellPhoneDB",
        "url": "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv",
        "dest": "cellphonedb/complex_input.csv",
        "role": "cpdb_complexes",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET  {url}", flush=True)
    subprocess.run(
        ["curl", "-fL", "--retry", "5", "--retry-delay", "8", "--retry-all-errors", "-o", str(tmp), url],
        check=True,
    )
    tmp.replace(dest)
    print(f"OK   {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in FILES:
        dest = args.datadir / spec["dest"]
        fetch(spec["url"], dest)
        rows.append(
            {
                "cohort": spec["cohort"],
                "file": spec["dest"],
                "role": spec["role"],
                "bytes": dest.stat().st_size,
                "sha256": sha256(dest),
                "url": spec["url"],
            }
        )
    manifest = args.datadir / "file_manifest.json"
    manifest.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"wrote {manifest}", flush=True)
    print("skipped_raw=HRA001033,EGAD00001005054,EGAD00001008703", flush=True)


if __name__ == "__main__":
    main()
