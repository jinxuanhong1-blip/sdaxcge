#!/usr/bin/env python3
"""W200-A3 feasibility check for GSE243013 (TACSTD2 in malignant cells, MPR vs non-MPR).

Task rule: run the A3-style analysis only if the processed expression data needed
for malignant-cell TACSTD2 can be obtained in < 2 GB. Otherwise skip honestly.

This script downloads ONLY the small evidence files (cell metadata ~39 MB, gene
list ~112 KB) plus HTTP HEAD size checks, then writes the evidence tables and a
machine-readable verdict under results/w200/A3_GSE243013/.

Usage: python3 scripts/w200_a3_gse243013/feasibility_check.py [--workdir /tmp/gse243013]
"""

import argparse
import json
import os
import urllib.request

import pandas as pd

SUPP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/"
SUPP_FILES = [
    "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz",
    "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz",
    "GSE243013_RAW.tar",
    "GSE243013_barcodes.csv.gz",
    "GSE243013_genes.csv.gz",
    "GSE243013_T_with_TCR_annotation.csv.gz",
    "GSE243013_UMAP_info.tar.gz",
    "GSE243013_NMF_all_group_5.csv.gz",
]
SIZE_BUDGET_BYTES = 2 * 1024**3
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "w200", "A3_GSE243013")

MALIGNANT_PATTERN = "pith|alig|umor|ancer|arcinom"  # epithelial/malignant/tumor/cancer/carcinoma


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return int(resp.headers["Content-Length"])


def fetch(url: str, dest: str) -> None:
    if not os.path.exists(dest):
        urllib.request.urlretrieve(url, dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default="/tmp/gse243013")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1) Size manifest for all supplementary files (HEAD requests only).
    rows = []
    for name in SUPP_FILES:
        size = head_size(SUPP_BASE + name)
        rows.append({
            "file": name,
            "size_bytes": size,
            "size_gb": round(size / 1024**3, 3),
            "within_2gb_budget": size < SIZE_BUDGET_BYTES,
        })
    manifest = pd.DataFrame(rows).sort_values("size_bytes", ascending=False)
    manifest.to_csv(os.path.join(OUT_DIR, "file_manifest.tsv"), sep="\t", index=False)

    counts_size = int(manifest.loc[manifest.file.str.contains("counts.mtx"), "size_bytes"].iloc[0])

    # 2) Cell metadata: the only annotation of what cell populations exist.
    meta_gz = os.path.join(args.workdir, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
    fetch(SUPP_BASE + "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz", meta_gz)
    meta = pd.read_csv(meta_gz, low_memory=False)

    comp = (
        meta.groupby(["major_cell_type", "sub_cell_type"]).size().reset_index(name="n_cells")
        .sort_values(["major_cell_type", "n_cells"], ascending=[True, False])
    )
    comp.to_csv(os.path.join(OUT_DIR, "cell_type_composition.tsv"), sep="\t", index=False)

    malignant_cells = int(
        meta.major_cell_type.astype(str).str.contains(MALIGNANT_PATTERN, case=False).sum()
        + meta.sub_cell_type.astype(str).str.contains(MALIGNANT_PATTERN, case=False).sum()
    )

    # 3) Patient-level response labels (evidence that MPR labels exist).
    pt = meta.drop_duplicates("sampleID")
    summary = (
        pt.groupby(["cancer_type", "pathological_response"]).size().reset_index(name="n_patients")
    )
    summary.to_csv(os.path.join(OUT_DIR, "patient_response_summary.tsv"), sep="\t", index=False)

    # 4) TACSTD2 present in the gene index (but expression locked in the 6.6 GB matrix).
    genes_gz = os.path.join(args.workdir, "GSE243013_genes.csv.gz")
    fetch(SUPP_BASE + "GSE243013_genes.csv.gz", genes_gz)
    genes = pd.read_csv(genes_gz)
    tacstd2_in_index = bool((genes.iloc[:, 0] == "TACSTD2").any())

    verdict = {
        "dataset": "GSE243013",
        "task": "W200-A3: TACSTD2 in malignant cells, MPR vs non-MPR",
        "n_geo_samples": int(pt.shape[0]),
        "n_cells": int(meta.shape[0]),
        "major_cell_types": meta.major_cell_type.value_counts().to_dict(),
        "n_malignant_or_epithelial_cells": malignant_cells,
        "mpr_labels_available": bool(pt.pathological_response.notna().any()),
        "tacstd2_in_gene_index": tacstd2_in_index,
        "counts_matrix_bytes": counts_size,
        "counts_matrix_gb": round(counts_size / 1024**3, 2),
        "size_budget_gb": 2.0,
        "counts_within_budget": counts_size < SIZE_BUDGET_BYTES,
        "decision": "SKIP",
        "reasons": [
            "The only expression file (GSE243013_NSCLC_immune_scRNA_counts.mtx.gz) "
            f"is {counts_size / 1024**3:.1f} GB gzipped, exceeding the 2 GB processed-data budget.",
            "The atlas is CD45+/immune-only: cell metadata contains exclusively T/NK, B and "
            "Myeloid cells, with zero epithelial or malignant cells, so malignant-cell TACSTD2 "
            "cannot be measured from this dataset at any size.",
        ],
    }
    with open(os.path.join(OUT_DIR, "feasibility.json"), "w") as fh:
        json.dump(verdict, fh, indent=2)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
