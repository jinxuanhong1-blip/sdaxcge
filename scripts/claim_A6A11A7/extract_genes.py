#!/usr/bin/env python3
"""Stream gene x cell UMI matrices and keep a small gene panel + library sizes."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

GENES = [
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "CLDN4",
    "CLDN7",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "NCAM1",
    "GZMB",
    "GZMK",
    "PRF1",
    "PTPRC",
    "MS4A1",
    "CD79A",
    "LYZ",
    "COL1A1",
    "PECAM1",
    "LGALS1",
    "LGALS3",
    "LGALS9",
    "NECTIN1",
    "NECTIN2",
    "NECTIN3",
    "NECTIN4",
    "PVR",
    "PVRL2",
    "TGFB1",
    "TGFB2",
    "TGFB3",
    "TGFBR1",
    "TGFBR2",
    "CD47",
    "SIRPA",
    "CD274",
    "PDCD1LG2",
    "PDCD1",
    "HAVCR2",
    "TIGIT",
]


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def extract(matrix_path: Path, out_prefix: Path) -> dict:
    wanted = set(GENES)
    found: dict[str, np.ndarray] = {}
    with open_text(matrix_path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        totals = np.zeros(n_cells, dtype=np.float64)
        n_genes = 0
        for line in handle:
            gene, sep, values = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(values, sep="\t", dtype=np.float64)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} != {n_cells}")
            totals += arr
            n_genes += 1
            if gene in wanted:
                found[gene] = arr
    missing = sorted(wanted - set(found))
    aliases = {"PVRL2": "NECTIN2"}
    for alias, canon in aliases.items():
        if alias in found and canon not in found:
            found[canon] = found[alias]
    expr = pd.DataFrame({g: found[g] for g in GENES if g in found}, index=cell_ids)
    expr.index.name = "cell_id"
    lib = pd.DataFrame({"n_umi": totals}, index=cell_ids)
    lib.index.name = "cell_id"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    expr.to_csv(out_prefix.with_name(out_prefix.name + ".genes.csv.gz"))
    lib.to_csv(out_prefix.with_name(out_prefix.name + ".lib.csv.gz"))
    meta = {
        "matrix": str(matrix_path),
        "n_cells": n_cells,
        "n_genes_scanned": n_genes,
        "genes_found": sorted(found),
        "genes_missing": missing,
        "median_umi": float(np.median(totals)),
    }
    out_prefix.with_suffix(".extract.json").write_text(json.dumps(meta, indent=2))
    return meta


if __name__ == "__main__":
    root = Path("/tmp/data")
    out = Path("/tmp/data/extracted")
    jobs = [
        (root / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", out / "GSE131907"),
        (root / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", out / "GSE207422"),
    ]
    for matrix, prefix in jobs:
        print("extracting", matrix)
        print(json.dumps(extract(matrix, prefix), indent=2))
