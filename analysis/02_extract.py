#!/usr/bin/env python3
"""Stream selected genes from both processed GSE131907 matrices.

Neither matrix is skipped for size. Genes are streamed; the full gene x cell
tables are not loaded into RAM. Output is a compact per-cell table.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SELECTED = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "NKX2-1",
    "NAPSA",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "NKG7",
    "CD79A",
    "MS4A1",
    "CD68",
    "LYZ",
    "FCGR3A",
    "KIT",
    "PECAM1",
    "COL1A1",
    "CD274",
    "PDCD1",
    "HAVCR2",
    "TIGIT",
    "LAG3",
    "CXCL9",
    "CXCL10",
    "STAT1",
    "B2M",
    "HLA-A",
    "HLA-B",
    "HLA-DRA",
]


def stream_selected_genes(matrix_path: Path, wanted: set[str]) -> tuple[list[str], dict[str, np.ndarray], int]:
    found: dict[str, np.ndarray] = {}
    n_genes = 0
    t0 = time.time()
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[0] != "Index":
            raise SystemExit(f"{matrix_path}: unexpected header start {header[0]!r}")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        print(f"{matrix_path.name}: {n_cells} cells, seeking {len(wanted)} genes", flush=True)
        for line in handle:
            n_genes += 1
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene in wanted:
                arr = np.fromstring(rest, sep="\t", dtype=np.float32)
                if arr.size != n_cells:
                    raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
                found[gene] = arr
                print(
                    f"  found {gene} ({len(found)}/{len(wanted)}) after {n_genes} rows "
                    f"[{time.time() - t0:.1f}s]",
                    flush=True,
                )
                if len(found) == len(wanted):
                    # Still count remaining genes so the inventory is honest.
                    for extra in handle:
                        if extra.strip():
                            n_genes += 1
                    break
            if n_genes % 2500 == 0:
                print(f"  scanned {n_genes} gene rows, found {len(found)} [{time.time() - t0:.1f}s]", flush=True)
    return cell_ids, found, n_genes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "extracted")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    wanted = set(SELECTED)
    matrices = {
        "umi": args.datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "log2tpm": args.datadir / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
    }
    for key, path in matrices.items():
        if not path.exists():
            raise SystemExit(f"missing {path}")

    umi_ids, umi, n_umi_genes = stream_selected_genes(matrices["umi"], wanted)
    tpm_ids, tpm, n_tpm_genes = stream_selected_genes(matrices["log2tpm"], wanted)
    if umi_ids != tpm_ids:
        raise SystemExit("UMI and log2TPM cell-ID order differ")

    missing_umi = sorted(wanted - set(umi))
    missing_tpm = sorted(wanted - set(tpm))
    print("UMI missing", missing_umi)
    print("log2TPM missing", missing_tpm)

    ann = pd.read_csv(
        args.datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
        dtype=str,
    )
    per_cell = ann.set_index("Index").reindex(umi_ids).reset_index()
    if per_cell["Sample"].isna().any():
        n_miss = int(per_cell["Sample"].isna().sum())
        raise SystemExit(f"{n_miss} matrix cell IDs missing from annotation")

    for gene in SELECTED:
        if gene in umi:
            per_cell[f"{gene}_umi"] = umi[gene]
        if gene in tpm:
            per_cell[f"{gene}_log2tpm"] = tpm[gene]

    out_table = args.outdir / "per_cell_selected_genes.csv.gz"
    per_cell.to_csv(out_table, index=False)
    inventory = {
        "n_cells": len(umi_ids),
        "n_genes_umi_matrix": n_umi_genes,
        "n_genes_log2tpm_matrix": n_tpm_genes,
        "genes_requested": SELECTED,
        "genes_found_umi": sorted(umi),
        "genes_found_log2tpm": sorted(tpm),
        "genes_missing_umi": missing_umi,
        "genes_missing_log2tpm": missing_tpm,
        "umi_bytes": matrices["umi"].stat().st_size,
        "log2tpm_bytes": matrices["log2tpm"].stat().st_size,
        "matrices_used": [p.name for p in matrices.values()],
        "skipped_for_size": [],
    }
    (args.outdir / "extract_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    print(json.dumps(inventory, indent=2))
    print(f"wrote {out_table} rows={len(per_cell)}", file=sys.stderr)


if __name__ == "__main__":
    main()
