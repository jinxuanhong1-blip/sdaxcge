#!/usr/bin/env python3
"""Stream TACSTD2/CLDN4/TJ-core genes + full library size from GSE131907 UMI matrix."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

IDENTITY = [
    "TACSTD2",
    "CLDN4",
    "ELF3",
    "EPCAM",
    "WT1",
    "CALB2",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "NKG7",
    "MS4A1",
    "CD68",
]


def load_tj(path: Path) -> list[str]:
    genes = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        genes.append(line.split()[0])
    return genes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/gse131907/paper_funnel_mpe"))
    ap.add_argument(
        "--tj-list",
        type=Path,
        default=Path("methods/paper_funnel_mpe_tacstd2_tj/data/tj_core_genes.txt"),
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tj = load_tj(args.tj_list)
    genes = list(dict.fromkeys(IDENTITY + tj))
    wanted = set(genes)
    matrix = args.datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    found: dict[str, np.ndarray] = {}
    print(f"streaming {len(genes)} genes + library size from {matrix}", flush=True)
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[0] != "Index":
            raise SystemExit(f"unexpected header start: {header[0]!r}")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        libsize = np.zeros(n_cells, dtype=np.float64)
        n_rows = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n_cells}")
            libsize += arr.astype(np.float64)
            n_rows += 1
            if gene in wanted and gene not in found:
                found[gene] = arr
            if n_rows % 5000 == 0:
                print(f"  rows {n_rows} found {len(found)}/{len(genes)}", flush=True)

    missing = [g for g in genes if g not in found]
    # PATJ is often annotated as INADL; keep INADL if PATJ absent.
    if "PATJ" in missing and "INADL" in found:
        missing = [g for g in missing if g != "PATJ"]
        genes = [g for g in genes if g != "PATJ"]
    if missing:
        raise SystemExit(f"missing genes: {missing}")

    counts = np.vstack([found[g] for g in genes]).astype(np.float32)
    out = args.outdir / "gene_counts.npz"
    np.savez_compressed(
        out,
        genes=np.array(genes),
        barcodes=np.array(cell_ids),
        counts=counts,
        libsize=libsize,
        tj_genes=np.array([g for g in tj if g in genes]),
    )
    meta = {
        "n_cells": n_cells,
        "n_gene_rows_scanned": n_rows,
        "genes": genes,
        "tj_genes": [g for g in tj if g in genes],
        "median_nCount_full": float(np.median(libsize)),
        "source": matrix.name,
        "normalization_downstream": "log1p(UMI / nCount_full * 10000)",
    }
    (args.outdir / "extract_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2), flush=True)
    print(f"wrote {out} ({out.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
