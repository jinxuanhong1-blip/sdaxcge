#!/usr/bin/env python3
"""Stream TACSTD2, CLDN4, TJ genes, and full-library sizes from GSE131907 UMI text.

Writes genes x cells float32 counts plus barcodes and libsize. Does not write
cell-level p-values.
"""
from __future__ import annotations

import json
from pathlib import Path

import gzip
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EXTRA = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "MUC1",
    "CDH1",
    "ELF3",
    "CLDN1",
    "CLDN3",
    "CLDN7",
    "CLDN8",
    "CLDN18",
    "OCLN",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "CGN",
    "CGNL1",
    "MARVELD2",
    "CRB3",
    "PARD3",
    "PARD6A",
    "PARD6B",
    "AMOT",
    "AMOTL1",
    "MAGI1",
    "MAGI3",
    "MPDZ",
    "JAM2",
    "JAM3",
    "INADL",
    "PATJ",
    "CD3D",
    "CD3E",
    "CD8A",
    "NKG7",
]


def wanted() -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    kegg = (DATA / "kegg_tj_no_cldn4.txt").read_text().splitlines()
    for g in kegg + EXTRA:
        g = g.strip()
        if not g or g.startswith("#") or g in seen:
            continue
        seen.add(g)
        genes.append(g)
    return genes


def main() -> None:
    datadir = Path("/tmp/gse131907")
    outdir = datadir / "maxeffect_subset"
    outdir.mkdir(parents=True, exist_ok=True)
    genes = wanted()
    need = set(genes)
    matrix = datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    size = matrix.stat().st_size
    if size != 408736818:
        raise SystemExit(f"UMI matrix size {size} != 408736818")
    print(f"streaming {len(need)} genes + libsize from {matrix}", flush=True)

    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[0] != "Index":
            raise SystemExit(f"unexpected header field {header[0]}")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        if n_cells != 208506:
            raise SystemExit(f"expected 208506 cells, got {n_cells}")
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
            if gene in need and gene not in found:
                found[gene] = arr
            if n_rows % 4000 == 0:
                print(f"  rows {n_rows} found {len(found)}", flush=True)
        print(f"scanned {n_rows} gene rows", flush=True)

    present = [g for g in genes if g in found]
    missing = [g for g in genes if g not in found]
    print(f"found {len(present)} missing {len(missing)}: {missing}", flush=True)
    mat = np.vstack([found[g] for g in present]).astype(np.float32)
    np.save(outdir / "counts.npy", mat)
    np.save(outdir / "libsize.npy", libsize.astype(np.float64))
    (outdir / "genes.txt").write_text("\n".join(present) + "\n")
    (outdir / "barcodes.txt").write_text("\n".join(cell_ids) + "\n")
    (outdir / "missing_genes.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    meta = {
        "n_cells": n_cells,
        "n_gene_rows_scanned": n_rows,
        "n_genes_found": len(present),
        "missing": missing,
        "median_nCount_full": float(np.median(libsize)),
        "source": "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "bytes": size,
    }
    (outdir / "extract_meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
