#!/usr/bin/env python3
"""Stream locked family genes + full-library sizes from GSE131907 UMI text.

Helper only. Primary analysis is R + Seurat CreateSeuratObject.
Does not require scipy: writes Matrix Market by hand.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gzip
import numpy as np


EXTRA = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "NKG7",
    "GNLY",
    "STAT1",
    "IRF1",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "OCLN",
    "TJP1",
]


def wanted_genes(families_path: Path) -> list[str]:
    fam = json.loads(families_path.read_text())
    genes: list[str] = []
    seen: set[str] = set()
    for key in ("IFN", "MHC_I_APM", "TJ_no_CLDN4", "chemokine"):
        for g in fam[key]:
            if g not in seen:
                seen.add(g)
                genes.append(g)
    for g in EXTRA:
        if g not in seen:
            seen.add(g)
            genes.append(g)
    return genes


def write_mtx(path: Path, rows: np.ndarray, cols: np.ndarray, data: np.ndarray, shape: tuple[int, int]) -> None:
    """1-based coordinate Matrix Market, genes x cells."""
    nnz = int(data.size)
    with path.open("w") as out:
        out.write("%%MatrixMarket matrix coordinate real general\n")
        out.write(f"{shape[0]} {shape[1]} {nnz}\n")
        if nnz:
            # MTX is 1-indexed
            np.savetxt(
                out,
                np.column_stack((rows + 1, cols + 1, data)),
                fmt="%d %d %.6g",
            )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/gse131907/subset"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    families = Path(__file__).resolve().parents[1] / "data" / "families.json"
    if not families.exists():
        families = args.datadir / "families.json"
    genes = wanted_genes(families)
    wanted = set(genes)
    matrix = args.datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    print(f"streaming {len(wanted)} genes + full libsize from {matrix}", flush=True)

    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
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
            libsize += arr
            n_rows += 1
            if gene in wanted:
                found[gene] = arr
        print(f"scanned {n_rows} gene rows", flush=True)

    missing = sorted(wanted - set(found))
    present = [g for g in genes if g in found]
    print(f"found {len(found)} missing {len(missing)}", flush=True)

    row_idx = []
    col_idx = []
    vals = []
    for i, gene in enumerate(present):
        arr = found[gene]
        nz = np.flatnonzero(arr)
        if nz.size == 0:
            continue
        row_idx.append(np.full(nz.size, i, dtype=np.int32))
        col_idx.append(nz.astype(np.int32))
        vals.append(arr[nz].astype(np.float32))
    if row_idx:
        rows = np.concatenate(row_idx)
        cols = np.concatenate(col_idx)
        data = np.concatenate(vals)
    else:
        rows = np.zeros(0, dtype=np.int32)
        cols = np.zeros(0, dtype=np.int32)
        data = np.zeros(0, dtype=np.float32)

    write_mtx(args.outdir / "matrix.mtx", rows, cols, data, (len(present), len(cell_ids)))
    (args.outdir / "features.tsv").write_text("\n".join(present) + "\n")
    (args.outdir / "barcodes.tsv").write_text("\n".join(cell_ids) + "\n")
    np.savetxt(
        args.outdir / "libsize.tsv",
        np.column_stack((np.array(cell_ids), libsize.astype(int))),
        fmt="%s\t%s",
        header="cell\tnCount_full",
        comments="",
    )
    (args.outdir / "missing_genes.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    meta = {
        "n_cells": n_cells,
        "n_gene_rows_scanned": n_rows,
        "n_genes_wanted": len(wanted),
        "n_genes_found": len(found),
        "missing": missing,
        "median_nCount_full": float(np.median(libsize)),
        "source": "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    }
    (args.outdir / "extract_meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
