#!/usr/bin/env python3
"""Stream GSE241934 MTX (IIT + REAL); keep panel genes + n_genes + total UMI.

Author metadata is joined by barcode / cellID. Does not load the full sparse matrix
into RAM (Real MTX is too large for a 16 GB box).
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pathlib import Path as _P
import sys as _sys
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from gene_sets import all_panel_genes

DATA = Path("/tmp/scrna_cytotrace")


def feature_row_map(path: Path, wanted: set[str]) -> dict[int, str]:
    """1-based MTX row -> gene symbol for genes in wanted."""
    out: dict[int, str] = {}
    seen: set[str] = set()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            symbols = [p for p in parts[:2] if p]
            hit = next((s for s in symbols if s in wanted), None)
            if hit and hit not in seen:
                out[i] = hit
                seen.add(hit)
    return out


def stream_mtx_panel(mtx: Path, n_cells: int, row_to_gene: dict[int, str]) -> dict[str, np.ndarray]:
    """Column-major MTX: keep only panel gene rows. n_genes comes from Seurat meta."""
    expr = {g: np.zeros(n_cells, dtype=np.float32) for g in row_to_gene.values()}
    wanted = row_to_gene
    n_hit = 0
    n_nz = 0
    with gzip.open(mtx, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            n_rows, n_cols, n_expected = (int(x) for x in line.split())
            print(f"  MTX {n_rows} x {n_cols} nnz={n_expected} panel_rows={len(wanted)}", flush=True)
            if n_cols != n_cells:
                print(f"  WARNING: barcodes={n_cells} MTX cols={n_cols}", flush=True)
            break
        for line in f:
            n_nz += 1
            r_s, c_s, v_s = line.split()
            gene = wanted.get(int(r_s))
            if gene is not None:
                c = int(c_s) - 1
                if 0 <= c < n_cells:
                    expr[gene][c] = float(v_s)
                    n_hit += 1
            if n_nz % 50_000_000 == 0:
                print(f"  nnz={n_nz} panel_hits={n_hit}", flush=True)
    print(f"  finished nnz={n_nz} panel_hits={n_hit}", flush=True)
    return expr


def align_meta(meta: pd.DataFrame, barcodes: list[str]) -> pd.DataFrame:
    md = meta.copy()
    if "cellID" in md.columns:
        md = md.set_index("cellID", drop=False)
    else:
        md = md.set_index(md.columns[0], drop=False)
    if barcodes[0] in md.index:
        return md.reindex(barcodes)
    # IIT: barcode AAACCTGAGCTATGCT-1 ; cellID P343_AAACCTGAGCTATGCT-1
    if "sampleID" in md.columns:
        pref = md["sampleID"].astype(str) + "_" + pd.Index(barcodes)
        # try constructing cellIDs from equal-length row order
        if len(md) == len(barcodes):
            print("  align by equal-length row order", flush=True)
            out = md.reset_index(drop=True)
            out["barcode"] = barcodes
            return out
    # strip sample prefix
    stripped = md.index.astype(str).str.replace(r"^[^_]+_", "", regex=True)
    if stripped.isin(barcodes).mean() > 0.8:
        md = md.set_index(stripped)
        return md.reindex(barcodes)
    if len(md) == len(barcodes):
        print("  fallback equal-length row order", flush=True)
        out = md.reset_index(drop=True)
        out["barcode"] = barcodes
        return out
    raise RuntimeError(f"cannot align meta {len(md)} vs barcodes {len(barcodes)}")


def extract_cohort(name: str, mtx: Path, features: Path, barcodes: Path, meta: Path, out: Path) -> None:
    print(f"=== {name} ===", flush=True)
    bc = pd.read_csv(barcodes, sep="\t", header=None)
    barcodes_list = bc.iloc[:, 0].astype(str).tolist()
    n_cells = len(barcodes_list)
    md = pd.read_csv(meta, sep="\t", dtype=str, low_memory=False)
    print(f"  meta={md.shape} barcodes={n_cells} cols={list(md.columns)}", flush=True)
    wanted = set(all_panel_genes())
    row_map = feature_row_map(features, wanted)
    print(f"  panel hits={len(row_map)} missing={sorted(wanted - set(row_map.values()))}", flush=True)
    expr = stream_mtx_panel(mtx, n_cells, row_map)
    aligned = align_meta(md, barcodes_list)
    aligned = aligned.reset_index(drop=True)
    aligned["barcode"] = barcodes_list
    # Author Seurat QC columns = same n_genes / depth CytoTRACE uses
    if "nCount_RNA" in aligned.columns:
        aligned["total_umi"] = pd.to_numeric(aligned["nCount_RNA"], errors="coerce")
    else:
        raise KeyError("nCount_RNA missing from author meta")
    if "nFeature_RNA" in aligned.columns:
        aligned["n_genes"] = pd.to_numeric(aligned["nFeature_RNA"], errors="coerce")
    else:
        raise KeyError("nFeature_RNA missing from author meta")
    aligned["cohort"] = name
    for g, arr in expr.items():
        aligned[g] = arr.astype(np.int32)
    aligned.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"  wrote {out} rows={len(aligned)}", flush=True)


def main() -> None:
    extract_cohort(
        "IIT_EGFRmut",
        DATA / "GSE241934_IIT_Matrix.mtx.gz",
        DATA / "GSE241934_IIT_features.tsv.gz",
        DATA / "GSE241934_IIT_barcodes.tsv.gz",
        DATA / "GSE241934_IIT_Meta.txt.gz",
        DATA / "gse241934_iit_panel_cells.tsv.gz",
    )
    extract_cohort(
        "REAL_WT",
        DATA / "GSE241934_Real_Matrix.mtx.gz",
        DATA / "GSE241934_RWC_features.tsv.gz",
        DATA / "GSE241934_RWC_barcodes.tsv.gz",
        DATA / "GSE241934_Real_Meta.txt.gz",
        DATA / "gse241934_real_panel_cells.tsv.gz",
    )


if __name__ == "__main__":
    sys.exit(main())
