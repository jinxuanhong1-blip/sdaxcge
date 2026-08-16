#!/usr/bin/env python3
"""Assemble per-cell gene-panel UMI and metadata from the streamed RDS dump."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PANEL_META_KEYS = (
    "orig.ident",
    "orig_ident",
    "sample",
    "sample_id",
    "patient",
    "donor",
    "tissue",
    "group",
    "cell_type",
    "celltype",
    "CellType",
    "garnett_cluster",
    "garnett_cell_type",
    "garnett_celltype",
    "cell_class",
    "seurat_clusters",
    "ident",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "nCount_SCT",
    "nFeature_SCT",
)


def parse_patient(text: str) -> str:
    m = re.search(r"(MRC0*\d+)", str(text).replace(" ", ""))
    return m.group(1) if m else "unknown"


def parse_tissue(text: str) -> str:
    t = str(text).upper()
    if "ANT" in t or "ADJACENT" in t:
        return "ANT"
    if "CD45" in t:
        return "ANT" if "ANT" in t else "Tumor"
    return "Tumor"


def pick_csc(arrays: list[dict], n_genes: int) -> dict | None:
    by_len = {}
    for a in arrays:
        by_len.setdefault(a["length"], []).append(a)
    # p is n_cells+1; i and x share nnz
    candidates = []
    for a in arrays:
        if a["dtype"] != "int32":
            continue
        # possible p
        n_cells = a["length"] - 1
        if n_cells < 1000:
            continue
        nnz_groups = [k for k, v in by_len.items() if k > n_cells and len(v) >= 1]
        for nnz in nnz_groups:
            ints = [x for x in by_len[nnz] if x["dtype"] == "int32"]
            reals = [x for x in by_len[nnz] if x["dtype"] == "float64"]
            if ints and reals:
                candidates.append(
                    {
                        "n_genes": n_genes,
                        "n_cells": n_cells,
                        "nnz": nnz,
                        "p": a,
                        "i": ints[0],
                        "x": reals[0],
                    }
                )
    if not candidates:
        return None
    # Prefer the largest nnz (counts usually denser than a subset assay)
    candidates.sort(key=lambda c: c["nnz"], reverse=True)
    return candidates[0]


def extract_rows(i_file: Path, p_file: Path, x_file: Path, gene_idx: dict[str, int], n_cells: int):
    i = np.load(i_file, mmap_mode="r")
    p = np.load(p_file, mmap_mode="r")
    x = np.load(x_file, mmap_mode="r")
    want = {idx: name for name, idx in gene_idx.items()}
    out = {name: np.zeros(n_cells, dtype=np.float32) for name in gene_idx}
    total = np.zeros(n_cells, dtype=np.float64)
    if len(p) != n_cells + 1:
        n_cells = len(p) - 1
        out = {name: np.zeros(n_cells, dtype=np.float32) for name in gene_idx}
        total = np.zeros(n_cells, dtype=np.float64)
    for col in range(n_cells):
        a, b = int(p[col]), int(p[col + 1])
        if b <= a:
            continue
        rows = np.asarray(i[a:b])
        vals = np.asarray(x[a:b], dtype=np.float64)
        total[col] = float(vals.sum())
        for row, val in zip(rows.tolist(), vals.tolist()):
            name = want.get(int(row))
            if name is not None:
                out[name][col] = float(val)
        if col and col % 50000 == 0:
            print(f"CSC {col}/{n_cells}", flush=True)
    out["total"] = total.astype(np.float32)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extracted", type=Path, default=Path("data/gse253013/extracted"))
    args = ap.parse_args()
    ex = args.extracted
    arrays = json.loads((ex / "arrays.json").read_text())
    gene_info = json.loads((ex / "gene_index.json").read_text())
    str_index = json.loads((ex / "str_index.json").read_text())
    gene_idx = gene_info.get("index") or {}
    n_genes = gene_info.get("n_genes_in_vector")
    if not gene_idx:
        raise SystemExit("no target genes in gene_index.json")

    csc = pick_csc(arrays, n_genes or 0)
    if csc is None:
        raise SystemExit(f"could not pair i/p/x arrays; see {ex/'arrays.json'}")
    print("CSC", {k: csc[k] if k in ("n_genes", "n_cells", "nnz") else csc[k]["file"] for k in csc})

    expr = extract_rows(
        Path(csc["i"]["file"]),
        Path(csc["p"]["file"]),
        Path(csc["x"]["file"]),
        {k: int(v) for k, v in gene_idx.items()},
        csc["n_cells"],
    )
    np.savez_compressed(ex / "gene_panel.npz", **expr)

    n_cells = len(next(iter(expr.values())))
    meta = pd.DataFrame({"cell_index": np.arange(n_cells)})

    # Attach string vectors whose length matches n_cells.
    used = []
    for rec in str_index:
        path = rec.get("file")
        if not path or not Path(path).exists():
            continue
        strings = Path(path).read_text().splitlines()
        if len(strings) != n_cells:
            continue
        col = rec.get("path") or Path(path).stem
        col = re.sub(r"[^A-Za-z0-9_.]+", "_", col)[-40:] or "meta"
        if col in meta.columns:
            col = col + "_2"
        meta[col] = strings
        used.append({"path": rec.get("path"), "col": col})

    # Patient / tissue from the most barcode-like or orig.ident-like column
    src_col = None
    for c in meta.columns:
        if meta[c].astype(str).str.contains("MRC", regex=False).mean() > 0.5:
            src_col = c
            break
    if src_col is None:
        for c in meta.columns:
            if c != "cell_index":
                src_col = c
                break
    if src_col:
        meta["orig.ident"] = meta[src_col]
        meta["patient"] = [parse_patient(x) for x in meta[src_col]]
        meta["tissue"] = [parse_tissue(x) for x in meta[src_col]]

    meta.to_csv(ex / "cell_metadata.tsv", sep="\t", index=False)
    (ex / "assemble_summary.json").write_text(
        json.dumps(
            {
                "n_cells": n_cells,
                "n_genes_extracted": len(gene_idx),
                "csc": {
                    "n_cells": csc["n_cells"],
                    "nnz": csc["nnz"],
                    "i": csc["i"]["file"],
                    "p": csc["p"]["file"],
                    "x": csc["x"]["file"],
                },
                "metadata_columns": used,
                "patients": meta["patient"].value_counts().to_dict() if "patient" in meta else {},
            },
            indent=2,
        )
        + "\n"
    )
    print("wrote", ex / "gene_panel.npz", "n_cells", n_cells)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
