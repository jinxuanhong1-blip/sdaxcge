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
    """Prefer RNA @counts (assays/[0]/counts/i), not SCT or @data."""
    by_path = {a["path"]: a for a in arrays}
    prefer = [
        ("assays/[0]/counts/i", "assays/[0]/counts/i/p", "assays/[0]/counts/i/p/Dim/Dimnames/x"),
        ("assays/[0]/counts/data/i", "assays/[0]/counts/data/i/p", "assays/[0]/counts/data/i/p/Dim/Dimnames/x"),
    ]
    for ip, pp, xp in prefer:
        if ip in by_path and pp in by_path and xp in by_path:
            p = by_path[pp]
            return {
                "n_genes": n_genes,
                "n_cells": p["length"] - 1,
                "nnz": by_path[ip]["length"],
                "p": p,
                "i": by_path[ip],
                "x": by_path[xp],
            }
    return None


def extract_rows(i_file: Path, p_file: Path, x_file: Path, gene_idx: dict[str, int], n_cells: int):
    i = np.load(i_file, mmap_mode="r")
    p = np.load(p_file, mmap_mode="r")
    x = np.load(x_file, mmap_mode="r")
    n_cells = len(p) - 1
    max_row = int(max(gene_idx.values())) if gene_idx else 0
    inv = np.full(max_row + 1, -1, dtype=np.int32)
    names = list(gene_idx)
    for slot, (name, idx) in enumerate(gene_idx.items()):
        if idx >= 0:
            inv[int(idx)] = slot
    out_mat = np.zeros((len(names), n_cells), dtype=np.float32)
    total = np.zeros(n_cells, dtype=np.float64)
    nnz = len(i)
    chunk = 4_000_000
    for start in range(0, nnz, chunk):
        stop = min(start + chunk, nnz)
        rows = np.asarray(i[start:stop])
        vals = np.asarray(x[start:stop], dtype=np.float64)
        # library size from this counts matrix
        cols_all = np.searchsorted(p, np.arange(start, stop, dtype=np.int64), side="right") - 1
        np.add.at(total, cols_all, vals)
        mapped = np.where((rows >= 0) & (rows <= max_row), inv[np.clip(rows, 0, max_row)], -1)
        hit = mapped >= 0
        if not np.any(hit):
            print(f"CSC {stop}/{nnz}", flush=True)
            continue
        out_mat[mapped[hit], cols_all[hit]] = vals[hit].astype(np.float32)
        print(f"CSC {stop}/{nnz}", flush=True)
    out = {name: out_mat[j] for j, name in enumerate(names)}
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
    used = []

    def add_str(path_suffix: str, col: str) -> None:
        for rec in str_index:
            if rec.get("length") != n_cells:
                continue
            if rec.get("path") == path_suffix or (rec.get("file") and path_suffix in str(rec.get("path"))):
                fp = rec.get("file")
                if fp and Path(fp).exists():
                    vals = Path(fp).read_text().splitlines()
                    if len(vals) == n_cells:
                        meta[col] = vals
                        used.append({"path": rec.get("path"), "col": col})
                        return

    # Author Seurat meta.data columns recovered from the RDS walk
    add_str("assays/[1]/counts/data/scale.data/key/var.features/meta.features/misc/[0]/[9]/[3]", "library_id")
    add_str("assays/[1]/counts/data/scale.data/key/var.features/meta.features/misc/[0]/[9]/[4]", "patient")
    add_str("assays/[1]/counts/data/scale.data/key/var.features/meta.features/misc/[0]/[9]/[5]", "tissue_raw")
    add_str("assays/[1]/counts/data/scale.data/key/var.features/meta.features/misc/[0]/[9]/[15]", "treatment")

    # Fallback: first 256379 STR with MRC* and T/NAT
    if "patient" not in meta.columns:
        for rec in str_index:
            fp = rec.get("file")
            if not fp or not Path(fp).exists():
                continue
            vals = Path(fp).read_text().splitlines()
            if len(vals) != n_cells:
                continue
            if sum(v.startswith("MRC") for v in vals[:200]) > 50:
                meta["patient"] = vals
                used.append({"path": rec.get("path"), "col": "patient"})
            if set(vals[:5000]).issubset({"T", "NAT", "ANT"}):
                meta["tissue_raw"] = vals
                used.append({"path": rec.get("path"), "col": "tissue_raw"})

    if "patient" in meta.columns:
        meta["orig.ident"] = meta["patient"]
        meta["patient"] = [parse_patient(x) for x in meta["patient"]]
    if "tissue_raw" in meta.columns:
        meta["tissue"] = [
            "ANT" if str(x).upper() in {"NAT", "ANT", "NORMAL"} else "Tumor" for x in meta["tissue_raw"]
        ]

    # Author coarse cell types: active.ident factor (1-based) + 10 levels
    levels_file = ex / "str_092_n10.txt"
    ident_arr = ex / "tmp_arrays" / "arr_0038_int32.npy"
    if levels_file.exists() and ident_arr.exists():
        levels = levels_file.read_text().splitlines()
        codes = np.load(ident_arr)
        mapped = []
        for c in codes:
            if 1 <= int(c) <= len(levels):
                mapped.append(levels[int(c) - 1])
            else:
                mapped.append("Unknown")
        meta["cell_type"] = mapped
        used.append({"path": "active.ident", "col": "cell_type"})

    ncount = ex / "tmp_arrays" / "arr_0014_f64.npy"
    nfeat = ex / "tmp_arrays" / "arr_0015_int32.npy"
    if ncount.exists() and np.load(ncount).shape[0] == n_cells:
        meta["nCount_RNA"] = np.load(ncount)
        used.append({"path": "meta.data nCount_RNA", "col": "nCount_RNA"})
    if nfeat.exists() and np.load(nfeat).shape[0] == n_cells:
        meta["nFeature_RNA"] = np.load(nfeat)
        used.append({"path": "meta.data nFeature_RNA", "col": "nFeature_RNA"})

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
