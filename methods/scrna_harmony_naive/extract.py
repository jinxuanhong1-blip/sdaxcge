#!/usr/bin/env python3
"""Extract the shared gene panel from each public naive tumor atlas.

GSE253013 uses the XDR RDS walker in extract_gse253013.py (the 9.3 GB
Seurat object does not fit in 16 GB RAM).
"""

from __future__ import annotations

import argparse
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from config import (
    DATA,
    EXTRACTED,
    FILES,
    GENE_PANEL,
    GSE131907_TUMOR_ORIGINS,
)


def _save_dataset(name: str, expr: dict[str, np.ndarray], meta: pd.DataFrame, present: list[str]) -> None:
    out = EXTRACTED / name
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "gene_panel.npz", **expr)
    meta.to_csv(out / "cell_metadata.tsv", sep="\t", index=False)
    (out / "gene_index.json").write_text(
        json.dumps(
            {
                "dataset": name,
                "present": present,
                "absent": [g for g in GENE_PANEL if g not in present],
                "n_cells": int(len(meta)),
                "TACSTD2": "TACSTD2" in present,
                "CLDN4": "CLDN4" in present,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"{name}: n={len(meta)} genes={present}", flush=True)


def stream_gene_rows(path: Path, wanted: set[str], sep: str = "\t") -> tuple[list[str], dict[str, np.ndarray]]:
    """Read a genes-by-cells text matrix; keep only wanted gene rows."""
    opener = gzip.open if str(path).endswith(".gz") else open
    found: dict[str, np.ndarray] = {}
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(sep)
        # GSE131907: first token is "Index". GSE148071: first token is a barcode.
        if header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n_cells = len(cell_ids)
        for line in handle:
            gene, _, rest = line.partition(sep)
            gene = gene.strip().strip('"').split(".")[0]
            if gene not in wanted:
                continue
            arr = np.fromstring(rest, sep=sep, dtype=np.float32)
            if arr.size != n_cells:
                # retry if the first field was not skipped correctly
                toks = line.rstrip("\n").split(sep)
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n_cells}")
            found[gene] = arr
            if len(found) == len(wanted):
                break
    return cell_ids, found


def extract_gse131907() -> None:
    annot = pd.read_csv(FILES["GSE131907_annot"], sep="\t")
    cell_ids, expr = stream_gene_rows(FILES["GSE131907_umi"], set(GENE_PANEL))
    # Align annotation to matrix columns
    annot = annot.set_index("Index", drop=False)
    missing = [c for c in cell_ids if c not in annot.index]
    if missing:
        # try Barcode_Sample
        annot["alt"] = annot["Barcode"].astype(str) + "_" + annot["Sample"].astype(str)
        alt = annot.set_index("alt", drop=False)
        if sum(c in alt.index for c in cell_ids[:100]) > 50:
            annot = alt
        else:
            raise SystemExit(f"GSE131907 barcode mismatch, e.g. {missing[:3]}")
    order = [c for c in cell_ids]
    annot = annot.loc[order].reset_index(drop=True)
    keep = annot["Sample_Origin"].isin(GSE131907_TUMOR_ORIGINS).to_numpy()
    annot = annot.loc[keep].copy()
    expr = {g: v[keep] for g, v in expr.items()}
    annot["dataset"] = "GSE131907"
    annot["donor"] = annot["Sample"].astype(str)
    annot["patient"] = annot["Sample"].astype(str)
    annot["tissue"] = "Tumor"
    annot["site"] = annot["Sample_Origin"].astype(str)
    annot["author_cell_type"] = annot["Cell_type"].astype(str)
    annot["author_subtype"] = annot["Cell_subtype"].astype(str)
    present = sorted(expr)
    lib = None
    # library size not in the gene panel; use n/a and later CP10k from panel sum
    _save_dataset("GSE131907", expr, annot, present)


def extract_gse148071() -> None:
    tar_path = FILES["GSE148071_tar"]
    untar = DATA / "GSE148071_files"
    untar.mkdir(parents=True, exist_ok=True)
    if not any(untar.glob("*_exp.txt.gz")):
        with tarfile.open(tar_path, "r") as tf:
            tf.extractall(untar)
    files = sorted(untar.glob("*_exp.txt.gz"))
    if not files:
        files = sorted(untar.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no GSE148071 exp files in {untar}")

    expr_all: dict[str, list[np.ndarray]] = {g: [] for g in GENE_PANEL}
    rows = []
    present_any: set[str] = set()
    for fp in files:
        # GSM4453576_P1_exp.txt.gz
        stem = fp.name.replace(".txt.gz", "")
        patient = stem.split("_")[1] if "_" in stem else stem
        cell_ids, found = stream_gene_rows(fp, set(GENE_PANEL))
        n = len(cell_ids)
        present_any.update(found)
        for g in GENE_PANEL:
            expr_all[g].append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        rec = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "dataset": "GSE148071",
                "donor": patient,
                "patient": patient,
                "tissue": "Tumor",
                "site": "tumor_biopsy",
                "author_cell_type": "",
                "author_subtype": "",
            }
        )
        rows.append(rec)
        print(f"  {fp.name} n={n} genes={sorted(found)}", flush=True)

    expr = {g: np.concatenate(v) for g, v in expr_all.items() if any(x.sum() for x in v) or g in present_any}
    # drop all-zero genes that were never present
    expr = {g: v for g, v in expr.items() if g in present_any}
    meta = pd.concat(rows, ignore_index=True)
    _save_dataset("GSE148071", expr, meta, sorted(present_any))


def extract_gse127465() -> None:
    genes = pd.read_csv(FILES["GSE127465_genes"], sep="\t", header=None)[0].astype(str)
    want_idx = {i: g for i, g in enumerate(genes) if g in set(GENE_PANEL)}
    meta = pd.read_csv(FILES["GSE127465_meta"], sep="\t")
    n_cells = len(meta)
    expr = {g: np.zeros(n_cells, dtype=np.float32) for g in want_idx.values()}

    # MatrixMarket: filename is 54773 x 41861 → rows=cells, cols=genes (1-based)
    opener = gzip.open if str(FILES["GSE127465_mtx"]).endswith(".gz") else open
    with opener(FILES["GSE127465_mtx"], "rt") as handle:
        header = handle.readline()
        while header.startswith("%"):
            header = handle.readline()
        nrow, ncol, nnz = [int(x) for x in header.split()[:3]]
        print(f"GSE127465 MTX {nrow} x {ncol} nnz={nnz}", flush=True)
        # decide orientation
        if nrow == n_cells and ncol == len(genes):
            cell_dim, gene_dim = 0, 1
        elif ncol == n_cells and nrow == len(genes):
            cell_dim, gene_dim = 1, 0
        else:
            raise SystemExit(f"MTX shape {nrow}x{ncol} vs cells={n_cells} genes={len(genes)}")
        # 1-based indices
        for i, line in enumerate(handle):
            a, b, val = line.split()[:3]
            r, c = int(a) - 1, int(b) - 1
            gene_i = r if gene_dim == 0 else c
            cell_i = r if cell_dim == 0 else c
            g = want_idx.get(gene_i)
            if g is not None:
                expr[g][cell_i] = float(val)
            if i and i % 5_000_000 == 0:
                print(f"  MTX {i}/{nnz}", flush=True)

    keep = meta["Tissue"].astype(str).str.lower().eq("tumor").to_numpy()
    meta = meta.loc[keep].copy()
    expr = {g: v[keep] for g, v in expr.items()}
    meta["dataset"] = "GSE127465"
    meta["donor"] = meta["Patient"].astype(str)
    meta["patient"] = meta["Patient"].astype(str)
    meta["tissue"] = "Tumor"
    meta["site"] = "tumor"
    meta["author_cell_type"] = meta["Major cell type"].astype(str)
    meta["author_subtype"] = meta["Minor subset"].astype(str)
    present = [g for g, v in expr.items() if float(np.max(v)) > 0 or g in set(want_idx.values())]
    _save_dataset("GSE127465", expr, meta.reset_index(drop=True), sorted(present))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["GSE131907", "GSE148071", "GSE127465"],
        help="text/mtx datasets (GSE253013 is a separate RDS walker)",
    )
    args = ap.parse_args()
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    if "GSE131907" in args.datasets:
        extract_gse131907()
    if "GSE148071" in args.datasets:
        extract_gse148071()
    if "GSE127465" in args.datasets:
        extract_gse127465()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
