#!/usr/bin/env python3
"""TISCH h5 + tissue-filter helpers."""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from config import EXCLUDE_TISSUE_TOKENS, TUMOR_TOKENS


def _as_str_array(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.array(out, dtype=object)


def load_tisch_genes(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, list[str], dict]:
    audit = {"path": str(path), "keys": [], "shape": None, "n_genes": None, "n_cells": None}
    with h5py.File(path, "r") as f:
        audit["keys"] = list(f.keys())
        grp = None
        if all(k in f for k in ("data", "indices", "indptr")):
            grp = f
        else:
            for key in f.keys():
                g = f[key]
                if isinstance(g, h5py.Group) and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix in {path}: {audit['keys']}")

        def _names(candidates):
            for c in candidates:
                node = grp.get(c) if hasattr(grp, "get") else None
                if node is None:
                    node = f.get(c)
                if node is None:
                    continue
                if isinstance(node, h5py.Dataset):
                    return _as_str_array(node[:])
                if isinstance(node, h5py.Group):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return _as_str_array(node[sub][:])
            return None

        genes = _names(["features", "gene_names", "genes", "rownames"])
        barcodes = _names(["barcodes", "cell_names", "colnames"])
        shape = tuple(int(x) for x in grp["shape"][:]) if "shape" in grp else None
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]

    if shape is None:
        raise ValueError(f"{path}: missing shape")
    n0, n1 = int(shape[0]), int(shape[1])
    if genes is not None and barcodes is not None and len(genes) == n1 and len(barcodes) == n0:
        n_genes, n_cells = n1, n0
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))

        def col(i: int) -> np.ndarray:
            return np.asarray(csc[:, i].todense()).ravel()

    else:
        n_genes, n_cells = n0, n1
        csc = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))

        def col(i: int) -> np.ndarray:
            return np.asarray(csc.getrow(i).todense()).ravel()

    if genes is None or barcodes is None:
        raise ValueError(f"{path}: missing gene or barcode names")
    if len(genes) != n_genes or len(barcodes) != n_cells:
        raise ValueError(
            f"{path}: dim mismatch genes={len(genes)} barcodes={len(barcodes)} "
            f"shape=({n_genes},{n_cells})"
        )
    audit["shape"] = [n_genes, n_cells]
    audit["n_genes"] = n_genes
    audit["n_cells"] = n_cells
    present = [g for g in wanted if g in set(genes)]
    audit["genes_present"] = present
    audit["genes_absent"] = [g for g in wanted if g not in set(genes)]
    gene_index = {g: int(np.where(genes == g)[0][0]) for g in present}
    table = {"barcode": barcodes}
    for g, i in gene_index.items():
        table[g] = col(i)
    return pd.DataFrame(table).set_index("barcode"), present, audit


SAMPLE_EXCLUDE_RE = re.compile(r"(^LUNG_N|Juxta|_Normal$|_NAT$|^D\d+N$)", re.IGNORECASE)
SAMPLE_FORCE_TUMOR_RE = re.compile(r"(^LUNG_T|_Tumor$|^D\d+T$)", re.IGNORECASE)


def is_tumor_tissue(val: str) -> bool:
    s = str(val).strip().lower()
    if not s or s in {"nan", "none", "na"}:
        return False
    if any(tok in s for tok in EXCLUDE_TISSUE_TOKENS):
        return False
    return any(tok in s for tok in TUMOR_TOKENS)


def keep_tumor_unit(sample_name: str | None, tissue_val: str | None) -> bool:
    if sample_name and SAMPLE_EXCLUDE_RE.search(str(sample_name)):
        return False
    if sample_name and SAMPLE_FORCE_TUMOR_RE.search(str(sample_name)):
        return True
    if tissue_val is None:
        return True
    return is_tumor_tissue(tissue_val)


def pick_unit_columns(meta: pd.DataFrame) -> tuple[str | None, str | None]:
    cols = {c.lower(): c for c in meta.columns}
    patient = cols.get("patient")
    sample = cols.get("sample")
    if sample is not None and meta[sample].nunique(dropna=True) >= 0.9 * len(meta):
        sample = None
    return patient, sample


def pick_tissue_col(meta: pd.DataFrame) -> str | None:
    for name in ("Source", "Tissue", "source", "tissue"):
        if name in meta.columns:
            return name
    return None
