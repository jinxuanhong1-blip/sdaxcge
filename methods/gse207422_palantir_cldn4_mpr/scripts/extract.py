#!/usr/bin/env python3
"""Two-pass extract of GSE207422 A3-malignant-like epithelium.

Pass 1: stream markers + library size; classify epithelial / A3-malignant-like.
Pass 2: stream the full UMI matrix and keep only A3-malignant-like cells as
a sparse cells×genes h5ad (counts layer). Uncompressed text is never written.

A3-malignant-like = epithelial (argmax lineage) AND zero UMI for
SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. Author CopyKAT barcodes are not on GEO.
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
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import A3_NORMAL_LUNG, LINEAGES, all_panel_genes  # noqa: E402

try:
    import anndata as ad
except ImportError as e:
    raise SystemExit("anndata is required") from e


PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}


def score_log1p(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score_log1p(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def stream_pass1(matrix: Path) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray, int]:
    wanted = set(all_panel_genes())
    t0 = time.time()
    rows: dict[str, np.ndarray] = {}
    ngenes = 0
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        ncells = len(cells)
        total = np.zeros(ncells, dtype=np.int64)
        n_nonzero_genes = np.zeros(ncells, dtype=np.int32)
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != ncells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {ncells}")
            total += arr
            n_nonzero_genes += arr > 0
            if gene in wanted:
                rows[gene] = arr
            if ngenes % 4000 == 0:
                print(
                    f"  pass1 {ngenes} genes, {len(rows)}/{len(wanted)} markers, {time.time() - t0:.0f}s",
                    flush=True,
                )
    print(
        f"pass1 done: {ngenes} genes × {ncells} cells in {time.time() - t0:.0f}s; "
        f"markers {len(rows)}/{len(wanted)}",
        flush=True,
    )
    return cells, rows, total, n_nonzero_genes, ngenes


def stream_pass2(matrix: Path, keep_idx: np.ndarray) -> tuple[list[str], csr_matrix]:
    t0 = time.time()
    n_keep = int(len(keep_idx))
    indptr = [0]
    indices: list[int] = []
    data: list[int] = []
    genes: list[str] = []
    ngenes = 0
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ncells = len(header) - 1
        for line in fh:
            ngenes += 1
            i = line.find("\t")
            gene = line[:i]
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != ncells:
                raise SystemExit(f"{gene}: {arr.size} values, expected {ncells}")
            sub = arr[keep_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                indices.extend(nz.tolist())
                data.extend(sub[nz].tolist())
            indptr.append(len(indices))
            genes.append(gene)
            if ngenes % 4000 == 0:
                print(
                    f"  pass2 {ngenes} genes, nnz={len(data)}, {time.time() - t0:.0f}s",
                    flush=True,
                )
    X = csr_matrix(
        (np.asarray(data, dtype=np.int32), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
        shape=(ngenes, n_keep),
        dtype=np.int32,
    )
    X = X.T.tocsr()
    print(
        f"pass2 done: {X.shape[0]} cells × {X.shape[1]} genes, nnz={X.nnz}, {time.time() - t0:.0f}s",
        flush=True,
    )
    return genes, X


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("data/GSE207422/a3_malignant_like.h5ad"),
    )
    args = ap.parse_args()
    matrix = args.workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_xlsx = args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    if not matrix.exists():
        raise SystemExit(f"missing {matrix}; run download.py")

    cells, expr, total, n_nonzero_genes, ngenes = stream_pass1(matrix)
    n = len(cells)
    lineage = assign_lineage(expr, n)
    is_epi = lineage == "epithelial"
    normal_umi = np.zeros(n, dtype=np.int64)
    present_normal = []
    for g in A3_NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g].astype(np.int64)
            present_normal.append(g)
    is_malig = is_epi & (normal_umi == 0)
    keep_idx = np.flatnonzero(is_malig)
    print(
        f"lineage epithelial={int(is_epi.sum())} A3-malignant-like={int(is_malig.sum())} "
        f"(normal panel present={present_normal})",
        flush=True,
    )

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells], dtype=object)
    meta = load_sample_meta(meta_xlsx) if meta_xlsx.exists() else pd.DataFrame()
    meta_map = meta.set_index("Sample") if len(meta) else None

    class_rows = []
    for sid in sorted(set(sample_ids.tolist())):
        m = sample_ids == sid
        class_rows.append(
            {
                "Sample": sid,
                "paper_group": PAPER_GROUP.get(sid, "NA"),
                "n_cells": int(m.sum()),
                "n_epithelial": int((m & is_epi).sum()),
                "n_malig_a3": int((m & is_malig).sum()),
            }
        )
    class_df = pd.DataFrame(class_rows)
    class_path = args.workdir / "classification_counts.tsv"
    class_df.to_csv(class_path, sep="\t", index=False)

    genes, X = stream_pass2(matrix, keep_idx)
    obs = pd.DataFrame(index=pd.Index(cells[keep_idx].astype(str), name="cell"))
    obs["barcode"] = cells[keep_idx].astype(str)
    obs["Sample"] = sample_ids[keep_idx]
    obs["paper_group"] = obs["Sample"].map(PAPER_GROUP).fillna("NA")
    obs["timing"] = np.where(obs["paper_group"].eq("TN"), "pre", "post")
    obs["total_umi"] = total[keep_idx]
    obs["n_genes"] = n_nonzero_genes[keep_idx]
    obs["lineage"] = lineage[keep_idx]
    if meta_map is not None:
        for col in ("path_response", "Resource"):
            if col in meta_map.columns:
                obs[col] = obs["Sample"].map(meta_map[col])

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.layers["counts"] = adata.X.copy()
    adata.uns["extract"] = {
        "dataset": "GSE207422",
        "definition": "A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3",
        "n_cells_matrix": int(n),
        "n_genes_matrix": int(ngenes),
        "n_epithelial": int(is_epi.sum()),
        "n_malig_a3": int(is_malig.sum()),
        "normal_panel_present": present_normal,
        "author_copykat_public": False,
        "dual_high_gate": False,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")
    print(f"wrote {args.out} {adata.shape}", flush=True)
    print(json.dumps(adata.uns["extract"], indent=2), flush=True)


if __name__ == "__main__":
    main()
