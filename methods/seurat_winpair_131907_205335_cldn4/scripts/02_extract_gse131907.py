#!/usr/bin/env python3
"""Stream GSE131907 UMI TSV to Matrix Market for the capped keep list.

Extract-only. Primary analysis is Seurat in 03_seurat_integrate.R.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

ROOT = Path(__file__).resolve().parents[1]
KEEP = ROOT / "extract" / "keep_GSE131907.tsv"
UMI = Path("/tmp/winpair_geo/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz")
OUT = ROOT / "extract" / "GSE131907"


def main() -> int:
    keep_ids: list[str] = []
    with KEEP.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bcol = header.index("barcode")
        for line in fh:
            keep_ids.append(line.rstrip("\n").split("\t")[bcol])
    keep_set = set(keep_ids)
    print(f"keep barcodes={len(keep_ids)} unique={len(keep_set)}", flush=True)

    with gzip.open(UMI, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE131907 cell IDs in UMI header")
        ordered = [cell_ids[i] for i in col_idx]
        print(f"header cells={len(cell_ids)} keep={len(col_idx)}", flush=True)
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz_total = 0
        for gi, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            vals = np.fromstring(raw[tab0 + 1 :], sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(f"row {gi} {gene}: {vals.size} != {len(cell_ids)}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 2000 == 0:
                print(f"  streamed {gi} genes nnz={nnz_total}", flush=True)

    if data:
        data_a = np.concatenate(data)
        indices_a = np.concatenate(indices)
    else:
        data_a = np.array([], dtype=np.float32)
        indices_a = np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, indices_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered)),
        dtype=np.float32,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    mmwrite(OUT / "matrix.mtx", mat)
    (OUT / "genes.tsv").write_text("\n".join(genes) + "\n")
    (OUT / "barcodes.tsv").write_text("\n".join(ordered) + "\n")
    print(f"wrote {OUT} genes={len(genes)} cells={len(ordered)} nnz={mat.nnz}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
