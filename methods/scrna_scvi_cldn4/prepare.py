#!/usr/bin/env python3
"""Stream GSE207422 genes×cells UMI TSV.gz to sparse h5ad. Dense text is never written."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def stream_umi(path: Path) -> ad.AnnData:
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n_cells = len(cells)
        print(f"GSE207422 cells={n_cells}", flush=True)
        genes: list[str] = []
        data_chunks: list[np.ndarray] = []
        idx_chunks: list[np.ndarray] = []
        indptr = [0]
        for i, line in enumerate(fh, start=1):
            tab = line.find("\t")
            gene = line[:tab]
            arr = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.int32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} != {n_cells}")
            nz = np.flatnonzero(arr)
            if nz.size:
                data_chunks.append(arr[nz])
                idx_chunks.append(nz.astype(np.int32, copy=False))
                indptr.append(indptr[-1] + int(nz.size))
            else:
                indptr.append(indptr[-1])
            genes.append(gene)
            if i % 4000 == 0:
                print(f"  genes={i} nnz={indptr[-1]}", flush=True)
    data = np.concatenate(data_chunks) if data_chunks else np.array([], dtype=np.int32)
    indices = np.concatenate(idx_chunks) if idx_chunks else np.array([], dtype=np.int32)
    X_genes = csr_matrix(
        (data, indices, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), n_cells),
        dtype=np.int32,
    )
    X = X_genes.T.tocsr()
    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    obs["barcode"] = cells
    obs["sample"] = [c.rsplit("_", 1)[0] for c in cells]
    obs["dataset"] = "GSE207422"
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    print(f"GSE207422 matrix {X.shape} nnz={X.nnz}", flush=True)
    return ad.AnnData(X=X, obs=obs, var=var)


def attach_meta(adata: ad.AnnData, xlsx: Path) -> None:
    meta = pd.read_excel(xlsx)
    meta = meta.dropna(subset=["Sample"]).copy()
    meta = meta[~meta["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(columns={"Sample": "sample"})
    meta["patient"] = meta["Patient"].astype(str)
    meta["timing"] = np.where(
        meta["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    meta["mpr"] = meta["Pathologic Response"].map(
        {"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR", "NE": "NE"}
    )
    meta["histology"] = meta["Pathology"].astype(str)
    keep = [
        "sample",
        "patient",
        "timing",
        "mpr",
        "histology",
        "Resource",
        "Pathologic Response",
        "RECIST",
        "PD1 Antibody",
    ]
    meta = meta[keep].drop_duplicates("sample")
    adata.obs = adata.obs.merge(meta, on="sample", how="left")
    adata.obs.index = adata.obs["barcode"].astype(str)
    adata.obs.index.name = "cell"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/scrna_scvi_cldn4/raw"))
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/scrna_scvi_cldn4/h5ad"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    umi = args.raw / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    xlsx = args.raw / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    adata = stream_umi(umi)
    attach_meta(adata, xlsx)
    out = args.outdir / "GSE207422.h5ad"
    adata.write_h5ad(out, compression="gzip")
    print(f"wrote {out} {adata.n_obs} x {adata.n_vars}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
