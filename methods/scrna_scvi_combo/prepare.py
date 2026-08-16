#!/usr/bin/env python3
"""Convert public GEO matrices to sparse h5ad (counts, sample metadata)."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def stream_gse207422_txt(path: Path) -> ad.AnnData:
    """Genes x cells TSV.gz -> cells x genes CSR. Does not materialize dense."""
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
    obs["cohort"] = "GSE207422"
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    print(f"GSE207422 matrix {X.shape} nnz={X.nnz}", flush=True)
    return ad.AnnData(X=X, obs=obs, var=var)


def attach_gse207422_sample_meta(adata: ad.AnnData, xlsx: Path) -> None:
    meta = pd.read_excel(xlsx)
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(columns={"Sample": "sample"})
    meta["patient"] = meta["Patient"].astype(str)
    meta["timing"] = meta["Resource"].map(
        {
            "Pre-treatment biopsy": "pre",
            "Post-treatment surgery": "post",
        }
    )
    meta["mpr"] = meta["Pathologic Response"].map(
        {"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR", "NE": "NE"}
    )
    meta["histology"] = meta["Pathology"].astype(str)
    meta["tissue"] = "tumor"
    keep = [
        "sample",
        "patient",
        "timing",
        "mpr",
        "histology",
        "tissue",
        "Resource",
        "Pathologic Response",
        "RECIST",
        "PD1 Antibody",
    ]
    meta = meta[keep].drop_duplicates("sample")
    adata.obs = adata.obs.merge(meta, on="sample", how="left")
    adata.obs.index = adata.obs["barcode"].astype(str)
    adata.obs.index.name = "cell"


def read_10x_mtx(mtx: Path, barcodes: Path, features: Path) -> ad.AnnData:
    import scanpy as sc

    adata = sc.read_mtx(mtx).T
    bc = pd.read_csv(barcodes, header=None, sep="\t")[0].astype(str).tolist()
    feat = pd.read_csv(features, header=None, sep="\t")
    genes = feat[0].astype(str).tolist()
    if adata.n_obs != len(bc) or adata.n_vars != len(genes):
        raise ValueError(
            f"shape mismatch mtx={adata.shape} barcodes={len(bc)} features={len(genes)}"
        )
    adata.obs_names = pd.Index(bc, name="cell")
    # collapse duplicate symbols by summing
    adata.var_names = pd.Index(genes, name="gene")
    if adata.var_names.duplicated().any():
        adata = adata[:, ~adata.var_names.duplicated()].copy()
    if hasattr(adata.X, "astype"):
        adata.X = adata.X.tocsr().astype(np.int32)
    return adata


def attach_gse241934_meta(adata: ad.AnnData, meta_path: Path, cohort: str) -> None:
    meta = pd.read_csv(meta_path, sep="\t")
    if "cellID" not in meta.columns:
        raise SystemExit(f"{meta_path} missing cellID")
    meta = meta.drop_duplicates("cellID")
    meta["barcode"] = meta["cellID"].astype(str)
    meta["sample"] = meta["sampleID"].astype(str)
    meta["patient"] = meta["sampleID"].astype(str)
    meta["timing"] = "post"
    meta["mpr"] = meta["Pathological Response"].map(
        {"MPR": "MPR", "pCR": "MPR", "non-MPR": "NMPR", "NMPR": "NMPR"}
    )
    meta["histology"] = meta["Histology"].astype(str)
    meta["tissue"] = "tumor"
    meta["author_major"] = meta["major.cell.type"].astype(str)
    keep = [
        "barcode",
        "sample",
        "patient",
        "timing",
        "mpr",
        "histology",
        "tissue",
        "author_major",
        "Pathological Response",
        "PD1",
        "EGFR",
    ]
    meta = meta[[c for c in keep if c in meta.columns]]
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs = adata.obs.merge(meta, on="barcode", how="left")
    adata.obs["dataset"] = "GSE241934"
    adata.obs["cohort"] = cohort
    adata.obs.index = adata.obs["barcode"].astype(str)
    adata.obs.index.name = "cell"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/scrna_scvi_combo"))
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/scrna_scvi_combo/h5ad"))
    ap.add_argument("--include-rwc", action="store_true")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    p207 = args.raw / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    m207 = args.raw / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    a207 = stream_gse207422_txt(p207)
    attach_gse207422_sample_meta(a207, m207)
    out207 = args.outdir / "GSE207422.h5ad"
    a207.write_h5ad(out207, compression="gzip")
    print(f"wrote {out207} {a207.n_obs} x {a207.n_vars}", flush=True)
    del a207

    p241 = args.raw / "GSE241934_IIT_Matrix.mtx.gz"
    a241 = read_10x_mtx(
        p241,
        args.raw / "GSE241934_IIT_barcodes.tsv.gz",
        args.raw / "GSE241934_IIT_features.tsv.gz",
    )
    attach_gse241934_meta(a241, args.raw / "GSE241934_IIT_Meta.txt.gz", "GSE241934_IIT")
    out241 = args.outdir / "GSE241934_IIT.h5ad"
    a241.write_h5ad(out241, compression="gzip")
    print(f"wrote {out241} {a241.n_obs} x {a241.n_vars}", flush=True)
    del a241

    if args.include_rwc:
        a_r = read_10x_mtx(
            args.raw / "GSE241934_Real_Matrix.mtx.gz",
            args.raw / "GSE241934_RWC_barcodes.tsv.gz",
            args.raw / "GSE241934_RWC_features.tsv.gz",
        )
        attach_gse241934_meta(
            a_r, args.raw / "GSE241934_Real_Meta.txt.gz", "GSE241934_RWC"
        )
        out_r = args.outdir / "GSE241934_RWC.h5ad"
        a_r.write_h5ad(out_r, compression="gzip")
        print(f"wrote {out_r} {a_r.n_obs} x {a_r.n_vars}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
