#!/usr/bin/env python3
"""Build GSE123902 marker-epithelial AnnData from GEO dense UMI CSVs.

Locked gate (PR #459): (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
Tumor (PRIMARY/METASTASIS) is the analysis object. Matched NORMAL epithelium
is kept only so the root can be AT2-like and not CLDN4-high.
Author 36.5 GB H5 is not used. Marker epithelium is not CNV-malignant.
"""
from __future__ import annotations

import argparse
import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
from scipy import sparse

EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")
NAME_RE = re.compile(
    r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz"
)


def parse_name(fname: str) -> dict:
    m = NAME_RE.match(fname)
    if not m:
        raise ValueError(fname)
    site = m.group(3)
    tissue = {"PRIMARY_TUMOUR": "PRIMARY", "METASTASIS": "METASTASIS", "NORMAL": "NORMAL"}[site]
    return {"gsm": m.group(1), "donor": m.group(2), "site": site, "tissue": tissue, "file": fname}


def ensure_extracted(tar_path: Path, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    files = sorted(dest.glob("GSM*_dense.csv.gz"))
    if files:
        return files
    print(f"extract {tar_path} -> {dest}", flush=True)
    with tarfile.open(tar_path) as tf:
        tf.extractall(dest)
    files = sorted(dest.glob("GSM*_dense.csv.gz"))
    if not files:
        raise SystemExit(f"no GSM*_dense.csv.gz in {dest}")
    return files


def stream_epithelial(path: Path) -> tuple[list[str], list[str], list[np.ndarray], dict]:
    meta = parse_name(path.name)
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
        genes = [g.strip() for g in header[1:]]
        upper = [g.upper() for g in genes]
        idx = {g: i for i, g in enumerate(upper)}
        epi_i = [idx[g] for g in EPI if g in idx]
        ptprc_i = idx.get("PTPRC")
        barcodes: list[str] = []
        rows: list[np.ndarray] = []
        n_cells = 0
        n_epi = 0
        for line in handle:
            parts = line.rstrip("\n").split(",")
            n_cells += 1
            vals = np.fromiter((float(x) for x in parts[1:]), dtype=np.float32, count=len(genes))
            if not any(vals[i] > 0 for i in epi_i):
                continue
            if ptprc_i is not None and vals[ptprc_i] != 0:
                continue
            barcodes.append(f"{meta['gsm']}_{parts[0]}")
            rows.append(vals)
            n_epi += 1
    info = {
        **meta,
        "n_cells": n_cells,
        "n_epithelial": n_epi,
        "n_genes": len(genes),
    }
    print(
        f"  {path.name}: cells={n_cells} epi={n_epi} genes={len(genes)} {meta['donor']} {meta['tissue']}",
        flush=True,
    )
    return genes, barcodes, rows, info


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/gse123902_slingshot"))
    p.add_argument("--out", type=Path, default=Path("/tmp/gse123902_slingshot/epithelium.h5ad"))
    args = p.parse_args()

    import anndata as ad
    import pandas as pd

    tar = args.data / "GSE123902_RAW.tar"
    if not tar.exists():
        raise SystemExit(f"missing {tar}; run download.py")
    csv_dir = args.data / "csv"
    files = ensure_extracted(tar, csv_dir)

    all_genes = None
    X_blocks = []
    obs_rows = []
    audit = []
    for path in files:
        genes, barcodes, rows, info = stream_epithelial(path)
        audit.append(info)
        if not rows:
            continue
        if all_genes is None:
            all_genes = genes
        elif genes != all_genes:
            raise SystemExit(f"gene order mismatch in {path.name}")
        X_blocks.append(sparse.csr_matrix(np.vstack(rows)))
        for bc in barcodes:
            obs_rows.append(
                {
                    "barcode": bc,
                    "gsm": info["gsm"],
                    "donor": info["donor"],
                    "tissue": info["tissue"],
                    "site": info["site"],
                    "file": info["file"],
                    "dataset": "GSE123902",
                    "is_tumor": info["tissue"] != "NORMAL",
                }
            )

    if not X_blocks:
        raise SystemExit("no marker-epithelial cells")
    X = sparse.vstack(X_blocks, format="csr")
    obs = pd.DataFrame(obs_rows)
    obs.index = obs["barcode"].astype(str)
    var = pd.DataFrame(index=pd.Index(all_genes, name="gene"))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["n_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    audit_df = pd.DataFrame(audit)
    audit_path = args.out.with_name("extract_audit.tsv")
    audit_df.to_csv(audit_path, sep="\t", index=False)
    print(
        f"wrote {args.out} cells={adata.n_obs} genes={adata.n_vars} "
        f"donors={adata.obs['donor'].nunique()} tumor={int(adata.obs['is_tumor'].sum())}",
        flush=True,
    )
    print(audit_df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
