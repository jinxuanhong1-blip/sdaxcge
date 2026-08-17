#!/usr/bin/env python3
"""Extract GSE189357 tumor epithelium (marker-malignant) to h5ad.

Gate (same as prior GSE189357 CLDN4 work): (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
Patient is the unit (TD1–TD9). No TACSTD2∩CLDN4 dual-high filter.
No normal-lung sample exists in this accession.
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

HERE = Path(__file__).resolve().parent
META = HERE.parent / "data" / "sample_metadata.tsv"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]


def _members(tf: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {m.name: m for m in tf.getmembers() if m.isfile()}


def _find(members: dict[str, tarfile.TarInfo], sample: str, tokens: tuple[str, ...]) -> str:
    hits = [
        n
        for n in members
        if f"_{sample}_" in n or n.endswith(f"_{sample}.mtx.gz")
        if any(t in n.lower() for t in tokens)
    ]
    # Prefer exact TD token surrounded by _ 
    hits = [n for n in hits if f"_{sample}_" in Path(n).name or f"_{sample}." in Path(n).name]
    if not hits:
        raise FileNotFoundError(f"{sample}: no file matching {tokens}")
    hits.sort(key=len)
    return hits[0]


def _read_tsv_gz(tf: tarfile.TarFile, member: tarfile.TarInfo) -> list[list[str]]:
    rows = []
    with gzip.GzipFile(fileobj=tf.extractfile(member)) as handle:
        for raw in handle:
            rows.append(raw.decode().rstrip("\n").split("\t"))
    return rows


def _read_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = _find(members, sample, ("features", "genes"))
    rows = _read_tsv_gz(tf, members[name])
    genes = []
    for p in rows:
        if len(p) >= 2 and not p[1].startswith("ENSG"):
            genes.append(p[1])
        else:
            genes.append(p[0])
    return genes


def _read_barcodes(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = _find(members, sample, ("barcodes",))
    return [r[0] for r in _read_tsv_gz(tf, members[name])]


def _mal_mask(mat: sparse.spmatrix, genes: list[str]) -> np.ndarray:
    upper = {g.upper(): i for i, g in enumerate(genes)}

    def col(g: str) -> np.ndarray:
        i = upper.get(g)
        if i is None:
            return np.zeros(mat.shape[1], dtype=float)
        return np.asarray(mat[i, :].todense()).ravel()

    epi = np.zeros(mat.shape[1], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    return epi & (col("PTPRC") == 0)


def extract(tar_path: Path) -> ad.AnnData:
    meta = pd.read_csv(META, sep="\t")
    meta = meta.set_index("patient")
    pieces: list[ad.AnnData] = []
    with tarfile.open(tar_path) as tf:
        members = _members(tf)
        print("tar members", len(members), flush=True)
        for sample in [f"TD{i}" for i in range(1, 10)]:
            print(f"reading {sample}", flush=True)
            genes = _read_features(tf, members, sample)
            barcodes = _read_barcodes(tf, members, sample)
            mtx_name = _find(members, sample, ("matrix.mtx",))
            with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as handle:
                mat = mmread(handle).tocsc()
            if mat.shape[0] != len(genes):
                if mat.shape[1] == len(genes):
                    mat = mat.T.tocsc()
                else:
                    raise RuntimeError(
                        f"{sample}: mtx {mat.shape} vs genes {len(genes)} barcodes {len(barcodes)}"
                    )
            if mat.shape[1] != len(barcodes):
                raise RuntimeError(f"{sample}: mtx cells {mat.shape[1]} vs barcodes {len(barcodes)}")
            mal = _mal_mask(mat, genes)
            n_mal = int(mal.sum())
            print(f"  {sample}: cells={mat.shape[1]} marker-malignant={n_mal}", flush=True)
            if n_mal == 0:
                continue
            X = mat[:, mal].T.tocsr()
            obs = pd.DataFrame(index=[f"{sample}_{b}" for b in np.asarray(barcodes)[mal]])
            obs["patient"] = sample
            obs["unit_id"] = sample
            obs["dataset"] = "GSE189357"
            obs["stage"] = str(meta.loc[sample, "stage"])
            obs["radiology"] = str(meta.loc[sample, "radiology"])
            obs["sex"] = str(meta.loc[sample, "sex"])
            obs["tissue"] = "TUMOR"
            obs["malig_def"] = "marker_malig"
            var = pd.DataFrame(index=pd.Index(genes, name="gene"))
            adata = ad.AnnData(X=X, obs=obs, var=var)
            adata.var_names_make_unique()
            pieces.append(adata)
            del mat, X
    if not pieces:
        raise SystemExit("no tumor epithelium extracted")
    # Harmonize genes
    genes = sorted(set().union(*[set(a.var_names) for a in pieces]))
    aligned = []
    for a in pieces:
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            pad = ad.AnnData(
                X=sparse.csr_matrix((a.n_obs, len(missing))),
                obs=a.obs.copy(),
                var=pd.DataFrame(index=missing),
            )
            a = ad.concat([a, pad], axis=1, merge="same")
        aligned.append(a[:, genes].copy())
    out = ad.concat(aligned, axis=0, merge="same")
    out.obs_names_make_unique()
    out.layers["counts"] = out.X.copy()
    out.uns["extraction"] = {
        "accession": "GSE189357",
        "gate": "(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0",
        "n_patients": int(out.obs["patient"].nunique()),
        "n_cells": int(out.n_obs),
        "dual_high": False,
    }
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tar", type=Path, default=Path("/tmp/gse189357/GSE189357_RAW.tar"))
    p.add_argument("--out", type=Path, default=Path("/tmp/gse189357/epithelium.h5ad"))
    args = p.parse_args()
    adata = extract(args.tar)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(
        f"wrote {args.out} cells={adata.n_obs} genes={adata.n_vars} "
        f"patients={adata.obs['patient'].nunique()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
