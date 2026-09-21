#!/usr/bin/env python3
"""Patient-batch scVI on the concordant-4 10% sample.

Same recipe as the locked scVI integration (n_latent=20, 2 layers, NB,
batch = dataset|unit, HVGs within dataset, seed 1), with CLDN4 and TACSTD2
forced into the model. The latent is the Milo graph. scANVI is not fit.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import numpy as np
import pandas as pd
import scanpy as sc
import scvi

HERE = Path(__file__).resolve().parent
PREP = HERE.parent / "scanvi_concordant4_cldn4" / "results" / "cache"
OUT = HERE / "results" / "cache"
DATASETS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]


def load() -> sc.AnnData:
    ads = []
    for ds in DATASETS:
        path = PREP / ds / "sub.h5ad"
        if not path.exists():
            raise SystemExit(f"missing {path}; run prepare.py first")
        ad = sc.read_h5ad(path)
        ad.obs_names_make_unique()
        ads.append(ad)
    adata = sc.concat(ads, join="inner", merge="same")
    adata.obs_names_make_unique()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.obs["batch"] = adata.obs["dataset"].astype(str) + "|" + adata.obs["unit_id"].astype(str)
    for col in ["dataset", "batch", "seed_class", "unit_id"]:
        adata.obs[col] = adata.obs[col].astype(str)
    return adata


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = 1
    adata = load()
    print(f"concat cells={adata.n_obs} genes={adata.n_vars} batches={adata.obs['batch'].nunique()}", flush=True)
    for gene in ("CLDN4", "TACSTD2"):
        if gene not in adata.var_names:
            raise SystemExit(f"{gene} dropped by the inner gene join")
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=2000,
        flavor="seurat_v3",
        batch_key="dataset",
        subset=False,
    )
    keep = adata.var["highly_variable"].to_numpy().copy()
    for gene in ("CLDN4", "TACSTD2"):
        keep[adata.var_names == gene] = True
    adata = adata[:, keep].copy()
    print(f"model genes={adata.n_vars} cells={adata.n_obs}", flush=True)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    model = scvi.model.SCVI(adata, n_latent=20, n_layers=2, gene_likelihood="nb")
    model.train(
        max_epochs=200,
        early_stopping=True,
        early_stopping_patience=15,
        batch_size=256,
        accelerator="cpu",
        devices=1,
        plan_kwargs={"lr": 1e-3},
    )
    latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
    epochs = int(len(model.history["elbo_train"]))
    obs = adata.obs[["dataset", "unit_id", "batch", "seed_class"]].copy()
    obs.insert(0, "cell_id", adata.obs_names.astype(str))
    np.save(OUT / "latent.npy", latent)
    obs.to_csv(OUT / "obs.tsv.gz", sep="\t", index=False)
    info = {
        "scvi_version": scvi.__version__,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_batches": int(adata.obs["batch"].nunique()),
        "n_latent": 20,
        "n_layers": 2,
        "gene_likelihood": "nb",
        "batch_key": "patient (dataset|unit_id)",
        "hvg_batch": "dataset",
        "forced_genes": ["CLDN4", "TACSTD2"],
        "seed": 1,
        "epochs_scvi": epochs,
        "early_stopping_patience": 15,
        "max_epochs": 200,
        "sample": "stable 10% of QC-pass cells; not the 220/140/40 cap",
    }
    (OUT / "train_info.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    main()
