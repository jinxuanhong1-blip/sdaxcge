#!/usr/bin/env python3
"""GEO GSE148071 raw UMI → TISCH-labeled epithelium AnnData.

Keeps TISCH2 Malignant + leftover epithelium (Alveolar / Basal / Epithelial).
Does not download the 1.1 GB TISCH expression.h5. Counts come from the public
GSE148071_RAW.tar (42 sample UMI matrices).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import TISCH_EPITHELIAL  # noqa: E402

SAMPLE_RE = re.compile(r"^(GSM\d+)_P(\d+)_exp\.txt\.gz$")


def _load_tisch(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t")
    if "Cell" not in meta.columns:
        raise SystemExit(f"TISCH table missing Cell: {list(meta.columns)}")
    meta = meta.copy()
    meta["lineage"] = meta["Celltype (major-lineage)"].astype(str).str.strip()
    meta["malignancy"] = meta["Celltype (malignancy)"].astype(str).str.strip()
    meta["patient"] = meta["Patient"].astype(str).str.strip()
    meta["gsm"] = meta["Sample"].astype(str).str.strip()
    keep = meta["lineage"].isin(TISCH_EPITHELIAL)
    out = meta.loc[keep].set_index("Cell", drop=False)
    print(
        f"TISCH cells={len(meta)} epithelial={len(out)} "
        f"lineages={out['lineage'].value_counts().to_dict()}",
        flush=True,
    )
    return out


def _read_sample_subset(
    handle,
    gsm: str,
    keep_ids: set[str],
) -> tuple[list[str], list[str], sparse.csr_matrix]:
    with gzip.open(handle, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = header
        col_idx = [i for i, b in enumerate(barcodes) if f"{gsm}@{b}" in keep_ids]
        if not col_idx:
            return [], [], sparse.csr_matrix((0, 0), dtype=np.float32)
        cells = [f"{gsm}@{barcodes[i]}" for i in col_idx]
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz = 0
        for gi, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            vals = np.fromstring(raw[tab0 + 1 :], sep="\t", dtype=np.float32)
            if vals.size != len(barcodes):
                raise SystemExit(
                    f"{gsm} row {gi} {gene}: expected {len(barcodes)} values, got {vals.size}"
                )
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz += int(nz.size)
            indptr.append(nnz)
            genes.append(gene)
            if gi % 10000 == 0:
                print(f"  {gsm} streamed {gi} genes nnz={nnz}", flush=True)
    X = sparse.csr_matrix(
        (
            np.concatenate(data) if data else np.array([], dtype=np.float32),
            np.concatenate(indices) if indices else np.array([], dtype=np.int32),
            np.asarray(indptr, dtype=np.int32),
        ),
        shape=(len(genes), len(cells)),
        dtype=np.float32,
    )
    return genes, cells, X.T.tocsr()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/tmp/gse148071_paga_data")
    ap.add_argument("--out", default="/tmp/gse148071_paga_data/epithelium.h5ad")
    args = ap.parse_args()

    dest = Path(args.data)
    tar_path = dest / "GSE148071_RAW.tar"
    tisch_path = dest / "NSCLC_GSE148071_CellMetainfo_table.tsv"
    if not tar_path.is_file() or not tisch_path.is_file():
        raise SystemExit("missing RAW.tar or TISCH CellMetainfo; run download.sh")

    tisch = _load_tisch(tisch_path)
    keep_ids = set(tisch.index.astype(str))

    blocks: list[tuple[list[str], list[str], sparse.csr_matrix]] = []
    sample_rows = []
    with tarfile.open(tar_path, "r") as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        for mem in sorted(members, key=lambda m: m.name):
            m = SAMPLE_RE.match(Path(mem.name).name)
            if not m:
                print(f"skip unexpected member {mem.name}", flush=True)
                continue
            gsm, pnum = m.group(1), m.group(2)
            print(f">> {mem.name} patient=P{pnum}", flush=True)
            handle = tf.extractfile(mem)
            if handle is None:
                raise SystemExit(f"cannot extract {mem.name}")
            genes, cells, X = _read_sample_subset(handle, gsm, keep_ids)
            print(f"   kept cells={len(cells)} genes={len(genes)} nnz={X.nnz}", flush=True)
            sample_rows.append(
                {
                    "gsm": gsm,
                    "patient": f"P{pnum}",
                    "n_tisch_epi": int((tisch["gsm"] == gsm).sum()),
                    "n_matrix_matched": len(cells),
                }
            )
            if cells:
                blocks.append((genes, cells, X))

    if not blocks:
        raise SystemExit("no epithelial cells matched GEO barcodes to TISCH")

    gene_sets = [set(g) for g, _, _ in blocks]
    genes = sorted(set.intersection(*gene_sets))
    if not genes:
        raise SystemExit("empty gene intersection across samples")
    print(f"gene intersection n={len(genes)} (from {len(gene_sets)} samples)", flush=True)

    Xs = []
    cells_all: list[str] = []
    for g_i, cells, X in blocks:
        idx = {g: i for i, g in enumerate(g_i)}
        take = [idx[g] for g in genes]
        Xs.append(X[:, take])
        cells_all.extend(cells)
    Xall = sparse.vstack(Xs, format="csr")

    obs = tisch.reindex(cells_all)
    missing = int(obs["lineage"].isna().sum())
    if missing:
        print(f"WARNING dropped {missing} cells missing TISCH rows", flush=True)
        keep = obs["lineage"].notna().to_numpy()
        Xall = Xall[keep]
        obs = obs.loc[keep]
        cells_all = [c for c, k in zip(cells_all, keep) if k]
    obs = obs.copy()
    obs["cell_id"] = cells_all
    obs["n_umi"] = np.asarray(Xall.sum(axis=1)).ravel()
    obs["n_genes"] = np.asarray((Xall > 0).sum(axis=1)).ravel()

    adata = AnnData(X=Xall, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.obs_names = pd.Index(cells_all, name="cell")
    adata.var_names_make_unique()
    # keep uns JSON-safe (no list-of-dicts) so h5ad write succeeds
    lineage_counts = {str(k): int(v) for k, v in adata.obs["lineage"].value_counts().items()}
    extract = {
        "accession": "GSE148071",
        "paper": "Wu et al. Nat Commun 2021 PMID 33953163",
        "n_geo_samples": 42,
        "n_tisch_cells": 82267,
        "n_tisch_epithelial": int(len(tisch)),
        "n_matched": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "counts": "GEO raw UMI (GSE148071_RAW.tar)",
        "labels": "TISCH2 CellMetainfo major-lineage",
    }
    adata.uns["extract"] = extract
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary_path = out.parent / "extract_summary.json"
    summary_path.write_text(
        json.dumps({**extract, "lineage_counts": lineage_counts, "samples": sample_rows}, indent=2)
    )
    adata.write_h5ad(out)
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_patients": int(adata.obs["patient"].nunique()),
                "out": str(out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
