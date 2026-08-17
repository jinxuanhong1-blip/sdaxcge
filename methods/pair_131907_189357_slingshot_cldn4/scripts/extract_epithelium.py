#!/usr/bin/env python3
"""Build a joint epithelial AnnData for GSE131907 + GSE189357.

ADDITIVE. Epithelium only. GSE148071 is not added. No dual-high gate.
GSE131907: author Epithelial cells on tLung / nLung / tL/B / mLN / mBrain
(PE unlabeled epithelium is dropped). nLung AT2 is kept for the Slingshot root.
GSE189357: marker epithelium (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
"""

from __future__ import annotations

import argparse
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

KEEP_ORIGINS_131907 = ("tLung", "nLung", "tL/B", "mLN", "mBrain")
EPI_TYPE_131907 = "Epithelial cells"
EPI_MARKERS_189357 = ("EPCAM", "KRT8", "KRT18", "KRT19")
CAP_PER_UNIT = 350
SEED = 1

HERE = Path(__file__).resolve().parents[1]
STAGES_PATH = HERE / "data" / "gse189357_stages.tsv"


def cap_barcodes(
    frame: pd.DataFrame,
    unit_col: str,
    cap: int,
    seed: int,
    protect: pd.Series | None = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    for _, sub in frame.groupby(unit_col, observed=True):
        if protect is not None:
            prot = sub.index[protect.loc[sub.index].to_numpy()]
        else:
            prot = sub.index[:0]
        prot = pd.Index(prot)
        rest = sub.index.difference(prot)
        n_rest = max(0, cap - len(prot))
        if len(rest) > n_rest:
            chosen = rng.choice(rest.to_numpy(), size=n_rest, replace=False)
            picked = prot.append(pd.Index(chosen))
        else:
            picked = prot.append(rest)
        keep_idx.extend(picked.tolist())
    return frame.loc[keep_idx].copy()


def stream_umi_subset(
    umi_path: Path,
    keep_ids: list[str],
    gene_cap: int | None = None,
) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE131907 cell IDs in UMI header")
        ordered_cells = [cell_ids[i] for i in col_idx]
        print(
            f"GSE131907 UMI header cells={len(cell_ids)} keep={len(col_idx)}",
            flush=True,
        )
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
            rest = raw[tab0 + 1 :]
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(
                    f"row {gi} {gene}: expected {len(cell_ids)} values, got {vals.size}"
                )
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 2000 == 0:
                print(f"  streamed {gi} genes, nnz={nnz_total}", flush=True)
            if gene_cap is not None and gi >= gene_cap:
                break
    if data:
        data_a = np.concatenate(data)
        indices_a = np.concatenate(indices)
    else:
        data_a = np.array([], dtype=np.float32)
        indices_a = np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, indices_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered_cells)),
        dtype=np.float32,
    )
    return ordered_cells, genes, mat


def extract_gse131907(data: Path, gene_cap: int | None, cap: int) -> "ad.AnnData":
    import anndata as ad

    ann_path = data / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
    epi = ann.loc[
        (ann["Cell_type"] == EPI_TYPE_131907)
        & (ann["Sample_Origin"].isin(KEEP_ORIGINS_131907))
    ].copy()
    epi["Index"] = epi["Index"].astype(str)
    epi = epi.set_index("Index", drop=False)
    protect = epi["Cell_subtype"].astype(str).eq("AT2") & epi["Sample_Origin"].eq("nLung")
    catalog = {
        "n_ann_rows": int(len(ann)),
        "n_epithelial_kept_origins": int(len(epi)),
        "by_origin": epi["Sample_Origin"].value_counts().to_dict(),
        "by_subtype": epi["Cell_subtype"].fillna("NA").value_counts().to_dict(),
        "n_samples": int(epi["Sample"].nunique()),
        "n_nLung_AT2": int(protect.sum()),
    }
    print(json.dumps({"GSE131907_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi, "Sample", cap, SEED, protect)
    keep_ids = epi["Index"].astype(str).tolist()
    cells, genes, mat = stream_umi_subset(umi_path, keep_ids, gene_cap=gene_cap)
    X = mat.T.tocsr()
    obs = epi.set_index("Index").reindex(cells)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.obs["dataset"] = "GSE131907"
    adata.obs["histology"] = "LUAD"
    adata.obs["unit_id"] = "GSE131907:" + adata.obs["Sample"].astype(str)
    adata.obs["patient_id"] = (
        adata.obs["Sample"]
        .astype(str)
        .str.extract(r"(?:LUNG_[NT]|EBUS_|BRONCHO_|EFFUSION_)?(\d+)", expand=False)
    )
    adata.obs["author_lineage"] = "Epithelial cells"
    adata.obs["author_subtype"] = adata.obs["Cell_subtype"].astype(str)
    adata.obs["tissue"] = adata.obs["Sample_Origin"].astype(str)
    adata.obs["is_normal_tissue"] = adata.obs["Sample_Origin"].eq("nLung")
    adata.obs["stage"] = np.where(adata.obs["Sample_Origin"].eq("nLung"), "nLung", "tumor")
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE131907 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def _read_10x_names(tf: tarfile.TarFile, members: dict, sample: str, kind: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_{kind}" in n)
    out = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            parts = line.decode().strip().split("\t")
            if kind.startswith("feature") or kind.startswith("gene"):
                out.append((parts[1] if len(parts) > 1 else parts[0]))
            else:
                out.append(parts[0])
    return out


def _extract_mtx(tf: tarfile.TarFile, members: dict, sample: str, dest: Path) -> Path:
    name = next(n for n in members if f"_{sample}_matrix.mtx" in n)
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"{sample}_matrix.mtx"
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as src, out.open("wb") as fh:
        while True:
            chunk = src.read(8 * 1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    return out


def _collapse_duplicate_genes(mat: sparse.spmatrix, genes: list[str]) -> tuple[sparse.csr_matrix, list[str]]:
    """Sum rows that share a gene symbol (10x ENSG collisions after symbol map)."""
    series = pd.Series(np.arange(len(genes)), index=pd.Index(genes, name="gene"))
    if not series.index.has_duplicates:
        return mat.tocsr(), genes
    groups = series.groupby(level=0).apply(lambda s: s.to_numpy())
    uniq = list(groups.index)
    rows = []
    for idxs in groups:
        if len(idxs) == 1:
            rows.append(mat[idxs[0]])
        else:
            rows.append(mat[idxs].sum(axis=0))
    stacked = sparse.vstack(rows, format="csr")
    return stacked, uniq


def extract_gse189357(data: Path, cap: int, scratch: Path) -> "ad.AnnData":
    import anndata as ad

    tar_path = data / "gse189357" / "GSE189357_RAW.tar"
    stages = pd.read_csv(STAGES_PATH, sep="\t").set_index("patient")
    pieces = []
    catalog_rows = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            print(f"GSE189357 reading {sample}", flush=True)
            genes = _read_10x_names(tf, members, sample, "features")
            barcodes = _read_10x_names(tf, members, sample, "barcodes")
            mtx_path = _extract_mtx(tf, members, sample, scratch)
            mat = mmread(mtx_path).tocsr()  # genes x cells
            if mat.shape[0] != len(genes) or mat.shape[1] != len(barcodes):
                raise SystemExit(
                    f"{sample}: mtx {mat.shape} vs genes {len(genes)} barcodes {len(barcodes)}"
                )
            mat, genes = _collapse_duplicate_genes(mat, genes)
            upper = {g.upper(): i for i, g in enumerate(genes)}

            def col(name: str) -> np.ndarray:
                i = upper.get(name.upper())
                if i is None:
                    return np.zeros(mat.shape[1], dtype=np.float32)
                return np.asarray(mat[i].todense()).ravel()

            epi_pos = np.zeros(mat.shape[1], dtype=bool)
            for g in EPI_MARKERS_189357:
                epi_pos |= col(g) > 0
            keep = epi_pos & (col("PTPRC") == 0)
            n_keep = int(keep.sum())
            catalog_rows.append(
                {
                    "patient": sample,
                    "n_cells": int(mat.shape[1]),
                    "n_marker_epithelium": n_keep,
                    "histology": str(stages.loc[sample, "histology"]),
                }
            )
            print(
                f"  {sample}: cells={mat.shape[1]} marker_epi={n_keep} "
                f"stage={stages.loc[sample, 'histology']}",
                flush=True,
            )
            if n_keep == 0:
                del mat
                continue
            sub = mat[:, keep].T.tocsr()
            obs = pd.DataFrame(
                {
                    "Index": [f"{sample}:{b}" for b, k in zip(barcodes, keep) if k],
                    "Sample": sample,
                    "patient": sample,
                    "barcode": [b for b, k in zip(barcodes, keep) if k],
                }
            ).set_index("Index")
            obs["dataset"] = "GSE189357"
            obs["histology"] = str(stages.loc[sample, "histology"])
            obs["Sample_Origin"] = "tumor"
            obs["Cell_type"] = "marker_epithelium"
            obs["Cell_subtype"] = "marker_epithelium"
            obs["unit_id"] = "GSE189357:" + sample
            obs["patient_id"] = sample
            obs["author_lineage"] = "marker_epithelium"
            obs["author_subtype"] = "marker_epithelium"
            obs["tissue"] = "tumor"
            obs["is_normal_tissue"] = False
            obs["stage"] = str(stages.loc[sample, "histology"])
            adata_s = ad.AnnData(
                X=sub,
                obs=obs,
                var=pd.DataFrame(index=pd.Index(genes, name="gene")),
            )
            adata_s = adata_s[
                cap_barcodes(adata_s.obs, "Sample", cap, SEED + 10 + int(sample[2:])).index
            ].copy()
            pieces.append(adata_s)
            del mat
    if not pieces:
        raise SystemExit("GSE189357: no marker-epithelial cells")
    print(json.dumps({"GSE189357_catalog": catalog_rows}, indent=2), flush=True)
    shared = pieces[0].var_names
    for p in pieces[1:]:
        shared = shared.intersection(p.var_names)
    pieces = [p[:, shared].copy() for p in pieces]
    out = ad.concat(pieces, axis=0, join="inner", merge="same")
    out.layers["counts"] = out.X.copy()
    print(f"GSE189357 written cells={out.n_obs} genes={out.n_vars}", flush=True)
    return out


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif pd.api.types.is_categorical_dtype(out[col]) or out[col].dtype == object:
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def concat_shared(a, b):
    import anndata as ad

    shared = a.var_names.intersection(b.var_names)
    if len(shared) < 5000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    a2 = a[:, shared].copy()
    b2 = b[:, shared].copy()
    a2.obs_names = "GSE131907:" + a2.obs_names.astype(str)
    b2.obs_names = "GSE189357:" + b2.obs_names.astype(str)
    cols = sorted(set(a2.obs.columns) | set(b2.obs.columns))
    for frame in (a2, b2):
        for c in cols:
            if c not in frame.obs.columns:
                frame.obs[c] = "NA"
        frame.obs = _sanitize_obs(frame.obs[cols])
    out = ad.concat([a2, b2], axis=0, join="inner", merge="same")
    out.obs = _sanitize_obs(out.obs)
    out.layers["counts"] = out.X.copy()
    return out, int(len(shared))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_131907_189357"))
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_131907_189357/epithelium.h5ad"))
    p.add_argument("--scratch", type=Path, default=Path("/tmp/pair_131907_189357/extract"))
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    p.add_argument("--gene-cap", type=int, default=None, help="debug only")
    args = p.parse_args()

    cache_a = args.data / "GSE131907_epithelium_cap.h5ad"
    cache_b = args.data / "GSE189357_epithelium_cap.h5ad"
    import anndata as ad

    if cache_a.is_file() and args.gene_cap is None:
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
    else:
        a = extract_gse131907(args.data, args.gene_cap, args.cap_per_unit)
        a.obs = _sanitize_obs(a.obs)
        a.write_h5ad(cache_a)
        print(f"cached {cache_a}", flush=True)
    if cache_b.is_file() and args.gene_cap is None:
        print(f"reuse {cache_b}", flush=True)
        b = ad.read_h5ad(cache_b)
    else:
        b = extract_gse189357(args.data, args.cap_per_unit, args.scratch)
        b.obs = _sanitize_obs(b.obs)
        b.write_h5ad(cache_b)
        print(f"cached {cache_b}", flush=True)
    joint, n_shared = concat_shared(a, b)
    del a, b
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes_shared": n_shared,
        "by_dataset": joint.obs["dataset"].value_counts().to_dict(),
        "by_origin": joint.obs["Sample_Origin"].astype(str).value_counts().to_dict(),
        "by_subtype": joint.obs["author_subtype"].astype(str).value_counts().to_dict(),
        "by_stage": joint.obs["stage"].astype(str).value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "dual_high": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
