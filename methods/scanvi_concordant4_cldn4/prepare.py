#!/usr/bin/env python3
"""Build full-unit CLDN4/T/NK tables and a stratified subsample h5ad per GSE.

The patient / donor / sample table uses every cell in the locked unit (no cap).
The h5ad is a QC'd stratified subsample (up to 220 malignant + 140 T/NK + 40
other cells per unit) for scVI/scANVI. Composition is never taken from that cap.
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from common import (
    AUTHOR_DATASETS,
    EPI_GENES,
    GSE131907_UNITS,
    MAL_CAP,
    OTHER_CAP,
    QC_MAX_MT,
    QC_MIN_GENES,
    QC_MIN_UMI,
    TNK_CAP,
    TNK_GENES,
)

HERE = Path(__file__).resolve().parent
CACHE = HERE / "results" / "cache"


def rng_for(unit_id: str) -> np.random.Generator:
    digest = hashlib.md5(unit_id.encode()).hexdigest()[:8]
    return np.random.default_rng(int(digest, 16))


def collapse_genes(mat: sparse.spmatrix, genes: np.ndarray):
    mat = sparse.csr_matrix(mat)
    genes = np.array([str(g).upper() for g in genes])
    order = np.argsort(genes, kind="mergesort")
    genes = genes[order]
    mat = mat[order]
    if len(genes) <= 1 or np.all(genes[1:] != genes[:-1]):
        return mat, genes
    change = np.nonzero(genes[1:] != genes[:-1])[0] + 1
    bounds = np.concatenate([[0], change, [len(genes)]])
    data, indices, indptr, names = [], [], [0], []
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        names.append(genes[b0])
        block = sparse.csr_matrix(mat[b0:b1].sum(axis=0))
        data.append(block.data.astype(np.float32, copy=False))
        indices.append(block.indices)
        indptr.append(indptr[-1] + block.nnz)
    out = sparse.csr_matrix(
        (np.concatenate(data) if data else np.array([], dtype=np.float32),
         np.concatenate(indices) if indices else np.array([], dtype=np.int32),
         indptr),
        shape=(len(names), mat.shape[1]),
    )
    return out, np.array(names)


def as_count_data(mat: sparse.csr_matrix) -> sparse.csr_matrix:
    if mat.nnz == 0:
        return mat
    rounded = np.rint(mat.data)
    if np.max(np.abs(mat.data - rounded)) < 1e-3:
        mat = mat.copy()
        mat.data = rounded.astype(np.float32)
    return mat


def gene_index(genes: np.ndarray) -> dict[str, int]:
    return {g: i for i, g in enumerate(genes)}


def gene_vec(mat: sparse.csr_matrix, gindex: dict[str, int], gene: str) -> np.ndarray:
    i = gindex.get(gene)
    if i is None:
        return np.zeros(mat.shape[1], dtype=np.float64)
    return np.asarray(mat[i].todense(), dtype=np.float64).ravel()


def marker_masks(mat: sparse.csr_matrix, gindex: dict[str, int]):
    epi = np.zeros(mat.shape[1], dtype=bool)
    for g in EPI_GENES:
        epi |= gene_vec(mat, gindex, g) > 0
    ptprc = gene_vec(mat, gindex, "PTPRC")
    mal = epi & (ptprc == 0)
    tnk = np.zeros(mat.shape[1], dtype=bool)
    for g in TNK_GENES:
        tnk |= gene_vec(mat, gindex, g) > 0
    tnk &= ~mal
    return mal, tnk


def qc_mask(mat: sparse.csr_matrix, genes: np.ndarray) -> np.ndarray:
    ncount = np.asarray(mat.sum(axis=0)).ravel()
    nfeat = np.asarray((mat > 0).sum(axis=0)).ravel()
    mt = np.array([g.startswith("MT-") for g in genes])
    if mt.any():
        mt_pct = np.asarray(mat[mt].sum(axis=0)).ravel() / np.maximum(ncount, 1) * 100
    else:
        mt_pct = np.zeros(mat.shape[1])
    return (nfeat >= QC_MIN_GENES) & (ncount >= QC_MIN_UMI) & (mt_pct < QC_MAX_MT)


def stratified_columns(mal: np.ndarray, tnk: np.ndarray, qc: np.ndarray, unit_id: str) -> np.ndarray:
    rng = rng_for(unit_id)
    ok = np.flatnonzero(qc)
    mal_i = ok[mal[ok]]
    tnk_i = ok[tnk[ok] & ~mal[ok]]
    oth_i = ok[~mal[ok] & ~tnk[ok]]

    def take(ix, k):
        if len(ix) <= k:
            return ix
        return rng.choice(ix, size=k, replace=False)

    chosen = np.concatenate([take(mal_i, MAL_CAP), take(tnk_i, TNK_CAP), take(oth_i, OTHER_CAP)])
    return np.unique(chosen.astype(int))


def cldn4_stats(umi: np.ndarray) -> tuple[float, float]:
    if len(umi) == 0:
        return float("nan"), float("nan")
    return float(100.0 * np.mean(umi > 0)), float(np.mean(np.log1p(umi)))


def unit_row(dataset, unit_id, unit_type, tissue, n_cells, n_mal, n_tnk, umi, malig_def) -> dict:
    pct, mean = cldn4_stats(umi)
    return {
        "dataset": dataset,
        "unit_id": unit_id,
        "patient_id": unit_id,
        "unit_type": unit_type,
        "tissue": tissue,
        "n_cells": int(n_cells),
        "n_malignant": int(n_mal),
        "n_tnk": int(n_tnk),
        "frac_tnk": float(n_tnk / n_cells) if n_cells else float("nan"),
        "mal_CLDN4_pct": pct,
        "mal_CLDN4_mean_log1p": mean,
        "malig_def": malig_def,
    }


def hstack_union(pieces: list[tuple[sparse.csr_matrix, np.ndarray]]) -> tuple[sparse.csr_matrix, np.ndarray]:
    genes0 = pieces[0][1]
    if all(len(g) == len(genes0) and np.array_equal(g, genes0) for _, g in pieces):
        return sparse.hstack([m for m, _ in pieces]).tocsr(), genes0
    print("gene panels differ; aligning to the union", flush=True)
    union = pd.Index(genes0.astype(str))
    for _, g in pieces[1:]:
        union = union.append(pd.Index(g.astype(str))).unique()
    genes = union.to_numpy()
    blocks = []
    for mat, g in pieces:
        indexer = pd.Index(g.astype(str)).get_indexer(genes)
        present = np.flatnonzero(indexer >= 0)
        mapped = sparse.lil_matrix((len(genes), mat.shape[1]), dtype=np.float32)
        mapped[present, :] = mat[indexer[present]]
        blocks.append(mapped.tocsr())
    return sparse.hstack(blocks).tocsr(), genes


def write_h5ad(path: Path, mat_gxc: sparse.csr_matrix, genes: np.ndarray, meta: pd.DataFrame) -> None:
    import anndata as ad

    path.parent.mkdir(parents=True, exist_ok=True)
    # Drop all-zero genes to keep the subsample small; CLDN4 is retained if present.
    keep = np.asarray(mat_gxc.sum(axis=1)).ravel() > 0
    if "CLDN4" in set(genes) and not keep[np.where(genes == "CLDN4")[0][0]]:
        keep[np.where(genes == "CLDN4")[0][0]] = True
    mat_gxc = mat_gxc[keep]
    genes = genes[keep]
    X = mat_gxc.T.tocsr().astype(np.float32)
    obs = meta.copy()
    obs.index = obs["cell_id"].astype(str)
    obs["scanvi_label"] = np.where(obs["dataset"].isin(list(AUTHOR_DATASETS)), obs["seed_class"], "Unknown")
    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = pd.Index(genes.astype(str))
    adata.var_names_make_unique()
    adata.write_h5ad(path)
    print(f"wrote {path} cells={adata.n_obs} genes={adata.n_vars}", flush=True)


def prepare_gse123902(raw: Path) -> None:
    out = CACHE / "GSE123902"
    out.mkdir(parents=True, exist_ok=True)
    tar_path = raw / "GSE123902_RAW.tar"
    if not (out / "_extracted").exists():
        print("extract GSE123902 tar", flush=True)
        with tarfile.open(tar_path) as tar:
            tar.extractall(out / "_extracted")
    locked = pd.read_csv(HERE / "data" / "GSE123902_marker_units.tsv", sep="\t")
    locked = locked[locked["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    locked = locked.sort_values(["patient", "tissue"]).drop_duplicates("patient")
    # Keep the 13 donors that entered the locked n=65 table.
    ref = pd.read_csv(HERE / "data" / "locked_patient_units.tsv", sep="\t")
    keep_ids = set(ref.loc[ref["dataset"] == "GSE123902", "unit_id"])
    locked = locked[locked["patient"].isin(keep_ids)]
    pieces = []
    metas = []
    units = []
    for rec in locked.itertuples(index=False):
        path = out / "_extracted" / rec.file
        print(f"GSE123902 {rec.patient} {path.name}", flush=True)
        df = pd.read_csv(path, index_col=0)
        mat, genes = collapse_genes(sparse.csr_matrix(df.to_numpy(dtype=np.float32).T), df.columns.to_numpy())
        mat = as_count_data(mat)
        gindex = gene_index(genes)
        if "CLDN4" not in gindex:
            raise SystemExit(f"CLDN4 missing in {path.name}")
        mal, tnk = marker_masks(mat, gindex)
        cld = gene_vec(mat, gindex, "CLDN4")
        units.append(
            unit_row("GSE123902", rec.patient, "donor", rec.tissue, mat.shape[1], int(mal.sum()), int(tnk.sum()), cld[mal], "marker_malig")
        )
        qc = qc_mask(mat, genes)
        cols = stratified_columns(mal, tnk, qc, f"GSE123902|{rec.patient}")
        sub = mat[:, cols]
        barcodes = df.index.to_numpy()[cols]
        meta = pd.DataFrame({
            "cell_id": [f"GSE123902|{rec.patient}|{b}" for b in barcodes],
            "dataset": "GSE123902",
            "unit_id": rec.patient,
            "patient_id": rec.patient,
            "unit_type": "donor",
            "barcode": barcodes.astype(str),
            "seed_class": np.where(mal[cols], "malignant", np.where(tnk[cols], "T/NK", "other")),
        })
        pieces.append((sub, genes))
        metas.append(meta)
        print(f"  cells={mat.shape[1]} mal={int(mal.sum())} tnk={int(tnk.sum())} sub={len(cols)} pct={units[-1]['mal_CLDN4_pct']:.2f}", flush=True)
        del df
    big, all_genes = hstack_union(pieces)
    meta = pd.concat(metas, ignore_index=True)
    units_df = pd.DataFrame(units)
    units_df.to_csv(out / "units.tsv", sep="\t", index=False)
    write_h5ad(out / "sub.h5ad", big, all_genes, meta)


def prepare_gse189357(raw: Path) -> None:
    out = CACHE / "GSE189357"
    out.mkdir(parents=True, exist_ok=True)
    tar_path = raw / "GSE189357_RAW.tar"
    if not (out / "_extracted").exists():
        print("extract GSE189357 tar", flush=True)
        with tarfile.open(tar_path) as tar:
            tar.extractall(out / "_extracted")
    root = out / "_extracted"
    ref = pd.read_csv(HERE / "data" / "locked_patient_units.tsv", sep="\t")
    patients = list(ref.loc[ref["dataset"] == "GSE189357", "unit_id"])
    pieces, metas, units = [], [], []
    for pat in patients:
        mtx = next(root.glob(f"*_{pat}_matrix.mtx.gz"))
        feat = next(root.glob(f"*_{pat}_features.tsv.gz"))
        bc = next(root.glob(f"*_{pat}_barcodes.tsv.gz"))
        print(f"GSE189357 {pat}", flush=True)
        features = pd.read_csv(feat, sep="\t", header=None)
        genes = features.iloc[:, 1].astype(str).to_numpy() if features.shape[1] > 1 else features.iloc[:, 0].astype(str).to_numpy()
        barcodes = pd.read_csv(bc, sep="\t", header=None).iloc[:, 0].astype(str).to_numpy()
        import gzip
        with gzip.open(mtx, "rb") as fh:
            mat = mmread(fh).tocsr()
        mat, genes = collapse_genes(mat, genes)
        mat = as_count_data(mat)
        if mat.shape[1] != len(barcodes):
            raise SystemExit(f"{pat} matrix cells {mat.shape[1]} != barcodes {len(barcodes)}")
        gindex = gene_index(genes)
        if "CLDN4" not in gindex:
            raise SystemExit(f"CLDN4 missing in {pat}")
        mal, tnk = marker_masks(mat, gindex)
        cld = gene_vec(mat, gindex, "CLDN4")
        units.append(unit_row("GSE189357", pat, "patient", "TUMOR", mat.shape[1], int(mal.sum()), int(tnk.sum()), cld[mal], "marker_malig"))
        qc = qc_mask(mat, genes)
        cols = stratified_columns(mal, tnk, qc, f"GSE189357|{pat}")
        sub = mat[:, cols]
        pieces.append((sub, genes))
        metas.append(pd.DataFrame({
            "cell_id": [f"GSE189357|{pat}|{b}" for b in barcodes[cols]],
            "dataset": "GSE189357",
            "unit_id": pat,
            "patient_id": pat,
            "unit_type": "patient",
            "barcode": barcodes[cols],
            "seed_class": np.where(mal[cols], "malignant", np.where(tnk[cols], "T/NK", "other")),
        }))
        print(f"  cells={mat.shape[1]} mal={int(mal.sum())} tnk={int(tnk.sum())} sub={len(cols)} pct={units[-1]['mal_CLDN4_pct']:.2f}", flush=True)
        del mat
    big, all_genes = hstack_union(pieces)
    pd.DataFrame(units).to_csv(out / "units.tsv", sep="\t", index=False)
    write_h5ad(out / "sub.h5ad", big, all_genes, pd.concat(metas, ignore_index=True))


def prepare_gse131907(raw: Path) -> None:
    out = CACHE / "GSE131907"
    out.mkdir(parents=True, exist_ok=True)
    ann_path = raw / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = raw / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    ann = ann[ann["Sample"].isin(GSE131907_UNITS)].copy()
    ann["author_malignant"] = ann["Cell_subtype"].eq("Malignant cells")
    ann["author_tnk"] = ann["Cell_type"].isin(["T lymphocytes", "NK cells"])
    print(f"GSE131907 annotation rows in 21 samples: {len(ann)}", flush=True)
    chosen_rows = []
    for sample, g in ann.groupby("Sample", sort=True):
        mal = g.index[g["author_malignant"].to_numpy()]
        tnk = g.index[g["author_tnk"].to_numpy() & ~g["author_malignant"].to_numpy()]
        oth = g.index[~g["author_malignant"].to_numpy() & ~g["author_tnk"].to_numpy()]
        rng = rng_for(f"GSE131907|{sample}")

        def take(ix, k):
            ix = np.asarray(ix)
            if len(ix) <= k:
                return ix
            return rng.choice(ix, size=k, replace=False)

        pick = np.concatenate([take(mal, MAL_CAP), take(tnk, TNK_CAP), take(oth, OTHER_CAP)])
        chosen_rows.append(g.loc[pick])
    sub_ann = pd.concat(chosen_rows, axis=0)
    sub_ids = set(sub_ann["Index"])
    mal_ids = set(ann.loc[ann["author_malignant"], "Index"])
    flags = out / "barcodes.tsv"
    with flags.open("w") as fh:
        seen = set()
        for idx in list(sub_ann["Index"]) + list(ann.loc[ann["author_malignant"], "Index"]):
            if idx in seen:
                continue
            seen.add(idx)
            fh.write(f"{idx} {1 if idx in sub_ids else 0} {1 if idx in mal_ids else 0}\n")
    exe = HERE / "src" / "extract_cols"
    if not exe.exists():
        subprocess.run(["gcc", "-O3", "-o", str(exe), str(HERE / "src" / "extract_cols.c"), "-lz"], check=True)
    prefix = out / "stream"
    print("stream GSE131907 matrix", flush=True)
    subprocess.run([str(exe), str(mat_path), str(flags), str(prefix)], check=True)
    genes = pd.read_csv(prefix.with_suffix(".genes.txt"), header=None)[0].astype(str).to_numpy()
    # suffix trick: prefix.genes.txt is not with_suffix of a path that already has dots.
    genes = pd.read_csv(str(prefix) + ".genes.txt", header=None)[0].astype(str).to_numpy()
    cells = pd.read_csv(str(prefix) + ".cells.txt", header=None)[0].astype(str).to_numpy()
    meta_txt = (str(prefix) + ".meta.txt")
    meta_kv = dict(line.rstrip("\n").split("\t") for line in open(meta_txt) if "\t" in line)
    n_genes, n_cells, nnz = int(meta_kv["n_genes"]), int(meta_kv["n_cells"]), int(meta_kv["nnz"])
    if int(meta_kv["cldn4_found"]) != 1:
        raise SystemExit("CLDN4 row was not found in GSE131907")
    rec = np.fromfile(str(prefix) + ".coo.bin", dtype=[("g", "<i4"), ("c", "<i4"), ("v", "<f4")])
    if len(rec) != nnz:
        print(f"warning nnz meta {nnz} file {len(rec)}", flush=True)
    mat = sparse.coo_matrix((rec["v"], (rec["g"], rec["c"])), shape=(n_genes, n_cells)).tocsr()
    del rec
    mat, genes = collapse_genes(mat, genes)
    mat = as_count_data(mat)
    ann_i = ann.drop_duplicates("Index").set_index("Index")
    sub_ann = sub_ann.drop_duplicates("Index")
    # cells.txt order is the subsample output order
    info = ann_i.loc[cells]
    mal = info["author_malignant"].to_numpy()
    tnk = info["author_tnk"].to_numpy() & ~mal
    qc = qc_mask(mat, genes)
    keep = qc
    print(f"GSE131907 subsample {n_cells} QC pass {int(keep.sum())}", flush=True)
    mat = mat[:, keep]
    cells_k = cells[keep]
    info = info.loc[cells_k]
    mal = mal[keep]
    tnk = tnk[keep]
    meta = pd.DataFrame({
        "cell_id": [f"GSE131907|{s}|{b}" for s, b in zip(info["Sample"], cells_k)],
        "dataset": "GSE131907",
        "unit_id": info["Sample"].to_numpy(),
        "patient_id": info["Sample"].to_numpy(),
        "unit_type": "sample",
        "barcode": cells_k,
        "seed_class": np.where(mal, "malignant", np.where(tnk, "T/NK", "other")),
    })
    cld = pd.read_csv(str(prefix) + ".cldn4.tsv", sep="\t", header=None, names=["Index", "umi"])
    cld = cld.drop_duplicates("Index")
    ann_m = ann.loc[ann["author_malignant"], ["Index", "Sample"]].merge(cld, on="Index", how="left")
    missing = int(ann_m["umi"].isna().sum())
    if missing:
        print(f"warning malignant cells missing from CLDN4 row: {missing}", flush=True)
    ann_m["umi"] = ann_m["umi"].fillna(0)
    units = []
    for sample, g in ann.groupby("Sample", sort=True):
        n_cells_s = len(g)
        n_mal = int(g["author_malignant"].sum())
        n_tnk = int(g["author_tnk"].sum())
        umi = ann_m.loc[ann_m["Sample"] == sample, "umi"].to_numpy(dtype=float)
        origin = g["Sample_Origin"].iloc[0]
        units.append(unit_row("GSE131907", sample, "sample", origin, n_cells_s, n_mal, n_tnk, umi, "author_malig"))
        print(f"  {sample} cells={n_cells_s} mal={n_mal} tnk={n_tnk} pct={units[-1]['mal_CLDN4_pct']:.2f}", flush=True)
    pd.DataFrame(units).to_csv(out / "units.tsv", sep="\t", index=False)
    write_h5ad(out / "sub.h5ad", mat, genes, meta)


def prepare_gse205335(raw: Path) -> None:
    out = CACHE / "GSE205335"
    if not (out / "units.tsv").exists() or not (out / "matrix.mtx").exists():
        print("R extract GSE205335", flush=True)
        subprocess.run(["Rscript", str(HERE / "extract_gse205335.R"), str(raw), str(out)], check=True)
    print("read GSE205335 subsample", flush=True)
    mat = mmread(out / "matrix.mtx").tocsr()
    genes = pd.read_csv(out / "features.tsv", header=None)[0].astype(str).to_numpy()
    barcodes = pd.read_csv(out / "barcodes.tsv", header=None)[0].astype(str).to_numpy()
    meta = pd.read_csv(out / "meta.tsv", sep="\t")
    if mat.shape[1] != len(barcodes):
        raise SystemExit(f"GSE205335 mtx {mat.shape[1]} vs barcodes {len(barcodes)}")
    meta = meta.set_index("cell_id").loc[barcodes].reset_index()
    mat, genes = collapse_genes(mat, genes)
    mat = as_count_data(mat)
    write_h5ad(out / "sub.h5ad", mat, genes, meta)


def compare_locked() -> None:
    parts = []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        p = CACHE / ds / "units.tsv"
        if p.exists():
            parts.append(pd.read_csv(p, sep="\t"))
    if not parts:
        return
    got = pd.concat(parts, ignore_index=True)
    locked = pd.read_csv(HERE / "data" / "locked_patient_units.tsv", sep="\t")
    m = locked.merge(got, on=["dataset", "unit_id"], suffixes=("_locked", "_new"))
    m["d_frac"] = m["frac_tnk_new"] - m["frac_tnk_locked"]
    m["d_pct"] = m["mal_CLDN4_pct_new"] - m["mal_CLDN4_pct_locked"]
    m["d_n"] = m["n_cells_new"] - m["n_cells_locked"]
    tab = HERE / "results" / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    cols = ["dataset", "unit_id", "n_cells_locked", "n_cells_new", "d_n", "frac_tnk_locked", "frac_tnk_new", "d_frac",
            "mal_CLDN4_pct_locked", "mal_CLDN4_pct_new", "d_pct", "n_malignant_locked", "n_malignant_new", "n_tnk_locked", "n_tnk_new"]
    m[cols].to_csv(tab / "locked_replication_check.tsv", sep="\t", index=False)
    print("locked check max |d_frac|", float(m["d_frac"].abs().max()), "max |d_pct|", float(m["d_pct"].abs().max()), "max |d_n|", float(m["d_n"].abs().max()), "n", len(m), flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--raw", type=Path, default=Path("/tmp/geo_c4"))
    p.add_argument("--dataset", default="all", choices=["all", "GSE123902", "GSE131907", "GSE205335", "GSE189357"])
    args = p.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    todo = ["GSE123902", "GSE189357", "GSE131907", "GSE205335"] if args.dataset == "all" else [args.dataset]
    fns = {
        "GSE123902": prepare_gse123902,
        "GSE189357": prepare_gse189357,
        "GSE131907": prepare_gse131907,
        "GSE205335": prepare_gse205335,
    }
    for ds in todo:
        fns[ds](args.raw)
    compare_locked()


if __name__ == "__main__":
    main()
