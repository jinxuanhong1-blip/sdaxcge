#!/usr/bin/env python3
"""Extract concordant-4 malignant-cell counts for one public cohort.

Observed pseudobulk uses every malignant cell in the locked unit (no QC),
matching the locked family OLS. The GRN subsample uses QC-pass malignant
cells from units with n_malignant >= 30, capped at 100 cells per unit.

Malignant definitions match the locked Seurat concordant-4 analysis:
  GSE123902, GSE189357: (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0
  GSE131907: author Cell_subtype == Malignant cells
  GSE205335 is extracted in R (RDS), not here.
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
LOCKED = ROOT / "data" / "locked_patient_units.tsv"
MARKER123 = ROOT / "data" / "GSE123902_marker_units.tsv"
QC_MIN_GENES = 200
QC_MIN_UMI = 500
QC_MAX_MT = 20
MIN_MAL_FOR_GRN = 30
CAP_PER_UNIT = 100
SEED = 1


def locked_units(dataset: str) -> pd.DataFrame:
    df = pd.read_csv(LOCKED, sep="\t")
    df = df.loc[df["dataset"] == dataset].copy()
    if df.empty:
        raise SystemExit(f"no locked units for {dataset}")
    return df


def collapse_rows(mat: sparse.csr_matrix, genes: np.ndarray):
    genes = np.asarray([str(g).upper() for g in genes])
    if pd.Index(genes).is_unique:
        return mat.tocsr(), genes
    uniq, inv = np.unique(genes, return_inverse=True)
    coo = mat.tocoo()
    out = sparse.coo_matrix(
        (coo.data, (inv[coo.row], coo.col)),
        shape=(len(uniq), mat.shape[1]),
    )
    return out.tocsr(), uniq


def row_vector(mat: sparse.spmatrix, index: dict, gene: str, n: int) -> np.ndarray:
    i = index.get(gene)
    if i is None:
        return np.zeros(n, dtype=np.float64)
    return np.asarray(mat.getrow(i).todense()).ravel()


def marker_masks(mat: sparse.csr_matrix, genes: np.ndarray):
    index = {g: i for i, g in enumerate(genes)}
    n = mat.shape[1]
    mal = (
        (row_vector(mat, index, "EPCAM", n) > 0)
        | (row_vector(mat, index, "KRT8", n) > 0)
        | (row_vector(mat, index, "KRT18", n) > 0)
        | (row_vector(mat, index, "KRT19", n) > 0)
    ) & (row_vector(mat, index, "PTPRC", n) == 0)
    tnk = (
        (row_vector(mat, index, "CD3D", n) > 0)
        | (row_vector(mat, index, "CD3E", n) > 0)
        | (row_vector(mat, index, "CD8A", n) > 0)
        | (row_vector(mat, index, "NKG7", n) > 0)
        | (row_vector(mat, index, "GNLY", n) > 0)
        | (row_vector(mat, index, "KLRD1", n) > 0)
    ) & (~mal)
    return mal, tnk


def qc_mask(mat: sparse.csr_matrix, genes: np.ndarray) -> np.ndarray:
    ncount = np.asarray(mat.sum(axis=0)).ravel()
    nfeat = np.asarray((mat > 0).sum(axis=0)).ravel()
    mt_idx = [i for i, g in enumerate(genes) if str(g).startswith("MT-")]
    if mt_idx:
        mt = np.asarray(mat[mt_idx].sum(axis=0)).ravel()
    else:
        mt = np.zeros(mat.shape[1])
    mt_pct = 100.0 * mt / np.maximum(ncount, 1.0)
    return (nfeat >= QC_MIN_GENES) & (ncount >= QC_MIN_UMI) & (mt_pct < QC_MAX_MT)


def cldn4_stats(mat: sparse.csr_matrix, genes: np.ndarray, mal: np.ndarray):
    index = {g: i for i, g in enumerate(genes)}
    x = row_vector(mat, index, "CLDN4", mat.shape[1])[mal]
    if x.size == 0:
        return np.nan, np.nan
    return float(100.0 * np.mean(x > 0)), float(np.mean(np.log1p(x)))


def take_cap(ix: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    ix = np.asarray(ix)
    if ix.size <= CAP_PER_UNIT:
        return np.sort(ix)
    return np.sort(rng.choice(ix, size=CAP_PER_UNIT, replace=False))


def align_hstack(pieces: list[tuple[np.ndarray, sparse.csr_matrix, np.ndarray]]):
    """pieces: (genes, csr genes x cells, unit_id per cell)."""
    all_genes = sorted(set().union(*[set(g.tolist()) for g, _, _ in pieces]))
    idx = {g: i for i, g in enumerate(all_genes)}
    blocks = []
    units = []
    for genes, mat, unit_ids in pieces:
        coo = mat.tocoo()
        mapped_r = np.fromiter((idx[genes[r]] for r in coo.row), dtype=np.int32, count=coo.nnz)
        mapped = sparse.coo_matrix(
            (coo.data, (mapped_r, coo.col)),
            shape=(len(all_genes), mat.shape[1]),
        ).tocsr()
        blocks.append(mapped)
        units.append(np.asarray(unit_ids))
    mat = sparse.hstack(blocks, format="csr")
    return np.asarray(all_genes), mat, np.concatenate(units)


def write_outputs(out: Path, dataset: str, units: pd.DataFrame, genes: list[str], pb: np.ndarray, grn_genes, grn_mat, grn_units):
    out.mkdir(parents=True, exist_ok=True)
    units.to_csv(out / "units.tsv", sep="\t", index=False)
    pb_df = pd.DataFrame(pb, index=genes, columns=units["unit_id"].tolist())
    pb_df.index.name = "gene"
    with gzip.open(out / "pseudobulk.tsv.gz", "wt") as fh:
        pb_df.to_csv(fh, sep="\t")
    spio.mmwrite(str(out / "grn_matrix.mtx"), grn_mat, field="real")
    pd.Series(grn_genes).to_csv(out / "grn_genes.tsv", index=False, header=False)
    pd.DataFrame({"unit_id": grn_units}).to_csv(out / "grn_cells.tsv", sep="\t", index=False)
    print(
        f"wrote {dataset}: units={len(units)} pb_genes={len(genes)} "
        f"grn={grn_mat.shape[0]}x{grn_mat.shape[1]} nnz={grn_mat.nnz}",
        flush=True,
    )


def dense_csv_to_csr(path: Path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [str(c).upper() for c in df.columns]
    if df.columns.has_duplicates:
        df = df.T.groupby(level=0).sum().T
    genes = df.columns.to_numpy()
    mat = sparse.csr_matrix(df.to_numpy(dtype=np.float64).T)
    del df
    return mat, genes


def extract_gse123902(geo: Path, out: Path):
    dataset = "GSE123902"
    locked = locked_units(dataset)
    marker = pd.read_csv(MARKER123, sep="\t")
    merged = locked.merge(
        marker,
        left_on=["unit_id", "tissue"],
        right_on=["patient", "tissue"],
        how="left",
        suffixes=("", "_marker"),
    )
    if merged["file"].isna().any():
        missing = merged.loc[merged["file"].isna(), "unit_id"].tolist()
        raise SystemExit(f"GSE123902 file map missing for {missing}")
    tar_path = geo / "GSE123902_RAW.tar"
    dest = geo / "gse123902"
    dest.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    unit_rows = []
    pb_pieces = []
    grn_pieces = []
    with tarfile.open(tar_path) as tar:
        for rec in merged.sort_values("unit_id").itertuples(index=False):
            fn = rec.file
            local = dest / fn
            if not local.exists():
                tar.extract(fn, path=dest)
            print(f"  {rec.unit_id} {fn}", flush=True)
            mat, genes = dense_csv_to_csr(local)
            mal, tnk = marker_masks(mat, genes)
            pct, mean = cldn4_stats(mat, genes, mal)
            qc = qc_mask(mat, genes)
            mal_qc = np.flatnonzero(mal & qc)
            n_mal = int(mal.sum())
            if n_mal >= MIN_MAL_FOR_GRN and mal_qc.size:
                keep = take_cap(mal_qc, rng)
                grn_pieces.append((genes, mat[:, keep], np.array([rec.unit_id] * len(keep))))
                n_grn = int(len(keep))
            else:
                n_grn = 0
            pb_pieces.append((genes, np.asarray(mat[:, mal].sum(axis=1)).ravel(), rec.unit_id))
            unit_rows.append(
                {
                    "dataset": dataset,
                    "unit_id": rec.unit_id,
                    "n_cells": int(mat.shape[1]),
                    "n_malignant": n_mal,
                    "n_tnk": int(tnk.sum()),
                    "n_malignant_qc": int(mal_qc.size),
                    "n_grn": n_grn,
                    "mal_CLDN4_pct": pct,
                    "mal_CLDN4_mean": mean,
                }
            )
            del mat
    all_genes = sorted(set().union(*[set(g.tolist()) for g, _, _ in pb_pieces]))
    gix = {g: i for i, g in enumerate(all_genes)}
    pb = np.zeros((len(all_genes), len(pb_pieces)), dtype=np.float64)
    for j, (genes, vec, _unit) in enumerate(pb_pieces):
        for g, v in zip(genes, vec):
            pb[gix[g], j] = v
    units = pd.DataFrame(unit_rows)
    grn_genes, grn_mat, grn_units = align_hstack(grn_pieces)
    write_outputs(out / dataset, dataset, units, all_genes, pb, grn_genes, grn_mat, grn_units)


def extract_gse189357(geo: Path, out: Path):
    dataset = "GSE189357"
    locked = locked_units(dataset).sort_values("unit_id")
    tar_path = geo / "GSE189357_RAW.tar"
    dest = geo / "gse189357"
    dest.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    unit_rows = []
    pb_pieces = []
    grn_pieces = []
    names = {m.name for m in tarfile.open(tar_path).getmembers()}
    with tarfile.open(tar_path) as tar:
        for unit in locked["unit_id"]:
            print(f"  {unit}", flush=True)
            paths = {}
            for kind, suffix in (("mtx", "_matrix.mtx.gz"), ("feat", "_features.tsv.gz"), ("bc", "_barcodes.tsv.gz")):
                hits = [n for n in names if n.endswith(f"_{unit}{suffix}")]
                if len(hits) != 1:
                    raise SystemExit(f"expected one {suffix} for {unit}, found {hits}")
                local = dest / hits[0]
                if not local.exists():
                    tar.extract(hits[0], path=dest)
                paths[kind] = local
            mat = spio.mmread(paths["mtx"]).tocsr()
            feat = pd.read_csv(paths["feat"], sep="\t", header=None)
            genes = feat.iloc[:, 1].astype(str).to_numpy() if feat.shape[1] > 1 else feat.iloc[:, 0].astype(str).to_numpy()
            mat, genes = collapse_rows(mat, genes)
            mal, tnk = marker_masks(mat, genes)
            pct, mean = cldn4_stats(mat, genes, mal)
            qc = qc_mask(mat, genes)
            mal_qc = np.flatnonzero(mal & qc)
            n_mal = int(mal.sum())
            if n_mal >= MIN_MAL_FOR_GRN and mal_qc.size:
                keep = take_cap(mal_qc, rng)
                grn_pieces.append((genes, mat[:, keep], np.array([unit] * len(keep))))
                n_grn = int(len(keep))
            else:
                n_grn = 0
            pb_pieces.append((genes, np.asarray(mat[:, mal].sum(axis=1)).ravel(), unit))
            unit_rows.append(
                {
                    "dataset": dataset,
                    "unit_id": unit,
                    "n_cells": int(mat.shape[1]),
                    "n_malignant": n_mal,
                    "n_tnk": int(tnk.sum()),
                    "n_malignant_qc": int((mal & qc).sum()),
                    "n_grn": n_grn,
                    "mal_CLDN4_pct": pct,
                    "mal_CLDN4_mean": mean,
                }
            )
            del mat
    # align pseudobulk genes
    all_genes = sorted(set().union(*[set(g.tolist()) for g, _, _ in pb_pieces]))
    idx = {g: i for i, g in enumerate(all_genes)}
    pb = np.zeros((len(all_genes), len(pb_pieces)), dtype=np.float64)
    for j, (genes, vec, _unit) in enumerate(pb_pieces):
        for g, v in zip(genes, vec):
            pb[idx[g], j] = v
    units = pd.DataFrame(unit_rows)
    grn_genes, grn_mat, grn_units = align_hstack(grn_pieces)
    write_outputs(out / dataset, dataset, units, all_genes, pb, grn_genes, grn_mat, grn_units)


def extract_gse131907(geo: Path, out: Path):
    dataset = "GSE131907"
    locked = locked_units(dataset)
    keep_samples = set(locked["unit_id"])
    ann_path = geo / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = geo / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cells = []
    with gzip.open(ann_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            sample = p[ix["Sample"]]
            if sample not in keep_samples:
                continue
            cells.append(
                {
                    "index": p[ix["Index"]],
                    "sample": sample,
                    "malignant": p[ix["Cell_subtype"]] == "Malignant cells",
                    "tnk": p[ix["Cell_type"]] in {"T lymphocytes", "NK cells"},
                }
            )
    by_sample_n = {}
    by_sample_tnk = {}
    by_sample_mal = {}
    for sample, sub in pd.DataFrame(cells).groupby("sample"):
        by_sample_n[sample] = int(len(sub))
        by_sample_tnk[sample] = int(sub["tnk"].sum())
        by_sample_mal[sample] = int(sub["malignant"].sum())
    mal_index = {c["index"] for c in cells if c["malignant"]}
    print(f"  annotation malignant in locked samples: {len(mal_index)}", flush=True)

    with gzip.open(mat_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    barcodes = header[1:]
    col_of = {b: i for i, b in enumerate(barcodes)}
    missing = [b for b in mal_index if b not in col_of]
    if missing:
        raise SystemExit(f"GSE131907 malignant barcodes missing from matrix: {len(missing)}")
    sample_of = {c["index"]: c["sample"] for c in cells if c["malignant"]}
    samples = sorted(keep_samples)
    sample_code = {s: i for i, s in enumerate(samples)}
    sel = np.array([col_of[b] for b in barcodes if b in mal_index], dtype=np.int64)
    # barcodes list is aligned with matrix columns; rebuild sel in barcode order
    mal_positions = [i for i, b in enumerate(barcodes) if b in mal_index]
    sel = np.asarray(mal_positions, dtype=np.int64)
    code = np.asarray([sample_code[sample_of[barcodes[i]]] for i in sel], dtype=np.int32)
    n_mal = len(sel)
    n_units = len(samples)
    print(f"  matrix malignant columns: {n_mal}", flush=True)

    ncount = np.zeros(n_mal, dtype=np.float64)
    nfeat = np.zeros(n_mal, dtype=np.int32)
    mt = np.zeros(n_mal, dtype=np.float64)
    cldn4 = np.zeros(n_mal, dtype=np.float64)
    genes = []
    pb_rows = []
    with gzip.open(mat_path, "rt") as fh:
        fh.readline()
        for gi, line in enumerate(fh):
            line = line.rstrip("\n")
            gene, rest = line.split("\t", 1)
            gene = gene.upper()
            arr = np.fromstring(rest, sep="\t", dtype=np.float64)
            if arr.size != len(barcodes):
                raise SystemExit(f"short row {gene}: {arr.size}")
            vals = arr[sel]
            pb_rows.append(np.bincount(code, weights=vals, minlength=n_units))
            ncount += vals
            nfeat += (vals > 0).astype(np.int32)
            if gene.startswith("MT-"):
                mt += vals
            if gene == "CLDN4":
                cldn4 = vals
            genes.append(gene)
            if gi and gi % 4000 == 0:
                print(f"  pass1 {gi} genes", flush=True)
    pb = np.vstack(pb_rows)
    if len(set(genes)) != len(genes):
        raise SystemExit("GSE131907 gene symbols are not unique; collapse is not implemented")
    genes = list(genes)
    mt_pct = 100.0 * mt / np.maximum(ncount, 1.0)
    qc = (nfeat >= QC_MIN_GENES) & (ncount >= QC_MIN_UMI) & (mt_pct < QC_MAX_MT)
    # per-unit stats from the CLDN4 vector and counts
    unit_rows = []
    rng = np.random.default_rng(SEED)
    chosen_local = []
    for s in samples:
        local = np.flatnonzero(code == sample_code[s])
        x = cldn4[local]
        n_m = int(local.size)
        qc_local = local[qc[local]]
        if n_m >= MIN_MAL_FOR_GRN and qc_local.size:
            keep = take_cap(qc_local, rng)
            chosen_local.append(keep)
            n_grn = int(len(keep))
        else:
            n_grn = 0
        unit_rows.append(
            {
                "dataset": dataset,
                "unit_id": s,
                "n_cells": by_sample_n[s],
                "n_malignant": n_m,
                "n_tnk": by_sample_tnk[s],
                "n_malignant_qc": int(qc[local].sum()),
                "n_grn": n_grn,
                "mal_CLDN4_pct": float(100.0 * np.mean(x > 0)) if n_m else np.nan,
                "mal_CLDN4_mean": float(np.mean(np.log1p(x))) if n_m else np.nan,
            }
        )
        if n_m != by_sample_mal[s]:
            raise SystemExit(f"{s} malignant annotation {by_sample_mal[s]} != matrix {n_m}")
    if not chosen_local:
        raise SystemExit("no GRN cells for GSE131907")
    chosen = np.concatenate(chosen_local)
    chosen_cols = sel[chosen]
    grn_units = np.asarray([samples[c] for c in code[chosen]])
    print(f"  pass2 extracting {len(chosen)} GRN cells", flush=True)
    rows_i, cols_j, data = [], [], []
    with gzip.open(mat_path, "rt") as fh:
        fh.readline()
        for gi, line in enumerate(fh):
            gene, rest = line.rstrip("\n").split("\t", 1)
            arr = np.fromstring(rest, sep="\t", dtype=np.float64)
            vals = arr[chosen_cols]
            nz = np.flatnonzero(vals)
            if nz.size:
                rows_i.append(np.full(nz.size, gi, dtype=np.int32))
                cols_j.append(nz.astype(np.int32))
                data.append(vals[nz])
            if gi and gi % 4000 == 0:
                print(f"  pass2 {gi} genes", flush=True)
    grn = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows_i), np.concatenate(cols_j))),
        shape=(len(genes) if isinstance(genes, list) else pb.shape[0], len(chosen)),
    ).tocsr()
    # If duplicate genes were collapsed in pb, grn rows were not. GSE131907 symbols
    # are unique in practice; collapse if needed.
    gene_arr = np.asarray(genes)
    if not pd.Index(gene_arr).is_unique:
        grn, gene_arr = collapse_rows(grn, gene_arr)
    units = pd.DataFrame(unit_rows)
    # column order of pb matches `samples`
    units = units.set_index("unit_id").loc[samples].reset_index()
    write_outputs(out / dataset, dataset, units, list(gene_arr), pb, gene_arr, grn, grn_units)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["GSE123902", "GSE131907", "GSE189357"])
    ap.add_argument("--geo", type=Path, default=Path("/tmp/geo_c4"))
    ap.add_argument("--out", type=Path, default=Path("/tmp/c4_work"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.dataset == "GSE123902":
        extract_gse123902(args.geo, args.out)
    elif args.dataset == "GSE189357":
        extract_gse189357(args.geo, args.out)
    else:
        extract_gse131907(args.geo, args.out)


if __name__ == "__main__":
    main()
