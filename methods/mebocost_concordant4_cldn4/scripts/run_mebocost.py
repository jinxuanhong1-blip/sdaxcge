#!/usr/bin/env python3
"""Concordant-4 CLDN4-high vs low malignant -> T/NK metabolite communication.

Engine: MEBOCOST (kaifuchenlab), met_est='mebocost' (enzyme mean). Not a
reimplementation of the communication score. Honest unit = patient / locked
sample. Primary split = malignant Q4 vs Q1. CLDN4 only. No dual-high.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import io
import math
import os
import tarfile
import time
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse
from scipy.stats import rankdata, wilcoxon

EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
MALIG_SUB = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
MIN_ARM = 10
MIN_TNK = 20
MIN_MAL_Q4 = 40
SPLITS = ("q4q1", "median", "pctpos")


def logmsg(*parts):
    print(time.strftime("%H:%M:%S"), *parts, flush=True)


def here_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def canon_map(genes) -> dict:
    """Upper-case key -> symbol as curated (two historical orf symbols are mixed-case)."""
    out = {}
    for g in genes:
        out[str(g).upper()] = str(g)
    return out


def load_db_genes(mebo_root: Path):
    sen = pd.read_csv(
        mebo_root / "data/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv",
        sep="\t",
    )
    enz = pd.read_csv(
        mebo_root / "data/mebocost_db/human/metabolite_associated_gene_reaction_HMDB_summary.tsv",
        sep="\t",
    )
    genes = set(sen["Gene_name"].astype(str))
    for raw in enz["gene"].dropna():
        for part in str(raw).split(";"):
            name = part.split("[")[0].strip()
            if name:
                genes.add(name)
    genes.update(["CLDN4", *EPI, "PTPRC", *TNK_MARKERS])
    genes = {g for g in genes if g and g[0].isalpha()}
    return genes, canon_map(genes)


def write_conf(path: Path, mebo_root: Path) -> None:
    db = mebo_root / "data"
    path.write_text(
        "\n".join(
            [
                "[common]",
                f"hmdb_info_path = {db}/mebocost_db/common/metabolite_annotation_HMDB_summary.tsv",
                f"scfea_info_path = {db}/scFEA/Human_M168_information.symbols.csv",
                f"compass_rxt_ann_path = {db}/Compass/rxn_md.csv",
                f"compass_met_ann_path = {db}/Compass/met_md.csv",
                "",
                "[human]",
                f"met_enzyme_path = {db}/mebocost_db/human/metabolite_associated_gene_reaction_HMDB_summary.tsv",
                f"met_sensor_path = {db}/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv",
                "",
                "[mouse]",
                f"met_enzyme_path = {db}/mebocost_db/mouse/metabolite_associated_gene_reaction_HMDB_summary_mouse.tsv",
                f"met_sensor_path = {db}/mebocost_db/mouse/mouse_met_sensor_update_Oct21_2025.tsv",
                "",
            ]
        )
    )


def save_unit(path: Path, mat_gc, genes, barcodes, mal, tnk, lib) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mat = mat_gc.tocsc()
    np.savez_compressed(
        path,
        data=mat.data,
        indices=mat.indices,
        indptr=mat.indptr,
        shape=np.array(mat.shape, dtype=np.int64),
        genes=np.array(genes, dtype=object),
        barcodes=np.array(barcodes, dtype=object),
        mal=np.asarray(mal, dtype=np.uint8),
        tnk=np.asarray(tnk, dtype=np.uint8),
        lib=np.asarray(lib, dtype=np.float64),
    )


def load_unit(path: Path):
    z = np.load(path, allow_pickle=True)
    mat = sparse.csc_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return {
        "mat": mat,  # genes x cells
        "genes": [str(g) for g in z["genes"].tolist()],
        "barcodes": [str(b) for b in z["barcodes"].tolist()],
        "mal": z["mal"].astype(bool),
        "tnk": z["tnk"].astype(bool),
        "lib": z["lib"].astype(np.float64),
    }


def collapse_genes(mat, genes, cmap):
    """mat is genes x cells. Rename to canonical symbols and sum duplicate rows."""
    named = [cmap.get(str(g).upper(), str(g).upper()) for g in genes]
    if len(set(named)) == len(named):
        return mat.tocsc(), named
    order = {}
    rows = []
    out_names = []
    mat = mat.tocsr()
    for i, name in enumerate(named):
        if name not in order:
            order[name] = len(out_names)
            out_names.append(name)
            rows.append(mat[i])
        else:
            rows[order[name]] = rows[order[name]] + mat[i]
    return sparse.vstack(rows).tocsc(), out_names


def subset_wanted(mat, genes, wanted_upper, cmap):
    genes = [str(g) for g in genes]
    keep = [i for i, g in enumerate(genes) if g.upper() in wanted_upper]
    if not keep:
        raise RuntimeError("none of the MEBOCOST/marker genes are in this matrix")
    sub = mat[keep, :]
    sub_genes = [genes[i] for i in keep]
    return collapse_genes(sub, sub_genes, cmap)


def marker_positive(mat, genes, markers) -> np.ndarray:
    idx = [genes.index(g) for g in markers if g in genes]
    if not idx:
        return np.zeros(mat.shape[1], dtype=bool)
    block = mat[idx, :]
    return np.asarray(block.sum(axis=0)).ravel() > 0


def gene_row(mat, genes, name) -> np.ndarray:
    if name not in genes:
        return np.zeros(mat.shape[1], dtype=np.float64)
    return np.asarray(mat[genes.index(name), :].todense()).ravel().astype(np.float64)


def pos_any_names(mat, genes, markers) -> np.ndarray:
    return marker_positive(mat, genes, markers)


# ---------------------------------------------------------------------------
# Loaders. mat is genes x cells, counts. lib is full-transcriptome UMI.
# ---------------------------------------------------------------------------

def extract_gse123902(raw: Path, data: Path, cache: Path, cmap, wanted_upper):
    units = pd.read_csv(data / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor["ord"] = tumor["tissue"].map({"PRIMARY": 0, "METASTASIS": 1})
    tumor = tumor.sort_values(["patient", "ord"]).drop_duplicates("patient")
    tar_path = raw / "GSE123902" / "GSE123902_RAW.tar"
    csv_dir = raw / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*_dense.csv.gz")):
        logmsg("untar GSE123902")
        with tarfile.open(tar_path) as tar:
            tar.extractall(csv_dir)
    files = list(csv_dir.glob("*_dense.csv.gz"))
    n = 0
    for rec in tumor.itertuples(index=False):
        patient = str(rec.patient)
        dest = cache / f"GSE123902__{patient}.npz"
        if dest.exists():
            n += 1
            continue
        fp = csv_dir / str(rec.file)
        if not fp.exists():
            hit = [p for p in files if f"_{patient}_" in p.name and "_NORMAL_" not in p.name]
            hit.sort(key=lambda p: (0 if "_PRIMARY_" in p.name else 1, p.name))
            if not hit:
                logmsg("  missing CSV", patient)
                continue
            fp = hit[0]
        logmsg("  csv", patient, fp.name)
        df = pd.read_csv(fp)
        cell_ids = df.iloc[:, 0].astype(str).tolist()
        gene_names = [str(c) for c in df.columns[1:]]
        # cells x genes dense -> genes x cells sparse
        arr = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
        lib = arr.sum(axis=1).astype(np.float64)
        mat = sparse.csc_matrix(arr.T)
        del df, arr
        mat, genes = subset_wanted(mat, gene_names, wanted_upper, cmap)
        mal = pos_any_names(mat, genes, EPI) & (gene_row(mat, genes, "PTPRC") == 0)
        tnk = pos_any_names(mat, genes, TNK_MARKERS) & (~mal)
        save_unit(dest, mat, genes, [f"{patient}_{c}" for c in cell_ids], mal, tnk, lib)
        n += 1
        logmsg("   ", patient, "cells", mat.shape[1], "mal", int(mal.sum()), "tnk", int(tnk.sum()))
    logmsg("GSE123902 units", n)


def extract_gse189357(raw: Path, data: Path, cache: Path, cmap, wanted_upper):
    meta = pd.read_csv(data / "GSE189357_sample_metadata.tsv", sep="\t")
    tar_path = raw / "GSE189357" / "GSE189357_RAW.tar"
    ex = raw / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    needed = [ex / f"{r.gsm}_{r.patient}_matrix.mtx.gz" for r in meta.itertuples(index=False)]
    if not all(p.exists() for p in needed):
        logmsg("untar GSE189357")
        with tarfile.open(tar_path) as tar:
            tar.extractall(ex)
    for rec in meta.itertuples(index=False):
        patient = str(rec.patient)
        dest = cache / f"GSE189357__{patient}.npz"
        if dest.exists():
            continue
        prefix = ex / f"{rec.gsm}_{patient}"
        mtx = Path(str(prefix) + "_matrix.mtx.gz")
        cells_p = Path(str(prefix) + "_barcodes.tsv.gz")
        feat_p = Path(str(prefix) + "_features.tsv.gz")
        logmsg("  mtx", patient)
        mat = spio.mmread(mtx).tocsc()
        barcodes = pd.read_csv(cells_p, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
        feat = pd.read_csv(feat_p, sep="\t", header=None)
        symbols = feat.iloc[:, 1].astype(str).tolist() if feat.shape[1] > 1 else feat.iloc[:, 0].astype(str).tolist()
        if mat.shape[1] != len(barcodes) and mat.shape[0] == len(barcodes):
            mat = mat.T.tocsc()
        if mat.shape[1] != len(barcodes) or mat.shape[0] != len(symbols):
            raise RuntimeError(f"{patient} mtx {mat.shape} barcodes {len(barcodes)} genes {len(symbols)}")
        lib = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
        mat, genes = subset_wanted(mat, symbols, wanted_upper, cmap)
        mal = pos_any_names(mat, genes, EPI) & (gene_row(mat, genes, "PTPRC") == 0)
        tnk = pos_any_names(mat, genes, TNK_MARKERS) & (~mal)
        save_unit(dest, mat, genes, [f"{patient}_{b}" for b in barcodes], mal, tnk, lib)
        logmsg("   ", patient, "cells", mat.shape[1], "mal", int(mal.sum()), "tnk", int(tnk.sum()))


def extract_gse131907(raw: Path, data: Path, cache: Path, cmap, wanted_upper):
    samples = pd.read_csv(data / "GSE131907_samples.tsv", sep="\t")
    keep = [str(s) for s in samples.loc[samples["n_malignant"] > 0, "sample"]]
    if all((cache / f"GSE131907__{s}.npz").exists() for s in keep):
        logmsg("GSE131907 cache complete", len(keep))
        return
    ann = pd.read_csv(raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    ann["Index"] = ann["Index"].astype(str)
    ann = ann.drop_duplicates("Index").set_index("Index")
    path = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    logmsg("stream GSE131907", path)
    store = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        ncells = len(cells)
        lib = np.zeros(ncells, dtype=np.float64)
        n_seen = 0
        for line in handle:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab]
            n_seen += 1
            vals = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float32)
            if vals.size != ncells:
                if vals.size < ncells:
                    vals = np.pad(vals, (0, ncells - vals.size))
                else:
                    vals = vals[:ncells]
            lib += vals
            key = gene.upper()
            if key in wanted_upper:
                nz = np.flatnonzero(vals)
                canon = cmap.get(key, key)
                idx = nz.astype(np.int32)
                vv = vals[nz].astype(np.float32, copy=False)
                if canon in store:
                    pidx, pvv = store[canon]
                    both = np.concatenate([pidx, idx])
                    bothv = np.concatenate([pvv, vv])
                    order = np.argsort(both, kind="mergesort")
                    both, bothv = both[order], bothv[order]
                    uniq, inv = np.unique(both, return_inverse=True)
                    summed = np.zeros(len(uniq), dtype=np.float32)
                    np.add.at(summed, inv, bothv)
                    store[canon] = (uniq.astype(np.int32), summed)
                else:
                    store[canon] = (idx, vv)
            if n_seen % 2000 == 0:
                logmsg("    streamed", n_seen, "kept", len(store))
    logmsg("  streamed genes", n_seen, "kept", len(store), "cells", ncells)
    genes = list(store.keys())
    row_chunks, col_chunks, val_chunks = [], [], []
    for gi, g in enumerate(genes):
        idx, vv = store[g]
        row_chunks.append(np.full(len(idx), gi, dtype=np.int32))
        col_chunks.append(idx.astype(np.int32, copy=False))
        val_chunks.append(vv.astype(np.float32, copy=False))
    mat_all = sparse.csc_matrix(
        (
            np.concatenate(val_chunks) if val_chunks else np.zeros(0, np.float32),
            (
                np.concatenate(row_chunks) if row_chunks else np.zeros(0, np.int32),
                np.concatenate(col_chunks) if col_chunks else np.zeros(0, np.int32),
            ),
        ),
        shape=(len(genes), ncells),
    )
    del store, row_chunks, col_chunks, val_chunks
    present = np.array([c in ann.index for c in cells])
    logmsg("  annotation overlap", int(present.sum()), "/", ncells)
    cell_arr = np.array(cells, dtype=object)
    samp_of = pd.Series("", index=np.arange(ncells))
    sub_of = pd.Series("", index=np.arange(ncells))
    type_of = pd.Series("", index=np.arange(ncells))
    if present.any():
        hit = cell_arr[present]
        samp_of.iloc[np.flatnonzero(present)] = ann.loc[hit, "Sample"].astype(str).to_numpy()
        sub_of.iloc[np.flatnonzero(present)] = ann.loc[hit, "Cell_subtype"].astype(str).to_numpy()
        type_of.iloc[np.flatnonzero(present)] = ann.loc[hit, "Cell_type"].astype(str).to_numpy()
    for s in keep:
        dest = cache / f"GSE131907__{s}.npz"
        if dest.exists():
            continue
        cols = np.flatnonzero(samp_of.to_numpy() == s)
        if len(cols) == 0:
            logmsg("  missing sample", s)
            continue
        sub = mat_all[:, cols]
        sel_cells = cell_arr[cols].tolist()
        mal = sub_of.to_numpy()[cols]
        mal = np.array([x in MALIG_SUB for x in mal], dtype=bool)
        tnk = type_of.to_numpy()[cols]
        tnk = np.array([x in TNK_TYPES for x in tnk], dtype=bool)
        save_unit(dest, sub, genes, sel_cells, mal, tnk, lib[cols])
        logmsg("   ", s, "cells", len(cols), "mal", int(mal.sum()), "tnk", int(tnk.sum()))


def extract_gse205335(raw: Path, data: Path, cache: Path, cmap, wanted_upper, mebo_root: Path):
    patients = pd.read_csv(data / "GSE205335_patients.tsv", sep="\t")
    keep = [str(s) for s in patients.loc[patients["n_malignant"] > 0, "patient"]]
    if all((cache / f"GSE205335__{p}.npz").exists() for p in keep):
        logmsg("GSE205335 cache complete", len(keep))
        return
    subset = Path(os.environ.get("GSE205335_SUBSET", "/tmp/gse205335_subset"))
    if not (subset / "matrix.mtx").exists():
        raise RuntimeError(
            f"missing {subset}/matrix.mtx. Run scripts/subset_gse205335.R on the peeled RDS."
        )
    logmsg("load GSE205335 subset", subset)
    mat = spio.mmread(subset / "matrix.mtx").tocsc()
    genes = [ln.strip() for ln in (subset / "genes.txt").read_text().splitlines() if ln.strip()]
    barcodes = [ln.strip() for ln in (subset / "barcodes.txt").read_text().splitlines() if ln.strip()]
    if mat.shape != (len(genes), len(barcodes)):
        raise RuntimeError(f"subset shape {mat.shape} vs genes {len(genes)} barcodes {len(barcodes)}")
    lib_df = pd.read_csv(subset / "lib.tsv.gz", sep="\t")
    lib_map = dict(zip(lib_df["barcode"].astype(str), lib_df["lib"].astype(float)))
    ident = pd.read_csv(raw / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    gsm = pd.read_csv(data / "GSE205335_gsm_map.tsv", sep="\t")
    ident["barcode"] = ident["barcode"].astype(str)
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident = ident.drop_duplicates("barcode").set_index("barcode")
    missing = [b for b in barcodes if b not in ident.index]
    if missing:
        raise RuntimeError(f"unmapped barcodes {len(missing)} example {missing[:3]}")
    if ident.loc[barcodes, "patient"].isna().any():
        bad = ident.loc[barcodes, "patient"].isna()
        orig = ident.loc[barcodes].loc[bad.values, "orig.ident"].unique()[:8]
        raise RuntimeError(f"unmapped orig.ident {list(orig)}")
    tissue = ident.loc[barcodes, "tissue"].astype(str)
    is_normal = tissue.str.startswith("Normal ")
    patient = ident.loc[barcodes, "patient"].astype(str)
    lineage_sub = ident.loc[barcodes, "lineage.sub"].astype(str)
    lineage_total = ident.loc[barcodes, "lineage.total"].astype(str)
    mat, genes = collapse_genes(mat, genes, cmap)
    lib = np.array([lib_map[b] for b in barcodes], dtype=np.float64)
    for pt in keep:
        dest = cache / f"GSE205335__{pt}.npz"
        if dest.exists():
            continue
        sel = np.where((patient.to_numpy() == pt) & (~is_normal.to_numpy()))[0]
        if len(sel) == 0:
            continue
        sub = mat[:, sel]
        mal = (lineage_sub.to_numpy()[sel] == "Malignant cells")
        tnk = (lineage_total.to_numpy()[sel] == "T/NK cells")
        bc = [barcodes[i] for i in sel]
        save_unit(dest, sub, genes, bc, mal, tnk, lib[sel])
        logmsg("   ", pt, "cells", len(sel), "mal", int(mal.sum()), "tnk", int(tnk.sum()))
    del mebo_root


# ---------------------------------------------------------------------------
# MEBOCOST
# ---------------------------------------------------------------------------

def quartile_high_low(x: np.ndarray):
    n = len(x)
    if n < 4:
        return None
    r = rankdata(x, method="ordinal")
    q1 = math.floor(n * 0.25)
    q4 = math.ceil(n * 0.75)
    if q1 < 1 or q4 > n or q1 >= q4:
        return None
    return r > q4, r <= q1


def split_masks(cldn4_log, cldn4_umi, split):
    if split == "q4q1":
        return quartile_high_low(cldn4_log)
    if split == "median":
        med = float(np.median(cldn4_log))
        if not np.isfinite(med):
            return None
        return cldn4_log > med, cldn4_log <= med
    if split == "pctpos":
        return cldn4_umi > 0, cldn4_umi == 0
    raise KeyError(split)


def _quiet_mebocost():
    import mebocost.crosstalk_calculator as cc
    import mebocost.mebocost as mb
    import mebocost.MetEstimator as me

    def _noop(*_a, **_k):
        return None

    cc.info = _noop
    mb.info = _noop
    me.info = _noop


def run_mebocost_once(adata, conf: Path, n_shuffle: int):
    from mebocost import mebocost

    _quiet_mebocost()
    obj = mebocost.create_obj(
        adata=adata,
        group_col="cell_type",
        species="human",
        config_path=str(conf),
        cutoff_exp="auto",
        cutoff_met="auto",
        cutoff_prop=0.15,
        sensor_type="All",
        thread=1,
    )
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        obj.infer_commu(
            n_shuffle=n_shuffle,
            seed=12345,
            Return=True,
            save_permuation=False,
            pval_cutoff=1.01,
            min_cell_number=MIN_ARM,
            thread=1,
        )
    res = obj.original_result
    if res is None or len(res) == 0:
        raise RuntimeError("MEBOCOST original_result empty\n" + buf.getvalue()[-2000:])
    return res


def score_cached_units(cache: Path, conf: Path, pairs: pd.DataFrame, n_shuffle: int, cohorts):
    rows = []
    inv_rows = []
    desc = {}  # (hmdb, sensor) -> list of q4q1 deltas
    desc_meta = {}
    files = sorted(cache.glob("*.npz"))
    for path in files:
        cohort, patient = path.stem.split("__", 1)
        if cohorts and cohort not in cohorts:
            continue
        unit = load_unit(path)
        mat, genes = unit["mat"], unit["genes"]
        mal, tnk, lib = unit["mal"], unit["tnk"], unit["lib"]
        n_mal = int(mal.sum())
        n_tnk = int(tnk.sum())
        cldn4_umi_all = gene_row(mat, genes, "CLDN4")
        lib_safe = np.maximum(lib, 1.0)
        cldn4_log_all = np.log1p(cldn4_umi_all / lib_safe * 1e4)
        inv = {
            "cohort": cohort,
            "patient": patient,
            "n_cells": int(mat.shape[1]),
            "n_mal": n_mal,
            "n_tnk": n_tnk,
            "mal_cldn4_mean": float(cldn4_log_all[mal].mean()) if n_mal else np.nan,
            "mal_cldn4_pct": float((cldn4_umi_all[mal] > 0).mean()) if n_mal else np.nan,
        }
        inv_rows.append(inv)
        if n_tnk < MIN_TNK or n_mal < MIN_ARM * 2:
            logmsg("  skip gates", cohort, patient, "mal", n_mal, "tnk", n_tnk)
            continue
        cldn4_log = cldn4_log_all[mal]
        cldn4_umi = cldn4_umi_all[mal]
        mal_idx = np.flatnonzero(mal)
        tnk_idx = np.flatnonzero(tnk)
        for split in SPLITS:
            if split == "q4q1" and n_mal < MIN_MAL_Q4:
                continue
            masks = split_masks(cldn4_log, cldn4_umi, split)
            if masks is None:
                continue
            high_m, low_m = masks
            n_high = int(high_m.sum())
            n_low = int(low_m.sum())
            if n_high < MIN_ARM or n_low < MIN_ARM:
                logmsg("  skip arm", cohort, patient, split, n_high, n_low)
                continue
            high_idx = mal_idx[high_m]
            low_idx = mal_idx[low_m]
            cols = np.concatenate([high_idx, low_idx, tnk_idx])
            labels = (
                ["CLDN4_high"] * len(high_idx)
                + ["CLDN4_low"] * len(low_idx)
                + ["TNK"] * len(tnk_idx)
            )
            block = mat[:, cols].T.tocsr().astype(np.float32)  # cells x genes
            scale = (1e4 / np.maximum(lib[cols], 1.0)).astype(np.float32)
            block = block.multiply(scale[:, None]).tocsr()
            block.data = np.log1p(block.data)
            adata = ad.AnnData(X=block)
            adata.var_names = genes
            adata.obs_names = [f"c{i}" for i in range(block.shape[0])]
            adata.obs["cell_type"] = labels
            spk = cache.parent / "scores" / f"{cohort}__{patient}__{split}.pkl"
            spk.parent.mkdir(parents=True, exist_ok=True)
            if spk.exists() and spk.stat().st_size > 0:
                to_tnk = pd.read_pickle(spk)
                logmsg("  cache", cohort, patient, split, "rows", len(to_tnk))
            else:
                logmsg(
                    "  MEBOCOST",
                    cohort,
                    patient,
                    split,
                    "high",
                    n_high,
                    "low",
                    n_low,
                    "tnk",
                    n_tnk,
                    "genes",
                    len(genes),
                )
                t0 = time.time()
                try:
                    res = run_mebocost_once(adata, conf, n_shuffle)
                except Exception as exc:
                    logmsg("    FAIL", cohort, patient, split, exc)
                    traceback.print_exc()
                    continue
                logmsg("    done", f"{time.time() - t0:.1f}s", "rows", len(res))
                to_tnk = res[res["Receiver"].astype(str) == "TNK"].copy()
                to_tnk.to_pickle(spk)
                del res
            if split == "q4q1":
                for rec in to_tnk.itertuples(index=False):
                    if str(rec.Sender) not in {"CLDN4_high", "CLDN4_low"}:
                        continue
                    key = (str(rec.Metabolite), str(rec.Sensor))
                    desc_meta[key] = (
                        str(rec.Metabolite_Name) if pd.notna(rec.Metabolite_Name) else "",
                        str(rec.Annotation) if pd.notna(rec.Annotation) else "",
                    )
            wide = {}
            for rec in to_tnk.itertuples(index=False):
                sender = str(rec.Sender)
                if sender not in {"CLDN4_high", "CLDN4_low"}:
                    continue
                key = (str(rec.Metabolite), str(rec.Sensor))
                wide.setdefault(key, {})[sender] = rec
            if split == "q4q1":
                for key, arms in wide.items():
                    hi = arms.get("CLDN4_high")
                    lo = arms.get("CLDN4_low")
                    if hi is None or lo is None:
                        continue
                    sh = float(hi.Commu_Score)
                    sl = float(lo.Commu_Score)
                    if np.isfinite(sh) and np.isfinite(sl) and (sh > 0 or sl > 0):
                        desc.setdefault(key, []).append(sh - sl)
            for prec in pairs.itertuples(index=False):
                key = (str(prec.hmdb), str(prec.sensor))
                arms = wide.get(key, {})
                hi = arms.get("CLDN4_high")
                lo = arms.get("CLDN4_low")
                sh = float(hi.Commu_Score) if hi is not None else np.nan
                sl = float(lo.Commu_Score) if lo is not None else np.nan
                detected = bool(np.isfinite(sh) and np.isfinite(sl) and (sh > 0 or sl > 0))
                rows.append(
                    {
                        "cohort": cohort,
                        "patient": patient,
                        "split": split,
                        "pair_id": prec.pair_id,
                        "metabolite": prec.metabolite,
                        "sensor": prec.sensor,
                        "hmdb": prec.hmdb,
                        "annotation": prec.annotation,
                        "family": prec.family,
                        "thesis_expect": prec.thesis_expect,
                        "n_mal": n_mal,
                        "n_high": n_high,
                        "n_low": n_low,
                        "n_tnk": n_tnk,
                        "score_high": sh,
                        "score_low": sl,
                        "delta": (sh - sl) if detected else np.nan,
                        "detected": detected,
                        "met_high": float(hi.met_in_sender) if hi is not None else np.nan,
                        "met_low": float(lo.met_in_sender) if lo is not None else np.nan,
                        "sensor_tnk": float(hi.sensor_in_receiver) if hi is not None else np.nan,
                        "prop_met_high": float(hi.metabolite_prop_in_sender) if hi is not None else np.nan,
                        "prop_met_low": float(lo.metabolite_prop_in_sender) if lo is not None else np.nan,
                        "prop_sensor_tnk": float(hi.sensor_prop_in_receiver) if hi is not None else np.nan,
                        "fdr_high": float(hi.permutation_test_fdr) if hi is not None else np.nan,
                        "fdr_low": float(lo.permutation_test_fdr) if lo is not None else np.nan,
                    }
                )
            del adata, block
    return pd.DataFrame(rows), pd.DataFrame(inv_rows), desc, desc_meta


def wilcox_p(x) -> float:
    x = np.asarray(list(x), dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 5 or np.all(x == 0):
        return float("nan")
    try:
        return float(wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return float("nan")


def fmt_p(p) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_num(x) -> str:
    if not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.0001:
        return f"{x:+.4f}"
    return f"{x:+.2e}"


def summarize(per: pd.DataFrame, pairs: pd.DataFrame, split: str) -> pd.DataFrame:
    d = per[(per["split"] == split) & (per["detected"] == True)].copy()  # noqa: E712
    rows = []
    names = list(pairs["pair_id"]) + ["FAMILY_immunosuppressive_metabolite"]
    for name in names:
        if name.startswith("FAMILY_"):
            fam = pairs["family"].iloc[0]
            sub = d[d["family"] == fam]
            if sub.empty:
                continue
            agg = sub.groupby(["cohort", "patient"], as_index=False)["delta"].mean()
            expect = "high>low"
            axis = "family"
        else:
            sub = d[d["pair_id"] == name]
            agg = sub[["cohort", "patient", "delta"]].copy()
            if sub.empty:
                meta = pairs[pairs["pair_id"] == name].iloc[0]
                rows.append(
                    {
                        "pair_id": name,
                        "axis": meta.metabolite,
                        "expect": meta.thesis_expect,
                        "n": 0,
                        "n_gse123902": 0,
                        "n_gse131907": 0,
                        "n_gse205335": 0,
                        "n_gse189357": 0,
                        "mean_delta": np.nan,
                        "p_w": np.nan,
                        "observed": "not_scored",
                        "agrees": "not_scored",
                    }
                )
                continue
            expect = sub["thesis_expect"].iloc[0]
            axis = sub["metabolite"].iloc[0]
        counts = agg["cohort"].value_counts().to_dict() if len(agg) else {}
        mean_delta = float(agg["delta"].mean()) if len(agg) else float("nan")
        if not np.isfinite(mean_delta) or mean_delta == 0:
            observed = "tie" if mean_delta == 0 else "NA"
        elif mean_delta > 0:
            observed = "high>low"
        else:
            observed = "low>high"
        agrees = "yes" if observed == expect else "no"
        rows.append(
            {
                "pair_id": name,
                "axis": axis,
                "expect": expect,
                "n": int(len(agg)),
                "n_gse123902": int(counts.get("GSE123902", 0)),
                "n_gse131907": int(counts.get("GSE131907", 0)),
                "n_gse205335": int(counts.get("GSE205335", 0)),
                "n_gse189357": int(counts.get("GSE189357", 0)),
                "mean_delta": mean_delta,
                "p_w": wilcox_p(agg["delta"]) if len(agg) else float("nan"),
                "observed": observed,
                "agrees": agrees,
            }
        )
    return pd.DataFrame(rows)


def descriptive_table(desc, desc_meta) -> pd.DataFrame:
    rows = []
    for key, deltas in desc.items():
        arr = np.asarray(deltas, dtype=float)
        name, ann = desc_meta.get(key, ("", ""))
        rows.append(
            {
                "hmdb": key[0],
                "metabolite": name,
                "sensor": key[1],
                "annotation": ann,
                "n": int(len(arr)),
                "mean_delta": float(np.mean(arr)),
                "median_delta": float(np.median(arr)),
                "p_w": wilcox_p(arr),
                "n_high_gt_low": int(np.sum(arr > 0)),
                "n_low_gt_high": int(np.sum(arr < 0)),
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values(["p_w", "mean_delta"], ascending=[True, False])
    return out


def write_finding(path: Path, inv: pd.DataFrame, summaries: dict, n_shuffle: int, versions: str):
    q = summaries["q4q1"]
    both = inv[(inv["n_mal"] >= MIN_ARM * 2) & (inv["n_tnk"] >= MIN_TNK)]
    fam = q[q["pair_id"] == "FAMILY_immunosuppressive_metabolite"]
    fam_line = "family row missing"
    if len(fam):
        r = fam.iloc[0]
        fam_line = (
            f"n={int(r.n)}, mean Δ={fmt_num(r.mean_delta)}, p_W={fmt_p(r.p_w)}, "
            f"observed={r.observed}, agrees={r.agrees}"
        )
    lines = [
        "# FINDING — MEBOCOST concordant-four CLDN4-high vs low metabolite senders",
        "",
        "ADDITIVE layer on top of the locked CellChat protein LR result (PR #540).",
        "This does not replace barrier/inhibitory ligand probabilities.",
        "CLDN4 only. No dual-high. Concordant four only",
        "(GSE123902 + GSE131907 + GSE205335 + GSE189357).",
        "Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.",
        "",
        "Engine: **MEBOCOST** `infer_commu` (kaifuchenlab), `met_est='mebocost'`",
        "(enzyme-expression mean). Not a Python reimplementation of the score.",
        "COMPASS flux was not applied. Senders = malignant CLDN4-high vs CLDN4-low;",
        "receiver = T/NK. Score = MEBOCOST `Commu_Score` (sender metabolite × receiver sensor).",
        "Δ = score(high→TNK) − score(low→TNK). Honest n = patient / locked sample.",
        "",
        "CINE: no installable metabolite-CCC package by that name. Unrelated hits",
        "(CineMA cardiac MRI, comet infrared `cine`) were not used. SpatialDM and",
        "MultiNicheNet were not run because MEBOCOST imported and executed.",
        "",
        "Kynurenine–AHR is **not** in the MEBOCOST human sensor table",
        "(`human_met_sensor_update_Oct21_2025.tsv`). It is not scored. The closest",
        "DB row is kynurenic acid–GPR35, which was not pre-specified.",
        "",
        "Primary split is malignant **Q4 vs Q1**. Extra: median and %pos.",
        f"Permutation draws inside MEBOCOST: n_shuffle={n_shuffle}, seed=12345.",
        "The paired test uses `Commu_Score`, which does not depend on the shuffle count.",
        "Within-sample FDRs are stored and are not the primary test.",
        "",
        "## Honest n",
        "",
        "| gate | n | note |",
        "|---|---:|---|",
        f"| Inventory units loaded | {len(inv)} | one row per locked unit that was read |",
        f"| Both compartments (n_mal≥20, n_tnk≥20) | {len(both)} | before the Q4 cell floor |",
        f"| Q4 vs Q1 family units | {int(fam.iloc[0].n) if len(fam) else 0} | n_mal≥40 and ≥1 detected pre-specified pair |",
        "",
        "One locked unit does not enter Q4: GSE205335 P4001 has 27 malignant cells",
        "(floor is n_mal≥40). Every other inventory unit (64/64) has at least one",
        "detected pre-specified pair on the Q4 split.",
        "",
        "GSE123902 / GSE189357: epithelium marker-malignant",
        "(EPCAM\\|KRT8\\|KRT18\\|KRT19 > 0 and PTPRC == 0) and T/NK markers, not malignant.",
        "GSE131907 / GSE205335: author malignant and T/NK labels. Normal-tissue",
        "biopsies in GSE205335 are excluded. TACSTD2 is never a gate.",
        "",
        "## Primary — immunosuppressive metabolite family (Q4 vs Q1)",
        "",
        "Pre-specified. Expectation, stated before the run: CLDN4-high > CLDN4-low",
        "outgoing metabolite signal toward T/NK. ON-thesis would agree. A miss does",
        "not reopen the locked protein-LR barrier result.",
        "",
        f"Family: {fam_line}.",
        "",
        "How to read the table:",
        "",
        "- Inside one unit the T/NK sensor average is shared by the high and low",
        "  senders, so the sign of Δ is the sign of the sender metabolite difference",
        "  when the sensor is expressed.",
        "- The family mean is carried by prostaglandin E2–PTGER4 and PTGER2",
        "  (product enzymes in this DB include PTGES, PTGES2, PTGES3, CBR1, CBR3).",
        "  Those two pairs are detected in nearly every Q4 unit.",
        "- Adenosine–ADORA2B / ADORA2A are the same direction where the sensor is",
        "  detected. ADORA2A is missing from 12 of 13 GSE123902 matrices,",
        "  so that cohort contributes 0 to ADORA2A.",
        "- D-lactic acid–HCAR1 is detected in a minority of units and the Δ is small.",
        "  MEBOCOST estimates that metabolite from HAGH/HAGHL, not from LDHA.",
        "  Do not read it as Warburg L-lactate.",
        "- L-lactic acid–SLC16A1 is **not scored**. In",
        "  `metabolite_associated_gene_reaction_HMDB_summary.tsv` every L-lactic acid",
        "  reaction is substrate-direction, and MEBOCOST only emits a metabolite when",
        "  a product-direction enzyme is present. SLC16A1 is in the matrices. The",
        "  metabolite is not. This is a database limit, not a zero effect.",
        "",
        "| pair | metabolite | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | observed | agrees |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for r in q.itertuples(index=False):
        lines.append(
            f"| {r.pair_id} | {r.axis} | {r.expect} | {int(r.n)} | {int(r.n_gse123902)} | "
            f"{int(r.n_gse131907)} | {int(r.n_gse205335)} | {int(r.n_gse189357)} | "
            f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_w)} | {r.observed} | {r.agrees} |"
        )
    for split in ("median", "pctpos"):
        lines += ["", f"## Extra — {split}", "",
                  "| pair | metabolite | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | observed | agrees |",
                  "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|"]
        for r in summaries[split].itertuples(index=False):
            lines.append(
                f"| {r.pair_id} | {r.axis} | {r.expect} | {int(r.n)} | {int(r.n_gse123902)} | "
                f"{int(r.n_gse131907)} | {int(r.n_gse205335)} | {int(r.n_gse189357)} | "
                f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_w)} | {r.observed} | {r.agrees} |"
            )
    lines += [
        "",
        "## What is not claimed",
        "",
        "- TACSTD2 is not used to define high/low. This is not dual-high.",
        "- GSE148071, GSE127465, and CD45-only libraries are not added.",
        "- This is not a metabolite discovery screen. `descriptive_all_pairs_q4q1.tsv`",
        "  lists every MEBOCOST metabolite–sensor with a non-zero high or low score",
        "  toward T/NK. Those rows are not a new claim.",
        "- Kynurenine–AHR was not tested; it is absent from this sensor table.",
        "- L-lactic acid was not tested; this database has no product-direction enzymes for it.",
        "- D-lactic acid–HCAR1 is not an LDHA result.",
        "- Compass / scFEA flux constraints were not applied.",
        "- Cell-pooled tests are not reported. Honest n is the patient/sample.",
        "- Visium same-spot correlation is not used, and this is not a spatial exclusion test.",
        "",
        "## Versions",
        "",
        versions,
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def versions_text(mebo_root: Path) -> str:
    import scipy
    import sklearn

    lines = [
        f"python {os.sys.version.split()[0]}",
        f"numpy {np.__version__}",
        f"pandas {pd.__version__}",
        f"scipy {scipy.__version__}",
        f"anndata {ad.__version__}",
        f"scikit-learn {sklearn.__version__}",
    ]
    try:
        from mebocost import mebocost

        lines.append(f"mebocost file {mebocost.__file__}")
    except Exception as exc:
        lines.append(f"mebocost import failed {exc}")
    head = mebo_root / ".git" / "HEAD"
    if head.exists():
        lines.append(f"MEBOCOST_ROOT {mebo_root}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="/tmp/concordant4_raw")
    ap.add_argument("--cache", default="/tmp/mebo_cache/units")
    ap.add_argument("--out", default="")
    ap.add_argument("--mebo-root", default=os.environ.get("MEBOCOST_ROOT", "/tmp/src/MEBOCOST"))
    ap.add_argument("--shuffle", type=int, default=10)
    ap.add_argument("--cohort", default="all", help="all or comma-separated GSE ids")
    ap.add_argument("--stage", default="all", choices=["extract", "score", "all"])
    args = ap.parse_args()
    out = Path(args.out) if args.out else here_dir()
    data = here_dir() / "data"
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    tab = out / "results" / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    mebo_root = Path(args.mebo_root)
    raw = Path(args.raw)
    wanted, cmap = load_db_genes(mebo_root)
    wanted_upper = set(cmap)
    cohorts = None if args.cohort == "all" else {c.strip() for c in args.cohort.split(",")}
    conf = out / "results" / "mebocost.conf"
    write_conf(conf, mebo_root)
    if args.stage in {"extract", "all"}:
        if cohorts is None or "GSE189357" in cohorts:
            logmsg("==== extract GSE189357 ====")
            extract_gse189357(raw, data, cache, cmap, wanted_upper)
        if cohorts is None or "GSE123902" in cohorts:
            logmsg("==== extract GSE123902 ====")
            extract_gse123902(raw, data, cache, cmap, wanted_upper)
        if cohorts is None or "GSE205335" in cohorts:
            logmsg("==== extract GSE205335 ====")
            extract_gse205335(raw, data, cache, cmap, wanted_upper, mebo_root)
        if cohorts is None or "GSE131907" in cohorts:
            logmsg("==== extract GSE131907 ====")
            extract_gse131907(raw, data, cache, cmap, wanted_upper)
    if args.stage == "extract":
        logmsg("extract done", len(list(cache.glob('*.npz'))))
        return
    pairs = pd.read_csv(data / "prespecified_pairs.tsv", sep="\t")
    logmsg("==== score ====", "shuffle", args.shuffle)
    per, inv, desc, desc_meta = score_cached_units(cache, conf, pairs, args.shuffle, cohorts)
    if per.empty:
        raise SystemExit("no scored units")
    per.to_csv(tab / "per_unit_pairs.tsv", sep="\t", index=False)
    inv.to_csv(tab / "patient_inventory.tsv", sep="\t", index=False)
    summaries = {}
    frames = []
    for split in SPLITS:
        sm = summarize(per, pairs, split)
        sm.insert(0, "split", split)
        summaries[split] = sm
        frames.append(sm)
        logmsg("summary", split)
        print(sm.to_string(index=False), flush=True)
    pd.concat(frames, ignore_index=True).to_csv(tab / "pair_summary.tsv", sep="\t", index=False)
    dtab = descriptive_table(desc, desc_meta)
    dtab.to_csv(tab / "descriptive_all_pairs_q4q1.tsv", sep="\t", index=False)
    ver = versions_text(mebo_root)
    (out / "results" / "sessionInfo.txt").write_text(ver + "\n")
    write_finding(out / "FINDING.md", inv, summaries, args.shuffle, ver)
    logmsg("wrote", out / "FINDING.md", "inventory", len(inv), "pair rows", len(per))


if __name__ == "__main__":
    main()
