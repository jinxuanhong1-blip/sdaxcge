#!/usr/bin/env python3
"""Stream lineage + LR genes from GSE131907, GSE148071, GSE205335.

Keeps malignant / malignant-like and T/NK cells only. Writes a compact
extract used by 03_run_ccc.py. TACSTD2 is not used as a gate.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import re
import shutil
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

HERE = Path(__file__).resolve().parents[1]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}


def log(msg: str) -> None:
    print(msg, flush=True)


def wanted_from_cfg(cfg: dict, extra: list[str]) -> set[str]:
    genes = set(extra)
    for vs in cfg["lineage_markers"].values():
        genes.update(vs)
    genes.update(cfg["normal_lung"])
    genes.update(cfg["state_genes"])
    return genes


def score_lineage(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], lineages: dict, n: int) -> np.ndarray:
    names = list(lineages)
    scores = np.vstack([score_lineage(expr, lineages[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def stream_dense(path: Path, keep: set[str], header_has_index: bool) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, int]:
    log(f"[stream] {path.name} keep={len(keep)}")
    t0 = time.time()
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header if h != ""]
        if header_has_index or (header and header[0] in {"", "Index", "index", "gene", "Gene", "GENE"}):
            cells = np.array(header[1:], dtype=object)
        else:
            cells = np.array(header, dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.strip().strip('"').split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                toks = line.rstrip("\n").split("\t")
                vals = np.asarray(toks[1:], dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"{path.name} {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"  {path.name}: scanned {n_genes} genes, kept {len(store)} [{time.time()-t0:.0f}s]")
    log(f"[stream] {path.name}: cells={n} genes={n_genes} kept={len(store)} [{time.time()-t0:.0f}s]")
    return cells, store, totals, n_genes


def parse_gse131907_series(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_gse205335_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_rds_genes(path: Path, wanted: set[str]):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            log(f"decompress {path.name}")
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        log("read GSE205335 RDS")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
            extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    del matrix, obj
    gc.collect()
    log(f"extracted {len(extracted)} / {len(wanted)} GSE205335 genes")
    return barcodes, extracted, library_umi, genes.tolist()


def pack_kept(
    dataset: str,
    cell_ids: np.ndarray,
    patients: np.ndarray,
    samples: np.ndarray,
    origins: np.ndarray,
    label_source: str,
    is_malig: np.ndarray,
    is_t: np.ndarray,
    is_nk: np.ndarray,
    totals: np.ndarray,
    store: dict[str, np.ndarray],
    gene_order: list[str],
) -> tuple[pd.DataFrame, np.ndarray]:
    keep = is_malig | is_t | is_nk
    idx = np.flatnonzero(keep)
    lib = np.maximum(totals[idx], 1.0)
    cldn = store["CLDN4"][idx] if "CLDN4" in store else np.zeros(idx.size, dtype=np.float32)
    cells = pd.DataFrame(
        {
            "cell_id": [f"{dataset}|{c}" for c in cell_ids[idx]],
            "dataset": dataset,
            "patient": patients[idx],
            "sample": samples[idx],
            "origin": origins[idx],
            "label_source": label_source,
            "is_malig": is_malig[idx],
            "is_t": is_t[idx],
            "is_nk": is_nk[idx],
            "is_tnk": is_t[idx] | is_nk[idx],
            "total_umi": totals[idx],
            "cldn4_log1p": np.log1p(cldn / lib * 1e4).astype(np.float32),
        }
    )
    mat = np.zeros((idx.size, len(gene_order)), dtype=np.float32)
    for j, g in enumerate(gene_order):
        if g in store:
            mat[:, j] = np.log1p(store[g][idx] / lib * 1e4)
    return cells, mat


def extract_gse131907(datadir: Path, keep: set[str], gene_order: list[str]) -> tuple[pd.DataFrame, np.ndarray, dict]:
    umi = datadir / "GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann_path = datadir / "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    series_path = datadir / "GSE131907/GSE131907_series_matrix.txt.gz"
    cells_ids, store, totals, n_genes = stream_dense(umi, keep, header_has_index=True)
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    id_col = "Index" if "Index" in ann.columns else ann.columns[0]
    ann = ann.set_index(id_col).reindex(cells_ids)
    if ann["Sample"].isna().any():
        raise SystemExit(f"GSE131907 annotation miss {int(ann['Sample'].isna().sum())}")
    sample_col = "Sample"
    origin_col = "Sample_Origin" if "Sample_Origin" in ann.columns else "Sample_Origion"
    type_col = "Cell_type"
    subtype_col = "Cell_subtype"
    series = parse_gse131907_series(series_path)
    title_to_patient = {}
    if "patient_id" in series.columns:
        for title, pid in zip(series["title"], series["patient_id"]):
            title_to_patient[str(title)] = str(pid)
    samples = ann[sample_col].to_numpy()
    origins = ann[origin_col].fillna("").to_numpy()
    patients = np.array(
        [f"GSE131907|{title_to_patient.get(s, s)}" for s in samples], dtype=object
    )
    tumor = np.isin(origins, list(TUMOR_ORIGINS))
    is_malig = (
        tumor
        & (ann[type_col].to_numpy() == "Epithelial cells")
        & np.isin(ann[subtype_col].fillna("").to_numpy(), list(MALIGNANT_SUBTYPES))
    )
    is_t = tumor & (ann[type_col].to_numpy() == "T lymphocytes")
    is_nk = tumor & (ann[type_col].to_numpy() == "NK cells")
    cells, mat = pack_kept(
        "GSE131907",
        cells_ids,
        patients,
        samples,
        origins,
        "author_Cell_type",
        is_malig,
        is_t,
        is_nk,
        totals,
        store,
        gene_order,
    )
    inv = {
        "n_cells_matrix": int(len(cells_ids)),
        "n_genes_matrix": int(n_genes),
        "n_tumor_cells": int(tumor.sum()),
        "n_malig_kept": int(is_malig.sum()),
        "n_t_kept": int(is_t.sum()),
        "n_nk_kept": int(is_nk.sum()),
        "n_patients_deposited": int(pd.Series(patients).nunique()),
        "genes_found": sorted(store),
    }
    return cells, mat, inv


def extract_gse148071(datadir: Path, keep: set[str], gene_order: list[str], cfg: dict) -> tuple[pd.DataFrame, np.ndarray, dict]:
    files = sorted((datadir / "GSE148071/GSE148071_files").rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit("GSE148071 exp matrices missing")
    parts_c, parts_m = [], []
    n_cells = 0
    n_mal = n_t = n_nk = 0
    for path in files:
        m = re.search(r"_(P\d+)_", path.name)
        pid = m.group(1) if m else path.name
        cell_ids, store, totals, _n = stream_dense(path, keep, header_has_index=False)
        n = len(cell_ids)
        n_cells += n
        lineage = assign_lineage(store, cfg["lineage_markers"], n)
        normal = np.zeros(n, dtype=bool)
        for g in cfg["normal_lung"]:
            if g in store:
                normal |= store[g] > 0
        is_malig = (lineage == "epithelial") & (~normal)
        is_t = lineage == "T"
        is_nk = lineage == "NK"
        n_mal += int(is_malig.sum())
        n_t += int(is_t.sum())
        n_nk += int(is_nk.sum())
        patients = np.array([f"GSE148071|{pid}"] * n, dtype=object)
        samples = np.array([pid] * n, dtype=object)
        origins = np.array(["biopsy"] * n, dtype=object)
        cells, mat = pack_kept(
            "GSE148071",
            cell_ids,
            patients,
            samples,
            origins,
            "A3_marker_argmax",
            is_malig,
            is_t,
            is_nk,
            totals,
            store,
            gene_order,
        )
        parts_c.append(cells)
        parts_m.append(mat)
        log(f"  {pid}: cells={n} malig={int(is_malig.sum())} T={int(is_t.sum())} NK={int(is_nk.sum())}")
    cells = pd.concat(parts_c, ignore_index=True)
    mat = np.vstack(parts_m) if parts_m else np.zeros((0, len(gene_order)), dtype=np.float32)
    inv = {
        "n_cells_matrix": int(n_cells),
        "n_files": int(len(files)),
        "n_malig_kept": n_mal,
        "n_t_kept": n_t,
        "n_nk_kept": n_nk,
        "n_patients_deposited": int(cells["patient"].nunique()) if len(cells) else 0,
    }
    return cells, mat, inv


def extract_gse205335(datadir: Path, keep: set[str], gene_order: list[str]) -> tuple[pd.DataFrame, np.ndarray, dict]:
    barcodes, store, totals, all_genes = load_rds_genes(
        datadir / "GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz", keep
    )
    ident = pd.read_csv(
        datadir / "GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t"
    )
    if "barcode" not in ident.columns:
        ident = ident.rename(columns={ident.columns[0]: "barcode"})
    indexed = ident.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    if len(missing):
        raise SystemExit(f"GSE205335 identity miss {len(missing)}")
    cells_ann = indexed.loc[barcodes]
    meta = parse_gse205335_soft(datadir / "GSE205335/GSE205335_family.soft.gz")
    if "orig.ident" in cells_ann.columns:
        cells_ann = cells_ann.merge(
            meta[["orig.ident", "patient"]].drop_duplicates(),
            on="orig.ident",
            how="left",
        )
    elif "patient" not in cells_ann.columns:
        raise SystemExit("GSE205335 identity table has no orig.ident/patient")
    if cells_ann["patient"].isna().any():
        nmiss = int(cells_ann["patient"].isna().sum())
        log(f"WARN GSE205335 {nmiss} cells without GEO patient; dropping")
        ok = cells_ann["patient"].notna().to_numpy()
        barcodes = barcodes[ok]
        totals = totals[ok]
        store = {g: v[ok] for g, v in store.items()}
        cells_ann = cells_ann.loc[ok]
    patients = np.array([f"GSE205335|{p}" for p in cells_ann["patient"].astype(str)], dtype=object)
    samples = (
        cells_ann["orig.ident"].astype(str).to_numpy()
        if "orig.ident" in cells_ann.columns
        else cells_ann["patient"].astype(str).to_numpy()
    )
    origins = (
        cells_ann["tissue"].astype(str).to_numpy()
        if "tissue" in cells_ann.columns
        else np.array(["tumor"] * len(barcodes), dtype=object)
    )
    is_malig = cells_ann["lineage.sub"].astype(str).to_numpy() == "Malignant cells"
    is_tnk = cells_ann["lineage.total"].astype(str).to_numpy() == "T/NK cells"
    # Split T vs NK when a finer label exists; otherwise both flags stay in T/NK.
    if "lineage.sub" in cells_ann.columns:
        sub = cells_ann["lineage.sub"].astype(str).to_numpy()
        is_t = is_tnk & ~np.isin(sub, ["NK cells", "NK"])
        is_nk = is_tnk & np.isin(sub, ["NK cells", "NK"])
        # If author T/NK is not split, keep all as T so is_tnk still works.
        if is_nk.sum() == 0:
            is_t = is_tnk
    else:
        is_t = is_tnk
        is_nk = np.zeros(len(barcodes), dtype=bool)
    cells, mat = pack_kept(
        "GSE205335",
        barcodes,
        patients,
        samples,
        origins,
        "author_lineage",
        is_malig,
        is_t,
        is_nk,
        totals,
        store,
        gene_order,
    )
    inv = {
        "n_cells_matrix": int(len(barcodes)),
        "n_genes_matrix": int(len(all_genes)),
        "n_malig_kept": int(is_malig.sum()),
        "n_t_kept": int(is_t.sum()),
        "n_nk_kept": int(is_nk.sum()),
        "n_patients_deposited": int(pd.Series(patients).nunique()),
        "genes_found": sorted(store),
    }
    return cells, mat, inv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    ap.add_argument("--outdir", type=Path, default=HERE / "data" / "extract")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    lr_genes = (HERE / "resources" / "lr_genes.txt").read_text().split()
    keep = wanted_from_cfg(cfg, lr_genes)
    gene_order = sorted(keep)
    log(f"wanted genes={len(gene_order)}")

    frames, mats, inventory = [], [], {}
    c, m, inv = extract_gse131907(args.datadir, keep, gene_order)
    frames.append(c)
    mats.append(m)
    inventory["GSE131907"] = inv
    log(f"GSE131907 kept {len(c)}")
    gc.collect()

    c, m, inv = extract_gse148071(args.datadir, keep, gene_order, cfg)
    frames.append(c)
    mats.append(m)
    inventory["GSE148071"] = inv
    log(f"GSE148071 kept {len(c)}")
    gc.collect()

    c, m, inv = extract_gse205335(args.datadir, keep, gene_order)
    frames.append(c)
    mats.append(m)
    inventory["GSE205335"] = inv
    log(f"GSE205335 kept {len(c)}")
    gc.collect()

    cells = pd.concat(frames, ignore_index=True)
    mat = np.vstack(mats)
    cells.to_parquet(args.outdir / "cells.parquet", index=False)
    np.savez_compressed(
        args.outdir / "log1p_cp10k.npz",
        log1p=mat,
        genes=np.array(gene_order, dtype=object),
        cell_id=cells["cell_id"].to_numpy(),
    )
    inventory["n_cells_kept"] = int(len(cells))
    inventory["n_genes"] = int(len(gene_order))
    inventory["n_malig"] = int(cells["is_malig"].sum())
    inventory["n_tnk"] = int(cells["is_tnk"].sum())
    (args.outdir / "inventory.json").write_text(json.dumps(inventory, indent=2))
    log(json.dumps({k: inventory[k] for k in ("n_cells_kept", "n_malig", "n_tnk")}, indent=2))


if __name__ == "__main__":
    main()
