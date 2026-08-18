#!/usr/bin/env python3
"""ADDITIVE CLDN4-only CellChat-style barrier call on the concordant four.

GSE123902 + GSE131907 + GSE205335 + GSE189357. Not a full-pool.
Do NOT add GSE148071 / GSE127465 / CD45-only.

Thesis already correct. Ligands stay:
  CLDN4-high → more barrier/inhibitory outgoing to T/NK
    (F11R, NECTIN2–TIGIT, CDH1, LGALS9)  — ON-thesis
  CLDN4-low / KD-like → more IFN / T-recruit outgoing
    (often NOT detected; CXCL9/10 sparse)

CLDN4-only gates (Q4 vs Q1 or %pos). No dual-high.
Primary readout: barrier/inhibitory family ΔP (high − low).
Honest n = units with both compartments.
PR #424 / #474 numbers are comparison rows only — not re-audited.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib import (  # noqa: E402
    COMPARISON_ROWS,
    EPI_MARKERS,
    EXPR_PROP,
    FAMILIES,
    GENE_ALIASES,
    KH,
    MALIG_SUBTYPES,
    MIN_CELLS_ARM,
    MIN_CELLS_COMP,
    MIN_MAL_Q4,
    TNK_MARKERS,
    TRIM,
    TUMOR_ORIGINS,
    apply_aliases,
    dersimonian_laird,
    fmt_num,
    fmt_p,
    highend_split,
    load_lr,
    score_outgoing,
    thesis_genes,
    unit_delta_stats,
)

LOCKED = ROOT / "data"
DB = ROOT / "db"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLORS = {
    "GSE123902": "#4C78A8",
    "GSE131907": "#F58518",
    "GSE205335": "#54A24B",
    "GSE189357": "#E45756",
}
N_KEYS = {
    "GSE123902": "n_gse123902",
    "GSE131907": "n_gse131907",
    "GSE205335": "n_gse205335",
    "GSE189357": "n_gse189357",
}

# PR #459 singles for the concordant four (T/NK Spearman; not re-audited).
GIVEN_SINGLES = [
    {"cohort": "GSE123902", "n": 13, "rho": -0.659, "p": 0.014, "score": "pct"},
    {"cohort": "GSE131907", "n": 21, "rho": -0.522, "p": 0.015, "score": "pct"},
    {"cohort": "GSE205335", "n": 22, "rho": -0.435, "p": 0.043, "score": "pct"},
    {"cohort": "GSE189357", "n": 9, "rho": -0.600, "p": 0.088, "score": "pct"},
]
GIVEN_LOCKED = {"GSE123902": 13, "GSE131907": 21, "GSE205335": 22, "GSE189357": 9}
PRIMARY_SPLIT = "q4q1"
SPLITS = ("q4q1", "pctpos", "median", "tertile")


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    extracted = apply_aliases(extracted)
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def load_locked() -> dict[str, pd.DataFrame]:
    out = {}
    d = pd.read_csv(LOCKED / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    el["unit"] = "GSE123902"
    el["patient_id"] = el["patient"].astype(str)
    out["GSE123902"] = el

    s = pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t")
    el2 = s[(s["n_malignant"] > 0) & (s["n_tnk"] > 0)].copy()
    el2["unit"] = "GSE131907"
    el2["patient_id"] = el2["sample"].astype(str)
    out["GSE131907"] = el2

    h = pd.read_csv(LOCKED / "GSE205335_patients.tsv", sep="\t")
    el3 = h[(h["n_malignant"] > 0) & (h["n_tnk"] > 0)].copy()
    el3["unit"] = "GSE205335"
    el3["patient_id"] = el3["patient"].astype(str)
    out["GSE205335"] = el3

    z = pd.read_csv(LOCKED / "GSE189357_marker_units.tsv", sep="\t")
    el4 = z[z["eligible"].astype(str).str.lower() == "true"].copy()
    el4["unit"] = "GSE189357"
    el4["patient_id"] = el4["patient"].astype(str)
    out["GSE189357"] = el4
    return out


def _marker_mal_tnk(extracted: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def col(g: str) -> np.ndarray:
        if g not in extracted:
            return np.zeros(n, dtype=float)
        return extracted[g]

    epi = np.zeros(n, dtype=bool)
    for g in EPI_MARKERS:
        epi |= col(g) > 0
    mal = epi & (col("PTPRC") == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in TNK_MARKERS:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    return mal, tnk


def _score_one_123902_csv(handle, wanted: set[str]):
    df = pd.read_csv(handle, index_col=0)
    df.columns = [str(c).upper() for c in df.columns]
    gene_set = set(df.columns)
    libs = df.sum(axis=1).to_numpy(dtype=np.float64)
    extracted = {g: df[g].to_numpy(dtype=np.float32) for g in wanted if g in df.columns}
    return extracted, libs, gene_set


def load_gse123902(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE123902 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    tar_path = raw / "GSE123902" / "GSE123902_RAW.tar"
    chosen: dict[str, str] = {}
    members = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            tissue = (
                "NORMAL"
                if "NORMAL" in name
                else (
                    "METASTASIS"
                    if "METASTASIS" in name
                    else ("PRIMARY" if "PRIMARY" in name else "OTHER")
                )
            )
            if tissue not in {"PRIMARY", "METASTASIS"}:
                continue
            if donor not in keep:
                continue
            members.append((donor, tissue, m, name))
        members.sort(key=lambda x: (x[0], 0 if x[1] == "PRIMARY" else 1))
        chunks = []
        gene_union: set[str] = set()
        for donor, tissue, m, name in members:
            if donor in chosen:
                continue
            chosen[donor] = tissue
            raw_f = gzip.GzipFile(fileobj=tf.extractfile(m))
            extracted, libs, genes = _score_one_123902_csv(raw_f, wanted)
            gene_union |= genes
            n = int(libs.size)
            chunks.append((donor, extracted, libs))
            print(f"  {name}: cells={n} stored={len(extracted)} donor={donor} {tissue}", flush=True)
    if not chunks:
        raise SystemExit("GSE123902: no locked tumor matrices read")
    units = []
    for donor, extracted, libs in chunks:
        n = int(libs.size)
        log_cp, pos = to_log_pos(extracted, libs)
        mal, tnk = _marker_mal_tnk(apply_aliases(extracted), n)
        cldn4 = log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))
        units.append(
            {
                "cohort": "GSE123902",
                "patient_id": donor,
                "n_cells": n,
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4,
                "log_cp": log_cp,
                "pos": pos,
                "genes": set(apply_aliases(extracted)),
                "malig_def": "marker_malignant (epithelium proxy; author labels absent)",
            }
        )
        print(f"  unit {donor}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())}", flush=True)
    return units, gene_union


def stream_131907(matrix_path: Path, wanted: set[str]):
    cache = Path("/tmp/concordant4_raw/GSE131907/wanted_stream_cache.npz")
    if cache.exists():
        print(f"  load stream cache {cache}", flush=True)
        z = np.load(cache, allow_pickle=True)
        cell_ids = z["cell_ids"].tolist()
        n_umi = z["n_umi"]
        n_streamed = int(z["n_streamed"])
        found = {g: z[f"g_{g}"] for g in z["genes"].tolist() if f"g_{g}" in z.files}
        print(f"  cache cells={len(cell_ids)} stored={len(found)} genes={n_streamed}", flush=True)
        return cell_ids, apply_aliases(found), n_umi, n_streamed
    found: dict[str, np.ndarray] = {}
    alias_wanted = set(wanted) | set(GENE_ALIASES_KEYS())
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in alias_wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    found = apply_aliases(found)
    payload = {
        "cell_ids": np.asarray(cell_ids, dtype=object),
        "n_umi": n_umi,
        "n_streamed": n_streamed,
        "genes": np.asarray(sorted(found), dtype=object),
    }
    for g, arr in found.items():
        payload[f"g_{g}"] = arr
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **payload)
    print(f"  wrote stream cache {cache}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def GENE_ALIASES_KEYS() -> set[str]:
    return set(GENE_ALIASES.keys()) | set(GENE_ALIASES.values())


def load_gse131907(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE131907 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    ann_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t")
    cell_ids, found, n_umi, n_genes = stream_131907(mat_path, wanted)
    cell_index = {c: i for i, c in enumerate(cell_ids)}
    key = "Index" if "Index" in ann.columns else ann.columns[0]
    ann = ann.copy()
    ann["cell_i"] = ann[key].map(cell_index)
    ann = ann[ann["cell_i"].notna()].copy()
    ann["cell_i"] = ann["cell_i"].astype(int)
    sample_col = "Sample" if "Sample" in ann.columns else "sample"
    origin_col = "Sample_Origin" if "Sample_Origin" in ann.columns else None
    type_col = "Cell_type" if "Cell_type" in ann.columns else None
    sub_col = "Cell_subtype" if "Cell_subtype" in ann.columns else None

    mal_mask = np.zeros(len(cell_ids), dtype=bool)
    tnk_mask = np.zeros(len(cell_ids), dtype=bool)
    sample_of = np.array([""] * len(cell_ids), dtype=object)
    for rec in ann.itertuples(index=False):
        i = int(rec.cell_i)
        sample_of[i] = str(getattr(rec, sample_col))
        subtype = str(getattr(rec, sub_col)) if sub_col else ""
        ctype = str(getattr(rec, type_col)) if type_col else ""
        origin = str(getattr(rec, origin_col)) if origin_col else ""
        if subtype in MALIG_SUBTYPES and (origin in TUMOR_ORIGINS or origin_col is None):
            mal_mask[i] = True
        if ctype in {"T lymphocytes", "NK cells"}:
            tnk_mask[i] = True

    cldn4_all_log, _ = to_log_pos(found, n_umi)
    cldn4 = cldn4_all_log.get("CLDN4", np.zeros(len(cell_ids), dtype=np.float32))
    units = []
    for sample in sorted(keep):
        idx = np.flatnonzero(sample_of == sample)
        if idx.size == 0:
            print(f"  MISSING sample {sample}", flush=True)
            continue
        mal = mal_mask[idx]
        tnk = tnk_mask[idx]
        extracted = {g: found[g][idx] for g in found}
        libs = n_umi[idx]
        log_u, pos_u = to_log_pos(extracted, libs)
        units.append(
            {
                "cohort": "GSE131907",
                "patient_id": sample,
                "n_cells": int(idx.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4[idx],
                "log_cp": log_u,
                "pos": pos_u,
                "genes": set(found),
                "malig_def": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
            }
        )
        print(
            f"  unit {sample}: n={int(idx.size)} mal={int(mal.sum())} tnk={int(tnk.sum())}",
            flush=True,
        )
    return units, set(found), n_genes


def _parse_205335_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records, current, descriptions, titles = [], None, [], []
    with opener(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
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
    return pd.DataFrame(records)


def load_gse205335(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE205335 ====", flush=True)
    import rdata

    keep = set(locked["patient_id"].astype(str))
    ident = pd.read_csv(raw / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    rds = raw / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        unzipped = Path(tmp) / "matrix.rds"
        print("  decompress RDS", flush=True)
        with gzip.open(rds, "rb") as src, unzipped.open("wb") as dest:
            shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        print("  read RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(unzipped)
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"  build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False).tocsr()
    library = np.asarray(matrix.sum(axis=0)).ravel()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(set(wanted) | GENE_ALIASES_KEYS()):
        row = name_to_row.get(gene)
        if row is not None:
            extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    del matrix, obj
    extracted = apply_aliases(extracted)
    print(f"  extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    indexed = ident.set_index("barcode")
    cells = indexed.loc[barcodes].reset_index()
    soft = _parse_205335_soft(raw / "GSE205335" / "GSE205335_family.soft.gz")
    if "platform" in soft.columns:
        read_end = soft["platform"].astype(str).str.extract(r"Single Cell ([35])'")[0]
        soft["orig.ident"] = (
            soft["description"].astype(str).str.replace("_", "-", regex=False) + "-" + read_end + "P"
        )
    if "patient" not in cells.columns and "orig.ident" in cells.columns and "patient" in soft.columns:
        cells = cells.merge(soft[["orig.ident", "patient"]], on="orig.ident", how="left")
    if "patient" not in cells.columns:
        raise SystemExit(f"GSE205335: cannot map patient IDs; columns={list(cells.columns)}")
    if "lineage.sub" not in cells.columns:
        raise SystemExit("GSE205335 identity table missing lineage.sub")
    mal_all = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    tnk_all = (
        cells["lineage.total"].eq("T/NK cells").to_numpy()
        if "lineage.total" in cells.columns
        else np.zeros(len(cells), dtype=bool)
    )
    patient = cells["patient"].astype(str).to_numpy()
    log_cp_all, _ = to_log_pos(extracted, library)
    cldn4_all = log_cp_all.get("CLDN4", np.zeros(len(barcodes), dtype=np.float32))
    units = []
    for pid in sorted(keep):
        idx = np.flatnonzero(patient == pid)
        if idx.size == 0:
            print(f"  MISSING patient {pid}", flush=True)
            continue
        mal = mal_all[idx]
        tnk = tnk_all[idx]
        extracted_u = {g: extracted[g][idx] for g in extracted}
        libs = library[idx]
        log_u, pos_u = to_log_pos(extracted_u, libs)
        units.append(
            {
                "cohort": "GSE205335",
                "patient_id": pid,
                "n_cells": int(idx.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4_all[idx],
                "log_cp": log_u,
                "pos": pos_u,
                "genes": set(extracted),
                "malig_def": "author lineage.sub == Malignant cells",
            }
        )
        print(
            f"  unit {pid}: n={int(idx.size)} mal={int(mal.sum())} tnk={int(tnk.sum())}",
            flush=True,
        )
    return units, set(genes.tolist())


def _member(members: dict[str, tarfile.TarInfo], sample: str, kind: str) -> str:
    for name in members:
        if f"_{sample}_{kind}" in Path(name).name:
            return name
    raise KeyError(f"no {kind} for {sample}")


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = _member(members, sample, "features")
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for line in handle:
            parts = line.decode().strip().split("\t")
            genes.append((parts[1] if len(parts) > 1 else parts[0]).upper())
    return genes


def _read_10x_nbarcodes(tf: tarfile.TarFile, members: dict, sample: str) -> int:
    name = _member(members, sample, "barcodes")
    n = 0
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for _ in handle:
            n += 1
    return n


def _stream_10x_wanted(tf, members, sample, genes, wanted) -> tuple[dict[str, np.ndarray], np.ndarray]:
    n_cells = _read_10x_nbarcodes(tf, members, sample)
    keep_set = set(wanted) | GENE_ALIASES_KEYS()
    keep_idx = {i for i, g in enumerate(genes) if g in keep_set}
    extracted = {genes[i]: np.zeros(n_cells, dtype=np.float32) for i in keep_idx}
    library = np.zeros(n_cells, dtype=np.float64)
    mtx_name = _member(members, sample, "matrix.mtx")
    with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as handle:
        for line in handle:
            if line.startswith(b"%"):
                continue
            break
        for line in handle:
            a, b, v = line.decode().split()
            gi = int(a) - 1
            ci = int(b) - 1
            val = float(v)
            library[ci] += val
            if gi in keep_idx:
                extracted[genes[gi]][ci] = val
    return apply_aliases(extracted), library


def load_gse189357(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE189357 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    tar_path = raw / "GSE189357" / "GSE189357_RAW.tar"
    units = []
    gene_union: set[str] = set()
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            if sample not in keep:
                print(f"  skip {sample} (not in locked extract)", flush=True)
                continue
            genes = _read_10x_features(tf, members, sample)
            gene_union |= set(genes)
            extracted, libs = _stream_10x_wanted(tf, members, sample, genes, wanted)
            n = int(libs.size)
            log_cp, pos = to_log_pos(extracted, libs)
            mal, tnk = _marker_mal_tnk(extracted, n)
            cldn4 = log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))
            units.append(
                {
                    "cohort": "GSE189357",
                    "patient_id": sample,
                    "n_cells": n,
                    "n_mal": int(mal.sum()),
                    "n_tnk": int(tnk.sum()),
                    "mal": mal,
                    "tnk": tnk,
                    "cldn4": cldn4,
                    "log_cp": log_cp,
                    "pos": pos,
                    "genes": set(extracted),
                    "malig_def": "marker_malignant (epithelium proxy; author labels thin)",
                }
            )
            print(
                f"  {sample}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())} stored={len(extracted)}",
                flush=True,
            )
    if not units:
        raise SystemExit("GSE189357: no locked tumor matrices read")
    return units, gene_union


def score_unit(unit: dict, lr: pd.DataFrame, mode: str) -> dict | None:
    mal_idx = np.flatnonzero(unit["mal"])
    tnk_idx = np.flatnonzero(unit["tnk"])
    # Honest n: both compartments present (malignant high/low + T/NK).
    if tnk_idx.size < MIN_CELLS_COMP:
        return None
    split = highend_split(unit["cldn4"], mal_idx, mode)
    if split is None:
        return None
    hi, lo = split
    hi_rows = score_outgoing(lr, unit["log_cp"], unit["pos"], hi, tnk_idx)
    lo_rows = score_outgoing(lr, unit["log_cp"], unit["pos"], lo, tnk_idx)
    if not hi_rows or not lo_rows:
        return None
    hi_map = {r["interaction_name"]: r for r in hi_rows}
    lo_map = {r["interaction_name"]: r for r in lo_rows}
    pair_rows = []
    for name in hi_map:
        h = hi_map[name]
        l = lo_map[name]
        detected = bool(h["detected"] or l["detected"])
        pair_rows.append(
            {
                "cohort": unit["cohort"],
                "patient_id": unit["patient_id"],
                "split": mode,
                "interaction_name": name,
                "ligand": h["ligand"],
                "receptor": h["receptor"],
                "prob_high": h["prob"],
                "prob_low": l["prob"],
                "delta": h["prob"] - l["prob"],
                "detected": detected,
                "detected_high": bool(h["detected"]),
                "detected_low": bool(l["detected"]),
                "ligand_prop_high": h["ligand_prop"],
                "ligand_prop_low": l["ligand_prop"],
                "receptor_prop": h["receptor_prop"],
                "n_mal_high": int(hi.size),
                "n_mal_low": int(lo.size),
                "n_tnk": int(tnk_idx.size),
                "n_mal": int(mal_idx.size),
            }
        )
    return {
        "cohort": unit["cohort"],
        "patient_id": unit["patient_id"],
        "split": mode,
        "n_mal": int(mal_idx.size),
        "n_mal_high": int(hi.size),
        "n_mal_low": int(lo.size),
        "n_tnk": int(tnk_idx.size),
        "pairs": pair_rows,
    }


def _empty_meta() -> dict:
    return {
        "n": 0,
        "n_123902": 0,
        "n_131907": 0,
        "n_205335": 0,
        "n_189357": 0,
        "mean": np.nan,
        "median": np.nan,
        "se": np.nan,
        "wilcoxon_p": np.nan,
        "meta_p": np.nan,
        "i2": np.nan,
        "n_pos": 0,
        "n_neg": 0,
    }


def meta_pair(df: pd.DataFrame) -> dict:
    """Patient-level ΔP meta. Honest n = units with a finite detected delta."""
    sub = df[df["detected"]].copy() if "detected" in df.columns else df.copy()
    if sub.empty:
        return _empty_meta()
    stats = unit_delta_stats(sub["delta"])
    if stats["n"] >= 2 and np.isfinite(stats["sd"]) and stats["sd"] > 0:
        v = np.full(stats["n"], stats["sd"] ** 2)
        dl = dersimonian_laird(sub["delta"].to_numpy(float), v)
        meta_p = dl["p"]
        i2 = dl["i2"]
        mean = dl["mean"]
        se = dl["se"]
    else:
        meta_p, i2, mean, se = stats["wilcoxon_p"], np.nan, stats["mean"], stats["se"]
    return {
        "n": stats["n"],
        "n_123902": int((sub["cohort"] == "GSE123902").sum()),
        "n_131907": int((sub["cohort"] == "GSE131907").sum()),
        "n_205335": int((sub["cohort"] == "GSE205335").sum()),
        "n_189357": int((sub["cohort"] == "GSE189357").sum()),
        "mean": mean,
        "median": stats["median"],
        "se": se,
        "wilcoxon_p": stats["wilcoxon_p"],
        "meta_p": meta_p,
        "i2": i2,
        "n_pos": int((sub["delta"] > 0).sum()),
        "n_neg": int((sub["delta"] < 0).sum()),
    }


def agrees(expect: str, mean_delta: float, n: int, p: float) -> str:
    if n < 3:
        return "thin"
    if not np.isfinite(mean_delta):
        return "undetected"
    if expect == "high>low":
        direction_ok = mean_delta > 0
    else:
        direction_ok = mean_delta < 0
    if direction_ok and np.isfinite(p) and p < 0.05:
        return "yes"
    if direction_ok:
        return "direction_only"
    if np.isfinite(p) and p < 0.05:
        return "opposite"
    return "no"


def family_table(long_df: pd.DataFrame, family_key: str, mode: str) -> pd.DataFrame:
    spec = FAMILIES[family_key]
    rows = []
    pair_long = []
    for rec in spec["pairs"]:
        name = rec["interaction_name"]
        sub = long_df[(long_df["interaction_name"] == name) & (long_df["split"] == mode)]
        m = meta_pair(sub)
        observed = (
            "high>low"
            if np.isfinite(m["mean"]) and m["mean"] > 0
            else ("low>high" if np.isfinite(m["mean"]) and m["mean"] < 0 else "flat/NA")
        )
        rows.append(
            {
                "family": family_key,
                "row": "pair",
                "interaction_name": name,
                "ligand": rec["ligand"],
                "receptor": rec["receptor"],
                "axis": rec["axis"],
                "thesis_expect": rec["thesis_expect"],
                "split": mode,
                "n_units": m["n"],
                "n_gse123902": m["n_123902"],
                "n_gse131907": m["n_131907"],
                "n_gse205335": m["n_205335"],
                "n_gse189357": m["n_189357"],
                "mean_delta": m["mean"],
                "median_delta": m["median"],
                "se": m["se"],
                "wilcoxon_p": m["wilcoxon_p"],
                "meta_p": m["meta_p"],
                "I2": m["i2"],
                "n_pos": m["n_pos"],
                "n_neg": m["n_neg"],
                "observed_direction": observed,
                "agrees_thesis": agrees(rec["thesis_expect"], m["mean"], m["n"], m["wilcoxon_p"]),
            }
        )
        if not sub.empty:
            pair_long.append(sub)

    if pair_long:
        cat = pd.concat(pair_long, ignore_index=True)
        det = cat[cat["detected"]].copy()
        if not det.empty:
            fam_unit = det.groupby(["cohort", "patient_id"], as_index=False).agg(
                delta=("delta", "mean"), n_pairs=("delta", "size")
            )
            fam_unit["detected"] = True
            m = meta_pair(fam_unit)
            observed = (
                "high>low"
                if np.isfinite(m["mean"]) and m["mean"] > 0
                else ("low>high" if np.isfinite(m["mean"]) and m["mean"] < 0 else "flat/NA")
            )
            rows.append(
                {
                    "family": family_key,
                    "row": "FAMILY_AGGREGATE",
                    "interaction_name": f"FAMILY:{family_key}",
                    "ligand": "|".join(sorted({r["ligand"] for r in spec["pairs"]})),
                    "receptor": "T/NK",
                    "axis": "family",
                    "thesis_expect": spec["thesis_expect"],
                    "split": mode,
                    "n_units": m["n"],
                    "n_gse123902": m["n_123902"],
                    "n_gse131907": m["n_131907"],
                    "n_gse205335": m["n_205335"],
                    "n_gse189357": m["n_189357"],
                    "mean_delta": m["mean"],
                    "median_delta": m["median"],
                    "se": m["se"],
                    "wilcoxon_p": m["wilcoxon_p"],
                    "meta_p": m["meta_p"],
                    "I2": m["i2"],
                    "n_pos": m["n_pos"],
                    "n_neg": m["n_neg"],
                    "observed_direction": observed,
                    "agrees_thesis": agrees(spec["thesis_expect"], m["mean"], m["n"], m["wilcoxon_p"]),
                }
            )
    return pd.DataFrame(rows)


def write_family_md_rows(df: pd.DataFrame) -> str:
    lines = [
        "| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for rec in df.itertuples(index=False):
        name = rec.interaction_name.replace("FAMILY:", "**FAMILY** ")
        lines.append(
            f"| {name} | {rec.axis} | {rec.thesis_expect} | {rec.n_units} | "
            f"{rec.n_gse123902} | {rec.n_gse131907} | {rec.n_gse205335} | {rec.n_gse189357} | "
            f"{fmt_num(rec.mean_delta)} | {fmt_p(rec.wilcoxon_p)} | "
            f"{rec.observed_direction} | {rec.agrees_thesis} |"
        )
    return "\n".join(lines)


def cxcl_line(fam_ifn: pd.DataFrame, long_df: pd.DataFrame, mode: str) -> str:
    names = ["CXCL9_CXCR3", "CXCL10_CXCR3"]
    bits = []
    for name in names:
        sub = long_df[(long_df["interaction_name"] == name) & (long_df["split"] == mode)]
        n_det = int(sub["detected"].sum()) if not sub.empty else 0
        row = fam_ifn[fam_ifn["interaction_name"] == name]
        if row.empty:
            bits.append(f"{name} undetected")
            continue
        r = row.iloc[0]
        if r.n_units < 3:
            bits.append(
                f"{name} undetected at expr_prop≥{EXPR_PROP:g} (honest n={int(r.n_units)}; "
                f"units with any detect flag={n_det})"
            )
        else:
            bits.append(
                f"{name} n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} p_W={fmt_p(r.wilcoxon_p)}"
            )
    return "CXCL9/10 (KD-like recruit arm): " + "; ".join(bits) + "."


def fig_given(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    names = [r["cohort"] for r in GIVEN_SINGLES]
    rhos = [r["rho"] for r in GIVEN_SINGLES]
    ns = [r["n"] for r in GIVEN_SINGLES]
    colors = [COHORT_COLORS[c] for c in names]
    ax.barh(names[::-1], rhos[::-1], color=colors[::-1])
    ax.axvline(0, color="k", lw=0.8)
    for i, (n, r) in enumerate(zip(ns[::-1], rhos[::-1])):
        ax.text(
            r - 0.02 if r < 0 else r + 0.02,
            i,
            f"n={n}  ρ={r:.3f}",
            va="center",
            ha="right" if r < 0 else "left",
            fontsize=8,
        )
    ax.set_xlabel("Spearman ρ (malignant CLDN4 %pos vs T/NK frac)")
    ax.set_title("PR #459 concordant-four singles — not re-audited")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_honest_n(cov: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    labels = [
        "locked four (given)",
        "GSE123902 locked",
        "GSE131907 locked",
        "GSE205335 locked",
        "GSE189357 locked",
        "Q4 vs Q1 LR (primary)",
        "%pos LR extra",
        "median LR extra",
        "tertile LR extra",
    ]
    vals = [
        sum(GIVEN_LOCKED.values()),
        int((cov["cohort"] == "GSE123902").sum()),
        int((cov["cohort"] == "GSE131907").sum()),
        int((cov["cohort"] == "GSE205335").sum()),
        int((cov["cohort"] == "GSE189357").sum()),
        int(cov["q4q1_ok"].sum()),
        int(cov["pctpos_ok"].sum()) if "pctpos_ok" in cov.columns else 0,
        int(cov["median_ok"].sum()),
        int(cov["tertile_ok"].sum()),
    ]
    ax.barh(labels[::-1], vals[::-1], color="#4C78A8")
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("honest n (patient / donor / locked sample)")
    ax.set_title("Honest n = units with both compartments")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_n_per_unit(cov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.2), sharey=False)
    for ax, cohort in zip(axes.ravel(), COHORTS):
        sub = cov[cov["cohort"] == cohort].sort_values("patient_id")
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["n_mal"], width=0.4, label="malignant / epithelium", color="#4C78A8")
        ax.bar(x + 0.2, sub["n_tnk"], width=0.4, label="T/NK", color="#F58518")
        ax.axhline(MIN_CELLS_ARM, color="k", ls="--", lw=0.7)
        ax.axhline(MIN_CELLS_COMP, color="0.4", ls=":", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient_id"], rotation=90, fontsize=5)
        ax.set_title(cohort)
        ax.set_ylabel("cells")
        ax.legend(fontsize=6)
    fig.suptitle("Per-unit malignant/epithelium and T/NK floors")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_bars(fam: pd.DataFrame, title: str, path: Path) -> None:
    plot = fam[fam["row"] == "pair"].copy()
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    y = np.arange(len(plot))
    colors = ["#54A24B" if d > 0 else "#E45756" for d in plot["mean_delta"].fillna(0)]
    ax.barh(y, plot["mean_delta"].fillna(0), color=colors)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["interaction_name"], fontsize=8)
    ax.set_xlabel("mean patient ΔP (CLDN4-high − CLDN4-low → T/NK)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_forest(long_df: pd.DataFrame, family_key: str, mode: str, path: Path) -> None:
    spec = FAMILIES[family_key]
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    ylabels, means, ses = [], [], []
    for rec in spec["pairs"]:
        sub = long_df[(long_df["interaction_name"] == rec["interaction_name"]) & (long_df["split"] == mode)]
        m = meta_pair(sub)
        ylabels.append(rec["interaction_name"])
        means.append(m["mean"] if np.isfinite(m["mean"]) else 0.0)
        ses.append(m["se"] if np.isfinite(m["se"]) else 0.0)
    y = np.arange(len(ylabels))
    ax.errorbar(means, y, xerr=ses, fmt="o", color="#4C78A8", capsize=3)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("mean ΔP ± SE (patient/donor unit)")
    ax.set_title(f"{spec['label']} — {mode}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_patient(long_df: pd.DataFrame, family_key: str, mode: str, path: Path) -> None:
    spec = FAMILIES[family_key]
    names = {r["interaction_name"] for r in spec["pairs"]}
    sub = long_df[
        (long_df["interaction_name"].isin(names)) & (long_df["split"] == mode) & (long_df["detected"])
    ]
    if sub.empty:
        return
    fam = sub.groupby(["cohort", "patient_id"], as_index=False).agg(delta=("delta", "mean"))
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.2), sharex=True)
    for ax, cohort in zip(axes.ravel(), COHORTS):
        part = fam[fam["cohort"] == cohort].sort_values("delta")
        y = np.arange(len(part))
        ax.scatter(part["delta"], y, color=COHORT_COLORS[cohort], s=22)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(part["patient_id"], fontsize=5)
        ax.set_title(cohort)
        ax.set_xlabel("family mean ΔP")
    fig.suptitle(f"Per-unit family score — {family_key} / {mode}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_cohort_family(long_df: pd.DataFrame, mode: str, path: Path) -> None:
    rows = []
    for family_key in FAMILIES:
        spec = FAMILIES[family_key]
        names = {r["interaction_name"] for r in spec["pairs"]}
        sub = long_df[
            (long_df["interaction_name"].isin(names))
            & (long_df["split"] == mode)
            & (long_df["detected"])
        ]
        if sub.empty:
            continue
        fam = sub.groupby(["cohort", "patient_id"], as_index=False).agg(delta=("delta", "mean"))
        for cohort, part in fam.groupby("cohort"):
            m = unit_delta_stats(part["delta"])
            rows.append({"family": family_key, "cohort": cohort, "n": m["n"], "mean": m["mean"], "se": m["se"]})
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), sharey=True)
    for ax, family_key, title in zip(
        axes,
        ["barrier_inhibitory", "ifn_recruit_mhci"],
        ["Barrier/inhibitory (expect +)", "IFN/recruit/MHC-I (expect −)"],
    ):
        part = df[df["family"] == family_key]
        x = np.arange(len(COHORTS))
        means = [
            float(part.loc[part["cohort"] == c, "mean"].iloc[0]) if (part["cohort"] == c).any() else 0
            for c in COHORTS
        ]
        ses = [
            float(part.loc[part["cohort"] == c, "se"].iloc[0]) if (part["cohort"] == c).any() else 0
            for c in COHORTS
        ]
        colors = [COHORT_COLORS[c] for c in COHORTS]
        ax.bar(x, means, yerr=ses, color=colors, capsize=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(COHORTS, rotation=20, fontsize=8)
        ax.set_title(title)
        ax.set_ylabel("mean family ΔP")
    fig.suptitle(f"Family ΔP by cohort ({mode}; patient unit)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_called_out(fam_bar: pd.DataFrame, path: Path) -> None:
    """Call out F11R, NECTIN2–TIGIT, CDH1, LGALS9."""
    axes_order = ["F11R", "NECTIN2-TIGIT", "CDH1", "LGALS9"]
    plot = fam_bar[fam_bar["row"] == "pair"].copy()
    rows = []
    for axis in axes_order:
        part = plot[plot["axis"] == axis]
        if part.empty:
            continue
        # mean of pair mean_deltas weighted by n
        w = part["n_units"].to_numpy(float)
        d = part["mean_delta"].to_numpy(float)
        mask = np.isfinite(d) & (w > 0)
        if not mask.any():
            continue
        rows.append(
            {
                "axis": axis,
                "mean": float(np.average(d[mask], weights=w[mask])),
                "n": int(w[mask].max()),
            }
        )
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    colors = ["#54A24B" if v > 0 else "#E45756" for v in df["mean"]]
    ax.barh(df["axis"][::-1], df["mean"][::-1], color=colors[::-1])
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("n-weighted mean ΔP (high − low → T/NK)")
    ax.set_title("Called-out barrier/inhibitory ligands (ON-thesis)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_comparison(fam_bar: pd.DataFrame, path: Path) -> None:
    agg = fam_bar[fam_bar["row"] == "FAMILY_AGGREGATE"]
    this_n = int(agg.iloc[0].n_units) if not agg.empty else 0
    this_d = float(agg.iloc[0].mean_delta) if not agg.empty else 0.0
    labels = [f"{r['source']}\n{r['row']}" for r in COMPARISON_ROWS] + ["this\nFAMILY barrier"]
    vals = [r["delta"] for r in COMPARISON_ROWS] + [this_d]
    colors = ["#9D9D9D", "#9D9D9D", "#54A24B"]
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    ax.bar(np.arange(len(labels)), vals, color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("ΔP (high − low)")
    ax.set_title("Comparison rows (PR #424/#474 not re-audited) vs four-set")
    for i, (v, n) in enumerate(zip(vals, [r["n"] for r in COMPARISON_ROWS] + [this_n])):
        ax.text(i, v + 0.004, f"n={n}\n{v:+.3f}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(
    cov: pd.DataFrame,
    fam_bar: pd.DataFrame,
    fam_ifn: pd.DataFrame,
    fam_bar_med: pd.DataFrame,
    long_df: pd.DataFrame,
    cxcl: str,
) -> str:
    n_q = int(cov["q4q1_ok"].sum())
    n_pct = int(cov["pctpos_ok"].sum()) if "pctpos_ok" in cov.columns else 0
    n_med = int(cov["median_ok"].sum())
    n_ter = int(cov["tertile_ok"].sum())
    counts = {c: int((cov["cohort"] == c).sum()) for c in COHORTS}
    bar_agg = fam_bar[fam_bar["row"] == "FAMILY_AGGREGATE"]
    ifn_agg = fam_ifn[fam_ifn["row"] == "FAMILY_AGGREGATE"]
    bar_txt = "NA"
    ifn_txt = "NA"
    if not bar_agg.empty:
        r = bar_agg.iloc[0]
        bar_txt = (
            f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} "
            f"p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"
        )
    if not ifn_agg.empty:
        r = ifn_agg.iloc[0]
        ifn_txt = (
            f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} "
            f"p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"
        )
    both = int(((cov["n_mal"] >= MIN_CELLS_ARM) & (cov["n_tnk"] >= MIN_CELLS_COMP)).sum())

    md = f"""# FINDING — Concordant-four CellChat barrier call (CLDN4-only)

ADDITIVE. **Thesis already correct. Ligands stay.**
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.

Patient/donor/locked-sample is the unit. Honest ligand n = units with
**both compartments** (malignant/epithelium high+low arms and T/NK).

Thesis (already correct; this folder reports the four-set pooled ligand call):

- Barrier/inhibitory outgoing **UP from CLDN4-high** (F11R, NECTIN2–TIGIT, CDH1, LGALS9) is **ON-thesis**.
- IFN/T-recruit outgoing UP from CLDN4-low is the KD-like arm (often **not** detected; CXCL9/10 sparse).
  Do not bury barrier-up-in-high as a recruit-up skip.

Primary split is malignant **Q4 vs Q1**. Extra: %pos and median.

## Honest n

Locked four n=65 (13+21+22+9) is **not** the ligand n.

| gate | n | note |
|---|---:|---|
| Locked four (given PR #459 singles) | 65 | 13+21+22+9; T/NK ρ not re-audited |
| Locked GSE123902 donors | {counts['GSE123902']} | marker-malignant / epithelium proxy; PRIMARY preferred |
| Locked GSE131907 samples | {counts['GSE131907']} | author Malignant cells / tS*; T lymphocytes + NK cells |
| Locked GSE205335 patients | {counts['GSE205335']} | author lineage.sub Malignant cells; lineage.total T/NK |
| Locked GSE189357 patients | {counts['GSE189357']} | marker-malignant / epithelium proxy; TD1–TD9 |
| Both compartments (n_mal≥{MIN_CELLS_ARM}, n_tnk≥{MIN_CELLS_COMP}) | {both} | honest inventory |
| Q4 vs Q1 LR | **{n_q}** | primary ligand n (n_mal≥{MIN_MAL_Q4}, ≥{MIN_CELLS_ARM}/arm, T/NK≥{MIN_CELLS_COMP}) |
| %pos LR extra | {n_pct} | CLDN4>0 vs =0 |
| Median-split LR extra | {n_med} | |
| Tertile extra | {n_ter} | exclusive arms |

GSE123902 / GSE189357 malignant labels are thin → epithelium marker-malignant
(EPCAM\\|KRT8\\|KRT18\\|KRT19 > 0 and PTPRC == 0) → T/NK.
GSE131907 / GSE205335 use author malignant → T/NK.
TACSTD2 is never a gate. PVRL2 is aliased to NECTIN2 on GSE131907.

## Primary — barrier / inhibitory family ΔP (high − low)

ON-thesis. Primary table: `results/family_barrier_inhibitory.tsv` (Q4 vs Q1).

{write_family_md_rows(fam_bar)}

Family aggregate: {bar_txt}.

Called-out axes (F11R, NECTIN2–TIGIT, CDH1, LGALS9) are the high-side barrier ligands.

## KD-like arm — IFN / T-recruit / MHC-I (expect low > high)

Do not file the barrier result as a recruit-up skip.

{cxcl}

Primary table: `results/family_ifn_recruit_mhci.tsv`.

{write_family_md_rows(fam_ifn)}

Family aggregate: {ifn_txt}.

## Extra: median-split barrier family (same units, different gate)

{write_family_md_rows(fam_bar_med)}

## Comparison rows (PR #424 / #474 — not re-audited)

The new result is the four-set pooled ligand call. These rows are context only.

| source | combo | row | split | n | ΔP | p | note |
|---|---|---|---|---:|---:|---|---|
| PR #424 | GSE131907+GSE205335 | NECTIN2–TIGIT | median | 39 | +0.094 | 3.6e-12 | given |
| PR #474 | GSE123902+GSE205335 | FAMILY barrier | median | 32 | +0.087 | 2.00e-08 | given |
| this | concordant four | FAMILY barrier | Q4 vs Q1 | {int(bar_agg.iloc[0].n_units) if not bar_agg.empty else 0} | {fmt_num(bar_agg.iloc[0].mean_delta) if not bar_agg.empty else 'NA'} | {fmt_p(bar_agg.iloc[0].wilcoxon_p) if not bar_agg.empty else 'NA'} | new |

## Extra figures

- `figures/fig_given_singles_rho.png` — PR #459 singles (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_called_out_barrier.png` — F11R, NECTIN2–TIGIT, CDH1, LGALS9
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`
- `figures/fig_family_by_cohort.png`
- `figures/fig_comparison_424_474.png` — comparison rows vs this call

## What is not claimed

- PR #424 / #474 numbers are comparison rows. They were not recomputed.
- PR #459 T/NK Spearmans are given. They were not re-audited.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071, GSE127465, and CD45-only libraries are not added.
- This is not a CellChat discovery screen and not a 7-pool.
- CellChat R / LIANA R were not run. Probability is Jin et al. 2021
  (10% truncated mean, Kh={KH}, expr_prop≥{EXPR_PROP}).

## Reproduce

```bash
python3 methods/concordant4_cellchat_barrier_cldn4/scripts/download.py
python3 methods/concordant4_cellchat_barrier_cldn4/scripts/analyze.py
```

Hill constants: trim={TRIM}, Kh={KH}, expr_prop={EXPR_PROP}.
Arm floor ≥{MIN_CELLS_ARM}; T/NK compartment floor ≥{MIN_CELLS_COMP}; Q4 vs Q1 also n_mal≥{MIN_MAL_Q4}.
"""
    return md


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/concordant4_raw"))
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    locked = load_locked()
    print("locked " + " ".join(f"{k} n={len(locked[k])}" for k in COHORTS), flush=True)
    wanted = thesis_genes()
    names = {r["interaction_name"] for fam in FAMILIES.values() for r in fam["pairs"]}
    lr = load_lr(DB, names, matrix_genes=None)
    print(f"pre-specified pairs in CellChatDB: {len(lr)} / {len(names)}", flush=True)
    missing = names - set(lr["interaction_name"])
    if missing:
        print("WARNING missing interactions:", sorted(missing), flush=True)

    u123, genes123 = load_gse123902(args.raw, wanted, locked["GSE123902"])
    u131, genes131, n_streamed = load_gse131907(args.raw, wanted, locked["GSE131907"])
    u205, genes205 = load_gse205335(args.raw, wanted, locked["GSE205335"])
    u189, genes189 = load_gse189357(args.raw, wanted, locked["GSE189357"])
    units = u123 + u131 + u205 + u189
    print(
        f"units loaded: {len(units)}  genes "
        f"123902={len(genes123)} 131907={len(genes131)} "
        f"205335={len(genes205)} 189357={len(genes189)}",
        flush=True,
    )

    coverage_rows = []
    long_rows = []
    for unit in units:
        rec = {
            "cohort": unit["cohort"],
            "patient_id": unit["patient_id"],
            "n_cells": unit["n_cells"],
            "n_mal": unit["n_mal"],
            "n_tnk": unit["n_tnk"],
            "malig_def": unit.get("malig_def", ""),
            "both_compartments": bool(
                unit["n_mal"] >= MIN_CELLS_ARM and unit["n_tnk"] >= MIN_CELLS_COMP
            ),
        }
        for mode in SPLITS:
            rec[f"{mode}_ok"] = False
            scored = score_unit(unit, lr, mode)
            if scored is None:
                continue
            rec[f"{mode}_ok"] = True
            rec[f"{mode}_n_high"] = scored["n_mal_high"]
            rec[f"{mode}_n_low"] = scored["n_mal_low"]
            long_rows.extend(scored["pairs"])
        coverage_rows.append(rec)

    cov = pd.DataFrame(coverage_rows)
    long_df = pd.DataFrame(long_rows)
    cov.to_csv(RESULTS / "patient_coverage.tsv", sep="\t", index=False)
    long_df.to_csv(RESULTS / "patient_lr_long.tsv", sep="\t", index=False)

    fam_bar = family_table(long_df, "barrier_inhibitory", PRIMARY_SPLIT)
    fam_ifn = family_table(long_df, "ifn_recruit_mhci", PRIMARY_SPLIT)
    fam_bar_med = family_table(long_df, "barrier_inhibitory", "median")
    fam_bar.to_csv(RESULTS / "family_barrier_inhibitory.tsv", sep="\t", index=False)
    fam_ifn.to_csv(RESULTS / "family_ifn_recruit_mhci.tsv", sep="\t", index=False)
    fam_bar_med.to_csv(RESULTS / "family_barrier_inhibitory_median.tsv", sep="\t", index=False)

    extras = [family_table(long_df, key, mode) for key in FAMILIES for mode in SPLITS]
    pd.concat(extras, ignore_index=True).to_csv(RESULTS / "family_all_splits.tsv", sep="\t", index=False)

    # Called-out four-axis table (primary split).
    call = fam_bar[fam_bar["axis"].isin(["F11R", "NECTIN2-TIGIT", "CDH1", "LGALS9", "family"])].copy()
    call.to_csv(RESULTS / "barrier_called_out.tsv", sep="\t", index=False)

    # Comparison table (given rows + this family).
    agg = fam_bar[fam_bar["row"] == "FAMILY_AGGREGATE"]
    comp = pd.DataFrame(COMPARISON_ROWS)
    if not agg.empty:
        r = agg.iloc[0]
        comp = pd.concat(
            [
                comp,
                pd.DataFrame(
                    [
                        {
                            "source": "this",
                            "combo": "GSE123902+GSE131907+GSE205335+GSE189357",
                            "row": "FAMILY barrier",
                            "split": PRIMARY_SPLIT,
                            "n": int(r.n_units),
                            "n_note": f"{int(r.n_gse123902)}/{int(r.n_gse131907)}/{int(r.n_gse205335)}/{int(r.n_gse189357)}",
                            "delta": float(r.mean_delta),
                            "p": float(r.wilcoxon_p),
                            "metric": "mean ΔP",
                            "note": "new four-set pooled ligand call",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    comp.to_csv(RESULTS / "comparison_pr424_pr474.tsv", sep="\t", index=False)

    cxcl = cxcl_line(fam_ifn, long_df, PRIMARY_SPLIT)

    fig_given(FIGURES / "fig_given_singles_rho")
    fig_honest_n(cov, FIGURES / "fig_honest_n")
    fig_n_per_unit(cov, FIGURES / "fig_n_per_unit")
    fig_family_bars(fam_bar, "Barrier/inhibitory outgoing ΔP (Q4 vs Q1)", FIGURES / "fig_family_barrier_bars")
    fig_family_bars(fam_ifn, "IFN/T-recruit/MHC-I outgoing ΔP (Q4 vs Q1)", FIGURES / "fig_family_ifn_bars")
    fig_family_forest(long_df, "barrier_inhibitory", PRIMARY_SPLIT, FIGURES / "fig_family_barrier_forest")
    fig_family_forest(long_df, "ifn_recruit_mhci", PRIMARY_SPLIT, FIGURES / "fig_family_ifn_forest")
    fig_family_patient(long_df, "barrier_inhibitory", PRIMARY_SPLIT, FIGURES / "fig_family_barrier_units")
    fig_family_patient(long_df, "ifn_recruit_mhci", PRIMARY_SPLIT, FIGURES / "fig_family_ifn_units")
    fig_cohort_family(long_df, PRIMARY_SPLIT, FIGURES / "fig_family_by_cohort")
    fig_called_out(fam_bar, FIGURES / "fig_called_out_barrier")
    fig_comparison(fam_bar, FIGURES / "fig_comparison_424_474")

    summary = {
        "combo": "GSE123902+GSE131907+GSE205335+GSE189357",
        "full_pool": False,
        "excluded": ["GSE148071", "GSE127465", "CD45-only"],
        "dual_high": False,
        "primary_split": PRIMARY_SPLIT,
        "n_locked": {k: int((cov.cohort == k).sum()) for k in COHORTS},
        "n_lr": {mode: int(cov[f"{mode}_ok"].sum()) for mode in SPLITS if f"{mode}_ok" in cov.columns},
        "family_tables": [
            "results/family_barrier_inhibitory.tsv",
            "results/family_ifn_recruit_mhci.tsv",
        ],
        "n_streamed_131907_genes": int(n_streamed),
        "re_audited_pr424_pr474": False,
        "re_audited_tnk_rho": False,
        "cxcl910": cxcl,
        "pvrl2_aliased_to_nectin2": True,
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    finding = write_finding(cov, fam_bar, fam_ifn, fam_bar_med, long_df, cxcl)
    (ROOT / "FINDING.md").write_text(finding)
    (RESULTS / "FINDING.md").write_text(finding)
    print("wrote", ROOT / "FINDING.md", flush=True)
    print("barrier family table:", RESULTS / "family_barrier_inhibitory.tsv", flush=True)
    print(fam_bar.to_string(index=False), flush=True)
    print(cxcl, flush=True)
    print(fam_ifn.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
