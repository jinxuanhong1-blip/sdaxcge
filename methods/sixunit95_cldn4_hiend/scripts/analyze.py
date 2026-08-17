#!/usr/bin/env python3
"""CLDN4-only high-end CellChat-style outgoing Mal → T/NK on the locked 6-unit combo.

Patient is the unit. The n=95 ρ=−0.260 Spearman row is taken as given and is
not re-audited. GSE253013 9.3 GB RDS is not downloaded.
"""
from __future__ import annotations

import argparse
import gzip
import json
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

from lib import (
    EXTRA,
    LINEAGE_MARKERS,
    MIN_CELLS_ARM,
    MIN_CELLS_COMP,
    MIN_MAL_Q4,
    NORMAL_LUNG,
    assign_patient_quartiles,
    dersimonian_laird,
    fmt_num,
    fmt_p,
    highend_split,
    ligand_class,
    load_lr,
    marker_lineage,
    marker_malignant,
    parse_symbols,
    score_outgoing,
    unit_delta_stats,
    wanted_genes,
)

ROOT = Path(__file__).resolve().parents[1]
LOCKED = ROOT / "data" / "locked"
DB = ROOT / "db"
GIVEN_RHO = -0.260
GIVEN_N = 95
GIVEN_P = 0.0195
GIVEN_I2 = 0.0

UNITS = [
    "GSE207422",
    "GSE205335",
    "GSE291670",
    "GSE253013",
    "GSE131907",
    "GSE325414",
]


def stream_gene_text(path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = np.array(header[1:], dtype=object)
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        all_genes: list[str] = []
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            all_genes.append(gene)
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, set(all_genes)


def stream_mtx(mtx_gz: Path, n_genes: int, n_cells: int, keep_idx: list[int]) -> sparse.csr_matrix:
    keep_set = set(keep_idx)
    pos = {g: i for i, g in enumerate(keep_idx)}
    rows, cols, vals = [], [], []
    with gzip.open(mtx_gz, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            break
        for line in handle:
            r, c, v = line.split()
            gi = int(r) - 1
            if gi in keep_set:
                rows.append(int(c) - 1)
                cols.append(pos[gi])
                vals.append(float(v))
    return sparse.coo_matrix((vals, (rows, cols)), shape=(n_cells, len(keep_idx))).tocsr()


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def load_locked_patients() -> dict[str, pd.DataFrame]:
    out = {}
    g = pd.read_csv(LOCKED / "GSE207422_drmref_patients.tsv", sep="\t")
    g["unit"] = "GSE207422"
    g["patient_id"] = g["sample"]
    g["cldn4_given"] = g["malig_CLDN4_mean"]
    out["GSE207422"] = g

    h = pd.read_csv(LOCKED / "GSE205335_patients.tsv", sep="\t")
    h["unit"] = "GSE205335"
    h["patient_id"] = h["patient"]
    h["cldn4_given"] = h["mal_CLDN4_mean"]
    out["GSE205335"] = h

    k = pd.read_csv(LOCKED / "GSE291670_patients.tsv", sep="\t")
    k["unit"] = "GSE291670"
    k["patient_id"] = k["sample"]
    k["cldn4_given"] = k["mal_CLDN4_mean_log1p_cp10k"]
    out["GSE291670"] = k

    a = pd.read_csv(LOCKED / "GSE253013_patients.tsv", sep="\t")
    a = a[(a["tissue"] == "Tumor") & (a["eligible_malig"])].copy()
    a["unit"] = "GSE253013"
    a["patient_id"] = a["patient"]
    a["cldn4_given"] = a["CLDN4_mean_log1p_cp10k"]
    out["GSE253013"] = a

    s = pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t")
    s = s[s["n_malignant"] >= 20].copy()
    s["unit"] = "GSE131907"
    s["patient_id"] = s["sample"]
    s["cldn4_given"] = s["mal_CLDN4_mean"]
    out["GSE131907"] = s

    d = pd.read_csv(LOCKED / "GSE325414_donors.csv")
    d["unit"] = "GSE325414"
    d["patient_id"] = d["donor"]
    d["cldn4_given"] = d["malignant_CLDN4_mean"]
    out["GSE325414"] = d

    n = sum(len(v) for v in out.values())
    if n != GIVEN_N:
        raise SystemExit(f"locked n={n} != given n={GIVEN_N}")
    return out


def load_gse207422(raw: Path, wanted: set[str]):
    print("==== GSE207422 ====", flush=True)
    locked = pd.read_csv(LOCKED / "GSE207422_drmref_patients.tsv", sep="\t")
    keep = set(locked["sample"])
    path = raw / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    cell_ids, extracted, n_umi, genes = stream_gene_text(path, wanted)
    sample = np.array(["_".join(str(bc).split("_")[:2]) for bc in cell_ids], dtype=object)
    log_cp, pos = to_log_pos(extracted, n_umi)
    lineage = marker_lineage(log_cp, len(cell_ids))
    mal = lineage == "Epithelial"  # DRMref barcodes are not public
    tnk = np.isin(lineage, ["T", "NK"])
    keep_mask = np.isin(sample, list(keep))
    print(f"  locked samples={len(keep)} epi={int((mal & keep_mask).sum())} tnk={int((tnk & keep_mask).sum())}", flush=True)
    return {
        "patient": sample,
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": genes,
        "keep": keep,
        "note": "DRMref barcodes not public; epithelial lineage is the malignant proxy on the same 12 locked samples",
    }


def load_gse205335(raw: Path, wanted: set[str]):
    print("==== GSE205335 ====", flush=True)
    import rdata

    ident = pd.read_csv(raw / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    rds = raw / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        unzipped = Path(tmp) / "matrix.rds"
        print("  decompress RDS", flush=True)
        with gzip.open(rds, "rb") as src, unzipped.open("wb") as dest:
            dest.write(src.read())
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(unzipped)
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False).tocsr()
    library = np.asarray(matrix.sum(axis=0)).ravel()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is not None:
            extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
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
    mal = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    tnk = (
        cells["lineage.total"].eq("T/NK cells").to_numpy()
        if "lineage.total" in cells.columns
        else np.zeros(len(cells), dtype=bool)
    )
    patient = cells["patient"].astype(str).to_numpy()
    locked = set(pd.read_csv(LOCKED / "GSE205335_patients.tsv", sep="\t")["patient"].astype(str))
    log_cp, pos = to_log_pos(extracted, library)
    print(f"  mal={int(mal.sum())} tnk={int(tnk.sum())} patients_in_matrix={pd.Series(patient).nunique()}", flush=True)
    return {
        "patient": np.asarray(patient, dtype=object),
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": set(genes.tolist()),
        "keep": locked,
        "note": "author lineage.sub == Malignant cells; lineage.total == T/NK cells",
    }


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


def load_gse291670(raw: Path, wanted: set[str]):
    print("==== GSE291670 ====", flush=True)
    import scipy.io

    tar_path = raw / "GSE291670" / "GSE291670_RAW.tar"
    extract_dir = raw / "GSE291670" / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)
    if not any(extract_dir.glob("*_matrix.mtx.gz")):
        with tarfile.open(tar_path) as tar:
            tar.extractall(extract_dir)
    samples = [
        "GSM8839599_MPR-1",
        "GSM8839600_MPR-2",
        "GSM8839601_MPR-3",
        "GSM8839602_Non-MPR-1",
        "GSM8839603_Non-MPR-2",
        "GSM8839604_Non-MPR-3",
    ]
    pats, mal_l, tnk_l = [], [], []
    extracted_all: dict[str, list[np.ndarray]] = {}
    libs = []
    genes_union: set[str] = set()
    for prefix in samples:
        matches = list(extract_dir.glob(f"{prefix}*_matrix.mtx.gz")) + list(extract_dir.glob(f"*{prefix}*_matrix.mtx.gz"))
        if not matches:
            # GEO RAW tar often uses GSM_xxx_matrix.mtx.gz
            matches = list(extract_dir.glob(f"{prefix.split('_')[0]}*_matrix.mtx.gz"))
        if not matches:
            raise FileNotFoundError(f"no MTX for {prefix} in {extract_dir}")
        mtx_path = matches[0]
        stem = str(mtx_path).replace("_matrix.mtx.gz", "")
        mtx = scipy.io.mmread(stem + "_matrix.mtx.gz").tocsr().astype(np.float32)
        with gzip.open(stem + "_barcodes.tsv.gz", "rt") as fh:
            barcodes = [line.strip().split("\t")[0] for line in fh]
        genes = []
        with gzip.open(stem + "_features.tsv.gz", "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                genes.append(parts[1] if len(parts) > 1 else parts[0])
        if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
            mtx = mtx.T.tocsr()
        genes_union.update(genes)
        name_to_row = {g: i for i, g in enumerate(genes)}
        library = np.asarray(mtx.sum(axis=0)).ravel()
        n = len(barcodes)
        sample_ext = {}
        for gene in wanted:
            row = name_to_row.get(gene)
            sample_ext[gene] = np.asarray(mtx.getrow(row).toarray()).ravel() if row is not None else np.zeros(n, np.float32)
        log_cp, _ = to_log_pos(sample_ext, library)
        lineage = marker_lineage(log_cp, n)
        mal = marker_malignant(lineage, log_cp)
        tnk = np.isin(lineage, ["T", "NK"])
        pats.append(np.array([prefix] * n, dtype=object))
        mal_l.append(mal)
        tnk_l.append(tnk)
        libs.append(library)
        for gene, vec in sample_ext.items():
            extracted_all.setdefault(gene, []).append(vec)
        print(f"  {prefix} cells={n} mal={int(mal.sum())} tnk={int(tnk.sum())}", flush=True)
    extracted = {g: np.concatenate(vs) for g, vs in extracted_all.items()}
    library = np.concatenate(libs)
    log_cp, pos = to_log_pos(extracted, library)
    locked = set(pd.read_csv(LOCKED / "GSE291670_patients.tsv", sep="\t")["sample"])
    return {
        "patient": np.concatenate(pats),
        "mal": np.concatenate(mal_l),
        "tnk": np.concatenate(tnk_l),
        "log_cp": log_cp,
        "pos": pos,
        "genes": genes_union,
        "keep": locked,
        "note": "marker malignant = epithelial minus high normal-lung score; no author labels on GEO",
    }


def load_gse131907(raw: Path, wanted: set[str]):
    print("==== GSE131907 ====", flush=True)
    annot = pd.read_csv(raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    path = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    cell_ids, extracted, n_umi, genes = stream_gene_text(path, wanted)
    # annotation index
    if "Index" in annot.columns:
        annot = annot.set_index("Index")
    elif annot.columns[0].lower() in {"cell", "barcode", "index"}:
        annot = annot.set_index(annot.columns[0])
    common = pd.Index(cell_ids).intersection(annot.index)
    if len(common) < 0.9 * len(cell_ids):
        # try Sample + barcode styles
        print(f"  warning annot overlap {len(common)} / {len(cell_ids)}", flush=True)
    order = pd.Index(cell_ids)
    per = annot.reindex(order)
    mal_sub = {"Malignant cells", "tS1", "tS2", "tS3"}
    mal = (per["Cell_type"].eq("Epithelial cells") & per["Cell_subtype"].isin(mal_sub)).fillna(False).to_numpy()
    tnk = per["Cell_type"].isin(["T lymphocytes", "NK cells"]).fillna(False).to_numpy()
    sample = per["Sample"].astype(str).to_numpy() if "Sample" in per.columns else np.array(cell_ids)
    locked = set(pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t").query("n_malignant >= 20")["sample"])
    log_cp, pos = to_log_pos(extracted, n_umi)
    print(f"  mal={int(mal.sum())} tnk={int(tnk.sum())} locked={len(locked)}", flush=True)
    return {
        "patient": np.asarray(sample, dtype=object),
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": genes,
        "keep": locked,
        "note": "author Epithelial ∩ {Malignant cells, tS1, tS2, tS3}; locked n_mal≥20 samples (n=21)",
    }


def load_gse325414(raw: Path, wanted: set[str]):
    print("==== GSE325414 ====", flush=True)
    feat = [line.strip().split("\t")[0] for line in gzip.open(raw / "GSE325414" / "GSE325414_features.tsv.gz", "rt")]
    barcodes = [line.strip().split("\t")[0] for line in gzip.open(raw / "GSE325414" / "GSE325414_barcodes.tsv.gz", "rt")]
    keep_names = [g for g in sorted(wanted) if g in feat]
    keep_idx = [feat.index(g) for g in keep_names]
    print(f"  streaming {len(keep_names)} genes from {len(feat)} x {len(barcodes)}", flush=True)
    X = stream_mtx(raw / "GSE325414" / "GSE325414_matrix.mtx.gz", len(feat), len(barcodes), keep_idx)
    meta = pd.read_csv(raw / "GSE325414" / "GSE325414_metadata_individual_cells.txt.gz", sep="\t")
    meta = meta.set_index("cell")
    common = pd.Index(barcodes).intersection(meta.index)
    bc_index = {b: i for i, b in enumerate(barcodes)}
    keep_cells = [bc_index[b] for b in common]
    X = X[keep_cells]
    meta = meta.loc[common]
    ncount = pd.to_numeric(meta["nCount_RNA"], errors="coerce").to_numpy()
    extracted = {g: np.asarray(X[:, j].todense()).ravel() for j, g in enumerate(keep_names)}
    tnk_major = {"Tcell_CD4", "Tcell_CD8", "NKcell", "Tcell_CD4-_CD8-", "Tcell_CD4+_CD8+"}
    mal_sub = {"EpithelialCells_TumorCells"}
    sub1 = meta["sub.pop.level1"].astype(str).to_numpy()
    major = meta["major.populations"].astype(str).to_numpy()
    mal = np.isin(sub1, list(mal_sub))
    tnk = np.isin(major, list(tnk_major)) | (sub1 == "Lymphoid")
    mal = np.isin(sub1, list(mal_sub))  # malignant wins
    log_cp, pos = to_log_pos(extracted, ncount)
    locked = set(pd.read_csv(LOCKED / "GSE325414_donors.csv")["donor"])
    print(f"  mal={int(mal.sum())} tnk={int(tnk.sum())} donors={meta['donor'].nunique()}", flush=True)
    return {
        "patient": meta["donor"].astype(str).to_numpy(),
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": set(feat),
        "keep": locked,
        "note": "author sub.pop.level1 == EpithelialCells_TumorCells; T/NK from major.populations",
    }


def score_bundle(unit: str, bundle: dict, lr: pd.DataFrame, modes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    patient = np.asarray(bundle["patient"], dtype=object)
    mal = np.asarray(bundle["mal"], dtype=bool)
    tnk = np.asarray(bundle["tnk"], dtype=bool)
    log_cp = bundle["log_cp"]
    pos = bundle["pos"]
    keep = set(bundle["keep"])
    cldn4 = log_cp.get("CLDN4")
    if cldn4 is None:
        raise SystemExit(f"{unit}: CLDN4 missing")
    n = patient.size
    idx_all = np.arange(n)
    coverage_rows = []
    pair_rows = []
    for pid in sorted(keep):
        mask = patient == pid
        mal_idx = idx_all[mask & mal]
        tnk_idx = idx_all[mask & tnk]
        rec = {
            "unit": unit,
            "patient_id": pid,
            "n_mal": int(mal_idx.size),
            "n_tnk": int(tnk_idx.size),
            "mean_cldn4_mal": float(cldn4[mal_idx].mean()) if mal_idx.size else np.nan,
            "eligible_allmal": int(mal_idx.size >= MIN_CELLS_COMP and tnk_idx.size >= MIN_CELLS_COMP),
        }
        # all-malignant outgoing (between-patient high-end uses this)
        if rec["eligible_allmal"]:
            scored = score_outgoing(lr, log_cp, pos, mal_idx, tnk_idx)
            for row in scored:
                row.update({"unit": unit, "patient_id": pid, "split": "all_mal", "arm": "all"})
                pair_rows.append(row)
        for mode in modes:
            split = highend_split(cldn4, mal_idx, mode)
            rec[f"eligible_{mode}"] = int(split is not None and tnk_idx.size >= MIN_CELLS_ARM)
            if split is None or tnk_idx.size < MIN_CELLS_ARM:
                rec[f"n_high_{mode}"] = 0
                rec[f"n_low_{mode}"] = 0
                continue
            hi, lo = split
            rec[f"n_high_{mode}"] = int(hi.size)
            rec[f"n_low_{mode}"] = int(lo.size)
            hi_rows = score_outgoing(lr, log_cp, pos, hi, tnk_idx)
            lo_rows = score_outgoing(lr, log_cp, pos, lo, tnk_idx)
            lo_map = {r["interaction_name"]: r for r in lo_rows}
            for hr in hi_rows:
                lr_name = hr["interaction_name"]
                lr_low = lo_map.get(lr_name)
                if lr_low is None:
                    continue
                pair_rows.append(
                    {
                        **{k: hr[k] for k in hr if k not in {"prob", "detected", "ligand_mean", "receptor_mean", "ligand_prop", "receptor_prop"}},
                        "unit": unit,
                        "patient_id": pid,
                        "split": mode,
                        "arm": "delta",
                        "prob_high": hr["prob"],
                        "prob_low": lr_low["prob"],
                        "delta": hr["prob"] - lr_low["prob"],
                        "detected_high": hr["detected"],
                        "detected_low": lr_low["detected"],
                        "detected_either": bool(hr["detected"] or lr_low["detected"]),
                    }
                )
        coverage_rows.append(rec)
        print(
            f"  {unit} {pid} mal={rec['n_mal']} tnk={rec['n_tnk']} "
            f"q4={rec.get('eligible_q4q1', 0)} med={rec.get('eligible_median', 0)}",
            flush=True,
        )
    return pd.DataFrame(coverage_rows), pd.DataFrame(pair_rows)


def meta_patient_deltas(pairs: pd.DataFrame, split: str) -> pd.DataFrame:
    sub = pairs[(pairs["split"] == split) & (pairs["arm"] == "delta") & (pairs["detected_either"])].copy()
    if sub.empty:
        return pd.DataFrame()
    rows = []
    for name, g in sub.groupby("interaction_name", sort=False):
        unit_stats = []
        for unit, ug in g.groupby("unit"):
            st = unit_delta_stats(ug["delta"].to_numpy())
            st["unit"] = unit
            if st["n"] >= 2 and np.isfinite(st["se"]) and st["se"] > 0:
                unit_stats.append(st)
        if not unit_stats:
            # still record patient-pooled Wilcoxon
            st_all = unit_delta_stats(g["delta"].to_numpy())
            first = g.iloc[0]
            rows.append(
                {
                    "interaction_name": name,
                    "pathway_name": first["pathway_name"],
                    "ligand": first["ligand"],
                    "receptor": first["receptor"],
                    "lr_class": first["lr_class"],
                    "split": split,
                    "k_units": 0,
                    "n_patients": st_all["n"],
                    "mean_delta": st_all["mean"],
                    "median_delta": st_all["median"],
                    "se": st_all["se"],
                    "p_meta": np.nan,
                    "p_wilcoxon_patients": st_all["wilcoxon_p"],
                    "i2": np.nan,
                    "units": "",
                }
            )
            continue
        meta = dersimonian_laird(
            np.array([s["mean"] for s in unit_stats]),
            np.array([s["se"] ** 2 for s in unit_stats]),
        )
        st_all = unit_delta_stats(g["delta"].to_numpy())
        first = g.iloc[0]
        rows.append(
            {
                "interaction_name": name,
                "pathway_name": first["pathway_name"],
                "ligand": first["ligand"],
                "receptor": first["receptor"],
                "lr_class": first["lr_class"],
                "split": split,
                "k_units": meta["k"],
                "n_patients": st_all["n"],
                "mean_delta": meta["mean"],
                "median_delta": st_all["median"],
                "se": meta["se"],
                "p_meta": meta["p"],
                "p_wilcoxon_patients": st_all["wilcoxon_p"],
                "i2": meta["i2"],
                "units": "+".join(sorted(s["unit"] for s in unit_stats)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["p_meta", "n_patients"], ascending=[True, False], na_position="last")


def meta_q4q1_between(pairs: pd.DataFrame, locked: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Between-patient high-end: Q4 vs Q1 of *given* CLDN4, patient-level all-mal P."""
    allmal = pairs[(pairs["split"] == "all_mal") & (pairs["detected"])].copy()
    if allmal.empty:
        return pd.DataFrame()
    qmap = {}
    for unit, df in locked.items():
        try:
            qs = assign_patient_quartiles(df["cldn4_given"])
        except ValueError:
            continue
        for pid, q in zip(df["patient_id"], qs.astype(str)):
            qmap[(unit, str(pid))] = q
    allmal["quartile"] = [qmap.get((u, str(p))) for u, p in zip(allmal["unit"], allmal["patient_id"])]
    rows = []
    for name, g in allmal.groupby("interaction_name", sort=False):
        unit_effects = []
        n_q1 = n_q4 = 0
        for unit, ug in g.groupby("unit"):
            a = ug.loc[ug["quartile"] == "Q1", "prob"]
            b = ug.loc[ug["quartile"] == "Q4", "prob"]
            if len(a) < 2 or len(b) < 2:
                continue
            n_q1 += int(len(a))
            n_q4 += int(len(b))
            u_stat, p = __import__("scipy").stats.mannwhitneyu(b, a, alternative="two-sided")
            r_rb = (2.0 * float(u_stat)) / (len(a) * len(b)) - 1.0
            delta = float(b.median() - a.median())
            # SE proxy from rank-biserial / n (descriptive)
            se = max(1e-6, abs(delta) / 1.96) if p < 1 else 0.1
            if len(a) + len(b) >= 4:
                se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))) if (a.var(ddof=1) + b.var(ddof=1)) > 0 else 0.05
            unit_effects.append({"unit": unit, "delta": delta, "se": se if se > 0 else 0.05, "r_rb": r_rb, "p": float(p), "n1": int(len(a)), "n4": int(len(b))})
        if not unit_effects:
            continue
        meta = dersimonian_laird(
            np.array([e["delta"] for e in unit_effects]),
            np.array([e["se"] ** 2 for e in unit_effects]),
        )
        first = g.iloc[0]
        rows.append(
            {
                "interaction_name": name,
                "pathway_name": first["pathway_name"],
                "ligand": first["ligand"],
                "receptor": first["receptor"],
                "lr_class": first.get("lr_class", ligand_class(first["ligand_genes"])),
                "split": "between_q4q1_given_cldn4",
                "k_units": meta["k"],
                "n_q1": n_q1,
                "n_q4": n_q4,
                "n_compared": n_q1 + n_q4,
                "mean_delta": meta["mean"],
                "se": meta["se"],
                "p_meta": meta["p"],
                "i2": meta["i2"],
                "units": "+".join(sorted(e["unit"] for e in unit_effects)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["p_meta", "n_compared"], ascending=[True, False], na_position="last")


def plot_n(coverage: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    units = [u for u in UNITS if u in set(coverage["unit"])]
    x = np.arange(len(units))
    locked_n = coverage.groupby("unit")["patient_id"].nunique().reindex(units).fillna(0)
    elig = coverage.groupby("unit")["eligible_q4q1"].sum().reindex(units).fillna(0) if "eligible_q4q1" in coverage else pd.Series(0, index=units)
    elig_m = coverage.groupby("unit")["eligible_median"].sum().reindex(units).fillna(0) if "eligible_median" in coverage else pd.Series(0, index=units)
    ax.bar(x - 0.25, locked_n, 0.24, label="locked n (given)", color="#9e9e9e")
    ax.bar(x, elig_m, 0.24, label="median-split eligible", color="#6a8aaa")
    ax.bar(x + 0.25, elig, 0.24, label="Q4 vs Q1 eligible", color="#b2182b")
    ax.set_xticks(x, units, rotation=25, ha="right")
    ax.set_ylabel("patients")
    ax.set_title("Honest n: locked combo vs CellChat-eligible patients")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(meta: pd.DataFrame, path: Path, title: str, n=12) -> None:
    show = meta.head(n).copy()
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * max(len(show), 1) + 1.6))
    if show.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No meta rows", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    y = np.arange(len(show))
    colors = np.where(show["mean_delta"] >= 0, "#b2182b", "#2166ac")
    ax.axvline(0, color="#444", lw=0.8)
    ax.barh(y, show["mean_delta"], color=colors, alpha=0.85)
    if "se" in show:
        ax.errorbar(show["mean_delta"], y, xerr=1.96 * show["se"].fillna(0), fmt="none", ecolor="#333", lw=0.8)
    labels = [
        f"{r.interaction_name}  n={int(r.n_patients) if 'n_patients' in r else int(r.get('n_compared', 0))}  p={fmt_p(r.p_meta)}"
        for r in show.itertuples()
    ]
    ax.set_yticks(y, labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("meta ΔP (CLDN4-high − low)")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_patient_deltas(pairs: pd.DataFrame, interaction: str, split: str, path: Path) -> None:
    sub = pairs[(pairs["split"] == split) & (pairs["arm"] == "delta") & (pairs["interaction_name"] == interaction)]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    if sub.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No patient deltas for {interaction}", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    units = [u for u in UNITS if u in set(sub["unit"])]
    data = [sub.loc[sub["unit"] == u, "delta"].to_numpy() for u in units]
    bp = ax.boxplot(data, tick_labels=units, patch_artist=True, widths=0.55)
    for patch in bp["boxes"]:
        patch.set_facecolor("#f4a582")
        patch.set_alpha(0.6)
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=14, zorder=3)
    ax.axhline(0, color="#444", lw=0.8)
    ax.set_ylabel("patient ΔP")
    ax.set_title(f"{interaction}  within-patient {split}  (patient is the unit)")
    ax.tick_params(axis="x", rotation=25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str) -> None:
    show = tbl.head(18).copy()
    fig, ax = plt.subplots(figsize=(8.8, 0.42 * max(len(show), 1) + 1.4))
    if show.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No differential pairs at the stated gates", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    y = np.arange(len(show))
    colors = np.where(show["mean_delta"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["mean_delta"], color=colors, alpha=0.85)
    ax.axvline(0, color="#444", lw=0.7)
    ax.set_yticks(
        y,
        [f"{r.interaction_name}  [{r.lr_class}]  n={int(r.n_patients)}  p={fmt_p(r.p_meta)}" for r in show.itertuples()],
        fontsize=8,
    )
    ax.invert_yaxis()
    ax.set_xlabel("meta mean patient ΔP")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(coverage, meta_q4, meta_med, meta_between, skip_notes, out: Path) -> None:
    n_locked = int(coverage["patient_id"].nunique()) if not coverage.empty else 0
    n_q4 = int(coverage["eligible_q4q1"].sum()) if "eligible_q4q1" in coverage else 0
    n_med = int(coverage["eligible_median"].sum()) if "eligible_median" in coverage else 0
    units_q4 = sorted(coverage.loc[coverage.get("eligible_q4q1", 0) == 1, "unit"].unique()) if not coverage.empty and "eligible_q4q1" in coverage else []

    def top_md(df: pd.DataFrame, n=12) -> str:
        if df is None or df.empty:
            return "_No pairs passed the detect gate with honest n._\n"
        lines = [
            "| pair | class | k | n_patients | mean ΔP | p_meta | p_Wilcoxon | I² | units |",
            "|---|---|---:|---:|---:|---|---|---:|---|",
        ]
        for r in df.head(n).itertuples():
            lines.append(
                f"| {r.interaction_name} | {r.lr_class} | {int(r.k_units)} | {int(r.n_patients)} | "
                f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_meta)} | {fmt_p(r.p_wilcoxon_patients)} | "
                f"{fmt_num(r.i2, 0) if np.isfinite(r.i2) else 'NA'} | {r.units} |"
            )
        return "\n".join(lines) + "\n"

    def top_between(df: pd.DataFrame, n=10) -> str:
        if df is None or df.empty:
            return "_No between-patient Q4 vs Q1 pairs._\n"
        lines = [
            "| pair | class | k | n_Q1/n_Q4 | mean ΔP | p_meta | I² | units |",
            "|---|---|---:|---|---:|---|---:|---|",
        ]
        for r in df.head(n).itertuples():
            lines.append(
                f"| {r.interaction_name} | {r.lr_class} | {int(r.k_units)} | {int(r.n_q1)}/{int(r.n_q4)} | "
                f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_meta)} | "
                f"{fmt_num(r.i2, 0) if np.isfinite(r.i2) else 'NA'} | {r.units} |"
            )
        return "\n".join(lines) + "\n"

    cov_lines = [
        "| unit | locked n | Q4-eligible | median-eligible | note |",
        "|---|---:|---:|---:|---|",
    ]
    for unit in UNITS:
        sub = coverage[coverage["unit"] == unit] if not coverage.empty else pd.DataFrame()
        locked_n = int(sub["patient_id"].nunique()) if not sub.empty else (9 if unit == "GSE253013" else 0)
        q4 = int(sub["eligible_q4q1"].sum()) if not sub.empty and "eligible_q4q1" in sub else 0
        med = int(sub["eligible_median"].sum()) if not sub.empty and "eligible_median" in sub else 0
        note = skip_notes.get(unit, "")
        cov_lines.append(f"| {unit} | {locked_n} | {q4} | {med} | {note} |")

    text = f"""# FINDING — CLDN4-only high-end CellChat on the strict-malignant 6-unit combo

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit.

The strict-malignant 6-unit mean-vs-T/NK Spearman is **taken as given** and is
**not re-audited**: n={GIVEN_N}, ρ={GIVEN_RHO}, p={GIVEN_P}, I²={GIVEN_I2:g}%
(PR #312 / `methods/strict_malig_tnk`). Units:
GSE207422 (DRMref, 12) + GSE205335 (author mal, 22) + GSE291670 (marker, 6)
+ GSE253013 (marker, 9) + GSE131907 (author mal n_mal≥20, 21) + GSE325414
(author mal, 25).

This folder adds CellChat-style **outgoing malignant → T/NK** on the
**high-end** (within-patient CLDN4 Q4 vs Q1; median split as sensitivity)
and meta-analyzes **patient ΔP** across units. Jin et al. 2021 Hill
probability on CellChatDB v2 protein pairs. CellChat R was not run.

GSE253013’s only public processed matrix is a **9.3 GB RDS** and was
**not downloaded**. That unit contributes locked n=9 to the given Spearman
row and **0** patients to the LR meta.

## Honest n

Locked combo n={GIVEN_N} (given). CellChat-eligible patients are fewer
because both CLDN4-high and CLDN4-low malignant arms plus T/NK must meet
the cell floor (Q4 vs Q1: n_mal≥{MIN_MAL_Q4} and ≥{MIN_CELLS_ARM}/arm;
median: ≥{2 * MIN_CELLS_ARM} malignant and ≥{MIN_CELLS_ARM} T/NK).

Within-patient Q4 vs Q1 eligible: **n={n_q4}**. Median-split eligible: **n={n_med}**.
Units with ≥1 Q4-eligible patient: {', '.join(units_q4) if units_q4 else 'none'}.

{chr(10).join(cov_lines)}

GSE207422 DRMref barcodes are not public. CellChat on that unit uses
epithelial lineage as the malignant proxy on the same 12 locked samples.

## Primary: within-patient Q4 vs Q1 ΔP (outgoing Mal → T/NK)

ΔP = P(CLDN4-high mal → same-patient T/NK) − P(CLDN4-low mal → T/NK).
Each patient is one delta. Unit means are random-effects pooled
(DerSimonian–Laird). p-values are descriptive.

{top_md(meta_q4)}

Full table: `results/lr_meta_q4q1.tsv`.

## Sensitivity: within-patient median split

{top_md(meta_med)}

Full table: `results/lr_meta_median.tsv`.

## Extra: between-patient high-end (given CLDN4 Q4 vs Q1)

Quartiles use the **given** malignant CLDN4 scores from the locked tables
(the same vectors as the n=95 ρ=−0.260 row). Per-patient P is scored from
all malignant cells → same-patient T/NK. This is not a re-audit of the
Spearman.

{top_between(meta_between)}

Full table: `results/lr_meta_between_q4q1.tsv`.

## Extra figures

- `figures/fig_honest_n.png` — locked vs CellChat-eligible n
- `figures/fig_extra_ligand_table.png` — top within-patient Q4 ΔP pairs
- `figures/fig_forest_q4q1.png` — meta forest
- `figures/fig_patient_delta_top.png` — patient ΔP by unit for the top pair
- `figures/fig_forest_between_q4q1.png` — between-patient high-end extra

## What is not claimed

- The n=95 ρ=−0.260 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE253013 LR is empty because the 9.3 GB RDS was not downloaded.
- Cell-pooled permutations (sibling CellChat PRs) are not the test.
- CellChat R visualizations were not generated.

## Reproduce

```bash
python3 methods/sixunit95_cldn4_hiend/scripts/download.py
python3 methods/sixunit95_cldn4_hiend/scripts/analyze.py
```
"""
    (out / "FINDING.md").write_text(text)
    (ROOT / "FINDING.md").write_text(text)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--raw", type=Path, default=Path("/tmp/sixunit95_raw"))
    p.add_argument("--out", type=Path, default=ROOT)
    args = p.parse_args()
    out = args.out
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    locked = load_locked_patients()
    print("locked n", {u: len(locked[u]) for u in UNITS}, "sum", sum(len(locked[u]) for u in UNITS), flush=True)

    # LR pairs: keep those present in the union of downloadable matrices.
    # Load DB against a large human symbol set from the interaction file itself
    # and filter per unit.
    inter = pd.read_csv(DB / "interaction_cellchatdb_v2_protein.csv")
    db_genes: set[str] = set(EXTRA) | set(NORMAL_LUNG)
    for vs in LINEAGE_MARKERS.values():
        db_genes.update(vs)
    for rec in inter.itertuples(index=False):
        db_genes.update(parse_symbols(getattr(rec, "ligand.symbol", None) or rec.ligand))
        db_genes.update(parse_symbols(getattr(rec, "receptor.symbol", None) or rec.receptor))
    wanted = db_genes

    loaders = {
        "GSE207422": load_gse207422,
        "GSE205335": load_gse205335,
        "GSE291670": load_gse291670,
        "GSE131907": load_gse131907,
        "GSE325414": load_gse325414,
    }
    skip_notes = {
        "GSE253013": "9.3 GB RDS not downloaded; LR n=0",
    }
    coverage_all = []
    pairs_all = []
    for unit, loader in loaders.items():
        bundle = loader(args.raw, wanted)
        skip_notes.setdefault(unit, bundle.get("note", ""))
        lr = load_lr(DB, bundle["genes"])
        print(f"  {unit} LR pairs in matrix: {len(lr)}", flush=True)
        cov, pairs = score_bundle(unit, bundle, lr, modes=["q4q1", "median", "tertile"])
        coverage_all.append(cov)
        pairs_all.append(pairs)
        del bundle

    coverage = pd.concat(coverage_all, ignore_index=True) if coverage_all else pd.DataFrame()
    # locked-only row for GSE253013
    miss = locked["GSE253013"][["unit", "patient_id"]].copy()
    miss["n_mal"] = np.nan
    miss["n_tnk"] = np.nan
    miss["mean_cldn4_mal"] = np.nan
    miss["eligible_allmal"] = 0
    miss["eligible_q4q1"] = 0
    miss["eligible_median"] = 0
    miss["eligible_tertile"] = 0
    coverage = pd.concat([coverage, miss], ignore_index=True)

    pairs = pd.concat(pairs_all, ignore_index=True) if pairs_all else pd.DataFrame()
    coverage.to_csv(out / "results" / "patient_coverage.tsv", sep="\t", index=False)
    if not pairs.empty:
        pairs.to_csv(out / "results" / "patient_lr_long.tsv.gz", sep="\t", index=False, compression="gzip")

    meta_q4 = meta_patient_deltas(pairs, "q4q1") if not pairs.empty else pd.DataFrame()
    meta_med = meta_patient_deltas(pairs, "median") if not pairs.empty else pd.DataFrame()
    meta_ter = meta_patient_deltas(pairs, "tertile") if not pairs.empty else pd.DataFrame()
    meta_between = meta_q4q1_between(pairs, locked) if not pairs.empty else pd.DataFrame()

    if not meta_q4.empty:
        meta_q4.to_csv(out / "results" / "lr_meta_q4q1.tsv", sep="\t", index=False)
    if not meta_med.empty:
        meta_med.to_csv(out / "results" / "lr_meta_median.tsv", sep="\t", index=False)
    if not meta_ter.empty:
        meta_ter.to_csv(out / "results" / "lr_meta_tertile.tsv", sep="\t", index=False)
    if not meta_between.empty:
        meta_between.to_csv(out / "results" / "lr_meta_between_q4q1.tsv", sep="\t", index=False)

    # unit-level delta table for the primary split
    if not pairs.empty:
        unit_rows = []
        sub = pairs[(pairs["split"] == "q4q1") & (pairs["arm"] == "delta") & (pairs["detected_either"])]
        for (unit, name), g in sub.groupby(["unit", "interaction_name"]):
            st = unit_delta_stats(g["delta"])
            first = g.iloc[0]
            unit_rows.append(
                {
                    "unit": unit,
                    "interaction_name": name,
                    "pathway_name": first["pathway_name"],
                    "lr_class": first["lr_class"],
                    **st,
                }
            )
        pd.DataFrame(unit_rows).to_csv(out / "results" / "lr_unit_q4q1.tsv", sep="\t", index=False)

    plot_n(coverage, out / "figures" / "fig_honest_n.png")
    plot_forest(meta_q4, out / "figures" / "fig_forest_q4q1.png", "Within-patient Q4 vs Q1 outgoing ΔP (RE meta)")
    plot_ligand_table(meta_q4, out / "figures" / "fig_extra_ligand_table.png", "Extra: top outgoing patient-ΔP pairs (Q4 vs Q1)")
    plot_forest(meta_between, out / "figures" / "fig_forest_between_q4q1.png", "Between-patient given-CLDN4 Q4 vs Q1 (extra)")
    if not meta_q4.empty:
        plot_patient_deltas(pairs, meta_q4.iloc[0]["interaction_name"], "q4q1", out / "figures" / "fig_patient_delta_top.png")
    if not meta_med.empty:
        plot_forest(meta_med, out / "figures" / "fig_forest_median.png", "Within-patient median-split outgoing ΔP (sensitivity)")

    summary = {
        "given_spearman": {"n": GIVEN_N, "rho": GIVEN_RHO, "p": GIVEN_P, "i2": GIVEN_I2, "re_audited": False},
        "locked_n": {u: int(len(locked[u])) for u in UNITS},
        "eligible_q4q1": int(coverage["eligible_q4q1"].sum()) if "eligible_q4q1" in coverage else 0,
        "eligible_median": int(coverage["eligible_median"].sum()) if "eligible_median" in coverage else 0,
        "n_meta_q4_pairs": int(len(meta_q4)),
        "n_meta_between_pairs": int(len(meta_between)),
        "skipped": skip_notes,
        "dual_high": False,
    }
    (out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(coverage, meta_q4, meta_med, meta_between, skip_notes, out / "results")
    print("done", json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
