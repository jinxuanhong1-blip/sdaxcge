#!/usr/bin/env python3
"""CLDN4-only ICI-adjacent scRNA trio: GSE207422 + GSE205335 (+ GSE179994 if usable).

Patient is the unit. No TACSTD2 gate. No dual-high.
Primary endpoints:
  1. Patient-level malignant CLDN4 (mean log1p CP10k) vs same-patient T/NK fraction.
  2. CellPhoneDB-style outgoing score from CLDN4-high vs CLDN4-low malignant cells → T/NK.

GSE179994 public processed expression is T-cell-only; it cannot enter either endpoint.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

from lib_stats import random_effects_dl, spearman, stouffer, wilcoxon_safe

HERE = Path(__file__).resolve().parent
MIN_MALIG = 10
MIN_TNK = 20
MIN_BIN = 10
EXPR_PROP = 0.10

PAPER_GROUP_207422 = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
    "mast": ["KIT"],
}
NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
NSCLC = {"ADC", "SQ"}


def log(msg: str) -> None:
    print(msg, flush=True)


def load_pairs(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def pair_genes(pairs: pd.DataFrame) -> set[str]:
    genes: set[str] = set()
    for col in ("ligand", "receptor"):
        for item in pairs[col].astype(str):
            genes.update(item.split("+"))
    return genes


def wanted_207422(pairs: pd.DataFrame) -> set[str]:
    keep = set(pair_genes(pairs))
    keep.update(["CLDN4", "TACSTD2", "PTPRC", "CD8A", "CD8B", "CXCL13"])
    for genes in LINEAGES.values():
        keep.update(genes)
    keep.update(NORMAL_LUNG)
    return keep


def wanted_205335(pairs: pd.DataFrame) -> set[str]:
    keep = set(pair_genes(pairs))
    keep.update(["CLDN4", "TACSTD2", "PTPRC", "CD8A", "CD8B", "CXCL13"])
    return keep


def score_lineage(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score_lineage(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def stream_207422(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    log(f"[GSE207422] stream {path}")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, rest = line.split("\t", 1)
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch for {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep:
                store[gene] = vals
    log(f"[GSE207422] cells={n} genes_in_file={n_genes} genes_kept={len(store)}")
    return cells, store, totals


def process_207422(datadir: Path, pairs: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    matrix = datadir / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    cells, expr, totals = stream_207422(matrix, wanted_207422(pairs))
    n = len(cells)
    sample = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    epi = lineage == "epithelial"
    normal_umi = np.zeros(n, dtype=np.float64)
    for g in NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]
    malig = epi & (normal_umi == 0)
    tnk = np.isin(lineage, ["T", "NK"])
    cldn4 = log1p_cp10k(expr.get("CLDN4", np.zeros(n)), totals)
    cd8 = np.zeros(n, dtype=bool)
    if "CD8A" in expr and "CD8B" in expr:
        cd8 = tnk & ((expr["CD8A"] > 0) | (expr["CD8B"] > 0))
    cxcl13 = np.zeros(n, dtype=bool)
    if "CXCL13" in expr:
        cxcl13 = tnk & (expr["CXCL13"] > 0)

    rows = []
    for s in sorted(set(sample)):
        m = sample == s
        n_mal = int(malig[m].sum())
        n_tnk = int(tnk[m].sum())
        n_all = int(m.sum())
        rows.append(
            {
                "cohort": "GSE207422",
                "patient": s.replace("BD_immune", "P"),
                "sample": s,
                "paper_group": PAPER_GROUP_207422.get(s, "NA"),
                "cancer_subtype": "NSCLC",
                "response": PAPER_GROUP_207422.get(s, "NA"),
                "ici_adjacent": PAPER_GROUP_207422.get(s) in {"MPR", "NMPR"},
                "n_cells": n_all,
                "n_malignant": n_mal,
                "n_tnk": n_tnk,
                "n_cd8": int(cd8[m].sum()),
                "n_cxcl13pos_tnk": int(cxcl13[m].sum()),
                "frac_tnk": n_tnk / n_all if n_all else np.nan,
                "frac_cd8": int(cd8[m].sum()) / n_all if n_all else np.nan,
                "malig_CLDN4_mean": float(cldn4[m & malig].mean()) if n_mal else np.nan,
                "malig_CLDN4_pct": float((expr.get("CLDN4", np.zeros(n))[m & malig] > 0).mean()) if n_mal else np.nan,
                "eligible_combo": int(n_mal >= MIN_MALIG and n_tnk >= MIN_TNK),
                "malig_def": "A3_marker_malignant_like",
                "immune_def": "marker_T_or_NK",
            }
        )
    patients = pd.DataFrame(rows)
    payload = {
        "cohort": "GSE207422",
        "cells": cells,
        "sample": sample,
        "expr": expr,
        "totals": totals,
        "malig": malig,
        "tnk": tnk,
        "cldn4": cldn4,
        "lineage": lineage,
    }
    return patients, payload


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
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
    read_end = metadata.get("platform", pd.Series([""] * len(metadata))).astype(str).str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].astype(str).str.replace("_", "-", regex=False) + "-" + read_end.fillna("") + "P"
    )
    return metadata.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


def load_205335_matrix(path: Path, keep: set[str]) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        # GEO supplement is often double-gzipped.
        if path.suffix == ".gz":
            inner = Path(tmp) / path.stem
            with gzip.open(path, "rb") as source, inner.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
            matrix_path = inner
            # If still gzip-magic, unwrap once more.
            with matrix_path.open("rb") as fh:
                magic = fh.read(2)
            if magic == b"\x1f\x8b":
                inner2 = Path(tmp) / "matrix.rds"
                with gzip.open(matrix_path, "rb") as source, inner2.open("wb") as dest:
                    shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
                matrix_path = inner2
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(keep):
        positions = np.flatnonzero(genes == gene)
        if len(positions) == 0:
            continue
        extracted[gene] = np.asarray(matrix.getrow(int(positions[0])).toarray()).ravel()
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    log(f"[GSE205335] cells={len(barcodes)} genes_kept={len(extracted)}")
    return extracted, library_umi, barcodes


def process_205335(datadir: Path, pairs: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ident_path = datadir / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    matrix_path = datadir / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    soft_path = datadir / "GSE205335" / "GSE205335_family.soft.gz"
    identities = pd.read_csv(ident_path, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("Cell-identity barcodes are not unique")
    expr, totals, barcodes = load_205335_matrix(matrix_path, wanted_205335(pairs))
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(f"Matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
    cells = indexed.loc[barcodes].reset_index()
    meta_csv = HERE / "resources" / "gse205335_sample_metadata.csv"
    if meta_csv.exists():
        meta = pd.read_csv(meta_csv)
    elif soft_path.exists() and soft_path.stat().st_size > 50_000:
        meta = parse_geo_soft(soft_path)
    else:
        raise FileNotFoundError("GSE205335 sample metadata CSV missing and family SOFT is a stub")
    keep_cols = [c for c in ["orig.ident", "gsm", "patient", "tissue", "recist", "platform", "cancer_subtype"] if c in meta.columns]
    cells = cells.merge(meta[keep_cols], on="orig.ident", how="left", validate="many_to_one")
    if cells["patient"].isna().any():
        n_miss = int(cells["patient"].isna().sum())
        raise ValueError(f"{n_miss} cells did not match GSE205335 sample metadata")
    n = len(cells)
    tissue = cells["tissue"].astype(str) if "tissue" in cells.columns else pd.Series(["NA"] * n)
    is_normal = tissue.str.contains("Normal", case=False, na=False).to_numpy()
    tumor = ~is_normal
    malig = cells["lineage.sub"].eq("Malignant cells").to_numpy() & tumor
    tnk = cells["lineage.total"].eq("T/NK cells").to_numpy() & tumor
    cldn4 = log1p_cp10k(expr.get("CLDN4", np.zeros(n)), totals)
    cd8 = np.zeros(n, dtype=bool)
    if "lineage.sub" in cells.columns:
        cd8 = cells["lineage.sub"].astype(str).str.contains("CD8", case=False, na=False).to_numpy()
    if not cd8.any() and "CD8A" in expr:
        cd8 = tnk & (expr["CD8A"] > 0)
    cxcl13 = np.zeros(n, dtype=bool)
    if "CXCL13" in expr:
        cxcl13 = tnk & (expr["CXCL13"] > 0)
    sample = cells["orig.ident"].astype(str).to_numpy()
    patient = cells["patient"].astype(str).to_numpy() if "patient" in cells.columns else sample
    subtype = cells["cancer_subtype"].astype(str).to_numpy() if "cancer_subtype" in cells.columns else np.array(["NA"] * n)
    recist = cells["recist"].astype(str).to_numpy() if "recist" in cells.columns else np.array(["NA"] * n)

    rows = []
    for p in pd.unique(patient):
        m = (patient == p) & tumor
        if not m.any():
            continue
        n_mal = int(malig[m].sum())
        n_tnk = int(tnk[m].sum())
        n_all = int(m.sum())
        sub = subtype[m][0] if n_all else "NA"
        rec = recist[m][0] if n_all else "NA"
        rows.append(
            {
                "cohort": "GSE205335",
                "patient": p,
                "sample": ",".join(sorted(set(sample[m]))),
                "paper_group": RESPONSE_MAP.get(rec, rec),
                "cancer_subtype": sub,
                "response": RESPONSE_MAP.get(rec, rec),
                "ici_adjacent": True,
                "n_cells": n_all,
                "n_malignant": n_mal,
                "n_tnk": n_tnk,
                "n_cd8": int(cd8[m].sum()),
                "n_cxcl13pos_tnk": int(cxcl13[m].sum()),
                "frac_tnk": n_tnk / n_all if n_all else np.nan,
                "frac_cd8": int(cd8[m].sum()) / n_all if n_all else np.nan,
                "malig_CLDN4_mean": float(cldn4[m & malig].mean()) if n_mal else np.nan,
                "malig_CLDN4_pct": float((expr.get("CLDN4", np.zeros(n))[m & malig] > 0).mean()) if n_mal else np.nan,
                "eligible_combo": int(n_mal >= MIN_MALIG and n_tnk >= MIN_TNK),
                "malig_def": "author_malignant",
                "immune_def": "author_TNK",
            }
        )
    patients = pd.DataFrame(rows)
    payload = {
        "cohort": "GSE205335",
        "cells": barcodes,
        "sample": patient,
        "expr": expr,
        "totals": totals,
        "malig": malig,
        "tnk": tnk,
        "cldn4": cldn4,
    }
    return patients, payload


def audit_179994(datadir: Path) -> dict:
    meta_path = datadir / "GSE179994" / "GSE179994_Tcell.metadata.tsv.gz"
    out = {
        "cohort": "GSE179994",
        "usable": False,
        "reason": (
            "Public processed expression is T-cell-only "
            "(GSE179994_all.Tcell.rawCounts.rds.gz, 421 MB). "
            "No all-cell UMI/h5/TISCH extract <2 GB. "
            "Malignant CLDN4 and a same-patient T/NK denominator are not estimable."
        ),
        "analysis_n_patients": 0,
        "tcell_rds_mb": 421.0,
        "tcell_rds_under_2gb": True,
        "tcell_only": True,
    }
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t")
        out["n_tcell_barcodes"] = int(len(meta))
        for col in meta.columns:
            if col.lower() in {"patient", "patient_id", "donor", "orig.ident"}:
                out["n_tcell_patients"] = int(meta[col].nunique())
            if "celltype" in col.lower() or col.lower() == "celltype":
                out["author_celltypes"] = meta[col].astype(str).value_counts().to_dict()
        if "celltype" in meta.columns:
            out["author_celltypes"] = meta["celltype"].astype(str).value_counts().to_dict()
        # patient-like columns
        for col in meta.columns:
            if "patient" in col.lower() or col.lower() in {"orig.ident", "sample", "donor"}:
                out[f"n_unique_{col}"] = int(meta[col].nunique())
    return out


def combo_rows(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = [
        ("GSE207422_postICI", (patients["cohort"].eq("GSE207422") & patients["ici_adjacent"] & patients["eligible_combo"].eq(1))),
        ("GSE207422_all_eligible", (patients["cohort"].eq("GSE207422") & patients["eligible_combo"].eq(1))),
        ("GSE205335_all", (patients["cohort"].eq("GSE205335") & patients["eligible_combo"].eq(1))),
        ("GSE205335_NSCLC", (patients["cohort"].eq("GSE205335") & patients["cancer_subtype"].isin(NSCLC) & patients["eligible_combo"].eq(1))),
        ("GSE179994", patients["cohort"].eq("GSE179994")),
    ]
    stats_by = {}
    for name, mask in subsets:
        sub = patients.loc[mask]
        rho, p, n = spearman(sub["malig_CLDN4_mean"], sub["frac_tnk"])
        rho_pct, p_pct, n_pct = spearman(sub["malig_CLDN4_pct"], sub["frac_tnk"])
        rho_cd8, p_cd8, n_cd8 = spearman(sub["malig_CLDN4_mean"], sub["frac_cd8"])
        stats_by[name] = {
            "subset": name,
            "k_cohorts": 1 if n else 0,
            "n_patients": n,
            "malig_score": "mean_log1p_cp10k",
            "immune": "frac_tnk",
            "rho": rho,
            "p": p,
            "rho_pctpos_vs_tnk": rho_pct,
            "p_pctpos_vs_tnk": p_pct,
            "rho_mean_vs_cd8": rho_cd8,
            "p_mean_vs_cd8": p_cd8,
            "I2": 0.0 if n else np.nan,
            "note": "single-cohort Spearman; patient is the unit",
        }
        rows.append(stats_by[name])

    def add_pool(name: str, members: list[str], note: str) -> None:
        rhos = [stats_by[m]["rho"] for m in members]
        ps = [stats_by[m]["p"] for m in members]
        ns = [stats_by[m]["n_patients"] for m in members]
        if any(n < 3 for n in ns):
            rows.append(
                {
                    "subset": name,
                    "k_cohorts": sum(n >= 3 for n in ns),
                    "n_patients": int(sum(ns)),
                    "malig_score": "mean_log1p_cp10k",
                    "immune": "frac_tnk",
                    "rho": np.nan,
                    "p": np.nan,
                    "rho_pctpos_vs_tnk": np.nan,
                    "p_pctpos_vs_tnk": np.nan,
                    "rho_mean_vs_cd8": np.nan,
                    "p_mean_vs_cd8": np.nan,
                    "I2": np.nan,
                    "note": "not pooled: a member has n<3",
                }
            )
            return
        re = random_effects_dl(rhos, ns)
        st = stouffer(rhos, ps, ns)
        rows.append(
            {
                "subset": name,
                "k_cohorts": re.get("k", 0),
                "n_patients": re.get("n_patients_total", int(sum(ns))),
                "malig_score": "mean_log1p_cp10k",
                "immune": "frac_tnk",
                "rho": re.get("pooled_rho", np.nan),
                "p": re.get("p", np.nan),
                "rho_pctpos_vs_tnk": random_effects_dl(
                    [stats_by[m]["rho_pctpos_vs_tnk"] for m in members],
                    [stats_by[m]["n_patients"] for m in members],
                ).get("pooled_rho", np.nan),
                "p_pctpos_vs_tnk": random_effects_dl(
                    [stats_by[m]["rho_pctpos_vs_tnk"] for m in members],
                    [stats_by[m]["n_patients"] for m in members],
                ).get("p", np.nan),
                "stouffer_z": st.get("z", np.nan),
                "stouffer_p": st.get("p", np.nan),
                "I2": re.get("I2", np.nan),
                "ci95_lo": re.get("ci95_rho", [np.nan, np.nan])[0],
                "ci95_hi": re.get("ci95_rho", [np.nan, np.nan])[1],
                "note": note,
            }
        )

    add_pool(
        "combo_207422post_205335all",
        ["GSE207422_postICI", "GSE205335_all"],
        "Primary ICI duo (179994 unusable). DL RE on Fisher-z. Includes 207422 post-ICI.",
    )
    add_pool(
        "combo_207422post_205335NSCLC",
        ["GSE207422_postICI", "GSE205335_NSCLC"],
        "Sensitivity: drop SCLC/NUT from 205335. DL RE on Fisher-z.",
    )
    add_pool(
        "combo_207422all_205335all",
        ["GSE207422_all_eligible", "GSE205335_all"],
        "Sensitivity: 207422 includes treatment-naive biopsies. DL RE on Fisher-z.",
    )
    return pd.DataFrame(rows)


def group_gene_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict[str, float], fracs: dict[str, float]) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(pairs: pd.DataFrame, logx: dict[str, np.ndarray], sender: np.ndarray, receiver: np.ndarray) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "pass_expr_prop": (l_frac >= EXPR_PROP) and (rec_f >= EXPR_PROP),
            }
        )
    return pd.DataFrame(rows)


def build_logx(expr: dict[str, np.ndarray], totals: np.ndarray, genes: set[str]) -> dict[str, np.ndarray]:
    logx = {}
    n = len(totals)
    for g in genes:
        if g in expr:
            logx[g] = log1p_cp10k(expr[g], totals)
        else:
            logx[g] = np.zeros(n, dtype=np.float32)
    return logx


def lr_outgoing(payload: dict, pairs: pd.DataFrame, eligible_patients: set[str] | None = None) -> pd.DataFrame:
    genes = pair_genes(pairs)
    logx = build_logx(payload["expr"], payload["totals"], genes)
    sample = payload["sample"]
    malig = payload["malig"]
    tnk = payload["tnk"]
    cldn4 = payload["cldn4"]
    rows = []
    for s in sorted(set(sample)):
        if eligible_patients is not None and s not in eligible_patients and s.replace("BD_immune", "P") not in eligible_patients:
            # still allow raw sample ids
            pass
        m = sample == s
        mal = m & malig
        rec = m & tnk
        if int(mal.sum()) < 2 * MIN_BIN or int(rec.sum()) < MIN_TNK:
            continue
        vals = cldn4[mal]
        thr = float(np.median(vals))
        high = mal & (cldn4 >= thr)
        low = mal & (cldn4 < thr)
        # if median ties leave one side short, use strict complementary split
        if int(high.sum()) < MIN_BIN or int(low.sum()) < MIN_BIN:
            order = np.where(mal)[0]
            ranked = order[np.argsort(cldn4[order])]
            half = len(ranked) // 2
            low_idx = ranked[:half]
            high_idx = ranked[half:]
            high = np.zeros(len(malig), dtype=bool)
            low = np.zeros(len(malig), dtype=bool)
            high[high_idx] = True
            low[low_idx] = True
        if int(high.sum()) < MIN_BIN or int(low.sum()) < MIN_BIN:
            continue
        hi = score_pairs(pairs, logx, high, rec)
        lo = score_pairs(pairs, logx, low, rec)
        if hi.empty or lo.empty:
            continue
        key = ["ligand", "receptor", "pathway"]
        merged = hi.merge(lo, on=key, suffixes=("_high", "_low"))
        merged["delta_high_minus_low"] = merged["cpdb_mean_score_high"] - merged["cpdb_mean_score_low"]
        merged["cohort"] = payload["cohort"]
        merged["patient"] = s.replace("BD_immune", "P") if payload["cohort"] == "GSE207422" else s
        merged["sample"] = s
        merged["n_malig_high"] = int(high.sum())
        merged["n_malig_low"] = int(low.sum())
        merged["n_tnk"] = int(rec.sum())
        merged["cldn4_threshold"] = thr
        rows.append(merged)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def summarize_lr(patient_lr: pd.DataFrame) -> pd.DataFrame:
    if patient_lr.empty:
        return pd.DataFrame()
    rows = []
    for (cohort, lig, rec, path), g in patient_lr.groupby(["cohort", "ligand", "receptor", "pathway"], sort=False):
        deltas = g["delta_high_minus_low"].to_numpy(dtype=float)
        highs = g["cpdb_mean_score_high"].to_numpy(dtype=float)
        lows = g["cpdb_mean_score_low"].to_numpy(dtype=float)
        n = int(len(deltas))
        p = wilcoxon_safe(highs, lows)
        pass_n = int((g["pass_expr_prop_high"] | g["pass_expr_prop_low"]).sum())
        rows.append(
            {
                "subset": cohort,
                "ligand": lig,
                "receptor": rec,
                "pathway": path,
                "n_patients": n,
                "n_pass_expr_prop": pass_n,
                "median_delta": float(np.median(deltas)),
                "mean_delta": float(np.mean(deltas)),
                "n_delta_neg": int((deltas < 0).sum()),
                "n_delta_pos": int((deltas > 0).sum()),
                "wilcoxon_p": p,
            }
        )
    # combined across usable cohorts
    for (lig, rec, path), g in patient_lr.groupby(["ligand", "receptor", "pathway"], sort=False):
        deltas = g["delta_high_minus_low"].to_numpy(dtype=float)
        highs = g["cpdb_mean_score_high"].to_numpy(dtype=float)
        lows = g["cpdb_mean_score_low"].to_numpy(dtype=float)
        n = int(len(deltas))
        p = wilcoxon_safe(highs, lows)
        rows.append(
            {
                "subset": "combo_" + "+".join(sorted(g["cohort"].unique())),
                "ligand": lig,
                "receptor": rec,
                "pathway": path,
                "n_patients": n,
                "n_pass_expr_prop": int((g["pass_expr_prop_high"] | g["pass_expr_prop_low"]).sum()),
                "median_delta": float(np.median(deltas)),
                "mean_delta": float(np.mean(deltas)),
                "n_delta_neg": int((deltas < 0).sum()),
                "n_delta_pos": int((deltas > 0).sum()),
                "wilcoxon_p": p,
            }
        )
    out = pd.DataFrame(rows)
    for subset, idx in out.groupby("subset").groups.items():
        pvals = out.loc[idx, "wilcoxon_p"].to_numpy(dtype=float)
        finite = np.isfinite(pvals)
        q = np.full(pvals.shape, np.nan)
        if finite.sum():
            q[finite] = multipletests(pvals[finite], method="fdr_bh")[1]
        out.loc[idx, "fdr"] = q
    return out.sort_values(["subset", "pathway", "median_delta"])


def plot_scatter(patients: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    colors = {"GSE207422": "#1f4e79", "GSE205335": "#c45911"}
    markers = {True: "o", False: "D"}
    for cohort, g in patients.groupby("cohort"):
        el = g[g["eligible_combo"] == 1]
        ax.scatter(
            el["malig_CLDN4_mean"],
            el["frac_tnk"],
            c=colors.get(cohort, "gray"),
            marker=markers.get(True, "o"),
            s=42,
            label=f"{cohort} n={len(el)}",
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title("CLDN4-only · ICI duo (179994 unusable)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_lr(lr_sum: pd.DataFrame, path: Path) -> None:
    combo = lr_sum[lr_sum["subset"].str.startswith("combo_")].copy()
    if combo.empty:
        return
    combo = combo.sort_values("median_delta")
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    y = np.arange(len(combo))
    colors = np.where(combo["median_delta"] < 0, "#7b2d26", "#1f4e79")
    ax.barh(y, combo["median_delta"], color=colors, height=0.7)
    ax.axvline(0, color="black", lw=0.6)
    labels = [f"{r.pathway} {r.ligand}–{r.receptor} (n={int(r.n_patients)})" for r in combo.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Median patient Δ (CLDN4-high − CLDN4-low) outgoing score")
    ax.set_title("Outgoing CLDN4-high malignant → T/NK")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if abs(x) < 0.001 and x != 0:
        return f"{x:.2e}"
    return f"{x:.{nd}f}"


def write_finding(
    path: Path,
    audit: dict,
    patients: pd.DataFrame,
    combo: pd.DataFrame,
    lr_sum: pd.DataFrame,
    lr_patients: pd.DataFrame,
) -> None:
    def row(name: str) -> pd.Series:
        hit = combo[combo["subset"] == name]
        return hit.iloc[0] if len(hit) else pd.Series(dtype=object)

    primary = row("combo_207422post_205335all")
    c422 = row("GSE207422_postICI")
    c335 = row("GSE205335_all")
    c335n = row("GSE205335_NSCLC")
    c422a = row("GSE207422_all_eligible")

    el = patients[patients["eligible_combo"] == 1]
    n422 = int(((el["cohort"] == "GSE207422") & el["ici_adjacent"]).sum())
    n335 = int((el["cohort"] == "GSE205335").sum())
    n335_nsclc = int(((el["cohort"] == "GSE205335") & el["cancer_subtype"].isin(NSCLC)).sum())

    lr_combo = lr_sum[lr_sum["subset"].str.startswith("combo_")].copy() if not lr_sum.empty else pd.DataFrame()
    lr_n = {}
    if not lr_patients.empty:
        for cohort, g in lr_patients.groupby("cohort"):
            lr_n[cohort] = int(g["patient"].nunique())
        lr_n["combo"] = int(lr_patients["patient"].nunique())

    lines = [
        "# ICI trio CLDN4-only: GSE207422 + GSE205335 + GSE179994",
        "",
        "Additive **CLDN4-only** (no dual-high, no TACSTD2 gate). Patient is the unit.",
        "Intended merge: GSE207422 (Hu *Genome Med* 2023, neoadjuvant PD-1 + chemo) +",
        "GSE205335 (palliative ICI biopsies) + GSE179994 (pembro + carboplatin + pemetrexed)",
        "if GSE179994 has a public processed all-cell matrix <2 GB.",
        "",
        "## GSE179994 is unusable",
        "",
        audit["reason"],
        "",
        f"- T-cell RDS exists and is under 2 GB ({audit.get('tcell_rds_mb', 421)} MB) but is **T cells only**.",
        f"- Analysis n (patients with malignant CLDN4 and a same-patient all-cell T/NK denominator) = **{audit.get('analysis_n_patients', 0)}**.",
        "- Do not cite T-cell barcode counts as n.",
        "- The trio therefore runs as the **ICI duo GSE207422 + GSE205335**. 207422 is included.",
        "",
        "## Honest n (patient-level malignant CLDN4 vs T/NK)",
        "",
        f"Gates: ≥{MIN_MALIG} malignant cells and ≥{MIN_TNK} T/NK cells. Score = malignant CLDN4 mean log1p(CP10k) vs same-patient T/NK fraction.",
        "",
        "| Cohort | Malignant def | Immune def | Eligible n | ICI-adjacent n | Note |",
        "|---|---|---|---:|---:|---|",
        f"| GSE207422 | A3 marker malignant-like (epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3) | marker T or NK | {int(c422a.get('n_patients', 0) or 0)} | {n422} | Author per-cell labels are not on GEO. TN biopsies exist; primary combo uses post-ICI only. |",
        f"| GSE205335 | author `Malignant cells` | author `T/NK cells` | {n335} | {n335} | RECIST only; no MPR. NSCLC (ADC+SQ) n={n335_nsclc}. |",
        f"| GSE179994 | none | none | **0** | **0** | T-cell-only public matrix. |",
        f"| **Primary combo** | mixed (as above) | mixed | **{int(primary.get('n_patients', 0) or 0)}** | **{n422 + n335}** | k=2 cohorts. 179994 not in the pool. |",
        "",
        "Cells are not replicates. Q4 vs Q1 is not done (would need n≥16 **per cohort** or a pre-specified pooled quartile; not claimed).",
        "",
        "## Combo table — malignant CLDN4 vs T/NK",
        "",
        "| Subset | k | N patients | ρ | p | I² | Note |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for rec in combo.itertuples():
        lines.append(
            f"| {rec.subset} | {int(rec.k_cohorts) if pd.notna(rec.k_cohorts) else 0} | {int(rec.n_patients) if pd.notna(rec.n_patients) else 0} | "
            f"{fmt(rec.rho)} | {fmt(rec.p)} | {fmt(getattr(rec, 'I2', np.nan), 0)} | {rec.note} |"
        )
    lines += [
        "",
        f"Primary row **combo_207422post_205335all**: ρ={fmt(primary.get('rho'))}, p={fmt(primary.get('p'))}, "
        f"N={int(primary.get('n_patients', 0) or 0)}, I²={fmt(primary.get('I2'), 0)}%. "
        f"GSE207422 post-ICI ρ={fmt(c422.get('rho'))} (n={int(c422.get('n_patients', 0) or 0)}); "
        f"GSE205335 ρ={fmt(c335.get('rho'))} (n={int(c335.get('n_patients', 0) or 0)}).",
        "",
        "p-values are descriptive. This is a two-cohort merge, not a 15-cohort claim.",
        "",
        "## CellChat-style outgoing CLDN4-high malignant → T/NK",
        "",
        "Not CellChat (R / CellChat were not run). Documented CellPhoneDB-style score:",
        "partner expression = min(subunit means) on log1p(CP10k); pair score = mean of the two partner means",
        "(Efremova 2020; Garcia-Alonso 2022). High/low = **within-patient median** of malignant CLDN4.",
        f"A patient enters if both bins have ≥{MIN_BIN} malignant cells and ≥{MIN_TNK} T/NK.",
        "Negative Δ = weaker outgoing score from the CLDN4-high state.",
        "",
        f"Paired n: GSE207422={lr_n.get('GSE207422', 0)}, GSE205335={lr_n.get('GSE205335', 0)}, "
        f"combo={lr_n.get('combo', 0)}. GSE179994 paired n=0.",
        "",
        "### Combined LR table (usable cohorts only)",
        "",
        "| Pathway | Pair | n patients | median Δ | n Δ<0 | Wilcoxon p | FDR |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    if lr_combo.empty:
        lines.append("| — | — | 0 | NA | 0 | NA | NA |")
    else:
        show = lr_combo.sort_values(["pathway", "median_delta"])
        for rec in show.itertuples():
            lines.append(
                f"| {rec.pathway} | {rec.ligand}–{rec.receptor} | {int(rec.n_patients)} | "
                f"{fmt(rec.median_delta)} | {int(rec.n_delta_neg)} | {fmt(rec.wilcoxon_p)} | {fmt(rec.fdr)} |"
            )
    lines += [
        "",
        "Per-cohort LR rows are in `tables/lr_outgoing_cldn4high_tnk.tsv`.",
        "Patient-level pair scores are in `tables/lr_patient_scores.tsv`.",
        "",
        "## What was not done",
        "",
        "- No dual-high / TACSTD2 gate.",
        "- No CellChat communication probability.",
        "- GSE179994 T-cell RDS was not downloaded (filename + author metadata already establish T-cell-only content).",
        "- Author CopyKAT / DRMref labels for GSE207422 are not public; A3 marker malignant-like is used as given.",
        "- GSE205335 has no MPR field; RECIST is not substituted for MPR.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `tables/combo_cldn4_tnk.tsv` | Combo + single-cohort Spearman / DL RE |",
        "| `tables/patient_cldn4_tnk.tsv` | Patient-level malignant CLDN4 and T/NK |",
        "| `tables/lr_outgoing_cldn4high_tnk.tsv` | LR table (per-cohort and combo) |",
        "| `tables/lr_patient_scores.tsv` | Per-patient outgoing scores |",
        "| `tables/gse179994_feasibility.json` | Why 179994 is out |",
        "| `figures/scatter_cldn4_tnk.png` | Patient scatter |",
        "| `figures/lr_outgoing_combo.png` | Combined LR Δ |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 -m pip install -r methods/ici_trio_207422_205335_179994_cldn4/requirements.txt",
        "python3 methods/ici_trio_207422_205335_179994_cldn4/download.py",
        "python3 methods/ici_trio_207422_205335_179994_cldn4/analyze.py",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    out = args.outdir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    pairs = load_pairs(HERE / "resources" / "lr_pairs.tsv")
    audit = audit_179994(args.datadir)
    (out / "tables" / "gse179994_feasibility.json").write_text(json.dumps(audit, indent=2))
    log(json.dumps(audit, indent=2))

    p422, pay422 = process_207422(args.datadir, pairs)
    p335, pay335 = process_205335(args.datadir, pairs)
    patients = pd.concat([p422, p335], ignore_index=True)
    patients.to_csv(out / "tables" / "patient_cldn4_tnk.tsv", sep="\t", index=False)

    combo = combo_rows(patients)
    combo.to_csv(out / "tables" / "combo_cldn4_tnk.tsv", sep="\t", index=False)
    log(combo.to_string(index=False))

    # LR: 207422 post-ICI patients only for the ICI-adjacent contrast; 205335 all eligible.
    post_ids = set(p422.loc[p422["ici_adjacent"] & p422["eligible_combo"].eq(1), "sample"])
    lr422 = lr_outgoing(pay422, pairs)
    if not lr422.empty:
        lr422 = lr422[lr422["sample"].isin(post_ids) | lr422["patient"].isin(set(p422.loc[p422["ici_adjacent"], "patient"]))]
        # keep only post-ICI
        post_patients = set(p422.loc[p422["ici_adjacent"], "patient"])
        lr422 = lr422[lr422["patient"].isin(post_patients)]
    lr335 = lr_outgoing(pay335, pairs)
    lr_patients = pd.concat([x for x in [lr422, lr335] if not x.empty], ignore_index=True) if (not lr422.empty or not lr335.empty) else pd.DataFrame()
    if not lr_patients.empty:
        lr_patients.to_csv(out / "tables" / "lr_patient_scores.tsv", sep="\t", index=False)
    lr_sum = summarize_lr(lr_patients)
    if not lr_sum.empty:
        lr_sum.to_csv(out / "tables" / "lr_outgoing_cldn4high_tnk.tsv", sep="\t", index=False)
        log(lr_sum[lr_sum["subset"].str.startswith("combo_")].to_string(index=False))

    plot_scatter(patients, out / "figures" / "scatter_cldn4_tnk")
    if not lr_sum.empty:
        plot_lr(lr_sum, out / "figures" / "lr_outgoing_combo")

    write_finding(out / "tables" / "FINDING_auto.md", audit, patients, combo, lr_sum, lr_patients)
    summary = {
        "gse179994_usable": False,
        "primary_combo_n": int(combo.loc[combo["subset"] == "combo_207422post_205335all", "n_patients"].iloc[0])
        if (combo["subset"] == "combo_207422post_205335all").any()
        else 0,
        "lr_paired_n": int(lr_patients["patient"].nunique()) if not lr_patients.empty else 0,
        "min_malig": MIN_MALIG,
        "min_tnk": MIN_TNK,
        "dual_high": False,
    }
    (out / "tables" / "summary.json").write_text(json.dumps(summary, indent=2))
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
