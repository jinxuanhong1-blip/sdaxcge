#!/usr/bin/env python3
"""ADDITIVE CLDN4-only pairwise combo: GSE123902 + GSE131907.

Patient-level malignant CLDN4 vs same-patient T/NK. No dual-high.
Not a pile: only these two public processed matrices (<2 GB).
If the two cohort rhos differ, run CellChat-style LR (Jin 2021 Hill).
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import trim_mean, wilcoxon

from lib_stats import (
    fisher_combine,
    random_effects_dl,
    spearman,
    spearman_ci,
    stouffer,
)

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("PAIR_123902_131907_DATA", "/tmp/pair_123902_131907"))
OUT = HERE / "results"
FIG = HERE / "figures"

MIN_MAL = 20
MIN_TNK = 20
MIN_N_SPEARMAN = 5
MIN_GROUP = 25
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
NBOOT = 100
SEED = 1

# CLDN4 is never used to assign compartments.
EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
MYE_MARKERS = ["LYZ", "CD14", "CSF1R", "AIF1"]
B_MARKERS = ["MS4A1", "CD79A"]
AUDIT_GENES = ["CLDN4", "TACSTD2", "PTPRC", "PECAM1", "COL1A1"]

TUMOR_ORIGINS_131907 = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES_131907 = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES_131907 = {"T lymphocytes", "NK cells"}

FOCUSED_PAIRS = {
    "MHC-I": {"HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"},
    "T-recruit": {"CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL4", "CCL5"},
    "checkpoint": {"CD274", "NECTIN2", "PVR", "PDCD1LG2"},
    "barrier": {"CDH1", "F11R", "CLDN4", "EPCAM"},
}


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand.symbol", None)) or parse_symbols(rec.ligand)
        recp = parse_symbols(getattr(rec, "receptor.symbol", None)) or parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        if any(g not in matrix_genes for g in lig + recp):
            continue
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


def lr_wanted_genes(db_dir: Path) -> set[str]:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    genes: set[str] = set()
    for rec in inter.itertuples(index=False):
        genes.update(parse_symbols(getattr(rec, "ligand.symbol", None)) or parse_symbols(rec.ligand))
        genes.update(parse_symbols(getattr(rec, "receptor.symbol", None)) or parse_symbols(rec.receptor))
    return genes


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def geom_mean_rows(mat: np.ndarray) -> np.ndarray:
    if mat.ndim == 1:
        return mat
    if mat.shape[0] == 1:
        return mat[0]
    out = np.exp(np.mean(np.log(np.clip(mat, 1e-12, None)), axis=0))
    out[np.any(mat <= 0, axis=0)] = 0.0
    return out


def marker_score(log1p: pd.DataFrame, markers: list[str]) -> np.ndarray:
    present = [m for m in markers if m in log1p.columns]
    if not present:
        return np.zeros(len(log1p), dtype=float)
    return log1p[present].mean(axis=1).to_numpy(dtype=float)


def assign_lineage(log1p: pd.DataFrame) -> np.ndarray:
    scores = {
        "epithelial": marker_score(log1p, EPI_MARKERS),
        "tnk": marker_score(log1p, TNK_MARKERS),
        "myeloid": marker_score(log1p, MYE_MARKERS),
        "b": marker_score(log1p, B_MARKERS),
    }
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    winner = np.argmax(mat, axis=0)
    top = mat.max(axis=0)
    second = np.partition(mat, -2, axis=0)[-2]
    labels = np.array(["other"] * mat.shape[1], dtype=object)
    keep = (top >= 0.12) & (top >= second * 1.15)
    labels[keep] = np.array(names, dtype=object)[winner[keep]]
    return labels


def log1p_cp10k_from_counts(counts: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    numi = counts.sum(axis=1).to_numpy(dtype=float)
    scale = np.where(numi > 0, 1e4 / numi, 0.0)
    expr = np.log1p(counts.to_numpy(dtype=float) * scale[:, None])
    return pd.DataFrame(expr, index=counts.index, columns=counts.columns), numi


# ---------------------------------------------------------------------------
# GSE123902
# ---------------------------------------------------------------------------

def parse_123902_name(fname: str) -> dict:
    m = re.match(r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz", fname)
    if not m:
        raise ValueError(fname)
    return {"gsm": m.group(1), "patient": m.group(2), "site": m.group(3), "file": fname}


def extract_123902(raw_tar: Path, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_tar) as tf:
        tf.extractall(dest)
    return sorted(dest.glob("GSM*_dense.csv.gz"))


def load_123902_sample(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """Return wanted-gene counts plus full-library nUMI (all ~16k genes)."""
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
        genes = header[1:]
        keep_idx = np.array([i for i, g in enumerate(genes) if g in set(wanted)], dtype=int)
        keep_names = [genes[i] for i in keep_idx]
        barcodes: list[str] = []
        kept_rows: list[np.ndarray] = []
        numi: list[float] = []
        for line in handle:
            parts = line.rstrip("\n").split(",")
            barcodes.append(parts[0])
            vals = np.fromstring(",".join(parts[1:]), sep=",", dtype=np.float32)
            if vals.size != len(genes):
                raise ValueError(f"{path.name}: {vals.size} != {len(genes)}")
            numi.append(float(vals.sum()))
            kept_rows.append(vals[keep_idx] if keep_idx.size else np.zeros(0, dtype=np.float32))
    counts = pd.DataFrame(kept_rows, index=barcodes, columns=keep_names)
    return counts, np.asarray(numi, dtype=float)


def analyze_123902(data_dir: Path, extra_genes: list[str] | None = None) -> dict:
    raw = data_dir / "gse123902" / "GSE123902_RAW.tar"
    csv_dir = data_dir / "gse123902"
    files = sorted(csv_dir.glob("GSM*_dense.csv.gz"))
    if not files:
        files = extract_123902(raw, csv_dir)
    wanted = sorted(set(EPI_MARKERS + TNK_MARKERS + MYE_MARKERS + B_MARKERS + AUDIT_GENES + (extra_genes or [])))
    cell_rows = []
    sample_rows = []
    for path in files:
        meta = parse_123902_name(path.name)
        counts, numi = load_123902_sample(path, wanted)
        scale = np.where(numi > 0, 1e4 / numi, 0.0)
        log1p = pd.DataFrame(
            np.log1p(counts.to_numpy(dtype=float) * scale[:, None]),
            index=counts.index,
            columns=counts.columns,
        )
        lineage = assign_lineage(log1p)
        lineage[numi <= 0] = "other"
        tumor = meta["site"] != "NORMAL"
        mal = (lineage == "epithelial") & tumor
        tnk = lineage == "tnk"
        rec = {
            **meta,
            "n_cells": int(len(counts)),
            "n_epithelial": int((lineage == "epithelial").sum()),
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(tnk.mean()) if len(counts) else float("nan"),
            "mal_CLDN4_mean": float(log1p.loc[mal, "CLDN4"].mean()) if mal.any() and "CLDN4" in log1p else float("nan"),
            "mal_CLDN4_pct": float((log1p.loc[mal, "CLDN4"] > 0).mean() * 100) if mal.any() and "CLDN4" in log1p else float("nan"),
            "tnk_CLDN4_mean": float(log1p.loc[tnk, "CLDN4"].mean()) if tnk.any() and "CLDN4" in log1p else float("nan"),
            "median_n_umi": float(np.median(numi)) if len(numi) else float("nan"),
            "n_empty": int((numi <= 0).sum()),
            "bytes": path.stat().st_size,
        }
        sample_rows.append(rec)
        tmp = pd.DataFrame(
            {
                "barcode": counts.index.astype(str),
                "patient": meta["patient"],
                "site": meta["site"],
                "gsm": meta["gsm"],
                "lineage": lineage,
                "malignant": mal,
                "tnk": tnk,
                "tumor": tumor,
                "CLDN4": log1p["CLDN4"].to_numpy() if "CLDN4" in log1p else np.nan,
                "n_umi": numi,
            }
        )
        cell_rows.append(tmp)
        print(f"123902 {meta['gsm']} {meta['patient']} {meta['site']} n={rec['n_cells']} mal={rec['n_malignant']} tnk={rec['n_tnk']}", flush=True)
    samples = pd.DataFrame(sample_rows)
    cells = pd.concat(cell_rows, ignore_index=True)
    # Patient unit: tumor samples only. One patient can have tumour + normal.
    tumor_cells = cells[cells["tumor"]]
    patients = []
    for pid, sub in tumor_cells.groupby("patient"):
        mal = sub["malignant"]
        tnk = sub["tnk"]
        patients.append(
            {
                "cohort": "GSE123902",
                "patient": pid,
                "n_tumor_samples": int(sub["gsm"].nunique()),
                "sites": ",".join(sorted(sub["site"].unique())),
                "n_cells": int(len(sub)),
                "n_malignant": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "frac_tnk": float(tnk.mean()),
                "mal_CLDN4_mean": float(sub.loc[mal, "CLDN4"].mean()) if mal.any() else float("nan"),
                "mal_CLDN4_pct": float((sub.loc[mal, "CLDN4"] > 0).mean() * 100) if mal.any() else float("nan"),
                "tnk_CLDN4_mean": float(sub.loc[tnk, "CLDN4"].mean()) if tnk.any() else float("nan"),
                "eligible": bool(mal.sum() >= MIN_MAL and tnk.sum() >= MIN_TNK),
            }
        )
    patients_df = pd.DataFrame(patients)
    elig = patients_df[patients_df["eligible"]]
    rho, p, n = spearman(elig["mal_CLDN4_mean"], elig["frac_tnk"]) if len(elig) >= 3 else (float("nan"), float("nan"), int(len(elig)))
    lo, hi = spearman_ci(rho, n)
    effect = {
        "cohort": "GSE123902",
        "assay": "scRNA GEO dense UMI CSVs; marker epithelial in tumor vs T/NK",
        "malignant_def": "marker epithelial (EPCAM/KRT*) in PRIMARY/METASTASIS; not CNV; author H5 36.5GB skipped",
        "n_catalog_tumor_patients": int(patients_df["patient"].nunique()),
        "n_eligible": n,
        "rho": rho,
        "p": p,
        "ci95": [lo, hi],
        "min_mal": MIN_MAL,
        "min_tnk": MIN_TNK,
    }
    return {
        "samples": samples,
        "patients": patients_df,
        "cells": cells,
        "effect": effect,
    }


# ---------------------------------------------------------------------------
# GSE131907
# ---------------------------------------------------------------------------

def parse_131907_series(path: Path) -> pd.DataFrame:
    titles: list[str] = []
    patients: list[str] = []
    origins: list[str] = []
    stages: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                if vals and vals[0].startswith("patient id:"):
                    patients = [v.split(": ", 1)[-1] for v in vals]
                elif vals and vals[0].startswith("tumor stage:"):
                    stages = [v.split(": ", 1)[-1] for v in vals]
                elif vals and vals[0].startswith("tissue origin"):
                    origins = [v.split(": ", 1)[-1] for v in vals]
    return pd.DataFrame({"sample": titles, "patient": patients, "origin": origins, "stage": stages})


def stream_131907(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
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
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 5000 == 0:
                print(f"  131907 stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"131907 stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def analyze_131907(data_dir: Path, extra_genes: list[str] | None = None) -> dict:
    ann = pd.read_csv(data_dir / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    series = parse_131907_series(data_dir / "gse131907" / "GSE131907_series_matrix.txt.gz")
    wanted = set(["CLDN4"] + (extra_genes or []))
    cell_ids, found, n_umi, n_genes = stream_131907(
        data_dir / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", wanted
    )
    expr = pd.DataFrame({"Index": cell_ids})
    if "CLDN4" not in found:
        raise RuntimeError("CLDN4 absent from GSE131907 UMI matrix")
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    cldn4 = np.log1p(found["CLDN4"].astype(float) * scale)
    expr["CLDN4"] = cldn4
    expr["n_umi"] = n_umi
    merged = ann.merge(expr, on="Index", how="inner")
    merged = merged.merge(series, left_on="Sample", right_on="sample", how="left")
    merged["malignant"] = merged["Cell_subtype"].isin(MALIG_SUBTYPES_131907)
    merged["tnk"] = merged["Cell_type"].isin(TNK_TYPES_131907)
    merged["tumor_origin"] = merged["Sample_Origin"].isin(TUMOR_ORIGINS_131907)
    tumor = merged[merged["tumor_origin"]].copy()
    patients = []
    for pid, sub in tumor.groupby("patient"):
        mal = sub["malignant"]
        tnk = sub["tnk"]
        patients.append(
            {
                "cohort": "GSE131907",
                "patient": pid,
                "n_tumor_samples": int(sub["Sample"].nunique()),
                "sites": ",".join(sorted(sub["Sample_Origin"].dropna().unique())),
                "stage": ",".join(sorted(sub["stage"].dropna().astype(str).unique())),
                "n_cells": int(len(sub)),
                "n_malignant": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "frac_tnk": float(tnk.mean()),
                "mal_CLDN4_mean": float(sub.loc[mal, "CLDN4"].mean()) if mal.any() else float("nan"),
                "mal_CLDN4_pct": float((sub.loc[mal, "CLDN4"] > 0).mean() * 100) if mal.any() else float("nan"),
                "tnk_CLDN4_mean": float(sub.loc[tnk, "CLDN4"].mean()) if tnk.any() else float("nan"),
                "eligible": bool(mal.sum() >= MIN_MAL and tnk.sum() >= MIN_TNK),
            }
        )
    patients_df = pd.DataFrame(patients)
    elig = patients_df[patients_df["eligible"]]
    rho, p, n = spearman(elig["mal_CLDN4_mean"], elig["frac_tnk"]) if len(elig) >= 3 else (float("nan"), float("nan"), int(len(elig)))
    lo, hi = spearman_ci(rho, n)
    effect = {
        "cohort": "GSE131907",
        "assay": "scRNA GEO raw UMI + author annotation",
        "malignant_def": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3} at tLung/tL-B/mLN/mBrain",
        "n_catalog_tumor_patients": int(patients_df["patient"].nunique()),
        "n_eligible": n,
        "rho": rho,
        "p": p,
        "ci95": [lo, hi],
        "min_mal": MIN_MAL,
        "min_tnk": MIN_TNK,
        "n_genes_streamed": n_genes,
        "n_cells_matrix": int(len(cell_ids)),
    }
    return {
        "ann": ann,
        "series": series,
        "merged": merged,
        "found": found,
        "n_umi": n_umi,
        "cell_ids": cell_ids,
        "patients": patients_df,
        "effect": effect,
    }


# ---------------------------------------------------------------------------
# Combo + difference rule
# ---------------------------------------------------------------------------

def cohorts_differ(a: dict, b: dict) -> tuple[bool, str]:
    if not math.isfinite(a["rho"]) or not math.isfinite(b["rho"]):
        return False, "one or both Spearman undefined"
    if a["n_eligible"] < MIN_N_SPEARMAN or b["n_eligible"] < MIN_N_SPEARMAN:
        return False, f"honest n below {MIN_N_SPEARMAN} on at least one arm"
    reasons = []
    if (a["rho"] > 0) != (b["rho"] > 0) and min(abs(a["rho"]), abs(b["rho"])) >= 0.10:
        reasons.append("opposite sign")
    if abs(a["rho"] - b["rho"]) >= 0.25:
        reasons.append(f"|Δρ|={abs(a['rho']-b['rho']):.3f}≥0.25")
    sig_a, sig_b = a["p"] < 0.05, b["p"] < 0.05
    if sig_a != sig_b and abs(a["rho"] - b["rho"]) >= 0.20:
        reasons.append("one arm p<0.05 and |Δρ|≥0.20")
    re = random_effects_dl([a["rho"], b["rho"]], [a["n_eligible"], b["n_eligible"]])
    if re.get("I2", 0) >= 50:
        reasons.append(f"I²={re['I2']:.1f}%≥50%")
    if reasons:
        return True, "; ".join(reasons)
    return False, "cohort rhos do not differ by the pre-specified rule"


def combo_table(effects: list[dict]) -> tuple[pd.DataFrame, dict]:
    rows = []
    for e in effects:
        lo, hi = e["ci95"]
        rows.append(
            {
                "cohort": e["cohort"],
                "assay": e["assay"],
                "n": e["n_eligible"],
                "rho": e["rho"],
                "ci95_lo": lo,
                "ci95_hi": hi,
                "p": e["p"],
                "malignant_def": e["malignant_def"],
            }
        )
    rhos = [e["rho"] for e in effects]
    ns = [e["n_eligible"] for e in effects]
    ps = [e["p"] for e in effects]
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    fi = fisher_combine(ps)
    rows.append(
        {
            "cohort": "RE_combo",
            "assay": "Fisher-z DerSimonian-Laird",
            "n": re.get("n_patients_total", sum(ns)),
            "rho": re.get("pooled_rho", float("nan")),
            "ci95_lo": re.get("ci95_rho", [float("nan"), float("nan")])[0],
            "ci95_hi": re.get("ci95_rho", [float("nan"), float("nan")])[1],
            "p": re.get("p", float("nan")),
            "malignant_def": f"I2={re.get('I2', float('nan')):.1f}%",
        }
    )
    return pd.DataFrame(rows), {"re": re, "stouffer": st, "fisher": fi}


def forest_plot(table: pd.DataFrame, path: Path, title: str) -> None:
    plot = table.copy()
    fig, ax = plt.subplots(figsize=(7.2, 2.6 + 0.35 * len(plot)))
    y = np.arange(len(plot))[::-1]
    for i, rec in enumerate(plot.itertuples(index=False)):
        yi = y[i]
        color = "#1f4e79" if rec.cohort == "RE_combo" else "#4c78a8"
        if math.isfinite(rec.rho) and math.isfinite(rec.ci95_lo):
            ax.plot([rec.ci95_lo, rec.ci95_hi], [yi, yi], color=color, lw=2)
            ax.plot(rec.rho, yi, "o", color=color, ms=7)
        ax.text(1.02, yi, f"n={int(rec.n)}" if math.isfinite(rec.n) else "n=?", transform=ax.get_yaxis_transform(), va="center", fontsize=8)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["cohort"])
    ax.set_xlabel("Spearman ρ, malignant CLDN4 vs T/NK fraction")
    ax.set_title(title)
    ax.set_xlim(-1.05, 1.05)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_plot(df: pd.DataFrame, title: str, path: Path) -> None:
    elig = df[df["eligible"]]
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(elig["mal_CLDN4_mean"], elig["frac_tnk"], c="#1f4e79", s=36)
    ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# CellChat-style (only if the pair differs)
# ---------------------------------------------------------------------------

def pair_class(ligand: str) -> str:
    tags = []
    if ligand in FOCUSED_PAIRS["MHC-I"] or ligand.startswith("HLA-"):
        tags.append("MHC-I")
    if ligand in FOCUSED_PAIRS["T-recruit"]:
        tags.append("T-recruit")
    if ligand in FOCUSED_PAIRS["checkpoint"]:
        tags.append("checkpoint")
    if ligand in FOCUSED_PAIRS["barrier"]:
        tags.append("barrier")
    return "|".join(tags) if tags else "other"


def cellchat_outgoing(mal_high: pd.DataFrame, mal_low: pd.DataFrame, tnk: pd.DataFrame, lr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    def pack(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        mat = df[genes].to_numpy(dtype=float)
        pos = mat > 0
        return mat.T, pos.T

    groups = ["Mal_high", "Mal_low", "TNK"]
    mats = [pack(mal_high), pack(mal_low), pack(tnk)]
    expr = np.concatenate([m[0] for m in mats], axis=1)
    pos = np.concatenate([m[1] for m in mats], axis=1)
    labels = np.array(
        ["Mal_high"] * len(mal_high) + ["Mal_low"] * len(mal_low) + ["TNK"] * len(tnk)
    )
    n_g = expr.shape[0]
    means = np.zeros((n_g, 3))
    props = np.zeros((n_g, 3))
    counts = {}
    for j, g in enumerate(groups):
        idx = np.flatnonzero(labels == g)
        counts[g] = int(idx.size)
        if idx.size == 0:
            continue
        sub = expr[:, idx]
        means[:, j] = sub[:, 0] if idx.size == 1 else trim_mean(sub, TRIM, axis=1)
        props[:, j] = pos[:, idx].mean(axis=1)
    gix = {g: i for i, g in enumerate(genes)}

    def complex_mu_pr(subunits):
        ix = [gix[g] for g in subunits]
        mu = geom_mean_rows(means[ix])
        pr = props[ix].min(axis=0) if len(ix) > 1 else props[ix[0]]
        return mu, pr

    rows = []
    for rec in lr.itertuples(index=False):
        if any(g not in gix for g in rec.ligand_genes + rec.receptor_genes):
            continue
        lig_mu, lig_pr = complex_mu_pr(rec.ligand_genes)
        rec_mu, rec_pr = complex_mu_pr(rec.receptor_genes)
        for src, si in (("Mal_high", 0), ("Mal_low", 1)):
            detected = (lig_pr[si] >= EXPR_PROP) and (rec_pr[2] >= EXPR_PROP)
            prob = hill_prob(float(lig_mu[si]), float(rec_mu[2])) if detected else 0.0
            rows.append(
                {
                    "interaction_name": rec.interaction_name,
                    "pathway_name": rec.pathway_name,
                    "ligand": rec.ligand,
                    "receptor": rec.receptor,
                    "source": src,
                    "target": "TNK",
                    "n_source": counts[src],
                    "n_target": counts["TNK"],
                    "detected": bool(detected),
                    "P": prob,
                    "lig_prop": float(lig_pr[si]),
                    "rec_prop": float(rec_pr[2]),
                    "pair_class": pair_class(rec.ligand),
                }
            )
    wide = pd.DataFrame(rows)
    if wide.empty:
        return wide
    hi = wide[wide["source"] == "Mal_high"].set_index("interaction_name")
    lo = wide[wide["source"] == "Mal_low"].set_index("interaction_name")
    out = hi[["pathway_name", "ligand", "receptor", "n_source", "n_target", "pair_class"]].copy()
    out = out.rename(columns={"n_source": "n_mal_high"})
    out["n_mal_low"] = lo["n_source"]
    out["P_high"] = hi["P"]
    out["P_low"] = lo["P"]
    out["delta_P"] = out["P_high"] - out["P_low"]
    out["detected_high"] = hi["detected"]
    out["detected_low"] = lo["detected"]
    out["lig_prop_high"] = hi["lig_prop"]
    out["lig_prop_low"] = lo["lig_prop"]
    return out.reset_index()


def patient_level_cellchat(cells: pd.DataFrame, gene_cols: list[str], lr: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """Within-patient median split of malignant CLDN4; outgoing to same-patient T/NK."""
    rows = []
    for pid, sub in cells.groupby("patient"):
        mal = sub[sub["malignant"]]
        tnk = sub[sub["tnk"]]
        if len(mal) < 2 * MIN_MAL or len(tnk) < MIN_TNK:
            continue
        cut = float(mal["CLDN4"].median())
        hi = mal[mal["CLDN4"] >= cut]
        lo = mal[mal["CLDN4"] < cut]
        if len(hi) < MIN_GROUP or len(lo) < MIN_GROUP:
            continue
        tab = cellchat_outgoing(hi, lo, tnk, lr, gene_cols)
        if tab.empty:
            continue
        tab["patient"] = pid
        tab["cohort"] = cohort
        tab["n_tnk"] = int(len(tnk))
        rows.append(tab)
    if not rows:
        return pd.DataFrame()
    allp = pd.concat(rows, ignore_index=True)
    # Wilcoxon on per-patient P for focused / detected pairs
    summary = []
    for name, g in allp.groupby("interaction_name"):
        if len(g) < 5:
            continue
        try:
            stat, p = wilcoxon(g["P_high"], g["P_low"], zero_method="wilcox", alternative="two-sided")
        except ValueError:
            stat, p = float("nan"), float("nan")
        summary.append(
            {
                "cohort": cohort,
                "interaction_name": name,
                "pathway_name": g["pathway_name"].iloc[0],
                "ligand": g["ligand"].iloc[0],
                "receptor": g["receptor"].iloc[0],
                "pair_class": g["pair_class"].iloc[0],
                "n_patients": int(len(g)),
                "mean_P_high": float(g["P_high"].mean()),
                "mean_P_low": float(g["P_low"].mean()),
                "mean_delta_P": float(g["delta_P"].mean()),
                "wilcoxon_stat": float(stat) if math.isfinite(stat) else float("nan"),
                "p": float(p) if math.isfinite(p) else float("nan"),
                "n_mal_high_median": float(g["n_mal_high"].median()),
                "n_mal_low_median": float(g["n_mal_low"].median()),
                "n_tnk_median": float(g["n_tnk"].median()),
            }
        )
    return pd.DataFrame(summary).sort_values(["p", "mean_delta_P"], na_position="last")


def cellchat_123902(data_dir: Path, cells_meta: pd.DataFrame, db_dir: Path) -> pd.DataFrame:
    wanted = sorted(lr_wanted_genes(db_dir) | set(AUDIT_GENES + EPI_MARKERS + TNK_MARKERS))
    csv_dir = data_dir / "gse123902"
    frames = []
    present_genes: set[str] | None = None
    for path in sorted(csv_dir.glob("GSM*_dense.csv.gz")):
        meta = parse_123902_name(path.name)
        if meta["site"] == "NORMAL":
            continue
        counts, numi = load_123902_sample(path, wanted)
        scale = np.where(numi > 0, 1e4 / numi, 0.0)
        log1p = pd.DataFrame(
            np.log1p(counts.to_numpy(dtype=float) * scale[:, None]),
            index=counts.index,
            columns=counts.columns,
        )
        if present_genes is None:
            present_genes = set(log1p.columns)
        log1p = log1p.reset_index().rename(columns={"index": "barcode"})
        log1p["barcode"] = log1p["barcode"].astype(str)
        log1p["patient"] = meta["patient"]
        log1p["gsm"] = meta["gsm"]
        frames.append(log1p)
    expr = pd.concat(frames, ignore_index=True)
    meta = cells_meta[cells_meta["tumor"]][["barcode", "patient", "malignant", "tnk", "CLDN4"]].copy()
    meta["barcode"] = meta["barcode"].astype(str)
    # barcodes may collide across samples; join on patient+barcode
    expr["barcode"] = expr["barcode"].astype(str)
    merged = expr.merge(meta.drop(columns=["CLDN4"]), on=["barcode", "patient"], how="inner")
    genes = [g for g in wanted if g in merged.columns]
    lr = load_lr(db_dir, set(genes))
    print(f"123902 CellChat genes={len(genes)} pairs={len(lr)} cells={len(merged)}", flush=True)
    return patient_level_cellchat(merged, genes, lr, "GSE123902")


def cellchat_131907(bundle: dict, extra_found: dict[str, np.ndarray], db_dir: Path) -> pd.DataFrame:
    merged = bundle["merged"]
    tumor = merged[merged["tumor_origin"]].copy()
    scale = np.where(bundle["n_umi"] > 0, 1e4 / bundle["n_umi"], 0.0)
    id_to_i = {cid: i for i, cid in enumerate(bundle["cell_ids"])}
    genes = ["CLDN4"]
    for g, arr in extra_found.items():
        if g == "CLDN4":
            continue
        tumor[g] = np.log1p(arr[[id_to_i[x] for x in tumor["Index"]]].astype(float) * scale[[id_to_i[x] for x in tumor["Index"]]])
        genes.append(g)
    # CLDN4 already on tumor
    lr = load_lr(db_dir, set(genes))
    print(f"131907 CellChat genes={len(genes)} pairs={len(lr)} cells={len(tumor)}", flush=True)
    return patient_level_cellchat(tumor, genes, lr, "GSE131907")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--skip-cellchat", action="store_true")
    parser.add_argument("--reuse-131907", action="store_true", help="Reuse results/gse131907_patients.tsv")
    args = parser.parse_args()
    ensure_dirs()
    db_dir = HERE / "db"

    print("=== GSE123902 patient-level ===", flush=True)
    g123 = analyze_123902(args.data)
    g123["samples"].to_csv(OUT / "gse123902_samples.tsv", sep="\t", index=False)
    g123["patients"].to_csv(OUT / "gse123902_patients.tsv", sep="\t", index=False)
    scatter_plot(g123["patients"], "GSE123902 malignant CLDN4 vs T/NK", FIG / "scatter_gse123902")

    print("=== GSE131907 patient-level ===", flush=True)
    reuse_path = OUT / "gse131907_patients.tsv"
    if args.reuse_131907 and reuse_path.exists():
        patients_df = pd.read_csv(reuse_path, sep="\t")
        elig = patients_df[patients_df["eligible"]]
        rho, p, n = spearman(elig["mal_CLDN4_mean"], elig["frac_tnk"]) if len(elig) >= 3 else (float("nan"), float("nan"), int(len(elig)))
        lo, hi = spearman_ci(rho, n)
        g131 = {
            "patients": patients_df,
            "effect": {
                "cohort": "GSE131907",
                "assay": "scRNA GEO raw UMI + author annotation",
                "malignant_def": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3} at tLung/tL-B/mLN/mBrain",
                "n_catalog_tumor_patients": int(patients_df["patient"].nunique()),
                "n_eligible": n,
                "rho": rho,
                "p": p,
                "ci95": [lo, hi],
                "min_mal": MIN_MAL,
                "min_tnk": MIN_TNK,
                "reused": True,
            },
        }
        print("reused", reuse_path, flush=True)
    else:
        g131 = analyze_131907(args.data)
        g131["patients"].to_csv(OUT / "gse131907_patients.tsv", sep="\t", index=False)
        g131["series"].to_csv(OUT / "gse131907_series.tsv", sep="\t", index=False)
    scatter_plot(g131["patients"], "GSE131907 malignant CLDN4 vs T/NK", FIG / "scatter_gse131907")

    effects = [g123["effect"], g131["effect"]]
    table, pooled = combo_table(effects)
    table.to_csv(OUT / "combo_cldn4_tnk.tsv", sep="\t", index=False)
    forest_plot(table, FIG / "forest_CLDN4_tnk", "Pairwise combo: malignant CLDN4 vs T/NK")

    differ, why = cohorts_differ(g123["effect"], g131["effect"])
    print(f"differ={differ} ({why})", flush=True)
    print(table.to_string(index=False), flush=True)

    cellchat_ran = False
    ligand_note = "not run; cohort rhos do not differ"
    if differ and not args.skip_cellchat:
        print("=== CellChat-style (pair differs) ===", flush=True)
        wanted = sorted(lr_wanted_genes(db_dir))
        # Re-stream 131907 with LR genes (CLDN4 already stored).
        extra_wanted = set(wanted) - set(g131["found"])
        if extra_wanted:
            print(f"re-stream 131907 for {len(extra_wanted)} LR genes", flush=True)
            _, extra_found, _, _ = stream_131907(
                args.data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
                extra_wanted,
            )
        else:
            extra_found = {}
        extra_found.update(g131["found"])
        tab123 = cellchat_123902(args.data, g123["cells"], db_dir)
        tab131 = cellchat_131907(g131, extra_found, db_dir)
        if not tab123.empty:
            tab123.to_csv(OUT / "ligand_table_gse123902.tsv", sep="\t", index=False)
        if not tab131.empty:
            tab131.to_csv(OUT / "ligand_table_gse131907.tsv", sep="\t", index=False)
        parts = [t for t in (tab123, tab131) if not t.empty]
        if parts:
            both = pd.concat(parts, ignore_index=True)
            both.to_csv(OUT / "ligand_table.tsv", sep="\t", index=False)
        cellchat_ran = True
        ligand_note = (
            f"patient-level within-patient median CLDN4 split; "
            f"GSE123902 n_pairs={0 if tab123.empty else len(tab123)}; "
            f"GSE131907 n_pairs={0 if tab131.empty else len(tab131)}"
        )

    summary = {
        "pair": ["GSE123902", "GSE131907"],
        "gene": "CLDN4",
        "dual_high": False,
        "pile": False,
        "effects": effects,
        "pooled": {
            "re": {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in pooled["re"].items()},
            "stouffer": pooled["stouffer"],
            "fisher": pooled["fisher"],
        },
        "differ": differ,
        "differ_reason": why,
        "cellchat_ran": cellchat_ran,
        "ligand_note": ligand_note,
        "gates": {
            "gse123902_raw_tar_mb": 90.4,
            "gse123902_author_h5_gb": 36.5,
            "gse123902_author_h5_skipped": True,
            "gse131907_umi_mb": 389.8,
            "gse131907_log2tpm_gb": 2.9,
            "gse131907_log2tpm_skipped": True,
            "tisch_gse123902": "404",
        },
        "min_mal": MIN_MAL,
        "min_tnk": MIN_TNK,
        "min_n_spearman": MIN_N_SPEARMAN,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps({"differ": differ, "why": why, "effects": effects}, indent=2, default=float))


if __name__ == "__main__":
    main()
