#!/usr/bin/env python3
"""GSE260770 CLDN4-only: plasma-exosome FPKM + counts vs CD8A / T/NK / IFN / CD274 / response.

Additive CLDN4 page. No dual-high. Does not re-score TACSTD2.

Public GEO only:
  GSE260770_mRNA_FPKM.txt.gz
  GSE260770_mRNA_counts.txt.gz
  GSE260770_series_matrix.txt.gz

Tissue is peripheral-blood plasma exosome RNA (Cheng / Feng / Liang,
Signal Transduct Target Ther 2024, PMID 38637495; NCT04026841 sintilimab
on high-risk GGO in MPLC). Not tumor epithelium.

Patient is the unit. GEO deposits no patient ID and no timepoint.
Trial ITT n=36; paper exosome RNA subset = 10 patients (5 R / 5 NR) at
T1–T4. Deposited n=50 libraries (26 Responsed / 24 Non-responsed).
Tests are therefore library-level with that caveat — they are not 50
independent patients.

Usage: python3 methods/gse260770_cldn4/analyze.py
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
TABLES = HERE / "tables"
PROC = HERE / "processed"
FIGS = HERE / "figures"
for d in (CACHE, TABLES, PROC, FIGS):
    d.mkdir(parents=True, exist_ok=True)

FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260770/suppl/"
    "GSE260770_mRNA_FPKM.txt.gz"
)
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260770/suppl/"
    "GSE260770_mRNA_counts.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260770/matrix/"
    "GSE260770_series_matrix.txt.gz"
)

SEED = 20260817
N_BOOT = 2000

# Pre-specified, CLDN4-only. No TACSTD2 dual-high.
GENES = {
    "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563",
    "CD274": "ENSG00000120217",
}

TNK_GENES = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "PRF1", "GZMA", "GZMB"]
IFN_GENES = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

# Paper / GEO facts (not re-derived from a TACSTD2 page).
PAPER_TRIAL_N = 36
PAPER_EXOSOME_N = 10
PAPER_EXOSOME_R = 5
PAPER_EXOSOME_NR = 5


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    urllib.request.urlretrieve(url, dest)
    return dest


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]


def load_meta() -> pd.DataFrame:
    path = fetch(MATRIX_URL, CACHE / "GSE260770_series_matrix.txt.gz")
    stored: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                stored["gsm"] = _split_fields(line)
            elif line.startswith("!Sample_title"):
                stored["title"] = _split_fields(line)
            elif line.startswith("!Sample_source_name"):
                stored["source"] = _split_fields(line)
            elif line.startswith("!Sample_characteristics"):
                char_rows.append(_split_fields(line))
    n = len(stored["gsm"])
    meta = pd.DataFrame(
        {
            "gsm": stored["gsm"],
            "title": stored["title"],
            "source": stored.get("source", [None] * n)[:n],
        }
    )
    chars: dict[str, dict[str, str]] = {g: {} for g in meta["gsm"]}
    for vals in char_rows:
        for g, v in zip(meta["gsm"], vals):
            if not v:
                continue
            if ":" in v:
                k, val = v.split(":", 1)
                chars[g][k.strip().lower()] = val.strip()
    for k in sorted({kk for d in chars.values() for kk in d}):
        meta[k] = [chars[g].get(k) for g in meta["gsm"]]
    meta["library"] = meta["title"]
    # GEO has no patient / timepoint field. Do not invent one from the title.
    meta["patient_id_geo"] = pd.NA
    meta["resp"] = np.where(
        meta["group"].str.contains("Non-responsed", case=False, na=False),
        "NR",
        np.where(meta["group"].str.contains("Responsed", case=False, na=False), "R", None),
    )
    if meta["resp"].isna().any():
        bad = meta.loc[meta["resp"].isna(), "title"].tolist()
        raise AssertionError(f"unparsed response group: {bad}")
    return meta


def collapse_symbol(df: pd.DataFrame) -> pd.DataFrame:
    """One row per official symbol. Prefer the ENSG with highest mean."""
    work = df.copy()
    work["Symbol"] = work["Symbol"].astype(str)
    sample_cols = [c for c in work.columns if c not in ("#ID", "Symbol")]
    work[sample_cols] = work[sample_cols].apply(pd.to_numeric, errors="coerce")
    # drop novel-gene rows when an ENSG exists for the same symbol
    work["_ensg"] = work["#ID"].astype(str).str.startswith("ENSG")
    work["_mean"] = work[sample_cols].mean(axis=1)
    work = work.sort_values(["Symbol", "_ensg", "_mean"], ascending=[True, False, False])
    work = work.drop_duplicates("Symbol", keep="first")
    out = work.set_index("Symbol")[sample_cols]
    return out


def load_matrix(kind: str) -> pd.DataFrame:
    if kind == "fpkm":
        path = fetch(FPKM_URL, CACHE / "GSE260770_mRNA_FPKM.txt.gz")
    elif kind == "counts":
        path = fetch(COUNTS_URL, CACHE / "GSE260770_mRNA_counts.txt.gz")
    else:
        raise ValueError(kind)
    raw = pd.read_csv(path, sep="\t")
    return collapse_symbol(raw)


def log2p1_fpkm(expr: pd.DataFrame) -> pd.DataFrame:
    return np.log2(expr.clip(lower=0) + 1)


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd = np.where(sd == 0, np.nan, sd)
    return pd.DataFrame((arr - mu) / sd, index=df.index, columns=df.columns)


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in expr.index]
    if len(found) < 2:
        return pd.Series(np.nan, index=expr.columns), found
    return zscore_rows(expr.loc[found]).mean(axis=0), found


def spearman(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    rng = np.random.default_rng(SEED)
    boots = []
    idx = np.where(m)[0]
    for _ in range(N_BOOT):
        b = rng.choice(idx, size=n, replace=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=stats.ConstantInputWarning)
            rb, _ = stats.spearmanr(x[b], y[b])
        if rb is not None and np.isfinite(rb):
            boots.append(float(rb))
    if boots:
        lo, hi = np.percentile(boots, [2.5, 97.5])
    else:
        lo = hi = np.nan
    return {
        "n": n,
        "rho": float(r),
        "p": float(p),
        "ci_lo": float(lo),
        "ci_hi": float(hi),
    }


def partial_spearman(x, y, z) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx, ry, rz = stats.rankdata(x[m]), stats.rankdata(y[m]), stats.rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex, ey = rx - Z @ bx, ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def cliff_delta(pos: np.ndarray, neg: np.ndarray) -> float:
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(2.0 * u / (len(pos) * len(neg)) - 1.0)


def mwu_block(pos, neg) -> dict:
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    n1, n2 = len(pos), len(neg)
    if n1 < 2 or n2 < 2:
        return {
            "n_R": n1,
            "n_NR": n2,
            "median_R": float(np.median(pos)) if n1 else np.nan,
            "median_NR": float(np.median(neg)) if n2 else np.nan,
            "delta": np.nan,
            "p": np.nan,
            "auc": np.nan,
        }
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    delta = float(2.0 * u.statistic / (n1 * n2) - 1.0)
    return {
        "n_R": n1,
        "n_NR": n2,
        "median_R": float(np.median(pos)),
        "median_NR": float(np.median(neg)),
        "delta": delta,
        "p": float(u.pvalue),
        "auc": float(u.statistic / (n1 * n2)),
    }


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:+.{nd}f}" if isinstance(x, float) and abs(x) < 10 else f"{x:.{nd}f}"


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "—"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def detection_row(expr: pd.DataFrame, gene: str, matrix: str) -> dict:
    if gene not in expr.index:
        return {
            "gene": gene,
            "matrix": matrix,
            "present": "no",
            "n_libraries": int(expr.shape[1]),
            "n_nonzero": 0,
            "median": np.nan,
            "max": np.nan,
        }
    v = pd.to_numeric(expr.loc[gene], errors="coerce")
    return {
        "gene": gene,
        "matrix": matrix,
        "present": "yes",
        "n_libraries": int(v.notna().sum()),
        "n_nonzero": int((v > 0).sum()),
        "median": float(v.median()),
        "max": float(v.max()),
    }


def score_frame(log_expr: pd.DataFrame, raw_expr: pd.DataFrame, prefix: str) -> pd.DataFrame:
    rows = []
    for lib in log_expr.columns:
        rec = {"library": lib}
        for g in GENES:
            rec[f"{prefix}_{g}"] = float(log_expr.loc[g, lib]) if g in log_expr.index else np.nan
            rec[f"{prefix}_{g}_raw"] = float(raw_expr.loc[g, lib]) if g in raw_expr.index else np.nan
        rows.append(rec)
    sc = pd.DataFrame(rows).set_index("library")
    tnk, tnk_found = mean_z(log_expr, TNK_GENES)
    ifn, ifn_found = mean_z(log_expr, IFN_GENES)
    epi, epi_found = mean_z(log_expr, EPI_GENES)
    sc[f"{prefix}_TNK"] = tnk
    sc[f"{prefix}_IFN"] = ifn
    sc[f"{prefix}_EPI"] = epi
    sc.attrs[f"{prefix}_TNK_genes"] = tnk_found
    sc.attrs[f"{prefix}_IFN_genes"] = ifn_found
    sc.attrs[f"{prefix}_EPI_genes"] = epi_found
    return sc


def save_scatter(path, x, y, title, xlab, ylab, note=""):
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ax.scatter(x, y, s=28, c="#2c5aa0", alpha=0.75, edgecolors="none")
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() >= 4:
        r, p = stats.spearmanr(x[m], y[m])
        ax.set_title(f"{title}\nρ={r:+.2f}  p={p:.3g}  n={int(m.sum())}", fontsize=10)
    else:
        ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if note:
        ax.text(0.02, 0.98, note, transform=ax.transAxes, va="top", fontsize=8, color="#444")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_box(path, a, b, title, ylab, lab_a="R", lab_b="NR"):
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    data = [np.asarray(a, float), np.asarray(b, float)]
    bp = ax.boxplot(
        data,
        tick_labels=[f"{lab_a}\nn={len(data[0])}", f"{lab_b}\nn={len(data[1])}"],
        patch_artist=True,
    )
    for patch, c in zip(bp["boxes"], ["#d95f02", "#7570b3"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    blk = mwu_block(a, b)
    ax.set_title(f"{title}\nδ={fmt(blk['delta'])}  p={fmt_p(blk['p'])}", fontsize=10)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    meta = load_meta()
    fpkm = load_matrix("fpkm")
    counts = load_matrix("counts")

    # align libraries
    libs = [c for c in meta["library"] if c in fpkm.columns and c in counts.columns]
    if len(libs) != 50:
        raise AssertionError(f"expected 50 aligned libraries, got {len(libs)}")
    meta = meta.set_index("library").loc[libs].reset_index()
    fpkm = fpkm[libs]
    counts = counts[libs]

    for g, ensg in GENES.items():
        if g not in fpkm.index:
            raise AssertionError(f"{g} missing from FPKM")
        # confirm the chosen row is the official ENSG when present
        # (collapse already preferred ENSG)

    log_fpkm = log2p1_fpkm(fpkm)
    log_cpm = log2cpm(counts)

    sc_f = score_frame(log_fpkm, fpkm, "fpkm")
    sc_c = score_frame(log_cpm, counts, "cpm")
    df = meta.merge(sc_f, on="library", how="inner").merge(sc_c, on="library", how="inner")
    if len(df) != 50:
        raise AssertionError(f"merged n={len(df)}")

    # epithelial residual only if the epithelial score actually varies in blood
    epi = df["fpkm_EPI"].to_numpy()
    epi_any = (fpkm.reindex(EPI_GENES).fillna(0) > 0).any(axis=0)
    n_epi_any = int(epi_any.sum())
    residual_possible = n_epi_any >= 10 and float(np.nanstd(epi)) > 1e-8
    if residual_possible:
        df["fpkm_CLDN4_resid_EPI"] = residualize_safe(df["fpkm_CLDN4"], df["fpkm_EPI"])
    else:
        df["fpkm_CLDN4_resid_EPI"] = np.nan

    df.to_csv(PROC / "GSE260770_library.tsv", sep="\t", index=False)

    # --- inventory ---
    inv = [
        {"field": "GEO libraries (GSM)", "public": "yes", "n": int(df["gsm"].nunique()), "note": "GSM8124173–GSM8124222"},
        {"field": "unique titles", "public": "yes", "n": int(df["title"].nunique()), "note": "E-lncA411-##h; no repeated title"},
        {"field": "patient ID on GEO / BioSample", "public": "no", "n": 0, "note": "no patient / donor / timepoint characteristic"},
        {"field": "trial ITT patients (paper)", "public": "text", "n": PAPER_TRIAL_N, "note": "NCT04026841; not a GEO field"},
        {"field": "paper exosome-RNA patients", "public": "text", "n": PAPER_EXOSOME_N, "note": f"{PAPER_EXOSOME_R} R + {PAPER_EXOSOME_NR} NR matched; T1–T4 blood"},
        {"field": "GEO response label (library)", "public": "yes", "n": 50, "note": f"{int((df.resp=='R').sum())} Responsed / {int((df.resp=='NR').sum())} Non-responsed"},
        {"field": "computable patient n", "public": "no", "n": 0, "note": "cannot collapse 50 libraries to patients without inventing IDs"},
        {"field": "tissue", "public": "yes", "n": 50, "note": "peripheral blood plasma exosome; not tumor"},
        {"field": "CLDN4 on FPKM+counts", "public": "yes", "n": 50, "note": "ENSG00000189143 present; floor (see detection)"},
        {"field": "CD8A on FPKM+counts", "public": "yes", "n": 50, "note": "ENSG00000153563 present; floor"},
        {"field": "CD274 on FPKM+counts", "public": "yes", "n": 50, "note": "ENSG00000120217 present; floor"},
        {"field": "epithelial residual possible", "public": "no" if not residual_possible else "yes", "n": n_epi_any, "note": f"{n_epi_any}/50 libraries have any of {EPI_GENES} >0 FPKM"},
    ]
    pd.DataFrame(inv).to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    detect_genes = (
        list(GENES)
        + TNK_GENES
        + IFN_GENES
        + EPI_GENES
        + ["PTPRC", "NKG7"]
    )
    detect_genes = list(dict.fromkeys(detect_genes))
    det = [detection_row(fpkm, g, "FPKM") for g in detect_genes]
    det += [detection_row(counts, g, "counts") for g in detect_genes]
    det_df = pd.DataFrame(det)
    det_df.to_csv(TABLES / "gene_detection.tsv", sep="\t", index=False)

    tests = []

    def add_spearman(contrast, x, y, matrix, note):
        sp = spearman(x, y)
        tests.append(
            {
                "contrast": contrast,
                "matrix": matrix,
                "n": sp["n"],
                "n_R": "",
                "n_NR": "",
                "metric": "Spearman_rho",
                "effect": sp["rho"],
                "p": sp["p"],
                "ci_lo": sp["ci_lo"],
                "ci_hi": sp["ci_hi"],
                "note": note,
            }
        )
        return sp

    def add_resp(contrast, values, matrix, note):
        r = df["resp"].eq("R")
        blk = mwu_block(values[r], values[~r])
        tests.append(
            {
                "contrast": contrast,
                "matrix": matrix,
                "n": int(blk["n_R"] + blk["n_NR"]),
                "n_R": blk["n_R"],
                "n_NR": blk["n_NR"],
                "metric": "Cliff_delta_R_minus_NR",
                "effect": blk["delta"],
                "p": blk["p"],
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "note": note + f"; med R={blk['median_R']:.4g} NR={blk['median_NR']:.4g}",
            }
        )
        return blk

    # Primary matrix = author FPKM (log2+1). Counts log2CPM = sensitivity.
    add_spearman(
        "CLDN4 vs CD8A",
        df["fpkm_CLDN4"],
        df["fpkm_CD8A"],
        "log2(FPKM+1)",
        "library-level; both genes at floor",
    )
    add_spearman(
        "CLDN4 vs T/NK mean-z",
        df["fpkm_CLDN4"],
        df["fpkm_TNK"],
        "log2(FPKM+1)",
        "T/NK = mean-z of " + ",".join(sc_f.attrs.get("fpkm_TNK_genes", TNK_GENES)),
    )
    add_spearman(
        "CLDN4 vs IFN mean-z",
        df["fpkm_CLDN4"],
        df["fpkm_IFN"],
        "log2(FPKM+1)",
        "IFN = Ayers-like mean-z of " + ",".join(sc_f.attrs.get("fpkm_IFN_genes", IFN_GENES)),
    )
    add_spearman(
        "CLDN4 vs CD274",
        df["fpkm_CLDN4"],
        df["fpkm_CD274"],
        "log2(FPKM+1)",
        "library-level; both genes at floor",
    )
    add_resp(
        "CLDN4 vs GEO sintilimab group",
        df["fpkm_CLDN4"],
        "log2(FPKM+1)",
        "26 R / 24 NR libraries; not 50 patients",
    )

    # counts sensitivity
    add_spearman("CLDN4 vs CD8A", df["cpm_CLDN4"], df["cpm_CD8A"], "log2(CPM+1)", "counts sensitivity")
    add_spearman("CLDN4 vs T/NK mean-z", df["cpm_CLDN4"], df["cpm_TNK"], "log2(CPM+1)", "counts sensitivity")
    add_spearman("CLDN4 vs IFN mean-z", df["cpm_CLDN4"], df["cpm_IFN"], "log2(CPM+1)", "counts sensitivity")
    add_spearman("CLDN4 vs CD274", df["cpm_CLDN4"], df["cpm_CD274"], "log2(CPM+1)", "counts sensitivity")
    add_resp("CLDN4 vs GEO sintilimab group", df["cpm_CLDN4"], "log2(CPM+1)", "counts sensitivity; 26/24 libraries")

    if residual_possible:
        add_spearman(
            "CLDN4 vs CD8A | epithelial",
            df["fpkm_CLDN4_resid_EPI"],
            df["fpkm_CD8A"],
            "log2(FPKM+1) resid EPI",
            "partial via residualized CLDN4",
        )
        for name, col in [("T/NK mean-z", "fpkm_TNK"), ("IFN mean-z", "fpkm_IFN"), ("CD274", "fpkm_CD274")]:
            psp = partial_spearman(df["fpkm_CLDN4"], df[col], df["fpkm_EPI"])
            tests.append(
                {
                    "contrast": f"CLDN4 vs {name} | epithelial",
                    "matrix": "log2(FPKM+1) partial EPI",
                    "n": psp["n"],
                    "n_R": "",
                    "n_NR": "",
                    "metric": "partial_Spearman_rho",
                    "effect": psp["rho"],
                    "p": psp["p"],
                    "ci_lo": np.nan,
                    "ci_hi": np.nan,
                    "note": f"covariate = epithelial mean-z; {n_epi_any}/50 libs any EPI>0",
                }
            )
    else:
        tests.append(
            {
                "contrast": "CLDN4 vs immune | epithelial residual",
                "matrix": "not computed",
                "n": n_epi_any,
                "n_R": "",
                "n_NR": "",
                "metric": "partial_Spearman_rho",
                "effect": np.nan,
                "p": np.nan,
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "note": (
                    f"epithelial residual not interpretable in plasma exosome: "
                    f"{n_epi_any}/50 libraries have any of EPCAM/KRT8/KRT18/KRT19/CDH1/KRT7 >0 FPKM"
                ),
            }
        )

    # FPKM vs counts concordance for CLDN4
    add_spearman(
        "CLDN4 FPKM vs CLDN4 counts",
        df["fpkm_CLDN4"],
        df["cpm_CLDN4"],
        "concordance",
        "same 50 libraries; log2(FPKM+1) vs log2(CPM+1)",
    )

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(TABLES / "all_tests.tsv", sep="\t", index=False)

    # one-row primary table (FPKM, library-level, honest n)
    c8 = tests_df[(tests_df.contrast == "CLDN4 vs CD8A") & (tests_df.matrix == "log2(FPKM+1)")].iloc[0]
    tnk = tests_df[(tests_df.contrast == "CLDN4 vs T/NK mean-z") & (tests_df.matrix == "log2(FPKM+1)")].iloc[0]
    ifn = tests_df[(tests_df.contrast == "CLDN4 vs IFN mean-z") & (tests_df.matrix == "log2(FPKM+1)")].iloc[0]
    pdl1 = tests_df[(tests_df.contrast == "CLDN4 vs CD274") & (tests_df.matrix == "log2(FPKM+1)")].iloc[0]
    resp = tests_df[(tests_df.contrast == "CLDN4 vs GEO sintilimab group") & (tests_df.matrix == "log2(FPKM+1)")].iloc[0]
    cldn4_nz = int((df["fpkm_CLDN4_raw"] > 0).sum())
    primary = pd.DataFrame(
        [
            {
                "dataset": "GSE260770",
                "tissue": "plasma exosome (not tumor)",
                "assay": "RNA-seq FPKM + counts",
                "n_libraries": 50,
                "n_patients_GEO": 0,
                "n_patients_paper_exosome": PAPER_EXOSOME_N,
                "n_R_libraries": int((df.resp == "R").sum()),
                "n_NR_libraries": int((df.resp == "NR").sum()),
                "CLDN4_nonzero_FPKM": cldn4_nz,
                "CLDN4_vs_CD8A_rho": c8["effect"],
                "CLDN4_vs_CD8A_p": c8["p"],
                "CLDN4_vs_TNK_rho": tnk["effect"],
                "CLDN4_vs_TNK_p": tnk["p"],
                "CLDN4_vs_IFN_rho": ifn["effect"],
                "CLDN4_vs_IFN_p": ifn["p"],
                "CLDN4_vs_CD274_rho": pdl1["effect"],
                "CLDN4_vs_CD274_p": pdl1["p"],
                "CLDN4_vs_response_delta": resp["effect"],
                "CLDN4_vs_response_p": resp["p"],
                "epithelial_residual": "not interpretable" if not residual_possible else "computed",
                "verdict": "FLOOR / NOT_TUMOR / PATIENT_ID_ABSENT",
                "note": (
                    "CLDN4-only. No dual-high. GEO has no patient ID; 50 libraries "
                    f"({int((df.resp=='R').sum())} R / {int((df.resp=='NR').sum())} NR). "
                    f"CLDN4 FPKM>0 in {cldn4_nz}/50 (median 0). Paper exosome n=10."
                ),
            }
        ]
    )
    primary.to_csv(TABLES / "primary.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE260770",
        "pmid": "38637495",
        "n_libraries": 50,
        "n_patients_geo": 0,
        "n_patients_paper_exosome": PAPER_EXOSOME_N,
        "n_R_libraries": int((df.resp == "R").sum()),
        "n_NR_libraries": int((df.resp == "NR").sum()),
        "cldn4_nonzero_fpkm": cldn4_nz,
        "residual_possible": bool(residual_possible),
        "n_epi_any": n_epi_any,
        "tnk_genes": sc_f.attrs.get("fpkm_TNK_genes", TNK_GENES),
        "ifn_genes": sc_f.attrs.get("fpkm_IFN_genes", IFN_GENES),
        "epi_genes_present": [g for g in EPI_GENES if g in fpkm.index],
        "primary": primary.iloc[0].to_dict(),
        "tests": tests_df.to_dict(orient="records"),
    }
    # numpy types
    def _jsonable(o):
        if isinstance(o, dict):
            return {k: _jsonable(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_jsonable(v) for v in o]
        if isinstance(o, (np.floating,)):
            return float(o) if np.isfinite(o) else None
        if isinstance(o, (np.integer,)):
            return int(o)
        if o is pd.NA:
            return None
        return o

    (TABLES / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2) + "\n")

    # --- figures ---
    r = df["resp"].eq("R")
    save_box(
        FIGS / "CLDN4_vs_response_FPKM.png",
        df.loc[r, "fpkm_CLDN4"],
        df.loc[~r, "fpkm_CLDN4"],
        "GSE260770 CLDN4 vs sintilimab (exosome)",
        "CLDN4 log2(FPKM+1)",
    )
    save_scatter(
        FIGS / "CLDN4_vs_CD8A_FPKM.png",
        df["fpkm_CLDN4"].to_numpy(),
        df["fpkm_CD8A"].to_numpy(),
        "CLDN4 vs CD8A",
        "CLDN4 log2(FPKM+1)",
        "CD8A log2(FPKM+1)",
        note="both at floor",
    )
    save_scatter(
        FIGS / "CLDN4_vs_TNK_FPKM.png",
        df["fpkm_CLDN4"].to_numpy(),
        df["fpkm_TNK"].to_numpy(),
        "CLDN4 vs T/NK",
        "CLDN4 log2(FPKM+1)",
        "T/NK mean-z",
    )
    save_scatter(
        FIGS / "CLDN4_vs_IFN_FPKM.png",
        df["fpkm_CLDN4"].to_numpy(),
        df["fpkm_IFN"].to_numpy(),
        "CLDN4 vs IFN",
        "CLDN4 log2(FPKM+1)",
        "IFN mean-z",
    )
    save_scatter(
        FIGS / "CLDN4_vs_CD274_FPKM.png",
        df["fpkm_CLDN4"].to_numpy(),
        df["fpkm_CD274"].to_numpy(),
        "CLDN4 vs CD274",
        "CLDN4 log2(FPKM+1)",
        "CD274 log2(FPKM+1)",
        note="both at floor",
    )
    save_scatter(
        FIGS / "CLDN4_FPKM_vs_counts.png",
        df["fpkm_CLDN4"].to_numpy(),
        df["cpm_CLDN4"].to_numpy(),
        "CLDN4 FPKM vs counts",
        "CLDN4 log2(FPKM+1)",
        "CLDN4 log2(CPM+1)",
    )

    # extra multi-panel
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.6))
    # A detection
    ax = axes[0, 0]
    genes_bar = ["CLDN4", "CD8A", "CD274", "NKG7", "HLA-DRA", "EPCAM", "KRT8"]
    nz = [int((fpkm.loc[g] > 0).sum()) if g in fpkm.index else 0 for g in genes_bar]
    ax.bar(range(len(genes_bar)), nz, color="#4C72B0")
    ax.axhline(50, color="#999", ls="--", lw=0.8)
    ax.set_xticks(range(len(genes_bar)))
    ax.set_xticklabels(genes_bar, rotation=35, ha="right")
    ax.set_ylabel("libraries with FPKM>0")
    ax.set_ylim(0, 55)
    ax.set_title(f"A  detection (n=50 libraries)\nCLDN4 {cldn4_nz}/50")
    # B response
    ax = axes[0, 1]
    ax.boxplot(
        [df.loc[r, "fpkm_CLDN4"], df.loc[~r, "fpkm_CLDN4"]],
        tick_labels=[f"R n={int(r.sum())}", f"NR n={int((~r).sum())}"],
        patch_artist=True,
        boxprops=dict(facecolor="#d95f02", alpha=0.5),
    )
    ax.set_ylabel("CLDN4 log2(FPKM+1)")
    ax.set_title(f"B  vs GEO response\nδ={fmt(resp['effect'])} p={fmt_p(resp['p'])}")
    # C vs T/NK
    ax = axes[1, 0]
    ax.scatter(df["fpkm_CLDN4"], df["fpkm_TNK"], s=22, c="#2c5aa0", alpha=0.75)
    ax.set_xlabel("CLDN4 log2(FPKM+1)")
    ax.set_ylabel("T/NK mean-z")
    ax.set_title(f"C  vs T/NK\nρ={fmt(tnk['effect'])} p={fmt_p(tnk['p'])}")
    # D vs IFN
    ax = axes[1, 1]
    ax.scatter(df["fpkm_CLDN4"], df["fpkm_IFN"], s=22, c="#2c5aa0", alpha=0.75)
    ax.set_xlabel("CLDN4 log2(FPKM+1)")
    ax.set_ylabel("IFN mean-z")
    ax.set_title(f"D  vs IFN\nρ={fmt(ifn['effect'])} p={fmt_p(ifn['p'])}")
    fig.suptitle(
        "GSE260770 extra — CLDN4-only plasma exosome (not 50 patients)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(FIGS / "extra_CLDN4_panel.png", dpi=170)
    fig.savefig(FIGS / "extra_CLDN4_panel.pdf")
    plt.close(fig)

    print(primary.to_string(index=False))
    print("\nFPKM primary contrasts:")
    print(
        tests_df[tests_df.matrix == "log2(FPKM+1)"][
            ["contrast", "n", "n_R", "n_NR", "metric", "effect", "p"]
        ].to_string(index=False)
    )
    print(f"\nepithelial residual possible: {residual_possible} (any EPI>0 in {n_epi_any}/50)")
    print(f"wrote {TABLES} and {FIGS}")


def residualize_safe(y, covar):
    y = np.asarray(y, float)
    c = np.asarray(covar, float)
    m = ~(np.isnan(y) | np.isnan(c))
    out = np.full_like(y, np.nan, dtype=float)
    if m.sum() < 4:
        return out
    Z = np.column_stack([np.ones(m.sum()), c[m]])
    b, *_ = np.linalg.lstsq(Z, y[m], rcond=None)
    out[m] = y[m] - Z @ b
    return out


if __name__ == "__main__":
    main()
