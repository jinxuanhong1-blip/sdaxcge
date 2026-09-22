#!/usr/bin/env python3
"""GSE131907 max-effect sweep for TJ delta and TACSTD2-high vs T depletion.

Winning rules are fixed here. They are not revised after seeing TACSTD2 output.

Honest unit: one number per patient. 208,506 cells are never the inferential n.
Sample-level Spearman p-values are saved as descriptive and cannot win.

TJ winner: among prespecified modules (>=3 genes, not the epithelial control),
the largest mean paired patient-level delta of a log1p(CP10k) score
(TACSTD2-high minus TACSTD2-low malignant/epithelial cells), requiring
mean delta > 0, median delta > 0, and n_patients >= 8.
The score scale is locked to natural-log log1p(count / library * 10000).

T winner: among patient-level Spearman tests with rho < 0 and n >= 8, the
smallest two-sided p for a TACSTD2 summary versus a T or CD8 fraction.
A second family is Q4 vs Q1 Mann-Whitney (Q4 lower than Q1, both arms >= 4).
T/NK is a reference outcome and cannot win. Naive sample-level p cannot win.
Searched p-values are descriptive, not a confirmatory test.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results"
TAB = OUT / "tables"
FIG = OUT / "figures"
SUB = Path("/tmp/gse131907/maxeffect_subset")

TUMOR_ORIGINS = ["tLung", "tL/B", "mLN", "PE", "mBrain"]
SITES = {
    "all": ["tLung", "tL/B", "mLN", "PE", "mBrain"],
    "no_mBrain": ["tLung", "tL/B", "mLN", "PE"],
    "tLung": ["tLung"],
    "met_lung": ["tL/B", "mLN", "PE"],
    "mBrain": ["mBrain"],
    "mLN": ["mLN"],
    "tLB": ["tL/B"],
}
CD8_SUB = {"Exhausted CD8+ T", "Naive CD8+ T", "Cytotoxic CD8+ T", "CD8 low T"}
CYTO_SUB = {"Cytotoxic CD8+ T"}
EXH_SUB = {"Exhausted CD8+ T"}
IMMUNE_TYPES = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
}
GATES = ("author", "ts", "broad", "epi")
MIN_N = 8
SPLITS = ("median", "q4q1", "pos")
MIN_ARMS = (10, 20)
SCORES = ("cell_mean", "pseudobulk")
POOLS = ("patient", "sample_then_patient")

# Prespecified TJ modules. CTRL_* is the epithelial coexpression control and
# cannot win the TJ contest.
MODULE_GENES = {
    "pre_core": [
        "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "TJP2", "TJP3",
        "F11R", "CGN", "MARVELD2", "CRB3",
    ],
    "pre_core_no_cldn4": [
        "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1", "TJP2", "TJP3",
        "F11R", "CGN", "MARVELD2", "CRB3",
    ],
    "claudins": ["CLDN1", "CLDN3", "CLDN4", "CLDN7"],
    "claudins_no_cldn4": ["CLDN1", "CLDN3", "CLDN7"],
    "structural": ["CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "CGN", "TJP1", "MARVELD2"],
    "structural_no_cldn4": ["CLDN3", "CLDN7", "OCLN", "F11R", "CGN", "TJP1", "MARVELD2"],
    "cldn_tjp": ["CLDN3", "CLDN4", "CLDN7", "TJP1", "TJP2", "TJP3"],
    "junction_adhesion": ["F11R", "CGN", "CGNL1", "OCLN", "MARVELD2", "CRB3", "TJP1"],
    "polarity": ["PARD3", "PARD6A", "PARD6B", "CRB3", "MPDZ", "INADL", "PATJ"],
    "CTRL_epi": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "CTRL_krt": ["KRT8", "KRT18", "KRT19"],
}

T_OUTCOMES = (
    "T_frac",
    "T_of_immune",
    "CD8_frac",
    "CD8_of_immune",
    "cyto_frac",
    "exh_frac",
)
REF_OUTCOMES = ("TNK_frac",)
PREDICTORS = ("tac_pct", "tac_pct2", "tac_mean")
DENOMS = ("gate_samples", "all_tumor_in_site")
MIN_MALS = (20, 50)

# Expected from the locked Seurat patient analysis (PR #541), author-malignant
# CLDN4 %pos vs T/NK, n_malignant >= 20.
CALIB_SAMPLE_RHO = -0.522077922077922
CALIB_PATIENT_RHO = -0.477922077922078


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return np.nan, np.nan, n
    xx, yy = x[m], y[m]
    if np.unique(xx).size < 2 or np.unique(yy).size < 2:
        return np.nan, np.nan, n
    res = stats.spearmanr(xx, yy, alternative="two-sided")
    return float(res.statistic), float(res.pvalue), n


def wilcoxon_p(deltas: np.ndarray) -> float:
    d = deltas[np.isfinite(deltas)]
    d = d[d != 0]
    if d.size < 1:
        return np.nan
    try:
        res = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox", method="auto")
    except ValueError:
        res = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox", method="asymptotic")
    return float(res.pvalue)


def mannwhitney(high: np.ndarray, low: np.ndarray) -> tuple[float, float]:
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    if high.size < 3 or low.size < 3:
        return np.nan, np.nan
    try:
        res = stats.mannwhitneyu(high, low, alternative="two-sided", method="auto")
    except ValueError:
        res = stats.mannwhitneyu(high, low, alternative="two-sided", method="asymptotic")
    n4, n1 = high.size, low.size
    # Rank-biserial from the Mann-Whitney U, matching the Seurat script.
    r_rb = (2.0 * float(res.statistic)) / (n4 * n1) - 1.0
    return r_rb, float(res.pvalue)


def assign_quartiles(x: np.ndarray) -> np.ndarray | None:
    if np.isfinite(x).sum() < 8:
        return None
    s = pd.Series(x)
    try:
        q = pd.qcut(s.rank(method="average"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    except ValueError:
        return None
    out = q.astype(str).to_numpy()
    if not {"Q1", "Q4"} <= set(out):
        return None
    return out


def parse_series(path: Path) -> pd.DataFrame:
    titles = patients = stages = origins = None
    with __import__("gzip").open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][len("!Sample_") :]
            vals = [v.strip('"') for v in parts[1:]]
            if key == "title":
                titles = vals
            elif key == "characteristics_ch1":
                if vals[0].startswith("patient id:"):
                    patients = [v.split(": ", 1)[1] for v in vals]
                elif vals[0].startswith("tumor stage:"):
                    stages = [v.split(": ", 1)[1] for v in vals]
                elif vals[0].startswith("tissue origin"):
                    origins = [v.split(": ", 1)[1] for v in vals]
    if not (titles and patients and stages and origins):
        raise SystemExit("series matrix missing sample characteristics")
    return pd.DataFrame(
        {
            "Sample": titles,
            "patient_id": patients,
            "tumor_stage": stages,
            "origin_geo": origins,
        }
    )


def load_kegg(present: set[str]) -> list[str]:
    genes = []
    for line in (DATA / "kegg_tj_no_cldn4.txt").read_text().splitlines():
        g = line.strip()
        if g and not g.startswith("#") and g in present and g not in ("TACSTD2", "CLDN4"):
            genes.append(g)
    return genes


def group_indices(ids: np.ndarray, mask: np.ndarray) -> list[tuple[str, np.ndarray]]:
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    labels = ids[idx]
    order = np.argsort(labels, kind="mergesort")
    idx = idx[order]
    labels = labels[order]
    cuts = np.flatnonzero(labels[1:] != labels[:-1]) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, idx.size]
    return [(labels[s], idx[s:e]) for s, e in zip(starts, ends)]


def split_arms(tac_log: np.ndarray, tac_umi: np.ndarray, mode: str, min_arm: int) -> tuple[np.ndarray, np.ndarray] | None:
    n = tac_log.size
    if mode == "median":
        if n < 2 * min_arm:
            return None
        order = np.argsort(tac_log, kind="mergesort")
        mid = n // 2
        if n % 2 == 0:
            low, high = order[:mid], order[mid:]
        else:
            low, high = order[:mid], order[mid + 1 :]
    elif mode == "q4q1":
        n1 = n // 4
        if n1 < min_arm:
            return None
        order = np.argsort(tac_log, kind="mergesort")
        low, high = order[:n1], order[-n1:]
    elif mode == "pos":
        high = np.flatnonzero(tac_umi > 0)
        low = np.flatnonzero(tac_umi <= 0)
        if high.size < min_arm or low.size < min_arm:
            return None
    else:
        raise ValueError(mode)
    if np.median(tac_log[high]) <= np.median(tac_log[low]) and np.mean(tac_log[high]) <= np.mean(tac_log[low]):
        return None
    return low, high


def arm_score(
    local_low: np.ndarray,
    local_high: np.ndarray,
    cells: np.ndarray,
    score: np.ndarray,
    gene_rows: np.ndarray | None,
    lib: np.ndarray,
    kind: str,
) -> float:
    if kind == "cell_mean":
        return float(score[cells[local_high]].mean() - score[cells[local_low]].mean())
    if gene_rows is None:
        raise ValueError("pseudobulk needs counts")
    def pb(local: np.ndarray) -> float:
        gi = cells[local]
        libsum = float(lib[gi].sum())
        if libsum <= 0:
            return np.nan
        tot = gene_rows[:, gi].sum(axis=1)
        return float(np.log1p(tot / libsum * 10000.0).mean())
    hi, lo = pb(local_high), pb(local_low)
    if not np.isfinite(hi) or not np.isfinite(lo):
        return np.nan
    return hi - lo


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.3e}"
    return f"{p:.4g}"


def fmt(x: float, d: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{d}f}"


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    genes = (SUB / "genes.txt").read_text().splitlines()
    gene_index = {g: i for i, g in enumerate(genes)}
    counts = np.load(SUB / "counts.npy")
    lib = np.load(SUB / "libsize.npy")
    barcodes = (SUB / "barcodes.txt").read_text().splitlines()
    if counts.shape != (len(genes), len(barcodes)):
        raise SystemExit(f"counts shape {counts.shape} != genes x barcodes")
    if lib.shape[0] != len(barcodes):
        raise SystemExit("libsize misaligned")
    if int((lib <= 0).sum()):
        raise SystemExit("nonpositive library size")

    lognorm = np.log1p(counts / lib[None, :] * 10000.0).astype(np.float32)

    ann = pd.read_csv(
        "/tmp/gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
        na_values=["NA"],
        keep_default_na=True,
    )
    series = parse_series(Path("/tmp/gse131907/GSE131907_series_matrix.txt.gz"))
    ann = ann.set_index("Index")
    if list(ann.index) != barcodes and not ann.index.isin(barcodes).all():
        missing = set(barcodes) - set(ann.index)
        raise SystemExit(f"annotation missing {len(missing)} barcodes")
    ann = ann.loc[barcodes]
    if ann["Sample"].isna().any():
        raise SystemExit("annotation Sample missing")
    smap = series.set_index("Sample")
    if not set(ann["Sample"].unique()) <= set(smap.index):
        raise SystemExit("annotation samples missing from series matrix")
    if not (ann["Sample_Origin"] == ann["Sample"].map(smap["origin_geo"])).all():
        raise SystemExit("annotation origin disagrees with GEO series matrix")

    patient = ann["Sample"].map(smap["patient_id"]).to_numpy()
    sample = ann["Sample"].to_numpy()
    origin = ann["Sample_Origin"].to_numpy()
    stage = ann["Sample"].map(smap["tumor_stage"]).to_numpy()
    cell_type = ann["Cell_type"].fillna("").to_numpy()
    subtype = ann["Cell_subtype"].fillna("").to_numpy()
    if patient.dtype.kind != "O":
        patient = patient.astype(str)

    n_cells = len(barcodes)
    if n_cells != 208506:
        raise SystemExit(f"n cells {n_cells}")
    author = subtype == "Malignant cells"
    ts = np.isin(subtype, ["tS1", "tS2", "tS3"])
    if int(author.sum()) != 24784:
        raise SystemExit(f"author malignant n={int(author.sum())}")
    if int(ts.sum()) != 6352:
        raise SystemExit(f"tS n={int(ts.sum())}")
    if int((author & (origin == "nLung")).sum()) != 0:
        raise SystemExit("nLung has author-malignant cells")
    if int((author & (origin == "tLung")).sum()) != 0:
        raise SystemExit("tLung has author-malignant cells; gate mix-up")
    if int((ts & (origin != "tLung")).sum()) != 0:
        raise SystemExit("tS cells outside tLung")

    is_tumor = np.isin(origin, TUMOR_ORIGINS)
    is_t = cell_type == "T lymphocytes"
    is_nk = cell_type == "NK cells"
    # A few author CD8 subtypes are Cell_type NK (231 cells). CD8 counts stay
    # inside T lymphocytes so they are not also part of the NK fraction.
    cd8_sub = np.isin(subtype, list(CD8_SUB))
    n_cd8_labeled_nk = int((cd8_sub & is_nk).sum())
    is_cd8 = cd8_sub & is_t
    is_cyto = np.isin(subtype, list(CYTO_SUB)) & is_t
    is_exh = np.isin(subtype, list(EXH_SUB)) & is_t
    is_immune = np.isin(cell_type, list(IMMUNE_TYPES))
    if n_cd8_labeled_nk != 231:
        raise SystemExit(f"unexpected CD8-labeled NK cells: {n_cd8_labeled_nk}")
    write_tsv(
        TAB / "n_honest.tsv",
        [
            {"piece": "cells_in_matrix", "n": int(n_cells), "note": "not an inferential n"},
            {"piece": "author_malignant_cells", "n": int(author.sum()), "note": "mBrain, tL/B, mLN only"},
            {"piece": "tS_cells", "n": int(ts.sum()), "note": "tLung only"},
            {"piece": "cd8_subtype_and_T", "n": int(is_cd8.sum()), "note": "CD8 subtypes inside T lymphocytes"},
            {"piece": "cd8_subtype_labeled_NK", "n": n_cd8_labeled_nk, "note": "excluded from CD8 fraction"},
        ],
    )

    tac = counts[gene_index["TACSTD2"]]
    tac_log = lognorm[gene_index["TACSTD2"]]
    cldn4 = counts[gene_index["CLDN4"]]

    present = set(genes)
    modules: dict[str, list[str]] = {}
    for name, glist in MODULE_GENES.items():
        use = [g for g in glist if g in present and g != "TACSTD2"]
        modules[name] = use
    modules["kegg_broad_no_cldn4"] = load_kegg(present)
    if len(modules["pre_core"]) < 8:
        raise SystemExit(f"pre_core missing genes: {modules['pre_core']}")
    if "CLDN4" not in modules["pre_core"]:
        raise SystemExit("CLDN4 missing from pre_core")
    if "TACSTD2" in modules["kegg_broad_no_cldn4"] or "CLDN4" in modules["kegg_broad_no_cldn4"]:
        raise SystemExit("TACSTD2 or CLDN4 leaked into the broad TJ set")

    # ------------------------------------------------------------------
    # Calibration: reproduce locked CLDN4 %pos vs T/NK (PR #541)
    # ------------------------------------------------------------------
    print("calibration", flush=True)
    cal_rows = []
    # sample rows
    samp_records = []
    for (pid, sid), idx in _group2(patient, sample, np.ones(n_cells, dtype=bool)):
        o = origin[idx][0]
        st = stage[idx][0]
        n = idx.size
        n_author = int(author[idx].sum())
        n_ts = int(ts[idx].sum())
        n_epi = int(((cell_type[idx] == "Epithelial cells") & is_tumor[idx]).sum())
        n_broad = n_author + n_ts
        rec = {
            "patient": pid,
            "sample": sid,
            "origin": o,
            "stage": st,
            "n_cells": n,
            "n_T": int(is_t[idx].sum()),
            "n_CD8": int(is_cd8[idx].sum()),
            "n_cyto": int(is_cyto[idx].sum()),
            "n_exh": int(is_exh[idx].sum()),
            "n_NK": int(is_nk[idx].sum()),
            "n_immune": int(is_immune[idx].sum()),
            "author_n": n_author,
            "author_tac_pos": int(((tac[idx] > 0) & author[idx]).sum()),
            "author_tac_ge2": int(((tac[idx] >= 2) & author[idx]).sum()),
            "author_tac_logsum": float(tac_log[idx][author[idx]].sum()) if n_author else 0.0,
            "author_cldn_pos": int(((cldn4[idx] > 0) & author[idx]).sum()),
            "ts_n": n_ts,
            "ts_tac_pos": int(((tac[idx] > 0) & ts[idx]).sum()),
            "ts_tac_ge2": int(((tac[idx] >= 2) & ts[idx]).sum()),
            "ts_tac_logsum": float(tac_log[idx][ts[idx]].sum()) if n_ts else 0.0,
            "broad_n": n_broad,
            "broad_tac_pos": int(((tac[idx] > 0) & (author[idx] | ts[idx])).sum()),
            "broad_tac_ge2": int(((tac[idx] >= 2) & (author[idx] | ts[idx])).sum()),
            "broad_tac_logsum": float(tac_log[idx][author[idx] | ts[idx]].sum()) if n_broad else 0.0,
            "epi_n": n_epi,
            "epi_tac_pos": int(((tac[idx] > 0) & (cell_type[idx] == "Epithelial cells")).sum()) if o in TUMOR_ORIGINS else int(((tac[idx] > 0) & (cell_type[idx] == "Epithelial cells") & is_tumor[idx]).sum()),
            "epi_tac_ge2": int(((tac[idx] >= 2) & (cell_type[idx] == "Epithelial cells") & is_tumor[idx]).sum()),
            "epi_tac_logsum": float(tac_log[idx][(cell_type[idx] == "Epithelial cells") & is_tumor[idx]].sum()) if n_epi else 0.0,
        }
        # epi counts must use tumor epithelial only. Recompute explicitly.
        epi_m = (cell_type[idx] == "Epithelial cells") & np.isin(origin[idx], TUMOR_ORIGINS)
        rec["epi_n"] = int(epi_m.sum())
        rec["epi_tac_pos"] = int((tac[idx][epi_m] > 0).sum()) if rec["epi_n"] else 0
        rec["epi_tac_ge2"] = int((tac[idx][epi_m] >= 2).sum()) if rec["epi_n"] else 0
        rec["epi_tac_logsum"] = float(tac_log[idx][epi_m].sum()) if rec["epi_n"] else 0.0
        samp_records.append(rec)
        cal_rows.append(rec)
    samples = pd.DataFrame(samp_records)
    tumor_samples = samples[samples["origin"].isin(TUMOR_ORIGINS)].copy()

    def frac_from(df: pd.DataFrame, num: str, den: str) -> np.ndarray:
        d = df[den].to_numpy(dtype=float)
        out = np.full(len(df), np.nan)
        ok = d > 0
        out[ok] = df[num].to_numpy(dtype=float)[ok] / d[ok]
        return out

    el_s = tumor_samples[tumor_samples["author_n"] >= 20].copy()
    el_s["cldn_pct"] = el_s["author_cldn_pos"] / el_s["author_n"]
    el_s["tnk"] = (el_s["n_T"] + el_s["n_NK"]) / el_s["n_cells"]
    rho_s, p_s, n_s = spearman(el_s["cldn_pct"].to_numpy(), el_s["tnk"].to_numpy())
    # patient pool of every tumor-origin sample
    pat = tumor_samples.groupby("patient", as_index=False).sum(numeric_only=True)
    # stage: one per patient
    stage_map = tumor_samples.groupby("patient")["stage"].agg(lambda s: sorted(set(s))[0] if len(set(s)) == 1 else "MIXED")
    pat["stage"] = pat["patient"].map(stage_map)
    el_p = pat[pat["author_n"] >= 20].copy()
    el_p["cldn_pct"] = el_p["author_cldn_pos"] / el_p["author_n"]
    el_p["tnk"] = (el_p["n_T"] + el_p["n_NK"]) / el_p["n_cells"]
    rho_p, p_p, n_p = spearman(el_p["cldn_pct"].to_numpy(), el_p["tnk"].to_numpy())
    print(f"  sample rho {rho_s:.6f} p {p_s:.4g} n {n_s}", flush=True)
    print(f"  patient rho {rho_p:.6f} p {p_p:.4g} n {n_p}", flush=True)
    if n_s != 21 or n_p != 21:
        raise SystemExit(f"calibration n sample/patient {n_s}/{n_p}, expected 21/21")
    if abs(rho_s - CALIB_SAMPLE_RHO) > 0.01 or abs(rho_p - CALIB_PATIENT_RHO) > 0.01:
        raise SystemExit(
            f"calibration rho mismatch sample {rho_s} vs {CALIB_SAMPLE_RHO}; "
            f"patient {rho_p} vs {CALIB_PATIENT_RHO}"
        )
    write_tsv(
        TAB / "calibration_cldn4_tnk.tsv",
        [
            {"unit": "sample_authorMalignant", "n": n_s, "rho_pct": rho_s, "p": p_s, "locked_rho": CALIB_SAMPLE_RHO},
            {"unit": "patient_authorMalignant_all_tumor_samples", "n": n_p, "rho_pct": rho_p, "p": p_p, "locked_rho": CALIB_PATIENT_RHO},
        ],
    )
    el_s.sort_values("patient").to_csv(TAB / "calibration_samples.tsv", sep="\t", index=False)
    el_p.sort_values("patient").to_csv(TAB / "calibration_patients.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    # TJ paired deltas
    # ------------------------------------------------------------------
    print("TJ sweep", flush=True)
    score_of = {}
    gene_rows_of = {}
    for name, glist in modules.items():
        ix = np.array([gene_index[g] for g in glist], dtype=int)
        score_of[name] = lognorm[ix].mean(axis=0).astype(np.float32)
        gene_rows_of[name] = counts[ix]

    def cell_mean_deltas(gi_hi: np.ndarray, gi_lo: np.ndarray) -> dict[str, float]:
        return {name: float(vec[gi_hi].mean() - vec[gi_lo].mean()) for name, vec in score_of.items()}

    def pb_from_sums(sum_hi: np.ndarray, sum_lo: np.ndarray, lib_hi: float, lib_lo: float, ix: np.ndarray) -> float:
        if lib_hi <= 0 or lib_lo <= 0:
            return np.nan
        hi = float(np.log1p(sum_hi[ix] / lib_hi * 10000.0).mean())
        lo = float(np.log1p(sum_lo[ix] / lib_lo * 10000.0).mean())
        return hi - lo

    gate_masks = {
        "author": author & is_tumor,
        "ts": ts & is_tumor,
        "broad": (author | ts) & is_tumor,
        "epi": (cell_type == "Epithelial cells") & is_tumor,
    }
    tj_rows: list[dict] = []
    # Store patient-level deltas for panels we will plot: filled after the sweep
    # by recomputing the chosen keys.
    for gate in GATES:
        base = gate_masks[gate]
        for site, origins_keep in SITES.items():
            site_mask = base & np.isin(origin, origins_keep)
            if int(site_mask.sum()) < 40:
                continue
            by_patient = group_indices(patient, site_mask)
            by_sample = group_indices(sample, site_mask)
            print(f"  {gate} {site} cells={int(site_mask.sum())} patients={len(by_patient)}", flush=True)
            # Arms do not depend on the module. Cache them, and cache per-arm
            # library size plus the sum of every extracted gene.
            arm_patient: dict[tuple, list] = {}
            for split in SPLITS:
                for min_arm in MIN_ARMS:
                    packed = []
                    for pid, cells in by_patient:
                        arms = split_arms(tac_log[cells], tac[cells], split, min_arm)
                        if arms is None:
                            continue
                        low, high = arms
                        gi_hi = cells[high]
                        gi_lo = cells[low]
                        packed.append(
                            (
                                str(pid),
                                cell_mean_deltas(gi_hi, gi_lo),
                                int(high.size),
                                int(low.size),
                                float(lib[gi_hi].sum()),
                                float(lib[gi_lo].sum()),
                                counts[:, gi_hi].sum(axis=1),
                                counts[:, gi_lo].sum(axis=1),
                            )
                        )
                    arm_patient[(split, min_arm)] = packed
            arm_sample: dict[tuple, dict] = {}
            for split in SPLITS:
                for min_arm in MIN_ARMS:
                    bags: dict[str, list] = {}
                    for _sid, cells in by_sample:
                        arms = split_arms(tac_log[cells], tac[cells], split, min_arm)
                        if arms is None:
                            continue
                        low, high = arms
                        pid = str(patient[cells[0]])
                        gi_hi = cells[high]
                        gi_lo = cells[low]
                        bags.setdefault(pid, []).append(
                            (
                                int(high.size),
                                int(low.size),
                                float(lib[gi_hi].sum()),
                                float(lib[gi_lo].sum()),
                                counts[:, gi_hi].sum(axis=1),
                                counts[:, gi_lo].sum(axis=1),
                                cells[high],
                                cells[low],
                            )
                        )
                    arm_sample[(split, min_arm)] = bags
            for module, score in score_of.items():
                ix = np.array([gene_index[g] for g in modules[module]], dtype=int)
                n_genes = int(ix.size)
                is_ctrl = module.startswith("CTRL_")
                for split in SPLITS:
                    for min_arm in MIN_ARMS:
                        for kind in SCORES:
                            for pool in POOLS:
                                deltas = []
                                n_hi = []
                                n_lo = []
                                if pool == "patient":
                                    for pid, cell_delta, nhi, nlo, lib_hi, lib_lo, sum_hi, sum_lo in arm_patient[(split, min_arm)]:
                                        if kind == "cell_mean":
                                            dlt = cell_delta[module]
                                        else:
                                            dlt = pb_from_sums(sum_hi, sum_lo, lib_hi, lib_lo, ix)
                                        if not np.isfinite(dlt):
                                            continue
                                        deltas.append(float(dlt))
                                        n_hi.append(nhi)
                                        n_lo.append(nlo)
                                else:
                                    for pid, items in arm_sample[(split, min_arm)].items():
                                        ds = []
                                        nhi = nlo = 0
                                        for item in items:
                                            if kind == "cell_mean":
                                                dlt = float(score[item[6]].mean() - score[item[7]].mean())
                                            else:
                                                dlt = pb_from_sums(item[4], item[5], item[2], item[3], ix)
                                            if np.isfinite(dlt):
                                                ds.append(float(dlt))
                                                nhi += item[0]
                                                nlo += item[1]
                                        if not ds:
                                            continue
                                        deltas.append(float(np.mean(ds)))
                                        n_hi.append(nhi)
                                        n_lo.append(nlo)
                                arr = np.asarray(deltas, dtype=float)
                                n_pat = int(arr.size)
                                mean_d = float(arr.mean()) if n_pat else np.nan
                                med_d = float(np.median(arr)) if n_pat else np.nan
                                frac_pos = float(np.mean(arr > 0)) if n_pat else np.nan
                                p_w = wilcoxon_p(arr) if n_pat else np.nan
                                eligible = (
                                    (not is_ctrl)
                                    and n_genes >= 3
                                    and n_pat >= MIN_N
                                    and np.isfinite(mean_d)
                                    and mean_d > 0
                                    and np.isfinite(med_d)
                                    and med_d > 0
                                )
                                tj_rows.append(
                                    {
                                        "gate": gate,
                                        "site": site,
                                        "pool": pool,
                                        "split": split,
                                        "min_arm": min_arm,
                                        "score": kind,
                                        "module": module,
                                        "n_genes": n_genes,
                                        "n_patients": n_pat,
                                        "mean_delta": mean_d,
                                        "median_delta": med_d,
                                        "frac_pos": frac_pos,
                                        "wilcoxon_p": p_w,
                                        "mean_n_high": float(np.mean(n_hi)) if n_hi else np.nan,
                                        "mean_n_low": float(np.mean(n_lo)) if n_lo else np.nan,
                                        "eligible": bool(eligible),
                                        "control": bool(is_ctrl),
                                    }
                                )
    tj = pd.DataFrame(tj_rows)
    tj.to_csv(TAB / "tj_sweep.tsv", sep="\t", index=False)
    elig = tj[tj["eligible"]].copy()
    if elig.empty:
        raise SystemExit("no eligible TJ panel")
    elig = elig.sort_values(
        ["mean_delta", "n_patients", "wilcoxon_p"],
        ascending=[False, False, True],
    )
    winner = elig.iloc[0].to_dict()
    print(
        f"  TJ winner {winner['module']} {winner['gate']} {winner['site']} "
        f"{winner['split']} {winner['score']} meanΔ={winner['mean_delta']:.4f} "
        f"n={winner['n_patients']} p={winner['wilcoxon_p']}",
        flush=True,
    )

    prespec_keys = [
        ("author", "all", "patient", "median", 10, "cell_mean", "pre_core"),
        ("author", "all", "patient", "q4q1", 10, "cell_mean", "pre_core"),
        ("author", "all", "patient", "median", 10, "cell_mean", "pre_core_no_cldn4"),
        ("ts", "tLung", "patient", "median", 10, "cell_mean", "pre_core"),
        ("ts", "tLung", "patient", "q4q1", 10, "cell_mean", "pre_core"),
        ("broad", "all", "patient", "median", 10, "cell_mean", "pre_core"),
    ]

    def lookup_tj(key: tuple) -> dict:
        g, site, pool, split, min_arm, kind, module = key
        hit = tj[
            (tj.gate == g)
            & (tj.site == site)
            & (tj.pool == pool)
            & (tj.split == split)
            & (tj.min_arm == min_arm)
            & (tj.score == kind)
            & (tj.module == module)
        ]
        if hit.empty:
            raise SystemExit(f"missing prespecified TJ row {key}")
        return hit.iloc[0].to_dict()

    prespec_tj = [lookup_tj(k) for k in prespec_keys]

    def same_setting(module: str, template: dict) -> dict:
        key = (
            template["gate"],
            template["site"],
            template["pool"],
            template["split"],
            int(template["min_arm"]),
            template["score"],
            module,
        )
        return lookup_tj(key)

    ctrl_at_winner = same_setting("CTRL_epi", winner)
    krt_at_winner = same_setting("CTRL_krt", winner)

    # Patient-level delta vectors for the figure / audit, plus residual on CTRL_epi.
    def patient_delta_table(spec: dict, tag: str) -> pd.DataFrame:
        gate = spec["gate"]
        site = spec["site"]
        module = spec["module"]
        split = spec["split"]
        min_arm = int(spec["min_arm"])
        kind = spec["score"]
        pool = spec["pool"]
        mask = gate_masks[gate] & np.isin(origin, SITES[site])
        score = score_of[module]
        grows = gene_rows_of[module]
        epi_score = score_of["CTRL_epi"]
        # residualize module score on epithelial control inside this mask
        xx = epi_score[mask]
        yy = score[mask]
        if np.std(xx) > 0:
            coef = np.polyfit(xx, yy, 1)
            resid = score - (coef[0] * epi_score + coef[1])
        else:
            resid = score.copy()
        rows = []
        if pool == "patient":
            groups = group_indices(patient, mask)
            bags: dict[str, list] = {}
            for pid, cells in groups:
                arms = split_arms(tac_log[cells], tac[cells], split, min_arm)
                if arms is None:
                    continue
                low, high = arms
                dlt = arm_score(low, high, cells, score, grows, lib, kind)
                dlt_r = arm_score(low, high, cells, resid, None if kind == "cell_mean" else grows, lib, "cell_mean" if kind == "cell_mean" else kind)
                # residual contrast is always cell-mean of residual scores
                dlt_r = float(resid[cells[high]].mean() - resid[cells[low]].mean())
                bags[pid] = [dlt, dlt_r, int(high.size), int(low.size), int(cells.size)]
            for pid, vals in bags.items():
                rows.append(
                    {
                        "tag": tag,
                        "patient": pid,
                        "delta": vals[0],
                        "delta_resid_epi": vals[1],
                        "n_high": vals[2],
                        "n_low": vals[3],
                        "n_gate_cells": vals[4],
                    }
                )
        else:
            acc: dict[str, list] = {}
            for sid, cells in group_indices(sample, mask):
                pid = str(patient[cells[0]])
                arms = split_arms(tac_log[cells], tac[cells], split, min_arm)
                if arms is None:
                    continue
                low, high = arms
                dlt = arm_score(low, high, cells, score, grows, lib, kind)
                dlt_r = float(resid[cells[high]].mean() - resid[cells[low]].mean())
                acc.setdefault(pid, []).append((dlt, dlt_r, int(high.size), int(low.size), int(cells.size)))
            for pid, items in acc.items():
                rows.append(
                    {
                        "tag": tag,
                        "patient": pid,
                        "delta": float(np.mean([it[0] for it in items])),
                        "delta_resid_epi": float(np.mean([it[1] for it in items])),
                        "n_high": int(np.sum([it[2] for it in items])),
                        "n_low": int(np.sum([it[3] for it in items])),
                        "n_gate_cells": int(np.sum([it[4] for it in items])),
                    }
                )
        return pd.DataFrame(rows)

    delta_frames = [
        patient_delta_table(prespec_tj[0], "prespec_author_median_pre_core"),
        patient_delta_table(prespec_tj[3], "prespec_ts_median_pre_core"),
        patient_delta_table(winner, "winner"),
    ]
    deltas_df = pd.concat(delta_frames, ignore_index=True)
    deltas_df.to_csv(TAB / "tj_patient_deltas.tsv", sep="\t", index=False)

    def resid_summary(tag: str) -> dict:
        sub = deltas_df[deltas_df.tag == tag]
        arr = sub["delta_resid_epi"].to_numpy(dtype=float)
        raw = sub["delta"].to_numpy(dtype=float)
        return {
            "tag": tag,
            "n_patients": int(len(sub)),
            "mean_delta": float(np.mean(raw)) if len(sub) else np.nan,
            "median_delta": float(np.median(raw)) if len(sub) else np.nan,
            "mean_delta_resid_epi": float(np.mean(arr)) if len(sub) else np.nan,
            "median_delta_resid_epi": float(np.median(arr)) if len(sub) else np.nan,
            "wilcoxon_p_resid": wilcoxon_p(arr) if len(sub) else np.nan,
            "frac_pos_resid": float(np.mean(arr > 0)) if len(sub) else np.nan,
        }

    resid_rows = [
        resid_summary("prespec_author_median_pre_core"),
        resid_summary("prespec_ts_median_pre_core"),
        resid_summary("winner"),
    ]
    # Recompute winner raw mean from the patient table and require it matches the sweep.
    wtab = deltas_df[deltas_df.tag == "winner"]["delta"].to_numpy(dtype=float)
    if abs(float(np.mean(wtab)) - float(winner["mean_delta"])) > 1e-6:
        raise SystemExit("winner delta table does not match sweep")
    write_tsv(TAB / "tj_residual_on_epi.tsv", resid_rows)

    # Single-gene maxima. Cell-mean, patient pool, min_arm 10, median and q4q1.
    print("single-gene TJ", flush=True)
    gene_menu = []
    for g in modules["pre_core"] + modules["structural"] + modules["polarity"] + ["CLDN18", "CDH1"]:
        if g in gene_index and g not in gene_menu and g != "TACSTD2":
            gene_menu.append(g)
    gene_rows_out = []
    for gname in gene_menu:
        gscore = lognorm[gene_index[gname]]
        grows = counts[gene_index[gname] : gene_index[gname] + 1]
        for gate, site in (("author", "all"), ("ts", "tLung"), ("broad", "all"), ("epi", "no_mBrain")):
            mask = gate_masks[gate] & np.isin(origin, SITES[site])
            for split in ("median", "q4q1"):
                ds = []
                for pid, cells in group_indices(patient, mask):
                    arms = split_arms(tac_log[cells], tac[cells], split, 10)
                    if arms is None:
                        continue
                    low, high = arms
                    dlt = float(gscore[cells[high]].mean() - gscore[cells[low]].mean())
                    ds.append(dlt)
                arr = np.asarray(ds, dtype=float)
                n_pat = int(arr.size)
                mean_d = float(arr.mean()) if n_pat else np.nan
                med_d = float(np.median(arr)) if n_pat else np.nan
                eligible = n_pat >= MIN_N and np.isfinite(mean_d) and mean_d > 0 and np.isfinite(med_d) and med_d > 0
                gene_rows_out.append(
                    {
                        "gene": gname,
                        "gate": gate,
                        "site": site,
                        "pool": "patient",
                        "split": split,
                        "min_arm": 10,
                        "score": "cell_mean",
                        "n_patients": n_pat,
                        "mean_delta": mean_d,
                        "median_delta": med_d,
                        "frac_pos": float(np.mean(arr > 0)) if n_pat else np.nan,
                        "wilcoxon_p": wilcoxon_p(arr) if n_pat else np.nan,
                        "eligible": bool(eligible),
                    }
                )
    gene_df = pd.DataFrame(gene_rows_out)
    gene_df.to_csv(TAB / "tj_gene_sweep.tsv", sep="\t", index=False)
    gene_el = gene_df[gene_df.eligible].sort_values(["mean_delta", "n_patients"], ascending=[False, False])
    gene_winner = gene_el.iloc[0].to_dict() if len(gene_el) else None

    # ------------------------------------------------------------------
    # TACSTD2 vs T depletion, patient unit
    # ------------------------------------------------------------------
    print("T depletion sweep", flush=True)
    gate_cols = {
        "author": ("author_n", "author_tac_pos", "author_tac_ge2", "author_tac_logsum"),
        "ts": ("ts_n", "ts_tac_pos", "ts_tac_ge2", "ts_tac_logsum"),
        "broad": ("broad_n", "broad_tac_pos", "broad_tac_ge2", "broad_tac_logsum"),
        "epi": ("epi_n", "epi_tac_pos", "epi_tac_ge2", "epi_tac_logsum"),
    }
    num_cols = {
        "T_frac": "n_T",
        "T_of_immune": "n_T",
        "CD8_frac": "n_CD8",
        "CD8_of_immune": "n_CD8",
        "cyto_frac": "n_cyto",
        "exh_frac": "n_exh",
        "TNK_frac": None,
    }
    t_rows = []
    patient_cache: dict[tuple, pd.DataFrame] = {}

    def build_patient_frame(gate: str, site: str, denom: str, min_mal: int) -> pd.DataFrame:
        key = (gate, site, denom, min_mal)
        if key in patient_cache:
            return patient_cache[key]
        n_col, pos_col, ge2_col, log_col = gate_cols[gate]
        site_df = tumor_samples[tumor_samples["origin"].isin(SITES[site])].copy()
        gate_tot = site_df.groupby("patient")[[n_col, pos_col, ge2_col, log_col]].sum()
        keep_patients = gate_tot.index[gate_tot[n_col] >= min_mal]
        if denom == "gate_samples":
            use = site_df[site_df[n_col] >= 1]
        elif denom == "all_tumor_in_site":
            use = site_df
        else:
            raise ValueError(denom)
        use = use[use["patient"].isin(keep_patients)]
        if use.empty:
            patient_cache[key] = pd.DataFrame()
            return patient_cache[key]
        agg = use.groupby("patient")[
            ["n_cells", "n_T", "n_CD8", "n_cyto", "n_exh", "n_NK", "n_immune", n_col]
        ].sum()
        pred = gate_tot.loc[agg.index]
        out = pd.DataFrame({"patient": agg.index.to_numpy()})
        out["n_gate"] = pred[n_col].to_numpy()
        out["tac_pct"] = pred[pos_col].to_numpy() / pred[n_col].to_numpy()
        out["tac_pct2"] = pred[ge2_col].to_numpy() / pred[n_col].to_numpy()
        out["tac_mean"] = pred[log_col].to_numpy() / pred[n_col].to_numpy()
        out["mal_frac"] = pred[n_col].to_numpy() / agg["n_cells"].to_numpy()
        out["n_cells"] = agg["n_cells"].to_numpy()
        out["T_frac"] = agg["n_T"].to_numpy() / agg["n_cells"].to_numpy()
        imm = agg["n_immune"].to_numpy().astype(float)
        out["T_of_immune"] = np.divide(agg["n_T"].to_numpy(), imm, out=np.full(len(agg), np.nan), where=imm > 0)
        out["CD8_frac"] = agg["n_CD8"].to_numpy() / agg["n_cells"].to_numpy()
        out["CD8_of_immune"] = np.divide(agg["n_CD8"].to_numpy(), imm, out=np.full(len(agg), np.nan), where=imm > 0)
        out["cyto_frac"] = agg["n_cyto"].to_numpy() / agg["n_cells"].to_numpy()
        out["exh_frac"] = agg["n_exh"].to_numpy() / agg["n_cells"].to_numpy()
        out["TNK_frac"] = (agg["n_T"].to_numpy() + agg["n_NK"].to_numpy()) / agg["n_cells"].to_numpy()
        out["stage"] = out["patient"].map(stage_map)
        patient_cache[key] = out.reset_index(drop=True)
        return patient_cache[key]

    for gate in GATES:
        for site in SITES:
            for denom in DENOMS:
                for min_mal in MIN_MALS:
                    frame = build_patient_frame(gate, site, denom, min_mal)
                    if frame.empty or len(frame) < 4:
                        continue
                    for pred_name in PREDICTORS:
                        x = frame[pred_name].to_numpy(dtype=float)
                        q = assign_quartiles(x)
                        for outcome in list(T_OUTCOMES) + list(REF_OUTCOMES):
                            y = frame[outcome].to_numpy(dtype=float)
                            rho, p, n = spearman(x, y)
                            row = {
                                "gate": gate,
                                "site": site,
                                "denom": denom,
                                "min_mal": min_mal,
                                "predictor": pred_name,
                                "outcome": outcome,
                                "test": "spearman",
                                "n_patients": n,
                                "rho": rho,
                                "p": p,
                                "n_q1": np.nan,
                                "n_q4": np.nan,
                                "median_q1": np.nan,
                                "median_q4": np.nan,
                                "delta_q4_minus_q1": np.nan,
                                "reference_only": outcome in REF_OUTCOMES,
                            }
                            can = (
                                outcome not in REF_OUTCOMES
                                and n >= MIN_N
                                and np.isfinite(rho)
                                and rho < 0
                                and np.isfinite(p)
                            )
                            row["eligible"] = bool(can)
                            t_rows.append(row)
                            if q is None:
                                continue
                            y1 = y[q == "Q1"]
                            y4 = y[q == "Q4"]
                            if np.isfinite(y1).sum() < 4 or np.isfinite(y4).sum() < 4:
                                continue
                            r_rb, p_mw = mannwhitney(y4[np.isfinite(y4)], y1[np.isfinite(y1)])
                            med1 = float(np.nanmedian(y1))
                            med4 = float(np.nanmedian(y4))
                            delta = med4 - med1
                            mw_row = {
                                "gate": gate,
                                "site": site,
                                "denom": denom,
                                "min_mal": min_mal,
                                "predictor": pred_name,
                                "outcome": outcome,
                                "test": "q4q1_mannwhitney",
                                "n_patients": n,
                                "rho": r_rb,
                                "p": p_mw,
                                "n_q1": int(np.isfinite(y1).sum()),
                                "n_q4": int(np.isfinite(y4).sum()),
                                "median_q1": med1,
                                "median_q4": med4,
                                "delta_q4_minus_q1": delta,
                                "reference_only": outcome in REF_OUTCOMES,
                                "eligible": bool(
                                    outcome not in REF_OUTCOMES
                                    and np.isfinite(p_mw)
                                    and delta < 0
                                    and np.isfinite(y1).sum() >= 4
                                    and np.isfinite(y4).sum() >= 4
                                ),
                            }
                            t_rows.append(mw_row)

    tdf = pd.DataFrame(t_rows)
    tdf.to_csv(TAB / "t_sweep.tsv", sep="\t", index=False)
    t_elig = tdf[tdf["eligible"]].copy()
    if t_elig.empty:
        t_winner = None
        t_winner_spear = None
        t_winner_mw = None
    else:
        t_winner_spear = (
            t_elig[t_elig.test == "spearman"].sort_values(["p", "n_patients"], ascending=[True, False]).head(1)
        )
        t_winner_mw = (
            t_elig[t_elig.test == "q4q1_mannwhitney"].sort_values(["p", "n_patients"], ascending=[True, False]).head(1)
        )
        t_winner_spear = t_winner_spear.iloc[0].to_dict() if len(t_winner_spear) else None
        t_winner_mw = t_winner_mw.iloc[0].to_dict() if len(t_winner_mw) else None
        cands = [c for c in (t_winner_spear, t_winner_mw) if c is not None]
        t_winner = sorted(cands, key=lambda r: (r["p"], -r["n_patients"]))[0]

    def t_prespec(gate, site, denom, min_mal, predictor, outcome) -> dict:
        frame = build_patient_frame(gate, site, denom, min_mal)
        x = frame[predictor].to_numpy(dtype=float) if len(frame) else np.array([])
        y = frame[outcome].to_numpy(dtype=float) if len(frame) else np.array([])
        rho, p, n = spearman(x, y)
        q = assign_quartiles(x) if len(frame) else None
        rec = {
            "gate": gate,
            "site": site,
            "denom": denom,
            "min_mal": min_mal,
            "predictor": predictor,
            "outcome": outcome,
            "n_patients": n,
            "rho": rho,
            "p_spearman": p,
            "p_q4q1": np.nan,
            "delta_q4_minus_q1": np.nan,
            "n_q1": np.nan,
            "n_q4": np.nan,
        }
        if q is not None:
            y1, y4 = y[q == "Q1"], y[q == "Q4"]
            if np.isfinite(y1).sum() >= 3 and np.isfinite(y4).sum() >= 3:
                _, p_mw = mannwhitney(y4[np.isfinite(y4)], y1[np.isfinite(y1)])
                rec["p_q4q1"] = p_mw
                rec["delta_q4_minus_q1"] = float(np.nanmedian(y4) - np.nanmedian(y1))
                rec["n_q1"] = int(np.isfinite(y1).sum())
                rec["n_q4"] = int(np.isfinite(y4).sum())
        return rec

    prespec_t = [
        t_prespec("author", "all", "gate_samples", 20, "tac_pct", "T_frac"),
        t_prespec("author", "all", "all_tumor_in_site", 20, "tac_pct", "T_frac"),
        t_prespec("author", "all", "gate_samples", 20, "tac_mean", "T_frac"),
        t_prespec("author", "all", "gate_samples", 20, "tac_pct", "CD8_frac"),
        t_prespec("ts", "tLung", "gate_samples", 20, "tac_pct", "T_frac"),
        t_prespec("ts", "tLung", "gate_samples", 20, "tac_mean", "T_frac"),
        t_prespec("ts", "tLung", "gate_samples", 20, "tac_pct", "CD8_frac"),
        t_prespec("broad", "all", "gate_samples", 20, "tac_pct", "T_frac"),
        t_prespec("epi", "all", "gate_samples", 20, "tac_pct", "T_frac"),
    ]
    write_tsv(TAB / "t_prespecified.tsv", prespec_t)

    def partial_spearman(x, y, z) -> tuple[float, float, int]:
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        n = int(m.sum())
        if n < 6:
            return np.nan, np.nan, n
        rx = stats.rankdata(x[m])
        ry = stats.rankdata(y[m])
        rz = stats.rankdata(z[m])
        def resid(a, b):
            design = np.column_stack([np.ones(len(b)), b])
            coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
            return a - design @ coef
        xr, yr = resid(rx, rz), resid(ry, rz)
        if np.std(xr) == 0 or np.std(yr) == 0:
            return np.nan, np.nan, n
        res = stats.pearsonr(xr, yr)
        return float(res.statistic), float(res.pvalue), n

    def panel_frame(spec: dict) -> pd.DataFrame:
        return build_patient_frame(spec["gate"], spec["site"], spec["denom"], int(spec["min_mal"]))

    robustness = []
    rng = np.random.default_rng(131907)
    if t_winner is not None:
        fr = panel_frame(t_winner)
        x = fr[t_winner["predictor"]].to_numpy(dtype=float)
        y = fr[t_winner["outcome"]].to_numpy(dtype=float)
        z = fr["mal_frac"].to_numpy(dtype=float)
        pr, pp, pn = partial_spearman(x, y, z)
        # permutation of this one panel, two-sided |rho|
        obs = float(t_winner["rho"]) if t_winner["test"] == "spearman" else np.nan
        p_perm = np.nan
        if t_winner["test"] == "spearman" and np.isfinite(obs):
            m = np.isfinite(x) & np.isfinite(y)
            xx, yy = x[m], y[m]
            exceed = 0
            n_perm = 4999
            for _ in range(n_perm):
                r, _ = stats.spearmanr(xx, rng.permutation(yy))
                if abs(float(r)) + 1e-15 >= abs(obs):
                    exceed += 1
            p_perm = (exceed + 1) / (n_perm + 1)
        fr.assign(panel="winner").to_csv(TAB / "t_winner_patients.tsv", sep="\t", index=False)
        loo_note = ""
        if t_winner["test"] == "spearman":
            y_w = fr[t_winner["outcome"]].to_numpy(dtype=float)
            x_w = fr[t_winner["predictor"]].to_numpy(dtype=float)
            n_w = fr["n_cells"].to_numpy(dtype=float)
            total_out = float(np.nansum(y_w * n_w))
            n_zero = int(np.sum(np.isfinite(y_w) & (y_w == 0)))
            loo_ps = []
            for i in range(len(fr)):
                keep = np.ones(len(fr), dtype=bool)
                keep[i] = False
                _, p_loo, _ = spearman(x_w[keep], y_w[keep])
                if np.isfinite(p_loo):
                    loo_ps.append(p_loo)
            loo_note = (
                f" Across these patients the outcome totals about {total_out:.0f} cells"
                f" ({n_zero} patients are zero)."
            )
            if loo_ps:
                loo_note += (
                    f" Leave-one-out Spearman p ranges from {fmt_p(float(min(loo_ps)))}"
                    f" to {fmt_p(float(max(loo_ps)))}."
                )
        robustness.append(
            {
                "panel": "winner",
                "test": t_winner["test"],
                "partial_rho_on_gate_fraction": pr,
                "partial_p": pp,
                "partial_n": pn,
                "permutation_p_this_panel": p_perm,
                "n_perm": 4999 if np.isfinite(p_perm) else 0,
            }
        )
    # partials for the two primary prespecified rows
    for label, spec in (
        ("prespec_author_pct_T", prespec_t[0]),
        ("prespec_ts_pct_T", prespec_t[4]),
    ):
        fr = panel_frame(
            {
                "gate": spec["gate"],
                "site": spec["site"],
                "denom": spec["denom"],
                "min_mal": spec["min_mal"],
            }
        )
        if fr.empty:
            continue
        pr, pp, pn = partial_spearman(
            fr[spec["predictor"]].to_numpy(dtype=float),
            fr[spec["outcome"]].to_numpy(dtype=float),
            fr["mal_frac"].to_numpy(dtype=float),
        )
        robustness.append(
            {
                "panel": label,
                "test": "spearman",
                "partial_rho_on_gate_fraction": pr,
                "partial_p": pp,
                "partial_n": pn,
                "permutation_p_this_panel": np.nan,
                "n_perm": 0,
            }
        )
        fr.assign(panel=label).to_csv(TAB / f"t_{label}_patients.tsv", sep="\t", index=False)
    write_tsv(TAB / "t_robustness.tsv", robustness)

    n_spear = int((tdf.test == "spearman").sum())
    n_mw = int((tdf.test == "q4q1_mannwhitney").sum())
    n_spear_elig = int(((tdf.test == "spearman") & tdf.eligible).sum())
    n_spear_p05 = int(((tdf.test == "spearman") & tdf.eligible & (tdf.p < 0.05)).sum())
    n_mw_p05 = int(((tdf.test == "q4q1_mannwhitney") & tdf.eligible & (tdf.p < 0.05)).sum())
    n_tj = int(len(tj))
    n_tj_elig = int(tj.eligible.sum())

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    print("figures", flush=True)
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), sharey=False)
    for ax, dfp, title, rho, p, n in (
        (axes[0], el_s, "Sample", rho_s, p_s, n_s),
        (axes[1], el_p, "Patient, tumor samples pooled", rho_p, p_p, n_p),
    ):
        ax.scatter(dfp["cldn_pct"] * 100, dfp["tnk"], s=36, c="#1b4f72", zorder=3)
        ax.set_xlabel("Author-malignant CLDN4 % positive")
        ax.set_ylabel("T/NK fraction")
        ax.set_title(f"{title}\nρ={rho:.3f}, p={fmt_p(p)}, n={n} patients" if "Patient" in title else f"{title}\nρ={rho:.3f}, p={fmt_p(p)}, n={n} samples")
    fig.suptitle("Calibration against the locked GSE131907 CLDN4 row", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "fig_calibration_cldn4_tnk.png", dpi=160)
    fig.savefig(FIG / "fig_calibration_cldn4_tnk.pdf")
    plt.close()

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    plot_tags = [
        ("prespec_author_median_pre_core", "Author malignant\npre-specified core"),
        ("prespec_ts_median_pre_core", "Primary tS\npre-specified core"),
        ("winner", "Largest mean Δ\n" + str(winner["module"])),
    ]
    # de-duplicate if winner is the prespecified author panel
    rng_plot = np.random.default_rng(1)
    for i, (tag, label) in enumerate(plot_tags):
        vals = deltas_df.loc[deltas_df.tag == tag, "delta"].to_numpy(dtype=float)
        jitter = rng_plot.uniform(-0.12, 0.12, size=vals.size)
        ax.scatter(np.full(vals.size, i) + jitter, vals, s=28, c="#1b4f72", zorder=3)
        if vals.size:
            ax.hlines(np.mean(vals), i - 0.25, i + 0.25, colors="#b03a2e", lw=2, zorder=4)
    ax.axhline(0, color="#666666", lw=1)
    ax.set_xticks(range(len(plot_tags)))
    ax.set_xticklabels([p[1] for p in plot_tags])
    ax.set_ylabel("Paired TJ Δ, log1p(CP10k)\nTACSTD2-high minus low")
    ax.set_title(
        f"Patient-level TJ Δ  ·  winner n={int(winner['n_patients'])} patients\n"
        "Red bar = mean. Cell counts are not the sample size."
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig_tj_patient_deltas.png", dpi=160)
    fig.savefig(FIG / "fig_tj_patient_deltas.pdf")
    plt.close()

    # Module bars at the winning setting, including controls.
    setting = tj[
        (tj.gate == winner["gate"])
        & (tj.site == winner["site"])
        & (tj.pool == winner["pool"])
        & (tj.split == winner["split"])
        & (tj.min_arm == winner["min_arm"])
        & (tj.score == winner["score"])
    ].copy()
    setting = setting[setting.n_patients >= MIN_N].sort_values("mean_delta", ascending=True)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    colors = ["#7f8c8d" if m.startswith("CTRL_") else "#1b4f72" for m in setting["module"]]
    ax.barh(setting["module"], setting["mean_delta"], color=colors)
    ax.axvline(0, color="#333333", lw=1)
    ax.set_xlabel("Mean paired patient Δ, log1p(CP10k)")
    ax.set_title(
        f"Same contrast as the maximum\n{winner['gate']}, {winner['site']}, {winner['split']}, {winner['score']}, n≥{MIN_N}"
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig_tj_modules_at_winner.png", dpi=160)
    fig.savefig(FIG / "fig_tj_modules_at_winner.pdf")
    plt.close()

    # T scatters: prespec author, prespec tS, winner
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
    scatter_specs = [
        (prespec_t[0], "Author malignant, pre-specified"),
        (prespec_t[4], "Primary tS, pre-specified"),
    ]
    if t_winner is not None and t_winner["test"] == "spearman":
        scatter_specs.append((t_winner, "Smallest depletion p"))
    elif t_winner is not None:
        # still show spearman winner if different
        if t_winner_spear is not None:
            scatter_specs.append((t_winner_spear, "Smallest Spearman p"))
        else:
            scatter_specs.append((prespec_t[2], "Author malignant, TACSTD2 mean"))
    else:
        scatter_specs.append((prespec_t[2], "Author malignant, TACSTD2 mean"))
    for ax, (spec, title) in zip(axes, scatter_specs):
        if "p_spearman" in spec:
            frame = build_patient_frame(spec["gate"], spec["site"], spec["denom"], int(spec["min_mal"]))
            xx = frame[spec["predictor"]].to_numpy(dtype=float)
            yy = frame[spec["outcome"]].to_numpy(dtype=float)
            rho, p, n = spec["rho"], spec["p_spearman"], spec["n_patients"]
            pred, outcome = spec["predictor"], spec["outcome"]
        else:
            frame = panel_frame(spec)
            xx = frame[spec["predictor"]].to_numpy(dtype=float)
            yy = frame[spec["outcome"]].to_numpy(dtype=float)
            rho, p, n = spec["rho"], spec["p"], spec["n_patients"]
            pred, outcome = spec["predictor"], spec["outcome"]
        ax.scatter(xx, yy, s=32, c="#1b4f72")
        ax.set_xlabel(pred)
        ax.set_ylabel(outcome.replace("_", " "))
        ax.set_title(f"{title}\nρ={fmt(rho)}, p={fmt_p(p)}, n={int(n)}")
    fig.suptitle("TACSTD2 versus the named T outcome, one point per patient", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "fig_tacstd2_vs_t.png", dpi=160)
    fig.savefig(FIG / "fig_tacstd2_vs_t.pdf")
    plt.close()

    # Honest n
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = [
        "Cells in matrix",
        "Calibration patients",
        "TJ winner patients",
        "T pre-specified\nauthor",
        "T pre-specified\nprimary tS",
    ]
    values = [
        208506,
        int(n_p),
        int(winner["n_patients"]),
        int(prespec_t[0]["n_patients"]),
        int(prespec_t[4]["n_patients"]),
    ]
    if t_winner is not None:
        labels.append("T depletion\nwinner")
        values.append(int(t_winner["n_patients"]))
    cols = ["#95a5a6", "#1b4f72", "#1b4f72", "#1b4f72", "#1b4f72"] + (["#b03a2e"] if t_winner is not None else [])
    ax.bar(range(len(values)), values, color=cols)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.set_title("Inferential n is patients, not cells")
    fig.tight_layout()
    fig.savefig(FIG / "fig_honest_n.png", dpi=160)
    fig.savefig(FIG / "fig_honest_n.pdf")
    plt.close()

    # ------------------------------------------------------------------
    # FINDING
    # ------------------------------------------------------------------
    def tj_line(r: dict) -> str:
        return (
            f"| {r['gate']} | {r['site']} | {r['pool']} | {r['split']} | {r['score']} | {r['module']} | "
            f"{int(r['n_genes'])} | {int(r['n_patients'])} | {fmt(float(r['mean_delta']))} | "
            f"{fmt(float(r['median_delta']))} | {fmt(float(r['frac_pos']))} | {fmt_p(float(r['wilcoxon_p']))} |"
        )

    def t_line(r: dict) -> str:
        return (
            f"| {r['gate']} | {r['site']} | {r['denom']} | {r['predictor']} | {r['outcome']} | "
            f"{int(r['n_patients'])} | {fmt(float(r['rho']))} | {fmt_p(float(r['p_spearman']))} | "
            f"{fmt(float(r['delta_q4_minus_q1']))} | {fmt_p(float(r['p_q4q1']))} |"
        )

    w_resid = resid_rows[2]
    top_tj = elig.drop_duplicates(
        subset=["gate", "site", "pool", "split", "score", "module", "n_patients", "mean_delta"]
    ).head(8)
    top_lines = "\n".join(tj_line(r) for _, r in top_tj.iterrows())
    pre_lines = "\n".join(tj_line(r) for r in prespec_tj)
    pre_t_lines = "\n".join(t_line(r) for r in prespec_t)

    if t_winner is None:
        t_head = "No patient-level panel with n≥8 had TACSTD2 higher and T (or CD8) lower."
        t_detail = ""
    else:
        tw = t_winner
        rob = robustness[0]
        bonf = 0.05 / max(n_spear, 1)
        t_head = (
            f"Smallest nominal depletion p is {fmt_p(float(tw['p']))} "
            f"({tw['test']}, {tw['gate']}, {tw['site']}, {tw['denom']}, "
            f"{tw['predictor']} vs {tw['outcome']}, n={int(tw['n_patients'])} patients"
            + (
                f", ρ={fmt(float(tw['rho']))}"
                if tw["test"] == "spearman"
                else f", Q4−Q1 median={fmt(float(tw['delta_q4_minus_q1']))}"
            )
            + ")."
        )
        perm_bit = ""
        if np.isfinite(rob["permutation_p_this_panel"]):
            perm_bit = f" Permutation p for this single panel is {fmt_p(float(rob['permutation_p_this_panel']))} (4999 shuffles)."
        partial_bit = (
            f" Partial Spearman on the gate-cell fraction is ρ={fmt(float(rob['partial_rho_on_gate_fraction']))}, "
            f"p={fmt_p(float(rob['partial_p']))} (n={int(rob['partial_n'])})."
        )
        hits = tdf[(tdf.test == "spearman") & tdf.eligible & (tdf.p < 0.05)]
        if len(hits):
            uniq = (
                hits.groupby(["site", "outcome", "predictor"], as_index=False)
                .agg(rho=("rho", "first"), n_patients=("n_patients", "first"), n_rows=("p", "size"))
            )
            uniq_txt = "; ".join(
                f"{r.outcome} at {r.site}, {r.predictor}, ρ={fmt(float(r.rho))}, n={int(r.n_patients)} ({int(r.n_rows)} duplicate rows)"
                for r in uniq.itertuples()
            )
        else:
            uniq_txt = "none"
        n_tfrac_p05 = int(
            (
                (tdf.test == "spearman")
                & (tdf.outcome == "T_frac")
                & (tdf.rho < 0)
                & (tdf.n_patients >= MIN_N)
                & (tdf.p < 0.05)
            ).sum()
        )
        partial_note = ""
        if np.isfinite(rob["partial_p"]) and float(rob["partial_p"]) >= 0.05:
            partial_note = " After the gate-fraction partial, that nominal p is no longer below 0.05."
        t_detail = (
            f"Spearman tests in the sweep: {n_spear}. "
            f"Depletion-direction and n≥8: {n_spear_elig}, of which {n_spear_p05} have nominal p<0.05. "
            f"Those p<0.05 rows collapse to: {uniq_txt}. "
            f"Total T fraction (T_frac) never reaches p<0.05 in this sweep (rows below 0.05: {n_tfrac_p05}). "
            f"Q4 vs Q1 tests: {n_mw}, of which {n_mw_p05} are eligible and p<0.05. "
            f"Bonferroni 0.05 line for the Spearman family: {fmt_p(bonf)}."
            f"{perm_bit}{partial_bit}{partial_note}{loo_note} "
            "The searched minimum is not a confirmatory p-value."
        )

    if gene_winner is None:
        gene_bit = "No single TJ gene had a positive paired delta at n≥8."
    else:
        gw = gene_winner
        gene_bit = (
            f"Largest single-gene mean Δ is {gw['gene']} "
            f"({gw['gate']}, {gw['site']}, {gw['split']}): "
            f"mean Δ={fmt(float(gw['mean_delta']))}, median Δ={fmt(float(gw['median_delta']))}, "
            f"n={int(gw['n_patients'])}, Wilcoxon p={fmt_p(float(gw['wilcoxon_p']))}. "
            "That row is not a multi-gene TJ module."
        )

    ctrl_sentence = (
        f"At the same contrast, the epithelial control EPCAM/KRT8/18/19 mean Δ is "
        f"{fmt(float(ctrl_at_winner['mean_delta']))} (n={int(ctrl_at_winner['n_patients'])}) and "
        f"KRT8/18/19 alone is {fmt(float(krt_at_winner['mean_delta']))}."
    )
    if float(ctrl_at_winner["mean_delta"]) >= float(winner["mean_delta"]):
        ctrl_sentence += " The control is at least as large as the winning TJ module, so this maximum is not TJ-specific."
    else:
        ctrl_sentence += " The winning TJ module is larger than that epithelial control on this contrast."

    resid_sentence = (
        f"After a cell-level linear residual of the winning score on EPCAM/KRT8/18/19, "
        f"the mean paired Δ is {fmt(float(w_resid['mean_delta_resid_epi']))} "
        f"(median {fmt(float(w_resid['median_delta_resid_epi']))}, "
        f"Wilcoxon p={fmt_p(float(w_resid['wilcoxon_p_resid']))}, n={int(w_resid['n_patients'])})."
    )

    pre0 = prespec_tj[0]
    pre_ts = prespec_tj[3]
    author_best = elig[elig.gate == "author"].iloc[0].to_dict()
    site_words = {
        "all": "all tumor origins",
        "no_mBrain": "tumor origins except brain metastasis",
        "tLung": "primary lung",
        "met_lung": "tL/B, mLN, and pleural effusion",
        "mBrain": "brain metastasis",
        "mLN": "metastatic lymph node",
        "tLB": "bronchus or EBUS tumor",
    }
    gate_words = {
        "author": "author Malignant cells",
        "ts": "primary tS1/tS2/tS3",
        "broad": "author malignant plus tS",
        "epi": "all tumor-site epithelial cells",
    }
    finding = f"""# FINDING — GSE131907 max effect, TJ Δ and TACSTD2 vs T

Kim et al., *Nat Commun* 2020, GSE131907. Public processed UMI matrix only (208,506 cells). The log2TPM text and EGA FASTQ were not used.

Inferential unit is the **patient**. Each TJ Δ is a within-patient contrast (TACSTD2-high minus TACSTD2-low cells in the gate). Each T association is one TACSTD2 summary and one T fraction per patient. Sample-level tests that repeat patients are not eligible to win. p-values from the sweep are descriptive.

The pipeline first reproduces the locked author-malignant CLDN4 % positive versus T/NK row (n=21): sample ρ={fmt(rho_s, 3)} (locked {fmt(CALIB_SAMPLE_RHO, 3)}), patient ρ={fmt(rho_p, 3)} (locked {fmt(CALIB_PATIENT_RHO, 3)}).

## TJ Δ

Score scale is locked: mean of log1p(UMI / full library × 10,000). TACSTD2 is not in any TJ score. A panel needs ≥3 genes, n≥8 patients, and both the mean and the median paired Δ > 0. Epithelial controls cannot win.

**Maximum mean paired Δ = {fmt(float(winner['mean_delta']))}** on `{winner['module']}` ({int(winner['n_genes'])} genes). Gate `{winner['gate']}` ({gate_words[str(winner['gate'])]}), site `{winner['site']}` ({site_words[str(winner['site'])]}), pool `{winner['pool']}`, split `{winner['split']}`, score `{winner['score']}`, min cells/arm {int(winner['min_arm'])}. n={int(winner['n_patients'])} patients. Median Δ={fmt(float(winner['median_delta']))}. Fraction of patients with Δ>0 = {fmt(float(winner['frac_pos']))}. Two-sided Wilcoxon p={fmt_p(float(winner['wilcoxon_p']))}.

Best author-malignant panel (metastatic cells only, not the broader epithelial gate): `{author_best['module']}`, site `{author_best['site']}` ({site_words[str(author_best['site'])]}), split `{author_best['split']}`, mean Δ={fmt(float(author_best['mean_delta']))}, n={int(author_best['n_patients'])}, Wilcoxon p={fmt_p(float(author_best['wilcoxon_p']))}. Author-malignant cells are absent from primary lung and from pleural effusion, so the author `no_mBrain` and author `met_lung` rows are the same cells.

{ctrl_sentence} {resid_sentence}

Pre-specified core (CLDN1/3/4/7, OCLN, TJP1/2/3, F11R, CGN, MARVELD2, CRB3), median split, cell-mean score, patient pool: author-malignant mean Δ={fmt(float(pre0['mean_delta']))} (n={int(pre0['n_patients'])}, p={fmt_p(float(pre0['wilcoxon_p']))}); primary tS mean Δ={fmt(float(pre_ts['mean_delta']))} (n={int(pre_ts['n_patients'])}, p={fmt_p(float(pre_ts['wilcoxon_p']))}).

Sweep size: {n_tj} module rows, {n_tj_elig} eligible. {gene_bit}

### Pre-specified TJ rows

| gate | site | pool | split | score | module | genes | n | mean Δ | median Δ | frac>0 | Wilcoxon p |
|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|
{pre_lines}

### Eight largest eligible mean Δ

| gate | site | pool | split | score | module | genes | n | mean Δ | median Δ | frac>0 | Wilcoxon p |
|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|
{top_lines}

## TACSTD2-high versus T depletion

Pre-specified primary: author-malignant cells, samples that contain those cells, TACSTD2 % (UMI>0) versus T-cell fraction, ≥20 malignant cells. ρ={fmt(float(prespec_t[0]['rho']))}, Spearman p={fmt_p(float(prespec_t[0]['p_spearman']))}, n={int(prespec_t[0]['n_patients'])}. Primary tS, same predictor versus T fraction: ρ={fmt(float(prespec_t[4]['rho']))}, p={fmt_p(float(prespec_t[4]['p_spearman']))}, n={int(prespec_t[4]['n_patients'])}.

{t_head} {t_detail}

T outcomes that can win: T fraction of the sample, T fraction of immune cells, CD8 fraction, CD8 fraction of immune cells, cytotoxic CD8 fraction, exhausted CD8 fraction. T/NK is reported in the sweep and cannot win. Denominators are either the samples that contain the gate cells, or every tumor sample in the site (including a capture with no gate cells). The second denominator is the patient pool used for the locked CLDN4 row.

### Pre-specified T rows

| gate | site | denominator | predictor | outcome | n | ρ | Spearman p | Q4−Q1 | Q4 vs Q1 p |
|---|---|---|---|---|---:|---:|---:|---:|---:|
{pre_t_lines}

## What this is not

- Not a cell-level Wilcoxon. The MPE and primary-tS cell-level TJ p-values in the paper-funnel pass are a different unit.
- Not a claim that TACSTD2 causes T-cell loss, and not a spatial exclusion test. This atlas is dissociated.
- Not a confirmatory p-value for the searched minimum.
- Not GSE148071, not the four-cohort pool, and not a mouse Tacstd2 result.
- Author `Malignant cells` are metastatic (mBrain, tL/B, mLN), not primary tLung. Primary tumor epithelium in this matrix is tS1/tS2/tS3.

## Files

- `results/tables/tj_sweep.tsv` — module sweep
- `results/tables/tj_patient_deltas.tsv` — one Δ per patient
- `results/tables/t_sweep.tsv` — T/CD8 sweep
- `results/tables/calibration_cldn4_tnk.tsv` — lock check
- `results/figures/fig_tj_patient_deltas.png`
- `results/figures/fig_tacstd2_vs_t.png`
- `results/figures/fig_honest_n.png`
"""
    (ROOT / "FINDING.md").write_text(finding)
    summary = {
        "calibration_sample_rho": rho_s,
        "calibration_patient_rho": rho_p,
        "tj_winner": {k: (float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int)) else v) for k, v in winner.items() if k != "eligible"},
        "tj_ctrl_epi_mean_delta": float(ctrl_at_winner["mean_delta"]),
        "tj_resid_mean_delta": float(w_resid["mean_delta_resid_epi"]),
        "tj_n_rows": n_tj,
        "tj_n_eligible": n_tj_elig,
        "gene_winner": None
        if gene_winner is None
        else {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in gene_winner.items()},
        "t_winner": None
        if t_winner is None
        else {k: (None if isinstance(v, float) and not np.isfinite(v) else (float(v) if isinstance(v, (np.floating, float)) else v)) for k, v in t_winner.items()},
        "t_n_spearman": n_spear,
        "t_n_spearman_p05_eligible": n_spear_p05,
        "prespec_author_T_rho": float(prespec_t[0]["rho"]),
        "prespec_author_T_p": float(prespec_t[0]["p_spearman"]),
        "prespec_ts_T_rho": float(prespec_t[4]["rho"]),
        "prespec_ts_T_p": float(prespec_t[4]["p_spearman"]),
        "n_cells_not_inferential": n_cells,
    }
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        return obj
    (OUT / "summary.json").write_text(json.dumps(_clean(summary), indent=2))
    print(json.dumps(_clean(summary), indent=2), flush=True)
    print("done", flush=True)


def _group2(patient: np.ndarray, sample: np.ndarray, mask: np.ndarray):
    idx = np.flatnonzero(mask)
    key = np.array([f"{patient[i]}\t{sample[i]}" for i in idx])
    order = np.argsort(key, kind="mergesort")
    idx = idx[order]
    key = key[order]
    cuts = np.flatnonzero(key[1:] != key[:-1]) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, idx.size]
    out = []
    for s, e in zip(starts, ends):
        pid, sid = key[s].split("\t")
        out.append(((pid, sid), idx[s:e]))
    return out


if __name__ == "__main__":
    main()
