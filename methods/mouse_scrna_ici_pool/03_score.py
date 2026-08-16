#!/usr/bin/env python3
"""Score Tacstd2/Cldn4 in public mouse-lung ICI scRNA (epithelial/tumor vs T/NK)."""
from __future__ import annotations

import gzip
import json
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "notes/mouse_scrna_ici_pool/raw/data"
OUT = ROOT / "results/mouse_scrna_ici_pool"
PANEL_PATH = Path(__file__).with_name("gene_panel.tsv")

# 10x Cell Ranger 3.0 mm10 (31,053 genes) — same table as GSE275877 features
MM10_31053 = DATA / "GSE275877_features.tsv.gz"


def load_panel():
    df = pd.read_csv(PANEL_PATH, sep="\t")
    roles = defaultdict(list)
    for _, r in df.iterrows():
        roles[r.role].append(r.symbol)
    return df["symbol"].tolist(), dict(roles)


def read_features(path: Path) -> list[str]:
    genes = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and not parts[1].isdigit():
                genes.append(parts[1])
            else:
                genes.append(parts[0])
    return genes


def norm_symbol(s: str) -> str:
    return s.replace(".", "-").upper()


ALIASES = {
    "Epcam": ["Tacstd1"],
    "Nkx2-1": ["Nkx2.1", "Ttf1"],
    "Klrb1c": ["Nk1.1"],
}


def gene_index(genes: list[str], wanted: list[str]) -> dict[str, int]:
    lut = {}
    for i, g in enumerate(genes):
        lut[norm_symbol(g)] = i
    out = {}
    for w in wanted:
        keys = [w] + ALIASES.get(w, [])
        for k in keys:
            i = lut.get(norm_symbol(k))
            if i is not None:
                out[w] = i
                break
    return out


def stream_mtx_panel(mtx_path: Path, idx: dict[str, int], n_genes: int, n_cells: int, min_umi: int | None):
    """One-pass MTX: per-cell UMI totals + panel counts. 10x is genes x cells, 1-indexed."""
    want = {i: name for name, i in idx.items()}
    umi = np.zeros(n_cells, dtype=np.int64)
    mat = {name: np.zeros(n_cells, dtype=np.float32) for name in idx}
    opener = gzip.open if str(mtx_path).endswith(".gz") else open
    with opener(mtx_path, "rt") as f:
        header_done = False
        for line in f:
            if line.startswith("%"):
                continue
            if not header_done:
                header_done = True
                continue
            a, b, c = line.split()
            gi = int(a) - 1
            ci = int(b) - 1
            v = float(c)
            if 0 <= ci < n_cells:
                umi[ci] += v
                if gi in want:
                    mat[want[gi]][ci] = v
    if min_umi is not None:
        keep = umi >= min_umi
    else:
        keep = umi > 0
    mat = {k: v[keep] for k, v in mat.items()}
    return mat, int(keep.sum()), int((~keep).sum())


def mtx_shape(path: Path) -> tuple[int, int, int]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            ng, nc, nz = map(int, line.split()[:3])
            return ng, nc, nz
    raise ValueError(path)


def extract_dense_tsv(path: Path, wanted: list[str]):
    """GSE157881/882: cells x genes, header gene symbols, already log-like."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        genes = header[1:]
        idx = gene_index(genes, wanted)
        cols = {name: [] for name in idx}
        n = 0
        for line in f:
            parts = line.rstrip("\n").split("\t")
            n += 1
            for name, j in idx.items():
                cols[name].append(float(parts[j + 1]))
    mat = {k: np.asarray(v, dtype=np.float32) for k, v in cols.items()}
    return mat, n, idx


def log1p_if_counts(mat: dict[str, np.ndarray], already_log: bool):
    if already_log:
        return mat
    return {k: np.log1p(v) for k, v in mat.items()}


def mean_score(mat, names):
    arrs = [mat[n] for n in names if n in mat]
    if not arrs:
        return None
    return np.mean(np.vstack(arrs), axis=0)


def _pos(mat, name, n):
    if name in mat:
        return mat[name] > 0
    return np.zeros(n, dtype=bool)


def classify(mat, roles, presort: str | None):
    """Epithelial/tumor = Epcam or Cdh1+Krt8 or SCLC Ascl1/Chga, not T-lineage.

    Do not use Sftpc/Scgb1a1 (ambient / normal AT2-club in whole-lung preps).
    Do not require Ptprc==0 (ambient CD45).
    """
    n = len(next(iter(mat.values())))
    if presort == "cd45neg":
        return np.array(["epithelial"] * n, dtype=object)
    tnk = (
        _pos(mat, "Cd3e", n)
        | _pos(mat, "Cd3d", n)
        | _pos(mat, "Cd8a", n)
        | _pos(mat, "Nkg7", n)
        | _pos(mat, "Ncr1", n)
    )
    if presort == "cd45pos":
        lab = np.array(["immune"] * n, dtype=object)
        lab[tnk] = "tnk"
        return lab
    if presort == "cd3":
        return np.array(["tnk"] * n, dtype=object)
    epi = (
        _pos(mat, "Epcam", n)
        | (_pos(mat, "Cdh1", n) & _pos(mat, "Krt8", n))
        | _pos(mat, "Ascl1", n)
        | _pos(mat, "Chga", n)
        | _pos(mat, "Insm1", n)
    )
    lab = np.array(["other"] * n, dtype=object)
    # Epcam/tumor markers win over ambient T-lineage UMIs
    lab[tnk] = "tnk"
    lab[epi] = "epithelial"
    return lab


def summarize_sample(mat, labels, already_log: bool):
    logmat = log1p_if_counts(mat, already_log)
    rows = []
    n = len(labels)
    for comp in ["epithelial", "tnk", "immune", "other"]:
        mask = labels == comp
        k = int(mask.sum())
        rec = {"compartment": comp, "n_cells": k, "frac": k / n if n else 0.0}
        for gene in ["Tacstd2", "Cldn4"]:
            if gene not in logmat:
                rec[f"{gene}_mean"] = np.nan
                rec[f"{gene}_pct"] = np.nan
                continue
            x = logmat[gene][mask] if k else np.array([])
            rec[f"{gene}_mean"] = float(x.mean()) if k else np.nan
            rec[f"{gene}_pct"] = float((mat[gene][mask] > 0).mean() * 100) if k else np.nan
        rows.append(rec)
    return rows


def welch_or_nan(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(t), float(p)


def mwu_or_nan(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(u), float(p)
    except ValueError:
        return np.nan, np.nan


SAMPLES = [
    # GSE157881 — HKP1 lung, RT not ICB; CD45± sorts
    dict(gse="GSE157881", gsm="GSM4777885", sample="0Gy_CD45neg", arm="control_0Gy",
         icb="no", sort="cd45neg", fmt="dense_tsv", already_log=True,
         path="GSM4777885_0Gy-CD45-minus_S15.txt.gz", leftover=False,
         note="HKP1 lung; radiation 0 Gy; CD45− epithelial/tumor sort"),
    dict(gse="GSE157881", gsm="GSM4777886", sample="0Gy_CD45pos", arm="control_0Gy",
         icb="no", sort="cd45pos", fmt="dense_tsv", already_log=True,
         path="GSM4777886_0Gy-CD45-plus_S14.txt.gz", leftover=False,
         note="HKP1 lung; radiation 0 Gy; CD45+ immune sort"),
    dict(gse="GSE157881", gsm="GSM4777887", sample="4Gy_CD45neg", arm="RT_4Gy",
         icb="no", sort="cd45neg", fmt="dense_tsv", already_log=True,
         path="GSM4777887_4Gy-CD45-minus_S17.txt.gz", leftover=False,
         note="HKP1 lung; 4 Gy RT; CD45−; ICB not in these libraries"),
    dict(gse="GSE157881", gsm="GSM4777888", sample="4Gy_CD45pos", arm="RT_4Gy",
         icb="no", sort="cd45pos", fmt="dense_tsv", already_log=True,
         path="GSM4777888_4Gy-CD45-plus_S16.txt.gz", leftover=False,
         note="HKP1 lung; 4 Gy RT; CD45+"),
    # GSE157882 — CD45+ only, club-cell depletion after 4 Gy
    dict(gse="GSE157882", gsm="GSM4777889", sample="PBS_CD45pos", arm="PBS_4Gy",
         icb="no", sort="cd45pos", fmt="dense_tsv", already_log=True,
         path="GSM4777889_PBS_CD45_plus_S6_2W.txt.gz", leftover=False,
         note="HKP1; 4 Gy + PBS; CD45+ only; no epithelial"),
    dict(gse="GSE157882", gsm="GSM4777890", sample="DT_CD45pos", arm="DT_club_deplete_4Gy",
         icb="no", sort="cd45pos", fmt="dense_tsv", already_log=True,
         path="GSM4777890_DT_CD45_plus_S8_2W.txt.gz", leftover=False,
         note="HKP1; 4 Gy + DT club-cell depletion; CD45+ only"),
    # GSE176091 — CA170/VISTA; CD45+ TILs only (2 vs 2)
    dict(gse="GSE176091", gsm="GSM5354872", sample="CD45_Control_1", arm="control",
         icb="control", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM5354872_CD45_Control_CD45_1_matrix.mtx.gz",
         features="GSM5354872_CD45_Control_CD45_1_features.tsv.gz", leftover=False,
         note="Vinyl-carbamate LUAD; CD45+ TILs; PBS"),
    dict(gse="GSE176091", gsm="GSM5354873", sample="CD45_Control_2", arm="control",
         icb="control", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM5354873_CD45_Control_CD45_2_matrix.mtx.gz",
         features="GSM5354873_CD45_Control_CD45_2_features.tsv.gz", leftover=False,
         note="CD45+ TILs; PBS replicate"),
    dict(gse="GSE176091", gsm="GSM5354876", sample="CD45_CA170_1", arm="CA170",
         icb="icb", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM5354876_CD45_CA170_1_matrix.mtx.gz",
         features="GSM5354876_CD45_CA170_1_features.tsv.gz", leftover=False,
         note="CA170 VISTA/PD-L1 antagonist; CD45+ TILs"),
    dict(gse="GSE176091", gsm="GSM5354877", sample="CD45_CA170_2", arm="CA170",
         icb="icb", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM5354877_CD45_CA170_2_matrix.mtx.gz",
         features="GSM5354877_CD45_CA170_2_features.tsv.gz", leftover=False,
         note="CA170 replicate; CD45+ TILs"),
    # GSE268525 — LLC1-sgLkb1 lung mets; ICI vs RT+ICI (no untreated)
    dict(gse="GSE268525", gsm="GSM8293873", sample="Flox_ICI", arm="ICI",
         icb="icb", sort=None, fmt="mtx", already_log=False,
         path="GSM8293873_NC_matrix.mtx.gz", features="GSM8293873_NC_features.tsv.gz",
         leftover=False, note="LLC1-sgLkb1 lung met; Sfrp2-flox; ICI only"),
    dict(gse="GSE268525", gsm="GSM8293874", sample="Flox_RI", arm="RT_ICI",
         icb="icb_combo", sort=None, fmt="mtx", already_log=False,
         path="GSM8293874_RI_matrix.mtx.gz", features="GSM8293874_RI_features.tsv.gz",
         leftover=False, note="LLC1-sgLkb1 lung met; Sfrp2-flox; RT+ICI"),
    dict(gse="GSE268525", gsm="GSM8293872", sample="CKO_RI", arm="Sfrp2cKO_RT_ICI",
         icb="icb_combo", sort=None, fmt="mtx", already_log=False,
         path="GSM8293872_CKORI_matrix.mtx.gz", features="GSM8293872_CKORI_features.tsv.gz",
         leftover=False, note="Sfrp2ΔCol1a2 + RT+ICI; not a no-ICI control"),
    # GSE303943 — CMT167R subcutaneous aPD-1-resistant; PKCi vs solvent (not ICB vs control)
    dict(gse="GSE303943", gsm="GSM9139390", sample="CMT167R_CON", arm="solvent",
         icb="no", sort=None, fmt="mtx", already_log=False,
         path="GSM9139390_CON_matrix.mtx.gz", features="GSM9139390_CON_features.tsv.gz",
         leftover=False, note="Subcutaneous CMT167R (lung line, aPD-1-resistant); solvent"),
    dict(gse="GSE303943", gsm="GSM9139391", sample="CMT167R_PKCi", arm="PKCi",
         icb="no", sort=None, fmt="mtx", already_log=False,
         path="GSM9139391_PKCi_matrix.mtx.gz", features="GSM9139391_PKCi_features.tsv.gz",
         leftover=False, note="Subcutaneous CMT167R; PKC inhibitor, not ICB"),
    # leftover: GSE133604 KP ± aPD-1 (mm10 31053)
    dict(gse="GSE133604", gsm="GSM3912860", sample="KP_Ctrl_IgG", arm="control",
         icb="control", sort=None, fmt="mtx", already_log=False, min_umi=200,
         path="GSM3912860_1_Ctrl_matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover KP lung; IgG; features borrowed from GSE275877 mm10-3.0 (31053)"),
    dict(gse="GSE133604", gsm="GSM3912861", sample="KP_Ctrl_aPD1", arm="aPD1",
         icb="icb", sort=None, fmt="mtx", already_log=False, min_umi=200,
         path="GSM3912861_2_CtrlplusPD1_matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover KP lung; anti-PD-1"),
    dict(gse="GSE133604", gsm="GSM3912862", sample="KP_Asf1aKO_IgG", arm="Asf1aKO",
         icb="control", sort=None, fmt="mtx", already_log=False, min_umi=200,
         path="GSM3912862_3_ko_matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover KP Asf1a-KO; IgG"),
    dict(gse="GSE133604", gsm="GSM3912863", sample="KP_Asf1aKO_aPD1", arm="Asf1aKO_aPD1",
         icb="icb", sort=None, fmt="mtx", already_log=False, min_umi=200,
         path="GSM3912863_4_KOplusPD1_matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover KP Asf1a-KO; anti-PD-1"),
    # leftover: GSE129297 SCLC ± aPD-1 (unfiltered 737k barcodes)
    dict(gse="GSE129297", gsm="GSM3704194", sample="SCLC_Ctrl", arm="control",
         icb="control", sort=None, fmt="mtx", already_log=False, min_umi=500,
         path="GSM3704194_Ctrl.matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover SCLC GEMM lung; control; unfiltered 10x, min UMI 500"),
    dict(gse="GSE129297", gsm="GSM3704195", sample="SCLC_aPD1", arm="aPD1",
         icb="icb", sort=None, fmt="mtx", already_log=False, min_umi=500,
         path="GSM3704195_PD1.matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover SCLC GEMM; anti-PD-1"),
    dict(gse="GSE129297", gsm="GSM3704196", sample="SCLC_YKL", arm="YKL",
         icb="no", sort=None, fmt="mtx", already_log=False, min_umi=500,
         path="GSM3704196_YKL.matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover SCLC; CDK7i YKL, not ICB"),
    dict(gse="GSE129297", gsm="GSM3704197", sample="SCLC_combo", arm="YKL_aPD1",
         icb="icb_combo", sort=None, fmt="mtx", already_log=False, min_umi=500,
         path="GSM3704197_combo.matrix.mtx.gz", features="GSE275877_features.tsv.gz",
         leftover=True, note="Leftover SCLC; YKL + anti-PD-1"),
    # leftover immune-only with aPD-1 vs control
    dict(gse="GSE232730", gsm="GSM7372335", sample="KPL_CD45_Ctrl", arm="control",
         icb="control", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM7372335_Control_matrix.mtx.gz", features="GSM7372335_Control_features.tsv.gz",
         leftover=True, note="Leftover KPL-3M subcutaneous; CD45+; no epithelial"),
    dict(gse="GSE232730", gsm="GSM7372336", sample="KPL_CD45_aPD1", arm="aPD1",
         icb="icb", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM7372336_PD1_matrix.mtx.gz", features="GSM7372336_PD1_features.tsv.gz",
         leftover=True, note="Leftover KPL-3M; CD45+; anti-PD-1"),
    dict(gse="GSE222158", gsm="GSM6915764", sample="CD45_Ctl", arm="control",
         icb="control", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM6915764_A_matrix.mtx.gz", features="GSM6915764_A_features.tsv.gz",
         leftover=True, note="Leftover FVB lung cancer; CD45+ control"),
    dict(gse="GSE222158", gsm="GSM6915765", sample="CD45_PD1", arm="aPD1",
         icb="icb", sort="cd45pos", fmt="mtx", already_log=False,
         path="GSM6915765_B_matrix.mtx.gz", features="GSM6915765_B_features.tsv.gz",
         leftover=True, note="Leftover; CD45+ anti-PD-1"),
]


def score_one(cfg, wanted, roles):
    path = DATA / cfg["path"]
    rec = {k: cfg.get(k) for k in
           ["gse", "gsm", "sample", "arm", "icb", "sort", "leftover", "note"]}
    if not path.exists():
        rec["status"] = "missing_file"
        rec["detail"] = str(path.name)
        return rec, []
    if cfg["fmt"] == "dense_tsv":
        mat, n, idx = extract_dense_tsv(path, wanted)
        rec["n_cells_raw"] = n
        rec["n_cells_used"] = n
        rec["genes_found"] = ",".join(sorted(idx))
        rec["Tacstd2_in_matrix"] = "Tacstd2" in idx
        rec["Cldn4_in_matrix"] = "Cldn4" in idx
    else:
        feat = DATA / cfg["features"]
        if not feat.exists():
            rec["status"] = "missing_features"
            return rec, []
        genes = read_features(feat)
        ng, nc, nz = mtx_shape(path)
        rec["mtx_shape"] = f"{ng}x{nc}"
        if ng != len(genes):
            rec["status"] = "features_mismatch"
            rec["detail"] = f"mtx {ng} genes vs features {len(genes)}"
            return rec, []
        idx = gene_index(genes, wanted)
        rec["genes_found"] = ",".join(sorted(idx))
        rec["Tacstd2_in_matrix"] = "Tacstd2" in idx
        rec["Cldn4_in_matrix"] = "Cldn4" in idx
        min_umi = cfg.get("min_umi")
        mat, n_use, n_drop = stream_mtx_panel(path, idx, ng, nc, min_umi)
        rec["n_cells_raw"] = nc
        rec["n_cells_used"] = n_use
        rec["n_cells_dropped"] = n_drop
    labels = classify(mat, roles, cfg.get("sort"))
    rec["n_epithelial"] = int((labels == "epithelial").sum())
    rec["n_tnk"] = int((labels == "tnk").sum())
    rec["pct_Epcam"] = float((_pos(mat, "Epcam", len(labels))).mean() * 100)
    rec["pct_Cd3e"] = float((_pos(mat, "Cd3e", len(labels))).mean() * 100)
    rec["status"] = "ok"
    rows = summarize_sample(mat, labels, cfg.get("already_log", False))
    for r in rows:
        r.update({k: rec[k] for k in ["gse", "gsm", "sample", "arm", "icb", "leftover"]})
    return rec, rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wanted, roles = load_panel()
    sample_rows = []
    cell_rows = []
    for cfg in SAMPLES:
        print("SCORE", cfg["gse"], cfg["sample"], flush=True)
        rec, rows = score_one(cfg, wanted, roles)
        sample_rows.append(rec)
        cell_rows.extend(rows)
        print(" ", rec.get("status"), rec.get("n_epithelial"), rec.get("n_tnk"),
              "Tacstd2", rec.get("Tacstd2_in_matrix"), flush=True)

    samp = pd.DataFrame(sample_rows)
    cells = pd.DataFrame(cell_rows)
    samp.to_csv(OUT / "sample_inventory.tsv", sep="\t", index=False)
    cells.to_csv(OUT / "compartment_scores.tsv", sep="\t", index=False)

    # contrasts
    contrasts = []

    def add_contrast(gse, name, a_samples, b_samples, compartment, gene, extra=""):
        sub = cells[(cells.gse == gse) & (cells.compartment == compartment)]
        a = sub[sub["sample"].isin(a_samples)][f"{gene}_mean"].astype(float).tolist()
        b = sub[sub["sample"].isin(b_samples)][f"{gene}_mean"].astype(float).tolist()
        na, nb = len(a), len(b)
        if na == 0 or nb == 0:
            return
        mean_a = float(np.nanmean(a)) if na else np.nan
        mean_b = float(np.nanmean(b)) if nb else np.nan
        delta = mean_b - mean_a
        t, p_t = welch_or_nan(b, a)
        u, p_u = mwu_or_nan(b, a)
        contrasts.append({
            "gse": gse, "contrast": name, "compartment": compartment, "gene": gene,
            "n_a": na, "n_b": nb, "mean_a": mean_a, "mean_b": mean_b,
            "delta_log": delta, "direction": ("up" if delta > 0 else "down" if delta < 0 else "tie"),
            "welch_p": p_t, "mwu_p": p_u, "note": extra,
        })

    # GSE157881 RT vs 0Gy in CD45− (epithelial) and CD45+ T/NK
    add_contrast("GSE157881", "4Gy_vs_0Gy", ["0Gy_CD45neg"], ["4Gy_CD45neg"],
                 "epithelial", "Tacstd2", "RT not ICB; n=1 vs 1 — p not computed")
    add_contrast("GSE157881", "4Gy_vs_0Gy", ["0Gy_CD45neg"], ["4Gy_CD45neg"],
                 "epithelial", "Cldn4", "RT not ICB; n=1 vs 1")
    add_contrast("GSE157881", "4Gy_vs_0Gy_CD45pos", ["0Gy_CD45pos"], ["4Gy_CD45pos"],
                 "tnk", "Tacstd2", "immune sort")
    add_contrast("GSE157881", "4Gy_vs_0Gy_CD45pos", ["0Gy_CD45pos"], ["4Gy_CD45pos"],
                 "tnk", "Cldn4", "immune sort")

    def add_sorted_pair(gse, name, epi_sample, tnk_sample, gene, extra=""):
        epi = cells[(cells.gse == gse) & (cells.sample == epi_sample) & (cells.compartment == "epithelial")]
        tnk = cells[(cells.gse == gse) & (cells.sample == tnk_sample) & (cells.compartment == "tnk")]
        if epi.empty or tnk.empty:
            return
        ea = float(epi[f"{gene}_mean"].iloc[0])
        tb = float(tnk[f"{gene}_mean"].iloc[0])
        contrasts.append({
            "gse": gse, "contrast": name, "compartment": "epithelial_minus_tnk",
            "gene": gene, "n_a": 1, "n_b": 1, "mean_a": tb, "mean_b": ea,
            "delta_log": ea - tb,
            "direction": "up" if ea > tb else "down" if ea < tb else "tie",
            "welch_p": np.nan, "mwu_p": np.nan,
            "note": extra,
        })

    add_sorted_pair("GSE157881", "CD45neg_vs_CD45pos_0Gy", "0Gy_CD45neg", "0Gy_CD45pos",
                    "Tacstd2", "author CD45 sort; n=1 library each")
    add_sorted_pair("GSE157881", "CD45neg_vs_CD45pos_0Gy", "0Gy_CD45neg", "0Gy_CD45pos",
                    "Cldn4", "author CD45 sort; n=1 library each")
    add_sorted_pair("GSE157881", "CD45neg_vs_CD45pos_4Gy", "4Gy_CD45neg", "4Gy_CD45pos",
                    "Tacstd2", "author CD45 sort; n=1 library each")
    add_sorted_pair("GSE157881", "CD45neg_vs_CD45pos_4Gy", "4Gy_CD45neg", "4Gy_CD45pos",
                    "Cldn4", "author CD45 sort; n=1 library each")

    add_contrast("GSE176091", "CA170_vs_control", ["CD45_Control_1", "CD45_Control_2"],
                 ["CD45_CA170_1", "CD45_CA170_2"], "tnk", "Tacstd2",
                 "CD45+ TILs; epithelial ABSENT")
    add_contrast("GSE176091", "CA170_vs_control", ["CD45_Control_1", "CD45_Control_2"],
                 ["CD45_CA170_1", "CD45_CA170_2"], "tnk", "Cldn4",
                 "CD45+ TILs; epithelial ABSENT")

    add_contrast("GSE268525", "RTICI_vs_ICI", ["Flox_ICI"], ["Flox_RI"],
                 "epithelial", "Tacstd2", "no untreated arm; n=1 vs 1")
    add_contrast("GSE268525", "RTICI_vs_ICI", ["Flox_ICI"], ["Flox_RI"],
                 "epithelial", "Cldn4", "no untreated arm; n=1 vs 1")
    add_contrast("GSE268525", "RTICI_vs_ICI", ["Flox_ICI"], ["Flox_RI"],
                 "tnk", "Tacstd2", "")
    add_contrast("GSE268525", "RTICI_vs_ICI", ["Flox_ICI"], ["Flox_RI"],
                 "tnk", "Cldn4", "")

    add_contrast("GSE303943", "PKCi_vs_solvent", ["CMT167R_CON"], ["CMT167R_PKCi"],
                 "epithelial", "Tacstd2", "subcutaneous; not ICB vs control; n=1 vs 1")
    add_contrast("GSE303943", "PKCi_vs_solvent", ["CMT167R_CON"], ["CMT167R_PKCi"],
                 "epithelial", "Cldn4", "subcutaneous; n=1 vs 1")
    add_contrast("GSE303943", "PKCi_vs_solvent", ["CMT167R_CON"], ["CMT167R_PKCi"],
                 "tnk", "Tacstd2", "")
    add_contrast("GSE303943", "PKCi_vs_solvent", ["CMT167R_CON"], ["CMT167R_PKCi"],
                 "tnk", "Cldn4", "")

    add_contrast("GSE133604", "KP_aPD1_vs_IgG", ["KP_Ctrl_IgG"], ["KP_Ctrl_aPD1"],
                 "epithelial", "Tacstd2", "leftover; n=1 vs 1")
    add_contrast("GSE133604", "KP_aPD1_vs_IgG", ["KP_Ctrl_IgG"], ["KP_Ctrl_aPD1"],
                 "epithelial", "Cldn4", "leftover; n=1 vs 1")
    add_contrast("GSE133604", "KP_aPD1_vs_IgG", ["KP_Ctrl_IgG"], ["KP_Ctrl_aPD1"],
                 "tnk", "Tacstd2", "leftover")
    add_contrast("GSE133604", "KP_aPD1_vs_IgG", ["KP_Ctrl_IgG"], ["KP_Ctrl_aPD1"],
                 "tnk", "Cldn4", "leftover")

    add_contrast("GSE129297", "SCLC_aPD1_vs_Ctrl", ["SCLC_Ctrl"], ["SCLC_aPD1"],
                 "epithelial", "Tacstd2", "leftover; n=1 vs 1")
    add_contrast("GSE129297", "SCLC_aPD1_vs_Ctrl", ["SCLC_Ctrl"], ["SCLC_aPD1"],
                 "epithelial", "Cldn4", "leftover; n=1 vs 1")
    add_contrast("GSE129297", "SCLC_aPD1_vs_Ctrl", ["SCLC_Ctrl"], ["SCLC_aPD1"],
                 "tnk", "Tacstd2", "leftover")
    add_contrast("GSE129297", "SCLC_aPD1_vs_Ctrl", ["SCLC_Ctrl"], ["SCLC_aPD1"],
                 "tnk", "Cldn4", "leftover")

    add_contrast("GSE232730", "aPD1_vs_control", ["KPL_CD45_Ctrl"], ["KPL_CD45_aPD1"],
                 "tnk", "Tacstd2", "CD45+ only; leftover; n=1 vs 1")
    add_contrast("GSE232730", "aPD1_vs_control", ["KPL_CD45_Ctrl"], ["KPL_CD45_aPD1"],
                 "tnk", "Cldn4", "CD45+ only")
    add_contrast("GSE222158", "aPD1_vs_control", ["CD45_Ctl"], ["CD45_PD1"],
                 "tnk", "Tacstd2", "CD45+ only; leftover; n=1 vs 1")
    add_contrast("GSE222158", "aPD1_vs_control", ["CD45_Ctl"], ["CD45_PD1"],
                 "tnk", "Cldn4", "CD45+ only")

    # paired epi vs T/NK within samples that have both
    both = []
    for (gse, sample), g in cells.groupby(["gse", "sample"]):
        epi = g[g.compartment == "epithelial"]
        tnk = g[g.compartment == "tnk"]
        if epi.empty or tnk.empty:
            continue
        if int(epi.n_cells.iloc[0]) < 20 or int(tnk.n_cells.iloc[0]) < 20:
            continue
        for gene in ["Tacstd2", "Cldn4"]:
            both.append({
                "gse": gse, "sample": sample, "gene": gene,
                "epi_mean": float(epi[f"{gene}_mean"].iloc[0]),
                "tnk_mean": float(tnk[f"{gene}_mean"].iloc[0]),
                "n_epi": int(epi.n_cells.iloc[0]),
                "n_tnk": int(tnk.n_cells.iloc[0]),
                "delta_epi_minus_tnk": float(epi[f"{gene}_mean"].iloc[0] - tnk[f"{gene}_mean"].iloc[0]),
            })
    both_df = pd.DataFrame(both)
    both_df.to_csv(OUT / "epi_vs_tnk_per_sample.tsv", sep="\t", index=False)

    if not both_df.empty:
        for gene, g in both_df.groupby("gene"):
            d = g["delta_epi_minus_tnk"].astype(float).values
            n_up = int((d > 0).sum())
            try:
                w, p = stats.wilcoxon(g["epi_mean"], g["tnk_mean"], alternative="two-sided")
            except ValueError:
                w, p = np.nan, np.nan
            contrasts.append({
                "gse": "POOL", "contrast": "epithelial_vs_tnk_paired",
                "compartment": "epithelial_minus_tnk", "gene": gene,
                "n_a": int(len(g)), "n_b": int(len(g)),
                "mean_a": float(g.tnk_mean.mean()), "mean_b": float(g.epi_mean.mean()),
                "delta_log": float(d.mean()),
                "direction": f"{n_up}/{len(g)} epi>tnk",
                "welch_p": np.nan, "mwu_p": float(p) if p == p else np.nan,
                "note": f"Wilcoxon signed-rank W={w}; samples with ≥20 cells each compartment",
            })

    cdf = pd.DataFrame(contrasts)
    cdf.to_csv(OUT / "contrasts.tsv", sep="\t", index=False)

    # combined direction table: ICB vs control in epithelial when both exist
    icb_rows = cdf[cdf.contrast.str.contains("aPD1_vs|CA170_vs|SCLC_aPD1")].copy()
    icb_rows.to_csv(OUT / "direction_table.tsv", sep="\t", index=False)

    samp.to_json(OUT / "sample_inventory.json", orient="records", indent=2)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
