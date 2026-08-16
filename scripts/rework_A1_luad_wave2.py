#!/usr/bin/env python3
"""REWORK A1 wave 2 — TCGA-LUAD only.

Question
--------
Prior LUAD-only A1 was weak after ABSOLUTE purity:
  |partial ρ| ≤ 0.22 on marker genes (PR #81)
  |partial ρ| ≤ 0.13 on xCell / MCP / GEP18 / ESTIMATE (PR #169)
OncoSG (PR #139) and TCGA-LUSC (PR #143 / #107) did support a negative
TACSTD2–immune association after purity. The user's pooled claim used
n ≈ 1031 (LUAD+LUSC, not LUAD alone).

This wave asks whether a *different residual* on LUAD can recover a
strong inverse:

  1. Double residual: TACSTD2 residualized on ABSOLUTE purity AND an
     epithelial fraction / KRT+EPCAM score, then correlated with
     CD8 / GEP18 / CYT / xCell CD8.
  2. ESTIMATE ImmuneScore residual of TACSTD2, then the same four
     immune readouts.

Do not pool LUSC. Do not invent numbers. Public Xena + ABSOLUTE only.

Outputs -> results/rework/A1_luad_wave2/
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SIG = DATA / "signatures"
OUT = ROOT / "results" / "rework" / "A1_luad_wave2"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

URLS = {
    "expression": {
        "file": "TCGA.LUAD.HiSeqV2.gz",
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUAD HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "purity": {
        "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
        "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "desc": "PanCanAtlas ABSOLUTE purity (GDC 4f277128-f793-4354-a13d-30cc7fe9f6b5)",
    },
    "estimate": {
        "file": "ESTIMATE_LUAD_RNAseqV2.txt",
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        "desc": "MD Anderson official ESTIMATE RNAseqV2 ImmuneScore / StromalScore / ESTIMATEScore",
    },
    "timer2": {
        "file": "infiltration_estimation_for_tcga.csv.gz",
        "url": "https://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz",
        "desc": "TIMER2.0 / immunedeconv precomputed TCGA infiltration (xCell CD8)",
    },
}

# Locked before looking at results.
# Strong inverse = OncoSG A1 after published purity (PR #139):
#   CD8A partial ρ = −0.309; GEP18 = −0.349; immune 8-gene = −0.318 (n=169).
STRONG_RHO = -0.30
# LUSC CD8 partial ρ = −0.244 (PR #143); prior LUAD marker |ρ| bound = 0.22.
LUSC_BOUND = -0.22
# Prior LUAD deconv primary max |ρ| = 0.126 (PR #169).
DECONV_BOUND = 0.13

PRIMARY_IMMUNE = ["CD8", "CYT", "GEP18_zmean", "xCell_T_cell_CD8"]
SENSITIVITY_IMMUNE = ["CD8B", "GEP18_ssGSEA", "GEP18_logmean", "ESTIMATE_ImmuneScore"]

# Yoshihara Nat Commun 2013 ESTIMATE TumorPurity transform.
ESTIMATE_PURITY_A = 0.6049872018
ESTIMATE_PURITY_B = 0.0001467884

GEP18_GENES = [
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(key: str) -> Path:
    meta = URLS[key]
    dest = DATA / meta["file"]
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"  cached {key}: {dest.name} ({dest.stat().st_size} bytes)")
        return dest
    DATA.mkdir(parents=True, exist_ok=True)
    print(f"downloading {key} -> {dest.name}")
    last_err = None
    for attempt in range(1, 5):
        try:
            r = requests.get(meta["url"], timeout=300)
            r.raise_for_status()
            dest.write_bytes(r.content)
            if dest.stat().st_size < 1000:
                raise RuntimeError(f"{key} download too small: {dest.stat().st_size} bytes")
            return dest
        except Exception as exc:
            last_err = exc
            print(f"  attempt {attempt} failed: {exc}")
    raise SystemExit(f"failed to download {key}: {last_err}")


def sample01(barcode: str) -> str | None:
    parts = str(barcode).replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3]) + "-01"


def bh_fdr(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    q = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    prev = 1.0
    out = np.empty(len(pv))
    for rank_from_end, j in enumerate(order[::-1]):
        rank = len(pv) - rank_from_end
        val = min(prev, pv[j] * len(pv) / rank)
        out[j] = val
        prev = val
    q[idx] = out
    return q


def fisher_z_ci(r: float, n: int, k_covariates: int):
    if r is None or not np.isfinite(r) or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = math.atanh(r)
    se = 1.0 / math.sqrt(n - k_covariates - 3)
    zcrit = stats.norm.ppf(0.975)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def residualize(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def rank_cols(Z: np.ndarray) -> np.ndarray:
    out = np.empty_like(Z, dtype=float)
    for j in range(Z.shape[1]):
        out[:, j] = stats.rankdata(Z[:, j])
    return out


def corr_p_from_r(r: float, n: int, k: int):
    if not np.isfinite(r) or n <= k + 2:
        return np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 2 - k
    if df <= 0:
        return np.nan
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    return float(2 * stats.t.sf(abs(t), df)), df


def spearman_pair(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def onesided_value(x, y, Z):
    """Spearman of value-space residual(x | Z) vs raw y."""
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(Z), axis=1)
    n = int(m.sum())
    k = Z.shape[1]
    if n < k + 6:
        return np.nan, np.nan, n, np.nan
    resid = residualize(x[m], Z[m])
    r, _ = stats.spearmanr(resid, y[m])
    r = float(np.clip(r, -0.999999, 0.999999))
    p, df = corr_p_from_r(r, n, k)
    return r, p, n, df


def onesided_rank(x, y, Z):
    """Pearson of rank residual(x | Z) vs rank(y). User-requested one-sided residual."""
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(Z), axis=1)
    n = int(m.sum())
    k = Z.shape[1]
    if n < k + 6:
        return np.nan, np.nan, n, np.nan
    resid = residualize(stats.rankdata(x[m]), rank_cols(Z[m]))
    r, _ = stats.pearsonr(resid, stats.rankdata(y[m]))
    r = float(np.clip(r, -0.999999, 0.999999))
    p, df = corr_p_from_r(r, n, k)
    return r, p, n, df


def twosided_rank(x, y, Z):
    """Standard multi-covariate partial Spearman (rank residuals on both sides)."""
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(Z), axis=1)
    n = int(m.sum())
    k = Z.shape[1]
    if n < k + 6:
        return np.nan, np.nan, n, np.nan
    Zr = rank_cols(Z[m])
    rx = residualize(stats.rankdata(x[m]), Zr)
    ry = residualize(stats.rankdata(y[m]), Zr)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    p, df = corr_p_from_r(r, n, k)
    return r, p, n, df


def ssgsea_one_set(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 2:
        raise SystemExit(f"ssGSEA gene set too small: {present}")
    n_genes = expr.shape[0]
    ranked = expr.rank(axis=0, method="average", ascending=True) * (10000.0 / n_genes)
    gene_set = set(present)
    scores = {}
    for sample in expr.columns:
        m = ranked[sample]
        order = m.sort_values(ascending=False).index
        m_ord = m.loc[order].to_numpy(float)
        hits = np.fromiter((g in gene_set for g in order), dtype=bool, count=len(order))
        w = np.abs(m_ord) ** tau
        w_hit = np.where(hits, w, 0.0)
        nhit = float(w_hit.sum())
        nmiss = float((~hits).sum())
        if nhit <= 0 or nmiss <= 0:
            scores[sample] = np.nan
            continue
        p_hit = np.cumsum(w_hit) / nhit
        p_miss = np.cumsum((~hits).astype(float)) / nmiss
        scores[sample] = float(np.sum(p_hit - p_miss))
    return pd.Series(scores, name="GEP18_ssGSEA")


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def call_rho(r: float) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    if r <= STRONG_RHO:
        return "STRONG_INVERSE"
    if r <= LUSC_BOUND:
        return "REACHES_LUSC_BOUND"
    if r < 0:
        return "WEAK_NEGATIVE"
    return "NOT_INVERSE"


def main():
    print("=== A1 LUAD wave 2: double residual ===")
    for k in URLS:
        download(k)

    print("Loading Xena LUAD expression ...")
    expr = pd.read_csv(DATA / URLS["expression"]["file"], sep="\t", index_col=0)
    expr = expr.loc[:, [c for c in expr.columns if str(c).endswith("-01")]]
    if expr.columns.duplicated().any():
        expr = expr.T.groupby(level=0).mean().T
    if "TACSTD2" not in expr.index:
        raise SystemExit("TACSTD2 missing from Xena HiSeqV2")
    tacstd2 = expr.loc["TACSTD2"].astype(float)
    n_expr = int(tacstd2.shape[0])
    print(f"  primary tumors with RNA: {n_expr}")

    epi_tab = pd.read_csv(SIG / "epithelial_krt_epcam.tsv", sep="\t")
    primary_epi = epi_tab.loc[epi_tab["panel"] == "primary", "hugo"].tolist()
    sens_epi = epi_tab.loc[epi_tab["panel"].isin(["primary", "sensitivity"]), "hugo"].tolist()
    epi_present = [g for g in primary_epi if g in expr.index]
    epi_missing = [g for g in primary_epi if g not in expr.index]
    if len(epi_present) < 2:
        raise SystemExit(f"KRT/EPCAM primary genes missing: {epi_missing}")
    krt_epcam = expr.loc[epi_present].mean(axis=0)
    krt_epcam.name = "KRT_EPCAM"
    sens_present = [g for g in sens_epi if g in expr.index]
    krt_epcam_cdh1 = expr.loc[sens_present].mean(axis=0)
    krt_epcam_cdh1.name = "KRT_EPCAM_CDH1"
    krt_all = [g for g in expr.index if re.fullmatch(r"KRT[0-9]+[A-Z]?", str(g))]
    all_krt_epcam_genes = ["EPCAM"] + krt_all if "EPCAM" in expr.index else krt_all
    all_krt_epcam = expr.loc[all_krt_epcam_genes].mean(axis=0)
    all_krt_epcam.name = "all_KRT_EPCAM"

    needed = ["TACSTD2", "CD8A", "CD8B", "GZMA", "PRF1"] + GEP18_GENES + primary_epi
    missing_core = [g for g in needed if g not in expr.index]
    if missing_core:
        raise SystemExit(f"core genes missing from Xena: {missing_core}")

    print("Loading ABSOLUTE ...")
    ab = pd.read_csv(DATA / URLS["purity"]["file"], sep="\t")
    ab = ab[ab["array"].astype(str).str.endswith("-01", na=False)][["array", "purity"]].dropna()
    purity = ab.groupby("array")["purity"].mean()

    print("Loading official ESTIMATE ...")
    est = pd.read_csv(DATA / URLS["estimate"]["file"], sep="\t")
    est.columns = [c.strip() for c in est.columns]
    idcol = est.columns[0]
    est["sample"] = est[idcol].map(lambda b: sample01(b) or (str(b) if str(b).endswith("-01") else None))
    est = est.dropna(subset=["sample"])
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if cl in ("immune_score", "immunescore"):
            rename[c] = "ESTIMATE_ImmuneScore"
        elif cl in ("stromal_score", "stromalscore"):
            rename[c] = "ESTIMATE_StromalScore"
        elif cl in ("estimate_score", "estimatescore"):
            rename[c] = "ESTIMATE_Score"
    est = est.rename(columns=rename)
    est_keep = [c for c in ("ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_Score") if c in est.columns]
    if "ESTIMATE_ImmuneScore" not in est_keep or "ESTIMATE_Score" not in est_keep:
        raise SystemExit(f"ESTIMATE columns not found: {list(est.columns)}")
    est = est.groupby("sample")[est_keep].mean()
    est["ESTIMATE_TumorPurity"] = np.cos(ESTIMATE_PURITY_A + ESTIMATE_PURITY_B * est["ESTIMATE_Score"])

    print("Loading TIMER2 xCell ...")
    timer = pd.read_csv(DATA / URLS["timer2"]["file"])
    timer["sample"] = timer["cell_type"].map(lambda b: sample01(b) or (str(b) if str(b).endswith("-01") else None))
    timer = timer.dropna(subset=["sample"])
    xcell_cols = [c for c in timer.columns if c.endswith("_XCELL")]
    timer_x = timer.groupby("sample")[xcell_cols].mean()

    def find_xcell(*needles: str) -> str | None:
        hits = []
        for c in timer_x.columns:
            cl = c.lower().replace("+", "").replace("-", "_").replace(" ", "_")
            if all(n.lower().replace(" ", "_") in cl for n in needles):
                hits.append(c)
        if not hits:
            return None
        # prefer the shortest / least-qualified name
        hits.sort(key=lambda s: (len(s), s))
        return hits[0]

    xcell_cd8_col = find_xcell("t_cell", "cd8") or find_xcell("cd8")
    # exact-ish CD8, not naive/memory if a plain CD8 exists
    plain = [c for c in timer_x.columns if re.search(r"T cell CD8\+_XCELL$", c) or c == "T cell CD8+_XCELL"]
    if plain:
        xcell_cd8_col = plain[0]
    xcell_epi_col = find_xcell("epithelial")
    xcell_ker_col = find_xcell("keratinocyte")
    print(f"  xCell CD8 column: {xcell_cd8_col}")
    print(f"  xCell epithelial column: {xcell_epi_col}")
    print(f"  xCell keratinocyte column: {xcell_ker_col}")
    if xcell_cd8_col is None:
        raise SystemExit(f"xCell CD8 column not found. Available: {list(timer_x.columns)}")

    print("Scoring GEP18 ...")
    gep_ssgsea = ssgsea_one_set(expr, GEP18_GENES, tau=0.25)
    z = expr.loc[GEP18_GENES].T
    gep_zmean = ((z - z.mean(axis=0)) / z.std(axis=0, ddof=0)).mean(axis=1)
    gep_zmean.name = "GEP18_zmean"
    gep_logmean = expr.loc[GEP18_GENES].mean(axis=0)
    gep_logmean.name = "GEP18_logmean"

    print("Merging analysis table ...")
    df = pd.DataFrame({"TACSTD2": tacstd2})
    df["CD8"] = expr.loc["CD8A"].astype(float)
    df["CD8B"] = expr.loc["CD8B"].astype(float)
    df["CYT"] = expr.loc[["GZMA", "PRF1"]].mean(axis=0)
    df = df.join(gep_zmean, how="left")
    df = df.join(gep_ssgsea, how="left")
    df = df.join(gep_logmean, how="left")
    df = df.join(purity.rename("ABSOLUTE_purity"), how="left")
    df = df.join(est, how="left")
    df = df.join(krt_epcam, how="left")
    df = df.join(krt_epcam_cdh1, how="left")
    df = df.join(all_krt_epcam, how="left")
    df["xCell_T_cell_CD8"] = timer_x[xcell_cd8_col].reindex(df.index)
    if xcell_epi_col:
        df["xCell_Epithelial"] = timer_x[xcell_epi_col].reindex(df.index)
    if xcell_ker_col:
        df["xCell_Keratinocytes"] = timer_x[xcell_ker_col].reindex(df.index)
    df.index.name = "sample"

    # Covariate models. Primary double residual is ABSOLUTE + KRT/EPCAM.
    models = {
        "unadjusted": {"cols": [], "k": 0, "role": "baseline"},
        "ABS_only": {"cols": ["ABSOLUTE_purity"], "k": 1, "role": "prior_A1_replication"},
        "KRT_only": {"cols": ["KRT_EPCAM"], "k": 1, "role": "epithelial_only"},
        "ABS_KRT": {"cols": ["ABSOLUTE_purity", "KRT_EPCAM"], "k": 2, "role": "primary_double_residual"},
        "ABS_KRT_CDH1": {"cols": ["ABSOLUTE_purity", "KRT_EPCAM_CDH1"], "k": 2, "role": "sensitivity_epithelial"},
        "ABS_allKRT": {"cols": ["ABSOLUTE_purity", "all_KRT_EPCAM"], "k": 2, "role": "sensitivity_epithelial"},
        "ABS_ESTPUR": {"cols": ["ABSOLUTE_purity", "ESTIMATE_TumorPurity"], "k": 2, "role": "sensitivity_rna_purity"},
        "IMMUNE": {"cols": ["ESTIMATE_ImmuneScore"], "k": 1, "role": "estimate_immunescore_residual"},
        "ABS_IMMUNE": {"cols": ["ABSOLUTE_purity", "ESTIMATE_ImmuneScore"], "k": 2, "role": "sensitivity_abs_plus_immune"},
    }
    if "xCell_Epithelial" in df.columns:
        models["ABS_xCellEpi"] = {
            "cols": ["ABSOLUTE_purity", "xCell_Epithelial"],
            "k": 2,
            "role": "sensitivity_xcell_epithelial",
        }

    immune_feats = [f for f in PRIMARY_IMMUNE + SENSITIVITY_IMMUNE if f in df.columns]

    rows = []
    for model_name, spec in models.items():
        cols = spec["cols"]
        k = spec["k"]
        for feat in immune_feats:
            x = df["TACSTD2"].to_numpy(float)
            y = df[feat].to_numpy(float)
            if k == 0:
                r, p, n = spearman_pair(x, y)
                rows.append(
                    {
                        "model": model_name,
                        "role": spec["role"],
                        "feature": feat,
                        "estimator": "unadjusted_spearman",
                        "rho": r,
                        "p": p,
                        "n": n,
                        "df": n - 2 if n >= 3 else np.nan,
                        "k_covariates": 0,
                        "ci_low": fisher_z_ci(r, n, 0)[0],
                        "ci_high": fisher_z_ci(r, n, 0)[1],
                        "call": call_rho(r),
                    }
                )
                continue
            Z = df[cols].to_numpy(float)
            for est_name, fn in (
                ("onesided_rank", onesided_rank),
                ("onesided_value", onesided_value),
                ("twosided_rank", twosided_rank),
            ):
                r, p, n, dfree = fn(x, y, Z)
                lo, hi = fisher_z_ci(r, n, k)
                rows.append(
                    {
                        "model": model_name,
                        "role": spec["role"],
                        "feature": feat,
                        "estimator": est_name,
                        "rho": r,
                        "p": p,
                        "n": n,
                        "df": dfree,
                        "k_covariates": k,
                        "ci_low": lo,
                        "ci_high": hi,
                        "call": call_rho(r),
                    }
                )
    res = pd.DataFrame(rows)
    # FDR within each model × estimator across the 4 primary immune features
    res["fdr_primary4"] = np.nan
    for (model, est_name), idx in res.groupby(["model", "estimator"]).groups.items():
        mask = res.index.isin(idx) & res["feature"].isin(PRIMARY_IMMUNE)
        if mask.sum():
            res.loc[mask, "fdr_primary4"] = bh_fdr(res.loc[mask, "p"])

    # Context correlations (not residual models)
    context_pairs = [
        ("TACSTD2", "ABSOLUTE_purity"),
        ("TACSTD2", "KRT_EPCAM"),
        ("TACSTD2", "KRT_EPCAM_CDH1"),
        ("TACSTD2", "all_KRT_EPCAM"),
        ("TACSTD2", "ESTIMATE_TumorPurity"),
        ("TACSTD2", "ESTIMATE_ImmuneScore"),
        ("ABSOLUTE_purity", "KRT_EPCAM"),
        ("ABSOLUTE_purity", "ESTIMATE_TumorPurity"),
        ("ABSOLUTE_purity", "ESTIMATE_ImmuneScore"),
        ("KRT_EPCAM", "ESTIMATE_TumorPurity"),
        ("KRT_EPCAM", "ESTIMATE_ImmuneScore"),
        ("ESTIMATE_TumorPurity", "ESTIMATE_ImmuneScore"),
        ("CD8", "KRT_EPCAM"),
        ("CYT", "KRT_EPCAM"),
        ("GEP18_zmean", "KRT_EPCAM"),
        ("xCell_T_cell_CD8", "KRT_EPCAM"),
        ("CD8", "ABSOLUTE_purity"),
        ("CYT", "ABSOLUTE_purity"),
        ("GEP18_zmean", "ABSOLUTE_purity"),
        ("xCell_T_cell_CD8", "ABSOLUTE_purity"),
        ("CD8", "ESTIMATE_ImmuneScore"),
        ("CYT", "ESTIMATE_ImmuneScore"),
        ("GEP18_zmean", "ESTIMATE_ImmuneScore"),
        ("xCell_T_cell_CD8", "ESTIMATE_ImmuneScore"),
    ]
    if "xCell_Epithelial" in df.columns:
        context_pairs += [
            ("TACSTD2", "xCell_Epithelial"),
            ("ABSOLUTE_purity", "xCell_Epithelial"),
            ("KRT_EPCAM", "xCell_Epithelial"),
        ]
    ctx_rows = []
    for a, b in context_pairs:
        if a not in df.columns or b not in df.columns:
            continue
        r, p, n = spearman_pair(df[a].to_numpy(float), df[b].to_numpy(float))
        ctx_rows.append({"a": a, "b": b, "spearman_rho": r, "p": p, "n": n})
    ctx = pd.DataFrame(ctx_rows)

    def ctx_rho(a, b):
        hit = ctx[(ctx.a == a) & (ctx.b == b)]
        if hit.empty:
            return np.nan, np.nan, 0
        r = hit.iloc[0]
        return float(r.spearman_rho), float(r.p), int(r.n)

    n_abs = int(df["ABSOLUTE_purity"].notna().sum())
    n_abs_krt = int((df["ABSOLUTE_purity"].notna() & df["KRT_EPCAM"].notna()).sum())
    n_imm = int(df["ESTIMATE_ImmuneScore"].notna().sum())
    n_xcell = int(df["xCell_T_cell_CD8"].notna().sum())

    # ---------- write tables ----------
    res.to_csv(OUT / "correlations.tsv", sep="\t", index=False)
    ctx.to_csv(OUT / "context_correlations.tsv", sep="\t", index=False)
    df.reset_index().to_csv(OUT / "sample_table.tsv", sep="\t", index=False)

    prim = res[
        res["feature"].isin(PRIMARY_IMMUNE)
        & res["model"].isin(["unadjusted", "ABS_only", "ABS_KRT", "IMMUNE"])
        & res["estimator"].isin(["unadjusted_spearman", "onesided_rank", "twosided_rank"])
    ].copy()
    prim.to_csv(OUT / "primary_results.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        [
            {
                "set": "KRT_EPCAM_primary",
                "n_official": len(primary_epi),
                "n_present": len(epi_present),
                "present": ",".join(epi_present),
                "missing": ",".join(epi_missing),
            },
            {
                "set": "KRT_EPCAM_CDH1",
                "n_official": len(sens_epi),
                "n_present": len(sens_present),
                "present": ",".join(sens_present),
                "missing": ",".join([g for g in sens_epi if g not in expr.index]),
            },
            {
                "set": "all_KRT_EPCAM",
                "n_official": len(all_krt_epcam_genes),
                "n_present": len(all_krt_epcam_genes),
                "present": ",".join(all_krt_epcam_genes),
                "missing": "",
            },
            {
                "set": "GEP18",
                "n_official": 18,
                "n_present": 18,
                "present": ",".join(GEP18_GENES),
                "missing": "",
            },
            {
                "set": "CYT",
                "n_official": 2,
                "n_present": 2,
                "present": "GZMA,PRF1",
                "missing": "",
            },
        ]
    )
    coverage.to_csv(OUT / "gene_coverage.tsv", sep="\t", index=False)

    # ---------- figures ----------
    order_feat = [f for f in PRIMARY_IMMUNE if f in df.columns]
    compare_models = [
        ("ABS_only", "onesided_rank", "ABS only (one-sided)"),
        ("ABS_KRT", "onesided_rank", "ABS+KRT/EPCAM (one-sided residual)"),
        ("ABS_KRT", "twosided_rank", "ABS+KRT/EPCAM (two-sided partial)"),
        ("IMMUNE", "onesided_rank", "ImmuneScore residual (one-sided)"),
        ("IMMUNE", "twosided_rank", "ImmuneScore residual (two-sided)"),
    ]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(len(order_feat))
    width = 0.16
    colors = ["#7f8c8d", "#c0392b", "#922b21", "#2980b9", "#1a5276"]
    for i, (model, est_name, label) in enumerate(compare_models):
        sub = res[(res.model == model) & (res.estimator == est_name)].set_index("feature")
        ys = [sub.loc[f, "rho"] if f in sub.index else np.nan for f in order_feat]
        los = [sub.loc[f, "ci_low"] if f in sub.index else np.nan for f in order_feat]
        his = [sub.loc[f, "ci_high"] if f in sub.index else np.nan for f in order_feat]
        xpos = x + (i - 2) * width
        ax.bar(xpos, ys, width=width, color=colors[i], label=label, edgecolor="none")
        yerr = np.vstack(
            [
                np.array(ys) - np.array(los),
                np.array(his) - np.array(ys),
            ]
        )
        ax.errorbar(xpos, ys, yerr=yerr, fmt="none", ecolor="black", elinewidth=0.7, capsize=2)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(STRONG_RHO, color="#8e44ad", ls="--", lw=1.0, label=f"OncoSG strong inverse ({STRONG_RHO:.2f})")
    ax.axhline(LUSC_BOUND, color="#d35400", ls=":", lw=1.0, label=f"LUSC / prior LUAD marker bound ({LUSC_BOUND:.2f})")
    ax.set_xticks(x)
    ax.set_xticklabels(["CD8A", "CYT (GZMA/PRF1)", "GEP18 z-mean", "xCell CD8"], fontsize=9)
    ax.set_ylabel("Spearman ρ (TACSTD2 vs immune)")
    ax.set_title("TCGA-LUAD only: does a second residual recover a strong TACSTD2–immune inverse?")
    ax.set_ylim(-0.45, 0.25)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "bar_residual_models.png", dpi=160)
    plt.close(fig)

    # Forest of primary double residual (both estimators)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for ax, est_name, title in (
        (axes[0], "onesided_rank", "One-sided: residual(TACSTD2 | ABS, KRT/EPCAM) vs raw immune"),
        (axes[1], "twosided_rank", "Two-sided partial Spearman | ABS + KRT/EPCAM"),
    ):
        sub = res[(res.model == "ABS_KRT") & (res.estimator == est_name) & (res.feature.isin(order_feat))]
        sub = sub.set_index("feature").loc[order_feat]
        y = np.arange(len(sub))
        cols = ["#c0392b" if r < 0 else "#2471a3" for r in sub["rho"]]
        ax.barh(y, sub["rho"], color=cols, edgecolor="none")
        ax.errorbar(
            sub["rho"],
            y,
            xerr=[sub["rho"] - sub["ci_low"], sub["ci_high"] - sub["rho"]],
            fmt="none",
            ecolor="black",
            elinewidth=0.8,
            capsize=2,
        )
        ax.axvline(0, color="k", lw=0.8)
        ax.axvline(STRONG_RHO, color="#8e44ad", ls="--", lw=0.9)
        ax.axvline(LUSC_BOUND, color="#d35400", ls=":", lw=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(order_feat, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("ρ")
        ax.set_title(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "forest_double_residual.png", dpi=160)
    plt.close(fig)

    # Residual vs immune scatters for primary double residual (value-space residual)
    m = df["ABSOLUTE_purity"].notna() & df["KRT_EPCAM"].notna()
    Z = df.loc[m, ["ABSOLUTE_purity", "KRT_EPCAM"]].to_numpy(float)
    tac_res = pd.Series(residualize(df.loc[m, "TACSTD2"].to_numpy(float), Z), index=df.index[m])
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.2), constrained_layout=True)
    lookup = res[(res.model == "ABS_KRT") & (res.estimator == "onesided_value")].set_index("feature")
    for ax, feat, ylab in zip(
        axes.ravel(),
        order_feat,
        ["CD8A log2(RSEM+1)", "CYT (mean GZMA, PRF1)", "GEP18 z-mean", "xCell CD8"],
    ):
        mm = tac_res.index.intersection(df.index[df[feat].notna()])
        ax.scatter(tac_res.loc[mm], df.loc[mm, feat], s=10, alpha=0.7, c="#1f618d", linewidths=0)
        r = lookup.loc[feat]
        ax.set_title(
            f"{feat}\nonesided-value ρ={r.rho:.3f} p={r.p:.2g} n={int(r.n)}",
            fontsize=9,
        )
        ax.set_xlabel("TACSTD2 residual | ABSOLUTE + KRT/EPCAM", fontsize=8)
        ax.set_ylabel(ylab, fontsize=8)
    fig.suptitle("TCGA-LUAD: value-space TACSTD2 double residual vs immune (raw Y)", fontsize=11)
    fig.savefig(FIG / "scatter_double_residual.png", dpi=150)
    plt.close(fig)

    # Context: TACSTD2 vs epithelial / purity / ImmuneScore
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    panels = [
        ("ABSOLUTE_purity", "ABSOLUTE purity"),
        ("KRT_EPCAM", "KRT/EPCAM score"),
        ("ESTIMATE_ImmuneScore", "ESTIMATE ImmuneScore"),
    ]
    for ax, (col, xlab) in zip(axes, panels):
        mm = df[col].notna()
        ax.scatter(df.loc[mm, col], df.loc[mm, "TACSTD2"], s=9, alpha=0.7, c="#6c3483", linewidths=0)
        r, p, n = ctx_rho("TACSTD2", col)
        ax.set_title(f"TACSTD2 vs {xlab}\nρ={r:.3f} p={fmt_p(p)} n={n}", fontsize=9)
        ax.set_xlabel(xlab, fontsize=8)
        ax.set_ylabel("TACSTD2 log2(RSEM+1)", fontsize=8)
    fig.savefig(FIG / "scatter_tacstd2_vs_covariates.png", dpi=150)
    plt.close(fig)

    # ---------- verdict ----------
    def pick(model, feat, estimator):
        hit = res[(res.model == model) & (res.feature == feat) & (res.estimator == estimator)]
        if hit.empty:
            return None
        return hit.iloc[0]

    def max_abs(model, estimator, feats=None):
        feats = feats or PRIMARY_IMMUNE
        sub = res[(res.model == model) & (res.estimator == estimator) & (res.feature.isin(feats))]
        if sub.empty:
            return np.nan
        return float(sub["rho"].abs().max())

    def min_rho(model, estimator, feats=None):
        feats = feats or PRIMARY_IMMUNE
        sub = res[(res.model == model) & (res.estimator == estimator) & (res.feature.isin(feats))]
        if sub.empty:
            return np.nan
        return float(sub["rho"].min())

    onesided_min = min_rho("ABS_KRT", "onesided_rank")
    twosided_min = min_rho("ABS_KRT", "twosided_rank")
    immune_one_min = min_rho("IMMUNE", "onesided_rank")
    immune_two_min = min_rho("IMMUNE", "twosided_rank")
    abs_one_min = min_rho("ABS_only", "onesided_rank")

    def verdict_from(min_r, label):
        if min_r <= STRONG_RHO:
            return (
                "REACHES_STRONG_INVERSE",
                f"{label} reaches the OncoSG-like strong inverse (most negative ρ = {min_r:.3f} ≤ {STRONG_RHO:.2f}).",
            )
        if min_r <= LUSC_BOUND:
            return (
                "REACHES_LUSC_NOT_ONCOSG",
                f"{label} reaches the LUSC / prior-LUAD-marker bound (most negative ρ = {min_r:.3f} ≤ {LUSC_BOUND:.2f}) "
                f"but not the OncoSG strong inverse ({STRONG_RHO:.2f}).",
            )
        return (
            "STILL_WEAK",
            f"{label} does not reach a strong inverse (most negative ρ = {min_r:.3f}; "
            f"|ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30).",
        )

    tag_one, txt_one = verdict_from(onesided_min, "User-requested one-sided double residual (TACSTD2 | ABS+KRT/EPCAM)")
    tag_two, txt_two = verdict_from(twosided_min, "Two-sided partial Spearman | ABS+KRT/EPCAM")
    tag_imm, txt_imm = verdict_from(immune_one_min, "ESTIMATE ImmuneScore residual of TACSTD2 (one-sided)")
    tag_imm2, txt_imm2 = verdict_from(immune_two_min, "Two-sided partial | ImmuneScore")

    # Overall: the claim is rescued only if the *fair* two-sided partial reaches strong,
    # or if the user-requested one-sided does AND two-sided is not a collapse.
    if tag_two == "REACHES_STRONG_INVERSE":
        overall_tag = "REACHES_STRONG_INVERSE"
    elif tag_one == "REACHES_STRONG_INVERSE" and tag_two != "STILL_WEAK":
        overall_tag = "ONESIDED_ONLY_STRONG"
    elif tag_two == "REACHES_LUSC_NOT_ONCOSG" or tag_one == "REACHES_LUSC_NOT_ONCOSG":
        overall_tag = "REACHES_LUSC_NOT_ONCOSG"
    else:
        overall_tag = "STILL_WEAK"

    rho_tp, p_tp, n_tp = ctx_rho("TACSTD2", "ABSOLUTE_purity")
    rho_tk, p_tk, n_tk = ctx_rho("TACSTD2", "KRT_EPCAM")
    rho_ti, p_ti, n_ti = ctx_rho("TACSTD2", "ESTIMATE_ImmuneScore")
    rho_ak, p_ak, n_ak = ctx_rho("ABSOLUTE_purity", "KRT_EPCAM")
    rho_ai, p_ai, n_ai = ctx_rho("ABSOLUTE_purity", "ESTIMATE_ImmuneScore")

    def row_md(model, feat, estimator):
        r = pick(model, feat, estimator)
        if r is None:
            return f"| `{feat}` | {model} / {estimator} | NA | NA | NA | NA | NA |"
        return (
            f"| `{feat}` | {fmt_rho(r.rho)} | {fmt_p(r.p)} | {fmt_p(r.fdr_primary4)} | "
            f"{int(r.n)} | {fmt_rho(r.ci_low)} to {fmt_rho(r.ci_high)} | {r.call} |"
        )

    def block(model, estimator):
        lines = [
            "| Feature | ρ | p | FDR (4) | n | 95% CI | Call |",
            "|---|---:|---:|---:|---:|---|---|",
        ]
        for feat in order_feat:
            lines.append(row_md(model, feat, estimator))
        return "\n".join(lines)

    def all_model_table(estimator):
        lines = [
            "| Model | CD8 | CYT | GEP18 z-mean | xCell CD8 | most negative |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for model in models:
            if model == "unadjusted" and estimator != "unadjusted_spearman":
                est_use = "unadjusted_spearman"
            else:
                est_use = estimator if model != "unadjusted" else "unadjusted_spearman"
            vals = []
            for feat in order_feat:
                r = pick(model, feat, est_use)
                vals.append(float(r.rho) if r is not None and np.isfinite(r.rho) else np.nan)
            arr = np.asarray(vals, dtype=float)
            mn = float(np.nanmin(arr)) if np.isfinite(arr).any() else np.nan
            lines.append(
                "| `"
                + model
                + "` | "
                + " | ".join(fmt_rho(v) for v in vals)
                + f" | {fmt_rho(mn)} |"
            )
        return "\n".join(lines)

    extra_note = ""
    if "xCell_Epithelial" not in df.columns:
        extra_note = (
            "TIMER2 xCell has **no Epithelial-cells column** in this freeze "
            f"(keratinocyte column: {xcell_ker_col or 'also absent'}). "
            "Epithelial fraction is therefore the KRT/EPCAM gene score and ESTIMATE TumorPurity, "
            "not an xCell epithelial fraction. That absence is reported, not filled in."
        )
    else:
        extra_note = f"TIMER2 xCell epithelial column used as sensitivity: `{xcell_epi_col}`."

    # Honest circularity flags
    circ_immune = (
        "Residualizing TACSTD2 on ESTIMATE ImmuneScore and then correlating with CD8/GEP18/CYT/xCell CD8 "
        "is valid (different features). Correlating that residual with ImmuneScore itself would be "
        "circular and is not used as a primary readout. ImmuneScore is strongly anti-correlated with "
        f"ABSOLUTE (ρ = {rho_ai:.3f}) by construction of leukocyte content."
    )

    report = f"""# REWORK A1 wave 2 — TCGA-LUAD double residual (ABSOLUTE + epithelium) and ImmuneScore residual

**Self-contained. Public data only. No fabricated numbers.**

**Cohort:** TCGA-LUAD primary tumors (`-01`) only. Not LUSC. Not OncoSG. Not the pooled n ≈ 1031.

**Question:** Prior LUAD-only A1 was weak after ABSOLUTE (|ρ| ≤ 0.22 markers, PR #81; |ρ| ≤ 0.13 xCell/MCP/GEP18/ESTIMATE, PR #169). OncoSG and LUSC *did* support a negative TACSTD2–immune association after purity. Can a second residual — epithelium or ESTIMATE ImmuneScore — recover the user's implied **strong inverse** on LUAD?

## Verdict: {overall_tag}

**User-requested one-sided double residual** (TACSTD2 residualized on ABSOLUTE **and** KRT/EPCAM, then Spearman vs raw immune): {txt_one}

**Two-sided partial** (both TACSTD2 and the immune feature residualized on ABSOLUTE + KRT/EPCAM): {txt_two}

**ESTIMATE ImmuneScore residual of TACSTD2** (one-sided): {txt_imm}

**Two-sided partial | ImmuneScore:** {txt_imm2}

Prior ABSOLUTE-only one-sided residual (replication of PR #81 / #169 style): most negative primary ρ = {abs_one_min:.3f}.

**Implied strong inverse, locked before results.** OncoSG A1 after published purity (PR #139, n=169) is the user's working example of a supported inverse: CD8A partial ρ = −0.309, GEP18 = −0.349. Thresholds:

| Call | Rule |
|---|---|
| `STRONG_INVERSE` | ρ ≤ **−0.30** (OncoSG CD8/GEP) |
| `REACHES_LUSC_BOUND` | −0.30 < ρ ≤ **−0.22** (LUSC CD8 / prior LUAD marker bound) |
| `WEAK_NEGATIVE` | −0.22 < ρ < 0 |
| `NOT_INVERSE` | ρ ≥ 0 |

n ≈ 1031 in the original claim is the **pooled TCGA LUAD+LUSC** ESTIMATE-complete set (PR #89), not LUAD-only. This file stays LUAD-only. Complete-case n is in every table cell (~{n_abs} with ABSOLUTE).

## Why one-sided and two-sided are both shown

The request is: residualize **TACSTD2** on (1) ABSOLUTE and (2) epithelial fraction / KRT/EPCAM, **then** correlate with CD8 / GEP18 / CYT / xCell CD8.

That is a **one-sided** residual. Immune scores still carry epithelial / purity composition. Because KRT/EPCAM is expected to anti-correlate with CD8 in bulk RNA, one-sided residualization can *inflate* a negative ρ. The **two-sided** rank residual (standard multi-covariate partial Spearman) is the fair test of “TACSTD2 vs immune beyond tumor content and epithelium.” Both are reported. The overall tag uses the two-sided result as the claim-level answer and states the one-sided number separately.

{circ_immune}

{extra_note}

## Analysis set

| Filter | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`) | {n_expr} |
| + ABSOLUTE purity | {n_abs} |
| + KRT/EPCAM score (same matrix) | {n_abs_krt} |
| + official ESTIMATE ImmuneScore | {n_imm} |
| + TIMER2 xCell CD8 | {n_xcell} |

Primary tumors only. One row per 15-character barcode. Replicate aliquots averaged.

## Definitions (locked)

### TACSTD2
UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` **log2(RSEM normalized_count + 1)**. Gene symbol `TACSTD2`. Not protein. Not ADC. Not ICI.

### Samples
TCGA-LUAD only. Sample type `01`. No LUSC, no OncoSG, no normals.

### Covariates
1. **ABSOLUTE purity** — PanCanAtlas GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5`, column `purity`.
2. **KRT/EPCAM score (primary epithelium)** — mean log2(RSEM+1) of `{', '.join(primary_epi)}`. Present: {', '.join(epi_present)}. Missing (not imputed): {', '.join(epi_missing) if epi_missing else 'none'}. **TACSTD2 is not in the score.**
3. **ESTIMATE TumorPurity (sensitivity epithelium)** — Yoshihara 2013 transform `cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` from the official MD Anderson LUAD RNAseqV2 table. This is RNA tumor fraction, not ABSOLUTE.
4. **ESTIMATE ImmuneScore** — official `Immune_score` from the same table. Used only as the requested TACSTD2 residual covariate, not as a primary immune *outcome*.
5. **xCell epithelial** — used only if a TIMER2 `*_XCELL` epithelial column exists; otherwise skipped (see above).

Sensitivity epithelium: `{', '.join(sens_epi)}` mean; and mean of EPCAM + every `KRT[0-9]+` gene in HiSeqV2 ({len(all_krt_epcam_genes)} genes).

### Immune outcomes (primary)
| Name | Definition |
|---|---|
| CD8 | `CD8A` log2(RSEM+1) on the same Xena matrix |
| CYT | Rooney *Cell* 2015 cytolytic score = mean(`GZMA`, `PRF1`) on log2 data |
| GEP18 | Ayers *JCI* 2017 18 genes, **unweighted within-LUAD z-mean** (histology-rework definition, PR #107) |
| xCell CD8 | TIMER2.0 immunedeconv column `{xcell_cd8_col}` |

Sensitivity outcomes (not in the verdict): `CD8B`, GEP18 ssGSEA (Barbie/GSVA tau=0.25, ranks 1..10000), GEP18 log-mean, official ESTIMATE ImmuneScore as a *Y* (only for ABS/KRT models).

### Estimators
- **onesided_rank (user-requested primary):** OLS of rank(TACSTD2) on rank(covariates); Pearson of that residual vs rank(immune). df = n − 2 − k.
- **onesided_value:** OLS of TACSTD2 on raw covariates; Spearman of residual vs raw immune. Same df.
- **twosided_rank:** both sides residualized on ranked covariates (partial Spearman). Same df.
- Unadjusted Spearman is the no-covariate baseline.
- 95% Fisher-z CI uses variance `1/(n − k − 3)`.
- BH-FDR across the 4 primary immune features, within each model × estimator.

## Covariate context (unadjusted Spearman)

| Pair | ρ | p | n |
|---|---:|---:|---:|
| TACSTD2 vs ABSOLUTE | {fmt_rho(rho_tp)} | {fmt_p(p_tp)} | {n_tp} |
| TACSTD2 vs KRT/EPCAM | {fmt_rho(rho_tk)} | {fmt_p(p_tk)} | {n_tk} |
| TACSTD2 vs ESTIMATE ImmuneScore | {fmt_rho(rho_ti)} | {fmt_p(p_ti)} | {n_ti} |
| ABSOLUTE vs KRT/EPCAM | {fmt_rho(rho_ak)} | {fmt_p(p_ak)} | {n_ak} |
| ABSOLUTE vs ImmuneScore | {fmt_rho(rho_ai)} | {fmt_p(p_ai)} | {n_ai} |

If TACSTD2 is nearly orthogonal to ABSOLUTE (as in PR #81 / #169), purity-only adjustment cannot create a strong inverse. A second residual can only help if TACSTD2 shares variance with epithelium (or ImmuneScore) that was *masking* an inverse with CD8/GEP/CYT.

## Primary results

### Unadjusted

{block("unadjusted", "unadjusted_spearman")}

### ABSOLUTE only (prior A1 replication, one-sided rank residual)

{block("ABS_only", "onesided_rank")}

### Double residual: ABSOLUTE + KRT/EPCAM — one-sided (requested)

{block("ABS_KRT", "onesided_rank")}

### Double residual: ABSOLUTE + KRT/EPCAM — two-sided partial (fair test)

{block("ABS_KRT", "twosided_rank")}

### ESTIMATE ImmuneScore residual of TACSTD2 — one-sided

{block("IMMUNE", "onesided_rank")}

### ESTIMATE ImmuneScore residual — two-sided partial

{block("IMMUNE", "twosided_rank")}

## All models, one-sided rank residual (primary four)

{all_model_table("onesided_rank")}

## All models, two-sided rank partial (primary four)

{all_model_table("twosided_rank")}

## Does |ρ| reach the implied strong inverse?

| Test | Most negative primary ρ | Strong (≤ −0.30)? | LUSC bound (≤ −0.22)? |
|---|---:|---|---|
| Unadjusted | {min_rho("unadjusted", "unadjusted_spearman"):.3f} | {str(min_rho("unadjusted", "unadjusted_spearman") <= STRONG_RHO)} | {str(min_rho("unadjusted", "unadjusted_spearman") <= LUSC_BOUND)} |
| ABS only, one-sided | {abs_one_min:.3f} | {str(abs_one_min <= STRONG_RHO)} | {str(abs_one_min <= LUSC_BOUND)} |
| ABS+KRT, one-sided | {onesided_min:.3f} | {str(onesided_min <= STRONG_RHO)} | {str(onesided_min <= LUSC_BOUND)} |
| ABS+KRT, two-sided | {twosided_min:.3f} | {str(twosided_min <= STRONG_RHO)} | {str(twosided_min <= LUSC_BOUND)} |
| ImmuneScore residual, one-sided | {immune_one_min:.3f} | {str(immune_one_min <= STRONG_RHO)} | {str(immune_one_min <= LUSC_BOUND)} |
| ImmuneScore residual, two-sided | {immune_two_min:.3f} | {str(immune_two_min <= STRONG_RHO)} | {str(immune_two_min <= LUSC_BOUND)} |

Largest |ρ| among primary four, ABS+KRT one-sided: {max_abs("ABS_KRT", "onesided_rank"):.3f}. Two-sided: {max_abs("ABS_KRT", "twosided_rank"):.3f}.

## Honest limits

1. **LUAD-only.** OncoSG (n=169, East-Asian surgical LUAD) and TCGA-LUSC are other cohorts. Their support does not transfer by residualization.
2. **n is not 1031.** Pooled LUAD+LUSC ESTIMATE n ≈ 1014–1034 (PR #89). This table is ~{n_abs} LUAD primaries with ABSOLUTE.
3. **KRT/EPCAM is a transcript score, not a cell fraction.** ESTIMATE TumorPurity is the published RNA tumor-fraction transform; ABSOLUTE is DNA. They are not interchangeable.
4. **xCell CD8 is TIMER2, not a from-scratch Xena xCell run.**
5. **GEP18 z-mean is not the Merck NanoString TIS.** ssGSEA is a sensitivity row.
6. **One-sided residual vs raw immune is not a partial correlation.** If it looks stronger than two-sided, composition in Y is the first explanation.
7. **No causality, no protein, no ICI outcome.**

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_luad_wave2.py
```

Downloads public tables into `data/` (gitignored) on first run. Signature lists are bundled under `data/signatures/`.

## Files

- `REPORT.md` — this writeup
- `DEFINITIONS.md` — locked definitions
- `correlations.tsv` — every model × estimator × feature
- `primary_results.tsv` — unadjusted / ABS / ABS+KRT / ImmuneScore × primary four
- `context_correlations.tsv` — TACSTD2 / purity / epithelium / ImmuneScore pairwise
- `sample_table.tsv` — per-sample values
- `gene_coverage.tsv`
- `summary.json` / `provenance.json`
- `figures/bar_residual_models.png`
- `figures/forest_double_residual.png`
- `figures/scatter_double_residual.png`
- `figures/scatter_tacstd2_vs_covariates.png`

## Data

- Expression: {URLS['expression']['url']}
- Purity: {URLS['purity']['url']}
- ESTIMATE: {URLS['estimate']['url']}
- TIMER2: {URLS['timer2']['url']}
"""

    (OUT / "REPORT.md").write_text(report)
    (OUT / "DEFINITIONS.md").write_text(
        report.split("## Primary results")[0]
        + "\n## Pointer\n\nNumbers live in `REPORT.md` and `correlations.tsv`. This file is the locked definition list.\n"
    )

    summary = {
        "task": "A1_luad_wave2",
        "cohort": "TCGA-LUAD primary tumors (-01)",
        "n_expression": n_expr,
        "n_absolute": n_abs,
        "n_abs_krt": n_abs_krt,
        "overall_verdict": overall_tag,
        "onesided_double_residual": {"tag": tag_one, "most_negative_rho": onesided_min, "text": txt_one},
        "twosided_double_residual": {"tag": tag_two, "most_negative_rho": twosided_min, "text": txt_two},
        "immunescore_residual_onesided": {"tag": tag_imm, "most_negative_rho": immune_one_min, "text": txt_imm},
        "immunescore_residual_twosided": {"tag": tag_imm2, "most_negative_rho": immune_two_min, "text": txt_imm2},
        "thresholds": {
            "strong_inverse": STRONG_RHO,
            "lusc_bound": LUSC_BOUND,
            "prior_deconv_abs_bound": DECONV_BOUND,
            "strong_inverse_source": "OncoSG A1 PR #139 CD8/GEP after purity (−0.309 / −0.349)",
        },
        "tacstd2_vs_absolute": {"rho": rho_tp, "p": p_tp, "n": n_tp},
        "tacstd2_vs_krt_epcam": {"rho": rho_tk, "p": p_tk, "n": n_tk},
        "tacstd2_vs_immunescore": {"rho": rho_ti, "p": p_ti, "n": n_ti},
        "absolute_vs_krt_epcam": {"rho": rho_ak, "p": p_ak, "n": n_ak},
        "xcell_cd8_column": xcell_cd8_col,
        "xcell_epithelial_column": xcell_epi_col,
        "krt_epcam_genes_present": epi_present,
        "primary_immune": PRIMARY_IMMUNE,
        "models": {k: v["cols"] for k, v in models.items()},
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    prov = {
        "claim": "A1 LUAD wave 2: TACSTD2 double residual on ABSOLUTE + KRT/EPCAM, and ImmuneScore residual, vs CD8/GEP18/CYT/xCell CD8",
        "inputs": {},
        "signatures": {
            "gep18": str(SIG / "gep18_genes.tsv"),
            "epithelial": str(SIG / "epithelial_krt_epcam.tsv"),
        },
        "randomness": "deterministic (no sampling)",
    }
    for k, v in URLS.items():
        p = DATA / v["file"]
        prov["inputs"][k] = {**v, "bytes": p.stat().st_size, "md5": md5(p)}
    with open(OUT / "provenance.json", "w") as f:
        json.dump(prov, f, indent=2)

    print("\n=== PRIMARY (ABS+KRT onesided_rank) ===")
    sub = res[(res.model == "ABS_KRT") & (res.estimator == "onesided_rank") & (res.feature.isin(PRIMARY_IMMUNE))]
    print(sub[["feature", "rho", "p", "n", "call"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n=== PRIMARY (ABS+KRT twosided_rank) ===")
    sub = res[(res.model == "ABS_KRT") & (res.estimator == "twosided_rank") & (res.feature.isin(PRIMARY_IMMUNE))]
    print(sub[["feature", "rho", "p", "n", "call"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n=== IMMUNE residual onesided_rank ===")
    sub = res[(res.model == "IMMUNE") & (res.estimator == "onesided_rank") & (res.feature.isin(PRIMARY_IMMUNE))]
    print(sub[["feature", "rho", "p", "n", "call"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nOverall verdict: {overall_tag}")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
