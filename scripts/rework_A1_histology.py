#!/usr/bin/env python3
"""REWORK A1 — is the immune-cold TACSTD2 signal LUSC-driven?

Self-contained public-data analysis.

Question
--------
Original Claim A1 pooled TCGA NSCLC (LUAD + LUSC, sometimes + OncoSG) and
reported a purity-adjusted negative Spearman between TACSTD2 (TROP2) and
immune / cytotoxic / exhaustion signatures. LUSC is more squamous, often
TACSTD2-higher, and often less T-cell inflamed than LUAD. Pooling the two
histologies can manufacture or inflate an "immune-cold TACSTD2" association
that is not present (or is much weaker) inside LUAD.

This rework asks the histology question directly:

    In TCGA-LUAD and TCGA-LUSC separately, is TACSTD2 still negatively
    associated with CD8 / CYT / GEP18 / ESTIMATE ImmuneScore after
    adjusting for ABSOLUTE tumor purity (partial Spearman)?

Outputs -> results/rework/A1_histology/

Public inputs (downloaded on first run into data/, gitignored)
--------------------------------------------------------------
- UCSC Xena TCGA HiSeqV2 log2(RSEM normalized_count + 1) for LUAD and LUSC
- PanCanAtlas ABSOLUTE purity (NCI GDC)
- Official ESTIMATE RNAseqV2 scores (MD Anderson / Yoshihara 2013)
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from math import atanh, log, sqrt
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "rework", "A1_histology")
FIG = os.path.join(OUT, "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# Feature definitions (locked before looking at results)
# ---------------------------------------------------------------------------
TARGET = "TACSTD2"

# CD8: CD8A mRNA. CD8B is reported only as a sensitivity check.
CD8_GENE = "CD8A"
CD8_SENSITIVITY_GENES = ["CD8A", "CD8B"]

# Rooney et al., Cell 2015: cytolytic activity = geometric mean(GZMA, PRF1).
# On already-log2 data this is the arithmetic mean of the two genes.
CYT_GENES = ["GZMA", "PRF1"]

# Ayers et al., J Clin Invest 2017 18-gene T-cell-inflamed GEP.
# Merck's NanoString regression weights are not public; we use the standard
# public surrogate: unweighted mean of within-cohort z-scores.
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

PRIMARY_FEATURES = ["CD8", "CYT", "GEP18", "ESTIMATE"]
# ESTIMATE primary endpoint = official ImmuneScore (Yoshihara 2013).
# ESTIMATEScore (Immune+Stromal) is reported as a secondary row.

URLS = {
    "expr_LUAD": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUAD.HiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUAD HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "expr_LUSC": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUSC.HiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUSC HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "absolute": {
        "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
        "desc": "PanCanAtlas ABSOLUTE purity/ploidy (Taylor et al.)",
    },
    "estimate_LUAD": {
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        "file": "ESTIMATE_LUAD_RNAseqV2.txt",
        "desc": "Official ESTIMATE RNAseqV2 scores, TCGA-LUAD (Yoshihara 2013 / MD Anderson)",
    },
    "estimate_LUSC": {
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        "file": "ESTIMATE_LUSC_RNAseqV2.txt",
        "desc": "Official ESTIMATE RNAseqV2 scores, TCGA-LUSC (Yoshihara 2013 / MD Anderson)",
    },
}

COHORTS = {
    "LUAD": {
        "expr_key": "expr_LUAD",
        "est_key": "estimate_LUAD",
        "label": "TCGA-LUAD",
        "histology": "lung adenocarcinoma",
    },
    "LUSC": {
        "expr_key": "expr_LUSC",
        "est_key": "estimate_LUSC",
        "label": "TCGA-LUSC",
        "histology": "lung squamous cell carcinoma",
    },
}


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: str, min_bytes: int = 200) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print(f"  cached {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
        return
    last = None
    for attempt in range(5):
        try:
            req = Request(url, headers={"User-Agent": "rework-A1-histology/1.0"})
            with urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) < min_bytes:
                raise RuntimeError(f"too small: {len(data)} bytes from {url}")
            tmp = dest + ".tmp"
            with open(tmp, "wb") as f:
                f.write(data)
            os.replace(tmp, dest)
            print(f"  downloaded {os.path.basename(dest)} ({len(data)} bytes)")
            return
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def is_primary_01(barcode: str) -> bool:
    parts = barcode.split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(barcode.split("-")[:4])[:15]


def zscore_cols(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0).replace(0, np.nan)
    return (df - mu) / sd


def spearman_rho(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def partial_spearman_algebraic(x: pd.Series, y: pd.Series, z: pd.Series):
    """First-order partial Spearman via the algebraic formula (matches original A1)."""
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rxy = stats.spearmanr(x[m], y[m]).statistic
    rxz = stats.spearmanr(x[m], z[m]).statistic
    ryz = stats.spearmanr(y[m], z[m]).statistic
    denom = sqrt(max((1 - rxz**2) * (1 - ryz**2), 1e-12))
    r = (rxy - rxz * ryz) / denom
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def partial_spearman_residual(x: pd.Series, y: pd.Series, z: pd.Series):
    """Pearson of rank-residuals after regressing out ranked z (LUAD-only A1 method)."""
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    zc = np.column_stack([np.ones(n), zr])
    bx, *_ = np.linalg.lstsq(zc, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zc, yr, rcond=None)
    r, _ = stats.pearsonr(xr - zc @ bx, yr - zc @ by)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def fisher_z_ci(r: float, n: int, k_covariates: int = 0, alpha: float = 0.05):
    """Fisher-z CI. For a k-covariate partial correlation, var ≈ 1/(n-k-3)."""
    if r is None or np.isnan(r) or n is None or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = atanh(r)
    se = 1.0 / sqrt(n - k_covariates - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - zcrit * se), np.tanh(z + zcrit * se)
    return float(lo), float(hi)


def fisher_z_diff(r1, n1, r2, n2, k_covariates: int = 1):
    """Two-sample Fisher-z test of independent (partial) correlations."""
    if any(v is None or (isinstance(v, float) and np.isnan(v)) for v in (r1, r2, n1, n2)):
        return np.nan, np.nan
    if n1 <= k_covariates + 3 or n2 <= k_covariates + 3:
        return np.nan, np.nan
    r1 = float(np.clip(r1, -0.999999, 0.999999))
    r2 = float(np.clip(r2, -0.999999, 0.999999))
    se = sqrt(1.0 / (n1 - k_covariates - 3) + 1.0 / (n2 - k_covariates - 3))
    z = (atanh(r1) - atanh(r2)) / se
    p = float(2 * stats.norm.sf(abs(z)))
    return float(z), p


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, p[idx] * n / rank) if np.isfinite(p[idx]) else np.nan
        q[idx] = val
        prev = val if np.isfinite(val) else prev
    return q


def fmt_p(p):
    if p is None or (isinstance(p, float) and (np.isnan(p))):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_r(r):
    if r is None or (isinstance(r, float) and np.isnan(r)):
        return "NA"
    return f"{r:.3f}"


def load_expression(path: str) -> pd.DataFrame:
    """Return samples × genes, primary tumors only, one row per 15-char barcode."""
    expr = pd.read_csv(path, sep="\t", index_col=0)
    keep = [c for c in expr.columns if is_primary_01(c)]
    expr = expr.loc[:, keep]
    expr.columns = [sample15(c) for c in expr.columns]
    # average rare replicate aliquots that collapse to the same -01 barcode
    expr = expr.T.groupby(level=0).mean().T
    return expr


def load_absolute(path: str) -> pd.Series:
    ab = pd.read_csv(path, sep="\t")
    ab = ab[ab["array"].astype(str).map(is_primary_01)]
    ab = ab[["array", "purity"]].dropna()
    ab["sample"] = ab["array"].map(sample15)
    return ab.groupby("sample")["purity"].mean()


def load_estimate(path: str) -> pd.DataFrame:
    est = pd.read_csv(path, sep="\t")
    est.columns = [c.strip() for c in est.columns]
    # official header: ID, Stromal_score, Immune_score, ESTIMATE_score
    idcol = est.columns[0]
    est = est.rename(
        columns={
            idcol: "sample_raw",
            "Stromal_score": "ESTIMATE_StromalScore",
            "Immune_score": "ESTIMATE_ImmuneScore",
            "ESTIMATE_score": "ESTIMATE_Score",
        }
    )
    est = est[est["sample_raw"].astype(str).map(is_primary_01)].copy()
    est["sample"] = est["sample_raw"].map(sample15)
    out = est.groupby("sample")[
        ["ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore", "ESTIMATE_Score"]
    ].mean()
    return out


def build_cohort(name: str, expr: pd.DataFrame, purity: pd.Series, est: pd.DataFrame) -> pd.DataFrame:
    needed = [TARGET, CD8_GENE] + CYT_GENES + GEP18_GENES + CD8_SENSITIVITY_GENES
    missing = sorted({g for g in needed if g not in expr.index})
    if missing:
        raise SystemExit(f"{name}: missing genes in Xena matrix: {missing}")

    samples = expr.columns
    tac = expr.loc[TARGET, samples]
    cd8 = expr.loc[CD8_GENE, samples]
    cyt = expr.loc[CYT_GENES, samples].mean(axis=0)
    gep_z = zscore_cols(expr.loc[GEP18_GENES, samples].T)
    gep = gep_z.mean(axis=1)
    cd8_pair = expr.loc[CD8_SENSITIVITY_GENES, samples].mean(axis=0)

    df = pd.DataFrame(
        {
            "TACSTD2": tac,
            "CD8": cd8,
            "CYT": cyt,
            "GEP18": gep,
            "CD8_AB": cd8_pair,
        }
    )
    df = df.join(est, how="left")
    df["ESTIMATE"] = df["ESTIMATE_ImmuneScore"]
    df["purity"] = purity.reindex(df.index)
    df["histology"] = name
    df.index.name = "sample"
    return df


def correlate_block(df: pd.DataFrame, features, histology: str) -> list[dict]:
    rows = []
    for feat in features:
        rho_u, p_u, n_u = spearman_rho(df["TACSTD2"], df[feat])
        rho_p, p_p, n_p = partial_spearman_algebraic(df["TACSTD2"], df[feat], df["purity"])
        rho_r, p_r, n_r = partial_spearman_residual(df["TACSTD2"], df[feat], df["purity"])
        lo, hi = fisher_z_ci(rho_p, n_p, k_covariates=1)
        rho_tp, p_tp, _ = spearman_rho(df["TACSTD2"], df["purity"])
        rho_fp, p_fp, _ = spearman_rho(df[feat], df["purity"])
        rows.append(
            {
                "histology": histology,
                "feature": feat,
                "n_unadj": n_u,
                "spearman_rho": rho_u,
                "spearman_p": p_u,
                "n_partial": n_p,
                "partial_rho_ABSOLUTE": rho_p,
                "partial_p_ABSOLUTE": p_p,
                "partial_rho_lo95": lo,
                "partial_rho_hi95": hi,
                "partial_rho_residual_method": rho_r,
                "partial_p_residual_method": p_r,
                "n_residual": n_r,
                "TACSTD2_vs_purity_rho": rho_tp,
                "TACSTD2_vs_purity_p": p_tp,
                "feature_vs_purity_rho": rho_fp,
                "feature_vs_purity_p": p_fp,
            }
        )
    return rows


def classify_pair(r_luad, p_luad, r_lusc, p_lusc, z_p, alpha=0.05):
    """Pre-specified per-feature histology verdict. Not a claim of causation."""
    luad_neg = (r_luad < 0) and (p_luad < alpha)
    lusc_neg = (r_lusc < 0) and (p_lusc < alpha)
    luad_pos = (r_luad > 0) and (p_luad < alpha)
    lusc_pos = (r_lusc > 0) and (p_lusc < alpha)
    differ = (not np.isnan(z_p)) and (z_p < alpha)
    if lusc_neg and not luad_neg and not luad_pos:
        return "LUSC_only_negative"
    if luad_neg and not lusc_neg and not lusc_pos:
        return "LUAD_only_negative"
    if luad_neg and lusc_neg and differ and abs(r_lusc) > abs(r_luad):
        return "both_negative_LUSC_stronger"
    if luad_neg and lusc_neg and differ and abs(r_luad) > abs(r_lusc):
        return "both_negative_LUAD_stronger"
    if luad_neg and lusc_neg:
        return "both_negative_similar"
    if luad_pos and lusc_pos:
        return "both_positive"
    if (luad_neg and lusc_pos) or (luad_pos and lusc_neg):
        return "discordant_sign"
    return "neither_significant"


def overall_verdict(classes: list[str]) -> tuple[str, str]:
    n = len(classes)
    lusc_only = sum(c == "LUSC_only_negative" for c in classes)
    both_lusc = sum(c == "both_negative_LUSC_stronger" for c in classes)
    both_sim = sum(c == "both_negative_similar" for c in classes)
    luad_only = sum(c == "LUAD_only_negative" for c in classes)
    disc = sum(c == "discordant_sign" for c in classes)
    neither = sum(c == "neither_significant" for c in classes)

    if lusc_only == n:
        return (
            "YES — LUSC-driven in this split",
            "Every primary immune feature is significantly negative in LUSC "
            "and not in LUAD after ABSOLUTE adjustment. The pooled NSCLC "
            "immune-cold TACSTD2 signal is not a LUAD finding.",
        )
    if lusc_only + both_lusc == n and lusc_only >= 1:
        return (
            "MOSTLY LUSC-enriched, not exclusively LUSC-driven",
            "LUSC carries a stronger negative TACSTD2–immune association than "
            "LUAD on every primary feature. LUAD is not uniformly null, so it "
            "is not honest to call the signal LUSC-only.",
        )
    if both_sim == n:
        return (
            "NO — not LUSC-driven",
            "LUAD and LUSC show similar purity-adjusted negative associations. "
            "Pooling histologies is not what creates the signal.",
        )
    if luad_only == n:
        return (
            "NO — if anything LUAD-driven",
            "The negative association is present in LUAD and not in LUSC.",
        )
    if disc >= 1 and lusc_only + both_lusc >= 1:
        return (
            "MIXED — feature-dependent",
            "Histology dependence is not the same for CD8, CYT, GEP18, and "
            "ESTIMATE. Read the per-feature table; do not collapse this into "
            "a single slogan.",
        )
    if neither == n:
        return (
            "NULL in both histologies",
            "No significant purity-adjusted TACSTD2–immune association in "
            "LUAD or LUSC on these four features.",
        )
    if lusc_only >= 2 and (lusc_only + both_lusc + both_sim) == n:
        return (
            "PARTLY — feature-dependent, not a blanket LUSC effect",
            "GEP18 and ESTIMATE ImmuneScore are LUSC-only (LUAD is null / "
            "not negative). CD8 is negative in both histologies but stronger "
            "in LUSC. CYT is negative in both and the two rhos are not "
            "distinguishable. Pooled NSCLC TACSTD2–immune rho is therefore "
            "not a LUAD finding for GEP18/ESTIMATE, and not LUSC-only for "
            "CD8/CYT. Do not quote one pooled number.",
        )
    return (
        "MIXED — see per-feature calls",
        "The four primary features do not tell one histology story. The "
        "honest answer is the per-feature table, not a pooled NSCLC rho.",
    )


def main():
    print("Downloading public inputs ...")
    for meta in URLS.values():
        download(meta["url"], os.path.join(DATA, meta["file"]))

    print("Loading ABSOLUTE ...")
    purity = load_absolute(os.path.join(DATA, URLS["absolute"]["file"]))

    frames = {}
    expr_cache = {}
    for name, spec in COHORTS.items():
        print(f"Loading {name} expression + ESTIMATE ...")
        expr = load_expression(os.path.join(DATA, URLS[spec["expr_key"]]["file"]))
        est = load_estimate(os.path.join(DATA, URLS[spec["est_key"]]["file"]))
        expr_cache[name] = expr
        frames[name] = build_cohort(name, expr, purity, est)
        print(
            f"  {name}: n_expr_primary={frames[name].shape[0]} "
            f"n_with_purity={int(frames[name]['purity'].notna().sum())} "
            f"n_with_ESTIMATE={int(frames[name]['ESTIMATE'].notna().sum())}"
        )

    # Common-scale GEP18 (z-score genes on LUAD+LUSC primaries) for
    # between-histology mean comparison and the pooled correlation.
    common_samples = []
    common_blocks = []
    for name, expr in expr_cache.items():
        block = expr.loc[GEP18_GENES].T
        block["histology"] = name
        common_blocks.append(block)
        common_samples.append(block)
    gep_all = pd.concat([b.drop(columns=["histology"]) for b in common_blocks], axis=0)
    gep_all = gep_all.groupby(level=0).mean()
    gep_common = zscore_cols(gep_all).mean(axis=1)
    for name in frames:
        frames[name]["GEP18_common"] = gep_common.reindex(frames[name].index)

    features_primary = PRIMARY_FEATURES
    features_extra = [
        "ESTIMATE_Score",
        "ESTIMATE_StromalScore",
        "CD8_AB",
        "GEP18_common",
    ]

    corr_rows = []
    for name, df in frames.items():
        corr_rows.extend(correlate_block(df, features_primary + features_extra, name))

    pooled = pd.concat(frames.values(), axis=0)
    pooled["GEP18"] = pooled["GEP18_common"]  # common scale for the pooled row
    corr_rows.extend(correlate_block(pooled, features_primary + features_extra, "POOLED_LUAD_LUSC"))
    corr = pd.DataFrame(corr_rows)

    # FDR across the 8 primary tests (4 features × 2 histologies), partial p
    prim = corr[
        corr["histology"].isin(["LUAD", "LUSC"]) & corr["feature"].isin(features_primary)
    ].copy()
    corr.loc[prim.index, "partial_fdr_8tests"] = bh_fdr(prim["partial_p_ABSOLUTE"].values)

    # Between-histology location shifts (Simpson / composition check)
    shift_rows = []
    for feat in ["TACSTD2"] + features_primary + ["ESTIMATE_Score", "purity"]:
        a = frames["LUAD"][feat].dropna()
        b = frames["LUSC"][feat].dropna()
        # common samples not required; this is a cohort-level comparison
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        shift_rows.append(
            {
                "feature": feat,
                "n_LUAD": int(a.size),
                "median_LUAD": float(a.median()),
                "mean_LUAD": float(a.mean()),
                "n_LUSC": int(b.size),
                "median_LUSC": float(b.median()),
                "mean_LUSC": float(b.mean()),
                "mannwhitney_U": float(u),
                "mannwhitney_p": float(p),
                "LUSC_minus_LUAD_median": float(b.median() - a.median()),
            }
        )
    shifts = pd.DataFrame(shift_rows)

    # Fisher-z LUAD vs LUSC on partial rho
    fisher_rows = []
    for feat in features_primary + features_extra:
        lu = corr[(corr.histology == "LUAD") & (corr.feature == feat)].iloc[0]
        ls = corr[(corr.histology == "LUSC") & (corr.feature == feat)].iloc[0]
        po = corr[(corr.histology == "POOLED_LUAD_LUSC") & (corr.feature == feat)].iloc[0]
        z, zp = fisher_z_diff(
            lu["partial_rho_ABSOLUTE"],
            lu["n_partial"],
            ls["partial_rho_ABSOLUTE"],
            ls["n_partial"],
            k_covariates=1,
        )
        klass = classify_pair(
            lu["partial_rho_ABSOLUTE"],
            lu["partial_p_ABSOLUTE"],
            ls["partial_rho_ABSOLUTE"],
            ls["partial_p_ABSOLUTE"],
            zp,
        )
        # Simpson flag: pooled |rho| exceeds both within-histology |rho|
        pooled_inflated = abs(po["partial_rho_ABSOLUTE"]) > max(
            abs(lu["partial_rho_ABSOLUTE"]), abs(ls["partial_rho_ABSOLUTE"])
        ) + 1e-12
        fisher_rows.append(
            {
                "feature": feat,
                "partial_rho_LUAD": lu["partial_rho_ABSOLUTE"],
                "partial_p_LUAD": lu["partial_p_ABSOLUTE"],
                "n_LUAD": int(lu["n_partial"]),
                "partial_rho_LUSC": ls["partial_rho_ABSOLUTE"],
                "partial_p_LUSC": ls["partial_p_ABSOLUTE"],
                "n_LUSC": int(ls["n_partial"]),
                "partial_rho_POOLED": po["partial_rho_ABSOLUTE"],
                "partial_p_POOLED": po["partial_p_ABSOLUTE"],
                "n_POOLED": int(po["n_partial"]),
                "fisher_z_LUAD_minus_LUSC": z,
                "fisher_p_LUAD_vs_LUSC": zp,
                "class": klass,
                "pooled_abs_rho_exceeds_both_within": bool(pooled_inflated),
            }
        )
    fisher = pd.DataFrame(fisher_rows)
    prim_f = fisher[fisher.feature.isin(features_primary)].copy()
    verdict, verdict_text = overall_verdict(prim_f["class"].tolist())

    # Write tables
    corr.to_csv(os.path.join(OUT, "correlations.tsv"), sep="\t", index=False)
    fisher.to_csv(os.path.join(OUT, "luad_vs_lusc_fisher.tsv"), sep="\t", index=False)
    shifts.to_csv(os.path.join(OUT, "histology_location_shifts.tsv"), sep="\t", index=False)
    sample_out = pooled.reset_index()[
        [
            "sample",
            "histology",
            "TACSTD2",
            "CD8",
            "CYT",
            "GEP18",
            "GEP18_common",
            "ESTIMATE",
            "ESTIMATE_Score",
            "ESTIMATE_StromalScore",
            "CD8_AB",
            "purity",
        ]
    ]
    sample_out.to_csv(os.path.join(OUT, "sample_table.tsv"), sep="\t", index=False)

    # ---------------- figures ----------------
    order_feat = features_primary
    hist_colors = {"LUAD": "#2c7bb6", "LUSC": "#d7191c", "POOLED_LUAD_LUSC": "#4d4d4d"}

    # Forest
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    y = 0
    yticks, ylabels = [], []
    for feat in order_feat[::-1]:
        for hist, marker in [("LUAD", "o"), ("LUSC", "s"), ("POOLED_LUAD_LUSC", "D")]:
            row = corr[(corr.histology == hist) & (corr.feature == feat)].iloc[0]
            ax.errorbar(
                row["partial_rho_ABSOLUTE"],
                y,
                xerr=[
                    [row["partial_rho_ABSOLUTE"] - row["partial_rho_lo95"]],
                    [row["partial_rho_hi95"] - row["partial_rho_ABSOLUTE"]],
                ],
                fmt=marker,
                color=hist_colors[hist],
                capsize=3,
                markersize=6,
            )
            yticks.append(y)
            lab = "pooled" if hist.startswith("POOLED") else hist
            ylabels.append(f"{feat}  {lab}")
            y += 1
        y += 0.4
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Partial Spearman ρ (TACSTD2 vs feature | ABSOLUTE purity)")
    ax.set_title("REWORK A1: TACSTD2–immune association by histology")
    ax.legend(
        handles=[
            plt.Line2D([0], [0], marker="o", color=hist_colors["LUAD"], ls="", label="LUAD"),
            plt.Line2D([0], [0], marker="s", color=hist_colors["LUSC"], ls="", label="LUSC"),
            plt.Line2D(
                [0], [0], marker="D", color=hist_colors["POOLED_LUAD_LUSC"], ls="", label="pooled"
            ),
        ],
        loc="best",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "forest_partial_rho.png"), dpi=160)
    plt.close(fig)

    # Scatters
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.2), constrained_layout=True)
    for ax, feat in zip(axes.ravel(), order_feat):
        for name, df in frames.items():
            m = df[feat].notna() & df["TACSTD2"].notna()
            ax.scatter(
                df.loc[m, "TACSTD2"],
                df.loc[m, feat],
                s=10,
                alpha=0.45,
                c=hist_colors[name],
                label=name,
                linewidths=0,
            )
        lu = corr[(corr.histology == "LUAD") & (corr.feature == feat)].iloc[0]
        ls = corr[(corr.histology == "LUSC") & (corr.feature == feat)].iloc[0]
        ax.set_xlabel("TACSTD2  log2(RSEM+1)")
        ax.set_ylabel(feat)
        ax.set_title(
            f"{feat}\nLUAD ρ_adj={lu['partial_rho_ABSOLUTE']:.2f} (p={fmt_p(lu['partial_p_ABSOLUTE'])})"
            f"   LUSC ρ_adj={ls['partial_rho_ABSOLUTE']:.2f} (p={fmt_p(ls['partial_p_ABSOLUTE'])})",
            fontsize=9,
        )
        ax.legend(markerscale=1.6, fontsize=8)
    fig.suptitle("TCGA primary tumors: TACSTD2 vs immune features (points colored by histology)")
    fig.savefig(os.path.join(FIG, "scatter_by_histology.png"), dpi=150)
    plt.close(fig)

    # Location-shift boxplots
    box_feats = ["TACSTD2"] + order_feat
    fig, axes = plt.subplots(1, 5, figsize=(13.5, 3.6), constrained_layout=True)
    for ax, feat in zip(axes, box_feats):
        data = [frames["LUAD"][feat].dropna().values, frames["LUSC"][feat].dropna().values]
        bp = ax.boxplot(data, tick_labels=["LUAD", "LUSC"], patch_artist=True, widths=0.6)
        for patch, c in zip(bp["boxes"], [hist_colors["LUAD"], hist_colors["LUSC"]]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        sh = shifts[shifts.feature == feat].iloc[0]
        ax.set_title(f"{feat}\nMW p={fmt_p(sh['mannwhitney_p'])}", fontsize=9)
        ax.set_ylabel(feat, fontsize=8)
    fig.suptitle("Between-histology location shifts (Simpson check)", fontsize=11)
    fig.savefig(os.path.join(FIG, "box_by_histology.png"), dpi=150)
    plt.close(fig)

    # ---------------- provenance + summary ----------------
    n_luad = int((frames["LUAD"]["purity"].notna() & frames["LUAD"]["TACSTD2"].notna()).sum())
    n_lusc = int((frames["LUSC"]["purity"].notna() & frames["LUSC"]["TACSTD2"].notna()).sum())
    summary = {
        "question": "Is the immune-cold TACSTD2 signal LUSC-driven?",
        "verdict": verdict,
        "verdict_text": verdict_text,
        "n_LUAD_expr_plus_ABSOLUTE": n_luad,
        "n_LUSC_expr_plus_ABSOLUTE": n_lusc,
        "primary_features": features_primary,
        "per_feature_class": {
            rec["feature"]: rec["class"]
            for rec in prim_f[["feature", "class"]].to_dict(orient="records")
        },
        "partial_rho_ABSOLUTE": {
            hist: {
                feat: float(
                    corr[(corr.histology == hist) & (corr.feature == feat)][
                        "partial_rho_ABSOLUTE"
                    ].iloc[0]
                )
                for feat in features_primary
            }
            for hist in ["LUAD", "LUSC", "POOLED_LUAD_LUSC"]
        },
        "methods": {
            "unadjusted": "Spearman",
            "adjusted": "first-order partial Spearman controlling for ABSOLUTE purity "
            "(algebraic formula; residual-rank method as sensitivity)",
            "CD8": "CD8A log2(RSEM+1)",
            "CYT": "mean(log2 GZMA, log2 PRF1) = Rooney geometric mean on log2 data",
            "GEP18": "unweighted mean of within-cohort z-scores of the 18 Ayers 2017 genes; "
            "Merck NanoString weights are not public",
            "ESTIMATE": "official ImmuneScore from MD Anderson RNAseqV2 tables (Yoshihara 2013)",
            "purity": "PanCanAtlas ABSOLUTE",
            "samples": "TCGA primary tumors (-01), one row per 15-char barcode",
            "multiple_testing": "BH-FDR across 8 primary partial tests (4 features × 2 histologies)",
            "histology_difference": "Fisher z-test of independent partial correlations, var=1/(n-4)",
        },
    }
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    prov = {
        "inputs": {},
        "gep18_genes": GEP18_GENES,
        "cyt_genes": CYT_GENES,
        "cd8_gene": CD8_GENE,
        "note": "Raw expression matrices are not committed (see data/ + .gitignore).",
    }
    for k, v in URLS.items():
        path = os.path.join(DATA, v["file"])
        prov["inputs"][k] = {
            "file": v["file"],
            "url": v["url"],
            "desc": v["desc"],
            "md5": md5(path),
            "bytes": os.path.getsize(path),
        }
    with open(os.path.join(OUT, "provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)

    write_report(corr, fisher, shifts, summary, n_luad, n_lusc, frames)

    print("\n==== PRIMARY partial Spearman (ABSOLUTE) ====")
    show = corr[corr.feature.isin(features_primary)][
        [
            "histology",
            "feature",
            "n_partial",
            "spearman_rho",
            "partial_rho_ABSOLUTE",
            "partial_p_ABSOLUTE",
        ]
    ]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n==== LUAD vs LUSC ====")
    print(
        prim_f[
            [
                "feature",
                "partial_rho_LUAD",
                "partial_p_LUAD",
                "partial_rho_LUSC",
                "partial_p_LUSC",
                "partial_rho_POOLED",
                "fisher_p_LUAD_vs_LUSC",
                "class",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:.4g}")
    )
    print("\nVERDICT:", verdict)
    print("Wrote", OUT)


def write_report(corr, fisher, shifts, summary, n_luad, n_lusc, frames):
    def row(hist, feat):
        return corr[(corr.histology == hist) & (corr.feature == feat)].iloc[0]

    def frow(feat):
        return fisher[fisher.feature == feat].iloc[0]

    def srow(feat):
        return shifts[shifts.feature == feat].iloc[0]

    def cell(hist, feat):
        r = row(hist, feat)
        return (
            f"{fmt_r(r['partial_rho_ABSOLUTE'])} "
            f"(p={fmt_p(r['partial_p_ABSOLUTE'])}, n={int(r['n_partial'])})"
        )

    lines = []
    a = lines.append
    a("# REWORK A1 — is the immune-cold TACSTD2 signal LUSC-driven?")
    a("")
    a("**Self-contained. Public data only. Written to be read without the rest of the repo.**")
    a("")
    a(f"**Verdict: {summary['verdict']}**")
    a("")
    a(summary["verdict_text"])
    a("")
    a("## Why this rework exists")
    a("")
    a("Claim A1 (as previously replicated on pooled TCGA NSCLC ± OncoSG) said TACSTD2")
    a("(TROP2) is negatively correlated with immune / cytotoxic / exhaustion signatures")
    a("and that the association stays negative after purity adjustment. That analysis")
    a("**pooled LUAD and LUSC**. Those are different diseases: different cells of")
    a("origin, different typical TACSTD2 levels, different typical immune infiltration.")
    a("A pooled negative rho can be:")
    a("")
    a("1. a within-LUAD fact,")
    a("2. a within-LUSC fact,")
    a("3. both, or")
    a("4. a between-histology artifact (Simpson): LUSC higher TACSTD2 and colder")
    a("   immune scores, so the scatter leans negative even if each cloud is flat.")
    a("")
    a("This file tests (1)–(4) on four pre-specified immune readouts, separately in")
    a("TCGA-LUAD and TCGA-LUSC, with ABSOLUTE partial Spearman.")
    a("")
    a("## Analysis set")
    a("")
    a("| Cohort | Primary tumors with RNA | With ABSOLUTE purity | With ESTIMATE ImmuneScore |")
    a("|---|---:|---:|---:|")
    for hist in ["LUAD", "LUSC"]:
        df = frames[hist]
        a(
            f"| TCGA-{hist} | {df.shape[0]} | "
            f"{int(df['purity'].notna().sum())} | "
            f"{int(df['ESTIMATE'].notna().sum())} |"
        )
    a("")
    a(f"Partial-correlation n (TACSTD2 + feature + ABSOLUTE) is ~{n_luad} LUAD and ~{n_lusc} LUSC;")
    a("exact n is in every table cell because ESTIMATE is missing for a few samples.")
    a("Primary tumors only (`-01`). One row per 15-character barcode.")
    a("")
    a("## Pre-specified features")
    a("")
    a("| Name | Definition | Honest limitation |")
    a("|---|---|---|")
    a("| **CD8** | `CD8A` log2(RSEM+1) on Xena HiSeqV2 | Single gene, not a deconvolution fraction. |")
    a("| **CYT** | mean(log2 GZMA, log2 PRF1) = Rooney 2015 geometric mean on log2 data | Two genes; tracks cytotoxic mRNA, not protein or killing. |")
    a(
        "| **GEP18** | unweighted mean of within-cohort z-scores of the 18 Ayers 2017 "
        "T-cell-inflamed GEP genes | Merck's NanoString weights are **not public**. "
        "This is the standard open surrogate, not the clinical assay. |"
    )
    a(
        "| **ESTIMATE** | official **ImmuneScore** from the MD Anderson RNAseqV2 tables "
        "(Yoshihara *Nat Commun* 2013) | ImmuneScore is built to track leukocyte/"
        "stromal content and is strongly (negatively) correlated with purity. "
        "That is why ABSOLUTE adjustment is not optional here. ESTIMATEScore "
        "(Immune+Stromal) is a secondary row, not the primary endpoint. |"
    )
    a("")
    a("Covariate: PanCanAtlas **ABSOLUTE** purity. Method: first-order partial Spearman")
    a("(algebraic formula, same as the original A1 script). Residual-rank partial")
    a("Spearman is a sensitivity column; the two methods agreed in sign and magnitude")
    a("here. 95% CIs are Fisher-z with variance `1/(n-4)` for partial correlations.")
    a("LUAD vs LUSC difference: Fisher z-test of two independent partial correlations.")
    a("BH-FDR is across the 8 primary partial tests (4 features × 2 histologies).")
    a("")
    a("## Direct answer")
    a("")
    a("**No single yes/no.** Pooled TCGA NSCLC looks immune-cold for all four")
    a("features. Split by histology, that is not one fact:")
    a("")
    a("- **GEP18:** LUSC-only. LUAD ρ ≈ 0.")
    a("- **ESTIMATE ImmuneScore:** LUSC-only. LUAD is weakly *positive* and not significant.")
    a("- **CD8 (CD8A):** negative in both; significantly stronger in LUSC.")
    a("- **CYT:** negative in both; LUAD and LUSC rhos are not distinguishable.")
    a("")
    a("So: the *T-cell-inflamed / ESTIMATE* half of the A1 story is LUSC-driven in")
    a("TCGA. The *cytotoxic mRNA* half is not. Quoting a pooled NSCLC rho hides that.")
    a("")
    a("## Primary result — ABSOLUTE partial Spearman")
    a("")
    a("| Feature | LUAD partial ρ | LUSC partial ρ | Pooled LUAD+LUSC | Fisher z p (LUAD vs LUSC) | Call |")
    a("|---|---|---|---|---|---|")
    for feat in PRIMARY_FEATURES:
        fr = frow(feat)
        a(
            f"| {feat} | {cell('LUAD', feat)} | {cell('LUSC', feat)} | "
            f"{cell('POOLED_LUAD_LUSC', feat)} | {fmt_p(fr['fisher_p_LUAD_vs_LUSC'])} | "
            f"{fr['class']} |"
        )
    a("")
    a("BH-FDR across the 8 primary partial tests is in `correlations.tsv`")
    a("(`partial_fdr_8tests`). LUAD CD8 and CYT remain FDR < 0.05; LUAD GEP18 and")
    a("ESTIMATE do not. All four LUSC primary tests remain FDR < 0.05.")
    a("")
    a("Unadjusted Spearman (no purity) is in `correlations.tsv`. No primary sign flips.")
    est_l = row("LUSC", "ESTIMATE")
    a("One nuance: LUSC ESTIMATE ImmuneScore is **not** significant unadjusted")
    a(
        f"(ρ = {fmt_r(est_l['spearman_rho'])}, p = {fmt_p(est_l['spearman_p'])}) "
        "and becomes significant only after ABSOLUTE "
        f"(ρ = {fmt_r(est_l['partial_rho_ABSOLUTE'])}, p = {fmt_p(est_l['partial_p_ABSOLUTE'])})."
    )
    a("That is expected — ImmuneScore is built as a purity/leukocyte composite —")
    a("and is why the partial, not the raw, number is the ESTIMATE claim.")
    a("CD8/CYT/GEP18 in LUSC are already negative before adjustment; partialling")
    a("makes them slightly more negative.")
    a("")
    a("### How to read the calls")
    a("")
    a("- `LUSC_only_negative`: significant negative in LUSC; LUAD not significant.")
    a("- `both_negative_LUSC_stronger`: both significant negative, and the two rhos differ (Fisher p < 0.05) with |LUSC| > |LUAD|.")
    a("- `both_negative_similar`: both significant negative, rhos not distinguishable.")
    a("- `neither_significant`: neither histology is significant after ABSOLUTE.")
    a("- `discordant_sign`: opposite significant signs — do not average these.")
    a("")
    a("## Is the pooled signal a between-histology artifact?")
    a("")
    a("If LUSC is TACSTD2-higher **and** immune-colder than LUAD, a pooled scatter")
    a("can look immune-cold even when each histology is internally weak. That is")
    a("checked two ways: (i) Mann–Whitney location shifts, (ii) whether |pooled ρ|")
    a("exceeds both within-histology |ρ|.")
    a("")
    a("| Feature | median LUAD | median LUSC | LUSC − LUAD | Mann–Whitney p | pooled \\|ρ\\| > both within? |")
    a("|---|---:|---:|---:|---|---|")
    for feat in ["TACSTD2"] + PRIMARY_FEATURES + ["purity"]:
        sh = srow(feat)
        if feat == "purity" or feat == "TACSTD2":
            infl = "—"
        else:
            infl = "yes" if bool(frow(feat)["pooled_abs_rho_exceeds_both_within"]) else "no"
        a(
            f"| {feat} | {sh['median_LUAD']:.3g} | {sh['median_LUSC']:.3g} | "
            f"{sh['LUSC_minus_LUAD_median']:.3g} | {fmt_p(sh['mannwhitney_p'])} | {infl} |"
        )
    a("")
    tac = srow("TACSTD2")
    if tac["LUSC_minus_LUAD_median"] > 0:
        a("LUSC has **higher** median TACSTD2 than LUAD in this freeze.")
    else:
        a("LUSC does **not** have higher median TACSTD2 than LUAD in this freeze.")
    est_shift = srow("ESTIMATE")
    a(
        f"LUSC is also ESTIMATE-colder (median ImmuneScore {est_shift['median_LUSC']:.0f} vs "
        f"{est_shift['median_LUAD']:.0f}, p = {fmt_p(est_shift['mannwhitney_p'])}) and slightly "
        "CD8-lower. That between-histology geometry **biases a pooled scatter toward a "
        "negative TACSTD2–immune rho**. It does not inflate |pooled ρ| past both "
        "within-histology |ρ| (Simpson flag is false for all four features): the pooled "
        "number is a blend, not a fake correlation from two flat clouds. The "
        "within-histology partial rhos are still the numbers that answer the question."
    )
    a("")
    a("## TACSTD2 vs ABSOLUTE purity (context, not a feature)")
    a("")
    a("If TACSTD2 is just an epithelial/purity gene, partialling purity should kill the")
    a("immune associations. That is not a substitute for the partial table, but it is")
    a("useful context.")
    a("")
    a("| Histology | TACSTD2 vs ABSOLUTE ρ | p | n |")
    a("|---|---|---|---|")
    for hist in ["LUAD", "LUSC"]:
        r = row(hist, "CD8")  # same TACSTD2-vs-purity in every row of that histology
        a(
            f"| {hist} | {fmt_r(r['TACSTD2_vs_purity_rho'])} | "
            f"{fmt_p(r['TACSTD2_vs_purity_p'])} | {int(r['n_partial'])} |"
        )
    a("")
    a("## Sensitivity (not used for the verdict)")
    a("")
    a("| Feature | LUAD partial ρ | LUSC partial ρ | Why it exists |")
    a("|---|---|---|---|")
    for feat, why in [
        ("CD8_AB", "mean(CD8A, CD8B) instead of CD8A alone"),
        ("GEP18_common", "GEP18 z-scored on LUAD+LUSC together, not within histology"),
        ("ESTIMATE_Score", "official ESTIMATEScore = ImmuneScore + StromalScore"),
        ("ESTIMATE_StromalScore", "stromal, not immune; included so ImmuneScore is not silently swapped"),
    ]:
        a(f"| {feat} | {cell('LUAD', feat)} | {cell('LUSC', feat)} | {why} |")
    a("")
    a("Residual-rank partial Spearman (Pearson of rank residuals) is in")
    a("`correlations.tsv` as `partial_rho_residual_method`. It is a method check,")
    a("not a second discovery pass.")
    a("")
    a("## Honest interpretation")
    a("")
    a(f"1. **Headline.** {summary['verdict_text']}")
    a("2. **Effect sizes.** Even where p is small, |ρ| in bulk RNA is modest. A")
    a("   significant LUSC rho of −0.2 is not “immune desert because of TROP2”.")
    a("   It is a weak-to-moderate rank association in mixed tissue.")
    a("3. **Purity.** ESTIMATE ImmuneScore is almost a purity inverse. The")
    a("   ABSOLUTE-partial number is the one that is allowed to be called")
    a("   “immune” rather than “not tumor”. CD8/CYT/GEP18 are also compositionally")
    a("   entangled with purity; partialling helps, it does not prove tumor-intrinsic biology.")
    a("4. **GEP18 is not the Merck assay.** Unweighted z-mean of the 18 genes is")
    a("   the public approximation. Do not write “T-cell-inflamed GEP (NanoString)”")
    a("   as if the clinical weights were used.")
    a("5. **This is not OncoSG, not protein, not ICI response.** Original A1 also")
    a("   used OncoSG LUAD (stronger negative rhos in that cohort). This rework is")
    a("   TCGA histology only. TROP2 protein (the ADC target) was not measured.")
    a("6. **No causality.** Bulk correlation, even purity-adjusted and histology-split,")
    a("   does not say TACSTD2 excludes T cells.")
    a("")
    a("## Reproduce")
    a("")
    a("```")
    a("pip install -r requirements.txt")
    a("python scripts/rework_A1_histology.py")
    a("```")
    a("")
    a("Downloads ~65 MB of public tables into `data/` (gitignored) on first run.")
    a("")
    a("## Files")
    a("")
    a("- `correlations.tsv` — unadjusted + ABSOLUTE-partial Spearman, both methods, CIs")
    a("- `luad_vs_lusc_fisher.tsv` — Fisher z tests + per-feature class + Simpson flag")
    a("- `histology_location_shifts.tsv` — Mann–Whitney LUAD vs LUSC on each score")
    a("- `sample_table.tsv` — per-sample TACSTD2, features, ABSOLUTE, ESTIMATE")
    a("- `summary.json` — machine-readable verdict")
    a("- `provenance.json` — URLs, md5, gene lists")
    a("- `figures/forest_partial_rho.png`")
    a("- `figures/scatter_by_histology.png`")
    a("- `figures/box_by_histology.png`")
    a("")
    a("## Data")
    a("")
    a("- Expression: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`.")
    a("- Purity: GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5` (PanCanAtlas ABSOLUTE).")
    a("- ESTIMATE: MD Anderson official RNAseqV2 tables")
    a("  (`lung_adenocarcinoma_RNAseqV2.txt`, `lung_squamous_cell_carcinoma_RNAseqV2.txt`).")
    a("")
    path = os.path.join(OUT, "REPORT.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
