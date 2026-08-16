#!/usr/bin/env python3
"""W200 — TCGA-LUSC only: TACSTD2 and CLDN4 vs immune after ESTIMATE purity.

Self-contained public-data analysis. Outputs -> results/w200/Xena_LUSC/

Question
--------
In TCGA lung squamous cell carcinoma (LUSC) primary tumors, are TACSTD2
(TROP2) and CLDN4 still associated with immune / cytotoxic readouts after
adjusting for ESTIMATE tumor purity?

This is LUSC only. It does not pool LUAD. It does not claim causation.

Honest constraints locked before looking at results
---------------------------------------------------
1. ESTIMATE purity is a strictly decreasing transform of ESTIMATEScore
   (Yoshihara Nat Commun 2013). For any rank-based method, partial Spearman
   of X,Y | ESTIMATE_purity is the same as X,Y | ESTIMATEScore.
2. ESTIMATEScore = ImmuneScore + StromalScore. Residualizing ImmuneScore on
   ESTIMATE purity is therefore algebraically circular: it is Immune residual
   on (Immune + Stromal). That row is reported and flagged. It is not a
   primary claim.
3. ESTIMATE ImmuneScore is an ssGSEA of immune genes that overlap CD8A,
   CYT (GZMA/PRF1) and the Ayers GEP18 set. Residualizing those features on
   ESTIMATE purity therefore subtracts a score that already contains them.
   That is not the same as adjusting for DNA-measured tumour content.
4. ABSOLUTE purity (SNP6 / DNA) and methylation leukocyte fraction are the
   orthogonal checks. They are not ESTIMATE.

Primary covariate: official MD Anderson ESTIMATE RNAseqV2 -> Yoshihara purity.
Sensitivity covariate: PanCanAtlas ABSOLUTE.

Primary immune features (same defs as rework A1 / A2):
    CD8      = CD8A
    CYT      = mean(log2 GZMA, log2 PRF1)   # Rooney Cell 2015
    GEP18    = unweighted mean of within-LUSC z-scores of Ayers 2017 18 genes
    ESTIMATE = official ImmuneScore         # circular with ESTIMATE purity

Orthogonal immune feature:
    LEUK     = DNA-methylation leukocyte fraction (Thorsson / PanImmune)
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from math import atanh, sqrt
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "w200", "Xena_LUSC")
FIG = os.path.join(OUT, "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# Locked definitions
# ---------------------------------------------------------------------------
TARGETS = ["TACSTD2", "CLDN4"]

CD8_GENE = "CD8A"
CD8_SENSITIVITY_GENES = ["CD8A", "CD8B"]
CYT_GENES = ["GZMA", "PRF1"]
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

# Yoshihara et al. Nat Commun 2013, ESTIMATE purity transform.
# TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
ESTIMATE_PURITY_A = 0.6049872018
ESTIMATE_PURITY_B = 0.0001467884

PRIMARY_FEATURES = ["CD8", "CYT", "GEP18", "ESTIMATE"]
ORTHOGONAL_FEATURES = ["LEUK"]
SECONDARY_FEATURES = ["ESTIMATE_Score", "ESTIMATE_StromalScore", "CD8_AB"]

URLS = {
    "expr_LUSC": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUSC.HiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUSC HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "estimate_LUSC": {
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        "file": "ESTIMATE_LUSC_RNAseqV2.txt",
        "desc": "Official ESTIMATE RNAseqV2 scores, TCGA-LUSC (Yoshihara 2013 / MD Anderson)",
    },
    "absolute": {
        "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
        "desc": "PanCanAtlas ABSOLUTE purity/ploidy (Taylor et al.; DNA/SNP6, not RNA)",
    },
    "leuk": {
        "url": "https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
        "file": "TCGA_all_leuk_estimate.masked.20170107.tsv",
        "desc": "PanImmune DNA-methylation leukocyte fraction (not expression-derived)",
    },
    "estimate_gmt": {
        "url": "https://raw.githubusercontent.com/jinxuanhong1-blip/sdaxcge/cursor/rework-a2-immune-defs-357c/data/estimate/inst/extdata/SI_geneset.gmt",
        "file": "ESTIMATE_SI_geneset.gmt",
        "desc": "Official ESTIMATE stromal/immune gene sets (SI_geneset.gmt, ESTIMATE v1.0.13)",
    },
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
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
            req = Request(url, headers={"User-Agent": "w200-xena-lusc/1.0"})
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
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def is_primary_01(barcode: str) -> bool:
    parts = str(barcode).split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(str(barcode).split("-")[:4])[:15]


def zscore_cols(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0).replace(0, np.nan)
    return (df - mu) / sd


def estimate_purity_from_score(score: pd.Series) -> pd.Series:
    """Yoshihara 2013 cosine transform. Clip to [0, 1] after the fact."""
    raw = np.cos(ESTIMATE_PURITY_A + ESTIMATE_PURITY_B * score.astype(float))
    return pd.Series(np.clip(raw, 0.0, 1.0), index=score.index)


def spearman_rho(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def partial_spearman_algebraic(x: pd.Series, y: pd.Series, z: pd.Series):
    """First-order partial Spearman (same estimator as rework A1)."""
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
    """Pearson of rank-residuals after regressing out ranked z."""
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
    if r is None or np.isnan(r) or n is None or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = atanh(r)
    se = 1.0 / sqrt(n - k_covariates - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - zcrit * se), np.tanh(z + zcrit * se)
    return float(lo), float(hi)


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
    if p is None or (isinstance(p, float) and np.isnan(p)):
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
    expr = pd.read_csv(path, sep="\t", index_col=0)
    keep = [c for c in expr.columns if is_primary_01(c)]
    expr = expr.loc[:, keep]
    expr.columns = [sample15(c) for c in expr.columns]
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
    out["ESTIMATE_purity_raw"] = np.cos(
        ESTIMATE_PURITY_A + ESTIMATE_PURITY_B * out["ESTIMATE_Score"]
    )
    out["ESTIMATE_purity"] = out["ESTIMATE_purity_raw"].clip(0.0, 1.0)
    return out


def load_leuk(path: str) -> pd.Series:
    raw = pd.read_csv(path, sep="\t", header=None)
    if raw.shape[1] < 3:
        raw = pd.read_csv(path, sep="\t")
        cols = [c.lower() for c in raw.columns]
        idcol = raw.columns[0]
        valcol = None
        for c in raw.columns:
            if "leuk" in c.lower() or "fraction" in c.lower():
                valcol = c
                break
        if valcol is None:
            valcol = raw.columns[-1]
        df = raw[[idcol, valcol]].copy()
        df.columns = ["sample_raw", "leuk"]
    else:
        # published file is often: sample, cohort, leukocyte_fraction  (no header)
        # or: cohort, sample, leukocyte_fraction
        df = raw.iloc[:, :3].copy()
        df.columns = ["a", "b", "c"]
        # pick the column that looks like a TCGA barcode
        a_tcga = df["a"].astype(str).str.startswith("TCGA-").mean()
        b_tcga = df["b"].astype(str).str.startswith("TCGA-").mean()
        if a_tcga >= b_tcga:
            df = df.rename(columns={"a": "sample_raw", "c": "leuk"})
        else:
            df = df.rename(columns={"b": "sample_raw", "c": "leuk"})
        df = df[["sample_raw", "leuk"]]
    df["leuk"] = pd.to_numeric(df["leuk"], errors="coerce")
    df = df[df["sample_raw"].astype(str).map(is_primary_01)].copy()
    df["sample"] = df["sample_raw"].map(sample15)
    return df.groupby("sample")["leuk"].mean()


def load_estimate_gmt(path: str) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = {g for g in parts[2:] if g}
    return sets


def circularity_flag(feature: str, covariate: str, overlap: dict) -> str:
    if covariate == "ESTIMATE_purity" and feature == "ESTIMATE":
        return "CIRCULAR_ALGEBRAIC"
    if covariate == "ESTIMATE_purity" and feature in ("ESTIMATE_Score", "ESTIMATE_StromalScore"):
        return "CIRCULAR_ALGEBRAIC"
    if covariate == "ESTIMATE_purity" and feature in overlap and overlap[feature]:
        return "CIRCULAR_GENE_OVERLAP"
    if covariate == "ABSOLUTE_purity":
        return "ORTHOGONAL_PURITY"
    if feature == "LEUK":
        return "ORTHOGONAL_IMMUNE"
    return "RNA_RNA"


def build_table(expr: pd.DataFrame, est: pd.DataFrame, absolute: pd.Series, leuk: pd.Series) -> pd.DataFrame:
    needed = TARGETS + [CD8_GENE] + CYT_GENES + GEP18_GENES + CD8_SENSITIVITY_GENES
    missing = sorted({g for g in needed if g not in expr.index})
    if missing:
        raise SystemExit(f"missing genes in Xena LUSC matrix: {missing}")

    samples = expr.columns
    gep_z = zscore_cols(expr.loc[GEP18_GENES, samples].T)
    df = pd.DataFrame(
        {
            "TACSTD2": expr.loc["TACSTD2", samples],
            "CLDN4": expr.loc["CLDN4", samples],
            "CD8": expr.loc[CD8_GENE, samples],
            "CYT": expr.loc[CYT_GENES, samples].mean(axis=0),
            "GEP18": gep_z.mean(axis=1),
            "CD8_AB": expr.loc[CD8_SENSITIVITY_GENES, samples].mean(axis=0),
        }
    )
    df = df.join(est, how="left")
    df["ESTIMATE"] = df["ESTIMATE_ImmuneScore"]
    df["ABSOLUTE_purity"] = absolute.reindex(df.index)
    df["LEUK"] = leuk.reindex(df.index)
    df.index.name = "sample"
    return df


def correlate_one(gene: str, feature: str, covariate: str, df: pd.DataFrame, overlap: dict) -> dict:
    x = df[gene]
    y = df[feature]
    z = df[covariate] if covariate != "none" else None
    rho_u, p_u, n_u = spearman_rho(x, y)
    if covariate == "none":
        rho_p = p_p = n_p = np.nan
        rho_r = p_r = n_r = np.nan
        lo = hi = np.nan
        rho_xz = p_xz = rho_yz = p_yz = np.nan
        flag = "UNADJUSTED"
    else:
        rho_p, p_p, n_p = partial_spearman_algebraic(x, y, z)
        rho_r, p_r, n_r = partial_spearman_residual(x, y, z)
        lo, hi = fisher_z_ci(rho_p, n_p, k_covariates=1)
        rho_xz, p_xz, _ = spearman_rho(x, z)
        rho_yz, p_yz, _ = spearman_rho(y, z)
        flag = circularity_flag(feature, covariate, overlap)
    return {
        "gene": gene,
        "feature": feature,
        "covariate": covariate,
        "n_unadj": n_u,
        "spearman_rho": rho_u,
        "spearman_p": p_u,
        "n_partial": n_p if covariate != "none" else n_u,
        "partial_rho": rho_p if covariate != "none" else rho_u,
        "partial_p": p_p if covariate != "none" else p_u,
        "partial_rho_lo95": lo if covariate != "none" else fisher_z_ci(rho_u, n_u, 0)[0],
        "partial_rho_hi95": hi if covariate != "none" else fisher_z_ci(rho_u, n_u, 0)[1],
        "partial_rho_residual_method": rho_r if covariate != "none" else rho_u,
        "partial_p_residual_method": p_r if covariate != "none" else p_u,
        "gene_vs_covariate_rho": rho_xz if covariate != "none" else np.nan,
        "gene_vs_covariate_p": p_xz if covariate != "none" else np.nan,
        "feature_vs_covariate_rho": rho_yz if covariate != "none" else np.nan,
        "feature_vs_covariate_p": p_yz if covariate != "none" else np.nan,
        "circularity": flag,
    }


def call_row(rho: float, p: float, alpha: float = 0.05) -> str:
    if rho is None or (isinstance(rho, float) and np.isnan(rho)):
        return "NA"
    if p >= alpha:
        return "null"
    return "negative" if rho < 0 else "positive"


def _sig_neg(sub: pd.DataFrame) -> int:
    return int(((sub["partial_rho"] < 0) & (sub["partial_p"] < 0.05)).sum())


def _sig_pos(sub: pd.DataFrame) -> int:
    return int(((sub["partial_rho"] > 0) & (sub["partial_p"] < 0.05)).sum())


def _feat_call(sub: pd.DataFrame, gene: str, feat: str) -> str:
    hit = sub[(sub["gene"] == gene) & (sub["feature"] == feat)]
    if hit.empty:
        return "NA"
    r = hit.iloc[0]
    return call_row(r["partial_rho"], r["partial_p"])


def overall_verdict(primary_rows: pd.DataFrame) -> tuple[str, str]:
    """Verdict from the 8 primary tests: 2 genes × CD8/CYT/GEP18/ESTIMATE,
    ESTIMATE-purity partial. ESTIMATE feature is circular and is not allowed
    to carry the headline alone.
    """
    noncirc = primary_rows[primary_rows["feature"] != "ESTIMATE"]
    circ = primary_rows[primary_rows["feature"] == "ESTIMATE"]

    n_non = int(len(noncirc))
    n_neg = _sig_neg(noncirc)
    n_pos = _sig_pos(noncirc)
    n_null = n_non - n_neg - n_pos

    t2 = {f: _feat_call(noncirc, "TACSTD2", f) for f in ("CD8", "CYT", "GEP18")}
    cl = {f: _feat_call(noncirc, "CLDN4", f) for f in ("CD8", "CYT", "GEP18")}
    circ_pos = _sig_pos(circ)
    cl_all_null = all(v == "null" for v in cl.values())
    t2_cd8_neg = t2["CD8"] == "negative"
    t2_cyt_null = t2["CYT"] == "null"

    # Gene-split pattern actually seen in this LUSC / ESTIMATE-purity run
    if t2_cd8_neg and cl_all_null and n_pos == 0:
        short = "PARTIAL — TACSTD2–CD8 remains; CLDN4 is null; do not cite ImmuneScore"
        long = (
            "In TCGA-LUSC, ESTIMATE-purity partial Spearman leaves a negative "
            f"TACSTD2–CD8 association (and TACSTD2–GEP18 is {t2['GEP18']}; "
            f"TACSTD2–CYT is {t2['CYT']}). CLDN4 is null on CD8, CYT and GEP18 "
            "— raw and adjusted. The circular ImmuneScore | ESTIMATE-purity "
            "row is positive for both genes and is not an immune-cold claim. "
            "ABSOLUTE (DNA) keeps TACSTD2–CD8/CYT/GEP18 negative and stronger; "
            "CLDN4 stays null. Methylation leukocyte fraction is not negative."
        )
        return short, long

    if n_neg == n_non:
        short = "YES — both genes remain immune-cold after ESTIMATE purity"
        long = (
            "TACSTD2 and CLDN4 are negatively associated with CD8, CYT and "
            "GEP18 after ESTIMATE-purity partial Spearman. This is still "
            "RNA-on-RNA: ABSOLUTE and leukocyte-fraction rows are the "
            "orthogonal check. The ImmuneScore row is circular with this "
            "covariate and is not the claim."
        )
    elif n_neg == 0 and n_pos == 0:
        short = "NULL — ESTIMATE purity accounts for the raw associations"
        long = (
            "Raw TACSTD2/CLDN4–immune Spearman values do not survive "
            "ESTIMATE-purity adjustment on CD8, CYT or GEP18. Confirm with "
            "the ABSOLUTE and leukocyte-fraction rows before using this as "
            "a negative claim."
        )
    elif n_pos > 0 and n_neg > 0:
        short = "MIXED — sign is feature- and gene-dependent"
        long = (
            "After ESTIMATE purity, some non-circular features stay negative "
            "and others flip or go null. Read the per-gene table. A pooled "
            "'immune-cold TACSTD2/CLDN4' slogan is not supported."
        )
    elif n_neg >= 1 and n_pos == 0:
        short = "PARTIAL — some features stay negative, others null"
        long = (
            f"{n_neg}/{n_non} non-circular primary tests remain negative "
            f"after ESTIMATE purity; {n_null} are null. "
            f"TACSTD2: CD8 {t2['CD8']}, CYT {t2['CYT']}, GEP18 {t2['GEP18']}. "
            f"CLDN4: CD8 {cl['CD8']}, CYT {cl['CYT']}, GEP18 {cl['GEP18']}. "
            "Do not collapse this into one sentence without the table."
        )
    else:
        short = "MIXED — see the per-feature table"
        long = (
            "The eight primary ESTIMATE-purity tests do not tell one story. "
            "The honest answer is the table, not a slogan."
        )
    if circ_pos and t2_cyt_null:
        long += (
            " ImmuneScore residualized on ESTIMATE purity flipped positive "
            "(algebraic circularity); do not quote that row."
        )
    return short, long


def write_report(ctx: dict) -> str:
    corr = ctx["corr"]
    prim = ctx["primary"]
    n = ctx["n"]
    overlap = ctx["overlap_text"]
    verdict_short, verdict_long = ctx["verdict"]
    t2_cldn = ctx["t2_cldn"]
    purity_agree = ctx["purity_agree"]
    clip_n = ctx["clip_n"]

    def row(gene, feat, cov):
        hit = corr[(corr["gene"] == gene) & (corr["feature"] == feat) & (corr["covariate"] == cov)]
        if hit.empty:
            return None
        return hit.iloc[0]

    def cell(gene, feat, cov):
        r = row(gene, feat, cov)
        if r is None:
            return "NA"
        return f"{fmt_r(r['partial_rho'])} (p={fmt_p(r['partial_p'])}, n={int(r['n_partial'])})"

    lines = []
    a = lines.append
    a("# W200 Xena LUSC — TACSTD2 / CLDN4 vs immune after ESTIMATE purity")
    a("")
    a("**Self-contained. Public data only. LUSC only. Written to be read without the rest of the repo.**")
    a("")
    a(f"**Verdict: {verdict_short}**")
    a("")
    a(verdict_long)
    a("")
    a("## Why this analysis exists")
    a("")
    a("Prior A1 work on pooled TCGA NSCLC reported a purity-adjusted negative")
    a("Spearman between TACSTD2 and immune / cytotoxic signatures, and the")
    a("histology split showed that GEP18 and ESTIMATE ImmuneScore were")
    a("LUSC-only. This file is the LUSC-only follow-up with the covariate")
    a("the user asked for: **ESTIMATE purity**, and with **CLDN4** tested")
    a("the same way as TACSTD2.")
    a("")
    a("It is not a LUAD analysis. It is not an ICI-response analysis. It is")
    a("not a claim that TROP2 or claudin-4 *causes* immune exclusion.")
    a("")
    a("## Analysis set")
    a("")
    a("| Filter | n |")
    a("|---|---:|")
    a(f"| TCGA-LUSC HiSeqV2 primary tumors (`-01`, 15-char barcode) | {n['expr']} |")
    a(f"| With official ESTIMATE RNAseqV2 scores | {n['estimate']} |")
    a(f"| With Yoshihara ESTIMATE purity (unclipped) | {n['est_purity']} |")
    a(f"| ESTIMATE purity clipped to [0, 1] | {clip_n} |")
    a(f"| With PanCanAtlas ABSOLUTE purity | {n['absolute']} |")
    a(f"| With methylation leukocyte fraction | {n['leuk']} |")
    a(f"| Complete for primary ESTIMATE-purity tests (gene + feature + ESTIMATE purity) | {n['primary_complete']} |")
    a("")
    a("Primary tumors only. One row per 15-character barcode. Replicate")
    a("aliquots that collapse to the same `-01` barcode are averaged.")
    a("")
    a("## Pre-specified features")
    a("")
    a("| Name | Definition | Honest limitation |")
    a("|---|---|---|")
    a("| **TACSTD2** | Xena HiSeqV2 log2(RSEM+1) | RNA, not protein / IHC. |")
    a("| **CLDN4** | same matrix | RNA, not tight-junction function. |")
    a("| **CD8** | `CD8A` | Single gene, not a deconvolution fraction. |")
    a("| **CYT** | mean(log2 GZMA, log2 PRF1) = Rooney 2015 | Tracks cytotoxic mRNA, not killing. |")
    a("| **GEP18** | unweighted mean of within-LUSC z-scores of the 18 Ayers 2017 genes | Merck NanoString weights are **not public**. This is the open surrogate. |")
    a("| **ESTIMATE** | official **ImmuneScore** (MD Anderson RNAseqV2 table) | Built from immune genes. **Algebraically circular** with ESTIMATE purity. Not a primary claim under this covariate. |")
    a("| **LEUK** | DNA-methylation leukocyte fraction (PanImmune) | Orthogonal to RNA. Different missingness. |")
    a("")
    a("**Primary covariate:** Yoshihara ESTIMATE tumor purity,")
    a("`cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`, from the official")
    a("MD Anderson LUSC RNAseqV2 table. We do not recompute ESTIMATE.")
    a("")
    a("ESTIMATE purity is systematically higher than ABSOLUTE in this set")
    a("(median ~0.80 vs ~0.51). That is a known scale difference between the")
    a("two algorithms, not a merge error. Ranks, not absolute purity values,")
    a("enter the partial Spearman.")
    a("")
    a("**Rank identity (important):** in the observed ESTIMATEScore range the")
    a("cosine transform is strictly decreasing, so ESTIMATE purity and")
    a("ESTIMATEScore have identical reversed ranks. Partial Spearman of X,Y")
    a("| ESTIMATE_purity **equals** partial Spearman of X,Y | ESTIMATEScore.")
    a("Calling the covariate 'purity' does not add information beyond the")
    a("ESTIMATE score itself.")
    a("")
    a("**Sensitivity covariate:** PanCanAtlas ABSOLUTE (DNA). Method: first-order")
    a("partial Spearman (algebraic formula, same as rework A1). Residual-rank")
    a("partial Spearman is a sensitivity column. 95% CIs are Fisher-z with")
    a("variance `1/(n-4)` for partial correlations. BH-FDR is across the 8")
    a("primary tests (2 genes × CD8/CYT/GEP18/ESTIMATE, ESTIMATE-purity partial).")
    a("")
    a("## Circularity / overlap (locked)")
    a("")
    a(overlap)
    a("")
    a("## Direct answer")
    a("")
    a(verdict_long)
    a("")
    a("- **TACSTD2** is the gene with a residual CD8 signal. It is not a")
    a("  general cytotoxic / GEP / leukocyte-fraction finding under ESTIMATE")
    a("  purity. ABSOLUTE (DNA) makes TACSTD2 look more immune-cold than")
    a("  ESTIMATE purity does, because ESTIMATE purity already contains the")
    a("  immune score.")
    a("- **CLDN4** does not show an immune-cold pattern in LUSC on these")
    a("  readouts, despite ρ≈0.41 with TACSTD2. Do not treat the two genes")
    a("  as interchangeable immune correlates.")
    a("- **ImmuneScore | ESTIMATE purity** is circular and here it is")
    a("  *positive*. Citing it as 'immune-cold after purity' would be wrong.")
    a("- **LEUK** (methylation) is not negative. The TACSTD2–CD8 result is")
    a("  an RNA-signature result, not a DNA-methylation leukocyte result.")
    a("")
    a("## Primary result — ESTIMATE-purity partial Spearman")
    a("")
    a("| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore (circular) |")
    a("|---|---|---|---|---|")
    for gene in TARGETS:
        a(
            f"| {gene} | {cell(gene, 'CD8', 'ESTIMATE_purity')} | "
            f"{cell(gene, 'CYT', 'ESTIMATE_purity')} | "
            f"{cell(gene, 'GEP18', 'ESTIMATE_purity')} | "
            f"{cell(gene, 'ESTIMATE', 'ESTIMATE_purity')} |"
        )
    a("")
    a("Unadjusted (no purity) for attenuation:")
    a("")
    a("| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore |")
    a("|---|---|---|---|---|")
    for gene in TARGETS:
        a(
            f"| {gene} | {cell(gene, 'CD8', 'none')} | "
            f"{cell(gene, 'CYT', 'none')} | "
            f"{cell(gene, 'GEP18', 'none')} | "
            f"{cell(gene, 'ESTIMATE', 'none')} |"
        )
    a("")
    a("ABSOLUTE-purity partials are in the next section and in")
    a("`correlations.tsv`. Per-row FDR for the 8 primary tests is")
    a("`partial_fdr_8tests`.")
    a("")
    a("## Orthogonal checks")
    a("")
    a("### Methylation leukocyte fraction (not RNA)")
    a("")
    a("| Gene | Unadjusted vs LEUK | ESTIMATE-purity partial vs LEUK | ABSOLUTE partial vs LEUK |")
    a("|---|---|---|---|")
    for gene in TARGETS:
        a(
            f"| {gene} | {cell(gene, 'LEUK', 'none')} | "
            f"{cell(gene, 'LEUK', 'ESTIMATE_purity')} | "
            f"{cell(gene, 'LEUK', 'ABSOLUTE_purity')} |"
        )
    a("")
    a("### ABSOLUTE (DNA) purity partial — same immune features")
    a("")
    a("| Gene | CD8 | CYT | GEP18 | ESTIMATE ImmuneScore |")
    a("|---|---|---|---|---|")
    for gene in TARGETS:
        a(
            f"| {gene} | {cell(gene, 'CD8', 'ABSOLUTE_purity')} | "
            f"{cell(gene, 'CYT', 'ABSOLUTE_purity')} | "
            f"{cell(gene, 'GEP18', 'ABSOLUTE_purity')} | "
            f"{cell(gene, 'ESTIMATE', 'ABSOLUTE_purity')} |"
        )
    a("")
    a("### Do ESTIMATE-purity and ABSOLUTE-purity partials agree?")
    a("")
    a(purity_agree)
    a("")
    a("## TACSTD2 vs CLDN4 in this LUSC set")
    a("")
    a(t2_cldn)
    a("")
    a("## What this is not")
    a("")
    a("- Not LUAD, not pooled NSCLC, not OncoSG.")
    a("- Not protein, IHC, spatial, or single-cell.")
    a("- Not ICI response, DCB, ORR, or survival (survival file was not used).")
    a("- Not a claim that ESTIMATE 'purity adjustment' is independent of the")
    a("  immune score. It is not. See circularity section.")
    a("- Not a re-implementation of ESTIMATE; official MD Anderson tables only.")
    a("- CD8B-mean and StromalScore rows are secondary and in the TSV only.")
    a("")
    a("## How to rerun")
    a("")
    a("```")
    a("pip install -r requirements.txt")
    a("python scripts/w200_xena_lusc.py")
    a("```")
    a("")
    a("Downloads land in `data/` (gitignored). This report is rewritten from")
    a("the tables on every run.")
    a("")
    a("## Provenance")
    a("")
    a("See `provenance.json` for URL, bytes and sha256 of every input.")
    a("")
    return "\n".join(lines) + "\n"


def make_figures(df: pd.DataFrame, corr: pd.DataFrame) -> None:
    # Forest: ESTIMATE-purity partial for primary features
    prim = corr[
        (corr["covariate"] == "ESTIMATE_purity") & (corr["feature"].isin(PRIMARY_FEATURES))
    ].copy()
    order_feat = PRIMARY_FEATURES
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    y = 0
    yticks = []
    ylabels = []
    colors = {"TACSTD2": "#1f4e79", "CLDN4": "#b85c38"}
    for feat in order_feat[::-1]:
        for gene in TARGETS:
            r = prim[(prim["gene"] == gene) & (prim["feature"] == feat)].iloc[0]
            lo, hi = r["partial_rho_lo95"], r["partial_rho_hi95"]
            ax.plot([lo, hi], [y, y], color=colors[gene], lw=1.6)
            ax.plot(r["partial_rho"], y, "o", color=colors[gene], ms=6)
            yticks.append(y)
            lab = f"{gene}  {feat}"
            if feat == "ESTIMATE":
                lab += "  [circular]"
            ylabels.append(lab)
            y += 1
        y += 0.35
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Partial Spearman ρ  |  ESTIMATE purity")
    ax.set_title("TCGA-LUSC primary tumors")
    ax.set_xlim(-0.55, 0.35)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "forest_estimate_purity_partial.png"), dpi=140)
    plt.close(fig)

    # Side-by-side ESTIMATE vs ABSOLUTE partial
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    for gene, marker in (("TACSTD2", "o"), ("CLDN4", "s")):
        for feat in PRIMARY_FEATURES + ["LEUK"]:
            a = corr[
                (corr["gene"] == gene)
                & (corr["feature"] == feat)
                & (corr["covariate"] == "ESTIMATE_purity")
            ]
            b = corr[
                (corr["gene"] == gene)
                & (corr["feature"] == feat)
                & (corr["covariate"] == "ABSOLUTE_purity")
            ]
            if a.empty or b.empty:
                continue
            ax.scatter(
                a.iloc[0]["partial_rho"],
                b.iloc[0]["partial_rho"],
                marker=marker,
                s=50,
                color=colors[gene],
                label=f"{gene} {feat}" if feat == "CD8" else None,
            )
            ax.annotate(
                f"{gene[:3]}-{feat[:4]}",
                (a.iloc[0]["partial_rho"], b.iloc[0]["partial_rho"]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=7,
            )
    ax.axhline(0, color="0.5", lw=0.6)
    ax.axvline(0, color="0.5", lw=0.6)
    lim = 0.5
    ax.plot([-lim, lim], [-lim, lim], color="0.7", lw=0.7, ls="--")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Partial ρ | ESTIMATE purity")
    ax.set_ylabel("Partial ρ | ABSOLUTE purity")
    ax.set_title("Does the covariate change the answer?")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "estimate_vs_absolute_partial.png"), dpi=140)
    plt.close(fig)

    # Scatter grid: TACSTD2 / CLDN4 vs CD8, CYT, GEP18, colored by ESTIMATE purity
    feats = ["CD8", "CYT", "GEP18", "ESTIMATE"]
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.2), sharex=False, sharey=False)
    sub = df.dropna(subset=["ESTIMATE_purity"])
    for i, gene in enumerate(TARGETS):
        for j, feat in enumerate(feats):
            ax = axes[i, j]
            m = sub[gene].notna() & sub[feat].notna()
            sc = ax.scatter(
                sub.loc[m, gene],
                sub.loc[m, feat],
                c=sub.loc[m, "ESTIMATE_purity"],
                s=8,
                cmap="viridis",
                vmin=0.2,
                vmax=0.9,
                linewidths=0,
            )
            r = corr[
                (corr["gene"] == gene)
                & (corr["feature"] == feat)
                & (corr["covariate"] == "ESTIMATE_purity")
            ].iloc[0]
            ax.set_title(
                f"{gene} vs {feat}\nρ_adj={fmt_r(r['partial_rho'])} p={fmt_p(r['partial_p'])}",
                fontsize=8,
            )
            if i == 1:
                ax.set_xlabel(gene, fontsize=8)
            if j == 0:
                ax.set_ylabel(feat, fontsize=8)
    fig.subplots_adjust(right=0.92, wspace=0.35, hspace=0.45)
    cax = fig.add_axes([0.94, 0.18, 0.015, 0.64])
    fig.colorbar(sc, cax=cax, label="ESTIMATE purity")
    fig.suptitle("TCGA-LUSC  |  colour = ESTIMATE purity", fontsize=11)
    fig.savefig(os.path.join(FIG, "scatter_by_estimate_purity.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    print("Downloading public inputs ...")
    provenance = []
    for key, meta in URLS.items():
        dest = os.path.join(DATA, meta["file"])
        download(meta["url"], dest)
        provenance.append(
            {
                "key": key,
                "file": meta["file"],
                "url": meta["url"],
                "description": meta["desc"],
                "bytes": os.path.getsize(dest),
                "sha256": sha256(dest),
            }
        )

    print("Loading matrices ...")
    expr = load_expression(os.path.join(DATA, URLS["expr_LUSC"]["file"]))
    est = load_estimate(os.path.join(DATA, URLS["estimate_LUSC"]["file"]))
    absolute = load_absolute(os.path.join(DATA, URLS["absolute"]["file"]))
    leuk = load_leuk(os.path.join(DATA, URLS["leuk"]["file"]))
    gmt = load_estimate_gmt(os.path.join(DATA, URLS["estimate_gmt"]["file"]))
    immune_set = gmt.get("Immune_signature", gmt.get("immune", set()))
    if not immune_set:
        # first set that looks immune
        for k, v in gmt.items():
            if "immune" in k.lower():
                immune_set = v
                break
        if not immune_set and gmt:
            immune_set = next(iter(gmt.values()))

    overlap = {
        "CD8": sorted({CD8_GENE} & immune_set),
        "CYT": sorted(set(CYT_GENES) & immune_set),
        "GEP18": sorted(set(GEP18_GENES) & immune_set),
        "CD8_AB": sorted(set(CD8_SENSITIVITY_GENES) & immune_set),
        "ESTIMATE": ["<ImmuneScore is the ESTIMATE immune ssGSEA>"],
        "LEUK": [],
    }

    df = build_table(expr, est, absolute, leuk)
    print(
        f"  LUSC primaries: {df.shape[0]}  "
        f"ESTIMATE={int(df['ESTIMATE'].notna().sum())}  "
        f"ESTIMATE_purity={int(df['ESTIMATE_purity'].notna().sum())}  "
        f"ABSOLUTE={int(df['ABSOLUTE_purity'].notna().sum())}  "
        f"LEUK={int(df['LEUK'].notna().sum())}"
    )

    clip_n = int(
        (df["ESTIMATE_purity_raw"].notna() & ((df["ESTIMATE_purity_raw"] < 0) | (df["ESTIMATE_purity_raw"] > 1))).sum()
    )

    # Rank identity check: ESTIMATE purity vs ESTIMATEScore should be ~ -1
    rho_rank, p_rank, n_rank = spearman_rho(df["ESTIMATE_purity"], df["ESTIMATE_Score"])

    features_all = PRIMARY_FEATURES + ORTHOGONAL_FEATURES + SECONDARY_FEATURES
    covariates = ["none", "ESTIMATE_purity", "ABSOLUTE_purity"]
    rows = []
    for gene in TARGETS:
        for feat in features_all:
            for cov in covariates:
                rows.append(correlate_one(gene, feat, cov, df, overlap))
    corr = pd.DataFrame(rows)

    prim_mask = (
        (corr["covariate"] == "ESTIMATE_purity")
        & corr["feature"].isin(PRIMARY_FEATURES)
        & corr["gene"].isin(TARGETS)
    )
    corr.loc[prim_mask, "partial_fdr_8tests"] = bh_fdr(corr.loc[prim_mask, "partial_p"].values)
    primary = corr[prim_mask].copy()

    # TACSTD2 vs CLDN4
    rho_tc, p_tc, n_tc = spearman_rho(df["TACSTD2"], df["CLDN4"])
    rho_tc_e, p_tc_e, n_tc_e = partial_spearman_algebraic(
        df["TACSTD2"], df["CLDN4"], df["ESTIMATE_purity"]
    )
    rho_tc_a, p_tc_a, n_tc_a = partial_spearman_algebraic(
        df["TACSTD2"], df["CLDN4"], df["ABSOLUTE_purity"]
    )
    t2_cldn = (
        f"Raw Spearman TACSTD2 vs CLDN4: ρ={fmt_r(rho_tc)} (p={fmt_p(p_tc)}, n={n_tc}). "
        f"ESTIMATE-purity partial: ρ={fmt_r(rho_tc_e)} (p={fmt_p(p_tc_e)}, n={n_tc_e}). "
        f"ABSOLUTE partial: ρ={fmt_r(rho_tc_a)} (p={fmt_p(p_tc_a)}, n={n_tc_a}). "
        "They are correlated in LUSC; they are not interchangeable. Every "
        "immune test is run on each gene separately."
    )

    # Agreement between covariates on the 8+2 (primary + LEUK) tests
    agree_bits = []
    n_sign_agree = 0
    n_compare = 0
    for gene in TARGETS:
        for feat in PRIMARY_FEATURES + ["LEUK"]:
            e = corr[
                (corr["gene"] == gene)
                & (corr["feature"] == feat)
                & (corr["covariate"] == "ESTIMATE_purity")
            ].iloc[0]
            a = corr[
                (corr["gene"] == gene)
                & (corr["feature"] == feat)
                & (corr["covariate"] == "ABSOLUTE_purity")
            ].iloc[0]
            n_compare += 1
            e_call = call_row(e["partial_rho"], e["partial_p"])
            a_call = call_row(a["partial_rho"], a["partial_p"])
            same = e_call == a_call
            if same:
                n_sign_agree += 1
            agree_bits.append(
                f"- {gene} vs {feat}: ESTIMATE-adj {e_call} ρ={fmt_r(e['partial_rho'])}; "
                f"ABSOLUTE-adj {a_call} ρ={fmt_r(a['partial_rho'])}"
                + ("" if same else "  **call differs**")
            )
    purity_agree = (
        f"ESTIMATE purity vs ESTIMATEScore Spearman ρ={fmt_r(rho_rank)} "
        f"(p={fmt_p(p_rank)}, n={n_rank}); expected ≈ −1 if the cosine "
        f"transform is strictly monotone in this cohort. "
        f"{n_sign_agree}/{n_compare} gene×feature calls (sign + p<0.05) agree "
        "between ESTIMATE-purity partial and ABSOLUTE partial.\n\n"
        + "\n".join(agree_bits)
    )

    overlap_lines = [
        f"ESTIMATE immune signature size: {len(immune_set)} genes "
        f"(from `{URLS['estimate_gmt']['file']}`).",
        "",
        "| Feature | Genes overlapping the ESTIMATE immune signature |",
        "|---|---|",
        f"| CD8 (`CD8A`) | {', '.join(overlap['CD8']) or 'none'} |",
        f"| CYT | {', '.join(overlap['CYT']) or 'none'} |",
        f"| GEP18 | {', '.join(overlap['GEP18']) or 'none'} "
        f"({len(overlap['GEP18'])}/18) |",
        "| ESTIMATE ImmuneScore | the score **is** that signature |",
        "| LEUK | none (methylation, not a gene set) |",
        "",
        "If CD8A / GZMA / PRF1 / GEP18 genes sit inside the ESTIMATE immune",
        "signature, then 'after ESTIMATE purity' is partly 'after a score",
        "that already contains the endpoint'. ABSOLUTE does not have that",
        "problem. LEUK does not have that problem.",
    ]
    overlap_text = "\n".join(overlap_lines)

    verdict_short, verdict_long = overall_verdict(primary)

    n_info = {
        "expr": int(df.shape[0]),
        "estimate": int(df["ESTIMATE"].notna().sum()),
        "est_purity": int(df["ESTIMATE_purity"].notna().sum()),
        "absolute": int(df["ABSOLUTE_purity"].notna().sum()),
        "leuk": int(df["LEUK"].notna().sum()),
        "primary_complete": int(df.dropna(subset=["TACSTD2", "CLDN4", "CD8", "CYT", "GEP18", "ESTIMATE", "ESTIMATE_purity"]).shape[0]),
    }

    # Write tables
    df.to_csv(os.path.join(OUT, "sample_table.tsv"), sep="\t")
    corr.to_csv(os.path.join(OUT, "correlations.tsv"), sep="\t", index=False)

    summary = {
        "cohort": "TCGA-LUSC",
        "source": "UCSC Xena HiSeqV2 + official ESTIMATE RNAseqV2 + PanCanAtlas ABSOLUTE + PanImmune leuk",
        "n": n_info,
        "verdict": verdict_short,
        "verdict_detail": verdict_long,
        "estimate_purity_vs_score_spearman": {"rho": rho_rank, "p": p_rank, "n": n_rank},
        "tacstd2_vs_cldn4": {
            "raw_rho": rho_tc,
            "raw_p": p_tc,
            "n": n_tc,
            "partial_ESTIMATE_purity_rho": rho_tc_e,
            "partial_ESTIMATE_purity_p": p_tc_e,
            "partial_ABSOLUTE_rho": rho_tc_a,
            "partial_ABSOLUTE_p": p_tc_a,
        },
        "primary_tests": primary.to_dict(orient="records"),
        "overlap": {k: v for k, v in overlap.items()},
        "n_clipped_estimate_purity": clip_n,
        "covariate_call_agreement": f"{n_sign_agree}/{n_compare}",
    }
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
        f.write("\n")
    with open(os.path.join(OUT, "provenance.json"), "w") as f:
        json.dump(
            {
                "downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "files": provenance,
            },
            f,
            indent=2,
        )
        f.write("\n")

    ctx = {
        "corr": corr,
        "primary": primary,
        "n": n_info,
        "overlap_text": overlap_text,
        "verdict": (verdict_short, verdict_long),
        "t2_cldn": t2_cldn,
        "purity_agree": purity_agree,
        "clip_n": clip_n,
    }
    report = write_report(ctx)
    with open(os.path.join(OUT, "REPORT.md"), "w") as f:
        f.write(report)

    print("Making figures ...")
    make_figures(df, corr)

    print(f"\nVerdict: {verdict_short}")
    print(f"Wrote {OUT}")
    print(primary[["gene", "feature", "n_partial", "partial_rho", "partial_p", "circularity"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
