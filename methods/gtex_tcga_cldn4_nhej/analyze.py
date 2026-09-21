#!/usr/bin/env python3
"""GTEx + TCGA pan-cancer: CLDN4 vs PRKDC, LIG4, STING1, HLA-A/B/C.

Public bulk RNA only. Joint matrix is the UCSC Xena Toil recomputation
(TcgaTargetGtex_rsem_gene_tpm): log2(RSEM TPM + 0.001), GENCODE v23,
same pipeline for TCGA and GTEx (Vivian et al., Nat Biotechnol 2017).

STING1 is TMEM173 on this GENCODE v23 freeze. Tables use the current
symbol STING1 and record the matrix symbol.

Primary family (pre-specified): Spearman CLDN4 vs each of the six genes
inside three lung cohorts, BH together (18 tests):
  TCGA-LUAD primary tumor, TCGA-LUSC primary tumor, GTEx lung.

Pan-cancer TCGA primaries and other GTEx tissues are consistency context,
with their own BH families. They are not pooled with each other or with
the lung tests into one correlation.

Sensitivity, not the primary estimand:
  * rank-residual partial Spearman on PTPRC (leukocyte) and on EPCAM
    (epithelial). Both covariates are same-matrix RNA, not DNA purity.
  * TCGA-LUAD / TCGA-LUSC legacy HiSeqV2 (log2 norm_count+1), the freeze
    used by earlier lung HLA tables. Separate 12-test BH. Not mixed into
    the Toil 18.

Lung adjacent normal (TCGA solid tissue normal) is reported separately.
It is not part of the 18.

Unit of analysis: one donor per cohort (TCGA patient; GTEx donor).
Replicate aliquots are averaged. No imputation.

Outputs -> methods/gtex_tcga_cldn4_nhej/
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
import warnings
from math import atanh, sqrt
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
TABLES = os.path.join(HERE, "tables")
FIG = os.path.join(HERE, "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# Canonical endpoint order. Companions are not in any FDR family.
ENDPOINTS = ["PRKDC", "LIG4", "STING1", "HLA-A", "HLA-B", "HLA-C"]
COMPANIONS = ["NHEJ_mean", "MHC_I"]
COVARIATES = ["PTPRC", "EPCAM"]
NEEDED = ["CLDN4", *ENDPOINTS, *COVARIATES]

ALIASES = {
    "CLDN4": ["CLDN4"],
    "PRKDC": ["PRKDC"],
    "LIG4": ["LIG4"],
    "STING1": ["STING1", "TMEM173"],
    "HLA-A": ["HLA-A"],
    "HLA-B": ["HLA-B"],
    "HLA-C": ["HLA-C"],
    "PTPRC": ["PTPRC"],
    "EPCAM": ["EPCAM"],
}

TCGA_ABBR = {
    "Acute Myeloid Leukemia": "LAML",
    "Adrenocortical Cancer": "ACC",
    "Bladder Urothelial Carcinoma": "BLCA",
    "Brain Lower Grade Glioma": "LGG",
    "Breast Invasive Carcinoma": "BRCA",
    "Cervical & Endocervical Cancer": "CESC",
    "Cholangiocarcinoma": "CHOL",
    "Colon Adenocarcinoma": "COAD",
    "Diffuse Large B-Cell Lymphoma": "DLBC",
    "Esophageal Carcinoma": "ESCA",
    "Glioblastoma Multiforme": "GBM",
    "Head & Neck Squamous Cell Carcinoma": "HNSC",
    "Kidney Chromophobe": "KICH",
    "Kidney Clear Cell Carcinoma": "KIRC",
    "Kidney Papillary Cell Carcinoma": "KIRP",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Lung Adenocarcinoma": "LUAD",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Mesothelioma": "MESO",
    "Ovarian Serous Cystadenocarcinoma": "OV",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Pheochromocytoma & Paraganglioma": "PCPG",
    "Prostate Adenocarcinoma": "PRAD",
    "Rectum Adenocarcinoma": "READ",
    "Sarcoma": "SARC",
    "Skin Cutaneous Melanoma": "SKCM",
    "Stomach Adenocarcinoma": "STAD",
    "Testicular Germ Cell Tumor": "TGCT",
    "Thymoma": "THYM",
    "Thyroid Carcinoma": "THCA",
    "Uterine Carcinosarcoma": "UCS",
    "Uterine Corpus Endometrioid Carcinoma": "UCEC",
    "Uveal Melanoma": "UVM",
}

LUNG_PRIMARY_IDS = ["TCGA-LUAD", "TCGA-LUSC", "GTEx-Lung"]
FLOOR = -9.9  # log2(0.001) is about -9.966; zeros are stored at that floor
MIN_N_TEST = 10
MIN_N_FDR = 30

TOIL_EXPR_URL = (
    "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz"
)
TOIL_PHENO_URL = (
    "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"
)
TOIL_PROBE_URL = (
    "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "probeMap/gencode.v23.annotation.gene.probemap"
)
HISEQ_URLS = {
    "LUAD": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
    "LUSC": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
}

COLORS = {
    "TCGA-LUAD": "#0072B2",
    "TCGA-LUSC": "#E69F00",
    "GTEx-Lung": "#009E73",
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: str, min_bytes: int) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print(f"  cached {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
        return
    last = None
    for attempt in range(5):
        try:
            req = Request(url, headers={"User-Agent": "gtex-tcga-cldn4-nhej/1.0"})
            tmp = dest + ".tmp"
            with urlopen(req, timeout=300) as r, open(tmp, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            if os.path.getsize(tmp) < min_bytes:
                raise RuntimeError(f"too small: {os.path.getsize(tmp)} bytes")
            os.replace(tmp, dest)
            print(f"  downloaded {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
            return
        except Exception as e:  # noqa: BLE001
            last = e
            print(f"  retry {attempt + 1} {url}: {e}")
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {dest}: {last}")


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg. Non-finite p-values stay non-finite and do not take a slot."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    qq = np.empty(n)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, pv[idx] * n / rank)
        qq[idx] = val
        prev = val
    q[np.flatnonzero(ok)] = qq
    return q


def fisher_z_ci(r: float, n: int, k_covariates: int = 0, alpha: float = 0.05):
    if r is None or not np.isfinite(r) or n is None or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = atanh(r)
    se = 1.0 / sqrt(n - k_covariates - 3)
    zcrit = float(stats.norm.ppf(1 - alpha / 2))
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def spearman_rho(x: np.ndarray, y: np.ndarray):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_TEST:
        return np.nan, np.nan, n
    xs, ys = x[m], y[m]
    if np.unique(xs).size < 2 or np.unique(ys).size < 2:
        return np.nan, np.nan, n
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=stats.ConstantInputWarning)
        rho, p = stats.spearmanr(xs, ys)
    return float(rho), float(p), n


def partial_rank_residual(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]):
    """Rank-residual Pearson. df = n - k - 2. k=1 matches the earlier ImmuneScore partial."""
    m = np.isfinite(x) & np.isfinite(y)
    for c in covariates:
        m = m & np.isfinite(c)
    n = int(m.sum())
    k = len(covariates)
    if n < k + MIN_N_TEST:
        return np.nan, np.nan, n
    xs, ys = x[m], y[m]
    if np.unique(xs).size < 2 or np.unique(ys).size < 2:
        return np.nan, np.nan, n
    xr = stats.rankdata(xs)
    yr = stats.rankdata(ys)
    cols = [np.ones(n)]
    for c in covariates:
        cr = stats.rankdata(c[m])
        if np.unique(cr).size < 2:
            return np.nan, np.nan, n
        cols.append(cr)
    z = np.column_stack(cols)
    bx, *_ = np.linalg.lstsq(z, xr, rcond=None)
    by, *_ = np.linalg.lstsq(z, yr, rcond=None)
    rx = xr - z @ bx
    ry = yr - z @ by
    if np.unique(rx).size < 2 or np.unique(ry).size < 2:
        return np.nan, np.nan, n
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - k - 2
    tstat = r * sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return r, p, n


def frac_floor(v: np.ndarray) -> float:
    x = v[np.isfinite(v)]
    if x.size == 0:
        return np.nan
    return float(np.mean(x < FLOOR))


def donor_id(sample: str, study: str) -> str:
    parts = str(sample).split("-")
    if study == "TCGA":
        return "-".join(parts[:3])
    if study == "GTEX":
        return "-".join(parts[:2])
    return str(sample)


def tcga_type_code(sample: str) -> str:
    parts = str(sample).split("-")
    if len(parts) < 4:
        return ""
    return parts[3][:2]


def self_check() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    rho, p, n = spearman_rho(x, x)
    assert n == 200 and abs(rho - 1) < 1e-9 and p < 1e-6
    z = rng.normal(size=400)
    # y shares only z with x, so the partial | z should be ~0
    x = z + rng.normal(size=400)
    y = z + rng.normal(size=400)
    r, p, n = partial_rank_residual(x, y, [z])
    assert n == 400 and abs(r) < 0.15
    q = bh_fdr(np.array([0.01, 0.04, np.nan, 0.5]))
    assert abs(q[0] - 0.03) < 1e-12
    assert abs(q[1] - 0.06) < 1e-12
    assert not np.isfinite(q[2])
    lo, hi = fisher_z_ci(0.0, 100, 0)
    assert lo < 0 < hi
    print("self-check ok")


def resolve_toil_ids(probe_path: str) -> dict[str, dict]:
    probe = pd.read_csv(probe_path, sep="\t")
    out = {}
    for canon, names in ALIASES.items():
        hit = probe[probe["gene"].isin(names)]
        if len(hit) != 1:
            raise SystemExit(f"probe map for {canon} ({names}) returned {len(hit)} rows")
        row = hit.iloc[0]
        out[canon] = {
            "ensembl": str(row["id"]),
            "matrix_symbol": str(row["gene"]),
            "chrom": str(row["chrom"]),
        }
    return out


def load_toil_genes(expr_path: str, idmap: dict[str, dict]) -> pd.DataFrame:
    want = {meta["ensembl"]: canon for canon, meta in idmap.items()}
    with gzip.open(expr_path, "rt") as f:
        samples = f.readline().rstrip("\n").split("\t")[1:]
        found: dict[str, np.ndarray] = {}
        n_rows = 0
        for line in f:
            n_rows += 1
            tab = line.find("\t")
            eid = line[:tab]
            canon = want.get(eid)
            if canon is None:
                continue
            vals = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float64)
            if vals.size != len(samples):
                raise SystemExit(f"{canon} parsed {vals.size} values, expected {len(samples)}")
            found[canon] = vals
            if len(found) == len(want):
                break
    missing = [g for g in NEEDED if g not in found]
    if missing:
        raise SystemExit(f"Toil matrix missing {missing} after scanning {n_rows} rows")
    print(f"  kept {len(found)} genes; scan stopped at row {n_rows}")
    df = pd.DataFrame(found, index=samples)
    df.index.name = "sample"
    return df


def load_symbol_genes(path: str) -> tuple[pd.DataFrame, dict[str, str]]:
    """HiSeqV2-style matrix: Hugo symbol in column 0. Rename aliases to canonical."""
    alias_to_canon = {}
    for canon, names in ALIASES.items():
        for name in names:
            alias_to_canon[name] = canon
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        samples = f.readline().rstrip("\n").split("\t")[1:]
        found: dict[str, np.ndarray] = {}
        used: dict[str, str] = {}
        for line in f:
            tab = line.find("\t")
            sym = line[:tab]
            canon = alias_to_canon.get(sym)
            if canon is None or canon in found:
                continue
            vals = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float64)
            if vals.size != len(samples):
                raise SystemExit(f"{sym} parsed {vals.size} values, expected {len(samples)}")
            found[canon] = vals
            used[canon] = sym
            if len(found) == len(NEEDED):
                break
    missing = [g for g in NEEDED if g not in found]
    if missing:
        raise SystemExit(f"{os.path.basename(path)} missing {missing}")
    df = pd.DataFrame(found, index=samples)
    df.index.name = "sample"
    return df, used


def add_companions(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["NHEJ_mean"] = out[["PRKDC", "LIG4"]].mean(axis=1)
    out["MHC_I"] = out[["HLA-A", "HLA-B", "HLA-C"]].mean(axis=1)
    return out


def collapse_donors(df: pd.DataFrame, study: str) -> pd.DataFrame:
    work = df.copy()
    work["donor"] = [donor_id(s, study) for s in work.index]
    value_cols = [c for c in work.columns if c != "donor"]
    n_ali = work.groupby("donor").size()
    means = work.groupby("donor")[value_cols].mean()
    means["n_aliquots"] = n_ali.astype(int)
    means.index.name = "donor"
    return means


def correlate_cohort(expr: pd.DataFrame, meta: dict) -> list[dict]:
    x = expr["CLDN4"].to_numpy(dtype=float)
    ptprc = expr["PTPRC"].to_numpy(dtype=float)
    epcam = expr["EPCAM"].to_numpy(dtype=float)
    r_ptprc, p_ptprc, n_ptprc = spearman_rho(x, ptprc)
    r_epcam, p_epcam, n_epcam = spearman_rho(x, epcam)
    rows = []
    for ep in ENDPOINTS + COMPANIONS:
        y = expr[ep].to_numpy(dtype=float)
        rho, p, n = spearman_rho(x, y)
        lo, hi = fisher_z_ci(rho, n, 0)
        r1, p1, n1 = partial_rank_residual(x, y, [ptprc])
        lo1, hi1 = fisher_z_ci(r1, n1, 1)
        r2, p2, n2 = partial_rank_residual(x, y, [epcam])
        lo2, hi2 = fisher_z_ci(r2, n2, 1)
        r3, p3, n3 = partial_rank_residual(x, y, [ptprc, epcam])
        lo3, hi3 = fisher_z_ci(r3, n3, 2)
        r_y_pt, p_y_pt, _ = spearman_rho(y, ptprc)
        r_y_ep, p_y_ep, _ = spearman_rho(y, epcam)
        rows.append(
            {
                **meta,
                "predictor": "CLDN4",
                "endpoint": ep,
                "in_primary_endpoint_set": ep in ENDPOINTS,
                "n": n,
                "rho": rho,
                "p": p,
                "rho_ci95_lo": lo,
                "rho_ci95_hi": hi,
                "n_partial_ptprc": n1,
                "rho_partial_ptprc": r1,
                "p_partial_ptprc": p1,
                "rho_partial_ptprc_ci95_lo": lo1,
                "rho_partial_ptprc_ci95_hi": hi1,
                "n_partial_epcam": n2,
                "rho_partial_epcam": r2,
                "p_partial_epcam": p2,
                "rho_partial_epcam_ci95_lo": lo2,
                "rho_partial_epcam_ci95_hi": hi2,
                "n_partial_ptprc_epcam": n3,
                "rho_partial_ptprc_epcam": r3,
                "p_partial_ptprc_epcam": p3,
                "rho_partial_ptprc_epcam_ci95_lo": lo3,
                "rho_partial_ptprc_epcam_ci95_hi": hi3,
                "rho_CLDN4_vs_PTPRC": r_ptprc,
                "p_CLDN4_vs_PTPRC": p_ptprc,
                "n_CLDN4_vs_PTPRC": n_ptprc,
                "rho_CLDN4_vs_EPCAM": r_epcam,
                "p_CLDN4_vs_EPCAM": p_epcam,
                "n_CLDN4_vs_EPCAM": n_epcam,
                "rho_endpoint_vs_PTPRC": r_y_pt,
                "p_endpoint_vs_PTPRC": p_y_pt,
                "rho_endpoint_vs_EPCAM": r_y_ep,
                "p_endpoint_vs_EPCAM": p_y_ep,
                "frac_cldn4_floor": frac_floor(x),
                "frac_endpoint_floor": frac_floor(y),
                "n_donors_collapsed": int(expr.shape[0]),
                "n_donors_with_replicates": int((expr["n_aliquots"] > 1).sum())
                if "n_aliquots" in expr.columns
                else 0,
            }
        )
    return rows


def assign_q(stats_df: pd.DataFrame) -> pd.DataFrame:
    df = stats_df.copy()
    df["q_lung"] = np.nan
    df["q_tcga"] = np.nan
    df["q_gtex"] = np.nan
    df["q_lung_adjacent"] = np.nan
    df["q_hiseqv2"] = np.nan
    df["q_lung_partial_ptprc"] = np.nan
    df["q_lung_partial_epcam"] = np.nan

    def _put(mask: pd.Series, pcol: str, qcol: str) -> None:
        idx = df.index[mask & df[pcol].notna() & (df["n"] >= MIN_N_FDR) & df["in_primary_endpoint_set"]]
        if pcol.startswith("p_partial"):
            # partial n can differ; require the partial n
            ncol = "n_partial_ptprc" if "ptprc" in pcol and "epcam" not in pcol else "n_partial_epcam"
            idx = df.index[mask & df[pcol].notna() & (df[ncol] >= MIN_N_FDR) & df["in_primary_endpoint_set"]]
        if len(idx) == 0:
            return
        df.loc[idx, qcol] = bh_fdr(df.loc[idx, pcol].to_numpy())

    toil = df["matrix"] == "Toil"
    _put(toil & df["cohort_id"].isin(LUNG_PRIMARY_IDS), "p", "q_lung")
    _put(toil & df["compartment"].isin(["TCGA_solid_primary", "TCGA_blood_primary", "TCGA_metastatic"]), "p", "q_tcga")
    _put(toil & (df["compartment"] == "GTEx_normal"), "p", "q_gtex")
    _put(toil & (df["compartment"] == "TCGA_adjacent_normal"), "p", "q_lung_adjacent")
    _put(df["matrix"] == "HiSeqV2", "p", "q_hiseqv2")
    _put(toil & df["cohort_id"].isin(LUNG_PRIMARY_IDS), "p_partial_ptprc", "q_lung_partial_ptprc")
    _put(toil & df["cohort_id"].isin(LUNG_PRIMARY_IDS), "p_partial_epcam", "q_lung_partial_epcam")
    return df


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def mag_word(rho: float) -> str:
    a = abs(rho)
    if a < 0.10:
        return "very small"
    if a < 0.20:
        return "small"
    if a < 0.30:
        return "modest"
    if a < 0.50:
        return "moderate"
    return "large"


def lookup(df: pd.DataFrame, cohort: str, endpoint: str, matrix: str = "Toil") -> pd.Series:
    hit = df[(df["cohort_id"] == cohort) & (df["endpoint"] == endpoint) & (df["matrix"] == matrix)]
    if len(hit) != 1:
        raise SystemExit(f"expected 1 row for {matrix} {cohort} {endpoint}, found {len(hit)}")
    return hit.iloc[0]


def lung_gene_paragraph(df: pd.DataFrame, gene: str) -> str:
    bits = []
    for cohort in LUNG_PRIMARY_IDS:
        r = lookup(df, cohort, gene)
        bits.append(
            f"{cohort} ρ={fmt_r(r['rho'])} (n={int(r['n'])}, p={fmt_p(r['p'])}, q={fmt_p(r['q_lung'])})"
        )
    rhos = [float(lookup(df, c, gene)["rho"]) for c in LUNG_PRIMARY_IDS]
    qs = [float(lookup(df, c, gene)["q_lung"]) for c in LUNG_PRIMARY_IDS]
    n_sig = sum(np.isfinite(q) and q < 0.05 for q in qs)
    signs = {np.sign(v) for v in rhos if np.isfinite(v) and v != 0}
    if n_sig == 3 and len(signs) == 1:
        word = "positive" if rhos[0] > 0 else "negative"
        lead = f"All three lung cohorts are {word} and survive the 18-test BH."
    elif n_sig == 0:
        lead = "None of the three survive the 18-test BH."
    else:
        lead = f"{n_sig} of 3 survive the 18-test BH."
    same = len(signs) == 1
    agree = "The sign agrees across the three cohorts." if same else "The sign does not agree across the three cohorts."
    return f"**{gene}.** {lead} {agree} " + "; ".join(bits) + "."


def partial_sentence(df: pd.DataFrame, gene: str) -> str:
    parts = []
    for cohort in LUNG_PRIMARY_IDS:
        r = lookup(df, cohort, gene)
        parts.append(
            f"{cohort} unadj {fmt_r(r['rho'])} → PTPRC {fmt_r(r['rho_partial_ptprc'])} "
            f"(q={fmt_p(r['q_lung_partial_ptprc'])}), EPCAM {fmt_r(r['rho_partial_epcam'])} "
            f"(q={fmt_p(r['q_lung_partial_epcam'])}), both {fmt_r(r['rho_partial_ptprc_epcam'])}"
        )
    return f"**{gene}.** " + "; ".join(parts) + "."


def sign_summary(df: pd.DataFrame, mask: pd.Series, qcol: str) -> pd.DataFrame:
    rows = []
    sub = df[mask & df["endpoint"].isin(ENDPOINTS)]
    for ep in ENDPOINTS:
        s = sub[sub["endpoint"] == ep]
        rho = s["rho"].to_numpy(dtype=float)
        q = s[qcol].to_numpy(dtype=float)
        ok = np.isfinite(rho)
        rows.append(
            {
                "endpoint": ep,
                "n_cohorts": int(ok.sum()),
                "median_rho": float(np.median(rho[ok])) if ok.any() else np.nan,
                "q25_rho": float(np.quantile(rho[ok], 0.25)) if ok.any() else np.nan,
                "q75_rho": float(np.quantile(rho[ok], 0.75)) if ok.any() else np.nan,
                "min_rho": float(np.min(rho[ok])) if ok.any() else np.nan,
                "max_rho": float(np.max(rho[ok])) if ok.any() else np.nan,
                "n_positive": int(np.sum(rho[ok] > 0)),
                "n_negative": int(np.sum(rho[ok] < 0)),
                "n_q_pos": int(np.sum((q < 0.05) & (rho > 0) & np.isfinite(q))),
                "n_q_neg": int(np.sum((q < 0.05) & (rho < 0) & np.isfinite(q))),
            }
        )
    return pd.DataFrame(rows)


def md_sign_table(summary: pd.DataFrame) -> list[str]:
    lines = [
        "| Endpoint | Cohorts | Median ρ | IQR | Positive | Exploratory q<0.05 and + | Exploratory q<0.05 and − |",
        "|---|---:|---:|---|---:|---:|---:|",
    ]
    for _, r in summary.iterrows():
        iqr = f"[{fmt_r(r['q25_rho'])}, {fmt_r(r['q75_rho'])}]"
        lines.append(
            f"| {r['endpoint']} | {int(r['n_cohorts'])} | {fmt_r(r['median_rho'])} | {iqr} | "
            f"{int(r['n_positive'])}/{int(r['n_cohorts'])} | {int(r['n_q_pos'])} | {int(r['n_q_neg'])} |"
        )
    return lines


def md_result_table(df: pd.DataFrame, cohorts: list[str], qcol: str, matrix: str = "Toil") -> list[str]:
    lines = [
        "| Cohort | Endpoint | n | ρ | 95% CI | p | q |",
        "|---|---|---:|---:|---|---:|---:|",
    ]
    for cohort in cohorts:
        for ep in ENDPOINTS:
            r = lookup(df, cohort, ep, matrix=matrix)
            ci = f"[{fmt_r(r['rho_ci95_lo'])}, {fmt_r(r['rho_ci95_hi'])}]"
            lines.append(
                f"| {cohort} | {ep} | {int(r['n'])} | {fmt_r(r['rho'])} | {ci} | "
                f"{fmt_p(r['p'])} | {fmt_p(r[qcol])} |"
            )
    return lines


def write_finding(
    stats_df: pd.DataFrame,
    counts: pd.DataFrame,
    idmap: dict,
    tcga_summary: pd.DataFrame,
    gtex_summary: pd.DataFrame,
) -> None:
    def n_of(cohort_id: str) -> int:
        hit = counts[counts["cohort_id"] == cohort_id]
        return int(hit.iloc[0]["n_donors"])

    luad_n = n_of("TCGA-LUAD")
    lusc_n = n_of("TCGA-LUSC")
    glung_n = n_of("GTEx-Lung")
    luad_norm = n_of("TCGA-LUAD-normal")
    lusc_norm = n_of("TCGA-LUSC-normal")
    h_luad = n_of("TCGA-LUAD-HiSeqV2")
    h_lusc = n_of("TCGA-LUSC-HiSeqV2")

    c_luad = counts[counts["cohort_id"] == "TCGA-LUAD"].iloc[0]
    c_lusc = counts[counts["cohort_id"] == "TCGA-LUSC"].iloc[0]
    c_glung = counts[counts["cohort_id"] == "GTEx-Lung"].iloc[0]

    sting_symbol = idmap["STING1"]["matrix_symbol"]
    sting_id = idmap["STING1"]["ensembl"]

    # Context: CLDN4 vs covariates in the three lung cohorts (same on every endpoint row).
    cov_lines = [
        "| Cohort | n | CLDN4 vs PTPRC ρ | p | CLDN4 vs EPCAM ρ | p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for cohort in LUNG_PRIMARY_IDS:
        r = lookup(stats_df, cohort, "PRKDC")
        cov_lines.append(
            f"| {cohort} | {int(r['n_CLDN4_vs_PTPRC'])} | {fmt_r(r['rho_CLDN4_vs_PTPRC'])} | "
            f"{fmt_p(r['p_CLDN4_vs_PTPRC'])} | {fmt_r(r['rho_CLDN4_vs_EPCAM'])} | "
            f"{fmt_p(r['p_CLDN4_vs_EPCAM'])} |"
        )

    # Floor QC
    floor_bits = []
    lung = stats_df[(stats_df["cohort_id"].isin(LUNG_PRIMARY_IDS)) & (stats_df["endpoint"].isin(ENDPOINTS))]
    for _, r in lung.iterrows():
        if r["frac_cldn4_floor"] > 0.10 or r["frac_endpoint_floor"] > 0.10:
            floor_bits.append(
                f"{r['cohort_id']} {r['endpoint']}: CLDN4 floor {r['frac_cldn4_floor']:.1%}, "
                f"endpoint floor {r['frac_endpoint_floor']:.1%}"
            )
    if floor_bits:
        floor_txt = "Lung cohorts with >10% of donors at the log2(TPM+0.001) floor: " + "; ".join(floor_bits) + "."
    else:
        floor_txt = (
            "No lung primary cohort has more than 10% of donors at the log2(TPM+0.001) floor "
            "for CLDN4 or for any of the six endpoints."
        )

    n_solid = int((counts["compartment"] == "TCGA_solid_primary").sum())
    n_gtex = int(((counts["compartment"] == "GTEx_normal") & (counts["n_donors"] >= MIN_N_FDR)).sum())

    lines: list[str] = []
    lines.append("# Finding — GTEx + TCGA pan-cancer: CLDN4 vs PRKDC / LIG4 / STING1 / HLA-A/B/C")
    lines.append("")
    lines.append(
        "**Additive public bulk RNA.** UCSC Xena Toil `TcgaTargetGtex_rsem_gene_tpm`, "
        "log2(RSEM TPM+0.001), GENCODE v23, one pipeline for TCGA and GTEx. "
        f"STING1 is **{sting_symbol}** on this freeze (`{sting_id}`). "
        "This folder does not re-audit CosMx exclusion, concordant-4, the keratin partials, or TISMO."
    )
    lines.append("")
    lines.append(
        "Primary question, pre-specified: within-cohort Spearman of **CLDN4** vs **PRKDC**, **LIG4**, "
        "**STING1**, **HLA-A**, **HLA-B**, and **HLA-C** in **TCGA-LUAD primary**, **TCGA-LUSC primary**, "
        "and **GTEx lung**. Those 18 tests are one BH family. Other TCGA primaries and other GTEx tissues "
        "are consistency context with their own BH families. Cohorts are not pooled into one correlation."
    )
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append("| Cohort | Samples before donor collapse | Donors (analysis n) | Donors with >1 aliquot |")
    lines.append("|---|---:|---:|---:|")
    lines.append(
        f"| TCGA-LUAD primary | {int(c_luad['n_samples'])} | **{luad_n}** | {int(c_luad['n_donors_with_replicates'])} |"
    )
    lines.append(
        f"| TCGA-LUSC primary | {int(c_lusc['n_samples'])} | **{lusc_n}** | {int(c_lusc['n_donors_with_replicates'])} |"
    )
    lines.append(
        f"| GTEx lung | {int(c_glung['n_samples'])} | **{glung_n}** | {int(c_glung['n_donors_with_replicates'])} |"
    )
    lines.append(f"| TCGA-LUAD solid-tissue normal |  | {luad_norm} |  |")
    lines.append(f"| TCGA-LUSC solid-tissue normal |  | {lusc_norm} |  |")
    lines.append(f"| TCGA-LUAD HiSeqV2 primary (sensitivity) |  | {h_luad} |  |")
    lines.append(f"| TCGA-LUSC HiSeqV2 primary (sensitivity) |  | {h_lusc} |  |")
    lines.append("")
    lines.append(
        "Analysis n is donors after averaging replicate aliquots. No imputation. "
        f"Pan-cancer context uses {n_solid} TCGA solid-primary cohorts plus LAML (blood primary) and "
        "SKCM metastatic (most SKCM RNA is metastatic; primary SKCM is its own solid row). "
        f"GTEx context uses {n_gtex} normal tissues with at least {MIN_N_FDR} donors. "
        "TARGET samples are not included. GTEx cell lines are not included."
    )
    lines.append("")
    lines.append("## Lung primary result")
    lines.append("")
    lines.append(
        "BH q is across these 18 tests only. A q on this table is not the pan-cancer q."
    )
    lines.append("")
    lines.extend(md_result_table(stats_df, LUNG_PRIMARY_IDS, "q_lung"))
    lines.append("")
    for gene in ENDPOINTS:
        lines.append(lung_gene_paragraph(stats_df, gene))
        lines.append("")
    lines.append(floor_txt)
    lines.append("")
    lines.append("## Composition sensitivity (lung 18, separate BH)")
    lines.append("")
    lines.append(
        "Partials are rank-residual Pearson on the same Toil matrix. "
        "**PTPRC** is a leukocyte marker. **EPCAM** is an epithelial marker. "
        "Neither covariate is inside the endpoint definitions. "
        "Residualizing RNA does not equal adjusting for DNA purity. "
        "An EPCAM partial that shrinks is shared epithelial variance, not by itself proof that the "
        "unadjusted association was only composition. q below is BH within the 18 partials of that covariate, "
        "not the unadjusted 18."
    )
    lines.append("")
    lines.append("CLDN4 versus the two covariates (context for the partials):")
    lines.append("")
    lines.extend(cov_lines)
    lines.append("")
    for gene in ENDPOINTS:
        lines.append(partial_sentence(stats_df, gene))
        lines.append("")
    lines.append("## HiSeqV2 lung freeze (sensitivity)")
    lines.append("")
    lines.append(
        "Legacy UCSC Xena HiSeqV2 log2(norm_count+1), primary tumors (`-01`), donor-averaged. "
        "This is the freeze earlier LUAD/LUSC HLA tables used. It is not the Toil TPM matrix. "
        "BH here is a separate 12-test family (2 cohorts × 6 genes)."
    )
    lines.append("")
    lines.extend(
        md_result_table(
            stats_df,
            ["TCGA-LUAD-HiSeqV2", "TCGA-LUSC-HiSeqV2"],
            "q_hiseqv2",
            matrix="HiSeqV2",
        )
    )
    lines.append("")
    # One comparison sentence per gene: Toil vs HiSeqV2 sign in LUAD and LUSC.
    agree_bits = []
    for gene in ENDPOINTS:
        flags = []
        for toil_id, hseq_id in (("TCGA-LUAD", "TCGA-LUAD-HiSeqV2"), ("TCGA-LUSC", "TCGA-LUSC-HiSeqV2")):
            a = float(lookup(stats_df, toil_id, gene)["rho"])
            b = float(lookup(stats_df, hseq_id, gene, matrix="HiSeqV2")["rho"])
            flags.append(np.sign(a) == np.sign(b) or abs(a) < 0.05 and abs(b) < 0.05)
        agree_bits.append(f"{gene} {'sign-agrees' if all(flags) else 'sign differs'} in both histologies")
    lines.append(
        "Toil vs HiSeqV2, LUAD and LUSC, calling |ρ|<0.05 on both a match: " + "; ".join(agree_bits) + "."
    )
    lines.append("")
    lines.append("## TCGA pan-cancer context")
    lines.append("")
    lines.append(
        f"Each solid-primary cohort is its own Spearman (n≥{MIN_N_FDR} to enter BH). "
        "Exploratory BH is across TCGA cohorts × the six genes, and it includes LUAD, LUSC, "
        "LAML blood, and SKCM metastatic as separate rows. Median and IQR below are **solid primaries only** "
        "(not LAML, not SKCM metastatic). This is a description of heterogeneity, not a single pan-cancer p-value."
    )
    lines.append("")
    solid_ids = set(counts.loc[counts["compartment"] == "TCGA_solid_primary", "cohort_id"])
    lines.extend(md_sign_table(tcga_summary))
    lines.append("")
    lines.append(
        f"Solid-primary cohort count in that table: {len(solid_ids)}. "
        "Full per-cohort ρ is in `tables/correlations.tsv` and `figures/heatmap_tcga.png`."
    )
    lines.append("")
    lines.append("## GTEx tissue context")
    lines.append("")
    lines.append(
        f"Normal tissues with at least {MIN_N_FDR} donors. Cell lines are excluded. "
        "Brain subregions stay separate. BH is within this GTEx family. Lung is included here as one tissue; "
        "its primary q is the lung-family q above, not this one."
    )
    lines.append("")
    lines.extend(md_sign_table(gtex_summary))
    lines.append("")
    lines.append("Full per-tissue ρ is in `tables/correlations.tsv` and `figures/heatmap_gtex.png`.")
    lines.append("")
    lines.append("## Lung adjacent normal")
    lines.append("")
    lines.append(
        "TCGA solid-tissue normal, not pooled with GTEx and not part of the 18. "
        "These donors overlap the tumor cohorts as paired normals; the correlation is within the normal samples, "
        "not a tumor-minus-normal delta. BH is a separate 12-test family."
    )
    lines.append("")
    lines.extend(
        md_result_table(
            stats_df,
            ["TCGA-LUAD-normal", "TCGA-LUSC-normal"],
            "q_lung_adjacent",
        )
    )
    lines.append("")
    lines.append("## Companions (not in any FDR family)")
    lines.append("")
    lines.append(
        "NHEJ_mean is the unweighted mean of PRKDC and LIG4 on the log2 TPM scale. "
        "MHC_I is the unweighted mean of HLA-A, HLA-B, and HLA-C. Neither is a gene-set enrichment."
    )
    lines.append("")
    lines.append("| Cohort | Companion | n | ρ | p |")
    lines.append("|---|---|---:|---:|---:|")
    for cohort in LUNG_PRIMARY_IDS:
        for ep in COMPANIONS:
            r = lookup(stats_df, cohort, ep)
            lines.append(
                f"| {cohort} | {ep} | {int(r['n'])} | {fmt_r(r['rho'])} | {fmt_p(r['p'])} |"
            )
    lines.append("")
    lines.append("## What this measurement is")
    lines.append("")
    lines.append(
        "These are within-cohort rank correlations on bulk RNA. A positive CLDN4–PRKDC or CLDN4–LIG4 ρ "
        "means the two transcripts move together across donors in that tissue or tumor type. "
        "A CLDN4–HLA or CLDN4–STING1 ρ mixes tumor-cell and immune-cell RNA. "
        "It is not a tumor-cell program, not a spatial neighborhood, and not an ICI endpoint. "
        "Negative HLA ρ is not spatial exclusion. Positive NHEJ ρ is not a repair mechanism."
    )
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "- **Toil matrix:** `TcgaTargetGtex_rsem_gene_tpm`, log2(TPM+0.001). "
        "Phenotype `TcgaTargetGTEX_phenotype.txt.gz`. "
        "Gene map `gencode.v23.annotation.gene.probemap` (Ensembl ID with version). "
        "Zeros are stored at log2(0.001) ≈ −9.966; Spearman treats that floor as a tie."
    )
    lines.append(
        "- **Cohorts:** TCGA `_sample_type == Primary Tumor` split by `detailed_category`, "
        "plus LAML blood primary and SKCM metastatic as labeled extra rows. "
        "GTEx `_sample_type == Normal Tissue` split by `detailed_category`. "
        "Lung adjacent: TCGA `_sample_type == Solid Tissue Normal` for LUAD and LUSC only."
    )
    lines.append(
        "- **Donor collapse:** TCGA patient = first three barcode fields; GTEx donor = first two. "
        "Mean of aliquots. Spearman is unchanged by a monotone transform, so log2 TPM and TPM agree."
    )
    lines.append(
        "- **Unadjusted:** Spearman, two-sided. **CI:** Fisher z, variance 1/(n−3). "
        f"Tests with n<{MIN_N_TEST} are not computed. BH uses cohorts with n≥{MIN_N_FDR}."
    )
    lines.append(
        "- **Partial:** rank-residual Pearson given PTPRC, EPCAM, or both. "
        "t degrees of freedom n−k−2. Fisher z variance 1/(n−k−3). "
        "Algebraic first-order partial Spearman was the formula in the earlier ImmuneScore lung table; "
        "rank-residual Pearson is the sensitivity used here because it extends to two covariates. "
        "On a single covariate the two agree closely when ties are mild."
    )
    lines.append(
        "- **FDR families, kept separate:** (1) lung Toil 18, (2) TCGA exploratory cohorts × 6, "
        "(3) GTEx tissues × 6, (4) lung adjacent 12, (5) HiSeqV2 lung 12, "
        "(6) lung PTPRC partials 18, (7) lung EPCAM partials 18. "
        "Companions are outside every family."
    )
    lines.append(
        "- **HiSeqV2:** `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`, "
        "primary `01` only, collapsed to the 15-character barcode by mean, then to the patient."
    )
    lines.append("")
    lines.append("## What is not done")
    lines.append("")
    lines.append(
        "- No pooled GTEx+TCGA correlation. No STAR-TPM pan-cancer re-run. "
        "No ESTIMATE or leukocyte-fraction residual (PTPRC/EPCAM are the pre-specified RNA covariates). "
        "No ICI outcome, protein, spatial statistic, or single-cell pseudobulk. "
        "No claim that this ranks CLDN4 against other surface genes."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/correlations.tsv` — every cohort, endpoint, partial, and q")
    lines.append("- `tables/sign_summary.tsv` — median ρ and sign counts")
    lines.append("- `tables/counts.tsv` — honest n per cohort")
    lines.append("- `tables/samples.tsv` — donor-level matrix the correlations were computed from")
    lines.append("- `tables/provenance.json` — URLs, sha256, Ensembl IDs")
    lines.append("- `figures/forest_lung.png` — lung unadjusted and PTPRC-partial ρ")
    lines.append("- `figures/heatmap_tcga.png` — TCGA per-cohort ρ")
    lines.append("- `figures/heatmap_gtex.png` — GTEx per-tissue ρ")
    lines.append("- Reproduce: `python3 methods/gtex_tcga_cldn4_nhej/analyze.py`")
    lines.append("")
    text = "\n".join(lines)
    path = os.path.join(HERE, "FINDING.md")
    with open(path, "w") as f:
        f.write(text)
    print(f"wrote {path}")


def forest_panel(ax, df: pd.DataFrame, rho_col: str, lo_col: str, hi_col: str, title: str) -> None:
    offsets = {"TCGA-LUAD": -0.22, "TCGA-LUSC": 0.0, "GTEx-Lung": 0.22}
    ybase = np.arange(len(ENDPOINTS))
    for cohort in LUNG_PRIMARY_IDS:
        ys, xs, xerr_lo, xerr_hi = [], [], [], []
        for i, ep in enumerate(ENDPOINTS):
            r = lookup(df, cohort, ep)
            rho = float(r[rho_col])
            lo = float(r[lo_col])
            hi = float(r[hi_col])
            if not np.isfinite(rho):
                continue
            ys.append(ybase[i] + offsets[cohort])
            xs.append(rho)
            xerr_lo.append(rho - lo if np.isfinite(lo) else 0)
            xerr_hi.append(hi - rho if np.isfinite(hi) else 0)
        ax.errorbar(
            xs,
            ys,
            xerr=[xerr_lo, xerr_hi],
            fmt="o",
            color=COLORS[cohort],
            label=cohort,
            ms=5,
            capsize=2,
            lw=1,
        )
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ybase)
    ax.set_yticklabels(ENDPOINTS)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman ρ (95% Fisher z CI)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, loc="best")


def draw_forest(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), sharey=True)
    forest_panel(axes[0], df, "rho", "rho_ci95_lo", "rho_ci95_hi", "Unadjusted")
    forest_panel(
        axes[1],
        df,
        "rho_partial_ptprc",
        "rho_partial_ptprc_ci95_lo",
        "rho_partial_ptprc_ci95_hi",
        "Partial | PTPRC",
    )
    fig.suptitle("CLDN4 in lung — Toil log2(TPM+0.001)", fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"forest_lung.{ext}"), dpi=160)
    plt.close(fig)


def draw_heatmap(df: pd.DataFrame, cohort_ids: list[str], qcol: str, title: str, filename: str) -> None:
    mat = np.full((len(cohort_ids), len(ENDPOINTS)), np.nan)
    annot = [[""] * len(ENDPOINTS) for _ in cohort_ids]
    ylabels = []
    for i, cid in enumerate(cohort_ids):
        n = int(lookup(df, cid, ENDPOINTS[0])["n"])
        ylabels.append(f"{cid} (n={n})")
        for j, ep in enumerate(ENDPOINTS):
            r = lookup(df, cid, ep)
            rho = float(r["rho"])
            mat[i, j] = rho
            star = ""
            q = r[qcol]
            if np.isfinite(q) and q < 0.05:
                star = "*"
            annot[i][j] = f"{rho:+.2f}{star}" if np.isfinite(rho) else ""
    finite = mat[np.isfinite(mat)]
    lim = 0.4 if finite.size == 0 else float(np.clip(np.ceil(np.quantile(np.abs(finite), 0.98) * 10) / 10, 0.3, 0.8))
    fig_h = max(4.5, 0.28 * len(cohort_ids) + 1.4)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(np.arange(len(ENDPOINTS)))
    ax.set_xticklabels(ENDPOINTS, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(cohort_ids)))
    ax.set_yticklabels(ylabels, fontsize=7)
    for i in range(len(cohort_ids)):
        for j in range(len(ENDPOINTS)):
            color = "white" if np.isfinite(mat[i, j]) and abs(mat[i, j]) > 0.55 * lim else "black"
            ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=6, color=color)
    ax.set_title(title + f"\n* exploratory/family q<0.05; color scale ±{lim:.1f}")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{filename}.{ext}"), dpi=160)
    plt.close(fig)


def prepare_toil() -> tuple[pd.DataFrame, dict, dict]:
    expr_path = os.path.join(DATA, "TcgaTargetGtex_rsem_gene_tpm.gz")
    pheno_path = os.path.join(DATA, "TcgaTargetGTEX_phenotype.txt.gz")
    probe_path = os.path.join(DATA, "gencode.v23.annotation.gene.probemap")
    download(TOIL_EXPR_URL, expr_path, min_bytes=1_200_000_000)
    download(TOIL_PHENO_URL, pheno_path, min_bytes=50_000)
    download(TOIL_PROBE_URL, probe_path, min_bytes=1_000_000)
    idmap = resolve_toil_ids(probe_path)
    print("gene map:")
    for k, v in idmap.items():
        print(f"  {k}: {v['matrix_symbol']} {v['ensembl']}")
    expr = load_toil_genes(expr_path, idmap)
    pheno = pd.read_csv(pheno_path, sep="\t", encoding="latin1").set_index("sample")
    if expr.index.duplicated().any() or pheno.index.duplicated().any():
        raise SystemExit("duplicate sample IDs in Toil expression or phenotype")
    both = expr.index.intersection(pheno.index)
    if len(both) < 0.99 * len(expr.index):
        raise SystemExit(f"expression/phenotype overlap {len(both)} / {len(expr.index)}")
    df = expr.join(pheno, how="inner")
    unknown = sorted(set(df.loc[df["_study"] == "TCGA", "detailed_category"]) - set(TCGA_ABBR))
    # phenotype includes normals and metastases under the same disease names; the abbr map covers diseases.
    if unknown:
        raise SystemExit(f"TCGA detailed_category missing from abbr map: {unknown}")
    prov = {
        "expr_url": TOIL_EXPR_URL,
        "pheno_url": TOIL_PHENO_URL,
        "probe_url": TOIL_PROBE_URL,
        "expr_sha256": sha256(expr_path),
        "pheno_sha256": sha256(pheno_path),
        "probe_sha256": sha256(probe_path),
        "expr_bytes": os.path.getsize(expr_path),
        "n_expr_samples": int(expr.shape[0]),
        "n_joined_samples": int(df.shape[0]),
        "genes": idmap,
        "value_scale": "log2(RSEM TPM + 0.001)",
        "floor_flag": FLOOR,
    }
    return df, idmap, prov


def cohort_from_subset(sub: pd.DataFrame, meta: dict, study: str) -> tuple[pd.DataFrame, dict]:
    sub = sub.copy()
    n_samples = int(sub.shape[0])
    if n_samples == 0:
        raise SystemExit(f"empty cohort {meta['cohort_id']}")
    base = add_companions(collapse_donors(sub[NEEDED], study))
    count = {
        "cohort_id": meta["cohort_id"],
        "matrix": meta["matrix"],
        "compartment": meta["compartment"],
        "in_lung_primary": meta["in_lung_primary"],
        "n_samples": n_samples,
        "n_donors": int(base.shape[0]),
        "n_donors_with_replicates": int((base["n_aliquots"] > 1).sum()),
        "label": meta["label"],
    }
    sample_out = base.copy()
    sample_out.insert(0, "cohort_id", meta["cohort_id"])
    sample_out.insert(1, "matrix", meta["matrix"])
    sample_out.insert(2, "compartment", meta["compartment"])
    return sample_out, count | {"_expr": base, "_meta": meta}


def build_toil_cohorts(df: pd.DataFrame) -> tuple[list[dict], list[pd.DataFrame]]:
    cohorts = []
    sample_frames = []
    tcga = df["_study"] == "TCGA"
    gtex = df["_study"] == "GTEX"

    # Solid primaries, including LUAD and LUSC.
    prim = tcga & (df["_sample_type"] == "Primary Tumor")
    for disease, sub in df.loc[prim].groupby("detailed_category"):
        abbr = TCGA_ABBR[disease]
        code_ok = sub.index.map(tcga_type_code) == "01"
        if not bool(code_ok.all()):
            bad = int((~code_ok).sum())
            print(f"  dropping {bad} {abbr} primary rows whose barcode is not type 01")
            sub = sub.loc[code_ok]
        meta = {
            "cohort_id": f"TCGA-{abbr}",
            "matrix": "Toil",
            "compartment": "TCGA_solid_primary",
            "in_lung_primary": abbr in {"LUAD", "LUSC"},
            "label": f"TCGA-{abbr} primary tumor",
            "disease": disease,
        }
        frame, packed = cohort_from_subset(sub, meta, "TCGA")
        sample_frames.append(frame.reset_index())
        cohorts.append(packed)

    # LAML blood primary.
    laml = tcga & (df["detailed_category"] == "Acute Myeloid Leukemia") & (
        df["_sample_type"] == "Primary Blood Derived Cancer - Peripheral Blood"
    )
    meta = {
        "cohort_id": "TCGA-LAML",
        "matrix": "Toil",
        "compartment": "TCGA_blood_primary",
        "in_lung_primary": False,
        "label": "TCGA-LAML blood primary",
        "disease": "Acute Myeloid Leukemia",
    }
    frame, packed = cohort_from_subset(df.loc[laml], meta, "TCGA")
    sample_frames.append(frame.reset_index())
    cohorts.append(packed)

    # SKCM metastatic (companion row; primary SKCM is already in the solid set).
    skcm_met = tcga & (df["detailed_category"] == "Skin Cutaneous Melanoma") & (df["_sample_type"] == "Metastatic")
    meta = {
        "cohort_id": "TCGA-SKCM-metastatic",
        "matrix": "Toil",
        "compartment": "TCGA_metastatic",
        "in_lung_primary": False,
        "label": "TCGA-SKCM metastatic",
        "disease": "Skin Cutaneous Melanoma",
    }
    frame, packed = cohort_from_subset(df.loc[skcm_met], meta, "TCGA")
    sample_frames.append(frame.reset_index())
    cohorts.append(packed)

    # Lung adjacent normals.
    for disease, abbr in (
        ("Lung Adenocarcinoma", "LUAD"),
        ("Lung Squamous Cell Carcinoma", "LUSC"),
    ):
        mask = (
            tcga
            & (df["detailed_category"] == disease)
            & (df["_sample_type"] == "Solid Tissue Normal")
        )
        meta = {
            "cohort_id": f"TCGA-{abbr}-normal",
            "matrix": "Toil",
            "compartment": "TCGA_adjacent_normal",
            "in_lung_primary": False,
            "label": f"TCGA-{abbr} solid tissue normal",
            "disease": disease,
        }
        frame, packed = cohort_from_subset(df.loc[mask], meta, "TCGA")
        sample_frames.append(frame.reset_index())
        cohorts.append(packed)

    # GTEx normal tissues. Drop cell lines via sample type.
    gtex_n = df.loc[gtex & (df["_sample_type"] == "Normal Tissue")]
    for tissue, sub in gtex_n.groupby("detailed_category"):
        safe = str(tissue)
        meta = {
            "cohort_id": f"GTEx-{safe}",
            "matrix": "Toil",
            "compartment": "GTEx_normal",
            "in_lung_primary": safe == "Lung",
            "label": f"GTEx {safe}",
            "disease": safe,
        }
        frame, packed = cohort_from_subset(sub, meta, "GTEX")
        sample_frames.append(frame.reset_index())
        cohorts.append(packed)
    ids = [c["cohort_id"] for c in cohorts]
    if len(ids) != len(set(ids)):
        dup = sorted({i for i in ids if ids.count(i) > 1})
        raise SystemExit(f"duplicate cohort ids: {dup}")
    return cohorts, sample_frames


def prepare_hiseqv2() -> tuple[list[dict], list[pd.DataFrame], dict]:
    packed_list = []
    frames = []
    prov = {}
    for abbr, url in HISEQ_URLS.items():
        dest = os.path.join(DATA, f"TCGA.{abbr}.HiSeqV2.gz")
        download(url, dest, min_bytes=1_000_000)
        expr, used = load_symbol_genes(dest)
        keep = [is_primary_barcode(s) for s in expr.index]
        expr = expr.loc[keep].copy()
        expr.index = ["-".join(str(s).replace(".", "-").split("-")[:4])[:15] for s in expr.index]
        expr = expr.groupby(level=0).mean()
        n_samples = int(expr.shape[0])
        base = collapse_donors(expr, "TCGA")
        base = add_companions(base)
        meta = {
            "cohort_id": f"TCGA-{abbr}-HiSeqV2",
            "matrix": "HiSeqV2",
            "compartment": "TCGA_solid_primary",
            "in_lung_primary": False,
            "label": f"TCGA-{abbr} HiSeqV2 primary",
            "disease": abbr,
        }
        frame = base.copy()
        frame.insert(0, "cohort_id", meta["cohort_id"])
        frame.insert(1, "matrix", meta["matrix"])
        frame.insert(2, "compartment", meta["compartment"])
        frames.append(frame.reset_index())
        packed_list.append(
            {
                "cohort_id": meta["cohort_id"],
                "matrix": "HiSeqV2",
                "compartment": "TCGA_solid_primary",
                "in_lung_primary": False,
                "n_samples": n_samples,
                "n_donors": int(base.shape[0]),
                "n_donors_with_replicates": int((base["n_aliquots"] > 1).sum()),
                "label": meta["label"],
                "_expr": base,
                "_meta": meta,
                "symbols_used": used,
            }
        )
        prov[abbr] = {"url": url, "sha256": sha256(dest), "bytes": os.path.getsize(dest), "symbols": used}
        print(f"  HiSeqV2 {abbr}: samples {n_samples}, donors {base.shape[0]}, symbols {used}")
    return packed_list, frames, prov


def is_primary_barcode(barcode: str) -> bool:
    parts = str(barcode).replace(".", "-").split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def heatmap_order(counts: pd.DataFrame, compartment: str, lung_first: list[str]) -> list[str]:
    sub = counts[counts["compartment"] == compartment].copy()
    if compartment == "GTEx_normal":
        sub = sub[sub["n_donors"] >= MIN_N_FDR]
    first = [c for c in lung_first if c in set(sub["cohort_id"])]
    rest = sub[~sub["cohort_id"].isin(first)].sort_values(["n_donors", "cohort_id"], ascending=[False, True])
    return first + rest["cohort_id"].tolist()


def main() -> None:
    self_check()
    print("loading Toil")
    toil, idmap, prov = prepare_toil()
    print("building cohorts")
    cohorts, sample_frames = build_toil_cohorts(toil)
    print("loading HiSeqV2 lung sensitivity")
    hiseq_cohorts, hiseq_frames, hiseq_prov = prepare_hiseqv2()
    cohorts.extend(hiseq_cohorts)
    sample_frames.extend(hiseq_frames)

    rows = []
    count_rows = []
    for packed in cohorts:
        meta = packed["_meta"]
        rows.extend(correlate_cohort(packed["_expr"], meta))
        count_rows.append({k: v for k, v in packed.items() if not k.startswith("_") and k != "symbols_used"})
    stats_df = assign_q(pd.DataFrame(rows))
    counts = pd.DataFrame(count_rows)

    # Sanity: lung filters actually kept the tumors.
    for cid, lo in (("TCGA-LUAD", 400), ("TCGA-LUSC", 400), ("GTEx-Lung", 200)):
        n = int(counts.loc[counts["cohort_id"] == cid, "n_donors"].iloc[0])
        if n < lo:
            raise SystemExit(f"{cid} n={n} is below the sanity floor {lo}")

    solid_mask = (stats_df["matrix"] == "Toil") & (stats_df["compartment"] == "TCGA_solid_primary")
    gtex_mask = (
        (stats_df["matrix"] == "Toil")
        & (stats_df["compartment"] == "GTEx_normal")
        & (stats_df["n"] >= MIN_N_FDR)
    )
    tcga_summary = sign_summary(stats_df, solid_mask, "q_tcga")
    tcga_summary.insert(0, "set", "TCGA_solid_primary")
    gtex_summary = sign_summary(stats_df, gtex_mask, "q_gtex")
    gtex_summary.insert(0, "set", "GTEx_normal_n_ge_30")
    summary = pd.concat([tcga_summary, gtex_summary], ignore_index=True)

    stats_df.to_csv(os.path.join(TABLES, "correlations.tsv"), sep="\t", index=False)
    counts.sort_values(["matrix", "compartment", "cohort_id"]).to_csv(
        os.path.join(TABLES, "counts.tsv"), sep="\t", index=False
    )
    summary.to_csv(os.path.join(TABLES, "sign_summary.tsv"), sep="\t", index=False)
    samples = pd.concat(sample_frames, ignore_index=True)
    samples.to_csv(os.path.join(TABLES, "samples.tsv"), sep="\t", index=False)

    prov["hiseqv2"] = hiseq_prov
    prov["n_correlation_rows"] = int(stats_df.shape[0])
    prov["n_sample_rows"] = int(samples.shape[0])
    prov["primary_family"] = "18 = 3 lung cohorts x 6 genes, Toil, BH"
    prov["generated_by"] = "methods/gtex_tcga_cldn4_nhej/analyze.py"
    with open(os.path.join(TABLES, "provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)

    draw_forest(stats_df)
    tcga_ids = heatmap_order(counts[counts["matrix"] == "Toil"], "TCGA_solid_primary", ["TCGA-LUAD", "TCGA-LUSC"])
    # Append labeled non-solid TCGA rows under the solid block.
    extra = ["TCGA-LAML", "TCGA-SKCM-metastatic"]
    draw_heatmap(
        stats_df,
        tcga_ids + extra,
        "q_tcga",
        "TCGA CLDN4 Spearman ρ (Toil). LUAD/LUSC at top; LAML blood and SKCM metastatic at bottom.",
        "heatmap_tcga",
    )
    gtex_ids = heatmap_order(counts, "GTEx_normal", ["GTEx-Lung"])
    draw_heatmap(
        stats_df,
        gtex_ids,
        "q_gtex",
        "GTEx normal tissues, CLDN4 Spearman ρ (Toil), n≥30 donors. Lung at top.",
        "heatmap_gtex",
    )
    write_finding(stats_df, counts, idmap, tcga_summary, gtex_summary)

    print("\n=== LUNG PRIMARY ===")
    show = stats_df[(stats_df["cohort_id"].isin(LUNG_PRIMARY_IDS)) & (stats_df["endpoint"].isin(ENDPOINTS))]
    cols = ["cohort_id", "endpoint", "n", "rho", "p", "q_lung", "rho_partial_ptprc", "rho_partial_epcam"]
    print(show[cols].to_string(index=False))
    print("\n=== SIGN SUMMARY ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
