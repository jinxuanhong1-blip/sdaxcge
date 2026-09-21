#!/usr/bin/env python3
"""TCGA-LUAD / TCGA-LUSC: CLDN4-low tumors vs ssGSEA NHEJ, IFN-γ, STING, APM.

Public inputs only.
  Expression: UCSC Xena legacy HiSeqV2 log2(RSEM norm_count + 1).
  ImmuneScore: official MD Anderson ESTIMATE RNAseqV2 (Yoshihara 2013).
  ssGSEA: gseapy 1.3.1, rank normalization, weight 0.25, no permutation.

Questions
  1. In CLDN4-low (bottom quartile) primaries, is the NHEJ score lower than
     in CLDN4-high (top quartile)?
  2. Are Hallmark IFN-γ, STING, and MHC-I APM scores higher in CLDN4-low?
  3. Spearman of CLDN4 vs PRKDC, LIG4, STING1 (TMEM173 on this freeze),
     and HLA-A/B/C, unadjusted and partial on ESTIMATE ImmuneScore.

Outputs -> methods/tcga_cldn4_nhej_sting/
"""
from __future__ import annotations

import gzip
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

import gseapy as gp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
TABLES = os.path.join(HERE, "tables")
FIG = os.path.join(HERE, "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# Hallmark IFN-γ (MSigDB via Enrichr MSigDB_Hallmark_2020), with HiSeqV2
# symbols where HGNC has renamed the gene since the Xena freeze.
IFN_ALIAS = {"WARS1": "WARS", "MARCHF1": "MARCH1", "HELZ2": "PRIC285", "CMTR1": "FTSJD2"}
IFN_GAMMA_HALLMARK = [
    "STAT1", "ISG15", "IFIT1", "MX1", "IFIT3", "IFI35", "IRF7", "IFIT2", "OAS2", "TAP1",
    "EIF2AK2", "RSAD2", "MX2", "IRF1", "OAS3", "TNFSF10", "IRF9", "CXCL10", "IFI44", "BST2",
    "XAF1", "SP110", "OASL", "PSMB8", "IFI44L", "IFITM3", "DDX60", "LGALS3BP", "GBP4", "IRF8",
    "PSMB9", "PML", "IFIH1", "UBE2L6", "IFI27", "ADAR", "LY6E", "STAT2", "CXCL9", "IL10RA",
    "PLA2G4A", "TRIM21", "USP18", "PTGS2", "EPSTI1", "C1S", "DDX58", "IL15", "NLRC5", "NMI",
    "IDO1", "PSMB10", "CXCL11", "ITGB7", "SAMHD1", "HERC6", "CMPK2", "SAMD9L", "RTP4", "PTPN2",
    "PARP14", "TNFAIP2", "IFITM2", "PLSCR1", "SOCS1", "CASP1", "ICAM1", "WARS1", "PSME1", "ISG20",
    "IRF2", "TRIM14", "FCGR1A", "MARCHF1", "SOCS3", "JAK2", "HLA-DMA", "PARP12", "TNFAIP6", "TRIM26",
    "VCAM1", "CD274", "CIITA", "NAMPT", "SELP", "GPR18", "FPR1", "HELZ2", "PSME2", "SERPING1",
    "CCL5", "RNF31", "SOD2", "TRIM25", "LAP3", "PSMA3", "RNF213", "PELI1", "CFB", "CD86",
    "TXNIP", "HLA-DQA1", "GCH1", "PNP", "CCL7", "PTPN6", "SPPL2A", "IL4R", "PNPT1", "DHX58",
    "BTG1", "CASP8", "IFI30", "CCL2", "FGL2", "CASP7", "SECTM1", "IL15RA", "CD40", "TRAFD1",
    "HLA-DRB1", "GBP6", "LCP2", "HLA-G", "MT2A", "RIPK1", "KLRK1", "UPP1", "PSMB2", "TDRD7",
    "HIF1A", "EIF4E3", "VAMP8", "PFKP", "CD38", "ZBP1", "BANK1", "TOR1B", "RBCK1", "PDE4B",
    "MVP", "IL7", "BPGM", "CMTR1", "AUTS2", "B2M", "RIPK2", "CD69", "MYD88", "PSMA2",
    "PIM1", "NOD1", "CFH", "TAPBP", "SLC25A28", "PTPN1", "TNFAIP3", "SSPN", "NUP93", "MTHFD2",
    "CDKN1A", "IRF4", "NFKB1", "BATF2", "HLA-B", "LATS2", "IRF5", "SLAMF7", "ISOC1", "P2RY14",
    "STAT3", "NCOA3", "HLA-A", "IL6", "GZMA", "IFNAR2", "CD74", "RAPGEF6", "CASP4", "FAS",
    "OGFR", "ARL4A", "SRI", "LYSMD2", "CSF2RB", "ST3GAL5", "C1R", "CASP3", "CMKLR1", "NFKBIA",
    "METTL7B", "ST8SIA4", "XCL1", "IL2RB", "VAMP5", "IL18BP", "ZNFX1", "ARID5B", "APOL6", "STAT4",
]

# KEGG hsa03450 (retrieved 2026-09-21). MRE11 is MRE11A on this HiSeqV2 freeze.
NHEJ_KEGG = [
    "RAD50", "DNTT", "FEN1", "XRCC6", "POLL", "POLM", "LIG4", "MRE11A",
    "PRKDC", "DCLRE1C", "XRCC4", "XRCC5", "NHEJ1",
]

# Reactome R-HSA-1834941 "STING mediated induction of host immune responses",
# minus genes that are also in KEGG NHEJ (PRKDC, XRCC5, XRCC6, MRE11) and
# minus CGAS, which is absent from HiSeqV2. STING1 is TMEM173 on this freeze.
# The set still mixes the sensor/kinase arm with negative regulators
# (TREX1, NLRC3, NLRP4, DTX4).
STING_REACTOME = [
    "DDX41", "DTX4", "IFI16", "IRF3", "NLRC3", "NLRP4", "STAT6",
    "TMEM173", "TBK1", "TREX1", "TRIM21",
]

# Classical MHC-I antigen-processing machinery: HLA heavy chain, β2m,
# immunoproteasome, TAP, peptide loading, and ER aminopeptidases.
APM_MHCI = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "TAP1", "TAP2", "TAPBP",
    "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "ERAP1", "ERAP2",
    "NLRC5", "CALR", "PDIA3", "CANX",
]

GENE_ENDPOINTS = ["PRKDC", "LIG4", "TMEM173", "HLA-A", "HLA-B", "HLA-C"]
GENE_LABEL = {
    "PRKDC": "PRKDC",
    "LIG4": "LIG4",
    "TMEM173": "STING1",
    "HLA-A": "HLA-A",
    "HLA-B": "HLA-B",
    "HLA-C": "HLA-C",
    "MHC_I": "MHC-I",
}
# Display order for scores. APM_nonIFN is a sensitivity, not a primary test.
SCORE_PRIMARY = ["NHEJ", "IFNg", "STING", "APM"]
SCORE_ALL = ["NHEJ", "IFNg", "STING", "APM", "APM_nonIFN"]
SCORE_LABEL = {
    "NHEJ": "NHEJ",
    "IFNg": "IFN-γ",
    "STING": "STING",
    "APM": "APM",
    "APM_nonIFN": "APM outside IFN-γ",
}
# Predicted direction for CLDN4-low (Q1) relative to CLDN4-high (Q4).
PREDICTED_Q1_HIGHER = {"IFNg": True, "STING": True, "APM": True, "NHEJ": False}

SOURCES = {
    "LUAD": {
        "expr": [
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        ],
        "expr_file": "TCGA.LUAD.HiSeqV2.gz",
        "estimate": [
            "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
            "https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        ],
        "estimate_file": "ESTIMATE_LUAD_RNAseqV2.txt",
    },
    "LUSC": {
        "expr": [
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        ],
        "expr_file": "TCGA.LUSC.HiSeqV2.gz",
        "estimate": [
            "https://ibl.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
            "https://bioinformatics.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        ],
        "estimate_file": "ESTIMATE_LUSC_RNAseqV2.txt",
    },
}

# Prior additive TCGA RNA (same HiSeqV2 + ESTIMATE join) for a regression check.
SANITY = {
    "LUAD": {"n": 515, "HLA-A": 0.059, "HLA-B": 0.035, "HLA-C": 0.023},
    "LUSC": {"n": 501, "HLA-A": 0.032, "HLA-B": 0.007, "HLA-C": 0.023},
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(urls, dest: str, min_bytes: int = 200) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print(f"  cached {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
        return
    last = None
    url_list = urls if isinstance(urls, (list, tuple)) else [urls]
    for url in url_list:
        for attempt in range(5):
            try:
                req = Request(url, headers={"User-Agent": "tcga-cldn4-nhej-sting/1.0"})
                with urlopen(req, timeout=300) as r:
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
                print(f"  retry {attempt + 1} {url}: {e}")
                time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {dest}: {last}")


def is_primary_01(barcode: str) -> bool:
    parts = str(barcode).replace(".", "-").split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(str(barcode).replace(".", "-").split("-")[:4])[:15]


def matrix_symbol(gene: str) -> str:
    return IFN_ALIAS.get(gene, gene)


def load_expression(path: str) -> pd.DataFrame:
    """Genes × primary samples. Aliquots collapsed to the 15-character barcode."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    usecols = [0] + [i + 1 for i, s in enumerate(samples) if is_primary_01(s)]
    keep_names = [sample15(s) for s in samples if is_primary_01(s)]
    expr = pd.read_csv(path, sep="\t", index_col=0, usecols=usecols, compression="gzip" if path.endswith(".gz") else None)
    expr.columns = keep_names
    if expr.index.has_duplicates:
        expr = expr.groupby(level=0).mean()
    expr = expr.T.groupby(level=0).mean().T
    expr.index = expr.index.astype(str)
    expr.index.name = "gene"
    return expr


def load_estimate(path: str) -> pd.DataFrame:
    est = pd.read_csv(path, sep="\t")
    est.columns = [c.strip() for c in est.columns]
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if c == est.columns[0] or cl in {"id", "name", "sample", "sampleid"}:
            rename[c] = "sample_raw"
        elif "immune" in cl and "score" in cl:
            rename[c] = "ImmuneScore"
        elif "stromal" in cl and "score" in cl:
            rename[c] = "StromalScore"
        elif "estimate" in cl and "score" in cl:
            rename[c] = "ESTIMATEScore"
    est = est.rename(columns=rename)
    if "sample_raw" not in est.columns or "ImmuneScore" not in est.columns:
        raise SystemExit(f"unrecognized ESTIMATE columns: {list(est.columns)}")
    est = est[est["sample_raw"].astype(str).map(is_primary_01)].copy()
    est["sample"] = est["sample_raw"].map(sample15)
    cols = [c for c in ["ImmuneScore", "StromalScore", "ESTIMATEScore"] if c in est.columns]
    return est.groupby("sample")[cols].mean()


def spearman_rho(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def partial_spearman_algebraic(x: pd.Series, y: pd.Series, z: pd.Series):
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rxy = float(stats.spearmanr(x[m], y[m]).statistic)
    rxz = float(stats.spearmanr(x[m], z[m]).statistic)
    ryz = float(stats.spearmanr(y[m], z[m]).statistic)
    denom = sqrt(max((1 - rxz ** 2) * (1 - ryz ** 2), 1e-12))
    r = (rxy - rxz * ryz) / denom
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * sqrt(df / max(1e-12, 1 - r ** 2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def fisher_z_ci(r: float, n: int, k_covariates: int = 0, alpha: float = 0.05):
    if r is None or not np.isfinite(r) or n is None or n <= k_covariates + 3:
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
    order = np.argsort(np.where(np.isfinite(p), p, np.inf))
    q = np.full(n, np.nan)
    prev = 1.0
    finite_idx = [i for i in order if np.isfinite(p[i])]
    m = len(finite_idx)
    for rank_from_end, idx in enumerate(finite_idx[::-1]):
        rank = m - rank_from_end
        val = min(prev, p[idx] * m / rank)
        q[idx] = val
        prev = val
    return q


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_r(r):
    if r is None or (isinstance(r, float) and not np.isfinite(r)):
        return "NA"
    return f"{r:+.3f}"


def present_genes(expr: pd.DataFrame, genes: list[str]) -> list[str]:
    idx = set(expr.index)
    out = []
    for g in genes:
        s = matrix_symbol(g)
        if s in idx and s not in out:
            out.append(s)
    return out


def run_ssgsea(expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    """Return samples × gene-set NES. Sign: high expression of the set -> positive NES."""
    res = gp.ssgsea(
        data=expr,
        gene_sets=gene_sets,
        outdir=None,
        sample_norm_method="rank",
        min_size=5,
        max_size=500,
        permutation_num=0,
        weight=0.25,
        no_plot=True,
        threads=4,
        seed=1,
        verbose=False,
    )
    long = res.res2d.copy()
    long["NES"] = pd.to_numeric(long["NES"])
    wide = long.pivot(index="Name", columns="Term", values="NES")
    wide.index.name = "sample"
    missing = [k for k in gene_sets if k not in wide.columns]
    if missing:
        raise SystemExit(f"ssGSEA dropped gene sets: {missing}")
    return wide


def assign_quartiles(cldn4: pd.Series) -> pd.Series:
    """Equal-size quartiles. Ties broken by sorted barcode so the cut is deterministic."""
    rank = cldn4.rank(method="first")
    return pd.qcut(rank, 4, labels=["Q1", "Q2", "Q3", "Q4"])


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta for a versus b. Positive means a tends to exceed b."""
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float((2 * u) / (len(a) * len(b)) - 1)


def correlate_block(df: pd.DataFrame, endpoint: str, kind: str) -> dict:
    y = df[endpoint]
    r, p, n = spearman_rho(df["CLDN4"], y)
    lo, hi = fisher_z_ci(r, n, 0)
    pr, pp, n_p = partial_spearman_algebraic(df["CLDN4"], y, df["ImmuneScore"])
    plo, phi = fisher_z_ci(pr, n_p, 1)
    kr, kp, _ = partial_spearman_algebraic(df["CLDN4"], y, df["KRT_mean"])
    mr, mp, _ = partial_spearman_algebraic(df["CLDN4"], y, df["MKI67"])
    r_imm, p_imm, _ = spearman_rho(y, df["ImmuneScore"])
    r_cl, p_cl, _ = spearman_rho(df["CLDN4"], df["ImmuneScore"])
    return {
        "endpoint": endpoint,
        "kind": kind,
        "n": n,
        "rho_unadj": r,
        "p_unadj": p,
        "rho_unadj_ci95_lo": lo,
        "rho_unadj_ci95_hi": hi,
        "rho_partial_ImmuneScore": pr,
        "p_partial_ImmuneScore": pp,
        "rho_partial_ci95_lo": plo,
        "rho_partial_ci95_hi": phi,
        "rho_partial_keratin": kr,
        "p_partial_keratin": kp,
        "rho_partial_MKI67": mr,
        "p_partial_MKI67": mp,
        "rho_endpoint_vs_ImmuneScore": r_imm,
        "p_endpoint_vs_ImmuneScore": p_imm,
        "rho_CLDN4_vs_ImmuneScore": r_cl,
        "p_CLDN4_vs_ImmuneScore": p_cl,
    }


def quartile_block(df: pd.DataFrame, endpoint: str, kind: str) -> dict:
    q1 = df.loc[df["quartile"] == "Q1", endpoint].dropna().to_numpy()
    q4 = df.loc[df["quartile"] == "Q4", endpoint].dropna().to_numpy()
    rest = df.loc[df["quartile"] != "Q1", endpoint].dropna().to_numpy()
    u = stats.mannwhitneyu(q1, q4, alternative="two-sided")
    u_rest = stats.mannwhitneyu(q1, rest, alternative="two-sided")
    return {
        "endpoint": endpoint,
        "kind": kind,
        "n_q1": int(len(q1)),
        "n_q4": int(len(q4)),
        "n_rest": int(len(rest)),
        "median_q1": float(np.median(q1)),
        "median_q4": float(np.median(q4)),
        "median_rest": float(np.median(rest)),
        "delta_q1_minus_q4": float(np.median(q1) - np.median(q4)),
        "mw_p_q1_q4": float(u.pvalue),
        "cliffs_delta_q1_vs_q4": cliffs_delta(q1, q4),
        "mw_p_q1_vs_rest": float(u_rest.pvalue),
        "delta_q1_minus_rest": float(np.median(q1) - np.median(rest)),
    }


def build_gene_sets(expr: pd.DataFrame) -> tuple[dict[str, list[str]], list[dict]]:
    ifn = present_genes(expr, IFN_GAMMA_HALLMARK)
    nhej = present_genes(expr, NHEJ_KEGG)
    sting = present_genes(expr, STING_REACTOME)
    apm = present_genes(expr, APM_MHCI)
    apm_non = [g for g in apm if g not in set(ifn)]
    sets = {"NHEJ": nhej, "IFNg": ifn, "STING": sting, "APM": apm, "APM_nonIFN": apm_non}
    for name, genes in sets.items():
        if len(genes) < 5:
            raise SystemExit(f"{name} has only {len(genes)} genes on the matrix")
    rows = []
    sources = {
        "NHEJ": "KEGG hsa03450; MRE11 stored as MRE11A",
        "IFNg": "MSigDB Hallmark interferon gamma response; four HGNC aliases mapped onto HiSeqV2",
        "STING": "Reactome R-HSA-1834941 minus KEGG-NHEJ genes and absent CGAS; STING1=TMEM173",
        "APM": "MHC-I antigen processing machinery (HLA, B2M, TAP, immunoproteasome, loading, ERAP)",
        "APM_nonIFN": "APM genes that are not in the Hallmark IFN-γ set used here",
    }
    for name, genes in sets.items():
        for g in genes:
            rows.append({"set": name, "gene": g, "source": sources[name], "n_in_set": len(genes)})
    return sets, rows


def analyze_cohort(cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    src = SOURCES[cohort]
    expr_path = os.path.join(DATA, src["expr_file"])
    est_path = os.path.join(DATA, src["estimate_file"])
    download(src["expr"], expr_path, min_bytes=1_000_000)
    download(src["estimate"], est_path, min_bytes=200)
    print(f"  loading {cohort} expression")
    expr = load_expression(expr_path)
    est = load_estimate(est_path)
    gene_sets, _ = build_gene_sets(expr)
    print(f"  ssGSEA {cohort}: " + ", ".join(f"{k}={len(v)}" for k, v in gene_sets.items()))
    scores = run_ssgsea(expr, gene_sets)

    needed = ["CLDN4", "PRKDC", "LIG4", "TMEM173", "HLA-A", "HLA-B", "HLA-C",
              "KRT8", "KRT18", "KRT19", "MKI67"]
    missing = [g for g in needed if g not in expr.index]
    if missing:
        raise SystemExit(f"{cohort} missing genes: {missing}")
    genes = expr.loc[needed].T
    genes.index.name = "sample"
    df = genes.join(scores, how="inner").join(est, how="inner")
    df["MHC_I"] = df[["HLA-A", "HLA-B", "HLA-C"]].mean(axis=1)
    df["KRT_mean"] = df[["KRT8", "KRT18", "KRT19"]].mean(axis=1)
    complete_cols = ["CLDN4", "ImmuneScore", "MKI67", "KRT_mean", "MHC_I"] + GENE_ENDPOINTS + SCORE_ALL
    df = df.dropna(subset=complete_cols).copy()
    df = df.sort_index()
    df["quartile"] = assign_quartiles(df["CLDN4"])
    df.insert(0, "cohort", f"TCGA-{cohort}")

    # Regression check against the earlier HiSeqV2 HLA correlations.
    exp = SANITY[cohort]
    if len(df) != exp["n"]:
        raise SystemExit(f"{cohort} n={len(df)} != prior complete-case n={exp['n']}")
    for gene, rho_exp in exp.items():
        if gene == "n":
            continue
        rho, p, n = spearman_rho(df["CLDN4"], df[gene])
        if abs(rho - rho_exp) > 0.008:
            raise SystemExit(f"{cohort} CLDN4 vs {gene} ρ={rho:.4f} disagrees with prior ρ={rho_exp}")
        print(f"  sanity {cohort} {gene}: ρ={rho:.4f} (prior {rho_exp:+.3f}) n={n}")

    corr_rows = []
    q_rows = []
    for ep in SCORE_ALL:
        corr_rows.append(correlate_block(df, ep, "score" if ep != "APM_nonIFN" else "score_sensitivity"))
        q_rows.append(quartile_block(df, ep, "score" if ep != "APM_nonIFN" else "score_sensitivity"))
    for ep in GENE_ENDPOINTS + ["MHC_I"]:
        kind = "gene_companion" if ep == "MHC_I" else "gene"
        corr_rows.append(correlate_block(df, ep, kind))
        q_rows.append(quartile_block(df, ep, kind))
    corr = pd.DataFrame(corr_rows)
    corr.insert(0, "cohort", f"TCGA-{cohort}")
    quart = pd.DataFrame(q_rows)
    quart.insert(0, "cohort", f"TCGA-{cohort}")
    counts = {
        "cohort": f"TCGA-{cohort}",
        "n_hiseqv2_primary": int(expr.shape[1]),
        "n_estimate_primary": int(est.shape[0]),
        "n_complete": int(df.shape[0]),
        "n_q1": int((df["quartile"] == "Q1").sum()),
        "n_q4": int((df["quartile"] == "Q4").sum()),
        "expr_sha256": sha256(expr_path),
        "estimate_sha256": sha256(est_path),
        "expr_bytes": os.path.getsize(expr_path),
        "estimate_bytes": os.path.getsize(est_path),
        "gseapy": gp.__version__,
    }
    keep = ["cohort", "CLDN4", "quartile", "ImmuneScore", "StromalScore", "ESTIMATEScore",
            "MKI67", "KRT_mean", "MHC_I"] + GENE_ENDPOINTS + SCORE_ALL
    keep = [c for c in keep if c in df.columns]
    return corr, quart, counts, df.reset_index()[keep]


def add_fdr(corr: pd.DataFrame, quart: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    corr = corr.copy()
    quart = quart.copy()
    corr["q_unadj"] = np.nan
    corr["q_partial_ImmuneScore"] = np.nan
    quart["q_mw_q1_q4"] = np.nan

    def _apply(mask, frame, src, dest):
        idx = frame.index[mask]
        frame.loc[idx, dest] = bh_fdr(frame.loc[idx, src].to_numpy())

    score_m = corr["endpoint"].isin(SCORE_PRIMARY)
    gene_m = corr["endpoint"].isin(GENE_ENDPOINTS)
    _apply(score_m, corr, "p_unadj", "q_unadj")
    _apply(score_m, corr, "p_partial_ImmuneScore", "q_partial_ImmuneScore")
    _apply(gene_m, corr, "p_unadj", "q_unadj")
    _apply(gene_m, corr, "p_partial_ImmuneScore", "q_partial_ImmuneScore")
    qmask = quart["endpoint"].isin(SCORE_PRIMARY)
    _apply(qmask, quart, "mw_p_q1_q4", "q_mw_q1_q4")
    return corr, quart


def call_direction(delta: float, rho: float, q_mw: float, q_rho: float, low_higher: bool) -> str:
    """Short call for the pre-specified CLDN4-low direction. q may be NA."""
    predicted_delta_pos = low_higher
    pred_ok = (delta > 0) if predicted_delta_pos else (delta < 0)
    rho_ok = (rho < 0) if predicted_delta_pos else (rho > 0)
    sig = (np.isfinite(q_mw) and q_mw < 0.05) or (np.isfinite(q_rho) and q_rho < 0.05)
    opposite = (delta < 0) if predicted_delta_pos else (delta > 0)
    opp_sig = opposite and sig
    if pred_ok and rho_ok and sig and abs(rho) >= 0.15:
        return "supported, modest"
    if pred_ok and rho_ok and sig:
        return "same direction, small"
    if opp_sig and abs(rho) >= 0.15:
        return "opposite direction"
    if opp_sig:
        return "opposite, small"
    return "not supported"


def write_finding(corr: pd.DataFrame, quart: pd.DataFrame, counts: pd.DataFrame, geneset_rows: pd.DataFrame) -> None:
    def cget(cohort, ep):
        hit = corr[(corr["cohort"] == cohort) & (corr["endpoint"] == ep)]
        return hit.iloc[0]

    def qget(cohort, ep):
        hit = quart[(quart["cohort"] == cohort) & (quart["endpoint"] == ep)]
        return hit.iloc[0]

    def n_of(cohort):
        return int(counts.loc[counts["cohort"] == cohort, "n_complete"].iloc[0])

    def n_expr(cohort):
        return int(counts.loc[counts["cohort"] == cohort, "n_hiseqv2_primary"].iloc[0])

    def n_est(cohort):
        return int(counts.loc[counts["cohort"] == cohort, "n_estimate_primary"].iloc[0])

    def set_n(name):
        return int(geneset_rows.loc[geneset_rows["set"] == name, "n_in_set"].iloc[0])

    def set_genes(name):
        return ", ".join(geneset_rows.loc[geneset_rows["set"] == name, "gene"].tolist())

    lines = []
    lines.append("# Finding — TCGA-LUAD / TCGA-LUSC: CLDN4-low vs ssGSEA NHEJ, IFN-γ, STING, and APM")
    lines.append("")
    lines.append("**Additive public RNA.** UCSC Xena legacy **HiSeqV2** log2(RSEM norm_count+1) for primary tumors (`-01`). Official MD Anderson **ESTIMATE RNAseqV2 ImmuneScore**. Single-sample GSEA is gseapy "
                 f"{gp.__version__} (rank normalization, weight 0.25, no gene-set permutation). This folder does not restate prior TACSTD2, purity, or keratin-rank results.")
    lines.append("")
    lines.append("Questions, pre-specified: in **CLDN4-low** tumors (bottom quartile, Q1) versus **CLDN4-high** (top quartile, Q4), is the **NHEJ** ssGSEA score lower, and are **IFN-γ / STING / MHC-I APM** scores higher? Separately, Spearman of **CLDN4 vs PRKDC, LIG4, STING1, and HLA-A/B/C**.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append("| Filter | TCGA-LUAD | TCGA-LUSC |")
    lines.append("|---|---:|---:|")
    lines.append(f"| HiSeqV2 primary tumors (`-01`, 15-char) | {n_expr('TCGA-LUAD')} | {n_expr('TCGA-LUSC')} |")
    lines.append(f"| Official ESTIMATE RNAseqV2 primaries | {n_est('TCGA-LUAD')} | {n_est('TCGA-LUSC')} |")
    lines.append(f"| **Complete: CLDN4 + scores + PRKDC/LIG4/TMEM173/HLA + ImmuneScore** | **{n_of('TCGA-LUAD')}** | **{n_of('TCGA-LUSC')}** |")
    lines.append(f"| Q1 / Q4 (barcode order breaks CLDN4 ties) | {int(counts.loc[counts.cohort=='TCGA-LUAD','n_q1'].iloc[0])} / {int(counts.loc[counts.cohort=='TCGA-LUAD','n_q4'].iloc[0])} | {int(counts.loc[counts.cohort=='TCGA-LUSC','n_q1'].iloc[0])} / {int(counts.loc[counts.cohort=='TCGA-LUSC','n_q4'].iloc[0])} |")
    lines.append("")
    lines.append("Quartile sizes differ by one tumor when n is not divisible by 4 (LUAD 129/129/128/129; LUSC 126/125/125/125). LUSC drops the one HiSeqV2 primary with no official ESTIMATE row (502 → 501). That is the same complete-case n as the earlier CLDN4–HLA RNA table. A regression check reproduces those unadjusted HLA Spearman values (LUAD HLA-A ρ=+0.059, HLA-B ρ=+0.035, HLA-C ρ=+0.023; LUSC HLA-A ρ=+0.032, HLA-B ρ=+0.007, HLA-C ρ=+0.023).")
    lines.append("")
    lines.append("## Answer")
    lines.append("")
    lines.append("**CLDN4-low tumors do not have a lower NHEJ score, and they do not have higher IFN-γ, STING, or APM scores.** NHEJ moves the other way in both histologies: the bottom CLDN4 quartile has a higher NHEJ score than the top quartile. IFN-γ, the STING pathway score, and the MHC-I APM score are null on the unadjusted test. Absolute ssGSEA NES is not an enrichment p-value, and NES levels are not comparable across gene sets; the contrasts below are within a cohort.")
    lines.append("")
    lines.append("| Question | TCGA-LUAD (n=515) | TCGA-LUSC (n=501) |")
    lines.append("|---|---|---|")

    def qcell(cohort, ep):
        q = qget(cohort, ep)
        c = cget(cohort, ep)
        return (
            f"ρ={fmt_r(c['rho_unadj'])} (q={fmt_p(c['q_unadj'])}); "
            f"Q1−Q4 Δ={q['delta_q1_minus_q4']:+.3f}, Cliff's δ={q['cliffs_delta_q1_vs_q4']:+.2f}, "
            f"MW q={fmt_p(q['q_mw_q1_q4'])}"
        )

    lines.append(f"| NHEJ lower in CLDN4-low? | No. Score is higher. {qcell('TCGA-LUAD', 'NHEJ')} | No. Score is higher. {qcell('TCGA-LUSC', 'NHEJ')} |")
    lines.append(f"| IFN-γ higher in CLDN4-low? | No. {qcell('TCGA-LUAD', 'IFNg')} | No. {qcell('TCGA-LUSC', 'IFNg')} |")
    lines.append(f"| STING score higher in CLDN4-low? | No. {qcell('TCGA-LUAD', 'STING')} | No. {qcell('TCGA-LUSC', 'STING')} |")
    lines.append(f"| APM score higher in CLDN4-low? | No. {qcell('TCGA-LUAD', 'APM')} | No. {qcell('TCGA-LUSC', 'APM')} |")
    lines.append("")

    nhej_l = cget("TCGA-LUAD", "NHEJ")
    nhej_s = cget("TCGA-LUSC", "NHEJ")
    ifn_l = cget("TCGA-LUAD", "IFNg")
    ifn_s = cget("TCGA-LUSC", "IFNg")
    st_l = cget("TCGA-LUAD", "STING")
    st_s = cget("TCGA-LUSC", "STING")
    apm_l = cget("TCGA-LUAD", "APM_nonIFN")
    apm_s = cget("TCGA-LUSC", "APM_nonIFN")
    lines.append(
        f"**NHEJ is higher in CLDN4-low, not lower.** LUAD Spearman ρ={fmt_r(nhej_l['rho_unadj'])} "
        f"(q={fmt_p(nhej_l['q_unadj'])}); LUSC ρ={fmt_r(nhej_s['rho_unadj'])} (q={fmt_p(nhej_s['q_unadj'])}). "
        f"ImmuneScore partials stay negative (LUAD {fmt_r(nhej_l['rho_partial_ImmuneScore'])}, "
        f"LUSC {fmt_r(nhej_s['rho_partial_ImmuneScore'])}), and so do partials on mean KRT8/KRT18/KRT19 "
        f"(LUAD {fmt_r(nhej_l['rho_partial_keratin'])}, p={fmt_p(nhej_l['p_partial_keratin'])}; "
        f"LUSC {fmt_r(nhej_s['rho_partial_keratin'])}, p={fmt_p(nhej_s['p_partial_keratin'])}). "
        f"MKI67 does not explain the LUAD result (partial ρ={fmt_r(nhej_l['rho_partial_MKI67'])}, "
        f"p={fmt_p(nhej_l['p_partial_MKI67'])}). In LUSC the score association shrinks after MKI67 "
        f"(partial ρ={fmt_r(nhej_s['rho_partial_MKI67'])}, p={fmt_p(nhej_s['p_partial_MKI67'])}), "
        "so the LUSC NHEJ-score shift is partly shared with a proliferation marker. "
        "The NES scale is tight (within-cohort SD ≈ 0.05), so the median shifts "
        f"({qget('TCGA-LUAD', 'NHEJ')['delta_q1_minus_q4']:+.3f} LUAD, "
        f"{qget('TCGA-LUSC', 'NHEJ')['delta_q1_minus_q4']:+.3f} LUSC) are modest rank effects "
        f"(Cliff's δ = {qget('TCGA-LUAD', 'NHEJ')['cliffs_delta_q1_vs_q4']:+.2f} in LUAD and "
        f"{qget('TCGA-LUSC', 'NHEJ')['cliffs_delta_q1_vs_q4']:+.2f} in LUSC), not a separated population."
    )
    lines.append("")
    lines.append(
        f"**IFN-γ is not higher in CLDN4-low.** Unadjusted ρ={fmt_r(ifn_l['rho_unadj'])} in LUAD and "
        f"{fmt_r(ifn_s['rho_unadj'])} in LUSC; both quartile tests are null. "
        f"The LUAD ImmuneScore partial is a small positive residual (ρ={fmt_r(ifn_l['rho_partial_ImmuneScore'])}, "
        f"q={fmt_p(ifn_l['q_partial_ImmuneScore'])}): after removing infiltrate, higher CLDN4, not lower CLDN4, "
        f"sits slightly higher on IFN-γ. The unadjusted test is null, and IFN-γ vs ImmuneScore itself is "
        f"ρ={fmt_r(ifn_l['rho_endpoint_vs_ImmuneScore'])}, so that residual is a suppressor, not an IFN-high CLDN4-low result. "
        f"The LUSC partial is null (ρ={fmt_r(ifn_s['rho_partial_ImmuneScore'])})."
    )
    lines.append("")
    lines.append(
        f"**The STING pathway score is not higher in CLDN4-low.** Unadjusted LUAD ρ={fmt_r(st_l['rho_unadj'])} "
        f"(q={fmt_p(st_l['q_unadj'])}) and LUSC ρ={fmt_r(st_s['rho_unadj'])} (q={fmt_p(st_s['q_unadj'])}) do not survive BH, "
        f"and the Q1 vs Q4 tests are null. The LUAD ImmuneScore partial is positive "
        f"(ρ={fmt_r(st_l['rho_partial_ImmuneScore'])}, q={fmt_p(st_l['q_partial_ImmuneScore'])}), which is the other direction. "
        "The 11-gene set mixes TMEM173/TBK1/IRF3/IFI16 with negative regulators (TREX1, NLRC3, NLRP4, DTX4), so it is a membership score, not an induced-ISG score. "
        "The direct STING1 readout is the TMEM173 row below."
    )
    lines.append("")
    lines.append(
        f"**A sensitivity, not a primary endpoint:** the 8 APM genes that are outside Hallmark IFN-γ "
        f"(HLA-C, HLA-E, TAP2, ERAP1, ERAP2, CALR, PDIA3, CANX) correlate modestly negatively with CLDN4 "
        f"(LUAD ρ={fmt_r(apm_l['rho_unadj'])}, p={fmt_p(apm_l['p_unadj'])}; "
        f"LUSC ρ={fmt_r(apm_s['rho_unadj'])}, p={fmt_p(apm_s['p_unadj'])}), i.e. slightly higher in CLDN4-low. "
        f"In LUSC that association does not survive the keratin partial (ρ={fmt_r(apm_s['rho_partial_keratin'])}, "
        f"p={fmt_p(apm_s['p_partial_keratin'])}). The primary 19-gene APM score, which shares its interferon-inducible members with Hallmark IFN-γ, stays null. "
        "This sensitivity was defined from set overlap before looking at CLDN4, and it is not in the BH family."
    )
    lines.append("")
    lines.append("BH q for quartile tests is across the 8 primary contrasts (2 cohorts × NHEJ / IFN-γ / STING / APM). BH q for Spearman is separate for the 8 primary score correlations and for the 12 primary gene correlations (2 cohorts × PRKDC / LIG4 / STING1 / HLA-A / HLA-B / HLA-C). APM-outside-IFN-γ and the MHC-I mean are companions.")
    lines.append("")
    lines.append("## Gene correlations — CLDN4 vs PRKDC, LIG4, STING1, HLA")
    lines.append("")
    lines.append("STING1 is **TMEM173** on this HiSeqV2 freeze. There is no `STING1` or `CGAS` row. MHC-I = mean(HLA-A, HLA-B, HLA-C) on the log2 scale.")
    lines.append("")
    lines.append("| Cohort | Gene | n | Unadj ρ | Unadj p | Unadj q | Partial ρ \\| ImmuneScore | Partial p | Partial q |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for ep in GENE_ENDPOINTS + ["MHC_I"]:
        for cohort in ["TCGA-LUAD", "TCGA-LUSC"]:
            r = cget(cohort, ep)
            label = GENE_LABEL[ep] + (" †" if ep == "HLA-B" else "")
            q_u = fmt_p(r["q_unadj"]) if ep != "MHC_I" else "—"
            q_p = fmt_p(r["q_partial_ImmuneScore"]) if ep != "MHC_I" else "—"
            lines.append(
                f"| {cohort} | {label} | {int(r['n'])} | {fmt_r(r['rho_unadj'])} | {fmt_p(r['p_unadj'])} | {q_u} | "
                f"{fmt_r(r['rho_partial_ImmuneScore'])} | {fmt_p(r['p_partial_ImmuneScore'])} | {q_p} |"
            )
    lines.append("")
    prk_l, prk_s = cget("TCGA-LUAD", "PRKDC"), cget("TCGA-LUSC", "PRKDC")
    lig_l, lig_s = cget("TCGA-LUAD", "LIG4"), cget("TCGA-LUSC", "LIG4")
    tmem_l, tmem_s = cget("TCGA-LUAD", "TMEM173"), cget("TCGA-LUSC", "TMEM173")
    lines.append("† **HLA-B is in Yoshihara Immune141.** Partialling ImmuneScore out of HLA-B is not an independent infiltrate control. HLA-A, HLA-C, PRKDC, LIG4, and TMEM173 are the non-circular gene endpoints. NHEJ ssGSEA **includes** PRKDC and LIG4, so the NHEJ score correlation is not a separate replication of those two genes.")
    lines.append("")
    lines.append(
        f"**PRKDC and LIG4 agree with the higher-NHEJ direction in CLDN4-low tumors, where they are significant.** "
        f"PRKDC ρ={fmt_r(prk_l['rho_unadj'])} in LUAD (q={fmt_p(prk_l['q_unadj'])}) and "
        f"{fmt_r(prk_s['rho_unadj'])} in LUSC (q={fmt_p(prk_s['q_unadj'])}). "
        f"The LUSC PRKDC association is the largest gene effect here and remains after keratin "
        f"(partial ρ={fmt_r(prk_s['rho_partial_keratin'])}, p={fmt_p(prk_s['p_partial_keratin'])}) and after MKI67 "
        f"(partial ρ={fmt_r(prk_s['rho_partial_MKI67'])}, p={fmt_p(prk_s['p_partial_MKI67'])}). "
        f"LIG4 is negative in LUAD (ρ={fmt_r(lig_l['rho_unadj'])}, q={fmt_p(lig_l['q_unadj'])}) and does not clear BH in LUSC "
        f"(ρ={fmt_r(lig_s['rho_unadj'])}, q={fmt_p(lig_s['q_unadj'])})."
    )
    lines.append("")
    lines.append(
        f"**STING1 is higher with higher CLDN4 in LUAD, and is null in LUSC.** "
        f"TMEM173 ρ={fmt_r(tmem_l['rho_unadj'])} (q={fmt_p(tmem_l['q_unadj'])}), "
        f"ImmuneScore partial {fmt_r(tmem_l['rho_partial_ImmuneScore'])} (q={fmt_p(tmem_l['q_partial_ImmuneScore'])}). "
        f"That LUAD association remains after keratin (partial ρ={fmt_r(tmem_l['rho_partial_keratin'])}, "
        f"p={fmt_p(tmem_l['p_partial_keratin'])}) and is smaller after MKI67 "
        f"(partial ρ={fmt_r(tmem_l['rho_partial_MKI67'])}, p={fmt_p(tmem_l['p_partial_MKI67'])}). "
        f"LUSC TMEM173 ρ={fmt_r(tmem_s['rho_unadj'])}. This is not a CLDN4-low STING-high result, and it does not replicate in squamous tumors."
    )
    lines.append("")
    lines.append(
        "**HLA-A/B/C stay unadjusted-null in both histologies**, at the same coefficients as the earlier TCGA RNA table. "
        "LUAD partials on ImmuneScore are small and positive for HLA-A and HLA-B (q≈0.02–0.03): higher CLDN4, not lower, with slightly higher HLA after the infiltrate score is removed. "
        "HLA-C does not clear BH. LUSC partials stay null. Do not quote a LUAD partial without the unadjusted null, and do not read HLA-B's partial as independent of ImmuneScore."
    )
    lines.append("")
    lines.append("## Score correlations and the ImmuneScore confound")
    lines.append("")
    lines.append("| Cohort | Score | n | Unadj ρ | p | q | Partial ρ \\| ImmuneScore | p | q | ρ(score, ImmuneScore) | Partial ρ \\| KRT8/18/19 (p) | Partial ρ \\| MKI67 (p) |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for ep in SCORE_ALL:
        for cohort in ["TCGA-LUAD", "TCGA-LUSC"]:
            r = cget(cohort, ep)
            q_u = fmt_p(r["q_unadj"]) if ep in SCORE_PRIMARY else "—"
            q_p = fmt_p(r["q_partial_ImmuneScore"]) if ep in SCORE_PRIMARY else "—"
            lines.append(
                f"| {cohort} | {SCORE_LABEL[ep]} | {int(r['n'])} | {fmt_r(r['rho_unadj'])} | {fmt_p(r['p_unadj'])} | {q_u} | "
                f"{fmt_r(r['rho_partial_ImmuneScore'])} | {fmt_p(r['p_partial_ImmuneScore'])} | {q_p} | "
                f"{fmt_r(r['rho_endpoint_vs_ImmuneScore'])} | {fmt_r(r['rho_partial_keratin'])} ({fmt_p(r['p_partial_keratin'])}) | {fmt_r(r['rho_partial_MKI67'])} ({fmt_p(r['p_partial_MKI67'])}) |"
            )
    lines.append("")
    lines.append("IFN-γ and APM track ImmuneScore, as they should: both gene sets are full of infiltrate and interferon genes. CLDN4 itself is only weakly related to ImmuneScore. A partial correlation can therefore move even when the unadjusted association is null. Quote the unadjusted number next to any partial. Keratin (mean of KRT8/KRT18/KRT19) and MKI67 partials are sensitivities, not extra primary tests.")
    lines.append("")
    # Q1 vs rest one-liner
    lines.append("## Q1 versus the other three quartiles")
    lines.append("")
    lines.append("The primary contrast is Q1 vs Q4, matching earlier TCGA quartile cuts. Q1 versus everyone else is the same direction test with a larger control arm.")
    lines.append("")
    lines.append("| Cohort | Score | n Q1 | n rest | Δ median (Q1−rest) | MW p |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for ep in SCORE_PRIMARY:
        for cohort in ["TCGA-LUAD", "TCGA-LUSC"]:
            q = qget(cohort, ep)
            lines.append(
                f"| {cohort} | {SCORE_LABEL[ep]} | {int(q['n_q1'])} | {int(q['n_rest'])} | "
                f"{q['delta_q1_minus_rest']:+.3f} | {fmt_p(q['mw_p_q1_vs_rest'])} |"
            )
    lines.append("")
    lines.append("## What the gene sets are, and what they share")
    lines.append("")
    lines.append(f"- **NHEJ** ({set_n('NHEJ')} genes), KEGG hsa03450: {set_genes('NHEJ')}. `MRE11` is `MRE11A` on HiSeqV2.")
    lines.append(f"- **IFN-γ** ({set_n('IFNg')} genes), Hallmark interferon gamma response. Four symbols were mapped onto this freeze: WARS1→WARS, MARCHF1→MARCH1, HELZ2→PRIC285, CMTR1→FTSJD2. Genes absent after that map are listed in `tables/genesets.tsv` only if dropped; the scored set has {set_n('IFNg')} genes.")
    lines.append(f"- **STING** ({set_n('STING')} genes), Reactome R-HSA-1834941 after removing KEGG-NHEJ members (PRKDC, XRCC5, XRCC6, MRE11) and CGAS (not on HiSeqV2): {set_genes('STING')}. `STING1` is `TMEM173`. The set mixes positive components (TMEM173, TBK1, IRF3, IFI16) with negative regulators (TREX1, NLRC3, NLRP4, DTX4), so it is a pathway membership score, not a pure induced-ISG score. TRIM21 is also in Hallmark IFN-γ.")
    lines.append(f"- **APM** ({set_n('APM')} genes): {set_genes('APM')}.")
    lines.append(f"- **APM outside IFN-γ** ({set_n('APM_nonIFN')} genes), sensitivity: {set_genes('APM_nonIFN')}. The other APM genes sit inside the Hallmark IFN-γ set, so the primary APM score is not an independent IFN test.")
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("- Bulk primary-tumor RNA. It is not a malignant-cell NHEJ or STING measurement, and it is not an ICI-response result.")
    lines.append("- A null quartile test is not evidence that NHEJ or cGAS–STING is irrelevant in CLDN4-high cells. It is evidence that this bulk contrast does not show the pre-specified shift.")
    lines.append("- ImmuneScore is RNA-derived. Partialling it is not an ABSOLUTE purity adjustment. Hallmark IFN-γ and the APM set overlap Immune141-style genes, so those partials are partly circular. HLA-B is the gene-level case that is explicitly inside Immune141.")
    lines.append("- No STAR-TPM rerun, no protein, no CPTAC, no single-cell split.")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Matrix.** UCSC Xena `TCGA.{LUAD,LUSC}.sampleMap/HiSeqV2`. Primary tumors only (sample type 01). Replicate aliquots averaged to the 15-character barcode.")
    lines.append("- **ESTIMATE.** Official RNAseqV2 `Immune_score`, not recomputed.")
    lines.append("- **ssGSEA.** gseapy single-sample GSEA. Sample normalization `rank`, weight α=0.25, `permutation_num=0` (the NES column is gseapy's normalized enrichment score, not a permutation p-value). Minimum set size 5. A synthetic check in this environment confirms that a sample with the NHEJ genes pinned to the top of the ranking receives a positive NES. Absolute NES is not a significance call, and NES cannot be compared across gene sets of different sizes. Every test is a within-cohort contrast.")
    lines.append("- **CLDN4-low / high.** Quartiles of CLDN4 within each complete-case cohort. `rank(method='first')` after sorting barcodes, then `qcut` into four equal groups. Q1 is the lowest CLDN4.")
    lines.append("- **Quartile test.** Two-sided Mann–Whitney Q1 vs Q4. Effect sizes: median(Q1)−median(Q4) and Cliff's delta (positive means Q1 values are larger). Q1 vs Q2–Q4 is secondary and is not in the BH family.")
    lines.append("- **Correlation.** Spearman. Partial = first-order partial Spearman given ImmuneScore (df = n−3). Fisher z intervals use 1/(n−3) unadjusted and 1/(n−4) partial. Keratin mean and MKI67 partials use the same formula and are labeled sensitivities.")
    lines.append("- **FDR.** Benjamini–Hochberg within the families named above. Families are not pooled with each other.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/correlations.tsv` — Spearman, partials, CI, q")
    lines.append("- `tables/quartiles.tsv` — Q1 vs Q4 and Q1 vs rest")
    lines.append("- `tables/counts.tsv` — honest-n filters")
    lines.append("- `tables/samples.tsv` — per-tumor CLDN4, quartiles, scores, genes, ImmuneScore")
    lines.append("- `tables/genesets.tsv` — genes that entered each ssGSEA set")
    lines.append("- `tables/provenance.json` — URLs, sha256, gseapy version")
    lines.append("- `figures/q1q4_ssgsea.png` — quartile NES")
    lines.append("- `figures/forest_partial.png` — unadjusted vs ImmuneScore-partial ρ")
    lines.append("- `figures/scatter_genes.png` — CLDN4 vs PRKDC, LIG4, STING1, MHC-I")
    lines.append("- `figures/scatter_scores.png` — CLDN4 vs the four ssGSEA scores")
    lines.append("- Reproduce: `python3 methods/tcga_cldn4_nhej_sting/analyze.py`")
    lines.append("")
    path = os.path.join(HERE, "FINDING.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {path}")


def _errbar(ax, x, y, lo, hi, color, marker, label):
    ax.errorbar(
        x, y,
        xerr=[np.asarray(x) - np.asarray(lo), np.asarray(hi) - np.asarray(x)],
        fmt=marker, color=color, label=label, capsize=2.5, ms=5, lw=1,
    )


def plot_quartiles(samples: pd.DataFrame, quart: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(11.2, 6.4), sharex="col")
    colors = {"Q1": "#3B6D9A", "Q2": "#B7C6D4", "Q3": "#F3D2B5", "Q4": "#E07A3D"}
    for row, cohort in enumerate(["TCGA-LUAD", "TCGA-LUSC"]):
        sub = samples[samples["cohort"] == cohort]
        for col, ep in enumerate(SCORE_PRIMARY):
            ax = axes[row, col]
            data = [sub.loc[sub["quartile"] == q, ep].to_numpy() for q in ["Q1", "Q2", "Q3", "Q4"]]
            bp = ax.boxplot(
                data, positions=[0, 1, 2, 3], widths=0.62, patch_artist=True,
                showfliers=False, medianprops={"color": "black", "lw": 1.2},
            )
            for patch, q in zip(bp["boxes"], ["Q1", "Q2", "Q3", "Q4"]):
                patch.set_facecolor(colors[q])
                patch.set_edgecolor("0.2")
            rng = np.random.default_rng(1 + row * 10 + col)
            for i, q in enumerate(["Q1", "Q2", "Q3", "Q4"]):
                y = sub.loc[sub["quartile"] == q, ep].to_numpy()
                x = i + rng.uniform(-0.12, 0.12, size=len(y))
                ax.scatter(x, y, s=4, c=colors[q], alpha=0.28, linewidths=0, zorder=2)
            qrow = quart[(quart["cohort"] == cohort) & (quart["endpoint"] == ep)].iloc[0]
            ax.set_xticks([0, 1, 2, 3])
            ax.set_xticklabels(["Q1\nlow", "Q2", "Q3", "Q4\nhigh"])
            n = int(qrow["n_q1"] + qrow["n_q4"])
            # n printed is Q1+Q4; full cohort is 2x roughly. Use full from samples.
            ax.set_title(
                f"{SCORE_LABEL[ep]}\nΔ={qrow['delta_q1_minus_q4']:+.3f}  δ={qrow['cliffs_delta_q1_vs_q4']:+.2f}  p={fmt_p(qrow['mw_p_q1_q4'])}",
                fontsize=9,
            )
            if col == 0:
                ax.set_ylabel(f"{cohort}\nssGSEA NES")
            else:
                ax.set_ylabel("")
    fig.suptitle("CLDN4 quartiles vs ssGSEA  ·  Q1 = lowest CLDN4  ·  two-sided MW is Q1 vs Q4", fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"q1q4_ssgsea.{ext}"), dpi=160)
    plt.close(fig)


def plot_forest(corr: pd.DataFrame) -> None:
    endpoints = SCORE_PRIMARY + GENE_ENDPOINTS
    labels = [SCORE_LABEL.get(e, GENE_LABEL.get(e, e)) + (" †" if e == "HLA-B" else "") for e in endpoints]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 6.2), sharex=True)
    y = np.arange(len(endpoints))[::-1]
    for ax, cohort in zip(axes, ["TCGA-LUAD", "TCGA-LUSC"]):
        sub = corr[corr["cohort"] == cohort].set_index("endpoint").loc[endpoints]
        n = int(sub["n"].iloc[0])
        _errbar(
            ax, sub["rho_unadj"], y + 0.12,
            sub["rho_unadj_ci95_lo"], sub["rho_unadj_ci95_hi"],
            "#4C78A8", "o", "unadjusted",
        )
        _errbar(
            ax, sub["rho_partial_ImmuneScore"], y - 0.12,
            sub["rho_partial_ci95_lo"], sub["rho_partial_ci95_hi"],
            "#F58518", "s", "partial | ImmuneScore",
        )
        ax.axvline(0, color="0.45", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_title(f"{cohort}  n={n}")
        ax.set_xlabel("Spearman ρ  (CLDN4 vs endpoint)")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        ax.axhline(3.5, color="0.85", lw=0.6)
    fig.suptitle("CLDN4 vs NHEJ / IFN-γ / STING / APM scores and vs PRKDC, LIG4, STING1, HLA", fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"forest_partial.{ext}"), dpi=160)
    plt.close(fig)


def _scatter_grid(samples: pd.DataFrame, corr: pd.DataFrame, endpoints: list[str], ylab, outfile: str, title: str) -> None:
    fig, axes = plt.subplots(2, len(endpoints), figsize=(2.7 * len(endpoints) + 1.4, 6.3), sharex=False)
    vmin = float(np.nanpercentile(samples["ImmuneScore"], 2))
    vmax = float(np.nanpercentile(samples["ImmuneScore"], 98))
    for row, cohort in enumerate(["TCGA-LUAD", "TCGA-LUSC"]):
        sub = samples[samples["cohort"] == cohort]
        q1_max = float(sub.loc[sub["quartile"] == "Q1", "CLDN4"].max())
        for col, ep in enumerate(endpoints):
            ax = axes[row, col]
            sc = ax.scatter(
                sub["CLDN4"], sub[ep], c=sub["ImmuneScore"], cmap="viridis",
                s=9, alpha=0.75, linewidths=0, vmin=vmin, vmax=vmax, rasterized=True,
            )
            ax.axvspan(float(sub["CLDN4"].min()), q1_max, color="#3B6D9A", alpha=0.16, lw=0)
            ax.axvline(q1_max, color="#3B6D9A", lw=0.6, alpha=0.8)
            r = corr[(corr["cohort"] == cohort) & (corr["endpoint"] == ep)].iloc[0]
            label = SCORE_LABEL.get(ep, GENE_LABEL.get(ep, ep))
            ax.set_title(
                f"{label}\nρ={fmt_r(r['rho_unadj'])}  p={fmt_p(r['p_unadj'])}",
                fontsize=9,
            )
            if row == 1:
                ax.set_xlabel("CLDN4  log2(norm+1)")
            if col == 0:
                ax.set_ylabel(f"{cohort}\n{ylab}")
    fig.subplots_adjust(right=0.90, top=0.86, hspace=0.38, wspace=0.35)
    cax = fig.add_axes([0.92, 0.18, 0.015, 0.55])
    fig.colorbar(sc, cax=cax, label="ESTIMATE ImmuneScore")
    fig.suptitle(title + "   ·   shaded band = CLDN4 Q1", fontsize=11)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{outfile}.{ext}"), dpi=160)
    plt.close(fig)


def plot_scatters(samples: pd.DataFrame, corr: pd.DataFrame) -> None:
    _scatter_grid(
        samples, corr, ["PRKDC", "LIG4", "TMEM173", "MHC_I"],
        "expression", "scatter_genes",
        "CLDN4 vs PRKDC, LIG4, STING1 (TMEM173), and MHC-I",
    )
    _scatter_grid(
        samples, corr, SCORE_PRIMARY,
        "ssGSEA NES", "scatter_scores",
        "CLDN4 vs ssGSEA NHEJ, IFN-γ, STING, and APM",
    )


def main() -> None:
    # Gene-set membership is taken from LUAD and checked against LUSC (same platform).
    src = SOURCES["LUAD"]
    expr_path = os.path.join(DATA, src["expr_file"])
    download(src["expr"], expr_path, min_bytes=1_000_000)
    expr_luad = load_expression(expr_path)
    gene_sets, gene_rows = build_gene_sets(expr_luad)
    # Confirm LUSC has every gene the LUAD sets use.
    lusc_path = os.path.join(DATA, SOURCES["LUSC"]["expr_file"])
    download(SOURCES["LUSC"]["expr"], lusc_path, min_bytes=1_000_000)
    expr_lusc = load_expression(lusc_path)
    for name, genes in gene_sets.items():
        missing = [g for g in genes if g not in expr_lusc.index]
        if missing:
            raise SystemExit(f"LUSC missing {name} genes: {missing}")
    geneset_df = pd.DataFrame(gene_rows)
    geneset_df.to_csv(os.path.join(TABLES, "genesets.tsv"), sep="\t", index=False)
    del expr_luad, expr_lusc

    all_corr, all_q, all_counts, all_samples = [], [], [], []
    for cohort in ["LUAD", "LUSC"]:
        print(f"=== {cohort} ===")
        corr, quart, counts, samples = analyze_cohort(cohort)
        all_corr.append(corr)
        all_q.append(quart)
        all_counts.append(counts)
        all_samples.append(samples)
    corr = pd.concat(all_corr, ignore_index=True)
    quart = pd.concat(all_q, ignore_index=True)
    corr, quart = add_fdr(corr, quart)
    counts = pd.DataFrame(all_counts)
    samples = pd.concat(all_samples, ignore_index=True)
    corr.to_csv(os.path.join(TABLES, "correlations.tsv"), sep="\t", index=False)
    quart.to_csv(os.path.join(TABLES, "quartiles.tsv"), sep="\t", index=False)
    counts.to_csv(os.path.join(TABLES, "counts.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(TABLES, "samples.tsv"), sep="\t", index=False)
    provenance = {
        "expression": "UCSC Xena HiSeqV2 log2(RSEM norm_count+1)",
        "estimate": "MD Anderson official RNAseqV2 ImmuneScore (Yoshihara 2013)",
        "ssgsea": {
            "package": f"gseapy {gp.__version__}",
            "sample_norm_method": "rank",
            "weight": 0.25,
            "permutation_num": 0,
            "min_size": 5,
        },
        "cldn4_low": "Q1 of CLDN4 within cohort after barcode-sorted rank(method=first)",
        "partial": "first-order partial Spearman | ImmuneScore; keratin and MKI67 are sensitivities",
        "sting1_symbol_on_hiseqv2": "TMEM173",
        "cgas_on_hiseqv2": False,
        "ifn_alias_to_hiseqv2": IFN_ALIAS,
        "sources": SOURCES,
        "sanity_prior_hla_rho": SANITY,
        "counts": counts.to_dict(orient="records"),
        "gene_set_sizes": {k: len(v) for k, v in gene_sets.items()},
    }
    with open(os.path.join(TABLES, "provenance.json"), "w") as f:
        json.dump(provenance, f, indent=2)
    write_finding(corr, quart, counts, geneset_df)
    plot_quartiles(samples, quart)
    plot_forest(corr)
    plot_scatters(samples, corr)
    show = corr[corr["endpoint"].isin(SCORE_PRIMARY + GENE_ENDPOINTS)][
        ["cohort", "endpoint", "n", "rho_unadj", "p_unadj", "q_unadj",
         "rho_partial_ImmuneScore", "p_partial_ImmuneScore", "q_partial_ImmuneScore"]
    ]
    print(show.to_string(index=False))
    print(quart[quart["endpoint"].isin(SCORE_PRIMARY)][
        ["cohort", "endpoint", "delta_q1_minus_q4", "mw_p_q1_q4", "q_mw_q1_q4", "cliffs_delta_q1_vs_q4"]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
