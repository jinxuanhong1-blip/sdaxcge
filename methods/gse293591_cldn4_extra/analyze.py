#!/usr/bin/env python3
"""GSE293591 leftover extra: CLDN4 vs CD274 and IFN beyond MHC-I.

Public BostonGene RNA-seq / IHC biomarker series (Kushnarev et al., Sci Rep 2025,
PMID 40753302; GEO GSE293591). MHC-I / APM on this matrix is treated as already
known (partial) and is scored only as a reference row. This leftover is CD274
and IFN genes/signatures after MHC-I genes are removed.

Honest n is the number of columns in the deposited TPM table that also have a
GEO diagnosis label. IHC PD-L1, ICI response, and FASTQ are not on GEO.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE293591_CLDN4_DATA", "/tmp/gse293591_cldn4"))
SEED = 20260817
N_BOOT = 2000

TPM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293591/"
    "suppl/GSE293591_TPM_all_samples.tsv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293591/"
    "matrix/GSE293591_series_matrix.txt.gz"
)

# Locked custom classical MHC-I / APM (same 21 genes as methods/gse285029_cldn4_gsea).
MHC_I_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
    "B2M", "TAP1", "TAP2", "TAPBP", "TAPBPL", "NLRC5",
    "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3", "IRF1",
]

# Compact MHC-I core (the "already known" partial slice).
MHC_I_CORE = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP"]

# Hallmark IFN-γ (MSigDB Hallmark 2020 / Liberzon Cell Syst 2015).
HALLMARK_IFNG = [
    "ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2",
    "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4",
    "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40",
    "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1",
    "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58",
    "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A",
    "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2",
    "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1",
    "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L",
    "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2",
    "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7",
    "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15",
    "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2",
    "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2",
    "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA",
    "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR",
    "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1",
    "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3",
    "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2",
    "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2",
    "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1",
    "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2",
    "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1",
    "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2",
    "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14",
    "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18",
    "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1",
]

# Ayers 6-gene IFN-γ (Ayers et al., JCI 2017). HLA-DRA is MHC-II, not MHC-I.
AYERS6 = ["IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"]

# Compact ISG / IFN panel with no classical MHC-I / APM gene.
ISG_NO_MHC = [
    "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1",
    "IFNG", "GBP1", "ISG15", "MX1", "IFIT1", "OAS2",
]

# Single-gene leftover panel (IFN beyond MHC-I + CD274).
IFN_GENES = [
    "CD274", "IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11",
    "IDO1", "ISG15", "MX1", "IFIT1", "OAS2", "GBP1",
]

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

LUNG_DX = {"Lung Adenocarcinoma", "Squamous Cell Carcinoma of Lung"}
GI_DX = {
    "Colorectal Adenocarcinoma",
    "Gastric Adenocarcinoma",
    "Esophageal Adenocarcinoma",
    "Pancreatic Adenocarcinoma",
    "Liver carcinoma and Cholangiocarcinoma",
}
BREAST_DX = {"Breast cancer"}

MHC_I_SET = set(MHC_I_APM)


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"download {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key = line.split("\t", 1)[0][1:]
            vals = next(csv.reader([line], delimiter="\t"))[1:]
            if key == "Sample_characteristics_ch1":
                if vals and vals[0].startswith("diagnosis:"):
                    fields["diagnosis"] = [v.replace("diagnosis: ", "") for v in vals]
                elif vals and vals[0].startswith("tissue preservation:"):
                    fields["preservation"] = [
                        v.replace("tissue preservation: ", "") for v in vals
                    ]
            else:
                fields[key] = vals
    n = len(fields["Sample_geo_accession"])
    sample_id = []
    for title, desc in zip(fields["Sample_title"], fields["Sample_description"]):
        # titles: "Sample_001, Breast cancer, biopsy"
        sid = title.split(",")[0].strip()
        if not sid.startswith("Sample_"):
            sid = desc.replace("Library name: ", "").strip()
        sample_id.append(sid)
    out = pd.DataFrame(
        {
            "sample_id": sample_id,
            "gsm": fields["Sample_geo_accession"],
            "title": fields["Sample_title"],
            "diagnosis": fields["diagnosis"],
            "preservation": fields["preservation"],
            "molecule": fields["Sample_molecule_ch1"],
        }
    )
    if len(out) != n:
        raise RuntimeError(f"series matrix parse length mismatch: {len(out)} vs {n}")
    if out["sample_id"].duplicated().any():
        raise RuntimeError("duplicate sample_id in series matrix")
    return out


def log2p1(x: pd.Series | pd.DataFrame) -> pd.Series | pd.DataFrame:
    return np.log2(np.clip(x, 0, None) + 1.0)


def zscore_mean(log_mat: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in log_mat.index]
    if not present:
        return pd.Series(np.nan, index=log_mat.columns), []
    z = log_mat.loc[present].apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1)
    return z.mean(axis=0), present


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def bootstrap_spearman_ci(
    x: pd.Series, y: pd.Series, n_boot: int = N_BOOT, seed: int = SEED
) -> tuple[float, float]:
    m = x.notna() & y.notna()
    xv = x[m].to_numpy()
    yv = y[m].to_numpy()
    n = len(xv)
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rhos[i] = stats.spearmanr(xv[idx], yv[idx]).statistic
    return float(np.nanpercentile(rhos, 2.5)), float(np.nanpercentile(rhos, 97.5))


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> tuple[float, float, int]:
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return float("nan"), float("nan"), n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    bx = np.polyfit(zr, xr, 1)
    by = np.polyfit(zr, yr, 1)
    rx = xr - (bx[0] * zr + bx[1])
    ry = yr - (by[0] * zr + by[1])
    rho, p = stats.spearmanr(rx, ry)
    return float(rho), float(p), n


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ: P(a>b) - P(a<b). Positive = a larger than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    # vectorized pairwise
    diff = a[:, None] - b[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def high_low_cut(x: pd.Series, y: pd.Series, how: str) -> dict:
    m = x.notna() & y.notna()
    xv = x[m]
    yv = y[m]
    if how == "median":
        lo_cut = hi_cut = float(xv.median())
        lo = yv[xv <= lo_cut]
        hi = yv[xv > hi_cut]
        # if ties pile on the median, split by rank so both sides are non-empty
        if lo.empty or hi.empty:
            r = xv.rank(method="first")
            mid = r.median()
            lo = yv[r <= mid]
            hi = yv[r > mid]
    elif how == "tertile":
        q1, q2 = xv.quantile([1 / 3, 2 / 3])
        lo = yv[xv <= q1]
        hi = yv[xv >= q2]
        lo_cut, hi_cut = float(q1), float(q2)
    elif how == "quartile":
        q1, q3 = xv.quantile([0.25, 0.75])
        lo = yv[xv <= q1]
        hi = yv[xv >= q3]
        lo_cut, hi_cut = float(q1), float(q3)
    else:
        raise ValueError(how)
    n_lo, n_hi = int(lo.size), int(hi.size)
    if n_lo < 3 or n_hi < 3:
        return {
            "cut": how,
            "n_low": n_lo,
            "n_high": n_hi,
            "cut_low": lo_cut,
            "cut_high": hi_cut,
            "mwu_p": float("nan"),
            "cliffs_delta": float("nan"),
            "median_low": float("nan"),
            "median_high": float("nan"),
        }
    u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    return {
        "cut": how,
        "n_low": n_lo,
        "n_high": n_hi,
        "cut_low": lo_cut,
        "cut_high": hi_cut,
        "mwu_p": float(p),
        "cliffs_delta": cliffs_delta(hi.to_numpy(), lo.to_numpy()),
        "median_low": float(lo.median()),
        "median_high": float(hi.median()),
        "U": float(u),
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        # rank from the end
        k = n - rank + 1
        val = p[i] * n / k
        prev = min(prev, val)
        q[i] = prev
    return [float(min(1.0, v)) for v in q]


def cohort_mask(meta: pd.DataFrame, name: str) -> pd.Series:
    if name == "all":
        return pd.Series(True, index=meta.index)
    if name == "FFPE":
        return meta["preservation"] == "FFPE"
    if name == "FF":
        return meta["preservation"] == "FF"
    if name == "breast":
        return meta["diagnosis"].isin(BREAST_DX)
    if name == "lung":
        return meta["diagnosis"].isin(LUNG_DX)
    if name == "LUAD":
        return meta["diagnosis"] == "Lung Adenocarcinoma"
    if name == "LUSC":
        return meta["diagnosis"] == "Squamous Cell Carcinoma of Lung"
    if name == "GI":
        return meta["diagnosis"].isin(GI_DX)
    raise KeyError(name)


def verdict(n: int, rho: float, p: float, *, leftover: bool) -> str:
    if n < 20 or not np.isfinite(rho):
        return "UNDERPOWERED"
    if leftover and n < 40:
        return "EXPLORATORY"
    if p >= 0.05:
        return "NO_EVIDENCE"
    if rho > 0:
        return "POSITIVE"
    return "NEGATIVE"


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    tpm_path = download(TPM_URL, DATA / "GSE293591_TPM_all_samples.tsv.gz")
    mtx_path = download(MATRIX_URL, DATA / "GSE293591_series_matrix.txt.gz")

    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    tpm.index = tpm.index.astype(str)
    if tpm.index.duplicated().any():
        tpm = tpm.groupby(level=0).mean()
    meta = parse_series_matrix(mtx_path).set_index("sample_id")

    shared = [c for c in tpm.columns if c in meta.index]
    missing_tpm = sorted(set(tpm.columns) - set(meta.index))
    missing_meta = sorted(set(meta.index) - set(tpm.columns))
    if missing_tpm or missing_meta:
        raise RuntimeError(
            f"sample ID mismatch TPM-only={missing_tpm[:5]} meta-only={missing_meta[:5]}"
        )
    tpm = tpm[shared]
    meta = meta.loc[shared].copy()
    log = log2p1(tpm)

    # signatures
    ifng_full_score, ifng_full_genes = zscore_mean(log, HALLMARK_IFNG)
    ifng_no_mhc_genes = [g for g in HALLMARK_IFNG if g not in MHC_I_SET]
    ifng_no_mhc_score, ifng_no_mhc_used = zscore_mean(log, ifng_no_mhc_genes)
    ifng_no_mhc_cd274_genes = [g for g in ifng_no_mhc_genes if g != "CD274"]
    ifng_no_mhc_cd274_score, ifng_no_mhc_cd274_used = zscore_mean(
        log, ifng_no_mhc_cd274_genes
    )
    ayers_score, ayers_used = zscore_mean(log, AYERS6)
    isg_score, isg_used = zscore_mean(log, ISG_NO_MHC)
    mhc_score, mhc_used = zscore_mean(log, MHC_I_APM)
    mhc_core_score, mhc_core_used = zscore_mean(log, MHC_I_CORE)
    epi_score, epi_used = zscore_mean(log, EPITHELIAL)

    scores = pd.DataFrame(
        {
            "CLDN4": log.loc["CLDN4"],
            "CD274": log.loc["CD274"],
            "TACSTD2": log.loc["TACSTD2"],
            "IFNG_gene": log.loc["IFNG"],
            "Hallmark_IFNg": ifng_full_score,
            "Hallmark_IFNg_noMHC": ifng_no_mhc_score,
            "Hallmark_IFNg_noMHC_noCD274": ifng_no_mhc_cd274_score,
            "Ayers6": ayers_score,
            "ISG_noMHC": isg_score,
            "MHC_I_APM": mhc_score,
            "MHC_I_core": mhc_core_score,
            "epithelial": epi_score,
        }
    )
    for g in IFN_GENES + MHC_I_CORE + ["TACSTD2"]:
        if g in log.index and g not in scores.columns:
            scores[g] = log.loc[g]
    scores = scores.join(meta)

    scores["group_lung"] = scores["diagnosis"].isin(LUNG_DX)
    scores["group_breast"] = scores["diagnosis"].isin(BREAST_DX)
    scores["group_GI"] = scores["diagnosis"].isin(GI_DX)

    # inventory
    inv_rows = [
        {"item": "GEO series samples / TPM columns", "n": int(tpm.shape[1]), "note": "deposited TPM"},
        {"item": "Genes on TPM matrix", "n": int(tpm.shape[0]), "note": "20062 protein-coding (author)"},
        {"item": "Unique GSM", "n": int(meta["gsm"].nunique()), "note": "all unique"},
        {"item": "Unique sample_id", "n": int(meta.index.nunique()), "note": "Sample_001..Sample_377 with gaps"},
        {"item": "FFPE / total RNA / target enrichment", "n": int((meta["preservation"] == "FFPE").sum()), "note": "Procrustes batch-corrected"},
        {"item": "FF / polyA", "n": int((meta["preservation"] == "FF").sum()), "note": "Cureline fresh-frozen"},
        {"item": "Breast cancer", "n": int((meta["diagnosis"] == "Breast cancer").sum()), "note": "GEO diagnosis"},
        {"item": "Lung Adenocarcinoma", "n": int((meta["diagnosis"] == "Lung Adenocarcinoma").sum()), "note": "GEO diagnosis"},
        {"item": "Squamous Cell Carcinoma of Lung", "n": int((meta["diagnosis"] == "Squamous Cell Carcinoma of Lung").sum()), "note": "GEO diagnosis"},
        {"item": "Lung (LUAD+LUSC)", "n": int(meta["diagnosis"].isin(LUNG_DX).sum()), "note": "does not include non-lung SCC"},
        {"item": "GI (CRC+gastric+eso+panc+liver/CCA)", "n": int(meta["diagnosis"].isin(GI_DX).sum()), "note": "locked GI bag"},
        {"item": "CLDN4 finite", "n": int(log.loc["CLDN4"].notna().sum()), "note": "all samples"},
        {"item": "CD274 finite", "n": int(log.loc["CD274"].notna().sum()), "note": "all samples"},
        {"item": "PD-L1 IHC on GEO", "n": 0, "note": "paper has IHC; not deposited"},
        {"item": "ICI / response labels on GEO", "n": 0, "note": "biomarker series, not an ICI trial"},
        {"item": "FASTQ / SRA on GEO", "n": 0, "note": "author withheld identifiable sequence"},
    ]
    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    dx_counts = (
        meta["diagnosis"].value_counts().rename_axis("diagnosis").reset_index(name="n")
    )
    dx_counts.to_csv(TABLES / "diagnosis_counts.tsv", sep="\t", index=False)

    # coverage
    cov_rows = []
    for name, genes, used in [
        ("Hallmark_IFNg", HALLMARK_IFNG, ifng_full_genes),
        ("Hallmark_IFNg_noMHC", ifng_no_mhc_genes, ifng_no_mhc_used),
        ("Hallmark_IFNg_noMHC_noCD274", ifng_no_mhc_cd274_genes, ifng_no_mhc_cd274_used),
        ("Ayers6", AYERS6, ayers_used),
        ("ISG_noMHC", ISG_NO_MHC, isg_used),
        ("MHC_I_APM", MHC_I_APM, mhc_used),
        ("MHC_I_core", MHC_I_CORE, mhc_core_used),
        ("epithelial", EPITHELIAL, epi_used),
    ]:
        cov_rows.append(
            {
                "signature": name,
                "n_listed": len(genes),
                "n_present": len(used),
                "missing": ",".join([g for g in genes if g not in log.index]) or "",
            }
        )
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    endpoints = [
        ("CD274", "CD274", True, "leftover: PD-L1 transcript"),
        ("Hallmark_IFNg_noMHC", "Hallmark_IFNg_noMHC", True, "leftover: Hallmark IFN-γ minus MHC-I/APM"),
        ("Hallmark_IFNg_noMHC_noCD274", "Hallmark_IFNg_noMHC_noCD274", True, "leftover: IFN-γ minus MHC-I and CD274"),
        ("Ayers6", "Ayers6", True, "leftover: Ayers 6-gene IFN-γ"),
        ("ISG_noMHC", "ISG_noMHC", True, "leftover: compact ISG, no MHC-I"),
        ("IFNG_gene", "IFNG_gene", True, "leftover: IFNG transcript"),
        ("MHC_I_APM", "MHC_I_APM", False, "already known (partial): 21-gene MHC-I/APM"),
        ("MHC_I_core", "MHC_I_core", False, "already known (partial): HLA-A/B/C + B2M + TAP1/2 + TAPBP"),
        ("Hallmark_IFNg", "Hallmark_IFNg", False, "reference: full Hallmark IFN-γ (includes MHC-I + CD274)"),
        ("TACSTD2", "TACSTD2", False, "companion only"),
        ("epithelial", "epithelial", False, "companion: epithelial mean-z"),
    ]

    cohorts = ["all", "FFPE", "breast", "lung", "GI"]

    spearman_rows = []
    cut_rows = []
    for cohort in cohorts:
        mask = cohort_mask(scores, cohort)
        sub = scores.loc[mask]
        n_coh = int(mask.sum())
        for ep, col, leftover, note in endpoints:
            rho, p, n = spearman(sub["CLDN4"], sub[col])
            lo, hi = bootstrap_spearman_ci(sub["CLDN4"], sub[col])
            rho_adj, p_adj, n_adj = partial_spearman(
                sub["CLDN4"], sub[col], sub["epithelial"]
            )
            spearman_rows.append(
                {
                    "cohort": cohort,
                    "n_cohort": n_coh,
                    "endpoint": ep,
                    "leftover": leftover,
                    "n": n,
                    "rho": rho,
                    "p": p,
                    "ci95_lo": lo,
                    "ci95_hi": hi,
                    "rho_adj_epithelial": rho_adj,
                    "p_adj_epithelial": p_adj,
                    "n_adj": n_adj,
                    "verdict": verdict(n, rho, p, leftover=leftover),
                    "note": note,
                }
            )
            for how in ("median", "tertile", "quartile"):
                cut = high_low_cut(sub["CLDN4"], sub[col], how)
                cut.update(
                    {
                        "cohort": cohort,
                        "n_cohort": n_coh,
                        "endpoint": ep,
                        "leftover": leftover,
                        "note": note,
                    }
                )
                cut_rows.append(cut)

    spearman_df = pd.DataFrame(spearman_rows)
    # BH within leftover endpoints, per cohort, continuous tests only
    spearman_df["q_bh_leftover"] = np.nan
    for cohort in cohorts:
        idx = spearman_df.index[
            (spearman_df["cohort"] == cohort) & (spearman_df["leftover"])
        ]
        spearman_df.loc[idx, "q_bh_leftover"] = bh_fdr(
            spearman_df.loc[idx, "p"].tolist()
        )
    spearman_df.to_csv(TABLES / "spearman_cldn4_vs_endpoints.tsv", sep="\t", index=False)

    cut_df = pd.DataFrame(cut_rows)
    cut_df.to_csv(TABLES / "extra_cuts_cldn4.tsv", sep="\t", index=False)

    # single-gene IFN leftover panel on all / lung
    gene_rows = []
    for cohort in ("all", "FFPE", "breast", "lung", "GI"):
        mask = cohort_mask(scores, cohort)
        sub = scores.loc[mask]
        for g in IFN_GENES + MHC_I_CORE:
            if g not in log.index:
                continue
            rho, p, n = spearman(sub["CLDN4"], log.loc[g, sub.index])
            rho_adj, p_adj, _ = partial_spearman(
                sub["CLDN4"], log.loc[g, sub.index], sub["epithelial"]
            )
            gene_rows.append(
                {
                    "cohort": cohort,
                    "gene": g,
                    "panel": "leftover_IFN_or_CD274" if g in IFN_GENES else "MHC_I_core_known",
                    "n": n,
                    "rho": rho,
                    "p": p,
                    "rho_adj_epithelial": rho_adj,
                    "p_adj_epithelial": p_adj,
                }
            )
    gene_df = pd.DataFrame(gene_rows)
    for cohort in gene_df["cohort"].unique():
        idx = gene_df.index[
            (gene_df["cohort"] == cohort) & (gene_df["panel"] == "leftover_IFN_or_CD274")
        ]
        gene_df.loc[idx, "q_bh"] = bh_fdr(gene_df.loc[idx, "p"].tolist())
    gene_df.to_csv(TABLES / "spearman_cldn4_vs_genes.tsv", sep="\t", index=False)

    # per-sample table (compact)
    keep_cols = [
        "gsm", "diagnosis", "preservation", "molecule",
        "CLDN4", "CD274", "TACSTD2", "IFNG_gene",
        "Hallmark_IFNg_noMHC", "Hallmark_IFNg_noMHC_noCD274",
        "Ayers6", "ISG_noMHC", "MHC_I_APM", "MHC_I_core", "epithelial",
    ]
    scores[keep_cols].to_csv(TABLES / "per_sample.tsv", sep="\t")

    # headline one-row table (all tumors, leftover + MHC-I reference)
    headline_eps = [
        "CD274",
        "Hallmark_IFNg_noMHC",
        "Hallmark_IFNg_noMHC_noCD274",
        "Ayers6",
        "ISG_noMHC",
        "MHC_I_APM",
    ]
    one = spearman_df[
        (spearman_df["cohort"] == "all") & (spearman_df["endpoint"].isin(headline_eps))
    ].copy()
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # summary json
    def pick(cohort: str, ep: str) -> dict:
        r = spearman_df[(spearman_df["cohort"] == cohort) & (spearman_df["endpoint"] == ep)].iloc[0]
        return {
            "n": int(r["n"]),
            "rho": float(r["rho"]),
            "p": float(r["p"]),
            "ci95": [float(r["ci95_lo"]), float(r["ci95_hi"])],
            "rho_adj": float(r["rho_adj_epithelial"]),
            "p_adj": float(r["p_adj_epithelial"]),
            "verdict": r["verdict"],
            "q_bh": (None if pd.isna(r["q_bh_leftover"]) else float(r["q_bh_leftover"])),
        }

    q4 = cut_df[(cut_df["cohort"] == "all") & (cut_df["cut"] == "quartile")]
    summary = {
        "dataset": "GSE293591",
        "pmid": "40753302",
        "n_tpm": int(tpm.shape[1]),
        "n_genes": int(tpm.shape[0]),
        "n_FFPE": int((meta["preservation"] == "FFPE").sum()),
        "n_FF": int((meta["preservation"] == "FF").sum()),
        "n_breast": int((meta["diagnosis"] == "Breast cancer").sum()),
        "n_lung": int(meta["diagnosis"].isin(LUNG_DX).sum()),
        "n_LUAD": int((meta["diagnosis"] == "Lung Adenocarcinoma").sum()),
        "n_LUSC": int((meta["diagnosis"] == "Squamous Cell Carcinoma of Lung").sum()),
        "n_GI": int(meta["diagnosis"].isin(GI_DX).sum()),
        "coverage": cov.to_dict(orient="records"),
        "all": {ep: pick("all", ep) for ep in headline_eps},
        "FFPE": {ep: pick("FFPE", ep) for ep in ["CD274", "Hallmark_IFNg_noMHC", "MHC_I_APM"]},
        "breast": {ep: pick("breast", ep) for ep in ["CD274", "Hallmark_IFNg_noMHC", "MHC_I_APM"]},
        "lung": {ep: pick("lung", ep) for ep in ["CD274", "Hallmark_IFNg_noMHC", "MHC_I_APM"]},
        "GI": {ep: pick("GI", ep) for ep in ["CD274", "Hallmark_IFNg_noMHC", "MHC_I_APM"]},
        "q4q1_all": q4.set_index("endpoint")[
            ["n_low", "n_high", "mwu_p", "cliffs_delta"]
        ].to_dict(orient="index"),
        "seed": SEED,
        "n_boot": N_BOOT,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # figures
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    def scatter(ax, x, y, title, xlabel, ylabel):
        ax.scatter(x, y, s=12, alpha=0.45, c="#1f4e79", edgecolors="none")
        rho, p, n = spearman(x, y)
        ax.set_title(f"{title}\nρ={rho:.3f} p={fmt_p(p)} n={n}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    scatter(
        axes[0],
        scores["CLDN4"],
        scores["CD274"],
        "All tumors",
        "CLDN4 log2(TPM+1)",
        "CD274 log2(TPM+1)",
    )
    lung = scores.loc[scores["diagnosis"].isin(LUNG_DX)]
    scatter(
        axes[1],
        lung["CLDN4"],
        lung["CD274"],
        "Lung (LUAD+LUSC)",
        "CLDN4 log2(TPM+1)",
        "CD274 log2(TPM+1)",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.png")
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.pdf")
    plt.close(fig)

    # forest leftover + MHC-I reference, all tumors
    forest_eps = [
        ("CD274", "CD274"),
        ("Hallmark_IFNg_noMHC", "IFN-γ no MHC-I"),
        ("Hallmark_IFNg_noMHC_noCD274", "IFN-γ no MHC-I / CD274"),
        ("Ayers6", "Ayers 6-gene"),
        ("ISG_noMHC", "ISG no MHC-I"),
        ("IFNG_gene", "IFNG"),
        ("MHC_I_APM", "MHC-I/APM (known)"),
        ("MHC_I_core", "MHC-I core (known)"),
    ]
    fsub = spearman_df[spearman_df["cohort"] == "all"].set_index("endpoint")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ys = np.arange(len(forest_eps))[::-1]
    for y, (ep, lab) in zip(ys, forest_eps):
        r = fsub.loc[ep]
        color = "#1f4e79" if r["leftover"] else "#7a7a7a"
        ax.plot([r["ci95_lo"], r["ci95_hi"]], [y, y], color=color, lw=1.6)
        ax.plot(r["rho"], y, "o", color=color, ms=6)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab in forest_eps])
    ax.set_xlabel("Spearman ρ vs CLDN4 (all tumors, n=365)")
    ax.set_title("GSE293591 leftover extra: CLDN4 vs CD274 / IFN beyond MHC-I")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_spearman_forest.png")
    fig.savefig(FIGURES / "fig2_spearman_forest.pdf")
    plt.close(fig)

    # Q4 vs Q1 violins
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.6))
    q1, q3 = scores["CLDN4"].quantile([0.25, 0.75])
    hi = scores["CLDN4"] >= q3
    lo = scores["CLDN4"] <= q1
    for ax, col, title in [
        (axes[0], "CD274", "CD274"),
        (axes[1], "Hallmark_IFNg_noMHC", "IFN-γ no MHC-I"),
        (axes[2], "MHC_I_APM", "MHC-I/APM (known)"),
    ]:
        data = [scores.loc[lo, col], scores.loc[hi, col]]
        parts = ax.violinplot(data, showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor("#1f4e79")
            pc.set_alpha(0.45)
        ax.scatter(
            np.repeat([1, 2], [lo.sum(), hi.sum()]),
            np.concatenate(data),
            s=6,
            alpha=0.25,
            c="0.2",
        )
        p = cut_df[
            (cut_df["cohort"] == "all")
            & (cut_df["cut"] == "quartile")
            & (cut_df["endpoint"] == col)
        ].iloc[0]
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"Q1\nn={int(lo.sum())}", f"Q4\nn={int(hi.sum())}"])
        ax.set_title(f"{title}\nδ={p['cliffs_delta']:+.2f} p={fmt_p(p['mwu_p'])}")
        ax.set_ylabel("score" if col != "CD274" else "log2(TPM+1)")
    fig.suptitle("CLDN4 Q4 vs Q1, all tumors", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_q4q1_violins.png")
    fig.savefig(FIGURES / "fig3_q4q1_violins.pdf")
    plt.close(fig)

    # subset forest for CD274 and IFN-no-MHC
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), sharex=True)
    for ax, ep, title in [
        (axes[0], "CD274", "CLDN4 vs CD274"),
        (axes[1], "Hallmark_IFNg_noMHC", "CLDN4 vs IFN-γ no MHC-I"),
    ]:
        rows = []
        for cohort in cohorts:
            r = spearman_df[
                (spearman_df["cohort"] == cohort) & (spearman_df["endpoint"] == ep)
            ].iloc[0]
            rows.append(r)
        ys = np.arange(len(rows))[::-1]
        for y, r in zip(ys, rows):
            ax.plot([r["ci95_lo"], r["ci95_hi"]], [y, y], color="#1f4e79", lw=1.6)
            ax.plot(r["rho"], y, "o", color="#1f4e79", ms=6)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r['cohort']} n={int(r['n'])}" for r in rows])
        ax.set_title(title)
        ax.set_xlabel("Spearman ρ")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_subset_forest.png")
    fig.savefig(FIGURES / "fig4_subset_forest.pdf")
    plt.close(fig)

    print("n_all", tpm.shape[1])
    print("diagnosis")
    print(dx_counts.to_string(index=False))
    print("\nheadline all:")
    print(
        spearman_df[spearman_df["cohort"] == "all"][
            ["endpoint", "n", "rho", "p", "verdict"]
        ].to_string(index=False)
    )
    print("wrote", TABLES)
    print("wrote", FIGURES)


if __name__ == "__main__":
    main()
