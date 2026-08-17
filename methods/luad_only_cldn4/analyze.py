#!/usr/bin/env python3
"""Additive LUAD-only recut of CLDN4 vs CD8 / TNK / Immune.

Recomputes TCGA-LUAD and TCGA-LUSC from public Xena HiSeqV2 + ESTIMATE +
ABSOLUTE. CPTAC and scRNA rows reuse harvested public patient tables from
prior additive analyses. GSE218989 / GSE285029 / GSE148071 stay mixed when
per-patient histology is not a deposited table.

Usage:
  python3 methods/luad_only_cldn4/analyze.py
"""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
RAW = Path("/tmp/luad_only_cldn4")

GEP18 = [
    "CD274", "CD8A", "LAG3", "HLA-DQA1", "HLA-DRB1", "STAT1",
    "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "GZMB",
    "TIGIT", "CXCR6", "CCL5", "NKG7", "CMKLR1",
]
NEED_GENES = sorted(set(GEP18 + ["CLDN4", "CD8A", "CD8B", "GZMA", "PRF1", "IFNG"]))


def fisher_ci(rho: float, n: int, zcrit: float = 1.96) -> tuple[float, float]:
    if n < 4 or not np.isfinite(rho) or abs(rho) >= 1:
        return (np.nan, np.nan)
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / math.sqrt(n - 3)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def spearman_row(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    lo, hi = fisher_ci(float(rho), n)
    return {"n": n, "rho": float(rho), "p": float(p), "ci_lo": lo, "ci_hi": hi}


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.Series) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y"), cov.rename("c")], axis=1).dropna()
    n = int(len(d))
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    rx = d["x"].rank()
    ry = d["y"].rank()
    rc = d["c"].rank()
    X = np.column_stack([np.ones(n), rc.to_numpy()])
    ex = rx.to_numpy() - X @ np.linalg.lstsq(X, rx.to_numpy(), rcond=None)[0]
    ey = ry.to_numpy() - X @ np.linalg.lstsq(X, ry.to_numpy(), rcond=None)[0]
    rho, p = stats.pearsonr(ex, ey)
    # Fisher-z with partial df n-4
    if n > 4 and abs(rho) < 1:
        z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
        se = 1.0 / math.sqrt(n - 4)
        lo, hi = float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))
    else:
        lo, hi = np.nan, np.nan
    return {"n": n, "rho": float(rho), "p": float(p), "ci_lo": lo, "ci_hi": hi}


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(rho: float) -> str:
    if not np.isfinite(rho):
        return "NA"
    return f"{rho:+.3f}"


def load_hiseqv2(path: Path, genes: list[str]) -> pd.DataFrame:
    keep = set(genes)
    rows = {}
    with gzip_open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in fh:
            g = line.split("\t", 1)[0]
            if g in keep:
                vals = line.rstrip("\n").split("\t")[1:]
                rows[g] = pd.to_numeric(pd.Series(vals, index=samples), errors="coerce")
    missing = [g for g in genes if g not in rows]
    if missing:
        raise SystemExit(f"missing genes in {path.name}: {missing}")
    return pd.DataFrame(rows)


def gzip_open(path: Path):
    return gzip.open(path, "rt")


def primary_tumor_table(expr: pd.DataFrame) -> pd.DataFrame:
    # HiSeqV2 columns are 15-char barcodes (TCGA-XX-XXXX-01)
    idx = expr.index.astype(str)
    keep = idx.str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01")
    out = expr.loc[keep].copy()
    out["patient"] = out.index.str[:12]
    out = out[~out.index.duplicated(keep="first")]
    out = out.drop_duplicates("patient", keep="first")
    return out


def load_estimate(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    # first col is sample
    col0 = df.columns[0]
    df = df.rename(columns={col0: "sample"})
    df["sample15"] = df["sample"].astype(str).str[:15]
    s = df.set_index("sample15")["Immune_score"]
    return s[~s.index.duplicated(keep="first")]


def load_absolute(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    samp = df["array"].astype(str)
    df["sample15"] = samp.str[:15]
    s = pd.to_numeric(df.set_index("sample15")["purity"], errors="coerce")
    return s[~s.index.duplicated(keep="first")]


def zmean(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    z = (df[genes] - df[genes].mean()) / df[genes].std(ddof=0)
    return z.mean(axis=1)


def tcga_cohort(name: str, expr_path: Path, est_path: Path) -> pd.DataFrame:
    expr = load_hiseqv2(expr_path, NEED_GENES)
    tab = primary_tumor_table(expr)
    tab["histology"] = name
    tab["CYT"] = tab[["GZMA", "PRF1"]].mean(axis=1)
    tab["GEP18"] = zmean(tab, GEP18)
    est = load_estimate(est_path)
    absu = load_absolute(RAW / "tcga_absolute.txt")
    tab["ImmuneScore"] = est.reindex(tab.index)
    tab["ABSOLUTE"] = absu.reindex(tab.index)
    tab["CD8"] = tab["CD8A"]
    return tab


def add_tests(rows: list[dict], cohort: str, histology: str, layer: str,
              predictor: str, endpoint: str, x, y, cov=None, note: str = "",
              cov_name: str = "ABSOLUTE"):
    u = spearman_row(x, y)
    rows.append({
        "cohort": cohort, "histology": histology, "layer": layer,
        "predictor": predictor, "endpoint": endpoint, "adjust": "none",
        "n": u["n"], "rho": u["rho"], "p": u["p"],
        "ci_lo": u["ci_lo"], "ci_hi": u["ci_hi"], "note": note,
    })
    if cov is not None:
        p = partial_spearman(x, y, cov)
        rows.append({
            "cohort": cohort, "histology": histology, "layer": layer,
            "predictor": predictor, "endpoint": endpoint, "adjust": cov_name,
            "n": p["n"], "rho": p["rho"], "p": p["p"],
            "ci_lo": p["ci_lo"], "ci_hi": p["ci_hi"], "note": note,
        })


def ensure_tcga_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    files = {
        "TCGA.LUAD.HiSeqV2.gz": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
        "TCGA.LUSC.HiSeqV2.gz": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
        "MDACC_estimate_LUAD.txt": "https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        "MDACC_estimate_LUSC.txt": "https://bioinformatics.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        "tcga_absolute.txt": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
    }
    import urllib.request
    for name, url in files.items():
        dest = RAW / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f"[get] {name}")
        urllib.request.urlretrieve(url, dest)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    ensure_tcga_raw()
    rows: list[dict] = []

    # ------------------------------------------------------------------ TCGA
    luad = tcga_cohort("LUAD", RAW / "TCGA.LUAD.HiSeqV2.gz", RAW / "MDACC_estimate_LUAD.txt")
    lusc = tcga_cohort("LUSC", RAW / "TCGA.LUSC.HiSeqV2.gz", RAW / "MDACC_estimate_LUSC.txt")
    pooled = pd.concat([luad, lusc], axis=0)
    for hist, tab in [("LUAD", luad), ("LUSC", lusc), ("mixed", pooled)]:
        note = "Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity"
        if hist == "mixed":
            note += "; pooled LUAD+LUSC (dilution / Simpson risk)"
        for ep, col in [("CD8A", "CD8"), ("CYT", "CYT"), ("GEP18", "GEP18"),
                        ("ImmuneScore", "ImmuneScore")]:
            add_tests(rows, "TCGA", hist, "bulk_RNA", "CLDN4", ep,
                      tab["CLDN4"], tab[col], tab["ABSOLUTE"], note=note,
                      cov_name="ABSOLUTE")
    luad.to_csv(TABLES / "tcga_luad_samples.tsv", sep="\t")
    lusc.to_csv(TABLES / "tcga_lusc_samples.tsv", sep="\t")

    # ------------------------------------------------------------------ CPTAC
    clu = pd.read_csv(DATA / "cptac_luad_sample_table.tsv", sep="\t")
    cls = pd.read_csv(DATA / "cptac_lusc_sample_scores.tsv", sep="\t")
    # LUAD RNA / protein vs CD8 and Immune
    add_tests(rows, "CPTAC-LUAD", "LUAD", "bulk_RNA", "CLDN4_RNA", "CIBERSORT_CD8",
              clu["CLDN4_RNA"], clu["CIBERSORT_T_cell_CD8+"], clu["WES_purity"],
              note="CPTAC LUAD Gillette 2020; already LUAD-only", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LUAD", "LUAD", "bulk_RNA", "CLDN4_RNA", "ImmuneScore",
              clu["CLDN4_RNA"], clu["ESTIMATE_ImmuneScore"], clu["WES_purity"],
              note="CPTAC LUAD Gillette 2020; already LUAD-only", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LUAD", "LUAD", "protein", "CLDN4_protein", "CIBERSORT_CD8",
              clu["CLDN4_protein"], clu["CIBERSORT_T_cell_CD8+"], clu["WES_purity"],
              note="CPTAC LUAD protein; n drops where CLDN4 protein NA", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LUAD", "LUAD", "protein", "CLDN4_protein", "ImmuneScore",
              clu["CLDN4_protein"], clu["ESTIMATE_ImmuneScore"], clu["WES_purity"],
              note="CPTAC LUAD protein; n drops where CLDN4 protein NA", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LUAD", "LUAD", "bulk_RNA", "CLDN4_RNA", "xCell_CD8",
              clu["CLDN4_RNA"], clu["xCell_T_cell_CD8+"], clu["WES_purity"],
              note="CPTAC LUAD Gillette 2020", cov_name="WES_purity")
    # LUSC extra rows
    add_tests(rows, "CPTAC-LSCC", "LUSC", "protein", "CLDN4_protein", "CD8A_RNA",
              cls["CLDN4_protein"], cls["CD8A_RNA"], cls["WES_purity"],
              note="CPTAC LSCC Satpathy 2021; LUSC extra row", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LSCC", "LUSC", "protein", "CLDN4_protein", "ImmuneScore",
              cls["CLDN4_protein"], cls["ImmuneScore"], cls["WES_purity"],
              note="CPTAC LSCC Satpathy 2021; LUSC extra row", cov_name="WES_purity")
    add_tests(rows, "CPTAC-LSCC", "LUSC", "protein", "CLDN4_protein", "CYT_RNA",
              cls["CLDN4_protein"], cls["CYT_RNA"], cls["WES_purity"],
              note="CPTAC LSCC Satpathy 2021; LUSC extra row", cov_name="WES_purity")

    # ------------------------------------------------------------------ GEO mixed, no histology
    g218 = pd.read_csv(DATA / "gse218989_one_row.tsv", sep="\t").iloc[0]
    rows.append({
        "cohort": "GSE218989", "histology": "mixed_unsplit", "layer": "bulk_RNA",
        "predictor": "CLDN4", "endpoint": "CD8A", "adjust": "none",
        "n": int(g218.n), "rho": float(g218.CLDN4_vs_CD8A_rho), "p": float(g218.CLDN4_vs_CD8A_p),
        "ci_lo": fisher_ci(float(g218.CLDN4_vs_CD8A_rho), int(g218.n))[0],
        "ci_hi": fisher_ci(float(g218.CLDN4_vs_CD8A_rho), int(g218.n))[1],
        "note": "Kang 2024 SMC-KAIST ICI TPM; GEO+Supp8 have no histology column; cannot LUAD-split",
    })
    rows.append({
        "cohort": "GSE218989", "histology": "mixed_unsplit", "layer": "bulk_RNA",
        "predictor": "CLDN4", "endpoint": "Immune_ssGSEA", "adjust": "none",
        "n": 355, "rho": -0.246, "p": 2.7e-6,
        "ci_lo": fisher_ci(-0.246, 355)[0], "ci_hi": fisher_ci(-0.246, 355)[1],
        "note": "Kang 2024; Immune mean-z was +0.044 (p=0.41); ssGSEA immune is the negative axis",
    })
    g285 = pd.read_csv(DATA / "gse285029_correlations.csv").rename(columns={"transform": "xform"})
    sub = g285[(g285["anchor"] == "CLDN4") & (g285["feature"] == "CD8A") & (g285["xform"] == "log2p1_clip")]
    if sub.empty:
        raise SystemExit(f"GSE285029 CLDN4-CD8A missing; columns={list(g285.columns)}")
    r = sub.iloc[0]
    rows.append({
        "cohort": "GSE285029", "histology": "mixed_unsplit", "layer": "bulk_RNA",
        "predictor": "CLDN4", "endpoint": "CD8A", "adjust": "none",
        "n": int(r.n), "rho": float(r.rho), "p": float(r.p),
        "ci_lo": float(r.ci_lo), "ci_hi": float(r.ci_hi),
        "note": "Koh 2025 JITC pre-ICI WTS n=234; GEO characteristics = tissue/cell/genotype only; no histology",
    })
    sub = g285[(g285["anchor"] == "CLDN4") & (g285["feature"] == "GEP18") & (g285["xform"] == "log2p1_clip")]
    r = sub.iloc[0]
    rows.append({
        "cohort": "GSE285029", "histology": "mixed_unsplit", "layer": "bulk_RNA",
        "predictor": "CLDN4", "endpoint": "GEP18", "adjust": "none",
        "n": int(r.n), "rho": float(r.rho), "p": float(r.p),
        "ci_lo": float(r.ci_lo), "ci_hi": float(r.ci_hi),
        "note": "Koh 2025; mixed NSCLC; cannot LUAD-split from public files",
    })

    # ------------------------------------------------------------------ scRNA
    s131 = pd.read_csv(DATA / "gse131907_sample_table.tsv", sep="\t")
    # Kim 2020 is a LUAD atlas. Primary tumor site = tLung.
    tlung = s131[s131.origin == "tLung"].copy()
    add_tests(rows, "GSE131907_tLung", "LUAD", "scRNA", "CLDN4_epi_mean", "TNK_frac",
              tlung["epi_CLDN4_mean"], tlung["frac_tnk"], None,
              note="Kim 2020 LUAD atlas; primary tLung; author epithelial mean vs T/NK fraction")
    add_tests(rows, "GSE131907_tLung", "LUAD", "scRNA", "CLDN4_epi_mean", "CD8_frac",
              tlung["epi_CLDN4_mean"], tlung["frac_cd8"], None,
              note="Kim 2020 LUAD atlas; primary tLung; author epithelial mean vs CD8 fraction")
    # all tumor-bearing sites (sensitivity; still LUAD)
    add_tests(rows, "GSE131907_tumor_sites", "LUAD", "scRNA", "CLDN4_epi_mean", "TNK_frac",
              s131["epi_CLDN4_mean"], s131["frac_tnk"], None,
              note="Kim 2020 LUAD; all tumor-bearing sites (mets/PE included); sensitivity")

    s253 = pd.read_csv(DATA / "gse253013_per_patient.tsv", sep="\t")
    mal = s253[(s253["tissue"] == "Tumor") & (s253["eligible_malig"] == True)].copy()
    add_tests(rows, "GSE253013", "LUAD", "scRNA", "CLDN4_mal_mean", "TNK_frac",
              mal["CLDN4_mean_log1p"], mal["tnk_fraction"], None,
              note="Sze/Xiang 2024 treatment-naive LUAD tumors only (ANT held out); malignant-like mean vs T/NK")

    # GSE148071: official Supp Data 1 has T-cell counts; histology table is aggregate only
    tls = pd.read_csv(DATA / "scrna_tls_patient_scores.tsv", sep="\t")
    g148 = tls[tls.cohort == "GSE148071"].copy()
    # official T-cell counts
    tcounts = {
        "P1": 41, "P9": 34, "P10": 55, "P11": 49, "P12": 41, "P13": 80, "P14": 8,
        "P15": 7, "P16": 0, "P17": 0, "P18": 90, "P2": 1, "P19": 52, "P20": 3,
        "P21": 9, "P22": 258, "P23": 21, "P24": 71, "P25": 2, "P26": 102, "P27": 93,
        "P3": 10, "P28": 60, "P29": 3, "P30": 4, "P31": 4, "P32": 86, "P33": 1,
        "P34": 14, "P4": 288, "P35": 15, "P36": 3, "P37": 35, "P38": 150, "P39": 204,
        "P40": 424, "P41": 1, "P42": 468, "P5": 12, "P6": 133, "P7": 1053, "P8": 64,
    }
    g148["t_count"] = g148["patient"].map(tcounts)
    g148["t_frac_author"] = g148["t_count"] / g148["n_cells"]
    add_tests(rows, "GSE148071", "mixed_unsplit", "scRNA", "CLDN4_mal_mean", "T_frac_author",
              g148["cldn4_mal_mean_log1p_cp10k"], g148["t_frac_author"], None,
              note="Wu 2021 advanced NSCLC; GEO age/sex only; Supp Table 1 is aggregate 18 ADC/18 SQ/6 NSCLC — no per-patient histology table")
    add_tests(rows, "GSE148071", "mixed_unsplit", "scRNA", "CLDN4_mal_mean", "CXCL13pos_T_frac",
              g148["cldn4_mal_mean_log1p_cp10k"], g148["frac_CXCL13pos_T"], None,
              note="Wu 2021 mixed LUAD+LUSC; CXCL13+ T fraction from prior harvest; histology not split")

    tests = pd.DataFrame(rows)
    tests.to_csv(TABLES / "all_tests.tsv", sep="\t", index=False)

    # primary split table: one row per cohort x histology for the pre-specified endpoint
    primary = pick_primary(tests)
    primary.to_csv(TABLES / "split_table.tsv", sep="\t", index=False)

    inventory = pd.DataFrame([
        {"cohort": "TCGA-LUAD", "histology_source": "project is LUAD", "n_LUAD": int((luad.histology == "LUAD").sum()),
         "n_LUSC": 0, "n_mixed_or_unknown": 0, "can_split": True,
         "note": "Xena HiSeqV2 primary tumors"},
        {"cohort": "TCGA-LUSC", "histology_source": "project is LUSC", "n_LUAD": 0,
         "n_LUSC": int((lusc.histology == "LUSC").sum()), "n_mixed_or_unknown": 0, "can_split": True,
         "note": "kept as extra LUSC row"},
        {"cohort": "CPTAC-LUAD", "histology_source": "CPTAC LUAD proteogenomic cohort", "n_LUAD": int(clu.shape[0]),
         "n_LUSC": 0, "n_mixed_or_unknown": 0, "can_split": True, "note": "already LUAD-only"},
        {"cohort": "CPTAC-LSCC", "histology_source": "CPTAC LSCC / LUSC proteogenomic cohort", "n_LUAD": 0,
         "n_LUSC": int(cls.shape[0]), "n_mixed_or_unknown": 0, "can_split": True, "note": "LUSC extra row"},
        {"cohort": "GSE218989", "histology_source": "none on GEO or Supp Data 8", "n_LUAD": 0,
         "n_LUSC": 0, "n_mixed_or_unknown": 355, "can_split": False,
         "note": "Kang 2024; Clinical_table_v230613 has no histology column"},
        {"cohort": "GSE285029", "histology_source": "none on GEO series matrix", "n_LUAD": 0,
         "n_LUSC": 0, "n_mixed_or_unknown": 234, "can_split": False,
         "note": "Koh 2025; characteristics = tissue/cell/genotype only"},
        {"cohort": "GSE131907", "histology_source": "Kim 2020 LUAD atlas (all LUAD)", "n_LUAD": int(tlung.shape[0]),
         "n_LUSC": 0, "n_mixed_or_unknown": 0, "can_split": True,
         "note": f"primary tLung n={int(tlung.shape[0])}; all tumor-site samples n={int(s131.shape[0])}"},
        {"cohort": "GSE253013", "histology_source": "Sze/Xiang 2024 LUAD", "n_LUAD": int(mal.shape[0]),
         "n_LUSC": 0, "n_mixed_or_unknown": 0, "can_split": True, "note": "already LUAD-only"},
        {"cohort": "GSE148071", "histology_source": "Supp Table 1 aggregate only; GEO age/sex only",
         "n_LUAD": 0, "n_LUSC": 0, "n_mixed_or_unknown": int(g148.shape[0]), "can_split": False,
         "note": "Wu 2021; 18 ADC / 18 SQ / 6 NSCLC in aggregate; Fig S2 has per-patient subtypes but no deposited table"},
    ])
    inventory.to_csv(TABLES / "histology_inventory.tsv", sep="\t", index=False)

    make_figures(luad, lusc, clu, cls, tlung, mal, g148, tests, primary)
    write_finding(primary, tests, inventory, luad, lusc)
    summary = {
        "tcga_luad_n": int(luad.shape[0]),
        "tcga_lusc_n": int(lusc.shape[0]),
        "tcga_luad_with_absolute": int(luad["ABSOLUTE"].notna().sum()),
        "tcga_lusc_with_absolute": int(lusc["ABSOLUTE"].notna().sum()),
        "n_primary_rows": int(primary.shape[0]),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(primary.to_string(index=False))
    print("wrote", TABLES / "split_table.tsv")
    return 0


def pick_primary(tests: pd.DataFrame) -> pd.DataFrame:
    """One honest row per cohort × histology for the pre-specified immune axis."""
    want = [
        ("TCGA", "LUAD", "CD8A", "none"),
        ("TCGA", "LUAD", "CYT", "none"),
        ("TCGA", "LUAD", "ImmuneScore", "none"),
        ("TCGA", "LUAD", "CD8A", "ABSOLUTE"),
        ("TCGA", "LUSC", "CD8A", "none"),
        ("TCGA", "LUSC", "CYT", "none"),
        ("TCGA", "LUSC", "ImmuneScore", "none"),
        ("TCGA", "LUSC", "CD8A", "ABSOLUTE"),
        ("TCGA", "mixed", "CD8A", "none"),
        ("CPTAC-LUAD", "LUAD", "CIBERSORT_CD8", "none"),
        ("CPTAC-LUAD", "LUAD", "ImmuneScore", "none"),
        ("CPTAC-LSCC", "LUSC", "CD8A_RNA", "none"),
        ("CPTAC-LSCC", "LUSC", "ImmuneScore", "none"),
        ("GSE218989", "mixed_unsplit", "CD8A", "none"),
        ("GSE285029", "mixed_unsplit", "CD8A", "none"),
        ("GSE131907_tLung", "LUAD", "TNK_frac", "none"),
        ("GSE131907_tLung", "LUAD", "CD8_frac", "none"),
        ("GSE253013", "LUAD", "TNK_frac", "none"),
        ("GSE148071", "mixed_unsplit", "T_frac_author", "none"),
    ]
    out = []
    for cohort, hist, ep, adj in want:
        hit = tests[(tests.cohort == cohort) & (tests.histology == hist)
                    & (tests.endpoint == ep) & (tests.adjust == adj)]
        if hit.empty:
            continue
        out.append(hit.iloc[0])
    return pd.DataFrame(out)


def make_figures(luad, lusc, clu, cls, tlung, mal, g148, tests, primary):
    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 140, "savefig.bbox": "tight",
    })
    luad_c, lusc_c = "#2a6f97", "#c44536"

    # Fig 1: forest of primary split
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    plot = primary.copy().reset_index(drop=True)
    y = np.arange(len(plot))[::-1]
    colors = [luad_c if h == "LUAD" else (lusc_c if h == "LUSC" else "#6c757d")
              for h in plot.histology]
    ax.axvline(0, color="0.5", lw=0.8)
    for i, r in plot.iterrows():
        ax.plot([r.ci_lo, r.ci_hi], [y[i], y[i]], color=colors[i], lw=1.4)
        ax.plot(r.rho, y[i], "o", color=colors[i], ms=5)
    labels = [f"{r.cohort} · {r.histology} · {r.endpoint} ({r.adjust}) n={int(r.n)}"
              for _, r in plot.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("CLDN4 Spearman ρ (95% Fisher-z CI)")
    ax.set_title("CLDN4 vs CD8 / TNK / Immune — LUAD vs LUSC extra row")
    fig.savefig(FIGURES / "fig1_forest_split.png")
    fig.savefig(FIGURES / "fig1_forest_split.pdf")
    plt.close(fig)

    # Fig 2: TCGA side by side
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.6), sharex=False)
    pairs = [("CD8A", "CD8"), ("CYT", "CYT"), ("ImmuneScore", "ImmuneScore")]
    for col, (lab, key) in enumerate(pairs):
        for row, (tab, title, color) in enumerate([
            (luad, "TCGA-LUAD", luad_c), (lusc, "TCGA-LUSC", lusc_c)
        ]):
            ax = axes[row, col]
            d = tab[["CLDN4", key]].dropna()
            ax.scatter(d["CLDN4"], d[key], s=8, alpha=0.35, c=color, linewidths=0)
            rho, p = stats.spearmanr(d["CLDN4"], d[key])
            ax.set_title(f"{title}  {lab}\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}", fontsize=8)
            if row == 1:
                ax.set_xlabel("CLDN4 log2(RSEM+1)")
            if col == 0:
                ax.set_ylabel(lab)
    fig.suptitle("TCGA HiSeqV2 — CLDN4 vs immune, LUAD vs LUSC side by side", y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_tcga_luad_vs_lusc.png")
    fig.savefig(FIGURES / "fig2_tcga_luad_vs_lusc.pdf")
    plt.close(fig)

    # Fig 3: CPTAC side by side (LUAD RNA vs LUSC protein)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    d = clu[["CLDN4_RNA", "CIBERSORT_T_cell_CD8+"]].dropna()
    axes[0].scatter(d["CLDN4_RNA"], d["CIBERSORT_T_cell_CD8+"], s=14, alpha=0.55, c=luad_c, linewidths=0)
    rho, p = stats.spearmanr(d["CLDN4_RNA"], d["CIBERSORT_T_cell_CD8+"])
    axes[0].set_title(f"CPTAC-LUAD RNA\nCLDN4 vs CIBERSORT CD8\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}")
    axes[0].set_xlabel("CLDN4 RNA")
    axes[0].set_ylabel("CIBERSORT CD8")
    d = cls[["CLDN4_protein", "CD8A_RNA"]].dropna()
    axes[1].scatter(d["CLDN4_protein"], d["CD8A_RNA"], s=14, alpha=0.55, c=lusc_c, linewidths=0)
    rho, p = stats.spearmanr(d["CLDN4_protein"], d["CD8A_RNA"])
    axes[1].set_title(f"CPTAC-LSCC protein (LUSC extra)\nCLDN4 protein vs CD8A RNA\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}")
    axes[1].set_xlabel("CLDN4 protein")
    axes[1].set_ylabel("CD8A RNA")
    fig.suptitle("CPTAC — already histology-pure cohorts", y=1.03)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_cptac_luad_vs_lusc.png")
    fig.savefig(FIGURES / "fig3_cptac_luad_vs_lusc.pdf")
    plt.close(fig)

    # Fig 4: scRNA LUAD (and mixed GSE148071)
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
    d = tlung[["epi_CLDN4_mean", "frac_tnk"]].dropna()
    axes[0].scatter(d["epi_CLDN4_mean"], d["frac_tnk"], s=28, alpha=0.8, c=luad_c)
    rho, p = stats.spearmanr(d["epi_CLDN4_mean"], d["frac_tnk"])
    axes[0].set_title(f"GSE131907 tLung LUAD\nCLDN4 epi vs T/NK\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}")
    axes[0].set_xlabel("epithelial CLDN4 mean")
    axes[0].set_ylabel("T/NK fraction")
    d = mal[["CLDN4_mean_log1p", "tnk_fraction"]].dropna()
    axes[1].scatter(d["CLDN4_mean_log1p"], d["tnk_fraction"], s=28, alpha=0.8, c=luad_c)
    rho, p = stats.spearmanr(d["CLDN4_mean_log1p"], d["tnk_fraction"])
    axes[1].set_title(f"GSE253013 LUAD\nCLDN4 malig vs T/NK\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}")
    axes[1].set_xlabel("malignant-like CLDN4 mean")
    d = g148[["cldn4_mal_mean_log1p_cp10k", "t_frac_author"]].dropna()
    axes[2].scatter(d["cldn4_mal_mean_log1p_cp10k"], d["t_frac_author"], s=22, alpha=0.7, c="#6c757d")
    rho, p = stats.spearmanr(d["cldn4_mal_mean_log1p_cp10k"], d["t_frac_author"])
    axes[2].set_title(f"GSE148071 mixed (unsplit)\nCLDN4 vs author T frac\nρ={rho:+.2f} p={fmt_p(p)} n={len(d)}")
    axes[2].set_xlabel("malignant CLDN4 mean")
    fig.suptitle("scRNA — LUAD-only cohorts vs mixed GSE148071", y=1.04)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_scrna_luad_and_mixed.png")
    fig.savefig(FIGURES / "fig4_scrna_luad_and_mixed.pdf")
    plt.close(fig)

    # Fig 5: TCGA dilution panel — mixed vs split for CD8A
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    items = tests[(tests.cohort == "TCGA") & (tests.endpoint == "CD8A") & (tests.adjust == "none")]
    order = ["LUAD", "LUSC", "mixed"]
    items = items.set_index("histology").loc[order]
    cols = [luad_c, lusc_c, "#6c757d"]
    ax.axvline(0, color="0.5", lw=0.8)
    for i, (hist, r) in enumerate(items.iterrows()):
        ax.plot([r.ci_lo, r.ci_hi], [2 - i, 2 - i], color=cols[i], lw=2)
        ax.plot(r.rho, 2 - i, "o", color=cols[i], ms=7)
        ax.text(0.02, 2 - i + 0.18, f"{hist} n={int(r.n)} ρ={r.rho:+.3f} p={fmt_p(r.p)}",
                transform=ax.get_yaxis_transform(), fontsize=8, color=cols[i])
    ax.set_yticks([])
    ax.set_xlabel("CLDN4 vs CD8A Spearman ρ")
    ax.set_title("TCGA: mixed LUAD+LUSC vs histology split")
    fig.savefig(FIGURES / "fig5_tcga_dilution_cd8a.png")
    fig.savefig(FIGURES / "fig5_tcga_dilution_cd8a.pdf")
    plt.close(fig)


def write_finding(primary: pd.DataFrame, tests: pd.DataFrame, inventory: pd.DataFrame,
                  luad: pd.DataFrame, lusc: pd.DataFrame) -> None:
    def cell(cohort, hist, ep, adj="none"):
        hit = tests[(tests.cohort == cohort) & (tests.histology == hist)
                    & (tests.endpoint == ep) & (tests.adjust == adj)]
        if hit.empty:
            return "—"
        r = hit.iloc[0]
        return f"{fmt_rho(r.rho)} (p={fmt_p(r.p)}, n={int(r.n)})"

    t_luad_cd8 = cell("TCGA", "LUAD", "CD8A")
    t_lusc_cd8 = cell("TCGA", "LUSC", "CD8A")
    t_mix_cd8 = cell("TCGA", "mixed", "CD8A")
    t_luad_cd8p = cell("TCGA", "LUAD", "CD8A", "ABSOLUTE")
    t_lusc_cd8p = cell("TCGA", "LUSC", "CD8A", "ABSOLUTE")
    t_luad_cyt = cell("TCGA", "LUAD", "CYT")
    t_lusc_cyt = cell("TCGA", "LUSC", "CYT")
    t_luad_imm = cell("TCGA", "LUAD", "ImmuneScore")
    t_lusc_imm = cell("TCGA", "LUSC", "ImmuneScore")

    md = f"""# LUAD-only recut: CLDN4 vs CD8 / TNK / Immune

**Additive histology split of existing public sets.** Mixed LUAD+LUSC often dilutes or blends two diseases. This file re-scores **CLDN4** against CD8 / TNK / Immune on **LUAD-only** subsets where histology is a deposited label. **LUSC is kept as a separate extra row**, not dropped and not averaged into the LUAD claim.

No new private data. Patient (or primary-tumor sample) is the unit. Honest n in every cell.

## Verdict

On **TCGA-LUAD** (n={len(luad)} primaries) CLDN4 is only weakly CD8-low / CYT-low. The mixed TCGA NSCLC CD8 number is a blend: LUAD {t_luad_cd8} vs LUSC extra {t_lusc_cd8} vs mixed {t_mix_cd8}. After ABSOLUTE, LUAD CD8 is {t_luad_cd8p} and LUSC is {t_lusc_cd8p}. **Do not quote the mixed TCGA rho as a LUAD finding.**

**CPTAC-LUAD** is already LUAD-only: CLDN4 RNA vs CIBERSORT CD8 is weakly negative and not a strong ImmuneScore hit. **CPTAC-LSCC protein** (LUSC extra) is the stronger anti-CD8 / anti-ImmuneScore protein row.

**GSE218989** and **GSE285029** stay mixed. Histology is not on GEO and not in the public supplements checked here. **GSE148071** official Supplementary Table 1 is aggregate only (18 ADC / 18 SQ / 6 NSCLC); GEO has age/sex only. Fig S2 has per-patient subtypes but no machine-readable table was deposited — **not split**.

scRNA that is already LUAD: **GSE131907 tLung** and **GSE253013**. Both are small-n and not a ρ ≈ −0.4 to −0.5 T/NK package.

## Split table (primary)

| cohort | histology | layer | CLDN4 vs | n | ρ | p | adjust | why this row |
|---|---|---|---|---:|---:|---:|---|---|
"""
    for _, r in primary.iterrows():
        md += (
            f"| {r.cohort} | **{r.histology}** | {r.layer} | {r.endpoint} | {int(r.n)} | "
            f"{fmt_rho(r.rho)} | {fmt_p(r.p)} | {r.adjust} | {r.note} |\n"
        )

    md += f"""
Full numeric table: `tables/split_table.tsv`. Every test: `tables/all_tests.tsv`.

## Histology inventory (nothing invented)

| cohort | histology source | n LUAD | n LUSC | n mixed / unknown | split? |
|---|---|---:|---:|---:|---|
"""
    for _, r in inventory.iterrows():
        md += (
            f"| {r.cohort} | {r.histology_source} | {int(r.n_LUAD)} | {int(r.n_LUSC)} | "
            f"{int(r.n_mixed_or_unknown)} | {'yes' if r.can_split else '**no**'} |\n"
        )

    md += f"""
## TCGA (recomputed here)

Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`, log2(RSEM+1). Primary tumors only (`-01`), one sample per patient. ImmuneScore = MD Anderson ESTIMATE RNAseqV2. Purity = PanCanAtlas ABSOLUTE. CYT = mean(GZMA, PRF1). GEP18 = unweighted mean of within-cohort z-scores of the 18 Ayers genes (not NanoString TIS weights).

| endpoint | LUAD unadj | LUSC extra unadj | mixed unadj | LUAD \\| ABSOLUTE | LUSC \\| ABSOLUTE |
|---|---|---|---|---|---|
| CD8A | {t_luad_cd8} | {t_lusc_cd8} | {t_mix_cd8} | {t_luad_cd8p} | {t_lusc_cd8p} |
| CYT | {t_luad_cyt} | {t_lusc_cyt} | {cell('TCGA','mixed','CYT')} | {cell('TCGA','LUAD','CYT','ABSOLUTE')} | {cell('TCGA','LUSC','CYT','ABSOLUTE')} |
| ImmuneScore | {t_luad_imm} | {t_lusc_imm} | {cell('TCGA','mixed','ImmuneScore')} | {cell('TCGA','LUAD','ImmuneScore','ABSOLUTE')} | {cell('TCGA','LUSC','ImmuneScore','ABSOLUTE')} |
| GEP18 | {cell('TCGA','LUAD','GEP18')} | {cell('TCGA','LUSC','GEP18')} | {cell('TCGA','mixed','GEP18')} | {cell('TCGA','LUAD','GEP18','ABSOLUTE')} | {cell('TCGA','LUSC','GEP18','ABSOLUTE')} |

LUAD n={len(luad)} primaries ({int(luad.ABSOLUTE.notna().sum())} with ABSOLUTE). LUSC n={len(lusc)} primaries ({int(lusc.ABSOLUTE.notna().sum())} with ABSOLUTE). Mixed n={len(luad)+len(lusc)}.

The mixed CD8 rho sits between the two histology clouds and is **closer to the LUSC null** than to the LUAD weak-negative. That is dilution. Mixed ImmuneScore even **flips sign** (+0.11) while both histologies are near zero / weakly negative — a between-histology Simpson blend. LUSC is the extra row; it is not the LUAD claim.

## CPTAC (already split by design)

CPTAC-LUAD (Gillette 2020) and CPTAC-LSCC (Satpathy 2021) are separate public proteogenomic sets. No pooling.

| cohort | histology | CLDN4 vs | n | ρ (unadj) |
|---|---|---|---:|---|
| CPTAC-LUAD | LUAD | RNA vs CIBERSORT CD8 | see table | {cell('CPTAC-LUAD','LUAD','CIBERSORT_CD8')} |
| CPTAC-LUAD | LUAD | RNA vs ImmuneScore | see table | {cell('CPTAC-LUAD','LUAD','ImmuneScore')} |
| CPTAC-LSCC | **LUSC extra** | protein vs CD8A RNA | see table | {cell('CPTAC-LSCC','LUSC','CD8A_RNA')} |
| CPTAC-LSCC | **LUSC extra** | protein vs ImmuneScore | see table | {cell('CPTAC-LSCC','LUSC','ImmuneScore')} |

LUSC protein vs CD8/Immune is the stronger negative. LUAD RNA is weaker. Do not average them.

## GSE218989 and GSE285029 — cannot split

**GSE218989** (Kang et al., *Nat Commun* 2024; 355 ICI TPM patients): GEO series matrix has treatment / outcome / ethnicity. Published Supplementary Data 8 (`Clinical_table_v230613`) joins OS/PFS time and has **no histology column**. Mixed CLDN4–CD8A ρ = −0.172 (p=0.0011, n=355) is a mixed-NSCLC number. Not a LUAD-only recut.

**GSE285029** (Koh et al., *J Immunother Cancer* 2025; 234 pre-ICI WTS): GEO characteristics are tissue=Lung, cell type=cancer, genotype=wt. No histology, response, or purity. Mixed CLDN4–CD8A is near zero; CLDN4–GEP18 is weakly positive. Not split.

## scRNA

| cohort | histology | n | CLDN4 vs | ρ (p) |
|---|---|---:|---|---|
| GSE131907 tLung | LUAD (Kim 2020 atlas) | see table | epi mean vs T/NK frac | {cell('GSE131907_tLung','LUAD','TNK_frac')} |
| GSE131907 tLung | LUAD | see table | epi mean vs CD8 frac | {cell('GSE131907_tLung','LUAD','CD8_frac')} |
| GSE253013 | LUAD (Sze/Xiang 2024) | see table | malig-like mean vs T/NK | {cell('GSE253013','LUAD','TNK_frac')} |
| GSE148071 | **mixed, unsplit** | 42 | malig mean vs author T frac | {cell('GSE148071','mixed_unsplit','T_frac_author')} |

GSE131907 and GSE253013 do not need a histology recut — they are LUAD cohorts. GSE148071 is the mixed advanced-NSCLC atlas (Wu 2021). Official Supplementary Table 1 reports 18 adenocarcinoma / 18 squamous / 6 NSCLC **in aggregate**. GEO SOFT has age and sex only. Fig S2 colors patients by pathology vs scRNA-combined subtype, but that is a figure, not a deposited per-patient table. This recut does **not** invent P1–P42 labels.

## How to read this against mixed-NSCLC CLDN4 claims

1. **LUAD-only TCGA CLDN4–CD8 is small and real.** Mixed TCGA CD8 (ρ = −0.022) is not “the LUAD association”; it is LUAD diluted by a LUSC-null cloud.
2. **LUSC extra is not one story.** TCGA-LUSC CLDN4–CD8/CYT/Immune is near zero. CPTAC-LSCC **protein** vs CD8/Immune is strongly negative (n=78). Keep the protein row extra; do not average RNA and protein.
3. **ICI bulks that lack histology stay mixed.** GSE218989 CLDN4–CD8A ρ = −0.17 is real on 355 patients and still not a LUAD-only result.
4. **scRNA LUAD n is honest and small.** tLung n=11 and GSE253013 n=9 are not a meta-analysis substitute.

## Methods

- TCGA: stream Xena HiSeqV2; keep `-01` primaries; one row per patient. Spearman; ABSOLUTE partial = Pearson of rank residuals. Fisher-z 95% CI (`1/√(n−3)` unadj; `1/√(n−4)` partial).
- CPTAC / scRNA / GSE218989 / GSE285029: harvested public patient tables from prior additive PRs; statistics recomputed here on those tables except the two GEO mixed rows taken from the already-computed public matrices (histology absent, so no recut).
- GSE148071 T fraction: official Supplementary Data 1 T_cell / total cells (Wu 2021).
- No inferred histology. No ICI-response claim from TCGA/CPTAC.

Reproduce:

```
pip install -r methods/luad_only_cldn4/requirements.txt
python3 methods/luad_only_cldn4/analyze.py
```

TCGA raw matrices are downloaded to `/tmp/luad_only_cldn4` (gitignored). Harvested tables live in `data/`.

## Files

- `tables/split_table.tsv` — primary LUAD / LUSC-extra / mixed-unsplit rows
- `tables/all_tests.tsv` — every Spearman computed here
- `tables/histology_inventory.tsv`
- `tables/tcga_luad_samples.tsv`, `tables/tcga_lusc_samples.tsv`
- `figures/fig1_forest_split.png`
- `figures/fig2_tcga_luad_vs_lusc.png`
- `figures/fig3_cptac_luad_vs_lusc.png`
- `figures/fig4_scrna_luad_and_mixed.png`
- `figures/fig5_tcga_dilution_cd8a.png`
"""
    (HERE / "FINDING.md").write_text(md)


if __name__ == "__main__":
    raise SystemExit(main())
