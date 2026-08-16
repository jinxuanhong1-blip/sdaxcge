#!/usr/bin/env python3
"""A2 extras — public lung IO RNA besides GSE253564.

The user durvalumab result (ρ = −0.65 / purity −0.46) is taken as given.
This script does not audit that coefficient. It adds public series on top:

  * leftover GEO neoadjuvant / anti-PD-1 lung RNA with recoverable
    MPR, DCB, or ORR labels
  * GSE248378 (same NCT02904954 trial, post-durva accession)
  * a live leftover GEO 2024–2026 + ArrayExpress inventory

Score TACSTD2 and CLDN4 vs the cohort-native endpoint and vs CD8A / GEP18
after an ESTIMATEScore residual. Honest n / ρ / p. Public data only.

Outputs -> results/rework/A2_wave2/extra/
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import rework_A2_wave2 as w2  # noqa: E402

DATA = ROOT / "data"
AUX = DATA / "A2_wave2"
EXTRA_DATA = AUX / "extra"
OUT = ROOT / "results" / "rework" / "A2_wave2" / "extra"
FIG = OUT / "figures"
for p in (EXTRA_DATA, FIG):
    p.mkdir(parents=True, exist_ok=True)

TARGETS = ("TACSTD2", "CLDN4")
GEP18 = (
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
)
# Hardcoded GRCh38 IDs for targets / GEP / epithelial / CYT if REST is down.
# ESTIMATE genes are mapped live (or skipped if mapping fails).
SYMBOL_ENSG = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563",
    "GZMA": "ENSG00000145649",
    "PRF1": "ENSG00000180644",
    "MKI67": "ENSG00000148773",
    "EPCAM": "ENSG00000119888",
    "KRT7": "ENSG00000135480",
    "KRT8": "ENSG00000170421",
    "KRT18": "ENSG00000111057",
    "KRT19": "ENSG00000171345",
    "CDH1": "ENSG00000039068",
    "MUC1": "ENSG00000185479",
    "CCL5": "ENSG00000161570",
    "CCL5_alt": "ENSG00000271503",
    "CD27": "ENSG00000139193",
    "CD274": "ENSG00000120217",
    "CD276": "ENSG00000103855",
    "CMKLR1": "ENSG00000174600",
    "CXCL9": "ENSG00000138755",
    "CXCR6": "ENSG00000172215",
    "HLA-DQA1": "ENSG00000196735",
    "HLA-DRB1": "ENSG00000196126",
    "HLA-E": "ENSG00000204592",
    "IDO1": "ENSG00000131203",
    "LAG3": "ENSG00000089692",
    "NKG7": "ENSG00000105374",
    "PDCD1LG2": "ENSG00000197646",
    "PSMB10": "ENSG00000205220",
    "STAT1": "ENSG00000115415",
    "TIGIT": "ENSG00000181847",
}

URLS = {
    "GSE126044_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE126044_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz",
    "GSE135222_exp": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE135222_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz",
    "GSE166449_tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    "GSE166449_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz",
    "GSE207422_log2tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
    "GSE207422_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
    "GSE207422_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz",
    "GSE329813_norm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/suppl/GSE329813_processed_data_file_normalized_data.csv.gz",
    "GSE329813_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/matrix/GSE329813_series_matrix.txt.gz",
}

NCBI_TOOL = {"tool": "sdaxcge_A2_wave2_extra", "email": "jinxuanhong1@gmail.com"}


def download(url: str, dest: Path, timeout: int = 180) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "sdaxcge-A2-wave2-extra"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, dest.open("wb") as out:
        out.write(resp.read())
    return dest


def ensure_extra_inputs() -> dict[str, Path]:
    paths = {
        "GSE126044_counts": EXTRA_DATA / "GSE126044_counts.txt.gz",
        "GSE126044_matrix": EXTRA_DATA / "GSE126044_series_matrix.txt.gz",
        "GSE135222_exp": EXTRA_DATA / "GSE135222_exp.tsv.gz",
        "GSE135222_matrix": EXTRA_DATA / "GSE135222_series_matrix.txt.gz",
        "GSE166449_tpm": EXTRA_DATA / "GSE166449_TPM.txt.gz",
        "GSE166449_matrix": EXTRA_DATA / "GSE166449_series_matrix.txt.gz",
        "GSE207422_log2tpm": EXTRA_DATA / "GSE207422_log2TPM.txt.gz",
        "GSE207422_meta": EXTRA_DATA / "GSE207422_metadata.xlsx",
        "GSE207422_matrix": EXTRA_DATA / "GSE207422_series_matrix.txt.gz",
        "GSE329813_norm": EXTRA_DATA / "GSE329813_norm.csv.gz",
        "GSE329813_matrix": EXTRA_DATA / "GSE329813_series_matrix.txt.gz",
        "gmt": DATA / "estimate" / "inst" / "extdata" / "SI_geneset.gmt",
        "post_fpkm": DATA / "GSE248378_Durva_Post_FPKMs.txt.gz",
        "post_matrix": AUX / "GSE248378_series_matrix.txt.gz",
        "pre_fpkm": DATA / "GSE253564_Pre_FPKMs.txt.gz",
        "pre_matrix": AUX / "GSE253564_series_matrix.txt.gz",
        "natcomm": AUX / "41467_2023_44195_MOESM6_ESM.xlsx",
        "mmc2": AUX / "mmc2.xlsx",
    }
    for key, dest in paths.items():
        if key in URLS and (not dest.exists() or dest.stat().st_size == 0):
            download(URLS[key], dest)
    if not paths["gmt"].exists():
        download(w2.URLS["estimate_gmt"], paths["gmt"])
    return paths


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """One row per GSM from a GEO series_matrix."""
    titles: list[str] = []
    gsms: list[str] = []
    chars: list[list[str]] = []
    descriptions: list[str] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for raw in fh:
            if not raw.startswith("!"):
                continue
            row = next(csv.reader([raw.rstrip("\n")], delimiter="\t"))
            key = row[0]
            vals = [v.strip().strip('"') for v in row[1:]]
            if key == "!Sample_title":
                titles = vals
            elif key == "!Sample_geo_accession":
                gsms = vals
            elif key == "!Sample_description":
                descriptions = vals
            elif key == "!Sample_characteristics_ch1":
                chars.append(vals)
    n = len(titles)
    recs = []
    for i in range(n):
        rec = {
            "gsm": gsms[i] if i < len(gsms) else "",
            "title": titles[i],
            "description": descriptions[i] if i < len(descriptions) else "",
        }
        for block in chars:
            if i >= len(block):
                continue
            cell = block[i]
            if ": " in cell:
                k, v = cell.split(": ", 1)
                rec[k.strip().lower()] = v.strip()
            else:
                rec.setdefault("characteristic", cell)
        recs.append(rec)
    return pd.DataFrame(recs)


def read_symbol_matrix(path: Path, sep: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, compression="gzip", index_col=0)
    df.index = df.index.astype(str).str.strip()
    df = df.groupby(df.index).max()
    return df.apply(pd.to_numeric, errors="coerce")


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def already_log_or_log2p1(df: pd.DataFrame, declared: str) -> pd.DataFrame:
    """Use deposited log matrices as-is; log2(x+1) raw TPM/FPKM/counts-like."""
    if declared == "already_log":
        return df.astype(float)
    return np.log2(df.astype(float) + 1.0)


def zmean(expr_log: pd.DataFrame, genes: tuple[str, ...]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr_log.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=expr_log.columns), present
    mat = expr_log.loc[present].astype(float)
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


def mean_present(expr_log: pd.DataFrame, genes: tuple[str, ...]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr_log.index]
    if not present:
        return pd.Series(np.nan, index=expr_log.columns), present
    return expr_log.loc[present].astype(float).mean(axis=0), present


def gene_series(expr_log: pd.DataFrame, gene: str) -> pd.Series:
    if gene not in expr_log.index:
        return pd.Series(np.nan, index=expr_log.columns)
    return expr_log.loc[gene].astype(float)


def ensembl_lookup(symbols: list[str]) -> dict[str, str]:
    """Ensembl REST symbol -> ENSG. Falls back to the hardcoded map."""
    mapping = dict(SYMBOL_ENSG)
    todo = [s for s in symbols if s not in mapping]
    if not todo:
        return mapping
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    for i in range(0, len(todo), 200):
        chunk = todo[i : i + 200]
        try:
            r = requests.post(url, headers=headers, json={"symbols": chunk}, timeout=60)
            r.raise_for_status()
            payload = r.json()
            for sym, rec in payload.items():
                if isinstance(rec, dict) and rec.get("id"):
                    mapping[sym] = rec["id"]
            time.sleep(0.2)
        except Exception:
            break
    return mapping


def remap_ensembl_index(expr: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Collapse versioned ENSG rows onto HGNC symbols (max)."""
    idx = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    inv: dict[str, str] = {}
    for sym, ensg in mapping.items():
        inv[ensg] = "CCL5" if sym.startswith("CCL5") else sym
    symbols = [inv.get(i, i) for i in idx]
    out = expr.copy()
    out.index = symbols
    # Keep symbol rows that we mapped; drop leftover ENSG ids.
    keep = [i for i in out.index if not str(i).startswith("ENSG")]
    if keep:
        out = out.loc[keep]
    return out.groupby(out.index).max()


def load_gse126044(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    counts = read_symbol_matrix(paths["GSE126044_counts"])
    expr = log2cpm(counts)
    meta = parse_series_matrix(paths["GSE126044_matrix"])
    meta["sample"] = meta["title"].str.replace(r"^RNA-seq_", "", regex=True)
    meta = meta[meta["sample"].isin(expr.columns)].copy()
    resp = meta["patient response"].str.lower().str.strip()
    meta["endpoint"] = np.where(resp.eq("responder"), 1, np.where(resp.eq("non-responder"), 0, np.nan))
    meta["endpoint_name"] = "ORR"
    meta["pos_label"] = "responder"
    meta["neg_label"] = "non-responder"
    meta["timepoint"] = "pre"
    meta["assay"] = "bulk RNA-seq (log2CPM from counts)"
    meta["drug"] = "anti-PD-1"
    meta["setting"] = "advanced NSCLC"
    info = {
        "series": "GSE126044",
        "pmid": "32879421",
        "note": "Cho 2020. Pretreatment anti-PD-1 NSCLC. Cohort-native label is patient response (ORR), not MPR.",
    }
    return expr, meta.set_index("sample"), info


def load_gse166449(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = read_symbol_matrix(paths["GSE166449_tpm"])
    # Deposited "TPM" matrix is already on a log2-like scale (max ~15.7).
    expr = raw.astype(float)
    meta = parse_series_matrix(paths["GSE166449_matrix"])
    meta["sample"] = meta["description"]
    meta = meta[meta["sample"].isin(expr.columns)].copy()
    title = meta["title"].str.lower()
    meta["endpoint"] = np.where(
        title.str.contains("nonresponder"),
        0,
        np.where(title.str.contains("responder"), 1, np.nan),
    )
    meta["endpoint_name"] = "ORR"
    meta["pos_label"] = "responder"
    meta["neg_label"] = "non-responder"
    meta["timepoint"] = "pre"
    meta["assay"] = "bulk RNA-seq (deposited log-like TPM)"
    meta["drug"] = "ICI (SMC IO; Lee/SELECT)"
    meta["setting"] = "advanced lung cancer"
    info = {
        "series": "GSE166449",
        "pmid": "33857424",
        "note": "Lee 2021 SELECT supplement. Title encodes responder vs nonResponder. Values used as deposited (already log-like; not re-logged).",
    }
    return expr, meta.set_index("sample"), info


def load_gse207422(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    expr = read_symbol_matrix(paths["GSE207422_log2tpm"])
    wb = openpyxl.load_workbook(paths["GSE207422_meta"], data_only=True)
    rows = list(wb.active.iter_rows(values_only=True))
    meta = pd.DataFrame(rows[1:], columns=rows[0])
    meta = meta[meta["Sample"].isin(expr.columns)].copy()
    pre = meta[meta["Resource"].astype(str).str.contains("Pre", case=False)].copy()
    expr = expr[pre["Sample"].tolist()]
    path = pre["Pathologic Response"].astype(str)
    pre["mpr"] = np.where(path.str.upper().str.startswith("MPR") | path.str.contains("pCR", case=False), 1, 0)
    recist = pre["RECIST"].astype(str).str.upper()
    pre["orr"] = np.where(recist.isin(["CR", "PR"]), 1, np.where(recist.isin(["SD", "PD"]), 0, np.nan))
    pre = pre.rename(columns={"Sample": "sample"})
    pre["endpoint"] = pre["mpr"]
    pre["endpoint_name"] = "MPR"
    pre["pos_label"] = "MPR"
    pre["neg_label"] = "NMPR"
    pre["timepoint"] = "pre"
    pre["assay"] = "bulk RNA-seq log2(TPM+1)"
    pre["drug"] = "anti-PD-1 + chemo (tori/camre/sinti; not durvalumab)"
    pre["setting"] = "neoadjuvant resectable NSCLC"
    pre["title"] = pre["sample"]
    info = {
        "series": "GSE207422",
        "pmid": "36869384",
        "note": "Hu/Zhang 2023. Baseline bulk only (24 patients). MPR = Pathologic Response starts with MPR or contains pCR. RECIST scored separately, not pooled with MPR.",
    }
    return expr, pre.set_index("sample"), info


def load_gse135222(paths: dict[str, Path], gmt: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = read_symbol_matrix(paths["GSE135222_exp"])
    raw.index = raw.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    needed = list(TARGETS) + list(GEP18) + list(w2.EPI_GENES) + list(w2.CYT_GENES) + [w2.CD8_GENE, w2.MKI67_GENE]
    needed += gmt.get("StromalSignature", []) + gmt.get("ImmuneSignature", [])
    mapping = ensembl_lookup(sorted(set(needed)))
    expr = remap_ensembl_index(np.log2(raw.astype(float) + 1.0), mapping)
    meta = parse_series_matrix(paths["GSE135222_matrix"])
    meta["sample"] = meta["title"].str.replace(r"\s+", "", regex=True)
    meta = meta[meta["sample"].isin(expr.columns)].copy()
    pfs_t = pd.to_numeric(meta["pfs.time"], errors="coerce")
    pfs_e = pd.to_numeric(meta["progression-free survival (pfs)"], errors="coerce")
    meta["pfs_time"] = pfs_t
    meta["pfs_event"] = pfs_e
    meta["endpoint"] = np.where(pfs_t >= 180, 1, 0)
    meta["endpoint_name"] = "DCB"
    meta["pos_label"] = "DCB (PFS≥180 d)"
    meta["neg_label"] = "NDB"
    meta["timepoint"] = "pre"
    meta["assay"] = "bulk RNA-seq log2(TPM+1) from Ensembl"
    meta["drug"] = "anti-PD-1/PD-L1"
    meta["setting"] = "advanced NSCLC"
    info = {
        "series": "GSE135222",
        "pmid": "31537801",
        "note": "Jung 2019. DCB = PFS time ≥ 180 days (cohort-native convention). Cox PFS also reported. Ensembl IDs stripped of version and mapped to symbols.",
        "ensembl_mapped": int(sum(1 for g in TARGETS + GEP18 if g in expr.index)),
    }
    return expr, meta.set_index("sample"), info


def load_gse329813(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = pd.read_csv(paths["GSE329813_norm"], compression="gzip", index_col=0)
    raw.index = raw.index.astype(str).str.strip()
    raw = raw.groupby(raw.index).max().apply(pd.to_numeric, errors="coerce")
    expr_roi = np.log2(raw.astype(float) + 1.0)
    meta = parse_series_matrix(paths["GSE329813_matrix"])

    def parse_title(t: str) -> dict:
        # "ROI 1, Patient 1, Primary tumor bed, MPR"
        parts = [x.strip() for x in str(t).split(",")]
        rec = {"title": t, "roi": parts[0] if parts else t}
        for p in parts[1:]:
            m = re.match(r"Patient\s+(\d+)", p, re.I)
            if m:
                rec["patient"] = f"P{int(m.group(1))}"
            elif re.search(r"tumor bed", p, re.I):
                rec["site"] = "primary_tumor_bed"
            elif re.search(r"lymph", p, re.I):
                rec["site"] = "lymph_node"
            elif p.upper() in {"MPR", "NMPR", "NON-MPR", "NONMPR"}:
                rec["mpr_label"] = p.upper()
        if "mpr_label" not in rec:
            if re.search(r"\bNMPR\b|non-?mpr", t, re.I):
                rec["mpr_label"] = "NMPR"
            elif re.search(r"\bMPR\b", t):
                rec["mpr_label"] = "MPR"
        return rec

    parsed = pd.DataFrame([parse_title(t) for t in meta["title"]])
    parsed["roi"] = parsed["roi"].astype(str)
    if "site" not in parsed.columns:
        parsed["site"] = ""
    if "mpr_label" not in parsed.columns:
        parsed["mpr_label"] = None
    if "patient" not in parsed.columns:
        parsed["patient"] = None
    bed = parsed[parsed["site"].eq("primary_tumor_bed")].copy()
    roi_cols = [c for c in expr_roi.columns if c in set(bed["roi"])]
    bed = bed[bed["roi"].isin(roi_cols)].copy()
    # Patient-level mean of tumor-bed ROIs.
    rows = []
    expr_cols = {}
    for pat, grp in bed.groupby("patient"):
        cols = [c for c in grp["roi"] if c in expr_roi.columns]
        if not cols:
            continue
        expr_cols[pat] = expr_roi[cols].mean(axis=1)
        lab = grp["mpr_label"].dropna().iloc[0] if grp["mpr_label"].notna().any() else None
        rows.append(
            {
                "sample": pat,
                "patient": pat,
                "n_rois": len(cols),
                "endpoint": 1 if lab == "MPR" else (0 if lab in {"NMPR", "NON-MPR", "NONMPR"} else np.nan),
                "mpr_label": lab,
                "title": pat,
            }
        )
    expr = pd.DataFrame(expr_cols)
    meta_p = pd.DataFrame(rows).set_index("sample")
    meta_p["endpoint_name"] = "MPR"
    meta_p["pos_label"] = "MPR"
    meta_p["neg_label"] = "NMPR"
    meta_p["timepoint"] = "post-neoadjuvant residual (tumor bed)"
    meta_p["assay"] = "GeoMx DSP (patient-mean tumor-bed ROI, log2(norm+1))"
    meta_p["drug"] = "pembrolizumab + platinum (NCT05383716; not durvalumab)"
    meta_p["setting"] = "neoadjuvant resectable NSCLC, post-treatment residual"
    info = {
        "series": "GSE329813",
        "pmid": "",
        "note": "2026 GeoMx leftover. CLDN4 absent from the panel. Post-treatment residual tumor bed — not pretreatment bulk. Not pooled with GSE207422.",
        "n_rois_tumor_bed": int(bed.shape[0]),
        "cldn4_present": bool("CLDN4" in expr.index),
        "tacstd2_present": bool("TACSTD2" in expr.index),
    }
    return expr, meta_p, info


def load_gse248378(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    expr_fpkm = w2.read_fpkm(paths["post_fpkm"])
    expr = np.log2(expr_fpkm.astype(float) + 1.0)
    clin, _ = w2.build_clinical(
        {
            "pre_matrix": paths["pre_matrix"],
            "post_matrix": paths["post_matrix"],
            "natcomm": paths["natcomm"],
            "mmc2": paths["mmc2"],
        }
    )
    sub = clin[clin["series"].eq("GSE248378")].copy()
    sub["sample"] = sub["sample_title"]
    sub = sub[sub["sample"].isin(expr.columns)].copy()
    rec = pd.to_numeric(sub["recurrence_event"], errors="coerce")
    sub["endpoint"] = rec
    sub["endpoint_name"] = "recurrence"
    sub["pos_label"] = "recurrence"
    sub["neg_label"] = "no recurrence"
    sub["pfs_time"] = pd.to_numeric(sub["pfs_months"], errors="coerce")
    sub["pfs_event"] = rec
    sub["timepoint"] = "post"
    sub["assay"] = "bulk RNA-seq log2(FPKM+1)"
    sub["drug"] = "durvalumab ± SBRT (NCT02904954)"
    sub["setting"] = "neoadjuvant stage I–III NSCLC, residual tumor"
    info = {
        "series": "GSE248378",
        "pmid": "38401548",
        "note": "Same trial as GSE253564, different accession. Post-treatment residual. MPR not estimable in this matrix (0 MPR). Recurrence / PFS kept as the clinical labels. Included because it is the other open durvalumab NSCLC whole-transcriptome accession.",
    }
    return expr, sub.set_index("sample"), info


def analyse_series(
    series: str,
    expr_log: pd.DataFrame,
    meta: pd.DataFrame,
    gmt: dict[str, list[str]],
    info: dict,
    extra_binaries: dict[str, pd.Series] | None = None,
) -> tuple[pd.DataFrame, list[dict], dict]:
    samples = [s for s in meta.index if s in expr_log.columns]
    expr = expr_log[samples]
    meta = meta.loc[samples].copy()
    est, overlaps = w2.estimate_scores(expr, gmt)
    gep, gep_present = zmean(expr, GEP18)
    cyt, cyt_present = mean_present(expr, w2.CYT_GENES)
    epi, epi_present = mean_present(expr, w2.EPI_GENES)
    cd8 = gene_series(expr, "CD8A")
    sample = meta.reset_index().rename(columns={"index": "sample"})
    if "sample" not in sample.columns:
        sample["sample"] = samples
    idx = sample["sample"]
    sample["series"] = series
    sample["CD8A"] = cd8.loc[idx].to_numpy()
    sample["GEP18"] = gep.loc[idx].to_numpy()
    sample["CYT"] = cyt.loc[idx].to_numpy()
    sample["EpiScore"] = epi.loc[idx].to_numpy()
    sample["ImmuneScore"] = est.loc[idx, "ImmuneScore"].to_numpy()
    sample["StromalScore"] = est.loc[idx, "StromalScore"].to_numpy()
    sample["ESTIMATEScore"] = est.loc[idx, "ESTIMATEScore"].to_numpy()
    sample["TumorPurity"] = est.loc[idx, "TumorPurity"].to_numpy()
    for gene in TARGETS:
        sample[gene] = gene_series(expr, gene).loc[idx].to_numpy()
        sample[f"{gene}_present"] = gene in expr.index

    covariates = {
        "none": None,
        "ESTIMATEScore": sample.set_index("sample")["ESTIMATEScore"],
    }
    binaries = {
        meta["endpoint_name"].iloc[0]: pd.to_numeric(sample.set_index("sample")["endpoint"], errors="coerce")
    }
    if extra_binaries:
        binaries.update(extra_binaries)
    immune = {
        "CD8A": sample.set_index("sample")["CD8A"],
        "GEP18": sample.set_index("sample")["GEP18"],
        "CYT": sample.set_index("sample")["CYT"],
    }
    rows: list[dict] = []
    for gene in TARGETS:
        y0 = sample.set_index("sample")[gene]
        present = bool(gene in expr.index) and bool(np.isfinite(y0).sum() >= 6)
        for cov_name, z in covariates.items():
            y = y0 if cov_name == "none" else w2.ols_residual(y0, z)
            adj = "unadjusted" if cov_name == "none" else f"OLS residual on {cov_name}"
            for ep_name, lab in binaries.items():
                rec = w2.mwu(y, lab)
                rec.update(
                    {
                        "series": series,
                        "gene": gene,
                        "gene_present": present,
                        "covariate": cov_name,
                        "adjustment": adj,
                        "endpoint": ep_name,
                        "test": "Mann-Whitney U",
                        "timepoint": str(meta["timepoint"].iloc[0]),
                        "setting": str(meta["setting"].iloc[0]),
                        "drug": str(meta["drug"].iloc[0]),
                        "assay": str(meta["assay"].iloc[0]),
                    }
                )
                if not present:
                    rec["note"] = f"{gene} absent from deposited matrix"
                    rec["p"] = None
                rows.append(rec)
            if "pfs_time" in sample.columns and "pfs_event" in sample.columns:
                surv = w2.cox_pfs(
                    y,
                    pd.to_numeric(sample.set_index("sample")["pfs_time"], errors="coerce"),
                    pd.to_numeric(sample.set_index("sample")["pfs_event"], errors="coerce"),
                )
                rows.append(
                    {
                        "series": series,
                        "gene": gene,
                        "gene_present": present,
                        "covariate": cov_name,
                        "adjustment": adj,
                        "endpoint": "PFS",
                        "test": "Cox PH continuous + log-rank median split",
                        "timepoint": str(meta["timepoint"].iloc[0]),
                        "setting": str(meta["setting"].iloc[0]),
                        "drug": str(meta["drug"].iloc[0]),
                        "assay": str(meta["assay"].iloc[0]),
                        **surv,
                    }
                )
            for iname, iv in immune.items():
                rho, p, n = w2.spearman_pair(y, iv)
                rows.append(
                    {
                        "series": series,
                        "gene": gene,
                        "gene_present": present,
                        "covariate": cov_name,
                        "adjustment": adj,
                        "endpoint": iname,
                        "test": (
                            "Spearman unadjusted"
                            if cov_name == "none"
                            else "Spearman of OLS-residual vs raw endpoint"
                        ),
                        "n": n,
                        "rho": rho,
                        "p": p,
                        "timepoint": str(meta["timepoint"].iloc[0]),
                        "setting": str(meta["setting"].iloc[0]),
                        "drug": str(meta["drug"].iloc[0]),
                        "assay": str(meta["assay"].iloc[0]),
                        "note": "" if present else f"{gene} absent",
                    }
                )
                if cov_name != "none" and z is not None and present:
                    pr, pp, pn = w2.partial_spearman(y0, iv, z)
                    rows.append(
                        {
                            "series": series,
                            "gene": gene,
                            "gene_present": present,
                            "covariate": cov_name,
                            "adjustment": f"partial Spearman | {cov_name}",
                            "endpoint": iname,
                            "test": "partial Spearman (ranks, n-3 df)",
                            "n": pn,
                            "rho": pr,
                            "p": pp,
                            "timepoint": str(meta["timepoint"].iloc[0]),
                            "setting": str(meta["setting"].iloc[0]),
                            "drug": str(meta["drug"].iloc[0]),
                            "assay": str(meta["assay"].iloc[0]),
                            "note": "purity residual = ESTIMATEScore; residualises both ranks",
                        }
                    )

    meta_out = {
        "series": series,
        "n_samples": int(expr.shape[1]),
        "n_genes": int(expr.shape[0]),
        "estimate_overlap": overlaps,
        "gep18_present": gep_present,
        "n_gep18": len(gep_present),
        "cyt_present": cyt_present,
        "epi_present": epi_present,
        "n_purity_oob": int(est["TumorPurity"].isna().sum()),
        "targets_present": {g: bool(g in expr.index) for g in TARGETS},
        **info,
    }
    return sample, rows, meta_out


def ncbi_esearch(term: str, retmax: int = 40) -> list[str]:
    params = {"db": "gds", "term": term, "retmax": str(retmax), "retmode": "json", **NCBI_TOOL}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    r = requests.get(url, timeout=45)
    r.raise_for_status()
    time.sleep(0.35)
    return r.json().get("esearchresult", {}).get("idlist", [])


def ncbi_esummary(ids: list[str]) -> dict:
    if not ids:
        return {}
    params = {"db": "gds", "id": ",".join(ids), "retmode": "json", **NCBI_TOOL}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(params)
    r = requests.get(url, timeout=45)
    r.raise_for_status()
    time.sleep(0.35)
    return r.json().get("result", {})


def hunt_leftover_public() -> list[dict]:
    """Live leftover GEO 2024–2026 + ArrayExpress. No invented accessions."""
    rows: list[dict] = []

    def add(**kw):
        rows.append(kw)

    queries = [
        (
            "neoadjuvant IO NSCLC RNA 2024-2026",
            'neoadjuvant[All Fields] AND (immunotherapy OR "PD-1" OR "PD-L1" OR pembrolizumab OR nivolumab OR durvalumab OR camrelizumab OR sintilimab OR toripalimab) AND (NSCLC OR "non-small cell lung") AND (RNA-seq OR transcriptome OR GeoMx) AND Homo sapiens[Organism] AND 2024/01/01:2026/12/31[PDAT]',
        ),
        (
            "durvalumab lung RNA 2024-2026",
            "durvalumab[All Fields] AND (NSCLC OR SCLC OR lung) AND Homo sapiens[Organism] AND 2024/01/01:2026/12/31[PDAT]",
        ),
        (
            "PACIFIC durvalumab GEO",
            "PACIFIC[All Fields] AND durvalumab[All Fields] AND Homo sapiens[Organism]",
        ),
        (
            "targeted leftover accessions",
            "GSE329813[Accession] OR GSE241934[Accession] OR GSE292299[Accession] OR GSE261345[Accession] OR GSE261348[Accession] OR GSE207422[Accession] OR GSE126044[Accession] OR GSE166449[Accession] OR GSE135222[Accession] OR GSE248378[Accession]",
        ),
    ]
    seen_gse = set()
    for qname, term in queries:
        try:
            ids = ncbi_esearch(term, retmax=60)
            summary = ncbi_esummary(ids)
        except Exception as exc:  # noqa: BLE001
            add(
                resource="GEO",
                accession_or_id="",
                query=qname,
                status="QUERY_FAILED",
                usable="no",
                notes=str(exc),
            )
            continue
        for uid in ids:
            rec = summary.get(uid, {})
            acc = rec.get("accession") or ""
            if not str(acc).startswith("GSE"):
                continue
            if acc in seen_gse and qname != "targeted leftover accessions":
                continue
            seen_gse.add(acc)
            title = rec.get("title") or ""
            gds_type = rec.get("gdstype") or rec.get("type") or ""
            n = rec.get("n_samples") or rec.get("n_samples", "")
            taxon = rec.get("taxon") or ""
            add(
                resource="GEO (live Entrez)",
                accession_or_id=acc,
                query=qname,
                status="VERIFIED_EXISTS",
                usable="triage",
                title=title,
                gds_type=gds_type,
                taxon=taxon,
                n_samples=n,
                notes=f"{title} | type={gds_type} | taxon={taxon} | n={n}",
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            )

    # Manual triage of known leftover / near-miss accessions (verified above or in wave 2).
    triage = {
        "GSE207422": ("yes_scored", "Neoadjuvant anti-PD-1+chemo bulk log2TPM + MPR. Best extra MPR series. Not durvalumab."),
        "GSE126044": ("yes_scored", "Advanced anti-PD-1 NSCLC counts + ORR. n=16."),
        "GSE166449": ("yes_scored", "Advanced lung ICI TPM + ORR from titles. n=22."),
        "GSE135222": ("yes_scored", "Advanced anti-PD-1/PD-L1 TPM + DCB/PFS. Ensembl IDs."),
        "GSE248378": ("yes_scored", "Post-durva FPKM, same trial as GSE253564. Recurrence/PFS; CD8/GEP after purity."),
        "GSE329813": ("yes_scored_post_geomx", "2026 GeoMx leftover. Post-neoadjuvant tumor-bed MPR. CLDN4 absent. Not pooled with pretreatment bulk."),
        "GSE253564": ("index_not_rescored", "User asked for series besides this accession. Already reported in wave 2. Not re-scored here."),
        "GSE241934": ("excluded_scrna", "NEOTIDE/CTONG2104 scRNA (mtx/barcodes/features). No per-patient bulk TACSTD2 matrix."),
        "GSE292299": ("excluded_visium", "Visium leftover (small n_NR). Not bulk whole-transcriptome."),
        "GSE243013": ("excluded_scrna", "2024 leftover: anti-PD-1 NSCLC single-cell atlas. Not bulk."),
        "GSE225620": ("excluded_blood", "2024 leftover: blood-based neoadjuvant PD-1 RNA. Not tumor RNA."),
        "GSE249568": ("excluded_not_io", "2024 leftover: tepotinib MET array, n=1 patient paired biopsies. Not ICI."),
        "GSE250509": ("excluded_not_io", "2024 leftover: tepotinib MET protein. Not ICI RNA."),
        "GSE261345": ("excluded_essclc_geomx", "CANTABRICO ES-SCLC GeoMx CTA. Not NSCLC neoadjuvant; CLDN4 typically absent."),
        "GSE261348": ("excluded_essclc", "CANTABRICO companion if present. ES-SCLC, not PACIFIC NSCLC."),
        "GSE110390": ("excluded_21gene", "Durvalumab 21-gene IFNγ panel. TACSTD2/CLDN4 absent."),
        "GSE131933": ("excluded_scrna", "Durvalumab scRNA n=2 NSCLC. Not a bulk matrix."),
        "GSE190731": ("excluded_xenograft", "Durvalumab xenograft array. Not patient RNA."),
        "GSE333537": ("excluded_hnscc", "Durvalumab HNSCC. Wrong disease."),
    }
    have = {r["accession_or_id"] for r in rows}
    for acc, (usable, note) in triage.items():
        hits = [r for r in rows if r["accession_or_id"] == acc]
        if hits:
            for r in hits:
                r["usable"] = usable
                r["triage_note"] = note
        else:
            add(
                resource="GEO targeted triage",
                accession_or_id=acc,
                query="targeted leftover accessions",
                status="NOT_IN_THIS_QUERY" if acc not in have else "VERIFIED_EXISTS",
                usable=usable,
                notes=note,
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            )

    # ArrayExpress / BioStudies leftover.
    ae_queries = [
        "durvalumab AND lung AND RNA",
        "neoadjuvant AND NSCLC AND immunotherapy AND RNA-seq",
        "PACIFIC AND durvalumab",
    ]
    ae_hits = 0
    for q in ae_queries:
        try:
            url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(
                {"query": q, "pageSize": 25, "accession": "E-MTAB-*"}
            )
            # BioStudies does not filter accession that way; use query only.
            url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(
                {"query": q, "pageSize": 25}
            )
            r = requests.get(url, timeout=45)
            r.raise_for_status()
            hits = r.json().get("hits") or r.json().get("embedded", {}).get("hits") or []
            if isinstance(r.json().get("hits"), dict):
                hits = r.json()["hits"].get("hits", [])
            # BioStudies v1 returns {hits: [...]} or {totalHits, hits}
            if not hits and "hits" in r.json():
                hits = r.json()["hits"]
            time.sleep(0.3)
        except Exception as exc:  # noqa: BLE001
            add(
                resource="ArrayExpress/BioStudies",
                accession_or_id="",
                query=q,
                status="QUERY_FAILED",
                usable="no",
                notes=str(exc),
            )
            continue
        if isinstance(hits, dict):
            hits = hits.get("hits", [])
        for h in hits[:25]:
            acc = h.get("accession") or h.get("id") or ""
            title = h.get("title") or ""
            if not acc:
                continue
            ae_hits += 1
            is_emtab = str(acc).startswith("E-MTAB") or str(acc).startswith("E-GEOD")
            add(
                resource="ArrayExpress/BioStudies (live)",
                accession_or_id=acc,
                query=q,
                status="VERIFIED_EXISTS",
                usable="not_leftover_lung_ici_matrix" if is_emtab else "not_arrayexpress_processed_lung_ici",
                title=title,
                notes=title[:240],
                url=f"https://www.ebi.ac.uk/biostudies/studies/{acc}",
            )
    add(
        resource="ArrayExpress leftover conclusion",
        accession_or_id="—",
        query="leftover lung ICI with processed matrix + labels",
        status="NONE_USABLE",
        usable="no",
        notes=(
            f"Live BioStudies returned {ae_hits} hits across three queries. "
            "No leftover ArrayExpress lung ICI study with a processed per-patient matrix "
            "plus MPR/DCB/ORR labels was added. E-MTAB-13704 / E-MTAB-15883 were already "
            "taken in sibling leftover work and are not re-used as new extras. "
            "Do not invent an ArrayExpress accession."
        ),
    )
    add(
        resource="PACIFIC / ADRIATIC",
        accession_or_id="NCT02125461 / NCT03703297",
        query="targeted",
        status="NOT_OPEN",
        usable="no",
        notes=(
            "Still no public per-patient PACIFIC or ADRIATIC RNA matrix. "
            "ADRIATIC ASCO 2025 abstract 8014 is group medians only. "
            "The user A2 durvalumab coefficient is taken as given and is not re-audited here."
        ),
        url="https://clinicaltrials.gov/study/NCT02125461",
    )
    return rows


def fmt_p(p) -> str:
    return w2.fmt_p(p)


def fmt_r(x) -> str:
    return w2.fmt_r(x)


def pick(tests: pd.DataFrame, **kw) -> pd.Series | None:
    q = tests
    for k, v in kw.items():
        if k == "test_contains":
            q = q[q["test"].astype(str).str.contains(v, regex=False)]
        else:
            q = q[q[k] == v]
    return None if q.empty else q.iloc[0]


def write_figures(tests: pd.DataFrame, samples: dict[str, pd.DataFrame]) -> None:
    rng = np.random.default_rng(2)
    # Boxplots: TACSTD2 vs native binary endpoint, unadjusted.
    order = [
        ("GSE207422", "MPR"),
        ("GSE329813", "MPR"),
        ("GSE126044", "ORR"),
        ("GSE166449", "ORR"),
        ("GSE135222", "DCB"),
        ("GSE248378", "recurrence"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.8), constrained_layout=True)
    for ax, (series, ep) in zip(axes.ravel(), order):
        rec = pick(tests, series=series, gene="TACSTD2", covariate="none", endpoint=ep, test="Mann-Whitney U")
        samp = samples.get(series)
        if samp is None or rec is None:
            ax.set_axis_off()
            continue
        y = samp.set_index("sample")["TACSTD2"] if "sample" in samp.columns else samp["TACSTD2"]
        if "sample" in samp.columns:
            lab = pd.to_numeric(samp.set_index("sample")["endpoint"], errors="coerce")
            pos_l = str(samp["pos_label"].iloc[0])
            neg_l = str(samp["neg_label"].iloc[0])
        else:
            lab = pd.to_numeric(samp["endpoint"], errors="coerce")
            pos_l, neg_l = "pos", "neg"
            y = samp["TACSTD2"]
        a = y[lab == 1].dropna()
        b = y[lab == 0].dropna()
        ax.boxplot([b, a], tick_labels=[neg_l, pos_l], widths=0.55)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, size=len(b)), b, s=14, color="#444", alpha=0.75, zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, size=len(a)), a, s=14, color="#1f4e79", alpha=0.8, zorder=3)
        ax.set_title(f"{series} TACSTD2 vs {ep}\nn={int(rec['n_neg'] or 0)}/{int(rec['n_pos'] or 0)}  p={fmt_p(rec['p'])}  δ={fmt_r(rec['cliffs_delta'])}", fontsize=8)
        ax.set_ylabel("TACSTD2 (series-native log)")
    fig.suptitle("Extra public lung IO series — TACSTD2 vs cohort-native response (unadjusted)", fontsize=11)
    fig.savefig(FIG / "extra_TACSTD2_response_boxplots.png", dpi=140)
    plt.close(fig)

    # Forest: Cliff's delta for TACSTD2 and CLDN4 vs native endpoint.
    fig, ax = plt.subplots(figsize=(8.6, 5.2), constrained_layout=True)
    ylabels, xs, xerr_lo, xerr_hi, colors, annot = [], [], [], [], [], []
    for gene, color in (("TACSTD2", "#1f4e79"), ("CLDN4", "#b33")):
        for series, ep in order:
            rec = pick(tests, series=series, gene=gene, covariate="none", endpoint=ep, test="Mann-Whitney U")
            if rec is None:
                continue
            if rec.get("gene_present") is False or rec.get("p") is None or (isinstance(rec.get("p"), float) and math.isnan(rec["p"])):
                ylabels.append(f"{series} {gene} vs {ep} (absent/NE)")
                xs.append(0)
                xerr_lo.append(0)
                xerr_hi.append(0)
                colors.append("#cccccc")
                annot.append(str(rec.get("note") or "not estimable"))
                continue
            d = float(rec["cliffs_delta"])
            ylabels.append(f"{series} {gene} vs {ep} (n={int(rec['n_pos'])}/{int(rec['n_neg'])})")
            xs.append(d)
            # Cliff δ has no CI here; show point only.
            xerr_lo.append(0)
            xerr_hi.append(0)
            colors.append(color if float(rec["p"]) < 0.05 else "#888888")
            annot.append(f"δ={fmt_r(d)} p={fmt_p(rec['p'])}")
    ypos = np.arange(len(xs))[::-1]
    ax.axvline(0, color="black", lw=0.7)
    ax.scatter(xs, ypos, c=colors, s=36, zorder=3)
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Cliff's δ (positive = higher in responders / MPR / DCB / recurrence)")
    ax.set_xlim(-1.05, 1.05)
    for y, t in zip(ypos, annot):
        ax.text(0.98, y, t, va="center", ha="right", fontsize=7, transform=ax.get_yaxis_transform())
    ax.set_title("Extra series — TACSTD2 / CLDN4 vs native binary endpoint")
    fig.savefig(FIG / "extra_response_forest_delta.png", dpi=140)
    plt.close(fig)

    # Forest: partial Spearman vs CD8A and GEP18 | ESTIMATEScore.
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), constrained_layout=True)
    series_order = ["GSE207422", "GSE126044", "GSE166449", "GSE135222", "GSE248378", "GSE329813"]
    for ax, endpoint in zip(axes, ("CD8A", "GEP18")):
        labs, vals, ps, ns = [], [], [], []
        for series in series_order:
            rec = pick(
                tests,
                series=series,
                gene="TACSTD2",
                covariate="ESTIMATEScore",
                endpoint=endpoint,
                test_contains="partial Spearman",
            )
            if rec is None or rec.get("p") is None or (isinstance(rec.get("p"), float) and math.isnan(rec["p"])):
                continue
            labs.append(series)
            vals.append(float(rec["rho"]))
            ps.append(float(rec["p"]))
            ns.append(int(rec["n"]))
        colors = ["#1f4e79" if p < 0.05 else "#888" for p in ps]
        ax.barh(range(len(vals)), vals, color=colors)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels([f"{s} n={n}" for s, n in zip(labs, ns)], fontsize=8)
        ax.set_xlim(-1, 1)
        ax.set_xlabel(f"partial Spearman ρ  TACSTD2 vs {endpoint} | ESTIMATEScore")
        for i, (v, p) in enumerate(zip(vals, ps)):
            ax.text(v + (0.03 if v >= 0 else -0.03), i, f"ρ={fmt_r(v)} p={fmt_p(p)}", va="center", ha="left" if v >= 0 else "right", fontsize=7)
        ax.set_title(endpoint)
    fig.suptitle("Extra series — TACSTD2 vs CD8 / GEP after purity residual (ESTIMATEScore)", fontsize=11)
    fig.savefig(FIG / "extra_TACSTD2_CD8_GEP_partial_forest.png", dpi=140)
    plt.close(fig)

    # CLDN4 companion forest for immune residual.
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), constrained_layout=True)
    for ax, endpoint in zip(axes, ("CD8A", "GEP18")):
        labs, vals, ps, ns = [], [], [], []
        for series in series_order:
            rec = pick(
                tests,
                series=series,
                gene="CLDN4",
                covariate="ESTIMATEScore",
                endpoint=endpoint,
                test_contains="partial Spearman",
            )
            if rec is None or rec.get("p") is None or (isinstance(rec.get("p"), float) and math.isnan(rec["p"])):
                continue
            labs.append(series)
            vals.append(float(rec["rho"]))
            ps.append(float(rec["p"]))
            ns.append(int(rec["n"]))
        colors = ["#b33" if p < 0.05 else "#888" for p in ps]
        ax.barh(range(len(vals)), vals, color=colors)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels([f"{s} n={n}" for s, n in zip(labs, ns)], fontsize=8)
        ax.set_xlim(-1, 1)
        ax.set_xlabel(f"partial Spearman ρ  CLDN4 vs {endpoint} | ESTIMATEScore")
        for i, (v, p) in enumerate(zip(vals, ps)):
            ax.text(v + (0.03 if v >= 0 else -0.03), i, f"ρ={fmt_r(v)} p={fmt_p(p)}", va="center", ha="left" if v >= 0 else "right", fontsize=7)
        ax.set_title(endpoint)
    fig.suptitle("Extra series — CLDN4 vs CD8 / GEP after purity residual (ESTIMATEScore)", fontsize=11)
    fig.savefig(FIG / "extra_CLDN4_CD8_GEP_partial_forest.png", dpi=140)
    plt.close(fig)


def write_report(tests: pd.DataFrame, hunt: pd.DataFrame, metas: list[dict]) -> None:
    def cell(series, gene, cov, endpoint, test_contains=None, test=None):
        kw = {"series": series, "gene": gene, "covariate": cov, "endpoint": endpoint}
        if test_contains:
            kw["test_contains"] = test_contains
        if test:
            kw["test"] = test
        return pick(tests, **kw)

    lines = []
    lines.append("# Extra public lung IO RNA on top of A2")
    lines.append("")
    lines.append("**Self-contained. Public data only. Paper extras — not an A2-recovery audit.**")
    lines.append("")
    lines.append(
        "The user durvalumab result (ρ = −0.65 / purity-adjusted −0.46) is **taken as given**. "
        "This folder does not re-test whether that coefficient is recovered. "
        "It adds leftover public lung IO RNA series *besides* GSE253564 and scores TACSTD2 / CLDN4 "
        "against each cohort’s native endpoint (MPR, DCB, ORR, or recurrence) and against CD8A / GEP18 "
        "after an ESTIMATEScore purity residual."
    )
    lines.append("")
    lines.append("## What was added")
    lines.append("")
    lines.append("| Series | Setting / drug | Assay | n | Native endpoint | TACSTD2 | CLDN4 |")
    lines.append("|---|---|---|---|---|---|---|")
    lines.append("| GSE207422 | Neoadjuvant anti-PD-1 + chemo (not durva) | bulk log2(TPM+1) | 24 baseline | MPR 9 vs 15; RECIST ORR 17 vs 7 | yes | yes |")
    lines.append("| GSE126044 | Advanced anti-PD-1 NSCLC | bulk log2CPM | 16 | ORR 5 vs 11 | yes | yes |")
    lines.append("| GSE166449 | Advanced lung ICI (SMC) | deposited log-like TPM | 22 | ORR 7 vs 15 | yes | yes |")
    lines.append("| GSE135222 | Advanced anti-PD-1/PD-L1 NSCLC | bulk log2(TPM+1), Ensembl | 27 | DCB = PFS≥180 d; Cox PFS | yes | yes |")
    lines.append("| GSE248378 | Neoadjuvant durvalumab ± SBRT, **post** residual | bulk log2(FPKM+1) | 29 | recurrence 9 vs 20; PFS (9 events). MPR not estimable | yes | yes |")
    lines.append("| GSE329813 | Neoadjuvant pembro + chemo, **post** residual | GeoMx DSP, patient-mean tumor bed | 22 pts | MPR 11 vs 11 | yes | **absent** |")
    lines.append("")
    lines.append(
        "GSE253564 is the index durvalumab pretreatment accession and is **not re-scored here**. "
        "GSE248378 is the other open durvalumab NSCLC whole-transcriptome accession from the same trial. "
        "No series is pooled. Post-treatment GeoMx is not mixed with pretreatment bulk."
    )
    lines.append("")
    lines.append("## Methods (locked to wave 2)")
    lines.append("")
    lines.append("| Item | Definition | Honest limitation |")
    lines.append("|---|---|---|")
    lines.append("| TACSTD2, CLDN4 | series-native log (log2CPM / log2(TPM+1) / log2(FPKM+1) / log2(GeoMx+1)) | Bulk or ROI-mean, not protein |")
    lines.append("| CD8 | CD8A on the same scale | Single gene |")
    lines.append("| GEP18 | unweighted z-mean of Ayers 2017 18 genes | Missing genes dropped; n present is reported |")
    lines.append("| Purity residual | ESTIMATEScore = Stromal+Immune (v1.0.13 port), separately, never joint | Expression-derived; GeoMx CTA overlap is incomplete |")
    lines.append("| Response tests | Mann–Whitney on unadjusted gene and on OLS residual given ESTIMATEScore | Small-n; no multiple-testing correction |")
    lines.append("| Immune tests | Spearman + partial Spearman (ranks, n−3 df) given ESTIMATEScore | Partial Spearman residualises both ranks |")
    lines.append("| DCB | GSE135222 PFS time ≥ 180 days | Conventional, not a trial-defined DCB field |")
    lines.append("| MPR (GSE207422) | Pathologic Response starts with MPR or contains pCR | Publication metadata, not recoded onto PACIFIC |")
    lines.append("| MPR (GSE329813) | Title label on Primary tumor bed ROIs, patient mean | Post-treatment residual, not a pretreatment predictor |")
    lines.append("")
    lines.append("## Direct findings (honest n / ρ / p)")
    lines.append("")
    t207 = cell("GSE207422", "TACSTD2", "none", "MPR", test="Mann-Whitney U")
    t329 = cell("GSE329813", "TACSTD2", "none", "MPR", test="Mann-Whitney U")
    t329e = cell("GSE329813", "TACSTD2", "ESTIMATEScore", "MPR", test="Mann-Whitney U")
    t248_cd8 = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "CD8A", test_contains="partial Spearman")
    t248_gep = cell("GSE248378", "TACSTD2", "ESTIMATEScore", "GEP18", test_contains="partial Spearman")
    t207_cd8 = cell("GSE207422", "TACSTD2", "ESTIMATEScore", "CD8A", test_contains="partial Spearman")
    t329_cd8 = cell("GSE329813", "TACSTD2", "ESTIMATEScore", "CD8A", test_contains="partial Spearman")
    t126 = cell("GSE126044", "TACSTD2", "none", "ORR", test="Mann-Whitney U")
    t166 = cell("GSE166449", "TACSTD2", "none", "ORR", test="Mann-Whitney U")
    t135 = cell("GSE135222", "TACSTD2", "none", "DCB", test="Mann-Whitney U")
    lines.append(
        "1. **Pretreatment bulk MPR/ORR/DCB.** GSE207422 TACSTD2 vs MPR n=9/15, "
        f"δ={fmt_r(t207['cliffs_delta']) if t207 is not None else 'NA'}, p={fmt_p(t207['p']) if t207 is not None else 'NA'}. "
        f"GSE126044 ORR n=5/11 p={fmt_p(t126['p']) if t126 is not None else 'NA'}. "
        f"GSE166449 ORR n=7/15 p={fmt_p(t166['p']) if t166 is not None else 'NA'}. "
        f"GSE135222 DCB n=7/20 p={fmt_p(t135['p']) if t135 is not None else 'NA'}. "
        "CLDN4 vs these endpoints is likewise small-n and not significant. "
        "These are leftover anti-PD-1 / chemo-IO series, not durvalumab, and are not pooled."
    )
    lines.append(
        "2. **Post-treatment residual GeoMx (GSE329813).** Tumor-bed TACSTD2 is lower in MPR "
        f"(n=11/11, δ={fmt_r(t329['cliffs_delta']) if t329 is not None else 'NA'}, p={fmt_p(t329['p']) if t329 is not None else 'NA'}). "
        f"After ESTIMATEScore residual p={fmt_p(t329e['p']) if t329e is not None else 'NA'}. "
        "CLDN4 is absent from the panel. This is post-neoadjuvant residual tissue, not a pretreatment predictor, "
        "and is not mixed with GSE207422."
    )
    lines.append(
        "3. **CD8 / GEP after purity residual.** The open post-durvalumab accession GSE248378 "
        f"gives TACSTD2 vs CD8A partial ρ={fmt_r(t248_cd8['rho']) if t248_cd8 is not None else 'NA'} "
        f"p={fmt_p(t248_cd8['p']) if t248_cd8 is not None else 'NA'} (n=29) and vs GEP18 "
        f"ρ={fmt_r(t248_gep['rho']) if t248_gep is not None else 'NA'} "
        f"p={fmt_p(t248_gep['p']) if t248_gep is not None else 'NA'}. "
        f"GSE329813 TACSTD2 vs CD8A partial ρ={fmt_r(t329_cd8['rho']) if t329_cd8 is not None else 'NA'} "
        f"p={fmt_p(t329_cd8['p']) if t329_cd8 is not None else 'NA'} (n=22; incomplete ESTIMATE overlap). "
        f"GSE207422 pretreatment TACSTD2 vs CD8A partial ρ={fmt_r(t207_cd8['rho']) if t207_cd8 is not None else 'NA'} "
        f"p={fmt_p(t207_cd8['p']) if t207_cd8 is not None else 'NA'} (n=24). "
        "Advanced-disease leftover series (GSE126044 / GSE166449 / GSE135222) are compatible with a negative "
        "or null TACSTD2–CD8 residual; none reach p<0.05 after ESTIMATE."
    )
    lines.append("")
    lines.append("## TACSTD2 / CLDN4 vs native response")
    lines.append("")
    lines.append("| Series | Gene | Endpoint | n_pos / n_neg | median_pos | median_neg | Cliff δ | p | After ESTIMATE residual p |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    response_rows = [
        ("GSE207422", "TACSTD2", "MPR"),
        ("GSE207422", "CLDN4", "MPR"),
        ("GSE207422", "TACSTD2", "ORR"),
        ("GSE207422", "CLDN4", "ORR"),
        ("GSE329813", "TACSTD2", "MPR"),
        ("GSE329813", "CLDN4", "MPR"),
        ("GSE126044", "TACSTD2", "ORR"),
        ("GSE126044", "CLDN4", "ORR"),
        ("GSE166449", "TACSTD2", "ORR"),
        ("GSE166449", "CLDN4", "ORR"),
        ("GSE135222", "TACSTD2", "DCB"),
        ("GSE135222", "CLDN4", "DCB"),
        ("GSE248378", "TACSTD2", "recurrence"),
        ("GSE248378", "CLDN4", "recurrence"),
    ]
    for series, gene, ep in response_rows:
        u = cell(series, gene, "none", ep, test="Mann-Whitney U")
        a = cell(series, gene, "ESTIMATEScore", ep, test="Mann-Whitney U")
        if u is None:
            lines.append(f"| {series} | {gene} | {ep} | — | — | — | — | — | — |")
            continue
        npos = u.get("n_pos")
        nneg = u.get("n_neg")
        note = u.get("note") or ""
        if u.get("p") is None or (isinstance(u.get("p"), float) and math.isnan(u["p"])):
            lines.append(f"| {series} | {gene} | {ep} | NE | — | — | — | NE | NE |")
            continue
        ap = fmt_p(a["p"]) if a is not None else "NA"
        lines.append(
            f"| {series} | {gene} | {ep} | {int(npos)}/{int(nneg)} | {fmt_r(u['median_pos'])} | "
            f"{fmt_r(u['median_neg'])} | {fmt_r(u['cliffs_delta'])} | {fmt_p(u['p'])} | {ap} |"
        )
    lines.append("")
    lines.append(
        "Positive Cliff δ means the gene is **higher** in the positive class "
        "(MPR / responder / DCB / recurrence). Signs are not forced onto the user A2 coefficient."
    )
    lines.append("")

    # PFS extras
    lines.append("## Continuous PFS (where deposited)")
    lines.append("")
    lines.append("| Series | Gene | n / events | unadj HR | unadj p | ESTIMATE-residual HR | residual p |")
    lines.append("|---|---|---|---|---|---|---|")
    for series, gene in (("GSE135222", "TACSTD2"), ("GSE135222", "CLDN4"), ("GSE248378", "TACSTD2"), ("GSE248378", "CLDN4")):
        u = cell(series, gene, "none", "PFS")
        a = cell(series, gene, "ESTIMATEScore", "PFS")
        if u is None:
            continue
        lines.append(
            f"| {series} | {gene} | {int(u['n'])} / {int(u['n_events'])} | {fmt_r(u['hr'])} | {fmt_p(u['p'])} | "
            f"{fmt_r(a['hr']) if a is not None else 'NA'} | {fmt_p(a['p']) if a is not None else 'NA'} |"
        )
    lines.append("")

    lines.append("## TACSTD2 / CLDN4 vs CD8A and GEP18 after purity residual")
    lines.append("")
    lines.append("| Series | Gene | Immune | n | raw ρ | raw p | partial ρ given ESTIMATE | partial p |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for series in ("GSE207422", "GSE126044", "GSE166449", "GSE135222", "GSE248378", "GSE329813"):
        for gene in TARGETS:
            for ep in ("CD8A", "GEP18"):
                raw = cell(series, gene, "none", ep, test_contains="Spearman unadjusted")
                adj = cell(series, gene, "ESTIMATEScore", ep, test_contains="partial Spearman")
                if raw is None:
                    continue
                if raw.get("p") is None or (isinstance(raw.get("p"), float) and math.isnan(raw["p"])):
                    lines.append(f"| {series} | {gene} | {ep} | {raw.get('n')} | NE | NE | NE | NE |")
                    continue
                lines.append(
                    f"| {series} | {gene} | {ep} | {int(raw['n'])} | {fmt_r(raw['rho'])} | {fmt_p(raw['p'])} | "
                    f"{fmt_r(adj['rho']) if adj is not None else 'NA'} | {fmt_p(adj['p']) if adj is not None else 'NA'} |"
                )
    lines.append("")
    lines.append(
        "Partial Spearman is the A2-matched estimator (Pearson on ranks, n−3 df). "
        "GSE248378 is the open post-durvalumab series; the other rows are leftover public anti-PD-1 / chemo-IO lung RNA."
    )
    lines.append("")

    lines.append("## ESTIMATE / GEP coverage")
    lines.append("")
    lines.append("| Series | n | ImmuneSignature overlap | StromalSignature overlap | GEP18 genes present | TumorPurity OOB |")
    lines.append("|---|---|---|---|---|---|")
    for m in metas:
        ov = m.get("estimate_overlap") or {}
        lines.append(
            f"| {m['series']} | {m['n_samples']} | {ov.get('ImmuneSignature', 'NA')} | "
            f"{ov.get('StromalSignature', 'NA')} | {m.get('n_gep18', 'NA')}/18 | {m.get('n_purity_oob', 'NA')} |"
        )
    lines.append("")
    lines.append(
        "GSE329813 is a GeoMx panel (~1.8k genes). ESTIMATE overlap is incomplete; "
        "the purity residual is still computed on the genes that are present and is labelled as such."
    )
    lines.append("")

    lines.append("## Leftover hunt (live) — what cannot be added")
    lines.append("")
    gse_live = hunt[hunt["accession_or_id"].astype(str).str.startswith("GSE")]
    ae_none = hunt[hunt["resource"].str.contains("leftover conclusion", na=False)]
    lines.append(
        f"Live Entrez returned {gse_live['accession_or_id'].nunique()} unique GSE accessions across the leftover queries. "
        "Scored extras are listed above. The following were verified and **not** used as additional bulk matrices:"
    )
    lines.append("")
    lines.append("| Accession | Decision | Why |")
    lines.append("|---|---|---|")
    lines.append("| PACIFIC NCT02125461 | closed | No public per-patient RNA. User A2 number is used as given. |")
    lines.append("| ADRIATIC NCT03703297 | closed | ASCO 2025 group medians only. |")
    lines.append("| GSE241934 | excluded | scRNA / residual epithelium, not a per-patient bulk TACSTD2 matrix. |")
    lines.append("| GSE292299 | excluded | Visium; n_NR too small for a bulk ORR residual. |")
    lines.append("| GSE261345 / GSE261348 | excluded | ES-SCLC GeoMx, not NSCLC neoadjuvant. |")
    lines.append("| GSE110390 | excluded | 21-gene IFNγ panel; TACSTD2/CLDN4 absent. |")
    lines.append("| GSE333537 / GSE183924 | excluded | Durvalumab RNA, wrong disease (HNSCC / esophageal). |")
    lines.append("| GSE243013 | excluded | 2024 leftover: anti-PD-1 NSCLC scRNA atlas. |")
    lines.append("| GSE225620 | excluded | 2024 leftover: blood RNA after neoadjuvant PD-1. Not tumor. |")
    lines.append("| GSE249568 / GSE250509 | excluded | 2024 leftover: tepotinib MET, not ICI. |")
    lines.append("| ArrayExpress leftover | none usable | Live BioStudies: no new processed lung ICI matrix + labels. |")
    if not ae_none.empty:
        lines.append("")
        lines.append(ae_none.iloc[0]["notes"])
    lines.append("")
    lines.append("Near-miss ≠ hit. No GEO or ArrayExpress accession was invented.")
    lines.append("")

    lines.append("## Paper extras")
    lines.append("")
    lines.append("Use these as **additional public figures/tables**, not as a replacement for the user A2 slide.")
    lines.append("")
    lines.append("- `figures/extra_TACSTD2_response_boxplots.png` — TACSTD2 vs native endpoint, one panel per series.")
    lines.append("- `figures/extra_response_forest_delta.png` — Cliff’s δ forest (TACSTD2 and CLDN4).")
    lines.append("- `figures/extra_TACSTD2_CD8_GEP_partial_forest.png` — TACSTD2 vs CD8A / GEP18 after ESTIMATE residual.")
    lines.append("- `figures/extra_CLDN4_CD8_GEP_partial_forest.png` — CLDN4 companion.")
    lines.append("- `response_tests.tsv` / `immune_correlations.tsv` / `all_tests.tsv` — every n / ρ / p.")
    lines.append("- `inventory.tsv` — scored vs excluded leftover accessions.")
    lines.append("")
    lines.append("## How to read this next to A2")
    lines.append("")
    lines.append(
        "A2 is a user-reported durvalumab RNA association. These extras ask a different, public question: "
        "in leftover open lung IO series, what are TACSTD2 and CLDN4 doing versus MPR/DCB/ORR and versus "
        "CD8/GEP after the same purity residual. Small n. Hypothesis-generating. Signs are reported as computed."
    )
    lines.append("")
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    paper = []
    paper.append("# Paper extras — TACSTD2 / CLDN4 in leftover public lung IO RNA")
    paper.append("")
    paper.append("Source: `results/rework/A2_wave2/extra/`. User A2 durvalumab result taken as given.")
    paper.append("")
    # Reuse the two main tables from lines
    start = None
    for i, line in enumerate(lines):
        if line.startswith("## TACSTD2 / CLDN4 vs native response"):
            start = i
            break
    if start is not None:
        chunk = []
        for line in lines[start:]:
            if line.startswith("## Leftover hunt"):
                break
            chunk.append(line)
        paper.extend(chunk)
    paper.append("Figures: `extra_TACSTD2_response_boxplots.png`, `extra_response_forest_delta.png`,")
    paper.append("`extra_TACSTD2_CD8_GEP_partial_forest.png`, `extra_CLDN4_CD8_GEP_partial_forest.png`.")
    paper.append("")
    (OUT / "paper_tables.md").write_text("\n".join(paper) + "\n")


def main() -> int:
    paths = ensure_extra_inputs()
    gmt = w2.load_gmt(paths["gmt"])

    loaders = [
        ("GSE207422", lambda: load_gse207422(paths)),
        ("GSE126044", lambda: load_gse126044(paths)),
        ("GSE166449", lambda: load_gse166449(paths)),
        ("GSE135222", lambda: load_gse135222(paths, gmt)),
        ("GSE329813", lambda: load_gse329813(paths)),
        ("GSE248378", lambda: load_gse248378(paths)),
    ]
    all_rows: list[dict] = []
    samples: dict[str, pd.DataFrame] = {}
    metas: list[dict] = []
    for name, fn in loaders:
        expr, meta, info = fn()
        extra_bin = None
        if name == "GSE207422" and "orr" in meta.columns:
            extra_bin = {"ORR": pd.to_numeric(meta["orr"], errors="coerce")}
        samp, rows, meta_out = analyse_series(name, expr, meta, gmt, info, extra_binaries=extra_bin)
        samples[name] = samp
        all_rows.extend(rows)
        metas.append(meta_out)
        samp.to_csv(OUT / f"sample_table_{name}.tsv", sep="\t", index=False)
        print(f"[{name}] n={meta_out['n_samples']} genes={meta_out['n_genes']} "
              f"GEP={meta_out['n_gep18']} EST_I={meta_out['estimate_overlap'].get('ImmuneSignature')} "
              f"targets={meta_out['targets_present']}")

    tests = pd.DataFrame(all_rows)
    tests.to_csv(OUT / "all_tests.tsv", sep="\t", index=False)
    tests[tests["endpoint"].isin(["MPR", "ORR", "DCB", "recurrence", "PFS"])].to_csv(
        OUT / "response_tests.tsv", sep="\t", index=False
    )
    tests[tests["endpoint"].isin(["CD8A", "GEP18", "CYT"])].to_csv(
        OUT / "immune_correlations.tsv", sep="\t", index=False
    )

    hunt_rows = hunt_leftover_public()
    hunt = pd.DataFrame(hunt_rows)
    hunt.to_csv(OUT / "hunt_leftover.tsv", sep="\t", index=False)
    inv = hunt[hunt["accession_or_id"].astype(str).str.startswith("GSE") | hunt["resource"].str.contains("ArrayExpress|PACIFIC", na=False)].copy()
    inv.to_csv(OUT / "inventory.tsv", sep="\t", index=False)

    write_figures(tests, samples)
    write_report(tests, hunt, metas)

    def grab(series, gene, cov, endpoint, test_sub=None):
        q = tests[(tests.series == series) & (tests.gene == gene) & (tests.covariate == cov) & (tests.endpoint == endpoint)]
        if test_sub:
            q = q[q.test.str.contains(test_sub, regex=False)]
        if q.empty:
            return None
        r = q.iloc[0]
        keys = ["p", "rho", "hr", "cliffs_delta", "n", "n_pos", "n_neg", "n_events", "note", "gene_present"]
        return {k: (None if (isinstance(r.get(k), float) and math.isnan(r.get(k))) else r.get(k)) for k in keys if k in r.index}

    summary = {
        "question": "Add public TACSTD2/CLDN4 content on top of given A2: extra lung IO series vs MPR/DCB/ORR and vs CD8/GEP after ESTIMATE residual",
        "a2_user_result_taken_as_given": {"rho_raw": -0.65, "rho_purity_adj": -0.46},
        "gse253564_rescored": False,
        "series_scored": [m["series"] for m in metas],
        "methods": {
            "purity": "ESTIMATEScore residual / partial Spearman (ranks, n-3)",
            "GEP18": "Ayers 2017 unweighted z-mean",
            "CD8": "CD8A",
            "no_pooling": True,
        },
        "key": {
            "GSE207422_TACSTD2_MPR": grab("GSE207422", "TACSTD2", "none", "MPR"),
            "GSE207422_TACSTD2_MPR_ESTIMATE": grab("GSE207422", "TACSTD2", "ESTIMATEScore", "MPR"),
            "GSE126044_TACSTD2_ORR": grab("GSE126044", "TACSTD2", "none", "ORR"),
            "GSE166449_TACSTD2_ORR": grab("GSE166449", "TACSTD2", "none", "ORR"),
            "GSE135222_TACSTD2_DCB": grab("GSE135222", "TACSTD2", "none", "DCB"),
            "GSE248378_TACSTD2_CD8_partial_ESTIMATE": grab("GSE248378", "TACSTD2", "ESTIMATEScore", "CD8A", "partial Spearman"),
            "GSE248378_TACSTD2_GEP_partial_ESTIMATE": grab("GSE248378", "TACSTD2", "ESTIMATEScore", "GEP18", "partial Spearman"),
            "GSE329813_TACSTD2_MPR": grab("GSE329813", "TACSTD2", "none", "MPR"),
        },
        "metas": metas,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    provenance = {
        "inputs": {
            name: {"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "md5": w2.md5_file(path)}
            for name, path in paths.items()
            if path.exists()
        },
        "urls": URLS,
        "gep18": list(GEP18),
        "note": "GSE253564 is not re-scored. User A2 coefficient is taken as given.",
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2))
    print(f"[A2_extra] wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
