#!/usr/bin/env python3
"""GSE81089 Swedish NSCLC RNA-seq: CLDN4 vs immune / IFN-γ / MHC-I.

Additive CLDN4-only slice. No TACSTD2∩CLDN4 dual-high. Unit is the
patient / array (one tumor GSM per patient). Histology codes are the
GEO overall_design integers (1=squamous, 2=AC unspecified, 3=large
cell/NOS). Series n=218 includes 19 paired normals and is not the
tumor or histology n.

Primary matrix is Cufflinks FPKM as deposited (log2(FPKM+1)).
featureCounts log2(CPM+1) is a same-run sensitivity, not a FASTQ
re-quant.

ImmuneScore = mean of gene-wise z-scores of the Yoshihara 2013
ESTIMATE ImmuneSignature. This is not the official R ssGSEA score.

Downloads stay under $GSE81089_CLDN4_DATA (default /tmp/gse81089_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import json
import os
import tarfile
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE81089_CLDN4_DATA", "/tmp/gse81089_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/"
    "matrix/GSE81089_series_matrix.txt.gz"
)
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/"
    "suppl/GSE81089_FPKM_cufflinks.tsv.gz"
)
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/"
    "suppl/GSE81089_readcounts_featurecounts.tsv.gz"
)
GTF_URL = (
    "https://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh37.73.gtf.gz"
)
ESTIMATE_TAR = "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"

# Ensembl v73 / GRCh37 IDs confirmed on the deposited FPKM matrix.
NAMED_ENSG = {
    "CLDN4": "ENSG00000189143",
    "TACSTD2": "ENSG00000184292",
    "CD8A": "ENSG00000153563",
    "CD274": "ENSG00000120217",
    "EPCAM": "ENSG00000119888",
    "CDH1": "ENSG00000039068",
    "KRT8": "ENSG00000170421",
    "KRT18": "ENSG00000111057",
    "KRT19": "ENSG00000171345",
    "KRT7": "ENSG00000135480",
    "IFNG": "ENSG00000111537",
    "STAT1": "ENSG00000115415",
    "CXCL9": "ENSG00000138755",
    "CXCL10": "ENSG00000169245",
    "IDO1": "ENSG00000131203",
    "HLA-DRA": "ENSG00000204287",
    "HLA-A": "ENSG00000206503",
    "HLA-B": "ENSG00000234745",
    "HLA-C": "ENSG00000204525",
    "B2M": "ENSG00000166710",
    "TAP1": "ENSG00000168394",
    "TAP2": "ENSG00000204267",
}

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
IFNG_AYERS = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]

# GEO overall_design (series record).
HISTOLOGY_MAP = {
    "1": "squamous_cell_cancer",
    "2": "AC_unspecified",
    "3": "large_cell_NOS",
}
STAGE_MAP = {
    "1": "1a",
    "2": "1b",
    "3": "2a",
    "4": "2b",
    "5": "3a",
    "6": "3b",
    "7": "IV",
}
SMOKING_MAP = {"1": "current", "2": "ex_gt_1year", "3": "never"}

COL_ALIAS = {"L608T_2122": "L608T", "L771T_1": "L771T"}


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse81089-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    return dest


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields.setdefault(key, []).append(vals)
    n = len(sample_fields["geo_accession"][0])
    rows = []
    for i in range(n):
        rec = {
            "gsm": sample_fields["geo_accession"][0][i],
            "title": sample_fields["title"][0][i],
            "source_name": sample_fields.get("source_name_ch1", [[""] * n])[0][i],
        }
        for block in sample_fields.get("characteristics_ch1", []):
            v = block[i]
            if ": " in v:
                k, val = v.split(": ", 1)
                rec[k.strip()] = val.strip()
            elif v:
                rec.setdefault("unparsed", v)
        rows.append(rec)
    meta = pd.DataFrame(rows).set_index("gsm")
    meta.index.name = "gsm"
    return meta, {k: " | ".join(v) for k, v in series.items()}


def parse_gtf_symbols(path: Path) -> pd.DataFrame:
    """Ensembl v73 GTF has exon/CDS rows only (no gene features)."""
    seen: dict[str, tuple[str, str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            gid = gname = biotype = ""
            for field in attr.strip().strip(";").split(";"):
                field = field.strip()
                if not field:
                    continue
                key, _, val = field.partition(" ")
                val = val.strip().strip('"')
                if key == "gene_id":
                    gid = val.split(".")[0]
                elif key == "gene_name":
                    gname = val
                elif key == "gene_biotype":
                    biotype = val
            if gid and gid not in seen:
                seen[gid] = (gname, biotype)
    return pd.DataFrame(
        [{"ensembl": k, "symbol": v[0], "biotype": v[1]} for k, v in seen.items()]
    )


def load_estimate_immune_genes(tar_path: Path) -> list[str]:
    with tarfile.open(tar_path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith("SI_geneset.gmt"))
        raw = tf.extractfile(member).read().decode("utf-8")
    immune = None
    for line in raw.splitlines():
        parts = line.strip().split("\t")
        if parts and parts[0] == "ImmuneSignature":
            immune = [g for g in parts[2:] if g]
    if not immune:
        raise RuntimeError("ImmuneSignature not found in SI_geneset.gmt")
    if len(immune) != 141:
        raise RuntimeError(f"expected 141 ImmuneSignature genes, got {len(immune)}")
    return immune


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 8:
        return dict(n=n, rho=np.nan, p=np.nan, ci_low=np.nan, ci_high=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            boots[i] = np.nan
            continue
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return dict(n=n, rho=float(rho), p=float(p), ci_low=float(lo), ci_high=float(hi))


def partial_spearman(x, y, covar):
    x = pd.Series(np.asarray(x, float))
    y = pd.Series(np.asarray(y, float))
    c = pd.Series(np.asarray(covar, float), name="cov")
    df = pd.concat([x.rename("x"), y.rename("y"), c], axis=1).dropna()
    n = int(df.shape[0])
    if n < 12:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(df["x"].values)
    yr = stats.rankdata(df["y"].values)
    cr = np.column_stack([np.ones(n), stats.rankdata(df["cov"].values)])
    bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
    rx = xr - cr @ bx
    ry = yr - cr @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = n - 3
    t = r * np.sqrt(dof / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def verdict(rho, p, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if not np.isfinite(rho) or not np.isfinite(p):
        return "NO_EVIDENCE"
    if rho < 0 and p < 0.05:
        return "HOLDS"
    if rho > 0 and p < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def score_mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(dtype=float), present
    z = zscore_rows(expr.loc[present])
    return z.mean(axis=0), present


def collapse_symbol(fpkm: pd.DataFrame, ens2sym: pd.Series) -> pd.DataFrame:
    """Max-mean collapse Ensembl rows that share a gene symbol."""
    df = fpkm.copy()
    df["_gene"] = df.index.map(lambda x: ens2sym.get(str(x), ""))
    df = df[df["_gene"].astype(str).str.len() > 0]
    means = df.drop(columns="_gene").astype(float).mean(axis=1)
    df["_mean"] = means
    df = df.sort_values("_mean", ascending=False)
    kept = df.loc[~df["_gene"].duplicated(keep="first")]
    out = kept.drop(columns=["_gene", "_mean"]).astype(float)
    out.index = kept["_gene"].values
    return out.sort_index()


def pair_row(name, x, y, epithelial, cohort, extra=None, matrix="FPKM_log2p1"):
    crude = spearman_ci(x.values, y.values)
    if epithelial is not None and len(epithelial):
        adj = partial_spearman(x.values, y.values, epithelial.reindex(x.index).values)
    else:
        adj = dict(n=crude["n"], rho_adj=np.nan, p_adj=np.nan)
    rec = {
        "cohort": cohort,
        "matrix": matrix,
        "pair": name,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "verdict_crude": verdict(crude["rho"], crude["p"], crude["n"]),
        "verdict_adj": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
    }
    if extra:
        rec.update(extra)
    return rec


def highlow(x: pd.Series, y: pd.Series, endpoint: str, cohort: str, matrix="FPKM_log2p1") -> dict:
    q = pd.qcut(x, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = y[q == "Q4"]
    b = y[q == "Q1"]
    if a.size < 5 or b.size < 5:
        return {
            "cohort": cohort,
            "matrix": matrix,
            "endpoint": endpoint,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
        }
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "cohort": cohort,
        "matrix": matrix,
        "endpoint": endpoint,
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
    }


def analyze_cohort(name, samples, gene_expr, immune_genes, matrix="FPKM_log2p1"):
    sub = gene_expr.loc[:, samples]
    cldn4 = sub.loc["CLDN4"]
    epi, epi_used = score_mean_z(sub, EPITHELIAL)
    immune, imm_used = score_mean_z(sub, immune_genes)
    ifng, ifn_used = score_mean_z(sub, IFNG_AYERS)
    mhc, mhc_used = score_mean_z(sub, MHC1)
    rows = []

    def add_gene(g, note="primary"):
        extra = {"n_genes_in_score": 1, "genes_used": g, "present": g in sub.index, "note": note}
        if g not in sub.index:
            rows.append(
                {
                    "cohort": name,
                    "matrix": matrix,
                    "pair": f"CLDN4 vs {g}",
                    "n": 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "rho_adj_epithelial": np.nan,
                    "p_adj_epithelial": np.nan,
                    "verdict_crude": "ABSENT",
                    "verdict_adj": "ABSENT",
                    **extra,
                }
            )
            return
        rows.append(pair_row(f"CLDN4 vs {g}", cldn4, sub.loc[g], epi if len(epi) else None, name, extra, matrix))

    add_gene("CD8A", "primary")
    add_gene("CD274", "primary")
    rows.append(
        pair_row(
            "CLDN4 vs ImmuneScore",
            cldn4,
            immune,
            epi if len(epi) else None,
            name,
            {"n_genes_in_score": len(imm_used), "genes_used": ",".join(imm_used), "present": True, "note": "primary"},
            matrix,
        )
    )
    rows.append(
        pair_row(
            "CLDN4 vs IFNg_Ayers",
            cldn4,
            ifng,
            epi if len(epi) else None,
            name,
            {"n_genes_in_score": len(ifn_used), "genes_used": ",".join(ifn_used), "present": True, "note": "primary"},
            matrix,
        )
    )
    rows.append(
        pair_row(
            "CLDN4 vs MHC-I",
            cldn4,
            mhc,
            epi if len(epi) else None,
            name,
            {"n_genes_in_score": len(mhc_used), "genes_used": ",".join(mhc_used), "present": True, "note": "primary"},
            matrix,
        )
    )
    if "TACSTD2" in sub.index:
        rows.append(
            pair_row(
                "CLDN4 vs TACSTD2",
                cldn4,
                sub.loc["TACSTD2"],
                epi if len(epi) else None,
                name,
                {"n_genes_in_score": 1, "genes_used": "TACSTD2", "present": True, "note": "companion"},
                matrix,
            )
        )
    if len(epi):
        rows.append(
            pair_row(
                "CLDN4 vs epithelial mean-z",
                cldn4,
                epi,
                None,
                name,
                {"n_genes_in_score": len(epi_used), "genes_used": ",".join(epi_used), "present": True, "note": "companion"},
                matrix,
            )
        )
    hl = []
    for endpoint, series in [
        ("CD8A", sub.loc["CD8A"] if "CD8A" in sub.index else None),
        ("CD274", sub.loc["CD274"] if "CD274" in sub.index else None),
        ("ImmuneScore", immune if len(immune) else None),
        ("IFNg_Ayers", ifng if len(ifng) else None),
        ("MHC-I", mhc if len(mhc) else None),
    ]:
        if series is not None and len(series):
            hl.append(highlow(cldn4, series, endpoint, name, matrix))
    info = {
        "cohort": name,
        "matrix": matrix,
        "n": int(sub.shape[1]),
        "immune_genes_used": len(imm_used),
        "immune_genes_listed": len(immune_genes),
        "ifng_genes_used": len(ifn_used),
        "mhc_genes_used": len(mhc_used),
        "epithelial_genes_used": len(epi_used),
        "cldn4_finite": int(cldn4.notna().sum()),
    }
    scores = {
        "CLDN4": cldn4,
        "CD8A": sub.loc["CD8A"] if "CD8A" in sub.index else pd.Series(dtype=float),
        "CD274": sub.loc["CD274"] if "CD274" in sub.index else pd.Series(dtype=float),
        "ImmuneScore": immune,
        "IFNg_Ayers": ifng,
        "MHC-I": mhc,
        "epithelial": epi,
        "TACSTD2": sub.loc["TACSTD2"] if "TACSTD2" in sub.index else pd.Series(dtype=float),
    }
    return rows, hl, info, scores


def clean_json(obj):
    if isinstance(obj, dict):
        return {k: clean_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_json(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return None if not np.isfinite(v) else v
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, pd.Series):
        return None
    return obj


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    # Reuse already-fetched copies when present (same filenames).
    for src, name in [
        (Path("/tmp/gse81089/GSE81089_series_matrix.txt.gz"), "GSE81089_series_matrix.txt.gz"),
        (Path("/tmp/gse81089/GSE81089_FPKM_cufflinks.tsv.gz"), "GSE81089_FPKM_cufflinks.tsv.gz"),
        (Path("/tmp/gse81089/GSE81089_readcounts_featurecounts.tsv.gz"), "GSE81089_readcounts_featurecounts.tsv.gz"),
        (Path("/tmp/gse81089/Homo_sapiens.GRCh37.73.gtf.gz"), "Homo_sapiens.GRCh37.73.gtf.gz"),
        (Path("/tmp/gse81089/estimate_1.0.13.tar.gz"), "estimate_1.0.13.tar.gz"),
    ]:
        dest = DATA / name
        if src.exists() and not dest.exists():
            dest.write_bytes(src.read_bytes())

    matrix_path = dl(DATA / "GSE81089_series_matrix.txt.gz", GEO_MATRIX)
    fpkm_path = dl(DATA / "GSE81089_FPKM_cufflinks.tsv.gz", FPKM_URL)
    counts_path = dl(DATA / "GSE81089_readcounts_featurecounts.tsv.gz", COUNTS_URL)
    gtf_path = dl(DATA / "Homo_sapiens.GRCh37.73.gtf.gz", GTF_URL)
    est_tar = dl(DATA / "estimate_1.0.13.tar.gz", ESTIMATE_TAR)

    immune_genes = load_estimate_immune_genes(est_tar)
    (DATA / "SI_geneset.ImmuneSignature.txt").write_text("\n".join(immune_genes) + "\n")

    meta, series = parse_series_meta(matrix_path)
    meta["sample_id"] = (
        meta["title"]
        .str.replace(r"^matched sample_", "", regex=True)
        .str.replace(r"^matched sample ", "", regex=True)
    )
    meta["tissue"] = np.where(
        meta["source_name"].str.contains("non-malignant", case=False, na=False)
        | meta["sample_id"].str.endswith("N"),
        "normal",
        "tumor",
    )
    hist_raw = meta.get("histology", pd.Series(index=meta.index, dtype=str)).astype(str)
    meta["histology_code"] = hist_raw.where(hist_raw.isin(["1", "2", "3"]), "")
    meta["histology_label"] = meta["histology_code"].map(HISTOLOGY_MAP).fillna("")
    stage_raw = meta.get("stage tnm", meta.get("stage_tnm", pd.Series(index=meta.index, dtype=str))).astype(str)
    meta["stage_code"] = stage_raw
    meta["stage_label"] = stage_raw.map(STAGE_MAP).fillna("")
    smoke_raw = meta.get("smoking", pd.Series(index=meta.index, dtype=str)).astype(str)
    meta["smoking_label"] = smoke_raw.map(SMOKING_MAP).fillna("")

    gtf = parse_gtf_symbols(gtf_path)
    ens2sym = gtf.set_index("ensembl")["symbol"]

    fpkm = pd.read_csv(fpkm_path, sep="\t", index_col=0)
    fpkm.index = fpkm.index.astype(str).str.split(".").str[0]
    fpkm.columns = [COL_ALIAS.get(c, c) for c in fpkm.columns]
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str).str.split(".").str[0]
    counts.columns = [COL_ALIAS.get(c, c) for c in counts.columns]

    # Align columns to GEO sample_id (patient/array unit).
    sid = meta["sample_id"].astype(str)
    if sid.duplicated().any():
        raise RuntimeError(f"duplicate sample_id: {sid[sid.duplicated()].tolist()}")
    id2gsm = pd.Series(meta.index, index=sid)
    shared = [c for c in fpkm.columns if c in id2gsm.index]
    if len(shared) != 218:
        raise RuntimeError(f"expected 218 FPKM columns matched to GEO, got {len(shared)}")
    fpkm = fpkm.loc[:, shared]
    fpkm.columns = id2gsm.loc[shared].values
    counts = counts.loc[:, shared]
    counts.columns = id2gsm.loc[shared].values
    meta = meta.loc[fpkm.columns]

    # Named-gene confirmation on deposited Ensembl IDs.
    confirm = []
    for gene, eid in NAMED_ENSG.items():
        confirm.append(
            {
                "gene": gene,
                "ensembl_v73": eid,
                "in_fpkm": eid in fpkm.index,
                "in_counts": eid in counts.index,
                "gtf_symbol": ens2sym.get(eid, ""),
                "mean_fpkm_all": float(fpkm.loc[eid].mean()) if eid in fpkm.index else np.nan,
            }
        )
    pd.DataFrame(confirm).to_csv(TABLES / "gene_id_confirm.tsv", sep="\t", index=False)
    for gene, eid in NAMED_ENSG.items():
        if eid not in fpkm.index:
            raise RuntimeError(f"{gene} {eid} missing from FPKM")

    # Symbol matrix: named IDs win; remaining symbols max-mean collapsed.
    named_expr = pd.DataFrame({g: fpkm.loc[eid] for g, eid in NAMED_ENSG.items()}).T
    collapsed = collapse_symbol(fpkm, ens2sym)
    for g in named_expr.index:
        collapsed.loc[g] = named_expr.loc[g]
    gene_fpkm = collapsed
    log_fpkm = np.log2(gene_fpkm.clip(lower=0) + 1.0)

    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    named_cpm = pd.DataFrame({g: cpm.loc[eid] for g, eid in NAMED_ENSG.items() if eid in cpm.index}).T
    collapsed_cpm = collapse_symbol(cpm, ens2sym)
    for g in named_cpm.index:
        collapsed_cpm.loc[g] = named_cpm.loc[g]
    log_cpm = np.log2(collapsed_cpm.clip(lower=0) + 1.0)

    tumors = meta.index[meta["tissue"] == "tumor"].tolist()
    normals = meta.index[meta["tissue"] == "normal"].tolist()
    luad = meta.index[(meta["tissue"] == "tumor") & (meta["histology_code"] == "2")].tolist()
    lusc = meta.index[(meta["tissue"] == "tumor") & (meta["histology_code"] == "1")].tolist()
    large = meta.index[(meta["tissue"] == "tumor") & (meta["histology_code"] == "3")].tolist()

    if len(meta) != 218:
        raise RuntimeError(f"expected 218 GSM, got {len(meta)}")
    if len(tumors) != 199:
        raise RuntimeError(f"expected 199 tumors, got {len(tumors)}")
    if len(normals) != 19:
        raise RuntimeError(f"expected 19 normals, got {len(normals)}")
    if len(luad) + len(lusc) + len(large) != 199:
        raise RuntimeError("histology codes do not partition tumors")

    hist_counts = (
        meta.loc[tumors, "histology_label"]
        .value_counts()
        .rename_axis("histology_label")
        .reset_index(name="n")
    )
    hist_counts["code"] = hist_counts["histology_label"].map(
        {v: k for k, v in HISTOLOGY_MAP.items()}
    )
    hist_counts.to_csv(TABLES / "histology_counts.tsv", sep="\t", index=False)

    sample = meta.copy()
    sample["patient_id"] = sample["sample_id"].str.replace(r"[TN]$", "", regex=True)
    sample["cldn4_log2fpkm"] = log_fpkm.loc["CLDN4"].reindex(sample.index).values
    sample["cd8a_log2fpkm"] = log_fpkm.loc["CD8A"].reindex(sample.index).values
    sample["cd274_log2fpkm"] = log_fpkm.loc["CD274"].reindex(sample.index).values
    sample["tacstd2_log2fpkm"] = log_fpkm.loc["TACSTD2"].reindex(sample.index).values
    sample["is_luad"] = sample.index.isin(luad)
    sample["is_lusc"] = sample.index.isin(lusc)
    sample["is_large_cell"] = sample.index.isin(large)
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    # Immune-gene coverage after symbol collapse.
    cov = pd.DataFrame(
        {
            "gene": immune_genes,
            "present_after_collapse": [g in gene_fpkm.index for g in immune_genes],
        }
    )
    cov.to_csv(TABLES / "immunescore_gene_coverage.tsv", sep="\t", index=False)

    pair_rows, hl_rows, infos = [], [], []
    score_store = {}
    for name, gsms in [
        ("NSCLC_tumor", tumors),
        ("LUAD_AC_unspecified", luad),
        ("LUSC_squamous", lusc),
        ("LargeCell_NOS_extra", large),
    ]:
        pr, hl, info, scores = analyze_cohort(name, gsms, log_fpkm, immune_genes, "FPKM_log2p1")
        pair_rows.extend(pr)
        hl_rows.extend(hl)
        infos.append(info)
        score_store[name] = scores
        # featureCounts sensitivity on the same samples
        pr2, hl2, info2, _ = analyze_cohort(name, gsms, log_cpm, immune_genes, "featureCounts_log2CPM")
        pair_rows.extend(pr2)
        hl_rows.extend(hl2)
        infos.append(info2)

    pairs_df = pd.DataFrame(pair_rows)
    pairs_df.to_csv(TABLES / "spearman_cldn4_vs_partners.tsv", sep="\t", index=False)
    hl_df = pd.DataFrame(hl_rows)
    hl_df.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    n_tumor = len(tumors)
    n_luad = len(luad)
    n_lusc = len(lusc)
    n_large = len(large)
    labels = pd.DataFrame(
        [
            {"field": "GSM in series / FPKM columns", "public": "yes", "n": 218, "note": "199 tumor + 19 paired normal; do not use 218 as a histology n"},
            {"field": "unique patients (tumor GSM)", "public": "yes", "n": n_tumor, "note": "one tumor array per patient; this is the NSCLC tumor n"},
            {"field": "paired normal arrays", "public": "yes", "n": 19, "note": "inventoried; not mixed into tumor Spearman"},
            {"field": "histology=2 AC unspecified (LUAD)", "public": "yes", "n": n_luad, "note": "GEO overall_design code 2; this is the LUAD n"},
            {"field": "histology=1 squamous cell cancer (LUSC)", "public": "yes", "n": n_lusc, "note": "GEO overall_design code 1; this is the LUSC n"},
            {"field": "histology=3 large cell / NOS", "public": "yes", "n": n_large, "note": "not folded into LUAD or LUSC"},
            {"field": "stage tnm (tumors)", "public": "yes", "n": int(meta.loc[tumors, "stage_label"].ne("").sum()), "note": "1=1a … 7=IV; not a covariate"},
            {"field": "OS / vital date / dead", "public": "yes", "n": n_tumor, "note": "deposited; not this claim"},
            {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected Uppsala 2006–2010 atlas, not an ICI series"},
            {"field": "Tumor % / ABSOLUTE purity", "public": "no", "n": 0, "note": "only public proxy is an RNA epithelial score"},
            {"field": "CLDN4 ENSG00000189143 finite (tumors)", "public": "yes", "n": int(log_fpkm.loc["CLDN4", tumors].notna().sum()), "note": "Cufflinks FPKM"},
            {"field": "CD8A finite (tumors)", "public": "yes", "n": int(log_fpkm.loc["CD8A", tumors].notna().sum()), "note": "ENSG00000153563"},
            {"field": "CD274 finite (tumors)", "public": "yes", "n": int(log_fpkm.loc["CD274", tumors].notna().sum()), "note": "ENSG00000120217"},
            {"field": "TACSTD2 finite (tumors, companion)", "public": "yes", "n": int(log_fpkm.loc["TACSTD2", tumors].notna().sum()), "note": "companion column only; no dual-high"},
            {"field": "primary pairwise n (NSCLC tumors)", "public": "yes", "n": n_tumor, "note": "honest tumor n; not series 218"},
            {"field": "ESTIMATE ImmuneSignature genes present", "public": "yes", "n": int(cov["present_after_collapse"].sum()), "note": f"of {len(immune_genes)} Yoshihara 2013 genes"},
            {"field": "Ayers IFN-γ 6-gene present", "public": "yes", "n": sum(g in gene_fpkm.index for g in IFNG_AYERS), "note": ",".join(IFNG_AYERS)},
            {"field": "MHC-I genes present", "public": "yes", "n": sum(g in gene_fpkm.index for g in MHC1), "note": ",".join(MHC1)},
        ]
    )
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # One-row table: NSCLC tumors + LUAD + LUSC (FPKM primary).
    primary_pairs = [
        "CLDN4 vs CD8A",
        "CLDN4 vs CD274",
        "CLDN4 vs ImmuneScore",
        "CLDN4 vs IFNg_Ayers",
        "CLDN4 vs MHC-I",
    ]
    one_rows = []
    for cohort, n, hist in [
        ("NSCLC_tumor", n_tumor, "NSCLC tumor (all histologies)"),
        ("LUAD_AC_unspecified", n_luad, "LUAD (AC unspecified, code 2)"),
        ("LUSC_squamous", n_lusc, "LUSC (squamous, code 1)"),
    ]:
        sub = pairs_df[(pairs_df.cohort == cohort) & (pairs_df.matrix == "FPKM_log2p1")]
        rec = {
            "dataset": "GSE81089",
            "cohort": cohort,
            "histology": hist,
            "n": n,
            "unit": "patient/array (tumor GSM)",
            "matrix": "Cufflinks FPKM log2(FPKM+1)",
            "CLDN4": "ENSG00000189143",
        }
        for pair, key in [
            ("CLDN4 vs CD8A", "CD8A"),
            ("CLDN4 vs CD274", "CD274"),
            ("CLDN4 vs ImmuneScore", "ImmuneScore"),
            ("CLDN4 vs IFNg_Ayers", "IFNg_Ayers"),
            ("CLDN4 vs MHC-I", "MHC1"),
        ]:
            row = sub[sub.pair == pair].iloc[0]
            rec[f"{key}_n"] = int(row["n"])
            rec[f"{key}_rho"] = float(row["rho"])
            rec[f"{key}_p"] = float(row["p"])
            rec[f"{key}_rho_adj"] = float(row["rho_adj_epithelial"])
            rec[f"{key}_p_adj"] = float(row["p_adj_epithelial"])
            rec[f"{key}_verdict_adj"] = row["verdict_adj"]
        one_rows.append(rec)
    one_df = pd.DataFrame(one_rows)
    one_df.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # Figures — extra, NSCLC tumors + histology split.
    nsclc = score_store["NSCLC_tumor"]
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.4))
    panels = [
        ("CD8A", nsclc["CD8A"], "CD8A  log2(FPKM+1)"),
        ("CD274", nsclc["CD274"], "CD274  log2(FPKM+1)"),
        ("ImmuneScore", nsclc["ImmuneScore"], "ImmuneScore  mean-z"),
        ("IFNg_Ayers", nsclc["IFNg_Ayers"], "IFN-γ Ayers6  mean-z"),
        ("MHC-I", nsclc["MHC-I"], "MHC-I  mean-z"),
        ("epithelial", nsclc["epithelial"], "epithelial  mean-z"),
    ]
    nsclc_pairs = pairs_df[(pairs_df.cohort == "NSCLC_tumor") & (pairs_df.matrix == "FPKM_log2p1")]
    x = nsclc["CLDN4"]
    for ax, (lab, y, ylab) in zip(axes.ravel(), panels):
        pair_name = "CLDN4 vs epithelial mean-z" if lab == "epithelial" else f"CLDN4 vs {lab}"
        row = nsclc_pairs[nsclc_pairs.pair == pair_name].iloc[0]
        ax.scatter(x.values, np.asarray(y.reindex(x.index), float), s=14, alpha=0.7, c="#2c5f8a", edgecolors="none")
        ax.set_xlabel("CLDN4  log2(FPKM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(f"n={int(row.n)}  ρ={row.rho:.3f}  p={fmt_p(row.p)}")
    fig.suptitle("GSE81089 NSCLC tumors (n=199; not series 218)", y=1.01, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_cldn4_vs_partners.pdf", bbox_inches="tight")
    plt.close(fig)

    plot = nsclc_pairs[nsclc_pairs.pair.isin(primary_pairs)].copy()
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    y_pos = np.arange(len(plot))
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        plot["rho"],
        y_pos,
        xerr=[plot["rho"] - plot["ci_low"], plot["ci_high"] - plot["rho"]],
        fmt="o",
        color="#2c5f8a",
        ecolor="#2c5f8a",
        capsize=3,
        ms=6,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot["pair"].str.replace("CLDN4 vs ", ""))
    ax.set_xlabel("Spearman ρ vs CLDN4 (bootstrap 95% CI)")
    ax.set_title(f"GSE81089 NSCLC tumors  n={n_tumor}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_forest.pdf")
    plt.close(fig)

    # Extra: LUAD vs LUSC forest.
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.4), sharex=True)
    for ax, cohort, title, n in [
        (axes[0], "LUAD_AC_unspecified", "LUAD  AC unspecified", n_luad),
        (axes[1], "LUSC_squamous", "LUSC  squamous", n_lusc),
    ]:
        sub = pairs_df[(pairs_df.cohort == cohort) & (pairs_df.matrix == "FPKM_log2p1") & (pairs_df.pair.isin(primary_pairs))]
        y_pos = np.arange(len(sub))
        ax.axvline(0, color="0.5", lw=0.8)
        ax.errorbar(
            sub["rho"],
            y_pos,
            xerr=[sub["rho"] - sub["ci_low"], sub["ci_high"] - sub["rho"]],
            fmt="o",
            color="#2c5f8a",
            ecolor="#2c5f8a",
            capsize=3,
            ms=6,
        )
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub["pair"].str.replace("CLDN4 vs ", ""))
        ax.set_title(f"{title}  n={n}")
        ax.set_xlabel("Spearman ρ vs CLDN4")
    fig.suptitle("GSE81089 histology split (honest n; not series 218)", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_forest_luad_lusc.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig3_forest_luad_lusc.pdf", bbox_inches="tight")
    plt.close(fig)

    # Extra: Q4 vs Q1 boxplots on NSCLC tumors.
    fig, axes = plt.subplots(1, 5, figsize=(13.0, 3.4))
    q = pd.qcut(nsclc["CLDN4"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    for ax, lab in zip(axes, ["CD8A", "CD274", "ImmuneScore", "IFNg_Ayers", "MHC-I"]):
        y = nsclc[lab].reindex(nsclc["CLDN4"].index)
        data = [y[q == qq].values for qq in ["Q1", "Q4"]]
        ax.boxplot(data, tick_labels=["Q1", "Q4"], widths=0.55)
        row = hl_df[(hl_df.cohort == "NSCLC_tumor") & (hl_df.matrix == "FPKM_log2p1") & (hl_df.endpoint == lab)].iloc[0]
        ax.set_title(f"{lab}\nMWU p={fmt_p(row.p)}")
        ax.set_ylabel(lab)
    fig.suptitle(f"GSE81089 CLDN4 Q4 vs Q1  NSCLC tumors n={n_tumor}", y=1.03, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_q4q1_boxplots.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig4_q4q1_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "dataset": "GSE81089",
        "pmid": "29282718",
        "title": series.get("title", ""),
        "platform": "GPL16791 Illumina HiSeq 2500",
        "processing": "Cufflinks 2.1.1 FPKM (Ensembl v73 GRCh37) as deposited; featureCounts 1.4.0-p1 sensitivity.",
        "unit": "patient/array (tumor GSM)",
        "n_gsm_series": 218,
        "n_tumor": n_tumor,
        "n_normal": 19,
        "n_luad": n_luad,
        "n_lusc": n_lusc,
        "n_large_cell": n_large,
        "cldn4": "ENSG00000189143",
        "tacstd2_companion": "ENSG00000184292",
        "no_dual_high": True,
        "immunescore": {
            "definition": "mean of gene-wise z-scores of ESTIMATE ImmuneSignature, z-scored within the analysis set",
            "source": "estimate 1.0.13 SI_geneset.gmt ImmuneSignature (Yoshihara 2013)",
            "n_genes_listed": 141,
            "n_genes_present": int(cov["present_after_collapse"].sum()),
            "not": "official R estimate::estimateScore ssGSEA",
        },
        "ifng": {"definition": "Ayers 2017 6-gene IFN-γ mean-z", "genes": IFNG_AYERS},
        "mhc1": {"definition": "classical MHC-I / APM mean-z", "genes": MHC1},
        "epithelial": {"genes": EPITHELIAL},
        "missing_public_labels": ["ICI", "tumor_percent", "ABSOLUTE purity"],
        "cohort_info": infos,
        "one_row": one_df.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(clean_json(summary), indent=2) + "\n")

    show = pairs_df[
        (pairs_df.matrix == "FPKM_log2p1")
        & (pairs_df.pair.isin(primary_pairs))
        & (pairs_df.cohort.isin(["NSCLC_tumor", "LUAD_AC_unspecified", "LUSC_squamous"]))
    ]
    print(show[["cohort", "pair", "n", "rho", "p", "rho_adj_epithelial", "p_adj_epithelial", "verdict_crude", "verdict_adj"]].to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
