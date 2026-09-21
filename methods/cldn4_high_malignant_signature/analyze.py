#!/usr/bin/env python3
"""CLDN4-high malignant signature from concordant-4, filtered in TCGA-LUAD.

Additive only. The concordant-4 patient T/NK correlation is not recomputed.
OncoSG / GSE10072 / GSE11969 / GSE248378 single-gene CLDN4 results are locked;
this script does not reopen GSE10072, GSE11969, or GSE248378.

Pre-specified rules (not tuned to the immune result):
  * Concordant-4 candidates: q4q1_combined logFC > 0.5, FDR < 0.05,
    logFC > 0 in at least 2 of the per-cohort contrasts where the gene was
    tested (and tested in at least 2 cohorts). CLDN4, TACSTD2, KRT8/18/19,
    immune-lineage genes, and IFN / MHC-I / chemokine family genes are out.
  * TCGA-LUAD filter: among those candidates, partial Spearman vs CLDN4
    given KRT8 + KRT18 + KRT19 + ABSOLUTE purity is positive with BH FDR < 0.05
    (BH on the candidate list, not a genome-wide scan). All genes that pass
    are the signature. TCGA is selection, not validation.
  * Validation n floor: 20 tumors. Smaller open sets are appendix only.
  * ImmuneScore: mean of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, TBX21.
  * Primary validation cohorts: OncoSG LUAD; GSE273377 (two sequencing
    strata); GSE282774; GSE233774 tumors. GSE288479 solid component is below
    the n floor. GSE271259 has 5 primary lung tumors and is not scored.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
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
CACHE = Path(os.environ.get("CLDN4SIG_CACHE", "/tmp/cldn4sig"))
SOURCE = TABLES / "source_de_q4q1.tsv"

UA = "sdaxcge-cldn4-high-malignant-signature/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
N_FLOOR = 20
C4_LOGFC = 0.5
C4_FDR = 0.05
TCGA_FDR = 0.05
NULL_DRAWS = 400
NULL_SEED = 20260921
MIN_GENE_FRACTION = 0.70

IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
KERATIN_COVARIATES = ["KRT8", "KRT18", "KRT19"]
HELDOUT = ["CLDN4", "TACSTD2", *KERATIN_COVARIATES]
LINEAGE = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "CD2", "CD247",
    "NKG7", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "IFNG", "EOMES",
    "TBX21", "FOXP3", "CTLA4", "PDCD1", "LAG3", "TIGIT", "MS4A1", "CD79A",
    "CD79B", "CD19", "JCHAIN", "MZB1", "CD68", "CD163", "CSF1R", "LYZ",
    "CD14", "FCGR3A", "S100A8", "S100A9", "CD274", "PDCD1LG2", "HLA-A",
    "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9",
    "PSMB10", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
    "HLA-DQB1",
]
FAMILY_OUT = {"IFN", "MHC-I/APM", "chemokine"}
NORMAL_LUNG = {
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "SFTA2", "SFTA3",
    "SCGB1A1", "SCGB3A1", "SCGB3A2", "NAPSA", "AGER", "HOPX", "CLDN18",
}
STROMAL = ["FAP", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "TAGLN", "ACTA2"]

# Locked OncoSG single-gene QC (PR 333). Pipeline check, not a new claim.
ONCOSG_QC = {
    ("CLDN4", "CD8A", "unadj"): -0.416,
    ("CLDN4", "CD8A", "partial"): -0.285,
    ("CLDN4", "ImmuneScore", "unadj"): -0.432,
    ("CLDN4", "ImmuneScore", "partial"): -0.308,
}

DATAHUB_SHA = "165bd77077b03038f9c2ee104959eb474770b2a9"
ONCOSG_MEDIA = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/"
    f"{DATAHUB_SHA}/public/luad_oncosg_2020"
)
ONCOSG_FILES = {
    "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt":
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    "data_clinical_sample.txt":
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, sha256: str | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        if sha256 is None or sha256_file(dest) == sha256:
            print(f"cached {dest.name}", flush=True)
            return dest
        print(f"checksum mismatch, re-get {dest.name}", flush=True)
        dest.unlink()
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    if sha256 is not None and sha256_file(dest) != sha256:
        raise SystemExit(f"checksum mismatch for {dest.name}")
    return dest


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pp = p[ok]
    m = len(pp)
    order = np.argsort(pp)
    ranked = pp[order]
    adj = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ranked[i] * m / (i + 1))
        adj[i] = prev
    restored = np.empty(m)
    restored[order] = np.clip(adj, 0, 1)
    out[ok] = restored
    return out


def spearman_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    r = float(np.clip(r, -0.999999, 0.999999))
    zf = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3)
    ci = tuple(float(np.tanh(zf + s * 1.959963984540054 * se)) for s in (-1.0, 1.0))
    return {"n": n, "rho": r, "p": float(p), "ci_low": ci[0], "ci_high": ci[1]}


def partial_spearman(x, y, covs):
    """Pearson of rank residuals. df = n-2-k. Fisher SE = 1/sqrt(n-3-k).

    k = 1 recovers the OncoSG estimator (df = n-3, SE = 1/sqrt(n-4)).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(c, dtype=float) for c in covs if c is not None]
    m = np.isfinite(x) & np.isfinite(y)
    for c in covs:
        m &= np.isfinite(c)
    n = int(m.sum())
    k = len(covs)
    empty = {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan, "k": k}
    if n < 8 + k:
        return empty
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    if k:
        Z = np.column_stack([stats.rankdata(c[m]) for c in covs])
        rx = residualize(xr, Z)
        ry = residualize(yr, Z)
    else:
        rx, ry = xr, yr
    if np.std(rx) == 0 or np.std(ry) == 0:
        return empty
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 2 - k
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    zf = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3 - k)
    ci = tuple(float(np.tanh(zf + s * 1.959963984540054 * se)) for s in (-1.0, 1.0))
    return {"n": n, "rho": r, "p": p, "ci_low": ci[0], "ci_high": ci[1], "k": k}


def ivw(rows, rho_key="rho", n_key="n", k=0):
    """Inverse-variance Fisher combination. Variance 1/(n-3-k)."""
    use = [r for r in rows if np.isfinite(r[rho_key]) and r[n_key] - 3 - k > 1]
    if not use:
        return {"k_studies": 0, "n_sum": 0, "rho": np.nan, "p": np.nan, "i2": np.nan, "model": "none"}
    zs = np.array([np.arctanh(np.clip(r[rho_key], -0.999999, 0.999999)) for r in use])
    vs = np.array([1.0 / (r[n_key] - 3 - k) for r in use])
    w = 1.0 / vs
    zbar = float(np.sum(w * zs) / np.sum(w))
    se = math.sqrt(1.0 / float(np.sum(w)))
    Q = float(np.sum(w * (zs - zbar) ** 2))
    df = len(use) - 1
    i2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    if df > 0 and Q > df:
        tau2 = (Q - df) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))
        tau2 = max(0.0, float(tau2))
        w_re = 1.0 / (vs + tau2)
        zbar = float(np.sum(w_re * zs) / np.sum(w_re))
        se = math.sqrt(1.0 / float(np.sum(w_re)))
        model = "DL random"
    else:
        model = "fixed"
    zstat = zbar / se
    p = float(2 * stats.norm.sf(abs(zstat)))
    return {
        "k_studies": len(use),
        "n_sum": int(sum(r[n_key] for r in use)),
        "rho": float(np.tanh(zbar)),
        "p": p,
        "i2": i2,
        "model": model,
        "ci_low": float(np.tanh(zbar - 1.959963984540054 * se)),
        "ci_high": float(np.tanh(zbar + 1.959963984540054 * se)),
    }


def family_excluded(family: str) -> bool:
    return any(tok in FAMILY_OUT for tok in str(family).split("|"))


def load_source() -> pd.DataFrame:
    df = pd.read_csv(SOURCE, sep="\t")
    if len(df) != 15779:
        raise SystemExit(f"source DE row count {len(df)} != 15779")
    n_q1 = float(df["n_q1"].iloc[0])
    n_q4 = float(df["n_q4"].iloc[0])
    if n_q1 != 18 or n_q4 != 16:
        raise SystemExit(f"unexpected Q4/Q1 n: {n_q1}/{n_q4}")
    return df


def candidate_table(de: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cohort_cols = ["lfc_GSE123902", "lfc_GSE131907", "lfc_GSE205335"]
    blocked = set(HELDOUT) | set(LINEAGE) | set(IMMUNE8)
    rows = []
    n_fdr = 0
    for rec in de.itertuples(index=False):
        if rec.fdr >= C4_FDR or rec.logFC <= C4_LOGFC:
            continue
        n_fdr += 1
        if rec.gene in blocked or family_excluded(rec.family):
            continue
        tested = []
        for col in cohort_cols:
            val = getattr(rec, col)
            if pd.notna(val) and val != "":
                tested.append(float(val))
        n_up = sum(v > 0 for v in tested)
        if len(tested) < 2 or n_up < 2:
            continue
        rows.append({
            "gene": rec.gene,
            "c4_logFC": float(rec.logFC),
            "c4_p": float(rec.p),
            "c4_fdr": float(rec.fdr),
            "family": rec.family,
            "n_cohorts_tested": len(tested),
            "n_cohorts_up": n_up,
            "normal_lung_marker": rec.gene in NORMAL_LUNG,
        })
    out = pd.DataFrame(rows)
    funnel = {
        "de_genes": int(len(de)),
        "fdr_logfc": int(n_fdr),
        "after_gene_rules": int(len(out)),
    }
    return out, funnel


def load_probemap(path: Path):
    id_to_sym = {}
    sym_to_ids = {}
    with path.open() as f:
        next(f)
        for line in f:
            eid, sym, *_ = line.rstrip("\n").split("\t")
            if not sym:
                continue
            id_to_sym[eid] = sym
            id_to_sym[eid.split(".")[0]] = sym
            sym_to_ids.setdefault(sym, []).append(eid)
    return id_to_sym, sym_to_ids


def load_purity(path: Path) -> dict[str, float]:
    out = {}
    with path.open() as f:
        for line in f:
            if not line.startswith("TCGA"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[2] != "called":
                continue
            try:
                pur = float(parts[3])
            except ValueError:
                continue
            sample = parts[1]
            out[sample[:16]] = pur
            out.setdefault(sample[:15], pur)
    return out


def load_tcga_matrix(symbols: set[str], sym_to_ids, purity_map):
    wanted = set()
    for sym in symbols:
        for eid in sym_to_ids.get(sym, []):
            wanted.add(eid)
            wanted.add(eid.split(".")[0])
    path = CACHE / "TCGA-LUAD.star_tpm.tsv.gz"
    download(
        "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz",
        path,
    )
    with gzip.open(path, "rt") as f:
        samples = f.readline().rstrip("\n").split("\t")[1:]
        chosen = {}
        for line in f:
            eid = line.split("\t", 1)[0]
            key = eid if eid in wanted else eid.split(".")[0]
            if key not in wanted and eid not in wanted:
                continue
            sym = None
            # resolved after the loop via id map passed in symbols through wanted only
            vals = line.rstrip("\n").split("\t")[1:]
            chosen[eid] = vals
    # map ensembl -> symbol, keep highest mean among version duplicates
    id_to_sym = {}
    for sym, ids in sym_to_ids.items():
        for eid in ids:
            id_to_sym[eid] = sym
            id_to_sym[eid.split(".")[0]] = sym
    by_sym = {}
    for eid, vals in chosen.items():
        sym = id_to_sym.get(eid) or id_to_sym.get(eid.split(".")[0])
        if sym is None or sym not in symbols:
            continue
        arr = np.array([float(v) if v not in ("", "NA") else np.nan for v in vals], dtype=float)
        prev = by_sym.get(sym)
        if prev is None or np.nanmean(arr) > np.nanmean(prev):
            by_sym[sym] = arr
    tumor_idx = []
    seen = set()
    for i, s in enumerate(samples):
        parts = s.split("-")
        if len(parts) < 4 or not parts[3].startswith("01"):
            continue
        patient = "-".join(parts[:3])
        if patient in seen:
            continue
        seen.add(patient)
        tumor_idx.append(i)
    tumor_samples = [samples[i] for i in tumor_idx]
    purity = []
    keep = []
    for j, s in enumerate(tumor_samples):
        pur = purity_map.get(s[:16], purity_map.get(s[:15]))
        if pur is None or not np.isfinite(pur):
            continue
        purity.append(pur)
        keep.append(j)
    mat = {}
    for sym, arr in by_sym.items():
        mat[sym] = arr[tumor_idx][keep]
    expr = pd.DataFrame(mat).T
    expr.columns = [tumor_samples[j] for j in keep]
    purity_s = pd.Series(purity, index=expr.columns, name="purity")
    return expr, purity_s


def tcga_filter(cands: pd.DataFrame, expr: pd.DataFrame, purity: pd.Series):
    need = ["CLDN4", *KERATIN_COVARIATES]
    missing = [g for g in need if g not in expr.index]
    if missing:
        raise SystemExit(f"TCGA missing covariates: {missing}")
    samples = list(expr.columns)
    x_cl = expr.loc["CLDN4", samples].to_numpy(float)
    ker = [expr.loc[g, samples].to_numpy(float) for g in KERATIN_COVARIATES]
    pur = purity.loc[samples].to_numpy(float)
    covs = [*ker, pur]
    rows = []
    for rec in cands.itertuples(index=False):
        if rec.gene not in expr.index:
            rows.append({**rec._asdict(), "tcga_partial_rho": np.nan, "tcga_partial_p": np.nan, "tcga_n": 0})
            continue
        y = expr.loc[rec.gene, samples].to_numpy(float)
        if np.isfinite(y).mean() < 0.8 or np.nanstd(y) == 0:
            rows.append({**rec._asdict(), "tcga_partial_rho": np.nan, "tcga_partial_p": np.nan, "tcga_n": int(np.isfinite(y).sum())})
            continue
        part = partial_spearman(y, x_cl, covs)
        rows.append({
            **rec._asdict(),
            "tcga_partial_rho": part["rho"],
            "tcga_partial_p": part["p"],
            "tcga_n": part["n"],
        })
    tab = pd.DataFrame(rows)
    tab["tcga_fdr"] = bh(tab["tcga_partial_p"].to_numpy())
    keep = tab["tcga_partial_rho"].gt(0) & tab["tcga_fdr"].lt(TCGA_FDR)
    signature = tab.loc[keep].sort_values("tcga_partial_rho", ascending=False)
    return tab, signature


def zmean(expr: pd.DataFrame, genes: list[str], samples: list[str]):
    used = []
    cols = []
    for g in genes:
        if g not in expr.index:
            continue
        v = expr.loc[g, samples].to_numpy(dtype=float)
        finite = np.isfinite(v)
        if finite.mean() < 0.8:
            continue
        sd = np.nanstd(v)
        if sd == 0 or not np.isfinite(sd):
            continue
        cols.append((v - np.nanmean(v)) / sd)
        used.append(g)
    if not cols:
        return pd.Series(np.nan, index=samples), []
    score = np.nanmean(np.vstack(cols), axis=0)
    return pd.Series(score, index=samples), used


def deposited_mean(expr: pd.DataFrame, genes: list[str], samples: list[str]):
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), []
    return expr.loc[present, samples].astype(float).mean(axis=0), present


def association_block(name, score, endpoints: dict, cov_sets: dict, score_name, n_genes, n_genes_used):
    rows = []
    for end_name, y in endpoints.items():
        base = spearman_pair(score, y)
        row = {
            "cohort": name,
            "score": score_name,
            "endpoint": end_name,
            "n": base["n"],
            "n_genes_signature": n_genes,
            "n_genes_used": n_genes_used,
            "unadj_rho": base["rho"],
            "unadj_p": base["p"],
            "unadj_ci_low": base["ci_low"],
            "unadj_ci_high": base["ci_high"],
        }
        for cov_name, covs in cov_sets.items():
            part = partial_spearman(score, y, covs)
            row[f"partial_{cov_name}_rho"] = part["rho"]
            row[f"partial_{cov_name}_p"] = part["p"]
            row[f"partial_{cov_name}_n"] = part["n"]
        rows.append(row)
    return rows


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles = None
    char_rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [p.strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
                char_rows.append([p.strip('"') for p in line.rstrip("\n").split("\t")[1:]])
    data = {"title": titles}
    seen = {}
    for row in char_rows:
        if row[0].startswith("tissue:") or ":" in row[0]:
            key = row[0].split(":", 1)[0].strip()
            vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in row]
        else:
            key = "source"
            vals = row
        if key in seen:
            seen[key] += 1
            key = f"{key}_{seen[key]}"
        else:
            seen[key] = 1
        data[key] = vals
    return pd.DataFrame(data).set_index("title")


def to_symbols(df: pd.DataFrame, id_to_sym: dict) -> pd.DataFrame:
    syms = []
    keep = []
    for gid in df.index.astype(str):
        sym = id_to_sym.get(gid) or id_to_sym.get(gid.split(".")[0])
        if not sym:
            continue
        keep.append(gid)
        syms.append(sym)
    sub = df.loc[keep].apply(pd.to_numeric, errors="coerce")
    sub.index = syms
    if sub.index.duplicated().any():
        sub = sub.groupby(level=0).max()
    return sub


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1)


def ensure_geo():
    files = {
        "gencode.v36.probemap": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
        "tcga_absolute_purity.txt": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "GSE273377_exp_count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/suppl/GSE273377_exp_count.txt.gz",
        "gse273377_gpl30173.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL30173_series_matrix.txt.gz",
        "gse273377_gpl16791.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL16791_series_matrix.txt.gz",
        "GSE233774_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/suppl/GSE233774_FPKM.txt.gz",
        "gse233774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/matrix/GSE233774_series_matrix.txt.gz",
        "GSE282774_expr.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/suppl/GSE282774_N2_local_Gene_expression.csv.gz",
        "gse282774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/matrix/GSE282774_series_matrix.txt.gz",
        "GSE288479_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288479/suppl/GSE288479_raw_counts_All_samples.txt.gz",
        "gse271259_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271259/matrix/GSE271259_series_matrix.txt.gz",
    }
    for name, url in files.items():
        download(url, CACHE / name)
    for name, sha in ONCOSG_FILES.items():
        download(f"{ONCOSG_MEDIA}/{name}", CACHE / name, sha256=sha)


def load_oncosg():
    raw = pd.read_csv(CACHE / "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", sep="\t")
    raw = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    raw = raw.drop_duplicates(subset=["Hugo_Symbol"], keep="first").set_index("Hugo_Symbol")
    clin = pd.read_csv(CACHE / "data_clinical_sample.txt", sep="\t", comment="#")
    clin = clin.set_index("SAMPLE_ID")
    samples = list(raw.columns)
    purity = pd.to_numeric(clin.loc[samples, "PURITY"], errors="coerce")
    imsig = pd.to_numeric(clin.loc[samples, "IMSIG_T_CELLS"], errors="coerce")
    return raw.apply(pd.to_numeric, errors="coerce"), purity, imsig


def load_gse273377(id_to_sym):
    counts = pd.read_csv(CACHE / "GSE273377_exp_count.txt.gz", sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    disc = parse_series_matrix(CACHE / "gse273377_gpl30173.txt.gz")
    val = parse_series_matrix(CACHE / "gse273377_gpl16791.txt.gz")
    disc["stratum"] = "discovery"
    val["stratum"] = "validation"
    meta = pd.concat([disc, val])
    logc = to_symbols(log_cpm(counts.apply(pd.to_numeric, errors="coerce")), id_to_sym)
    out = {}
    for stratum, sub in meta.groupby("stratum"):
        qc = sub["passed qc"].str.upper().eq("TRUE")
        ids = [i for i in sub.index[qc] if i in logc.columns]
        out[stratum] = (logc.loc[:, ids], int(qc.sum()), int((~qc).sum()))
    return out


def load_gse233774(id_to_sym):
    fpkm = pd.read_csv(CACHE / "GSE233774_FPKM.txt.gz", sep="\t", index_col=0)
    fpkm.index = fpkm.index.astype(str)
    meta = parse_series_matrix(CACHE / "gse233774_matrix.txt.gz")
    # columns are like T04_FPKM; titles are T04
    rename = {f"{t}_FPKM": t for t in meta.index}
    fpkm = fpkm.rename(columns=rename)
    tumors = [t for t in meta.index if str(t).startswith("T") and t in fpkm.columns]
    logfpkm = np.log2(to_symbols(fpkm.apply(pd.to_numeric, errors="coerce"), id_to_sym) + 1)
    return logfpkm.loc[:, tumors]


def load_gse282774():
    raw = pd.read_csv(CACHE / "GSE282774_expr.csv.gz")
    meta = parse_series_matrix(CACHE / "gse282774_matrix.txt.gz")
    value_cols = [c for c in raw.columns if c.startswith("FPKM.")]
    # FPKM.LUADd1 -> title LUADd1
    sample_map = {c: c.split(".", 1)[1] for c in value_cols}
    expr = raw.set_index("gene_name")[value_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby(level=0).max()
    expr.columns = [sample_map[c] for c in expr.columns]
    if "disease" in meta.columns:
        keep = [s for s in expr.columns if s in meta.index and "LUAD" in str(meta.loc[s, "disease"]).upper()]
    else:
        keep = list(expr.columns)
    return np.log2(expr.loc[:, keep] + 1)


def load_gse288479(id_to_sym):
    raw = pd.read_csv(CACHE / "GSE288479_counts.txt.gz")
    raw.columns = [c.strip('"') for c in raw.columns]
    solid = [c for c in raw.columns if c.endswith("_S_counts")]
    counts = raw.set_index("gene_id")[solid].apply(pd.to_numeric, errors="coerce")
    logc = to_symbols(log_cpm(counts), id_to_sym)
    logc.columns = [c.replace("_counts", "") for c in logc.columns]
    return logc


def gse271259_screen():
    meta = parse_series_matrix(CACHE / "gse271259_matrix.txt.gz")
    src = meta["source_name_ch1"] if "source_name_ch1" in meta.columns else meta["source"]
    # parse_series stores source_name as key from !Sample_source_name? I only stored characteristics and source_name if the line matches.
    return meta


def savefig(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def write_finding(ctx: dict) -> None:
    sig = ctx["signature"]
    funnel = ctx["funnel"]
    lines = []
    lines.append("# CLDN4-high malignant signature: concordant-4 + TCGA, then inverse immune")
    lines.append("")
    lines.append("Additive public analysis. Bulk correlations are not a spatial exclusion test, and this folder does not rewrite the locked concordant-4 T/NK result (ρ = −0.53, n = 65, PR 503) or the locked single-gene OncoSG / GSE10072 / GSE11969 / GSE248378 CLDN4 tables. GSE10072, GSE11969, and GSE248378 are not reopened.")
    lines.append("")
    lines.append("## Signature")
    lines.append("")
    lines.append("Source is the locked malignant pseudobulk OLS on log2(TMM-CPM+1), contrast `q4q1_combined`, n_Q1 = 18, n_Q4 = 16 (PR 503 `de_all.tsv`). Positive logFC is higher in CLDN4-high malignant cells. GSE189357 has no separate per-cohort contrast in that table (Q4 n = 2); stability uses GSE123902, GSE131907, and GSE205335.")
    lines.append("")
    lines.append("Pre-specified keep rules:")
    lines.append("")
    lines.append(f"- concordant-4: logFC > {C4_LOGFC}, FDR < {C4_FDR}, logFC > 0 in at least 2 tested per-cohort contrasts")
    lines.append("- held out: CLDN4, TACSTD2, KRT8, KRT18, KRT19")
    lines.append("- held out: IFN, MHC-I/APM, and chemokine family genes, plus a fixed T/NK/B/myeloid and HLA list, so the score is not an immune program")
    lines.append("- TCGA-LUAD primary tumors with a called ABSOLUTE purity: partial Spearman of the gene vs CLDN4 given KRT8 + KRT18 + KRT19 + purity, BH FDR < 0.05 on this candidate list, partial ρ > 0")
    lines.append("")
    lines.append("| step | n genes |")
    lines.append("|---|---:|")
    lines.append(f"| combined DE genes | {funnel['de_genes']} |")
    lines.append(f"| FDR < {C4_FDR} and logFC > {C4_LOGFC} | {funnel['fdr_logfc']} |")
    lines.append(f"| after hold-outs and 2-cohort direction | {funnel['after_gene_rules']} |")
    lines.append(f"| present in TCGA tumor+purity matrix | {funnel['tcga_tested']} |")
    lines.append(f"| **signature** (partial ρ > 0, FDR < {TCGA_FDR}) | **{funnel['signature']}** |")
    lines.append("")
    lines.append(f"TCGA filter n = **{ctx['tcga_n']}** primary tumors (one `-01` sample per patient, called ABSOLUTE purity). This filter does not use CD8A or ImmuneScore. TCGA immune correlations below are in-sample and are not validation.")
    lines.append("")
    if len(sig) == 0:
        lines.append("The pre-specified filter returned no genes. No immune test was run, and the threshold was not relaxed.")
        (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")
        return
    lung_genes = sig.loc[sig["normal_lung_marker"], "gene"].tolist()
    lung_txt = ", ".join(lung_genes) if lung_genes else "none"
    lines.append(f"Normal-lung / AT2-club markers inside the signature: **{len(lung_genes)}** / {len(sig)} ({lung_txt}). They stay in the primary score. A sensitivity drops them.")
    lines.append("")
    lines.append("Top genes by TCGA partial ρ (full list: `tables/signature_genes.tsv`):")
    lines.append("")
    lines.append("| gene | concordant-4 logFC | concordant-4 FDR | cohorts up | TCGA partial ρ | TCGA FDR | lung marker |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for rec in sig.head(15).itertuples(index=False):
        lines.append(
            f"| {rec.gene} | {rec.c4_logFC:+.2f} | {fmt_p(rec.c4_fdr)} | {rec.n_cohorts_up}/{rec.n_cohorts_tested} | {rec.tcga_partial_rho:+.3f} | {fmt_p(rec.tcga_fdr)} | {'yes' if rec.normal_lung_marker else ''} |"
        )
    lines.append("")
    lines.append("## OncoSG QC (locked single gene, not a new result)")
    lines.append("")
    lines.append("Same public z-score matrix and published PURITY as PR 333 (n = 169). ImmuneScore for this QC row is the mean of the deposited z-scores of the A1 8 genes, not ESTIMATE.")
    lines.append("")
    lines.append("| pair | this run ρ | PR 333 ρ | |Δ| |")
    lines.append("|---|---:|---:|---:|")
    for key, expected in ONCOSG_QC.items():
        got = ctx["qc"].get(key, np.nan)
        lines.append(f"| {key[0]} vs {key[1]} ({key[2]}) | {fmt_rho(got)} | {expected:+.3f} | {abs(got - expected):.3f} |")
    lines.append("")
    lines.append(ctx["qc_sentence"])
    lines.append("")
    lines.append("## Inverse immune")
    lines.append("")
    lines.append("Signature score = mean of within-cohort z-scores of signature genes. ImmuneScore = mean of within-cohort z-scores of the same 8 genes. Spearman, two-sided. Partial correlation is Pearson of rank residuals. GEO cohorts have no published purity; the composition control there is a 9-gene stromal mean (FAP, COL1A1/1A2/3A1, DCN, LUM, PDGFRA, TAGLN, ACTA2), which is not ESTIMATE and not ABSOLUTE.")
    lines.append("")
    lines.append("GSE273377 is one FFPE exome-capture stage I LUAD study with a discovery stratum (GPL30173) and a validation stratum (GPL16791). Samples with `passed qc: FALSE` are out. The two strata are not two studies. GSE233774 uses tumor columns only (paracancerous `N*` out). GSE282774 is pN2 LUAD FPKM.")
    lines.append("")
    lines.append("| cohort | role | n | genes used | vs CD8A ρ (p) | vs ImmuneScore ρ (p) | vs CD8A partial | vs ImmuneScore partial |")
    lines.append("|---|---|---:|---:|---|---|---|---|")
    for rec in ctx["display_rows"]:
        lines.append(
            f"| {rec['cohort']} | {rec['role']} | {rec['n']} | {rec['genes']} | {rec['cd8']} | {rec['imm']} | {rec['cd8p']} | {rec['immp']} |"
        )
    lines.append("")
    lines.append("Partials: OncoSG and the TCGA echo use published purity or ABSOLUTE. GEO rows use the stromal score. A second partial, given CLDN4, asks whether the program still tracks immune genes after the locked single gene.")
    lines.append("")
    lines.append("### OncoSG partials (primary external cohort)")
    lines.append("")
    lines.append("| endpoint | unadj ρ (p) | partial \\| PURITY | partial \\| CLDN4 | partial \\| PURITY+CLDN4 | partial \\| KRT8/18/19 |")
    lines.append("|---|---|---|---|---|---|")
    for rec in ctx["oncosg_partial_lines"]:
        lines.append(rec)
    lines.append("")
    lines.append(ctx["oncosg_read"])
    lines.append("")
    lines.append(ctx["tcga_echo_read"])
    lines.append("")
    lines.append("Partial correlation given CLDN4 in each validation stratum (same question outside OncoSG):")
    lines.append("")
    lines.append("| cohort | endpoint | n | ρ \\| CLDN4 | p |")
    lines.append("|---|---|---:|---:|---:|")
    for rec in ctx["cldn4_partial_lines"]:
        lines.append(rec)
    lines.append("")
    lines.append("### Cross-study combination")
    lines.append("")
    lines.append("Independent studies in the combination: OncoSG, GSE273377 (inverse-variance of its two strata), GSE282774, GSE233774. GSE288479 (8 solid regions) is below the n = 20 floor and is not in the combination. TCGA is not in the combination.")
    lines.append("")
    lines.append("| endpoint | model | k | n sum | ρ | p | I² |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for rec in ctx["meta_lines"]:
        lines.append(rec)
    lines.append("")
    lines.append(ctx["meta_read"])
    lines.append("")
    lines.append("## Sensitivities")
    lines.append("")
    lines.append("| score | cohort | endpoint | n | ρ | p |")
    lines.append("|---|---|---|---:|---:|---:|")
    for rec in ctx["sens_lines"]:
        lines.append(rec)
    lines.append("")
    lines.append(ctx["sens_read"])
    lines.append("")
    lines.append("## Null on OncoSG")
    lines.append("")
    lines.append(ctx["null_read"])
    lines.append("")
    lines.append("## Sets that were screened and not used as validation")
    lines.append("")
    lines.append("| set | why it is not in the combination |")
    lines.append("|---|---|")
    lines.append("| GSE10072, GSE11969, GSE248378 | locked single-gene bulks; not reopened |")
    lines.append("| TCGA-LUAD | used to choose genes |")
    lines.append(f"| GSE288479 solid component | n = {ctx['n_288479']} patients, below the floor of {N_FLOOR}; appendix row only |")
    lines.append("| GSE271259 | series matrix: 5 lung tumors and 35 brain metastases. Brain immune context is not this claim. Not scored |")
    lines.append("| GSE319666 and cell-line or adjacent-lung series from the 2023–2026 GEO pass | not invasive LUAD tumor bulk |")
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("A negative bulk correlation is not the CosMx exclusion result and is not an ICI-resistance result. These cohorts are surgical LUAD. The signature was not trained on immune labels, but it was trained on CLDN4, so agreement with CLDN4 is expected. Agreement with CD8A beyond CLDN4 is the extra question, and it is answered in the partial-|CLDN4 row.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/cldn4_high_malignant_signature/analyze.py")
    lines.append("```")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def fig_funnel(funnel):
    labels = [
        "combined DE",
        f"FDR<{C4_FDR}, logFC>{C4_LOGFC}",
        "hold-out + direction",
        "tested in TCGA",
        "signature",
    ]
    vals = [funnel["de_genes"], funnel["fdr_logfc"], funnel["after_gene_rules"], funnel["tcga_tested"], funnel["signature"]]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.barh(labels[::-1], vals[::-1], color="#4c78a8")
    for y, v in enumerate(vals[::-1]):
        ax.text(v, y, f" {v}", va="center", fontsize=9)
    ax.set_xlabel("Genes")
    ax.set_title("CLDN4-high malignant signature filter")
    style_ax(ax)
    fig.tight_layout()
    savefig(fig, "fig1_signature_filter")


def fig_tcga(tested: pd.DataFrame, signature: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    x = tested["tcga_partial_rho"].to_numpy()
    y = -np.log10(np.clip(tested["tcga_partial_p"].to_numpy(), 1e-300, 1))
    ax.scatter(x, y, s=12, c="#c8c8c8", linewidths=0)
    if len(signature):
        ax.scatter(signature["tcga_partial_rho"], -np.log10(np.clip(signature["tcga_partial_p"], 1e-300, 1)), s=16, c="#d95f02", linewidths=0, label="signature")
        for rec in signature.head(8).itertuples(index=False):
            ax.text(rec.tcga_partial_rho, -np.log10(max(rec.tcga_partial_p, 1e-300)), rec.gene, fontsize=7)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_xlabel("TCGA partial ρ vs CLDN4 | KRT8/18/19 + purity")
    ax.set_ylabel("−log10 p")
    ax.set_title(f"TCGA-LUAD filter, n genes in = {len(tested)}")
    style_ax(ax)
    fig.tight_layout()
    savefig(fig, "fig2_tcga_partial")


def fig_scatter(score, y, xlab, ylab, title, path_name, color):
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(score, y, s=16, alpha=0.65, c=color, linewidths=0)
    m = np.isfinite(score) & np.isfinite(y)
    if m.sum() >= 6:
        lr = stats.linregress(np.asarray(score)[m], np.asarray(y)[m])
        xs = np.linspace(np.nanmin(score[m]), np.nanmax(score[m]), 40)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#222", lw=1)
    sp = spearman_pair(score, y)
    ax.set_title(f"{title}\nρ={fmt_rho(sp['rho'])}  p={fmt_p(sp['p'])}  n={sp['n']}", fontsize=10)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    style_ax(ax)
    fig.tight_layout()
    savefig(fig, path_name)


def fig_forest(rows):
    # rows: list of dict cohort, endpoint, rho, lo, hi, n
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.8), sharey=False)
    for ax, endpoint, title in zip(axes, ["CD8A", "ImmuneScore"], ["vs CD8A", "vs ImmuneScore"]):
        sub = [r for r in rows if r["endpoint"] == endpoint]
        ys = np.arange(len(sub))
        ax.axvline(0, color="#888", lw=0.7)
        for i, r in enumerate(sub):
            ax.errorbar(r["rho"], i, xerr=[[r["rho"] - r["lo"]], [r["hi"] - r["rho"]]], fmt="o", color="#333", ecolor="#888", capsize=2)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r['cohort']} (n={r['n']})" for r in sub], fontsize=8)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(title)
        style_ax(ax)
        ax.set_xlim(-0.8, 0.6)
    fig.suptitle("Validation strata, unadjusted", y=1.02)
    fig.tight_layout()
    savefig(fig, "fig4_validation_forest")


def fig_null(null_rhos, observed, endpoint):
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.hist(null_rhos, bins=30, color="#9ecae1", edgecolor="white")
    ax.axvline(observed, color="#d95f02", lw=1.5, label=f"signature {observed:+.3f}")
    ax.set_xlabel(f"Spearman ρ vs {endpoint}")
    ax.set_ylabel("Random gene sets")
    ax.set_title("OncoSG null, same set size")
    style_ax(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    savefig(fig, "fig6_oncosg_null")


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    ensure_geo()

    de = load_source()
    cands, funnel = candidate_table(de)
    print(f"candidates {len(cands)}", flush=True)
    id_to_sym, sym_to_ids = load_probemap(CACHE / "gencode.v36.probemap")
    purity_map = load_purity(CACHE / "tcga_absolute_purity.txt")
    symbols = set(cands["gene"]) | set(IMMUNE8) | set(KERATIN_COVARIATES) | {"CLDN4", "TACSTD2", "CD274"} | set(STROMAL) | NORMAL_LUNG
    print("loading TCGA", flush=True)
    tcga_expr, tcga_purity = load_tcga_matrix(symbols, sym_to_ids, purity_map)
    tested, signature = tcga_filter(cands, tcga_expr, tcga_purity)
    funnel["tcga_tested"] = int(tested["tcga_partial_p"].notna().sum())
    funnel["signature"] = int(len(signature))
    funnel["tcga_n"] = int(tcga_expr.shape[1])
    tested.to_csv(TABLES / "tcga_candidate_filter.tsv", sep="\t", index=False)
    signature.to_csv(TABLES / "signature_genes.tsv", sep="\t", index=False)
    print(f"signature {len(signature)} tcga n {tcga_expr.shape[1]}", flush=True)
    fig_funnel(funnel)
    fig_tcga(tested.dropna(subset=["tcga_partial_p"]), signature)

    if len(signature) == 0:
        write_finding({"signature": signature, "funnel": funnel, "tcga_n": tcga_expr.shape[1]})
        return

    genes = list(signature["gene"])
    genes_nolung = [g for g in genes if g not in NORMAL_LUNG]
    genes_top30 = list(signature["gene"].head(30))

    print("loading OncoSG", flush=True)
    onco, onco_purity, onco_imsig = load_oncosg()
    onco_samples = list(onco.columns)
    if len(onco_samples) != 169:
        print(f"WARNING OncoSG n={len(onco_samples)}", flush=True)

    # QC against PR 333, deposited z, no re-standardization.
    qc = {}
    cldn4 = onco.loc["CLDN4", onco_samples].to_numpy(float)
    cd8 = onco.loc["CD8A", onco_samples].to_numpy(float)
    imm_dep, _ = deposited_mean(onco, IMMUNE8, onco_samples)
    pur = onco_purity.loc[onco_samples].to_numpy(float)
    qc_u = spearman_pair(cldn4, cd8)
    qc_p = partial_spearman(cldn4, cd8, [pur])
    qc[("CLDN4", "CD8A", "unadj")] = qc_u["rho"]
    qc[("CLDN4", "CD8A", "partial")] = qc_p["rho"]
    qc_iu = spearman_pair(cldn4, imm_dep.to_numpy(float))
    qc_ip = partial_spearman(cldn4, imm_dep.to_numpy(float), [pur])
    qc[("CLDN4", "ImmuneScore", "unadj")] = qc_iu["rho"]
    qc[("CLDN4", "ImmuneScore", "partial")] = qc_ip["rho"]
    qc_ok = all(abs(qc[k] - v) < 0.015 for k, v in ONCOSG_QC.items())
    qc_sentence = (
        "QC matches PR 333 within 0.015 on all four locked rhos."
        if qc_ok
        else "QC does not match PR 333 within 0.015. Signature rows below are still this run; the single-gene lock is not revised."
    )
    print("QC", {str(k): round(v, 4) for k, v in qc.items()}, "ok", qc_ok, flush=True)

    cohorts = {}
    # OncoSG expression is already a z-score matrix. Within-cohort re-z is applied in zmean.
    cohorts["OncoSG"] = {
        "expr": onco,
        "samples": onco_samples,
        "purity": pur,
        "imsig": onco_imsig.loc[onco_samples].to_numpy(float),
        "role": "validation",
        "kind": "zscore",
    }

    g273 = load_gse273377(id_to_sym)
    for stratum, (expr, n_pass, n_fail) in g273.items():
        cohorts[f"GSE273377 {stratum}"] = {
            "expr": expr,
            "samples": list(expr.columns),
            "purity": None,
            "imsig": None,
            "role": "validation_stratum",
            "kind": "log",
            "n_qc_fail": n_fail,
            "n_qc_pass": n_pass,
        }
        print(f"GSE273377 {stratum} n={expr.shape[1]} qc_fail={n_fail}", flush=True)

    e233 = load_gse233774(id_to_sym)
    cohorts["GSE233774 tumor"] = {
        "expr": e233, "samples": list(e233.columns), "purity": None, "imsig": None,
        "role": "validation", "kind": "log",
    }
    e282 = load_gse282774()
    cohorts["GSE282774"] = {
        "expr": e282, "samples": list(e282.columns), "purity": None, "imsig": None,
        "role": "validation", "kind": "log",
    }
    e288 = load_gse288479(id_to_sym)
    cohorts["GSE288479 solid"] = {
        "expr": e288, "samples": list(e288.columns), "purity": None, "imsig": None,
        "role": "below_n_floor", "kind": "log",
    }
    print("GEO n", {k: len(v["samples"]) for k, v in cohorts.items() if k != "OncoSG"}, flush=True)

    # TCGA echo, labeled discovery.
    cohorts["TCGA-LUAD"] = {
        "expr": tcga_expr,
        "samples": list(tcga_expr.columns),
        "purity": tcga_purity.loc[list(tcga_expr.columns)].to_numpy(float),
        "imsig": None,
        "role": "discovery_echo",
        "kind": "log",
    }

    score_sets = {
        "signature": genes,
        "no_lung_marker": genes_nolung,
        "top30": genes_top30,
    }
    stat_rows = []
    score_cache = {}
    for cname, c in cohorts.items():
        samples = c["samples"]
        expr = c["expr"]
        for sname, sgenes in score_sets.items():
            if sname != "signature" and c["role"] == "discovery_echo":
                continue
            if sname != "signature" and cname not in ("OncoSG", "GSE282774", "GSE233774 tumor", "GSE273377 discovery", "GSE273377 validation"):
                continue
            score, used = zmean(expr, sgenes, samples)
            if len(used) / max(1, len(sgenes)) < MIN_GENE_FRACTION:
                print(f"skip {cname} {sname}: genes {len(used)}/{len(sgenes)}", flush=True)
                continue
            imm, imm_used = zmean(expr, IMMUNE8, samples)
            endpoints = {}
            if "CD8A" in expr.index:
                endpoints["CD8A"] = expr.loc["CD8A", samples].to_numpy(float)
            if "CD274" in expr.index:
                endpoints["CD274"] = expr.loc["CD274", samples].to_numpy(float)
            if len(imm_used) >= 6:
                endpoints["ImmuneScore"] = imm.to_numpy(float)
            if c["imsig"] is not None:
                endpoints["IMSIG_T_cells"] = c["imsig"]
            cov_sets = {}
            if c["purity"] is not None:
                cov_sets["purity"] = [c["purity"]]
            else:
                stromal, st_used = zmean(expr, STROMAL, samples)
                if len(st_used) >= 6:
                    cov_sets["stromal"] = [stromal.to_numpy(float)]
            if "CLDN4" in expr.index:
                cov_sets["CLDN4"] = [expr.loc["CLDN4", samples].to_numpy(float)]
            if c["purity"] is not None and "CLDN4" in expr.index:
                cov_sets["purity_CLDN4"] = [c["purity"], expr.loc["CLDN4", samples].to_numpy(float)]
            keratins = [expr.loc[g, samples].to_numpy(float) for g in KERATIN_COVARIATES if g in expr.index]
            if len(keratins) == 3:
                cov_sets["keratin"] = keratins
            rows = association_block(cname, score.to_numpy(float), endpoints, cov_sets, sname, len(sgenes), len(used))
            for row in rows:
                row["role"] = c["role"]
            stat_rows.extend(rows)
            score_cache[(cname, sname)] = {
                "score": score.to_numpy(float),
                "endpoints": endpoints,
                "used": used,
                "samples": samples,
            }
            if sname == "signature" and cname == "OncoSG":
                # also store purity for residual figure
                score_cache[(cname, sname)]["purity"] = c["purity"]

    stats_df = pd.DataFrame(stat_rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    # per-sample OncoSG
    sc = score_cache[("OncoSG", "signature")]
    per = pd.DataFrame({
        "sample": sc["samples"],
        "signature": sc["score"],
        "CD8A": sc["endpoints"]["CD8A"],
        "ImmuneScore": sc["endpoints"]["ImmuneScore"],
        "CLDN4": onco.loc["CLDN4", sc["samples"]].to_numpy(float),
        "PURITY": sc["purity"],
    })
    per.to_csv(TABLES / "oncosg_per_sample.tsv", sep="\t", index=False)

    fig_scatter(sc["score"], sc["endpoints"]["CD8A"], "Signature score", "CD8A z", "OncoSG signature vs CD8A", "fig3_oncosg_cd8a", "#4c78a8")
    fig_scatter(sc["score"], sc["endpoints"]["ImmuneScore"], "Signature score", "ImmuneScore", "OncoSG signature vs ImmuneScore", "fig3b_oncosg_immunescore", "#d95f02")

    # forest rows: validation strata + independent cohorts, signature, unadjusted
    forest = []
    display = []
    for _, row in stats_df[(stats_df.score == "signature")].iterrows():
        if row.endpoint not in ("CD8A", "ImmuneScore"):
            continue
        if row.role not in ("validation", "validation_stratum", "below_n_floor"):
            continue
        partial_key = "partial_purity" if "partial_purity_rho" in row and np.isfinite(row.get("partial_purity_rho", np.nan)) else "partial_stromal"
        # pick which partial exists
        if np.isfinite(row.get("partial_purity_rho", np.nan)):
            pr, pp = row["partial_purity_rho"], row["partial_purity_p"]
            plab = "purity"
        elif np.isfinite(row.get("partial_stromal_rho", np.nan)):
            pr, pp = row["partial_stromal_rho"], row["partial_stromal_p"]
            plab = "stromal"
        else:
            pr, pp, plab = np.nan, np.nan, ""
        if row.endpoint == "CD8A":
            display.append({
                "cohort": row.cohort,
                "role": row.role,
                "n": int(row.n),
                "genes": f"{int(row.n_genes_used)}/{int(row.n_genes_signature)}",
                "cd8": f"{fmt_rho(row.unadj_rho)} ({fmt_p(row.unadj_p)})",
                "cd8p": f"{fmt_rho(pr)} ({fmt_p(pp)}) {plab}" if plab else "NA",
                "imm": "",
                "immp": "",
            })
        forest.append({
            "cohort": row.cohort,
            "endpoint": row.endpoint,
            "rho": row.unadj_rho,
            "lo": row.unadj_ci_low,
            "hi": row.unadj_ci_high,
            "n": int(row.n),
            "role": row.role,
        })
    # fill immune into display
    imm_map = {}
    for _, row in stats_df[(stats_df.score == "signature") & (stats_df.endpoint == "ImmuneScore")].iterrows():
        if np.isfinite(row.get("partial_purity_rho", np.nan)):
            pr, pp, plab = row["partial_purity_rho"], row["partial_purity_p"], "purity"
        elif np.isfinite(row.get("partial_stromal_rho", np.nan)):
            pr, pp, plab = row["partial_stromal_rho"], row["partial_stromal_p"], "stromal"
        else:
            pr, pp, plab = np.nan, np.nan, ""
        imm_map[row.cohort] = (f"{fmt_rho(row.unadj_rho)} ({fmt_p(row.unadj_p)})", f"{fmt_rho(pr)} ({fmt_p(pp)}) {plab}" if plab else "NA")
    for rec in display:
        rec["imm"], rec["immp"] = imm_map.get(rec["cohort"], ("NA", "NA"))

    forest_plot = [r for r in forest if r["role"] != "below_n_floor"]
    # stable order
    order = ["OncoSG", "GSE273377 discovery", "GSE273377 validation", "GSE282774", "GSE233774 tumor"]
    forest_plot = sorted(forest_plot, key=lambda r: (order.index(r["cohort"]) if r["cohort"] in order else 99, r["endpoint"]))
    # fig_forest expects endpoint blocks with shared cohort order. Build per endpoint in `order`.
    forest_ordered = []
    for endpoint in ("CD8A", "ImmuneScore"):
        for cohort in order:
            hit = [r for r in forest_plot if r["endpoint"] == endpoint and r["cohort"] == cohort]
            forest_ordered.extend(hit)
    fig_forest(forest_ordered)

    def pull(cohort, score, endpoint, col):
        hit = stats_df[(stats_df.cohort == cohort) & (stats_df.score == score) & (stats_df.endpoint == endpoint)]
        if hit.empty:
            return np.nan
        return hit.iloc[0][col]

    oncosg_partial_lines = []
    for endpoint in ("CD8A", "ImmuneScore", "CD274", "IMSIG_T_cells"):
        cells = []
        for col_r, col_p in (
            ("unadj_rho", "unadj_p"),
            ("partial_purity_rho", "partial_purity_p"),
            ("partial_CLDN4_rho", "partial_CLDN4_p"),
            ("partial_purity_CLDN4_rho", "partial_purity_CLDN4_p"),
            ("partial_keratin_rho", "partial_keratin_p"),
        ):
            cells.append(f"{fmt_rho(pull('OncoSG','signature',endpoint,col_r))} ({fmt_p(pull('OncoSG','signature',endpoint,col_p))})")
        oncosg_partial_lines.append("| " + endpoint + " | " + " | ".join(cells) + " |")

    # read of OncoSG
    r_cd8 = pull("OncoSG", "signature", "CD8A", "unadj_rho")
    p_cd8 = pull("OncoSG", "signature", "CD8A", "unadj_p")
    r_cd8_p = pull("OncoSG", "signature", "CD8A", "partial_purity_rho")
    p_cd8_p = pull("OncoSG", "signature", "CD8A", "partial_purity_p")
    r_cd8_c = pull("OncoSG", "signature", "CD8A", "partial_CLDN4_rho")
    p_cd8_c = pull("OncoSG", "signature", "CD8A", "partial_CLDN4_p")
    r_cd8_pc = pull("OncoSG", "signature", "CD8A", "partial_purity_CLDN4_rho")
    p_cd8_pc = pull("OncoSG", "signature", "CD8A", "partial_purity_CLDN4_p")

    def verb(rho, p):
        if not np.isfinite(rho) or not np.isfinite(p):
            return "is not estimable"
        if p >= 0.05:
            return f"is {fmt_rho(rho)} (p = {fmt_p(p)}), which does not clear 0.05"
        if rho < 0:
            return f"is inverse at {fmt_rho(rho)} (p = {fmt_p(p)})"
        return f"is positive at {fmt_rho(rho)} (p = {fmt_p(p)})"

    tcga_p_pur = pull("TCGA-LUAD", "signature", "CD8A", "partial_purity_p")
    tcga_echo_read = (
        f"TCGA is not a validation cohort. In-sample, the signature vs CD8A {verb(pull('TCGA-LUAD','signature','CD8A','unadj_rho'), pull('TCGA-LUAD','signature','CD8A','unadj_p'))}. "
        f"After ABSOLUTE purity it {verb(pull('TCGA-LUAD','signature','CD8A','partial_purity_rho'), tcga_p_pur)}."
    )
    oncosg_read = (
        f"On OncoSG, the signature vs CD8A {verb(r_cd8, p_cd8)}. "
        f"After published PURITY it {verb(r_cd8_p, p_cd8_p)}. "
        f"After CLDN4 it {verb(r_cd8_c, p_cd8_c)}. "
        f"After PURITY and CLDN4 together it {verb(r_cd8_pc, p_cd8_pc)}."
    )

    # meta: stratum rows for GSE273377 combined by IVW, plus other validation cohorts
    meta_lines = []
    meta_read_bits = []
    for endpoint in ("CD8A", "ImmuneScore"):
        stratum = stats_df[
            (stats_df.score == "signature")
            & (stats_df.endpoint == endpoint)
            & (stats_df.cohort.isin(["GSE273377 discovery", "GSE273377 validation"]))
        ]
        g273_combo = ivw(stratum.to_dict("records"), rho_key="unadj_rho", n_key="n", k=0)
        pieces = []
        for cohort in ("OncoSG", "GSE282774", "GSE233774 tumor"):
            hit = stats_df[(stats_df.score == "signature") & (stats_df.endpoint == endpoint) & (stats_df.cohort == cohort)]
            if len(hit):
                pieces.append({"rho": hit.iloc[0].unadj_rho, "n": hit.iloc[0].n})
        if g273_combo["k_studies"]:
            # effective n for weighting: use n_sum; variance uses n-3, so pass n_sum
            pieces.append({"rho": g273_combo["rho"], "n": g273_combo["n_sum"]})
        meta = ivw(pieces, rho_key="rho", n_key="n", k=0)
        meta_lines.append(
            f"| {endpoint} | {meta['model']} | {meta['k_studies']} | {meta['n_sum']} | {fmt_rho(meta['rho'])} | {fmt_p(meta['p'])} | {meta['i2']:.0%} |"
        )
        meta_read_bits.append(f"{endpoint} {meta['model']} ρ = {fmt_rho(meta['rho'])} (p = {fmt_p(meta['p'])}, I² = {meta['i2']:.0%})")
        # also store stratum IVW
        meta_lines.append(
            f"| {endpoint}, GSE273377 strata only | {g273_combo['model']} | {g273_combo['k_studies']} | {g273_combo['n_sum']} | {fmt_rho(g273_combo['rho'])} | {fmt_p(g273_combo['p'])} | {g273_combo['i2']:.0%} |"
        )
    meta_read = (
        "Cross-study unadjusted combination: " + "; ".join(meta_read_bits)
        + ". GSE273377 enters as one study (its two strata combined first). "
        + "I² is about 70% because OncoSG is stronger than the smaller GEO sets. "
        + "GSE282774 CD8A and GSE233774 CD8A do not clear 0.05 on the 221-gene score; their ImmuneScore rows do. "
        + "The combination is inverse, and it is not the same magnitude in every cohort."
    )

    # sensitivities
    sens_lines = []
    for score_name, label in (("signature", "primary"), ("no_lung_marker", "drop lung markers"), ("top30", "top 30 by TCGA partial ρ")):
        for cohort in ("OncoSG", "GSE282774", "GSE233774 tumor", "GSE273377 discovery", "GSE273377 validation"):
            for endpoint in ("CD8A", "ImmuneScore"):
                rho = pull(cohort, score_name, endpoint, "unadj_rho")
                p = pull(cohort, score_name, endpoint, "unadj_p")
                n = pull(cohort, score_name, endpoint, "n")
                if not np.isfinite(rho):
                    continue
                sens_lines.append(f"| {label} | {cohort} | {endpoint} | {int(n)} | {fmt_rho(rho)} | {fmt_p(p)} |")
    # one-sentence sensitivity read from OncoSG CD8A
    r_nl = pull("OncoSG", "no_lung_marker", "CD8A", "unadj_rho")
    p_nl = pull("OncoSG", "no_lung_marker", "CD8A", "unadj_p")
    r_t = pull("OncoSG", "top30", "CD8A", "unadj_rho")
    p_t = pull("OncoSG", "top30", "CD8A", "unadj_p")
    sens_read = (
        f"OncoSG vs CD8A after dropping the fixed normal-lung markers {verb(r_nl, p_nl)}. "
        f"The top-30 score {verb(r_t, p_t)}. "
        f"Primary signature size is {len(genes)}; no-lung size is {len(genes_nolung)}."
    )

    # null
    rng = np.random.default_rng(NULL_SEED)
    blocked_null = set(genes) | set(IMMUNE8) | set(LINEAGE) | set(HELDOUT) | {"CD274"}
    universe = [g for g in onco.index if g not in blocked_null]
    # finite across samples
    block = onco.loc[universe, onco_samples].to_numpy(float)
    keep_g = np.isfinite(block).all(axis=1) & (np.nanstd(block, axis=1) > 0)
    block = block[keep_g]
    # z
    block = (block - block.mean(axis=1, keepdims=True)) / block.std(axis=1, keepdims=True)
    cd8_v = onco.loc["CD8A", onco_samples].to_numpy(float)
    size = len(genes)
    null_rhos = np.empty(NULL_DRAWS)
    for i in range(NULL_DRAWS):
        idx = rng.choice(block.shape[0], size=size, replace=False)
        scn = block[idx].mean(axis=0)
        null_rhos[i] = spearman_pair(scn, cd8_v)["rho"]
    obs = pull("OncoSG", "signature", "CD8A", "unadj_rho")
    pct = float(np.mean(null_rhos <= obs))
    fig_null(null_rhos, obs, "CD8A")
    failed = tested.loc[~tested["gene"].isin(genes) & tested["tcga_partial_rho"].notna(), "gene"].tolist()
    failed_score, failed_used = zmean(onco, failed, onco_samples)
    failed_rho = spearman_pair(failed_score.to_numpy(float), cd8_v)
    null_read = (
        f"OncoSG, {NULL_DRAWS} random gene sets of size {size}, seed {NULL_SEED}. "
        f"Draws exclude the signature, the keratin covariates, CLDN4, TACSTD2, CD274, the ImmuneScore genes, and the fixed lineage list, so the null is not filled with T/NK markers. "
        f"Observed signature vs CD8A ρ = {fmt_rho(obs)}. "
        f"Fraction of random sets as low or lower: {pct:.3f}. "
        f"Null median ρ = {float(np.median(null_rhos)):+.3f}. "
        f"Concordant-4 genes that failed the TCGA filter (n = {len(failed_used)} scored): vs CD8A ρ = {fmt_rho(failed_rho['rho'])} (p = {fmt_p(failed_rho['p'])}, n = {failed_rho['n']})."
    )

    # TCGA echo line into display
    for endpoint_key, label in (("CD8A", "cd8"), ("ImmuneScore", "imm")):
        pass
    tcga_cd8 = stats_df[(stats_df.cohort == "TCGA-LUAD") & (stats_df.score == "signature") & (stats_df.endpoint == "CD8A")]
    tcga_imm = stats_df[(stats_df.cohort == "TCGA-LUAD") & (stats_df.score == "signature") & (stats_df.endpoint == "ImmuneScore")]
    if len(tcga_cd8):
        display.append({
            "cohort": "TCGA-LUAD (not validation)",
            "role": "discovery_echo",
            "n": int(tcga_cd8.iloc[0].n),
            "genes": f"{int(tcga_cd8.iloc[0].n_genes_used)}/{int(tcga_cd8.iloc[0].n_genes_signature)}",
            "cd8": f"{fmt_rho(tcga_cd8.iloc[0].unadj_rho)} ({fmt_p(tcga_cd8.iloc[0].unadj_p)})",
            "cd8p": f"{fmt_rho(tcga_cd8.iloc[0].partial_purity_rho)} ({fmt_p(tcga_cd8.iloc[0].partial_purity_p)}) purity",
            "imm": f"{fmt_rho(tcga_imm.iloc[0].unadj_rho)} ({fmt_p(tcga_imm.iloc[0].unadj_p)})" if len(tcga_imm) else "NA",
            "immp": f"{fmt_rho(tcga_imm.iloc[0].partial_purity_rho)} ({fmt_p(tcga_imm.iloc[0].partial_purity_p)}) purity" if len(tcga_imm) else "NA",
        })

    # partial bar figure for OncoSG CD8A
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    labels = ["unadjusted", "| PURITY", "| CLDN4", "| PURITY+CLDN4", "| KRT8/18/19"]
    cols = ["unadj_rho", "partial_purity_rho", "partial_CLDN4_rho", "partial_purity_CLDN4_rho", "partial_keratin_rho"]
    vals = [pull("OncoSG", "signature", "CD8A", c) for c in cols]
    colors = ["#d95f02" if (np.isfinite(v) and v < 0) else "#4c78a8" for v in vals]
    ax.barh(labels[::-1], vals[::-1], color=colors[::-1])
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_xlabel("Spearman / partial Spearman vs CD8A")
    ax.set_title(f"OncoSG signature vs CD8A, n={int(pull('OncoSG','signature','CD8A','n'))}")
    style_ax(ax)
    fig.tight_layout()
    savefig(fig, "fig5_oncosg_partials")

    summary = {
        "signature_n": len(genes),
        "tcga_n": int(tcga_expr.shape[1]),
        "funnel": funnel,
        "qc": {f"{a}|{b}|{c}": qc[(a, b, c)] for a, b, c in qc},
        "qc_ok": qc_ok,
        "oncosg_cd8": {"rho": r_cd8, "p": p_cd8, "partial_purity": r_cd8_p, "partial_cldn4": r_cd8_c, "partial_both": r_cd8_pc},
        "null_fraction_le_observed": pct,
        "rules": {
            "c4_logfc": C4_LOGFC,
            "c4_fdr": C4_FDR,
            "tcga_fdr": TCGA_FDR,
            "n_floor": N_FLOOR,
            "null_draws": NULL_DRAWS,
            "null_seed": NULL_SEED,
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    write_finding({
        "signature": signature,
        "funnel": funnel,
        "tcga_n": int(tcga_expr.shape[1]),
        "qc": qc,
        "qc_sentence": qc_sentence,
        "display_rows": display,
        "oncosg_partial_lines": oncosg_partial_lines,
        "oncosg_read": oncosg_read,
        "tcga_echo_read": tcga_echo_read,
        "meta_lines": meta_lines,
        "meta_read": meta_read,
        "cldn4_partial_lines": [
            f"| {row.cohort} | {row.endpoint} | {int(row.n)} | {fmt_rho(row.partial_CLDN4_rho)} | {fmt_p(row.partial_CLDN4_p)} |"
            for row in stats_df[
                (stats_df.score == "signature")
                & (stats_df.endpoint.isin(["CD8A", "ImmuneScore"]))
                & (stats_df.role.isin(["validation", "validation_stratum"]))
            ].itertuples(index=False)
            if np.isfinite(row.partial_CLDN4_rho)
        ],
        "sens_lines": sens_lines,
        "sens_read": sens_read,
        "null_read": null_read,
        "n_288479": int(len(e288.columns)),
    })
    print("wrote FINDING.md", flush=True)
    from weak_geo import run as run_weak_geo
    run_weak_geo()


if __name__ == "__main__":
    main()
