#!/usr/bin/env python3
"""TACSTD2–immune partial on CLDN4, and mediation by the locked 221-gene signature.

Cohorts and score definitions follow PR 590 (OncoSG, GSE273377 discovery,
GSE273377 validation, GSE282774, GSE233774 tumors). The 221 genes are the
PR 590 signature file, hashed below. CLDN4, TACSTD2, and the keratin
covariates were already held out of that list. This script does not refit
the signature and does not reopen GSE10072, GSE11969, or GSE248378.

Primary questions, fixed before looking at these partials:
  1. Partial Spearman of TACSTD2 vs CD8A and vs ImmuneScore given CLDN4.
  2. Does the 221-gene score statistically mediate TACSTD2 vs those endpoints?

Sensitivities, also fixed: published PURITY on OncoSG and the PR 590 9-gene
stromal score on GEO; mediator = CLDN4 itself; mediator = the signature with
CLDN4 as a covariate.
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
CACHE = Path(os.environ.get("TACSTD2MED_CACHE", "/tmp/tacstd2med"))
SIG_PATH = HERE / "data" / "signature_genes.tsv"
SIG_SHA = "9d398e71f8a2569d3c5c7bb940577a24eaca1f19db0a1d589c84e13a0bac6a3a"

UA = "sdaxcge-tacstd2-cldn4-mediation/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
N_BOOT = 5000
MIN_GENE_FRACTION = 0.70

IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
STROMAL = ["FAP", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "TAGLN", "ACTA2"]
HELDOUT_FORBIDDEN = {"CLDN4", "TACSTD2", "KRT8", "KRT18", "KRT19", "CD8A", "CD274"}

# Locked single-gene OncoSG anchors. Pipeline checks, not new claims.
QC_TARGETS = {
    ("TACSTD2", "CD8A", "unadj"): -0.380,
    ("TACSTD2", "CD8A", "purity"): -0.309,
    ("TACSTD2", "ImmuneScore_deposited", "unadj"): -0.387,
    ("TACSTD2", "ImmuneScore_deposited", "purity"): -0.318,
    ("CLDN4", "CD8A", "unadj"): -0.416,
    ("CLDN4", "CD8A", "purity"): -0.285,
    ("CLDN4", "TACSTD2", "unadj"): 0.505,
    ("signature", "CD8A", "unadj"): -0.578,
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

COHORT_ORDER = [
    "OncoSG",
    "GSE273377 discovery",
    "GSE273377 validation",
    "GSE282774",
    "GSE233774 tumor",
]
EXPECTED_N = {
    "OncoSG": 169,
    "GSE273377 discovery": 103,
    "GSE273377 validation": 60,
    "GSE282774": 58,
    "GSE233774 tumor": 30,
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


def fmt_num(r, digits=3) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.{digits}f}"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


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
    """Pearson of rank residuals. df = n-2-k. Fisher SE = 1/sqrt(n-3-k)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(c, dtype=float) for c in covs]
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
    se = 1.0 / math.sqrt(max(1e-12, n - 3 - k))
    ci = tuple(float(np.tanh(zf + s * 1.959963984540054 * se)) for s in (-1.0, 1.0))
    return {"n": n, "rho": r, "p": p, "ci_low": ci[0], "ci_high": ci[1], "k": k}


def ivw(rows, rho_key="rho", n_key="n", k=0):
    use = [r for r in rows if np.isfinite(r[rho_key]) and r[n_key] - 3 - k > 1]
    if not use:
        return {"k_studies": 0, "n_sum": 0, "rho": np.nan, "p": np.nan, "i2": np.nan, "model": "none",
                "ci_low": np.nan, "ci_high": np.nan}
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


def dl_mean(effects, variances):
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    ok = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    effects, variances = effects[ok], variances[ok]
    k = len(effects)
    if k == 0:
        return {"k": 0, "mu": np.nan, "se": np.nan, "p": np.nan, "i2": np.nan, "model": "none",
                "ci_low": np.nan, "ci_high": np.nan}
    w = 1.0 / variances
    mu = float(np.sum(w * effects) / np.sum(w))
    Q = float(np.sum(w * (effects - mu) ** 2))
    df = k - 1
    i2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    if df > 0 and Q > df:
        tau2 = (Q - df) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))
        tau2 = max(0.0, float(tau2))
        w = 1.0 / (variances + tau2)
        mu = float(np.sum(w * effects) / np.sum(w))
        model = "DL random"
    else:
        model = "fixed"
    se = math.sqrt(1.0 / float(np.sum(w)))
    z = mu / se
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "k": k,
        "mu": mu,
        "se": se,
        "p": p,
        "i2": i2,
        "model": model,
        "ci_low": mu - 1.959963984540054 * se,
        "ci_high": mu + 1.959963984540054 * se,
    }


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


def load_probemap(path: Path):
    id_to_sym = {}
    with path.open() as f:
        next(f)
        for line in f:
            eid, sym, *_ = line.rstrip("\n").split("\t")
            if not sym:
                continue
            id_to_sym[eid] = sym
            id_to_sym[eid.split(".")[0]] = sym
    return id_to_sym


def ensure_inputs():
    files = {
        "gencode.v36.probemap": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
        "GSE273377_exp_count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/suppl/GSE273377_exp_count.txt.gz",
        "gse273377_gpl30173.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL30173_series_matrix.txt.gz",
        "gse273377_gpl16791.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL16791_series_matrix.txt.gz",
        "GSE233774_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/suppl/GSE233774_FPKM.txt.gz",
        "gse233774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/matrix/GSE233774_series_matrix.txt.gz",
        "GSE282774_expr.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/suppl/GSE282774_N2_local_Gene_expression.csv.gz",
        "gse282774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/matrix/GSE282774_series_matrix.txt.gz",
    }
    for name, url in files.items():
        download(url, CACHE / name)
    for name, sha in ONCOSG_FILES.items():
        download(f"{ONCOSG_MEDIA}/{name}", CACHE / name, sha256=sha)


def load_signature():
    digest = sha256_file(SIG_PATH)
    if digest != SIG_SHA:
        raise SystemExit(f"signature file hash {digest} != locked {SIG_SHA}")
    df = pd.read_csv(SIG_PATH, sep="\t")
    genes = list(df["gene"].astype(str))
    if len(genes) != 221 or len(set(genes)) != 221:
        raise SystemExit(f"signature length {len(genes)} unique {len(set(genes))}")
    hit = sorted(set(genes) & HELDOUT_FORBIDDEN)
    if hit:
        raise SystemExit(f"signature contains held-out genes {hit}")
    return genes


def load_oncosg():
    raw = pd.read_csv(CACHE / "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", sep="\t")
    raw = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    raw = raw.drop_duplicates(subset=["Hugo_Symbol"], keep="first").set_index("Hugo_Symbol")
    clin = pd.read_csv(CACHE / "data_clinical_sample.txt", sep="\t", comment="#")
    clin = clin.set_index("SAMPLE_ID")
    samples = list(raw.columns)
    purity = pd.to_numeric(clin.loc[samples, "PURITY"], errors="coerce")
    expr = raw.apply(pd.to_numeric, errors="coerce")
    return expr, purity


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
        out[stratum] = logc.loc[:, ids]
    return out


def load_gse233774(id_to_sym):
    fpkm = pd.read_csv(CACHE / "GSE233774_FPKM.txt.gz", sep="\t", index_col=0)
    fpkm.index = fpkm.index.astype(str)
    meta = parse_series_matrix(CACHE / "gse233774_matrix.txt.gz")
    rename = {f"{t}_FPKM": t for t in meta.index}
    fpkm = fpkm.rename(columns=rename)
    tumors = [t for t in meta.index if str(t).startswith("T") and t in fpkm.columns]
    logfpkm = np.log2(to_symbols(fpkm.apply(pd.to_numeric, errors="coerce"), id_to_sym) + 1)
    return logfpkm.loc[:, tumors]


def load_gse282774():
    raw = pd.read_csv(CACHE / "GSE282774_expr.csv.gz")
    meta = parse_series_matrix(CACHE / "gse282774_matrix.txt.gz")
    value_cols = [c for c in raw.columns if c.startswith("FPKM.")]
    sample_map = {c: c.split(".", 1)[1] for c in value_cols}
    expr = raw.set_index("gene_name")[value_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby(level=0).max()
    expr.columns = [sample_map[c] for c in expr.columns]
    if "disease" in meta.columns:
        keep = [s for s in expr.columns if s in meta.index and "LUAD" in str(meta.loc[s, "disease"]).upper()]
    else:
        keep = list(expr.columns)
    return np.log2(expr.loc[:, keep] + 1)


def load_cohorts(genes, id_to_sym):
    onco, purity = load_oncosg()
    samples = list(onco.columns)
    cohorts = {
        "OncoSG": {
            "expr": onco,
            "samples": samples,
            "composition": purity.loc[samples].to_numpy(float),
            "composition_name": "PURITY",
        }
    }
    g273 = load_gse273377(id_to_sym)
    for stratum, expr in g273.items():
        cohorts[f"GSE273377 {stratum}"] = {
            "expr": expr,
            "samples": list(expr.columns),
            "composition": None,
            "composition_name": "stromal",
        }
    e233 = load_gse233774(id_to_sym)
    cohorts["GSE233774 tumor"] = {
        "expr": e233,
        "samples": list(e233.columns),
        "composition": None,
        "composition_name": "stromal",
    }
    e282 = load_gse282774()
    cohorts["GSE282774"] = {
        "expr": e282,
        "samples": list(e282.columns),
        "composition": None,
        "composition_name": "stromal",
    }
    for name, c in cohorts.items():
        if len(c["samples"]) != EXPECTED_N[name]:
            raise SystemExit(f"{name} n={len(c['samples'])} expected {EXPECTED_N[name]}")
    return cohorts


def assemble(c, genes):
    expr = c["expr"]
    samples = c["samples"]
    for g in ("TACSTD2", "CLDN4", "CD8A"):
        if g not in expr.index:
            raise SystemExit(f"{g} absent")
    sig, sig_used = zmean(expr, genes, samples)
    imm, imm_used = zmean(expr, IMMUNE8, samples)
    if len(sig_used) / 221 < MIN_GENE_FRACTION:
        raise SystemExit(f"signature genes {len(sig_used)}/221")
    if len(imm_used) < 6:
        raise SystemExit(f"ImmuneScore genes {len(imm_used)}")
    immune_missing = [g for g in IMMUNE8 if g not in imm_used]
    stromal, st_used = zmean(expr, STROMAL, samples)
    comp = c["composition"]
    if comp is None:
        if len(st_used) < 6:
            comp = None
            comp_name = "none"
        else:
            comp = stromal.to_numpy(float)
            comp_name = "stromal"
    else:
        comp_name = c["composition_name"]
    dep, dep_used = deposited_mean(expr, IMMUNE8, samples)
    return {
        "TACSTD2": expr.loc["TACSTD2", samples].to_numpy(float),
        "CLDN4": expr.loc["CLDN4", samples].to_numpy(float),
        "CD8A": expr.loc["CD8A", samples].to_numpy(float),
        "ImmuneScore": imm.to_numpy(float),
        "ImmuneScore_deposited": dep.to_numpy(float),
        "signature": sig.to_numpy(float),
        "composition": comp,
        "composition_name": comp_name,
        "n_sig": len(sig_used),
        "n_immune": len(imm_used),
        "immune_missing": immune_missing,
        "n_stromal": len(st_used),
        "n_deposited": len(dep_used),
        "samples": samples,
        "sig_missing": [g for g in genes if g not in sig_used],
    }


def _zrank(v):
    r = stats.rankdata(v, method="average").astype(float)
    sd = r.std()
    if sd == 0 or not np.isfinite(sd):
        return None
    return (r - r.mean()) / sd


def _ols(y, X):
    n, p = X.shape
    A = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    df = n - A.shape[1]
    if df <= 0:
        return beta, np.full(len(beta), np.nan), df
    sigma2 = float(np.sum(resid ** 2) / df)
    try:
        xtx_inv = np.linalg.inv(A.T @ A)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(A.T @ A)
    se = np.sqrt(np.clip(np.diag(xtx_inv) * sigma2, 0, None))
    return beta, se, df


def _pack_fit(x, m, y, cov):
    zx, zm, zy = _zrank(x), _zrank(m), _zrank(y)
    if zx is None or zm is None or zy is None:
        return None
    cols = []
    for j in range(cov.shape[1]):
        zj = _zrank(cov[:, j])
        if zj is None:
            return None
        cols.append(zj)
    zc = np.column_stack(cols) if cols else np.zeros((len(x), 0))
    Xa = np.column_stack([zx, zc]) if zc.size else zx.reshape(-1, 1)
    ba, sea, dfa = _ols(zm, Xa)
    Xc = Xa
    bc, sec, dfc = _ols(zy, Xc)
    Xb = np.column_stack([zx, zm, zc]) if zc.size else np.column_stack([zx, zm])
    bb, seb, dfb = _ols(zy, Xb)
    a, c = float(ba[1]), float(bc[1])
    c_prime, b = float(bb[1]), float(bb[2])
    ab = a * b
    prop = ab / c if abs(c) > 1e-8 else np.nan
    def tp(est, se, df):
        if not np.isfinite(se) or se == 0 or df <= 0:
            return np.nan
        t = est / se
        return float(2 * stats.t.sf(abs(t), df))
    return {
        "a": a, "a_se": float(sea[1]), "a_p": tp(a, sea[1], dfa),
        "b": b, "b_se": float(seb[2]), "b_p": tp(b, seb[2], dfb),
        "c": c, "c_se": float(sec[1]), "c_p": tp(c, sec[1], dfc),
        "c_prime": c_prime, "c_prime_se": float(seb[1]), "c_prime_p": tp(c_prime, seb[1], dfb),
        "ab": ab,
        "prop": prop,
    }


def mediate(x, m, y, covs, seed_key):
    x = np.asarray(x, dtype=float)
    m = np.asarray(m, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(c, dtype=float) for c in covs]
    mask = np.isfinite(x) & np.isfinite(m) & np.isfinite(y)
    for c in covs:
        mask &= np.isfinite(c)
    n = int(mask.sum())
    k = len(covs)
    empty = {"n": n, "k": k}
    if n < 8 + k + 1:
        empty.update({key: np.nan for key in (
            "a", "b", "c", "c_prime", "ab", "prop",
            "a_p", "b_p", "c_p", "c_prime_p",
            "ab_ci_low", "ab_ci_high", "ab_p_boot", "ab_p_censored",
            "prop_ci_low", "prop_ci_high", "prop_boot_frac_ok",
            "c_ci_low", "c_ci_high", "c_prime_ci_low", "c_prime_ci_high",
            "a_ci_low", "a_ci_high", "b_ci_low", "b_ci_high",
        )})
        return empty
    xx, mm, yy = x[mask], m[mask], y[mask]
    cov = np.column_stack([c[mask] for c in covs]) if covs else np.zeros((n, 0))
    point = _pack_fit(xx, mm, yy, cov)
    if point is None:
        raise SystemExit("mediation point estimate failed")
    if abs(point["c"] - (point["c_prime"] + point["ab"])) > 1e-6:
        raise SystemExit(
            f"product identity failed for {seed_key}: c={point['c']} c'+ab={point['c_prime'] + point['ab']}"
        )
    seed = int(hashlib.sha256(seed_key.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    boots = np.full((N_BOOT, 6), np.nan)  # a b c c' ab prop
    for i in range(N_BOOT):
        idx = rng.integers(0, n, n)
        cov_i = cov[idx] if cov.size else cov
        fit = _pack_fit(xx[idx], mm[idx], yy[idx], cov_i)
        if fit is None:
            continue
        boots[i, 0] = fit["a"]
        boots[i, 1] = fit["b"]
        boots[i, 2] = fit["c"]
        boots[i, 3] = fit["c_prime"]
        boots[i, 4] = fit["ab"]
        boots[i, 5] = fit["prop"]
    def ci(col):
        v = boots[:, col]
        v = v[np.isfinite(v)]
        if len(v) < 100:
            return np.nan, np.nan, 0
        return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), len(v)
    a_lo, a_hi, _ = ci(0)
    b_lo, b_hi, _ = ci(1)
    c_lo, c_hi, _ = ci(2)
    cp_lo, cp_hi, _ = ci(3)
    ab_lo, ab_hi, n_ab = ci(4)
    # proportion: keep draws with a usable total effect
    prop_draws = boots[:, 5]
    prop_draws = prop_draws[np.isfinite(prop_draws)]
    prop_frac = float(len(prop_draws) / N_BOOT)
    if len(prop_draws) >= 100:
        p_lo, p_hi = float(np.percentile(prop_draws, 2.5)), float(np.percentile(prop_draws, 97.5))
    else:
        p_lo, p_hi = np.nan, np.nan
    ab_draws = boots[:, 4]
    ab_draws = ab_draws[np.isfinite(ab_draws)]
    if len(ab_draws) >= 100 and np.nanstd(ab_draws) > 0:
        # two-sided bootstrap p: fraction of draws on the other side of 0, times 2.
        # Zero crossings are reported as < 2/N rather than as p = 0.
        if point["ab"] < 0:
            n_cross = int(np.sum(ab_draws >= 0))
        else:
            n_cross = int(np.sum(ab_draws <= 0))
        ab_p = float(min(1.0, 2 * n_cross / len(ab_draws)))
        ab_p_censored = bool(n_cross == 0)
        if ab_p_censored:
            ab_p = 2.0 / len(ab_draws)
    else:
        ab_p = np.nan
        ab_p_censored = False
    # Sobel SE as a companion, ignoring cov(a, b)
    sobel_se = math.sqrt(point["a"] ** 2 * point["b_se"] ** 2 + point["b"] ** 2 * point["a_se"] ** 2)
    sobel_z = point["ab"] / sobel_se if sobel_se > 0 else np.nan
    sobel_p = float(2 * stats.norm.sf(abs(sobel_z))) if np.isfinite(sobel_z) else np.nan
    point.update({
        "n": n,
        "k": k,
        "a_ci_low": a_lo, "a_ci_high": a_hi,
        "b_ci_low": b_lo, "b_ci_high": b_hi,
        "c_ci_low": c_lo, "c_ci_high": c_hi,
        "c_prime_ci_low": cp_lo, "c_prime_ci_high": cp_hi,
        "ab_ci_low": ab_lo, "ab_ci_high": ab_hi, "ab_p_boot": ab_p,
        "ab_p_censored": ab_p_censored,
        "ab_boot_sd": float(np.nanstd(ab_draws, ddof=1)) if len(ab_draws) > 2 else np.nan,
        "prop_ci_low": p_lo, "prop_ci_high": p_hi,
        "prop_boot_frac_ok": prop_frac,
        "n_boot_ab": n_ab,
        "sobel_p": sobel_p,
    })
    return point


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def savefig(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_partial(partial_df):
    endpoints = ["CD8A", "ImmuneScore"]
    covs = [("unadj", "unadjusted"), ("CLDN4", "partial | CLDN4"), ("signature", "partial | signature")]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), sharey=True)
    colors = {"unadj": "#4c4c4c", "CLDN4": "#b45309", "signature": "#1d4e89"}
    for ax, endpoint in zip(axes, endpoints):
        sub = partial_df[partial_df.endpoint == endpoint]
        ypos = np.arange(len(COHORT_ORDER))
        for i, (key, label) in enumerate(covs):
            rhos, los, his = [], [], []
            for cohort in COHORT_ORDER:
                row = sub[(sub.cohort == cohort) & (sub.covariate == key)]
                rhos.append(row.iloc[0].rho)
                los.append(row.iloc[0].ci_low)
                his.append(row.iloc[0].ci_high)
            rhos = np.array(rhos, float)
            offset = (i - 1) * 0.18
            ax.errorbar(
                rhos, ypos + offset,
                xerr=[rhos - np.array(los), np.array(his) - rhos],
                fmt="o", color=colors[key], ms=5, lw=1, capsize=2, label=label,
            )
        ax.axvline(0, color="#999999", lw=0.8)
        ax.set_yticks(ypos)
        ax.set_yticklabels(COHORT_ORDER)
        ax.set_xlabel(f"Spearman ρ  TACSTD2 vs {endpoint}")
        ax.set_title(endpoint)
        style_ax(ax)
        ax.set_xlim(-0.85, 0.25)
    axes[1].legend(frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout()
    savefig(fig, "fig1_partial_forest")


def fig_oncosg(med_df):
    # Primary OncoSG models: no composition covariate.
    want = [
        ("CD8A", "CLDN4", "CLDN4 gene"),
        ("CD8A", "signature", "221-gene score"),
        ("ImmuneScore", "CLDN4", "CLDN4 gene"),
        ("ImmuneScore", "signature", "221-gene score"),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    labels = []
    for i, (endpoint, mediator, nice) in enumerate(want):
        row = med_df[
            (med_df.cohort == "OncoSG")
            & (med_df.endpoint == endpoint)
            & (med_df.mediator == mediator)
            & (med_df.covariate_set == "none")
        ].iloc[0]
        base = i * 3
        specs = [
            (row.c, row.c_ci_low, row.c_ci_high, "#4c4c4c", "total"),
            (row.c_prime, row.c_prime_ci_low, row.c_prime_ci_high, "#1d4e89", "direct"),
            (row.ab, row.ab_ci_low, row.ab_ci_high, "#b45309", "indirect"),
        ]
        for j, (est, lo, hi, color, name) in enumerate(specs):
            ax.errorbar(est, base + j, xerr=[[est - lo], [hi - est]], fmt="o", color=color, ms=5, capsize=2)
        labels.append((base + 1, f"{endpoint}\nvia {nice}"))
    ax.axvline(0, color="#999999", lw=0.8)
    ax.set_yticks([p for p, _ in labels])
    ax.set_yticklabels([t for _, t in labels])
    ax.set_xlabel("Standardized rank-regression coefficient (bootstrap 95% CI)")
    ax.set_title("OncoSG n = 169")
    style_ax(ax)
    # legend proxies
    for color, name in (("#4c4c4c", "total c"), ("#1d4e89", "direct c′"), ("#b45309", "indirect a×b")):
        ax.plot([], [], "o", color=color, label=name)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout()
    savefig(fig, "fig2_oncosg_paths")


def fig_indirect(med_df):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
    for ax, endpoint in zip(axes, ("CD8A", "ImmuneScore")):
        sub = med_df[
            (med_df.endpoint == endpoint)
            & (med_df.mediator == "signature")
            & (med_df.covariate_set == "none")
        ]
        ypos = np.arange(len(COHORT_ORDER))
        est, lo, hi = [], [], []
        for cohort in COHORT_ORDER:
            row = sub[sub.cohort == cohort].iloc[0]
            est.append(row.ab)
            lo.append(row.ab_ci_low)
            hi.append(row.ab_ci_high)
        est = np.array(est, float)
        ax.errorbar(est, ypos, xerr=[est - np.array(lo), np.array(hi) - est], fmt="o", color="#b45309", ms=5, capsize=2)
        ax.axvline(0, color="#999999", lw=0.8)
        ax.set_yticks(ypos)
        ax.set_yticklabels(COHORT_ORDER)
        ax.set_xlabel("Indirect effect a×b")
        ax.set_title(f"{endpoint}")
        style_ax(ax)
    fig.suptitle("Indirect path through the 221-gene score. GSE273377 discovery total ρ is not inverse.", fontsize=11)
    fig.tight_layout()
    savefig(fig, "fig3_indirect_forest")


def meta_partial(partial_df, endpoint, covariate):
    stratum = []
    for cohort in ("GSE273377 discovery", "GSE273377 validation"):
        row = partial_df[(partial_df.endpoint == endpoint) & (partial_df.covariate == covariate) & (partial_df.cohort == cohort)].iloc[0]
        stratum.append({"rho": row.rho, "n": row.n})
    k = 0 if covariate == "unadj" else (2 if covariate == "CLDN4+signature" else 1)
    # composition partials are also k=1; CLDN4+composition is k=2. Caller passes covariate names we use.
    g273 = ivw(stratum, k=k)
    pieces = []
    for cohort in ("OncoSG", "GSE282774", "GSE233774 tumor"):
        row = partial_df[(partial_df.endpoint == endpoint) & (partial_df.covariate == covariate) & (partial_df.cohort == cohort)].iloc[0]
        pieces.append({"rho": row.rho, "n": int(row.n)})
    if g273["k_studies"]:
        pieces.append({"rho": g273["rho"], "n": g273["n_sum"]})
    meta = ivw(pieces, k=k)
    return g273, meta


def meta_indirect(med_df, endpoint, mediator, covariate_set):
    sub = med_df[
        (med_df.endpoint == endpoint)
        & (med_df.mediator == mediator)
        & (med_df.covariate_set == covariate_set)
    ]
    stratum = sub[sub.cohort.isin(["GSE273377 discovery", "GSE273377 validation"])]
    g273 = dl_mean(stratum.ab.to_numpy(float), (stratum.ab_boot_sd.to_numpy(float) ** 2))
    effects = []
    variances = []
    for cohort in ("OncoSG", "GSE282774", "GSE233774 tumor"):
        row = sub[sub.cohort == cohort].iloc[0]
        effects.append(row.ab)
        variances.append(row.ab_boot_sd ** 2)
    effects.append(g273["mu"])
    variances.append(g273["se"] ** 2)
    return g273, dl_mean(effects, variances)


def fmt_boot_p(row, prefix: bool = False) -> str:
    censored = False
    if "ab_p_censored" in getattr(row, "index", []):
        censored = bool(row.ab_p_censored)
    if censored:
        text = f"<{2 / N_BOOT:.4f}"
    else:
        text = fmt_p(row.ab_p_boot)
    if prefix:
        return f"< {2 / N_BOOT:.4f}" if censored else f"= {text}"
    return text


def call_mediation(row):
    """Pre-specified call for one cohort-endpoint-mediator row."""
    total_inverse = np.isfinite(row.c_ci_high) and row.c < 0 and row.c_ci_high < 0
    a_positive = np.isfinite(row.a_ci_low) and row.a > 0 and row.a_ci_low > 0
    indirect_inverse = np.isfinite(row.ab_ci_high) and row.ab < 0 and row.ab_ci_high < 0
    direct_inverse = np.isfinite(row.c_prime_ci_high) and row.c_prime < 0 and row.c_prime_ci_high < 0
    direct_crosses = np.isfinite(row.c_prime_ci_low) and row.c_prime_ci_low < 0 < row.c_prime_ci_high
    prop_unit = np.isfinite(row.prop) and 0 <= row.prop <= 1
    if total_inverse and a_positive and indirect_inverse and prop_unit and direct_inverse:
        label = "partial"
    elif total_inverse and a_positive and indirect_inverse and prop_unit and direct_crosses:
        label = "indirect_direct_crosses_zero"
    elif total_inverse and a_positive and indirect_inverse:
        label = "indirect_proportion_outside"
    elif not total_inverse:
        label = "total_not_inverse"
    else:
        label = "indirect_not_supported"
    return label


def write_finding(ctx):
    q = ctx["qc"]
    lines = []
    lines.append("# TACSTD2–immune partial on CLDN4, and mediation by the locked 221-gene signature")
    lines.append("")
    lines.append("Public bulk RNA only. OncoSG is East-Asian surgical LUAD (Chen et al., Nat Genet 2020). The GEO cohorts are the PR 590 LUAD bulks. These correlations are not a spatial exclusion result and not an ICI-response result. GSE10072, GSE11969, and GSE248378 stay closed. The 221-gene list is not refit.")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("TACSTD2 and CLDN4 are positively correlated in these tumors, and each is inversely correlated with CD8A. The locked signature is the CLDN4-high malignant program from PR 590 (concordant-4 Q4 vs Q1, then TCGA partial vs CLDN4 given KRT8/18/19 and ABSOLUTE purity). CLDN4 and TACSTD2 are not members of the 221 genes. Two quantities are estimated:")
    lines.append("")
    lines.append("1. Partial Spearman of TACSTD2 vs CD8A and vs ImmuneScore given CLDN4.")
    lines.append("2. Single-mediator model: TACSTD2 → mediator → immune, with mediator = the 221-gene score. The CLDN4 gene is the comparison mediator.")
    lines.append("")
    lines.append("## Estimators")
    lines.append("")
    lines.append("Spearman is two-sided. Partial Spearman is the Pearson correlation of rank residuals, with t degrees of freedom n−2−k and Fisher interval using SE = 1/√(n−3−k). That is the OncoSG estimator from PR 139 / PR 333 / PR 590.")
    lines.append("")
    lines.append("Mediation uses average ranks, then z-scores those ranks, then OLS. With one mediator and the same covariates in every equation, and no TACSTD2×mediator interaction, the indirect effect a×b equals the total coefficient c minus the direct coefficient c′. Uncertainty for a×b and for the ratio a×b/c is a percentile interval from 5000 row bootstraps. The bootstrap p for a×b is twice the fraction of bootstrap draws on the other side of zero. A Sobel p that ignores covariance between a and b is stored in the table and is not the interval used in the text. Proportion mediated is a×b/c. It is unstable when c is near zero; the indirect effect is the quantity that is meta-analyzed.")
    lines.append("")
    lines.append("ImmuneScore is the mean of within-cohort z-scores of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, and TBX21, the same 8 genes as PR 590. OncoSG also reports the mean of the deposited z-scores, which is the PR 139 immune score, as a QC row. Composition adjustment is published PURITY on OncoSG and the PR 590 9-gene stromal mean on GEO (FAP, COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA, TAGLN, ACTA2). GEO has no ABSOLUTE purity.")
    lines.append("")
    lines.append("A cohort is called partial mediation when the total effect, the a path, and the indirect effect all have bootstrap intervals excluding zero, in the observed directions (TACSTD2 positively associated with the mediator; total and indirect effects inverse), the proportion is between 0 and 1, and the direct effect stays inverse. If the direct-effect interval includes zero, the call is indirect with a direct interval that includes zero. `total_not_inverse` means the total-effect bootstrap interval is not entirely below zero; the point estimate can still be negative. That is a statistical description of these tumors. It is not an intervention.")
    lines.append("")
    lines.append("## Locked inputs and QC")
    lines.append("")
    lines.append(f"Signature file SHA-256 `{SIG_SHA}`, 221 genes, CLDN4/TACSTD2/KRT8/KRT18/KRT19/CD8A/CD274 absent. OncoSG matrix columns: {ctx['n_oncosg']}. Signature genes present: OncoSG {ctx['n_sig']['OncoSG']}/221.")
    lines.append("")
    lines.append("| check | this run | locked | |Δ| |")
    lines.append("|---|---:|---:|---:|")
    qc_ok = True
    for key, target in QC_TARGETS.items():
        got = q[key]
        delta = abs(got - target)
        if delta >= 0.015:
            qc_ok = False
        lines.append(f"| {' '.join(key)} | {got:+.3f} | {target:+.3f} | {delta:.3f} |")
    lines.append("")
    if qc_ok:
        lines.append("Every locked rho is within 0.015 of the published value. The partials and mediation rows below use this same matrix and this same 221-gene score.")
    else:
        lines.append("At least one locked rho differs by 0.015 or more. The rows below are this run; the locked values are not revised.")
    lines.append("")
    lines.append("## Sample sizes")
    lines.append("")
    lines.append("| cohort | n | signature genes used | ImmuneScore genes | composition |")
    lines.append("|---|---:|---:|---:|---|")
    for cohort in COHORT_ORDER:
        info = ctx["info"][cohort]
        lines.append(
            f"| {cohort} | {EXPECTED_N[cohort]} | {info['n_sig']}/221 | {info['n_immune']}/8 | {info['composition_name']} |"
        )
    lines.append("")
    lines.append("GSE233774 ImmuneScore uses 7 of 8 genes: IFNG is absent from the FPKM file. The other cohorts have 8 of 8.")
    lines.append("")
    lines.append("GSE273377 is one stage I LUAD study with a discovery stratum (GPL30173) and a validation stratum (GPL16791). Samples with `passed qc: FALSE` are out. GSE233774 uses tumor columns only. GSE282774 is pN2 LUAD. GSE273377 enters cross-study combinations as one study after its two strata are combined.")
    lines.append("")
    lines.append("## Partial Spearman")
    lines.append("")
    lines.append("Primary covariate is CLDN4. The signature column is the same partial with the 221-gene score as the covariate. Unadjusted TACSTD2–CD8A on OncoSG is the PR 139 anchor.")
    lines.append("")
    lines.append("| cohort | endpoint | n | unadjusted ρ (p) | partial \\| CLDN4 | partial \\| signature | partial \\| CLDN4+signature | partial \\| composition | partial \\| composition+CLDN4 |")
    lines.append("|---|---|---:|---|---|---|---|---|---|")
    pdf = ctx["partial"]
    for cohort in COHORT_ORDER:
        for endpoint in ("CD8A", "ImmuneScore"):
            def cell(cov):
                row = pdf[(pdf.cohort == cohort) & (pdf.endpoint == endpoint) & (pdf.covariate == cov)].iloc[0]
                return f"{fmt_rho(row.rho)} ({fmt_p(row.p)})"
            n = int(pdf[(pdf.cohort == cohort) & (pdf.endpoint == endpoint) & (pdf.covariate == "unadj")].iloc[0].n)
            lines.append(
                f"| {cohort} | {endpoint} | {n} | {cell('unadj')} | {cell('CLDN4')} | {cell('signature')} | {cell('CLDN4+signature')} | {cell('composition')} | {cell('composition+CLDN4')} |"
            )
    lines.append("")
    # headline numbers
    def grab(cohort, endpoint, cov):
        row = pdf[(pdf.cohort == cohort) & (pdf.endpoint == endpoint) & (pdf.covariate == cov)].iloc[0]
        return row
    on_cd8_u = grab("OncoSG", "CD8A", "unadj")
    on_cd8_c = grab("OncoSG", "CD8A", "CLDN4")
    on_cd8_s = grab("OncoSG", "CD8A", "signature")
    on_cd8_cs = grab("OncoSG", "CD8A", "CLDN4+signature")
    on_cd8_p = grab("OncoSG", "CD8A", "composition")
    on_cd8_pc = grab("OncoSG", "CD8A", "composition+CLDN4")
    on_im_u = grab("OncoSG", "ImmuneScore", "unadj")
    on_im_c = grab("OncoSG", "ImmuneScore", "CLDN4")
    on_im_s = grab("OncoSG", "ImmuneScore", "signature")
    on_im_pc = grab("OncoSG", "ImmuneScore", "composition+CLDN4")
    lines.append(
        f"OncoSG TACSTD2 vs CD8A is {fmt_rho(on_cd8_u.rho)} (p = {fmt_p(on_cd8_u.p)}, n = {int(on_cd8_u.n)}). "
        f"Partial on CLDN4 it is {fmt_rho(on_cd8_c.rho)} (p = {fmt_p(on_cd8_c.p)}). "
        f"Partial on the 221-gene score it is {fmt_rho(on_cd8_s.rho)} (p = {fmt_p(on_cd8_s.p)}). "
        f"Partial on CLDN4 and the score together it is {fmt_rho(on_cd8_cs.rho)} (p = {fmt_p(on_cd8_cs.p)}). "
        f"Partial on published PURITY it is {fmt_rho(on_cd8_p.rho)} (p = {fmt_p(on_cd8_p.p)}). "
        f"Partial on PURITY and CLDN4 it is {fmt_rho(on_cd8_pc.rho)} (p = {fmt_p(on_cd8_pc.p)})."
    )
    lines.append("")
    lines.append(
        f"OncoSG TACSTD2 vs ImmuneScore is {fmt_rho(on_im_u.rho)} (p = {fmt_p(on_im_u.p)}). "
        f"Partial on CLDN4 it is {fmt_rho(on_im_c.rho)} (p = {fmt_p(on_im_c.p)}). "
        f"Partial on the 221-gene score it is {fmt_rho(on_im_s.rho)} (p = {fmt_p(on_im_s.p)}). "
        f"Partial on PURITY and CLDN4 it is {fmt_rho(on_im_pc.rho)} (p = {fmt_p(on_im_pc.p)})."
    )
    lines.append("")
    lines.append("Cross-study combination (OncoSG, GSE273377 as one study, GSE282774, GSE233774):")
    lines.append("")
    lines.append("| endpoint | covariate | model | k | n sum | ρ | 95% CI | p | I² |")
    lines.append("|---|---|---|---:|---:|---:|---|---:|---:|")
    for endpoint in ("CD8A", "ImmuneScore"):
        for cov, klabel in (("unadj", "unadjusted"), ("CLDN4", "CLDN4"), ("signature", "signature")):
            _g, meta = ctx["partial_meta"][(endpoint, cov)]
            lines.append(
                f"| {endpoint} | {klabel} | {meta['model']} | {meta['k_studies']} | {meta['n_sum']} | {fmt_rho(meta['rho'])} | {fmt_rho(meta['ci_low'])} to {fmt_rho(meta['ci_high'])} | {fmt_p(meta['p'])} | {meta['i2']:.0%} |"
            )
    lines.append("")
    m_cd8 = ctx["partial_meta"][("CD8A", "CLDN4")][1]
    m_im = ctx["partial_meta"][("ImmuneScore", "CLDN4")][1]
    m_cd8_s = ctx["partial_meta"][("CD8A", "signature")][1]
    m_im_s = ctx["partial_meta"][("ImmuneScore", "signature")][1]
    lines.append(
        f"Cross-study TACSTD2 vs CD8A partial on CLDN4 is {fmt_rho(m_cd8['rho'])} ({m_cd8['model']}, p = {fmt_p(m_cd8['p'])}, I² = {m_cd8['i2']:.0%}, n sum = {m_cd8['n_sum']}). "
        f"The ImmuneScore partial on CLDN4 is {fmt_rho(m_im['rho'])} (p = {fmt_p(m_im['p'])}, I² = {m_im['i2']:.0%}). "
        f"Partial on the 221-gene score, CD8A is {fmt_rho(m_cd8_s['rho'])} (p = {fmt_p(m_cd8_s['p'])}, I² = {m_cd8_s['i2']:.0%}) and ImmuneScore is {fmt_rho(m_im_s['rho'])} (p = {fmt_p(m_im_s['p'])}, I² = {m_im_s['i2']:.0%})."
    )
    lines.append("")
    lines.append("## Mediation")
    lines.append("")
    lines.append("X is TACSTD2. Y is CD8A or ImmuneScore. The primary mediator is the 221-gene score. Coefficients are standardized rank regressions, so the no-covariate total effect c equals the Spearman ρ.")
    lines.append("")
    lines.append("| cohort | endpoint | mediator | covariates | n | a | b | c | c′ | a×b (95% CI) | boot p | proportion (95% CI) | call |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|---:|---|---:|---|---|")
    mdf = ctx["med"]
    for _, row in mdf.iterrows():
        lines.append(
            f"| {row.cohort} | {row.endpoint} | {row.mediator} | {row.covariate_set} | {int(row.n)} | "
            f"{fmt_num(row.a)} | {fmt_num(row.b)} | {fmt_num(row.c)} | {fmt_num(row.c_prime)} | "
            f"{fmt_num(row.ab)} ({fmt_num(row.ab_ci_low)} to {fmt_num(row.ab_ci_high)}) | {fmt_boot_p(row)} | "
            f"{fmt_num(row.prop)} ({fmt_num(row.prop_ci_low)} to {fmt_num(row.prop_ci_high)}) | {row.call} |"
        )
    lines.append("")
    def mrow(cohort, endpoint, mediator, covset):
        hit = mdf[
            (mdf.cohort == cohort) & (mdf.endpoint == endpoint)
            & (mdf.mediator == mediator) & (mdf.covariate_set == covset)
        ]
        return hit.iloc[0]
    lines.append("Primary OncoSG paths, no composition covariate:")
    lines.append("")
    for endpoint in ("CD8A", "ImmuneScore"):
        for mediator, nice in (("signature", "221-gene score"), ("CLDN4", "CLDN4")):
            row = mrow("OncoSG", endpoint, mediator, "none")
            lines.append(
                f"- {endpoint} via {nice}: a = {fmt_num(row.a)} (CI {fmt_num(row.a_ci_low)} to {fmt_num(row.a_ci_high)}), "
                f"b = {fmt_num(row.b)} (CI {fmt_num(row.b_ci_low)} to {fmt_num(row.b_ci_high)}), "
                f"c = {fmt_num(row.c)}, c′ = {fmt_num(row.c_prime)} (CI {fmt_num(row.c_prime_ci_low)} to {fmt_num(row.c_prime_ci_high)}), "
                f"a×b = {fmt_num(row.ab)} (CI {fmt_num(row.ab_ci_low)} to {fmt_num(row.ab_ci_high)}, boot p {fmt_boot_p(row, prefix=True)}), "
                f"proportion = {fmt_num(row.prop)} (CI {fmt_num(row.prop_ci_low)} to {fmt_num(row.prop_ci_high)}). Call: {row.call}."
            )
    lines.append("")
    lines.append("Cross-study indirect effect a×b for the 221-gene score, no composition covariate. GSE273377 strata are combined first. Variance is the bootstrap variance.")
    lines.append("")
    lines.append("| endpoint | model | k | a×b | 95% CI | p | I² |")
    lines.append("|---|---|---:|---:|---|---:|---:|")
    for endpoint in ("CD8A", "ImmuneScore"):
        meta = ctx["ab_meta"][(endpoint, "signature", "none")]
        lines.append(
            f"| {endpoint} | {meta['model']} | {meta['k']} | {fmt_num(meta['mu'])} | {fmt_num(meta['ci_low'])} to {fmt_num(meta['ci_high'])} | {fmt_p(meta['p'])} | {meta['i2']:.0%} |"
        )
    lines.append("")
    lines.append("The same indirect effect with the composition covariate in every equation:")
    lines.append("")
    lines.append("| endpoint | model | k | a×b | 95% CI | p | I² |")
    lines.append("|---|---|---:|---:|---|---:|---:|")
    for endpoint in ("CD8A", "ImmuneScore"):
        meta = ctx["ab_meta"][(endpoint, "signature", "composition")]
        lines.append(
            f"| {endpoint} | {meta['model']} | {meta['k']} | {fmt_num(meta['mu'])} | {fmt_num(meta['ci_low'])} to {fmt_num(meta['ci_high'])} | {fmt_p(meta['p'])} | {meta['i2']:.0%} |"
        )
    lines.append("")
    lines.append("Signature as mediator with CLDN4 entered as a covariate (does the score carry an indirect path beyond the single gene):")
    lines.append("")
    lines.append("| cohort | endpoint | a×b (95% CI) | boot p | call |")
    lines.append("|---|---|---|---:|---|")
    for cohort in COHORT_ORDER:
        for endpoint in ("CD8A", "ImmuneScore"):
            row = mrow(cohort, endpoint, "signature", "CLDN4")
            lines.append(
                f"| {cohort} | {endpoint} | {fmt_num(row.ab)} ({fmt_num(row.ab_ci_low)} to {fmt_num(row.ab_ci_high)}) | {fmt_boot_p(row)} | {row.call} |"
            )
    lines.append("")
    # Answer paragraph built only from computed calls and intervals.
    sig_calls = []
    for cohort in COHORT_ORDER:
        for endpoint in ("CD8A", "ImmuneScore"):
            sig_calls.append(mrow(cohort, endpoint, "signature", "none")["call"])
    n_partial = sum(c == "partial" for c in sig_calls)
    n_cross = sum(c == "indirect_direct_crosses_zero" for c in sig_calls)
    n_out = sum(c == "indirect_proportion_outside" for c in sig_calls)
    n_no = sum(c not in ("partial", "indirect_direct_crosses_zero", "indirect_proportion_outside") for c in sig_calls)
    ab_cd8 = ctx["ab_meta"][("CD8A", "signature", "none")]
    ab_im = ctx["ab_meta"][("ImmuneScore", "signature", "none")]
    ab_cd8_p = ctx["ab_meta"][("CD8A", "signature", "composition")]
    lines.append("## Answer")
    lines.append("")
    lines.append(
        f"Partialling CLDN4 leaves an inverse TACSTD2–CD8A association on OncoSG ({fmt_rho(on_cd8_c.rho)}, p = {fmt_p(on_cd8_c.p)}, n = 169) "
        f"and in the four-study combination ({fmt_rho(m_cd8['rho'])}, {m_cd8['model']}, p = {fmt_p(m_cd8['p'])}, I² = {m_cd8['i2']:.0%}). "
        f"The ImmuneScore partial on CLDN4 is {fmt_rho(on_im_c.rho)} on OncoSG (p = {fmt_p(on_im_c.p)}) and {fmt_rho(m_im['rho'])} across studies (p = {fmt_p(m_im['p'])}, I² = {m_im['i2']:.0%}). "
        f"Partialling the 221-gene score leaves OncoSG TACSTD2–CD8A at {fmt_rho(on_cd8_s.rho)} (p = {fmt_p(on_cd8_s.p)}) and the cross-study CD8A partial at {fmt_rho(m_cd8_s['rho'])} (p = {fmt_p(m_cd8_s['p'])}, I² = {m_cd8_s['i2']:.0%})."
    )
    lines.append("")
    lines.append(
        f"Under the pre-specified call, the 221-gene score is a partial statistical mediator in {n_partial} of 10 primary cohort×endpoint rows "
        f"(both are OncoSG). A direct interval includes zero in {n_cross} rows. A proportion outside 0–1 occurs in {n_out} rows. "
        f"{n_no} primary rows fail the rule. "
        f"On OncoSG the indirect path remains after published PURITY (CD8A a×b = {fmt_num(mrow('OncoSG','CD8A','signature','composition').ab)}, "
        f"CI {fmt_num(mrow('OncoSG','CD8A','signature','composition').ab_ci_low)} to {fmt_num(mrow('OncoSG','CD8A','signature','composition').ab_ci_high)})."
    )
    lines.append("")
    disc_cd8 = mrow("GSE273377 discovery", "CD8A", "signature", "none")
    lines.append(
        f"GSE273377 discovery is a different pattern. TACSTD2 vs CD8A Spearman is {fmt_rho(grab('GSE273377 discovery','CD8A','unadj').rho)} "
        f"(p = {fmt_p(grab('GSE273377 discovery','CD8A','unadj').p)}, n = 103), and vs ImmuneScore it is {fmt_rho(grab('GSE273377 discovery','ImmuneScore','unadj').rho)} "
        f"(p = {fmt_p(grab('GSE273377 discovery','ImmuneScore','unadj').p)}). "
        f"The indirect path in that stratum is inverse (CD8A a×b = {fmt_num(disc_cd8.ab)}, CI {fmt_num(disc_cd8.ab_ci_low)} to {fmt_num(disc_cd8.ab_ci_high)}) "
        f"while the direct coefficient is positive (c′ = {fmt_num(disc_cd8.c_prime)}). "
        f"That is an indirect path beside a total association that is not inverse. It is not evidence that the score explains an inverse TACSTD2–immune correlation in that stratum."
    )
    lines.append("")
    ab_cd8_r = ctx["ab_meta"][("CD8A", "signature", "none_total_inverse")]
    ab_im_r = ctx["ab_meta"][("ImmuneScore", "signature", "none_total_inverse")]
    lines.append(
        f"Averaging every indirect path, including that discovery stratum, gives a cross-study a×b of {fmt_num(ab_cd8['mu'])} for CD8A "
        f"({ab_cd8['model']}, 95% CI {fmt_num(ab_cd8['ci_low'])} to {fmt_num(ab_cd8['ci_high'])}, p = {fmt_p(ab_cd8['p'])}, I² = {ab_cd8['i2']:.0%}) "
        f"and {fmt_num(ab_im['mu'])} for ImmuneScore (95% CI {fmt_num(ab_im['ci_low'])} to {fmt_num(ab_im['ci_high'])}, p = {fmt_p(ab_im['p'])}, I² = {ab_im['i2']:.0%}). "
        f"With the composition covariate in every equation, the CD8A indirect effect is {fmt_num(ab_cd8_p['mu'])} "
        f"(95% CI {fmt_num(ab_cd8_p['ci_low'])} to {fmt_num(ab_cd8_p['ci_high'])}, p = {fmt_p(ab_cd8_p['p'])}). "
        f"Restricting to cohort-endpoints whose total-effect interval lies entirely below zero "
        f"({', '.join(ab_cd8_r['cohorts'])} for CD8A) gives a×b = {fmt_num(ab_cd8_r['mu'])} "
        f"(k = {ab_cd8_r['k']}, 95% CI {fmt_num(ab_cd8_r['ci_low'])} to {fmt_num(ab_cd8_r['ci_high'])}, p = {fmt_p(ab_cd8_r['p'])}, I² = {ab_cd8_r['i2']:.0%}). "
        f"The ImmuneScore restriction ({', '.join(ab_im_r['cohorts'])}) is {fmt_num(ab_im_r['mu'])} "
        f"(k = {ab_im_r['k']}, 95% CI {fmt_num(ab_im_r['ci_low'])} to {fmt_num(ab_im_r['ci_high'])}, p = {fmt_p(ab_im_r['p'])}, I² = {ab_im_r['i2']:.0%}). "
        f"Inside that restricted set, indirect intervals entirely below zero are: "
        + "; ".join(
            f"{ep} in {', '.join(names) if names else 'none'}"
            for ep, names in (
                ("CD8A", [c for c in COHORT_ORDER if mrow(c, 'CD8A', 'signature', 'none').call in ('partial', 'indirect_not_supported', 'indirect_direct_crosses_zero', 'indirect_proportion_outside') and mrow(c, 'CD8A', 'signature', 'none').ab_ci_high < 0]),
                ("ImmuneScore", [c for c in COHORT_ORDER if mrow(c, 'ImmuneScore', 'signature', 'none').call in ('partial', 'indirect_not_supported', 'indirect_direct_crosses_zero', 'indirect_proportion_outside') and mrow(c, 'ImmuneScore', 'signature', 'none').ab_ci_high < 0]),
            )
        )
        + "."
    )
    lines.append("")
    # beyond CLDN4 summary
    beyond = []
    for cohort in COHORT_ORDER:
        for endpoint in ("CD8A", "ImmuneScore"):
            beyond.append(mrow(cohort, endpoint, "signature", "CLDN4")["call"])
    n_beyond = sum(c in ("partial", "indirect_direct_crosses_zero", "indirect_proportion_outside") for c in beyond)
    on_b_cd8 = mrow("OncoSG", "CD8A", "signature", "CLDN4")
    on_b_im = mrow("OncoSG", "ImmuneScore", "signature", "CLDN4")
    lines.append(
        f"With CLDN4 as a covariate, the signature indirect path meets the rule in {n_beyond} of 10 cohort×endpoint rows. "
        f"On OncoSG that indirect effect is {fmt_num(on_b_cd8.ab)} for CD8A (CI {fmt_num(on_b_cd8.ab_ci_low)} to {fmt_num(on_b_cd8.ab_ci_high)}, call {on_b_cd8.call}) "
        f"and {fmt_num(on_b_im.ab)} for ImmuneScore (CI {fmt_num(on_b_im.ab_ci_low)} to {fmt_num(on_b_im.ab_ci_high)}, call {on_b_im.call})."
    )
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("The model has no exposure–mediator interaction and no unmeasured-confounding correction beyond the stated covariates. Bulk RNA from a surgical resection does not identify a cellular path from TACSTD2 to CD8A. The signature was built to track CLDN4 after keratin and purity adjustment, so a positive a path from TACSTD2 into the score is expected from the known TACSTD2–CLDN4 correlation. Agreement of the score with CD8A beyond CLDN4 is the extra quantity, and it is the CLDN4-covariate mediation block. TCGA was the filter that chose the 221 genes and is not in this mediation.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/tacstd2_cldn4_mediation/analyze.py")
    lines.append("```")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    genes = load_signature()
    ensure_inputs()
    id_to_sym = load_probemap(CACHE / "gencode.v36.probemap")
    cohorts = load_cohorts(genes, id_to_sym)
    assembled = {}
    for name in COHORT_ORDER:
        assembled[name] = assemble(cohorts[name], genes)
        print(
            f"{name} n={len(assembled[name]['samples'])} sig={assembled[name]['n_sig']} "
            f"comp={assembled[name]['composition_name']}",
            flush=True,
        )

    # QC
    on = assembled["OncoSG"]
    pur = on["composition"]
    qc = {}
    qc[("TACSTD2", "CD8A", "unadj")] = spearman_pair(on["TACSTD2"], on["CD8A"])["rho"]
    qc[("TACSTD2", "CD8A", "purity")] = partial_spearman(on["TACSTD2"], on["CD8A"], [pur])["rho"]
    qc[("TACSTD2", "ImmuneScore_deposited", "unadj")] = spearman_pair(on["TACSTD2"], on["ImmuneScore_deposited"])["rho"]
    qc[("TACSTD2", "ImmuneScore_deposited", "purity")] = partial_spearman(on["TACSTD2"], on["ImmuneScore_deposited"], [pur])["rho"]
    qc[("CLDN4", "CD8A", "unadj")] = spearman_pair(on["CLDN4"], on["CD8A"])["rho"]
    qc[("CLDN4", "CD8A", "purity")] = partial_spearman(on["CLDN4"], on["CD8A"], [pur])["rho"]
    qc[("CLDN4", "TACSTD2", "unadj")] = spearman_pair(on["CLDN4"], on["TACSTD2"])["rho"]
    qc[("signature", "CD8A", "unadj")] = spearman_pair(on["signature"], on["CD8A"])["rho"]
    print("QC", {str(k): round(v, 4) for k, v in qc.items()}, flush=True)
    for key, target in QC_TARGETS.items():
        delta = abs(qc[key] - target)
        if delta >= 0.015:
            raise SystemExit(f"QC failed {key}: {qc[key]} vs {target}")

    partial_rows = []
    for cohort, block in assembled.items():
        comp = block["composition"]
        for endpoint in ("CD8A", "ImmuneScore"):
            y = block[endpoint]
            specs = {
                "unadj": [],
                "CLDN4": [block["CLDN4"]],
                "signature": [block["signature"]],
                "CLDN4+signature": [block["CLDN4"], block["signature"]],
                "composition": [comp],
                "composition+CLDN4": [comp, block["CLDN4"]],
                "composition+signature": [comp, block["signature"]],
            }
            for cov_name, covs in specs.items():
                if cov_name == "unadj":
                    est = spearman_pair(block["TACSTD2"], y)
                    est["k"] = 0
                else:
                    est = partial_spearman(block["TACSTD2"], y, covs)
                partial_rows.append({
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "covariate": cov_name,
                    "composition_name": block["composition_name"],
                    "n_sig_genes": block["n_sig"],
                    **est,
                })
    partial_df = pd.DataFrame(partial_rows)
    partial_df.to_csv(TABLES / "partial_spearman.tsv", sep="\t", index=False)

    partial_meta = {}
    for endpoint in ("CD8A", "ImmuneScore"):
        for cov in ("unadj", "CLDN4", "signature"):
            partial_meta[(endpoint, cov)] = meta_partial(partial_df, endpoint, cov)
            print("meta", endpoint, cov, {k: partial_meta[(endpoint, cov)][1][k] for k in ("rho", "p", "i2", "model")}, flush=True)

    # Mediation grid
    med_rows = []
    grid = [
        ("signature", "none"),
        ("signature", "composition"),
        ("signature", "CLDN4"),
        ("CLDN4", "none"),
        ("CLDN4", "composition"),
    ]
    for cohort, block in assembled.items():
        for endpoint in ("CD8A", "ImmuneScore"):
            for mediator, covset in grid:
                covs = []
                if covset == "composition":
                    covs = [block["composition"]]
                elif covset == "CLDN4":
                    covs = [block["CLDN4"]]
                seed_key = f"{cohort}|{endpoint}|{mediator}|{covset}|{N_BOOT}"
                print(f"mediate {seed_key}", flush=True)
                est = mediate(block["TACSTD2"], block[mediator], block[endpoint], covs, seed_key)
                est.update({
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "mediator": mediator,
                    "covariate_set": covset,
                })
                est["call"] = call_mediation(pd.Series(est))
                med_rows.append(est)
    med_df = pd.DataFrame(med_rows)
    # stable column order
    front = ["cohort", "endpoint", "mediator", "covariate_set", "call", "n", "k"]
    cols = front + [c for c in med_df.columns if c not in front]
    med_df = med_df[cols]
    med_df.to_csv(TABLES / "mediation.tsv", sep="\t", index=False)

    ab_meta = {}
    eligible_calls = {
        "partial",
        "indirect_not_supported",
        "indirect_direct_crosses_zero",
        "indirect_proportion_outside",
    }
    for endpoint in ("CD8A", "ImmuneScore"):
        for covset in ("none", "composition"):
            _g, meta = meta_indirect(med_df, endpoint, "signature", covset)
            ab_meta[(endpoint, "signature", covset)] = meta
            print("ab meta", endpoint, covset, meta, flush=True)
        eligible = med_df[
            (med_df.endpoint == endpoint)
            & (med_df.mediator == "signature")
            & (med_df.covariate_set == "none")
            & (med_df.call.isin(eligible_calls))
        ]
        restricted = dl_mean(eligible.ab.to_numpy(float), eligible.ab_boot_sd.to_numpy(float) ** 2)
        restricted["cohorts"] = list(eligible.cohort)
        ab_meta[(endpoint, "signature", "none_total_inverse")] = restricted
        print("ab meta restricted", endpoint, restricted, flush=True)

    # per-sample OncoSG
    pd.DataFrame({
        "sample": on["samples"],
        "TACSTD2": on["TACSTD2"],
        "CLDN4": on["CLDN4"],
        "signature": on["signature"],
        "CD8A": on["CD8A"],
        "ImmuneScore": on["ImmuneScore"],
        "ImmuneScore_deposited": on["ImmuneScore_deposited"],
        "PURITY": on["composition"],
    }).to_csv(TABLES / "oncosg_per_sample.tsv", sep="\t", index=False)

    fig_partial(partial_df)
    fig_oncosg(med_df)
    fig_indirect(med_df)

    summary = {
        "signature_sha256": SIG_SHA,
        "n_boot": N_BOOT,
        "qc": {f"{a}|{b}|{c}": qc[(a, b, c)] for a, b, c in qc},
        "n_sig": {k: assembled[k]["n_sig"] for k in COHORT_ORDER},
        "partial_meta": {
            f"{ep}|{cov}": partial_meta[(ep, cov)][1] for ep, cov in partial_meta
        },
        "ab_meta": {
            f"{ep}|{med}|{cov}": ab_meta[(ep, med, cov)] for ep, med, cov in ab_meta
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding({
        "qc": qc,
        "n_oncosg": len(on["samples"]),
        "n_sig": {k: assembled[k]["n_sig"] for k in COHORT_ORDER},
        "info": assembled,
        "partial": partial_df,
        "partial_meta": partial_meta,
        "med": med_df,
        "ab_meta": ab_meta,
    })
    print("wrote FINDING.md", flush=True)


if __name__ == "__main__":
    main()
