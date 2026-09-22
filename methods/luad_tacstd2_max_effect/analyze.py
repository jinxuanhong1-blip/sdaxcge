#!/usr/bin/env python3
"""OncoSG + GEO LUAD max-effect sweep for TACSTD2.

Pre-declared before any maximum was read. The script does not edit the grid
after computing it.

Story-fit (required to enter the maximum)
------------------------------------------
- TACSTD2 vs CD8 / ImmuneScore: Spearman rho < 0 (TACSTD2-high, immune-lower).
- TJ / junction GSEA: NES > 0 in the TACSTD2-high direction.
Opposite-sign results are stored and are not eligible for the maximum.

Rho grid
--------
Cohorts (strict LUAD, headline n >= 50): OncoSG, GSE31210, GSE72094,
GSE68465, GSE273377 discovery, GSE273377 validation, GSE282774.
Computed but ineligible for the headline: GSE19804 (metadata says lung
cancer, not LUAD) and GSE233774 tumors (n = 30).
Closed, not downloaded: GSE10072, GSE11969, GSE248378.

CD8 endpoints: CD8A, CD8B, z-mean of CD8A+CD8B, cytotoxic z-mean of
CD8A, GZMA, GZMB, PRF1, NKG7, GNLY (CD8A required, at least 4 of 6).
ImmuneScore endpoints: immune8 z-mean (CD8A, GZMA, GZMB, IFNG, EOMES,
CXCL9, CXCL10, TBX21; at least 6 of 8) and Yoshihara Immune141 ssGSEA
(GEO only; OncoSG is a z-score deposit, so ESTIMATE is not run there).
Ayers GEP18 is computed as a companion and is not an ImmuneScore endpoint.

Subsets: all samples; high tumor content (OncoSG published PURITY at or
above the median; GEO stromal score at or below the median). Arrays also
use a cosine-purity subset when at least 50 samples have TumorPurity in
[0, 1]. Low-purity-only subsets are not in the grid.

Covariates: none; composition (OncoSG published PURITY, array ESTIMATE
StromalScore, RNA-seq 9-gene stromal z-mean); array cosine TumorPurity
and ESTIMATEScore, eligible only for CD8 endpoints (those axes include
the immune program).

Objective: maximum |rho| among story-fit cells, and maximum |DerSimonian–
Laird meta rho| among specs with at least 3 negative cohorts.

NES grid
--------
Gene sets are the name-filtered TJ/junction collection in
data/junction_sets.gmt (A8 MSigDB lists plus Enrichr KEGG 2021, GO 2023,
Reactome 2022, Hallmark 2020). Include names that match tight junction,
apical junction, occluding junction, bicellular/tricellular tight junction,
or cell-cell junction. Exclude gap junction, neuromuscular, keratin,
desmosome, focal adhesion, adherens, septate, synapse. Set size 15–500.
Contrasts: Welch t Q4 vs Q1, outer tercile, outer decile (arm floors
below), and Spearman of each gene vs TACSTD2. TACSTD2 is removed from
the ranked list. GSEA is weighted KS (p = 1), 1000 gene-set permutations,
seed 42. Headline cohorts are strict LUAD with n >= 50.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
CACHE = Path("/tmp/luad_max_cache")
RNA = CACHE / "rnaseq"
GMT_PATH = DATA / "junction_sets.gmt"

UA = "sdaxcge-luad-tacstd2-max-effect/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
DATAHUB = (
    "165bd77077b03038f9c2ee104959eb474770b2a9"
)
ONCOSG_URL = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/{DATAHUB}"
    "/public/luad_oncosg_2020"
)
ONCOSG_SHA = {
    "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt":
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    "data_clinical_sample.txt":
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
}

EST_A = 0.6049872018
EST_B = 0.0001467884
SEED = 42
NPERM = 1000
MIN_SET = 15
MAX_SET = 500
MIN_N = 50
MIN_ARM_Q = 12
MIN_ARM_DECILE = 15

IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
CYT6 = ["CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY"]
GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
STROMAL9 = ["FAP", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "TAGLN", "ACTA2"]

INC_SET = re.compile(
    r"tight[\s_-]*junction|apical[\s_-]*junction|occluding[\s_-]*junction|"
    r"bicellular[\s_-]*tight|tricellular[\s_-]*tight|cell[\s_-]*cell[\s_-]*junction",
    re.I,
)
EXC_SET = re.compile(
    r"gap[\s_-]*junction|neuromuscular|keratin|desmosom|focal[\s_-]*adhesion|"
    r"adherens|septate|synapse",
    re.I,
)

# Locked anchors from PR 139 / PR 736 / PR 235. Fail the run if OncoSG or
# GSE31210 CD8A leaves this band (loader bug). Other cohorts only warn.
ANCHOR = {
    ("OncoSG", "CD8A", "all", "none"): -0.380,
    ("OncoSG", "CD8A", "all", "composition"): -0.309,
    ("OncoSG", "immune8_deposited_mean", "all", "none"): -0.387,
    ("OncoSG", "immune8_deposited_mean", "all", "composition"): -0.318,
    ("GSE31210", "CD8A", "all", "none"): -0.289,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int = 1000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"cached {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"GET {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    if dest.stat().st_size < min_bytes:
        raise SystemExit(f"download too small: {dest}")
    return dest


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), np.asarray(Z, dtype=float)])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def partial_spearman(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8 or np.nanstd(z[m]) == 0:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m]).reshape(-1, 1)
    rx = residualize(xr, zr)
    ry = residualize(yr, zr)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan, n
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def zmean(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = []
    cols = []
    for g in genes:
        if g not in expr.index:
            continue
        v = expr.loc[g].to_numpy(dtype=float)
        if np.isfinite(v).mean() < 0.8:
            continue
        sd = np.nanstd(v)
        if not np.isfinite(sd) or sd == 0:
            continue
        cols.append((v - np.nanmean(v)) / sd)
        present.append(g)
    if not cols:
        return pd.Series(np.nan, index=expr.columns), []
    return pd.Series(np.nanmean(np.vstack(cols), axis=0), index=expr.columns), present


def deposited_mean(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), []
    return expr.loc[present].astype(float).mean(axis=0), present


def ssgsea(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> tuple[pd.Series, int]:
    present = [g for g in genes if g in expr.index]
    if len(present) < 10:
        return pd.Series(np.nan, index=expr.columns), len(present)
    Xdf = expr
    finite = np.isfinite(Xdf.to_numpy(dtype=np.float64)).all(axis=1)
    if not bool(finite.all()):
        Xdf = Xdf.loc[finite]
        present = [g for g in present if g in Xdf.index]
        if len(present) < 10:
            return pd.Series(np.nan, index=expr.columns), len(present)
    X = Xdf.to_numpy(dtype=np.float64)
    n_genes, n_s = X.shape
    ranks = stats.rankdata(X, method="average", axis=0) * (10000.0 / n_genes)
    hitset = np.zeros(n_genes, dtype=bool)
    idx = {g: i for i, g in enumerate(Xdf.index)}
    for g in present:
        hitset[idx[g]] = True
    scores = np.empty(n_s)
    for j in range(n_s):
        order = np.argsort(-ranks[:, j], kind="mergesort")
        hit = hitset[order]
        w = np.abs(ranks[order, j]) ** tau
        w_hit = np.where(hit, w, 0.0)
        nh = float(w_hit.sum())
        nm = float((~hit).sum())
        if nh <= 0 or nm <= 0:
            scores[j] = np.nan
            continue
        ph = np.cumsum(w_hit) / nh
        pm = np.cumsum(~hit) / nm
        scores[j] = float(np.sum(ph - pm))
    out = pd.Series(scores, index=Xdf.columns)
    return out.reindex(expr.columns), len(present)


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q
    pv = p[valid]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    qv = np.empty(n)
    qv[order] = adj
    q[valid] = qv
    return q


def dl_meta(rhos, ns, k_cov: int):
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    m = np.isfinite(rhos) & (ns - 3 - k_cov > 1)
    rhos, ns = rhos[m], ns[m]
    k = int(len(rhos))
    if k == 0:
        return None
    z = np.arctanh(np.clip(rhos, -0.999, 0.999))
    v = 1.0 / (ns - 3 - k_cov)
    w = 1.0 / v
    zf = float(np.sum(w * z) / np.sum(w))
    if k == 1:
        mu, se, i2, model = zf, math.sqrt(1.0 / float(np.sum(w))), 0.0, "fixed"
    else:
        Q = float(np.sum(w * (z - zf) ** 2))
        df = k - 1
        i2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
        c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
        tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
        w2 = 1.0 / (v + tau2)
        mu = float(np.sum(w2 * z) / np.sum(w2))
        se = math.sqrt(1.0 / float(np.sum(w2)))
        model = "DL random"
    zstat = mu / se
    p = float(2 * stats.norm.sf(abs(zstat)))
    return {
        "k": k,
        "rho": float(np.tanh(mu)),
        "se_z": se,
        "p": p,
        "i2": i2,
        "model": model,
        "ci_low": float(np.tanh(mu - 1.959963984540054 * se)),
        "ci_high": float(np.tanh(mu + 1.959963984540054 * se)),
        "n_sum": int(ns.sum()),
    }


def fmt_rho(r):
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------


def _keep_set(name: str, genes: list[str]) -> bool:
    if not INC_SET.search(name) or EXC_SET.search(name):
        return False
    genes = [g for g in genes if g and g != "TACSTD2"]
    return MIN_SET <= len(set(genes)) <= MAX_SET


def _parse_enrichr_line(line: str):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        return None
    name = parts[0].strip()
    genes = [g.strip() for g in parts[2:] if g.strip()]
    return name, genes


def build_gmt() -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    a8 = CACHE / "a8_sets.json"
    if not a8.exists():
        raise SystemExit(f"missing {a8}")
    obj = json.loads(a8.read_text())
    raw = obj["sets"] if isinstance(obj, dict) and "sets" in obj else obj
    for name, genes in raw.items():
        genes = [str(g) for g in genes]
        if _keep_set(name, genes):
            sets[f"A8::{name}"] = sorted(set(genes) - {"TACSTD2"})
    enrichr_files = {
        "ENRICHR_KEGG2021": RNA / "enrichr_kegg.txt",
        "ENRICHR_GOBP2023": RNA / "enrichr_gobp.txt",
        "ENRICHR_GOCC2023": RNA / "enrichr_gocc.txt",
        "ENRICHR_REACTOME2022": RNA / "enrichr_reactome.txt",
        "ENRICHR_HALLMARK2020": RNA / "enrichr_hallmark.txt",
    }
    for src, path in enrichr_files.items():
        if not path.exists():
            raise SystemExit(f"missing {path}")
        for line in path.read_text().splitlines():
            parsed = _parse_enrichr_line(line)
            if not parsed:
                continue
            name, genes = parsed
            if _keep_set(name, genes):
                sets[f"{src}::{name}"] = sorted(set(genes) - {"TACSTD2"})
    if len(sets) < 8:
        raise SystemExit(f"junction library too small: {len(sets)}")
    GMT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with GMT_PATH.open("w") as f:
        for name in sorted(sets):
            f.write(name + "\t" + "\t".join(sets[name]) + "\n")
    print(f"junction sets: {len(sets)}", flush=True)
    for name in sorted(sets):
        print(f"  {name} n={len(sets[name])}", flush=True)
    return sets


def load_gmt() -> dict[str, list[str]]:
    sets = {}
    for line in GMT_PATH.read_text().splitlines():
        parts = line.split("\t")
        sets[parts[0]] = parts[1:]
    return sets


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def parse_geo_matrix(path: Path):
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    meta_rows = {}
    samples = None
    expr_start = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_title"):
            meta_rows["title"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_source_name_ch1"):
            meta_rows["source"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            key = None
            for v in vals:
                if v and ":" in v:
                    key = v.split(":", 1)[0].strip().lower()
                    break
            if key is None:
                key = f"char_{len(meta_rows)}"
            cleaned = []
            for v in vals:
                cleaned.append(v.split(":", 1)[1].strip() if ":" in v else v)
            k = key
            n = 2
            while k in meta_rows:
                k = f"{key}_{n}"
                n += 1
            meta_rows[k] = cleaned
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            expr_start = i
            break
    if samples is None or expr_start is None:
        raise SystemExit(f"could not parse {path}")
    meta = pd.DataFrame(meta_rows, index=samples)
    header = [x.strip().strip('"') for x in lines[expr_start].rstrip("\n").split("\t")]
    rows = []
    idx = []
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        idx.append(parts[0].strip().strip('"'))
        rows.append([float(x) if x not in ("", "NA", "null", "null") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def load_gpl_annot(path: Path) -> pd.Series:
    with gzip.open(path, "rt", errors="replace") as f:
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t"):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    cands = [c for c in df.columns if "symbol" in c.lower()]
    raw = df.set_index(id_col)[cands[0]].astype(str)
    raw = raw.replace({"nan": np.nan, "None": np.nan, "": np.nan, "---": np.nan}).dropna()
    first = raw.map(lambda x: str(x).split("///")[0].strip().split()[0])
    first = first[first.str.len() > 0]
    return first


def load_gpl15048(path: Path) -> pd.Series:
    ids, syms = [], []
    started = False
    with path.open() as f:
        for line in f:
            if line.startswith("ID\t"):
                started = True
                continue
            if not started or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or not parts[3].strip():
                continue
            sym = parts[3].strip().split()[0]
            if sym in ("---", "NA", "null"):
                continue
            ids.append(parts[0])
            syms.append(sym)
    return pd.Series(syms, index=ids)


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.to_numpy()]
    out.index = pick.index
    out = out.groupby(level=0).mean()
    return out


def ensure_inputs():
    CACHE.mkdir(parents=True, exist_ok=True)
    RNA.mkdir(parents=True, exist_ok=True)
    download(f"{ONCOSG_URL}/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", CACHE / "data_mrna.txt", 10000)
    download(f"{ONCOSG_URL}/data_clinical_sample.txt", CACHE / "data_clinical_sample.txt", 1000)
    for name, expect in ONCOSG_SHA.items():
        src = CACHE / ("data_mrna.txt" if name.startswith("data_mrna") else name)
        got = sha256_file(src)
        if got != expect:
            raise SystemExit(f"sha256 mismatch {name}: {got}")
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz",
        CACHE / "GSE31210_series_matrix.txt.gz", 1_000_000,
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz",
        CACHE / "GSE72094_series_matrix.txt.gz", 1_000_000,
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68465/matrix/GSE68465_series_matrix.txt.gz",
        CACHE / "GSE68465_series_matrix.txt.gz", 1_000_000,
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz",
        CACHE / "GSE19804_series_matrix.txt.gz", 1_000_000,
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        CACHE / "GPL570.annot.gz", 100_000,
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
        CACHE / "GPL96.annot.gz", 100_000,
    )
    if not (CACHE / "GPL15048_platform.txt").exists():
        download(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&view=data&form=text",
            CACHE / "GPL15048_platform.txt", 100_000,
        )
    pairs = {
        "gencode.v36.probemap": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
        "GSE273377_exp_count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/suppl/GSE273377_exp_count.txt.gz",
        "gse273377_gpl30173.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL30173_series_matrix.txt.gz",
        "gse273377_gpl16791.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273377/matrix/GSE273377-GPL16791_series_matrix.txt.gz",
        "GSE233774_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/suppl/GSE233774_FPKM.txt.gz",
        "gse233774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/matrix/GSE233774_series_matrix.txt.gz",
        "GSE282774_expr.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/suppl/GSE282774_N2_local_Gene_expression.csv.gz",
        "gse282774_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282774/matrix/GSE282774_series_matrix.txt.gz",
        "enrichr_kegg.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
        "enrichr_gobp.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023",
        "enrichr_gocc.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2023",
        "enrichr_reactome.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022",
        "enrichr_hallmark.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020",
    }
    for name, url in pairs.items():
        download(url, RNA / name, 1000)
    if not (CACHE / "a8_sets.json").exists():
        raise SystemExit("a8_sets.json missing from cache; place the PR 736 gene-set JSON at /tmp/luad_max_cache/a8_sets.json")


def load_oncosg():
    raw = pd.read_csv(CACHE / "data_mrna.txt", sep="\t")
    raw = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    raw = raw.drop_duplicates(subset=["Hugo_Symbol"], keep="first").set_index("Hugo_Symbol")
    expr = raw.apply(pd.to_numeric, errors="coerce")
    clin = pd.read_csv(CACHE / "data_clinical_sample.txt", sep="\t", comment="#").set_index("SAMPLE_ID")
    purity = pd.to_numeric(clin.loc[list(expr.columns), "PURITY"], errors="coerce")
    purity.index = expr.columns
    return expr, purity


def load_probemap(path: Path) -> dict[str, str]:
    id_to_sym = {}
    with path.open() as f:
        next(f)
        for line in f:
            eid, sym, *_ = line.rstrip("\n").split("\t")
            if sym:
                id_to_sym[eid] = sym
                id_to_sym[eid.split(".")[0]] = sym
    return id_to_sym


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


def parse_series_titles(path: Path) -> pd.DataFrame:
    titles = None
    char_rows = []
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
                char_rows.append([p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    data = {"title": titles}
    seen = {}
    for row in char_rows:
        if ":" in row[0]:
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


def load_gse273377(id_to_sym):
    counts = pd.read_csv(RNA / "GSE273377_exp_count.txt.gz", sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    disc = parse_series_titles(RNA / "gse273377_gpl30173.txt.gz")
    val = parse_series_titles(RNA / "gse273377_gpl16791.txt.gz")
    disc["stratum"] = "discovery"
    val["stratum"] = "validation"
    meta = pd.concat([disc, val])
    logc = to_symbols(log_cpm(counts.apply(pd.to_numeric, errors="coerce")), id_to_sym)
    out = {}
    for stratum, sub in meta.groupby("stratum"):
        qc = sub["passed qc"].astype(str).str.upper().eq("TRUE")
        ids = [i for i in sub.index[qc] if i in logc.columns]
        out[stratum] = logc.loc[:, ids]
    return out


def load_gse233774(id_to_sym):
    fpkm = pd.read_csv(RNA / "GSE233774_FPKM.txt.gz", sep="\t", index_col=0)
    fpkm.index = fpkm.index.astype(str)
    meta = parse_series_titles(RNA / "gse233774_matrix.txt.gz")
    rename = {f"{t}_FPKM": t for t in meta.index}
    fpkm = fpkm.rename(columns=rename)
    tumors = [t for t in meta.index if str(t).startswith("T") and t in fpkm.columns]
    logfpkm = np.log2(to_symbols(fpkm.apply(pd.to_numeric, errors="coerce"), id_to_sym) + 1)
    return logfpkm.loc[:, tumors]


def load_gse282774():
    raw = pd.read_csv(RNA / "GSE282774_expr.csv.gz")
    meta = parse_series_titles(RNA / "gse282774_matrix.txt.gz")
    value_cols = [c for c in raw.columns if str(c).startswith("FPKM.")]
    sample_map = {c: c.split(".", 1)[1] for c in value_cols}
    expr = raw.set_index("gene_name")[value_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby(level=0).max()
    expr.columns = [sample_map[c] for c in expr.columns]
    if "disease" in meta.columns:
        keep = [s for s in expr.columns if s in meta.index and "LUAD" in str(meta.loc[s, "disease"]).upper()]
    else:
        keep = list(expr.columns)
    return np.log2(expr.loc[:, keep] + 1)


# ---------------------------------------------------------------------------
# Cohort assembly
# ---------------------------------------------------------------------------


def add_endpoint(bucket, name, family, series, genes_used, genes_required):
    bucket.append({
        "endpoint": name,
        "family": family,
        "values": series.astype(float),
        "n_genes": len(genes_used),
        "genes_used": ",".join(genes_used),
        "n_required": len(genes_required),
    })


def prepare_endpoints(expr: pd.DataFrame, do_estimate: bool, est_genes: dict) -> tuple[list[dict], dict]:
    endpoints = []
    extra = {}
    if "CD8A" in expr.index:
        add_endpoint(endpoints, "CD8A", "CD8", expr.loc["CD8A"], ["CD8A"], ["CD8A"])
    if "CD8B" in expr.index:
        add_endpoint(endpoints, "CD8B", "CD8", expr.loc["CD8B"], ["CD8B"], ["CD8B"])
    if "CD8A" in expr.index and "CD8B" in expr.index:
        s, used = zmean(expr, ["CD8A", "CD8B"])
        add_endpoint(endpoints, "CD8A_CD8B_zmean", "CD8", s, used, ["CD8A", "CD8B"])
    cyt, cyt_used = zmean(expr, CYT6)
    if "CD8A" in cyt_used and len(cyt_used) >= 4:
        add_endpoint(endpoints, "cytotoxic6_zmean", "CD8", cyt, cyt_used, CYT6)
    imm, imm_used = zmean(expr, IMMUNE8)
    if len(imm_used) >= 6:
        add_endpoint(endpoints, "immune8_zmean", "ImmuneScore", imm, imm_used, IMMUNE8)
    dep, dep_used = deposited_mean(expr, IMMUNE8)
    if len(dep_used) >= 6:
        add_endpoint(endpoints, "immune8_deposited_mean", "ImmuneScoreAnchor", dep, dep_used, IMMUNE8)
    gep, gep_used = zmean(expr, GEP18)
    if len(gep_used) >= 12:
        add_endpoint(endpoints, "ayers_gep18_zmean", "GEP18_companion", gep, gep_used, GEP18)
    if do_estimate:
        immune_s, n_imm = ssgsea(expr, est_genes["immune"])
        stromal_s, n_str = ssgsea(expr, est_genes["stromal"])
        extra["estimate_immune"] = immune_s
        extra["estimate_stromal"] = stromal_s
        extra["estimate_immune_n"] = n_imm
        extra["estimate_stromal_n"] = n_str
        extra["estimate_score"] = stromal_s + immune_s
        if n_imm >= 100:
            endpoints.append({
                "endpoint": "estimate_immune141",
                "family": "ImmuneScore",
                "values": immune_s.astype(float),
                "n_genes": int(n_imm),
                "genes_used": f"Yoshihara_Immune141_present_{n_imm}_of_{len(est_genes['immune'])}",
                "n_required": len(est_genes["immune"]),
            })
    return endpoints, extra


def high_mask(score: pd.Series, higher_is_tumor: bool) -> pd.Series:
    v = score.astype(float)
    med = np.nanmedian(v.to_numpy(dtype=float))
    if higher_is_tumor:
        return v >= med
    return v <= med


def cosine_from_score(est_score: pd.Series) -> pd.Series:
    return pd.Series(np.cos(EST_A + EST_B * est_score.to_numpy(dtype=float)), index=est_score.index)


# ---------------------------------------------------------------------------
# Rho sweep
# ---------------------------------------------------------------------------


def rho_rows_for_cohort(cohort: dict) -> list[dict]:
    tac = cohort["tac"]
    rows = []
    subsets = {"all": pd.Series(True, index=tac.index)}
    if cohort["tumor_axis"] is not None:
        subsets["high_tumor"] = high_mask(cohort["tumor_axis"], cohort["tumor_axis_high"])
    if cohort.get("cosine") is not None:
        cos = cohort["cosine"]
        in_range = (cos >= 0) & (cos <= 1) & np.isfinite(cos)
        if int(in_range.sum()) >= MIN_N:
            med = np.nanmedian(cos[in_range].to_numpy(dtype=float))
            subsets["high_cosine_purity"] = in_range & (cos >= med)
    covariates = [("none", None, "none")]
    if cohort["composition"] is not None:
        covariates.append(("composition", cohort["composition"], cohort["composition_name"]))
    if cohort.get("cosine") is not None:
        covariates.append(("cosine_purity", cohort["cosine"], "ESTIMATE_cosine_TumorPurity"))
    if cohort.get("estimate_score") is not None:
        covariates.append(("estimate_score", cohort["estimate_score"], "ESTIMATEScore_stromal_plus_immune"))
    for ep in cohort["endpoints"]:
        for subset_name, mask in subsets.items():
            for cov_class, cov, cov_name in covariates:
                if cov_class in ("cosine_purity", "estimate_score") and ep["family"] != "CD8":
                    eligible_cov = False
                else:
                    eligible_cov = True
                idx = tac.index[mask.reindex(tac.index).fillna(False).to_numpy()]
                x = tac.loc[idx].to_numpy(dtype=float)
                y = ep["values"].reindex(idx).to_numpy(dtype=float)
                if cov is None:
                    r, p, n = spearman_pair(x, y)
                else:
                    z = cov.reindex(idx).to_numpy(dtype=float)
                    r, p, n = partial_spearman(x, y, z)
                strict = bool(cohort["strict_luad"]) and cohort["n"] >= MIN_N and n >= MIN_N
                family_ok = ep["family"] in ("CD8", "ImmuneScore")
                story = bool(np.isfinite(r) and r < 0)
                rows.append({
                    "cohort": cohort["name"],
                    "strict_luad": cohort["strict_luad"],
                    "n_cohort": cohort["n"],
                    "platform": cohort["platform"],
                    "endpoint": ep["endpoint"],
                    "endpoint_family": ep["family"],
                    "n_genes": ep["n_genes"],
                    "genes_used": ep["genes_used"],
                    "subset": subset_name,
                    "covariate_class": cov_class,
                    "covariate_name": cov_name,
                    "n": n,
                    "rho": r,
                    "p": p,
                    "story_negative": story,
                    "primary_family": bool(strict and family_ok and eligible_cov),
                    "eligible_max": bool(strict and family_ok and eligible_cov and story),
                })
    return rows


# ---------------------------------------------------------------------------
# GSEA
# ---------------------------------------------------------------------------


def split_extremes(x: pd.Series, low_q: float, high_q: float):
    lo = x.quantile(low_q)
    hi = x.quantile(high_q)
    low = x.index[x <= lo]
    high = x.index[x >= hi]
    both = high.intersection(low)
    if len(both):
        high = high.difference(both)
        low = low.difference(both)
    return high, low


def welch_rank(expr: pd.DataFrame, high, low) -> pd.Series:
    eh = expr.loc[:, high].to_numpy(dtype=np.float64)
    el = expr.loc[:, low].to_numpy(dtype=np.float64)
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    ok &= (vh + vl) > 1e-8
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    tstat = np.full(expr.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    s = pd.Series(tstat, index=expr.index).replace([np.inf, -np.inf], np.nan).dropna()
    s = s.drop(labels=["TACSTD2"], errors="ignore")
    return s.sort_values(ascending=False)


def spearman_rank(expr: pd.DataFrame, tac: pd.Series) -> pd.Series:
    samples = [s for s in tac.index if s in expr.columns and np.isfinite(tac.loc[s])]
    X = expr.loc[:, samples].to_numpy(dtype=np.float64)
    finite_gene = np.isfinite(X).all(axis=1)
    X = X[finite_gene]
    genes = expr.index.to_numpy()[finite_gene]
    # drop TACSTD2
    keep = genes != "TACSTD2"
    X = X[keep]
    genes = genes[keep]
    xr = stats.rankdata(X, method="average", axis=1)
    tr = stats.rankdata(tac.loc[samples].to_numpy(dtype=float), method="average")
    xr_c = xr - xr.mean(axis=1, keepdims=True)
    tr_c = tr - tr.mean()
    num = xr_c @ tr_c
    den = np.sqrt((xr_c ** 2).sum(axis=1) * float((tr_c ** 2).sum()))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = num / den
    rho = np.where(np.isfinite(rho), rho, np.nan)
    s = pd.Series(rho, index=genes).dropna().sort_values(ascending=False)
    return s


def enrichment_walk(abs_s: np.ndarray, hit: np.ndarray):
    n_hit = int(hit.sum())
    n_miss = int((~hit).sum())
    if n_hit == 0 or n_miss == 0:
        return 0.0, 0
    hit_w = np.where(hit, abs_s, 0.0)
    denom = hit_w.sum()
    if denom == 0:
        hit_w = hit.astype(float) / n_hit
    else:
        hit_w = hit_w / denom
    step_miss = (~hit).astype(float) / n_miss
    walk = np.cumsum(hit_w - step_miss)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), i


def null_es(abs_s: np.ndarray, n_hit: int, rng: np.random.Generator) -> np.ndarray:
    n = abs_s.size
    n_miss = n - n_hit
    out = np.empty(NPERM, dtype=float)
    ar = np.arange(n_hit)
    for i in range(NPERM):
        pos = np.sort(rng.choice(n, size=n_hit, replace=False))
        weights = abs_s[pos]
        denom = weights.sum()
        if denom == 0:
            weights = np.full(n_hit, 1.0 / n_hit)
        else:
            weights = weights / denom
        phit = np.cumsum(weights)
        pmiss = (pos - ar) / n_miss
        walk = phit - pmiss
        walk_before = (phit - weights) - pmiss
        cand = np.concatenate([walk, walk_before])
        out[i] = float(cand[int(np.argmax(np.abs(cand)))])
    return out


def gsea_table(rank: pd.Series, gene_sets: dict[str, list[str]], rng: np.random.Generator) -> pd.DataFrame:
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    abs_s = np.abs(scores)
    n = len(genes)
    gene_pos = {g: i for i, g in enumerate(genes)}
    prepared = []
    for term in sorted(gene_sets):
        members = gene_sets[term]
        idx = [gene_pos[g] for g in members if g in gene_pos]
        idx = np.unique(np.asarray(idx, dtype=int))
        if len(idx) < MIN_SET or len(idx) > MAX_SET:
            continue
        hit = np.zeros(n, dtype=bool)
        hit[idx] = True
        es, lead_i = enrichment_walk(abs_s, hit)
        prepared.append((term, idx, hit, es, lead_i))
    null_by_size = {}
    rows = []
    for term, idx, hit, es, lead_i in prepared:
        n_hit = int(len(idx))
        if n_hit not in null_by_size:
            null_by_size[n_hit] = null_es(abs_s, n_hit, rng)
        null = null_by_size[n_hit]
        if es >= 0:
            pos = null[null >= 0]
            nes = float(es / pos.mean()) if len(pos) and pos.mean() != 0 else 0.0
            nom_p = float((np.sum(null >= es) + 1) / (NPERM + 1))
            lead = genes[: lead_i + 1][hit[: lead_i + 1]]
        else:
            neg = null[null < 0]
            nes = float(es / abs(neg.mean())) if len(neg) and neg.mean() != 0 else 0.0
            nom_p = float((np.sum(null <= es) + 1) / (NPERM + 1))
            lead = genes[lead_i:][hit[lead_i:]]
        rows.append({
            "term": term,
            "es": es,
            "nes": nes,
            "nom_p": nom_p,
            "n_set_in_rank": n_hit,
            "lead_genes": ",".join(list(lead[:20])),
        })
    return pd.DataFrame(rows)


def gsea_for_cohort(cohort: dict, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    expr = cohort["expr"]
    tac = cohort["tac"]
    frames = []
    contrasts = []
    high, low = split_extremes(tac, 0.25, 0.75)
    if len(high) >= MIN_ARM_Q and len(low) >= MIN_ARM_Q:
        contrasts.append(("q4_vs_q1_welch", welch_rank(expr, high, low), len(high), len(low)))
    high, low = split_extremes(tac, 1 / 3, 2 / 3)
    if len(high) >= MIN_ARM_Q and len(low) >= MIN_ARM_Q:
        contrasts.append(("tercile_welch", welch_rank(expr, high, low), len(high), len(low)))
    high, low = split_extremes(tac, 0.10, 0.90)
    if len(high) >= MIN_ARM_DECILE and len(low) >= MIN_ARM_DECILE:
        contrasts.append(("decile_welch", welch_rank(expr, high, low), len(high), len(low)))
    contrasts.append(("spearman_vs_tacstd2", spearman_rank(expr, tac), cohort["n"], cohort["n"]))
    for contrast, rank, n_high, n_low in contrasts:
        print(f"  GSEA {cohort['name']} {contrast} genes={len(rank)}", flush=True)
        rng = np.random.default_rng(SEED)
        tab = gsea_table(rank, gene_sets, rng)
        if tab.empty:
            continue
        tab.insert(0, "cohort", cohort["name"])
        tab.insert(1, "contrast", contrast)
        tab["n_high"] = n_high
        tab["n_low"] = n_low
        tab["strict_luad"] = cohort["strict_luad"]
        tab["n_cohort"] = cohort["n"]
        frames.append(tab)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Cohort build
# ---------------------------------------------------------------------------


def finish_cohort(name, expr, tac_ok, strict, platform, composition, composition_name,
                  tumor_axis, tumor_axis_high, tumor_axis_name, do_estimate, est_genes,
                  cosine=None, estimate_score=None) -> dict:
    expr = expr.loc[:, tac_ok]
    if "TACSTD2" not in expr.index:
        raise SystemExit(f"{name} missing TACSTD2")
    tac = expr.loc["TACSTD2"].astype(float)
    endpoints, extra = prepare_endpoints(expr, do_estimate, est_genes)
    if do_estimate and "estimate_stromal" in extra:
        if composition is None:
            composition = extra["estimate_stromal"]
        if tumor_axis is None:
            tumor_axis = extra["estimate_stromal"]
            tumor_axis_high = False
            tumor_axis_name = "ESTIMATE_StromalScore"
        if estimate_score is None:
            estimate_score = extra["estimate_score"]
    if cosine is None and estimate_score is not None:
        cosine = cosine_from_score(estimate_score)
    print(
        f"{name}: n={expr.shape[1]} genes={expr.shape[0]} endpoints={len(endpoints)} strict={strict}",
        flush=True,
    )
    return {
        "name": name,
        "expr": expr,
        "tac": tac,
        "n": int(expr.shape[1]),
        "strict_luad": strict,
        "platform": platform,
        "endpoints": endpoints,
        "composition": None if composition is None else composition.reindex(expr.columns).astype(float),
        "composition_name": composition_name,
        "tumor_axis": None if tumor_axis is None else tumor_axis.reindex(expr.columns).astype(float),
        "tumor_axis_high": tumor_axis_high,
        "tumor_axis_name": tumor_axis_name,
        "cosine": None if cosine is None else cosine.reindex(expr.columns).astype(float),
        "estimate_score": None if estimate_score is None else estimate_score.reindex(expr.columns).astype(float),
        "estimate_immune_n": extra.get("estimate_immune_n"),
        "estimate_stromal_n": extra.get("estimate_stromal_n"),
    }


def build_cohorts(est_genes) -> list[dict]:
    cohorts = []
    print("loading OncoSG", flush=True)
    expr, purity = load_oncosg()
    cohorts.append(finish_cohort(
        "OncoSG", expr, list(expr.columns), True, "cBioPortal z-score RSEM",
        purity, "published_PURITY", purity, True, "published_PURITY",
        False, est_genes,
    ))

    print("loading GPL annotations", flush=True)
    gpl570 = load_gpl_annot(CACHE / "GPL570.annot.gz")
    gpl96 = load_gpl_annot(CACHE / "GPL96.annot.gz")
    gpl15048 = load_gpl15048(CACHE / "GPL15048_platform.txt")

    def geo_array(name, matrix, probe2gene, mask, strict, platform):
        print(f"loading {name}", flush=True)
        meta, probes = parse_geo_matrix(matrix)
        keep = mask(meta)
        print(f"  {name} tumor samples {int(keep.sum())} / {len(keep)}", flush=True)
        gene = collapse_maxmean(probes.loc[:, meta.index[keep]], probe2gene)
        return finish_cohort(
            name, gene, list(gene.columns), strict, platform,
            None, "ESTIMATE_StromalScore", None, False, "ESTIMATE_StromalScore",
            True, est_genes,
        )

    def mask_31210(meta):
        col = "tissue" if "tissue" in meta.columns else "source"
        return meta[col].astype(str).str.lower().str.contains("primary lung tumor")

    def mask_72094(meta):
        s = meta["source"].astype(str).str.lower() if "source" in meta.columns else meta.iloc[:, 0].astype(str)
        return s.str.contains("adenocarcinoma")

    def mask_68465(meta):
        return meta["disease_state"].astype(str).str.contains("Adenocarcinoma", case=False, na=False)

    def mask_19804(meta):
        t = meta["tissue"].astype(str).str.lower()
        return t.str.contains("lung cancer") & ~t.str.contains("adjacent") & ~t.str.contains("normal")

    cohorts.append(geo_array(
        "GSE31210", CACHE / "GSE31210_series_matrix.txt.gz", gpl570, mask_31210,
        True, "GPL570 MAS5",
    ))
    cohorts.append(geo_array(
        "GSE72094", CACHE / "GSE72094_series_matrix.txt.gz", gpl15048, mask_72094,
        True, "GPL15048 HuRSTA",
    ))
    cohorts.append(geo_array(
        "GSE68465", CACHE / "GSE68465_series_matrix.txt.gz", gpl96, mask_68465,
        True, "GPL96 U133A",
    ))
    cohorts.append(geo_array(
        "GSE19804", CACHE / "GSE19804_series_matrix.txt.gz", gpl570, mask_19804,
        False, "GPL570 MAS5",
    ))

    print("loading RNA-seq GEO", flush=True)
    id_to_sym = load_probemap(RNA / "gencode.v36.probemap")
    g273 = load_gse273377(id_to_sym)
    for stratum, expr in g273.items():
        strom, used = zmean(expr, STROMAL9)
        if len(used) < 6:
            raise SystemExit(f"GSE273377 {stratum} stromal genes {len(used)}")
        cohorts.append(finish_cohort(
            f"GSE273377_{stratum}", expr, list(expr.columns), True, "RNA-seq log2 CPM",
            strom, "stromal9_zmean", strom, False, "stromal9_zmean",
            True, est_genes,
        ))
    e233 = load_gse233774(id_to_sym)
    strom, used = zmean(e233, STROMAL9)
    cohorts.append(finish_cohort(
        "GSE233774_tumor", e233, list(e233.columns), True, "RNA-seq log2 FPKM",
        strom, "stromal9_zmean", strom, False, "stromal9_zmean",
        True, est_genes,
    ))
    e282 = load_gse282774()
    strom, used = zmean(e282, STROMAL9)
    cohorts.append(finish_cohort(
        "GSE282774", e282, list(e282.columns), True, "RNA-seq log2 FPKM",
        strom, "stromal9_zmean", strom, False, "stromal9_zmean",
        True, est_genes,
    ))
    return cohorts


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def check_anchors(rho: pd.DataFrame):
    fails = []
    for (cohort, endpoint, subset, cov), expect in ANCHOR.items():
        sub = rho[
            (rho.cohort == cohort) & (rho.endpoint == endpoint)
            & (rho.subset == subset) & (rho.covariate_class == cov)
        ]
        if sub.empty:
            fails.append(f"missing anchor {cohort} {endpoint} {subset} {cov}")
            continue
        got = float(sub.iloc[0]["rho"])
        delta = abs(got - expect)
        print(f"anchor {cohort} {endpoint} {cov}: got {got:.4f} expected {expect:.3f} |d|={delta:.4f}", flush=True)
        if delta > 0.02:
            fails.append(f"{cohort} {endpoint} {cov} got {got:.4f} expected {expect:.3f}")
    if fails:
        raise SystemExit("anchor mismatch:\n" + "\n".join(fails))


def reference_warnings(rho: pd.DataFrame):
    refs = {
        ("GSE72094", "CD8A", "all", "none"): -0.227,
        ("GSE68465", "CD8A", "all", "none"): -0.152,
        ("GSE19804", "CD8A", "all", "none"): -0.024,
        ("GSE273377_validation", "CD8A", "all", "none"): -0.423,
        ("GSE282774", "CD8A", "all", "none"): -0.206,
        ("GSE233774_tumor", "CD8A", "all", "none"): -0.437,
    }
    # GSE273377 discovery CD8A was +0.039 in PR 719 (immune8-era single gene).
    rows = []
    for key, expect in refs.items():
        cohort, endpoint, subset, cov = key
        sub = rho[
            (rho.cohort == cohort) & (rho.endpoint == endpoint)
            & (rho.subset == subset) & (rho.covariate_class == cov)
        ]
        if sub.empty:
            print(f"ref missing {key}", flush=True)
            continue
        got = float(sub.iloc[0]["rho"])
        rows.append({"cohort": cohort, "endpoint": endpoint, "published_rho": expect, "this_run_rho": got, "abs_delta": abs(got - expect)})
        print(f"ref {cohort} {endpoint}: got {got:.4f} prior {expect:.3f}", flush=True)
    return pd.DataFrame(rows)


def meta_table(rho: pd.DataFrame) -> pd.DataFrame:
    elig = rho[rho.eligible_max].copy()
    rows = []
    keys = ["endpoint", "endpoint_family", "subset", "covariate_class"]
    for key, sub in elig.groupby(keys, dropna=False):
        k_cov = 0 if key[3] == "none" else 1
        meta = dl_meta(sub.rho.to_numpy(), sub.n.to_numpy(), k_cov)
        if meta is None or meta["k"] < 3:
            continue
        rec = dict(zip(keys, key))
        rec.update(meta)
        rec["cohorts"] = ",".join(sub.sort_values("cohort").cohort.tolist())
        rec["min_rho"] = float(sub.rho.min())
        rec["max_rho"] = float(sub.rho.max())
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["abs_rho"] = out["rho"].abs()
    return out.sort_values("abs_rho", ascending=False)


def consistent_nes(nes: pd.DataFrame, min_k: int = 5) -> pd.DataFrame:
    """Term x contrast positive in every strict LUAD cohort with n >= 50."""
    strict = nes[(nes.strict_luad.astype(bool)) & (nes.n_cohort >= MIN_N)].copy()
    rows = []
    for (term, contrast), g in strict.groupby(["term", "contrast"], dropna=False):
        if len(g) < min_k or (g.nes <= 0).any():
            continue
        best = g.loc[g.nes.idxmax()]
        worst = g.loc[g.nes.idxmin()]
        rows.append({
            "term": term,
            "contrast": contrast,
            "k": int(len(g)),
            "min_nes": float(g.nes.min()),
            "median_nes": float(g.nes.median()),
            "max_nes": float(g.nes.max()),
            "best_cohort": best.cohort,
            "best_nes": float(best.nes),
            "best_p": float(best.nom_p),
            "best_q": float(best.q_primary) if np.isfinite(best.q_primary) else np.nan,
            "worst_cohort": worst.cohort,
            "worst_nes": float(worst.nes),
            "cohorts": ",".join(sorted(g.cohort)),
        })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["max_nes", "min_nes"], ascending=False)


def write_finding(rho, meta, nes, summary):
    cell = summary["cell_max"]
    meta_w = summary["meta_max"]
    nes_w = summary["nes_max"]
    nes_locked = summary["nes_locked_max"]
    lines = []
    lines.append("# OncoSG + GEO LUAD: maximum story-fit TACSTD2 |ρ| and TJ/junction NES")
    lines.append("")
    lines.append("Public bulk LUAD only. Every ρ and NES below is written by `analyze.py` from the downloaded matrices. The maximum is the extreme of a pre-declared grid, so its p-value is smaller than a single pre-specified test. Sweep q is Benjamini–Hochberg across that grid.")
    lines.append("")
    lines.append("Story-fit rule, fixed before the maximum was read: keep TACSTD2–CD8/ImmuneScore correlations only when ρ is negative, and keep TJ/junction enrichments only when NES is positive in TACSTD2-high. GSE19804 (tissue labeled lung cancer) and GSE233774 (n = 30) are computed and excluded from the headline. GSE10072, GSE11969, and GSE248378 stay closed.")
    lines.append("")
    lines.append("## Maximum |ρ|")
    lines.append("")
    lines.append(
        f"Eligible cells: {summary['n_eligible_rho']} of {summary['n_primary_rho']} primary-family tests "
        f"({summary['n_rho_rows']} rows stored, including companions and ineligible cohorts). "
        f"Eligible |ρ| median {summary['elig_abs_median']:.3f}, 95th percentile {summary['elig_abs_p95']:.3f}."
    )
    lines.append("")
    lines.append(
        f"**Search maximum:** {cell['cohort']}, {cell['endpoint']} vs TACSTD2, subset `{cell['subset']}`, "
        f"covariate `{cell['covariate_name']}` ({cell['covariate_class']}). "
        f"n = {int(cell['n'])}. Spearman ρ = **{cell['rho']:.3f}** (p = {fmt_p(cell['p'])}, "
        f"primary-family q = {fmt_p(cell['q_primary'])}). "
        f"Genes in the endpoint: {cell['n_genes']} ({cell['genes_used']})."
    )
    lines.append("")
    lines.append(summary["cell_context"])
    lines.append("")
    if meta_w:
        lines.append(
            f"**Cross-cohort maximum:** {meta_w['endpoint']}, subset `{meta_w['subset']}`, "
            f"covariate class `{meta_w['covariate_class']}`. "
            f"Pooled ρ = **{meta_w['rho']:.3f}** ({meta_w['model']}, 95% CI {meta_w['ci_low']:.3f} to {meta_w['ci_high']:.3f}, "
            f"p = {fmt_p(meta_w['p'])}, I² = {meta_w['i2']:.0%}, k = {int(meta_w['k'])}, n sum = {int(meta_w['n_sum'])}). "
            f"Cohorts: {meta_w['cohorts']}."
        )
        lines.append("")
    lines.append("Locked single-gene reference, recomputed in this run. The partial column is the composition covariate: published PURITY on OncoSG, ESTIMATE StromalScore on arrays, and the 9-gene stromal z-mean on RNA-seq GEO.")
    lines.append("")
    lines.append("| cohort | n | CD8A unadjusted ρ | CD8A composition partial ρ | immune8 unadjusted ρ | immune8 composition partial ρ |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for cohort in ["OncoSG", "GSE31210", "GSE72094", "GSE68465", "GSE273377_discovery", "GSE273377_validation", "GSE282774", "GSE19804", "GSE233774_tumor"]:
        def grab(ep, subset, cov):
            sub = rho[(rho.cohort == cohort) & (rho.endpoint == ep) & (rho.subset == subset) & (rho.covariate_class == cov)]
            if sub.empty:
                return "NA", "NA"
            r = sub.iloc[0]
            return str(int(r.n)), fmt_rho(r.rho)
        n1, a = grab("CD8A", "all", "none")
        _, b = grab("CD8A", "all", "composition")
        # OncoSG immune anchor is the deposited mean; others use z-mean
        ep_i = "immune8_deposited_mean" if cohort == "OncoSG" else "immune8_zmean"
        _, c = grab(ep_i, "all", "none")
        _, d = grab(ep_i, "all", "composition")
        lines.append(f"| {cohort} | {n1} | {a} | {b} | {c} | {d} |")
    lines.append("")
    est_bits = []
    for cohort in ["GSE31210", "GSE72094", "GSE68465"]:
        hit = rho[
            (rho.cohort == cohort) & (rho.endpoint == "CD8A")
            & (rho.subset == "all") & (rho.covariate_class == "estimate_score")
        ]
        if len(hit):
            est_bits.append(f"{cohort} {float(hit.iloc[0].rho):.3f}")
    lines.append(
        "Array partials on ESTIMATEScore (stromal ssGSEA + immune ssGSEA), the GSE31210 covariate in PR 736: "
        + "; ".join(est_bits)
        + ". Cosine TumorPurity matches those three partials exactly in this run."
    )
    lines.append("")
    lines.append(
        f"CD8A, all samples, unadjusted, pooled across cohorts where that correlation is negative: "
        f"ρ = {summary['cd8a_meta_rho']:.3f} (k = {summary['cd8a_meta_k']}, I² = {summary['cd8a_meta_i2']:.0%}). "
        f"immune8 z-mean on the same slice: ρ = {summary['imm8_meta_rho']:.3f} (k = {summary['imm8_meta_k']}, I² = {summary['imm8_meta_i2']:.0%})."
    )
    lines.append("")
    lines.append(summary["discovery_sentence"])
    lines.append("")
    lines.append(summary["gep_sentence"])
    lines.append("")
    lines.append("## Maximum TJ/junction NES in TACSTD2-high")
    lines.append("")
    lines.append(
        f"Library: {summary['n_sets']} name-filtered sets (`data/junction_sets.gmt`). "
        f"Tests in the strict-LUAD n≥50 family: {summary['n_nes_primary']}. "
        f"Positive NES (story-fit): {summary['n_nes_positive']}."
    )
    lines.append("")
    if nes_w:
        lines.append(
            f"**Positive-NES maximum:** {nes_w['cohort']}, `{nes_w['term']}`, contrast `{nes_w['contrast']}`. "
            f"NES = **{nes_w['nes']:.3f}** (ES = {nes_w['es']:.3f}, nominal p = {fmt_p(nes_w['nom_p'])}, "
            f"family q = {fmt_p(nes_w['q_primary'])}, genes in rank = {int(nes_w['n_set_in_rank'])}, "
            f"arms {int(nes_w['n_high'])}/{int(nes_w['n_low'])}). "
            f"Leading edge: {nes_w['lead_genes']}."
        )
        lines.append("")
        lines.append(summary["nes_extreme_context"])
        lines.append("")
    if nes_locked:
        lines.append(
            f"**Locked A8 set with the largest positive NES** (KEGG tight junction, GO tight-junction organization, "
            f"GO bicellular tight-junction assembly, Hallmark apical junction): "
            f"{nes_locked['cohort']}, `{nes_locked['term']}`, `{nes_locked['contrast']}`, "
            f"NES = **{nes_locked['nes']:.3f}** (nominal p = {fmt_p(nes_locked['nom_p'])}, "
            f"q within the four-set family = {fmt_p(nes_locked['q_locked'])})."
        )
        lines.append("")
    cons = summary.get("consistent_top")
    if cons:
        lines.append(
            f"**Positive in all {int(cons['k'])} strict LUAD cohorts:** `{cons['term']}`, `{cons['contrast']}`, "
            f"largest NES = **{cons['best_nes']:.3f}** in {cons['best_cohort']} "
            f"(nominal p = {fmt_p(cons['best_p'])}, family q = {fmt_p(cons['best_q'])}). "
            f"The smallest NES in that specification is {cons['min_nes']:.3f} in {cons['worst_cohort']}. "
            f"The specification with the highest floor is `{summary['consistent_floor']['term']}`, "
            f"`{summary['consistent_floor']['contrast']}`, minimum NES **{summary['consistent_floor']['min_nes']:.3f}** "
            f"({summary['consistent_floor']['worst_cohort']}) and maximum {summary['consistent_floor']['max_nes']:.3f} "
            f"({summary['consistent_floor']['best_cohort']})."
        )
        lines.append("")
    lines.append("A8 sets at the Q4 vs Q1 Welch contrast (positive NES is the story direction):")
    lines.append("")
    lines.append("| cohort | KEGG tight junction | GO TJ organization | GO bicellular TJ assembly | Hallmark apical junction |")
    lines.append("|---|---:|---:|---:|---:|")
    a8_terms = [
        "A8::KEGG_TIGHT_JUNCTION",
        "A8::GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "A8::GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
        "A8::HALLMARK_APICAL_JUNCTION",
    ]
    order = ["OncoSG", "GSE31210", "GSE72094", "GSE68465", "GSE273377_discovery", "GSE273377_validation", "GSE282774"]
    sub = nes[(nes.contrast == "q4_vs_q1_welch") & (nes.term.isin(a8_terms))]
    for cohort in order:
        vals = []
        for term in a8_terms:
            hit = sub[(sub.cohort == cohort) & (sub.term == term)]
            vals.append(fmt_rho(float(hit.iloc[0].nes)) if len(hit) else "NA")
        lines.append(f"| {cohort} | " + " | ".join(vals) + " |")
    lines.append("")
    lines.append(summary["nes_direction_sentence"])
    lines.append("")
    lines.append("## What this is")
    lines.append("")
    lines.append("These are bulk RNA associations in surgical or resected LUAD. They are not a spatial exclusion measurement and not an ICI-response result. OncoSG remains the East-Asian z-score matrix with n = 169 (portal RNA list 181 is not the matrix n). ESTIMATE is not applied to OncoSG z-scores. Selecting the grid maximum inflates |effect| relative to the locked CD8A and KEGG rows above.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r methods/luad_tacstd2_max_effect/requirements.txt")
    lines.append("python3 methods/luad_tacstd2_max_effect/analyze.py")
    lines.append("```")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (HERE / "FINDING.md").write_text(text)
    return text


def plots(rho, meta, nes, summary):
    FIGURES.mkdir(parents=True, exist_ok=True)
    elig = rho[rho.eligible_max].copy()
    elig["abs_rho"] = elig.rho.abs()
    # fig1 distribution
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.hist(elig.abs_rho, bins=30, color="#1f4e79", edgecolor="white")
    cell = summary["cell_max"]
    ax.axvline(abs(cell["rho"]), color="#b03a2e", lw=1.5, label=f"maximum |ρ| {abs(cell['rho']):.3f}")
    ax.axvline(summary["elig_abs_median"], color="#7f8c8d", lw=1, ls="--", label=f"median {summary['elig_abs_median']:.3f}")
    ax.set_xlabel("|Spearman ρ| among story-fit cells")
    ax.set_ylabel("cells")
    ax.set_title("TACSTD2 vs CD8 / ImmuneScore, negative ρ only")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_abs_rho_distribution.png", dpi=140)
    plt.close(fig)

    # fig2 top cells
    top = elig.sort_values("abs_rho", ascending=False).head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    labels = [f"{r.cohort} | {r.endpoint} | {r.subset} | {r.covariate_class}" for r in top.itertuples()]
    ax.barh(np.arange(len(top)), top.rho, color="#1f4e79")
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Spearman ρ")
    ax.set_title("Largest story-fit |ρ| cells")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_top_rho_cells.png", dpi=140)
    plt.close(fig)

    # fig3 meta winner forest if present
    if summary["meta_max"]:
        mw = summary["meta_max"]
        sub = rho[
            (rho.endpoint == mw["endpoint"]) & (rho.subset == mw["subset"])
            & (rho.covariate_class == mw["covariate_class"])
            & (rho.strict_luad) & (rho.n >= MIN_N)
        ].sort_values("cohort")
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        y = np.arange(len(sub))
        colors = ["#1f4e79" if r < 0 else "#7f8c8d" for r in sub.rho]
        ax.barh(y, sub.rho, color=colors)
        ax.axvline(mw["rho"], color="#b03a2e", lw=1, label=f"pooled ρ {mw['rho']:.3f}")
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{c} (n={int(n)})" for c, n in zip(sub.cohort, sub.n)], fontsize=8)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(f"{mw['endpoint']} | {mw['subset']} | {mw['covariate_class']}")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / "fig3_meta_winner_forest.png", dpi=140)
        plt.close(fig)

    # fig4 A8 NES heatmap-like bars
    a8_terms = [
        "A8::KEGG_TIGHT_JUNCTION",
        "A8::GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "A8::GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
        "A8::HALLMARK_APICAL_JUNCTION",
    ]
    short = {
        "A8::KEGG_TIGHT_JUNCTION": "KEGG TJ",
        "A8::GOBP_TIGHT_JUNCTION_ORGANIZATION": "GO TJ org.",
        "A8::GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY": "GO bicell. TJ",
        "A8::HALLMARK_APICAL_JUNCTION": "Hallmark apical",
    }
    order = ["OncoSG", "GSE31210", "GSE72094", "GSE68465", "GSE273377_discovery", "GSE273377_validation", "GSE282774"]
    sub = nes[(nes.contrast == "q4_vs_q1_welch") & (nes.term.isin(a8_terms)) & (nes.cohort.isin(order))]
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    x = np.arange(len(order))
    width = 0.18
    for i, term in enumerate(a8_terms):
        vals = []
        for c in order:
            hit = sub[(sub.cohort == c) & (sub.term == term)]
            vals.append(float(hit.iloc[0].nes) if len(hit) else np.nan)
        ax.bar(x + (i - 1.5) * width, vals, width=width, label=short[term])
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("NES (TACSTD2 Q4 vs Q1)")
    ax.set_title("Locked TJ / junction sets")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_a8_nes_q4q1.png", dpi=140)
    plt.close(fig)

    # fig5: specification positive in every strict cohort, highest floor
    cons = consistent_nes(nes, min_k=5)
    if len(cons):
        floor = cons.sort_values(["min_nes", "max_nes"], ascending=False).iloc[0]
        sub = nes[(nes.term == floor.term) & (nes.contrast == floor.contrast) & (nes.strict_luad.astype(bool)) & (nes.n_cohort >= MIN_N)]
        sub = sub.sort_values("nes")
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
        colors = ["#1f4e79" if v > 0 else "#7f8c8d" for v in sub.nes]
        ax.barh(np.arange(len(sub)), sub.nes, color=colors)
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels(sub.cohort, fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("NES")
        label = floor.term.split("::", 1)[-1]
        ax.set_title(f"{label}\n{floor.contrast} (positive in all {int(floor.k)} strict cohorts)")
        fig.tight_layout()
        fig.savefig(FIGURES / "fig5_consistent_nes.png", dpi=140)
        plt.close(fig)


def _rho_cell(rho, cohort, endpoint, subset, cov):
    sub = rho[
        (rho.cohort == cohort) & (rho.endpoint == endpoint)
        & (rho.subset == subset) & (rho.covariate_class == cov)
    ]
    if sub.empty:
        return None
    return sub.iloc[0]


def cell_context(rho) -> str:
    base = _rho_cell(rho, "GSE273377_validation", "CD8A", "all", "none")
    comp = _rho_cell(rho, "GSE273377_validation", "CD8A", "all", "composition")
    est = _rho_cell(rho, "GSE273377_validation", "CD8A", "all", "estimate_score")
    cos = _rho_cell(rho, "GSE273377_validation", "CD8A", "all", "cosine_purity")
    onco = _rho_cell(rho, "OncoSG", "cytotoxic6_zmean", "high_tumor", "none")
    same = (
        abs(float(est.rho) - float(cos.rho)) < 1e-12
        if est is not None and cos is not None else False
    )
    same_txt = (
        "Cosine TumorPurity and ESTIMATEScore give the same partial ρ in this cohort, "
        "so the cosine transform is monotone with ESTIMATEScore here and is the same adjustment."
        if same else
        "Cosine TumorPurity and ESTIMATEScore partials differ in this cohort."
    )
    return (
        f"In that same cohort the unadjusted CD8A ρ is {base.rho:.3f} and the stromal9 partial is {comp.rho:.3f}. "
        f"The search maximum is {abs(float(est.rho) - float(base.rho)):.3f} above the unadjusted CD8A ρ. "
        f"{same_txt} "
        f"OncoSG high-tumor cytotoxic6 (CD8A, GZMA, GZMB, PRF1, NKG7, GNLY z-mean), no further covariate, "
        f"is ρ = {onco.rho:.3f} (n = {int(onco.n)}, q = {fmt_p(onco.q_primary)})."
    )


def discovery_sentence(rho) -> str:
    bits = []
    for ep in ["CD8A", "cytotoxic6_zmean", "immune8_zmean", "estimate_immune141"]:
        row = _rho_cell(rho, "GSE273377_discovery", ep, "all", "none")
        if row is not None:
            bits.append(f"{ep} {row.rho:+.3f}")
    return (
        "GSE273377 discovery (n = 103) is strict LUAD and is left out of the negative pools where its sign is positive. "
        "All-sample unadjusted ρ: " + "; ".join(bits) + ". "
        "The pooled ρ values are therefore pooled over the cohorts that stay negative."
    )


def gep_sentence(rho) -> str:
    g = rho[(rho.endpoint == "ayers_gep18_zmean") & (rho.strict_luad) & (rho.n >= MIN_N) & (rho.rho < 0)]
    if g.empty:
        return "Ayers GEP18 has no story-fit cell."
    top = g.loc[g.rho.abs().idxmax()]
    return (
        f"Ayers GEP18 is a companion, outside the CD8/ImmuneScore maximum. "
        f"Its largest negative ρ on strict LUAD with n ≥ 50 is {top.rho:.3f} "
        f"({top.cohort}, {top.subset}, {top.covariate_class}, n = {int(top.n)})."
    )


def nes_extreme_context(nes) -> str:
    term = "A8::HALLMARK_APICAL_JUNCTION"
    contrast = "tercile_welch"
    sub = nes[(nes.term == term) & (nes.contrast == contrast) & (nes.strict_luad.astype(bool))]
    bits = [f"{r.cohort} {r.nes:+.2f}" for r in sub.sort_values("cohort").itertuples()]
    kegg = nes[(nes.term == "A8::KEGG_TIGHT_JUNCTION") & (nes.contrast == "q4_vs_q1_welch")]
    kegg_bits = []
    for cohort in ["OncoSG", "GSE31210"]:
        hit = kegg[kegg.cohort == cohort]
        if len(hit):
            kegg_bits.append(f"{cohort} {float(hit.iloc[0].nes):+.2f}")
    return (
        "That Hallmark apical-junction specification across strict LUAD: " + "; ".join(bits) + ". "
        "OncoSG and GSE31210 are negative, so this cell is the positive-NES extreme and does not carry a cross-cohort apical-junction claim. "
        "Locked KEGG tight junction, Q4 vs Q1 Welch, recomputed here: " + "; ".join(kegg_bits) +
        " (PR 736 reported OncoSG −1.00 and GSE31210 +2.00 under the same ranking statistic and an independent permutation stream)."
    )


def direction_sentence(nes: pd.DataFrame) -> str:
    a8 = [
        "A8::KEGG_TIGHT_JUNCTION",
        "A8::GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "A8::GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
        "A8::HALLMARK_APICAL_JUNCTION",
    ]
    sub = nes[(nes.contrast == "q4_vs_q1_welch") & (nes.term.isin(a8)) & (nes.cohort == "OncoSG")]
    if sub.empty:
        return "OncoSG Q4 vs Q1 A8 rows were not produced."
    pos = sub[sub.nes > 0]
    neg = sub[sub.nes <= 0]
    bits = [f"{r.term.split('::',1)[1]} NES {r.nes:.2f}" for r in sub.itertuples()]
    return (
        "OncoSG Q4 vs Q1 for the four locked sets: "
        + "; ".join(bits)
        + f". {len(pos)} of 4 are positive in TACSTD2-high. "
        + (f"Sets at or below zero: {', '.join(neg.term.str.split('::').str[-1])}." if len(neg) else "All four are positive.")
    )


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    ensure_inputs()
    if GMT_PATH.exists():
        gene_sets = load_gmt()
        print(f"loaded GMT {len(gene_sets)}", flush=True)
    else:
        gene_sets = build_gmt()
    est = pd.read_csv(DATA / "estimate_yoshihara_2013.tsv", sep="\t")
    est_genes = {
        "stromal": est.loc[est["set"] == "Stromal141_UP", "hugo"].astype(str).tolist(),
        "immune": est.loc[est["set"] == "Immune141_UP", "hugo"].astype(str).tolist(),
    }
    cohorts = build_cohorts(est_genes)
    rho_rows = []
    nes_frames = []
    coverage = []
    for c in cohorts:
        rho_rows.extend(rho_rows_for_cohort(c))
        coverage.append({
            "cohort": c["name"],
            "n": c["n"],
            "strict_luad": c["strict_luad"],
            "platform": c["platform"],
            "composition": c["composition_name"],
            "tumor_axis": c["tumor_axis_name"],
            "estimate_immune_n": c["estimate_immune_n"],
            "estimate_stromal_n": c["estimate_stromal_n"],
            "endpoints": ",".join(ep["endpoint"] for ep in c["endpoints"]),
        })
        nes_frames.append(gsea_for_cohort(c, gene_sets))
        # free the matrix
        c["expr"] = None
    rho = pd.DataFrame(rho_rows)
    check_anchors(rho)
    refs = reference_warnings(rho)
    primary = rho["primary_family"].to_numpy()
    q = np.full(len(rho), np.nan)
    q[primary] = bh(rho.loc[primary, "p"].to_numpy())
    rho["q_primary"] = q
    rho.to_csv(TABLES / "rho_grid.tsv", sep="\t", index=False)
    refs.to_csv(TABLES / "published_reference_delta.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(TABLES / "cohort_coverage.tsv", sep="\t", index=False)

    meta = meta_table(rho)
    if len(meta):
        # BH across meta specs is descriptive; the cell-level q is the sweep q
        meta.to_csv(TABLES / "meta_grid.tsv", sep="\t", index=False)
    nes = pd.concat(nes_frames, ignore_index=True) if nes_frames else pd.DataFrame()
    nes_primary = (
        nes["strict_luad"].astype(bool) & (nes["n_cohort"] >= MIN_N) & np.isfinite(nes["nes"])
    )
    nes["primary_family"] = nes_primary
    nes["story_positive"] = nes["nes"] > 0
    nes["eligible_max"] = nes_primary & (nes["nes"] > 0)
    qn = np.full(len(nes), np.nan)
    if nes_primary.any():
        qn[nes_primary.to_numpy()] = bh(nes.loc[nes_primary, "nom_p"].to_numpy())
    nes["q_primary"] = qn
    locked_mask = nes.term.str.startswith("A8::") & nes_primary
    ql = np.full(len(nes), np.nan)
    if locked_mask.any():
        ql[locked_mask.to_numpy()] = bh(nes.loc[locked_mask, "nom_p"].to_numpy())
    nes["q_locked"] = ql
    nes.to_csv(TABLES / "nes_grid.tsv", sep="\t", index=False)

    elig = rho[rho.eligible_max]
    cell = elig.loc[elig.rho.abs().idxmax()].to_dict()
    meta_w = meta.iloc[0].to_dict() if len(meta) else None

    def one_meta(endpoint):
        sub = rho[(rho.endpoint == endpoint) & (rho.subset == "all") & (rho.covariate_class == "none") & (rho.eligible_max)]
        m = dl_meta(sub.rho.to_numpy(), sub.n.to_numpy(), 0)
        if m is None:
            return np.nan, 0, np.nan
        return m["rho"], m["k"], m["i2"]

    cd8_rho, cd8_k, cd8_i2 = one_meta("CD8A")
    # immune8 z-mean; OncoSG deposited mean is the anchor, z-mean is the cross-cohort endpoint
    imm_rho, imm_k, imm_i2 = one_meta("immune8_zmean")

    nes_elig = nes[nes.eligible_max]
    nes_w = nes_elig.loc[nes_elig.nes.idxmax()].to_dict() if len(nes_elig) else None
    locked = nes[nes.eligible_max & nes.term.str.startswith("A8::")]
    nes_locked = locked.loc[locked.nes.idxmax()].to_dict() if len(locked) else None

    abs_vals = elig.rho.abs()
    summary = {
        "n_rho_rows": int(len(rho)),
        "n_primary_rho": int(rho.primary_family.sum()),
        "n_eligible_rho": int(rho.eligible_max.sum()),
        "elig_abs_median": float(abs_vals.median()),
        "elig_abs_p95": float(np.quantile(abs_vals, 0.95)),
        "cell_max": cell,
        "meta_max": meta_w,
        "cd8a_meta_rho": cd8_rho,
        "cd8a_meta_k": cd8_k,
        "cd8a_meta_i2": cd8_i2,
        "imm8_meta_rho": imm_rho,
        "imm8_meta_k": imm_k,
        "imm8_meta_i2": imm_i2,
        "n_sets": len(gene_sets),
        "n_nes_primary": int(nes.primary_family.sum()),
        "n_nes_positive": int(nes.eligible_max.sum()),
        "nes_max": nes_w,
        "nes_locked_max": nes_locked,
        "nes_direction_sentence": direction_sentence(nes),
        "cell_context": cell_context(rho),
        "discovery_sentence": discovery_sentence(rho),
        "gep_sentence": gep_sentence(rho),
        "nes_extreme_context": nes_extreme_context(nes),
        "anchors_checked": {f"{a}|{b}|{c}|{d}": v for (a, b, c, d), v in ANCHOR.items()},
    }
    cons = consistent_nes(nes, min_k=5)
    if len(cons):
        cons.to_csv(TABLES / "nes_consistent.tsv", sep="\t", index=False)
        summary["consistent_top"] = cons.iloc[0].to_dict()
        floor = cons.sort_values(["min_nes", "max_nes"], ascending=False).iloc[0]
        summary["consistent_floor"] = floor.to_dict()
    else:
        summary["consistent_top"] = None
        summary["consistent_floor"] = None
    # json-safe
    def conv(o):
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    (TABLES / "summary.json").write_text(json.dumps(conv(summary), indent=2))
    write_finding(rho, meta, nes, summary)
    plots(rho, meta, nes, summary)
    print("CELL MAX", cell["cohort"], cell["endpoint"], cell["subset"], cell["covariate_class"], cell["rho"], flush=True)
    if meta_w:
        print("META MAX", meta_w["endpoint"], meta_w["subset"], meta_w["covariate_class"], meta_w["rho"], meta_w["k"], flush=True)
    if nes_w:
        print("NES MAX", nes_w["cohort"], nes_w["term"], nes_w["contrast"], nes_w["nes"], flush=True)
    if nes_locked:
        print("NES LOCKED", nes_locked["cohort"], nes_locked["term"], nes_locked["contrast"], nes_locked["nes"], flush=True)
    print("done", flush=True)


def report_only():
    """Rewrite FINDING and figures from tables already produced by main()."""
    rho = pd.read_csv(TABLES / "rho_grid.tsv", sep="\t")
    nes = pd.read_csv(TABLES / "nes_grid.tsv", sep="\t")
    meta = pd.read_csv(TABLES / "meta_grid.tsv", sep="\t") if (TABLES / "meta_grid.tsv").exists() else pd.DataFrame()
    summary = json.loads((TABLES / "summary.json").read_text())
    summary["cell_context"] = cell_context(rho)
    summary["discovery_sentence"] = discovery_sentence(rho)
    summary["gep_sentence"] = gep_sentence(rho)
    summary["nes_extreme_context"] = nes_extreme_context(nes)
    summary["nes_direction_sentence"] = direction_sentence(nes)
    cons = consistent_nes(nes, min_k=5)
    cons.to_csv(TABLES / "nes_consistent.tsv", sep="\t", index=False)
    summary["consistent_top"] = cons.iloc[0].to_dict()
    floor = cons.sort_values(["min_nes", "max_nes"], ascending=False).iloc[0]
    summary["consistent_floor"] = floor.to_dict()
    write_finding(rho, meta, nes, summary)
    plots(rho, meta, nes, summary)
    print("report refreshed", flush=True)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--report-only":
        report_only()
    else:
        main()
