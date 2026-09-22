#!/usr/bin/env python3
"""PAPER FUNNEL — Bulk LUAD (OncoSG + GEO GSE31210): TACSTD2-high vs low
DEG/GSEA for TJ, and TACSTD2/CLDN4 vs ImmuneScore/CD8.

PPT source: User PPT 2026-08-17 (金炫宏), claims A1 and A8.
Rule: never fabricate. Report pass / fail / partial from computed numbers only.

Cohorts (public only)
---------------------
- OncoSG LUAD (cBioPortal luad_oncosg_2020; Chen et al. Nat Genet 2020).
  Public deposit = all-sample z-scores. Honest n = 169 (not portal RNA list 181).
  ImmuneScore = A1 8-gene T-cell effector mean z (not ESTIMATE; raw RSEM absent).
  Purity = published clinical PURITY.
- GSE31210 Japanese stage I–II LUAD tumors (Okayama et al. Cancer Res 2012;
  GPL570). Honest n = 226 primary lung tumors. ImmuneScore = Yoshihara 2013
  ESTIMATE Immune141 ssGSEA. Purity covariate = ESTIMATEScore (rank residual;
  cosine TumorPurity wraps outside [0,1] for most arrays).

Locked design (set before looking at NES / ρ)
---------------------------------------------
- Split (GSEA primary): TACSTD2 top vs bottom quartile.
- Ranking: Welch t (high − low).
- GSEA: preranked weighted KS p=1, 1000 gene-set permutations, seed=42.
  Primary sets = A8 12-set list (keratin/TJ/EMT/apical junction).
- Immune: Spearman of TACSTD2 and CLDN4 vs CD8A and ImmuneScore; partial
  Spearman residualising on the cohort purity covariate.
- PPT call: A1 = TACSTD2 vs CD8/Immune remains negative after purity;
  A8 = TJ/keratin up and Hallmark EMT down in TACSTD2-high.

Outputs -> methods/paper_funnel_luad_tacstd2_tj/{tables,figures,FINDING.md}
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
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
CACHE = Path("/tmp/paper_funnel_luad_tacstd2_tj")
SET_JSON = HERE / "data" / "a8_sets.json"
EST_TSV = HERE / "data" / "estimate_yoshihara_2013.tsv"

UA = (
    "sdaxcge-paper-funnel-luad-tacstd2-tj/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

DATAHUB_SHA = "165bd77077b03038f9c2ee104959eb474770b2a9"
STUDY = "luad_oncosg_2020"
MEDIA = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/"
    f"{DATAHUB_SHA}/public/{STUDY}"
)
ONCOSG_FILES = {
    "expression": (
        "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    ),
    "clinical_sample": (
        "data_clinical_sample.txt",
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
    ),
}

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
    "matrix/GSE31210_series_matrix.txt.gz"
)
GPL570_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
)

IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
FOCAL_TJ = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3", "OCLN", "TJP1", "CGN", "CRB3"]
SEED = 42
NPERM = 1000
MIN_SIZE = 8
MAX_SIZE = 500
EST_A = 0.6049872018
EST_B = 0.0001467884


# ---------------------------------------------------------------------------
# utilities
# ---------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int = 1000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > min_bytes:
        print(f"cached {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"GET {url}\n -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
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
    return f"{r:.3f}"


def fmt_nes(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.2f}"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_pair(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def partial_spearman(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m]).reshape(-1, 1)
    rx = residualize(xr, zr)
    ry = residualize(yr, zr)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def mean_present(expr: pd.DataFrame, genes: list[str], samples) -> tuple[pd.Series, list, list]:
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), present, missing
    return expr.loc[present, samples].astype(float).mean(axis=0), present, missing


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    q = p.copy()
    valid = p.dropna()
    if valid.empty:
        return q
    n = len(valid)
    order = valid.sort_values().index
    adj = valid.loc[order] * n / np.arange(1, n + 1)
    adj = adj[::-1].cummin()[::-1].clip(upper=1.0)
    q.loc[order] = adj
    return q


# ---------------------------------------------------------------------------
# GSEA (same statistic as A8 rework)
# ---------------------------------------------------------------------------


def quartile_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    low = x.index[x <= q1]
    high = x.index[x >= q3]
    return high, low


def rank_high_vs_low(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.Series:
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
    return s.sort_values(ascending=False)


def enrichment_walk(abs_s: np.ndarray, hit: np.ndarray) -> tuple[float, np.ndarray, int]:
    n_hit = int(hit.sum())
    n_miss = int((~hit).sum())
    if n_hit == 0 or n_miss == 0:
        return 0.0, np.zeros(len(hit)), 0
    hit_w = np.where(hit, abs_s, 0.0)
    denom = hit_w.sum()
    if denom == 0:
        hit_w = hit.astype(float) / n_hit
    else:
        hit_w = hit_w / denom
    step_miss = (~hit).astype(float) / n_miss
    walk = np.cumsum(hit_w - step_miss)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), walk, i


def es_from_hits(hit_idx: np.ndarray, abs_s: np.ndarray, n: int) -> float:
    hit_idx = np.sort(np.asarray(hit_idx, dtype=int))
    n_hit = hit_idx.size
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0
    weights = abs_s[hit_idx]
    denom = weights.sum()
    if denom == 0:
        weights = np.full(n_hit, 1.0 / n_hit)
    else:
        weights = weights / denom
    phit = np.cumsum(weights)
    pmiss = (hit_idx - np.arange(n_hit)) / n_miss
    walk = phit - pmiss
    walk_before = (phit - weights) - pmiss
    cand = np.concatenate([walk, walk_before])
    return float(cand[int(np.argmax(np.abs(cand)))])


def gsea_prerank(rank: pd.Series, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    abs_s = np.abs(scores)
    n = len(genes)
    gene_pos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(SEED)
    rows = []
    for term, members in gene_sets.items():
        idx = np.unique(np.array([gene_pos[g] for g in members if g in gene_pos], dtype=int))
        if len(idx) < MIN_SIZE or len(idx) > MAX_SIZE:
            continue
        hit = np.zeros(n, dtype=bool)
        hit[idx] = True
        es, walk, lead_i = enrichment_walk(abs_s, hit)
        if es >= 0:
            lead = genes[: lead_i + 1][hit[: lead_i + 1]]
        else:
            lead = genes[lead_i:][hit[lead_i:]]
        n_hit = int(len(idx))
        null_es = np.empty(NPERM, dtype=float)
        deck = np.arange(n)
        for i in range(NPERM):
            rng.shuffle(deck)
            null_es[i] = es_from_hits(deck[:n_hit], abs_s, n)
        if es >= 0:
            pos = null_es[null_es >= 0]
            nes = float(es / pos.mean()) if len(pos) and pos.mean() != 0 else 0.0
            nom_p = float((np.sum(null_es >= es) + 1) / (NPERM + 1))
        else:
            neg = null_es[null_es < 0]
            nes = float(es / abs(neg.mean())) if len(neg) and neg.mean() != 0 else 0.0
            nom_p = float((np.sum(null_es <= es) + 1) / (NPERM + 1))
        rows.append(
            {
                "term": term,
                "es": es,
                "nes": nes,
                "nom_p": nom_p,
                "n_set_in_rank": n_hit,
                "mean_t": float(scores[hit].mean()),
                "lead_genes": ",".join(lead[:25]),
                "n_lead": int(len(lead)),
            }
        )
    return pd.DataFrame(rows)


def deg_table(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.DataFrame:
    eh = expr.loc[:, high].to_numpy(dtype=np.float64)
    el = expr.loc[:, low].to_numpy(dtype=np.float64)
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    # OncoSG z-scores / Affy log-ish: report mean difference as delta
    delta = mh - ml
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (mh - ml) / se
        # Welch-Satterthwaite df
        num = (vh / nh + vl / nl) ** 2
        den = (vh / nh) ** 2 / (nh - 1) + (vl / nl) ** 2 / (nl - 1)
        df = num / np.where(den > 0, den, np.nan)
        p = 2 * stats.t.sf(np.abs(t), df)
    out = pd.DataFrame(
        {
            "gene": expr.index,
            "mean_high": mh,
            "mean_low": ml,
            "delta_high_minus_low": delta,
            "welch_t": t,
            "welch_p": p,
            "n_high": nh,
            "n_low": nl,
        }
    )
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["welch_t"])
    out["fdr_bh"] = bh_fdr(out["welch_p"])
    return out.sort_values("welch_t", ascending=False)


# ---------------------------------------------------------------------------
# OncoSG loader
# ---------------------------------------------------------------------------


def load_oncosg() -> tuple[pd.DataFrame, pd.Series, pd.Series, dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    paths = {}
    hashes = {}
    for key, (fname, expect_sha) in ONCOSG_FILES.items():
        dest = CACHE / fname
        download(f"{MEDIA}/{fname}", dest, min_bytes=5000)
        got = sha256_file(dest)
        if expect_sha and got != expect_sha:
            raise SystemExit(f"sha256 mismatch {fname}: got {got}")
        paths[key] = dest
        hashes[key] = got

    raw = pd.read_csv(paths["expression"], sep="\t")
    expr = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr = expr.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
    expr = expr.set_index("Hugo_Symbol")
    samples = list(expr.columns)

    clin = pd.read_csv(paths["clinical_sample"], sep="\t", comment="#").set_index(
        "SAMPLE_ID", drop=False
    )
    missing = set(samples) - set(clin.index)
    if missing:
        raise SystemExit(f"OncoSG expression samples missing clinical: {missing}")

    purity = pd.to_numeric(clin.loc[samples, "PURITY"], errors="coerce")
    if int(purity.notna().sum()) != len(samples):
        raise SystemExit("OncoSG PURITY incomplete on matrix samples")

    immune, imm_present, imm_missing = mean_present(expr, IMMUNE8, samples)
    meta = {
        "cohort": "OncoSG_LUAD",
        "n": len(samples),
        "expression_sha256": hashes["expression"],
        "clinical_sha256": hashes["clinical_sample"],
        "datahub_sha": DATAHUB_SHA,
        "matrix": "z-score RSEM (all-sample ref)",
        "immune_score": "A1_8gene_mean_z",
        "immune8_present": imm_present,
        "immune8_missing": imm_missing,
        "purity": "published PURITY",
        "note": "ESTIMATE skipped — public deposit is z-scores only",
    }
    return expr.astype(float), purity, immune, meta


# ---------------------------------------------------------------------------
# GSE31210 loader
# ---------------------------------------------------------------------------


def ssgsea(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 10:
        return pd.Series(np.nan, index=expr.columns)
    n_genes = expr.shape[0]
    ranked = expr.rank(axis=0, method="average", ascending=True) * (10000.0 / n_genes)
    gene_set = set(present)
    scores = {}
    for sample in expr.columns:
        m = ranked[sample]
        order = m.sort_values(ascending=False).index
        m_ord = m.loc[order].to_numpy(float)
        hits = np.fromiter((g in gene_set for g in order), dtype=bool, count=len(order))
        w = np.abs(m_ord) ** tau
        w_hit = np.where(hits, w, 0.0)
        nhit = float(w_hit.sum())
        nmiss = float((~hits).sum())
        if nhit <= 0 or nmiss <= 0:
            scores[sample] = np.nan
            continue
        p_hit = np.cumsum(w_hit) / nhit
        p_miss = np.cumsum((~hits).astype(float)) / nmiss
        scores[sample] = float(np.sum(p_hit - p_miss))
    return pd.Series(scores)


def parse_geo_matrix(path: Path):
    meta_rows = {}
    expr_start = None
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    samples = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics"):
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
                if ":" in v:
                    cleaned.append(v.split(":", 1)[1].strip())
                else:
                    cleaned.append(v)
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
        raise SystemExit(f"could not parse GEO matrix {path}")
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
        rows.append(
            [float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]]
        )
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def load_gpl_annot(path: Path) -> pd.Series:
    with gzip.open(path, "rt", errors="replace") as f:
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith("ID "):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    cands = [c for c in df.columns if "symbol" in c.lower()]
    if not cands:
        raise SystemExit(f"no symbol column in {path}")
    raw = df.set_index(id_col)[cands[0]].astype(str)
    raw = raw.replace({"nan": np.nan, "None": np.nan, "": np.nan}).dropna()
    n_sym = raw.map(lambda x: len([p.strip() for p in str(x).split("///") if p.strip()]))
    first = raw.map(lambda x: str(x).split("///")[0].strip())
    first = first[first.str.len() > 0]
    # unique-mapped only for named targets; first-symbol for ESTIMATE coverage
    return first, first[n_sym == 1]


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out


def load_gse31210() -> tuple[pd.DataFrame, pd.Series, pd.Series, dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    mat = download(GSE_MATRIX, CACHE / "GSE31210_series_matrix.txt.gz", min_bytes=100_000)
    annot = download(GPL570_ANNOT, CACHE / "GPL570.annot.gz", min_bytes=100_000)
    meta, probe_expr = parse_geo_matrix(mat)
    tissue_col = None
    for c in meta.columns:
        if "tissue" in c.lower():
            tissue_col = c
            break
    if tissue_col is None:
        raise SystemExit(f"no tissue column in GSE31210 meta: {list(meta.columns)}")
    tumor = meta.index[meta[tissue_col].astype(str).str.lower().str.contains("primary lung tumor")]
    probe_expr = probe_expr.loc[:, tumor]
    first_map, unique_map = load_gpl_annot(annot)
    # collapse for analysis genes with unique map where available; else first
    gene_expr = collapse_maxmean(probe_expr, first_map)

    est = pd.read_csv(EST_TSV, sep="\t")
    stromal = est.loc[est["set"] == "Stromal141_UP", "hugo"].tolist()
    immune_genes = est.loc[est["set"] == "Immune141_UP", "hugo"].tolist()
    strom = ssgsea(gene_expr, stromal)
    imm = ssgsea(gene_expr, immune_genes)
    est_score = strom + imm
    # purity covariate = ESTIMATEScore (monotone with cosine purity on this matrix)
    purity = est_score.copy()
    purity.name = "ESTIMATEScore"

    n_pur_ok = int(((np.cos(EST_A + EST_B * est_score.to_numpy(float)) >= 0) &
                    (np.cos(EST_A + EST_B * est_score.to_numpy(float)) <= 1)).sum())

    meta_out = {
        "cohort": "GSE31210_LUAD",
        "n": int(len(tumor)),
        "n_matrix": int(meta.shape[0]),
        "matrix_sha256": sha256_file(mat),
        "annot_sha256": sha256_file(annot),
        "matrix": "Affymetrix GPL570 MAS5 series matrix",
        "immune_score": "ESTIMATE_Immune141_ssGSEA",
        "estimate_immune_present": int(sum(g in gene_expr.index for g in immune_genes)),
        "estimate_stromal_present": int(sum(g in gene_expr.index for g in stromal)),
        "purity": "ESTIMATEScore (rank residual; cosine TumorPurity often outside [0,1])",
        "n_cosine_purity_in_0_1": n_pur_ok,
        "tissue_rule": "primary lung tumor",
    }
    return gene_expr.astype(float), purity, imm, meta_out


# ---------------------------------------------------------------------------
# per-cohort analysis
# ---------------------------------------------------------------------------


def analyze_cohort(
    name: str,
    expr: pd.DataFrame,
    purity: pd.Series,
    immune: pd.Series,
    gene_sets: dict[str, list[str]],
    set_meta: dict,
    meta: dict,
) -> dict:
    for g in ["TACSTD2", "CLDN4", "CD8A"]:
        if g not in expr.index:
            raise SystemExit(f"{name}: {g} missing")

    samples = list(expr.columns)
    tac = expr.loc["TACSTD2", samples].astype(float)
    cldn4 = expr.loc["CLDN4", samples].astype(float)
    cd8a = expr.loc["CD8A", samples].astype(float)
    imm = immune.loc[samples].astype(float)
    pur = purity.loc[samples].astype(float)

    high, low = quartile_split(tac)
    rank = rank_high_vs_low(expr, high, low)
    # drop TACSTD2 from ranking for GSEA (avoid trivial self-hit bias in leading edge)
    rank_gsea = rank.drop(labels=["TACSTD2"], errors="ignore")

    gsea = gsea_prerank(rank_gsea, gene_sets)
    gsea["cohort"] = name
    gsea["primary"] = gsea["term"].map(lambda t: bool(set_meta.get(t, {}).get("primary")))
    # BH within primary sets; also BH across all tested
    prim = gsea["primary"] == True  # noqa: E712
    gsea["fdr_bh_primary"] = np.nan
    gsea.loc[prim, "fdr_bh_primary"] = bh_fdr(gsea.loc[prim, "nom_p"]).values
    gsea["fdr_bh_all"] = bh_fdr(gsea["nom_p"])

    deg = deg_table(expr, high, low)
    deg["cohort"] = name
    focal = deg[deg["gene"].isin(FOCAL_TJ + ["TACSTD2", "CD8A", "CD274"])].copy()
    focal["cohort"] = name

    # immune correlations
    rows = []
    pairs = [
        ("TACSTD2", tac.to_numpy(float), "CD8A", cd8a.to_numpy(float)),
        ("TACSTD2", tac.to_numpy(float), "ImmuneScore", imm.to_numpy(float)),
        ("CLDN4", cldn4.to_numpy(float), "CD8A", cd8a.to_numpy(float)),
        ("CLDN4", cldn4.to_numpy(float), "ImmuneScore", imm.to_numpy(float)),
        ("TACSTD2", tac.to_numpy(float), "CLDN4", cldn4.to_numpy(float)),
        ("CD8A", cd8a.to_numpy(float), "ImmuneScore", imm.to_numpy(float)),
        ("TACSTD2", tac.to_numpy(float), "PURITY", pur.to_numpy(float)),
        ("CLDN4", cldn4.to_numpy(float), "PURITY", pur.to_numpy(float)),
    ]
    for pred, x, end, y in pairs:
        ru, pu, nu = spearman_pair(x, y)
        if end == "PURITY":
            rp, pp, np_ = np.nan, np.nan, nu
        else:
            rp, pp, np_ = partial_spearman(x, y, pur.to_numpy(float))
        rows.append(
            {
                "cohort": name,
                "predictor": pred,
                "endpoint": end,
                "n_unadjusted": nu,
                "unadj_rho": ru,
                "unadj_p": pu,
                "n_partial": np_,
                "partial_rho": rp,
                "partial_p": pp,
                "purity_covariate": meta["purity"],
                "immune_definition": meta["immune_score"],
            }
        )
    corr = pd.DataFrame(rows)

    per = pd.DataFrame(
        {
            "sample": samples,
            "TACSTD2": tac.values,
            "CLDN4": cldn4.values,
            "CD8A": cd8a.values,
            "ImmuneScore": imm.values,
            "purity_covariate": pur.values,
            "tacstd2_quartile": [
                "Q4" if s in high else ("Q1" if s in low else "mid") for s in samples
            ],
        }
    )
    per["cohort"] = name

    return {
        "name": name,
        "meta": meta,
        "n": len(samples),
        "n_q4": int(len(high)),
        "n_q1": int(len(low)),
        "gsea": gsea,
        "deg": deg,
        "focal": focal,
        "corr": corr,
        "per": per,
        "rank_head": rank.head(50).rename("welch_t").reset_index().rename(
            columns={"index": "gene", "Hugo_Symbol": "gene"}
        ),
        "rank_tail": rank.tail(50).rename("welch_t").reset_index().rename(
            columns={"index": "gene", "Hugo_Symbol": "gene"}
        ),
    }


# ---------------------------------------------------------------------------
# PPT corroboration (honest)
# ---------------------------------------------------------------------------


def ppt_verdicts(results: list[dict]) -> pd.DataFrame:
    """Map computed numbers onto PPT A1 / A8. Never invent."""
    rows = []
    for r in results:
        name = r["name"]
        corr = r["corr"]
        gsea = r["gsea"]
        prim = gsea[gsea["primary"] == True]  # noqa: E712

        def get_corr(pred, end):
            hit = corr[(corr.predictor == pred) & (corr.endpoint == end)]
            return hit.iloc[0] if len(hit) else None

        def get_gsea(term):
            hit = prim[prim.term == term]
            return hit.iloc[0] if len(hit) else None

        # A1: TACSTD2 vs CD8 / Immune negative after purity
        t2_cd8 = get_corr("TACSTD2", "CD8A")
        t2_imm = get_corr("TACSTD2", "ImmuneScore")
        a1_cd8_pass = (
            t2_cd8 is not None
            and t2_cd8.partial_rho < 0
            and t2_cd8.partial_p < 0.05
        )
        a1_imm_pass = (
            t2_imm is not None
            and t2_imm.partial_rho < 0
            and t2_imm.partial_p < 0.05
        )
        if a1_cd8_pass and a1_imm_pass:
            a1_call = "PASS"
        elif (t2_cd8 is not None and t2_cd8.unadj_rho < 0 and t2_cd8.unadj_p < 0.05) or (
            t2_imm is not None and t2_imm.unadj_rho < 0 and t2_imm.unadj_p < 0.05
        ):
            # unadjusted negative but partial fails
            a1_call = "PARTIAL"
        else:
            a1_call = "FAIL"

        rows.append(
            {
                "cohort": name,
                "ppt_claim": "A1",
                "ppt_text": "TACSTD2 vs immune negative after purity (CD8A + ImmuneScore)",
                "call": a1_call,
                "detail": (
                    f"TACSTD2–CD8A unadj ρ={fmt_rho(t2_cd8.unadj_rho)} p={fmt_p(t2_cd8.unadj_p)}; "
                    f"partial ρ={fmt_rho(t2_cd8.partial_rho)} p={fmt_p(t2_cd8.partial_p)}. "
                    f"TACSTD2–ImmuneScore unadj ρ={fmt_rho(t2_imm.unadj_rho)} p={fmt_p(t2_imm.unadj_p)}; "
                    f"partial ρ={fmt_rho(t2_imm.partial_rho)} p={fmt_p(t2_imm.partial_p)}."
                ),
                "n": r["n"],
            }
        )

        # A8: TJ/keratin UP; Hallmark EMT DOWN
        tj = get_gsea("KEGG_TIGHT_JUNCTION")
        apj = get_gsea("HALLMARK_APICAL_JUNCTION")
        krt = get_gsea("GOBP_KERATINIZATION")
        emt = get_gsea("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

        def up_ok(row):
            return (
                row is not None
                and row.nes > 0
                and (row.fdr_bh_primary < 0.05 if np.isfinite(row.fdr_bh_primary) else False)
            )

        def down_ok(row):
            return (
                row is not None
                and row.nes < 0
                and (row.fdr_bh_primary < 0.05 if np.isfinite(row.fdr_bh_primary) else False)
            )

        tj_up = up_ok(tj) or up_ok(apj)
        krt_up = up_ok(krt)
        emt_down = down_ok(emt)
        if tj_up and krt_up and emt_down:
            a8_call = "PASS"
        elif tj_up or krt_up:
            a8_call = "PARTIAL"  # barrier up without EMT-down (or only one of TJ/krt)
        else:
            a8_call = "FAIL"

        rows.append(
            {
                "cohort": name,
                "ppt_claim": "A8",
                "ppt_text": "TACSTD2-high: keratin/TJ up; Hallmark EMT down (GSEA)",
                "call": a8_call,
                "detail": (
                    f"KEGG_TJ NES={fmt_nes(tj.nes if tj is not None else np.nan)} "
                    f"FDR_prim={fmt_p(tj.fdr_bh_primary if tj is not None else np.nan)}; "
                    f"HALLMARK_APICAL_JUNCTION NES={fmt_nes(apj.nes if apj is not None else np.nan)} "
                    f"FDR_prim={fmt_p(apj.fdr_bh_primary if apj is not None else np.nan)}; "
                    f"GOBP_KERATINIZATION NES={fmt_nes(krt.nes if krt is not None else np.nan)} "
                    f"FDR_prim={fmt_p(krt.fdr_bh_primary if krt is not None else np.nan)}; "
                    f"HALLMARK_EMT NES={fmt_nes(emt.nes if emt is not None else np.nan)} "
                    f"FDR_prim={fmt_p(emt.fdr_bh_primary if emt is not None else np.nan)}."
                ),
                "n": r["n"],
            }
        )

        # companion: CLDN4 vs CD8/Immune (not a numbered PPT A1/A8 line, but requested)
        c4_cd8 = get_corr("CLDN4", "CD8A")
        c4_imm = get_corr("CLDN4", "ImmuneScore")
        c4_pass = (
            c4_cd8 is not None
            and c4_cd8.unadj_rho < 0
            and c4_cd8.unadj_p < 0.05
            and c4_imm is not None
            and c4_imm.unadj_rho < 0
            and c4_imm.unadj_p < 0.05
        )
        c4_partial_pass = (
            c4_cd8 is not None
            and c4_cd8.partial_rho < 0
            and c4_cd8.partial_p < 0.05
            and c4_imm is not None
            and c4_imm.partial_rho < 0
            and c4_imm.partial_p < 0.05
        )
        if c4_partial_pass:
            c4_call = "PASS"
        elif c4_pass:
            c4_call = "PARTIAL"
        else:
            c4_call = "FAIL"
        rows.append(
            {
                "cohort": name,
                "ppt_claim": "CLDN4_immune_companion",
                "ppt_text": "CLDN4 vs CD8A / ImmuneScore negative (requested companion)",
                "call": c4_call,
                "detail": (
                    f"CLDN4–CD8A unadj ρ={fmt_rho(c4_cd8.unadj_rho)} p={fmt_p(c4_cd8.unadj_p)}; "
                    f"partial ρ={fmt_rho(c4_cd8.partial_rho)} p={fmt_p(c4_cd8.partial_p)}. "
                    f"CLDN4–ImmuneScore unadj ρ={fmt_rho(c4_imm.unadj_rho)} p={fmt_p(c4_imm.unadj_p)}; "
                    f"partial ρ={fmt_rho(c4_imm.partial_rho)} p={fmt_p(c4_imm.partial_p)}."
                ),
                "n": r["n"],
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# figures + FINDING
# ---------------------------------------------------------------------------


def make_figures(results: list[dict], verdicts: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    # 1. GSEA NES primary
    fig, axes = plt.subplots(1, len(results), figsize=(5.2 * len(results), 5.5), sharey=True)
    if len(results) == 1:
        axes = [axes]
    for ax, r in zip(axes, results):
        g = r["gsea"]
        prim = g[g["primary"] == True].copy()  # noqa: E712
        prim = prim.sort_values("nes")
        colors = ["#1b7837" if n > 0 else "#c51b7d" for n in prim["nes"]]
        ax.barh(prim["term"], prim["nes"], color=colors, alpha=0.85)
        ax.axvline(0, color="#333", lw=0.8)
        ax.set_xlabel("NES (TACSTD2 Q4 vs Q1)")
        ax.set_title(f"{r['name']}  n={r['n']} ({r['n_q4']} vs {r['n_q1']})")
        for i, (_, row) in enumerate(prim.iterrows()):
            star = ""
            if np.isfinite(row.fdr_bh_primary):
                if row.fdr_bh_primary < 0.01:
                    star = "***"
                elif row.fdr_bh_primary < 0.05:
                    star = "**"
                elif row.fdr_bh_primary < 0.25:
                    star = "*"
            if star:
                x = row.nes + (0.05 if row.nes >= 0 else -0.05)
                ax.text(x, i, star, va="center", ha="left" if row.nes >= 0 else "right", fontsize=8)
    fig.suptitle("Primary TJ / keratin / EMT GSEA — paper funnel", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_gsea_primary_nes.png", dpi=150)
    plt.close(fig)

    # 2. Immune forest
    fig, axes = plt.subplots(1, len(results), figsize=(5.0 * len(results), 4.2), sharey=True)
    if len(results) == 1:
        axes = [axes]
    labels = [
        ("TACSTD2", "CD8A"),
        ("TACSTD2", "ImmuneScore"),
        ("CLDN4", "CD8A"),
        ("CLDN4", "ImmuneScore"),
    ]
    for ax, r in zip(axes, results):
        corr = r["corr"]
        ylabels = []
        unadj = []
        partial = []
        for pred, end in labels:
            hit = corr[(corr.predictor == pred) & (corr.endpoint == end)].iloc[0]
            ylabels.append(f"{pred} vs {end}")
            unadj.append(hit.unadj_rho)
            partial.append(hit.partial_rho)
        y = np.arange(len(ylabels))
        ax.axvline(0, color="#333", lw=0.8)
        ax.scatter(unadj, y - 0.12, c="#4c78a8", s=50, label="unadjusted", zorder=3)
        ax.scatter(partial, y + 0.12, c="#f58518", s=50, label="partial | purity", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(ylabels)
        ax.set_xlabel("Spearman ρ")
        ax.set_title(f"{r['name']}  n={r['n']}")
        ax.set_xlim(-0.7, 0.3)
        ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("TACSTD2 / CLDN4 vs CD8A / ImmuneScore", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_immune_spearman_forest.png", dpi=150)
    plt.close(fig)

    # 3. scatter grid OncoSG / GSE31210 TACSTD2 vs CD8A and CLDN4 vs CD8A
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.5))
    for col, r in enumerate(results):
        per = r["per"]
        for row, (xcol, ycol, color) in enumerate(
            [("TACSTD2", "CD8A", "#4c78a8"), ("CLDN4", "CD8A", "#54a24b")]
        ):
            ax = axes[row, col]
            x = per[xcol].to_numpy(float)
            y = per[ycol].to_numpy(float)
            ax.scatter(x, y, s=12, alpha=0.55, c=color, linewidths=0)
            rho, p, n = spearman_pair(x, y)
            ax.set_xlabel(xcol)
            ax.set_ylabel(ycol)
            ax.set_title(f"{r['name']}: ρ={fmt_rho(rho)} p={fmt_p(p)} n={n}", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_scatter_vs_cd8a.png", dpi=150)
    plt.close(fig)

    # 4. verdict heatmap-ish table as bars
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    call_map = {"PASS": 1, "PARTIAL": 0.5, "FAIL": 0}
    color_map = {"PASS": "#1b7837", "PARTIAL": "#e08214", "FAIL": "#c51b7d"}
    v = verdicts.copy()
    v["y"] = v["ppt_claim"] + " | " + v["cohort"]
    ax.barh(v["y"], [call_map[c] for c in v["call"]], color=[color_map[c] for c in v["call"]])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("call (PASS=1, PARTIAL=0.5, FAIL=0)")
    ax.set_title("PPT corroboration (computed; not fabricated)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_ppt_verdicts.png", dpi=150)
    plt.close(fig)


def write_finding(results: list[dict], verdicts: pd.DataFrame) -> None:
    by = {r["name"]: r for r in results}

    def corr_cell(cohort, pred, end, which="unadj"):
        c = by[cohort]["corr"]
        hit = c[(c.predictor == pred) & (c.endpoint == end)].iloc[0]
        if which == "unadj":
            return f"{fmt_rho(hit.unadj_rho)} ({fmt_p(hit.unadj_p)})"
        return f"{fmt_rho(hit.partial_rho)} ({fmt_p(hit.partial_p)})"

    def gsea_cell(cohort, term):
        g = by[cohort]["gsea"]
        hit = g[(g.term == term) & (g.primary == True)]  # noqa: E712
        if hit.empty:
            return "NA"
        row = hit.iloc[0]
        return f"NES={fmt_nes(row.nes)} FDR={fmt_p(row.fdr_bh_primary)}"

    lines = []
    lines.append("# PAPER FUNNEL — Bulk LUAD TACSTD2-high vs low (OncoSG + GSE31210)")
    lines.append("")
    lines.append("**Rule:** never fabricate. Every number below was computed in `analyze.py`.")
    lines.append("PPT claims A1 / A8 are corroborated or not from these numbers; empty/null is reported as such.")
    lines.append("")
    lines.append("## Cohorts")
    lines.append("")
    lines.append("| cohort | n | Q4 vs Q1 | matrix | ImmuneScore | purity covariate |")
    lines.append("|---|---:|---|---|---|---|")
    for r in results:
        m = r["meta"]
        lines.append(
            f"| {r['name']} | **{r['n']}** | {r['n_q4']} vs {r['n_q1']} | {m['matrix']} | "
            f"{m['immune_score']} | {m['purity']} |"
        )
    lines.append("")
    lines.append(
        "OncoSG portal RNA list is 181; public z-score matrix has **169** columns. "
        "Do not write n=181. ESTIMATE is **skipped** on OncoSG (z-scores only)."
    )
    lines.append("")
    lines.append("## PPT corroboration")
    lines.append("")
    lines.append("| cohort | claim | call | detail |")
    lines.append("|---|---|---|---|")
    for _, row in verdicts.iterrows():
        lines.append(
            f"| {row.cohort} | {row.ppt_claim} | **{row.call}** | {row.detail} |"
        )
    lines.append("")
    lines.append("### Call rules (locked)")
    lines.append("")
    lines.append(
        "- **A1 PASS:** TACSTD2 vs CD8A **and** vs ImmuneScore both have partial ρ < 0 and partial p < 0.05."
    )
    lines.append(
        "- **A1 PARTIAL:** unadjusted Spearmans are negative and significant for at least one endpoint, "
        "but purity-adjusted tests do not both pass."
    )
    lines.append(
        "- **A8 PASS:** KEGG TJ or Hallmark apical junction UP (NES>0, primary FDR<0.05) **and** "
        "GOBP keratinization UP **and** Hallmark EMT DOWN."
    )
    lines.append(
        "- **A8 PARTIAL:** keratin and/or TJ/apical UP at primary FDR<0.05, but the full triad "
        "(TJ-or-apical UP + keratin UP + Hallmark EMT DOWN) does not hold "
        "(includes cases where Hallmark apical junction is DOWN, or EMT fails)."
    )
    lines.append(
        "- **CLDN4 companion PASS:** CLDN4 vs CD8A and vs ImmuneScore both negative after purity; "
        "PARTIAL = unadjusted only."
    )
    lines.append("")
    lines.append("## TACSTD2 / CLDN4 vs ImmuneScore / CD8A")
    lines.append("")
    lines.append(
        "| cohort | TACSTD2–CD8A ρ (p) | TACSTD2–Immune ρ (p) | "
        "TACSTD2–CD8A partial (p) | TACSTD2–Immune partial (p) | "
        "CLDN4–CD8A ρ (p) | CLDN4–Immune ρ (p) |"
    )
    lines.append("|---|---|---|---|---|---|---|")
    for r in results:
        n = r["name"]
        lines.append(
            f"| {n} | {corr_cell(n,'TACSTD2','CD8A')} | {corr_cell(n,'TACSTD2','ImmuneScore')} | "
            f"{corr_cell(n,'TACSTD2','CD8A','partial')} | {corr_cell(n,'TACSTD2','ImmuneScore','partial')} | "
            f"{corr_cell(n,'CLDN4','CD8A')} | {corr_cell(n,'CLDN4','ImmuneScore')} |"
        )
    lines.append("")
    lines.append(
        "ImmuneScore definitions differ by cohort (A1 8-gene mean z on OncoSG; "
        "ESTIMATE Immune141 on GSE31210). Do not pool the two ImmuneScore ρ values."
    )
    lines.append("")
    lines.append("## GSEA primary sets (TACSTD2 Q4 vs Q1, Welch prerank)")
    lines.append("")
    lines.append("| set | want | OncoSG | GSE31210 |")
    lines.append("|---|---|---|---|")
    want = {
        "KEGG_TIGHT_JUNCTION": "UP",
        "HALLMARK_APICAL_JUNCTION": "UP",
        "GOBP_KERATINIZATION": "UP",
        "GOBP_TIGHT_JUNCTION_ORGANIZATION": "UP",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "DOWN",
        "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION": "DOWN",
    }
    for term, w in want.items():
        lines.append(
            f"| {term} | {w} | {gsea_cell('OncoSG_LUAD', term)} | {gsea_cell('GSE31210_LUAD', term)} |"
        )
    lines.append("")
    lines.append("BH-FDR is within the 12 primary sets per cohort. Full table: `tables/gsea_primary.tsv`.")
    lines.append("")
    lines.append("## Focal TJ genes (Welch delta, Q4 − Q1)")
    lines.append("")
    lines.append(
        "Delta units differ: OncoSG = z-score mean difference; GSE31210 = MAS5 intensity mean difference "
        "(not log2FC). Do not compare absolute delta magnitudes across cohorts."
    )
    lines.append("")
    lines.append("| gene | OncoSG delta / FDR | GSE31210 delta / FDR |")
    lines.append("|---|---|---|")
    for g in FOCAL_TJ:
        cells = []
        for r in results:
            hit = r["focal"][r["focal"].gene == g]
            if hit.empty:
                cells.append("absent")
            else:
                row = hit.iloc[0]
                cells.append(
                    f"{row.delta_high_minus_low:+.3f} / {fmt_p(row.fdr_bh)}"
                )
        lines.append(f"| {g} | {cells[0]} | {cells[1]} |")
    lines.append("")
    lines.append("## What this is / is not")
    lines.append("")
    lines.append("- East-Asian surgical LUAD (OncoSG + GSE31210). **Not** an ICI-response cohort.")
    lines.append("- Does **not** re-audit TCGA-LUAD A8 (see PR 117): LUAD Hallmark EMT can be UP.")
    lines.append("- Does **not** invent ESTIMATE on OncoSG z-scores.")
    lines.append("- Does **not** merge private 8KL scRNA.")
    lines.append("- PPT numbers are corroborated or contradicted by **these** public matrices only.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r methods/paper_funnel_luad_tacstd2_tj/requirements.txt")
    lines.append("python3 methods/paper_funnel_luad_tacstd2_tj/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/immune_correlations.tsv`, `gsea_primary.tsv`, `gsea_all.tsv`")
    lines.append("- `tables/deg_focal.tsv`, `ppt_verdicts.tsv`, `per_sample.tsv`, `summary.json`")
    lines.append("- `figures/fig1_gsea_primary_nes.png` … `fig4_ppt_verdicts.png`")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    with SET_JSON.open() as f:
        set_blob = json.load(f)
    gene_sets = set_blob["sets"]
    set_meta = set_blob["meta"]
    # strip TACSTD2 from every set if present
    for k, members in list(gene_sets.items()):
        gene_sets[k] = [g for g in members if g != "TACSTD2"]

    print("loading OncoSG…", flush=True)
    expr_o, pur_o, imm_o, meta_o = load_oncosg()
    print(f"  OncoSG n={meta_o['n']}", flush=True)

    print("loading GSE31210…", flush=True)
    expr_g, pur_g, imm_g, meta_g = load_gse31210()
    print(f"  GSE31210 n={meta_g['n']}", flush=True)

    results = []
    print("analyzing OncoSG (DEG + GSEA + immune)…", flush=True)
    results.append(
        analyze_cohort("OncoSG_LUAD", expr_o, pur_o, imm_o, gene_sets, set_meta, meta_o)
    )
    print("analyzing GSE31210 (DEG + GSEA + immune)…", flush=True)
    results.append(
        analyze_cohort("GSE31210_LUAD", expr_g, pur_g, imm_g, gene_sets, set_meta, meta_g)
    )

    corr_all = pd.concat([r["corr"] for r in results], ignore_index=True)
    gsea_all = pd.concat([r["gsea"] for r in results], ignore_index=True)
    gsea_prim = gsea_all[gsea_all["primary"] == True].copy()  # noqa: E712
    focal_all = pd.concat([r["focal"] for r in results], ignore_index=True)
    per_all = pd.concat([r["per"] for r in results], ignore_index=True)
    verdicts = ppt_verdicts(results)

    # DEG top tables (per cohort, top 200 up/down by |t|)
    deg_tops = []
    for r in results:
        d = r["deg"]
        top = pd.concat([d.head(100), d.tail(100)])
        top = top.drop_duplicates(subset=["gene"])
        deg_tops.append(top)
    deg_top = pd.concat(deg_tops, ignore_index=True)

    corr_all.to_csv(TABLES / "immune_correlations.tsv", sep="\t", index=False)
    gsea_prim.to_csv(TABLES / "gsea_primary.tsv", sep="\t", index=False)
    gsea_all.to_csv(TABLES / "gsea_all.tsv", sep="\t", index=False)
    focal_all.to_csv(TABLES / "deg_focal.tsv", sep="\t", index=False)
    deg_top.to_csv(TABLES / "deg_top200_per_tail.tsv", sep="\t", index=False)
    per_all.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)
    verdicts.to_csv(TABLES / "ppt_verdicts.tsv", sep="\t", index=False)

    summary = {
        "task": "PAPER_FUNNEL_bulk_LUAD_TACSTD2_TJ_immune",
        "ppt_source": "User PPT 2026-08-17 (金炫宏); claims A1, A8 + CLDN4 companion",
        "rule": "never fabricate",
        "cohorts": [r["meta"] | {"n_q4": r["n_q4"], "n_q1": r["n_q1"]} for r in results],
        "gsea": {
            "split": "TACSTD2 quartile",
            "rank": "Welch t high-minus-low",
            "nperm": NPERM,
            "seed": SEED,
            "primary_sets": sorted([k for k, v in set_meta.items() if v.get("primary")]),
        },
        "verdicts": verdicts.to_dict(orient="records"),
        "immune_correlations": corr_all.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    make_figures(results, verdicts)
    write_finding(results, verdicts)

    print("\n=== PPT verdicts ===", flush=True)
    print(verdicts.to_string(index=False), flush=True)
    print("\nWrote FINDING.md + tables/ + figures/", flush=True)


if __name__ == "__main__":
    main()
