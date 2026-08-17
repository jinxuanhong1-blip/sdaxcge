#!/usr/bin/env python3
"""GSE68465 LUAD: CLDN4 vs CD274 and HLA-A/B/C after ESTIMATE.

Additive only. New CD274 / classical MHC-I cut.
CLDN4 vs CD8A after ESTIMATE TumorPurity is already known (PR 298:
n=443, unadj ρ=−0.177, partial ρ=−0.092, p=0.054) and is not re-cut.

Public GEO series matrix (Shedden / Director's Challenge; GPL96 U133A)
+ official GPL96 annotation + Yoshihara 2013 ESTIMATE lists.

CD274 (Entrez 29126) is not on U133A. Do not invent a PD-L1 ρ.
PDCD1LG2 (PD-L2) is on the chip and is not substituted.
"""

from __future__ import annotations

import gzip
import io
import json
import math
import os
import urllib.request
from datetime import datetime, timezone
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
DATA = Path(os.environ.get("GSE68465_CLDN4_CD274_DATA", "/tmp/gse68465_cldn4_cd274"))
SIG = HERE / "estimate_yoshihara_2013.tsv"

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68465/"
    "matrix/GSE68465_series_matrix.txt.gz"
)
GPL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"
UA = "sdaxcge-gse68465-cldn4-cd274/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"

EST_A = 0.6049872018
EST_B = 0.0001467884
BOOT_SEED = 20260817
N_BOOT = 2000

# New-cut named genes. CD8A is taken as given (PR 298) and is not a target.
CLAIM_GENES = ["CD274", "HLA-A", "HLA-B", "HLA-C"]
ALIASES = {
    "CLDN4": ["CLDN4"],
    "CD274": ["CD274", "PDCD1LG1"],
    "HLA-A": ["HLA-A"],
    "HLA-B": ["HLA-B"],
    "HLA-C": ["HLA-C"],
    "CD8A": ["CD8A"],
}
ENTREZ = {
    "CLDN4": "1364",
    "CD274": "29126",
    "HLA-A": "3105",
    "HLA-B": "3106",
    "HLA-C": "3107",
    "CD8A": "925",
    "PDCD1LG2": "80380",
}
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
}

# Given from PR 298; not recomputed as a claim.
CD8_GIVEN = {
    "pr": 298,
    "n": 443,
    "unadj_rho": -0.177,
    "unadj_p": 1.7e-4,
    "partial_rho": -0.092,
    "partial_p": 0.054,
    "covariate": "ESTIMATE TumorPurity",
}


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10_000:
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
    sign = "+" if r > 0 else "−" if r < 0 else ""
    return f"{sign}{abs(r):.3f}"


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
        return "—"
    return f"{fmt_rho(lo)} to {fmt_rho(hi)}"


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
    zr = stats.rankdata(z[m])
    Z = zr.reshape(-1, 1)
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def bootstrap_spearman(x, y, n_boot=N_BOOT, seed=BOOT_SEED):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[m], y[m]
    n = int(len(xx))
    if n < 8:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(xx[idx]).size < 3 or np.unique(yy[idx]).size < 3:
            rhos[i] = np.nan
            continue
        r, _ = stats.spearmanr(xx[idx], yy[idx])
        rhos[i] = r
    lo, hi = np.nanquantile(rhos, [0.025, 0.975])
    return float(lo), float(hi)


def bootstrap_partial(x, y, z, n_boot=N_BOOT, seed=BOOT_SEED):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    xx, yy, zz = x[m], y[m], z[m]
    n = int(len(xx))
    if n < 8:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        r, _, _ = partial_spearman(xx[idx], yy[idx], zz[idx])
        rhos[i] = r
    lo, hi = np.nanquantile(rhos, [0.025, 0.975])
    return float(lo), float(hi)


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
        if line.startswith("!Sample_title"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["title"] = vals
        if line.startswith("!Sample_source_name"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["source"] = vals
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
        rows.append([float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def parse_soft_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(
        i
        for i, l in enumerate(lines)
        if l.startswith("!platform_table_begin") or l.startswith("ID\t") or l.startswith("ID ")
    )
    if lines[start].startswith("!platform_table_begin"):
        start += 1
    end = next(
        (i for i, l in enumerate(lines) if l.startswith("!platform_table_end")),
        len(lines),
    )
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def first_symbol(s: str) -> str:
    if not isinstance(s, str) or not s or s == "nan":
        return ""
    return s.split("///")[0].strip()


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out, pick


def search_annot_exact(annot: pd.DataFrame, aliases: list[str], entrez: str | None = None) -> pd.DataFrame:
    """Exact first-symbol or exact Entrez. No substring aliases (PDL1 ⊂ SPDL1)."""
    cols = [c for c in ["ID", "Gene symbol", "Gene title", "Gene ID"] if c in annot.columns]
    first = annot["Gene symbol"].map(first_symbol)
    mask = first.isin(aliases)
    if entrez:
        mask = mask | (annot["Gene ID"].fillna("").astype(str) == str(entrez))
    return annot.loc[mask, cols].copy()


def verdict(rho, p, n, present: bool) -> str:
    if not present:
        return "ABSENT"
    if n < 40:
        return "UNDERPOWERED"
    if not np.isfinite(rho) or not np.isfinite(p):
        return "NO_EVIDENCE"
    if p < 0.05 and rho < 0:
        return "NEGATIVE"
    if p < 0.05 and rho > 0:
        return "POSITIVE"
    return "NO_EVIDENCE"


def rec_row(predictor, endpoint, x, y, z, covar_name, note="", present=True):
    if not present or x is None or y is None:
        return {
            "predictor": predictor,
            "endpoint": endpoint,
            "present": False,
            "n_unadjusted": 0,
            "unadj_rho": np.nan,
            "unadj_p": np.nan,
            "unadj_ci_lo": np.nan,
            "unadj_ci_hi": np.nan,
            "n_partial": 0,
            "partial_rho": np.nan,
            "partial_p": np.nan,
            "partial_ci_lo": np.nan,
            "partial_ci_hi": np.nan,
            "covariate": covar_name,
            "verdict_unadj": "ABSENT",
            "verdict_partial": "ABSENT",
            "note": note or "CD274 not on GPL96 U133A",
        }
    ru, pu, nu = spearman_pair(x, y)
    lo, hi = bootstrap_spearman(x, y)
    if z is None:
        rp, pp, np_ = np.nan, np.nan, nu
        plo, phi = np.nan, np.nan
    else:
        rp, pp, np_ = partial_spearman(x, y, z)
        plo, phi = bootstrap_partial(x, y, z)
    return {
        "predictor": predictor,
        "endpoint": endpoint,
        "present": True,
        "n_unadjusted": nu,
        "unadj_rho": ru,
        "unadj_p": pu,
        "unadj_ci_lo": lo,
        "unadj_ci_hi": hi,
        "n_partial": np_,
        "partial_rho": rp,
        "partial_p": pp,
        "partial_ci_lo": plo,
        "partial_ci_hi": phi,
        "covariate": covar_name,
        "verdict_unadj": verdict(ru, pu, nu, True),
        "verdict_partial": verdict(rp, pp, np_, True) if z is not None else "",
        "note": note,
    }


def mwu_q4q1(x: np.ndarray, y: np.ndarray | None, present: bool) -> dict:
    if not present or y is None:
        return {
            "present": False,
            "n_Q4": 0,
            "n_Q1": 0,
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
            "verdict": "ABSENT",
        }
    xs = pd.Series(x)
    ys = pd.Series(y)
    q = pd.qcut(xs, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    a = ys[q == "Q4"].dropna()
    b = ys[q == "Q1"].dropna()
    if a.size < 4 or b.size < 4:
        return {
            "present": True,
            "n_Q4": int(a.size),
            "n_Q1": int(b.size),
            "median_Q4": np.nan,
            "median_Q1": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial_Q4_minus_Q1": np.nan,
            "verdict": "UNDERPOWERED",
        }
    U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
    rbc = 2 * U / (a.size * b.size) - 1
    return {
        "present": True,
        "n_Q4": int(a.size),
        "n_Q1": int(b.size),
        "median_Q4": float(np.median(a)),
        "median_Q1": float(np.median(b)),
        "U": float(U),
        "p": float(p_mw),
        "rank_biserial_Q4_minus_Q1": float(rbc),
        "verdict": verdict(float(rbc), float(p_mw), int(a.size + b.size), True),
    }


def write_finding(stats_df, cov, probe_used, probe_rule, n_tumor, n_matrix, n_normal, hl_df) -> None:
    est_name = "ESTIMATEScore (Stromal+Immune)"

    def pick(end, covar=est_name):
        hit = stats_df[(stats_df.predictor == "CLDN4") & (stats_df.endpoint == end) & (stats_df.covariate == covar)]
        if hit.empty:
            hit = stats_df[(stats_df.predictor == "CLDN4") & (stats_df.endpoint == end)]
        return hit.iloc[0] if len(hit) else None

    c4_pd = pick("CD274")
    c4_a = pick("HLA-A")
    c4_b = pick("HLA-B")
    c4_c = pick("HLA-C")
    c4_mhc = pick("MHC1_HLAabc")
    c4_imm = pick("ImmuneScore")

    def cell(row, which="unadj"):
        if row is None or (which == "unadj" and not bool(row.present)):
            return "ABSENT (n=0)"
        if which == "unadj":
            return f"{fmt_rho(row.unadj_rho)} ({fmt_p(row.unadj_p)})"
        if not bool(row.present):
            return "ABSENT (n=0)"
        return f"{fmt_rho(row.partial_rho)} ({fmt_p(row.partial_p)})"

    def pair_row(label, row):
        if row is None or not bool(row.present):
            return f"| CLDN4 vs **{label}** | **0** | — | — | — | — | — | **ABSENT** |"
        return (
            f"| CLDN4 vs {label} | **{int(row.n_unadjusted)}** | {fmt_rho(row.unadj_rho)} | "
            f"{fmt_ci(row.unadj_ci_lo, row.unadj_ci_hi)} | {fmt_p(row.unadj_p)} | "
            f"{fmt_rho(row.partial_rho)} ({fmt_p(row.partial_p)}) | "
            f"{fmt_ci(row.partial_ci_lo, row.partial_ci_hi)} | **{row.verdict_unadj}** |"
        )

    def hl_line(name):
        hit = hl_df[hl_df.endpoint == name]
        if hit.empty:
            return f"| {name} | — | — | — | **ABSENT** |"
        r = hit.iloc[0]
        if not bool(r.present):
            return f"| {name} | — | — | — | **ABSENT** |"
        return (
            f"| {name} | {int(r.n_Q4)} vs {int(r.n_Q1)} | "
            f"{fmt_rho(r.rank_biserial_Q4_minus_Q1)} | {fmt_p(r.p)} | {r.verdict} |"
        )

    n_pur_ok = cov.get("n_purity_in_0_1", "NA")
    n_pur_out = cov.get("n_purity_outside_0_1", "NA")
    strom_n = cov["estimate_stromal_present"]
    imm_n = cov["estimate_immune_present"]

    md = f"""# GSE68465 LUAD — CLDN4 vs CD274 and HLA-A/B/C (partial on ESTIMATE)

**Additive only. New CD274 / classical MHC-I cut.** CLDN4 vs CD8A after ESTIMATE TumorPurity is **already known** ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298): n=443, unadj ρ=−0.177, partial ρ=−0.092, p=0.054) and is **not** re-cut. This folder only measures CLDN4 against **CD274** and **HLA-A / HLA-B / HLA-C**.

Public Director's Challenge LUAD microarray (Shedden et al., *Nat Med* 2008; GEO [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465); GPL96 Affymetrix Human Genome U133A). Unit is the **tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {n_matrix} | GEO `GSE68465_series_matrix.txt.gz` (22,283 probes) |
| LUAD tumors | **{n_tumor}** | `disease_state` contains `Adenocarcinoma` |
| Normal arrays | {n_normal} | dropped |
| CLDN4 finite (`201428_at`) | {n_tumor} | official GPL96 = CLDN4 / Entrez 1364 |
| **CD274 finite (Entrez 29126)** | **0** | no GPL96 probe; Plus-2 `223834_at` / `227458_at` are not on U133A |
| HLA-A / HLA-B / HLA-C finite | {n_tumor} | first-symbol max-mean; 3/3 present |
| ESTIMATEScore | {n_tumor} | Yoshihara Stromal141+Immune141 ssGSEA ({strom_n}/141 + {imm_n}/141) |
| ESTIMATE TumorPurity in [0, 1] | {n_pur_ok} | cosine map; {n_pur_out} tumors wrap outside [0, 1] |
| ICI labels | 0 | surgical / multi-site prognostic series |
| **Primary pairwise n (CLDN4 + CD274)** | **0** | ABSENT — do not invent a PD-L1 ρ |
| **Primary pairwise n (CLDN4 + HLA-A/B/C)** | **{n_tumor}** | this is the n used for HLA / MHC-I |

Primary HLA tests use **n={n_tumor}** tumors. That is the same tumor rule as the sibling CD8 cut ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298)); it is not a new cohort hunt.

PDCD1LG2 (PD-L2, `220049_s_at`, Entrez 80380) **is** on U133A. It is not CD274 and is **not** substituted.

## One-row table

| dataset | n tumors | CLDN4 | CD274 | HLA-A/B/C | CLDN4–CD274 ρ (p) | CLDN4–HLA-A ρ (p) | CLDN4–HLA-B ρ (p) | CLDN4–HLA-C ρ (p) | CLDN4–CD274 partial \\| ESTIMATE (p) | CLDN4–HLA-A partial \\| ESTIMATE (p) | CLDN4–HLA-B partial \\| ESTIMATE (p) | CLDN4–HLA-C partial \\| ESTIMATE (p) |
|---|---:|---|---|---|---|---|---|---|---|---|---|---|
| GSE68465 LUAD tumors | {n_tumor} | `{probe_used.get('CLDN4', '201428_at')}` | **ABSENT** (Entrez 29126 not on U133A) | 3/3 | **ABSENT (n=0)** | {cell(c4_a)} | {cell(c4_b)} | {cell(c4_c)} | **ABSENT (n=0)** | {cell(c4_a, "partial")} | {cell(c4_b, "partial")} | {cell(c4_c, "partial")} |

MHC-I mean-z (HLA-A/B/C only): unadj {cell(c4_mhc)}; partial \\| ESTIMATE {cell(c4_mhc, "partial")}.

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## New CD274 cut

MAS5 as deposited. Tests use native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. MHC-I mean-z = gene-wise z of **HLA-A/B/C only** (3/3; not B2M/TAP). Partial = Pearson of rank residuals on ESTIMATEScore; df = n − 3. On this matrix TumorPurity is a monotone cosine of ESTIMATEScore (when the angle stays in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore for tumors with a defined cosine.

| pair | n | ρ | 95% CI | p | ρ \\| ESTIMATE (p) | 95% CI partial | verdict |
|---|---:|---:|---|---:|---|---|---|
{pair_row("**CD274**", c4_pd)}
{pair_row("HLA-A", c4_a)}
{pair_row("HLA-B †", c4_b)}
{pair_row("HLA-C", c4_c)}
{pair_row("MHC-I (HLA-A/B/C mean-z)", c4_mhc)}

† **HLA-B is in Yoshihara Immune141.** Residualizing ESTIMATEScore (stromal + immune) out of HLA-B is partly circular. HLA-A and HLA-C are not in Immune141. The HLA-B row is reported because it was requested; it is not an independent infiltrate control.

CLDN4 Q4 vs Q1 on the same endpoints:

| endpoint | n Q4 vs Q1 | rank-biserial (Q4−Q1) | MWU p | verdict |
|---|---|---:|---:|---|
{hl_line("CD274")}
{hl_line("HLA-A")}
{hl_line("HLA-B")}
{hl_line("HLA-C")}
{hl_line("MHC1_HLAabc")}

There is no public CD274 (PD-L1) measurement on this U133A series. Classical MHC-I is present at honest n={n_tumor}. Antigen-presentation / PD-L1 cut only; not an ICI-response test. CLDN4–CD8A is not restated.

## ESTIMATE

R `estimate` is not required. Scores follow the sibling public implementation ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298)):

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

**Partial covariate is ESTIMATEScore**, not ImmuneScore. ImmuneScore is one addend of ESTIMATEScore and is collinear with it. The cosine is below 0 for **{n_pur_out} / {n_tumor}** tumors and is **not** used as a 0–1 fraction. The [0, 1] subset (n={n_pur_ok}) is only the least-impure tail and is **not** the analysis n.

CLDN4 vs ImmuneScore is reported only as a context row (same matrix as the sibling cut). It is not a new CD8 test.

| predictor | endpoint | n | ρ | p | partial ρ \\| ESTIMATE | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | ImmuneScore | {int(c4_imm.n_unadjusted) if c4_imm is not None else "NA"} | {fmt_rho(c4_imm.unadj_rho) if c4_imm is not None else "NA"} | {fmt_p(c4_imm.unadj_p) if c4_imm is not None else "NA"} | {fmt_rho(c4_imm.partial_rho) if c4_imm is not None else "NA"} | {fmt_p(c4_imm.partial_p) if c4_imm is not None else "NA"} |

## Given CD8 residual (not re-run)

From [PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298), GSE68465 CLDN4 vs CD8A after the same ESTIMATE TumorPurity residual is **partial ρ = −0.092, p = 0.054, n = 443** (unadj −0.177, p = 1.7e-4). That CD8 table is not re-downloaded or re-fit here.

## Named-gene collapse

First-symbol max-mean (same as the sibling CD8 cut). Exact first-symbol / Entrez only — no substring aliases (`PDL1` hits `SPDL1`).

| gene | probe | rule |
|---|---|---|
| CLDN4 | `{probe_used.get('CLDN4', 'NA')}` | {probe_rule.get('CLDN4', 'NA')} |
| CD274 | **ABSENT** | no GPL96 hit for CD274 / PDCD1LG1 / Entrez 29126 |
| HLA-A | `{probe_used.get('HLA-A', 'NA')}` | {probe_rule.get('HLA-A', 'NA')} |
| HLA-B | `{probe_used.get('HLA-B', 'NA')}` | {probe_rule.get('HLA-B', 'NA')} |
| HLA-C | `{probe_used.get('HLA-C', 'NA')}` | {probe_rule.get('HLA-C', 'NA')} |

ESTIMATE gene collapse uses the same first-symbol mapping.

## What is not done

- No re-cut of CLDN4 vs CD8A / TACSTD2 (sibling folder / PR 298 / PR 235).
- No PD-L1 RNA claim (CD274 is not on the chip).
- No PD-L2 (PDCD1LG2) substitute for CD274.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_confirm.tsv`, `coverage.tsv`, `highlow_cldn4_cd274_hla.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png` — ABSENT call
- `figures/fig2_cldn4_vs_hla.png`
- `figures/fig3_spearman_forest.png`
- `figures/fig4_correlation_heatmap.png`

```bash
python3 -m pip install -r methods/gse68465_cldn4_cd274/requirements.txt
python3 methods/gse68465_cldn4_cd274/analyze.py
```

Downloads (not committed) go to `$GSE68465_CLDN4_CD274_DATA` (default `/tmp/gse68465_cldn4_cd274`).
"""
    (HERE / "FINDING.md").write_text(md)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    est_tab = pd.read_csv(SIG, sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()
    if len(stromal) != 141 or len(immune) != 141:
        raise SystemExit(f"Yoshihara lists expected 141+141, got {len(stromal)}+{len(immune)}")

    matrix_path = download(GSE_MATRIX, DATA / "GSE68465_series_matrix.txt.gz")
    annot_path = download(GPL_ANNOT, DATA / "GPL96.annot.gz")

    print("parse matrix", flush=True)
    meta, probes = parse_geo_matrix(matrix_path)
    with gzip.open(annot_path, "rt", errors="replace") as fh:
        annot = parse_soft_table(fh.read())
    annot["symbol"] = annot["Gene symbol"].map(first_symbol)
    first_map = annot.drop_duplicates("ID").set_index("ID")["symbol"]
    first_map = first_map[first_map.astype(str).str.len() > 0]

    if "disease_state" not in meta.columns:
        raise SystemExit(f"GSE68465 disease_state missing: {list(meta.columns)}")
    disease = meta["disease_state"].astype(str)
    is_tumor = disease.str.contains("Adenocarcinoma", case=False, na=False)
    n_matrix = int(len(meta))
    n_tumor = int(is_tumor.sum())
    n_normal = int((~is_tumor).sum())
    print(f"arrays={n_matrix} tumors={n_tumor} held_out={n_normal}", flush=True)
    print(disease.value_counts().to_string(), flush=True)
    if n_tumor != 443:
        print(f"WARNING: tumor n={n_tumor} (sibling n=443)", flush=True)

    probes_t = probes.loc[:, meta.index[is_tumor]]
    genes, probe_pick = collapse_maxmean(probes_t, first_map)
    print(f"genes after max-mean collapse: {genes.shape[0]}", flush=True)

    # Probe confirmation: exact first-symbol + Entrez (no substring aliases).
    confirm_rows = []
    resolved = {}
    probe_used = {}
    probe_rule = {}
    for gene, aliases in ALIASES.items():
        hits = search_annot_exact(annot, aliases, ENTREZ.get(gene))
        chosen = str(probe_pick.loc[gene]) if gene in probe_pick.index else ""
        if gene in genes.index:
            resolved[gene] = gene
            probe_used[gene] = chosen
            probe_rule[gene] = "first-symbol max-mean"
        else:
            resolved[gene] = None
            if gene == "CD274":
                probe_used[gene] = "ABSENT"
                probe_rule[gene] = "no GPL96 hit for CD274 / PDCD1LG1 / Entrez 29126"
        if hits.empty:
            confirm_rows.append({
                "gene": gene,
                "named_probe": NAMED_PROBES.get(gene, ""),
                "probe": "",
                "gpl96_symbol": "",
                "entrez": ENTREZ.get(gene, ""),
                "title": "",
                "in_matrix": False,
                "is_named": False,
                "chosen_for_collapse": False,
                "note": "no GPL96 hit for locked aliases",
            })
            continue
        for _, r in hits.iterrows():
            confirm_rows.append({
                "gene": gene,
                "named_probe": NAMED_PROBES.get(gene, ""),
                "probe": r["ID"],
                "gpl96_symbol": r.get("Gene symbol", ""),
                "entrez": r.get("Gene ID", ""),
                "title": r.get("Gene title", ""),
                "in_matrix": r["ID"] in probes.index,
                "is_named": r["ID"] == NAMED_PROBES.get(gene, ""),
                "chosen_for_collapse": chosen != "" and r["ID"] == chosen,
                "note": "",
            })
    pdl2 = search_annot_exact(annot, ["PDCD1LG2"], ENTREZ["PDCD1LG2"])
    if pdl2.empty:
        confirm_rows.append({
            "gene": "PDCD1LG2",
            "named_probe": "",
            "probe": "",
            "gpl96_symbol": "",
            "entrez": ENTREZ["PDCD1LG2"],
            "title": "",
            "in_matrix": False,
            "is_named": False,
            "chosen_for_collapse": False,
            "note": "PD-L2; related ligand, not CD274; not used as a substitute",
        })
    else:
        for _, r in pdl2.iterrows():
            confirm_rows.append({
                "gene": "PDCD1LG2",
                "named_probe": "",
                "probe": r["ID"],
                "gpl96_symbol": r.get("Gene symbol", ""),
                "entrez": r.get("Gene ID", ""),
                "title": r.get("Gene title", ""),
                "in_matrix": r["ID"] in probes.index,
                "is_named": False,
                "chosen_for_collapse": False,
                "note": "PD-L2; related ligand, not CD274; not used as a substitute",
            })
    confirm = pd.DataFrame(confirm_rows)
    confirm.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    if resolved["CLDN4"] is None:
        raise SystemExit("CLDN4 missing after collapse")
    for g in ("HLA-A", "HLA-B", "HLA-C"):
        if resolved[g] is None:
            raise SystemExit(f"{g} missing after collapse")
        print(f"{g}: {probe_used[g]} ({probe_rule[g]})", flush=True)
    print(f"CD274: ABSENT ({probe_rule.get('CD274', 'no GPL96 hit')})", flush=True)

    est_strom = ssgsea(genes, stromal)
    est_imm = ssgsea(genes, immune)
    estscore_s = est_strom + est_imm
    pur_s = pd.Series(np.cos(EST_A + EST_B * estscore_s.to_numpy(float)), index=genes.columns)
    strom_n = int(sum(g in genes.index for g in stromal))
    imm_n = int(sum(g in genes.index for g in immune))
    print(f"ESTIMATE coverage stromal={strom_n}/141 immune={imm_n}/141", flush=True)

    cldn4 = genes.loc["CLDN4"].to_numpy(float)
    hla_a = genes.loc["HLA-A"].to_numpy(float)
    hla_b = genes.loc["HLA-B"].to_numpy(float)
    hla_c = genes.loc["HLA-C"].to_numpy(float)
    imm = est_imm.to_numpy(float)
    strom = est_strom.to_numpy(float)
    estscore = estscore_s.to_numpy(float)
    pur = pur_s.to_numpy(float)
    in01 = np.isfinite(pur) & (pur >= 0) & (pur <= 1)

    hla_z = []
    for name, v in (("HLA-A", hla_a), ("HLA-B", hla_b), ("HLA-C", hla_c)):
        s = pd.Series(v, dtype=float)
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        hla_z.append(((s - s.mean()) / sd).to_numpy(float))
    mhc = np.mean(np.vstack(hla_z), axis=0) if hla_z else np.full(n_tumor, np.nan)

    est_name = "ESTIMATEScore (Stromal+Immune)"
    pur_name = "ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141)"
    note_eq = "rank-equivalent to TumorPurity residual when cosine is monotone"

    rows = []
    endpoints = [
        ("CD274", None, False),
        ("HLA-A", hla_a, True),
        ("HLA-B", hla_b, True),
        ("HLA-C", hla_c, True),
        ("MHC1_HLAabc", mhc, True),
    ]
    for end, y, present in endpoints:
        note = ""
        if end == "CD274":
            note = "ABSENT on GPL96 U133A; Plus-2 223834_at / 227458_at not on this chip; PDCD1LG2 not substituted"
        elif end == "HLA-B":
            note = "HLA-B is in Yoshihara Immune141; ESTIMATE residual is partly circular"
        rows.append(rec_row("CLDN4", end, cldn4, y, estscore, est_name, note=note or note_eq, present=present))
        rows.append(rec_row("CLDN4", end, cldn4, y, pur, pur_name, note="sensitivity; TumorPurity cosine", present=present))

    rows.append(rec_row(
        "CLDN4", "ImmuneScore", cldn4, imm, estscore, est_name,
        note="context only; ImmuneScore is an addend of ESTIMATEScore; not a CD8 re-cut",
    ))
    rows.append(rec_row("CLDN4", "StromalScore", cldn4, strom, None, ""))
    rows.append(rec_row("CLDN4", "TumorPurity", cldn4, pur, None, ""))
    rows.append(rec_row("CLDN4", "ESTIMATEScore", cldn4, estscore, None, ""))
    rows.append(rec_row("MHC1_HLAabc", "ImmuneScore", mhc, imm, None, "", note="context"))

    if int(in01.sum()) >= 8:
        for end, y, present in (("HLA-A", hla_a, True), ("HLA-B", hla_b, True), ("HLA-C", hla_c, True), ("MHC1_HLAabc", mhc, True)):
            rows.append(rec_row(
                "CLDN4", end, cldn4[in01], y[in01], estscore[in01],
                "ESTIMATEScore; TumorPurity in [0,1]",
                note=f"sensitivity; n={int(in01.sum())} of {n_tumor}",
                present=present,
            ))

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    hl_rows = []
    for end, y, present in endpoints:
        rec = {"endpoint": end}
        rec.update(mwu_q4q1(cldn4, y, present))
        hl_rows.append(rec)
    hl_df = pd.DataFrame(hl_rows)
    hl_df.to_csv(TABLES / "highlow_cldn4_cd274_hla.tsv", sep="\t", index=False)

    per = pd.DataFrame({
        "sample": probes_t.columns,
        "disease_state": disease.loc[probes_t.columns].to_numpy(),
        "CLDN4": cldn4,
        "CD274": np.nan,
        "HLA-A": hla_a,
        "HLA-B": hla_b,
        "HLA-C": hla_c,
        "MHC1_HLAabc_meanz": mhc,
        "ESTIMATE_StromalScore": strom,
        "ESTIMATE_ImmuneScore": imm,
        "ESTIMATE_Score": estscore,
        "ESTIMATE_TumorPurity": pur,
        "purity_in_0_1": in01,
    })
    for col in ("sex", "age", "disease_stage", "histologic_grade", "smoking_history"):
        if col in meta.columns:
            per[col] = meta.loc[probes_t.columns, col].to_numpy()
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    cov = {
        "accession": "GSE68465",
        "title": "Director's Challenge LUAD (Jacob / Shedden multi-site)",
        "platform": "GPL96 U133A",
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_normal": n_normal,
        "tumor_rule": "disease_state contains Adenocarcinoma",
        "n_genes_firstsymbol": int(genes.shape[0]),
        "targets_present": "CLDN4,HLA-A,HLA-B,HLA-C",
        "cd274_present": False,
        "estimate_stromal_present": strom_n,
        "estimate_immune_present": imm_n,
        "mhc1_genes_used": "HLA-A,HLA-B,HLA-C",
        "n_purity_in_0_1": int(in01.sum()),
        "n_purity_outside_0_1": int((~in01).sum()),
        "partial_covariate": est_name,
        "cd8_cut": "taken_as_given_PR298",
        "pdl2_substituted": False,
    }
    pd.DataFrame([cov]).to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    def g(end, col, covar=est_name):
        hit = stats_df[(stats_df.predictor == "CLDN4") & (stats_df.endpoint == end) & (stats_df.covariate == covar)]
        if hit.empty:
            hit = stats_df[(stats_df.predictor == "CLDN4") & (stats_df.endpoint == end)]
        if hit.empty:
            return np.nan
        return hit.iloc[0][col]

    one = pd.DataFrame([{
        "dataset": "GSE68465 LUAD tumors",
        "histology": "LUAD",
        "platform": "GPL96 U133A",
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_normal": n_normal,
        "n_CLDN4_CD274": 0,
        "n_CLDN4_HLAABC": n_tumor,
        "CLDN4_probe": probe_used.get("CLDN4", ""),
        "CD274_probe": "ABSENT",
        "CD274_present": False,
        "HLA_A_probe": probe_used.get("HLA-A", ""),
        "HLA_B_probe": probe_used.get("HLA-B", ""),
        "HLA_C_probe": probe_used.get("HLA-C", ""),
        "CLDN4_vs_CD274_rho": np.nan,
        "CLDN4_vs_CD274_p": np.nan,
        "CLDN4_vs_HLA_A_rho": g("HLA-A", "unadj_rho"),
        "CLDN4_vs_HLA_A_p": g("HLA-A", "unadj_p"),
        "CLDN4_vs_HLA_B_rho": g("HLA-B", "unadj_rho"),
        "CLDN4_vs_HLA_B_p": g("HLA-B", "unadj_p"),
        "CLDN4_vs_HLA_C_rho": g("HLA-C", "unadj_rho"),
        "CLDN4_vs_HLA_C_p": g("HLA-C", "unadj_p"),
        "CLDN4_vs_MHC1_HLAabc_rho": g("MHC1_HLAabc", "unadj_rho"),
        "CLDN4_vs_MHC1_HLAabc_p": g("MHC1_HLAabc", "unadj_p"),
        "CLDN4_vs_CD274_partial_rho": np.nan,
        "CLDN4_vs_CD274_partial_p": np.nan,
        "CLDN4_vs_HLA_A_partial_rho": g("HLA-A", "partial_rho"),
        "CLDN4_vs_HLA_A_partial_p": g("HLA-A", "partial_p"),
        "CLDN4_vs_HLA_B_partial_rho": g("HLA-B", "partial_rho"),
        "CLDN4_vs_HLA_B_partial_p": g("HLA-B", "partial_p"),
        "CLDN4_vs_HLA_C_partial_rho": g("HLA-C", "partial_rho"),
        "CLDN4_vs_HLA_C_partial_p": g("HLA-C", "partial_p"),
        "CLDN4_vs_MHC1_HLAabc_partial_rho": g("MHC1_HLAabc", "partial_rho"),
        "CLDN4_vs_MHC1_HLAabc_partial_p": g("MHC1_HLAabc", "partial_p"),
        "partial_covariate": est_name,
        "ICI_labels": "not_deposited",
        "cd8_recut": "no",
        "cd8_given_pr298_partial_rho": CD8_GIVEN["partial_rho"],
        "cd8_given_pr298_partial_p": CD8_GIVEN["partial_p"],
        "note_cd274": "ABSENT on GPL96 U133A; PDCD1LG2 not substituted",
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "task": "gse68465_cldn4_cd274",
        "additive": True,
        "new_cut_only": True,
        "cd8_taken_as_given": CD8_GIVEN,
        "coverage": cov,
        "probe_used": probe_used,
        "probe_rule": probe_rule,
        "one_row": one.iloc[0].to_dict(),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "bootstrap": {"n": N_BOOT, "seed": BOOT_SEED},
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    write_finding(stats_df, cov, probe_used, probe_rule, n_tumor, n_matrix, n_normal, hl_df)

    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.axis("off")
    ax.text(
        0.5, 0.5,
        "GSE68465 / GPL96 U133A\nCD274 (PD-L1, Entrez 29126) is ABSENT\n"
        "no official first-symbol or Entrez probe\n"
        f"Plus-2 223834_at / 227458_at not on this chip\n"
        f"PDCD1LG2 (PD-L2) present and not substituted\n"
        f"HLA / MHC-I tests use n={n_tumor} LUAD tumors",
        ha="center", va="center", fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.png", bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.pdf", bbox_inches="tight")
    plt.close(fig)

    def scatter(ax, x, y, xlab, ylab, rho, p, n, color):
        ax.scatter(x, y, s=10, alpha=0.45, c=color, linewidths=0)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_rho(rho)}  p={fmt_p(p)}  n={n}", fontsize=10)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    for ax, gene, vec, color in zip(
        axes,
        ("HLA-A", "HLA-B", "HLA-C"),
        (hla_a, hla_b, hla_c),
        ("#E45756", "#54A24B", "#F58518"),
    ):
        scatter(
            ax, cldn4, vec, "CLDN4 (max-mean)", f"{gene} (max-mean)",
            g(gene, "unadj_rho"), g(gene, "unadj_p"), n_tumor, color,
        )
    fig.suptitle(f"GSE68465 LUAD tumors n={n_tumor} — new cut (not CD8)", y=1.03)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_hla.png", bbox_inches="tight")
    fig.savefig(FIGURES / "fig2_cldn4_vs_hla.pdf", bbox_inches="tight")
    plt.close(fig)

    forest = [
        ("CLDN4 vs CD274", np.nan, "#4C78A8", False),
        ("CLDN4 vs HLA-A", g("HLA-A", "unadj_rho"), "#E45756", True),
        ("CLDN4 vs HLA-B", g("HLA-B", "unadj_rho"), "#54A24B", True),
        ("CLDN4 vs HLA-C", g("HLA-C", "unadj_rho"), "#F58518", True),
        ("CLDN4 vs MHC-I mean-z", g("MHC1_HLAabc", "unadj_rho"), "#B279A2", True),
        ("CLDN4 vs HLA-A | ESTIMATE", g("HLA-A", "partial_rho"), "#f4a582", True),
        ("CLDN4 vs HLA-B | ESTIMATE", g("HLA-B", "partial_rho"), "#a1d99b", True),
        ("CLDN4 vs HLA-C | ESTIMATE", g("HLA-C", "partial_rho"), "#fdae6b", True),
        ("CLDN4 vs MHC-I | ESTIMATE", g("MHC1_HLAabc", "partial_rho"), "#d4b9da", True),
    ]
    fig, ax = plt.subplots(figsize=(7.8, 5.4))
    y = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=1)
    for i, (lab, rho, c, present) in enumerate(forest):
        if not present or not np.isfinite(rho):
            ax.text(0.02, i, "ABSENT", color="#922b21", va="center", fontsize=8)
            continue
        ax.plot(rho, i, "o", color=c, ms=7)
        ax.plot([0, rho], [i, i], color=c, lw=2)
    ax.set_yticks(y, [r[0] for r in forest], fontsize=9)
    ax.set_xlabel(f"Spearman ρ (n={n_tumor}; partial = rank residual on ESTIMATEScore)")
    ax.set_title("GSE68465 LUAD tumors — new CD274/HLA cut (CD274 absent on U133A)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_spearman_forest.png", bbox_inches="tight")
    fig.savefig(FIGURES / "fig3_spearman_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    mat = pd.DataFrame({
        "CLDN4": cldn4,
        "HLA-A": hla_a,
        "HLA-B": hla_b,
        "HLA-C": hla_c,
        "MHC1_HLAabc": mhc,
        "ImmuneScore": imm,
        "ESTIMATEScore": estscore,
    })
    cols = list(mat.columns)
    R = np.zeros((len(cols), len(cols)))
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            r, _p, _n = spearman_pair(mat[a].to_numpy(float), mat[b].to_numpy(float))
            R[i, j] = r
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(R, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right")
    ax.set_yticks(range(len(cols)), cols)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, label="Spearman ρ")
    ax.set_title(f"GSE68465 LUAD tumors n={n_tumor} (CD274 absent)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_correlation_heatmap.png", bbox_inches="tight")
    fig.savefig(FIGURES / "fig4_correlation_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    print(one.T.to_string())
    print(stats_df[["predictor", "endpoint", "n_unadjusted", "unadj_rho", "unadj_p", "partial_rho", "partial_p", "verdict_unadj"]].to_string(index=False))
    print("wrote", TABLES)
    print("wrote", FIGURES)
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
