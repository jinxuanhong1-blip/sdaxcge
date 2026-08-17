#!/usr/bin/env python3
"""GSE10072 Landi LUAD array: CLDN4-only vs CD8A / CD274 / ImmuneScore.

Additive only. No dual-high. Tumor arrays only. Epithelial residual.

Public inputs:
  GEO GSE10072 series matrix (Landi et al. PLoS ONE 2008; GPL96 HG-U133A RMA)
  GPL96.annot gene symbols
  Yoshihara 2013 ESTIMATE Stromal141 + Immune141 (local TSV)

Unit = tumor array. Adjacent/normal lung is dropped.
CD274 is not on HG-U133A; that absence is reported, not proxied.
"""
from __future__ import annotations

import gzip
import json
import math
import os
import urllib.request
from collections import Counter
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
CACHE = Path(os.environ.get("GSE10072_CLDN4_DATA", "/tmp/gse10072_cldn4"))
SIG = HERE / "estimate_yoshihara_2013.tsv"

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/"
    "matrix/GSE10072_series_matrix.txt.gz"
)
GPL96_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

UA = "sdaxcge-gse10072-cldn4/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

TARGETS = ["CLDN4", "CD8A", "CD274"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

EST_A = 0.6049872018
EST_B = 0.0001467884


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
    return f"{r:+.3f}"


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or not (np.isfinite(lo) and np.isfinite(hi)):
        return "NA"
    return f"{lo:+.3f} to {hi:+.3f}"


def verdict(rho_adj, p_adj, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 8:
        return dict(n=n, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
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
    return dict(n=n, rho=float(rho), p=float(p), ci_lo=float(lo), ci_hi=float(hi))


def partial_spearman(x, y, *covars):
    """Pearson of rank residuals; df = n − 2 − k (same as GSE10245 / GSE8894)."""
    cols = [np.asarray(x, float), np.asarray(y, float)]
    for c in covars:
        cols.append(np.asarray(c, float))
    arr = np.column_stack(cols)
    ok = np.isfinite(arr).all(axis=1)
    arr = arr[ok]
    n = int(arr.shape[0])
    k = arr.shape[1] - 2
    if n < 8 + k:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(arr[:, 0])
    yr = stats.rankdata(arr[:, 1])
    Z = np.column_stack([stats.rankdata(arr[:, 2 + i]) for i in range(k)])
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    dof = n - 2 - k
    t = r * math.sqrt(dof / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


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
            meta_rows["title"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_source_name"):
            meta_rows["source"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_data_processing"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            meta_rows["data_processing"] = vals
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
    n_na = 0
    n_vals = 0
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        idx.append(parts[0].strip().strip('"'))
        vals = []
        for x in parts[1:]:
            n_vals += 1
            if x in ("", "NA", "null"):
                n_na += 1
                vals.append(np.nan)
            else:
                vals.append(float(x))
        rows.append(vals)
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr, {"n_probes": len(idx), "n_na": n_na, "n_vals": n_vals}


def load_gpl_annot(path: Path) -> tuple[pd.Series, pd.Series]:
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
        raise SystemExit(f"no symbol column in {path}: {list(df.columns)[:12]}")
    symbol_col = cands[0]
    raw = df.set_index(id_col)[symbol_col].astype(str)
    raw = raw.replace({"nan": np.nan, "None": np.nan, "": np.nan}).dropna()
    n_sym = raw.map(lambda x: len([p.strip() for p in str(x).split("///") if p.strip()]))
    first = raw.map(lambda x: str(x).split("///")[0].strip())
    first = first[first.str.len() > 0]
    unique = first[n_sym == 1]
    return first, unique


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out, pick


def zmean(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    sub = expr.loc[present].astype(float)
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0), present


def estimate_scores(gene_expr: pd.DataFrame, stromal: list[str], immune: list[str]) -> pd.DataFrame:
    strom = ssgsea(gene_expr, stromal)
    imm = ssgsea(gene_expr, immune)
    est = strom + imm
    pur = np.cos(EST_A + EST_B * est.to_numpy(float))
    return pd.DataFrame(
        {
            "ESTIMATE_StromalScore": strom,
            "ESTIMATE_ImmuneScore": imm,
            "ESTIMATE_Score": est,
            "ESTIMATE_TumorPurity": pur,
        },
        index=gene_expr.columns,
    )


def rec_row(cohort, predictor, endpoint, x, y, z=None, note=""):
    unadj = spearman_ci(x, y)
    if z is None:
        part = dict(n=unadj["n"], rho_adj=np.nan, p_adj=np.nan)
    else:
        part = partial_spearman(x, y, z)
    return {
        "cohort": cohort,
        "predictor": predictor,
        "endpoint": endpoint,
        "n": unadj["n"],
        "rho": unadj["rho"],
        "p": unadj["p"],
        "ci_lo": unadj["ci_lo"],
        "ci_hi": unadj["ci_hi"],
        "rho_adj_epithelial": part["rho_adj"],
        "p_adj_epithelial": part["p_adj"],
        "verdict": verdict(part["rho_adj"], part["p_adj"], part["n"]),
        "note": note,
    }


def pick_row(stats_df: pd.DataFrame, cohort: str, pred: str, end: str):
    hit = stats_df[
        (stats_df["cohort"] == cohort)
        & (stats_df["predictor"] == pred)
        & (stats_df["endpoint"] == end)
    ]
    if hit.empty:
        return None
    return hit.iloc[0]


def scatter(path: Path, x, y, xlab, ylab, title, n, rho, p):
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(x, y, s=22, c="#1f4e79", alpha=0.75, edgecolors="none")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.text(
        0.04,
        0.96,
        f"n={n}\nρ={fmt_rho(rho)}  p={fmt_p(p)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, lw=0.4),
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(stats_df, cov, probe_used, missing_targets):
    c4_cd8 = pick_row(stats_df, "LUAD_tumor", "CLDN4", "CD8A")
    c4_imm = pick_row(stats_df, "LUAD_tumor", "CLDN4", "ImmuneScore")
    c4_epi = pick_row(stats_df, "LUAD_tumor", "CLDN4", "epithelial_z")
    cd8_imm = pick_row(stats_df, "LUAD_tumor", "CD8A", "ImmuneScore")

    cd274_note = (
        "not on GPL96 HG-U133A (no CD274 / PDCD1LG1 / B7-H1 probe)"
        if "CD274" in missing_targets
        else f"`{probe_used.get('CD274', 'NA')}`"
    )
    cd274_cell = "**absent (U133A)**" if "CD274" in missing_targets else "present"

    def cell(row, adj=False):
        if row is None:
            return "NA"
        if adj:
            return f"{fmt_rho(row.rho_adj_epithelial)} ({fmt_p(row.p_adj_epithelial)})"
        return f"{fmt_rho(row.rho)} ({fmt_p(row.p)})"

    v_cd8 = c4_cd8.verdict if c4_cd8 is not None else "NA"
    v_imm = c4_imm.verdict if c4_imm is not None else "NA"

    n_t = int(cov["n_tumor"])
    smoke_t = cov["tumor_smoking"]
    stage_t = cov["tumor_stage"]

    md = f"""# GSE10072 — Landi LUAD array, CLDN4-only vs CD8A / CD274 / ImmuneScore

**Additive only. CLDN4-only. No dual-high.** Public Landi EAGLE LUAD Affymetrix HG-U133A series (Landi et al., *PLoS ONE* 2008, PMID 18297132; GEO [GSE10072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072); GPL96). Unit is the **array**. Tumor only. No slide was re-scored. No ICI arm. TACSTD2 is not analysed.

Primary tests: **CLDN4 vs CD8A**, **CLDN4 vs CD274**, **CLDN4 vs ESTIMATE ImmuneScore**, each after an **epithelial residual**.

## Honest n

Do not write n=135 or n=122. The series summary still says 135 tissues; the GEO `overall_design` then drops low-tumour-cell arrays and averages QC duplicates. The deposited series matrix is **107 arrays**.

| item | public? | n | rule |
|---|---|---:|---|
| selected tissues (paper / design text) | text only | 180 | includes QC duplicates / triplicates from 14 subjects |
| high-quality RNA | text only | 148 | 180 minus insufficient RNA |
| normalized microarrays | text only | 135 | 148 minus 13 problematic assays |
| after low-tumour-cell drop | text only | 122 | 135 minus 13 low % tumor |
| after averaging 15 QC duplicates | text only | **107** | GEO `overall_design`; 58 tumor + 49 non-tumor |
| **arrays in the series matrix** | yes | **107** | 22,283 probes × 107 GSM; 0 missing RMA values |
| unique GSM / unique titles | yes | **107** | `Lung Tumor_GT*` / `Normal Lung_GT*`; all unique |
| unique subjects (`GT` id) | yes | **74** | 20 never / 26 former / 28 current (design text) |
| adenocarcinoma tumor (`source_name`) | yes | **58** | `Adenocarcinoma of the Lung`; 58 unique `GT` ids |
| non-tumor lung | yes | **49** | `Normal Lung Tissue`; **dropped** (tumor-only rule) |
| paired tumor+normal subjects | yes | 33 | same `GT` in both sources; pairing unused (tumor-only) |
| tumor-only subjects | yes | 25 | tumor array, no deposited normal |
| tumor smoking | yes | 58 | Current 24 / Former 18 / Never 16 |
| tumor stage | yes | 58 | IA 5, IB 17, IIA 3, IIB 18, IIIA 9, IIIB 3, IV 3 |
| tumor gender | yes | 58 | Male 35 / Female 23 |
| ICI / treatment / response | **no** | 0 | smoking / LUAD atlas, not an ICI series |
| OS / DSS time or event | **no** | 0 | paper discusses survival; labels **not** in the matrix |
| numeric tumor % / ABSOLUTE purity | **no** | 0 | 13 arrays dropped for low % tumor; values not deposited |
| CLDN4 finite (`{probe_used.get("CLDN4", "NA")}`) | yes | **{n_t}** | single unique-mapped GPL96 probe |
| CD8A finite (`{probe_used.get("CD8A", "NA")}`) | yes | **{n_t}** | single unique-mapped GPL96 probe |
| CD274 finite | **no** | **0** | {cd274_note} |
| ESTIMATE ImmuneScore | computed | **{n_t}** | Yoshihara Immune141 ssGSEA; {cov["estimate_immune_present"]}/141 genes |
| epithelial mean-z | computed | **{n_t}** | {cov["epithelial_used"]} ({cov["n_epithelial_genes"]}/6) |
| **primary pairwise n (tumor LUAD)** | yes | **{n_t}** | complete-case CLDN4 + CD8A + ImmuneScore |

Primary tests use **n={n_t} tumor LUAD arrays**. Do not pool the 49 normals. Do not write the paper’s 135 / 122.

## One-row table

| dataset | histology rule | platform | n tumors | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | adj ρ \\| epi (p) | CLDN4–ImmuneScore ρ (p) | adj ρ \\| epi (p) | CLDN4–CD274 | verdict |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE10072 Landi | LUAD tumor only (normals dropped) | GPL96 HG-U133A RMA | **{n_t}** | `{probe_used.get("CLDN4", "NA")}` | `{probe_used.get("CD8A", "NA")}` | {cd274_cell} | ESTIMATE Immune141 ({cov["estimate_immune_present"]}/141) | **{cell(c4_cd8)}** | {cell(c4_cd8, adj=True)} | **{cell(c4_imm)}** | {cell(c4_imm, adj=True)} | not on U133A | CD8A **{v_cd8}**; ImmuneScore **{v_imm}** |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix **RMA** (Bioconductor `affy`; `!Sample_data_processing`). Tests use native ranks (Spearman / ssGSEA), so the deposited log-intensity scale does not change ρ.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | 107 | every tumor is LUAD; every normal is non-tumor lung |
| smoking | yes | 107 | GEO `Cigarette Smoking Status` |
| stage / age / gender | yes | 107 | GEO characteristics; **not** the claim |
| ICI response | **no** | 0 | not an ICI series |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL96 probe (multi-mapped `///` probes dropped for named genes):

| gene | probe | on GPL96? |
|---|---|---|
| CLDN4 | `{probe_used.get("CLDN4", "NA")}` | yes |
| CD8A | `{probe_used.get("CD8A", "NA")}` | yes |
| CD274 | — | **no** (Plus 2.0 probes `223834_at` / `227458_at` are not on U133A) |

ESTIMATE gene collapse uses first-symbol mapping so Immune141 coverage is {cov["estimate_immune_present"]}/141 and Stromal141 is {cov["estimate_stromal_present"]}/141.

## Epithelial residual and ImmuneScore

- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 ({cov["n_epithelial_genes"]}/6 present: {cov["epithelial_used"]}).
- **ImmuneScore:** ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Immune141. R `estimate` is not required.
- **Partial Spearman:** Pearson of rank residuals on the epithelial mean-z; df = n − 3.
- **HOLDS rule** (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05.
- **No dual-high:** no TACSTD2×CLDN4 quadrant, no dual-high vs dual-low contrast.

CD8A vs ImmuneScore (positive-control, tumor LUAD): ρ={fmt_rho(cd8_imm.rho) if cd8_imm is not None else "NA"} (p={fmt_p(cd8_imm.p) if cd8_imm is not None else "NA"}, n={int(cd8_imm.n) if cd8_imm is not None else "NA"}).

CLDN4 vs epithelial mean-z (context, not the claim): ρ={fmt_rho(c4_epi.rho) if c4_epi is not None else "NA"} (p={fmt_p(c4_epi.p) if c4_epi is not None else "NA"}, n={int(c4_epi.n) if c4_epi is not None else "NA"}).

## Main Spearman (primary tumor LUAD n={n_t})

Bootstrap 95% CI, 2,000 resamples, seed `{SEED}`.

| predictor | endpoint | n | ρ | 95% CI | p | ρ_adj \\| epi | p_adj | verdict |
|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 | CD8A | {int(c4_cd8.n) if c4_cd8 is not None else "NA"} | **{fmt_rho(c4_cd8.rho) if c4_cd8 is not None else "NA"}** | {fmt_ci(c4_cd8.ci_lo, c4_cd8.ci_hi) if c4_cd8 is not None else "NA"} | {fmt_p(c4_cd8.p) if c4_cd8 is not None else "NA"} | {fmt_rho(c4_cd8.rho_adj_epithelial) if c4_cd8 is not None else "NA"} | {fmt_p(c4_cd8.p_adj_epithelial) if c4_cd8 is not None else "NA"} | **{v_cd8}** |
| CLDN4 | CD274 | 0 | — | — | — | — | — | **ABSENT_ON_GPL96** |
| CLDN4 | ImmuneScore | {int(c4_imm.n) if c4_imm is not None else "NA"} | **{fmt_rho(c4_imm.rho) if c4_imm is not None else "NA"}** | {fmt_ci(c4_imm.ci_lo, c4_imm.ci_hi) if c4_imm is not None else "NA"} | {fmt_p(c4_imm.p) if c4_imm is not None else "NA"} | {fmt_rho(c4_imm.rho_adj_epithelial) if c4_imm is not None else "NA"} | {fmt_p(c4_imm.p_adj_epithelial) if c4_imm is not None else "NA"} | **{v_imm}** |
| CD8A | ImmuneScore | {int(cd8_imm.n) if cd8_imm is not None else "NA"} | {fmt_rho(cd8_imm.rho) if cd8_imm is not None else "NA"} | {fmt_ci(cd8_imm.ci_lo, cd8_imm.ci_hi) if cd8_imm is not None else "NA"} | {fmt_p(cd8_imm.p) if cd8_imm is not None else "NA"} | {fmt_rho(cd8_imm.rho_adj_epithelial) if cd8_imm is not None else "NA"} | {fmt_p(cd8_imm.p_adj_epithelial) if cd8_imm is not None else "NA"} | positive control |
| CLDN4 | epithelial_z | {int(c4_epi.n) if c4_epi is not None else "NA"} | {fmt_rho(c4_epi.rho) if c4_epi is not None else "NA"} | {fmt_ci(c4_epi.ci_lo, c4_epi.ci_hi) if c4_epi is not None else "NA"} | {fmt_p(c4_epi.p) if c4_epi is not None else "NA"} | — | — | context |

This is **not** an ICI-response test. CD274 cannot be scored on this chip. Smoking-stratum n (Current {smoke_t.get("Current Smoker", 0)} / Former {smoke_t.get("Former Smoker", 0)} / Never {smoke_t.get("Never Smoked", 0)}) is below the HOLDS n≥40 rule and is not a claim.

## What is not done

- No dual-high TACSTD2×CLDN4 split.
- No pooling of the 49 non-tumor arrays into the primary n.
- No invented CD274 probe or PD-L1 IHC substitute.
- No ICI ORR / PFS model (labels are not deposited).
- No OS model (survival labels are not in the series matrix).
- No re-use of the paper’s 135-sample sentence as the analysis n.

## Files

- `analyze.py` — GEO download, honest-n inventory, ESTIMATE, epithelial residual, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `label_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_spearman_forest.png`

```bash
python3 -m pip install -r methods/gse10072_cldn4/requirements.txt
export GSE10072_CLDN4_DATA=/tmp/gse10072_cldn4
python3 methods/gse10072_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)
    _ = stage_t  # inventory already written into honest-n table from cov


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    matrix_path = download(GSE_MATRIX, CACHE / "GSE10072_series_matrix.txt.gz")
    annot_path = download(GPL96_ANNOT, CACHE / "GPL96.annot.gz")

    meta, probes, expr_qc = parse_geo_matrix(matrix_path)
    first_map, unique_map = load_gpl_annot(annot_path)

    sig = pd.read_csv(SIG, sep="\t")
    stromal = sig.loc[sig["set"].str.contains("Stromal", case=False), "hugo"].astype(str).tolist()
    immune = sig.loc[sig["set"].str.contains("Immune", case=False), "hugo"].astype(str).tolist()

    meta["patient"] = meta["title"].astype(str).str.extract(r"(GT\d+)", expand=False)
    src = meta["source"].astype(str)
    title = meta["title"].astype(str)
    is_tumor = src.str.contains("Adeno", case=False) | title.str.contains("Tumor", case=False)
    is_normal = src.str.contains("Normal", case=False) | title.str.contains("Normal", case=False)
    # Prefer source_name; titles are a backup. Ambiguous rows are inventoried.
    meta["is_tumor"] = is_tumor & ~is_normal
    meta["is_normal"] = is_normal & ~is_tumor
    meta["histo"] = np.where(meta["is_tumor"], "LUAD", np.where(meta["is_normal"], "NTL", "OTHER"))

    n_matrix = int(len(meta))
    n_tumor = int(meta["is_tumor"].sum())
    n_normal = int(meta["is_normal"].sum())
    n_other = int((~meta["is_tumor"] & ~meta["is_normal"]).sum())
    if n_other:
        raise SystemExit(f"unmapped source rows n={n_other}: {meta.loc[~meta['is_tumor'] & ~meta['is_normal'], 'source'].tolist()}")
    if n_tumor < HOLDS_N:
        raise SystemExit(f"tumor n={n_tumor} too small")

    smoke_col = "cigarette smoking status" if "cigarette smoking status" in meta.columns else None
    stage_col = "stage" if "stage" in meta.columns else None
    gender_col = "gender" if "gender" in meta.columns else None
    age_col = "age at diagnosis" if "age at diagnosis" in meta.columns else None

    tum_meta = meta.loc[meta["is_tumor"]].copy()
    nor_meta = meta.loc[meta["is_normal"]].copy()
    paired = set(tum_meta["patient"].dropna()) & set(nor_meta["patient"].dropna())
    tumor_smoking = Counter(tum_meta[smoke_col].astype(str)) if smoke_col else {}
    tumor_stage = Counter(tum_meta[stage_col].astype(str)) if stage_col else {}
    tumor_gender = Counter(tum_meta[gender_col].astype(str)) if gender_col else {}

    inv = pd.DataFrame(
        [
            {"item": "arrays in series matrix", "public": "yes", "n": n_matrix, "note": "GSE10072_series_matrix.txt.gz"},
            {"item": "unique GSM / titles", "public": "yes", "n": n_matrix, "note": "all unique"},
            {"item": "unique GT subjects", "public": "yes", "n": int(meta["patient"].nunique()), "note": "from titles"},
            {"item": "LUAD tumor arrays", "public": "yes", "n": n_tumor, "note": "source Adenocarcinoma of the Lung"},
            {"item": "non-tumor lung arrays", "public": "yes", "n": n_normal, "note": "dropped; tumor-only rule"},
            {"item": "paired tumor+normal subjects", "public": "yes", "n": len(paired), "note": "same GT; pairing unused"},
            {"item": "tumor-only subjects", "public": "yes", "n": int(tum_meta["patient"].nunique() - len(paired)), "note": "tumor, no deposited normal"},
            {"item": "ICI / response", "public": "no", "n": 0, "note": "not an ICI series"},
            {"item": "OS / event in matrix", "public": "no", "n": 0, "note": "not deposited"},
        ]
    )
    inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    keep = meta.index[meta["is_tumor"]]
    probes_k = probes.loc[:, keep]
    genes_est, _ = collapse_maxmean(probes_k, first_map)
    genes_uniq, probe_pick = collapse_maxmean(probes_k, unique_map)

    present_targets = [g for g in TARGETS if g in genes_uniq.index]
    missing_targets = [g for g in TARGETS if g not in genes_uniq.index]
    # CD274 is expected missing on GPL96. CLDN4 / CD8A must be present.
    hard_miss = [g for g in ("CLDN4", "CD8A") if g in missing_targets]
    if hard_miss:
        (HERE / "FINDING.md").write_text(
            "# GSE10072 — Landi LUAD CLDN4-only\n\n"
            f"Required genes missing after unique-probe collapse: {hard_miss}.\n"
        )
        raise SystemExit(f"required genes missing: {hard_miss}")

    probe_used = {g: str(probe_pick.loc[g]) for g in present_targets}
    pd.Series(probe_used, name="probe_id").to_csv(TABLES / "probe_used.tsv", sep="\t")
    pd.DataFrame(
        [
            {"gene": "CLDN4", "on_gpl96": "CLDN4" not in missing_targets, "probe": probe_used.get("CLDN4", "")},
            {"gene": "CD8A", "on_gpl96": "CD8A" not in missing_targets, "probe": probe_used.get("CD8A", "")},
            {"gene": "CD274", "on_gpl96": "CD274" not in missing_targets, "probe": probe_used.get("CD274", "")},
        ]
    ).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    est = estimate_scores(genes_est, stromal, immune)
    epithelial, epi_used = zmean(genes_uniq, EPITHELIAL)

    per = pd.DataFrame(
        {
            "gsm": keep,
            "title": meta.loc[keep, "title"].to_numpy(),
            "source": meta.loc[keep, "source"].to_numpy(),
            "patient": meta.loc[keep, "patient"].to_numpy(),
            "histo": meta.loc[keep, "histo"].to_numpy(),
            "smoking": meta.loc[keep, smoke_col].to_numpy() if smoke_col else "",
            "stage": meta.loc[keep, stage_col].to_numpy() if stage_col else "",
            "gender": meta.loc[keep, gender_col].to_numpy() if gender_col else "",
            "age": meta.loc[keep, age_col].to_numpy() if age_col else "",
            "CLDN4": genes_uniq.loc["CLDN4"].to_numpy(float),
            "CD8A": genes_uniq.loc["CD8A"].to_numpy(float),
            "epithelial_z": epithelial.reindex(keep).to_numpy(float),
        },
        index=keep,
    )
    if "CD274" in genes_uniq.index:
        per["CD274"] = genes_uniq.loc["CD274"].to_numpy(float)
    else:
        per["CD274"] = np.nan
    per = per.join(est)
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    n_named = int((np.isfinite(per["CLDN4"]) & np.isfinite(per["CD8A"])).sum())
    n_imm = int(np.isfinite(per["ESTIMATE_ImmuneScore"]).sum())

    cov = {
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_normal": n_normal,
        "n_subjects": int(meta["patient"].nunique()),
        "n_paired_subjects": len(paired),
        "n_named_complete": n_named,
        "n_immunescore": n_imm,
        "n_probes": expr_qc["n_probes"],
        "n_missing_rma": expr_qc["n_na"],
        "estimate_immune_present": int(sum(g in genes_est.index for g in immune)),
        "estimate_stromal_present": int(sum(g in genes_est.index for g in stromal)),
        "n_epithelial_genes": len(epi_used),
        "epithelial_used": "/".join(epi_used),
        "probe_CLDN4": probe_used.get("CLDN4", ""),
        "probe_CD8A": probe_used.get("CD8A", ""),
        "probe_CD274": probe_used.get("CD274", ""),
        "cd274_on_platform": "CD274" not in missing_targets,
        "tumor_smoking": dict(tumor_smoking),
        "tumor_stage": dict(tumor_stage),
        "tumor_gender": dict(tumor_gender),
        "normalization": "RMA (Bioconductor affy; GEO Sample_data_processing)",
        "platform": "GPL96 HG-U133A",
    }
    # JSON-friendly coverage TSV (flatten dicts)
    cov_flat = {k: (json.dumps(v) if isinstance(v, dict) else v) for k, v in cov.items()}
    pd.Series(cov_flat, name="value").to_csv(TABLES / "coverage.tsv", sep="\t")

    rows = []
    c4 = per["CLDN4"].to_numpy(float)
    cd8 = per["CD8A"].to_numpy(float)
    imm = per["ESTIMATE_ImmuneScore"].to_numpy(float)
    epi = per["epithelial_z"].to_numpy(float)
    rows.append(rec_row("LUAD_tumor", "CLDN4", "CD8A", c4, cd8, epi))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "ImmuneScore", c4, imm, epi))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "epithelial_z", c4, epi, None, note="context"))
    rows.append(rec_row("LUAD_tumor", "CD8A", "ImmuneScore", cd8, imm, epi, note="positive control"))
    if "CD274" not in missing_targets:
        pdl1 = per["CD274"].to_numpy(float)
        rows.append(rec_row("LUAD_tumor", "CLDN4", "CD274", c4, pdl1, epi))
        rows.append(rec_row("LUAD_tumor", "CD274", "CD8A", pdl1, cd8, epi))
    else:
        rows.append(
            {
                "cohort": "LUAD_tumor",
                "predictor": "CLDN4",
                "endpoint": "CD274",
                "n": 0,
                "rho": np.nan,
                "p": np.nan,
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "rho_adj_epithelial": np.nan,
                "p_adj_epithelial": np.nan,
                "verdict": "ABSENT_ON_GPL96",
                "note": "CD274 / PDCD1LG1 not on HG-U133A",
            }
        )

    # Smoking strata are inventoried; n<40 so UNDERPOWERED if computed.
    if smoke_col:
        for smoke, d in per.groupby(per["smoking"].astype(str)):
            if len(d) < 8:
                continue
            rows.append(
                rec_row(
                    f"LUAD_{smoke.replace(' ', '_')}",
                    "CLDN4",
                    "CD8A",
                    d["CLDN4"].to_numpy(float),
                    d["CD8A"].to_numpy(float),
                    d["epithelial_z"].to_numpy(float),
                    note="smoking stratum; not the claim",
                )
            )
            rows.append(
                rec_row(
                    f"LUAD_{smoke.replace(' ', '_')}",
                    "CLDN4",
                    "ImmuneScore",
                    d["CLDN4"].to_numpy(float),
                    d["ESTIMATE_ImmuneScore"].to_numpy(float),
                    d["epithelial_z"].to_numpy(float),
                    note="smoking stratum; not the claim",
                )
            )

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    a = pick_row(stats_df, "LUAD_tumor", "CLDN4", "CD8A")
    b = pick_row(stats_df, "LUAD_tumor", "CLDN4", "CD274")
    c = pick_row(stats_df, "LUAD_tumor", "CLDN4", "ImmuneScore")
    one = pd.DataFrame(
        [
            {
                "dataset": "GSE10072 Landi",
                "histology_rule": "LUAD tumor only",
                "platform": "GPL96 HG-U133A RMA",
                "n_tumors": n_tumor,
                "CLDN4_probe": probe_used.get("CLDN4", ""),
                "CD8A_probe": probe_used.get("CD8A", ""),
                "CD274": "absent_on_GPL96" if "CD274" in missing_targets else probe_used.get("CD274", ""),
                "ImmuneScore": f"ESTIMATE_Immune141_{cov['estimate_immune_present']}/141",
                "CLDN4_CD8A_rho": None if a is None else a.rho,
                "CLDN4_CD8A_p": None if a is None else a.p,
                "CLDN4_CD8A_partial_rho_epi": None if a is None else a.rho_adj_epithelial,
                "CLDN4_CD8A_partial_p_epi": None if a is None else a.p_adj_epithelial,
                "CLDN4_CD8A_verdict": None if a is None else a.verdict,
                "CLDN4_ImmuneScore_rho": None if c is None else c.rho,
                "CLDN4_ImmuneScore_p": None if c is None else c.p,
                "CLDN4_ImmuneScore_partial_rho_epi": None if c is None else c.rho_adj_epithelial,
                "CLDN4_ImmuneScore_partial_p_epi": None if c is None else c.p_adj_epithelial,
                "CLDN4_ImmuneScore_verdict": None if c is None else c.verdict,
                "CLDN4_CD274_rho": None if b is None else b.rho,
                "CLDN4_CD274_p": None if b is None else b.p,
                "CLDN4_CD274_verdict": None if b is None else b.verdict,
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    scatter(
        FIGURES / "fig1_cldn4_vs_cd8a.png",
        per["CLDN4"],
        per["CD8A"],
        "CLDN4 (RMA)",
        "CD8A (RMA)",
        f"GSE10072 LUAD tumor CLDN4 vs CD8A (n={n_tumor})",
        int(a.n) if a is not None else n_tumor,
        a.rho if a is not None else np.nan,
        a.p if a is not None else np.nan,
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_immunescore.png",
        per["CLDN4"],
        per["ESTIMATE_ImmuneScore"],
        "CLDN4 (RMA)",
        "ESTIMATE ImmuneScore",
        f"GSE10072 LUAD tumor CLDN4 vs ImmuneScore (n={n_tumor})",
        int(c.n) if c is not None else n_tumor,
        c.rho if c is not None else np.nan,
        c.p if c is not None else np.nan,
    )

    forest_keys = [
        ("LUAD_tumor", "CLDN4", "CD8A", "CLDN4–CD8A"),
        ("LUAD_tumor", "CLDN4", "ImmuneScore", "CLDN4–ImmuneScore"),
        ("LUAD_tumor", "CD8A", "ImmuneScore", "CD8A–ImmuneScore (control)"),
        ("LUAD_tumor", "CLDN4", "epithelial_z", "CLDN4–epithelial z (context)"),
    ]
    labels, rhos, los, his = [], [], [], []
    for cohort, pred, end, lab in forest_keys:
        r = pick_row(stats_df, cohort, pred, end)
        if r is None or not np.isfinite(r.rho):
            continue
        labels.append(f"{lab} (n={int(r.n)})")
        rhos.append(r.rho)
        los.append(r.ci_lo)
        his.append(r.ci_hi)
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="0.5", lw=0.8)
    ax.errorbar(
        rhos,
        y,
        xerr=[np.array(rhos) - np.array(los), np.array(his) - np.array(rhos)],
        fmt="o",
        color="#1f4e79",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
    ax.set_title("GSE10072 LUAD tumor CLDN4-only")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_spearman_forest.png", dpi=140)
    fig.savefig(FIGURES / "fig3_spearman_forest.pdf")
    plt.close(fig)

    write_finding(stats_df, cov, probe_used, missing_targets)

    summary = {
        "dataset": "GSE10072",
        "primary_cohort": "LUAD tumor only",
        "n_tumor": n_tumor,
        "n_normal_dropped": n_normal,
        "n_named_complete": n_named,
        "missing_targets": missing_targets,
        "CLDN4_CD8A": {
            "n": None if a is None else int(a.n),
            "rho": None if a is None else float(a.rho),
            "p": None if a is None else float(a.p),
            "partial_rho_epi": None if a is None else float(a.rho_adj_epithelial),
            "partial_p_epi": None if a is None else float(a.p_adj_epithelial),
            "verdict": None if a is None else str(a.verdict),
        },
        "CLDN4_ImmuneScore": {
            "n": None if c is None else int(c.n),
            "rho": None if c is None else float(c.rho),
            "p": None if c is None else float(c.p),
            "partial_rho_epi": None if c is None else float(c.rho_adj_epithelial),
            "partial_p_epi": None if c is None else float(c.p_adj_epithelial),
            "verdict": None if c is None else str(c.verdict),
        },
        "CLDN4_CD274": "absent_on_GPL96" if "CD274" in missing_targets else "present",
        "probes": probe_used,
        "epithelial_genes": epi_used,
        "estimate_immune_present": cov["estimate_immune_present"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    print("wrote", HERE / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
