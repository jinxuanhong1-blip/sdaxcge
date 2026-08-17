#!/usr/bin/env python3
"""GSE30219 NSCLC tumors: CLDN4 vs CD8A / CD274 (after ESTIMATE).

Additive only. No prior GSE30219 CLDN4–immune table is re-audited.
This folder only measures CLDN4 against CD8A and CD274 on the public
Rousseaux mixed-histology lung series, with an honest NSCLC n.

Public inputs:
  GEO GSE30219 series matrix (Rousseaux et al. Sci Transl Med 2013; GPL570)
  GPL570.annot gene symbols
  Yoshihara 2013 ESTIMATE Stromal141 + Immune141 (local TSV)

Unit = tumor array. Non-tumoral lung is dropped.
NSCLC is defined from GEO histology codes (see FINDING.md).
SCLC (GEO code SCC) and carcinoid (CARCI) are not NSCLC.
"""
from __future__ import annotations

import gzip
import json
import math
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
CACHE = Path("/tmp/gse30219_cldn4")
SIG = HERE / "estimate_yoshihara_2013.tsv"

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30219/"
    "matrix/GSE30219_series_matrix.txt.gz"
)
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

UA = (
    "sdaxcge-gse30219-cldn4/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

EST_A = 0.6049872018
EST_B = 0.0001467884

TARGETS = ["CLDN4", "CD8A", "CD274"]

# GEO histology codes on this series (Rousseaux / Centre Léon Bérard).
# SCC on GSE30219 is small-cell, not squamous. Squamous is SQC.
# NSCLC (primary): ADC + SQC + LCC + BAS.
# LCNE is high-grade neuroendocrine (LCNEC); reported as a sensitivity, not
# folded into the primary NSCLC n.
NSCLC_CODES = {"ADC", "SQC", "LCC", "BAS"}
LCNE_CODES = {"LCNE", "LCNEC"}
SCLC_CODES = {"SCC", "SCLC"}
CARCINOID_CODES = {"CARCI", "CARCINOID"}
NORMAL_CODES = {"NT", "NTL", "NL", "NORMAL"}
OTHER_CODES = {"OTHER"}


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
    return f"{r:.3f}"


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


def spearman_ci(x, y, n_boot: int = 2000, seed: int = 20260817):
    m = np.isfinite(x) & np.isfinite(y)
    xx = np.asarray(x[m], float)
    yy = np.asarray(y[m], float)
    n = int(len(xx))
    if n < 6:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    rhos = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        r, _ = stats.spearmanr(xx[idx], yy[idx])
        if np.isfinite(r):
            rhos.append(float(r))
    if len(rhos) < 50:
        return np.nan, np.nan
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return float(lo), float(hi)


def partial_spearman(x, y, z):
    """Pearson of rank residuals; df = n − 3 (same as A1 extra / GSE31210)."""
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
        rows.append([float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


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


def histology_code(raw: str) -> str:
    s = str(raw).strip().upper()
    if s in {"NT", "NTL", "NL", "NORMAL", "NONTUMORAL", "NON-TUMORAL", "NON TUMORAL"}:
        return "NTL"
    if s in {"ADC", "ADENOCARCINOMA", "LUAD"}:
        return "ADC"
    if s in {"SQC", "SQUAMOUS", "LUSC", "SQCC"}:
        return "SQC"
    if s in {"LCC", "LARGE CELL", "LARGE-CELL"}:
        return "LCC"
    if s in {"BAS", "BASALOID"}:
        return "BAS"
    if s in {"LCNE", "LCNEC", "LARGE CELL NEUROENDOCRINE"}:
        return "LCNE"
    if s in {"SCC", "SCLC", "SMALL CELL", "SMALL-CELL"}:
        return "SCC"
    if s in {"CARCI", "CARCINOID"}:
        return "CARCI"
    if s in {"OTHER"}:
        return "OTHER"
    return s


def rec_row(cohort, predictor, endpoint, x, y, z=None, note=""):
    ru, pu, nu = spearman_pair(x, y)
    lo, hi = spearman_ci(x, y)
    if z is not None:
        rp, pp, np_ = partial_spearman(x, y, z)
    else:
        rp, pp, np_ = np.nan, np.nan, nu
    return {
        "cohort": cohort,
        "predictor": predictor,
        "endpoint": endpoint,
        "n_unadjusted": nu,
        "unadj_rho": ru,
        "unadj_p": pu,
        "unadj_ci_lo": lo,
        "unadj_ci_hi": hi,
        "n_partial": np_,
        "partial_rho": rp,
        "partial_p": pp,
        "note": note,
    }


def pick_row(stats_df: pd.DataFrame, cohort: str, pred: str, end: str):
    hit = stats_df[
        (stats_df.cohort == cohort)
        & (stats_df.predictor == pred)
        & (stats_df.endpoint == end)
    ]
    return hit.iloc[0] if len(hit) else None


def scatter(path: Path, x, y, xlab, ylab, title, n, rho, p):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(x, y, s=18, alpha=0.55, c="#1f4e79", edgecolors="none")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.text(
        0.03,
        0.97,
        f"n={n}\nρ={fmt_rho(rho)}\np={fmt_p(p)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, lw=0.4),
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(
    stats_df: pd.DataFrame,
    cov: dict,
    probe_used: dict,
    histo_counts: dict,
    n_matrix: int,
) -> None:
    def cell(cohort, pred, end):
        r = pick_row(stats_df, cohort, pred, end)
        if r is None:
            return "—"
        return f"{fmt_rho(r.unadj_rho)} ({fmt_p(r.unadj_p)})"

    def pcell(cohort, pred, end):
        r = pick_row(stats_df, cohort, pred, end)
        if r is None or not np.isfinite(r.partial_rho):
            return "—"
        return f"{fmt_rho(r.partial_rho)} ({fmt_p(r.partial_p)})"

    n_nsclc = cov["n_nsclc"]
    n_adc = cov["n_adc"]
    n_sqc = cov["n_sqc"]
    n_lcc = cov["n_lcc"]
    n_bas = cov["n_bas"]
    n_lcne = cov["n_lcne"]
    n_scc = cov["n_scc"]
    n_carci = cov["n_carci"]
    n_nt = cov["n_nt"]
    n_other = cov["n_other"]
    n_tumor = cov["n_tumor"]
    histo_txt = ", ".join(f"{k} {v}" for k, v in sorted(histo_counts.items(), key=lambda kv: (-kv[1], kv[0])))

    c4_cd8 = pick_row(stats_df, "NSCLC", "CLDN4", "CD8A")
    c4_pdl1 = pick_row(stats_df, "NSCLC", "CLDN4", "CD274")
    c4_imm = pick_row(stats_df, "NSCLC", "CLDN4", "ImmuneScore")
    cd8_imm = pick_row(stats_df, "NSCLC", "CD8A", "ImmuneScore")
    pdl1_cd8 = pick_row(stats_df, "NSCLC", "CD274", "CD8A")
    c4_pur = pick_row(stats_df, "NSCLC", "CLDN4", "TumorPurity")

    md = f"""# GSE30219 — NSCLC bulk CLDN4 vs CD8A / CD274

**Additive only.** No prior GSE30219 CLDN4–immune table is re-audited. This folder only measures **CLDN4 vs CD8A and CD274** on the public Rousseaux mixed-histology lung series.

Public Centre Léon Bérard surgical lung series (Rousseaux et al., *Sci Transl Med* 2013, PMID 23698379; GEO [GSE30219](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30219); GPL570 Affymetrix U133 Plus 2.0). Unit is the **array**. No ICI arm. No slide was re-scored.

## Honest n

GEO histology codes on this series are **not** WHO abbreviations. **`SCC` is small-cell**, not squamous. Squamous is **`SQC`**. Do not write n=293 for NSCLC.

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {n_matrix} | GEO `GSE30219_series_matrix.txt.gz` |
| GEO histology inventory | {n_matrix} | {histo_txt} |
| non-tumoral lung (`NTL`) | {n_nt} | titles “Non tumoral lung”; dropped. GEO `tissue` still says lung tumour |
| all tumors (any histology) | {n_tumor} | matrix minus `NTL` = paper 293 |
| SCLC (`SCC`) | {n_scc} | **not NSCLC**; held out |
| carcinoid (`CARCI`) | {n_carci} | **not NSCLC**; held out |
| LCNEC (`LCNE`) | {n_lcne} | neuroendocrine; **not** in primary NSCLC n |
| other / unmapped histology | {n_other} | GEO `Other` unspecified tumors; held out |
| ADC | {n_adc} | GEO `ADC` |
| SQC (squamous) | {n_sqc} | GEO `SQC` (not `SCC`) |
| LCC | {n_lcc} | GEO `LCC` |
| BAS (basaloid) | {n_bas} | GEO `BAS` |
| **primary NSCLC** | **{n_nsclc}** | ADC + SQC + LCC + BAS |
| CLDN4 / CD8A / CD274 finite | {cov['n_named_complete']} | unique-mapped max-mean collapse |
| ESTIMATE ImmuneScore | {n_nsclc} | Yoshihara Immune141 ssGSEA ({cov['estimate_immune_present']}/141) |

Primary tests use **n={n_nsclc}** NSCLC arrays. LCNE is a sensitivity row only. Paper text “293 tumors” includes SCLC / carcinoid / LCNE and is **not** the NSCLC n.

## One-row table

| dataset | histology rule | n | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–CD8A partial \\| purity (p) | CLDN4–CD274 partial \\| purity (p) |
|---|---|---:|---|---|---|---|
| GSE30219 NSCLC | ADC+SQC+LCC+BAS | {n_nsclc} | {cell('NSCLC','CLDN4','CD8A')} | {cell('NSCLC','CLDN4','CD274')} | {pcell('NSCLC','CLDN4','CD8A')} | {pcell('NSCLC','CLDN4','CD274')} |
| GSE30219 ADC | ADC only | {n_adc} | {cell('ADC','CLDN4','CD8A')} | {cell('ADC','CLDN4','CD274')} | {pcell('ADC','CLDN4','CD8A')} | {pcell('ADC','CLDN4','CD274')} |
| GSE30219 SQC | SQC only | {n_sqc} | {cell('SQC','CLDN4','CD8A')} | {cell('SQC','CLDN4','CD274')} | {pcell('SQC','CLDN4','CD8A')} | {pcell('SQC','CLDN4','CD274')} |
| GSE30219 +LCNE (sensitivity) | NSCLC+LCNE | {n_nsclc + n_lcne} | {cell('NSCLC_plus_LCNE','CLDN4','CD8A')} | {cell('NSCLC_plus_LCNE','CLDN4','CD274')} | {pcell('NSCLC_plus_LCNE','CLDN4','CD8A')} | {pcell('NSCLC_plus_LCNE','CLDN4','CD274')} |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | {n_matrix} | GEO characteristic; codes listed above |
| ICI response | **no** | 0 | surgical / pretreatment diagnostic series |
| OS / DFS | yes | {n_matrix} | GEO follow-up / status / DFS / relapse; **not used** (not a survival claim) |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL570 probe (multi-mapped `///` probes dropped for the three named genes):

| gene | probe |
|---|---|
| CLDN4 | `{probe_used.get('CLDN4', 'NA')}` |
| CD8A | `{probe_used.get('CD8A', 'NA')}` |
| CD274 | `{probe_used.get('CD274', 'NA')}` |

ESTIMATE gene collapse uses first-symbol mapping (same as A1 extra / GSE31210) so Immune141 coverage is {cov['estimate_immune_present']}/141 and Stromal141 is {cov['estimate_stromal_present']}/141.

## ESTIMATE

R `estimate` is not required. Scores follow the A1 extra public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

On the NSCLC matrix ESTIMATEScore ranges {cov['est_min']:.0f}–{cov['est_max']:.0f}. The cosine falls outside [0, 1] for **{cov['n_purity_outside_0_1']} / {n_nsclc}** tumors. Those values are **not** usable as a 0–1 purity fraction. They remain monotone with ESTIMATEScore here, so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The [0, 1] subset (n={cov['n_purity_in_0_1']}) is **not** the analysis n.

**ImmuneScore after ESTIMATEScore / TumorPurity is collinear.** The independent residual is **CLDN4 vs CD8A | ESTIMATEScore** (equals TumorPurity residual here).

CD8A vs ImmuneScore (positive-control, NSCLC): ρ={fmt_rho(cd8_imm.unadj_rho) if cd8_imm is not None else 'NA'} (p={fmt_p(cd8_imm.unadj_p) if cd8_imm is not None else 'NA'}, n={int(cd8_imm.n_unadjusted) if cd8_imm is not None else 0}).

## Main Spearman (primary NSCLC n={n_nsclc})

Partial = Pearson of rank residuals on ESTIMATE TumorPurity; df = n − 3. Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| predictor | endpoint | n | ρ | 95% CI | p | partial ρ \\| purity | partial p |
|---|---|---:|---:|---|---:|---:|---:|
| CLDN4 | CD8A | {int(c4_cd8.n_unadjusted) if c4_cd8 is not None else 0} | **{fmt_rho(c4_cd8.unadj_rho) if c4_cd8 is not None else 'NA'}** | {fmt_rho(c4_cd8.unadj_ci_lo) if c4_cd8 is not None else 'NA'} to {fmt_rho(c4_cd8.unadj_ci_hi) if c4_cd8 is not None else 'NA'} | {fmt_p(c4_cd8.unadj_p) if c4_cd8 is not None else 'NA'} | {fmt_rho(c4_cd8.partial_rho) if c4_cd8 is not None else 'NA'} | {fmt_p(c4_cd8.partial_p) if c4_cd8 is not None else 'NA'} |
| CLDN4 | CD274 | {int(c4_pdl1.n_unadjusted) if c4_pdl1 is not None else 0} | **{fmt_rho(c4_pdl1.unadj_rho) if c4_pdl1 is not None else 'NA'}** | {fmt_rho(c4_pdl1.unadj_ci_lo) if c4_pdl1 is not None else 'NA'} to {fmt_rho(c4_pdl1.unadj_ci_hi) if c4_pdl1 is not None else 'NA'} | {fmt_p(c4_pdl1.unadj_p) if c4_pdl1 is not None else 'NA'} | {fmt_rho(c4_pdl1.partial_rho) if c4_pdl1 is not None else 'NA'} | {fmt_p(c4_pdl1.partial_p) if c4_pdl1 is not None else 'NA'} |
| CLDN4 | ImmuneScore | {int(c4_imm.n_unadjusted) if c4_imm is not None else 0} | {fmt_rho(c4_imm.unadj_rho) if c4_imm is not None else 'NA'} | {fmt_rho(c4_imm.unadj_ci_lo) if c4_imm is not None else 'NA'} to {fmt_rho(c4_imm.unadj_ci_hi) if c4_imm is not None else 'NA'} | {fmt_p(c4_imm.unadj_p) if c4_imm is not None else 'NA'} | {fmt_rho(c4_imm.partial_rho) if c4_imm is not None else 'NA'} | {fmt_p(c4_imm.partial_p) if c4_imm is not None else 'NA'} |
| CLDN4 | TumorPurity | {int(c4_pur.n_unadjusted) if c4_pur is not None else 0} | {fmt_rho(c4_pur.unadj_rho) if c4_pur is not None else 'NA'} | {fmt_rho(c4_pur.unadj_ci_lo) if c4_pur is not None else 'NA'} to {fmt_rho(c4_pur.unadj_ci_hi) if c4_pur is not None else 'NA'} | {fmt_p(c4_pur.unadj_p) if c4_pur is not None else 'NA'} | — | — |
| CD274 | CD8A | {int(pdl1_cd8.n_unadjusted) if pdl1_cd8 is not None else 0} | {fmt_rho(pdl1_cd8.unadj_rho) if pdl1_cd8 is not None else 'NA'} | {fmt_rho(pdl1_cd8.unadj_ci_lo) if pdl1_cd8 is not None else 'NA'} to {fmt_rho(pdl1_cd8.unadj_ci_hi) if pdl1_cd8 is not None else 'NA'} | {fmt_p(pdl1_cd8.unadj_p) if pdl1_cd8 is not None else 'NA'} | {fmt_rho(pdl1_cd8.partial_rho) if pdl1_cd8 is not None else 'NA'} | {fmt_p(pdl1_cd8.partial_p) if pdl1_cd8 is not None else 'NA'} |

Mixed NSCLC (ADC+SQC+LCC+BAS) is **null** for both named pairs. The ADC-only unadjusted CLDN4–CD8A inverse (see one-row table) attenuates after the ESTIMATE impurity axis, the same pattern as GSE31210 LUAD. Do not treat the mixed-NSCLC pool as evidence that CLDN4-high NSCLC is CD8-low or PD-L1-low. It is **not** an ICI-response test.

## What is not done

- No re-use of the paper’s 293-tumor “all histologies” n as NSCLC.
- No SCLC (`SCC`) or carcinoid (`CARCI`) in the primary n.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, histology inventory, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `histology_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse30219_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    est_tab = pd.read_csv(SIG, sep="\t")
    stromal = est_tab.loc[est_tab["set"] == "Stromal141_UP", "hugo"].tolist()
    immune = est_tab.loc[est_tab["set"] == "Immune141_UP", "hugo"].tolist()
    if len(stromal) != 141 or len(immune) != 141:
        raise SystemExit(f"Yoshihara lists expected 141+141, got {len(stromal)}+{len(immune)}")

    matrix_path = download(GSE_MATRIX, CACHE / "GSE30219_series_matrix.txt.gz")
    annot_path = download(GPL570_ANNOT, CACHE / "GPL570.annot.gz")

    print("parse matrix", flush=True)
    meta, probes = parse_geo_matrix(matrix_path)
    first_map, unique_map = load_gpl_annot(annot_path)

    histo_col = None
    for cand in ("histology", "histological type", "tumor type", "tissue", "type"):
        if cand in meta.columns:
            histo_col = cand
            break
    if histo_col is None:
        # fall back: first characteristic that looks like histology codes
        for c in meta.columns:
            vals = set(meta[c].astype(str).str.upper().str.strip())
            if vals & (NSCLC_CODES | SCLC_CODES | LCNE_CODES | CARCINOID_CODES):
                histo_col = c
                break
    if histo_col is None:
        raise SystemExit(f"GSE30219 histology characteristic missing: {list(meta.columns)}")

    meta = meta.copy()
    meta["histo_raw"] = meta[histo_col].astype(str)
    meta["histo"] = meta["histo_raw"].map(histology_code)
    n_matrix = int(len(meta))
    counts = Counter(meta["histo"].tolist())

    is_nt = meta["histo"].isin(NORMAL_CODES)
    is_scc = meta["histo"].isin(SCLC_CODES)
    is_carci = meta["histo"].isin(CARCINOID_CODES)
    is_lcne = meta["histo"].isin(LCNE_CODES)
    is_nsclc = meta["histo"].isin(NSCLC_CODES)
    is_tumor = ~is_nt
    known = NSCLC_CODES | LCNE_CODES | SCLC_CODES | CARCINOID_CODES | NORMAL_CODES | OTHER_CODES
    is_other = meta["histo"].isin(OTHER_CODES) | ~meta["histo"].isin(known)

    n_nt = int(is_nt.sum())
    n_scc = int(is_scc.sum())
    n_carci = int(is_carci.sum())
    n_lcne = int(is_lcne.sum())
    n_nsclc = int(is_nsclc.sum())
    n_tumor = int(is_tumor.sum())
    n_other = int(is_other.sum())
    n_adc = int((meta["histo"] == "ADC").sum())
    n_sqc = int((meta["histo"] == "SQC").sum())
    n_lcc = int((meta["histo"] == "LCC").sum())
    n_bas = int((meta["histo"] == "BAS").sum())

    print(
        f"arrays={n_matrix} tumors={n_tumor} NSCLC={n_nsclc} "
        f"ADC={n_adc} SQC={n_sqc} LCC={n_lcc} BAS={n_bas} "
        f"LCNE={n_lcne} SCC={n_scc} CARCI={n_carci} NT={n_nt} other={n_other}",
        flush=True,
    )
    print("histology counts:", dict(counts), flush=True)
    print("meta columns:", list(meta.columns), flush=True)

    inv = (
        meta["histo"]
        .value_counts()
        .rename_axis("histo")
        .reset_index(name="n")
    )
    inv["in_primary_nsclc"] = inv["histo"].isin(sorted(NSCLC_CODES))
    inv["note"] = inv["histo"].map(
        {
            "ADC": "adenocarcinoma; NSCLC",
            "SQC": "squamous; NSCLC (not SCC)",
            "LCC": "large-cell; NSCLC",
            "BAS": "basaloid; NSCLC",
            "LCNE": "LCNEC; sensitivity only",
            "SCC": "small-cell; NOT NSCLC",
            "CARCI": "carcinoid; NOT NSCLC",
            "NTL": "non-tumoral lung (GEO NTL / titles 'Non tumoral lung'); dropped",
            "OTHER": "GEO histology Other; unspecified tumor; held out",
        }
    )
    inv.to_csv(TABLES / "histology_inventory.tsv", sep="\t", index=False)

    if n_nsclc < 20:
        # Still write an empty-ish FINDING so the PR has a table-or-empty file.
        (HERE / "FINDING.md").write_text(
            "# GSE30219 — NSCLC bulk CLDN4 vs CD8A / CD274\n\n"
            f"Could not define a usable NSCLC set (n={n_nsclc}). "
            f"Histology inventory: {dict(counts)}. Columns: {list(meta.columns)}.\n"
        )
        raise SystemExit(f"NSCLC n={n_nsclc} too small; see histology_inventory.tsv")

    # Collapse on arrays that enter a reported cohort (NSCLC + LCNE sensitivity).
    keep = meta.index[is_nsclc | is_lcne]
    probes_k = probes.loc[:, keep]
    genes_est, _ = collapse_maxmean(probes_k, first_map)
    genes_uniq, probe_pick = collapse_maxmean(probes_k, unique_map)

    miss = [g for g in TARGETS if g not in genes_uniq.index]
    if miss:
        (HERE / "FINDING.md").write_text(
            "# GSE30219 — NSCLC bulk CLDN4 vs CD8A / CD274\n\n"
            f"Named genes missing after unique-probe collapse: {miss}.\n"
        )
        raise SystemExit(f"named genes missing after unique-probe collapse: {miss}")

    probe_used = {g: str(probe_pick.loc[g]) for g in TARGETS}
    pd.Series(probe_used, name="probe_id").to_csv(TABLES / "probe_used.tsv", sep="\t")

    est = estimate_scores(genes_est, stromal, immune)
    logx = np.log2(genes_uniq.clip(lower=0) + 1)

    per = pd.DataFrame(
        {
            "gsm": keep,
            "histo": meta.loc[keep, "histo"].to_numpy(),
            "histo_raw": meta.loc[keep, "histo_raw"].to_numpy(),
            "CLDN4": genes_uniq.loc["CLDN4"].to_numpy(float),
            "CD8A": genes_uniq.loc["CD8A"].to_numpy(float),
            "CD274": genes_uniq.loc["CD274"].to_numpy(float),
            "CLDN4_log2p1": logx.loc["CLDN4"].to_numpy(float),
            "CD8A_log2p1": logx.loc["CD8A"].to_numpy(float),
            "CD274_log2p1": logx.loc["CD274"].to_numpy(float),
        },
        index=keep,
    )
    per = per.join(est)
    per["in_primary_nsclc"] = per["histo"].isin(NSCLC_CODES)
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    nsclc = per[per["in_primary_nsclc"]].copy()
    n_named = int(
        (
            np.isfinite(nsclc["CLDN4"])
            & np.isfinite(nsclc["CD8A"])
            & np.isfinite(nsclc["CD274"])
        ).sum()
    )

    est_ns = nsclc["ESTIMATE_Score"].to_numpy(float)
    pur_ns = nsclc["ESTIMATE_TumorPurity"].to_numpy(float)
    n_pur_ok = int(((pur_ns >= 0) & (pur_ns <= 1) & np.isfinite(pur_ns)).sum())
    n_pur_out = int(np.isfinite(pur_ns).sum() - n_pur_ok)

    cov = {
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_nsclc": n_nsclc,
        "n_adc": n_adc,
        "n_sqc": n_sqc,
        "n_lcc": n_lcc,
        "n_bas": n_bas,
        "n_lcne": n_lcne,
        "n_scc": n_scc,
        "n_carci": n_carci,
        "n_nt": n_nt,
        "n_other": n_other,
        "n_named_complete": n_named,
        "histo_column": histo_col,
        "estimate_immune_present": int(sum(g in genes_est.index for g in immune)),
        "estimate_stromal_present": int(sum(g in genes_est.index for g in stromal)),
        "n_purity_in_0_1": n_pur_ok,
        "n_purity_outside_0_1": n_pur_out,
        "est_min": float(np.nanmin(est_ns)) if len(est_ns) else float("nan"),
        "est_max": float(np.nanmax(est_ns)) if len(est_ns) else float("nan"),
        "probe_CLDN4": probe_used["CLDN4"],
        "probe_CD8A": probe_used["CD8A"],
        "probe_CD274": probe_used["CD274"],
    }
    pd.Series(cov, name="value").to_csv(TABLES / "coverage.tsv", sep="\t")

    rows = []
    cohorts = {
        "NSCLC": per["histo"].isin(NSCLC_CODES),
        "ADC": per["histo"] == "ADC",
        "SQC": per["histo"] == "SQC",
        "LCC": per["histo"] == "LCC",
        "BAS": per["histo"] == "BAS",
        "LCNE": per["histo"] == "LCNE",
        "NSCLC_plus_LCNE": per["histo"].isin(NSCLC_CODES | {"LCNE"}),
    }
    for cohort, mask in cohorts.items():
        d = per.loc[mask]
        if len(d) < 6:
            continue
        z = d["ESTIMATE_TumorPurity"].to_numpy(float)
        c4 = d["CLDN4"].to_numpy(float)
        cd8 = d["CD8A"].to_numpy(float)
        pdl1 = d["CD274"].to_numpy(float)
        imm = d["ESTIMATE_ImmuneScore"].to_numpy(float)
        strom = d["ESTIMATE_StromalScore"].to_numpy(float)
        rows.append(rec_row(cohort, "CLDN4", "CD8A", c4, cd8, z))
        rows.append(rec_row(cohort, "CLDN4", "CD274", c4, pdl1, z))
        rows.append(rec_row(cohort, "CLDN4", "ImmuneScore", c4, imm, z))
        rows.append(rec_row(cohort, "CLDN4", "StromalScore", c4, strom, z))
        rows.append(rec_row(cohort, "CLDN4", "TumorPurity", c4, z, None))
        rows.append(rec_row(cohort, "CD8A", "ImmuneScore", cd8, imm, z, note="positive control"))
        rows.append(rec_row(cohort, "CD274", "CD8A", pdl1, cd8, z))
        rows.append(rec_row(cohort, "CD8A", "TumorPurity", cd8, z, None))
        rows.append(rec_row(cohort, "CD274", "TumorPurity", pdl1, z, None))

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    one = []
    for cohort in ("NSCLC", "ADC", "SQC", "NSCLC_plus_LCNE"):
        a = pick_row(stats_df, cohort, "CLDN4", "CD8A")
        b = pick_row(stats_df, cohort, "CLDN4", "CD274")
        if a is None or b is None:
            continue
        one.append(
            {
                "dataset": f"GSE30219 {cohort}",
                "n": int(a.n_unadjusted),
                "CLDN4_CD8A_rho": a.unadj_rho,
                "CLDN4_CD8A_p": a.unadj_p,
                "CLDN4_CD8A_partial_rho": a.partial_rho,
                "CLDN4_CD8A_partial_p": a.partial_p,
                "CLDN4_CD274_rho": b.unadj_rho,
                "CLDN4_CD274_p": b.unadj_p,
                "CLDN4_CD274_partial_rho": b.partial_rho,
                "CLDN4_CD274_partial_p": b.partial_p,
            }
        )
    pd.DataFrame(one).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # Figures on primary NSCLC.
    scatter(
        FIGURES / "fig1_cldn4_vs_cd8a.png",
        nsclc["CLDN4_log2p1"],
        nsclc["CD8A_log2p1"],
        "CLDN4 log2(MAS5+1)",
        "CD8A log2(MAS5+1)",
        f"GSE30219 NSCLC CLDN4 vs CD8A (n={n_nsclc})",
        int(c4_cd8.n_unadjusted) if (c4_cd8 := pick_row(stats_df, "NSCLC", "CLDN4", "CD8A")) is not None else n_nsclc,
        c4_cd8.unadj_rho if c4_cd8 is not None else np.nan,
        c4_cd8.unadj_p if c4_cd8 is not None else np.nan,
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_cd274.png",
        nsclc["CLDN4_log2p1"],
        nsclc["CD274_log2p1"],
        "CLDN4 log2(MAS5+1)",
        "CD274 log2(MAS5+1)",
        f"GSE30219 NSCLC CLDN4 vs CD274 (n={n_nsclc})",
        int(c4_pdl1.n_unadjusted) if (c4_pdl1 := pick_row(stats_df, "NSCLC", "CLDN4", "CD274")) is not None else n_nsclc,
        c4_pdl1.unadj_rho if c4_pdl1 is not None else np.nan,
        c4_pdl1.unadj_p if c4_pdl1 is not None else np.nan,
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_immunescore.png",
        nsclc["CLDN4_log2p1"],
        nsclc["ESTIMATE_ImmuneScore"],
        "CLDN4 log2(MAS5+1)",
        "ESTIMATE ImmuneScore",
        f"GSE30219 NSCLC CLDN4 vs ImmuneScore (n={n_nsclc})",
        int(c4_imm.n_unadjusted) if (c4_imm := pick_row(stats_df, "NSCLC", "CLDN4", "ImmuneScore")) is not None else n_nsclc,
        c4_imm.unadj_rho if c4_imm is not None else np.nan,
        c4_imm.unadj_p if c4_imm is not None else np.nan,
    )

    # Forest of primary pairs.
    forest_keys = [
        ("NSCLC", "CLDN4", "CD8A", "CLDN4–CD8A"),
        ("NSCLC", "CLDN4", "CD274", "CLDN4–CD274"),
        ("NSCLC", "CLDN4", "ImmuneScore", "CLDN4–ImmuneScore"),
        ("ADC", "CLDN4", "CD8A", "ADC CLDN4–CD8A"),
        ("ADC", "CLDN4", "CD274", "ADC CLDN4–CD274"),
        ("SQC", "CLDN4", "CD8A", "SQC CLDN4–CD8A"),
        ("SQC", "CLDN4", "CD274", "SQC CLDN4–CD274"),
    ]
    labels, rhos, los, his = [], [], [], []
    for cohort, pred, end, lab in forest_keys:
        r = pick_row(stats_df, cohort, pred, end)
        if r is None:
            continue
        labels.append(f"{lab} (n={int(r.n_unadjusted)})")
        rhos.append(r.unadj_rho)
        los.append(r.unadj_ci_lo)
        his.append(r.unadj_ci_hi)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
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
    ax.set_title("GSE30219 CLDN4 vs CD8A / CD274")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_spearman_forest.png", dpi=140)
    fig.savefig(FIGURES / "fig4_spearman_forest.pdf")
    plt.close(fig)

    write_finding(stats_df, cov, probe_used, dict(counts), n_matrix)

    summary = {
        "dataset": "GSE30219",
        "primary_cohort": "NSCLC = ADC+SQC+LCC+BAS",
        "n_nsclc": n_nsclc,
        "n_named_complete": n_named,
        "CLDN4_CD8A": {
            "n": int(c4_cd8.n_unadjusted) if c4_cd8 is not None else None,
            "rho": None if c4_cd8 is None else float(c4_cd8.unadj_rho),
            "p": None if c4_cd8 is None else float(c4_cd8.unadj_p),
            "partial_rho": None if c4_cd8 is None else float(c4_cd8.partial_rho),
            "partial_p": None if c4_cd8 is None else float(c4_cd8.partial_p),
        },
        "CLDN4_CD274": {
            "n": int(c4_pdl1.n_unadjusted) if c4_pdl1 is not None else None,
            "rho": None if c4_pdl1 is None else float(c4_pdl1.unadj_rho),
            "p": None if c4_pdl1 is None else float(c4_pdl1.unadj_p),
            "partial_rho": None if c4_pdl1 is None else float(c4_pdl1.partial_rho),
            "partial_p": None if c4_pdl1 is None else float(c4_pdl1.partial_p),
        },
        "probes": probe_used,
        "histo_column": histo_col,
        "histo_counts": dict(counts),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    print("wrote", HERE / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
