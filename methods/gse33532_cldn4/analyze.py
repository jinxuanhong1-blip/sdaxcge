#!/usr/bin/env python3
"""GSE33532 LUAD tumors: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only. No dual-high. Tumor only. Epithelial residual.

Public Meister / Thoraxklinik Heidelberg Affymetrix U133 Plus 2.0 series
(Meister et al., J Bioinformatics Research Studies 2014; GEO GSE33532).
Four tumor sub-sites (A–D) plus matched normal lung from 20 stage I–II
NSCLC patients (n=100 arrays). GEO histology is adeno / mixed / squamous.

Primary cohort is LUAD tumor arrays. Matched normal is dropped. Mixed and
squamous are inventoried and are not LUAD. Independent n is the patient
(10 LUAD), not the 40 arrays and not the 100-array series.

Downloads stay under $GSE33532_CLDN4_DATA (default /tmp/gse33532_cldn4)
and are not committed.
"""
from __future__ import annotations

import gzip
import json
import math
import os
import re
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
DATA = Path(os.environ.get("GSE33532_CLDN4_DATA", "/tmp/gse33532_cldn4"))
SIG = HERE / "estimate_yoshihara_2013.tsv"

SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33532/"
    "matrix/GSE33532_series_matrix.txt.gz"
)
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

UA = (
    "sdaxcge-gse33532-cldn4/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

# Named Plus2 probes (same as GSE19188 / GSE4573). 227458_at is not on
# the deposited GSE33532 matrix (filtered 25,906 / 54,675 probe sets).
NAMED_PROBES = {
    "CLDN4": "201428_at",
    "CD8A": "205758_at",
    "CD274": "223834_at",
}

EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
TARGETS = ["CLDN4", "CD8A", "CD274"]


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


def fmt_rho_signed(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[ok], y[ok]
    n = int(xx.size)
    if n < 6 or np.unique(xx).size < 3 or np.unique(yy).size < 3:
        return dict(n=n, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    rho, p = stats.spearmanr(xx, yy)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(xx[idx]).size < 3 or np.unique(yy[idx]).size < 3:
            continue
        r, _ = stats.spearmanr(xx[idx], yy[idx])
        if np.isfinite(r):
            boots.append(float(r))
    if len(boots) < 50:
        lo = hi = np.nan
    else:
        lo, hi = np.percentile(boots, [2.5, 97.5])
    return dict(n=n, rho=float(rho), p=float(p), ci_lo=float(lo), ci_hi=float(hi))


def partial_spearman(x, y, z):
    """Pearson of rank residuals; df = n − 3."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m]).reshape(-1, 1)
    rx = residualize(xr, zr)
    ry = residualize(yr, zr)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return dict(n=n, rho_adj=r, p_adj=p)


def verdict(rho_adj, p_adj, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


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


def parse_title(title: str) -> tuple[str, str]:
    m = re.search(
        r"patient\s+(\S+),\s+(tumor sub-sample ([ABCD])|matched normal lung)",
        str(title),
        flags=re.I,
    )
    if not m:
        return "NA", "NA"
    patient = m.group(1)
    site = "N" if "normal" in str(title).lower() else m.group(3)
    return patient, site


def histo_code(raw: str) -> str:
    s = str(raw).strip().lower()
    if s in {"adeno", "adenocarcinoma", "luad", "adc"}:
        return "LUAD"
    if s in {"squamous", "scc", "sqc", "lusc"}:
        return "LUSC"
    if s in {"mixed", "adenosquamous"}:
        return "MIXED"
    if "normal" in s:
        return "NORMAL"
    return s.upper() or "NA"


def pair_row(cohort, unit, n_patients, predictor, endpoint, x, y, z) -> dict:
    crude = spearman_ci(x, y)
    adj = partial_spearman(x, y, z)
    n_for_verdict = int(crude["n"])
    return {
        "cohort": cohort,
        "unit": unit,
        "n_patients": int(n_patients),
        "predictor": predictor,
        "endpoint": endpoint,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_lo": crude["ci_lo"],
        "ci_hi": crude["ci_hi"],
        "rho_adj_epithelial": adj["rho_adj"],
        "p_adj_epithelial": adj["p_adj"],
        "n_adj": adj["n"],
        "verdict": verdict(adj["rho_adj"], adj["p_adj"], n_for_verdict),
    }


def pick_row(df: pd.DataFrame, cohort: str, endpoint: str, unit: str | None = None):
    m = (df["cohort"] == cohort) & (df["predictor"] == "CLDN4") & (df["endpoint"] == endpoint)
    if unit is not None:
        m = m & (df["unit"] == unit)
    hit = df.loc[m]
    if hit.empty:
        return None
    return hit.iloc[0]


def scatter(path: Path, x, y, xlab, ylab, title, n, rho, p, patients=None):
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(x, y, s=22, c="#1f4e79", alpha=0.75, edgecolors="none")
    if patients is not None:
        tmp = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float), "p": patients})
        means = tmp.groupby("p", sort=False).mean(numeric_only=True)
        ax.scatter(
            means["x"],
            means["y"],
            s=48,
            facecolors="none",
            edgecolors="#c0392b",
            linewidths=1.2,
            label="patient mean",
            zorder=3,
        )
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nn={n}  ρ={fmt_rho(rho)}  p={fmt_p(p)}")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(stats_df, cov, probe_used, epi_used, counts):
    def row(cohort, endpoint, unit):
        return pick_row(stats_df, cohort, endpoint, unit)

    luad_a_cd8 = row("LUAD", "CD8A", "array")
    luad_a_pdl1 = row("LUAD", "CD274", "array")
    luad_a_imm = row("LUAD", "ImmuneScore", "array")
    luad_p_cd8 = row("LUAD", "CD8A", "patient")
    luad_p_pdl1 = row("LUAD", "CD274", "patient")
    luad_p_imm = row("LUAD", "ImmuneScore", "patient")
    luad_a_epi = row("LUAD", "epithelial_mean_z", "array")
    nsclc_a_cd8 = row("all_tumor", "CD8A", "array")
    nsclc_a_pdl1 = row("all_tumor", "CD274", "array")
    nsclc_p_cd8 = row("all_tumor", "CD8A", "patient")
    nsclc_p_pdl1 = row("all_tumor", "CD274", "patient")
    ctrl = pick_row(stats_df, "LUAD", "ImmuneScore", "array")
    # CD8A vs ImmuneScore is stored with predictor CD8A
    cd8_imm = stats_df[
        (stats_df["cohort"] == "LUAD")
        & (stats_df["predictor"] == "CD8A")
        & (stats_df["endpoint"] == "ImmuneScore")
        & (stats_df["unit"] == "array")
    ]
    cd8_imm = None if cd8_imm.empty else cd8_imm.iloc[0]

    def pair_md(r):
        if r is None:
            return "NA", "NA", "NA"
        return (
            f"{fmt_rho_signed(r.rho)} ({fmt_p(r.p)})",
            f"{fmt_rho_signed(r.rho_adj_epithelial)} ({fmt_p(r.p_adj_epithelial)})",
            r.verdict,
        )

    a8, a8a, a8v = pair_md(luad_a_cd8)
    a27, a27a, a27v = pair_md(luad_a_pdl1)
    ai, aia, aiv = pair_md(luad_a_imm)
    p8, p8a, p8v = pair_md(luad_p_cd8)
    p27, p27a, p27v = pair_md(luad_p_pdl1)
    pi, pia, piv = pair_md(luad_p_imm)

    n_luad_arr = int(cov["n_luad_tumor_arrays"])
    n_luad_pt = int(cov["n_luad_patients"])
    n_probe = int(cov["n_probes_deposited"])

    text = f"""# GSE33532 — LUAD array CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** No dual-high. No prior GSE33532 CLDN4–immune table is re-audited. This folder only measures **CLDN4 vs CD8A, CD274, and ESTIMATE ImmuneScore** on public LUAD **tumor** arrays, after an epithelial residual.

Public Meister / Thoraxklinik Heidelberg surgical series (Meister et al., *J Bioinformatics Research Studies* 2014, 1(1):1; GEO [GSE33532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33532); GPL570 Affymetrix U133 Plus 2.0, RMA-PLM). Unit for independence is the **patient**. No ICI arm. No slide was re-scored.

## Honest n

Do not write n=100. The series is 20 patients × (4 tumor sub-sites + 1 matched normal). GEO `histology` is **adeno / mixed / squamous / normal lung**. Mixed is not LUAD.

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 100 | GEO `GSE33532_series_matrix.txt.gz` |
| deposited probe sets | {n_probe} | filtered GPL570 (full Plus2 is 54,675); not a re-annotation |
| unique GSM | 100 | all unique |
| patients (GEO design) | 20 | `overall_design`: 4 sites A–D + matched normal from 20 patients |
| matched normal lung | 20 | titles “matched normal lung”; **dropped** |
| all tumor arrays | 80 | 20 patients × 4 sites; **not** the LUAD n |
| LUAD tumor arrays (`histology: adeno`) | **{n_luad_arr}** | 10 patients × 4 sites |
| LUAD patients | **{n_luad_pt}** | **independent n** |
| mixed tumor arrays | {int(cov["n_mixed_tumor_arrays"])} | 6 patients; not LUAD |
| squamous tumor arrays | {int(cov["n_lusc_tumor_arrays"])} | 4 patients; not LUAD |
| ICI / response | 0 | resected early-stage atlas; not deposited |
| OS / DFS | 0 | not a GEO characteristic |
| tumor % / ABSOLUTE purity | 0 | only public proxy used here is an RNA epithelial mean-z |
| CLDN4 finite (`{probe_used["CLDN4"]}`) | {n_luad_arr} | named Plus2 probe |
| CD8A finite (`{probe_used["CD8A"]}`) | {n_luad_arr} | named Plus2 probe |
| CD274 finite (`{probe_used["CD274"]}`) | {n_luad_arr} | named Plus2 probe; `227458_at` **absent** from the deposited matrix |
| ESTIMATE ImmuneScore | {n_luad_arr} | Yoshihara Immune141 ssGSEA ({int(cov["estimate_immune_present"])}/141 genes on this matrix) |
| epithelial mean-z | {n_luad_arr} | {len(epi_used)}/{len(EPI_GENES)} of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |

Primary tests use **{n_luad_pt} LUAD patients** (mean of the four tumor sites). Array-level n={n_luad_arr} is the same 10 tumors counted four times and is **not** an independent n. HOLDS (PR 229 / GSE4573) requires n≥40, ρ_adj<0, p_adj<0.05 — the independent LUAD n is **UNDERPOWERED**.

## One-row table

| dataset | histology rule | n arrays | n patients | CLDN4–CD8A ρ (p) | adj ρ \\| epi (p) | CLDN4–CD274 ρ (p) | adj ρ \\| epi (p) | CLDN4–ImmuneScore ρ (p) | adj ρ \\| epi (p) | verdict |
|---|---|---:|---:|---|---|---|---|---|---|---|
| GSE33532 LUAD (patient mean) | GEO `adeno` tumors only | {n_luad_arr} | **{n_luad_pt}** | {p8} | {p8a} | {p27} | {p27a} | {pi} | {pia} | **{p8v}** |
| GSE33532 LUAD (array; not independent) | same 10 tumors × 4 sites | {n_luad_arr} | {n_luad_pt} | {a8} | {a8a} | {a27} | {a27a} | {ai} | {aia} | {a8v} (pseudo-replicated) |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

Deposited series-matrix is RMA-PLM log2 intensity on a **filtered** 25,906-probe GPL570 subset. Spearman is rank-based.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | 100 | GEO `adeno` 40 / `mixed` 24 / `squamous` 16 / `normal lung` 20 |
| tumor stage | yes | 80 tumors | 1A / 1B / 2A / 2B; unused here |
| gender / age | yes | 100 | unused here |
| ICI response | **no** | 0 | surgical / pretreatment diagnostic series |
| ESTIMATE published scores | **no** | 0 | ImmuneScore computed here from Yoshihara lists |

Named Plus2 probes (forced; max-mean collapse is not used for the three named genes):

| gene | probe | on deposited matrix |
|---|---|---|
| CLDN4 | `{probe_used["CLDN4"]}` | yes |
| CD8A | `{probe_used["CD8A"]}` | yes |
| CD274 | `{probe_used["CD274"]}` | yes (`227458_at` missing) |

ESTIMATE gene collapse uses first-symbol mapping so Immune141 coverage is {int(cov["estimate_immune_present"])}/141 on this filtered matrix.

## Epithelial residual

Partial Spearman residualises ranks on the pan-epithelial mean-z ({", ".join(EPI_GENES)}; **{len(epi_used)}/{len(EPI_GENES)}** present: {", ".join(epi_used)}). No dual-high TACSTD2×CLDN4 split.

CLDN4 vs epithelial mean-z (LUAD arrays): ρ={fmt_rho_signed(luad_a_epi.rho) if luad_a_epi is not None else "NA"} (p={fmt_p(luad_a_epi.p) if luad_a_epi is not None else "NA"}, n={n_luad_arr}).

CD8A vs ImmuneScore (positive-control, LUAD arrays): ρ={fmt_rho_signed(cd8_imm.rho) if cd8_imm is not None else "NA"} (p={fmt_p(cd8_imm.p) if cd8_imm is not None else "NA"}, n={n_luad_arr}).

## Main Spearman (primary = LUAD patient mean, n={n_luad_pt})

Partial = Pearson of rank residuals on epithelial mean-z; df = n − 3. Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| predictor | endpoint | unit | n | ρ | 95% CI | p | ρ_adj \\| epi | p_adj | verdict |
|---|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 | CD8A | patient | {int(luad_p_cd8.n) if luad_p_cd8 is not None else n_luad_pt} | **{fmt_rho(luad_p_cd8.rho) if luad_p_cd8 is not None else "NA"}** | {fmt_rho(luad_p_cd8.ci_lo) if luad_p_cd8 is not None else "NA"} to {fmt_rho(luad_p_cd8.ci_hi) if luad_p_cd8 is not None else "NA"} | {fmt_p(luad_p_cd8.p) if luad_p_cd8 is not None else "NA"} | {fmt_rho(luad_p_cd8.rho_adj_epithelial) if luad_p_cd8 is not None else "NA"} | {fmt_p(luad_p_cd8.p_adj_epithelial) if luad_p_cd8 is not None else "NA"} | **{p8v}** |
| CLDN4 | CD274 | patient | {int(luad_p_pdl1.n) if luad_p_pdl1 is not None else n_luad_pt} | **{fmt_rho(luad_p_pdl1.rho) if luad_p_pdl1 is not None else "NA"}** | {fmt_rho(luad_p_pdl1.ci_lo) if luad_p_pdl1 is not None else "NA"} to {fmt_rho(luad_p_pdl1.ci_hi) if luad_p_pdl1 is not None else "NA"} | {fmt_p(luad_p_pdl1.p) if luad_p_pdl1 is not None else "NA"} | {fmt_rho(luad_p_pdl1.rho_adj_epithelial) if luad_p_pdl1 is not None else "NA"} | {fmt_p(luad_p_pdl1.p_adj_epithelial) if luad_p_pdl1 is not None else "NA"} | **{p27v}** |
| CLDN4 | ImmuneScore | patient | {int(luad_p_imm.n) if luad_p_imm is not None else n_luad_pt} | **{fmt_rho(luad_p_imm.rho) if luad_p_imm is not None else "NA"}** | {fmt_rho(luad_p_imm.ci_lo) if luad_p_imm is not None else "NA"} to {fmt_rho(luad_p_imm.ci_hi) if luad_p_imm is not None else "NA"} | {fmt_p(luad_p_imm.p) if luad_p_imm is not None else "NA"} | {fmt_rho(luad_p_imm.rho_adj_epithelial) if luad_p_imm is not None else "NA"} | {fmt_p(luad_p_imm.p_adj_epithelial) if luad_p_imm is not None else "NA"} | **{piv}** |
| CLDN4 | CD8A | array | {int(luad_a_cd8.n) if luad_a_cd8 is not None else n_luad_arr} | {fmt_rho(luad_a_cd8.rho) if luad_a_cd8 is not None else "NA"} | {fmt_rho(luad_a_cd8.ci_lo) if luad_a_cd8 is not None else "NA"} to {fmt_rho(luad_a_cd8.ci_hi) if luad_a_cd8 is not None else "NA"} | {fmt_p(luad_a_cd8.p) if luad_a_cd8 is not None else "NA"} | {fmt_rho(luad_a_cd8.rho_adj_epithelial) if luad_a_cd8 is not None else "NA"} | {fmt_p(luad_a_cd8.p_adj_epithelial) if luad_a_cd8 is not None else "NA"} | {a8v} |
| CLDN4 | CD274 | array | {int(luad_a_pdl1.n) if luad_a_pdl1 is not None else n_luad_arr} | {fmt_rho(luad_a_pdl1.rho) if luad_a_pdl1 is not None else "NA"} | {fmt_rho(luad_a_pdl1.ci_lo) if luad_a_pdl1 is not None else "NA"} to {fmt_rho(luad_a_pdl1.ci_hi) if luad_a_pdl1 is not None else "NA"} | {fmt_p(luad_a_pdl1.p) if luad_a_pdl1 is not None else "NA"} | {fmt_rho(luad_a_pdl1.rho_adj_epithelial) if luad_a_pdl1 is not None else "NA"} | {fmt_p(luad_a_pdl1.p_adj_epithelial) if luad_a_pdl1 is not None else "NA"} | {a27v} |
| CLDN4 | ImmuneScore | array | {int(luad_a_imm.n) if luad_a_imm is not None else n_luad_arr} | {fmt_rho(luad_a_imm.rho) if luad_a_imm is not None else "NA"} | {fmt_rho(luad_a_imm.ci_lo) if luad_a_imm is not None else "NA"} to {fmt_rho(luad_a_imm.ci_hi) if luad_a_imm is not None else "NA"} | {fmt_p(luad_a_imm.p) if luad_a_imm is not None else "NA"} | {fmt_rho(luad_a_imm.rho_adj_epithelial) if luad_a_imm is not None else "NA"} | {fmt_p(luad_a_imm.p_adj_epithelial) if luad_a_imm is not None else "NA"} | {aiv} |

Do not treat the 40 LUAD arrays as 40 patients. The independent test is n={n_luad_pt} and is **UNDERPOWERED** for HOLDS. This is **not** an ICI-response test.

## Sensitivity (not the claim)

All-tumor NSCLC (LUAD+MIXED+LUSC) is 80 arrays / 20 patients. Still not n=100.

| pair | unit | n | ρ (p) | adj ρ \\| epi (p) | verdict |
|---|---|---:|---|---|---|
| CLDN4 vs CD8A | patient | {int(nsclc_p_cd8.n) if nsclc_p_cd8 is not None else 20} | {fmt_rho_signed(nsclc_p_cd8.rho) if nsclc_p_cd8 is not None else "NA"} ({fmt_p(nsclc_p_cd8.p) if nsclc_p_cd8 is not None else "NA"}) | {fmt_rho_signed(nsclc_p_cd8.rho_adj_epithelial) if nsclc_p_cd8 is not None else "NA"} ({fmt_p(nsclc_p_cd8.p_adj_epithelial) if nsclc_p_cd8 is not None else "NA"}) | {nsclc_p_cd8.verdict if nsclc_p_cd8 is not None else "NA"} |
| CLDN4 vs CD274 | patient | {int(nsclc_p_pdl1.n) if nsclc_p_pdl1 is not None else 20} | {fmt_rho_signed(nsclc_p_pdl1.rho) if nsclc_p_pdl1 is not None else "NA"} ({fmt_p(nsclc_p_pdl1.p) if nsclc_p_pdl1 is not None else "NA"}) | {fmt_rho_signed(nsclc_p_pdl1.rho_adj_epithelial) if nsclc_p_pdl1 is not None else "NA"} ({fmt_p(nsclc_p_pdl1.p_adj_epithelial) if nsclc_p_pdl1 is not None else "NA"}) | {nsclc_p_pdl1.verdict if nsclc_p_pdl1 is not None else "NA"} |
| CLDN4 vs CD8A | array | {int(nsclc_a_cd8.n) if nsclc_a_cd8 is not None else 80} | {fmt_rho_signed(nsclc_a_cd8.rho) if nsclc_a_cd8 is not None else "NA"} ({fmt_p(nsclc_a_cd8.p) if nsclc_a_cd8 is not None else "NA"}) | {fmt_rho_signed(nsclc_a_cd8.rho_adj_epithelial) if nsclc_a_cd8 is not None else "NA"} ({fmt_p(nsclc_a_cd8.p_adj_epithelial) if nsclc_a_cd8 is not None else "NA"}) | {nsclc_a_cd8.verdict if nsclc_a_cd8 is not None else "NA"} |
| CLDN4 vs CD274 | array | {int(nsclc_a_pdl1.n) if nsclc_a_pdl1 is not None else 80} | {fmt_rho_signed(nsclc_a_pdl1.rho) if nsclc_a_pdl1 is not None else "NA"} ({fmt_p(nsclc_a_pdl1.p) if nsclc_a_pdl1 is not None else "NA"}) | {fmt_rho_signed(nsclc_a_pdl1.rho_adj_epithelial) if nsclc_a_pdl1 is not None else "NA"} ({fmt_p(nsclc_a_pdl1.p_adj_epithelial) if nsclc_a_pdl1 is not None else "NA"}) | {nsclc_a_pdl1.verdict if nsclc_a_pdl1 is not None else "NA"} |

## What is not done

- No n=100 (tumor+normal mix) and no n=80 as a LUAD n.
- No mixed or squamous arrays in the primary LUAD row.
- No TACSTD2 / CLDN4 dual-high split.
- No ICI ORR / PFS model (labels are not deposited).
- No survival model (labels are not deposited).

## Files

- `analyze.py` — GEO download, honest-n inventory, ESTIMATE ImmuneScore, epithelial residual, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `patient_means.tsv`, `probe_used.tsv`, `coverage.tsv`, `histology_inventory.tsv`, `label_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse33532_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = download(GSE_MATRIX, DATA / "GSE33532_series_matrix.txt.gz")
    annot_path = download(GPL570_ANNOT, DATA / "GPL570.annot.gz")
    # reuse already-downloaded cache names if present
    if (DATA / "matrix.gz").exists() and not matrix_path.exists():
        matrix_path = DATA / "matrix.gz"
    if (DATA / "gpl.gz").exists() and annot_path.stat().st_size < 10_000:
        annot_path = DATA / "gpl.gz"

    # Prefer the files we already fetched under short names.
    if (DATA / "matrix.gz").exists():
        matrix_path = DATA / "matrix.gz"
    if (DATA / "gpl.gz").exists():
        annot_path = DATA / "gpl.gz"

    sig = pd.read_csv(SIG, sep="\t")
    stromal = sig.loc[sig["set"].str.contains("Stromal", case=False), "hugo"].astype(str).tolist()
    immune = sig.loc[sig["set"].str.contains("Immune", case=False), "hugo"].astype(str).tolist()

    meta, probes = parse_geo_matrix(matrix_path)
    first_map, unique_map = load_gpl_annot(annot_path)

    parsed = meta["title"].map(parse_title)
    meta["patient"] = [p for p, _ in parsed]
    meta["site"] = [s for _, s in parsed]
    meta["histo"] = meta["histology"].map(histo_code)
    meta["is_tumor"] = meta["site"].ne("N") & meta["histo"].ne("NORMAL")
    meta["is_normal"] = ~meta["is_tumor"]
    meta["is_luad"] = meta["is_tumor"] & meta["histo"].eq("LUAD")

    n_matrix = int(probes.shape[1])
    n_probes = int(probes.shape[0])
    n_tumor = int(meta["is_tumor"].sum())
    n_normal = int(meta["is_normal"].sum())
    n_luad_arr = int(meta["is_luad"].sum())
    n_mixed_arr = int((meta["is_tumor"] & meta["histo"].eq("MIXED")).sum())
    n_lusc_arr = int((meta["is_tumor"] & meta["histo"].eq("LUSC")).sum())
    n_patients = int(meta["patient"].nunique())
    n_luad_pt = int(meta.loc[meta["is_luad"], "patient"].nunique())
    n_mixed_pt = int(meta.loc[meta["is_tumor"] & meta["histo"].eq("MIXED"), "patient"].nunique())
    n_lusc_pt = int(meta.loc[meta["is_tumor"] & meta["histo"].eq("LUSC"), "patient"].nunique())

    counts = Counter(meta["histo"].tolist())
    inv = (
        meta.groupby(["histo", "is_tumor"], dropna=False)
        .agg(n_arrays=("title", "size"), n_patients=("patient", "nunique"))
        .reset_index()
    )
    inv.to_csv(TABLES / "histology_inventory.tsv", sep="\t", index=False)

    if n_luad_arr < 8 or n_luad_pt < 4:
        raise SystemExit(
            f"LUAD n too small: arrays={n_luad_arr} patients={n_luad_pt}; "
            f"histo={dict(counts)}"
        )

    # Collapse on all arrays (scores need a common gene space); tests subset later.
    genes_est, _ = collapse_maxmean(probes, first_map)
    genes_uniq, probe_pick = collapse_maxmean(probes, unique_map)

    miss = [g for g in TARGETS if g not in genes_uniq.index and NAMED_PROBES[g] not in probes.index]
    if miss:
        raise SystemExit(f"named genes missing: {miss}")

    # Force named probes for the three headline genes.
    probe_used = {}
    for g, pid in NAMED_PROBES.items():
        if pid not in probes.index:
            raise SystemExit(f"named probe {pid} for {g} missing from deposited matrix")
        genes_uniq.loc[g] = probes.loc[pid].astype(float)
        probe_used[g] = pid
    pd.Series(probe_used, name="probe_id").to_csv(TABLES / "probe_used.tsv", sep="\t")

    imm_score = ssgsea(genes_est, immune)
    epi_present = [g for g in EPI_GENES if g in genes_est.index]
    if len(epi_present) < 3:
        raise SystemExit(f"epithelial genes missing: present={epi_present}")
    z = genes_est.loc[epi_present]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    epi_score = z.mean(axis=0)

    per = pd.DataFrame(
        {
            "gsm": meta.index,
            "patient": meta["patient"].to_numpy(),
            "site": meta["site"].to_numpy(),
            "histo": meta["histo"].to_numpy(),
            "histo_raw": meta["histology"].to_numpy(),
            "source": meta["source"].to_numpy(),
            "stage": meta["tumor_stage"].to_numpy() if "tumor_stage" in meta.columns else "",
            "is_tumor": meta["is_tumor"].to_numpy(),
            "is_luad": meta["is_luad"].to_numpy(),
            "CLDN4": genes_uniq.loc["CLDN4"].to_numpy(float),
            "CD8A": genes_uniq.loc["CD8A"].to_numpy(float),
            "CD274": genes_uniq.loc["CD274"].to_numpy(float),
            "ImmuneScore": imm_score.reindex(meta.index).to_numpy(float),
            "epithelial_mean_z": epi_score.reindex(meta.index).to_numpy(float),
        },
        index=meta.index,
    )
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    tumor = per[per["is_tumor"]].copy()
    luad = per[per["is_luad"]].copy()

    def patient_means(df: pd.DataFrame) -> pd.DataFrame:
        cols = ["CLDN4", "CD8A", "CD274", "ImmuneScore", "epithelial_mean_z"]
        out = df.groupby("patient", sort=True)[cols].mean()
        out["histo"] = df.groupby("patient")["histo"].first()
        out["n_sites"] = df.groupby("patient").size()
        return out

    luad_pt = patient_means(luad)
    tumor_pt = patient_means(tumor)
    luad_pt.to_csv(TABLES / "patient_means.tsv", sep="\t")

    n_named = int(
        (
            np.isfinite(luad["CLDN4"])
            & np.isfinite(luad["CD8A"])
            & np.isfinite(luad["CD274"])
        ).sum()
    )

    cov = {
        "n_matrix": n_matrix,
        "n_probes_deposited": n_probes,
        "n_patients": n_patients,
        "n_tumor_arrays": n_tumor,
        "n_normal_arrays": n_normal,
        "n_luad_tumor_arrays": n_luad_arr,
        "n_luad_patients": n_luad_pt,
        "n_mixed_tumor_arrays": n_mixed_arr,
        "n_mixed_patients": n_mixed_pt,
        "n_lusc_tumor_arrays": n_lusc_arr,
        "n_lusc_patients": n_lusc_pt,
        "n_named_complete_luad_arrays": n_named,
        "estimate_immune_present": int(sum(g in genes_est.index for g in immune)),
        "estimate_stromal_present": int(sum(g in genes_est.index for g in stromal)),
        "epithelial_present": len(epi_present),
        "probe_CLDN4": probe_used["CLDN4"],
        "probe_CD8A": probe_used["CD8A"],
        "probe_CD274": probe_used["CD274"],
        "cd274_227458_at_present": bool("227458_at" in probes.index),
    }
    pd.Series(cov, name="value").to_csv(TABLES / "coverage.tsv", sep="\t")

    labels = pd.DataFrame(
        [
            {"field": "arrays in series matrix", "public": "yes", "n": n_matrix, "note": "20 patients × 5 arrays; do not use as the test n"},
            {"field": "deposited probe sets", "public": "yes", "n": n_probes, "note": "filtered GPL570 RMA-PLM; full Plus2 is 54,675"},
            {"field": "patients", "public": "yes", "n": n_patients, "note": "GEO overall_design"},
            {"field": "matched normal arrays", "public": "yes", "n": n_normal, "note": "dropped; not mixed into tumor tests"},
            {"field": "all tumor arrays", "public": "yes", "n": n_tumor, "note": "4 sites × 20 patients; not the LUAD n"},
            {"field": "LUAD tumor arrays", "public": "yes", "n": n_luad_arr, "note": "GEO histology: adeno; 4 sites × 10 patients"},
            {"field": "LUAD patients (independent n)", "public": "yes", "n": n_luad_pt, "note": "PRIMARY n; patient-mean of 4 tumor sites"},
            {"field": "mixed tumor arrays", "public": "yes", "n": n_mixed_arr, "note": "not LUAD"},
            {"field": "squamous tumor arrays", "public": "yes", "n": n_lusc_arr, "note": "not LUAD"},
            {"field": "ICI / treatment", "public": "no", "n": 0, "note": "resected atlas, not an ICI series"},
            {"field": "OS / DFS", "public": "no", "n": 0, "note": "not deposited"},
            {"field": "tumor % / ABSOLUTE purity", "public": "no", "n": 0, "note": "RNA epithelial mean-z is the public proxy"},
            {"field": "CLDN4 finite on LUAD tumors", "public": "yes", "n": int(luad["CLDN4"].notna().sum()), "note": probe_used["CLDN4"]},
            {"field": "CD8A finite on LUAD tumors", "public": "yes", "n": int(luad["CD8A"].notna().sum()), "note": probe_used["CD8A"]},
            {"field": "CD274 finite on LUAD tumors", "public": "yes", "n": int(luad["CD274"].notna().sum()), "note": f"{probe_used['CD274']}; 227458_at absent"},
        ]
    )
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    rows = []
    cohorts = {
        "LUAD": (luad, luad_pt, n_luad_pt),
        "all_tumor": (tumor, tumor_pt, int(tumor["patient"].nunique())),
        "MIXED": (
            tumor[tumor["histo"] == "MIXED"],
            patient_means(tumor[tumor["histo"] == "MIXED"]),
            n_mixed_pt,
        ),
        "LUSC": (
            tumor[tumor["histo"] == "LUSC"],
            patient_means(tumor[tumor["histo"] == "LUSC"]),
            n_lusc_pt,
        ),
    }
    endpoints = ["CD8A", "CD274", "ImmuneScore", "epithelial_mean_z"]
    for cohort, (arr, pt, n_pt) in cohorts.items():
        if len(arr) < 6:
            continue
        z_arr = arr["epithelial_mean_z"].to_numpy(float)
        z_pt = pt["epithelial_mean_z"].to_numpy(float)
        c4_arr = arr["CLDN4"].to_numpy(float)
        c4_pt = pt["CLDN4"].to_numpy(float)
        for end in endpoints:
            rows.append(
                pair_row(cohort, "array", n_pt, "CLDN4", end, c4_arr, arr[end].to_numpy(float), z_arr)
            )
            if len(pt) >= 6:
                rows.append(
                    pair_row(cohort, "patient", n_pt, "CLDN4", end, c4_pt, pt[end].to_numpy(float), z_pt)
                )
        rows.append(
            pair_row(
                cohort,
                "array",
                n_pt,
                "CD8A",
                "ImmuneScore",
                arr["CD8A"].to_numpy(float),
                arr["ImmuneScore"].to_numpy(float),
                z_arr,
            )
        )
        if len(pt) >= 6:
            rows.append(
                pair_row(
                    cohort,
                    "patient",
                    n_pt,
                    "CD8A",
                    "ImmuneScore",
                    pt["CD8A"].to_numpy(float),
                    pt["ImmuneScore"].to_numpy(float),
                    z_pt,
                )
            )

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    one = []
    for cohort, unit in (("LUAD", "patient"), ("LUAD", "array"), ("all_tumor", "patient"), ("all_tumor", "array")):
        a = pick_row(stats_df, cohort, "CD8A", unit)
        b = pick_row(stats_df, cohort, "CD274", unit)
        c = pick_row(stats_df, cohort, "ImmuneScore", unit)
        if a is None or b is None or c is None:
            continue
        one.append(
            {
                "dataset": f"GSE33532 {cohort} {unit}",
                "n_arrays": n_luad_arr if cohort == "LUAD" else n_tumor,
                "n_patients": int(a.n_patients),
                "unit": unit,
                "CLDN4_CD8A_rho": a.rho,
                "CLDN4_CD8A_p": a.p,
                "CLDN4_CD8A_partial_rho": a.rho_adj_epithelial,
                "CLDN4_CD8A_partial_p": a.p_adj_epithelial,
                "CLDN4_CD274_rho": b.rho,
                "CLDN4_CD274_p": b.p,
                "CLDN4_CD274_partial_rho": b.rho_adj_epithelial,
                "CLDN4_CD274_partial_p": b.p_adj_epithelial,
                "CLDN4_ImmuneScore_rho": c.rho,
                "CLDN4_ImmuneScore_p": c.p,
                "CLDN4_ImmuneScore_partial_rho": c.rho_adj_epithelial,
                "CLDN4_ImmuneScore_partial_p": c.p_adj_epithelial,
                "verdict_CD8A": a.verdict,
                "verdict_CD274": b.verdict,
                "verdict_ImmuneScore": c.verdict,
            }
        )
    pd.DataFrame(one).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    c4_cd8 = pick_row(stats_df, "LUAD", "CD8A", "array")
    c4_pdl1 = pick_row(stats_df, "LUAD", "CD274", "array")
    c4_imm = pick_row(stats_df, "LUAD", "ImmuneScore", "array")

    scatter(
        FIGURES / "fig1_cldn4_vs_cd8a.png",
        luad["CLDN4"],
        luad["CD8A"],
        "CLDN4 (RMA-PLM log2)",
        "CD8A (RMA-PLM log2)",
        f"GSE33532 LUAD CLDN4 vs CD8A ({n_luad_arr} arrays / {n_luad_pt} patients)",
        int(c4_cd8.n) if c4_cd8 is not None else n_luad_arr,
        c4_cd8.rho if c4_cd8 is not None else np.nan,
        c4_cd8.p if c4_cd8 is not None else np.nan,
        patients=luad["patient"],
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_cd274.png",
        luad["CLDN4"],
        luad["CD274"],
        "CLDN4 (RMA-PLM log2)",
        "CD274 (RMA-PLM log2)",
        f"GSE33532 LUAD CLDN4 vs CD274 ({n_luad_arr} arrays / {n_luad_pt} patients)",
        int(c4_pdl1.n) if c4_pdl1 is not None else n_luad_arr,
        c4_pdl1.rho if c4_pdl1 is not None else np.nan,
        c4_pdl1.p if c4_pdl1 is not None else np.nan,
        patients=luad["patient"],
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_immunescore.png",
        luad["CLDN4"],
        luad["ImmuneScore"],
        "CLDN4 (RMA-PLM log2)",
        "ESTIMATE ImmuneScore",
        f"GSE33532 LUAD CLDN4 vs ImmuneScore ({n_luad_arr} arrays / {n_luad_pt} patients)",
        int(c4_imm.n) if c4_imm is not None else n_luad_arr,
        c4_imm.rho if c4_imm is not None else np.nan,
        c4_imm.p if c4_imm is not None else np.nan,
        patients=luad["patient"],
    )

    forest_keys = [
        ("LUAD", "patient", "CD8A", "LUAD patient CLDN4–CD8A"),
        ("LUAD", "patient", "CD274", "LUAD patient CLDN4–CD274"),
        ("LUAD", "patient", "ImmuneScore", "LUAD patient CLDN4–ImmuneScore"),
        ("LUAD", "array", "CD8A", "LUAD array CLDN4–CD8A"),
        ("LUAD", "array", "CD274", "LUAD array CLDN4–CD274"),
        ("LUAD", "array", "ImmuneScore", "LUAD array CLDN4–ImmuneScore"),
    ]
    labels_f, rhos, los, his = [], [], [], []
    for cohort, unit, end, lab in forest_keys:
        r = pick_row(stats_df, cohort, end, unit)
        if r is None or not np.isfinite(r.rho):
            continue
        labels_f.append(f"{lab} (n={int(r.n)})")
        rhos.append(r.rho)
        los.append(r.ci_lo)
        his.append(r.ci_hi)
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    y = np.arange(len(labels_f))[::-1]
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
    ax.set_yticklabels(labels_f)
    ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
    ax.set_title("GSE33532 LUAD CLDN4 vs CD8A / CD274 / ImmuneScore")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_spearman_forest.png", dpi=140)
    fig.savefig(FIGURES / "fig4_spearman_forest.pdf")
    plt.close(fig)

    write_finding(stats_df, cov, probe_used, epi_present, dict(counts))

    p_cd8 = pick_row(stats_df, "LUAD", "CD8A", "patient")
    p_pdl1 = pick_row(stats_df, "LUAD", "CD274", "patient")
    p_imm = pick_row(stats_df, "LUAD", "ImmuneScore", "patient")
    summary = {
        "dataset": "GSE33532",
        "primary_cohort": "LUAD tumor, patient-mean of 4 sites",
        "n_luad_arrays": n_luad_arr,
        "n_luad_patients": n_luad_pt,
        "n_named_complete": n_named,
        "CLDN4_CD8A_patient": {
            "n": None if p_cd8 is None else int(p_cd8.n),
            "rho": None if p_cd8 is None else float(p_cd8.rho),
            "p": None if p_cd8 is None else float(p_cd8.p),
            "partial_rho": None if p_cd8 is None else float(p_cd8.rho_adj_epithelial),
            "partial_p": None if p_cd8 is None else float(p_cd8.p_adj_epithelial),
            "verdict": None if p_cd8 is None else p_cd8.verdict,
        },
        "CLDN4_CD274_patient": {
            "n": None if p_pdl1 is None else int(p_pdl1.n),
            "rho": None if p_pdl1 is None else float(p_pdl1.rho),
            "p": None if p_pdl1 is None else float(p_pdl1.p),
            "partial_rho": None if p_pdl1 is None else float(p_pdl1.rho_adj_epithelial),
            "partial_p": None if p_pdl1 is None else float(p_pdl1.p_adj_epithelial),
            "verdict": None if p_pdl1 is None else p_pdl1.verdict,
        },
        "CLDN4_ImmuneScore_patient": {
            "n": None if p_imm is None else int(p_imm.n),
            "rho": None if p_imm is None else float(p_imm.rho),
            "p": None if p_imm is None else float(p_imm.p),
            "partial_rho": None if p_imm is None else float(p_imm.rho_adj_epithelial),
            "partial_p": None if p_imm is None else float(p_imm.p_adj_epithelial),
            "verdict": None if p_imm is None else p_imm.verdict,
        },
        "probes": probe_used,
        "histo_counts": dict(counts),
        "epithelial_genes": epi_present,
        "estimate_immune_present": cov["estimate_immune_present"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    print("wrote", HERE / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
