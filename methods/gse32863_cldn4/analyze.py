#!/usr/bin/env python3
"""GSE32863 LUAD Illumina WG-6 v3.0: CLDN4 vs CD8A / CD274 / ImmuneScore.

Additive CLDN4-only public series-matrix slice. No dual-high. No TACSTD2 claim.

Honest n is the deposited tumor arrays after matched adjacent lung is dropped.
GEO series summary text says “60 tumors”; overall design and the matrix are
58 LUAD + 58 adjacent non-tumor (116 arrays). Do not write n=60, n=59
(methylation companion GSE32867), or n=116.

Downloads stay under $GSE32863_CLDN4_DATA (default /tmp/gse32863_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
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
DATA = Path(os.environ.get("GSE32863_CLDN4_DATA", "/tmp/gse32863_cldn4"))
SIG = HERE / "estimate_yoshihara_2013.tsv"
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32863/"
    "matrix/GSE32863_series_matrix.txt.gz"
)
GPL_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/"
    "annot/GPL6884.annot.gz"
)

UA = "sdaxcge-gse32863-cldn4/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"

# Named Illumina HumanWG-6 v3.0 probes (GEO GPL6884 annot). Collapse still uses max-mean.
NAMED_PROBES = {
    "CLDN4": "ILMN_2132458",
    "CD8A": "ILMN_1768482",  # max-mean will pick among the three CD8A probes
    "CD274": "ILMN_1701914",
}

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
PRIMARY = ["CD8A", "CD274", "ImmuneScore"]
TARGETS = ["CLDN4", "CD8A", "CD274"]


def dl(dest: Path, url: str) -> Path:
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


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_characteristics"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                char_rows.append(vals)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"

    # GEO puts mixed keys in later characteristic slots. Parse by prefix.
    n = len(meta)
    parsed: dict[str, list[str]] = {}
    for vals in char_rows:
        if len(vals) != n:
            continue
        for i, v in enumerate(vals):
            if not v or ":" not in v:
                continue
            key, val = v.split(":", 1)
            key = key.strip().lower()
            val = val.strip()
            parsed.setdefault(key, [""] * n)
            # keep first non-empty for this key
            if parsed[key][i] == "":
                parsed[key][i] = val
    for key, vals in parsed.items():
        col = re.sub(r"[^a-z0-9]+", "_", key).strip("_")
        meta[col] = vals
    return meta, {k: " | ".join(v) for k, v in series.items()}


def first_symbol(s: str) -> str:
    if not isinstance(s, str) or not s or s == "nan":
        return ""
    return s.split(" /// ")[0].strip()


def load_gpl_annot(path: Path) -> pd.Series:
    with gzip.open(path, "rt", errors="replace") as f:
        text = f.read()
    start = text.index("!platform_table_begin")
    end = text.index("!platform_table_end")
    body = text[start:].split("\n", 1)[1]
    body = body[: body.index("!platform_table_end")]
    df = pd.read_csv(io.StringIO(body), sep="\t", dtype=str)
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    if "Gene symbol" not in df.columns:
        raise SystemExit(f"no Gene symbol column in {path}: {list(df.columns)[:12]}")
    raw = df.set_index(id_col)["Gene symbol"].map(first_symbol)
    return raw[raw.str.len() > 0]


def collapse_maxmean(probe_expr: pd.DataFrame, id2gene: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = probe_expr.index.map(lambda x: id2gene.get(str(x), ""))
    df = probe_expr.copy()
    df["_gene"] = genes
    df = df[df["_gene"].astype(str).str.len() > 0]
    means = df.drop(columns="_gene").astype(float).mean(axis=1)
    df["_mean"] = means
    df = df.sort_values("_mean", ascending=False)
    kept = df.loc[~df["_gene"].duplicated(keep="first")]
    gene_expr = kept.drop(columns=["_gene", "_mean"]).astype(float)
    gene_expr.index = kept["_gene"].values
    gene_expr = gene_expr.sort_index()
    tmp = df.reset_index()
    probe_col = tmp.columns[0]
    probe_audit = (
        tmp.rename(columns={probe_col: "probe", "_gene": "gene", "_mean": "mean"})
        [["gene", "probe", "mean"]]
        .sort_values(["gene", "mean"], ascending=[True, False])
    )
    return gene_expr, probe_audit


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


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def partial_spearman(x, y, *covars):
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


def verdict(rho_adj, p_adj, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho_adj) and rho_adj < 0 and p_adj < 0.05:
        return "HOLDS"
    if np.isfinite(rho_adj) and rho_adj > 0 and p_adj < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


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


def zmean(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    sub = expr.loc[present].astype(float)
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0), present


def tissue_class(row) -> str:
    bits = " ".join(
        str(row.get(c, ""))
        for c in ("source_name_ch1", "tissue", "type", "title")
    ).lower()
    if "adenocarcinoma" in bits or "lung tumor" in bits:
        return "tumor"
    if "non-tumor" in bits or "nontumor" in bits or "normal lung" in bits or "adjacent" in bits:
        return "normal"
    return "other"


def patient_id(title: str) -> str:
    m = re.search(r"(\d{2}L\d+[A-Z]?)_([NT])", str(title))
    if m:
        return m.group(1)
    m = re.search(r"(\d{2}L\d+)", str(title))
    return m.group(1) if m else ""


def pair_row(cldn4, y, epithelial, subset, partner, source):
    crude = spearman_ci(cldn4.values, y.values)
    adj_epi = partial_spearman(cldn4.values, y.values, epithelial.values)
    return {
        "subset": subset,
        "partner": partner,
        "source": source,
        "n": crude["n"],
        "rho": crude["rho"],
        "p": crude["p"],
        "ci_low": crude["ci_low"],
        "ci_high": crude["ci_high"],
        "rho_adj_epithelial": adj_epi["rho_adj"],
        "p_adj_epithelial": adj_epi["p_adj"],
        "verdict_epithelial": verdict(adj_epi["rho_adj"], adj_epi["p_adj"], adj_epi["n"]),
    }


def pick(corr: pd.DataFrame, subset, partner, source="collapsed_maxmean"):
    hit = corr[(corr.subset == subset) & (corr.partner == partner) & (corr.source == source)]
    return hit.iloc[0].to_dict() if len(hit) else {}


def md_pair(row) -> str:
    return (
        f"{fmt_rho(row.get('rho'))} ({fmt_p(row.get('p'))})"
        if row
        else "NA"
    )


def md_adj(row) -> str:
    return (
        f"{fmt_rho(row.get('rho_adj_epithelial'))} ({fmt_p(row.get('p_adj_epithelial'))})"
        if row
        else "NA"
    )


def write_finding(corr, cov, probe_used, epi_used, immune_n, stromal_n, meta_t, series):
    r8 = pick(corr, "tumor_LUAD", "CD8A")
    r274 = pick(corr, "tumor_LUAD", "CD274")
    rimm = pick(corr, "tumor_LUAD", "ImmuneScore")
    repi = pick(corr, "tumor_LUAD", "epithelial_z", "context")
    rctrl = pick(corr, "tumor_LUAD", "ImmuneScore", "CD8A_vs_ImmuneScore_control")
    n = int(cov["n_tumor"])

    smoking = Counter(meta_t.get("smoking_status", pd.Series(dtype=str)).fillna(""))
    stage = Counter(meta_t.get("stage", pd.Series(dtype=str)).fillna(""))
    gender = Counter(meta_t.get("gender", pd.Series(dtype=str)).fillna(""))
    eth = Counter(meta_t.get("ethnicity", pd.Series(dtype=str)).fillna(""))
    kras = Counter(meta_t.get("kras_mutation_status", pd.Series(dtype=str)).fillna(""))
    egfr = Counter(meta_t.get("egfr_mutation_status", pd.Series(dtype=str)).fillna(""))
    source = Counter(meta_t.get("tissue_source", pd.Series(dtype=str)).fillna(""))
    batch = Counter(meta_t.get("batch", pd.Series(dtype=str)).fillna(""))

    def cnt(c: Counter) -> str:
        return ", ".join(f"{k or 'NA'} {v}" for k, v in c.most_common())

    text = f"""# GSE32863 LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** Public Selamat / Laird-Offringa resected LUAD Illumina HumanWG-6 v3.0 series (Selamat et al., *Genome Res* 2012, PMID 22613842; GEO [GSE32863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32863); GPL6884). Unit is the **tumor array**. Matched adjacent non-tumor arrays are dropped. No slide was re-scored. No ICI arm. No dual-high. TACSTD2 is not a claim.

Primary question: **CLDN4 vs CD8A**, **CLDN4 vs CD274**, and **CLDN4 vs ESTIMATE ImmuneScore**, with an **epithelial residual**.

## Honest n

Do **not** write n=60. The GEO series *summary* says “60 lung adenocarcinoma tumors”; the *overall design* and the deposited series matrix are **58 LUAD + 58 adjacent non-tumor**. Do **not** write n=59 (that is the methylation companion [GSE32867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32867) / paper Infinium 27k set). Do **not** write n=116 (that pool includes matched normal).

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **116** | 48,803 probes × 116 GSM; log2 + lumi RSN as deposited |
| Unique GSM / unique titles | yes | **116** | all unique |
| Adjacent non-tumor lung | yes | **58** | `source_name_ch1` / `tissue: Normal lung`; **dropped** |
| **Tumor LUAD arrays (primary n)** | yes | **{n}** | `source_name_ch1` = Lung adenocarcinoma; `tissue: Lung tumor` |
| Unique patients (title `xxL*_T`) | yes | **{int(cov['n_patients'])}** | one tumor array per patient |
| Matched pairs (tumor + adjacent) | yes | **{int(cov['n_pairs'])}** | pairing is inventory only; tests are tumor-only |
| Stage | yes | {n} | deposited; **not** the claim |
| Smoking / pack-years | yes | {n} | deposited; **not** the claim |
| KRAS / EGFR / LKB1 | yes (tumor slot) | {n} | deposited; **not** the claim |
| Recurrence yes/no | yes | {n} | deposited; **not** a survival claim (no time) |
| OS / DSS time or event | no | 0 | not a GEO characteristic |
| ICI / treatment | no | 0 | 2012 methylation/expression atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | 0 | protocol is macrodissection; values not deposited |
| CLDN4 finite (`{probe_used['CLDN4']}`) | yes | **{n}** | single GPL6884 probe; max-mean collapse |
| CD8A finite (max-mean of 3 probes) | yes | **{n}** | collapse chose `{probe_used['CD8A']}` |
| CD274 finite (`{probe_used['CD274']}`) | yes | **{n}** | single GPL6884 probe |
| Epithelial mean-z | yes | **{n}** | {len(epi_used)}/{len(EPITHELIAL)} of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |
| ESTIMATE ImmuneScore | computed | **{n}** | Yoshihara Immune141 ssGSEA ({immune_n}/141 genes) |
| **Primary pairwise n** | yes | **{n}** | complete-case tumor CLDN4 + CD8A + CD274 + ImmuneScore |

The computable public tumor n is **{n} LUAD arrays**. That is the n in the one-row table.

Tumor-only clinical inventory (not tested): smoking {cnt(smoking)}; stage {cnt(stage)}; gender {cnt(gender)}; ethnicity {cnt(eth)}; KRAS {cnt(kras)}; EGFR {cnt(egfr)}; tissue source {cnt(source)}; expression batch {cnt(batch)}.

## One-row table

| dataset | histology | platform | n tumor | CLDN4–CD8A ρ (p) | ρ\\|epi (p) | CLDN4–CD274 ρ (p) | ρ\\|epi (p) | CLDN4–ImmuneScore ρ (p) | ρ\\|epi (p) | verdict |
|---|---|---|---:|---|---|---|---|---|---|---|
| GSE32863 Selamat | LUAD (normals dropped) | GPL6884 HumanWG-6 v3.0 log2-RSN | **{n}** | **{md_pair(r8)}** | {md_adj(r8)} | **{md_pair(r274)}** | {md_adj(r274)} | **{md_pair(rimm)}** | {md_adj(rimm)} | **{r8.get('verdict_epithelial', 'NA')} / {r274.get('verdict_epithelial', 'NA')} / {rimm.get('verdict_epithelial', 'NA')}** |

Full numbers: `tables/one_row.tsv`, `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary)

Deposited log2 + Robust Spline Normalization (lumi). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z ({'/'.join(epi_used)}; **{len(epi_used)}/{len(EPITHELIAL)}** present). HOLDS rule (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05. No dual-high cut.

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | verdict |
|---|---|---:|---:|---|---:|---|---|
| CLDN4 vs **CD8A** | tumor LUAD | **{n}** | **{fmt_rho(r8.get('rho'))}** | {fmt_ci(r8.get('ci_low'), r8.get('ci_high'))} | {fmt_p(r8.get('p'))} | {md_adj(r8)} | **{r8.get('verdict_epithelial', 'NA')}** |
| CLDN4 vs **CD274** | tumor LUAD | **{n}** | **{fmt_rho(r274.get('rho'))}** | {fmt_ci(r274.get('ci_low'), r274.get('ci_high'))} | {fmt_p(r274.get('p'))} | {md_adj(r274)} | **{r274.get('verdict_epithelial', 'NA')}** |
| CLDN4 vs **ImmuneScore** | tumor LUAD | **{n}** | **{fmt_rho(rimm.get('rho'))}** | {fmt_ci(rimm.get('ci_low'), rimm.get('ci_high'))} | {fmt_p(rimm.get('p'))} | {md_adj(rimm)} | **{rimm.get('verdict_epithelial', 'NA')}** |
| CLDN4 vs epithelial mean-z | tumor LUAD | {n} | {fmt_rho(repi.get('rho'))} | {fmt_ci(repi.get('ci_low'), repi.get('ci_high'))} | {fmt_p(repi.get('p'))} | — | context |
| CD8A vs ImmuneScore | tumor LUAD | {n} | {fmt_rho(rctrl.get('rho'))} | {fmt_ci(rctrl.get('ci_low'), rctrl.get('ci_high'))} | {fmt_p(rctrl.get('p'))} | — | positive control |

ImmuneScore is Yoshihara 2013 Immune141 ssGSEA (Barbie/GSVA, τ=0.25; {immune_n}/141 genes present after first-symbol max-mean collapse). Stromal141 coverage is {stromal_n}/141 (not used in the residual). There is no deposited numeric purity; the residual is the **epithelial mean-z**, not ESTIMATE TumorPurity.

Named-probe sensitivity (CLDN4 `{NAMED_PROBES['CLDN4']}` vs CD8A `{NAMED_PROBES['CD8A']}` / CD274 `{NAMED_PROBES['CD274']}`) is in `tables/spearman.tsv`.

### Q4 vs Q1 (descriptive)

See `tables/highlow_cldn4.tsv`. Quartiles are descriptive only and are not a dual-high claim.

## Matrix and labels (nothing invented)

| field | public? | n | what is there |
|---|---|---:|---|
| tissue / source | yes | 116 | 58 tumor / 58 adjacent; tests use tumor only |
| histology | yes | {n} | LUAD only on the tumor arrays |
| ICI response | **no** | 0 | surgical methylation/expression series |
| OS time | **no** | 0 | recurrence yes/no is deposited; unused |
| ESTIMATE published scores | **no** | 0 | ImmuneScore computed here from Yoshihara lists |

Max-mean unique-symbol GPL6884 collapse for the named genes:

| gene | probe used |
|---|---|
| CLDN4 | `{probe_used['CLDN4']}` |
| CD8A | `{probe_used['CD8A']}` |
| CD274 | `{probe_used['CD274']}` |

## What is not done

- No dual-high (no TACSTD2∩CLDN4 high intersection).
- No pooling of matched normal into the immune Spearman.
- No re-use of the series-summary “60 tumors” or the methylation-paper n=59.
- No ICI ORR / PFS model (labels are not deposited).
- No OS model (no time-to-event).

## Files

- `analyze.py` — GEO download, tumor filter, epithelial residual, ImmuneScore, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `spearman.tsv`, `label_inventory.tsv`, `per_sample.tsv`, `probe_confirm.tsv`, `gene_coverage.tsv`, `highlow_cldn4.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse32863_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)
    return text


def scatter(path: Path, x, y, xlab, ylab, title, n, rho, p):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(x, y, s=26, alpha=0.82, c="#2c5f8a", edgecolors="none")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nn={n}  ρ={fmt_rho(rho)}  p={fmt_p(p)}")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE32863_series_matrix.txt.gz", GEO_MATRIX)
    annot_path = dl(DATA / "GPL6884.annot.gz", GPL_ANNOT)

    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        raw = fh.read()
    meta, series = parse_series_meta(matrix_path)
    lines = raw.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    probe_expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    probe_expr.index = probe_expr.index.astype(str).str.strip('"')
    probe_expr.columns = [c.strip('"') for c in probe_expr.columns]
    probe_expr = probe_expr.astype(float)
    probe_expr = probe_expr.loc[:, [c for c in probe_expr.columns if c in meta.index]]
    meta = meta.loc[probe_expr.columns].copy()
    meta["tissue_class"] = meta.apply(tissue_class, axis=1)
    meta["patient"] = meta["title"].map(patient_id)

    counts = Counter(meta["tissue_class"])
    if counts.get("other", 0):
        raise SystemExit(f"unexpected tissue_class: {dict(counts)}")
    n_matrix = int(meta.shape[0])
    n_tumor = int((meta["tissue_class"] == "tumor").sum())
    n_normal = int((meta["tissue_class"] == "normal").sum())
    if n_tumor != 58 or n_normal != 58 or n_matrix != 116:
        raise SystemExit(
            f"honest-n check failed: matrix={n_matrix} tumor={n_tumor} normal={n_normal}; "
            "expected 116 / 58 / 58. See label_inventory."
        )

    tumor_idx = meta.index[meta["tissue_class"] == "tumor"]
    meta_t = meta.loc[tumor_idx].copy()
    n_patients = int(meta_t["patient"].nunique())
    pair_patients = set(meta.loc[meta["tissue_class"] == "tumor", "patient"]) & set(
        meta.loc[meta["tissue_class"] == "normal", "patient"]
    )
    n_pairs = len(pair_patients)
    if n_patients != n_tumor:
        print(f"WARN unique tumor patients={n_patients} vs tumor arrays={n_tumor}", flush=True)

    id2gene = load_gpl_annot(annot_path)
    # Collapse on tumor arrays only (primary universe).
    probes_t = probe_expr.loc[:, tumor_idx]
    gene_expr, probe_audit = collapse_maxmean(probes_t, id2gene)

    miss = [g for g in TARGETS if g not in gene_expr.index]
    if miss:
        raise SystemExit(f"named genes missing after max-mean collapse: {miss}")

    probe_used = {}
    for g in TARGETS:
        probe_used[g] = str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0])

    sig = pd.read_csv(SIG, sep="\t")
    stromal = sig.loc[sig["set"].str.contains("Stromal", case=False), "hugo"].astype(str).tolist()
    immune = sig.loc[sig["set"].str.contains("Immune", case=False), "hugo"].astype(str).tolist()
    immune_score = ssgsea(gene_expr, immune)
    stromal_score = ssgsea(gene_expr, stromal)
    immune_n = int(sum(g in gene_expr.index for g in immune))
    stromal_n = int(sum(g in gene_expr.index for g in stromal))

    epithelial, epi_used = zmean(gene_expr, EPITHELIAL)
    if len(epi_used) < 4:
        raise SystemExit(f"epithelial panel too thin: {epi_used}")

    n_nan = int(probes_t.isna().sum().sum())
    n_named = int(
        (
            np.isfinite(gene_expr.loc["CLDN4"])
            & np.isfinite(gene_expr.loc["CD8A"])
            & np.isfinite(gene_expr.loc["CD274"])
            & np.isfinite(immune_score)
        ).sum()
    )
    if n_named != n_tumor:
        print(f"WARN complete-case named n={n_named} vs tumor n={n_tumor}", flush=True)

    per = pd.DataFrame(
        {
            "gsm": tumor_idx,
            "patient": meta_t["patient"].to_numpy(),
            "title": meta_t["title"].to_numpy(),
            "tissue_class": meta_t["tissue_class"].to_numpy(),
            "tissue": meta_t.get("tissue", pd.Series("", index=tumor_idx)).to_numpy(),
            "source_name": meta_t.get("source_name_ch1", pd.Series("", index=tumor_idx)).to_numpy(),
            "batch": meta_t.get("batch", pd.Series("", index=tumor_idx)).to_numpy(),
            "tissue_source": meta_t.get("tissue_source", pd.Series("", index=tumor_idx)).to_numpy(),
            "age": meta_t.get("age", pd.Series("", index=tumor_idx)).to_numpy(),
            "gender": meta_t.get("gender", pd.Series("", index=tumor_idx)).to_numpy(),
            "ethnicity": meta_t.get("ethnicity", pd.Series("", index=tumor_idx)).to_numpy(),
            "smoking_status": meta_t.get("smoking_status", pd.Series("", index=tumor_idx)).to_numpy(),
            "pack_years": meta_t.get("pack_years", pd.Series("", index=tumor_idx)).to_numpy(),
            "stage": meta_t.get("stage", pd.Series("", index=tumor_idx)).to_numpy(),
            "recurrence": meta_t.get("recurrence", pd.Series("", index=tumor_idx)).to_numpy(),
            "kras_mutation_status": meta_t.get("kras_mutation_status", pd.Series("", index=tumor_idx)).to_numpy(),
            "egfr_mutation_status": meta_t.get("egfr_mutation_status", pd.Series("", index=tumor_idx)).to_numpy(),
            "lkb1_mutation_status": meta_t.get("lkb1_mutation_status", pd.Series("", index=tumor_idx)).to_numpy(),
            "CLDN4": gene_expr.loc["CLDN4"].to_numpy(float),
            "CD8A": gene_expr.loc["CD8A"].to_numpy(float),
            "CD274": gene_expr.loc["CD274"].to_numpy(float),
            "ImmuneScore": immune_score.reindex(tumor_idx).to_numpy(float),
            "StromalScore": stromal_score.reindex(tumor_idx).to_numpy(float),
            "epithelial_z": epithelial.reindex(tumor_idx).to_numpy(float),
        },
        index=tumor_idx,
    )
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    cov = {
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_normal": n_normal,
        "n_patients": n_patients,
        "n_pairs": n_pairs,
        "n_named_complete": n_named,
        "n_probes": int(probe_expr.shape[0]),
        "n_genes_collapsed": int(gene_expr.shape[0]),
        "n_missing_tumor_values": n_nan,
        "estimate_immune_present": immune_n,
        "estimate_stromal_present": stromal_n,
        "epithelial_genes": ",".join(epi_used),
        "probe_CLDN4": probe_used["CLDN4"],
        "probe_CD8A": probe_used["CD8A"],
        "probe_CD274": probe_used["CD274"],
        "processing": "log2 + lumi Robust Spline Normalization (deposited)",
    }
    pd.Series(cov, name="value").to_csv(TABLES / "coverage.tsv", sep="\t")

    labels = pd.DataFrame(
        [
            {"item": "arrays in series matrix", "public": "yes", "n": n_matrix,
             "note": "do not use as the immune n; includes matched normal"},
            {"item": "adjacent non-tumor lung", "public": "yes", "n": n_normal,
             "note": "dropped; source_name_ch1 / tissue: Normal lung"},
            {"item": "tumor LUAD arrays (primary n)", "public": "yes", "n": n_tumor,
             "note": "honest tumor n; not the series-summary '60'"},
            {"item": "unique tumor patients", "public": "yes", "n": n_patients,
             "note": "one tumor array per title patient id"},
            {"item": "matched pairs", "public": "yes", "n": n_pairs,
             "note": "inventory only; Spearman is tumor-only"},
            {"item": "series-summary '60 tumors'", "public": "text only", "n": 60,
             "note": "DO NOT USE; contradicts overall design and matrix"},
            {"item": "methylation companion GSE32867 / paper 59 pairs", "public": "other series", "n": 59,
             "note": "DO NOT USE; different assay"},
            {"item": "ICI / OS time", "public": "no", "n": 0,
             "note": "recurrence yes/no is deposited but unused"},
            {"item": "numeric tumor percent", "public": "no", "n": 0,
             "note": "macrodissection mentioned; values not deposited"},
            {"item": "primary pairwise n", "public": "yes", "n": n_named,
             "note": "complete-case CLDN4 + CD8A + CD274 + ImmuneScore"},
        ]
    )
    labels.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    cov_rows = []
    for g in TARGETS + EPITHELIAL:
        cov_rows.append(
            {
                "gene": g,
                "present": int(g in gene_expr.index),
                "collapse_probe": (
                    str(probe_audit.loc[probe_audit["gene"] == g, "probe"].iloc[0])
                    if g in gene_expr.index
                    else ""
                ),
                "n_probes_mapped": int((probe_audit["gene"] == g).sum()) if g in set(probe_audit["gene"]) else 0,
            }
        )
    pd.DataFrame(cov_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"gene": g, "named_probe": NAMED_PROBES[g], "collapse_probe": probe_used[g]} for g in TARGETS]
    ).to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    rows = []
    cldn4 = per["CLDN4"]
    epi = per["epithelial_z"]
    for partner, series_y in (
        ("CD8A", per["CD8A"]),
        ("CD274", per["CD274"]),
        ("ImmuneScore", per["ImmuneScore"]),
    ):
        rows.append(pair_row(cldn4, series_y, epi, "tumor_LUAD", partner, "collapsed_maxmean"))

    # named-probe sensitivity
    rows.append(
        pair_row(
            probes_t.loc[NAMED_PROBES["CLDN4"]],
            probes_t.loc[NAMED_PROBES["CD8A"]],
            epi,
            "tumor_LUAD",
            "CD8A",
            f"named_{NAMED_PROBES['CLDN4']}_vs_{NAMED_PROBES['CD8A']}",
        )
    )
    rows.append(
        pair_row(
            probes_t.loc[NAMED_PROBES["CLDN4"]],
            probes_t.loc[NAMED_PROBES["CD274"]],
            epi,
            "tumor_LUAD",
            "CD274",
            f"named_{NAMED_PROBES['CLDN4']}_vs_{NAMED_PROBES['CD274']}",
        )
    )

    # context / controls
    crude_epi = spearman_ci(cldn4.values, epi.values)
    rows.append(
        {
            "subset": "tumor_LUAD",
            "partner": "epithelial_z",
            "source": "context",
            "n": crude_epi["n"],
            "rho": crude_epi["rho"],
            "p": crude_epi["p"],
            "ci_low": crude_epi["ci_low"],
            "ci_high": crude_epi["ci_high"],
            "rho_adj_epithelial": np.nan,
            "p_adj_epithelial": np.nan,
            "verdict_epithelial": "COMPANION",
        }
    )
    crude_ctrl = spearman_ci(per["CD8A"].values, per["ImmuneScore"].values)
    rows.append(
        {
            "subset": "tumor_LUAD",
            "partner": "ImmuneScore",
            "source": "CD8A_vs_ImmuneScore_control",
            "n": crude_ctrl["n"],
            "rho": crude_ctrl["rho"],
            "p": crude_ctrl["p"],
            "ci_low": crude_ctrl["ci_low"],
            "ci_high": crude_ctrl["ci_high"],
            "rho_adj_epithelial": np.nan,
            "p_adj_epithelial": np.nan,
            "verdict_epithelial": "CONTROL",
        }
    )
    crude_pdl1_cd8 = spearman_ci(per["CD274"].values, per["CD8A"].values)
    rows.append(
        {
            "subset": "tumor_LUAD",
            "partner": "CD8A",
            "source": "CD274_vs_CD8A_context",
            "n": crude_pdl1_cd8["n"],
            "rho": crude_pdl1_cd8["rho"],
            "p": crude_pdl1_cd8["p"],
            "ci_low": crude_pdl1_cd8["ci_low"],
            "ci_high": crude_pdl1_cd8["ci_high"],
            "rho_adj_epithelial": np.nan,
            "p_adj_epithelial": np.nan,
            "verdict_epithelial": "COMPANION",
        }
    )

    corr = pd.DataFrame(rows)
    corr.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)

    r8 = pick(corr, "tumor_LUAD", "CD8A")
    r274 = pick(corr, "tumor_LUAD", "CD274")
    rimm = pick(corr, "tumor_LUAD", "ImmuneScore")
    one = pd.DataFrame(
        [
            {
                "dataset": "GSE32863 Selamat LUAD",
                "histology": "LUAD",
                "platform": "GPL6884 Illumina HumanWG-6 v3.0 log2-RSN",
                "n_tumor": n_tumor,
                "n_normal_dropped": n_normal,
                "CLDN4_probe": probe_used["CLDN4"],
                "CD8A_probe": probe_used["CD8A"],
                "CD274_probe": probe_used["CD274"],
                "CLDN4_CD8A_rho": r8.get("rho"),
                "CLDN4_CD8A_p": r8.get("p"),
                "CLDN4_CD8A_partial_epi_rho": r8.get("rho_adj_epithelial"),
                "CLDN4_CD8A_partial_epi_p": r8.get("p_adj_epithelial"),
                "CLDN4_CD8A_verdict": r8.get("verdict_epithelial"),
                "CLDN4_CD274_rho": r274.get("rho"),
                "CLDN4_CD274_p": r274.get("p"),
                "CLDN4_CD274_partial_epi_rho": r274.get("rho_adj_epithelial"),
                "CLDN4_CD274_partial_epi_p": r274.get("p_adj_epithelial"),
                "CLDN4_CD274_verdict": r274.get("verdict_epithelial"),
                "CLDN4_ImmuneScore_rho": rimm.get("rho"),
                "CLDN4_ImmuneScore_p": rimm.get("p"),
                "CLDN4_ImmuneScore_partial_epi_rho": rimm.get("rho_adj_epithelial"),
                "CLDN4_ImmuneScore_partial_epi_p": rimm.get("p_adj_epithelial"),
                "CLDN4_ImmuneScore_verdict": rimm.get("verdict_epithelial"),
                "epithelial_genes": ",".join(epi_used),
                "Immune141_present": immune_n,
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    hl_rows = []
    q = pd.qcut(cldn4, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    if q.nunique() == 4:
        for partner in ("CD8A", "CD274", "ImmuneScore"):
            a = per[partner][q == "Q4"]
            b = per[partner][q == "Q1"]
            U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
            rbc = 2 * U / (a.size * b.size) - 1
            hl_rows.append(
                {
                    "subset": "tumor_LUAD",
                    "endpoint": partner,
                    "n_Q4": int(a.size),
                    "n_Q1": int(b.size),
                    "median_Q4": float(np.median(a)),
                    "median_Q1": float(np.median(b)),
                    "U": float(U),
                    "p": float(p_mw),
                    "rank_biserial_Q4_minus_Q1": float(rbc),
                }
            )
    pd.DataFrame(hl_rows).to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    scatter(
        FIGURES / "fig1_cldn4_vs_cd8a.png",
        per["CLDN4"],
        per["CD8A"],
        "CLDN4  log2-RSN",
        "CD8A  log2-RSN",
        "GSE32863 LUAD tumor  CLDN4 vs CD8A",
        int(r8.get("n", n_tumor)),
        r8.get("rho"),
        r8.get("p"),
    )
    scatter(
        FIGURES / "fig2_cldn4_vs_cd274.png",
        per["CLDN4"],
        per["CD274"],
        "CLDN4  log2-RSN",
        "CD274  log2-RSN",
        "GSE32863 LUAD tumor  CLDN4 vs CD274",
        int(r274.get("n", n_tumor)),
        r274.get("rho"),
        r274.get("p"),
    )
    scatter(
        FIGURES / "fig3_cldn4_vs_immunescore.png",
        per["CLDN4"],
        per["ImmuneScore"],
        "CLDN4  log2-RSN",
        "ESTIMATE ImmuneScore",
        "GSE32863 LUAD tumor  CLDN4 vs ImmuneScore",
        int(rimm.get("n", n_tumor)),
        rimm.get("rho"),
        rimm.get("p"),
    )

    plot = corr[(corr.source == "collapsed_maxmean") & (corr.partner.isin(PRIMARY))].copy()
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    y_pos = np.arange(len(plot))[::-1]
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
    ax.set_yticklabels([f"CLDN4 vs {r.partner}  n={int(r.n)}" for r in plot.itertuples()])
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    ax.set_title("GSE32863 LUAD tumor  CLDN4 vs CD8A / CD274 / ImmuneScore")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_spearman_forest.png", dpi=160)
    fig.savefig(FIGURES / "fig4_spearman_forest.pdf")
    plt.close(fig)

    write_finding(corr, cov, probe_used, epi_used, immune_n, stromal_n, meta_t, series)

    summary = {
        "dataset": "GSE32863",
        "pmid": "22613842",
        "platform": "GPL6884 Illumina HumanWG-6 v3.0",
        "processing": "log2 + lumi Robust Spline Normalization as deposited. Spearman is rank-based.",
        "claim": "CLDN4-only vs CD8A / CD274 / ImmuneScore; epithelial residual; tumor arrays only",
        "dual_high": False,
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_normal_dropped": n_normal,
        "n_patients": n_patients,
        "do_not_use_n": {
            "series_summary_60": 60,
            "methylation_GSE32867_59": 59,
            "matrix_including_normal_116": 116,
        },
        "epithelial_genes_used": epi_used,
        "Immune141_present": immune_n,
        "Stromal141_present": stromal_n,
        "probes": probe_used,
        "holds_rule": "n>=40, residual ρ<0, residual p<0.05 (same as PR 229 / GSE4573 / GSE10245)",
        "primary": {
            "CLDN4_CD8A": r8,
            "CLDN4_CD274": r274,
            "CLDN4_ImmuneScore": rimm,
        },
        "missing_public_labels": ["OS_time", "ICI"],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps({k: summary[k] for k in ("dataset", "n_tumor", "n_normal_dropped", "probes")}, indent=2))
    print(
        corr[corr.source.isin(["collapsed_maxmean", "context", "CD8A_vs_ImmuneScore_control"])][
            [
                "subset",
                "partner",
                "source",
                "n",
                "rho",
                "p",
                "rho_adj_epithelial",
                "p_adj_epithelial",
                "verdict_epithelial",
            ]
        ].to_string(index=False)
    )
    print("[done]", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
