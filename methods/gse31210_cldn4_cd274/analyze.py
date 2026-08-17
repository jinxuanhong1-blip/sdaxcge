#!/usr/bin/env python3
"""GSE31210 LUAD tumors: CLDN4 vs CD274 and HLA-A/B/C after ESTIMATE.

Additive only. New cut. CLDN4–CD8A ρ on this series is already known
(methods/gse31210_cldn4_immune) and is not re-cut here.

Public inputs:
  GEO GSE31210 series matrix (Okayama et al. Cancer Res 2012; GPL570)
  GPL570.annot gene symbols
  Yoshihara 2013 ESTIMATE Stromal141 + Immune141 (local TSV)

Unit = primary lung tumor array. Normals are dropped.
"""
from __future__ import annotations

import gzip
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
CACHE = Path("/tmp/gse31210_cldn4_cd274")
SIG = HERE / "estimate_yoshihara_2013.tsv"

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
    "matrix/GSE31210_series_matrix.txt.gz"
)
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

UA = (
    "sdaxcge-gse31210-cldn4-cd274/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

EST_A = 0.6049872018
EST_B = 0.0001467884

# New-cut named genes only. CD8A is taken as given and is not a target.
TARGETS = ["CLDN4", "CD274", "HLA-A", "HLA-B", "HLA-C"]
HLA_GENES = ["HLA-A", "HLA-B", "HLA-C"]


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


def partial_spearman(x, y, z):
    """Pearson of rank residuals; df = n − 3 (same as sibling GSE31210 cut)."""
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
    """Barbie/GSVA ssGSEA; ranks scaled 1..10000 (sibling / estimate-like)."""
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
    """Return (meta samples x characteristics, expression probe x sample)."""
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
    """Return (first-symbol map, unique-only map). Multi-mapped probes dropped in unique."""
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


def rec_row(predictor, endpoint, x, y, z, covar_name, note=""):
    ru, pu, nu = spearman_pair(x, y)
    rp, pp, np_ = partial_spearman(x, y, z) if z is not None else (np.nan, np.nan, nu)
    return {
        "predictor": predictor,
        "endpoint": endpoint,
        "n_unadjusted": nu,
        "unadj_rho": ru,
        "unadj_p": pu,
        "n_partial": np_,
        "partial_rho": rp,
        "partial_p": pp,
        "covariate": covar_name,
        "note": note,
    }


def pick_named_gene(genes_uniq, genes_first, probe_uniq, probe_first, gene: str):
    """Prefer unique-mapped max-mean; fall back to first-symbol and record why."""
    if gene in genes_uniq.index:
        return (
            genes_uniq.loc[gene].to_numpy(float),
            str(probe_uniq.loc[gene]),
            "unique-mapped max-mean",
        )
    if gene in genes_first.index:
        return (
            genes_first.loc[gene].to_numpy(float),
            str(probe_first.loc[gene]),
            "first-symbol max-mean (no unique-mapped probe)",
        )
    raise SystemExit(f"{gene} missing after unique and first-symbol collapse")


def write_finding(stats_df, cov, probe_used, probe_rule, n_tumor, n_matrix) -> None:
    def pick(pred, end, covar=None):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        if covar is not None:
            hit = hit[hit.covariate == covar]
        return hit.iloc[0] if len(hit) else None

    est_name = "ESTIMATEScore (Stromal+Immune)"
    c4_pd = pick("CLDN4", "CD274", est_name)
    c4_a = pick("CLDN4", "HLA-A", est_name)
    c4_b = pick("CLDN4", "HLA-B", est_name)
    c4_c = pick("CLDN4", "HLA-C", est_name)
    c4_mhc = pick("CLDN4", "MHC1_HLAabc", est_name)
    c4_imm = pick("CLDN4", "ImmuneScore", est_name)

    n_pur_ok = cov.get("n_purity_in_0_1", "NA")
    n_pur_out = cov.get("n_purity_outside_0_1", "NA")

    def cell(row, which="unadj"):
        if row is None:
            return "NA"
        if which == "unadj":
            return f"{fmt_rho(row.unadj_rho)} ({fmt_p(row.unadj_p)})"
        return f"{fmt_rho(row.partial_rho)} ({fmt_p(row.partial_p)})"

    md = f"""# GSE31210 — CLDN4 vs CD274 and HLA-A/B/C (partial on ESTIMATE)

**Additive only. New cut.** CLDN4–CD8A ρ on this series is **already known** (`methods/gse31210_cldn4_immune`) and is not re-cut. This folder only measures CLDN4 against **CD274** and **HLA-A / HLA-B / HLA-C**, with a rank residual on ESTIMATEScore.

Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570 Affymetrix U133 Plus 2.0). Unit is the **primary lung tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {n_matrix} | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **{n_tumor}** | `tissue: primary lung tumor` |
| adjacent / non-tumor | {n_matrix - n_tumor} | dropped |
| CLDN4 / CD274 / HLA-A / HLA-B / HLA-C | {n_tumor} | all five genes present after collapse |
| ESTIMATEScore | {n_tumor} | Yoshihara Stromal141+Immune141 ssGSEA ({cov['estimate_stromal_present']}/141 + {cov['estimate_immune_present']}/141) |
| ESTIMATE TumorPurity in [0, 1] | {n_pur_ok} | cosine map; {n_pur_out} tumors wrap outside [0, 1] |

Primary tests use **n={n_tumor}** tumors. That is the same tumor rule as the sibling CD8 cut; it is not a new cohort hunt.

## One-row table

| dataset | n | CLDN4–CD274 ρ (p) | CLDN4–HLA-A ρ (p) | CLDN4–HLA-B ρ (p) | CLDN4–HLA-C ρ (p) | CLDN4–CD274 partial \\| ESTIMATE (p) | CLDN4–HLA-A partial \\| ESTIMATE (p) | CLDN4–HLA-B partial \\| ESTIMATE (p) | CLDN4–HLA-C partial \\| ESTIMATE (p) |
|---|---:|---|---|---|---|---|---|---|---|
| GSE31210 LUAD tumors | {n_tumor} | {cell(c4_pd)} | {cell(c4_a)} | {cell(c4_b)} | {cell(c4_c)} | {cell(c4_pd, "partial")} | {cell(c4_a, "partial")} | {cell(c4_b, "partial")} | {cell(c4_c, "partial")} |

MHC-I mean-z (HLA-A/B/C, {cov.get('mhc1_genes_used', 'NA')}): unadj {cell(c4_mhc)}; partial \\| ESTIMATE {cell(c4_mhc, "partial")}.

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use the native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| tissue | yes | {n_matrix} | primary lung tumor {n_tumor} / other {n_matrix - n_tumor} |
| histology | yes | {n_tumor} | LUAD (series is stage I–II lung adenocarcinoma) |
| ICI response | **no** | 0 | treatment-naive surgical cohort |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Named-gene collapse (unique-mapped max-mean when a unique probe exists; otherwise first-symbol max-mean):

| gene | probe | rule |
|---|---|---|
| CLDN4 | `{probe_used.get('CLDN4', 'NA')}` | {probe_rule.get('CLDN4', 'NA')} |
| CD274 | `{probe_used.get('CD274', 'NA')}` | {probe_rule.get('CD274', 'NA')} |
| HLA-A | `{probe_used.get('HLA-A', 'NA')}` | {probe_rule.get('HLA-A', 'NA')} |
| HLA-B | `{probe_used.get('HLA-B', 'NA')}` | {probe_rule.get('HLA-B', 'NA')} |
| HLA-C | `{probe_used.get('HLA-C', 'NA')}` | {probe_rule.get('HLA-C', 'NA')} |

ESTIMATE gene collapse uses first-symbol mapping (same as the sibling cut).

## ESTIMATE

R `estimate` is not required. Scores follow the sibling public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

**Partial covariate is ESTIMATEScore**, not ImmuneScore. ImmuneScore is one addend of ESTIMATEScore and is collinear with it. On this matrix TumorPurity remains monotone with ESTIMATEScore (angles still in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The cosine is below 0 for **{n_pur_out} / {n_tumor}** tumors and is **not** used as a 0–1 fraction. The [0, 1] subset (n={n_pur_ok}) is only the least-impure tail and is **not** the analysis n.

CLDN4 vs ImmuneScore is reported only as a context row (same matrix as the sibling cut). It is not a new CD8 test.

## Main Spearman (n={n_tumor})

Partial = Pearson of rank residuals on ESTIMATEScore; df = n − 3.

| predictor | endpoint | n | ρ | p | partial ρ \\| ESTIMATE | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | CD274 | {int(c4_pd.n_unadjusted)} | **{fmt_rho(c4_pd.unadj_rho)}** | {fmt_p(c4_pd.unadj_p)} | {fmt_rho(c4_pd.partial_rho)} | {fmt_p(c4_pd.partial_p)} |
| CLDN4 | HLA-A | {int(c4_a.n_unadjusted)} | **{fmt_rho(c4_a.unadj_rho)}** | {fmt_p(c4_a.unadj_p)} | {fmt_rho(c4_a.partial_rho)} | {fmt_p(c4_a.partial_p)} |
| CLDN4 | HLA-B | {int(c4_b.n_unadjusted)} | **{fmt_rho(c4_b.unadj_rho)}** | {fmt_p(c4_b.unadj_p)} | {fmt_rho(c4_b.partial_rho)} | {fmt_p(c4_b.partial_p)} |
| CLDN4 | HLA-C | {int(c4_c.n_unadjusted)} | **{fmt_rho(c4_c.unadj_rho)}** | {fmt_p(c4_c.unadj_p)} | {fmt_rho(c4_c.partial_rho)} | {fmt_p(c4_c.partial_p)} |
| CLDN4 | MHC-I (HLA-A/B/C mean-z) | {int(c4_mhc.n_unadjusted)} | {fmt_rho(c4_mhc.unadj_rho)} | {fmt_p(c4_mhc.unadj_p)} | {fmt_rho(c4_mhc.partial_rho)} | {fmt_p(c4_mhc.partial_p)} |
| CLDN4 | ImmuneScore | {int(c4_imm.n_unadjusted)} | {fmt_rho(c4_imm.unadj_rho)} | {fmt_p(c4_imm.unadj_p)} | {fmt_rho(c4_imm.partial_rho)} | {fmt_p(c4_imm.partial_p)} |

This is an antigen-presentation / PD-L1 cut on the same n={n_tumor} LUAD tumor matrix. It is not an ICI-response test. CLDN4–CD8A is not restated.

## What is not done

- No re-cut of CLDN4 vs CD8A / TACSTD2 (sibling folder).
- No OncoSG or CPTAC re-run.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png`
- `figures/fig2_cldn4_vs_hla.png`
- `figures/fig3_spearman_forest.png`
- `figures/fig4_correlation_heatmap.png`

```bash
python3 methods/gse31210_cldn4_cd274/analyze.py
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

    matrix_path = download(GSE_MATRIX, CACHE / "GSE31210_series_matrix.txt.gz")
    annot_path = download(GPL570_ANNOT, CACHE / "GPL570.annot.gz")

    print("parse matrix", flush=True)
    meta, probes = parse_geo_matrix(matrix_path)
    first_map, unique_map = load_gpl_annot(annot_path)

    if "tissue" not in meta.columns:
        raise SystemExit(f"GSE31210 tissue characteristic missing: {list(meta.columns)}")
    tissue = meta["tissue"].astype(str)
    is_tumor = tissue.str.lower().str.contains("primary lung tumor")
    n_matrix = int(len(meta))
    n_tumor = int(is_tumor.sum())
    print(f"arrays={n_matrix} tumors={n_tumor}", flush=True)
    if n_tumor != 226:
        print(f"WARNING: tumor n={n_tumor} (catalog n=226)", flush=True)

    probes_t = probes.loc[:, meta.index[is_tumor]]
    genes_est, _ = collapse_maxmean(probes_t, first_map)
    genes_uniq, probe_uniq = collapse_maxmean(probes_t, unique_map)
    _, probe_first = collapse_maxmean(probes_t, first_map)

    probe_used = {}
    probe_rule = {}
    vecs = {}
    for g in TARGETS:
        v, pid, rule = pick_named_gene(genes_uniq, genes_est, probe_uniq, probe_first, g)
        vecs[g] = v
        probe_used[g] = pid
        probe_rule[g] = rule
        print(f"{g}: {pid} ({rule})", flush=True)

    pd.DataFrame(
        [{"gene": g, "probe_id": probe_used[g], "rule": probe_rule[g]} for g in TARGETS]
    ).to_csv(TABLES / "probe_used.tsv", sep="\t", index=False)

    est = estimate_scores(genes_est, stromal, immune)
    log_parts = {}
    for g in TARGETS:
        raw = np.clip(vecs[g], a_min=0, a_max=None)
        log_parts[g] = np.log2(raw + 1)
    logx = pd.DataFrame(log_parts, index=probes_t.columns)

    cldn4 = vecs["CLDN4"]
    cd274 = vecs["CD274"]
    hla_a = vecs["HLA-A"]
    hla_b = vecs["HLA-B"]
    hla_c = vecs["HLA-C"]
    imm = est["ESTIMATE_ImmuneScore"].to_numpy(float)
    strom = est["ESTIMATE_StromalScore"].to_numpy(float)
    estscore = est["ESTIMATE_Score"].to_numpy(float)
    pur = est["ESTIMATE_TumorPurity"].to_numpy(float)
    in01 = (pur >= 0) & (pur <= 1)

    hla_z = []
    hla_used = []
    for name, v in (("HLA-A", hla_a), ("HLA-B", hla_b), ("HLA-C", hla_c)):
        s = pd.Series(v, dtype=float)
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        hla_z.append(((s - s.mean()) / sd).to_numpy(float))
        hla_used.append(name)
    mhc = np.mean(np.vstack(hla_z), axis=0) if hla_z else np.full(n_tumor, np.nan)

    rows = []
    est_name = "ESTIMATEScore (Stromal+Immune)"
    pur_name = "ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141)"
    note_eq = "rank-equivalent to TumorPurity residual on this cohort"

    for end, y in (
        ("CD274", cd274),
        ("HLA-A", hla_a),
        ("HLA-B", hla_b),
        ("HLA-C", hla_c),
        ("MHC1_HLAabc", mhc),
    ):
        rows.append(rec_row("CLDN4", end, cldn4, y, estscore, est_name, note=note_eq))
        rows.append(rec_row("CLDN4", end, cldn4, y, pur, pur_name,
                            note="rank-equivalent to ESTIMATEScore residual on this cohort"))

    rows.append(rec_row("CLDN4", "ImmuneScore", cldn4, imm, estscore, est_name,
                        note="context only; ImmuneScore is an addend of ESTIMATEScore; not a CD8 re-cut"))
    rows.append(rec_row("CLDN4", "StromalScore", cldn4, strom, None, ""))
    rows.append(rec_row("CLDN4", "TumorPurity", cldn4, pur, None, ""))
    rows.append(rec_row("CLDN4", "ESTIMATEScore", cldn4, estscore, None, ""))
    rows.append(rec_row("CD274", "ImmuneScore", cd274, imm, None, "", note="context"))
    rows.append(rec_row("MHC1_HLAabc", "ImmuneScore", mhc, imm, None, "", note="context"))
    rows.append(rec_row("CD274", "MHC1_HLAabc", cd274, mhc, estscore, est_name, note="context"))

    if int(in01.sum()) >= 8:
        for end, y in (("CD274", cd274), ("HLA-A", hla_a), ("HLA-B", hla_b), ("HLA-C", hla_c)):
            rows.append(rec_row(
                "CLDN4", end, cldn4[in01], y[in01], estscore[in01],
                "ESTIMATEScore; TumorPurity in [0,1]",
                note=f"sensitivity; n={int(in01.sum())} of {n_tumor}",
            ))

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    per = pd.DataFrame({
        "sample": probes_t.columns,
        "tissue": tissue.loc[probes_t.columns].to_numpy(),
        "CLDN4": cldn4,
        "CD274": cd274,
        "HLA-A": hla_a,
        "HLA-B": hla_b,
        "HLA-C": hla_c,
        "CLDN4_log2": logx["CLDN4"].to_numpy(float),
        "CD274_log2": logx["CD274"].to_numpy(float),
        "HLA-A_log2": logx["HLA-A"].to_numpy(float),
        "HLA-B_log2": logx["HLA-B"].to_numpy(float),
        "HLA-C_log2": logx["HLA-C"].to_numpy(float),
        "MHC1_HLAabc_meanz": mhc,
        "ESTIMATE_StromalScore": strom,
        "ESTIMATE_ImmuneScore": imm,
        "ESTIMATE_Score": estscore,
        "ESTIMATE_TumorPurity": pur,
        "purity_in_0_1": in01,
    })
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    cov = {
        "accession": "GSE31210",
        "platform": "GPL570 U133 Plus 2.0",
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "tumor_rule": "tissue == primary lung tumor",
        "n_genes_unique": int(genes_uniq.shape[0]),
        "n_genes_firstsymbol": int(genes_est.shape[0]),
        "targets_present": ",".join(TARGETS),
        "estimate_stromal_present": int(sum(g in genes_est.index for g in stromal)),
        "estimate_immune_present": int(sum(g in genes_est.index for g in immune)),
        "mhc1_genes_used": ",".join(hla_used),
        "n_purity_in_0_1": int(in01.sum()),
        "n_purity_outside_0_1": int((~in01).sum()),
        "unit_plot": "log2(MAS5+1)",
        "unit_estimate": "ssGSEA on MAS5 linear ranks (rank-invariant)",
        "partial_covariate": est_name,
        "cd8_cut": "taken_as_given_sibling_gse31210_cldn4_immune",
    }
    pd.DataFrame([cov]).to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    def g(pred, end, col, covar=est_name):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)
                       & (stats_df.covariate == covar)]
        if hit.empty:
            hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        return hit.iloc[0][col] if len(hit) else np.nan

    one = pd.DataFrame([{
        "dataset": "GSE31210 LUAD tumors",
        "n": n_tumor,
        "CLDN4_probe": probe_used["CLDN4"],
        "CD274_probe": probe_used["CD274"],
        "HLA_A_probe": probe_used["HLA-A"],
        "HLA_B_probe": probe_used["HLA-B"],
        "HLA_C_probe": probe_used["HLA-C"],
        "CLDN4_vs_CD274_rho": g("CLDN4", "CD274", "unadj_rho"),
        "CLDN4_vs_CD274_p": g("CLDN4", "CD274", "unadj_p"),
        "CLDN4_vs_HLA_A_rho": g("CLDN4", "HLA-A", "unadj_rho"),
        "CLDN4_vs_HLA_A_p": g("CLDN4", "HLA-A", "unadj_p"),
        "CLDN4_vs_HLA_B_rho": g("CLDN4", "HLA-B", "unadj_rho"),
        "CLDN4_vs_HLA_B_p": g("CLDN4", "HLA-B", "unadj_p"),
        "CLDN4_vs_HLA_C_rho": g("CLDN4", "HLA-C", "unadj_rho"),
        "CLDN4_vs_HLA_C_p": g("CLDN4", "HLA-C", "unadj_p"),
        "CLDN4_vs_MHC1_HLAabc_rho": g("CLDN4", "MHC1_HLAabc", "unadj_rho"),
        "CLDN4_vs_MHC1_HLAabc_p": g("CLDN4", "MHC1_HLAabc", "unadj_p"),
        "CLDN4_vs_CD274_partial_rho": g("CLDN4", "CD274", "partial_rho"),
        "CLDN4_vs_CD274_partial_p": g("CLDN4", "CD274", "partial_p"),
        "CLDN4_vs_HLA_A_partial_rho": g("CLDN4", "HLA-A", "partial_rho"),
        "CLDN4_vs_HLA_A_partial_p": g("CLDN4", "HLA-A", "partial_p"),
        "CLDN4_vs_HLA_B_partial_rho": g("CLDN4", "HLA-B", "partial_rho"),
        "CLDN4_vs_HLA_B_partial_p": g("CLDN4", "HLA-B", "partial_p"),
        "CLDN4_vs_HLA_C_partial_rho": g("CLDN4", "HLA-C", "partial_rho"),
        "CLDN4_vs_HLA_C_partial_p": g("CLDN4", "HLA-C", "partial_p"),
        "CLDN4_vs_MHC1_HLAabc_partial_rho": g("CLDN4", "MHC1_HLAabc", "partial_rho"),
        "CLDN4_vs_MHC1_HLAabc_partial_p": g("CLDN4", "MHC1_HLAabc", "partial_p"),
        "partial_covariate": est_name,
        "ICI_labels": "not_deposited",
        "cd8_recut": "no",
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "question": "GSE31210 LUAD n=226: CLDN4 vs CD274 and HLA-A/B/C after ESTIMATE",
        "additive": True,
        "new_cut_only": True,
        "cd8_taken_as_given": True,
        "coverage": cov,
        "probe_used": probe_used,
        "probe_rule": probe_rule,
        "one_row": one.iloc[0].to_dict(),
        "stats": stats_df.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    write_finding(stats_df, cov, probe_used, probe_rule, n_tumor, n_matrix)

    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    def scatter(ax, x, y, xlab, ylab, rho, p, n, color):
        ax.scatter(x, y, s=12, alpha=0.45, c=color, linewidths=0)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_rho(rho)}  p={fmt_p(p)}  n={n}", fontsize=10)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(ax, logx["CLDN4"], logx["CD274"], "CLDN4 log2(MAS5+1)", "CD274 log2(MAS5+1)",
            g("CLDN4", "CD274", "unadj_rho"), g("CLDN4", "CD274", "unadj_p"),
            n_tumor, "#4C78A8")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd274.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    for ax, gene, color in zip(
        axes, HLA_GENES, ("#E45756", "#54A24B", "#F58518")
    ):
        scatter(ax, logx["CLDN4"], logx[gene], "CLDN4 log2(MAS5+1)", f"{gene} log2(MAS5+1)",
                g("CLDN4", gene, "unadj_rho"), g("CLDN4", gene, "unadj_p"),
                n_tumor, color)
    fig.suptitle(f"GSE31210 LUAD tumors n={n_tumor}", y=1.03)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_hla.png", bbox_inches="tight")
    plt.close(fig)

    forest = [
        ("CLDN4 vs CD274", g("CLDN4", "CD274", "unadj_rho"), "#4C78A8"),
        ("CLDN4 vs HLA-A", g("CLDN4", "HLA-A", "unadj_rho"), "#E45756"),
        ("CLDN4 vs HLA-B", g("CLDN4", "HLA-B", "unadj_rho"), "#54A24B"),
        ("CLDN4 vs HLA-C", g("CLDN4", "HLA-C", "unadj_rho"), "#F58518"),
        ("CLDN4 vs MHC-I mean-z", g("CLDN4", "MHC1_HLAabc", "unadj_rho"), "#B279A2"),
        ("CLDN4 vs CD274 | ESTIMATE", g("CLDN4", "CD274", "partial_rho"), "#9ecae1"),
        ("CLDN4 vs HLA-A | ESTIMATE", g("CLDN4", "HLA-A", "partial_rho"), "#f4a582"),
        ("CLDN4 vs HLA-B | ESTIMATE", g("CLDN4", "HLA-B", "partial_rho"), "#a1d99b"),
        ("CLDN4 vs HLA-C | ESTIMATE", g("CLDN4", "HLA-C", "partial_rho"), "#fdae6b"),
        ("CLDN4 vs MHC-I | ESTIMATE", g("CLDN4", "MHC1_HLAabc", "partial_rho"), "#d4b9da"),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    y = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=1)
    for i, (_lab, rho, c) in enumerate(forest):
        ax.plot(rho, i, "o", color=c, ms=7)
        ax.plot([0, rho], [i, i], color=c, lw=2)
    ax.set_yticks(y, [r[0] for r in forest], fontsize=9)
    ax.set_xlabel(f"Spearman ρ (n={n_tumor}; partial = rank residual on ESTIMATEScore)")
    ax.set_title("GSE31210 LUAD tumors — new cut (not CD8)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_spearman_forest.png", bbox_inches="tight")
    plt.close(fig)

    mat = pd.DataFrame({
        "CLDN4": logx["CLDN4"].to_numpy(float),
        "CD274": logx["CD274"].to_numpy(float),
        "HLA-A": logx["HLA-A"].to_numpy(float),
        "HLA-B": logx["HLA-B"].to_numpy(float),
        "HLA-C": logx["HLA-C"].to_numpy(float),
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
    ax.set_title(f"GSE31210 LUAD tumors n={n_tumor}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_correlation_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    print(one.T.to_string())
    print("wrote", TABLES)
    print("wrote", FIGURES)
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
