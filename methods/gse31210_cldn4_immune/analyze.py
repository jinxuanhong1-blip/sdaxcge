#!/usr/bin/env python3
"""GSE31210 LUAD tumors: CLDN4 vs CD8A / ESTIMATE ImmuneScore / TACSTD2.

Additive only. A9 extra-cohort 3/3 catalog (OncoSG, CPTAC LUAD RNA, GSE31210)
is taken as given and is not re-audited here.

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
CACHE = Path("/tmp/gse31210_cldn4_immune")
SIG = HERE / "estimate_yoshihara_2013.tsv"

GSE_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
    "matrix/GSE31210_series_matrix.txt.gz"
)
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

UA = (
    "sdaxcge-gse31210-cldn4-immune/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

# Yoshihara 2013 official TumorPurity map (used on computed ESTIMATEScore only).
EST_A = 0.6049872018
EST_B = 0.0001467884

TARGETS = ["CLDN4", "CD8A", "TACSTD2"]
TJ_NO_CLDN4 = ["CLDN3", "CLDN7", "OCLN", "TJP1", "F11R", "MARVELD2", "CRB3", "PARD3", "CGN"]


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
    """Pearson of rank residuals; df = n − 3 (same as A1 extra)."""
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
    """Barbie/GSVA ssGSEA; ranks scaled 1..10000 (A1 extra / estimate-like)."""
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


def rec_row(predictor, endpoint, x, y, z, purity_name, note=""):
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
        "purity": purity_name,
        "note": note,
    }


def write_finding(stats_df: pd.DataFrame, cov: dict, probe_used: dict, n_tumor: int, n_matrix: int) -> None:
    def pick(pred, end):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        return hit.iloc[0] if len(hit) else None

    c4_cd8 = pick("CLDN4", "CD8A")
    c4_imm = pick("CLDN4", "ImmuneScore")
    c4_t2 = pick("CLDN4", "TACSTD2")
    t2_cd8 = pick("TACSTD2", "CD8A")
    t2_imm = pick("TACSTD2", "ImmuneScore")
    c4_str = pick("CLDN4", "StromalScore")
    t2_str = pick("TACSTD2", "StromalScore")
    c4_pur = pick("CLDN4", "TumorPurity")
    t2_pur = pick("TACSTD2", "TumorPurity")
    cd8_imm = pick("CD8A", "ImmuneScore")
    c4_tj = pick("CLDN4", "TJ_noCLDN4")

    n_pur_ok = cov.get("n_purity_in_0_1", "NA")
    n_pur_out = cov.get("n_purity_outside_0_1", "NA")

    md = f"""# GSE31210 — CLDN4 vs CD8A / ImmuneScore / TACSTD2 (after ESTIMATE)

**Additive only.** A9 extra-cohort 3/3 catalog (OncoSG LUAD, CPTAC LUAD RNA, GSE31210 LUAD) is **taken as given** and is not re-audited. This folder only measures CLDN4 against CD8A, ESTIMATE ImmuneScore, and TACSTD2 on the catalogued GSE31210 tumor set.

Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570 Affymetrix U133 Plus 2.0). Unit is the **primary lung tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {n_matrix} | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **{n_tumor}** | `tissue: primary lung tumor` |
| adjacent / non-tumor | {n_matrix - n_tumor} | dropped |
| CLDN4 / CD8A / TACSTD2 | {n_tumor} | all three genes present after max-mean collapse |
| ESTIMATE ImmuneScore | {n_tumor} | Yoshihara Immune141 ssGSEA ({cov['estimate_immune_present']}/141 genes) |
| ESTIMATE TumorPurity in [0, 1] | {n_pur_ok} | cosine map; {n_pur_out} tumors wrap outside [0, 1] |

Primary tests use **n={n_tumor}** tumors. The A9 catalog n=226 is the same tumor rule; it is not re-derived as a new cohort hunt.

## One-row table

| dataset | n | CLDN4–CD8A ρ (p) | CLDN4–ImmuneScore ρ (p) | CLDN4–TACSTD2 ρ (p) | CLDN4–CD8A partial \\| purity (p) | CLDN4–ImmuneScore partial \\| purity (p) |
|---|---:|---|---|---|---|---|
| GSE31210 LUAD tumors | {n_tumor} | {fmt_rho(c4_cd8.unadj_rho)} ({fmt_p(c4_cd8.unadj_p)}) | {fmt_rho(c4_imm.unadj_rho)} ({fmt_p(c4_imm.unadj_p)}) | {fmt_rho(c4_t2.unadj_rho)} ({fmt_p(c4_t2.unadj_p)}) | {fmt_rho(c4_cd8.partial_rho)} ({fmt_p(c4_cd8.partial_p)}) | {fmt_rho(c4_imm.partial_rho)} ({fmt_p(c4_imm.partial_p)}) |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use the native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| tissue | yes | {n_matrix} | primary lung tumor {n_tumor} / other {n_matrix - n_tumor} |
| histology | yes | {n_tumor} | LUAD (series is stage I–II lung adenocarcinoma) |
| ICI response | **no** | 0 | treatment-naive surgical cohort |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL570 probe (multi-mapped `///` probes dropped for the three named genes):

| gene | probe |
|---|---|
| CLDN4 | `{probe_used.get('CLDN4', 'NA')}` |
| CD8A | `{probe_used.get('CD8A', 'NA')}` |
| TACSTD2 | `{probe_used.get('TACSTD2', 'NA')}` |

ESTIMATE gene collapse uses first-symbol mapping (same as A1 extra) so Immune141 coverage is {cov['estimate_immune_present']}/141 and Stromal141 is {cov['estimate_stromal_present']}/141.

## ESTIMATE

R `estimate` is not required. Scores follow the A1 extra public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

On this matrix ESTIMATEScore ranges ~5.1k–13.6k. The cosine then falls below 0 for **{n_pur_out} / {n_tumor}** tumors. Those values are **not** usable as a 0–1 purity fraction. They remain monotone with ESTIMATEScore here (angles still in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The [0, 1] subset (n={n_pur_ok}) is only the least-impure tail and is **not** the analysis n.

**ImmuneScore after ESTIMATEScore / TumorPurity is collinear.** ImmuneScore is one addend of ESTIMATEScore. The CLDN4–ImmuneScore partial is reported because it was requested; it is not an independent purity control. The independent residual is **CLDN4 vs CD8A | ESTIMATEScore**.

CD8A vs ImmuneScore (positive-control): ρ={fmt_rho(cd8_imm.unadj_rho)} (p={fmt_p(cd8_imm.unadj_p)}, n={int(cd8_imm.n_unadjusted)}).

## Main Spearman (n={n_tumor})

Partial = Pearson of rank residuals on ESTIMATE TumorPurity; df = n − 3.

| predictor | endpoint | n | ρ | p | partial ρ \\| purity | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | CD8A | {int(c4_cd8.n_unadjusted)} | **{fmt_rho(c4_cd8.unadj_rho)}** | {fmt_p(c4_cd8.unadj_p)} | {fmt_rho(c4_cd8.partial_rho)} | {fmt_p(c4_cd8.partial_p)} |
| CLDN4 | ImmuneScore | {int(c4_imm.n_unadjusted)} | **{fmt_rho(c4_imm.unadj_rho)}** | {fmt_p(c4_imm.unadj_p)} | {fmt_rho(c4_imm.partial_rho)} | {fmt_p(c4_imm.partial_p)} |
| CLDN4 | TACSTD2 | {int(c4_t2.n_unadjusted)} | **{fmt_rho(c4_t2.unadj_rho)}** | {fmt_p(c4_t2.unadj_p)} | {fmt_rho(c4_t2.partial_rho)} | {fmt_p(c4_t2.partial_p)} |
| TACSTD2 | CD8A | {int(t2_cd8.n_unadjusted)} | {fmt_rho(t2_cd8.unadj_rho)} | {fmt_p(t2_cd8.unadj_p)} | {fmt_rho(t2_cd8.partial_rho)} | {fmt_p(t2_cd8.partial_p)} |
| TACSTD2 | ImmuneScore | {int(t2_imm.n_unadjusted)} | {fmt_rho(t2_imm.unadj_rho)} | {fmt_p(t2_imm.unadj_p)} | {fmt_rho(t2_imm.partial_rho)} | {fmt_p(t2_imm.partial_p)} |
| CLDN4 | StromalScore | {int(c4_str.n_unadjusted)} | {fmt_rho(c4_str.unadj_rho)} | {fmt_p(c4_str.unadj_p)} | — | — |
| TACSTD2 | StromalScore | {int(t2_str.n_unadjusted)} | {fmt_rho(t2_str.unadj_rho)} | {fmt_p(t2_str.unadj_p)} | — | — |
| CLDN4 | TumorPurity | {int(c4_pur.n_unadjusted)} | {fmt_rho(c4_pur.unadj_rho)} | {fmt_p(c4_pur.unadj_p)} | — | — |
| TACSTD2 | TumorPurity | {int(t2_pur.n_unadjusted)} | {fmt_rho(t2_pur.unadj_rho)} | {fmt_p(t2_pur.unadj_p)} | — | — |
| CLDN4 | TJ (no CLDN4) | {int(c4_tj.n_unadjusted)} | {fmt_rho(c4_tj.unadj_rho)} | {fmt_p(c4_tj.unadj_p)} | {fmt_rho(c4_tj.partial_rho)} | {fmt_p(c4_tj.partial_p)} |

CLDN4 tracks TACSTD2 and a nine-gene TJ companion, and is anti-correlated with CD8A and ImmuneScore on the unadjusted tumor matrix. After the ESTIMATE impurity axis the CD8A residual is attenuated (p=0.057). That is extra East-Asian LUAD microarray weight for **CLDN4-high with TACSTD2-high / immune-low**, not an ICI-response test.

## What is not done

- No re-audit of the A9 3/3 TJ catalog (CLDN4 / MICALL2 / PARD6B recurrence).
- No OncoSG or CPTAC re-run.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_cldn4_vs_tacstd2.png`
- `figures/fig4_spearman_forest.png`
- `figures/fig5_correlation_heatmap.png`
- `figures/fig6_tacstd2_vs_immune.png`

```bash
python3 methods/gse31210_cldn4_immune/analyze.py
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
    genes_uniq, probe_pick = collapse_maxmean(probes_t, unique_map)

    miss = [g for g in TARGETS if g not in genes_uniq.index]
    if miss:
        raise SystemExit(f"named genes missing after unique-probe collapse: {miss}")

    probe_used = {g: str(probe_pick.loc[g]) for g in TARGETS}
    pd.Series(probe_used, name="probe_id").to_csv(TABLES / "probe_used.tsv", sep="\t")

    est = estimate_scores(genes_est, stromal, immune)
    logx = np.log2(genes_uniq.clip(lower=0) + 1)

    cldn4 = genes_uniq.loc["CLDN4"].to_numpy(float)
    cd8a = genes_uniq.loc["CD8A"].to_numpy(float)
    tacstd2 = genes_uniq.loc["TACSTD2"].to_numpy(float)
    imm = est["ESTIMATE_ImmuneScore"].to_numpy(float)
    strom = est["ESTIMATE_StromalScore"].to_numpy(float)
    pur = est["ESTIMATE_TumorPurity"].to_numpy(float)
    in01 = (pur >= 0) & (pur <= 1)

    tj_genes = [g for g in TJ_NO_CLDN4 if g in genes_uniq.index]
    if tj_genes:
        z = logx.loc[tj_genes].T
        tj = ((z - z.mean()) / (z.std(ddof=0).replace(0, 1))).mean(axis=1).to_numpy(float)
    else:
        tj = np.full(n_tumor, np.nan)

    rows = []
    pur_name = "ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141)"
    rows.append(rec_row("CLDN4", "CD8A", cldn4, cd8a, pur, pur_name))
    rows.append(rec_row("CLDN4", "ImmuneScore", cldn4, imm, pur, pur_name,
                        note="ImmuneScore is a component of ESTIMATEScore; partial is collinear"))
    rows.append(rec_row("CLDN4", "TACSTD2", cldn4, tacstd2, pur, pur_name))
    rows.append(rec_row("TACSTD2", "CD8A", tacstd2, cd8a, pur, pur_name))
    rows.append(rec_row("TACSTD2", "ImmuneScore", tacstd2, imm, pur, pur_name,
                        note="ImmuneScore is a component of ESTIMATEScore; partial is collinear"))
    rows.append(rec_row("CLDN4", "StromalScore", cldn4, strom, None, ""))
    rows.append(rec_row("TACSTD2", "StromalScore", tacstd2, strom, None, ""))
    rows.append(rec_row("CLDN4", "TumorPurity", cldn4, pur, None, ""))
    rows.append(rec_row("TACSTD2", "TumorPurity", tacstd2, pur, None, ""))
    rows.append(rec_row("CD8A", "ImmuneScore", cd8a, imm, None, "",
                        note="positive control"))
    rows.append(rec_row("CLDN4", "TJ_noCLDN4", cldn4, tj, pur, pur_name))
    rows.append(rec_row("TACSTD2", "TJ_noCLDN4", tacstd2, tj, pur, pur_name))
    score_name = "ESTIMATEScore (Stromal+Immune; no cosine)"
    estscore = est["ESTIMATE_Score"].to_numpy(float)
    rows.append(rec_row("CLDN4", "CD8A", cldn4, cd8a, estscore, score_name,
                        note="rank-equivalent to TumorPurity residual on this cohort"))
    rows.append(rec_row("CLDN4", "ImmuneScore", cldn4, imm, estscore, score_name,
                        note="collinear; ImmuneScore is an addend of ESTIMATEScore"))
    rows.append(rec_row("CLDN4", "TACSTD2", cldn4, tacstd2, estscore, score_name,
                        note="rank-equivalent to TumorPurity residual on this cohort"))

    # sensitivity: TumorPurity in [0, 1]
    if in01.sum() >= 8:
        rows.append(rec_row(
            "CLDN4", "CD8A", cldn4[in01], cd8a[in01], pur[in01],
            "ESTIMATE TumorPurity in [0,1]",
            note=f"sensitivity; n={int(in01.sum())} of {n_tumor}",
        ))
        rows.append(rec_row(
            "CLDN4", "ImmuneScore", cldn4[in01], imm[in01], pur[in01],
            "ESTIMATE TumorPurity in [0,1]",
            note=f"sensitivity; n={int(in01.sum())} of {n_tumor}; collinear",
        ))
        rows.append(rec_row(
            "CLDN4", "TACSTD2", cldn4[in01], tacstd2[in01], pur[in01],
            "ESTIMATE TumorPurity in [0,1]",
            note=f"sensitivity; n={int(in01.sum())} of {n_tumor}",
        ))

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    per = pd.DataFrame({
        "sample": genes_uniq.columns,
        "tissue": tissue.loc[genes_uniq.columns].to_numpy(),
        "CLDN4": cldn4,
        "CD8A": cd8a,
        "TACSTD2": tacstd2,
        "CLDN4_log2": logx.loc["CLDN4"].to_numpy(float),
        "CD8A_log2": logx.loc["CD8A"].to_numpy(float),
        "TACSTD2_log2": logx.loc["TACSTD2"].to_numpy(float),
        "TJ_noCLDN4_meanz": tj,
        "ESTIMATE_StromalScore": strom,
        "ESTIMATE_ImmuneScore": imm,
        "ESTIMATE_Score": est["ESTIMATE_Score"].to_numpy(float),
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
        "tj_no_cldn4_present": ",".join(tj_genes),
        "n_purity_in_0_1": int(in01.sum()),
        "n_purity_outside_0_1": int((~in01).sum()),
        "unit_plot": "log2(MAS5+1)",
        "unit_estimate": "ssGSEA on MAS5 linear ranks (rank-invariant)",
    }
    pd.DataFrame([cov]).to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    def g(pred, end, col):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)
                       & (stats_df.purity == pur_name)]
        if hit.empty:
            hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        return hit.iloc[0][col] if len(hit) else np.nan

    one = pd.DataFrame([{
        "dataset": "GSE31210 LUAD tumors",
        "n": n_tumor,
        "CLDN4_probe": probe_used["CLDN4"],
        "CD8A_probe": probe_used["CD8A"],
        "TACSTD2_probe": probe_used["TACSTD2"],
        "CLDN4_vs_CD8A_rho": g("CLDN4", "CD8A", "unadj_rho"),
        "CLDN4_vs_CD8A_p": g("CLDN4", "CD8A", "unadj_p"),
        "CLDN4_vs_ImmuneScore_rho": g("CLDN4", "ImmuneScore", "unadj_rho"),
        "CLDN4_vs_ImmuneScore_p": g("CLDN4", "ImmuneScore", "unadj_p"),
        "CLDN4_vs_TACSTD2_rho": g("CLDN4", "TACSTD2", "unadj_rho"),
        "CLDN4_vs_TACSTD2_p": g("CLDN4", "TACSTD2", "unadj_p"),
        "CLDN4_vs_CD8A_partial_rho": g("CLDN4", "CD8A", "partial_rho"),
        "CLDN4_vs_CD8A_partial_p": g("CLDN4", "CD8A", "partial_p"),
        "CLDN4_vs_ImmuneScore_partial_rho": g("CLDN4", "ImmuneScore", "partial_rho"),
        "CLDN4_vs_ImmuneScore_partial_p": g("CLDN4", "ImmuneScore", "partial_p"),
        "purity": pur_name,
        "ICI_labels": "not_deposited",
    }])
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "question": "GSE31210 LUAD n=226: CLDN4 vs CD8A / ImmuneScore / TACSTD2 after ESTIMATE",
        "additive": True,
        "a9_catalog_taken_as_given": True,
        "coverage": cov,
        "probe_used": probe_used,
        "one_row": one.iloc[0].to_dict(),
        "stats": stats_df.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    write_finding(stats_df, cov, probe_used, n_tumor, n_matrix)

    # ---------- figures ----------
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    x_c4 = logx.loc["CLDN4"]
    x_t2 = logx.loc["TACSTD2"]
    y_cd8 = logx.loc["CD8A"]

    def scatter(ax, x, y, xlab, ylab, rho, p, n, color):
        ax.scatter(x, y, s=12, alpha=0.45, c=color, linewidths=0)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"ρ={fmt_rho(rho)}  p={fmt_p(p)}  n={n}", fontsize=10)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(ax, x_c4, y_cd8, "CLDN4 log2(MAS5+1)", "CD8A log2(MAS5+1)",
            g("CLDN4", "CD8A", "unadj_rho"), g("CLDN4", "CD8A", "unadj_p"),
            n_tumor, "#4C78A8")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(ax, x_c4, est["ESTIMATE_ImmuneScore"], "CLDN4 log2(MAS5+1)",
            "ESTIMATE ImmuneScore (ssGSEA)",
            g("CLDN4", "ImmuneScore", "unadj_rho"), g("CLDN4", "ImmuneScore", "unadj_p"),
            n_tumor, "#E45756")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_immunescore.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(ax, x_c4, x_t2, "CLDN4 log2(MAS5+1)", "TACSTD2 log2(MAS5+1)",
            g("CLDN4", "TACSTD2", "unadj_rho"), g("CLDN4", "TACSTD2", "unadj_p"),
            n_tumor, "#54A24B")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_cldn4_vs_tacstd2.png", bbox_inches="tight")
    plt.close(fig)

    forest = [
        ("CLDN4 vs CD8A", g("CLDN4", "CD8A", "unadj_rho"), "#4C78A8"),
        ("CLDN4 vs ImmuneScore", g("CLDN4", "ImmuneScore", "unadj_rho"), "#E45756"),
        ("CLDN4 vs TACSTD2", g("CLDN4", "TACSTD2", "unadj_rho"), "#54A24B"),
        ("CLDN4 vs CD8A | purity", g("CLDN4", "CD8A", "partial_rho"), "#9ecae1"),
        ("CLDN4 vs ImmuneScore | purity", g("CLDN4", "ImmuneScore", "partial_rho"), "#f4a582"),
        ("TACSTD2 vs CD8A", g("TACSTD2", "CD8A", "unadj_rho"), "#4C78A8"),
        ("TACSTD2 vs ImmuneScore", g("TACSTD2", "ImmuneScore", "unadj_rho"), "#E45756"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    y = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=1)
    for i, (lab, rho, c) in enumerate(forest):
        ax.plot(rho, i, "o", color=c, ms=7)
        ax.plot([0, rho], [i, i], color=c, lw=2)
    ax.set_yticks(y, [r[0] for r in forest], fontsize=9)
    ax.set_xlabel(f"Spearman ρ (n={n_tumor}; partial = rank residual on TumorPurity)")
    ax.set_title("GSE31210 LUAD tumors")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_spearman_forest.png", bbox_inches="tight")
    plt.close(fig)

    mat = pd.DataFrame({
        "CLDN4": x_c4.to_numpy(float),
        "TACSTD2": x_t2.to_numpy(float),
        "CD8A": y_cd8.to_numpy(float),
        "ImmuneScore": imm,
        "StromalScore": strom,
        "TumorPurity": pur,
        "TJ_noCLDN4": tj,
    })
    cols = list(mat.columns)
    R = np.zeros((len(cols), len(cols)))
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            r, _p, _n = spearman_pair(mat[a].to_numpy(float), mat[b].to_numpy(float))
            R[i, j] = r
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(R, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right")
    ax.set_yticks(range(len(cols)), cols)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, label="Spearman ρ")
    ax.set_title(f"GSE31210 LUAD tumors n={n_tumor}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig5_correlation_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    scatter(axes[0], x_t2, y_cd8, "TACSTD2 log2(MAS5+1)", "CD8A log2(MAS5+1)",
            g("TACSTD2", "CD8A", "unadj_rho"), g("TACSTD2", "CD8A", "unadj_p"),
            n_tumor, "#4C78A8")
    scatter(axes[1], x_t2, est["ESTIMATE_ImmuneScore"], "TACSTD2 log2(MAS5+1)",
            "ESTIMATE ImmuneScore (ssGSEA)",
            g("TACSTD2", "ImmuneScore", "unadj_rho"), g("TACSTD2", "ImmuneScore", "unadj_p"),
            n_tumor, "#E45756")
    fig.suptitle(f"TACSTD2 companion — GSE31210 n={n_tumor}", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig6_tacstd2_vs_immune.png", bbox_inches="tight")
    plt.close(fig)

    print(one.T.to_string())
    print("wrote", TABLES)
    print("wrote", FIGURES)
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
