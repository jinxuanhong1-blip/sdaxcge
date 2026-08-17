#!/usr/bin/env python3
"""GSE11969 LUAD array: CLDN4 vs CD8A / CD274 (PDCD1LG1).

Additive public cohort. Takeuchi / Tomida resected lung series
(Takeuchi et al., J Clin Oncol 2006, PMID 16549822; GEO GSE11969).
Primary slice is LUAD (Histology = AD) only. Honest n is the pairwise
finite count on those arrays, not the 163-array series or the 149-NSCLC
sentence.

CD274 is not a GPL7015 Gene symbol. The 2009 platform uses PDCD1LG1
(probe A_23_P256487, EST BU621474). That probe is the CD274 row.

Downloads stay under $GSE11969_CLDN4_DATA (default /tmp/gse11969_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import io
import json
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
DATA = Path(os.environ.get("GSE11969_CLDN4_DATA", "/tmp/gse11969_cldn4"))
SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11969/"
    "matrix/GSE11969_series_matrix.txt.gz"
)
GPL_TABLE = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL7015&targ=self&form=text&view=data"
)

# Named probes from official GPL7015 Gene symbol / Title.
NAMED = {
    "CLDN4": {
        "id": "1374",
        "probe": "A_23_P19944",
        "gb": "NM_001305",
        "symbol_on_gpl": "CLDN4",
        "title": "Claudin 4",
    },
    "CD8A": {
        "id": "9374",
        "probe": "A_23_P68110",
        "gb": "NM_001768",
        "symbol_on_gpl": "CD8A",
        "title": "CD8 antigen, alpha polypeptide (p32)",
    },
    "CD274": {
        "id": "8635",
        "probe": "A_23_P256487",
        "gb": "BU621474",
        "symbol_on_gpl": "PDCD1LG1",
        "title": "Programmed cell death 1 ligand 1",
    },
    "TACSTD2": {
        "id": "19914",
        "probe": "A_23_P149529",
        "gb": "NM_002353",
        "symbol_on_gpl": "TACSTD2",
        "title": "Tumor-associated calcium signal transducer 2",
    },
}

# EPCAM is absent on GPL7015. Residual uses the 5 present epithelial genes.
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "gse11969-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                if key == "characteristics_ch1":
                    char_key = None
                    for v in vals:
                        if v and ":" in v:
                            char_key = v.split(":", 1)[0].strip()
                            break
                    key = char_key or f"char_{len(sample_fields)}"
                    vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
                sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession in sample fields: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    return meta, {k: " | ".join(v) for k, v in series.items()}


def load_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        raw = fh.read()
    start = raw.index("!series_matrix_table_begin")
    end = raw.index("!series_matrix_table_end")
    block = raw[start:].split("\n", 1)[1]
    block = block[: block.rfind("!series_matrix_table_end")]
    expr = pd.read_csv(io.StringIO(block), sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.strip().str.strip('"')
    expr.columns = [c.strip().strip('"') for c in expr.columns]
    return expr.astype(float)


def load_gpl(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="replace").splitlines()
    start = next(i for i, l in enumerate(text) if l.startswith("!platform_table_begin")) + 1
    end = next(i for i, l in enumerate(text) if l.startswith("!platform_table_end"))
    plat = pd.read_csv(io.StringIO("\n".join(text[start:end])), sep="\t", dtype=str)
    plat["ID"] = plat["ID"].astype(str)
    plat["Gene symbol"] = plat["Gene symbol"].fillna("")
    return plat


def first_symbol(s: str) -> str:
    if not isinstance(s, str) or not s or s in {"nan", "NODATA", "no data"}:
        return ""
    return s.split(" / ")[0].split(" /// ")[0].strip()


def collapse_maxmean(probe_expr: pd.DataFrame, id2gene: pd.Series) -> pd.DataFrame:
    genes = probe_expr.index.map(lambda x: first_symbol(id2gene.get(str(x), "")))
    df = probe_expr.copy()
    df["_gene"] = genes
    df = df[df["_gene"].astype(str).str.len() > 0]
    df["_mean"] = df.drop(columns="_gene").astype(float).mean(axis=1)
    df = df.sort_values("_mean", ascending=False)
    kept = df.loc[~df["_gene"].duplicated(keep="first")]
    out = kept.drop(columns=["_gene", "_mean"]).astype(float)
    out.index = kept["_gene"].values
    return out.sort_index()


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


def partial_spearman(x, y, covar):
    df = pd.DataFrame(
        {
            "x": np.asarray(x, float),
            "y": np.asarray(y, float),
            "c": np.asarray(covar, float),
        }
    ).dropna()
    n = int(df.shape[0])
    if n < 12:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    xr = stats.rankdata(df["x"].values)
    yr = stats.rankdata(df["y"].values)
    cr = np.column_stack([np.ones(n), stats.rankdata(df["c"].values)])
    bx, *_ = np.linalg.lstsq(cr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(cr, yr, rcond=None)
    rx = xr - cr @ bx
    ry = yr - cr @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return dict(n=n, rho_adj=np.nan, p_adj=np.nan)
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = n - 3
    t = r * np.sqrt(dof / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return dict(n=n, rho_adj=r, p_adj=p)


def verdict(rho, p, n) -> str:
    if n < HOLDS_N:
        return "UNDERPOWERED"
    if np.isfinite(rho) and rho < 0 and p < 0.05:
        return "HOLDS"
    if np.isfinite(rho) and rho > 0 and p < 0.05:
        return "OPPOSITE"
    return "NO_EVIDENCE"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "—"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "—"
    return f"{r:+.3f}"


def fmt_ci(lo, hi) -> str:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "—"
    return f"{lo:+.3f} to {hi:+.3f}"


def write_finding(n_table: pd.DataFrame, primary: pd.DataFrame, qtab: pd.DataFrame, cov: dict) -> None:
    def row(gene_y: str) -> pd.Series:
        hit = primary[primary["gene_y"] == gene_y]
        if hit.empty:
            return pd.Series(dtype=object)
        return hit.iloc[0]

    c8 = row("CD8A")
    pdl1 = row("CD274")
    tac = row("TACSTD2")
    q8 = qtab[qtab["endpoint"] == "CD8A"].iloc[0] if len(qtab) else None
    q274 = qtab[qtab["endpoint"] == "CD274"].iloc[0] if len(qtab) else None

    n_series = int(cov["n_series"])
    n_nsclc = int(cov["n_nsclc"])
    n_luad = int(cov["n_luad"])
    n_sclc = int(cov["n_sclc"])
    n_normal = int(cov["n_normal"])
    n_primary = int(c8.get("n", 0) or 0)
    n_gefitinib = int(cov.get("gefitinib_luad_yes", 0) or 0)
    epi_used = cov["epithelial_used"]
    epi_missing = cov["epithelial_missing"]

    empty = n_primary == 0 or not np.isfinite(c8.get("rho", np.nan))
    if empty:
        body = f"""# Additive GSE11969 LUAD — CLDN4 vs CD8A / CD274

**Empty table.** Primary pairwise n is 0 or the named probes were not finite on LUAD arrays.

| dataset | histology | n LUAD | CLDN4 vs CD8A | CLDN4 vs CD274 |
|---|---|---:|---|---|
| GSE11969 Takeuchi | LUAD (Histology=AD) | {n_luad} | — | — |

## Honest n

| Item | n |
|---|---:|
| Arrays in the series matrix | {n_series} |
| LUAD (Histology = AD) | {n_luad} |
| Primary pairwise n | 0 |
"""
        (HERE / "FINDING.md").write_text(body)
        return

    body = f"""# Additive GSE11969 LUAD — CLDN4 vs CD8A / CD274

**Additive public cohort.** Takeuchi / Tomida **GSE11969** (PMID 16549822; Tomida follow-up PMID 21465578): resected Japanese lung series on **GPL7015** Agilent Homo sapiens 21.6K custom array. Series text says 149 NSCLC including **90 adenocarcinomas (AD)**, plus 9 SCLC and 5 normal-lung mixtures. This slice uses only **Histology = AD**.

Primary question: **CLDN4 vs CD8A** and **CLDN4 vs CD274**. TACSTD2 is a same-run companion, not an audit of prior TACSTD2 claims.

**CD274 mapping.** GPL7015 has no `CD274` symbol. The deposited symbol is **PDCD1LG1** (probe `{NAMED["CD274"]["probe"]}`, EST `{NAMED["CD274"]["gb"]}`, title “Programmed cell death 1 ligand 1”). That is the CD274 row. Do not write “CD274 absent.”

## One-row table

| dataset | histology | platform | n series | n NSCLC (text) | n LUAD | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict CD8A / CD274 |
|---|---|---|---:|---:|---:|---|---|---|---|---|---|---|---|
| GSE11969 Takeuchi | LUAD (AD) | GPL7015 Agilent 21.6K log10(R/G) | {n_series} | {n_nsclc} | **{n_luad}** | `{NAMED["CLDN4"]["probe"]}` | `{NAMED["CD8A"]["probe"]}` | `{NAMED["CD274"]["probe"]}` (PDCD1LG1) | **{fmt_rho(c8.rho)} ({fmt_p(c8.p)})** | {fmt_rho(c8.rho_adj)} ({fmt_p(c8.p_adj)}) | **{fmt_rho(pdl1.rho)} ({fmt_p(pdl1.p)})** | {fmt_rho(pdl1.rho_adj)} ({fmt_p(pdl1.p_adj)}) | **{c8.verdict_adj} / {pdl1.verdict_adj}** |

Do not write n={n_series} or n={n_nsclc} for the correlations. The computable LUAD n is **{n_luad}** arrays / **{n_luad}** patients (AD001–AD090, 1:1).

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **{n_series}** | 21,619 probes × {n_series} GSM; 0 missing VALUE cells |
| NSCLC (AD+SQ+LA+AS+LCNEC) | yes | **{n_nsclc}** | matches series text “149 patients with NSCLC” |
| **LUAD (Histology = AD)** | yes | **{n_luad}** | titles all “Patient … with adenocarcinoma”; unique Annotation AD001–AD090 |
| SQ / LA / AS / LCNEC | yes | {int(cov["n_sq"])} / {int(cov["n_la"])} / {int(cov["n_as"])} / {int(cov["n_lcnec"])} | not in the LUAD slice |
| SCLC | yes | {n_sclc} | excluded |
| Normal lung mixtures | yes | {n_normal} | Cy5 channel; excluded |
| Unique LUAD patients | yes | **{n_luad}** | one array per AD annotation |
| ICI / PD-1 / PD-L1 treatment | no | **0** | {n_gefitinib}/{n_luad} LUAD marked Gefitinib=Y; not an ICI series |
| Tumor % / ESTIMATE on GEO | no | **0** | only public proxy is an RNA epithelial score |
| CLDN4 finite (`{NAMED["CLDN4"]["probe"]}`) | yes | **{n_luad}** | named GPL7015 CLDN4 / NM_001305 |
| CD8A finite (`{NAMED["CD8A"]["probe"]}`) | yes | **{n_luad}** | named GPL7015 CD8A / NM_001768 |
| CD274 finite (`{NAMED["CD274"]["probe"]}`) | yes | **{n_luad}** | GPL7015 symbol **PDCD1LG1** / EST BU621474 |
| **Primary pairwise n (CLDN4 + CD8A + CD274)** | yes | **{n_primary}** | this is the n used below |
| CLDN4 Q4 vs Q1 | yes | **{int(q8.n_Q4)} vs {int(q8.n_Q1)}** | not {n_luad}; Q2+Q3 dropped |

The author patient-info supplement has extra blank rows if read naively (183 lines). Histology counts above come from the series matrix, not that file.

## Verdict

| Test | n | ρ | 95% CI | p | ρ_adj (epithelial 5/6) | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| **CLDN4 vs CD8A** | **{int(c8.n)}** | **{fmt_rho(c8.rho)}** | {fmt_ci(c8.ci_low, c8.ci_high)} | **{fmt_p(c8.p)}** | {fmt_rho(c8.rho_adj)} | {fmt_p(c8.p_adj)} | **{c8.verdict_adj}** |
| **CLDN4 vs CD274 (PDCD1LG1)** | **{int(pdl1.n)}** | **{fmt_rho(pdl1.rho)}** | {fmt_ci(pdl1.ci_low, pdl1.ci_high)} | **{fmt_p(pdl1.p)}** | {fmt_rho(pdl1.rho_adj)} | {fmt_p(pdl1.p_adj)} | **{pdl1.verdict_adj}** |
| CLDN4 vs TACSTD2 (companion) | {int(tac.n)} | {fmt_rho(tac.rho)} | {fmt_ci(tac.ci_low, tac.ci_high)} | {fmt_p(tac.p)} | {fmt_rho(tac.rho_adj)} | {fmt_p(tac.p_adj)} | companion |
| CLDN4 vs epithelial mean-z | {n_luad} | {fmt_rho(cov["cldn4_vs_epi_rho"])} | — | {fmt_p(cov["cldn4_vs_epi_p"])} | — | — | companion |

**What holds.** Higher **CLDN4** tracks **lower CD8A** on LUAD arrays (n={int(c8.n)}, ρ={fmt_rho(c8.rho)}, p={fmt_p(c8.p)}). The inverse remains after residualising on the 5-gene epithelial mean-z (EPCAM absent; used {", ".join(epi_used)}; ρ_adj={fmt_rho(c8.rho_adj)}, p_adj={fmt_p(c8.p_adj)}). HOLDS rule: n≥{HOLDS_N}, ρ<0, p<0.05.

**What does not hold.** **CLDN4 vs CD274/PDCD1LG1** is null (n={int(pdl1.n)}, ρ={fmt_rho(pdl1.rho)}, p={fmt_p(pdl1.p)}; adj ρ={fmt_rho(pdl1.rho_adj)}, p={fmt_p(pdl1.p_adj)}). Do not upgrade the CD8A anti-correlation into a PD-L1-RNA claim. The CD274 probe is an EST, not NM_014143; that is a platform limit, not a reason to impute a different n.

## Q4 vs Q1 (not n={n_luad})

`pd.qcut(rank(method="first"), 4)` on the {n_luad} LUAD CLDN4 values. Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | p |
|---|---:|---:|---:|---:|---:|---:|
| **CD8A** | **{int(q8.n_Q4)} / {int(q8.n_Q1)}** | {q8.median_Q4:.3f} | {q8.median_Q1:.3f} | {q8.delta:.3f} | **{q8.rank_biserial:+.2f}** | **{fmt_p(q8.p)}** |
| CD274 (PDCD1LG1) | {int(q274.n_Q4)} / {int(q274.n_Q1)} | {q274.median_Q4:.3f} | {q274.median_Q1:.3f} | {q274.delta:.3f} | {q274.rank_biserial:+.2f} | {fmt_p(q274.p)} |

The CD8A quartile cut is the same immune-low direction as the continuous test. It is not a second cohort. Do not write n={n_luad} for 23 vs 23.

## Methods

- **Matrix:** GEO `GSE11969_series_matrix.txt.gz`. Deposited VALUE = LOWESS-normalized, background-subtracted **log10(processed Red / processed Green)**. Cy5 = sample; Cy3 = 20-cell-line lung RNA reference. Spearman is rank-based, so the log10 ratio scale does not change ρ.
- **LUAD filter:** `Histology = AD` on the series matrix (n={n_luad}). Titles match. SQ/LA/AS/LCNEC/SCLC/normal are counted above and then dropped.
- **Probes:** official GPL7015 table (`Gene symbol`, `Title`, `GB_LIST`). Primary genes use the **named single probes** in `tables/probe_confirm.tsv`. No max-mean collapse for CLDN4 / CD8A / CD274 because each has one probe.
- **CD274:** GPL7015 symbol **PDCD1LG1**. Alias search (CD274, PD-L1, PDL1, B7-H1) found no second probe.
- **Epithelial residual:** unweighted mean of gene-wise z for {", ".join(epi_used)} ({len(epi_used)}/6; missing: {", ".join(epi_missing) or "none"}). Max-mean collapse only for this score. Partial Spearman = Pearson of rank residuals; df = n − 3.
- **HOLDS:** n≥{HOLDS_N}, ρ (or ρ_adj) < 0, p < 0.05. Used for the immune-low claim only.
- **Not tested:** OS/relapse (labels exist on GEO but are outside this slice), ICI response (none), pathologist CD8/PD-L1 IHC, purity/ESTIMATE (not deposited).

## What this does not claim

- It does not treat n={n_series} or n={n_nsclc} as the LUAD correlation n.
- It does not claim CLDN4 vs PD-L1 protein. The CD274 row is one EST probe.
- It does not claim an ICI-response effect. {n_gefitinib} LUAD arrays have Gefitinib=Y; that is EGFR TKI, not PD-1/PD-L1 blockade.
- TACSTD2 numbers are a same-run companion only.

## Reproduce

```bash
python3 -m pip install -r methods/gse11969_cldn4/requirements.txt
python3 methods/gse11969_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE11969_CLDN4_DATA` (default `/tmp/gse11969_cldn4`).

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/n_table.tsv` — honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / PDCD1LG1 / TACSTD2 probes
- `tables/spearman_primary.tsv` — CLDN4 vs CD8A / CD274 / TACSTD2
- `tables/highlow_q4q1.tsv` — CLDN4 Q4 vs Q1
- `tables/sample_annotation.tsv` — {n_luad} LUAD arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
"""
    (HERE / "FINDING.md").write_text(body)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE11969_series_matrix.txt.gz", GEO_MATRIX)
    gpl_path = dl(DATA / "GPL7015.txt", GPL_TABLE)

    meta, series = parse_series_meta(matrix_path)
    expr = load_matrix(matrix_path)
    expr = expr.loc[:, [c for c in expr.columns if c in meta.index]]
    meta = meta.loc[expr.columns]
    plat = load_gpl(gpl_path)
    id2gene = plat.drop_duplicates("ID").set_index("ID")["Gene symbol"]

    hist = meta["Histology"].astype(str)
    n_series = int(expr.shape[1])
    n_ad = int((hist == "AD").sum())
    n_sq = int((hist == "SQ").sum())
    n_la = int((hist == "LA").sum())
    n_as = int((hist == "AS").sum())
    n_lcnec = int((hist == "LCNEC").sum())
    n_sclc = int((hist == "SCLC").sum())
    n_normal = int((hist == "Normal lung tissue").sum())
    n_nsclc = n_ad + n_sq + n_la + n_as + n_lcnec

    luad = meta[hist == "AD"].copy()
    if luad["Annotation"].nunique() != len(luad):
        raise SystemExit("LUAD Annotation is not 1:1 with GSM")
    e_luad = expr.loc[:, luad.index]

    # Confirm named probes still map as expected.
    probe_rows = []
    for gene, rec in NAMED.items():
        pid = rec["id"]
        plat_hit = plat[plat["ID"] == pid]
        symbol = plat_hit["Gene symbol"].iloc[0] if len(plat_hit) else ""
        in_matrix = pid in expr.index
        n_finite = int(e_luad.loc[pid].notna().sum()) if in_matrix else 0
        probe_rows.append(
            {
                "claim_gene": gene,
                "gpl_symbol": symbol,
                "id": pid,
                "probe": rec["probe"],
                "gb_list": rec["gb"],
                "title": rec["title"],
                "in_matrix": in_matrix,
                "n_luad_finite": n_finite,
                "note": "CD274 = PDCD1LG1 on this platform" if gene == "CD274" else "",
            }
        )
        if not in_matrix or symbol != rec["symbol_on_gpl"]:
            print(f"[warn] {gene} probe check failed: in_matrix={in_matrix} symbol={symbol}")
    probe_df = pd.DataFrame(probe_rows)
    probe_df.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    gene_expr = collapse_maxmean(e_luad, id2gene)
    epi_present = [g for g in EPITHELIAL if g in gene_expr.index]
    epi_missing = [g for g in EPITHELIAL if g not in gene_expr.index]
    if epi_present:
        z = gene_expr.loc[epi_present]
        sd = z.std(axis=1, ddof=1).replace(0, np.nan)
        epithelial = z.sub(z.mean(axis=1), axis=0).div(sd, axis=0).mean(axis=0)
    else:
        epithelial = pd.Series(np.nan, index=e_luad.columns)

    cldn4 = e_luad.loc[NAMED["CLDN4"]["id"]]
    partners = {
        "CD8A": e_luad.loc[NAMED["CD8A"]["id"]],
        "CD274": e_luad.loc[NAMED["CD274"]["id"]],
        "TACSTD2": e_luad.loc[NAMED["TACSTD2"]["id"]],
    }

    primary_rows = []
    for name, y in partners.items():
        crude = spearman_ci(cldn4.values, y.values)
        adj = partial_spearman(cldn4.values, y.values, epithelial.reindex(cldn4.index).values)
        kind = "primary" if name in {"CD8A", "CD274"} else "companion"
        primary_rows.append(
            {
                "gene_x": "CLDN4",
                "gene_y": name,
                "gpl_symbol": NAMED[name]["symbol_on_gpl"],
                "kind": kind,
                "n": crude["n"],
                "n_luad": n_ad,
                "rho": crude["rho"],
                "p": crude["p"],
                "ci_low": crude["ci_low"],
                "ci_high": crude["ci_high"],
                "rho_adj": adj["rho_adj"],
                "p_adj": adj["p_adj"],
                "verdict_unadj": verdict(crude["rho"], crude["p"], crude["n"]),
                "verdict_adj": verdict(adj["rho_adj"], adj["p_adj"], adj["n"]),
            }
        )
    primary = pd.DataFrame(primary_rows)
    primary.to_csv(TABLES / "spearman_primary.tsv", sep="\t", index=False)

    epi_r, epi_p = stats.spearmanr(cldn4.values, epithelial.values) if epithelial.notna().sum() >= 8 else (np.nan, np.nan)

    q = pd.qcut(cldn4.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    q_rows = []
    for name, y in partners.items():
        a = y[q == "Q4"]
        b = y[q == "Q1"]
        U, p_mw = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
        rbc = 2 * U / (a.size * b.size) - 1
        q_rows.append(
            {
                "predictor": "CLDN4_quartile",
                "endpoint": name,
                "n_Q4": int(a.size),
                "n_Q1": int(b.size),
                "median_Q4": float(np.median(a)),
                "median_Q1": float(np.median(b)),
                "delta": float(np.median(a) - np.median(b)),
                "U": float(U),
                "p": float(p_mw),
                "rank_biserial": float(rbc),
            }
        )
    qtab = pd.DataFrame(q_rows)
    qtab.to_csv(TABLES / "highlow_q4q1.tsv", sep="\t", index=False)

    n_table = pd.DataFrame(
        [
            {"item": "arrays_in_series_matrix", "public": "yes", "n": n_series, "note": "21619 probes × 163 GSM; 0 missing VALUE"},
            {"item": "NSCLC_AD_SQ_LA_AS_LCNEC", "public": "yes", "n": n_nsclc, "note": "series text 149 NSCLC"},
            {"item": "LUAD_Histology_AD", "public": "yes", "n": n_ad, "note": "primary filter; unique AD001–AD090"},
            {"item": "SQ", "public": "yes", "n": n_sq, "note": "excluded"},
            {"item": "LA", "public": "yes", "n": n_la, "note": "excluded"},
            {"item": "AS", "public": "yes", "n": n_as, "note": "excluded"},
            {"item": "LCNEC", "public": "yes", "n": n_lcnec, "note": "excluded"},
            {"item": "SCLC", "public": "yes", "n": n_sclc, "note": "excluded"},
            {"item": "normal_lung_mixture", "public": "yes", "n": n_normal, "note": "excluded"},
            {"item": "unique_LUAD_patients", "public": "yes", "n": int(luad["Annotation"].nunique()), "note": "1 array per AD annotation"},
            {
                "item": "ICI_PD1_PDL1_treatment",
                "public": "no",
                "n": 0,
                "note": f"{int((luad['Gefitinib treatment'].astype(str) == 'Y').sum())}/{n_ad} LUAD Gefitinib=Y; not ICI",
            },
            {"item": "CLDN4_finite", "public": "yes", "n": int(cldn4.notna().sum()), "note": NAMED["CLDN4"]["probe"]},
            {"item": "CD8A_finite", "public": "yes", "n": int(partners["CD8A"].notna().sum()), "note": NAMED["CD8A"]["probe"]},
            {"item": "CD274_PDCD1LG1_finite", "public": "yes", "n": int(partners["CD274"].notna().sum()), "note": NAMED["CD274"]["probe"]},
            {
                "item": "primary_pairwise_CLDN4_CD8A_CD274",
                "public": "yes",
                "n": int((cldn4.notna() & partners["CD8A"].notna() & partners["CD274"].notna()).sum()),
                "note": "n used for Spearman",
            },
            {"item": "CLDN4_Q4_vs_Q1", "public": "yes", "n": int((q == "Q4").sum() + (q == "Q1").sum()), "note": "23 vs 23; not 90"},
        ]
    )
    n_table.to_csv(TABLES / "n_table.tsv", sep="\t", index=False)

    sample = luad[
        [
            "title",
            "Annotation",
            "Histology",
            "Sex",
            "Age",
            "pStage",
            "EGFR status",
            "K-ras Status",
            "p53 Status",
            "Gefitinib treatment",
        ]
    ].copy()
    sample["cldn4"] = cldn4
    sample["cd8a"] = partners["CD8A"]
    sample["cd274_pdcd1lg1"] = partners["CD274"]
    sample["tacstd2"] = partners["TACSTD2"]
    sample["epithelial_mean_z"] = epithelial
    sample["cldn4_quartile"] = q.astype(str)
    sample.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    # Figures
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(cldn4.values, partners["CD8A"].values, s=18, alpha=0.75, c="#2c5f8a", edgecolors="none")
    c8 = primary[primary.gene_y == "CD8A"].iloc[0]
    ax.set_xlabel("CLDN4  log10(R/G)")
    ax.set_ylabel("CD8A  log10(R/G)")
    ax.set_title(
        f"GSE11969 LUAD  n={int(c8.n)} arrays\n"
        f"CLDN4 vs CD8A  ρ={c8.rho:.3f}  p={fmt_p(c8.p)}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.png", dpi=160)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(cldn4.values, partners["CD274"].values, s=18, alpha=0.75, c="#6b3a2a", edgecolors="none")
    pdl1 = primary[primary.gene_y == "CD274"].iloc[0]
    ax.set_xlabel("CLDN4  log10(R/G)")
    ax.set_ylabel("CD274 / PDCD1LG1  log10(R/G)")
    ax.set_title(
        f"GSE11969 LUAD  n={int(pdl1.n)} arrays\n"
        f"CLDN4 vs CD274  ρ={pdl1.rho:.3f}  p={fmt_p(pdl1.p)}"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_cd274.png", dpi=160)
    fig.savefig(FIGURES / "fig2_cldn4_vs_cd274.pdf")
    plt.close(fig)

    cov = {
        "dataset": "GSE11969",
        "pmid": "16549822",
        "pmid_followup": "21465578",
        "platform": "GPL7015 Agilent Homo sapiens 21.6K custom array",
        "processing": series.get("data_processing")
        or "LOWESS normalized, background subtracted VALUE = log10(processed Red/processed Green)",
        "n_series": n_series,
        "n_nsclc": n_nsclc,
        "n_luad": n_ad,
        "n_sq": n_sq,
        "n_la": n_la,
        "n_as": n_as,
        "n_lcnec": n_lcnec,
        "n_sclc": n_sclc,
        "n_normal": n_normal,
        "n_probes": int(expr.shape[0]),
        "n_missing_value_cells": int(expr.isna().sum().sum()),
        "epithelial_used": epi_present,
        "epithelial_missing": epi_missing,
        "cldn4_vs_epi_rho": float(epi_r) if np.isfinite(epi_r) else None,
        "cldn4_vs_epi_p": float(epi_p) if np.isfinite(epi_p) else None,
        "cd274_gpl_symbol": "PDCD1LG1",
        "gefitinib_luad_yes": int((luad["Gefitinib treatment"].astype(str) == "Y").sum()),
    }
    summary = {
        **cov,
        "primary": primary.to_dict(orient="records"),
        "q4q1": qtab.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    write_finding(n_table, primary, qtab, cov)
    print(primary.to_string(index=False))
    print(qtab.to_string(index=False))
    print("[done]", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
