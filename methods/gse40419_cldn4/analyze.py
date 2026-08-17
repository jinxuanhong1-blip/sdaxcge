#!/usr/bin/env python3
"""GSE40419 Korean LUAD RNA-seq (Seo): CLDN4-only vs CD8A / CD274 / ImmuneScore / IFN.

Additive public processed-matrix slice. No FASTQ. No dual-high.
Tumor only. Epithelial residual. Honest n from the deposited RPKM file.

The GEO series matrix has 0 expression rows. Expression is the author
RPKM supplementary file (GSE40419_LC-87_RPKM_expression.txt.gz, ~10 MB).

Downloads stay under $GSE40419_CLDN4_DATA (default /tmp/gse40419_cldn4)
and are not committed.
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
DATA = Path(os.environ.get("GSE40419_CLDN4_DATA", "/tmp/gse40419_cldn4"))
SIG = HERE / "estimate_yoshihara_2013.tsv"

SEED = 20260817
N_BOOT = 2000
HOLDS_N = 40

UA = (
    "sdaxcge-gse40419-cldn4/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40419/"
    "matrix/GSE40419_series_matrix.txt.gz"
)
GEO_RPKM = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40419/"
    "suppl/GSE40419_LC-87_RPKM_expression.txt.gz"
)

EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
# Ayers et al. J Clin Invest 2017 IFNG 6-gene (PMID 28240696).
AYERS_IFNG6 = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
# Compact type-I ISG core (sensitivity; not the primary IFN row).
ISG_TYPE1 = [
    "ISG15",
    "MX1",
    "OAS1",
    "OAS2",
    "IFIT1",
    "IFIT3",
    "STAT1",
    "IRF7",
    "IFI27",
    "IFI44L",
]
NAMED = ["CLDN4", "CD8A", "CD274", "TACSTD2"]
PRIMARY_ENDS = ["CD8A", "CD274", "ImmuneScore", "IFN_Ayers6"]


def download(url: str, dest: Path, min_bytes: int = 1000) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
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


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
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
        if np.unique(xx[idx]).size < 3 or np.unique(yy[idx]).size < 3:
            continue
        r, _ = stats.spearmanr(xx[idx], yy[idx])
        if np.isfinite(r):
            rhos.append(float(r))
    if len(rhos) < 50:
        return np.nan, np.nan
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return float(lo), float(hi)


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
        return np.nan, np.nan, n
    xr = stats.rankdata(arr[:, 0])
    yr = stats.rankdata(arr[:, 1])
    Z = np.column_stack([stats.rankdata(arr[:, 2 + i]) for i in range(k)])
    rx = residualize(xr, Z)
    ry = residualize(yr, Z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan, n
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    dof = n - 2 - k
    t = r * math.sqrt(dof / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), dof))
    return r, p, n


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


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    z = expr.loc[present].apply(lambda s: (s - s.mean()) / (s.std(ddof=0) or np.nan), axis=1)
    return z.mean(axis=0), present


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict, int]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[str]] = {}
    n_expr_rows = 0
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                if line.startswith("ID_REF") or line.startswith('"ID_REF"'):
                    continue
                if line.strip():
                    n_expr_rows += 1
                continue
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                if key in sample_fields:
                    n = 2
                    k2 = f"{key}_{n}"
                    while k2 in sample_fields:
                        n += 1
                        k2 = f"{key}_{n}"
                    sample_fields[k2] = vals
                else:
                    sample_fields[key] = vals
    meta = pd.DataFrame(sample_fields)
    if "geo_accession" not in meta.columns:
        raise KeyError(f"no geo_accession: {list(meta.columns)}")
    meta.index = meta["geo_accession"].astype(str)
    meta.index.name = "gsm"
    return meta, {k: " | ".join(v) for k, v in series.items()}, n_expr_rows


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
        "n": nu,
        "rho": ru,
        "p": pu,
        "ci_low": lo,
        "ci_high": hi,
        "n_partial": np_,
        "rho_adj_epi": rp,
        "p_adj_epi": pp,
        "verdict": verdict(rp, pp, np_) if note == "primary" else note,
        "note": note,
    }


def pick_row(df: pd.DataFrame, cohort: str, pred: str, end: str):
    hit = df[(df.cohort == cohort) & (df.predictor == pred) & (df.endpoint == end)]
    return hit.iloc[0] if len(hit) else None


def scatter(ax, x, y, xlab, ylab, n, rho, p):
    ax.scatter(x, y, s=22, alpha=0.65, c="#1f4e79", edgecolors="none")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
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


def write_finding(stats_df, cov, highlow, inventory) -> None:
    def cell(cohort, pred, end, adj=False):
        r = pick_row(stats_df, cohort, pred, end)
        if r is None:
            return "—"
        if adj:
            return f"{fmt_rho(r.rho_adj_epi)} ({fmt_p(r.p_adj_epi)})"
        return f"{fmt_rho(r.rho)} ({fmt_p(r.p)})"

    def row(cohort, pred, end):
        return pick_row(stats_df, cohort, pred, end)

    c8 = row("LUAD_tumor", "CLDN4", "CD8A")
    pdl1 = row("LUAD_tumor", "CLDN4", "CD274")
    imm = row("LUAD_tumor", "CLDN4", "ImmuneScore")
    ifn = row("LUAD_tumor", "CLDN4", "IFN_Ayers6")
    isg = row("LUAD_tumor", "CLDN4", "IFN_ISG_type1")
    cd8_nr = row("LUAD_tumor", "CLDN4", "CD8A_NR027353")
    epi = row("LUAD_tumor", "CLDN4", "epithelial_z")
    tac_cd8 = row("LUAD_tumor", "TACSTD2", "CD8A")
    tac_pdl1 = row("LUAD_tumor", "TACSTD2", "CD274")
    c4_tac = row("LUAD_tumor", "CLDN4", "TACSTD2")
    cd8_imm = row("LUAD_tumor", "CD8A", "ImmuneScore")
    ifn_cd8 = row("LUAD_tumor", "IFN_Ayers6", "CD8A")

    n = cov["n_tumor"]
    overall = "NO_EVIDENCE"
    holds_any = any(
        r is not None and r.verdict == "HOLDS" for r in (c8, pdl1, imm, ifn)
    )
    opp_any = any(
        r is not None and r.verdict == "OPPOSITE" for r in (c8, pdl1, imm, ifn)
    )
    if holds_any and not opp_any:
        overall = "HOLDS (subset of endpoints)"
    elif opp_any and not holds_any:
        overall = "OPPOSITE (subset of endpoints)"
    elif holds_any and opp_any:
        overall = "MIXED"
    elif all(r is not None and r.verdict == "NO_EVIDENCE" for r in (c8, pdl1, imm, ifn)):
        overall = "NO_EVIDENCE"

    def hl(end):
        hit = highlow[highlow.endpoint == end]
        return hit.iloc[0] if len(hit) else None

    md = f"""# GSE40419 Korean LUAD RNA-seq (Seo) — CLDN4-only vs CD8A / CD274 / ImmuneScore / IFN

**Additive only. CLDN4-only. No dual-high.** Public Seo Korean lung adenocarcinoma RNA-seq (Seo et al., *Genome Res* 2012, PMID 22975805; GEO [GSE40419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419)). Unit is the **tumor RNA-seq library**. No slide was re-scored. No ICI arm. TACSTD2 is a same-run companion only and is **never a gate**.

The GEO series matrix has **0 expression rows** (`!Sample_data_row_count = 0`). Expression is the author processed RPKM file `GSE40419_LC-87_RPKM_expression.txt.gz` (~10 MB, under 2 GB). No FASTQ was used.

Primary question: **CLDN4 vs CD8A / CD274 / ImmuneScore / IFN** on **tumors only**, with an **epithelial residual**.

## Honest n

Do not write n=200 (paper clinical cohort before RNA-seq subset). Do not write n=164 (87 tumors + 77 adjacent normal). Do not write an expression n from the empty series matrix.

| Item | Public? | n | Note |
|---|---|---:|---|
| GSM in the series matrix | yes | **{cov['n_gsm']}** | 87 tumor + 77 adjacent normal titles |
| Expression rows in the series matrix | yes | **{cov['n_matrix_expr_rows']}** | empty table; `data_row_count=0` |
| Genes × libraries in author RPKM | yes | **{cov['n_rpkm_rows']} × {cov['n_rpkm_libs']}** | RefSeq rows; 87 tumor + 77 `*_nor` |
| Unique HUGO after max-mean | yes | **{cov['n_hugo']}** | isoform collapse on the tumor matrix |
| Tumors (`Lung cancer cells`) | yes | **{cov['n_tumor']}** | titles `LC_C*` / `LC_S*`; all unique |
| Adjacent normal (`*_nor`) | yes | **{cov['n_normal']}** | **dropped** (tumor-only slice) |
| Tumors with a paired `*_nor` | yes | **{cov['n_paired']}** | pairing not used |
| Unique tumor titles / patients | yes | **{cov['n_tumor']}** | 1 title = 1 library; no duplicate-patient sentence |
| Paper text “200 cancer patients” | paper only | 200 | broader clinical set; **not** the RNA-seq n |
| Stage on tumors | yes | **{cov['n_tumor_stage']}** | GEO `Stage`; {cov['n_tumor_stage_na']} tumor NA |
| Smoking on tumors | yes | **{cov['n_tumor_smoke']}** | never / smoker / current / NA |
| OS / ICI / response | no | **0** | surgical atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | **0** | not deposited |
| CLDN4 finite | yes | **{cov['n_cldn4']}** | single RefSeq row `{cov['cldn4_acc']}` |
| CD8A finite | yes | **{cov['n_cd8a']}** | protein-coding `{cov['cd8a_acc']}` (max-mean among NM_); {cov['n_cd8a_iso']} isoforms on the locus |
| CD274 finite | yes | **{cov['n_cd274']}** | single RefSeq row `{cov['cd274_acc']}` |
| ImmuneScore (Yoshihara Immune141) | computed | **{cov['n_tumor']}** | ssGSEA; {cov['immune_present']}/141 genes present |
| IFN Ayers-6 | computed | **{cov['n_tumor']}** | {cov['ayers_present']}/6 genes present |
| Epithelial mean-z (6/6) | computed | **{cov['n_tumor']}** | EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |
| **Primary pairwise n (tumor LUAD)** | yes | **{n}** | complete-case CLDN4 + CD8A + CD274 + scores |

The computable public n is **{n} tumor libraries**. Use 164 only when counting GSM metadata. Use 77 only for the unused adjacent-normal columns.

## One-row table

| dataset | histology | platform | n tumors | CLDN4–CD8A ρ (p) | epi adj ρ (p) | CLDN4–CD274 ρ (p) | epi adj ρ (p) | CLDN4–ImmuneScore ρ (p) | epi adj ρ (p) | CLDN4–IFN ρ (p) | epi adj ρ (p) | verdict | OS / ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE40419 Seo | LUAD | HiSeq 2000 author RPKM | **{n}** | {cell('LUAD_tumor','CLDN4','CD8A')} | {cell('LUAD_tumor','CLDN4','CD8A', True)} | {cell('LUAD_tumor','CLDN4','CD274')} | {cell('LUAD_tumor','CLDN4','CD274', True)} | {cell('LUAD_tumor','CLDN4','ImmuneScore')} | {cell('LUAD_tumor','CLDN4','ImmuneScore', True)} | {cell('LUAD_tumor','CLDN4','IFN_Ayers6')} | {cell('LUAD_tumor','CLDN4','IFN_Ayers6', True)} | **{overall}** | not deposited |

Full numbers: `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore / IFN (primary)

Author RPKM as deposited. Spearman is rank-based (log2 does not change ρ). Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **{cov['epi_present']}/6** present). HOLDS rule (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05.

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | verdict |
|---|---|---:|---:|---|---:|---|---|
| CLDN4 vs **CD8A** | LUAD tumor | **{c8['n']}** | **{fmt_rho(c8['rho'])}** | {fmt_ci(c8['ci_low'], c8['ci_high'])} | {fmt_p(c8['p'])} | {fmt_rho(c8['rho_adj_epi'])} ({fmt_p(c8['p_adj_epi'])}) | **{c8['verdict']}** |
| CLDN4 vs **CD274** | LUAD tumor | **{pdl1['n']}** | **{fmt_rho(pdl1['rho'])}** | {fmt_ci(pdl1['ci_low'], pdl1['ci_high'])} | {fmt_p(pdl1['p'])} | {fmt_rho(pdl1['rho_adj_epi'])} ({fmt_p(pdl1['p_adj_epi'])}) | **{pdl1['verdict']}** |
| CLDN4 vs **ImmuneScore** | LUAD tumor | **{imm['n']}** | **{fmt_rho(imm['rho'])}** | {fmt_ci(imm['ci_low'], imm['ci_high'])} | {fmt_p(imm['p'])} | {fmt_rho(imm['rho_adj_epi'])} ({fmt_p(imm['p_adj_epi'])}) | **{imm['verdict']}** |
| CLDN4 vs **IFN (Ayers-6)** | LUAD tumor | **{ifn['n']}** | **{fmt_rho(ifn['rho'])}** | {fmt_ci(ifn['ci_low'], ifn['ci_high'])} | {fmt_p(ifn['p'])} | {fmt_rho(ifn['rho_adj_epi'])} ({fmt_p(ifn['p_adj_epi'])}) | **{ifn['verdict']}** |
| CLDN4 vs IFN type-I ISG (sensitivity) | LUAD tumor | {isg['n']} | {fmt_rho(isg['rho'])} | {fmt_ci(isg['ci_low'], isg['ci_high'])} | {fmt_p(isg['p'])} | {fmt_rho(isg['rho_adj_epi'])} ({fmt_p(isg['p_adj_epi'])}) | {isg['verdict']} |
| CLDN4 vs CD8A `NR_027353` (sensitivity) | LUAD tumor | {cd8_nr['n']} | {fmt_rho(cd8_nr['rho'])} | {fmt_ci(cd8_nr['ci_low'], cd8_nr['ci_high'])} | {fmt_p(cd8_nr['p'])} | {fmt_rho(cd8_nr['rho_adj_epi'])} ({fmt_p(cd8_nr['p_adj_epi'])}) | {cd8_nr['verdict']} |

Positive-control (same tumors): CD8A vs ImmuneScore ρ={fmt_rho(cd8_imm['rho'])} (p={fmt_p(cd8_imm['p'])}); IFN Ayers-6 vs CD8A ρ={fmt_rho(ifn_cd8['rho'])} (p={fmt_p(ifn_cd8['p'])}).

CLDN4 vs epithelial mean-z: ρ={fmt_rho(epi['rho'])} (p={fmt_p(epi['p'])}, n={epi['n']}). The residual asks whether any immune association is more than epithelial content.

**What is inverse unadjusted.** CLDN4 vs CD8A and vs ImmuneScore are negative at n=87 before the residual (CIs exclude 0). Q4 vs Q1 (22 vs 22) matches that direction for CD8A and ImmuneScore. That is a coarsened unadjusted test, not a residual claim.

**What does not HOLDS.** After the specified epithelial residual, every primary pair is NS (CD8A adj p={fmt_p(c8['p_adj_epi'])}; ImmuneScore adj p={fmt_p(imm['p_adj_epi'])}; CD274 and IFN already NS unadjusted). HOLDS requires ρ_adj<0 and p_adj<0.05. Do not upgrade the unadjusted inverses into an epithelial-independent immune-low claim.

### Q4 vs Q1 (descriptive)

Quartiles are `pd.qcut(rank(method="first"), 4)` on the {n} tumors with finite CLDN4. Honest n is the quartile arms, not {n}.

| endpoint | n Q4 vs Q1 | MWU p | rank-biserial (Q4>Q1) |
|---|---|---:|---:|
| CD8A | {int(hl('CD8A')['n_q4'])} vs {int(hl('CD8A')['n_q1'])} | {fmt_p(hl('CD8A')['p'])} | {hl('CD8A')['rank_biserial']:+.3f} |
| CD274 | {int(hl('CD274')['n_q4'])} vs {int(hl('CD274')['n_q1'])} | {fmt_p(hl('CD274')['p'])} | {hl('CD274')['rank_biserial']:+.3f} |
| ImmuneScore | {int(hl('ImmuneScore')['n_q4'])} vs {int(hl('ImmuneScore')['n_q1'])} | {fmt_p(hl('ImmuneScore')['p'])} | {hl('ImmuneScore')['rank_biserial']:+.3f} |
| IFN Ayers-6 | {int(hl('IFN_Ayers6')['n_q4'])} vs {int(hl('IFN_Ayers6')['n_q1'])} | {fmt_p(hl('IFN_Ayers6')['p'])} | {hl('IFN_Ayers6')['rank_biserial']:+.3f} |

## Context (not the claim)

This is a surgical Korean LUAD transcriptome atlas (Seo 2012). It is **not** an ICI-response series. Adjacent-normal columns exist in the RPKM file and were not used. Smoking and stage are deposited for tumors and were not residualised (the specified residual is epithelial). No dual-high (CLDN4-high AND TACSTD2-high) split was run.

## TACSTD2 companion (not the claim)

| pair | subset | n | ρ (p) | ρ_adj epi (p) | verdict |
|---|---|---:|---|---|---|
| TACSTD2 vs CD8A | LUAD tumor | {tac_cd8['n']} | {fmt_rho(tac_cd8['rho'])} ({fmt_p(tac_cd8['p'])}) | {fmt_rho(tac_cd8['rho_adj_epi'])} ({fmt_p(tac_cd8['p_adj_epi'])}) | {tac_cd8['verdict']} |
| TACSTD2 vs CD274 | LUAD tumor | {tac_pdl1['n']} | {fmt_rho(tac_pdl1['rho'])} ({fmt_p(tac_pdl1['p'])}) | {fmt_rho(tac_pdl1['rho_adj_epi'])} ({fmt_p(tac_pdl1['p_adj_epi'])}) | {tac_pdl1['verdict']} |
| CLDN4 vs TACSTD2 | LUAD tumor | {c4_tac['n']} | {fmt_rho(c4_tac['rho'])} ({fmt_p(c4_tac['p'])}) | {fmt_rho(c4_tac['rho_adj_epi'])} ({fmt_p(c4_tac['p_adj_epi'])}) | {c4_tac['verdict']} |

Do not treat CLDN4 as interchangeable with TACSTD2. Dual-high was not computed.

## Methods

- **Matrix:** GEO supplementary `GSE40419_LC-87_RPKM_expression.txt.gz` (author RPKM, NCBI build 37.1 / GSNAP). The series matrix expression table is empty and was used only for sample characteristics.
- **Tumor rule:** GEO `source_name = Lung cancer cells` (titles without `_nor`). Adjacent normal dropped.
- **Gene collapse:** max-mean of RefSeq rows that share a HUGO symbol, computed on the tumor matrix. Named genes (CLDN4 / CD8A / CD274 / TACSTD2) take the highest-mean **NM_** accession when one exists.
- **CD8** = protein-coding `CD8A` `NM_001768` (highest-mean NM_ on the locus). The non-coding `NR_027353` isoform is a sensitivity row only. **CD274** = PD-L1 RNA (not protein).
- **ImmuneScore:** Yoshihara 2013 Immune141 ssGSEA (Barbie/GSVA, τ=0.25); ranks scaled to 1…10000. Stromal141 is computed only to form ESTIMATEScore for the coverage audit; the specified residual is epithelial, not ESTIMATE purity.
- **IFN (primary):** Ayers 2017 IFNG 6-gene mean of gene-wise z (`IFNG`, `STAT1`, `CXCL9`, `CXCL10`, `IDO1`, `HLA-DRA`).
- **IFN (sensitivity):** 10-gene type-I ISG mean-z (`ISG15`, `MX1`, `OAS1`, `OAS2`, `IFIT1`, `IFIT3`, `STAT1`, `IRF7`, `IFI27`, `IFI44L`).
- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7.
- **Partial Spearman:** Pearson of rank residuals; df = n − 2 − k.
- **Not done:** no FASTQ / recount. No dual-high. No ICI / OS model. No gene-set fishing beyond the pre-specified IFN lists.

## Files

- `tables/spearman.tsv` — all n / ρ / p / CI / epithelial residuals
- `tables/label_inventory.tsv` — honest n
- `tables/sample_annotation.tsv` — per-tumor genes + scores + GEO labels
- `tables/gene_coverage.tsv` / `highlow_cldn4.tsv` / `isoform_used.tsv` (named + set genes only) / `summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`
- `figures/fig2_cldn4_vs_immunescore_ifn.png`
- `figures/fig3_forest.png`

## Reproduce

```bash
python3 -m pip install -r methods/gse40419_cldn4/requirements.txt
export GSE40419_CLDN4_DATA=/tmp/gse40419_cldn4
python3 methods/gse40419_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)
    print("wrote FINDING.md", flush=True)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = download(GEO_MATRIX, DATA / "GSE40419_series_matrix.txt.gz", min_bytes=1000)
    rpkm_path = download(GEO_RPKM, DATA / "GSE40419_LC-87_RPKM_expression.txt.gz", min_bytes=100_000)

    meta, series, n_matrix_expr = parse_series_meta(matrix_path)
    print(f"series matrix GSM={len(meta)} expr_rows={n_matrix_expr}", flush=True)

    # characteristics: age / gender / smoking / Stage (order in the SOFT dump)
    char_cols = [c for c in meta.columns if c.startswith("characteristics_ch1")]
    rename = {}
    for c in char_cols:
        vals = meta[c].astype(str)
        key = None
        for v in vals:
            if ":" in v:
                key = v.split(":", 1)[0].strip().lower().replace(" ", "_")
                break
        if key:
            rename[c] = key
            meta[c] = vals.map(lambda x: x.split(":", 1)[1].strip() if ":" in str(x) else str(x))
    meta = meta.rename(columns=rename)
    title_col = "title" if "title" in meta.columns else "title_1"
    source_col = "source_name_ch1" if "source_name_ch1" in meta.columns else "source_name_ch1_1"
    meta["title"] = meta[title_col].astype(str)
    meta["source"] = meta[source_col].astype(str)
    meta["is_tumor"] = meta["source"].str.contains("cancer", case=False, na=False) & ~meta[
        "title"
    ].str.endswith("_nor")
    meta["is_normal"] = meta["title"].str.endswith("_nor") | meta["source"].str.contains(
        "normal", case=False, na=False
    )

    raw = pd.read_csv(rpkm_path, sep="\t", compression="gzip", low_memory=False)
    if "gene" not in raw.columns:
        raise SystemExit(f"unexpected RPKM columns: {list(raw.columns)[:8]}")
    lib_cols = [c for c in raw.columns if c.startswith("LC_")]
    tumor_cols = [c for c in lib_cols if not c.endswith("_nor")]
    normal_cols = [c for c in lib_cols if c.endswith("_nor")]
    print(f"RPKM rows={len(raw)} libs={len(lib_cols)} tumor={len(tumor_cols)} normal={len(normal_cols)}", flush=True)

    # collapse isoforms on the tumor matrix (max-mean). Named immune genes
    # prefer a protein-coding NM_ accession when one exists (CD8A NR_027353
    # is a non-coding transcript with a slightly higher mean).
    tum = raw[["gene", "accession"] + tumor_cols].copy()
    tum[tumor_cols] = tum[tumor_cols].apply(pd.to_numeric, errors="coerce")
    tum["_mean"] = tum[tumor_cols].mean(axis=1)
    tum["_nm"] = tum["accession"].astype(str).str.startswith("NM_")
    tum = tum.sort_values(["_nm", "_mean"], ascending=[False, False])
    named_pref = tum[tum["gene"].isin(NAMED)].copy()
    named_pref = named_pref.loc[~named_pref["gene"].duplicated(keep="first")]
    rest = tum.loc[~tum["gene"].isin(NAMED)].sort_values("_mean", ascending=False)
    rest = rest.loc[~rest["gene"].duplicated(keep="first")]
    kept = pd.concat([named_pref, rest], axis=0)
    gene_expr = kept.set_index("gene")[tumor_cols].astype(float)
    isoform_used = kept[["gene", "accession", "_mean"]].rename(columns={"_mean": "tumor_mean_rpkm"})
    cd8a_nr = tum.loc[tum["accession"] == "NR_027353", tumor_cols]
    cd8a_nr_s = (
        cd8a_nr.iloc[0].astype(float)
        if len(cd8a_nr)
        else pd.Series(np.nan, index=tumor_cols)
    )

    # map GEO labels onto RPKM titles
    title2meta = meta.set_index("title")
    missing_titles = [c for c in tumor_cols if c not in title2meta.index]
    if missing_titles:
        print(f"WARNING titles not in GEO meta: {missing_titles[:8]}", flush=True)

    sig = pd.read_csv(SIG, sep="\t")
    stromal = sig.loc[sig["set"] == "Stromal141_UP", "hugo"].astype(str).tolist()
    immune = sig.loc[sig["set"] == "Immune141_UP", "hugo"].astype(str).tolist()
    immune_score = ssgsea(gene_expr, immune)
    stromal_score = ssgsea(gene_expr, stromal)
    epi_z, epi_present = mean_z(gene_expr, EPITHELIAL)
    ifn6, ayers_present = mean_z(gene_expr, AYERS_IFNG6)
    isg, isg_present = mean_z(gene_expr, ISG_TYPE1)

    def g(name: str) -> pd.Series:
        if name not in gene_expr.index:
            return pd.Series(np.nan, index=gene_expr.columns)
        return gene_expr.loc[name]

    cldn4 = g("CLDN4")
    cd8a = g("CD8A")
    cd274 = g("CD274")
    tac = g("TACSTD2")
    log_c4 = np.log2(cldn4 + 1)
    log_cd8 = np.log2(cd8a + 1)
    log_pdl1 = np.log2(cd274 + 1)

    rows = []
    z = epi_z.to_numpy(float)
    rows.append(rec_row("LUAD_tumor", "CLDN4", "CD8A", cldn4, cd8a, z, "primary"))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "CD274", cldn4, cd274, z, "primary"))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "ImmuneScore", cldn4, immune_score, z, "primary"))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "IFN_Ayers6", cldn4, ifn6, z, "primary"))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "IFN_ISG_type1", cldn4, isg, z, "sensitivity"))
    rows.append(
        rec_row(
            "LUAD_tumor",
            "CLDN4",
            "CD8A_NR027353",
            cldn4,
            cd8a_nr_s.reindex(cldn4.index),
            z,
            "sensitivity",
        )
    )
    rows.append(rec_row("LUAD_tumor", "CLDN4", "epithelial_z", cldn4, epi_z, None, "context"))
    rows.append(rec_row("LUAD_tumor", "TACSTD2", "CD8A", tac, cd8a, z, "companion"))
    rows.append(rec_row("LUAD_tumor", "TACSTD2", "CD274", tac, cd274, z, "companion"))
    rows.append(rec_row("LUAD_tumor", "TACSTD2", "ImmuneScore", tac, immune_score, z, "companion"))
    rows.append(rec_row("LUAD_tumor", "TACSTD2", "IFN_Ayers6", tac, ifn6, z, "companion"))
    rows.append(rec_row("LUAD_tumor", "CLDN4", "TACSTD2", cldn4, tac, z, "companion"))
    rows.append(rec_row("LUAD_tumor", "CD8A", "ImmuneScore", cd8a, immune_score, None, "positive_control"))
    rows.append(rec_row("LUAD_tumor", "IFN_Ayers6", "CD8A", ifn6, cd8a, None, "positive_control"))
    rows.append(rec_row("LUAD_tumor", "IFN_Ayers6", "ImmuneScore", ifn6, immune_score, None, "positive_control"))
    stats_df = pd.DataFrame(rows)

    # Q4 vs Q1
    ranks = cldn4.rank(method="first")
    q = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    hl_rows = []
    for end_name, series in [
        ("CD8A", cd8a),
        ("CD274", cd274),
        ("ImmuneScore", immune_score),
        ("IFN_Ayers6", ifn6),
        ("TACSTD2", tac),
    ]:
        a = series[q == "Q4"].to_numpy(float)
        b = series[q == "Q1"].to_numpy(float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        n4, n1 = int(a.size), int(b.size)
        if n4 >= 3 and n1 >= 3:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            r = 2.0 * float(u) / (n4 * n1) - 1.0
        else:
            u, p, r = np.nan, np.nan, np.nan
        hl_rows.append(
            {
                "endpoint": end_name,
                "n_q4": n4,
                "n_q1": n1,
                "median_q4": float(np.median(a)) if n4 else np.nan,
                "median_q1": float(np.median(b)) if n1 else np.nan,
                "delta_q4_minus_q1": (float(np.median(a) - np.median(b)) if n4 and n1 else np.nan),
                "U": float(u) if np.isfinite(u) else np.nan,
                "p": p,
                "rank_biserial": r,
            }
        )
    highlow = pd.DataFrame(hl_rows)

    # annotation
    ann = pd.DataFrame(
        {
            "title": tumor_cols,
            "CLDN4": cldn4.values,
            "CD8A": cd8a.values,
            "CD274": cd274.values,
            "TACSTD2": tac.values,
            "ImmuneScore": immune_score.reindex(tumor_cols).values,
            "StromalScore": stromal_score.reindex(tumor_cols).values,
            "IFN_Ayers6": ifn6.reindex(tumor_cols).values,
            "IFN_ISG_type1": isg.reindex(tumor_cols).values,
            "epithelial_z": epi_z.reindex(tumor_cols).values,
            "CLDN4_quartile": q.astype(str).values,
            "log2_CLDN4": log_c4.values,
            "log2_CD8A": log_cd8.values,
            "log2_CD274": log_pdl1.values,
        }
    )
    for col in ["geo_accession", "source", "age_at_diagnosis", "gender", "smoking_status", "stage"]:
        src_col = col if col != "stage" else ("stage" if "stage" in title2meta.columns else None)
        if col == "stage":
            src_col = "stage" if "stage" in title2meta.columns else None
        if col == "geo_accession":
            src_col = "geo_accession"
        if src_col and src_col in title2meta.columns:
            ann[col] = [title2meta.loc[t, src_col] if t in title2meta.index else np.nan for t in tumor_cols]
        elif col == "stage" and "Stage" in title2meta.columns:
            ann["stage"] = [title2meta.loc[t, "Stage"] if t in title2meta.index else np.nan for t in tumor_cols]

    # coverage / honest n
    tum_meta = meta.loc[meta["is_tumor"]].copy()
    nor_meta = meta.loc[meta["is_normal"]].copy()
    stage_col = "stage" if "stage" in tum_meta.columns else None
    smoke_col = "smoking_status" if "smoking_status" in tum_meta.columns else None

    def n_non_na(s):
        if s is None:
            return 0
        return int((~s.astype(str).isin(["", "NA", "nan", "None"])).sum())

    n_cd8a_iso = int((raw["gene"] == "CD8A").sum())
    n_krt18_iso = int((raw["gene"] == "KRT18").sum())
    acc = isoform_used.set_index("gene")["accession"]

    cov = {
        "n_gsm": int(len(meta)),
        "n_matrix_expr_rows": int(n_matrix_expr),
        "n_rpkm_rows": int(len(raw)),
        "n_rpkm_libs": int(len(lib_cols)),
        "n_hugo": int(gene_expr.shape[0]),
        "n_tumor": int(len(tumor_cols)),
        "n_normal": int(len(normal_cols)),
        "n_paired": int(sum((t + "_nor") in set(normal_cols) for t in tumor_cols)),
        "n_tumor_stage": n_non_na(tum_meta[stage_col]) if stage_col else 0,
        "n_tumor_stage_na": int(len(tum_meta) - (n_non_na(tum_meta[stage_col]) if stage_col else 0)),
        "n_tumor_smoke": n_non_na(tum_meta[smoke_col]) if smoke_col else 0,
        "n_cldn4": int(np.isfinite(cldn4).sum()),
        "n_cd8a": int(np.isfinite(cd8a).sum()),
        "n_cd274": int(np.isfinite(cd274).sum()),
        "n_cd8a_iso": n_cd8a_iso,
        "n_krt18_iso": n_krt18_iso,
        "cldn4_acc": str(acc.get("CLDN4", "NA")),
        "cd8a_acc": str(acc.get("CD8A", "NA")),
        "cd274_acc": str(acc.get("CD274", "NA")),
        "immune_present": int(sum(g in gene_expr.index for g in immune)),
        "stromal_present": int(sum(g in gene_expr.index for g in stromal)),
        "ayers_present": int(len(ayers_present)),
        "isg_present": int(len(isg_present)),
        "epi_present": int(len(epi_present)),
        "series_title": series.get("title", "GSE40419"),
        "series_design": series.get("overall_design", ""),
    }

    inventory = pd.DataFrame(
        [
            {"item": "GSM in series matrix", "n": cov["n_gsm"], "note": "87 tumor + 77 adjacent normal"},
            {"item": "series matrix expression rows", "n": cov["n_matrix_expr_rows"], "note": "empty; do not use"},
            {"item": "RPKM gene rows", "n": cov["n_rpkm_rows"], "note": "RefSeq; before collapse"},
            {"item": "RPKM libraries", "n": cov["n_rpkm_libs"], "note": "87 tumor + 77 normal"},
            {"item": "HUGO after max-mean (tumor)", "n": cov["n_hugo"], "note": "analysis matrix"},
            {"item": "tumors (primary n)", "n": cov["n_tumor"], "note": "Lung cancer cells; unique titles"},
            {"item": "adjacent normal (dropped)", "n": cov["n_normal"], "note": "tumor-only slice"},
            {"item": "tumors with paired normal", "n": cov["n_paired"], "note": "pairing unused"},
            {"item": "paper clinical n=200", "n": 200, "note": "NOT the RNA-seq n; do not write"},
            {"item": "CLDN4 finite", "n": cov["n_cldn4"], "note": cov["cldn4_acc"]},
            {"item": "CD8A finite", "n": cov["n_cd8a"], "note": cov["cd8a_acc"]},
            {"item": "CD274 finite", "n": cov["n_cd274"], "note": cov["cd274_acc"]},
            {"item": "Immune141 present", "n": cov["immune_present"], "note": "of 141"},
            {"item": "Ayers IFNG-6 present", "n": cov["ayers_present"], "note": "of 6"},
            {"item": "epithelial 6 present", "n": cov["epi_present"], "note": "of 6"},
            {"item": "OS / ICI labels", "n": 0, "note": "not deposited"},
        ]
    )

    coverage = pd.DataFrame(
        [
            {"set": "epithelial", "n_list": 6, "n_present": cov["epi_present"], "genes": ",".join(epi_present)},
            {"set": "Ayers_IFNG6", "n_list": 6, "n_present": cov["ayers_present"], "genes": ",".join(ayers_present)},
            {"set": "ISG_type1", "n_list": 10, "n_present": cov["isg_present"], "genes": ",".join(isg_present)},
            {"set": "Immune141", "n_list": 141, "n_present": cov["immune_present"], "genes": ""},
            {"set": "Stromal141", "n_list": 141, "n_present": cov["stromal_present"], "genes": ""},
        ]
    )

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))
    r_cd8 = pick_row(stats_df, "LUAD_tumor", "CLDN4", "CD8A")
    r_pdl1 = pick_row(stats_df, "LUAD_tumor", "CLDN4", "CD274")
    scatter(
        axes[0],
        log_c4,
        log_cd8,
        "CLDN4 log2(RPKM+1)",
        "CD8A log2(RPKM+1)",
        int(r_cd8["n"]),
        r_cd8["rho"],
        r_cd8["p"],
    )
    scatter(
        axes[1],
        log_c4,
        log_pdl1,
        "CLDN4 log2(RPKM+1)",
        "CD274 log2(RPKM+1)",
        int(r_pdl1["n"]),
        r_pdl1["rho"],
        r_pdl1["p"],
    )
    fig.suptitle(f"GSE40419 Seo LUAD tumors, n={len(tumor_cols)} (CLDN4-only)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.png", dpi=140)
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a_cd274.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))
    scatter(
        axes[0],
        log_c4,
        immune_score.reindex(tumor_cols),
        "CLDN4 log2(RPKM+1)",
        "ImmuneScore (ssGSEA)",
        cldn4.size,
        pick_row(stats_df, "LUAD_tumor", "CLDN4", "ImmuneScore")["rho"],
        pick_row(stats_df, "LUAD_tumor", "CLDN4", "ImmuneScore")["p"],
    )
    scatter(
        axes[1],
        log_c4,
        ifn6.reindex(tumor_cols),
        "CLDN4 log2(RPKM+1)",
        "IFN Ayers-6 (mean-z)",
        cldn4.size,
        pick_row(stats_df, "LUAD_tumor", "CLDN4", "IFN_Ayers6")["rho"],
        pick_row(stats_df, "LUAD_tumor", "CLDN4", "IFN_Ayers6")["p"],
    )
    fig.suptitle(f"GSE40419 Seo LUAD tumors, n={len(tumor_cols)} (CLDN4-only)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_immunescore_ifn.png", dpi=140)
    fig.savefig(FIGURES / "fig2_cldn4_vs_immunescore_ifn.pdf")
    plt.close(fig)

    forest_ends = ["CD8A", "CD274", "ImmuneScore", "IFN_Ayers6"]
    labels = []
    xs, xlo, xhi = [], [], []
    for end in forest_ends:
        r = pick_row(stats_df, "LUAD_tumor", "CLDN4", end)
        labels.append(f"CLDN4 vs {end}\nunadjusted")
        xs.append(r["rho"])
        xlo.append(r["ci_low"])
        xhi.append(r["ci_high"])
        labels.append(f"CLDN4 vs {end}\nepi residual")
        xs.append(r["rho_adj_epi"])
        xlo.append(np.nan)
        xhi.append(np.nan)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    y = np.arange(len(labels))[::-1]
    for i, (lab, x, lo, hi) in enumerate(zip(labels, xs, xlo, xhi)):
        yy = y[i]
        color = "#1f4e79" if "unadjusted" in lab else "#b85c38"
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [yy, yy], color=color, lw=1.6)
        ax.plot(x, yy, "o", color=color, ms=6)
    ax.axvline(0, color="#666", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI on unadjusted)")
    ax.set_title(f"GSE40419 CLDN4-only, n={len(tumor_cols)} tumors")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_forest.png", dpi=140)
    fig.savefig(FIGURES / "fig3_forest.pdf")
    plt.close(fig)

    stats_df.to_csv(TABLES / "spearman.tsv", sep="\t", index=False)
    inventory.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)
    ann.to_csv(TABLES / "sample_annotation.tsv", sep="\t", index=False)
    coverage.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)
    highlow.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)
    keep_genes = set(NAMED + EPITHELIAL + AYERS_IFNG6 + ISG_TYPE1)
    named_iso = isoform_used[isoform_used["gene"].isin(keep_genes)].sort_values("gene")
    named_iso.to_csv(TABLES / "isoform_used.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE40419",
        "citation": "Seo et al. Genome Res 2012 PMID 22975805",
        "n_tumor": cov["n_tumor"],
        "n_normal_dropped": cov["n_normal"],
        "n_matrix_expr_rows": cov["n_matrix_expr_rows"],
        "primary": stats_df[stats_df.note == "primary"][
            ["endpoint", "n", "rho", "p", "rho_adj_epi", "p_adj_epi", "verdict"]
        ].to_dict(orient="records"),
        "coverage": cov,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    write_finding(stats_df, cov, highlow, inventory)
    print(stats_df[stats_df.note.isin(["primary", "positive_control"])].to_string(index=False), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
