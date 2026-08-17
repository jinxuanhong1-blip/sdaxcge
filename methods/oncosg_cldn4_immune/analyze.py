#!/usr/bin/env python3
"""OncoSG LUAD public RNA: CLDN4 vs CD8A / ImmuneScore / CD274.

Additive only. A9 extra-cohort 3/3 catalog (OncoSG, CPTAC LUAD RNA, GSE31210)
is taken as given. TACSTD2 vs CD8/GEP/immune on this same matrix is already
in PR 139 and is treated as given. This folder only measures CLDN4.

Public cBioPortal study luad_oncosg_2020 (Chen et al., Nat Genet 2020).
The open deposit is all-sample z-scores. ESTIMATE ImmuneScore / xCell /
MCP-counter are not computed. ImmuneScore here is the A1 8-gene T-cell
effector mean z (same definition as PR 139), plus the published IMSIG
T-cell column as a deposited immune score.

Honest n = 169 public z-score columns (portal RNA list is 181).
"""

from __future__ import annotations

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
CACHE = Path("/tmp/oncosg_cldn4_immune")

DATAHUB_SHA = "165bd77077b03038f9c2ee104959eb474770b2a9"
STUDY = "luad_oncosg_2020"
MEDIA = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/"
    f"{DATAHUB_SHA}/public/{STUDY}"
)
FILES = {
    "expression": (
        "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    ),
    "clinical_sample": (
        "data_clinical_sample.txt",
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
    ),
    "clinical_patient": (
        "data_clinical_patient.txt",
        None,
    ),
}

UA = "sdaxcge-oncosg-cldn4-immune/1.0 (+https://github.com/jinxuanhong1-blip/sdaxcge)"

# Portal RNA list IDs that are not columns in the public z-score file (PR 139).
RNA_LIST_NOT_IN_MATRIX = [
    "A008",
    "A114",
    "A122",
    "A136",
    "A139",
    "A147",
    "A184",
    "A302",
    "A435",
    "A484",
    "A489",
    "A507",
]

# Same A1 8-gene immune / T-cell effector list as PR 139. Not ESTIMATE.
IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
TJ_NO_CLDN4 = ["CLDN3", "CLDN7", "OCLN", "TJP1", "F11R", "MARVELD2", "CRB3", "PARD3", "CGN"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"cached {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    print(f"GET {url}\n -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=180) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    return dest


def ensure_file(key: str) -> Path:
    name, expected = FILES[key]
    dest = download(f"{MEDIA}/{name}", CACHE / name)
    if expected is not None:
        got = sha256_file(dest)
        if got != expected:
            raise SystemExit(f"checksum mismatch for {name}: {got}")
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
    """Pearson of rank residuals; df = n − 3. Fisher-z 95% CI, SE = 1/√(n−4)."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return {
            "n": n,
            "rho": np.nan,
            "p": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }
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
    zf = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 4)
    ci = tuple(float(np.tanh(zf + s * 1.959963984540054 * se)) for s in (-1.0, 1.0))
    return {"n": n, "rho": r, "p": p, "ci_low": ci[0], "ci_high": ci[1]}


def rec_row(predictor, endpoint, x, y, z, purity_name, note=""):
    r, p, n = spearman_pair(x, y)
    out = {
        "predictor": predictor,
        "endpoint": endpoint,
        "n": n,
        "unadj_rho": r,
        "unadj_p": p,
        "partial_rho": np.nan,
        "partial_p": np.nan,
        "partial_ci_low": np.nan,
        "partial_ci_high": np.nan,
        "purity": purity_name,
        "note": note,
    }
    if z is not None:
        part = partial_spearman(x, y, z)
        out["n"] = part["n"]
        out["partial_rho"] = part["rho"]
        out["partial_p"] = part["p"]
        out["partial_ci_low"] = part["ci_low"]
        out["partial_ci_high"] = part["ci_high"]
    return out


def mean_present(expr: pd.DataFrame, genes: list[str], samples: list[str]):
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), present, missing
    score = expr.loc[present, samples].astype(float).mean(axis=0)
    return score, present, missing


def write_finding(stats_df: pd.DataFrame, cov: dict) -> None:
    def g(pred, end, col, purity=None):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        if purity is not None:
            hit = hit[hit.purity == purity]
        if hit.empty:
            return np.nan
        return hit.iloc[0][col]

    n = cov["n_matrix"]
    pur = "published PURITY"

    def cell(pred, end):
        return f"{fmt_rho(g(pred, end, 'unadj_rho'))} ({fmt_p(g(pred, end, 'unadj_p'))})"

    def pcell(pred, end):
        return f"{fmt_rho(g(pred, end, 'partial_rho'))} ({fmt_p(g(pred, end, 'partial_p'))})"

    def row_md(pred, end, label=None):
        lab = label or f"{pred} vs {end}"
        return (
            f"| {lab} | {int(g(pred, end, 'n'))} | "
            f"{fmt_rho(g(pred, end, 'unadj_rho'))} | {fmt_p(g(pred, end, 'unadj_p'))} | "
            f"{fmt_rho(g(pred, end, 'partial_rho'))} | {fmt_p(g(pred, end, 'partial_p'))} |"
        )

    md = f"""# OncoSG LUAD — CLDN4 vs CD8A / ImmuneScore / CD274

**Additive only.** A9 extra-cohort 3/3 catalog (OncoSG LUAD, CPTAC LUAD RNA, GSE31210 LUAD) is **taken as given** and is not re-audited. TACSTD2 vs CD8 / GEP / immune on this same public matrix is already in [PR 139](https://github.com/jinxuanhong1-blip/sdaxcge/pull/139) and is treated as given. This folder only measures **CLDN4**.

Public East-Asian surgical LUAD (Chen et al., *Nat Genet* 2020; cBioPortal [`luad_oncosg_2020`](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)). Unit is the **public z-score column**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| cBioPortal RNA sample list `luad_oncosg_2020_rna_seq_v2_mrna` | 181 | portal description |
| public z-score matrix columns | **{n}** | `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` |
| RNA-list IDs absent from the public matrix | 12 | {", ".join(RNA_LIST_NOT_IN_MATRIX)} |
| clinical samples | {cov["clinical_samples"]} | `data_clinical_sample.txt` |
| clinical PURITY non-NA | {cov["clinical_purity_non_na"]} | includes many RNA-absent tumors |
| matrix samples with PURITY | **{n}** | all {n} columns have published PURITY |
| matrix samples with IMSIG T cells | **{n}** | published clinical column |
| CLDN4 / CD8A / CD274 finite | **{n}** | 0 missing z-scores |
| ONCOTREE LUAD on matrix | **{n}** | all columns are Lung Adenocarcinoma / Tumor |
| ESTIMATE ImmuneScore | **0** | skipped — deposit is z-scores, no raw RSEM |
| ICI / PD-1 labels | **0** | not an ICI-response cohort |

Primary tests use **n={n}**. Do not write n=181. The 12 portal-listed RNA IDs are not in the open expression table and cannot enter a CLDN4 correlation.

## One-row table

| dataset | n | CLDN4–CD8A ρ (p) | CLDN4–ImmuneScore ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–CD8A partial \\| PURITY (p) | CLDN4–ImmuneScore partial \\| PURITY (p) | CLDN4–CD274 partial \\| PURITY (p) |
|---|---:|---|---|---|---|---|---|
| OncoSG LUAD public RNA | {n} | {cell("CLDN4", "CD8A")} | {cell("CLDN4", "ImmuneScore")} | {cell("CLDN4", "CD274")} | {pcell("CLDN4", "CD8A")} | {pcell("CLDN4", "ImmuneScore")} | {pcell("CLDN4", "CD274")} |

ImmuneScore = A1 8-gene T-cell effector mean z (`CD8A GZMA GZMB IFNG EOMES CXCL9 CXCL10 TBX21`; 8/8 present). This is **not** Yoshihara ESTIMATE ImmuneScore. Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## ImmuneScore (honest skip of ESTIMATE)

cBioPortal exposes **only z-score** mRNA for this study (`luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores`). Datahub has no `data_mrna_seq_v2_rsem.txt`.

| score | status | why |
|---|---|---|
| ESTIMATE ImmuneScore / StromalScore / TumorPurity | **skipped** | official ssGSEA ranks genes *within a sample* on raw/log expression; ranking across-sample z-scores is a different transform; Affymetrix cosine purity is invalid on z-scores |
| xCell / MCP-counter | **skipped** | need counts or TPM |
| ImmuneScore used here | A1 8-gene mean z | same gene list as PR 139; Spearman is rank-based |
| published immune column | IMSIG T cells | deposited on `data_clinical_sample.txt`; not recomputed from genes |

Do not quote an ESTIMATE ImmuneScore for OncoSG from these public files.

**Covariate:** published sample-level `PURITY` (OncoSG clinical attribute). This is **not** TCGA ABSOLUTE and **not** ESTIMATE cosine purity.

**Statistic:** Spearman on deposited z-scores. Partial = Pearson of average-rank residuals after published PURITY; two-sided t, df = n−3; Fisher-z 95% CI with SE = 1/√(n−4). Same estimator as PR 139.

## Main Spearman (n={n})

Partial residualises on published PURITY.

| pair | n | ρ | p | partial ρ \\| PURITY | partial p |
|---|---:|---:|---:|---:|---:|
{row_md("CLDN4", "CD8A")}
{row_md("CLDN4", "ImmuneScore", "CLDN4 vs ImmuneScore (A1 8-gene)")}
{row_md("CLDN4", "CD274")}
{row_md("CLDN4", "IMSIG_T_CELLS", "CLDN4 vs IMSIG T cells (published)")}
{row_md("CLDN4", "TACSTD2")}
{row_md("CLDN4", "PURITY")}
{row_md("CLDN4", "TJ_noCLDN4", "CLDN4 vs TJ (no CLDN4)")}
{row_md("CD8A", "ImmuneScore", "CD8A vs ImmuneScore (positive control)")}
{row_md("CD274", "ImmuneScore", "CD274 vs ImmuneScore")}
{row_md("TACSTD2", "CD8A", "TACSTD2 vs CD8A (companion; PR 139)")}

CLDN4 vs PURITY Spearman ρ = {fmt_rho(g("CLDN4", "PURITY", "unadj_rho"))} (p = {fmt_p(g("CLDN4", "PURITY", "unadj_p"))}). Partialling purity shrinks the immune |ρ| but does not flip the CD8A / ImmuneScore / IMSIG T-cell signs.

CD8A vs ImmuneScore ρ = {fmt_rho(g("CD8A", "ImmuneScore", "unadj_rho"))} is a positive control and is partly by construction (CD8A is one of the eight ImmuneScore genes).

CD274 (PD-L1 transcript) is **not** a CD8 / ImmuneScore substitute. On this matrix CLDN4–CD274 is weaker than CLDN4–CD8A.

## Extra scatter

Unadjusted scatters for the three requested pairs, plus extra residual and companion panels:

- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_cldn4_vs_cd274.png`
- `figures/fig4_spearman_forest.png`
- `figures/fig5_correlation_heatmap.png`
- `figures/fig6_extra_residual_scatter.png` — rank residuals after PURITY (CD8A / ImmuneScore / CD274)
- `figures/fig7_extra_imsig_tacstd2_scatter.png` — CLDN4 vs IMSIG T cells and vs TACSTD2
- `figures/fig8_extra_purity_control_scatter.png` — CLDN4 vs PURITY; CD8A vs ImmuneScore control

## What is not done

- No re-audit of the A9 3/3 catalog or of PR 139 TACSTD2 numbers.
- No ESTIMATE / xCell / MCP on z-scores.
- No ICI ORR / PFS model (labels are not deposited).
- No imputation of the 12 portal-listed RNA IDs that lack a public z-score column.
- This is East-Asian surgical LUAD, not an ICI-response cohort. A negative bulk correlation is not evidence that CLDN4-high tumors fail checkpoint blockade.

## Files

- `analyze.py` — Datahub download, complete-case n, Spearman, extra scatters
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `coverage.tsv`, `gene_coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png` … `fig8_extra_purity_control_scatter.png`

```bash
python3 methods/oncosg_cldn4_immune/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)


def scatter(ax, x, y, xlab, ylab, rho, p, n, color):
    ax.scatter(x, y, s=14, alpha=0.55, c=color, linewidths=0)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() >= 6:
        lr = stats.linregress(np.asarray(x)[m], np.asarray(y)[m])
        xs = np.linspace(np.nanmin(np.asarray(x)[m]), np.nanmax(np.asarray(x)[m]), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#333333", lw=1.2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"ρ={fmt_rho(rho)}  p={fmt_p(p)}  n={n}", fontsize=10)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    expr_path = ensure_file("expression")
    clin_path = ensure_file("clinical_sample")
    pat_path = ensure_file("clinical_patient")

    raw = pd.read_csv(expr_path, sep="\t")
    if "Hugo_Symbol" not in raw.columns:
        raise SystemExit(f"unexpected expression header: {list(raw.columns)[:5]}")
    n_hugo_rows = int(len(raw))
    n_hugo_dups = int(raw["Hugo_Symbol"].duplicated().sum())
    expr = raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr = expr.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
    expr = expr.set_index("Hugo_Symbol")
    samples = list(expr.columns)

    clin_full = pd.read_csv(clin_path, sep="\t", comment="#")
    pat = pd.read_csv(pat_path, sep="\t", comment="#")
    clin = clin_full.set_index("SAMPLE_ID", drop=False)
    missing_clin = set(samples) - set(clin.index)
    if missing_clin:
        raise SystemExit(f"expression samples missing from clinical: {missing_clin}")

    for g in ["CLDN4", "CD8A", "CD274", "TACSTD2"]:
        if g not in expr.index:
            raise SystemExit(f"{g} missing from public z-score matrix")

    cldn4 = expr.loc["CLDN4", samples].astype(float)
    cd8a = expr.loc["CD8A", samples].astype(float)
    cd274 = expr.loc["CD274", samples].astype(float)
    tacstd2 = expr.loc["TACSTD2", samples].astype(float)
    purity = pd.to_numeric(clin.loc[samples, "PURITY"], errors="coerce")
    imsig_t = pd.to_numeric(clin.loc[samples, "IMSIG_T_CELLS"], errors="coerce")

    immune, imm_present, imm_missing = mean_present(expr, IMMUNE8, samples)
    tj, tj_present, tj_missing = mean_present(expr, TJ_NO_CLDN4, samples)

    n = len(samples)
    if n != 169:
        print(f"WARNING: matrix n={n} (expected 169)", flush=True)
    if int(purity.notna().sum()) != n:
        raise SystemExit(f"PURITY incomplete on matrix: {int(purity.notna().sum())}/{n}")
    for name, s in [("CLDN4", cldn4), ("CD8A", cd8a), ("CD274", cd274)]:
        if int(s.notna().sum()) != n:
            raise SystemExit(f"{name} missing values: {int(s.isna().sum())}")

    x = cldn4.to_numpy(float)
    y_cd8 = cd8a.to_numpy(float)
    y_imm = immune.to_numpy(float)
    y_pd = cd274.to_numpy(float)
    y_t2 = tacstd2.to_numpy(float)
    y_ims = imsig_t.to_numpy(float)
    y_tj = tj.to_numpy(float)
    z = purity.to_numpy(float)

    pur_name = "published PURITY"
    rows = [
        rec_row("CLDN4", "CD8A", x, y_cd8, z, pur_name),
        rec_row(
            "CLDN4",
            "ImmuneScore",
            x,
            y_imm,
            z,
            pur_name,
            note="A1 8-gene mean z; not ESTIMATE",
        ),
        rec_row("CLDN4", "CD274", x, y_pd, z, pur_name),
        rec_row(
            "CLDN4",
            "IMSIG_T_CELLS",
            x,
            y_ims,
            z,
            pur_name,
            note="published clinical IMSIG T-cell column",
        ),
        rec_row("CLDN4", "TACSTD2", x, y_t2, z, pur_name, note="companion; PR 139 is the TACSTD2 claim"),
        rec_row("CLDN4", "PURITY", x, z, None, ""),
        rec_row("CLDN4", "TJ_noCLDN4", x, y_tj, z, pur_name),
        rec_row("CD8A", "ImmuneScore", y_cd8, y_imm, None, "", note="positive control"),
        rec_row("CD274", "ImmuneScore", y_pd, y_imm, z, pur_name),
        rec_row("CD274", "CD8A", y_pd, y_cd8, z, pur_name),
        rec_row("TACSTD2", "CD8A", y_t2, y_cd8, z, pur_name, note="companion; should match PR 139 unadj/partial"),
        rec_row("TACSTD2", "ImmuneScore", y_t2, y_imm, z, pur_name, note="companion"),
        rec_row("TACSTD2", "CD274", y_t2, y_pd, z, pur_name, note="companion"),
        rec_row("TACSTD2", "PURITY", y_t2, z, None, ""),
    ]
    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)

    per = pd.DataFrame(
        {
            "sample_id": samples,
            "patient_id": clin.loc[samples, "PATIENT_ID"].to_numpy(),
            "CLDN4_z": x,
            "CD8A_z": y_cd8,
            "CD274_z": y_pd,
            "TACSTD2_z": y_t2,
            "ImmuneScore_A1_8gene": y_imm,
            "IMSIG_T_CELLS": y_ims,
            "TJ_noCLDN4_meanz": y_tj,
            "PURITY": z,
            "ONCOTREE_CODE": clin.loc[samples, "ONCOTREE_CODE"].to_numpy(),
            "CANCER_TYPE_DETAILED": clin.loc[samples, "CANCER_TYPE_DETAILED"].to_numpy(),
            "SAMPLE_CLASS": clin.loc[samples, "SAMPLE_CLASS"].to_numpy(),
        }
    )
    if "PATIENT_ID" in pat.columns:
        pat_i = pat.set_index("PATIENT_ID")
        for col in ["SEX", "SMOKING_STATUS", "STAGE", "ETHNICITY", "COHORT"]:
            if col in pat_i.columns:
                per[col] = per["patient_id"].map(pat_i[col])
    per.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)

    gene_cov = [
        {
            "feature": "CLDN4",
            "n_genes_defined": 1,
            "n_genes_present": 1,
            "genes_present": "CLDN4",
            "genes_missing": "",
            "source": "public z-score matrix",
        },
        {
            "feature": "CD8A",
            "n_genes_defined": 1,
            "n_genes_present": 1,
            "genes_present": "CD8A",
            "genes_missing": "",
            "source": "public z-score matrix",
        },
        {
            "feature": "CD274",
            "n_genes_defined": 1,
            "n_genes_present": 1,
            "genes_present": "CD274",
            "genes_missing": "",
            "source": "public z-score matrix",
        },
        {
            "feature": "ImmuneScore_A1_8gene",
            "n_genes_defined": len(IMMUNE8),
            "n_genes_present": len(imm_present),
            "genes_present": ",".join(imm_present),
            "genes_missing": ",".join(imm_missing),
            "source": "mean of present gene z-scores; not ESTIMATE",
        },
        {
            "feature": "TJ_noCLDN4",
            "n_genes_defined": len(TJ_NO_CLDN4),
            "n_genes_present": len(tj_present),
            "genes_present": ",".join(tj_present),
            "genes_missing": ",".join(tj_missing),
            "source": "mean of present gene z-scores",
        },
        {
            "feature": "ESTIMATE_ImmuneScore",
            "n_genes_defined": 141,
            "n_genes_present": 0,
            "genes_present": "",
            "genes_missing": "entire SI_geneset; no raw RSEM on cBioPortal",
            "source": "SKIP — public deposit is z-scores only",
        },
        {
            "feature": "IMSIG_T_CELLS",
            "n_genes_defined": 0,
            "n_genes_present": 0,
            "genes_present": "",
            "genes_missing": "",
            "source": "data_clinical_sample.txt published IMSIG column",
        },
    ]
    pd.DataFrame(gene_cov).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    cov = {
        "study": STUDY,
        "datahub_sha": DATAHUB_SHA,
        "n_matrix": n,
        "cbioportal_rna_sample_list_n": 181,
        "rna_list_not_in_public_zscore_matrix": ",".join(RNA_LIST_NOT_IN_MATRIX),
        "clinical_samples": int(len(clin_full)),
        "clinical_purity_non_na": int(pd.to_numeric(clin_full["PURITY"], errors="coerce").notna().sum()),
        "clinical_imsig_t_non_na": int(pd.to_numeric(clin_full["IMSIG_T_CELLS"], errors="coerce").notna().sum()),
        "n_purity_on_matrix": int(purity.notna().sum()),
        "n_imsig_t_on_matrix": int(imsig_t.notna().sum()),
        "n_cldn4": int(cldn4.notna().sum()),
        "n_cd8a": int(cd8a.notna().sum()),
        "n_cd274": int(cd274.notna().sum()),
        "genes_in_matrix": int(expr.shape[0]),
        "hugo_rows_before_dedup": n_hugo_rows,
        "hugo_duplicate_rows_dropped": n_hugo_dups,
        "oncotree_luad": int((clin.loc[samples, "ONCOTREE_CODE"] == "LUAD").sum()),
        "immune8_present": ",".join(imm_present),
        "immune8_missing": ",".join(imm_missing),
        "tj_present": ",".join(tj_present),
        "tj_missing": ",".join(tj_missing),
        "ESTIMATE_status": "skipped_no_raw_RSEM",
        "ICI_labels": "not_deposited",
        "expression_unit": "cBioPortal all-sample z-score",
        "purity_field": "published PURITY",
    }
    pd.DataFrame([cov]).to_csv(TABLES / "coverage.tsv", sep="\t", index=False)

    def g(pred, end, col):
        hit = stats_df[(stats_df.predictor == pred) & (stats_df.endpoint == end)]
        return hit.iloc[0][col] if len(hit) else np.nan

    one = pd.DataFrame(
        [
            {
                "dataset": "OncoSG LUAD public RNA (luad_oncosg_2020)",
                "n": n,
                "n_portal_rna_list": 181,
                "ImmuneScore_definition": "A1 8-gene T-cell effector mean z (not ESTIMATE)",
                "CLDN4_vs_CD8A_rho": g("CLDN4", "CD8A", "unadj_rho"),
                "CLDN4_vs_CD8A_p": g("CLDN4", "CD8A", "unadj_p"),
                "CLDN4_vs_ImmuneScore_rho": g("CLDN4", "ImmuneScore", "unadj_rho"),
                "CLDN4_vs_ImmuneScore_p": g("CLDN4", "ImmuneScore", "unadj_p"),
                "CLDN4_vs_CD274_rho": g("CLDN4", "CD274", "unadj_rho"),
                "CLDN4_vs_CD274_p": g("CLDN4", "CD274", "unadj_p"),
                "CLDN4_vs_CD8A_partial_rho": g("CLDN4", "CD8A", "partial_rho"),
                "CLDN4_vs_CD8A_partial_p": g("CLDN4", "CD8A", "partial_p"),
                "CLDN4_vs_ImmuneScore_partial_rho": g("CLDN4", "ImmuneScore", "partial_rho"),
                "CLDN4_vs_ImmuneScore_partial_p": g("CLDN4", "ImmuneScore", "partial_p"),
                "CLDN4_vs_CD274_partial_rho": g("CLDN4", "CD274", "partial_rho"),
                "CLDN4_vs_CD274_partial_p": g("CLDN4", "CD274", "partial_p"),
                "CLDN4_vs_IMSIG_T_rho": g("CLDN4", "IMSIG_T_CELLS", "unadj_rho"),
                "CLDN4_vs_IMSIG_T_p": g("CLDN4", "IMSIG_T_CELLS", "unadj_p"),
                "purity": pur_name,
                "ESTIMATE": "skipped_no_raw_RSEM",
                "ICI_labels": "not_deposited",
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    summary = {
        "question": "OncoSG LUAD public RNA: CLDN4 vs CD8A / ImmuneScore / CD274",
        "additive": True,
        "a9_catalog_taken_as_given": True,
        "honest_n": n,
        "portal_rna_list_n": 181,
        "ImmuneScore": "A1 8-gene mean z; ESTIMATE skipped (z-scores only)",
        "coverage": cov,
        "one_row": one.iloc[0].to_dict(),
        "stats": stats_df.to_dict(orient="records"),
        "cannot_compute": [
            "ESTIMATE ImmuneScore / StromalScore / TumorPurity (no raw RSEM; only z-scores)",
            "xCell (needs counts/TPM + spillover calibration)",
            "MCP-counter as published (marker means on log2 TPM, not z-scores)",
        ],
        "datahub_sha": DATAHUB_SHA,
        "expression_sha256": FILES["expression"][1],
        "clinical_sha256": FILES["clinical_sample"][1],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    write_finding(stats_df, cov)

    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(
        ax,
        x,
        y_cd8,
        "CLDN4 z-score",
        "CD8A z-score",
        g("CLDN4", "CD8A", "unadj_rho"),
        g("CLDN4", "CD8A", "unadj_p"),
        n,
        "#4C78A8",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_cldn4_vs_cd8a.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(
        ax,
        x,
        y_imm,
        "CLDN4 z-score",
        "ImmuneScore (A1 8-gene mean z)",
        g("CLDN4", "ImmuneScore", "unadj_rho"),
        g("CLDN4", "ImmuneScore", "unadj_p"),
        n,
        "#E45756",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_cldn4_vs_immunescore.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    scatter(
        ax,
        x,
        y_pd,
        "CLDN4 z-score",
        "CD274 z-score",
        g("CLDN4", "CD274", "unadj_rho"),
        g("CLDN4", "CD274", "unadj_p"),
        n,
        "#F58518",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_cldn4_vs_cd274.png", bbox_inches="tight")
    plt.close(fig)

    forest = [
        ("CLDN4 vs CD8A", g("CLDN4", "CD8A", "unadj_rho"), "#4C78A8"),
        ("CLDN4 vs ImmuneScore", g("CLDN4", "ImmuneScore", "unadj_rho"), "#E45756"),
        ("CLDN4 vs CD274", g("CLDN4", "CD274", "unadj_rho"), "#F58518"),
        ("CLDN4 vs IMSIG T cells", g("CLDN4", "IMSIG_T_CELLS", "unadj_rho"), "#B279A2"),
        ("CLDN4 vs CD8A | PURITY", g("CLDN4", "CD8A", "partial_rho"), "#9ecae1"),
        ("CLDN4 vs ImmuneScore | PURITY", g("CLDN4", "ImmuneScore", "partial_rho"), "#f4a582"),
        ("CLDN4 vs CD274 | PURITY", g("CLDN4", "CD274", "partial_rho"), "#fdb462"),
        ("CLDN4 vs TACSTD2", g("CLDN4", "TACSTD2", "unadj_rho"), "#54A24B"),
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    y = np.arange(len(forest))
    ax.axvline(0, color="0.5", lw=1)
    for i, (_lab, rho, c) in enumerate(forest):
        ax.plot(rho, i, "o", color=c, ms=7)
        ax.plot([0, rho], [i, i], color=c, lw=2)
    ax.set_yticks(y, [r[0] for r in forest], fontsize=9)
    ax.set_xlabel(f"Spearman ρ (n={n}; partial = rank residual on published PURITY)")
    ax.set_title("OncoSG LUAD public RNA")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_spearman_forest.png", bbox_inches="tight")
    plt.close(fig)

    mat = pd.DataFrame(
        {
            "CLDN4": x,
            "CD8A": y_cd8,
            "ImmuneScore": y_imm,
            "CD274": y_pd,
            "IMSIG_T": y_ims,
            "TACSTD2": y_t2,
            "PURITY": z,
            "TJ_noCLDN4": y_tj,
        }
    )
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
    ax.set_title(f"OncoSG LUAD public RNA n={n}")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig5_correlation_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    # Extra scatter: purity rank residuals for the three requested endpoints.
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.7), sharey=True)
    for ax, yv, title, feat, color in zip(
        axes,
        [y_cd8, y_imm, y_pd],
        ["CD8A", "ImmuneScore (A1 8-gene)", "CD274"],
        ["CD8A", "ImmuneScore", "CD274"],
        ["#4C78A8", "#E45756", "#F58518"],
    ):
        xr = residualize(stats.rankdata(x), stats.rankdata(z).reshape(-1, 1))
        yr = residualize(stats.rankdata(yv), stats.rankdata(z).reshape(-1, 1))
        ax.scatter(xr, yr, s=14, alpha=0.55, c=color, linewidths=0)
        lr = stats.linregress(xr, yr)
        xs = np.linspace(xr.min(), xr.max(), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#333333", lw=1.2)
        ax.set_title(
            f"{title}\npartial ρ={fmt_rho(g('CLDN4', feat, 'partial_rho'))}  "
            f"p={fmt_p(g('CLDN4', feat, 'partial_p'))}"
        )
        ax.set_xlabel("CLDN4 rank residual")
    axes[0].set_ylabel("Feature rank residual")
    fig.suptitle(f"Extra scatter — OncoSG rank residuals after PURITY n={n}", y=1.03, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig6_extra_residual_scatter.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    scatter(
        axes[0],
        x,
        y_ims,
        "CLDN4 z-score",
        "IMSIG T cells (published)",
        g("CLDN4", "IMSIG_T_CELLS", "unadj_rho"),
        g("CLDN4", "IMSIG_T_CELLS", "unadj_p"),
        n,
        "#B279A2",
    )
    scatter(
        axes[1],
        x,
        y_t2,
        "CLDN4 z-score",
        "TACSTD2 z-score",
        g("CLDN4", "TACSTD2", "unadj_rho"),
        g("CLDN4", "TACSTD2", "unadj_p"),
        n,
        "#54A24B",
    )
    fig.suptitle(f"Extra scatter — IMSIG T / TACSTD2 companion n={n}", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig7_extra_imsig_tacstd2_scatter.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    scatter(
        axes[0],
        x,
        z,
        "CLDN4 z-score",
        "published PURITY",
        g("CLDN4", "PURITY", "unadj_rho"),
        g("CLDN4", "PURITY", "unadj_p"),
        n,
        "#72B7B2",
    )
    scatter(
        axes[1],
        y_cd8,
        y_imm,
        "CD8A z-score",
        "ImmuneScore (A1 8-gene mean z)",
        g("CD8A", "ImmuneScore", "unadj_rho"),
        g("CD8A", "ImmuneScore", "unadj_p"),
        n,
        "#4C78A8",
    )
    fig.suptitle(f"Extra scatter — purity and CD8A–ImmuneScore control n={n}", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig8_extra_purity_control_scatter.png", bbox_inches="tight")
    plt.close(fig)

    print(one.T.to_string())
    print("wrote", TABLES)
    print("wrote", FIGURES)
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()
