#!/usr/bin/env python3
"""DepMap 24Q4 lung lines: CLDN4 RNA vs NHEJ and cGAS–STING genes.

Expression is the primary cut (CCLE/DepMap log2(TPM+1)). CRISPR Chronos
co-dependency is reported on the overlapping lung lines. CLDN4 Chronos is
expected to be nearly flat (prior cut: 0/123 lung lines < −0.5); partner
profiles are summarized with SD so a flat correlation is not read as biology.

This script does not restate the already reported CLDN4–Hallmark IFN-γ
Spearman. STAT1 is flagged because it sits inside that IFN-γ gene set.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

NHEJ = ["PRKDC", "XRCC4", "LIG4", "TP53BP1"]
STING = ["STING1", "CGAS", "STAT1"]
SENSOR = ["STING1", "CGAS"]  # STAT1 is downstream IFN, not the sensor
PANEL = NHEJ + STING
AXIS = {"NHEJ": NHEJ, "STING": STING, "SENSOR": SENSOR}
BOOT_N = 5000
BOOT_SEED = 0


def cohort_tag(row: pd.Series) -> str:
    disease = str(row.get("OncotreePrimaryDisease") or "")
    subtype = str(row.get("OncotreeSubtype") or "")
    if disease == "Non-Small Cell Lung Cancer":
        if subtype == "Lung Adenocarcinoma":
            return "LUAD"
        if subtype == "Lung Squamous Cell Carcinoma":
            return "LUSC"
        return "other_NSCLC"
    if disease == "Lung Neuroendocrine Tumor" or subtype == "Small Cell Lung Cancer":
        return "SCLC_NET"
    return "other_lung"


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    z = df[use].apply(zscore, axis=0)
    return z.mean(axis=1), use


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(a))
    if n < 8:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(a["x"], a["y"])
    return {"n": n, "rho": float(rho), "p": float(p)}


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = BOOT_N, seed: int = BOOT_SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def partial_spearman(x: pd.Series, y: pd.Series, covar: pd.Series) -> dict:
    """Spearman partial correlation: rank, residualize on covar, then Pearson."""
    a = pd.concat(
        [x.rename("x"), y.rename("y"), covar.rename("c")], axis=1
    ).dropna()
    n = int(len(a))
    if n < 10:
        return {"n": n, "rho": None, "p": None}
    ranks = a.rank(method="average")
    def resid(col: str) -> np.ndarray:
        slope, intercept, *_ = stats.linregress(ranks["c"], ranks[col])
        return ranks[col].to_numpy() - (intercept + slope * ranks["c"].to_numpy())
    rx, ry = resid("x"), resid("y")
    if np.std(rx) == 0 or np.std(ry) == 0:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.pearsonr(rx, ry)
    return {"n": n, "rho": float(rho), "p": float(p)}


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's δ for x (Q4) versus y (Q1). Positive means Q4 is larger."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Vectorized pairwise compare, chunked if needed. n≈50 is small.
    diff = x[:, None] - y[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def bh_fdr(pvals: list[float | None]) -> list[float | None]:
    idx = [i for i, p in enumerate(pvals) if p is not None and np.isfinite(p)]
    out: list[float | None] = [None] * len(pvals)
    if not idx:
        return out
    p = np.array([pvals[i] for i in idx], dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, j in enumerate(order[::-1], start=1):
        k = n - rank + 1
        val = min(prev, p[j] * n / k)
        q[j] = val
        prev = val
    for j, i in enumerate(idx):
        out[i] = float(min(1.0, q[j]))
    return out


def attach_fdr(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    df["q"] = np.nan
    if df.empty:
        return df
    grouper = group_cols[0] if len(group_cols) == 1 else group_cols
    for _, idx in df.groupby(grouper, dropna=False).groups.items():
        rows = list(idx)
        q = bh_fdr([None if pd.isna(df.at[i, "p"]) else float(df.at[i, "p"]) for i in rows])
        for i, qi in zip(rows, q):
            df.at[i, "q"] = np.nan if qi is None else qi
    return df


def fmt_p(p: float | None) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(rho: float | None) -> str:
    if rho is None or (isinstance(rho, float) and not np.isfinite(rho)):
        return "NA"
    return f"{rho:+.3f}"


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return f"[{lo:+.2f}, {hi:+.2f}]"


def add_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for name, genes in AXIS.items():
        score, used = signature(out, genes)
        out[f"sig_{name}"] = score
        out.attrs[f"{name}_used"] = used
    return out


def corr_block(
    frame: pd.DataFrame,
    cut: str,
    pairs: list[tuple[str, str, str, str]],
    with_ci: bool,
) -> list[dict]:
    """pairs: (family, predictor_label, predictor_col, target_label/col or special)."""
    rows = []
    for family, pred_label, pred_col, target in pairs:
        rec = spearman(frame[pred_col], frame[target])
        rec.update(
            {
                "cut": cut,
                "family": family,
                "predictor": pred_label,
                "target": target,
                "ci95_low": None,
                "ci95_high": None,
            }
        )
        if with_ci and rec["n"] and rec["n"] >= 8 and rec["rho"] is not None:
            a = pd.concat([frame[pred_col], frame[target]], axis=1).dropna()
            # Stable seed per pair so CIs do not depend on row order in the file.
            seed = BOOT_SEED + (sum(ord(c) for c in f"{cut}|{pred_label}|{target}") % 10000)
            lo, hi = bootstrap_spearman_ci(a.iloc[:, 0].to_numpy(), a.iloc[:, 1].to_numpy(), seed=seed)
            rec["ci95_low"] = lo
            rec["ci95_high"] = hi
        rows.append(rec)
    return rows


def expr_pairs() -> list[tuple[str, str, str, str]]:
    pairs = []
    for g in PANEL:
        pairs.append(("CLDN4_vs_gene", "CLDN4_RNA", "CLDN4", g))
    pairs.append(("CLDN4_vs_score", "CLDN4_RNA", "CLDN4", "sig_NHEJ"))
    pairs.append(("CLDN4_vs_score", "CLDN4_RNA", "CLDN4", "sig_STING"))
    pairs.append(("CLDN4_vs_score", "CLDN4_RNA", "CLDN4", "sig_SENSOR"))
    pairs.append(("score_vs_score", "sig_NHEJ", "sig_NHEJ", "sig_STING"))
    pairs.append(("score_vs_score", "sig_NHEJ", "sig_NHEJ", "sig_SENSOR"))
    return pairs


def chronos_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def write_md_table(df: pd.DataFrame, cols: list[str], headers: list[str]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if c in {"rho"}:
                cells.append(fmt_rho(None if pd.isna(v) else float(v)))
            elif c in {"p", "q"}:
                cells.append(fmt_p(None if pd.isna(v) else float(v)))
            elif c == "ci":
                cells.append(v)
            elif c == "n":
                cells.append(str(int(v)) if pd.notna(v) else "NA")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def call_phrase(rho: float, q: float) -> str:
    if rho is None or q is None or not np.isfinite(rho) or not np.isfinite(q):
        return "not estimated"
    if q < 0.05 and rho >= 0.15:
        return "positive at q<0.05"
    if q < 0.05 and rho <= -0.15:
        return "negative at q<0.05"
    if q < 0.05:
        return "q<0.05 but |ρ|<0.15"
    return "not q<0.05"


def finding_text(
    stats: dict,
    corr: pd.DataFrame,
    crispr_corr: pd.DataFrame,
    q4: pd.DataFrame,
    sd: pd.DataFrame,
    prkdc_lung: pd.DataFrame,
) -> str:
    lung = corr[(corr["cut"] == "lung") & (corr["family"] == "CLDN4_vs_gene")].copy()
    nsclc = corr[(corr["cut"] == "NSCLC") & (corr["family"] == "CLDN4_vs_gene")].copy()

    def row_line(frame: pd.DataFrame) -> str:
        bits = []
        for _, r in frame.iterrows():
            bits.append(
                f"{r['target']} ρ={fmt_rho(r['rho'])} (p={fmt_p(r['p'])}, q={fmt_p(r['q'])}, {fmt_ci(r['ci95_low'], r['ci95_high'])})"
            )
        return "; ".join(bits)

    def score_line(cut: str, target: str) -> str:
        r = corr[(corr["cut"] == cut) & (corr["target"] == target) & (corr["predictor"] == "CLDN4_RNA")].iloc[0]
        return f"ρ={fmt_rho(r['rho'])}, p={fmt_p(r['p'])}, q={fmt_p(r['q'])}, 95% CI {fmt_ci(r['ci95_low'], r['ci95_high'])}, n={int(r['n'])}"

    def pick(cut: str, target: str, predictor: str = "CLDN4_RNA") -> pd.Series:
        hit = corr[(corr["cut"] == cut) & (corr["target"] == target) & (corr["predictor"] == predictor)]
        return hit.iloc[0]

    def rho_q(cut: str, target: str, predictor: str = "CLDN4_RNA") -> str:
        r = pick(cut, target, predictor)
        return f"{fmt_rho(r['rho'])} ({fmt_p(r['q'])})"

    table_cuts = ["lung", "NSCLC", "LUAD", "LUSC", "SCLC_NET", "other_NSCLC"]
    n_by_cut = {c: int(pick(c, "STING1")["n"]) for c in table_cuts}
    header = "| gene | " + " | ".join(f"{c} n={n_by_cut[c]}" for c in table_cuts) + " |"
    sep = "|" + "|".join(["---"] * (1 + len(table_cuts))) + "|"
    body_rows = []
    for gene in PANEL + ["sig_NHEJ", "sig_STING", "sig_SENSOR"]:
        label = {
            "sig_NHEJ": "NHEJ mean-z",
            "sig_STING": "STING mean-z",
            "sig_SENSOR": "SENSOR mean-z",
        }.get(gene, gene)
        body_rows.append("| " + label + " | " + " | ".join(rho_q(c, gene) for c in table_cuts) + " |")
    gene_table = "\n".join([header, sep, *body_rows])

    nhej_sting = corr[(corr["cut"] == "lung") & (corr["target"] == "sig_STING") & (corr["predictor"] == "sig_NHEJ")].iloc[0]
    nhej_sensor = corr[(corr["cut"] == "lung") & (corr["target"] == "sig_SENSOR") & (corr["predictor"] == "sig_NHEJ")].iloc[0]

    def calls(frame: pd.DataFrame) -> str:
        bits = []
        for _, r in frame.iterrows():
            bits.append(f"{r['target']} {call_phrase(float(r['rho']), float(r['q']))} (ρ={fmt_rho(r['rho'])}, q={fmt_p(r['q'])})")
        return "; ".join(bits)

    cl = sd[sd["gene"] == "CLDN4"].iloc[0]
    crispr_genes = ", ".join(sd["gene"].tolist())
    n_prkdc = int(len(prkdc_lung))
    if n_prkdc:
        prkdc_lines = "; ".join(
            f"{r.CellLineName} ({r.group}, Chronos {r.PRKDC_screen:+.3f})"
            for r in prkdc_lung.itertuples()
        )
    else:
        prkdc_lines = "none"

    def crispr_bits(predictor: str) -> str:
        sub = crispr_corr[(crispr_corr["cut"] == "lung_CRISPR") & (crispr_corr["predictor"] == predictor)]
        bits = []
        for _, r in sub.iterrows():
            bits.append(f"{r['target']} ρ={fmt_rho(r['rho'])} (n={int(r['n'])}, p={fmt_p(r['p'])}, q={fmt_p(r['q'])})")
        return "; ".join(bits)

    q4_bits = []
    for _, r in q4[q4["cut"] == "lung"].iterrows():
        q4_bits.append(
            f"{r['target']}: Q4 median {r['q4_median']:+.3f} vs Q1 {r['q1_median']:+.3f} "
            f"(Cliff δ={r['cliffs_delta']:+.2f}, p={fmt_p(r['p'])}, n={int(r['n_q4'])}/{int(r['n_q1'])})"
        )

    return f"""# DepMap 24Q4 lung: CLDN4 expression vs NHEJ and cGAS–STING

**Additive public evidence.** Cultured lung lines only. This page is CLDN4 RNA against four NHEJ genes (`PRKDC`, `XRCC4`, `LIG4`, `TP53BP1`) and three cGAS–STING-axis genes (`STING1`, `CGAS`, `STAT1`), plus Chronos co-dependency on the CRISPR overlap. It does not restate the already reported all-lung CLDN4–Hallmark IFN-γ Spearman.

**Data.** DepMap Public 24Q4 ([10.6084/m9.figshare.27993248](https://doi.org/10.6084/m9.figshare.27993248)). RNA = `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, log2(TPM+1), the CCLE/DepMap expression matrix. Dependency = `CRISPRGeneEffect.csv` Chronos (more negative = stronger dependency). Expression columns use the HGNC symbols above (`STING1`, not `TMEM173`; `CGAS`, not `MB21D1`). Integrated Chronos is missing **PRKDC** (no `PRKDC` / `XRCC7` / Entrez 5591 column). PRKDC remains only in `ScreenGeneEffect.csv`, and only on Humagne-CD screens.

**Cohort.** `OncotreeLineage == Lung` and `ModelType == Cell Line`: **{stats['n_lung_models']}** models in `Model.csv` → **{stats['n_lung_rna']}** with finite RNA for CLDN4 and the seven partners. NSCLC **{stats['n_nsclc']}** (LUAD {stats['n_luad']} / LUSC {stats['n_lusc']} / other NSCLC {stats['n_other_nsclc']}). SCLC+NET **{stats['n_sclc']}**. Other lung **{stats['n_other']}**. CRISPR overlap with finite integrated Chronos for CLDN4 and every partner that exists in `CRISPRGeneEffect` ({crispr_genes}): **{stats['n_crispr_lung']}** lung lines (**{stats['n_crispr_nsclc']}** NSCLC). Humagne-CD screens with a finite PRKDC gene effect: **{n_prkdc}** lung lines. No immune infiltrate and no IFN treatment in the dish.

Scores are the mean of gene-wise z-scores computed **inside each cut**. NHEJ = PRKDC, XRCC4, LIG4, TP53BP1. STING = STING1, CGAS, STAT1. SENSOR = STING1 and CGAS only, because STAT1 is an interferon-stimulated gene and a Hallmark IFN-γ member. Spearman is two-sided, complete cases. BH *q* is within the family on that cut (seven genes; three scores; Chronos pairs). Bootstrap 95% CI: {BOOT_N} resamples. A gene is called only when q<0.05 and |ρ|≥0.15.

## CLDN4 RNA vs each gene

Each cell is Spearman ρ (BH *q* within that cut’s seven-gene family, or within that cut’s three scores). Master file, with p and bootstrap CIs: `tables/expr_correlations.tsv`.

{gene_table}

Bootstrap 95% CIs on the two large cuts: lung PRKDC {fmt_ci(pick('lung','PRKDC')['ci95_low'], pick('lung','PRKDC')['ci95_high'])}, lung STING1 {fmt_ci(pick('lung','STING1')['ci95_low'], pick('lung','STING1')['ci95_high'])}, NSCLC PRKDC {fmt_ci(pick('NSCLC','PRKDC')['ci95_low'], pick('NSCLC','PRKDC')['ci95_high'])}, NSCLC STING1 {fmt_ci(pick('NSCLC','STING1')['ci95_low'], pick('NSCLC','STING1')['ci95_high'])}. Score CIs: lung NHEJ {score_line('lung', 'sig_NHEJ')}; lung STING {score_line('lung', 'sig_STING')}; lung SENSOR {score_line('lung', 'sig_SENSOR')}; NSCLC NHEJ {score_line('NSCLC', 'sig_NHEJ')}; NSCLC STING {score_line('NSCLC', 'sig_STING')}; NSCLC SENSOR {score_line('NSCLC', 'sig_SENSOR')}.

NHEJ mean-z vs STING mean-z on all lung lines: ρ={fmt_rho(nhej_sting['rho'])}, p={fmt_p(nhej_sting['p'])}, q={fmt_p(nhej_sting['q'])}, 95% CI {fmt_ci(nhej_sting['ci95_low'], nhej_sting['ci95_high'])}, n={int(nhej_sting['n'])}. NHEJ vs SENSOR: ρ={fmt_rho(nhej_sensor['rho'])}, p={fmt_p(nhej_sensor['p'])}, q={fmt_p(nhej_sensor['q'])}, 95% CI {fmt_ci(nhej_sensor['ci95_low'], nhej_sensor['ci95_high'])}. On NSCLC the same NHEJ-vs-STING score correlation is ρ={fmt_rho(pick('NSCLC','sig_STING','sig_NHEJ')['rho'])}, q={fmt_p(pick('NSCLC','sig_STING','sig_NHEJ')['q'])}.

Partial Spearman (rank-residual) on all lung lines: CLDN4 vs STING score given NHEJ score ρ={fmt_rho(stats['partial_cldn4_sting_given_nhej']['rho'])} (p={fmt_p(stats['partial_cldn4_sting_given_nhej']['p'])}, n={stats['partial_cldn4_sting_given_nhej']['n']}); CLDN4 vs NHEJ score given STING score ρ={fmt_rho(stats['partial_cldn4_nhej_given_sting']['rho'])} (p={fmt_p(stats['partial_cldn4_nhej_given_sting']['p'])}, n={stats['partial_cldn4_nhej_given_sting']['n']}); CLDN4 vs SENSOR score given NHEJ score ρ={fmt_rho(stats['partial_cldn4_sensor_given_nhej']['rho'])} (p={fmt_p(stats['partial_cldn4_sensor_given_nhej']['p'])}, n={stats['partial_cldn4_sensor_given_nhej']['n']}).

CLDN4 RNA quartile split (top vs bottom quartile) on the lung RNA cohort: {'; '.join(q4_bits)}.

`sig_*` columns in `tables/lung_lines.tsv` are z-scored on the all-lung cohort. Correlations in the table above recompute z inside each cut.

## What the expression cut shows

Calls use q<0.05 and |ρ|≥0.15 inside the family for that cut.

The gene that stays positive with CLDN4 RNA in the mixed lung set, in NSCLC, and inside LUAD is **STING1** (lung ρ={fmt_rho(pick('lung','STING1')['rho'])}, LUAD ρ={fmt_rho(pick('LUAD','STING1')['rho'])}, q={fmt_p(pick('LUAD','STING1')['q'])}, n={int(pick('LUAD','STING1')['n'])}). SCLC+NET is the same direction (ρ={fmt_rho(pick('SCLC_NET','STING1')['rho'])}, q={fmt_p(pick('SCLC_NET','STING1')['q'])}, n={int(pick('SCLC_NET','STING1')['n'])}). **CGAS** is not called in any of those cuts. **STAT1** is not called in the mixed lung set or in NSCLC. STAT1 is one Hallmark IFN-γ gene; this page does not recompute that score, and a null STAT1 rank is not a revision of it.

The SENSOR mean-z (STING1+CGAS) is called on the mixed lung set and is not called in NSCLC or LUAD. CGAS dilutes STING1. The three-gene STING score is the same story: called on all-lung, not called in NSCLC.

**PRKDC** RNA is negatively associated with CLDN4 on the mixed lung set (ρ={fmt_rho(pick('lung','PRKDC')['rho'])}, q={fmt_p(pick('lung','PRKDC')['q'])}) and on mixed NSCLC (ρ={fmt_rho(pick('NSCLC','PRKDC')['rho'])}, q={fmt_p(pick('NSCLC','PRKDC')['q'])}). It is not a LUAD association (ρ={fmt_rho(pick('LUAD','PRKDC')['rho'])}, q={fmt_p(pick('LUAD','PRKDC')['q'])}, n={int(pick('LUAD','PRKDC')['n'])}). The NSCLC number is carried by the non-LUAD slices (other NSCLC ρ={fmt_rho(pick('other_NSCLC','PRKDC')['rho'])}, q={fmt_p(pick('other_NSCLC','PRKDC')['q'])}, n={int(pick('other_NSCLC','PRKDC')['n'])}; LUSC ρ={fmt_rho(pick('LUSC','PRKDC')['rho'])}, q={fmt_p(pick('LUSC','PRKDC')['q'])}, n={int(pick('LUSC','PRKDC')['n'])}). **XRCC4** and **LIG4** are not called. **TP53BP1** meets the call rule only on the mixed lung set and not in NSCLC. The NHEJ mean-z is therefore not a four-gene block. It is called on mixed lung and mixed NSCLC, not on LUAD (ρ={fmt_rho(pick('LUAD','sig_NHEJ')['rho'])}, q={fmt_p(pick('LUAD','sig_NHEJ')['q'])}). LUSC NHEJ mean-z is ρ={fmt_rho(pick('LUSC','sig_NHEJ')['rho'])}, q={fmt_p(pick('LUSC','sig_NHEJ')['q'])}, n={int(pick('LUSC','sig_NHEJ')['n'])}.

NHEJ mean-z and STING mean-z are uncorrelated on all lung lines (ρ={fmt_rho(nhej_sting['rho'])}). On NSCLC they are weakly positive (ρ={fmt_rho(pick('NSCLC','sig_STING','sig_NHEJ')['rho'])}, q={fmt_p(pick('NSCLC','sig_STING','sig_NHEJ')['q'])}). That is not an NHEJ-up / STING-down pair.

Lung gene calls: {calls(lung)}.

NSCLC gene calls: {calls(nsclc)}.

## CRISPR co-dependency

Chronos on the lung CRISPR∩panel set (n={stats['n_crispr_lung']}). Profile SD is the standard deviation of the gene effect across those lines.

| gene | n | median Chronos | SD | fraction < −0.5 |
|---|---:|---:|---:|---:|
{chr(10).join(
    f"| {r.gene} | {int(r.n)} | {r.median:+.3f} | {r.sd:.3f} | {r.frac_lt_m0_5:.3f} |"
    for r in sd.itertuples()
)}

CLDN4 Chronos median {cl['median']:+.3f}, SD {cl['sd']:.3f}, fraction < −0.5 = {cl['frac_lt_m0_5']:.3f}. That profile is too tight for a co-dependency scan to mean selective CLDN4 essentiality. The CLDN4-Chronos correlations are still reported, as the available co-dependency numbers:

- CLDN4 Chronos vs partner Chronos, lung: {crispr_bits('CLDN4_Chronos')}
- CLDN4 RNA vs partner Chronos, lung (expression versus dependency, not co-dependency): {crispr_bits('CLDN4_RNA')}

NHEJ×STING Chronos pairs use the integrated genes that exist (XRCC4, LIG4, TP53BP1 × STING1, CGAS, STAT1; BH within that 3×3 family). On the lung CRISPR cut the minimum *q* in that 3×3 is {fmt_p(float(crispr_corr[(crispr_corr['cut']=='lung_CRISPR') & (crispr_corr['family']=='NHEJ_Chronos_vs_STING_Chronos')]['q'].min()))}. Table: `tables/crispr_nhej_vs_sting.tsv`. Read each ρ next to the Chronos SD above. A flat profile cannot support a selective co-dependency.

**PRKDC Chronos is not in the integrated matrix.** Humagne-CD `ScreenGeneEffect` has a finite PRKDC value for {n_prkdc} lung cell lines: {prkdc_lines}. n={n_prkdc} is too small for a Spearman. Those screen-level values are listed in `tables/prkdc_humagne_lung.tsv` and are not mixed into the integrated co-dependency table. Guide maps: Avana and KY have no PRKDC guides; Humagne-CD has two guides marked `UsedByChronos=True`.

## What this is not

- Not a tumor immune-exclusion result. These are cell lines.
- Not an IFN-stimulated or ICI result. STAT1 here is basal RNA.
- Not a restatement of the all-lung CLDN4–Hallmark IFN-γ Spearman. STAT1 is one gene inside that axis.
- Not evidence that knocking out CLDN4 changes NHEJ or STING. CLDN4 Chronos in these lines has little spread.
- Not a DNA-PK–versus–cGAS mechanism result. PRKDC dependency is almost absent from this release, and the expression ranks are reported as ranks.

Scatter and forest: `figures/fig_cldn4_nhej_sting.png`.

## How to rerun

```bash
python3 -m pip install -r methods/depmap_nhej_sting/requirements.txt
python3 methods/depmap_nhej_sting/download.py --outdir data/depmap_nhej_sting
python3 methods/depmap_nhej_sting/analyze.py --data data/depmap_nhej_sting --outdir methods/depmap_nhej_sting
```
"""


def make_figure(
    lung: pd.DataFrame,
    corr: pd.DataFrame,
    chronos_mat: pd.DataFrame,
    path: Path,
) -> None:
    fig = plt.figure(figsize=(13.6, 8.8), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    genes = list(reversed(PANEL))
    for offset, cut, color in ((-0.15, "lung", "#1f4e79"), (0.15, "NSCLC", "#c47b2b")):
        sub = corr[(corr["cut"] == cut) & (corr["family"] == "CLDN4_vs_gene")].set_index("target")
        ys = np.arange(len(genes)) + offset
        rhos = [sub.loc[g, "rho"] for g in genes]
        lo = [sub.loc[g, "rho"] - sub.loc[g, "ci95_low"] for g in genes]
        hi = [sub.loc[g, "ci95_high"] - sub.loc[g, "rho"] for g in genes]
        ax.errorbar(rhos, ys, xerr=[lo, hi], fmt="o", color=color, label=f"{cut} n={int(sub.loc[genes[0], 'n'])}", ms=5, lw=1)
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes)
    ax.set_xlabel("Spearman ρ, CLDN4 RNA (95% bootstrap CI)")
    ax.set_title("CLDN4 expression vs panel")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(-0.55, 0.7)

    def scatter(ax, frame, ycol, title, ylabel):
        colors = {
            "LUAD": "#1f4e79",
            "LUSC": "#b85c38",
            "other_NSCLC": "#7aa0c4",
            "SCLC_NET": "#2f6b4f",
            "other_lung": "#888888",
        }
        for name, sub in frame.groupby("group"):
            ax.scatter(sub["CLDN4"], sub[ycol], s=16, alpha=0.8, c=colors.get(name, "#444"), label=f"{name} ({len(sub)})", linewidths=0)
        ax.set_xlabel("CLDN4 log2(TPM+1)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(
            frameon=False,
            fontsize=6.5,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
        )

    ax = fig.add_subplot(gs[0, 1])
    scatter(ax, lung, "PRKDC", "Lung lines, PRKDC RNA", "PRKDC log2(TPM+1)")

    ax = fig.add_subplot(gs[1, 0])
    scatter(ax, lung, "STING1", "Lung lines, STING1 RNA", "STING1 log2(TPM+1)")

    ax = fig.add_subplot(gs[1, 1])
    labels = list(chronos_mat.index)
    mat = chronos_mat.loc[labels, labels].to_numpy(dtype=float)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="equal")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(len(labels)):
        for j in range(len(labels)):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=6, color="#111" if abs(val) < 0.35 else "white")
    ax.set_title(
        f"Integrated Chronos Spearman, lung n={int(chronos_mat.attrs.get('n', 0))}\n(PRKDC absent from CRISPRGeneEffect)"
    )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="ρ")

    fig.suptitle("DepMap 24Q4 lung cell lines: CLDN4, NHEJ, cGAS–STING", fontsize=12)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_nhej_sting")
    ap.add_argument("--outdir", default="methods/depmap_nhej_sting")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(data / "Model.csv")
    expr = chronos_numeric(pd.read_csv(data / "expression_panel.csv"), ["CLDN4"] + PANEL)
    crispr = chronos_numeric(pd.read_csv(data / "crispr_panel.csv"), ["CLDN4"] + PANEL)
    manifest = json.loads((data / "download_manifest.json").read_text())

    meta_cols = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "ModelType",
    ]
    meta_cols = [c for c in meta_cols if c in model.columns]
    df = expr.merge(model[meta_cols], on="ModelID", how="left")
    lung_models = model[(model["OncotreeLineage"] == "Lung") & (model["ModelType"] == "Cell Line")]
    lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"] == "Cell Line")].copy()
    need = ["CLDN4"] + PANEL
    n_lung_rna_any = int(len(lung))
    lung = lung.dropna(subset=need).copy()
    lung["group"] = lung.apply(cohort_tag, axis=1)
    lung = add_scores(lung)

    cohorts = {
        "lung": lung,
        "NSCLC": lung[lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy(),
        "LUAD": lung[lung["group"] == "LUAD"].copy(),
        "LUSC": lung[lung["group"] == "LUSC"].copy(),
        "SCLC_NET": lung[lung["group"] == "SCLC_NET"].copy(),
        "other_NSCLC": lung[lung["group"] == "other_NSCLC"].copy(),
    }
    # Recompute scores inside each smaller cut so z is not inherited from all-lung.
    for name, frame in list(cohorts.items()):
        if name == "lung":
            continue
        cohorts[name] = add_scores(frame)

    corr_rows = []
    for name, frame in cohorts.items():
        with_ci = name in {"lung", "NSCLC"}
        corr_rows.extend(corr_block(frame, name, expr_pairs(), with_ci=with_ci))
    corr = pd.DataFrame(corr_rows)
    corr = attach_fdr(corr, ["cut", "family"])

    # Quartile contrast on lung and NSCLC, using within-cut scores already stored.
    q_rows = []
    for name in ("lung", "NSCLC"):
        frame = cohorts[name]
        lo_thr = frame["CLDN4"].quantile(0.25)
        hi_thr = frame["CLDN4"].quantile(0.75)
        q1 = frame[frame["CLDN4"] <= lo_thr]
        q4 = frame[frame["CLDN4"] >= hi_thr]
        for target in ["sig_NHEJ", "sig_STING", "sig_SENSOR"] + PANEL:
            x = q4[target].dropna().to_numpy()
            y = q1[target].dropna().to_numpy()
            if len(x) < 5 or len(y) < 5:
                p = None
                delta = None
            else:
                p = float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
                delta = cliffs_delta(x, y)
            q_rows.append(
                {
                    "cut": name,
                    "target": target,
                    "n_q4": int(len(x)),
                    "n_q1": int(len(y)),
                    "cldn4_q4_min": float(hi_thr),
                    "cldn4_q1_max": float(lo_thr),
                    "q4_median": float(np.median(x)) if len(x) else None,
                    "q1_median": float(np.median(y)) if len(y) else None,
                    "cliffs_delta": delta,
                    "p": p,
                }
            )
    q4 = pd.DataFrame(q_rows)
    q4 = attach_fdr(q4, ["cut"])

    partial = {
        "partial_cldn4_sting_given_nhej": partial_spearman(lung["CLDN4"], lung["sig_STING"], lung["sig_NHEJ"]),
        "partial_cldn4_nhej_given_sting": partial_spearman(lung["CLDN4"], lung["sig_NHEJ"], lung["sig_STING"]),
        "partial_cldn4_sensor_given_nhej": partial_spearman(lung["CLDN4"], lung["sig_SENSOR"], lung["sig_NHEJ"]),
    }

    # Integrated Chronos. PRKDC is not in CRISPRGeneEffect; do not require it.
    crispr_genes = [g for g in ["CLDN4"] + PANEL if g in crispr.columns]
    partner_genes = [g for g in PANEL if g in crispr.columns]
    nhej_crispr = [g for g in NHEJ if g in crispr.columns]
    sting_crispr = [g for g in STING if g in crispr.columns]
    chron = crispr.rename(columns={g: f"{g}_Chronos" for g in crispr_genes})
    merged = lung.merge(chron, on="ModelID", how="left")
    chron_cols = [f"{g}_Chronos" for g in crispr_genes]
    crispr_lung = merged.dropna(subset=chron_cols).copy()
    crispr_nsclc = crispr_lung[crispr_lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy()

    sd_rows = []
    for g in crispr_genes:
        s = crispr_lung[f"{g}_Chronos"]
        sd_rows.append(
            {
                "gene": g,
                "n": int(s.notna().sum()),
                "median": float(s.median()),
                "sd": float(s.std(ddof=0)),
                "iqr_low": float(s.quantile(0.25)),
                "iqr_high": float(s.quantile(0.75)),
                "frac_lt_m0_5": float((s < -0.5).mean()),
                "min": float(s.min()),
                "max": float(s.max()),
            }
        )
    sd = pd.DataFrame(sd_rows)

    crispr_rows = []
    for cut_name, frame in (("lung_CRISPR", crispr_lung), ("NSCLC_CRISPR", crispr_nsclc)):
        pairs = []
        for g in partner_genes:
            pairs.append(("CLDN4_Chronos_vs_Chronos", "CLDN4_Chronos", "CLDN4_Chronos", f"{g}_Chronos"))
            pairs.append(("CLDN4_RNA_vs_Chronos", "CLDN4_RNA", "CLDN4", f"{g}_Chronos"))
        raw = corr_block(frame, cut_name, pairs, with_ci=(cut_name == "lung_CRISPR"))
        for rec, g in zip(raw[0::2], partner_genes):
            rec["target"] = g
        for rec, g in zip(raw[1::2], partner_genes):
            rec["target"] = g
        crispr_rows.extend(raw)
        cross = []
        for a in nhej_crispr:
            for b in sting_crispr:
                rec = spearman(frame[f"{a}_Chronos"], frame[f"{b}_Chronos"])
                rec.update(
                    {
                        "cut": cut_name,
                        "family": "NHEJ_Chronos_vs_STING_Chronos",
                        "predictor": a,
                        "target": b,
                        "ci95_low": None,
                        "ci95_high": None,
                    }
                )
                cross.append(rec)
        crispr_rows.extend(cross)
    crispr_corr = pd.DataFrame(crispr_rows)
    crispr_corr = attach_fdr(crispr_corr, ["cut", "family"])

    mat_src = pd.DataFrame({g: crispr_lung[f"{g}_Chronos"] for g in crispr_genes})
    rho_mat = mat_src.corr(method="spearman")
    rho_mat.attrs["n"] = int(len(crispr_lung))

    screen = pd.read_csv(data / "screen_prkdc.csv")
    smap = pd.read_csv(data / "CRISPRScreenMap.csv")
    if "PRKDC" not in screen.columns and "PRKDC_screen" in screen.columns:
        screen = screen.rename(columns={"PRKDC_screen": "PRKDC"})
    screen["PRKDC_screen"] = pd.to_numeric(screen["PRKDC"], errors="coerce")
    prkdc_joined = screen.dropna(subset=["PRKDC_screen"]).merge(smap, on="ScreenID", how="left")
    prkdc_joined = prkdc_joined.merge(model[meta_cols], on="ModelID", how="left")
    prkdc_lung = prkdc_joined[
        (prkdc_joined["OncotreeLineage"] == "Lung") & (prkdc_joined["ModelType"] == "Cell Line")
    ].copy()
    prkdc_lung["group"] = prkdc_lung.apply(cohort_tag, axis=1)
    prkdc_keep = [
        c
        for c in [
            "ScreenID",
            "ModelID",
            "CellLineName",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "group",
            "PRKDC_screen",
        ]
        if c in prkdc_lung.columns
    ]
    prkdc_lung[prkdc_keep].to_csv(tabdir / "prkdc_humagne_lung.tsv", sep="\t", index=False)

    # Line table.
    line_cols = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "group",
        "CLDN4",
        *PANEL,
        "sig_NHEJ",
        "sig_STING",
        "sig_SENSOR",
    ]
    line_cols = [c for c in line_cols if c in lung.columns]
    lines = lung[line_cols].merge(
        crispr_lung[["ModelID"] + chron_cols],
        on="ModelID",
        how="left",
    )
    lines.to_csv(tabdir / "lung_lines.tsv", sep="\t", index=False)
    corr.to_csv(tabdir / "expr_correlations.tsv", sep="\t", index=False)
    q4.to_csv(tabdir / "cldn4_quartile.tsv", sep="\t", index=False)
    crispr_corr.to_csv(tabdir / "crispr_correlations.tsv", sep="\t", index=False)
    sd.to_csv(tabdir / "chronos_spread.tsv", sep="\t", index=False)
    cross = crispr_corr[crispr_corr["family"] == "NHEJ_Chronos_vs_STING_Chronos"].copy()
    cross.to_csv(tabdir / "crispr_nhej_vs_sting.tsv", sep="\t", index=False)

    counts = pd.DataFrame(
        [
            {"cohort": "lung_models", "n": int(len(lung_models))},
            {"cohort": "lung_rna_any_model_match", "n": n_lung_rna_any},
            {"cohort": "lung_rna_complete_panel", "n": int(len(lung))},
            {"cohort": "NSCLC", "n": int(len(cohorts["NSCLC"]))},
            {"cohort": "LUAD", "n": int(len(cohorts["LUAD"]))},
            {"cohort": "LUSC", "n": int(len(cohorts["LUSC"]))},
            {"cohort": "SCLC_NET", "n": int(len(cohorts["SCLC_NET"]))},
            {"cohort": "other_lung", "n": int((lung["group"] == "other_lung").sum())},
            {"cohort": "other_NSCLC", "n": int((lung["group"] == "other_NSCLC").sum())},
            {"cohort": "lung_CRISPR_complete_panel", "n": int(len(crispr_lung))},
            {"cohort": "NSCLC_CRISPR_complete_panel", "n": int(len(crispr_nsclc))},
            {"cohort": "lung_PRKDC_Humagne_screen", "n": int(len(prkdc_lung))},
        ]
    )
    counts.to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    key = {
        "release": manifest.get("release"),
        "doi": manifest.get("doi"),
        "n_lung_models": int(len(lung_models)),
        "n_lung_rna_any": n_lung_rna_any,
        "n_lung_rna": int(len(lung)),
        "n_nsclc": int(len(cohorts["NSCLC"])),
        "n_luad": int(len(cohorts["LUAD"])),
        "n_lusc": int(len(cohorts["LUSC"])),
        "n_sclc": int(len(cohorts["SCLC_NET"])),
        "n_other": int((lung["group"] == "other_lung").sum()),
        "n_other_nsclc": int((lung["group"] == "other_NSCLC").sum()),
        "n_crispr_lung": int(len(crispr_lung)),
        "n_crispr_nsclc": int(len(crispr_nsclc)),
        "expression_missing": manifest.get("expression_extract", {}).get("missing_genes"),
        "crispr_missing": manifest.get("crispr_extract", {}).get("missing_genes"),
        "n_prkdc_humagne_lung": int(len(prkdc_lung)),
        "crispr_genes": crispr_genes,
        **partial,
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    make_figure(lung, corr, rho_mat, figdir / "fig_cldn4_nhej_sting.png")
    text = finding_text(key, corr, crispr_corr, q4, sd, prkdc_lung)
    (outdir / "FINDING.md").write_text(text)
    print(text)
    print("--- counts ---")
    print(counts.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
