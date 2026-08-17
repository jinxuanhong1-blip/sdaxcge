#!/usr/bin/env python3
"""ADDITIVE extra — prerank GSEA on GSE126044 CLDN4 Q4 vs Q1.

Public Cho 2020 anti-PD-1 NSCLC counts (n=16). Single-gene CLDN4 vs
ESTIMATE ImmuneScore is already known (Spearman r=-0.524, p=0.037, n=16)
and is recorded, not re-tested.

Engine is the same as methods/cldn4_ko_gsea and scrna_pseudobulk_gsea_meta
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Hallmark IFN-α is secondary only (not in the BH-FDR headline).

Positive NES = enriched at the CLDN4-high (Q4) end of the rank.
Honest n is the actual Q4 / Q1 membership after the quartile cut.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from gsea_core import (  # noqa: E402
    NPERM,
    SEED,
    bh_fdr,
    gsea_prerank,
    quartile_split,
    rank_high_vs_low,
    rank_spearman_vs_target,
)

DATA = HERE / "data"
FIG = HERE / "figures"
TAB = HERE / "tables"
SET_JSON = HERE / "gene_sets.json"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

COUNTS_FILE = "GSE126044_counts.txt.gz"
CLIN_FILE = "GSE126044_clinical.csv"
CONTRAST = "GSE126044_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE126044"

# Known from methods/bulk_immune (not re-tested)
KNOWN_IMMUNESCORE_R = -0.5235294117647059
KNOWN_IMMUNESCORE_P = 0.037412648130149384
KNOWN_IMMUNESCORE_N = 16

HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "KEGG_TIGHT_JUNCTION",
]
SECONDARY = ["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
ALL_SETS = HEADLINE + SECONDARY

HEADLINE_LABEL = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Hallmark IFN-α",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
    "KEGG_TIGHT_JUNCTION": "KEGG tight junction",
}

EXCEL_DATE_SYMBOLS = {"1-Mar", "2-Mar", "1-Sep", "2-Sep"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (p != p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def fmt_fc(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def load_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in ALL_SETS}


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def load_counts() -> pd.DataFrame:
    raw = pd.read_csv(DATA / COUNTS_FILE, sep="\t", index_col=0)
    raw.index = raw.index.astype(str).str.strip()
    raw = raw.loc[~raw.index.isin(EXCEL_DATE_SYMBOLS)]
    raw = raw.loc[raw.index.notna() & (raw.index != "") & (raw.index != "nan")]
    raw = raw.groupby(raw.index).sum()
    return raw.astype(float)


def load_clinical() -> pd.DataFrame:
    clin = pd.read_csv(DATA / CLIN_FILE)
    clin["sample_id"] = clin["sample_id"].astype(str)
    return clin.set_index("sample_id")


def rank_log2fc(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.Series:
    mh = expr.loc[:, high].mean(axis=1)
    ml = expr.loc[:, low].mean(axis=1)
    s = (mh - ml).replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def run_rank(
    rank: pd.Series,
    sets: dict[str, list[str]],
    rank_metric: str,
    n_q4: int,
    n_q1: int,
    n_used: int,
) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = CONTRAST
    g["accession"] = ACCESSION
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_q4"] = n_q4
    g["n_q1"] = n_q1
    g["n_used_for_rank"] = n_used
    g["n_cohort"] = 16
    g["headline"] = g["term"].isin(HEADLINE)
    head = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    g["direction"] = np.where(g["nes"] > 0, "UP_in_CLDN4_Q4", "DOWN_in_CLDN4_Q4")
    return g


def write_figures(
    headline: pd.DataFrame,
    cldn4: pd.Series,
    high: pd.Index,
    low: pd.Index,
    clin: pd.DataFrame,
) -> list[str]:
    h = headline[headline["rank_metric"] == "Welch_t_Q4_minus_Q1"].copy()
    h["label"] = h["term"].map(HEADLINE_LABEL)
    order = [HEADLINE_LABEL[t] for t in HEADLINE]
    h = h.set_index("label").reindex(order).reset_index()
    colors = ["#c0392b" if n > 0 else "#2980b9" for n in h["nes"]]

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    y = np.arange(len(h))
    ax.barh(y, h["nes"], color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(h["label"])
    ax.invert_yaxis()
    ax.set_xlabel("NES (positive = up in CLDN4 Q4)")
    ax.set_title(f"GSE126044 CLDN4 Q4 vs Q1 — headline prerank GSEA (n={len(high)} vs {len(low)})")
    for i, row in h.iterrows():
        txt = f"NES {row['nes']:+.2f}  FDR {fmt_p(row['fdr'])}"
        x = row["nes"]
        ax.text(
            x + (0.04 if x >= 0 else -0.04),
            i,
            txt,
            va="center",
            ha="left" if x >= 0 else "right",
            fontsize=8,
        )
    fig.tight_layout()
    png = FIG / "fig_headline_nes.png"
    pdf = FIG / "fig_headline_nes.pdf"
    fig.savefig(png, dpi=160)
    fig.savefig(pdf)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    rng = np.random.default_rng(SEED)
    for label, idx, color in (("Q1", low, "#2980b9"), ("Q4", high, "#c0392b")):
        vals = cldn4.loc[idx].to_numpy()
        axes[0].scatter(
            rng.normal({"Q1": 0, "Q4": 1}[label], 0.04, size=len(vals)),
            vals,
            color=color,
            edgecolor="black",
            linewidth=0.4,
            s=36,
            zorder=3,
            label=f"{label} n={len(idx)}",
        )
        axes[0].hlines(np.median(vals), {"Q1": -0.18, "Q4": 0.82}[label],
                       {"Q1": 0.18, "Q4": 1.18}[label], color="black", lw=1.2)
    axes[0].set_xticks([0, 1])
    axes[0].set_xticklabels(["CLDN4 Q1", "CLDN4 Q4"])
    axes[0].set_ylabel("CLDN4 log2(CPM+1)")
    axes[0].set_title("Honest n after quartile cut")
    axes[0].legend(frameon=False, fontsize=8)
    # annotate response / FFPE
    for i, sid in enumerate(list(low) + list(high)):
        x = 0 if sid in set(low) else 1
        resp = clin.loc[sid, "patient_response"] if sid in clin.index else "?"
        st = clin.loc[sid, "sample_type"] if sid in clin.index else "?"
        tag = ("R" if resp == "responder" else "NR") + ("/FFPE" if st == "FFPE" else "")
        axes[0].annotate(
            tag,
            (x, cldn4.loc[sid]),
            textcoords="offset points",
            xytext=(8, 0),
            fontsize=6,
            color="#555555",
        )
    axes[1].barh(y, h["nes"], color=colors, edgecolor="black", linewidth=0.4)
    axes[1].axvline(0, color="black", lw=0.7)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(h["label"])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("NES")
    axes[1].set_title("Headline NES")
    fig.tight_layout()
    png2 = FIG / "fig_qc_and_nes.png"
    pdf2 = FIG / "fig_qc_and_nes.pdf"
    fig.savefig(png2, dpi=160)
    fig.savefig(pdf2)
    plt.close(fig)
    return [
        str(png.relative_to(HERE.parent.parent)),
        str(png2.relative_to(HERE.parent.parent)),
    ]


def nes_md_rows(gsea: pd.DataFrame, rank_metric: str) -> str:
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(HEADLINE))].copy()
    sub["ord"] = sub["term"].map({t: i for i, t in enumerate(HEADLINE)})
    sub = sub.sort_values("ord")
    lines = []
    for _, r in sub.iterrows():
        lines.append(
            f"| {HEADLINE_LABEL[r['term']]} | `{r['term']}` | {fmt_nes(r['nes'])} | "
            f"{fmt_p(r['fdr'])} | {fmt_p(r['nom_p'])} | {int(r['n_set_in_rank'])} |"
        )
    return "\n".join(lines)


def call_headline(gsea: pd.DataFrame, rank_metric: str, terms=None) -> str:
    terms = list(terms) if terms is not None else list(HEADLINE)
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(terms))]
    bits = []
    for term in terms:
        r = sub[sub["term"] == term].iloc[0]
        bits.append(f"{HEADLINE_LABEL[term]} NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])}")
    return "; ".join(bits)


def write_finding(
    gsea: pd.DataFrame,
    sample_tab: pd.DataFrame,
    figures: list[str],
    n_ranked: int,
    n_q4: int,
    n_q1: int,
    q4_ids: list[str],
    q1_ids: list[str],
    cldn4_q4: float,
    cldn4_q1: float,
    n_q4_nr: int,
    n_q4_ffpe: int,
    n_q1_r: int,
) -> None:
    ifna = gsea[
        (gsea["rank_metric"] == "Welch_t_Q4_minus_Q1")
        & (gsea["term"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE")
    ]
    ifna_txt = "not ranked"
    if not ifna.empty:
        r = ifna.iloc[0]
        ifna_txt = (
            f"NES {fmt_nes(r['nes'])} nom p {fmt_p(r['nom_p'])} "
            f"(n={int(r['n_set_in_rank'])}; not in headline FDR)"
        )

    lead = {}
    for term in HEADLINE:
        r = gsea[(gsea["rank_metric"] == "Welch_t_Q4_minus_Q1") & (gsea["term"] == term)].iloc[0]
        lead[term] = r["lead_genes"]

    q4_str = ", ".join(q4_ids)
    q1_str = ", ".join(q1_ids)

    md = f"""# FINDING — GSE126044 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive only.** Public Cho 2020 anti-PD-1 NSCLC bulk RNA-seq. This folder does **not** audit or retract any slide. Single-gene CLDN4 vs ESTIMATE ImmuneScore is **already known** (Spearman r = **−0.524**, p = **0.037**, n = **16**; `methods/bulk_immune`) and is recorded, not re-tested. The extra is prerank GSEA on **CLDN4 Q4 vs Q1**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `scrna_pseudobulk_gsea_meta` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Honest n.** Cohort n = **16**. After the quartile cut, CLDN4 Q4 n = **{n_q4}**, Q1 n = **{n_q1}**. That is the contrast. Gene-set permutation on a prerank does **not** become a 16-sample test.

---

## 一句话 / TL;DR

| Contrast | Honest n | Known CLDN4 vs ImmuneScore | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE126044 CLDN4 Q4 vs Q1 | {n_q4} vs {n_q1} (cohort 16) | Spearman r=−0.524, p=0.037, n=16 (given; not re-tested) | {call_headline(gsea, "Welch_t_Q4_minus_Q1")} |

Primary rank is Welch *t* on log2(CPM+1). Q4 is **{n_q4_nr}/{n_q4} non-responder** and includes **{n_q4_ffpe} FFPE**; Q1 is **{n_q1_r}/{n_q1} responder**, all fresh. Bulk biopsy, so IFN / MHC NES can be infiltrate.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO counts only. No FASTQ / SRA. |
| File | `GSE126044_counts.txt.gz` |
| Cohort | 16 pre-treatment NSCLC biopsies, anti-PD-1 (Cho 2020, PMID 32879421) |
| Labels | 5 responder / 11 non-responder; 11 fresh / 5 FFPE (all FFPE are NR) |
| Target | CLDN4 log2(CPM+1) |
| Contrast | Q4 vs Q1 via `gsea_core.quartile_split` (≤P25 vs ≥P75) |
| Honest n | Q4 = {n_q4} ({q4_str}); Q1 = {n_q1} ({q1_str}) |
| Sign | t>0 / log2FC>0 = **up in CLDN4 Q4** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | mean log2(CPM+1) difference Q4−Q1; Spearman ρ vs continuous CLDN4 (n=16) |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | ImmuneScore re-test; TACSTD2 median A8 GSEA rerun; FASTQ |

Series: [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044). Title on GEO: *Genome-wide identification of differentially methylated promoters and enhancers associated with response to anti-PD-1 therapy in non-small cell lung cancer*.

---

## Known single-gene vs ImmuneScore (not re-tested)

From `methods/bulk_immune/results/demo/GSE126044_correlation_spearman.tsv` (same 16 samples, same CLDN4):

| Target | Score | n | Spearman r | p | FDR (within family) |
|---|---|---:|---:|---:|---:|
| CLDN4 | `estimate:ImmuneScore` | {KNOWN_IMMUNESCORE_N} | −0.524 | 0.037 | 0.142 |

This extra does not recompute ESTIMATE / MCP-counter / xCell.

---

## Headline NES (Welch *t* rank, Q4 minus Q1)

{n_ranked} genes ranked. CLDN4 itself: median log2(CPM+1) Q4 = **{cldn4_q4:.3f}**, Q1 = **{cldn4_q1:.3f}** (the split gene; not a GSEA claim).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "Welch_t_Q4_minus_Q1")}

Hallmark IFN-α (secondary, not in BH): {ifna_txt}.

Leading-edge (first 25, Welch *t* rank):

- IFN-γ: `{lead["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}`
- MHC-I / APM: `{lead["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}`
- KEGG TJ: `{lead["KEGG_TIGHT_JUNCTION"]}`

---

## Sensitivity: log2FC rank (Q4 − Q1 mean)

Same three sets, rank = mean log2(CPM+1) difference instead of Welch *t*. Still n = {n_q4} vs {n_q1}.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "log2FC_Q4_minus_Q1")}

---

## Sensitivity: Spearman vs continuous CLDN4 (uses all 16)

Rank = per-gene Spearman ρ vs CLDN4. This uses the full cohort (n=16) and is **not** a Q4 vs Q1 test. Reported so the extreme-quartile NES is not the only cut.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "Spearman_vs_CLDN4")}

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse126044_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `sample_table.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-test of CLDN4 vs ImmuneScore / CD8 / MCP NK (already in `methods/bulk_immune`).
- Not the A8 TACSTD2 **median** GSEA on this same series (that was 8 vs 8, TJ/keratin/EMT).
- Not FASTQ / salmon / DESeq2 from SRA.
- Not sample-permutation GSEA. Honest contrast n = **{n_q4} vs {n_q1}**; the engine is gene-set permutation on a prerank.
- Not a response GSEA. Q4 happens to be {n_q4_nr}/{n_q4} NR and includes {n_q4_ffpe} FFPE; that confounding is recorded, not adjusted.

---

## 中文摘要

只补公开 **GSE126044**（Cho 2020，抗 PD-1 NSCLC，n=16）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。CLDN4 与 ESTIMATE ImmuneScore 的单基因相关已经知道（Spearman r=−0.524，p=0.037，n=16），这里不重测。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 诚实 n：队列 16；四分位切割后 Q4 = {n_q4}，Q1 = {n_q1}。这不是 16 对 16 的样本置换检验。
- Q4 为 {n_q4_nr}/{n_q4} 非应答，含 {n_q4_ffpe} 例 FFPE；Q1 为 {n_q1_r}/{n_q1} 应答、全部新鲜组织。
- Headline（Welch *t* 秩）：{call_headline(gsea, "Welch_t_Q4_minus_Q1")}。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    counts = load_counts()
    clin = load_clinical()
    common = [c for c in counts.columns if c in clin.index]
    if len(common) != 16:
        raise SystemExit(f"expected 16 overlapping samples, got {len(common)}: {common}")
    counts = counts.loc[:, common]
    expr = log2cpm(counts)
    if "CLDN4" not in expr.index:
        raise SystemExit("CLDN4 missing from GSE126044 counts")
    cldn4 = expr.loc["CLDN4"]
    high, low = quartile_split(cldn4)
    n_q4, n_q1 = int(len(high)), int(len(low))
    log(f"cohort n={len(common)}; CLDN4 Q4 n={n_q4} {list(high)}; Q1 n={n_q1} {list(low)}")

    sample_tab = clin.loc[common].copy()
    sample_tab["CLDN4_log2cpm"] = cldn4.reindex(common)
    sample_tab["TACSTD2_log2cpm"] = expr.loc["TACSTD2"].reindex(common) if "TACSTD2" in expr.index else np.nan
    sample_tab["cldn4_arm"] = "mid"
    sample_tab.loc[list(high), "cldn4_arm"] = "Q4"
    sample_tab.loc[list(low), "cldn4_arm"] = "Q1"
    sample_tab = sample_tab.reset_index().rename(columns={"index": "sample_id"})

    rank_t = rank_high_vs_low(expr, high, low)
    rank_fc = rank_log2fc(expr, high, low)
    rank_rho = rank_spearman_vs_target(expr, cldn4)
    log(f"Welch t rank {len(rank_t)}; log2FC rank {len(rank_fc)}; Spearman rank {len(rank_rho)}")

    g_t = run_rank(rank_t, sets, "Welch_t_Q4_minus_Q1", n_q4, n_q1, n_q4 + n_q1)
    g_fc = run_rank(rank_fc, sets, "log2FC_Q4_minus_Q1", n_q4, n_q1, n_q4 + n_q1)
    g_rho = run_rank(rank_rho, sets, "Spearman_vs_CLDN4", n_q4, n_q1, 16)
    gsea = pd.concat([g_t, g_fc, g_rho], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    q4_ids = list(high)
    q1_ids = list(low)
    n_q4_nr = int((sample_tab.loc[sample_tab["cldn4_arm"] == "Q4", "patient_response"] == "non-responder").sum())
    n_q4_ffpe = int((sample_tab.loc[sample_tab["cldn4_arm"] == "Q4", "sample_type"] == "FFPE").sum())
    n_q1_r = int((sample_tab.loc[sample_tab["cldn4_arm"] == "Q1", "patient_response"] == "responder").sum())

    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "NSCLC pre-treatment biopsy, anti-PD-1 (Cho 2020)",
                "n_cohort": 16,
                "n_q4": n_q4,
                "n_q1": n_q1,
                "q4_samples": ",".join(q4_ids),
                "q1_samples": ",".join(q1_ids),
                "n_responder": 5,
                "n_nonresponder": 11,
                "n_fresh": 11,
                "n_ffpe": 5,
                "n_q4_nonresponder": n_q4_nr,
                "n_q4_ffpe": n_q4_ffpe,
                "n_q1_responder": n_q1_r,
                "rank_metric_primary": "Welch t on log2(CPM+1), CLDN4 Q4 minus Q1",
                "n_genes_ranked": int(rank_t.shape[0]),
                "counts_file": COUNTS_FILE,
                "known_immunescore_spearman_r": KNOWN_IMMUNESCORE_R,
                "known_immunescore_p": KNOWN_IMMUNESCORE_P,
                "known_immunescore_n": KNOWN_IMMUNESCORE_N,
                "known_immunescore_note": "given; not re-tested",
                "series_title": "anti-PD-1 NSCLC RNA-seq (Cho 2020)",
                "not": "ImmuneScore re-test; TACSTD2 median A8 GSEA; FASTQ; sample-permutation GSEA",
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    sample_tab.to_csv(TAB / "sample_table.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    rank_t.rename("welch_t").to_csv(TAB / "rank_welch_t.tsv", sep="\t")
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline, cldn4, high, low, clin)
    write_finding(
        gsea,
        sample_tab,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        q4_ids,
        q1_ids,
        float(cldn4.loc[high].median()),
        float(cldn4.loc[low].median()),
        n_q4_nr,
        n_q4_ffpe,
        n_q1_r,
    )

    print("\n=== HEADLINE NES (Welch t, Q4 vs Q1) ===")
    show = headline[headline["rank_metric"] == "Welch_t_Q4_minus_Q1"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print(f"\nHonest n: Q4={n_q4} vs Q1={n_q1} (cohort 16)")
    print(sample_tab.loc[sample_tab["cldn4_arm"].isin(["Q1", "Q4"]),
                         ["sample_id", "cldn4_arm", "CLDN4_log2cpm", "patient_response", "sample_type"]]
          .sort_values(["cldn4_arm", "CLDN4_log2cpm"])
          .to_string(index=False))


if __name__ == "__main__":
    main()
