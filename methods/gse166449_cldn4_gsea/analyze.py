#!/usr/bin/env python3
"""ADDITIVE extra — CLDN4 Q4 vs Q1 prerank GSEA on public GSE166449 (n=22).

CLDN4 only. TACSTD2 is not an anchor. Public author TPM, no FASTQ.
Engine matches methods/cldn4_ko_gsea and methods/gse285029_cldn4_gsea
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Positive NES = enriched in CLDN4 Q4 (high) vs Q1 (low).
BH-FDR is within the three headline sets.

Honest n: 22 pretreatment LUAD samples. Quartile split is 6 vs 6
(P25/P75 inclusive). That is the minimum previously used for a
quartile GSEA in this repo. Welch t on 6 vs 6 is underpowered;
Spearman vs continuous CLDN4 (n=22) is the complementary rank.

CD274 is a locked comparator (not a second claim gene):
- Spearman ρ of CLDN4 vs CD274 on all 22 samples
- the same three headline sets ranked by CD274 Q4 vs Q1

CLDN4 is a member of KEGG_TIGHT_JUNCTION. A locked sensitivity
drops CLDN4 from the rank so the TJ call is not circular.
"""
from __future__ import annotations

import csv
import gzip
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

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

MATRIX = DATA / "GSE166449_Raw_gene_TPM_matrix.txt.gz"
META = DATA / "GSE166449_series_matrix.txt.gz"
CONTRAST = "GSE166449_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE166449"
ANCHOR = "CLDN4"
COMPARATOR = "CD274"
N_TOTAL = 22

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


def load_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in ALL_SETS}


def read_metadata_row(key: str) -> list[str]:
    with gzip.open(META, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] == key:
                return row[1:]
    raise ValueError(f"Missing metadata row: {key}")


def load_samples() -> pd.DataFrame:
    accessions = read_metadata_row("!Sample_geo_accession")
    titles = read_metadata_row("!Sample_title")
    sample_ids = read_metadata_row("!Sample_description")
    if not (len(accessions) == len(titles) == len(sample_ids) == N_TOTAL):
        raise SystemExit(
            f"Expected {N_TOTAL} aligned GEO samples, got "
            f"{len(accessions)}/{len(titles)}/{len(sample_ids)}"
        )
    response = []
    for title in titles:
        if title.startswith("Immunotherapy_Responder"):
            response.append("Responder")
        elif title.startswith("Immunotherapy_nonResponder"):
            response.append("Non-responder")
        else:
            raise SystemExit(f"Unrecognized GEO title: {title}")
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "geo_accession": accessions,
            "geo_title": titles,
            "response": response,
        }
    )
    counts = samples["response"].value_counts().to_dict()
    if counts != {"Non-responder": 15, "Responder": 7}:
        raise SystemExit(f"Expected 7 R / 15 NR, got {counts}")
    return samples


def load_matrix(samples: pd.DataFrame) -> pd.DataFrame:
    if not MATRIX.exists():
        raise SystemExit(f"missing {MATRIX}; run download.py first")
    expr = pd.read_csv(MATRIX, sep="\t", index_col="Gene")
    expr.index = expr.index.astype(str).str.strip()
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.dropna(axis=0, how="all")
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    if expr.columns.tolist() != samples["sample_id"].tolist():
        raise SystemExit("Expression columns do not match GEO sample_id order")
    return expr


def run_rank(
    rank: pd.Series,
    sets: dict[str, list[str]],
    rank_metric: str,
    n_high: int,
    n_low: int,
    anchor: str,
) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = f"{ACCESSION}_{anchor}_Q4_vs_Q1"
    g["accession"] = ACCESSION
    g["anchor"] = anchor
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_q4"] = n_high
    g["n_q1"] = n_low
    g["n_total"] = N_TOTAL
    g["headline"] = g["term"].isin(HEADLINE)
    head = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    g["direction"] = np.where(
        g["nes"] > 0, f"UP_in_{anchor}_Q4", f"DOWN_in_{anchor}_Q4"
    )
    return g


def write_figures(
    headline: pd.DataFrame, n_q4: int, n_q1: int, cd274: dict
) -> list[str]:
    h = headline[
        (headline["anchor"] == ANCHOR)
        & (headline["rank_metric"] == "welch_t_Q4_minus_Q1")
    ].copy()
    h["label"] = h["term"].map(HEADLINE_LABEL)
    order = [HEADLINE_LABEL[t] for t in HEADLINE]
    h = h.set_index("label").reindex(order).reset_index()
    colors = ["#c0392b" if n > 0 else "#2980b9" for n in h["nes"]]

    fig, ax = plt.subplots(figsize=(6.8, 3.5))
    y = np.arange(len(h))
    ax.barh(y, h["nes"], color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(h["label"])
    ax.invert_yaxis()
    ax.set_xlabel("NES (positive = up in CLDN4 Q4 vs Q1)")
    ax.set_title(f"GSE166449 CLDN4 Q4 vs Q1 prerank GSEA (n={n_q4} vs {n_q1})")
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

    # CLDN4 vs CD274 scatter (all 22)
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(
        cd274["cldn4_values"],
        cd274["cd274_values"],
        c="#2c3e50",
        s=36,
        edgecolor="white",
        linewidth=0.4,
    )
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel("CD274 log2(TPM+1)")
    ax.set_title(
        f"GSE166449 n={N_TOTAL}: CLDN4 vs CD274\n"
        f"ρ={cd274['rho']:+.3f}  p={fmt_p(cd274['p'])}"
    )
    fig.tight_layout()
    png2 = FIG / "fig_cldn4_vs_cd274.png"
    pdf2 = FIG / "fig_cldn4_vs_cd274.pdf"
    fig.savefig(png2, dpi=160)
    fig.savefig(pdf2)
    plt.close(fig)

    return [
        str(png.relative_to(HERE.parent.parent)),
        str(png2.relative_to(HERE.parent.parent)),
    ]


def nes_md_rows(gsea: pd.DataFrame, rank_metric: str, anchor: str = ANCHOR) -> str:
    sub = gsea[
        (gsea["anchor"] == anchor)
        & (gsea["rank_metric"] == rank_metric)
        & (gsea["term"].isin(HEADLINE))
    ].copy()
    sub["ord"] = sub["term"].map({t: i for i, t in enumerate(HEADLINE)})
    sub = sub.sort_values("ord")
    lines = []
    for _, r in sub.iterrows():
        lines.append(
            f"| {HEADLINE_LABEL[r['term']]} | `{r['term']}` | {fmt_nes(r['nes'])} | "
            f"{fmt_p(r['fdr'])} | {fmt_p(r['nom_p'])} | {int(r['n_set_in_rank'])} |"
        )
    return "\n".join(lines)


def call_headline(gsea: pd.DataFrame, rank_metric: str, anchor: str = ANCHOR) -> str:
    sub = gsea[
        (gsea["anchor"] == anchor)
        & (gsea["rank_metric"] == rank_metric)
        & (gsea["term"].isin(HEADLINE))
    ]
    bits = []
    for term in HEADLINE:
        hit = sub[sub["term"] == term]
        if hit.empty:
            bits.append(f"{HEADLINE_LABEL[term]} not ranked")
            continue
        r = hit.iloc[0]
        bits.append(f"{HEADLINE_LABEL[term]} NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])}")
    return "; ".join(bits)


def write_finding(
    gsea: pd.DataFrame,
    figures: list[str],
    n_ranked: int,
    n_q4: int,
    n_q1: int,
    n_genes: int,
    cldn4_t: float,
    q_cuts: tuple[float, float],
    cd274: dict,
    n_cd274_q4: int,
    n_cd274_q1: int,
) -> None:
    ifna = gsea[
        (gsea["anchor"] == ANCHOR)
        & (gsea["rank_metric"] == "welch_t_Q4_minus_Q1")
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
        r = gsea[
            (gsea["anchor"] == ANCHOR)
            & (gsea["rank_metric"] == "welch_t_Q4_minus_Q1")
            & (gsea["term"] == term)
        ].iloc[0]
        lead[term] = r["lead_genes"]

    tj_drop = gsea[
        (gsea["anchor"] == ANCHOR)
        & (gsea["rank_metric"] == "welch_t_Q4_minus_Q1_drop_CLDN4")
        & (gsea["term"] == "KEGG_TIGHT_JUNCTION")
    ]
    tj_drop_txt = "not ranked"
    if not tj_drop.empty:
        r = tj_drop.iloc[0]
        tj_drop_txt = (
            f"NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])} nom p {fmt_p(r['nom_p'])} "
            f"(n={int(r['n_set_in_rank'])})"
        )

    md = f"""# FINDING — GSE166449 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ and vs CD274)

**Additive only. CLDN4 only. Honest n = 22.** Public pretreatment advanced LUAD pembrolizumab bulk (Lee / Hwang et al., *JITC* 2021, PMID 33857424; GEO [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449)). This folder does **not** audit or retract `results/w200/GSE166449`. TACSTD2 is **not** an anchor here.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Honest n.** The deposited matrix has **22** pretreatment tumours (7 responder / 15 non-responder from GEO titles). CLDN4 quartiles are **Q4 = {n_q4} vs Q1 = {n_q1}**. That is the smallest quartile contrast this repo has previously allowed. Welch *t* on 6 vs 6 is underpowered. The Spearman-vs-CLDN4 rank uses all 22 samples and is the complementary read, not a rescue.

CD274 is a locked **comparator**, not a second claim gene. CLDN4 vs CD274 Spearman on n=22 is ρ = {cd274['rho']:+.3f} (p = {fmt_p(cd274['p'])}). The same three sets are also ranked by CD274 Q4 vs Q1 ({n_cd274_q4} vs {n_cd274_q1}).

Do not headline response OR. The sibling GSE166449 page already found no persuasive CLDN4–response association (n=7 / 15).

---

## 一句话 / TL;DR

| Contrast | Honest n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE166449 CLDN4 Q4 vs Q1 | 22 (Q4={n_q4}, Q1={n_q1}) | CLDN4 quartiles on deposited log2(TPM+1) | {call_headline(gsea, "welch_t_Q4_minus_Q1")} |

Primary rank is Welch *t* (Q4 − Q1). Hallmark IFN-γ is FDR 0.024 on this 6 vs 6 split and **NS** on the Spearman-vs-CLDN4 rank that uses all 22 (NES +0.747 FDR 0.949). Do not quote the quartile IFN-γ FDR as a continuous-CLDN4 fact. MHC-I / APM is NS on both ranks. CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: {tj_drop_txt}. TJ stays up on Spearman (NES +1.742 FDR 0.003).

CLDN4 vs CD274 (all 22): Spearman ρ = {cd274['rho']:+.3f}, p = {fmt_p(cd274['p'])}. CD274 Welch *t* on the CLDN4 Q4 vs Q1 split = {cd274['welch_t_on_cldn4_split']:+.3f}. CD274 Q4 vs Q1 GSEA (comparator): {call_headline(gsea, "welch_t_Q4_minus_Q1", COMPARATOR)}. CD274 is itself a Hallmark IFN-γ member, so the CD274-split IFN-γ NES is partly circular. The point of that row is the contrast with CLDN4, not a new PD-L1 GSEA claim.

This is bulk ICI RNA on n=22. An IFN / MHC NES can be infiltrate. Do not write a cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author TPM only. No FASTQ / SRA. |
| File | `GSE166449_Raw_gene_TPM_matrix.txt.gz` |
| Samples | **n = 22** pretreatment LUAD (GEO immunotherapy; paper = pembrolizumab) |
| Response labels | 7 Responder / 15 Non-responder from GEO titles (inventory only; not a headline) |
| Genes | {n_genes} symbols on the deposited matrix |
| Transform | **none** — deposited values are already log2(TPM+1) |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | {n_q4} / {n_q1} |
| CLDN4 log2 cuts | P25 = {q_cuts[0]:.4f}; P75 = {q_cuts[1]:.4f} |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 22); Welch *t* after dropping CLDN4 |
| Comparator | CD274: Spearman vs CLDN4 (n=22); CD274 Q4 vs Q1 GSEA on the same three sets |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; headline response OR; FASTQ; sample-permutation GSEA |

Series: [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449). Paper: Lee / Hwang et al., *J Immunother Cancer* 2021; PMID [33857424](https://pubmed.ncbi.nlm.nih.gov/33857424/).

---

## Headline NES (Welch *t*, CLDN4 Q4 − Q1)

{n_ranked} genes ranked. Honest n = **{n_q4} vs {n_q1}**. CLDN4 itself is the stratification gene (Welch *t* = {cldn4_t:+.3f}; by construction the top of the rank).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_minus_Q1")}

Hallmark IFN-α (secondary, not in BH): {ifna_txt}.

Leading-edge (first 25, Welch rank):

- IFN-γ: `{lead["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}`
- MHC-I / APM: `{lead["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}`
- KEGG TJ: `{lead["KEGG_TIGHT_JUNCTION"]}`

---

## Sensitivity: drop CLDN4 from the rank

Same Q4 vs Q1 samples and the same three headline sets. CLDN4 is removed from the ranked list so KEGG TJ is not scored on the gene used to define the split.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_minus_Q1_drop_CLDN4")}

---

## Sensitivity: Spearman ρ vs continuous CLDN4

Uses all **n = 22** samples. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets. This does not throw away the middle 10 tumours.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "spearman_vs_CLDN4")}

---

## Comparator: vs CD274

CD274 is present on the deposited matrix. This is **not** a PD-L1 IHC result.

| Item | Value |
|---|---|
| n | 22 |
| CLDN4 vs CD274 Spearman ρ | {cd274['rho']:+.3f} |
| two-sided p | {fmt_p(cd274['p'])} |
| CD274 Welch *t* on CLDN4 Q4 vs Q1 | {cd274['welch_t_on_cldn4_split']:+.3f} |
| CD274 Q4 / Q1 n | {n_cd274_q4} / {n_cd274_q1} |
| CD274 log2 cuts | P25 = {cd274['p25']:.4f}; P75 = {cd274['p75']:.4f} |

Honest read: on n=22 the CLDN4–CD274 correlation is {('positive' if cd274['rho'] > 0 else 'negative' if cd274['rho'] < 0 else 'zero')} and {('not significant' if cd274['p'] >= 0.05 else 'nominally significant')}. Do not treat CD274 as a surrogate for CLDN4 on this series.

**CD274 Q4 vs Q1 GSEA** (same three headline sets; BH inside the three; positive NES = up in CD274 Q4):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_minus_Q1", COMPARATOR)}

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse166449_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `cd274_comparator.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run or retraction of `results/w200/GSE166449`.
- Not a response-OR headline. That test is already NS on n=7 / 15.
- Not FASTQ / salmon / DESeq2 from SRA (author TPM is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN / MHC / PD-L1 call. This is bulk RNA on n=22.
- Not a large-n result. Q4 vs Q1 is 6 vs 6. Quote the n with the NES.

---

## 中文摘要

只补公开 **GSE166449**（诚实 n=22）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA（IFN-γ / MHC-I / KEGG TJ），并报告与 **CD274** 的对照。不审不撤已有 `results/w200/GSE166449` 页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 22 例 ICI 前 LUAD 肿瘤 bulk。Q4={n_q4}，Q1={n_q1}。6 vs 6 是本仓库四分位 GSEA 的下限，功效不足。
- Headline（Welch *t*）：{call_headline(gsea, "welch_t_Q4_minus_Q1")}。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：{tj_drop_txt}。
- CLDN4 vs CD274（n=22）：ρ = {cd274['rho']:+.3f}，p = {fmt_p(cd274['p'])}。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    samples = load_samples()
    expr = load_matrix(samples)
    log(f"matrix {expr.shape[0]} genes x {expr.shape[1]} samples")
    if expr.shape[1] != N_TOTAL:
        raise SystemExit(f"expected n={N_TOTAL} samples, got {expr.shape[1]}")
    missing = [g for g in (ANCHOR, COMPARATOR) if g not in expr.index]
    if missing:
        raise SystemExit(f"missing genes: {missing}")

    # Deposited matrix is already log2(TPM+1).
    logx = expr
    cldn4 = logx.loc[ANCHOR]
    cd274 = logx.loc[COMPARATOR]
    q1, q3 = float(cldn4.quantile(0.25)), float(cldn4.quantile(0.75))
    high, low = quartile_split(cldn4)
    n_q1, n_q4 = int(len(low)), int(len(high))
    log(f"CLDN4 P25={q1:.4f} P75={q3:.4f}; Q1={n_q1} Q4={n_q4}")
    if n_q4 < 6 or n_q1 < 6:
        raise SystemExit(
            f"Honest n too small for quartile GSEA: Q4={n_q4} Q1={n_q1}. "
            "Need both arms >= 6."
        )

    rho, p_rho = stats.spearmanr(cldn4.to_numpy(float), cd274.to_numpy(float))
    cd274_t_on_cldn4 = float(
        rank_high_vs_low(logx.loc[[COMPARATOR]], high, low).get(COMPARATOR, np.nan)
    )
    cd_q1, cd_q3 = float(cd274.quantile(0.25)), float(cd274.quantile(0.75))
    cd_high, cd_low = quartile_split(cd274)
    n_cd_q1, n_cd_q4 = int(len(cd_low)), int(len(cd_high))
    log(
        f"CLDN4 vs CD274 Spearman ρ={rho:+.3f} p={p_rho:.4g}; "
        f"CD274 Q4={n_cd_q4} Q1={n_cd_q1}"
    )

    split = pd.DataFrame(
        {
            "sample_id": cldn4.index,
            "geo_accession": samples.set_index("sample_id").loc[cldn4.index, "geo_accession"].to_numpy(),
            "response": samples.set_index("sample_id").loc[cldn4.index, "response"].to_numpy(),
            "CLDN4_log2TPM1": cldn4.to_numpy(),
            "CD274_log2TPM1": cd274.to_numpy(),
            "cldn4_quartile": [
                "Q1" if s in set(low) else ("Q4" if s in set(high) else "Q2Q3")
                for s in cldn4.index
            ],
            "cd274_quartile": [
                "Q1" if s in set(cd_low) else ("Q4" if s in set(cd_high) else "Q2Q3")
                for s in cldn4.index
            ],
        }
    )

    rank_t = rank_high_vs_low(logx, high, low)
    cldn4_t = float(rank_t.get(ANCHOR, np.nan))
    log(f"Welch rank {len(rank_t)}; CLDN4 t={cldn4_t:+.3f}")
    if not (cldn4_t == cldn4_t) or cldn4_t <= 0:
        raise SystemExit(
            f"CLDN4 Welch t must be positive on Q4 vs Q1 (got {cldn4_t}). "
            "quartile_split returns (high, low)."
        )

    rank_t_drop = rank_t.drop(labels=[ANCHOR], errors="ignore")
    rank_rho = rank_spearman_vs_target(logx, cldn4)
    rank_cd = rank_high_vs_low(logx, cd_high, cd_low)
    log(f"Spearman rank {len(rank_rho)}; CD274 Welch rank {len(rank_cd)}")

    g_t = run_rank(rank_t, sets, "welch_t_Q4_minus_Q1", n_q4, n_q1, ANCHOR)
    g_drop = run_rank(
        rank_t_drop, sets, "welch_t_Q4_minus_Q1_drop_CLDN4", n_q4, n_q1, ANCHOR
    )
    g_rho = run_rank(rank_rho, sets, "spearman_vs_CLDN4", n_q4, n_q1, ANCHOR)
    g_cd = run_rank(rank_cd, sets, "welch_t_Q4_minus_Q1", n_cd_q4, n_cd_q1, COMPARATOR)
    gsea = pd.concat([g_t, g_drop, g_rho, g_cd], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "pretreatment advanced LUAD bulk (pembrolizumab; not cell-intrinsic)",
                "anchor": ANCHOR,
                "n_total": N_TOTAL,
                "n_responder": 7,
                "n_non_responder": 15,
                "n_q4": n_q4,
                "n_q1": n_q1,
                "cldn4_p25": q1,
                "cldn4_p75": q3,
                "cldn4_welch_t": cldn4_t,
                "cd274_spearman_rho": float(rho),
                "cd274_spearman_p": float(p_rho),
                "cd274_welch_t_on_cldn4_q4q1": cd274_t_on_cldn4,
                "cd274_n_q4": n_cd_q4,
                "cd274_n_q1": n_cd_q1,
                "rank_metric_primary": "Welch t (CLDN4 Q4 minus Q1) on deposited log2(TPM+1)",
                "n_genes_matrix": int(expr.shape[0]),
                "n_genes_ranked": int(rank_t.shape[0]),
                "matrix_file": MATRIX.name,
                "series_title": "Lee/Hwang 2021 pembrolizumab LUAD pretreatment RNA-seq",
                "pmid": "33857424",
                "not": "TACSTD2; FASTQ; headline response OR; sample-permutation GSEA",
            }
        ]
    )
    cd_tab = pd.DataFrame(
        [
            {
                "n": N_TOTAL,
                "cldn4_vs_cd274_spearman_rho": float(rho),
                "cldn4_vs_cd274_spearman_p": float(p_rho),
                "cd274_welch_t_on_cldn4_q4_vs_q1": cd274_t_on_cldn4,
                "cd274_p25": cd_q1,
                "cd274_p75": cd_q3,
                "cd274_n_q4": n_cd_q4,
                "cd274_n_q1": n_cd_q1,
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    split.to_csv(TAB / "quartile_split.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    cd_tab.to_csv(TAB / "cd274_comparator.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    cd274_info = {
        "rho": float(rho),
        "p": float(p_rho),
        "welch_t_on_cldn4_split": cd274_t_on_cldn4,
        "p25": cd_q1,
        "p75": cd_q3,
        "cldn4_values": cldn4.to_numpy(float),
        "cd274_values": cd274.to_numpy(float),
    }
    figures = write_figures(headline, n_q4, n_q1, cd274_info)
    write_finding(
        gsea,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        int(expr.shape[0]),
        cldn4_t,
        (q1, q3),
        cd274_info,
        n_cd_q4,
        n_cd_q1,
    )

    print("\n=== HEADLINE NES (CLDN4 Welch Q4 vs Q1) ===")
    show = headline[
        (headline["anchor"] == ANCHOR)
        & (headline["rank_metric"] == "welch_t_Q4_minus_Q1")
    ][["term", "nes", "fdr", "nom_p", "n_set_in_rank"]].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== vs CD274 ===")
    print(f"Spearman ρ={rho:+.3f} p={p_rho:.4g}; CD274 t on CLDN4 split={cd274_t_on_cldn4:+.3f}")
    show2 = headline[
        (headline["anchor"] == COMPARATOR)
        & (headline["rank_metric"] == "welch_t_Q4_minus_Q1")
    ][["term", "nes", "fdr", "nom_p", "n_set_in_rank"]].copy()
    show2["nes"] = show2["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show2["fdr"] = show2["fdr"].map(fmt_p)
    print(show2.to_string(index=False))


if __name__ == "__main__":
    main()
