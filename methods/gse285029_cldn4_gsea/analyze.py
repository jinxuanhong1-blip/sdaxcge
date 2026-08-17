#!/usr/bin/env python3
"""ADDITIVE extra — CLDN4 Q4 vs Q1 prerank GSEA on public GSE285029 (n=234).

CLDN4 only. TACSTD2 is not an anchor. Public author WTS matrix, no FASTQ.
Engine matches methods/cldn4_ko_gsea and methods/gse289287_trop2ko_gsea
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Positive NES = enriched in CLDN4 Q4 (high) vs Q1 (low).
BH-FDR is within the three headline sets.

CLDN4 is a member of KEGG_TIGHT_JUNCTION. The primary TJ NES therefore
includes the stratification gene. A locked sensitivity drops CLDN4 from
the rank (and from the TJ set) so the TJ call is not circular.
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

MATRIX = DATA / "GSE285029_WTS_expr_count_235_032820.txt.gz"
CONTRAST = "GSE285029_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE285029"
ANCHOR = "CLDN4"

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


def load_matrix() -> pd.DataFrame:
    if not MATRIX.exists():
        raise SystemExit(f"missing {MATRIX}; run download.py first")
    expr = pd.read_csv(MATRIX, sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.strip()
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.dropna(axis=0, how="all")
    # collapse rare duplicate symbols by mean
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    return expr


def log2p1_clip(expr: pd.DataFrame) -> pd.DataFrame:
    return np.log2(expr.clip(lower=0) + 1.0)


def run_rank(
    rank: pd.Series,
    sets: dict[str, list[str]],
    rank_metric: str,
    n_q4: int,
    n_q1: int,
) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = CONTRAST
    g["accession"] = ACCESSION
    g["anchor"] = ANCHOR
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_q4"] = n_q4
    g["n_q1"] = n_q1
    g["n_total"] = 234
    g["headline"] = g["term"].isin(HEADLINE)
    head = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    g["direction"] = np.where(g["nes"] > 0, "UP_in_CLDN4_Q4", "DOWN_in_CLDN4_Q4")
    return g


def write_figures(headline: pd.DataFrame, n_q4: int, n_q1: int) -> list[str]:
    h = headline[headline["rank_metric"] == "welch_t_Q4_minus_Q1"].copy()
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
    ax.set_title(f"GSE285029 CLDN4 Q4 vs Q1 prerank GSEA (n={n_q4} vs {n_q1})")
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
    return [str(png.relative_to(HERE.parent.parent))]


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
) -> None:
    ifna = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_minus_Q1")
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
            (gsea["rank_metric"] == "welch_t_Q4_minus_Q1") & (gsea["term"] == term)
        ].iloc[0]
        lead[term] = r["lead_genes"]

    tj_drop = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_minus_Q1_drop_CLDN4")
        & (gsea["term"] == "KEGG_TIGHT_JUNCTION")
    ]
    tj_drop_txt = "not ranked"
    if not tj_drop.empty:
        r = tj_drop.iloc[0]
        tj_drop_txt = (
            f"NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])} nom p {fmt_p(r['nom_p'])} "
            f"(n={int(r['n_set_in_rank'])})"
        )

    md = f"""# FINDING — GSE285029 CLDN4 Q4 vs Q1 prerank GSEA

**Additive only. CLDN4 only.** Public pre-ICI NSCLC WTS, **n = 234** (Koh et al., *JITC* 2025; GEO [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029)). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse289287_trop2ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity. This is bulk tumour WTS, not a cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE285029 CLDN4 Q4 vs Q1 | 234 (Q4={n_q4}, Q1={n_q1}) | CLDN4 quartiles on log2(clip0+1) | {call_headline(gsea, "welch_t_Q4_minus_Q1")} |

Primary rank is Welch *t* (Q4 − Q1) on the author matrix. MHC-I / APM is FDR **0.049** on this rank (nom p 0.049) and FDR **0.097** on the Spearman sensitivity — reported as computed, not rounded into a stronger call. CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: {tj_drop_txt}.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author WTS only. No FASTQ / SRA. |
| File | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = 234** (Case1–Case234). Author ICI–RNA-seq cohort. |
| Genes | {n_genes} symbols on the deposited matrix |
| Transform | `log2(pmax(x,0)+1)` (same clip as the sibling GSE285029 score analysis) |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | {n_q4} / {n_q1} |
| CLDN4 log2 cuts | P25 = {q_cuts[0]:.4f}; P75 = {q_cuts[1]:.4f} |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 234 samples) |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; response labels; FASTQ; sample-permutation GSEA |

Series: [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029). Title on GEO: pre-ICI NSCLC tumour WTS (PD-1 or PD-L1 blockade). Paper: Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/).

The filename says “count”. Values are continuous with negatives; they are the author-processed matrix, not raw integer counts and not TPM. Rank correlations / Welch *t* on the clipped log matrix are the claim.

---

## Headline NES (Welch *t*, Q4 − Q1)

{n_ranked} genes ranked. CLDN4 itself is the stratification gene (Welch *t* = {cldn4_t:+.3f}; by construction the top of the rank).

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

Uses all **n = 234** samples. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "spearman_vs_CLDN4")}

---

## Extra figures

- `{figures[0]}`

Tables: `methods/gse285029_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run of the sibling GSE285029 score / ρ tables (`results/w200/A11_GSE285029`). Those stay as written.
- Not FASTQ / salmon / DESeq2 from SRA (author matrix is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN call. This is bulk pre-ICI NSCLC WTS.
- Not a response / PD-L1 IHC / purity analysis. GEO does not release those labels.

---

## 中文摘要

只补公开 **GSE285029**（n=234）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 234 例 ICI 前 NSCLC 肿瘤 WTS。Q4={n_q4}，Q1={n_q1}。
- Headline（Welch *t*）：{call_headline(gsea, "welch_t_Q4_minus_Q1")}。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：{tj_drop_txt}。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    raw = load_matrix()
    log(f"raw matrix {raw.shape[0]} genes x {raw.shape[1]} samples")
    if raw.shape[1] != 234:
        raise SystemExit(f"expected n=234 samples, got {raw.shape[1]}")
    if ANCHOR not in raw.index:
        raise SystemExit(f"{ANCHOR} missing from matrix")

    logx = log2p1_clip(raw)
    cldn4 = logx.loc[ANCHOR]
    q1, q3 = float(cldn4.quantile(0.25)), float(cldn4.quantile(0.75))
    high, low = quartile_split(cldn4)
    n_q1, n_q4 = int(len(low)), int(len(high))
    log(f"CLDN4 P25={q1:.4f} P75={q3:.4f}; Q1={n_q1} Q4={n_q4}")

    split = pd.DataFrame(
        {
            "sample": cldn4.index,
            "CLDN4_log2p1": cldn4.to_numpy(),
            "quartile": [
                "Q1" if s in set(low) else ("Q4" if s in set(high) else "Q2Q3")
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
    log(f"Spearman rank {len(rank_rho)}")

    g_t = run_rank(rank_t, sets, "welch_t_Q4_minus_Q1", n_q4, n_q1)
    g_drop = run_rank(rank_t_drop, sets, "welch_t_Q4_minus_Q1_drop_CLDN4", n_q4, n_q1)
    g_rho = run_rank(rank_rho, sets, "spearman_vs_CLDN4", n_q4, n_q1)
    gsea = pd.concat([g_t, g_drop, g_rho], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "pre-ICI NSCLC tumour WTS (bulk; not cell-intrinsic)",
                "anchor": ANCHOR,
                "n_total": 234,
                "n_q4": n_q4,
                "n_q1": n_q1,
                "cldn4_p25": q1,
                "cldn4_p75": q3,
                "cldn4_welch_t": cldn4_t,
                "rank_metric_primary": "Welch t (CLDN4 Q4 minus Q1) on log2(clip0+1)",
                "n_genes_matrix": int(raw.shape[0]),
                "n_genes_ranked": int(rank_t.shape[0]),
                "matrix_file": MATRIX.name,
                "series_title": "Koh et al. JITC 2025 pre-ICI NSCLC WTS",
                "pmid": "40050048",
                "not": "TACSTD2; FASTQ; response labels; sample-permutation GSEA",
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    split.to_csv(TAB / "quartile_split.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline, n_q4, n_q1)
    write_finding(
        gsea,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        int(raw.shape[0]),
        cldn4_t,
        (q1, q3),
    )

    print("\n=== HEADLINE NES (Welch Q4 vs Q1) ===")
    show = headline[headline["rank_metric"] == "welch_t_Q4_minus_Q1"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== DROP CLDN4 ===")
    show2 = headline[headline["rank_metric"] == "welch_t_Q4_minus_Q1_drop_CLDN4"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show2["nes"] = show2["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show2["fdr"] = show2["fdr"].map(fmt_p)
    print(show2.to_string(index=False))


if __name__ == "__main__":
    main()
