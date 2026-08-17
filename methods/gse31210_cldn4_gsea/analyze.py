#!/usr/bin/env python3
"""ADDITIVE extra — CLDN4 Q4 vs Q1 prerank GSEA on public GSE31210 (n=226).

CLDN4 only. TACSTD2 is not an anchor. CD8A ρ is already known in
methods/gse31210_cldn4_immune and is not re-headlined here.

Public GEO series-matrix MAS5 + GPL570.annot. Unique-mapped max-mean collapse.
Engine matches methods/cldn4_ko_gsea and methods/gse285029_cldn4_gsea
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Positive NES = enriched in CLDN4 Q4 (high) vs Q1 (low).
BH-FDR is within the three headline sets.

CLDN4 is a member of KEGG_TIGHT_JUNCTION. The primary TJ NES therefore
includes the stratification gene. A locked sensitivity drops CLDN4 from
the rank (and from the TJ set) so the TJ call is not circular.
"""
from __future__ import annotations

import gzip
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

MATRIX = DATA / "GSE31210_series_matrix.txt.gz"
ANNOT = DATA / "GPL570.annot.gz"
CONTRAST = "GSE31210_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE31210"
ANCHOR = "CLDN4"
N_TUMOR_EXPECTED = 226
N_MATRIX_EXPECTED = 246

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


def load_gpl_annot(path: Path) -> pd.Series:
    """Unique-mapped first-symbol map. Multi-mapped /// probes dropped."""
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
    return unique


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out, pick


def load_tumor_matrix() -> tuple[pd.DataFrame, pd.Series, dict]:
    if not MATRIX.exists():
        raise SystemExit(f"missing {MATRIX}; run download.py first")
    if not ANNOT.exists():
        raise SystemExit(f"missing {ANNOT}; run download.py first")
    log("parse series matrix")
    meta, probes = parse_geo_matrix(MATRIX)
    unique_map = load_gpl_annot(ANNOT)
    if "tissue" not in meta.columns:
        raise SystemExit(f"GSE31210 tissue characteristic missing: {list(meta.columns)}")
    tissue = meta["tissue"].astype(str)
    is_tumor = tissue.str.lower().str.contains("primary lung tumor")
    n_matrix = int(len(meta))
    n_tumor = int(is_tumor.sum())
    log(f"arrays={n_matrix} tumors={n_tumor}")
    if n_matrix != N_MATRIX_EXPECTED or n_tumor != N_TUMOR_EXPECTED:
        raise SystemExit(
            f"honest n mismatch: arrays={n_matrix} (expected {N_MATRIX_EXPECTED}), "
            f"tumors={n_tumor} (expected {N_TUMOR_EXPECTED})"
        )
    probes_t = probes.loc[:, meta.index[is_tumor]]
    genes, probe_pick = collapse_maxmean(probes_t, unique_map)
    if ANCHOR not in genes.index:
        raise SystemExit(f"{ANCHOR} missing after unique-probe collapse")
    probe_used = {ANCHOR: str(probe_pick.loc[ANCHOR])}
    inv = {
        "n_matrix": n_matrix,
        "n_tumor": n_tumor,
        "n_nontumor": n_matrix - n_tumor,
        "n_genes_unique": int(genes.shape[0]),
        "cldn4_probe": probe_used[ANCHOR],
        "tumor_rule": "tissue == primary lung tumor",
    }
    return genes, probe_pick, inv


def run_rank(
    rank: pd.Series,
    sets: dict[str, list[str]],
    rank_metric: str,
    n_q4: int,
    n_q1: int,
    n_total: int,
) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = CONTRAST
    g["accession"] = ACCESSION
    g["anchor"] = ANCHOR
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_q4"] = n_q4
    g["n_q1"] = n_q1
    g["n_total"] = n_total
    g["headline"] = g["term"].isin(HEADLINE)
    head = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    g["direction"] = np.where(g["nes"] > 0, "UP_in_CLDN4_Q4", "DOWN_in_CLDN4_Q4")
    return g


def write_figures(headline: pd.DataFrame, n_q4: int, n_q1: int, n_total: int) -> list[str]:
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
    ax.set_title(f"GSE31210 CLDN4 Q4 vs Q1 prerank GSEA (n={n_total}; {n_q4} vs {n_q1})")
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
    n_total: int,
    n_matrix: int,
    cldn4_t: float,
    q_cuts: tuple[float, float],
    cldn4_probe: str,
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

    def pick(rank_metric: str, term: str):
        hit = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"] == term)]
        return hit.iloc[0] if len(hit) else None

    welch_ifn = pick("welch_t_Q4_minus_Q1", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
    welch_mhc = pick("welch_t_Q4_minus_Q1", "CUSTOM_MHC_I_ANTIGEN_PRESENTATION")
    welch_tj = pick("welch_t_Q4_minus_Q1", "KEGG_TIGHT_JUNCTION")
    rho_mhc = pick("spearman_vs_CLDN4", "CUSTOM_MHC_I_ANTIGEN_PRESENTATION")
    tj_drop = pick("welch_t_Q4_minus_Q1_drop_CLDN4", "KEGG_TIGHT_JUNCTION")
    tj_drop_txt = "not ranked"
    if tj_drop is not None:
        tj_drop_txt = (
            f"NES {fmt_nes(tj_drop['nes'])} FDR {fmt_p(tj_drop['fdr'])} "
            f"nom p {fmt_p(tj_drop['nom_p'])} (n={int(tj_drop['n_set_in_rank'])})"
        )
    ifn_txt = (
        f"NES {fmt_nes(welch_ifn['nes'])} FDR {fmt_p(welch_ifn['fdr'])}"
        if welch_ifn is not None
        else "not ranked"
    )
    mhc_txt = (
        f"FDR {fmt_p(welch_mhc['fdr'])}"
        if welch_mhc is not None
        else "not ranked"
    )
    mhc_rho_txt = (
        f"FDR {fmt_p(rho_mhc['fdr'])}"
        if rho_mhc is not None
        else "NA"
    )
    tj_txt = (
        f"NES {fmt_nes(welch_tj['nes'])} FDR {fmt_p(welch_tj['fdr'])}"
        if welch_tj is not None
        else "not ranked"
    )

    md = f"""# FINDING — GSE31210 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive extra. CLDN4 only.** Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

**CD8 ρ is already known** in `methods/gse31210_cldn4_immune` (CLDN4 vs CD8A ρ=−0.341, n=226) and is **not** re-tested or re-headlined. This folder only adds prerank GSEA of IFN-γ / MHC-I / KEGG TJ on the same tumor set.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

No ICI labels. Bulk MAS5 microarray, not a cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE31210 CLDN4 Q4 vs Q1 | {n_total} (Q4={n_q4}, Q1={n_q1}) | CLDN4 quartiles on log2(MAS5+1) | {call_headline(gsea, "welch_t_Q4_minus_Q1")} |

Primary rank is Welch *t* (Q4 − Q1) on unique-mapped max-mean collapse.

**Honest read:** Hallmark IFN-γ is **down** in CLDN4 Q4 ({ifn_txt}). That is the same direction as the already-known CLDN4–CD8A anti-correlation, but it is a program NES, not a CD8 re-test. MHC-I / APM is the same sign and **NS** ({mhc_txt}; Spearman sensitivity {mhc_rho_txt}) — reported as computed, not rounded into a call. KEGG TJ is up ({tj_txt}). CLDN4 is a member of that set, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: {tj_drop_txt}.

---

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | {n_matrix} | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **{n_total}** | `tissue: primary lung tumor` |
| adjacent / non-tumor | {n_matrix - n_total} | dropped |
| CLDN4 Q4 / Q1 | {n_q4} / {n_q1} | ≤P25 vs ≥P75; ties stay in the tail |
| genes after unique-mapped max-mean | {n_genes} | multi-mapped `///` probes dropped |
| genes ranked (Welch) | {n_ranked} | finite Welch *t*, Q4 minus Q1 |

Same tumor rule as `methods/gse31210_cldn4_immune`. Not a new cohort hunt.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL / FASTQ / SRA. |
| File | `GSE31210_series_matrix.txt.gz` + `GPL570.annot.gz` |
| Samples | **n = {n_total}** primary lung tumors (of {n_matrix} arrays) |
| Platform | GPL570 Affymetrix U133 Plus 2.0 |
| Collapse | unique-mapped max-mean (same as sibling immune folder) |
| CLDN4 probe | `{cldn4_probe}` |
| Transform | `log2(MAS5+1)` |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | {n_q4} / {n_q1} |
| CLDN4 log2 cuts | P25 = {q_cuts[0]:.4f}; P75 = {q_cuts[1]:.4f} |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all {n_total} tumors) |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | CD8 ρ (already known); TACSTD2 ranks; ICI labels; sample-permutation GSEA |

Series: [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210). Okayama et al., *Cancer Res* 2012; PMID [22261853](https://pubmed.ncbi.nlm.nih.gov/22261853/). Treatment-naive surgical LUAD. No ICI arm.

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

Uses all **n = {n_total}** tumors. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "spearman_vs_CLDN4")}

---

## Extra figures

- `{figures[0]}`

Tables: `methods/gse31210_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `probe_used.tsv`.

---

## What this is not

- Not a re-run of the sibling CD8 / ImmuneScore / TACSTD2 ρ tables (`methods/gse31210_cldn4_immune`). Those stay as written. CD8 ρ is already known.
- Not a TACSTD2 analysis. CLDN4 only.
- Not CEL / RMA / FASTQ from SRA (series-matrix MAS5 is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN call. This is bulk LUAD microarray.
- Not an ICI-response analysis. Labels are not deposited.

---

## 中文摘要

只补公开 **GSE31210**（诚实 n=226 原发肺腺癌）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。CD8 ρ 已在 `gse31210_cldn4_immune` 给出，这里不再当主结果。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 246 张芯片，**226** 例 `primary lung tumor`，20 例非肿瘤丢掉。Q4={n_q4}，Q1={n_q1}。
- Headline（Welch *t*）：{call_headline(gsea, "welch_t_Q4_minus_Q1")}。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：{tj_drop_txt}。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    raw, probe_pick, inv = load_tumor_matrix()
    n_total = int(inv["n_tumor"])
    n_matrix = int(inv["n_matrix"])
    log(f"tumor matrix {raw.shape[0]} genes x {raw.shape[1]} samples")
    if raw.shape[1] != n_total:
        raise SystemExit(f"expected n={n_total} tumors, got {raw.shape[1]}")

    logx = np.log2(raw.clip(lower=0) + 1.0)
    cldn4 = logx.loc[ANCHOR]
    q1, q3 = float(cldn4.quantile(0.25)), float(cldn4.quantile(0.75))
    high, low = quartile_split(cldn4)
    n_q1, n_q4 = int(len(low)), int(len(high))
    log(f"CLDN4 P25={q1:.4f} P75={q3:.4f}; Q1={n_q1} Q4={n_q4}")

    split = pd.DataFrame(
        {
            "sample": cldn4.index,
            "CLDN4_log2_MAS5p1": cldn4.to_numpy(),
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

    g_t = run_rank(rank_t, sets, "welch_t_Q4_minus_Q1", n_q4, n_q1, n_total)
    g_drop = run_rank(rank_t_drop, sets, "welch_t_Q4_minus_Q1_drop_CLDN4", n_q4, n_q1, n_total)
    g_rho = run_rank(rank_rho, sets, "spearman_vs_CLDN4", n_q4, n_q1, n_total)
    gsea = pd.concat([g_t, g_drop, g_rho], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    inv_row = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "stage I-II LUAD primary tumor microarray (bulk; not cell-intrinsic)",
                "platform": "GPL570 U133 Plus 2.0",
                "anchor": ANCHOR,
                "cldn4_probe": inv["cldn4_probe"],
                "n_matrix": n_matrix,
                "n_total": n_total,
                "n_nontumor_dropped": n_matrix - n_total,
                "n_q4": n_q4,
                "n_q1": n_q1,
                "cldn4_p25": q1,
                "cldn4_p75": q3,
                "cldn4_welch_t": cldn4_t,
                "rank_metric_primary": "Welch t (CLDN4 Q4 minus Q1) on log2(MAS5+1)",
                "n_genes_matrix": int(raw.shape[0]),
                "n_genes_ranked": int(rank_t.shape[0]),
                "matrix_file": MATRIX.name,
                "annot_file": ANNOT.name,
                "series_title": "Okayama et al. Cancer Res 2012 Japanese stage I-II LUAD",
                "pmid": "22261853",
                "cd8_rho_already_known": "methods/gse31210_cldn4_immune CLDN4 vs CD8A rho=-0.341 n=226",
                "not": "CD8 re-test; TACSTD2; CEL/RMA; ICI labels; sample-permutation GSEA",
            }
        ]
    )

    probe_tab = pd.DataFrame(
        [{"gene": ANCHOR, "probe_id": inv["cldn4_probe"], "collapse": "unique-mapped max-mean"}]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    split.to_csv(TAB / "quartile_split.tsv", sep="\t", index=False)
    inv_row.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    probe_tab.to_csv(TAB / "probe_used.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline, n_q4, n_q1, n_total)
    write_finding(
        gsea,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        int(raw.shape[0]),
        n_total,
        n_matrix,
        cldn4_t,
        (q1, q3),
        inv["cldn4_probe"],
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
    print("\n=== SPEARMAN vs CLDN4 ===")
    show3 = headline[headline["rank_metric"] == "spearman_vs_CLDN4"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show3["nes"] = show3["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show3["fdr"] = show3["fdr"].map(fmt_p)
    print(show3.to_string(index=False))


if __name__ == "__main__":
    main()
