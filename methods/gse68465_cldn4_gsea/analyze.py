#!/usr/bin/env python3
"""ADDITIVE extra — CLDN4 Q4 vs Q1 prerank GSEA on public GSE68465 LUAD.

CLDN4 only. TACSTD2 is not an anchor. Public GEO series matrix (MAS5) +
official GPL96 annotation. No CEL reprocess. No ICI labels.

Engine matches methods/cldn4_ko_gsea / methods/gse285029_cldn4_gsea
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Positive NES = enriched in CLDN4 Q4 (high) vs Q1 (low).
BH-FDR is within the three headline sets.

Honest n: GEO deposits 462 arrays (443 LUAD + 19 Normal). The GSEA
contrast is the quartile arms of the 443 tumours, not n=462 and not
the paper's 442.

CLDN4 is a member of KEGG_TIGHT_JUNCTION. The primary TJ NES therefore
includes the stratification gene. A locked sensitivity drops CLDN4 from
the rank so the TJ call is not circular.

U133A does not carry every current MSigDB symbol. n in rank is the
intersection after a small documented HGNC→U133A alias map
(WARS1→WARS, MARCHF1→MARCH1, PALS1→MPP5).
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

MATRIX = DATA / "GSE68465_series_matrix.txt.gz"
ANNOT = DATA / "GPL96.annot.gz"
CONTRAST = "GSE68465_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE68465"
ANCHOR = "CLDN4"
NAMED_CLDN4 = "201428_at"

# Current MSigDB / HGNC symbols that official GPL96 still lists under
# the U133A-era name. Applied only when matching sets to the collapsed
# matrix. Not a fishing remap.
U133A_ALIASES = {
    "WARS1": "WARS",
    "MARCHF1": "MARCH1",
    "PALS1": "MPP5",
}

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


def parse_geo_matrix(path: Path):
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
            cleaned = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
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


def load_gpl_annot(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as f:
        skip = 0
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith("ID "):
                skip = i
                break
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out


def load_sets_raw() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in ALL_SETS}


def remap_sets(raw: dict[str, list[str]], present: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Map current symbols onto the collapsed U133A matrix. Honest intersection."""
    mapped: dict[str, list[str]] = {}
    rows = []
    for term, members in raw.items():
        seen = []
        for g in members:
            hit = g if g in present else U133A_ALIASES.get(g)
            if hit in present and hit not in seen:
                seen.append(hit)
                rows.append(
                    {
                        "term": term,
                        "query": g,
                        "matrix_symbol": hit,
                        "status": "direct" if g in present else "alias",
                    }
                )
            else:
                rows.append(
                    {
                        "term": term,
                        "query": g,
                        "matrix_symbol": "",
                        "status": "absent_on_U133A",
                    }
                )
        mapped[term] = seen
    return mapped, pd.DataFrame(rows)


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
    ax.set_title(f"GSE68465 CLDN4 Q4 vs Q1 prerank GSEA (n={n_q4} vs {n_q1})")
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
    n_tumor: int,
    n_matrix: int,
    n_normal: int,
    n_genes: int,
    cldn4_t: float,
    q_cuts: tuple[float, float],
    coverage: pd.DataFrame,
    used_probe: str,
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

    cov_bits = []
    for term in HEADLINE + SECONDARY:
        sub = coverage[coverage["term"] == term]
        n_dir = int((sub["status"] == "direct").sum())
        n_alias = int((sub["status"] == "alias").sum())
        n_abs = int((sub["status"] == "absent_on_U133A").sum())
        cov_bits.append(
            f"{HEADLINE_LABEL[term]} {n_dir + n_alias}/{len(sub)} on U133A "
            f"({n_dir} direct, {n_alias} alias, {n_abs} absent)"
        )
    cov_txt = "; ".join(cov_bits)

    md = f"""# FINDING — GSE68465 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive only. CLDN4 only.** Public Director's Challenge LUAD microarray (Shedden et al., *Nat Med* 2008, PMID [18641660](https://pubmed.ncbi.nlm.nih.gov/18641660/); GEO [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465)). This folder does **not** audit or retract `methods/gse68465_cldn4` (CLDN4 vs CD8A after ESTIMATE purity). TACSTD2 is **not** an anchor here. No slide was re-scored. Unit is the **tumour array**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

GEO deposits **462** arrays (**443** `Lung Adenocarcinoma` + **19** `Normal`). The paper text says 442 LUAD. Honest tumour n here is the GEO filter: **443**. The GSEA contrast is the quartile arms (**{n_q4} vs {n_q1}**), not n=443 and not n=462. No ICI labels. This is bulk U133A — an IFN / MHC NES can be infiltrate, not a tumour-cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | Honest n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE68465 CLDN4 Q4 vs Q1 | **{n_q4} vs {n_q1}** of 443 LUAD | Welch *t* on log2(MAS5) | {call_headline(gsea, "welch_t_Q4_minus_Q1")} |

Primary rank is Welch *t* (Q4 − Q1) on log2(MAS5). CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: {tj_drop_txt}. U133A coverage is incomplete; n in rank is the intersection.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL / MAS5 reprocess. |
| File | `GSE68465_series_matrix.txt.gz` |
| Platform | GPL96 Affymetrix HG-U133A |
| Arrays | **n = {n_matrix}** (22,283 probes) |
| Tumours used | **n = {n_tumor}** `disease_state = Lung Adenocarcinoma` |
| Held out | **{n_normal}** `Normal` (Stratagene) |
| Paper vs GEO | Shedden text = 442 LUAD; GEO filter = **443**. This page uses 443. |
| Genes after max-mean | {n_genes} (first `///` HUGO from official GPL96) |
| Transform | `log2(MAS5)` (deposited values are MAS5 intensity, min>0) |
| Anchor | **CLDN4 only** (`{used_probe}`; named probe `{NAMED_CLDN4}`) |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| Honest GSEA n | Q4 = {n_q4}, Q1 = {n_q1} (not 443, not 462) |
| CLDN4 log2 cuts | P25 = {q_cuts[0]:.4f}; P75 = {q_cuts[1]:.4f} |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all {n_tumor} LUAD); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| U133A aliases | `WARS1→WARS`, `MARCHF1→MARCH1`, `PALS1→MPP5` (documented HGNC updates only) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; ICI response; CEL reprocess; sample-permutation GSEA; CD8 / ESTIMATE re-run |

Series: [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465). Title on GEO: *caArray_jacob-00182: gene expression–based survival prediction in lung adenocarcinoma*. Paper: Shedden et al., *Nat Med* 2008; PMID [18641660](https://pubmed.ncbi.nlm.nih.gov/18641660/).

U133A intersection (after aliases): {cov_txt}. `CD274` and `NLRC5` are **absent** on this platform.

CLDN4 itself is a KEGG TJ member and is the split gene. Its Welch *t* on this contrast is **{cldn4_t:+.2f}**. TJ NES on the full rank is therefore partly circular. The drop-CLDN4 row is the non-circular companion.

---

## Headline NES (Welch *t*, Q4 − Q1)

{n_ranked} genes ranked. **{n_q4} Q4 vs {n_q1} Q1** of 443 LUAD. CLDN4 itself is the stratification gene (Welch *t* = {cldn4_t:+.3f}; by construction the top of the rank).

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

Uses all **n = {n_tumor}** LUAD arrays. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets. This is **not** the GSEA n — the quartile contrast remains {n_q4} vs {n_q1}.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "spearman_vs_CLDN4")}

---

## Extra figures

- `{figures[0]}`

Tables: `methods/gse68465_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `geneset_coverage.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run or retraction of `methods/gse68465_cldn4` (CD8A / ESTIMATE purity).
- Not n=462 and not the paper's 442. Tumour n = 443; GSEA n = {n_q4} vs {n_q1}.
- Not CEL / RMA / fRMA from RAW.tar. The deposited MAS5 series matrix is used as deposited.
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN/MHC call. This is bulk U133A surgical LUAD.
- Not an ICI / RECIST / PD-L1 analysis. This series has no ICI labels.
- Not a claim that CLDN4 *induces* IFN or MHC. Association only.

---

## 中文摘要

只补公开 **GSE68465** Director's Challenge LUAD 芯片上 **CLDN4 Q4 vs Q1** 的 prerank GSEA（IFN / MHC / TJ），不审不撤已有 `gse68465_cldn4`（CD8 / 纯度）页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- GEO 462 张芯片：443 LUAD + 19 Normal。论文写 442；本页用 GEO 过滤 **n=443**。GSEA 对比是四分位臂 **{n_q4} vs {n_q1}**，不要把 443 写成 GSEA 的 n。
- Headline（Welch *t*，log2 MAS5）：{call_headline(gsea, "welch_t_Q4_minus_Q1")}。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：{tj_drop_txt}。
- U133A 覆盖不全（{cov_txt}）。`CD274` / `NLRC5` 不在此平台。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    if not MATRIX.exists() or not ANNOT.exists():
        raise SystemExit("missing GEO files; run download.py first")

    raw_sets = load_sets_raw()
    for t, genes in raw_sets.items():
        log(f"set {t}: {len(genes)} current symbols")

    gpl = load_gpl_annot(ANNOT)
    id_col = "ID" if "ID" in gpl.columns else gpl.columns[0]
    sym_col = "Gene symbol"
    p2g = gpl.set_index(id_col)[sym_col].astype(str)
    p2g = p2g.replace({"nan": np.nan, "None": np.nan, "": np.nan}).dropna()
    p2g = p2g.map(lambda x: str(x).split("///")[0].strip())
    p2g = p2g[p2g.str.len() > 0]

    log("parsing GSE68465 series matrix")
    meta, probes = parse_geo_matrix(MATRIX)
    n_matrix = int(meta.shape[0])
    log(f"matrix samples={n_matrix} probes={probes.shape[0]}")
    disease = meta["disease_state"].astype(str)
    keep = disease.str.contains("Adenocarcinoma", case=False, na=False)
    n_tumor = int(keep.sum())
    n_normal = int((~keep).sum())
    log(f"LUAD n={n_tumor}; held-out Normal n={n_normal}")
    if n_tumor != 443:
        raise SystemExit(f"expected 443 LUAD arrays, got {n_tumor}")
    if NAMED_CLDN4 not in probes.index:
        raise SystemExit(f"named CLDN4 probe {NAMED_CLDN4} missing")

    probes_t = probes.loc[:, meta.index[keep]]
    genes = collapse_maxmean(probes_t, p2g)
    log(f"genes after max-mean collapse: {genes.shape[0]}")
    if ANCHOR not in genes.index:
        raise SystemExit(f"{ANCHOR} missing after collapse")

    used_probe = ""
    gvals = genes.loc[ANCHOR]
    cands = p2g[p2g == ANCHOR].index.intersection(probes_t.index)
    for pr in cands:
        if np.allclose(probes_t.loc[pr].to_numpy(float), gvals.to_numpy(float), equal_nan=True):
            used_probe = str(pr)
            break
    log(f"CLDN4 max-mean probe={used_probe} named={NAMED_CLDN4}")

    mapped, coverage = remap_sets(raw_sets, set(genes.index.astype(str)))
    for t, members in mapped.items():
        log(f"U133A {t}: {len(members)} / {len(raw_sets[t])}")

    # MAS5 intensity → log2. All values > 0 on this deposit.
    if float(np.nanmin(genes.to_numpy())) <= 0:
        raise SystemExit("unexpected non-positive MAS5; refuse silent clip")
    logx = np.log2(genes.astype(float))
    cldn4 = logx.loc[ANCHOR]
    q1, q3 = float(cldn4.quantile(0.25)), float(cldn4.quantile(0.75))
    high, low = quartile_split(cldn4)
    n_q1, n_q4 = int(len(low)), int(len(high))
    log(f"CLDN4 P25={q1:.4f} P75={q3:.4f}; Q1={n_q1} Q4={n_q4} of {n_tumor}")
    if n_q4 + n_q1 >= n_tumor:
        raise SystemExit("quartile tails ate the middle; check ties")

    split = pd.DataFrame(
        {
            "sample": cldn4.index,
            "CLDN4_log2_MAS5": cldn4.to_numpy(),
            "quartile": [
                "Q1" if s in set(low) else ("Q4" if s in set(high) else "Q2Q3")
                for s in cldn4.index
            ],
            "site": meta.loc[cldn4.index, "source"].to_numpy() if "source" in meta.columns else "",
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

    g_t = run_rank(rank_t, mapped, "welch_t_Q4_minus_Q1", n_q4, n_q1, n_tumor)
    g_drop = run_rank(rank_t_drop, mapped, "welch_t_Q4_minus_Q1_drop_CLDN4", n_q4, n_q1, n_tumor)
    g_rho = run_rank(rank_rho, mapped, "spearman_vs_CLDN4", n_q4, n_q1, n_tumor)
    gsea = pd.concat([g_t, g_drop, g_rho], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "surgical LUAD bulk U133A (not cell-intrinsic; no ICI labels)",
                "anchor": ANCHOR,
                "n_matrix": n_matrix,
                "n_normal_held_out": n_normal,
                "n_tumor": n_tumor,
                "n_q4": n_q4,
                "n_q1": n_q1,
                "paper_n_luad": 442,
                "honest_gsea_n": f"{n_q4} vs {n_q1}",
                "cldn4_p25": q1,
                "cldn4_p75": q3,
                "cldn4_welch_t": cldn4_t,
                "cldn4_maxmean_probe": used_probe,
                "cldn4_named_probe": NAMED_CLDN4,
                "rank_metric_primary": "Welch t (CLDN4 Q4 minus Q1) on log2(MAS5)",
                "n_genes_collapsed": int(genes.shape[0]),
                "n_genes_ranked": int(rank_t.shape[0]),
                "matrix_file": MATRIX.name,
                "series_title": "Director's Challenge LUAD (Shedden / jacob-00182)",
                "pmid": "18641660",
                "u133a_aliases": "WARS1->WARS; MARCHF1->MARCH1; PALS1->MPP5",
                "not": "TACSTD2; CEL reprocess; ICI labels; sample-permutation GSEA; CD8/ESTIMATE re-run",
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    split.to_csv(TAB / "quartile_split.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    coverage.to_csv(TAB / "geneset_coverage.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline, n_q4, n_q1)
    write_finding(
        gsea,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        n_tumor,
        n_matrix,
        n_normal,
        int(genes.shape[0]),
        cldn4_t,
        (q1, q3),
        coverage,
        used_probe,
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
    print(f"\nhonest n: {n_q4} vs {n_q1} of {n_tumor} LUAD (matrix {n_matrix}, normals {n_normal})")


if __name__ == "__main__":
    main()
