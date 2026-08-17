#!/usr/bin/env python3
"""ADDITIVE leftover — GSE253564 CLDN4-only Q4 vs Q1 GSEA (IFN / MHC / TJ) and vs CD274.

Public pre-treatment FPKM only. No FASTQ / SRA. Does not re-run the TACSTD2
keratin/TJ leftover in a8_ici_gsea or the MPR/PFS leftover in opus_geo_leftover.

Engine matches methods/cldn4_ko_gsea (gsea_core.py): weighted KS p=1, 1000
gene-set permutations, seed=42. Positive NES = enriched in CLDN4 Q4.
"""
from __future__ import annotations

import hashlib
import json
import time
import urllib.request
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
DATA.mkdir(parents=True, exist_ok=True)

ACCESSION = "GSE253564"
CONTRAST = "GSE253564_CLDN4_Q4_vs_Q1"
FPKM_NAME = "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/"
    "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
)

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

MHC_GENES = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (p != p)):
        return "NA"
    if abs(p) < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def fmt_rho(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:+.3f}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_fpkm() -> Path:
    dest = DATA / FPKM_NAME
    if dest.exists() and dest.stat().st_size > 1000:
        log(f"using cached {dest.name} ({dest.stat().st_size} bytes)")
        return dest
    log(f"download {FPKM_URL}")
    req = urllib.request.Request(FPKM_URL, headers={"User-Agent": "gse253564-cldn4-gsea/1.0"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as fh, tmp.open("wb") as out:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.rename(dest)
    log(f"saved {dest.name} ({dest.stat().st_size} bytes) sha256={sha256(dest)}")
    return dest


def collapse_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    expr.index = expr.index.str.strip().str.strip('"')
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.index.duplicated().any():
        means = expr.mean(axis=1)
        keep_idx = means.groupby(level=0).idxmax()
        expr = expr.loc[keep_idx]
    return expr.loc[~expr.index.duplicated(keep="first")]


def load_expr(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t")
    gene_col = raw.columns[0]
    drop = [c for c in raw.columns if c.lower().replace(".", "_") in {"entrez_id", "entrez"}]
    expr = raw.set_index(gene_col).drop(columns=drop, errors="ignore")
    expr = collapse_symbols(expr)
    expr = np.log2(expr.clip(lower=0) + 1.0)
    return expr


def mean_z(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns)
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def spearman_pair(x: pd.Series, y: pd.Series) -> dict:
    aligned = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(aligned))
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(aligned["x"], aligned["y"])
    return {"n": n, "rho": float(rho), "p": float(p)}


def mwu_q4_q1(values: pd.Series, high: pd.Index, low: pd.Index) -> dict:
    a = values.reindex(high).dropna()
    b = values.reindex(low).dropna()
    if len(a) < 3 or len(b) < 3:
        return {
            "n_q4": int(len(a)),
            "n_q1": int(len(b)),
            "median_q4": float(a.median()) if len(a) else np.nan,
            "median_q1": float(b.median()) if len(b) else np.nan,
            "p": np.nan,
        }
    u, p = stats.mannwhitneyu(a.to_numpy(), b.to_numpy(), alternative="two-sided")
    return {
        "n_q4": int(len(a)),
        "n_q1": int(len(b)),
        "median_q4": float(a.median()),
        "median_q1": float(b.median()),
        "U": float(u),
        "p": float(p),
    }


def load_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in ALL_SETS}


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


def write_figures(
    headline: pd.DataFrame,
    cldn4: pd.Series,
    cd274: pd.Series,
    high: pd.Index,
    low: pd.Index,
) -> list[str]:
    h = headline[headline["rank_metric"] == "welch_t_Q4_vs_Q1"].copy()
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
    ax.set_title(f"GSE253564 leftover — CLDN4 Q4 vs Q1 prerank GSEA (n={len(high)} vs {len(low)})")
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

    aligned = pd.concat([cldn4.rename("CLDN4"), cd274.rename("CD274")], axis=1).dropna()
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    axes[0].scatter(aligned["CLDN4"], aligned["CD274"], s=28, c="#2c3e50", edgecolors="white", linewidths=0.4)
    if len(aligned) >= 3:
        slope, intercept, *_ = stats.linregress(aligned["CLDN4"], aligned["CD274"])
        xs = np.linspace(aligned["CLDN4"].min(), aligned["CLDN4"].max(), 50)
        axes[0].plot(xs, slope * xs + intercept, color="#c0392b", lw=1.0)
    rho_row = spearman_pair(cldn4, cd274)
    axes[0].set_xlabel("CLDN4 log2(FPKM+1)")
    axes[0].set_ylabel("CD274 log2(FPKM+1)")
    axes[0].set_title(f"CLDN4 vs CD274  ρ={fmt_rho(rho_row['rho'])}  p={fmt_p(rho_row['p'])}  n={rho_row['n']}")

    axes[1].barh(y, h["nes"], color=colors, edgecolor="black", linewidth=0.4)
    axes[1].axvline(0, color="black", lw=0.7)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(h["label"])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("NES")
    axes[1].set_title("Headline NES (Welch t, Q4 vs Q1)")
    fig.tight_layout()
    png2 = FIG / "fig_cd274_and_nes.png"
    pdf2 = FIG / "fig_cd274_and_nes.pdf"
    fig.savefig(png2, dpi=160)
    fig.savefig(pdf2)
    plt.close(fig)
    return [str(png.relative_to(HERE.parent.parent)), str(png2.relative_to(HERE.parent.parent))]


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


def call_headline(gsea: pd.DataFrame, rank_metric: str) -> str:
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(HEADLINE))]
    bits = []
    for term in HEADLINE:
        r = sub[sub["term"] == term].iloc[0]
        bits.append(f"{HEADLINE_LABEL[term]} NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])}")
    return "; ".join(bits)


def write_finding(
    gsea: pd.DataFrame,
    corr: pd.DataFrame,
    inventory: dict,
    figures: list[str],
) -> None:
    n_total = int(inventory["n_total"])
    n_q4 = int(inventory["n_q4"])
    n_q1 = int(inventory["n_q1"])
    n_ranked = int(inventory["n_genes_ranked"])
    q1_cut = float(inventory["q1_cut"])
    q3_cut = float(inventory["q3_cut"])
    cldn4_t = float(inventory["cldn4_welch_t"])
    cldn4_lfc = float(inventory["cldn4_log2fc_q4_minus_q1"])

    cd = corr[corr["pair"] == "CLDN4 vs CD274"].iloc[0]
    ifn = corr[corr["pair"] == "CLDN4 vs Hallmark IFN-γ score"].iloc[0]
    mhc = corr[corr["pair"] == "CLDN4 vs MHC-I score"].iloc[0]
    cd_q = corr[corr["pair"] == "CD274 Q4 vs Q1 MWU"].iloc[0]

    ifna = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_vs_Q1")
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
        r = gsea[(gsea["rank_metric"] == "welch_t_Q4_vs_Q1") & (gsea["term"] == term)].iloc[0]
        lead[term] = r["lead_genes"]

    drop_txt = "not run"
    drop = gsea[gsea["rank_metric"] == "welch_t_Q4_vs_Q1_drop_CLDN4"]
    if not drop.empty:
        drop_txt = call_headline(gsea, "welch_t_Q4_vs_Q1_drop_CLDN4")

    spear_txt = "not run"
    spear = gsea[gsea["rank_metric"] == "spearman_vs_CLDN4"]
    if not spear.empty:
        spear_txt = call_headline(gsea, "spearman_vs_CLDN4")

    md = f"""# FINDING — GSE253564 leftover CLDN4 Q4 vs Q1 GSEA (IFN / MHC / TJ) and vs CD274

**Additive only.** Public pre-treatment FPKM from the leftover neoadjuvant **durvalumab ± SBRT** series. This folder does **not** audit or retract the TACSTD2 keratin/TJ leftover (`a8_ici_gsea`) or the MPR/PFS leftover (`opus_geo_leftover`). **CLDN4 only** — TACSTD2 is not used to split samples.

Prerank engine is the same as `methods/cldn4_ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (high minus low). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv` and `tables/correlations.tsv`.

---

## 一句话 / TL;DR

| Contrast | Honest n | CLDN4 Q4 vs Q1 | IFN-γ / MHC-I / KEGG TJ NES (FDR) | CLDN4 vs CD274 |
|---|---|---|---|---|
| GSE253564 pre-treatment leftover | **{n_q4} vs {n_q1}** of {n_total} | Welch *t* {cldn4_t:+.2f}; mean log2FC {cldn4_lfc:+.3f} | {call_headline(gsea, "welch_t_Q4_vs_Q1")} | ρ = {fmt_rho(cd["rho"])} p = {fmt_p(cd["p"])} (n={int(cd["n"])}) |

GEO deposits **{n_total}** pre-treatment tumours. The GSEA contrast is the quartile arms (**{n_q4} vs {n_q1}**), not n={n_total}. The CD274 Spearman uses every sample with both genes (n={int(cd["n"])}).

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author FPKM only. No FASTQ / SRA. |
| File | `{FPKM_NAME}` |
| Scale | log2(FPKM+1); duplicate symbols collapsed to the highest-mean locus |
| Split | CLDN4 quartiles: Q1 ≤ {q1_cut:.3f}, Q4 ≥ {q3_cut:.3f} (log2(FPKM+1)) |
| Honest n | GEO n={n_total}; Q4 vs Q1 = **{n_q4} vs {n_q1}** |
| Sign | positive = up in CLDN4 Q4 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (n={n_total}); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| CD274 | Spearman on all complete pairs; MWU of CD274 in Q4 vs Q1 |
| Out of scope | TACSTD2 splits; GSE248378 post-treatment; MPR/PFS re-test; Hallmark EMT / keratin |

Series: [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564). Title on GEO: *Features associated with enhanced proliferation define lung tumors that respond to dual immune checkpoint and radiation therapy*. Platform GPL24676. PMID 38401548.

CLDN4 is a **KEGG tight-junction member**. Q4 vs Q1 is defined on CLDN4, so the TJ NES is partly circular. The drop-CLDN4 sensitivity is the non-circular TJ number.

---

## Headline NES (Welch *t*, Q4 vs Q1)

{n_ranked} genes ranked. CLDN4 itself: Welch *t* **{cldn4_t:+.2f}**, mean log2(FPKM+1) difference Q4−Q1 **{cldn4_lfc:+.3f}** (the split worked).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_vs_Q1")}

Hallmark IFN-α (secondary, not in BH): {ifna_txt}.

Leading-edge (first 25, Welch rank):

- IFN-γ: `{lead["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}`
- MHC-I / APM: `{lead["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}`
- KEGG TJ: `{lead["KEGG_TIGHT_JUNCTION"]}`

---

## Sensitivity

Same three headline sets.

| Rank | n used | Result |
|---|---|---|
| Spearman ρ vs continuous CLDN4 | {n_total} | {spear_txt} |
| Welch *t* after dropping CLDN4 from the rank | {n_q4} vs {n_q1} | {drop_txt} |

---

## CLDN4 vs CD274 (and companion scores)

Spearman, two-sided, complete cases. IFN-γ and MHC-I scores are the mean of within-cohort gene-wise z-scores (Hallmark IFN-γ genes present; MHC-I = HLA-A/B/C + B2M).

| Pair | n | ρ | p | Note |
|---|---:|---:|---:|---|
| CLDN4 vs CD274 | {int(cd["n"])} | {fmt_rho(cd["rho"])} | {fmt_p(cd["p"])} | primary leftover correlation |
| CLDN4 vs Hallmark IFN-γ score | {int(ifn["n"])} | {fmt_rho(ifn["rho"])} | {fmt_p(ifn["p"])} | companion, not a GSEA substitute |
| CLDN4 vs MHC-I score (HLA-A/B/C + B2M) | {int(mhc["n"])} | {fmt_rho(mhc["rho"])} | {fmt_p(mhc["p"])} | companion |
| CD274 Q4 vs Q1 (MWU) | {int(cd_q["n_q4"])} vs {int(cd_q["n_q1"])} | med {cd_q["median_q4"]:.3f} vs {cd_q["median_q1"]:.3f} | {fmt_p(cd_q["p"])} | same quartile arms as GSEA |

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse253564_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `correlations.tsv`, `inventory.tsv`.

---

## What this is not

- Not a TACSTD2 analysis and not a dual-high gate.
- Not a re-run of A8 leftover keratin/TJ/EMT (that split was TACSTD2).
- Not an MPR / PFS / recurrence test (already in `opus_geo_leftover`; CLDN4 vs MPR was null).
- Not GSE248378 (post-treatment residual tumours; 0 MPR in that matrix).
- Not sample-permutation GSEA. Quartile n = {n_q4} vs {n_q1}; the engine is gene-set permutation on a prerank.
- Not a claim that CLDN4 *induces* IFN or PD-L1. Association only. Bulk FFPE mixes epithelium and infiltrate.

---

## 中文摘要

只补公开 **GSE253564** 治疗前 leftover durvalumab ± SBRT bulk 的 **CLDN4 Q4 vs Q1** prerank GSEA（IFN / MHC / TJ）以及 CLDN4 vs CD274，不审不撤已有页。不用 TACSTD2 分层。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- GEO n={n_total}。四分位对比是 **{n_q4} vs {n_q1}**，不要把 n={n_total} 写成 GSEA 的 n。
- Headline（Welch *t*）：{call_headline(gsea, "welch_t_Q4_vs_Q1")}。
- CLDN4 vs CD274：ρ = {fmt_rho(cd["rho"])}，p = {fmt_p(cd["p"])}，n={int(cd["n"])}。
- CLDN4 在 KEGG TJ 里；去掉 CLDN4 后再跑的 TJ 才是非循环数字。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    path = download_fpkm()
    expr = load_expr(path)
    if "CLDN4" not in expr.index:
        raise SystemExit("CLDN4 absent from GSE253564 FPKM")
    if "CD274" not in expr.index:
        raise SystemExit("CD274 absent from GSE253564 FPKM")

    cldn4 = expr.loc["CLDN4"]
    cd274 = expr.loc["CD274"]
    high, low = quartile_split(cldn4)
    n_total = int(cldn4.notna().sum())
    n_q4, n_q1 = int(len(high)), int(len(low))
    q1_cut = float(cldn4.quantile(0.25))
    q3_cut = float(cldn4.quantile(0.75))
    log(f"n_total={n_total} Q4={n_q4} Q1={n_q1} cuts=({q1_cut:.3f}, {q3_cut:.3f})")
    if n_q4 < 6 or n_q1 < 6:
        raise SystemExit(f"quartile arms too small: {n_q4} vs {n_q1}")

    sets = load_sets()
    rank_t = rank_high_vs_low(expr, high, low)
    rank_t_drop = rank_t.drop(labels=["CLDN4"], errors="ignore")
    rank_rho = rank_spearman_vs_target(expr, cldn4)

    gsea_parts = [
        run_rank(rank_t, sets, "welch_t_Q4_vs_Q1", n_q4, n_q1, n_total),
        run_rank(rank_t_drop, sets, "welch_t_Q4_vs_Q1_drop_CLDN4", n_q4, n_q1, n_total),
        run_rank(rank_rho, sets, "spearman_vs_CLDN4", n_total, n_total, n_total),
    ]
    gsea = pd.concat(gsea_parts, ignore_index=True)

    ifng_genes = sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]
    ifng_score = mean_z(expr, ifng_genes)
    mhc_score = mean_z(expr, MHC_GENES)

    pairs = []
    for name, y in [
        ("CLDN4 vs CD274", cd274),
        ("CLDN4 vs Hallmark IFN-γ score", ifng_score),
        ("CLDN4 vs MHC-I score", mhc_score),
    ]:
        row = spearman_pair(cldn4, y)
        row.update({"pair": name, "test": "spearman"})
        pairs.append(row)
    mwu = mwu_q4_q1(cd274, high, low)
    mwu.update({"pair": "CD274 Q4 vs Q1 MWU", "test": "mannwhitneyu", "n": mwu["n_q4"] + mwu["n_q1"], "rho": np.nan})
    corr = pd.DataFrame(pairs + [mwu])

    cldn4_t = float(rank_t.loc["CLDN4"]) if "CLDN4" in rank_t.index else np.nan
    cldn4_lfc = float(expr.loc["CLDN4", high].mean() - expr.loc["CLDN4", low].mean())

    inventory = {
        "accession": ACCESSION,
        "file": FPKM_NAME,
        "sha256": sha256(path),
        "bytes": int(path.stat().st_size),
        "url": FPKM_URL,
        "n_total": n_total,
        "n_q4": n_q4,
        "n_q1": n_q1,
        "q1_cut": q1_cut,
        "q3_cut": q3_cut,
        "n_genes_matrix": int(expr.shape[0]),
        "n_genes_ranked": int(rank_t.shape[0]),
        "cldn4_present": True,
        "cd274_present": True,
        "cldn4_welch_t": cldn4_t,
        "cldn4_log2fc_q4_minus_q1": cldn4_lfc,
        "engine": f"gsea_core.py weighted KS p=1 nperm={NPERM} seed={SEED}",
        "note": "CLDN4 only; leftover Durva pre-treatment bulk; Q4 vs Q1 GSEA + CD274 Spearman",
    }

    headline = gsea[gsea["term"].isin(HEADLINE)].copy()
    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    corr.to_csv(TAB / "correlations.tsv", sep="\t", index=False)
    pd.DataFrame([inventory]).to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    (TAB / "summary.json").write_text(json.dumps(inventory, indent=2) + "\n")

    sample_tbl = pd.DataFrame(
        {
            "sample": cldn4.index,
            "CLDN4_log2fpkm1": cldn4.to_numpy(),
            "CD274_log2fpkm1": cd274.reindex(cldn4.index).to_numpy(),
            "IFN_score": ifng_score.reindex(cldn4.index).to_numpy(),
            "MHC_score": mhc_score.reindex(cldn4.index).to_numpy(),
            "quartile": [
                "Q4" if s in set(high) else ("Q1" if s in set(low) else "Q2Q3") for s in cldn4.index
            ],
        }
    )
    sample_tbl.to_csv(TAB / "sample_table.tsv", sep="\t", index=False)

    figures = write_figures(headline, cldn4, cd274, high, low)
    write_finding(gsea, corr, inventory, figures)
    log("done")
    print(headline[headline["rank_metric"] == "welch_t_Q4_vs_Q1"][["term", "nes", "fdr", "nom_p", "n_set_in_rank"]].to_string(index=False))
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
