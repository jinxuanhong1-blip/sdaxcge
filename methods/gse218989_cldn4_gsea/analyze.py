#!/usr/bin/env python3
"""ADDITIVE extra — GSE218989 CLDN4 Q4 vs Q1 prerank GSEA.

Public SMC-KAIST NSCLC ICI bulk TPM (Kang et al., Nat Commun 2024).
Headline sets: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.

Engine matches methods/cldn4_ko_gsea (gsea_core.py): weighted KS p=1,
1000 gene-set permutations, seed=42. Positive NES = enriched in CLDN4 Q4.

Response OR is already known NS from methods/gse218989_cldn4_ici
(GEO MWU p=0.40, AUC=0.47). It is recorded in a note table only —
do not headline it.
"""
from __future__ import annotations

import gzip
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
CACHE = Path("/tmp/gse218989")
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

TPM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/suppl/"
    "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218989/matrix/"
    "GSE218989_series_matrix.txt.gz"
)
TPM_FILE = "GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"
MATRIX_FILE = "GSE218989_series_matrix.txt.gz"

CONTRAST = "GSE218989_CLDN4_Q4_vs_Q1"
ACCESSION = "GSE218989"

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
    "KEGG_TIGHT_JUNCTION_NO_CLDN4": "KEGG TJ (no CLDN4)",
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


def fmt_or(x) -> str:
    if x is None or (isinstance(x, float) and (x != x)):
        return "NA"
    return f"{x:.2f}"


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    log(f"download {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def load_tpm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df.index = df.index.astype(str)
    df = df[~df.index.duplicated(keep="first")]
    return df.astype(float)


def load_geo_meta(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    titles = None
    cols: dict[str, list[str]] = {}
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            cols["gsm"] = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in line.split("\t")[1:]]
            key = vals[0].split(":", 1)[0].strip()
            cols[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return pd.DataFrame({"patient": titles, **cols})


def load_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in ALL_SETS}


def fisher_or(a: int, b: int, c: int, d: int) -> dict:
    """OR for 2x2 [[a,b],[c,d]] = Q4 R/NR vs Q1 R/NR. Haldane-Anscombe if a zero cell."""
    table = np.array([[a, b], [c, d]], dtype=float)
    if (table == 0).any():
        table = table + 0.5
    or_hat = float((table[0, 0] * table[1, 1]) / (table[0, 1] * table[1, 0]))
    se = float(np.sqrt(1 / table[0, 0] + 1 / table[0, 1] + 1 / table[1, 0] + 1 / table[1, 1]))
    lo = float(np.exp(np.log(or_hat) - 1.96 * se))
    hi = float(np.exp(np.log(or_hat) + 1.96 * se))
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return {"or": or_hat, "or_lo": lo, "or_hi": hi, "p_fisher": float(p)}


def run_rank(
    rank: pd.Series,
    sets: dict[str, list[str]],
    rank_metric: str,
    n_q4: int,
    n_q1: int,
    headline_terms: list[str] | None = None,
) -> pd.DataFrame:
    headline_terms = list(headline_terms) if headline_terms is not None else list(HEADLINE)
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = CONTRAST
    g["accession"] = ACCESSION
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_q4"] = n_q4
    g["n_q1"] = n_q1
    g["headline"] = g["term"].isin(headline_terms)
    g["fdr"] = np.nan
    head = g["term"].isin(headline_terms)
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    g["direction"] = np.where(g["nes"] > 0, "UP_in_CLDN4_Q4", "DOWN_in_CLDN4_Q4")
    return g


def write_figures(headline: pd.DataFrame) -> list[str]:
    h = headline[headline["rank_metric"] == "welch_t_Q4_minus_Q1"].copy()
    h = h[h["term"].isin(HEADLINE)]
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
    ax.set_xlabel("NES (positive = enriched in CLDN4 Q4)")
    ax.set_title("GSE218989 CLDN4 Q4 vs Q1 — headline prerank GSEA")
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

    # two-panel: primary NES + Spearman sensitivity
    s = headline[headline["rank_metric"] == "spearman_vs_CLDN4"].copy()
    s = s[s["term"].isin(HEADLINE)]
    s["label"] = s["term"].map(HEADLINE_LABEL)
    s = s.set_index("label").reindex(order).reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5), sharey=True)
    for ax, df, title in (
        (axes[0], h, "Primary: Welch t Q4−Q1"),
        (axes[1], s, "Sensitivity: Spearman vs CLDN4"),
    ):
        cols = ["#c0392b" if n > 0 else "#2980b9" for n in df["nes"]]
        ax.barh(y, df["nes"], color=cols, edgecolor="black", linewidth=0.4)
        ax.axvline(0, color="black", lw=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(order)
        ax.invert_yaxis()
        ax.set_xlabel("NES")
        ax.set_title(title)
        for i, row in df.iterrows():
            ax.text(
                row["nes"] + (0.03 if row["nes"] >= 0 else -0.03),
                i,
                f"{row['nes']:+.2f}",
                va="center",
                ha="left" if row["nes"] >= 0 else "right",
                fontsize=8,
            )
    fig.tight_layout()
    png2 = FIG / "fig_nes_primary_vs_spearman.png"
    pdf2 = FIG / "fig_nes_primary_vs_spearman.pdf"
    fig.savefig(png2, dpi=160)
    fig.savefig(pdf2)
    plt.close(fig)
    return [
        str(png.relative_to(HERE.parent.parent)),
        str(png2.relative_to(HERE.parent.parent)),
    ]


def nes_md_rows(gsea: pd.DataFrame, rank_metric: str, terms=None) -> str:
    terms = list(terms) if terms is not None else list(HEADLINE)
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(terms))].copy()
    sub["ord"] = sub["term"].map({t: i for i, t in enumerate(terms)})
    sub = sub.sort_values("ord")
    lines = []
    for _, r in sub.iterrows():
        label = HEADLINE_LABEL.get(r["term"], r["term"])
        fdr = r["fdr"] if pd.notna(r["fdr"]) else r["nom_p"]
        lines.append(
            f"| {label} | `{r['term']}` | {fmt_nes(r['nes'])} | "
            f"{fmt_p(fdr)} | {fmt_p(r['nom_p'])} | {int(r['n_set_in_rank'])} |"
        )
    return "\n".join(lines)


def call_headline(gsea: pd.DataFrame, rank_metric: str, terms=None) -> str:
    terms = list(terms) if terms is not None else list(HEADLINE)
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(terms))]
    bits = []
    for term in terms:
        hit = sub[sub["term"] == term]
        if hit.empty:
            bits.append(f"{HEADLINE_LABEL.get(term, term)} NES NA")
            continue
        r = hit.iloc[0]
        fdr = r["fdr"] if pd.notna(r["fdr"]) else r["nom_p"]
        bits.append(f"{HEADLINE_LABEL[term]} NES {fmt_nes(r['nes'])} FDR {fmt_p(fdr)}")
    return "; ".join(bits)


def honest_read(gsea: pd.DataFrame) -> str:
    """One-sentence honest read of the three headline NES. No spin."""
    sub = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_minus_Q1") & (gsea["term"].isin(HEADLINE))
    ]
    bits = []
    for term in HEADLINE:
        r = sub[sub["term"] == term].iloc[0]
        sig = pd.notna(r["fdr"]) and r["fdr"] < 0.05
        if not sig:
            bits.append(f"{HEADLINE_LABEL[term]} is NS (NES {fmt_nes(r['nes'])}, FDR {fmt_p(r['fdr'])})")
        elif r["nes"] > 0:
            bits.append(f"{HEADLINE_LABEL[term]} is up in Q4 (NES {fmt_nes(r['nes'])}, FDR {fmt_p(r['fdr'])})")
        else:
            bits.append(f"{HEADLINE_LABEL[term]} is down in Q4 (NES {fmt_nes(r['nes'])}, FDR {fmt_p(r['fdr'])})")
    return "; ".join(bits) + "."


def write_finding(
    gsea: pd.DataFrame,
    figures: list[str],
    n_ranked: int,
    n_q4: int,
    n_q1: int,
    n_genes: int,
    n_patients: int,
    q4_cut: float,
    q1_cut: float,
    or_row: dict,
    cldn4_t: float,
) -> None:
    lead = {}
    for term in HEADLINE:
        r = gsea[
            (gsea["rank_metric"] == "welch_t_Q4_minus_Q1") & (gsea["term"] == term)
        ].iloc[0]
        lead[term] = r["lead_genes"]

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

    tj_drop = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_minus_Q1_drop_CLDN4")
        & (gsea["term"] == "KEGG_TIGHT_JUNCTION")
    ]
    tj_drop_txt = "not ranked"
    if not tj_drop.empty:
        r = tj_drop.iloc[0]
        tj_drop_txt = (
            f"NES {fmt_nes(r['nes'])} nom p {fmt_p(r['nom_p'])} "
            f"(n={int(r['n_set_in_rank'])}; not in headline FDR)"
        )

    tj_noc = gsea[
        (gsea["rank_metric"] == "welch_t_Q4_minus_Q1")
        & (gsea["term"] == "KEGG_TIGHT_JUNCTION_NO_CLDN4")
    ]
    tj_noc_txt = "not ranked"
    if not tj_noc.empty:
        r = tj_noc.iloc[0]
        tj_noc_txt = (
            f"NES {fmt_nes(r['nes'])} nom p {fmt_p(r['nom_p'])} "
            f"(n={int(r['n_set_in_rank'])}; not in headline FDR)"
        )

    md = f"""# FINDING — GSE218989 CLDN4 Q4 vs Q1 prerank GSEA (IFN-γ / MHC-I / TJ)

**Additive only.** Public SMC–KAIST PD-1/PD-L1 NSCLC bulk TPM (Kang et al., *Nat Commun* 2024, PMID 38744958; GEO [GSE218989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218989)). This folder does **not** audit or retract `methods/gse218989_cldn4_ici`. No slide was re-scored. Unit is the **patient**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `gse289287_trop2ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Do not headline response OR.** CLDN4 vs GEO ICI response is already known NS in `methods/gse218989_cldn4_ici` (MWU p=0.40, AUC=0.47; CD8A / Ayers IFN-γ do separate responders). The Q4 vs Q1 2×2 is recorded in `tables/response_or_note.tsv` only.

---

## 一句话 / TL;DR

| Contrast | n (Q4 vs Q1) | Rank | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE218989 CLDN4 Q4 vs Q1 | {n_q4} vs {n_q1} | Welch *t* on log2(TPM+1) | {call_headline(gsea, "welch_t_Q4_minus_Q1")} |

Honest read: {honest_read(gsea)} Histology is **not deposited**. This is bulk ICI RNA — an IFN / MHC NES can be infiltrate, not a tumour-cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author TPM only. No FASTQ / SRA. |
| File | `{TPM_FILE}` |
| n | {n_patients} patients × {n_genes} protein-coding genes |
| Split | CLDN4 `log2(TPM+1)` quartiles: Q1 ≤ {q1_cut:.4f}, Q4 ≥ {q4_cut:.4f} |
| Contrast | CLDN4 Q4 (n={n_q4}) vs Q1 (n={n_q1}) |
| Sign | Welch *t* > 0 = **higher in CLDN4 Q4** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all {n_patients}); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α; KEGG TJ without CLDN4 (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | Headline response OR; LUAD vs LUSC (not deposited); ADC+ICI |

CLDN4 itself is a KEGG TJ member and is the split gene. Its Welch *t* on this contrast is **{cldn4_t:+.2f}**. TJ NES on the full rank is therefore partly circular. The drop-CLDN4 / TJ-minus-CLDN4 rows are the non-circular companions.

---

## Headline NES (Welch *t*, Q4 minus Q1)

{n_ranked} genes ranked. {n_q4} Q4 vs {n_q1} Q1.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_minus_Q1")}

Hallmark IFN-α (secondary, not in BH): {ifna_txt}.

Leading-edge (first 25, Welch rank):

- IFN-γ: `{lead["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}`
- MHC-I / APM: `{lead["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}`
- KEGG TJ: `{lead["KEGG_TIGHT_JUNCTION"]}`

---

## Sensitivity

Same three headline sets.

**Spearman ρ vs continuous CLDN4** (all {n_patients} patients; not a quartile cut):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "spearman_vs_CLDN4")}

**Welch *t* after dropping CLDN4 from the rank** (same Q4/Q1 patients):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "welch_t_Q4_minus_Q1_drop_CLDN4")}

KEGG TJ minus CLDN4 on the primary rank (not in BH): {tj_noc_txt}.

KEGG TJ on the drop-CLDN4 rank (not in BH, same as the TJ row above): {tj_drop_txt}.

---

## Response OR — already known NS; not a headline

From `methods/gse218989_cldn4_ici`: CLDN4 vs GEO response MWU p=0.40, AUC=0.47 (168 R / 187 NR). CD8A p=0.0027 and Ayers IFN-γ p=0.00070 do separate responders, so the endpoint is not dead.

Q4 vs Q1 2×2 on the same GEO labels (note only): Q4 {or_row['q4_r']} R / {or_row['q4_nr']} NR vs Q1 {or_row['q1_r']} R / {or_row['q1_nr']} NR; OR {fmt_or(or_row['or'])} (95% CI {fmt_or(or_row['or_lo'])}–{fmt_or(or_row['or_hi'])}), Fisher p={fmt_p(or_row['p_fisher'])}. That is the already-known null, not a new claim. Table: `tables/response_or_note.tsv`.

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse218989_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_n.tsv`, `response_or_note.tsv`.

---

## What this is not

- Not a re-run or retraction of `methods/gse218989_cldn4_ici`.
- Not a response-OR headline. That test is already NS.
- Not FASTQ / salmon / DESeq2 from SRA (author TPM is used as deposited).
- Not a LUAD vs LUSC split — histology is not on GEO or Supplementary Data 8.
- Not an ADC+ICI trial. ICI-only cohort.
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a tumour-cell-intrinsic IFN/MHC call. Bulk ICI RNA.

---

## 中文摘要

只补公开 **GSE218989** CLDN4 四分位 Q4 vs Q1 的 prerank GSEA（IFN-γ / MHC-I / KEGG TJ），不审不撤已有 `gse218989_cldn4_ici` 页。**不要把缓解 OR 当标题**——那边已经是 NS（GEO MWU p=0.40）。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- {n_q4} Q4 vs {n_q1} Q1；秩 = Welch *t*（log2(TPM+1)）。
- Headline：{call_headline(gsea, "welch_t_Q4_minus_Q1")}。
- CLDN4 在 KEGG TJ 里，又是分组基因，TJ 全秩 NES 有循环成分；去 CLDN4 的伴随结果见上表。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes; CLDN4 in set={('CLDN4' in genes)}")

    tpm_path = download(TPM_URL, CACHE / TPM_FILE)
    mtx_path = download(MATRIX_URL, CACHE / MATRIX_FILE)
    tpm = load_tpm(tpm_path)
    geo = load_geo_meta(mtx_path).set_index("patient")
    log(f"TPM {tpm.shape[0]} genes × {tpm.shape[1]} patients")

    patients = list(tpm.columns)
    assert len(patients) == len(set(patients))
    assert set(patients) == set(geo.index)
    assert "CLDN4" in tpm.index

    logx = np.log2(tpm + 1.0)
    cldn4 = logx.loc["CLDN4"]
    q4, q1 = quartile_split(cldn4)
    n_q4, n_q1 = int(len(q4)), int(len(q1))
    q1_cut = float(cldn4.quantile(0.25))
    q4_cut = float(cldn4.quantile(0.75))
    log(f"CLDN4 Q4 n={n_q4} (≥{q4_cut:.4f}); Q1 n={n_q1} (≤{q1_cut:.4f})")

    rank_t = rank_high_vs_low(logx, q4, q1)
    rank_rho = rank_spearman_vs_target(logx, cldn4)
    rank_t_drop = rank_t.drop(labels=["CLDN4"], errors="ignore")
    cldn4_t = float(rank_t.loc["CLDN4"]) if "CLDN4" in rank_t.index else float("nan")
    log(f"Welch rank {len(rank_t)}; Spearman rank {len(rank_rho)}; CLDN4 t={cldn4_t:+.2f}")

    sets_plus = dict(sets)
    sets_plus["KEGG_TIGHT_JUNCTION_NO_CLDN4"] = [
        g for g in sets["KEGG_TIGHT_JUNCTION"] if g != "CLDN4"
    ]

    g_t = run_rank(rank_t, sets_plus, "welch_t_Q4_minus_Q1", n_q4, n_q1)
    g_rho = run_rank(rank_rho, sets, "spearman_vs_CLDN4", n_q4, n_q1)
    g_drop = run_rank(
        rank_t_drop, sets, "welch_t_Q4_minus_Q1_drop_CLDN4", n_q4, n_q1
    )
    gsea = pd.concat([g_t, g_rho, g_drop], ignore_index=True)

    # Headline table = 3 sets × primary + spearman (drop-CLDN4 is sensitivity, kept in all)
    headline = gsea[gsea["term"].isin(HEADLINE)].copy()

    # Response OR note (GEO labels). Not a headline.
    geo_r = geo.loc[patients, "treatment outcome"]
    y = geo_r.map({"Responder": 1, "Non-responder": 0})
    q4_r = int(((y.loc[q4] == 1)).sum())
    q4_nr = int(((y.loc[q4] == 0)).sum())
    q1_r = int(((y.loc[q1] == 1)).sum())
    q1_nr = int(((y.loc[q1] == 0)).sum())
    or_stats = fisher_or(q4_r, q4_nr, q1_r, q1_nr)
    # continuous MWU echo of the known result
    mwu = stats.mannwhitneyu(
        cldn4.loc[y[y == 1].index],
        cldn4.loc[y[y == 0].index],
        alternative="two-sided",
    )
    or_row = {
        "source": "GEO_treatment_outcome",
        "note": "already_known_NS_do_not_headline",
        "known_from": "methods/gse218989_cldn4_ici MWU p=0.40 AUC=0.47",
        "n_q4": n_q4,
        "n_q1": n_q1,
        "q4_r": q4_r,
        "q4_nr": q4_nr,
        "q1_r": q1_r,
        "q1_nr": q1_nr,
        "or": or_stats["or"],
        "or_lo": or_stats["or_lo"],
        "or_hi": or_stats["or_hi"],
        "p_fisher": or_stats["p_fisher"],
        "continuous_MWU_p": float(mwu.pvalue),
        "n_GEO_R": int((y == 1).sum()),
        "n_GEO_NR": int((y == 0).sum()),
    }

    quart = pd.DataFrame(
        [
            {
                "gene": "CLDN4",
                "transform": "log2(TPM+1)",
                "n_patients": int(len(patients)),
                "q1_cut": q1_cut,
                "q4_cut": q4_cut,
                "n_q1": n_q1,
                "n_q4": n_q4,
                "n_middle": int(len(patients) - n_q1 - n_q4),
                "cldn4_welch_t_Q4_minus_Q1": cldn4_t,
                "cldn4_in_KEGG_TJ": True,
            }
        ]
    )
    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "NSCLC ICI bulk TPM (SMC-KAIST); histology not deposited",
                "n_patients": int(len(patients)),
                "n_genes": int(tpm.shape[0]),
                "n_q4": n_q4,
                "n_q1": n_q1,
                "rank_metric_primary": "Welch t on log2(TPM+1), CLDN4 Q4 minus Q1",
                "n_genes_ranked": int(rank_t.shape[0]),
                "tpm_file": TPM_FILE,
                "series_title": "SMC-KAIST PD-1/PD-L1 NSCLC (Kang 2024)",
                "not": "response-OR headline; FASTQ; LUAD/LUSC split; ADC+ICI",
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    quart.to_csv(TAB / "quartile_n.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    pd.DataFrame([or_row]).to_csv(TAB / "response_or_note.tsv", sep="\t", index=False)
    rank_t.rename("welch_t").to_csv(TAB / "rank_welch_t.tsv", sep="\t")
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline)
    write_finding(
        gsea,
        figures,
        int(rank_t.shape[0]),
        n_q4,
        n_q1,
        int(tpm.shape[0]),
        int(tpm.shape[1]),
        q4_cut,
        q1_cut,
        or_row,
        cldn4_t,
    )

    print("\n=== HEADLINE NES (Welch t Q4-Q1) ===")
    show = headline[headline["rank_metric"] == "welch_t_Q4_minus_Q1"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== RESPONSE OR (note only) ===")
    print(
        f"Q4 {q4_r}/{q4_nr} vs Q1 {q1_r}/{q1_nr}  OR {or_row['or']:.2f} "
        f"p={or_row['p_fisher']:.3g}  continuous MWU p={or_row['continuous_MWU_p']:.3g}"
    )


if __name__ == "__main__":
    main()
