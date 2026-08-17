#!/usr/bin/env python3
"""ADDITIVE extra — prerank GSEA on public GSE289287 T-47D Trop-2 KO xenografts.

Uses the author DESeq2 table Trop2KO tumors vs WT (log2FC>0 / stat>0 = up in KO).
Engine is the same as methods/cldn4_ko_gsea and scrna_pseudobulk_gsea_meta
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Headline sets requested: Hallmark IFN-γ, MHC-I/APM, KEGG tight junction.
Hallmark IFN-α is run as a secondary set only (not in the BH-FDR headline).

Not lung. Not SKB264. Not FASTQ / SRA. CLDN4 padj ~0.47 is the known
single-gene result and is recorded, not re-litigated.
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
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

DATA = HERE / "data"
FIG = HERE / "figures"
TAB = HERE / "tables"
SET_JSON = HERE / "gene_sets.json"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

DESEQ_FILE = "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz"
CONTRAST = "GSE289287_T47D_Trop2KO_xenograft"
ACCESSION = "GSE289287"

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

# Author table columns: WT animals 2808/2810/2812; Trop-2 KO 2807/2815/2817/2818
WT_SAMPLES = ["OV.2808.RNA", "OV.2810.RNA", "OV.2812.RNA"]
KO_SAMPLES = ["OV.2807.RNA", "OV.2815.RNA", "OV.2817.RNA", "OV.2818.RNA"]
WT_GSM = ["GSM8788420", "GSM8788421", "GSM8788419"]
KO_GSM = ["GSM8788425", "GSM8788422", "GSM8788423", "GSM8788424"]


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


def load_deseq() -> pd.DataFrame:
    df = pd.read_csv(DATA / DESEQ_FILE, sep="\t")
    df["symbol"] = df["Feature_name"].astype(str).str.strip()
    df = df[df["symbol"].notna() & (df["symbol"] != "") & (df["symbol"] != "nan")]
    df["stat"] = pd.to_numeric(df["stat"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["baseMean"] = pd.to_numeric(df["baseMean"], errors="coerce")
    # one row per symbol already; keep highest |stat| if a collision appears
    df = df.assign(_abs=df["stat"].abs()).sort_values("_abs", ascending=False)
    df = df.drop_duplicates("symbol", keep="first")
    return df


def gene_row(df: pd.DataFrame, symbol: str) -> dict:
    hit = df[df["symbol"].str.upper() == symbol.upper()]
    if hit.empty:
        return {"symbol": symbol, "present": False}
    r = hit.iloc[0]
    return {
        "symbol": symbol,
        "present": True,
        "ensembl": r.get("Ensembl_Id"),
        "log2FC": float(r["log2FoldChange"]),
        "stat": float(r["stat"]),
        "pvalue": float(r["pvalue"]) if pd.notna(r["pvalue"]) else np.nan,
        "padj": float(r["padj"]) if pd.notna(r["padj"]) else np.nan,
        "baseMean": float(r["baseMean"]) if pd.notna(r["baseMean"]) else np.nan,
    }


def run_rank(rank: pd.Series, sets: dict[str, list[str]], rank_metric: str) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    g["contrast"] = CONTRAST
    g["accession"] = ACCESSION
    g["rank_metric"] = rank_metric
    g["n_genes_ranked"] = int(rank.shape[0])
    g["n_ko"] = 4
    g["n_wt"] = 3
    g["headline"] = g["term"].isin(HEADLINE)
    head = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if head.any():
        g.loc[head, "fdr"] = bh_fdr(g.loc[head, "nom_p"]).to_numpy()
    # secondary FDR is the nominal p (single extra set); leave NaN and report nom_p
    g["direction"] = np.where(g["nes"] > 0, "UP_in_Trop2KO", "DOWN_in_Trop2KO")
    return g


def write_figures(headline: pd.DataFrame, genes: pd.DataFrame) -> list[str]:
    h = headline[headline["rank_metric"] == "DESeq2_Wald_stat"].copy()
    h["label"] = h["term"].map(HEADLINE_LABEL)
    order = [HEADLINE_LABEL[t] for t in HEADLINE]
    h = h.set_index("label").reindex(order).reset_index()
    colors = ["#c0392b" if n > 0 else "#2980b9" for n in h["nes"]]

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(h))
    ax.barh(y, h["nes"], color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(h["label"])
    ax.invert_yaxis()
    ax.set_xlabel("NES (positive = up in Trop-2 KO)")
    ax.set_title("GSE289287 T-47D Trop-2 KO xenograft — headline prerank GSEA")
    for i, row in h.iterrows():
        txt = f"NES {row['nes']:+.2f}  FDR {fmt_p(row['fdr'])}"
        x = row["nes"]
        ax.text(x + (0.04 if x >= 0 else -0.04), i, txt, va="center",
                ha="left" if x >= 0 else "right", fontsize=8)
    fig.tight_layout()
    png = FIG / "fig_headline_nes.png"
    pdf = FIG / "fig_headline_nes.pdf"
    fig.savefig(png, dpi=160)
    fig.savefig(pdf)
    plt.close(fig)

    # small two-panel: TACSTD2 / CLDN4 + NES
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    gdf = genes.set_index("symbol")
    names = ["TACSTD2", "CLDN4"]
    fcs = [gdf.loc[n, "log2FC"] if n in gdf.index else np.nan for n in names]
    padjs = [gdf.loc[n, "padj"] if n in gdf.index else np.nan for n in names]
    axes[0].bar(names, fcs, color=["#8e44ad", "#7f8c8d"], edgecolor="black", linewidth=0.4)
    axes[0].axhline(0, color="black", lw=0.7)
    axes[0].set_ylabel("author DESeq2 log2FC (KO vs WT)")
    axes[0].set_title("Perturbation QC / known CLDN4")
    for i, (fc, p) in enumerate(zip(fcs, padjs)):
        axes[0].text(i, fc, f"{fc:+.2f}\npadj {fmt_p(p)}", ha="center",
                     va="bottom" if fc >= 0 else "top", fontsize=8)
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


def call_headline(gsea: pd.DataFrame, rank_metric: str, terms=None) -> str:
    terms = list(terms) if terms is not None else list(HEADLINE)
    sub = gsea[(gsea["rank_metric"] == rank_metric) & (gsea["term"].isin(terms))]
    bits = []
    for term in terms:
        r = sub[sub["term"] == term].iloc[0]
        bits.append(
            f"{HEADLINE_LABEL[term]} NES {fmt_nes(r['nes'])} FDR {fmt_p(r['fdr'])}"
        )
    return "; ".join(bits)


def write_finding(gsea: pd.DataFrame, genes: pd.DataFrame, figures: list[str], n_ranked: int) -> None:
    tac = genes[genes["symbol"] == "TACSTD2"].iloc[0].to_dict()
    cld = genes[genes["symbol"] == "CLDN4"].iloc[0].to_dict()
    ifna = gsea[
        (gsea["rank_metric"] == "DESeq2_Wald_stat")
        & (gsea["term"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE")
    ]
    ifna_txt = "not ranked"
    if not ifna.empty:
        r = ifna.iloc[0]
        ifna_txt = f"NES {fmt_nes(r['nes'])} nom p {fmt_p(r['nom_p'])} (n={int(r['n_set_in_rank'])}; not in headline FDR)"

    lead = {}
    for term in HEADLINE:
        r = gsea[
            (gsea["rank_metric"] == "DESeq2_Wald_stat") & (gsea["term"] == term)
        ].iloc[0]
        lead[term] = r["lead_genes"]

    md = f"""# FINDING — GSE289287 T-47D Trop-2 KO xenograft prerank GSEA

**Additive only.** Public author DESeq2 on T-47D **Trop-2 KO vs WT xenografts**. This folder does **not** audit or retract any slide. It is **not** lung and it is **not** SKB264.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `scrna_pseudobulk_gsea_meta` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **Trop-2 KO** end of the rank (KO minus WT). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

CLDN4 single-gene padj **{fmt_p(cld['padj'])}** is the known result (log2FC {fmt_fc(cld['log2FC'])}); it is recorded, not treated as a new claim.

---

## 一句话 / TL;DR

| Contrast | n (KO vs WT) | TACSTD2 log2FC (padj) | CLDN4 log2FC (padj) | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|---|
| GSE289287 T-47D Trop-2 KO xenograft | 4 vs 3 | {fmt_fc(tac['log2FC'])} ({fmt_p(tac['padj'])}) | {fmt_fc(cld['log2FC'])} ({fmt_p(cld['padj'])}) | {call_headline(gsea, 'DESeq2_Wald_stat')} |

Primary rank is the author DESeq2 **Wald statistic**. Hosts are **NRG** females (mammary fat pad), so any IFN / APM movement is **tumour-cell-intrinsic**, not adaptive infiltrate.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author DESeq2 only. No FASTQ / SRA. |
| File | `GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz` |
| Contrast | Trop-2 KO xenografts (GSM8788422/23/24/25; animals 2815/2817/2818/2807) vs WT (GSM8788419/20/21; animals 2812/2808/2810) |
| Sign | log2FC>0 and Wald stat>0 = **up in Trop-2 KO** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | DESeq2 Wald `stat` |
| Rank (sensitivity) | DESeq2 `log2FoldChange` |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | DSG2 KO tables; in-vitro T-47D (no Trop-2 KO cell DESeq2 on GEO); lung; SKB264 |

Series: [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287). Title on GEO: *Trop-2 governs anti-metastatic desmosomal integrity* (Vacek / Souček Lab). Platform GPL24676, human T-47D in NRG mammary fat pad.

---

## Headline NES (Wald-stat rank)

{n_ranked} genes ranked. TACSTD2 itself: log2FC **{fmt_fc(tac['log2FC'])}**, Wald **{fmt_fc(tac['stat'])}**, padj **{fmt_p(tac['padj'])}** (KO worked). CLDN4: log2FC **{fmt_fc(cld['log2FC'])}**, padj **{fmt_p(cld['padj'])}** (known, ns).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "DESeq2_Wald_stat")}

Hallmark IFN-α (secondary, not in BH): {ifna_txt}.

Leading-edge (first 25, Wald rank):

- IFN-γ: `{lead["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}`
- MHC-I / APM: `{lead["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}`
- KEGG TJ: `{lead["KEGG_TIGHT_JUNCTION"]}`

---

## Sensitivity: log2FC rank

Same three sets, rank = author `log2FoldChange` instead of Wald stat.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_md_rows(gsea, "DESeq2_log2FC")}

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse289287_trop2ko_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `key_genes.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-run of any existing slide or of SKB264.
- Not FASTQ / salmon / DESeq2 from SRA (author table is used as deposited).
- Not lung. T-47D is luminal breast; hosts are NRG mice.
- Not the DSG2-KO arms of the same series.
- Not sample-permutation GSEA. n = 4 vs 3; the engine is gene-set permutation on a prerank.

---

## 中文摘要

只补公开 **GSE289287** T-47D Trop-2 KO 移植瘤的 author DESeq2 prerank GSEA，不审不撤已有页。不是肺癌，不是 SKB264。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = Trop-2 KO 端富集。FDR 只在三个 headline set 内做 BH。

- 4 KO vs 3 WT（NRG 乳腺脂肪垫）。TACSTD2 log2FC {fmt_fc(tac['log2FC'])}，padj {fmt_p(tac['padj'])}（敲除成立）。
- CLDN4 log2FC {fmt_fc(cld['log2FC'])}，padj {fmt_p(cld['padj'])}（已知，单基因不显著）。
- Headline（Wald 秩）：{call_headline(gsea, "DESeq2_Wald_stat")}。
"""
    (HERE / "FINDING.md").write_text(md)
    log(f"wrote {HERE / 'FINDING.md'}")


def main() -> None:
    sets = load_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    df = load_deseq()
    log(f"DESeq2 rows after symbol collapse: {len(df)}")

    rank_stat = (
        df.dropna(subset=["stat"])
        .set_index("symbol")["stat"]
        .astype(float)
        .sort_values(ascending=False)
    )
    rank_lfc = (
        df.dropna(subset=["log2FoldChange"])
        .set_index("symbol")["log2FoldChange"]
        .astype(float)
        .sort_values(ascending=False)
    )
    log(f"Wald rank {len(rank_stat)}; log2FC rank {len(rank_lfc)}")

    g_stat = run_rank(rank_stat, sets, "DESeq2_Wald_stat")
    g_lfc = run_rank(rank_lfc, sets, "DESeq2_log2FC")
    gsea = pd.concat([g_stat, g_lfc], ignore_index=True)
    headline = gsea[gsea["headline"]].copy()

    key_syms = ["TACSTD2", "CLDN4", "CLDN7", "B2M", "HLA-A", "HLA-B", "HLA-C",
                "STAT1", "ISG15", "IFI44L", "MX1"]
    genes = pd.DataFrame([gene_row(df, s) for s in key_syms])

    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": ACCESSION,
                "species": "human",
                "tissue": "T-47D luminal breast xenograft, NRG mammary fat pad (not lung)",
                "n_ko": 4,
                "n_wt": 3,
                "ko_gsm": ",".join(KO_GSM),
                "wt_gsm": ",".join(WT_GSM),
                "ko_animals": "2815,2817,2818,2807",
                "wt_animals": "2812,2808,2810",
                "rank_metric_primary": "author DESeq2 Wald stat (KO minus WT)",
                "n_genes_ranked": int(rank_stat.shape[0]),
                "deseq_file": DESEQ_FILE,
                "series_title": "Trop-2 governs anti-metastatic desmosomal integrity",
                "not": "lung; SKB264; FASTQ; DSG2 KO arms",
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    genes.to_csv(TAB / "key_genes.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    figures = write_figures(headline, genes)
    write_finding(gsea, genes, figures, int(rank_stat.shape[0]))

    print("\n=== HEADLINE NES (Wald) ===")
    show = headline[headline["rank_metric"] == "DESeq2_Wald_stat"][
        ["term", "nes", "fdr", "nom_p", "n_set_in_rank"]
    ].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== KEY GENES ===")
    print(genes[["symbol", "log2FC", "stat", "pvalue", "padj"]].to_string(index=False))


if __name__ == "__main__":
    main()
