#!/usr/bin/env python3
"""ADDITIVE extra — prerank GSEA on GSE334497 4T1 Trop2 KO tumors (5 vs 5).

Same engine as methods/scrna_pseudobulk_gsea_meta / scripts/cldn4_ko_gsea
(gsea_core.py): weighted KS p=1, 1000 gene-set permutations, seed=42.

Positive NES = the set is enriched at the Trop2-KO end of the rank
(Welch t, KO minus WT, on log2(norm+1)).

Headline sets (user-requested): IFN-γ, MHC-I/APM, KEGG tight junction, keratin.

Cldn4 gene-level p=0.32 is taken as given. This script reports Cldn4 log2FC
as an effect size and does not re-test that gene.

Not lung. Not SKB264. Not a CLDN4 KO. Additive only — does not audit or
retract any existing slide.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

DATA = ROOT / "methods" / "gse334497_trop2ko_gsea" / "data"
OUT = ROOT / "methods" / "gse334497_trop2ko_gsea"
FIG = OUT / "figures"
TAB = OUT / "tables"
SET_JSON = ROOT / "data" / "genesets" / "a8_sets.json"
MM_GMT = ROOT / "data" / "genesets" / "mh.all.v2023.2.Mm.symbols.gmt"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# Library names from GEO series matrix Sample_description (Wu et al., JITC 2026)
KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["control170", "RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R"]

# User-requested headline sets. FDR is BH among these four only.
HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "KEGG_TIGHT_JUNCTION",
    "GOBP_KERATINIZATION",
]
HEADLINE_LABEL = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
    "KEGG_TIGHT_JUNCTION": "KEGG tight junction",
    "GOBP_KERATINIZATION": "GO keratinization",
}
MOUSE_MHC_OVERRIDE = {
    "HLA-A": ["H2-K1", "H2-D1"],
    "HLA-B": ["H2-K1", "H2-D1"],
    "HLA-C": ["H2-D1"],
    "HLA-E": ["H2-T23"],
    "HLA-F": ["H2-Q7", "H2-Q6"],
    "HLA-G": ["H2-Q10"],
    "ERAP2": ["Lnpep"],
}

# Given. Do not re-audit this gene-level test.
CLDN4_P_GIVEN = 0.32
CONTRAST = "GSE334497_4T1_Trop2KO"
PSEUDO = 1.0


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fmt_p(p) -> str:
    if p is None or p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x) -> str:
    if x is None or x != x:
        return "NA"
    return f"{x:+.3f}"


def fmt_fc(x) -> str:
    if x is None or x != x:
        return "NA"
    return f"{x:+.3f}"


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > 2:
                sets[p[0]] = [g for g in p[2:] if g]
    return sets


def to_mouse(symbols: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for s in symbols:
        cands = MOUSE_MHC_OVERRIDE[s] if s in MOUSE_MHC_OVERRIDE else [s[0] + s[1:].lower()]
        for m in cands:
            if m and m not in seen:
                seen.add(m)
                out.append(m)
    return out


def load_mouse_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    human = {k: list(payload["sets"][k]) for k in HEADLINE}
    mm_hallmark = read_gmt(MM_GMT)
    mouse: dict[str, list[str]] = {}
    for k in HEADLINE:
        if k in mm_hallmark:
            mouse[k] = mm_hallmark[k]
        else:
            mouse[k] = to_mouse(human[k])
    return mouse


def load_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(DATA / "GSE334497_normalized_counts.csv.gz", index_col=0)
    raw.index = raw.index.astype(str)
    for c in KO + WT:
        if c not in raw.columns:
            raise SystemExit(f"missing library {c} in counts; have {list(raw.columns)}")
    annot = pd.read_csv(DATA / "ensembl_to_symbol.tsv", sep="\t")
    ens2sym = dict(zip(annot["ensembl"].astype(str), annot["symbol"].astype(str)))
    df = raw[KO + WT].copy()
    df["symbol"] = [ens2sym.get(i, "") for i in df.index]
    df = df[df["symbol"].astype(str).str.len() > 0]
    df["_mean"] = df[KO + WT].mean(axis=1)
    df = (
        df.sort_values("_mean", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")
    )
    counts = df[KO + WT].astype(float)
    lg = np.log2(counts + PSEUDO)
    return counts, lg


def welch_t_rank(lg: pd.DataFrame) -> pd.Series:
    """Signed Welch t, KO minus WT. Positive = up in Trop2 KO."""
    eh = lg[KO].to_numpy(dtype=np.float64)
    el = lg[WT].to_numpy(dtype=np.float64)
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    ok &= (vh + vl) > 1e-8
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    tstat = np.full(lg.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    s = pd.Series(tstat, index=lg.index).replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def gene_log2fc(lg: pd.DataFrame, symbol: str) -> dict:
    if symbol not in lg.index:
        return {"symbol": symbol, "present": False, "log2FC": np.nan}
    ko = lg.loc[symbol, KO].astype(float)
    wt = lg.loc[symbol, WT].astype(float)
    return {
        "symbol": symbol,
        "present": True,
        "log2FC": float(ko.mean() - wt.mean()),
        "mean_log2_KO": float(ko.mean()),
        "mean_log2_WT": float(wt.mean()),
        "values_KO": ",".join(f"{v:.4f}" for v in ko.to_numpy()),
        "values_WT": ",".join(f"{v:.4f}" for v in wt.to_numpy()),
    }


def run_one(rank: pd.Series, sets: dict[str, list[str]], meta: dict) -> pd.DataFrame:
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    if g.empty:
        raise SystemExit("no gene set passed min_size")
    g["fdr"] = bh_fdr(g["nom_p"])
    g["contrast"] = CONTRAST
    g["accession"] = meta["accession"]
    g["n_loss"] = meta["n_loss"]
    g["n_wt"] = meta["n_wt"]
    g["rank_metric"] = meta["rank_metric"]
    g["n_genes_ranked"] = int(rank.shape[0])
    g["label"] = g["term"].map(HEADLINE_LABEL)
    return g


def write_figures(gsea: pd.DataFrame, cldn4: dict, tacstd2: dict) -> list[str]:
    written: list[str] = []
    sub = gsea.set_index("term").reindex(HEADLINE)
    nes = sub["nes"].astype(float).to_numpy()
    fdr = sub["fdr"].astype(float).to_numpy()
    labels = [HEADLINE_LABEL[t] for t in HEADLINE]
    colors = ["#4C78A8" if v >= 0 else "#E45756" for v in nes]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax = axes[0]
    y = np.arange(len(HEADLINE))
    ax.barh(y, nes, color=colors, zorder=3, height=0.62)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("NES (positive = up in Trop2 KO)")
    ax.set_title(
        "GSE334497 4T1 Trop2 KO vs WT — prerank GSEA\n"
        f"* BH-FDR < 0.05 among {len(HEADLINE)} headline sets; "
        f"{NPERM} gene-set perms, seed={SEED}"
    )
    xmax = max(2.2, float(np.nanmax(np.abs(nes))) + 0.35)
    ax.set_xlim(-xmax, xmax)
    for i, (nv, q, nset) in enumerate(zip(nes, fdr, sub["n_set_in_rank"].to_numpy())):
        star = "*" if pd.notna(q) and q < 0.05 else ""
        side = 0.08 if nv >= 0 else -0.08
        ax.text(
            nv + side,
            i,
            f"{nv:+.2f}{star}  FDR {fmt_p(q)}  n={int(nset)}",
            va="center",
            ha="left" if nv >= 0 else "right",
            fontsize=8,
        )
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, ls=":", alpha=0.5)

    ax = axes[1]
    genes = ["Tacstd2", "Cldn4"]
    vals = [tacstd2.get("log2FC", np.nan), cldn4.get("log2FC", np.nan)]
    bar_c = ["#4C78A8" if (v == v and v >= 0) else "#E45756" for v in vals]
    ax.barh(np.arange(2), [0 if v != v else v for v in vals], color=bar_c, height=0.55, zorder=3)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Tacstd2 (KO QC)", "Cldn4 (effect size)"])
    ax.invert_yaxis()
    ax.set_xlabel("log2FC (KO − WT), log2(norm+1)")
    ax.set_title("Gene-level log2FC (not a re-test of Cldn4 p)")
    for i, (g, v) in enumerate(zip(genes, vals)):
        if v == v:
            extra = "  given p=0.32, not re-tested" if g == "Cldn4" else "  perturbation QC"
            ax.text(
                v + (0.12 if v >= 0 else -0.12),
                i,
                f"{v:+.2f}{extra}",
                va="center",
                ha="left" if v >= 0 else "right",
                fontsize=8,
            )
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, ls=":", alpha=0.5)

    fig.tight_layout()
    p1 = FIG / "fig_headline_nes.png"
    fig.savefig(p1, dpi=160)
    fig.savefig(FIG / "fig_headline_nes.pdf")
    plt.close(fig)
    written.append(str(p1.relative_to(ROOT)))

    # Compact NES-only bar (the "NES figure" for the PR)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(HEADLINE))
    ax.bar(x, nes, color=colors, zorder=3, width=0.62)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("NES (positive = up in Trop2 KO)")
    ax.set_title(
        "GSE334497 4T1 Trop2 KO 5 vs 5 — headline prerank NES\n"
        f"Cldn4 log2FC {fmt_fc(cldn4.get('log2FC'))} (p=0.32 given, not re-tested)"
    )
    for i, (nv, q) in enumerate(zip(nes, fdr)):
        star = "*" if pd.notna(q) and q < 0.05 else ""
        ax.text(
            i,
            nv + (0.08 if nv >= 0 else -0.08),
            f"{nv:+.2f}{star}\nFDR {fmt_p(q)}",
            ha="center",
            va="bottom" if nv >= 0 else "top",
            fontsize=8,
        )
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()
    p2 = FIG / "fig_nes_bars.png"
    fig.savefig(p2, dpi=160)
    fig.savefig(FIG / "fig_nes_bars.pdf")
    plt.close(fig)
    written.append(str(p2.relative_to(ROOT)))
    return written


def write_finding(gsea: pd.DataFrame, cldn4: dict, tacstd2: dict, figures: list[str], n_ranked: int) -> None:
    def row(term: str) -> pd.Series:
        hit = gsea[gsea["term"] == term]
        if hit.empty:
            return pd.Series({"nes": np.nan, "fdr": np.nan, "nom_p": np.nan, "n_set_in_rank": 0, "es": np.nan, "lead_genes": ""})
        return hit.iloc[0]

    def nes_line() -> str:
        parts = []
        for t in HEADLINE:
            r = row(t)
            parts.append(
                f"| {HEADLINE_LABEL[t]} | `{t}` | {fmt_nes(r['nes'])} | {fmt_p(r['fdr'])} | {fmt_p(r['nom_p'])} | {int(r['n_set_in_rank']) if pd.notna(r['n_set_in_rank']) else 0} |"
            )
        return "\n".join(parts)

    ifn = row("HALLMARK_INTERFERON_GAMMA_RESPONSE")
    mhc = row("CUSTOM_MHC_I_ANTIGEN_PRESENTATION")
    tj = row("KEGG_TIGHT_JUNCTION")
    krt = row("GOBP_KERATINIZATION")

    def dir_call(r: pd.Series) -> str:
        if pd.isna(r.get("nes")):
            return "not tested"
        side = "UP" if r["nes"] > 0 else "DOWN"
        sig = "FDR<0.05" if pd.notna(r.get("fdr")) and r["fdr"] < 0.05 else "FDR≥0.05"
        return f"{side} NES {fmt_nes(r['nes'])} ({sig})"

    md = f"""# FINDING — GSE334497 4T1 Trop2 KO prerank GSEA

**Additive only.** This folder does **not** audit or retract any slide. It is **not** lung and it is **not** SKB264. It is **not** a CLDN4 knockout: the perturbation is CRISPR **Trop2 (Tacstd2)** in 4T1 mammary tumors.

Prerank engine is the same as `methods/scrna_pseudobulk_gsea_meta` (`scripts/gse334497_trop2ko_gsea/gsea_core.py`): weighted KS *p*=1, **{NPERM}** gene-set permutations, seed=**{SEED}**. Positive NES = the set is enriched at the **Trop2-KO** end of the rank (Welch *t*, KO minus WT). BH-FDR is within the **four headline sets**. Numbers below are written from `tables/gsea_headline.tsv` and `tables/cldn4_log2fc.tsv`.

**Cldn4 gene-level *p* = 0.32 is given and is not re-tested.** Only the log2FC effect size is reported.

---

## 一句话 / TL;DR

| Contrast | n | Tacstd2 log2FC (QC) | Cldn4 log2FC | Cldn4 *p* | IFN-γ | MHC-I/APM | KEGG TJ | keratin |
|---|---|---|---|---|---|---|---|---|
| GSE334497 4T1 Trop2 KO vs WT | **5 vs 5** | {fmt_fc(tacstd2.get("log2FC"))} | {fmt_fc(cldn4.get("log2FC"))} | **0.32 (given)** | {dir_call(ifn)} | {dir_call(mhc)} | {dir_call(tj)} | {dir_call(krt)} |

RNA is from **frozen whole-tumor sections** in immunocompetent BALB/c hosts. An IFN / MHC NES on this matrix can be stroma or infiltrate; it is not a tumor-cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public processed file only: `GSE334497_normalized_counts.csv.gz` (GRCm38 Ensembl). No FASTQ / SRA. |
| Model | Mouse **4T1** mammary tumors, CRISPR Trop2 KO vs WT, 3 weeks in BALB/c (Wu *et al.*, *JITC* 2026). Breast, not lung. |
| Libraries | KO: KO162, KO164, KO165, KO172, RESUB-KO163R. WT: control170, RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R. |
| Transform | log2(normalized count + 1); Ensembl → official symbol, collapse by max mean. |
| Rank | Signed Welch *t*, KO minus WT. Positive = up in Trop2 KO. |
| Engine | Same prerank as `gsea_core.py` (weighted KS *p*=1, 1000 gene-set perms, seed=42, min size 8) |
| Headline sets | Hallmark IFN-γ (MSigDB `mh.all.v2023.2.Mm`); custom MHC-I/APM (human → mouse, HLA→H2); KEGG tight junction; GO keratinization |
| FDR | BH inside the **4** headline sets |
| Cldn4 gene test | **Not re-run.** Given *p* = 0.32. This folder reports log2FC only. |

Primary genes: mouse Cldn4 `ENSMUSG00000047501`; Tacstd2 `ENSMUSG00000051397` (perturbation QC only; not a GSEA-set member).

---

## Headline NES

GEO: [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497). {n_ranked} genes ranked.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line()}

| Gene | log2FC (KO − WT) | Note |
|---|---:|---|
| Tacstd2 | {fmt_fc(tacstd2.get("log2FC"))} | Perturbation QC. Knockout worked. |
| Cldn4 | {fmt_fc(cldn4.get("log2FC"))} | Effect size only. Gene-level *p* = **0.32 (given; not re-tested)**. |

Leading-edge (first 25, from the walk):

| Set | Lead genes |
|---|---|
| Hallmark IFN-γ | `{row("HALLMARK_INTERFERON_GAMMA_RESPONSE").get("lead_genes", "")}` |
| MHC-I / APM | `{row("CUSTOM_MHC_I_ANTIGEN_PRESENTATION").get("lead_genes", "")}` |
| KEGG tight junction | `{row("KEGG_TIGHT_JUNCTION").get("lead_genes", "")}` |
| GO keratinization | `{row("GOBP_KERATINIZATION").get("lead_genes", "")}` |

---

## What this is and is not

**Is**

- An additive prerank GSEA on the public 4T1 Trop2 KO 5-vs-5 tumor matrix.
- A report of NES / FDR for IFN-γ, MHC-I/APM, KEGG TJ, and keratin, plus Cldn4 log2FC.

**Is not**

- A CLDN4 knockdown or knockout.
- A lung model or SKB264.
- A re-audit of the given Cldn4 gene-level *p* = 0.32.
- Tumor-cell-intrinsic IFN/APM (bulk immunocompetent tumor RNA).
- Sample-permutation GSEA. The null is gene-set permutation on a prerank.

---

## Figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/gse334497_trop2ko_gsea/tables/gsea_headline.tsv`, `gsea_prerank_all.tsv`, `cldn4_log2fc.tsv`, `inventory.tsv`.

---

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/gse334497_trop2ko_gsea/analyze.py
```

---

## 中文摘要

只补 GSE334497 **4T1 Trop2 KO 肿瘤 5 vs 5** 的 prerank GSEA，不审不撤已有页。不是肺，不是 SKB264，也不是 CLDN4 KO。引擎与 `scrna_pseudobulk_gsea_meta` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = Trop2 KO 端富集。Cldn4 单基因 *p*=0.32 **按给定，不重测**，只报 log2FC。

- Tacstd2 log2FC {fmt_fc(tacstd2.get("log2FC"))}（敲除成立）。
- Cldn4 log2FC {fmt_fc(cldn4.get("log2FC"))}（*p*=0.32 给定）。
- IFN-γ：{dir_call(ifn)}。
- MHC-I/APM：{dir_call(mhc)}。
- KEGG 紧密连接：{dir_call(tj)}。
- GO 角化：{dir_call(krt)}。
- 全瘤冰冻切片、免疫健全宿主，IFN/MHC NES 可以来自间质/浸润，不能写成肿瘤内在。
"""
    (OUT / "FINDING.md").write_text(md)
    log(f"wrote {OUT / 'FINDING.md'}")


def main() -> None:
    sets = load_mouse_sets()
    for t, genes in sets.items():
        log(f"set {t}: {len(genes)} genes")

    counts, lg = load_matrix()
    log(f"matrix symbols={lg.shape[0]} libraries={list(lg.columns)}")

    rank = welch_t_rank(lg)
    log(f"ranked {len(rank)} genes (Welch t, KO minus WT)")

    tacstd2 = gene_log2fc(lg, "Tacstd2")
    cldn4 = gene_log2fc(lg, "Cldn4")
    cldn4["gene_p_given"] = CLDN4_P_GIVEN
    cldn4["gene_p_note"] = "given; not re-tested"
    log(f"Tacstd2 log2FC={fmt_fc(tacstd2.get('log2FC'))}")
    log(f"Cldn4 log2FC={fmt_fc(cldn4.get('log2FC'))} (p={CLDN4_P_GIVEN} given, not re-tested)")

    meta = {
        "accession": "GSE334497",
        "contrast": "4T1 Trop2 CRISPR KO vs WT tumors (bulk frozen section)",
        "species": "mouse",
        "tissue": "4T1 mammary tumor, immunocompetent BALB/c (not lung)",
        "n_loss": 5,
        "n_wt": 5,
        "n_note": "5 vs 5 biological tumors; author-normalized counts",
        "rank_metric": "Welch t on log2(norm+1), KO minus WT",
        "perturb_gene": "Tacstd2",
    }
    gsea = run_one(rank, sets, meta)
    headline = gsea[gsea["term"].isin(HEADLINE)].copy()

    targets = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": "GSE334497",
                **cldn4,
                "n_KO": 5,
                "n_WT": 5,
                "tissue": meta["tissue"],
            }
        ]
    )
    qc = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                "accession": "GSE334497",
                **tacstd2,
                "role": "perturbation_QC",
                "n_KO": 5,
                "n_WT": 5,
            }
        ]
    )
    inv = pd.DataFrame(
        [
            {
                "contrast": CONTRAST,
                **meta,
                "n_genes_ranked": int(rank.shape[0]),
                "n_symbols_mapped": int(lg.shape[0]),
                "cldn4_p_given": CLDN4_P_GIVEN,
            }
        ]
    )

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    targets.to_csv(TAB / "cldn4_log2fc.tsv", sep="\t", index=False)
    qc.to_csv(TAB / "tacstd2_log2fc.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    rank.rename("welch_t_KO_minus_WT").to_csv(TAB / "rank_welch_t.tsv", sep="\t")
    log(f"wrote tables under {TAB}")

    figures = write_figures(gsea, cldn4, tacstd2)
    log("figures: " + ", ".join(figures))
    write_finding(gsea, cldn4, tacstd2, figures, int(rank.shape[0]))

    print("\n=== HEADLINE NES ===")
    show = headline[["term", "nes", "fdr", "nom_p", "n_set_in_rank"]].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== CLDN4 (effect size only; p=0.32 given) ===")
    print(f"log2FC={fmt_fc(cldn4.get('log2FC'))}  given_p={CLDN4_P_GIVEN}")
    print(f"Tacstd2 log2FC={fmt_fc(tacstd2.get('log2FC'))}")


if __name__ == "__main__":
    main()
