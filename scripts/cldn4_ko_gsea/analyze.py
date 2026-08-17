#!/usr/bin/env python3
"""ADDITIVE extra — prerank GSEA on public CLDN4 / Cldn4 loss matrices.

Same engine as methods/scrna_pseudobulk_gsea_meta (gsea_core.py): weighted KS
p=1, 1000 gene-set permutations, seed=42. Positive NES = enriched after
CLDN4 / Cldn4 loss (KO or siRNA vs WT / CLDN4-high).

Does not re-run or retract any existing slide. Not SKB264. GSE50927 is
mouse whole lung (ventilator-injury series), not lung cancer.
"""
from __future__ import annotations

import gzip
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

DATA = ROOT / "methods" / "cldn4_ko_gsea" / "data"
OUT = ROOT / "methods" / "cldn4_ko_gsea"
FIG = OUT / "figures"
TAB = OUT / "tables"
SET_JSON = ROOT / "data" / "genesets" / "a8_sets.json"
MM_GMT = ROOT / "data" / "genesets" / "mh.all.v2023.2.Mm.symbols.gmt"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "KEGG_TIGHT_JUNCTION",
    "GOBP_KERATINIZATION",
]
HEADLINE_LABEL = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Hallmark IFN-α",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
    "KEGG_TIGHT_JUNCTION": "KEGG tight junction",
    "GOBP_KERATINIZATION": "GO keratinization",
}
IFN_MHC = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
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
FPKM_PSEUDO = 0.5
ALIAS_22493 = {"G1P2": "ISG15", "G1P3": "IFI6"}


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


def fmt_num(x, nd=3) -> str:
    if x is None or x != x:
        return "NA"
    return f"{float(x):.{nd}f}"


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


def collapse_max_mean(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    df["_m"] = df[cols].mean(axis=1)
    return (
        df.sort_values("_m", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")[cols]
    )


def load_human_mouse_sets() -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    payload = json.loads(SET_JSON.read_text())
    human = {k: list(payload["sets"][k]) for k in HEADLINE}
    mm_hallmark = read_gmt(MM_GMT)
    mouse = {}
    for k in HEADLINE:
        if k in mm_hallmark:
            mouse[k] = mm_hallmark[k]
        else:
            mouse[k] = to_mouse(human[k])
    return human, mouse


def load_gse207704() -> dict[str, object]:
    df = pd.read_csv(DATA / "GSE207704_CLDN4_RNAseq.txt.gz", sep="\t", low_memory=False)
    cols = {
        "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
        "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
        "T47D_WT": "T47D_WT_FPKM (fpkm)",
    }
    df = df.rename(columns={"gene_short_name": "symbol"})
    df = df[df["symbol"].notna() & (df["symbol"].astype(str).str.len() > 0)]
    for c in cols.values():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    fpkm = collapse_max_mean(df[["symbol", *cols.values()]], list(cols.values()))
    fpkm = fpkm.rename(columns={v: k for k, v in cols.items()})
    lg = np.log2(fpkm + FPKM_PSEUDO)
    ranks = {
        "GSE207704_T47D": (lg["T47D_KO"] - lg["T47D_WT"]).dropna(),
        "GSE207704_MCF7": (lg["MCF7_KO"] - lg["MCF7_WT"]).dropna(),
    }
    ranks["GSE207704_mean_both_lines"] = pd.concat(
        [ranks["GSE207704_T47D"], ranks["GSE207704_MCF7"]], axis=1
    ).mean(axis=1)
    meta = {
        "GSE207704_T47D": {
            "accession": "GSE207704",
            "contrast": "T47D CLDN4-/- vs WT (breast RNA-seq, group-mean FPKM)",
            "species": "human",
            "tissue": "T47D breast cancer cell line (in vitro)",
            "n_loss": 2,
            "n_wt": 2,
            "n_note": "2 vs 2 biological replicates deposited; GEO file collapses them to one FPKM per group",
            "rank_metric": f"log2((KO+{FPKM_PSEUDO})/(WT+{FPKM_PSEUDO})) from group-mean FPKM",
            "perturb_gene": "CLDN4",
        },
        "GSE207704_MCF7": {
            "accession": "GSE207704",
            "contrast": "MCF7 CLDN4-/- vs WT (breast RNA-seq, group-mean FPKM)",
            "species": "human",
            "tissue": "MCF7 breast cancer cell line (in vitro)",
            "n_loss": 2,
            "n_wt": 2,
            "n_note": "2 vs 2 biological replicates deposited; GEO file collapses them to one FPKM per group",
            "rank_metric": f"log2((KO+{FPKM_PSEUDO})/(WT+{FPKM_PSEUDO})) from group-mean FPKM",
            "perturb_gene": "CLDN4",
        },
        "GSE207704_mean_both_lines": {
            "accession": "GSE207704",
            "contrast": "mean log2FC of T47D and MCF7 CLDN4-/- vs WT",
            "species": "human",
            "tissue": "T47D + MCF7 breast lines (descriptive mean rank)",
            "n_loss": 4,
            "n_wt": 4,
            "n_note": "mean of two group-mean log2FC ranks; not an independent sample",
            "rank_metric": "mean of T47D and MCF7 group-mean log2FC",
            "perturb_gene": "CLDN4",
        },
    }
    target = {}
    for key, line in [
        ("GSE207704_T47D", "T47D"),
        ("GSE207704_MCF7", "MCF7"),
        ("GSE207704_mean_both_lines", None),
    ]:
        if "CLDN4" not in fpkm.index:
            target[key] = {"symbol": "CLDN4", "log2FC": np.nan}
            continue
        if line is None:
            target[key] = {
                "symbol": "CLDN4",
                "log2FC": float(ranks[key].loc["CLDN4"]),
                "value_loss": float((fpkm.loc["CLDN4", "T47D_KO"] + fpkm.loc["CLDN4", "MCF7_KO"]) / 2),
                "value_wt": float((fpkm.loc["CLDN4", "T47D_WT"] + fpkm.loc["CLDN4", "MCF7_WT"]) / 2),
                "value_unit": "mean FPKM across lines",
            }
        else:
            target[key] = {
                "symbol": "CLDN4",
                "log2FC": float(ranks[key].loc["CLDN4"]),
                "value_loss": float(fpkm.loc["CLDN4", f"{line}_KO"]),
                "value_wt": float(fpkm.loc["CLDN4", f"{line}_WT"]),
                "value_unit": "group-mean FPKM",
            }
    return {"ranks": ranks, "meta": meta, "target": target, "fpkm": fpkm}


def load_gse50927() -> dict[str, object]:
    df = pd.read_csv(DATA / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    df = df[df["Marker.Symbol"].notna()]
    df = (
        df.assign(_abs=df["logFC"].abs())
        .sort_values(["logCPM", "_abs"], ascending=False)
        .drop_duplicates("Marker.Symbol")
        .set_index("Marker.Symbol")
    )
    rank = df["logFC"].astype(float).dropna()
    key = "GSE50927_naive_lung"
    meta = {
        key: {
            "accession": "GSE50927",
            "contrast": "Cldn4 KO vs WT, naive whole lung (no VILI)",
            "species": "mouse",
            "tissue": "mouse whole lung (not lung cancer; mixed-cell bulk)",
            "n_loss": 1,
            "n_wt": 1,
            "n_note": "n=1 vs 1; author edgeR table; no count matrix on GEO",
            "rank_metric": "author edgeR logFC (KO minus WT); Cldn4 logFC = -6.06 confirms sign",
            "perturb_gene": "Cldn4",
        }
    }
    target = {
        key: {
            "symbol": "Cldn4",
            "log2FC": float(df.loc["Cldn4", "logFC"]) if "Cldn4" in df.index else np.nan,
            "value_loss": np.nan,
            "value_wt": np.nan,
            "value_unit": "author edgeR logFC (no FPKM deposited)",
            "author_FDR": float(df.loc["Cldn4", "FDR"]) if "Cldn4" in df.index else np.nan,
            "author_logCPM": float(df.loc["Cldn4", "logCPM"]) if "Cldn4" in df.index else np.nan,
        }
    }
    return {"ranks": {key: rank}, "meta": meta, "target": target}


def _gse22493_symbol(orf: str, desc: str) -> str:
    s = (orf or "").strip().upper()
    if not s:
        s = (desc or "").split("--")[0].strip().upper()
    return ALIAS_22493.get(s, s)


def load_gse22493() -> dict[str, object]:
    plat = pd.read_csv(DATA / "GPL10555_platform.tsv.gz", sep="\t", dtype=str)
    plat["ID"] = plat["ID"].astype(str)
    rows, hdr = [], None
    with gzip.open(DATA / "GSE22493_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                hdr = next(fh).rstrip("\n").replace('"', "").split("\t")
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if hdr is not None:
                rows.append(line.rstrip("\n").replace('"', "").split("\t"))
    mat = pd.DataFrame(rows, columns=hdr)
    for c in mat.columns[1:]:
        mat[c] = pd.to_numeric(mat[c], errors="coerce")
    mat["ID_REF"] = mat["ID_REF"].astype(str)
    m = mat.merge(plat, left_on="ID_REF", right_on="ID", how="left")
    m["symbol"] = [
        _gse22493_symbol(o, d)
        for o, d in zip(m["ORF"].fillna(""), m["DESCRIPTION"].fillna(""))
    ]
    mapped = int((m["symbol"].str.len() > 0).sum())
    n_probe = len(m)
    gsms = ["GSM558700", "GSM558701", "GSM558702"]
    keep = m[m["symbol"].str.len() > 0].copy()
    gene = keep.groupby("symbol")[gsms].mean()
    rank = gene.mean(axis=1, skipna=True).dropna()
    key = "GSE22493_SKOV3"
    meta = {
        key: {
            "accession": "GSE22493",
            "contrast": "SKOV-3 CLDN4 siRNA vs CLDN4-overexpressing control",
            "species": "human",
            "tissue": "SKOV-3 ovarian line (in vitro two-colour array)",
            "n_loss": 3,
            "n_wt": 3,
            "n_note": "3 two-colour arrays; VALUE = author log2(KD / CLDN4-OE control), not parental WT",
            "rank_metric": "mean deposited log2(KD/OE) across 3 arrays after ORF/DESCRIPTION symbol map",
            "perturb_gene": "CLDN4",
            "n_probes": n_probe,
            "n_probes_mapped": mapped,
            "n_genes": int(rank.shape[0]),
            "tacstd2_on_platform": bool("TACSTD2" in rank.index),
        }
    }
    cldn = {
        "symbol": "CLDN4",
        "log2FC": float(rank.loc["CLDN4"]) if "CLDN4" in rank.index else np.nan,
        "value_unit": "mean log2(KD/OE) across arrays",
        "per_array": {},
    }
    if "CLDN4" in gene.index:
        for a in gsms:
            cldn["per_array"][a] = float(gene.loc["CLDN4", a]) if pd.notna(gene.loc["CLDN4", a]) else np.nan
    return {
        "ranks": {key: rank},
        "meta": meta,
        "target": {key: cldn},
        "mapping": {
            "included": True,
            "n_probes": n_probe,
            "n_probes_mapped": mapped,
            "n_genes": int(rank.shape[0]),
            "method": "GPL10555 ORF, else DESCRIPTION prefix before '--'",
            "tacstd2_on_platform": bool("TACSTD2" in rank.index),
        },
    }


def run_one(rank: pd.Series, sets: dict[str, list[str]], key: str, meta: dict) -> pd.DataFrame:
    rank = rank.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    g = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
    if g.empty:
        g = pd.DataFrame(columns=["term", "es", "nes", "nom_p", "n_set_in_rank"])
    g = g.copy()
    g["contrast"] = key
    g["accession"] = meta["accession"]
    g["species"] = meta["species"]
    g["n_loss"] = meta["n_loss"]
    g["n_wt"] = meta["n_wt"]
    g["n_genes_ranked"] = int(rank.shape[0])
    g["headline"] = g["term"].isin(HEADLINE)
    g["fdr"] = np.nan
    if not g.empty:
        g["fdr"] = bh_fdr(g["nom_p"])
    g["ifn_mhc_up"] = np.where(
        g["term"].isin(IFN_MHC),
        np.where(g["nes"] > 0, "UP", np.where(g["nes"] < 0, "DOWN", "flat")),
        "",
    )
    return g


def write_figures(gsea: pd.DataFrame, targets: pd.DataFrame) -> list[str]:
    written = []
    order_contrast = [
        "GSE207704_T47D",
        "GSE207704_MCF7",
        "GSE207704_mean_both_lines",
        "GSE50927_naive_lung",
        "GSE22493_SKOV3",
    ]
    short = {
        "GSE207704_T47D": "T47D KO",
        "GSE207704_MCF7": "MCF7 KO",
        "GSE207704_mean_both_lines": "Breast mean",
        "GSE50927_naive_lung": "Mouse lung KO",
        "GSE22493_SKOV3": "SKOV-3 siRNA",
    }
    prim = gsea[gsea["term"].isin(HEADLINE)].copy()
    prim["contrast"] = pd.Categorical(prim["contrast"], order_contrast, ordered=True)
    prim["term"] = pd.Categorical(prim["term"], HEADLINE, ordered=True)
    prim = prim.sort_values(["contrast", "term"])

    # Figure 1 — headline NES bars
    fig, ax = plt.subplots(figsize=(11.2, 5.6))
    x = np.arange(len(HEADLINE))
    n_c = len(order_contrast)
    width = 0.15
    colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#B279A2"]
    for i, c in enumerate(order_contrast):
        sub = prim[prim["contrast"] == c].set_index("term").reindex(HEADLINE)
        nes = sub["nes"].astype(float).to_numpy()
        ax.bar(x + (i - (n_c - 1) / 2) * width, nes, width=width, color=colors[i], label=short[c], zorder=3)
        for j, (nv, q) in enumerate(zip(nes, sub["fdr"].to_numpy())):
            if pd.notna(q) and q < 0.05 and pd.notna(nv):
                ax.text(x[j] + (i - (n_c - 1) / 2) * width, nv + (0.08 if nv >= 0 else -0.08),
                        "*", ha="center", va="bottom" if nv >= 0 else "top", fontsize=9, color=colors[i])
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([HEADLINE_LABEL[t] for t in HEADLINE], rotation=18, ha="right")
    ax.set_ylabel("NES (positive = up after CLDN4 / Cldn4 loss)")
    ax.set_title("CLDN4-loss prerank GSEA — headline sets\n* = BH-FDR < 0.05 within the 5 headline sets; gene-set permutation, seed=42")
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper right")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()
    p1 = FIG / "fig_headline_nes.png"
    fig.savefig(p1, dpi=160)
    fig.savefig(FIG / "fig_headline_nes.pdf")
    plt.close(fig)
    written.append(str(p1.relative_to(ROOT)))

    # Figure 2 — CLDN4 log2FC + IFN/MHC direction
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    ax = axes[0]
    t = targets.set_index("contrast").reindex(order_contrast)
    vals = t["log2FC"].astype(float)
    ax.barh(np.arange(len(order_contrast)), vals.fillna(0),
            color=["#E45756" if v < 0 else "#4C78A8" for v in vals.fillna(0)])
    ax.set_yticks(np.arange(len(order_contrast)))
    ax.set_yticklabels([short[c] for c in order_contrast])
    ax.axvline(0, color="k", lw=0.7)
    ax.invert_yaxis()
    ax.set_xlabel("CLDN4 / Cldn4 log2FC (loss vs WT / OE)")
    ax.set_title("Target transcript after CLDN4 / Cldn4 loss")
    for i, (c, v) in enumerate(vals.items()):
        if pd.notna(v):
            nL, nW = int(t.loc[c, "n_loss"]), int(t.loc[c, "n_wt"])
            ax.text(v + (0.12 if v >= 0 else -0.12), i, f"{v:+.2f}  n={nL} vs {nW}",
                    va="center", ha="left" if v >= 0 else "right", fontsize=8)

    ax = axes[1]
    heat = prim.pivot(index="term", columns="contrast", values="nes").reindex(index=HEADLINE, columns=order_contrast)
    qmat = prim.pivot(index="term", columns="contrast", values="fdr").reindex(index=HEADLINE, columns=order_contrast)
    im = ax.imshow(heat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
    ax.set_xticks(range(len(order_contrast)))
    ax.set_xticklabels([short[c] for c in order_contrast], rotation=25, ha="right")
    ax.set_yticks(range(len(HEADLINE)))
    ax.set_yticklabels([HEADLINE_LABEL[t] for t in HEADLINE])
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            v = heat.iloc[i, j]
            q = qmat.iloc[i, j]
            if pd.isna(v):
                txt = "NA"
            else:
                star = "*" if pd.notna(q) and q < 0.05 else ""
                txt = f"{v:+.2f}{star}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5,
                    color="white" if pd.notna(v) and abs(v) > 1.2 else "black")
    ax.set_title("Headline NES (loss-high rank)\n* FDR<0.05")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="NES")
    fig.tight_layout()
    p2 = FIG / "fig_cldn4_log2fc_and_nes_heatmap.png"
    fig.savefig(p2, dpi=160)
    fig.savefig(FIG / "fig_cldn4_log2fc_and_nes_heatmap.pdf")
    plt.close(fig)
    written.append(str(p2.relative_to(ROOT)))
    return written


def ifn_mhc_call(rows: pd.DataFrame) -> str:
    sub = rows[rows["term"].isin(IFN_MHC)]
    if sub.empty:
        return "not tested"
    ups = int((sub["nes"] > 0).sum())
    dns = int((sub["nes"] < 0).sum())
    sig_up = int(((sub["nes"] > 0) & (sub["fdr"] < 0.05)).sum())
    sig_dn = int(((sub["nes"] < 0) & (sub["fdr"] < 0.05)).sum())
    names_up = [HEADLINE_LABEL[t] for t in sub.loc[sub["nes"] > 0, "term"]]
    names_dn = [HEADLINE_LABEL[t] for t in sub.loc[sub["nes"] < 0, "term"]]
    bits = []
    if names_up:
        bits.append("UP: " + ", ".join(names_up) + (f" ({sig_up} FDR<0.05)" if sig_up else " (FDR≥0.05)"))
    if names_dn:
        bits.append("DOWN: " + ", ".join(names_dn) + (f" ({sig_dn} FDR<0.05)" if sig_dn else " (FDR≥0.05)"))
    return "; ".join(bits) if bits else "flat"


def write_finding(gsea: pd.DataFrame, targets: pd.DataFrame, mapping_22493: dict, figures: list[str]) -> None:
    def row(contrast: str, term: str) -> pd.Series:
        hit = gsea[(gsea["contrast"] == contrast) & (gsea["term"] == term)]
        if hit.empty:
            return pd.Series({"nes": np.nan, "fdr": np.nan, "nom_p": np.nan, "n_set_in_rank": 0, "es": np.nan})
        return hit.iloc[0]

    def nes_line(contrast: str) -> str:
        parts = []
        for t in HEADLINE:
            r = row(contrast, t)
            parts.append(
                f"| {HEADLINE_LABEL[t]} | `{t}` | {fmt_nes(r['nes'])} | {fmt_p(r['fdr'])} | {fmt_p(r['nom_p'])} | {int(r['n_set_in_rank']) if pd.notna(r['n_set_in_rank']) else 0} |"
            )
        return "\n".join(parts)

    def tgt(contrast: str) -> pd.Series:
        hit = targets[targets["contrast"] == contrast]
        return hit.iloc[0] if not hit.empty else pd.Series()

    t47 = tgt("GSE207704_T47D")
    mcf = tgt("GSE207704_MCF7")
    mn = tgt("GSE207704_mean_both_lines")
    lu = tgt("GSE50927_naive_lung")
    sk = tgt("GSE22493_SKOV3")

    def call(contrast: str) -> str:
        return ifn_mhc_call(gsea[gsea["contrast"] == contrast])

    per_array = sk.get("per_array", {}) if isinstance(sk.get("per_array"), dict) else {}
    arr_txt = ", ".join(
        f"{k}={v:+.3f}" if pd.notna(v) else f"{k}=NA" for k, v in per_array.items()
    ) if per_array else "see table"

    md = f"""# FINDING — public CLDN4 / Cldn4-loss prerank GSEA

**Additive only.** User thesis is taken as given: CLDN4 loss shares the IFN / MHC-I axis; CLDN4-high is the tight-junction barrier. This folder does **not** audit or retract any slide. It is **not** SKB264 and it is **not** a lung-cancer KO (GSE50927 is mouse **whole lung**).

Prerank engine is the same as `methods/scrna_pseudobulk_gsea_meta` (`scripts/cldn4_ko_gsea/gsea_core.py`): weighted KS *p*=1, **{NPERM}** gene-set permutations, seed=**{SEED}**. Positive NES = the set is enriched at the **CLDN4-loss** end of the rank (KO / siRNA minus WT / CLDN4-high). BH-FDR is within the **five headline sets** per contrast. Numbers below are written from `tables/gsea_headline.tsv` and `tables/cldn4_log2fc.tsv`.

---

## 一句话 / TL;DR

| Contrast | n (loss vs WT) | CLDN4 / Cldn4 log2FC | IFN / MHC-I after CLDN4 loss |
|---|---|---|---|
| GSE207704 T47D CLDN4-/- | 2 vs 2 (group-mean FPKM) | {fmt_fc(t47.get("log2FC"))} | {call("GSE207704_T47D")} |
| GSE207704 MCF7 CLDN4-/- | 2 vs 2 (group-mean FPKM) | {fmt_fc(mcf.get("log2FC"))} | {call("GSE207704_MCF7")} |
| GSE207704 mean of both lines | descriptive mean rank | {fmt_fc(mn.get("log2FC"))} | {call("GSE207704_mean_both_lines")} |
| GSE50927 Cldn4 KO naive **whole lung** | **1 vs 1** | {fmt_fc(lu.get("log2FC"))} | {call("GSE50927_naive_lung")} |
| GSE22493 SKOV-3 CLDN4 siRNA vs OE | 3 arrays | {fmt_fc(sk.get("log2FC"))} | {call("GSE22493_SKOV3")} |

GSE22493 **was included**: GPL10555 ORF (else DESCRIPTION prefix) maps **{mapping_22493["n_probes_mapped"]} / {mapping_22493["n_probes"]}** probes to **{mapping_22493["n_genes"]}** gene symbols. TACSTD2 is **not** on the platform.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrices | Public processed files only. No FASTQ / SRA. |
| Engine | Same prerank as `gsea_core.py` in `scrna_pseudobulk_gsea_meta` |
| Rank | log2FC of CLDN4-loss minus WT (or CLDN4-OE). Positive = up after loss. |
| Headline sets | Hallmark IFN-γ, Hallmark IFN-α, custom MHC-I/APM (21 genes), KEGG tight junction, GO keratinization |
| Mouse Hallmark IFN | MSigDB `mh.all.v2023.2.Mm` (not title-cased human lists) |
| Mouse MHC-I / TJ / keratin | human set → mouse symbols (HLA → classical H2; ERAP2 → Lnpep) |
| FDR | BH inside the 5 headline sets, per contrast |
| GSE207704 | cuffdiff **group-mean FPKM** (replicates already collapsed on GEO) |
| GSE50927 | author edgeR table, **naive whole lung**, not a tumour |
| GSE22493 | series-matrix VALUE; control arm is CLDN4-**overexpressing**, not parental |

Primary genes: CLDN4 `ENSG00000189143`, mouse Cldn4 `ENSMUSG00000047501`. TACSTD2 `ENSG00000184292` is recorded when present; it is not a GSEA-set member.

---

## 1. GSE207704 — CLDN4 KO T47D and MCF7 (breast RNA-seq)

GEO: [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704). File: `GSE207704_CLDN4_RNAseq.txt.gz` (FPKM). Design on the record: T47D WT GSM6310640–41, T47D CLDN4-/- GSM6310642–43, MCF7 WT GSM6310644–45, MCF7 CLDN4-/- GSM6310646–47 (2 vs 2 per line). The open processed table has **one FPKM column per group**, so replicate-level *t* statistics cannot be computed.

CLDN4 itself (highest-FPKM locus after symbol collapse; pseudocount {FPKM_PSEUDO}):

| Line | FPKM WT | FPKM KO | log2FC KO vs WT | n |
|---|---:|---:|---:|---|
| T47D | {fmt_num(t47.get("value_wt"))} | {fmt_num(t47.get("value_loss"))} | {fmt_fc(t47.get("log2FC"))} | 2 vs 2, collapsed |
| MCF7 | {fmt_num(mcf.get("value_wt"))} | {fmt_num(mcf.get("value_loss"))} | {fmt_fc(mcf.get("log2FC"))} | 2 vs 2, collapsed |

### T47D CLDN4-/- vs WT

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line("GSE207704_T47D")}

IFN / MHC-I after CLDN4 loss: **{call("GSE207704_T47D")}**.

### MCF7 CLDN4-/- vs WT

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line("GSE207704_MCF7")}

IFN / MHC-I after CLDN4 loss: **{call("GSE207704_MCF7")}**.

### Mean of both lines (descriptive)

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line("GSE207704_mean_both_lines")}

MHC-I / APM on this cuffdiff table is sparse (classical HLA-A/B, B2M, TAP1/2, PSMB8/9 are absent as `gene_short_name`). `n_set_in_rank` is the number that were actually ranked.

---

## 2. GSE50927 — Cldn4 KO mouse **whole lung** (not lung cancer)

GEO: [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). Author edgeR table `GSE50927_Cldn4lungWTvsKOgenes.csv.gz`. Naive (no ventilator) Cldn4 KO vs WT. **n = 1 vs 1**. Mixed-cell whole lung on a mixed 129S6/C57BL/6/BALB/c background. This is **not** a lung-tumour series.

Cldn4 author logFC = **{fmt_fc(lu.get("log2FC"))}** (logCPM {fmt_num(lu.get("author_logCPM"))}; submitter edgeR FDR {fmt_p(lu.get("author_FDR"))}). Sign is KO minus WT.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line("GSE50927_naive_lung")}

IFN / MHC-I after Cldn4 loss: **{call("GSE50927_naive_lung")}**.

Author edgeR *P* / FDR on an unreplicated design assume a dispersion; they are **not** used as the GSEA null. The rank is the deposited logFC.

---

## 3. GSE22493 — SKOV-3 CLDN4 siRNA (mapped, included)

GEO: [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493), platform [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555). Series-matrix VALUE columns GSM558700–702.

Mapping: **{mapping_22493["method"]}**. **{mapping_22493["n_probes_mapped"]} / {mapping_22493["n_probes"]}** probes → **{mapping_22493["n_genes"]}** symbols. TACSTD2 on platform: **{mapping_22493["tacstd2_on_platform"]}**.

CLDN4 mean log2(KD / OE) = **{fmt_fc(sk.get("log2FC"))}** (per array: {arr_txt}). Control arm is CLDN4-overexpressing SKOV-3, not parental WT. n = 3 arrays.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
{nes_line("GSE22493_SKOV3")}

IFN / MHC-I after CLDN4 loss: **{call("GSE22493_SKOV3")}**.

---

## Extra figures

- `{figures[0]}`
- `{figures[1]}`

Tables: `methods/cldn4_ko_gsea/tables/gsea_headline.tsv`, `gsea_prerank_all.tsv`, `cldn4_log2fc.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-run of any existing slide or of SKB264.
- Not FASTQ / salmon / DESeq2 from SRA.
- Not a lung-cancer Cldn4 KO. GSE50927 is whole lung.
- Not sample-permutation GSEA. n on GSE207704 and GSE50927 is too small for that; the engine is gene-set permutation on a prerank, same as the scRNA pseudobulk extra.

---

## 中文摘要

只补公开 CLDN4 / Cldn4 缺失的 prerank GSEA，不审不撤已有页。引擎与 `scrna_pseudobulk_gsea_meta` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 缺失端富集。

- **GSE207704** 乳腺 T47D / MCF7 CRISPR KO，GEO 只有组均 FPKM（2 vs 2 已合并）。CLDN4 log2FC T47D {fmt_fc(t47.get("log2FC"))}、MCF7 {fmt_fc(mcf.get("log2FC"))}。IFN/MHC：T47D {call("GSE207704_T47D")}；MCF7 {call("GSE207704_MCF7")}。
- **GSE50927** 是小鼠**全肺** Cldn4 KO，不是肺癌。n=1 vs 1。Cldn4 logFC {fmt_fc(lu.get("log2FC"))}。IFN/MHC：{call("GSE50927_naive_lung")}。
- **GSE22493** SKOV-3 siRNA 可用 GPL10555 映射到基因符号（{mapping_22493["n_genes"]} 个），已纳入。对照是 CLDN4 过表达而非亲本。CLDN4 log2(KD/OE) {fmt_fc(sk.get("log2FC"))}。IFN/MHC：{call("GSE22493_SKOV3")}。TACSTD2 不在芯片上。
"""
    (OUT / "FINDING.md").write_text(md)
    log(f"wrote {OUT / 'FINDING.md'}")


def main() -> None:
    human_sets, mouse_sets = load_human_mouse_sets()
    for name, sets in [("human", human_sets), ("mouse", mouse_sets)]:
        for t, genes in sets.items():
            log(f"set {name} {t}: {len(genes)} genes")

    d704 = load_gse207704()
    d509 = load_gse50927()
    d224 = load_gse22493()
    log(f"GSE22493 mapped {d224['mapping']['n_probes_mapped']}/{d224['mapping']['n_probes']} probes -> {d224['mapping']['n_genes']} symbols")

    packs = [
        (d704, human_sets),
        (d509, mouse_sets),
        (d224, human_sets),
    ]
    gsea_rows = []
    target_rows = []
    inv_rows = []
    for pack, sets in packs:
        for key, rank in pack["ranks"].items():
            meta = pack["meta"][key]
            log(f"GSEA {key}: {len(rank)} ranked genes, n={meta['n_loss']} vs {meta['n_wt']}")
            g = run_one(rank, sets, key, meta)
            gsea_rows.append(g)
            tgt = dict(pack["target"][key])
            per_array = tgt.pop("per_array", None)
            tgt["contrast"] = key
            tgt["accession"] = meta["accession"]
            tgt["n_loss"] = meta["n_loss"]
            tgt["n_wt"] = meta["n_wt"]
            tgt["n_note"] = meta["n_note"]
            tgt["tissue"] = meta["tissue"]
            if per_array:
                tgt["per_array"] = per_array
            target_rows.append(tgt)
            inv_rows.append(
                {
                    "contrast": key,
                    **{k: meta[k] for k in ("accession", "contrast", "species", "tissue", "n_loss", "n_wt", "n_note", "rank_metric", "perturb_gene")},
                    "n_genes_ranked": int(rank.shape[0]),
                }
            )

    gsea = pd.concat(gsea_rows, ignore_index=True)
    headline = gsea[gsea["term"].isin(HEADLINE)].copy()
    targets = pd.DataFrame(target_rows)
    # keep per_array as JSON string for the TSV
    if "per_array" in targets.columns:
        targets["per_array"] = targets["per_array"].apply(
            lambda x: json.dumps(x) if isinstance(x, dict) else x
        )
    inv = pd.DataFrame(inv_rows)

    gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
    headline.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    targets.to_csv(TAB / "cldn4_log2fc.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "inventory.tsv", sep="\t", index=False)
    log(f"wrote tables under {TAB}")

    # restore per_array dict for FINDING
    if "per_array" in targets.columns:
        targets["per_array"] = targets["per_array"].apply(
            lambda x: json.loads(x) if isinstance(x, str) and x.startswith("{") else (x if isinstance(x, dict) else {})
        )

    figures = write_figures(gsea, targets)
    log("figures: " + ", ".join(figures))
    write_finding(gsea, targets, d224["mapping"], figures)

    print("\n=== HEADLINE NES ===")
    show = headline[["contrast", "term", "nes", "fdr", "nom_p", "n_set_in_rank"]].copy()
    show["nes"] = show["nes"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
    show["fdr"] = show["fdr"].map(fmt_p)
    print(show.to_string(index=False))
    print("\n=== CLDN4 / Cldn4 ===")
    print(targets[["contrast", "symbol", "log2FC", "n_loss", "n_wt"]].to_string(index=False))


if __name__ == "__main__":
    main()
