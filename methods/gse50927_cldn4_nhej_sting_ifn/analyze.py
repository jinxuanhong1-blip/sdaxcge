#!/usr/bin/env python3
"""GSE50927 mouse VILI — Cldn4 NHEJ / STING / IFN gene panels.

Additive public mouse. Whole-lung ventilator injury, not a tumour series.
Public processed EdgeR tables only (no SRA / FASTQ). Cldn4-high = WT
(Cldn4-intact). Cldn4-low = Cldn4 KO. KOlow / KOhigh are BAL-protein injury
strata, not Cldn4-expression strata.

Author logFC is the second group minus the first (KO minus WT, or VILI minus
naive). Positive NES on a genotype contrast means the set sits at the
Cldn4-loss end. Compact means are also reported as Cldn4-high minus Cldn4-low
(sign flip of the deposited logFC) so they read the same way as the earlier
IFN/MHC/TJ folder.

Gene-set permutation and a Wilcoxon on deposited logFC are not mouse-level
tests. Honest genotype n is 1 GSM vs 1 GSM.
"""
from __future__ import annotations

import json
import os
import time
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("GSE50927_NHEJ_DATA", "/tmp/gse50927_nhej"))
TAB = HERE / "tables"
FIG = HERE / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927"
FILES = {
    "naive": f"{FTP}/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "vili_wt": f"{FTP}/suppl/GSE50927_VILIwtGenes.csv.gz",
    "vili_kohi": f"{FTP}/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "vili_kolo": f"{FTP}/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
}

CONTRASTS = {
    "naive_KO_vs_WT": {
        "file": "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        "key": "naive",
        "role": "primary",
        "label": "Cldn4 KO vs WT, naive whole lung (no VILI)",
        "logfc_is": "KO minus WT",
        "flip_to_high_minus_low": True,
        "n_high": 1,
        "n_low": 1,
    },
    "VILI_WT_vs_naive_WT": {
        "file": "GSE50927_VILIwtGenes.csv.gz",
        "key": "vili_wt",
        "role": "companion_induction",
        "label": "WT VILI vs WT no VILI (Cldn4 induction, not a genotype split)",
        "logfc_is": "VILI minus naive (both WT)",
        "flip_to_high_minus_low": False,
        "n_high": 1,
        "n_low": 1,
    },
    "VILI_KOhigh_vs_WT": {
        "file": "GSE50927_VILIwtkohiGenes.csv.gz",
        "key": "vili_kohi",
        "role": "companion_injury",
        "label": "Cldn4 KO VILIhigh vs WT VILI",
        "logfc_is": "KO VILIhigh minus WT VILI",
        "flip_to_high_minus_low": True,
        "n_high": 1,
        "n_low": 1,
    },
    "VILI_KOlow_vs_WT": {
        "file": "GSE50927_VILIwtkoloGenes.csv.gz",
        "key": "vili_kolo",
        "role": "companion_injury",
        "label": "Cldn4 KO VILIlow vs WT VILI",
        "logfc_is": "KO VILIlow minus WT VILI",
        "flip_to_high_minus_low": True,
        "n_high": 1,
        "n_low": 1,
    },
}

# Human or uppercase symbols that do not mouse-case onto this mm9 table.
ALIAS = {
    "MRE11": "Mre11a",
    "MRE11A": "Mre11a",
    "TP53BP1": "Trp53bp1",
    "H2AX": "H2afx",
    "CGAS": "Mb21d1",
    "MB21D1": "Mb21d1",
    "STING1": "Tmem173",
    "TMEM173": "Tmem173",
}

# Official Reactome STING (R-HSA-1834941) includes the DNA-PK complex.
# Those three genes are scored in NHEJ, so a second set drops them.
NHEJ_IN_STING = {"PRKDC", "XRCC5", "XRCC6"}

GSEA_SETS = [
    "NHEJ_CORE",
    "KEGG_NHEJ_mmu03450",
    "GO_NHEJ_0006303",
    "REACTOME_NHEJ_R-HSA-5693571",
    "STING_CORE",
    "REACTOME_STING_R-HSA-1834941",
    "REACTOME_STING_NO_NHEJ",
    "REACTOME_CYTOSOLIC_DNA_R-HSA-3134975",
    "KEGG_CYTOSOLIC_DNA_SENSING_mmu04623",
    "IFN_COMPACT",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

COMPACT_SETS = ["NHEJ_CORE", "STING_CORE", "IFN_COMPACT"]

FOCAL = [
    "Cldn4",
    "Xrcc6", "Xrcc5", "Prkdc", "Dclre1c", "Lig4", "Xrcc4", "Nhej1", "Poll", "Polm",
    "Fen1", "Dntt", "Mre11a", "Rad50", "Trp53bp1",
    "Mb21d1", "Tmem173", "Tbk1", "Ikbke", "Irf3", "Ddx41", "Trex1",
    "Ifnb1", "Stat1", "Irf7", "Isg15", "Cxcl10", "Cxcl9", "Ifng",
    "Tnf", "Il1b", "Egr1",
]

SET_LABEL = {
    "NHEJ_CORE": "NHEJ core",
    "KEGG_NHEJ_mmu03450": "KEGG NHEJ",
    "GO_NHEJ_0006303": "GO NHEJ",
    "REACTOME_NHEJ_R-HSA-5693571": "Reactome NHEJ",
    "STING_CORE": "STING core",
    "REACTOME_STING_R-HSA-1834941": "Reactome STING",
    "REACTOME_STING_NO_NHEJ": "Reactome STING minus DNA-PK",
    "REACTOME_CYTOSOLIC_DNA_R-HSA-3134975": "Reactome cytosolic DNA",
    "KEGG_CYTOSOLIC_DNA_SENSING_mmu04623": "KEGG cytosolic DNA-sensing",
    "IFN_COMPACT": "IFN compact",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Hallmark IFN-α",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    log(f"download {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def load_edger(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    scol = "Marker.Symbol" if "Marker.Symbol" in df.columns else "GeneSymbol"
    df = df.rename(columns={scol: "symbol"})
    df["symbol"] = df["symbol"].astype(str)
    df = df[df["symbol"].notna() & (df["symbol"] != "") & (df["symbol"] != "nan")]
    for c in ("logFC", "logCPM", "PValue", "FDR"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.sort_values("logCPM", ascending=False).drop_duplicates("symbol")
    return df.set_index("symbol", drop=False)


def to_mouse(symbol: str, universe: set[str]) -> str | None:
    if symbol in universe:
        return symbol
    key = symbol.upper()
    if key in ALIAS and ALIAS[key] in universe:
        return ALIAS[key]
    if not symbol:
        return None
    mapped = symbol[0].upper() + symbol[1:].lower()
    if mapped in universe:
        return mapped
    return None


def map_set(genes: list[str], universe: set[str]) -> tuple[list[str], list[str]]:
    used: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for g in genes:
        m = to_mouse(g, universe)
        if m is None:
            missing.append(g)
            continue
        if m not in seen:
            seen.add(m)
            used.append(m)
    return used, missing


def load_sets(universe: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    raw = json.loads((HERE / "genesets" / "panels.json").read_text())
    raw["REACTOME_STING_NO_NHEJ"] = [
        g for g in raw["REACTOME_STING_R-HSA-1834941"] if g.upper() not in NHEJ_IN_STING
    ]
    sets: dict[str, list[str]] = {}
    rows = []
    for name in GSEA_SETS:
        used, missing = map_set(raw[name], universe)
        sets[name] = used
        rows.append(
            {
                "set": name,
                "label": SET_LABEL[name],
                "n_listed": len(raw[name]),
                "n_in_rank": len(used),
                "n_missing": len(missing),
                "missing": ",".join(missing),
                "symbols_used": ",".join(used) if name in COMPACT_SETS or len(used) <= 20 else "",
            }
        )
    return sets, pd.DataFrame(rows)


def compact_score(df: pd.DataFrame, genes: list[str], flip: bool, min_cpm: float = 0.0) -> dict:
    """Mean of deposited logFC. Genes with logCPM < 0 are dropped when at least
    five genes remain, matching the earlier GSE50927 IFN compact (Ifng / Cxcl11).
    """
    present = [g for g in genes if g in df.index]
    sub = df.loc[present]
    detectable = sub[sub["logCPM"] >= min_cpm]
    used = detectable if len(detectable) >= 5 else sub
    fc = used["logFC"].astype(float)
    signed = -fc if flip else fc
    bg_fc = df.loc[~df.index.isin(used.index), "logFC"].astype(float)
    bg = -bg_fc if flip else bg_fc
    w_stat = w_p = mwu_p = np.nan
    if len(signed) >= 6 and np.any(signed != 0):
        try:
            w_stat, w_p = wilcoxon(signed.to_numpy(), zero_method="wilcox", alternative="two-sided")
        except ValueError:
            pass
    if len(signed) >= 5 and len(bg) >= 20:
        _, mwu_p = mannwhitneyu(signed.to_numpy(), bg.to_numpy(), alternative="two-sided")
    return {
        "n_listed": len(genes),
        "n_present": int(len(sub)),
        "n_low_cpm_dropped": int(len(sub) - len(used)),
        "n_used": int(len(signed)),
        "mean_logfc_deposited": float(fc.mean()) if len(fc) else np.nan,
        "mean_logfc_high_minus_low": float(signed.mean()) if len(signed) else np.nan,
        "median_logfc_high_minus_low": float(signed.median()) if len(signed) else np.nan,
        "n_up_in_high": int((signed > 0).sum()),
        "n_down_in_high": int((signed < 0).sum()),
        "wilcoxon_stat": float(w_stat) if pd.notna(w_stat) else np.nan,
        "wilcoxon_p": float(w_p) if pd.notna(w_p) else np.nan,
        "mwu_vs_background_p": float(mwu_p) if pd.notna(mwu_p) else np.nan,
    }


def one_gene(df: pd.DataFrame, symbol: str) -> dict:
    if symbol not in df.index:
        return {"symbol": symbol, "logFC": np.nan, "logCPM": np.nan, "FDR": np.nan}
    r = df.loc[symbol]
    return {
        "symbol": symbol,
        "logFC": float(r["logFC"]),
        "logCPM": float(r["logCPM"]),
        "FDR": float(r["FDR"]),
    }


def plot_compact(compact: pd.DataFrame) -> None:
    order = list(CONTRASTS)
    panels = COMPACT_SETS
    colors = {"NHEJ_CORE": "#3C6E9F", "STING_CORE": "#D4762C", "IFN_COMPACT": "#6B4C9A"}
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    x = np.arange(len(order))
    width = 0.24
    for i, panel in enumerate(panels):
        sub = compact[compact["set"] == panel].set_index("contrast").loc[order]
        ax.bar(
            x + (i - 1) * width,
            sub["mean_logfc_high_minus_low"],
            width=width,
            color=colors[panel],
            label=SET_LABEL[panel],
        )
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["naive KO vs WT\n(primary)", "WT VILI vs naive\n(induction)", "KO VILIhigh\nvs WT VILI", "KO VILIlow\nvs WT VILI"],
        fontsize=8,
    )
    ax.set_ylabel("Mean panel logFC\n(Cldn4-high − low, or VILI − naive)")
    ax.set_title("GSE50927 gene panels (1 vs 1 GSM; not a mouse-level test)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_compact_panels.png", dpi=160)
    fig.savefig(FIG / "fig1_compact_panels.pdf")
    plt.close(fig)


def plot_nes(gsea: pd.DataFrame) -> None:
    show = [
        "KEGG_NHEJ_mmu03450",
        "NHEJ_CORE",
        "STING_CORE",
        "REACTOME_STING_NO_NHEJ",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    ]
    order = list(CONTRASTS)
    mat = (
        gsea[gsea["term"].isin(show)]
        .pivot(index="term", columns="contrast", values="nes")
        .reindex(index=show, columns=order)
    )
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["naive\nKO vs WT", "WT VILI\nvs naive", "KO VILIhigh\nvs WT", "KO VILIlow\nvs WT"], fontsize=8)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels([SET_LABEL[s] for s in show], fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.to_numpy()[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="black")
    ax.set_title("Prerank NES on deposited logFC (positive = up in KO or in VILI)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="NES")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_gsea_nes.png", dpi=160)
    fig.savefig(FIG / "fig2_gsea_nes.pdf")
    plt.close(fig)


def plot_genes(genes: pd.DataFrame) -> None:
    show = [
        "Cldn4",
        "Xrcc6", "Xrcc5", "Prkdc", "Lig4", "Xrcc4", "Nhej1", "Dclre1c",
        "Mb21d1", "Tmem173", "Tbk1", "Irf3", "Trex1",
        "Ifnb1", "Stat1", "Isg15", "Cxcl10",
    ]
    order = list(CONTRASTS)
    mat = (
        genes[genes["symbol"].isin(show)]
        .pivot(index="symbol", columns="contrast", values="logFC")
        .reindex(index=show, columns=order)
    )
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    lim = np.nanmax(np.abs(mat.to_numpy()))
    lim = max(2.5, float(lim))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["naive\nKO−WT", "WT VILI\n− naive", "VILIhigh\nKO−WT", "VILIlow\nKO−WT"], fontsize=8)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show, fontsize=8)
    ax.set_title("Deposited logFC (KO − WT, except the VILI induction column)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="logFC")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_focal_logfc.png", dpi=160)
    fig.savefig(FIG / "fig3_focal_logfc.pdf")
    plt.close(fig)


def main() -> None:
    frames = {}
    for key, url in FILES.items():
        path = fetch(url, DATA / Path(url).name)
        frames[key] = load_edger(path)
        log(f"{key} genes {len(frames[key])}")

    universe = set(frames["naive"].index)
    sets, coverage = load_sets(universe)
    coverage.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)
    log("coverage\n" + coverage[["set", "n_listed", "n_in_rank", "n_missing"]].to_string(index=False))

    gsea_rows = []
    compact_rows = []
    gene_rows = []
    cldn_rows = []
    for name, meta in CONTRASTS.items():
        df = frames[meta["key"]]
        rank = df["logFC"].replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
        log(f"GSEA {name} n_rank={len(rank)}")
        res = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED, min_size=6)
        res["contrast"] = name
        res["role"] = meta["role"]
        res["fdr"] = bh_fdr(res["nom_p"])
        gsea_rows.append(res)
        for panel in COMPACT_SETS:
            row = compact_score(df, sets[panel], flip=meta["flip_to_high_minus_low"])
            row.update({"contrast": name, "role": meta["role"], "set": panel, "label": SET_LABEL[panel]})
            compact_rows.append(row)
        for symbol in FOCAL:
            row = one_gene(df, symbol)
            row.update({"contrast": name, "role": meta["role"]})
            gene_rows.append(row)
        cldn = one_gene(df, "Cldn4")
        cldn_rows.append({"contrast": name, "role": meta["role"], "logfc_is": meta["logfc_is"], **{f"Cldn4_{k}": v for k, v in cldn.items() if k != "symbol"}})

    gsea = pd.concat(gsea_rows, ignore_index=True)
    gsea["label"] = gsea["term"].map(SET_LABEL)
    compact = pd.DataFrame(compact_rows)
    genes = pd.DataFrame(gene_rows)
    gsea.to_csv(TAB / "gsea_panels.tsv", sep="\t", index=False)
    compact.to_csv(TAB / "compact_panel_scores.tsv", sep="\t", index=False)
    genes.to_csv(TAB / "gene_level.tsv", sep="\t", index=False)
    pd.DataFrame(cldn_rows).to_csv(TAB / "cldn4_logfc.tsv", sep="\t", index=False)

    # one row per contrast: Cldn4 plus the three compact means and headline NES
    one = []
    for name, meta in CONTRASTS.items():
        row = {
            "contrast": name,
            "role": meta["role"],
            "n_high_vs_low": f"{meta['n_high']} vs {meta['n_low']}",
            "logfc_is": meta["logfc_is"],
            "Cldn4_logFC_deposited": one_gene(frames[meta["key"]], "Cldn4")["logFC"],
        }
        for panel in COMPACT_SETS:
            sub = compact[(compact["contrast"] == name) & (compact["set"] == panel)].iloc[0]
            row[f"{panel}_mean_high_minus_low"] = sub["mean_logfc_high_minus_low"]
            row[f"{panel}_wilcoxon_p"] = sub["wilcoxon_p"]
            row[f"{panel}_n_used"] = sub["n_used"]
        for term in (
            "KEGG_NHEJ_mmu03450",
            "STING_CORE",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        ):
            sub = gsea[(gsea["contrast"] == name) & (gsea["term"] == term)]
            if sub.empty:
                continue
            row[f"{term}_NES"] = float(sub["nes"].iloc[0])
            row[f"{term}_FDR"] = float(sub["fdr"].iloc[0])
        one.append(row)
    one_df = pd.DataFrame(one)
    one_df.to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    plot_compact(compact)
    plot_nes(gsea)
    plot_genes(genes)

    summary = {
        "dataset": "GSE50927",
        "organism": "Mus musculus",
        "context": "whole-lung VILI, Cldn4 KO, not cancer",
        "n_gsm_per_arm": "1 vs 1",
        "nperm": NPERM,
        "seed": SEED,
        "contrasts": one,
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    log("wrote tables and figures")
    print(one_df.to_string(index=False))
    print(compact[["contrast", "set", "n_used", "mean_logfc_deposited", "mean_logfc_high_minus_low", "n_up_in_high", "n_down_in_high", "wilcoxon_p", "mwu_vs_background_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
