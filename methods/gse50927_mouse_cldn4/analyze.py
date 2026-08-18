#!/usr/bin/env python3
"""GSE50927 mouse whole-lung Cldn4 KO — Cldn4-high vs low IFN / MHC / TJ.

Additive public mouse. Cldn4-only. Injury, not tumour. Public processed
EdgeR tables only (no SRA / FASTQ). Thesis is taken as given: Cldn4-high
epithelium is barrier / IFN-low. This folder does not audit that as failed.

Cldn4-high = WT (Cldn4-intact). Cldn4-low = KO (Cldn4-null).
KOlow / KOhigh are injury-severity strata, not Cldn4-expression strata.
"""
from __future__ import annotations

import gzip
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
DATA = Path(os.environ.get("GSE50927_MOUSE_CLDN4_DATA", "/tmp/gse50927_mouse_cldn4"))
TAB = HERE / "tables"
FIG = HERE / "figures"
SETDIR = HERE / "genesets"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927"
FILES = {
    "matrix": f"{FTP}/matrix/GSE50927_series_matrix.txt.gz",
    "naive": f"{FTP}/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "vili_wt": f"{FTP}/suppl/GSE50927_VILIwtGenes.csv.gz",
    "vili_kohi": f"{FTP}/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "vili_kolo": f"{FTP}/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
}

# Author EdgeR logFC is second group minus first (KO minus WT, or VILI minus naive).
# Confirmed vs Kage et al. 2014 Table 3 (Ccl2 log2FC +3.0, FDR 3.6e-24 on KOhigh).
CONTRASTS = {
    "naive_KO_vs_WT": {
        "file": "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        "key": "naive",
        "role": "primary",
        "label": "Cldn4 KO vs WT, naive whole lung (no VILI)",
        "high_arm": "WT no VILI (Cldn4-intact)",
        "low_arm": "Cldn4 KO no VILI (Cldn4-null)",
        "logfc_is": "KO minus WT",
        "n_high": 1,
        "n_low": 1,
        "n_note": "1 GSM vs 1 GSM; design text says duplicates; 2 SRA runs/GSM; no count matrix",
        "is_genotype_split": True,
    },
    "VILI_WT_vs_naive_WT": {
        "file": "GSE50927_VILIwtGenes.csv.gz",
        "key": "vili_wt",
        "role": "companion_induction",
        "label": "WT VILI vs WT no VILI (Cldn4 induction, not high-vs-low)",
        "high_arm": "WT VILI",
        "low_arm": "WT no VILI",
        "logfc_is": "VILI minus naive (both WT)",
        "n_high": 1,
        "n_low": 1,
        "n_note": "1 GSM vs 1 GSM; not a Cldn4 genotype split",
        "is_genotype_split": False,
    },
    "VILI_KOhigh_vs_WT": {
        "file": "GSE50927_VILIwtkohiGenes.csv.gz",
        "key": "vili_kohi",
        "role": "companion_injury",
        "label": "Cldn4 KO VILIhigh vs WT VILI",
        "high_arm": "WT VILI (Cldn4-intact + injury)",
        "low_arm": "Cldn4 KO VILIhigh (Cldn4-null + higher BAL leak)",
        "logfc_is": "KO VILIhigh minus WT VILI",
        "n_high": 1,
        "n_low": 1,
        "n_note": "KOhigh = injury severity, not Cldn4 expression; 1 vs 1 GSM",
        "is_genotype_split": True,
    },
    "VILI_KOlow_vs_WT": {
        "file": "GSE50927_VILIwtkoloGenes.csv.gz",
        "key": "vili_kolo",
        "role": "companion_injury",
        "label": "Cldn4 KO VILIlow vs WT VILI",
        "high_arm": "WT VILI (Cldn4-intact + injury)",
        "low_arm": "Cldn4 KO VILIlow (Cldn4-null + WT-like BAL leak)",
        "logfc_is": "KO VILIlow minus WT VILI",
        "n_high": 1,
        "n_low": 1,
        "n_note": "KOlow = injury severity, not Cldn4 expression; 1 vs 1 GSM",
        "is_genotype_split": True,
    },
}

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

# Compact pre-specified epithelial / thesis panels (mouse symbols).
# Cldn4 is the splitter / KO gene and is excluded from TJ_COMPACT.
COMPACT = {
    "IFN_COMPACT": [
        "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15", "Mx1",
        "Oas1a", "Oas2", "Oasl1", "Rsad2", "Ifit1", "Ifit2", "Ifit3",
        "Ifi44", "Usp18", "Gbp2", "Gbp3", "Cxcl9", "Cxcl10", "Cxcl11",
        "Ccl5", "Ifng", "Eif2ak2", "Bst2", "Dhx58",
    ],
    "MHC_COMPACT": [
        "H2-K1", "H2-D1", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T23",
        "B2m", "Tap1", "Tap2", "Tapbp", "Nlrc5", "Ciita",
        "Psmb8", "Psmb9", "Psmb10", "Erap1", "Cd274",
    ],
    "TJ_COMPACT": [
        "Cldn3", "Cldn5", "Cldn7", "Cldn18", "Tjp1", "Tjp2",
        "Ocln", "Marveld2", "F11r", "Cdh1",
    ],
    "EPI_COMPACT": [
        "Epcam", "Krt8", "Krt18", "Krt19", "Cdh1", "Aqp5", "Hopx", "Sftpc",
    ],
}
FOCAL = [
    "Cldn4", "Cldn3", "Cldn18", "Cldn7", "Tjp1", "Ocln", "Cdh1",
    "Stat1", "Irf1", "Isg15", "Cxcl9", "Cxcl10", "Ifng",
    "B2m", "H2-K1", "H2-D1", "Nlrc5", "Cd274",
    "Tnf", "Il1b", "Il6", "Egr1", "Tacstd2",
]


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
    human = json.loads((SETDIR / "custom_human_sets.json").read_text())
    mm = read_gmt(SETDIR / "mh.all.v2023.2.Mm.symbols.gmt")
    mouse: dict[str, list[str]] = {}
    for k in HEADLINE:
        if k in mm:
            mouse[k] = mm[k]
        else:
            mouse[k] = to_mouse(human[k])
    mouse.update(COMPACT)
    return mouse


def load_edger(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    scol = "Marker.Symbol" if "Marker.Symbol" in df.columns else "GeneSymbol"
    df = df.rename(columns={scol: "symbol"})
    df["symbol"] = df["symbol"].astype(str)
    df = df[df["symbol"].notna() & (df["symbol"] != "") & (df["symbol"] != "nan")]
    for c in ("logFC", "logCPM", "PValue", "FDR"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.sort_values("logCPM", ascending=False).drop_duplicates("symbol")
    df = df.set_index("symbol", drop=False)
    return df


def rank_from_edger(df: pd.DataFrame) -> pd.Series:
    s = df["logFC"].replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def compact_score(df: pd.DataFrame, genes: list[str], min_cpm: float = 0.0) -> dict:
    present = [g for g in genes if g in df.index]
    sub = df.loc[present].copy()
    detectable = sub[sub["logCPM"] >= min_cpm]
    used = detectable if len(detectable) >= 5 else sub
    fc = used["logFC"].astype(float)
    bg = df.loc[~df.index.isin(used.index), "logFC"].astype(float)
    n_up = int((fc > 0).sum())
    n_dn = int((fc < 0).sum())
    w_stat = w_p = mwu_stat = mwu_p = np.nan
    if len(fc) >= 6 and np.any(fc != 0):
        try:
            w_stat, w_p = wilcoxon(fc.values, zero_method="wilcox", alternative="two-sided")
        except ValueError:
            pass
    if len(fc) >= 5 and len(bg) >= 20:
        mwu_stat, mwu_p = mannwhitneyu(fc.values, bg.values, alternative="two-sided")
    return {
        "n_listed": len(genes),
        "n_present": len(present),
        "n_used": int(len(used)),
        "n_low_cpm_dropped": int(len(sub) - len(used)),
        "mean_logFC_KO_minus_WT": float(fc.mean()) if len(fc) else np.nan,
        "median_logFC_KO_minus_WT": float(fc.median()) if len(fc) else np.nan,
        "mean_logFC_Cldn4high_minus_low": float(-fc.mean()) if len(fc) else np.nan,
        "n_up_in_KO": n_up,
        "n_down_in_KO": n_dn,
        "wilcoxon_stat": float(w_stat) if w_stat == w_stat else np.nan,
        "wilcoxon_p": float(w_p) if w_p == w_p else np.nan,
        "mwu_vs_background_stat": float(mwu_stat) if mwu_stat == mwu_stat else np.nan,
        "mwu_vs_background_p": float(mwu_p) if mwu_p == mwu_p else np.nan,
        "genes_used": ",".join(used.index.astype(str)),
    }


def thesis_call(gsea_row: pd.Series | None, compact: dict, is_genotype: bool) -> str:
    """Direction relative to thesis: Cldn4-high is IFN-low / barrier.

    On a KO-minus-WT rank, IFN/MHC UP after loss SUPPORTS IFN-low in Cldn4-high.
    This is injury, n=1 vs 1 — never HOLDS, never FAILED.
    """
    if not is_genotype:
        return "COMPANION (not a Cldn4-high vs low split)"
    if gsea_row is None or gsea_row.empty:
        return "ADDITIVE / UNDERPOWERED (n=1 vs 1)"
    nes = float(gsea_row["nes"])
    fdr = float(gsea_row["fdr"])
    mean_hi = compact.get("mean_logFC_Cldn4high_minus_low", np.nan)
    if nes > 0 and fdr < 0.05:
        return "SUPPORTS (IFN/MHC up after Cldn4 loss; Cldn4-high is IFN-low)"
    if nes > 0:
        return "DIRECTION SUPPORTS (NES>0, FDR≥0.05; n=1 vs 1)"
    if mean_hi == mean_hi and mean_hi < 0:
        return "DIRECTION SUPPORTS on compact mean; GSEA NES≤0"
    return "ADDITIVE / INCONCLUSIVE at n=1 vs 1 (not a fail)"


def parse_series_matrix(path: Path) -> pd.DataFrame:
    # GEO characteristics lines are ragged; titles/GSM are reliable.
    titles = gsms = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]]
    meta = {
        "GSM1232580": ("WT", "none", "Cldn4-intact / Cldn4-high", "SRR988118,SRR988119"),
        "GSM1232581": ("Cldn4 KO", "none", "Cldn4-null / Cldn4-low", "SRR988120,SRR988121"),
        "GSM1232582": ("WT", "VILI 40 cmH2O 2h", "Cldn4-intact / Cldn4-high", "SRR988122,SRR988123"),
        "GSM1232583": ("Cldn4 KO", "VILI 40 cmH2O 2h; injury similar to WT", "Cldn4-null; KOlow = injury not Cldn4", "SRR988124,SRR988125"),
        "GSM1232584": ("Cldn4 KO", "VILI 40 cmH2O 2h; injury higher than WT", "Cldn4-null; KOhigh = injury not Cldn4", "SRR988126,SRR988127"),
    }
    rows = []
    for gsm, title in zip(gsms, titles):
        genotype, vili, role, runs = meta[gsm]
        rows.append(
            {
                "gsm": gsm,
                "title": title,
                "genotype": genotype,
                "vili": vili,
                "cldn4_class": role,
                "tissue": "whole lung",
                "strain": "mixed 129S6, C57BL/6, BALB/c",
                "sra_runs": runs,
                "n_sra_runs": 2,
                "processed_count_columns": 0,
            }
        )
    return pd.DataFrame(rows)


def download_all() -> None:
    for key, url in FILES.items():
        dest = DATA / Path(url).name
        fetch(url, dest)
        # also copy into DATA with expected names
        log(f"ok {dest.name} ({dest.stat().st_size} bytes)")


def main() -> None:
    download_all()
    mouse_sets = load_mouse_sets()
    headline_sets = {k: mouse_sets[k] for k in HEADLINE}

    coverage_rows = []
    gene_rows = []
    compact_rows = []
    gsea_rows = []
    one_rows = []
    cldn4_rows = []

    for cid, meta in CONTRASTS.items():
        df = load_edger(DATA / meta["file"])
        rank = rank_from_edger(df)
        log(f"{cid}: {len(rank)} genes ranked")

        gsea = gsea_prerank(rank, headline_sets, nperm=NPERM, seed=SEED)
        gsea["fdr"] = bh_fdr(gsea["nom_p"])
        gsea["contrast"] = cid
        gsea["role"] = meta["role"]
        gsea["n_high"] = meta["n_high"]
        gsea["n_low"] = meta["n_low"]
        gsea["n_genes_ranked"] = int(len(rank))
        gsea["label"] = gsea["term"].map(HEADLINE_LABEL)
        ifn_block = gsea[gsea["term"].isin(IFN_MHC)]
        gsea_rows.append(gsea)

        for set_name, genes in COMPACT.items():
            sc = compact_score(df, genes)
            sc.update(
                {
                    "contrast": cid,
                    "role": meta["role"],
                    "set": set_name,
                    "is_genotype_split": meta["is_genotype_split"],
                }
            )
            compact_rows.append(sc)

        for g in sorted(set(FOCAL + COMPACT["IFN_COMPACT"] + COMPACT["MHC_COMPACT"] + COMPACT["TJ_COMPACT"] + COMPACT["EPI_COMPACT"])):
            if g not in df.index:
                gene_rows.append(
                    {
                        "contrast": cid,
                        "symbol": g,
                        "present": False,
                    }
                )
                continue
            r = df.loc[g]
            gene_rows.append(
                {
                    "contrast": cid,
                    "symbol": g,
                    "present": True,
                    "logFC_KO_minus_WT_or_deposited": float(r["logFC"]),
                    "logFC_Cldn4high_minus_low": float(-r["logFC"]) if meta["is_genotype_split"] else np.nan,
                    "logCPM": float(r["logCPM"]),
                    "PValue": float(r["PValue"]),
                    "FDR": float(r["FDR"]),
                    "panel": (
                        "Cldn4"
                        if g == "Cldn4"
                        else "IFN"
                        if g in COMPACT["IFN_COMPACT"]
                        else "MHC"
                        if g in COMPACT["MHC_COMPACT"]
                        else "TJ"
                        if g in COMPACT["TJ_COMPACT"]
                        else "EPI"
                        if g in COMPACT["EPI_COMPACT"]
                        else "focal"
                    ),
                }
            )

        for set_name, genes in {**headline_sets, **COMPACT}.items():
            present = [g for g in genes if g in df.index]
            coverage_rows.append(
                {
                    "contrast": cid,
                    "set": set_name,
                    "n_listed": len(genes),
                    "n_present": len(present),
                    "frac_present": len(present) / len(genes) if genes else np.nan,
                }
            )

        cldn = df.loc["Cldn4"] if "Cldn4" in df.index else None
        ifn_g = gsea[gsea["term"] == "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
        ifn_a = gsea[gsea["term"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]
        mhc = gsea[gsea["term"] == "CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]
        tj = gsea[gsea["term"] == "KEGG_TIGHT_JUNCTION"]
        ifn_c = [r for r in compact_rows if r["contrast"] == cid and r["set"] == "IFN_COMPACT"][0]
        mhc_c = [r for r in compact_rows if r["contrast"] == cid and r["set"] == "MHC_COMPACT"][0]
        tj_c = [r for r in compact_rows if r["contrast"] == cid and r["set"] == "TJ_COMPACT"][0]
        gsea_ifn = ifn_g.iloc[0] if len(ifn_g) else None
        call = thesis_call(gsea_ifn, ifn_c, meta["is_genotype_split"])

        cldn4_rows.append(
            {
                "contrast": cid,
                "role": meta["role"],
                "Cldn4_present": cldn is not None,
                "Cldn4_logFC": float(cldn["logFC"]) if cldn is not None else np.nan,
                "Cldn4_logCPM": float(cldn["logCPM"]) if cldn is not None else np.nan,
                "Cldn4_FDR": float(cldn["FDR"]) if cldn is not None else np.nan,
            }
        )
        one_rows.append(
            {
                "dataset": "GSE50927",
                "contrast": cid,
                "role": meta["role"],
                "histology": "mouse whole lung VILI / epithelium (not tumour)",
                "platform": "GPL13112 Illumina HiSeq 2000",
                "n_Cldn4high": meta["n_high"],
                "n_Cldn4low": meta["n_low"],
                "honest_n": "1 vs 1 GSM (design says duplicates; no processed count matrix)",
                "Cldn4_logFC_deposited": float(cldn["logFC"]) if cldn is not None else np.nan,
                "IFNg_NES": float(ifn_g["nes"].iloc[0]) if len(ifn_g) else np.nan,
                "IFNg_FDR": float(ifn_g["fdr"].iloc[0]) if len(ifn_g) else np.nan,
                "IFNa_NES": float(ifn_a["nes"].iloc[0]) if len(ifn_a) else np.nan,
                "IFNa_FDR": float(ifn_a["fdr"].iloc[0]) if len(ifn_a) else np.nan,
                "MHC_NES": float(mhc["nes"].iloc[0]) if len(mhc) else np.nan,
                "MHC_FDR": float(mhc["fdr"].iloc[0]) if len(mhc) else np.nan,
                "TJ_NES": float(tj["nes"].iloc[0]) if len(tj) else np.nan,
                "TJ_FDR": float(tj["fdr"].iloc[0]) if len(tj) else np.nan,
                "IFN_compact_mean_logFC_high_minus_low": ifn_c["mean_logFC_Cldn4high_minus_low"]
                if meta["is_genotype_split"]
                else np.nan,
                "MHC_compact_mean_logFC_high_minus_low": mhc_c["mean_logFC_Cldn4high_minus_low"]
                if meta["is_genotype_split"]
                else np.nan,
                "TJ_compact_mean_logFC_high_minus_low": tj_c["mean_logFC_Cldn4high_minus_low"]
                if meta["is_genotype_split"]
                else np.nan,
                "verdict": call,
                "ICI": "no",
            }
        )

    gsea_df = pd.concat(gsea_rows, ignore_index=True)
    compact_df = pd.DataFrame(compact_rows)
    gene_df = pd.DataFrame(gene_rows)
    cov_df = pd.DataFrame(coverage_rows)
    one_df = pd.DataFrame(one_rows)
    cldn4_df = pd.DataFrame(cldn4_rows)

    ann = parse_series_matrix(DATA / "GSE50927_series_matrix.txt.gz")
    related = pd.DataFrame(
        [
            {
                "accession": "GSE50927",
                "relation": "this series",
                "usable": "yes",
                "note": "only GEO series linked from PMID 25106430; 4 author EdgeR tables",
            },
            {
                "accession": "E-GEOD-50927",
                "relation": "ArrayExpress mirror",
                "usable": "duplicate",
                "note": "same five GSM; not a second matrix",
            },
            {
                "accession": "Wray 2009 PMID 19376879",
                "relation": "prior Cldn4 VILI physiology (CPE peptide / siRNA)",
                "usable": "no",
                "note": "no GEO / ArrayExpress processed series",
            },
            {
                "accession": "hyperoxia arm in PMID 25106430",
                "relation": "same paper, qRT-PCR only",
                "usable": "no",
                "note": "no deposited processed matrix",
            },
        ]
    )

    inventory = pd.DataFrame(
        [
            {"field": "GEO series", "public": "yes", "n": 1, "note": "GSE50927 only (PMID 25106430)"},
            {"field": "GSM records", "public": "yes", "n": 5, "note": "one GSM per condition; do not use 5 as pairwise n"},
            {"field": "BioSamples", "public": "yes", "n": 5, "note": "SAMN02358107–111"},
            {"field": "SRA experiments", "public": "yes", "n": 5, "note": "SRX352050–054"},
            {"field": "SRA runs", "public": "yes", "n": 10, "note": "2 FASTQ runs per GSM; not a processed count matrix"},
            {"field": "series-matrix expression rows", "public": "no", "n": 0, "note": "Sample_data_row_count=0"},
            {"field": "author EdgeR tables", "public": "yes", "n": 4, "note": "naive KO vs WT; WT VILI; KOhigh; KOlow"},
            {"field": "Cldn4 in EdgeR tables", "public": "yes", "n": 4, "note": "naive logFC −6.061; VILI induction +3.959"},
            {"field": "sample-level Cldn4 Q4 vs Q1", "public": "no", "n": 0, "note": "no per-sample matrix; cannot Spearman"},
            {"field": "primary genotype n (naive)", "public": "yes", "n": "1 vs 1", "note": "design text says duplicates; honest n is 1 vs 1 GSM"},
            {"field": "KOlow / KOhigh meaning", "public": "yes", "n": "NA", "note": "BAL protein injury strata, NOT Cldn4 expression"},
            {"field": "tumour / ICI / LUAD", "public": "no", "n": 0, "note": "injury whole lung; mixed 129S6/C57BL/6/BALB/c"},
            {"field": "related processed mouse series", "public": "no", "n": 0, "note": "Wray 2009 and hyperoxia arm have no matrix"},
        ]
    )

    gsea_df.to_csv(TAB / "gsea_headline.tsv", sep="\t", index=False)
    compact_df.to_csv(TAB / "compact_panel_scores.tsv", sep="\t", index=False)
    gene_df.to_csv(TAB / "gene_level.tsv", sep="\t", index=False)
    cov_df.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)
    one_df.to_csv(TAB / "one_row.tsv", sep="\t", index=False)
    cldn4_df.to_csv(TAB / "cldn4_logfc.tsv", sep="\t", index=False)
    ann.to_csv(TAB / "sample_annotation.tsv", sep="\t", index=False)
    related.to_csv(TAB / "related_series.tsv", sep="\t", index=False)
    inventory.to_csv(TAB / "label_inventory.tsv", sep="\t", index=False)

    summary = {
        "accession": "GSE50927",
        "species": "Mus musculus",
        "tissue": "whole lung",
        "context": "Cldn4 KO ± VILI (injury, not tumour)",
        "primary_contrast": "naive_KO_vs_WT",
        "honest_n_primary": "1 vs 1 GSM",
        "engine": f"prerank GSEA weighted KS p=1, {NPERM} gene-set permutations, seed={SEED}",
        "thesis": "Cldn4-high epithelium is barrier / IFN-low (taken as given; not audited as failed)",
        "one_row": one_df.to_dict(orient="records"),
        "cldn4": cldn4_df.to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # figures
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    plot = gsea_df.copy()
    terms = HEADLINE
    x = np.arange(len(CONTRASTS))
    width = 0.15
    colors = ["#c0392b", "#e67e22", "#2980b9", "#27ae60", "#7f8c8d"]
    for i, term in enumerate(terms):
        sub = plot[plot["term"] == term].set_index("contrast").reindex(CONTRASTS)
        ax.bar(x + (i - 2) * width, sub["nes"].to_numpy(), width, label=HEADLINE_LABEL[term], color=colors[i])
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(["naive KO vs WT\n(primary)", "WT VILI vs\nnaive WT", "KO VILIhigh\nvs WT VILI", "KO VILIlow\nvs WT VILI"], fontsize=8)
    ax.set_ylabel("NES (positive = up after Cldn4 loss / in VILI)")
    ax.set_title("GSE50927 mouse lung — prerank GSEA (n=1 vs 1 per contrast)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_gsea_nes.png", dpi=140)
    fig.savefig(FIG / "fig1_gsea_nes.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    for ax, set_name, title in zip(
        axes,
        ["IFN_COMPACT", "MHC_COMPACT", "TJ_COMPACT"],
        ["IFN compact", "MHC compact", "TJ compact (no Cldn4)"],
    ):
        sub = compact_df[compact_df["set"] == set_name]
        vals = []
        labs = []
        for cid in CONTRASTS:
            row = sub[sub["contrast"] == cid].iloc[0]
            if CONTRASTS[cid]["is_genotype_split"]:
                vals.append(row["mean_logFC_Cldn4high_minus_low"])
                labs.append(cid.replace("_", "\n"))
        ax.bar(range(len(vals)), vals, color=["#2c5f8a", "#c0392b", "#8e44ad"])
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(["naive", "VILI KOhigh", "VILI KOlow"], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("mean logFC  Cldn4-high − low" if ax is axes[0] else "")
    fig.suptitle("GSE50927 compact panels — negative = lower in Cldn4-intact WT (thesis IFN-low)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_compact_high_minus_low.png", dpi=140)
    fig.savefig(FIG / "fig2_compact_high_minus_low.pdf")
    plt.close(fig)

    log("wrote tables and figures")


if __name__ == "__main__":
    main()
