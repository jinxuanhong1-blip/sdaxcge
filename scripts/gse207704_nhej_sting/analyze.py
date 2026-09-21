#!/usr/bin/env python3
"""GSE207704 breast CLDN4 KO: c-NHEJ, STING core, Hallmark IFN / APM.

Public object: GEO cufflinks group-mean FPKM (T47D and MCF7, CLDN4-/- vs WT).
Replicates are already pooled, so there is no gene-level sample FDR.

log2FC = log2((KO + 0.5) / (WT + 0.5)) on the cufflinks locus with the
highest mean FPKM for that symbol (same collapse and pseudocount as the
CLDN4-loss prerank). Positive = higher after CLDN4 KO.

Thesis checked here, not assumed: KD -> NHEJ down, STING/IFN up.
"""
from __future__ import annotations

import json
import math
import os
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse207704_nhej_sting"
DATA = OUT / "data"
TAB = OUT / "tables"
FIG = OUT / "figures"
SET_JSON = DATA / "ifn_apm_sets.json"
RAW_GZ = DATA / "GSE207704_CLDN4_RNAseq.txt.gz"
GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/"
    "GSE207704_CLDN4_RNAseq.txt.gz"
)

PSEUDO = 0.5
WEAK = 0.25
LOW_FPKM = 1.0

COLS = {
    "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
    "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
    "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
    "T47D_WT": "T47D_WT_FPKM (fpkm)",
}

# Official symbols. CGAS is MB21D1 in this cufflinks annotation.
NHEJ = ["PRKDC", "XRCC4", "LIG4", "RIF1", "TP53BP1", "XRCC5", "XRCC6"]
STING = ["CGAS", "STING1", "TBK1", "IRF3", "STAT1"]
ALIASES = {"CGAS": ["MB21D1", "C6orf150"], "STING1": ["TMEM173", "STING"]}
ENTREZ = {
    "CGAS": "115004",
    "MB21D1": "115004",
    "STING1": "340061",
    "TMEM173": "340061",
    "TBK1": "29110",
    "IRF3": "3661",
    "STAT1": "6772",
    "PRKDC": "5591",
    "XRCC4": "7518",
    "LIG4": "3981",
    "RIF1": "55183",
    "TP53BP1": "7158",
    "XRCC5": "7520",
    "XRCC6": "2547",
    "CLDN4": "1364",
    "TACSTD2": "4070",
}

SET_ORDER = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
]
SET_LABEL = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Hallmark IFN-α",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
}

# Locked to the earlier prerank (same rank, engine, seed, set order).
LOCKED_ES = {
    ("T47D", "HALLMARK_INTERFERON_GAMMA_RESPONSE"): -0.4537188687805774,
    ("T47D", "HALLMARK_INTERFERON_ALPHA_RESPONSE"): -0.5870919845851025,
    ("T47D", "CUSTOM_MHC_I_ANTIGEN_PRESENTATION"): -0.5923068171480215,
    ("MCF7", "HALLMARK_INTERFERON_GAMMA_RESPONSE"): -0.297865653250094,
    ("MCF7", "HALLMARK_INTERFERON_ALPHA_RESPONSE"): -0.3042719321685066,
    ("MCF7", "CUSTOM_MHC_I_ANTIGEN_PRESENTATION"): 0.4128969206378983,
}


def log2fc(ko: float, wt: float) -> float:
    return float(math.log2((ko + PSEUDO) / (wt + PSEUDO)))


def direction(value: float, max_fpkm: float) -> str:
    if max_fpkm < LOW_FPKM:
        return "LOW"
    if value > WEAK:
        return "UP"
    if value < -WEAK:
        return "DOWN"
    return "FLAT"


def consensus(d_mcf7: str, d_t47d: str) -> str:
    if "ABSENT" in (d_mcf7, d_t47d):
        return "ABSENT"
    dirs = {d_mcf7, d_t47d}
    if dirs <= {"LOW"}:
        return "LOW"
    if "UP" in dirs and "DOWN" in dirs:
        return "DISCORDANT"
    if "UP" in dirs:
        return "UP"
    if "DOWN" in dirs:
        return "DOWN"
    return "FLAT"


def vs_thesis(group: str, call: str, d_mcf7: str = "", d_t47d: str = "") -> str:
    """Thesis: NHEJ genes down after KO; STING / IFN / APM genes up after KO."""
    if call == "ABSENT":
        return "unmeasured"
    expects_down = group == "cNHEJ"
    both = {d_mcf7, d_t47d}
    if expects_down:
        if call == "DOWN" and both == {"DOWN"}:
            return "agrees (down in both)"
        if call == "DOWN":
            return "down in one line"
        if call == "UP" and both == {"UP"}:
            return "opposite (up in both)"
        if call == "UP":
            return "not decreased (up in one line)"
        if call == "DISCORDANT":
            return "mixed"
        return "not decreased"
    if call == "UP" and both == {"UP"}:
        return "agrees (up in both)"
    if call == "UP":
        return "up in one line"
    if call == "DOWN":
        return "opposite (down)"
    if call == "DISCORDANT":
        return "mixed"
    return "not increased"


def norm_entrez(value: object) -> str:
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    if text in {"", "nan", "None"}:
        return ""
    return text


def ensure_matrix() -> None:
    if RAW_GZ.exists() and RAW_GZ.stat().st_size > 0:
        return
    DATA.mkdir(parents=True, exist_ok=True)
    print(f"downloading {GEO_URL}")
    urllib.request.urlretrieve(GEO_URL, RAW_GZ)


def load_fpkm() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(RAW_GZ, sep="\t", low_memory=False)
    raw = raw.rename(columns={"gene_short_name": "symbol"})
    entrez_col = [c for c in raw.columns if str(c).startswith("gene_id")][-1]
    raw["entrez"] = raw[entrez_col].map(norm_entrez)
    raw = raw[raw["symbol"].notna() & (raw["symbol"].astype(str).str.len() > 0)].copy()
    raw["symbol"] = raw["symbol"].astype(str).str.strip()
    for col in COLS.values():
        raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
    raw["_mean"] = raw[list(COLS.values())].mean(axis=1)
    primary = (
        raw.sort_values("_mean", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")
    )
    primary = primary.rename(columns={v: k for k, v in COLS.items()})
    return raw, primary


def resolve(symbol: str, primary: pd.DataFrame, entrez_to_symbol: dict[str, str]) -> tuple[str | None, str]:
    if symbol in primary.index:
        return symbol, "symbol"
    for alias in ALIASES.get(symbol, []):
        if alias in primary.index:
            return alias, f"alias:{alias}"
    entrez = ENTREZ.get(symbol, "")
    mapped = entrez_to_symbol.get(entrez, "")
    if mapped and mapped in primary.index:
        return mapped, f"entrez:{entrez}:{mapped}"
    return None, "absent"


def gene_rows(raw: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    entrez_to_symbol = {}
    for symbol, entrez in (
        raw.loc[raw["entrez"] != "", ["symbol", "entrez"]]
        .drop_duplicates()
        .itertuples(index=False)
    ):
        entrez_to_symbol.setdefault(entrez, symbol)

    value_cols = list(COLS)
    rows = []
    panel = [("cNHEJ", g) for g in NHEJ] + [("STING", g) for g in STING]
    panel += [("QC", "CLDN4"), ("QC", "TACSTD2")]
    for group, symbol in panel:
        mapped, how = resolve(symbol, primary, entrez_to_symbol)
        loci = raw[raw["symbol"] == (mapped or symbol)]
        if mapped is None:
            # alias may be the cufflinks name
            for alias in ALIASES.get(symbol, []):
                loci = raw[raw["symbol"] == alias]
                if len(loci):
                    break
        rec = {
            "group": group,
            "symbol": symbol,
            "mapped_symbol": mapped or "",
            "entrez": ENTREZ.get(symbol, ""),
            "mapping": how,
            "n_loci": int((raw["symbol"] == mapped).sum()) if mapped else 0,
            "consensus": "ABSENT",
            "vs_thesis": vs_thesis(group, "ABSENT"),
            "note": "symbol and Entrez absent from GSE207704_CLDN4_RNAseq.txt",
        }
        if mapped is None:
            rows.append(rec)
            continue
        hit = primary.loc[mapped]
        m_wt, m_ko = float(hit["MCF7_WT"]), float(hit["MCF7_KO"])
        t_wt, t_ko = float(hit["T47D_WT"]), float(hit["T47D_KO"])
        lfc_m, lfc_t = log2fc(m_ko, m_wt), log2fc(t_ko, t_wt)
        d_m = direction(lfc_m, max(m_wt, m_ko))
        d_t = direction(lfc_t, max(t_wt, t_ko))
        call = consensus(d_m, d_t)
        same = raw[raw["symbol"] == mapped]
        sum_fpkm = same[list(COLS.values())].sum()
        sum_m = log2fc(float(sum_fpkm[COLS["MCF7_KO"]]), float(sum_fpkm[COLS["MCF7_WT"]]))
        sum_t = log2fc(float(sum_fpkm[COLS["T47D_KO"]]), float(sum_fpkm[COLS["T47D_WT"]]))
        rec.update(
            {
                "MCF7_WT": m_wt,
                "MCF7_KO": m_ko,
                "T47D_WT": t_wt,
                "T47D_KO": t_ko,
                "log2FC_MCF7": lfc_m,
                "log2FC_T47D": lfc_t,
                "log2FC_mean": 0.5 * (lfc_m + lfc_t),
                "dir_MCF7": d_m,
                "dir_T47D": d_t,
                "consensus": call,
                "vs_thesis": vs_thesis(group, call, d_m, d_t),
                "sum_loci_log2FC_MCF7": sum_m,
                "sum_loci_log2FC_T47D": sum_t,
                "sum_loci_log2FC_mean": 0.5 * (sum_m + sum_t),
                "note": "descriptive; group-mean FPKM, no gene-level sample FDR",
            }
        )
        rows.append(rec)
    return pd.DataFrame(rows)


def rank_table(primary: pd.DataFrame) -> pd.DataFrame:
    lg = np.log2(primary[list(COLS)] + PSEUDO)
    out = pd.DataFrame(
        {
            "MCF7": lg["MCF7_KO"] - lg["MCF7_WT"],
            "T47D": lg["T47D_KO"] - lg["T47D_WT"],
        }
    )
    out["mean"] = out.mean(axis=1)
    return out


def sign_row(values: pd.Series) -> dict:
    values = values.dropna()
    n = int(len(values))
    n_pos = int((values > 0).sum())
    n_neg = int((values < 0).sum())
    n_zero = int((values == 0).sum())
    # Sign test ignores zeros. Two-sided.
    n_sign = n_pos + n_neg
    if n_sign == 0:
        sign_p = np.nan
    else:
        sign_p = float(binomtest(n_pos, n_sign, 0.5, alternative="two-sided").pvalue)
    if n >= 1 and np.any(values != 0):
        wil_p = float(wilcoxon(values.to_numpy(), alternative="two-sided", zero_method="wilcox").pvalue)
    else:
        wil_p = np.nan
    return {
        "n": n,
        "mean_log2FC": float(values.mean()) if n else np.nan,
        "median_log2FC": float(values.median()) if n else np.nan,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_zero": n_zero,
        "sign_p": sign_p,
        "wilcoxon_p": wil_p,
    }


def load_sets() -> dict[str, list[str]]:
    payload = json.loads(SET_JSON.read_text())
    return {k: list(payload["sets"][k]) for k in SET_ORDER}


def run_gsea(ranks: pd.DataFrame, sets: dict[str, list[str]]) -> pd.DataFrame:
    frames = []
    for contrast in ["T47D", "MCF7", "mean"]:
        rank = ranks[contrast].dropna().sort_values(ascending=False)
        scored = gsea_prerank(rank, sets, nperm=NPERM, seed=SEED)
        scored["fdr"] = bh_fdr(scored["nom_p"])
        scored["contrast"] = contrast
        scored["n_genes_ranked"] = int(rank.shape[0])
        scored["fdr_family"] = "BH within Hallmark IFN-γ, Hallmark IFN-α, and MHC-I/APM"
        frames.append(scored)
    out = pd.concat(frames, ignore_index=True)
    for _, row in out.iterrows():
        key = (row["contrast"], row["term"])
        if key in LOCKED_ES:
            if abs(float(row["es"]) - LOCKED_ES[key]) > 1e-8:
                raise SystemExit(f"ES drift for {key}: {row['es']} vs {LOCKED_ES[key]}")
    return out


def set_gene_table(ranks: pd.DataFrame, primary: pd.DataFrame, sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for term, genes in sets.items():
        for gene in genes:
            if gene not in ranks.index:
                rows.append(
                    {
                        "term": term,
                        "symbol": gene,
                        "in_matrix": 0,
                        "log2FC_MCF7": np.nan,
                        "log2FC_T47D": np.nan,
                        "log2FC_mean": np.nan,
                    }
                )
                continue
            mx = float(primary.loc[gene, list(COLS)].max())
            rows.append(
                {
                    "term": term,
                    "symbol": gene,
                    "in_matrix": 1,
                    "max_fpkm": mx,
                    "log2FC_MCF7": float(ranks.loc[gene, "MCF7"]),
                    "log2FC_T47D": float(ranks.loc[gene, "T47D"]),
                    "log2FC_mean": float(ranks.loc[gene, "mean"]),
                    "dir_MCF7": direction(float(ranks.loc[gene, "MCF7"]), mx),
                    "dir_T47D": direction(float(ranks.loc[gene, "T47D"]), mx),
                }
            )
    return pd.DataFrame(rows)


def fmt(value: float, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}f}"


def thesis_rows(genes: pd.DataFrame, gsea: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, label, expect in [
        ("cNHEJ", "c-NHEJ (7 genes)", "down"),
        ("STING", "STING core (5 genes)", "up"),
    ]:
        sub = genes[genes["group"] == group]
        measured = sub[sub["consensus"] != "ABSENT"]
        for contrast, col in [("T47D", "log2FC_T47D"), ("MCF7", "log2FC_MCF7"), ("mean", "log2FC_mean")]:
            stats_row = sign_row(measured[col])
            n_up = int((measured[f"dir_{contrast}"] == "UP").sum()) if contrast != "mean" else int((measured["consensus"] == "UP").sum())
            n_down = int((measured[f"dir_{contrast}"] == "DOWN").sum()) if contrast != "mean" else int((measured["consensus"] == "DOWN").sum())
            if contrast == "mean":
                detail = (
                    f"mean log2FC {fmt(stats_row['mean_log2FC'])}; "
                    f"{stats_row['n_pos']}/{stats_row['n']} genes > 0; "
                    f"consensus UP {int((measured['consensus']=='UP').sum())}, "
                    f"DOWN {int((measured['consensus']=='DOWN').sum())}, "
                    f"FLAT {int((measured['consensus']=='FLAT').sum())}, "
                    f"ABSENT {int((sub['consensus']=='ABSENT').sum())}"
                )
            else:
                detail = (
                    f"mean log2FC {fmt(stats_row['mean_log2FC'])}; "
                    f"UP {n_up}, DOWN {n_down}, "
                    f"of {int(stats_row['n'])} measured"
                )
            if expect == "down":
                call = "not decreased" if stats_row["n_neg"] < stats_row["n"] and stats_row["mean_log2FC"] > 0 else "see genes"
                if int((measured["consensus"] == "DOWN").sum()) == 0 and stats_row["mean_log2FC"] > 0:
                    call = "not decreased"
            else:
                n_absent = int((sub["consensus"] == "ABSENT").sum())
                if int((measured["consensus"] == "UP").sum()) == 0:
                    call = "not increased" + ("; STING1 absent" if n_absent else "")
                else:
                    call = "see genes"
            rows.append(
                {
                    "axis": label,
                    "thesis_after_KO": expect,
                    "contrast": contrast,
                    "call": call if contrast == "mean" else call,
                    "detail": detail,
                    "sign_p": stats_row["sign_p"],
                    "wilcoxon_p": stats_row["wilcoxon_p"],
                }
            )
    for _, row in gsea.iterrows():
        nes = float(row["nes"])
        fdr = float(row["fdr"])
        if nes > 0 and fdr < 0.05:
            call = "increased"
        elif nes < 0 and fdr < 0.05:
            call = "decreased (opposite of up)"
        elif nes > 0:
            call = "not increased (NES>0, FDR>=0.05)"
        else:
            call = "not increased (NES<0, FDR>=0.05)"
        rows.append(
            {
                "axis": SET_LABEL[row["term"]],
                "thesis_after_KO": "up",
                "contrast": row["contrast"],
                "call": call,
                "detail": (
                    f"NES {fmt(nes)} nominal p {fmt(float(row['nom_p']), 4)} "
                    f"FDR {fmt(fdr, 4)}; {int(row['n_set_in_rank'])} genes in rank; "
                    f"mean log2FC {fmt(float(row['mean_stat']))}"
                ),
                "sign_p": np.nan,
                "wilcoxon_p": np.nan,
            }
        )
    return pd.DataFrame(rows)


def write_finding(genes: pd.DataFrame, gsea: pd.DataFrame, sets: dict[str, list[str]]) -> None:
    def g(symbol: str) -> pd.Series:
        return genes.loc[genes["symbol"] == symbol].iloc[0]

    cldn4 = g("CLDN4")
    nhej = genes[genes["group"] == "cNHEJ"]
    sting = genes[genes["group"] == "STING"]

    def gene_md(frame: pd.DataFrame) -> str:
        lines = [
            "| Gene | Mapped | MCF7 log2FC | T47D log2FC | Mean | MCF7 | T47D | Both lines | vs thesis |",
            "|---|---|---:|---:|---:|---|---|---|---|",
        ]
        for _, r in frame.iterrows():
            if r["consensus"] == "ABSENT":
                lines.append(
                    f"| {r['symbol']} | absent | NA | NA | NA | ABSENT | ABSENT | ABSENT | {r['vs_thesis']} |"
                )
                continue
            mapped = r["symbol"] if r["mapped_symbol"] == r["symbol"] else f"{r['mapped_symbol']}"
            lines.append(
                "| {sym} | {mapped} | {m} | {t} | {mn} | {dm} | {dt} | {c} | {v} |".format(
                    sym=r["symbol"],
                    mapped=mapped,
                    m=fmt(r["log2FC_MCF7"]),
                    t=fmt(r["log2FC_T47D"]),
                    mn=fmt(r["log2FC_mean"]),
                    dm=r["dir_MCF7"],
                    dt=r["dir_T47D"],
                    c=r["consensus"],
                    v=r["vs_thesis"],
                )
            )
        return "\n".join(lines)

    def gsea_md() -> str:
        lines = [
            "| Set | Contrast | NES | nominal p | BH-FDR | Genes in rank | Mean log2FC | Call vs IFN/APM up |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
        order = {"T47D": 0, "MCF7": 1, "mean": 2}
        view = gsea.copy()
        view["_o"] = view["contrast"].map(order)
        view["_t"] = view["term"].map({k: i for i, k in enumerate(SET_ORDER)})
        view = view.sort_values(["_t", "_o"])
        for _, r in view.iterrows():
            nes, fdr = float(r["nes"]), float(r["fdr"])
            if nes < 0 and fdr < 0.05:
                call = "decreased"
            elif nes > 0 and fdr < 0.05:
                call = "increased"
            elif nes < 0:
                call = "not increased"
            else:
                call = "not increased"
            contrast = r["contrast"] if r["contrast"] != "mean" else "mean of the two lines"
            lines.append(
                f"| {SET_LABEL[r['term']]} | {contrast} | {fmt(nes)} | {fmt(float(r['nom_p']), 4)} | "
                f"{fmt(fdr, 4)} | {int(r['n_set_in_rank'])} | {fmt(float(r['mean_stat']))} | {call} |"
            )
        return "\n".join(lines)

    nhej_mean = sign_row(nhej["log2FC_mean"])
    sting_measured = sting[sting["consensus"] != "ABSENT"]
    sting_mean = sign_row(sting_measured["log2FC_mean"])
    prkdc = g("PRKDC")
    lig4 = g("LIG4")

    ifna = sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
    ifng = sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]
    apm = sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]

    # lead genes from T47D negative IFN sets
    def leads(term: str, contrast: str) -> str:
        hit = gsea[(gsea["term"] == term) & (gsea["contrast"] == contrast)].iloc[0]
        return str(hit["lead_genes"])

    text = f"""# GSE207704 breast CLDN4 KO — c-NHEJ, STING, Hallmark IFN / APM

**FINAL discordant cancer-line KD.** The cufflinks directions below were re-checked from the open SRA runs (kallisto on Ensembl 90, both replicates). No c-NHEJ set was down. T47D Hallmark IFN stayed down (IFN-α NES −1.745, FDR 0.0077). MCF7 IFN NES became weakly positive and was not FDR < 0.05, and the c-NHEJ panel there was up. STING1/TMEM173 is quantified and essentially unexpressed. The sweep, the sets, and the count tables are in `SWEEP.md`.

Public cufflinks table (T47D and MCF7, CLDN4 CRISPR KO vs parental WT). The question is the direction against **KD → NHEJ down, STING/IFN up**.

c-NHEJ transcripts stay flat to slightly higher (0/7 consensus down). The four measured STING-core genes stay flat, and STING1 is absent from the deposit. Hallmark IFN-α and IFN-γ are lower after knockout in T47D (both BH-FDR < 0.05) and are not higher in MCF7. MHC-I/APM is mostly missing from the file; the eight genes that remain are not higher together.

## Design

Murakami et al., Breast Cancer Research 2023 (GEO [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704)). Two breast lines, CLDN4−/− vs WT. GEO lists 2 biological replicates per genotype (GSM6310640–GSM6310647). The only GEO expression file collapses each genotype to one FPKM. Ranks below are those group means, so this file has no gene-level sample FDR. The open SRA runs were re-quantified separately; that result is in `SWEEP.md`.

log2FC = log2((KO + 0.5) / (WT + 0.5)). One symbol can have more than one cufflinks locus; the reported value is the locus with the higher mean FPKM. A second column in `gene_panel.tsv` sums loci. UP / DOWN requires |log2FC| > 0.25 and max FPKM ≥ 1; otherwise FLAT. Consensus UP requires one line UP and the other not DOWN (same rule in the other direction).

CLDN4 QC, higher-mean locus: MCF7 {fmt(cldn4['log2FC_MCF7'])} (FPKM {fmt(cldn4['MCF7_WT'], 1)} → {fmt(cldn4['MCF7_KO'], 1)}), T47D {fmt(cldn4['log2FC_T47D'])} (FPKM {fmt(cldn4['T47D_WT'], 1)} → {fmt(cldn4['T47D_KO'], 1)}). Both lines are DOWN. Residual mRNA is about {100 * cldn4['MCF7_KO'] / cldn4['MCF7_WT']:.0f}% (MCF7) and {100 * cldn4['T47D_KO'] / cldn4['T47D_WT']:.0f}% (T47D).

## Scorecard

| Axis | After CLDN4 KO, thesis says | What this file shows |
|---|---|---|
| c-NHEJ (PRKDC, XRCC4, LIG4, RIF1, TP53BP1, XRCC5, XRCC6) | down | **Not decreased.** 7/7 genes present. Both-line mean log2FC {fmt(nhej_mean['mean_log2FC'])}; {nhej_mean['n_pos']}/7 positive. Consensus DOWN: 0. Consensus UP: LIG4, PRKDC (T47D only). |
| STING core (CGAS, STING1, TBK1, IRF3, STAT1) | up | **Not increased.** STING1 / TMEM173 absent. The other four are FLAT (\\|log2FC\\| ≤ 0.25). Their mean log2FC is {fmt(sting_mean['mean_log2FC'])}. |
| Hallmark IFN-α / IFN-γ | up | **Not increased.** T47D NES is negative (IFN-α {fmt(float(gsea[(gsea.term.str.contains('ALPHA')) & (gsea.contrast=='T47D')].nes.iloc[0]))}, IFN-γ {fmt(float(gsea[(gsea.term.str.contains('GAMMA')) & (gsea.contrast=='T47D')].nes.iloc[0]))}; both BH-FDR < 0.05). MCF7 NES is negative and not FDR < 0.05. |
| MHC-I / APM (21-gene list; Hallmark has no APM set) | up | **Not increased.** 8/21 genes are in the file. NES T47D {fmt(float(gsea[(gsea.term.str.contains('MHC')) & (gsea.contrast=='T47D')].nes.iloc[0]))}, MCF7 {fmt(float(gsea[(gsea.term.str.contains('MHC')) & (gsea.contrast=='MCF7')].nes.iloc[0]))}, neither FDR < 0.05. |

## 1. c-NHEJ

{gene_md(nhej)}

Call rule on the both-line mean: {nhej_mean['n_pos']} positive, {nhej_mean['n_neg']} negative. Two-sided sign test p = {fmt(nhej_mean['sign_p'], 3)}. Wilcoxon signed-rank p = {fmt(nhej_mean['wilcoxon_p'], 3)}; that smaller p is the single negative gene (XRCC6, mean log2FC {fmt(g('XRCC6')['log2FC_mean'])}) having the smallest magnitude. These p-values are gene concordance on one collapsed rank, not replicate tests. They do not establish an NHEJ increase. They do show the panel is not the predicted decrease.

The only |log2FC| > 0.25 calls are T47D-limited and point up: LIG4 {fmt(lig4['log2FC_T47D'])} (MCF7 {fmt(lig4['log2FC_MCF7'])}) and PRKDC {fmt(prkdc['log2FC_T47D'])} (MCF7 {fmt(prkdc['log2FC_MCF7'])}, just short of the −0.25 line). XRCC4, RIF1, TP53BP1, XRCC5, and XRCC6 stay FLAT in both lines.

PRKDC and TBK1 each have a second, lower-FPKM cufflinks fragment of the same locus. Summing fragments does not create a both-line NHEJ decrease: PRKDC sum-of-loci log2FC is MCF7 {fmt(prkdc['sum_loci_log2FC_MCF7'])}, T47D {fmt(prkdc['sum_loci_log2FC_T47D'])}. Primary numbers stay on the higher-mean locus so they match the earlier rank.

c-NHEJ has 7 genes, under the prerank minimum of 8, so there is no NES for this panel.

## 2. STING pathway

{gene_md(sting)}

CGAS is in the file as **MB21D1** (Entrez 115004): MCF7 {fmt(g('CGAS')['log2FC_MCF7'])}, T47D {fmt(g('CGAS')['log2FC_T47D'])}. STING1 and TMEM173 (Entrez 340061) are not in the table, so the receptor is unmeasured rather than a zero. TBK1, IRF3, and STAT1 are expressed (FPKM tens) and FLAT. Of the four measured genes, {sting_mean['n_pos']}/4 have a positive both-line mean, sign-test p = {fmt(sting_mean['sign_p'], 3)}, and every |log2FC| is ≤ 0.25. That is not STING up.

## 3. Hallmark IFN and APM

Prerank GSEA, weighted KS p = 1, {NPERM} gene-set permutations, seed {SEED}, same engine as the CLDN4-loss prerank. Positive NES = the set sits at the CLDN4-KO end of the rank. BH-FDR is within these three sets for that contrast. Enrichment scores match the earlier GSE207704 rows (same rank and the same set order). The earlier FDR was BH across five headline sets, so the FDR column here is not the same number even though NES is.

{gsea_md()}

Coverage in this cufflinks annotation: Hallmark IFN-α {int(gsea[(gsea.term=='HALLMARK_INTERFERON_ALPHA_RESPONSE') & (gsea.contrast=='T47D')].n_set_in_rank.iloc[0])}/{len(ifna)}, Hallmark IFN-γ {int(gsea[(gsea.term=='HALLMARK_INTERFERON_GAMMA_RESPONSE') & (gsea.contrast=='T47D')].n_set_in_rank.iloc[0])}/{len(ifng)}, MHC-I/APM {int(gsea[(gsea.term=='CUSTOM_MHC_I_ANTIGEN_PRESENTATION') & (gsea.contrast=='T47D')].n_set_in_rank.iloc[0])}/{len(apm)}. Genes missing from the deposit cannot move the rank. The T47D IFN decrease is carried by genes that are present. MCF7 IFN-α mean log2FC is {fmt(float(gsea[(gsea.term=='HALLMARK_INTERFERON_ALPHA_RESPONSE') & (gsea.contrast=='MCF7')].mean_stat.iloc[0]))} while its NES is negative: the set is not shifted up. Leading edge, T47D IFN-α (KO-low end): {leads('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'T47D')}. T47D IFN-γ: {leads('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'T47D')}.

APM genes absent here include HLA-A, HLA-B, HLA-E, HLA-F, HLA-G, B2M, TAP1, TAP2, TAPBP, NLRC5, PSMB8, PSMB9, and ERAP2. Of the eight that are present, ERAP1 is down in both lines (MCF7 −0.729, T47D −0.386). IRF1 is down in T47D (−0.477) and flat in MCF7 (−0.154). HLA-C is up in MCF7 (+1.287) and down in T47D (−0.499). PSMB10 is up in MCF7 only (+0.357). No APM gene is up in both lines.

The “mean of the two lines” rank averages the two group-mean log2FC vectors. It is a summary, not a third cohort.

## Reading

Relative to KD → NHEJ down and STING/IFN up, GSE207704’s breast CLDN4 KO is not a supporting public example. NHEJ is flat to slightly higher (LIG4 and PRKDC up in T47D only). STING core is flat where it is measured. Hallmark interferon moves down in T47D. APM cannot be scored for classical MHC-I genes because they are absent, and the genes that remain are not up.

Numbers in this file are descriptive directions on pooled FPKM. The replicate-level kallisto and PyDESeq2 result is in `SWEEP.md`.

## Reproduce

```bash
python3 scripts/gse207704_nhej_sting/analyze.py
```

Inputs: `methods/gse207704_nhej_sting/data/GSE207704_CLDN4_RNAseq.txt.gz` (downloaded from GEO if missing) and `ifn_apm_sets.json`. Figure: `methods/gse207704_nhej_sting/figures/fig_nhej_sting_ifn.png`.
"""
    # The f-string above uses gsea.term.str which needs the dataframe name in scope.
    # It is evaluated here, so `gsea` must support .term — pandas does via attribute access.
    (OUT / "FINDING.md").write_text(text)


def plot(genes: pd.DataFrame, gsea: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    panel = genes[genes["group"].isin(["cNHEJ", "STING"])].copy()
    # Plot top-to-bottom in the requested order.
    order = list(reversed(NHEJ + STING))
    panel = panel.set_index("symbol").loc[order]
    y = np.arange(len(order))
    fig, axes = plt.subplots(
        1, 2, figsize=(10.2, 5.8), gridspec_kw={"width_ratios": [1.35, 1.0]}
    )
    ax = axes[0]
    h = 0.36
    present = panel["consensus"].to_numpy() != "ABSENT"
    mcf = panel["log2FC_MCF7"].to_numpy(dtype=float)
    t47 = panel["log2FC_T47D"].to_numpy(dtype=float)
    # STING block is the first 5 rows of `order` (bottom of the axis).
    ax.axhspan(-0.55, 4.55, color="#F4F1EA", zorder=0)
    ax.axhspan(4.55, len(order) - 0.45, color="#F7F7F7", zorder=0)
    ax.barh(y[present] + h / 2, t47[present], height=h, color="#3C5488", label="T47D", zorder=2)
    ax.barh(y[present] - h / 2, mcf[present], height=h, color="#E64B35", label="MCF7", zorder=2)
    for i, sym in enumerate(order):
        if panel.loc[sym, "consensus"] == "ABSENT":
            ax.plot([0], [y[i]], marker="x", color="#666666", markersize=7, zorder=3)
            ax.text(0.04, y[i], "not in deposit", va="center", ha="left", fontsize=7, color="#444444")
    ax.axvline(0, color="#222222", lw=0.6)
    ax.axvline(WEAK, color="#888888", lw=0.6, ls="--")
    ax.axvline(-WEAK, color="#888888", lw=0.6, ls="--")
    ax.set_yticks(y)
    labels = []
    for sym in order:
        if sym == "CGAS":
            labels.append("CGAS (MB21D1)")
        elif sym == "STING1":
            labels.append("STING1 (absent)")
        else:
            labels.append(sym)
    ax.set_yticklabels(labels)
    for tick, sym in zip(ax.get_yticklabels(), order):
        tick.set_color("#3C5488" if sym in NHEJ else "#8C5A2A")
    ax.set_xlabel("log2FC  CLDN4 KO vs WT")
    ax.set_title("c-NHEJ and STING core")
    ax.legend(frameon=False, loc="lower left")
    ax.set_xlim(-1.15, 1.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    terms = SET_ORDER
    contrasts = ["T47D", "MCF7"]
    colors = {"T47D": "#3C5488", "MCF7": "#E64B35"}
    ypos = np.arange(len(terms))
    for j, contrast in enumerate(contrasts):
        vals = []
        fdrs = []
        for term in terms:
            hit = gsea[(gsea["term"] == term) & (gsea["contrast"] == contrast)].iloc[0]
            vals.append(float(hit["nes"]))
            fdrs.append(float(hit["fdr"]))
        offset = 0.16 if j == 0 else -0.16
        ax.barh(
            ypos + offset,
            vals,
            height=0.3,
            color=colors[contrast],
            label=contrast,
            zorder=2,
        )
        for i, (nes, fdr) in enumerate(zip(vals, fdrs)):
            if fdr < 0.05:
                x = nes - 0.06 if nes < 0 else nes + 0.06
                ax.text(
                    x,
                    ypos[i] + offset,
                    "*",
                    va="center",
                    ha="right" if nes < 0 else "left",
                    fontsize=13,
                    color=colors[contrast],
                )
    ax.axvline(0, color="#222222", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels([SET_LABEL[t] for t in terms])
    ax.set_title("Hallmark IFN and APM")
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-2.25, 1.45)
    ax.set_xlabel("NES  (positive = higher after KO)\n* BH-FDR < 0.05 within these 3 sets")

    fig.suptitle("GSE207704  CLDN4 KO vs WT   group-mean FPKM", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "fig_nhej_sting_ifn.png", dpi=160)
    fig.savefig(FIG / "fig_nhej_sting_ifn.pdf")
    plt.close(fig)


def main() -> None:
    ensure_matrix()
    raw, primary = load_fpkm()
    if "CLDN4" not in primary.index:
        raise SystemExit("CLDN4 missing from matrix")
    ranks = rank_table(primary)
    # Lock CLDN4 log2FC to the earlier prerank.
    if abs(float(ranks.loc["CLDN4", "MCF7"]) - (-0.749921130859537)) > 1e-9:
        raise SystemExit("MCF7 CLDN4 log2FC drifted")
    if abs(float(ranks.loc["CLDN4", "T47D"]) - (-1.0385143134411772)) > 1e-9:
        raise SystemExit("T47D CLDN4 log2FC drifted")
    genes = gene_rows(raw, primary)
    sets = load_sets()
    gsea = run_gsea(ranks, sets)
    members = set_gene_table(ranks, primary, sets)
    TAB.mkdir(parents=True, exist_ok=True)
    genes.to_csv(TAB / "gene_panel.tsv", sep="\t", index=False, float_format="%.6g")
    gsea.to_csv(TAB / "gsea_ifn_apm.tsv", sep="\t", index=False, float_format="%.6g")
    members.to_csv(TAB / "ifn_apm_genes.tsv", sep="\t", index=False, float_format="%.6g")
    thesis = thesis_rows(genes, gsea)
    thesis.to_csv(TAB / "thesis_direction.tsv", sep="\t", index=False, float_format="%.6g")
    write_finding(genes, gsea, sets)
    plot(genes, gsea)
    print(f"wrote {OUT}")
    print(genes.loc[genes["group"] != "QC", ["symbol", "log2FC_MCF7", "log2FC_T47D", "consensus", "vs_thesis"]].to_string(index=False))


if __name__ == "__main__":
    main()
