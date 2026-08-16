#!/usr/bin/env python3
"""Score CLDN4/TJ and IFN/MHC-I/APM in public TROP2-ADC / TROP2-loss analogs.

Public only. IMMU132 = sacituzumab govitecan (SG), not SKB264/sac-TMT.
TROP2 KO series are TROP2-loss analogs, not ADC treatment.
"""

from __future__ import annotations

import json
import sys
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "c_public_trop2_adc_analogs"
RAW = OUT / "raw"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import APM, C4_IFN_MHCI, IFN, JUNCTION, MOUSE_ENSEMBL, TARGETS, all_sets

PSEUDO = 1.0
DPI = 160


def welch_rows(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    n_a = np.isfinite(a).sum(axis=1)
    n_b = np.isfinite(b).sum(axis=1)
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    var_a = np.nanvar(a, axis=1, ddof=1)
    var_b = np.nanvar(b, axis=1, ddof=1)
    log2fc = mean_a - mean_b
    se2 = var_a / np.maximum(n_a, 1) + var_b / np.maximum(n_b, 1)
    tvals = log2fc / np.sqrt(se2)
    denom = (var_a / np.maximum(n_a, 1)) ** 2 / np.maximum(n_a - 1, 1) + (
        var_b / np.maximum(n_b, 1)
    ) ** 2 / np.maximum(n_b - 1, 1)
    dfree = se2**2 / denom
    pvals = 2 * stats.t.sf(np.abs(tvals), dfree)
    bad = (n_a < 2) | (n_b < 2) | ~np.isfinite(se2) | (se2 == 0)
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "n_treat": n_a,
            "n_ctrl": n_b,
            "mean_treat_log2": mean_a,
            "mean_ctrl_log2": mean_b,
            "log2FC": log2fc,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out


def paired_rows(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    d = a - b
    n = np.isfinite(d).sum(axis=1)
    mean_d = np.nanmean(d, axis=1)
    sd = np.nanstd(d, axis=1, ddof=1)
    tvals = mean_d / (sd / np.sqrt(np.maximum(n, 1)))
    pvals = 2 * stats.t.sf(np.abs(tvals), n - 1)
    bad = (n < 3) | ~np.isfinite(sd) | (sd == 0)
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "n_pairs": n,
            "mean_treat_log2": np.nanmean(a, axis=1),
            "mean_ctrl_log2": np.nanmean(b, axis=1),
            "log2FC": mean_d,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out


def geneset_mw(de: pd.DataFrame, genes: list[str], expressed: set[str] | None = None) -> dict:
    present = [g for g in genes if g in set(de["gene"])]
    set_fc = de.loc[de["gene"].isin(present), "log2FC"].dropna()
    if expressed is None:
        bg = de.loc[~de["gene"].isin(present), "log2FC"].dropna()
    else:
        bg = de.loc[de["gene"].isin(expressed) & ~de["gene"].isin(present), "log2FC"].dropna()
    rec = {
        "n_in_set": len(genes),
        "n_present": len(present),
        "n_tested": int(len(set_fc)),
        "median_log2FC_set": float(set_fc.median()) if len(set_fc) else None,
        "median_log2FC_bg": float(bg.median()) if len(bg) else None,
        "n_up": int((set_fc > 0).sum()) if len(set_fc) else 0,
        "n_down": int((set_fc < 0).sum()) if len(set_fc) else 0,
    }
    if len(set_fc) >= 1:
        rec["binom_p"] = float(stats.binomtest(int((set_fc > 0).sum()), n=int(len(set_fc)), p=0.5).pvalue)
    else:
        rec["binom_p"] = None
    if len(set_fc) < 3 or len(bg) < 20:
        rec.update({"mw_p": None, "direction": "untestable"})
        return rec
    u, p = stats.mannwhitneyu(set_fc, bg, alternative="two-sided")
    med_s, med_b = float(set_fc.median()), float(bg.median())
    if p >= 0.05:
        direction = "null"
    elif med_s > med_b:
        direction = "up"
    else:
        direction = "down"
    rec.update(
        {
            "delta_median": med_s - med_b,
            "mw_U": float(u),
            "mw_p": float(p),
            "direction": direction,
        }
    )
    return rec


def zmean_score(logx: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in logx.index]
    sub = logx.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
    return z.mean(axis=0)


def score_unpaired(scores: pd.Series, treat: list[str], ctrl: list[str]) -> dict:
    a = scores[treat].to_numpy(dtype=float)
    b = scores[ctrl].to_numpy(dtype=float)
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return {
        "n_treat": int(len(a)),
        "n_ctrl": int(len(b)),
        "mean_treat": float(np.mean(a)),
        "mean_ctrl": float(np.mean(b)),
        "delta": float(np.mean(a) - np.mean(b)),
        "t": float(t),
        "p": float(p),
    }


def score_paired(scores: pd.Series, treat: list[str], ctrl: list[str]) -> dict:
    a = scores[treat].to_numpy(dtype=float)
    b = scores[ctrl].to_numpy(dtype=float)
    d = a - b
    t, p = stats.ttest_rel(a, b)
    try:
        wp = float(stats.wilcoxon(d).pvalue)
    except ValueError:
        wp = None
    return {
        "n_pairs": int(len(d)),
        "mean_treat": float(np.mean(a)),
        "mean_ctrl": float(np.mean(b)),
        "delta": float(np.mean(d)),
        "t": float(t),
        "p": float(p),
        "wilcoxon_p": wp,
    }


def key_gene_rows(de: pd.DataFrame, expr: pd.DataFrame | None, cols_treat: list[str], cols_ctrl: list[str]) -> pd.DataFrame:
    genes = TARGETS + C4_IFN_MHCI + ["CLDN1", "CLDN7", "TJP1", "OCLN", "F11R", "B2M", "TAP1", "TAP2"]
    seen = []
    rows = []
    for g in genes:
        if g in seen or g not in set(de["gene"]):
            continue
        seen.append(g)
        r = de.loc[de["gene"] == g].iloc[0]
        rec = {
            "gene": g,
            "log2FC": float(r["log2FC"]),
            "p": None if pd.isna(r["p"]) else float(r["p"]),
            "q": None if pd.isna(r["q"]) else float(r["q"]),
        }
        if expr is not None and g in expr.index:
            rec["mean_treat"] = float(expr.loc[g, cols_treat].mean())
            rec["mean_ctrl"] = float(expr.loc[g, cols_ctrl].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def collapse_human_fpkm(path: Path, encoding: str | None = None) -> pd.DataFrame:
    kw = {"sep": "\t"}
    if encoding:
        kw["encoding"] = encoding
    df = pd.read_csv(path, **kw)
    df["gene_name"] = df["gene_name"].astype(str).str.strip()
    sample_cols = [c for c in df.columns if c not in {
        "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
        "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family",
    }]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[df["gene_name"].notna() & (df["gene_name"] != "") & (df["gene_name"] != "nan")]
    if "gene_biotype" in df.columns:
        coding = df[df["gene_biotype"].astype(str) == "protein_coding"]
        rest = df[df["gene_biotype"].astype(str) != "protein_coding"]
        df = pd.concat([coding, rest], ignore_index=True)
    df["mean_all"] = df[sample_cols].mean(axis=1)
    df = df.sort_values("mean_all", ascending=False).drop_duplicates("gene_name", keep="first")
    return df.set_index("gene_name")[sample_cols]


def summarize_contrast(name: str, de: pd.DataFrame, logx: pd.DataFrame, treat, ctrl, paired: bool, expr=None) -> dict:
    expressed = set(logx.index[logx.mean(axis=1) > np.log2(PSEUDO + 1)].astype(str))
    sets = {}
    for sname, genes in all_sets().items():
        sets[sname] = geneset_mw(de, genes, expressed)
    scores = {}
    for sname, genes in all_sets().items():
        sc = zmean_score(logx, genes)
        scores[sname] = score_paired(sc, treat, ctrl) if paired else score_unpaired(sc, treat, ctrl)
    key = key_gene_rows(de, expr, treat, ctrl)
    return {
        "contrast": name,
        "paired": paired,
        "n_treat": len(treat),
        "n_ctrl": len(ctrl),
        "sets": sets,
        "scores": scores,
        "key_genes": key.to_dict(orient="records"),
    }


def load_gse235812() -> pd.DataFrame | None:
    tar_path = RAW / "GSE235812" / "GSE235812_RAW.tar"
    if not tar_path.exists():
        return None
    genes_keep = set(TARGETS + C4_IFN_MHCI + APM + IFN + JUNCTION)
    cols = {}
    with tarfile.open(tar_path, "r") as tar:
        for m in tar.getmembers():
            if not m.name.endswith(".tsv.gz"):
                continue
            gsm = Path(m.name).name.split("_")[0]
            f = tar.extractfile(m)
            if f is None or m.size == 0:
                continue
            try:
                df = pd.read_csv(f, sep="\t", compression="gzip", usecols=["gene_name", "RPKM"])
            except (ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
                continue
            if "gene_name" not in df.columns or "RPKM" not in df.columns:
                continue
            df = df[df["gene_name"].isin(genes_keep)].drop_duplicates("gene_name")
            cols[gsm] = df.set_index("gene_name")["RPKM"]
    mat = pd.DataFrame(cols)
    # series order: GSM7509978 = BCX-006.P0 ... alternating P0/P1
    return mat


def style_ax(ax, title: str) -> None:
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def box_groups(ax, logx: pd.DataFrame, gene: str, groups: dict[str, list[str]], colors) -> None:
    data, labs = [], []
    for lab, cols in groups.items():
        if gene not in logx.index:
            continue
        data.append(logx.loc[gene, cols].to_numpy(dtype=float))
        labs.append(f"{lab}\nn={len(cols)}")
    if not data:
        ax.text(0.5, 0.5, f"{gene} absent", ha="center")
        return
    bp = ax.boxplot(data, tick_labels=labs, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    for i, d in enumerate(data, start=1):
        ax.scatter(np.full(len(d), i) + np.random.default_rng(0).uniform(-0.08, 0.08, len(d)), d, s=18, c="k", zorder=3)
    ax.set_ylabel(f"log2({gene}+1)")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)

    catalog = []
    summaries = []
    key_all = []

    # ----- GSE312098 CX-1 CRC IMMU132 (SG analog), n=3 vs 3 -----
    fp = collapse_human_fpkm(RAW / "GSE312098_gene_fpkm.txt.gz", encoding="utf-16")
    g312_ctrl = ["X_1", "X_2", "X_3"]
    g312_immu = ["X_4", "X_5", "X_6"]
    log312 = np.log2(fp[g312_ctrl + g312_immu] + PSEUDO)
    de312 = welch_rows(log312, g312_immu, g312_ctrl)
    de312.to_csv(TABLES / "de_GSE312098_IMMU132_vs_ctrl.tsv.gz", sep="\t", index=False, compression="gzip")
    s312 = summarize_contrast("GSE312098_CX1_IMMU132_vs_ctrl", de312, log312, g312_immu, g312_ctrl, False, fp)
    s312.update({
        "accession": "GSE312098",
        "agent": "IMMU132 (sacituzumab govitecan)",
        "is_skb264": False,
        "is_adc_treatment": True,
        "tissue": "CRC cell line CX-1",
        "lung": False,
        "label": "non-lung analog",
        "time": "2 days",
        "note": "SG analog, not SKB264. FPKM, n=3 vs 3.",
    })
    summaries.append(s312)
    catalog.append({**{k: s312[k] for k in ("accession", "agent", "is_adc_treatment", "tissue", "lung", "label", "note")}, "usable_adc_vs_ctrl": True})

    # ----- GSE311016 CRC PDX IMMU132, n=5 paired -----
    fp311 = collapse_human_fpkm(RAW / "GSE311016_gene_fpkm.txt.gz", encoding="utf-16")
    pairs = ["36", "82", "83", "114", "196"]
    g311_ctrl = [f"C_{p}" for p in pairs]
    g311_immu = [f"T_{p}" for p in pairs]
    log311 = np.log2(fp311[g311_ctrl + g311_immu] + PSEUDO)
    de311 = paired_rows(log311, g311_immu, g311_ctrl)
    de311.to_csv(TABLES / "de_GSE311016_PDX_IMMU132_vs_ctrl.tsv.gz", sep="\t", index=False, compression="gzip")
    s311 = summarize_contrast("GSE311016_CRC_PDX_IMMU132_vs_ctrl", de311, log311, g311_immu, g311_ctrl, True, fp311)
    s311.update({
        "accession": "GSE311016",
        "agent": "IMMU132 (sacituzumab govitecan)",
        "is_skb264": False,
        "is_adc_treatment": True,
        "tissue": "CRC PDX (5 models, day 29)",
        "lung": False,
        "label": "non-lung analog",
        "time": "29 days",
        "note": "SG analog, not SKB264. Paired PDX FPKM, n=5 pairs.",
    })
    summaries.append(s311)
    catalog.append({**{k: s311[k] for k in ("accession", "agent", "is_adc_treatment", "tissue", "lung", "label", "note")}, "usable_adc_vs_ctrl": True})

    # ----- GSE304294 KYSE30 ESCC IMMU132, n=2 vs 3 -----
    fp304 = collapse_human_fpkm(RAW / "GSE304294_gene_fpkm.txt.gz")
    g304_ctrl = ["OX1_1", "OX1_2", "OX1_3"]
    g304_immu = ["OX2_1", "OX2_2"]
    log304 = np.log2(fp304[g304_ctrl + g304_immu] + PSEUDO)
    de304 = welch_rows(log304, g304_immu, g304_ctrl)
    de304.to_csv(TABLES / "de_GSE304294_IMMU132_vs_ctrl.tsv.gz", sep="\t", index=False, compression="gzip")
    s304 = summarize_contrast("GSE304294_KYSE30_IMMU132_vs_ctrl", de304, log304, g304_immu, g304_ctrl, False, fp304)
    s304.update({
        "accession": "GSE304294",
        "agent": "IMMU132 (sacituzumab govitecan)",
        "is_skb264": False,
        "is_adc_treatment": True,
        "tissue": "ESCC cell line KYSE30",
        "lung": False,
        "label": "non-lung analog",
        "time": "1 day",
        "note": "SG analog, not SKB264. IMMU arm n=2 vs control n=3. Underpowered.",
    })
    summaries.append(s304)
    catalog.append({**{k: s304[k] for k in ("accession", "agent", "is_adc_treatment", "tissue", "lung", "label", "note")}, "usable_adc_vs_ctrl": True})

    # ----- GSE334497 4T1 Trop2 KO vs WT, NOT ADC -----
    m334 = pd.read_csv(RAW / "GSE334497_normalized_counts.csv.gz", index_col=0)
    wt = ["control170", "RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R"]
    ko = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
    # remap to human symbols
    rows = {}
    for hum, eid in MOUSE_ENSEMBL.items():
        if eid in m334.index:
            rows[hum] = m334.loc[eid]
    expr334 = pd.DataFrame(rows).T[wt + ko]
    log334 = np.log2(expr334 + PSEUDO)
    de334 = welch_rows(log334, ko, wt)
    # background from all genes on Ensembl IDs then attach symbols for set tests
    log334_all = np.log2(m334[wt + ko] + PSEUDO)
    de334_all = welch_rows(log334_all, ko, wt)
    # replace gene with symbol where mapped
    inv = {v: k for k, v in MOUSE_ENSEMBL.items()}
    de334_all["symbol"] = de334_all["gene"].map(inv)
    de334_sym = de334_all.dropna(subset=["symbol"]).copy()
    de334_sym["gene"] = de334_sym["symbol"]
    de334.to_csv(TABLES / "de_GSE334497_Trop2KO_vs_WT_symbols.tsv", sep="\t", index=False)
    s334 = summarize_contrast("GSE334497_4T1_Trop2KO_vs_WT", de334_sym, log334, ko, wt, False, expr334)
    # set MW against full mouse background
    for sname, genes in all_sets().items():
        s334["sets"][sname] = geneset_mw(de334_all.assign(gene=de334_all["symbol"].fillna(de334_all["gene"])), genes)
    s334.update({
        "accession": "GSE334497",
        "agent": "none (Trop2 knockout)",
        "is_skb264": False,
        "is_adc_treatment": False,
        "tissue": "4T1 mouse TNBC tumors",
        "lung": False,
        "label": "non-lung TROP2-loss analog (NOT ADC)",
        "time": "3 weeks tumor growth",
        "note": "TROP2 KO, not sacituzumab / not SKB264. Mouse 4T1, n=5 vs 5.",
    })
    summaries.append(s334)
    catalog.append({**{k: s334[k] for k in ("accession", "agent", "is_adc_treatment", "tissue", "lung", "label", "note")}, "usable_adc_vs_ctrl": False})

    # ----- GSE289287 T-47D Trop2 KO xenografts, author DESeq2 + normCounts -----
    d289 = pd.read_csv(RAW / "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz", sep="\t")
    d289["gene"] = d289["Feature_name"].astype(str)
    d289 = d289.sort_values("baseMean", ascending=False).drop_duplicates("gene", keep="first")
    wt289 = ["OV.2808.RNA_normCounts", "OV.2810.RNA_normCounts", "OV.2812.RNA_normCounts"]
    ko289 = ["OV.2807.RNA_normCounts", "OV.2815.RNA_normCounts", "OV.2817.RNA_normCounts", "OV.2818.RNA_normCounts"]
    expr289 = d289.set_index("gene")[wt289 + ko289]
    log289 = np.log2(expr289 + PSEUDO)
    de289_w = welch_rows(log289, ko289, wt289)
    # also keep author DESeq2
    author = d289.set_index("gene")[["log2FoldChange", "pvalue", "padj", "baseMean"]].rename(
        columns={"log2FoldChange": "author_log2FC", "pvalue": "author_p", "padj": "author_q"}
    )
    de289_w = de289_w.merge(author, left_on="gene", right_index=True, how="left")
    de289_w.to_csv(TABLES / "de_GSE289287_Trop2KO_tumors_vs_WT.tsv.gz", sep="\t", index=False, compression="gzip")
    s289 = summarize_contrast("GSE289287_T47D_Trop2KO_xenograft", de289_w, log289, ko289, wt289, False, expr289)
    s289.update({
        "accession": "GSE289287",
        "agent": "none (TACSTD2 knockout)",
        "is_skb264": False,
        "is_adc_treatment": False,
        "tissue": "T-47D luminal breast xenografts",
        "lung": False,
        "label": "non-lung TROP2-loss analog (NOT ADC; not TNBC)",
        "time": "xenograft",
        "note": "TROP2 KO, not sacituzumab / not SKB264. Author DESeq2 + Welch on normCounts. n=4 KO vs 3 WT.",
        "author_TACSTD2_log2FC": float(author.loc["TACSTD2", "author_log2FC"]) if "TACSTD2" in author.index else None,
        "author_CLDN4_log2FC": float(author.loc["CLDN4", "author_log2FC"]) if "CLDN4" in author.index else None,
        "author_CLDN4_p": float(author.loc["CLDN4", "author_p"]) if "CLDN4" in author.index else None,
    })
    summaries.append(s289)
    catalog.append({**{k: s289[k] for k in ("accession", "agent", "is_adc_treatment", "tissue", "lung", "label", "note")}, "usable_adc_vs_ctrl": False})

    # ----- GSE235812 untreated breast P0/P1: correlation only -----
    g235 = load_gse235812()
    corr235 = None
    if g235 is not None:
        # GSM order from series matrix
        gsm = [f"GSM{i}" for i in range(7509978, 7510014)]
        g235 = g235.reindex(columns=[c for c in gsm if c in g235.columns])
        titles = []
        for i, c in enumerate(g235.columns):
            pair = i // 2
            passage = "P0" if i % 2 == 0 else "P1"
            titles.append(f"pair{pair:02d}_{passage}")
        g235.columns = titles
        log235 = np.log2(g235.astype(float) + PSEUDO)
        log235.to_csv(TABLES / "GSE235812_geneset_log2rpkm.tsv", sep="\t")
        if "TACSTD2" in log235.index:
            rows = []
            for g in ["CLDN4"] + [x for x in JUNCTION if x != "CLDN4"] + C4_IFN_MHCI + ["B2M", "TAP1"]:
                if g not in log235.index:
                    continue
                rho, p = stats.spearmanr(log235.loc["TACSTD2"], log235.loc[g])
                rows.append({"gene": g, "spearman_vs_TACSTD2": float(rho), "p": float(p), "n": int(log235.shape[1])})
            corr235 = pd.DataFrame(rows)
            corr235["q"] = multipletests(corr235["p"], method="fdr_bh")[1]
            corr235.to_csv(TABLES / "GSE235812_TACSTD2_spearman.tsv", sep="\t", index=False)
        catalog.append({
            "accession": "GSE235812",
            "agent": "none (untreated breast tumor vs PDX)",
            "is_adc_treatment": False,
            "tissue": "breast tumor P0 and matching PDX P1",
            "lung": False,
            "label": "non-lung untreated analog (NOT ADC)",
            "note": "Paper used SG in other experiments; deposited RNA is untreated P0/P1. No ADC vs control.",
            "usable_adc_vs_ctrl": False,
        })
    else:
        catalog.append({
            "accession": "GSE235812",
            "agent": "none",
            "is_adc_treatment": False,
            "usable_adc_vs_ctrl": False,
            "note": "RAW.tar missing; cannot score.",
        })

    # catalog leftovers (honest cannot)
    catalog.extend([
        {
            "accession": "GSE245459",
            "agent": "none (shTACSTD2)",
            "is_adc_treatment": False,
            "tissue": "SKOV3 ovarian",
            "lung": False,
            "label": "already covered as C4 TACSTD2-KD analog",
            "note": "Not ADC. See results/w200/C4_GSE245459/. Not re-analysed here.",
            "usable_adc_vs_ctrl": False,
        },
        {
            "accession": "GSE302284",
            "agent": "none in deposited omics (osimertinib / residual tissue)",
            "is_adc_treatment": False,
            "tissue": "EGFR-mutant NSCLC",
            "lung": True,
            "label": "not an SG treatment contrast",
            "note": "Paper mentions TROP2 CAR-T / SG context; GEO samples are LN/tumor and osi vs vehicle. Sibling slice owns this accession. Not scored as ADC vs control.",
            "usable_adc_vs_ctrl": False,
        },
        {
            "accession": "GSE278664",
            "agent": "none (prexasertib trial baseline)",
            "is_adc_treatment": False,
            "tissue": "HGSOC biopsies",
            "lung": False,
            "label": "SG paper, deposited RNA is CHK1i trial",
            "note": "Cannot test SG vs control. Treatment field is Prexasertib.",
            "usable_adc_vs_ctrl": False,
        },
        {
            "accession": "public GEO search",
            "agent": "sac-TMT / SKB264 / MK-2870",
            "is_adc_treatment": False,
            "tissue": None,
            "lung": None,
            "label": "no public expression series found",
            "note": "esearch SKB264 / sac-TMT / tirumotecan / datopotamab returned 0 GSE records (2026-08-16).",
            "usable_adc_vs_ctrl": False,
        },
    ])

    # key-gene table across ADC contrasts
    for s in summaries:
        for r in s["key_genes"]:
            key_all.append({
                "contrast": s["contrast"],
                "accession": s["accession"],
                "is_adc_treatment": s["is_adc_treatment"],
                **r,
            })
    key_df = pd.DataFrame(key_all)
    key_df.to_csv(TABLES / "key_genes_all_contrasts.tsv", sep="\t", index=False)

    set_rows = []
    for s in summaries:
        for sname, rec in s["sets"].items():
            set_rows.append({"contrast": s["contrast"], "accession": s["accession"], "is_adc_treatment": s["is_adc_treatment"], "set": sname, **rec})
        for sname, rec in s["scores"].items():
            set_rows.append({"contrast": s["contrast"], "accession": s["accession"], "is_adc_treatment": s["is_adc_treatment"], "set": sname + "_zscore", **rec})
    pd.DataFrame(set_rows).to_csv(TABLES / "geneset_stats.tsv", sep="\t", index=False)
    pd.DataFrame(catalog).to_csv(TABLES / "catalog.tsv", sep="\t", index=False)

    # compact JSON (drop huge key lists already in TSV)
    slim = []
    for s in summaries:
        slim.append({
            "contrast": s["contrast"],
            "accession": s["accession"],
            "agent": s["agent"],
            "is_skb264": False,
            "is_adc_treatment": s["is_adc_treatment"],
            "tissue": s["tissue"],
            "label": s["label"],
            "n_treat": s["n_treat"],
            "n_ctrl": s["n_ctrl"],
            "note": s["note"],
            "sets": s["sets"],
            "scores": s["scores"],
            "key_genes": s["key_genes"],
        })
    stats_out = {
        "disclaimer": "Public analogs only. None of these series is SKB264 / sac-TMT. IMMU132 is sacituzumab govitecan.",
        "catalog": catalog,
        "contrasts": slim,
        "gse235812_spearman": None if corr235 is None else corr235.to_dict(orient="records"),
    }
    (OUT / "key_stats.json").write_text(json.dumps(stats_out, indent=2, default=str) + "\n")

    # ================= figures =================
    plt.rcParams.update({"font.size": 9, "figure.facecolor": "white"})

    # fig_extra1 heatmap of key genes, ADC contrasts only
    adc_contrasts = [s for s in summaries if s["is_adc_treatment"]]
    genes_hm = ["TACSTD2", "CLDN4", "CLDN1", "CLDN7", "TJP1", "OCLN"] + C4_IFN_MHCI + ["B2M", "TAP1"]
    mat = []
    annot = []
    ylabs = []
    for s in adc_contrasts:
        kd = {r["gene"]: r for r in s["key_genes"]}
        row, an = [], []
        for g in genes_hm:
            if g in kd and kd[g]["log2FC"] is not None:
                row.append(kd[g]["log2FC"])
                p = kd[g]["p"]
                star = "" if p is None or p >= 0.05 else ("*" if p >= 0.01 else "**")
                an.append(f"{kd[g]['log2FC']:+.2f}{star}")
            else:
                row.append(np.nan)
                an.append("")
        mat.append(row)
        annot.append(an)
        ylabs.append(f"{s['accession']}\n{s['tissue']}\nn={s['n_treat']} vs {s['n_ctrl']}")
    fig, ax = plt.subplots(figsize=(12.5, 3.8))
    arr = np.array(mat, dtype=float)
    vmax = np.nanmax(np.abs(arr))
    vmax = 1.5 if not np.isfinite(vmax) or vmax < 0.5 else min(2.5, vmax)
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(genes_hm)))
    ax.set_xticklabels(genes_hm, rotation=45, ha="right")
    ax.set_yticks(range(len(ylabs)))
    ax.set_yticklabels(ylabs)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if annot[i][j]:
                ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.02, label="log2FC (ADC vs control)")
    ax.set_title("EXTRA analog Fig 1. IMMU132 (sacituzumab govitecan) vs control — not SKB264\n* p<0.05  ** p<0.01  Welch or paired t on log2(FPKM+1)")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra1_adc_keygene_heatmap.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra1_adc_keygene_heatmap.pdf")
    plt.close(fig)

    # fig_extra2 GSE312098 boxes
    fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.2))
    groups312 = {"ctrl": g312_ctrl, "IMMU132": g312_immu}
    cols = ["#9aa0a6", "#d97706"]
    for ax, g in zip(axes.flat, ["TACSTD2", "CLDN4", "CLDN7", "ISG15", "IFI27", "OAS2", "HLA-A", "B2M"]):
        box_groups(ax, log312, g, groups312, cols)
        if g in set(de312["gene"]):
            r = de312.loc[de312["gene"] == g].iloc[0]
            style_ax(ax, f"{g}  log2FC={r['log2FC']:+.2f}  p={r['p']:.3g}")
        else:
            style_ax(ax, g)
    fig.suptitle("EXTRA analog Fig 2. GSE312098 CX-1 CRC, IMMU132 vs control (n=3 vs 3), 2 d\nPublic SG analog. Not SKB264. Non-lung.", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra2_gse312098_cx1_immu132.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra2_gse312098_cx1_immu132.pdf")
    plt.close(fig)

    # fig_extra3 paired PDX
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.6))
    for ax, g in zip(axes, ["TACSTD2", "CLDN4", "ISG15", "HLA-A"]):
        if g not in log311.index:
            continue
        x0 = log311.loc[g, g311_ctrl].to_numpy(dtype=float)
        x1 = log311.loc[g, g311_immu].to_numpy(dtype=float)
        for a, b in zip(x0, x1):
            ax.plot([0, 1], [a, b], color="#888", lw=0.8)
            ax.scatter([0, 1], [a, b], c=["#9aa0a6", "#d97706"], s=28, zorder=3)
        r = de311.loc[de311["gene"] == g].iloc[0]
        style_ax(ax, f"{g}\nlog2FC={r['log2FC']:+.2f} p={r['p']:.3g}")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["ctrl", "IMMU132"])
        ax.set_ylabel("log2(FPKM+1)")
    fig.suptitle("EXTRA analog Fig 3. GSE311016 CRC PDX IMMU132 vs control, 5 paired models, day 29\nPublic SG analog. Not SKB264. Non-lung.", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra3_gse311016_crc_pdx.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra3_gse311016_crc_pdx.pdf")
    plt.close(fig)

    # fig_extra4 ESCC n=2
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.8))
    groups304 = {"ctrl n=3": g304_ctrl, "IMMU132 n=2": g304_immu}
    for ax, g in zip(axes.flat, ["TACSTD2", "CLDN4", "CLDN1", "ISG15", "IFIT1", "HLA-A"]):
        box_groups(ax, log304, g, groups304, cols)
        r = de304.loc[de304["gene"] == g].iloc[0]
        style_ax(ax, f"{g}  log2FC={r['log2FC']:+.2f}  p={r['p']:.3g}")
    fig.suptitle("EXTRA analog Fig 4. GSE304294 KYSE30 ESCC IMMU132 vs control\nHonest n=2 vs 3, 1 day. Public SG analog. Not SKB264. Underpowered.", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra4_gse304294_escc.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra4_gse304294_escc.pdf")
    plt.close(fig)

    # fig_extra5 TROP2-loss (NOT ADC)
    fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.2))
    groups334 = {"WT": wt, "Trop2 KO": ko}
    for ax, g in zip(axes[0], ["TACSTD2", "CLDN4", "CLDN1", "ISG15"]):
        box_groups(ax, log334, g, groups334, ["#6b8f71", "#6a4c93"])
        if g in set(de334["gene"]):
            r = de334.loc[de334["gene"] == g].iloc[0]
            style_ax(ax, f"4T1 {g}  log2FC={r['log2FC']:+.2f} p={r['p']:.3g}")
        else:
            style_ax(ax, f"4T1 {g}")
    groups289 = {"WT": wt289, "Trop2 KO": ko289}
    for ax, g in zip(axes[1], ["TACSTD2", "CLDN4", "ISG15", "HLA-A"]):
        box_groups(ax, log289, g, groups289, ["#6b8f71", "#6a4c93"])
        r = de289_w.loc[de289_w["gene"] == g].iloc[0]
        ap = r["author_p"] if "author_p" in r.index else np.nan
        style_ax(ax, f"T-47D {g}  log2FC={r['log2FC']:+.2f} p={r['p']:.3g}\nauthor DESeq2 p={ap:.3g}")
    fig.suptitle("EXTRA analog Fig 5. TROP2-loss analogs — NOT ADC, not SKB264\nTop: GSE334497 4T1 mouse TNBC KO vs WT (n=5 vs 5). Bottom: GSE289287 T-47D luminal xenografts (n=4 vs 3).", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra5_trop2loss_not_adc.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra5_trop2loss_not_adc.pdf")
    plt.close(fig)

    # fig_extra6 untreated correlation
    if g235 is not None and "TACSTD2" in log235.index and "CLDN4" in log235.index:
        fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.5))
        p0 = [c for c in log235.columns if c.endswith("_P0")]
        p1 = [c for c in log235.columns if c.endswith("_P1")]
        for ax, cols, lab in (
            (axes[0], list(log235.columns), "all n=36"),
            (axes[1], p0, "patient P0 n=18"),
            (axes[2], p1, "PDX P1 n=18"),
        ):
            x = log235.loc["TACSTD2", cols].to_numpy(dtype=float)
            y = log235.loc["CLDN4", cols].to_numpy(dtype=float)
            ax.scatter(x, y, s=22, c="#3b6d9a")
            rho, p = stats.spearmanr(x, y)
            style_ax(ax, f"{lab}\nρ={rho:.2f} p={p:.3g}")
            ax.set_xlabel("log2(TACSTD2 RPKM+1)")
            ax.set_ylabel("log2(CLDN4 RPKM+1)")
        fig.suptitle("EXTRA analog Fig 6. GSE235812 untreated breast tumor/PDX — no ADC contrast\nNot SKB264. Not sacituzumab-treated RNA.", fontsize=11)
        fig.tight_layout()
        fig.savefig(FIGS / "fig_extra6_gse235812_untreated_corr.png", dpi=DPI)
        fig.savefig(FIGS / "fig_extra6_gse235812_untreated_corr.pdf")
        plt.close(fig)

    # fig_extra7 set-score forest for ADC only
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    y = 0
    yticks, ylabs = [], []
    for s in adc_contrasts:
        for sname, color in (("JUNCTION", "#3b6d9a"), ("C4_IFN_MHCI", "#c0392b"), ("APM", "#d97706"), ("IFN", "#6a4c93")):
            rec = s["sets"][sname]
            med = rec.get("median_log2FC_set")
            p = rec.get("mw_p")
            if med is None:
                continue
            ax.scatter(med, y, c=color, s=40, zorder=3)
            ax.plot([0, med], [y, y], color=color, lw=1.4)
            labp = "NA" if p is None else f"{p:.3g}"
            yticks.append(y)
            ylabs.append(f"{s['accession']} {sname}  med={med:+.2f} MW p={labp}")
            y -= 1
        y -= 0.4
    ax.axvline(0, color="#444", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabs, fontsize=8)
    ax.set_xlabel("median log2FC of set (ADC vs control)")
    style_ax(ax, "EXTRA analog Fig 7. Gene-set median log2FC after IMMU132 — not SKB264")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra7_adc_geneset_forest.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra7_adc_geneset_forest.pdf")
    plt.close(fig)

    print("wrote", OUT)
    for s in summaries:
        c4 = next((r for r in s["key_genes"] if r["gene"] == "CLDN4"), None)
        print(s["contrast"], "CLDN4", c4, "C4set", s["sets"]["C4_IFN_MHCI"]["direction"], s["sets"]["C4_IFN_MHCI"].get("mw_p"))


if __name__ == "__main__":
    main()
