#!/usr/bin/env python3
"""Co-perturbation and dual-readout sweep: TACSTD2 and CLDN4 versus an IFN/immune score.

Every log2FC and p-value is computed from a public GEO matrix (or taken from an
author table when the deposit is a contrast, not a sample matrix). No value is filled
in when a gene or a sample label is missing.

Call rule, fixed before looking at the results
----------------------------------------------
Partial dependence here means the immune/IFN score and the two epithelial genes move
as a pair.

Genetic loss or TROP2-ADC (perturb gene A, read partner B and the IFN score):
  perturbation_failed  A is not down (log2FC > -0.5)
  null                 partner is not significant, or the IFN score is not
  supports             A is down, partner is down (p<0.05 and |log2FC|>=0.25),
                       and the IFN score moves (p<0.05 and |delta|>=0.25)
  opposite             A is down, partner is up by that same threshold, and the
                       IFN score moves

IFN or dsRNA stimulus (both genes are readouts; the score must rise):
  stimulus_failed      IFN score does not increase (delta<=0 or p>=0.05)
  supports             score rises, and both genes move the same direction
                       (each p<0.05 and |log2FC|>=0.25)
  opposite             score rises, and the two genes move in opposite directions
                       at that threshold
  null                 score rises, but the pair does not co-move

Immune treatment uses the stimulus rule but does not require the score to rise.
A flat score is null.

Unreplicated contrasts (n=1) are reported and are not given a three-way call.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RAW = Path("/tmp/sweep/raw")
OUT = Path("/workspace/results/tacstd2_cldn4_partial")
OUT.mkdir(parents=True, exist_ok=True)

HUMAN_IFN = [
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "IFI6", "IFI27", "IFI44", "IFI44L", "STAT1", "IRF1", "IRF7", "B2M",
    "HLA-A", "HLA-B", "HLA-C", "TAP1", "PSMB8", "PSMB9", "CD274", "CXCL10", "CXCL9",
]
MOUSE_IFN = [
    "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2", "Oas3",
    "Ifi27l2a", "Ifi44", "Stat1", "Irf1", "Irf7", "B2m", "H2-K1", "H2-D1",
    "Tap1", "Psmb8", "Psmb9", "Cd274", "Cxcl10", "Cxcl9",
]
MIN_IFN = 8
LFC_FLOOR = 0.25
P_CUT = 0.05
PERT_DOWN = -0.5


def welch(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.nanstd(a) == 0 and np.nanstd(b) == 0:
        return 1.0 if np.nanmean(a) == np.nanmean(b) else np.nan
    return float(stats.ttest_ind(a, b, equal_var=False, alternative="two-sided").pvalue)


def paired(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 2:
        return np.nan
    return float(stats.ttest_rel(a, b, alternative="two-sided").pvalue)


def onesample(deltas):
    d = np.asarray(deltas, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) < 2:
        return np.nan
    if np.nanstd(d) == 0:
        return np.nan
    return float(stats.ttest_1samp(d, 0.0, alternative="two-sided").pvalue)


def partial_spearman(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]
    n = len(x)
    if n < 8:
        return n, np.nan, np.nan
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return n, float(r), float(p)


def strip_ens(s):
    return str(s).split(".")[0]


def log2p1(v):
    return np.log2(np.asarray(v, dtype=float) + 1.0)


def gene_row_symbol(df, symbol, symbol_col):
    hit = df.loc[df[symbol_col] == symbol]
    if hit.empty:
        return None
    return hit.iloc[0]


def mean_log_groups(values_by_sample, treat_cols, ctrl_cols, already_log):
    def arr(cols):
        raw = np.array([float(values_by_sample[c]) for c in cols], dtype=float)
        return raw if already_log else log2p1(raw)

    a = arr(treat_cols)
    b = arr(ctrl_cols)
    return float(np.mean(a) - np.mean(b)), welch(a, b), a, b


def score_matrix(log_df, genes):
    """log_df: genes x samples, already log2."""
    present = [g for g in genes if g in log_df.index]
    if len(present) < MIN_IFN:
        return None, present
    return log_df.loc[present].mean(axis=0), present


def classify_genetic(pert_lfc, pert_p, partner_lfc, partner_p, ifn_delta, ifn_p):
    if pert_lfc is None or not np.isfinite(pert_lfc) or pert_lfc > PERT_DOWN:
        return "perturbation_failed"
    partner_sig = (
        partner_lfc is not None
        and np.isfinite(partner_lfc)
        and partner_p is not None
        and np.isfinite(partner_p)
        and partner_p < P_CUT
        and abs(partner_lfc) >= LFC_FLOOR
    )
    ifn_sig = (
        ifn_delta is not None
        and np.isfinite(ifn_delta)
        and ifn_p is not None
        and np.isfinite(ifn_p)
        and ifn_p < P_CUT
        and abs(ifn_delta) >= LFC_FLOOR
    )
    if not partner_sig or not ifn_sig:
        return "null"
    if partner_lfc < 0:
        return "supports"
    if partner_lfc > 0:
        return "opposite"
    return "null"


def classify_stimulus(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p, require_ifn_up):
    score_up = (
        ifn_delta is not None
        and np.isfinite(ifn_delta)
        and ifn_p is not None
        and np.isfinite(ifn_p)
        and ifn_p < P_CUT
        and ifn_delta >= LFC_FLOOR
    )
    score_moved = (
        ifn_delta is not None
        and np.isfinite(ifn_delta)
        and ifn_p is not None
        and np.isfinite(ifn_p)
        and ifn_p < P_CUT
        and abs(ifn_delta) >= LFC_FLOOR
    )
    if require_ifn_up and not score_up:
        return "stimulus_failed"
    if not score_moved:
        return "null"

    def moved(lfc, p):
        return (
            lfc is not None
            and p is not None
            and np.isfinite(lfc)
            and np.isfinite(p)
            and p < P_CUT
            and abs(lfc) >= LFC_FLOOR
        )

    tm, cm = moved(t_lfc, t_p), moved(c_lfc, c_p)
    if tm and cm:
        return "supports" if t_lfc * c_lfc > 0 else "opposite"
    return "null"


def ifn_direction(delta):
    if delta is None or not np.isfinite(delta) or abs(delta) < LFC_FLOOR:
        return "flat"
    return "up" if delta > 0 else "down"


def empty_row(**kwargs):
    base = {
        "accession": "",
        "contrast": "",
        "design_class": "",
        "species": "",
        "n_treat": "",
        "n_ctrl": "",
        "test": "",
        "tacstd2_log2fc": np.nan,
        "tacstd2_p": np.nan,
        "cldn4_log2fc": np.nan,
        "cldn4_p": np.nan,
        "ifn_score_delta": np.nan,
        "ifn_score_p": np.nan,
        "ifn_n_genes": "",
        "ifn_direction": "",
        "partial_n": "",
        "partial_rho_tacstd2_ifn_given_cldn4": np.nan,
        "partial_p_tacstd2": np.nan,
        "partial_rho_cldn4_ifn_given_tacstd2": np.nan,
        "partial_p_cldn4": np.nan,
        "call": "",
        "note": "",
    }
    base.update(kwargs)
    return base


def attach_partial(row, tac_s, cld_s, ifn_s):
    n1, r1, p1 = partial_spearman(tac_s, ifn_s, cld_s)
    n2, r2, p2 = partial_spearman(cld_s, ifn_s, tac_s)
    row["partial_n"] = int(n1)
    row["partial_rho_tacstd2_ifn_given_cldn4"] = r1
    row["partial_p_tacstd2"] = p1
    row["partial_rho_cldn4_ifn_given_tacstd2"] = r2
    row["partial_p_cldn4"] = p2
    return row


def contrast_from_logs(log_df, treat, ctrl, ifn_genes, species, meta, require_ifn_up, genetic_gene):
    """log_df indexed by gene symbol, columns are samples, values already log2."""
    missing = [g for g in ("TACSTD2", "CLDN4", "Tacstd2", "Cldn4") if False]
    t_sym = "TACSTD2" if species == "human" else "Tacstd2"
    c_sym = "CLDN4" if species == "human" else "Cldn4"
    if t_sym not in log_df.index or c_sym not in log_df.index:
        row = empty_row(**meta)
        row["call"] = "not_evaluable"
        row["note"] = f"missing {t_sym if t_sym not in log_df.index else c_sym}"
        return row, []
    tac = log_df.loc[t_sym, treat + ctrl].astype(float)
    cld = log_df.loc[c_sym, treat + ctrl].astype(float)
    # if duplicate index, mean them
    if isinstance(tac, pd.DataFrame):
        tac = tac.mean(axis=0)
    if isinstance(cld, pd.DataFrame):
        cld = cld.mean(axis=0)
    t_lfc = float(tac[treat].mean() - tac[ctrl].mean())
    c_lfc = float(cld[treat].mean() - cld[ctrl].mean())
    t_p = welch(tac[treat], tac[ctrl])
    c_p = welch(cld[treat], cld[ctrl])
    score, present = score_matrix(log_df, ifn_genes)
    gene_rows = []
    if score is None:
        row = empty_row(**meta, tacstd2_log2fc=t_lfc, tacstd2_p=t_p, cldn4_log2fc=c_lfc, cldn4_p=c_p)
        row["call"] = "not_evaluable"
        row["ifn_n_genes"] = len(present)
        row["note"] = f"IFN panel present {len(present)} < {MIN_IFN}"
        return row, gene_rows
    ifn_delta = float(score[treat].mean() - score[ctrl].mean())
    ifn_p = welch(score[treat], score[ctrl])
    for g in present:
        gvec = log_df.loc[g]
        if isinstance(gvec, pd.DataFrame):
            gvec = gvec.mean(axis=0)
        glfc = float(gvec[treat].mean() - gvec[ctrl].mean())
        gene_rows.append({"accession": meta["accession"], "contrast": meta["contrast"], "gene": g, "log2fc": glfc})
    if genetic_gene == "TACSTD2":
        call = classify_genetic(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p)
    elif genetic_gene == "CLDN4":
        call = classify_genetic(c_lfc, c_p, t_lfc, t_p, ifn_delta, ifn_p)
    else:
        call = classify_stimulus(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p, require_ifn_up)
    row = empty_row(
        **meta,
        n_treat=len(treat),
        n_ctrl=len(ctrl),
        tacstd2_log2fc=t_lfc,
        tacstd2_p=t_p,
        cldn4_log2fc=c_lfc,
        cldn4_p=c_p,
        ifn_score_delta=ifn_delta,
        ifn_score_p=ifn_p,
        ifn_n_genes=len(present),
        ifn_direction=ifn_direction(ifn_delta),
        call=call,
    )
    attach_partial(row, tac[treat + ctrl], cld[treat + ctrl], score[treat + ctrl])
    return row, gene_rows


def matrix_from_symbol_frame(df, symbol_col, sample_cols, already_log):
    sub = df[[symbol_col] + sample_cols].copy()
    sub[symbol_col] = sub[symbol_col].astype(str)
    for c in sample_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    vals = sub.groupby(symbol_col, as_index=True)[sample_cols].mean()
    if already_log:
        return vals
    return np.log2(vals + 1.0)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_gse334497():
    df = pd.read_csv(RAW / "GSE334497_normalized_counts.csv.gz", index_col=0)
    df.index = df.index.map(strip_ens)
    # column names encode KO vs the remaining samples
    ko = [c for c in df.columns if "KO" in str(c).upper()]
    wt = [c for c in df.columns if c not in ko]
    keep = {
        "ENSMUSG00000051397": "Tacstd2",
        "ENSMUSG00000047501": "Cldn4",
        "ENSMUSG00000035692": "Isg15",
        "ENSMUSG00000034459": "Ifit1",
        "ENSMUSG00000045932": "Ifit2",
        "ENSMUSG00000074896": "Ifit3",
        "ENSMUSG00000000386": "Mx1",
        "ENSMUSG00000023341": "Mx2",
        "ENSMUSG00000052776": "Oas1a",
        "ENSMUSG00000032690": "Oas2",
        "ENSMUSG00000032661": "Oas3",
        "ENSMUSG00000079017": "Ifi27l2a",
        "ENSMUSG00000028037": "Ifi44",
        "ENSMUSG00000026104": "Stat1",
        "ENSMUSG00000018899": "Irf1",
        "ENSMUSG00000025498": "Irf7",
        "ENSMUSG00000060802": "B2m",
        "ENSMUSG00000061232": "H2-K1",
        "ENSMUSG00000073411": "H2-D1",
        "ENSMUSG00000037321": "Tap1",
        "ENSMUSG00000024338": "Psmb8",
        "ENSMUSG00000096727": "Psmb9",
        "ENSMUSG00000016496": "Cd274",
        "ENSMUSG00000034855": "Cxcl10",
        "ENSMUSG00000029417": "Cxcl9",
    }
    present = df.loc[df.index.intersection(keep.keys())].copy()
    present.index = present.index.map(keep)
    log_df = np.log2(present.astype(float) + 1.0)
    # QC: KO columns must be the low Tacstd2 group
    if log_df.loc["Tacstd2", ko].mean() >= log_df.loc["Tacstd2", wt].mean():
        raise SystemExit("GSE334497 KO columns are not Tacstd2-low; mapping refused")
    meta = dict(
        accession="GSE334497",
        contrast="4T1 Trop2 KO vs WT tumor",
        design_class="genetic_ko",
        species="mouse",
        test="Welch on log2(norm+1)",
        note="KO columns are those whose names contain KO (5); remaining 5 are WT/control. Checked: Tacstd2 is lower in the KO columns.",
    )
    return contrast_from_logs(log_df, ko, wt, MOUSE_IFN, "mouse", meta, False, "TACSTD2")


def load_gse289287():
    df = pd.read_csv(RAW / "GSE289287_Trop2KO.tsv.gz", sep="\t")
    wt = ["OV.2808.RNA_normCounts", "OV.2810.RNA_normCounts", "OV.2812.RNA_normCounts"]
    ko = ["OV.2807.RNA_normCounts", "OV.2815.RNA_normCounts", "OV.2817.RNA_normCounts", "OV.2818.RNA_normCounts"]
    sym = {
        "TACSTD2": "TACSTD2", "CLDN4": "CLDN4", "ISG15": "ISG15", "IFIT1": "IFIT1",
        "IFIT2": "IFIT2", "IFIT3": "IFIT3", "MX1": "MX1", "MX2": "MX2", "OAS1": "OAS1",
        "OAS2": "OAS2", "OAS3": "OAS3", "IFI6": "IFI6", "IFI27": "IFI27", "IFI44": "IFI44",
        "IFI44L": "IFI44L", "STAT1": "STAT1", "IRF1": "IRF1", "IRF7": "IRF7", "B2M": "B2M",
        "HLA-A": "HLA-A", "HLA-B": "HLA-B", "HLA-C": "HLA-C", "TAP1": "TAP1",
        "PSMB8": "PSMB8", "PSMB9": "PSMB9", "CD274": "CD274", "CXCL10": "CXCL10", "CXCL9": "CXCL9",
    }
    sub = df[df["Feature_name"].isin(sym)].copy()
    log_df = matrix_from_symbol_frame(sub, "Feature_name", wt + ko, already_log=False)
    meta = dict(
        accession="GSE289287",
        contrast="T-47D Trop-2 KO vs WT xenograft",
        design_class="genetic_ko",
        species="human",
        test="Welch on log2(author normCounts+1)",
        note="WT animals 2808/2810/2812; KO animals 2807/2815/2817/2818, from GSM titles. Author DESeq2 columns are in the same file and are not the p used for the call.",
    )
    row, genes = contrast_from_logs(log_df, ko, wt, HUMAN_IFN, "human", meta, False, "TACSTD2")
    # attach author statistics for the two genes
    for gene, prefix in (("TACSTD2", "tacstd2"), ("CLDN4", "cldn4")):
        hit = df.loc[df["Feature_name"] == gene]
        if not hit.empty:
            row["note"] += (
                f" Author DESeq2 {gene} log2FC={float(hit.iloc[0]['log2FoldChange']):.4f},"
                f" padj={float(hit.iloc[0]['padj']):.3g}."
            )
    return row, genes


def load_gse245459():
    df = pd.read_csv(RAW / "GSE245459_fpkm.anno.txt.gz", sep="\t", usecols=lambda c: c in {
        "GeneName", "shNC1", "shNC2", "shNC3", "shNCDDP1", "shNCDDP2", "shNCDDP3",
        "sh1", "sh2", "sh3", "shDDP1", "shDDP2", "shDDP3",
    })
    rows = []
    genes = []
    contrasts = [
        ("SKOV3 shTACSTD2 vs shNC, no drug", ["sh1", "sh2", "sh3"], ["shNC1", "shNC2", "shNC3"]),
        ("SKOV3 shTACSTD2 vs shNC, both on cisplatin", ["shDDP1", "shDDP2", "shDDP3"], ["shNCDDP1", "shNCDDP2", "shNCDDP3"]),
    ]
    for name, treat, ctrl in contrasts:
        log_df = matrix_from_symbol_frame(df, "GeneName", treat + ctrl, already_log=False)
        meta = dict(
            accession="GSE245459",
            contrast=name,
            design_class="genetic_kd",
            species="human",
            test="Welch on log2(FPKM+1)",
            note="Sample columns are labeled in the FPKM file. Cisplatin arms are a second contrast, not pooled with the untreated arm.",
        )
        row, g = contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, False, "TACSTD2")
        rows.append(row)
        genes.extend(g)
    return rows, genes


def load_gse207704():
    df = pd.read_csv(RAW / "GSE207704_CLDN4_RNAseq.txt.gz", sep="\t")
    cols = {
        "MCF7_CLDN4KO_FPKM (fpkm)": "MCF7_KO",
        "MCF7_WT_FPKM (fpkm)": "MCF7_WT",
        "T47D_CLDN4KO_FPKM (fpkm)": "T47D_KO",
        "T47D_WT_FPKM (fpkm)": "T47D_WT",
    }
    sub = df[["gene_short_name"] + list(cols)].rename(columns=cols)
    log_df = matrix_from_symbol_frame(sub, "gene_short_name", list(cols.values()), already_log=False)
    # one collapsed FPKM per line; test is across the two lines
    lines = [("MCF7_KO", "MCF7_WT"), ("T47D_KO", "T47D_WT")]
    t_sym, c_sym = "TACSTD2", "CLDN4"
    if t_sym not in log_df.index or c_sym not in log_df.index:
        meta = dict(accession="GSE207704", contrast="CLDN4-/- vs WT, MCF7 and T47D", design_class="genetic_ko", species="human", test="two cell lines", note="gene missing")
        row = empty_row(**meta, call="not_evaluable")
        return row, []
    t_deltas, c_deltas, ifn_deltas = [], [], []
    score, present = score_matrix(log_df, HUMAN_IFN)
    gene_rows = []
    for ko, wt in lines:
        t_deltas.append(float(log_df.loc[t_sym, ko] - log_df.loc[t_sym, wt]))
        c_deltas.append(float(log_df.loc[c_sym, ko] - log_df.loc[c_sym, wt]))
        if score is not None:
            ifn_deltas.append(float(score[ko] - score[wt]))
    t_lfc, c_lfc = float(np.mean(t_deltas)), float(np.mean(c_deltas))
    t_p, c_p = onesample(t_deltas), onesample(c_deltas)
    ifn_delta = float(np.mean(ifn_deltas)) if ifn_deltas else np.nan
    ifn_p = onesample(ifn_deltas) if ifn_deltas else np.nan
    for g in present:
        gl = [float(log_df.loc[g, ko] - log_df.loc[g, wt]) for ko, wt in lines]
        gene_rows.append({"accession": "GSE207704", "contrast": "CLDN4-/- vs WT, two lines", "gene": g, "log2fc": float(np.mean(gl))})
    call = classify_genetic(c_lfc, c_p, t_lfc, t_p, ifn_delta, ifn_p)
    note = (
        "Deposit is one FPKM per genotype per line; the two GEO replicates are already collapsed. "
        f"n=2 lines. TACSTD2 line log2FC={t_deltas[0]:.3f},{t_deltas[1]:.3f}. "
        f"CLDN4 line log2FC={c_deltas[0]:.3f},{c_deltas[1]:.3f}. "
        f"IFN-score line delta={ifn_deltas[0]:.3f},{ifn_deltas[1]:.3f}."
    )
    row = empty_row(
        accession="GSE207704",
        contrast="CLDN4-/- vs WT, MCF7 and T47D",
        design_class="genetic_ko",
        species="human",
        n_treat=2,
        n_ctrl=2,
        test="one-sample t of two line-level log2(FPKM+1) differences",
        tacstd2_log2fc=t_lfc,
        tacstd2_p=t_p,
        cldn4_log2fc=c_lfc,
        cldn4_p=c_p,
        ifn_score_delta=ifn_delta,
        ifn_score_p=ifn_p,
        ifn_n_genes=len(present),
        ifn_direction=ifn_direction(ifn_delta),
        call=call,
        note=note,
    )
    return row, gene_rows


def load_gse50927():
    df = pd.read_csv(RAW / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    want = set(MOUSE_IFN + ["Tacstd2", "Cldn4"])
    sub = df[df["Marker.Symbol"].isin(want)].copy()
    def grab(sym):
        hit = sub.loc[sub["Marker.Symbol"] == sym]
        if hit.empty:
            return None
        return hit.iloc[0]
    rows_out = []
    cld = grab("Cldn4")
    tac = grab("Tacstd2")
    ifn_lfcs = []
    gene_rows = []
    for g in MOUSE_IFN:
        hit = grab(g)
        if hit is None:
            continue
        ifn_lfcs.append(float(hit["logFC"]))
        gene_rows.append({"accession": "GSE50927", "contrast": "Cldn4 KO vs WT lung, no VILI", "gene": g, "log2fc": float(hit["logFC"])})
    ifn_delta = float(np.mean(ifn_lfcs)) if ifn_lfcs else np.nan
    ifn_p = onesample(ifn_lfcs) if len(ifn_lfcs) >= 2 else np.nan
    note = (
        "Author edgeR table GSE50927_Cldn4lungWTvsKOgenes. GEO has one WT and one KO library without VILI, "
        "so this is unreplicated and is not given a supports/null/opposite call. "
        "edgeR logFC is used as deposited. The IFN value is the mean of those gene logFCs, not a sample-level test."
    )
    if cld is not None:
        note += f" Cldn4 logFC={float(cld['logFC']):.4f}, author P={float(cld['PValue']):.3g}, FDR={float(cld['FDR']):.3g}."
    row = empty_row(
        accession="GSE50927",
        contrast="Cldn4 KO lung vs WT, no VILI",
        design_class="unreplicated",
        species="mouse",
        n_treat=1,
        n_ctrl=1,
        test="author edgeR logFC; no replicate p used for a call",
        tacstd2_log2fc=float(tac["logFC"]) if tac is not None else np.nan,
        tacstd2_p=float(tac["PValue"]) if tac is not None else np.nan,
        cldn4_log2fc=float(cld["logFC"]) if cld is not None else np.nan,
        cldn4_p=float(cld["PValue"]) if cld is not None else np.nan,
        ifn_score_delta=ifn_delta,
        ifn_score_p=np.nan,
        ifn_n_genes=len(ifn_lfcs),
        ifn_direction=ifn_direction(ifn_delta),
        call="unreplicated",
        note=note + f" Mean of deposited IFN-gene logFC={ifn_delta:.4f} across {len(ifn_lfcs)} genes. That gene-wise mean is not a sample-level p-value and is not used as one.",
    )
    return row, gene_rows


def load_gse304294():
    df = pd.read_csv(RAW / "GSE304294_gene_fpkm.txt.gz", sep="\t")
    # OX2 is the only block with 2 replicates and matches IMMU132 (the only n=2 GSM group).
    # OX1 is the n=3 block placed immediately before it, matching Control submitted before IMMU132.
    ctrl = ["OX1_1", "OX1_2", "OX1_3"]
    treat = ["OX2_1", "OX2_2"]
    log_df = matrix_from_symbol_frame(df, "gene_name", ctrl + treat, already_log=False)
    meta = dict(
        accession="GSE304294",
        contrast="KYSE30 IMMU132 vs control, 1 day",
        design_class="trop2_adc",
        species="human",
        test="Welch on log2(FPKM+1)",
        note="OX2 (n=2) is IMMU132 because it is the only n=2 block and IMMU132 is the only n=2 group. OX1 is the preceding n=3 block (Control in submission order). IACS and combination are not this contrast.",
    )
    return contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, False, "TACSTD2")


def load_gse311016():
    df = pd.read_csv(RAW / "GSE311016_gene_fpkm.txt.gz", sep="\t", encoding="utf-16")
    df.columns = [str(c).replace("\ufeff", "") for c in df.columns]
    pairs = ["114", "36", "82", "83", "196"]
    treat = [f"T_{p}" for p in pairs]
    ctrl = [f"C_{p}" for p in pairs]
    for c in treat + ctrl:
        if c not in df.columns:
            raise SystemExit(f"GSE311016 missing {c}: {list(df.columns)[:15]}")
    idcol = "gene_name"
    sub = df[[idcol] + treat + ctrl].copy()
    sub[idcol] = sub[idcol].astype(str)
    log_df = np.log2(sub.groupby(idcol)[treat + ctrl].mean(numeric_only=True) + 1.0)
    t_sym, c_sym = "TACSTD2", "CLDN4"
    meta_base = dict(
        accession="GSE311016",
        contrast="CRC PDX IMMU132 vs paired control, day 29",
        design_class="trop2_adc",
        species="human",
        test="paired t on log2(FPKM+1), 5 PDX pairs",
    )
    if t_sym not in log_df.index or c_sym not in log_df.index:
        row = empty_row(**meta_base, call="not_evaluable", note="gene missing")
        return row, []
    def paired_lfc(sym):
        a = log_df.loc[sym, treat].astype(float).to_numpy()
        b = log_df.loc[sym, ctrl].astype(float).to_numpy()
        return float(np.mean(a - b)), paired(a, b), a, b
    t_lfc, t_p, ta, tb = paired_lfc(t_sym)
    c_lfc, c_p, ca, cb = paired_lfc(c_sym)
    score, present = score_matrix(log_df, HUMAN_IFN)
    gene_rows = []
    if score is None:
        row = empty_row(**meta_base, call="not_evaluable", ifn_n_genes=len(present), note="IFN panel short")
        return row, gene_rows
    sa = score[treat].to_numpy()
    sb = score[ctrl].to_numpy()
    ifn_delta = float(np.mean(sa - sb))
    ifn_p = paired(sa, sb)
    for g in present:
        gene_rows.append({
            "accession": "GSE311016",
            "contrast": meta_base["contrast"],
            "gene": g,
            "log2fc": float(np.mean(log_df.loc[g, treat].to_numpy() - log_df.loc[g, ctrl].to_numpy())),
        })
    call = classify_genetic(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p)
    row = empty_row(
        **meta_base,
        n_treat=5,
        n_ctrl=5,
        tacstd2_log2fc=t_lfc,
        tacstd2_p=t_p,
        cldn4_log2fc=c_lfc,
        cldn4_p=c_p,
        ifn_score_delta=ifn_delta,
        ifn_score_p=ifn_p,
        ifn_n_genes=len(present),
        ifn_direction=ifn_direction(ifn_delta),
        call=call,
        note="Human genes in PDX RNA-seq. T_/C_ suffixes match PDX ids 114, 36, 82, 83, 196. Paired by PDX. Immune infiltrate is not in the human alignment.",
    )
    return row, gene_rows


def load_gse274940():
    df = pd.read_csv(RAW / "GSE274940_raw_counts.csv.gz")
    df["symbol"] = df["ALIAS"].astype(str)
    samples = ["WT1", "WT2", "WT3", "KO1", "KO2", "KO3"]
    counts = df.groupby("symbol")[samples].sum(numeric_only=True)
    lib = counts.sum(axis=0)
    log_df = np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)
    meta = dict(
        accession="GSE274940",
        contrast="EpH4 Cldn-null vs WT",
        design_class="genetic_ko",
        species="mouse",
        test="Welch on log2(CPM+1)",
        note="Pan-claudin null line, not a Cldn4-only knockout. Call is perturbation_failed when Cldn4 itself is not down.",
    )
    return contrast_from_logs(log_df, ["KO1", "KO2", "KO3"], ["WT1", "WT2", "WT3"], MOUSE_IFN, "mouse", meta, False, "CLDN4")


def load_fpkm_or_count_named(path, symbol_col, sample_map, already_log, raw_counts, meta, ifn_genes, species, treat, ctrl, genetic, require_ifn_up, sep=","):
    df = pd.read_csv(path, sep=sep)
    if symbol_col not in df.columns:
        # try first columns
        raise SystemExit(f"{path.name} missing {symbol_col}: {list(df.columns)[:8]}")
    cols = treat + ctrl
    if raw_counts:
        # CPM from all genes, then log
        num = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        lib = num.sum(axis=0)
        logged = np.log2(num.div(lib, axis=1) * 1e6 + 1.0)
        logged[symbol_col] = df[symbol_col].astype(str).values
        log_df = logged.groupby(symbol_col)[cols].mean()
    else:
        log_df = matrix_from_symbol_frame(df, symbol_col, cols, already_log=already_log)
    return contrast_from_logs(log_df, treat, ctrl, ifn_genes, species, meta, require_ifn_up, genetic)


def load_gse156295():
    df = pd.read_csv(RAW / "GSE156295_count_hg38.txt.gz", sep="\t")
    out_rows, out_genes = [], []
    contrasts = [
        ("A549 type I IFN vs mock", ["A549-IFN-1", "A549-IFN-2"], ["A549-Mock-1", "A549-Mock-2"]),
        ("HTBE type I IFN vs mock", ["HTBE-IFN-1", "HTBE-IFN-2"], ["HTBE-Mock-1", "HTBE-Mock-2"]),
    ]
    for name, treat, ctrl in contrasts:
        num = df[treat + ctrl].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        lib = num.sum(axis=0)
        logged = np.log2(num.div(lib, axis=1) * 1e6 + 1.0)
        logged.insert(0, "Gene", df["Gene"].astype(str).values)
        log_df = logged.groupby("Gene")[treat + ctrl].mean()
        meta = dict(
            accession="GSE156295",
            contrast=name,
            design_class="ifn_stimulus",
            species="human",
            test="Welch on log2(CPM+1)",
            note="Columns are named in the count file. n=2 vs 2.",
        )
        row, genes = contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, True, None)
        out_rows.append(row)
        out_genes.extend(genes)
    return out_rows, out_genes


def load_gse184456():
    df = pd.read_csv(RAW / "GSE184456_htseq_rawCounts.txt.gz", sep="\t")
    df["ens"] = df["ENSEMBL_ID"].map(strip_ens)
    # map needed ids
    want = {
        "ENSG00000184292": "TACSTD2", "ENSG00000189143": "CLDN4",
        "ENSG00000187608": "ISG15", "ENSG00000185745": "IFIT1", "ENSG00000119922": "IFIT2",
        "ENSG00000119917": "IFIT3", "ENSG00000157601": "MX1", "ENSG00000183486": "MX2",
        "ENSG00000089127": "OAS1", "ENSG00000111335": "OAS2", "ENSG00000111331": "OAS3",
        "ENSG00000126709": "IFI6", "ENSG00000165949": "IFI27", "ENSG00000137965": "IFI44",
        "ENSG00000137959": "IFI44L", "ENSG00000115415": "STAT1", "ENSG00000125347": "IRF1",
        "ENSG00000185507": "IRF7", "ENSG00000166710": "B2M", "ENSG00000206503": "HLA-A",
        "ENSG00000234745": "HLA-B", "ENSG00000204525": "HLA-C", "ENSG00000168394": "TAP1",
        "ENSG00000204264": "PSMB8", "ENSG00000240065": "PSMB9", "ENSG00000120217": "CD274",
        "ENSG00000169245": "CXCL10", "ENSG00000138755": "CXCL9",
    }
    sub = df[df["ens"].isin(want)].copy()
    ut = ["JB-01-A1-TIGK-UT-1_S1", "JB-02-A2-TIGK-UT-2_S2", "JB-03-A3-TIGK-UT-3_S3"]
    b = ["JB-04-B1-TIGK-IFNB-1_S4", "JB-05-B2-TIGK-IFNB-2_S5", "JB-06-B3-TIGK-IFNB-3_S6"]
    l = ["JB-07-C1-TIGK-IFNL-1_S7", "JB-08-C2-TIGK-IFNL-2_S8", "JB-09-C3-TIGK-IFNL-3_S9"]
    # CPM using all genes
    all_cols = ut + b + l
    num_all = df[all_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib = num_all.sum(axis=0)
    sub_num = sub[all_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    logged = np.log2(sub_num.div(lib, axis=1) * 1e6 + 1.0)
    logged.insert(0, "symbol", sub["ens"].map(want).values)
    log_df = logged.groupby("symbol")[all_cols].mean()
    rows, genes = [], []
    for name, treat in (("TIGK IFN-β vs unstimulated", b), ("TIGK IFN-λ vs unstimulated", l)):
        meta = dict(
            accession="GSE184456",
            contrast=name,
            design_class="ifn_stimulus",
            species="human",
            test="Welch on log2(CPM+1)",
            note="Gingival epithelial keratinocytes. Column names encode UT, IFNB, IFNL.",
        )
        row, g = contrast_from_logs(log_df, treat, ut, HUMAN_IFN, "human", meta, True, None)
        rows.append(row)
        genes.extend(g)
    return rows, genes


def load_gse241852():
    df = pd.read_csv(RAW / "GSE241852_RPKM.txt.gz", sep="\t")
    neo_nt = [c for c in df.columns if "NEO_NT" in c]
    neo_pic = [c for c in df.columns if "NEO_PolyIC" in c]
    ko_nt = [c for c in df.columns if "GRHL2KO_NT" in c]
    ko_pic = [c for c in df.columns if "GRHL2KO_PolyIC" in c]
    rows, genes = [], []
    for name, treat, ctrl in (
        ("MCF10A NEO Poly(I:C) vs untreated", neo_pic, neo_nt),
        ("MCF10A GRHL2KO Poly(I:C) vs untreated", ko_pic, ko_nt),
    ):
        log_df = matrix_from_symbol_frame(df, "gene_symbol", treat + ctrl, already_log=False)
        meta = dict(
            accession="GSE241852",
            contrast=name,
            design_class="ifn_stimulus",
            species="human",
            test="Welch on log2(RPKM+1)",
            note="Poly(I:C) is a dsRNA stimulus, not recombinant IFN. RPKM file is linear.",
        )
        row, g = contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, True, None)
        rows.append(row)
        genes.extend(g)
    return rows, genes


def load_gse150255():
    df = pd.read_csv(RAW / "GSE150255.txt.gz", sep="\t")
    rows = []
    genes = []
    pairs = [
        ("A549 IFN-γ 24h vs 0h", "A54924h_FPKM", "A5490h_FPKM"),
        ("HCC827 IFN-γ 24h vs 0h", "HCC82724h_FPKM", "HCC8270h_FPKM"),
    ]
    for name, tcol, ccol in pairs:
        sub = df[["gene_name", tcol, ccol]].copy()
        sub["gene_name"] = sub["gene_name"].astype(str)
        log_df = np.log2(sub.groupby("gene_name")[[tcol, ccol]].mean(numeric_only=True) + 1.0)
        t_sym, c_sym = "TACSTD2", "CLDN4"
        if t_sym not in log_df.index or c_sym not in log_df.index:
            rows.append(empty_row(
                accession="GSE150255",
                contrast=name,
                design_class="unreplicated",
                species="human",
                n_treat=1,
                n_ctrl=1,
                test="not dual readout",
                call="not_evaluable",
                note=f"{t_sym if t_sym not in log_df.index else c_sym} is absent from the deposited gene_name column, so both genes are not read out.",
            ))
            continue
        t_lfc = float(log_df.loc[t_sym, tcol] - log_df.loc[t_sym, ccol])
        c_lfc = float(log_df.loc[c_sym, tcol] - log_df.loc[c_sym, ccol])
        present = [g for g in HUMAN_IFN if g in log_df.index]
        deltas = [float(log_df.loc[g, tcol] - log_df.loc[g, ccol]) for g in present]
        ifn_delta = float(np.mean(deltas)) if deltas else np.nan
        for g, d in zip(present, deltas):
            genes.append({"accession": "GSE150255", "contrast": name, "gene": g, "log2fc": d})
        rows.append(empty_row(
            accession="GSE150255",
            contrast=name,
            design_class="unreplicated",
            species="human",
            n_treat=1,
            n_ctrl=1,
            test="single library per time point; no p-value",
            tacstd2_log2fc=t_lfc,
            cldn4_log2fc=c_lfc,
            ifn_score_delta=ifn_delta,
            ifn_n_genes=len(present),
            ifn_direction=ifn_direction(ifn_delta),
            call="unreplicated",
            note="Series title is an IFN-γ time course. Each time point is one library. 24h versus 0h only. No replicate test.",
        ))
    return rows, genes


def load_gse215771():
    # barcode -> WT group, from family SOFT supplementary filenames
    wt_ctrl = ["AGGCTCATCG-ATCATACCGC", "AGAAGACCTA-GATAGGTTGC", "AACATCTCGA-TATTCGCCAG"]
    wt_ifn = ["TCAACGCGTA-ATGACAGCAC", "CGCAACATGC-ATACGAAGCA", "GTAACGTCAC-ACCAACTAAG"]
    # read header only
    header = pd.read_csv(RAW / "GSE215771_Expression_Values.csv.gz", nrows=0)
    cols = list(header.columns)
    def pick(barcodes):
        chosen = []
        for b in barcodes:
            hits = [c for c in cols if b in c and c.endswith("TPM")]
            if len(hits) != 1:
                raise SystemExit(f"GSE215771 barcode {b} TPM hits {hits}")
            chosen.append(hits[0])
        return chosen
    treat = pick(wt_ifn)
    ctrl = pick(wt_ctrl)
    df = pd.read_csv(RAW / "GSE215771_Expression_Values.csv.gz", usecols=["Name"] + treat + ctrl)
    log_df = matrix_from_symbol_frame(df, "Name", treat + ctrl, already_log=False)
    meta = dict(
        accession="GSE215771",
        contrast="A549 WT IFN-γ vs control",
        design_class="ifn_stimulus",
        species="human",
        test="Welch on log2(TPM+1)",
        note="WT arms only. Barcodes mapped from GEO supplementary filenames in the family SOFT. IRF1KO and NF2KO are not in this contrast.",
    )
    return contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, True, None)


def load_gse306855():
    df = pd.read_excel(RAW / "GSE306855.xlsx", sheet_name="PR3940_gene_expression", usecols=[
        "Gene ID",
        "PR3940_01_a.TPM", "PR3940_02_a.TPM", "PR3940_03_a.TPM",
        "PR3940_04_a.TPM", "PR3940_05_a.TPM", "PR3940_06_a.TPM",
    ])
    df["ens"] = df["Gene ID"].map(strip_ens)
    want = {
        "ENSG00000184292": "TACSTD2", "ENSG00000189143": "CLDN4",
        "ENSG00000187608": "ISG15", "ENSG00000185745": "IFIT1", "ENSG00000119922": "IFIT2",
        "ENSG00000119917": "IFIT3", "ENSG00000157601": "MX1", "ENSG00000183486": "MX2",
        "ENSG00000089127": "OAS1", "ENSG00000111335": "OAS2", "ENSG00000111331": "OAS3",
        "ENSG00000126709": "IFI6", "ENSG00000165949": "IFI27", "ENSG00000137965": "IFI44",
        "ENSG00000137959": "IFI44L", "ENSG00000115415": "STAT1", "ENSG00000125347": "IRF1",
        "ENSG00000185507": "IRF7", "ENSG00000166710": "B2M", "ENSG00000206503": "HLA-A",
        "ENSG00000234745": "HLA-B", "ENSG00000204525": "HLA-C", "ENSG00000168394": "TAP1",
        "ENSG00000204264": "PSMB8", "ENSG00000240065": "PSMB9", "ENSG00000120217": "CD274",
        "ENSG00000169245": "CXCL10", "ENSG00000138755": "CXCL9",
    }
    sub = df[df["ens"].isin(want)].copy()
    tpm_cols = [f"PR3940_0{i}_a.TPM" for i in range(1, 7)]
    sub = sub.groupby("ens")[tpm_cols].mean(numeric_only=True)
    sub.index = sub.index.map(want)
    # submission order is control rep1-3 then IFNγ rep1-3. Confirm with ISG15.
    ctrl = tpm_cols[:3]
    treat = tpm_cols[3:]
    isg = np.log2(sub.loc["ISG15", treat].astype(float) + 1).mean() - np.log2(sub.loc["ISG15", ctrl].astype(float) + 1).mean()
    note = "Columns 01–03 then 04–06 follow GEO submission order (control rep1–3, IFNγ rep1–3)."
    if isg <= 0.5:
        # try the reverse before refusing
        isg_rev = np.log2(sub.loc["ISG15", ctrl].astype(float) + 1).mean() - np.log2(sub.loc["ISG15", treat].astype(float) + 1).mean()
        if isg_rev > 0.5:
            treat, ctrl = ctrl, treat
            note = "Submission-order mapping failed the ISG15 check; columns were swapped so the IFN-γ arm is the ISG15-high arm."
            isg = isg_rev
        else:
            row = empty_row(
                accession="GSE306855",
                contrast="HSC-2 IFN-γ vs PBS",
                design_class="ifn_stimulus",
                species="human",
                call="not_evaluable",
                note=f"ISG15 did not identify the IFN arm (log2FC {isg:.3f}). Contrast refused.",
            )
            return row, []
    note += f" ISG15 log2FC={isg:.3f} on the arm used as IFN-γ."
    log_df = np.log2(sub.astype(float) + 1.0)
    meta = dict(
        accession="GSE306855",
        contrast="HSC-2 IFN-γ vs PBS",
        design_class="ifn_stimulus",
        species="human",
        test="Welch on log2(TPM+1)",
        note=note,
    )
    return contrast_from_logs(log_df, treat, ctrl, HUMAN_IFN, "human", meta, True, None)


def load_gse190899():
    df = pd.read_csv(RAW / "GSE190899_Raw_gene_counts_matrix.txt.gz", sep="\t")
    df["ens"] = df["Geneid"].map(strip_ens)
    want = {
        "ENSG00000184292": "TACSTD2", "ENSG00000189143": "CLDN4",
        "ENSG00000187608": "ISG15", "ENSG00000185745": "IFIT1", "ENSG00000119922": "IFIT2",
        "ENSG00000119917": "IFIT3", "ENSG00000157601": "MX1", "ENSG00000183486": "MX2",
        "ENSG00000089127": "OAS1", "ENSG00000111335": "OAS2", "ENSG00000111331": "OAS3",
        "ENSG00000126709": "IFI6", "ENSG00000165949": "IFI27", "ENSG00000137965": "IFI44",
        "ENSG00000137959": "IFI44L", "ENSG00000115415": "STAT1", "ENSG00000125347": "IRF1",
        "ENSG00000185507": "IRF7", "ENSG00000166710": "B2M", "ENSG00000206503": "HLA-A",
        "ENSG00000234745": "HLA-B", "ENSG00000204525": "HLA-C", "ENSG00000168394": "TAP1",
        "ENSG00000204264": "PSMB8", "ENSG00000240065": "PSMB9", "ENSG00000120217": "CD274",
        "ENSG00000169245": "CXCL10", "ENSG00000138755": "CXCL9",
    }
    sample_cols = [c for c in df.columns if c not in ("Geneid", "Chr", "Start", "End", "Strand", "Length", "ens")]
    num_all = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib = num_all.sum(axis=0)
    sub = df[df["ens"].isin(want)].copy()
    sub_num = sub[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    logged = np.log2(sub_num.div(lib, axis=1) * 1e6 + 1.0)
    logged.insert(0, "symbol", sub["ens"].map(want).values)
    log_df = logged.groupby("symbol")[sample_cols].mean()
    lines = {
        "005-ileum": {
            "N": ["5INA", "5INB", "5INC", "5IND", "5INE", "5INF"],
            "B": ["5IBA", "5IBB", "5IBC"],
            "L": ["5ILA", "5ILB", "5ILC"],
        },
        "005-sigmoid": {
            "N": ["5SNA", "5SNB", "5SNC"],
            "B": ["5SBA", "5SBB", "5SBC"],
            "G": ["5SGA", "5SGB", "5SGC"],
            "L": ["5SLA", "5SLC"],
        },
        "007-sigmoid": {
            "N": ["7SNA", "7SNB", "7SNC"],
            "B": ["7SBA", "7SBB", "7SBC"],
            "G": ["7SGA", "7SGB", "7SGC"],
            "L": ["7SLA", "7SLB", "7SLC"],
        },
    }
    # ISG15 QC of the letter code
    isg_notes = []
    code_ok = True
    for line, arms in lines.items():
        for letter in ("B", "G", "L"):
            if letter not in arms:
                continue
            d = float(log_df.loc["ISG15", arms[letter]].mean() - log_df.loc["ISG15", arms["N"]].mean())
            isg_notes.append(f"{line} {letter} ISG15 log2FC={d:.2f}")
            if d <= 0.5:
                code_ok = False
    label = {"B": "IFN-β", "G": "IFN-γ", "L": "IFN-λ"}
    rows, genes = [], []
    if not code_ok:
        rows.append(empty_row(
            accession="GSE190899",
            contrast="IEC organoid IFN vs untreated",
            design_class="ifn_stimulus",
            species="human",
            call="not_evaluable",
            note="Letter code failed the ISG15 check. " + "; ".join(isg_notes),
        ))
        return rows, genes
    for letter, nice in label.items():
        line_t, line_c, line_ifn = [], [], []
        used = []
        present = [g for g in HUMAN_IFN if g in log_df.index]
        per_gene = {g: [] for g in present}
        for line, arms in lines.items():
            if letter not in arms:
                continue
            tr, ct = arms[letter], arms["N"]
            line_t.append(float(log_df.loc["TACSTD2", tr].mean() - log_df.loc["TACSTD2", ct].mean()))
            line_c.append(float(log_df.loc["CLDN4", tr].mean() - log_df.loc["CLDN4", ct].mean()))
            score, present = score_matrix(log_df, HUMAN_IFN)
            line_ifn.append(float(score[tr].mean() - score[ct].mean()))
            for g in present:
                per_gene[g].append(float(log_df.loc[g, tr].mean() - log_df.loc[g, ct].mean()))
            used.append(line)
        t_lfc, c_lfc = float(np.mean(line_t)), float(np.mean(line_c))
        ifn_delta = float(np.mean(line_ifn))
        t_p, c_p, ifn_p = onesample(line_t), onesample(line_c), onesample(line_ifn)
        call = classify_stimulus(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p, True)
        for g, vals in per_gene.items():
            genes.append({"accession": "GSE190899", "contrast": f"IEC organoids {nice} vs untreated", "gene": g, "log2fc": float(np.mean(vals))})
        # partial correlation on the samples in this contrast
        sample_cols_used = []
        for line, arms in lines.items():
            if letter not in arms:
                continue
            sample_cols_used.extend(arms[letter] + arms["N"])
        row = empty_row(
            accession="GSE190899",
            contrast=f"IEC organoids {nice} vs untreated",
            design_class="ifn_stimulus",
            species="human",
            n_treat=len(line_t),
            n_ctrl=len(line_t),
            test="one-sample t of line-level log2(CPM+1) differences",
            tacstd2_log2fc=t_lfc,
            tacstd2_p=t_p,
            cldn4_log2fc=c_lfc,
            cldn4_p=c_p,
            ifn_score_delta=ifn_delta,
            ifn_score_p=ifn_p,
            ifn_n_genes=len(present),
            ifn_direction=ifn_direction(ifn_delta),
            call=call,
            note=(
                "n is the number of organoid lines, not libraries. "
                f"Lines: {', '.join(used)}. "
                f"TACSTD2 line log2FC={', '.join(f'{x:.3f}' for x in line_t)}. "
                f"CLDN4 line log2FC={', '.join(f'{x:.3f}' for x in line_c)}. "
                "Codes: subject 005 ileum (5I), 005 sigmoid (5S), 007 sigmoid (7S); "
                "B/G/L/N = IFN-β / IFN-γ / IFN-λ / untreated, matching the series summary and confirmed by ISG15. "
                + "; ".join(isg_notes)
            ),
        )
        attach_partial(
            row,
            log_df.loc["TACSTD2", sample_cols_used],
            log_df.loc["CLDN4", sample_cols_used],
            score_matrix(log_df, HUMAN_IFN)[0][sample_cols_used],
        )
        rows.append(row)
    return rows, genes


def load_gse239485():
    df = pd.read_excel(RAW / "GSE239485.xlsx", sheet_name="DataNorm")
    sample_cols = [c for c in df.columns if str(c).startswith(("C_", "D_", "T_"))]
    # detect log: Tacstd2 or a high-abundance gene
    sub = df[["gene_name"] + sample_cols].copy()
    log_candidate = sub.groupby("gene_name")[sample_cols].mean(numeric_only=True)
    p99 = float(np.nanpercentile(log_candidate.to_numpy(), 99))
    already = p99 < 40
    if not already:
        log_df = np.log2(log_candidate.clip(lower=0) + 1.0)
        scale_note = f"99th percentile of the matrix is {p99:.2f}, so values were treated as linear and log2(x+1) was applied."
    else:
        log_df = log_candidate
        scale_note = f"99th percentile of the matrix is {p99:.2f}, so DataNorm was treated as already log2. The reported delta is a difference of those values."
    c_cols = [c for c in sample_cols if str(c).startswith("C_")]
    d_cols = [c for c in sample_cols if str(c).startswith("D_")]
    t_cols = [c for c in sample_cols if str(c).startswith("T_")]
    rows, genes = [], []
    for name, treat, note in (
        ("Poly(I:C)+anti-PD-1 vs vehicle", d_cols, "D_ columns match GSM titles 'Poly I:C and anti-PD-1'. C_ columns are vehicle. Tumor model is not renamed beyond the GSM titles."),
        ("Poly(I:C)+anti-PD-1+anti-C5aR1 vs vehicle", t_cols, "T_ columns match GSM titles 'Poly I:C, anti-PD-1 and anti-C5aR1'."),
    ):
        meta = dict(
            accession="GSE239485",
            contrast=name,
            design_class="immune_treatment",
            species="mouse",
            test="Welch on the normalized matrix",
            note=scale_note + " " + note,
        )
        row, g = contrast_from_logs(log_df, treat, c_cols, MOUSE_IFN, "mouse", meta, False, None)
        rows.append(row)
        genes.extend(g)
    return rows, genes


def load_gse15212():
    # probe id -> symbol for the genes we need
    want = set(HUMAN_IFN + ["TACSTD2", "CLDN4"])
    id_to_symbol = {}
    with gzip.open(RAW / "GPL4133.annot.gz", "rt", errors="replace") as fh:
        started = False
        for ln in fh:
            if ln.startswith("!platform_table_begin"):
                started = True
                header = next(fh).rstrip("\n").split("\t")
                sym_i = header.index("Gene symbol")
                continue
            if not started:
                continue
            if ln.startswith("!platform_table_end"):
                break
            parts = ln.rstrip("\n").split("\t")
            if len(parts) <= sym_i:
                continue
            sym = parts[sym_i]
            if sym in want:
                id_to_symbol.setdefault(sym, []).append(parts[0])
    wanted_ids = {i for ids in id_to_symbol.values() for i in ids}
    # sample titles
    titles = accessions = None
    with gzip.open(RAW / "GSE15212_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_geo_accession"):
                accessions = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!series_matrix_table_begin"):
                break
    if titles is None or accessions is None:
        raise SystemExit("GSE15212 sample header missing")
    treat_gsm = [a for a, t in zip(accessions, titles) if t.startswith("TACSTD2.") and ".72h." in t]
    ctrl_gsm = [a for a, t in zip(accessions, titles) if t.startswith("Neg.2.72h.")]
    # stream values
    values = {i: None for i in wanted_ids}
    with gzip.open(RAW / "GSE15212_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("!series_matrix_table_begin"):
                header = next(fh).rstrip("\n").split("\t")
                header = [h.strip().strip('"') for h in header]
                idx = {g: header.index(g) for g in treat_gsm + ctrl_gsm}
                continue
            if not ln or ln.startswith("!"):
                if ln.startswith("!series_matrix_table_end"):
                    break
                continue
            key = ln.split("\t", 1)[0].strip().strip('"')
            if key not in values:
                continue
            parts = ln.rstrip("\n").split("\t")
            values[key] = {g: float(parts[idx[g]].strip().strip('"')) for g in idx}
    # average probes, values are already log2
    def gene_vector(sym):
        ids = id_to_symbol.get(sym, [])
        mats = []
        for i in ids:
            if values[i] is None:
                continue
            mats.append(pd.Series(values[i]))
        if not mats:
            return None
        return pd.concat(mats, axis=1).mean(axis=1)
    tac = gene_vector("TACSTD2")
    cld = gene_vector("CLDN4")
    if tac is None or cld is None:
        return empty_row(accession="GSE15212", contrast="SW480 TACSTD2 siRNA vs neg, 72h", call="not_evaluable", note="probe missing"), []
    present_genes = []
    gene_mat = {}
    for g in HUMAN_IFN:
        vec = gene_vector(g)
        if vec is not None:
            present_genes.append(g)
            gene_mat[g] = vec
    if len(present_genes) < MIN_IFN:
        return empty_row(accession="GSE15212", call="not_evaluable", note="IFN probes short"), []
    score = pd.DataFrame(gene_mat).mean(axis=1)
    t_lfc = float(tac[treat_gsm].mean() - tac[ctrl_gsm].mean())
    c_lfc = float(cld[treat_gsm].mean() - cld[ctrl_gsm].mean())
    t_p = welch(tac[treat_gsm], tac[ctrl_gsm])
    c_p = welch(cld[treat_gsm], cld[ctrl_gsm])
    ifn_delta = float(score[treat_gsm].mean() - score[ctrl_gsm].mean())
    ifn_p = welch(score[treat_gsm], score[ctrl_gsm])
    gene_rows = []
    for g in present_genes:
        gene_rows.append({
            "accession": "GSE15212",
            "contrast": "SW480 TACSTD2 siRNA vs negative siRNA, 72h",
            "gene": g,
            "log2fc": float(gene_mat[g][treat_gsm].mean() - gene_mat[g][ctrl_gsm].mean()),
        })
    call = classify_genetic(t_lfc, t_p, c_lfc, c_p, ifn_delta, ifn_p)
    row = empty_row(
        accession="GSE15212",
        contrast="SW480 TACSTD2 siRNA vs negative siRNA, 72h",
        design_class="genetic_kd",
        species="human",
        n_treat=len(treat_gsm),
        n_ctrl=len(ctrl_gsm),
        test="Welch on log2 quantile-normalized intensities (already log2)",
        tacstd2_log2fc=t_lfc,
        tacstd2_p=t_p,
        cldn4_log2fc=c_lfc,
        cldn4_p=c_p,
        ifn_score_delta=ifn_delta,
        ifn_score_p=ifn_p,
        ifn_n_genes=len(present_genes),
        ifn_direction=ifn_direction(ifn_delta),
        call=call,
        note=(
            "Both TACSTD2 siRNAs (.1 and .4) at 72h versus Neg.2 72h. "
            f"n_treat={len(treat_gsm)} (titles: {[t for a,t in zip(accessions, titles) if a in treat_gsm]}). "
            "GPL4133 probes averaged when a symbol has more than one. Platform processing: log2 and quantile normalization, so no second log was applied."
        ),
    )
    attach_partial(row, tac[treat_gsm + ctrl_gsm], cld[treat_gsm + ctrl_gsm], score[treat_gsm + ctrl_gsm])
    return row, gene_rows


def main():
    rows = []
    genes = []
    loaders = [
        load_gse334497,
        load_gse289287,
        load_gse245459,
        load_gse207704,
        load_gse50927,
        load_gse304294,
        load_gse311016,
        load_gse274940,
        load_gse15212,
        load_gse156295,
        load_gse184456,
        load_gse241852,
        load_gse150255,
        load_gse215771,
        load_gse306855,
        load_gse190899,
        load_gse239485,
    ]
    for fn in loaders:
        print("running", fn.__name__, flush=True)
        result = fn()
        if isinstance(result, tuple) and result and isinstance(result[0], list):
            rows.extend(result[0])
            genes.extend(result[1])
        else:
            row, g = result
            rows.append(row)
            genes.extend(g)
    table = pd.DataFrame(rows)

    def residual_call(r):
        n = r["partial_n"]
        if n == "" or n is None or (isinstance(n, float) and not np.isfinite(n)) or float(n) < 8:
            return ""
        sig = []
        for rho_key, p_key in (
            ("partial_rho_tacstd2_ifn_given_cldn4", "partial_p_tacstd2"),
            ("partial_rho_cldn4_ifn_given_tacstd2", "partial_p_cldn4"),
        ):
            rho, p = r[rho_key], r[p_key]
            if rho is None or p is None:
                continue
            if not (np.isfinite(rho) and np.isfinite(p)):
                continue
            if p < P_CUT:
                sig.append(float(rho))
        if not sig:
            return "null"
        if all(v > 0 for v in sig):
            return "supports"
        if all(v < 0 for v in sig):
            return "opposite"
        return "opposite"

    table["residual_partial_call"] = table.apply(residual_call, axis=1)
    # literal "null" must survive later reads
    table.to_csv(OUT / "evidence_table.tsv", sep="\t", index=False, na_rep="")
    pd.DataFrame(genes).to_csv(OUT / "ifn_gene_log2fc.tsv", sep="\t", index=False)
    # printable summary
    cols = ["accession", "contrast", "n_treat", "n_ctrl", "tacstd2_log2fc", "tacstd2_p", "cldn4_log2fc", "cldn4_p", "ifn_score_delta", "ifn_score_p", "ifn_direction", "call"]
    with pd.option_context("display.max_rows", 100, "display.max_colwidth", 60, "display.width", 200):
        print(table[cols].to_string(index=False))
    print("wrote", OUT / "evidence_table.tsv")


if __name__ == "__main__":
    main()
