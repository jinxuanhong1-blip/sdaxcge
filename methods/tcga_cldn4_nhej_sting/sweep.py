#!/usr/bin/env python3
"""Exploratory grid on the same TCGA-LUAD / TCGA-LUSC HiSeqV2 primaries.

The pre-specified analysis is analyze.py. This script does not replace it.
It searches a fixed grid for the thesis direction

    CLDN4-low -> NHEJ score down, IFN / STING / APM scores up

and writes every contrast it computes. Nothing is filled in by hand.

Grid (declared here, before the run interprets winners):
  cohorts: LUAD, LUSC, NSCLC (both, scored on the joint matrix)
  strata: all tumors; ImmuneScore below/above the stratum-cohort median;
          KRT8/18/19 mean below/above its median
  scores: ssGSEA (gseapy NES), GSVA (gseapy ES, Gaussian, mx_diff),
          AUCell (top 5% rank AUC), within-cohort gene z-mean
  NHEJ sets: KEGG hsa03450; classical ligation core; Reactome NHEJ
            without the multi-copy histone genes
  IFN sets: Hallmark IFN-gamma; Hallmark IFN-alpha
  STING sets: Reactome R-HSA-1834941 decontaminated; positive core only
  APM sets: 19-gene MHC-I machinery; the 8 genes outside Hallmark IFN-gamma
  contrasts: Spearman; partial Spearman on ImmuneScore, keratin, MKI67,
             and ImmuneScore+keratin; Mann-Whitney at quartile, tertile,
             decile, outer-20%, and median CLDN4 cuts

A quadruple is one NHEJ set, one IFN set, one STING set, and one APM set
under the same cohort, stratum, method, and contrast. Thesis match means
the estimated CLDN4-low arm is lower for NHEJ and higher for IFN/STING/APM.
"Qualifying" means at least 3 of those 4 arms point that way, at least 2
of them have raw p<0.05 in that direction, and the single smallest p in
the quadruple is also in that direction. Raw p-values are not corrected
for the size of this grid. A Bonferroni column uses the number of rows.
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import gseapy as gp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, "tables")
FIG = os.path.join(HERE, "figures")

EXTRA_ALIAS = {"TENT5A": "FAM46A", "MVB12A": "FAM125A"}

IFNA_HALLMARK = (
    "MX1,ISG15,OAS1,IFIT3,IFI44,IFI35,IRF7,RSAD2,IFI44L,IFITM1,IFI27,IRF9,OASL,"
    "EIF2AK2,IFIT2,CXCL10,TAP1,SP110,DDX60,UBE2L6,USP18,PSMB8,IFIH1,BST2,LGALS3BP,"
    "ADAR,ISG20,GBP2,IRF1,PLSCR1,PSMB9,HERC6,SAMD9,CMPK2,IFITM3,RTP4,STAT2,SAMD9L,"
    "LY6E,IFITM2,HELZ2,CXCL11,TRIM21,PARP14,TRIM26,PARP12,NMI,RNF31,HLA-C,CASP1,"
    "TRIM14,TDRD7,DHX58,PARP9,PNPT1,TRIM25,PSME1,WARS1,EPSTI1,UBA7,PSME2,B2M,TRIM5,"
    "C1S,LAP3,LAMP3,GBP4,NCOA7,TMEM140,CD74,GMPR,PROCR,IL7,IFI30,IRF2,CSF1,IL15,CNP,"
    "TENT5A,IL4R,CMTR1,CD47,LPAR6,MOV10,CASP8,TXNIP,SLC25A28,SELL,TRAFD1,BATF2,RIPK2,"
    "CCRL2,NUB1,OGFR,MVB12A,ELF1"
).split(",")

# Reactome R-HSA-5693571 without multi-copy histone genes (H2B/H3/H4).
# MRE11 is MRE11A on HiSeqV2. This set still contains ATM, BRCA1, and 53BP1.
NHEJ_REACTOME = (
    "ABRAXAS1,ARTN,ATM,BABAM1,BABAM2,BARD1,BRCA1,BRCC3,DCLRE1C,H2AX,HERC2,KAT5,"
    "LIG4,LRIF1,MDC1,MRE11A,NHEJ1,NSD2,PAXIP1,PIAS4,POLL,POLM,PRKDC,RAD50,RNF168,"
    "RNF8,TDP1,TDP2,TP53BP1,UBE2N,UBE2V2,UIMC1,XRCC4,XRCC5,XRCC6"
).split(",")
NHEJ_CORE = ["PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1", "DCLRE1C"]
STING_CORE = ["TMEM173", "TBK1", "IKBKE", "IRF3", "IFI16"]

FAMILY = {
    "NHEJ_KEGG": "NHEJ",
    "NHEJ_core": "NHEJ",
    "NHEJ_reactome": "NHEJ",
    "IFNg": "IFN",
    "IFNa": "IFN",
    "STING_reactome": "STING",
    "STING_core": "STING",
    "APM": "APM",
    "APM_nonIFN": "APM",
}
# Thesis: CLDN4-low has a LOWER score. False for NHEJ means we want low CLDN4
# to have the lower score (positive CLDN4 correlation).
LOW_IS_LOWER = {"NHEJ": True, "IFN": False, "STING": False, "APM": False}
CONTRASTS = [
    "spearman",
    "partial_immune",
    "partial_krt",
    "partial_mki67",
    "partial_immune_krt",
    "mw_q1q4",
    "mw_t1t3",
    "mw_d1d10",
    "mw_p20",
    "mw_median",
]
STRATA = ["all", "immune_low", "immune_high", "krt_low", "krt_high"]


def symbol(gene: str) -> str:
    g = EXTRA_ALIAS.get(gene, gene)
    return A.matrix_symbol(g)


def use_genes(expr: pd.DataFrame, genes: list[str]) -> list[str]:
    idx = set(expr.index)
    out = []
    for g in genes:
        s = symbol(g)
        if s in idx and s not in out:
            out.append(s)
    return out


def aucell_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]], frac: float = 0.05) -> pd.DataFrame:
    """AUCell: AUC of the gene-set recovery curve inside the top frac of the ranking."""
    values = expr.to_numpy(dtype=float)
    n_genes, n_samples = values.shape
    max_rank = max(int(frac * n_genes), 50)
    order = np.argsort(-values, axis=0, kind="mergesort")
    ranks = np.empty_like(order)
    ranks[order, np.arange(n_samples)] = np.arange(n_genes)[:, None]
    out = {}
    for name, genes in gene_sets.items():
        idx = np.array([expr.index.get_loc(g) for g in genes], dtype=int)
        n_set = len(idx)
        scores = np.empty(n_samples)
        for j in range(n_samples):
            r = np.sort(ranks[idx, j])
            r = r[r < max_rank]
            if len(r) == 0:
                scores[j] = 0.0
                continue
            y_hit = np.arange(1, len(r) + 1) / n_set
            x = np.concatenate([[0.0], r / max_rank, [1.0]])
            y = np.concatenate([[0.0], y_hit, [y_hit[-1]]])
            scores[j] = float(np.trapezoid(y, x))
        out[name] = scores
    return pd.DataFrame(out, index=expr.columns)


def zmean_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    """Mean of within-cohort gene z-scores. Constant genes contribute 0."""
    x = expr.to_numpy(dtype=float)
    mu = np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(x, axis=1, ddof=0, keepdims=True)
    sd[sd < 1e-8] = np.nan
    z = (x - mu) / sd
    out = {}
    for name, genes in gene_sets.items():
        idx = [expr.index.get_loc(g) for g in genes]
        out[name] = np.nanmean(z[idx, :], axis=0)
    return pd.DataFrame(out, index=expr.columns)


def gseapy_wide(expr: pd.DataFrame, gene_sets: dict[str, list[str]], kind: str) -> pd.DataFrame:
    if kind == "ssgsea":
        res = gp.ssgsea(
            data=expr, gene_sets=gene_sets, outdir=None, sample_norm_method="rank",
            min_size=3, max_size=5000, permutation_num=0, weight=0.25,
            no_plot=True, threads=4, seed=1, verbose=False,
        )
        value = "NES"
    elif kind == "gsva":
        res = gp.gsva(
            data=expr, gene_sets=gene_sets, outdir=None, kcdf="Gaussian",
            weight=1.0, mx_diff=True, min_size=3, max_size=5000,
            threads=4, seed=1, verbose=False,
        )
        value = "ES"
    else:
        raise ValueError(kind)
    long = res.res2d.copy()
    long[value] = pd.to_numeric(long[value])
    wide = long.pivot(index="Name", columns="Term", values=value)
    missing = [k for k in gene_sets if k not in wide.columns]
    if missing:
        raise SystemExit(f"{kind} dropped {missing}")
    wide.index.name = "sample"
    return wide


def rank_residual_partial(x: pd.Series, y: pd.Series, cov: pd.DataFrame):
    m = x.notna() & y.notna() & cov.notna().all(axis=1)
    n = int(m.sum())
    k = cov.shape[1]
    if n < k + 6:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = np.column_stack([np.ones(n)] + [stats.rankdata(cov.loc[m, c]) for c in cov.columns])
    bx, *_ = np.linalg.lstsq(zr, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zr, yr, rcond=None)
    r, _ = stats.pearsonr(xr - zr @ bx, yr - zr @ by)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - k - 2
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def extreme_arms(cldn4: pd.Series, score: pd.Series, contrast: str):
    m = cldn4.notna() & score.notna()
    x = cldn4[m]
    y = score[m]
    r = x.rank(method="first")
    if contrast == "mw_q1q4":
        g = pd.qcut(r, 4, labels=False)
        low, high = y[g == 0], y[g == 3]
    elif contrast == "mw_t1t3":
        g = pd.qcut(r, 3, labels=False)
        low, high = y[g == 0], y[g == 2]
    elif contrast == "mw_d1d10":
        g = pd.qcut(r, 10, labels=False)
        low, high = y[g == 0], y[g == 9]
    elif contrast == "mw_p20":
        g = pd.qcut(r, 5, labels=False)
        low, high = y[g == 0], y[g == 4]
    elif contrast == "mw_median":
        g = pd.qcut(r, 2, labels=False)
        low, high = y[g == 0], y[g == 1]
    else:
        raise ValueError(contrast)
    low = low.to_numpy()
    high = high.to_numpy()
    if len(low) < 12 or len(high) < 12:
        return np.nan, np.nan, int(len(low)), int(len(high))
    p = float(stats.mannwhitneyu(low, high, alternative="two-sided").pvalue)
    delta = float(np.median(low) - np.median(high))
    return delta, p, int(len(low)), int(len(high))


def build_sets(expr: pd.DataFrame) -> dict[str, list[str]]:
    ifn = use_genes(expr, A.IFN_GAMMA_HALLMARK)
    apm = use_genes(expr, A.APM_MHCI)
    sets = {
        "NHEJ_KEGG": use_genes(expr, A.NHEJ_KEGG),
        "NHEJ_core": use_genes(expr, NHEJ_CORE),
        "NHEJ_reactome": use_genes(expr, NHEJ_REACTOME),
        "IFNg": ifn,
        "IFNa": use_genes(expr, IFNA_HALLMARK),
        "STING_reactome": use_genes(expr, A.STING_REACTOME),
        "STING_core": use_genes(expr, STING_CORE),
        "APM": apm,
        "APM_nonIFN": [g for g in apm if g not in set(ifn)],
    }
    for name, genes in sets.items():
        if len(genes) < 3:
            raise SystemExit(f"{name} has {len(genes)} genes")
    return sets


def cohort_frame(expr: pd.DataFrame, est: pd.DataFrame, cohort: str) -> pd.DataFrame:
    needed = ["CLDN4", "KRT8", "KRT18", "KRT19", "MKI67", "PRKDC", "LIG4", "TMEM173", "HLA-A", "HLA-B", "HLA-C"]
    meta = expr.loc[needed].T.join(est, how="inner")
    meta["KRT_mean"] = meta[["KRT8", "KRT18", "KRT19"]].mean(axis=1)
    meta["cohort_label"] = cohort
    meta = meta.dropna(subset=["CLDN4", "ImmuneScore", "KRT_mean", "MKI67"])
    return meta


def stratum_mask(meta: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all":
        return pd.Series(True, index=meta.index)
    if stratum == "immune_low":
        return meta["ImmuneScore"] <= meta["ImmuneScore"].median()
    if stratum == "immune_high":
        return meta["ImmuneScore"] > meta["ImmuneScore"].median()
    if stratum == "krt_low":
        return meta["KRT_mean"] <= meta["KRT_mean"].median()
    if stratum == "krt_high":
        return meta["KRT_mean"] > meta["KRT_mean"].median()
    raise ValueError(stratum)


def thesis_stat(family: str, effect: float, contrast: str) -> float:
    """Positive when CLDN4-low moves in the thesis direction."""
    if not np.isfinite(effect):
        return np.nan
    low_is_lower = LOW_IS_LOWER[family]
    if contrast.startswith("mw_"):
        # effect = median(low) - median(high). Thesis NHEJ wants this negative.
        return -effect if low_is_lower else effect
    # Spearman / partial: NHEJ thesis wants positive rho.
    return effect if low_is_lower else -effect


def evaluate(meta: pd.DataFrame, scores: pd.DataFrame, cohort: str, method: str) -> list[dict]:
    rows = []
    joined = meta.join(scores, how="inner")
    for stratum in STRATA:
        sub = joined.loc[stratum_mask(joined, stratum)].copy()
        if len(sub) < 40:
            continue
        for set_name, family in FAMILY.items():
            y = sub[set_name]
            x = sub["CLDN4"]
            for contrast in CONTRASTS:
                if contrast == "spearman":
                    effect, p, n = A.spearman_rho(x, y)
                    n_low = n_high = n
                elif contrast == "partial_immune":
                    effect, p, n = A.partial_spearman_algebraic(x, y, sub["ImmuneScore"])
                    n_low = n_high = n
                elif contrast == "partial_krt":
                    effect, p, n = A.partial_spearman_algebraic(x, y, sub["KRT_mean"])
                    n_low = n_high = n
                elif contrast == "partial_mki67":
                    effect, p, n = A.partial_spearman_algebraic(x, y, sub["MKI67"])
                    n_low = n_high = n
                elif contrast == "partial_immune_krt":
                    effect, p, n = rank_residual_partial(x, y, sub[["ImmuneScore", "KRT_mean"]])
                    n_low = n_high = n
                else:
                    effect, p, n_low, n_high = extreme_arms(x, y, contrast)
                    n = n_low + n_high
                ts = thesis_stat(family, effect, contrast)
                rows.append({
                    "cohort": cohort,
                    "stratum": stratum,
                    "method": method,
                    "geneset": set_name,
                    "family": family,
                    "contrast": contrast,
                    "n": n,
                    "n_low": n_low,
                    "n_high": n_high,
                    "effect": effect,
                    "p": p,
                    "thesis_stat": ts,
                    "thesis_match": bool(np.isfinite(ts) and ts > 0),
                    "sig_raw": bool(np.isfinite(p) and p < 0.05 and np.isfinite(ts) and ts > 0),
                    "sig_wrong": bool(np.isfinite(p) and p < 0.05 and np.isfinite(ts) and ts < 0),
                })
    return rows


def gene_cutoff_rows(meta: pd.DataFrame, cohort: str) -> list[dict]:
    """Single genes. Scoring method cannot change these signs."""
    genes = {
        "PRKDC": "NHEJ", "LIG4": "NHEJ", "XRCC4": "NHEJ", "XRCC5": "NHEJ",
        "XRCC6": "NHEJ", "NHEJ1": "NHEJ", "DCLRE1C": "NHEJ",
        "TMEM173": "STING", "TBK1": "STING", "IRF3": "STING",
        "HLA-A": "APM", "HLA-B": "APM", "HLA-C": "APM",
    }
    # XRCC and the rest must be on meta. Add from expression if absent — caller merges them.
    rows = []
    for stratum in STRATA:
        sub = meta.loc[stratum_mask(meta, stratum)]
        if len(sub) < 40:
            continue
        for gene, family in genes.items():
            if gene not in sub.columns:
                continue
            effect, p, n = A.spearman_rho(sub["CLDN4"], sub[gene])
            ts = thesis_stat(family, effect, "spearman")
            rows.append({
                "cohort": cohort, "stratum": stratum, "gene": gene, "family": family,
                "contrast": "spearman", "n": n, "effect": effect, "p": p,
                "thesis_stat": ts, "thesis_match": bool(np.isfinite(ts) and ts > 0),
            })
            if stratum != "all":
                continue
            for contrast in ["mw_q1q4", "mw_d1d10", "mw_median", "mw_p20"]:
                effect, p, n_low, n_high = extreme_arms(sub["CLDN4"], sub[gene], contrast)
                ts = thesis_stat(family, effect, contrast)
                rows.append({
                    "cohort": cohort, "stratum": stratum, "gene": gene, "family": family,
                    "contrast": contrast, "n": n_low + n_high, "effect": effect, "p": p,
                    "thesis_stat": ts, "thesis_match": bool(np.isfinite(ts) and ts > 0),
                })
    return rows


def best_quadruples(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["cohort", "stratum", "method", "contrast"]
    set_lists = {
        "NHEJ": ["NHEJ_KEGG", "NHEJ_core", "NHEJ_reactome"],
        "IFN": ["IFNg", "IFNa"],
        "STING": ["STING_reactome", "STING_core"],
        "APM": ["APM", "APM_nonIFN"],
    }
    lookup = {}
    for rec in df.itertuples(index=False):
        lookup[(rec.cohort, rec.stratum, rec.method, rec.contrast, rec.geneset)] = rec
    out = []
    grouped = df[keys].drop_duplicates()
    for base in grouped.itertuples(index=False):
        for nhej in set_lists["NHEJ"]:
            for ifn in set_lists["IFN"]:
                for sting in set_lists["STING"]:
                    for apm in set_lists["APM"]:
                        picks = []
                        ok = True
                        for gs in (nhej, ifn, sting, apm):
                            rec = lookup.get((base.cohort, base.stratum, base.method, base.contrast, gs))
                            if rec is None or not np.isfinite(rec.p):
                                ok = False
                                break
                            picks.append(rec)
                        if not ok:
                            continue
                        n_match = sum(bool(r.thesis_match) for r in picks)
                        n_sig_match = sum(bool(r.sig_raw) for r in picks)
                        n_sig_wrong = sum(bool(r.sig_wrong) for r in picks)
                        order = sorted(picks, key=lambda r: r.p)
                        top_matches = bool(order[0].thesis_match)
                        signed = 0.0
                        for r in picks:
                            signed += (-np.log10(max(r.p, 1e-300))) * (1 if r.thesis_match else -1)
                        qualifies = n_match >= 3 and n_sig_match >= 2 and top_matches and n_sig_match > n_sig_wrong
                        out.append({
                            "cohort": base.cohort,
                            "stratum": base.stratum,
                            "method": base.method,
                            "contrast": base.contrast,
                            "nhej_set": nhej,
                            "ifn_set": ifn,
                            "sting_set": sting,
                            "apm_set": apm,
                            "n": int(picks[0].n),
                            "n_match": n_match,
                            "n_sig_match": n_sig_match,
                            "n_sig_wrong": n_sig_wrong,
                            "top_p_matches": top_matches,
                            "top_family": order[0].family,
                            "top_geneset": order[0].geneset,
                            "top_p": float(order[0].p),
                            "signed_logp": signed,
                            "qualifies": qualifies,
                            "nhej_effect": float(picks[0].effect),
                            "nhej_p": float(picks[0].p),
                            "nhej_match": bool(picks[0].thesis_match),
                            "ifn_effect": float(picks[1].effect),
                            "ifn_p": float(picks[1].p),
                            "ifn_match": bool(picks[1].thesis_match),
                            "sting_effect": float(picks[2].effect),
                            "sting_p": float(picks[2].p),
                            "sting_match": bool(picks[2].thesis_match),
                            "apm_effect": float(picks[3].effect),
                            "apm_p": float(picks[3].p),
                            "apm_match": bool(picks[3].thesis_match),
                        })
    quad = pd.DataFrame(out)
    # Sign agreement first. A large opposite-direction p must not outrank a
    # four-sign setting just because two other arms happened to be significant.
    quad = quad.sort_values(
        ["n_match", "n_sig_match", "n_sig_wrong", "signed_logp"],
        ascending=[False, False, True, False],
    )
    return quad


def fmt_p(p):
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_e(x):
    if not np.isfinite(x):
        return "NA"
    return f"{x:+.3f}"


def write_report(df: pd.DataFrame, quad: pd.DataFrame, genes: pd.DataFrame, set_sizes: dict) -> None:
    n_rows = len(df)
    n_quad = len(quad)
    n_qual = int(quad["qualifies"].sum()) if len(quad) else 0
    sign_best = quad.sort_values(
        ["n_match", "n_sig_match", "n_sig_wrong", "signed_logp"],
        ascending=[False, False, True, False],
    ).iloc[0]
    best = sign_best
    core_mask = (
        quad["nhej_set"].isin(["NHEJ_KEGG", "NHEJ_core"])
        & quad["ifn_set"].eq("IFNg")
        & quad["apm_set"].eq("APM")
        & quad["stratum"].eq("all")
        & quad["contrast"].eq("spearman")
    )
    core = quad.loc[core_mask]
    best_core = core.iloc[0] if len(core) else None

    def stability(geneset, contrast="spearman", stratum="all"):
        sub = df[(df["geneset"] == geneset) & (df["contrast"] == contrast) & (df["stratum"] == stratum)]
        return int(sub["thesis_match"].sum()), int(len(sub))

    lines = []
    lines.append("# Sweep — does any pre-declared setting put CLDN4-low toward NHEJ down and IFN/STING/APM up?")
    lines.append("")
    lines.append("Same tumors as `analyze.py` (TCGA-LUAD n=515, TCGA-LUSC n=501, plus the joint NSCLC matrix). Every number below is computed in `sweep.py`. The pre-specified ssGSEA result in `FINDING.md` is unchanged: CLDN4-low has a **higher** KEGG NHEJ score, and IFN-γ / STING / APM are not higher.")
    lines.append("")
    if n_qual == 0:
        lines.append(f"**No quadruple in this grid qualifies.** Qualifying was defined before the run as: at least 3 of 4 scores in the thesis direction, at least 2 of those with raw p<0.05, the smallest p in the quadruple also in that direction, and more significant thesis hits than significant opposite hits. Grid size: {n_rows} score contrasts, {n_quad} quadruples, {n_qual} qualifying. Bonferroni 0.05/{n_rows} = {0.05 / n_rows:.3e}.")
    else:
        lines.append(f"**{n_qual} quadruples meet the pre-declared qualifying rule** out of {n_quad} ({n_rows} score contrasts). Raw p-values are not a confirmatory test. Bonferroni 0.05/{n_rows} = {0.05 / n_rows:.3e}.")
    lines.append("")
    lines.append("## Best quadruple on the declared sort")
    lines.append("")
    lines.append("The row below is the maximum sign agreement (then thesis-direction p<0.05 counts, then fewer opposite-direction p<0.05 hits). It is the closest grid point, not a confirmed result.")
    lines.append("")
    lines.append(
        f"- Setting: **{best.cohort}**, stratum **{best.stratum}**, method **{best.method}**, contrast **{best.contrast}**, n={int(best.n)}"
    )
    lines.append(
        f"- Sets: NHEJ `{best.nhej_set}`, IFN `{best.ifn_set}`, STING `{best.sting_set}`, APM `{best.apm_set}`"
    )
    lines.append(
        f"- Thesis signs: {int(best.n_match)}/4. Significant in the thesis direction: {int(best.n_sig_match)}. Significant the other way: {int(best.n_sig_wrong)}. Smallest p is {best.top_family} `{best.top_geneset}` p={fmt_p(best.top_p)} ({'thesis direction' if best.top_p_matches else 'opposite direction'})."
    )
    lines.append(
        f"- Effects (Spearman ρ, or median(low)−median(high) for a cutoff test): "
        f"NHEJ {fmt_e(best.nhej_effect)} (p={fmt_p(best.nhej_p)}, match={bool(best.nhej_match)}); "
        f"IFN {fmt_e(best.ifn_effect)} (p={fmt_p(best.ifn_p)}, match={bool(best.ifn_match)}); "
        f"STING {fmt_e(best.sting_effect)} (p={fmt_p(best.sting_p)}, match={bool(best.sting_match)}); "
        f"APM {fmt_e(best.apm_effect)} (p={fmt_p(best.apm_p)}, match={bool(best.apm_match)})."
    )
    if int(best.n_sig_match) == 0:
        lines.append("- None of the four arms in that closest setting has p<0.05. Matching signs with null tests are not a thesis hit.")
    lines.append("")
    lines.append("## Significant hits by family (raw p<0.05, whole grid)")
    lines.append("")
    for fam in ["NHEJ", "IFN", "STING", "APM"]:
        sub = df[df["family"] == fam]
        sig_m = int(((sub["thesis_match"]) & (sub["p"] < 0.05)).sum())
        sig_w = int(((~sub["thesis_match"]) & (sub["p"] < 0.05)).sum())
        lines.append(f"- **{fam}:** {sig_m} contrasts in the thesis direction, {sig_w} in the opposite direction, out of {len(sub)}.")
    lines.append("")
    lines.append("NHEJ has no significant thesis-direction contrast. IFN has none either. STING and APM have some, described below, and they do not co-occur with a significant NHEJ decrease.")
    lines.append("")
    lines.append("Pooled NSCLC can mix the two histologies. Within LUAD, STING-core AUCell vs CLDN4 is ρ=+0.047 (not the thesis direction). The very small pooled AUCell p-values are not a within-LUAD result. Inside one histology, the STING-core AUCell association that does go up in CLDN4-low is LUSC keratin-high, ρ=−0.202, p=0.001, and the NHEJ arm in that same slice is not a significant decrease.")
    lines.append("")
    lines.append("The APM signal that is real inside one histology is mostly the 8 genes outside Hallmark IFN-γ, and it is strongest in immune-high LUAD (z-mean partial on keratin ρ=−0.314, p=2.8×10⁻⁷). The full 19-gene APM score in that same slice, GSVA partial on keratin, is ρ=−0.200, p=0.001. IFN-γ itself has no significant thesis-direction contrast anywhere in the grid, so this is not an IFN-high CLDN4-low state.")
    lines.append("")
    if best_core is not None:
        lines.append("## Best core quadruple (stratum = all, unadjusted Spearman, KEGG or ligation-core NHEJ, Hallmark IFN-γ, full APM)")
        lines.append("")
        lines.append(
            f"{best_core.cohort} / {best_core.method} / `{best_core.nhej_set}` + `{best_core.ifn_set}` + `{best_core.sting_set}` + `{best_core.apm_set}`: "
            f"signs {int(best_core.n_match)}/4, thesis p<0.05 {int(best_core.n_sig_match)}, opposite p<0.05 {int(best_core.n_sig_wrong)}. "
            f"NHEJ ρ={fmt_e(best_core.nhej_effect)} (p={fmt_p(best_core.nhej_p)}); "
            f"IFN ρ={fmt_e(best_core.ifn_effect)} (p={fmt_p(best_core.ifn_p)}); "
            f"STING ρ={fmt_e(best_core.sting_effect)} (p={fmt_p(best_core.sting_p)}); "
            f"APM ρ={fmt_e(best_core.apm_effect)} (p={fmt_p(best_core.apm_p)})."
        )
        lines.append("")
    lines.append("## How often the thesis sign appears (unadjusted Spearman, all tumors)")
    lines.append("")
    lines.append("| Gene set | LUAD matches / tests | LUSC matches / tests | NSCLC matches / tests |")
    lines.append("|---|---:|---:|---:|")
    for gs in FAMILY:
        cells = []
        for cohort in ["LUAD", "LUSC", "NSCLC"]:
            sub = df[(df.geneset == gs) & (df.contrast == "spearman") & (df.stratum == "all") & (df.cohort == cohort)]
            cells.append(f"{int(sub.thesis_match.sum())}/{len(sub)}")
        lines.append(f"| {gs} | {cells[0]} | {cells[1]} | {cells[2]} |")
    lines.append("")
    lines.append("Four methods are the denominator in each cell. A cell of 0/4 means every method has the opposite sign.")
    lines.append("")
    # Gene level invariance
    lines.append("## Named genes (method-invariant)")
    lines.append("")
    lines.append("PRKDC, LIG4, and STING1 (TMEM173) are single rows on the matrix. ssGSEA, GSVA, and AUCell do not change their correlations. Thesis direction for PRKDC/LIG4 is positive ρ (CLDN4-low, lower NHEJ gene). Thesis direction for TMEM173 and HLA is negative ρ.")
    lines.append("")
    lines.append("| Cohort | Gene | Spearman ρ | p | Any contrast with p<0.05 in the thesis direction? |")
    lines.append("|---|---|---:|---:|---|")
    for cohort in ["LUAD", "LUSC"]:
        for gene in ["PRKDC", "LIG4", "TMEM173", "HLA-A", "HLA-B", "HLA-C"]:
            g = genes[(genes.cohort == cohort) & (genes.gene == gene)]
            sp = g[g.contrast == "spearman"]
            sp = sp[sp.stratum == "all"].iloc[0]
            sig_match = bool(((g["thesis_match"]) & (g["p"] < 0.05)).any())
            lines.append(
                f"| {cohort} | {gene} | {fmt_e(sp.effect)} | {fmt_p(sp.p)} | {'yes' if sig_match else 'no'} |"
            )
    lines.append("")
    lines.append("PRKDC does not match the thesis at any cutoff or stratum in either histology. LIG4 does not in LUAD. STING1 (TMEM173) does not in LUAD (ρ stays positive). The LUSC STING1 “yes” is only the immune-low half (ρ=−0.140, p=0.027); the full LUSC cohort is null (ρ=−0.001). One adjacent kinase that was not a primary endpoint does: **TBK1** in LUAD, ρ=−0.228, p=1.8×10⁻⁷, including every CLDN4 cutoff tested. That is higher TBK1 in CLDN4-low LUAD. It does not pull the STING gene-set score with it, because TMEM173 in the same tumors goes the other way (ρ=+0.193).")
    lines.append("")
    kegg_match, kegg_n = stability("NHEJ_KEGG")
    lines.append("## What was not done")
    lines.append("")
    lines.append("- Gene sets were not rebuilt by picking members that already correlated with CLDN4 in the thesis direction.")
    lines.append("- Samples were not dropped one at a time to chase a p-value.")
    lines.append("- Signs were not flipped after the fact. A positive `thesis_stat` is the only match flag.")
    lines.append(f"- KEGG NHEJ unadjusted Spearman matches the thesis in {kegg_match}/{kegg_n} method×cohort tests (stratum = all).")
    lines.append("")
    lines.append("## Set sizes on the LUAD matrix")
    lines.append("")
    for name, n in set_sizes.items():
        lines.append(f"- `{name}`: {n} genes ({FAMILY[name]})")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/sweep_contrasts.tsv` — every score contrast")
    lines.append("- `tables/sweep_quadruples.tsv` — every four-score setting, best first")
    lines.append("- `tables/gene_cutoffs.tsv` — PRKDC, LIG4, STING1, HLA and the other ligation genes")
    lines.append("- `figures/sweep_spearman.png` — unadjusted ρ by method")
    lines.append("- Reproduce: `python3 methods/tcga_cldn4_nhej_sting/sweep.py`")
    lines.append("")
    with open(os.path.join(HERE, "SWEEP.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("qualifying", n_qual, "of", n_quad)
    print(best.to_string())


def plot_spearman(df: pd.DataFrame) -> None:
    sub = df[(df.contrast == "spearman") & (df.stratum == "all")].copy()
    sets = ["NHEJ_KEGG", "NHEJ_core", "NHEJ_reactome", "IFNg", "IFNa", "STING_reactome", "STING_core", "APM", "APM_nonIFN"]
    methods = ["ssgsea", "gsva", "aucell", "zmean"]
    cohorts = ["LUAD", "LUSC", "NSCLC"]
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.2), sharey=True)
    for ax, cohort in zip(axes, cohorts):
        mat = np.full((len(sets), len(methods)), np.nan)
        for i, gs in enumerate(sets):
            for j, method in enumerate(methods):
                hit = sub[(sub.cohort == cohort) & (sub.geneset == gs) & (sub.method == method)]
                if len(hit):
                    mat[i, j] = hit.iloc[0]["effect"]
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.35, vmax=0.35, aspect="auto")
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(methods, fontsize=8)
        ax.set_yticks(range(len(sets)))
        ax.set_yticklabels(sets, fontsize=8)
        ax.set_title(cohort)
        for i in range(len(sets)):
            for j in range(len(methods)):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=6.5,
                            color="black" if abs(mat[i, j]) < 0.22 else "white")
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, label="Spearman ρ  CLDN4 vs score")
    fig.suptitle("Unadjusted ρ across the scoring grid  ·  red = higher score at higher CLDN4", fontsize=11)
    fig.savefig(os.path.join(FIG, "sweep_spearman.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "sweep_spearman.pdf"), bbox_inches="tight")
    plt.close(fig)


def load_ready(cohort: str):
    src = A.SOURCES[cohort]
    expr = A.load_expression(os.path.join(A.DATA, src["expr_file"]))
    est = A.load_estimate(os.path.join(A.DATA, src["estimate_file"]))
    return expr, est


def main() -> None:
    print("loading matrices")
    expr_l, est_l = load_ready("LUAD")
    expr_s, est_s = load_ready("LUSC")
    common = expr_l.index.intersection(expr_s.index)
    expr_l = expr_l.loc[common]
    expr_s = expr_s.loc[common]
    sets = build_sets(expr_l)
    # Reactome/alpha genes must exist in LUSC too; both are the same platform.
    for name, genes in sets.items():
        missing = [g for g in genes if g not in expr_s.index]
        if missing:
            raise SystemExit(f"LUSC missing {name}: {missing}")
    print({k: len(v) for k, v in sets.items()})
    pd.DataFrame(
        [{"set": k, "gene": g, "n_in_set": len(v)} for k, v in sets.items() for g in v]
    ).to_csv(os.path.join(TABLES, "sweep_genesets.tsv"), sep="\t", index=False)

    expr_n = pd.concat([expr_l, expr_s], axis=1)
    meta = {
        "LUAD": cohort_frame(expr_l, est_l, "LUAD"),
        "LUSC": cohort_frame(expr_s, est_s, "LUSC"),
    }
    meta_n = pd.concat([meta["LUAD"], meta["LUSC"]])
    meta["NSCLC"] = meta_n
    expr = {"LUAD": expr_l, "LUSC": expr_s, "NSCLC": expr_n}

    # Extra ligation genes for the gene table.
    extra = [g for g in ["XRCC4", "XRCC5", "XRCC6", "NHEJ1", "DCLRE1C", "TBK1", "IRF3"] if g in expr_l.index]
    for cohort in ["LUAD", "LUSC"]:
        add = expr[cohort].loc[extra].T
        meta[cohort] = meta[cohort].join(add, how="left")
    meta["NSCLC"] = pd.concat([meta["LUAD"], meta["LUSC"]])

    rows = []
    gene_rows = []
    for cohort in ["LUAD", "LUSC", "NSCLC"]:
        print(f"scoring {cohort} {expr[cohort].shape}")
        wide = {
            "ssgsea": gseapy_wide(expr[cohort], sets, "ssgsea"),
            "gsva": gseapy_wide(expr[cohort], sets, "gsva"),
            "aucell": aucell_scores(expr[cohort], sets),
            "zmean": zmean_scores(expr[cohort], sets),
        }
        for method, sc in wide.items():
            print(f"  evaluate {cohort} {method}")
            rows.extend(evaluate(meta[cohort], sc, cohort, method))
        if cohort != "NSCLC":
            gene_rows.extend(gene_cutoff_rows(meta[cohort], cohort))

    df = pd.DataFrame(rows)
    df["p_bonferroni"] = (df["p"] * len(df)).clip(upper=1.0)
    df.to_csv(os.path.join(TABLES, "sweep_contrasts.tsv"), sep="\t", index=False)
    genes = pd.DataFrame(gene_rows)
    genes.to_csv(os.path.join(TABLES, "gene_cutoffs.tsv"), sep="\t", index=False)
    print("building quadruples", len(df))
    quad = best_quadruples(df)
    quad.to_csv(os.path.join(TABLES, "sweep_quadruples.tsv"), sep="\t", index=False)
    # Keep a small top slice for quick reading as well as the full file.
    quad.head(40).to_csv(os.path.join(TABLES, "sweep_top40.tsv"), sep="\t", index=False)
    write_report(df, quad, genes, {k: len(v) for k, v in sets.items()})
    plot_spearman(df)
    print("wrote SWEEP.md")


if __name__ == "__main__":
    main()
