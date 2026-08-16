#!/usr/bin/env python3
"""
Purity-adjusted TACSTD2 (TROP2) vs immune-signature correlations in open
durvalumab NSCLC trial RNA-seq (GEO).

Datasets (verified accessions, publication-level FPKM matrices from GEO suppl):
  - GSE253564: PRE-treatment tumor RNA-seq (n=32), randomized phase II
    neoadjuvant durvalumab +/- SBRT in stages I-III NSCLC (PMID 38401548).
  - GSE248378: POST-treatment resected tumors (n=29), same trial
    (PMIDs 38114518, 38401548).

Purity: ESTIMATE (Yoshihara 2013) — official 141-gene Stromal/Immune sets from
the `estimate` R package v1.0.13 (R-Forge), ssGSEA (alpha=0.25) on within-sample
gene ranks, ESTIMATEScore = Stromal + Immune, purity = cos(0.6049872018 +
0.0001467884 * ESTIMATEScore). NOTE: the purity formula was calibrated on
Affymetrix arrays; applied to RNA-seq FPKM it is an approximation (common
practice, but flagged honestly).

Immune scores:
  - ESTIMATE ImmuneScore (ssGSEA)
  - Ayers 6-gene IFN-gamma signature (mean log2 FPKM)
  - Ayers 18-gene T-cell-inflamed GEP (unweighted mean log2 FPKM; the original
    is a weighted, housekeeping-normalized NanoString score — approximation)
  - Rooney cytolytic activity CYT (geometric mean of GZMA, PRF1)
  - CD8A single gene

Statistics: Spearman rho (unadjusted) and purity-adjusted partial Spearman
(rank-transform, residualize on purity ranks, Pearson on residuals; p from
t-distribution with n-3 df). Small n — treat as hypothesis-level evidence only.

Caveat: ESTIMATE purity is derived from immune+stromal scores, so adjusting
immune-vs-TACSTD2 correlations for ESTIMATE purity is partially circular
(especially for the ImmuneScore row itself). This mirrors TIMER-style partial
correlation practice; rows using non-ESTIMATE immune scores are the primary
readout.
"""

import gzip
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "hunt_pacific")
OUT = os.path.join(BASE, "results", "hunt_pacific")
os.makedirs(OUT, exist_ok=True)

GMT = os.path.join(DATA, "SI_geneset.gmt")

IFNG6 = ["IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA"]
GEP18 = ["CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
         "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
         "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"]
CYT = ["GZMA", "PRF1"]


def load_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_fpkm(path, drop_cols=()):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "gene"})
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    df = df.set_index("gene")
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    # collapse duplicate gene symbols by max expression
    df = df.groupby(level=0).max()
    return df


def ssgsea(expr_log, gene_set, alpha=0.25):
    """Barbie-style ssGSEA on within-sample ranks (as in the estimate R pkg)."""
    scores = {}
    genes = expr_log.index
    inset = genes.isin(gene_set)
    n = len(genes)
    for col in expr_log.columns:
        r = stats.rankdata(expr_log[col].values)  # ascending ranks
        order = np.argsort(-r)  # descending by expression
        inset_ord = inset[order]
        rw = np.abs(r[order]) ** alpha
        pin_num = np.where(inset_ord, rw, 0.0).cumsum()
        pin = pin_num / pin_num[-1]
        pout = np.where(~inset_ord, 1.0, 0.0).cumsum() / (n - inset.sum())
        scores[col] = float((pin - pout).sum())
    return pd.Series(scores)


def sig_mean(expr_log, genes, name, log):
    present = [g for g in genes if g in expr_log.index]
    missing = sorted(set(genes) - set(present))
    log.append(f"  {name}: {len(present)}/{len(genes)} genes present"
               + (f" (missing: {', '.join(missing)})" if missing else ""))
    return expr_log.loc[present].mean(axis=0)


def partial_spearman(x, y, z):
    """Partial Spearman rho of x,y controlling z (residualize ranks)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    def resid(a, b):
        b1 = np.column_stack([np.ones_like(b), b])
        beta, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ beta
    ex, ey = resid(rx, rz), resid(ry, rz)
    r, _ = stats.pearsonr(ex, ey)
    n = len(x)
    dof = n - 3
    t = r * np.sqrt(dof / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    return r, p, n


def analyze(tag, fpkm_path, drop_cols=()):
    log = [f"== {tag} =="]
    fpkm = load_fpkm(fpkm_path, drop_cols=drop_cols)
    log.append(f"  matrix: {fpkm.shape[0]} genes x {fpkm.shape[1]} samples")
    if "TACSTD2" not in fpkm.index:
        log.append("  TACSTD2 not found — SKIPPING")
        return None, log
    expr = np.log2(fpkm + 1)

    sets = load_gmt(GMT)
    for k in sets:
        n_present = expr.index.isin(sets[k]).sum()
        log.append(f"  ESTIMATE {k}: {n_present}/{len(sets[k])} genes present")
    stromal = ssgsea(expr, sets["StromalSignature"])
    immune = ssgsea(expr, sets["ImmuneSignature"])
    est = stromal + immune
    # Yoshihara purity formula (Affymetrix-calibrated; out-of-[0,1] values on
    # RNA-seq indicate calibration mismatch — reported for reference only).
    purity = np.cos(0.6049872018 + 0.0001467884 * est)
    log.append(f"  ESTIMATE purity (uncalibrated for RNA-seq): median "
               f"{purity.median():.3f} (range {purity.min():.3f}-{purity.max():.3f})")

    scores = pd.DataFrame({
        "TACSTD2_log2FPKM": expr.loc["TACSTD2"],
        "ESTIMATE_StromalScore": stromal,
        "ESTIMATE_ImmuneScore": immune,
        "ESTIMATE_Score": est,
        "ESTIMATE_purity": purity,
        "IFNG6_Ayers": sig_mean(expr, IFNG6, "IFNG6_Ayers", log),
        "GEP18_TcellInflamed": sig_mean(expr, GEP18, "GEP18", log),
        "CYT_Rooney": sig_mean(expr, CYT, "CYT", log),
        "CD8A": expr.loc["CD8A"] if "CD8A" in expr.index else np.nan,
    })
    scores.index.name = "sample"
    scores.to_csv(os.path.join(OUT, f"scores_{tag}.csv"))

    rows = []
    x = scores["TACSTD2_log2FPKM"].values
    # Partial Spearman is invariant to monotone transforms of the covariate,
    # so adjust on the raw ESTIMATE score (purity is a monotone decreasing
    # function of it); this avoids the Affymetrix cos-formula calibration and
    # clipping issues on RNA-seq. Sign flipped so covariate is purity-directed.
    z = -scores["ESTIMATE_Score"].values
    # Sensitivity covariate: stromal score only (purity proxy that excludes
    # the immune component, avoiding over-adjustment/circularity for immune
    # outcomes, at the cost of being a weaker purity proxy).
    zs = -scores["ESTIMATE_StromalScore"].values
    for name in ["ESTIMATE_ImmuneScore", "IFNG6_Ayers",
                 "GEP18_TcellInflamed", "CYT_Rooney", "CD8A"]:
        y = scores[name].values
        rho, p = stats.spearmanr(x, y)
        prho, pp, n = partial_spearman(x, y, z)
        srho, sp, _ = partial_spearman(x, y, zs)
        rows.append({"dataset": tag, "immune_score": name, "n": n,
                     "spearman_rho": round(rho, 4), "spearman_p": round(p, 4),
                     "purity_adj_rho": round(prho, 4),
                     "purity_adj_p": round(pp, 4),
                     "stromal_adj_rho": round(srho, 4),
                     "stromal_adj_p": round(sp, 4)})
    rho_pur, p_pur = stats.spearmanr(x, z)
    rows.append({"dataset": tag, "immune_score": "purity_proxy(-ESTIMATEScore)",
                 "n": len(x), "spearman_rho": round(rho_pur, 4),
                 "spearman_p": round(p_pur, 4),
                 "purity_adj_rho": np.nan, "purity_adj_p": np.nan,
                 "stromal_adj_rho": np.nan, "stromal_adj_p": np.nan})
    res = pd.DataFrame(rows)
    for r in rows:
        log.append(f"  TACSTD2 vs {r['immune_score']}: rho={r['spearman_rho']}"
                   f" (p={r['spearman_p']}), purity-adj rho="
                   f"{r['purity_adj_rho']} (p={r['purity_adj_p']}), "
                   f"stromal-adj rho={r['stromal_adj_rho']} "
                   f"(p={r['stromal_adj_p']})")
    return res, log


def main():
    all_res, all_log = [], []
    jobs = [
        ("GSE253564_pretreatment",
         os.path.join(DATA, "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"),
         ("Entrez.ID",)),
        ("GSE248378_post_durvalumab",
         os.path.join(DATA, "GSE248378_Durva_Post_FPKMs.txt.gz"), ()),
    ]
    for tag, path, drop in jobs:
        res, log = analyze(tag, path, drop_cols=drop)
        all_log.extend(log)
        if res is not None:
            all_res.append(res)
    out = pd.concat(all_res, ignore_index=True)
    out.to_csv(os.path.join(OUT, "tacstd2_immune_correlations.csv"), index=False)
    with open(os.path.join(OUT, "analysis_log.txt"), "w") as fh:
        fh.write("\n".join(all_log) + "\n")
    print("\n".join(all_log))
    print("\nWrote:", os.path.join(OUT, "tacstd2_immune_correlations.csv"))


if __name__ == "__main__":
    main()
