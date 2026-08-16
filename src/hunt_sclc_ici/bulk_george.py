"""Bulk RNA-seq analysis: TACSTD2 / CLDN4 vs the SCLC immune (SCLC-I) axis.

Cohort: George et al. 2015 (Nature), 81 primary SCLC tumors, RNA-seq FPKM.
Source: cBioPortal datahub study `sclc_ucologne_2015` (public).
This is the discovery cohort used by Gay et al. 2021 (Cancer Cell) to define
SCLC-A/N/P/I. It is fully public, so results here are REAL (not simulated).

Honesty notes:
- Subtype calls here are a transparent, reproducible approximation of the Gay
  NMF classification, not the original NMF labels (those labels are not public
  per-sample for George). We sanity-check our subtype proportions against the
  proportions reported by Gay et al. for this exact cohort.
- Signature scores are unweighted mean z-scores (or the exact Rooney CYT
  formula), not proprietary weighted scores. Stated openly.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from . import signatures as sig

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "results" / "hunt_sclc_ici" / "data"
OUT = ROOT / "results" / "hunt_sclc_ici"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
LOGS = OUT / "logs"
for d in (TABLES, FIGS, LOGS):
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler(LOGS / "bulk_george.log", mode="w"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("bulk_george")

# Gay et al. 2021 reported subtype proportions for the George cohort (n=81)
GAY_GEORGE_PROPORTIONS = {"SCLC-A": 0.36, "SCLC-N": 0.31, "SCLC-I": 0.17, "SCLC-P": 0.16}


def load_expression() -> pd.DataFrame:
    """Return log2(FPKM+1), genes x samples, duplicate symbols collapsed."""
    fp = DATA / "data_mrna_seq_fpkm.txt"
    raw = pd.read_csv(fp, sep="\t")
    raw = raw.drop(columns=[c for c in ("Entrez_Gene_Id",) if c in raw.columns])
    raw = raw[raw["Hugo_Symbol"].notna()]
    # collapse duplicate gene symbols by max mean expression (keep most expressed)
    raw["_mean"] = raw.drop(columns=["Hugo_Symbol"]).mean(axis=1, numeric_only=True)
    raw = raw.sort_values("_mean", ascending=False).drop_duplicates("Hugo_Symbol")
    raw = raw.drop(columns=["_mean"]).set_index("Hugo_Symbol")
    raw = raw.apply(pd.to_numeric, errors="coerce")
    expr = np.log2(raw + 1.0)
    log.info("Loaded expression: %d genes x %d samples", expr.shape[0], expr.shape[1])
    return expr


def zscore_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """z-score each gene (row) across samples."""
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def mean_z_signature(z: pd.DataFrame, genes: list[str], name: str) -> pd.Series:
    present = [g for g in genes if g in z.index]
    missing = [g for g in genes if g not in z.index]
    if missing:
        log.info("Signature %s: %d/%d genes present; missing: %s",
                 name, len(present), len(genes), ",".join(missing))
    return z.loc[present].mean(axis=0).rename(name)


def cyt_score(expr: pd.DataFrame) -> pd.Series:
    """Rooney 2015 CYT = geometric mean of GZMA, PRF1 (computed in log2 space)."""
    present = [g for g in sig.CYT_GENES if g in expr.index]
    return expr.loc[present].mean(axis=0).rename("CYT_log2")


def _unit(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=1)


def assign_subtype(z: pd.DataFrame, inflamed_score: pd.Series) -> pd.DataFrame:
    """Assign SCLC-A/N/P/I per sample by an argmax over four standardized axes.

    Rule (transparent approximation of Gay NMF; NOT the original NMF labels):
      subtype = argmax over unit-variance-standardized {ASCL1, NEUROD1, POU2F3,
                inflamed-GEP} axis scores.

    Rationale: Gay et al. define SCLC-I as tumors with LOW ASCL1/NEUROD1/POU2F3
    and HIGH inflamed signature. Adding the inflamed GEP as a fourth competing
    axis (rather than an arbitrary triple-low threshold) reproduces that concept
    reproducibly and yields proportions close to the published George values.
    We also record a secondary "triple-low" call for sensitivity analysis.
    """
    tfs = ["ASCL1", "NEUROD1", "POU2F3"]
    zt = z.loc[tfs].T  # samples x TFs, already per-gene z-scored
    axes = pd.DataFrame({
        "SCLC-A": _unit(zt["ASCL1"]),
        "SCLC-N": _unit(zt["NEUROD1"]),
        "SCLC-P": _unit(zt["POU2F3"]),
        "SCLC-I": _unit(inflamed_score.reindex(zt.index)),
    })
    subtype = axes.idxmax(axis=1)

    out = pd.DataFrame(index=zt.index)
    out["z_ASCL1"] = zt["ASCL1"]
    out["z_NEUROD1"] = zt["NEUROD1"]
    out["z_POU2F3"] = zt["POU2F3"]
    out["axis_A"] = axes["SCLC-A"]
    out["axis_N"] = axes["SCLC-N"]
    out["axis_P"] = axes["SCLC-P"]
    out["axis_I"] = axes["SCLC-I"]
    out["max_TF_z"] = zt.max(axis=1)
    out["subtype"] = subtype
    # secondary sensitivity call: triple-low (all TF z < 0)
    out["triple_low"] = zt.max(axis=1) < 0
    return out


def bh(pvals: list[float]) -> np.ndarray:
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    if ok.sum():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def correlation_table(expr: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Spearman correlation of each gene of interest with each immune/axis score."""
    rows = []
    axis_cols = [c for c in scores.columns if c not in ("subtype",)]
    for gene in sig.GENES_OF_INTEREST:
        if gene not in expr.index:
            log.warning("Gene of interest %s not in matrix", gene)
            continue
        gvec = expr.loc[gene]
        for col in axis_cols:
            svec = scores[col]
            common = gvec.index.intersection(svec.dropna().index)
            x = gvec.loc[common].astype(float)
            y = svec.loc[common].astype(float)
            if len(common) < 5 or y.nunique() < 3:
                continue
            rho, p = stats.spearmanr(x, y)
            rows.append({"gene": gene, "axis": col, "n": len(common),
                         "spearman_rho": rho, "p_value": p})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value_BH"] = bh(df["p_value"].tolist())
        df = df.sort_values(["gene", "p_value"]).reset_index(drop=True)
    return df


def group_comparison(expr: pd.DataFrame, subtype: pd.Series) -> pd.DataFrame:
    """Compare gene-of-interest expression across subtypes and SCLC-I vs rest."""
    rows = []
    groups = ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]
    for gene in sig.GENES_OF_INTEREST + sig.CONTEXT_SURFACE_TARGETS:
        if gene not in expr.index:
            continue
        gvec = expr.loc[subtype.index].loc[:, :] if False else expr.loc[gene, subtype.index]
        gvec = gvec.astype(float)
        by = {g: gvec[subtype == g].values for g in groups if (subtype == g).any()}
        rec = {"gene": gene}
        for g in groups:
            v = by.get(g, np.array([]))
            rec[f"median_{g}"] = float(np.median(v)) if v.size else np.nan
            rec[f"n_{g}"] = int(v.size)
        # Kruskal-Wallis across available groups
        vals = [v for v in by.values() if v.size > 0]
        if len(vals) >= 3 and all(len(v) >= 2 for v in vals):
            try:
                rec["kruskal_H"], rec["kruskal_p"] = stats.kruskal(*vals)
            except ValueError:
                rec["kruskal_H"], rec["kruskal_p"] = np.nan, np.nan
        # SCLC-I vs rest (Mann-Whitney)
        i_vals = gvec[subtype == "SCLC-I"].values
        rest_vals = gvec[subtype != "SCLC-I"].values
        if i_vals.size >= 3 and rest_vals.size >= 3:
            u, p = stats.mannwhitneyu(i_vals, rest_vals, alternative="two-sided")
            rec["SCLC_I_vs_rest_U"] = u
            rec["SCLC_I_vs_rest_p"] = p
            rec["SCLC_I_median"] = float(np.median(i_vals))
            rec["rest_median"] = float(np.median(rest_vals))
            rec["log2FC_I_minus_rest"] = float(np.median(i_vals) - np.median(rest_vals))
            # rank-biserial effect size
            rec["cliffs_delta"] = 2 * u / (i_vals.size * rest_vals.size) - 1
        rows.append(rec)
    df = pd.DataFrame(rows)
    if "SCLC_I_vs_rest_p" in df:
        df["SCLC_I_vs_rest_q_BH"] = bh(df["SCLC_I_vs_rest_p"].tolist())
    return df


def run() -> dict:
    expr = load_expression()
    z = zscore_genes(expr)

    # scores
    scores = pd.DataFrame(index=expr.columns)
    for name, genes in sig.MEANZ_SIGNATURES.items():
        scores[name] = mean_z_signature(z, genes, name)
    scores["CYT_log2"] = cyt_score(expr)
    scores["CD8A_log2"] = expr.loc["CD8A"] if "CD8A" in expr.index else np.nan
    scores["NE_score"] = (mean_z_signature(z, sig.NE_UP, "NE_UP")
                          - mean_z_signature(z, sig.NON_NE_UP, "NON_NE_UP"))

    subt = assign_subtype(z, scores["GEP_Ayers_Tcell_inflamed"])
    scores = scores.join(subt[["axis_A", "axis_N", "axis_P", "axis_I",
                               "max_TF_z", "subtype", "triple_low"]])

    # add genes of interest raw log2 for convenience
    for g in sig.GENES_OF_INTEREST + ["ASCL1", "NEUROD1", "POU2F3", "YAP1"]:
        if g in expr.index:
            scores[f"log2_{g}"] = expr.loc[g]

    scores.index.name = "sample_id"
    scores.to_csv(TABLES / "george_sample_scores.tsv", sep="\t")

    # subtype proportions vs Gay
    props = subt["subtype"].value_counts(normalize=True).to_dict()
    prop_df = pd.DataFrame({
        "subtype": ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"],
        "n": [int((subt["subtype"] == s).sum()) for s in ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]],
        "proportion_this_analysis": [props.get(s, 0.0) for s in ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]],
        "proportion_Gay2021": [GAY_GEORGE_PROPORTIONS[s] for s in ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]],
    })
    prop_df.to_csv(TABLES / "george_subtype_proportions.tsv", sep="\t", index=False)
    log.info("Subtype proportions (this analysis vs Gay 2021):\n%s", prop_df.to_string(index=False))

    corr = correlation_table(expr, scores)
    corr.to_csv(TABLES / "george_correlations.tsv", sep="\t", index=False)
    log.info("Correlation table:\n%s", corr.to_string(index=False))

    grp = group_comparison(expr, subt["subtype"])
    grp.to_csv(TABLES / "george_group_comparison.tsv", sep="\t", index=False)
    log.info("Group comparison (subset):\n%s",
             grp[[c for c in grp.columns if c in
                  ("gene", "SCLC_I_median", "rest_median",
                   "log2FC_I_minus_rest", "SCLC_I_vs_rest_p", "SCLC_I_vs_rest_q_BH")]].to_string(index=False))

    summary = {
        "cohort": "George et al. 2015 SCLC (cBioPortal sclc_ucologne_2015)",
        "n_samples": int(expr.shape[1]),
        "n_genes": int(expr.shape[0]),
        "subtype_proportions": {s: props.get(s, 0.0) for s in ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]},
    }
    (TABLES / "george_summary.json").write_text(json.dumps(summary, indent=2))
    return {"expr": expr, "scores": scores, "subtype": subt["subtype"],
            "corr": corr, "grp": grp}


if __name__ == "__main__":
    run()
