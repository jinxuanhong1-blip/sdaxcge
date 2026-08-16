"""Finish leftover tables from 03, then ICI outcomes + histology splits.

Does NOT re-run the 1000-bootstrap Spearman table. Uses the already-written
correlations_by_cohort.csv plus a cheap re-score of the pickled expression
objects for CIBERSORT / ICI / histology analyses.
"""

from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
from scipy import stats

from opus_tls_lib import (
    FOCUS_GENES, PROC_DIR, RESULTS_DIR, SIGNATURES, TLS_SIGNATURES,
    CONTROL_SIGNATURES, bh_fdr, fmt_p, mann_whitney_effect, meta_fisher_z,
    partial_spearman, signature_scores, spearman_ci,
)
from importlib.machinery import SourceFileLoader

core = SourceFileLoader(
    "core", os.path.join(os.path.dirname(__file__), "03_core_association.py")
).load_module()

TABLES = os.path.join(RESULTS_DIR, "tables")
PRIMARY = ["TLS_Cabrita", "TLS_12chemokine", "TLS_imprint", "B_cell",
           "Plasma_cell", "Tfh", "T_cell_CD8", "IFNg_Ayers"]


def load_scored() -> dict:
    cohorts = core.load_cohorts()
    return {n: core.annotate_scores(c) for n, c in cohorts.items()}


# ---------------------------------------------------------------------------
# leftover from 03
# ---------------------------------------------------------------------------

def write_leftovers(cohorts: dict) -> None:
    corr = pd.read_csv(os.path.join(TABLES, "correlations_by_cohort.csv"))
    core.orthogonal_table(cohorts).to_csv(
        os.path.join(TABLES, "orthogonal_readouts_tcga.csv"), index=False)
    core.meta_table(corr).to_csv(os.path.join(TABLES, "meta_analysis.csv"), index=False)

    # histology-aware meta: LUAD-like vs LUSC-like vs ICI
    rows = []
    groups = {
        "LUAD_like": ["TCGA-LUAD", "GSE72094"],
        "LUSC_like": ["TCGA-LUSC"],
        "mixed_NSCLC_atlas": ["GSE81089"],
        "ICI_all": ["GSE135222", "GSE126044", "GSE207422"],
        "large_atlas": ["TCGA-LUAD", "TCGA-LUSC", "GSE81089", "GSE72094"],
    }
    for gname, members in groups.items():
        for gene in FOCUS_GENES:
            for sig in PRIMARY:
                for stat, ncol in [("rho", "n"), ("rho_adj_epithelial", "n_adj_epithelial"),
                                   ("rho_adj_purity", "n_adj_purity")]:
                    sub = corr[(corr.gene == gene) & (corr.signature == sig)
                               & (corr.cohort.isin(members))].dropna(subset=[stat])
                    if len(sub) < 1:
                        continue
                    m = meta_fisher_z(sub[stat].tolist(), sub[ncol].fillna(sub["n"]).tolist())
                    if not m:
                        continue
                    rows.append({"group": gname, "gene": gene, "signature": sig,
                                 "statistic": stat, "cohorts": ",".join(sub.cohort),
                                 **{k: (v if not isinstance(v, tuple) else v) for k, v in m.items()}})
    df = pd.DataFrame(rows)
    for c in ["ci_fe", "ci_re"]:
        if c in df.columns:
            df[c + "_low"] = df[c].map(lambda t: t[0] if isinstance(t, tuple) else np.nan)
            df[c + "_high"] = df[c].map(lambda t: t[1] if isinstance(t, tuple) else np.nan)
            df = df.drop(columns=c)
    df.to_csv(os.path.join(TABLES, "meta_by_histology_group.csv"), index=False)


# ---------------------------------------------------------------------------
# decision table: where does TACSTD2-high = immune-low survive purity?
# ---------------------------------------------------------------------------

def decision_table(corr: pd.DataFrame) -> pd.DataFrame:
    """One row per cohort x signature. Uses ABSOLUTE purity when present,
    otherwise the epithelial-score partial correlation."""
    rows = []
    for (coh, gene, sig), g in corr.groupby(["cohort", "gene", "signature"]):
        r = g.iloc[0]
        if pd.notna(r.get("rho_adj_purity")):
            method, rho, p, n = "ABSOLUTE_purity", r["rho_adj_purity"], r["p_adj_purity"], r["n_adj_purity"]
        else:
            method, rho, p, n = "epithelial_score", r["rho_adj_epithelial"], r["p_adj_epithelial"], r["n_adj_epithelial"]
        holds = bool(np.isfinite(rho) and rho < 0 and p < 0.05 and n >= 40)
        opposite = bool(np.isfinite(rho) and rho > 0 and p < 0.05 and n >= 40)
        rows.append({
            "cohort": coh, "kind": r["kind"], "gene": gene, "signature": sig,
            "n": int(n) if np.isfinite(n) else int(r["n"]),
            "adjustment": method,
            "rho_crude": r["rho"], "p_crude": r["p"],
            "rho_adj": rho, "p_adj": p,
            "holds_after_purity": holds,
            "opposite_after_purity": opposite,
            "verdict": ("HOLDS" if holds else ("OPPOSITE" if opposite else
                        ("UNDERPOWERED" if float(r["n"]) < 40 else "NO_EVIDENCE"))),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# histology split of GSE81089 (independent mixed NSCLC atlas)
# ---------------------------------------------------------------------------

def histology_split(cohorts: dict) -> pd.DataFrame:
    rows = []
    c = cohorts["GSE81089"]
    expr, scores, ph = c["expr"], c["scores"], c["pheno"]
    for hist, mask in [
        ("AC", ph["histology"] == "AC"),
        ("SqCC", ph["histology"] == "SqCC"),
        ("LCC", ph["histology"] == "LCC"),
    ]:
        idx = ph.index[mask]
        if len(idx) < 15:
            continue
        for g in FOCUS_GENES:
            x = expr.loc[g, idx]
            for sig in PRIMARY:
                y = scores.loc[idx, sig]
                rho, p, ci, n = spearman_ci(x.values, y.values, n_boot=800)
                r_ep, p_ep, n_ep = partial_spearman(
                    x.values, y.values, ph.loc[idx, ["epithelial_score"]].values)
                rows.append({"cohort": "GSE81089", "stratum": hist, "gene": g,
                             "signature": sig, "n": n, "rho": rho, "p": p,
                             "ci_low": ci[0] if isinstance(ci, tuple) else np.nan,
                             "ci_high": ci[1] if isinstance(ci, tuple) else np.nan,
                             "rho_adj_epithelial": r_ep, "p_adj_epithelial": p_ep})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GSE207422 restricted to pre-treatment biopsies
# ---------------------------------------------------------------------------

def gse207422_pre(cohorts: dict) -> pd.DataFrame:
    c = cohorts["GSE207422"]
    ph = c["pheno"]
    pre = ph.index[ph["timepoint"] == "pre"]
    rows = []
    for g in FOCUS_GENES:
        x = c["expr"].loc[g, pre]
        for sig in PRIMARY:
            y = c["scores"].loc[pre, sig]
            rho, p, ci, n = spearman_ci(x.values, y.values, n_boot=800)
            r_ep, p_ep, n_ep = partial_spearman(
                x.values, y.values, ph.loc[pre, ["epithelial_score"]].values)
            rows.append({"subset": "pre_only", "gene": g, "signature": sig, "n": n,
                         "rho": rho, "p": p,
                         "ci_low": ci[0] if isinstance(ci, tuple) else np.nan,
                         "ci_high": ci[1] if isinstance(ci, tuple) else np.nan,
                         "rho_adj_epithelial": r_ep, "p_adj_epithelial": p_ep})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# ICI clinical endpoints (honest: small n)
# ---------------------------------------------------------------------------

def ici_outcomes(cohorts: dict) -> pd.DataFrame:
    rows = []

    # GSE135222: PFS Cox + DCB
    c = cohorts["GSE135222"]
    ph, expr = c["pheno"], c["expr"]
    try:
        from lifelines import CoxPHFitter
        for g in FOCUS_GENES:
            df = pd.DataFrame({
                "y": expr.loc[g, ph.index],
                "T": ph["pfs_time"],
                "E": ph["pfs_event"],
                "epithelial": ph["epithelial_score"],
            }).dropna()
            if df["E"].sum() < 5:
                continue
            for label, cols in [("unadjusted", ["y", "T", "E"]),
                                ("+epithelial", ["y", "epithelial", "T", "E"])]:
                m = CoxPHFitter()
                m.fit(df[cols], duration_col="T", event_col="E")
                s = m.summary.loc["y"]
                rows.append({"cohort": "GSE135222", "endpoint": "PFS_Cox",
                             "gene": g, "model": label, "n": len(df),
                             "n_events": int(df["E"].sum()),
                             "hr": float(s["exp(coef)"]),
                             "hr_low": float(s["exp(coef) lower 95%"]),
                             "hr_high": float(s["exp(coef) upper 95%"]),
                             "p": float(s["p"])})
    except Exception as e:
        rows.append({"cohort": "GSE135222", "endpoint": "PFS_Cox",
                     "note": f"cox failed: {e}"})

    for g in FOCUS_GENES:
        dcb = ph["dcb_6mo"]
        hi = expr.loc[g, ph.index][dcb == 1]
        lo = expr.loc[g, ph.index][dcb == 0]
        r = mann_whitney_effect(hi, lo)
        rows.append({"cohort": "GSE135222", "endpoint": "DCB_6mo",
                     "gene": g, "model": "MW_DCB_vs_NDB", **r})

    # GSE126044: responder vs non-responder
    c = cohorts["GSE126044"]
    ph, expr = c["pheno"], c["expr"]
    for g in FOCUS_GENES:
        r = mann_whitney_effect(expr.loc[g, ph.index][ph["responder"] == 1],
                                expr.loc[g, ph.index][ph["responder"] == 0])
        rows.append({"cohort": "GSE126044", "endpoint": "RECIST_responder",
                     "gene": g, "model": "MW_R_vs_NR", **r})
        # also: is TLS itself associated with response? (positive control)
        for sig in ["TLS_Cabrita", "B_cell", "T_cell_CD8"]:
            r2 = mann_whitney_effect(c["scores"].loc[ph.index, sig][ph["responder"] == 1],
                                     c["scores"].loc[ph.index, sig][ph["responder"] == 0])
            rows.append({"cohort": "GSE126044", "endpoint": "RECIST_responder",
                         "gene": sig, "model": "MW_R_vs_NR_signature", **r2})

    # GSE207422: MPR, pre-treatment only
    c = cohorts["GSE207422"]
    ph, expr = c["pheno"], c["expr"]
    pre = ph["timepoint"] == "pre"
    for g in FOCUS_GENES:
        r = mann_whitney_effect(expr.loc[g, ph.index][pre & (ph["mpr"] == 1)],
                                expr.loc[g, ph.index][pre & (ph["mpr"] == 0)])
        rows.append({"cohort": "GSE207422", "endpoint": "MPR_pre",
                     "gene": g, "model": "MW_MPR_vs_NMPR", **r})
    return pd.DataFrame(rows)


def cache_slim(cohorts: dict) -> None:
    slim = {n: {"scores": c["scores"], "pheno": c["pheno"], "kind": c["kind"],
                "platform": c["platform"], "note": c.get("note", ""),
                "focus_expr": c["expr"].reindex(FOCUS_GENES).dropna(how="all")}
            for n, c in cohorts.items()}
    with open(os.path.join(PROC_DIR, "scored_cohorts.pkl"), "wb") as fh:
        pickle.dump(slim, fh, protocol=4)


def main() -> None:
    print("[load] scoring cohorts ...", flush=True)
    cohorts = load_scored()
    cache_slim(cohorts)

    print("[leftovers] orthogonal + meta ...", flush=True)
    write_leftovers(cohorts)

    corr = pd.read_csv(os.path.join(TABLES, "correlations_by_cohort.csv"))
    dec = decision_table(corr)
    dec.to_csv(os.path.join(TABLES, "where_it_holds.csv"), index=False)
    print("\n=== TACSTD2 holds after purity (n>=40, rho<0, p<0.05) ===")
    h = dec[(dec.gene == "TACSTD2") & (dec.signature.isin(PRIMARY))]
    print(h[["cohort", "signature", "n", "adjustment", "rho_adj", "p_adj", "verdict"]]
          .sort_values(["cohort", "signature"]).to_string(index=False))

    print("[histology] GSE81089 ...", flush=True)
    histology_split(cohorts).to_csv(os.path.join(TABLES, "gse81089_by_histology.csv"), index=False)

    print("[GSE207422] pre-only ...", flush=True)
    gse207422_pre(cohorts).to_csv(os.path.join(TABLES, "gse207422_pre_only.csv"), index=False)

    print("[ICI] outcomes ...", flush=True)
    ici_outcomes(cohorts).to_csv(os.path.join(TABLES, "ici_outcomes.csv"), index=False)
    print("[done]")


if __name__ == "__main__":
    main()
