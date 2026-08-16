"""Core test of the Bessede hypothesis: TACSTD2 (TROP2) / CLDN4 high == TLS low.

The analysis is built around one obvious confounder: TACSTD2 and CLDN4 are
epithelial (tumour-cell) genes, so in bulk RNA anything immune will correlate
negatively with them simply because tumour cells and immune cells compete for
the same library. Every association is therefore reported three ways:

  1. crude Spearman
  2. partial Spearman given tumour purity (ABSOLUTE in TCGA, pan-epithelial
     score elsewhere; leukocyte fraction from methylation as a second control)
  3. position of the gene in the transcriptome-wide null of the same statistic
     (is TACSTD2 special, or just an average epithelial gene?)

Outputs land in results/opus_tls/tables/.
"""

from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
from scipy import stats

from opus_tls_lib import (
    COMPARATOR_GENES, CONTROL_SIGNATURES, FOCUS_GENES, PROC_DIR, RESULTS_DIR,
    SIGNATURES, TLS_SIGNATURES, bh_fdr, genomewide_rank, inverse_normal,
    mann_whitney_effect, meta_fisher_z, ols_table, partial_spearman,
    rank_signature_scores, signature_scores, spearman_ci, zscore_rows,
)

TABLES = os.path.join(RESULTS_DIR, "tables")
PRIMARY_SIGS = ["TLS_Cabrita", "B_cell", "Plasma_cell"]


def load_cohorts() -> dict:
    with open(os.path.join(PROC_DIR, "tcga.pkl"), "rb") as fh:
        tcga = pickle.load(fh)
    with open(os.path.join(PROC_DIR, "geo_cohorts.pkl"), "rb") as fh:
        geo = pickle.load(fh)
    cohorts = {}
    for proj, d in tcga.items():
        cohorts[proj] = {"expr": d["expr"], "pheno": d["pheno"], "kind": "atlas",
                         "platform": "RNA-seq (Xena GDC log2 TPM+1)",
                         "note": "TCGA primary tumours", "expr_normal": d["expr_normal"]}
    cohorts.update(geo)
    return cohorts


def annotate_scores(c: dict) -> dict:
    scores, used = signature_scores(c["expr"])
    c["scores"] = scores
    c["scores_rank"] = rank_signature_scores(c["expr"])
    c["sig_genes_used"] = used
    ph = c["pheno"].reindex(scores.index)
    ph["epithelial_score"] = scores["Epithelial"]
    if "purity_absolute" in ph.columns:
        ph["purity"] = ph["purity_absolute"]
    else:
        # rank-scaled epithelial score as a tumour-content proxy where ABSOLUTE
        # purity is unavailable
        ph["purity"] = ph["epithelial_score"].rank(pct=True)
    c["pheno"] = ph
    return c


# ---------------------------------------------------------------------------
# 1. crude and adjusted correlations
# ---------------------------------------------------------------------------

def correlation_table(cohorts: dict) -> pd.DataFrame:
    rows = []
    for name, c in cohorts.items():
        expr, scores, ph = c["expr"], c["scores"], c["pheno"]
        genes = [g for g in FOCUS_GENES + COMPARATOR_GENES if g in expr.index]
        for g in genes:
            x = expr.loc[g, scores.index]
            for sig in scores.columns:
                y = scores[sig]
                rho, p, ci, n = spearman_ci(x.values, y.values, n_boot=1000)
                row = {"cohort": name, "kind": c["kind"], "gene": g, "signature": sig,
                       "rho": rho, "p": p, "ci_low": ci[0] if isinstance(ci, tuple) else np.nan,
                       "ci_high": ci[1] if isinstance(ci, tuple) else np.nan, "n": n}
                # adjusted for the pan-epithelial tumour-content proxy
                r_ep, p_ep, n_ep = partial_spearman(
                    x.values, y.values, ph[["epithelial_score"]].values)
                row.update({"rho_adj_epithelial": r_ep, "p_adj_epithelial": p_ep, "n_adj_epithelial": n_ep})
                # adjusted for ABSOLUTE purity where available
                if "purity_absolute" in ph.columns:
                    r_pu, p_pu, n_pu = partial_spearman(
                        x.values, y.values, ph[["purity_absolute"]].values)
                    row.update({"rho_adj_purity": r_pu, "p_adj_purity": p_pu, "n_adj_purity": n_pu})
                if "leukocyte_fraction" in ph.columns:
                    r_lf, p_lf, n_lf = partial_spearman(
                        x.values, y.values, ph[["purity_absolute", "leukocyte_fraction"]].values)
                    row.update({"rho_adj_purity_leuk": r_lf, "p_adj_purity_leuk": p_lf,
                                "n_adj_purity_leuk": n_lf})
                rows.append(row)
    df = pd.DataFrame(rows)
    for col in [c for c in df.columns if c.startswith("p")]:
        df["q" + col[1:]] = np.nan
    # FDR within cohort x gene family (all signatures tested for a given gene)
    for (coh, gene), idx in df.groupby(["cohort", "gene"]).groups.items():
        for pcol in [c for c in df.columns if c.startswith("p_") or c == "p"]:
            qcol = "q" + pcol[1:]
            df.loc[idx, qcol] = bh_fdr(df.loc[idx, pcol].values)
    return df


# ---------------------------------------------------------------------------
# 2. confounding structure: how much of the crude signal is tumour content?
# ---------------------------------------------------------------------------

def confounding_table(cohorts: dict) -> pd.DataFrame:
    rows = []
    for name, c in cohorts.items():
        expr, scores, ph = c["expr"], c["scores"], c["pheno"]
        for g in FOCUS_GENES:
            if g not in expr.index:
                continue
            x = expr.loc[g, scores.index]
            rec = {"cohort": name, "gene": g, "n": len(x)}
            rec["rho_gene_vs_epithelial"] = stats.spearmanr(x, scores["Epithelial"])[0]
            if "purity_absolute" in ph.columns:
                ok = ph["purity_absolute"].notna()
                rec["rho_gene_vs_purity"] = stats.spearmanr(x[ok], ph.loc[ok, "purity_absolute"])[0]
                rec["rho_epithelial_vs_purity"] = stats.spearmanr(
                    scores.loc[ok, "Epithelial"], ph.loc[ok, "purity_absolute"])[0]
                for sig in PRIMARY_SIGS:
                    rec[f"rho_{sig}_vs_purity"] = stats.spearmanr(
                        scores.loc[ok, sig], ph.loc[ok, "purity_absolute"])[0]
            for sig in PRIMARY_SIGS:
                rec[f"rho_{sig}_vs_epithelial"] = stats.spearmanr(scores[sig], scores["Epithelial"])[0]
            rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. transcriptome-wide null
# ---------------------------------------------------------------------------

def transcriptome_null(cohorts: dict, min_expr_frac: float = 0.5) -> tuple[pd.DataFrame, dict]:
    """Correlate every expressed gene with each primary signature, then locate
    TACSTD2/CLDN4 in that distribution (crude and purity-adjusted)."""
    rows, dists = [], {}
    for name, c in cohorts.items():
        expr, scores, ph = c["expr"], c["scores"], c["pheno"]
        if expr.shape[1] < 40:
            continue  # null distributions are meaningless in n<40 cohorts
        expressed = (expr > 0).mean(axis=1) >= min_expr_frac
        E = expr.loc[expressed, scores.index]
        R = E.T.rank()  # samples x genes ranks
        for sig in PRIMARY_SIGS + ["T_cell_CD8", "Myeloid"]:
            y = scores[sig]
            yr = stats.rankdata(y.values)
            rho_all = _fast_corr(R.values, yr)
            rho_all = pd.Series(rho_all, index=E.index)
            dists[(name, sig, "crude")] = rho_all
            # purity-adjusted version
            cov = ph["purity_absolute"] if "purity_absolute" in ph.columns else ph["epithelial_score"]
            covname = "ABSOLUTE purity" if "purity_absolute" in ph.columns else "epithelial score"
            ok = cov.notna().values
            cr = stats.rankdata(cov[ok].values)
            X = np.column_stack([np.ones(ok.sum()), cr])
            Rr = R.values[ok]
            Rr = Rr - X @ np.linalg.lstsq(X, Rr, rcond=None)[0]
            yy = stats.rankdata(y.values[ok])
            yy = yy - X @ np.linalg.lstsq(X, yy, rcond=None)[0]
            rho_adj = pd.Series(_fast_corr(Rr, yy), index=E.index)
            dists[(name, sig, "adjusted")] = rho_adj
            for g in FOCUS_GENES + ["EPCAM", "KRT8", "CDH1", "MUC1"]:
                if g not in rho_all.index:
                    continue
                rec = {"cohort": name, "signature": sig, "gene": g,
                       "n_samples": E.shape[1], "n_genes": E.shape[0],
                       "adjustment": covname,
                       "rho_crude": rho_all[g], "rho_adjusted": rho_adj[g]}
                rec.update({f"crude_{k}": v for k, v in genomewide_rank(rho_all[g], rho_all).items()})
                rec.update({f"adj_{k}": v for k, v in genomewide_rank(rho_adj[g], rho_adj).items()})
                rows.append(rec)
    return pd.DataFrame(rows), dists


def _fast_corr(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson correlation of every column of X with y (X already rank/residual)."""
    Xc = X - X.mean(axis=0)
    yc = y - y.mean()
    num = Xc.T @ yc
    den = np.sqrt((Xc ** 2).sum(axis=0) * (yc ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


# ---------------------------------------------------------------------------
# 4. quartile contrasts, purity-stratified
# ---------------------------------------------------------------------------

def quartile_table(cohorts: dict) -> pd.DataFrame:
    rows = []
    for name, c in cohorts.items():
        expr, scores, ph = c["expr"], c["scores"], c["pheno"]
        if expr.shape[1] < 40:
            continue
        for g in FOCUS_GENES:
            if g not in expr.index:
                continue
            x = expr.loc[g, scores.index]
            grp = pd.qcut(x, 4, labels=["Q1", "Q2", "Q3", "Q4"])
            for sig in PRIMARY_SIGS + CONTROL_SIGNATURES:
                y = scores[sig]
                res = mann_whitney_effect(y[grp == "Q4"], y[grp == "Q1"])
                rows.append({"cohort": name, "gene": g, "signature": sig, "stratum": "all",
                             **res})
                # stratified within tertiles of tumour purity
                pur = ph["purity"]
                tert = pd.qcut(pur.rank(method="first"), 3, labels=["P1", "P2", "P3"])
                zs, ws = [], []
                for t in ["P1", "P2", "P3"]:
                    sel = tert == t
                    xs = x[sel]
                    if xs.notna().sum() < 20:
                        continue
                    gq = pd.qcut(xs, 4, labels=["Q1", "Q2", "Q3", "Q4"])
                    r = mann_whitney_effect(y[sel][gq == "Q4"], y[sel][gq == "Q1"])
                    rows.append({"cohort": name, "gene": g, "signature": sig,
                                 "stratum": f"purity_{t}", **r})
                    if np.isfinite(r["p"]) and r["n1"] and r["n2"]:
                        z = stats.norm.isf(r["p"] / 2) * np.sign(r["rbc"])
                        zs.append(z)
                        ws.append(np.sqrt(r["n1"] + r["n2"]))
                if len(zs) >= 2:
                    z_st = float(np.sum(np.array(ws) * np.array(zs)) / np.sqrt(np.sum(np.array(ws) ** 2)))
                    rows.append({"cohort": name, "gene": g, "signature": sig,
                                 "stratum": "purity_stratified_stouffer",
                                 "p": float(2 * stats.norm.sf(abs(z_st))), "U": np.nan,
                                 "rbc": np.nan, "hl_shift": np.nan, "n1": np.nan, "n2": np.nan,
                                 "z": z_st})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. multivariable models
# ---------------------------------------------------------------------------

def model_table(cohorts: dict) -> pd.DataFrame:
    out = []
    for name, c in cohorts.items():
        expr, scores, ph = c["expr"], c["scores"], c["pheno"]
        if expr.shape[1] < 60:
            continue
        for g in FOCUS_GENES:
            if g not in expr.index:
                continue
            x = inverse_normal(expr.loc[g, scores.index])
            for sig in PRIMARY_SIGS:
                y = inverse_normal(scores[sig])
                X = pd.DataFrame({g: x}, index=scores.index)
                out.append(ols_table(y, X, label=f"{name}|{sig}|unadjusted"))
                Xc = X.copy()
                if "purity_absolute" in ph.columns and ph["purity_absolute"].notna().sum() > 50:
                    Xc["purity"] = ph["purity_absolute"]
                else:
                    Xc["epithelial"] = inverse_normal(ph["epithelial_score"])
                out.append(ols_table(y, Xc, label=f"{name}|{sig}|+purity"))
                Xf = Xc.copy()
                if "age" in ph.columns:
                    Xf["age"] = pd.to_numeric(ph["age"], errors="coerce")
                if "sex" in ph.columns:
                    Xf["male"] = (ph["sex"].astype(str).str.lower().str[0] == "m").astype(float)
                if "stage" in ph.columns:
                    Xf["stage"] = pd.to_numeric(ph["stage"], errors="coerce")
                if "leukocyte_fraction" in ph.columns:
                    Xf["leukocyte_fraction"] = ph["leukocyte_fraction"]
                out.append(ols_table(y, Xf, label=f"{name}|{sig}|+purity+clinical"))
    res = [d for d in out if len(d)]
    return pd.concat(res, ignore_index=True) if res else pd.DataFrame()


# ---------------------------------------------------------------------------
# 6. orthogonal (non-signature) immune read-outs in TCGA
# ---------------------------------------------------------------------------

def orthogonal_table(cohorts: dict) -> pd.DataFrame:
    rows = []
    targets = ["cib_B_naive", "cib_B_memory", "cib_B_total", "cib_plasma", "cib_Tfh",
               "cib_CD8", "leukocyte_fraction", "tmb_nonsilent"]
    for name, c in cohorts.items():
        ph = c["pheno"]
        if "cib_plasma" not in ph.columns:
            continue
        expr, scores = c["expr"], c["scores"]
        for g in FOCUS_GENES:
            if g not in expr.index:
                continue
            x = expr.loc[g, scores.index]
            for t in targets:
                if t not in ph.columns:
                    continue
                y = ph[t]
                rho, p, ci, n = spearman_ci(x.values, y.values, n_boot=1000)
                r_pu, p_pu, n_pu = partial_spearman(
                    x.values, y.values, ph[["purity_absolute"]].values)
                rows.append({"cohort": name, "gene": g, "readout": t, "rho": rho, "p": p,
                             "ci_low": ci[0] if isinstance(ci, tuple) else np.nan,
                             "ci_high": ci[1] if isinstance(ci, tuple) else np.nan, "n": n,
                             "rho_adj_purity": r_pu, "p_adj_purity": p_pu, "n_adj": n_pu})
    df = pd.DataFrame(rows)
    if len(df):
        df["q"] = bh_fdr(df["p"].values)
        df["q_adj_purity"] = bh_fdr(df["p_adj_purity"].values)
    return df


# ---------------------------------------------------------------------------
# 7. cross-cohort meta-analysis
# ---------------------------------------------------------------------------

def meta_table(corr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in FOCUS_GENES:
        for sig in TLS_SIGNATURES + CONTROL_SIGNATURES:
            for stat, ncol in [("rho", "n"), ("rho_adj_epithelial", "n_adj_epithelial")]:
                sub = corr[(corr.gene == gene) & (corr.signature == sig)].dropna(subset=[stat])
                if sub.empty:
                    continue
                m = meta_fisher_z(sub[stat].tolist(), sub[ncol].tolist())
                if not m:
                    continue
                rows.append({"gene": gene, "signature": sig, "statistic": stat,
                             "cohorts": ",".join(sub.cohort), **{
                                 k: (v if not isinstance(v, tuple) else v) for k, v in m.items()}})
    df = pd.DataFrame(rows)
    for c in ["ci_fe", "ci_re"]:
        if c in df.columns:
            df[c + "_low"] = df[c].map(lambda t: t[0])
            df[c + "_high"] = df[c].map(lambda t: t[1])
            df = df.drop(columns=c)
    return df


def main() -> None:
    cohorts = load_cohorts()
    for name in list(cohorts):
        cohorts[name] = annotate_scores(cohorts[name])
        print(f"[scores] {name}: n={cohorts[name]['expr'].shape[1]}")

    cov = pd.DataFrame([
        {"cohort": n, "signature": s, "n_genes_in_set": len(SIGNATURES[s]),
         "n_genes_found": len(g), "genes_found": ";".join(g)}
        for n, c in cohorts.items() for s, g in c["sig_genes_used"].items()
    ])
    cov.to_csv(os.path.join(TABLES, "signature_gene_coverage.csv"), index=False)

    print("[1/7] correlations ...", flush=True)
    corr = correlation_table(cohorts)
    corr.to_csv(os.path.join(TABLES, "correlations_by_cohort.csv"), index=False)

    print("[2/7] confounding structure ...", flush=True)
    confounding_table(cohorts).to_csv(os.path.join(TABLES, "confounding_structure.csv"), index=False)

    print("[3/7] transcriptome-wide null ...", flush=True)
    tn, dists = transcriptome_null(cohorts)
    tn.to_csv(os.path.join(TABLES, "transcriptome_wide_null.csv"), index=False)
    with open(os.path.join(PROC_DIR, "null_distributions.pkl"), "wb") as fh:
        pickle.dump(dists, fh, protocol=4)

    print("[4/7] quartile contrasts ...", flush=True)
    quartile_table(cohorts).to_csv(os.path.join(TABLES, "quartile_contrasts.csv"), index=False)

    print("[5/7] multivariable models ...", flush=True)
    model_table(cohorts).to_csv(os.path.join(TABLES, "multivariable_models.csv"), index=False)

    print("[6/7] orthogonal readouts ...", flush=True)
    orthogonal_table(cohorts).to_csv(os.path.join(TABLES, "orthogonal_readouts_tcga.csv"), index=False)

    print("[7/7] meta-analysis ...", flush=True)
    meta_table(corr).to_csv(os.path.join(TABLES, "meta_analysis.csv"), index=False)

    # cache scored cohorts for the outcome / figure scripts
    slim = {n: {"scores": c["scores"], "scores_rank": c["scores_rank"], "pheno": c["pheno"],
                "kind": c["kind"], "platform": c["platform"], "note": c.get("note", ""),
                "focus_expr": c["expr"].reindex(FOCUS_GENES + COMPARATOR_GENES).dropna(how="all")}
            for n, c in cohorts.items()}
    with open(os.path.join(PROC_DIR, "scored_cohorts.pkl"), "wb") as fh:
        pickle.dump(slim, fh, protocol=4)
    print("[done] tables written to", TABLES)


if __name__ == "__main__":
    main()
