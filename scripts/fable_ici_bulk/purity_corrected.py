#!/usr/bin/env python3
"""Purity-corrected, tumor-intrinsic re-analysis of TACSTD2/CLDN4 vs ICI response.

Tests the biological hypothesis (per user direction):
  * tumor-intrinsic TACSTD2 (TROP2) is HIGHER in non-responders (NR / NMPR)
  * tumor TACSTD2 ANTI-correlates with CD8 / NK cytotoxic signatures

For every usable bulk ICI lung cohort we:
  1. score CD8 / NK / cytotoxic and an epithelial (tumor-content / purity) signature,
  2. compute tumor-intrinsic target = z(target) - epithelial_score,
  3. correlate target vs immune signatures with raw AND partial (purity-corrected)
     Spearman controlling for epithelial content,
  4. run a directional Mann-Whitney (NR>R) on raw and tumor-intrinsic target.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from geo_utils import parse_series_matrix
import signatures as S

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

DATA = "/tmp/ici_bulk_data"
RES = "/workspace/results/fable_ici_bulk"
os.makedirs(os.path.join(RES, "tables"), exist_ok=True)
os.makedirs(os.path.join(RES, "figures"), exist_ok=True)

corr_rows = []
dir_rows = []
survival_rows = []


# ---------------------------------------------------------------------------
# Loaders: return (expr_log_symbol_df, labels_dict)
# ---------------------------------------------------------------------------
def load_gse126044():
    gse = "GSE126044"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    resp = {}
    for s in meta["samples"]:
        col = s["title"].replace("RNA-seq_", "")
        resp[col] = s.get("patient response")
    counts = pd.read_csv(f"{DATA}/{gse}_counts.txt.gz", sep="\t", index_col=0)
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1e6
    expr = np.log2(cpm + 1)
    return gse, expr, {"type": "binary", "resp": resp,
                       "responder": "responder", "nonresponder": "non-responder"}


def load_gse166449():
    gse = "GSE166449"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    resp = {}
    for s in meta["samples"]:
        col = s.get("description")
        t = s["title"].lower()
        resp[col] = "non-responder" if "nonresponder" in t else (
            "responder" if "responder" in t else None)
    tpm = pd.read_csv(f"{DATA}/{gse}_TPM.txt.gz", sep="\t", index_col=0)
    expr = np.log2(tpm + 1)
    return gse, expr, {"type": "binary", "resp": resp,
                       "responder": "responder", "nonresponder": "non-responder"}


def load_gse207422():
    gse = "GSE207422"
    md = pd.read_excel(f"{DATA}/{gse}_metadata.xlsx").dropna(
        subset=["Sample", "Pathologic Response"])
    md = md[md["Pathologic Response"].astype(str).str.contains("MPR")]
    resp = {}
    for _, r in md.iterrows():
        v = str(r["Pathologic Response"])
        resp[r["Sample"]] = "non-responder" if v.startswith("NMPR") else "responder"
    expr = pd.read_csv(f"{DATA}/{gse}_log2TPM.txt.gz", sep="\t", index_col=0)
    return gse, expr, {"type": "binary", "resp": resp,
                       "responder": "responder", "nonresponder": "non-responder"}


def load_gse253564():
    # Neoadjuvant durvalumab +/- radiation NSCLC, PRE-treatment resected tumors, FPKM.
    # GEO provides no responder label (outcome is in the paper), so correlation-only.
    gse = "GSE253564"
    df = pd.read_csv(f"{DATA}/{gse}_pre_FPKM.txt.gz", sep="\t", index_col=0)
    if "Entrez.ID" in df.columns:
        df = df.drop(columns=["Entrez.ID"])
    df = df[~df.index.duplicated(keep="first")]
    expr = np.log2(df.clip(lower=0) + 1)
    return gse, expr, {"type": "corr_only"}


def load_gse248378():
    # Neoadjuvant durvalumab +/- radiation NSCLC, POST-treatment resected tumors, FPKM.
    gse = "GSE248378"
    df = pd.read_csv(f"{DATA}/{gse}_post_FPKM.txt.gz", sep="\t", index_col=0)
    if "Entrez.ID" in df.columns:
        df = df.drop(columns=["Entrez.ID"])
    df = df[~df.index.duplicated(keep="first")]
    expr = np.log2(df.clip(lower=0) + 1)
    return gse, expr, {"type": "corr_only"}


def load_gse135222():
    gse = "GSE135222"
    meta = parse_series_matrix(f"{DATA}/{gse}_series_matrix.txt.gz")
    surv = {}
    for s in meta["samples"]:
        col = s["title"].replace(" ", "")
        e = s.get("progression-free survival (pfs)")
        t = s.get("pfs.time")
        if e is not None and t is not None:
            surv[col] = (int(e), float(t))
    tpm = pd.read_csv(f"{DATA}/{gse}_exp.tsv.gz", sep="\t", index_col=0)
    tpm.index = [i.split(".")[0] for i in tpm.index]
    tpm = tpm[tpm.index.isin(S.ENSG2SYMBOL)]
    tpm.index = [S.ENSG2SYMBOL[i] for i in tpm.index]
    expr = np.log2(tpm + 1)
    return gse, expr, {"type": "survival", "surv": surv}


# ---------------------------------------------------------------------------
def analyze_cohort(gse, expr, labels):
    # signature scores
    scores = {}
    found = {}
    for name, genes in S.SIGNATURES.items():
        sc, fg = S.signature_score(expr, genes)
        scores[name] = sc
        found[name] = fg
    epi = scores["Epithelial"]
    leuk = scores["Leukocyte"]

    for gene in S.TARGET_GENES:
        if gene not in expr.index:
            continue
        target = expr.loc[gene]
        tz = pd.Series(S.zscore(target.values), index=target.index)
        intrinsic = tz - epi  # tumor-intrinsic (purity-adjusted) target

        # correlations with immune signatures
        for sig in ["CD8_Tcell", "NK_cell", "Cytotoxic", "CD8_NK"]:
            sv = scores[sig]
            common = target.index
            r_raw, p_raw, n = S.spearman(target.reindex(common).values,
                                         sv.reindex(common).values)
            r_par, p_par, _ = S.partial_spearman(
                target.reindex(common).values, sv.reindex(common).values,
                [epi.reindex(common).values])
            r_leuk, p_leuk, _ = S.partial_spearman(
                target.reindex(common).values, sv.reindex(common).values,
                [leuk.reindex(common).values])
            corr_rows.append({
                "dataset": gse, "gene": gene, "immune_signature": sig,
                "n": n, "spearman_raw": r_raw, "p_raw": p_raw,
                "partial_epi": r_par, "p_partial_epi": p_par,
                "partial_leuk": r_leuk, "p_partial_leuk": p_leuk,
                "epithelial_genes_found": len(found["Epithelial"]),
                "signature_genes_found": len(found[sig]),
            })

        # directional response test (binary cohorts)
        if labels["type"] == "binary":
            resp = labels["resp"]
            r_lab, n_lab = labels["responder"], labels["nonresponder"]
            nr = [c for c in target.index if resp.get(c) == n_lab]
            r = [c for c in target.index if resp.get(c) == r_lab]
            for measure, series in [("raw_log_expr", target),
                                    ("tumor_intrinsic", intrinsic)]:
                res = S.mannwhitney_directional(series[nr].values,
                                                series[r].values)
                res.update({"dataset": gse, "gene": gene, "measure": measure,
                            "hypothesis": "NR > R"})
                dir_rows.append(res)
            _scatter(gse, gene, target, scores["CD8_NK"], resp,
                     r_lab, n_lab)

        elif labels["type"] == "survival":
            surv = labels["surv"]
            cols = [c for c in target.index if c in surv]
            df = pd.DataFrame({
                "time": [surv[c][1] for c in cols],
                "event": [surv[c][0] for c in cols],
                "raw": [target[c] for c in cols],
                "intrinsic": [intrinsic[c] for c in cols],
            })
            for measure in ["raw", "intrinsic"]:
                cph = CoxPHFitter()
                cph.fit(df[["time", "event", measure]], "time", "event")
                survival_rows.append({
                    "dataset": gse, "gene": gene, "measure": measure,
                    "n": len(df), "n_events": int(df["event"].sum()),
                    "cox_HR": float(np.exp(cph.params_[measure])),
                    "cox_CI_low": float(np.exp(cph.confidence_intervals_.loc[measure].iloc[0])),
                    "cox_CI_high": float(np.exp(cph.confidence_intervals_.loc[measure].iloc[1])),
                    "cox_p": float(cph.summary.loc[measure, "p"]),
                    "note": "HR>1 => higher expr worse PFS",
                })
            _scatter(gse, gene, target, scores["CD8_NK"], None, None, None)

        else:  # corr_only cohorts (durvalumab NSCLC, no GEO response label)
            _scatter(gse, gene, target, scores["CD8_NK"], None, None, None)


def _scatter(gse, gene, target, cd8nk, resp, r_lab, n_lab):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    x = cd8nk.reindex(target.index).values
    y = target.values
    if resp is not None:
        for lab, col in [(r_lab, "#2e8b57"), (n_lab, "#d1495b")]:
            idx = [i for i, c in enumerate(target.index) if resp.get(c) == lab]
            ax.scatter(np.asarray(x)[idx], np.asarray(y)[idx], s=30, color=col,
                       label=lab, alpha=0.85, edgecolor="k", linewidth=0.3)
        ax.legend(fontsize=8, title="response")
    else:
        ax.scatter(x, y, s=30, color="#555", alpha=0.8, edgecolor="k", linewidth=0.3)
    r_raw, p_raw, n = S.spearman(x, y)
    ax.set_xlabel("CD8/NK cytotoxic signature (mean z)")
    ax.set_ylabel(f"{gene} expression (log)")
    ax.set_title(f"{gse}: {gene} vs CD8/NK\nraw Spearman rho={r_raw:.2f}, p={p_raw:.3f} (n={n})")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "figures", f"{gse}_{gene}_vs_CD8NK.png"), dpi=140)
    plt.close(fig)


def main():
    for loader in [load_gse126044, load_gse166449, load_gse207422, load_gse135222,
                   load_gse253564, load_gse248378]:
        gse, expr, labels = loader()
        analyze_cohort(gse, expr, labels)

    cdf = pd.DataFrame(corr_rows)
    ddf = pd.DataFrame(dir_rows)
    sdf = pd.DataFrame(survival_rows)
    cdf.to_csv(f"{RES}/tables/purity_corrected_correlations.csv", index=False)
    ddf.to_csv(f"{RES}/tables/purity_corrected_directional.csv", index=False)
    sdf.to_csv(f"{RES}/tables/purity_corrected_survival.csv", index=False)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print("=== Purity-corrected correlations (target vs immune) ===")
    print(cdf.to_string(index=False))
    print("\n=== Directional response (NR > R), raw vs tumor-intrinsic ===")
    print(ddf[["dataset", "gene", "measure", "n_high", "n_low", "median_high",
               "median_low", "auc_high_gt_low", "p_one_sided_greater",
               "p_two_sided"]].to_string(index=False))
    print("\n=== Survival (GSE135222 PFS), raw vs tumor-intrinsic ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
