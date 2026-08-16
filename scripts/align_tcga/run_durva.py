"""
USER-ALIGN: public durvalumab NSCLC RNA hunt + same Spearman / purity-partial
method on every open matrix we could actually download.

Cohorts (GEO, open processed matrices):
  * GSE253564 — Altorki et al. Nat Commun 2024, *pre-treatment* FPKM (n=32).
  * GSE248378 — same trial, *post-treatment / resected* FPKM (n=29).

Neither matrix ships a published purity column. We compute ESTIMATE scores
with a Python port of Yoshihara et al. 2013 (validated against the official
sample_estimate.gct) and use `purity_proxy = -ESTIMATEScore` as the covariate.

No public matrix was found that reports the user's cited ρ = -0.65 / -0.46;
this script records what *is* computable from the open files.
"""
from __future__ import annotations

import gzip
import os
import pandas as pd

import common as C
from estimate_score import estimate_scores

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "results", "align_tcga", "data")
OUT = os.path.join(HERE, "..", "..", "results", "align_tcga")

COHORTS = {
    "GSE253564_pre": {
        "file": "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
        "label": "GSE253564_durvalumab_NSCLC_pretreatment",
    },
    "GSE248378_post": {
        "file": "GSE248378_Durva_Post_FPKMs.txt.gz",
        "label": "GSE248378_durvalumab_NSCLC_posttreatment",
    },
}


def load_fpkm(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    # drop annotation columns (e.g. numeric Entrez.ID in GSE253564)
    drop = [c for c in df.columns
            if str(c).lower() in {"entrez.id", "entrez", "entrezid",
                                  "gene_id", "description"}]
    df = df.drop(columns=drop, errors="ignore")
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    # log2(FPKM+1) — Spearman is rank-based so this does not change single-gene
    # ρ, but it stabilizes z-scored signature scores.
    return np_log2p1(df)


def np_log2p1(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        __import__("numpy").log2(df.astype(float).clip(lower=0) + 1.0),
        index=df.index, columns=df.columns,
    )


def compute_features(expr: pd.DataFrame) -> pd.DataFrame:
    if C.TARGET_GENE not in expr.index:
        raise SystemExit(f"{C.TARGET_GENE} missing in {expr.shape}")
    feats = {C.TARGET_GENE: expr.loc[C.TARGET_GENE].astype(float)}
    for name, genes in C.SIGNATURES.items():
        score, _ = C.signature_score(expr, genes)
        feats[name] = score
    avail = set(expr.index)
    for g in C.SINGLE_IMMUNE_GENES:
        r = C.resolve_symbol(g, avail)
        if r is not None:
            feats["gene_" + g] = expr.loc[r].astype(float)
    for g in C.TJ_GENES:
        r = C.resolve_symbol(g, avail)
        if r is not None:
            feats["TJ_" + g] = expr.loc[r].astype(float)
    return pd.DataFrame(feats)


def run_one(tag: str, expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    feat = compute_features(expr)
    est = estimate_scores(expr)
    # align
    common = feat.index.intersection(est.index)
    feat = feat.loc[common]
    est = est.loc[common]
    x = feat[C.TARGET_GENE].to_numpy(float)
    z = est["purity_proxy"].to_numpy(float)
    rows = []
    for col in [c for c in feat.columns if c != C.TARGET_GENE]:
        y = feat[col].to_numpy(float)
        r, p, n = C.spearman(x, y)
        pr, pp, pn = C.partial_spearman(x, y, z)
        rows.append({
            "cohort": tag, "feature": col,
            "class": "TJ" if col.startswith("TJ_") else "immune",
            "spearman_rho": r, "spearman_p": p, "n": n,
            "partial_rho_ESTIMATE": pr, "partial_p_ESTIMATE": pp,
            "partial_n_ESTIMATE": pn,
        })
    res = pd.DataFrame(rows)
    imm = res["class"] == "immune"
    res["spearman_q"] = float("nan")
    res["partial_q_ESTIMATE"] = float("nan")
    idx = list(res.index[imm])
    for pcol, qcol in [("spearman_p", "spearman_q"),
                       ("partial_p_ESTIMATE", "partial_q_ESTIMATE")]:
        qs = C.bh_fdr(res.loc[idx, pcol].tolist())
        for i, qi in zip(idx, qs):
            res.at[i, qcol] = qi
    # TACSTD2 vs ESTIMATE score / purity proxy
    vs = []
    for name, series in [("ESTIMATEScore", est["ESTIMATEScore"]),
                         ("purity_proxy", est["purity_proxy"]),
                         ("affy_calibrated_purity", est["affy_calibrated_purity"])]:
        r, p, n = C.spearman(x, series.to_numpy(float))
        vs.append({"cohort": tag, "purity_feature": name,
                   "spearman_rho": r, "spearman_p": p, "n": n,
                   "n_common_genes": int(est["n_common_genes"].iloc[0]),
                   "n_stromal_overlap": int(est["n_stromal_overlap"].iloc[0]),
                   "n_immune_overlap": int(est["n_immune_overlap"].iloc[0])})
    return res, pd.DataFrame(vs)


def main() -> None:
    all_res, all_vs = [], []
    for key, meta in COHORTS.items():
        path = os.path.join(DATA, meta["file"])
        if not os.path.exists(path):
            print(f"[durva] MISSING {path} — skip {key}")
            continue
        expr = load_fpkm(path)
        print(f"[durva] {key}: genes={expr.shape[0]} samples={expr.shape[1]}")
        res, vs = run_one(meta["label"], expr)
        all_res.append(res)
        all_vs.append(vs)
        sigs = list(C.SIGNATURES)
        print(res[res.feature.isin(sigs)][
            ["feature", "spearman_rho", "spearman_p",
             "partial_rho_ESTIMATE", "partial_p_ESTIMATE", "n"]
        ].to_string(index=False))

    if not all_res:
        raise SystemExit("no durvalumab matrices found")
    res = pd.concat(all_res, ignore_index=True)
    vs = pd.concat(all_vs, ignore_index=True)
    for c in res.columns:
        if c.endswith("_rho"):
            res[c] = res[c].round(4)
    for c in vs.columns:
        if c.endswith("_rho"):
            vs[c] = vs[c].round(4)
    res.to_csv(os.path.join(OUT, "durva_tacstd2_correlations.csv"), index=False)
    vs.to_csv(os.path.join(OUT, "durva_tacstd2_vs_estimate.csv"), index=False)
    print("[durva] wrote durva_tacstd2_correlations.csv and "
          "durva_tacstd2_vs_estimate.csv")


if __name__ == "__main__":
    main()
