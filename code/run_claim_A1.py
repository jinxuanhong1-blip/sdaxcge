#!/usr/bin/env python3
"""Claim A1 replication.

TACSTD2 (TROP2) vs immune / cytotoxic / exhaustion signatures in TCGA + OncoSG
NSCLC, using purity-adjusted (partial) Spearman correlation.

Public data only:
  - Expression + OncoSG purity: cBioPortal REST API
  - TCGA tumor purity: TCGA PanCanAtlas ABSOLUTE supplemental table (NCI GDC)

Outputs -> results/claim_A1/
"""
import csv
import io
import json
import os
import time
from math import sqrt, log

import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLAIM_DIR = os.path.join(ROOT, "claims", "A1")
DATA_DIR = os.path.join(ROOT, "data")
OUT_DIR = os.path.join(ROOT, "results", "claim_A1")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(CLAIM_DIR, "config.json")) as f:
    CFG = json.load(f)
with open(os.path.join(CLAIM_DIR, "signatures.json")) as f:
    SIG = json.load(f)

API = CFG["cbioportal_api"]
PROFILE_SUFFIX = CFG["expression_profile_suffix"]
TARGET = SIG["target_gene"]


def _get(url, **kw):
    last = None
    for i in range(4):
        try:
            r = requests.get(url, timeout=120, **kw)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa
            last = e
            time.sleep(2 ** i)
    raise last


def _post(url, **kw):
    last = None
    for i in range(4):
        try:
            r = requests.post(url, timeout=180, **kw)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa
            last = e
            time.sleep(2 ** i)
    raise last


def all_signature_genes():
    genes = {TARGET}
    for s in SIG["signatures"].values():
        genes.update(s["genes"])
    genes.update(SIG["cytotoxic_rooney_CYT"]["genes"])
    genes.update(SIG["individual_gene_checks"])
    return sorted(genes)


def resolve_entrez(genes):
    r = _post(API + "/genes/fetch", json=genes,
              params={"geneIdType": "HUGO_GENE_SYMBOL"})
    return {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in r.json()}


def fetch_expression(study_id, sample_list, entrez_map):
    """Return DataFrame indexed by sampleId, columns = hugo symbols (z-scores)."""
    profile = study_id + PROFILE_SUFFIX
    entrez_ids = sorted(set(entrez_map.values()))
    payload = {"entrezGeneIds": entrez_ids, "sampleListId": sample_list}
    r = _post(API + f"/molecular-profiles/{profile}/molecular-data/fetch",
              params={"projection": "SUMMARY"}, json=payload)
    rows = r.json()
    ent2hugo = {v: k for k, v in entrez_map.items()}
    recs = {}
    for d in rows:
        v = d.get("value")
        if v is None:
            continue
        sid = d["sampleId"]
        recs.setdefault(sid, {})[ent2hugo[d["entrezGeneId"]]] = float(v)
    df = pd.DataFrame.from_dict(recs, orient="index")
    df.index.name = "sampleId"
    return df


def load_tcga_purity():
    path = os.path.join(DATA_DIR, "abs_purity.txt")
    if not os.path.exists(path):
        url = CFG["purity_sources"]["TCGA_ABSOLUTE"]["url"]
        _get(url)  # ensure reachable
        with open(path, "wb") as fh:
            fh.write(_get(url).content)
    pur = {}
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            bc = row["array"]           # e.g. TCGA-05-4384-01
            val = row["purity"]
            if val in ("", "NA", None):
                continue
            try:
                pur[bc[:15]] = float(val)
            except ValueError:
                continue
    return pur


def load_oncosg_purity(study_id):
    r = _get(API + f"/studies/{study_id}/clinical-data",
             params={"clinicalDataType": "SAMPLE", "attributeId": "PURITY",
                     "projection": "DETAILED"})
    pur = {}
    for d in r.json():
        try:
            pur[d["sampleId"]] = float(d["value"])
        except (ValueError, KeyError):
            continue
    return pur


def purity_for_cohort(cohort, sample_ids):
    src = cohort["purity_source"]
    if src == "TCGA_ABSOLUTE":
        tab = load_tcga_purity()
        return {sid: tab[sid[:15]] for sid in sample_ids if sid[:15] in tab}
    if src == "cbioportal_clinical_PURITY":
        tab = load_oncosg_purity(cohort["study_id"])
        return {sid: tab[sid] for sid in sample_ids if sid in tab}
    raise ValueError(src)


def zscore(series):
    s = series.astype(float)
    sd = s.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return s * 0.0
    return (s - s.mean()) / sd


def build_signature_scores(expr):
    """expr: samples x genes z-scores. Returns DataFrame of signature scores."""
    scores = pd.DataFrame(index=expr.index)
    for name, spec in SIG["signatures"].items():
        genes = [g for g in spec["genes"] if g in expr.columns]
        # values are already per-gene z-scores from cBioPortal; mean across genes
        scores[name] = expr[genes].mean(axis=1)
    return scores


def spearman_rho(x, y):
    m = x.notna() & y.notna()
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(x[m], y[m])
    return rho, p, int(m.sum())


def partial_spearman(x, y, z):
    """First-order partial Spearman correlation of x,y controlling for z."""
    m = x.notna() & y.notna() & z.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    xr, yr, zr = x[m], y[m], z[m]
    rxy = stats.spearmanr(xr, yr).statistic
    rxz = stats.spearmanr(xr, zr).statistic
    ryz = stats.spearmanr(yr, zr).statistic
    denom = sqrt(max((1 - rxz ** 2) * (1 - ryz ** 2), 1e-12))
    r = (rxy - rxz * ryz) / denom
    r = max(min(r, 0.999999), -0.999999)
    df = n - 3
    if df <= 0:
        return r, np.nan, n
    t = r * sqrt(df / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), df)
    return r, p, n


def fisher_meta(rhos, ns):
    """Fixed-effect meta-analysis of partial correlations via Fisher z.
    variance approx 1/(n-3)."""
    zs, ws = [], []
    for r, n in zip(rhos, ns):
        if r is None or np.isnan(r) or n is None or n <= 3:
            continue
        r = max(min(r, 0.999999), -0.999999)
        zs.append(0.5 * log((1 + r) / (1 - r)))
        ws.append(n - 3)
    if not zs:
        return np.nan, np.nan, 0
    zs = np.array(zs); ws = np.array(ws, dtype=float)
    zbar = np.sum(ws * zs) / np.sum(ws)
    se = sqrt(1.0 / np.sum(ws))
    zstat = zbar / se
    p = 2 * stats.norm.sf(abs(zstat))
    rho = (np.exp(2 * zbar) - 1) / (np.exp(2 * zbar) + 1)
    return rho, p, int(np.sum(ws)) + 3 * len(zs)


def main():
    genes = all_signature_genes()
    entrez_map = resolve_entrez(genes)
    print(f"Resolved {len(entrez_map)}/{len(genes)} genes")

    per_cohort = {}      # cohort_id -> dict(target series, scores df, purity series)
    cohort_meta = []
    for cohort in CFG["cohorts"]:
        cid = cohort["id"]
        expr = fetch_expression(cohort["study_id"], cohort["sample_list"], entrez_map)
        scores = build_signature_scores(expr)
        pur_map = purity_for_cohort(cohort, list(expr.index))
        purity = pd.Series(pur_map, name="purity").reindex(expr.index)
        target = expr[TARGET] if TARGET in expr.columns else pd.Series(index=expr.index, dtype=float)
        per_cohort[cid] = {
            "target": target, "scores": scores, "purity": purity,
            "expr": expr, "histology": cohort["histology"],
        }
        n_expr = int(target.notna().sum())
        n_pur = int((target.notna() & purity.notna()).sum())
        cohort_meta.append({"cohort": cid, "study_id": cohort["study_id"],
                            "histology": cohort["histology"],
                            "n_expression": n_expr, "n_expr_and_purity": n_pur})
        print(f"{cid}: n_expr={n_expr} n_expr+purity={n_pur}")

    sig_names = list(SIG["signatures"].keys())

    # ---- per-cohort correlation table ----
    rows = []
    for cid, d in per_cohort.items():
        for sig in sig_names:
            rho_u, p_u, n_u = spearman_rho(d["target"], d["scores"][sig])
            rho_p, p_p, n_p = partial_spearman(d["target"], d["scores"][sig], d["purity"])
            rows.append({
                "cohort": cid, "histology": d["histology"], "signature": sig,
                "spearman_rho": rho_u, "spearman_p": p_u, "n_unadj": n_u,
                "partial_rho_purityadj": rho_p, "partial_p_purityadj": p_p,
                "n_partial": n_p,
            })
    per_df = pd.DataFrame(rows)
    per_df.to_csv(os.path.join(OUT_DIR, "per_cohort_correlations.csv"), index=False)

    # ---- pooled (concatenated) analyses ----
    def concat(cids):
        t = pd.concat([per_cohort[c]["target"] for c in cids])
        p = pd.concat([per_cohort[c]["purity"] for c in cids])
        sc = {s: pd.concat([per_cohort[c]["scores"][s] for c in cids]) for s in sig_names}
        return t, p, sc

    pooled_rows = []
    for pool_name, cids in CFG["pooled_definitions"].items():
        t, p, sc = concat(cids)
        for sig in sig_names:
            rho_u, p_u, n_u = spearman_rho(t, sc[sig])
            rho_p, p_p, n_p = partial_spearman(t, sc[sig], p)
            pooled_rows.append({
                "pool": pool_name, "cohorts": "+".join(cids), "signature": sig,
                "spearman_rho": rho_u, "spearman_p": p_u, "n_unadj": n_u,
                "partial_rho_purityadj": rho_p, "partial_p_purityadj": p_p,
                "n_partial": n_p,
            })
    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df.to_csv(os.path.join(OUT_DIR, "pooled_correlations.csv"), index=False)

    # ---- meta-analysis across the 3 cohorts (partial rho) ----
    meta_rows = []
    for sig in sig_names:
        rr = per_df[per_df.signature == sig]
        rho_m, p_m, n_m = fisher_meta(rr["partial_rho_purityadj"].tolist(),
                                      rr["n_partial"].tolist())
        meta_rows.append({"signature": sig, "meta_partial_rho": rho_m,
                          "meta_p": p_m, "meta_total_n": int(rr["n_partial"].sum())})
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(os.path.join(OUT_DIR, "meta_analysis.csv"), index=False)

    # ---- individual gene checks (pooled TCGA+OncoSG, purity-adjusted) ----
    t, p, _ = concat(CFG["pooled_definitions"]["TCGA_plus_OncoSG"])
    gene_rows = []
    for g in SIG["individual_gene_checks"]:
        gexpr = pd.concat([per_cohort[c]["expr"][g] for c in CFG["pooled_definitions"]["TCGA_plus_OncoSG"]
                           if g in per_cohort[c]["expr"].columns])
        rho_u, p_u, n_u = spearman_rho(t, gexpr)
        rho_p, p_p, n_p = partial_spearman(t, gexpr, p)
        gene_rows.append({"gene": g, "spearman_rho": rho_u, "spearman_p": p_u,
                          "partial_rho_purityadj": rho_p, "partial_p_purityadj": p_p,
                          "n_partial": n_p})
    gene_df = pd.DataFrame(gene_rows)
    gene_df.to_csv(os.path.join(OUT_DIR, "individual_gene_checks.csv"), index=False)

    # ---- summary json ----
    tp = pooled_df[pooled_df.pool == "TCGA_plus_OncoSG"].set_index("signature")
    summary = {
        "target_gene": TARGET,
        "cohorts": cohort_meta,
        "primary_pool": "TCGA_plus_OncoSG",
        "primary_pool_n_partial": {s: int(tp.loc[s, "n_partial"]) for s in sig_names},
        "primary_pool_partial_rho": {s: float(tp.loc[s, "partial_rho_purityadj"]) for s in sig_names},
        "primary_pool_partial_p": {s: float(tp.loc[s, "partial_p_purityadj"]) for s in sig_names},
        "all_partial_negative_after_purity": bool(
            all(tp.loc[s, "partial_rho_purityadj"] < 0 for s in sig_names)),
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # save merged per-sample data for reproducibility / figures
    for cid, d in per_cohort.items():
        m = d["scores"].copy()
        m.insert(0, TARGET, d["target"])
        m["purity"] = d["purity"]
        m.to_csv(os.path.join(OUT_DIR, f"sample_data_{cid}.csv"))

    print("\n== PRIMARY POOL (TCGA+OncoSG), purity-adjusted partial Spearman ==")
    print(tp[["cohorts", "partial_rho_purityadj", "partial_p_purityadj", "n_partial"]])
    print("\nDone. Outputs in", OUT_DIR)


if __name__ == "__main__":
    main()
