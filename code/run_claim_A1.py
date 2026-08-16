#!/usr/bin/env python3
"""Claim A1 replication.

TACSTD2 (TROP2) vs immune / cytotoxic / exhaustion signatures in TCGA
(LUAD+LUSC; PanCanAtlas and Firehose) + OncoSG, using purity-adjusted
(partial) Spearman correlation under ESTIMATE, ABSOLUTE, CPE, and
ESTIMATE-or-ABSOLUTE.

Public data only:
  - Expression + OncoSG purity: cBioPortal REST API
  - TCGA ABSOLUTE: PanCanAtlas ABSOLUTE supplemental table (NCI GDC)
  - TCGA ESTIMATE / Aran ABSOLUTE / CPE: Aran et al. Nat Commun 2015 Supp Data 1

Outputs -> results/claim_A1/
"""
import csv
import json
import os
import time
from math import log, sqrt

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
CLAIMED_N = CFG["claimed_n"]
N_TOL = CFG["n_match_tolerance"]


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


def load_gdc_absolute():
    path = os.path.join(DATA_DIR, "abs_purity.txt")
    src = CFG["purity_sources"]["ABSOLUTE_GDC"]
    if not os.path.exists(path):
        with open(path, "wb") as fh:
            fh.write(_get(src["url"]).content)
    pur = {}
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            val = row["purity"]
            if val in ("", "NA", None):
                continue
            try:
                pur[row["array"][:15]] = float(val)
            except ValueError:
                continue
    return pd.Series(pur, name="ABSOLUTE_GDC")


def load_aran_table():
    path = os.path.join(DATA_DIR, "aran_purity.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            "data/aran_purity.tsv missing; parse Aran Nat Commun 2015 Supp Data 1")
    df = pd.read_csv(path, sep="\t")
    df["bc15"] = df["sample_id"].astype(str).str[:15]
    df["vial"] = df["sample_id"].astype(str).str[15:16]
    df = df.sort_values(["bc15", "vial"])
    out = {}
    for col in ["ESTIMATE", "ABSOLUTE", "CPE"]:
        sub = df[df[col].notna()].drop_duplicates("bc15", keep="first")
        out[col] = pd.Series(sub[col].values, index=sub["bc15"], name=col)
    return out


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
    return pd.Series(pur, name="OncoSG_clinical_PURITY")


def build_signature_scores(expr):
    scores = pd.DataFrame(index=expr.index)
    for name, spec in SIG["signatures"].items():
        genes = [g for g in spec["genes"] if g in expr.columns]
        scores[name] = expr[genes].mean(axis=1)
    return scores


def spearman_rho(x, y):
    m = x.notna() & y.notna()
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def partial_spearman(x, y, z):
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
        return float(r), np.nan, n
    t = r * sqrt(df / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), df)
    return float(r), float(p), n


def fisher_meta(rhos, ns):
    zs, ws = [], []
    for r, n in zip(rhos, ns):
        if r is None or np.isnan(r) or n is None or n <= 3:
            continue
        r = max(min(float(r), 0.999999), -0.999999)
        zs.append(0.5 * log((1 + r) / (1 - r)))
        ws.append(n - 3)
    if not zs:
        return np.nan, np.nan, 0
    zs = np.array(zs)
    ws = np.array(ws, dtype=float)
    zbar = np.sum(ws * zs) / np.sum(ws)
    se = sqrt(1.0 / np.sum(ws))
    p = 2 * stats.norm.sf(abs(zbar / se))
    rho = (np.exp(2 * zbar) - 1) / (np.exp(2 * zbar) + 1)
    return float(rho), float(p), int(np.sum(ws)) + 3 * len(zs)


def attach_purity(sample_ids, atlas, tables):
    """Return DataFrame of purity columns aligned to sample_ids."""
    idx = pd.Index(sample_ids, name="sampleId")
    bc = pd.Series([s[:15] if str(s).startswith("TCGA") else s for s in idx],
                   index=idx)
    out = pd.DataFrame(index=idx)
    if atlas == "OncoSG":
        out["OncoSG_clinical_PURITY"] = tables["oncosg"].reindex(idx)
        return out
    out["ESTIMATE_Aran"] = tables["aran"]["ESTIMATE"].reindex(bc.values).values
    out["ABSOLUTE_Aran"] = tables["aran"]["ABSOLUTE"].reindex(bc.values).values
    out["CPE_Aran"] = tables["aran"]["CPE"].reindex(bc.values).values
    out["ABSOLUTE_GDC"] = tables["gdc"].reindex(bc.values).values
    out["ESTIMATE_or_ABSOLUTE"] = out["ESTIMATE_Aran"].fillna(out["ABSOLUTE_GDC"])
    return out


def method_purity_series(purity_df, atlas, method):
    spec = CFG["purity_methods"][method]
    if atlas == "OncoSG":
        src = spec["oncosg"]
        if src is None:
            return pd.Series(np.nan, index=purity_df.index, name=method)
        return purity_df[src].rename(method)
    tcga = spec["tcga"]
    if isinstance(tcga, list):
        s = purity_df[tcga[0]]
        for c in tcga[1:]:
            s = s.fillna(purity_df[c])
        return s.rename(method)
    return purity_df[tcga].rename(method)


def main():
    genes = all_signature_genes()
    entrez_map = resolve_entrez(genes)
    print(f"Resolved {len(entrez_map)}/{len(genes)} genes")

    tables = {
        "gdc": load_gdc_absolute(),
        "aran": load_aran_table(),
        "oncosg": load_oncosg_purity("luad_oncosg_2020"),
    }
    print(f"GDC ABSOLUTE barcodes: {tables['gdc'].shape[0]}")
    print(f"Aran ESTIMATE barcodes: {tables['aran']['ESTIMATE'].shape[0]}")
    print(f"OncoSG clinical PURITY: {tables['oncosg'].shape[0]}")

    per_cohort = {}
    cohort_meta = []
    for cohort in CFG["cohorts"]:
        cid = cohort["id"]
        expr = fetch_expression(cohort["study_id"], cohort["sample_list"], entrez_map)
        scores = build_signature_scores(expr)
        purity_df = attach_purity(list(expr.index), cohort["atlas"], tables)
        target = (expr[TARGET] if TARGET in expr.columns
                  else pd.Series(index=expr.index, dtype=float))
        per_cohort[cid] = {
            "target": target, "scores": scores, "purity_df": purity_df,
            "expr": expr, "histology": cohort["histology"],
            "atlas": cohort["atlas"],
        }
        n_expr = int(target.notna().sum())
        rec = {"cohort": cid, "study_id": cohort["study_id"],
               "histology": cohort["histology"], "atlas": cohort["atlas"],
               "n_expression": n_expr}
        for col in purity_df.columns:
            rec[f"n_{col}"] = int((target.notna() & purity_df[col].notna()).sum())
        cohort_meta.append(rec)
        print(f"{cid}: n_expr={n_expr} " +
              " ".join(f"{c}={rec.get('n_'+c,'')}" for c in purity_df.columns))

    sig_names = list(SIG["signatures"].keys())
    methods = list(CFG["purity_methods"].keys())

    # ---- per-cohort x method ----
    rows = []
    for cid, d in per_cohort.items():
        for method in methods:
            pur = method_purity_series(d["purity_df"], d["atlas"], method)
            for sig in sig_names:
                rho_u, p_u, n_u = spearman_rho(d["target"], d["scores"][sig])
                rho_p, p_p, n_p = partial_spearman(d["target"], d["scores"][sig], pur)
                rows.append({
                    "cohort": cid, "atlas": d["atlas"], "histology": d["histology"],
                    "purity_method": method, "signature": sig,
                    "spearman_rho": rho_u, "spearman_p": p_u, "n_unadj": n_u,
                    "partial_rho": rho_p, "partial_p": p_p, "n_partial": n_p,
                })
    per_df = pd.DataFrame(rows)
    per_df.to_csv(os.path.join(OUT_DIR, "per_cohort_correlations.csv"), index=False)

    # ---- pooled x method ----
    def concat(cids, method):
        t = pd.concat([per_cohort[c]["target"] for c in cids])
        p = pd.concat([method_purity_series(per_cohort[c]["purity_df"],
                                            per_cohort[c]["atlas"], method)
                       for c in cids])
        sc = {s: pd.concat([per_cohort[c]["scores"][s] for c in cids])
              for s in sig_names}
        return t, p, sc

    pooled_rows = []
    for pool_name, cids in CFG["pooled_definitions"].items():
        for method in methods:
            t, p, sc = concat(cids, method)
            for sig in sig_names:
                rho_u, p_u, n_u = spearman_rho(t, sc[sig])
                rho_p, p_p, n_p = partial_spearman(t, sc[sig], p)
                pooled_rows.append({
                    "pool": pool_name, "cohorts": "+".join(cids),
                    "purity_method": method, "signature": sig,
                    "spearman_rho": rho_u, "spearman_p": p_u, "n_unadj": n_u,
                    "partial_rho": rho_p, "partial_p": p_p, "n_partial": n_p,
                    "n_vs_claimed": n_p - CLAIMED_N,
                    "n_within_tolerance": bool(abs(n_p - CLAIMED_N) <= N_TOL),
                    "sign_match": bool(rho_p < 0) if not np.isnan(rho_p) else False,
                })
    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df.to_csv(os.path.join(OUT_DIR, "pooled_correlations.csv"), index=False)

    # ---- meta-analysis: PanCan LUAD, LUSC, OncoSG (where purity exists) ----
    meta_rows = []
    meta_cohorts = ["TCGA_LUAD_PanCan", "TCGA_LUSC_PanCan", "OncoSG_LUAD"]
    for method in methods:
        for sig in sig_names:
            rr = per_df[(per_df.signature == sig) & (per_df.purity_method == method)
                        & (per_df.cohort.isin(meta_cohorts)) & (per_df.n_partial >= 6)]
            rho_m, p_m, _ = fisher_meta(rr["partial_rho"].tolist(),
                                        rr["n_partial"].tolist())
            meta_rows.append({
                "purity_method": method, "signature": sig,
                "meta_partial_rho": rho_m, "meta_p": p_m,
                "meta_total_n": int(rr["n_partial"].sum()),
                "n_cohorts": int(len(rr)),
            })
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(os.path.join(OUT_DIR, "meta_analysis.csv"), index=False)

    # ---- individual gene checks: PanCan+OncoSG, ESTIMATE_or_ABSOLUTE ----
    t, p, _ = concat(CFG["pooled_definitions"]["PanCan_plus_OncoSG"],
                     "ESTIMATE_or_ABSOLUTE")
    gene_rows = []
    for g in SIG["individual_gene_checks"]:
        gexpr = pd.concat([
            per_cohort[c]["expr"][g]
            for c in CFG["pooled_definitions"]["PanCan_plus_OncoSG"]
            if g in per_cohort[c]["expr"].columns
        ])
        rho_u, p_u, n_u = spearman_rho(t, gexpr)
        rho_p, p_p, n_p = partial_spearman(t, gexpr, p)
        gene_rows.append({
            "gene": g, "purity_method": "ESTIMATE_or_ABSOLUTE",
            "spearman_rho": rho_u, "spearman_p": p_u,
            "partial_rho": rho_p, "partial_p": p_p, "n_partial": n_p,
        })
    gene_df = pd.DataFrame(gene_rows)
    gene_df.to_csv(os.path.join(OUT_DIR, "individual_gene_checks.csv"), index=False)

    # ---- match / mismatch vs user claim ----
    # Claim: still negative after purity, n≈1031, LUAD+LUSC+OncoSG
    primary_pools = ["PanCan_plus_OncoSG", "Firehose_plus_OncoSG",
                     "PanCan_LUAD_LUSC", "Firehose_LUAD_LUSC"]
    match_rows = []
    for _, r in pooled_df.iterrows():
        if r["pool"] not in primary_pools:
            continue
        n_ok = bool(r["n_within_tolerance"])
        sign_ok = bool(r["sign_match"])
        sig_ok = (not np.isnan(r["partial_p"])) and (r["partial_p"] < 0.05)
        oncosg_dropped = (
            "OncoSG" in r["pool"] and r["purity_method"] == "ESTIMATE"
        )
        if oncosg_dropped and sign_ok:
            verdict = "MATCH_sign_OncoSG_dropped_no_public_ESTIMATE"
        elif n_ok and sign_ok:
            verdict = "MATCH"
        elif sign_ok and not n_ok:
            verdict = "PARTIAL_MATCH_sign_ok_n_mismatch"
        elif n_ok and not sign_ok:
            verdict = "MISMATCH_n_ok_sign_fail"
        else:
            verdict = "MISMATCH"
        match_rows.append({
            "pool": r["pool"], "purity_method": r["purity_method"],
            "signature": r["signature"],
            "n": int(r["n_partial"]), "claimed_n": CLAIMED_N,
            "n_delta": int(r["n_vs_claimed"]),
            "partial_rho": r["partial_rho"], "partial_p": r["partial_p"],
            "sign_negative": sign_ok, "p_lt_0.05": bool(sig_ok),
            "n_within_pm50_of_1031": n_ok,
            "verdict": verdict,
        })
    match_df = pd.DataFrame(match_rows)
    match_df.to_csv(os.path.join(OUT_DIR, "match_mismatch.csv"), index=False)

    # closest n to 1031
    closest = (pooled_df.assign(abs_delta=lambda d: (d.n_partial - CLAIMED_N).abs())
               .sort_values(["abs_delta", "pool"]))
    closest_rec = closest.iloc[0]

    # primary user-requested slice: LUAD+LUSC+OncoSG, ESTIMATE/ABSOLUTE
    primary = pooled_df[(pooled_df.pool == "PanCan_plus_OncoSG")
                        & (pooled_df.purity_method == "ESTIMATE_or_ABSOLUTE")]
    primary_abs = pooled_df[(pooled_df.pool == "PanCan_plus_OncoSG")
                            & (pooled_df.purity_method == "ABSOLUTE")]
    primary_est_tcga = pooled_df[(pooled_df.pool == "PanCan_LUAD_LUSC")
                                 & (pooled_df.purity_method == "ESTIMATE")]
    fh_est = pooled_df[(pooled_df.pool == "Firehose_LUAD_LUSC")
                       & (pooled_df.purity_method == "ESTIMATE")]

    summary = {
        "target_gene": TARGET,
        "claimed_n": CLAIMED_N,
        "cohorts": cohort_meta,
        "primary_request": "LUAD+LUSC+OncoSG, ESTIMATE/ABSOLUTE purity-partial Spearman",
        "primary_pool": "PanCan_plus_OncoSG",
        "primary_purity_method": "ESTIMATE_or_ABSOLUTE",
        "primary": {
            s: {
                "n": int(primary.set_index("signature").loc[s, "n_partial"]),
                "partial_rho": float(primary.set_index("signature").loc[s, "partial_rho"]),
                "partial_p": float(primary.set_index("signature").loc[s, "partial_p"]),
            } for s in sig_names
        },
        "absolute_only_PanCan_plus_OncoSG": {
            s: {
                "n": int(primary_abs.set_index("signature").loc[s, "n_partial"]),
                "partial_rho": float(primary_abs.set_index("signature").loc[s, "partial_rho"]),
                "partial_p": float(primary_abs.set_index("signature").loc[s, "partial_p"]),
            } for s in sig_names
        },
        "estimate_only_PanCan_LUAD_LUSC": {
            s: {
                "n": int(primary_est_tcga.set_index("signature").loc[s, "n_partial"]),
                "partial_rho": float(primary_est_tcga.set_index("signature").loc[s, "partial_rho"]),
                "partial_p": float(primary_est_tcga.set_index("signature").loc[s, "partial_p"]),
            } for s in sig_names
        },
        "estimate_only_Firehose_LUAD_LUSC": {
            s: {
                "n": int(fh_est.set_index("signature").loc[s, "n_partial"]),
                "partial_rho": float(fh_est.set_index("signature").loc[s, "partial_rho"]),
                "partial_p": float(fh_est.set_index("signature").loc[s, "partial_p"]),
            } for s in sig_names
        },
        "closest_n_to_1031": {
            "pool": closest_rec["pool"],
            "purity_method": closest_rec["purity_method"],
            "signature": closest_rec["signature"],
            "n": int(closest_rec["n_partial"]),
            "n_delta": int(closest_rec["n_vs_claimed"]),
        },
        "all_primary_partial_negative": bool(
            all(primary.set_index("signature").loc[s, "partial_rho"] < 0
                for s in sig_names)),
        "n_match_primary_vs_1031": bool(
            abs(int(primary.iloc[0]["n_partial"]) - CLAIMED_N) <= N_TOL),
        "verdict": (
            "MATCH_sign_MISMATCH_n"
            if all(primary.set_index("signature").loc[s, "partial_rho"] < 0
                   for s in sig_names)
            and abs(int(primary.iloc[0]["n_partial"]) - CLAIMED_N) > N_TOL
            else "see match_mismatch.csv"
        ),
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    for cid, d in per_cohort.items():
        m = d["scores"].copy()
        m.insert(0, TARGET, d["target"])
        m = m.join(d["purity_df"])
        m.to_csv(os.path.join(OUT_DIR, f"sample_data_{cid}.csv"))

    print("\n== PRIMARY: PanCan LUAD+LUSC+OncoSG, ESTIMATE_or_ABSOLUTE ==")
    print(primary[["signature", "partial_rho", "partial_p", "n_partial"]].to_string(index=False))
    print("\n== ABSOLUTE only, same pool ==")
    print(primary_abs[["signature", "partial_rho", "partial_p", "n_partial"]].to_string(index=False))
    print("\n== ESTIMATE only, PanCan LUAD+LUSC (no OncoSG ESTIMATE) ==")
    print(primary_est_tcga[["signature", "partial_rho", "partial_p", "n_partial"]].to_string(index=False))
    print("\n== ESTIMATE only, Firehose LUAD+LUSC (closest published ESTIMATE n) ==")
    print(fh_est[["signature", "partial_rho", "partial_p", "n_partial"]].to_string(index=False))
    print("\nClosest n to 1031:", summary["closest_n_to_1031"])
    print("Done.", OUT_DIR)


if __name__ == "__main__":
    main()
