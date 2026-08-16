#!/usr/bin/env python3
"""B5 analog — KIRC/RCC: CLDN4 vs ICI response on open public expression cohorts.

Pre-specified in results/w200/B5_KIRC/ANALYSIS_PLAN.md (committed before this
script was run against outcomes). Re-run:

    python3 scripts/w200/B5_KIRC/analyze.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

ROOT = Path("/workspace")
DATA = ROOT / "data"
OUT = ROOT / "results" / "w200" / "B5_KIRC"
OUT.mkdir(parents=True, exist_ok=True)

RESP_ORR = {"CR", "PR", "CRPR"}
PD_ORR = {"PD"}
EVAL_ORR = {"CR", "PR", "CRPR", "SD", "PD"}
REF_GENES = [
    "CLDN4",
    "TACSTD2",
    "CLDN3",
    "CLDN7",
    "EPCAM",
    "CDH1",
    "CD8A",
    "GZMB",
    "PDCD1",
    "CD274",
    "PBRM1",
    "ACTB",
    "CA9",
    "PAX8",
]


def zscore(x: pd.Series) -> pd.Series:
    s = float(x.std(ddof=0))
    if s == 0 or not np.isfinite(s):
        return x * np.nan
    return (x - x.mean()) / s


def logit_or(y: np.ndarray, z: np.ndarray, label: str) -> dict:
    """Logistic OR per +1 SD of z. y is 0/1."""
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(y) & np.isfinite(z)
    y, z = y[ok], z[ok]
    n = int(y.size)
    n1 = int(y.sum())
    n0 = n - n1
    if n1 < 2 or n0 < 2 or np.nanstd(z) == 0:
        return {
            "label": label,
            "n": n,
            "n_resp": n1,
            "n_non": n0,
            "or_per_sd": np.nan,
            "or_lo": np.nan,
            "or_hi": np.nan,
            "p": np.nan,
            "beta": np.nan,
            "se": np.nan,
            "converged": False,
        }
    X = sm.add_constant(z)
    try:
        fit = sm.Logit(y, X).fit(disp=False, maxiter=100)
        beta = float(fit.params[1])
        se = float(fit.bse[1])
        ci = fit.conf_int()[1]
        return {
            "label": label,
            "n": n,
            "n_resp": n1,
            "n_non": n0,
            "or_per_sd": float(np.exp(beta)),
            "or_lo": float(np.exp(ci[0])),
            "or_hi": float(np.exp(ci[1])),
            "p": float(fit.pvalues[1]),
            "beta": beta,
            "se": se,
            "converged": bool(fit.mle_retvals.get("converged", True)),
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "label": label,
            "n": n,
            "n_resp": n1,
            "n_non": n0,
            "or_per_sd": np.nan,
            "or_lo": np.nan,
            "or_hi": np.nan,
            "p": np.nan,
            "beta": np.nan,
            "se": np.nan,
            "converged": False,
            "error": str(exc),
        }


def mwu_auc(x_resp: np.ndarray, x_non: np.ndarray) -> dict:
    x_resp = np.asarray(x_resp, dtype=float)
    x_non = np.asarray(x_non, dtype=float)
    x_resp = x_resp[np.isfinite(x_resp)]
    x_non = x_non[np.isfinite(x_non)]
    if x_resp.size < 2 or x_non.size < 2:
        return {"U": np.nan, "p": np.nan, "auc": np.nan, "median_resp": np.nan, "median_non": np.nan}
    res = stats.mannwhitneyu(x_resp, x_non, alternative="two-sided")
    auc = float(res.statistic) / (x_resp.size * x_non.size)
    return {
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "auc": auc,
        "median_resp": float(np.median(x_resp)),
        "median_non": float(np.median(x_non)),
        "mean_resp": float(np.mean(x_resp)),
        "mean_non": float(np.mean(x_non)),
    }


def fisher_or(table: np.ndarray) -> dict:
    """table = [[a,b],[c,d]] rows = high/low, cols = resp/non."""
    table = np.asarray(table, dtype=int)
    if table.shape != (2, 2) or table.sum() == 0:
        return {"or": np.nan, "or_lo": np.nan, "or_hi": np.nan, "p": np.nan}
    oddsr, p = stats.fisher_exact(table, alternative="two-sided")
    try:
        t = Table2x2(table)
        ci = t.oddsratio_confint()
        or_mle = float(t.oddsratio)
        return {
            "or": or_mle,
            "or_lo": float(ci[0]),
            "or_hi": float(ci[1]),
            "p": float(p),
            "fisher_or": float(oddsr),
        }
    except Exception:  # noqa: BLE001
        return {"or": float(oddsr), "or_lo": np.nan, "or_hi": np.nan, "p": float(p)}


def cox_hr(time, event, z, label: str) -> dict:
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(time) & np.isfinite(event) & np.isfinite(z) & (time > 0)
    time, event, z = time[ok], event[ok], z[ok]
    n = int(time.size)
    n_evt = int(event.sum())
    if n < 10 or n_evt < 5 or np.nanstd(z) == 0:
        return {
            "label": label,
            "n": n,
            "n_event": n_evt,
            "hr_per_sd": np.nan,
            "hr_lo": np.nan,
            "hr_hi": np.nan,
            "p": np.nan,
        }
    try:
        fit = PHReg(time, z, status=event).fit(disp=0)
        beta = float(np.asarray(fit.params).ravel()[0])
        se = float(np.asarray(fit.bse).ravel()[0])
        zstat = beta / se
        p = float(2 * stats.norm.sf(abs(zstat)))
        return {
            "label": label,
            "n": n,
            "n_event": n_evt,
            "hr_per_sd": float(np.exp(beta)),
            "hr_lo": float(np.exp(beta - 1.96 * se)),
            "hr_hi": float(np.exp(beta + 1.96 * se)),
            "p": p,
            "beta": beta,
            "se": se,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "label": label,
            "n": n,
            "n_event": n_evt,
            "hr_per_sd": np.nan,
            "hr_lo": np.nan,
            "hr_hi": np.nan,
            "p": np.nan,
            "error": str(exc),
        }


def power_wald(n: int, p_event: float, or_true: float, alpha: float = 0.05) -> float:
    """Approximate Wald power for logistic OR per +1 SD (var(z)=1)."""
    info = n * p_event * (1.0 - p_event)
    if info <= 0:
        return float("nan")
    se = 1.0 / math.sqrt(info)
    mu = abs(math.log(or_true)) / se
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(stats.norm.sf(zcrit - mu) + stats.norm.cdf(-zcrit - mu))


def load_braun():
    clin = pd.read_excel(
        DATA / "braun2020" / "braun2020_supp.xlsx",
        sheet_name="S1_Clinical_and_Immune_Data",
        header=1,
    )
    clin = clin[clin["RNA_ID"].notna()].copy()
    clin["RNA_ID"] = clin["RNA_ID"].astype(str)
    expr = pd.read_csv(DATA / "braun2020" / "S4A_rna.tsv", sep="\t", skiprows=1, index_col=0)
    expr.columns = expr.columns.astype(str)
    missing = set(clin["RNA_ID"]) - set(expr.columns)
    if missing:
        raise SystemExit(f"RNA_ID not in expression matrix: {missing}")
    genes = {g: expr.loc[g, clin["RNA_ID"]].astype(float).to_numpy() for g in REF_GENES if g in expr.index}
    clin = clin.reset_index(drop=True)
    for g, vals in genes.items():
        clin[g] = vals
    return clin, expr


def infer_event_coding(pfs, cnsr, orr) -> str:
    """Return 'cnsr1_is_event' or 'cnsr1_is_censored' from PD vs long-PFS pattern."""
    df = pd.DataFrame({"pfs": pfs, "cnsr": cnsr, "orr": orr})
    df = df[df["pfs"].notna() & df["cnsr"].notna()]
    pd_rate = float((df.loc[df["orr"] == "PD", "cnsr"] == 1).mean()) if (df["orr"] == "PD").any() else np.nan
    # If almost all PD have CNSR=1, CNSR=1 is the event flag (BMS-style).
    if pd_rate >= 0.8:
        return "cnsr1_is_event"
    if pd_rate <= 0.2:
        return "cnsr1_is_censored"
    # fallback: more events than censored among PD
    return "cnsr1_is_event" if pd_rate >= 0.5 else "cnsr1_is_censored"


def main() -> None:
    clin, expr = load_braun()

    # --- primary population ---
    nivo = clin[clin["Arm"] == "NIVOLUMAB"].copy()
    prim = nivo[nivo["ORR"].isin(RESP_ORR | PD_ORR)].copy()
    prim["y"] = prim["ORR"].isin(RESP_ORR).astype(int)
    prim["z_CLDN4"] = zscore(prim["CLDN4"])
    primary = logit_or(prim["y"].to_numpy(), prim["z_CLDN4"].to_numpy(), "primary_nivo_CRPR_vs_PD")
    primary_mwu = mwu_auc(
        prim.loc[prim["y"] == 1, "CLDN4"].to_numpy(),
        prim.loc[prim["y"] == 0, "CLDN4"].to_numpy(),
    )

    # median / tertile splits
    med = float(prim["CLDN4"].median())
    prim["high"] = (prim["CLDN4"] >= med).astype(int)
    # rows: high, low; cols: resp, non
    tab_med = np.array(
        [
            [
                int(((prim["high"] == 1) & (prim["y"] == 1)).sum()),
                int(((prim["high"] == 1) & (prim["y"] == 0)).sum()),
            ],
            [
                int(((prim["high"] == 0) & (prim["y"] == 1)).sum()),
                int(((prim["high"] == 0) & (prim["y"] == 0)).sum()),
            ],
        ]
    )
    med_or = fisher_or(tab_med)
    q1, q2 = prim["CLDN4"].quantile([1 / 3, 2 / 3])
    tert_hi = prim[prim["CLDN4"] >= q2]
    tert_lo = prim[prim["CLDN4"] <= q1]
    tab_ter = np.array(
        [
            [int((tert_hi["y"] == 1).sum()), int((tert_hi["y"] == 0).sum())],
            [int((tert_lo["y"] == 1).sum()), int((tert_lo["y"] == 0).sum())],
        ]
    )
    ter_or = fisher_or(tab_ter)

    # alternative endpoints
    eval_nivo = nivo[nivo["ORR"].isin(EVAL_ORR)].copy()
    eval_nivo["y"] = eval_nivo["ORR"].isin(RESP_ORR).astype(int)
    eval_nivo["z_CLDN4"] = zscore(eval_nivo["CLDN4"])
    alt_orr = logit_or(eval_nivo["y"].to_numpy(), eval_nivo["z_CLDN4"].to_numpy(), "nivo_CRPR_vs_SDPD")

    cb = nivo[nivo["Benefit"].isin(["CB", "NCB"])].copy()
    cb["y"] = (cb["Benefit"] == "CB").astype(int)
    cb["z_CLDN4"] = zscore(cb["CLDN4"])
    alt_cb = logit_or(cb["y"].to_numpy(), cb["z_CLDN4"].to_numpy(), "nivo_CB_vs_NCB")

    # multivariable
    mv = prim.copy()
    mv["purity"] = pd.to_numeric(mv["Purity"], errors="coerce")
    mv_ok = mv.dropna(subset=["z_CLDN4", "y", "purity", "Cohort"])
    mv_res = {"n": int(len(mv_ok)), "note": "adjusted for cohort + purity"}
    if len(mv_ok) >= 30 and mv_ok["y"].nunique() == 2:
        X = sm.add_constant(
            pd.concat(
                [mv_ok[["z_CLDN4", "purity"]], pd.get_dummies(mv_ok["Cohort"], drop_first=True, dtype=float)],
                axis=1,
            )
        )
        try:
            fit = sm.Logit(mv_ok["y"].astype(float), X).fit(disp=False, maxiter=100)
            mv_res.update(
                {
                    "or_per_sd": float(np.exp(fit.params["z_CLDN4"])),
                    "or_lo": float(np.exp(fit.conf_int().loc["z_CLDN4", 0])),
                    "or_hi": float(np.exp(fit.conf_int().loc["z_CLDN4", 1])),
                    "p": float(fit.pvalues["z_CLDN4"]),
                    "n_with_purity": int(len(mv_ok)),
                    "n_purity_missing_in_primary": int(prim["Purity"].isna().sum())
                    if "Purity" in prim
                    else None,
                }
            )
        except Exception as exc:  # noqa: BLE001
            mv_res["error"] = str(exc)

    mv2 = prim.copy()
    mv2["MSKCC"] = mv2["MSKCC"].astype(str)
    mv2_ok = mv2[mv2["MSKCC"].isin(["FAVORABLE", "INTERMEDIATE", "POOR"])].copy()
    mv2_res = {"n": int(len(mv2_ok)), "note": "adjusted for MSKCC"}
    if len(mv2_ok) >= 30:
        X = sm.add_constant(
            pd.concat(
                [zscore(mv2_ok["CLDN4"]).rename("z_CLDN4"), pd.get_dummies(mv2_ok["MSKCC"], drop_first=True, dtype=float)],
                axis=1,
            )
        )
        try:
            fit = sm.Logit(mv2_ok["y"].astype(float), X).fit(disp=False, maxiter=100)
            mv2_res.update(
                {
                    "or_per_sd": float(np.exp(fit.params["z_CLDN4"])),
                    "or_lo": float(np.exp(fit.conf_int().loc["z_CLDN4", 0])),
                    "or_hi": float(np.exp(fit.conf_int().loc["z_CLDN4", 1])),
                    "p": float(fit.pvalues["z_CLDN4"]),
                }
            )
        except Exception as exc:  # noqa: BLE001
            mv2_res["error"] = str(exc)

    # everolimus control
    eve = clin[clin["Arm"] == "EVEROLIMUS"].copy()
    eve_p = eve[eve["ORR"].isin(RESP_ORR | PD_ORR)].copy()
    eve_p["y"] = eve_p["ORR"].isin(RESP_ORR).astype(int)
    eve_p["z_CLDN4"] = zscore(eve_p["CLDN4"])
    eve_logit = logit_or(eve_p["y"].to_numpy(), eve_p["z_CLDN4"].to_numpy(), "everolimus_CRPR_vs_PD")

    # survival coding
    coding = infer_event_coding(nivo["PFS"], nivo["PFS_CNSR"], nivo["ORR"])
    if coding == "cnsr1_is_event":
        nivo_event_pfs = nivo["PFS_CNSR"].astype(float)
        nivo_event_os = nivo["OS_CNSR"].astype(float)
        eve_event_pfs = eve["PFS_CNSR"].astype(float)
        eve_event_os = eve["OS_CNSR"].astype(float)
    else:
        nivo_event_pfs = 1 - nivo["PFS_CNSR"].astype(float)
        nivo_event_os = 1 - nivo["OS_CNSR"].astype(float)
        eve_event_pfs = 1 - eve["PFS_CNSR"].astype(float)
        eve_event_os = 1 - eve["OS_CNSR"].astype(float)

    cox_nivo_pfs = cox_hr(nivo["PFS"], nivo_event_pfs, zscore(nivo["CLDN4"]), "nivo_PFS")
    cox_nivo_os = cox_hr(nivo["OS"], nivo_event_os, zscore(nivo["CLDN4"]), "nivo_OS")
    cox_eve_pfs = cox_hr(eve["PFS"], eve_event_pfs, zscore(eve["CLDN4"]), "everolimus_PFS")
    cox_eve_os = cox_hr(eve["OS"], eve_event_os, zscore(eve["CLDN4"]), "everolimus_OS")

    # reference genes on primary endpoint
    ref_rows = []
    for g in REF_GENES:
        if g not in prim.columns:
            continue
        r = logit_or(prim["y"].to_numpy(), zscore(prim[g]).to_numpy(), g)
        r["gene"] = g
        ref_rows.append(r)
    ref_df = pd.DataFrame(ref_rows)

    # correlations for interpretation
    corr = {}
    for g in ["Purity", "CD8A", "GZMB", "TACSTD2", "EPCAM", "CLDN3", "CLDN7"]:
        if g in nivo.columns:
            a = pd.to_numeric(nivo["CLDN4"], errors="coerce")
            b = pd.to_numeric(nivo[g], errors="coerce")
            ok = a.notna() & b.notna()
            if ok.sum() >= 10:
                rho, p = stats.spearmanr(a[ok], b[ok])
                corr[g] = {"rho": float(rho), "p": float(p), "n": int(ok.sum())}

    # --- transcriptome-wide logistic on primary samples ---
    sample_ids = prim["RNA_ID"].astype(str).tolist()
    y_prim = prim["y"].to_numpy(dtype=float)
    # restrict expression to primary samples
    mat = expr.loc[:, sample_ids].apply(pd.to_numeric, errors="coerce")
    gw_or = []
    gw_p = []
    gw_genes = []
    for gene, row in mat.iterrows():
        x = row.to_numpy(dtype=float)
        if np.isfinite(x).sum() < 20 or np.nanstd(x) == 0:
            continue
        z = (x - np.nanmean(x)) / np.nanstd(x)
        # fast IRLS
        X = np.column_stack([np.ones(z.size), z])
        beta = np.zeros(2)
        ok_fit = True
        hess = None
        for _ in range(20):
            eta = np.clip(X @ beta, -30, 30)
            p = 1.0 / (1.0 + np.exp(-eta))
            w = np.clip(p * (1.0 - p), 1e-12, None)
            WX = X * w[:, None]
            hess = X.T @ WX
            try:
                delta = np.linalg.solve(hess, X.T @ (y_prim - p))
            except np.linalg.LinAlgError:
                ok_fit = False
                break
            beta = beta + delta
            if np.max(np.abs(delta)) < 1e-7:
                break
        if not ok_fit or hess is None:
            continue
        try:
            se = math.sqrt(np.linalg.inv(hess)[1, 1])
        except np.linalg.LinAlgError:
            continue
        if not np.isfinite(se) or se == 0:
            continue
        pval = float(2 * stats.norm.sf(abs(beta[1] / se)))
        gw_genes.append(str(gene))
        gw_or.append(float(np.exp(beta[1])))
        gw_p.append(pval)

    gw_p = np.asarray(gw_p)
    gw_or = np.asarray(gw_or)
    gw_genes = np.asarray(gw_genes)
    if "CLDN4" not in set(gw_genes):
        cldn4_rank = None
        cldn4_fdr = None
        n_tested = int(gw_p.size)
    else:
        order = np.argsort(gw_p)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, order.size + 1)
        idx = int(np.where(gw_genes == "CLDN4")[0][0])
        cldn4_rank = int(ranks[idx])
        _, fdr, _, _ = multipletests(gw_p, method="fdr_bh")
        cldn4_fdr = float(fdr[idx])
        n_tested = int(gw_p.size)

    # top 15 by p (for honesty: show what actually ranks)
    top_idx = np.argsort(gw_p)[:15]
    gw_top = pd.DataFrame(
        {"gene": gw_genes[top_idx], "or_per_sd": gw_or[top_idx], "p": gw_p[top_idx]}
    )
    if n_tested:
        _, fdr_all, _, _ = multipletests(gw_p, method="fdr_bh")
        gw_top["fdr"] = fdr_all[top_idx]
        n_p05 = int((gw_p < 0.05).sum())
        n_fdr05 = int((fdr_all < 0.05).sum())
    else:
        n_p05 = 0
        n_fdr05 = 0

    # --- JAVELIN ---
    jav_clin = pd.read_csv(DATA / "javelin101" / "S11_clinical.tsv", sep="\t", skiprows=1)
    # expression: extract CLDN4 line only
    jav_expr_path = DATA / "javelin101" / "S13_tpm.tsv"
    jav_cldn4 = None
    jav_header = None
    with open(jav_expr_path) as fh:
        next(fh)  # title
        jav_header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.startswith("CLDN4\t"):
                parts = line.rstrip("\n").split("\t")
                jav_cldn4 = np.array([float(x) if x not in ("", "NA") else np.nan for x in parts[1:]])
                break
    if jav_cldn4 is None:
        raise SystemExit("CLDN4 not found in JAVELIN TPM matrix")
    jav_ids = jav_header[1:]
    jav = pd.DataFrame({"ID": jav_ids, "CLDN4": jav_cldn4})
    jav_clin["ID"] = jav_clin["ID"].astype(str)
    jav = jav.merge(jav_clin, on="ID", how="inner")
    # infer CNSR coding: CDISC typically 1=censored. Check event rate.
    # If PFS_P_CNSR==1 is common among short PFS, could be event. Compare mean PFS.
    mean_pfs_c1 = float(jav.loc[jav["PFS_P_CNSR"] == 1, "PFS_P"].mean())
    mean_pfs_c0 = float(jav.loc[jav["PFS_P_CNSR"] == 0, "PFS_P"].mean())
    # longer time with CNSR=1 => 1 is censored (CDISC)
    jav_cnsr1_censored = mean_pfs_c1 > mean_pfs_c0
    jav["event"] = (1 - jav["PFS_P_CNSR"].astype(float)) if jav_cnsr1_censored else jav["PFS_P_CNSR"].astype(float)
    jav_res = {}
    for arm, name in [
        ("Avelumab+Axitinib", "javelin_ave_axi_PFS"),
        ("Sunitinib", "javelin_sunitinib_PFS"),
    ]:
        sub = jav[jav["TRT01P"] == arm]
        jav_res[name] = cox_hr(sub["PFS_P"], sub["event"], zscore(sub["CLDN4"]), name)
    # interaction
    both = jav[jav["TRT01P"].isin(["Avelumab+Axitinib", "Sunitinib"])].copy()
    both["z"] = zscore(both["CLDN4"])
    both["ici"] = (both["TRT01P"] == "Avelumab+Axitinib").astype(float)
    both["zx"] = both["z"] * both["ici"]
    inter = {"n": int(len(both))}
    try:
        fit = PHReg(both["PFS_P"].astype(float), both[["z", "ici", "zx"]], status=both["event"].astype(float)).fit(disp=0)
        params = np.asarray(fit.params).ravel()
        bse = np.asarray(fit.bse).ravel()
        # column order z, ici, zx
        zstat = params[2] / bse[2]
        inter.update(
            {
                "hr_interaction": float(np.exp(params[2])),
                "p_interaction": float(2 * stats.norm.sf(abs(zstat))),
                "beta_interaction": float(params[2]),
            }
        )
    except Exception as exc:  # noqa: BLE001
        inter["error"] = str(exc)

    # --- GSE67501 ---
    import gzip

    gse_path = DATA / "gse67501" / "GSE67501_series_matrix.txt.gz"
    titles = recist = None
    table_lines = []
    with gzip.open(gse_path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            if "complete response (cr) or partial response" in line.lower():
                recist = [x.strip().strip('"').split(": ")[-1] for x in line.rstrip().split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break
        header = next(fh)
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            table_lines.append(line)
    probes = {}
    for line in table_lines:
        pid, *vals = line.rstrip().split("\t")
        pid = pid.strip('"')
        if pid in {"ILMN_2132458"}:  # CLDN4
            probes[pid] = np.array([float(v.strip('"')) for v in vals])
    gse = pd.DataFrame(
        {
            "title": titles,
            "recist": recist,
            "CLDN4": probes["ILMN_2132458"],
        }
    )
    # paper titles: anti-PD-1_Response_* vs anti-PD-1_No_Response_*
    gse["binary_paper"] = ~gse["title"].str.contains("No_Response")
    gse["y_crpr_vs_nr"] = gse["recist"].isin(["CR", "PR"]).astype(int)
    gse_keep = gse[gse["recist"].isin(["CR", "PR", "NR"])]
    gse_primary = logit_or(
        gse_keep["y_crpr_vs_nr"].to_numpy(),
        zscore(gse_keep["CLDN4"]).to_numpy(),
        "GSE67501_CRPR_vs_NR",
    )
    gse_paper = logit_or(
        gse["binary_paper"].astype(int).to_numpy(),
        zscore(gse["CLDN4"]).to_numpy(),
        "GSE67501_paper_response_vs_noresponse",
    )
    gse_mwu = mwu_auc(
        gse.loc[gse["binary_paper"], "CLDN4"].to_numpy(),
        gse.loc[~gse["binary_paper"], "CLDN4"].to_numpy(),
    )

    # power
    n = primary["n"]
    p_event = primary["n_resp"] / n if n else np.nan
    power = {
        "n": n,
        "n_resp": primary["n_resp"],
        "n_pd": primary["n_non"],
        "event_rate": p_event,
        "power_OR_0.42": power_wald(n, p_event, 0.42),
        "power_OR_0.50": power_wald(n, p_event, 0.50),
        "power_OR_0.67": power_wald(n, p_event, 0.67),
        "min_OR_80pct_power": None,
    }
    # OR closest to 1 (below 1) that still has ≥80% power
    detectable = None
    for cand in np.linspace(0.99, 0.20, 80):
        if power_wald(n, p_event, float(cand)) >= 0.80:
            detectable = float(cand)
            break
    power["or_below_1_with_80pct_power"] = detectable
    power["or_above_1_with_80pct_power"] = (float(round(1.0 / detectable, 3)) if detectable else None)
    power.pop("min_OR_80pct_power", None)
    power.pop("max_OR_80pct_power_protective", None)
    power.pop("min_OR_80pct_power_harmful", None)

    # save patient-level (small, reproducible)
    prim_out = prim[
        ["SUBJID", "RNA_ID", "Cohort", "Arm", "ORR", "Benefit", "PFS", "PFS_CNSR", "OS", "OS_CNSR", "Purity", "MSKCC", "CLDN4", "z_CLDN4", "y"]
    ].copy()
    prim_out.to_csv(OUT / "braun_primary_patients.tsv", sep="\t", index=False)
    nivo[
        ["SUBJID", "RNA_ID", "Cohort", "Arm", "ORR", "Benefit", "PFS", "PFS_CNSR", "OS", "OS_CNSR", "Purity", "MSKCC", "CLDN4"]
        + [g for g in REF_GENES if g in nivo.columns and g != "CLDN4"]
    ].to_csv(OUT / "braun_nivo_patients.tsv", sep="\t", index=False)
    gse.to_csv(OUT / "gse67501_patients.tsv", sep="\t", index=False)
    jav[["ID", "TRT01P", "PFS_P", "PFS_P_CNSR", "event", "CLDN4", "AGE", "SEX", "PDL1FL"]].to_csv(
        OUT / "javelin_patients.tsv", sep="\t", index=False
    )
    ref_df.to_csv(OUT / "reference_genes_primary.tsv", sep="\t", index=False)
    gw_top.to_csv(OUT / "genomewide_top15.tsv", sep="\t", index=False)

    # figures
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    data = [prim.loc[prim["y"] == 0, "CLDN4"], prim.loc[prim["y"] == 1, "CLDN4"]]
    bp = ax.boxplot(
        data,
        tick_labels=["PD (n={})".format(primary["n_non"]), "CR/PR (n={})".format(primary["n_resp"])],
        patch_artist=True,
    )
    for patch, c in zip(bp["boxes"], ["#c0392b", "#1e8449"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    rng = np.random.default_rng(1)
    for i, d in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(d)), d, s=14, c="k", alpha=0.55, zorder=3)
    ax.set_ylabel("CLDN4 (Braun normalized log2)")
    ax.set_title("Braun CheckMate nivo — primary endpoint")
    ax.axhline(med, ls="--", c="gray", lw=0.8, label="median")
    fig.tight_layout()
    fig.savefig(OUT / "fig_braun_primary_boxplot.png", dpi=160)
    plt.close(fig)

    # forest of ORs
    forest_items = [
        ("PRIMARY nivo CR/PR vs PD", primary),
        ("nivo CR/PR vs SD+PD", alt_orr),
        ("nivo CB vs NCB", alt_cb),
        ("everolimus CR/PR vs PD", eve_logit),
        ("GSE67501 CR/PR vs NR", gse_primary),
        ("GSE67501 paper binary", gse_paper),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ys = []
    labs = []
    k = 0
    for lab, r in forest_items:
        if r.get("or_per_sd") is None or not np.isfinite(r.get("or_per_sd", np.nan)):
            continue
        k += 1
        ys.append(k)
        labs.append(f"{lab}  n={r['n']} ({r['n_resp']}/{r['n_non']})")
        ax.errorbar(
            r["or_per_sd"],
            k,
            xerr=[[r["or_per_sd"] - r["or_lo"]], [r["or_hi"] - r["or_per_sd"]]],
            fmt="o",
            color="#1a5276",
            capsize=3,
        )
    ax.axvline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels(labs)
    ax.set_xlabel("Odds ratio per +1 SD CLDN4 (response)")
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_title("CLDN4 vs response — pre-specified logistic models")
    fig.tight_layout()
    fig.savefig(OUT / "fig_or_forest.png", dpi=160)
    plt.close(fig)

    summary = {
        "claim": "B5 analog, KIRC/RCC: CLDN4-high associated with worse ICI response (user meta OR=0.42)",
        "verdict": None,  # filled below
        "primary": primary,
        "primary_mwu": primary_mwu,
        "median_split": {
            **med_or,
            "table_high_resp_non": tab_med[0].tolist(),
            "table_low_resp_non": tab_med[1].tolist(),
            "median": med,
        },
        "tertile_split_hi_vs_lo": {
            **ter_or,
            "table_hi_resp_non": tab_ter[0].tolist(),
            "table_lo_resp_non": tab_ter[1].tolist(),
        },
        "alt_CRPR_vs_SDPD": alt_orr,
        "alt_CB_vs_NCB": alt_cb,
        "multivariable_cohort_purity": mv_res,
        "multivariable_MSKCC": mv2_res,
        "everolimus_control": eve_logit,
        "cox": {
            "nivo_PFS": cox_nivo_pfs,
            "nivo_OS": cox_nivo_os,
            "everolimus_PFS": cox_eve_pfs,
            "everolimus_OS": cox_eve_os,
            "braun_cnsr_coding": coding,
        },
        "javelin": {
            **jav_res,
            "interaction": inter,
            "cnsr1_is_censored": jav_cnsr1_censored,
            "mean_PFS_CNSR1": mean_pfs_c1,
            "mean_PFS_CNSR0": mean_pfs_c0,
            "n_ave_axi": int((jav["TRT01P"] == "Avelumab+Axitinib").sum()),
            "n_sunitinib": int((jav["TRT01P"] == "Sunitinib").sum()),
            "note": "JAVELIN public supplement has PFS only; no RECIST response, so no OR.",
        },
        "gse67501": {
            "primary_CRPR_vs_NR": gse_primary,
            "paper_binary": gse_paper,
            "mwu_paper_binary": gse_mwu,
            "n_total": int(len(gse)),
            "recist_counts": gse["recist"].value_counts().to_dict(),
            "probe": "ILMN_2132458",
        },
        "genomewide": {
            "n_genes_tested": n_tested,
            "cldn4_p_rank": cldn4_rank,
            "cldn4_fdr": cldn4_fdr,
            "n_genes_p_lt_0.05": n_p05,
            "n_genes_fdr_lt_0.05": n_fdr05,
        },
        "reference_genes": ref_df[["gene", "n", "n_resp", "n_non", "or_per_sd", "or_lo", "or_hi", "p"]].to_dict(
            orient="records"
        ),
        "cldn4_correlations_nivo": corr,
        "power": power,
        "inventory": {
            "braun_rna_nivo": int(len(nivo)),
            "braun_rna_everolimus": int(len(eve)),
            "braun_primary_n": int(len(prim)),
            "javelin_rna": int(len(jav)),
            "gse67501": int(len(gse)),
        },
        "sources": {
            "braun": "https://www.nature.com/articles/s41591-020-0839-y  (Supplementary Table 1 + 4A)",
            "javelin": "https://www.nature.com/articles/s41591-020-1044-8  (Supplementary Table 11 + 13)",
            "gse67501": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501  (Ascierto 2016 CIR)",
        },
    }

    # honest verdict
    p = primary["p"]
    or_ = primary["or_per_sd"]
    if not np.isfinite(p):
        verdict = "INCONCLUSIVE — primary model did not converge"
    elif p >= 0.05:
        verdict = (
            "NOT SUPPORTED in open KIRC ICI RNA. Primary Braun nivolumab CR/PR vs PD: "
            f"OR={or_:.3f} (95% CI {primary['or_lo']:.3f}–{primary['or_hi']:.3f}), "
            f"p={p:.3g}, n={primary['n']} ({primary['n_resp']} resp / {primary['n_non']} PD). "
            "Null. Does not reproduce user B5 OR=0.42 in RCC."
        )
    elif or_ < 1:
        verdict = (
            "DIRECTIONALLY CONSISTENT with B5 (CLDN4-high worse) in Braun nivo, "
            f"OR={or_:.3f}, p={p:.3g}, n={primary['n']}. "
            "See genome-wide rank and JAVELIN/GSE67501 before treating as confirmation."
        )
    else:
        verdict = (
            "OPPOSITE DIRECTION to B5 (CLDN4-high better) in Braun nivo, "
            f"OR={or_:.3f}, p={p:.3g}, n={primary['n']}. "
            "Does not support the user claim in RCC."
        )
    summary["verdict"] = verdict

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            if not np.isfinite(o):
                return None
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    with open(OUT / "summary.json", "w") as fh:
        json.dump(_clean(summary), fh, indent=2)
        fh.write("\n")

    print(json.dumps(_clean({"verdict": verdict, "primary": primary, "genomewide": summary["genomewide"], "power": power}), indent=2))


if __name__ == "__main__":
    main()
