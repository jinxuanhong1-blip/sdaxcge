#!/usr/bin/env python3
"""TCGA: share of the keratin-adjusted TACSTD2–immune correlation accounted for by CLDN4.

Pre-specified before looking at the mediation output
----------------------------------------------------
Expression: UCSC Xena GDC STAR log2(TPM+1), GENCODE v36. Primary solid tumor
only (sample type 01). Replicate aliquots averaged to the patient.

Cohorts, in table order:
  lung + keratin funnel (8): LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD
  keratin funnel (7): LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD

Keratin covariates: KRT8, KRT18, KRT19.
CD8 score: mean of CD8A and CD8B.
CD3 score: mean of CD3D, CD3E, and CD3G.
Association: partial Spearman. Pooled with DerSimonian–Laird on Fisher z.

The keratin-only model is fit first. CLDN4 is then added. CLDN7 and EPCAM are
added one at a time as negative controls. A further model adds CLDN4 on top
of CLDN7 and EPCAM.

Mediation percentage uses the pooled correlations:
    100 * (rho_keratin - rho_keratin_plus_covariate) / rho_keratin
The interval is the 2.5 and 97.5 percentiles of that quantity across patient
bootstraps done inside each cohort. Rank-OLS proportion mediated (a*b/c on
average ranks) is reported beside it. The path is TACSTD2 → CLDN4 → CD8;
CD3 gets the same arithmetic.

LUSC is also fit with KRT5+KRT6A+KRT14 in place of KRT8/18/19. That row is
not pooled.

This is a bulk-RNA decomposition. It does not measure spatial exclusion,
a knockdown, or an ICI-treated cohort.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bootstrap import (
    KERATIN,
    MEDIATORS,
    OUTCOMES,
    PREDICTOR,
    RHO_MODELS,
    SQUAMOUS,
    SQUAMOUS_MODELS,
    bootstrap_design,
    bootstrap_primary,
)
from io_tcga import extract_cohort, load_cached, load_probemap
from stats import (
    attenuation_fraction,
    bh_fdr,
    coef_given,
    covariate_matrix,
    fisher_ci,
    partial_spearman,
    product_mediation,
    random_effects_meta,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "tcga_trop2_cldn4_mediation")
TABLES = os.path.join(OUT, "tables")
FIGS = os.path.join(OUT, "figures")
CACHE = os.environ.get("TCGA_MED_CACHE", "/tmp/tcga_trop2_cldn4_mediation")

COHORTS = ["LUAD", "LUSC", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
FUNNEL7 = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
POOLS = {"lung_funnel8": COHORTS, "keratin_funnel7": FUNNEL7}
POOL_LABEL = {
    "lung_funnel8": "Pooled, 8 cohorts",
    "keratin_funnel7": "Pooled, 7-cohort funnel",
}

GENE_SYMBOLS = [
    "TACSTD2", "CLDN4", "CLDN7", "EPCAM",
    "KRT8", "KRT18", "KRT19", "KRT5", "KRT6A", "KRT14",
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G",
]
PRIMARY_GENES = [
    "TACSTD2", "CLDN4", "CLDN7", "EPCAM",
    "KRT8", "KRT18", "KRT19",
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G",
]
ADDED = {
    "keratin_CLDN4": "CLDN4",
    "keratin_CLDN7": "CLDN7",
    "keratin_EPCAM": "EPCAM",
    "keratin_CLDN7_EPCAM": "CLDN7+EPCAM",
    "keratin_CLDN7_EPCAM_CLDN4": "CLDN7+EPCAM+CLDN4",
}
SHARE_MODELS = ["keratin_CLDN4", "keratin_CLDN7", "keratin_EPCAM"]
CONTROL_MODELS = ["keratin_CLDN7", "keratin_EPCAM"]

# Three-decimal TACSTD2–CD8 partial rho | KRT8+KRT18+KRT19 printed by
# scripts/tcga_trop2_cldn_keratin/run_analysis.py. Reproduction check only.
PRIOR_CD8_KERATIN_RHO = {
    "LUAD": -0.104,
    "LUSC": -0.221,
    "BRCA": -0.039,
    "CESC": 0.118,
    "KIRC": -0.047,
    "STAD": -0.079,
    "BLCA": -0.158,
    "PAAD": 0.029,
}
PRIOR_N = {
    "LUAD": 516,
    "LUSC": 501,
    "BRCA": 1095,
    "CESC": 304,
    "KIRC": 533,
    "STAD": 412,
    "BLCA": 406,
    "PAAD": 178,
}

COLOR = {
    "keratin": "#222222",
    "keratin_CLDN4": "#0072B2",
    "keratin_CLDN7": "#D55E00",
    "keratin_EPCAM": "#009E73",
    "CLDN4": "#0072B2",
    "CLDN7": "#D55E00",
    "EPCAM": "#009E73",
}
MODEL_LABEL = {
    "keratin": "KRT8/18/19",
    "keratin_CLDN4": "+ CLDN4",
    "keratin_CLDN7": "+ CLDN7",
    "keratin_EPCAM": "+ EPCAM",
}
Q_FAMILIES = []
for _outcome in OUTCOMES:
    Q_FAMILIES.append((_outcome, ["keratin"]))
    Q_FAMILIES.append((_outcome, ["keratin_CLDN4"]))
    Q_FAMILIES.append((_outcome, ["keratin_CLDN7", "keratin_EPCAM"]))
    Q_FAMILIES.append((_outcome, ["keratin_CLDN7_EPCAM"]))
    Q_FAMILIES.append((_outcome, ["keratin_CLDN7_EPCAM_CLDN4"]))


def fmt_num(x, digits=3):
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(x):
    if not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def fmt_pct(x):
    """Format a fraction as a percentage with one decimal."""
    if not np.isfinite(x):
        return "NA"
    return f"{100.0 * x:.1f}%"


def fmt_i2(x):
    if not np.isfinite(x):
        return "NA"
    return f"{100.0 * x:.1f}%"


def fmt_ci(lo, hi, percent=False):
    if percent:
        return f"{fmt_pct(lo)} to {fmt_pct(hi)}"
    return f"{fmt_num(lo)} to {fmt_num(hi)}"


def percentile_ci(draws, alpha=0.05):
    d = np.asarray(draws, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < 50:
        return np.nan, np.nan, int(d.size)
    lo, hi = np.percentile(d, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi), int(d.size)


def md_table(headers, rows):
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def frame_to_arrays(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    data = {}
    for gene in ["TACSTD2", "CLDN4", "CLDN7", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT6A", "KRT14"]:
        data[gene] = frame[gene].to_numpy(dtype=float)
    data["CD8_score"] = frame[["CD8A", "CD8B"]].mean(axis=1).to_numpy(dtype=float)
    data["CD3_score"] = frame[["CD3D", "CD3E", "CD3G"]].mean(axis=1).to_numpy(dtype=float)
    return data


def load_cohorts(mapping):
    id_to_gene = {ens: sym for sym, ens in mapping.items()}
    provenance = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(extract_cohort, cohort, id_to_gene, CACHE): cohort for cohort in COHORTS
        }
        by_cohort = {}
        for fut in as_completed(futures):
            cohort = futures[fut]
            by_cohort[cohort] = fut.result()
            print(f"[ready] {cohort} n={by_cohort[cohort]['n_patients']} cache={by_cohort[cohort]['from_cache']}", flush=True)
    matrices = {}
    counts = []
    for cohort in COHORTS:
        frame = load_cached(cohort, CACHE)
        n_extract = int(frame.shape[0])
        before = frame[PRIMARY_GENES + ["KRT5", "KRT6A", "KRT14"]]
        n_na = int(before.isna().any(axis=1).sum())
        kept = frame.dropna(subset=PRIMARY_GENES + ["KRT5", "KRT6A", "KRT14"]).copy()
        matrices[cohort] = frame_to_arrays(kept)
        counts.append({
            "cohort": cohort,
            "n_extracted": n_extract,
            "n_complete": int(kept.shape[0]),
            "n_dropped_nonfinite": n_na,
            "n_01_columns": by_cohort[cohort]["n_01_columns"],
            "prior_n": PRIOR_N[cohort],
            "sha256": by_cohort[cohort]["sha256"],
            "url": by_cohort[cohort]["url"],
            "from_cache": by_cohort[cohort]["from_cache"],
        })
        provenance.append(by_cohort[cohort])
    return matrices, pd.DataFrame(counts), provenance


def point_partials(matrices):
    rows = []
    rho = {cohort: {outcome: {} for outcome in OUTCOMES} for cohort in COHORTS}
    for cohort in COHORTS:
        data = matrices[cohort]
        n = int(data["TACSTD2"].size)
        for outcome in OUTCOMES:
            y = data[outcome]
            for model, cov_names in RHO_MODELS.items():
                covariates = [data[name] for name in cov_names]
                r, p, nobs = partial_spearman(data["TACSTD2"], y, covariates)
                lo, hi = fisher_ci(r, nobs, len(cov_names))
                rho[cohort][outcome][model] = r
                rows.append({
                    "cohort": cohort,
                    "outcome": outcome,
                    "model": model,
                    "k": len(cov_names),
                    "n": nobs,
                    "rho": r,
                    "p": p,
                    "ci_low": lo,
                    "ci_high": hi,
                    "q": np.nan,
                })
                if nobs != n:
                    raise RuntimeError(f"{cohort} {outcome} {model} dropped rows after the complete-case filter")
    table = pd.DataFrame(rows)
    for outcome, models in Q_FAMILIES:
        mask = (table["outcome"] == outcome) & (table["model"].isin(models))
        table.loc[mask, "q"] = bh_fdr(table.loc[mask, "p"].to_numpy())
    return table, rho


def point_paths_and_pm(matrices):
    path_rows = []
    pm_rows = []
    pm = {cohort: {outcome: {} for outcome in OUTCOMES} for cohort in COHORTS}
    for cohort in COHORTS:
        data = matrices[cohort]
        keratin = [data[name] for name in KERATIN]
        for mediator in MEDIATORS:
            a_rho, a_p, a_n = partial_spearman(data["TACSTD2"], data[mediator], keratin)
            a_lo, a_hi = fisher_ci(a_rho, a_n, len(KERATIN))
            path_rows.append({
                "cohort": cohort,
                "outcome": "",
                "mediator": mediator,
                "path": "a",
                "rho": a_rho,
                "p": a_p,
                "n": a_n,
                "k": len(KERATIN),
                "ci_low": a_lo,
                "ci_high": a_hi,
            })
            for outcome in OUTCOMES:
                b_cov = keratin + [data["TACSTD2"]]
                b_rho, b_p, b_n = partial_spearman(data[mediator], data[outcome], b_cov)
                b_lo, b_hi = fisher_ci(b_rho, b_n, len(KERATIN) + 1)
                path_rows.append({
                    "cohort": cohort,
                    "outcome": outcome,
                    "mediator": mediator,
                    "path": "b",
                    "rho": b_rho,
                    "p": b_p,
                    "n": b_n,
                    "k": len(KERATIN) + 1,
                    "ci_low": b_lo,
                    "ci_high": b_hi,
                })
                fit = product_mediation(data["TACSTD2"], data[mediator], data[outcome], keratin, rank=True)
                if abs(fit["identity_gap"]) > 1e-6:
                    raise RuntimeError(f"mediation identity failed for {cohort} {outcome} {mediator}")
                pm[cohort][outcome][mediator] = fit["proportion"]
                pm_rows.append({
                    "cohort": cohort,
                    "outcome": outcome,
                    "mediator": mediator,
                    "a": fit["a"],
                    "b": fit["b"],
                    "c": fit["c"],
                    "c_prime": fit["c_prime"],
                    "indirect": fit["indirect"],
                    "proportion": fit["proportion"],
                    "identity_gap": fit["identity_gap"],
                    "n": fit["n"],
                    "k": fit["k"],
                })
    path = pd.DataFrame(path_rows)
    for path_name, k_expected in (("a", 3), ("b", 4)):
        for mediator in MEDIATORS:
            for outcome in ([""] if path_name == "a" else OUTCOMES):
                mask = (path["path"] == path_name) & (path["mediator"] == mediator)
                if path_name == "b":
                    mask &= path["outcome"] == outcome
                block = path.loc[mask]
                if block.empty:
                    continue
                q = bh_fdr(block["p"].to_numpy())
                path.loc[block.index, "q"] = q
    return path, pd.DataFrame(pm_rows), pm


def independent_incremental(data):
    """Rank-OLS incremental fraction, separate from the bootstrap implementation."""
    names = KERATIN + ["TACSTD2", "CLDN4", "CLDN7", "EPCAM", "CD8_score", "CD3_score"]
    ranked = {name: stats.rankdata(data[name]).astype(float) for name in names}
    out = {}
    Zk = covariate_matrix([ranked[name] for name in KERATIN], len(ranked["TACSTD2"]))
    Zc = covariate_matrix([ranked[name] for name in KERATIN + ["CLDN7", "EPCAM"]], len(ranked["TACSTD2"]))
    Zf = covariate_matrix([ranked[name] for name in KERATIN + ["CLDN7", "EPCAM", "CLDN4"]], len(ranked["TACSTD2"]))
    x = ranked["TACSTD2"]
    for outcome in OUTCOMES:
        y = ranked[outcome]
        c = coef_given(y, x, Zk)
        c_ctrl = coef_given(y, x, Zc)
        c_full = coef_given(y, x, Zf)
        out[outcome] = (c_ctrl - c_full) / c if abs(c) > 1e-12 else np.nan
    return out


def run_bootstraps(matrices, n_boot, seed):
    boots = {}
    for i, cohort in enumerate(COHORTS):
        cohort_seed = seed + 1000 * (i + 1)
        print(f"[boot] {cohort} B={n_boot} seed={cohort_seed}", flush=True)
        boots[cohort] = bootstrap_primary(matrices[cohort], n_boot, cohort_seed)
        inc = independent_incremental(matrices[cohort])
        for outcome in OUTCOMES:
            if abs(inc[outcome] - boots[cohort]["full_incremental"][outcome]) > 1e-6:
                raise RuntimeError(f"incremental mismatch in {cohort} {outcome}")
    print(f"[boot] LUSC squamous B={n_boot}", flush=True)
    squamous = bootstrap_design(
        data=matrices["LUSC"],
        models=SQUAMOUS_MODELS,
        outcomes=OUTCOMES,
        mediators=MEDIATORS,
        predictor=PREDICTOR,
        base_covariates=SQUAMOUS,
        n_boot=n_boot,
        seed=seed + 9000,
        incremental_pair=(SQUAMOUS + ["CLDN7", "EPCAM"], SQUAMOUS + ["CLDN7", "EPCAM", "CLDN4"]),
    )
    return boots, squamous


def agree_with_points(boots, rho, pm):
    worst = 0.0
    for cohort in COHORTS:
        for outcome in OUTCOMES:
            for model in RHO_MODELS:
                gap = abs(boots[cohort]["full_rho"][outcome][model] - rho[cohort][outcome][model])
                worst = max(worst, gap)
            for mediator in MEDIATORS:
                gap = abs(boots[cohort]["full_pm"][outcome][mediator] - pm[cohort][outcome][mediator])
                worst = max(worst, gap)
    if worst > 1e-6:
        raise RuntimeError(f"bootstrap full-sample metrics disagree with partial Spearman / rank-OLS by {worst}")
    return worst


def meta_draws(boots, cohorts, ns, outcome, model, k):
    n_boot = boots[cohorts[0]]["n_boot"]
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        rhos = [boots[c]["rho"][outcome][model][b] for c in cohorts]
        if not np.all(np.isfinite(rhos)):
            draws[b] = np.nan
            continue
        draws[b] = random_effects_meta(rhos, ns, k)["rho"]
    return draws


def pm_draws(boots, cohorts, ns, outcome, mediator):
    n_boot = boots[cohorts[0]]["n_boot"]
    weight = np.asarray(ns, dtype=float)
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        vals = np.array([boots[c]["pm"][outcome][mediator][b] for c in cohorts], dtype=float)
        if not np.isfinite(vals).all():
            draws[b] = np.nan
            continue
        draws[b] = float(np.sum(weight * vals) / np.sum(weight))
    return draws


def incremental_draws(boots, cohorts, ns, outcome):
    n_boot = boots[cohorts[0]]["n_boot"]
    weight = np.asarray(ns, dtype=float)
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        vals = np.array([boots[c]["incremental"][outcome][b] for c in cohorts], dtype=float)
        if not np.isfinite(vals).all():
            draws[b] = np.nan
            continue
        draws[b] = float(np.sum(weight * vals) / np.sum(weight))
    return draws


def assemble(partial_table, rho, pm, boots, counts):
    ns = {row.cohort: int(row.n_complete) for row in counts.itertuples()}
    pooled_rows = []
    share_rows = []
    pm_pool_rows = []
    diff_rows = []
    incremental_rows = []
    draw_cache = {}
    for pool_name, cohorts in POOLS.items():
        nvec = [ns[c] for c in cohorts]
        for outcome in OUTCOMES:
            for model, cov_names in RHO_MODELS.items():
                rhos = [rho[c][outcome][model] for c in cohorts]
                meta = random_effects_meta(rhos, nvec, len(cov_names))
                if meta["n_cohorts"] != len(cohorts):
                    raise RuntimeError(f"meta dropped cohorts for {pool_name} {outcome} {model}")
                draws = meta_draws(boots, cohorts, nvec, outcome, model, len(cov_names))
                draw_cache[(pool_name, outcome, model)] = draws
                pooled_rows.append({
                    "pool": pool_name,
                    "outcome": outcome,
                    "model": model,
                    "k": len(cov_names),
                    "n_cohorts": meta["n_cohorts"],
                    "n_patients": int(sum(nvec)),
                    "n_negative": int(sum(r < 0 for r in rhos)),
                    "rho": meta["rho"],
                    "ci_low": meta["ci_low"],
                    "ci_high": meta["ci_high"],
                    "p": meta["p"],
                    "I2": meta["I2"],
                    "tau2": meta["tau2"],
                    "Q": meta["Q"],
                })
            before = draw_cache[(pool_name, outcome, "keratin")]
            rho_before = next(
                row["rho"] for row in pooled_rows
                if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == "keratin"
            )
            for model in list(ADDED):
                after = draw_cache[(pool_name, outcome, model)]
                rho_after = next(
                    row["rho"] for row in pooled_rows
                    if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == model
                )
                delta = after - before
                with np.errstate(divide="ignore", invalid="ignore"):
                    share = (before - after) / before
                share = np.where(np.isfinite(before) & (np.abs(before) >= 1e-8), share, np.nan)
                d_lo, d_hi, d_n = percentile_ci(delta)
                s_lo, s_hi, s_n = percentile_ci(share)
                near = float(np.mean(np.isfinite(before) & (np.abs(before) < 0.02)))
                share_rows.append({
                    "pool": pool_name,
                    "outcome": outcome,
                    "model": model,
                    "added": ADDED[model],
                    "rho_before": rho_before,
                    "rho_after": rho_after,
                    "delta_rho": rho_after - rho_before,
                    "delta_lo": d_lo,
                    "delta_hi": d_hi,
                    "share": attenuation_fraction(rho_before, rho_after),
                    "share_lo": s_lo,
                    "share_hi": s_hi,
                    "near_zero_frac": near,
                    "n_boot_share": s_n,
                    "n_boot_delta": d_n,
                })
            for left, right in (("keratin_CLDN4", "keratin_CLDN7"), ("keratin_CLDN4", "keratin_EPCAM")):
                s_left = next(row for row in share_rows if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == left)
                s_right = next(row for row in share_rows if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == right)
                # Rebuild drawwise shares so the contrast uses paired resamples.
                a_draws = draw_cache[(pool_name, outcome, left)]
                b_draws = draw_cache[(pool_name, outcome, right)]
                with np.errstate(divide="ignore", invalid="ignore"):
                    share_l = (before - a_draws) / before
                    share_r = (before - b_draws) / before
                ok = np.isfinite(before) & (np.abs(before) >= 1e-8)
                contrast = np.where(ok, share_l - share_r, np.nan)
                c_lo, c_hi, c_n = percentile_ci(contrast)
                diff_rows.append({
                    "pool": pool_name,
                    "outcome": outcome,
                    "contrast": f"{ADDED[left]} minus {ADDED[right]}",
                    "share_difference": s_left["share"] - s_right["share"],
                    "ci_low": c_lo,
                    "ci_high": c_hi,
                    "n_boot": c_n,
                })
            for mediator in MEDIATORS:
                point_vals = [pm[c][outcome][mediator] for c in cohorts]
                point = float(np.sum(np.asarray(nvec) * np.asarray(point_vals)) / np.sum(nvec))
                draws = pm_draws(boots, cohorts, nvec, outcome, mediator)
                lo, hi, n_ok = percentile_ci(draws)
                pm_pool_rows.append({
                    "pool": pool_name,
                    "outcome": outcome,
                    "mediator": mediator,
                    "proportion": point,
                    "ci_low": lo,
                    "ci_high": hi,
                    "n_boot": n_ok,
                    "weight": "n",
                })
            inc_point_vals = []
            for cohort in cohorts:
                inc_point_vals.append(boots[cohort]["full_incremental"][outcome])
            inc_point = float(np.sum(np.asarray(nvec) * np.asarray(inc_point_vals, dtype=float)) / np.sum(nvec))
            inc_draw = incremental_draws(boots, cohorts, nvec, outcome)
            # Rho-scale incremental share of the original keratin correlation.
            rho_ctrl = next(row["rho"] for row in pooled_rows if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == "keratin_CLDN7_EPCAM")
            rho_full = next(row["rho"] for row in pooled_rows if row["pool"] == pool_name and row["outcome"] == outcome and row["model"] == "keratin_CLDN7_EPCAM_CLDN4")
            ctrl_draws = draw_cache[(pool_name, outcome, "keratin_CLDN7_EPCAM")]
            full_draws = draw_cache[(pool_name, outcome, "keratin_CLDN7_EPCAM_CLDN4")]
            with np.errstate(divide="ignore", invalid="ignore"):
                inc_share = (ctrl_draws - full_draws) / before
            inc_share = np.where(np.isfinite(before) & (np.abs(before) >= 1e-8), inc_share, np.nan)
            i_lo, i_hi, i_n = percentile_ci(inc_share)
            r_lo, r_hi, r_n = percentile_ci(inc_draw)
            incremental_rows.append({
                "pool": pool_name,
                "outcome": outcome,
                "rho_keratin": rho_before,
                "rho_controls": rho_ctrl,
                "rho_controls_plus_CLDN4": rho_full,
                # Further share of the original keratin correlation removed by CLDN4.
                "share_of_keratin_rho": float((rho_ctrl - rho_full) / rho_before) if abs(rho_before) >= 1e-8 else np.nan,
                "share_lo": i_lo,
                "share_hi": i_hi,
                "n_boot_share": i_n,
                "rank_ols_proportion": inc_point,
                "rank_ols_lo": r_lo,
                "rank_ols_hi": r_hi,
                "n_boot_rank": r_n,
            })
    return {
        "pooled": pd.DataFrame(pooled_rows),
        "share": pd.DataFrame(share_rows),
        "pm_pool": pd.DataFrame(pm_pool_rows),
        "diff": pd.DataFrame(diff_rows),
        "incremental": pd.DataFrame(incremental_rows),
        "ns": ns,
        "draw_cache": draw_cache,
    }


def cohort_share_table(rho, boots, ns):
    rows = []
    for cohort in COHORTS:
        for outcome in OUTCOMES:
            before_point = rho[cohort][outcome]["keratin"]
            before_draws = boots[cohort]["rho"][outcome]["keratin"]
            for model, added in ADDED.items():
                after_point = rho[cohort][outcome][model]
                after_draws = boots[cohort]["rho"][outcome][model]
                delta = after_draws - before_draws
                with np.errstate(divide="ignore", invalid="ignore"):
                    share = (before_draws - after_draws) / before_draws
                share = np.where(np.isfinite(before_draws) & (np.abs(before_draws) >= 1e-8), share, np.nan)
                d_lo, d_hi, _ = percentile_ci(delta)
                s_lo, s_hi, _ = percentile_ci(share)
                rows.append({
                    "cohort": cohort,
                    "outcome": outcome,
                    "model": model,
                    "added": added,
                    "n": ns[cohort],
                    "rho_before": before_point,
                    "rho_after": after_point,
                    "delta_rho": after_point - before_point,
                    "delta_lo": d_lo,
                    "delta_hi": d_hi,
                    "share": attenuation_fraction(before_point, after_point),
                    "share_lo": s_lo,
                    "share_hi": s_hi,
                    "small_denominator": abs(before_point) < 0.05,
                })
    return pd.DataFrame(rows)


def calibration_table(rho, counts):
    rows = []
    for cohort in COHORTS:
        observed = rho[cohort]["CD8_score"]["keratin"]
        prior = PRIOR_CD8_KERATIN_RHO[cohort]
        n = int(counts.loc[counts["cohort"] == cohort, "n_complete"].iloc[0])
        rows.append({
            "cohort": cohort,
            "n_complete": n,
            "prior_n": PRIOR_N[cohort],
            "rho": observed,
            "prior_rho_3dp": prior,
            "abs_diff": abs(observed - prior),
        })
    return pd.DataFrame(rows)


def squamous_table(matrices, squamous_boot):
    data = matrices["LUSC"]
    rows = []
    n = int(data["TACSTD2"].size)
    for outcome in OUTCOMES:
        rhos = {}
        for model, cov_names in SQUAMOUS_MODELS.items():
            r, p, nobs = partial_spearman(data["TACSTD2"], data[outcome], [data[name] for name in cov_names])
            if abs(r - squamous_boot["full_rho"][outcome][model]) > 1e-6:
                raise RuntimeError("LUSC squamous full-sample rho mismatch")
            lo, hi = fisher_ci(r, nobs, len(cov_names))
            rhos[model] = r
            rows.append({
                "outcome": outcome,
                "model": model,
                "k": len(cov_names),
                "n": nobs,
                "rho": r,
                "p": p,
                "ci_low": lo,
                "ci_high": hi,
                "share_of_squamous_rho": np.nan if model == "squamous" else attenuation_fraction(rhos["squamous"], r),
            })
        before = squamous_boot["rho"][outcome]["squamous"]
        for model in ["squamous_CLDN4", "squamous_CLDN7", "squamous_EPCAM"]:
            after = squamous_boot["rho"][outcome][model]
            with np.errstate(divide="ignore", invalid="ignore"):
                share = (before - after) / before
            share = np.where(np.isfinite(before) & (np.abs(before) >= 1e-8), share, np.nan)
            lo, hi, n_ok = percentile_ci(share)
            for row in rows:
                if row["outcome"] == outcome and row["model"] == model:
                    row["share_lo"] = lo
                    row["share_hi"] = hi
                    row["n_boot"] = n_ok
        fit = product_mediation(data["TACSTD2"], data["CLDN4"], data[outcome], [data[name] for name in SQUAMOUS], rank=True)
        for row in rows:
            if row["outcome"] == outcome and row["model"] == "squamous_CLDN4":
                row["rank_ols_proportion"] = fit["proportion"]
                pm_draws_local = squamous_boot["pm"][outcome]["CLDN4"]
                plo, phi, _ = percentile_ci(pm_draws_local)
                row["rank_ols_lo"] = plo
                row["rank_ols_hi"] = phi
    return pd.DataFrame(rows)


def lookup(frame, **kwargs):
    mask = np.ones(len(frame), dtype=bool)
    for key, value in kwargs.items():
        mask &= frame[key].to_numpy() == value
    hit = frame.loc[mask]
    if len(hit) != 1:
        raise KeyError(f"expected one row for {kwargs}, found {len(hit)}")
    return hit.iloc[0]


def plot_rho_forest(path, partial_table, pooled, outcome, ylabel_ns):
    cohorts_top_first = ["lung_funnel8", "keratin_funnel7"] + COHORTS
    labels = []
    for key in cohorts_top_first:
        if key in POOL_LABEL:
            labels.append(POOL_LABEL[key])
        else:
            labels.append(f"{key} (n={ylabel_ns[key]})")
    models = ["keratin", "keratin_CLDN4", "keratin_CLDN7", "keratin_EPCAM"]
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    offsets = np.linspace(-0.28, 0.28, len(models))
    for off, model in zip(offsets, models):
        xs, los, his, ys = [], [], [], []
        for i, key in enumerate(cohorts_top_first):
            if key in POOLS:
                row = lookup(pooled, pool=key, outcome=outcome, model=model)
            else:
                row = lookup(partial_table, cohort=key, outcome=outcome, model=model)
            xs.append(row["rho"])
            los.append(max(0.0, float(row["rho"] - row["ci_low"])))
            his.append(max(0.0, float(row["ci_high"] - row["rho"])))
            ys.append(y[i] + off)
        ax.errorbar(
            xs, ys, xerr=[los, his], fmt="o", ms=4.2, lw=1.0, capsize=0,
            color=COLOR[model], label=MODEL_LABEL[model],
        )
    ax.axvline(0, color="0.45", lw=0.8)
    ax.axhline(y[1] - 0.55, color="0.75", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    outcome_label = "CD8 score" if outcome == "CD8_score" else "CD3 score"
    ax.set_xlabel(f"Partial Spearman ρ, TACSTD2 vs {outcome_label}")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_title("Keratin covariates, then one added gene")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_delta_forest(path, cohort_share, share, ns):
    keys = ["lung_funnel8", "keratin_funnel7"] + COHORTS
    labels = []
    for key in keys:
        labels.append(POOL_LABEL[key] if key in POOL_LABEL else f"{key} (n={ns[key]})")
    y = np.arange(len(labels))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 6.6), sharey=True)
    for ax, outcome, title in zip(axes, OUTCOMES, ["CD8 score", "CD3 score"]):
        offsets = np.linspace(-0.22, 0.22, len(SHARE_MODELS))
        for off, model in zip(offsets, SHARE_MODELS):
            for i, key in enumerate(keys):
                if key in POOLS:
                    row = lookup(share, pool=key, outcome=outcome, model=model)
                else:
                    row = lookup(cohort_share, cohort=key, outcome=outcome, model=model)
                yy = y[i] + off
                if np.isfinite(row["delta_lo"]) and np.isfinite(row["delta_hi"]):
                    ax.plot([row["delta_lo"], row["delta_hi"]], [yy, yy], color=COLOR[model], lw=1.15)
                ax.plot(row["delta_rho"], yy, "o", ms=4.0, color=COLOR[model],
                        label=MODEL_LABEL[model] if i == 0 else None)
        ax.axvline(0, color="0.45", lw=0.8)
        ax.axhline(y[1] - 0.55, color="0.75", lw=0.6)
        ax.set_title(title)
        ax.set_xlabel("Δρ = ρ(added) − ρ(keratin)")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels)
    axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle("Change in the TACSTD2 partial correlation", y=1.01)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_share(path, share):
    """Pooled mediation percentages. Clip only the axis, and mark clipped intervals."""
    order = [
        ("keratin_CLDN4", "CLDN4"),
        ("keratin_CLDN7", "CLDN7"),
        ("keratin_EPCAM", "EPCAM"),
        ("keratin_CLDN7_EPCAM_CLDN4", "CLDN4 after CLDN7+EPCAM"),
    ]
    # The last model is not a one-step share of the keratin rho in `share`
    # under that label's usual meaning. Plot the three one-at-a-time shares.
    order = order[:3]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharey=True)
    y = np.arange(len(order))[::-1]
    xlim = [-40, 120]
    for ax, outcome, title in zip(axes, OUTCOMES, ["CD8 score", "CD3 score"]):
        for pool, marker in (("lung_funnel8", "o"), ("keratin_funnel7", "D")):
            color = "#222222" if pool == "lung_funnel8" else "#7a7a7a"
            first = True
            for i, (model, _label) in enumerate(order):
                row = lookup(share, pool=pool, outcome=outcome, model=model)
                yy = y[i] + (-0.12 if pool == "lung_funnel8" else 0.12)
                point = 100 * row["share"]
                lo = 100 * row["share_lo"]
                hi = 100 * row["share_hi"]
                if np.isfinite(lo) and np.isfinite(hi):
                    draw_lo = max(lo, xlim[0])
                    draw_hi = min(hi, xlim[1])
                    if draw_hi >= draw_lo:
                        ax.plot([draw_lo, draw_hi], [yy, yy], color=color, lw=1.2)
                    if lo < xlim[0]:
                        ax.plot(xlim[0] + 1.5, yy, marker="<", color=color, ms=4)
                    if hi > xlim[1]:
                        ax.plot(xlim[1] - 1.5, yy, marker=">", color=color, ms=4)
                ax.plot(point, yy, marker, color=color, ms=5.5,
                        label=POOL_LABEL[pool] if first else None)
                first = False
        ax.axvline(0, color="0.45", lw=0.8)
        ax.set_xlim(*xlim)
        ax.set_title(title)
        ax.set_xlabel("Percent of pooled keratin ρ accounted for")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([label for _model, label in order])
    axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def contains_zero(lo, hi):
    return np.isfinite(lo) and np.isfinite(hi) and lo <= 0 <= hi


def interval_sentence(name, estimate, lo, hi, percent=True):
    shown = fmt_ci(lo, hi, percent=percent)
    point = fmt_pct(estimate) if percent else fmt_num(estimate)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return f"The bootstrap interval for {name} is undefined (point {point})."
    if contains_zero(lo, hi):
        return f"The bootstrap interval for {name} contains zero (point {point}, interval {shown})."
    side = "above zero" if lo > 0 else "below zero"
    return f"The bootstrap interval for {name} lies {side} (point {point}, interval {shown})."


def cohort_sign_sentence(rho, outcome, cohorts):
    negative = [c for c in cohorts if rho[c][outcome]["keratin"] < 0]
    positive = [c for c in cohorts if rho[c][outcome]["keratin"] > 0]
    zero = [c for c in cohorts if rho[c][outcome]["keratin"] == 0]
    return (
        f"Keratin-adjusted TACSTD2–{outcome.replace('_score', '')} ρ is negative in "
        f"{', '.join(negative) if negative else 'no cohort'} "
        f"({len(negative)}/{len(cohorts)}) and positive in "
        f"{', '.join(positive) if positive else 'no cohort'}"
        + (f" and zero in {', '.join(zero)}" if zero else "")
        + "."
    )


def build_report(ctx):
    pooled = ctx["assembled"]["pooled"]
    share = ctx["assembled"]["share"]
    diff = ctx["assembled"]["diff"]
    incremental = ctx["assembled"]["incremental"]
    pm_pool = ctx["assembled"]["pm_pool"]
    partial = ctx["partial"]
    cal = ctx["calibration"]
    path = ctx["path"]
    counts = ctx["counts"]
    mapping = ctx["mapping"]
    squamous = ctx["squamous"]
    rho = ctx["rho"]

    def prow(pool, outcome, model):
        return lookup(pooled, pool=pool, outcome=outcome, model=model)

    def srow(pool, outcome, model):
        return lookup(share, pool=pool, outcome=outcome, model=model)

    cd8_k = prow("lung_funnel8", "CD8_score", "keratin")
    cd8_4 = prow("lung_funnel8", "CD8_score", "keratin_CLDN4")
    sh8 = srow("lung_funnel8", "CD8_score", "keratin_CLDN4")
    sh7c = srow("lung_funnel8", "CD8_score", "keratin_CLDN7")
    she = srow("lung_funnel8", "CD8_score", "keratin_EPCAM")
    cd3_k = prow("lung_funnel8", "CD3_score", "keratin")
    cd3_4 = prow("lung_funnel8", "CD3_score", "keratin_CLDN4")
    sh3 = srow("lung_funnel8", "CD3_score", "keratin_CLDN4")
    f_cd8_k = prow("keratin_funnel7", "CD8_score", "keratin")
    f_cd8_4 = prow("keratin_funnel7", "CD8_score", "keratin_CLDN4")
    fsh = srow("keratin_funnel7", "CD8_score", "keratin_CLDN4")
    inc8 = lookup(incremental, pool="lung_funnel8", outcome="CD8_score")
    d47 = lookup(diff, pool="lung_funnel8", outcome="CD8_score", contrast="CLDN4 minus CLDN7")
    d4e = lookup(diff, pool="lung_funnel8", outcome="CD8_score", contrast="CLDN4 minus EPCAM")
    max_cal = float(cal["abs_diff"].max())
    n8 = int(cd8_k["n_patients"])
    n7 = int(f_cd8_k["n_patients"])

    answer = (
        f"In the eight lung and keratin-funnel cohorts ({n8} primary tumors), "
        f"the DerSimonian–Laird pooled partial Spearman correlation of TACSTD2 with the CD8 score "
        f"after KRT8+KRT18+KRT19 is {fmt_num(cd8_k['rho'])} "
        f"(95% CI {fmt_ci(cd8_k['ci_low'], cd8_k['ci_high'])}, p={fmt_p(cd8_k['p'])}, "
        f"I²={fmt_i2(cd8_k['I2'])}, {int(cd8_k['n_negative'])}/8 cohorts negative). "
        f"Adding CLDN4 changes that pooled correlation to {fmt_num(cd8_4['rho'])} "
        f"(95% CI {fmt_ci(cd8_4['ci_low'], cd8_4['ci_high'])}, p={fmt_p(cd8_4['p'])}, "
        f"I²={fmt_i2(cd8_4['I2'])}, {int(cd8_4['n_negative'])}/8 negative). "
        f"The share accounted for, (ρ_before − ρ_after) / ρ_before, is {fmt_pct(sh8['share'])} "
        f"(patient-bootstrap 95% percentile interval {fmt_ci(sh8['share_lo'], sh8['share_hi'], percent=True)}). "
        f"Adding CLDN7 instead accounts for {fmt_pct(sh7c['share'])} "
        f"(interval {fmt_ci(sh7c['share_lo'], sh7c['share_hi'], percent=True)}). "
        f"Adding EPCAM instead accounts for {fmt_pct(she['share'])} "
        f"(interval {fmt_ci(she['share_lo'], she['share_hi'], percent=True)})."
    )
    answer2 = (
        f"CLDN4 minus CLDN7 is {fmt_pct(d47['share_difference'])} "
        f"(interval {fmt_ci(d47['ci_low'], d47['ci_high'], percent=True)}). "
        f"CLDN4 minus EPCAM is {fmt_pct(d4e['share_difference'])} "
        f"(interval {fmt_ci(d4e['ci_low'], d4e['ci_high'], percent=True)}). "
        f"After CLDN7 and EPCAM are already in the model, adding CLDN4 accounts for a further "
        f"{fmt_pct(inc8['share_of_keratin_rho'])} of the original keratin-adjusted CD8 correlation "
        f"(interval {fmt_ci(inc8['share_lo'], inc8['share_hi'], percent=True)}). "
        f"On the CD3 score the eight-cohort keratin ρ is {fmt_num(cd3_k['rho'])} "
        f"(95% CI {fmt_ci(cd3_k['ci_low'], cd3_k['ci_high'])}, p={fmt_p(cd3_k['p'])}, I²={fmt_i2(cd3_k['I2'])}) "
        f"and the CLDN4 model is {fmt_num(cd3_4['rho'])} "
        f"(95% CI {fmt_ci(cd3_4['ci_low'], cd3_4['ci_high'])}, p={fmt_p(cd3_4['p'])}, I²={fmt_i2(cd3_4['I2'])}), "
        f"a CLDN4 share of {fmt_pct(sh3['share'])} "
        f"(interval {fmt_ci(sh3['share_lo'], sh3['share_hi'], percent=True)})."
    )
    answer3 = (
        f"In the seven-cohort keratin funnel (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD; {n7} tumors) "
        f"the CD8 pooled ρ is {fmt_num(f_cd8_k['rho'])} after keratin "
        f"(95% CI {fmt_ci(f_cd8_k['ci_low'], f_cd8_k['ci_high'])}, p={fmt_p(f_cd8_k['p'])}, "
        f"I²={fmt_i2(f_cd8_k['I2'])}, {int(f_cd8_k['n_negative'])}/7 negative) and "
        f"{fmt_num(f_cd8_4['rho'])} after adding CLDN4 "
        f"(95% CI {fmt_ci(f_cd8_4['ci_low'], f_cd8_4['ci_high'])}, p={fmt_p(f_cd8_4['p'])}, "
        f"I²={fmt_i2(f_cd8_4['I2'])}). "
        f"The CLDN4 share is {fmt_pct(fsh['share'])} "
        f"(interval {fmt_ci(fsh['share_lo'], fsh['share_hi'], percent=True)})."
    )
    cd8_cldn4 = ctx["cohort_share"]
    cd8_cldn4 = cd8_cldn4[(cd8_cldn4["outcome"] == "CD8_score") & (cd8_cldn4["model"] == "keratin_CLDN4")]
    attenuated = []
    strengthened = []
    positive_grew = []
    for _, row in cd8_cldn4.iterrows():
        if row["rho_before"] < 0 and row["delta_rho"] > 0:
            attenuated.append(f"{row['cohort']} ({fmt_num(row['rho_before'])} to {fmt_num(row['rho_after'])})")
        elif row["rho_before"] < 0 and row["delta_rho"] < 0:
            strengthened.append(f"{row['cohort']} ({fmt_num(row['rho_before'])} to {fmt_num(row['rho_after'])})")
        elif row["rho_before"] > 0 and row["delta_rho"] > 0:
            positive_grew.append(f"{row['cohort']} ({fmt_num(row['rho_before'])} to {fmt_num(row['rho_after'])})")
    most_neg = cd8_cldn4.loc[cd8_cldn4["rho_before"].idxmin()]
    answer4 = (
        "That eight-cohort percentage describes the two pooled correlations. "
        f"Adding CLDN4 shrinks a negative CD8 correlation in {', '.join(attenuated) if attenuated else 'no cohort'}. "
        f"It makes a negative correlation more negative in {', '.join(strengthened) if strengthened else 'no cohort'}. "
        f"It makes a positive correlation more positive in {', '.join(positive_grew) if positive_grew else 'no cohort'}. "
        f"The most negative keratin-adjusted correlation is {most_neg['cohort']} "
        f"({fmt_num(most_neg['rho_before'])} to {fmt_num(most_neg['rho_after'])} after CLDN4)."
    )
    sq = ctx["squamous"]
    sq_base = sq[(sq["outcome"] == "CD8_score") & (sq["model"] == "squamous")].iloc[0]
    sq_cldn4 = sq[(sq["outcome"] == "CD8_score") & (sq["model"] == "squamous_CLDN4")].iloc[0]
    answer5 = (
        f"In LUSC only, KRT5+KRT6A+KRT14 in place of KRT8/18/19 gives TACSTD2–CD8 ρ "
        f"{fmt_num(sq_base['rho'])} and a CLDN4 share of {fmt_pct(sq_cldn4['share_of_squamous_rho'])} "
        f"(interval {fmt_ci(sq_cldn4['share_lo'], sq_cldn4['share_hi'], percent=True)}). "
        "That row is outside both pools."
    )
    zh = (
        f"八队列（{n8} 例原发灶）：TACSTD2 与 CD8 分数在 KRT8+KRT18+KRT19 后的合并偏 Spearman ρ 为 "
        f"{fmt_num(cd8_k['rho'])}（95% CI {fmt_ci(cd8_k['ci_low'], cd8_k['ci_high'])}），"
        f"加入 CLDN4 后为 {fmt_num(cd8_4['rho'])}。"
        f"CLDN4 解释的合并相关比例为 {fmt_pct(sh8['share'])}（bootstrap 95% 区间 "
        f"{fmt_ci(sh8['share_lo'], sh8['share_hi'], percent=True)}）。"
        f"CLDN7 为 {fmt_pct(sh7c['share'])}，EPCAM 为 {fmt_pct(she['share'])}。"
        f"CD3 上 CLDN4 的比例为 {fmt_pct(sh3['share'])}。"
        f"七队列角蛋白漏斗（{n7} 例）CD8 的 CLDN4 比例为 {fmt_pct(fsh['share'])}。"
        f"CLDN4 与 CLDN7 的比例之差为 {fmt_pct(d47['share_difference'])}（区间 {fmt_ci(d47['ci_low'], d47['ci_high'], percent=True)}）。"
        f"在 CLDN7 与 EPCAM 已进入模型后，CLDN4 再解释的比例为 {fmt_pct(inc8['share_of_keratin_rho'])} "
        f"（区间 {fmt_ci(inc8['share_lo'], inc8['share_hi'], percent=True)}）。"
        f"这是批量 RNA 上的统计分解。"
    )

    pooled_headers = ["pool", "outcome", "model", "k", "N", "negative cohorts", "pooled ρ", "95% CI", "p", "I²"]
    pooled_rows = []
    for _, row in pooled.iterrows():
        pooled_rows.append([
            row["pool"], row["outcome"], row["model"], int(row["k"]), int(row["n_patients"]),
            f"{int(row['n_negative'])}/{int(row['n_cohorts'])}",
            fmt_num(row["rho"]), fmt_ci(row["ci_low"], row["ci_high"]), fmt_p(row["p"]), fmt_i2(row["I2"]),
        ])

    share_headers = [
        "pool", "outcome", "added", "ρ before", "ρ after", "Δρ", "Δρ bootstrap CI",
        "mediation %", "mediation % bootstrap CI", "fraction of bootstraps with |ρ_before|<0.02",
    ]
    share_md = []
    for _, row in share.iterrows():
        if row["model"] not in SHARE_MODELS and row["model"] != "keratin_CLDN7_EPCAM":
            continue
        share_md.append([
            row["pool"], row["outcome"], row["added"],
            fmt_num(row["rho_before"]), fmt_num(row["rho_after"]), fmt_num(row["delta_rho"]),
            fmt_ci(row["delta_lo"], row["delta_hi"]),
            fmt_pct(row["share"]), fmt_ci(row["share_lo"], row["share_hi"], percent=True),
            fmt_num(row["near_zero_frac"], 3),
        ])

    diff_headers = ["pool", "outcome", "contrast", "mediation % difference", "bootstrap CI"]
    diff_md = []
    for _, row in diff.iterrows():
        diff_md.append([
            row["pool"], row["outcome"], row["contrast"],
            fmt_pct(row["share_difference"]), fmt_ci(row["ci_low"], row["ci_high"], percent=True),
        ])

    inc_headers = [
        "pool", "outcome", "ρ keratin", "ρ +CLDN7+EPCAM", "ρ +CLDN7+EPCAM+CLDN4",
        "further share from CLDN4", "bootstrap CI", "rank-OLS incremental proportion", "rank-OLS CI",
    ]
    inc_md = []
    for _, row in incremental.iterrows():
        inc_md.append([
            row["pool"], row["outcome"],
            fmt_num(row["rho_keratin"]), fmt_num(row["rho_controls"]), fmt_num(row["rho_controls_plus_CLDN4"]),
            fmt_pct(row["share_of_keratin_rho"]), fmt_ci(row["share_lo"], row["share_hi"], percent=True),
            fmt_pct(row["rank_ols_proportion"]), fmt_ci(row["rank_ols_lo"], row["rank_ols_hi"], percent=True),
        ])

    pm_headers = ["pool", "outcome", "mediator", "n-weighted rank-OLS proportion", "bootstrap CI"]
    pm_md = []
    for _, row in pm_pool.iterrows():
        pm_md.append([
            row["pool"], row["outcome"], row["mediator"],
            fmt_pct(row["proportion"]), fmt_ci(row["ci_low"], row["ci_high"], percent=True),
        ])

    def cohort_rho_md(outcome):
        headers = ["cohort", "n", "ρ keratin", "p", "q", "ρ +CLDN4", "p", "q", "ρ +CLDN7", "ρ +EPCAM"]
        rows = []
        for cohort in COHORTS:
            cells = [cohort, str(int(lookup(partial, cohort=cohort, outcome=outcome, model="keratin")["n"]))]
            for model in ["keratin", "keratin_CLDN4", "keratin_CLDN7", "keratin_EPCAM"]:
                row = lookup(partial, cohort=cohort, outcome=outcome, model=model)
                cells.append(fmt_num(row["rho"]))
                if model in ("keratin", "keratin_CLDN4"):
                    cells.extend([fmt_p(row["p"]), fmt_p(row["q"])])
            rows.append(cells)
        return headers, rows

    def cohort_share_md(outcome):
        headers = ["cohort", "added", "ρ before", "ρ after", "Δρ", "Δρ CI", "mediation %", "mediation % CI", "small ρ"]
        rows = []
        sub = ctx["cohort_share"]
        sub = sub[(sub["outcome"] == outcome) & (sub["model"].isin(SHARE_MODELS))]
        for _, row in sub.iterrows():
            rows.append([
                row["cohort"], row["added"], fmt_num(row["rho_before"]), fmt_num(row["rho_after"]),
                fmt_num(row["delta_rho"]), fmt_ci(row["delta_lo"], row["delta_hi"]),
                fmt_pct(row["share"]), fmt_ci(row["share_lo"], row["share_hi"], percent=True),
                "yes" if row["small_denominator"] else "",
            ])
        return headers, rows

    # Path summaries: pool a and b with the same meta function, printed from stored cohort path rows.
    path_pool_rows = []
    ns = ctx["assembled"]["ns"]
    for mediator in MEDIATORS:
        for pool_name, cohorts in POOLS.items():
            block = path[(path["path"] == "a") & (path["mediator"] == mediator)]
            ordered_a = [float(block.loc[block["cohort"] == c, "rho"].iloc[0]) for c in cohorts]
            meta = random_effects_meta(ordered_a, [ns[c] for c in cohorts], 3)
            path_pool_rows.append([
                pool_name, "a: TACSTD2–" + mediator + " | keratin",
                fmt_num(meta["rho"]), fmt_ci(meta["ci_low"], meta["ci_high"]), fmt_p(meta["p"]), fmt_i2(meta["I2"]),
                f"{int(sum(r > 0 for r in ordered_a))}/{len(cohorts)} positive",
            ])
            for outcome in OUTCOMES:
                block_b = path[
                    (path["path"] == "b") & (path["mediator"] == mediator)
                    & (path["outcome"] == outcome) & (path["cohort"].isin(cohorts))
                ]
                # Align to cohort order.
                ordered = [float(block_b.loc[block_b["cohort"] == c, "rho"].iloc[0]) for c in cohorts]
                meta_b = random_effects_meta(ordered, [ns[c] for c in cohorts], 4)
                label = "CD8" if outcome == "CD8_score" else "CD3"
                path_pool_rows.append([
                    pool_name, f"b: {mediator}–{label} | keratin+TACSTD2",
                    fmt_num(meta_b["rho"]), fmt_ci(meta_b["ci_low"], meta_b["ci_high"]),
                    fmt_p(meta_b["p"]), fmt_i2(meta_b["I2"]),
                    f"{int(sum(r < 0 for r in ordered))}/{len(cohorts)} negative",
                ])

    path_headers = ["cohort", "mediator", "a ρ", "a p", "b ρ CD8", "b p CD8", "b ρ CD3", "b p CD3"]
    path_md = []
    for cohort in COHORTS:
        for mediator in MEDIATORS:
            a = path[(path["cohort"] == cohort) & (path["path"] == "a") & (path["mediator"] == mediator)].iloc[0]
            b8 = path[(path["cohort"] == cohort) & (path["path"] == "b") & (path["mediator"] == mediator) & (path["outcome"] == "CD8_score")].iloc[0]
            b3 = path[(path["cohort"] == cohort) & (path["path"] == "b") & (path["mediator"] == mediator) & (path["outcome"] == "CD3_score")].iloc[0]
            path_md.append([
                cohort, mediator, fmt_num(a["rho"]), fmt_p(a["p"]),
                fmt_num(b8["rho"]), fmt_p(b8["p"]), fmt_num(b3["rho"]), fmt_p(b3["p"]),
            ])

    sq_headers = ["outcome", "model", "n", "ρ", "p", "95% CI", "share of squamous ρ", "share CI", "rank-OLS CLDN4 proportion"]
    sq_md = []
    for _, row in squamous.iterrows():
        sq_md.append([
            row["outcome"], row["model"], int(row["n"]), fmt_num(row["rho"]), fmt_p(row["p"]),
            fmt_ci(row["ci_low"], row["ci_high"]),
            fmt_pct(row["share_of_squamous_rho"]) if np.isfinite(row["share_of_squamous_rho"]) else "",
            fmt_ci(row["share_lo"], row["share_hi"], percent=True) if "share_lo" in row and np.isfinite(row.get("share_lo", np.nan)) else "",
            fmt_pct(row["rank_ols_proportion"]) if "rank_ols_proportion" in row and np.isfinite(row.get("rank_ols_proportion", np.nan)) else "",
        ])

    cal_headers = ["cohort", "n", "prior n", "ρ", "prior ρ (3 d.p.)", "|difference|"]
    cal_md = []
    for _, row in cal.iterrows():
        cal_md.append([
            row["cohort"], int(row["n_complete"]), int(row["prior_n"]),
            fmt_num(row["rho"], 4), fmt_num(row["prior_rho_3dp"], 3), f"{row['abs_diff']:.4f}",
        ])

    cd8_h, cd8_rows = cohort_rho_md("CD8_score")
    cd3_h, cd3_rows = cohort_rho_md("CD3_score")
    csh, csrows = cohort_share_md("CD8_score")
    c3h, c3rows = cohort_share_md("CD3_score")

    intervals = [
        interval_sentence("the eight-cohort CD8 CLDN4 mediation percentage", sh8["share"], sh8["share_lo"], sh8["share_hi"]),
        interval_sentence("the eight-cohort CD8 CLDN7 mediation percentage", sh7c["share"], sh7c["share_lo"], sh7c["share_hi"]),
        interval_sentence("the eight-cohort CD8 EPCAM mediation percentage", she["share"], she["share_lo"], she["share_hi"]),
        interval_sentence("the CD8 CLDN4−CLDN7 difference", d47["share_difference"], d47["ci_low"], d47["ci_high"]),
        interval_sentence("the CD8 CLDN4−EPCAM difference", d4e["share_difference"], d4e["ci_low"], d4e["ci_high"]),
        interval_sentence(
            "the further CD8 share from CLDN4 after CLDN7 and EPCAM",
            inc8["share_of_keratin_rho"], inc8["share_lo"], inc8["share_hi"],
        ),
        interval_sentence("the eight-cohort CD3 CLDN4 mediation percentage", sh3["share"], sh3["share_lo"], sh3["share_hi"]),
        interval_sentence("the seven-cohort CD8 CLDN4 mediation percentage", fsh["share"], fsh["share_lo"], fsh["share_hi"]),
    ]
    if cd8_k["I2"] >= 0.5:
        intervals.append(
            f"I² for the eight-cohort keratin-adjusted CD8 correlation is {fmt_i2(cd8_k['I2'])}, "
            f"so that pooled ρ averages cohorts that do not share one effect size."
        )
    intervals.append(cohort_sign_sentence(rho, "CD8_score", COHORTS))
    intervals.append(cohort_sign_sentence(rho, "CD3_score", COHORTS))
    if sh8["near_zero_frac"] >= 0.05:
        intervals.append(
            f"In {fmt_pct(sh8['near_zero_frac'])} of eight-cohort CD8 bootstrap draws the keratin pooled ρ "
            f"has absolute value below 0.02, so the mediation-percentage interval is wide."
        )
    if max_cal > 0.0008:
        cal_sentence = (
            f"The largest absolute difference from the earlier three-decimal TACSTD2–CD8 keratin correlations "
            f"is {max_cal:.4f}. Those earlier figures are a reproduction check; they are not inputs."
        )
    else:
        cal_sentence = (
            f"The keratin-only TACSTD2–CD8 correlations match the earlier three-decimal values "
            f"within {max_cal:.4f}. That check does not use those values as inputs."
        )

    gene_rows = [[sym, mapping[sym]] for sym in GENE_SYMBOLS]
    count_rows = []
    for _, row in counts.iterrows():
        count_rows.append([
            row["cohort"], int(row["n_extracted"]), int(row["n_complete"]),
            int(row["n_dropped_nonfinite"]), int(row["n_01_columns"]) if pd.notna(row["n_01_columns"]) else "NA",
        ])

    lines = [
        "# TCGA: share of the keratin-adjusted TACSTD2–immune correlation accounted for by CLDN4",
        "",
        "Numbers in this file are written by `scripts/tcga_trop2_cldn4_mediation/run_analysis.py`.",
        "",
        "## Answer",
        "",
        answer,
        "",
        answer2,
        "",
        answer3,
        "",
        answer4,
        "",
        answer5,
        "",
        "A positive mediation percentage means the later pooled correlation is closer to zero than the keratin-only correlation. "
        f"The eight-cohort keratin-only CD8 correlation is {'negative' if cd8_k['rho'] < 0 else 'positive' if cd8_k['rho'] > 0 else 'zero'}.",
        "",
        "## 中文摘要",
        "",
        zh,
        "",
        "## Design",
        "",
        "The expression values are UCSC Xena GDC STAR log2(TPM+1), GENCODE v36, primary solid tumor (sample type 01), with replicate aliquots averaged to the patient. "
        "The eight-cohort pool is LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD. "
        "The seven-cohort pool is the locked keratin funnel: LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD. "
        "Keratin covariates are KRT8, KRT18, and KRT19. The CD8 score is the mean of CD8A and CD8B. The CD3 score is the mean of CD3D, CD3E, and CD3G.",
        "",
        "Partial Spearman is the Pearson correlation of average-rank residuals after regression on an intercept and the rank-transformed covariates. "
        "The t test uses df = n − 2 − k. Fisher intervals use variance 1/(n − 3 − k). "
        "Cohorts are combined by DerSimonian–Laird random effects on Fisher z.",
        "",
        "Mediation percentage = 100 × (ρ_keratin − ρ_keratin+covariate) / ρ_keratin, using the two pooled correlations. "
        "The bootstrap resamples patients inside each cohort, refits every cohort, re-pools, and takes the 2.5 and 97.5 percentiles. "
        f"Draws: {ctx['n_boot']}. Master seed: {ctx['seed']}. Cohort c uses seed {ctx['seed']}+1000×(index+1).",
        "",
        "The rank-OLS proportion is a×b/c on the same average ranks, with keratins in every equation, then an n-weighted mean across cohorts. "
        "On this scale c − c′ = a×b. CLDN7 and EPCAM are fit the same way, one mediator at a time. "
        "The incremental row puts CLDN7 and EPCAM in the covariate set and then adds CLDN4. "
        "Its further share is (ρ_controls − ρ_controls+CLDN4) / ρ_keratin.",
        "",
        "Benjamini–Hochberg q-values are computed inside each pre-specified family: keratin tests for one outcome; keratin+CLDN4 tests for one outcome; the two negative-control tests for one outcome. Families are not mixed.",
        "",
        "Cohort rows with |keratin ρ| < 0.05 are marked `small ρ` because that cohort's percentage divides by a small number. Those cohorts stay in every pool.",
        "",
        "## Pooled correlations",
        "",
        md_table(pooled_headers, pooled_rows),
        "",
        "![TACSTD2 vs CD8 partial correlation](figures/fig_cd8_partial_rho.png)",
        "",
        "![TACSTD2 vs CD3 partial correlation](figures/fig_cd3_partial_rho.png)",
        "",
        "Fisher-z intervals are drawn for each cohort. Pooled rows use the DerSimonian–Laird interval.",
        "",
        "## Mediation percentage",
        "",
        md_table(share_headers, share_md),
        "",
        "Δρ = ρ(after adding the covariate) − ρ(keratin). A positive Δρ moves the correlation upward. "
        "When the keratin correlation is negative, a positive Δρ moves it toward zero.",
        "",
        "![Pooled mediation percentage](figures/fig_mediation_share.png)",
        "",
        "The mediation-percentage axis is fixed at −40% to 120% so the three covariates can be compared. An arrow marks an interval that continues past the axis. The table holds the unclipped interval.",
        "",
        "![Change in partial correlation](figures/fig_delta_rho.png)",
        "",
        "### CLDN4 compared with the controls",
        "",
        md_table(diff_headers, diff_md),
        "",
        "### CLDN4 after CLDN7 and EPCAM are in the model",
        "",
        md_table(inc_headers, inc_md),
        "",
        "### Rank-OLS proportion mediated",
        "",
        "On average ranks, the cohort proportion is a×b/c. The n-weighted mean of those proportions is a different summary from the pooled-ρ percentage. "
        "A cohort with a small total coefficient, PAAD in particular, can dominate the weighted mean. The cohort rows are the ones to read.",
        "",
        md_table(pm_headers, pm_md),
        "",
        md_table(
            ["cohort", "outcome", "mediator", "n", "proportion", "total coefficient c"],
            [
                [
                    lookup(ctx["pm_table"], cohort=cohort, outcome=outcome, mediator=mediator)["cohort"],
                    outcome,
                    mediator,
                    int(lookup(ctx["pm_table"], cohort=cohort, outcome=outcome, mediator=mediator)["n"]),
                    fmt_pct(lookup(ctx["pm_table"], cohort=cohort, outcome=outcome, mediator=mediator)["proportion"]),
                    fmt_num(lookup(ctx["pm_table"], cohort=cohort, outcome=outcome, mediator=mediator)["c"], 3),
                ]
                for outcome in OUTCOMES
                for cohort in COHORTS
                for mediator in MEDIATORS
            ],
        ),
        "",
        "## Cohort correlations",
        "",
        "### CD8",
        "",
        md_table(cd8_h, cd8_rows),
        "",
        md_table(csh, csrows),
        "",
        "### CD3",
        "",
        md_table(cd3_h, cd3_rows),
        "",
        md_table(c3h, c3rows),
        "",
        "## Paths",
        "",
        "Path a is the partial Spearman correlation of TACSTD2 with the candidate mediator after KRT8/18/19. "
        "Path b is the partial Spearman correlation of the mediator with the immune score after KRT8/18/19 and TACSTD2.",
        "",
        md_table(["pool", "path", "pooled ρ", "95% CI", "p", "I²", "sign count"], path_pool_rows),
        "",
        md_table(path_headers, path_md),
        "",
        "## LUSC with squamous keratins",
        "",
        "KRT5+KRT6A+KRT14 replace KRT8/18/19 in LUSC only. This row is not entered into either pool.",
        "",
        md_table(sq_headers, sq_md),
        "",
        "## Reproduction check",
        "",
        "The earlier keratin script reported TACSTD2–CD8 partial correlations at three decimals. "
        "The comparison below is an audit. Those printed values are not covariates, weights, or filters.",
        "",
        cal_sentence,
        "",
        md_table(cal_headers, cal_md),
        "",
        "## Intervals",
        "",
        "\n".join(intervals),
        "",
        "## What the percentage is",
        "",
        "The percentage is an accounting split of a bulk-tumor partial correlation. "
        "Patients were not randomized to TACSTD2 or CLDN4. "
        "The samples are untreated primary tumors in TCGA, so the split does not estimate a change after checkpoint blockade. "
        "A keratin residual is a covariate adjustment for KRT8, KRT18, and KRT19, and the LUSC squamous row is the corresponding adjustment for KRT5, KRT6A, and KRT14. "
        "The split does not locate immune cells relative to tumor cells.",
        "",
        "The surface-gene screen that ranked CLDN4 against other membrane genes is a different analysis and is not rerun here.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 scripts/tcga_trop2_cldn4_mediation/test_stats.py",
        "python3 scripts/tcga_trop2_cldn4_mediation/run_analysis.py",
        "```",
        "",
        "The script streams each cohort from the Xena GDC hub, caches the patient-level gene table under `TCGA_MED_CACHE` (default `/tmp/tcga_trop2_cldn4_mediation`), and rewrites this file, the TSV tables, and the figures. "
        f"This run used {ctx['n_boot']} bootstrap draws and master seed {ctx['seed']}.",
        "",
        "## Gene IDs",
        "",
        md_table(["symbol", "ensembl"], gene_rows),
        "",
        "## Sample counts",
        "",
        md_table(["cohort", "extracted patients", "complete", "dropped", "type-01 columns"], count_rows),
        "",
    ]
    return "\n".join(lines) + "\n"


def headline_dict(ctx):
    share = ctx["assembled"]["share"]
    pooled = ctx["assembled"]["pooled"]
    incremental = ctx["assembled"]["incremental"]
    diff = ctx["assembled"]["diff"]
    out = {"n_boot": ctx["n_boot"], "seed": ctx["seed"], "calibration_max_abs_diff": float(ctx["calibration"]["abs_diff"].max())}
    for pool in POOLS:
        for outcome in OUTCOMES:
            for model in ["keratin", "keratin_CLDN4", "keratin_CLDN7", "keratin_EPCAM"]:
                row = lookup(pooled, pool=pool, outcome=outcome, model=model)
                key = f"{pool}|{outcome}|{model}"
                out[key] = {
                    "rho": row["rho"], "ci_low": row["ci_low"], "ci_high": row["ci_high"],
                    "p": row["p"], "I2": row["I2"], "n_negative": int(row["n_negative"]),
                    "n_patients": int(row["n_patients"]),
                }
            for model in SHARE_MODELS:
                row = lookup(share, pool=pool, outcome=outcome, model=model)
                out[f"{pool}|{outcome}|share|{row['added']}"] = {
                    "share": row["share"], "share_lo": row["share_lo"], "share_hi": row["share_hi"],
                    "delta_rho": row["delta_rho"], "rho_before": row["rho_before"], "rho_after": row["rho_after"],
                }
            inc = lookup(incremental, pool=pool, outcome=outcome)
            out[f"{pool}|{outcome}|incremental_CLDN4"] = {
                "share": inc["share_of_keratin_rho"], "share_lo": inc["share_lo"], "share_hi": inc["share_hi"],
            }
        for contrast in diff.loc[diff["pool"] == pool, "contrast"].unique():
            for outcome in OUTCOMES:
                row = lookup(diff, pool=pool, outcome=outcome, contrast=contrast)
                out[f"{pool}|{outcome}|diff|{contrast}"] = {
                    "difference": row["share_difference"], "ci_low": row["ci_low"], "ci_high": row["ci_high"],
                }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--boot", type=int, default=int(os.environ.get("TCGA_MED_BOOT", "2000")))
    parser.add_argument("--seed", type=int, default=20260921)
    args = parser.parse_args()
    if args.boot < 50:
        print(f"[warn] bootstrap count {args.boot} is below 50; percentile intervals will be marked undefined", flush=True)
    os.makedirs(TABLES, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    print("[map] probemap", flush=True)
    mapping = load_probemap(CACHE, GENE_SYMBOLS)
    matrices, counts, provenance = load_cohorts(mapping)
    for cohort in COHORTS:
        n = int(counts.loc[counts["cohort"] == cohort, "n_complete"].iloc[0])
        if n != PRIOR_N[cohort]:
            print(f"[warn] {cohort} complete n={n}, earlier extract n={PRIOR_N[cohort]}", flush=True)
    partial, rho = point_partials(matrices)
    path, pm_table, pm = point_paths_and_pm(matrices)
    boots, squamous_boot = run_bootstraps(matrices, args.boot, args.seed)
    worst = agree_with_points(boots, rho, pm)
    print(f"[agree] max |full-sample gap| {worst:.3e}", flush=True)
    assembled = assemble(partial, rho, pm, boots, counts)
    cohort_share = cohort_share_table(rho, boots, assembled["ns"])
    cal = calibration_table(rho, counts)
    squamous = squamous_table(matrices, squamous_boot)
    print(f"[calibration] max |diff| vs prior 3 d.p. CD8 rho = {cal['abs_diff'].max():.6f}", flush=True)

    partial.to_csv(os.path.join(TABLES, "partial_rho_by_cohort.tsv"), sep="\t", index=False, float_format="%.8g")
    assembled["pooled"].to_csv(os.path.join(TABLES, "pooled_rho.tsv"), sep="\t", index=False, float_format="%.8g")
    assembled["share"].to_csv(os.path.join(TABLES, "mediation_share_pooled.tsv"), sep="\t", index=False, float_format="%.8g")
    assembled["diff"].to_csv(os.path.join(TABLES, "mediation_share_differences.tsv"), sep="\t", index=False, float_format="%.8g")
    assembled["incremental"].to_csv(os.path.join(TABLES, "incremental_cldn4.tsv"), sep="\t", index=False, float_format="%.8g")
    assembled["pm_pool"].to_csv(os.path.join(TABLES, "rank_ols_pooled.tsv"), sep="\t", index=False, float_format="%.8g")
    cohort_share.to_csv(os.path.join(TABLES, "mediation_share_by_cohort.tsv"), sep="\t", index=False, float_format="%.8g")
    path.to_csv(os.path.join(TABLES, "paths_by_cohort.tsv"), sep="\t", index=False, float_format="%.8g")
    pm_table.to_csv(os.path.join(TABLES, "rank_ols_by_cohort.tsv"), sep="\t", index=False, float_format="%.8g")
    squamous.to_csv(os.path.join(TABLES, "lusc_squamous.tsv"), sep="\t", index=False, float_format="%.8g")
    cal.to_csv(os.path.join(TABLES, "calibration_cd8_keratin.tsv"), sep="\t", index=False, float_format="%.8g")
    counts.to_csv(os.path.join(TABLES, "sample_counts.tsv"), sep="\t", index=False)

    plot_rho_forest(os.path.join(FIGS, "fig_cd8_partial_rho.png"), partial, assembled["pooled"], "CD8_score", assembled["ns"])
    plot_rho_forest(os.path.join(FIGS, "fig_cd3_partial_rho.png"), partial, assembled["pooled"], "CD3_score", assembled["ns"])
    plot_delta_forest(os.path.join(FIGS, "fig_delta_rho.png"), cohort_share, assembled["share"], assembled["ns"])
    plot_share(os.path.join(FIGS, "fig_mediation_share.png"), assembled["share"])

    ctx = {
        "assembled": assembled,
        "partial": partial,
        "calibration": cal,
        "path": path,
        "counts": counts,
        "mapping": mapping,
        "squamous": squamous,
        "rho": rho,
        "cohort_share": cohort_share,
        "pm_table": pm_table,
        "n_boot": args.boot,
        "seed": args.seed,
    }
    report = build_report(ctx)
    with open(os.path.join(OUT, "RESULTS.md"), "w", encoding="utf-8") as handle:
        handle.write(report)
    head = headline_dict(ctx)
    with open(os.path.join(TABLES, "headline.json"), "w", encoding="utf-8") as handle:
        json.dump(head, handle, indent=2, sort_keys=True)
        handle.write("\n")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    prov = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_boot": args.boot,
        "seed": args.seed,
        "cohort_seed": {cohort: args.seed + 1000 * (i + 1) for i, cohort in enumerate(COHORTS)},
        "lusc_squamous_seed": args.seed + 9000,
        "gene_ids": mapping,
        "agreement_max_abs_gap": worst,
        "calibration_max_abs_diff": float(cal["abs_diff"].max()),
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": stats.__version__ if hasattr(stats, "__version__") else "scipy",
        "script_sha256": {
            name: file_sha256(os.path.join(script_dir, name))
            for name in ["run_analysis.py", "stats.py", "bootstrap.py", "io_tcga.py"]
        },
        "extracts": provenance,
        "definition": {
            "cd8_score": "mean(CD8A, CD8B)",
            "cd3_score": "mean(CD3D, CD3E, CD3G)",
            "keratin": KERATIN,
            "mediation_percent": "(rho_keratin - rho_added) / rho_keratin",
            "pools": {name: cohorts for name, cohorts in POOLS.items()},
        },
    }
    # scipy version
    import scipy
    prov["scipy"] = scipy.__version__
    with open(os.path.join(TABLES, "provenance.json"), "w", encoding="utf-8") as handle:
        json.dump(prov, handle, indent=2)
        handle.write("\n")
    print(f"[done] {os.path.join(OUT, 'RESULTS.md')}", flush=True)
    if float(cal["abs_diff"].max()) > 0.0008:
        print("[fail] keratin-only CD8 correlations do not reproduce the earlier three-decimal values", flush=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
