#!/usr/bin/env python3
"""CPTAC LUAD/LSCC protein: TROP2 → CLDN4 → CD8A/MHC-I, with EPCAM as control.

Public TMT freeze v1.2, tumor protein only. LUAD and LSCC are not pooled.
No log2 value is imputed. CLDN4 is used only where it was quantified.

Estimands are fixed in this file before the cohort fits are read:

- Exposure X: TACSTD2 protein (TROP2).
- Mediator under test M1: CLDN4 protein.
- Control protein M2: EPCAM protein, entered in the same equations.
- Outcomes Y: CD8A protein, and the MHC-I score (mean of HLA-A/B/C z-scores,
  the same score as the earlier protein page on this freeze).
- Primary sample: tumors with X, CLDN4, EPCAM, and Y all quantified.
- Single-mediator model: z-rank OLS, M ~ X and Y ~ X + M.
  Indirect effect ab = a*b. Total effect c from Y ~ X. Direct effect c'.
  On this scale a and c equal Spearman ρ.
- Parallel model: M1 ~ X, M2 ~ X, Y ~ X + M1 + M2.
  Indirect effects a1*b1 and a2*b2. Algebra: c = c' + a1*b1 + a2*b2.
- Uncertainty: percentile bootstrap, 2,000 resamples, ranks recomputed
  inside each resample. Seed 20260921, one stream per estimand label.
- Partial Spearman: Pearson correlation of rank residuals. This is the
  "partial CLDN4" / "partial EPCAM" estimand. It is not the regression
  coefficient b.
- Primary q: Benjamini–Hochberg on the 16 indirect-effect bootstrap
  p-values (2 cohorts × 2 outcomes × CLDN4/EPCAM × single/parallel).
- Partial-family q: BH on the partial Spearman tests listed in PARTIAL_SPECS.
- WES purity enters a sensitivity model as an extra rank covariate in every
  equation. It is not part of the primary call.

Mediation rule (parallel model, each cohort × outcome). All five must hold:

1. a for CLDN4 has bootstrap CI entirely above 0 (TROP2 and CLDN4 move together).
2. b for CLDN4 has bootstrap CI entirely below 0.
3. CLDN4 indirect effect ab has bootstrap CI entirely below 0.
4. Total effect c has bootstrap CI entirely below 0.
5. EPCAM indirect effect in the same parallel model does not have its CI
   entirely below 0.

Meeting the rule is an associational statement on a treatment-naive
resection proteome. It is not an experiment and not an ICI result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genes import (  # noqa: E402
    CONTROL,
    ENSEMBL,
    EXPOSURE,
    MEDIATOR,
    MHC1_GENES,
    OUTCOMES,
    REQUIRED,
)

SEED = 20260921
N_BOOT = 2000
MIN_N = 20
MIN_BOOT = 100
EXPECTED_N = {"LUAD": 110, "LSCC": 108}

# Point estimates from methods/cptac_ifn_mhc_protein/tables/summary.json
# on the same freeze files. Used only as a reproduction check.
PRIOR = [
    ("LUAD", "CLDN4", "CD8A", 79, -0.09790652385589094, 0.3906679722862911),
    ("LUAD", "CLDN4", "MHC1", 79, 0.04946445959104188, 0.6650821699051745),
    ("LUAD", "TACSTD2", "CD8A", 110, -0.03836079615896129, 0.6907172403595156),
    ("LUAD", "TACSTD2", "MHC1", 110, -0.08329012916168879, 0.38699507924129223),
    ("LSCC", "CLDN4", "CD8A", 78, -0.4441002035938745, 4.64372740784413e-05),
    ("LSCC", "CLDN4", "MHC1", 78, -0.4097042198308021, 0.0001953844621326853),
    ("LSCC", "TACSTD2", "CD8A", 108, -0.10335629293015898, 0.28711792033709344),
    ("LSCC", "TACSTD2", "MHC1", 108, 0.0291996532243467, 0.7641888603082562),
]

PARTIAL_SPECS = [
    ("TACSTD2", "CD8A", ("CLDN4",)),
    ("TACSTD2", "CD8A", ("EPCAM",)),
    ("TACSTD2", "CD8A", ("CLDN4", "EPCAM")),
    ("TACSTD2", "MHC1", ("CLDN4",)),
    ("TACSTD2", "MHC1", ("EPCAM",)),
    ("TACSTD2", "MHC1", ("CLDN4", "EPCAM")),
    ("CLDN4", "CD8A", ("TACSTD2",)),
    ("CLDN4", "CD8A", ("EPCAM",)),
    ("CLDN4", "CD8A", ("TACSTD2", "EPCAM")),
    ("CLDN4", "MHC1", ("TACSTD2",)),
    ("CLDN4", "MHC1", ("EPCAM",)),
    ("CLDN4", "MHC1", ("TACSTD2", "EPCAM")),
    ("EPCAM", "CD8A", ("CLDN4",)),
    ("EPCAM", "CD8A", ("TACSTD2",)),
    ("EPCAM", "CD8A", ("CLDN4", "TACSTD2")),
    ("EPCAM", "MHC1", ("CLDN4",)),
    ("EPCAM", "MHC1", ("TACSTD2",)),
    ("EPCAM", "MHC1", ("CLDN4", "TACSTD2")),
]

ASSOC_PAIRS = [
    ("TACSTD2", "CD8A"),
    ("TACSTD2", "MHC1"),
    ("CLDN4", "CD8A"),
    ("CLDN4", "MHC1"),
    ("EPCAM", "CD8A"),
    ("EPCAM", "MHC1"),
    ("TACSTD2", "CLDN4"),
    ("TACSTD2", "EPCAM"),
    ("CLDN4", "EPCAM"),
]

TREAT_KEYS = (
    "treat",
    "therap",
    "neoadj",
    "ici",
    "immunotherap",
    "pd-1",
    "pd1",
    "pembro",
    "nivol",
    "atezo",
    "chemo",
    "drug",
)


def rng_for(*labels: str) -> np.random.Generator:
    payload = "\n".join(labels).encode()
    digest = hashlib.sha256(payload).digest()
    extra = int.from_bytes(digest[:8], "little") % (2**32 - 1)
    return np.random.default_rng(np.random.SeedSequence([SEED, extra]))


def zrank(v: np.ndarray) -> np.ndarray:
    r = stats.rankdata(v, method="average").astype(float)
    sd = float(r.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        raise ValueError("constant")
    return (r - r.mean()) / sd


def ols(y: np.ndarray, x_cols: list[np.ndarray]) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), *x_cols])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return coef[1:]


def resid(y: np.ndarray, cov: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), cov])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def spearman_stat(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    res = stats.spearmanr(x, y)
    rho = float(getattr(res, "statistic", res[0]))
    p = float(getattr(res, "pvalue", res[1]))
    return rho, p


def pearson_stat(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    res = stats.pearsonr(x, y)
    rho = float(getattr(res, "statistic", res[0]))
    p = float(getattr(res, "pvalue", res[1]))
    return rho, p


def partial_spearman(y: np.ndarray, x: np.ndarray, covs: list[np.ndarray]) -> tuple[float, float]:
    ry = stats.rankdata(y, method="average").astype(float)
    rx = stats.rankdata(x, method="average").astype(float)
    rc = np.column_stack([stats.rankdata(c, method="average").astype(float) for c in covs])
    yr = resid(ry, rc)
    xr = resid(rx, rc)
    if float(np.std(xr)) == 0 or float(np.std(yr)) == 0:
        return np.nan, np.nan
    return pearson_stat(xr, yr)


def vif_max(cols: list[np.ndarray]) -> float:
    vals = []
    for i, y in enumerate(cols):
        others = [cols[j] for j in range(len(cols)) if j != i]
        design = np.column_stack([np.ones(len(y)), *others])
        coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        pred = design @ coef
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        ss_res = float(np.sum((y - pred) ** 2))
        if ss_tot == 0:
            return np.inf
        r2 = 1.0 - ss_res / ss_tot
        if r2 >= 1:
            return np.inf
        vals.append(1.0 / (1.0 - r2))
    return float(max(vals))


def fit_shared(x: np.ndarray, m1: np.ndarray, m2: np.ndarray, y: np.ndarray) -> dict[str, float]:
    zx, z1, z2, zy = zrank(x), zrank(m1), zrank(m2), zrank(y)
    a1 = float(ols(z1, [zx])[0])
    a2 = float(ols(z2, [zx])[0])
    c = float(ols(zy, [zx])[0])
    c1, b1 = (float(v) for v in ols(zy, [zx, z1]))
    c2, b2 = (float(v) for v in ols(zy, [zx, z2]))
    cp, bp1, bp2 = (float(v) for v in ols(zy, [zx, z1, z2]))
    out = {
        "a_cldn4": a1,
        "a_epcam": a2,
        "b_cldn4_single": b1,
        "b_epcam_single": b2,
        "c_prime_cldn4_single": c1,
        "c_prime_epcam_single": c2,
        "ab_cldn4_single": a1 * b1,
        "ab_epcam_single": a2 * b2,
        "b_cldn4_parallel": bp1,
        "b_epcam_parallel": bp2,
        "c_prime_parallel": cp,
        "ab_cldn4_parallel": a1 * bp1,
        "ab_epcam_parallel": a2 * bp2,
        "c": c,
        "vif_max": vif_max([zx, z1, z2]),
    }
    if abs(out["ab_cldn4_single"] - (c - c1)) > 1e-8:
        raise AssertionError("single CLDN4 identity failed")
    if abs(out["ab_epcam_single"] - (c - c2)) > 1e-8:
        raise AssertionError("single EPCAM identity failed")
    if abs(out["ab_cldn4_parallel"] + out["ab_epcam_parallel"] - (c - cp)) > 1e-8:
        raise AssertionError("parallel identity failed")
    return out


def fit_single(x: np.ndarray, m: np.ndarray, y: np.ndarray) -> dict[str, float]:
    zx, zm, zy = zrank(x), zrank(m), zrank(y)
    a = float(ols(zm, [zx])[0])
    c = float(ols(zy, [zx])[0])
    cp, b = (float(v) for v in ols(zy, [zx, zm]))
    out = {"a": a, "b": b, "c": c, "c_prime": cp, "ab": a * b}
    if abs(out["ab"] - (c - cp)) > 1e-8:
        raise AssertionError("single identity failed")
    return out


def fit_shared_purity(
    x: np.ndarray, m1: np.ndarray, m2: np.ndarray, y: np.ndarray, pur: np.ndarray
) -> dict[str, float]:
    zx, z1, z2, zy, zp = zrank(x), zrank(m1), zrank(m2), zrank(y), zrank(pur)
    a1 = float(ols(z1, [zx, zp])[0])
    a2 = float(ols(z2, [zx, zp])[0])
    c = float(ols(zy, [zx, zp])[0])
    c1, b1, _g1 = (float(v) for v in ols(zy, [zx, z1, zp]))
    c2, b2, _g2 = (float(v) for v in ols(zy, [zx, z2, zp]))
    cp, bp1, bp2, _gp = (float(v) for v in ols(zy, [zx, z1, z2, zp]))
    out = {
        "a_cldn4": a1,
        "a_epcam": a2,
        "b_cldn4_single": b1,
        "b_epcam_single": b2,
        "c_prime_cldn4_single": c1,
        "c_prime_epcam_single": c2,
        "ab_cldn4_single": a1 * b1,
        "ab_epcam_single": a2 * b2,
        "b_cldn4_parallel": bp1,
        "b_epcam_parallel": bp2,
        "c_prime_parallel": cp,
        "ab_cldn4_parallel": a1 * bp1,
        "ab_epcam_parallel": a2 * bp2,
        "c": c,
    }
    if abs(out["ab_cldn4_single"] - (c - c1)) > 1e-8:
        raise AssertionError("purity single CLDN4 identity failed")
    if abs(out["ab_epcam_single"] - (c - c2)) > 1e-8:
        raise AssertionError("purity single EPCAM identity failed")
    if abs(out["ab_cldn4_parallel"] + out["ab_epcam_parallel"] - (c - cp)) > 1e-8:
        raise AssertionError("purity parallel identity failed")
    return out


def bootstrap(point_fn, n: int, rng: np.random.Generator) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    point = point_fn(np.arange(n))
    keys = list(point)
    store = {k: [] for k in keys}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        try:
            est = point_fn(idx)
        except ValueError:
            continue
        for k in keys:
            v = est[k]
            if np.isfinite(v):
                store[k].append(float(v))
    samples = {k: np.asarray(v, dtype=float) for k, v in store.items()}
    return point, samples


def summarize_boot(point: dict[str, float], samples: dict[str, np.ndarray], keys: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for k in keys:
        arr = samples[k]
        out[k] = float(point[k])
        if len(arr) < MIN_BOOT:
            out[f"{k}_ci_lo"] = np.nan
            out[f"{k}_ci_hi"] = np.nan
            out[f"{k}_p"] = np.nan
            out[f"{k}_p_floor"] = np.nan
            out[f"{k}_nboot"] = float(len(arr))
            continue
        out[f"{k}_ci_lo"] = float(np.percentile(arr, 2.5))
        out[f"{k}_ci_hi"] = float(np.percentile(arr, 97.5))
        p = 2.0 * min(float(np.mean(arr <= 0)), float(np.mean(arr >= 0)))
        p = min(p, 1.0)
        floor = 0.0
        if p == 0:
            p = 1.0 / len(arr)
            floor = 1.0
        out[f"{k}_p"] = p
        out[f"{k}_p_floor"] = floor
        out[f"{k}_nboot"] = float(len(arr))
    return out


def bh_q(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    m = len(idx)
    if m == 0:
        return out.tolist()
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out.tolist()


def ci_side(lo: float, hi: float) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "na"
    if hi < 0:
        return "below"
    if lo > 0:
        return "above"
    return "crosses"


def self_check() -> None:
    rng = np.random.default_rng(1)
    n = 160
    x = rng.normal(size=n)
    x[::7] = x[0]  # ties
    m1 = 0.7 * x + rng.normal(size=n)
    m2 = 0.35 * x + 0.25 * m1 + rng.normal(size=n)
    y = 0.15 * x - 0.55 * m1 + 0.2 * m2 + rng.normal(size=n)
    pur = 0.2 * x + rng.normal(size=n)
    single = fit_single(x, m1, y)
    rho, _p = spearman_stat(x, m1)
    if abs(single["a"] - rho) > 1e-8:
        raise AssertionError(f"a {single['a']} != spearman {rho}")
    rho_c, _ = spearman_stat(x, y)
    if abs(single["c"] - rho_c) > 1e-8:
        raise AssertionError("c != spearman")
    fit_shared(x, m1, m2, y)
    fit_shared_purity(x, m1, m2, y, pur)
    # Partial Spearman is invariant to monotone transforms of one covariate-free pair.
    r1, _ = partial_spearman(y, m1, [x])
    if not np.isfinite(r1):
        raise AssertionError("partial spearman failed")
    print("self-check ok", flush=True)


def find_row(index: pd.Index, prefixes: list[str]) -> str | None:
    hits: list[str] = []
    for prefix in prefixes:
        for i in index:
            s = str(i)
            if s == prefix or s.startswith(prefix + ".") or s.startswith(prefix + "|"):
                hits.append(s)
    hits = list(dict.fromkeys(hits))
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {prefixes}: {hits}")
    return hits[0]


def protein_path(data: Path, cohort: str) -> Path:
    return (
        data
        / cohort
        / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )


def phenotype_path(data: Path, cohort: str) -> Path:
    return data / cohort / f"{cohort}_phenotype.txt"


def mean_z(series_list: list[pd.Series], min_genes: int) -> pd.Series:
    rows = []
    for s in series_list:
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        rows.append((s - s.mean()) / sd)
    if not rows:
        return pd.Series(dtype=float)
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    return score.where(n >= min_genes)


def load_cohort(data: Path, cohort: str) -> dict:
    mat = pd.read_csv(protein_path(data, cohort), sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    mat = mat.apply(pd.to_numeric, errors="coerce")
    ph = pd.read_csv(phenotype_path(data, cohort), sep="\t")
    if "idx" not in ph.columns:
        raise ValueError(f"{cohort} phenotype has no idx")
    ph = ph.set_index(ph["idx"].astype(str))
    if set(mat.columns.astype(str)) != set(ph.index.astype(str)):
        raise ValueError(f"{cohort} protein columns do not match phenotype idx")
    mat.columns = mat.columns.astype(str)
    ph = ph.reindex(mat.columns)
    n_tumors = int(mat.shape[1])
    if n_tumors != EXPECTED_N[cohort]:
        raise ValueError(f"{cohort} n_tumors={n_tumors}, expected {EXPECTED_N[cohort]}")
    treat_cols = [c for c in ph.columns if any(k in str(c).lower() for k in TREAT_KEYS)]
    series: dict[str, pd.Series] = {}
    rows_used: dict[str, str | None] = {}
    coverage = []
    symbols = list(ENSEMBL)
    for symbol in symbols:
        row = find_row(mat.index, ENSEMBL[symbol])
        rows_used[symbol] = row
        if row is None:
            s = pd.Series(np.nan, index=mat.columns, name=symbol)
        else:
            s = pd.to_numeric(mat.loc[row], errors="coerce")
            s.index = s.index.astype(str)
            s.name = symbol
        series[symbol] = s
        n_q = int(s.notna().sum()) if row is not None else 0
        coverage.append(
            {
                "cohort": cohort,
                "symbol": symbol,
                "ensembl_query": ";".join(ENSEMBL[symbol]),
                "row": row or "",
                "present": row is not None,
                "n_tumors": n_tumors,
                "n_quantified": n_q,
                "n_missing": n_tumors - n_q,
                "pct_missing": 100.0 * (n_tumors - n_q) / n_tumors,
            }
        )
    missing_required = [g for g in REQUIRED if rows_used[g] is None]
    if missing_required:
        raise ValueError(f"{cohort} missing required rows: {missing_required}")
    mhc_vecs = []
    mhc_members = []
    for g in MHC1_GENES:
        if int(series[g].notna().sum()) >= MIN_N:
            mhc_members.append(g)
            mhc_vecs.append(series[g])
    if len(mhc_members) < 3:
        raise ValueError(f"{cohort} MHC-I genes quantified: {mhc_members}")
    series["MHC1"] = mean_z(mhc_vecs, 3)
    series["MHC1"].name = "MHC1"
    wes = pd.to_numeric(ph["WES_purity"], errors="coerce")
    wes.index = wes.index.astype(str)
    series["WES_purity"] = wes.reindex(mat.columns)
    return {
        "cohort": cohort,
        "n_tumors": n_tumors,
        "n_pheno_cols": int(ph.shape[1]),
        "treat_cols": treat_cols,
        "series": series,
        "rows_used": rows_used,
        "coverage": coverage,
        "mhc_members": mhc_members,
    }


def align(bundle: dict, names: list[str]) -> pd.DataFrame:
    df = pd.concat([bundle["series"][n].rename(n) for n in names], axis=1)
    return df.dropna()


def boot_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> dict[str, float]:
    def point_fn(idx: np.ndarray) -> dict[str, float]:
        rho, _p = spearman_stat(x[idx], y[idx])
        if not np.isfinite(rho):
            raise ValueError("nonfinite rho")
        return {"rho": float(rho)}

    point, samples = bootstrap(point_fn, len(x), rng)
    rho, p = spearman_stat(x, y)
    sm = summarize_boot(point, samples, ["rho"])
    sm["rho"] = float(rho)
    sm["p_analytic"] = float(p)
    return sm


def mediation_rows_from_summary(
    cohort: str,
    outcome: str,
    n: int,
    sample_set: str,
    summary: dict[str, float],
    model: str,
    mediator: str,
    analytic_a_p: float,
    analytic_c_p: float,
) -> dict:
    if model == "single" and mediator == "CLDN4":
        prefix = "cldn4_single"
        a_key, b_key, cp_key, ab_key = (
            "a_cldn4",
            "b_cldn4_single",
            "c_prime_cldn4_single",
            "ab_cldn4_single",
        )
    elif model == "single" and mediator == "EPCAM":
        a_key, b_key, cp_key, ab_key = (
            "a_epcam",
            "b_epcam_single",
            "c_prime_epcam_single",
            "ab_epcam_single",
        )
    elif model == "parallel" and mediator == "CLDN4":
        a_key, b_key, cp_key, ab_key = (
            "a_cldn4",
            "b_cldn4_parallel",
            "c_prime_parallel",
            "ab_cldn4_parallel",
        )
    elif model == "parallel" and mediator == "EPCAM":
        a_key, b_key, cp_key, ab_key = (
            "a_epcam",
            "b_epcam_parallel",
            "c_prime_parallel",
            "ab_epcam_parallel",
        )
    else:
        raise ValueError((model, mediator))
    c = summary["c"]
    ab = summary[ab_key]
    prop = ab / c if np.isfinite(c) and c != 0 else np.nan
    return {
        "cohort": cohort,
        "outcome": outcome,
        "mediator": mediator,
        "model": model,
        "sample_set": sample_set,
        "n": n,
        "a": summary[a_key],
        "a_ci_lo": summary[f"{a_key}_ci_lo"],
        "a_ci_hi": summary[f"{a_key}_ci_hi"],
        "a_p_analytic": analytic_a_p,
        "b": summary[b_key],
        "b_ci_lo": summary[f"{b_key}_ci_lo"],
        "b_ci_hi": summary[f"{b_key}_ci_hi"],
        "b_p": summary[f"{b_key}_p"],
        "c": c,
        "c_ci_lo": summary["c_ci_lo"],
        "c_ci_hi": summary["c_ci_hi"],
        "c_p_analytic": analytic_c_p,
        "c_prime": summary[cp_key],
        "c_prime_ci_lo": summary[f"{cp_key}_ci_lo"],
        "c_prime_ci_hi": summary[f"{cp_key}_ci_hi"],
        "ab": ab,
        "ab_ci_lo": summary[f"{ab_key}_ci_lo"],
        "ab_ci_hi": summary[f"{ab_key}_ci_hi"],
        "ab_p": summary[f"{ab_key}_p"],
        "ab_p_floor": summary[f"{ab_key}_p_floor"],
        "ab_nboot": summary[f"{ab_key}_nboot"],
        "proportion_ab_over_c": prop,
        "vif_max": summary.get("vif_max", np.nan),
    }


SHARED_KEYS = [
    "a_cldn4",
    "a_epcam",
    "b_cldn4_single",
    "b_epcam_single",
    "c_prime_cldn4_single",
    "c_prime_epcam_single",
    "ab_cldn4_single",
    "ab_epcam_single",
    "b_cldn4_parallel",
    "b_epcam_parallel",
    "c_prime_parallel",
    "ab_cldn4_parallel",
    "ab_epcam_parallel",
    "c",
]


def run_shared(bundle: dict, outcome: str) -> tuple[list[dict], dict]:
    cohort = bundle["cohort"]
    df = align(bundle, ["TACSTD2", "CLDN4", "EPCAM", outcome])
    n = int(len(df))
    if n < MIN_N:
        raise ValueError(f"{cohort} {outcome} shared n={n}")
    x = df["TACSTD2"].to_numpy()
    m1 = df["CLDN4"].to_numpy()
    m2 = df["EPCAM"].to_numpy()
    y = df[outcome].to_numpy()
    rng = rng_for("shared", cohort, outcome)

    def point_fn(idx: np.ndarray) -> dict[str, float]:
        est = fit_shared(x[idx], m1[idx], m2[idx], y[idx])
        est.pop("vif_max")
        return est

    point, samples = bootstrap(point_fn, n, rng)
    point["vif_max"] = fit_shared(x, m1, m2, y)["vif_max"]
    samples["vif_max"] = np.asarray([point["vif_max"]])
    summary = summarize_boot(point, samples, SHARED_KEYS + ["vif_max"])
    summary["vif_max"] = float(point["vif_max"])
    _ra, pa = spearman_stat(x, m1)
    _rb, pb = spearman_stat(x, m2)
    _rc, pc = spearman_stat(x, y)
    if abs(summary["a_cldn4"] - _ra) > 1e-8 or abs(summary["a_epcam"] - _rb) > 1e-8:
        raise AssertionError("shared a diverged from Spearman")
    if abs(summary["c"] - _rc) > 1e-8:
        raise AssertionError("shared c diverged from Spearman")
    rows = []
    for model, mediator, ap in (
        ("single", "CLDN4", pa),
        ("single", "EPCAM", pb),
        ("parallel", "CLDN4", pa),
        ("parallel", "EPCAM", pb),
    ):
        rows.append(
            mediation_rows_from_summary(
                cohort, outcome, n, "shared", summary, model, mediator, ap, pc
            )
        )
    meta = {"n": n, "index": df.index.astype(str).tolist(), "vif_max": float(point["vif_max"])}
    return rows, meta


def run_purity(bundle: dict, outcome: str) -> list[dict]:
    cohort = bundle["cohort"]
    df = align(bundle, ["TACSTD2", "CLDN4", "EPCAM", outcome, "WES_purity"])
    n = int(len(df))
    x = df["TACSTD2"].to_numpy()
    m1 = df["CLDN4"].to_numpy()
    m2 = df["EPCAM"].to_numpy()
    y = df[outcome].to_numpy()
    pur = df["WES_purity"].to_numpy()
    rng = rng_for("purity", cohort, outcome)

    def point_fn(idx: np.ndarray) -> dict[str, float]:
        return fit_shared_purity(x[idx], m1[idx], m2[idx], y[idx], pur[idx])

    point, samples = bootstrap(point_fn, n, rng)
    summary = summarize_boot(point, samples, SHARED_KEYS)
    rows = []
    for model, mediator in (
        ("single", "CLDN4"),
        ("single", "EPCAM"),
        ("parallel", "CLDN4"),
        ("parallel", "EPCAM"),
    ):
        rec = mediation_rows_from_summary(
            cohort, outcome, n, "shared_wes_purity", summary, model, mediator, np.nan, np.nan
        )
        rows.append(rec)
    return rows


def run_mediator_only(bundle: dict, outcome: str, mediator: str) -> dict:
    cohort = bundle["cohort"]
    df = align(bundle, ["TACSTD2", mediator, outcome])
    n = int(len(df))
    x = df["TACSTD2"].to_numpy()
    m = df[mediator].to_numpy()
    y = df[outcome].to_numpy()
    rng = rng_for("mediator_only", cohort, outcome, mediator)

    def point_fn(idx: np.ndarray) -> dict[str, float]:
        est = fit_single(x[idx], m[idx], y[idx])
        return est

    point, samples = bootstrap(point_fn, n, rng)
    summary_s = summarize_boot(point, samples, ["a", "b", "c", "c_prime", "ab"])
    _ra, pa = spearman_stat(x, m)
    _rc, pc = spearman_stat(x, y)
    rec = {
        "cohort": cohort,
        "outcome": outcome,
        "mediator": mediator,
        "model": "single",
        "sample_set": "mediator_complete",
        "n": n,
        "a": summary_s["a"],
        "a_ci_lo": summary_s["a_ci_lo"],
        "a_ci_hi": summary_s["a_ci_hi"],
        "a_p_analytic": pa,
        "b": summary_s["b"],
        "b_ci_lo": summary_s["b_ci_lo"],
        "b_ci_hi": summary_s["b_ci_hi"],
        "b_p": summary_s["b_p"],
        "c": summary_s["c"],
        "c_ci_lo": summary_s["c_ci_lo"],
        "c_ci_hi": summary_s["c_ci_hi"],
        "c_p_analytic": pc,
        "c_prime": summary_s["c_prime"],
        "c_prime_ci_lo": summary_s["c_prime_ci_lo"],
        "c_prime_ci_hi": summary_s["c_prime_ci_hi"],
        "ab": summary_s["ab"],
        "ab_ci_lo": summary_s["ab_ci_lo"],
        "ab_ci_hi": summary_s["ab_ci_hi"],
        "ab_p": summary_s["ab_p"],
        "ab_p_floor": summary_s["ab_p_floor"],
        "ab_nboot": summary_s["ab_nboot"],
        "proportion_ab_over_c": (
            summary_s["ab"] / summary_s["c"] if summary_s["c"] not in (0, np.nan) else np.nan
        ),
        "vif_max": np.nan,
        "same_as_shared_n": np.nan,
    }
    if abs(rec["a"] - _ra) > 1e-8:
        raise AssertionError("mediator-only a diverged from Spearman")
    return rec


def run_partial(bundle: dict) -> list[dict]:
    cohort = bundle["cohort"]
    rows = []
    for focal, outcome, cov_names in PARTIAL_SPECS:
        names = [focal, outcome, *cov_names]
        df = align(bundle, list(dict.fromkeys(names)))
        n = int(len(df))
        y = df[outcome].to_numpy()
        x = df[focal].to_numpy()
        covs = [df[c].to_numpy() for c in cov_names]
        rho, p = partial_spearman(y, x, covs)
        rng = rng_for("partial", cohort, focal, outcome, ",".join(cov_names))

        def point_fn(idx: np.ndarray, y=y, x=x, covs=covs) -> dict[str, float]:
            r, _pp = partial_spearman(y[idx], x[idx], [c[idx] for c in covs])
            if not np.isfinite(r):
                raise ValueError("nonfinite partial")
            return {"rho": float(r)}

        _point, samples = bootstrap(point_fn, n, rng)
        arr = samples["rho"]
        if len(arr) < MIN_BOOT:
            lo, hi = np.nan, np.nan
        else:
            lo = float(np.percentile(arr, 2.5))
            hi = float(np.percentile(arr, 97.5))
        rows.append(
            {
                "cohort": cohort,
                "focal": focal,
                "outcome": outcome,
                "covariates": "+".join(cov_names),
                "n": n,
                "rho": rho,
                "p_analytic": p,
                "ci_lo": lo,
                "ci_hi": hi,
                "nboot": int(len(arr)),
            }
        )
    return rows


def run_assoc(bundle: dict) -> list[dict]:
    cohort = bundle["cohort"]
    rows = []
    for a, b in ASSOC_PAIRS:
        df = align(bundle, [a, b])
        sm = boot_spearman(
            df[a].to_numpy(),
            df[b].to_numpy(),
            rng_for("assoc", cohort, a, b),
        )
        rows.append(
            {
                "cohort": cohort,
                "x": a,
                "y": b,
                "n": int(len(df)),
                "rho": sm["rho"],
                "p_analytic": sm["p_analytic"],
                "ci_lo": sm["rho_ci_lo"],
                "ci_hi": sm["rho_ci_hi"],
                "nboot": int(sm["rho_nboot"]),
            }
        )
    return rows


def reproduction(assoc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, x, y, n_prior, rho_prior, p_prior in PRIOR:
        hit = assoc[(assoc["cohort"] == cohort) & (assoc["x"] == x) & (assoc["y"] == y)]
        if len(hit) != 1:
            raise AssertionError(f"missing assoc {cohort} {x} {y}")
        r = hit.iloc[0]
        rows.append(
            {
                "cohort": cohort,
                "x": x,
                "y": y,
                "n": int(r["n"]),
                "n_prior": n_prior,
                "rho": float(r["rho"]),
                "rho_prior": rho_prior,
                "abs_delta_rho": abs(float(r["rho"]) - rho_prior),
                "p_analytic": float(r["p_analytic"]),
                "p_prior": p_prior,
                "abs_delta_p": abs(float(r["p_analytic"]) - p_prior),
                "match": (
                    int(r["n"]) == n_prior
                    and abs(float(r["rho"]) - rho_prior) < 1e-8
                    and abs(float(r["p_analytic"]) - p_prior) < 1e-8
                ),
            }
        )
    return pd.DataFrame(rows)


def indirect_call(row: pd.Series) -> str:
    ab = ci_side(float(row["ab_ci_lo"]), float(row["ab_ci_hi"]))
    a = ci_side(float(row["a_ci_lo"]), float(row["a_ci_hi"]))
    b = ci_side(float(row["b_ci_lo"]), float(row["b_ci_hi"]))
    c = ci_side(float(row["c_ci_lo"]), float(row["c_ci_hi"]))
    if ab == "below" and a == "above" and b == "below":
        ab_txt = "indirect inverse (a CI > 0, b CI < 0, ab CI < 0)"
    elif ab == "below":
        ab_txt = "ab CI below 0; a and b are not both TROP2-up and immune-down"
    elif ab == "above":
        ab_txt = "ab CI above 0"
    elif ab == "crosses":
        ab_txt = "ab CI crosses 0"
    else:
        ab_txt = "ab not estimated"
    c_map = {
        "below": "total CI below 0",
        "above": "total CI above 0",
        "crosses": "total CI crosses 0",
        "na": "total not estimated",
    }
    return f"{ab_txt}; {c_map[c]}"


def apply_rule(primary: pd.DataFrame) -> pd.DataFrame:
    primary = primary.copy()
    primary["call"] = primary.apply(indirect_call, axis=1)
    primary["mechanism_rule_met"] = False
    for cohort in primary["cohort"].unique():
        for outcome in OUTCOMES:
            cldn = primary[
                (primary["cohort"] == cohort)
                & (primary["outcome"] == outcome)
                & (primary["mediator"] == "CLDN4")
                & (primary["model"] == "parallel")
                & (primary["sample_set"] == "shared")
            ]
            epc = primary[
                (primary["cohort"] == cohort)
                & (primary["outcome"] == outcome)
                & (primary["mediator"] == "EPCAM")
                & (primary["model"] == "parallel")
                & (primary["sample_set"] == "shared")
            ]
            if len(cldn) != 1 or len(epc) != 1:
                continue
            cr, er = cldn.iloc[0], epc.iloc[0]
            met = (
                ci_side(float(cr["a_ci_lo"]), float(cr["a_ci_hi"])) == "above"
                and ci_side(float(cr["b_ci_lo"]), float(cr["b_ci_hi"])) == "below"
                and ci_side(float(cr["ab_ci_lo"]), float(cr["ab_ci_hi"])) == "below"
                and ci_side(float(cr["c_ci_lo"]), float(cr["c_ci_hi"])) == "below"
                and ci_side(float(er["ab_ci_lo"]), float(er["ab_ci_hi"])) != "below"
            )
            primary.loc[cldn.index, "mechanism_rule_met"] = met
    return primary


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_r(r: float) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    # Keep a nonzero bound from rounding to 0.000 when it decides which side of 0 the CI is on.
    if 0 < abs(float(r)) < 5e-4:
        return f"{float(r):+.1e}"
    return f"{float(r):+.3f}"


def fmt_ci(lo: float, hi: float) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return f"{fmt_r(lo)} to {fmt_r(hi)}"


def md_table(df: pd.DataFrame, cols: list[tuple[str, str, callable]]) -> str:
    header = "| " + " | ".join(h for _k, h, _f in cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for _, row in df.iterrows():
        cells = []
        for key, _h, fn in cols:
            cells.append(str(fn(row[key])))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def forest(ax, labels, est, lo, hi, colors, xlabel: str) -> None:
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#666666", lw=0.8)
    for yi, e, a, b, col in zip(y, est, lo, hi, colors):
        if not (np.isfinite(e) and np.isfinite(a) and np.isfinite(b)):
            continue
        ax.plot([a, b], [yi, yi], color=col, lw=1.6, solid_capstyle="round")
        ax.plot([e], [yi], "o", color=col, ms=4.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel(xlabel, fontsize=8)
    style(ax)


MEDIATOR_COLOR = {"CLDN4": "#0072B2", "EPCAM": "#E69F00"}


def plot_missingness(cov: pd.DataFrame, path: Path) -> None:
    symbols = ["TACSTD2", "CLDN4", "EPCAM", "CD8A", "HLA-A", "HLA-B", "HLA-C", "B2M"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    cohorts = ["LUAD", "LSCC"]
    x = np.arange(len(symbols))
    width = 0.36
    for i, cohort in enumerate(cohorts):
        sub = cov[cov["cohort"] == cohort].set_index("symbol")
        vals = [float(sub.loc[s, "pct_missing"]) if s in sub.index else np.nan for s in symbols]
        ax.bar(x + (i - 0.5) * width, vals, width=width, label=cohort, color=["#4C78A8", "#F58518"][i])
    ax.set_xticks(x)
    ax.set_xticklabels(symbols, rotation=30, ha="right")
    ax.set_ylabel("% tumors missing")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Protein quantification, CPTAC tumor TMT")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_indirect(primary: pd.DataFrame, path: Path) -> None:
    df = primary[primary["sample_set"] == "shared"].copy()
    order = []
    for cohort in ("LSCC", "LUAD"):
        for outcome in ("CD8A", "MHC1"):
            for model in ("single", "parallel"):
                for mediator in ("CLDN4", "EPCAM"):
                    order.append((cohort, outcome, model, mediator))
    labels, est, lo, hi, colors = [], [], [], [], []
    for cohort, outcome, model, mediator in order:
        hit = df[
            (df["cohort"] == cohort)
            & (df["outcome"] == outcome)
            & (df["model"] == model)
            & (df["mediator"] == mediator)
        ].iloc[0]
        labels.append(f"{cohort} {outcome} {model} {mediator} (n={int(hit['n'])})")
        est.append(float(hit["ab"]))
        lo.append(float(hit["ab_ci_lo"]))
        hi.append(float(hit["ab_ci_hi"]))
        colors.append(MEDIATOR_COLOR[mediator])
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    forest(ax, labels, est, lo, hi, colors, "Indirect effect ab, z-rank OLS (95% bootstrap CI)")
    ax.set_title("TROP2 → mediator → CD8A / MHC-I, shared complete cases")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_partial(partial: pd.DataFrame, assoc: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), sharex=True)
    specs = [
        ("TACSTD2", "total", None),
        ("TACSTD2", "CLDN4", ("CLDN4",)),
        ("TACSTD2", "EPCAM", ("EPCAM",)),
        ("TACSTD2", "CLDN4+EPCAM", ("CLDN4", "EPCAM")),
        ("CLDN4", "total", None),
        ("CLDN4", "EPCAM", ("EPCAM",)),
        ("CLDN4", "TACSTD2", ("TACSTD2",)),
        ("EPCAM", "total", None),
        ("EPCAM", "CLDN4", ("CLDN4",)),
    ]
    for ax, cohort, outcome in (
        (axes[0, 0], "LUAD", "CD8A"),
        (axes[0, 1], "LSCC", "CD8A"),
        (axes[1, 0], "LUAD", "MHC1"),
        (axes[1, 1], "LSCC", "MHC1"),
    ):
        labels, est, lo, hi, colors = [], [], [], [], []
        for focal, tag, cov in specs:
            if cov is None:
                hit = assoc[(assoc["cohort"] == cohort) & (assoc["x"] == focal) & (assoc["y"] == outcome)]
                if len(hit) != 1:
                    continue
                r = hit.iloc[0]
                rho, a, b, n = float(r["rho"]), float(r["ci_lo"]), float(r["ci_hi"]), int(r["n"])
            else:
                key = "+".join(cov)
                hit = partial[
                    (partial["cohort"] == cohort)
                    & (partial["focal"] == focal)
                    & (partial["outcome"] == outcome)
                    & (partial["covariates"] == key)
                ]
                if len(hit) != 1:
                    continue
                r = hit.iloc[0]
                rho, a, b, n = float(r["rho"]), float(r["ci_lo"]), float(r["ci_hi"]), int(r["n"])
            labels.append(f"{focal} | {tag} (n={n})")
            est.append(rho)
            lo.append(a)
            hi.append(b)
            colors.append("#333333" if tag == "total" else MEDIATOR_COLOR.get(tag.split("+")[0], "#009E73"))
        forest(ax, labels, est, lo, hi, colors, "Spearman ρ (95% bootstrap CI)")
        ax.set_title(f"{cohort} {outcome}", fontsize=10)
    fig.suptitle("Total and partial protein associations", fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_scatter(bundles: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 6.0), sharex=False, sharey=False)
    genes = ["TACSTD2", "CLDN4", "EPCAM"]
    for row, cohort in enumerate(("LUAD", "LSCC")):
        b = bundles[cohort]
        cldn = b["series"]["CLDN4"]
        y = b["series"]["CD8A"]
        for col, gene in enumerate(genes):
            ax = axes[row, col]
            x = b["series"][gene]
            df = pd.concat([x.rename("x"), y.rename("y"), cldn.rename("cldn4")], axis=1)
            present = df.dropna(subset=["x", "y", "cldn4"])
            missing = df[df["cldn4"].isna()].dropna(subset=["x", "y"])
            ax.scatter(present["x"], present["y"], s=12, c="#0072B2", alpha=0.8, linewidths=0, label="CLDN4 quantified")
            if len(missing):
                ax.scatter(missing["x"], missing["y"], s=12, c="#B0B0B0", alpha=0.9, linewidths=0, label="CLDN4 missing")
            both = df.dropna(subset=["x", "y"])
            rho, p = spearman_stat(both["x"].to_numpy(), both["y"].to_numpy())
            ax.set_title(f"{cohort} {gene} vs CD8A\nρ={rho:+.3f}, n={len(both)}, p={fmt_p(p)}", fontsize=8)
            if row == 1:
                ax.set_xlabel(f"{gene} protein, log2", fontsize=8)
            if col == 0:
                ax.set_ylabel("CD8A protein, log2", fontsize=8)
            style(ax)
            if row == 0 and col == 2:
                ax.legend(frameon=False, fontsize=6, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def side_phrase(lo: float, hi: float) -> str:
    s = ci_side(lo, hi)
    return {"below": "entirely below 0", "above": "entirely above 0", "crosses": "crosses 0", "na": "not estimated"}[s]


def one(df: pd.DataFrame, **kw) -> pd.Series:
    m = df
    for k, v in kw.items():
        m = m[m[k] == v]
    if len(m) != 1:
        raise AssertionError(f"expected 1 row for {kw}, found {len(m)}")
    return m.iloc[0]


def rule_bits(cr: pd.Series, er: pd.Series) -> tuple[list[str], list[str]]:
    checks = [
        ("CLDN4 a CI entirely above 0", ci_side(float(cr["a_ci_lo"]), float(cr["a_ci_hi"])) == "above"),
        ("CLDN4 b CI entirely below 0", ci_side(float(cr["b_ci_lo"]), float(cr["b_ci_hi"])) == "below"),
        ("CLDN4 ab CI entirely below 0", ci_side(float(cr["ab_ci_lo"]), float(cr["ab_ci_hi"])) == "below"),
        ("total c CI entirely below 0", ci_side(float(cr["c_ci_lo"]), float(cr["c_ci_hi"])) == "below"),
        (
            "EPCAM ab CI not entirely below 0",
            ci_side(float(er["ab_ci_lo"]), float(er["ab_ci_hi"])) != "below",
        ),
    ]
    met = [name for name, ok in checks if ok]
    unmet = [name for name, ok in checks if not ok]
    return met, unmet


def coef_sentence(label: str, row: pd.Series, key: str) -> str:
    return (
        f"{label} {fmt_r(float(row[key]))} "
        f"(95% CI {fmt_ci(float(row[key + '_ci_lo']), float(row[key + '_ci_hi']))})"
    )


def assoc_phrase(r: pd.Series, label: str) -> str:
    return (
        f"{label} ρ={fmt_r(float(r['rho']))} "
        f"(n={int(r['n'])}, 95% CI {fmt_ci(float(r['ci_lo']), float(r['ci_hi']))}, "
        f"{side_phrase(float(r['ci_lo']), float(r['ci_hi']))})"
    )


def partial_phrase(r: pd.Series, label: str) -> str:
    return (
        f"{label} partial ρ={fmt_r(float(r['rho']))} "
        f"(n={int(r['n'])}, 95% CI {fmt_ci(float(r['ci_lo']), float(r['ci_hi']))}, "
        f"q={fmt_p(float(r['q_partial']))}, "
        f"{side_phrase(float(r['ci_lo']), float(r['ci_hi']))})"
    )


def pattern_lines(
    assoc: pd.DataFrame,
    partial: pd.DataFrame,
    primary: pd.DataFrame,
    purity: pd.DataFrame,
) -> list[str]:
    """Prose that only restates CI sides and formatted coefficients.

    The sentences below are templates. If a CI side is not the one the sentence
    names, stop instead of writing a mismatched paragraph.
    """
    def need(cond: bool, msg: str) -> None:
        if not cond:
            raise SystemExit(f"finding template does not match the tables: {msg}")

    for c in ("LUAD", "LSCC"):
        for y in OUTCOMES:
            r = one(assoc, cohort=c, x="TACSTD2", y=y)
            need(ci_side(float(r["ci_lo"]), float(r["ci_hi"])) == "crosses", f"all-tumor TROP2 {c} {y}")
            for gene in ("CLDN4", "EPCAM"):
                g = one(assoc, cohort="LSCC", x=gene, y=y)
                need(ci_side(float(g["ci_lo"]), float(g["ci_hi"])) == "below", f"LSCC {gene} {y}")
    for y in OUTCOMES:
        for focal, cov in (
            ("CLDN4", "EPCAM"),
            ("CLDN4", "TACSTD2"),
            ("EPCAM", "CLDN4"),
        ):
            r = one(partial, cohort="LSCC", focal=focal, outcome=y, covariates=cov)
            need(ci_side(float(r["ci_lo"]), float(r["ci_hi"])) == "below", f"LSCC partial {focal}|{cov} {y}")
    need(
        ci_side(*[float(one(assoc, cohort="LUAD", x="TACSTD2", y="CLDN4")[k]) for k in ("ci_lo", "ci_hi")]) == "above",
        "LUAD TROP2-CLDN4",
    )
    need(
        ci_side(*[float(one(assoc, cohort="LUAD", x="TACSTD2", y="EPCAM")[k]) for k in ("ci_lo", "ci_hi")]) == "above",
        "LUAD TROP2-EPCAM",
    )
    for y in ("CLDN4", "EPCAM"):
        r = one(assoc, cohort="LSCC", x="TACSTD2", y=y)
        need(ci_side(float(r["ci_lo"]), float(r["ci_hi"])) == "crosses", f"LSCC TROP2-{y}")
    shared = primary[primary["sample_set"] == "shared"]
    need(
        int(shared.apply(lambda r: ci_side(float(r["ab_ci_lo"]), float(r["ab_ci_hi"])) == "crosses", axis=1).sum())
        == len(shared),
        "not every shared ab CI crosses 0",
    )
    pur_cross = int(
        purity.apply(lambda r: ci_side(float(r["ab_ci_lo"]), float(r["ab_ci_hi"])) == "crosses", axis=1).sum()
    )
    need(pur_cross == len(purity), "not every purity ab CI crosses 0")

    lines: list[str] = ["## What the paths show", ""]
    trop = [
        assoc_phrase(one(assoc, cohort=c, x="TACSTD2", y=y), f"{c} TROP2 vs {y}")
        for c in ("LUAD", "LSCC")
        for y in OUTCOMES
    ]
    lines.append(
        "TROP2 versus CD8A and MHC-I, using every tumor in which both proteins were quantified: "
        + "; ".join(trop)
        + "."
    )
    lines.append("")
    lsc_bits = []
    for y, name in (("CD8A", "CD8A"), ("MHC1", "MHC-I")):
        lsc_bits.append(assoc_phrase(one(assoc, cohort="LSCC", x="CLDN4", y=y), f"CLDN4 vs {name}"))
        lsc_bits.append(assoc_phrase(one(assoc, cohort="LSCC", x="EPCAM", y=y), f"EPCAM vs {name}"))
    lines.append("In LSCC the inverse CD8A and MHC-I associations are with CLDN4 and with EPCAM. " + "; ".join(lsc_bits) + ".")
    lines.append("")
    hold = []
    for y, name in (("CD8A", "CD8A"), ("MHC1", "MHC-I")):
        hold.append(
            partial_phrase(
                one(partial, cohort="LSCC", focal="CLDN4", outcome=y, covariates="EPCAM"),
                f"CLDN4 vs {name} given EPCAM",
            )
        )
        hold.append(
            partial_phrase(
                one(partial, cohort="LSCC", focal="CLDN4", outcome=y, covariates="TACSTD2"),
                f"CLDN4 vs {name} given TROP2",
            )
        )
        hold.append(
            partial_phrase(
                one(partial, cohort="LSCC", focal="EPCAM", outcome=y, covariates="CLDN4"),
                f"EPCAM vs {name} given CLDN4",
            )
        )
    lines.append(
        "Each of those LSCC associations remains after the other epithelial protein, or TROP2, is held constant. "
        "Partial Spearman on the CLDN4-quantified tumors: "
        + "; ".join(hold)
        + "."
    )
    lines.append("")
    a_bits = [
        assoc_phrase(one(assoc, cohort="LUAD", x="TACSTD2", y="CLDN4"), "LUAD TROP2 vs CLDN4"),
        assoc_phrase(one(assoc, cohort="LUAD", x="TACSTD2", y="EPCAM"), "LUAD TROP2 vs EPCAM"),
        assoc_phrase(one(assoc, cohort="LSCC", x="TACSTD2", y="CLDN4"), "LSCC TROP2 vs CLDN4"),
        assoc_phrase(one(assoc, cohort="LSCC", x="TACSTD2", y="EPCAM"), "LSCC TROP2 vs EPCAM"),
        assoc_phrase(one(assoc, cohort="LUAD", x="CLDN4", y="EPCAM"), "LUAD CLDN4 vs EPCAM"),
        assoc_phrase(one(assoc, cohort="LSCC", x="CLDN4", y="EPCAM"), "LSCC CLDN4 vs EPCAM"),
    ]
    vif_luad = float(one(primary, cohort="LUAD", outcome="CD8A", mediator="CLDN4", model="parallel")["vif_max"])
    vif_lscc = float(one(primary, cohort="LSCC", outcome="CD8A", mediator="CLDN4", model="parallel")["vif_max"])
    lines.append(
        "The a path, TROP2 to the epithelial protein, is histology-specific. "
        + "; ".join(a_bits)
        + f". Max VIF of the three rank predictors is {vif_luad:.2f} in LUAD and {vif_lscc:.2f} in LSCC."
    )
    lines.append("")
    shared = primary[primary["sample_set"] == "shared"]
    ab_cross = int(
        shared.apply(lambda r: ci_side(float(r["ab_ci_lo"]), float(r["ab_ci_hi"])) == "crosses", axis=1).sum()
    )
    n_ab = int(len(shared))

    def _side_count(frame: pd.DataFrame, key: str, side: str) -> tuple[int, int]:
        m = frame.apply(lambda r: ci_side(float(r[f"{key}_ci_lo"]), float(r[f"{key}_ci_hi"])) == side, axis=1)
        return int(m.sum()), int(len(frame))

    lscc = shared[shared["cohort"] == "LSCC"]
    luad_cl = shared[(shared["cohort"] == "LUAD") & (shared["mediator"] == "CLDN4")]
    lscc_a_below, lscc_a_n = _side_count(lscc[lscc["mediator"] == "CLDN4"], "a", "crosses")
    lscc_b_below, lscc_b_n = _side_count(lscc, "b", "below")
    luad_a_above, luad_a_n = _side_count(luad_cl, "a", "above")
    luad_b_cross, luad_b_n = _side_count(luad_cl, "b", "crosses")
    lines.append(
        f"Shared-sample indirect effects: {ab_cross} of {n_ab} ab confidence intervals cross 0. "
        f"LSCC b intervals entirely below 0: {lscc_b_below} of {lscc_b_n}. "
        f"LSCC TROP2→CLDN4 a intervals that cross 0: {lscc_a_below} of {lscc_a_n}. "
        f"LUAD CLDN4-mediator a intervals entirely above 0: {luad_a_above} of {luad_a_n}. "
        f"LUAD CLDN4-mediator b intervals that cross 0: {luad_b_cross} of {luad_b_n}. "
        f"Primary q on the 16 indirect tests has minimum {fmt_p(float(primary['q_primary'].min()))}."
    )
    lines.append("")
    luad_m = one(primary, cohort="LUAD", outcome="MHC1", mediator="CLDN4", model="single")
    luad_part = one(partial, cohort="LUAD", focal="TACSTD2", outcome="MHC1", covariates="CLDN4")
    luad_both = one(partial, cohort="LUAD", focal="TACSTD2", outcome="MHC1", covariates="CLDN4+EPCAM")
    pur_c = one(purity, cohort="LUAD", outcome="MHC1", mediator="CLDN4", model="parallel")
    lines.append(
        "The one shared-sample total effect whose interval lies entirely below 0 is LUAD TROP2 vs MHC-I "
        f"on the CLDN4-quantified tumors: c={fmt_r(float(luad_m['c']))} "
        f"(n={int(luad_m['n'])}, 95% CI {fmt_ci(float(luad_m['c_ci_lo']), float(luad_m['c_ci_hi']))}). "
        "The all-tumor TROP2–MHC-I interval crosses 0, so this is the CLDN4-complete subset, not the "
        f"full LUAD series. {partial_phrase(luad_part, 'TROP2 vs MHC-I given CLDN4')}. "
        f"{partial_phrase(luad_both, 'The same contrast given CLDN4 and EPCAM')}. "
        f"The CLDN4 indirect effect on that subset is ab={fmt_r(float(luad_m['ab']))} "
        f"(95% CI {fmt_ci(float(luad_m['ab_ci_lo']), float(luad_m['ab_ci_hi']))}), "
        f"ab/c={float(luad_m['proportion_ab_over_c']):+.3f}. "
        "The product is positive because the CLDN4 b coefficient is positive. "
        "The inverse total association on this subset is the direct path. "
        "WES-purity rank adjustment of the parallel model on this subset "
        f"(n={int(pur_c['n'])}) moves the total-effect interval to "
        f"{fmt_ci(float(pur_c['c_ci_lo']), float(pur_c['c_ci_hi']))} "
        f"({side_phrase(float(pur_c['c_ci_lo']), float(pur_c['c_ci_hi']))}) and the direct-effect interval to "
        f"{fmt_ci(float(pur_c['c_prime_ci_lo']), float(pur_c['c_prime_ci_hi']))} "
        f"({side_phrase(float(pur_c['c_prime_ci_lo']), float(pur_c['c_prime_ci_hi']))})."
    )
    lines.append("")
    # Purity b sides for LSCC, stated from the table rather than assumed.
    b_notes = []
    for outcome, name in (("CD8A", "CD8A"), ("MHC1", "MHC-I")):
        for mediator in ("CLDN4", "EPCAM"):
            r = one(purity, cohort="LSCC", outcome=outcome, mediator=mediator, model="parallel")
            b_notes.append(
                f"{name} {mediator} b={fmt_r(float(r['b']))} "
                f"(CI {fmt_ci(float(r['b_ci_lo']), float(r['b_ci_hi']))}, "
                f"{side_phrase(float(r['b_ci_lo']), float(r['b_ci_hi']))})"
            )
    pur_ab_cross = int(
        purity.apply(lambda r: ci_side(float(r["ab_ci_lo"]), float(r["ab_ci_hi"])) == "crosses", axis=1).sum()
    )
    lines.append(
        f"WES-purity rank adjustment: {pur_ab_cross} of {len(purity)} indirect-effect intervals cross 0. "
        "LSCC parallel-model b coefficients under that adjustment: "
        + "; ".join(b_notes)
        + ". Unadjusted partial Spearman had already kept both LSCC epithelial proteins inverse; "
        "the purity-adjusted regression coefficient is a different estimand and is where EPCAM and CLDN4 separate."
    )
    lines.append("")
    n_met = int(
        primary.loc[primary["model"].eq("parallel") & primary["mediator"].eq("CLDN4"), "mechanism_rule_met"].sum()
    )
    lines.append(
        f"The pre-specified mediation rule, which requires an inverse total TROP2 association and an "
        f"inverse CLDN4 indirect effect that EPCAM does not copy, is met for {n_met} of 4 cohort–outcome pairs."
    )
    lines.append("")
    return lines


def write_finding(
    path: Path,
    cov: pd.DataFrame,
    assoc: pd.DataFrame,
    partial: pd.DataFrame,
    primary: pd.DataFrame,
    purity: pd.DataFrame,
    mediator_only: pd.DataFrame,
    repro: pd.DataFrame,
    cohort_meta: dict,
) -> None:
    lines: list[str] = []
    lines.append("# Finding — CPTAC protein: TROP2, partial CLDN4, CD8A/MHC-I, EPCAM control")
    lines.append("")
    lines.append(
        "Public CPTAC TMT freeze v1.2 tumor protein. **LUAD** (Gillette *Cell* 2020, n=110) and "
        "**LSCC** (Satpathy *Cell* 2021, n=108) kept separate. TROP2 is the protein product of "
        "TACSTD2. No log2 abundance is imputed. Coefficients below are written by `analyze.py` "
        "from `tables/`; displayed values are rounded to 3 decimals."
    )
    lines.append("")
    # Headline from CI sides, all-tumor TROP2 totals and parallel indirect effects.
    trop_totals = []
    for cohort in ("LUAD", "LSCC"):
        for outcome in OUTCOMES:
            r = one(assoc, cohort=cohort, x="TACSTD2", y=outcome)
            trop_totals.append(ci_side(float(r["ci_lo"]), float(r["ci_hi"])))
    par = primary[(primary["model"] == "parallel") & (primary["sample_set"] == "shared")]
    cldn_ab = [
        ci_side(float(r.ab_ci_lo), float(r.ab_ci_hi))
        for r in par[par["mediator"] == "CLDN4"].itertuples()
    ]
    ep_ab = [
        ci_side(float(r.ab_ci_lo), float(r.ab_ci_hi))
        for r in par[par["mediator"] == "EPCAM"].itertuples()
    ]
    lines.append(
        f"All-tumor TROP2 versus CD8A and versus MHC-I: "
        f"{sum(s == 'below' for s in trop_totals)} of {len(trop_totals)} bootstrap CIs lie entirely below 0, "
        f"{sum(s == 'above' for s in trop_totals)} entirely above 0, "
        f"{sum(s == 'crosses' for s in trop_totals)} cross 0. "
        f"Parallel-model indirect effects on the shared complete cases: "
        f"CLDN4 mediator, {sum(s == 'below' for s in cldn_ab)} of {len(cldn_ab)} entirely below 0; "
        f"EPCAM mediator, {sum(s == 'below' for s in ep_ab)} of {len(ep_ab)} entirely below 0."
    )
    lines.append("")
    lines.extend(pattern_lines(assoc, partial, primary, purity))
    lines.append("## Treatment-naive")
    lines.append("")
    for cohort, meta in cohort_meta.items():
        lines.append(
            f"{cohort} phenotype has {meta['n_pheno_cols']} columns. "
            f"Therapy-like column names found: {meta['treat_cols'] or 'none'}."
        )
    lines.append("")
    lines.append(
        "Treatment-naive is the published cohort definition: prospectively collected, previously "
        "untreated surgical resections. There is no ICI arm, no on-treatment biopsy, and no response "
        "label in these files. A path coefficient here is not acquired resistance and not a drug effect."
    )
    lines.append("")
    lines.append("## Estimands")
    lines.append("")
    lines.append(
        "Exposure X = TACSTD2 protein. Mediator under test M1 = CLDN4 protein. Control protein "
        "M2 = EPCAM, a second epithelial surface protein run through the same equations. Outcomes "
        "are CD8A protein and the MHC-I score used on the earlier page for this freeze: the mean of "
        "per-gene z-scores (ddof=0) of HLA-A, HLA-B, and HLA-C, requiring all three. B2M is not in "
        "that score. Ranks use average ties. Single-mediator z-rank OLS is M ~ X and Y ~ X + M. "
        "The product ab equals c − c′ on that scale, and a and c equal Spearman ρ. The parallel "
        "model is M1 ~ X, M2 ~ X, and Y ~ X + M1 + M2, so each indirect effect is adjusted for the "
        "other epithelial protein at the b step. Bootstrap: 2,000 resamples, seed `20260921`, ranks "
        "recomputed inside each resample, percentile 95% CI. Bootstrap p for ab is two-sided on the "
        "resampled products; a value of 0 is stored as 1/2000 and flagged. Primary q is "
        "Benjamini–Hochberg across the 16 shared-sample indirect tests (2 cohorts × 2 outcomes × "
        "2 mediators × single and parallel). Partial Spearman is the Pearson correlation of rank "
        "residuals and has its own BH family (the covariate list in `PARTIAL_SPECS`). WES purity "
        "is a rank covariate in a sensitivity refit, not part of the primary rule."
    )
    lines.append("")
    lines.append(
        "The mediation rule, applied to the parallel model on tumors with TACSTD2, CLDN4, EPCAM, "
        "and the outcome all quantified, requires all of: CLDN4 a CI entirely above 0, CLDN4 b CI "
        "entirely below 0, CLDN4 ab CI entirely below 0, total c CI entirely below 0, and EPCAM ab "
        "CI not entirely below 0. Where a total-effect CI includes 0, ab/c is left in "
        "`tables/mediation_shared.tsv` and is not read as a fraction mediated."
    )
    lines.append("")
    lines.append("## Reproduction of the eight earlier totals")
    lines.append("")
    lines.append(
        "Pairwise-complete Spearman ρ and analytic p for TACSTD2 and CLDN4 versus CD8A and MHC-I "
        "were compared with `summary.json` from the earlier protein page, which used these same "
        "freeze files. All eight match within 1e-8 on n, ρ, and p. That page's q and its "
        "WES-purity partial correlations are not re-adjudicated here."
    )
    lines.append("")
    show = repro.copy()
    lines.append(
        md_table(
            show,
            [
                ("cohort", "Cohort", str),
                ("x", "X", str),
                ("y", "Y", str),
                ("n", "n", lambda v: str(int(v))),
                ("rho", "ρ", fmt_r),
                ("abs_delta_rho", "|Δρ|", lambda v: f"{float(v):.1e}"),
                ("match", "Match", lambda v: "yes" if bool(v) else "no"),
            ],
        )
    )
    lines.append("")
    lines.append("## Partial quantification")
    lines.append("")
    lines.append(
        "CLDN4 is the incompletely quantified protein. Mediation uses complete cases only. "
        "Tumors missing CLDN4 stay in the all-tumor TROP2 and EPCAM scatters and are dropped "
        "from every path that includes CLDN4. The B2M row is absent from both matrices, which "
        "is why the MHC-I score stays HLA-A/B/C."
    )
    lines.append("")
    cov_show = cov[cov["symbol"].isin(["TACSTD2", "CLDN4", "EPCAM", "CD8A", "HLA-A", "HLA-B", "HLA-C", "B2M"])]
    lines.append(
        md_table(
            cov_show,
            [
                ("cohort", "Cohort", str),
                ("symbol", "Protein", str),
                ("row", "Matrix row", str),
                ("n_quantified", "Quantified", lambda v: str(int(v))),
                ("n_missing", "Missing", lambda v: str(int(v))),
                ("pct_missing", "% missing", lambda v: f"{float(v):.1f}"),
            ],
        )
    )
    lines.append("")
    lines.append("Shared complete-case n (TACSTD2, CLDN4, EPCAM, and the outcome):")
    lines.append("")
    for cohort in ("LUAD", "LSCC"):
        for outcome in OUTCOMES:
            r = one(primary, cohort=cohort, outcome=outcome, mediator="CLDN4", model="parallel")
            lines.append(f"- {cohort} {outcome}: n={int(r['n'])}")
    lines.append("")
    lines.append("## All-tumor Spearman")
    lines.append("")
    lines.append(
        "These are pairwise-complete. They are not the mediation-sample total effects. "
        "Analytic Spearman p is the usual t approximation. The CI is this page's bootstrap. "
        "No FDR is applied to this descriptive table."
    )
    lines.append("")
    lines.append(
        md_table(
            assoc,
            [
                ("cohort", "Cohort", str),
                ("x", "X", str),
                ("y", "Y", str),
                ("n", "n", lambda v: str(int(v))),
                ("rho", "ρ", fmt_r),
                ("ci_lo", "CI low", fmt_r),
                ("ci_hi", "CI high", fmt_r),
                ("p_analytic", "p", fmt_p),
            ],
        )
    )
    lines.append("")
    lines.append("## Mediation on the shared sample")
    lines.append("")
    lines.append(
        "a is the standardized rank coefficient, equal to Spearman ρ for the single-mediator "
        "and parallel a paths (those a paths are M ~ X, not adjusted for the other mediator). "
        "b and c′ are standardized rank regression coefficients. ab is their product. "
        "q is the 16-test primary BH on ab."
    )
    lines.append("")
    prim_show = primary[primary["sample_set"] == "shared"].copy()
    lines.append(
        md_table(
            prim_show,
            [
                ("cohort", "Cohort", str),
                ("outcome", "Y", str),
                ("mediator", "Mediator", str),
                ("model", "Model", str),
                ("n", "n", lambda v: str(int(v))),
                ("a", "a", fmt_r),
                ("b", "b", fmt_r),
                ("c", "c", fmt_r),
                ("c_prime", "c′", fmt_r),
                ("ab", "ab", fmt_r),
                ("ab_ci_lo", "ab CI low", fmt_r),
                ("ab_ci_hi", "ab CI high", fmt_r),
                ("ab_p", "ab p", fmt_p),
                ("q_primary", "q", fmt_p),
            ],
        )
    )
    lines.append("")
    lines.append("Parallel-model detail and the mediation rule:")
    lines.append("")
    for cohort in ("LSCC", "LUAD"):
        for outcome in OUTCOMES:
            cr = one(primary, cohort=cohort, outcome=outcome, mediator="CLDN4", model="parallel")
            er = one(primary, cohort=cohort, outcome=outcome, mediator="EPCAM", model="parallel")
            sr = one(primary, cohort=cohort, outcome=outcome, mediator="CLDN4", model="single")
            se = one(primary, cohort=cohort, outcome=outcome, mediator="EPCAM", model="single")
            met, unmet = rule_bits(cr, er)
            lines.append(
                f"**{cohort} {outcome}, n={int(cr['n'])}.** "
                f"Total c {fmt_r(float(cr['c']))} (CI {fmt_ci(float(cr['c_ci_lo']), float(cr['c_ci_hi']))}, "
                f"analytic p={fmt_p(float(cr['c_p_analytic']))}), which {side_phrase(float(cr['c_ci_lo']), float(cr['c_ci_hi']))}. "
                f"CLDN4 a {fmt_r(float(cr['a']))} (CI {fmt_ci(float(cr['a_ci_lo']), float(cr['a_ci_hi']))}). "
                f"Single-mediator CLDN4 b {fmt_r(float(sr['b']))} "
                f"(CI {fmt_ci(float(sr['b_ci_lo']), float(sr['b_ci_hi']))}), "
                f"c′ {fmt_r(float(sr['c_prime']))} "
                f"(CI {fmt_ci(float(sr['c_prime_ci_lo']), float(sr['c_prime_ci_hi']))}), "
                f"ab {fmt_r(float(sr['ab']))} "
                f"(CI {fmt_ci(float(sr['ab_ci_lo']), float(sr['ab_ci_hi']))}, "
                f"p={fmt_p(float(sr['ab_p']))}, q={fmt_p(float(sr['q_primary']))}). "
                f"Single-mediator EPCAM ab {fmt_r(float(se['ab']))} "
                f"(CI {fmt_ci(float(se['ab_ci_lo']), float(se['ab_ci_hi']))}, "
                f"q={fmt_p(float(se['q_primary']))}). "
                f"Parallel CLDN4 b {fmt_r(float(cr['b']))} "
                f"(CI {fmt_ci(float(cr['b_ci_lo']), float(cr['b_ci_hi']))}), "
                f"ab {fmt_r(float(cr['ab']))} "
                f"(CI {fmt_ci(float(cr['ab_ci_lo']), float(cr['ab_ci_hi']))}, "
                f"p={fmt_p(float(cr['ab_p']))}, q={fmt_p(float(cr['q_primary']))}). "
                f"Parallel EPCAM b {fmt_r(float(er['b']))} "
                f"(CI {fmt_ci(float(er['b_ci_lo']), float(er['b_ci_hi']))}), "
                f"ab {fmt_r(float(er['ab']))} "
                f"(CI {fmt_ci(float(er['ab_ci_lo']), float(er['ab_ci_hi']))}, "
                f"q={fmt_p(float(er['q_primary']))}). "
                f"Parallel direct c′ {fmt_r(float(cr['c_prime']))} "
                f"(CI {fmt_ci(float(cr['c_prime_ci_lo']), float(cr['c_prime_ci_hi']))}). "
                f"Max VIF of TACSTD2, CLDN4, and EPCAM ranks in the parallel outcome model: "
                f"{float(cr['vif_max']):.2f}. "
                f"Rule clauses met: {', '.join(met) if met else 'none'}. "
                f"Rule clauses not met: {', '.join(unmet) if unmet else 'none'}."
            )
            lines.append("")
    lines.append("## Partial Spearman")
    lines.append("")
    lines.append(
        "Partial ρ is the correlation of rank residuals, with analytic p from that Pearson "
        "correlation and a bootstrap CI. q is BH inside this covariate family only. "
        "A partial CI that excludes 0 is not promoted into the mediation rule."
    )
    lines.append("")
    lines.append(
        md_table(
            partial,
            [
                ("cohort", "Cohort", str),
                ("focal", "Focal", str),
                ("outcome", "Y", str),
                ("covariates", "Given", str),
                ("n", "n", lambda v: str(int(v))),
                ("rho", "partial ρ", fmt_r),
                ("ci_lo", "CI low", fmt_r),
                ("ci_hi", "CI high", fmt_r),
                ("p_analytic", "p", fmt_p),
                ("q_partial", "q", fmt_p),
            ],
        )
    )
    lines.append("")
    lines.append("Contrasts that carry the partial-CLDN4 question:")
    lines.append("")
    for cohort in ("LSCC", "LUAD"):
        for outcome in OUTCOMES:
            tot = one(assoc, cohort=cohort, x="TACSTD2", y=outcome)
            p_cl = one(partial, cohort=cohort, focal="TACSTD2", outcome=outcome, covariates="CLDN4")
            p_ep = one(partial, cohort=cohort, focal="TACSTD2", outcome=outcome, covariates="EPCAM")
            p_both = one(partial, cohort=cohort, focal="TACSTD2", outcome=outcome, covariates="CLDN4+EPCAM")
            c_ep = one(partial, cohort=cohort, focal="CLDN4", outcome=outcome, covariates="EPCAM")
            c_tx = one(partial, cohort=cohort, focal="CLDN4", outcome=outcome, covariates="TACSTD2")
            lines.append(
                f"**{cohort} {outcome}.** All-tumor TROP2 ρ={fmt_r(float(tot['rho']))} "
                f"(n={int(tot['n'])}, CI {fmt_ci(float(tot['ci_lo']), float(tot['ci_hi']))}). "
                f"TROP2 given CLDN4: partial ρ={fmt_r(float(p_cl['rho']))} "
                f"(n={int(p_cl['n'])}, CI {fmt_ci(float(p_cl['ci_lo']), float(p_cl['ci_hi']))}, "
                f"q={fmt_p(float(p_cl['q_partial']))}). "
                f"TROP2 given EPCAM: partial ρ={fmt_r(float(p_ep['rho']))} "
                f"(n={int(p_ep['n'])}, CI {fmt_ci(float(p_ep['ci_lo']), float(p_ep['ci_hi']))}, "
                f"q={fmt_p(float(p_ep['q_partial']))}). "
                f"TROP2 given CLDN4 and EPCAM: partial ρ={fmt_r(float(p_both['rho']))} "
                f"(n={int(p_both['n'])}, CI {fmt_ci(float(p_both['ci_lo']), float(p_both['ci_hi']))}, "
                f"q={fmt_p(float(p_both['q_partial']))}). "
                f"CLDN4 given EPCAM: partial ρ={fmt_r(float(c_ep['rho']))} "
                f"(n={int(c_ep['n'])}, CI {fmt_ci(float(c_ep['ci_lo']), float(c_ep['ci_hi']))}, "
                f"q={fmt_p(float(c_ep['q_partial']))}). "
                f"CLDN4 given TROP2: partial ρ={fmt_r(float(c_tx['rho']))} "
                f"(n={int(c_tx['n'])}, CI {fmt_ci(float(c_tx['ci_lo']), float(c_tx['ci_hi']))}, "
                f"q={fmt_p(float(c_tx['q_partial']))})."
            )
            lines.append("")
    lines.append("## Sensitivity")
    lines.append("")
    lines.append(
        "Mediator-complete fits drop the requirement that the other epithelial protein was "
        "quantified. WES-purity fits add purity rank to every equation on tumors that also "
        "have a purity value. q in the purity file is BH across those purity indirect tests "
        "only. Neither family replaces the shared-sample primary rule."
    )
    lines.append("")
    lines.append("Mediator-complete single-mediator models:")
    lines.append("")
    lines.append(
        md_table(
            mediator_only,
            [
                ("cohort", "Cohort", str),
                ("outcome", "Y", str),
                ("mediator", "Mediator", str),
                ("n", "n", lambda v: str(int(v))),
                ("a", "a", fmt_r),
                ("b", "b", fmt_r),
                ("c", "c", fmt_r),
                ("ab", "ab", fmt_r),
                ("ab_ci_lo", "ab CI low", fmt_r),
                ("ab_ci_hi", "ab CI high", fmt_r),
                ("ab_p", "ab p", fmt_p),
            ],
        )
    )
    lines.append("")
    lines.append("WES-purity rank-adjusted paths:")
    lines.append("")
    lines.append(
        md_table(
            purity,
            [
                ("cohort", "Cohort", str),
                ("outcome", "Y", str),
                ("mediator", "Mediator", str),
                ("model", "Model", str),
                ("n", "n", lambda v: str(int(v))),
                ("a", "a", fmt_r),
                ("b", "b", fmt_r),
                ("c", "c", fmt_r),
                ("ab", "ab", fmt_r),
                ("ab_ci_lo", "ab CI low", fmt_r),
                ("ab_ci_hi", "ab CI high", fmt_r),
                ("ab_p", "ab p", fmt_p),
                ("q_purity", "q", fmt_p),
            ],
        )
    )
    lines.append("")
    lines.append("## What this page does not do")
    lines.append("")
    lines.append(
        "It does not impute CLDN4. It does not pool LUAD with LSCC. It does not substitute RNA "
        "for a missing protein. It does not use ImmuneScore, GEP, or CIBERSORT as the outcome; "
        "those are already on other pages. It does not restate private PDX numbers. It does not "
        "turn a cross-sectional product of rank coefficients into a causal demonstration that "
        "TROP2 changes CD8 or MHC through CLDN4."
    )
    lines.append("")
    lines.append("Figures: `figures/fig1_missingness.png`, `figures/fig2_indirect_forest.png`, `figures/fig3_partial_forest.png`, `figures/fig4_cd8a_scatter.png`.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    self_check()
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data/cptac_trop2_cldn4_mediation")
    p.add_argument("--outdir", default=str(here))
    args = p.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    fig = out / "figures"
    tab = out / "tables"
    fig.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)

    bundles = {c: load_cohort(data, c) for c in ("LUAD", "LSCC")}
    cov = pd.DataFrame([row for b in bundles.values() for row in b["coverage"]])
    assoc = pd.DataFrame([row for b in bundles.values() for row in run_assoc(b)])
    repro = reproduction(assoc)
    if not bool(repro["match"].all()):
        repro.to_csv(tab / "reproduction_check.tsv", sep="\t", index=False)
        print(repro.to_string(index=False), file=sys.stderr)
        raise SystemExit("reproduction check failed; MHC score or gene rows do not match the earlier page")

    partial = pd.DataFrame([row for b in bundles.values() for row in run_partial(b)])
    partial["q_partial"] = bh_q(partial["p_analytic"].tolist())

    primary_rows = []
    shared_meta = {}
    purity_rows = []
    med_only_rows = []
    for cohort, bundle in bundles.items():
        for outcome in OUTCOMES:
            rows, meta = run_shared(bundle, outcome)
            primary_rows.extend(rows)
            shared_meta[(cohort, outcome)] = meta
            purity_rows.extend(run_purity(bundle, outcome))
            for mediator in ("CLDN4", "EPCAM"):
                med_only_rows.append(run_mediator_only(bundle, outcome, mediator))
    primary = apply_rule(pd.DataFrame(primary_rows))
    primary["q_primary"] = bh_q(primary["ab_p"].tolist())
    # apply_rule copied the frame before q; q was added after. mechanism flag does not depend on q.
    purity = pd.DataFrame(purity_rows)
    purity["q_purity"] = bh_q(purity["ab_p"].tolist())
    mediator_only = pd.DataFrame(med_only_rows)

    cov.to_csv(tab / "coverage.tsv", sep="\t", index=False)
    assoc.to_csv(tab / "associations.tsv", sep="\t", index=False)
    repro.to_csv(tab / "reproduction_check.tsv", sep="\t", index=False)
    partial.to_csv(tab / "partial_spearman.tsv", sep="\t", index=False)
    primary.to_csv(tab / "mediation_shared.tsv", sep="\t", index=False)
    purity.to_csv(tab / "mediation_purity.tsv", sep="\t", index=False)
    mediator_only.to_csv(tab / "mediation_mediator_complete.tsv", sep="\t", index=False)

    plot_missingness(cov, fig / "fig1_missingness.png")
    plot_indirect(primary, fig / "fig2_indirect_forest.png")
    plot_partial(partial, assoc, fig / "fig3_partial_forest.png")
    plot_scatter(bundles, fig / "fig4_cd8a_scatter.png")

    cohort_meta = {
        c: {"n_pheno_cols": bundles[c]["n_pheno_cols"], "treat_cols": bundles[c]["treat_cols"]}
        for c in bundles
    }
    write_finding(
        out / "FINDING.md",
        cov,
        assoc,
        partial,
        primary,
        purity,
        mediator_only,
        repro,
        cohort_meta,
    )
    summary = {
        "seed": SEED,
        "n_boot": N_BOOT,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": stats.__version__ if hasattr(stats, "__version__") else "see scipy",
        "mechanism_rule_met_n": int(
            primary.loc[primary["model"].eq("parallel") & primary["mediator"].eq("CLDN4"), "mechanism_rule_met"].sum()
        ),
        "shared_n": {f"{c}_{o}": shared_meta[(c, o)]["n"] for c in bundles for o in OUTCOMES},
        "treat_cols": {c: bundles[c]["treat_cols"] for c in bundles},
        "reproduction_all_match": bool(repro["match"].all()),
    }
    # scipy version
    import scipy
    summary["scipy"] = scipy.__version__
    (tab / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(primary[["cohort", "outcome", "mediator", "model", "n", "a", "b", "c", "ab", "ab_ci_lo", "ab_ci_hi", "q_primary", "mechanism_rule_met"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
