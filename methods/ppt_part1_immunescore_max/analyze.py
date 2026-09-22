#!/usr/bin/env python3
"""PPT Part 1. Maximize |ρ| of the locked CLDN4-high signature and of TACSTD2
versus ImmuneScore.

ImmuneScore is the PR 590 definition: the within-cohort mean of z-scores of
CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, TBX21. ESTIMATE immune ssGSEA
is a labeled sensitivity and is not renamed ImmuneScore.

The 221 genes and their order are the PR 590 signature. Prefix length is the
only size knob. No gene is added or dropped for its correlation with
ImmuneScore. TACSTD2 and CLDN4 are not in the list. An equal-weight score is
the mean of the within-cohort z-score of a signature score and the
within-cohort z-score of TACSTD2. Those weights are not fit to ImmuneScore.

Cohorts are OncoSG, GSE273377 discovery, GSE273377 validation, GSE282774, and
GSE233774 tumors. A spec is eligible only when all five strata have a finite
Spearman, so a cohort cannot be dropped to raise |ρ|. GSE10072, GSE11969, and
GSE248378 are not opened.

Covariate modes, all pre-specified: unadjusted; partial Spearman on published
PURITY (OncoSG) or ESTIMATE StromalScore (GEO); partial on KRT8+KRT18+KRT19;
partial on that purity covariate plus the three keratins; gene-level residual
on the purity covariate, then Spearman with the raw ImmuneScore.

The selected spec is chosen on these cohorts. Its meta p-value is the p-value
of a maximized |ρ|. The pre-specified baseline is the 221-gene z-mean,
unadjusted. A GEO-only selection is applied unchanged to OncoSG.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
METHODS = HERE.parent
sys.path.insert(0, str(METHODS / "cldn4_sig_max_rho"))

import max_effect as M  # noqa: E402
import analyze as A  # noqa: E402

TABLES = HERE / "tables"
FIGURES = HERE / "figures"
STUDIES = M.STUDIES
REQUIRED = M.REQUIRED
K_COV = {
    "none": 0,
    "purity": 1,
    "keratin": 3,
    "purity_keratin": 4,
    "residual": 0,
}
IMMUNE8_EXPECTED = {
    "OncoSG": -0.600,
    "GSE273377 discovery": -0.557,
    "GSE273377 validation": -0.482,
    "GSE282774": -0.604,
    "GSE233774 tumor": -0.601,
}


def self_test() -> None:
    rng = np.random.default_rng(20260922)
    score = rng.normal(size=(6, 48))
    y = rng.normal(size=48)
    rho = spearman_rows(score, y)
    for i in range(score.shape[0]):
        ref = A.spearman_pair(score[i], y)["rho"]
        if abs(rho[i] - ref) > 1e-8:
            raise SystemExit(f"spearman mismatch {rho[i]} vs {ref}")
    covs = [rng.normal(size=48), rng.normal(size=48)]
    rho_p = partial_rows(score, y, covs)
    for i in range(score.shape[0]):
        ref = A.partial_spearman(score[i], y, covs)["rho"]
        if abs(rho_p[i] - ref) > 1e-8:
            raise SystemExit(f"partial mismatch {rho_p[i]} vs {ref}")
    M.self_test()
    print("correlation self-test ok", flush=True)


def spearman_rows(score: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson of average ranks. Matches analyze.spearman_pair on finite rows."""
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    if score.ndim == 1:
        score = score.reshape(1, -1)
    out = np.full(score.shape[0], np.nan)
    finite_y = np.isfinite(y)
    if int(finite_y.sum()) < 6:
        return out
    y = y[finite_y]
    score = score[:, finite_y]
    ok = np.isfinite(score).all(axis=1) & (np.std(score, axis=1) > 0)
    if not np.any(ok) or np.std(y) == 0:
        return out
    rs = stats.rankdata(score[ok], axis=1)
    ry = stats.rankdata(y)
    rs = rs - rs.mean(axis=1, keepdims=True)
    ry = ry - ry.mean()
    denom = np.sqrt((rs ** 2).sum(axis=1) * (ry ** 2).sum())
    out[ok] = (rs @ ry) / denom
    return out


def partial_rows(score: np.ndarray, y: np.ndarray, covs: list[np.ndarray]) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    if score.ndim == 1:
        score = score.reshape(1, -1)
    covs = [np.asarray(c, dtype=float) for c in covs]
    mask = np.isfinite(y)
    for cov in covs:
        mask &= np.isfinite(cov)
    out = np.full(score.shape[0], np.nan)
    n = int(mask.sum())
    k = len(covs)
    if n < 8 + k:
        return out
    y = y[mask]
    cov_m = [cov[mask] for cov in covs]
    score = score[:, mask]
    ok = np.isfinite(score).all(axis=1) & (np.std(score, axis=1) > 0)
    if not np.any(ok):
        return out
    rs = stats.rankdata(score[ok], axis=1)
    ry = stats.rankdata(y)
    rc = np.column_stack([stats.rankdata(cov) for cov in cov_m])
    design = np.column_stack([np.ones(n), rc])
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    ry_res = ry - design @ by
    bs, *_ = np.linalg.lstsq(design, rs.T, rcond=None)
    rs_res = rs - (design @ bs).T
    rs_res = rs_res - rs_res.mean(axis=1, keepdims=True)
    ry_res = ry_res - ry_res.mean()
    denom = np.sqrt((rs_res ** 2).sum(axis=1) * (ry_res ** 2).sum())
    good = denom > 0
    vals = np.full(rs.shape[0], np.nan)
    vals[good] = (rs_res[good] @ ry_res) / denom[good]
    out[ok] = vals
    return out


def equal_weight(score: np.ndarray, tac: np.ndarray) -> np.ndarray:
    """Mean of within-mask z-scores. Weights are not fit to the endpoint."""
    score = np.asarray(score, dtype=float)
    tac = np.asarray(tac, dtype=float)
    if score.ndim == 1:
        score = score.reshape(1, -1)
    out = np.full_like(score, np.nan)
    mask = np.isfinite(tac)
    if int(mask.sum()) < 8 or np.nanstd(tac) == 0:
        return out
    tac_z = np.full(tac.shape, np.nan)
    tac_z[mask] = (tac[mask] - np.mean(tac[mask])) / np.std(tac[mask])
    for i in range(score.shape[0]):
        row = score[i]
        ok = mask & np.isfinite(row)
        if int(ok.sum()) < 8 or np.std(row[ok]) == 0:
            continue
        z = (row[ok] - np.mean(row[ok])) / np.std(row[ok])
        out[i, ok] = 0.5 * (z + tac_z[ok])
    return out


def gene_vector(expr: pd.DataFrame, gene: str) -> np.ndarray | None:
    if gene not in expr.index:
        return None
    return expr.loc[gene].to_numpy(dtype=float)


def estimate_immune(expr: pd.DataFrame, immune_genes: list[str], common: set[str]) -> np.ndarray:
    keep = [gene for gene in expr.index if gene in common]
    return M.ssgsea_fast(expr.loc[keep], immune_genes, alpha=M.PRIMARY_ALPHA)


def pack_scores(zmeans: dict, ss: dict) -> tuple[list[tuple], np.ndarray]:
    keys = []
    rows = []
    for size, values in zmeans.items():
        keys.append(("z-mean", None, int(size)))
        rows.append(np.asarray(values, dtype=float))
    for (alpha, size), values in ss.items():
        keys.append(("ssGSEA", float(alpha), int(size)))
        rows.append(np.asarray(values, dtype=float))
    if not rows:
        raise SystemExit("no signature scores")
    return keys, np.vstack(rows)


def correlate_block(score: np.ndarray, y: np.ndarray, covs: list[np.ndarray] | None) -> np.ndarray:
    if covs:
        return partial_rows(score, y, covs)
    return spearman_rows(score, y)


def exact_stat(score: np.ndarray, y: np.ndarray, covs: list[np.ndarray] | None) -> dict:
    if covs:
        return A.partial_spearman(score, y, covs)
    return A.spearman_pair(score, y)


def meta_from_cohorts(stats: dict[str, dict], k: int, cohorts: tuple[str, ...] | None = None) -> dict:
    use_cohorts = REQUIRED if cohorts is None else cohorts
    by_study: dict[str, list[dict]] = {}
    for cohort in use_cohorts:
        stat = stats.get(cohort)
        if stat is None or not np.isfinite(stat["rho"]):
            continue
        study = "GSE273377" if cohort.startswith("GSE273377") else cohort
        by_study.setdefault(study, []).append({"rho": stat["rho"], "n": stat["n"]})
    pieces = []
    for study in STUDIES:
        group = by_study.get(study, [])
        if not group:
            continue
        if len(group) == 1:
            pieces.append(group[0])
        else:
            combined = A.ivw(group, k=k)
            if combined["k_studies"]:
                pieces.append({"rho": combined["rho"], "n": combined["n_sum"]})
    return A.ivw(pieces, k=k)


def eligible(stats: dict[str, dict]) -> bool:
    return all(cohort in stats and np.isfinite(stats[cohort]["rho"]) for cohort in REQUIRED)


def spec_key(feature: str, method: str, alpha, size: int, covariate: str, immune: str) -> tuple:
    alpha_key = None if alpha is None or (isinstance(alpha, float) and not np.isfinite(alpha)) else float(alpha)
    return (feature, method, alpha_key, int(size), covariate, immune)


def load_inputs():
    M.download_inputs()
    signature = M.load_ranked_signature()
    ranked = list(signature["gene"])
    stromal_genes, immune_genes, common = M.load_estimate_sets()
    cohorts = M.build_cohorts(stromal_genes, common)
    return signature, ranked, immune_genes, common, cohorts


def score_cohort(cohort: M.Cohort, ranked: list[str], immune_genes, common):
    expr = cohort.expr
    samples = list(expr.columns)
    immune8, immune_used = A.zmean(expr, A.IMMUNE8, samples)
    endpoints = {
        "immune8": immune8.to_numpy(dtype=float),
        "estimate": estimate_immune(expr, immune_genes, common),
        "CD8A": gene_vector(expr, "CD8A"),
    }
    tac = gene_vector(expr, "TACSTD2")
    if tac is None or endpoints["CD8A"] is None:
        raise SystemExit(f"{cohort.name} is missing TACSTD2 or CD8A")
    keratins = []
    for gene in A.KERATIN_COVARIATES:
        values = gene_vector(expr, gene)
        if values is None:
            keratins = None
            break
        keratins.append(values)
    zmeans, ss = M.score_bank(expr, ranked)
    keys, matrix = pack_scores(zmeans, ss)
    finite = np.isfinite(cohort.purity)
    resid_expr = M.residualize(expr.loc[:, np.asarray(samples)[finite]], cohort.purity[finite])
    resid_tac = gene_vector(resid_expr, "TACSTD2")
    if resid_tac is None:
        raise SystemExit(f"{cohort.name} residual matrix has no TACSTD2")
    rz, rss = M.score_bank(resid_expr, ranked)
    rkeys, rmatrix = pack_scores(rz, rss)
    if rkeys != keys:
        raise SystemExit(f"{cohort.name} residual score keys differ from raw keys")
    resid_samples = list(resid_expr.columns)
    resid_index = [samples.index(sample) for sample in resid_samples]
    return {
        "name": cohort.name,
        "study": cohort.study,
        "n": len(samples),
        "samples": samples,
        "purity": np.asarray(cohort.purity, dtype=float),
        "purity_name": cohort.purity_name,
        "endpoints": endpoints,
        "n_immune8": len(immune_used),
        "immune_used": list(immune_used),
        "resid_tac": resid_tac,
        "tac": tac,
        "keratins": keratins,
        "keys": keys,
        "matrix": matrix,
        "resid_matrix": rmatrix,
        "resid_index": np.asarray(resid_index, dtype=int),
        "n_signature": M.present_count(expr, ranked, 221),
        "expr": expr,
    }


def covariate_sets(block: dict, index: np.ndarray | None) -> dict[str, list[np.ndarray] | None]:
    def take(values):
        values = np.asarray(values, dtype=float)
        return values if index is None else values[index]

    purity = take(block["purity"])
    out = {
        "none": None,
        "purity": [purity],
        "residual": None,
    }
    keratins = block["keratins"]
    if keratins is None:
        out["keratin"] = None
        out["purity_keratin"] = None
        out["keratin_missing"] = True
    else:
        kvecs = [take(values) for values in keratins]
        out["keratin"] = kvecs
        out["purity_keratin"] = [purity, *kvecs]
        out["keratin_missing"] = False
    return out


def accumulate(store, feature, keys, rhos, n, covariate, immune, cohort):
    for (method, alpha, size), rho in zip(keys, rhos):
        key = spec_key(feature, method, alpha, size, covariate, immune)
        store.setdefault(key, {})[cohort] = {"rho": float(rho) if np.isfinite(rho) else np.nan, "n": int(n)}


def finite_mask(y: np.ndarray, covs: list[np.ndarray] | None, extra: np.ndarray | None = None) -> np.ndarray:
    mask = np.isfinite(np.asarray(y, dtype=float))
    if covs:
        for cov in covs:
            mask &= np.isfinite(np.asarray(cov, dtype=float))
    if extra is not None:
        mask &= np.isfinite(extra)
    return mask


def slice_covs(covs: list[np.ndarray] | None, mask: np.ndarray) -> list[np.ndarray] | None:
    if not covs:
        return None
    return [np.asarray(cov, dtype=float)[mask] for cov in covs]


def sweep_aligned(store, feature, keys, score, y, covs, covariate, immune, cohort) -> None:
    """`score` and `y` are already restricted to the samples that enter the test."""
    rhos = correlate_block(score, y, covs)
    n = int(len(np.asarray(y)))
    accumulate(store, feature, keys, rhos, n, covariate, immune, cohort)
    if np.ndim(score) == 2 and not np.isfinite(score).all():
        n_vec = np.isfinite(score).sum(axis=1)
        for (method, alpha, size), n_i in zip(keys, n_vec):
            key = spec_key(feature, method, alpha, size, covariate, immune)
            store[key][cohort]["n"] = int(n_i)
            if n_i < 8:
                store[key][cohort]["rho"] = np.nan


def sweep_cohort(block: dict, store: dict) -> None:
    keys = block["keys"]
    name = block["name"]
    covsets = covariate_sets(block, None)
    n_raw = block["n"]
    for immune_name, y in block["endpoints"].items():
        if immune_name == "CD8A":
            continue
        for covariate in ("none", "purity", "keratin", "purity_keratin"):
            if covariate in ("keratin", "purity_keratin") and covsets["keratin_missing"]:
                continue
            cov = covsets[covariate]
            mask = finite_mask(y, cov)
            if int(mask.sum()) < 8 + K_COV[covariate]:
                continue
            y_m = np.asarray(y, dtype=float)[mask]
            cov_m = slice_covs(cov, mask)
            sweep_aligned(store, "signature", keys, block["matrix"][:, mask], y_m, cov_m, covariate, immune_name, name)
            tac_mask = mask & np.isfinite(block["tac"])
            if int(tac_mask.sum()) < 8 + K_COV[covariate]:
                continue
            y_t = np.asarray(y, dtype=float)[tac_mask]
            cov_t = slice_covs(cov, tac_mask)
            tac_m = block["tac"][tac_mask]
            sweep_aligned(store, "tacstd2", [("gene", None, 1)], tac_m, y_t, cov_t, covariate, immune_name, name)
            combo = equal_weight(block["matrix"][:, tac_mask], tac_m)
            sweep_aligned(store, "signature_plus_tacstd2", keys, combo, y_t, cov_t, covariate, immune_name, name)
    # Residual scores align to finite-purity samples. Endpoint stays raw.
    rindex = block["resid_index"]
    y_resid = {ename: np.asarray(values)[rindex] for ename, values in block["endpoints"].items()}
    tac_resid = np.asarray(block["resid_tac"], dtype=float)
    n_resid = int(len(rindex))
    for immune_name, y in y_resid.items():
        if immune_name == "CD8A":
            continue
        mask = finite_mask(y, None)
        if int(mask.sum()) < 8:
            continue
        sweep_aligned(
            store, "signature", keys, block["resid_matrix"][:, mask], y[mask], None, "residual", immune_name, name
        )
        tac_mask = mask & np.isfinite(tac_resid)
        if int(tac_mask.sum()) < 8:
            continue
        sweep_aligned(
            store, "tacstd2", [("gene", None, 1)], tac_resid[tac_mask], y[tac_mask], None, "residual", immune_name, name
        )
        combo = equal_weight(block["resid_matrix"][:, tac_mask], tac_resid[tac_mask])
        sweep_aligned(
            store, "signature_plus_tacstd2", keys, combo, y[tac_mask], None, "residual", immune_name, name
        )
    print(f"scored {name} n={n_raw} residual_n={n_resid} immune8_genes={block['n_immune8']}", flush=True)


def analysis_n(y: np.ndarray, covs: list[np.ndarray] | None) -> int:
    mask = np.isfinite(y)
    if covs:
        for cov in covs:
            mask &= np.isfinite(cov)
    return int(mask.sum())


def build_meta_table(store: dict) -> pd.DataFrame:
    rows = []
    for key, stats in store.items():
        feature, method, alpha, size, covariate, immune = key
        k = K_COV[covariate]
        full = meta_from_cohorts(stats, k)
        geo_cohorts = tuple(cohort for cohort in REQUIRED if cohort != "OncoSG")
        geo = meta_from_cohorts(stats, k, geo_cohorts)
        ok = eligible(stats)
        signs = [stats[cohort]["rho"] for cohort in REQUIRED if cohort in stats and np.isfinite(stats[cohort]["rho"])]
        rows.append({
            "feature": feature,
            "method": method,
            "alpha": np.nan if alpha is None else alpha,
            "size": size,
            "covariate": covariate,
            "immune": immune,
            "eligible": ok,
            "n_negative": int(sum(rho < 0 for rho in signs)),
            "n_strata": len(signs),
            "meta_rho": full["rho"],
            "meta_p": full["p"],
            "meta_i2": full["i2"],
            "meta_model": full["model"],
            "meta_n": full["n_sum"],
            "meta_k": full["k_studies"],
            "meta_ci_low": full.get("ci_low", np.nan),
            "meta_ci_high": full.get("ci_high", np.nan),
            "geo_rho": geo["rho"],
            "geo_p": geo["p"],
            "geo_i2": geo["i2"],
            "geo_model": geo["model"],
            "geo_n": geo["n_sum"],
            "abs_meta": abs(full["rho"]) if np.isfinite(full["rho"]) else np.nan,
            "abs_geo": abs(geo["rho"]) if np.isfinite(geo["rho"]) else np.nan,
        })
    return pd.DataFrame(rows)


def pick_max(metas: pd.DataFrame, mask: pd.Series, column: str = "abs_meta") -> pd.Series:
    hit = metas[mask & metas.eligible & np.isfinite(metas[column])].copy()
    if hit.empty:
        raise SystemExit("no eligible spec in this slice")
    hit["immune_rank"] = np.where(hit.immune == "immune8", 0, 1)
    hit["cov_rank"] = hit.covariate.map({"none": 0, "purity": 1, "keratin": 2, "purity_keratin": 3, "residual": 4}).fillna(9)
    hit["alpha_rank"] = np.where(hit.method == "z-mean", 0, np.where(np.isclose(hit.alpha.fillna(-1), M.PRIMARY_ALPHA), 0, 1))
    hit = hit.sort_values(
        [column, "immune_rank", "cov_rank", "alpha_rank", "size"],
        ascending=[False, True, True, True, False],
    )
    return hit.iloc[0]


def checksum(blocks: list[dict], ranked: list[str]) -> pd.DataFrame:
    """Match published PR 590 / PR 693 cells before the maximum is interpreted."""
    records = []
    by_name = {block["name"]: block for block in blocks}
    onco = by_name["OncoSG"]
    idx = {(method, alpha, size): i for i, (method, alpha, size) in enumerate(onco["keys"])}
    score = onco["matrix"][idx[("z-mean", None, 221)]]
    cd8 = exact_stat(score, onco["endpoints"]["CD8A"], None)
    imm = exact_stat(score, onco["endpoints"]["immune8"], None)
    checks = [
        ("OncoSG", "z-mean", 221, "CD8A", cd8["rho"], -0.5779846849982597, 1e-9),
        ("OncoSG", "z-mean", 221, "immune8", imm["rho"], -0.615, 0.0015),
    ]
    for cohort, expected in IMMUNE8_EXPECTED.items():
        block = by_name[cohort]
        bidx = {(method, alpha, size): i for i, (method, alpha, size) in enumerate(block["keys"])}
        values = block["matrix"][bidx[("ssGSEA", 0.75, 163)]]
        stat = exact_stat(values, block["endpoints"]["immune8"], None)
        checks.append((cohort, "ssGSEA α=0.75", 163, "immune8", stat["rho"], expected, 0.0015))
    ok = True
    for cohort, method, size, endpoint, observed, expected, tol in checks:
        match = abs(observed - expected) <= tol
        ok = ok and match
        records.append({
            "cohort": cohort,
            "method": method,
            "size": size,
            "endpoint": endpoint,
            "observed": observed,
            "expected": expected,
            "abs_diff": abs(observed - expected),
            "ok": match,
        })
        print(f"checksum {cohort} {method} {endpoint} observed {observed:+.6f} expected {expected:+.6f} ok={match}", flush=True)
    out = pd.DataFrame(records)
    if not ok:
        out.to_csv(TABLES / "locked_checks.tsv", sep="\t", index=False)
        raise SystemExit("locked ImmuneScore / CD8A cells did not match PR 590 / PR 693")
    # silence unused ranked in case of future gene-count checks
    _ = ranked
    return out


def slice_mask(metas: pd.DataFrame, feature: str | None, immune: str | None) -> pd.Series:
    mask = metas.eligible.copy()
    if feature is not None:
        mask &= metas.feature == feature
    if immune is not None:
        mask &= metas.immune == immune
    return mask


def top_rows(metas: pd.DataFrame, mask: pd.Series, n: int = 8) -> pd.DataFrame:
    hit = metas[mask & np.isfinite(metas.abs_meta)].sort_values(
        ["abs_meta", "size"], ascending=[False, False]
    )
    return hit.head(n)


def cohort_stats_for(store, row) -> dict[str, dict]:
    alpha = None if row.method == "z-mean" or not np.isfinite(row.alpha) else float(row.alpha)
    key = spec_key(row.feature, row.method, alpha, int(row["size"]), row.covariate, row.immune)
    return store[key]


def relabel_exact(blocks, row) -> list[dict]:
    """Recompute winner cells with the project Spearman so printed p-values match."""
    by_name = {block["name"]: block for block in blocks}
    stats = []
    for cohort in REQUIRED:
        block = by_name[cohort]
        if row.covariate == "residual":
            index = block["resid_index"]
            matrix = block["resid_matrix"]
            keys = block["keys"]
            y = block["endpoints"][row.immune][index]
            tac = np.asarray(block["resid_tac"], dtype=float)
            purity = block["purity"][index]
            keratins = None if block["keratins"] is None else [values[index] for values in block["keratins"]]
        else:
            matrix = block["matrix"]
            keys = block["keys"]
            y = block["endpoints"][row.immune]
            tac = block["tac"]
            purity = block["purity"]
            keratins = block["keratins"]
        covs = None
        if row.covariate == "purity":
            covs = [purity]
        elif row.covariate == "keratin":
            covs = list(keratins)
        elif row.covariate == "purity_keratin":
            covs = [purity, *keratins]
        if row.feature == "tacstd2":
            mask = finite_mask(y, covs, tac)
            stat = exact_stat(tac[mask], np.asarray(y)[mask], slice_covs(covs, mask))
        else:
            alpha = None if row.method == "z-mean" else float(row.alpha)
            idx = {(method, alpha_i, size): i for i, (method, alpha_i, size) in enumerate(keys)}
            score = matrix[idx[(row.method, alpha, int(row["size"]))]]
            if row.feature == "signature_plus_tacstd2":
                mask = finite_mask(y, covs, tac)
                mask &= np.isfinite(score)
                sliced = equal_weight(score[mask], tac[mask])[0]
                stat = exact_stat(sliced, np.asarray(y)[mask], slice_covs(covs, mask))
            else:
                mask = finite_mask(y, covs, score)
                stat = exact_stat(score[mask], np.asarray(y)[mask], slice_covs(covs, mask))
        stats.append({
            "cohort": cohort,
            "n": stat["n"],
            "rho": stat["rho"],
            "p": stat["p"],
            "ci_low": stat["ci_low"],
            "ci_high": stat["ci_high"],
        })
    return stats


def spec_label(row) -> str:
    if row.feature == "tacstd2":
        feature = "TACSTD2"
    elif row.feature == "signature_plus_tacstd2":
        feature = f"equal-weight signature + TACSTD2, {row.method}"
        if row.method == "ssGSEA":
            feature += f" α={float(row.alpha):g}"
        feature += f", size {int(row['size'])}"
    else:
        feature = row.method
        if row.method == "ssGSEA":
            feature += f" α={float(row.alpha):g}"
        feature += f", size {int(row['size'])}"
    cov = {
        "none": "unadjusted",
        "purity": "partial | purity/stroma",
        "keratin": "partial | KRT8/18/19",
        "purity_keratin": "partial | purity/stroma + KRT8/18/19",
        "residual": "gene residual | purity/stroma",
    }[row.covariate]
    immune = "ImmuneScore" if row.immune == "immune8" else "ESTIMATE immune ssGSEA"
    return f"{feature}, {cov}, vs {immune}"


def fmt_rho(r) -> str:
    return A.fmt_rho(r)


def fmt_p(p) -> str:
    return A.fmt_p(p)


def write_figures(reported: dict[str, list[dict]], metas_for_curve: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    order = list(REQUIRED)
    colors = {
        "signature": "#1b4f72",
        "tacstd2": "#922b21",
        "signature_plus_tacstd2": "#b9770e",
    }
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.6), sharex=False)
    panels = [
        ("signature", "221-gene signature"),
        ("tacstd2", "TACSTD2"),
        ("signature_plus_tacstd2", "Equal-weight sum"),
    ]
    for ax, (key, title) in zip(axes, panels):
        rows = reported[key]
        meta = reported[key + "_meta"]
        ys = np.arange(len(order) + 1)
        ax.axvline(0, color="#888888", lw=0.7)
        for y, cohort in enumerate(order):
            row = next(item for item in rows if item["cohort"] == cohort)
            ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=colors[key], lw=1.4)
            ax.plot(row["rho"], y, "o", color=colors[key], ms=5.5)
        ax.plot(meta["meta_rho"], len(order), "D", color=colors[key], ms=6)
        ax.plot([meta["meta_ci_low"], meta["meta_ci_high"]], [len(order), len(order)], color=colors[key], lw=1.6)
        ax.set_yticks(list(ys))
        ax.set_yticklabels([*order, "meta"])
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Spearman ρ")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-0.85, 0.35)
        ax.text(
            0.02,
            0.02,
            f"ρ={meta['meta_rho']:+.3f}\nI²={meta['meta_i2']:.0%}",
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
        )
    fig.suptitle("ImmuneScore (8-gene), maximum |meta ρ| inside each feature", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_immunescore_forests.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_immunescore_forests.pdf", bbox_inches="tight")
    plt.close(fig)

    styles = (
        ("z-mean", np.nan, "#1b4f72", "z-mean"),
        ("ssGSEA", 0.0, "#7f8c8d", "ssGSEA α=0"),
        ("ssGSEA", 0.25, "#148f77", "ssGSEA α=0.25"),
        ("ssGSEA", 0.75, "#b9770e", "ssGSEA α=0.75"),
        ("ssGSEA", 1.0, "#922b21", "ssGSEA α=1"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharey=True)
    for ax, covariate, title in (
        (axes[0], "none", "Unadjusted"),
        (axes[1], "keratin", "Partial | KRT8/18/19"),
    ):
        curve = metas_for_curve[
            (metas_for_curve.feature == "signature")
            & (metas_for_curve.immune == "immune8")
            & (metas_for_curve.covariate == covariate)
            & (metas_for_curve.eligible)
        ]
        ax.axhline(0, color="#888888", lw=0.6)
        for method, alpha, color, label in styles:
            if method == "z-mean":
                hit = curve[curve.method == "z-mean"].sort_values("size")
            else:
                hit = curve[(curve.method == "ssGSEA") & np.isclose(curve.alpha, alpha)].sort_values("size")
            ax.plot(hit["size"], hit["meta_rho"], color=color, lw=1.2, label=label)
        ax.set_xlabel("Prefix size (locked order)")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Meta Spearman ρ vs ImmuneScore")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Signature vs ImmuneScore, all five strata required", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_meta_rho_vs_size.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig2_meta_rho_vs_size.pdf", bbox_inches="tight")
    plt.close(fig)


def write_finding(ctx: dict) -> None:
    sig = ctx["signature_winner"]
    tac = ctx["tacstd2_winner"]
    combo = ctx["combo_winner"]
    base = ctx["baseline"]
    held = ctx["heldout"]
    lines = []
    lines.append("# PPT Part 1 — maximum |ρ| vs ImmuneScore")
    lines.append("")
    lines.append("Public LUAD bulk only. The 221 genes and their order are the PR 590 CLDN4-high malignant signature (`methods/cldn4_high_malignant_signature/tables/signature_genes.tsv`, SHA-256 `" + ctx["sha"] + "`). CLDN4 and TACSTD2 are held out of that list. Prefix length is the only size knob. No gene was added or removed because of ImmuneScore. GSE10072, GSE11969, and GSE248378 stay closed. Bulk ρ is not a spatial exclusion result and does not replace the locked CosMx ratios or the concordant-4 malignant CLDN4 % vs T/NK result.")
    lines.append("")
    lines.append("ImmuneScore is the mean of within-cohort z-scores of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, and TBX21. None of these genes are in the 221. ESTIMATE immune ssGSEA (α=0.25, common-gene filter, same engine as the stromal covariate) is reported separately and is not called ImmuneScore.")
    lines.append("")
    lines.append("## Objective")
    lines.append("")
    lines.append("Primary number: absolute DerSimonian–Laird meta-analytic Spearman versus ImmuneScore. Studies are OncoSG, GSE273377 (discovery and validation inverse-variance combined first), GSE282774, and GSE233774 tumors. Eligible specs have a finite correlation in all five strata.")
    lines.append("")
    lines.append("Three feature classes are maximized separately, then the largest of the three is the Part 1 maximum:")
    lines.append("")
    lines.append("1. Signature score (z-mean or ssGSEA; α in {0, 0.25, 0.75, 1}; prefix size 5–221).")
    lines.append("2. TACSTD2 itself. Spearman is invariant to a monotone transform, so the gene has no scale knob.")
    lines.append("3. Equal-weight mean of the within-cohort z-score of a signature score and the within-cohort z-score of TACSTD2. Weights are not fit to ImmuneScore.")
    lines.append("")
    lines.append("Covariate modes inside each class: unadjusted; partial on published PURITY or ESTIMATE StromalScore; partial on KRT8, KRT18, and KRT19; partial on purity/stroma plus those keratins; gene-level residual on purity/stroma. The p-value on a selected row is the p-value of a maximized |ρ|.")
    lines.append("")
    lines.append("## Check against published cells")
    lines.append("")
    lines.append(f"All {ctx['n_locked']} locked cells matched, including OncoSG z-mean size 221 versus CD8A at the stored PR 590 precision and the PR 693 ImmuneScore cohort rhos at ssGSEA α=0.75, size 163 (tolerance 0.0015).")
    lines.append("")
    lines.append(ctx["immune_gene_line"])
    lines.append("")
    lines.append("## Part 1 maximum")
    lines.append("")
    overall = ctx["overall"]
    lines.append(f"Largest |meta ρ| versus ImmuneScore is **{fmt_rho(overall.meta_rho)}** ({overall.meta_model}, p={fmt_p(overall.meta_p)}, I²={overall.meta_i2:.0%}, n sum={int(overall.meta_n)}, k={int(overall.meta_k)}). Spec: **{spec_label(overall)}**.")
    lines.append("")
    lines.append(cohort_sentence(ctx["overall_cohorts"]))
    lines.append("")
    lines.append("## Signature alone")
    lines.append("")
    lines.append(f"Signature-only maximum: **{fmt_rho(sig.meta_rho)}** ({sig.meta_model}, p={fmt_p(sig.meta_p)}, I²={sig.meta_i2:.0%}, n sum={int(sig.meta_n)}). Spec: {spec_label(sig)}.")
    lines.append("")
    lines.append(cohort_sentence(ctx["signature_cohorts"]))
    lines.append("")
    lines.append(f"Pre-specified baseline (z-mean, size 221, unadjusted): meta ρ {fmt_rho(base.meta_rho)} ({base.meta_model}, p={fmt_p(base.meta_p)}, I²={base.meta_i2:.0%}, n sum={int(base.meta_n)}).")
    lines.append("")
    lines.append("Same signature score across covariate modes:")
    lines.append("")
    lines.append(ctx["signature_covariate_lines"])
    lines.append("")
    lines.append(ctx["keratin_note"])
    lines.append("")
    lines.append("## TACSTD2 alone")
    lines.append("")
    lines.append(f"TACSTD2 maximum: **{fmt_rho(tac.meta_rho)}** ({tac.meta_model}, p={fmt_p(tac.meta_p)}, I²={tac.meta_i2:.0%}, n sum={int(tac.meta_n)}). Spec: {spec_label(tac)}.")
    lines.append("")
    lines.append(cohort_sentence(ctx["tacstd2_cohorts"]))
    lines.append("")
    lines.append(f"Pre-specified TACSTD2 unadjusted ImmuneScore meta ρ is {fmt_rho(ctx['tac_unadj'].meta_rho)} ({ctx['tac_unadj'].meta_model}, p={fmt_p(ctx['tac_unadj'].meta_p)}, I²={ctx['tac_unadj'].meta_i2:.0%}).")
    lines.append("")
    lines.append(ctx["tacstd2_sign_line"])
    lines.append("")
    lines.append("## Equal-weight signature + TACSTD2")
    lines.append("")
    lines.append(f"Equal-weight maximum: **{fmt_rho(combo.meta_rho)}** ({combo.meta_model}, p={fmt_p(combo.meta_p)}, I²={combo.meta_i2:.0%}, n sum={int(combo.meta_n)}). Spec: {spec_label(combo)}.")
    lines.append("")
    lines.append(cohort_sentence(ctx["combo_cohorts"]))
    lines.append("")
    lines.append(ctx["combo_vs_signature"])
    lines.append("")
    lines.append("## Held-out OncoSG")
    lines.append("")
    lines.append(f"The GEO-only maximum inside the ImmuneScore search (OncoSG not used to pick the spec) is {spec_label(ctx['geo_winner'])}. GEO meta ρ {fmt_rho(ctx['geo_winner'].geo_rho)} ({ctx['geo_winner'].geo_model}, p={fmt_p(ctx['geo_winner'].geo_p)}, I²={ctx['geo_winner'].geo_i2:.0%}, n sum={int(ctx['geo_winner'].geo_n)}).")
    lines.append("")
    lines.append(f"Applied to OncoSG: ρ={fmt_rho(held['rho'])} (p={fmt_p(held['p'])}, n={int(held['n'])}).")
    lines.append("")
    lines.append("Frozen GEO spec, all five strata (OncoSG was not used to choose it): " + cohort_sentence(ctx["heldout_cohorts"]))
    lines.append("")
    frozen = ctx["heldout_meta"]
    lines.append(
        f"That frozen spec has four-study meta ρ {fmt_rho(frozen['rho'])} ({frozen['model']}, p={fmt_p(frozen['p'])}, I²={frozen['i2']:.0%}, n sum={int(frozen['n_sum'])})."
    )
    lines.append("")
    lines.append("## ESTIMATE immune ssGSEA")
    lines.append("")
    lines.append(ctx["estimate_line"])
    lines.append("")
    lines.append("## What the 221 genes are")
    lines.append("")
    lines.append(ctx["junction_line"])
    lines.append("")
    lines.append("## Plateau")
    lines.append("")
    lines.append("Largest eligible ImmuneScore |meta ρ| values:")
    lines.append("")
    lines.append(ctx["plateau_md"])
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/ppt_part1_immunescore_max/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("Matrices download to `/tmp/cldn4sig` and are not committed. Purity covariate: OncoSG uses the cBioPortal published PURITY column. GEO uses ESTIMATE stromal ssGSEA (α=0.25). Partial correlation is Pearson of rank residuals.")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines))


def cohort_sentence(rows: list[dict]) -> str:
    parts = []
    for row in rows:
        parts.append(f"{row['cohort']} {fmt_rho(row['rho'])} (p={fmt_p(row['p'])}, n={int(row['n'])})")
    return "; ".join(parts) + "."


def md_table(frame: pd.DataFrame) -> str:
    show = frame.copy()
    lines = ["| feature | spec | covariate | meta ρ | I² | model |", "|---|---|---|---:|---:|---|"]
    for rec in show.itertuples(index=False):
        method = rec.method if rec.feature == "tacstd2" else (
            rec.method if rec.method == "z-mean" else f"ssGSEA α={float(rec.alpha):g}, size {int(rec.size)}"
        )
        if rec.feature != "tacstd2" and rec.method == "z-mean":
            method = f"z-mean, size {int(rec.size)}"
        lines.append(
            f"| {rec.feature} | {method} | {rec.covariate} | {rec.meta_rho:+.3f} | {rec.meta_i2:.0%} | {rec.meta_model} |"
        )
    return "\n".join(lines)


def covariate_lines(metas: pd.DataFrame, winner: pd.Series) -> str:
    hit = metas[
        (metas.feature == winner.feature)
        & (metas.method == winner.method)
        & (metas["size"] == winner["size"])
        & (metas.immune == "immune8")
        & (np.isclose(metas.alpha.fillna(-1), -1 if winner.method == "z-mean" else float(winner.alpha)))
    ].sort_values("abs_meta", ascending=False)
    parts = []
    for rec in hit.itertuples(index=False):
        parts.append(f"{rec.covariate}: {fmt_rho(rec.meta_rho)} ({rec.meta_model}, I²={rec.meta_i2:.0%})")
    return "; ".join(parts) + "."


def main() -> None:
    self_test()
    if "--self-test" in sys.argv:
        return
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    signature, ranked, immune_genes, common, cohorts = load_inputs()
    blocks = [score_cohort(cohort, ranked, immune_genes, common) for cohort in cohorts]
    locked = checksum(blocks, ranked)
    locked.to_csv(TABLES / "locked_checks.tsv", sep="\t", index=False)
    store: dict = {}
    for block in blocks:
        sweep_cohort(block, store)
    metas = build_meta_table(store)
    print(f"meta rows {len(metas)} eligible {int(metas.eligible.sum())}", flush=True)

    immune8 = metas.immune == "immune8"
    signature_winner = pick_max(metas, slice_mask(metas, "signature", "immune8"))
    tacstd2_winner = pick_max(metas, slice_mask(metas, "tacstd2", "immune8"))
    combo_winner = pick_max(metas, slice_mask(metas, "signature_plus_tacstd2", "immune8"))
    overall = pick_max(metas, immune8)
    geo_winner = pick_max(metas, immune8, "abs_geo")
    baseline = metas[
        (metas.feature == "signature")
        & (metas.method == "z-mean")
        & (metas["size"] == 221)
        & (metas.covariate == "none")
        & (metas.immune == "immune8")
    ].iloc[0]
    tac_unadj = metas[
        (metas.feature == "tacstd2")
        & (metas.covariate == "none")
        & (metas.immune == "immune8")
    ].iloc[0]
    estimate_winner = pick_max(metas, metas.immune == "estimate")

    def pack(row):
        cohorts_exact = relabel_exact(blocks, row)
        # Confirm the vectorized meta matches the exact Spearman to 1e-6 per cohort.
        alpha = None if row.method == "z-mean" else float(row.alpha)
        key = spec_key(row.feature, row.method, alpha, int(row["size"]), row.covariate, row.immune)
        rough = store[key]
        for item in cohorts_exact:
            if abs(item["rho"] - rough[item["cohort"]]["rho"]) > 1e-6:
                raise SystemExit(
                    f"exact rho {item['rho']} != grid {rough[item['cohort']]['rho']} for {item['cohort']} {key}"
                )
        meta = meta_from_cohorts({item["cohort"]: item for item in cohorts_exact}, K_COV[row.covariate])
        return cohorts_exact, meta

    sig_c, sig_m = pack(signature_winner)
    tac_c, tac_m = pack(tacstd2_winner)
    combo_c, combo_m = pack(combo_winner)
    all_c, all_m = pack(overall)
    held_cohorts, held_meta = pack(geo_winner)
    held = next(item for item in held_cohorts if item["cohort"] == "OncoSG")

    # Overwrite printed meta with the exact-Spearman meta so the file matches the cohort rows.
    def apply_meta(row, meta):
        row = row.copy()
        for field in ("rho", "p", "i2", "model", "n_sum", "k_studies", "ci_low", "ci_high"):
            src = {"rho": "meta_rho", "p": "meta_p", "i2": "meta_i2", "model": "meta_model", "n_sum": "meta_n", "k_studies": "meta_k", "ci_low": "meta_ci_low", "ci_high": "meta_ci_high"}[field]
            row[src] = meta[field] if field != "n_sum" else meta["n_sum"]
        return row

    signature_winner = apply_meta(signature_winner, sig_m)
    tacstd2_winner = apply_meta(tacstd2_winner, tac_m)
    combo_winner = apply_meta(combo_winner, combo_m)
    overall = apply_meta(overall, all_m)

    flat = [row for row in tac_c if abs(row["rho"]) < 0.05 or row["p"] > 0.1]
    supporting = [row for row in tac_c if row["rho"] < 0 and row not in flat]
    positive = [row for row in tac_c if row["rho"] > 0 and row not in flat]
    tac_sign = (
        f"At the TACSTD2 maximum, {len(supporting)}/5 strata are inverse. "
        + (
            "Flat strata: " + ", ".join(f"{row['cohort']} {fmt_rho(row['rho'])} (p={fmt_p(row['p'])})" for row in flat) + ". "
            if flat else ""
        )
        + (
            "Positive strata: " + ", ".join(f"{row['cohort']} {fmt_rho(row['rho'])} (p={fmt_p(row['p'])})" for row in positive) + ". "
            if positive else ""
        )
        + "Flat and positive strata stay in the meta. TACSTD2 is not a concordant cold result on these five strata."
    )
    immune_bits = []
    for block in blocks:
        missing = [gene for gene in A.IMMUNE8 if gene not in block["immune_used"]]
        if missing:
            immune_bits.append(f"{block['name']} uses {block['n_immune8']}/8 ({', '.join(missing)} absent)")
        else:
            immune_bits.append(f"{block['name']} uses 8/8")
    immune_gene_line = "ImmuneScore gene coverage: " + "; ".join(immune_bits) + "."
    same_score = metas[
        (metas.feature == signature_winner.feature)
        & (metas.method == signature_winner.method)
        & (metas["size"] == signature_winner["size"])
        & (metas.immune == "immune8")
        & np.isclose(metas.alpha.fillna(-1), float(signature_winner.alpha) if signature_winner.method != "z-mean" else -1)
        & (metas.covariate == "none")
    ].iloc[0]
    keratin_note = (
        f"Unadjusted correlation at this same score is meta ρ {fmt_rho(same_score.meta_rho)} "
        f"({same_score.meta_model}, I²={same_score.meta_i2:.0%}). "
    )
    if (
        signature_winner.method == "ssGSEA"
        and int(signature_winner["size"]) == 163
        and abs(float(signature_winner.alpha) - 0.75) < 1e-9
    ):
        keratin_note += "That unadjusted score is the PR 693 ImmuneScore maximum. "
    if signature_winner.covariate == "keratin" and abs(signature_winner.meta_rho) > abs(same_score.meta_rho):
        keratin_note += "Partial correlation on KRT8/18/19 is the covariate mode that increases |ρ| over that unadjusted score. "
    purity_same = metas[
        (metas.feature == signature_winner.feature)
        & (metas.method == signature_winner.method)
        & (metas["size"] == signature_winner["size"])
        & (metas.immune == "immune8")
        & np.isclose(metas.alpha.fillna(-1), float(signature_winner.alpha) if signature_winner.method != "z-mean" else -1)
        & (metas.covariate == "purity")
    ].iloc[0]
    if abs(purity_same.meta_rho) < abs(same_score.meta_rho):
        keratin_note += "Partial correlation on purity/stroma decreases |ρ|."

    if abs(combo_winner.meta_rho) > abs(signature_winner.meta_rho) + 1e-6:
        combo_line = (
            f"Equal-weight with TACSTD2 raises |meta ρ| relative to the signature-only maximum "
            f"({fmt_rho(combo_winner.meta_rho)} vs {fmt_rho(signature_winner.meta_rho)})."
        )
    else:
        combo_line = (
            f"Equal-weight with TACSTD2 does not raise |meta ρ| above the signature-only maximum "
            f"({fmt_rho(combo_winner.meta_rho)} vs {fmt_rho(signature_winner.meta_rho)})."
        )

    est_same = (
        f"The largest |meta ρ| versus ESTIMATE immune ssGSEA is {fmt_rho(estimate_winner.meta_rho)} "
        f"({estimate_winner.meta_model}, p={fmt_p(estimate_winner.meta_p)}, I²={estimate_winner.meta_i2:.0%}, "
        f"n sum={int(estimate_winner.meta_n)}). Spec: {spec_label(estimate_winner)}. "
        "This row is not relabeled ImmuneScore."
    )

    n_tj = int((signature.family == "TJ").sum())
    top = signature.iloc[0]
    classical = [
        "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "TJP2", "TJP3", "F11R",
        "CDH1", "EPCAM", "TACSTD2", "CGN", "MARVELD2", "CRB3", "PARD3", "CXADR",
    ]
    present_classical = [gene for gene in classical if (signature.gene == gene).any()]
    immune_overlap = sorted(set(signature.gene) & set(A.IMMUNE8))
    estimate_overlap = sorted(set(signature.gene) & set(immune_genes))
    junction = (
        f"The list is the CLDN4-tracking malignant program: concordant-4 CLDN4-high versus CLDN4-low genes "
        f"that still track CLDN4 in TCGA-LUAD after KRT8, KRT18, KRT19, and ABSOLUTE purity. "
        f"Top gene {top.gene} (family {top.family}, TCGA partial ρ vs CLDN4 {float(top.tcga_partial_rho):+.3f}). "
        f"Genes labeled TJ in the signature table: {n_tj}/221 ({', '.join(signature.loc[signature.family == 'TJ', 'gene'])}). "
        f"CLDN4 in the list: {'yes' if (signature.gene == 'CLDN4').any() else 'no'}. "
        f"TACSTD2 in the list: {'yes' if (signature.gene == 'TACSTD2').any() else 'no'}. "
        f"Overlap with ImmuneScore genes: {immune_overlap or 'none'}. "
        f"Overlap with the ESTIMATE immune set: {estimate_overlap or 'none'}. "
        f"Classical junction-panel genes present: {', '.join(present_classical) if present_classical else 'none'}. "
        "This score is not a KEGG tight-junction module."
    )

    plateau = top_rows(metas, immune8 & metas.eligible, 8)
    reported = {
        "signature": sig_c,
        "signature_meta": signature_winner,
        "tacstd2": tac_c,
        "tacstd2_meta": tacstd2_winner,
        "signature_plus_tacstd2": combo_c,
        "signature_plus_tacstd2_meta": combo_winner,
    }
    write_figures(reported, metas)

    ctx = {
        "sha": signature.attrs["sha256"],
        "n_locked": int(len(locked)),
        "overall": overall,
        "overall_cohorts": all_c,
        "signature_winner": signature_winner,
        "signature_cohorts": sig_c,
        "tacstd2_winner": tacstd2_winner,
        "tacstd2_cohorts": tac_c,
        "combo_winner": combo_winner,
        "combo_cohorts": combo_c,
        "baseline": baseline,
        "tac_unadj": tac_unadj,
        "tacstd2_sign_line": tac_sign,
        "combo_vs_signature": combo_line,
        "geo_winner": geo_winner,
        "heldout": held,
        "heldout_cohorts": held_cohorts,
        "heldout_meta": held_meta,
        "immune_gene_line": immune_gene_line,
        "estimate_line": est_same,
        "junction_line": junction,
        "signature_covariate_lines": covariate_lines(metas, signature_winner),
        "keratin_note": keratin_note,
        "plateau_md": md_table(plateau),
    }
    write_finding(ctx)

    def row_dict(row):
        out = {}
        for key, value in row.to_dict().items():
            if isinstance(value, (np.floating, float)):
                value = float(value)
                out[key] = value if np.isfinite(value) else None
            elif isinstance(value, (np.integer, int)) and not isinstance(value, bool):
                out[key] = int(value)
            elif isinstance(value, (np.bool_, bool)):
                out[key] = bool(value)
            else:
                out[key] = value
        return out

    summary = {
        "signature_sha256": signature.attrs["sha256"],
        "overall": row_dict(overall),
        "signature": row_dict(signature_winner),
        "tacstd2": row_dict(tacstd2_winner),
        "combo": row_dict(combo_winner),
        "baseline_meta_rho": float(baseline.meta_rho),
        "tacstd2_unadjusted_meta_rho": float(tac_unadj.meta_rho),
        "estimate": row_dict(estimate_winner),
        "heldout_oncosg": held,
        "geo_winner": row_dict(geo_winner),
        "n_locked_ok": int(locked.ok.sum()),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    plateau.to_csv(TABLES / "plateau_immune8.tsv", sep="\t", index=False)
    curve = metas[
        (metas.feature == "signature")
        & (metas.immune == "immune8")
        & (metas.covariate == "none")
    ]
    curve.to_csv(TABLES / "signature_unadjusted_curve.tsv", sep="\t", index=False)
    tac_grid = metas[(metas.feature == "tacstd2") & (metas.immune == "immune8")]
    tac_grid.to_csv(TABLES / "tacstd2_immune8.tsv", sep="\t", index=False)
    cohort_rows = []
    for label, rows in (("signature", sig_c), ("tacstd2", tac_c), ("signature_plus_tacstd2", combo_c), ("overall", all_c)):
        for row in rows:
            cohort_rows.append({"which": label, **row})
    pd.DataFrame(cohort_rows).to_csv(TABLES / "reported_cohorts.tsv", sep="\t", index=False)
    print(json.dumps({
        "overall": spec_label(overall),
        "overall_rho": float(overall.meta_rho),
        "signature_rho": float(signature_winner.meta_rho),
        "tacstd2_rho": float(tacstd2_winner.meta_rho),
        "combo_rho": float(combo_winner.meta_rho),
        "baseline_rho": float(baseline.meta_rho),
        "estimate_rho": float(estimate_winner.meta_rho),
        "heldout_oncosg": held["rho"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
