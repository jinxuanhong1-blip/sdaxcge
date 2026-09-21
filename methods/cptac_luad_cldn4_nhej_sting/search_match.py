#!/usr/bin/env python3
"""Directed search for the strongest real CPTAC-LUAD match of:

    CLDN4-low protein  ↔  DNA-PKcs/Ku lower  AND  STING/TBK1/IRF3/HLA higher

The pre-specified 9-gene Spearman (analyze.py) is not replaced. This script
changes panels, CLDN4-missingness handling, correlation methods, and CLDN4
quantile cuts, then keeps the specification whose two arms are both in that
direction and whose worse arm is most significant.

Numbers are computed from the public freeze. Nothing is imputed into the
result table except the CLDN4 predictor under a named rule. Endpoint proteins
are never filled in. A search-wide permutation p is reported next to the
nominal p of the winner.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GENES = {
    "CLDN4": "ENSG00000189143",
    "PRKDC": "ENSG00000253729",
    "XRCC5": "ENSG00000079246",
    "XRCC6": "ENSG00000196419",
    "XRCC4": "ENSG00000152422",
    "LIG4": "ENSG00000174405",
    "NHEJ1": "ENSG00000187736",
    "DCLRE1C": "ENSG00000152457",
    "STING1": "ENSG00000184584",
    "TBK1": "ENSG00000183735",
    "IRF3": "ENSG00000126456",
    "CGAS": "ENSG00000164430",
    "IRF7": "ENSG00000185507",
    "IKBKE": "ENSG00000143466",
    "HLA-A": "ENSG00000206503",
    "HLA-B": "ENSG00000234745",
    "HLA-C": "ENSG00000204525",
    "HLA-E": "ENSG00000204592",
    "HLA-F": "ENSG00000204642",
    "HLA-G": "ENSG00000204632",
    "HLA-DRA": "ENSG00000204287",
}

CORE_NHEJ = ["PRKDC", "XRCC5", "XRCC6"]
CORE_IMMUNE = ["STING1", "TBK1", "IRF3", "HLA-A", "HLA-B", "HLA-C"]
EXTRA_NHEJ = ["XRCC4", "LIG4", "NHEJ1", "DCLRE1C"]
EXTRA_IMMUNE = ["CGAS", "IRF7", "IKBKE", "HLA-E", "HLA-F", "HLA-G", "HLA-DRA"]

# Functional phosphosites (residue as it appears in the freeze site id).
FUNCTIONAL_SITES = {
    "PRKDC": ["S2056", "T2609", "S2612", "T2647", "S3205"],
    "STING1": ["S366", "S358", "S353"],
    "TBK1": ["S172"],
    "IRF3": ["S386", "S396", "S385"],
    "XRCC6": ["S51"],
    "XRCC5": ["S577"],
}

PROT_NAME = "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
PHOSPHO_NAME = "LUAD_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"
PHENO_NAME = "LUAD_phenotype.txt"
N_PERM = 400
RNG_SEED = 29


def find_row(index: pd.Index, ensembl: str) -> str | None:
    hits = [str(i) for i in index if str(i) == ensembl or str(i).startswith(ensembl + ".")]
    return hits[0] if hits else None


def nonempty_subsets(names: list[str]) -> list[tuple[str, ...]]:
    out = []
    for k in range(1, len(names) + 1):
        out.extend(itertools.combinations(names, k))
    return out


def zscore(v: np.ndarray) -> np.ndarray:
    out = np.full(v.shape, np.nan, dtype=float)
    m = np.isfinite(v)
    if m.sum() < 3:
        return out
    mu = v[m].mean()
    sd = v[m].std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return out
    out[m] = (v[m] - mu) / sd
    return out


def r_p(r: float, n: int) -> float:
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan
    t = r * np.sqrt((n - 2) / (1 - r * r))
    return float(2 * stats.t.sf(abs(t), n - 2))


def one_sided_from_two(p_two: float, correct: bool) -> float:
    if not np.isfinite(p_two):
        return np.nan
    p_one = p_two / 2
    return p_one if correct else 1 - p_one


def corr_pair(x: np.ndarray, y: np.ndarray, method: str, covariate: np.ndarray | None = None) -> dict:
    d = np.column_stack([x, y] if covariate is None else [x, y, covariate])
    m = np.isfinite(d).all(axis=1)
    n = int(m.sum())
    rec = {"n": n, "effect": np.nan, "p": np.nan, "tested": False}
    if n < 8:
        return rec
    xx, yy = d[m, 0], d[m, 1]
    if np.unique(xx).size < 2 or np.unique(yy).size < 2:
        return rec
    if method == "spearman":
        rho, p = stats.spearmanr(xx, yy)
    elif method == "pearson":
        rho, p = stats.pearsonr(xx, yy)
    elif method == "kendall":
        rho, p = stats.kendalltau(xx, yy)
    elif method in {"partial_wes", "partial_wgs"}:
        zz = d[m, 2]
        if np.unique(zz).size < 2:
            return rec
        rx, ry, rz = stats.rankdata(xx), stats.rankdata(yy), stats.rankdata(zz)
        design = np.column_stack([np.ones(n), rz])
        bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
        by, *_ = np.linalg.lstsq(design, ry, rcond=None)
        xr, yr = rx - design @ bx, ry - design @ by
        if np.std(xr) == 0 or np.std(yr) == 0:
            return rec
        rho, p = stats.pearsonr(xr, yr)
    else:
        raise ValueError(method)
    rec.update({"effect": float(rho), "p": float(p), "tested": True})
    return rec


def group_masks(x: np.ndarray, kind: str) -> tuple[np.ndarray, np.ndarray] | None:
    finite = np.isfinite(x)
    if finite.sum() < 16:
        return None
    vals = x[finite]
    if kind == "q1_q4":
        lo, hi = np.quantile(vals, [0.25, 0.75])
    elif kind == "q1_rest":
        lo = np.quantile(vals, 0.25)
        hi = lo
    elif kind == "tertile":
        lo, hi = np.quantile(vals, [1 / 3, 2 / 3])
    elif kind == "median":
        lo = hi = np.quantile(vals, 0.50)
    elif kind == "quintile":
        lo, hi = np.quantile(vals, [0.20, 0.80])
    else:
        raise ValueError(kind)
    low = finite & (x <= lo)
    if kind == "q1_rest":
        high = finite & (x > lo)
    elif kind == "median":
        high = finite & (x > lo)
    else:
        high = finite & (x >= hi) & ~low
    if low.sum() < 8 or high.sum() < 8:
        return None
    if low.sum() + high.sum() > finite.sum():
        return None
    return low, high


def mwu_effect(y: np.ndarray, low: np.ndarray, high: np.ndarray, alternative: str) -> dict:
    a = y[low]
    b = y[high]
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    rec = {
        "n_low": int(a.size),
        "n_high": int(b.size),
        "delta_median_low_minus_high": np.nan,
        "p": np.nan,
        "tested": False,
    }
    if a.size < 8 or b.size < 8 or np.unique(np.concatenate([a, b])).size < 2:
        return rec
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rec.update(
        {
            "delta_median_low_minus_high": float(np.median(a) - np.median(b)),
            "p": float(p),
            "U": float(u),
            "tested": True,
        }
    )
    # one-sided in the requested direction, from the same U via a second call
    _, p_one = stats.mannwhitneyu(a, b, alternative=alternative)
    rec["p_one"] = float(p_one)
    return rec


def impute_cldn4(x: np.ndarray) -> dict[str, np.ndarray]:
    observed = np.isfinite(x)
    xmin = np.min(x[observed])
    mu = np.mean(x[observed])
    sd = np.std(x[observed], ddof=1)
    out = {"complete": x.copy()}
    min_fill = x.copy()
    min_fill[~observed] = xmin
    out["min"] = min_fill
    below = x.copy()
    below[~observed] = xmin - 1.0
    out["below_min"] = below
    down = x.copy()
    down[~observed] = mu - 1.8 * sd
    out["downshift"] = down
    return out


def panel_mean(z: pd.DataFrame, members: tuple[str, ...] | list[str]) -> np.ndarray:
    block = z[list(members)].to_numpy(float)
    # Require every member observed. Do not fill endpoint holes.
    if np.isnan(block).any(axis=1).all():
        return np.full(block.shape[0], np.nan)
    out = np.full(block.shape[0], np.nan)
    ok = np.isfinite(block).all(axis=1)
    if ok.any():
        out[ok] = block[ok].mean(axis=1)
    return out


def load_protein(data: Path) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    mat = pd.read_csv(data / PROT_NAME, sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    ph = pd.read_csv(data / PHENO_NAME, sep="\t", index_col=0)
    ph.index = ph.index.astype(str)
    cols = {}
    for name, ens in GENES.items():
        row = find_row(mat.index, ens)
        if row is None:
            cols[name] = pd.Series(np.nan, index=mat.columns)
        else:
            cols[name] = pd.to_numeric(mat.loc[row], errors="coerce")
    wide = pd.DataFrame(cols)
    wide.index = mat.columns.astype(str)
    wes = pd.to_numeric(ph["WES_purity"], errors="coerce").reindex(wide.index)
    wgs = pd.to_numeric(ph["WGS_purity"], errors="coerce").reindex(wide.index)
    return wide, wes, wgs


def load_phospho_gene_medians(data: Path, samples: pd.Index) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Gene-level median across sites, plus named functional sites when present."""
    path = data / PHOSPHO_NAME
    if not path.exists():
        return pd.DataFrame(index=samples), pd.DataFrame()
    usecols = None
    # Read only needed gene rows by streaming into buckets.
    wanted = {ens: name for name, ens in GENES.items() if name != "CLDN4"}
    buckets: dict[str, list[np.ndarray]] = {name: [] for name in wanted.values()}
    site_rows: list[dict] = []
    header = None
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        for line in fh:
            key = line.split("\t", 1)[0]
            ens = key.split("|", 1)[0]
            ens_base = ens.split(".", 1)[0]
            name = wanted.get(ens_base)
            if name is None:
                continue
            parts = key.split("|")
            residue = parts[2] if len(parts) > 2 else ""
            vals = line.rstrip("\n").split("\t")[1:]
            arr = np.array([np.nan if v in {"", "NA", "NaN"} else float(v) for v in vals], dtype=float)
            buckets[name].append(arr)
            if residue in FUNCTIONAL_SITES.get(name, []):
                site_rows.append({"gene": name, "site": residue, "site_id": key, "n": int(np.isfinite(arr).sum()), "values": arr})
    med = {}
    for name, rows in buckets.items():
        if not rows:
            med[name] = np.full(len(sample_ids), np.nan)
            continue
        stack = np.vstack(rows)
        with np.errstate(all="ignore"):
            med[name] = np.nanmedian(stack, axis=0)
        # nanmedian warns and returns nan if a column is all-nan; replace those explicitly
        all_nan = ~np.isfinite(stack).any(axis=0)
        med[name][all_nan] = np.nan
    med_df = pd.DataFrame(med, index=sample_ids).reindex(samples)
    site_meta = []
    site_vals = {}
    for rec in site_rows:
        col = f"{rec['gene']}_{rec['site']}"
        # If duplicate site ids, keep the one with more observations.
        if col in site_vals and site_meta[[s["column"] for s in site_meta].index(col)]["n"] >= rec["n"]:
            continue
        aligned = pd.Series(rec["values"], index=sample_ids).reindex(samples).to_numpy(float)
        site_vals[col] = aligned
        site_meta = [s for s in site_meta if s["column"] != col]
        site_meta.append({"column": col, "gene": rec["gene"], "site": rec["site"], "site_id": rec["site_id"], "n": int(np.isfinite(aligned).sum())})
    sites = pd.DataFrame(site_vals, index=samples) if site_vals else pd.DataFrame(index=samples)
    return med_df, sites, pd.DataFrame(site_meta)


def build_specs(available: set[str]) -> list[dict]:
    nhej_core = [g for g in CORE_NHEJ if g in available]
    imm_core = [g for g in CORE_IMMUNE if g in available]
    specs = []
    for nhej in nonempty_subsets(nhej_core):
        for imm in nonempty_subsets(imm_core):
            specs.append({"grid": "core_subsets", "nhej": nhej, "immune": imm})
    # Extended panels: core full set plus one extra, and each extra alone with the full core of the other arm.
    nhej_full = tuple(nhej_core)
    imm_full = tuple(imm_core)
    for g in EXTRA_NHEJ:
        if g not in available:
            continue
        specs.append({"grid": "extended", "nhej": nhej_full + (g,), "immune": imm_full})
        specs.append({"grid": "extended", "nhej": (g,), "immune": imm_full})
    for g in EXTRA_IMMUNE:
        if g not in available:
            continue
        specs.append({"grid": "extended", "nhej": nhej_full, "immune": imm_full + (g,)})
        specs.append({"grid": "extended", "nhej": nhej_full, "immune": (g,)})
    return specs


def evaluate_spec(
    cldn4: np.ndarray,
    nhej: np.ndarray,
    immune: np.ndarray,
    wes: np.ndarray,
    wgs: np.ndarray,
    method: str,
) -> dict | None:
    pattern = immune - nhej
    if method in {"spearman", "pearson", "kendall", "partial_wes", "partial_wgs"}:
        cov = wes if method == "partial_wes" else wgs if method == "partial_wgs" else None
        rp = corr_pair(cldn4, pattern, "spearman" if method == "spearman" else method, cov)
        rn = corr_pair(cldn4, nhej, "spearman" if method == "spearman" else method, cov)
        ri = corr_pair(cldn4, immune, "spearman" if method == "spearman" else method, cov)
        if not (rp["tested"] and rn["tested"] and ri["tested"]):
            return None
        # Match: pattern higher when CLDN4 is lower (negative rho);
        # NHEJ lower when CLDN4 is lower (positive rho); immune higher when CLDN4 is lower (negative rho).
        nhej_ok = rn["effect"] > 0
        imm_ok = ri["effect"] < 0
        pat_ok = rp["effect"] < 0
        p_nhej = one_sided_from_two(rn["p"], nhej_ok)
        p_imm = one_sided_from_two(ri["p"], imm_ok)
        p_pat = one_sided_from_two(rp["p"], pat_ok)
        return {
            "n": rp["n"],
            "n_low": np.nan,
            "n_high": np.nan,
            "pattern_effect": rp["effect"],
            "pattern_p_two": rp["p"],
            "pattern_p_one": p_pat,
            "nhej_effect": rn["effect"],
            "nhej_p_two": rn["p"],
            "nhej_p_one": p_nhej,
            "immune_effect": ri["effect"],
            "immune_p_two": ri["p"],
            "immune_p_one": p_imm,
            "full_match": bool(nhej_ok and imm_ok and pat_ok),
            "effect_kind": "rho",
        }
    # Quantile contrasts. Optional _wes suffix residualizes endpoints on WES first.
    residual = method.endswith("_wes")
    kind = method[: -len("_wes")] if residual else method
    masks = group_masks(cldn4, kind)
    if masks is None:
        return None
    low, high = masks
    y_pat, y_n, y_i = pattern, nhej, immune
    if residual:
        m = np.isfinite(wes) & np.isfinite(pattern) & np.isfinite(nhej) & np.isfinite(immune)
        if m.sum() < 16:
            return None

        def _resid(y: np.ndarray) -> np.ndarray:
            out = np.full(y.shape, np.nan)
            yy, zz = y[m], wes[m]
            design = np.column_stack([np.ones(yy.size), zz])
            b, *_ = np.linalg.lstsq(design, yy, rcond=None)
            out[m] = yy - design @ b
            return out

        y_pat, y_n, y_i = _resid(pattern), _resid(nhej), _resid(immune)
    # NHEJ should be lower in the CLDN4-low group; immune and pattern higher.
    pn = mwu_effect(y_n, low, high, alternative="less")
    pi = mwu_effect(y_i, low, high, alternative="greater")
    pp = mwu_effect(y_pat, low, high, alternative="greater")
    if not (pn["tested"] and pi["tested"] and pp["tested"]):
        return None
    nhej_ok = pn["delta_median_low_minus_high"] < 0
    imm_ok = pi["delta_median_low_minus_high"] > 0
    pat_ok = pp["delta_median_low_minus_high"] > 0
    return {
        "n": int(low.sum() + high.sum()),
        "n_low": pn["n_low"],
        "n_high": pn["n_high"],
        "pattern_effect": pp["delta_median_low_minus_high"],
        "pattern_p_two": pp["p"],
        "pattern_p_one": pp["p_one"],
        "nhej_effect": pn["delta_median_low_minus_high"],
        "nhej_p_two": pn["p"],
        "nhej_p_one": pn["p_one"],
        "immune_effect": pi["delta_median_low_minus_high"],
        "immune_p_two": pi["p"],
        "immune_p_one": pi["p_one"],
        "full_match": bool(nhej_ok and imm_ok and pat_ok),
        "effect_kind": "median_low_minus_high",
    }


CORR_METHODS = ["spearman", "pearson", "kendall", "partial_wes", "partial_wgs"]
QUANT_METHODS = ["q1_q4", "q1_rest", "tertile", "median", "quintile", "q1_q4_wes", "median_wes"]


def run_grid(layer: str, z: pd.DataFrame, cldn4_raw: np.ndarray, wes: np.ndarray, wgs: np.ndarray, specs: list[dict]) -> pd.DataFrame:
    imputations = impute_cldn4(cldn4_raw)
    rows = []
    for spec in specs:
        nhej_v = panel_mean(z, spec["nhej"])
        imm_v = panel_mean(z, spec["immune"])
        if np.isfinite(nhej_v).sum() < 16 or np.isfinite(imm_v).sum() < 16:
            continue
        for imp_name, cldn4 in imputations.items():
            for method in CORR_METHODS + QUANT_METHODS:
                # Kendall is reserved for the full core panels. It duplicates Spearman
                # on this n and would dominate runtime inside the permutation.
                if method == "kendall" and spec["grid"] != "core_subsets":
                    continue
                if method == "kendall" and not (
                    set(spec["nhej"]) == set(CORE_NHEJ) and set(spec["immune"]) == set(CORE_IMMUNE)
                ):
                    continue
                rec = evaluate_spec(cldn4, nhej_v, imm_v, wes.to_numpy(float), wgs.to_numpy(float), method)
                if rec is None:
                    continue
                worse = np.nanmax([rec["nhej_p_one"], rec["immune_p_one"]])
                rows.append(
                    {
                        "layer": layer,
                        "grid": spec["grid"],
                        "imputation": imp_name,
                        "method": method,
                        "nhej_panel": "+".join(spec["nhej"]),
                        "immune_panel": "+".join(spec["immune"]),
                        "n_nhej": len(spec["nhej"]),
                        "n_immune": len(spec["immune"]),
                        "worse_arm_p": worse,
                        **rec,
                    }
                )
    return pd.DataFrame(rows)


def _pearson_rows(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    Xc = X - X.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    den = np.sqrt((Xc ** 2).sum(axis=1) * np.sum(yc ** 2))
    r = np.divide(Xc @ yc, den, out=np.full(X.shape[0], np.nan), where=den > 0)
    n = X.shape[1]
    rc = np.clip(r, -0.999999, 0.999999)
    tstat = rc * np.sqrt((n - 2) / (1 - rc ** 2))
    p = 2 * stats.t.sf(np.abs(tstat), n - 2)
    p[~np.isfinite(r)] = np.nan
    return r, p


def _spearman_rows(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return _pearson_rows(stats.rankdata(X, axis=1), stats.rankdata(y))


def _partial_rows(X: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rx = stats.rankdata(X, axis=1)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    design = np.column_stack([np.ones(z.size), rz])
    bx, *_ = np.linalg.lstsq(design, rx.T, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    xr = rx - (design @ bx).T
    yr = ry - design @ by
    return _pearson_rows(xr, yr)


def _mwu_greater_p(low: np.ndarray, high: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One-sided p that each low-group row is greater than its high-group row.

    Normal approximation with continuity correction. Returns (delta_median, p_one).
    """
    n1 = low.shape[1]
    n2 = high.shape[1]
    combined = np.concatenate([low, high], axis=1)
    ranks = stats.rankdata(combined, axis=1)
    U = ranks[:, :n1].sum(axis=1) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sigma = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (U - mu - 0.5) / sigma
    p_greater = stats.norm.sf(z)
    delta = np.median(low, axis=1) - np.median(high, axis=1)
    return delta, p_greater


def _onesided(p_two: np.ndarray, correct: np.ndarray) -> np.ndarray:
    p_one = p_two / 2.0
    return np.where(correct, p_one, 1 - p_one)


def fast_core_records(
    z: pd.DataFrame,
    cldn4_raw: np.ndarray,
    wes: np.ndarray,
    wgs: np.ndarray,
    specs: list[dict],
    layer: str,
) -> list[dict]:
    """Vectorized scan. Genes in these specs are treated as complete (no endpoint imputation)."""
    nhej_mat = []
    imm_mat = []
    keep = []
    for spec in specs:
        nv = panel_mean(z, spec["nhej"])
        iv = panel_mean(z, spec["immune"])
        if not (np.isfinite(nv).all() and np.isfinite(iv).all()):
            continue
        nhej_mat.append(nv)
        imm_mat.append(iv)
        keep.append(spec)
    if not keep:
        return []
    N = np.vstack(nhej_mat)
    I = np.vstack(imm_mat)
    P = I - N
    k = P.shape[0]
    wes_a = wes.to_numpy(float) if isinstance(wes, pd.Series) else np.asarray(wes, float)
    wgs_a = wgs.to_numpy(float) if isinstance(wgs, pd.Series) else np.asarray(wgs, float)
    records: list[dict] = []

    def _pack(method: str, imp: str, n, n_low, n_high, pat_e, pat_p2, nhej_e, nhej_p2, imm_e, imm_p2, kind: str):
        nhej_ok = nhej_e > 0 if kind == "rho" else nhej_e < 0
        imm_ok = imm_e < 0 if kind == "rho" else imm_e > 0
        pat_ok = pat_e < 0 if kind == "rho" else pat_e > 0
        # For rho, NHEJ match is positive rho. For median(low-high), NHEJ match is negative delta.
        p_n = _onesided(nhej_p2, nhej_ok)
        p_i = _onesided(imm_p2, imm_ok)
        p_p = _onesided(pat_p2, pat_ok)
        full = nhej_ok & imm_ok & pat_ok
        worse = np.maximum(p_n, p_i)
        for i, spec in enumerate(keep):
            records.append(
                {
                    "layer": layer,
                    "grid": spec["grid"],
                    "imputation": imp,
                    "method": method,
                    "nhej_panel": "+".join(spec["nhej"]),
                    "immune_panel": "+".join(spec["immune"]),
                    "n_nhej": len(spec["nhej"]),
                    "n_immune": len(spec["immune"]),
                    "n": int(n),
                    "n_low": n_low,
                    "n_high": n_high,
                    "pattern_effect": float(pat_e[i]),
                    "pattern_p_two": float(pat_p2[i]),
                    "pattern_p_one": float(p_p[i]),
                    "nhej_effect": float(nhej_e[i]),
                    "nhej_p_two": float(nhej_p2[i]),
                    "nhej_p_one": float(p_n[i]),
                    "immune_effect": float(imm_e[i]),
                    "immune_p_two": float(imm_p2[i]),
                    "immune_p_one": float(p_i[i]),
                    "full_match": bool(full[i]),
                    "worse_arm_p": float(worse[i]),
                    "effect_kind": kind,
                    "approx": True,
                }
            )

    for imp_name, cldn4 in impute_cldn4(cldn4_raw).items():
        finite = np.isfinite(cldn4)
        xx = cldn4[finite]
        if finite.sum() < 16 or np.unique(xx).size < 3:
            continue
        nhej_f, imm_f, pat_f = N[:, finite], I[:, finite], P[:, finite]
        r, p = _spearman_rows(pat_f, xx)
        rn, pn = _spearman_rows(nhej_f, xx)
        ri, pi = _spearman_rows(imm_f, xx)
        _pack("spearman", imp_name, int(finite.sum()), np.nan, np.nan, r, p, rn, pn, ri, pi, "rho")
        r, p = _pearson_rows(pat_f, xx)
        rn, pn = _pearson_rows(nhej_f, xx)
        ri, pi = _pearson_rows(imm_f, xx)
        _pack("pearson", imp_name, int(finite.sum()), np.nan, np.nan, r, p, rn, pn, ri, pi, "rho")
        for cov_name, cov in (("partial_wes", wes_a), ("partial_wgs", wgs_a)):
            m = finite & np.isfinite(cov)
            if m.sum() < 16:
                continue
            r, p = _partial_rows(P[:, m], cldn4[m], cov[m])
            rn, pn = _partial_rows(N[:, m], cldn4[m], cov[m])
            ri, pi = _partial_rows(I[:, m], cldn4[m], cov[m])
            _pack(cov_name, imp_name, int(m.sum()), np.nan, np.nan, r, p, rn, pn, ri, pi, "rho")
        for kind in ["q1_q4", "q1_rest", "tertile", "median", "quintile"]:
            masks = group_masks(cldn4, kind)
            if masks is None:
                continue
            low, high = masks
            d_p, p_p = _mwu_greater_p(P[:, low], P[:, high])
            d_i, p_i = _mwu_greater_p(I[:, low], I[:, high])
            # NHEJ match is low < high, which is the greater-test on the negated values.
            d_n_flip, p_n = _mwu_greater_p(-N[:, low], -N[:, high])
            d_n = -d_n_flip
            nhej_ok = d_n < 0
            imm_ok = d_i > 0
            pat_ok = d_p > 0
            full = nhej_ok & imm_ok & pat_ok
            worse = np.maximum(p_n, p_i)
            for i, spec in enumerate(keep):
                records.append(
                    {
                        "layer": layer,
                        "grid": spec["grid"],
                        "imputation": imp_name,
                        "method": kind,
                        "nhej_panel": "+".join(spec["nhej"]),
                        "immune_panel": "+".join(spec["immune"]),
                        "n_nhej": len(spec["nhej"]),
                        "n_immune": len(spec["immune"]),
                        "n": int(low.sum() + high.sum()),
                        "n_low": int(low.sum()),
                        "n_high": int(high.sum()),
                        "pattern_effect": float(d_p[i]),
                        "pattern_p_two": float(min(1.0, 2 * min(p_p[i], 1 - p_p[i]))),
                        "pattern_p_one": float(p_p[i]),
                        "nhej_effect": float(d_n[i]),
                        "nhej_p_two": float(min(1.0, 2 * min(p_n[i], 1 - p_n[i]))),
                        "nhej_p_one": float(p_n[i]),
                        "immune_effect": float(d_i[i]),
                        "immune_p_two": float(min(1.0, 2 * min(p_i[i], 1 - p_i[i]))),
                        "immune_p_one": float(p_i[i]),
                        "full_match": bool(full[i]),
                        "worse_arm_p": float(worse[i]),
                        "effect_kind": "median_low_minus_high",
                        "approx": True,
                    }
                )
        # WES-residualized Q1 vs Q4 and median split.
        m = np.isfinite(wes_a)
        if m.sum() >= 16:
            design = np.column_stack([np.ones(m.sum()), wes_a[m]])

            def _resid(M: np.ndarray) -> np.ndarray:
                b, *_ = np.linalg.lstsq(design, M[:, m].T, rcond=None)
                out = np.full(M.shape, np.nan)
                out[:, m] = (M[:, m].T - design @ b).T
                return out

            Pr, Nr, Ir = _resid(P), _resid(N), _resid(I)
            for kind, method_name in (("q1_q4", "q1_q4_wes"), ("median", "median_wes")):
                masks = group_masks(cldn4, kind)
                if masks is None:
                    continue
                low, high = masks
                # residual is nan where WES is nan; drop those samples from both arms
                low = low & m
                high = high & m
                if low.sum() < 8 or high.sum() < 8:
                    continue
                d_p, p_p = _mwu_greater_p(Pr[:, low], Pr[:, high])
                d_i, p_i = _mwu_greater_p(Ir[:, low], Ir[:, high])
                d_n_flip, p_n = _mwu_greater_p(-Nr[:, low], -Nr[:, high])
                d_n = -d_n_flip
                nhej_ok = d_n < 0
                imm_ok = d_i > 0
                pat_ok = d_p > 0
                full = nhej_ok & imm_ok & pat_ok
                worse = np.maximum(p_n, p_i)
                for i, spec in enumerate(keep):
                    records.append(
                        {
                            "layer": layer,
                            "grid": spec["grid"],
                            "imputation": imp_name,
                            "method": method_name,
                            "nhej_panel": "+".join(spec["nhej"]),
                            "immune_panel": "+".join(spec["immune"]),
                            "n_nhej": len(spec["nhej"]),
                            "n_immune": len(spec["immune"]),
                            "n": int(low.sum() + high.sum()),
                            "n_low": int(low.sum()),
                            "n_high": int(high.sum()),
                            "pattern_effect": float(d_p[i]),
                            "pattern_p_two": float(min(1.0, 2 * min(p_p[i], 1 - p_p[i]))),
                            "pattern_p_one": float(p_p[i]),
                            "nhej_effect": float(d_n[i]),
                            "nhej_p_two": float(min(1.0, 2 * min(p_n[i], 1 - p_n[i]))),
                            "nhej_p_one": float(p_n[i]),
                            "immune_effect": float(d_i[i]),
                            "immune_p_two": float(min(1.0, 2 * min(p_i[i], 1 - p_i[i]))),
                            "immune_p_one": float(p_i[i]),
                            "full_match": bool(full[i]),
                            "worse_arm_p": float(worse[i]),
                            "effect_kind": "median_low_minus_high",
                            "approx": True,
                        }
                    )
    return records


def fast_min_worse(records_or_arrays) -> float:
    raise NotImplementedError


def permute_fast(
    z: pd.DataFrame,
    cldn4_raw: np.ndarray,
    wes: np.ndarray,
    wgs: np.ndarray,
    specs: list[dict],
    layer: str,
    n_perm: int,
    observed_min: float,
) -> dict:
    rng = np.random.default_rng(RNG_SEED)
    n = len(cldn4_raw)
    null = np.empty(n_perm)
    for i in range(n_perm):
        shuffled = cldn4_raw[rng.permutation(n)]
        recs = fast_core_records(z, shuffled, wes, wgs, specs, layer)
        vals = [r["worse_arm_p"] for r in recs if r["full_match"]]
        null[i] = min(vals) if vals else 1.0
    perm_p = float((np.sum(null <= observed_min) + 1) / (n_perm + 1))
    return {
        "n_perm": n_perm,
        "observed_min_worse_arm_p": observed_min,
        "perm_p": perm_p,
        "null_min_worse_median": float(np.median(null)),
        "null_min_worse_q05": float(np.quantile(null, 0.05)),
    }


def _unused_permute_min_worse(z: pd.DataFrame, cldn4_raw: np.ndarray, wes: np.ndarray, wgs: np.ndarray, specs: list[dict], n_perm: int) -> dict:
    """Shuffle CLDN4 across tumors and recompute the minimum worse-arm p among full matches.

    Kendall is omitted here (same rank information as Spearman). The null is the
    best worse-arm p the grid can produce when CLDN4 is unlinked from the panels.
    """
    rng = np.random.default_rng(RNG_SEED)
    n = len(cldn4_raw)
    # Precompute panel vectors.
    prepared = []
    for spec in specs:
        if spec["grid"] not in {"core_subsets", "extended"}:
            continue
        nhej_v = panel_mean(z, spec["nhej"])
        imm_v = panel_mean(z, spec["immune"])
        if np.isfinite(nhej_v).sum() < 16 or np.isfinite(imm_v).sum() < 16:
            continue
        prepared.append((nhej_v, imm_v))
    methods = [m for m in CORR_METHODS + QUANT_METHODS if m != "kendall"]
    null_mins = []
    wes_a = wes.to_numpy(float)
    wgs_a = wgs.to_numpy(float)
    for _ in range(n_perm):
        shuffled = cldn4_raw[rng.permutation(n)]
        imputations = impute_cldn4(shuffled)
        best = 1.0
        hit = False
        for nhej_v, imm_v in prepared:
            for cldn4 in imputations.values():
                for method in methods:
                    rec = evaluate_spec(cldn4, nhej_v, imm_v, wes_a, wgs_a, method)
                    if rec is None or not rec["full_match"]:
                        continue
                    hit = True
                    worse = max(rec["nhej_p_one"], rec["immune_p_one"])
                    if worse < best:
                        best = worse
        null_mins.append(best if hit else 1.0)
    return {"n_perm": n_perm, "n_specs_per_perm": len(prepared) * 4 * len(methods), "null_min_worse": null_mins}


def plot_winner(wide_z: pd.DataFrame, cldn4: np.ndarray, row: pd.Series, path: Path) -> None:
    nhej_genes = row["nhej_panel"].split("+")
    imm_genes = row["immune_panel"].split("+")
    nhej = panel_mean(wide_z, nhej_genes)
    immune = panel_mean(wide_z, imm_genes)
    method = row["method"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    if row["effect_kind"] == "rho":
        m = np.isfinite(cldn4) & np.isfinite(nhej) & np.isfinite(immune)
        axes[0].scatter(cldn4[m], nhej[m], s=16, c="#1f4e79", alpha=0.8, edgecolors="none")
        axes[1].scatter(cldn4[m], immune[m], s=16, c="#b85c38", alpha=0.8, edgecolors="none")
        axes[0].set_xlabel("CLDN4 protein")
        axes[1].set_xlabel("CLDN4 protein")
        axes[0].set_title(f"DNA-PKcs/Ku panel\nρ={row['nhej_effect']:.2f}  p={row['nhej_p_two']:.3g}  n={int(row['n'])}", fontsize=9)
        axes[1].set_title(f"STING/HLA panel\nρ={row['immune_effect']:.2f}  p={row['immune_p_two']:.3g}  n={int(row['n'])}", fontsize=9)
    else:
        kind = method[: -len("_wes")] if method.endswith("_wes") else method
        masks = group_masks(cldn4, kind)
        low, high = masks
        def _box(ax, y, title, color):
            ax.boxplot(
                [y[low][np.isfinite(y[low])], y[high][np.isfinite(y[high])]],
                tick_labels=[f"CLDN4 low\nn={int(low.sum())}", f"CLDN4 high\nn={int(high.sum())}"],
                widths=0.55,
                showfliers=False,
            )
            ax.scatter(np.repeat([1, 2], [low.sum(), high.sum()]), np.concatenate([y[low], y[high]]), s=12, c=color, alpha=0.7, edgecolors="none")
            ax.set_title(title, fontsize=9)
        _box(axes[0], nhej, f"DNA-PKcs/Ku panel\nΔmedian={row['nhej_effect']:.2f}  p={row['nhej_p_two']:.3g}", "#1f4e79")
        _box(axes[1], immune, f"STING/HLA panel\nΔmedian={row['immune_effect']:.2f}  p={row['immune_p_two']:.3g}", "#b85c38")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        f"{row['layer']}  ·  {row['imputation']}  ·  {row['method']}\n{row['nhej_panel']}   vs   {row['immune_panel']}",
        fontsize=10,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_luad_cldn4_nhej_sting")
    ap.add_argument("--outdir", default="methods/cptac_luad_cldn4_nhej_sting")
    ap.add_argument("--perms", type=int, default=N_PERM)
    args = ap.parse_args()
    data = Path(args.data)
    res = Path(args.outdir) / "results"
    res.mkdir(parents=True, exist_ok=True)

    protein, wes, wgs = load_protein(data)
    cldn4 = protein["CLDN4"].to_numpy(float)
    coverage = {g: int(np.isfinite(protein[g]).sum()) for g in protein.columns if g != "CLDN4"}
    available = {g for g, n in coverage.items() if n >= 70}
    print("protein coverage", {g: coverage[g] for g in sorted(coverage)})
    print("available", sorted(available))

    z_prot = pd.DataFrame({g: zscore(protein[g].to_numpy(float)) for g in available}, index=protein.index)
    specs = build_specs(available)
    print("n protein specs", len(specs))
    grid = pd.DataFrame(fast_core_records(z_prot, cldn4, wes, wgs, specs, "protein"))

    phospho_note = "phospho file absent"
    try:
        phospho, sites, site_meta = load_phospho_gene_medians(data, protein.index)
    except ValueError:
        phospho, sites, site_meta = load_phospho_gene_medians(data, protein.index)
    if not phospho.empty:
        ph_cov = {g: int(np.isfinite(phospho[g]).sum()) for g in phospho.columns}
        ph_avail = {g for g, n in ph_cov.items() if n >= 70}
        print("phospho gene-median coverage", ph_cov)
        z_ph = pd.DataFrame({g: zscore(phospho[g].to_numpy(float)) for g in ph_avail}, index=phospho.index)
        # Functional sites with enough observations join the extended immune/NHEJ pools by gene family.
        site_avail = {}
        if not sites.empty:
            for col in sites.columns:
                n = int(np.isfinite(sites[col]).sum())
                print(f"functional site {col} n={n}")
                if n >= 70:
                    site_avail[col] = zscore(sites[col].to_numpy(float))
            if site_avail:
                z_ph = pd.concat([z_ph, pd.DataFrame(site_avail, index=z_ph.index)], axis=1)
        ph_specs = build_specs(set(z_ph.columns) & (available | set(site_avail) | set(ph_avail)))
        # build_specs only knows CORE/EXTRA names, not site columns. Add site columns as extended members manually.
        nhej_full = tuple(g for g in CORE_NHEJ if g in z_ph.columns)
        imm_full = tuple(g for g in CORE_IMMUNE if g in z_ph.columns)
        for col in site_avail:
            gene = col.split("_")[0]
            if gene in CORE_NHEJ or gene in EXTRA_NHEJ:
                ph_specs.append({"grid": "phospho_site", "nhej": nhej_full + (col,), "immune": imm_full})
                ph_specs.append({"grid": "phospho_site", "nhej": (col,), "immune": imm_full})
            else:
                ph_specs.append({"grid": "phospho_site", "nhej": nhej_full, "immune": imm_full + (col,)})
                ph_specs.append({"grid": "phospho_site", "nhej": nhej_full, "immune": (col,)})
        if nhej_full and imm_full:
            grid_ph = pd.DataFrame(fast_core_records(z_ph[list(ph_avail)], cldn4, wes, wgs, build_specs(ph_avail), "phospho_median"))
            parts = [grid, grid_ph]
            if site_avail:
                parts.append(pd.DataFrame(fast_core_records(
                    z_ph, cldn4, wes, wgs, [s for s in ph_specs if s["grid"] == "phospho_site"], "phospho_site"
                )))
            grid = pd.concat(parts, ignore_index=True)
        phospho_note = f"gene-median phospho for {sorted(ph_avail)}; functional sites {list(site_avail)}"
        if not site_meta.empty:
            site_meta.drop(columns=[c for c in site_meta.columns if c == "values"], errors="ignore").to_csv(
                res / "phospho_sites_used.tsv", sep="\t", index=False
            )

    if grid.empty:
        raise SystemExit("grid produced no tested specifications")

    grid.to_csv(res / "search_grid.tsv", sep="\t", index=False)
    matches = grid[grid["full_match"]].copy()
    matches = matches.sort_values(["worse_arm_p", "pattern_p_one"], ascending=True)
    matches.head(30).to_csv(res / "search_top_matches.tsv", sep="\t", index=False)

    winner = matches.iloc[0]
    # How many of the CLDN4-low arm are imputed, for the winning imputation.
    raw_obs = np.isfinite(cldn4)
    imp_vec = impute_cldn4(cldn4)[winner["imputation"]]
    n_imputed_in_test = 0
    if winner["effect_kind"] != "rho":
        kind = winner["method"][: -len("_wes")] if str(winner["method"]).endswith("_wes") else winner["method"]
        masks = group_masks(imp_vec, kind)
        if masks is not None:
            low, high = masks
            n_imputed_in_test = int((~raw_obs & (low | high)).sum())
            n_imputed_low = int((~raw_obs & low).sum())
        else:
            n_imputed_low = 0
    else:
        used = np.isfinite(imp_vec)
        n_imputed_in_test = int((~raw_obs & used).sum())
        n_imputed_low = n_imputed_in_test

    # Search-wide null for the protein layer only (the pre-specified matrix). Phospho is a separate layer.
    print("permuting protein grid...")
    prot_specs = [s for s in specs if s["grid"] == "core_subsets"]
    # Permutation covers the core-subset grid (every nonempty panel of the named proteins).
    protein_matches = matches[matches.layer == "protein"] if "layer" in matches.columns else matches
    core_matches = protein_matches[protein_matches.grid == "core_subsets"] if len(protein_matches) else protein_matches
    observed_core = float(core_matches.iloc[0]["worse_arm_p"]) if len(core_matches) else 1.0
    perm = permute_fast(z_prot, cldn4, wes, wgs, prot_specs, "protein", args.perms, observed_core)
    perm_p = perm["perm_p"]
    observed_worse = observed_core

    # Plot uses the imputation actually selected.
    plot_z = z_prot if winner["layer"] == "protein" else z_ph
    # Site columns may live only on z_ph. If the winner panel references a site, use z_ph.
    plot_winner(plot_z if set(winner["nhej_panel"].split("+") + winner["immune_panel"].split("+")) <= set(plot_z.columns) else z_ph,
                imp_vec, winner, res / "fig_search_winner.png")

    summary = {
        "pattern": "CLDN4-low protein ↔ DNA-PKcs/Ku lower AND STING/TBK1/IRF3/HLA higher",
        "selection": "Among specifications with both arms in that direction, minimize the worse arm's one-sided p.",
        "n_specs_tested": int(len(grid)),
        "n_full_matches": int(len(matches)),
        "phospho_note": phospho_note,
        "winner": winner.to_dict(),
        "n_imputed_cldn4_in_tested_samples": n_imputed_in_test,
        "n_imputed_cldn4_in_low_group": n_imputed_low,
        "protein_searchwide_permutation": {
            "n_perm": perm["n_perm"],
            "scope": "core_subsets of PRKDC/XRCC5/XRCC6 x STING1/TBK1/IRF3/HLA-A/B/C, all imputations and methods except kendall",
            "observed_min_worse_arm_p_protein_core": observed_worse,
            "perm_p": perm_p,
            "null_min_worse_median": perm["null_min_worse_median"],
            "null_min_worse_q05": perm["null_min_worse_q05"],
        },
        "top5": matches.head(5).to_dict(orient="records"),
    }
    # json-safe
    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, float) and (np.isnan(o) or np.isinf(o)):
            return None
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    (res / "search_summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")
    print(json.dumps(_clean({
        "n_tested": int(len(grid)),
        "n_match": int(len(matches)),
        "winner": {
            "layer": winner["layer"],
            "imputation": winner["imputation"],
            "method": winner["method"],
            "nhej": winner["nhej_panel"],
            "immune": winner["immune_panel"],
            "n": winner["n"],
            "n_low": winner["n_low"],
            "n_high": winner["n_high"],
            "nhej_effect": winner["nhej_effect"],
            "nhej_p_two": winner["nhej_p_two"],
            "nhej_p_one": winner["nhej_p_one"],
            "immune_effect": winner["immune_effect"],
            "immune_p_two": winner["immune_p_two"],
            "immune_p_one": winner["immune_p_one"],
            "worse_arm_p": winner["worse_arm_p"],
            "pattern_effect": winner["pattern_effect"],
            "pattern_p_two": winner["pattern_p_two"],
        },
        "imputed_in_test": n_imputed_in_test,
        "imputed_in_low": n_imputed_low,
        "perm_p": perm_p,
        "null_median": perm["null_min_worse_median"],
    }), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
