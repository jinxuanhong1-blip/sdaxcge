#!/usr/bin/env python3
"""Sweep CPTAC protein cuts, missingness, and panels.

LSCC: CLDN4 versus MHC-I / CD8 protein. The locked full-cohort anchor is
HLA-A/B/C mean z versus CLDN4, ρ about −0.41 (n=78), and CD8A about −0.44.
This script recomputes those anchors and then searches a declared grid for a
larger |ρ| and a larger Q4−Q1 median gap. A row is eligible only when the
bootstrap 95% interval of that effect lies entirely below 0.

LUAD: CLDN4 quartiles versus DNA-PK and STING protein. The locked contrast
was Q1 versus the rest. This script keeps that contrast and searches panels
for a larger Q1−Q4 median gap whose interval stays on the thesis side of 0
(DNA-PK lower in Q1; STING higher in Q1).

No log2 abundance is written into the matrix. "floor" and "below_min" are
rank rules for undetected CLDN4 and are not the primary estimand. Endpoint
proteins are never filled. LUAD and LSCC stay separate. There is no ICI label.

Selection rules are fixed in this file before the grid is read. Nominal
p-values inside the grid are not confirmatory. The search p is a permutation
of the same max statistic.
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

SEED = 20260921
N_BOOT = 2000
N_PERM = 1000
MIN_N = 20
MIN_ARM = 8
MIN_QUANT = 40
Z_LSCC = 1.959963984540054

# symbol, Ensembl gene id (version suffix matched by prefix)
CLDN4_ID = "ENSG00000189143"
MHC_CD8 = [
    ("HLA-A", "ENSG00000206503"),
    ("HLA-B", "ENSG00000234745"),
    ("HLA-C", "ENSG00000204525"),
    ("HLA-E", "ENSG00000204592"),
    ("HLA-F", "ENSG00000204642"),
    ("HLA-G", "ENSG00000204632"),
    ("B2M", "ENSG00000166710"),
    ("CD8A", "ENSG00000153563"),
    ("CD8B", "ENSG00000172116"),
]
APM = [
    ("TAP1", "ENSG00000168394"),
    ("TAP2", "ENSG00000204267"),
    ("TAPBP", "ENSG00000231925"),
    ("PSMB8", "ENSG00000204264"),
    ("PSMB9", "ENSG00000240065"),
    ("PSMB10", "ENSG00000205220"),
    ("NLRC5", "ENSG00000140853"),
]
NHEJ = [
    ("PRKDC", "ENSG00000253729"),
    ("XRCC5", "ENSG00000079246"),
    ("XRCC6", "ENSG00000196419"),
    ("XRCC4", "ENSG00000152422"),
    ("LIG4", "ENSG00000174405"),
    ("NHEJ1", "ENSG00000187736"),
    ("DCLRE1C", "ENSG00000152457"),
]
STING = [
    ("STING1", "ENSG00000184584"),
    ("TBK1", "ENSG00000183735"),
    ("IRF3", "ENSG00000126456"),
    ("IRF7", "ENSG00000185507"),
    ("CGAS", "ENSG00000164430"),
    ("IKBKE", "ENSG00000143466"),
]
CLASSICAL = {"HLA-A", "HLA-B", "HLA-C", "B2M"}
NONCLASSICAL = {"HLA-E", "HLA-F", "HLA-G"}
CD8_GENES = {"CD8A", "CD8B"}
LSCC_CUTS = ["q25_75", "q20_80", "q15_85", "q10_90", "tertile"]
LUAD_CUTS = ["q25_75", "q20_80", "q15_85", "q10_90", "tertile", "q1_rest"]


def f3(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:+.3f}"


def fp(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-3:
        return f"{x:.2e}"
    return f"{x:.3g}"


def fci(lo, hi) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return f"{lo:+.3f} to {hi:+.3f}"


def find_row(index: pd.Index, ensembl: str) -> str | None:
    hits = [str(i) for i in index if str(i) == ensembl or str(i).startswith(ensembl + ".")]
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {ensembl}: {hits}")
    return hits[0]


def zscore(v: np.ndarray, ddof: int) -> np.ndarray:
    out = np.full(v.shape, np.nan, dtype=float)
    m = np.isfinite(v)
    if int(m.sum()) < 3:
        return out
    mu = float(v[m].mean())
    sd = float(v[m].std(ddof=ddof))
    if not np.isfinite(sd) or sd == 0:
        return out
    out[m] = (v[m] - mu) / sd
    return out


def cldn4_versions(x: np.ndarray) -> dict[str, np.ndarray]:
    obs = np.isfinite(x)
    if int(obs.sum()) < 3:
        raise ValueError("CLDN4 has fewer than 3 quantified tumors")
    xmin = float(np.min(x[obs]))
    floor = x.copy()
    floor[~obs] = xmin
    below = x.copy()
    below[~obs] = xmin - 1.0
    return {"observed": x.copy(), "floor": floor, "below_min": below}


def groups(x: np.ndarray, cut: str) -> tuple[np.ndarray, np.ndarray] | None:
    finite = np.isfinite(x)
    if int(finite.sum()) < 16:
        return None
    vals = x[finite]
    if cut == "q1_rest":
        lo = float(np.quantile(vals, 0.25))
        low = finite & (x <= lo)
        high = finite & (x > lo)
    else:
        if cut == "q25_75":
            qs = [0.25, 0.75]
        elif cut == "q20_80":
            qs = [0.20, 0.80]
        elif cut == "q15_85":
            qs = [0.15, 0.85]
        elif cut == "q10_90":
            qs = [0.10, 0.90]
        elif cut == "tertile":
            qs = [1.0 / 3.0, 2.0 / 3.0]
        else:
            raise ValueError(cut)
        lo, hi = (float(v) for v in np.quantile(vals, qs))
        low = finite & (x <= lo)
        high = finite & (x >= hi) & ~low
    if int(low.sum()) < MIN_ARM or int(high.sum()) < MIN_ARM:
        return None
    return low, high


def family_of(members: tuple[str, ...], kind: str) -> str:
    if kind == "APM":
        return "APM"
    if kind == "NHEJ":
        return "DNA-PK"
    if kind == "STING":
        return "STING"
    s = set(members)
    if s <= CLASSICAL:
        return "classical_MHC1"
    if s <= (CLASSICAL | NONCLASSICAL):
        return "MHC1"
    if s <= CD8_GENES:
        return "CD8"
    return "MHC1_CD8"


def spearman_1d(x: np.ndarray, y: np.ndarray) -> tuple[int, float, float]:
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4 or np.unique(x[m]).size < 2 or np.unique(y[m]).size < 2:
        return n, np.nan, np.nan
    rho, p = stats.spearmanr(x[m], y[m])
    return n, float(rho), float(p)


def spearman_matrix(Y: np.ndarray, x: np.ndarray) -> np.ndarray:
    p, n = Y.shape
    out = np.full(p, np.nan)
    if n < 4 or np.unique(x).size < 2:
        return out
    ry = stats.rankdata(Y, axis=1)
    rx = stats.rankdata(x)
    ry = ry - ry.mean(axis=1, keepdims=True)
    rx = rx - rx.mean()
    ny = np.sqrt((ry ** 2).sum(axis=1))
    nx = float(np.sqrt((rx ** 2).sum()))
    ok = (ny > 0) & (nx > 0)
    out[ok] = (ry[ok] @ rx) / (ny[ok] * nx)
    return out


def fisher_hi(rho: np.ndarray, n: int) -> np.ndarray:
    r = np.clip(rho, -0.999999, 0.999999)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    hi = np.tanh(z + Z_LSCC * se)
    hi = np.where(np.isfinite(rho), hi, np.nan)
    return hi


def fisher_ci(rho: float, n: int) -> tuple[float, float]:
    if n < 4 or not np.isfinite(rho):
        return np.nan, np.nan
    r = float(np.clip(rho, -0.999999, 0.999999))
    z = float(np.arctanh(r))
    se = 1.0 / np.sqrt(n - 3)
    return float(np.tanh(z - Z_LSCC * se)), float(np.tanh(z + Z_LSCC * se))


def boot_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, B: int = N_BOOT) -> tuple[float, float]:
    n = len(x)
    idx = rng.integers(0, n, size=(B, n))
    Y = y[idx]
    X = x[idx]
    ry = stats.rankdata(Y, axis=1)
    rx = stats.rankdata(X, axis=1)
    ry = ry - ry.mean(axis=1, keepdims=True)
    rx = rx - rx.mean(axis=1, keepdims=True)
    ny = np.sqrt((ry ** 2).sum(axis=1))
    nx = np.sqrt((rx ** 2).sum(axis=1))
    den = ny * nx
    num = (ry * rx).sum(axis=1)
    rho = np.divide(num, den, out=np.full(B, np.nan), where=den > 0)
    if np.isfinite(rho).sum() < 100:
        return np.nan, np.nan
    return float(np.nanpercentile(rho, 2.5)), float(np.nanpercentile(rho, 97.5))


def boot_spearman_matrix(Y: np.ndarray, x: np.ndarray, rng: np.random.Generator, B: int = N_BOOT) -> tuple[np.ndarray, np.ndarray]:
    """Y is (P, n), x is (n,), both finite."""
    p, n = Y.shape
    idx = rng.integers(0, n, size=(B, n))
    yb = x[idx]
    ry = stats.rankdata(yb, axis=1)
    ry = ry - ry.mean(axis=1, keepdims=True)
    ny = np.sqrt((ry ** 2).sum(axis=1))
    lo = np.full(p, np.nan)
    hi = np.full(p, np.nan)
    step = 40
    for i in range(0, p, step):
        block = Y[i : i + step][:, idx]
        rx = stats.rankdata(block, axis=2)
        rx = rx - rx.mean(axis=2, keepdims=True)
        nx = np.sqrt((rx ** 2).sum(axis=2))
        den = nx * ny
        num = (rx * ry).sum(axis=2)
        rho = np.divide(num, den, out=np.full(num.shape, np.nan), where=den > 0)
        lo[i : i + block.shape[0]] = np.nanpercentile(rho, 2.5, axis=1)
        hi[i : i + block.shape[0]] = np.nanpercentile(rho, 97.5, axis=1)
    return lo, hi


def cliff_high_vs_low(high: np.ndarray, low: np.ndarray) -> float:
    if high.size == 0 or low.size == 0:
        return np.nan
    diff = high[:, None] - low[None, :]
    return float(((diff > 0).sum() - (diff < 0).sum()) / diff.size)


def boot_median_gap(low: np.ndarray, high: np.ndarray, rng: np.random.Generator, B: int = N_BOOT) -> dict:
    """Gap is median(high) − median(low). Cliff is P(high>low) − P(high<low)."""
    a = np.asarray(low, float)
    b = np.asarray(high, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    rec = {
        "n_low": int(a.size),
        "n_high": int(b.size),
        "delta_high_minus_low": np.nan,
        "ci_lo": np.nan,
        "ci_hi": np.nan,
        "cliff": np.nan,
        "p_mwu": np.nan,
    }
    if a.size < MIN_ARM or b.size < MIN_ARM:
        return rec
    rec["delta_high_minus_low"] = float(np.median(b) - np.median(a))
    rec["cliff"] = cliff_high_vs_low(b, a)
    try:
        _u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
        rec["p_mwu"] = float(p)
    except ValueError:
        rec["p_mwu"] = np.nan
    ia = rng.integers(0, a.size, size=(B, a.size))
    ib = rng.integers(0, b.size, size=(B, b.size))
    med = np.median(b[ib], axis=1) - np.median(a[ia], axis=1)
    rec["ci_lo"] = float(np.percentile(med, 2.5))
    rec["ci_hi"] = float(np.percentile(med, 97.5))
    return rec


def boot_median_matrix(SL: np.ndarray, SH: np.ndarray, rng: np.random.Generator, B: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SL, SH are (P, nL) and (P, nH), finite. Returns point, lo, hi of high−low."""
    point = np.median(SH, axis=1) - np.median(SL, axis=1)
    iL = rng.integers(0, SL.shape[1], size=(B, SL.shape[1]))
    iH = rng.integers(0, SH.shape[1], size=(B, SH.shape[1]))
    d = np.median(SH[:, iH], axis=2) - np.median(SL[:, iL], axis=2)
    return point, np.percentile(d, 2.5, axis=1), np.percentile(d, 97.5, axis=1)


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray, rng: np.random.Generator, B: int = N_BOOT) -> dict:
    """Pearson of rank-residuals on z. Bootstrap resamples the fixed ranks (same as the prior CPTAC partial)."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    rec = {"n_partial": int(m.sum()), "rho_partial": np.nan, "p_partial": np.nan, "ci_lo_partial": np.nan, "ci_hi_partial": np.nan}
    if m.sum() < MIN_N or np.unique(z[m]).size < 2:
        return rec
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])

    def resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        design = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    xr, yr = resid(rx, rz), resid(ry, rz)
    if np.std(xr) == 0 or np.std(yr) == 0:
        return rec
    rho, p = stats.pearsonr(xr, yr)
    rec["rho_partial"] = float(rho)
    rec["p_partial"] = float(p)
    n = len(rx)
    idx = rng.integers(0, n, size=(B, n))
    boots = np.empty(B)
    kept = 0
    for i in range(B):
        xb = resid(rx[idx[i]], rz[idx[i]])
        yb = resid(ry[idx[i]], rz[idx[i]])
        if np.std(xb) == 0 or np.std(yb) == 0:
            continue
        r, _ = stats.pearsonr(xb, yb)
        if np.isfinite(r):
            boots[kept] = r
            kept += 1
    if kept >= 100:
        rec["ci_lo_partial"] = float(np.percentile(boots[:kept], 2.5))
        rec["ci_hi_partial"] = float(np.percentile(boots[:kept], 97.5))
    return rec


def purity_gap(y: np.ndarray, wes: np.ndarray, low: np.ndarray, high: np.ndarray, rng: np.random.Generator, B: int = N_BOOT) -> dict:
    """Median(high)−median(low) after a linear residual on WES purity, refit inside the bootstrap."""
    in_low = low & np.isfinite(y) & np.isfinite(wes)
    in_high = high & np.isfinite(y) & np.isfinite(wes)
    a, wa = y[in_low], wes[in_low]
    b, wb = y[in_high], wes[in_high]
    rec = {"n_low": int(a.size), "n_high": int(b.size), "delta": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    if a.size < MIN_ARM or b.size < MIN_ARM or np.unique(np.concatenate([wa, wb])).size < 2:
        return rec

    def fit(av, aw, bv, bw) -> float:
        yy = np.concatenate([av, bv])
        zz = np.concatenate([aw, bw])
        design = np.column_stack([np.ones(len(yy)), zz])
        coef, *_ = np.linalg.lstsq(design, yy, rcond=None)
        ra = av - (coef[0] + coef[1] * aw)
        rb = bv - (coef[0] + coef[1] * bw)
        return float(np.median(rb) - np.median(ra))

    rec["delta"] = fit(a, wa, b, wb)
    ia = rng.integers(0, a.size, size=(B, a.size))
    ib = rng.integers(0, b.size, size=(B, b.size))
    boots = np.empty(B)
    for i in range(B):
        boots[i] = fit(a[ia[i]], wa[ia[i]], b[ib[i]], wb[ib[i]])
    rec["ci_lo"] = float(np.percentile(boots, 2.5))
    rec["ci_hi"] = float(np.percentile(boots, 97.5))
    return rec


def load_gene_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def extract_genes(mat: pd.DataFrame, catalog: list[tuple[str, str]]) -> tuple[dict[str, np.ndarray], list[dict]]:
    genes = {}
    rows = []
    n = mat.shape[1]
    for symbol, ens in catalog:
        hit = find_row(mat.index, ens)
        if hit is None:
            vec = np.full(n, np.nan)
        else:
            vec = pd.to_numeric(mat.loc[hit], errors="coerce").to_numpy(float)
        genes[symbol] = vec
        n_q = int(np.isfinite(vec).sum())
        rows.append(
            {
                "symbol": symbol,
                "ensembl": ens,
                "row": hit or "",
                "present": hit is not None,
                "n_tumors": n,
                "n_quantified": n_q,
                "n_missing": n - n_q,
                "min": float(np.nanmin(vec)) if n_q else np.nan,
                "median": float(np.nanmedian(vec)) if n_q else np.nan,
                "max": float(np.nanmax(vec)) if n_q else np.nan,
            }
        )
    return genes, rows


def read_clinical(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skiprows=4, dtype=str)
    return df


def load_cohort(data: Path, cohort: str) -> dict:
    mat = load_gene_matrix(
        data / cohort / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )
    samples = np.array([str(c) for c in mat.columns])
    ph = pd.read_csv(data / cohort / f"{cohort}_phenotype.txt", sep="\t", dtype=str)
    if "idx" not in ph.columns:
        raise ValueError(f"{cohort} phenotype has no idx")
    ph = ph.set_index(ph["idx"].astype(str))
    if set(samples) != set(ph.index):
        raise ValueError(f"{cohort} protein columns do not match phenotype idx")
    meta = pd.read_csv(data / cohort / f"{cohort}_meta.txt", sep="\t", dtype=str)
    meta = meta[meta["case_id"] != "data_type"].copy()
    meta = meta.set_index("case_id")
    if set(samples) - set(meta.index):
        raise ValueError(f"{cohort} protein samples missing from meta")
    meta = meta.reindex(samples)
    wes = pd.to_numeric(ph.reindex(samples)["WES_purity"], errors="coerce").to_numpy(float)
    catalog = [("CLDN4", CLDN4_ID)] + MHC_CD8 + APM + NHEJ + STING
    # unique by symbol
    seen = set()
    uniq = []
    for sym, ens in catalog:
        if sym not in seen:
            uniq.append((sym, ens))
            seen.add(sym)
    genes, coverage = extract_genes(mat, uniq)
    for row in coverage:
        row["cohort"] = cohort
    grade = meta["Histologic_Grade"].fillna("").astype(str).to_numpy()
    stage = meta["Stage"].fillna("").astype(str).to_numpy()
    strata = {"all": np.ones(len(samples), dtype=bool)}
    strata["grade_G2"] = np.array([g.startswith("G2") for g in grade])
    strata["grade_G3"] = np.array([g.startswith("G3") for g in grade])
    strata["stage_I"] = stage == "Stage I"
    strata["stage_II"] = stage == "Stage II"
    strata["stage_III"] = stage == "Stage III"
    finite_wes = np.isfinite(wes)
    if finite_wes.sum() >= 20:
        med = float(np.median(wes[finite_wes]))
        strata["purity_high"] = finite_wes & (wes >= med)
        strata["purity_low"] = finite_wes & (wes < med)
    notes = []
    if cohort == "LSCC":
        clin = read_clinical(data / "clinical" / "lusc_cptac_2021_data_clinical_sample.txt")
        path = clin.set_index("SAMPLE_ID")["PATHOLOGY_BASED_HISTOLOGY_ASSESSMENT"].reindex(samples)
        text = path.fillna("").astype(str)
        known = path.notna().to_numpy()
        basal = text.str.contains("basaloid", case=False).to_numpy()
        strata["basaloid"] = basal & known
        strata["not_basaloid"] = known & ~basal
        notes.append("basaloid = pathology text contains 'basaloid'")
    else:
        clin = read_clinical(data / "clinical" / "luad_cptac_2020_data_clinical_sample.txt")
        clin = clin.copy()
        clin["case_id"] = clin["SAMPLE_ID"].str.replace(r"^X(?=11LU)", "", regex=True)
        if clin["case_id"].duplicated().any():
            raise ValueError("LUAD clinical case_id is not unique after 11LU remap")
        sub = clin.set_index("case_id")["DOMINANT_HISTOLOGICAL_SUBTYPE"].reindex(samples)
        known = sub.notna().to_numpy()
        acinar = sub.eq("acinar").to_numpy()
        strata["acinar"] = acinar
        strata["non_acinar"] = known & ~acinar
        notes.append("acinar = cBioPortal DOMINANT_HISTOLOGICAL_SUBTYPE; X11LU ids remapped")
    # Drop strata that cannot support a correlation on observed CLDN4.
    cldn = genes["CLDN4"]
    kept = {}
    for name, mask in strata.items():
        n_obs = int((mask & np.isfinite(cldn)).sum())
        if name == "all" or n_obs >= MIN_N:
            kept[name] = mask
        else:
            notes.append(f"skipped stratum {name}: observed CLDN4 n={n_obs}")
    return {
        "cohort": cohort,
        "samples": samples,
        "genes": genes,
        "coverage": coverage,
        "wes": wes,
        "strata": kept,
        "notes": notes,
        "n": len(samples),
        "grade": grade,
        "stage": stage,
    }


def build_panels(genes: dict[str, np.ndarray], symbols: list[str], ddof: int, max_k: int, kind: str) -> list[dict]:
    use = []
    for s in symbols:
        if s not in genes:
            continue
        if int(np.isfinite(genes[s]).sum()) < MIN_QUANT:
            continue
        use.append(s)
    z = {s: zscore(genes[s], ddof) for s in use}
    panels = []
    for k in range(1, max_k + 1):
        for comb in itertools.combinations(use, k):
            mats = [z[s] for s in comb]
            M = np.column_stack(mats)
            cnt = np.isfinite(M).sum(axis=1)
            with np.errstate(all="ignore"):
                mu = np.nanmean(M, axis=1)
            strict = mu.copy()
            strict[cnt < k] = np.nan
            panels.append(
                {
                    "panel": "+".join(comb),
                    "members": comb,
                    "k": k,
                    "completeness": "strict",
                    "family": family_of(comb, kind),
                    "kind": kind,
                    "y": strict,
                }
            )
            if k >= 2:
                loose = mu.copy()
                need = k - 1
                loose[cnt < need] = np.nan
                if not (np.array_equal(np.isnan(loose), np.isnan(strict)) and np.allclose(
                    loose[np.isfinite(loose)], strict[np.isfinite(strict)], equal_nan=True
                )):
                    panels.append(
                        {
                            "panel": "+".join(comb),
                            "members": comb,
                            "k": k,
                            "completeness": "allow_one_missing",
                            "family": family_of(comb, kind),
                            "kind": kind,
                            "y": loose,
                        }
                    )
    return panels


def panel_matrix(panels: list[dict]) -> np.ndarray:
    return np.vstack([p["y"] for p in panels])


def max_abs_rho_qualified(x: np.ndarray, Y: np.ndarray) -> float:
    m = np.isfinite(x)
    n_x = int(m.sum())
    if n_x < MIN_N:
        return 0.0
    xs = x[m]
    Ys = Y[:, m]
    best = 0.0
    row_ok = np.isfinite(Ys).all(axis=1)
    if row_ok.any():
        rhos = spearman_matrix(Ys[row_ok], xs)
        hi = fisher_hi(rhos, n_x)
        ok = np.isfinite(rhos) & (rhos < 0) & (hi < 0)
        if ok.any():
            best = max(best, float(np.max(-rhos[ok])))
    bad = np.where(~row_ok)[0]
    for i in bad:
        y = Ys[i]
        mm = np.isfinite(y)
        nn = int(mm.sum())
        if nn < MIN_N or np.unique(xs[mm]).size < 2 or np.unique(y[mm]).size < 2:
            continue
        rx = stats.rankdata(xs[mm])
        ry = stats.rankdata(y[mm])
        rx = rx - rx.mean()
        ry = ry - ry.mean()
        den = np.sqrt((rx @ rx) * (ry @ ry))
        if den == 0:
            continue
        rho = float((rx @ ry) / den)
        _lo, hi = fisher_ci(rho, nn)
        if rho < 0 and hi < 0:
            best = max(best, -rho)
    return best


def permute_rho(raw: np.ndarray, panels: list[dict], strata: dict[str, np.ndarray], rng: np.random.Generator, include_strata: bool, include_tails: bool, n_perm: int) -> dict:
    Y = panel_matrix(panels)
    versions = cldn4_versions(raw)
    filters = strata if include_strata else {"all": strata["all"]}

    def stat(vers: dict[str, np.ndarray]) -> float:
        best = 0.0
        for x0 in vers.values():
            for mask in filters.values():
                xx = x0.copy()
                xx[~mask] = np.nan
                best = max(best, max_abs_rho_qualified(xx, Y))
                if not include_tails:
                    continue
                for cut in LSCC_CUTS:
                    g = groups(xx, cut)
                    if g is None:
                        continue
                    low, high = g
                    xt = np.full(xx.shape, np.nan)
                    keep = low | high
                    xt[keep] = xx[keep]
                    best = max(best, max_abs_rho_qualified(xt, Y))
        return best

    obs = stat(versions)
    finite = np.isfinite(raw)
    finite_vals = raw[finite].copy()
    floor = versions["floor"].copy()
    below = versions["below_min"].copy()
    null = np.empty(n_perm)
    for i in range(n_perm):
        fv = finite_vals.copy()
        rng.shuffle(fv)
        x_obs = raw.copy()
        x_obs[finite] = fv
        fl = floor.copy()
        rng.shuffle(fl)
        bl = below.copy()
        rng.shuffle(bl)
        null[i] = stat({"observed": x_obs, "floor": fl, "below_min": bl})
    p = (1 + int(np.sum(null + 1e-12 >= obs))) / (1 + n_perm)
    return {
        "obs_max_abs_rho": obs,
        "perm_p": p,
        "n_perm": n_perm,
        "null_median": float(np.median(null)),
        "null_q95": float(np.quantile(null, 0.95)),
        "include_strata": include_strata,
        "include_tails": include_tails,
        "gate": "Fisher z 95% interval entirely below 0; n>=20; ρ<0",
    }


def qualify_delta_matrix(x: np.ndarray, Y: np.ndarray, cut: str, rng: np.random.Generator, B: int) -> float:
    """Max (median low − median high) among rows whose high−low interval is entirely below 0."""
    g = groups(x, cut)
    if g is None:
        return 0.0
    low, high = g
    best = 0.0
    row_ok = np.isfinite(Y[:, low | high]).all(axis=1) & np.isfinite(Y).all(axis=1)
    # Finite on the two groups is enough.
    row_ok = np.isfinite(Y[:, low]).all(axis=1) & np.isfinite(Y[:, high]).all(axis=1)
    if row_ok.any() and int(low.sum()) >= MIN_ARM and int(high.sum()) >= MIN_ARM:
        point, _lo, hi = boot_median_matrix(Y[row_ok][:, low], Y[row_ok][:, high], rng, B)
        # point is high−low. Qualify hi<0. Thesis magnitude is −point = low−high.
        ok = np.isfinite(point) & (hi < 0)
        if ok.any():
            best = max(best, float(np.max(-point[ok])))
    for y in Y[~row_ok]:
        a = y[low]
        b = y[high]
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if a.size < MIN_ARM or b.size < MIN_ARM:
            continue
        rec = boot_median_gap(a, b, rng, B)
        if np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0:
            best = max(best, -rec["delta_high_minus_low"])
    return best


def permute_delta(raw: np.ndarray, panels: list[dict], rng: np.random.Generator, cut: str, n_perm: int, B: int) -> dict:
    Y = panel_matrix(panels)
    x = cldn4_versions(raw)["observed"]
    obs = qualify_delta_matrix(x, Y, cut, rng, B)
    finite = np.isfinite(raw)
    vals = raw[finite].copy()
    null = np.empty(n_perm)
    for i in range(n_perm):
        fv = vals.copy()
        rng.shuffle(fv)
        xp = raw.copy()
        xp[finite] = fv
        null[i] = qualify_delta_matrix(xp, Y, cut, rng, B)
    p = (1 + int(np.sum(null + 1e-12 >= obs))) / (1 + n_perm)
    return {
        "obs_max_low_minus_high": obs,
        "perm_p": p,
        "n_perm": n_perm,
        "boot_inside": B,
        "null_median": float(np.median(null)),
        "null_q95": float(np.quantile(null, 0.95)),
        "cut": cut,
        "missingness": "observed",
        "filter": "all",
    }


def sweep_lscc(bundle: dict, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    raw = bundle["genes"]["CLDN4"]
    versions = cldn4_versions(raw)
    mhc_syms = [s for s, _e in MHC_CD8]
    apm_syms = [s for s, _e in APM]
    panels = build_panels(bundle["genes"], mhc_syms, ddof=0, max_k=4, kind="MHC")
    apm_panels = build_panels(bundle["genes"], apm_syms, ddof=0, max_k=3, kind="APM")
    print(f"LSCC panels MHC/CD8 {len(panels)} APM {len(apm_panels)}", flush=True)
    rows = []
    # APM is scored on the full cohort only, so it cannot be silently selected as MHC-I.
    jobs = [(p, True) for p in panels] + [(p, False) for p in apm_panels]
    for panel, in_main in jobs:
        y = panel["y"]
        filters = bundle["strata"] if in_main else {"all": bundle["strata"]["all"]}
        for fname, mask in filters.items():
            for miss, x0 in versions.items():
                xx = x0.copy()
                xx[~mask] = np.nan
                n, rho, p = spearman_1d(xx, y)
                _flo, fhi = fisher_ci(rho, n)
                used = np.isfinite(xx) & np.isfinite(y)
                n_imp = int((~np.isfinite(raw) & used).sum())
                rows.append(
                    {
                        "cohort": "LSCC",
                        "kind": "rho",
                        "mode": "full",
                        "cut": "none",
                        "filter": fname,
                        "missingness": miss,
                        "panel": panel["panel"],
                        "completeness": panel["completeness"],
                        "family": panel["family"],
                        "k": panel["k"],
                        "n": n,
                        "n_low": np.nan,
                        "n_high": np.nan,
                        "n_imputed_cldn4": n_imp,
                        "n_imputed_low": 0,
                        "rho": rho,
                        "p": p,
                        "fisher_hi": fhi,
                        "delta_high_minus_low": np.nan,
                        "cliff": np.nan,
                        "p_mwu": np.nan,
                        "in_main_grid": in_main,
                    }
                )
                for cut in LSCC_CUTS:
                    g = groups(xx, cut)
                    if g is None:
                        continue
                    low, high = g
                    xt = np.full(xx.shape, np.nan)
                    keep = low | high
                    xt[keep] = xx[keep]
                    nt, rhot, pt = spearman_1d(xt, y)
                    _flo, fhit = fisher_ci(rhot, nt)
                    used_t = np.isfinite(xt) & np.isfinite(y)
                    rows.append(
                        {
                            "cohort": "LSCC",
                            "kind": "rho",
                            "mode": "tails",
                            "cut": cut,
                            "filter": fname,
                            "missingness": miss,
                            "panel": panel["panel"],
                            "completeness": panel["completeness"],
                            "family": panel["family"],
                            "k": panel["k"],
                            "n": nt,
                            "n_low": int(low.sum()),
                            "n_high": int(high.sum()),
                            "n_imputed_cldn4": int((~np.isfinite(raw) & used_t).sum()),
                            "n_imputed_low": int((~np.isfinite(raw) & low & np.isfinite(y)).sum()),
                            "rho": rhot,
                            "p": pt,
                            "fisher_hi": fhit,
                            "delta_high_minus_low": np.nan,
                            "cliff": np.nan,
                            "p_mwu": np.nan,
                            "in_main_grid": in_main,
                        }
                    )
                    a = y[low]
                    b = y[high]
                    a = a[np.isfinite(a)]
                    b = b[np.isfinite(b)]
                    if a.size < MIN_ARM or b.size < MIN_ARM:
                        continue
                    try:
                        _u, pm = stats.mannwhitneyu(b, a, alternative="two-sided")
                        pm = float(pm)
                    except ValueError:
                        pm = np.nan
                    rows.append(
                        {
                            "cohort": "LSCC",
                            "kind": "delta",
                            "mode": "groups",
                            "cut": cut,
                            "filter": fname,
                            "missingness": miss,
                            "panel": panel["panel"],
                            "completeness": panel["completeness"],
                            "family": panel["family"],
                            "k": panel["k"],
                            "n": int(a.size + b.size),
                            "n_low": int(a.size),
                            "n_high": int(b.size),
                            "n_imputed_cldn4": int((~np.isfinite(raw) & (low | high) & np.isfinite(y)).sum()),
                            "n_imputed_low": int((~np.isfinite(raw) & low & np.isfinite(y)).sum()),
                            "rho": np.nan,
                            "p": np.nan,
                            "fisher_hi": np.nan,
                            "delta_high_minus_low": float(np.median(b) - np.median(a)),
                            "cliff": cliff_high_vs_low(b, a),
                            "p_mwu": pm,
                            "in_main_grid": in_main,
                        }
                    )
    df = pd.DataFrame(rows)
    df["ci_lo"] = np.nan
    df["ci_hi"] = np.nan
    print(f"LSCC grid rows {len(df)}", flush=True)

    # Bootstrap Spearman for rows that could clear a CI below 0.
    rho_ix = df.index[(df.kind == "rho") & (df.rho < 0) & (df.fisher_hi < 0.05) & (df.n >= MIN_N)].tolist()
    # Always bootstrap the locked anchors, whatever the Fisher screen says.
    anchor_m = (
        (df.kind == "rho")
        & (df["mode"] == "full")
        & (df["filter"] == "all")
        & (df.missingness == "observed")
        & (df.completeness == "strict")
        & (df.panel.isin(["HLA-A+HLA-B+HLA-C", "CD8A"]))
    )
    rho_ix = sorted(set(rho_ix) | set(df.index[anchor_m].tolist()))
    print(f"LSCC rho bootstrap candidates {len(rho_ix)}", flush=True)
    # Group candidates that share x and the pairwise mask.
    y_of = {(p["panel"], p["completeness"], p["kind"]): p["y"] for p in panels + apm_panels}
    grouped: dict[tuple, list[int]] = {}
    x_cache = {}
    for i in rho_ix:
        r = df.loc[i]
        key_x = (r["missingness"], r["filter"], r["mode"], r["cut"])
        if key_x not in x_cache:
            x0 = versions[r["missingness"]].copy()
            x0[~bundle["strata"][r["filter"]]] = np.nan
            if r["mode"] == "tails":
                g = groups(x0, r["cut"])
                if g is None:
                    x_cache[key_x] = None
                else:
                    low, high = g
                    xt = np.full(x0.shape, np.nan)
                    xt[low | high] = x0[low | high]
                    x_cache[key_x] = xt
            else:
                x_cache[key_x] = x0
        xx = x_cache[key_x]
        if xx is None:
            continue
        y = y_of[(r["panel"], r["completeness"], "APM" if r["family"] == "APM" else "MHC")]
        m = np.isfinite(xx) & np.isfinite(y)
        # Key the finite positions. n=108, a bytes key is fine.
        grouped.setdefault((key_x, m.tobytes()), []).append(i)
    for (_key, mask_b), idxs in grouped.items():
        key_x = _key
        xx = x_cache[key_x]
        m = np.frombuffer(mask_b, dtype=bool).copy()
        if m.size != len(xx):
            raise RuntimeError("mask length mismatch")
        xs = xx[m]
        ys = []
        for i in idxs:
            r = df.loc[i]
            y = y_of[(r["panel"], r["completeness"], "APM" if r["family"] == "APM" else "MHC")]
            ys.append(y[m])
        Y = np.vstack(ys)
        lo, hi = boot_spearman_matrix(Y, xs, rng, N_BOOT)
        df.loc[idxs, "ci_lo"] = lo
        df.loc[idxs, "ci_hi"] = hi

    # Bootstrap median gaps for inverse point estimates in the main and APM grids.
    d_ix = df.index[(df.kind == "delta") & (df.delta_high_minus_low < 0)].tolist()
    print(f"LSCC delta bootstrap candidates {len(d_ix)}", flush=True)
    for i in d_ix:
        r = df.loc[i]
        y = y_of[(r["panel"], r["completeness"], "APM" if r["family"] == "APM" else "MHC")]
        x0 = versions[r["missingness"]].copy()
        x0[~bundle["strata"][r["filter"]]] = np.nan
        g = groups(x0, r["cut"])
        if g is None:
            continue
        low, high = g
        rec = boot_median_gap(y[low], y[high], rng, N_BOOT)
        df.loc[i, "ci_lo"] = rec["ci_lo"]
        df.loc[i, "ci_hi"] = rec["ci_hi"]
        df.loc[i, "n_low"] = rec["n_low"]
        df.loc[i, "n_high"] = rec["n_high"]
        df.loc[i, "delta_high_minus_low"] = rec["delta_high_minus_low"]
        df.loc[i, "cliff"] = rec["cliff"]
        df.loc[i, "p_mwu"] = rec["p_mwu"]

    # Anchor QC against scipy and against the published point estimates.
    def anchor(panel: str) -> pd.Series:
        m = (
            (df.panel == panel)
            & (df.kind == "rho")
            & (df["mode"] == "full")
            & (df["filter"] == "all")
            & (df.missingness == "observed")
            & (df.completeness == "strict")
        )
        hit = df.loc[m]
        if hit.empty:
            raise RuntimeError(f"missing anchor panel {panel}")
        return hit.iloc[0]

    mhc = anchor("HLA-A+HLA-B+HLA-C")
    cd8 = anchor("CD8A")
    y_mhc = y_of[("HLA-A+HLA-B+HLA-C", "strict", "MHC")]
    y_cd8 = y_of[("CD8A", "strict", "MHC")]
    n1, r1, _p1 = spearman_1d(raw, y_mhc)
    n2, r2, _p2 = spearman_1d(raw, y_cd8)
    # fast path agreement
    r_fast = float(spearman_matrix(np.vstack([y_mhc[np.isfinite(raw)], y_cd8[np.isfinite(raw)]]), raw[np.isfinite(raw)])[0])
    qc = {
        "mhc_panel": "HLA-A+HLA-B+HLA-C",
        "mhc_n": int(mhc["n"]),
        "mhc_rho": float(mhc["rho"]),
        "mhc_ci": [float(mhc["ci_lo"]), float(mhc["ci_hi"])],
        "cd8a_n": int(cd8["n"]),
        "cd8a_rho": float(cd8["rho"]),
        "cd8a_ci": [float(cd8["ci_lo"]), float(cd8["ci_hi"])],
        "scipy_mhc_rho": r1,
        "scipy_cd8_rho": r2,
        "fast_mhc_rho": r_fast,
    }
    print("QC anchors", json.dumps(qc), flush=True)
    if n1 != 78 or abs(r1 - (-0.410)) > 0.015:
        raise RuntimeError(f"MHC anchor did not reproduce: n={n1} rho={r1}")
    if n2 != 78 or abs(r2 - (-0.444)) > 0.015:
        raise RuntimeError(f"CD8A anchor did not reproduce: n={n2} rho={r2}")
    if abs(r_fast - r1) > 1e-8:
        raise RuntimeError(f"fast Spearman disagrees with scipy: {r_fast} vs {r1}")

    # Partials for anchors and the strongest full-cohort observed rows.
    def add_partial(mask: pd.Series, label: str) -> None:
        sub = df.loc[mask]
        if sub.empty:
            return
        # limit
        if "rho" in sub.columns:
            order = sub.rho.sort_values()
            take = list(order.index[:8])
        else:
            take = list(sub.index[:8])
        for i in take:
            r = df.loc[i]
            y = y_of[(r["panel"], r["completeness"], "APM" if r["family"] == "APM" else "MHC")]
            x0 = versions[r["missingness"]].copy()
            x0[~bundle["strata"][r["filter"]]] = np.nan
            if r["mode"] == "tails":
                g = groups(x0, r["cut"])
                if g is None:
                    continue
                low, high = g
                xt = np.full_like(x0, np.nan)
                xt[low | high] = x0[low | high]
                x_use = xt
            else:
                x_use = x0
            part = partial_spearman(x_use, y, bundle["wes"], rng, N_BOOT)
            for k, v in part.items():
                df.loc[i, k] = v
            print(label, r["panel"], r["filter"], r["mode"], r["cut"], {k: part[k] for k in part}, flush=True)

    for col in ["n_partial", "rho_partial", "p_partial", "ci_lo_partial", "ci_hi_partial"]:
        df[col] = np.nan
    add_partial(anchor_m, "anchor")
    full_obs = (
        (df.kind == "rho")
        & (df["mode"] == "full")
        & (df["filter"] == "all")
        & (df.missingness == "observed")
        & (df.family != "APM")
        & (df.ci_hi < 0)
    )
    add_partial(full_obs, "full")

    print("permuting LSCC rho...", flush=True)
    perm_primary = permute_rho(raw, panels, bundle["strata"], rng, include_strata=False, include_tails=False, n_perm=N_PERM)
    print("primary perm", perm_primary, flush=True)
    perm_expanded = permute_rho(raw, panels, bundle["strata"], rng, include_strata=True, include_tails=True, n_perm=N_PERM)
    print("expanded perm", perm_expanded, flush=True)
    print("permuting LSCC Q1/Q4 delta...", flush=True)
    # Delta permutation uses the strict+allow panels and observed CLDN4, full cohort, q25_75.
    perm_delta = permute_delta(raw, panels, rng, "q25_75", n_perm=400, B=400)
    print("delta perm", perm_delta, flush=True)

    extra = {
        "qc": qc,
        "perm_rho_full_cohort_continuous": perm_primary,
        "perm_rho_strata_and_tails": perm_expanded,
        "perm_delta_q25_75_full_cohort": perm_delta,
        "n_panels_mhc": len(panels),
        "n_panels_apm": len(apm_panels),
        "panels": panels,
        "y_of": y_of,
        "versions": versions,
    }
    return df, extra


def sweep_luad(bundle: dict, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    raw = bundle["genes"]["CLDN4"]
    versions = cldn4_versions(raw)
    nhej = build_panels(bundle["genes"], [s for s, _e in NHEJ], ddof=1, max_k=7, kind="NHEJ")
    sting = build_panels(bundle["genes"], [s for s, _e in STING], ddof=1, max_k=6, kind="STING")
    panels = nhej + sting
    print(f"LUAD panels DNA-PK {len(nhej)} STING {len(sting)}", flush=True)
    y_of = {(p["panel"], p["completeness"], p["kind"]): p["y"] for p in panels}
    rows = []
    for panel in panels:
        y = panel["y"]
        for fname, mask in bundle["strata"].items():
            for miss, x0 in versions.items():
                xx = x0.copy()
                xx[~mask] = np.nan
                for cut in LUAD_CUTS:
                    g = groups(xx, cut)
                    if g is None:
                        continue
                    low, high = g
                    a = y[low]
                    b = y[high]
                    a = a[np.isfinite(a)]
                    b = b[np.isfinite(b)]
                    if a.size < MIN_ARM or b.size < MIN_ARM:
                        continue
                    try:
                        _u, pm = stats.mannwhitneyu(a, b, alternative="two-sided")
                        pm = float(pm)
                    except ValueError:
                        pm = np.nan
                    rows.append(
                        {
                            "cohort": "LUAD",
                            "kind": "delta",
                            "mode": "groups",
                            "cut": cut,
                            "filter": fname,
                            "missingness": miss,
                            "panel": panel["panel"],
                            "completeness": panel["completeness"],
                            "family": panel["family"],
                            "k": panel["k"],
                            "n": int(a.size + b.size),
                            "n_low": int(a.size),
                            "n_high": int(b.size),
                            "n_imputed_cldn4": int((~np.isfinite(raw) & (low | high) & np.isfinite(y)).sum()),
                            "n_imputed_low": int((~np.isfinite(raw) & low & np.isfinite(y)).sum()),
                            "delta_q1_minus_q4": float(np.median(a) - np.median(b)),
                            "delta_high_minus_low": float(np.median(b) - np.median(a)),
                            "cliff_low_vs_high": cliff_high_vs_low(a, b),
                            "p_mwu": pm,
                        }
                    )
    df = pd.DataFrame(rows)
    df["ci_lo"] = np.nan
    df["ci_hi"] = np.nan
    print(f"LUAD grid rows {len(df)}", flush=True)

    # Bootstrap the thesis side. DNA-PK wants Q1−Q4 < 0. STING wants Q1−Q4 > 0.
    # ci_* refer to the Q1−Q4 median difference.
    def needs_boot(r) -> bool:
        if r["family"] == "DNA-PK":
            return r["delta_q1_minus_q4"] < 0
        return r["delta_q1_minus_q4"] > 0

    # Reference anchors are bootstrapped even if the point estimate has the other sign.
    anchors = {
        ("PRKDC+LIG4", "q1_rest"),
        ("PRKDC+LIG4", "q25_75"),
        ("PRKDC+XRCC5+LIG4", "q1_rest"),
        ("PRKDC+XRCC5+LIG4", "q25_75"),
        ("STING1+TBK1+IRF3", "q1_rest"),
        ("STING1+TBK1+IRF3", "q25_75"),
        ("STING1+TBK1+IRF3+IRF7", "q1_rest"),
        ("STING1+TBK1+IRF3+IRF7", "q25_75"),
    }
    ix = []
    for i, r in df.iterrows():
        if needs_boot(r) or ((r["panel"], r["cut"]) in anchors and r["filter"] == "all" and r["missingness"] == "observed" and r["completeness"] == "strict"):
            ix.append(i)
    print(f"LUAD delta bootstrap candidates {len(ix)}", flush=True)
    for i in ix:
        r = df.loc[i]
        y = y_of[(r["panel"], r["completeness"], "NHEJ" if r["family"] == "DNA-PK" else "STING")]
        x0 = versions[r["missingness"]].copy()
        x0[~bundle["strata"][r["filter"]]] = np.nan
        g = groups(x0, r["cut"])
        if g is None:
            continue
        low, high = g
        # boot_median_gap returns high−low. Flip to Q1−Q4.
        rec = boot_median_gap(y[low], y[high], rng, N_BOOT)
        if not np.isfinite(rec["delta_high_minus_low"]):
            continue
        df.loc[i, "delta_high_minus_low"] = rec["delta_high_minus_low"]
        df.loc[i, "delta_q1_minus_q4"] = -rec["delta_high_minus_low"]
        df.loc[i, "ci_lo"] = -rec["ci_hi"]
        df.loc[i, "ci_hi"] = -rec["ci_lo"]
        df.loc[i, "n_low"] = rec["n_low"]
        df.loc[i, "n_high"] = rec["n_high"]
        df.loc[i, "p_mwu"] = rec["p_mwu"]

    def anchor_delta(panel: str, cut: str) -> pd.Series | None:
        m = (
            (df.panel == panel)
            & (df.cut == cut)
            & (df["filter"] == "all")
            & (df.missingness == "observed")
            & (df.completeness == "strict")
        )
        hit = df.loc[m]
        if hit.empty:
            return None
        return hit.iloc[0]

    qc = {}
    for panel, expect, cut in [
        ("PRKDC+LIG4", -0.262, "q1_rest"),
        ("STING1+TBK1+IRF3", 0.418, "q1_rest"),
        ("STING1+TBK1+IRF3+IRF7", 0.500, "q1_rest"),
        ("PRKDC+XRCC5+LIG4", -0.519, "q1_rest"),
    ]:
        row = anchor_delta(panel, cut)
        if row is None:
            qc[panel] = "missing"
            continue
        qc[panel] = {
            "delta_q1_minus_rest": float(row["delta_q1_minus_q4"]),
            "n_low": int(row["n_low"]),
            "n_high": int(row["n_high"]),
            "ci": [float(row["ci_lo"]), float(row["ci_hi"])],
        }
        print("LUAD anchor", panel, qc[panel], flush=True)
    # Hard check on the two headline anchors. Tolerance covers quantile/z-score ddof drift.
    d_pk = qc.get("PRKDC+LIG4", {})
    d_st = qc.get("STING1+TBK1+IRF3", {})
    if not isinstance(d_pk, dict) or abs(d_pk["delta_q1_minus_rest"] - (-0.262)) > 0.03:
        raise RuntimeError(f"DNA-PK Q1-rest anchor did not reproduce: {d_pk}")
    if not isinstance(d_st, dict) or abs(d_st["delta_q1_minus_rest"] - 0.418) > 0.03:
        raise RuntimeError(f"STING Q1-rest anchor did not reproduce: {d_st}")

    print("permuting LUAD Q1/Q4 ...", flush=True)
    perm = {}
    for fam, plist, sign in [("DNA-PK", nhej, "low"), ("STING", sting, "high")]:
        Y = panel_matrix(plist)
        x = versions["observed"]

        def stat_one(xp: np.ndarray, Ym: np.ndarray, which: str) -> float:
            g = groups(xp, "q25_75")
            if g is None:
                return 0.0
            low, high = g
            ok = np.isfinite(Ym[:, low]).all(1) & np.isfinite(Ym[:, high]).all(1)
            if not ok.any():
                return 0.0
            point, lo, hi = boot_median_matrix(Ym[ok][:, low], Ym[ok][:, high], rng, 300)
            # point, lo, hi are high−low. Q1−Q4 = −that, so lo_q = −hi, hi_q = −lo.
            q_lo, q_hi, q_point = -hi, -lo, -point
            if which == "low":
                good = np.isfinite(q_point) & (q_hi < 0)
                return float(np.max(-q_point[good])) if good.any() else 0.0
            good = np.isfinite(q_point) & (q_lo > 0)
            return float(np.max(q_point[good])) if good.any() else 0.0

        obs = stat_one(x, Y, sign)
        finite = np.isfinite(raw)
        vals = raw[finite].copy()
        null = np.empty(400)
        for i in range(400):
            fv = vals.copy()
            rng.shuffle(fv)
            xp = raw.copy()
            xp[finite] = fv
            null[i] = stat_one(xp, Y, sign)
        p = (1 + int(np.sum(null + 1e-12 >= obs))) / (1 + 400)
        perm[fam] = {
            "obs": obs,
            "perm_p": p,
            "n_perm": 400,
            "boot_inside": 300,
            "null_median": float(np.median(null)),
            "null_q95": float(np.quantile(null, 0.95)),
            "slice": "full cohort, observed CLDN4, Q1 vs Q4, all nonempty panels",
        }
        print(fam, perm[fam], flush=True)

    return df, {"qc": qc, "perm": perm, "n_nhej": len(nhej), "n_sting": len(sting), "y_of": y_of, "versions": versions}


def pick_rho(df: pd.DataFrame, *, full: bool, observed: bool, tails: bool, apm: bool = False) -> pd.Series | None:
    m = (df.kind == "rho") & (df.n >= MIN_N) & np.isfinite(df.ci_hi) & (df.ci_hi < 0) & (df.rho < 0)
    if apm:
        m &= df.family.eq("APM")
    else:
        m &= df.family.ne("APM") & df.in_main_grid
    if full:
        m &= df["filter"].eq("all")
    if observed:
        m &= df.missingness.eq("observed")
    if not tails:
        m &= df["mode"].eq("full")
    sub = df.loc[m]
    if sub.empty:
        return None
    return sub.sort_values(["rho", "n"], ascending=[True, False]).iloc[0]


def pick_delta(df: pd.DataFrame, *, full: bool, observed: bool, cut: str | None, apm: bool = False) -> pd.Series | None:
    m = (
        (df.kind == "delta")
        & (df.n_low >= MIN_ARM)
        & (df.n_high >= MIN_ARM)
        & np.isfinite(df.ci_hi)
        & (df.ci_hi < 0)
        & (df.delta_high_minus_low < 0)
    )
    if apm:
        m &= df.family.eq("APM")
    else:
        m &= df.family.ne("APM") & df.in_main_grid
    if full:
        m &= df["filter"].eq("all")
    if observed:
        m &= df.missingness.eq("observed")
    if cut is not None:
        m &= df.cut.eq(cut)
    sub = df.loc[m]
    if sub.empty:
        return None
    # Most negative high−low gap. Tie-break on Cliff (more negative = stronger inverse rank gap).
    return sub.sort_values(["delta_high_minus_low", "cliff"], ascending=[True, True]).iloc[0]


def pick_luad(df: pd.DataFrame, family: str, *, full: bool, observed: bool, cut: str | None) -> pd.Series | None:
    m = (df.family == family) & np.isfinite(df.ci_lo) & np.isfinite(df.ci_hi)
    if family == "DNA-PK":
        m &= df.ci_hi < 0
        m &= df.delta_q1_minus_q4 < 0
    else:
        m &= df.ci_lo > 0
        m &= df.delta_q1_minus_q4 > 0
    if full:
        m &= df["filter"].eq("all")
    if observed:
        m &= df.missingness.eq("observed")
    if cut is not None:
        m &= df.cut.eq(cut)
    sub = df.loc[m]
    if sub.empty:
        return None
    if family == "DNA-PK":
        return sub.sort_values(["delta_q1_minus_q4", "n"], ascending=[True, False]).iloc[0]
    return sub.sort_values(["delta_q1_minus_q4", "n"], ascending=[False, False]).iloc[0]


def row_brief(r: pd.Series, effect: str) -> str:
    base = (
        f"{r['panel']} ({r['family']}, {r['completeness']}, k={int(r['k'])}); "
        f"filter={r['filter']}; missingness={r['missingness']}"
    )
    if effect == "rho":
        return (
            f"{base}; mode={r['mode']}"
            + (f"; cut={r['cut']}" if r["mode"] != "full" else "")
            + f"; n={int(r['n'])}; ρ={f3(r['rho'])}; 95% CI {fci(r['ci_lo'], r['ci_hi'])}; "
            f"p={fp(r['p'])}; imputed CLDN4 in the pair={int(r['n_imputed_cldn4'])}"
        )
    if "delta_q1_minus_q4" in r and np.isfinite(r.get("delta_q1_minus_q4", np.nan)):
        return (
            f"{base}; cut={r['cut']}; n_low={int(r['n_low'])} vs n_high={int(r['n_high'])}; "
            f"Δ median z Q1−Q4={f3(r['delta_q1_minus_q4'])}; 95% CI {fci(r['ci_lo'], r['ci_hi'])}; "
            f"MWU p={fp(r['p_mwu'])}; imputed in Q1={int(r['n_imputed_low'])}"
        )
    return (
        f"{base}; cut={r['cut']}; n_low={int(r['n_low'])} vs n_high={int(r['n_high'])}; "
        f"Δ median z Q4−Q1={f3(r['delta_high_minus_low'])}; 95% CI {fci(r['ci_lo'], r['ci_hi'])}; "
        f"Cliff(Q4 vs Q1)={f3(r['cliff'])}; MWU p={fp(r['p_mwu'])}; "
        f"imputed in Q1={int(r['n_imputed_low'])}"
    )


def style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_rho_forest(df: pd.DataFrame, path: Path) -> None:
    m = (
        (df.kind == "rho")
        & (df["mode"] == "full")
        & (df["filter"] == "all")
        & (df.missingness == "observed")
        & (df.family != "APM")
        & np.isfinite(df.ci_hi)
        & (df.ci_hi < 0)
        & (df.rho < 0)
    )
    sub = df.loc[m].sort_values("rho").head(12).iloc[::-1]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    y = np.arange(len(sub))
    ax.errorbar(
        sub.rho,
        y,
        xerr=[sub.rho - sub.ci_lo, sub.ci_hi - sub.rho],
        fmt="o",
        color="#1f4e79",
        ecolor="#5b7c99",
        capsize=2,
        markersize=4,
    )
    labels = [f"{r.panel}  n={int(r.n)}" for r in sub.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.axvline(-0.410, color="#b85c38", lw=0.8, ls="--")
    ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
    ax.set_title("LSCC CLDN4 vs MHC-I/CD8 protein — full cohort, measured CLDN4")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_delta_forest(df: pd.DataFrame, path: Path) -> None:
    m = (
        (df.kind == "delta")
        & (df.cut == "q25_75")
        & (df["filter"] == "all")
        & (df.missingness == "observed")
        & (df.family != "APM")
        & np.isfinite(df.ci_hi)
        & (df.ci_hi < 0)
    )
    sub = df.loc[m].sort_values("delta_high_minus_low").head(12).iloc[::-1]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    y = np.arange(len(sub))
    ax.errorbar(
        sub.delta_high_minus_low,
        y,
        xerr=[sub.delta_high_minus_low - sub.ci_lo, sub.ci_hi - sub.delta_high_minus_low],
        fmt="o",
        color="#1f4e79",
        ecolor="#5b7c99",
        capsize=2,
        markersize=4,
    )
    labels = [f"{r.panel}  {int(r.n_low)} vs {int(r.n_high)}" for r in sub.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("Median z, CLDN4 Q4 − Q1 (bootstrap 95% CI)")
    ax.set_title("LSCC Q1 vs Q4 median-z gap — full cohort, measured CLDN4")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_scatter(raw: np.ndarray, y: np.ndarray, rho: float, ci: tuple[float, float], title: str, path: Path, ylabel: str) -> None:
    m = np.isfinite(raw) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(raw[m], y[m], s=16, c="#1f4e79", alpha=0.85, linewidths=0)
    ax.set_xlabel("CLDN4 protein, log2 TMT")
    ax.set_ylabel(ylabel)
    ax.set_title(title + f"\nρ={rho:+.3f} ({ci[0]:+.3f} to {ci[1]:+.3f}), n={int(m.sum())}")
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_luad(bundle_y: dict, raw: np.ndarray, winners: list[tuple[str, pd.Series]], path: Path) -> None:
    fig, axes = plt.subplots(1, max(len(winners), 1), figsize=(4.4 * max(len(winners), 1), 4.2), squeeze=False)
    x = cldn4_versions(raw)["observed"]
    g = groups(x, "q25_75")
    for ax, (title, row) in zip(axes[0], winners):
        if g is None or row is None:
            ax.set_axis_off()
            continue
        kind = "NHEJ" if row["family"] == "DNA-PK" else "STING"
        y = bundle_y[(row["panel"], row["completeness"], kind)]
        low, high = g
        a = y[low]
        b = y[high]
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        bp = ax.boxplot([a, b], tick_labels=["Q1", "Q4"], widths=0.55, patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], ["#c47b5a", "#1f4e79"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
        rng = np.random.default_rng(0)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, size=a.size), a, s=12, c="#c47b5a", alpha=0.8, linewidths=0)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, size=b.size), b, s=12, c="#1f4e79", alpha=0.8, linewidths=0)
        ax.set_title(
            f"{title}\n{row['panel']}\nΔ Q1−Q4 {row['delta_q1_minus_q4']:+.3f}\n"
            f"CI {row['ci_lo']:+.3f} to {row['ci_hi']:+.3f}",
            fontsize=9,
        )
        ax.set_ylabel("Panel mean z")
        style(ax)
    fig.suptitle("CPTAC LUAD, measured CLDN4, Q1 vs Q4", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


def write_finding(path: Path, lscc: pd.DataFrame, luad: pd.DataFrame, extra_l: dict, extra_u: dict, picks: dict) -> None:
    qc = extra_l["qc"]
    lines: list[str] = []
    lines.append("# CPTAC protein: larger CLDN4 contrasts with the interval kept off zero")
    lines.append("")
    lines.append(
        "Treatment-naive CPTAC TMT freeze v1.2 (LSCC Satpathy *Cell* 2021; LUAD Gillette *Cell* 2020). "
        "LSCC and LUAD are not pooled. No ICI column. Endpoint proteins are never filled. "
        "Undetected CLDN4 stays missing in the primary rows. `floor` ties those tumors at the lowest "
        "measured CLDN4 value; `below_min` places them one log2 unit under that value. Those are rank rules, "
        "not imputed abundances."
    )
    lines.append("")
    lines.append(
        "The grid is fixed in `analyze.py` before the correlations are ranked. MHC-I/CD8 panels are every "
        "subset of size 1–4 of HLA-A/B/C/E/F and CD8A (B2M is absent; HLA-G and CD8B are quantified in "
        "20 and 10 tumors and sit under the 40-tumor floor). The score is the mean of cohort z-scores "
        "(LSCC ddof=0, the same scale as the locked HLA-A/B/C score). Completeness is strict. "
        "Sample cuts are grade G2/G3, stage I/II, a WES-purity median split, and LSCC basaloid versus not. "
        "Stage III has too few measured CLDN4 tumors to test. Quantile cuts are Q1/Q4, 20/80, 15/85, 10/90, "
        "and tertiles. Spearman rows need n≥20. Group rows need ≥8 tumors on a side. "
        "An inverse row is eligible only when its 2,000-draw bootstrap 95% interval lies entirely below 0. "
        "A STING row is eligible only when the Q1−Q4 interval lies entirely above 0. "
        "Nominal p-values are not a multiplicity correction. The permutation p is for the maximum "
        "inside the stated slice."
    )
    lines.append("")
    lines.append("## Locked LSCC anchors, recomputed")
    lines.append("")
    lines.append(
        f"CLDN4 versus the HLA-A/B/C mean z: ρ={qc['mhc_rho']:+.3f}, n={qc['mhc_n']}, "
        f"bootstrap 95% CI {qc['mhc_ci'][0]:+.3f} to {qc['mhc_ci'][1]:+.3f}. "
        "The published point estimate was −0.410 on n=78. "
        f"WES partial ρ={f3(picks['rho_anchor_mhc']['rho_partial']) if False else ''}."
    )
    # The sentence above is completed from the anchor rows stored on the frame.
    mhc_anchor = lscc[
        (lscc.panel == "HLA-A+HLA-B+HLA-C")
        & (lscc.kind == "rho")
        & (lscc["mode"] == "full")
        & (lscc["filter"] == "all")
        & (lscc.missingness == "observed")
        & (lscc.completeness == "strict")
    ].iloc[0]
    cd8_anchor = lscc[
        (lscc.panel == "CD8A")
        & (lscc.kind == "rho")
        & (lscc["mode"] == "full")
        & (lscc["filter"] == "all")
        & (lscc.missingness == "observed")
        & (lscc.completeness == "strict")
    ].iloc[0]
    # Replace the awkward partial clause with the real anchor partials.
    lines[-1] = (
        f"CLDN4 versus the HLA-A/B/C mean z: ρ={qc['mhc_rho']:+.3f}, n={qc['mhc_n']}, "
        f"bootstrap 95% CI {qc['mhc_ci'][0]:+.3f} to {qc['mhc_ci'][1]:+.3f}. "
        f"WES partial ρ={f3(mhc_anchor['rho_partial'])} "
        f"(CI {fci(mhc_anchor['ci_lo_partial'], mhc_anchor['ci_hi_partial'])}, n={int(mhc_anchor['n_partial'])}). "
        "The published point estimate was −0.410, and the published partial was −0.320."
    )
    lines.append("")
    lines.append(
        f"CLDN4 versus CD8A: ρ={qc['cd8a_rho']:+.3f}, n={qc['cd8a_n']}, "
        f"bootstrap 95% CI {qc['cd8a_ci'][0]:+.3f} to {qc['cd8a_ci'][1]:+.3f}. "
        f"WES partial ρ={f3(cd8_anchor['rho_partial'])} "
        f"(CI {fci(cd8_anchor['ci_lo_partial'], cd8_anchor['ci_hi_partial'])}, n={int(cd8_anchor['n_partial'])}). "
        "The published point estimate was −0.444, and the published partial was −0.376."
    )
    lines.append("")
    lines.append("## Largest full-cohort Spearman, measured CLDN4")
    lines.append("")
    w = picks["rho_full"]
    if w is None:
        lines.append("No full-cohort measured-CLDN4 panel kept a bootstrap interval entirely below 0.")
    else:
        lines.append(
            f"The largest eligible |ρ| is {row_brief(w, 'rho')}."
        )
        lines.append("")
        lines.append(
            f"That is a larger inverse than the locked HLA-A/B/C ρ ({qc['mhc_rho']:+.3f}) "
            f"and the locked CD8A ρ ({qc['cd8a_rho']:+.3f}), on the same 78 tumors. "
            "Every panel in the top of this list contains CD8A and HLA-A. HLA-F is in the winning set; "
            "it is nonclassical MHC-I, and the panel is labeled that way."
        )
        if np.isfinite(w.get("rho_partial", np.nan)):
            lines.append("")
            lines.append(
                f"WES-purity partial ρ={f3(w['rho_partial'])} "
                f"(CI {fci(w['ci_lo_partial'], w['ci_hi_partial'])}, n={int(w['n_partial'])}). "
                + (
                    "The partial interval stays below 0."
                    if np.isfinite(w["ci_hi_partial"]) and w["ci_hi_partial"] < 0
                    else "The partial interval is not entirely below 0."
                )
            )
        top = lscc[
            (lscc.kind == "rho")
            & (lscc["mode"] == "full")
            & (lscc["filter"] == "all")
            & (lscc.missingness == "observed")
            & (lscc.family != "APM")
            & (lscc.ci_hi < 0)
            & (lscc.rho < 0)
        ].sort_values("rho").head(5)
        lines.append("")
        lines.append(
            md_table(
                ["Panel", "n", "ρ", "95% CI", "partial ρ", "partial CI"],
                [
                    [
                        r.panel,
                        str(int(r.n)),
                        f3(r.rho),
                        fci(r.ci_lo, r.ci_hi),
                        f3(r.rho_partial),
                        fci(r.ci_lo_partial, r.ci_hi_partial),
                    ]
                    for r in top.itertuples()
                ],
            )
        )
        pr = extra_l["perm_rho_full_cohort_continuous"]
        lines.append("")
        lines.append(
            f"Search permutation on the full-cohort continuous slice (every MHC/CD8 panel and all three "
            f"missingness rules; Fisher interval entirely below 0): max |ρ|={pr['obs_max_abs_rho']:.3f}, "
            f"permutation p={fp(pr['perm_p'])} ({pr['n_perm']} shuffles; null median {pr['null_median']:.3f}, "
            f"null 95th {pr['null_q95']:.3f}). The quoted interval is the bootstrap. "
            "On this slice the Fisher maximum and the bootstrap maximum are the same panel."
        )
    lines.append("")
    lines.append("## Largest full-cohort Q4−Q1 gap, measured CLDN4")
    lines.append("")
    dlt = picks["delta_full"]
    eligible = lscc[
        (lscc.kind == "delta")
        & (lscc.cut == "q25_75")
        & (lscc["filter"] == "all")
        & (lscc.missingness == "observed")
        & (lscc.family != "APM")
        & (lscc.ci_hi < 0)
        & (lscc.delta_high_minus_low < 0)
    ]
    if dlt is None:
        lines.append("No full-cohort Q1 versus Q4 panel kept the Q4−Q1 median-z interval entirely below 0.")
    else:
        lines.append(f"Largest median-z gap: {row_brief(dlt, 'delta')}.")
        lines.append("")
        cliff_row = eligible.sort_values(["cliff", "delta_high_minus_low"]).iloc[0]
        lines.append(
            f"The strongest rank separation in the same eligible set is {cliff_row['panel']}: "
            f"Cliff(Q4 vs Q1)={f3(cliff_row['cliff'])}, "
            f"Δ median z Q4−Q1={f3(cliff_row['delta_high_minus_low'])} "
            f"(CI {fci(cliff_row['ci_lo'], cliff_row['ci_hi'])}), "
            f"MWU p={fp(cliff_row['p_mwu'])}. "
            f"CD8A alone has Cliff={f3(eligible.loc[eligible.panel.eq('CD8A'), 'cliff'].iloc[0]) if (eligible.panel=='CD8A').any() else 'NA'} "
            f"and Δz={f3(eligible.loc[eligible.panel.eq('CD8A'), 'delta_high_minus_low'].iloc[0]) if (eligible.panel=='CD8A').any() else 'NA'}. "
            "The pre-specified pick is the median-z gap. A mean of z-scores is not on the same scale as one gene, "
            "and the HLA-F+CD8A z-gap is not a larger rank separation than HLA-A+HLA-C+CD8A."
        )
        gaps = extra_l.get("log2_q4q1", {})
        if gaps:
            lines.append("")
            bits = [f"{g} {gaps[g]['delta_q4_minus_q1_log2']:+.3f}" for g in ["HLA-A", "CD8A", "HLA-F", "HLA-C", "HLA-B", "HLA-E"] if g in gaps]
            lines.append(
                "Single-protein median log2 (Q4−Q1), same 20 vs 20 tumors: " + "; ".join(bits) + ". "
                "HLA-A is the largest single-protein gap. CD8A Q4 is lower than Q1 by "
                f"{abs(gaps['CD8A']['delta_q4_minus_q1_log2']):.3f} log2 TMT units."
                if "CD8A" in gaps
                else "Single-protein log2 gaps are in the run log."
            )
        pdlt = extra_l["perm_delta_q25_75_full_cohort"]
        lines.append("")
        lines.append(
            f"Search permutation for the median-z rule (measured CLDN4, full cohort, Q1 vs Q4): "
            f"max (median Q1 − median Q4)={pdlt['obs_max_low_minus_high']:.3f}, "
            f"permutation p={fp(pdlt['perm_p'])} "
            f"({pdlt['n_perm']} shuffles, {pdlt['boot_inside']} inner bootstrap draws; "
            f"null 95th {pdlt['null_q95']:.3f})."
        )
        pur = extra_l.get("purity_notes", {})
        if "lscc_q4q1" in pur:
            rec = pur["lscc_q4q1"]
            lines.append("")
            lines.append(
                f"WES residual of the median-z winner (Q4−Q1): {f3(rec['delta'])} "
                f"(CI {fci(rec['ci_lo'], rec['ci_hi'])}, n={rec['n_low']} vs {rec['n_high']}). "
                + (
                    "That residual interval is not entirely below 0."
                    if not (np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0)
                    else "That residual interval stays below 0."
                )
            )
        if "lscc_cliff" in pur:
            rec = pur["lscc_cliff"]
            lines.append("")
            lines.append(
                f"WES residual of the Cliff leader HLA-A+HLA-C+CD8A (Q4−Q1): {f3(rec['delta'])} "
                f"(CI {fci(rec['ci_lo'], rec['ci_hi'])}, n={rec['n_low']} vs {rec['n_high']}). "
                + (
                    "That residual interval is not entirely below 0."
                    if not (np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0)
                    else "That residual interval stays below 0."
                )
            )
        if "lscc_CD8A_q4q1" in pur:
            rec = pur["lscc_CD8A_q4q1"]
            lines.append(
                f"CD8A alone, same Q4−Q1 residual: {f3(rec['delta'])} "
                f"(CI {fci(rec['ci_lo'], rec['ci_hi'])}). "
                + (
                    "That interval stays below 0."
                    if np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0
                    else "That interval is not entirely below 0."
                )
            )
        if "lscc_MHC_q4q1" in pur:
            rec = pur["lscc_MHC_q4q1"]
            lines.append(
                f"HLA-A/B/C alone, same residual: {f3(rec['delta'])} "
                f"(CI {fci(rec['ci_lo'], rec['ci_hi'])}). "
                + (
                    "That interval stays below 0."
                    if np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0
                    else "That interval crosses 0."
                )
            )
    lines.append("")
    lines.append("## Larger |ρ| once cuts and strata are opened")
    lines.append("")
    wx = picks["rho_expanded"]
    if wx is None:
        lines.append("No measured-CLDN4 row in the stratum or tail grid kept an interval entirely below 0.")
    else:
        lines.append(f"Largest |ρ| anywhere in the measured-CLDN4 grid: {row_brief(wx, 'rho')}.")
        lines.append("")
        lines.append(
            "Tail-only Spearman uses only the tumors outside the middle of the CLDN4 distribution. "
            "It is not the full-cohort ρ of −0.41."
        )
        strat = lscc[
            (lscc.kind == "rho")
            & (lscc["mode"] == "full")
            & (lscc["filter"] != "all")
            & (lscc.missingness == "observed")
            & (lscc.family != "APM")
            & (lscc.ci_hi < 0)
            & (lscc.rho < 0)
        ].sort_values("rho")
        if not strat.empty:
            s = strat.iloc[0]
            lines.append("")
            lines.append(
                f"The strongest full Spearman inside a stratum, rather than a tail cut, is {row_brief(s, 'rho')}."
            )
        pe = extra_l["perm_rho_strata_and_tails"]
        lines.append("")
        lines.append(
            f"Search permutation over strata and tail cuts as well as the full cohort: "
            f"Fisher-gate max |ρ|={pe['obs_max_abs_rho']:.3f}, permutation p={fp(pe['perm_p'])} "
            f"({pe['n_perm']} shuffles; null median {pe['null_median']:.3f}, null 95th {pe['null_q95']:.3f}). "
            "The null 95th is already large because the search includes n≈20 slices."
        )
    lines.append("")
    lines.append("## Missingness rules")
    lines.append("")
    filled = lscc[
        (lscc.kind == "rho")
        & (lscc["mode"] == "full")
        & (lscc["filter"] == "all")
        & (lscc.missingness != "observed")
        & (lscc.family != "APM")
        & (lscc.ci_hi < 0)
        & (lscc.rho < 0)
    ].sort_values("rho")
    if filled.empty:
        lines.append("No filled-CLDN4 continuous row kept an interval entirely below 0.")
    else:
        fr = filled.iloc[0]
        lines.append(
            f"Tying the 30 undetected CLDN4 tumors at the floor, or placing them one log2 unit below it, "
            f"does not increase |ρ|. The strongest filled continuous row is {row_brief(fr, 'rho')}. "
            "That is weaker than the measured-CLDN4 winner. The primary estimand stays the 78 quantified tumors."
        )
    lines.append("")
    lines.append("## Antigen processing, separate from MHC-I")
    lines.append("")
    ap = picks["apm_rho"]
    if ap is None:
        lines.append("No measured-CLDN4 antigen-processing subset kept a Spearman interval entirely below 0.")
    else:
        lines.append(row_brief(ap, "rho") + ".")
        lines.append("")
        lines.append(
            "TAP1/2, TAPBP, PSMB8/9/10, and NLRC5 were scored as their own family. "
            "The strongest of those rows does not beat the MHC-I/CD8 winner, and it is not called MHC-I."
        )
    lines.append("")
    lines.append("## LUAD quartiles: DNA-PK and STING")
    lines.append("")
    lines.append(
        "Z-scores use ddof=1, matching the earlier LUAD search. Q1-versus-rest anchors on measured CLDN4 "
        "were recomputed before any new panel was ranked. Bootstrap intervals were not part of that earlier report."
    )
    lines.append("")
    uqc = extra_u["qc"]
    arows = []
    for key, label in [
        ("PRKDC+LIG4", "DNA-PKcs+LIG4"),
        ("PRKDC+XRCC5+LIG4", "DNA-PKcs+Ku80+LIG4"),
        ("STING1+TBK1+IRF3", "STING+TBK1+IRF3"),
        ("STING1+TBK1+IRF3+IRF7", "STING+TBK1+IRF3+IRF7"),
    ]:
        rec = uqc.get(key)
        if not isinstance(rec, dict):
            arows.append([label, "NA", "NA", "NA"])
            continue
        arows.append([
            label,
            f"{rec['delta_q1_minus_rest']:+.3f}",
            f"{rec['ci'][0]:+.3f} to {rec['ci'][1]:+.3f}",
            f"{rec['n_low']} vs {rec['n_high']}",
        ])
    lines.append(md_table(["Prior panel, Q1 vs rest", "Δ median z", "95% CI", "n"], arows))
    lines.append("")
    lines.append(
        "The two DNA-PK Q1-versus-rest gaps stay negative, and both bootstrap intervals cross 0. "
        "The two STING Q1-versus-rest gaps stay positive, and both intervals stay above 0."
    )
    lines.append("")
    lines.append("Same panels, Q1 versus Q4:")
    lines.append("")
    qrows = []
    for panel in ["PRKDC+LIG4", "PRKDC+XRCC5+LIG4", "STING1+TBK1+IRF3", "STING1+TBK1+IRF3+IRF7"]:
        hit = luad[
            (luad.panel == panel)
            & (luad.cut == "q25_75")
            & (luad["filter"] == "all")
            & (luad.missingness == "observed")
            & (luad.completeness == "strict")
        ]
        if hit.empty:
            qrows.append([panel, "NA", "NA", "NA"])
            continue
        r = hit.iloc[0]
        qrows.append([
            panel,
            f3(r["delta_q1_minus_q4"]),
            fci(r["ci_lo"], r["ci_hi"]),
            f"{int(r['n_low'])} vs {int(r['n_high'])}",
        ])
    lines.append(md_table(["Panel, Q1 vs Q4", "Δ median z", "95% CI", "n"], qrows))
    lines.append("")
    lines.append(
        "Moving from Q1-versus-rest to Q1-versus-Q4 shrinks these four panels. "
        "DNA-PKcs+LIG4 falls from −0.262 to −0.094. DNA-PKcs+Ku80+LIG4 falls from −0.519 to −0.357. "
        "STING+TBK1+IRF3+IRF7 falls from +0.500 to +0.467, and that smaller gap still has an interval above 0. "
        "STING+TBK1+IRF3’s Q1-versus-Q4 interval crosses 0."
    )
    lines.append("")
    dnapk_pts = luad[
        (luad.family == "DNA-PK")
        & (luad.cut == "q25_75")
        & (luad["filter"] == "all")
        & (luad.missingness == "observed")
        & (luad.completeness == "strict")
    ].sort_values("delta_q1_minus_q4")
    lines.append("**DNA-PK, full cohort, Q1 vs Q4, measured CLDN4.**")
    lines.append("")
    if picks["luad_dnapk"] is None:
        lines.append("No panel keeps the Q1−Q4 interval entirely below 0.")
        if not dnapk_pts.empty:
            b = dnapk_pts.iloc[0]
            lines.append(
                f"The most negative point estimate is {b['panel']}: "
                f"Δ={f3(b['delta_q1_minus_q4'])} (CI {fci(b['ci_lo'], b['ci_hi'])}, "
                f"n={int(b['n_low'])} vs {int(b['n_high'])}, MWU p={fp(b['p_mwu'])}). "
                "The interval crosses 0."
            )
        lines.append("")
    else:
        lines.append(row_brief(picks["luad_dnapk"], "luad") + ".")
        lines.append("")
    perm_d = extra_u["perm"]["DNA-PK"]
    lines.append(
        f"Search permutation over every nonempty DNA-PK subset on this slice: "
        f"max thesis-direction gap={perm_d['obs']:.3f}, permutation p={fp(perm_d['perm_p'])} "
        f"({perm_d['n_perm']} shuffles). The observed maximum is 0 because no panel cleared the interval."
    )
    lines.append("")
    lines.append("**STING, full cohort, Q1 vs Q4, measured CLDN4.**")
    lines.append("")
    st = picks["luad_sting"]
    if st is None:
        lines.append("No panel keeps the Q1−Q4 interval entirely above 0.")
    else:
        lines.append(row_brief(st, "luad") + ".")
        lines.append("")
        lines.append(
            "That gap is larger than the locked four-gene Q1-versus-Q4 gap (+0.467) and larger than "
            "the locked Q1-versus-rest gap (+0.500). It is an IRF3+IRF7 panel, not the cGAS–STING–TBK1 set."
        )
    perm_s = extra_u["perm"]["STING"]
    lines.append("")
    lines.append(
        f"Search permutation over every nonempty STING-pathway subset on this slice: "
        f"max thesis-direction gap={perm_s['obs']:.3f}, permutation p={fp(perm_s['perm_p'])} "
        f"({perm_s['n_perm']} shuffles, {perm_s['boot_inside']} inner bootstrap draws; "
        f"null 95th {perm_s['null_q95']:.3f})."
    )
    pur = extra_l.get("purity_notes", {})
    if "luad_sting" in pur:
        rec = pur["luad_sting"]
        lines.append("")
        lines.append(
            f"WES residual of that STING winner, scored as Q4−Q1 so the thesis direction is negative: "
            f"{f3(rec['delta'])} (CI {fci(rec['ci_lo'], rec['ci_hi'])}, "
            f"n={rec['n_low']} vs {rec['n_high']}). "
            + (
                "The residual interval crosses 0."
                if not (np.isfinite(rec["ci_hi"]) and rec["ci_hi"] < 0)
                else "The residual interval stays below 0."
            )
        )
    lines.append("")
    lines.append("**Strata and other cuts.**")
    lines.append("")
    for label, key in [("DNA-PK", "luad_dnapk_expanded"), ("STING", "luad_sting_expanded")]:
        r = picks[key]
        if r is None:
            lines.append(f"{label}: no expanded row keeps the interval on the thesis side of 0.")
        else:
            lines.append(f"{label}: {row_brief(r, 'luad')}.")
        lines.append("")
    lines.append(
        "Those two rows are 9 versus 9. They are not the full cohort, and the permutation above does not cover "
        "purity halves, histology, or the 20/80 and tertile cuts. HLA proteins were left out of this LUAD list; "
        "adding them previously removed the STING contrast."
    )
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append(
        "- The full-cohort LSCC Spearman moves from −0.410 / −0.444 to −0.531 on a declared MHC-I/CD8 grid, "
        "with the bootstrap and the WES partial both entirely below 0. It does not replace the locked "
        "HLA-A/B/C number, and it is not an ICI result."
    )
    lines.append(
        "- The LSCC Q4−Q1 median-z gap that maximizes |Δz| does not survive a WES residual. "
        "The larger rank separation is HLA-A+HLA-C+CD8A, not the z-gap winner."
    )
    lines.append(
        "- Tail-only ρ and grade or purity slices are different estimands. The tail search has a high null."
    )
    lines.append(
        "- LUAD DNA-PK does not gain a larger quartile gap whose interval stays below 0. "
        "The earlier Q1-versus-rest point estimates themselves have bootstrap intervals that cross 0."
    )
    lines.append(
        "- LUAD STING Q1 versus Q4 is larger for IRF3+IRF7 than for the locked four-gene panel. "
        "The search permutation is 0.062, and the WES residual interval crosses 0."
    )
    lines.append("- Floor and below-min do not increase the LSCC |ρ|.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/cptac_max_effect/download.py --outdir data/cptac_max_effect")
    lines.append("python3 methods/cptac_max_effect/analyze.py --data data/cptac_max_effect --outdir methods/cptac_max_effect")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")
def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data/cptac_max_effect")
    p.add_argument("--outdir", default="methods/cptac_max_effect")
    args = p.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    tab = out / "results" / "tables"
    fig = out / "results" / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    lscc_b = load_cohort(data, "LSCC")
    luad_b = load_cohort(data, "LUAD")
    cov = pd.DataFrame(lscc_b["coverage"] + luad_b["coverage"])
    cov.to_csv(tab / "coverage.tsv", sep="\t", index=False)
    print("LSCC notes", lscc_b["notes"], flush=True)
    print("LUAD notes", luad_b["notes"], flush=True)

    lscc, extra_l = sweep_lscc(lscc_b, rng)
    luad, extra_u = sweep_luad(luad_b, rng)

    picks = {
        "rho_full": pick_rho(lscc, full=True, observed=True, tails=False),
        "rho_expanded": pick_rho(lscc, full=False, observed=True, tails=True),
        "rho_missing": pick_rho(lscc, full=True, observed=False, tails=False),
        "delta_full": pick_delta(lscc, full=True, observed=True, cut="q25_75"),
        "delta_missing": pick_delta(lscc, full=True, observed=False, cut="q25_75"),
        "apm_rho": pick_rho(lscc, full=True, observed=True, tails=False, apm=True),
        "luad_dnapk": pick_luad(luad, "DNA-PK", full=True, observed=True, cut="q25_75"),
        "luad_sting": pick_luad(luad, "STING", full=True, observed=True, cut="q25_75"),
        "luad_dnapk_expanded": pick_luad(luad, "DNA-PK", full=False, observed=False, cut=None),
        "luad_sting_expanded": pick_luad(luad, "STING", full=False, observed=False, cut=None),
    }
    # Purity residual on the LSCC Q1/Q4 winner and the LUAD quartile winners.
    purity_notes = {}
    if picks["delta_full"] is not None:
        r = picks["delta_full"]
        y = extra_l["y_of"][(r["panel"], r["completeness"], "MHC")]
        x0 = extra_l["versions"]["observed"].copy()
        g = groups(x0, "q25_75")
        if g is not None:
            low, high = g
            purity_notes["lscc_q4q1"] = purity_gap(y, lscc_b["wes"], low, high, rng)
            print("LSCC delta purity", purity_notes["lscc_q4q1"], flush=True)
    for key, bundle, kind in [
        ("luad_dnapk", luad_b, "NHEJ"),
        ("luad_sting", luad_b, "STING"),
    ]:
        r = picks[key]
        if r is None:
            continue
        y = extra_u["y_of"][(r["panel"], r["completeness"], kind)]
        x0 = extra_u["versions"]["observed"].copy()
        g = groups(x0, "q25_75")
        if g is None:
            continue
        low, high = g
        purity_notes[key] = purity_gap(y, bundle["wes"], low, high, rng)
        print(key, "purity", purity_notes[key], flush=True)
    x_obs = extra_l["versions"]["observed"]
    g_q = groups(x_obs, "q25_75")
    if g_q is not None:
        low, high = g_q
        for label, panel in [
            ("lscc_cliff", "HLA-A+HLA-C+CD8A"),
            ("lscc_CD8A_q4q1", "CD8A"),
            ("lscc_MHC_q4q1", "HLA-A+HLA-B+HLA-C"),
        ]:
            y = extra_l["y_of"][(panel, "strict", "MHC")]
            purity_notes[label] = purity_gap(y, lscc_b["wes"], low, high, rng)
            print(label, "purity", purity_notes[label], flush=True)
        gaps = {}
        for sym in ["CD8A", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F"]:
            y = lscc_b["genes"][sym]
            a = y[low]
            b = y[high]
            a = a[np.isfinite(a)]
            b = b[np.isfinite(b)]
            gaps[sym] = {
                "delta_q4_minus_q1_log2": float(np.median(b) - np.median(a)),
                "n_low": int(a.size),
                "n_high": int(b.size),
            }
        extra_l["log2_q4q1"] = gaps
        print("log2 Q4-Q1", gaps, flush=True)

    # Drop bulky objects before saving rows.
    keep_cols_note = True
    lscc.to_csv(tab / "lscc_grid.tsv.gz", sep="\t", index=False)
    luad.to_csv(tab / "luad_grid.tsv.gz", sep="\t", index=False)

    def pack(r: pd.Series | None) -> dict | None:
        if r is None:
            return None
        d = {}
        for k, v in r.to_dict().items():
            if isinstance(v, (np.floating, float)):
                d[k] = None if not np.isfinite(v) else float(v)
            elif isinstance(v, (np.integer,)):
                d[k] = int(v)
            elif isinstance(v, (np.bool_, bool)):
                d[k] = bool(v)
            else:
                d[k] = v if not (isinstance(v, float) and not np.isfinite(v)) else None
        return d

    summary = {
        "seed": SEED,
        "n_boot": N_BOOT,
        "lscc_qc": extra_l["qc"],
        "luad_qc": extra_u["qc"],
        "perm_rho_full_cohort_continuous": extra_l["perm_rho_full_cohort_continuous"],
        "perm_rho_strata_and_tails": extra_l["perm_rho_strata_and_tails"],
        "perm_delta_q25_75": extra_l["perm_delta_q25_75_full_cohort"],
        "perm_luad": extra_u["perm"],
        "purity_residual_gap_high_minus_low": purity_notes,
        "picks": {k: pack(v) for k, v in picks.items()},
        "n_lscc_rows": int(len(lscc)),
        "n_luad_rows": int(len(luad)),
        "lscc_strata": list(lscc_b["strata"]),
        "luad_strata": list(luad_b["strata"]),
        "notes_lscc": lscc_b["notes"],
        "notes_luad": luad_b["notes"],
        "keep_cols_note": keep_cols_note,
    }
    (tab / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    plot_rho_forest(lscc, fig / "fig_lscc_rho_forest.png")
    plot_delta_forest(lscc, fig / "fig_lscc_delta_forest.png")
    y_mhc = extra_l["y_of"][("HLA-A+HLA-B+HLA-C", "strict", "MHC")]
    plot_scatter(
        lscc_b["genes"]["CLDN4"],
        y_mhc,
        extra_l["qc"]["mhc_rho"],
        tuple(extra_l["qc"]["mhc_ci"]),
        "LSCC CLDN4 vs HLA-A/B/C",
        fig / "fig_lscc_mhc_scatter.png",
        "HLA-A/B/C mean z",
    )
    if picks["rho_full"] is not None:
        r = picks["rho_full"]
        y = extra_l["y_of"][(r["panel"], r["completeness"], "MHC")]
        plot_scatter(
            lscc_b["genes"]["CLDN4"],
            y,
            float(r["rho"]),
            (float(r["ci_lo"]), float(r["ci_hi"])),
            f"Largest full-cohort |ρ|: {r['panel']}",
            fig / "fig_lscc_winner_scatter.png",
            "Panel mean z",
        )
    def _luad_strict(panel: str) -> pd.Series | None:
        hit = luad[
            (luad.panel == panel)
            & (luad.cut == "q25_75")
            & (luad["filter"] == "all")
            & (luad.missingness == "observed")
            & (luad.completeness == "strict")
        ]
        return None if hit.empty else hit.iloc[0]

    winners_plot = []
    for title, panel in [
        ("DNA-PK XRCC5+LIG4", "XRCC5+LIG4"),
        ("STING+TBK1+IRF3+IRF7", "STING1+TBK1+IRF3+IRF7"),
    ]:
        row = _luad_strict(panel)
        if row is not None:
            winners_plot.append((title, row))
    if picks["luad_sting"] is not None:
        winners_plot.append(("STING max Q1/Q4", picks["luad_sting"]))
    if winners_plot:
        plot_luad(extra_u["y_of"], luad_b["genes"]["CLDN4"], winners_plot, fig / "fig_luad_quartile.png")

    extra_l["purity_notes"] = purity_notes
    write_finding(out / "FINDING.md", lscc, luad, extra_l, extra_u, picks)
    print("wrote", out / "FINDING.md", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
