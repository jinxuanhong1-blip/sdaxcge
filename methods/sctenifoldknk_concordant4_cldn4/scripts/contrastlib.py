"""Shared contrast helpers. TMM and family OLS follow the locked concordant-4 script."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
LOCKED_UNITS = ROOT / "data" / "locked_patient_units.tsv"
DATASETS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_LEVELS = DATASETS


def load_pseudobulk(work: Path) -> dict[str, pd.DataFrame]:
    frames = {}
    for ds in DATASETS:
        path = work / ds / "pseudobulk.tsv.gz"
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index = df.index.astype(str).str.upper()
        if df.index.has_duplicates:
            df = df.groupby(level=0).sum()
        df = df.apply(pd.to_numeric, errors="raise")
        frames[ds] = df
    return frames


def load_units(work: Path) -> pd.DataFrame:
    parts = [pd.read_csv(work / ds / "units.tsv", sep="\t") for ds in DATASETS]
    return pd.concat(parts, ignore_index=True)


def align_zero_fill(frames: dict[str, pd.DataFrame]):
    genes = sorted(set().union(*[set(df.index) for df in frames.values()]))
    columns = []
    blocks = []
    panel = {ds: set(df.index) for ds, df in frames.items()}
    for ds in DATASETS:
        df = frames[ds]
        for unit in df.columns:
            columns.append({"dataset": ds, "unit_id": str(unit), "col": f"{ds}|{unit}"})
            blocks.append(df[unit].reindex(genes).fillna(0.0).to_numpy(dtype=np.float64))
    mat = np.column_stack(blocks)
    meta = pd.DataFrame(columns)
    return genes, mat, meta, panel


def tmm_factors(counts: np.ndarray) -> np.ndarray:
    """Port of the locked tmm_factors(). counts is genes x samples."""
    counts = np.asarray(counts, dtype=np.float64)
    lib = counts.sum(axis=0)
    lib = np.where(lib == 0, 1.0, lib)
    rel = counts / lib
    q75 = np.empty(counts.shape[1], dtype=np.float64)
    for j in range(counts.shape[1]):
        x = rel[:, j]
        x = x[x > 0]
        q75[j] = 0.0 if x.size == 0 else float(np.quantile(x, 0.75))
    ref = int(np.argmin(np.abs(q75 - q75.mean())))
    sf = np.ones(counts.shape[1], dtype=np.float64)
    for j in range(counts.shape[1]):
        if j == ref:
            continue
        keep = (counts[:, j] > 0) & (counts[:, ref] > 0)
        if int(keep.sum()) < 20:
            sf[j] = 1.0
            continue
        cj = counts[keep, j]
        cr = counts[keep, ref]
        m = np.log2((cj / lib[j]) / (cr / lib[ref]))
        a = 0.5 * np.log2((cj / lib[j]) * (cr / lib[ref]))
        n = m.size

        def trim_idx(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
            order = np.argsort(values, kind="mergesort")
            start = int(np.ceil(lo * n))
            end = max(start, int(np.floor(hi * n)))
            return order[start - 1 : end]

        keep_m = trim_idx(m, 0.3, 0.7)
        keep_a = trim_idx(a, 0.05, 0.95)
        keep2 = np.intersect1d(keep_m, keep_a)
        if keep2.size < 10:
            keep2 = np.arange(n)
        weights = 1.0 / (1.0 / np.maximum(cj[keep2], 1.0) + 1.0 / np.maximum(cr[keep2], 1.0))
        sf[j] = float(2.0 ** np.average(m[keep2], weights=weights))
    return sf / sf.mean()


def log2_tmm_cpm(counts: np.ndarray) -> np.ndarray:
    sf = tmm_factors(counts)
    lib = counts.sum(axis=0)
    norm = (lib / sf) / 1e6
    norm = np.where(norm == 0, 1.0, norm)
    cpm = counts / norm
    return np.log2(cpm + 1.0), sf


def ols_q4(y: np.ndarray, dataset: np.ndarray, is_q4: np.ndarray):
    """lm(y ~ cohort + Q4). Coefficient is Q4 minus Q1, cohort-adjusted."""
    y = np.asarray(y, dtype=np.float64)
    dataset = np.asarray(dataset)
    is_q4 = np.asarray(is_q4, dtype=np.float64)
    dummies = []
    for level in COHORT_LEVELS[1:]:
        dummies.append((dataset == level).astype(np.float64))
    x = np.column_stack([np.ones(y.shape[0]), *dummies, is_q4]) if dummies else np.column_stack([np.ones(y.shape[0]), is_q4])
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    df = int(y.shape[0] - rank)
    if df <= 0:
        return float(beta[-1]), np.nan, df, int(rank)
    sigma2 = float(np.sum(resid ** 2) / df)
    xtx = x.T @ x
    try:
        cov = sigma2 * np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        cov = sigma2 * np.linalg.pinv(xtx)
    se = float(np.sqrt(max(cov[-1, -1], 0.0)))
    if se == 0 or not np.isfinite(se):
        return float(beta[-1]), np.nan, df, int(rank)
    tstat = float(beta[-1] / se)
    p = float(2 * stats.t.sf(abs(tstat), df))
    return float(beta[-1]), p, df, int(rank)


def ols_matrix(logcpm: np.ndarray, dataset: np.ndarray, is_q4: np.ndarray):
    """Vectorized lm for genes x samples. Returns logFC, p per gene."""
    y = np.asarray(logcpm, dtype=np.float64)
    dataset = np.asarray(dataset)
    is_q4 = np.asarray(is_q4, dtype=np.float64)
    dummies = [(dataset == level).astype(np.float64) for level in COHORT_LEVELS[1:]]
    x = np.column_stack([np.ones(y.shape[1]), *dummies, is_q4])
    beta, _, rank, _ = np.linalg.lstsq(x, y.T, rcond=None)
    fitted = (x @ beta).T
    resid = y - fitted
    df = int(y.shape[1] - rank)
    sigma2 = np.sum(resid ** 2, axis=1) / df
    xtx_inv = np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[-1, -1], 0.0))
    logfc = beta[-1]
    tstat = np.divide(logfc, se, out=np.full_like(logfc, np.nan), where=se > 0)
    p = 2 * stats.t.sf(np.abs(tstat), df)
    return logfc, p, df
