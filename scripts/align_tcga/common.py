"""
Shared helpers for the USER-ALIGN TACSTD2 / immune-signature replication.

All statistics are computed from public open-access matrices only.
No values are hard-coded or fabricated; every number written to results/
is produced by these functions from the downloaded data.
"""
from __future__ import annotations

import gzip
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Gene signatures (immune / cytotoxic / exhaustion).
# Literature-based, standard sets. Symbols follow the TCGA HiSeqV2 (hg19/RSEM)
# annotation; VSIR is stored under its legacy symbol C10orf54.
# ---------------------------------------------------------------------------
SIGNATURES: dict[str, list[str]] = {
    # Cytolytic activity, Rooney et al. Cell 2015 (geometric mean of GZMA/PRF1).
    "Cytolytic_CYT": ["GZMA", "PRF1"],
    # CD8 T-cell markers.
    "CD8_Tcell": ["CD8A", "CD8B"],
    # Broad cytotoxic effector program.
    "Cytotoxic_effector": ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY",
                            "NKG7", "KLRD1", "KLRK1"],
    # T-cell exhaustion / co-inhibitory checkpoints.
    "Exhaustion": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA",
                   "VSIR"],
    # IFN-gamma 6-gene signature, Ayers et al. JCI 2017.
    "IFNgamma_6gene": ["IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"],
    # T-cell-inflamed GEP (18-gene), Ayers et al. JCI 2017.
    "Tcell_inflamed_GEP": ["CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1",
                           "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E",
                           "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10",
                           "STAT1", "TIGIT"],
}

# Individual immune / exhaustion genes reported one-by-one.
SINGLE_IMMUNE_GENES = ["CD8A", "GZMA", "GZMB", "PRF1", "IFNG",
                       "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
                       "CD274", "GZMK", "CXCL9"]

# Tight-junction "intersection" genes referenced by the user.
TJ_GENES = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]

# Legacy-symbol remaps (query symbol -> symbol in the matrix).
SYMBOL_ALIASES = {"VSIR": "C10orf54"}

TARGET_GENE = "TACSTD2"


def load_xena_matrix(path: str) -> pd.DataFrame:
    """Load a Xena HiSeqV2 gz matrix: genes x samples, values log2(norm+1)."""
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def resolve_symbol(sym: str, available: set[str]) -> str | None:
    if sym in available:
        return sym
    alias = SYMBOL_ALIASES.get(sym)
    if alias and alias in available:
        return alias
    return None


def signature_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean of per-gene z-scored log2 expression across samples (columns)."""
    avail = set(expr.index)
    rows = []
    used = []
    for g in genes:
        r = resolve_symbol(g, avail)
        if r is None:
            continue
        used.append(g)
        v = expr.loc[r].astype(float)
        z = (v - v.mean()) / v.std(ddof=0)
        rows.append(z)
    if not rows:
        return pd.Series(dtype=float), used
    mat = pd.concat(rows, axis=1)
    return mat.mean(axis=1), used


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return np.nan, np.nan, len(x)
    r, p = stats.spearmanr(x, y)
    return float(r), float(p), int(len(x))


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR. NaNs stay NaN."""
    q = [np.nan] * len(pvals)
    finite_idx = [i for i, p in enumerate(pvals) if np.isfinite(p)]
    if not finite_idx:
        return q
    finite_idx.sort(key=lambda i: pvals[i])
    m = len(finite_idx)
    prev = 1.0
    for rank, i in enumerate(reversed(finite_idx), start=0):
        k = m - rank  # 1-based rank of this p among finite values
        prev = min(prev, pvals[i] * m / k)
        q[i] = min(prev, 1.0)
    return q


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray
                     ) -> tuple[float, float, int]:
    """First-order partial Spearman correlation of x,y controlling for z.

    Computed as the partial Pearson correlation on the rank-transformed
    variables. p-value from a t-test with df = n - 3.
    """
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]
    n = len(x)
    if n < 5:
        return np.nan, np.nan, n
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    denom = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    if denom == 0:
        return np.nan, np.nan, n
    pr = (rxy - rxz * ryz) / denom
    pr = max(min(pr, 1.0), -1.0)
    df = n - 3
    if df <= 0 or abs(pr) >= 1:
        return float(pr), np.nan, n
    t = pr * np.sqrt(df / (1 - pr**2))
    p = 2 * stats.t.sf(abs(t), df)
    return float(pr), float(p), n
