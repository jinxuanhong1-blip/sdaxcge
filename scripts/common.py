"""Shared helpers for the Claim A10 analysis.

Claim A10 (see claims/claim_A10.md):
    TF module ELF3/GRHL1/KLF4/TFAP2A  <->  target program TACSTD2/CLDN4,
    in a state where NKX2-1 is down.  Human + mouse lung.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
INTERIM = os.path.join(ROOT, "data", "interim")
OUT = os.path.join(ROOT, "results", "claim_A10")
TABLES = os.path.join(OUT, "tables")
FIGURES = os.path.join(OUT, "figures")
LOGS = os.path.join(OUT, "logs")

for _d in (RAW, INTERIM, TABLES, FIGURES, LOGS):
    os.makedirs(_d, exist_ok=True)

RNG_SEED = 20240610

MODULE_H = ["ELF3", "GRHL1", "KLF4", "TFAP2A"]
TARGET_H = ["TACSTD2", "CLDN4"]
NKX_H = "NKX2-1"
CLAIM_GENES_H = MODULE_H + TARGET_H + [NKX_H]

MODULE_M = ["Elf3", "Grhl1", "Klf4", "Tfap2a"]
TARGET_M = ["Tacstd2", "Cldn4"]
NKX_M = "Nkx2-1"
CLAIM_GENES_M = MODULE_M + TARGET_M + [NKX_M]

# Lineage / composition markers used to test whether bulk associations are just
# "how much epithelium (and which kind) is in this piece of lung".
MARKERS_H = {
    "EPCAM": "pan-epithelial",
    "SFTPC": "AT2",
    "SFTPB": "AT2/secretory",
    "AGER": "AT1",
    "SCGB1A1": "club/secretory airway",
    "SCGB3A2": "club/secretory airway",
    "KRT5": "basal airway",
    "TP63": "basal airway",
    "FOXJ1": "ciliated airway",
    "MUC5B": "mucous airway",
    "PTPRC": "immune",
    "COL1A1": "fibroblast",
    "PECAM1": "endothelium",
}
MARKERS_M = {
    "Epcam": "pan-epithelial",
    "Sftpc": "AT2",
    "Sftpb": "AT2/secretory",
    "Ager": "AT1",
    "Scgb1a1": "club/secretory airway",
    "Scgb3a2": "club/secretory airway",
    "Krt5": "basal airway",
    "Trp63": "basal airway",
    "Foxj1": "ciliated airway",
    "Muc5b": "mucous airway",
    "Ptprc": "immune",
    "Col1a1": "fibroblast",
    "Pecam1": "endothelium",
}


# ----------------------------------------------------------------------------
# IO
# ----------------------------------------------------------------------------
def read_recount3_annotation(path: str) -> pd.DataFrame:
    """gene_id -> gene_name / gene_type / length from a recount3 gene_sums GTF."""
    rows = []
    pat_id = re.compile(r'gene_id "([^"]+)"')
    pat_name = re.compile(r'gene_name "([^"]+)"')
    pat_type = re.compile(r'gene_type "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = f[8]
            gid = pat_id.search(attr)
            if gid is None:
                continue
            name = pat_name.search(attr)
            gtype = pat_type.search(attr)
            try:
                length = float(f[5])
            except ValueError:
                length = np.nan
            rows.append((gid.group(1), name.group(1) if name else gid.group(1),
                         gtype.group(1) if gtype else "NA", length))
    return pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_type", "length"])


def read_recount3_counts(path: str) -> pd.DataFrame:
    """recount3 gene_sums file -> genes x samples DataFrame of coverage counts."""
    return pd.read_csv(path, sep="\t", comment="#", index_col=0)


def read_recount3_metadata(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


# ----------------------------------------------------------------------------
# Normalisation
# ----------------------------------------------------------------------------
def to_log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Library-size normalise then log2(x+1).

    recount3 'gene sums' are base-level coverage sums, so absolute values are
    length-scaled. That is irrelevant here: every downstream statistic is a
    rank correlation or a within-gene z-score across samples, both of which are
    invariant to a per-gene constant factor.
    """
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def collapse_to_symbols(mat: pd.DataFrame, annot: pd.DataFrame) -> pd.DataFrame:
    """Map gene_id rows to gene symbols, keeping the highest-expressed duplicate."""
    id2name = annot.set_index("gene_id")["gene_name"]
    sym = mat.index.map(id2name)
    keep = pd.notna(sym)
    mat = mat.loc[keep].copy()
    mat.index = pd.Index(sym[keep], name="gene_name")
    order = mat.mean(axis=1).sort_values(ascending=False).index
    mat = mat.loc[order]
    return mat[~mat.index.duplicated(keep="first")]


def expressed_genes(logcpm: pd.DataFrame, min_frac: float = 0.25,
                    min_value: float = 1.0) -> pd.Index:
    """Genes detected above `min_value` log2CPM in at least `min_frac` of samples."""
    frac = (logcpm > min_value).mean(axis=1)
    return logcpm.index[frac >= min_frac]


# ----------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------
def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return q


def spearman_pairs(mat: pd.DataFrame, genes_a, genes_b=None) -> pd.DataFrame:
    """All pairwise Spearman correlations between two gene lists (rows of mat)."""
    genes_b = genes_a if genes_b is None else genes_b
    have_a = [g for g in genes_a if g in mat.index]
    have_b = [g for g in genes_b if g in mat.index]
    rows = []
    seen = set()
    for ga in have_a:
        for gb in have_b:
            if ga == gb:
                continue
            key = tuple(sorted((ga, gb)))
            if key in seen:
                continue
            seen.add(key)
            x = mat.loc[ga].to_numpy(dtype=float)
            y = mat.loc[gb].to_numpy(dtype=float)
            if np.nanstd(x) == 0 or np.nanstd(y) == 0:
                rho, p = np.nan, np.nan
            else:
                rho, p = stats.spearmanr(x, y)
            rows.append({"gene_a": key[0], "gene_b": key[1], "n": x.size,
                         "spearman_rho": rho, "p_value": p})
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = bh_fdr(df["p_value"].to_numpy())
    return df


def background_pair_distribution(mat: pd.DataFrame, universe: pd.Index,
                                 n_pairs: int = 200_000, seed: int = RNG_SEED) -> np.ndarray:
    """Spearman rho for random gene pairs, to calibrate 'is this correlation big?'."""
    rng = np.random.default_rng(seed)
    uni = np.asarray(universe)
    sub = mat.loc[uni]
    ranks = sub.rank(axis=1).to_numpy(dtype=np.float32)
    ranks -= ranks.mean(axis=1, keepdims=True)
    norm = np.sqrt((ranks ** 2).sum(axis=1, keepdims=True))
    norm[norm == 0] = np.nan
    ranks /= norm
    i = rng.integers(0, len(uni), size=n_pairs)
    j = rng.integers(0, len(uni), size=n_pairs)
    ok = i != j
    return np.nansum(ranks[i[ok]] * ranks[j[ok]], axis=1)


def signed_percentile(rho: float, background: np.ndarray) -> float:
    """Percentile of rho within the background distribution (0-100)."""
    if not np.isfinite(rho):
        return np.nan
    return float((background < rho).mean() * 100.0)


def expression_bins(mat: pd.DataFrame, universe: pd.Index, n_bins: int = 20) -> pd.Series:
    mean_expr = mat.loc[universe].mean(axis=1)
    return pd.qcut(mean_expr.rank(method="first"), n_bins, labels=False)


def matched_random_sets(bins: pd.Series, genes: list[str], n_sets: int,
                        seed: int = RNG_SEED) -> list[list[str]]:
    """Random gene sets matched to `genes` on mean-expression bin."""
    rng = np.random.default_rng(seed)
    by_bin = {b: np.asarray(idx) for b, idx in bins.groupby(bins).groups.items()}
    target_bins = [bins[g] for g in genes if g in bins.index]
    out = []
    for _ in range(n_sets):
        pick = []
        for b in target_bins:
            pool = by_bin[b]
            pick.append(str(rng.choice(pool)))
        out.append(pick)
    return out


def mean_pairwise_rho(mat: pd.DataFrame, genes: list[str]) -> float:
    have = [g for g in genes if g in mat.index]
    if len(have) < 2:
        return np.nan
    r = mat.loc[have].T.rank().corr(method="pearson").to_numpy()
    iu = np.triu_indices_from(r, k=1)
    return float(np.nanmean(r[iu]))


def module_score(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean of per-gene z-scores across samples (genes missing are dropped)."""
    have = [g for g in genes if g in mat.index]
    sub = mat.loc[have]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def partial_spearman(x: np.ndarray, y: np.ndarray, covars: np.ndarray):
    """Spearman correlation of x,y after linearly removing rank(covars)."""
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    C = np.column_stack([stats.rankdata(covars[:, k]) for k in range(covars.shape[1])])
    C = np.column_stack([np.ones(len(rx)), C])
    beta_x, *_ = np.linalg.lstsq(C, rx, rcond=None)
    beta_y, *_ = np.linalg.lstsq(C, ry, rcond=None)
    ex = rx - C @ beta_x
    ey = ry - C @ beta_y
    r, _ = stats.pearsonr(ex, ey)
    n = len(rx)
    k = C.shape[1] - 1
    dof = n - k - 2
    if dof <= 0 or abs(r) >= 1:
        return float(r), np.nan
    t = r * np.sqrt(dof / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    return float(r), float(p)


def fisher_ci(rho: float, n: int, alpha: float = 0.05):
    """Fisher z CI for a (Spearman) correlation."""
    if not np.isfinite(rho) or n < 5 or abs(rho) >= 1:
        return (np.nan, np.nan)
    z = np.arctanh(rho)
    se = 1.06 / np.sqrt(n - 3)  # Bonett-Wright SE for Spearman
    zc = stats.norm.ppf(1 - alpha / 2)
    return (float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se)))


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n_a, n_b = len(a), len(b)
    if n_a == 0 or n_b == 0:
        return np.nan
    # rank-based O(n log n) formulation
    all_v = np.concatenate([a, b])
    r = stats.rankdata(all_v)
    r_a = r[:n_a].sum()
    u_a = r_a - n_a * (n_a + 1) / 2
    return float(2 * u_a / (n_a * n_b) - 1)


def write_table(df: pd.DataFrame, name: str) -> str:
    path = os.path.join(TABLES, name)
    df.to_csv(path, index=False)
    return path


def write_json(obj, name: str) -> str:
    path = os.path.join(TABLES, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    return path


@dataclass
class Verdict:
    prediction: str
    dataset: str
    statistic: str
    value: float
    direction_ok: bool
    note: str = ""
