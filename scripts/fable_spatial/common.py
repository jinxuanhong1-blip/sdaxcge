"""Shared utilities for the fable_spatial slice.

Marker sets, Visium loading, spot-neighborhood construction and the
statistics used across datasets (Spearman / partial Spearman on ranks,
naive within-section permutation, cross-section meta-tests).
"""

import gzip
import io
import os
import tarfile

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results", "fable_spatial")

# Target genes
TARGETS = ["TACSTD2", "CLDN4"]

# Epithelial score: excludes TACSTD2/CLDN4 and their close paralogs (CLDN3/7)
# to avoid circularity with the target genes.
EPI_MARKERS = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "CDH1", "ELF3", "SFN"]

# Broad immune score (lineage markers across T/NK/B/myeloid)
IMMUNE_MARKERS = [
    "PTPRC", "CD3E", "CD3D", "CD2", "TRAC", "IL7R", "CD8A", "CCL5", "NKG7",
    "GZMB", "GZMA", "GZMK", "PRF1", "CD79A", "MS4A1", "LYZ", "C1QA", "C1QB",
    "CD68", "AIF1", "FCER1G", "CD14", "CD52",
]

# Effector / cytotoxic T-cell axis (used for the GeoMx stromal immune score too)
TEFF_MARKERS = [
    "CD8A", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "IFNG",
    "CXCL9", "CXCL10", "CXCL13", "CD3E", "CD3D", "CD2", "TRAC",
]


def ensure_results():
    os.makedirs(RESULTS, exist_ok=True)
    return RESULTS


# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------

def spearman(x, y):
    r, p = stats.spearmanr(x, y)
    return float(r), float(p)


def partial_spearman(x, y, covars):
    """Rank-based partial correlation of x and y given covariate matrix."""
    def _resid(v):
        ranks = stats.rankdata(v).astype(float)
        C = np.column_stack([np.ones(len(ranks))] + [stats.rankdata(c) for c in covars])
        beta, *_ = np.linalg.lstsq(C, ranks, rcond=None)
        return ranks - C @ beta

    rx, ry = _resid(x), _resid(y)
    n = len(x)
    r, _ = stats.pearsonr(rx, ry)
    k = len(covars)
    dof = n - 2 - k
    if dof <= 0 or abs(r) >= 1:
        return float(r), np.nan
    t = r * np.sqrt(dof / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    return float(r), float(p)


def perm_pvalue(x, y, stat_fn, n_perm=1000, seed=0):
    """Naive permutation p for association between x and y (shuffles y)."""
    rng = np.random.default_rng(seed)
    obs = stat_fn(x, y)
    null = np.empty(n_perm)
    y = np.asarray(y)
    for i in range(n_perm):
        null[i] = stat_fn(x, rng.permutation(y))
    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (n_perm + 1)
    return float(obs), float(p)


def stouffer(pvals, signs):
    """Two-sided Stouffer combination of two-sided p-values with effect signs."""
    pvals = np.asarray(pvals, dtype=float)
    signs = np.asarray(signs, dtype=float)
    ok = np.isfinite(pvals) & np.isfinite(signs) & (pvals > 0)
    pvals, signs = pvals[ok], signs[ok]
    if len(pvals) == 0:
        return np.nan, np.nan
    z = stats.norm.isf(pvals / 2) * np.sign(signs)
    zc = z.sum() / np.sqrt(len(z))
    return float(zc), float(2 * stats.norm.sf(abs(zc)))


def wilcoxon_rhos(rhos):
    """Wilcoxon signed-rank test that per-section correlations differ from 0."""
    rhos = np.asarray([r for r in rhos if np.isfinite(r)])
    if len(rhos) < 5:
        return np.nan, np.nan
    try:
        w, p = stats.wilcoxon(rhos)
        return float(w), float(p)
    except ValueError:
        return np.nan, np.nan


# ---------------------------------------------------------------------------
# Visium loading
# ---------------------------------------------------------------------------

def load_10x_h5(path):
    import h5py
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        # stored as genes x barcodes in CSC (columns = barcodes)
        X = sp.csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]),
            shape=g["shape"][:],
        ).T.tocsr()  # -> barcodes x genes
        barcodes = [b.decode() for b in g["barcodes"][:]]
        names = [n.decode() for n in g["features/name"][:]]
    return X, barcodes, names


def load_mtx_triplet(mtx, features, barcodes):
    X = sp.csr_matrix(_read_mtx(mtx)).T  # genes x spots -> spots x genes
    feats = pd.read_csv(features, sep="\t", header=None, compression="gzip")
    bcs = pd.read_csv(barcodes, sep="\t", header=None, compression="gzip")[0].tolist()
    names = feats[1].tolist()
    return X, bcs, names


def _read_mtx(path):
    from scipy.io import mmread
    with gzip.open(path, "rb") as f:
        return mmread(f)


def load_positions_from_tar(tar_path):
    """Extract tissue positions from a 10x `spatial.tar` archive."""
    with tarfile.open(tar_path) as t:
        member = None
        for m in t.getmembers():
            base = os.path.basename(m.name)
            if base in ("tissue_positions.csv", "tissue_positions_list.csv"):
                member = m
                break
        if member is None:
            raise FileNotFoundError(f"no tissue positions in {tar_path}")
        raw = t.extractfile(member).read()
    first = raw.split(b"\n", 1)[0]
    header = 0 if b"barcode" in first else None
    df = pd.read_csv(io.BytesIO(raw), header=header)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    return df


def load_positions_csv_gz(path):
    df = pd.read_csv(path, header=None, compression="gzip")
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return df


# ---------------------------------------------------------------------------
# core per-section analysis
# ---------------------------------------------------------------------------

def collapse_duplicate_genes(X, names):
    """Sum counts of duplicated gene symbols; returns X (spots x genes), unique names."""
    names = pd.Index(names)
    if names.is_unique:
        return X, names.tolist()
    order = pd.factorize(names)[0]
    n_unique = order.max() + 1
    M = sp.csr_matrix(
        (np.ones(len(order)), (np.arange(len(order)), order)),
        shape=(len(order), n_unique),
    )
    Xu = X @ M
    uniq = pd.unique(names)
    return sp.csr_matrix(Xu), list(uniq)


def score(logX, names_idx, genes):
    """Mean z-scored log-normalized expression over available marker genes."""
    cols = [names_idx[g] for g in genes if g in names_idx]
    if not cols:
        return None, []
    sub = np.asarray(logX[:, cols].todense())
    mu, sd = sub.mean(0), sub.std(0)
    sd[sd == 0] = 1.0
    return ((sub - mu) / sd).mean(1), [g for g in genes if g in names_idx]


def analyze_section(X, barcodes, names, pos, section, min_counts=250, k_ring=6, n_perm=1000):
    """TACSTD2/CLDN4 (epithelial spots) vs neighborhood immune score.

    Returns (per-gene result rows, section QC dict).
    """
    from scipy.spatial import cKDTree

    pos = pos[pos.in_tissue == 1].set_index("barcode")
    keep = [i for i, b in enumerate(barcodes) if b in pos.index]
    X = X[keep]
    bcs = [barcodes[i] for i in keep]
    pos = pos.loc[bcs]

    X, names = collapse_duplicate_genes(X, names)

    tot = np.asarray(X.sum(1)).ravel()
    ok = tot >= min_counts
    X, pos, tot = X[ok], pos.iloc[ok], tot[ok]
    n_spots = X.shape[0]
    if n_spots < 100:
        return [], {"section": section, "n_spots": int(n_spots), "note": "too few spots"}

    # CP10K + log1p
    Xn = X.multiply(1e4 / tot[:, None]).tocsr()
    logX = Xn.copy()
    logX.data = np.log1p(logX.data)

    names_idx = {g: i for i, g in enumerate(names)}
    epi, _ = score(logX, names_idx, EPI_MARKERS)
    imm, _ = score(logX, names_idx, IMMUNE_MARKERS)
    teff, _ = score(logX, names_idx, TEFF_MARKERS)

    # neighborhood graph from pixel coordinates
    xy = pos[["pxl_row", "pxl_col"]].to_numpy(float)
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=k_ring + 1)
    # exclude self (col 0); drop neighbors farther than 2x median ring distance
    ring_d = np.median(d[:, 1])
    nbr_imm = np.full(n_spots, np.nan)
    nbr_teff = np.full(n_spots, np.nan)
    for i in range(n_spots):
        js = [j for j, dist in zip(idx[i, 1:], d[i, 1:]) if dist <= 2.0 * ring_d]
        if len(js) >= 3:
            nbr_imm[i] = imm[js].mean()
            nbr_teff[i] = teff[js].mean()

    epi_mask = (epi > np.median(epi)) & np.isfinite(nbr_imm)
    n_epi = int(epi_mask.sum())
    qc = {"section": section, "n_spots": int(n_spots), "n_epi_spots": n_epi,
          "median_counts": float(np.median(tot))}
    if n_epi < 50:
        qc["note"] = "too few epithelial spots"
        return [], qc

    rows = []
    for gene in TARGETS:
        if gene not in names_idx:
            continue
        expr = np.asarray(logX[:, names_idx[gene]].todense()).ravel()
        e = expr[epi_mask]
        pct = float((e > 0).mean())
        for nb_name, nb in [("immune_broad", nbr_imm), ("t_effector", nbr_teff)]:
            nbv = nb[epi_mask]
            rho, p = spearman(e, nbv)
            prho, pp = partial_spearman(e, nbv, [epi[epi_mask], imm[epi_mask]])
            _, permp = perm_pvalue(e, nbv, lambda a, b: stats.spearmanr(a, b)[0],
                                   n_perm=n_perm, seed=7)
            # hot vs cold neighborhood tertiles
            q1, q2 = np.quantile(nbv, [1 / 3, 2 / 3])
            lo, hi = e[nbv <= q1], e[nbv >= q2]
            try:
                u, up = stats.mannwhitneyu(hi, lo, alternative="two-sided")
            except ValueError:
                u, up = np.nan, np.nan
            delta = float(np.median(hi) - np.median(lo))
            rows.append(dict(section=section, gene=gene, neighborhood=nb_name,
                             n_epi_spots=n_epi, pct_detected=round(pct, 3),
                             spearman_rho=round(rho, 4), spearman_p=p,
                             partial_rho=round(prho, 4), partial_p=pp,
                             perm_p=permp,
                             hot_vs_cold_median_diff=round(delta, 4),
                             mannwhitney_p=up))
    return rows, qc
