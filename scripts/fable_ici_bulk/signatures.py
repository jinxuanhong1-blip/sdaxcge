"""Gene signatures and purity-aware statistics for the TACSTD2/CLDN4 ICI analysis.

Rationale: TACSTD2 (TROP2) and CLDN4 are epithelial/tumor-cell genes. In bulk
tumor RNA their apparent level is diluted by immune/stromal infiltrate, so a raw
comparison between responders (often immune-hot, low purity) and non-responders
confounds tumor-intrinsic expression with tumor purity. We therefore:
  (1) build cytotoxic (CD8 / NK) and epithelial/purity signatures,
  (2) estimate tumor-intrinsic TACSTD2/CLDN4 by adjusting for epithelial content,
  (3) use partial Spearman (controlling epithelial/purity) for correlations.
"""
import numpy as np
from scipy import stats

# --- immune cytotoxic signatures -------------------------------------------
CD8_TCELL = ["CD8A", "CD8B", "GZMK", "EOMES"]
NK_CELL = ["KLRD1", "KLRF1", "NCR1", "NKG7", "GNLY", "KLRC1", "KIR2DL4", "KLRK1"]
CYTOTOXIC = ["GZMA", "GZMB", "GZMH", "PRF1", "IFNG", "NKG7", "GNLY"]
CD8_NK = sorted(set(CD8_TCELL + NK_CELL))

# --- epithelial / tumor-content signature (purity proxy) --------------------
# CLDN4 and TACSTD2 are deliberately excluded to avoid circularity.
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "ELF3",
              "CLDN3", "CLDN7", "KRT17"]

# --- leukocyte / immune-content signature (inverse purity) ------------------
LEUKOCYTE = ["PTPRC", "CD3D", "CD3E", "CD2", "CD52", "CD53", "CD48", "LCP2"]

SIGNATURES = {
    "CD8_Tcell": CD8_TCELL,
    "NK_cell": NK_CELL,
    "Cytotoxic": CYTOTOXIC,
    "CD8_NK": CD8_NK,
    "Epithelial": EPITHELIAL,
    "Leukocyte": LEUKOCYTE,
}

TARGET_GENES = ["TACSTD2", "CLDN4"]
ENSG = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}

# Stable (unversioned) Ensembl IDs for all target + signature genes, used to map
# the Ensembl-indexed GSE135222 matrix back to HGNC symbols.
SYMBOL2ENSG = {
    "TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563", "CD8B": "ENSG00000172116",
    "GZMK": "ENSG00000113088", "EOMES": "ENSG00000163508",
    "KLRD1": "ENSG00000134539", "KLRF1": "ENSG00000150045",
    "NCR1": "ENSG00000189430", "NKG7": "ENSG00000105374",
    "GNLY": "ENSG00000115523", "KLRC1": "ENSG00000134545",
    "KIR2DL4": "ENSG00000189013", "KLRK1": "ENSG00000213809",
    "GZMA": "ENSG00000145649", "GZMB": "ENSG00000100453",
    "GZMH": "ENSG00000100450", "PRF1": "ENSG00000180644",
    "IFNG": "ENSG00000111537", "EPCAM": "ENSG00000119888",
    "KRT8": "ENSG00000170421", "KRT18": "ENSG00000111057",
    "KRT19": "ENSG00000171345", "KRT7": "ENSG00000135480",
    "CDH1": "ENSG00000039068", "ELF3": "ENSG00000163435",
    "CLDN3": "ENSG00000165215", "CLDN7": "ENSG00000181885",
    "KRT17": "ENSG00000128422", "PTPRC": "ENSG00000081237",
    "CD3D": "ENSG00000167286", "CD3E": "ENSG00000198851",
    "CD2": "ENSG00000116824", "CD52": "ENSG00000169442",
    "CD53": "ENSG00000143119", "CD48": "ENSG00000117091",
    "LCP2": "ENSG00000043462",
}
ENSG2SYMBOL = {v: k for k, v in SYMBOL2ENSG.items()}


def zscore(v):
    v = np.asarray(v, float)
    sd = np.nanstd(v)
    if sd == 0 or np.isnan(sd):
        return np.zeros_like(v)
    return (v - np.nanmean(v)) / sd


def signature_score(expr_df, genes):
    """Mean z-score across available signature genes (rows=genes, cols=samples).

    expr_df is assumed already on a log-like scale. Returns a Series over samples
    plus the list of genes actually found.
    """
    found = [g for g in genes if g in expr_df.index]
    if not found:
        return None, found
    z = expr_df.loc[found].apply(lambda r: zscore(r.values), axis=1,
                                 result_type="expand")
    z.columns = expr_df.columns
    return z.mean(axis=0), found


def _rank(a):
    return stats.rankdata(a)


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def partial_spearman(x, y, covars):
    """Spearman partial correlation of x,y controlling for covars (list of arrays).

    Computed as Pearson correlation of the residuals of rank(x) and rank(y)
    after linear regression on the ranked covariates. p-value from a t-test with
    df = n - 2 - k (k = number of covariates).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    C = [np.asarray(c, float) for c in covars]
    m = ~(np.isnan(x) | np.isnan(y))
    for c in C:
        m &= ~np.isnan(c)
    x, y = x[m], y[m]
    C = [c[m] for c in C]
    n = m.sum()
    k = len(C)
    if n - 2 - k < 1:
        return float("nan"), float("nan"), int(n)
    rx, ry = _rank(x), _rank(y)
    Z = np.column_stack([np.ones(n)] + [_rank(c) for c in C])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    ex = rx - Z @ bx
    ey = ry - Z @ by
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan"), float("nan"), int(n)
    r = np.corrcoef(ex, ey)[0, 1]
    df = n - 2 - k
    t = r * np.sqrt(df / max(1e-12, (1 - r * r)))
    p = 2 * stats.t.sf(abs(t), df)
    return float(r), float(p), int(n)


def mannwhitney_directional(high_group, low_group):
    """Two-sided + one-sided (greater) Mann-Whitney; effect = AUC(high>low).

    Returns dict. AUC>0.5 means values in `high_group` tend to exceed `low_group`.
    """
    a = np.asarray(high_group, float)
    b = np.asarray(low_group, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    U, p2 = stats.mannwhitneyu(a, b, alternative="two-sided")
    _, p_greater = stats.mannwhitneyu(a, b, alternative="greater")
    auc = U / (len(a) * len(b))
    return {
        "n_high": len(a), "n_low": len(b),
        "median_high": float(np.median(a)), "median_low": float(np.median(b)),
        "auc_high_gt_low": float(auc), "rank_biserial": float(2 * auc - 1),
        "p_two_sided": float(p2), "p_one_sided_greater": float(p_greater),
    }
