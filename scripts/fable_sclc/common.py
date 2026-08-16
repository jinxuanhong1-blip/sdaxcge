"""Shared helpers for the SCLC TACSTD2/CLDN4 immunotherapy dataset slice.

Data directory is controlled by the SCLC_DATA_DIR environment variable
(default: /tmp/sclc_data). Only small derived tables/figures are written
into the repository (results/fable_sclc/).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(os.environ.get("SCLC_DATA_DIR", "/tmp/sclc_data"))
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "fable_sclc"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

TARGET_GENES = ["TACSTD2", "CLDN4"]

# SCLC subtype-defining transcription factors (Gay et al., Cancer Cell 2021)
SUBTYPE_TFS = ["ASCL1", "NEUROD1", "POU2F3"]

# T-effector signature used in IMpower studies (Fehrenbacher/Gay-style)
TEFF_GENES = ["CD8A", "GZMA", "GZMB", "PRF1", "IFNG", "CXCL9"]

# Cytolytic activity (Rooney et al. 2015)
CYT_GENES = ["GZMA", "PRF1"]

# MHC class I / antigen presentation machinery
APM_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"]

# Additional immune-context genes fetched for correlation analyses
EXTRA_IMMUNE_GENES = [
    "CD8B", "CD3D", "CD3E", "CD2", "GZMK", "CXCL10", "CXCL11", "IDO1",
    "STAT1", "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
    "TIGIT", "HLA-DRA", "HLA-E", "CIITA", "CD68", "CD163", "MRC1", "CSF1R",
]

# Neuroendocrine / epithelial context genes
CONTEXT_GENES = ["YAP1", "INSM1", "CHGA", "SYP", "DLL3", "GRP", "REST",
                 "VIM", "EPCAM", "KRT8", "KRT18", "NCAM1", "MKI67"]

PANEL_GENES = sorted(set(
    TARGET_GENES + SUBTYPE_TFS + TEFF_GENES + CYT_GENES + APM_GENES
    + EXTRA_IMMUNE_GENES + CONTEXT_GENES
))


def ensure_dirs():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def zscore(series_or_df, axis=0):
    x = series_or_df
    mu = x.mean(axis=axis)
    sd = x.std(axis=axis, ddof=1)
    return (x - mu) / sd


def _bimodal_threshold(values: pd.Series) -> float:
    """1D two-means threshold (midpoint of the two cluster centroids).

    ASCL1/NEUROD1/POU2F3 are essentially bimodal in SCLC (expressed vs
    silent), so a two-cluster split separates 'high' from 'low' samples.
    """
    from scipy.cluster.vq import kmeans2
    x = values.to_numpy(dtype=float)
    centroids, _ = kmeans2(x, k=2, seed=0, minit="++")
    return float(centroids.mean())


def assign_sclc_subtype(log_expr: pd.DataFrame) -> pd.Series:
    """Marker-based SCLC subtype approximation (samples x genes input).

    Approximation of the NMF classification of Gay et al. 2021 (subtypes
    dominated by ASCL1/NEUROD1/POU2F3; SCLC-I = low expression of all
    three TFs): each TF is dichotomized at a per-cohort bimodal (two-means)
    threshold. Samples with all three TFs in the low mode are SCLC-I;
    otherwise the subtype is the high-mode TF with the largest z-score.

    This is an approximation of the published NMF classification and is
    labelled 'subtype_approx' in all outputs.
    """
    z = zscore(log_expr[SUBTYPE_TFS])
    thresholds = {g: _bimodal_threshold(log_expr[g]) for g in SUBTYPE_TFS}
    high = pd.DataFrame({g: log_expr[g] > thresholds[g] for g in SUBTYPE_TFS})
    mapping = {"ASCL1": "SCLC-A", "NEUROD1": "SCLC-N", "POU2F3": "SCLC-P"}
    labels = []
    for i in log_expr.index:
        if not high.loc[i].any():
            labels.append("SCLC-I")
        else:
            cand = z.loc[i].where(high.loc[i])
            labels.append(mapping[cand.idxmax()])
    return pd.Series(labels, index=log_expr.index, name="subtype_approx")


def signature_score(log_expr: pd.DataFrame, genes) -> pd.Series:
    """Mean of per-gene z-scores over available signature genes."""
    genes = [g for g in genes if g in log_expr.columns]
    return zscore(log_expr[genes]).mean(axis=1)


def nmf_subtypes(log_expr_full: pd.DataFrame, n_top_genes=1250, seed=0):
    """NMF-based SCLC subtype assignment, replicating Gay et al. 2021.

    Input: samples x genes log2 expression (full transcriptome).
    Procedure: select the most variable genes, run NMF with k=4, assign
    each sample to its maximum-weight factor, then label factors by their
    ASCL1/NEUROD1/POU2F3 gene loadings (highest-loading factor per TF);
    the remaining factor is the TF-low, inflamed class SCLC-I.
    Returns (labels, sample_factor_weights).
    """
    from scipy.stats import spearmanr
    from sklearn.decomposition import NMF

    expr = log_expr_full.loc[:, log_expr_full.std() > 0]
    var = expr.var().nlargest(n_top_genes).index
    x = expr[var]
    model = NMF(n_components=4, init="nndsvda", max_iter=2000,
                random_state=seed)
    w = model.fit_transform(x.to_numpy())          # samples x factors
    h = model.components_                          # factors x genes
    # transfer factor scale into W so argmax across factors is meaningful
    w = w * np.linalg.norm(h, axis=1)

    # label factors by correlating factor weights with TF expression,
    # assigning greedily (strongest correlation first, unique factors)
    cors = []
    for tf, label in [("ASCL1", "SCLC-A"), ("NEUROD1", "SCLC-N"),
                      ("POU2F3", "SCLC-P")]:
        for k in range(4):
            rho = spearmanr(w[:, k], log_expr_full[tf]).statistic
            cors.append((rho, k, label))
    mapping, used_factors, used_labels = {}, set(), set()
    for rho, k, label in sorted(cors, reverse=True):
        if k not in used_factors and label not in used_labels:
            mapping[k] = label
            used_factors.add(k)
            used_labels.add(label)
    for k in range(4):
        mapping.setdefault(k, "SCLC-I")
    if sorted(mapping.values()) != ["SCLC-A", "SCLC-I", "SCLC-N", "SCLC-P"]:
        raise RuntimeError(f"ambiguous NMF factor labelling: {mapping}")

    labels = pd.Series([mapping[i] for i in w.argmax(axis=1)],
                       index=log_expr_full.index, name="subtype_nmf")
    weights = pd.DataFrame(
        w, index=log_expr_full.index,
        columns=[f"factor_{mapping[k]}" for k in range(4)])
    return labels, weights


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out
