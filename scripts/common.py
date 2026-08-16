"""
Shared utilities for CLAIM B6 spatial colocalization analysis.

Claim B6: In public Visium / GeoMx spatial data, spots/regions high in the
epithelial-tumor markers CLDN4 / TACSTD2 anti-colocalize (spatially exclude)
T cells and B cells.

Design of the test (kept deliberately simple and honest):
  * Per capture area / section (Visium) or across AOIs (GeoMx) we build two
    signature scores per spatial unit:
       - epithelial score   = CLDN4 / TACSTD2
       - T-cell score        = T lymphocyte markers
       - B-cell score        = B lymphocyte markers
  * "Colocalization" between two signatures is measured as the Spearman rank
    correlation of the two scores ACROSS the spatial units of one sample.
    Spots are physical locations, so a negative correlation across spots means
    the two programs occupy different spots => spatial anti-colocalization.
  * Negative rho supports the claim; positive rho contradicts it.
  * We also report a top-tertile contrast (immune score in CLDN4/TACSTD2-high
    spots vs the rest) with a rank-biserial effect size, which is easier to
    interpret than rho alone.

Everything is reported per-sample and then aggregated with a sign test, so a
reader can see exactly how consistent (or not) the effect is.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

# ---- marker panels -------------------------------------------------------
# The claim names these two genes specifically; this is the primary signature.
EPI_MARKERS = ["CLDN4", "TACSTD2"]

# Broader epithelial panel used only as a sensitivity check.
EPI_PANEL_BROAD = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]

T_MARKERS = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "TRAC", "CD247", "IL7R"]
B_MARKERS = ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22"]


def present(genes, var_names) -> list[str]:
    s = set(var_names)
    return [g for g in genes if g in s]


def score_signature(adata, genes, name):
    """Background-corrected signature score (scanpy score_genes) on log-norm data.

    Returns the list of genes actually used (intersection with the panel).
    Writes adata.obs[name]. If no gene present, writes NaN.
    """
    used = present(genes, adata.var_names)
    if len(used) == 0:
        adata.obs[name] = np.nan
        return used
    sc.tl.score_genes(adata, used, score_name=name, use_raw=False)
    return used


def rank_biserial(high, rest):
    """Rank-biserial correlation from Mann-Whitney U (effect size in [-1, 1]).

    Positive => 'high' group tends to exceed 'rest'. For the claim we expect the
    immune score to be LOWER in epithelial-high spots => negative value.
    """
    high = np.asarray(high, float)
    rest = np.asarray(rest, float)
    n1, n2 = len(high), len(rest)
    if n1 == 0 or n2 == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(high, rest, alternative="two-sided")
    rb = 2.0 * u / (n1 * n2) - 1.0
    return rb, p


def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(a[m], b[m])
    return float(rho), float(p), int(m.sum())


def analyze_sample(adata, sample_id, extra=None):
    """Given a log-normalized AnnData for ONE sample, compute claim-B6 metrics.

    Returns a dict row. `extra` is merged in (dataset, tissue, platform, ...).
    """
    row = dict(sample=sample_id, n_spots=int(adata.n_obs))
    if extra:
        row.update(extra)

    epi_used = score_signature(adata, EPI_MARKERS, "epi_score")
    epi_broad_used = score_signature(adata, EPI_PANEL_BROAD, "epi_broad_score")
    t_used = score_signature(adata, T_MARKERS, "t_score")
    b_used = score_signature(adata, B_MARKERS, "b_score")

    # combined lymphocyte score
    tb_genes = present(T_MARKERS + B_MARKERS, adata.var_names)
    if tb_genes:
        sc.tl.score_genes(adata, tb_genes, score_name="tb_score", use_raw=False)
    else:
        adata.obs["tb_score"] = np.nan

    row["epi_markers_used"] = ",".join(epi_used)
    row["t_markers_used"] = ",".join(t_used)
    row["b_markers_used"] = ",".join(b_used)
    row["cldn4_detected_frac"] = _det_frac(adata, "CLDN4")
    row["tacstd2_detected_frac"] = _det_frac(adata, "TACSTD2")

    epi = adata.obs["epi_score"].values
    for tgt, col in [("T", "t_score"), ("B", "b_score"), ("TB", "tb_score"),
                     ("epi_broad", "epi_broad_score")]:
        rho, p, n = spearman(epi, adata.obs[col].values)
        row[f"rho_epi_vs_{tgt}"] = rho
        row[f"p_epi_vs_{tgt}"] = p
    row["n_used"] = n

    # top-tertile contrast on the combined lymphocyte score
    if np.isfinite(epi).sum() >= 30 and adata.obs["tb_score"].notna().any():
        thr = np.nanquantile(epi, 2 / 3)
        hi = adata.obs["tb_score"].values[epi >= thr]
        lo = adata.obs["tb_score"].values[epi < thr]
        rb, p = rank_biserial(hi, lo)
        row["rb_TB_in_epihigh"] = rb
        row["p_rb_TB_in_epihigh"] = p
        row["n_epihigh"] = int(np.sum(epi >= thr))
    else:
        row["rb_TB_in_epihigh"] = np.nan
        row["p_rb_TB_in_epihigh"] = np.nan
        row["n_epihigh"] = 0
    return row


def _det_frac(adata, gene):
    if gene not in adata.var_names:
        return 0.0
    x = adata[:, gene].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    return float(np.mean(x > 0))


def sign_test_summary(rhos):
    """Two-sided sign test that the median correlation differs from 0."""
    rhos = np.asarray([r for r in rhos if np.isfinite(r)], float)
    n = len(rhos)
    neg = int(np.sum(rhos < 0))
    pos = int(np.sum(rhos > 0))
    if n == 0:
        return dict(n=0, n_neg=0, n_pos=0, median_rho=np.nan, sign_p=np.nan)
    k = min(neg, pos)
    # two-sided sign test via binomial
    p = float(stats.binomtest(k, n, 0.5).pvalue) if n > 0 else np.nan
    # if all same sign, binomtest(min,n) gives correct small p
    return dict(n=n, n_neg=neg, n_pos=pos, median_rho=float(np.median(rhos)),
                mean_rho=float(np.mean(rhos)), sign_p=p)
