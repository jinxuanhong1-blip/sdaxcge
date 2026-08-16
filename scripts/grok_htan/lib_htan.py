"""Shared helpers for HTAN/HCA TACSTD2/CLDN4 vs immune analyses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr

# Human gene IDs (do not invent accessions; these are official HGNC/Ensembl).
TACSTD2_IDS = {
    "TACSTD2", "ENSG00000184292", "TROP2", "GA733-1", "M1S1",
    "ENSMUSG00000051397",  # mouse Tacstd2 (verified in HTAN Visium var)
}
CLDN4_IDS = {
    "CLDN4", "ENSG00000189143", "CPE-R", "CPER",
    "ENSMUSG00000047501",  # mouse Cldn4 (verified in HTAN Visium var)
}

# Compact immune / exclusion modules used across sc and spatial.
# Coverage is reported; missing genes are skipped (never fabricated).
SIGNATURES = {
    "CD8_Tcell": ["CD8A", "CD8B", "CD8B1", "CD3E", "CD3D"],
    "CYT": ["GZMA", "PRF1"],
    "Cytotoxic_effector": ["GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG"],
    "IFNG_6gene": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"],
    "TLS_12chemokine": [
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18",
        "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
    ],
    "TGFB_exclusion": ["TGFB1", "TGFBR2", "ACTA2", "COL1A1", "COL1A2", "FN1", "SMAD3"],
    "Immune_general": ["PTPRC", "CD3E", "CD8A", "CD4", "NKG7", "MS4A1", "CD79A", "CD68", "LYZ"],
    "B_cell": ["MS4A1", "CD79A", "CD19", "CD79B"],
    "Myeloid": ["CD68", "LYZ", "CSF1R", "ITGAM", "CD14"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
}

EPITHELIAL_TOKENS = (
    "epithelial", "alveolar", "club", "goblet", "ciliated", "basal",
    "ionocyte", "tuft", "pneumocyte", "at1", "at2", "atii", "secretory",
    "tumor", "malignant", "cancer", "nsclc", "sclc", "luad", "lusc",
    "neuroendocrine", "keratinocyte", "airway",
)
IMMUNE_TOKENS = (
    "t cell", "t-cell", "cd4", "cd8", "nk", "b cell", "b-cell", "plasma",
    "macrophage", "monocyte", "dendritic", "neutrophil", "mast", "basophil",
    "lymphocyte", "myeloid", "microglia", "ilc", "treg", "immune",
)
CD8_TOKENS = ("cd8", "cytotoxic t")
T_TOKENS = ("t cell", "t-cell", "cd4", "cd8", "treg", "regulatory t", "thymocyte")


def repo_paths():
    root = Path(__file__).resolve().parents[2]
    return {
        "root": root,
        "notes": root / "notes" / "grok_htan",
        "scripts": root / "scripts" / "grok_htan",
        "results": root / "results" / "grok_htan",
        "tables": root / "results" / "grok_htan" / "tables",
        "figures": root / "results" / "grok_htan" / "figures",
        "data": Path("/tmp/grok_htan_data"),
    }


def ensure_dirs():
    p = repo_paths()
    for k in ("notes", "results", "tables", "figures", "data"):
        p[k].mkdir(parents=True, exist_ok=True)
    return p


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return out
    pv = p[mask]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n, dtype=float)
    tmp[order] = q
    out[mask] = tmp
    return out


def spearman_safe(x, y, min_n=8):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < min_n:
        return {"n": n, "rho": np.nan, "p": np.nan}
    if np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def mwu_safe(a, b, min_n=5):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < min_n or b.size < min_n:
        return {
            "n_a": int(a.size),
            "n_b": int(b.size),
            "median_a": float(np.median(a)) if a.size else np.nan,
            "median_b": float(np.median(b)) if b.size else np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial": np.nan,
        }
    res = mannwhitneyu(a, b, alternative="two-sided")
    n1, n2 = a.size, b.size
    rbc = 1.0 - (2.0 * res.statistic) / (n1 * n2)
    return {
        "n_a": int(n1),
        "n_b": int(n2),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "rank_biserial": float(rbc),
    }


def _to_dense_1d(x):
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    return np.ravel(x).astype(float)


def extract_vector(adata, idx):
    return _to_dense_1d(adata.X[:, idx])


def resolve_gene_index(var, wanted):
    """Return first matching var position for a set of gene aliases."""
    wanted_u = {str(x).upper() for x in wanted}
    index_vals = pd.Index(var.index.astype(str))
    for i, g in enumerate(index_vals):
        if g.upper() in wanted_u or g.split(".")[0].upper() in wanted_u:
            return i, g
    for col in ("feature_name", "gene_symbols", "gene_symbol", "symbol", "name"):
        if col in var.columns:
            s = var[col].map(lambda v: "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
            for i, g in enumerate(s):
                if not g or g.lower() in ("nan", "none", "<na>"):
                    continue
                if g.upper() in wanted_u:
                    return i, g
    return None, None


def resolve_many(var, genes):
    found = {}
    missing = []
    for g in genes:
        idx, name = resolve_gene_index(var, {g})
        if idx is None:
            missing.append(g)
        else:
            found[g] = (idx, name)
    return found, missing


def mean_z_score(mat):
    """mat: cells x genes dense. Per-gene z, then mean. Constant genes -> 0."""
    mat = np.asarray(mat, dtype=float)
    if mat.ndim == 1:
        mat = mat[:, None]
    mu = np.nanmean(mat, axis=0)
    sd = np.nanstd(mat, axis=0)
    sd = np.where(sd == 0, np.nan, sd)
    z = (mat - mu) / sd
    z = np.where(np.isfinite(z), z, 0.0)
    return np.mean(z, axis=1)


def signature_scores(adata, signatures=SIGNATURES):
    rows = {}
    coverage = {}
    for name, genes in signatures.items():
        found, missing = resolve_many(adata.var, genes)
        coverage[name] = {
            "n_requested": len(genes),
            "n_found": len(found),
            "found": sorted(found),
            "missing": missing,
        }
        if not found:
            rows[name] = np.full(adata.n_obs, np.nan)
            continue
        cols = [found[g][0] for g in found]
        mat = adata.X[:, cols]
        if sparse.issparse(mat):
            mat = mat.toarray()
        rows[name] = mean_z_score(np.asarray(mat))
    return pd.DataFrame(rows, index=adata.obs_names), coverage


def classify_compartment(cell_type):
    s = str(cell_type).lower()
    if any(t in s for t in EPITHELIAL_TOKENS):
        return "epithelial"
    if any(t in s for t in IMMUNE_TOKENS):
        return "immune"
    return "other"


def is_cd8(cell_type):
    s = str(cell_type).lower()
    return any(t in s for t in CD8_TOKENS)


def is_tcell(cell_type):
    s = str(cell_type).lower()
    return any(t in s for t in T_TOKENS)


def pick_sample_col(obs):
    for c in (
        "sample", "sample_id", "orig.ident", "donor_id", "donor",
        "biosample_id", "library_id", "patient", "HTAN_Biospecimen_ID",
        "specimen", "orig_ident",
    ):
        if c in obs.columns and obs[c].nunique() >= 2:
            return c
    # fallback: any column with 3-200 levels
    for c in obs.columns:
        n = obs[c].nunique(dropna=True)
        if 3 <= n <= 400 and obs[c].dtype == object:
            return c
    return None


def pick_celltype_col(obs):
    for c in (
        "cell_type", "celltype", "CellType", "annotation", "ann_finest_level",
        "ann_level_2", "ann_level_3", "ann_coarse_for_GWAS_and_modeling",
        "author_cell_type", "cell_type_original", "majority_voting",
    ):
        if c in obs.columns:
            return c
    return None


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=str) + "\n")
