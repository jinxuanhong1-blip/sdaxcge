"""Shared configuration, gene sets, and helper functions.

Analysis slice: TACSTD2 (TROP2) and CLDN4 (Claudin-4) vs immune microenvironment /
immunotherapy response in EGFR/ALK-mutant lung cancer.

All paths are resolved relative to the repository root so scripts can be run from
anywhere.  Raw inputs live under data/fable_egfr_alk/ (git-ignored, reproducible via
the download commands recorded in notes/); all committed outputs go under
results/fable_egfr_alk/.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "fable_egfr_alk")
RESULTS = os.path.join(ROOT, "results", "fable_egfr_alk")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")
PROCESSED = os.path.join(RESULTS, "processed")

for _d in (RESULTS, TABLES, FIGURES, PROCESSED):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------------------
# Genes of interest
# ---------------------------------------------------------------------------
# Target genes (ADC / TJ targets under study)
TARGETS = {
    "TACSTD2": "ENSG00000184292",  # TROP2
    "CLDN4": "ENSG00000189143",    # Claudin-4
}

# Driver genes
DRIVERS = {
    "EGFR": "ENSG00000146648",
    "ALK": "ENSG00000171094",
}

# Immune marker genes -> Ensembl (unversioned)
IMMUNE_GENES = {
    "CD8A": "ENSG00000153563",
    "CD8B": "ENSG00000172116",
    "GZMA": "ENSG00000145649",
    "GZMB": "ENSG00000100453",
    "PRF1": "ENSG00000180644",
    "IFNG": "ENSG00000111537",
    "TBX21": "ENSG00000073861",
    "CXCL9": "ENSG00000138755",
    "CXCL10": "ENSG00000169245",
    "CXCL11": "ENSG00000169248",
    "STAT1": "ENSG00000115415",
    "IDO1": "ENSG00000131203",
    "HLA-DRA": "ENSG00000204287",
    "CD274": "ENSG00000120217",   # PD-L1
    "PDCD1": "ENSG00000188389",   # PD-1
    "CTLA4": "ENSG00000163599",
    "LAG3": "ENSG00000089692",
    "HAVCR2": "ENSG00000135077",  # TIM-3
    "TIGIT": "ENSG00000181847",
    "FOXP3": "ENSG00000049768",
    "CD4": "ENSG00000010610",
    "CD3D": "ENSG00000167286",
    "CD3E": "ENSG00000198851",
    "GZMK": "ENSG00000113088",
    "NKG7": "ENSG00000105374",
    "CCL5": "ENSG00000271503",
    "CCR5": "ENSG00000160791",
    "ITGAX": "ENSG00000140678",   # CD11c / myeloid
    "CD68": "ENSG00000129226",
    "CD163": "ENSG00000177575",   # M2 macrophage
    "MRC1": "ENSG00000260314",    # CD206 / M2
}

# ---------------------------------------------------------------------------
# Immune signatures (built from IMMUNE_GENES symbols)
# ---------------------------------------------------------------------------
# T-cell inflamed / IFN-gamma (Ayers et al. 2017, JCI) 18-gene GEP -- subset that
# is available in our marker panel.
GEP_TCELL_INFLAMED = [
    "CD8A", "CD274", "CD27", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
    "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT", "CCL5",
]  # restricted to available symbols at scoring time
GEP_IFNG_6 = ["IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"]

# Cytolytic activity (Rooney et al. 2015): geometric mean of GZMA & PRF1.
CYT_GENES = ["GZMA", "PRF1"]

# Broad effector CD8 / cytotoxic panel
CD8_EFFECTOR = ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "IFNG"]

# Inhibitory checkpoint panel
CHECKPOINTS = ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"]

# Broad immune infiltration panel (T + APC + myeloid markers)
BROAD_IMMUNE = [
    "CD8A", "CD8B", "CD3D", "CD3E", "CD4", "GZMA", "GZMB", "PRF1", "IFNG",
    "CXCL9", "CXCL10", "STAT1", "HLA-DRA", "CD274", "CD68", "ITGAX",
]

SIGNATURES = {
    "IFNG_6gene": GEP_IFNG_6,
    "CYT": CYT_GENES,
    "CD8_effector": CD8_EFFECTOR,
    "Checkpoint": CHECKPOINTS,
    "Broad_immune": BROAD_IMMUNE,
}

RANDOM_SEED = 20260816


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def strip_ens_version(idx: pd.Index) -> pd.Index:
    """Remove the ``.NN`` version suffix from Ensembl gene IDs."""
    return idx.to_series().str.replace(r"\.\d+$", "", regex=True).values


def zscore(x: pd.Series) -> pd.Series:
    sd = x.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return x * 0.0
    return (x - x.mean()) / sd


def signature_score(expr: pd.DataFrame, genes) -> pd.Series:
    """Mean of per-sample z-scored genes (samples x genes input).

    ``expr`` must have samples as the index and gene *symbols* as columns.
    Missing genes are silently skipped.  For CYT we use the log2 of the
    geometric mean convention already satisfied by mean-of-z on log-scale data.
    """
    avail = [g for g in genes if g in expr.columns]
    if not avail:
        return pd.Series(np.nan, index=expr.index)
    z = expr[avail].apply(zscore, axis=0)
    return z.mean(axis=1)


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Return BH-FDR adjusted q-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q
