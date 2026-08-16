"""xCell (Aran, Hu & Butte, Genome Biology 2017).

Reimplementation of ``xCell::xCellAnalysis`` against the resource files shipped
in ``xCell.data`` (489 signatures over 64 cell types, the calibration values
``fv`` and the spillover matrix ``K``), which
``scripts/00_fetch_resources.py`` extracts from the package .rda.

The pipeline, in the original's order:

1. restrict to xCell's 10,808 gene universe (errors out below 5,000 genes);
2. rank each sample;
3. raw ssGSEA (``ssgsea.norm = FALSE``) over the 489 signatures;
4. subtract each signature's minimum across samples;
5. average the signatures belonging to the same cell type;
6. calibrate: ``((score - min) / 5000) ** fv[,2] / (fv[,3] * 2)``;
7. spillover compensation: non-negative constrained least squares against
   ``K * alpha`` with a unit diagonal, alpha = 0.5;
8. microenvironment aggregates.

Steps 4 and 6 subtract a **per-dataset minimum**, so xCell scores are
explicitly relative within the set of samples analysed together. Running xCell
twice on two halves of a cohort gives two incomparable score sets; this is the
single most important caveat for cross-cohort or batch-split designs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .ssgsea import ssgsea

__all__ = ["load_signatures", "xcell_analysis", "IMMUNE_TYPES", "STROMA_TYPES"]

IMMUNE_TYPES = [
    "B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Eosinophils",
    "Macrophages", "Monocytes", "Mast cells", "Neutrophils", "NK cells",
]
STROMA_TYPES = ["Adipocytes", "Endothelial cells", "Fibroblasts"]


def load_signatures(path: str | Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def _nnls_compensate(K: np.ndarray, y: np.ndarray) -> np.ndarray:
    """``pracma::lsqlincon(K, y, lb = 0)`` -- least squares with x >= 0."""
    from scipy.optimize import lsq_linear

    res = lsq_linear(K, y, bounds=(0, np.inf), method="trf", tol=1e-12)
    return res.x


def xcell_analysis(
    expr: pd.DataFrame,
    signatures: dict[str, list[str]],
    genes: list[str],
    spill_K: pd.DataFrame,
    fit_values: pd.DataFrame,
    alpha: float = 0.5,
    scale: bool = True,
    cell_types_use: list[str] | None = None,
    return_raw: bool = False,
):
    """Run the xCell pipeline. Returns cell types x samples.

    ``cell_types_use`` restricts the compensation step to the cell types you
    expect to be present. The original documentation recommends this because the
    spillover step over-compensates when it is asked to explain a tumour biopsy
    with all 64 types (including neurons, hepatocytes and osteoblasts).
    """
    shared = [g for g in expr.index if g in set(genes)]
    if len(shared) < 5000:
        raise ValueError(
            f"only {len(shared)} of xCell's {len(genes)} genes present; "
            "xCell refuses to run below 5,000 genes"
        )
    work = expr.loc[shared]

    # The original ranks before calling GSVA; ssGSEA re-ranks internally, which
    # is idempotent, but we keep the step so the tie structure matches exactly.
    ranked = work.rank(axis=0, method="average")

    raw = ssgsea(ranked, signatures, alpha=0.25, normalize="none", tie_method="average_int")
    raw = raw.sub(raw.min(axis=1), axis=0)

    cell_type = [name.split("%")[0] for name in raw.index]
    by_type = raw.groupby(cell_type, sort=True).mean()

    rows = [r for r in by_type.index if r in fit_values.index]
    tscores = by_type.loc[rows]
    tscores = tscores.sub(tscores.min(axis=1), axis=0) / 5000.0
    tscores[tscores < 0] = 0.0

    fv = fit_values.loc[rows]
    power = fv.iloc[:, 1].astype(float)
    divisor = fv.iloc[:, 2].astype(float) if scale else pd.Series(1.0, index=rows)
    tscores = tscores.pow(power, axis=0).div(divisor * 2.0, axis=0)

    if cell_types_use is not None:
        missing = set(cell_types_use) - set(spill_K.index)
        if missing:
            raise ValueError(f"unknown cell types: {sorted(missing)}")
        tscores = tscores.loc[[c for c in cell_types_use if c in tscores.index]]

    use = [r for r in tscores.index if r in spill_K.index]
    K = spill_K.loc[use, use].to_numpy(dtype=float) * alpha
    np.fill_diagonal(K, 1.0)

    adjusted = np.column_stack(
        [_nnls_compensate(K, tscores[col].to_numpy(dtype=float)) for col in tscores.columns]
    )
    adjusted[adjusted < 0] = 0.0
    out = pd.DataFrame(adjusted, index=use, columns=tscores.columns)

    if cell_types_use is None:
        immune = [c for c in IMMUNE_TYPES if c in out.index]
        stroma = [c for c in STROMA_TYPES if c in out.index]
        if len(immune) == len(IMMUNE_TYPES) and len(stroma) == len(STROMA_TYPES):
            immune_score = out.loc[immune].sum() / 1.5
            stroma_score = out.loc[stroma].sum() / 2.0
            out.loc["ImmuneScore"] = immune_score
            out.loc["StromaScore"] = stroma_score
            out.loc["MicroenvironmentScore"] = immune_score + stroma_score

    out.attrs["n_genes_used"] = len(shared)
    if return_raw:
        return out, raw
    return out
