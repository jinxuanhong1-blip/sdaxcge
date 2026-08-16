"""MCP-counter (Becht et al., Genome Biology 2016).

The estimator is deliberately simple: for each of the ten populations, the score
is the arithmetic mean of the **log2** expression of its transcriptomic markers.

Consequences that the playbook leans on:

* Scores are in log2 expression units, so they are comparable *across samples*
  for one population and **not** comparable *across populations* (fibroblast and
  CD8 markers have different baseline expression). Never read an MCP-counter
  table as a cell-fraction table.
* Because it is a mean of logs, MCP-counter is sensitive to library-size
  normalisation and to which markers are missing. The number of markers found
  is returned alongside the scores.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["load_markers", "mcp_counter", "MCP_POPULATIONS"]

MCP_POPULATIONS = [
    "T cells",
    "CD8 T cells",
    "Cytotoxic lymphocytes",
    "B lineage",
    "NK cells",
    "Monocytic lineage",
    "Myeloid dendritic cells",
    "Neutrophils",
    "Endothelial cells",
    "Fibroblasts",
]


def load_markers(path: str | Path, feature: str = "HUGO symbols") -> dict[str, list[str]]:
    """Load the MCP-counter marker table (``Signatures/genes.txt``)."""
    table = pd.read_csv(path, sep="\t", dtype=str)
    table.columns = [c.strip().strip('"') for c in table.columns]
    for col in table.columns:
        table[col] = table[col].str.strip().str.strip('"')
    if feature not in table.columns:
        raise ValueError(f"{feature!r} not in {list(table.columns)}")
    markers: dict[str, list[str]] = {}
    for population, block in table.groupby("Cell population", sort=False):
        genes = [g for g in block[feature].dropna().unique() if g]
        markers[population] = genes
    return markers


def mcp_counter(
    expr_log2: pd.DataFrame,
    markers: dict[str, list[str]],
    min_markers: int = 1,
) -> pd.DataFrame:
    """MCP-counter scores.

    Parameters
    ----------
    expr_log2
        genes x samples, **log2 scale**. Passing linear TPM here inflates the
        score of whichever population happens to contain the highest-expressed
        marker and is a silent error, so the scale is asserted.
    min_markers
        populations with fewer markers found are returned as NaN rather than
        as a mean over one or two genes.

    Returns
    -------
    DataFrame
        samples x populations. ``attrs["n_markers_found"]`` records coverage.
    """
    arr = expr_log2.to_numpy(dtype=float)
    if np.nanmax(arr) > 50:
        raise ValueError(
            "expression does not look log2 transformed (max > 50); MCP-counter "
            "is defined on log2 expression"
        )

    scores: dict[str, pd.Series] = {}
    found: dict[str, int] = {}
    for population, genes in markers.items():
        present = [g for g in dict.fromkeys(genes) if g in expr_log2.index]
        found[population] = len(present)
        if len(present) < min_markers:
            scores[population] = pd.Series(np.nan, index=expr_log2.columns)
        else:
            scores[population] = expr_log2.loc[present].mean(axis=0, skipna=True)

    out = pd.DataFrame(scores)
    out.attrs["n_markers_found"] = found
    out.attrs["n_markers_total"] = {k: len(set(v)) for k, v in markers.items()}
    return out
