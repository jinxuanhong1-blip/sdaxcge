"""Documented NicheNet-v2 ligand activity + MultiNicheNet-style ranks.

This is a Python reimplementation of the published scoring logic
(Browaeys et al. 2020; Bonte / Browaeys MultiNicheNet vignette).
It is **not** a run of R `nichenetr` or `multinichenetr`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

from lib_stats import minmax_scale


def ligand_activity(
    lt: pd.DataFrame,
    geneset: set[str],
    background: list[str],
    ligands: list[str],
) -> pd.DataFrame:
    """Pearson / AUROC / AUPR of prior target scores vs gene-set membership.

    `lt` is targets × ligands (NicheNet-v2 nsga2r matrix).
    """
    genes = [g for g in background if g in lt.index]
    if len(genes) < 20 or not ligands:
        return pd.DataFrame()
    y = np.array([1 if g in geneset else 0 for g in genes], dtype=int)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    if not use:
        return pd.DataFrame()
    X = lt.loc[genes, use].to_numpy(dtype=float)
    rows = []
    for j, lig in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s, y)
        try:
            auroc = float(roc_auc_score(y, s))
        except ValueError:
            auroc = np.nan
        try:
            aupr = float(average_precision_score(y, s))
        except ValueError:
            aupr = np.nan
        rows.append(
            {
                "ligand": lig,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "auroc": auroc,
                "aupr": aupr,
                "n_background": len(genes),
                "n_geneset": int(y.sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def ligand_activity_continuous(
    lt: pd.DataFrame,
    target_score: pd.Series,
    ligands: list[str],
) -> pd.DataFrame:
    """Pearson of prior columns vs a continuous receiver score (e.g. DE t)."""
    genes = [g for g in target_score.index if g in lt.index and np.isfinite(target_score[g])]
    if len(genes) < 20:
        return pd.DataFrame()
    y = target_score.reindex(genes).to_numpy(dtype=float)
    if np.allclose(y, y[0]):
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    if not use:
        return pd.DataFrame()
    X = lt.loc[genes, use].to_numpy(dtype=float)
    rows = []
    for j, lig in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s, y)
        rows.append(
            {
                "ligand": lig,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "n_background": len(genes),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def receptors_for(lr: pd.DataFrame, ligand: str) -> list[str]:
    return sorted(set(lr.loc[lr["from"] == ligand, "to"].astype(str)))


def prioritization_table(
    ligand_de: pd.DataFrame,
    activity: pd.DataFrame,
    receptor_score: pd.Series,
    activity_col: str = "pearson",
) -> pd.DataFrame:
    """Equal-weight MultiNicheNet-style rank (documented, not the R package).

    Features (min-max scaled, then averaged):
      1. sender ligand Δ (CLDN4-high − low), signed
      2. ligand activity (prior vs receiver program / DE)
      3. receptor expression in same-patient T/NK
      4. fraction of patients with ligand expressed in high senders
    """
    tab = ligand_de.merge(activity, on="ligand", how="inner")
    tab["receptor_score"] = tab["ligand"].map(receptor_score)
    parts = {
        "scaled_delta": minmax_scale(tab["median_delta"]),
        "scaled_activity": minmax_scale(tab[activity_col]),
        "scaled_receptor": minmax_scale(tab["receptor_score"]),
        "scaled_frac": minmax_scale(tab["frac_patients_expressed_high"]),
    }
    for name, s in parts.items():
        tab[name] = s
    tab["prioritization_score"] = tab[list(parts)].mean(axis=1)
    tab = tab.sort_values("prioritization_score", ascending=False).reset_index(drop=True)
    tab["rank"] = np.arange(1, len(tab) + 1)
    return tab
