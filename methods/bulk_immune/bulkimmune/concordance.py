"""Cross-method and cross-cohort concordance for the playbook's pre-specified axes.

A deconvolution playbook that never asks whether MCP-counter CD8, xCell CD8,
TIDE CD8 and the Ayers effector signature even rank samples the same way is
unfinished. These are the axes we actually use when we talk about TACSTD2 /
CLDN4 vs “exclusion” or “inflamed”.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

__all__ = ["AXES", "HEADLINE_PAIRS", "pairwise_spearman", "sign_concordance"]

# Short names as they appear AFTER the "method:" prefix strip, plus a few
# fully-prefixed names when two methods share a short name (CD8, Fibroblasts).
AXES: dict[str, list[str]] = {
    "CD8_cytotoxic": [
        "mcp:CD8 T cells",
        "xcell:CD8+ T-cells",
        "tide:CD8",
        "tide:CTL",
        "ssgsea:CD8_Teffector_Mariathasan",
        "ssgsea:CYT",
        "ssgsea:IFNG_Ayers6",
        "quantiseq_style:T.cells.CD8",
        "tip:Step4_CD8_T_cell",
        "tip:Step7_killing_of_cancer_cells",
    ],
    "IFN_inflamed": [
        "ssgsea:IFNG_Ayers6",
        "ssgsea:TcellInflamed_GEP18",
        "ssgsea:ExpandedImmune_Ayers18",
        "ssgsea:HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "tide:IFNG",
        "estimate:ImmuneScore",
    ],
    "exclusion_stroma": [
        "tide:Exclusion",
        "tide:CAF",
        "tide:MDSC",
        "mcp:Fibroblasts",
        "xcell:Fibroblasts",
        "estimate:StromalScore",
        "ssgsea:HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "ssgsea:HALLMARK_TGF_BETA_SIGNALING",
    ],
    "TLS_Bcell": [
        "ssgsea:TLS_12chemokine_Coppola",
        "ssgsea:TLS_Cabrita9",
        "ssgsea:TLS_Cabrita9_noY",
        "ssgsea:Bcell_follicular",
        "xcell:B-cells",
        "mcp:B lineage",
    ],
}

# Target–score pairs whose *sign* we track across cohorts.
HEADLINE_PAIRS: list[tuple[str, str]] = [
    ("TACSTD2", "xcell:CD8+ T-cells"),
    ("TACSTD2", "mcp:CD8 T cells"),
    ("TACSTD2", "tide:CD8"),
    ("TACSTD2", "ssgsea:CD8_Teffector_Mariathasan"),
    ("TACSTD2", "ssgsea:CYT"),
    ("TACSTD2", "tide:Exclusion"),
    ("TACSTD2", "estimate:ImmuneScore"),
    ("CLDN4", "xcell:CD8+ T-cells"),
    ("CLDN4", "mcp:CD8 T cells"),
    ("CLDN4", "tide:Exclusion"),
    ("CLDN4", "estimate:StromalScore"),
    ("CLDN4", "tide:MDSC"),
]


def pairwise_spearman(
    scores: pd.DataFrame,
    columns: Iterable[str],
    min_n: int = 8,
) -> pd.DataFrame:
    """Upper-triangle Spearman table among ``columns`` that exist in ``scores``."""
    cols = [c for c in columns if c in scores.columns]
    rows = []
    for i, a in enumerate(cols):
        for b in cols[i + 1 :]:
            x = scores[a].to_numpy(dtype=float)
            y = scores[b].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            n = int(mask.sum())
            if n < min_n or np.unique(x[mask]).size < 2 or np.unique(y[mask]).size < 2:
                r = p = np.nan
            else:
                r, p = spearmanr(x[mask], y[mask])
            rows.append({"a": a, "b": b, "n": n, "spearman_r": float(r) if r == r else np.nan, "p": float(p) if p == p else np.nan})
    return pd.DataFrame(rows)


def sign_concordance(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """For each headline pair, the Spearman r/p in every named cohort table.

    ``tables`` maps cohort name -> correlation_spearman.tsv-like frame with
    columns target, score, spearman_r, p, p_adj (p_adj optional).
    """
    rows = []
    for target, score in HEADLINE_PAIRS:
        rec: dict = {"target": target, "score": score}
        signs = []
        for cohort, frame in tables.items():
            hit = frame[(frame["target"] == target) & (frame["score"] == score)]
            if hit.empty:
                rec[f"{cohort}_r"] = np.nan
                rec[f"{cohort}_p"] = np.nan
                rec[f"{cohort}_fdr"] = np.nan
                continue
            r = float(hit.iloc[0]["spearman_r"])
            rec[f"{cohort}_r"] = r
            rec[f"{cohort}_p"] = float(hit.iloc[0]["p"])
            rec[f"{cohort}_fdr"] = float(hit.iloc[0]["p_adj"]) if "p_adj" in hit.columns else np.nan
            if np.isfinite(r) and r != 0:
                signs.append(np.sign(r))
        rec["n_cohorts_with_r"] = len(signs)
        rec["n_negative"] = int(sum(s < 0 for s in signs))
        rec["n_positive"] = int(sum(s > 0 for s in signs))
        rec["sign_agree"] = bool(len(signs) >= 2 and len(set(signs)) == 1)
        rows.append(rec)
    return pd.DataFrame(rows)
