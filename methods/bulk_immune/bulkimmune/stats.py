"""Association of TACSTD2 / CLDN4 with immune scores, with an honest multiple-testing plan.

The analysis this module is built for:

    For each target gene (TACSTD2, CLDN4) and each immune score, report a
    rank correlation, a partial rank correlation after purity/batch, and
    (when a clinical endpoint exists) a two-group or survival test of the
    score itself.

Multiple testing is *not* "BH on the whole table". The table mixes three
families that were not equally likely a priori:

    A.  Exclusion / CAF / fibroblast / TIDE Exclusion -- the hypothesis that
        TACSTD2/CLDN4 mark an excluded epithelium.
    B.  Effector / IFN / TLS / CD8 -- the competing "hot tumour" hypothesis.
    C.  Everything else (the rest of CIBERSORT, xCell, MCP-counter, TIP steps).

Family-wise BH-FDR is applied inside each family. A result in family C that
does not survive its own FDR is a hypothesis, not a finding. Cross-family
comparisons ("Exclusion FDR = 0.01, B-cell FDR = 0.04, therefore Exclusion
wins") are not licensed by this procedure.

Spearman is the default. Immune fractions are zero-inflated and signatures
are skewed; Pearson on raw scores is almost always the wrong first look.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

__all__ = [
    "spearman_table",
    "partial_spearman",
    "two_group",
    "cox_ph",
    "adjust_pvalues",
    "DEFAULT_FAMILIES",
]


DEFAULT_FAMILIES = {
    "exclusion": [
        "TIDE", "Exclusion", "Dysfunction", "CAF", "MDSC", "TAM M2",
        "Fibroblasts", "Endothelial cells", "StromalScore", "StromaScore",
        "Step5_infiltration_of_immune_cells",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_HYPOXIA",
        "HALLMARK_ANGIOGENESIS",
    ],
    "inflamed": [
        "IFNG", "CD8", "CTL", "ImmuneScore", "CYT",
        "IFNG_Ayers6", "IFNG_Ayers_preliminary10", "ExpandedImmune_Ayers18",
        "TcellInflamed_GEP18", "CD8_Teffector_Mariathasan",
        "TLS_12chemokine_Coppola", "TLS_Cabrita9", "TLS_Cabrita9_noY",
        "Bcell_follicular",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_ALLOGRAFT_REJECTION",
        "CD8 T cells", "T cells", "Cytotoxic lymphocytes", "B lineage",
        "CD8+ T-cells", "CD4+ T-cells", "B-cells", "NK cells",
        "T.cells.CD8", "T.cells.CD4", "NK.cells", "B.cells",
        "Step4_trafficking_of_immune_cells",
        "Step4_CD8_T_cell", "Step4_T_cell", "Step4_NK_cell", "Step4_B_cell",
        "Step7_killing_of_cancer_cells",
        "Cytolytic_CYT_Rooney",
    ],
}


def _safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 5:
        return np.nan, np.nan, n
    if np.unique(x[mask]).size < 2 or np.unique(y[mask]).size < 2:
        return np.nan, np.nan, n
    r, p = spearmanr(x[mask], y[mask])
    return float(r), float(p), n


def spearman_table(
    targets: pd.DataFrame,
    scores: pd.DataFrame,
    min_n: int = 8,
) -> pd.DataFrame:
    """Spearman correlation of every target column with every score column."""
    common = targets.index.intersection(scores.index)
    T = targets.loc[common]
    S = scores.loc[common]
    rows = []
    for tcol in T.columns:
        x = T[tcol].to_numpy(dtype=float)
        for scol in S.columns:
            r, p, n = _safe_spearman(x, S[scol].to_numpy(dtype=float))
            if n < min_n:
                r, p = np.nan, np.nan
            rows.append(
                {
                    "target": tcol,
                    "score": scol,
                    "n": n,
                    "spearman_r": r,
                    "p": p,
                }
            )
    return pd.DataFrame(rows)


def partial_spearman(
    x: pd.Series,
    y: pd.Series,
    covariates: pd.DataFrame,
) -> dict[str, float]:
    """Spearman correlation of OLS residuals of rank(x) and rank(y) on covariates.

    This is the rank analogue of a partial Pearson correlation. It answers
    "do TACSTD2 and ImmuneScore still co-vary after linear removal of purity
    and batch", which is the question the playbook asks. It is *not* a
    model-based causal estimate.
    """
    frame = pd.concat(
        [x.rename("_x"), y.rename("_y"), covariates],
        axis=1,
    ).dropna()
    if frame.shape[0] < 8:
        return {"n": float(frame.shape[0]), "partial_r": np.nan, "p": np.nan}

    def _rank_resid(col: str) -> np.ndarray:
        rnk = frame[col].rank(method="average").to_numpy(dtype=float)
        cov = frame.drop(columns=["_x", "_y"])
        parts = []
        for c in cov.columns:
            s = cov[c]
            if pd.api.types.is_numeric_dtype(s):
                parts.append(s.astype(float).to_numpy()[:, None])
            else:
                parts.append(pd.get_dummies(s, drop_first=True).to_numpy(dtype=float))
        A = np.column_stack([np.ones(len(frame))] + parts) if parts else np.ones((len(frame), 1))
        beta, *_ = np.linalg.lstsq(A, rnk, rcond=None)
        return rnk - A @ beta

    rx, ry = _rank_resid("_x"), _rank_resid("_y")
    r, p, n = _safe_spearman(rx, ry)
    return {"n": float(n), "partial_r": r, "p": p}


def two_group(
    scores: pd.DataFrame,
    group: pd.Series,
    positive: str | None = None,
) -> pd.DataFrame:
    """Mann-Whitney U of each score between two groups (e.g. ICI responder)."""
    group = group.loc[scores.index]
    levels = list(pd.unique(group.dropna()))
    if positive is not None:
        a, b = positive, [lv for lv in levels if lv != positive][0]
    else:
        a, b = levels[0], levels[1]
    rows = []
    for col in scores.columns:
        xa = scores.loc[group == a, col].dropna().to_numpy()
        xb = scores.loc[group == b, col].dropna().to_numpy()
        if xa.size < 3 or xb.size < 3:
            stat, p, delta = np.nan, np.nan, np.nan
        else:
            stat, p = mannwhitneyu(xa, xb, alternative="two-sided")
            # Cliff's delta: (P(a>b) - P(a<b))
            # U = n_a n_b P(a>b) + 0.5 n_a n_b P(tie)
            delta = (2.0 * stat) / (xa.size * xb.size) - 1.0
        rows.append(
            {
                "score": col,
                "group_a": a,
                "n_a": int(xa.size),
                "median_a": float(np.median(xa)) if xa.size else np.nan,
                "group_b": b,
                "n_b": int(xb.size),
                "median_b": float(np.median(xb)) if xb.size else np.nan,
                "U": stat,
                "cliffs_delta": delta,
                "p": p,
            }
        )
    return pd.DataFrame(rows)


def cox_ph(
    scores: pd.DataFrame,
    time: pd.Series,
    event: pd.Series,
) -> pd.DataFrame:
    """Univariable Cox PH per score (Breslow ties), via a log-rank-style score test fallback.

    Uses ``lifelines`` when installed; otherwise reports a log-rank p-value on
    a median split so the demo still produces a real survival statistic
    without a heavy dependency. The method used is recorded in the table.
    """
    time = time.loc[scores.index].astype(float)
    event = event.loc[scores.index].astype(int)
    rows = []
    try:
        from lifelines import CoxPHFitter

        have_lifelines = True
    except ImportError:
        have_lifelines = False

    from scipy.stats import chi2

    for col in scores.columns:
        frame = pd.DataFrame(
            {"score": scores[col], "time": time, "event": event}
        ).dropna()
        if frame.shape[0] < 10 or frame["event"].sum() < 3:
            rows.append(
                {"score": col, "n": frame.shape[0], "events": int(frame["event"].sum()),
                 "hr": np.nan, "p": np.nan, "method": "skipped"}
            )
            continue
        if have_lifelines:
            cph = CoxPHFitter()
            try:
                cph.fit(frame, duration_col="time", event_col="event")
                rows.append(
                    {
                        "score": col,
                        "n": frame.shape[0],
                        "events": int(frame["event"].sum()),
                        "hr": float(np.exp(cph.params_["score"])),
                        "p": float(cph.summary.loc["score", "p"]),
                        "method": "coxph_lifelines",
                    }
                )
                continue
            except Exception:
                pass
        # Median-split log-rank (Mantel-Haenszel)
        hi = frame["score"] >= frame["score"].median()
        p = _logrank_p(frame.loc[hi, "time"], frame.loc[hi, "event"],
                       frame.loc[~hi, "time"], frame.loc[~hi, "event"])
        rows.append(
            {
                "score": col,
                "n": frame.shape[0],
                "events": int(frame["event"].sum()),
                "hr": np.nan,
                "p": p,
                "method": "logrank_median_split",
            }
        )
    return pd.DataFrame(rows)


def _logrank_p(t1, e1, t2, e2) -> float:
    from scipy.stats import chi2

    t1, e1 = np.asarray(t1, float), np.asarray(e1, int)
    t2, e2 = np.asarray(t2, float), np.asarray(e2, int)
    times = np.unique(np.concatenate([t1, t2]))
    o1 = o2 = e_exp1 = var = 0.0
    n1, n2 = float(t1.size), float(t2.size)
    i1 = i2 = 0
    order1, order2 = np.argsort(t1), np.argsort(t2)
    t1, e1 = t1[order1], e1[order1]
    t2, e2 = t2[order2], e2[order2]
    for t in times:
        d1 = 0.0
        while i1 < t1.size and t1[i1] == t:
            d1 += e1[i1]
            i1 += 1
        d2 = 0.0
        while i2 < t2.size and t2[i2] == t:
            d2 += e2[i2]
            i2 += 1
        d = d1 + d2
        n = n1 + n2
        if n > 1 and d > 0:
            e_exp1 += d * (n1 / n)
            var += d * (n - d) * n1 * n2 / (n ** 2 * (n - 1))
            o1 += d1
            o2 += d2
        n1 -= float(np.sum(t1 == t))
        n2 -= float(np.sum(t2 == t))
    if var <= 0:
        return np.nan
    stat = (o1 - e_exp1) ** 2 / var
    return float(chi2.sf(stat, 1))


def assign_family(score: str, families: dict[str, Sequence[str]] | None = None) -> str:
    families = families or DEFAULT_FAMILIES
    for name, members in families.items():
        if score in members:
            return name
        # prefix match for TIP cell-type scores and CIBERSORT LM22 names
        for m in members:
            if score == m or score.startswith(m + ".") or score.endswith(m):
                return name
    return "other"


def adjust_pvalues(
    table: pd.DataFrame,
    p_col: str = "p",
    group_cols: Sequence[str] = ("target",),
    family_col: str | None = "family",
    method: str = "fdr_bh",
) -> pd.DataFrame:
    """BH-FDR within (target x family), leaving NaN p-values untouched."""
    out = table.copy()
    if family_col and family_col not in out.columns:
        out[family_col] = out["score"].map(lambda s: assign_family(str(s)))
    out["p_adj"] = np.nan
    grouper = list(group_cols) + ([family_col] if family_col else [])
    if not grouper:
        mask = out[p_col].notna()
        if mask.any():
            _, padj, _, _ = multipletests(out.loc[mask, p_col].to_numpy(), method=method)
            out.loc[mask, "p_adj"] = padj
            out.loc[mask, "rejected_0.05"] = padj < 0.05
        return out
    for _, idx in out.groupby(grouper, dropna=False).groups.items():
        sub = out.loc[idx, p_col]
        mask = sub.notna()
        if mask.sum() == 0:
            continue
        rejected, padj, _, _ = multipletests(sub[mask].to_numpy(), method=method)
        out.loc[sub.index[mask], "p_adj"] = padj
        out.loc[sub.index[mask], "rejected_0.05"] = rejected
    return out
