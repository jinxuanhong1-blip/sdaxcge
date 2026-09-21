#!/usr/bin/env python3
"""Augur-style cell-type prioritization with the patient as the CV unit.

An L2 logistic regression predicts within-cohort CLDN4 Q4 vs Q1 from the
shared activity panel. Folds are grouped by patient/sample/donor. The
reported score is the AUC of the unit-mean out-of-fold probability, with
the cell-level AUC kept as the classical Augur number. CLDN4 is held out.
Malignant cells also drop the epithelial gate genes.

A random forest was the first classifier. Its out-of-fold scores sat below
0.5 under within-cohort label shuffles, so the reported model is the
logistic regression, whose null is centered near 0.5. The forest table is
kept as a sensitivity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

EXTRACT = Path("/tmp/geo_c4/augur_extract/augur_panel.npz")
META = Path("/tmp/geo_c4/augur_extract/augur_panel.tsv")
ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"

SEED = 1
N_REPEATS = 12
N_PERM = 99
N_PER_UNIT = 12
N_GENES = 40
HOLDOUT = {"CLDN4"}
MALIGNANT_EXTRA = {"CLDN4", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN3", "CLDN7", "TACSTD2"}
CELL_TYPES = ("T", "NK", "myeloid", "B", "malignant")


def log_cpm(expr: np.ndarray, lib: np.ndarray) -> np.ndarray:
    lib = np.clip(lib.astype(np.float64), 1.0, None)
    return np.log1p(expr.astype(np.float64) / lib[:, None] * 1e4)


def residualize(train: np.ndarray, test: np.ndarray, train_cohort: np.ndarray, test_cohort: np.ndarray):
    tr = train.copy()
    te = test.copy()
    for coh in np.unique(train_cohort):
        m = train_cohort == coh
        mu = tr[m].mean(axis=0)
        tr[m] -= mu
        te[test_cohort == coh] -= mu
    return tr, te


def top_genes(train: np.ndarray, k: int) -> np.ndarray:
    v = train.var(axis=0)
    order = np.argsort(-v)
    keep = [i for i in order if v[i] > 0][:k]
    return np.array(keep, dtype=int)


def unit_auc_from_oof(units: np.ndarray, y: np.ndarray, proba: np.ndarray) -> tuple[float, float]:
    # mean probability per unit; label is constant within unit
    rows = []
    for u in pd.unique(units):
        m = units == u
        if not np.isfinite(proba[m]).any():
            continue
        rows.append((float(np.nanmean(proba[m])), int(y[m][0])))
    if len(rows) < 6:
        return float("nan"), float("nan")
    p = np.array([r[0] for r in rows])
    yy = np.array([r[1] for r in rows])
    if len(set(yy.tolist())) < 2:
        return float("nan"), float("nan")
    cell_ok = np.isfinite(proba)
    cell = float(roc_auc_score(y[cell_ok], proba[cell_ok])) if cell_ok.sum() > 10 and len(set(y[cell_ok].tolist())) > 1 else float("nan")
    return float(roc_auc_score(yy, p)), cell


def one_auc(x, y, units, cohorts, rng, gene_idx) -> tuple[float, float]:
    x = x[:, gene_idx]
    # equal cells per unit
    chosen = []
    for u in pd.unique(units):
        idx = np.flatnonzero(units == u)
        if len(idx) < 8:
            continue
        take = idx if len(idx) <= N_PER_UNIT else rng.choice(idx, N_PER_UNIT, replace=False)
        chosen.append(take)
    if not chosen:
        return float("nan"), float("nan")
    sel = np.concatenate(chosen)
    x, y, units, cohorts = x[sel], y[sel], units[sel], cohorts[sel]
    n_splits = 5
    if len(np.unique(units)) < n_splits or min(np.bincount(y)) < n_splits:
        n_splits = int(min(np.bincount(y)))
    if n_splits < 2:
        return float("nan"), float("nan")
    try:
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=int(rng.integers(0, 1_000_000)))
        splits = list(cv.split(x, y, groups=units))
    except ValueError:
        return float("nan"), float("nan")
    proba = np.full(len(y), np.nan)
    for train, test in splits:
        if len(set(y[train].tolist())) < 2 or len(set(y[test].tolist())) < 2:
            continue
        tr, te = residualize(x[train], x[test], cohorts[train], cohorts[test])
        keep = top_genes(tr, N_GENES)
        if len(keep) < 5:
            continue
        clf = LogisticRegression(
            C=0.2,
            max_iter=500,
            solver="lbfgs",
        )
        clf.fit(tr[:, keep], y[train])
        proba[test] = clf.predict_proba(te[:, keep])[:, 1]
    return unit_auc_from_oof(units, y, proba)


def run_type(df: pd.DataFrame, x: np.ndarray, genes: list[str], cell_type: str, rng: np.random.Generator) -> dict:
    sub = df[df.cell_type == cell_type]
    if cell_type == "malignant":
        hold = MALIGNANT_EXTRA
        hold_name = "CLDN4_and_epithelial_gate"
    else:
        hold = HOLDOUT
        hold_name = "CLDN4"
    gene_idx = np.array([i for i, g in enumerate(genes) if g not in hold], dtype=int)
    y = (sub.quartile.to_numpy() == "Q4").astype(int)
    units = sub.unit_id.to_numpy()
    cohorts = sub.dataset.to_numpy()
    # require both classes and enough units
    n_q1 = len(set(units[y == 0].tolist()))
    n_q4 = len(set(units[y == 1].tolist()))
    base = {
        "cell_type": cell_type,
        "holdout": hold_name,
        "n_cells": int(len(sub)),
        "n_units_q1": n_q1,
        "n_units_q4": n_q4,
        "datasets": ",".join(sorted(set(sub.dataset.tolist()))),
    }
    if n_q1 < 6 or n_q4 < 6:
        base.update({"unit_auc": float("nan"), "cell_auc": float("nan"), "p_perm": float("nan"), "note": "too few units"})
        return base
    xx = x[sub.index.to_numpy()]
    unit_scores, cell_scores = [], []
    for _ in range(N_REPEATS):
        ua, ca = one_auc(xx, y, units, cohorts, rng, gene_idx)
        unit_scores.append(ua)
        cell_scores.append(ca)
    unit_auc = float(np.nanmean(unit_scores))
    cell_auc = float(np.nanmean(cell_scores))
    # null: shuffle quartile by unit within cohort, one subsample repeated
    null = np.empty(N_PERM)
    # map unit -> cohort
    unit_cohort = {}
    unit_label = {}
    for u, c, lab in zip(units, cohorts, y):
        unit_cohort[u] = c
        unit_label[u] = int(lab)
    cohort_units = {}
    for u, c in unit_cohort.items():
        cohort_units.setdefault(c, []).append(u)
    for p_i in range(N_PERM):
        new_label = dict(unit_label)
        for c, us in cohort_units.items():
            labs = [unit_label[u] for u in us]
            shuffled = rng.permutation(labs)
            for u, lab in zip(us, shuffled):
                new_label[u] = int(lab)
        y_perm = np.array([new_label[u] for u in units])
        if len(set(y_perm.tolist())) < 2:
            null[p_i] = 0.5
            continue
        ua, _ = one_auc(xx, y_perm, units, cohorts, rng, gene_idx)
        null[p_i] = ua
    valid = null[np.isfinite(null)]
    if np.isfinite(unit_auc) and len(valid):
        # Augur priority is one-sided: higher AUC than the label-shuffle null.
        p_perm = float((1 + np.sum(valid >= unit_auc)) / (1 + len(valid)))
        p_two = float((1 + np.sum(np.abs(valid - 0.5) >= abs(unit_auc - 0.5))) / (1 + len(valid)))
    else:
        p_perm = float("nan")
        p_two = float("nan")
    base.update(
        {
            "unit_auc": unit_auc,
            "cell_auc": cell_auc,
            "unit_auc_sd": float(np.nanstd(unit_scores, ddof=1)),
            "p_perm": p_perm,
            "p_two_sided": p_two,
            "null_mean": float(np.nanmean(valid)) if len(valid) else float("nan"),
            "n_repeats": N_REPEATS,
            "n_perm": int(len(valid)),
            "note": "",
        }
    )
    return base


def main() -> None:
    if not EXTRACT.exists():
        raise SystemExit(f"missing {EXTRACT}")
    blob = np.load(EXTRACT)
    expr = blob["expr"]
    genes = [g.upper() for g in blob["genes"].tolist()]
    meta = pd.read_csv(META, sep="\t")
    if len(meta) != expr.shape[0]:
        raise SystemExit("meta/expr length mismatch")
    x = log_cpm(expr, meta.lib_size.to_numpy())
    meta = meta.reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    rows = []
    for ct in CELL_TYPES:
        print(f"Augur {ct}", flush=True)
        rec = run_type(meta, x, genes, ct, rng)
        print(
            f"  units Q1/Q4 {rec['n_units_q1']}/{rec['n_units_q4']} "
            f"unit_auc={rec.get('unit_auc')} cell_auc={rec.get('cell_auc')} p={rec.get('p_perm')}",
            flush=True,
        )
        rows.append(rec)
    # BH on the finite permutation p-values
    pvals = [r.get("p_perm", float("nan")) for r in rows]
    order = np.argsort([p if np.isfinite(p) else 2 for p in pvals])
    finite_idx = [i for i in order if np.isfinite(pvals[i])]
    m = len(finite_idx)
    q = {i: float("nan") for i in range(len(rows))}
    prev = 1.0
    for rank_from_end, i in enumerate(finite_idx[::-1]):
        orig_rank = m - rank_from_end
        val = pvals[i] * m / orig_rank
        prev = min(prev, val)
        q[i] = float(min(prev, 1.0))
    for i, r in enumerate(rows):
        r["q_bh"] = q[i]
    out = TAB / "augur_priority.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    datasets = sorted(meta.dataset.unique())
    summary = {
        "datasets": datasets,
        "n_cells": int(len(meta)),
        "n_units": int(meta.groupby(["dataset", "unit_id"]).ngroups),
        "n_genes_panel": len(genes),
        "primary_score": "unit-mean OOF probability AUC; folds grouped by unit",
        "cldn4_held_out": True,
    }
    (TAB / "augur_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
