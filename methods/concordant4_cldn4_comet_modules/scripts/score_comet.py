#!/usr/bin/env python3
"""COMET-style XL-mHG plus within-unit Q4 vs Q1 deltas for one cohort.

CLDN4-high / CLDN4-low are the top and bottom quartiles of malignant-cell
log1p(CP10k CLDN4), assigned inside each locked unit (ties.method = first).
A unit is used only when n_malignant >= 40, at least 10 malignant cells have
CLDN4 > 0, and the Q4 mean exceeds the Q1 mean.

The patient-level mean of (Q4 - Q1) is the consistency unit. Cell-level AUC
and the XL-mHG statistic describe the same contrast on pooled cells. Cell
counts are not an n.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import rankdata

# xlmhg 2.5.4 was built against NumPy 2; its test module still mentions np.float.
if not hasattr(np, "float"):
    np.float = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "bool"):
    np.bool = np.bool_  # type: ignore[attr-defined]
import xlmhg.mhg_cython as xlmhg_cython

MIN_MAL = 40
MIN_CLDN4_POS = 10
MIN_ARM = 10
TOL = np.longdouble(1e-12)


def load_cohort(path: Path):
    if (path / "counts.npz").exists():
        X = sparse.load_npz(path / "counts.npz").tocsr()
        genes = pd.read_csv(path / "genes.tsv", header=None)[0].astype(str).to_numpy()
        lib = np.load(path / "libsize.npy")
    elif (path / "counts.mtx").exists():
        X = mmread(path / "counts.mtx").tocsr()
        genes = pd.read_csv(path / "features.tsv", header=None)[0].astype(str).to_numpy()
        lib = pd.read_csv(path / "libsize.tsv", sep="\t")["libsize"].to_numpy(dtype=np.float64)
        genes_u = np.array([g.strip().upper() for g in genes], dtype=object)
        if len(set(genes_u.tolist())) != len(genes_u):
            uniq, inv = np.unique(genes_u, return_inverse=True)
            coo = X.tocoo()
            X = sparse.csr_matrix((coo.data, (inv[coo.row], coo.col)), shape=(len(uniq), X.shape[1]))
            X.sum_duplicates()
            genes = uniq
        else:
            genes = genes_u
        X = X.astype(np.float32)
    else:
        raise SystemExit(f"no counts in {path}")
    meta = pd.read_csv(path / "meta.tsv", sep="\t")
    if X.shape[1] != len(meta) or X.shape[1] != len(lib):
        raise SystemExit(f"shape mismatch {X.shape} meta {len(meta)} lib {len(lib)}")
    if len(genes) != X.shape[0]:
        raise SystemExit("gene length mismatch")
    return X.tocsr(), genes, meta, np.asarray(lib, dtype=np.float64)


def ranks_first(x: np.ndarray) -> np.ndarray:
    """1..n ranks. Ties keep the earlier cell at the smaller rank (R ties.method='first')."""
    order = np.lexsort((np.arange(x.size), x))
    ranks = np.empty(x.size, dtype=np.int32)
    ranks[order] = np.arange(1, x.size + 1)
    return ranks


def xlmhg_stat(expr: np.ndarray, positive: np.ndarray) -> tuple[float, int]:
    """XL-mHG statistic for `positive` cells enriched at high `expr`.

    Ties are broken against the positive class (positives sort later), so equal
    expression does not count as evidence of enrichment.
    """
    n = expr.size
    if n > 65535:
        raise RuntimeError(f"XL-mHG list longer than 65535 ({n})")
    order = np.lexsort((positive.astype(np.int8), -expr))
    ix = np.flatnonzero(positive[order]).astype(np.uint16)
    k = int(ix.size)
    if k == 0 or k == n:
        return 1.0, 0
    stat, cutoff = xlmhg_cython.get_xlmhg_stat(ix, n, k, 1, n, TOL)
    return float(stat), int(cutoff)


def score_cohort(cohort_dir: Path, out_dir: Path) -> None:
    cohort = cohort_dir.name
    print(f"==== score {cohort} ====", flush=True)
    X, genes, meta, lib = load_cohort(cohort_dir)
    gene_index = {g: i for i, g in enumerate(genes.tolist())}
    if "CLDN4" not in gene_index:
        raise SystemExit(f"{cohort} has no CLDN4")
    scale = (1e4 / np.maximum(lib, 1.0)).astype(np.float64)
    cldn4_counts = np.asarray(X[gene_index["CLDN4"]].todense()).ravel()
    cldn4_log = np.log1p(cldn4_counts * scale)

    arm = np.full(X.shape[1], "unused", dtype=object)
    unit_rows = []
    units = meta["unit_id"].astype(str).to_numpy()
    for uid in pd.unique(units):
        idx = np.flatnonzero(units == uid)
        n = int(idx.size)
        n_pos = int((cldn4_counts[idx] > 0).sum())
        rec = {
            "cohort": cohort,
            "unit_id": uid,
            "n_malignant": n,
            "n_cldn4_pos": n_pos,
            "used": False,
            "reason": "",
            "n_high": 0,
            "n_low": 0,
            "cldn4_delta": np.nan,
        }
        if n < MIN_MAL:
            rec["reason"] = "n_malignant<40"
            unit_rows.append(rec)
            continue
        if n_pos < MIN_CLDN4_POS:
            rec["reason"] = "cldn4_pos<10"
            unit_rows.append(rec)
            continue
        x = cldn4_log[idx]
        r = ranks_first(x)
        q1 = int(np.floor(n * 0.25))
        q4 = int(np.ceil(n * 0.75))
        high = r > q4
        low = r <= q1
        if int(high.sum()) < MIN_ARM or int(low.sum()) < MIN_ARM:
            rec["reason"] = "arm<10"
            unit_rows.append(rec)
            continue
        delta = float(x[high].mean() - x[low].mean())
        if not (delta > 0):
            rec["reason"] = "cldn4_q4_not_above_q1"
            unit_rows.append(rec)
            continue
        arm[idx[high]] = "high"
        arm[idx[low]] = "low"
        arm[idx[~(high | low)]] = "mid"
        rec.update(used=True, reason="ok", n_high=int(high.sum()), n_low=int(low.sum()), cldn4_delta=delta)
        unit_rows.append(rec)
    unit_df = pd.DataFrame(unit_rows)
    used_units = unit_df.loc[unit_df["used"], "unit_id"].astype(str).tolist()
    print(unit_df.groupby("reason").size().to_string(), flush=True)
    print(f"  units used {len(used_units)} high {(arm=='high').sum()} low {(arm=='low').sum()}", flush=True)
    if len(used_units) < 5:
        raise SystemExit(f"{cohort} has fewer than 5 usable units")

    high_m = arm == "high"
    low_m = arm == "low"
    keep = high_m | low_m
    # Column-scale then log1p on nonzeros. Library size is the full-gene UMI sum.
    logX = X.dot(sparse.diags(scale)).tocsr()
    logX.data = np.log1p(logX.data.astype(np.float64))

    n_genes = X.shape[0]
    n_units = len(used_units)
    deltas = np.zeros((n_genes, n_units), dtype=np.float32)
    for ui, uid in enumerate(used_units):
        in_unit = units == uid
        h = in_unit & high_m
        l = in_unit & low_m
        deltas[:, ui] = (
            np.asarray(logX[:, h].mean(axis=1)).ravel() - np.asarray(logX[:, l].mean(axis=1)).ravel()
        )
    patient_mean = deltas.mean(axis=1)
    sign_up = (deltas > 0).mean(axis=1)
    sign_down = (deltas < 0).mean(axis=1)

    logA = logX[:, keep].tocsr()
    is_high = high_m[keep]
    n_high = int(is_high.sum())
    n_low = int((~is_high).sum())
    expr = np.zeros(logA.shape[1], dtype=np.float64)
    indptr, indices, data = logA.indptr, logA.indices, logA.data
    auc = np.full(n_genes, 0.5, dtype=np.float64)
    cell_delta = np.zeros(n_genes, dtype=np.float64)
    pct_high = np.zeros(n_genes, dtype=np.float64)
    pct_low = np.zeros(n_genes, dtype=np.float64)
    stat_up = np.ones(n_genes, dtype=np.float64)
    stat_down = np.ones(n_genes, dtype=np.float64)
    cut_up = np.zeros(n_genes, dtype=np.int32)
    cut_down = np.zeros(n_genes, dtype=np.int32)
    # Detection uses raw counts on the same cells.
    cntA = X[:, keep].tocsr()
    print(f"  XL-mHG on {n_genes} genes x {logA.shape[1]} arm cells", flush=True)
    for i in range(n_genes):
        a, b = indptr[i], indptr[i + 1]
        ca, cb = cntA.indptr[i], cntA.indptr[i + 1]
        if cb > ca:
            cols = cntA.indices[ca:cb]
            pct_high[i] = np.count_nonzero(is_high[cols]) / n_high
            pct_low[i] = np.count_nonzero(~is_high[cols]) / n_low
        if b - a < 20:
            continue
        expr[:] = 0.0
        expr[indices[a:b]] = data[a:b]
        ranks = rankdata(expr, method="average")
        r_high = ranks[is_high].sum()
        u = r_high - n_high * (n_high + 1) / 2.0
        auc[i] = u / (n_high * n_low)
        cell_delta[i] = expr[is_high].mean() - expr[~is_high].mean()
        stat_up[i], cut_up[i] = xlmhg_stat(expr, is_high)
        stat_down[i], cut_down[i] = xlmhg_stat(expr, ~is_high)
        if (i + 1) % 2000 == 0:
            print(f"    {i+1}/{n_genes}", flush=True)

    cldn4_i = gene_index["CLDN4"]
    if not (patient_mean[cldn4_i] > 0 and sign_up[cldn4_i] == 1.0):
        raise SystemExit(
            f"CLDN4 positive control failed: mean delta {patient_mean[cldn4_i]} "
            f"sign_up {sign_up[cldn4_i]}"
        )
    print(
        f"  CLDN4 control patient_delta={patient_mean[cldn4_i]:.3f} "
        f"auc={auc[cldn4_i]:.3f} cell_delta={cell_delta[cldn4_i]:.3f}",
        flush=True,
    )

    stats = pd.DataFrame({
        "gene": genes,
        "cohort": cohort,
        "n_units": n_units,
        "n_high": n_high,
        "n_low": n_low,
        "patient_delta_mean": patient_mean,
        "sign_frac_up": sign_up,
        "sign_frac_down": sign_down,
        "cell_delta": cell_delta,
        "auc": auc,
        "pct_high": pct_high,
        "pct_low": pct_low,
        "xlmhg_stat_up": stat_up,
        "xlmhg_cutoff_up": cut_up,
        "xlmhg_stat_down": stat_down,
        "xlmhg_cutoff_down": cut_down,
    })
    out_dir.mkdir(parents=True, exist_ok=True)
    stats.to_csv(out_dir / f"{cohort}_gene_stats.tsv.gz", sep="\t", index=False, compression="gzip")
    unit_df.to_csv(out_dir / f"{cohort}_units_used.tsv", sep="\t", index=False)
    np.savez_compressed(
        out_dir / f"{cohort}_patient_delta.npz",
        deltas=deltas,
        genes=np.array(genes, dtype=object),
        units=np.array(used_units, dtype=object),
    )
    labels = pd.DataFrame({
        "barcode": meta["barcode"].astype(str),
        "unit_id": units,
        "arm": arm,
    })
    labels.to_csv(out_dir / f"{cohort}_arm_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = {
        "cohort": cohort,
        "n_units_used": int(n_units),
        "n_units_input": int(unit_df.shape[0]),
        "n_high_cells": n_high,
        "n_low_cells": n_low,
        "n_genes": int(n_genes),
        "cldn4_patient_delta": float(patient_mean[cldn4_i]),
        "cldn4_auc": float(auc[cldn4_i]),
        "cldn4_cell_delta": float(cell_delta[cldn4_i]),
        "cldn4_sign_frac_up": float(sign_up[cldn4_i]),
    }
    (out_dir / f"{cohort}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--malig", default="/tmp/concordant4_malig")
    ap.add_argument("--results", default=str(Path(__file__).resolve().parents[1] / "results" / "tables"))
    ap.add_argument("--cohort", default="all")
    args = ap.parse_args()
    malig = Path(args.malig)
    out = Path(args.results)
    cohorts = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"] if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        score_cohort(malig / c, out)


if __name__ == "__main__":
    main()
