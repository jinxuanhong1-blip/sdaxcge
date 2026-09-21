#!/usr/bin/env python3
"""Exhaustive direction sweep on GSE207422.

Thesis to hunt: CLDN4-low cells have higher IFN, higher STING, and lower NHEJ.
Delta is always CLDN4-high minus CLDN4-low (or a Spearman rho of CLDN4 vs the
module). Thesis signs are IFN delta < 0, STING delta < 0, NHEJ delta > 0.

A specification is a hit only if all three signs match and the contrast is
reliable:
- paired within-patient deltas: n >= 4 patients, and each module has a strict
  majority of patients on the thesis side of zero
- one correlation across samples: n >= 8 and |rho| >= 0.3 on every module
Pooled cell-level tests are recorded and cannot be a hit.

This file lists the specifications before any sweep result is interpreted.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import assign_lineage, log1p_cp10k, module_vector
from gene_sets import A3_NORMAL_LUNG, MODULES

HERE = Path(__file__).resolve().parents[1]
MIN_CELLS = 20
MIN_TAIL = 8
MIN_MATCH = 8
MIN_N_PAIRED = 4
MIN_N_CORR = 8
MIN_ABS_RHO = 0.3

def _is_cyto_ribo(gene: str) -> bool:
    """Cytosolic ribosomal-protein genes. Pseudogenes (RPL7P*, RPS3AP*) are out."""
    if gene in {"RPSA", "RPLP0", "RPLP1", "RPLP2"}:
        return True
    if gene.startswith("RPS") or gene.startswith("RPL"):
        rest = gene[3:]
        if rest[:1].isdigit() and "P" not in rest:
            return True
    return False


MODULES_3 = ("ifn", "sting", "nhej")
SCORE_SRC = {"ifn": "ifn_effector", "sting": "sting", "nhej": "nhej"}


def stream_aux(matrix: Path, cells_expected: np.ndarray, cache: Path) -> dict[str, np.ndarray]:
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if len(z["cells"]) == len(cells_expected) and np.array_equal(z["cells"], cells_expected):
            print(f"aux cache {cache}", flush=True)
            return {k: z[k] for k in z.files}
    subtype_genes = ["KRT5", "KRT6A", "KRT7", "KRT17", "TP63", "SOX2"]
    want = set(subtype_genes)
    n = len(cells_expected)
    ribo = np.zeros(n, dtype=np.int64)
    mito = np.zeros(n, dtype=np.int64)
    n_ribo_genes = 0
    n_mito_genes = 0
    rows: dict[str, np.ndarray] = {}
    with gzip.open(matrix, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        if not np.array_equal(cells, cells_expected):
            raise SystemExit("UMI cell order does not match the extracted matrix")
        for line in fh:
            i = line.find("\t")
            gene = line[:i]
            is_ribo = _is_cyto_ribo(gene)
            is_mito = gene.startswith("MT-")
            keep = gene in want
            if not (is_ribo or is_mito or keep):
                continue
            arr = np.fromstring(line[i + 1 :], sep="\t", dtype=np.int32)
            if arr.size != n:
                raise SystemExit(f"{gene}: {arr.size} values, expected {n}")
            if is_ribo:
                ribo += arr.astype(np.int64)
                n_ribo_genes += 1
            elif is_mito:
                mito += arr.astype(np.int64)
                n_mito_genes += 1
            else:
                rows[gene] = arr
    payload = {
        "cells": cells,
        "ribo": ribo,
        "mito": mito,
        "n_ribo_genes": np.array([n_ribo_genes]),
        "n_mito_genes": np.array([n_mito_genes]),
    }
    for g in subtype_genes:
        payload[g] = rows.get(g, np.zeros(n, dtype=np.int32))
        payload[f"{g}_present"] = np.array([g in rows])
    np.savez_compressed(cache, **payload)
    print(
        f"aux streamed ribo_genes={n_ribo_genes} mito_genes={n_mito_genes} "
        f"subtype_present={[g for g in subtype_genes if g in rows]}",
        flush=True,
    )
    return {k: payload[k] for k in payload}


def majority_thesis(deltas: np.ndarray, thesis_negative: bool) -> tuple[bool, int, int, int]:
    d = deltas[np.isfinite(deltas)]
    n_neg = int(np.sum(d < 0))
    n_pos = int(np.sum(d > 0))
    n_tie = int(np.sum(d == 0))
    if len(d) == 0 or np.median(d) == 0:
        return False, n_neg, n_pos, n_tie
    if thesis_negative:
        ok = (np.median(d) < 0) and (n_neg > n_pos)
    else:
        ok = (np.median(d) > 0) and (n_pos > n_neg)
    return ok, n_neg, n_pos, n_tie


def judge_paired(spec: str, family: str, patient_df: pd.DataFrame) -> dict:
    row = {
        "spec": spec,
        "family": family,
        "kind": "paired",
        "n": int(len(patient_df)),
        "hit": False,
        "joint_sign": False,
        "reliable": False,
    }
    for name, thesis_neg in (("ifn", True), ("sting", True), ("nhej", False)):
        if patient_df.empty or name not in patient_df.columns:
            row[f"{name}_median"] = np.nan
            row[f"{name}_n_thesis"] = 0
            row[f"{name}_n_opposite"] = 0
            row[f"{name}_sign_ok"] = False
            continue
        d = patient_df[name].to_numpy(dtype=float)
        ok, n_neg, n_pos, _n_tie = majority_thesis(d, thesis_neg)
        row[f"{name}_median"] = float(np.median(d[np.isfinite(d)])) if np.isfinite(d).any() else np.nan
        row[f"{name}_n_thesis"] = n_neg if thesis_neg else n_pos
        row[f"{name}_n_opposite"] = n_pos if thesis_neg else n_neg
        row[f"{name}_sign_ok"] = bool(ok)
    row["joint_sign"] = bool(row["ifn_sign_ok"] and row["sting_sign_ok"] and row["nhej_sign_ok"])
    row["reliable"] = bool(row["n"] >= MIN_N_PAIRED and row["joint_sign"])
    row["hit"] = bool(row["reliable"])
    return row


def judge_corr(spec: str, family: str, rho: dict[str, float], n: int) -> dict:
    row = {"spec": spec, "family": family, "kind": "correlation", "n": int(n), "hit": False, "joint_sign": False, "reliable": False}
    signs = []
    for name, thesis_neg in (("ifn", True), ("sting", True), ("nhej", False)):
        r = rho.get(name, np.nan)
        row[f"{name}_median"] = r
        row[f"{name}_n_thesis"] = np.nan
        row[f"{name}_n_opposite"] = np.nan
        if not np.isfinite(r) or r == 0:
            row[f"{name}_sign_ok"] = False
            signs.append(False)
            continue
        sign_ok = (r < 0) if thesis_neg else (r > 0)
        mag_ok = abs(r) >= MIN_ABS_RHO
        row[f"{name}_sign_ok"] = bool(sign_ok and mag_ok)
        signs.append(sign_ok)
    row["joint_sign"] = bool(all(signs) and all(np.isfinite([rho.get(k, np.nan) for k in MODULES_3])))
    row["reliable"] = bool(n >= MIN_N_CORR and row["ifn_sign_ok"] and row["sting_sign_ok"] and row["nhej_sign_ok"])
    row["hit"] = bool(row["reliable"])
    return row


def split_masks(values: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray] | None:
    v = np.asarray(values, dtype=float)
    if mode == "detected":
        hi, lo = v > 0, v == 0
    elif mode == "q4q1":
        q75, q25 = np.quantile(v, 0.75), np.quantile(v, 0.25)
        if not np.isfinite(q75) or q75 <= q25:
            return None
        hi, lo = v >= q75, v <= q25
    elif mode == "median":
        med = np.median(v)
        if np.all(v == med):
            return None
        hi, lo = v > med, v <= med
    elif mode == "tertile":
        hi_c, lo_c = np.quantile(v, 2 / 3), np.quantile(v, 1 / 3)
        if hi_c <= lo_c:
            return None
        hi, lo = v >= hi_c, v <= lo_c
    elif mode == "decile":
        hi_c, lo_c = np.quantile(v, 0.90), np.quantile(v, 0.10)
        if hi_c <= lo_c:
            return None
        hi, lo = v >= hi_c, v <= lo_c
    elif mode == "top20":
        hi_c, lo_c = np.quantile(v, 0.80), np.quantile(v, 0.20)
        if hi_c <= lo_c:
            return None
        hi, lo = v >= hi_c, v <= lo_c
    elif mode == "pos_median":
        m = v > 0
        if m.sum() < MIN_CELLS:
            return None
        med = np.median(v[m])
        hi, lo = (v > med) & m, (v <= med) & m & (v > 0)
        if med <= 0:
            return None
    elif mode == "pos_q4q1":
        m = v > 0
        if m.sum() < MIN_CELLS:
            return None
        q75, q25 = np.quantile(v[m], 0.75), np.quantile(v[m], 0.25)
        if q75 <= q25:
            return None
        hi, lo = (v >= q75) & m, (v <= q25) & m
    else:
        raise ValueError(mode)
    if hi.sum() < MIN_TAIL or lo.sum() < MIN_TAIL:
        return None
    if np.any(hi & lo):
        return None
    return hi, lo


def patient_split_deltas(
    df: pd.DataFrame,
    mask: np.ndarray,
    splitter: str,
    score_cols: dict[str, str],
    mode: str,
    min_cells: int = MIN_CELLS,
) -> pd.DataFrame:
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < min_cells:
            continue
        parts = split_masks(g[splitter].to_numpy(), mode)
        if parts is None:
            continue
        hi, lo = parts
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": int(hi.sum()), "n_low": int(lo.sum())}
        for out_name, col in score_cols.items():
            rec[out_name] = float(g[col].to_numpy()[hi].mean() - g[col].to_numpy()[lo].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_double_low(df: pd.DataFrame, mask: np.ndarray, score_cols: dict[str, str]) -> pd.DataFrame:
    """Delta = double-high minus double-low. Thesis: double-low has higher IFN/STING."""
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        c4 = g["CLDN4"].to_numpy()
        tac = g["TACSTD2"].to_numpy()
        cmed, tmed = np.median(c4), np.median(tac)
        if np.all(c4 == cmed) or np.all(tac == tmed):
            continue
        hi = (c4 > cmed) & (tac > tmed)
        lo = (c4 <= cmed) & (tac <= tmed)
        if hi.sum() < MIN_TAIL or lo.sum() < MIN_TAIL:
            continue
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": int(hi.sum()), "n_low": int(lo.sum())}
        for out_name, col in score_cols.items():
            arr = g[col].to_numpy()
            rec[out_name] = float(arr[hi].mean() - arr[lo].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_caliper(df: pd.DataFrame, mask: np.ndarray, score_cols: dict[str, str], caliper: float) -> pd.DataFrame:
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        parts = split_masks(g["CLDN4"].to_numpy(), "q4q1")
        if parts is None:
            continue
        hi, lo = parts
        umi = g["total_umi"].to_numpy(dtype=float)
        hi_i = np.flatnonzero(hi)
        lo_i = np.flatnonzero(lo)
        lh = np.log1p(umi[hi_i])
        ll = np.log1p(umi[lo_i])
        used = np.zeros(len(lo_i), dtype=bool)
        keep_h: list[int] = []
        keep_l: list[int] = []
        for i in np.argsort(lh):
            d = np.abs(ll - lh[i])
            d[used] = np.inf
            j = int(np.argmin(d))
            if not np.isfinite(d[j]) or d[j] > caliper:
                continue
            used[j] = True
            keep_h.append(int(hi_i[i]))
            keep_l.append(int(lo_i[j]))
        if len(keep_h) < MIN_MATCH:
            continue
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": len(keep_h), "n_low": len(keep_l)}
        for out_name, col in score_cols.items():
            arr = g[col].to_numpy()
            rec[out_name] = float(arr[keep_h].mean() - arr[keep_l].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_umi_bins(df: pd.DataFrame, mask: np.ndarray, score_cols: dict[str, str], n_bins: int) -> pd.DataFrame:
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        umi = g["total_umi"].to_numpy(dtype=float)
        cldn = g["CLDN4"].to_numpy(dtype=float)
        try:
            bins = np.asarray(pd.qcut(np.asarray(umi), n_bins, duplicates="drop", labels=False))
        except ValueError:
            continue
        acc = {k: [] for k in score_cols}
        for b in np.unique(bins[~pd.isna(bins)]):
            m = bins == b
            if m.sum() < 2 * MIN_TAIL:
                continue
            med = np.median(cldn[m])
            if np.all(cldn[m] == med):
                continue
            hi = m & (cldn > med)
            lo = m & (cldn <= med)
            if hi.sum() < MIN_TAIL or lo.sum() < MIN_TAIL:
                continue
            for out_name, col in score_cols.items():
                arr = g[col].to_numpy()
                acc[out_name].append(float(arr[hi].mean() - arr[lo].mean()))
        if any(len(v) < 2 for v in acc.values()):
            continue
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": int(sum(len(v) for v in acc.values()) / 3), "n_low": int(sum(len(v) for v in acc.values()) / 3)}
        for out_name, vals in acc.items():
            rec[out_name] = float(np.mean(vals))
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_resid_split(
    df: pd.DataFrame,
    mask: np.ndarray,
    covariates: list[str],
    score_cols: dict[str, str],
) -> pd.DataFrame:
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        parts = split_masks(g["CLDN4"].to_numpy(), "q4q1")
        if parts is None:
            continue
        hi, lo = parts
        X = np.column_stack([g[c].to_numpy(dtype=float) for c in covariates])
        if not np.isfinite(X).all():
            continue
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": int(hi.sum()), "n_low": int(lo.sum())}
        A = np.column_stack([np.ones(len(g)), X])
        for out_name, col in score_cols.items():
            y = g[col].to_numpy(dtype=float)
            if not np.isfinite(y).all():
                rec[out_name] = np.nan
                continue
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            resid = y - A @ coef
            rec[out_name] = float(resid[hi].mean() - resid[lo].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_spearman(df: pd.DataFrame, mask: np.ndarray, score_cols: dict[str, str]) -> pd.DataFrame:
    rows = []
    sub = df.loc[mask]
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        x = g["CLDN4"].to_numpy(dtype=float)
        if np.unique(x).size < 2:
            continue
        rec = {"Sample": sid, "n_cells": int(len(g)), "n_high": int(len(g)), "n_low": int(len(g))}
        ok = True
        for out_name, col in score_cols.items():
            y = g[col].to_numpy(dtype=float)
            if np.unique(y).size < 2:
                ok = False
                break
            rec[out_name] = float(stats.spearmanr(x, y).statistic)
        if ok:
            rows.append(rec)
    return pd.DataFrame(rows)


def mean_spearman(df: pd.DataFrame, mask: np.ndarray) -> tuple[dict[str, float], int]:
    """One rho per module from patient-mean CLDN4 vs patient-mean module."""
    sub = df.loc[mask]
    rows = []
    for sid, g in sub.groupby("Sample", sort=True):
        if len(g) < MIN_CELLS:
            continue
        rows.append({
            "CLDN4": float(g["CLDN4"].mean()),
            "ifn": float(g["ifn"].mean()),
            "sting": float(g["sting"].mean()),
            "nhej": float(g["nhej"].mean()),
        })
    tab = pd.DataFrame(rows)
    if len(tab) < 4:
        return {k: np.nan for k in MODULES_3}, int(len(tab))
    out = {}
    for name in MODULES_3:
        if tab["CLDN4"].nunique() < 2 or tab[name].nunique() < 2:
            out[name] = np.nan
        else:
            out[name] = float(stats.spearmanr(tab["CLDN4"], tab[name]).statistic)
    return out, int(len(tab))


def bulk_partial(x, y, z) -> float:
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    if len(x) < MIN_N_CORR or np.unique(z).size < 2 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan

    def resid(a, cov):
        A = np.column_stack([np.ones(len(cov)), stats.rankdata(cov)])
        coef, *_ = np.linalg.lstsq(A, stats.rankdata(a), rcond=None)
        return stats.rankdata(a) - A @ coef

    rx, ry = resid(x, z), resid(y, z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(stats.pearsonr(rx, ry)[0])


def add_paired(results, long_rows, spec, family, patient_df):
    results.append(judge_paired(spec, family, patient_df))
    if len(patient_df):
        for _, r in patient_df.iterrows():
            long_rows.append({
                "spec": spec,
                "Sample": r["Sample"],
                "ifn": r.get("ifn", np.nan),
                "sting": r.get("sting", np.nan),
                "nhej": r.get("nhej", np.nan),
                "n_high": r.get("n_high", np.nan),
                "n_low": r.get("n_low", np.nan),
            })


def main() -> None:
    ap_workdir = Path("data/GSE207422")
    matrix = ap_workdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    extracted_path = ap_workdir / "extracted_nhej_sting_ifn.npz"
    cache = ap_workdir / "aux_ribo_subtype_v2.npz"
    extracted = np.load(extracted_path, allow_pickle=True)
    cells = extracted["cells"]
    aux = stream_aux(matrix, cells, cache)
    genes = [str(g) for g in extracted["genes"]]
    mat = extracted["mat"]
    total = extracted["total"].astype(np.float64)
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)
    ribo = aux["ribo"].astype(np.float64)
    nonribo_lib = np.maximum(total - ribo, 0.0)
    ribo_frac = np.divide(ribo, total, out=np.zeros(n), where=total > 0)

    lineage = assign_lineage(expr, n)
    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_epi = lineage == "epithelial"
    is_malig = is_epi & (a3_normal == 0)

    cldn4 = log1p_cp10k(expr["CLDN4"], total)
    tacstd2 = log1p_cp10k(expr["TACSTD2"], total)
    cldn4_nr = log1p_cp10k(expr["CLDN4"], nonribo_lib)
    tac_nr = log1p_cp10k(expr["TACSTD2"], nonribo_lib)

    raw_scores = {}
    nr_scores = {}
    for name, src in SCORE_SRC.items():
        sc, _used = module_vector(expr, MODULES[src], total)
        raw_scores[name] = sc
        nr_scores[name] = module_vector(expr, MODULES[src], nonribo_lib)[0]
    prolif, _p = module_vector(expr, MODULES["prolif"], total)

    def subtype_score(symbols: list[str], lib: np.ndarray) -> np.ndarray:
        cols = []
        for g in symbols:
            if g in expr and float(expr[g].sum()) > 0:
                cols.append(log1p_cp10k(expr[g], lib))
            elif g in aux and bool(np.asarray(aux[f"{g}_present"])[0]) and float(aux[g].sum()) > 0:
                cols.append(log1p_cp10k(aux[g], lib))
        if not cols:
            return np.full(n, np.nan)
        return np.mean(np.vstack(cols), axis=0)

    squamous = subtype_score(["KRT5", "KRT17", "TP63"], total)
    adeno = subtype_score(["KRT7"], total)
    keratin = subtype_score(["KRT8", "KRT18", "KRT19"], total)

    meta = pd.read_csv(HERE / "tables" / "sample_metadata.tsv", sep="\t")
    meta_map = meta.set_index("Sample")
    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])

    df = pd.DataFrame({
        "Sample": sample_ids,
        "total_umi": total,
        "ribo": ribo,
        "ribo_frac": ribo_frac,
        "log_umi": np.log1p(total),
        "log_ribo": np.log1p(ribo),
        "CLDN4": cldn4,
        "TACSTD2": tacstd2,
        "CLDN4_nr": cldn4_nr,
        "TACSTD2_nr": tac_nr,
        "ifn": raw_scores["ifn"],
        "sting": raw_scores["sting"],
        "nhej": raw_scores["nhej"],
        "ifn_nr": nr_scores["ifn"],
        "sting_nr": nr_scores["sting"],
        "nhej_nr": nr_scores["nhej"],
        "prolif": prolif,
        "squamous": squamous,
        "adeno": adeno,
        "keratin": keratin,
        "is_malig": is_malig,
        "is_epi": is_epi,
    })
    df["timing"] = df["Sample"].map(meta_map["timing"])
    df["paper_group"] = df["Sample"].map(meta_map["paper_group"])
    df["pathology"] = df["Sample"].map(meta_map["Pathology"])
    df["recist"] = df["Sample"].map(meta_map["RECIST"])
    df = df.reset_index(drop=True)

    base_scores = {"ifn": "ifn", "sting": "sting", "nhej": "nhej"}
    nr_score_cols = {"ifn": "ifn_nr", "sting": "sting_nr", "nhej": "nhej_nr"}

    post = (df["timing"] == "post").to_numpy()
    pre = (df["timing"] == "pre").to_numpy()
    a3 = df["is_malig"].to_numpy()
    epi = df["is_epi"].to_numpy()

    results: list[dict] = []
    long_rows: list[dict] = []

    def run_split(spec, family, mask, mode, scores=base_scores, splitter="CLDN4"):
        add_paired(results, long_rows, spec, family, patient_split_deltas(df, mask, splitter, scores, mode))

    # Q / tail splits on post A3 malignant cells.
    for mode in ("q4q1", "median", "tertile", "decile", "top20", "detected", "pos_median", "pos_q4q1"):
        run_split(f"post_a3_{mode}", "q_split", post & a3, mode)
    run_split("post_epi_q4q1", "q_split", post & epi, "q4q1")
    run_split("post_epi_median", "q_split", post & epi, "median")
    run_split("pre_a3_q4q1", "pre_vs_post", pre & a3, "q4q1")
    run_split("pre_a3_median", "pre_vs_post", pre & a3, "median")
    run_split("all_a3_q4q1", "pre_vs_post", a3, "q4q1")

    # Response and histology strata. Post surgery labels only.
    run_split("post_a3_nmpr_q4q1", "response", post & a3 & df["paper_group"].eq("NMPR").to_numpy(), "q4q1")
    run_split("post_a3_mpr_q4q1", "response", post & a3 & df["paper_group"].eq("MPR").to_numpy(), "q4q1")
    run_split("post_a3_squamous_q4q1", "response", post & a3 & df["pathology"].eq("Squamous").to_numpy(), "q4q1")
    run_split("post_a3_adeno_q4q1", "response", post & a3 & df["pathology"].eq("Adeno").to_numpy(), "q4q1")
    run_split("post_a3_recist_pr_q4q1", "response", post & a3 & df["recist"].eq("PR").to_numpy(), "q4q1")
    run_split("post_a3_recist_sd_q4q1", "response", post & a3 & df["recist"].eq("SD").to_numpy(), "q4q1")
    run_split("post_a3_nmpr_median", "response", post & a3 & df["paper_group"].eq("NMPR").to_numpy(), "median")
    run_split("post_a3_squamous_median", "response", post & a3 & df["pathology"].eq("Squamous").to_numpy(), "median")

    # Stricter UMI matching.
    add_paired(results, long_rows, "post_a3_caliper_logumi_0.10", "umi_match", patient_caliper(df, post & a3, base_scores, 0.10))
    add_paired(results, long_rows, "post_a3_caliper_logumi_0.25", "umi_match", patient_caliper(df, post & a3, base_scores, 0.25))
    add_paired(results, long_rows, "post_a3_caliper_logumi_0.50", "umi_match", patient_caliper(df, post & a3, base_scores, 0.50))
    add_paired(results, long_rows, "post_a3_umi_quintile", "umi_match", patient_umi_bins(df, post & a3, base_scores, 5))
    add_paired(results, long_rows, "post_a3_umi_decile", "umi_match", patient_umi_bins(df, post & a3, base_scores, 10))
    add_paired(
        results, long_rows, "post_a3_resid_logumi", "umi_match",
        patient_resid_split(df, post & a3, ["log_umi"], base_scores),
    )

    # Ribosomal depth.
    run_split("post_a3_nonribo_lib_q4q1", "ribo", post & a3, "q4q1", nr_score_cols, splitter="CLDN4_nr")
    run_split("post_a3_nonribo_scores_on_raw_cldn4", "ribo", post & a3, "q4q1", nr_score_cols, splitter="CLDN4")
    add_paired(
        results, long_rows, "post_a3_resid_ribo_frac", "ribo",
        patient_resid_split(df, post & a3, ["ribo_frac"], base_scores),
    )
    add_paired(
        results, long_rows, "post_a3_resid_logumi_and_ribo", "ribo",
        patient_resid_split(df, post & a3, ["log_umi", "ribo_frac"], base_scores),
    )
    low_ribo = df["ribo_frac"] <= df.loc[post & a3, "ribo_frac"].median()
    run_split("post_a3_drop_high_ribo", "ribo", post & a3 & low_ribo.to_numpy(), "q4q1")

    # Malignant subtypes.
    noncycling = df["prolif"].fillna(0).to_numpy() <= 0
    cycling = df["prolif"].fillna(0).to_numpy() > 0
    run_split("post_a3_noncycling_q4q1", "subtype", post & a3 & noncycling, "q4q1")
    run_split("post_a3_cycling_q4q1", "subtype", post & a3 & cycling, "q4q1")
    kerat_hi = np.zeros(len(df), dtype=bool)
    kerat_lo = np.zeros(len(df), dtype=bool)
    epcam_hi = np.zeros(len(df), dtype=bool)
    # Within-patient median splits for subtype membership.
    for sid, idx in df.loc[post & a3].groupby("Sample").groups.items():
        idx = np.asarray(list(idx))
        if len(idx) < MIN_CELLS:
            continue
        kmed = np.median(df.loc[idx, "keratin"])
        kerat_hi[idx] = df.loc[idx, "keratin"].to_numpy() > kmed
        kerat_lo[idx] = df.loc[idx, "keratin"].to_numpy() <= kmed
        if "EPCAM" in expr:
            ep = log1p_cp10k(expr["EPCAM"], total)[idx]
            em = np.median(ep)
            epcam_hi[idx] = ep > em
    run_split("post_a3_keratin_high_q4q1", "subtype", post & a3 & kerat_hi, "q4q1")
    run_split("post_a3_keratin_low_q4q1", "subtype", post & a3 & kerat_lo, "q4q1")
    run_split("post_a3_epcam_high_q4q1", "subtype", post & a3 & epcam_hi, "q4q1")
    sq_like = np.isfinite(df["squamous"]) & np.isfinite(df["adeno"]) & (df["squamous"] > df["adeno"]) & (df["squamous"] > 0.2)
    ad_like = np.isfinite(df["squamous"]) & np.isfinite(df["adeno"]) & (df["adeno"] > df["squamous"]) & (df["adeno"] > 0.2)
    run_split("post_a3_squamous_like_q4q1", "subtype", post & a3 & sq_like.to_numpy(), "q4q1")
    run_split("post_a3_adeno_like_q4q1", "subtype", post & a3 & ad_like.to_numpy(), "q4q1")

    # TACSTD2-low ∩ CLDN4-low versus double-high.
    add_paired(results, long_rows, "post_a3_double_low_vs_double_high", "double_low", patient_double_low(df, post & a3, base_scores))
    add_paired(results, long_rows, "post_epi_double_low_vs_double_high", "double_low", patient_double_low(df, post & epi, base_scores))
    tac_low = np.zeros(len(df), dtype=bool)
    for sid, idx in df.loc[post & a3].groupby("Sample").groups.items():
        idx = np.asarray(list(idx))
        if len(idx) < MIN_CELLS:
            continue
        tac_low[idx] = df.loc[idx, "TACSTD2"].to_numpy() <= np.median(df.loc[idx, "TACSTD2"])
    run_split("post_a3_within_tacstd2_low_q4q1", "double_low", post & a3 & tac_low, "q4q1")

    # Continuous correlations. Within-patient rhos use the paired majority rule.
    add_paired(results, long_rows, "post_a3_within_spearman", "continuous", patient_spearman(df, post & a3, base_scores))
    add_paired(results, long_rows, "post_epi_within_spearman", "continuous", patient_spearman(df, post & epi, base_scores))
    add_paired(results, long_rows, "pre_a3_within_spearman", "continuous", patient_spearman(df, pre & a3, base_scores))
    add_paired(results, long_rows, "post_a3_nmpr_within_spearman", "continuous", patient_spearman(df, post & a3 & df["paper_group"].eq("NMPR").to_numpy(), base_scores))
    for spec, mask in (
        ("post_a3_patientmean_spearman", post & a3),
        ("pre_a3_patientmean_spearman", pre & a3),
        ("post_a3_nmpr_patientmean_spearman", post & a3 & df["paper_group"].eq("NMPR").to_numpy()),
        ("post_a3_squamous_patientmean_spearman", post & a3 & df["pathology"].eq("Squamous").to_numpy()),
        ("all_a3_patientmean_spearman", a3),
    ):
        rho, n_pat = mean_spearman(df, mask)
        results.append(judge_corr(spec, "continuous", rho, n_pat))

    # Bulk baseline. One rho per stratum. Response is pathologic MPR vs NMPR.
    bulk = pd.read_csv(HERE / "tables" / "bulk_per_sample.tsv", sep="\t")
    bulk["MPR"] = bulk["Pathologic Response"].astype(str).str.upper().str.startswith("MPR")
    bulk_map = {"ifn": "ifn_effector", "sting": "sting", "nhej": "nhej"}

    def bulk_rho(sub: pd.DataFrame) -> tuple[dict[str, float], int]:
        out = {}
        n_ok = int(len(sub))
        for name, col in bulk_map.items():
            if n_ok < 4 or sub["CLDN4"].nunique() < 2 or sub[col].nunique() < 2:
                out[name] = np.nan
            else:
                out[name] = float(stats.spearmanr(sub["CLDN4"], sub[col]).statistic)
        return out, n_ok

    for spec, sub in (
        ("bulk_pre_all", bulk),
        ("bulk_pre_nmpr", bulk[~bulk["MPR"]]),
        ("bulk_pre_mpr", bulk[bulk["MPR"]]),
        ("bulk_pre_squamous", bulk[bulk["Pathology"].astype(str).str.contains("Squamous", case=False)]),
        ("bulk_pre_adeno", bulk[bulk["Pathology"].astype(str).str.contains("Adeno", case=False)]),
    ):
        rho, n_pat = bulk_rho(sub)
        results.append(judge_corr(spec, "bulk", rho, n_pat))
        # Quartile of samples: one unpaired contrast, stored as a 1-row paired judge only if both tails have >= 4.
        if len(sub) >= 8 and sub["CLDN4"].nunique() > 1:
            q75, q25 = sub["CLDN4"].quantile(0.75), sub["CLDN4"].quantile(0.25)
            hi, lo = sub[sub["CLDN4"] >= q75], sub[sub["CLDN4"] <= q25]
            if len(hi) >= 4 and len(lo) >= 4 and q75 > q25:
                rec = {"Sample": "bulk_samples", "n_high": int(len(hi)), "n_low": int(len(lo))}
                for name, col in bulk_map.items():
                    rec[name] = float(hi[col].mean() - lo[col].mean())
                # Single contrast: majority rule is just the sign. Require n_samples in the stratum >= 8
                # by duplicating the sign across a note, and use judge_paired only if we expand.
                # Treat as correlation-style magnitude using Cohen-free sign plus arm size via judge on one row
                # which fails majority (1 vs 0 is a majority). n=1 fails MIN_N_PAIRED. Record with kind bulk_split
                # and hit only if all signs match and each arm has >= 4 samples (n_patients reported as n_high+n_low).
                signs = {
                    "ifn": rec["ifn"] < 0,
                    "sting": rec["sting"] < 0,
                    "nhej": rec["nhej"] > 0,
                }
                results.append({
                    "spec": spec + "_sample_q4q1",
                    "family": "bulk",
                    "kind": "bulk_split",
                    "n": int(len(hi) + len(lo)),
                    "ifn_median": rec["ifn"],
                    "sting_median": rec["sting"],
                    "nhej_median": rec["nhej"],
                    "ifn_n_thesis": int(signs["ifn"]),
                    "sting_n_thesis": int(signs["sting"]),
                    "nhej_n_thesis": int(signs["nhej"]),
                    "ifn_n_opposite": int(not signs["ifn"]),
                    "sting_n_opposite": int(not signs["sting"]),
                    "nhej_n_opposite": int(not signs["nhej"]),
                    "ifn_sign_ok": signs["ifn"],
                    "sting_sign_ok": signs["sting"],
                    "nhej_sign_ok": signs["nhej"],
                    "joint_sign": all(signs.values()),
                    "reliable": bool(all(signs.values()) and len(hi) >= 4 and len(lo) >= 4),
                    "hit": bool(all(signs.values()) and len(hi) >= 4 and len(lo) >= 4),
                })
        rho_p = {}
        for name, col in bulk_map.items():
            rho_p[name] = bulk_partial(sub["CLDN4"].to_numpy(), sub[col].to_numpy(), sub["EPCAM"].to_numpy())
        results.append(judge_corr(spec + "_partial_epcam", "bulk", rho_p, int(len(sub))))

    # Bulk double-low vs double-high, sample means.
    if len(bulk) >= 8:
        cmed, tmed = bulk["CLDN4"].median(), bulk["TACSTD2"].median()
        hi = bulk[(bulk["CLDN4"] > cmed) & (bulk["TACSTD2"] > tmed)]
        lo = bulk[(bulk["CLDN4"] <= cmed) & (bulk["TACSTD2"] <= tmed)]
        if len(hi) >= 4 and len(lo) >= 4:
            signs = {}
            meds = {}
            for name, col in bulk_map.items():
                meds[name] = float(hi[col].mean() - lo[col].mean())
            signs = {"ifn": meds["ifn"] < 0, "sting": meds["sting"] < 0, "nhej": meds["nhej"] > 0}
            results.append({
                "spec": "bulk_double_low_vs_double_high",
                "family": "double_low",
                "kind": "bulk_split",
                "n": int(len(hi) + len(lo)),
                "ifn_median": meds["ifn"],
                "sting_median": meds["sting"],
                "nhej_median": meds["nhej"],
                "ifn_n_thesis": int(signs["ifn"]),
                "sting_n_thesis": int(signs["sting"]),
                "nhej_n_thesis": int(signs["nhej"]),
                "ifn_n_opposite": int(not signs["ifn"]),
                "sting_n_opposite": int(not signs["sting"]),
                "nhej_n_opposite": int(not signs["nhej"]),
                "ifn_sign_ok": signs["ifn"],
                "sting_sign_ok": signs["sting"],
                "nhej_sign_ok": signs["nhej"],
                "joint_sign": all(signs.values()),
                "reliable": bool(all(signs.values())),
                "hit": bool(all(signs.values())),
            })

    res = pd.DataFrame(results)
    res["n_modules_thesis_sign"] = (
        res["ifn_sign_ok"].astype(bool).astype(int)
        + res["sting_sign_ok"].astype(bool).astype(int)
        + res["nhej_sign_ok"].astype(bool).astype(int)
    )
    # For correlations, sign_ok already includes the |rho| floor. Also store raw joint sign
    # without the magnitude floor so a tiny flip is visible and not called a hit.
    res = res.sort_values(["hit", "joint_sign", "n_modules_thesis_sign", "spec"], ascending=[False, False, False, True])
    tabdir = HERE / "tables"
    figdir = HERE / "figures"
    tabdir.mkdir(exist_ok=True)
    figdir.mkdir(exist_ok=True)
    res.to_csv(tabdir / "sweep_specs.tsv", sep="\t", index=False)
    pd.DataFrame(long_rows).to_csv(tabdir / "sweep_patient_deltas.tsv", sep="\t", index=False)

    hits = res[res["hit"]].copy()
    joint = res[res["joint_sign"]].copy()
    aux_note = {
        "n_ribo_genes": int(np.asarray(aux["n_ribo_genes"])[0]),
        "n_mito_genes": int(np.asarray(aux["n_mito_genes"])[0]),
        "subtype_genes_present": {
            g: bool(np.asarray(aux[f"{g}_present"])[0]) for g in ("KRT5", "KRT6A", "KRT7", "KRT17", "TP63", "SOX2")
        },
        "median_ribo_frac_post_a3": float(df.loc[post & a3, "ribo_frac"].median()),
    }
    verdict = {
        "dataset": "GSE207422",
        "thesis": "CLDN4-low → IFN↑ STING↑ NHEJ↓",
        "n_specs": int(len(res)),
        "n_joint_sign": int(joint.shape[0]),
        "n_hits": int(hits.shape[0]),
        "hit_specs": hits["spec"].tolist(),
        "joint_sign_specs": joint["spec"].tolist(),
        "final_status": "FINAL_DISCORDANT" if hits.empty else "THESIS_ALIGNED_HIT",
        "rule": {
            "paired": f"n>={MIN_N_PAIRED} and strict majority on each module",
            "correlation": f"n>={MIN_N_CORR} and |rho|>={MIN_ABS_RHO} on each module",
            "bulk_split": "both arms n>=4 and all three mean differences in the thesis direction",
        },
        "aux": aux_note,
    }
    (tabdir / "sweep_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2), flush=True)
    show = res[["spec", "family", "kind", "n", "ifn_median", "sting_median", "nhej_median", "ifn_sign_ok", "sting_sign_ok", "nhej_sign_ok", "joint_sign", "hit"]]
    print(show.to_string(index=False), flush=True)

    # Sign map. Thesis direction is the lower row position.
    plot = res.sort_values(["family", "spec"]).reset_index(drop=True)
    fig_h = max(6.0, 0.28 * len(plot) + 1.2)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    colors = {True: "#1b7f4e", False: "#b23a48"}
    for i, r in plot.iterrows():
        for j, name in enumerate(("ifn_sign_ok", "sting_sign_ok", "nhej_sign_ok")):
            ax.scatter(j, i, s=36, color=colors[bool(r[name])], zorder=3)
        if bool(r["hit"]):
            ax.scatter(3.15, i, marker="*", s=70, color="#1b7f4e", zorder=3)
    ax.set_yticks(np.arange(len(plot)), plot["spec"], fontsize=7)
    ax.set_xticks([0, 1, 2, 3.15], ["IFN thesis", "STING thesis", "NHEJ thesis", "hit"])
    ax.set_xlim(-0.5, 3.6)
    ax.invert_yaxis()
    ax.set_title(f"GSE207422 sweep ({verdict['final_status']})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(figdir / "fig_sweep_signs.png", dpi=140, bbox_inches="tight")
    fig.savefig(figdir / "fig_sweep_signs.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
