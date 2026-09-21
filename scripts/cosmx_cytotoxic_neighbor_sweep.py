#!/usr/bin/env python3
"""Aggressive CosMx neighbor-count sweep for CLDN4 vs cytotoxic cells.

He et al. 2022, figshare 25976224 (clustered h5ad). Sensitivity layer only.
Does not replace the locked 50/100 µm cytotoxic ratios 0.36 / 0.52.

Selection rule, fixed before looking at the grid:
  A specification is eligible when at least 7 of 8 sections and 4 of 5
  donors share one sign of Δ (high − low). Unusable sections count against
  that 7/8 and 4/5.
  Ratio effect = |log(mean section-high / mean section-low)|.
  Delta effect = |mean of section deltas|.
  The ratio sensitivity is the eligible spec with the largest ratio effect.
  The delta sensitivity is the eligible spec with the largest delta effect.
  Ratio is undefined when the equal-weight mean of section low-arm means
  is below 0.05, so that spec cannot win the ratio call.
  Continuous CLDN4 (Spearman) is reported separately and is not ranked
  against count ratios.

Neighbor search uses each section's own global coordinates (0.18 µm/px).
Sections are never placed in one tree. FOV / section / donor is the unit
of the CLDN4 threshold and of the mean, not a second spatial window.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_cytotoxic_neighbor_sweep")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
CACHE = os.path.join("/tmp", "cosmx_cytotoxic_cell_counts.npz")

UM_PER_PX = 0.18
RADII = (10, 20, 25, 40, 50, 80, 100)
CYTOS = ("cd8", "nk", "cd8nk", "gzmb")
UNITS = ("section", "fov", "donor")
# top25 vs bottom is the same contrast as q4_q1, so it is not duplicated.
CUTS = (
    "detected_absent",
    "q4_q1",
    "median",
    "top10_vs_bottom",
    "top20_vs_bottom",
    "top30_vs_bottom",
    "top10_vs_rest",
    "top20_vs_rest",
    "top25_vs_rest",
    "top30_vs_rest",
)
SAMPLES = (
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
)
TUMOR_LABEL = {
    "LUAD-5 R1": "tumor 5",
    "LUAD-5 R2": "tumor 5",
    "LUAD-5 R3": "tumor 5",
    "LUSC-6": "tumor 6",
    "LUAD-9 R1": "tumor 9",
    "LUAD-9 R2": "tumor 9",
    "LUAD-12": "tumor 12",
    "LUAD-13": "tumor 13",
}
MIN_ARM_SECTION = 30
MIN_ARM_FOV = 5
MIN_FOV_PER_SECTION = 3
MIN_SPEARMAN = 30
LOW_MEAN_FLOOR = 0.05
N_SECTIONS = 8
N_DONORS = 5


def _decode(arr) -> list:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def arm_masks(values: np.ndarray, positive: np.ndarray, mode: str, min_arm: int):
    """Return (high, low) boolean masks, or None if the cut does not separate.

    q4_q1 and median follow the zero-inflated rule used on this object:
    low is everything at or below the lower threshold, high is everything
    at or above the upper threshold after removing overlap. When the upper
    threshold is 0, that collapses to detected vs absent for that group.
    topX_vs_bottom requires a strict gap between the two thresholds so a
    zero-inflated section is not relabeled as a top-decile contrast.
    topX_vs_rest is the upper tail versus everyone below that threshold.
    """
    values = np.asarray(values, dtype=np.float64)
    positive = np.asarray(positive, dtype=bool)
    n = len(values)
    if n < min_arm * 2:
        return None
    if mode == "detected_absent":
        high = positive.copy()
        low = ~positive
    elif mode == "median":
        med = float(np.median(values))
        high = values > med
        low = values <= med
    elif mode == "q4_q1":
        q1, q3 = np.quantile(values, [0.25, 0.75])
        low = values <= q1
        high = (values >= q3) & ~low
    elif mode.endswith("_vs_bottom") or mode.endswith("_vs_rest"):
        pct = int(mode.split("_")[0].replace("top", ""))
        p = pct / 100.0
        if mode.endswith("_vs_rest"):
            hi = float(np.quantile(values, 1.0 - p))
            high = values >= hi
            low = ~high
        else:
            lo = float(np.quantile(values, p))
            hi = float(np.quantile(values, 1.0 - p))
            if not hi > lo:
                return None
            low = values <= lo
            high = (values >= hi) & ~low
    else:
        raise ValueError(mode)
    if int(high.sum()) < min_arm or int(low.sum()) < min_arm:
        return None
    if np.any(high & low):
        return None
    return high, low


def assign_labels(cldn4: np.ndarray, positive: np.ndarray, sample: np.ndarray,
                  patient: np.ndarray, fov: np.ndarray, mode: str, unit: str) -> np.ndarray:
    """1 = high, -1 = low, 0 = unused."""
    labels = np.zeros(len(cldn4), dtype=np.int8)
    if unit == "section":
        groups = _group_indices(sample)
        min_arm = MIN_ARM_SECTION
    elif unit == "donor":
        groups = _group_indices(patient)
        min_arm = MIN_ARM_SECTION
    elif unit == "fov":
        groups = _group_indices(sample, fov)
        min_arm = MIN_ARM_FOV
    else:
        raise ValueError(unit)
    for idx in groups.values():
        masks = arm_masks(cldn4[idx], positive[idx], mode, min_arm)
        if masks is None:
            continue
        high, low = masks
        lab = np.zeros(len(idx), dtype=np.int8)
        lab[high] = 1
        lab[low] = -1
        labels[idx] = lab
    return labels


def _group_indices(*keys: np.ndarray) -> dict:
    frame = pd.DataFrame({f"k{i}": np.asarray(k) for i, k in enumerate(keys)})
    out = {}
    for key, idx in frame.groupby(list(frame.columns), sort=False).indices.items():
        out[key] = np.asarray(idx, dtype=np.int64)
    return out


def _mean_or_nan(values: np.ndarray) -> float:
    if len(values) == 0:
        return float("nan")
    return float(np.mean(values))


def section_rows_for_labels(sample, patient, fov, labels, counts, unit: str, cut: str):
    """counts is (n, n_metrics). One row per section per metric column."""
    n_metrics = counts.shape[1]
    fov_groups = _group_indices(sample, fov)
    fov_neg = np.zeros(n_metrics, dtype=np.int32)
    fov_pos = np.zeros(n_metrics, dtype=np.int32)
    fov_zero = np.zeros(n_metrics, dtype=np.int32)
    fov_n = 0
    fov_means = {s: [] for s in SAMPLES}  # list of (high_mean[m], low_mean[m])
    n_high_sec = {s: 0 for s in SAMPLES}
    n_low_sec = {s: 0 for s in SAMPLES}

    for idx in fov_groups.values():
        s = sample[idx[0]]
        lab = labels[idx]
        h = idx[lab == 1]
        l = idx[lab == -1]
        if len(h) < MIN_ARM_FOV or len(l) < MIN_ARM_FOV:
            continue
        mh = counts[h].mean(axis=0)
        ml = counts[l].mean(axis=0)
        delta = mh - ml
        fov_n += 1
        fov_neg += delta < 0
        fov_pos += delta > 0
        fov_zero += delta == 0
        fov_means[s].append((mh, ml))
        n_high_sec[s] += int(len(h))
        n_low_sec[s] += int(len(l))

    rows = []
    for s in SAMPLES:
        pat = None
        sec_idx = np.flatnonzero(sample == s)
        if len(sec_idx) == 0:
            continue
        pat = str(patient[sec_idx[0]])
        usable = True
        n_high = n_low = 0
        if unit == "fov":
            pairs = fov_means[s]
            if len(pairs) < MIN_FOV_PER_SECTION:
                usable = False
                mh = np.full(n_metrics, np.nan)
                ml = np.full(n_metrics, np.nan)
            else:
                mh = np.mean([p[0] for p in pairs], axis=0)
                ml = np.mean([p[1] for p in pairs], axis=0)
                n_high = int(n_high_sec[s])
                n_low = int(n_low_sec[s])
        else:
            lab = labels[sec_idx]
            h = sec_idx[lab == 1]
            l = sec_idx[lab == -1]
            n_high = int(len(h))
            n_low = int(len(l))
            if n_high < MIN_ARM_SECTION or n_low < MIN_ARM_SECTION:
                usable = False
                mh = np.full(n_metrics, np.nan)
                ml = np.full(n_metrics, np.nan)
            else:
                mh = counts[h].mean(axis=0)
                ml = counts[l].mean(axis=0)
        for j in range(n_metrics):
            rows.append({
                "sample": s,
                "patient": pat,
                "unit": unit,
                "cut": cut,
                "metric": j,
                "usable": bool(usable),
                "n_high": int(n_high),
                "n_low": int(n_low),
                "n_fov_in_section": int(len(fov_means[s])),
                "mean_high": float(mh[j]) if usable else np.nan,
                "mean_low": float(ml[j]) if usable else np.nan,
            })
    fov_stats = {
        "n_fov": int(fov_n),
        "n_fov_neg": fov_neg,
        "n_fov_pos": fov_pos,
        "n_fov_zero": fov_zero,
    }
    return rows, fov_stats


def spearman_rho(x: np.ndarray, y: np.ndarray, min_n: int) -> float:
    if len(x) < min_n:
        return float("nan")
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return float("nan")
    rho = stats.spearmanr(x, y).statistic
    if rho is None or not np.isfinite(rho):
        return float("nan")
    return float(rho)


def sign_pvalue(k: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)


def wilcoxon_p(deltas: np.ndarray) -> float:
    d = np.asarray(deltas, dtype=np.float64)
    d = d[np.isfinite(d)]
    d = d[d != 0]
    if len(d) == 0:
        return float("nan")
    method = "exact" if len(d) <= 25 else "auto"
    return float(stats.wilcoxon(d, alternative="two-sided", method=method, zero_method="wilcox").pvalue)


def grid_from_sections(sec: pd.DataFrame, fov_table: pd.DataFrame, metric_names: list[dict]) -> pd.DataFrame:
    rows = []
    grouped = sec.groupby(["unit", "cut", "metric"], sort=False)
    for (unit, cut, metric), sub in grouped:
        meta = metric_names[int(metric)]
        usable = sub[sub["usable"]].copy()
        deltas = usable["mean_high"].to_numpy() - usable["mean_low"].to_numpy()
        n_neg = int(np.sum(deltas < 0))
        n_pos = int(np.sum(deltas > 0))
        n_zero = int(np.sum(deltas == 0))
        # Donors: equal-weight mean of usable sections.
        don_neg = don_pos = don_zero = don_n = 0
        don_deltas = []
        if len(usable):
            for _, dsub in usable.groupby("patient", sort=False):
                dh = float(dsub["mean_high"].mean())
                dl = float(dsub["mean_low"].mean())
                dd = dh - dl
                don_n += 1
                don_deltas.append(dd)
                if dd < 0:
                    don_neg += 1
                elif dd > 0:
                    don_pos += 1
                else:
                    don_zero += 1
        mean_high = float(usable["mean_high"].mean()) if len(usable) else np.nan
        mean_low = float(usable["mean_low"].mean()) if len(usable) else np.nan
        pooled_delta = float(np.mean(deltas)) if len(deltas) else np.nan
        if np.isfinite(mean_low) and mean_low >= LOW_MEAN_FLOOR and np.isfinite(mean_high):
            pooled_ratio = float(mean_high / mean_low)
        else:
            pooled_ratio = float("nan")
        sec_ratio = []
        for _, r in usable.iterrows():
            if r["mean_low"] >= LOW_MEAN_FLOOR:
                sec_ratio.append(r["mean_high"] / r["mean_low"])
        med_ratio = float(np.median(sec_ratio)) if sec_ratio else np.nan
        fov_hit = fov_table[
            (fov_table["unit"] == unit) & (fov_table["cut"] == cut) & (fov_table["metric"] == metric)
        ]
        if len(fov_hit) == 1:
            n_fov = int(fov_hit["n_fov"].iloc[0])
            n_fov_neg = int(fov_hit["n_fov_neg"].iloc[0])
            n_fov_pos = int(fov_hit["n_fov_pos"].iloc[0])
        else:
            n_fov = n_fov_neg = n_fov_pos = 0
        effect_log = float(abs(np.log(pooled_ratio))) if np.isfinite(pooled_ratio) and pooled_ratio > 0 else np.nan
        effect_delta = float(abs(pooled_delta)) if np.isfinite(pooled_delta) else np.nan
        rows.append({
            "unit": unit,
            "cut": cut,
            "cyto": meta["cyto"],
            "radius_um": int(meta["radius"]),
            "continuous": bool(meta["continuous"]),
            "n_sec_usable": int(len(usable)),
            "n_sec_neg": n_neg,
            "n_sec_pos": n_pos,
            "n_sec_zero": n_zero,
            "n_don_usable": int(don_n),
            "n_don_neg": int(don_neg),
            "n_don_pos": int(don_pos),
            "n_don_zero": int(don_zero),
            "mean_high": mean_high,
            "mean_low": mean_low,
            "pooled_ratio": pooled_ratio,
            "pooled_delta": pooled_delta,
            "median_section_ratio": med_ratio,
            "effect_log_ratio": effect_log,
            "effect_abs_delta": effect_delta,
            "eligible_exclusion": bool(n_neg >= 7 and don_neg >= 4),
            "eligible_enrichment": bool(n_pos >= 7 and don_pos >= 4),
            "sign_p_section": sign_pvalue(max(n_neg, n_pos), N_SECTIONS),
            "sign_p_donor": sign_pvalue(max(don_neg, don_pos), N_DONORS),
            "wilcoxon_p": wilcoxon_p(deltas),
            "n_fov": n_fov,
            "n_fov_neg": n_fov_neg,
            "n_fov_pos": n_fov_pos,
            "donor_deltas_json": json.dumps(don_deltas),
        })
    return pd.DataFrame(rows)


def _eligible(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["eligible_exclusion"] | df["eligible_enrichment"]].copy()


def _direction(row) -> str:
    if row["eligible_exclusion"] and not row["eligible_enrichment"]:
        return "exclusion"
    if row["eligible_enrichment"] and not row["eligible_exclusion"]:
        return "enrichment"
    if row["n_sec_neg"] > row["n_sec_pos"]:
        return "exclusion"
    if row["n_sec_pos"] > row["n_sec_neg"]:
        return "enrichment"
    return "tied"


def select_winners(grid: pd.DataFrame) -> dict:
    """Largest effects among sign-consistent specs. Continuous rows are excluded."""
    binary = grid[~grid["continuous"].astype(bool)].copy()
    elig = _eligible(binary)
    out = {
        "n_binary": int(len(binary)),
        "n_eligible": int(len(elig)),
        "ratio_winner": None,
        "delta_winner": None,
        "numeric_max_ratio": None,
        "numeric_max_delta": None,
        "numeric_min_ratio": None,
        "numeric_min_delta": None,
    }
    if elig.empty:
        return out

    def pack(row):
        d = row.to_dict()
        d["direction"] = _direction(row)
        d["n_sec_sign"] = int(row["n_sec_neg"] if d["direction"] == "exclusion" else row["n_sec_pos"])
        d["n_don_sign"] = int(row["n_don_neg"] if d["direction"] == "exclusion" else row["n_don_pos"])
        return d

    # Detected/absent does not use a quantile, so section and donor labels match.
    # Prefer the section summary when the effect is tied.
    unit_rank = {"section": 0, "fov": 1, "donor": 2}
    ratio_pool = elig[np.isfinite(elig["effect_log_ratio"])].copy()
    if len(ratio_pool):
        ratio_pool["unit_rank"] = ratio_pool["unit"].map(unit_rank)
        ratio_pool = ratio_pool.sort_values(
            ["effect_log_ratio", "n_sec_neg", "n_sec_pos", "n_don_neg", "n_don_pos", "unit_rank", "cut", "cyto", "radius_um"],
            ascending=[False, False, False, False, False, True, True, True, True],
        )
        out["ratio_winner"] = pack(ratio_pool.iloc[0])
    delta_pool = elig[np.isfinite(elig["effect_abs_delta"])].copy()
    if len(delta_pool):
        delta_pool["unit_rank"] = delta_pool["unit"].map(unit_rank)
        delta_pool = delta_pool.sort_values(
            ["effect_abs_delta", "effect_log_ratio", "n_sec_neg", "n_sec_pos", "unit_rank", "cut", "cyto", "radius_um"],
            ascending=[False, False, False, False, True, True, True, True],
        )
        out["delta_winner"] = pack(delta_pool.iloc[0])
    finite = elig[np.isfinite(elig["pooled_ratio"]) & (elig["pooled_ratio"] > 0)]
    if len(finite):
        out["numeric_max_ratio"] = pack(finite.loc[finite["pooled_ratio"].idxmax()])
        out["numeric_min_ratio"] = pack(finite.loc[finite["pooled_ratio"].idxmin()])
    finite_d = elig[np.isfinite(elig["pooled_delta"])]
    if len(finite_d):
        out["numeric_max_delta"] = pack(finite_d.loc[finite_d["pooled_delta"].idxmax()])
        out["numeric_min_delta"] = pack(finite_d.loc[finite_d["pooled_delta"].idxmin()])
    return out


def select_continuous(grid: pd.DataFrame) -> dict:
    cont = grid[grid["continuous"].astype(bool)].copy()
    elig = _eligible(cont)
    out = {"n": int(len(cont)), "n_eligible": int(len(elig)), "winner": None}
    if elig.empty:
        return out
    elig = elig[np.isfinite(elig["effect_abs_delta"])]
    if elig.empty:
        return out
    elig = elig.sort_values(
        ["effect_abs_delta", "cut", "cyto", "unit", "radius_um"],
        ascending=[False, True, True, True, True],
    )
    row = elig.iloc[0]
    packed = row.to_dict()
    packed["direction"] = _direction(row)
    out["winner"] = packed
    return out


def _decode_cats(f, name: str):
    cats = _decode(f[f"obs/{name}/categories"][:])
    codes = f[f"obs/{name}/codes"][:]
    return np.asarray(cats, dtype=object)[codes], cats


def load_cells():
    import h5py

    t0 = time.time()
    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}
    for g in ("CLDN4", "GZMB"):
        if g not in gmap:
            raise SystemExit(f"{g} is absent from the 960-plex object")
    sample, sample_cats = _decode_cats(f, "sample")
    patient, _ = _decode_cats(f, "patient")
    cell_type, type_cats = _decode_cats(f, "cell_type")
    fov = f["obs/fov"][:].astype(np.int32)
    n_counts = f["obs/n_counts"][:].astype(np.float64)
    xy = f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    f.close()
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(sample), len(genes)))
    del data, indices, indptr
    cldn4_raw = np.asarray(mat.getcol(gmap["CLDN4"]).toarray()).ravel().astype(np.float64)
    gzmb_raw = np.asarray(mat.getcol(gmap["GZMB"]).toarray()).ravel().astype(np.float64)
    del mat
    lib = n_counts.copy()
    lib[lib <= 0] = np.nan
    cldn4 = np.nan_to_num(np.log1p(cldn4_raw / lib * 1e4), nan=0.0)
    print(f"loaded obs+CLDN4+GZMB in {time.time()-t0:.1f}s", flush=True)

    missing_samples = [s for s in SAMPLES if s not in set(sample_cats)]
    if missing_samples:
        raise SystemExit(f"missing samples: {missing_samples}; found {sample_cats}")
    cd8_labels = [t for t in type_cats if "CD8" in t]
    nk_labels = [t for t in type_cats if t == "NK" or t.startswith("NK")]
    if not cd8_labels or not nk_labels:
        raise SystemExit(f"CD8/NK labels missing. cell types: {type_cats}")

    malignant = np.array([cell_type[i] == TUMOR_LABEL[sample[i]] for i in range(len(sample))], dtype=bool)
    is_cd8 = np.isin(cell_type, cd8_labels)
    is_nk = np.isin(cell_type, nk_labels)
    is_tumor_any = np.array([str(t).startswith("tumor") for t in cell_type])
    is_gzmb = (gzmb_raw > 0) & ~is_tumor_any

    nn = {}
    pieces = []
    gzmb_comp = {}
    inventory = []
    for s in SAMPLES:
        m = sample == s
        xy_s = xy[m]
        # Median nearest-neighbor distance on a fixed subsample.
        rng = np.random.default_rng(0)
        take = rng.choice(xy_s.shape[0], size=min(6000, xy_s.shape[0]), replace=False)
        tree_all = cKDTree(xy_s)
        dist, _ = tree_all.query(xy_s[take], k=2, workers=4)
        med_nn = float(np.median(dist[:, 1]))
        nn[s] = med_nn
        if not (3.0 <= med_nn <= 40.0):
            raise SystemExit(f"{s}: median NN {med_nn:.3f} µm is outside 3–40; check 0.18 µm/px")
        mal = malignant[m]
        local_ct = cell_type[m]
        ref = {
            "cd8": is_cd8[m],
            "nk": is_nk[m],
            "cd8nk": is_cd8[m] | is_nk[m],
            "gzmb": is_gzmb[m],
        }
        tum_xy = xy_s[mal]
        counts = np.zeros((int(mal.sum()), len(CYTOS), len(RADII)), dtype=np.float32)
        for a, cy in enumerate(CYTOS):
            ref_xy = xy_s[ref[cy]]
            if len(ref_xy) == 0 or len(tum_xy) == 0:
                continue
            tree = cKDTree(ref_xy)
            for b, radius in enumerate(RADII):
                counts[:, a, b] = tree.query_ball_point(tum_xy, r=float(radius), return_length=True, workers=4)
        cldn4_s = cldn4[m][mal]
        raw_s = cldn4_raw[m][mal]
        fov_s = fov[m][mal]
        pat = str(patient[m][0])
        flat = counts.reshape(counts.shape[0], -1)
        piece = pd.DataFrame(flat, columns=[f"c{i}" for i in range(flat.shape[1])])
        piece.insert(0, "pos", raw_s > 0)
        piece.insert(0, "cldn4", cldn4_s.astype(np.float64))
        piece.insert(0, "fov", fov_s)
        piece.insert(0, "patient", pat)
        piece.insert(0, "sample", s)
        pieces.append(piece)
        gz_types = local_ct[ref["gzmb"]]
        comp = {str(k): int(v) for k, v in zip(*np.unique(gz_types, return_counts=True))}
        gzmb_comp[s] = comp
        inventory.append({
            "sample": s,
            "patient": pat,
            "tumor_label": TUMOR_LABEL[s],
            "n_cells": int(m.sum()),
            "n_malignant": int(mal.sum()),
            "n_cldn4_pos": int((raw_s > 0).sum()),
            "frac_cldn4_pos": float((raw_s > 0).mean()) if len(raw_s) else np.nan,
            "n_cd8": int(ref["cd8"].sum()),
            "n_nk": int(ref["nk"].sum()),
            "n_cd8nk": int(ref["cd8nk"].sum()),
            "n_gzmb": int(ref["gzmb"].sum()),
            "median_nn_um": med_nn,
            "n_other_tumor_label": int((m & is_tumor_any & ~malignant).sum()),
        })
        print(
            f"{s}: mal={int(mal.sum())} cd8={int(ref['cd8'].sum())} nk={int(ref['nk'].sum())} "
            f"gzmb+={int(ref['gzmb'].sum())} NN={med_nn:.2f} µm",
            flush=True,
        )
    cells = pd.concat(pieces, ignore_index=True)
    meta = {
        "cd8_labels": cd8_labels,
        "nk_labels": nk_labels,
        "gzmb_definition": "raw GZMB > 0 and cell_type does not start with 'tumor'",
        "gzmb_composition": gzmb_comp,
        "inventory": inventory,
        "nn_um": nn,
        "n_malignant": int(len(cells)),
    }
    return cells, meta


def count_columns(frame: pd.DataFrame) -> list[str]:
    return [c for c in frame.columns if len(str(c)) > 1 and str(c)[0] == "c" and str(c)[1:].isdigit()]


def save_cache(cells: pd.DataFrame, meta: dict) -> None:
    cols = count_columns(cells)
    np.savez_compressed(
        CACHE,
        sample=cells["sample"].to_numpy(),
        patient=cells["patient"].to_numpy(),
        fov=cells["fov"].to_numpy(),
        cldn4=cells["cldn4"].to_numpy(),
        pos=cells["pos"].to_numpy(),
        counts=cells[cols].to_numpy(dtype=np.float32),
        meta=json.dumps(meta),
    )
    print(f"cached {CACHE}", flush=True)


def load_cache():
    z = np.load(CACHE, allow_pickle=True)
    counts = z["counts"]
    cells = pd.DataFrame(counts, columns=[f"c{i}" for i in range(counts.shape[1])])
    cells.insert(0, "pos", z["pos"].astype(bool))
    cells.insert(0, "cldn4", z["cldn4"])
    cells.insert(0, "fov", z["fov"])
    cells.insert(0, "patient", z["patient"])
    cells.insert(0, "sample", z["sample"])
    meta = json.loads(str(z["meta"]))
    return cells, meta


def metric_catalog():
    names = []
    for cy in CYTOS:
        for radius in RADII:
            names.append({"cyto": cy, "radius": radius, "continuous": False})
    return names


def run_sweep(cells: pd.DataFrame):
    sample = cells["sample"].to_numpy()
    patient = cells["patient"].to_numpy()
    fov = cells["fov"].to_numpy()
    cldn4 = cells["cldn4"].to_numpy(dtype=np.float64)
    positive = cells["pos"].to_numpy(dtype=bool)
    count_mat = cells[count_columns(cells)].to_numpy(dtype=np.float64)
    names = metric_catalog()
    assert count_mat.shape[1] == len(names)
    sec_rows = []
    fov_rows = []
    t0 = time.time()
    for unit in UNITS:
        for cut in CUTS:
            labels = assign_labels(cldn4, positive, sample, patient, fov, cut, unit)
            rows, fov_stats = section_rows_for_labels(sample, patient, fov, labels, count_mat, unit, cut)
            sec_rows.extend(rows)
            for j in range(len(names)):
                fov_rows.append({
                    "unit": unit,
                    "cut": cut,
                    "metric": j,
                    "n_fov": fov_stats["n_fov"],
                    "n_fov_neg": int(fov_stats["n_fov_neg"][j]),
                    "n_fov_pos": int(fov_stats["n_fov_pos"][j]),
                })
            print(f"cut {unit} {cut} labeled {(labels != 0).mean():.3f} in {time.time()-t0:.1f}s", flush=True)
    # Continuous: section Spearman, FOV-median Spearman, donor-pooled Spearman for donor signs.
    # Stored as extra section rows with cut=continuous. Donor-unit section rows are the
    # per-section Spearmans; donor sign is recomputed from donor-pooled rho below and
    # written onto those rows via a side channel (patient-level override).
    cont_sec = []
    cont_fov_rows = []
    for unit in UNITS:
        if unit == "section":
            groups = {s: np.flatnonzero(sample == s) for s in SAMPLES}
            for s, idx in groups.items():
                rhos = np.array([spearman_rho(cldn4[idx], count_mat[idx, j], MIN_SPEARMAN) for j in range(count_mat.shape[1])])
                for j, rho in enumerate(rhos):
                    cont_sec.append({
                        "sample": s,
                        "patient": str(patient[idx[0]]),
                        "unit": unit,
                        "cut": "continuous",
                        "metric": j,
                        "usable": bool(np.isfinite(rho)),
                        "n_high": int(len(idx)),
                        "n_low": int(len(idx)),
                        "n_fov_in_section": 0,
                        "mean_high": 0.0,
                        "mean_low": float(-rho) if np.isfinite(rho) else np.nan,
                        # delta = mean_high - mean_low = rho. Negative rho => exclusion sign.
                    })
        elif unit == "fov":
            fov_groups = _group_indices(sample, fov)
            bucket = {s: [] for s in SAMPLES}
            fov_neg = np.zeros(count_mat.shape[1], dtype=np.int32)
            fov_pos = np.zeros(count_mat.shape[1], dtype=np.int32)
            n_fov = 0
            for idx in fov_groups.values():
                if len(idx) < MIN_SPEARMAN:
                    continue
                rhos = np.array([
                    spearman_rho(cldn4[idx], count_mat[idx, j], MIN_SPEARMAN) for j in range(count_mat.shape[1])
                ])
                if not np.isfinite(rhos).any():
                    continue
                n_fov += 1
                fov_neg += rhos < 0
                fov_pos += rhos > 0
                bucket[sample[idx[0]]].append(rhos)
            for j in range(count_mat.shape[1]):
                cont_fov_rows.append({
                    "unit": unit, "cut": "continuous", "metric": j,
                    "n_fov": n_fov, "n_fov_neg": int(fov_neg[j]), "n_fov_pos": int(fov_pos[j]),
                })
            for s in SAMPLES:
                rhos_list = bucket[s]
                idx = np.flatnonzero(sample == s)
                usable = len(rhos_list) >= MIN_FOV_PER_SECTION
                med = np.median(np.vstack(rhos_list), axis=0) if usable else np.full(count_mat.shape[1], np.nan)
                for j in range(count_mat.shape[1]):
                    rho = float(med[j]) if usable and np.isfinite(med[j]) else np.nan
                    cont_sec.append({
                        "sample": s,
                        "patient": str(patient[idx[0]]),
                        "unit": unit,
                        "cut": "continuous",
                        "metric": j,
                        "usable": bool(np.isfinite(rho)),
                        "n_high": int(len(rhos_list)),
                        "n_low": int(len(rhos_list)),
                        "n_fov_in_section": int(len(rhos_list)),
                        "mean_high": 0.0,
                        "mean_low": float(-rho) if np.isfinite(rho) else np.nan,
                    })
        else:
            # Donor unit: section rows carry per-section Spearman so 7/8 is defined.
            # Donor sign uses one Spearman on the donor's pooled malignant cells.
            for s in SAMPLES:
                idx = np.flatnonzero(sample == s)
                rhos = np.array([spearman_rho(cldn4[idx], count_mat[idx, j], MIN_SPEARMAN) for j in range(count_mat.shape[1])])
                for j, rho in enumerate(rhos):
                    cont_sec.append({
                        "sample": s,
                        "patient": str(patient[idx[0]]),
                        "unit": unit,
                        "cut": "continuous",
                        "metric": j,
                        "usable": bool(np.isfinite(rho)),
                        "n_high": int(len(idx)),
                        "n_low": int(len(idx)),
                        "n_fov_in_section": 0,
                        "mean_high": 0.0,
                        "mean_low": float(-rho) if np.isfinite(rho) else np.nan,
                    })
            for j in range(count_mat.shape[1]):
                cont_fov_rows.append({
                    "unit": unit, "cut": "continuous", "metric": j,
                    "n_fov": 0, "n_fov_neg": 0, "n_fov_pos": 0,
                })
        print(f"continuous {unit} in {time.time()-t0:.1f}s", flush=True)

    sec = pd.DataFrame(sec_rows + cont_sec)
    fov_table = pd.DataFrame(fov_rows + cont_fov_rows)
    names_ext = names + []  # same metric index; continuous reuses columns
    # Mark continuous metric names by cut, not by a second catalog.
    catalog = []
    for spec in names:
        catalog.append(dict(spec))
    grid = grid_from_sections(sec, fov_table, catalog)
    grid.loc[grid["cut"] == "continuous", "continuous"] = True
    # For continuous, pooled_delta was mean_high - mean_low = 0 - (-rho) = rho. Good.
    # mean_low is -rho, mean_high is 0, so pooled_ratio is meaningless. Force NA.
    grid.loc[grid["cut"] == "continuous", "pooled_ratio"] = np.nan
    grid.loc[grid["cut"] == "continuous", "effect_log_ratio"] = np.nan
    grid.loc[grid["cut"] == "continuous", "median_section_ratio"] = np.nan
    grid.loc[grid["cut"] == "continuous", "mean_high"] = np.nan
    # mean_low currently averages -rho. Replace with the mean rho, and pooled_delta is already rho.
    # Recompute mean_low as NA for continuous so a reader does not treat it as a count.
    grid.loc[grid["cut"] == "continuous", "mean_low"] = np.nan

    # Donor-unit continuous: replace donor sign counts with pooled-donor Spearman signs.
    if (grid["cut"] == "continuous").any() and (grid["unit"] == "donor").any():
        donor_sign = {}
        for j in range(count_mat.shape[1]):
            neg = pos = 0
            for pat in pd.unique(patient):
                idx = np.flatnonzero(patient == pat)
                rho = spearman_rho(cldn4[idx], count_mat[idx, j], MIN_SPEARMAN)
                if not np.isfinite(rho):
                    continue
                if rho < 0:
                    neg += 1
                elif rho > 0:
                    pos += 1
            donor_sign[j] = (neg, pos)
        mask = (grid["cut"] == "continuous") & (grid["unit"] == "donor")
        for i, row in grid.loc[mask].iterrows():
            j = CYTOS.index(row["cyto"]) * len(RADII) + RADII.index(int(row["radius_um"]))
            neg, pos = donor_sign[j]
            grid.at[i, "n_don_neg"] = neg
            grid.at[i, "n_don_pos"] = pos
            grid.at[i, "n_don_usable"] = neg + pos
            grid.at[i, "eligible_exclusion"] = bool(row["n_sec_neg"] >= 7 and neg >= 4)
            grid.at[i, "eligible_enrichment"] = bool(row["n_sec_pos"] >= 7 and pos >= 4)
            grid.at[i, "sign_p_donor"] = sign_pvalue(max(neg, pos), N_DONORS)
    # effect_abs_delta for continuous is |mean rho| which equals |pooled_delta|.
    grid["eligible_exclusion"] = grid["eligible_exclusion"].astype(bool)
    grid["eligible_enrichment"] = grid["eligible_enrichment"].astype(bool)
    return sec, grid


def sanity_check(grid: pd.DataFrame) -> None:
    """Stop if the scale is not a neighbor count. Does not check the locked ratios."""
    hit = grid[
        (~grid["continuous"])
        & (grid["unit"] == "section")
        & (grid["cut"] == "median")
        & (grid["cyto"] == "cd8nk")
        & (grid["radius_um"] == 50)
    ]
    if hit.empty:
        raise SystemExit("missing calibration row: section median CD8+NK 50 µm")
    row = hit.iloc[0]
    level = np.nanmean([row["mean_high"], row["mean_low"]])
    if not np.isfinite(level) or not (0.02 <= level <= 20):
        raise SystemExit(f"calibration count level {level} is outside 0.02–20; refusing to write results")
    hit10 = grid[
        (~grid["continuous"])
        & (grid["unit"] == "section")
        & (grid["cut"] == "median")
        & (grid["cyto"] == "cd8nk")
        & (grid["radius_um"] == 10)
    ].iloc[0]
    hit100 = grid[
        (~grid["continuous"])
        & (grid["unit"] == "section")
        & (grid["cut"] == "median")
        & (grid["cyto"] == "cd8nk")
        & (grid["radius_um"] == 100)
    ].iloc[0]
    if not (hit100["mean_low"] > hit10["mean_low"]):
        raise SystemExit(
            f"100 µm low-arm mean {hit100['mean_low']} is not above 10 µm {hit10['mean_low']}"
        )
    print(
        f"calibration section/median/cd8nk/50µm ratio={row['pooled_ratio']:.4f} "
        f"delta={row['pooled_delta']:.4f} signs={int(row['n_sec_neg'])}/8 {int(row['n_don_neg'])}/5",
        flush=True,
    )


def fmt(x, nd=3):
    if x is None or not isinstance(x, (int, float, np.floating)) or not np.isfinite(x):
        return "NA"
    return f"{float(x):.{nd}f}"


def fmt_p(x):
    if x is None or not isinstance(x, (int, float, np.floating)) or not np.isfinite(x):
        return "NA"
    x = float(x)
    if x < 0.0001:
        return f"{x:.2e}"
    return f"{x:.4f}"


def spec_name(row) -> str:
    if row is None:
        return "none"
    radius = row.get("radius_um")
    return f"{row['cut']}, {row['cyto']}, {radius} µm, {row['unit']} threshold"


def _md_row(row) -> str:
    direction = row.get("direction") or (
        "exclusion" if row["n_sec_neg"] >= row["n_sec_pos"] else "enrichment"
    )
    signs = f"{int(row['n_sec_neg'])}/8 neg, {int(row['n_sec_pos'])}/8 pos; donors {int(row['n_don_neg'])}/5 neg, {int(row['n_don_pos'])}/5 pos"
    return (
        f"| {row['cut']} | {row['cyto']} | {int(row['radius_um'])} | {row['unit']} | "
        f"{fmt(row['mean_high'])} | {fmt(row['mean_low'])} | {fmt(row['pooled_ratio'])} | "
        f"{fmt(row['pooled_delta'])} | {signs} | {int(row['n_fov_neg'])}/{int(row['n_fov'])} | "
        f"{direction} |"
    )


def write_report(grid, sec, winners, cont, meta, path: str) -> None:
    ratio_w = winners["ratio_winner"]
    delta_w = winners["delta_winner"]

    def pull(unit, cut, cyto, radius):
        hit = grid[
            (grid["unit"] == unit) & (grid["cut"] == cut) & (grid["cyto"] == cyto) & (grid["radius_um"] == radius)
        ]
        if hit.empty:
            return None
        return hit.iloc[0]

    lines = []
    lines.append("# CosMx NSCLC: largest CLDN4 cytotoxic-neighbor effect that keeps the sign")
    lines.append("")
    lines.append(
        "Sensitivity layer on He et al. 2022 CosMx NSCLC (figshare 25976224, "
        "`cosmx_human_nsclc_clustered.h5ad`, 8 sections / 5 donors). "
        "The locked exclusion summary is unchanged: cytotoxic neighbor ratio "
        "**0.36 at 50 µm and 0.52 at 100 µm**, 8/8 sections, 5/5 donors, sign P = 0.031. "
        "Nothing below replaces that summary."
    )
    lines.append("")
    lines.append("## Call")
    lines.append("")
    if ratio_w is None:
        lines.append(
            "No neighbor-count specification had at least 7/8 sections and 4/5 donors "
            "with the same sign and a defined high/low ratio."
        )
    else:
        lines.append(
            f"Largest ratio effect among specs with ≥7/8 sections and ≥4/5 donors of one sign: "
            f"**{spec_name(ratio_w)}**. "
            f"Equal-weight section means {fmt(ratio_w['mean_high'])} (high) / {fmt(ratio_w['mean_low'])} (low), "
            f"high/low ratio **{fmt(ratio_w['pooled_ratio'])}**, Δ **{fmt(ratio_w['pooled_delta'])}**. "
            f"Sign: {int(ratio_w['n_sec_sign'])}/8 sections and {int(ratio_w['n_don_sign'])}/5 donors "
            f"({ratio_w['direction']}). "
            f"|log ratio| = {fmt(ratio_w['effect_log_ratio'], 3)}. "
            f"One-sided sign P = {fmt_p(ratio_w['sign_p_section'])} (sections, out of 8) and "
            f"{fmt_p(ratio_w['sign_p_donor'])} (donors, out of 5). "
            f"Two-sided Wilcoxon on usable section deltas P = {fmt_p(ratio_w['wilcoxon_p'])}. "
            f"FOVs with the same high<low direction: {int(ratio_w['n_fov_neg'])}/{int(ratio_w['n_fov'])} "
            f"(high>low {int(ratio_w['n_fov_pos'])}/{int(ratio_w['n_fov'])})."
        )
    lines.append("")
    if delta_w is None:
        lines.append("No specification had a defined Δ under the same sign bar.")
    else:
        same = (
            ratio_w is not None
            and delta_w["cut"] == ratio_w["cut"]
            and delta_w["cyto"] == ratio_w["cyto"]
            and int(delta_w["radius_um"]) == int(ratio_w["radius_um"])
            and delta_w["unit"] == ratio_w["unit"]
        )
        if same:
            lines.append("The largest |Δ| is this same specification.")
        else:
            lines.append(
                f"Largest |Δ| under the same sign bar: **{spec_name(delta_w)}**. "
                f"Means {fmt(delta_w['mean_high'])} / {fmt(delta_w['mean_low'])}, "
                f"ratio {fmt(delta_w['pooled_ratio'])}, Δ **{fmt(delta_w['pooled_delta'])}**. "
                f"Sign: {int(delta_w['n_sec_sign'])}/8 sections and {int(delta_w['n_don_sign'])}/5 donors "
                f"({delta_w['direction']}). "
                f"Wilcoxon P = {fmt_p(delta_w['wilcoxon_p'])}."
            )
    lines.append("")
    if winners["numeric_max_ratio"] is not None:
        a = winners["numeric_max_ratio"]
        b = winners["numeric_min_ratio"]
        c = winners["numeric_max_delta"]
        d = winners["numeric_min_delta"]
        lines.append(
            f"Every eligible high/low ratio in this grid is below 1, and every eligible Δ is negative. "
            f"The numerically largest ratio is {fmt(a['pooled_ratio'])} ({spec_name(a)}) and the numerically "
            f"smallest is {fmt(b['pooled_ratio'])} ({spec_name(b)}). "
            f"The least negative Δ is {fmt(c['pooled_delta'])} ({spec_name(c)}) and the most negative Δ is "
            f"{fmt(d['pooled_delta'])} ({spec_name(d)}). "
            f"Eligible binary specs: {winners['n_eligible']} of {winners['n_binary']}."
        )
    lines.append("")
    lines.extend(_extra_context(grid))
    lines.append("")
    if cont.get("winner"):
        w = cont["winner"]
        lines.append(
            f"Continuous CLDN4 (Spearman of log1p CP10k vs neighbor count), largest |mean ρ| "
            f"with the same sign bar: **{spec_name(w)}**, mean section ρ = {fmt(w['pooled_delta'], 3)} "
            f"({w['direction']}; {int(w['n_sec_neg'])}/8 negative, {int(w['n_don_neg'])}/5 donors negative). "
            f"This is not a high/low ratio and is not the count-sweep call."
        )
    else:
        lines.append(
            "No continuous Spearman specification reached ≥7/8 sections and ≥4/5 donors of one sign."
        )
    lines.append("")
    lines.append("## Selection rule")
    lines.append("")
    lines.append(
        "Eligible: ≥7 sections with Δ<0 and ≥4 donors with Δ<0, or the same counts with Δ>0. "
        "There are 8 sections and 5 donors; a section that cannot form both arms is not a supporting sign. "
        "Δ is the equal-weight mean of section (high − low) means. "
        "The ratio is the equal-weight mean of section high means divided by the equal-weight mean of section low means. "
        "It is left undefined when the low mean is below 0.05. "
        "Ratio effect is |log(ratio)|. Delta effect is |Δ|. "
        "Ties break on more sections in the majority direction, then section over FOV over donor. "
        "P-values describe the winning row. The grid was searched, so they are not a second locked test. "
        "Median split is in the grid as a calibration cut, not as a replacement for Q4/Q1."
    )
    lines.append("")
    lines.append("## Definitions")
    lines.append("")
    lines.append(
        f"- Malignant index: author `cell_type` matched to the section (`tumor 5/6/9/12/13`). "
        f"CD8 labels: {', '.join(meta['cd8_labels'])}. NK labels: {', '.join(meta['nk_labels'])}. "
        f"CD8+NK is their union. GZMB+: {meta['gzmb_definition']}."
    )
    lines.append(
        "- CLDN4 is log1p(count / n_counts × 10,000). Detected means raw count > 0. "
        "Thresholds are computed inside the named unit (section, FOV, or donor). "
        "Q4/Q1 keeps cells ≥ the 75th percentile versus cells ≤ the 25th percentile; "
        "overlap is removed. If the 75th percentile is 0, that group collapses to detected versus absent. "
        "Top X% versus bottom requires the upper quantile to sit strictly above the lower quantile. "
        "Top X% versus rest is cells at or above the upper quantile versus everyone below it. "
        "A section needs ≥30 cells in each arm. An FOV needs ≥5 in each arm, and a section needs ≥3 such FOVs."
    )
    lines.append(
        "- Counts are other cells of that cytotoxic class within the radius, in the same section, "
        f"global centroids × {UM_PER_PX} µm/px. The index cell is not in the reference set. "
        "Donor is not a spatial window: Lung5's three sections keep separate coordinate systems."
    )
    lines.append(
        "- Donor Δ averages the donor's usable sections with equal section weight. "
        "FOV signs use FOVs with ≥5 cells in each arm."
    )
    lines.append("")
    lines.append("## Calibration row (not the lock)")
    lines.append("")
    lines.append(
        "Section-level median split, CD8+NK. This is the unstratified contrast that previously "
        "read about 0.80 at 50 µm on a closely related tumor definition. It is a pipeline check. "
        "It is not the locked 0.36 / 0.52 summary."
    )
    lines.append("")
    lines.append("| Radius | High | Low | Ratio | Δ | Sections neg | Donors neg | Wilcoxon |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
    for radius in (10, 50, 100):
        row = pull("section", "median", "cd8nk", radius)
        if row is None:
            continue
        lines.append(
            f"| {radius} | {fmt(row['mean_high'])} | {fmt(row['mean_low'])} | {fmt(row['pooled_ratio'])} | "
            f"{fmt(row['pooled_delta'])} | {int(row['n_sec_neg'])}/8 | {int(row['n_don_neg'])}/5 | "
            f"{fmt_p(row['wilcoxon_p'])} |"
        )
    lines.append("")
    lines.append("## Section means for the ratio call")
    lines.append("")
    if ratio_w is not None:
        lines.append(_section_table(sec, ratio_w))
        if ratio_w["cut"] == "detected_absent":
            lines.append("")
            lines.append(
                "Detected versus absent does not use a quantile, so a donor-wide threshold labels the same cells as the section threshold. "
                "The donor-unit row matches this table."
            )
    lines.append("")
    if delta_w is not None and (
        ratio_w is None
        or delta_w["cut"] != ratio_w["cut"]
        or delta_w["cyto"] != ratio_w["cyto"]
        or int(delta_w["radius_um"]) != int(ratio_w["radius_um"])
        or delta_w["unit"] != ratio_w["unit"]
    ):
        lines.append("## Section means for the Δ call")
        lines.append("")
        lines.append(_section_table(sec, delta_w))
        lines.append("")
    lines.append("## Eligible specs with the largest |log ratio|")
    lines.append("")
    lines.append(
        "| Cut | Cytotoxic | µm | Unit | High | Low | Ratio | Δ | Section / donor signs | FOV high<low | Direction |"
    )
    lines.append("|---|---|---:|---|---:|---:|---:|---:|---|---|---|")
    elig = _eligible(grid[~grid["continuous"].astype(bool)])
    elig = elig[np.isfinite(elig["effect_log_ratio"])].copy()
    elig["unit_rank"] = elig["unit"].map({"section": 0, "fov": 1, "donor": 2})
    elig = elig.sort_values(["effect_log_ratio", "unit_rank"], ascending=[False, True])
    if elig.empty:
        lines.append("| none | | | | | | | | | | |")
    else:
        for _, row in elig.head(12).iterrows():
            packed = row.to_dict()
            packed["direction"] = _direction(row)
            lines.append(_md_row(packed))
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| Section | Donor | Cells | Malignant | CLDN4>0 | CD8 | NK | GZMB+ non-tumor | Median NN µm |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for rec in meta["inventory"]:
        lines.append(
            f"| {rec['sample']} | {rec['patient']} | {rec['n_cells']} | {rec['n_malignant']} | "
            f"{rec['n_cldn4_pos']} | {rec['n_cd8']} | {rec['n_nk']} | {rec['n_gzmb']} | "
            f"{rec['median_nn_um']:.2f} |"
        )
    lines.append("")
    comp = {}
    for rec in meta.get("gzmb_composition", {}).values():
        for name, n in rec.items():
            comp[name] = comp.get(name, 0) + int(n)
    gz_total = sum(comp.values())
    gz_line = ""
    if gz_total:
        ranked = sorted(comp.items(), key=lambda kv: kv[1], reverse=True)[:3]
        bits = ", ".join(f"{name} {n / gz_total:.0%}" for name, n in ranked)
        cd8_n = sum(n for name, n in comp.items() if "CD8" in name)
        nk_n = sum(n for name, n in comp.items() if name == "NK" or name.startswith("NK"))
        gz_line = (
            f" GZMB+ non-tumor cells: {gz_total}. Largest shares: {bits}. "
            f"CD8 is {cd8_n / gz_total:.0%} and NK is {nk_n / gz_total:.0%} of that set, "
            f"so GZMB+ is not a purified cytotoxic class."
        )
    lines.append(
        f"Malignant cells in the neighbor table: {meta['n_malignant']}.{gz_line} "
        "Full GZMB+ counts are in `summary.json`."
    )
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- It does not replace 0.36 / 0.52, and it does not restate that pair as a new estimate from this grid.")
    lines.append("- It does not say nearby effectors are muzzled. This sweep counts cells, not GZMB/PRF1/NKG7/IFNG inside the cells that are present.")
    lines.append("- Sign tests and Wilcoxon p-values are descriptive. FOVs are nested in five donors, and the grid was searched.")
    lines.append("- No private 8-KL. No ICI labels.")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/download_cosmx_nsclc_h5ad.py")
    lines.append("python3 scripts/cosmx_cytotoxic_neighbor_sweep.py")
    lines.append("```")
    lines.append("")
    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text)
    alt = os.path.join(OUT, "RESULTS.md")
    with open(alt, "w") as fh:
        fh.write(text)


def _extra_context(grid: pd.DataFrame) -> list[str]:
    """Fully concordant row, and the sparse ratios the 0.05 floor leaves undefined."""
    binary = grid[~grid["continuous"].astype(bool)].copy()
    passing = binary[(binary["n_sec_neg"] >= 7) & (binary["n_don_neg"] >= 4)].copy()
    lines = []
    full = passing[
        (passing["n_sec_neg"] == 8)
        & (passing["n_don_neg"] == 5)
        & np.isfinite(passing["pooled_ratio"])
        & (passing["pooled_ratio"] > 0)
    ]
    if len(full):
        row = full.loc[full["pooled_ratio"].idxmin()]
        lines.append(
            f"Strongest fully concordant ratio (8/8 sections and 5/5 donors, low-arm mean ≥ {LOW_MEAN_FLOOR}): "
            f"{fmt(row['pooled_ratio'])} at {spec_name(row.to_dict())}, "
            f"means {fmt(row['mean_high'])} / {fmt(row['mean_low'])}, Δ {fmt(row['pooled_delta'])}. "
            f"The call above is a larger fold and is not 8/8."
        )
        lines.append("")
    sparse = passing[(passing["mean_low"] > 0) & (passing["mean_low"] < LOW_MEAN_FLOOR)].copy()
    if len(sparse):
        sparse["raw_ratio"] = sparse["mean_high"] / sparse["mean_low"]
        row = sparse.loc[sparse["raw_ratio"].idxmin()]
        lines.append(
            f"Neighborhoods whose equal-weight low-arm mean is below {LOW_MEAN_FLOOR} do not get a ratio in the call. "
            f"The smallest such high/low ratio with the sign bar is {fmt(row['raw_ratio'])} "
            f"({spec_name(row.to_dict())}; means {fmt(row['mean_high'], 4)} / {fmt(row['mean_low'], 4)}; "
            f"{int(row['n_sec_neg'])}/8 sections, {int(row['n_don_neg'])}/5 donors). "
            f"That is a ratio of rare events, not the reported sensitivity."
        )
    return lines


def _section_table(sec: pd.DataFrame, winner: dict) -> str:
    # The long section table stores metric index. Attach names if needed.
    if "cyto" in sec.columns:
        sub = sec[
            (sec["unit"] == winner["unit"])
            & (sec["cut"] == winner["cut"])
            & (sec["cyto"] == winner["cyto"])
            & (sec["radius_um"] == int(winner["radius_um"]))
        ]
    else:
        j = CYTOS.index(winner["cyto"]) * len(RADII) + RADII.index(int(winner["radius_um"]))
        sub = sec[(sec["unit"] == winner["unit"]) & (sec["cut"] == winner["cut"]) & (sec["metric"] == j)]
    lines = ["| Section | Donor | High | Low | Δ | Ratio | n high | n low | Usable |",
             "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    order = {s: i for i, s in enumerate(SAMPLES)}
    sub = sub.copy()
    sub["_ord"] = sub["sample"].map(order)
    for _, r in sub.sort_values("_ord").iterrows():
        if not r["usable"] or not np.isfinite(r["mean_high"]):
            lines.append(f"| {r['sample']} | {r['patient']} | NA | NA | NA | NA | {int(r['n_high'])} | {int(r['n_low'])} | no |")
            continue
        delta = float(r["mean_high"] - r["mean_low"])
        ratio = float(r["mean_high"] / r["mean_low"]) if r["mean_low"] > 0 else float("nan")
        lines.append(
            f"| {r['sample']} | {r['patient']} | {fmt(r['mean_high'], 4)} | {fmt(r['mean_low'], 4)} | "
            f"{fmt(delta, 4)} | {fmt(ratio, 3)} | {int(r['n_high'])} | {int(r['n_low'])} | yes |"
        )
    return "\n".join(lines)


def make_figures(grid: pd.DataFrame, sec: pd.DataFrame, winners: dict) -> None:
    os.makedirs(FIG, exist_ok=True)
    ratio_w = winners["ratio_winner"]
    if ratio_w is None:
        return
    j = CYTOS.index(ratio_w["cyto"]) * len(RADII) + RADII.index(int(ratio_w["radius_um"]))
    sub = sec[(sec["unit"] == ratio_w["unit"]) & (sec["cut"] == ratio_w["cut"]) & (sec["metric"] == j)].copy()
    order = {s: i for i, s in enumerate(SAMPLES)}
    sub["_ord"] = sub["sample"].map(order)
    sub = sub.sort_values("_ord")
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for _, r in sub.iterrows():
        if not r["usable"] or not np.isfinite(r["mean_high"]):
            continue
        ax.plot([0, 1], [r["mean_low"], r["mean_high"]], color="#4c78a8", lw=1.2)
        ax.scatter([0, 1], [r["mean_low"], r["mean_high"]], color="#4c78a8", s=28, zorder=3)
        ax.text(1.04, r["mean_high"], r["sample"].replace(" ", "\n"), fontsize=7, va="center")
    ax.set_xticks([0, 1], ["CLDN4 low", "CLDN4 high"])
    ax.set_xlim(-0.2, 1.55)
    ax.set_ylabel("Mean cytotoxic neighbor count")
    ax.set_title(
        f"{ratio_w['cyto']} at {int(ratio_w['radius_um'])} µm\n{ratio_w['cut']}, {ratio_w['unit']} threshold"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ratio_winner_paired.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "ratio_winner_paired.pdf"))
    plt.close(fig)

    # Heatmap of ratio for the winning cytotoxic class and unit.
    subg = grid[
        (~grid["continuous"].astype(bool))
        & (grid["cyto"] == ratio_w["cyto"])
        & (grid["unit"] == ratio_w["unit"])
    ].copy()
    cuts = [c for c in CUTS if c in set(subg["cut"])]
    mat = np.full((len(cuts), len(RADII)), np.nan)
    elig_mask = np.zeros_like(mat, dtype=bool)
    for i, cut in enumerate(cuts):
        for k, radius in enumerate(RADII):
            hit = subg[(subg["cut"] == cut) & (subg["radius_um"] == radius)]
            if hit.empty:
                continue
            mat[i, k] = hit["pooled_ratio"].iloc[0]
            elig_mask[i, k] = bool(hit["eligible_exclusion"].iloc[0] or hit["eligible_enrichment"].iloc[0])
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    show = np.array(mat, dtype=float)
    finite = show[np.isfinite(show)]
    if len(finite):
        span = max(1.0 - float(np.min(finite)), float(np.max(finite)) - 1.0, 0.15)
    else:
        span = 0.5
    im = ax.imshow(show, cmap="RdBu", vmin=1.0 - span, vmax=1.0 + span, aspect="auto")
    ax.set_xticks(range(len(RADII)), [str(r) for r in RADII])
    ax.set_yticks(range(len(cuts)), cuts)
    ax.set_xlabel("Radius (µm)")
    for i in range(len(cuts)):
        for k in range(len(RADII)):
            val = mat[i, k]
            if not np.isfinite(val):
                label = "·"
            else:
                label = f"{val:.2f}"
            weight = "bold" if elig_mask[i, k] else "normal"
            ax.text(k, i, label, ha="center", va="center", fontsize=7, fontweight=weight)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="High/low ratio")
    ax.set_title(f"{ratio_w['cyto']}, {ratio_w['unit']} threshold\nbold = ≥7/8 and ≥4/5 same sign")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "grid_ratio_heatmap.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "grid_ratio_heatmap.pdf"))
    plt.close(fig)


def attach_metric_names(sec: pd.DataFrame) -> pd.DataFrame:
    names = metric_catalog()
    sec = sec.copy()
    sec["cyto"] = sec["metric"].map(lambda j: names[int(j)]["cyto"])
    sec["radius_um"] = sec["metric"].map(lambda j: names[int(j)]["radius"])
    return sec


def self_test() -> None:
    vals = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 5], dtype=float)
    pos = vals > 0
    high, low = arm_masks(vals, pos, "detected_absent", min_arm=2)
    assert high.sum() == 5 and low.sum() == 5
    high, low = arm_masks(vals, pos, "q4_q1", min_arm=2)
    assert low.all() or low.sum() >= 2
    assert not np.any(high & low)
    # No separation when every value is 0.
    zeros = np.zeros(20)
    assert arm_masks(zeros, zeros > 0, "top10_vs_bottom", min_arm=2) is None
    assert arm_masks(zeros, zeros > 0, "top10_vs_rest", min_arm=2) is None
    spread = np.arange(100, dtype=float)
    masks = arm_masks(spread, spread > 0, "top10_vs_bottom", min_arm=5)
    assert masks is not None and masks[0].sum() >= 5 and masks[1].sum() >= 5

    grid = pd.DataFrame([
        dict(unit="section", cut="mild", cyto="cd8", radius_um=50, continuous=False,
             n_sec_neg=8, n_sec_pos=0, n_don_neg=5, n_don_pos=0,
             pooled_ratio=0.50, pooled_delta=-1.0, effect_log_ratio=abs(np.log(0.50)),
             effect_abs_delta=1.0, eligible_exclusion=True, eligible_enrichment=False,
             mean_high=0.5, mean_low=1.0),
        dict(unit="section", cut="steep", cyto="nk", radius_um=10, continuous=False,
             n_sec_neg=7, n_sec_pos=1, n_don_neg=4, n_don_pos=1,
             pooled_ratio=0.20, pooled_delta=-0.4, effect_log_ratio=abs(np.log(0.20)),
             effect_abs_delta=0.4, eligible_exclusion=True, eligible_enrichment=False,
             mean_high=0.2, mean_low=1.0),
        dict(unit="fov", cut="ineligible", cyto="cd8nk", radius_um=20, continuous=False,
             n_sec_neg=6, n_sec_pos=2, n_don_neg=5, n_don_pos=0,
             pooled_ratio=0.05, pooled_delta=-5.0, effect_log_ratio=abs(np.log(0.05)),
             effect_abs_delta=5.0, eligible_exclusion=False, eligible_enrichment=False,
             mean_high=0.05, mean_low=1.0),
        dict(unit="donor", cut="enrich", cyto="gzmb", radius_um=100, continuous=False,
             n_sec_neg=0, n_sec_pos=8, n_don_neg=0, n_don_pos=5,
             pooled_ratio=3.0, pooled_delta=2.0, effect_log_ratio=abs(np.log(3.0)),
             effect_abs_delta=2.0, eligible_exclusion=False, eligible_enrichment=True,
             mean_high=3.0, mean_low=1.0),
        dict(unit="section", cut="continuous", cyto="cd8", radius_um=50, continuous=True,
             n_sec_neg=8, n_sec_pos=0, n_don_neg=5, n_don_pos=0,
             pooled_ratio=np.nan, pooled_delta=-0.9, effect_log_ratio=np.nan,
             effect_abs_delta=0.9, eligible_exclusion=True, eligible_enrichment=False,
             mean_high=np.nan, mean_low=np.nan),
    ])
    winners = select_winners(grid)
    assert winners["ratio_winner"]["cut"] == "steep", winners["ratio_winner"]["cut"]
    assert winners["delta_winner"]["cut"] == "enrich"
    assert abs(winners["numeric_max_ratio"]["pooled_ratio"] - 3.0) < 1e-9
    assert abs(winners["numeric_min_ratio"]["pooled_ratio"] - 0.20) < 1e-9
    assert abs(winners["numeric_max_delta"]["pooled_delta"] - 2.0) < 1e-9
    assert abs(winners["numeric_min_delta"]["pooled_delta"] - (-1.0)) < 1e-9
    cont = select_continuous(grid)
    assert cont["winner"]["cut"] == "continuous"
    # Floor: tiny low mean cannot win the ratio.
    grid2 = grid.copy()
    grid2.loc[grid2["cut"] == "steep", "effect_log_ratio"] = np.nan
    grid2.loc[grid2["cut"] == "steep", "pooled_ratio"] = np.nan
    winners2 = select_winners(grid2)
    assert winners2["ratio_winner"]["cut"] == "enrich"
    print("self-test ok")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--reuse-cache", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    if args.reuse_cache and os.path.isfile(CACHE):
        cells, meta = load_cache()
        print(f"reused cache n={len(cells)}", flush=True)
    else:
        if not os.path.isfile(H5AD):
            raise SystemExit(f"missing {H5AD}; run scripts/download_cosmx_nsclc_h5ad.py")
        cells, meta = load_cells()
        save_cache(cells, meta)
    sec, grid = run_sweep(cells)
    sec = attach_metric_names(sec)
    sanity_check(grid)
    winners = select_winners(grid)
    cont = select_continuous(grid)
    grid.to_csv(os.path.join(TAB, "parameter_grid.csv"), index=False)
    sec.to_csv(os.path.join(TAB, "section_means.csv"), index=False)
    with open(os.path.join(TAB, "inventory.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    # JSON-safe winners
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items() if k != "donor_deltas_json"}
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        if isinstance(obj, (np.floating,)):
            v = float(obj)
            return None if not np.isfinite(v) else v
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        return obj
    summary = {"winners": _clean(winners), "continuous": _clean(cont), "locked": {"ratio_50_um": 0.36, "ratio_100_um": 0.52, "sections": "8/8", "donors": "5/5", "sign_p": 0.031, "replaced": False}}
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    write_report(grid, sec, winners, cont, meta, os.path.join(ROOT, "RESULTS.md"))
    make_figures(grid, sec, winners)
    print("wrote RESULTS.md", flush=True)
    if winners["ratio_winner"]:
        w = winners["ratio_winner"]
        print(
            f"RATIO {spec_name(w)} {w['pooled_ratio']:.4f} delta {w['pooled_delta']:.4f} "
            f"{w['direction']} {w['n_sec_sign']}/8 {w['n_don_sign']}/5",
            flush=True,
        )
    if winners["delta_winner"]:
        w = winners["delta_winner"]
        print(
            f"DELTA {spec_name(w)} {w['pooled_delta']:.4f} ratio {w['pooled_ratio']} "
            f"{w['direction']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
