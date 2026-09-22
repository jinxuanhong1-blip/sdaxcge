#!/usr/bin/env python3
"""Smallest TACSTD2-high vs low immune and CD8+NK ratio that stays 8/8 and 5/5.

He et al. 2022 CosMx NSCLC, figshare 25976224,
cosmx_human_nsclc_clustered.h5ad. Index cells are the author tumor label
matched to the section. Coordinates are obsm/spatial times 0.18 µm/pixel.
Trees are built per section. The index cell is excluded from the ball.

Immune labels and the CD8+NK labels are fixed to the set used for the
published TACSTD2 neighbor note. They are not swapped after seeing which
label is coldest. The grid below (radii, cuts, estimands) is fixed in this
file. Ranking happens only after every spec is scored.

A spec is eligible only when the high arm is strictly lower in every
section and every patient, each arm has at least 30 finite malignant cells
in every section, the high-arm mean is above 0 in every section, and the
low arm clears the estimand floor below. Patient values are the unweighted
mean of that patient's section means.

Fraction (frac) drops cells with no neighbor, matching the published immune
fraction. Count keeps those cells as zero. The zero-filled fraction (frac0)
is eligible only when the dropped-empty fraction for the same masks is also
strictly lower in 8/8 and 5/5 and clears the fraction floor. That stops an
empty-neighborhood effect from standing in for a compositional one.

Calibration against the published median and detected contrasts is an
assert. If it fails, this script exits and does not write a headline.
"""

from __future__ import annotations

import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_tacstd2_max_effect")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

UM_PER_PX = 0.18
MIN_N = 30
FRAC_FLOOR = 0.005
COUNT_FLOOR = 0.001
COUNT_MIN_LOW_EVENTS = 20.0
COUNT_MIN_HIGH_EVENTS = 5.0

RADII = tuple(range(4, 26)) + (28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100)

SAMPLES = [
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
]
PATIENT = {
    "LUAD-5 R1": "Lung5",
    "LUAD-5 R2": "Lung5",
    "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9",
    "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}
TUMOR = {
    "LUAD-5 R1": "tumor 5",
    "LUAD-5 R2": "tumor 5",
    "LUAD-5 R3": "tumor 5",
    "LUSC-6": "tumor 6",
    "LUAD-9 R1": "tumor 9",
    "LUAD-9 R2": "tumor 9",
    "LUAD-12": "tumor 12",
    "LUAD-13": "tumor 13",
}
BROAD = (
    "B-cell",
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
    "Treg",
    "mDC",
    "macrophage",
    "mast",
    "monocyte",
    "neutrophil",
    "pDC",
    "plasmablast",
)
CD8 = ("T CD8 memory", "T CD8 naive")
CD8NK = ("T CD8 memory", "T CD8 naive", "NK")
CLASSES = {
    "immune": BROAD,
    "cd8nk": CD8NK,
    "cd8": CD8,
}
PRIMARY_CLASSES = ("immune", "cd8nk")
ESTIMANDS = ("frac", "count", "frac0")
PATIENT_COLOR = {
    "Lung5": "#1b9e77",
    "Lung6": "#d95f02",
    "Lung9": "#7570b3",
    "Lung12": "#e7298a",
    "Lung13": "#66a61e",
}
SHORT = {
    "LUAD-5 R1": "5-R1",
    "LUAD-5 R2": "5-R2",
    "LUAD-5 R3": "5-R3",
    "LUSC-6": "6",
    "LUAD-9 R1": "9-R1",
    "LUAD-9 R2": "9-R2",
    "LUAD-12": "12",
    "LUAD-13": "13",
}

# Published equal-weight means this pipeline has to reproduce before ranking.
CALIBRATION = [
    ("log_median", "immune", "frac", 10, 0.0389631200782216, 0.07510678178070865, 0.5187696657271683),
    ("log_median", "immune", "frac", 20, 0.06137493243867926, 0.09551769859878653, 0.6425503685602719),
    ("log_median", "cd8nk", "count", 10, 0.004342177231251601, 0.009534437192328247, 0.4554204032877236),
    ("log_median", "cd8nk", "count", 20, 0.03940746073144875, 0.05880262090783151, 0.6701650389566283),
    ("detected", "immune", "count", 10, 0.03486319633932945, 0.08864670741665898, 0.39328247326169463),
    ("log_median", "immune", "count", 50, 5.300324394399313, 6.233482901632886, 0.8502990187092471),
    ("detected", "immune", "frac", 10, 0.03821405850050407, 0.07798484689242138, 0.490019023224084),
    ("log_q4q1", "immune", "frac", 10, 0.043094840611556315, 0.07798484689242138, 0.5526053115293645),
]
PUBLISHED_POS = {
    "LUAD-5 R1": 0.6626114256881762,
    "LUAD-5 R2": 0.6336943172937367,
    "LUAD-5 R3": 0.5907632290945699,
    "LUSC-6": 0.3724545289745604,
    "LUAD-9 R1": 0.26844895018526144,
    "LUAD-9 R2": 0.2682343269109153,
    "LUAD-12": 0.32402693381507636,
    "LUAD-13": 0.3477910103726469,
}


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def stable_seed(*parts) -> int:
    acc = 2166136261
    for part in parts:
        for ch in str(part):
            acc ^= ord(ch)
            acc = (acc * 16777619) & 0xFFFFFFFF
    return int(acc)


def ball_counts(tree: cKDTree, xy: np.ndarray, radius: float) -> np.ndarray:
    if tree.n == 0 or len(xy) == 0:
        return np.zeros(len(xy), dtype=np.float64)
    out = tree.query_ball_point(xy, r=float(radius), return_length=True, workers=-1)
    return np.asarray(out, dtype=np.float64)


def cut_catalog():
    """Fixed cut list. Names are the only switch in make_masks."""
    cuts = ["detected", "ge2", "ge3", "ge4", "ge5", "ge6", "ge8", "ge10", "ge12", "ge15", "ge20"]
    cuts += [f"pos_q{q}" for q in (50, 60, 70, 75, 80, 90, 95)]
    cuts += ["log_median", "log_q4q1"]
    cuts += ["log_q90_vs_q10", "log_q80_vs_q20", "log_q75_vs_q25", "log_q70_vs_q30"]
    cuts += ["log_q90_vs_rest", "log_q80_vs_rest", "log_q75_vs_rest", "log_q70_vs_rest"]
    cuts += ["log_pos_q75", "log_pos_q90", "log_pos_q95"]
    cuts += ["rank10", "rank20", "rank30"]
    cuts += ["fov_median"]
    return tuple(cuts)


CUTS = cut_catalog()


def make_masks(name: str, raw: np.ndarray, logv: np.ndarray, fov: np.ndarray, sample: str):
    """Return high, low, threshold. Cells in neither arm are the middle."""
    n = len(raw)
    high = np.zeros(n, dtype=bool)
    low = np.zeros(n, dtype=bool)
    thr = np.nan
    if name == "detected":
        high = raw > 0
        low = ~high
        thr = 0.0
    elif name.startswith("ge"):
        k = float(name[2:])
        high = raw >= k
        low = raw <= 0
        high = high & ~low
        thr = k
    elif name.startswith("pos_q"):
        q = int(name.split("q")[1]) / 100.0
        pos = raw > 0
        if int(pos.sum()) < MIN_N:
            return high, low, thr
        thr = float(np.quantile(raw[pos].astype(np.float64), q))
        high = raw >= thr
        low = raw <= 0
        high = high & ~low
    elif name == "log_median":
        thr = float(np.median(logv))
        high = logv > thr
        low = logv <= thr
    elif name == "log_q4q1":
        q1, q3 = np.quantile(logv, [0.25, 0.75])
        thr = float(q3)
        low = logv <= q1
        high = (logv >= q3) & ~low
    elif name.startswith("log_q") and "_vs_q" in name:
        # log_q90_vs_q10
        left, right = name[len("log_q") :].split("_vs_q")
        hi_q = int(left) / 100.0
        lo_q = int(right) / 100.0
        q_lo, q_hi = np.quantile(logv, [lo_q, hi_q])
        thr = float(q_hi)
        low = logv <= q_lo
        high = (logv >= q_hi) & ~low
    elif name.startswith("log_q") and name.endswith("_vs_rest"):
        hi_q = int(name[len("log_q") :].split("_")[0]) / 100.0
        thr = float(np.quantile(logv, hi_q))
        high = logv >= thr
        low = ~high
    elif name.startswith("log_pos_q"):
        q = int(name.split("q")[1]) / 100.0
        pos = logv > 0
        if int(pos.sum()) < MIN_N:
            return high, low, thr
        thr = float(np.quantile(logv[pos].astype(np.float64), q))
        high = logv >= thr
        low = logv <= 0
        high = high & ~low
    elif name.startswith("rank"):
        frac = int(name[4:]) / 100.0
        n_arm = int(np.floor(frac * n))
        if n_arm < MIN_N or n_arm * 2 > n:
            return high, low, thr
        rng = np.random.default_rng(stable_seed(sample, name))
        jitter = rng.random(n) * 1e-8
        order = np.argsort(logv + jitter, kind="mergesort")
        low[order[:n_arm]] = True
        high[order[-n_arm:]] = True
        thr = float(logv[high].min()) if high.any() else np.nan
    elif name == "fov_median":
        for fv in np.unique(fov):
            m = fov == fv
            if int(m.sum()) < 2:
                continue
            med = float(np.median(logv[m]))
            high[m] = logv[m] > med
            low[m] = logv[m] <= med
        thr = np.nan
    else:
        raise KeyError(name)
    if np.any(high & low):
        raise RuntimeError(f"{name} produced overlapping arms")
    return high, low, thr


def load_sections():
    import h5py

    t0 = time.time()
    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gi = genes.index("TACSTD2")
    cats = _decode(f["obs/cell_type/categories"][:])
    missing = sorted(set(BROAD) - set(cats))
    if missing:
        raise SystemExit(f"immune labels missing: {missing}")
    codes = f["obs/cell_type/codes"][:]
    cell_type = np.asarray(cats, dtype=object)[codes]
    sample_cats = _decode(f["obs/sample/categories"][:])
    sample = np.asarray(sample_cats, dtype=object)[f["obs/sample/codes"][:]]
    fov = f["obs/fov"][:].astype(np.int32)
    n_counts = f["obs/n_counts"][:].astype(np.float64)
    xy = f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    f.close()
    hits = np.flatnonzero(indices == gi)
    rows = np.searchsorted(indptr, hits, side="right") - 1
    raw = np.zeros(cell_type.shape[0], dtype=np.float64)
    np.add.at(raw, rows, np.asarray(data[hits], dtype=np.float64))
    del indices, data, indptr, hits, rows
    lib = n_counts.copy()
    lib[lib <= 0] = np.nan
    logv = np.nan_to_num(np.log1p(raw / lib * 1e4), nan=0.0)
    print(f"loaded TACSTD2 in {time.time() - t0:.1f}s", flush=True)

    sections = []
    inventory = []
    for s in SAMPLES:
        t1 = time.time()
        m = sample == s
        xy_s = np.ascontiguousarray(xy[m])
        ct = cell_type[m]
        mal = ct == TUMOR[s]
        if int(mal.sum()) < MIN_N * 2:
            raise SystemExit(f"{s} has too few matched tumor cells")
        rng = np.random.default_rng(0)
        take = rng.choice(xy_s.shape[0], size=min(6000, xy_s.shape[0]), replace=False)
        dist, _ = cKDTree(xy_s).query(xy_s[take], k=2, workers=-1)
        med_nn = float(np.median(dist[:, 1]))
        if not (3.0 <= med_nn <= 40.0):
            raise SystemExit(f"{s}: median NN {med_nn:.3f} µm is outside 3–40; check 0.18 µm/px")
        tum_xy = np.ascontiguousarray(xy_s[mal])
        trees = {"all": cKDTree(xy_s)}
        for cname, labels in CLASSES.items():
            trees[cname] = cKDTree(xy_s[np.isin(ct, list(labels))])
        raw_m = raw[m][mal]
        log_m = logv[m][mal]
        fov_m = fov[m][mal]
        pos_frac = float((raw_m > 0).mean())
        if abs(pos_frac - PUBLISHED_POS[s]) > 1e-6:
            raise SystemExit(f"{s} TACSTD2>0 fraction {pos_frac} != published {PUBLISHED_POS[s]}")
        denom = {}
        counts = {cname: {} for cname in CLASSES}
        for radius in RADII:
            n_all = ball_counts(trees["all"], tum_xy, radius)
            if float(n_all.min()) < 1:
                raise SystemExit(f"{s} {radius} µm: a tumor cell did not find itself")
            denom[radius] = n_all - 1.0
            for cname in CLASSES:
                counts[cname][radius] = ball_counts(trees[cname], tum_xy, radius)
        sec = {
            "sample": s,
            "patient": PATIENT[s],
            "raw": raw_m,
            "log": log_m,
            "fov": fov_m,
            "denom": denom,
            "counts": counts,
        }
        sections.append(sec)
        inventory.append(
            {
                "sample": s,
                "patient": PATIENT[s],
                "tumor_label": TUMOR[s],
                "n_cells": int(m.sum()),
                "n_tumor": int(mal.sum()),
                "n_immune": int(np.isin(ct, list(BROAD)).sum()),
                "n_cd8": int(np.isin(ct, list(CD8)).sum()),
                "n_cd8nk": int(np.isin(ct, list(CD8NK)).sum()),
                "tacstd2_pos_frac": pos_frac,
                "tacstd2_median_log": float(np.median(log_m)),
                "median_split_is_detected": bool(np.median(log_m) == 0.0),
                "median_nn_um": med_nn,
                "n_zero_library": int((n_counts[m] <= 0).sum()),
            }
        )
        print(
            f"{s}: tumor={int(mal.sum())} pos={pos_frac:.3f} NN={med_nn:.2f} µm in {time.time() - t1:.1f}s",
            flush=True,
        )
    if sum(r["n_tumor"] for r in inventory) != 295877:
        raise SystemExit("tumor-cell total is not 295,877")
    return sections, inventory


def outcome_vector(sec, class_name: str, radius: int, estimand: str) -> np.ndarray:
    counts = sec["counts"][class_name][radius]
    denom = sec["denom"][radius]
    if estimand == "count":
        return counts
    out = np.full(counts.shape, np.nan, dtype=np.float64)
    ok = denom > 0
    out[ok] = counts[ok] / denom[ok]
    if estimand == "frac0":
        out[~ok] = 0.0
    elif estimand != "frac":
        raise KeyError(estimand)
    return out


def arm_stats(y: np.ndarray, high: np.ndarray, low: np.ndarray):
    yh = y[high]
    yl = y[low]
    yh = yh[np.isfinite(yh)]
    yl = yl[np.isfinite(yl)]
    if yh.size < MIN_N or yl.size < MIN_N:
        return None
    return {
        "n_high": int(yh.size),
        "n_low": int(yl.size),
        "high_mean": float(yh.mean()),
        "low_mean": float(yl.mean()),
        "high_sum": float(yh.sum()),
        "low_sum": float(yl.sum()),
    }


def wilcoxon_p(high, low) -> float:
    d = np.asarray(high, float) - np.asarray(low, float)
    d = d[np.isfinite(d) & (d != 0)]
    if d.size == 0:
        return float("nan")
    method = "exact" if d.size <= 25 else "auto"
    return float(stats.wilcoxon(d, alternative="two-sided", method=method, zero_method="wilcox").pvalue)


def sign_p(n_neg: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return float(stats.binomtest(int(n_neg), int(n), 0.5, alternative="greater").pvalue)


def score(sections):
    """Score every cut x radius x class x estimand. Return a summary frame."""
    # key -> list of per-section dicts, in sample order, or None if any section fails n
    bucket = {}
    t0 = time.time()
    for sec in sections:
        raw, logv, fov = sec["raw"], sec["log"], sec["fov"]
        sample = sec["sample"]
        vectors = {}
        for radius in RADII:
            for class_name in CLASSES:
                counts = sec["counts"][class_name][radius]
                denom = sec["denom"][radius]
                vectors[(radius, class_name, "count")] = counts
                frac = np.full(counts.shape, np.nan, dtype=np.float64)
                ok = denom > 0
                frac[ok] = counts[ok] / denom[ok]
                vectors[(radius, class_name, "frac")] = frac
                filled = frac.copy()
                filled[~ok] = 0.0
                vectors[(radius, class_name, "frac0")] = filled
        for cut in CUTS:
            high, low, thr = make_masks(cut, raw, logv, fov, sample)
            if int(high.sum()) < MIN_N or int(low.sum()) < MIN_N:
                for radius in RADII:
                    for class_name in CLASSES:
                        for estimand in ESTIMANDS:
                            bucket[(cut, radius, class_name, estimand)] = None
                continue
            for radius in RADII:
                for class_name in CLASSES:
                    for estimand in ESTIMANDS:
                        key = (cut, radius, class_name, estimand)
                        if bucket.get(key, "missing") is None:
                            continue
                        y = vectors[(radius, class_name, estimand)]
                        st = arm_stats(y, high, low)
                        if st is None:
                            bucket[key] = None
                            continue
                        st.update(
                            {
                                "sample": sample,
                                "patient": sec["patient"],
                                "threshold": thr,
                                "cut": cut,
                                "radius_um": int(radius),
                                "class_name": class_name,
                                "estimand": estimand,
                            }
                        )
                        bucket.setdefault(key, []).append(st)
        print(f"scored cuts for {sample} at {time.time() - t0:.1f}s", flush=True)

    rows = []
    detail_keep = []
    for key, parts in bucket.items():
        if not parts or len(parts) != 8:
            continue
        parts = sorted(parts, key=lambda r: SAMPLES.index(r["sample"]))
        hi = np.array([p["high_mean"] for p in parts], dtype=np.float64)
        lo = np.array([p["low_mean"] for p in parts], dtype=np.float64)
        if not np.isfinite(hi).all() or not np.isfinite(lo).all():
            continue
        sec_neg = int(np.sum(hi < lo))
        sec_pos = int(np.sum(hi > lo))
        tmp = pd.DataFrame({"patient": [p["patient"] for p in parts], "hi": hi, "lo": lo})
        ph, pl = [], []
        for _, g in tmp.groupby("patient", sort=False):
            ph.append(float(g["hi"].mean()))
            pl.append(float(g["lo"].mean()))
        ph = np.asarray(ph)
        pl = np.asarray(pl)
        pat_neg = int(np.sum(ph < pl))
        mean_hi = float(hi.mean())
        mean_lo = float(lo.mean())
        ratio = float(mean_hi / mean_lo) if mean_lo > 0 else float("nan")
        delta = float(mean_hi - mean_lo)
        with np.errstate(divide="ignore", invalid="ignore"):
            sec_ratios = hi / lo
        finite_ratios = sec_ratios[np.isfinite(sec_ratios)]
        n_high = np.array([p["n_high"] for p in parts])
        n_low = np.array([p["n_low"] for p in parts])
        high_events = np.array([p["high_sum"] for p in parts])
        low_events = np.array([p["low_sum"] for p in parts])
        pooled_hi = float(np.sum([p["high_sum"] for p in parts]) / np.sum(n_high))
        pooled_lo = float(np.sum([p["low_sum"] for p in parts]) / np.sum(n_low))
        pooled_ratio = float(pooled_hi / pooled_lo) if pooled_lo > 0 else float("nan")
        cut, radius, class_name, estimand = key
        full_sign = sec_neg == 8 and pat_neg == 5 and len(ph) == 5
        min_high = float(hi.min())
        min_low = float(lo.min())
        row = {
            "cut": cut,
            "radius_um": int(radius),
            "class_name": class_name,
            "estimand": estimand,
            "mean_high": mean_hi,
            "mean_low": mean_lo,
            "delta": delta,
            "ratio": ratio,
            "pooled_ratio": pooled_ratio,
            "mean_section_ratio": float(finite_ratios.mean()) if finite_ratios.size else float("nan"),
            "max_section_ratio": float(finite_ratios.max()) if finite_ratios.size else float("nan"),
            "min_section_ratio": float(finite_ratios.min()) if finite_ratios.size else float("nan"),
            "sections_high_lt_low": f"{sec_neg}/8",
            "patients_high_lt_low": f"{pat_neg}/5",
            "full_sign": bool(full_sign),
            "min_high": min_high,
            "min_low": min_low,
            "min_high_events": float(high_events.min()),
            "min_low_events": float(low_events.min()),
            "n_high_min": int(n_high.min()),
            "n_low_min": int(n_low.min()),
            "section_sign_p": sign_p(sec_neg, 8),
            "patient_sign_p": sign_p(pat_neg, 5),
            "section_wilcoxon_p": wilcoxon_p(hi, lo),
            "patient_wilcoxon_p": wilcoxon_p(ph, pl),
        }
        rows.append(row)
        if full_sign or (cut in {"log_median", "detected", "log_q4q1"} and radius in {10, 20, 50}):
            for p, h, l in zip(parts, hi, lo):
                detail_keep.append(
                    {
                        **{k: p[k] for k in ("sample", "patient", "cut", "radius_um", "class_name", "estimand", "threshold", "n_high", "n_low")},
                        "high_mean": float(h),
                        "low_mean": float(l),
                        "delta": float(h - l),
                        "ratio": float(h / l) if l > 0 else float("nan"),
                    }
                )
    res = pd.DataFrame(rows)
    detail = pd.DataFrame(detail_keep)
    print(f"scored {len(res)} complete specs in {time.time() - t0:.1f}s", flush=True)
    return res, detail


def apply_eligibility(res: pd.DataFrame) -> pd.DataFrame:
    res = res.copy()
    floor_ok = []
    for rec in res.itertuples(index=False):
        if rec.estimand in {"frac", "frac0"}:
            ok = rec.min_low >= FRAC_FLOOR and rec.min_high > 0
        elif rec.estimand == "count":
            ok = (
                rec.min_low >= COUNT_FLOOR
                and rec.min_high > 0
                and rec.min_low_events >= COUNT_MIN_LOW_EVENTS
                and rec.min_high_events >= COUNT_MIN_HIGH_EVENTS
            )
        else:
            ok = False
        floor_ok.append(bool(ok))
    res["floor_ok"] = floor_ok
    # frac0 also needs the empty-dropped fraction, same cut/radius/class, to be fully concordant
    frac = res[res["estimand"] == "frac"][
        ["cut", "radius_um", "class_name", "full_sign", "floor_ok", "min_high", "min_low"]
    ].rename(
        columns={
            "full_sign": "contact_full_sign",
            "floor_ok": "contact_floor_ok",
            "min_high": "contact_min_high",
            "min_low": "contact_min_low",
        }
    )
    res = res.merge(frac, on=["cut", "radius_um", "class_name"], how="left")
    eligible = res["full_sign"] & res["floor_ok"] & np.isfinite(res["ratio"]) & (res["ratio"] > 0)
    is_frac0 = res["estimand"] == "frac0"
    eligible = eligible & np.where(
        is_frac0,
        res["contact_full_sign"].fillna(False) & res["contact_floor_ok"].fillna(False),
        True,
    )
    res["eligible"] = eligible
    return res


def check_calibration(res: pd.DataFrame) -> None:
    for cut, class_name, estimand, radius, high, low, ratio in CALIBRATION:
        hit = res[
            (res["cut"] == cut)
            & (res["class_name"] == class_name)
            & (res["estimand"] == estimand)
            & (res["radius_um"] == radius)
        ]
        if hit.empty:
            raise SystemExit(f"calibration spec missing: {cut} {class_name} {estimand} {radius}")
        row = hit.iloc[0]
        for got, exp, label in (
            (row["mean_high"], high, "high"),
            (row["mean_low"], low, "low"),
            (row["ratio"], ratio, "ratio"),
        ):
            if abs(float(got) - float(exp)) > 5e-4:
                raise SystemExit(
                    f"calibration {cut} {class_name} {estimand} r={radius} {label}: got {got} expected {exp}"
                )
        if int(row["sections_high_lt_low"].split("/")[0]) != 8 or int(row["patients_high_lt_low"].split("/")[0]) != 5:
            raise SystemExit(f"calibration sign mismatch for {cut} {class_name} {estimand} r={radius}")
    print("calibration matched published means within 5e-4", flush=True)


def pick_ratio(df: pd.DataFrame) -> pd.Series | None:
    if df.empty:
        return None
    ordered = df.sort_values(["ratio", "delta", "radius_um", "cut", "estimand"], kind="mergesort")
    return ordered.iloc[0]


def pick_delta(df: pd.DataFrame) -> pd.Series | None:
    if df.empty:
        return None
    ordered = df.sort_values(["delta", "ratio", "radius_um", "cut", "estimand"], kind="mergesort")
    return ordered.iloc[0]


def section_detail_for(sections, spec: pd.Series) -> pd.DataFrame:
    rows = []
    for sec in sections:
        high, low, thr = make_masks(spec["cut"], sec["raw"], sec["log"], sec["fov"], sec["sample"])
        y = outcome_vector(sec, spec["class_name"], int(spec["radius_um"]), spec["estimand"])
        st = arm_stats(y, high, low)
        if st is None:
            raise SystemExit(f"winner lost its section stats: {sec['sample']} {spec['cut']}")
        # degree, to show whether the fraction drop is only fewer neighbors
        deg = sec["denom"][int(spec["radius_um"])]
        deg_h = float(deg[high].mean()) if high.any() else float("nan")
        deg_l = float(deg[low].mean()) if low.any() else float("nan")
        rows.append(
            {
                "sample": sec["sample"],
                "patient": sec["patient"],
                "threshold": thr,
                "n_high": st["n_high"],
                "n_low": st["n_low"],
                "high_mean": st["high_mean"],
                "low_mean": st["low_mean"],
                "delta": st["high_mean"] - st["low_mean"],
                "ratio": st["high_mean"] / st["low_mean"] if st["low_mean"] > 0 else float("nan"),
                "high_events": st["high_sum"],
                "low_events": st["low_sum"],
                "degree_high": deg_h,
                "degree_low": deg_l,
            }
        )
    return pd.DataFrame(rows)


def fov_test(sections, spec: pd.Series) -> dict:
    frames = []
    for sec in sections:
        high, low, _thr = make_masks(spec["cut"], sec["raw"], sec["log"], sec["fov"], sec["sample"])
        y = outcome_vector(sec, spec["class_name"], int(spec["radius_um"]), spec["estimand"])
        arm = np.full(len(y), "mid", dtype=object)
        arm[high] = "high"
        arm[low] = "low"
        use = (arm != "mid") & np.isfinite(y)
        frames.append(
            pd.DataFrame(
                {
                    "sample": np.full(int(use.sum()), sec["sample"], dtype=object),
                    "fov": sec["fov"][use],
                    "arm": arm[use],
                    "y": y[use],
                }
            )
        )
    cells = pd.concat(frames, ignore_index=True)
    gmu = cells.groupby(["sample", "fov", "arm"])["y"].mean()
    gct = cells.groupby(["sample", "fov", "arm"]).size()
    wide = gct.unstack("arm")
    if "high" not in wide.columns or "low" not in wide.columns:
        return {"n_fovs": 0, "fovs_high_lt_low": "0/0", "fov_p": float("nan"), "fov_median_delta": float("nan")}
    ok = wide.index[(wide["high"] >= 8) & (wide["low"] >= 8)]
    him = gmu.xs("high", level="arm")
    lom = gmu.xs("low", level="arm")
    both = him.index.intersection(lom.index).intersection(ok)
    hv = him.loc[both].to_numpy()
    lv = lom.loc[both].to_numpy()
    n_neg = int(np.sum(hv < lv))
    return {
        "n_fovs": int(len(both)),
        "fovs_high_lt_low": f"{n_neg}/{len(both)}",
        "fov_median_delta": float(np.median(hv - lv)) if len(both) else float("nan"),
        "fov_p": wilcoxon_p(hv, lv),
    }


def _fmt(x, digits=3):
    return f"{float(x):.{digits}f}"


def _fmt_p(p):
    p = float(p)
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def plot_paired(detail: pd.DataFrame, spec: pd.Series, path_stem: str, ylabel: str, title: str):
    fig, ax = plt.subplots(figsize=(5.4, 4.5))
    for _, row in detail.iterrows():
        color = PATIENT_COLOR[row["patient"]]
        ax.plot([0, 1], [row["low_mean"], row["high_mean"]], color=color, lw=1.6)
        ax.scatter([0, 1], [row["low_mean"], row["high_mean"]], color=color, s=32, zorder=3)
    ax.set_xticks([0, 1], ["TACSTD2 low", "TACSTD2 high"])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.25, 1.25)
    handles = [
        plt.Line2D([0], [0], color=PATIENT_COLOR[p], lw=2, label=p)
        for p in ("Lung5", "Lung6", "Lung9", "Lung12", "Lung13")
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path_stem + ".png", dpi=160)
    fig.savefig(path_stem + ".pdf")
    plt.close(fig)


def plot_ratio_curves(res: pd.DataFrame, winners: dict):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    show_cuts = ["detected", "log_median", "log_pos_q90", "ge5"]
    for ax, class_name, title in zip(axes, ("immune", "cd8nk"), ("Immune fraction", "CD8+NK count")):
        estimand = "frac" if class_name == "immune" else "count"
        sub = res[(res["class_name"] == class_name) & (res["estimand"] == estimand)]
        for cut in show_cuts:
            part = sub[sub["cut"] == cut].sort_values("radius_um")
            if part.empty:
                continue
            ax.plot(part["radius_um"], part["ratio"], lw=1.2, label=cut)
            elig = part["eligible"].to_numpy()
            ax.scatter(part.loc[elig, "radius_um"], part.loc[elig, "ratio"], s=14, zorder=3)
        win = winners.get((class_name, "ratio"))
        if win is not None:
            ax.scatter([win["radius_um"]], [win["ratio"]], marker="*", s=90, color="black", zorder=4, label="ratio call")
        ax.set_xlabel("Radius (µm)")
        ax.set_title(title)
        ax.set_ylim(0, 1.05)
        ax.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("High / low ratio")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ratio_vs_radius.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "ratio_vs_radius.pdf"))
    plt.close(fig)


def spec_phrase(spec: pd.Series) -> str:
    return (
        f"`{spec['cut']}`, {spec['class_name']}, {spec['estimand']}, {int(spec['radius_um'])} µm"
    )


def write_results(res, eligible, inventory, calls, details, fovs):
    n_scored = int(len(res))
    n_eligible = int(len(eligible))
    n_primary = int(eligible["class_name"].isin(PRIMARY_CLASSES).sum())
    lines = []
    a = lines.append
    a("# CosMx NSCLC: smallest TACSTD2-high immune and CD8+NK ratio at 8/8 and 5/5")
    a("")
    a(
        "He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; "
        "295,877 patient-matched tumor cells, 8 sections, 5 patients). "
        "This note searches the TACSTD2 high-versus-low neighbor contrast. "
        "It does **not** replace the locked CLDN4 cytotoxic ratios 0.36 at 50 µm and 0.52 at 100 µm, "
        "and it does **not** say that effector cells which are present have lower GZMB, PRF1, NKG7, or IFNG. "
        "No private 8-KL."
    )
    a("")
    a("## What was held fixed")
    a("")
    a(
        "Index cells are the author tumor label matched to the section (`tumor 5/6/9/12/13`). "
        f"Centroids are global pixels × {UM_PER_PX} µm. The index cell is excluded. "
        "Immune is the same 14 labels as the previous TACSTD2 note "
        "(B, plasmablast, CD4 naive/memory, CD8 naive/memory, Treg, NK, mDC, pDC, monocyte, macrophage, mast, neutrophil). "
        "CD8+NK is T CD8 memory, T CD8 naive, and NK. CD8 alone is scored and is not allowed to replace CD8+NK. "
        "The ratio is the unweighted mean of the eight section high means divided by the unweighted mean of the eight section low means. "
        "A patient value is the unweighted mean of that patient's sections."
    )
    a("")
    a(
        f"Radii were {', '.join(str(r) for r in RADII)} µm. "
        f"Cuts were {', '.join(CUTS)}. "
        "`detected` and `geK` compare raw TACSTD2 count > 0 or ≥ K with count = 0. "
        "`pos_q*` is the numpy quantile of positive raw counts versus count = 0. "
        "`log_*` uses log1p(count / n_counts × 10,000) inside the section. "
        "`rank10/20/30` takes that exact fraction of cells from each tail; ties are broken by a fixed per-section jitter (seed from the section name and the cut), not by file order. "
        "`fov_median` is the within-FOV log median. "
        "Estimands: `frac` (neighbor fraction, cells with no neighbor dropped), `count` (neighbor count, empty balls kept as zero), "
        f"`frac0` (empty balls kept as fraction 0). A fraction spec needs every low-arm section mean ≥ {FRAC_FLOOR} and every high-arm section mean > 0. "
        f"A count spec needs every low-arm section mean ≥ {COUNT_FLOOR}, every low arm summing to ≥ {COUNT_MIN_LOW_EVENTS:.0f} events, and every high arm summing to ≥ {COUNT_MIN_HIGH_EVENTS:.0f} events. "
        "`frac0` is eligible only when `frac` for the same masks is also 8/8 and 5/5 and clears the fraction floor. "
        f"{n_scored} specs had both arms at n ≥ {MIN_N} in every section. {n_eligible} of those were eligible; {n_primary} of the eligible specs are immune or CD8+NK."
    )
    a("")
    a("## Calibration")
    a("")
    a("These equal-weight means match the previous TACSTD2 note within 5×10⁻⁴. Ranking started only after that check.")
    a("")
    a("| Cut | Class | Estimand | µm | High | Low | Ratio |")
    a("|---|---|---|---:|---:|---:|---:|")
    for cut, class_name, estimand, radius, high, low, ratio in CALIBRATION:
        a(f"| {cut} | {class_name} | {estimand} | {radius} | {high:.4f} | {low:.4f} | {ratio:.3f} |")
    a("")
    a("Where the log median is 0 (Lung6, Lung9, Lung12, Lung13), `log_median` is the detected-versus-absent split.")
    a("")
    a("## Calls")
    a("")
    a("Section and patient sign tests sit on the 8/8 and 5/5 floors (sign p = 0.0039 and 0.031; Wilcoxon p = 0.0078 and 0.0625) whenever every unit has the same sign, so those p-values do not choose the spec. FOV p-values are nominal: FOVs sit inside five patients, and the grid was searched.")
    a("")
    a("| Call | Spec | High | Low | Ratio | Δ | Weakest section ratio | Pooled ratio | FOVs lower |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---|")
    for key, label in (
        (("immune", "ratio"), "Smallest immune ratio"),
        (("immune", "delta"), "Largest immune \\|Δ\\|"),
        (("cd8nk", "ratio"), "Smallest CD8+NK ratio"),
        (("cd8nk", "delta"), "Largest CD8+NK \\|Δ\\|"),
    ):
        spec = calls[key]
        fov = fovs[key]
        a(
            f"| {label} | {spec_phrase(spec)} | {spec['mean_high']:.4f} | {spec['mean_low']:.4f} | "
            f"**{spec['ratio']:.3f}** | **{spec['delta']:.4f}** | {spec['max_section_ratio']:.4f} | "
            f"{spec['pooled_ratio']:.3f} | {fov['fovs_high_lt_low']} |"
        )
    a("")
    a("Per estimand, still immune or CD8+NK, still eligible:")
    a("")
    a("| Class | Estimand | Smallest ratio | That Δ | Largest \\|Δ\\| | That ratio |")
    a("|---|---|---|---:|---|---:|")
    for class_name in PRIMARY_CLASSES:
        for estimand in ESTIMANDS:
            sub = eligible[(eligible["class_name"] == class_name) & (eligible["estimand"] == estimand)]
            rw = pick_ratio(sub)
            dw = pick_delta(sub)
            if rw is None:
                a(f"| {class_name} | {estimand} | none eligible |  | none eligible |  |")
                continue
            a(
                f"| {class_name} | {estimand} | {spec_phrase(rw)} = {rw['ratio']:.3f} | {rw['delta']:.4f} | "
                f"{spec_phrase(dw)} = {dw['delta']:.4f} | {dw['ratio']:.3f} |"
            )
    a("")
    cd8_elig = eligible[eligible["class_name"] == "cd8"]
    cd8_r = pick_ratio(cd8_elig)
    cd8_d = pick_delta(cd8_elig)
    if cd8_r is None:
        a("CD8 alone had no eligible spec. It does not replace the CD8+NK call.")
    else:
        a(
            f"CD8 alone, smallest eligible ratio {cd8_r['ratio']:.3f} ({spec_phrase(cd8_r)}), "
            f"largest |Δ| {cd8_d['delta']:.4f} ({spec_phrase(cd8_d)}). "
            "Reported so a CD8-only cut cannot be mistaken for the CD8+NK call."
        )
    a("")
    a("## How the calls sit next to the published contrasts")
    a("")
    imm_ratio = calls[("immune", "ratio")]
    if imm_ratio["estimand"] == "frac0":
        contact = res[
            (res["cut"] == imm_ratio["cut"])
            & (res["class_name"] == imm_ratio["class_name"])
            & (res["radius_um"] == int(imm_ratio["radius_um"]))
            & (res["estimand"] == "frac")
        ]
        counted = res[
            (res["cut"] == imm_ratio["cut"])
            & (res["class_name"] == imm_ratio["class_name"])
            & (res["radius_um"] == int(imm_ratio["radius_um"]))
            & (res["estimand"] == "count")
        ]
        if not contact.empty:
            c = contact.iloc[0]
            a(
                f"The immune ratio call is the zero-filled fraction: empty balls count as 0. "
                f"On the same masks, the fraction among cells that have a neighbor is **{c['ratio']:.3f}** "
                f"(high {c['mean_high']:.4f}, low {c['mean_low']:.4f}, Δ {c['delta']:.4f}, "
                f"{'eligible' if bool(c['eligible']) else 'not eligible'}). "
                f"Quietest contact-fraction high arm: n = {int(c['n_high_min'])}, "
                f"section mean {c['min_high']:.4f}."
            )
            a("")
        if not counted.empty:
            c = counted.iloc[0]
            a(
                f"The neighbor-count ratio on those same masks is {c['ratio']:.3f} "
                f"(Δ {c['delta']:.4f}; quietest high arm {c['min_high_events']:.1f} immune neighbors). "
                + (
                    "That count spec clears the event floor."
                    if bool(c["eligible"])
                    else "That count spec does not clear the event floor, so it is not a call."
                )
            )
        a("")
    pub = res[
        (res["cut"] == "log_median")
        & (res["class_name"] == "cd8nk")
        & (res["estimand"] == "count")
        & (res["radius_um"] == 10)
    ].iloc[0]
    a(
        f"The published CD8+NK count ratio at 10 µm, log median, is {pub['ratio']:.3f} "
        f"(Δ {pub['delta']:.4f}), 8/8 and 5/5. It is not eligible: the quietest low arm has "
        f"{pub['min_low_events']:.0f} events and the quietest high arm has {pub['min_high_events']:.0f}, "
        f"against floors of {COUNT_MIN_LOW_EVENTS:.0f} and {COUNT_MIN_HIGH_EVENTS:.0f}. "
        f"The eligible CD8+NK ratio is larger than that published ratio. "
        f"No cut in this grid is a smaller CD8+NK ratio that also clears the event floor."
    )
    a("")
    for key, label in (
        (("immune", "delta"), "Immune"),
        (("cd8nk", "delta"), "CD8+NK"),
    ):
        detail = details[key]
        weakest = detail.loc[detail["ratio"].idxmax()]
        a(
            f"{label} absolute-drop call, weakest section: {weakest['sample']} ratio {weakest['ratio']:.4f} "
            f"(Δ {weakest['delta']:.4f}, high {weakest['high_mean']:.4f}, low {weakest['low_mean']:.4f}). "
            f"The equal-weight Δ is not the drop in every section."
        )
    a("")

    def block(title, key, ylabel):
        spec = calls[key]
        detail = details[key]
        fov = fovs[key]
        a(f"## {title}")
        a("")
        a(
            f"{spec_phrase(spec)}. Equal-weight high {spec['mean_high']:.4f}, low {spec['mean_low']:.4f}, "
            f"ratio {spec['ratio']:.3f}, Δ {spec['delta']:.4f}. "
            f"Sections {spec['sections_high_lt_low']}, patients {spec['patients_high_lt_low']}. "
            f"Weakest section ratio {spec['max_section_ratio']:.4f}. "
            f"Cell-weighted ratio {spec['pooled_ratio']:.3f}. "
            f"Smallest high-arm section mean {spec['min_high']:.4f}; smallest low-arm section mean {spec['min_low']:.4f}. "
            f"FOVs with at least 8 cells in each arm: {fov['fovs_high_lt_low']} lower, "
            f"median Δ {fov['fov_median_delta']:.4f}, nominal Wilcoxon p = {_fmt_p(fov['fov_p'])}. "
            f"Mean degree {float(detail['degree_high'].mean()):.2f} (high) and {float(detail['degree_low'].mean()):.2f} (low)."
        )
        a("")
        a("| Section | Patient | n high | n low | High | Low | Ratio | Δ | High events | Low events | Degree high | Degree low |")
        a("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in detail.iterrows():
            a(
                f"| {r['sample']} | {r['patient']} | {int(r['n_high'])} | {int(r['n_low'])} | "
                f"{r['high_mean']:.4f} | {r['low_mean']:.4f} | {r['ratio']:.3f} | {r['delta']:.4f} | "
                f"{r['high_events']:.1f} | {r['low_events']:.1f} | {r['degree_high']:.2f} | {r['degree_low']:.2f} |"
            )
        a("")

    block("Smallest immune ratio", ("immune", "ratio"), "immune")
    block("Largest immune absolute drop", ("immune", "delta"), "immune")
    block("Smallest CD8+NK ratio", ("cd8nk", "ratio"), "cd8nk")
    block("Largest CD8+NK absolute drop", ("cd8nk", "delta"), "cd8nk")

    a("## Concordant specs that missed the floor")
    a("")
    below = res[res["full_sign"] & ~res["eligible"] & np.isfinite(res["ratio"]) & res["class_name"].isin(PRIMARY_CLASSES)]
    below = below.sort_values("ratio")
    if below.empty:
        a("No immune or CD8+NK spec was 8/8 and 5/5 and still ineligible.")
    else:
        b = below.iloc[0]
        a(
            f"The smallest immune or CD8+NK ratio with a strict 8/8 and 5/5 sign, ignoring the floor and the frac0 contact rule, "
            f"is **{b['ratio']:.3f}** ({spec_phrase(b)}; high {b['mean_high']:.4f}, low {b['mean_low']:.4f}, "
            f"min low {b['min_low']:.4f}, min high events {b['min_high_events']:.1f}, min low events {b['min_low_events']:.1f}). "
            "It is not a call. The floor is what keeps a near-empty arm from manufacturing a ratio."
        )
        a("")
        a("The five smallest such non-calls:")
        a("")
        a("| Cut | Class | Estimand | µm | Ratio | Δ | Min low | Min high events | Min low events |")
        a("|---|---|---|---:|---:|---:|---:|---:|---:|")
        for _, r in below.head(5).iterrows():
            a(
                f"| {r['cut']} | {r['class_name']} | {r['estimand']} | {int(r['radius_um'])} | {r['ratio']:.3f} | "
                f"{r['delta']:.4f} | {r['min_low']:.4f} | {r['min_high_events']:.1f} | {r['min_low_events']:.1f} |"
            )
    a("")
    a("## Ten smallest eligible immune or CD8+NK ratios")
    a("")
    primary = eligible[eligible["class_name"].isin(PRIMARY_CLASSES)].sort_values(["ratio", "delta", "radius_um"])
    a("| Cut | Class | Estimand | µm | Ratio | Δ | Max section ratio | Pooled ratio | Min n high |")
    a("|---|---|---|---:|---:|---:|---:|---:|---:|")
    for _, r in primary.head(10).iterrows():
        a(
            f"| {r['cut']} | {r['class_name']} | {r['estimand']} | {int(r['radius_um'])} | {r['ratio']:.3f} | "
            f"{r['delta']:.4f} | {r['max_section_ratio']:.3f} | {r['pooled_ratio']:.3f} | {int(r['n_high_min'])} |"
        )
    a("")
    a("## Ten largest eligible absolute drops, immune or CD8+NK")
    a("")
    by_delta = eligible[eligible["class_name"].isin(PRIMARY_CLASSES)].sort_values(["delta", "ratio"])
    a("| Cut | Class | Estimand | µm | Δ | Ratio | Max section ratio |")
    a("|---|---|---|---:|---:|---:|---:|")
    for _, r in by_delta.head(10).iterrows():
        a(
            f"| {r['cut']} | {r['class_name']} | {r['estimand']} | {int(r['radius_um'])} | {r['delta']:.4f} | "
            f"{r['ratio']:.3f} | {r['max_section_ratio']:.3f} |"
        )
    a("")
    a("Count Δ grows with the area of the ball. A larger count Δ at a long radius is a larger absolute drop, not automatically a stronger fold. The fold call is the ratio table.")
    a("")
    a("## What this does not claim")
    a("")
    a("- Not a re-estimate of the locked CLDN4 50/100 µm cytotoxic ratios 0.36 / 0.52.")
    a("- Not muzzling. Neighbor counts are not GZMB, PRF1, NKG7, or IFNG levels inside the effector cells that remain.")
    a("- Not a smaller section-level p-value than the 8/8 Wilcoxon floor.")
    a("- FOV p-values are nominal.")
    a("- CD8 alone is not the CD8+NK result.")
    a("- No ICI labels. No private 8-KL. No Visium same-spot correlation written as exclusion.")
    a("")
    a("```bash")
    a("python3 scripts/download_cosmx_nsclc_h5ad.py")
    a("python3 scripts/cosmx_tacstd2_max_effect.py")
    a("```")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    return text


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    sections, inventory = load_sections()
    pd.DataFrame(inventory).to_csv(os.path.join(TAB, "inventory.csv"), index=False)
    res, rough_detail = score(sections)
    check_calibration(res)
    res = apply_eligibility(res)
    res.to_csv(os.path.join(TAB, "parameter_grid.csv"), index=False)
    eligible = res[res["eligible"]].copy()
    eligible.to_csv(os.path.join(TAB, "eligible.csv"), index=False)
    rough_detail.to_csv(os.path.join(TAB, "sign_concordant_sections.csv"), index=False)

    calls = {}
    details = {}
    fovs = {}
    for class_name in PRIMARY_CLASSES:
        pool = eligible[eligible["class_name"] == class_name]
        if pool.empty:
            raise SystemExit(f"no eligible spec for {class_name}")
        calls[(class_name, "ratio")] = pick_ratio(pool)
        calls[(class_name, "delta")] = pick_delta(pool)
    for key, spec in calls.items():
        print(
            "CALL",
            key,
            spec_phrase(spec),
            "ratio",
            round(float(spec["ratio"]), 4),
            "delta",
            round(float(spec["delta"]), 4),
            flush=True,
        )
        details[key] = section_detail_for(sections, spec)
        # The recomputed equal-weight ratio has to match the grid row.
        recomputed = float(details[key]["high_mean"].mean() / details[key]["low_mean"].mean())
        if abs(recomputed - float(spec["ratio"])) > 5e-4:
            raise SystemExit(f"recomputed ratio {recomputed} != grid {spec['ratio']} for {key}")
        fovs[key] = fov_test(sections, spec)
        details[key].to_csv(os.path.join(TAB, f"{key[0]}_{key[1]}_by_section.csv"), index=False)

    ylabels = {
        ("immune", "ratio"): "Section mean",
        ("immune", "delta"): "Section mean",
        ("cd8nk", "ratio"): "Section mean",
        ("cd8nk", "delta"): "Section mean",
    }
    for key, spec in calls.items():
        estimand = spec["estimand"]
        noun = {"frac": "fraction", "count": "count", "frac0": "zero-filled fraction"}[estimand]
        kind = "smallest ratio" if key[1] == "ratio" else "largest |Δ|"
        plot_paired(
            details[key],
            spec,
            os.path.join(FIG, f"{key[0]}_{key[1]}_paired"),
            f"Section mean {spec['class_name']} {noun}",
            f"{kind}: {spec['cut']}, {int(spec['radius_um'])} µm\nratio {spec['ratio']:.3f}, Δ {spec['delta']:.3f}",
        )
    plot_ratio_curves(res, calls)

    summary = {
        "n_specs_scored": int(len(res)),
        "n_eligible": int(eligible.shape[0]),
        "n_eligible_primary": int(eligible["class_name"].isin(PRIMARY_CLASSES).sum()),
        "floors": {
            "fraction_min_low": FRAC_FLOOR,
            "count_min_low": COUNT_FLOOR,
            "count_min_low_events": COUNT_MIN_LOW_EVENTS,
            "count_min_high_events": COUNT_MIN_HIGH_EVENTS,
            "min_n": MIN_N,
        },
        "calls": {f"{a}_{b}": calls[(a, b)].to_dict() for a, b in calls},
        "fovs": {f"{a}_{b}": fovs[(a, b)] for a, b in fovs},
        "inventory": inventory,
        "radii_um": list(RADII),
        "cuts": list(CUTS),
        "classes": {k: list(v) for k, v in CLASSES.items()},
    }
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    write_results(res, eligible, inventory, calls, details, fovs)
    print(f"wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
