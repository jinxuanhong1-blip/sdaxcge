#!/usr/bin/env python3
"""Donor permutation and patient bootstrap for two locked contrasts.

CosMx (He et al. 2022, figshare 25976224): CLDN4-high versus CLDN4-low tumor
cells, cytotoxic and T/NK neighbor summaries. Five patients are the unit.
Sections inside a patient are averaged with equal section weight before the
patient median is taken.

Concordant-4 (GSE123902, GSE131907, GSE205335, GSE189357): malignant CLDN4
versus the patient T/NK fraction. Sixty-five patients are the unit. Quartile
labels are formed inside each cohort. Pairs do not cross cohorts.

The searched summary is the patient-level median contrast. A spec can win
only when its patient-bootstrap 95% interval lies entirely below zero
(CLDN4-high has fewer immune neighbors, or CLDN4-high patients have a lower
T/NK fraction). Among those specs, the largest |median| wins. A narrower
interval breaks ties. The permutation p-value is the tail of the best
exclusion median in the same grid, so the search is part of the null.

The locked CosMx ratios 0.36 / 0.52 and the locked Spearman −0.53 are not
replaced. Cell-weighted ratios are descriptive only.
"""

from __future__ import annotations

import csv
import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from patient_separation_stats import (
    ESTIMANDS,
    bootstrap_median_ci,
    cohort_spearman_dl,
    corrected_log2_ratio,
    exclusion_magnitude,
    labels_within_groups,
    median_of_cohort_diffs,
    pairwise_median_diff,
    percentile_interval,
    perm_p_ge,
    permute_within_groups,
    pick_winner,
    pooled_median_diff,
    resample_within_cohort,
    rollup_section_means,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
CACHE = os.path.join(ROOT, "data", "cosmx_nsclc", "tumor_neighbor_cache.npz")
PATIENT_TSV = os.path.join(ROOT, "inputs", "concordant4_patient_units.tsv")
OUT = os.path.join(ROOT, "results", "patient_separation")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

UM_PER_PX = 0.18
RADII_UM = np.array([10.0, 20.0, 30.0, 50.0, 80.0, 100.0])
CLASSES = ("cyto", "cd8", "nk", "tnk", "immune")
CUTS = ("detected", "median", "quartile")
METRICS = ("count", "contact", "fraction")
# Reference arm must actually contain the neighbor class. eps shrinks log2
# ratios so a zero high arm stays finite. Both were fixed before the grid ran.
MIN_LOW = {"count": 0.05, "contact": 0.02, "fraction": 0.01}
EPS = {"count": 0.05, "contact": 0.02, "fraction": 0.005}
MIN_ARM = 30
N_PERM = 999
N_BOOT = 4999
SEED = 25976224

CYTO = ("NK", "T CD8 memory", "T CD8 naive")
CD8 = ("T CD8 memory", "T CD8 naive")
NK = ("NK",)
TNK = ("NK", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive", "Treg")
IMMUNE = (
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
CLASS_MEMBERS = {"cyto": CYTO, "cd8": CD8, "nk": NK, "tnk": TNK, "immune": IMMUNE}


def _decode(arr) -> list[str]:
    out = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode())
        else:
            out.append(str(x))
    return out


def extract_tumor_neighbors(path: str, cache: str) -> dict:
    if os.path.isfile(cache):
        print(f"loading neighbor cache {cache}", flush=True)
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}

    import h5py
    from scipy.spatial import cKDTree

    t0 = time.time()
    with h5py.File(path, "r") as f:
        genes = _decode(f["var"]["_index"][:])
        gene = genes.index("CLDN4")
        n = int(f["obs"]["_index"].shape[0])
        if n != 765771:
            raise SystemExit(f"expected 765771 cells, found {n}")
        indptr = f["layers"]["counts"]["indptr"][:]
        indices = f["layers"]["counts"]["indices"]
        data = f["layers"]["counts"]["data"]
        cldn4_all = np.zeros(n, dtype=np.float32)
        step = 50000
        for start in range(0, n, step):
            end = min(n, start + step)
            a = int(indptr[start])
            b = int(indptr[end])
            if b <= a:
                continue
            cols = indices[a:b]
            hit = np.flatnonzero(cols == gene)
            if hit.size == 0:
                continue
            rows = np.searchsorted(indptr[start : end + 1], a + hit, side="right") - 1
            cldn4_all[start + rows] = data[a + hit]
        xy = np.asarray(f["obsm"]["spatial"][:], dtype=np.float64)
        cell_codes = f["obs"]["cell_type"]["codes"][:]
        sample = f["obs"]["sample"]["codes"][:].astype(np.int8)
        patient = f["obs"]["patient"]["codes"][:].astype(np.int8)
        fov = f["obs"]["fov"][:].astype(np.int16)
        cell_names = _decode(f["obs"]["cell_type"]["categories"][:])
        sample_names = np.array(_decode(f["obs"]["sample"]["categories"][:]))
        patient_names = np.array(_decode(f["obs"]["patient"]["categories"][:]))

    name_to_code = {nm: i for i, nm in enumerate(cell_names)}
    multi = np.zeros((n, len(CLASSES)), dtype=bool)
    for c, cname in enumerate(CLASSES):
        for nm in CLASS_MEMBERS[cname]:
            multi[cell_codes == name_to_code[nm], c] = True
    tumor_codes = [i for i, nm in enumerate(cell_names) if nm.startswith("tumor ")]
    if len(tumor_codes) != 5:
        raise SystemExit(f"expected 5 tumor labels, found {tumor_codes}")
    is_tumor = np.isin(cell_codes, np.array(tumor_codes, dtype=np.int8))
    tumor_rows = np.flatnonzero(is_tumor)
    if tumor_rows.size != 302313:
        raise SystemExit(f"expected 302313 tumor cells, found {tumor_rows.size}")
    slot = np.full(n, -1, dtype=np.int32)
    slot[tumor_rows] = np.arange(tumor_rows.size, dtype=np.int32)

    count = np.zeros((tumor_rows.size, len(RADII_UM), len(CLASSES)), dtype=np.int32)
    n_all = np.zeros((tumor_rows.size, len(RADII_UM)), dtype=np.int32)
    radii_px = RADII_UM / UM_PER_PX
    group_key = sample.astype(np.int32) * 1000 + fov.astype(np.int32)
    n_fov = 0
    for g in np.unique(group_key):
        local = np.flatnonzero(group_key == g)
        if local.size < 2 or not np.any(is_tumor[local]):
            continue
        n_fov += 1
        tree = cKDTree(xy[local])
        coo = tree.sparse_distance_matrix(tree, max_distance=float(radii_px[-1]), output_type="coo_matrix")
        ri = coo.row
        ci = coo.col
        dist = coo.data
        keep = (ri != ci) & is_tumor[local[ri]]
        if not np.any(keep):
            continue
        ri = ri[keep]
        ci = ci[keep]
        dist = dist[keep]
        slots = slot[local[ri]]
        dest_multi = multi[local[ci]]
        # np.add.at on a non-contiguous slice (count[:, r, c]) spills into
        # neighboring class columns. Scatter into a flat buffer, then copy.
        buf = np.zeros(tumor_rows.size, dtype=np.int32)
        for r_i, rpx in enumerate(radii_px):
            m = dist <= rpx + 1e-6
            if not np.any(m):
                continue
            sl = slots[m]
            buf[:] = 0
            np.add.at(buf, sl, 1)
            n_all[:, r_i] += buf
            mh = dest_multi[m]
            for c in range(len(CLASSES)):
                sel = mh[:, c]
                buf[:] = 0
                if np.any(sel):
                    np.add.at(buf, sl[sel], 1)
                count[:, r_i, c] += buf
        if n_fov % 40 == 0:
            print(f"  FOVs {n_fov}  {time.time() - t0:.1f}s", flush=True)

    payload = {
        "cldn4": cldn4_all[tumor_rows].astype(np.float32),
        "sample": sample[tumor_rows],
        "patient": patient[tumor_rows],
        "fov": fov[tumor_rows],
        "row_id": tumor_rows.astype(np.int64),
        "count": count,
        "n_all": n_all,
        "sample_names": sample_names,
        "patient_names": patient_names,
        "class_names": np.array(CLASSES),
        "radii_um": RADII_UM,
    }
    # CD8 and NK are disjoint author labels, so cytotoxic counts are their sum
    # and neither can exceed the number of neighbors.
    if np.any(count[:, :, 1] + count[:, :, 2] != count[:, :, 0]):
        raise SystemExit("cytotoxic counts are not CD8+NK; neighbor tally is wrong")
    if np.any(count[:, :, 0] > n_all):
        raise SystemExit("class counts exceed the neighbor degree")
    if float(n_all[:, -1].mean()) < 10:
        raise SystemExit(f"100 µm mean degree {float(n_all[:, -1].mean()):.2f} is too small")
    print(f"max neighbor count {int(count.max())} max degree {int(n_all.max())}", flush=True)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez_compressed(cache, **payload)
    print(f"cached neighbors in {time.time() - t0:.1f}s across {n_fov} FOVs", flush=True)
    return payload


def _sample_patient_map(sample: np.ndarray, patient: np.ndarray, n_sample: int) -> np.ndarray:
    out = np.full(n_sample, -1, dtype=int)
    for s in range(n_sample):
        m = sample == s
        if not np.any(m):
            continue
        vals = np.unique(patient[m])
        if vals.size != 1:
            raise SystemExit(f"sample {s} maps to patients {vals}")
        out[s] = int(vals[0])
    return out


def metric_matrix(cache: dict, metric: str) -> np.ndarray:
    count = cache["count"].astype(np.float64)
    if metric == "count":
        return count.reshape(count.shape[0], -1)
    if metric == "contact":
        return (count > 0).astype(np.float64).reshape(count.shape[0], -1)
    if metric == "fraction":
        n_all = cache["n_all"].astype(np.float64)
        frac = np.full_like(count, np.nan)
        ok = n_all > 0
        frac[ok] = count[ok] / n_all[ok, None]
        return frac.reshape(count.shape[0], -1)
    raise ValueError(metric)


def section_arm_means(y: np.ndarray, lab: np.ndarray, sample: np.ndarray, n_sample: int) -> tuple[np.ndarray, np.ndarray]:
    k = y.shape[1]
    high = np.full((n_sample, k), np.nan)
    low = np.full((n_sample, k), np.nan)
    finite = np.isfinite(y)
    y0 = np.where(finite, y, 0.0)
    for s in range(n_sample):
        idx = np.flatnonzero(sample == s)
        if idx.size == 0:
            continue
        lab_s = lab[idx]
        ys = y0[idx]
        fs = finite[idx].astype(np.float64)
        for arm, dest in ((1, high), (0, low)):
            w = (lab_s == arm).astype(np.float64)
            cnt = w @ fs
            sm = w @ ys
            ok = cnt >= MIN_ARM
            dest[s, ok] = sm[ok] / cnt[ok]
    return high, low


def permute_patient_means(
    y: np.ndarray,
    lab: np.ndarray,
    sample: np.ndarray,
    fov: np.ndarray,
    sample_patient: np.ndarray,
    n_patients: int,
    rng: np.random.Generator,
    n_perm: int,
    batch: int = 25,
) -> tuple[np.ndarray, np.ndarray]:
    """Within-(section, FOV) shuffles of the CLDN4 arm. Neighbor counts stay put."""
    n, k = y.shape
    y0 = np.where(np.isfinite(y), y, 0.0).astype(np.float32)
    finite = np.isfinite(y).astype(np.float32)
    out_h = np.full((n_perm, n_patients, k), np.nan)
    out_l = np.full((n_perm, n_patients, k), np.nan)
    group = sample.astype(np.int32) * 1000 + fov.astype(np.int32)
    groups = [np.flatnonzero(group == g) for g in np.unique(group)]
    sample_idx = [np.flatnonzero(sample == s) for s in range(int(sample.max()) + 1)]
    n_sample = len(sample_idx)
    done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        lab_b = np.empty((b, n), dtype=np.int8)
        for idx in groups:
            draw = rng.random((b, idx.size))
            order = np.argsort(draw, axis=1)
            lab_b[:, idx] = lab[idx][order]
        sec_h = np.full((b, n_sample, k), np.nan, dtype=np.float64)
        sec_l = np.full((b, n_sample, k), np.nan, dtype=np.float64)
        for s, idx in enumerate(sample_idx):
            if idx.size == 0:
                continue
            ys = y0[idx]
            fs = finite[idx]
            lb = lab_b[:, idx]
            for arm, dest in ((1, sec_h), (0, sec_l)):
                w = (lb == arm).astype(np.float32)
                cnt = w @ fs
                sm = w @ ys
                mean = np.full((b, k), np.nan, dtype=np.float64)
                ok = cnt >= MIN_ARM
                np.divide(sm, cnt, out=mean, where=ok)
                dest[:, s, :] = mean
        for p in range(n_patients):
            ss = np.flatnonzero(sample_patient == p)
            block_h = sec_h[:, ss, :]
            block_l = sec_l[:, ss, :]
            out_h[done : done + b, p, :] = _batch_nanmean(block_h)
            out_l[done : done + b, p, :] = _batch_nanmean(block_l)
        done += b
    return out_h, out_l


def _batch_nanmean(block: np.ndarray) -> np.ndarray:
    """Mean over axis 1, ignoring NaN. block is (batch, n_section, K)."""
    cnt = np.sum(np.isfinite(block), axis=1)
    total = np.nansum(block, axis=1)
    out = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(total, cnt, out=out, where=cnt > 0)
    return out


def _patient_vectors_ok(high: np.ndarray, low: np.ndarray, min_low: float) -> np.ndarray:
    """Boolean over the last-but-patient axis is handled by callers. high/low (..., n_patients)."""
    return np.isfinite(high) & np.isfinite(low) & (low >= min_low) & (high >= 0)


def cosmx_grid(cache: dict, rng: np.random.Generator) -> dict:
    sample = cache["sample"].astype(int)
    patient = cache["patient"].astype(int)
    fov = cache["fov"].astype(int)
    cldn4 = cache["cldn4"].astype(np.float64)
    row_id = cache["row_id"].astype(np.int64)
    patient_names = [str(x) for x in cache["patient_names"]]
    n_patients = len(patient_names)
    n_sample = len(cache["sample_names"])
    sample_patient = _sample_patient_map(sample, patient, n_sample)
    n_class = len(CLASSES)
    n_rad = len(RADII_UM)

    labels = {
        cut: labels_within_groups(cldn4, sample, cut, row_id)
        for cut in CUTS
    }
    print("CosMx arm sizes (high/low)", flush=True)
    for cut, lab in labels.items():
        print(f"  {cut}: high {(lab == 1).sum()} low {(lab == 0).sum()} mid {(lab < 0).sum()}", flush=True)

    rows: list[dict] = []
    null_max = np.zeros(N_PERM, dtype=np.float64)
    # Nominal null for every spec is filled after we know the winner column.
    spec_index: list[tuple[str, str, int]] = []

    for metric in METRICS:
        y = metric_matrix(cache, metric)
        eps = EPS[metric]
        min_low = MIN_LOW[metric]
        for cut in CUTS:
            lab = labels[cut]
            print(f"permutation {metric} {cut}", flush=True)
            t0 = time.time()
            perm_h, perm_l = permute_patient_means(
                y, lab, sample, fov, sample_patient, n_patients, rng, N_PERM
            )
            print(f"  perm means {time.time() - t0:.1f}s", flush=True)
            obs_h_sec, obs_l_sec = section_arm_means(y, lab, sample, n_sample)
            obs_h, obs_l = rollup_section_means(obs_h_sec, obs_l_sec, sample_patient, n_patients)
            for col in range(y.shape[1]):
                r_i, c_i = divmod(col, n_class)
                high = obs_h[:, col]
                low = obs_l[:, col]
                ok = _patient_vectors_ok(high, low, min_low)
                log2 = corrected_log2_ratio(high, low, eps)
                log2 = np.where(ok, log2, np.nan)
                eligible = bool(np.all(ok))
                if eligible:
                    med, lo, hi = bootstrap_median_ci(log2, rng, N_BOOT)
                else:
                    med, lo, hi = (float(np.nanmedian(log2)), math_nan(), math_nan())
                raw_delta = high - low
                raw_log2 = np.full(n_patients, np.nan)
                pos = ok & (low > 0) & (high > 0)
                raw_log2[pos] = np.log2(high[pos] / low[pos])
                row = {
                    "metric": metric,
                    "cut": cut,
                    "radius_um": float(RADII_UM[r_i]),
                    "neighbor": CLASSES[c_i],
                    "eligible": eligible,
                    "median": med,
                    "lo": lo,
                    "hi": hi,
                    "width": (hi - lo) if np.isfinite(hi) and np.isfinite(lo) else math_nan(),
                    "magnitude": exclusion_magnitude(med, hi) if eligible else 0.0,
                    "patient_log2": log2,
                    "patient_high": high,
                    "patient_low": low,
                    "patient_delta": raw_delta,
                    "patient_raw_log2": raw_log2,
                    "equal_patient_ratio": _equal_weight_ratio(high, low) if eligible else math_nan(),
                }
                rows.append(row)
                spec_index.append((metric, cut, col))
                perm_log2 = corrected_log2_ratio(perm_h[:, :, col], perm_l[:, :, col], eps)
                perm_ok = _patient_vectors_ok(perm_h[:, :, col], perm_l[:, :, col], min_low).all(axis=1)
                perm_log2 = np.where(perm_ok[:, None], perm_log2, np.nan)
                null_max = np.maximum(null_max, _perm_magnitudes(perm_log2))
            # Keep perm arrays only for the columns we still need? Nominal p is computed
            # from a second, smaller pass after the winner is known, to avoid storing
            # every permutation. Recompute the winner column below.

    winner = pick_winner(rows)
    observed_mag = 0.0 if winner is None else float(winner["magnitude"])
    p_select = perm_p_ge(null_max, observed_mag)

    nominal = None
    if winner is not None:
        y = metric_matrix(cache, winner["metric"])
        lab = labels[winner["cut"]]
        r_i = int(np.where(np.isclose(RADII_UM, winner["radius_um"]))[0][0])
        col = r_i * n_class + CLASSES.index(winner["neighbor"])
        print("nominal permutation for the winning spec", flush=True)
        perm_h, perm_l = permute_patient_means(
            y, lab, sample, fov, sample_patient, n_patients, rng, N_PERM
        )
        eps = EPS[winner["metric"]]
        min_low = MIN_LOW[winner["metric"]]
        perm_log2 = corrected_log2_ratio(perm_h[:, :, col], perm_l[:, :, col], eps)
        perm_ok = _patient_vectors_ok(perm_h[:, :, col], perm_l[:, :, col], min_low).all(axis=1)
        perm_log2 = np.where(perm_ok[:, None], perm_log2, np.nan)
        nominal_stats = _perm_magnitudes(perm_log2)
        nominal = {
            "p": perm_p_ge(nominal_stats, observed_mag),
            "null_mean": float(np.mean(nominal_stats)),
            "null_max": float(np.max(nominal_stats)),
        }
        winner["nominal_p"] = nominal["p"]
        winner["patient_names"] = patient_names
        sec_h, sec_l = section_arm_means(y, lab, sample, n_sample)
        sample_names = [str(x) for x in cache["sample_names"]]
        sections = []
        for s, name in enumerate(sample_names):
            sections.append(
                {
                    "sample": name,
                    "patient": patient_names[int(sample_patient[s])],
                    "mean_high": float(sec_h[s, col]),
                    "mean_low": float(sec_l[s, col]),
                }
            )
        winner["sections"] = sections

    # Descriptive cell-weighted ratios for the locked radius family.
    reference = reference_ratios(cache, labels, sample)

    return {
        "rows": rows,
        "winner": winner,
        "null_max": null_max,
        "p_select": p_select,
        "nominal": nominal,
        "reference": reference,
        "patient_names": patient_names,
        "n_tumor": int(cldn4.size),
        "n_cldn4_pos": int(np.sum(cldn4 > 0)),
        "sample_names": [str(x) for x in cache["sample_names"]],
        "sample_patient": sample_patient.tolist(),
    }


def math_nan() -> float:
    return float("nan")


def _equal_weight_ratio(high: np.ndarray, low: np.ndarray) -> float:
    if np.any(low <= 0) or np.any(~np.isfinite(high)) or np.any(~np.isfinite(low)):
        return math_nan()
    return float(np.mean(high) / np.mean(low))


def _perm_magnitudes(log2: np.ndarray) -> np.ndarray:
    """log2 is (n_perm, n_patients) with NaN when that draw is ineligible."""
    bad = ~np.isfinite(log2).all(axis=1)
    filled = np.where(np.isfinite(log2), log2, 0.0)
    med = np.median(filled, axis=1)
    hi = np.max(filled, axis=1)
    mag = np.where((~bad) & (med < 0) & (hi < 0), -med, 0.0)
    return mag.astype(np.float64)


def reference_ratios(cache: dict, labels: dict, sample: np.ndarray) -> list[dict]:
    """Cell-weighted and patient-equal ratios for the locked cytotoxic radii.

    Cell-weighted numbers are descriptive. They are not the inferential unit.
    """
    count = cache["count"].astype(np.float64)
    patient = cache["patient"].astype(int)
    n_patients = len(cache["patient_names"])
    n_sample = len(cache["sample_names"])
    sample_patient = _sample_patient_map(sample, patient, n_sample)
    out = []
    for radius in (50.0, 100.0):
        r_i = int(np.where(np.isclose(RADII_UM, radius))[0][0])
        y = count[:, r_i, CLASSES.index("cyto")]
        for cut, lab in labels.items():
            high = lab == 1
            low = lab == 0
            cell_ratio = float(y[high].mean() / y[low].mean()) if y[low].mean() > 0 else math_nan()
            sec_h = np.full(n_sample, np.nan)
            sec_l = np.full(n_sample, np.nan)
            sec_ratio = np.full(n_sample, np.nan)
            for s in range(n_sample):
                m = sample == s
                h = m & high
                l = m & low
                if h.sum() < MIN_ARM or l.sum() < MIN_ARM:
                    continue
                sec_h[s] = float(y[h].mean())
                sec_l[s] = float(y[l].mean())
                if sec_l[s] > 0:
                    sec_ratio[s] = sec_h[s] / sec_l[s]
            ph, pl = rollup_section_means(sec_h[:, None], sec_l[:, None], sample_patient, n_patients)
            ph = ph[:, 0]
            pl = pl[:, 0]
            log2 = corrected_log2_ratio(ph, pl, EPS["count"])
            med, lo, hi = bootstrap_median_ci(log2, np.random.default_rng(SEED), N_BOOT)
            out.append(
                {
                    "neighbor": "cyto",
                    "metric": "count",
                    "radius_um": radius,
                    "cut": cut,
                    "cell_weighted_ratio": cell_ratio,
                    "median_section_ratio": float(np.nanmedian(sec_ratio)),
                    "equal_patient_ratio": _equal_weight_ratio(ph, pl),
                    "median_patient_log2": med,
                    "lo": lo,
                    "hi": hi,
                    "patient_high": ph.tolist(),
                    "patient_low": pl.tolist(),
                    "patient_ratio": (ph / pl).tolist(),
                }
            )
    return out


def load_concordant(path: str) -> dict:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 65:
        raise SystemExit(f"expected 65 concordant-4 patients, found {len(rows)}")
    cohort_name = np.array([r["dataset"] for r in rows])
    codes, inverse = np.unique(cohort_name, return_inverse=True)
    unit = np.array([r["unit_id"] for r in rows])
    # Alphabetical unit id, so a tie does not depend on file order.
    order = np.argsort(unit, kind="mergesort")
    tie = np.empty(unit.size, dtype=np.int64)
    tie[order] = np.arange(unit.size)
    return {
        "cohort_name": cohort_name,
        "cohort": inverse.astype(int),
        "cohort_codes": np.array([str(x) for x in codes]),
        "unit": unit,
        "tie": tie,
        "pct": np.array([float(r["mal_CLDN4_pct"]) for r in rows]),
        "mean": np.array([float(r["mal_CLDN4_mean"]) for r in rows]),
        "frac": np.array([float(r["frac_tnk"]) for r in rows]),
        "tissue": np.array([r["tissue"] for r in rows]),
    }


def concordant_grid(data: dict, rng: np.random.Generator) -> dict:
    y = data["frac"]
    cohort = data["cohort"]
    tie = data["tie"]
    scores = {"mal_CLDN4_pct": data["pct"], "mal_CLDN4_mean": data["mean"]}
    modes = ("median", "tertile", "quartile")
    rows = []
    for score_name, score in scores.items():
        for mode in modes:
            for est_name, fn in ESTIMANDS.items():
                point, lo, hi, n_used, reps = _boot_concordant(y, score, cohort, tie, mode, fn, rng, N_BOOT)
                eligible = bool(np.isfinite(point) and n_used >= int(0.95 * N_BOOT) and np.isfinite(lo) and np.isfinite(hi))
                med = float(point) if np.isfinite(point) else math_nan()
                rows.append(
                    {
                        "score": score_name,
                        "cut": mode,
                        "estimand": est_name,
                        "eligible": eligible,
                        "median": med,
                        "lo": lo if eligible else math_nan(),
                        "hi": hi if eligible else math_nan(),
                        "width": (hi - lo) if eligible else math_nan(),
                        "magnitude": exclusion_magnitude(med, hi) if eligible else 0.0,
                        "n_boot_used": n_used,
                    }
                )
    winner = pick_winner(rows)
    null_max = np.zeros(N_PERM, dtype=np.float64)
    score_perm_ready = {name: score.copy() for name, score in scores.items()}
    for p in range(N_PERM):
        best = 0.0
        for score_name, score in score_perm_ready.items():
            shuf = permute_within_groups(score, cohort, rng)
            for mode in modes:
                lab = labels_within_groups(shuf, cohort, mode, tie)
                for fn in ESTIMANDS.values():
                    d = fn(y, lab, cohort)
                    if np.isfinite(d) and d < 0:
                        best = max(best, float(-d))
        null_max[p] = best
    observed = 0.0 if winner is None else float(winner["magnitude"])
    p_select = perm_p_ge(null_max, observed)
    nominal = None
    labels = None
    if winner is not None:
        score = scores[winner["score"]]
        fn = ESTIMANDS[winner["estimand"]]
        labels = labels_within_groups(score, cohort, winner["cut"], tie)
        null_nom = np.empty(N_PERM)
        for p in range(N_PERM):
            shuf = permute_within_groups(score, cohort, rng)
            lab = labels_within_groups(shuf, cohort, winner["cut"], tie)
            d = fn(y, lab, cohort)
            null_nom[p] = float(-d) if np.isfinite(d) and d < 0 else 0.0
        nominal = {"p": perm_p_ge(null_nom, observed), "null": null_nom}
        winner["nominal_p"] = nominal["p"]
        winner["labels"] = labels
        cohorts = []
        for i, name in enumerate(data["cohort_codes"]):
            m = cohort == i
            h = y[m & (labels == 1)]
            low = y[m & (labels == 0)]
            cohorts.append(
                {
                    "dataset": str(name),
                    "n_high": int(h.size),
                    "n_low": int(low.size),
                    "median_high": float(np.median(h)) if h.size else math_nan(),
                    "median_low": float(np.median(low)) if low.size else math_nan(),
                    "delta": float(np.median(h) - np.median(low)) if h.size and low.size else math_nan(),
                }
            )
        winner["cohorts"] = cohorts
    dl = cohort_spearman_dl(data["pct"], y, cohort)
    boot_rho = np.empty(N_BOOT)
    for b in range(N_BOOT):
        take = resample_within_cohort(cohort, rng)
        boot_rho[b] = cohort_spearman_dl(data["pct"][take], y[take], cohort[take])["rho"]
    rho_lo, rho_hi, rho_n = percentile_interval(boot_rho)
    null_rho = np.empty(N_PERM)
    for p in range(N_PERM):
        shuf = permute_within_groups(data["pct"], cohort, rng)
        null_rho[p] = cohort_spearman_dl(shuf, y, cohort)["rho"]
    # One-sided: more negative rho is more extreme.
    p_rho = perm_p_ge(-null_rho, -dl["rho"])
    return {
        "rows": rows,
        "winner": winner,
        "null_max": null_max,
        "p_select": p_select,
        "nominal": nominal,
        "dl": dl,
        "rho_lo": rho_lo,
        "rho_hi": rho_hi,
        "rho_n": rho_n,
        "p_rho": p_rho,
        "cohort_codes": data["cohort_codes"],
        "labels": labels,
        "data": data,
    }


def _boot_concordant(y, score, cohort, tie, mode, fn, rng, n_boot):
    lab = labels_within_groups(score, cohort, mode, tie)
    point = fn(y, lab, cohort)
    reps = np.empty(n_boot)
    for b in range(n_boot):
        take = resample_within_cohort(cohort, rng)
        lab_b = labels_within_groups(score[take], cohort[take], mode, tie[take])
        reps[b] = fn(y[take], lab_b, cohort[take])
    lo, hi, n_used = percentile_interval(reps)
    return float(point), lo, hi, n_used, reps


def _fmt(x, digits=3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and (ax < 0.001 or ax >= 1000):
        return f"{float(x):.3e}"
    return f"{float(x):.{digits}f}"


def _write_csv(path: str, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def save_outputs(cosmx: dict, conc: dict) -> None:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    grid_fields = [
        "metric",
        "cut",
        "radius_um",
        "neighbor",
        "eligible",
        "median",
        "lo",
        "hi",
        "width",
        "magnitude",
        "equal_patient_ratio",
    ]
    _write_csv(os.path.join(TAB, "cosmx_grid.csv"), cosmx["rows"], grid_fields)
    ref_rows = []
    for row in cosmx["reference"]:
        ref_rows.append({k: v for k, v in row.items() if k not in {"patient_high", "patient_low", "patient_ratio"}})
    _write_csv(
        os.path.join(TAB, "cosmx_reference_ratios.csv"),
        ref_rows,
        [
            "neighbor",
            "metric",
            "radius_um",
            "cut",
            "cell_weighted_ratio",
            "median_section_ratio",
            "equal_patient_ratio",
            "median_patient_log2",
            "lo",
            "hi",
        ],
    )
    if cosmx["winner"] is not None:
        w = cosmx["winner"]
        patient_rows = []
        for i, name in enumerate(cosmx["patient_names"]):
            patient_rows.append(
                {
                    "patient": name,
                    "mean_high": w["patient_high"][i],
                    "mean_low": w["patient_low"][i],
                    "delta_count_or_rate": w["patient_delta"][i],
                    "raw_log2": w["patient_raw_log2"][i],
                    "corrected_log2": w["patient_log2"][i],
                }
            )
        _write_csv(
            os.path.join(TAB, "cosmx_winner_patients.csv"),
            patient_rows,
            ["patient", "mean_high", "mean_low", "delta_count_or_rate", "raw_log2", "corrected_log2"],
        )
        if w.get("sections"):
            _write_csv(
                os.path.join(TAB, "cosmx_winner_sections.csv"),
                w["sections"],
                ["sample", "patient", "mean_high", "mean_low"],
            )
    c_fields = ["score", "cut", "estimand", "eligible", "median", "lo", "hi", "width", "magnitude", "n_boot_used"]
    _write_csv(os.path.join(TAB, "concordant4_grid.csv"), conc["rows"], c_fields)
    if conc["labels"] is not None:
        data = conc["data"]
        arm_rows = []
        for i in range(data["frac"].size):
            arm = conc["labels"][i]
            arm_rows.append(
                {
                    "dataset": data["cohort_name"][i],
                    "unit_id": data["unit"][i],
                    "tissue": data["tissue"][i],
                    "arm": {1: "high", 0: "low", -1: "mid"}[int(arm)],
                    "mal_CLDN4_pct": data["pct"][i],
                    "mal_CLDN4_mean": data["mean"][i],
                    "frac_tnk": data["frac"][i],
                }
            )
        _write_csv(
            os.path.join(TAB, "concordant4_winner_patients.csv"),
            arm_rows,
            ["dataset", "unit_id", "tissue", "arm", "mal_CLDN4_pct", "mal_CLDN4_mean", "frac_tnk"],
        )
        if conc["winner"] is not None and conc["winner"].get("cohorts"):
            _write_csv(
                os.path.join(TAB, "concordant4_winner_cohorts.csv"),
                conc["winner"]["cohorts"],
                ["dataset", "n_high", "n_low", "median_high", "median_low", "delta"],
            )

    def _jsonable(row: dict) -> dict:
        out = {}
        for k, v in row.items():
            if isinstance(v, np.ndarray):
                out[k] = [None if not np.isfinite(x) else float(x) for x in v.tolist()]
            elif isinstance(v, (np.floating, float)):
                out[k] = None if not np.isfinite(v) else float(v)
            elif isinstance(v, (np.integer, int)):
                out[k] = int(v)
            elif isinstance(v, (np.bool_, bool)):
                out[k] = bool(v)
            elif isinstance(v, dict):
                out[k] = _jsonable(v)
            elif isinstance(v, list):
                out[k] = [_jsonable(x) if isinstance(x, dict) else x for x in v]
            else:
                out[k] = v
        return out

    winner_c = None if cosmx["winner"] is None else _jsonable(cosmx["winner"])
    winner_k = None if conc["winner"] is None else _jsonable({k: v for k, v in conc["winner"].items() if k != "labels"})
    payload = {
        "seed": SEED,
        "n_perm": N_PERM,
        "n_boot": N_BOOT,
        "min_arm": MIN_ARM,
        "min_low": MIN_LOW,
        "eps": EPS,
        "cosmx_p_select": float(cosmx["p_select"]),
        "cosmx_nominal": None
        if cosmx["nominal"] is None
        else {key: float(val) for key, val in cosmx["nominal"].items()},
        "cosmx_winner": winner_c,
        "cosmx_n_tumor": cosmx["n_tumor"],
        "cosmx_n_cldn4_pos": cosmx["n_cldn4_pos"],
        "cosmx_null_max_quantiles": _quantiles(cosmx["null_max"]),
        "concordant4_p_select": float(conc["p_select"]),
        "concordant4_nominal_p": None if conc["nominal"] is None else float(conc["nominal"]["p"]),
        "concordant4_winner": winner_k,
        "concordant4_dl_rho": conc["dl"]["rho"],
        "concordant4_dl_p": conc["dl"]["p"],
        "concordant4_dl_i2": conc["dl"]["i2"],
        "concordant4_rho_ci": [conc["rho_lo"], conc["rho_hi"]],
        "concordant4_rho_perm_p": conc["p_rho"],
        "concordant4_cohort_rho": {
            str(conc["cohort_codes"][i]): {"rho": conc["dl"]["cohort_rho"][i], "n": conc["dl"]["cohort_n"][i]}
            for i in range(len(conc["cohort_codes"]))
        },
    }
    with open(os.path.join(OUT, "stats.json"), "w") as handle:
        json.dump(payload, handle, indent=2)
    np.save(os.path.join(TAB, "cosmx_null_max.npy"), cosmx["null_max"])
    np.save(os.path.join(TAB, "concordant4_null_max.npy"), conc["null_max"])
    _figures(cosmx, conc)
    text = _results_markdown(cosmx, conc)
    for path in (os.path.join(OUT, "RESULTS.md"), os.path.join(ROOT, "RESULTS.md")):
        with open(path, "w") as handle:
            handle.write(text)


def _quantiles(x: np.ndarray) -> dict:
    qs = np.quantile(x, [0.5, 0.9, 0.95, 1.0])
    return {"p50": float(qs[0]), "p90": float(qs[1]), "p95": float(qs[2]), "max": float(qs[3])}


def _figures(cosmx: dict, conc: dict) -> None:
    _forest_cosmx(cosmx)
    _grid_scatter(cosmx["rows"], os.path.join(FIG, "cosmx_grid.png"), "CosMx patient-median log2 ratio")
    _null_hist(
        cosmx["null_max"],
        0.0 if cosmx["winner"] is None else cosmx["winner"]["magnitude"],
        os.path.join(FIG, "cosmx_null.png"),
        "CosMx max |median log2| under within-FOV label swaps",
    )
    _forest_concordant(conc)
    _grid_scatter(conc["rows"], os.path.join(FIG, "concordant4_grid.png"), "Concordant-4 patient median Δ T/NK")
    _null_hist(
        conc["null_max"],
        0.0 if conc["winner"] is None else conc["winner"]["magnitude"],
        os.path.join(FIG, "concordant4_null.png"),
        "Concordant-4 max |median Δ| under within-cohort CLDN4 swaps",
    )


def _forest_cosmx(cosmx: dict) -> None:
    w = cosmx["winner"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if w is None:
        ax.text(0.5, 0.5, "No CosMx spec had a patient interval entirely below 0", ha="center")
    else:
        y = np.arange(len(cosmx["patient_names"]))
        ax.axvline(0, color="#444444", lw=0.8)
        ax.plot(w["patient_log2"], y, "o", color="#1f4e79", ms=8, label="Patient corrected log2 ratio")
        ax.plot([w["lo"], w["hi"]], [-1, -1], color="#8fa6bf", lw=3, solid_capstyle="butt")
        ax.plot([w["median"]], [-1], "D", color="#c45911", ms=7, label="Median and 95% patient interval")
        ax.set_yticks(list(y) + [-1])
        ax.set_yticklabels(list(cosmx["patient_names"]) + ["median"])
        ax.set_xlabel("log2 ( (CLDN4-high + ε) / (CLDN4-low + ε) )")
        ax.set_title(
            f"{w['neighbor']} {w['metric']} at {w['radius_um']:.0f} µm, {w['cut']} cut"
        )
        ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "cosmx_patient_forest.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "cosmx_patient_forest.pdf"))
    plt.close(fig)


def _grid_scatter(rows: list[dict], path: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for row in rows:
        if not row.get("eligible", False) or not np.isfinite(row["median"]):
            continue
        clears = row["hi"] < 0 and row["median"] < 0
        ax.scatter(
            row["width"],
            abs(row["median"]),
            s=28,
            c="#1f4e79" if clears else "#b0b0b0",
            zorder=3 if clears else 1,
        )
    winner = pick_winner(rows)
    if winner is not None:
        ax.scatter([winner["width"]], [abs(winner["median"])], s=90, facecolors="none", edgecolors="#c45911", lw=1.8, zorder=4)
    ax.set_xlabel("Patient-bootstrap 95% interval width")
    ax.set_ylabel("|median contrast|")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def _null_hist(null: np.ndarray, observed: float, path: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.hist(null, bins=30, color="#8fa6bf", edgecolor="white")
    ax.axvline(observed, color="#c45911", lw=1.8, label=f"Observed {_fmt(observed)}")
    ax.set_xlabel("Best exclusion |median| in the grid")
    ax.set_ylabel("Permutations")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def _forest_concordant(conc: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    w = conc["winner"]
    data = conc["data"]
    if w is None or conc["labels"] is None:
        ax.text(0.5, 0.5, "No concordant-4 spec had a patient interval entirely below 0", ha="center")
    else:
        lab = conc["labels"]
        rng = np.random.default_rng(1)
        colors = {"high": "#1f4e79", "low": "#c45911", "mid": "#b0b0b0"}
        for arm, name in ((1, "high"), (0, "low")):
            m = lab == arm
            y = data["frac"][m]
            x = np.full(y.size, 1 if arm == 1 else 0, dtype=float)
            x = x + rng.uniform(-0.12, 0.12, size=y.size)
            ax.scatter(x, y, s=18, c=colors[name], alpha=0.85, label=f"CLDN4 {name}")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4 low", "CLDN4 high"])
        ax.set_ylabel("Patient T/NK fraction")
        ax.set_title(
            f"{w['score']} {w['cut']} / {w['estimand']}\n"
            f"median Δ {_fmt(w['median'])}  [{_fmt(w['lo'])}, {_fmt(w['hi'])}]"
        )
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "concordant4_patients.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "concordant4_patients.pdf"))
    plt.close(fig)


def _results_markdown(cosmx: dict, conc: dict) -> str:
    w = cosmx["winner"]
    k = conc["winner"]
    lines = []
    lines.append("# Patient bootstrap of CosMx exclusion and concordant-4 CLDN4–T/NK")
    lines.append("")
    lines.append(
        "Additive layer on the locked public results. The CosMx cytotoxic ratios "
        "0.36 / 0.52 (50 / 100 µm, 8/8 sections, 5/5 patients, sign P = 0.031) are not recomputed "
        "as a replacement, and the concordant-4 Spearman is not replaced. "
        "Inference here uses patients only."
    )
    lines.append("")
    lines.append("## How a spec was allowed to win")
    lines.append("")
    lines.append(
        "The grid below was fixed before either outcome was ranked. "
        "A spec is eligible only when every patient has a usable contrast. "
        "The patient-bootstrap 95% interval of the median contrast must lie entirely below zero "
        "(CLDN4-high is the lower arm). Among those specs, the largest absolute median wins. "
        "If two medians match, the narrower interval wins."
    )
    lines.append("")
    lines.append(
        "With five CosMx patients the percentile bootstrap of the median is exactly the range "
        "of the five patient values: the upper bound is the weakest patient. "
        "That is the whole-interval clearance, not a large-sample standard error."
    )
    lines.append("")
    lines.append(
        f"CosMx log2 ratios use a fixed continuity correction ε "
        f"(count {EPS['count']}, contact rate {EPS['contact']}, neighbor fraction {EPS['fraction']}). "
        f"The reference arm must be at least "
        f"{MIN_LOW['count']} / {MIN_LOW['contact']} / {MIN_LOW['fraction']} "
        "for those three metrics. Sections need at least "
        f"{MIN_ARM} tumor cells in each arm. "
        f"Permutations: {N_PERM}. Bootstrap draws for concordant-4 and for the five-patient identity check: {N_BOOT}. "
        f"Seed {SEED}."
    )
    lines.append("")
    lines.append("## CosMx")
    lines.append("")
    lines.append(
        f"Object: figshare 25976224, {cosmx['n_tumor']} author tumor cells, "
        f"{cosmx['n_cldn4_pos']} with CLDN4 count > 0, 8 sections, 5 patients. "
        "Neighbors are counted inside the same FOV. Coordinates are CosMx pixels at 0.18 µm/pixel. "
        "Cytotoxic cells are author NK plus CD8 T cells. T/NK adds CD4 and Treg. "
        "Immune adds B, myeloid, mast, and plasmablast. "
        "CLDN4 cuts, inside each section: detected versus undetected, top half versus bottom half, "
        "top quartile versus bottom quartile. Ties break on the cell index."
    )
    lines.append("")
    lines.append(
        "Patient value = unweighted mean of that patient's section means, then "
        "log2((high + ε) / (low + ε)). The reported median is the median of the five patients. "
        "Within-FOV permutation swaps CLDN4 arms inside each section's FOV and leaves coordinates "
        "and neighbor counts where they are."
    )
    lines.append("")
    if w is None:
        lines.append("No spec had all five patients below zero.")
    else:
        lines.append(
            f"Winning spec: **{w['neighbor']} {w['metric']}**, **{w['radius_um']:.0f} µm**, "
            f"**{w['cut']}** CLDN4 cut."
        )
        lines.append("")
        lines.append(
            f"Median corrected log2 ratio **{_fmt(w['median'])}** "
            f"(95% patient interval **[{_fmt(w['lo'])}, {_fmt(w['hi'])}]**, "
            f"width {_fmt(w['width'])}). "
            f"Equal-patient-weight ratio of means **{_fmt(w['equal_patient_ratio'])}**."
        )
        lines.append("")
        lines.append(
            f"Selection-adjusted permutation p = **{_fmt(cosmx['p_select'], 4)}** "
            f"({N_PERM} within-FOV shuffles; the null statistic is the best |median| in the same "
            f"{len(cosmx['rows'])}-spec grid). "
            f"Nominal permutation p for this spec alone = **{_fmt(w.get('nominal_p'), 4)}**. "
            "The nominal p does not pay for the search."
        )
        lines.append("")
        lines.append("| Patient | High-arm mean | Low-arm mean | High − low | Corrected log2 |")
        lines.append("|---|---:|---:|---:|---:|")
        for i, name in enumerate(cosmx["patient_names"]):
            lines.append(
                f"| {name} | {_fmt(w['patient_high'][i])} | {_fmt(w['patient_low'][i])} | "
                f"{_fmt(w['patient_delta'][i])} | {_fmt(w['patient_log2'][i])} |"
            )
        lines.append("")
        lines.append(
            "Patient means are the unweighted average of section means. "
            "For a fraction, the mean is taken over tumor cells that have at least one neighbor inside the radius."
        )
        lines.append("")
        if w.get("sections"):
            n_down = sum(
                np.isfinite(s["mean_high"]) and np.isfinite(s["mean_low"]) and s["mean_low"] > 0 and s["mean_high"] < s["mean_low"]
                for s in w["sections"]
            )
            lines.append(f"Section means for the same spec ({n_down}/8 with high < low):")
            lines.append("")
            lines.append("| Section | Patient | High | Low | High/low |")
            lines.append("|---|---|---:|---:|---:|")
            for s in w["sections"]:
                ratio = s["mean_high"] / s["mean_low"] if s["mean_low"] not in (0, None) and np.isfinite(s["mean_low"]) and s["mean_low"] != 0 else math_nan()
                lines.append(
                    f"| {s['sample']} | {s['patient']} | {_fmt(s['mean_high'])} | {_fmt(s['mean_low'])} | {_fmt(ratio)} |"
                )
    lines.append("")
    lines.append("### Locked-family reference (cytotoxic neighbor counts)")
    lines.append("")
    lines.append(
        "These rows are the 50 µm and 100 µm cytotoxic counts under every CLDN4 cut. "
        "The cell-weighted ratio pools tumor cells and is not a patient-level estimate. "
        "It is printed so it can be compared with the locked 0.36 / 0.52 ratios without substituting for them."
    )
    lines.append("")
    lines.append("| Radius | Cut | Cell-weighted ratio | Median section ratio | Equal-patient ratio | Patient median log2 | Patient interval |")
    lines.append("|---:|---|---:|---:|---:|---:|---|")
    for row in cosmx["reference"]:
            lines.append(
            f"| {row['radius_um']:.0f} | {row['cut']} | {_fmt(row['cell_weighted_ratio'])} | "
            f"{_fmt(row['median_section_ratio'])} | {_fmt(row['equal_patient_ratio'])} | "
            f"{_fmt(row['median_patient_log2'])} | [{_fmt(row['lo'])}, {_fmt(row['hi'])}] |"
        )
    lines.append("")
    lines.append(
        "Under this author-label, within-FOV, equal-patient definition, the 50 µm cytotoxic "
        "median split has equal-patient ratio 0.840 and a patient interval that just clears zero. "
        "The 100 µm intervals include zero. Those ratios are not the locked 0.36 / 0.52 summary."
    )
    lines.append("")
    top = sorted(
        [r for r in cosmx["rows"] if r["magnitude"] > 0],
        key=lambda r: (-r["magnitude"], r["width"]),
    )[:8]
    lines.append("### Next exclusion specs by |median|")
    lines.append("")
    lines.append("| Neighbor | Metric | µm | Cut | Median log2 | Interval | Width | Equal-patient ratio |")
    lines.append("|---|---|---:|---|---:|---|---:|---:|")
    for row in top:
        lines.append(
            f"| {row['neighbor']} | {row['metric']} | {row['radius_um']:.0f} | {row['cut']} | "
            f"{_fmt(row['median'])} | [{_fmt(row['lo'])}, {_fmt(row['hi'])}] | {_fmt(row['width'])} | "
            f"{_fmt(row['equal_patient_ratio'])} |"
        )
    lines.append("")
    lines.append("## Concordant-4")
    lines.append("")
    lines.append(
        "Patient table from the four-cohort integration (65 units: GSE123902 donors, GSE131907 samples, "
        "GSE205335 patients, GSE189357 patients). Outcome is `frac_tnk = n_tnk / n_cells`. "
        "CLDN4 scores are malignant % positive and malignant mean. "
        "Cuts are formed inside each cohort (median halves, tertile extremes, quartile extremes). "
        "Estimands: median of within-cohort pairwise high−low differences; pooled difference of medians; "
        "median of the four cohort differences. The bootstrap resamples patients inside each cohort "
        "and rebuilds the cut."
    )
    lines.append("")
    dl = conc["dl"]
    lines.append(
        f"Locked estimand, recomputed from the same table: DerSimonian–Laird Spearman "
        f"**ρ = {_fmt(dl['rho'], 3)}** (two-sided p = {_fmt(dl['p'], 4)}, I² = {_fmt(100 * dl['i2'], 1)}%). "
        f"Analytical 95% interval on the Fisher-z scale: **[{_fmt(dl['ci_lo'])}, {_fmt(dl['ci_hi'])}]**. "
        f"Stratified patient-bootstrap interval: [{_fmt(conc['rho_lo'])}, {_fmt(conc['rho_hi'])}]. "
        f"The bootstrap lower tail is wider because a resample of the n=9 cohort can sit on the correlation boundary. "
        f"Within-cohort permutation p (one-sided, more negative) = **{_fmt(conc['p_rho'], 4)}**."
    )
    lines.append("")
    if k is None:
        lines.append("No median-Δ spec had a patient-bootstrap interval entirely below zero.")
    else:
        lines.append(
            f"Winning median-Δ spec: **{k['score']}**, **{k['cut']}** cut, **{k['estimand']}**."
        )
        lines.append("")
        lines.append(
            f"Median Δ (CLDN4-high − CLDN4-low T/NK fraction) **{_fmt(k['median'])}** "
            f"(95% stratified patient-bootstrap interval **[{_fmt(k['lo'])}, {_fmt(k['hi'])}]**, "
            f"width {_fmt(k['width'])})."
        )
        lines.append("")
        lines.append(
            f"Selection-adjusted permutation p = **{_fmt(conc['p_select'], 4)}** "
            f"(within-cohort CLDN4 shuffles; null statistic = best |median Δ| among the "
            f"{len(conc['rows'])} specs). "
            f"Nominal permutation p for this spec alone = **{_fmt(k.get('nominal_p'), 4)}**."
        )
        if k.get("cohorts"):
            lines.append("")
            lines.append("| Cohort | n high | n low | Median T/NK high | Median T/NK low | Δ |")
            lines.append("|---|---:|---:|---:|---:|---:|")
            for c in k["cohorts"]:
                lines.append(
                    f"| {c['dataset']} | {c['n_high']} | {c['n_low']} | {_fmt(c['median_high'])} | "
                    f"{_fmt(c['median_low'])} | {_fmt(c['delta'])} |"
                )
            lines.append("")
            lines.append(
                "The winning Δ is the median of those four cohort differences. "
                "GSE189357 has nine patients, so the quartile arms are two and two."
            )
    lines.append("")
    lines.append("| Score | Cut | Estimand | Median Δ | 95% interval | Width | Eligible |")
    lines.append("|---|---|---|---:|---|---:|---|")
    for row in conc["rows"]:
        lines.append(
            f"| {row['score']} | {row['cut']} | {row['estimand']} | {_fmt(row['median'])} | "
            f"[{_fmt(row['lo'])}, {_fmt(row['hi'])}] | {_fmt(row['width'])} | {row['eligible']} |"
        )
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append(
        "Sections are not five extra patients. Cells are not the sample size. "
        "The search picked the contrast, so the interval on the winner is the interval for that selected contrast; "
        "the permutation p-value is the one that includes the search. "
        "A Visium spot correlation is not in this file. Public mouse Harmony objects are not merged with the private KL cohort. "
        "Nearby effector cells are not claimed to be muzzled."
    )
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/download_cosmx_nsclc_h5ad.py")
    lines.append("python3 scripts/run_patient_separation.py")
    lines.append("```")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    print("concordant-4", flush=True)
    conc = concordant_grid(load_concordant(PATIENT_TSV), rng)
    print("concordant winner", None if conc["winner"] is None else {k: conc["winner"][k] for k in ("score", "cut", "estimand", "median", "lo", "hi", "magnitude")}, flush=True)
    print("cosmx extract", flush=True)
    cache = extract_tumor_neighbors(H5AD, CACHE)
    print("cosmx grid", flush=True)
    cosmx = cosmx_grid(cache, rng)
    print("cosmx winner", None if cosmx["winner"] is None else {k: cosmx["winner"][k] for k in ("neighbor", "metric", "radius_um", "cut", "median", "lo", "hi", "magnitude")}, flush=True)
    save_outputs(cosmx, conc)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
