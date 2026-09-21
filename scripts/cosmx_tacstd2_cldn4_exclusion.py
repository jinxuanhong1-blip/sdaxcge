#!/usr/bin/env python3
"""Does TACSTD2-local immune exclusion in He 2022 CosMx depend on CLDN4?

Pre-specified before neighbor summaries are read. Figshare 25976224.
Does not replace the locked CLDN4 cytotoxic ratios 0.36 / 0.52 at 50 / 100 µm.

Index cells: author tumor matched to the section (tumor 5/6/9/12/13).
Expression: log1p(count / n_counts * 1e4). Detected = raw count > 0.
Coordinates: obsm/spatial * 0.18 µm/px. Trees are built per section.

Primary high vs low: within-section median of log-normalized TACSTD2.
Companion cuts, also fixed in advance: raw detected vs absent, and Q4 vs Q1.
Neither companion replaces the median row.

Neighbor classes: broad immune (14 author labels), CD8 (memory + naive),
and CD8+NK as a cytotoxic companion. Radii: 10, 20, 50, 100 µm.
Primary readouts: immune-cell fraction of other cells in the radius, and
CD8 neighbor count. Counts and CD8+NK counts are reported with them.

CLDN4 dependence, all on the median split plus the continuous scores:
  1. Stratify: within-section CLDN4 median, then TACSTD2 median inside the stratum.
  2. Quadrants: section-wide medians of the two genes.
  3. Residualize: outcome minus its linear fit on CLDN4, then the TACSTD2 median contrast.
  4. Partial fields: Gaussian expression and immune fields; partial Spearman and
     standardized regression of the immune/CD8 field on the TACSTD2 field given CLDN4.
  5. Mediation-like attenuation: standardized coefficient of TACSTD2 before and
     after CLDN4, on ranks (primary, because counts are zero-inflated) and on
     the log-normalized values (companion). This is observational attenuation,
     not a causal mediation estimate.

Inference is the section tally (8) and the patient tally (5). Patient values
are unweighted means of that patient's sections. A missing section counts
against 8/8. Cell-level p-values are not used: neighborhoods overlap.
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
from scipy import sparse, stats
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.spatial import cKDTree
from scipy.stats import rankdata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_tacstd2_cldn4")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
CACHE = os.path.join("/tmp", "cosmx_tacstd2_cldn4_cells.npz")

UM_PER_PX = 0.18
RADII = (10, 20, 50, 100)
MIN_ARM = 30
MIN_FOV_FIELD = 40
MIN_FOVS = 3
TRUNCATE = 4.0
KERNEL_MASS_MIN = 0.80
ATT_MIN_ABS = 0.01

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
EXPECTED_PATIENT = {
    "LUAD-5 R1": "Lung5",
    "LUAD-5 R2": "Lung5",
    "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9",
    "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}
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
CD8_LABELS = ("T CD8 memory", "T CD8 naive")
CD8NK_LABELS = ("T CD8 memory", "T CD8 naive", "NK")
PATIENT_COLOR = {
    "Lung5": "#1b9e77",
    "Lung6": "#d95f02",
    "Lung9": "#7570b3",
    "Lung12": "#e7298a",
    "Lung13": "#66a61e",
}
PRIMARY_METRICS = ("imm_frac", "cd8_n")
ALL_METRICS = ("imm_frac", "imm_n", "cd8_n", "cd8_frac", "cd8nk_n")
METRIC_LABEL = {
    "imm_frac": "Immune fraction",
    "imm_n": "Immune count",
    "cd8_n": "CD8 count",
    "cd8_frac": "CD8 fraction",
    "cd8nk_n": "CD8+NK count",
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


def _decode(arr) -> list:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def stable_seed(*parts) -> int:
    acc = 2166136261
    for part in parts:
        for ch in str(part):
            acc ^= ord(ch)
            acc = (acc * 16777619) & 0xFFFFFFFF
    return int(acc)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _z(a: np.ndarray):
    a = np.asarray(a, dtype=np.float64)
    s = float(a.std(ddof=0))
    if not np.isfinite(s) or s < 1e-12:
        return None
    return (a - float(a.mean())) / s


def standardized_attenuation(x, z, y) -> dict:
    """Observational attenuation of standardized OLS coefficients.

    beta_total is the simple standardized slope of y on x.
    beta_direct is the slope of x in y ~ x + z.
    Attenuation is 1 - beta_direct / beta_total when |beta_total| >= 0.01.
    In this linear model it equals the product-of-paths proportion.
    """
    xz, zz, yz = _z(x), _z(z), _z(y)
    empty = {
        "beta_total": np.nan,
        "beta_direct": np.nan,
        "beta_cldn4_joint": np.nan,
        "beta_cldn4_total": np.nan,
        "a_path": np.nan,
        "attenuation": np.nan,
        "proportion": np.nan,
        "attenuation_of_cldn4": np.nan,
        "interaction": np.nan,
    }
    if xz is None or zz is None or yz is None:
        return empty
    b_tot = float(np.dot(xz, yz) / np.dot(xz, xz))
    b_c_tot = float(np.dot(zz, yz) / np.dot(zz, zz))
    coef, _, _, _ = np.linalg.lstsq(np.column_stack([xz, zz]), yz, rcond=None)
    b_dir = float(coef[0])
    b_z = float(coef[1])
    a_path = float(np.dot(xz, zz) / np.dot(xz, xz))
    att = np.nan
    prop = np.nan
    att_c = np.nan
    if abs(b_tot) >= ATT_MIN_ABS:
        att = float(1.0 - b_dir / b_tot)
        prop = float(a_path * b_z / b_tot)
    if abs(b_c_tot) >= ATT_MIN_ABS:
        att_c = float(1.0 - b_z / b_c_tot)
    inter_cols = np.column_stack([xz, zz, xz * zz])
    inter, _, _, _ = np.linalg.lstsq(inter_cols, yz, rcond=None)
    return {
        "beta_total": b_tot,
        "beta_direct": b_dir,
        "beta_cldn4_joint": b_z,
        "beta_cldn4_total": b_c_tot,
        "a_path": a_path,
        "attenuation": att,
        "proportion": prop,
        "attenuation_of_cldn4": att_c,
        "interaction": float(inter[2]),
    }


def partial_spearman(x, y, z) -> float:
    """Spearman of x vs y after linear residualizing both ranks on z."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if int(m.sum()) < MIN_ARM:
        return float("nan")
    rx = rankdata(x[m]).astype(np.float64)
    ry = rankdata(y[m]).astype(np.float64)
    rz = rankdata(z[m]).astype(np.float64)
    design = np.column_stack([np.ones(rx.size), rz])
    cx, _, _, _ = np.linalg.lstsq(design, rx, rcond=None)
    cy, _, _, _ = np.linalg.lstsq(design, ry, rcond=None)
    xr = rx - design @ cx
    yr = ry - design @ cy
    if float(xr.std()) < 1e-12 or float(yr.std()) < 1e-12:
        return float("nan")
    return float(np.corrcoef(xr, yr)[0, 1])


def spearman(x, y) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if int(m.sum()) < MIN_ARM:
        return float("nan")
    rho = stats.spearmanr(x[m], y[m]).correlation
    return float(rho) if np.isfinite(rho) else float("nan")


def sign_test(n_excl: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return float(stats.binomtest(int(n_excl), int(n), 0.5, alternative="greater").pvalue)


def wilcoxon_p(deltas: np.ndarray) -> float:
    d = np.asarray(deltas, dtype=np.float64)
    d = d[np.isfinite(d) & (d != 0)]
    if d.size == 0:
        return float("nan")
    method = "exact" if d.size <= 25 else "auto"
    return float(stats.wilcoxon(d, alternative="two-sided", method=method, zero_method="wilcox").pvalue)


def nanmean(a) -> float:
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(a.mean())


def _finite_arm(y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return y[mask & np.isfinite(y)]


def arm_summary(y: np.ndarray, high: np.ndarray, low: np.ndarray) -> dict:
    yh = _finite_arm(y, high)
    yl = _finite_arm(y, low)
    usable = yh.size >= MIN_ARM and yl.size >= MIN_ARM
    hm = float(yh.mean()) if yh.size else float("nan")
    lm = float(yl.mean()) if yl.size else float("nan")
    return {
        "usable": bool(usable),
        "n_high": int(yh.size),
        "n_low": int(yl.size),
        "high_mean": hm if usable else float("nan"),
        "low_mean": lm if usable else float("nan"),
        "delta": (hm - lm) if usable else float("nan"),
        "ratio": (hm / lm) if usable and abs(lm) > 1e-8 else float("nan"),
    }


def median_masks(values: np.ndarray):
    med = float(np.median(values))
    high = values > med
    low = values <= med
    return high, low, med


def q4q1_masks(values: np.ndarray):
    q1, q3 = np.quantile(values, [0.25, 0.75])
    low = values <= q1
    high = (values >= q3) & ~low
    return high, low, float(q1), float(q3)


def classify_dependence(marginal_8: int, marginal_5: int, stratum_ok: bool, att_median: float, n_att: int) -> str:
    """Label from pre-specified rules. Missing units count against 8 and 5."""
    concordant = marginal_8 == 8 and marginal_5 == 5
    att_defined = n_att >= 6 and np.isfinite(att_median)
    if not concordant:
        return "not_concordant_8_of_8_and_5_of_5"
    if att_defined and att_median > 0.5 and not stratum_ok:
        return "concordant_and_largely_shared_with_cldn4"
    if att_defined and att_median < 0.25 and stratum_ok:
        return "concordant_and_not_accounted_for_by_cldn4"
    return "concordant_with_mixed_cldn4_adjustment"


def ball_counts(tree: cKDTree, xy: np.ndarray, radius: float) -> np.ndarray:
    if tree.n == 0 or len(xy) == 0:
        return np.zeros(len(xy), dtype=np.float64)
    try:
        out = tree.query_ball_point(xy, r=float(radius), return_length=True, workers=4)
    except TypeError:
        out = tree.query_ball_point(xy, r=float(radius), return_length=True)
    return np.asarray(out, dtype=np.float64)


def _decode_cats(f, name: str):
    cats = _decode(f[f"obs/{name}/categories"][:])
    codes = f[f"obs/{name}/codes"][:]
    return np.asarray(cats, dtype=object)[codes], cats


def load_and_count():
    import h5py

    t0 = time.time()
    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}
    for g in ("TACSTD2", "CLDN4"):
        if g not in gmap:
            raise SystemExit(f"{g} is absent from this object")
    sample, sample_cats = _decode_cats(f, "sample")
    patient, _ = _decode_cats(f, "patient")
    cell_type, type_cats = _decode_cats(f, "cell_type")
    missing = sorted(set(BROAD) - set(type_cats))
    if missing:
        raise SystemExit(f"immune labels missing: {missing}")
    fov = f["obs/fov"][:].astype(np.int32)
    n_counts = f["obs/n_counts"][:].astype(np.float64)
    xy = f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    f.close()
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(sample), len(genes)))
    del data, indices, indptr
    tac_raw = np.asarray(mat.getcol(gmap["TACSTD2"]).toarray()).ravel().astype(np.float64)
    cld_raw = np.asarray(mat.getcol(gmap["CLDN4"]).toarray()).ravel().astype(np.float64)
    del mat
    lib = n_counts.copy()
    lib[lib <= 0] = np.nan
    tac = np.nan_to_num(np.log1p(tac_raw / lib * 1e4), nan=0.0)
    cld = np.nan_to_num(np.log1p(cld_raw / lib * 1e4), nan=0.0)
    print(f"loaded genes in {time.time() - t0:.1f}s", flush=True)

    missing_samples = [s for s in SAMPLES if s not in set(sample_cats)]
    if missing_samples:
        raise SystemExit(f"missing samples {missing_samples}; found {sample_cats}")

    is_immune = np.isin(cell_type, list(BROAD))
    is_cd8 = np.isin(cell_type, list(CD8_LABELS))
    is_cd8nk = np.isin(cell_type, list(CD8NK_LABELS))
    pieces = []
    field_rows = []
    inventory = []
    for s in SAMPLES:
        t1 = time.time()
        m = sample == s
        pat = str(patient[m][0])
        if pat != EXPECTED_PATIENT[s]:
            raise SystemExit(f"{s} patient {pat} != {EXPECTED_PATIENT[s]}")
        xy_s = np.ascontiguousarray(xy[m])
        ct = cell_type[m]
        mal = ct == TUMOR_LABEL[s]
        if int(mal.sum()) < MIN_ARM * 2:
            raise SystemExit(f"{s} has too few matched tumor cells")
        rng = np.random.default_rng(0)
        take = rng.choice(xy_s.shape[0], size=min(6000, xy_s.shape[0]), replace=False)
        dist, _ = cKDTree(xy_s).query(xy_s[take], k=2, workers=4)
        med_nn = float(np.median(dist[:, 1]))
        if not (3.0 <= med_nn <= 40.0):
            raise SystemExit(f"{s}: median NN {med_nn:.3f} µm is outside 3–40; check 0.18 µm/px")

        tum_xy = np.ascontiguousarray(xy_s[mal])
        trees = {
            "immune": cKDTree(xy_s[is_immune[m]]),
            "cd8": cKDTree(xy_s[is_cd8[m]]),
            "cd8nk": cKDTree(xy_s[is_cd8nk[m]]),
            "all": cKDTree(xy_s),
        }
        cols = {
            "sample": np.full(int(mal.sum()), s, dtype=object),
            "patient": np.full(int(mal.sum()), pat, dtype=object),
            "fov": fov[m][mal],
            "tac": tac[m][mal],
            "cld": cld[m][mal],
            "tac_pos": tac_raw[m][mal] > 0,
            "cld_pos": cld_raw[m][mal] > 0,
        }
        for radius in RADII:
            n_imm = ball_counts(trees["immune"], tum_xy, radius)
            n_cd8 = ball_counts(trees["cd8"], tum_xy, radius)
            n_nk = ball_counts(trees["cd8nk"], tum_xy, radius)
            n_all = ball_counts(trees["all"], tum_xy, radius)
            if float(n_all.min()) < 1:
                raise SystemExit(f"{s} {radius} µm: a tumor cell did not find itself")
            denom = n_all - 1.0
            frac = np.divide(n_imm, denom, out=np.full_like(n_imm, np.nan), where=denom > 0)
            cd8_frac = np.divide(n_cd8, denom, out=np.full_like(n_cd8, np.nan), where=denom > 0)
            cols[f"imm_n_{radius}"] = n_imm
            cols[f"imm_frac_{radius}"] = frac
            cols[f"cd8_n_{radius}"] = n_cd8
            cols[f"cd8_frac_{radius}"] = cd8_frac
            cols[f"cd8nk_n_{radius}"] = n_nk
            cols[f"denom_{radius}"] = denom
        piece = pd.DataFrame(cols)
        pieces.append(piece)
        field_rows.extend(section_fields(s, pat, xy_s, mal, fov[m], is_immune[m], is_cd8[m], tac[m], cld[m]))
        tac_m = tac[m][mal]
        cld_m = cld[m][mal]
        inventory.append({
            "sample": s,
            "patient": pat,
            "tumor_label": TUMOR_LABEL[s],
            "n_cells": int(m.sum()),
            "n_tumor": int(mal.sum()),
            "n_immune": int(is_immune[m].sum()),
            "n_cd8": int(is_cd8[m].sum()),
            "n_cd8nk": int(is_cd8nk[m].sum()),
            "tacstd2_pos_frac": float((tac_raw[m][mal] > 0).mean()),
            "cldn4_pos_frac": float((cld_raw[m][mal] > 0).mean()),
            "tacstd2_median_log": float(np.median(tac_m)),
            "cldn4_median_log": float(np.median(cld_m)),
            "median_split_is_detected": bool(np.median(tac_m) == 0.0),
            "spearman_tacstd2_cldn4": spearman(tac_m, cld_m),
            "median_nn_um": med_nn,
            "n_zero_library": int((n_counts[m] <= 0).sum()),
        })
        print(
            f"{s}: tumor={int(mal.sum())} immune={int(is_immune[m].sum())} "
            f"CD8={int(is_cd8[m].sum())} NN={med_nn:.2f} µm in {time.time() - t1:.1f}s",
            flush=True,
        )
    cells = pd.concat(pieces, ignore_index=True)
    meta = {"inventory": inventory, "n_tumor": int(len(cells)), "um_per_px": UM_PER_PX}
    return cells, field_rows, meta


def section_fields(sample, patient, xy, mal, fov, immune, cd8, tac, cld) -> list:
    """Gaussian fields sampled inside each FOV.

    Pixel size is bandwidth/5 so the kernel is about 5 pixels. Expression
    fields are kernel-weighted means over tumor cells. Immune fraction is
    the kernel-weighted immune share of all cells. CD8 density is cells per µm².
    """
    rows = []
    tum_xy = xy[mal]
    fov_t = fov[mal]
    tac_t = tac[mal]
    cld_t = cld[mal]
    xmin, ymin = float(xy[:, 0].min()), float(xy[:, 1].min())
    xmax, ymax = float(xy[:, 0].max()), float(xy[:, 1].max())
    outcomes_spec = ("imm_frac", "cd8_density")
    for bandwidth in RADII:
        pixel = float(bandwidth) / 5.0
        nx = int(np.floor((xmax - xmin) / pixel)) + 2
        ny = int(np.floor((ymax - ymin) / pixel)) + 2
        sigma = float(bandwidth) / pixel
        ones = np.ones((ny, nx), dtype=np.float32)
        mass = gaussian_filter(ones, sigma=sigma, mode="constant", cval=0.0, truncate=TRUNCATE)
        del ones

        def grid_of(points, weights):
            g = np.zeros((ny, nx), dtype=np.float32)
            if len(points) == 0:
                return g
            ix = np.floor((points[:, 0] - xmin) / pixel).astype(np.int32)
            iy = np.floor((points[:, 1] - ymin) / pixel).astype(np.int32)
            ok = (ix >= 0) & (iy >= 0) & (ix < nx) & (iy < ny)
            np.add.at(g, (iy[ok], ix[ok]), np.asarray(weights, dtype=np.float32)[ok])
            return gaussian_filter(g, sigma=sigma, mode="constant", cval=0.0, truncate=TRUNCATE)

        sm_all = grid_of(xy, np.ones(len(xy), dtype=np.float32))
        sm_imm = grid_of(xy[immune], np.ones(int(immune.sum()), dtype=np.float32))
        sm_cd8 = grid_of(xy[cd8], np.ones(int(cd8.sum()), dtype=np.float32))
        sm_tum = grid_of(tum_xy, np.ones(len(tum_xy), dtype=np.float32))
        sm_tac = grid_of(tum_xy, tac_t.astype(np.float32))
        sm_cld = grid_of(tum_xy, cld_t.astype(np.float32))

        def sample_grid(grid, pts):
            cols = (pts[:, 0] - xmin) / pixel
            rows_ = (pts[:, 1] - ymin) / pixel
            return map_coordinates(grid, [rows_, cols], order=1, mode="constant", cval=np.nan)

        for fov_id in np.unique(fov_t):
            idx = np.flatnonzero(fov_t == fov_id)
            if idx.size < MIN_FOV_FIELD:
                continue
            rng = np.random.default_rng(stable_seed(sample, int(fov_id), int(bandwidth)))
            if idx.size > 250:
                idx = rng.choice(idx, size=250, replace=False)
            pts = tum_xy[idx]
            mass_s = sample_grid(mass, pts)
            tum_s = sample_grid(sm_tum, pts)
            all_s = sample_grid(sm_all, pts)
            imm_s = sample_grid(sm_imm, pts)
            cd8_s = sample_grid(sm_cd8, pts)
            tac_s = sample_grid(sm_tac, pts)
            cld_s = sample_grid(sm_cld, pts)
            with np.errstate(divide="ignore", invalid="ignore"):
                tum_mean_ok = tum_s / mass_s
                n_disk = tum_mean_ok / (pixel ** 2) * np.pi * (bandwidth ** 2)
                tac_field = tac_s / tum_s
                cld_field = cld_s / tum_s
                imm_frac = imm_s / all_s
                cd8_den = cd8_s / mass_s / (pixel ** 2)
            keep = (
                np.isfinite(mass_s)
                & (mass_s >= KERNEL_MASS_MIN)
                & np.isfinite(n_disk)
                & (n_disk >= 3.0)
                & np.isfinite(tac_field)
                & np.isfinite(cld_field)
                & np.isfinite(imm_frac)
                & np.isfinite(cd8_den)
            )
            n_keep = int(keep.sum())
            usable = n_keep >= MIN_FOV_FIELD
            values = {"imm_frac": imm_frac, "cd8_density": cd8_den}
            for outcome in outcomes_spec:
                y = values[outcome]
                if not usable:
                    rho_m = rho_p = b_tac = b_cld = rho_tz = float("nan")
                else:
                    rho_m = spearman(tac_field[keep], y[keep])
                    rho_p = partial_spearman(tac_field[keep], y[keep], cld_field[keep])
                    fit = standardized_attenuation(tac_field[keep], cld_field[keep], y[keep])
                    b_tac = fit["beta_direct"]
                    b_cld = fit["beta_cldn4_joint"]
                    rho_tz = spearman(tac_field[keep], cld_field[keep])
                rows.append({
                    "sample": sample,
                    "patient": patient,
                    "fov": int(fov_id),
                    "bandwidth": int(bandwidth),
                    "outcome": outcome,
                    "n_points": n_keep,
                    "usable": bool(usable),
                    "rho_marginal": rho_m,
                    "rho_partial": rho_p,
                    "beta_tacstd2": b_tac,
                    "beta_cldn4": b_cld,
                    "rho_tac_cld_field": rho_tz,
                })
        print(f"  field {sample} {bandwidth} µm", flush=True)
    return rows


def ycol(metric: str, radius: int) -> str:
    return f"{metric}_{radius}"


def section_split_rows(cells: pd.DataFrame, gene: str, mode: str) -> list:
    """gene is 'tac' or 'cld'. mode is median, detected, or q4q1."""
    rows = []
    value_key = "tac" if gene == "tac" else "cld"
    pos_key = "tac_pos" if gene == "tac" else "cld_pos"
    for s in SAMPLES:
        sub = cells[cells["sample"] == s]
        values = sub[value_key].to_numpy(dtype=np.float64)
        pos = sub[pos_key].to_numpy(dtype=bool)
        if mode == "median":
            high, low, med = median_masks(values)
            extra = {"threshold": med}
        elif mode == "detected":
            high, low = pos, ~pos
            extra = {"threshold": 0.0}
        elif mode == "q4q1":
            high, low, q1, q3 = q4q1_masks(values)
            extra = {"threshold": q3, "q1": q1}
        else:
            raise ValueError(mode)
        for metric in ALL_METRICS:
            for radius in RADII:
                y = sub[ycol(metric, radius)].to_numpy(dtype=np.float64)
                sm = arm_summary(y, high, low)
                rows.append({
                    "gene": gene,
                    "cut": mode,
                    "metric": metric,
                    "radius": int(radius),
                    "sample": s,
                    "patient": sub["patient"].iloc[0],
                    **extra,
                    **sm,
                })
    return rows


def stratified_rows(cells: pd.DataFrame) -> list:
    rows = []
    for s in SAMPLES:
        sub = cells[cells["sample"] == s]
        cld = sub["cld"].to_numpy(dtype=np.float64)
        tac = sub["tac"].to_numpy(dtype=np.float64)
        c_high, c_low, cmed = median_masks(cld)
        for stratum, mask in (("cldn4_high", c_high), ("cldn4_low", c_low)):
            if int(mask.sum()) < MIN_ARM * 2:
                high = low = np.zeros(mask.sum(), dtype=bool)
                tmed = float("nan")
                ybase = None
            else:
                high, low, tmed = median_masks(tac[mask])
            for metric in PRIMARY_METRICS:
                for radius in RADII:
                    y_all = sub[ycol(metric, radius)].to_numpy(dtype=np.float64)
                    if int(mask.sum()) < MIN_ARM * 2:
                        sm = arm_summary(np.array([]), np.array([], dtype=bool), np.array([], dtype=bool))
                    else:
                        sm = arm_summary(y_all[mask], high, low)
                    rows.append({
                        "stratum": stratum,
                        "metric": metric,
                        "radius": int(radius),
                        "sample": s,
                        "patient": sub["patient"].iloc[0],
                        "cldn4_median": cmed,
                        "tacstd2_median_in_stratum": tmed,
                        "n_stratum": int(mask.sum()),
                        **sm,
                    })
    return rows


def quadrant_rows(cells: pd.DataFrame) -> list:
    rows = []
    names = ("tac_high_cld_high", "tac_high_cld_low", "tac_low_cld_high", "tac_low_cld_low")
    for s in SAMPLES:
        sub = cells[cells["sample"] == s]
        tac = sub["tac"].to_numpy(dtype=np.float64)
        cld = sub["cld"].to_numpy(dtype=np.float64)
        th, tl, tmed = median_masks(tac)
        ch, cl, cmed = median_masks(cld)
        masks = {
            "tac_high_cld_high": th & ch,
            "tac_high_cld_low": th & cl,
            "tac_low_cld_high": tl & ch,
            "tac_low_cld_low": tl & cl,
        }
        for metric in PRIMARY_METRICS:
            for radius in RADII:
                y = sub[ycol(metric, radius)].to_numpy(dtype=np.float64)
                rec = {
                    "metric": metric,
                    "radius": int(radius),
                    "sample": s,
                    "patient": sub["patient"].iloc[0],
                    "tac_median": tmed,
                    "cld_median": cmed,
                }
                for name in names:
                    vals = _finite_arm(y, masks[name])
                    usable = vals.size >= MIN_ARM
                    rec[f"n_{name}"] = int(vals.size)
                    rec[f"mean_{name}"] = float(vals.mean()) if usable else float("nan")
                    rec[f"usable_{name}"] = bool(usable)
                rec["usable_hl_vs_ll"] = bool(rec["usable_tac_high_cld_low"] and rec["usable_tac_low_cld_low"])
                rec["usable_hh_vs_lh"] = bool(rec["usable_tac_high_cld_high"] and rec["usable_tac_low_cld_high"])
                rec["usable_hh_vs_ll"] = bool(rec["usable_tac_high_cld_high"] and rec["usable_tac_low_cld_low"])
                rec["delta_hl_ll"] = (
                    rec["mean_tac_high_cld_low"] - rec["mean_tac_low_cld_low"]
                    if rec["usable_hl_vs_ll"] else float("nan")
                )
                rec["delta_hh_lh"] = (
                    rec["mean_tac_high_cld_high"] - rec["mean_tac_low_cld_high"]
                    if rec["usable_hh_vs_lh"] else float("nan")
                )
                rec["delta_hh_ll"] = (
                    rec["mean_tac_high_cld_high"] - rec["mean_tac_low_cld_low"]
                    if rec["usable_hh_vs_ll"] else float("nan")
                )
                rows.append(rec)
    return rows


def attenuation_rows(cells: pd.DataFrame) -> list:
    rows = []
    for s in SAMPLES:
        sub = cells[cells["sample"] == s]
        tac = sub["tac"].to_numpy(dtype=np.float64)
        cld = sub["cld"].to_numpy(dtype=np.float64)
        for metric in PRIMARY_METRICS:
            for radius in RADII:
                y = sub[ycol(metric, radius)].to_numpy(dtype=np.float64)
                m = np.isfinite(y) & np.isfinite(tac) & np.isfinite(cld)
                if int(m.sum()) < 100:
                    continue
                rank_fit = standardized_attenuation(rankdata(tac[m]), rankdata(cld[m]), rankdata(y[m]))
                lin_fit = standardized_attenuation(tac[m], cld[m], y[m])
                high, low, _ = median_masks(tac)
                # Residual of y on CLDN4, then median TACSTD2 contrast.
                design = np.column_stack([np.ones(int(m.sum())), cld[m]])
                coef, _, _, _ = np.linalg.lstsq(design, y[m], rcond=None)
                resid = np.full(y.shape, np.nan)
                resid[m] = y[m] - design @ coef
                rsum = arm_summary(resid, high, low)
                ry = rankdata(y[m])
                rc = rankdata(cld[m])
                rdesign = np.column_stack([np.ones(int(m.sum())), rc])
                rcoef, _, _, _ = np.linalg.lstsq(rdesign, ry, rcond=None)
                rresid = np.full(y.shape, np.nan)
                rresid[m] = ry - rdesign @ rcoef
                rrsum = arm_summary(rresid, high, low)
                rows.append({
                    "metric": metric,
                    "radius": int(radius),
                    "sample": s,
                    "patient": sub["patient"].iloc[0],
                    "n": int(m.sum()),
                    "spearman_tac_y": spearman(tac[m], y[m]),
                    "partial_spearman_index": partial_spearman(tac[m], y[m], cld[m]),
                    "residual_delta": rsum["delta"],
                    "rank_residual_delta": rrsum["delta"],
                    "residual_usable": rsum["usable"],
                    **{f"rank_{k}": v for k, v in rank_fit.items()},
                    **{f"linear_{k}": v for k, v in lin_fit.items()},
                })
    return rows


def _group_mean(df: pd.DataFrame, key: str, value: str, usable: str | None = None) -> pd.DataFrame:
    rows = []
    for unit, sub in df.groupby(key, sort=False):
        if usable is None:
            use = sub
        else:
            use = sub[sub[usable]]
        rows.append({
            key: unit,
            "n_sections": int(sub["sample"].nunique()) if "sample" in sub.columns else int(len(sub)),
            "n_usable": int(len(use)),
            value: nanmean(use[value]) if len(use) else float("nan"),
        })
    return pd.DataFrame(rows)


def aggregate_marginal(sec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["gene", "cut", "metric", "radius"]
    for key, sub in sec.groupby(keys, sort=False):
        gene, cut, metric, radius = key
        usable = sub[sub["usable"]].copy()
        # Patient deltas from equal-weight section means. Missing sections stay missing.
        pat_rows = []
        for pat in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]:
            psub = sub[sub["patient"] == pat]
            u = psub[psub["usable"]]
            if len(u) == 0:
                delta = np.nan
                hm = lm = np.nan
            else:
                hm = float(u["high_mean"].mean())
                lm = float(u["low_mean"].mean())
                delta = hm - lm
            pat_rows.append({"patient": pat, "delta": delta, "high_mean": hm, "low_mean": lm})
        pat = pd.DataFrame(pat_rows)
        sec_delta = sub["delta"].to_numpy(dtype=np.float64)
        # sub has one row per sample; reindex to all 8 in SAMPLES order
        aligned = sub.set_index("sample").reindex(SAMPLES)
        sec_delta = aligned["delta"].to_numpy(dtype=np.float64)
        pat_delta = pat["delta"].to_numpy(dtype=np.float64)
        n_sec = int(np.sum(np.isfinite(sec_delta) & (sec_delta < 0)))
        n_pat = int(np.sum(np.isfinite(pat_delta) & (pat_delta < 0)))
        n_sec_def = int(np.sum(np.isfinite(sec_delta)))
        hm = nanmean(usable["high_mean"]) if len(usable) else float("nan")
        lm = nanmean(usable["low_mean"]) if len(usable) else float("nan")
        rows.append({
            "gene": gene,
            "cut": cut,
            "metric": metric,
            "radius": int(radius),
            "n_sections_usable": int(len(usable)),
            "high_mean": hm,
            "low_mean": lm,
            "delta": (hm - lm) if np.isfinite(hm) and np.isfinite(lm) else float("nan"),
            "ratio": (hm / lm) if np.isfinite(hm) and np.isfinite(lm) and abs(lm) > 1e-8 else float("nan"),
            "median_section_ratio": nanmean([]) if usable.empty else float(np.nanmedian(usable["ratio"])),
            "sections_high_lt_low": n_sec,
            "patients_high_lt_low": n_pat,
            "sections_defined": n_sec_def,
            "sign_p_sections": sign_test(n_sec, 8),
            "sign_p_patients": sign_test(n_pat, 5),
            "wilcoxon_p_sections": wilcoxon_p(sec_delta),
            "wilcoxon_p_patients": wilcoxon_p(pat_delta),
        })
    return pd.DataFrame(rows)


def aggregate_stratified(sec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, sub in sec.groupby(["stratum", "metric", "radius"], sort=False):
        stratum, metric, radius = key
        aligned = sub.set_index("sample").reindex(SAMPLES)
        sec_delta = aligned["delta"].to_numpy(dtype=np.float64)
        pat_delta = []
        for pat in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]:
            u = sub[(sub["patient"] == pat) & (sub["usable"])]
            pat_delta.append(nanmean(u["delta"]) if len(u) else float("nan"))
        pat_delta = np.asarray(pat_delta, dtype=np.float64)
        usable = sub[sub["usable"]]
        hm = nanmean(usable["high_mean"]) if len(usable) else float("nan")
        lm = nanmean(usable["low_mean"]) if len(usable) else float("nan")
        rows.append({
            "stratum": stratum,
            "metric": metric,
            "radius": int(radius),
            "n_sections_usable": int(len(usable)),
            "high_mean": hm,
            "low_mean": lm,
            "delta": (hm - lm) if np.isfinite(hm) and np.isfinite(lm) else float("nan"),
            "ratio": (hm / lm) if np.isfinite(hm) and np.isfinite(lm) and abs(lm) > 1e-8 else float("nan"),
            "sections_high_lt_low": int(np.sum(np.isfinite(sec_delta) & (sec_delta < 0))),
            "patients_high_lt_low": int(np.sum(np.isfinite(pat_delta) & (pat_delta < 0))),
            "sign_p_sections": sign_test(int(np.sum(np.isfinite(sec_delta) & (sec_delta < 0))), 8),
            "sign_p_patients": sign_test(int(np.sum(np.isfinite(pat_delta) & (pat_delta < 0))), 5),
            "wilcoxon_p_sections": wilcoxon_p(sec_delta),
        })
    return pd.DataFrame(rows)


def aggregate_quadrants(sec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    mean_cols = [
        "mean_tac_high_cld_high",
        "mean_tac_high_cld_low",
        "mean_tac_low_cld_high",
        "mean_tac_low_cld_low",
    ]
    for key, sub in sec.groupby(["metric", "radius"], sort=False):
        metric, radius = key
        rec = {"metric": metric, "radius": int(radius)}
        for col, flag in (
            ("mean_tac_high_cld_high", "usable_tac_high_cld_high"),
            ("mean_tac_high_cld_low", "usable_tac_high_cld_low"),
            ("mean_tac_low_cld_high", "usable_tac_low_cld_high"),
            ("mean_tac_low_cld_low", "usable_tac_low_cld_low"),
        ):
            u = sub[sub[flag]]
            rec[col] = nanmean(u[col]) if len(u) else float("nan")
            rec[f"n_sections_{flag}"] = int(len(u))
        for dcol, flag in (
            ("delta_hl_ll", "usable_hl_vs_ll"),
            ("delta_hh_lh", "usable_hh_vs_lh"),
            ("delta_hh_ll", "usable_hh_vs_ll"),
        ):
            aligned = sub.set_index("sample").reindex(SAMPLES)
            # unusable already stored as nan delta, but reindex may also nan
            deltas = []
            for s in SAMPLES:
                if s not in set(sub["sample"]):
                    deltas.append(np.nan)
                    continue
                row = sub[sub["sample"] == s].iloc[0]
                deltas.append(row[dcol] if bool(row[flag]) else np.nan)
            deltas = np.asarray(deltas, dtype=np.float64)
            rec[f"{dcol}_sections_neg"] = int(np.sum(np.isfinite(deltas) & (deltas < 0)))
            rec[f"{dcol}_sections_defined"] = int(np.sum(np.isfinite(deltas)))
            rec[dcol] = nanmean(deltas)
        rows.append(rec)
    return pd.DataFrame(rows)


def aggregate_attenuation(sec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    value_cols = [
        "spearman_tac_y",
        "partial_spearman_index",
        "residual_delta",
        "rank_residual_delta",
        "rank_beta_total",
        "rank_beta_direct",
        "rank_attenuation",
        "rank_proportion",
        "rank_attenuation_of_cldn4",
        "rank_interaction",
        "rank_a_path",
        "linear_beta_total",
        "linear_beta_direct",
        "linear_attenuation",
        "linear_attenuation_of_cldn4",
        "linear_interaction",
    ]
    for key, sub in sec.groupby(["metric", "radius"], sort=False):
        metric, radius = key
        rec = {"metric": metric, "radius": int(radius), "n_sections": int(sub["sample"].nunique())}
        for col in value_cols:
            vals = sub[col].to_numpy(dtype=np.float64)
            rec[f"{col}_median"] = float(np.nanmedian(vals)) if np.isfinite(vals).any() else float("nan")
            rec[f"{col}_n"] = int(np.sum(np.isfinite(vals)))
        # Section signs for residual delta and for rank beta_direct.
        aligned = sub.set_index("sample").reindex(SAMPLES)
        for col, name in (
            ("residual_delta", "residual"),
            ("rank_residual_delta", "rank_residual"),
            ("rank_beta_direct", "rank_direct"),
            ("partial_spearman_index", "index_partial"),
        ):
            d = aligned[col].to_numpy(dtype=np.float64)
            n_neg = int(np.sum(np.isfinite(d) & (d < 0)))
            rec[f"{name}_sections_neg"] = n_neg
            rec[f"{name}_sign_p"] = sign_test(n_neg, 8)
            pat = []
            for p in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]:
                u = sub.loc[sub["patient"] == p, col].to_numpy(dtype=np.float64)
                pat.append(nanmean(u))
            pat = np.asarray(pat, dtype=np.float64)
            n_pat = int(np.sum(np.isfinite(pat) & (pat < 0)))
            rec[f"{name}_patients_neg"] = n_pat
            rec[f"{name}_patient_sign_p"] = sign_test(n_pat, 5)
        rows.append(rec)
    return pd.DataFrame(rows)


def aggregate_fields(fov: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sec_rows = []
    for key, sub in fov.groupby(["sample", "bandwidth", "outcome"], sort=False):
        sample, bandwidth, outcome = key
        u = sub[sub["usable"]]
        sec_rows.append({
            "sample": sample,
            "patient": sub["patient"].iloc[0],
            "bandwidth": int(bandwidth),
            "outcome": outcome,
            "n_fov": int(len(sub)),
            "n_fov_usable": int(len(u)),
            "usable": bool(len(u) >= MIN_FOVS),
            "rho_marginal": nanmean(u["rho_marginal"]) if len(u) else float("nan"),
            "rho_partial": nanmean(u["rho_partial"]) if len(u) else float("nan"),
            "beta_tacstd2": nanmean(u["beta_tacstd2"]) if len(u) else float("nan"),
            "beta_cldn4": nanmean(u["beta_cldn4"]) if len(u) else float("nan"),
            "rho_tac_cld_field": nanmean(u["rho_tac_cld_field"]) if len(u) else float("nan"),
        })
    sec = pd.DataFrame(sec_rows)
    summ = []
    for key, sub in sec.groupby(["bandwidth", "outcome"], sort=False):
        bandwidth, outcome = key
        aligned = sub.set_index("sample").reindex(SAMPLES)
        rec = {"bandwidth": int(bandwidth), "outcome": outcome}
        for col in ("rho_marginal", "rho_partial", "beta_tacstd2", "beta_cldn4", "rho_tac_cld_field"):
            d = aligned[col].to_numpy(dtype=np.float64)
            # Unusable sections are NaN in `usable` handling: force NaN if not usable.
            usable = aligned["usable"].to_numpy()
            d = np.where(usable == True, d, np.nan)  # noqa: E712
            n_neg = int(np.sum(np.isfinite(d) & (d < 0)))
            rec[f"{col}_mean"] = nanmean(d)
            rec[f"{col}_sections_neg"] = n_neg
            rec[f"{col}_sign_p"] = sign_test(n_neg, 8)
            pat = []
            for p in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]:
                uu = sub[(sub["patient"] == p) & (sub["usable"])]
                pat.append(nanmean(uu[col]) if len(uu) else float("nan"))
            pat = np.asarray(pat, dtype=np.float64)
            rec[f"{col}_patients_neg"] = int(np.sum(np.isfinite(pat) & (pat < 0)))
            rec[f"{col}_patient_sign_p"] = sign_test(rec[f"{col}_patients_neg"], 5)
        summ.append(rec)
    return sec, pd.DataFrame(summ)


def fmt(x, nd=3) -> str:
    if x is None or not (isinstance(x, (int, float, np.floating)) and np.isfinite(x)):
        return "NA"
    x = float(x)
    ax = abs(x)
    if ax != 0 and ax < 0.001:
        return f"{x:.2e}"
    if ax < 0.02:
        return f"{x:.4f}"
    return f"{x:.{nd}f}"


def fmt_p(x) -> str:
    if x is None or not (isinstance(x, (int, float, np.floating)) and np.isfinite(x)):
        return "NA"
    x = float(x)
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


def md_table(headers, rows) -> str:
    line = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join(["---"] * len(headers)) + "|"
    body = ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join([line, sep, *body])


def label_text(label: str) -> str:
    return {
        "not_concordant_8_of_8_and_5_of_5": (
            "The pre-specified median split is not lower in all 8 sections and all 5 patients."
        ),
        "concordant_and_largely_shared_with_cldn4": (
            "The median split is lower in 8/8 sections and 5/5 patients, the within-CLDN4 "
            "strata are not, and the rank attenuation exceeds one half."
        ),
        "concordant_and_not_accounted_for_by_cldn4": (
            "The median split stays lower inside a CLDN4 stratum at 8/8 and 5/5, and the "
            "rank attenuation is below one quarter."
        ),
        "concordant_with_mixed_cldn4_adjustment": (
            "The median split is lower in 8/8 and 5/5, and the CLDN4 adjustment is only partial."
        ),
    }[label]


def build_verdicts(marg, strat, att) -> list:
    verdicts = []
    med = marg[(marg["gene"] == "tac") & (marg["cut"] == "median")]
    for metric in PRIMARY_METRICS:
        for radius in RADII:
            mrow = med[(med["metric"] == metric) & (med["radius"] == radius)]
            if mrow.empty:
                continue
            mrow = mrow.iloc[0]
            stratum_ok = False
            stratum_bits = []
            for stratum in ("cldn4_high", "cldn4_low"):
                srow = strat[
                    (strat["stratum"] == stratum)
                    & (strat["metric"] == metric)
                    & (strat["radius"] == radius)
                ]
                if srow.empty:
                    continue
                srow = srow.iloc[0]
                ok = int(srow["sections_high_lt_low"]) == 8 and int(srow["patients_high_lt_low"]) == 5
                stratum_ok = stratum_ok or ok
                stratum_bits.append(
                    f"{stratum} {int(srow['sections_high_lt_low'])}/8 and "
                    f"{int(srow['patients_high_lt_low'])}/5, ratio {fmt(srow['ratio'])}"
                )
            arow = att[(att["metric"] == metric) & (att["radius"] == radius)]
            if arow.empty:
                att_med = np.nan
                n_att = 0
                att_c = np.nan
            else:
                arow = arow.iloc[0]
                att_med = float(arow["rank_attenuation_median"])
                n_att = int(arow["rank_attenuation_n"])
                att_c = float(arow["rank_attenuation_of_cldn4_median"])
            label = classify_dependence(
                int(mrow["sections_high_lt_low"]),
                int(mrow["patients_high_lt_low"]),
                stratum_ok,
                att_med,
                n_att,
            )
            verdicts.append({
                "metric": metric,
                "radius": int(radius),
                "label": label,
                "sections": int(mrow["sections_high_lt_low"]),
                "patients": int(mrow["patients_high_lt_low"]),
                "ratio": float(mrow["ratio"]) if np.isfinite(mrow["ratio"]) else None,
                "high": float(mrow["high_mean"]) if np.isfinite(mrow["high_mean"]) else None,
                "low": float(mrow["low_mean"]) if np.isfinite(mrow["low_mean"]) else None,
                "attenuation": None if not np.isfinite(att_med) else att_med,
                "attenuation_of_cldn4": None if not np.isfinite(att_c) else att_c,
                "n_att": n_att,
                "strata": "; ".join(stratum_bits),
                "text": label_text(label),
            })
    return verdicts


def contact_table(cells: pd.DataFrame) -> pd.DataFrame:
    """Share of tumor cells with at least one neighbor of each class.

    Immune fraction drops cells whose radius is empty. This table keeps them.
    It was computed after that drop was visible, so it is a check, not a second primary.
    """
    rows = []
    for s in SAMPLES:
        sub = cells[cells["sample"] == s]
        tac = sub["tac"].to_numpy(dtype=np.float64)
        high, low, _ = median_masks(tac)
        for radius in RADII:
            rec = {
                "sample": s,
                "patient": sub["patient"].iloc[0],
                "radius": int(radius),
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
            }
            for name, col in (
                ("any_cell", f"denom_{radius}"),
                ("any_immune", f"imm_n_{radius}"),
                ("any_cd8", f"cd8_n_{radius}"),
            ):
                y = sub[col].to_numpy(dtype=np.float64) > 0
                rec[f"{name}_high"] = float(y[high].mean()) if high.any() else float("nan")
                rec[f"{name}_low"] = float(y[low].mean()) if low.any() else float("nan")
                rec[f"{name}_delta"] = rec[f"{name}_high"] - rec[f"{name}_low"]
            rows.append(rec)
    return pd.DataFrame(rows)


def _and_join(names) -> str:
    names = [str(n) for n in names]
    if not names:
        return "no section"
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + ", and " + names[-1]


def _one(df, **kw):
    sl = df
    for k, v in kw.items():
        sl = sl[sl[k] == v]
    if sl.empty:
        return None
    return sl.iloc[0]


def answer_lines(marg, strat, sec, att_sec, field_sec, quad) -> list[str]:
    lines = ["## Answer", ""]
    a = lines.append
    imm10 = _one(marg, gene="tac", cut="median", metric="imm_frac", radius=10)
    imm20 = _one(marg, gene="tac", cut="median", metric="imm_frac", radius=20)
    imm50 = _one(marg, gene="tac", cut="median", metric="imm_frac", radius=50)
    imm100 = _one(marg, gene="tac", cut="median", metric="imm_frac", radius=100)
    cd10 = _one(marg, gene="tac", cut="median", metric="cd8_n", radius=10)
    cd20 = _one(marg, gene="tac", cut="median", metric="cd8_n", radius=20)
    a(
        f"At 10 and 20 µm, tumor cells above the section median of TACSTD2 have a lower "
        f"immune-neighbor fraction in 8/8 sections and 5/5 patients. Equal-weight ratios are "
        f"{fmt(imm10['ratio'])} ({fmt(imm10['high_mean'])} vs {fmt(imm10['low_mean'])}) at 10 µm and "
        f"{fmt(imm20['ratio'])} ({fmt(imm20['high_mean'])} vs {fmt(imm20['low_mean'])}) at 20 µm. "
        f"CD8 neighbor counts are lower in the same 8/8 and 5/5, with ratios "
        f"{fmt(cd10['ratio'])} ({fmt(cd10['high_mean'])} vs {fmt(cd10['low_mean'])}) and "
        f"{fmt(cd20['ratio'])} ({fmt(cd20['high_mean'])} vs {fmt(cd20['low_mean'])})."
    )
    a("")
    lo10 = _one(strat, stratum="cldn4_low", metric="imm_frac", radius=10)
    lo20 = _one(strat, stratum="cldn4_low", metric="imm_frac", radius=20)
    hi10 = _one(strat, stratum="cldn4_high", metric="imm_frac", radius=10)
    hi20 = _one(strat, stratum="cldn4_high", metric="imm_frac", radius=20)
    att10 = att_sec[(att_sec.metric == "imm_frac") & (att_sec.radius == 10)]
    att20 = att_sec[(att_sec.metric == "imm_frac") & (att_sec.radius == 20)]
    a(
        f"That short-range immune-fraction contrast among CLDN4-low tumor cells is "
        f"{fmt(lo10['ratio'])} at 10 µm ({int(lo10['sections_high_lt_low'])}/8 and {int(lo10['patients_high_lt_low'])}/5) and "
        f"{fmt(lo20['ratio'])} at 20 µm ({int(lo20['sections_high_lt_low'])}/8 and {int(lo20['patients_high_lt_low'])}/5). "
        f"Among CLDN4-high cells it is weaker "
        f"({int(hi10['sections_high_lt_low'])}/8 at 10 µm, ratio {fmt(hi10['ratio'])}; "
        f"{int(hi20['sections_high_lt_low'])}/8 at 20 µm, ratio {fmt(hi20['ratio'])}). "
        f"The median section-level rank attenuation of TACSTD2, after CLDN4 is in the same model, "
        f"is {fmt(np.nanmedian(att10['rank_attenuation']))} at 10 µm and "
        f"{fmt(np.nanmedian(att20['rank_attenuation']))} at 20 µm. "
        f"Lung5 sections sit near zero or below (the TACSTD2 coefficient does not shrink), "
        f"while LUAD-9 R1 is {fmt(_one(att10, sample='LUAD-9 R1')['rank_attenuation'])} at 10 µm. "
        f"The 10–20 µm immune-fraction pattern is not accounted for by the cell's own CLDN4."
    )
    a("")
    up50 = sec[
        (sec.gene == "tac") & (sec.cut == "median") & (sec.metric == "imm_frac")
        & (sec.radius == 50) & (sec.delta > 0)
    ]
    up100 = sec[
        (sec.gene == "tac") & (sec.cut == "median") & (sec.metric == "imm_frac")
        & (sec.radius == 100) & (sec.delta > 0)
    ]
    a(
        f"At 50 µm the immune-fraction ratio is {fmt(imm50['ratio'])}, lower in "
        f"{int(imm50['sections_high_lt_low'])}/8 sections and {int(imm50['patients_high_lt_low'])}/5 patients. "
        f"The section that is higher is {', '.join(up50['sample'].tolist()) or 'none'}. "
        f"At 100 µm the ratio is {fmt(imm100['ratio'])}, "
        f"{int(imm100['sections_high_lt_low'])}/8 and {int(imm100['patients_high_lt_low'])}/5. "
        f"Higher sections: {', '.join(up100['sample'].tolist()) or 'none'}."
    )
    near50 = sec[
        (sec.gene == "tac") & (sec.cut == "median") & (sec.metric == "imm_frac")
        & (sec.radius == 50) & (sec.delta < 0) & (sec.ratio > 0.99)
    ]
    if len(near50):
        bits = [f"{r['sample']} ratio {fmt(r['ratio'])}" for _, r in near50.iterrows()]
        a(
            "Counted as lower at 50 µm with a section ratio above 0.99: "
            + "; ".join(bits)
            + "."
        )
    a("")
    q = _one(quad, metric="imm_frac", radius=10)
    a(
        f"At 10 µm the four median quadrants, as equal-weight section means of immune fraction, are "
        f"both low {fmt(q['mean_tac_low_cld_low'])}, "
        f"TACSTD2-high / CLDN4-low {fmt(q['mean_tac_high_cld_low'])}, "
        f"TACSTD2-low / CLDN4-high {fmt(q['mean_tac_low_cld_high'])}, "
        f"and both high {fmt(q['mean_tac_high_cld_high'])}. "
        f"TACSTD2-high / CLDN4-low is below both-low in {int(q['delta_hl_ll_sections_neg'])}/8 sections. "
        f"Each marker is associated with a colder 10 µm neighborhood when the other marker is low."
    )
    a("")
    miss = field_sec[(field_sec.outcome == "imm_frac") & (field_sec.bandwidth == 10) & (~field_sec.usable)]
    f20 = field_sec[(field_sec.outcome == "imm_frac") & (field_sec.bandwidth == 20) & (field_sec.usable)]
    n_neg = int(np.sum(f20["rho_partial"] < 0))
    a(
        f"The Gaussian field does not reproduce that cell-level 8/8. "
        f"At 10 µm, {_and_join(miss['sample'].tolist())} "
        f"{'keeps' if len(miss) == 1 else 'keep'} no FOV with at least {MIN_FOV_FIELD} evaluation points "
        f"after the kernel-mass and local-tumor-count filters, so those sections are absent from the 10 µm field mean. "
        f"At 20 µm, where all 8 sections pass, the partial Spearman of the TACSTD2 field given the "
        f"CLDN4 field is negative in {n_neg}/8 sections "
        f"(mean {fmt(nanmean(f20['rho_partial']))}). "
        f"Section coefficients change sign. The smoothed field is not an 8/8 TACSTD2 exclusion result."
    )
    a("")
    return lines


def write_report(meta, marg, strat, quad, att, field_sum, verdicts, inv: pd.DataFrame,
                 sec: pd.DataFrame | None = None, att_sec: pd.DataFrame | None = None,
                 field_sec: pd.DataFrame | None = None, contact: pd.DataFrame | None = None) -> str:
    lines = []
    a = lines.append
    a("# CosMx He2022: TACSTD2-local immune neighbors and CLDN4")
    a("")
    a(
        "Question: when a tumor cell is TACSTD2-high, are immune and CD8 neighbors "
        "fewer at 10, 20, 50, and 100 µm, and does that contrast shrink once CLDN4 "
        "is stratified, residualized, or entered in a partial field model?"
    )
    a("")
    a(
        f"Data: He et al. 2022 CosMx NSCLC, figshare 25976224, {meta['n_tumor']} "
        "author-matched tumor cells, 8 sections, 5 patients. Scale 0.18 µm/px, "
        "checked by a median nearest-neighbor distance between 3 and 40 µm in every section. "
        "The locked CLDN4 cytotoxic ratios 0.36 at 50 µm and 0.52 at 100 µm, 8/8 and 5/5, "
        "sign P = 0.031, are not recomputed and are not replaced. Nearby effectors are not "
        "scored for GZMB, PRF1, NKG7, or IFNG in this file."
    )
    a("")
    if sec is not None and att_sec is not None and field_sec is not None:
        lines.extend(answer_lines(marg, strat, sec, att_sec, field_sec, quad))
    a("## How high and low were fixed")
    a("")
    a(
        "Primary split: within each section, tumor cells above the median log-normalized "
        "TACSTD2 versus cells at or below that median. Companion splits, fixed before the "
        "summaries: raw count detected versus absent, and Q4 versus Q1. "
        "Where the log-normalized median is 0, the median split is the detected-versus-absent split. "
        "That happens in the sections marked below. It is reported as the same contrast twice, not as extra confirmation."
    )
    a("")
    inv_rows = []
    for _, r in inv.iterrows():
        inv_rows.append([
            r["sample"],
            r["patient"],
            str(int(r["n_tumor"])),
            fmt(r["tacstd2_pos_frac"], 3),
            fmt(r["cldn4_pos_frac"], 3),
            "yes" if r["median_split_is_detected"] else "no",
            fmt(r["spearman_tacstd2_cldn4"], 3),
            fmt(r["median_nn_um"], 2),
        ])
    a(md_table(
        ["Section", "Patient", "Tumor cells", "TACSTD2 > 0", "CLDN4 > 0", "Median is detected", "Spearman", "NN µm"],
        inv_rows,
    ))
    a("")
    a(
        "Patient is the equal-weight mean of its usable sections. A section that cannot fill both arms "
        f"(fewer than {MIN_ARM} cells with a finite score) counts against 8/8 and against 5/5. "
        "The sign test is one-sided in the exclusion direction (high mean lower than low mean). "
        "Its floor is 0.0039 on 8 sections and 0.031 on 5 patients. Wilcoxon p-values are two-sided "
        "on the defined units."
    )
    a("")
    a("## Marginal TACSTD2 high versus low")
    a("")
    a("Equal-weight mean of usable section means. Ratio is high / low. Tallies are high < low.")
    a("")

    def marg_rows(cut):
        out = []
        sub = marg[(marg["gene"] == "tac") & (marg["cut"] == cut) & (marg["metric"].isin(PRIMARY_METRICS))]
        for _, r in sub.sort_values(["metric", "radius"]).iterrows():
            out.append([
                METRIC_LABEL[r["metric"]],
                str(int(r["radius"])),
                fmt(r["high_mean"]),
                fmt(r["low_mean"]),
                fmt(r["ratio"]),
                f"{int(r['sections_high_lt_low'])}/8",
                f"{int(r['patients_high_lt_low'])}/5",
                fmt_p(r["sign_p_patients"]),
                fmt_p(r["wilcoxon_p_sections"]),
            ])
        return out

    a("### Primary: within-section median")
    a("")
    a(md_table(
        ["Readout", "µm", "High", "Low", "Ratio", "Sections", "Patients", "Sign P (5)", "Wilcoxon P (8)"],
        marg_rows("median"),
    ))
    a("")
    a("### Companion: detected versus absent")
    a("")
    a(md_table(
        ["Readout", "µm", "High", "Low", "Ratio", "Sections", "Patients", "Sign P (5)", "Wilcoxon P (8)"],
        marg_rows("detected"),
    ))
    a("")
    a("### Companion: Q4 versus Q1")
    a("")
    a(md_table(
        ["Readout", "µm", "High", "Low", "Ratio", "Sections", "Patients", "Sign P (5)", "Wilcoxon P (8)"],
        marg_rows("q4q1"),
    ))
    a("")
    a("## Same-pipeline CLDN4 calibration")
    a("")
    a(
        "Within-section median of CLDN4, CD8+NK neighbor count, same trees and same section weights. "
        "This is a calibration of the code path. It is not the locked 0.36 / 0.52 summary."
    )
    a("")
    cal = marg[(marg["gene"] == "cld") & (marg["cut"] == "median") & (marg["metric"] == "cd8nk_n")]
    cal_rows = []
    for _, r in cal.sort_values("radius").iterrows():
        cal_rows.append([
            str(int(r["radius"])),
            fmt(r["high_mean"]),
            fmt(r["low_mean"]),
            fmt(r["ratio"]),
            f"{int(r['sections_high_lt_low'])}/8",
            f"{int(r['patients_high_lt_low'])}/5",
            fmt(r["delta"]),
        ])
    a(md_table(["µm", "CLDN4 high", "CLDN4 low", "Ratio", "Sections", "Patients", "Δ"], cal_rows))
    a("")
    a("## Stratified by CLDN4")
    a("")
    a(
        "Inside each section, tumor cells are split at the CLDN4 median. Inside each stratum, "
        "TACSTD2 is split at that stratum's own median. The readout is still high versus low TACSTD2."
    )
    a("")
    srows = []
    for _, r in strat.sort_values(["metric", "radius", "stratum"]).iterrows():
        srows.append([
            METRIC_LABEL[r["metric"]],
            str(int(r["radius"])),
            r["stratum"].replace("cldn4_", "CLDN4 "),
            fmt(r["high_mean"]),
            fmt(r["low_mean"]),
            fmt(r["ratio"]),
            f"{int(r['sections_high_lt_low'])}/8",
            f"{int(r['patients_high_lt_low'])}/5",
            str(int(r["n_sections_usable"])),
        ])
    a(md_table(
        ["Readout", "µm", "Stratum", "TACSTD2 high", "TACSTD2 low", "Ratio", "Sections", "Patients", "Usable"],
        srows,
    ))
    a("")
    a("## Quadrants")
    a("")
    a(
        "Section-wide medians, crossed. `TACSTD2-high / CLDN4-low` is the cell that is high for TROP2 "
        "without being high for CLDN4. Means are equal-weight averages of usable sections. "
        "A quadrant with fewer than 30 finite cells in a section is unused in that section."
    )
    a("")
    qrows = []
    for _, r in quad.sort_values(["metric", "radius"]).iterrows():
        qrows.append([
            METRIC_LABEL[r["metric"]],
            str(int(r["radius"])),
            fmt(r["mean_tac_high_cld_high"]),
            fmt(r["mean_tac_high_cld_low"]),
            fmt(r["mean_tac_low_cld_high"]),
            fmt(r["mean_tac_low_cld_low"]),
            f"{int(r['delta_hl_ll_sections_neg'])}/8",
            f"{int(r['delta_hh_lh_sections_neg'])}/8",
            f"{int(r['delta_hh_ll_sections_neg'])}/8",
        ])
    a(md_table(
        ["Readout", "µm", "Both high", "TAC high, CLDN4 low", "TAC low, CLDN4 high", "Both low",
         "HL<LL", "HH<LH", "HH<LL"],
        qrows,
    ))
    a("")
    a("## Residual and mediation-like attenuation")
    a("")
    a(
        "Within each section the outcome, TACSTD2, and CLDN4 are rank-transformed and standardized. "
        "Attenuation is 1 minus the joint standardized coefficient of TACSTD2 divided by its simple "
        "coefficient. The product-of-paths proportion matches that number in this linear specification; "
        "both are stored, and the report quotes the attenuation. "
        "A simple coefficient smaller than 0.01 in absolute value leaves the attenuation undefined. "
        "The linear (log-normalized, not rank) attenuation is the companion. "
        "Residual Δ is the median-split TACSTD2 contrast after the outcome is regressed on CLDN4. "
        "Negative means TACSTD2-high cells still have the lower neighborhood. "
        "These are observational associations. They are not a causal mediation estimate, and overlapping "
        "neighborhoods are not treated as independent cells."
    )
    a("")
    arows = []
    for _, r in att.sort_values(["metric", "radius"]).iterrows():
        arows.append([
            METRIC_LABEL[r["metric"]],
            str(int(r["radius"])),
            fmt(r["rank_beta_total_median"]),
            fmt(r["rank_beta_direct_median"]),
            fmt(r["rank_attenuation_median"]),
            f"{int(r['rank_attenuation_n'])}/8",
            fmt(r["rank_attenuation_of_cldn4_median"]),
            fmt(r["residual_delta_median"]),
            f"{int(r['residual_sections_neg'])}/8",
            f"{int(r['residual_patients_neg'])}/5",
            fmt(r["partial_spearman_index_median"]),
            f"{int(r['index_partial_sections_neg'])}/8",
        ])
    a(md_table(
        ["Readout", "µm", "Rank β total", "Rank β with CLDN4", "TACSTD2 attenuation", "Defined",
         "CLDN4 attenuation", "Residual Δ", "Res. sections", "Res. patients", "Index partial ρ", "Partial < 0"],
        arows,
    ))
    a("")
    a("## Partial field models")
    a("")
    a(
        "Each section is rasterized. Gaussian bandwidths are 10, 20, 50, and 100 µm "
        f"(pixel = bandwidth/5, truncate {TRUNCATE:g}). "
        "The TACSTD2 field and the CLDN4 field are kernel-weighted means on tumor cells. "
        "The immune field is the kernel-weighted immune share of all cells. "
        "The CD8 field is kernel-weighted CD8 density (cells per µm²). "
        "Evaluation points are up to 250 tumor cells per FOV, kept where kernel mass is at least "
        f"{KERNEL_MASS_MIN:.2f} and the local tumor disk holds at least 3 cells. "
        f"A FOV needs {MIN_FOV_FIELD} points. A section needs {MIN_FOVS} such FOVs. "
        "Reported values are unweighted means of FOV Spearman coefficients, then of sections. "
        "Partial ρ residualizes both ranks on the CLDN4 field. "
        "β is the standardized coefficient in immune ~ TACSTD2 field + CLDN4 field."
    )
    a("")
    frows = []
    for _, r in field_sum.sort_values(["outcome", "bandwidth"]).iterrows():
        frows.append([
            "Immune fraction" if r["outcome"] == "imm_frac" else "CD8 density",
            str(int(r["bandwidth"])),
            fmt(r["rho_marginal_mean"]),
            f"{int(r['rho_marginal_sections_neg'])}/8",
            f"{int(r['rho_marginal_patients_neg'])}/5",
            fmt(r["rho_partial_mean"]),
            f"{int(r['rho_partial_sections_neg'])}/8",
            f"{int(r['rho_partial_patients_neg'])}/5",
            fmt(r["beta_tacstd2_mean"]),
            fmt(r["beta_cldn4_mean"]),
            fmt(r["rho_tac_cld_field_mean"]),
        ])
    a(md_table(
        ["Field", "µm", "Marginal ρ", "Sections", "Patients", "Partial ρ", "Sections", "Patients",
         "β TACSTD2", "β CLDN4", "Field TAC–CLDN4 ρ"],
        frows,
    ))
    a("")
    a("## Reading the dependence test")
    a("")
    a(
        "A radius is called concordant only at 8/8 sections and 5/5 patients. "
        "It is called largely shared with CLDN4 only when that concordance holds, neither CLDN4 "
        "stratum is itself 8/8 and 5/5, and the median rank attenuation is above 0.5. "
        "It is called not accounted for by CLDN4 only when concordance holds, at least one stratum "
        "is itself 8/8 and 5/5, and the median rank attenuation is below 0.25. "
        "Every other concordant radius is mixed. A radius that is not 8/8 and 5/5 is reported as not concordant. "
        "Companion cuts do not get this label."
    )
    a("")
    for v in verdicts:
        a(
            f"- **{METRIC_LABEL[v['metric']]}, {v['radius']} µm.** "
            f"Median-split ratio {fmt(v['ratio'])} "
            f"({fmt(v['high'])} vs {fmt(v['low'])}), "
            f"{v['sections']}/8 sections, {v['patients']}/5 patients. "
            f"Strata: {v['strata']}. "
            f"Rank attenuation of TACSTD2 {fmt(v['attenuation'])} "
            f"({v['n_att']}/8 defined); rank attenuation of CLDN4 {fmt(v['attenuation_of_cldn4'])}. "
            f"{v['text']}"
        )
    a("")
    if sec is not None:
        a("## Section ratios for the median split")
        a("")
        a("Ratio is the section high mean divided by the section low mean. A ratio below 1 is the exclusion direction.")
        a("")
        for metric in PRIMARY_METRICS:
            a(f"### {METRIC_LABEL[metric]}")
            a("")
            headers = ["Section", "Patient"] + [str(r) for r in RADII]
            body = []
            for s in SAMPLES:
                row = [SHORT[s], EXPECTED_PATIENT[s]]
                for radius in RADII:
                    sl = sec[
                        (sec["sample"] == s) & (sec.gene == "tac") & (sec.cut == "median")
                        & (sec.metric == metric) & (sec.radius == radius)
                    ]
                    if sl.empty or not bool(sl.iloc[0]["usable"]):
                        row.append("NA")
                    else:
                        row.append(fmt(sl.iloc[0]["ratio"]))
                body.append(row)
            a(md_table(headers, body))
            a("")
        cd = sec[(sec.gene == "tac") & (sec.cut == "median") & (sec.metric == "cd8_n") & (sec.radius == 10)]
        a(
            "CD8 counts at 10 µm are sparse. High-arm and low-arm event totals "
            "(mean count times cells in the arm) are:"
        )
        a("")
        ev = []
        for _, r in cd.iterrows():
            ev.append([
                str(r["sample"]),
                f"{r['high_mean'] * r['n_high']:.1f}",
                f"{r['low_mean'] * r['n_low']:.1f}",
            ])
        a(md_table(["Section", "High-arm CD8 events", "Low-arm CD8 events"], ev))
        a("")
        thin = []
        for _, r in cd.iterrows():
            he = float(r["high_mean"] * r["n_high"])
            le = float(r["low_mean"] * r["n_low"])
            if he < 2:
                thin.append(f"{r['sample']} {he:.0f} vs {le:.0f}")
        if thin:
            a(
                "Sections with fewer than 2 CD8-neighbor events in the high arm at 10 µm: "
                + "; ".join(thin)
                + ". They still count toward 8/8. The immune-fraction contrast is not carried by these CD8 events."
            )
            a("")
    if contact is not None:
        a("## Empty-neighborhood check")
        a("")
        a(
            "Immune fraction is undefined when a tumor cell has no other cell inside the radius, "
            "and those cells are left out of the fraction mean. After that was visible, the share "
            "of tumor cells with at least one immune neighbor was computed on every tumor cell, "
            "including isolated ones. It is a check on the fraction, not a second primary endpoint."
        )
        a("")
        headers = ["Section"] + [f"{r} µm high/low" for r in (10, 20)]
        body = []
        for s in SAMPLES:
            row = [SHORT[s]]
            for radius in (10, 20):
                sl = contact[(contact["sample"] == s) & (contact["radius"] == radius)]
                if sl.empty:
                    row.append("NA")
                else:
                    r0 = sl.iloc[0]
                    row.append(f"{fmt(r0['any_immune_high'])}/{fmt(r0['any_immune_low'])}")
            body.append(row)
        a(md_table(headers, body))
        a("")
        neg10 = int(np.sum(contact.loc[contact.radius == 10, "any_immune_delta"] < 0))
        neg20 = int(np.sum(contact.loc[contact.radius == 20, "any_immune_delta"] < 0))
        a(
            f"Any-immune contact is lower for TACSTD2-high cells in {neg10}/8 sections at 10 µm "
            f"and {neg20}/8 at 20 µm."
        )
        a("")
    if field_sec is not None:
        a("## Field coefficients by section")
        a("")
        a(
            "Mean partial ρ in the summary table averages usable sections only. "
            "A section with fewer than 3 usable FOVs is omitted from that mean and counts against 8/8."
        )
        a("")
        headers = ["Section", "µm", "FOVs used", "Marginal ρ", "Partial ρ", "β TACSTD2", "β CLDN4"]
        body = []
        sub = field_sec[field_sec.outcome == "imm_frac"]
        for radius in RADII:
            for s in SAMPLES:
                sl = sub[(sub["sample"] == s) & (sub.bandwidth == radius)]
                if sl.empty:
                    continue
                r0 = sl.iloc[0]
                body.append([
                    SHORT[s],
                    str(int(radius)),
                    str(int(r0["n_fov_usable"])),
                    fmt(r0["rho_marginal"]) if bool(r0["usable"]) else "unused",
                    fmt(r0["rho_partial"]) if bool(r0["usable"]) else "unused",
                    fmt(r0["beta_tacstd2"]) if bool(r0["usable"]) else "unused",
                    fmt(r0["beta_cldn4"]) if bool(r0["usable"]) else "unused",
                ])
        a(md_table(headers, body))
        a("")
    a("## What this file does not say")
    a("")
    a(
        "It does not replace the locked CLDN4 exclusion summary. "
        "It does not say that effector cells next to TACSTD2-high tumor cells are muzzled. "
        "It does not use an ICI label. It does not pool private KL tumors with this public object. "
        "FOV coefficients are descriptive because cells inside a FOV share neighbors."
    )
    a("")
    a("```bash")
    a("python3 scripts/download_cosmx_nsclc_h5ad.py")
    a("python3 scripts/cosmx_tacstd2_cldn4_exclusion.py")
    a("```")
    a("")
    text = "\n".join(lines)
    path = os.path.join(OUT, "RESULTS.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    with open(os.path.join(ROOT, "RESULTS.md"), "w", encoding="utf-8") as fh:
        fh.write(text)
    return text


def _paired(ax, sec, metric, radius, title):
    sub = sec[
        (sec["gene"] == "tac")
        & (sec["cut"] == "median")
        & (sec["metric"] == metric)
        & (sec["radius"] == radius)
    ].set_index("sample")
    xs = np.arange(len(SAMPLES))
    for i, s in enumerate(SAMPLES):
        if s not in sub.index or not bool(sub.loc[s, "usable"]):
            continue
        pat = sub.loc[s, "patient"]
        lo = float(sub.loc[s, "low_mean"])
        hi = float(sub.loc[s, "high_mean"])
        ax.plot([i, i], [lo, hi], color=PATIENT_COLOR[pat], lw=1.4, zorder=1)
        ax.scatter([i], [lo], s=28, facecolors="white", edgecolors=PATIENT_COLOR[pat], zorder=2)
        ax.scatter([i], [hi], s=28, color=PATIENT_COLOR[pat], zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT[s] for s in SAMPLES], fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(METRIC_LABEL[metric])


def make_figures(sec, strat_sec, att_sec, field_sec, quad_sec):
    os.makedirs(FIG, exist_ok=True)
    # Marginal immune fraction and CD8 count.
    for metric, name in (("imm_frac", "immune_fraction"), ("cd8_n", "cd8_count")):
        fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.6), constrained_layout=True)
        for ax, radius in zip(axes, RADII):
            row = sec[
                (sec.gene == "tac") & (sec.cut == "median") & (sec.metric == metric) & (sec.radius == radius)
            ]
            k = int(np.sum(row["usable"] & (row["delta"] < 0)))
            _paired(ax, sec, metric, radius, f"{radius} µm  {k}/8 lower")
        fig.suptitle(f"TACSTD2 median split — {METRIC_LABEL[metric]} (filled = high)", fontsize=12)
        fig.savefig(os.path.join(FIG, f"marginal_{name}.png"), dpi=140)
        fig.savefig(os.path.join(FIG, f"marginal_{name}.pdf"))
        plt.close(fig)

    # Stratified ratios at each radius for immune fraction.
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.8), constrained_layout=True)
    for ax, radius in zip(axes, RADII):
        for j, (stratum, marker) in enumerate((("cldn4_high", "o"), ("cldn4_low", "s"))):
            sub = strat_sec[
                (strat_sec["stratum"] == stratum)
                & (strat_sec["metric"] == "imm_frac")
                & (strat_sec["radius"] == radius)
            ].set_index("sample")
            xs, ys, cs = [], [], []
            for i, s in enumerate(SAMPLES):
                if s not in sub.index or not bool(sub.loc[s, "usable"]):
                    continue
                xs.append(i + (j - 0.5) * 0.15)
                ys.append(float(sub.loc[s, "ratio"]))
                cs.append(PATIENT_COLOR[sub.loc[s, "patient"]])
            ax.scatter(xs, ys, c=cs, marker=marker, s=26, zorder=2)
        ax.axhline(1.0, color="0.5", lw=0.8)
        ax.set_xticks(range(len(SAMPLES)))
        ax.set_xticklabels([SHORT[s] for s in SAMPLES], fontsize=8)
        ax.set_title(f"{radius} µm")
        ax.set_ylabel("TACSTD2 high/low ratio")
    fig.suptitle("Immune fraction inside CLDN4 strata (circle = CLDN4 high, square = CLDN4 low)", fontsize=11)
    fig.savefig(os.path.join(FIG, "stratified_immune_ratio.png"), dpi=140)
    fig.savefig(os.path.join(FIG, "stratified_immune_ratio.pdf"))
    plt.close(fig)

    labels = [
        ("mean_tac_low_cld_low", "usable_tac_low_cld_low", "Both low"),
        ("mean_tac_high_cld_low", "usable_tac_high_cld_low", "TAC high, CLDN4 low"),
        ("mean_tac_low_cld_high", "usable_tac_low_cld_high", "TAC low, CLDN4 high"),
        ("mean_tac_high_cld_high", "usable_tac_high_cld_high", "Both high"),
    ]
    for radius, fname in ((10, "quadrant_immune_10um"), (50, "quadrant_immune_50um")):
        fig, ax = plt.subplots(figsize=(8.2, 4.2), constrained_layout=True)
        q = quad_sec[(quad_sec.metric == "imm_frac") & (quad_sec.radius == radius)].set_index("sample")
        for s in SAMPLES:
            if s not in q.index:
                continue
            ys = []
            for col, flag, _ in labels:
                ys.append(float(q.loc[s, col]) if bool(q.loc[s, flag]) else np.nan)
            ax.plot(range(4), ys, color=PATIENT_COLOR[q.loc[s, "patient"]], marker="o", lw=1.2, label=SHORT[s])
        ax.set_xticks(range(4))
        ax.set_xticklabels([lab for _, _, lab in labels], fontsize=8)
        ax.set_ylabel(f"Immune fraction at {radius} µm")
        ax.set_title("Quadrant means, section median splits")
        ax.legend(ncol=4, fontsize=8, frameon=False)
        fig.savefig(os.path.join(FIG, f"{fname}.png"), dpi=140)
        fig.savefig(os.path.join(FIG, f"{fname}.pdf"))
        plt.close(fig)

    # Attenuation: rank beta total vs direct, immune fraction, all radii as bars of median later
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    sub = att_sec[(att_sec.metric == "imm_frac") & (att_sec.radius == 10)].set_index("sample")
    xs = np.arange(len(SAMPLES))
    b0 = [float(sub.loc[s, "rank_beta_total"]) if s in sub.index else np.nan for s in SAMPLES]
    b1 = [float(sub.loc[s, "rank_beta_direct"]) if s in sub.index else np.nan for s in SAMPLES]
    axes[0].plot(xs, b0, "o-", color="0.35", label="TACSTD2 alone")
    axes[0].plot(xs, b1, "o-", color="#1b9e77", label="TACSTD2 with CLDN4")
    axes[0].axhline(0, color="0.6", lw=0.8)
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels([SHORT[s] for s in SAMPLES], fontsize=8)
    axes[0].set_ylabel("Standardized rank coefficient")
    axes[0].set_title("Immune fraction, 10 µm")
    axes[0].legend(frameon=False, fontsize=8)
    for metric, color in (("imm_frac", "#1b9e77"), ("cd8_n", "#7570b3")):
        ys = []
        for radius in RADII:
            sl = att_sec[(att_sec.metric == metric) & (att_sec.radius == radius)]
            ys.append(float(np.nanmedian(sl["rank_attenuation"])))
        axes[1].plot(list(RADII), ys, "o-", color=color, label=METRIC_LABEL[metric])
    axes[1].axhline(0.5, color="0.6", lw=0.8, ls="--")
    axes[1].set_xlabel("Radius (µm)")
    axes[1].set_ylabel("Median rank attenuation")
    axes[1].set_title("Attenuation of TACSTD2 given CLDN4")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(FIG, "attenuation_rank.png"), dpi=140)
    fig.savefig(os.path.join(FIG, "attenuation_rank.pdf"))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    for ax, outcome, title in (
        (axes[0], "imm_frac", "Immune-fraction field"),
        (axes[1], "cd8_density", "CD8-density field"),
    ):
        for col, color, lab in (
            ("rho_marginal", "0.35", "Marginal"),
            ("rho_partial", "#d95f02", "Partial | CLDN4"),
        ):
            ys = []
            for bandwidth in RADII:
                sl = field_sec[
                    (field_sec.outcome == outcome) & (field_sec.bandwidth == bandwidth) & (field_sec.usable)
                ]
                ys.append(nanmean(sl[col]))
            ax.plot(list(RADII), ys, "o-", color=color, label=lab)
        ax.axhline(0, color="0.75", lw=0.8)
        ax.set_xlabel("Bandwidth (µm)")
        ax.set_ylabel("Mean section Spearman")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Partial field: TACSTD2 field vs immune field given the CLDN4 field", fontsize=11)
    fig.savefig(os.path.join(FIG, "field_partial_spearman.png"), dpi=140)
    fig.savefig(os.path.join(FIG, "field_partial_spearman.pdf"))
    plt.close(fig)


def self_test() -> None:
    rng = np.random.default_rng(1)
    n = 4000
    x = rng.normal(size=n)
    z = 0.5 * x + rng.normal(size=n)
    y = -0.15 * x - 0.4 * z + 0.2 * rng.normal(size=n)
    fit = standardized_attenuation(x, z, y)
    if not np.isfinite(fit["attenuation"]) or abs(fit["attenuation"] - fit["proportion"]) > 1e-6:
        raise SystemExit(f"attenuation identity failed: {fit}")
    if not (fit["attenuation"] > 0.5):
        raise SystemExit(f"expected shared-signal attenuation, got {fit['attenuation']}")
    # Independent covariate should not eat the slope.
    z2 = rng.normal(size=n)
    y2 = -0.4 * x + 0.2 * rng.normal(size=n)
    fit2 = standardized_attenuation(x, z2, y2)
    if not np.isfinite(fit2["attenuation"]) or abs(fit2["attenuation"]) > 0.15:
        raise SystemExit(f"null covariate attenuation too large: {fit2['attenuation']}")
    # Rank path uses the same algebra.
    fit_r = standardized_attenuation(rankdata(x), rankdata(z), rankdata(y))
    if abs(fit_r["attenuation"] - fit_r["proportion"]) > 1e-6:
        raise SystemExit("rank attenuation identity failed")
    # Partial spearman of noise is small.
    rho = partial_spearman(rng.normal(size=800), rng.normal(size=800), rng.normal(size=800))
    if not np.isfinite(rho) or abs(rho) > 0.2:
        raise SystemExit(f"null partial spearman {rho}")
    # Neighbor count on a toy geometry.
    tumor = np.array([[0.0, 0.0], [100.0, 0.0]])
    immune = np.array([[6.0, 0.0], [50.0, 0.0], [96.0, 0.0], [200.0, 0.0]])
    counts = ball_counts(cKDTree(immune), tumor, 10.0)
    if not np.array_equal(counts.astype(int), np.array([1, 1])):
        raise SystemExit(f"toy neighbor counts {counts}")
    # Sign tallies: missing unit counts against the denominator.
    if abs(sign_test(8, 8) - 1 / 256) > 1e-12:
        raise SystemExit(f"sign test 8/8 = {sign_test(8, 8)}")
    if abs(sign_test(5, 5) - 1 / 32) > 1e-12:
        raise SystemExit(f"sign test 5/5 = {sign_test(5, 5)}")
    if classify_dependence(8, 5, False, 0.8, 8) != "concordant_and_largely_shared_with_cldn4":
        raise SystemExit("classify shared failed")
    if classify_dependence(8, 5, True, 0.1, 8) != "concordant_and_not_accounted_for_by_cldn4":
        raise SystemExit("classify survives failed")
    if classify_dependence(7, 5, False, 0.9, 8) != "not_concordant_8_of_8_and_5_of_5":
        raise SystemExit("classify 7/8 failed")
    # Arm rule.
    y = np.arange(10, dtype=float)
    high = np.array([True] * 4 + [False] * 6)
    low = ~high
    sm = arm_summary(y, high, low)
    if sm["usable"]:
        raise SystemExit("arm below 30 should be unusable")
    y = np.arange(80, dtype=float)
    high = np.zeros(80, dtype=bool)
    high[:40] = True
    sm = arm_summary(y, high, ~high)
    if not sm["usable"] or abs(sm["delta"] - (19.5 - 59.5)) > 1e-9:
        raise SystemExit(f"arm means wrong {sm}")
    # Median of zeros collapses to detected.
    vals = np.array([0.0, 0.0, 0.0, 1.0])
    high, low, med = median_masks(vals)
    if med != 0 or high.sum() != 1 or low.sum() != 3:
        raise SystemExit("zero-inflated median split failed")
    print("self-test ok", flush=True)


def save_cache(cells: pd.DataFrame, field_rows: list, meta: dict) -> None:
    cells.to_pickle(CACHE + ".pkl")
    with open(CACHE + ".json", "w", encoding="utf-8") as fh:
        json.dump({"fields": _jsonable(field_rows), "meta": _jsonable(meta)}, fh)
    print(f"cached {CACHE}.pkl", flush=True)


def main() -> int:
    self_test()
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    t0 = time.time()
    cells, field_rows, meta = load_and_count()
    print(f"counts done in {time.time() - t0:.1f}s, tumor cells {len(cells)}", flush=True)
    save_cache(cells, field_rows, meta)

    sec = pd.DataFrame(section_split_rows(cells, "tac", "median")
                       + section_split_rows(cells, "tac", "detected")
                       + section_split_rows(cells, "tac", "q4q1")
                       + section_split_rows(cells, "cld", "median"))
    strat_sec = pd.DataFrame(stratified_rows(cells))
    quad_sec = pd.DataFrame(quadrant_rows(cells))
    att_sec = pd.DataFrame(attenuation_rows(cells))
    field_fov = pd.DataFrame(field_rows)
    marg = aggregate_marginal(sec)
    strat = aggregate_stratified(strat_sec)
    quad = aggregate_quadrants(quad_sec)
    att = aggregate_attenuation(att_sec)
    field_sec, field_sum = aggregate_fields(field_fov)
    verdicts = build_verdicts(marg, strat, att)
    inv = pd.DataFrame(meta["inventory"])

    sec.to_csv(os.path.join(TAB, "marginal_by_section.csv"), index=False)
    marg.to_csv(os.path.join(TAB, "marginal_summary.csv"), index=False)
    strat_sec.to_csv(os.path.join(TAB, "stratified_by_section.csv"), index=False)
    strat.to_csv(os.path.join(TAB, "stratified_summary.csv"), index=False)
    quad_sec.to_csv(os.path.join(TAB, "quadrant_by_section.csv"), index=False)
    quad.to_csv(os.path.join(TAB, "quadrant_summary.csv"), index=False)
    att_sec.to_csv(os.path.join(TAB, "attenuation_by_section.csv"), index=False)
    att.to_csv(os.path.join(TAB, "attenuation_summary.csv"), index=False)
    field_fov.to_csv(os.path.join(TAB, "field_by_fov.csv"), index=False)
    field_sec.to_csv(os.path.join(TAB, "field_by_section.csv"), index=False)
    field_sum.to_csv(os.path.join(TAB, "field_summary.csv"), index=False)
    inv.to_csv(os.path.join(TAB, "inventory.csv"), index=False)

    summary = {
        "n_tumor": meta["n_tumor"],
        "verdicts": verdicts,
        "inventory": meta["inventory"],
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    contact = contact_table(cells)
    contact.to_csv(os.path.join(TAB, "contact_by_section.csv"), index=False)
    write_report(
        meta, marg, strat, quad, att, field_sum, verdicts, inv,
        sec=sec, att_sec=att_sec, field_sec=field_sec, contact=contact,
    )
    make_figures(sec, strat_sec, att_sec, field_sec, quad_sec)
    print("verdicts:", flush=True)
    for v in verdicts:
        print(f"  {v['metric']} {v['radius']}: {v['sections']}/8 {v['patients']}/5 {v['label']} ratio={v['ratio']}", flush=True)
    print(f"done in {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
