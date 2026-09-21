#!/usr/bin/env python3
"""CosMx-native cell-type density fields: regress tumor CLDN4 on predicted CD8.

He et al. 2022 CosMx SMI NSCLC (8 sections / 5 patients) has no matched Visium
companion, so Tangram / RCTD / cell2location are not run. Author cell-type
points are mapped onto edge-corrected Gaussian density fields. The primary
regression is within-section Spearman of tumor-cell CLDN4 versus the CD8
density at that cell.

Pre-specified before looking at the association:
  - Primary predictor: author CD8 intensity (T CD8 naive + T CD8 memory),
    cells per (100 µm)^2, bandwidth 50 µm, kernel-mass mask >= 0.80.
  - Primary statistic: unweighted mean of the 8 within-section Spearman ρ.
  - Companion contrast: CLDN4-high (count >= section 75th percentile and >= 1)
    versus CLDN4-undetected (count == 0); mean of 8 median differences.
  - Secondary: CD8 fraction of the leave-one-out density field (RCTD-style
    proportion), bandwidth 100 µm, CD8+NK intensity, partial Spearman given
    local tumor density, and the same contrasts at mask >= 0.20.
  - Inference: 999 within-FOV permutations of CLDN4. Cell-level Spearman
    p-values are not used (spatial dependence).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import binomtest, norm, rankdata, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx_nsclc_cells.csv.gz"
GENES = ROOT / "data" / "cosmx_960_genes.csv"
OUT = ROOT / "results" / "cosmx_density_field"
FIG = OUT / "figures"
TAB = OUT / "tables"

UM_PER_MM = 1000.0
PIXEL_UM = 2.0
BANDWIDTHS = (50.0, 100.0)
PRIMARY_H = 50.0
PRIMARY_MASK = 0.80
SENS_MASK = 0.20
TRUNCATE = 4.0
N_PERM = 999
SEED = 20260921
PER_100 = 100.0 ** 2  # µm² → cells per (100 µm)²

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]

SUPER_OF = {}
for _name in ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13"):
    SUPER_OF[_name] = "tumor"
SUPER_OF.update(
    {
        "epithelial": "epithelial",
        "fibroblast": "fibroblast",
        "endothelial": "endothelial",
        "T CD8 naive": "CD8",
        "T CD8 memory": "CD8",
        "T CD4 naive": "CD4",
        "T CD4 memory": "CD4",
        "Treg": "Treg",
        "NK": "NK",
        "B-cell": "B",
        "plasmablast": "plasmablast",
        "macrophage": "macrophage",
        "monocyte": "monocyte",
        "neutrophil": "neutrophil",
        "mDC": "mDC",
        "pDC": "pDC",
        "mast": "mast",
    }
)
SUPER_ORDER = [
    "CD8",
    "CD4",
    "Treg",
    "NK",
    "macrophage",
    "monocyte",
    "neutrophil",
    "mDC",
    "pDC",
    "B",
    "plasmablast",
    "mast",
    "fibroblast",
    "endothelial",
    "epithelial",
    "tumor",
]


def _grid_spec(window_xy: np.ndarray, pixel: float):
    xmin = float(window_xy[:, 0].min())
    ymin = float(window_xy[:, 1].min())
    xmax = float(window_xy[:, 0].max())
    ymax = float(window_xy[:, 1].max())
    nx = max(1, int(np.ceil(max(xmax - xmin, pixel) / pixel)))
    ny = max(1, int(np.ceil(max(ymax - ymin, pixel) / pixel)))
    return xmin, ymin, xmin + nx * pixel, ymin + ny * pixel, nx, ny


def _kernel_center(sigma: float, truncate: float = TRUNCATE) -> float:
    rad = int(np.ceil(truncate * sigma)) + 1
    g = np.zeros((2 * rad + 1, 2 * rad + 1), dtype=np.float64)
    g[rad, rad] = 1.0
    f = gaussian_filter(g, sigma=sigma, mode="constant", cval=0.0, truncate=truncate)
    return float(f[rad, rad])


def kde_filtered_counts(
    query_xy: np.ndarray,
    sources: dict[str, np.ndarray],
    window_xy: np.ndarray,
    h_um: float,
    pixel: float = PIXEL_UM,
) -> tuple[dict[str, np.ndarray], np.ndarray, float]:
    """Sample edge-normalization denominator and filtered source counts.

    Returns filtered-count fields (kernel sums to 1), the kernel mass inside
    the rectangular window at each query, and the discrete kernel center k0.
    Intensity in cells per µm² is filtered_count / kernel_mass / pixel².
    """
    xmin, ymin, _xmax, _ymax, nx, ny = _grid_spec(window_xy, pixel)
    sigma = h_um / pixel
    den_grid = gaussian_filter(
        np.ones((ny, nx), dtype=np.float64),
        sigma=sigma,
        mode="constant",
        cval=0.0,
        truncate=TRUNCATE,
    )
    ix_f = np.clip((query_xy[:, 0] - xmin) / pixel - 0.5, 0, nx - 1)
    iy_f = np.clip((query_xy[:, 1] - ymin) / pixel - 0.5, 0, ny - 1)
    coords = np.vstack([iy_f, ix_f])
    den_q = map_coordinates(den_grid, coords, order=1, mode="nearest")
    out: dict[str, np.ndarray] = {}
    for name, xy in sources.items():
        grid = np.zeros((ny, nx), dtype=np.float64)
        if len(xy):
            ix = np.clip(np.floor((xy[:, 0] - xmin) / pixel).astype(np.int32), 0, nx - 1)
            iy = np.clip(np.floor((xy[:, 1] - ymin) / pixel).astype(np.int32), 0, ny - 1)
            np.add.at(grid, (iy, ix), 1.0)
        num = gaussian_filter(grid, sigma=sigma, mode="constant", cval=0.0, truncate=TRUNCATE)
        out[name] = map_coordinates(num, coords, order=1, mode="nearest")
    return out, den_q, _kernel_center(sigma)


def intensities_from_counts(
    counts: dict[str, np.ndarray], den_q: np.ndarray, k0: float, pixel: float
) -> dict[str, np.ndarray]:
    """Convert filtered counts to cells per µm², with leave-one-out all/tumor."""
    den = np.clip(den_q, 1e-6, None)
    area = pixel * pixel
    lam = {k: v / den / area for k, v in counts.items()}
    num_all = np.zeros(den_q.shape, dtype=np.float64)
    for v in counts.values():
        num_all = num_all + v
    lam["all"] = num_all / den / area
    lam["all_loo"] = np.maximum(num_all - k0, 0.0) / den / area
    if "tumor" in counts:
        lam["tumor_loo"] = np.maximum(counts["tumor"] - k0, 0.0) / den / area
    return lam


def exact_intensity(query, points, h, xmin, xmax, ymin, ymax) -> np.ndarray:
    """Continuous edge-corrected Gaussian intensity (cells per µm²)."""
    if len(points) == 0:
        return np.zeros(len(query), dtype=np.float64)
    dx = query[:, None, 0] - points[None, :, 0]
    dy = query[:, None, 1] - points[None, :, 1]
    dens = np.exp(-(dx * dx + dy * dy) / (2 * h * h)).sum(axis=1) / (2 * np.pi * h * h)
    mx = norm.cdf((xmax - query[:, 0]) / h) - norm.cdf((xmin - query[:, 0]) / h)
    my = norm.cdf((ymax - query[:, 1]) / h) - norm.cdf((ymin - query[:, 1]) / h)
    return dens / np.clip(mx * my, 1e-6, None)


def self_check() -> None:
    rng = np.random.default_rng(0)
    window = np.column_stack(
        [rng.uniform(0.0, 980.0, 600), rng.uniform(0.0, 650.0, 600)]
    )
    gx0, gy0, gx1, gy1, _nx, _ny = _grid_spec(window, PIXEL_UM)
    sources_xy = np.column_stack(
        [
            rng.uniform(gx0 + 30, gx0 + 0.45 * (gx1 - gx0), 70),
            rng.uniform(gy0 + 30, gy1 - 30, 70),
        ]
    )
    query = np.column_stack(
        [
            rng.uniform(gx0 + 20, gx1 - 20, 250),
            rng.uniform(gy0 + 20, gy1 - 20, 250),
        ]
    )
    h = 50.0
    counts, den_q, k0 = kde_filtered_counts(
        query, {"cd8": sources_xy}, window, h, PIXEL_UM
    )
    grid = intensities_from_counts(counts, den_q, k0, PIXEL_UM)["cd8"]
    exact = exact_intensity(query, sources_xy, h, gx0, gx1, gy0, gy1)
    rel = np.abs(grid - exact) / np.maximum(exact, 1e-12)
    interior = (
        (query[:, 0] > gx0 + 2 * h)
        & (query[:, 0] < gx1 - 2 * h)
        & (query[:, 1] > gy0 + 2 * h)
        & (query[:, 1] < gy1 - 2 * h)
        & (exact > np.median(exact))
    )
    if interior.sum() < 10:
        raise RuntimeError("self-check interior set too small")
    med = float(np.median(rel[interior]))
    corr = float(np.corrcoef(grid[interior], exact[interior])[0, 1])
    if med > 0.03 or corr < 0.995:
        raise RuntimeError(f"kernel self-check failed med_rel={med:.4f} corr={corr:.4f}")

    # Planted exclusion: CD8 on the left, CLDN4 increases to the right.
    cd8 = np.column_stack([rng.uniform(20, 180, 80), rng.uniform(30, 600, 80)])
    tumor = np.column_stack([rng.uniform(40, 900, 400), rng.uniform(30, 600, 400)])
    both = np.vstack([cd8, tumor])
    counts, den_q, k0 = kde_filtered_counts(tumor, {"CD8": cd8}, both, 50.0, PIXEL_UM)
    lam = intensities_from_counts(counts, den_q, k0, PIXEL_UM)["CD8"] * PER_100
    rho = fast_spearman(tumor[:, 0], lam)
    if not np.isfinite(rho) or rho > -0.4:
        raise RuntimeError(f"planted CD8 exclusion self-check failed rho={rho}")
    print(f"self-check ok  kernel med_rel={med:.4f} corr={corr:.4f} planted rho={rho:.3f}", flush=True)


def fast_spearman(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    a = np.asarray(a)[m]
    b = np.asarray(b)[m]
    if a.size < 30:
        return float("nan")
    if np.unique(a).size < 2 or np.unique(b).size < 2:
        return float("nan")
    ra = rankdata(a)
    rb = rankdata(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt(np.dot(ra, ra) * np.dot(rb, rb))
    if denom == 0:
        return float("nan")
    return float(np.dot(ra, rb) / denom)


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    r_xy = fast_spearman(x, y)
    r_xz = fast_spearman(x, z)
    r_yz = fast_spearman(y, z)
    if not (np.isfinite(r_xy) and np.isfinite(r_xz) and np.isfinite(r_yz)):
        return float("nan")
    den = np.sqrt(max(0.0, 1.0 - r_xz * r_xz) * max(0.0, 1.0 - r_yz * r_yz))
    if den < 0.05:
        return float("nan")
    return float((r_xy - r_xz * r_yz) / den)


def load_cells() -> pd.DataFrame:
    df = pd.read_csv(DATA)
    missing = sorted(set(df["cell_type"]) - set(SUPER_OF))
    if missing:
        raise RuntimeError(f"unmapped author cell types: {missing}")
    genes = pd.read_csv(GENES)
    gene_set = set(genes.iloc[:, 0].astype(str))
    if "CLDN4" not in gene_set:
        raise RuntimeError("CLDN4 is absent from data/cosmx_960_genes.csv; stopping.")
    if "CLDN4" not in df.columns:
        raise RuntimeError("CLDN4 column missing from the Giotto cell table; stopping.")
    df["super"] = df["cell_type"].map(SUPER_OF)
    df["x_um"] = df["x"].to_numpy(dtype=np.float64) * UM_PER_MM
    df["y_um"] = df["y"].to_numpy(dtype=np.float64) * UM_PER_MM
    df["is_tumor"] = df["super"].eq("tumor")
    df["is_cd8"] = df["super"].eq("CD8")
    return df


def field_tumor_table(df: pd.DataFrame) -> pd.DataFrame:
    """Evaluate density fields at every author tumor cell, within FOV."""
    tumor = df.loc[df["is_tumor"], ["cell_ID", "sample", "patient", "fov", "cell_type", "CLDN4"]].copy()
    n = len(tumor)
    store = {
        "den50": np.zeros(n),
        "den100": np.zeros(n),
    }
    for h in BANDWIDTHS:
        tag = f"{int(h)}"
        for name in ("cd8", "nk", "tumor", "tumor_loo", "all", "all_loo"):
            store[f"{name}{tag}"] = np.zeros(n)
        for name in SUPER_ORDER:
            store[f"lam{tag}_{name}"] = np.zeros(n)

    # positional index into `tumor`
    tumor = tumor.reset_index(drop=True)
    key_to_pos = {}
    for (sample, fov), sub in tumor.groupby(["sample", "fov"], sort=False):
        key_to_pos[(sample, int(fov))] = sub.index.to_numpy()

    n_fov = 0
    for (sample, fov), fov_df in df.groupby(["sample", "fov"], sort=False):
        n_fov += 1
        pos = key_to_pos.get((sample, int(fov)))
        if pos is None or pos.size == 0:
            continue
        window = fov_df[["x_um", "y_um"]].to_numpy(dtype=np.float64)
        if window.shape[0] < 5:
            continue
        q = window[fov_df["is_tumor"].to_numpy()]
        if q.shape[0] != pos.size:
            raise RuntimeError(f"tumor count mismatch in {sample} FOV {fov}")
        q_ids = fov_df.loc[fov_df["is_tumor"], "cell_ID"].to_numpy()
        if not np.array_equal(q_ids, tumor.loc[pos, "cell_ID"].to_numpy()):
            raise RuntimeError(f"tumor order mismatch in {sample} FOV {fov}")
        sources = {}
        supers = fov_df["super"].to_numpy()
        xy = window
        for name in SUPER_ORDER:
            sources[name] = xy[supers == name]
        for h in BANDWIDTHS:
            counts, den_q, k0 = kde_filtered_counts(q, sources, window, h, PIXEL_UM)
            lam = intensities_from_counts(counts, den_q, k0, PIXEL_UM)
            tag = f"{int(h)}"
            store[f"den{tag}"][pos] = den_q
            store[f"cd8{tag}"][pos] = lam["CD8"] * PER_100
            store[f"nk{tag}"][pos] = lam["NK"] * PER_100
            store[f"tumor{tag}"][pos] = lam["tumor"] * PER_100
            store[f"tumor_loo{tag}"][pos] = lam["tumor_loo"] * PER_100
            store[f"all{tag}"][pos] = lam["all"] * PER_100
            store[f"all_loo{tag}"][pos] = lam["all_loo"] * PER_100
            for name in SUPER_ORDER:
                store[f"lam{tag}_{name}"][pos] = lam[name] * PER_100
        if n_fov % 40 == 0:
            print(f"  fields {n_fov} FOVs", flush=True)

    for k, v in store.items():
        tumor[k] = v
    for tag in ("50", "100"):
        tumor[f"cyto{tag}"] = tumor[f"cd8{tag}"] + tumor[f"nk{tag}"]
        all_loo = tumor[f"all_loo{tag}"].to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = tumor[f"cd8{tag}"].to_numpy() / all_loo
        frac[~np.isfinite(frac) | (all_loo <= 0)] = np.nan
        tumor[f"p_cd8_{tag}"] = frac
    print(f"density fields evaluated at {len(tumor)} tumor cells in {n_fov} FOVs", flush=True)
    return tumor


def _section_codes(sample: pd.Series) -> np.ndarray:
    lookup = {s: i for i, s in enumerate(SAMPLES)}
    return sample.map(lookup).to_numpy(dtype=np.int32)


def section_rho(y: np.ndarray, x: np.ndarray, codes: np.ndarray) -> np.ndarray:
    out = np.full(len(SAMPLES), np.nan)
    for s in range(len(SAMPLES)):
        m = codes == s
        if int(m.sum()) < 30:
            continue
        out[s] = fast_spearman(y[m], x[m])
    return out


def section_partial(y: np.ndarray, x: np.ndarray, z: np.ndarray, codes: np.ndarray) -> np.ndarray:
    out = np.full(len(SAMPLES), np.nan)
    for s in range(len(SAMPLES)):
        m = codes == s
        if int(m.sum()) < 30:
            continue
        out[s] = partial_spearman(y[m], x[m], z[m])
    return out


def section_contrast(y: np.ndarray, x: np.ndarray, codes: np.ndarray) -> dict[str, np.ndarray]:
    """High = CLDN4 >= max(section p75, 1); low = CLDN4 == 0. Delta = median high − median low."""
    delta = np.full(len(SAMPLES), np.nan)
    med_high = np.full(len(SAMPLES), np.nan)
    med_low = np.full(len(SAMPLES), np.nan)
    n_high = np.zeros(len(SAMPLES), dtype=int)
    n_low = np.zeros(len(SAMPLES), dtype=int)
    thr = np.full(len(SAMPLES), np.nan)
    for s in range(len(SAMPLES)):
        m = np.flatnonzero(codes == s)
        if m.size < 30:
            continue
        ok = np.isfinite(x[m]) & np.isfinite(y[m])
        m = m[ok]
        if m.size < 30:
            continue
        c = y[m]
        p75 = float(np.quantile(c, 0.75))
        threshold = p75 if p75 >= 1.0 else 1.0
        thr[s] = threshold
        high = c >= threshold
        low = c == 0
        n_high[s] = int(high.sum())
        n_low[s] = int(low.sum())
        if high.sum() < 20 or low.sum() < 20:
            continue
        med_high[s] = float(np.median(x[m][high]))
        med_low[s] = float(np.median(x[m][low]))
        delta[s] = med_high[s] - med_low[s]
    return {
        "delta": delta,
        "med_high": med_high,
        "med_low": med_low,
        "n_high": n_high,
        "n_low": n_low,
        "thr": thr,
    }


def _mean(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(a.mean())


def patient_means(section_delta: np.ndarray) -> np.ndarray:
    patients = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
    sec_pat = [s.split("_")[0] for s in SAMPLES]
    out = []
    for p in patients:
        vals = [section_delta[i] for i, sp in enumerate(sec_pat) if sp == p]
        out.append(float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan"))
    return np.array(out, dtype=np.float64)


def _wilcoxon_p(delta: np.ndarray) -> float | None:
    d = delta[np.isfinite(delta)]
    d = d[d != 0]
    if d.size < 5:
        return None
    try:
        return float(wilcoxon(d, alternative="two-sided", method="exact", zero_method="wilcox").pvalue)
    except ValueError:
        return float(wilcoxon(d, alternative="two-sided", zero_method="wilcox").pvalue)


def _two_sided_perm_p(obs: float, null: np.ndarray) -> float:
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    return float((1 + np.sum(np.abs(null) >= abs(obs))) / (null.size + 1))


def _one_sided_less_p(obs: float, null: np.ndarray) -> float:
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    return float((1 + np.sum(null <= obs)) / (null.size + 1))


def permutation_test(frame: pd.DataFrame, h_tag: str, n_perm: int, seed: int) -> dict:
    y = frame["CLDN4"].to_numpy(dtype=np.float64)
    codes = _section_codes(frame["sample"])
    fov_key = frame["sample"].astype(str) + "|" + frame["fov"].astype(str)
    _fov_labels, fov_codes = np.unique(fov_key.to_numpy(), return_inverse=True)
    groups = [np.flatnonzero(fov_codes == f) for f in range(fov_codes.max() + 1)]
    predictors = {
        "lam": frame[f"cd8{h_tag}"].to_numpy(dtype=np.float64),
        "frac": frame[f"p_cd8_{h_tag}"].to_numpy(dtype=np.float64),
        "cyto": frame[f"cyto{h_tag}"].to_numpy(dtype=np.float64),
    }
    z = frame[f"tumor_loo{h_tag}"].to_numpy(dtype=np.float64)
    obs = {}
    for key, x in predictors.items():
        rho = section_rho(y, x, codes)
        contr = section_contrast(y, x, codes)
        part = section_partial(y, x, z, codes) if key in ("lam", "frac") else None
        obs[key] = {"rho": rho, "contrast": contr, "partial": part}

    null_rho = {k: np.full((n_perm, len(SAMPLES)), np.nan) for k in predictors}
    null_delta = {k: np.full((n_perm, len(SAMPLES)), np.nan) for k in predictors}
    null_partial = {k: np.full((n_perm, len(SAMPLES)), np.nan) for k in ("lam", "frac")}
    rng = np.random.default_rng(seed)
    y_perm = np.empty_like(y)
    for b in range(n_perm):
        y_perm[:] = y
        for g in groups:
            if g.size > 1:
                y_perm[g] = rng.permutation(y[g])
        for key, x in predictors.items():
            null_rho[key][b] = section_rho(y_perm, x, codes)
            null_delta[key][b] = section_contrast(y_perm, x, codes)["delta"]
            if key in null_partial:
                null_partial[key][b] = section_partial(y_perm, x, z, codes)
        if (b + 1) % 250 == 0:
            print(f"  perm {h_tag} {b + 1}/{n_perm}", flush=True)

    def pack(key: str) -> dict:
        rho = obs[key]["rho"]
        delta = obs[key]["contrast"]["delta"]
        mean_rho = _mean(rho)
        mean_delta = _mean(delta)
        null_mean_rho = np.nanmean(null_rho[key], axis=1)
        null_mean_delta = np.nanmean(null_delta[key], axis=1)
        pat = patient_means(delta)
        null_pat = np.array([_mean(patient_means(null_delta[key][b])) for b in range(n_perm)])
        packed = {
            "rho_by_section": rho,
            "delta_by_section": delta,
            "med_high": obs[key]["contrast"]["med_high"],
            "med_low": obs[key]["contrast"]["med_low"],
            "n_high": obs[key]["contrast"]["n_high"],
            "n_low": obs[key]["contrast"]["n_low"],
            "thr": obs[key]["contrast"]["thr"],
            "mean_rho": mean_rho,
            "mean_delta": mean_delta,
            "rho_p_two": _two_sided_perm_p(mean_rho, null_mean_rho),
            "rho_p_less": _one_sided_less_p(mean_rho, null_mean_rho),
            "delta_p_two": _two_sided_perm_p(mean_delta, null_mean_delta),
            "delta_p_less": _one_sided_less_p(mean_delta, null_mean_delta),
            "n_sections_rho_neg": int(np.sum(rho < 0)),
            "n_sections_rho_finite": int(np.sum(np.isfinite(rho))),
            "n_sections_delta_neg": int(np.sum(delta < 0)),
            "n_sections_delta_finite": int(np.sum(np.isfinite(delta))),
            "patient_delta": pat,
            "patient_mean_delta": _mean(pat),
            "patient_delta_p_two": _two_sided_perm_p(_mean(pat), null_pat),
            "n_patients_delta_neg": int(np.sum(pat < 0)),
            "n_patients_finite": int(np.sum(np.isfinite(pat))),
            "wilcoxon_delta_p": _wilcoxon_p(delta),
            "wilcoxon_rho_p": _wilcoxon_p(rho),
            "null_rho": null_rho[key],
            "null_delta": null_delta[key],
            "null_mean_rho": null_mean_rho,
            "null_mean_delta": null_mean_delta,
        }
        if key in null_partial:
            part = obs[key]["partial"]
            mean_part = _mean(part)
            null_mean_part = np.nanmean(null_partial[key], axis=1)
            packed["partial_by_section"] = part
            packed["mean_partial"] = mean_part
            packed["partial_p_two"] = _two_sided_perm_p(mean_part, null_mean_part)
            packed["partial_p_less"] = _one_sided_less_p(mean_part, null_mean_part)
            packed["n_sections_partial_neg"] = int(np.sum(part < 0))
            packed["null_mean_partial"] = null_mean_part
        n_neg = packed["n_patients_delta_neg"]
        n_pat = packed["n_patients_finite"]
        if n_pat > 0:
            packed["patient_sign_p_two"] = float(binomtest(n_neg, n_pat, 0.5, alternative="two-sided").pvalue)
        else:
            packed["patient_sign_p_two"] = float("nan")
        return packed

    return {k: pack(k) for k in predictors}


def type_spearman(frame: pd.DataFrame, h_tag: str) -> pd.DataFrame:
    y = frame["CLDN4"].to_numpy(dtype=np.float64)
    codes = _section_codes(frame["sample"])
    rows = []
    for name in SUPER_ORDER:
        col = f"tumor_loo{h_tag}" if name == "tumor" else f"lam{h_tag}_{name}"
        rho = section_rho(y, frame[col].to_numpy(dtype=np.float64), codes)
        rec = {"cell_type": name, "mean_rho": _mean(rho)}
        for i, sample in enumerate(SAMPLES):
            rec[sample] = float(rho[i]) if np.isfinite(rho[i]) else None
        rows.append(rec)
    return pd.DataFrame(rows)


def fov_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, fov), sub in frame.groupby(["sample", "fov"], sort=False):
        if len(sub) < 30:
            continue
        y = sub["CLDN4"].to_numpy(dtype=np.float64)
        rows.append(
            {
                "sample": sample,
                "patient": sub["patient"].iloc[0],
                "fov": int(fov),
                "n_tumor": int(len(sub)),
                "rho_lam50": fast_spearman(y, sub["cd850"].to_numpy(dtype=np.float64)),
                "rho_p50": fast_spearman(y, sub["p_cd8_50"].to_numpy(dtype=np.float64)),
            }
        )
    return pd.DataFrame(rows)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items() if k not in ("null_rho", "null_delta", "null_mean_rho", "null_mean_delta", "null_mean_partial")}
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return None if not np.isfinite(v) else v
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None:
        return None
    return obj


def write_section_csv(path: Path, packed_lam: dict, packed_frac: dict, packed_cyto: dict, n_mask: np.ndarray) -> None:
    rows = []
    for i, sample in enumerate(SAMPLES):
        rows.append(
            {
                "sample": sample,
                "patient": sample.split("_")[0],
                "n_tumor_mask": int(n_mask[i]),
                "cldn4_high_threshold": packed_lam["thr"][i],
                "n_high": int(packed_lam["n_high"][i]),
                "n_low": int(packed_lam["n_low"][i]),
                "median_cd8_density_high": packed_lam["med_high"][i],
                "median_cd8_density_low": packed_lam["med_low"][i],
                "delta_cd8_density": packed_lam["delta_by_section"][i],
                "rho_cd8_density": packed_lam["rho_by_section"][i],
                "partial_rho_cd8_density_given_tumor": packed_lam.get("partial_by_section", [np.nan] * 8)[i],
                "median_cd8_fraction_high": packed_frac["med_high"][i],
                "median_cd8_fraction_low": packed_frac["med_low"][i],
                "delta_cd8_fraction": packed_frac["delta_by_section"][i],
                "rho_cd8_fraction": packed_frac["rho_by_section"][i],
                "partial_rho_cd8_fraction_given_tumor": packed_frac.get("partial_by_section", [np.nan] * 8)[i],
                "rho_cd8nk_density": packed_cyto["rho_by_section"][i],
                "delta_cd8nk_density": packed_cyto["delta_by_section"][i],
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "figure.dpi": 140,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def savefig(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIG / f"{name}.png")
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)


def plot_example(df: pd.DataFrame, fov_stats: pd.DataFrame) -> dict:
    usable = fov_stats.dropna(subset=["rho_lam50"])
    usable = usable[(usable["n_tumor"] >= 80)]
    # CD8 count filter needs the parent table
    counts = (
        df.groupby(["sample", "fov"])["is_cd8"].sum().rename("n_cd8").reset_index()
    )
    usable = usable.merge(counts, on=["sample", "fov"], how="left")
    usable = usable[usable["n_cd8"] >= 20]
    if usable.empty:
        return {}
    med = float(usable["rho_lam50"].median())
    usable = usable.assign(dist=(usable["rho_lam50"] - med).abs())
    pick = usable.sort_values(["dist", "sample", "fov"]).iloc[0]
    sample, fov = pick["sample"], int(pick["fov"])
    sub = df[(df["sample"] == sample) & (df["fov"] == fov)]
    window = sub[["x_um", "y_um"]].to_numpy(dtype=np.float64)
    tumor = sub[sub["is_tumor"]]
    q = tumor[["x_um", "y_um"]].to_numpy(dtype=np.float64)
    sources = {}
    xy = window
    supers = sub["super"].to_numpy()
    for name in ("CD8", "tumor"):
        sources[name] = xy[supers == name]
    counts, den_q, k0 = kde_filtered_counts(q, sources, window, PRIMARY_H, PIXEL_UM)
    # Rebuild the image on the grid for display.
    xmin, ymin, xmax, ymax, nx, ny = _grid_spec(window, PIXEL_UM)
    sigma = PRIMARY_H / PIXEL_UM
    grid = np.zeros((ny, nx), dtype=np.float64)
    cd8_xy = sources["CD8"]
    if len(cd8_xy):
        ix = np.clip(np.floor((cd8_xy[:, 0] - xmin) / PIXEL_UM).astype(np.int32), 0, nx - 1)
        iy = np.clip(np.floor((cd8_xy[:, 1] - ymin) / PIXEL_UM).astype(np.int32), 0, ny - 1)
        np.add.at(grid, (iy, ix), 1.0)
    num = gaussian_filter(grid, sigma=sigma, mode="constant", cval=0.0, truncate=TRUNCATE)
    den = gaussian_filter(np.ones_like(grid), sigma=sigma, mode="constant", cval=0.0, truncate=TRUNCATE)
    field = (num / np.clip(den, 1e-6, None)) / (PIXEL_UM ** 2) * PER_100
    field = np.where(den >= PRIMARY_MASK, field, np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    cldn = tumor["CLDN4"].to_numpy()
    vmax = max(1.0, float(np.quantile(cldn, 0.95)))
    ax = axes[0]
    sc = ax.scatter(
        tumor["x_um"],
        tumor["y_um"],
        c=cldn,
        s=6,
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )
    fig.colorbar(sc, ax=ax, label="Tumor CLDN4 count", fraction=0.046, pad=0.02)
    ax.set_title("Author tumor cells")
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")

    ax = axes[1]
    im = ax.imshow(
        field,
        origin="lower",
        extent=[xmin, xmax, ymin, ymax],
        cmap="magma",
        interpolation="nearest",
    )
    ax.scatter(cd8_xy[:, 0], cd8_xy[:, 1], s=8, c="cyan", linewidths=0, label="CD8", rasterized=True)
    fig.colorbar(im, ax=ax, label="Predicted CD8 density\n(cells / (100 µm)²)", fraction=0.046, pad=0.02)
    ax.set_title(f"CD8 density field, h = {int(PRIMARY_H)} µm")
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.legend(loc="upper right", frameon=False, markerscale=1.5)
    fig.suptitle(
        f"{sample} FOV {fov}  ·  within-FOV Spearman {pick['rho_lam50']:.2f} "
        f"(median FOV {med:.2f}; typical, not the extreme)",
        fontsize=11,
    )
    savefig(fig, "example_fov_cd8_density")
    return {
        "sample": sample,
        "fov": fov,
        "rho_lam50": float(pick["rho_lam50"]),
        "median_fov_rho": med,
        "n_tumor": int(pick["n_tumor"]),
        "n_cd8": int(pick["n_cd8"]),
    }


def plot_paired(primary: dict) -> None:
    """Per-section high−low contrast. Absolute levels differ by orders of magnitude across sections."""
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.3), constrained_layout=True)
    y = np.arange(len(SAMPLES))
    for ax, key, title, xlab in (
        (
            axes[0],
            "lam",
            "CD8 density",
            "Median density, CLDN4-high − undetected\n(cells / (100 µm)²)",
        ),
        (
            axes[1],
            "frac",
            "CD8 fraction of the density field",
            "Median fraction, CLDN4-high − undetected",
        ),
    ):
        delta = np.asarray(primary[key]["delta_by_section"], dtype=float)
        colors = ["#E45756" if np.isfinite(v) and v < 0 else "#4C78A8" for v in delta]
        ax.axvline(0, color="0.35", lw=0.8)
        ax.barh(y, delta, color=colors, height=0.72)
        ax.set_yticks(y)
        ax.set_yticklabels(SAMPLES)
        ax.invert_yaxis()
        ax.set_xlabel(xlab)
        ax.set_title(title)
    fig.suptitle(
        f"h = {int(PRIMARY_H)} µm  ·  red: CLDN4-high has less predicted CD8",
        fontsize=11,
    )
    savefig(fig, "section_paired_cd8")


def plot_forest(packed: dict, name: str, title: str) -> None:
    rho = packed["rho_by_section"]
    null = packed["null_rho"]
    lo = np.nanpercentile(null, 2.5, axis=0)
    hi = np.nanpercentile(null, 97.5, axis=0)
    y = np.arange(len(SAMPLES))
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.hlines(y, lo, hi, color="0.65", lw=2)
    ax.plot(rho, y, "o", color="#E45756", ms=6)
    ax.set_yticks(y)
    ax.set_yticklabels(SAMPLES)
    ax.set_xlabel("Within-section Spearman ρ (CLDN4, predicted CD8)")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.text(
        0.0,
        -0.18,
        "Bars: 2.5–97.5% null from within-FOV CLDN4 permutations, not confidence intervals.",
        transform=ax.transAxes,
        fontsize=8,
        color="0.35",
    )
    savefig(fig, name)


def plot_heatmap(type_df: pd.DataFrame) -> None:
    mat = type_df[SAMPLES].to_numpy(dtype=float)
    v = float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 0.2
    v = max(v, 0.15)
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(type_df)))
    ax.set_yticklabels(type_df["cell_type"])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6, color="black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ with tumor CLDN4")
    ax.set_title(f"Tumor CLDN4 vs each cell-type density field ({int(PRIMARY_H)} µm)")
    savefig(fig, "type_density_heatmap")


def plot_deciles(frame: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    colors = {
        "Lung5": "#4C78A8",
        "Lung6": "#F58518",
        "Lung9": "#54A24B",
        "Lung12": "#B279A2",
        "Lung13": "#72B7B2",
    }
    curves = []
    for sample in SAMPLES:
        sub = frame[frame["sample"] == sample]
        x = sub["cd850"].to_numpy(dtype=np.float64)
        y = sub["CLDN4"].to_numpy(dtype=np.float64)
        if np.unique(x).size < 10:
            continue
        pct = rankdata(x, method="average") / len(x)
        bins = np.minimum((pct * 10).astype(int), 9)
        means = np.array([y[bins == b].mean() if np.any(bins == b) else np.nan for b in range(10)])
        curves.append(means)
        ax.plot(
            np.arange(1, 11),
            means,
            color=colors[sample.split("_")[0]],
            lw=1.2,
            alpha=0.9,
            label=sample,
        )
    if curves:
        med = np.nanmedian(np.vstack(curves), axis=0)
        ax.plot(np.arange(1, 11), med, color="black", lw=2.2, label="Median of sections", zorder=3)
    ax.set_xticks(range(1, 11))
    ax.set_xlabel("Decile of predicted CD8 density (1 = lowest)")
    ax.set_ylabel("Mean tumor CLDN4 count")
    ax.set_title(f"CLDN4 vs predicted CD8 density, h = {int(PRIMARY_H)} µm")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    savefig(fig, "cldn4_vs_predicted_cd8_deciles")


def plot_null(packed: dict) -> None:
    obs = packed["mean_rho"]
    null = packed["null_mean_rho"]
    fig, ax = plt.subplots(figsize=(5.6, 3.6), constrained_layout=True)
    ax.hist(null, bins=30, color="#9ecae1", edgecolor="white")
    ax.axvline(obs, color="#E45756", lw=2, label=f"Observed mean ρ = {obs:.3f}")
    ax.axvline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("Mean of 8 section Spearman ρ")
    ax.set_ylabel("Permutations")
    ax.set_title("Null: CLDN4 shuffled within FOV")
    ax.legend(frameon=False, fontsize=8)
    savefig(fig, "permutation_null_mean_rho")


def n_mask_by_section(frame: pd.DataFrame) -> np.ndarray:
    codes = _section_codes(frame["sample"])
    return np.array([(codes == s).sum() for s in range(len(SAMPLES))], dtype=int)


def main() -> None:
    apply_style()
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    self_check()
    print("loading cells", flush=True)
    df = load_cells()
    digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
    inv = {
        "n_cells": int(len(df)),
        "n_tumor": int(df["is_tumor"].sum()),
        "n_cd8": int(df["is_cd8"].sum()),
        "n_fov": int(df.groupby(["sample", "fov"]).ngroups),
        "n_patients": int(df["patient"].nunique()),
        "sha256_cells": digest,
        "cldn4_on_panel": True,
        "source": (
            "NanoString public Giotto object for He et al. 2022 CosMx NSCLC "
            "(author cell types, CLDN4 counts, millimetre coordinates). "
            "Same export as the Ripley analysis table."
        ),
    }
    print(json.dumps({k: inv[k] for k in ("n_cells", "n_tumor", "n_cd8", "n_fov")}, indent=2), flush=True)
    tumor = field_tumor_table(df)
    # Fraction coherence: raw type densities should sum to the all-cell density.
    part_sum = np.zeros(len(tumor))
    for name in SUPER_ORDER:
        part_sum += tumor["lam50_" + name].to_numpy()
    ratio = part_sum / np.clip(tumor["all50"].to_numpy(), 1e-8, None)
    interior = tumor["den50"].to_numpy() >= PRIMARY_MASK
    med_ratio = float(np.median(ratio[interior]))
    print(f"partition/all density ratio median (mask 0.8) = {med_ratio:.4f}", flush=True)
    if not 0.95 <= med_ratio <= 1.05:
        raise RuntimeError(f"type densities do not sum to the total field (median ratio {med_ratio})")

    primary = tumor[tumor["den50"] >= PRIMARY_MASK].copy()
    secondary = tumor[tumor["den100"] >= PRIMARY_MASK].copy()
    sens = tumor[tumor["den50"] >= SENS_MASK].copy()
    print(
        f"primary mask n={len(primary)} / {len(tumor)}; h100 n={len(secondary)}; sens n={len(sens)}",
        flush=True,
    )

    print("permutation, 50 µm", flush=True)
    perm50 = permutation_test(primary, "50", N_PERM, SEED)
    print("permutation, 100 µm", flush=True)
    perm100 = permutation_test(secondary, "100", N_PERM, SEED + 1)
    print("sensitivity summary (mask 0.2, no extra permutation grid)", flush=True)
    # Observed-only sensitivity using the same contrast helpers via a 0-perm call is wasteful.
    # Compute directly.
    y = sens["CLDN4"].to_numpy(dtype=np.float64)
    codes = _section_codes(sens["sample"])
    sens_rho = section_rho(y, sens["cd850"].to_numpy(dtype=np.float64), codes)
    sens_delta = section_contrast(y, sens["cd850"].to_numpy(dtype=np.float64), codes)["delta"]
    sens_partial = section_partial(
        y,
        sens["cd850"].to_numpy(dtype=np.float64),
        sens["tumor_loo50"].to_numpy(dtype=np.float64),
        codes,
    )

    types = type_spearman(primary, "50")
    types.to_csv(TAB / "type_spearman_h50.csv", index=False)
    fovs = fov_table(primary)
    fovs.to_csv(TAB / "fov_spearman_h50.csv", index=False)
    write_section_csv(
        TAB / "section_h50.csv",
        perm50["lam"],
        perm50["frac"],
        perm50["cyto"],
        n_mask_by_section(primary),
    )
    write_section_csv(
        TAB / "section_h100.csv",
        perm100["lam"],
        perm100["frac"],
        perm100["cyto"],
        n_mask_by_section(secondary),
    )

    slim_cols = [
        "cell_ID",
        "sample",
        "patient",
        "fov",
        "cell_type",
        "CLDN4",
        "den50",
        "cd850",
        "p_cd8_50",
        "cyto50",
        "tumor_loo50",
        "all_loo50",
        "den100",
        "cd8100",
        "p_cd8_100",
        "cyto100",
        "tumor_loo100",
    ]
    tumor[slim_cols].to_csv(TAB / "tumor_cell_density.csv.gz", index=False)

    example = plot_example(df, fovs)
    plot_paired(perm50)
    plot_forest(
        perm50["lam"],
        "spearman_forest_h50",
        f"CLDN4 vs CD8 density, h = {int(PRIMARY_H)} µm",
    )
    plot_forest(
        perm100["lam"],
        "spearman_forest_h100",
        "CLDN4 vs CD8 density, h = 100 µm",
    )
    plot_heatmap(types)
    plot_deciles(primary)
    plot_null(perm50["lam"])

    summary = {
        "visium_companion": False,
        "visium_gate": (
            "He et al. Nat Biotechnol 2022 and the public CosMx NSCLC release "
            "(SMI flat files plus the Giotto object) contain these 8 sections only. "
            "No matched Visium, Visium HD, or other spot-level companion is in that release. "
            "Tangram, RCTD, and cell2location were not run. Unrelated public Visium LUAD "
            "cohorts are not a companion for these sections."
        ),
        "method": {
            "bandwidths_um": list(BANDWIDTHS),
            "primary_bandwidth_um": PRIMARY_H,
            "pixel_um": PIXEL_UM,
            "primary_mask_kernel_mass": PRIMARY_MASK,
            "sensitivity_mask": SENS_MASK,
            "n_perm": N_PERM,
            "seed": SEED,
            "cd8_labels": ["T CD8 memory", "T CD8 naive"],
            "tumor_labels": "author labels starting with 'tumor'",
            "high_definition": "CLDN4 >= section 75th percentile, and at least 1 count",
            "low_definition": "CLDN4 == 0",
            "coordinates": "Giotto millimetres × 1000 = µm; within-FOV only",
            "primary_predictor": "edge-corrected Gaussian CD8 intensity, cells per (100 µm)^2",
            "fraction_predictor": "CD8 intensity / leave-one-out total intensity",
        },
        "inventory": inv,
        "n_primary_mask": int(len(primary)),
        "n_h100_mask": int(len(secondary)),
        "n_sens_mask": int(len(sens)),
        "partition_over_all_median_ratio": med_ratio,
        "primary_h50_density": _jsonable(perm50["lam"]),
        "primary_h50_fraction": _jsonable(perm50["frac"]),
        "primary_h50_cyto": _jsonable(perm50["cyto"]),
        "h100_density": _jsonable(perm100["lam"]),
        "h100_fraction": _jsonable(perm100["frac"]),
        "sensitivity_mask_0.2": {
            "mean_rho": _mean(sens_rho),
            "rho_by_section": sens_rho,
            "mean_delta": _mean(sens_delta),
            "delta_by_section": sens_delta,
            "mean_partial": _mean(sens_partial),
            "n_sections_rho_neg": int(np.sum(sens_rho < 0)),
            "n_sections_delta_neg": int(np.sum(sens_delta < 0)),
        },
        "example_fov": example,
        "type_mean_rho": {r.cell_type: r.mean_rho for r in types.itertuples()},
    }
    # numpy arrays inside sensitivity
    summary["sensitivity_mask_0.2"] = _jsonable(summary["sensitivity_mask_0.2"])
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lam = perm50["lam"]
    print(
        "PRIMARY mean section Spearman",
        round(lam["mean_rho"], 4),
        "two-sided perm p",
        lam["rho_p_two"],
        "sections neg",
        f"{lam['n_sections_rho_neg']}/{lam['n_sections_rho_finite']}",
        flush=True,
    )
    print(
        "PRIMARY mean delta density",
        round(lam["mean_delta"], 4),
        "two-sided perm p",
        lam["delta_p_two"],
        "sections neg",
        f"{lam['n_sections_delta_neg']}/{lam['n_sections_delta_finite']}",
        "patients neg",
        f"{lam['n_patients_delta_neg']}/{lam['n_patients_finite']}",
        flush=True,
    )
    print("partial mean", lam.get("mean_partial"), "p", lam.get("partial_p_two"), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
