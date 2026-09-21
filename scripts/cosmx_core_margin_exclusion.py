#!/usr/bin/env python3
"""CosMx NSCLC: CLDN4 exclusion stratified by tumor core vs margin.

Official public object: CellCharter clustered CosMx human NSCLC
(figshare 25976224, He et al. 2022; 8 sections / 5 patients).

Question
    Is CLDN4-high exclusion of cytotoxic cells stronger in tumor core
    than at the margin, where core/margin is distance to immune-rich niches?

Pre-specified before outcome inspection
    Immune-rich niches (author CellCharter `niche`, not generic stroma):
        immune, lymphoid structure, macrophages, myeloid-enriched stroma,
        neutrophils, plasmablast-enriched stroma.
    Distance: Euclidean µm from each tumor cell to the nearest cell in an
    immune-rich niche (0 if the tumor cell itself sits in one).
    Margin: distance <= 50 µm. Core: distance >= 150 µm.
    Tumor: author cell_type starting with "tumor ".
    Cytotoxic: NK + T CD8 memory + T CD8 naive.
    CLDN4-high: log-normalized CLDN4 above the within-section tumor median.
    Primary metric: mean cytotoxic-neighbor count at 50 µm.
    Exclusion: mean_high - mean_low < 0 (ratio high/low < 1).
    Stronger in core: (delta_core - delta_margin) < 0, with delta_core < 0.
    Test unit: section (n=8) paired Wilcoxon; patient (n=5) sign test.
    Null: shuffle CLDN4-high labels within FOV (999 perms).
    Secondary: 100 µm; cytotoxic fraction; contact rate; distance tertiles;
    niche-label core (tumor interior) vs margin (tumor-stroma boundary);
    immune-rich definition that also includes generic stroma.
    Infiltration-depth AUC (0–200 µm trapezoid of the 50 µm cytotoxic count
    versus distance) is reported and is claimed only if it is lower for
    CLDN4-high than CLDN4-low. The opposite direction is not claimed.

Pixel size 0.18 µm/px (CosMx NSCLC SMI). No private 8-KL. No Visium.
"""

from __future__ import annotations

import json
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "data", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_core_margin")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

UM_PER_PX = 0.18
MARGIN_UM = 50.0
CORE_UM = 150.0
RADII = (50.0, 100.0)
PRIMARY_R = 50.0
MIN_ARM = 40
N_PERM = 999
SEED = 25976224
AUC_MAX_UM = 200.0
AUC_BIN_UM = 20.0
MIN_BIN = 30

IMMUNE_RICH = (
    "immune",
    "lymphoid structure",
    "macrophages",
    "myeloid-enriched stroma",
    "neutrophils",
    "plasmablast-enriched stroma",
)
IMMUNE_RICH_PLUS_STROMA = IMMUNE_RICH + ("stroma",)
CYTOTOXIC = ("NK", "T CD8 memory", "T CD8 naive")
SAMPLE_ORDER = (
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
)


def _decode(arr) -> list[str]:
    out = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode())
        else:
            out.append(str(x))
    return out


def load_obs_and_cldn4(path: str) -> dict:
    with h5py.File(path, "r") as f:
        genes = _decode(f["var"]["_index"][:])
        if "CLDN4" not in genes:
            raise SystemExit("CLDN4 is absent from the CosMx panel. Stopping.")
        gene = genes.index("CLDN4")
        n = int(f["obs"]["_index"].shape[0])
        if n != 765771:
            raise SystemExit(f"expected 765771 cells, found {n}")
        indptr = f["X"]["indptr"][:]
        cldn4 = np.zeros(n, dtype=np.float32)
        indices = f["X"]["indices"]
        data = f["X"]["data"]
        step = 40000
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
            cldn4[start + rows] = data[a + hit]
        xy = np.asarray(f["obsm"]["spatial"][:], dtype=np.float64) * UM_PER_PX
        obs = {
            "cell_type": f["obs"]["cell_type"]["codes"][:],
            "niche": f["obs"]["niche"]["codes"][:],
            "sample": f["obs"]["sample"]["codes"][:],
            "patient": f["obs"]["patient"]["codes"][:],
            "fov": f["obs"]["fov"][:],
        }
        cats = {
            "cell_type": _decode(f["obs"]["cell_type"]["categories"][:]),
            "niche": _decode(f["obs"]["niche"]["categories"][:]),
            "sample": _decode(f["obs"]["sample"]["categories"][:]),
            "patient": _decode(f["obs"]["patient"]["categories"][:]),
        }
    return {"n": n, "cldn4": cldn4, "xy": xy, "obs": obs, "cats": cats, "gene_index": gene}


def codes_for(cats: list[str], names: tuple[str, ...] | list[str]) -> np.ndarray:
    missing = [nm for nm in names if nm not in cats]
    if missing:
        raise SystemExit(f"missing categories {missing} in {cats}")
    return np.array([cats.index(nm) for nm in names], dtype=np.int16)


def tumor_codes(cats: list[str]) -> np.ndarray:
    ids = [i for i, nm in enumerate(cats) if nm.startswith("tumor ")]
    if len(ids) != 5:
        raise SystemExit(f"expected 5 tumor labels, found {[cats[i] for i in ids]}")
    return np.array(ids, dtype=np.int16)


def arm_means(values: np.ndarray, high: np.ndarray, mask: np.ndarray) -> dict:
    h = mask & high
    low = mask & ~high
    nh = int(h.sum())
    nl = int(low.sum())
    out = {
        "n_high": nh,
        "n_low": nl,
        "usable": nh >= MIN_ARM and nl >= MIN_ARM,
    }
    if nh == 0 or nl == 0:
        out.update(mean_high=np.nan, mean_low=np.nan, delta=np.nan, ratio=np.nan)
        return out
    mh = float(values[h].mean())
    ml = float(values[low].mean())
    out["mean_high"] = mh
    out["mean_low"] = ml
    out["delta"] = mh - ml
    out["ratio"] = (mh / ml) if ml >= 0.05 else np.nan
    return out


def contact_table(contact: np.ndarray, high: np.ndarray, mask: np.ndarray) -> dict:
    """2x2 of any-cytotoxic-neighbor vs CLDN4 arm. Haldane 0.5 used for OR."""
    h = mask & high
    low = mask & ~high
    a = float((contact & h).sum())
    b = float((~contact & h).sum())
    c = float((contact & low).sum())
    d = float((~contact & low).sum())
    or_ = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return {
        "a_high_pos": a,
        "b_high_neg": b,
        "c_low_pos": c,
        "d_low_neg": d,
        "rate_high": a / (a + b) if (a + b) else np.nan,
        "rate_low": c / (c + d) if (c + d) else np.nan,
        "or": or_,
        "log_or": float(np.log(or_)),
    }


def mh_or(rows: list[dict]) -> float:
    num = 0.0
    den = 0.0
    for r in rows:
        if "a_high_pos" in r:
            a, b, c, d = r["a_high_pos"], r["b_high_neg"], r["c_low_pos"], r["d_low_neg"]
        else:
            a, b, c, d = r["a"], r["b"], r["c"], r["d"]
        n = a + b + c + d
        if n <= 0:
            continue
        num += a * d / n
        den += b * c / n
    if den <= 0:
        return np.nan
    return float(num / den)


def wilcoxon_pair(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    diff = a[m] - b[m]
    out = {"n": int(m.sum()), "n_a_lt_b": int(np.sum(diff < 0)), "n_a_gt_b": int(np.sum(diff > 0))}
    if m.sum() < 5 or np.allclose(diff, 0):
        out["p_two"] = np.nan
        out["stat"] = np.nan
        return out
    res = stats.wilcoxon(a[m], b[m], alternative="two-sided", zero_method="wilcox")
    out["stat"] = float(res.statistic)
    out["p_two"] = float(res.pvalue)
    return out


def sign_test(n_success: int, n: int, alternative: str = "two-sided") -> float:
    if n <= 0:
        return float("nan")
    return float(stats.binomtest(n_success, n, 0.5, alternative=alternative).pvalue)


def depth_auc(distance: np.ndarray, values: np.ndarray, high: np.ndarray) -> dict:
    edges = np.arange(0.0, AUC_MAX_UM + AUC_BIN_UM, AUC_BIN_UM)
    centers = []
    y_high = []
    y_low = []
    n_high = []
    n_low = []
    for lo, hi_edge in zip(edges[:-1], edges[1:]):
        m = (distance >= lo) & (distance < hi_edge)
        nh = int((m & high).sum())
        nl = int((m & ~high).sum())
        n_high.append(nh)
        n_low.append(nl)
        if nh < MIN_BIN or nl < MIN_BIN:
            continue
        centers.append(0.5 * (lo + hi_edge))
        y_high.append(float(values[m & high].mean()))
        y_low.append(float(values[m & ~high].mean()))
    if len(centers) < 4:
        return {"auc_high": np.nan, "auc_low": np.nan, "n_bins": len(centers)}
    x = np.asarray(centers)
    auc_h = float(np.trapezoid(y_high, x))
    auc_l = float(np.trapezoid(y_low, x))
    return {
        "auc_high": auc_h,
        "auc_low": auc_l,
        "delta_low_minus_high": auc_l - auc_h,
        "n_bins": len(centers),
        "centers": centers,
        "y_high": y_high,
        "y_low": y_low,
        "n_high": n_high,
        "n_low": n_low,
    }


def summarize_section(
    count: np.ndarray,
    frac: np.ndarray,
    contact: np.ndarray,
    high: np.ndarray,
    strata: dict[str, np.ndarray],
) -> dict:
    out = {}
    for name, mask in strata.items():
        out[name] = {
            "count": arm_means(count, high, mask),
            "fraction": arm_means(frac, high, mask),
            "contact": contact_table(contact, high, mask),
        }
    return out


def interaction_from_summary(summary: dict, stratum_a: str, stratum_b: str, metric: str) -> dict:
    """a minus b. For core minus margin, negative delta means a larger drop in core."""
    A = summary[stratum_a][metric]
    B = summary[stratum_b][metric]
    if not (A["usable"] and B["usable"]):
        return {"usable": False, "delta_a": np.nan, "delta_b": np.nan, "interaction": np.nan,
                "ratio_a": np.nan, "ratio_b": np.nan, "ratio_interaction": np.nan}
    return {
        "usable": True,
        "delta_a": A["delta"],
        "delta_b": B["delta"],
        "interaction": A["delta"] - B["delta"],
        "ratio_a": A["ratio"],
        "ratio_b": B["ratio"],
        "ratio_interaction": (
            A["ratio"] - B["ratio"]
            if np.isfinite(A["ratio"]) and np.isfinite(B["ratio"])
            else np.nan
        ),
        "mean_high_a": A["mean_high"],
        "mean_low_a": A["mean_low"],
        "mean_high_b": B["mean_high"],
        "mean_low_b": B["mean_low"],
        "n_high_a": A["n_high"],
        "n_low_a": A["n_low"],
        "n_high_b": B["n_high"],
        "n_low_b": B["n_low"],
    }


def prepare_sample(bundle: dict, sample_code: int, rich_ids: np.ndarray) -> dict | None:
    obs = bundle["obs"]
    cats = bundle["cats"]
    idx = np.flatnonzero(obs["sample"] == sample_code)
    if idx.size == 0:
        return None
    ct = obs["cell_type"][idx]
    niche = obs["niche"][idx]
    xy = bundle["xy"][idx]
    tumor_ids = tumor_codes(cats["cell_type"])
    cyto_ids = codes_for(cats["cell_type"], CYTOTOXIC)
    is_tumor = np.isin(ct, tumor_ids)
    is_cyto = np.isin(ct, cyto_ids)
    is_rich = np.isin(niche, rich_ids)
    if is_tumor.sum() < 100 or is_rich.sum() < 20 or is_cyto.sum() < 20:
        return None
    tum_xy = xy[is_tumor]
    tree_rich = cKDTree(xy[is_rich])
    dist, _ = tree_rich.query(tum_xy, k=1, workers=4)
    dist = np.asarray(dist, dtype=np.float64)
    tum_in_rich = is_rich[is_tumor]
    dist[tum_in_rich] = 0.0
    tree_cyto = cKDTree(xy[is_cyto])
    tree_all = cKDTree(xy)
    expr = bundle["cldn4"][idx][is_tumor]
    med = float(np.median(expr))
    high = expr > med
    fov = obs["fov"][idx][is_tumor]
    niche_tum = niche[is_tumor]
    counts = {}
    fracs = {}
    contacts = {}
    for radius in RADII:
        n_cyto = np.asarray(
            tree_cyto.query_ball_point(tum_xy, r=radius, return_length=True, workers=4),
            dtype=np.float64,
        )
        n_all = np.asarray(
            tree_all.query_ball_point(tum_xy, r=radius, return_length=True, workers=4),
            dtype=np.float64,
        )
        n_all = np.clip(n_all - 1.0, 0.0, None)
        frac = np.divide(n_cyto, n_all, out=np.zeros_like(n_cyto), where=n_all > 0)
        counts[radius] = n_cyto
        fracs[radius] = frac
        contacts[radius] = n_cyto > 0
    return {
        "dist": dist,
        "high": high,
        "fov": fov,
        "expr": expr,
        "median_cldn4": med,
        "frac_high": float(high.mean()),
        "n_tumor": int(is_tumor.sum()),
        "n_cyto": int(is_cyto.sum()),
        "n_rich": int(is_rich.sum()),
        "niche": niche_tum,
        "counts": counts,
        "fracs": fracs,
        "contacts": contacts,
        "xy": tum_xy,
        "in_rich": tum_in_rich,
    }


def strata_distance(dist: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "all": np.ones(dist.shape[0], dtype=bool),
        "margin": dist <= MARGIN_UM,
        "core": dist >= CORE_UM,
        "buffer": (dist > MARGIN_UM) & (dist < CORE_UM),
    }


def strata_tertile(dist: np.ndarray) -> dict[str, np.ndarray]:
    q1, q2 = np.quantile(dist, [1 / 3, 2 / 3])
    return {
        "near": dist <= q1,
        "mid": (dist > q1) & (dist <= q2),
        "far": dist > q2,
        "q_near": float(q1),
        "q_far": float(q2),
    }


def strata_niche(niche: np.ndarray, niche_cats: list[str]) -> dict[str, np.ndarray]:
    interior = niche_cats.index("tumor interior")
    boundary = niche_cats.index("tumor-stroma boundary")
    return {
        "niche_core": niche == interior,
        "niche_margin": niche == boundary,
    }


def permute_high(high: np.ndarray, fov: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = high.copy()
    for fv in np.unique(fov):
        m = fov == fv
        out[m] = rng.permutation(high[m])
    return out


def section_interaction(prep: dict, mask_core: np.ndarray, mask_margin: np.ndarray, radius: float) -> dict:
    high = prep["high"]
    sm = summarize_section(
        prep["counts"][radius],
        prep["fracs"][radius],
        prep["contacts"][radius],
        high,
        {"core": mask_core, "margin": mask_margin, "all": np.ones(high.shape[0], dtype=bool)},
    )
    inter_c = interaction_from_summary(sm, "core", "margin", "count")
    inter_f = interaction_from_summary(sm, "core", "margin", "fraction")
    return {"summary": sm, "count": inter_c, "fraction": inter_f}


def nanmean(xs: list[float]) -> float:
    arr = np.asarray(xs, dtype=float)
    if np.isfinite(arr).sum() == 0:
        return float("nan")
    return float(np.nanmean(arr))


def perm_p(obs: float, null: np.ndarray) -> dict:
    """Tails are versus the permutation distribution, which is not centered at 0.

    Absolute gaps are larger where counts are higher, so a label shuffle already
    expects a more negative delta at the margin. p_two is recentered on the
    null median. p_one_core_stronger is the lower tail (core minus margin).
    """
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return {
            "p_one_core_stronger": np.nan,
            "p_one_margin_stronger": np.nan,
            "p_two": np.nan,
            "n": int(null.size),
        }
    p_low = (1 + np.sum(null <= obs)) / (1 + null.size)
    p_high = (1 + np.sum(null >= obs)) / (1 + null.size)
    center = float(np.median(null))
    p_two = (1 + np.sum(np.abs(null - center) >= abs(obs - center))) / (1 + null.size)
    return {
        "p_one_core_stronger": float(p_low),
        "p_one_margin_stronger": float(p_high),
        "p_two": float(p_two),
        "n": int(null.size),
        "null_mean": float(null.mean()),
        "null_median": center,
        "null_sd": float(null.std()),
    }


def run_definition(bundle: dict, rich_names: tuple[str, ...], label: str) -> dict:
    cats = bundle["cats"]
    rich_ids = codes_for(cats["niche"], rich_names)
    sample_cats = cats["sample"]
    patient_cats = cats["patient"]
    niche_cats = cats["niche"]
    rows = []
    auc_rows = []
    curve_rows = []
    preps = {}
    print(f"\n=== definition {label} ===", flush=True)
    for sname in SAMPLE_ORDER:
        code = sample_cats.index(sname)
        print(f"  {sname}: query neighbors", flush=True)
        prep = prepare_sample(bundle, code, rich_ids)
        if prep is None:
            print(f"  {sname}: skipped (too few cells)", flush=True)
            continue
        preps[sname] = prep
        dist = prep["dist"]
        st = strata_distance(dist)
        tert = strata_tertile(dist)
        nich = strata_niche(prep["niche"], niche_cats)
        patient = patient_cats[int(bundle["obs"]["patient"][bundle["obs"]["sample"] == code][0])]
        for radius in RADII:
            summary = summarize_section(
                prep["counts"][radius],
                prep["fracs"][radius],
                prep["contacts"][radius],
                prep["high"],
                {
                    "all": st["all"],
                    "margin": st["margin"],
                    "core": st["core"],
                    "buffer": st["buffer"],
                    "near": tert["near"],
                    "far": tert["far"],
                    "niche_core": nich["niche_core"],
                    "niche_margin": nich["niche_margin"],
                },
            )
            for stratum, block in summary.items():
                c = block["count"]
                fr = block["fraction"]
                ct = block["contact"]
                rows.append({
                    "definition": label,
                    "sample": sname,
                    "patient": patient,
                    "radius_um": radius,
                    "stratum": stratum,
                    "n_tumor": prep["n_tumor"],
                    "n_high_arm": c["n_high"],
                    "n_low_arm": c["n_low"],
                    "usable": c["usable"],
                    "count_high": c["mean_high"],
                    "count_low": c["mean_low"],
                    "count_delta": c["delta"],
                    "count_ratio": c["ratio"],
                    "frac_high": fr["mean_high"],
                    "frac_low": fr["mean_low"],
                    "frac_delta": fr["delta"],
                    "frac_ratio": fr["ratio"],
                    "contact_high": ct["rate_high"],
                    "contact_low": ct["rate_low"],
                    "contact_or": ct["or"],
                    "a": ct["a_high_pos"],
                    "b": ct["b_high_neg"],
                    "c": ct["c_low_pos"],
                    "d": ct["d_low_neg"],
                    "cldn4_median": prep["median_cldn4"],
                    "frac_cldn4_high": prep["frac_high"],
                    "dist_median": float(np.median(dist)),
                    "dist_p90": float(np.quantile(dist, 0.9)),
                    "n_core_cells": int(st["core"].sum()),
                    "n_margin_cells": int(st["margin"].sum()),
                    "tertile_near_um": tert["q_near"],
                    "tertile_far_um": tert["q_far"],
                })
        auc = depth_auc(dist, prep["counts"][PRIMARY_R], prep["high"])
        auc_rows.append({
            "definition": label,
            "sample": sname,
            "patient": patient,
            "auc_high": auc.get("auc_high", np.nan),
            "auc_low": auc.get("auc_low", np.nan),
            "delta_low_minus_high": auc.get("delta_low_minus_high", np.nan),
            "n_bins": auc.get("n_bins", 0),
        })
        if "centers" in auc:
            for x, yh, yl in zip(auc["centers"], auc["y_high"], auc["y_low"]):
                curve_rows.append({
                    "definition": label,
                    "sample": sname,
                    "patient": patient,
                    "distance_um": x,
                    "count_high": yh,
                    "count_low": yl,
                })
        print(
            f"    tumor={prep['n_tumor']} core={int(st['core'].sum())} "
            f"margin={int(st['margin'].sum())} median_d={np.median(dist):.1f}",
            flush=True,
        )
    return {
        "label": label,
        "rows": rows,
        "auc_rows": auc_rows,
        "curve_rows": curve_rows,
        "preps": preps,
        "rich_names": rich_names,
    }


def fov_permutation(preps: dict, radius: float) -> dict:
    """FOV-blocked shuffle of CLDN4 labels. Statistic is the mean section interaction."""
    rng = np.random.default_rng(SEED)
    names = list(preps)
    obs_delta = []
    obs_ratio = []
    for name in names:
        prep = preps[name]
        dist = prep["dist"]
        inter = section_interaction(prep, dist >= CORE_UM, dist <= MARGIN_UM, radius)
        obs_delta.append(inter["count"]["interaction"])
        obs_ratio.append(inter["count"]["ratio_interaction"])
    obs_delta_mean = nanmean(obs_delta)
    obs_ratio_mean = nanmean(obs_ratio)
    null_delta = np.empty(N_PERM)
    null_ratio = np.empty(N_PERM)
    print(f"  permutations n={N_PERM} radius={radius}", flush=True)
    for p in range(N_PERM):
        dvals = []
        rvals = []
        for name in names:
            prep = preps[name]
            high0 = prep["high"]
            prep["high"] = permute_high(high0, prep["fov"], rng)
            inter = section_interaction(prep, prep["dist"] >= CORE_UM, prep["dist"] <= MARGIN_UM, radius)
            dvals.append(inter["count"]["interaction"])
            rvals.append(inter["count"]["ratio_interaction"])
            prep["high"] = high0
        null_delta[p] = nanmean(dvals)
        null_ratio[p] = nanmean(rvals)
        if (p + 1) % 250 == 0:
            print(f"    perm {p+1}", flush=True)
    return {
        "radius_um": radius,
        "obs_interaction_delta_mean": obs_delta_mean,
        "obs_interaction_ratio_mean": obs_ratio_mean,
        "section_interaction_delta": obs_delta,
        "section_interaction_ratio": obs_ratio,
        "delta": perm_p(obs_delta_mean, null_delta),
        "ratio": perm_p(obs_ratio_mean, null_ratio),
    }


def paired_tests(df: pd.DataFrame, definition: str, radius: float, stratum_a: str, stratum_b: str) -> dict:
    sub = df[(df.definition == definition) & (df.radius_um == radius)]
    a = sub[sub.stratum == stratum_a].set_index("sample")
    b = sub[sub.stratum == stratum_b].set_index("sample")
    common = a.index.intersection(b.index)
    a = a.loc[common]
    b = b.loc[common]
    usable = a["usable"] & b["usable"]
    a = a.loc[usable]
    b = b.loc[usable]

    def pack(col_delta: str, col_ratio: str) -> dict:
        wa = wilcoxon_pair(a[col_delta].to_numpy(), np.zeros(len(a)))
        # Wilcoxon of delta vs 0 is wilcoxon on the delta values themselves
        delta = a[col_delta].to_numpy(dtype=float)
        delta_b = b[col_delta].to_numpy(dtype=float)
        inter = delta - delta_b
        m = np.isfinite(inter)
        w_inter = {"n": int(m.sum()), "p_two": np.nan, "n_core_stronger": int(np.sum(inter[m] < 0))}
        if m.sum() >= 5 and not np.allclose(inter[m], 0):
            res = stats.wilcoxon(inter[m], alternative="two-sided", zero_method="wilcox")
            w_inter["p_two"] = float(res.pvalue)
            w_inter["stat"] = float(res.statistic)
        w_a = {"n": int(np.isfinite(delta).sum()), "p_two": np.nan,
               "n_exclusion": int(np.sum(delta[np.isfinite(delta)] < 0))}
        dfin = delta[np.isfinite(delta)]
        if dfin.size >= 5 and not np.allclose(dfin, 0):
            res = stats.wilcoxon(dfin, alternative="two-sided", zero_method="wilcox")
            w_a["p_two"] = float(res.pvalue)
        dfin_b = delta_b[np.isfinite(delta_b)]
        w_b = {"n": int(dfin_b.size), "p_two": np.nan,
               "n_exclusion": int(np.sum(dfin_b < 0))}
        if dfin_b.size >= 5 and not np.allclose(dfin_b, 0):
            res = stats.wilcoxon(dfin_b, alternative="two-sided", zero_method="wilcox")
            w_b["p_two"] = float(res.pvalue)
        ratio_a = a[col_ratio].to_numpy(dtype=float)
        ratio_b = b[col_ratio].to_numpy(dtype=float)
        rinter = ratio_a - ratio_b
        rm = np.isfinite(rinter)
        w_r = {"n": int(rm.sum()), "p_two": np.nan, "n_core_stronger": int(np.sum(rinter[rm] < 0))}
        if rm.sum() >= 5 and not np.allclose(rinter[rm], 0):
            res = stats.wilcoxon(rinter[rm], alternative="two-sided", zero_method="wilcox")
            w_r["p_two"] = float(res.pvalue)
        # patient: unweighted mean of sections
        pat_rows = []
        for patient, g_idx in a.groupby("patient").groups.items():
            sa = a.loc[g_idx]
            sb = b.loc[sa.index]
            pat_rows.append({
                "patient": patient,
                "delta_a": float(sa[col_delta].mean()),
                "delta_b": float(sb[col_delta].mean()),
                "interaction": float((sa[col_delta] - sb[col_delta]).mean()),
                "ratio_a": float(np.nanmean(sa[col_ratio])),
                "ratio_b": float(np.nanmean(sb[col_ratio])),
            })
        pat = pd.DataFrame(pat_rows)
        n_pat = int(len(pat))
        n_pat_excl_a = int(np.sum(pat["delta_a"] < 0)) if n_pat else 0
        n_pat_excl_b = int(np.sum(pat["delta_b"] < 0)) if n_pat else 0
        n_pat_core = int(np.sum(pat["interaction"] < 0)) if n_pat else 0
        return {
            "n_sections": int(len(a)),
            "stratum_a": {
                "mean_delta": float(np.nanmean(delta)) if len(delta) else np.nan,
                "mean_ratio": float(np.nanmean(ratio_a)) if len(ratio_a) else np.nan,
                "wilcoxon_delta_vs_0": w_a,
                "patients_exclusion": f"{n_pat_excl_a}/{n_pat}",
                "patient_sign_p_one": sign_test(n_pat_excl_a, n_pat, "greater"),
            },
            "stratum_b": {
                "mean_delta": float(np.nanmean(delta_b)) if len(delta_b) else np.nan,
                "mean_ratio": float(np.nanmean(ratio_b)) if len(ratio_b) else np.nan,
                "wilcoxon_delta_vs_0": w_b,
                "patients_exclusion": f"{n_pat_excl_b}/{n_pat}",
                "patient_sign_p_one": sign_test(n_pat_excl_b, n_pat, "greater"),
            },
            "interaction_delta": w_inter,
            "interaction_ratio": w_r,
            "mean_interaction_delta": float(np.nanmean(inter)) if m.any() else np.nan,
            "mean_interaction_ratio": float(np.nanmean(rinter)) if rm.any() else np.nan,
            "patients_core_stronger_delta": f"{n_pat_core}/{n_pat}",
            "patient_sign_p_one_core_stronger": sign_test(n_pat_core, n_pat, "greater"),
            "patient_sign_p_two_core_stronger": sign_test(n_pat_core, n_pat, "two-sided"),
            "patients": pat.to_dict(orient="records"),
        }

    tests = {
        "count": pack("count_delta", "count_ratio"),
        "fraction": pack("frac_delta", "frac_ratio"),
    }
    # contact OR paired
    if len(a):
        log_or_a = np.log(a["contact_or"].to_numpy(dtype=float))
        log_or_b = np.log(b["contact_or"].to_numpy(dtype=float))
        inter = log_or_a - log_or_b
        m = np.isfinite(inter)
        w = {"n": int(m.sum()), "p_two": np.nan, "n_core_stronger": int(np.sum(inter[m] < 0))}
        if m.sum() >= 5 and not np.allclose(inter[m], 0):
            res = stats.wilcoxon(inter[m], alternative="two-sided", zero_method="wilcox")
            w["p_two"] = float(res.pvalue)
        tests["contact_log_or"] = {
            "mh_or_a": mh_or(a.to_dict(orient="records")),
            "mh_or_b": mh_or(b.to_dict(orient="records")),
            "mean_log_or_a": float(np.nanmean(log_or_a)),
            "mean_log_or_b": float(np.nanmean(log_or_b)),
            "interaction": w,
            "mean_interaction_log_or": float(np.nanmean(inter)) if m.any() else np.nan,
        }
    tests["samples"] = list(a.index)
    return tests


def auc_tests(auc_df: pd.DataFrame, definition: str) -> dict:
    sub = auc_df[auc_df.definition == definition].copy()
    h = sub["auc_high"].to_numpy(dtype=float)
    low = sub["auc_low"].to_numpy(dtype=float)
    m = np.isfinite(h) & np.isfinite(low)
    diff = h[m] - low[m]  # positive => high has MORE cytotoxic area (opposite of exclusion)
    out = {
        "n": int(m.sum()),
        "mean_auc_high": float(np.mean(h[m])) if m.any() else np.nan,
        "mean_auc_low": float(np.mean(low[m])) if m.any() else np.nan,
        "n_high_lt_low": int(np.sum(diff < 0)),
        "n_high_gt_low": int(np.sum(diff > 0)),
        "p_two": np.nan,
    }
    if m.sum() >= 5 and not np.allclose(diff, 0):
        res = stats.wilcoxon(h[m], low[m], alternative="two-sided", zero_method="wilcox")
        out["p_two"] = float(res.pvalue)
    # patient
    pat_hi = []
    for patient, g in sub.groupby("patient"):
        hh = g["auc_high"].to_numpy(dtype=float)
        ll = g["auc_low"].to_numpy(dtype=float)
        mm = np.isfinite(hh) & np.isfinite(ll)
        if mm.sum() == 0:
            continue
        pat_hi.append(float(np.mean(hh[mm] - ll[mm])))
    n_pat = len(pat_hi)
    n_opp = int(np.sum(np.asarray(pat_hi) > 0))  # high AUC greater = opposite
    n_same = int(np.sum(np.asarray(pat_hi) < 0))
    out["patients_high_gt_low"] = f"{n_opp}/{n_pat}"
    out["patients_high_lt_low"] = f"{n_same}/{n_pat}"
    # claim only if high AUC is lower (exclusion direction) in a majority and p<0.05
    exclusion_direction = out["n_high_lt_low"] >= 6 and np.isfinite(out["p_two"]) and out["p_two"] < 0.05
    opposite = out["n_high_gt_low"] > out["n_high_lt_low"]
    out["claim"] = bool(exclusion_direction and not opposite)
    out["opposite"] = bool(opposite or (m.any() and out["mean_auc_high"] > out["mean_auc_low"]))
    return out


def niche_composition(bundle: dict) -> pd.DataFrame:
    cats = bundle["cats"]
    niche = bundle["obs"]["niche"]
    ct = bundle["obs"]["cell_type"]
    tumor_ids = set(tumor_codes(cats["cell_type"]).tolist())
    cyto_ids = set(codes_for(cats["cell_type"], CYTOTOXIC).tolist())
    immune_names = {
        "B-cell", "NK", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive",
        "Treg", "mDC", "macrophage", "mast", "monocyte", "neutrophil", "pDC", "plasmablast",
    }
    immune_ids = set(codes_for(cats["cell_type"], tuple(immune_names)).tolist())
    rows = []
    for i, name in enumerate(cats["niche"]):
        m = niche == i
        n = int(m.sum())
        rows.append({
            "niche": name,
            "n": n,
            "tumor_frac": float(np.isin(ct[m], list(tumor_ids)).mean()) if n else np.nan,
            "cytotoxic_frac": float(np.isin(ct[m], list(cyto_ids)).mean()) if n else np.nan,
            "immune_frac": float(np.isin(ct[m], list(immune_ids)).mean()) if n else np.nan,
            "immune_rich_primary": name in IMMUNE_RICH,
        })
    return pd.DataFrame(rows)


def plot_paired(df: pd.DataFrame, path: str) -> None:
    sub = df[(df.definition == "immune_rich") & (df.radius_um == PRIMARY_R)]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), constrained_layout=True)
    for ax, radius, title in (
        (axes[0], 50.0, "50 µm cytotoxic neighbors"),
        (axes[1], 100.0, "100 µm cytotoxic neighbors"),
    ):
        s = df[(df.definition == "immune_rich") & (np.isclose(df.radius_um, radius))]
        core = s[s.stratum == "core"].set_index("sample")
        margin = s[s.stratum == "margin"].set_index("sample")
        samples = [nm for nm in SAMPLE_ORDER if nm in core.index and nm in margin.index]
        for i, nm in enumerate(samples):
            if not (bool(core.loc[nm, "usable"]) and bool(margin.loc[nm, "usable"])):
                continue
            y0 = margin.loc[nm, "count_ratio"]
            y1 = core.loc[nm, "count_ratio"]
            if not (np.isfinite(y0) and np.isfinite(y1)):
                continue
            ax.plot([0, 1], [y0, y1], color="#4d4d4d", lw=1.0, zorder=1)
            ax.scatter([0], [y0], color="#e08214", s=36, zorder=2)
            ax.scatter([1], [y1], color="#542788", s=36, zorder=2)
        ax.axhline(1.0, color="#999999", ls="--", lw=0.8)
        ax.set_xticks([0, 1], ["Margin\n(≤50 µm)", "Core\n(≥150 µm)"])
        ax.set_ylabel("CLDN4-high / CLDN4-low\nmean cytotoxic count")
        ax.set_title(title)
        ax.set_xlim(-0.35, 1.35)
    fig.suptitle("Section-paired cytotoxic-neighbor ratio", fontsize=12)
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def plot_deltas(df: pd.DataFrame, path: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    s = df[(df.definition == "immune_rich") & (np.isclose(df.radius_um, PRIMARY_R))]
    core = s[s.stratum == "core"].set_index("sample")
    margin = s[s.stratum == "margin"].set_index("sample")
    labels = []
    y_m, y_c = [], []
    for nm in SAMPLE_ORDER:
        if nm not in core.index:
            continue
        if not (bool(core.loc[nm, "usable"]) and bool(margin.loc[nm, "usable"])):
            continue
        labels.append(nm.replace(" ", "\n"))
        y_m.append(margin.loc[nm, "count_delta"])
        y_c.append(core.loc[nm, "count_delta"])
    x = np.arange(len(labels))
    ax.axhline(0, color="#999999", lw=0.8)
    ax.scatter(x - 0.12, y_m, color="#e08214", s=40, label="Margin ≤50 µm", zorder=2)
    ax.scatter(x + 0.12, y_c, color="#542788", s=40, label="Core ≥150 µm", zorder=2)
    for i in range(len(labels)):
        ax.plot([x[i] - 0.12, x[i] + 0.12], [y_m[i], y_c[i]], color="#4d4d4d", lw=0.8)
    ax.set_xticks(x, labels, fontsize=8)
    ax.set_ylabel("Δ cytotoxic neighbors (high − low)")
    ax.set_title("50 µm exclusion gap by section")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def plot_curves(curve: pd.DataFrame, auc_df: pd.DataFrame, path: str) -> None:
    sub = curve[curve.definition == "immune_rich"]
    if sub.empty:
        return
    g = sub.groupby("distance_um")[["count_high", "count_low"]].mean()
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), constrained_layout=True)
    ax = axes[0]
    ax.plot(g.index, g["count_high"], color="#b2182b", marker="o", ms=4, label="CLDN4-high")
    ax.plot(g.index, g["count_low"], color="#2166ac", marker="o", ms=4, label="CLDN4-low")
    ax.set_xlabel("Distance to immune-rich niche (µm)")
    ax.set_ylabel("Mean cytotoxic neighbors in 50 µm")
    ax.set_title("Section-mean depth curve")
    ax.legend(frameon=False)
    ax = axes[1]
    auc = auc_df[auc_df.definition == "immune_rich"].set_index("sample")
    for nm in SAMPLE_ORDER:
        if nm not in auc.index:
            continue
        y0 = float(auc.loc[nm, "auc_low"])
        y1 = float(auc.loc[nm, "auc_high"])
        if not (np.isfinite(y0) and np.isfinite(y1)):
            continue
        ax.plot([0, 1], [y0, y1], color="#4d4d4d", lw=1.0)
        ax.scatter([0], [y0], color="#2166ac", s=36, zorder=2)
        ax.scatter([1], [y1], color="#b2182b", s=36, zorder=2)
    ax.set_xticks([0, 1], ["CLDN4-low", "CLDN4-high"])
    ax.set_ylabel("AUC of 50 µm cytotoxic count\nfrom 0 to 200 µm")
    ax.set_title("Section AUC (not a primary claim)")
    ax.set_xlim(-0.35, 1.35)
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def plot_map(prep: dict, sample: str, path: str) -> None:
    xy = prep["xy"]
    dist = prep["dist"]
    rng = np.random.default_rng(1)
    # plot a deterministic subset if huge
    n = len(xy)
    if n > 80000:
        take = rng.choice(n, 80000, replace=False)
    else:
        take = np.arange(n)
    fig, ax = plt.subplots(figsize=(6.2, 5.4), constrained_layout=True)
    d = dist[take]
    cols = np.full(take.shape[0], "#d9d9d9")
    cols[d <= MARGIN_UM] = "#e08214"
    cols[d >= CORE_UM] = "#542788"
    ax.scatter(xy[take, 0], xy[take, 1], c=cols, s=0.35, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.set_title(f"{sample} tumor cells\norange margin ≤50 µm, purple core ≥150 µm")
    ax.set_xlabel("µm")
    ax.set_ylabel("µm")
    # scale bar
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    bar_x = x0 + 0.06 * (x1 - x0)
    bar_y = y0 + 0.06 * (y1 - y0)
    ax.plot([bar_x, bar_x + 500], [bar_y, bar_y], color="black", lw=2)
    ax.text(bar_x, bar_y + 0.03 * (y1 - y0), "500 µm", fontsize=8, ha="left")
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def json_clean(obj):
    if isinstance(obj, dict):
        return {k: json_clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_clean(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if not np.isfinite(v) else v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, float):
        return None if not np.isfinite(obj) else obj
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def main() -> None:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    if not os.path.isfile(H5):
        raise SystemExit(f"missing {H5}. Download figshare file 46841842.")
    print("loading obs + CLDN4", flush=True)
    bundle = load_obs_and_cldn4(H5)
    print(f"cells={bundle['n']} CLDN4>0 {(bundle['cldn4']>0).mean():.3f}", flush=True)
    comp = niche_composition(bundle)
    comp.to_csv(os.path.join(TAB, "niche_composition.csv"), index=False)

    primary = run_definition(bundle, IMMUNE_RICH, "immune_rich")
    plus = run_definition(bundle, IMMUNE_RICH_PLUS_STROMA, "immune_rich_plus_stroma")
    # drop bulky coordinates from the stroma sensitivity
    plus["preps"] = {}

    print("FOV-blocked permutation (primary definition, 50 µm)", flush=True)
    perm50 = fov_permutation(primary["preps"], 50.0)
    print("FOV-blocked permutation (primary definition, 100 µm)", flush=True)
    perm100 = fov_permutation(primary["preps"], 100.0)

    rows = primary["rows"] + plus["rows"]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "section_stratum_metrics.csv"), index=False)
    auc_df = pd.DataFrame(primary["auc_rows"] + plus["auc_rows"])
    auc_df.to_csv(os.path.join(TAB, "infiltration_depth_auc.csv"), index=False)
    curve = pd.DataFrame(primary["curve_rows"] + plus["curve_rows"])
    curve.to_csv(os.path.join(TAB, "depth_curves.csv"), index=False)

    tests = {
        "primary_core_vs_margin_50": paired_tests(df, "immune_rich", 50.0, "core", "margin"),
        "primary_core_vs_margin_100": paired_tests(df, "immune_rich", 100.0, "core", "margin"),
        "tertile_far_vs_near_50": paired_tests(df, "immune_rich", 50.0, "far", "near"),
        "tertile_far_vs_near_100": paired_tests(df, "immune_rich", 100.0, "far", "near"),
        "niche_label_50": paired_tests(df, "immune_rich", 50.0, "niche_core", "niche_margin"),
        "niche_label_100": paired_tests(df, "immune_rich", 100.0, "niche_core", "niche_margin"),
        "plus_stroma_50": paired_tests(df, "immune_rich_plus_stroma", 50.0, "core", "margin"),
        "unstratified_50": None,
    }
    # unstratified reference: delta vs 0 on stratum all
    all50 = df[(df.definition == "immune_rich") & (np.isclose(df.radius_um, 50)) & (df.stratum == "all")]
    all100 = df[(df.definition == "immune_rich") & (np.isclose(df.radius_um, 100)) & (df.stratum == "all")]
    def unstrat(block: pd.DataFrame) -> dict:
        delta = block["count_delta"].to_numpy(dtype=float)
        ratio = block["count_ratio"].to_numpy(dtype=float)
        m = np.isfinite(delta)
        w = {"n": int(m.sum()), "p_two": np.nan, "n_exclusion": int(np.sum(delta[m] < 0))}
        if m.sum() >= 5 and not np.allclose(delta[m], 0):
            res = stats.wilcoxon(delta[m], alternative="two-sided", zero_method="wilcox")
            w["p_two"] = float(res.pvalue)
        # patient
        pat = block.groupby("patient")["count_delta"].mean()
        n_pat = int(len(pat))
        n_ex = int(np.sum(pat < 0))
        return {
            "mean_ratio": float(np.nanmean(ratio)),
            "mean_delta": float(np.nanmean(delta)),
            "section_exclusion": f"{int(np.sum(delta[m] < 0))}/{int(m.sum())}",
            "wilcoxon": w,
            "patient_exclusion": f"{n_ex}/{n_pat}",
            "patient_sign_p_one": sign_test(n_ex, n_pat, "greater"),
            "ratios": {r["sample"]: r["count_ratio"] for _, r in block.iterrows()},
            "deltas": {r["sample"]: r["count_delta"] for _, r in block.iterrows()},
        }
    tests["unstratified_50"] = unstrat(all50)
    tests["unstratified_100"] = unstrat(all100)
    tests["auc_primary"] = auc_tests(auc_df, "immune_rich")
    tests["auc_plus_stroma"] = auc_tests(auc_df, "immune_rich_plus_stroma")
    tests["permutation_50"] = {k: v for k, v in perm50.items() if k not in ("section_interaction_delta", "section_interaction_ratio")}
    tests["permutation_100"] = {k: v for k, v in perm100.items() if k not in ("section_interaction_delta", "section_interaction_ratio")}
    tests["permutation_50"]["section_interaction_delta"] = perm50["section_interaction_delta"]
    tests["permutation_100"]["section_interaction_delta"] = perm100["section_interaction_delta"]

    with open(os.path.join(OUT, "stats.json"), "w") as fh:
        json.dump(json_clean(tests), fh, indent=2)

    plot_paired(df, os.path.join(FIG, "paired_ratio_core_margin.png"))
    plot_deltas(df, os.path.join(FIG, "section_delta_50um.png"))
    plot_curves(curve, auc_df, os.path.join(FIG, "depth_curve_50um.png"))
    # map: sample with the largest balanced core and margin
    best, best_score = None, -1
    for name, prep in primary["preps"].items():
        n_c = int((prep["dist"] >= CORE_UM).sum())
        n_m = int((prep["dist"] <= MARGIN_UM).sum())
        score = min(n_c, n_m)
        if score > best_score:
            best, best_score = name, score
    if best:
        plot_map(primary["preps"][best], best, os.path.join(FIG, "map_core_margin.png"))
        with open(os.path.join(OUT, "map_sample.txt"), "w") as fh:
            fh.write(best + "\n")
    print("wrote", OUT, flush=True)
    # console headline
    p = tests["primary_core_vs_margin_50"]["count"]
    print("PRIMARY 50 µm count", json.dumps(json_clean({
        "core": p["stratum_a"],
        "margin": p["stratum_b"],
        "interaction_delta": p["interaction_delta"],
        "interaction_ratio": p["interaction_ratio"],
        "mean_interaction_delta": p["mean_interaction_delta"],
        "mean_interaction_ratio": p["mean_interaction_ratio"],
        "patients": p["patients_core_stronger_delta"],
        "patient_p": p["patient_sign_p_one_core_stronger"],
    }), indent=2))
    print("AUC", json.dumps(json_clean(tests["auc_primary"]), indent=2))
    print("UNSTRAT", json.dumps(json_clean(tests["unstratified_50"]), indent=2)[:2000])


if __name__ == "__main__":
    main()
