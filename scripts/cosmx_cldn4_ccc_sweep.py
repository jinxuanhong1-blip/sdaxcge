#!/usr/bin/env python3
"""Sensitivity sweep for CosMx CLDN4-high vs CD8 barrier / exclusion effect sizes.

The pre-specified v1 run (40 µm, tertiles, COMMOT 50 µm) is not overwritten.
This script changes the CLDN4 split, neighborhood radius, and FOV filter, and
reports effect sizes on a probability and fold scale. COMMOT, SpatialDM, and
squidpy are then re-run at the two configurations selected by the rules below.

Selection rules (applied after the grid exists, not tuned by hand):

- Eligible: 8 slides and at least 80 FOVs.
- Immune-cold winner: most negative median slide difference in the fraction of
  tumor cells with a CD8 cell within the radius, among eligible rows with the
  high arm lower on at least 6 slides. Endpoint name ``cd8_frac``.
- Barrier winner: largest mean of the three slide-median interface log2 folds
  (CDH1, ICAM1, MIF), among eligible rows with each of those folds > 1 and the
  high arm higher on at least 6 slides. Both class means must be at least 0.05
  so a near-zero denominator cannot win.
- Gap and count filters compete in that same grid. They are not a second,
  hidden screen.

Wilcoxon p-values on a selected cell of this grid are descriptive.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

import cosmx_cldn4_spatial_ccc as ccc

ROOT = Path(__file__).resolve().parents[1]
OUT = ccc.OUT / "sweep"
CACHE = OUT / "cache"

SPLITS = {
    "tertile": (1.0 / 3.0, 2.0 / 3.0),
    "quartile": (0.25, 0.75),
    "quintile": (0.20, 0.80),
    "decile": (0.10, 0.90),
}
RADII = [15.0, 20.0, 25.0, 40.0, 50.0, 80.0, 100.0]
# min CLDN4-high, min CLDN4-low, min CD8, min mean-CLDN4 gap (log1p units)
FILTERS = [
    (10, 10, 5, 0.0),
    (10, 10, 5, 0.25),
    (10, 10, 5, 0.50),
    (20, 20, 10, 0.0),
    (20, 20, 10, 0.25),
    (40, 40, 15, 0.0),
]
LIGANDS = [
    "CDH1",
    "ICAM1",
    "MIF",
    "CD274",
    "PDCD1LG2",
    "LGALS9",
    "TGFB1",
    "CXCL9",
    "CXCL16",
    "CCL5",
    "HLA-A",
]
BARRIER_FOCUS = ["CDH1", "ICAM1", "MIF"]
MIN_FOVS = 80
MIN_LEVEL = 0.05


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def quantile_split(values: np.ndarray, lo_q: float, hi_q: float) -> np.ndarray:
    lab = np.zeros(len(values), dtype=np.int8)
    if len(values) < 6:
        return lab
    ranks = pd.Series(values).rank(method="first").to_numpy()
    a, b = np.quantile(ranks, [lo_q, hi_q])
    lab[ranks <= a] = -1
    lab[ranks >= b] = 1
    return lab


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 5 or len(b) < 5:
        return np.nan
    va = float(np.var(a, ddof=1))
    vb = float(np.var(b, ddof=1))
    sp2 = ((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2)
    if not np.isfinite(sp2) or sp2 <= 1e-12:
        return np.nan
    return float((a.mean() - b.mean()) / np.sqrt(sp2))


def _mean(x: np.ndarray) -> float:
    return float(x.mean()) if len(x) else np.nan


def _rate(x: np.ndarray) -> float:
    return float((x > 0).mean()) if len(x) else np.nan


def score_sample(df: pd.DataFrame, X: np.ndarray, gmap: dict[str, int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    lig_idx = {g: gmap[g] for g in LIGANDS + ["IFNG", "CLDN4"] if g in gmap}
    fov = df["fov"].to_numpy()
    tumor = df["is_tumor"].to_numpy()
    cd8 = df["is_cd8"].to_numpy()
    immune = df["compartment"].eq("immune").to_numpy()
    xy = df[["x_um", "y_um"]].to_numpy(float)
    cldn = df["cldn4"].to_numpy(float)
    excl_rows = []
    lig_rows = []
    for f in np.unique(fov):
        pos = np.flatnonzero(fov == f)
        tmask = tumor[pos]
        cmask = cd8[pos]
        imask = immune[pos]
        if int(tmask.sum()) < 6 or int(cmask.sum()) < 5:
            continue
        xy_f = xy[pos]
        cldn_f = cldn[pos]
        tree_cd8 = cKDTree(xy_f[cmask])
        tree_imm = cKDTree(xy_f[imask]) if imask.any() else None
        d_cd8, _ = tree_cd8.query(xy_f[tmask], k=1)
        d_imm = tree_imm.query(xy_f[tmask], k=1)[0] if tree_imm is not None else np.full(int(tmask.sum()), np.nan)
        expr = {g: X[pos, j] for g, j in lig_idx.items()}
        tumor_local = np.flatnonzero(tmask)
        cd8_local = np.flatnonzero(cmask)
        for split, (lo_q, hi_q) in SPLITS.items():
            labs = np.zeros(len(pos), dtype=np.int8)
            labs[tmask] = quantile_split(cldn_f[tmask], lo_q, hi_q)
            hi = labs == 1
            lo = labs == -1
            n_hi = int(hi.sum())
            n_lo = int(lo.sum())
            n_cd = int(cmask.sum())
            if n_hi < 5 or n_lo < 5:
                continue
            gap = float(cldn_f[hi].mean() - cldn_f[lo].mean())
            # Map tumor-local nearest distances back onto hi/lo masks.
            hi_t = hi[tmask]
            lo_t = lo[tmask]
            nn_hi = float(d_cd8[hi_t].mean()) if hi_t.any() else np.nan
            nn_lo = float(d_cd8[lo_t].mean()) if lo_t.any() else np.nan
            # IFNG on CD8 near each tumor class: distance from CD8 to that class.
            d_to_hi = d_to_lo = None
            if hi.any() and cmask.any():
                d_to_hi, _ = cKDTree(xy_f[hi]).query(xy_f[cmask], k=1)
            if lo.any() and cmask.any():
                d_to_lo, _ = cKDTree(xy_f[lo]).query(xy_f[cmask], k=1)
            ifng = expr.get("IFNG")
            rank_cache = {}
            for g in LIGANDS:
                if g not in expr:
                    continue
                vals = expr[g][tmask]
                rank_cache[g] = rankdata(vals, method="average") / max(len(vals), 1)
            for radius in RADII:
                frac_cd8_hi = float((d_cd8[hi_t] <= radius).mean())
                frac_cd8_lo = float((d_cd8[lo_t] <= radius).mean())
                frac_imm_hi = float((d_imm[hi_t] <= radius).mean()) if np.isfinite(d_imm).any() else np.nan
                frac_imm_lo = float((d_imm[lo_t] <= radius).mean()) if np.isfinite(d_imm).any() else np.nan
                excl_rows.append(
                    {
                        "sample": df["sample"].iloc[0],
                        "patient": df["patient"].iloc[0],
                        "fov": int(f),
                        "split": split,
                        "radius": radius,
                        "n_hi": n_hi,
                        "n_lo": n_lo,
                        "n_cd8": n_cd,
                        "gap": gap,
                        "frac_cd8_hi": frac_cd8_hi,
                        "frac_cd8_lo": frac_cd8_lo,
                        "frac_imm_hi": frac_imm_hi,
                        "frac_imm_lo": frac_imm_lo,
                        "nn_hi": nn_hi,
                        "nn_lo": nn_lo,
                    }
                )
                hi_if = hi.copy()
                lo_if = lo.copy()
                hi_if[tmask] = hi[tmask] & (d_cd8 <= radius)
                lo_if[tmask] = lo[tmask] & (d_cd8 <= radius)
                for g in LIGANDS:
                    if g not in expr:
                        continue
                    vec = expr[g]
                    a = vec[hi]
                    b = vec[lo]
                    ai = vec[hi_if]
                    bi = vec[lo_if]
                    rh = rank_cache[g][hi[tmask]]
                    rl = rank_cache[g][lo[tmask]]
                    lig_rows.append(
                        {
                            "sample": df["sample"].iloc[0],
                            "patient": df["patient"].iloc[0],
                            "fov": int(f),
                            "split": split,
                            "radius": radius,
                            "ligand": g,
                            "n_hi": n_hi,
                            "n_lo": n_lo,
                            "n_cd8": n_cd,
                            "gap": gap,
                            "mean_hi": _mean(a),
                            "mean_lo": _mean(b),
                            "det_hi": _rate(a),
                            "det_lo": _rate(b),
                            "rank_hi": _mean(rh),
                            "rank_lo": _mean(rl),
                            "if_n_hi": int(hi_if.sum()),
                            "if_n_lo": int(lo_if.sum()),
                            "if_mean_hi": _mean(ai),
                            "if_mean_lo": _mean(bi),
                            "if_det_hi": _rate(ai),
                            "if_det_lo": _rate(bi),
                            "if_d": cohens_d(ai, bi),
                            "class_d": cohens_d(a, b),
                        }
                    )
                if ifng is not None and d_to_hi is not None and d_to_lo is not None:
                    a = ifng[cmask][d_to_hi <= radius]
                    b = ifng[cmask][d_to_lo <= radius]
                    lig_rows.append(
                        {
                            "sample": df["sample"].iloc[0],
                            "patient": df["patient"].iloc[0],
                            "fov": int(f),
                            "split": split,
                            "radius": radius,
                            "ligand": "IFNG_on_cd8",
                            "n_hi": n_hi,
                            "n_lo": n_lo,
                            "n_cd8": n_cd,
                            "gap": gap,
                            "mean_hi": _mean(a),
                            "mean_lo": _mean(b),
                            "det_hi": _rate(a),
                            "det_lo": _rate(b),
                            "rank_hi": np.nan,
                            "rank_lo": np.nan,
                            "if_n_hi": int(len(a)),
                            "if_n_lo": int(len(b)),
                            "if_mean_hi": _mean(a),
                            "if_mean_lo": _mean(b),
                            "if_det_hi": _rate(a),
                            "if_det_lo": _rate(b),
                            "if_d": cohens_d(a, b),
                            "class_d": cohens_d(a, b),
                        }
                    )
    return pd.DataFrame(excl_rows), pd.DataFrame(lig_rows)


def cache_sample(sample: str, df: pd.DataFrame, X: np.ndarray, genes: list[str]) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "patient", "fov", "x_um", "y_um", "is_tumor", "is_cd8", "cldn4", "compartment"]
    df.loc[:, cols].to_pickle(CACHE / f"{sample}.df.pkl")
    np.savez_compressed(CACHE / f"{sample}.X.npz", X=X, genes=np.array(genes))


def load_cache(sample: str) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    df = pd.read_pickle(CACHE / f"{sample}.df.pkl")
    blob = np.load(CACHE / f"{sample}.X.npz", allow_pickle=False)
    genes = [str(g) for g in blob["genes"].tolist()]
    return df, blob["X"], genes


def assign_split(df: pd.DataFrame, lo_q: float, hi_q: float) -> pd.DataFrame:
    df = df.copy()
    n = len(df)
    tert = np.zeros(n, dtype=int)
    nn_cd8 = np.full(n, np.nan)
    nn_tumor = np.full(n, np.nan)
    nn_hi = np.full(n, np.nan)
    nn_lo = np.full(n, np.nan)
    fov = df["fov"].to_numpy()
    tumor = df["is_tumor"].to_numpy()
    cd8 = df["is_cd8"].to_numpy()
    xy = df[["x_um", "y_um"]].to_numpy(float)
    cldn = df["cldn4"].to_numpy(float)
    for f in np.unique(fov):
        pos = np.flatnonzero(fov == f)
        tmask = tumor[pos]
        cmask = cd8[pos]
        if tmask.sum() >= 6:
            tert[pos[tmask]] = quantile_split(cldn[pos][tmask], lo_q, hi_q)
        hi = tmask & (tert[pos] == 1)
        lo = tmask & (tert[pos] == -1)
        if cmask.any() and tmask.any():
            d_t, _ = cKDTree(xy[pos][cmask]).query(xy[pos][tmask], k=1)
            nn_cd8[pos[tmask]] = d_t
            d_c, _ = cKDTree(xy[pos][tmask]).query(xy[pos][cmask], k=1)
            nn_tumor[pos[cmask]] = d_c
        if cmask.any() and hi.any():
            d_c, _ = cKDTree(xy[pos][hi]).query(xy[pos][cmask], k=1)
            nn_hi[pos[cmask]] = d_c
        if cmask.any() and lo.any():
            d_c, _ = cKDTree(xy[pos][lo]).query(xy[pos][cmask], k=1)
            nn_lo[pos[cmask]] = d_c
    df["cldn4_tertile"] = tert
    df["nn_cd8_um"] = nn_cd8
    df["nn_tumor_um"] = nn_tumor
    df["nn_tumor_hi_um"] = nn_hi
    df["nn_tumor_lo_um"] = nn_lo
    return df


def _wilcox(slide_delta: pd.Series) -> dict:
    w = ccc.wilcoxon_vec(slide_delta.to_numpy(float))
    return w


def _patient_signs(slide: pd.DataFrame, col: str) -> tuple[int, int, int]:
    pat = slide.groupby("patient")[col].mean()
    return int(len(pat)), int((pat > 0).sum()), int((pat < 0).sum())


def summarize_grid(excl: pd.DataFrame, lig: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add_row(split, radius, filt, endpoint, ligand, slide, hi, lo, extra=None):
        min_hi, min_lo, min_cd, min_gap = filt
        delta = slide[hi] - slide[lo]
        w = _wilcox(delta)
        n_pat, n_pat_pos, n_pat_neg = _patient_signs(slide.assign(delta=delta.to_numpy()), "delta")
        med_hi = float(slide[hi].median())
        med_lo = float(slide[lo].median())
        fold = med_hi / med_lo if med_lo > 0 and np.isfinite(med_hi) else np.nan
        rec = {
            "split": split,
            "radius_um": radius,
            "min_hi": min_hi,
            "min_lo": min_lo,
            "min_cd8": min_cd,
            "min_gap": min_gap,
            "endpoint": endpoint,
            "ligand": ligand,
            "n_fov": int(extra or 0),
            "n_slides": int(w["n"]),
            "median_delta": float(delta.median()) if len(delta) else np.nan,
            "median_hi": med_hi,
            "median_lo": med_lo,
            "fold": float(fold) if np.isfinite(fold) else np.nan,
            "n_pos": int(w["n_pos"]),
            "n_neg": int(w["n_neg"]),
            "p": w["p"],
            "n_patients": n_pat,
            "n_patients_pos": n_pat_pos,
            "n_patients_neg": n_pat_neg,
        }
        rows.append(rec)

    for filt in FILTERS:
        min_hi, min_lo, min_cd, min_gap = filt
        ex = excl[(excl["n_hi"] >= min_hi) & (excl["n_lo"] >= min_lo) & (excl["n_cd8"] >= min_cd) & (excl["gap"] >= min_gap)]
        if ex.empty:
            continue
        for (split, radius), g in ex.groupby(["split", "radius"]):
            slide = g.groupby(["sample", "patient"], as_index=False)[
                ["frac_cd8_hi", "frac_cd8_lo", "frac_imm_hi", "frac_imm_lo", "nn_hi", "nn_lo"]
            ].mean()
            n_fov = int(len(g))
            add_row(split, radius, filt, "cd8_frac", "CD8", slide, "frac_cd8_hi", "frac_cd8_lo", n_fov)
            add_row(split, radius, filt, "immune_frac", "immune", slide, "frac_imm_hi", "frac_imm_lo", n_fov)
            add_row(split, radius, filt, "nn_um", "CD8", slide, "nn_hi", "nn_lo", n_fov)
        lg = lig[(lig["n_hi"] >= min_hi) & (lig["n_lo"] >= min_lo) & (lig["n_cd8"] >= min_cd) & (lig["gap"] >= min_gap)]
        if lg.empty:
            continue
        for (split, radius, ligand), g in lg.groupby(["split", "radius", "ligand"]):
            use = g
            slide = use.groupby(["sample", "patient"], as_index=False)[
                ["mean_hi", "mean_lo", "det_hi", "det_lo", "rank_hi", "rank_lo"]
            ].mean()
            add_row(split, radius, filt, "class_logmean", ligand, slide, "mean_hi", "mean_lo", len(use))
            add_row(split, radius, filt, "class_detection", ligand, slide, "det_hi", "det_lo", len(use))
            add_row(split, radius, filt, "class_rank", ligand, slide, "rank_hi", "rank_lo", len(use))
            iface = g[(g["if_n_hi"] >= 5) & (g["if_n_lo"] >= 5)]
            if iface.empty:
                continue
            slide_i = iface.groupby(["sample", "patient"], as_index=False)[
                ["if_mean_hi", "if_mean_lo", "if_det_hi", "if_det_lo", "if_d"]
            ].mean()
            add_row(split, radius, filt, "interface_logmean", ligand, slide_i, "if_mean_hi", "if_mean_lo", len(iface))
            add_row(split, radius, filt, "interface_detection", ligand, slide_i, "if_det_hi", "if_det_lo", len(iface))
            # Cohen's d is already a difference; store it as delta with hi=d, lo=0 so fold stays blank.
            dslide = iface.groupby(["sample", "patient"], as_index=False)["if_d"].median()
            fake = dslide.rename(columns={"if_d": "d"}).assign(zero=0.0)
            add_row(split, radius, filt, "interface_cohens_d", ligand, fake, "d", "zero", len(iface))
    return pd.DataFrame(rows)


def eligible(tab: pd.DataFrame) -> pd.DataFrame:
    return tab[(tab["n_slides"] == 8) & (tab["n_fov"] >= MIN_FOVS)].copy()


def choose_exclusion(tab: pd.DataFrame) -> pd.Series | None:
    sub = eligible(tab)
    sub = sub[(sub["endpoint"] == "cd8_frac") & (sub["n_neg"] >= 6)]
    if sub.empty:
        return None
    sub = sub.sort_values(["median_delta", "p"], ascending=[True, True])
    return sub.iloc[0]


def choose_barrier(tab: pd.DataFrame) -> pd.Series | None:
    sub = eligible(tab)
    sub = sub[(sub["endpoint"] == "interface_logmean") & (sub["ligand"].isin(BARRIER_FOCUS))]
    sub = sub[(sub["median_hi"] >= MIN_LEVEL) & (sub["median_lo"] >= MIN_LEVEL) & (sub["fold"] > 1) & (sub["n_pos"] >= 6)]
    if sub.empty:
        return None
    keys = ["split", "radius_um", "min_hi", "min_lo", "min_cd8", "min_gap"]
    parts = []
    for key, g in sub.groupby(keys):
        if set(BARRIER_FOCUS) - set(g["ligand"]):
            continue
        folds = [float(g.loc[g["ligand"] == lig, "fold"].iloc[0]) for lig in BARRIER_FOCUS]
        parts.append((*key, float(np.mean(np.log2(folds))), int(g["n_fov"].min())))
    if not parts:
        return None
    scored = pd.DataFrame(parts, columns=keys + ["mean_log2_fold", "n_fov"])
    scored = scored.sort_values(["mean_log2_fold", "n_fov"], ascending=[False, False])
    best = scored.iloc[0]
    hit = sub
    for k in keys:
        hit = hit[hit[k] == best[k]]
    # Representative row: CDH1, with the joint score attached.
    row = hit[hit["ligand"] == "CDH1"].iloc[0].copy()
    row["mean_log2_fold"] = float(best["mean_log2_fold"])
    return row


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{float(p):.3f}"


def fmt_num(x, nd=3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and ax < 1e-3:
        return f"{float(x):.2e}"
    return f"{float(x):.{nd}f}"


def config_phrase(r: pd.Series) -> str:
    return (
        f"split={r['split']}, radius={float(r['radius_um']):.0f} µm, "
        f"FOV filter ≥{int(r['min_hi'])}/{int(r['min_lo'])}/{int(r['min_cd8'])} "
        f"(high/low/CD8), CLDN4 gap ≥{float(r['min_gap']):.2f}"
    )


def top_table(tab: pd.DataFrame, endpoint: str, ligand: str | None, thesis_neg: bool, n=12) -> pd.DataFrame:
    sub = tab[tab["endpoint"] == endpoint].copy()
    if ligand is not None:
        sub = sub[sub["ligand"] == ligand]
    sub = eligible(sub)
    if sub.empty:
        return sub
    sub = sub.sort_values("median_delta", ascending=thesis_neg)
    return sub.head(n)


def write_sweep_md(path: Path, tab: pd.DataFrame, excl_win: pd.Series | None, bar_win: pd.Series | None, ccc_note: str) -> None:
    lines = []
    lines.append("# CosMx CLDN4 sweep — effect sizes for immune-cold and barrier ligands")
    lines.append("")
    lines.append("This file is a parameter sweep. The pre-specified v1 analysis (tertiles, 40 µm proximity, SpatialDM `l=30` µm, COMMOT `dis_thr=50` µm) is unchanged and is still the primary CCC report in `RESULTS.md`. Nothing here replaces the locked contact odds ratio or the Ripley g(r) result.")
    lines.append("")
    lines.append("Effect sizes below are computed from the same QC cells, epithelial tumor gate, and CD8 gate as v1. The CLDN4 split, the radius, and the FOV filter change. A negative CD8-fraction delta means CLDN4-high tumor cells less often have a CD8 cell inside the radius (immune-cold). A ligand fold above 1 means higher expression on the CLDN4-high arm.")
    lines.append("")
    lines.append("Wilcoxon tests are on the 8 slide means. p-values on a cell that was chosen from this grid are descriptive of the sweep.")
    lines.append("")
    lines.append("## Selection rules")
    lines.append("")
    lines.append("- Eligible: 8 slides and at least 80 FOVs.")
    lines.append("- Immune-cold: most negative median slide Δ in the fraction of tumor cells with a CD8 neighbor within the radius, with the high arm lower on at least 6 slides.")
    lines.append("- Barrier ligands: largest mean log2 fold of interface expression for CDH1, ICAM1, and MIF together, each fold > 1, each high arm higher on at least 6 slides, and both medians at least 0.05 so a tiny denominator cannot win.")
    lines.append("- CDH1 is also one of the epithelial markers used to call tumor. ICAM1 and MIF are not.")
    lines.append("")

    def block(title, r, interpretation):
        lines.append(f"## {title}")
        lines.append("")
        if r is None:
            lines.append("No eligible configuration met the rule.")
            lines.append("")
            return
        lines.append(interpretation)
        lines.append("")
        lines.append(f"- Setting: {config_phrase(r)}")
        lines.append(f"- FOVs used: {int(r['n_fov'])}")
        lines.append(
            f"- Levels: high {fmt_num(r['median_hi'])}, low {fmt_num(r['median_lo'])}, "
            f"median slide Δ {fmt_num(r['median_delta'])}, fold {fmt_num(r['fold'], 2)}"
        )
        lines.append(
            f"- Slides: {int(r['n_pos'])}/8 Δ>0, {int(r['n_neg'])}/8 Δ<0, Wilcoxon p={fmt_p(r['p'])}; "
            f"tissues {int(r['n_patients_pos'])}/{int(r['n_patients'])} Δ>0 and {int(r['n_patients_neg'])}/{int(r['n_patients'])} Δ<0"
        )
        lines.append("")

    if excl_win is not None:
        pp = 100 * float(excl_win["median_delta"])
        block(
            "Strongest immune-cold panel",
            excl_win,
            f"CLDN4-high tumor cells are less often within the radius of a CD8 cell than CLDN4-low tumor cells "
            f"in this setting (median slide difference {pp:.1f} percentage points).",
        )
    else:
        block("Strongest immune-cold panel", None, "")

    if bar_win is not None:
        focus = tab[
            (tab["endpoint"] == "interface_logmean")
            & (tab["split"] == bar_win["split"])
            & (np.isclose(tab["radius_um"], float(bar_win["radius_um"])))
            & (tab["min_hi"] == int(bar_win["min_hi"]))
            & (tab["min_lo"] == int(bar_win["min_lo"]))
            & (tab["min_cd8"] == int(bar_win["min_cd8"]))
            & (np.isclose(tab["min_gap"], float(bar_win["min_gap"])))
            & (tab["ligand"].isin(BARRIER_FOCUS))
        ]
        bits = []
        for lig in BARRIER_FOCUS:
            hit = focus[focus["ligand"] == lig]
            if hit.empty:
                continue
            h = hit.iloc[0]
            bits.append(
                f"{lig} high {fmt_num(h['median_hi'])} vs low {fmt_num(h['median_lo'])} "
                f"(fold {fmt_num(h['fold'], 2)}, Δ {fmt_num(h['median_delta'])}, "
                f"{int(h['n_pos'])}/8 slides Δ>0, p={fmt_p(h['p'])}, "
                f"tissues {int(h['n_patients_pos'])}/{int(h['n_patients'])})"
            )
        block(
            "Strongest barrier-ligand panel",
            bar_win,
            "Interface means are log1p median-library expression on tumor cells that have a CD8 cell within the radius. "
            + " ".join(bits),
        )
        dsub = tab[
            (tab["endpoint"] == "interface_cohens_d")
            & (tab["split"] == bar_win["split"])
            & (np.isclose(tab["radius_um"].astype(float), float(bar_win["radius_um"])))
            & (tab["min_hi"] == int(bar_win["min_hi"]))
            & (np.isclose(tab["min_gap"].astype(float), float(bar_win["min_gap"])))
            & (tab["ligand"].isin(BARRIER_FOCUS))
        ]
        if len(dsub):
            dbits = [
                f"{r.ligand} d={fmt_num(r.median_delta, 2)} ({int(r.n_pos)}/8, p={fmt_p(r.p)})"
                for r in dsub.sort_values("ligand").itertuples(index=False)
            ]
            lines.append(
                "Cell-level Cohen's d on the same interface cells (median across slides of the per-FOV d) is a standardized effect, not a fold. "
                + "; ".join(dbits)
                + ". These d values are modest. The fold is a shift in the mean, and the two cell distributions still overlap."
            )
            lines.append("")
        det = tab[
            (tab["endpoint"] == "interface_detection")
            & (tab["split"] == bar_win["split"])
            & (np.isclose(tab["radius_um"].astype(float), float(bar_win["radius_um"])))
            & (tab["min_hi"] == int(bar_win["min_hi"]))
            & (np.isclose(tab["min_gap"].astype(float), float(bar_win["min_gap"])))
            & (tab["ligand"].isin(BARRIER_FOCUS))
        ]
        if len(det):
            dbits = []
            for r in det.sort_values("ligand").itertuples(index=False):
                dbits.append(
                    f"{r.ligand} {100*float(r.median_hi):.1f}% vs {100*float(r.median_lo):.1f}% "
                    f"({100*float(r.median_delta):+.1f} pp, {int(r.n_pos)}/8, p={fmt_p(r.p)}, "
                    f"tissues {int(r.n_patients_pos)}/{int(r.n_patients)})"
                )
            lines.append("Detection rate on those same interface cells: " + "; ".join(dbits) + ".")
            lines.append("")
    else:
        block("Strongest barrier-ligand panel", None, "")

    lines.append("## Immune-compartment neighborhood")
    lines.append("")
    lines.append("The CD8-fraction rule above is the neighborhood of the CCC receiver. The broader immune-compartment call (the same epithelial / immune / stromal argmax used to define tumor) is a larger immune-cold contrast. It is reported here and is not a ligand–receptor score.")
    lines.append("")
    imm = eligible(tab)
    imm = imm[(imm["endpoint"] == "immune_frac") & (imm["n_neg"] >= 6)].sort_values("median_delta")
    if len(imm):
        r = imm.iloc[0]
        pp = 100 * float(r["median_delta"])
        lines.append(
            f"Strongest eligible row: {config_phrase(r)}. "
            f"High {fmt_num(r['median_hi'])} vs low {fmt_num(r['median_lo'])} "
            f"({pp:.1f} percentage points, fold {fmt_num(r['fold'], 2)}), "
            f"{int(r['n_neg'])}/8 slides Δ<0, Wilcoxon p={fmt_p(r['p'])}, "
            f"tissues {int(r['n_patients_neg'])}/{int(r['n_patients'])} Δ<0, {int(r['n_fov'])} FOVs."
        )
        base_imm = imm[(imm["min_hi"] == 10) & (imm["min_cd8"] == 5) & (np.isclose(imm["min_gap"], 0.0))]
        if len(base_imm):
            b = base_imm.iloc[0]
            lines.append("")
            lines.append(
                f"Without a CLDN4-gap filter and with the loosest cell floor (≥10/10/5), the strongest radius/split is "
                f"{b['split']} at {float(b['radius_um']):.0f} µm: high {fmt_num(b['median_hi'])} vs low {fmt_num(b['median_lo'])} "
                f"({100*float(b['median_delta']):.1f} percentage points), {int(b['n_neg'])}/8 slides, p={fmt_p(b['p'])}, "
                f"tissues {int(b['n_patients_neg'])}/{int(b['n_patients'])}, {int(b['n_fov'])} FOVs."
            )
        lines.append("")
    cd8_base = tab[
        (tab["endpoint"] == "cd8_frac")
        & (tab["split"] == "tertile")
        & (np.isclose(tab["radius_um"], 40.0))
        & (tab["min_hi"] == 10)
        & (tab["min_cd8"] == 5)
        & (np.isclose(tab["min_gap"], 0.0))
    ]
    if len(cd8_base):
        b = cd8_base.iloc[0]
        lines.append(
            f"The same tertile / 40 µm CD8 contrast with no gap filter is high {fmt_num(b['median_hi'])} vs low {fmt_num(b['median_lo'])} "
            f"({100*float(b['median_delta']):.1f} percentage points), {int(b['n_neg'])}/8, p={fmt_p(b['p'])}, "
            f"tissues {int(b['n_patients_neg'])}/{int(b['n_patients'])}, {int(b['n_fov'])} FOVs. "
            f"The gap filter changes this by a fraction of a percentage point."
        )
        lines.append("")

    lines.append("## Same settings, the other ligands")
    lines.append("")
    lines.append("These rows use the barrier-winning setting when one exists, so ligands that do not support the thesis stay visible. Interface log-mean, CLDN4-high minus CLDN4-low.")
    lines.append("")
    if bar_win is not None:
        others = tab[
            (tab["endpoint"] == "interface_logmean")
            & (tab["split"] == bar_win["split"])
            & (np.isclose(tab["radius_um"].astype(float), float(bar_win["radius_um"])))
            & (tab["min_hi"] == int(bar_win["min_hi"]))
            & (np.isclose(tab["min_gap"].astype(float), float(bar_win["min_gap"])))
            & (tab["min_cd8"] == int(bar_win["min_cd8"]))
        ].sort_values("ligand")
        for r in others.itertuples(index=False):
            lines.append(
                f"- {r.ligand}: high {fmt_num(r.median_hi)} vs low {fmt_num(r.median_lo)}, "
                f"fold {fmt_num(r.fold, 2)}, Δ {fmt_num(r.median_delta)}, "
                f"{int(r.n_pos)}/8 Δ>0, p={fmt_p(r.p)}, "
                f"tissues {int(r.n_patients_pos)}/{int(r.n_patients)} Δ>0"
            )
        lines.append("")
    lines.append("IFNG_on_cd8 is IFNG on CD8 cells sitting within the radius of that tumor class. A flat contrast is the not-muzzling check.")
    lines.append("")

    lines.append("## How much of the grid points the same way")
    lines.append("")
    elig = eligible(tab)
    cd8 = elig[elig["endpoint"] == "cd8_frac"]
    if len(cd8):
        cold = cd8[cd8["median_delta"] < 0]
        sig = cold[cold["p"] <= 0.05]
        lines.append(
            f"Eligible CD8-fraction rows: {len(cd8)}. "
            f"Median Δ < 0 (CLDN4-high has fewer nearby CD8): {len(cold)}. "
            f"Of those, Wilcoxon p≤0.05: {len(sig)}."
        )
    base = tab[
        (tab["endpoint"] == "cd8_frac")
        & (tab["min_hi"] == 10)
        & (tab["min_cd8"] == 5)
        & (np.isclose(tab["min_gap"], 0.0))
    ]
    if len(base):
        lines.append("")
        lines.append("Base filter only (≥10/10/5, no CLDN4-gap filter), median slide Δ in CD8 fraction (negative = immune-cold):")
        lines.append("")
        pivot = base.pivot_table(index="split", columns="radius_um", values="median_delta")
        lines.append(_md_table(pivot))
        lines.append("")
    lines.append("Full grid: `sweep_table.tsv`. One row per split, radius, FOV filter, endpoint, and ligand.")
    lines.append("")
    lines.append("## CCC tools at the selected settings")
    lines.append("")
    lines.append(ccc_note)
    lines.append("")
    lines.append("SpatialDM length scale is set so the RBF weight exp(−d² / 2l²) falls through 0.15 near the selected radius (`l = radius / sqrt(−2 ln 0.15)`). COMMOT `dis_thr` is that same radius. Squidpy `gr.ligrec` is restricted to the same radius. Caps, permutation count (1,000), and the CellChat pair list match v1.")
    lines.append("")
    lines.append("A Moran I delta stays on a correlation scale. A COMMOT per-sender delta stays on a transport-mass scale and will look small even when the high/low fold is not. The fraction of senders with any positive transport is the transport score on a probability scale.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def _md_table(pivot: pd.DataFrame) -> str:
    cols = list(pivot.columns)
    header = "| split | " + " | ".join(f"{float(c):.0f} µm" for c in cols) + " |"
    sep = "|---|" + "|".join(["---:"] * len(cols)) + "|"
    body = []
    for idx, row in pivot.iterrows():
        cells = " | ".join(fmt_num(row[c], 3) for c in cols)
        body.append(f"| {idx} | {cells} |")
    return "\n".join([header, sep, *body])


def run_fast() -> tuple[pd.DataFrame, pd.DataFrame]:
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = ccc.pair_frame()
    panel = ccc.read_panel(ccc.find_file(ccc.SAMPLES[0], "exprMat", ccc.DATA))
    keep = [g for g in ccc.gene_universe(pairs) if g in set(panel)]
    excl_parts = []
    lig_parts = []
    for sample in ccc.SAMPLES:
        stamp = OUT / f"{sample}.excl.tsv.gz"
        lstamp = OUT / f"{sample}.lig.tsv.gz"
        if stamp.exists() and lstamp.exists() and (CACHE / f"{sample}.df.pkl").exists():
            log(f"{sample}: resume sweep rows")
            excl_parts.append(pd.read_csv(stamp, sep="\t"))
            lig_parts.append(pd.read_csv(lstamp, sep="\t"))
            continue
        df, X, genes = ccc.load_sample(sample, ccc.DATA, keep, panel)
        cache_sample(sample, df, X, genes)
        gmap = {g: i for i, g in enumerate(genes)}
        log(f"{sample}: scoring splits and radii")
        ex, lg = score_sample(df, X, gmap)
        ex.to_csv(stamp, sep="\t", index=False)
        lg.to_csv(lstamp, sep="\t", index=False)
        excl_parts.append(ex)
        lig_parts.append(lg)
        del df, X
    return pd.concat(excl_parts, ignore_index=True), pd.concat(lig_parts, ignore_index=True)


def ccc_summary(fov: pd.DataFrame) -> pd.DataFrame:
    fov = ccc.add_deltas(fov)
    rows = []
    for pid, g in fov.groupby("pair_id"):
        slide = g.groupby(["sample", "patient"], as_index=False).mean(numeric_only=True)
        meta = g.iloc[0]
        for ep, hi, lo in (
            ("expr", "expr_hi", "expr_lo"),
            ("sdm_prox", "I_hi_prox", "I_lo_prox"),
            ("sdm_geom", "I_hi_geom", "I_lo_geom"),
            ("commot", "commot_per_sender_hi", "commot_per_sender_lo"),
            ("commot_supply", "commot_per_supply_hi", "commot_per_supply_lo"),
            ("commot_frac_pos", "commot_frac_pos_hi", "commot_frac_pos_lo"),
        ):
            if hi not in slide.columns:
                continue
            delta = slide[hi] - slide[lo]
            w = ccc.wilcoxon_vec(delta.to_numpy(float))
            med_hi = float(slide[hi].median())
            med_lo = float(slide[lo].median())
            fold = med_hi / med_lo if np.isfinite(med_hi) and np.isfinite(med_lo) and med_lo > 0 else np.nan
            n_pat, n_pos, n_neg = _patient_signs(slide.assign(delta=delta.to_numpy()), "delta")
            rows.append(
                {
                    "pair_id": pid,
                    "ligand": meta["ligand"],
                    "family": meta["family"],
                    "endpoint": ep,
                    "n_fov": int(g.groupby(["sample", "fov"]).ngroups),
                    "n_slides": int(w["n"]),
                    "median_hi": med_hi,
                    "median_lo": med_lo,
                    "median_delta": float(delta.median()) if len(delta) else np.nan,
                    "fold": fold,
                    "n_pos": int(w["n_pos"]),
                    "n_neg": int(w["n_neg"]),
                    "p": w["p"],
                    "n_patients": n_pat,
                    "n_patients_pos": n_pos,
                    "n_patients_neg": n_neg,
                }
            )
    return pd.DataFrame(rows)


def filter_fovs(df: pd.DataFrame, X: np.ndarray, min_hi: int, min_lo: int, min_cd: int, min_gap: float) -> tuple[pd.DataFrame, np.ndarray]:
    keep = []
    for f, g in df.groupby("fov"):
        hi = g["is_tumor"] & g["cldn4_tertile"].eq(1)
        lo = g["is_tumor"] & g["cldn4_tertile"].eq(-1)
        if int(hi.sum()) < min_hi or int(lo.sum()) < min_lo or int(g["is_cd8"].sum()) < min_cd:
            continue
        gap = float(g.loc[hi, "cldn4"].mean() - g.loc[lo, "cldn4"].mean())
        if gap < min_gap:
            continue
        keep.append(f)
    mask = df["fov"].isin(keep).to_numpy()
    return df.loc[mask].reset_index(drop=True), X[mask]


def run_one_ccc(name: str, split: str, radius: float, filt: tuple, pairs: pd.DataFrame) -> pd.DataFrame:
    lo_q, hi_q = SPLITS[split]
    min_hi, min_lo, min_cd, min_gap = filt
    # RBF weight exp(-d^2 / (2 l^2)) = 0.15 at d = radius.
    sdm_l = float(radius) / np.sqrt(-2.0 * np.log(0.15))
    saved = (ccc.PROX_UM, ccc.COMMOT_DIS_UM, ccc.SDM_L_UM)
    ccc.PROX_UM = float(radius)
    ccc.COMMOT_DIS_UM = float(radius)
    ccc.SDM_L_UM = sdm_l
    log(f"CCC {name}: {split} radius={radius} l={sdm_l:.2f}")
    fov_parts = []
    sq_parts = []
    try:
        for sample in ccc.SAMPLES:
            ck = OUT / "cache" / f"{name}.{sample}.fov.tsv.gz"
            sk = OUT / "cache" / f"{name}.{sample}.squidpy.tsv"
            if ck.exists() and sk.exists():
                log(f"{name} {sample}: resume")
                fov_parts.append(pd.read_csv(ck, sep="\t"))
                sq_parts.append(pd.read_csv(sk, sep="\t"))
                continue
            df, X, genes = load_cache(sample)
            gmap = {g: i for i, g in enumerate(genes)}
            df = assign_split(df, lo_q, hi_q)
            df, X = filter_fovs(df, X, min_hi, min_lo, min_cd, min_gap)
            log(f"{name} {sample}: FOVs after filter={df['fov'].nunique()} cells={len(df)}")
            fov, _geo = ccc.fov_records(df, X, gmap, pairs, None, skip_commot=False)
            squid = ccc.squidpy_slide(df, X, gmap, pairs)
            fov.to_csv(ck, sep="\t", index=False)
            pd.DataFrame(squid).to_csv(sk, sep="\t", index=False)
            fov_parts.append(fov)
            sq_parts.append(pd.DataFrame(squid))
            del df, X
    finally:
        ccc.PROX_UM, ccc.COMMOT_DIS_UM, ccc.SDM_L_UM = saved
    fov = pd.concat(fov_parts, ignore_index=True) if fov_parts else pd.DataFrame()
    summary = ccc_summary(fov) if len(fov) else pd.DataFrame()
    summary.insert(0, "config", name)
    summary.insert(1, "split", split)
    summary.insert(2, "radius_um", radius)
    summary.insert(3, "sdm_l_um", sdm_l)
    squid = pd.concat(sq_parts, ignore_index=True) if sq_parts else pd.DataFrame()
    if len(squid):
        squid.to_csv(OUT / f"{name}.squidpy_slide.tsv", sep="\t", index=False)
        fam = fov[["pair_id", "ligand", "family"]].drop_duplicates() if len(fov) else pd.DataFrame(columns=["pair_id", "ligand", "family"])
        clean = []
        for pid, g in squid.groupby("pair_id"):
            slide = g.groupby(["sample", "patient"], as_index=False)[["mean_hi", "mean_lo", "delta_mean"]].mean()
            w = ccc.wilcoxon_vec(slide["delta_mean"].to_numpy(float))
            med_hi = float(slide["mean_hi"].median())
            med_lo = float(slide["mean_lo"].median())
            n_pat, n_pos, n_neg = _patient_signs(slide, "delta_mean")
            meta = fam[fam["pair_id"] == pid]
            clean.append(
                {
                    "config": name,
                    "split": split,
                    "radius_um": radius,
                    "sdm_l_um": sdm_l,
                    "pair_id": pid,
                    "ligand": meta["ligand"].iloc[0] if len(meta) else "",
                    "family": meta["family"].iloc[0] if len(meta) else "",
                    "endpoint": "squidpy",
                    "n_fov": np.nan,
                    "n_slides": int(w["n"]),
                    "median_hi": med_hi,
                    "median_lo": med_lo,
                    "median_delta": float(slide["delta_mean"].median()),
                    "fold": (med_hi / med_lo) if med_lo > 0 else np.nan,
                    "n_pos": int(w["n_pos"]),
                    "n_neg": int(w["n_neg"]),
                    "p": w["p"],
                    "n_patients": n_pat,
                    "n_patients_pos": n_pos,
                    "n_patients_neg": n_neg,
                }
            )
        summary = pd.concat([summary, pd.DataFrame(clean)], ignore_index=True)
    summary.to_csv(OUT / f"{name}.pair_summary.tsv", sep="\t", index=False)
    return summary


def append_ccc_section(summary: pd.DataFrame) -> str:
    if summary is None or summary.empty:
        return "CCC re-run did not produce a table."
    lines = []
    for config, g in summary.groupby("config", sort=False):
        split = g["split"].iloc[0]
        radius = float(g["radius_um"].iloc[0])
        lines.append(f"### {config} ({split}, {radius:.0f} µm)")
        lines.append("")
        show = g[g["family"].isin(["barrier", "mhc_control", "effector_reverse"])]
        focus_pairs = {
            "CDH1|CDH1",
            "CDH1|ITGA2_ITGB1",
            "ICAM1|ITGAL",
            "ICAM1|ITGAL_ITGB2",
            "MIF|CD74_CD44",
            "LGALS9|HAVCR2",
            "CD274|PDCD1",
            "PDCD1LG2|PDCD1",
            "TGFB1|TGFBR1_TGFBR2",
            "HLA-A|CD8A",
            "IFNG|IFNGR1_IFNGR2",
        }
        show = g[g["pair_id"].isin(focus_pairs) & g["endpoint"].isin(["expr", "commot", "commot_frac_pos", "commot_supply", "sdm_prox", "squidpy"])]
        show = show.sort_values(["pair_id", "endpoint"])
        for r in show.itertuples(index=False):
            lines.append(
                f"- {r.pair_id} {r.endpoint}: high {fmt_num(r.median_hi)} vs low {fmt_num(r.median_lo)}, "
                f"fold {fmt_num(r.fold, 2)}, Δ {fmt_num(r.median_delta)}, "
                f"{int(r.n_pos)}/{int(r.n_slides)} slides Δ>0, p={fmt_p(r.p)}, "
                f"tissues {int(r.n_patients_pos)}/{int(r.n_patients)} Δ>0"
            )
        lines.append("")
    return "\n".join(lines)


def plot_base_heatmap(tab: pd.DataFrame, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = tab[
        (tab["endpoint"] == "cd8_frac")
        & (tab["min_hi"] == 10)
        & (tab["min_cd8"] == 5)
        & (np.isclose(tab["min_gap"], 0.0))
    ]
    if base.empty:
        return
    pivot = base.pivot_table(index="split", columns="radius_um", values="median_delta")
    pivot = pivot.reindex(index=list(SPLITS), columns=RADII)
    fig, ax = plt.subplots(figsize=(8.2, 3.4))
    im = ax.imshow(pivot.to_numpy(float), cmap="RdBu", vmin=-0.15, vmax=0.15, aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{c:.0f}" for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(list(pivot.index))
    ax.set_xlabel("Radius (µm)")
    ax.set_title("CD8 within radius, CLDN4-high minus CLDN4-low")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.to_numpy(float)[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{100 * val:.1f} pp", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-ccc", action="store_true")
    args = ap.parse_args()
    excl, lig = run_fast()
    log(f"fov rows exclusion={len(excl)} ligand={len(lig)}")
    tab = summarize_grid(excl, lig)
    tab.to_csv(OUT / "sweep_table.tsv", sep="\t", index=False)
    log(f"sweep rows {len(tab)}")
    excl_win = choose_exclusion(tab)
    bar_win = choose_barrier(tab)
    wins = []
    if excl_win is not None:
        wins.append(("immune_cold", excl_win))
    if bar_win is not None:
        # Skip a second CCC run when the barrier winner is the same split and radius.
        same = (
            excl_win is not None
            and excl_win["split"] == bar_win["split"]
            and np.isclose(float(excl_win["radius_um"]), float(bar_win["radius_um"]))
            and int(excl_win["min_hi"]) == int(bar_win["min_hi"])
            and int(excl_win["min_lo"]) == int(bar_win["min_lo"])
            and int(excl_win["min_cd8"]) == int(bar_win["min_cd8"])
            and np.isclose(float(excl_win["min_gap"]), float(bar_win["min_gap"]))
        )
        if not same:
            wins.append(("barrier", bar_win))
    pd.DataFrame(
        [
            {
                "role": role,
                "split": r["split"],
                "radius_um": r["radius_um"],
                "min_hi": r["min_hi"],
                "min_lo": r["min_lo"],
                "min_cd8": r["min_cd8"],
                "min_gap": r["min_gap"],
                "median_delta": r["median_delta"],
                "fold": r.get("fold", np.nan),
                "p": r["p"],
                "n_fov": r["n_fov"],
                "n_pos": r["n_pos"],
                "n_neg": r["n_neg"],
            }
            for role, r in (("immune_cold", excl_win), ("barrier", bar_win))
            if r is not None
        ]
    ).to_csv(OUT / "selected_configs.tsv", sep="\t", index=False)
    plot_base_heatmap(tab, OUT / "figures" / "cd8_frac_delta_heatmap.png")
    ccc_note = "CCC re-run skipped."
    if not args.skip_ccc and wins:
        pairs = ccc.pair_frame()
        # Confirm against the bundled CellChat table once.
        pairs = ccc.confirm_pairs_in_cellchat(pairs)
        pieces = []
        seen = set()
        for role, r in wins:
            filt = (int(r["min_hi"]), int(r["min_lo"]), int(r["min_cd8"]), float(r["min_gap"]))
            key = (r["split"], float(r["radius_um"]), filt)
            if key in seen:
                continue
            seen.add(key)
            name = f"{role}_{r['split']}_{int(float(r['radius_um']))}um_g{int(round(float(r['min_gap']) * 100))}"
            pieces.append(run_one_ccc(name, r["split"], float(r["radius_um"]), filt, pairs))
        if pieces:
            all_sum = pd.concat(pieces, ignore_index=True)
            all_sum.to_csv(OUT / "ccc_selected_summary.tsv", sep="\t", index=False)
            ccc_note = append_ccc_section(all_sum)
    elif args.skip_ccc:
        ccc_note = "CCC re-run not started in this process (`--skip-ccc`)."
    write_sweep_md(OUT / "SWEEP.md", tab, excl_win, bar_win, ccc_note)
    # Root copy so the sweep is visible next to RESULTS.md
    text = (OUT / "SWEEP.md").read_text()
    (ROOT / "SWEEP.md").write_text(text)
    log(f"wrote {OUT / 'SWEEP.md'}")
    if excl_win is not None:
        log(f"immune-cold {config_phrase(excl_win)} delta={excl_win['median_delta']:.4f} p={excl_win['p']}")
    if bar_win is not None:
        log(f"barrier {config_phrase(bar_win)} fold={bar_win['fold']:.3f}")


if __name__ == "__main__":
    main()
