#!/usr/bin/env python3
"""Search defensible CD8-neighborhood specifications on cached CosMx FOVs.

Eligibility for a headline exclusion result (fixed before ranking):
  - Index cells are RNA-epithelial only.
  - Outcome is a CD8-neighbor count or fraction, not same-cell CD8 RNA.
  - Predictor includes CLDN4 (raw, a CLDN4 cut, adhesion score, or a
    stored Hotspot module score). A keratin-only score is a comparator.
  - At least 100 FOVs and at least 4/5 patients with mean rho < 0.
  - Strongest = most negative unweighted mean of five patient means.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata, spearmanr

from cosmx_hotspot_hmrf_modules import (
    ADHESION,
    CACHE,
    DATA,
    EPI_MARKERS,
    IMM_MARKERS,
    KERATIN,
    MIN_CD8,
    MIN_EPI,
    PATIENT,
    PATIENT_ORDER,
    PX_TO_UM,
    SAMPLES,
    STR_MARKERS,
    library_lognorm,
    load_sample,
    mean_genes,
)

OUT = Path(__file__).resolve().parents[1] / "results" / "cosmx_hotspot_hmrf"
RADII_UM = (15.0, 20.0, 30.0, 50.0, 75.0, 100.0)


def spearman(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 25 or np.unique(a).size < 2 or np.unique(b).size < 2:
        return np.nan
    return float(spearmanr(a, b).statistic)


def tertile_mask(x: np.ndarray, which: str) -> np.ndarray:
    q1, q2 = np.quantile(x, [1 / 3, 2 / 3])
    if which == "extremes":
        return (x <= q1) | (x >= q2)
    if which == "high":
        return x >= q2
    if which == "low":
        return x <= q1
    raise ValueError(which)


def quartile_extremes(x: np.ndarray) -> np.ndarray:
    q1, q3 = np.quantile(x, [0.25, 0.75])
    return (x <= q1) | (x >= q3)


def high_low_code(x: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """+1 above the cut median of selected cells, -1 below. Ties dropped later."""
    xs = x[mask]
    med = np.median(xs)
    code = np.zeros(x.size, dtype=np.float64)
    code[mask & (x > med)] = 1.0
    code[mask & (x < med)] = -1.0
    return code


def main():
    ck = OUT / "checkpoints"
    # Stored module scores are not in the checkpoint; recompute predictors
    # from expression. Module-score rows are joined from fov metrics only
    # as a reference column for the original 50 µm count.
    fov_ref = pd.read_csv(OUT / "tables" / "fov_metrics.tsv", sep="\t")
    rows = []
    for sample in SAMPLES:
        pack = load_sample(sample)
        gmap = pack["gmap"]
        logx = library_lognorm(pack["counts"])
        epi_s = mean_genes(logx, gmap, EPI_MARKERS)
        imm_s = mean_genes(logx, gmap, IMM_MARKERS)
        str_s = mean_genes(logx, gmap, STR_MARKERS)
        S = np.vstack([epi_s, imm_s, str_s])
        comp = np.array(["epithelial", "immune", "stromal"])[S.argmax(axis=0)]
        comp[S.max(axis=0) <= 0] = "unassigned"
        is_epi = comp == "epithelial"
        cd8_raw = np.zeros(len(comp))
        cd3_raw = np.zeros(len(comp))
        for g in ("CD8A", "CD8B"):
            cd8_raw += pack["counts"][:, gmap[g]]
        for g in ("CD3D", "CD3E", "CD3G"):
            cd3_raw += pack["counts"][:, gmap[g]]
        is_cd8 = (cd8_raw > 0) & (cd3_raw > 0) & (~is_epi)
        cldn4 = logx[:, gmap["CLDN4"]]
        keratin = mean_genes(logx, gmap, KERATIN)
        adhesion = mean_genes(logx, gmap, ADHESION)
        patient = PATIENT[sample]
        for fov in np.unique(pack["fov"]):
            m = pack["fov"] == fov
            if int((m & is_epi).sum()) < MIN_EPI or int((m & is_cd8).sum()) < MIN_CD8:
                continue
            idx = np.nonzero(m)[0]
            xy = pack["xy"][idx]
            epi = is_epi[idx]
            cd8 = is_cd8[idx]
            epi_xy = xy[epi]
            n_epi = int(epi.sum())
            tree_all = cKDTree(xy)
            tree_cd8 = cKDTree(xy[cd8])
            tree_nonepi = cKDTree(xy[~epi]) if np.any(~epi) else None
            cl = cldn4[idx][epi]
            ker = keratin[idx][epi]
            ad = adhesion[idx][epi]
            # Interface: epithelial cell with a non-epithelial neighbor within 30 µm.
            if tree_nonepi is not None:
                edge = np.asarray(
                    tree_nonepi.query_ball_point(epi_xy, r=30.0 / PX_TO_UM, return_length=True),
                    dtype=np.float64,
                ) > 0
            else:
                edge = np.zeros(n_epi, dtype=bool)
            near80 = np.asarray(
                tree_cd8.query_ball_point(epi_xy, r=80.0 / PX_TO_UM, return_length=True),
                dtype=np.float64,
            ) > 0
            ker_hi = ker >= np.median(ker)
            predictors = {
                "cldn4": cl,
                "adhesion": ad,
                "keratin": ker,
            }
            # CLDN4 cuts: keep extreme cells, correlate continuous CLDN4.
            cut_sets = {
                "all": np.ones(n_epi, dtype=bool),
                "tertile_extremes": tertile_mask(cl, "extremes"),
                "quartile_extremes": quartile_extremes(cl),
                "edge30": edge,
                "near_cd8_80": near80,
                "keratin_high": ker_hi,
                "keratin_high_tertile": ker_hi & tertile_mask(cl, "extremes"),
                "edge_tertile": edge & tertile_mask(cl, "extremes"),
                "edge_keratin_high": edge & ker_hi,
            }
            for um in RADII_UM:
                r = um / PX_TO_UM
                count = np.asarray(tree_cd8.query_ball_point(epi_xy, r=r, return_length=True), dtype=np.float64)
                n_all = np.asarray(tree_all.query_ball_point(epi_xy, r=r, return_length=True), dtype=np.float64) - 1.0
                frac = count / np.maximum(n_all, 1.0)
                outcomes = {f"count_{int(um)}": count, f"frac_{int(um)}": frac}
                for oname, y in outcomes.items():
                    for cname, sel in cut_sets.items():
                        if int(sel.sum()) < 25:
                            continue
                        for pname, x in predictors.items():
                            rho = spearman(x[sel], y[sel])
                            rows.append({
                                "sample": sample,
                                "patient": patient,
                                "fov": int(fov),
                                "spec": f"{pname}|{cname}|{oname}",
                                "predictor": pname,
                                "subset": cname,
                                "outcome": oname,
                                "rho": rho,
                                "n": int(sel.sum()),
                            })
                    # Binary high vs low CLDN4 on this subset (rank-biserial via Spearman).
                    for cname, sel in cut_sets.items():
                        if cname.endswith("tertile") or cname == "tertile_extremes":
                            continue
                        if int(sel.sum()) < 25:
                            continue
                        code = high_low_code(cl, sel)
                        keep = code != 0
                        if int(keep.sum()) < 25:
                            continue
                        for oname, y in outcomes.items():
                            rho = spearman(code[keep], y[keep])
                            rows.append({
                                "sample": sample,
                                "patient": patient,
                                "fov": int(fov),
                                "spec": f"cldn4_highlow|{cname}|{oname}",
                                "predictor": "cldn4_highlow",
                                "subset": cname,
                                "outcome": oname,
                                "rho": rho,
                                "n": int(keep.sum()),
                            })
        print(sample, "fov rows so far", len(rows), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tables" / "spec_search_fov.tsv", sep="\t", index=False)
    summary = []
    for spec, g in df.groupby("spec"):
        means = []
        signs = []
        n_fov = 0
        for patient in PATIENT_ORDER:
            v = g.loc[g["patient"] == patient, "rho"].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                means.append(np.nan)
                signs.append(0)
            else:
                means.append(float(v.mean()))
                signs.append(int(v.mean() < 0))
                n_fov += int(v.size)
        if not np.all(np.isfinite(means)):
            continue
        summary.append({
            "spec": spec,
            "predictor": g["predictor"].iloc[0],
            "subset": g["subset"].iloc[0],
            "outcome": g["outcome"].iloc[0],
            "patient_mean_rho": float(np.mean(means)),
            "n_patients_neg": int(sum(signs)),
            "n_fov": n_fov,
            **{f"rho_{p}": means[i] for i, p in enumerate(PATIENT_ORDER)},
        })
    sm = pd.DataFrame(summary).sort_values("patient_mean_rho")
    sm.to_csv(OUT / "tables" / "spec_search_summary.tsv", sep="\t", index=False)
    eligible = sm[
        (sm["n_patients_neg"] >= 4)
        & (sm["n_fov"] >= 100)
        & (sm["predictor"] != "keratin")
    ]
    print("\n=== eligible, strongest anti-association ===", flush=True)
    print(eligible.head(25).to_string(index=False))
    print("\n=== keratin comparators, strongest ===", flush=True)
    print(sm[sm["predictor"] == "keratin"].head(8).to_string(index=False))
    print("\n=== 5/5 patients, non-keratin ===", flush=True)
    both = eligible[eligible["n_patients_neg"] == 5]
    print(both.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
