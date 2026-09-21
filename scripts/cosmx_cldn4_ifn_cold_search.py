#!/usr/bin/env python3
"""Search IFN-neighborhood specifications for CLDN4-high → IFN-low / immune-cold.

Varies gene set, aggregation (mean vs sum vs count of positive neighbors),
radius, CLDN4 cutoff, and FOV vs section unit. Writes every spec. Does not
invent statistics. The broad equal-weight mean from the first pass is included
so a sign flip has to come from the data.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_cldn4_ifn_sting", "tables")

UM_PER_PX = 0.18
RADII = (20, 40, 60, 100, 150)
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
IMMUNE = {
    "B-cell", "NK", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive",
    "Treg", "mDC", "macrophage", "mast", "monocyte", "neutrophil", "pDC", "plasmablast",
}
CYTO = {"NK", "T CD8 memory", "T CD8 naive"}
CHEMO = ["IFNG", "CXCL9", "CXCL10", "CCL5"]
ISG = ["STAT1", "MX1", "IFIT1", "IFITM3", "OAS3", "IDO1"]
BROAD = [
    "STAT1", "JAK1", "JAK2", "IFNGR1", "IFNGR2", "IFNAR1", "IFNAR2",
    "IFNG", "IFNB1", "IFNA1", "IFIT1", "IFITM1", "IFITM3", "MX1",
    "OAS1", "OAS2", "OAS3", "OASL", "CXCL9", "CXCL10", "CCL5",
    "CD274", "IDO1", "BST2", "TAP1", "TAP2",
]
SINGLES = ["IFNG", "CXCL9", "CXCL10", "CCL5", "CD274", "IDO1", "STAT1", "MX1", "CXCL8", "TNF"]


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def load():
    import h5py

    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}

    def categorical(name):
        cats = _decode(f[f"obs/{name}/categories"][:])
        codes = f[f"obs/{name}/codes"][:]
        return np.asarray(cats, dtype=object)[codes]

    obs = {
        "sample": categorical("sample"),
        "cell_type": categorical("cell_type"),
        "fov": f["obs/fov"][:].astype(np.int32),
        "n_counts": f["obs/n_counts"][:].astype(np.float64),
        "xy": f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX,
    }
    counts = sparse.csr_matrix(
        (f["layers/counts/data"][:], f["layers/counts/indices"][:], f["layers/counts/indptr"][:]),
        shape=(len(obs["sample"]), len(genes)),
    )
    f.close()
    need = sorted(set(BROAD + CHEMO + ISG + SINGLES + ["CLDN4"]))
    lib = obs["n_counts"].copy()
    lib[lib <= 0] = np.nan
    loge = {}
    for g in need:
        col = np.asarray(counts.getcol(gmap[g]).todense()).ravel().astype(np.float64)
        loge[g] = np.nan_to_num(np.log1p(col / lib * 1e4))
    del counts
    return obs, loge


def adjacency(xy, radius):
    tree = cKDTree(xy)
    coo = tree.sparse_distance_matrix(tree, max_distance=radius, output_type="coo_matrix")
    n = xy.shape[0]
    mask = coo.row != coo.col
    A = sparse.csr_matrix(
        (np.ones(int(mask.sum()), dtype=np.float32), (coo.row[mask], coo.col[mask])),
        shape=(n, n),
    )
    A.sum_duplicates()
    return A


def wilcoxon(high, low):
    d = np.asarray(high, float) - np.asarray(low, float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return dict(n=0, n_neg=0, n_pos=0, median_delta=np.nan, p=np.nan)
    n_neg = int((d < 0).sum())
    n_pos = int((d > 0).sum())
    usable = d[d != 0]
    p = np.nan
    if len(usable) >= 1:
        method = "exact" if len(usable) <= 25 else "auto"
        p = float(stats.wilcoxon(usable, alternative="two-sided", method=method, zero_method="wilcox").pvalue)
    return dict(n=int(len(d)), n_neg=n_neg, n_pos=n_pos, median_delta=float(np.median(d)), p=p)


def arms_from_values(vals, mode):
    """Return high, low boolean masks over vals (1d)."""
    high = np.zeros(len(vals), dtype=bool)
    low = np.zeros(len(vals), dtype=bool)
    if len(vals) < 20:
        return high, low
    if mode == "q4_q1":
        q1, q3 = np.quantile(vals, [0.25, 0.75])
        high = vals >= q3
        low = vals <= q1
        high &= ~low
    elif mode == "median":
        q2 = np.median(vals)
        high = vals > q2
        low = vals <= q2
    elif mode == "pos_vs_zero":
        high = vals > 0
        low = vals == 0
    elif mode == "top20_vs_zero":
        q80 = np.quantile(vals, 0.80)
        high = vals >= q80
        low = vals == 0
        high &= ~low
    elif mode == "top10_vs_bot50":
        q50, q90 = np.quantile(vals, [0.50, 0.90])
        high = vals >= q90
        low = vals <= q50
        high &= ~low
    else:
        raise KeyError(mode)
    return high, low


def main():
    os.makedirs(OUT, exist_ok=True)
    print("loading", flush=True)
    obs, loge = load()
    sample = obs["sample"]
    cell_type = obs["cell_type"]
    malignant = np.array([cell_type[i] == TUMOR[sample[i]] for i in range(len(sample))])
    immune = np.isin(cell_type, list(IMMUNE))
    cyto = np.isin(cell_type, list(CYTO))

    # Per malignant cell, per radius, a dict of score arrays stored as columns
    # Build a long cell table after all radii: one row per malignant cell.
    pieces = []
    for s in SAMPLES:
        local = np.flatnonzero(sample == s)
        mal_local = malignant[local]
        xy = obs["xy"][local]
        imm = immune[local].astype(np.float32)
        cyt = cyto[local].astype(np.float32)
        cldn4 = loge["CLDN4"][local]
        fov = obs["fov"][local]
        print(f"{s} n={len(local)} malignant={int(mal_local.sum())}", flush=True)
        base = {
            "sample": np.full(int(mal_local.sum()), s, dtype=object),
            "patient": np.full(int(mal_local.sum()), PATIENT[s], dtype=object),
            "fov": fov[mal_local],
            "cldn4": cldn4[mal_local],
        }
        # gene matrices on the local index
        chemo = np.column_stack([loge[g][local] for g in CHEMO])
        isg = np.column_stack([loge[g][local] for g in ISG])
        broad = np.column_stack([loge[g][local] for g in BROAD])
        singles = {g: loge[g][local] for g in SINGLES}
        chemo_sum = chemo.sum(axis=1)
        isg_sum = isg.sum(axis=1)
        broad_sum = broad.sum(axis=1)
        chemo_pos = (chemo > 0).any(axis=1).astype(np.float32)
        for radius in RADII:
            print(f"  radius {radius}", flush=True)
            A = adjacency(xy, float(radius))
            deg = np.asarray(A.sum(axis=1)).ravel()
            n_imm = np.asarray(A @ imm).ravel()
            n_cyt = np.asarray(A @ cyt).ravel()
            sum_chemo = np.asarray(A @ chemo_sum).ravel()
            sum_isg = np.asarray(A @ isg_sum).ravel()
            sum_broad = np.asarray(A @ broad_sum).ravel()
            n_chemo_pos = np.asarray(A @ chemo_pos).ravel()
            # immune-only chemokine sum: zero non-immune columns
            A_imm = A @ sparse.diags(imm)
            sum_chemo_imm = np.asarray(A_imm @ chemo_sum).ravel()
            n_chemo_pos_imm = np.asarray(A_imm @ chemo_pos).ravel()
            scores = {
                "n_immune": n_imm,
                "n_cyto": n_cyt,
                "frac_immune": np.divide(n_imm, deg, out=np.zeros_like(n_imm), where=deg > 0),
                "sum_chemokine": sum_chemo,
                "mean_chemokine": np.divide(sum_chemo, deg, out=np.zeros_like(sum_chemo), where=deg > 0),
                "n_chemokine_pos": n_chemo_pos,
                "sum_chemokine_immune": sum_chemo_imm,
                "n_chemokine_pos_immune": n_chemo_pos_imm,
                "sum_isg": sum_isg,
                "mean_isg": np.divide(sum_isg, deg, out=np.zeros_like(sum_isg), where=deg > 0),
                "sum_broad": sum_broad,
                "mean_broad": np.divide(sum_broad, deg, out=np.zeros_like(sum_broad), where=deg > 0),
            }
            for g, vec in singles.items():
                sm = np.asarray(A @ vec).ravel()
                pos = (vec > 0).astype(np.float32)
                npos = np.asarray(A @ pos).ravel()
                scores[f"sum_{g}"] = sm
                scores[f"mean_{g}"] = np.divide(sm, deg, out=np.zeros_like(sm), where=deg > 0)
                scores[f"npos_{g}"] = npos
            block = {k: v[mal_local] for k, v in scores.items()}
            block = {f"r{radius}_{k}": v for k, v in block.items()}
            base.update(block)
            del A, A_imm
        pieces.append(pd.DataFrame(base))
    cells = pd.concat(pieces, ignore_index=True)
    print("cells", len(cells), "cols", cells.shape[1], flush=True)
    cells.to_pickle("/tmp/ifn_cold_cell_scores.pkl")

    score_names = [c[len("r20_"):] for c in cells.columns if c.startswith("r20_")]
    score_cols = [c for c in cells.columns if c.startswith("r")]
    cuts = ("q4_q1", "median", "pos_vs_zero", "top20_vs_zero", "top10_vs_bot50", "fov_q4_q1")
    rows = []
    # Precompute arm labels per sample (section cuts) and per fov
    for mode in cuts:
        arm = np.full(len(cells), "mid", dtype=object)
        if mode == "fov_q4_q1":
            for (s, fov), idx in cells.groupby(["sample", "fov"]).groups.items():
                idx = np.asarray(list(idx))
                hi, lo = arms_from_values(cells.loc[idx, "cldn4"].to_numpy(), "q4_q1")
                # arms_from_values requires len>=20; groupby index may be smaller
                if len(idx) < 20:
                    continue
                arm[idx[hi]] = "high"
                arm[idx[lo]] = "low"
        else:
            for s, idx in cells.groupby("sample").groups.items():
                idx = np.asarray(list(idx))
                hi, lo = arms_from_values(cells.loc[idx, "cldn4"].to_numpy(), mode)
                arm[idx[hi]] = "high"
                arm[idx[lo]] = "low"
        cells["_arm"] = arm
        use = cells[cells["_arm"].isin(["high", "low"])]
        sec_n = use.groupby(["sample", "_arm"]).size()
        sec_mean = use.groupby(["sample", "_arm"])[score_cols].mean()
        fov_n = use.groupby(["sample", "fov", "_arm"]).size()
        fov_mean = use.groupby(["sample", "fov", "_arm"])[score_cols].mean()
        for radius in RADII:
            for score in score_names:
                col = f"r{radius}_{score}"
                sec_hi, sec_lo, sec_names = [], [], []
                for s in SAMPLES:
                    try:
                        nh = int(sec_n.loc[(s, "high")])
                        nl = int(sec_n.loc[(s, "low")])
                        hv = float(sec_mean.loc[(s, "high"), col])
                        lv = float(sec_mean.loc[(s, "low"), col])
                    except KeyError:
                        continue
                    if nh < 30 or nl < 30:
                        continue
                    sec_hi.append(hv)
                    sec_lo.append(lv)
                    sec_names.append(s)
                st = wilcoxon(sec_hi, sec_lo)
                pat_hi, pat_lo = [], []
                if sec_names:
                    tmp = pd.DataFrame({"sample": sec_names, "hi": sec_hi, "lo": sec_lo})
                    tmp["patient"] = tmp["sample"].map(PATIENT)
                    for _, g in tmp.groupby("patient"):
                        pat_hi.append(float(g["hi"].mean()))
                        pat_lo.append(float(g["lo"].mean()))
                pt = wilcoxon(pat_hi, pat_lo)
                fov_hi, fov_lo = [], []
                if len(fov_n):
                    wide_n = fov_n.unstack("_arm")
                    if "high" in wide_n.columns and "low" in wide_n.columns:
                        ok = (wide_n["high"] >= 8) & (wide_n["low"] >= 8)
                        ok_idx = wide_n.index[ok.fillna(False)]
                        hi_m = fov_mean[col].xs("high", level="_arm")
                        lo_m = fov_mean[col].xs("low", level="_arm")
                        both = hi_m.index.intersection(lo_m.index).intersection(ok_idx)
                        fov_hi = hi_m.loc[both].to_numpy(dtype=float)
                        fov_lo = lo_m.loc[both].to_numpy(dtype=float)
                ft = wilcoxon(fov_hi, fov_lo)
                rows.append(
                    {
                        "score": score,
                        "radius_um": radius,
                        "cutoff": mode,
                        "n_sections": st["n"],
                        "sections_high_lt_low": f"{st['n_neg']}/{st['n']}" if st["n"] else "NA",
                        "section_median_delta": st["median_delta"],
                        "section_p": st["p"],
                        "patients_high_lt_low": f"{pt['n_neg']}/{pt['n']}" if pt["n"] else "NA",
                        "patient_median_delta": pt["median_delta"],
                        "patient_p": pt["p"],
                        "n_fovs": ft["n"],
                        "fovs_high_lt_low": f"{ft['n_neg']}/{ft['n']}" if ft["n"] else "NA",
                        "fov_median_delta": ft["median_delta"],
                        "fov_p": ft["p"],
                        "section_mean_high": float(np.mean(sec_hi)) if sec_hi else np.nan,
                        "section_mean_low": float(np.mean(sec_lo)) if sec_lo else np.nan,
                    }
                )
        print(f"cutoff {mode} done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "ifn_cold_search.csv"), index=False)
    # rank exclusion-direction specs
    neg = res[res["section_median_delta"] < 0].copy()
    neg["section_p_rank"] = neg["section_p"].fillna(1)
    neg["fov_p_rank"] = neg["fov_p"].fillna(1)
    neg = neg.sort_values(["section_p_rank", "fov_p_rank", "section_median_delta"])
    neg.head(40).to_csv(os.path.join(OUT, "ifn_cold_search_top_negative.csv"), index=False)
    print("--- top negative by section p then fov p ---", flush=True)
    cols = ["score", "radius_um", "cutoff", "sections_high_lt_low", "section_median_delta", "section_p", "patients_high_lt_low", "patient_p", "n_fovs", "fovs_high_lt_low", "fov_median_delta", "fov_p"]
    print(neg[cols].head(25).to_string(index=False))
    print("--- most negative section median delta among 8/8 ---", flush=True)
    full = neg[neg["sections_high_lt_low"] == "8/8"].sort_values("section_median_delta")
    print(full[cols].head(15).to_string(index=False))
    print("n specs", len(res), "n negative section delta", len(neg), flush=True)


if __name__ == "__main__":
    main()
