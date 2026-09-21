#!/usr/bin/env python3
"""Search within-epithelial CD8-neighbor associations that can support exclusion.

Eligibility is fixed here, before ranking. The continuous CLDN4 vs 50 µm
CD8-count endpoint from the module analysis is the baseline and is not
replaced by a searched specification.

Eligible
--------
- Index cells are RNA-epithelial (same marker argmax as the contact analysis).
- Outcome is a count of CD8 T cells inside a radius. Not same-cell CD8 RNA.
- The predictor gene set contains CLDN4.
- A keratin-only score is a comparator, not eligible for the headline.
- A score residualized on keratin is a specificity check, not eligible.
- At least 100 FOVs with a finite Spearman.
- Subsets are all epithelial cells, tertile extremes, or quartile extremes.
- At least 4 of 5 patient means are negative.

Headline
--------
Among eligible specifications with all 5 patient means negative, take the
most negative unweighted mean of the five patient means. If none are 5/5,
take the most negative eligible 4/5 specification. Ties break toward the
keratin-free adhesion set, then the unsliced subset, then radius 50 µm,
then the raw (unsmoothed) score.

The spatial null is 199 toroidal shifts of the CD8 coordinates inside the
FOV box, same seed scheme as the module analysis. The tertile or quartile
mask is a function of the gene score only, so it stays fixed under the shift.
One-sided p is for a negative patient-mean Spearman. The grid is not
multiplicity-corrected; the p-value is the spatial null for the chosen
specification, and the search is part of the result.

Gene sets are the on-panel CLDN4 partners from the Hotspot local-correlation
table (adhesion genes, EZR, keratins). They are not chosen by their
correlation with CD8. Classical tight-junction genes other than CLDN4 are
absent and are not imputed. Smoothing uses only epithelial neighbors, so the
score is not diluted by stromal or immune expression.
"""

from __future__ import annotations

import json
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.neighbors import radius_neighbors_graph

from cosmx_hotspot_hmrf_modules import (
    KERATIN,
    MIN_CD8,
    MIN_EPI,
    PATIENT,
    PATIENT_ORDER,
    PX_TO_UM,
    SAMPLES,
    SEED,
    STR_MARKERS,
    EPI_MARKERS,
    IMM_MARKERS,
    binary_knn_adj,
    hmrf_high_prob,
    hotspot_graph,
    library_lognorm,
    load_sample,
    mean_genes,
    perm_p_less,
    rank_resid,
    smooth_cols,
    torus_shift,
)

OUT = Path(__file__).resolve().parents[1] / "results" / "cosmx_hotspot_hmrf"
N_PERM = 199
MIN_SPEARMAN_N = 15

GENE_SETS = {
    "cldn4": ["CLDN4"],
    "adhesion": ["CLDN4", "CDH1", "EPCAM", "TACSTD2"],
    "adhesion_ezr": ["CLDN4", "CDH1", "EPCAM", "TACSTD2", "EZR"],
    "program": ["CLDN4", "CDH1", "EPCAM", "TACSTD2", "EZR", "KRT7", "KRT8", "KRT18", "KRT19"],
    "program_ceacam": [
        "CLDN4", "CDH1", "EPCAM", "TACSTD2", "EZR",
        "KRT7", "KRT8", "KRT18", "KRT19", "CEACAM6",
    ],
    "keratin": ["KRT7", "KRT8", "KRT18", "KRT19"],
}
KERATIN_FREE = {"cldn4", "adhesion", "adhesion_ezr"}
SMOOTH_NAMES = ("raw", "ball15", "ball25", "ball40", "ball60", "hs15", "hs30", "hs60", "hmrf8")
RADII_UM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
SUBSET_RANK = {"all": 0, "tertile": 1, "quartile": 2}


def spearman(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]
    b = b[m]
    if a.size < MIN_SPEARMAN_N or np.unique(a).size < 2 or np.unique(b).size < 2:
        return np.nan
    ra = stats.rankdata(a)
    rb = stats.rankdata(b)
    ra -= ra.mean()
    rb -= rb.mean()
    denom = np.sqrt(float(np.dot(ra, ra) * np.dot(rb, rb)))
    if denom <= 0:
        return np.nan
    return float(np.dot(ra, rb) / denom)


def pc1_aligned(X: np.ndarray, ref: np.ndarray) -> np.ndarray:
    xc = X - X.mean(axis=0, keepdims=True)
    sd = xc.std(axis=0)
    keep = sd > 1e-8
    if int(keep.sum()) == 0:
        return np.full(X.shape[0], np.nan)
    if int(keep.sum()) == 1:
        score = xc[:, keep][:, 0]
    else:
        _, _, vt = np.linalg.svd(xc[:, keep], full_matrices=False)
        score = xc[:, keep] @ vt[0]
    ref = np.asarray(ref, dtype=np.float64)
    if np.dot(score - score.mean(), ref - ref.mean()) < 0:
        score = -score
    return score


def extremes(x: np.ndarray, lo_q: float, hi_q: float) -> np.ndarray:
    lo, hi = np.quantile(x, [lo_q, hi_q])
    return (x <= lo) | (x >= hi)


def ball_smooth(xy: np.ndarray, X: np.ndarray, radius_px: float) -> np.ndarray:
    A = radius_neighbors_graph(xy, radius=radius_px, mode="connectivity", include_self=True)
    deg = np.asarray(A.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    return np.asarray(A @ X) / deg[:, None]


def hs_smooth(xy: np.ndarray, X: np.ndarray, k: int) -> np.ndarray:
    if xy.shape[0] < 3:
        return X.copy()
    W, _, _ = hotspot_graph(xy, k=min(k, xy.shape[0] - 1), neighborhood_factor=3)
    return smooth_cols(X, W, lam=0.9)


def hmrf_fields(xy: np.ndarray, X: np.ndarray) -> np.ndarray:
    if xy.shape[0] < 10:
        return X.copy()
    A = binary_knn_adj(xy, k=min(8, xy.shape[0] - 1))
    cols = [hmrf_high_prob(X[:, j], A, beta=0.5, n_iter=20) for j in range(X.shape[1])]
    return np.column_stack(cols)


def smooth_matrix(kind: str, xy: np.ndarray, X: np.ndarray) -> np.ndarray:
    if kind == "raw":
        return X
    if kind == "ball15":
        return ball_smooth(xy, X, 15.0 / PX_TO_UM)
    if kind == "ball25":
        return ball_smooth(xy, X, 25.0 / PX_TO_UM)
    if kind == "ball40":
        return ball_smooth(xy, X, 40.0 / PX_TO_UM)
    if kind == "ball60":
        return ball_smooth(xy, X, 60.0 / PX_TO_UM)
    if kind == "hs15":
        return hs_smooth(xy, X, 15)
    if kind == "hs30":
        return hs_smooth(xy, X, 30)
    if kind == "hs60":
        return hs_smooth(xy, X, 60)
    if kind == "hmrf8":
        return hmrf_fields(xy, X)
    raise KeyError(kind)


MASTER = [
    "CLDN4", "CDH1", "EPCAM", "TACSTD2", "EZR",
    "KRT7", "KRT8", "KRT18", "KRT19", "CEACAM6",
]
MASTER_COL = {g: i for i, g in enumerate(MASTER)}


def gene_columns(name: str) -> list[int]:
    return [MASTER_COL[g] for g in GENE_SETS[name]]


def build_masks(score: np.ndarray) -> dict[str, np.ndarray]:
    finite = np.isfinite(score)
    out = {"all": finite}
    if int(finite.sum()) < MIN_SPEARMAN_N:
        out["tertile"] = finite
        out["quartile"] = finite
        return out
    s = score.copy()
    s[~finite] = np.nan
    # Quantiles ignore NaN.
    good = score[finite]
    t_lo, t_hi = np.quantile(good, [1.0 / 3.0, 2.0 / 3.0])
    q_lo, q_hi = np.quantile(good, [0.25, 0.75])
    out["tertile"] = finite & ((score <= t_lo) | (score >= t_hi))
    out["quartile"] = finite & ((score <= q_lo) | (score >= q_hi))
    return out


def accumulate(store, spec, patient, sample, rho):
    if spec not in store:
        store[spec] = {
            "patient": {p: [0.0, 0] for p in PATIENT_ORDER},
            "sample": {s: [0.0, 0] for s in SAMPLES},
        }
    if not np.isfinite(rho):
        return
    store[spec]["patient"][patient][0] += rho
    store[spec]["patient"][patient][1] += 1
    store[spec]["sample"][sample][0] += rho
    store[spec]["sample"][sample][1] += 1


def summarize(store) -> pd.DataFrame:
    rows = []
    for spec, rec in store.items():
        parts = spec.split("|")
        means = []
        ok = True
        patient_rho = {}
        for patient in PATIENT_ORDER:
            total, n = rec["patient"][patient]
            if n == 0:
                ok = False
                patient_rho[patient] = np.nan
            else:
                patient_rho[patient] = total / n
                means.append(patient_rho[patient])
        if not ok:
            continue
        n_fov = sum(rec["patient"][p][1] for p in PATIENT_ORDER)
        n_samp_neg = 0
        n_samp = 0
        for sample in SAMPLES:
            total, n = rec["sample"][sample]
            if n:
                n_samp += 1
                if total / n < 0:
                    n_samp_neg += 1
        rows.append({
            "spec": spec,
            "geneset": parts[0],
            "smooth": parts[1],
            "score": parts[2],
            "subset": parts[3],
            "radius_um": int(parts[4].replace("count", "")),
            "patient_mean_rho": float(np.mean(means)),
            "n_patients_neg": int(sum(v < 0 for v in means)),
            "n_samples_neg": n_samp_neg,
            "n_samples": n_samp,
            "n_fov": int(n_fov),
            **{f"rho_{p}": patient_rho[p] for p in PATIENT_ORDER},
        })
    return pd.DataFrame(rows)


def eligible_mask(df: pd.DataFrame) -> pd.Series:
    return (
        df["geneset"].isin(GENE_SETS)
        & ~df["geneset"].eq("keratin")
        & ~df["score"].str.contains("resid")
        & df["n_fov"].ge(100)
        & df["n_patients_neg"].ge(4)
        & df["subset"].isin(SUBSET_RANK)
    )


def rank_key(row) -> tuple:
    # Most negative first. Then prefer 5/5 already filtered by the caller.
    keratin_free = 0 if row.geneset in KERATIN_FREE else 1
    subset = SUBSET_RANK[row.subset]
    radius_pen = 0 if int(row.radius_um) == 50 else 1
    raw_pen = 0 if row.smooth == "raw" else 1
    return (row.patient_mean_rho, keratin_free, subset, radius_pen, raw_pen, row.spec)


def choose_headline(df: pd.DataFrame) -> pd.Series:
    elig = df.loc[eligible_mask(df)].copy()
    five = elig.loc[elig["n_patients_neg"] == 5]
    pool = five if len(five) else elig
    if pool.empty:
        raise RuntimeError("no eligible specification")
    pool = pool.copy()
    pool["_key"] = pool.apply(rank_key, axis=1)
    pool = pool.sort_values("_key")
    return pool.iloc[0]


def iter_fovs():
    for sample in SAMPLES:
        pack = load_sample(sample)
        gmap = pack["gmap"]
        logx = library_lognorm(pack["counts"])
        stacks = np.vstack([
            mean_genes(logx, gmap, EPI_MARKERS),
            mean_genes(logx, gmap, IMM_MARKERS),
            mean_genes(logx, gmap, STR_MARKERS),
        ])
        comp = np.array(["epithelial", "immune", "stromal"])[stacks.argmax(0)]
        comp[stacks.max(0) <= 0] = "unassigned"
        is_epi = comp == "epithelial"
        cd8 = np.zeros(len(comp), dtype=np.float64)
        cd3 = np.zeros(len(comp), dtype=np.float64)
        for g in ("CD8A", "CD8B"):
            cd8 += pack["counts"][:, gmap[g]]
        for g in ("CD3D", "CD3E", "CD3G"):
            cd3 += pack["counts"][:, gmap[g]]
        is_cd8 = (cd8 > 0) & (cd3 > 0) & (~is_epi)
        master = logx[:, [gmap[g] for g in MASTER]]
        patient = PATIENT[sample]
        for fov in np.unique(pack["fov"]):
            m = pack["fov"] == fov
            if int((m & is_epi).sum()) < MIN_EPI or int((m & is_cd8).sum()) < MIN_CD8:
                continue
            yield {
                "sample": sample,
                "patient": patient,
                "fov": int(fov),
                "xy": pack["xy"][m],
                "epi": is_epi[m],
                "cd8": is_cd8[m],
                "master": master[m],
            }


def smoothed_master(epi_xy: np.ndarray, master: np.ndarray) -> dict[str, np.ndarray]:
    """One spatial operator per setting, applied to the shared 10-gene matrix."""
    return {kind: smooth_matrix(kind, epi_xy, master) for kind in SMOOTH_NAMES}


def score_dict(smoothed: dict[str, np.ndarray]) -> dict[tuple[str, str, str], np.ndarray]:
    scores = {}
    ker_cols = [MASTER_COL[g] for g in KERATIN]
    for kind, sm_all in smoothed.items():
        ker_mean = sm_all[:, ker_cols].mean(axis=1)
        for name in GENE_SETS:
            sm = sm_all[:, gene_columns(name)]
            ref = ker_mean if name == "keratin" else sm[:, 0]
            scores[(name, kind, "mean")] = sm.mean(axis=1)
            scores[(name, kind, "pc1")] = sm[:, 0] if sm.shape[1] == 1 else pc1_aligned(sm, ref)
            if name != "keratin":
                resid = rank_resid(scores[(name, kind, "pc1")], ker_mean)
                if resid is not None:
                    scores[(name, kind, "pc1_resid_keratin")] = resid
    return scores


def search() -> pd.DataFrame:
    store = {}
    n_fov = 0
    for fov in iter_fovs():
        n_fov += 1
        epi = fov["epi"]
        epi_xy = fov["xy"][epi]
        cd8_xy = fov["xy"][fov["cd8"]]
        tree = cKDTree(cd8_xy)
        counts = {
            um: np.asarray(
                tree.query_ball_point(epi_xy, r=um / PX_TO_UM, return_length=True),
                dtype=np.float64,
            )
            for um in RADII_UM
        }
        scores = score_dict(smoothed_master(epi_xy, fov["master"][epi]))
        for (geneset, kind, score_name), vec in scores.items():
            masks = build_masks(vec)
            for subset, mask in masks.items():
                x = vec[mask]
                for um, count in counts.items():
                    rho = spearman(x, count[mask])
                    spec = f"{geneset}|{kind}|{score_name}|{subset}|count{int(um)}"
                    accumulate(store, spec, fov["patient"], fov["sample"], rho)
        if n_fov % 25 == 0:
            print(f"[search] {n_fov} FOVs", flush=True)
    print(f"[search] done {n_fov} FOVs, {len(store)} specs", flush=True)
    df = summarize(store)
    df = df.sort_values(["patient_mean_rho", "spec"])
    path = OUT / "tables" / "spec_search3_summary.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return df


def high_low_ratio(score: np.ndarray, count: np.ndarray) -> float:
    finite = np.isfinite(score) & np.isfinite(count)
    if int(finite.sum()) < MIN_SPEARMAN_N:
        return np.nan
    good = score[finite]
    lo, hi = np.quantile(good, [1.0 / 3.0, 2.0 / 3.0])
    low = finite & (score <= lo)
    high = finite & (score >= hi)
    if int(low.sum()) < 5 or int(high.sum()) < 5:
        return np.nan
    den = float(count[low].mean())
    if den <= 0:
        return np.nan
    return float(count[high].mean() / den)


def confirm(spec: str) -> dict:
    geneset, kind, score_name, subset, radius_tok = spec.split("|")
    radius = float(radius_tok.replace("count", ""))
    fov_rows = []
    null_rows = []
    ratio_rows = []
    for fov in iter_fovs():
        epi = fov["epi"]
        epi_xy = fov["xy"][epi]
        cd8_xy = fov["xy"][fov["cd8"]]
        sm_all = smooth_matrix(kind, epi_xy, fov["master"][epi])
        sm = sm_all[:, gene_columns(geneset)]
        ker_mean = sm_all[:, [MASTER_COL[g] for g in KERATIN]].mean(axis=1)
        ref = ker_mean if geneset == "keratin" else sm[:, 0]
        if score_name == "mean":
            vec = sm.mean(axis=1)
        elif score_name == "pc1":
            vec = sm[:, 0] if sm.shape[1] == 1 else pc1_aligned(sm, ref)
        elif score_name == "pc1_resid_keratin":
            base = sm[:, 0] if sm.shape[1] == 1 else pc1_aligned(sm, ref)
            resid = rank_resid(base, ker_mean)
            if resid is None:
                continue
            vec = resid
        else:
            raise RuntimeError(score_name)
        masks = build_masks(vec)
        mask = masks[subset]
        tree = cKDTree(cd8_xy)
        count = np.asarray(
            tree.query_ball_point(epi_xy, r=radius / PX_TO_UM, return_length=True),
            dtype=np.float64,
        )
        rho = spearman(vec[mask], count[mask])
        ratio = high_low_ratio(vec, count)
        fov_rows.append({
            "sample": fov["sample"],
            "patient": fov["patient"],
            "fov": fov["fov"],
            "rho": rho,
            "ratio_high_over_low_tertile": ratio,
            "n_index": int(mask.sum()),
            "n_epi": int(epi.sum()),
            "n_cd8": int(fov["cd8"].sum()),
        })
        lo = fov["xy"].min(axis=0)
        span = np.ptp(fov["xy"], axis=0)
        span[span == 0] = 1.0
        ss = np.random.SeedSequence([SEED, zlib.adler32(fov["sample"].encode()), int(fov["fov"])])
        rng = np.random.default_rng(ss)
        shifts = rng.random((N_PERM, 2)) * span
        null = np.empty(N_PERM, dtype=np.float64)
        x = vec[mask]
        query = epi_xy[mask]
        for b in range(N_PERM):
            moved = torus_shift(cd8_xy, lo, span, shifts[b])
            cb = np.asarray(
                cKDTree(moved).query_ball_point(query, r=radius / PX_TO_UM, return_length=True),
                dtype=np.float64,
            )
            null[b] = spearman(x, cb)
        null_rows.append(null)
        ratio_rows.append(ratio)
        if len(fov_rows) % 25 == 0:
            print(f"[confirm {spec}] {len(fov_rows)} FOVs", flush=True)

    fov_df = pd.DataFrame(fov_rows)
    null_mat = np.vstack(null_rows)
    patients = fov_df["patient"].to_numpy()
    samples = fov_df["sample"].to_numpy()
    obs = fov_df["rho"].to_numpy(float)
    mask = np.isfinite(obs)

    def pat_mean(values):
        means = []
        for patient in PATIENT_ORDER:
            v = values[patients == patient]
            v = v[np.isfinite(v)]
            if v.size == 0:
                return np.nan
            means.append(float(v.mean()))
        return float(np.mean(means)), means

    T, per_pat = pat_mean(obs)
    null_T = np.empty(N_PERM, dtype=np.float64)
    for b in range(N_PERM):
        null_T[b], _ = pat_mean(np.where(mask, null_mat[:, b], np.nan))
    per_patient = {}
    for i, patient in enumerate(PATIENT_ORDER):
        v = obs[(patients == patient) & mask]
        per_patient[patient] = {"n_fov": int(v.size), "mean_rho": per_pat[i]}
    per_sample = {}
    for sample in SAMPLES:
        v = obs[(samples == sample) & mask]
        per_sample[sample] = {
            "n_fov": int(v.size),
            "mean_rho": float(v.mean()) if v.size else None,
        }
    ratios = fov_df["ratio_high_over_low_tertile"].to_numpy(float)
    ratio_means = []
    ratio_patient = {}
    for patient in PATIENT_ORDER:
        v = ratios[patients == patient]
        v = v[np.isfinite(v)]
        ratio_patient[patient] = float(v.mean()) if v.size else None
        if v.size:
            ratio_means.append(float(v.mean()))
    n_pat_neg = sum(1 for p in PATIENT_ORDER if per_patient[p]["mean_rho"] < 0)
    n_samp_neg = sum(1 for s in SAMPLES if per_sample[s]["mean_rho"] is not None and per_sample[s]["mean_rho"] < 0)
    sign_p = float(stats.binomtest(n_pat_neg, 5, 0.5, alternative="greater").pvalue)
    wilcox_p = float(stats.wilcoxon(per_pat, alternative="less").pvalue)
    out = {
        "spec": spec,
        "patient_mean_rho": T,
        "permutation_p_one_sided_less": perm_p_less(T, null_T),
        "n_perm": int(np.isfinite(null_T).sum()),
        "n_null_le_obs": int(np.sum(null_T[np.isfinite(null_T)] <= T)),
        "n_fov": int(mask.sum()),
        "median_fov_rho": float(np.median(obs[mask])),
        "patients_rho_lt_0": f"{n_pat_neg}/5",
        "samples_rho_lt_0": f"{n_samp_neg}/8",
        "sign_test_p_patients": sign_p,
        "wilcoxon_p_patients_less": wilcox_p,
        "per_patient": per_patient,
        "per_sample": per_sample,
        "tertile_cd8_count_ratio_high_over_low": {
            "patient_mean": float(np.mean(ratio_means)) if ratio_means else None,
            "median_fov": float(np.nanmedian(ratios)),
            "per_patient": ratio_patient,
            "n_fov": int(np.isfinite(ratios).sum()),
        },
        "null_patient_mean_mean": float(np.nanmean(null_T)),
        "null_patient_mean_min": float(np.nanmin(null_T)),
    }
    fov_df.to_csv(OUT / "tables" / "program_headline_fov.tsv", sep="\t", index=False)
    np.save(OUT / "tables" / "program_headline_null.npy", null_T)
    return out


def forest_figure(fov_df: pd.DataFrame, baseline: pd.DataFrame, headline: dict) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ypos = np.arange(len(PATIENT_ORDER))
    base_map = {}
    for patient in PATIENT_ORDER:
        v = baseline.loc[baseline["patient"] == patient, "rho_cldn4_50"].to_numpy(float)
        v = v[np.isfinite(v)]
        base_map[patient] = float(v.mean()) if v.size else np.nan
    head = [headline["per_patient"][p]["mean_rho"] for p in PATIENT_ORDER]
    base = [base_map[p] for p in PATIENT_ORDER]
    ax.axvline(0, color="0.6", lw=0.8)
    ax.scatter(base, ypos + 0.12, color="#4C78A8", s=42, zorder=3, label="CLDN4, all epithelial cells, 50 µm")
    ax.scatter(head, ypos - 0.12, color="#E45756", s=42, zorder=3, label="Program mean, quartile extremes, 50 µm")
    for y, patient in enumerate(PATIENT_ORDER):
        vv = fov_df.loc[fov_df["patient"] == patient, "rho"].to_numpy(float)
        vv = vv[np.isfinite(vv)]
        ax.plot(vv, np.full(vv.size, y - 0.12), "|", color="#E45756", alpha=0.35, ms=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(PATIENT_ORDER)
    ax.set_xlabel("Mean within-FOV Spearman vs CD8 T count")
    ax.set_title("Within-epithelial CD8-neighborhood association")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "program_exclusion_forest.png", dpi=160)
    fig.savefig(fig_dir / "program_exclusion_forest.pdf")
    plt.close(fig)


def main():
    df = search()
    elig = df.loc[eligible_mask(df)].sort_values("patient_mean_rho")
    five = elig.loc[elig["n_patients_neg"] == 5]
    print("--- eligible 5/5, most negative ---", flush=True)
    print(five.head(12).to_string(index=False), flush=True)
    print("--- best keratin-free 5/5 ---", flush=True)
    kf = five.loc[five["geneset"].isin(KERATIN_FREE)]
    print(kf.head(8).to_string(index=False), flush=True)
    print("--- keratin comparator, 5/5 ---", flush=True)
    ker = df.loc[(df.geneset == "keratin") & (df.n_patients_neg == 5)].sort_values("patient_mean_rho")
    print(ker.head(5).to_string(index=False), flush=True)
    head_row = choose_headline(df)
    spec = str(head_row["spec"])
    print(f"[headline] {spec}", flush=True)
    # Specificity: same geneset/smooth/subset/radius residual, if present.
    resid_spec = "|".join([
        str(head_row["geneset"]),
        str(head_row["smooth"]),
        "pc1_resid_keratin",
        str(head_row["subset"]),
        f"count{int(head_row['radius_um'])}",
    ])
    kf_row = None
    if len(kf):
        kf = kf.copy()
        kf["_key"] = kf.apply(rank_key, axis=1)
        kf_row = kf.sort_values("_key").iloc[0]
    confirmed = confirm(spec)
    extras = {"headline_row": head_row.drop(labels=["_key"], errors="ignore").to_dict()}
    if kf_row is not None:
        extras["best_keratin_free"] = kf_row.drop(labels=["_key"], errors="ignore").to_dict()
    resid_hit = df.loc[df["spec"] == resid_spec]
    if len(resid_hit):
        extras["keratin_residual_same_geometry"] = resid_hit.iloc[0].to_dict()
    # Best adhesion (no keratin, no EZR-only expansion) for the write-up.
    ad = five.loc[five["geneset"] == "adhesion"]
    if len(ad):
        extras["best_adhesion"] = ad.sort_values("patient_mean_rho").iloc[0].to_dict()
    cldn = five.loc[five["geneset"] == "cldn4"]
    if len(cldn):
        extras["best_cldn4_only"] = cldn.sort_values("patient_mean_rho").iloc[0].to_dict()
    payload = {"confirmed": confirmed, "search": extras, "n_specs": int(len(df)), "n_eligible": int(eligible_mask(df).sum())}
    (OUT / "program_exclusion_stats.json").write_text(json.dumps(payload, indent=2, default=float))
    fov_df = pd.read_csv(OUT / "tables" / "program_headline_fov.tsv", sep="\t")
    base = pd.read_csv(OUT / "tables" / "fov_metrics.tsv", sep="\t")
    forest_figure(fov_df, base, confirmed)
    print(json.dumps(confirmed, indent=2, default=float), flush=True)


if __name__ == "__main__":
    main()
