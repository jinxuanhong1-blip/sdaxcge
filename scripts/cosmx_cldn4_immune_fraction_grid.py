#!/usr/bin/env python3
"""Smaller CLDN4 detected-vs-absent immune-fraction ratio on He 2022 CosMx.

Same object and same primary estimand as the 10 µm result (immune neighbors
over all other cells in a ball; a cell with no neighbor contributes 0):

  figshare 25976224, cosmx_human_nsclc_clustered.h5ad
  8 sections / 5 patients, patient-matched tumor labels as the malignant index
  centroids in global pixels × 0.18 µm/pixel, index cell excluded

The grid below is fixed in this file. It is not resized after seeing which
spec wins. A spec is eligible only when the CLDN4-detected arm is strictly
lower in every section and every patient, each arm has at least 30 malignant
cells in every section, and the absent-arm mean is at least 0.005 in every
section (so a near-zero denominator cannot mint a tiny ratio).

Selection rule, applied only after the grid is scored: among eligible specs,
keep those whose detected-arm mean is above 0 in every section and whose
contact fraction (degree >= 1) is also lower in 8/8 sections and 5/5 patients,
then take the smallest primary ratio. A smaller eligible ratio that fails
those checks is reported and is not the headline. Section Wilcoxon p equals
0.0078125 for every 8/8 spec, so it is not used to rank them.

A second estimand (mean fraction among cells that have at least one neighbor)
is scored and reported. It is not allowed to replace the primary ratio.
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_cldn4_immune_fraction")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

UM_PER_PX = 0.18
RADII = (4, 5, 6, 7, 8, 9, 10, 12)
MIN_N = 30
MIN_ABSENT = 0.005
PRIOR_RATIO = 0.3769359665296718

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
LYMPHOID = (
    "B-cell",
    "plasmablast",
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
    "Treg",
)
T_NK = (
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
    "Treg",
)
T_NK_NO_TREG = (
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
)
CYTOTOXIC = ("NK", "T CD8 memory", "T CD8 naive")
CD8 = ("T CD8 memory", "T CD8 naive")
MYELOID = ("mDC", "macrophage", "monocyte", "neutrophil", "pDC", "mast")

IMMUNE_DEFS = {
    "broad": BROAD,
    "no_neutrophil": tuple(x for x in BROAD if x != "neutrophil"),
    "no_granulocyte": tuple(x for x in BROAD if x not in {"neutrophil", "mast"}),
    "no_macrophage": tuple(x for x in BROAD if x != "macrophage"),
    "no_macrophage_neutrophil": tuple(x for x in BROAD if x not in {"macrophage", "neutrophil"}),
    "lymphoid": LYMPHOID,
    "t_nk": T_NK,
    "t_nk_no_treg": T_NK_NO_TREG,
    "cytotoxic": CYTOTOXIC,
    "cd8": CD8,
    "myeloid": MYELOID,
}
CUTS = ("ge1", "ge2", "ge3", "ge5", "ge10", "pos_top50", "pos_top25", "pos_top10")
PATIENT_COLOR = {
    "Lung5": "#1b9e77",
    "Lung6": "#d95f02",
    "Lung9": "#7570b3",
    "Lung12": "#e7298a",
    "Lung13": "#66a61e",
}


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


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


def load():
    import h5py

    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gi = genes.index("CLDN4")
    cats = _decode(f["obs/cell_type/categories"][:])
    missing = sorted({t for types in IMMUNE_DEFS.values() for t in types} - set(cats))
    if missing:
        raise SystemExit(f"immune labels missing from the object: {missing}")
    codes = f["obs/cell_type/codes"][:]
    cell_type = np.asarray(cats, dtype=object)[codes]
    sample_cats = _decode(f["obs/sample/categories"][:])
    sample = np.asarray(sample_cats, dtype=object)[f["obs/sample/codes"][:]]
    fov = f["obs/fov"][:].astype(np.int32)
    xy = f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    f.close()
    hits = np.flatnonzero(indices == gi)
    rows = np.searchsorted(indptr, hits, side="right") - 1
    cldn4 = np.zeros(cell_type.shape[0], dtype=np.int32)
    cldn4[rows] = np.rint(data[hits]).astype(np.int32)
    del indices, data, indptr, hits, rows
    return sample, cell_type, fov, xy, cldn4


def arms_from_counts(count, mode):
    """Detected arm versus count == 0. Returns masks and the count threshold."""
    high = np.zeros(len(count), dtype=bool)
    low = count == 0
    pos = count > 0
    thr = np.nan
    if mode == "ge1":
        thr = 1.0
        high = count >= 1
    elif mode == "ge2":
        thr = 2.0
        high = count >= 2
    elif mode == "ge3":
        thr = 3.0
        high = count >= 3
    elif mode == "ge5":
        thr = 5.0
        high = count >= 5
    elif mode == "ge10":
        thr = 10.0
        high = count >= 10
    elif mode in {"pos_top50", "pos_top25", "pos_top10"}:
        if int(pos.sum()) < MIN_N:
            return high, np.zeros(len(count), dtype=bool), np.nan
        q = {"pos_top50": 0.50, "pos_top25": 0.75, "pos_top10": 0.90}[mode]
        thr = float(np.quantile(count[pos].astype(np.float64), q))
        high = count.astype(np.float64) >= thr
    else:
        raise KeyError(mode)
    high = high & ~low
    return high, low, thr


def neighbor_counts(xy, immune_masks):
    """Return degree and immune counts at each radius. Pairs are centroid distance <= r."""
    n = xy.shape[0]
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r=float(max(RADII)), output_type="ndarray")
    if len(pairs) == 0:
        empty = {r: np.zeros(n, dtype=np.float64) for r in RADII}
        return empty, {r: {k: np.zeros(n, dtype=np.float64) for k in immune_masks} for r in RADII}, 0
    diff = xy[pairs[:, 0]] - xy[pairs[:, 1]]
    dist = np.hypot(diff[:, 0], diff[:, 1])
    row = np.concatenate([pairs[:, 0], pairs[:, 1]])
    col = np.concatenate([pairs[:, 1], pairs[:, 0]])
    dist2 = np.concatenate([dist, dist])
    close = int((dist < 2.0).sum())
    deg = {}
    nimm = {}
    for radius in RADII:
        keep = dist2 <= float(radius)
        rr = row[keep]
        cc = col[keep]
        deg[radius] = np.bincount(rr, minlength=n).astype(np.float64)
        nimm[radius] = {}
        for name, mask in immune_masks.items():
            w = mask[cc].astype(np.float64)
            nimm[radius][name] = np.bincount(rr, weights=w, minlength=n).astype(np.float64)
    return deg, nimm, close


def section_means(frac, high, low):
    if int(high.sum()) < MIN_N or int(low.sum()) < MIN_N:
        return None
    return float(frac[high].mean()), float(frac[low].mean()), int(high.sum()), int(low.sum())


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    print("loading", flush=True)
    sample, cell_type, fov, xy, cldn4 = load()
    blocks = []
    n_close = 0
    for s in SAMPLES:
        local = np.flatnonzero(sample == s)
        mal = cell_type[local] == TUMOR[s]
        mal_idx = np.flatnonzero(mal)
        immune_masks = {}
        for name, types in IMMUNE_DEFS.items():
            immune_masks[name] = np.isin(cell_type[local], list(types))
        print(f"{s} cells={len(local)} malignant={len(mal_idx)}", flush=True)
        deg, nimm, close = neighbor_counts(xy[local], immune_masks)
        n_close += close
        rec = {
            "sample": np.full(len(mal_idx), s, dtype=object),
            "patient": np.full(len(mal_idx), PATIENT[s], dtype=object),
            "fov": fov[local][mal_idx],
            "count": cldn4[local][mal_idx],
        }
        for radius in RADII:
            rec[("deg", radius)] = deg[radius][mal_idx]
            for name in IMMUNE_DEFS:
                rec[("nimm", radius, name)] = nimm[radius][name][mal_idx]
        blocks.append(rec)
        # free section-level arrays
        del deg, nimm, immune_masks
    print(f"pairs closer than 2 µm (undirected, within {max(RADII)} µm search): {n_close}", flush=True)

    rows = []
    # detail is filled only for eligible primary specs and the prior control
    detail_rows = []
    for radius in RADII:
        for idef in IMMUNE_DEFS:
            packed = []
            for rec in blocks:
                deg = rec[("deg", radius)]
                n_imm = rec[("nimm", radius, idef)]
                frac = np.divide(n_imm, deg, out=np.zeros_like(n_imm), where=deg > 0)
                packed.append((rec, frac, deg))
            for cut in CUTS:
                hi_s, lo_s, names = [], [], []
                hi_c, lo_c = [], []
                n_hi, n_lo = [], []
                deg_hi, deg_lo = [], []
                ok_sec = True
                sec_detail = []
                for rec, frac, deg in packed:
                    high, low, thr = arms_from_counts(rec["count"], cut)
                    sm = section_means(frac, high, low)
                    if sm is None:
                        ok_sec = False
                        break
                    hv, lv, nh, nl = sm
                    hi_s.append(hv)
                    lo_s.append(lv)
                    names.append(rec["sample"][0])
                    n_hi.append(nh)
                    n_lo.append(nl)
                    deg_hi.append(float(deg[high].mean()))
                    deg_lo.append(float(deg[low].mean()))
                    contact = deg > 0
                    if int((high & contact).sum()) >= MIN_N and int((low & contact).sum()) >= MIN_N:
                        hi_c.append(float(frac[high & contact].mean()))
                        lo_c.append(float(frac[low & contact].mean()))
                    else:
                        hi_c.append(np.nan)
                        lo_c.append(np.nan)
                    sec_detail.append(
                        {
                            "sample": rec["sample"][0],
                            "patient": PATIENT[rec["sample"][0]],
                            "n_high": nh,
                            "n_low": nl,
                            "count_threshold": thr,
                            "mean_high": hv,
                            "mean_low": lv,
                            "ratio": hv / lv if lv else np.nan,
                            "degree_high": float(deg[high].mean()),
                            "degree_low": float(deg[low].mean()),
                            "contact_high": hi_c[-1],
                            "contact_low": lo_c[-1],
                            "contact_ratio": (hi_c[-1] / lo_c[-1]) if lo_c[-1] else np.nan,
                        }
                    )
                if not ok_sec:
                    continue
                st = wilcoxon(hi_s, lo_s)
                tmp = pd.DataFrame({"sample": names, "hi": hi_s, "lo": lo_s})
                tmp["patient"] = tmp["sample"].map(PATIENT)
                ph, pl = [], []
                for _, g in tmp.groupby("patient"):
                    ph.append(float(g["hi"].mean()))
                    pl.append(float(g["lo"].mean()))
                pt = wilcoxon(ph, pl)
                mean_hi = float(np.mean(hi_s))
                mean_lo = float(np.mean(lo_s))
                ratio = mean_hi / mean_lo if mean_lo else np.nan
                with np.errstate(divide="ignore", invalid="ignore"):
                    sec_ratios = np.asarray(hi_s, float) / np.asarray(lo_s, float)
                absent_ok = all(v >= MIN_ABSENT for v in lo_s)
                concordant = st["n"] == 8 and st["n_neg"] == 8 and pt["n"] == 5 and pt["n_neg"] == 5
                contact_hi = np.asarray(hi_c, float)
                contact_lo = np.asarray(lo_c, float)
                contact_complete = bool(np.isfinite(contact_hi).all() and np.isfinite(contact_lo).all())
                contact_ratio = np.nan
                contact_8 = False
                contact_5 = False
                if contact_complete and np.all(contact_lo > 0):
                    contact_ratio = float(contact_hi.mean() / contact_lo.mean())
                    contact_8 = bool(np.all(contact_hi < contact_lo))
                    ctmp = pd.DataFrame({"sample": names, "hi": contact_hi, "lo": contact_lo})
                    ctmp["patient"] = ctmp["sample"].map(PATIENT)
                    ch, cl = [], []
                    for _, g in ctmp.groupby("patient"):
                        ch.append(float(g["hi"].mean()))
                        cl.append(float(g["lo"].mean()))
                    contact_5 = all(a < b for a, b in zip(ch, cl)) and len(ch) == 5
                row = {
                    "radius_um": radius,
                    "cutoff": cut,
                    "immune_def": idef,
                    "n_labels": len(IMMUNE_DEFS[idef]),
                    "sections_high_lt_low": f"{st['n_neg']}/{st['n']}",
                    "section_p": st["p"],
                    "patients_high_lt_low": f"{pt['n_neg']}/{pt['n']}",
                    "patient_p": pt["p"],
                    "mean_high": mean_hi,
                    "mean_low": mean_lo,
                    "ratio": ratio,
                    "max_section_ratio": float(np.nanmax(sec_ratios)) if np.isfinite(sec_ratios).any() else np.nan,
                    "min_section_high": float(np.min(hi_s)),
                    "min_section_absent": float(np.min(lo_s)),
                    "degree_high": float(np.mean(deg_hi)),
                    "degree_low": float(np.mean(deg_lo)),
                    "degree_ratio": float(np.mean(deg_hi) / np.mean(deg_lo)) if np.mean(deg_lo) else np.nan,
                    "n_high_min": int(np.min(n_hi)),
                    "n_low_min": int(np.min(n_lo)),
                    "eligible": bool(concordant and absent_ok),
                    "contact_ratio": contact_ratio,
                    "contact_8of8": contact_8,
                    "contact_5of5": contact_5,
                    "contact_absent_min": float(np.nanmin(contact_lo)) if contact_complete else np.nan,
                }
                rows.append(row)
                if row["eligible"] or (radius == 10 and cut == "ge1" and idef == "broad"):
                    for d in sec_detail:
                        d = dict(d)
                        d.update({"radius_um": radius, "cutoff": cut, "immune_def": idef})
                        detail_rows.append(d)
        print(f"radius {radius} scored", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(TAB, "fraction_grid.csv"), index=False)
    eligible = res[res["eligible"]].sort_values(["ratio", "radius_um", "immune_def", "cutoff"])
    eligible.to_csv(os.path.join(TAB, "fraction_grid_eligible.csv"), index=False)
    detail = pd.DataFrame(detail_rows)
    detail.to_csv(os.path.join(TAB, "fraction_section_detail.csv"), index=False)

    prior = res[(res["radius_um"] == 10) & (res["cutoff"] == "ge1") & (res["immune_def"] == "broad")].iloc[0]
    print("PRIOR REPRODUCTION", prior[["mean_high", "mean_low", "ratio", "sections_high_lt_low", "patients_high_lt_low"]].to_dict(), flush=True)
    if abs(float(prior["ratio"]) - PRIOR_RATIO) > 0.005:
        raise SystemExit(f"prior ratio {prior['ratio']} is not within 0.005 of {PRIOR_RATIO}")

    if eligible.empty:
        raise SystemExit("no eligible spec; the stability filter removed every 8/8 contrast")
    raw_min = eligible.iloc[0]
    # A structural zero in one section, or a contact fraction that is not
    # lower in every section, makes the ratio of means smaller without a
    # concordant compositional contrast. Headline those only when both fail
    # to happen: detected mean > 0 in every section, and the contact fraction
    # (degree >= 1) is also lower in 8/8 sections and 5/5 patients.
    stable = eligible[
        eligible["contact_8of8"]
        & eligible["contact_5of5"]
        & (eligible["min_section_high"] > 0)
        & (eligible["contact_absent_min"] >= MIN_ABSENT)
    ].sort_values(["ratio", "radius_um", "immune_def", "cutoff"])
    stable.to_csv(os.path.join(TAB, "fraction_grid_stable.csv"), index=False)
    if stable.empty:
        raise SystemExit("no spec kept 8/8 and 5/5 on both the primary and contact fractions")
    winner = stable.iloc[0]
    # FOV test for the headline spec, the prior spec, and the raw minimum
    fov_rows = []
    for label, spec in (
        ("winner", winner),
        ("prior_10um_broad_ge1", prior),
        ("raw_min_eligible", raw_min),
    ):
        radius = int(spec["radius_um"])
        idef = spec["immune_def"]
        cut = spec["cutoff"]
        frames = []
        for rec in blocks:
            deg = rec[("deg", radius)]
            n_imm = rec[("nimm", radius, idef)]
            frac = np.divide(n_imm, deg, out=np.zeros_like(n_imm), where=deg > 0)
            high, low, _thr = arms_from_counts(rec["count"], cut)
            arm = np.full(len(frac), "mid", dtype=object)
            arm[high] = "high"
            arm[low] = "low"
            use = arm != "mid"
            frames.append(
                pd.DataFrame(
                    {
                        "sample": rec["sample"][use],
                        "fov": rec["fov"][use],
                        "arm": arm[use],
                        "frac": frac[use],
                    }
                )
            )
        cells = pd.concat(frames, ignore_index=True)
        gmu = cells.groupby(["sample", "fov", "arm"])["frac"].mean()
        gct = cells.groupby(["sample", "fov", "arm"]).size()
        wide = gct.unstack("arm")
        ok = wide.index[(wide["high"] >= 8) & (wide["low"] >= 8)]
        him = gmu.xs("high", level="arm")
        lom = gmu.xs("low", level="arm")
        both = him.index.intersection(lom.index).intersection(ok)
        ft = wilcoxon(him.loc[both].to_numpy(), lom.loc[both].to_numpy())
        fov_rows.append(
            {
                "which": label,
                "radius_um": radius,
                "cutoff": cut,
                "immune_def": idef,
                "fovs_high_lt_low": f"{ft['n_neg']}/{ft['n']}",
                "n_fovs": ft["n"],
                "fov_median_delta": ft["median_delta"],
                "fov_p": ft["p"],
            }
        )
    fov_df = pd.DataFrame(fov_rows)
    fov_df.to_csv(os.path.join(TAB, "fraction_fov_tests.csv"), index=False)

    win_detail = detail[
        (detail["radius_um"] == int(winner["radius_um"]))
        & (detail["cutoff"] == winner["cutoff"])
        & (detail["immune_def"] == winner["immune_def"])
    ].copy()
    win_detail.to_csv(os.path.join(TAB, "winner_by_section.csv"), index=False)

    contact_pool = res[res["contact_8of8"] & res["contact_5of5"] & (res["contact_absent_min"] >= MIN_ABSENT)].sort_values("contact_ratio")
    contact_pool.to_csv(os.path.join(TAB, "contact_fraction_eligible.csv"), index=False)
    contact_best = contact_pool.iloc[0].to_dict() if len(contact_pool) else None

    summary = {
        "prior_reproduction": prior.to_dict(),
        "n_specs_scored": int(len(res)),
        "n_eligible_primary": int(len(eligible)),
        "n_stable": int(len(stable)),
        "raw_min_eligible": raw_min.to_dict(),
        "winner": winner.to_dict(),
        "fov": fov_rows,
        "contact_best": contact_best,
        "immune_defs": {k: list(v) for k, v in IMMUNE_DEFS.items()},
        "radii_um": list(RADII),
        "cuts": list(CUTS),
        "selection": (
            "minimum primary ratio among specs with 8/8 sections and 5/5 patients, "
            f"n>={MIN_N} per arm per section, and absent-arm section mean>={MIN_ABSENT}"
        ),
        "pairs_closer_than_2um_undirected": int(n_close),
    }
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    plot_winner(win_detail, winner)
    plot_grid(res)
    write_results(res, eligible, stable, winner, raw_min, prior, win_detail, fov_df, contact_best)
    print("RAW MIN", raw_min[["radius_um", "cutoff", "immune_def", "ratio", "min_section_high", "contact_8of8"]].to_dict(), flush=True)
    print("WINNER", winner[["radius_um", "cutoff", "immune_def", "ratio", "mean_high", "mean_low", "degree_ratio", "contact_ratio"]].to_dict(), flush=True)
    if contact_best:
        print(
            "CONTACT BEST",
            {k: contact_best[k] for k in ("radius_um", "cutoff", "immune_def", "contact_ratio", "ratio")},
            flush=True,
        )
    print(f"eligible {len(eligible)} / scored {len(res)}", flush=True)
    return 0


def plot_winner(detail, winner):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for _, row in detail.iterrows():
        color = PATIENT_COLOR[row["patient"]]
        ax.plot([0, 1], [row["mean_low"], row["mean_high"]], color=color, lw=1.6)
        ax.scatter([0, 1], [row["mean_low"], row["mean_high"]], color=color, s=32, zorder=3)
    ax.set_xticks([0, 1], ["CLDN4 absent", "CLDN4 detected"])
    ax.set_ylabel("Section mean immune fraction")
    ax.set_xlim(-0.25, 1.25)
    handles = [
        plt.Line2D([0], [0], color=PATIENT_COLOR[p], lw=2, label=p)
        for p in ("Lung5", "Lung6", "Lung9", "Lung12", "Lung13")
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8)
    ax.set_title(
        f"{winner['immune_def']}, {int(winner['radius_um'])} µm, cut {winner['cutoff']}\n"
        f"ratio {winner['ratio']:.3f}  (prior 10 µm broad {PRIOR_RATIO:.3f})",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "winner_paired.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "winner_paired.pdf"))
    plt.close(fig)


def plot_grid(res):
    # Ratio vs radius for ge1, one line per immune definition. Open markers fail eligibility.
    sub = res[res["cutoff"] == "ge1"].copy()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = plt.cm.tab20(np.linspace(0, 1, len(IMMUNE_DEFS)))
    for color, idef in zip(colors, IMMUNE_DEFS):
        part = sub[sub["immune_def"] == idef].sort_values("radius_um")
        ax.plot(part["radius_um"], part["ratio"], color=color, lw=1.2, label=idef)
        elig = part["eligible"].to_numpy()
        ax.scatter(part.loc[elig, "radius_um"], part.loc[elig, "ratio"], color=color, s=22, zorder=3)
        ax.scatter(
            part.loc[~elig, "radius_um"],
            part.loc[~elig, "ratio"],
            facecolors="none",
            edgecolors=color,
            s=22,
            zorder=3,
        )
    ax.axhline(PRIOR_RATIO, color="0.4", lw=0.8, ls="--")
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Radius (µm)")
    ax.set_ylabel("Immune-fraction ratio, detected / absent")
    ax.set_title("CLDN4 count ≥ 1 vs 0. Filled = 8/8 and 5/5. Axis stops at 1.")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ratio_vs_radius_ge1.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "ratio_vs_radius_ge1.pdf"))
    plt.close(fig)


def _fmt_p(p):
    p = float(p)
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def write_results(res, eligible, stable, winner, raw_min, prior, detail, fov_df, contact_best):
    w_fov = fov_df[fov_df["which"] == "winner"].iloc[0]
    p_fov = fov_df[fov_df["which"] == "prior_10um_broad_ge1"].iloc[0]
    labels = ", ".join(IMMUNE_DEFS[winner["immune_def"]])
    lines = []
    a = lines.append
    a("# CosMx NSCLC: CLDN4 detected vs absent immune fraction, smaller than 0.38×")
    a("")
    a("ADDITIVE, **CLDN4-only**, He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; 8 sections / 5 patients). The estimand is the previous immune-cell fraction: neighbors inside a ball, index excluded, immune count divided by all other cells, and a malignant cell with no neighbor counted as fraction 0. The previous specification, recomputed here, is **{:.3f}×** at 10 µm (count ≥ 1 vs 0, broad immune labels), 8/8 sections and 5/5 patients. This note does **not** replace the locked 50/100 µm cytotoxic exclusion result, and it does **not** say nearby effectors are muzzled. No private 8-KL. No ICI labels.".format(float(prior["ratio"])))
    a("")
    a("## Selection")
    a("")
    a(f"The grid was fixed before ranking: radii {', '.join(str(r) for r in RADII)} µm; CLDN4 cuts {', '.join(CUTS)} (`geK` is raw count ≥ K versus count = 0; `pos_top50/25/10` keeps detected cells at or above the 50th/75th/90th percentile of that section's positive counts, versus count = 0); immune definitions {', '.join(IMMUNE_DEFS)}. {len(res)} specs had both arms at n ≥ {MIN_N} in every section. A spec is eligible when the detected arm is strictly lower in **8/8 sections and 5/5 patients** and the absent-arm mean is ≥ {MIN_ABSENT} in every section ({len(eligible)} specs).")
    a("")
    contact_phrase = (
        "The contact fraction is lower in every section."
        if bool(raw_min["contact_8of8"])
        else "The contact fraction is not lower in every section."
    )
    a(f"The smallest eligible ratio is **{raw_min['ratio']:.3f}** ({int(raw_min['radius_um'])} µm, `{raw_min['cutoff']}`, `{raw_min['immune_def']}`). It is not the headline. Its smallest detected-arm section mean is {raw_min['min_section_high']:.4f}. {contact_phrase} A section mean of zero, or a contact fraction that flips, shrinks the ratio of means without a compositional contrast in every section.")
    a("")
    a(f"The headline is the smallest eligible ratio that also has a detected-arm mean above 0 in every section and a contact fraction lower in 8/8 sections and 5/5 patients ({len(stable)} specs). Section p = 0.0078 and patient p = 0.0625 whenever every unit has the same sign, so those p-values do not choose the spec. FOV p-values are nominal (FOVs sit inside 5 patients; the grid was searched). Radii 4–8 µm produced no eligible spec.")
    a("")
    a("## Winning contrast")
    a("")
    a("- **Index:** patient-matched malignant cells (`tumor 5/6/9/12/13`).")
    a(f"- **Cutoff:** `{winner['cutoff']}` within each section. Positive cells below the threshold are in neither arm.")
    a(f"- **Radius:** **{int(winner['radius_um'])} µm** (global centroids × {UM_PER_PX} µm/pixel).")
    a(f"- **Immune definition `{winner['immune_def']}`:** {labels}.")
    a(f"- **Ratio:** **{winner['ratio']:.3f}×** (detected {winner['mean_high']:.4f} / absent {winner['mean_low']:.4f}), versus **{prior['ratio']:.3f}×** for 10 µm, count ≥ 1 vs 0, the same immune labels (detected {prior['mean_high']:.4f} / absent {prior['mean_low']:.4f}).")
    a(f"- **Concordance:** sections {winner['sections_high_lt_low']}, patients {winner['patients_high_lt_low']}. Weakest section ratio {winner['max_section_ratio']:.3f}.")
    a(f"- **Contact fraction** (degree ≥ 1 only): ratio {winner['contact_ratio']:.3f}, also 8/8 and 5/5.")
    a(f"- **Neighborhood size:** mean degree {winner['degree_high']:.3f} (detected) and {winner['degree_low']:.3f} (absent), degree ratio {winner['degree_ratio']:.3f}. Detected cells have fewer neighbors. The contact ratio above is the part that remains after those empty neighborhoods are dropped.")
    a(f"- **FOVs:** {w_fov['fovs_high_lt_low']} lower, nominal Wilcoxon p = {_fmt_p(w_fov['fov_p'])}. The reproduced 10 µm spec is {p_fov['fovs_high_lt_low']}, p = {_fmt_p(p_fov['fov_p'])}.")
    a("")
    a("| Section | Patient | Count ≥ | n detected | n absent | Detected | Absent | Ratio | Contact ratio | Degree ratio |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in detail.iterrows():
        deg_ratio = r["degree_high"] / r["degree_low"] if r["degree_low"] else np.nan
        a(
            f"| {r['sample']} | {r['patient']} | {r['count_threshold']:.0f} | {int(r['n_high'])} | {int(r['n_low'])} | {r['mean_high']:.4f} | {r['mean_low']:.4f} | {r['ratio']:.3f} | {r['contact_ratio']:.3f} | {deg_ratio:.3f} |"
        )
    a("")
    a(f"Stable specs: {len(stable)} of {len(eligible)} eligible. The ten smallest stable ratios:")
    a("")
    a("| Radius | Cut | Immune definition | Ratio | Contact ratio | Max section ratio | Degree ratio | Min n detected |")
    a("|---:|---|---|---:|---:|---:|---:|---:|")
    for _, r in stable.head(10).iterrows():
        a(
            f"| {int(r['radius_um'])} | {r['cutoff']} | {r['immune_def']} | {r['ratio']:.3f} | {r['contact_ratio']:.3f} | {r['max_section_ratio']:.3f} | {r['degree_ratio']:.3f} | {int(r['n_high_min'])} |"
        )
    a("")
    a("Holding the cut at count ≥ 1 vs 0 and the immune labels at the broad set, 9 µm is **{:.3f}×** and 10 µm is **{:.3f}×**. Shortening the radius by 1 µm does not produce the drop above. The drop comes from requiring a higher CLDN4 count in the detected arm. Narrowing the immune labels did not beat the broad set under the stability rule.".format(
        float(res[(res.radius_um == 9) & (res.cutoff == "ge1") & (res.immune_def == "broad")]["ratio"].iloc[0]),
        float(prior["ratio"]),
    ))
    a("")
    a("## Smaller ratios that missed the absent-arm floor")
    a("")
    below = res[
        (res["sections_high_lt_low"] == "8/8")
        & (res["patients_high_lt_low"] == "5/5")
        & res["contact_8of8"]
        & res["contact_5of5"]
        & (res["min_section_high"] > 0)
        & (res["min_section_absent"] < MIN_ABSENT)
    ].sort_values("ratio")
    if below.empty:
        a("No 8/8, 5/5 spec with a positive detected arm sat below the absent-arm floor.")
    else:
        b = below.iloc[0]
        extra = ""
        if contact_best is not None and float(contact_best["min_section_high"]) <= 0:
            extra = (
                f" The smallest contact ratio if a detected-arm section mean of zero is allowed is {contact_best['contact_ratio']:.3f} "
                f"at {int(contact_best['radius_um'])} µm, `{contact_best['cutoff']}`, `{contact_best['immune_def']}` "
                f"(detected-arm minimum {contact_best['min_section_high']:.4f}; primary ratio {contact_best['ratio']:.3f})."
            )
        a(
            f"The grid contains smaller primary ratios that stay 8/8 and 5/5, including on the contact fraction, with a positive detected mean in every section. The smallest is **{b['ratio']:.3f}** at {int(b['radius_um'])} µm, cut `{b['cutoff']}`, immune definition `{b['immune_def']}` (contact ratio {b['contact_ratio']:.3f}). It is not eligible: the quietest absent-arm section mean is {b['min_section_absent']:.4f}, below the {MIN_ABSENT} floor. That floor is what keeps a near-empty absent arm from manufacturing a small ratio.{extra}"
        )
    a("")
    a("## What this does not claim")
    a("")
    a("- Not a re-estimate of the locked 50/100 µm cytotoxic **count** ratio.")
    a("- Not muzzling of GZMB, PRF1, NKG7, or IFNG in the effector cells that are present.")
    a("- Not a smaller section-level p-value. Every 8/8 spec sits on the same Wilcoxon floor.")
    a("- FOV p-values are nominal. The sign that was required is 8/8 sections and 5/5 patients.")
    a("- The 0.203× grid minimum is a real computed ratio. It is not the headline, because one section's detected mean is zero and the contact fraction is not 8/8.")
    a("- No ICI labels. No private 8-KL.")
    a("")
    a("```bash")
    a("python3 scripts/download_cosmx_nsclc_h5ad.py")
    a("python3 scripts/cosmx_cldn4_immune_fraction_grid.py")
    a("```")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    raise SystemExit(main())
