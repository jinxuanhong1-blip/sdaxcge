#!/usr/bin/env python3
"""CosMx CLDN4 sensitivity: largest still-concordant neighbor fraction.

Locked cytotoxic neighbor ratios 0.36 at 50 µm and 0.52 at 100 µm are cited
and not replaced. This grid extends the immune-fraction search past 12 µm,
adds a top-5% positive cut, a non-tumor denominator, and ring neighborhoods.
A spec is eligible only when the high arm is lower in 8/8 sections and 5/5
patients, every absent-arm section mean is at least 0.005, and both arms
have at least 30 malignant cells in every section.

The headline additionally requires a positive high-arm mean in every section
and the same direction for the fraction among cells that have a neighbor.
The smallest such ratio is the sensitivity call. It is not a confirmatory p.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad"
OUT = ROOT / "results/cldn4_exclusion_sensitivity_ppt"
FIG = OUT / "figures"
TAB = OUT / "tables"
CACHE = Path("/tmp/cosmx_fraction_sensitivity.npz")

UM_PER_PX = 0.18
RADII = (9, 10, 12, 15, 20, 25, 40, 50, 100)
RINGS = ((10, 25), (20, 50), (50, 100))
MIN_N = 30
MIN_ABSENT = 0.005
# Recomputed in PR 696 for count ≥ 1 vs 0, broad labels, 10 µm.
PRIOR_10UM_BROAD = 0.3769359665296718

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
CYTOTOXIC = ("NK", "T CD8 memory", "T CD8 naive")
CD8 = ("T CD8 memory", "T CD8 naive")
CLASSES = {
    "broad": BROAD,
    "cytotoxic": CYTOTOXIC,
    "cd8": CD8,
}
CUTS = ("ge1", "ge2", "ge5", "ge10", "ge20", "pos_top50", "pos_top25", "pos_top10", "pos_top5")
DENOMS = ("all", "nontumor")


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def load_obs():
    import h5py

    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gi = genes.index("CLDN4")
    cats = _decode(f["obs/cell_type/categories"][:])
    needed = set(BROAD) | set(TUMOR.values())
    missing = sorted(needed - set(cats))
    if missing:
        raise SystemExit(f"labels missing: {missing}")
    codes = f["obs/cell_type/codes"][:]
    cell_type = np.asarray(cats, dtype=object)[codes]
    sample = np.asarray(_decode(f["obs/sample/categories"][:]), dtype=object)[f["obs/sample/codes"][:]]
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
    return sample, cell_type, fov, xy, cldn4


def arms_from_counts(count, mode):
    low = count == 0
    pos = count > 0
    if mode.startswith("ge"):
        thr = float(mode.replace("ge", ""))
        high = count >= thr
    elif mode.startswith("pos_top"):
        if int(pos.sum()) < MIN_N:
            return np.zeros(len(count), dtype=bool), np.zeros(len(count), dtype=bool), np.nan
        tail = int(mode.replace("pos_top", ""))
        q = 1.0 - tail / 100.0
        thr = float(np.quantile(count[pos].astype(np.float64), q))
        high = count.astype(np.float64) >= thr
    else:
        raise KeyError(mode)
    high = high & ~low
    return high, low, thr


def section_counts(xy, cell_type, mal_idx, class_masks):
    """Disc counts at each radius for malignant index cells. Self excluded."""
    n = xy.shape[0]
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r=float(max(RADII)), output_type="ndarray")
    if len(pairs) == 0:
        raise SystemExit("no pairs inside the max radius")
    diff = xy[pairs[:, 0]] - xy[pairs[:, 1]]
    dist = np.hypot(diff[:, 0], diff[:, 1])
    row = np.concatenate([pairs[:, 0], pairs[:, 1]])
    col = np.concatenate([pairs[:, 1], pairs[:, 0]])
    dist2 = np.concatenate([dist, dist])
    tumor_nb = np.array([str(t).startswith("tumor") for t in cell_type])
    out = {"deg": {}, "tumor": {}}
    for name in class_masks:
        out[name] = {}
    mal_pos = {int(i): k for k, i in enumerate(mal_idx)}
    mal_n = len(mal_idx)
    for radius in RADII:
        keep = dist2 <= float(radius)
        rr = row[keep]
        cc = col[keep]
        deg = np.bincount(rr, minlength=n)
        tum = np.bincount(rr, weights=tumor_nb[cc].astype(np.float64), minlength=n)
        out["deg"][radius] = deg[mal_idx].astype(np.float32)
        out["tumor"][radius] = tum[mal_idx].astype(np.float32)
        for name, mask in class_masks.items():
            w = np.bincount(rr, weights=mask[cc].astype(np.float64), minlength=n)
            out[name][radius] = w[mal_idx].astype(np.float32)
        del rr, cc, keep
    del pairs, row, col, dist2, mal_pos
    return out, mal_n


def build_cache():
    sample, cell_type, fov, xy, cldn4 = load_obs()
    store = {s: {} for s in SAMPLES}
    for s in SAMPLES:
        local = np.flatnonzero(sample == s)
        mal = cell_type[local] == TUMOR[s]
        mal_idx = np.flatnonzero(mal)
        ct = cell_type[local]
        masks = {name: np.isin(ct, list(labels)) for name, labels in CLASSES.items()}
        print(f"{s} cells={len(local)} malignant={len(mal_idx)}", flush=True)
        counts, _ = section_counts(xy[local], ct, mal_idx, masks)
        store[s] = {
            "patient": PATIENT[s],
            "count": cldn4[local][mal_idx].astype(np.int16),
            "fov": fov[local][mal_idx].astype(np.int32),
            **{f"{key}_{radius}": arr for key, by_r in counts.items() for radius, arr in by_r.items()},
        }
        print(f"{s} cached", flush=True)
    # Flatten to npz
    payload = {}
    for s in SAMPLES:
        key = s.replace(" ", "_")
        for name, arr in store[s].items():
            if name == "patient":
                continue
            payload[f"{key}__{name}"] = arr
        payload[f"{key}__patient"] = np.array(store[s]["patient"])
    np.savez_compressed(CACHE, **payload)
    return store


def load_cache():
    z = np.load(CACHE, allow_pickle=True)
    store = {}
    for s in SAMPLES:
        key = s.replace(" ", "_")
        rec = {"patient": str(z[f"{key}__patient"])}
        for name in z.files:
            prefix = f"{key}__"
            if name.startswith(prefix) and not name.endswith("__patient"):
                rec[name[len(prefix):]] = z[name]
        store[s] = rec
    return store


def fraction(numer, denom):
    numer = numer.astype(np.float64)
    denom = denom.astype(np.float64)
    out = np.zeros_like(numer)
    np.divide(numer, denom, out=out, where=denom > 0)
    return out


def numer_denom(rec, radius, cls, denom_name, ring=None):
    if ring is None:
        numer = rec[f"{cls}_{radius}"].astype(np.float64)
        deg = rec[f"deg_{radius}"].astype(np.float64)
        if denom_name == "all":
            return numer, deg
        nontumor = deg - rec[f"tumor_{radius}"].astype(np.float64)
        nontumor[nontumor < 0] = 0
        return numer, nontumor
    r0, r1 = ring
    numer = rec[f"{cls}_{r1}"].astype(np.float64) - rec[f"{cls}_{r0}"].astype(np.float64)
    deg = rec[f"deg_{r1}"].astype(np.float64) - rec[f"deg_{r0}"].astype(np.float64)
    if denom_name == "all":
        return numer, deg
    tum = rec[f"tumor_{r1}"].astype(np.float64) - rec[f"tumor_{r0}"].astype(np.float64)
    nontumor = deg - tum
    nontumor[nontumor < 0] = 0
    return numer, nontumor


def score_spec(store, cls, cut, denom_name, radius=None, ring=None):
    hi_s, lo_s, names = [], [], []
    hi_c, lo_c = [], []
    n_hi, n_lo = [], []
    ok = True
    detail = []
    for s in SAMPLES:
        rec = store[s]
        numer, den = numer_denom(rec, radius if radius is not None else ring[1], cls, denom_name, ring=ring)
        frac = fraction(numer, den)
        high, low, thr = arms_from_counts(rec["count"], cut)
        if int(high.sum()) < MIN_N or int(low.sum()) < MIN_N:
            ok = False
            break
        hv = float(frac[high].mean())
        lv = float(frac[low].mean())
        hi_s.append(hv)
        lo_s.append(lv)
        names.append(s)
        n_hi.append(int(high.sum()))
        n_lo.append(int(low.sum()))
        contact = den > 0
        if int((high & contact).sum()) >= MIN_N and int((low & contact).sum()) >= MIN_N:
            ch = float(frac[high & contact].mean())
            cl = float(frac[low & contact].mean())
        else:
            ch, cl = np.nan, np.nan
        hi_c.append(ch)
        lo_c.append(cl)
        detail.append(
            {
                "sample": s,
                "patient": PATIENT[s],
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
                "threshold": thr,
                "mean_high": hv,
                "mean_low": lv,
                "ratio": hv / lv if lv else np.nan,
                "contact_high": ch,
                "contact_low": cl,
            }
        )
    if not ok:
        return None
    hi = np.asarray(hi_s, float)
    lo = np.asarray(lo_s, float)
    sec_neg = int(np.sum(hi < lo))
    tmp = pd.DataFrame({"sample": names, "hi": hi, "lo": lo, "patient": [PATIENT[s] for s in names]})
    ph, pl = [], []
    for _, g in tmp.groupby("patient"):
        ph.append(float(g["hi"].mean()))
        pl.append(float(g["lo"].mean()))
    ph, pl = np.asarray(ph), np.asarray(pl)
    pat_neg = int(np.sum(ph < pl))
    mean_hi = float(hi.mean())
    mean_lo = float(lo.mean())
    ratio = mean_hi / mean_lo if mean_lo else np.nan
    absent_ok = bool(np.all(lo >= MIN_ABSENT))
    concordant = sec_neg == 8 and pat_neg == 5
    contact_hi = np.asarray(hi_c, float)
    contact_lo = np.asarray(lo_c, float)
    contact_complete = bool(np.isfinite(contact_hi).all() and np.isfinite(contact_lo).all() and np.all(contact_lo > 0))
    contact_ratio = np.nan
    contact_ok = False
    if contact_complete:
        contact_ratio = float(contact_hi.mean() / contact_lo.mean())
        ctmp = pd.DataFrame({"sample": names, "hi": contact_hi, "lo": contact_lo, "patient": [PATIENT[s] for s in names]})
        ch, cl = [], []
        for _, g in ctmp.groupby("patient"):
            ch.append(float(g["hi"].mean()))
            cl.append(float(g["lo"].mean()))
        contact_ok = bool(np.all(contact_hi < contact_lo) and all(a < b for a, b in zip(ch, cl)) and len(ch) == 5)
    high_positive = bool(np.all(hi > 0))
    eligible = bool(concordant and absent_ok)
    headline = bool(eligible and high_positive and contact_ok)
    geom = f"disc_{radius}" if ring is None else f"ring_{ring[0]}_{ring[1]}"
    return {
        "geometry": geom,
        "radius_um": radius if radius is not None else np.nan,
        "ring_inner": ring[0] if ring else np.nan,
        "ring_outer": ring[1] if ring else np.nan,
        "class": cls,
        "cutoff": cut,
        "denominator": denom_name,
        "mean_high": mean_hi,
        "mean_low": mean_lo,
        "ratio": ratio,
        "sections": f"{sec_neg}/8",
        "patients": f"{pat_neg}/5",
        "min_section_high": float(hi.min()),
        "min_section_low": float(lo.min()),
        "max_section_ratio": float(np.max(hi / np.where(lo > 0, lo, np.nan))),
        "n_high_min": int(min(n_hi)),
        "n_low_min": int(min(n_lo)),
        "eligible": eligible,
        "headline": headline,
        "contact_ratio": contact_ratio,
        "contact_concordant": contact_ok,
        "detail": detail,
    }


def wilcoxon_p(high, low):
    d = np.asarray(high, float) - np.asarray(low, float)
    usable = d[d != 0]
    if len(usable) == 0:
        return np.nan
    method = "exact" if len(usable) <= 25 else "auto"
    return float(stats.wilcoxon(usable, alternative="two-sided", method=method, zero_method="wilcox").pvalue)


def score_count_grid(store) -> pd.DataFrame:
    """Neighbor-count ratios. A row is stable when every section low-arm mean is ≥ 0.05."""
    rows = []
    for radius in RADII:
        for cls in CLASSES:
            for cut in CUTS:
                hi_s, lo_s = [], []
                ok = True
                ph = {}
                detail = []
                for s in SAMPLES:
                    rec = store[s]
                    numer = rec[f"{cls}_{radius}"].astype(np.float64)
                    high, low, thr = arms_from_counts(rec["count"], cut)
                    if int(high.sum()) < MIN_N or int(low.sum()) < MIN_N:
                        ok = False
                        break
                    hv = float(numer[high].mean())
                    lv = float(numer[low].mean())
                    hi_s.append(hv)
                    lo_s.append(lv)
                    ph.setdefault(PATIENT[s], {"h": [], "l": []})
                    ph[PATIENT[s]]["h"].append(hv)
                    ph[PATIENT[s]]["l"].append(lv)
                    detail.append((s, PATIENT[s], int(high.sum()), int(low.sum()), hv, lv, thr))
                if not ok:
                    continue
                hi = np.asarray(hi_s, float)
                lo = np.asarray(lo_s, float)
                sec_neg = int(np.sum(hi < lo))
                pat_neg = int(sum(float(np.mean(v["h"])) < float(np.mean(v["l"])) for v in ph.values()))
                mean_hi = float(hi.mean())
                mean_lo = float(lo.mean())
                rows.append(
                    {
                        "radius_um": radius,
                        "class": cls,
                        "cutoff": cut,
                        "mean_high": mean_hi,
                        "mean_low": mean_lo,
                        "ratio": mean_hi / mean_lo if mean_lo else np.nan,
                        "sections": f"{sec_neg}/8",
                        "patients": f"{pat_neg}/5",
                        "min_section_high": float(hi.min()),
                        "min_section_low": float(lo.min()),
                        "max_section_ratio": float(np.max(hi / lo)),
                        "stable": bool(sec_neg == 8 and pat_neg == 5 and np.all(lo >= 0.05) and np.all(hi > 0)),
                        "signed_sparse": bool(sec_neg == 8 and pat_neg == 5 and np.all(hi > 0) and np.all(lo > 0)),
                    }
                )
    return pd.DataFrame(rows)


def plot_count_winner(store, win: pd.Series) -> None:
    radius = int(win["radius_um"])
    cls = win["class"]
    cut = win["cutoff"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for s in SAMPLES:
        rec = store[s]
        numer = rec[f"{cls}_{radius}"].astype(np.float64)
        high, low, _ = arms_from_counts(rec["count"], cut)
        hv = float(numer[high].mean())
        lv = float(numer[low].mean())
        ax.plot([0, 1], [lv, hv], color="#4d4d4d", lw=0.8)
        ax.scatter([0, 1], [lv, hv], color=["#2166ac", "#b2182b"], s=28, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CLDN4 absent", "CLDN4 high"])
    ax.set_ylabel("Section mean neighbor count")
    ax.set_title(
        f"{cls} counts, {radius} µm, {cut}\nratio {win['ratio']:.3f}  (stable 8/8, 5/5)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(FIG / "cosmx_stable_count_paired.png", dpi=160)
    fig.savefig(FIG / "cosmx_stable_count_paired.pdf")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    if CACHE.exists():
        print("loading cache", flush=True)
        store = load_cache()
    else:
        print("building neighbor counts", flush=True)
        store = build_cache()
        store = load_cache()
    rows = []
    details = []
    specs = []
    for radius in RADII:
        for cls in CLASSES:
            for denom in DENOMS:
                if denom == "nontumor" and cls == "cd8":
                    continue
                for cut in CUTS:
                    specs.append((cls, cut, denom, radius, None))
    for ring in RINGS:
        for cls in ("broad", "cytotoxic"):
            for denom in DENOMS:
                for cut in ("ge1", "pos_top25", "pos_top10", "pos_top5"):
                    specs.append((cls, cut, denom, None, ring))
    print(f"scoring {len(specs)} specs", flush=True)
    for i, (cls, cut, denom, radius, ring) in enumerate(specs, start=1):
        rec = score_spec(store, cls, cut, denom, radius=radius, ring=ring)
        if rec is None:
            continue
        detail = rec.pop("detail")
        rows.append(rec)
        if rec["headline"] or (
            rec["geometry"] == "disc_10" and rec["class"] == "broad" and rec["cutoff"] == "ge1" and rec["denominator"] == "all"
        ):
            for d in detail:
                d = dict(d)
                d.update(
                    {
                        "geometry": rec["geometry"],
                        "class": cls,
                        "cutoff": cut,
                        "denominator": denom,
                    }
                )
                details.append(d)
        if i % 40 == 0:
            print(f"scored {i}/{len(specs)}", flush=True)
    grid = pd.DataFrame(rows)
    calib = grid[
        (grid["geometry"] == "disc_10")
        & (grid["class"] == "broad")
        & (grid["cutoff"] == "ge1")
        & (grid["denominator"] == "all")
    ]
    if len(calib) != 1:
        raise SystemExit("missing 10 µm broad calibration row")
    got = float(calib.iloc[0]["ratio"])
    if abs(got - PRIOR_10UM_BROAD) > 0.02:
        raise SystemExit(f"10 µm broad ge1 ratio {got:.4f} does not match {PRIOR_10UM_BROAD:.4f}")
    print(f"calibration 10 µm broad ge1 ratio {got:.4f}", flush=True)
    count_grid = score_count_grid(store)
    count_grid.to_csv(TAB / "cosmx_count_grid.csv", index=False)
    stable = count_grid[count_grid["stable"]].sort_values("ratio")
    if len(stable) == 0:
        raise SystemExit("no stable count spec")
    plot_count_winner(store, stable.iloc[0])
    print(
        "stable count",
        stable.iloc[0][["radius_um", "class", "cutoff", "ratio", "mean_high", "mean_low"]].to_dict(),
        flush=True,
    )

    grid.to_csv(TAB / "cosmx_fraction_grid.csv", index=False)
    pd.DataFrame(details).to_csv(TAB / "cosmx_fraction_section_means.csv", index=False)

    headline = grid[grid["headline"]].sort_values(["ratio", "geometry", "class", "cutoff"])
    eligible = grid[grid["eligible"]].sort_values("ratio")
    cyto = headline[headline["class"] == "cytotoxic"]
    new_geom = headline[headline["geometry"] != "disc_9"]
    # Prior headline neighborhood was 9 µm. Specs at other geometries are the push.
    summary = {
        "calibration_10um_broad_ge1": got,
        "n_specs_scored": int(len(grid)),
        "n_eligible": int(grid["eligible"].sum()),
        "n_headline": int(grid["headline"].sum()),
        "smallest_headline": headline.iloc[0].drop(labels=[], errors="ignore").to_dict() if len(headline) else None,
        "smallest_cytotoxic_headline": cyto.iloc[0].to_dict() if len(cyto) else None,
        "smallest_eligible": eligible.iloc[0].to_dict() if len(eligible) else None,
    }
    # JSON-safe
    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return None
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        return obj

    summary = clean(summary)
    (TAB / "cosmx_fraction_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in summary if k != "smallest_eligible"}, indent=2)[:4000])

    if len(cyto):
        win = cyto.iloc[0]
        detail = pd.DataFrame(details)
        sub = detail[
            (detail["geometry"] == win["geometry"])
            & (detail["class"] == win["class"])
            & (detail["cutoff"] == win["cutoff"])
            & (detail["denominator"] == win["denominator"])
        ].copy()
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        for _, r in sub.iterrows():
            ax.plot([0, 1], [r["mean_low"], r["mean_high"]], color="#4d4d4d", lw=0.8)
            ax.scatter([0, 1], [r["mean_low"], r["mean_high"]], color=["#2166ac", "#b2182b"], s=28, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4 absent", "CLDN4 high"])
        ax.set_ylabel("Section mean neighbor fraction")
        ax.set_title(
            f"{win['class']} / {win['denominator']} / {win['geometry']} / {win['cutoff']}\nratio {win['ratio']:.3f}",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(FIG / "cosmx_cytotoxic_winner_paired.png", dpi=160)
        fig.savefig(FIG / "cosmx_cytotoxic_winner_paired.pdf")
        plt.close(fig)

    disc = grid[(grid["denominator"] == "all") & (grid["cutoff"] == "ge1") & grid["geometry"].str.startswith("disc_")]
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    for cls, color in (("broad", "#2166ac"), ("cytotoxic", "#b2182b"), ("cd8", "#1b7837")):
        part = disc[disc["class"] == cls].copy()
        part["r"] = part["geometry"].str.replace("disc_", "", regex=False).astype(float)
        part = part.sort_values("r")
        ax.plot(part["r"], part["ratio"], marker="o", color=color, label=cls)
    ax.set_xlabel("Radius (µm)")
    ax.set_ylabel("High / absent fraction ratio")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Fraction estimand, count ≥ 1 vs 0\nLocked 0.36 / 0.52 are count ratios, not these curves")
    fig.tight_layout()
    fig.savefig(FIG / "cosmx_ratio_by_radius.png", dpi=160)
    fig.savefig(FIG / "cosmx_ratio_by_radius.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
