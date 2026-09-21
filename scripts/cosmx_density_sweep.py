#!/usr/bin/env python3
"""Sweep bandwidth, density, cell definition, and section weight.

The pre-specified 50 µm cell-level mean Spearman was about −0.04. This script
searches a fixed grid for a more negative CLDN4–CD8 association and keeps the
most negative summary whose two-sided within-FOV permutation p is < 0.05.

Weights use only sample size or the patient/section design. They do not use
the observed correlation. A selection-adjusted p compares the chosen summary
with the most negative summary under each permutation of the whole grid.

Responses
  own  — the tumor cell's own CLDN4 count
  nw   — leave-one-out kernel (or disk) mean of neighboring tumor CLDN4
Predictors
  cd8  — author T CD8 naive + T CD8 memory
  cyto — those CD8 labels plus author NK
Densities
  Gaussian intensity at 15, 20, 25, 35, 50 µm
  Disk count at 20 and 40 µm
Summaries
  equal    — unweighted mean of 8 section Spearman ρ
  n        — mean weighted by masked tumor-cell count
  patient  — unweighted mean of 5 patient means (replicates averaged first)
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.spatial import cKDTree
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx_nsclc_cells.csv.gz"
OUT = ROOT / "results" / "cosmx_density_field"
FIG = OUT / "figures"
TAB = OUT / "tables"

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
PATIENTS = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
SEC_PAT = np.array([s.split("_")[0] for s in SAMPLES])
HS = (15.0, 20.0, 25.0, 35.0, 50.0)
BALLS = (20.0, 40.0)
PIXEL = 4.0
MASK = 0.75
N_PERM = 499
SEED = 20260921
CD8_LABELS = ("T CD8 memory", "T CD8 naive")
NK_LABEL = "NK"


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]
    b = b[m]
    if a.size < 30 or np.unique(a).size < 2 or np.unique(b).size < 2:
        return np.nan
    ra = rankdata(a)
    rb = rankdata(b)
    ra -= ra.mean()
    rb -= rb.mean()
    den = np.sqrt(np.dot(ra, ra) * np.dot(rb, rb))
    if den == 0:
        return np.nan
    return float(np.dot(ra, rb) / den)


def kernel_center(sigma: float) -> float:
    rad = int(np.ceil(4.0 * sigma)) + 1
    g = np.zeros((2 * rad + 1, 2 * rad + 1))
    g[rad, rad] = 1.0
    f = gaussian_filter(g, sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
    return float(f[rad, rad])


def section_rho(y: np.ndarray, x: np.ndarray, codes: np.ndarray) -> np.ndarray:
    out = np.full(len(SAMPLES), np.nan)
    for s in range(len(SAMPLES)):
        m = codes == s
        if int(m.sum()) < 30:
            continue
        out[s] = spearman(y[m], x[m])
    return out


def pack_summaries(rho: np.ndarray, n: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(rho) & (n > 0)
    equal = float(np.mean(rho[ok])) if ok.any() else np.nan
    nwt = float(np.sum(rho[ok] * n[ok]) / np.sum(n[ok])) if ok.any() else np.nan
    pm = []
    for p in PATIENTS:
        m = (SEC_PAT == p) & ok
        if m.any():
            pm.append(float(np.mean(rho[m])))
    patient = float(np.mean(pm)) if pm else np.nan
    return {
        "equal": equal,
        "n": nwt,
        "patient": patient,
        "n_neg": int(np.sum(rho < 0)),
        "n_finite": int(np.sum(np.isfinite(rho))),
        "rho": rho,
    }


def build_fovs(df: pd.DataFrame) -> list[dict]:
    fovs = []
    for (sample, fov), g in df.groupby(["sample", "fov"], sort=False):
        xy = np.column_stack(
            [g["x_um"].to_numpy(np.float64), g["y_um"].to_numpy(np.float64)]
        )
        types = g["cell_type"].to_numpy()
        cl_all = g["CLDN4"].to_numpy(np.float64)
        is_t = np.fromiter((str(t).startswith("tumor") for t in types), dtype=bool, count=len(types))
        if int(is_t.sum()) < 40 or len(g) < 40:
            continue
        tumor_xy = xy[is_t]
        cl = cl_all[is_t]
        xmin, ymin = xy.min(axis=0)
        xmax, ymax = xy.max(axis=0)
        nx = max(1, int(np.ceil(max(float(xmax - xmin), PIXEL) / PIXEL)))
        ny = max(1, int(np.ceil(max(float(ymax - ymin), PIXEL) / PIXEL)))
        ix = np.clip(np.floor((tumor_xy[:, 0] - xmin) / PIXEL).astype(np.int32), 0, nx - 1)
        iy = np.clip(np.floor((tumor_xy[:, 1] - ymin) / PIXEL).astype(np.int32), 0, ny - 1)
        coords = np.vstack(
            [
                np.clip((tumor_xy[:, 1] - ymin) / PIXEL - 0.5, 0, ny - 1),
                np.clip((tumor_xy[:, 0] - xmin) / PIXEL - 0.5, 0, nx - 1),
            ]
        )
        cd8 = np.isin(types, CD8_LABELS)
        cyto = cd8 | (types == NK_LABEL)

        def count_grid(mask: np.ndarray) -> np.ndarray:
            grid = np.zeros((ny, nx), dtype=np.float64)
            if int(mask.sum()) == 0:
                return grid
            sx = np.clip(np.floor((xy[mask, 0] - xmin) / PIXEL).astype(np.int32), 0, nx - 1)
            sy = np.clip(np.floor((xy[mask, 1] - ymin) / PIXEL).astype(np.int32), 0, ny - 1)
            np.add.at(grid, (sy, sx), 1.0)
            return grid

        cd8_grid = count_grid(cd8)
        cyto_grid = count_grid(cyto)
        tum_grid = np.zeros((ny, nx), dtype=np.float64)
        np.add.at(tum_grid, (iy, ix), 1.0)
        gauss = {}
        for h in HS:
            sigma = h / PIXEL
            den = gaussian_filter(np.ones((ny, nx)), sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
            den_q = map_coordinates(den, coords, order=1, mode="nearest")
            keep = den_q >= MASK
            if int(keep.sum()) < 25:
                continue
            k0 = kernel_center(sigma)
            num_cd8 = gaussian_filter(cd8_grid, sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
            num_cy = gaussian_filter(cyto_grid, sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
            num_n = gaussian_filter(tum_grid, sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
            cd8_q = map_coordinates(num_cd8, coords, order=1, mode="nearest") / np.clip(den_q, 1e-6, None)
            cy_q = map_coordinates(num_cy, coords, order=1, mode="nearest") / np.clip(den_q, 1e-6, None)
            nn_q = map_coordinates(num_n, coords, order=1, mode="nearest")
            nn_loo = np.maximum(nn_q - k0, 0.0)
            gauss[h] = {
                "keep": np.flatnonzero(keep).astype(np.int32),
                "sigma": sigma,
                "k0": k0,
                "cd8": cd8_q[keep].astype(np.float64),
                "cyto": cy_q[keep].astype(np.float64),
                "nn_loo": nn_loo[keep].astype(np.float64),
            }
        balls = {}
        # Border inset so the disk is mostly inside the FOV rectangle.
        x0, y0 = float(tumor_xy[:, 0].min()), float(tumor_xy[:, 1].min())
        x1, y1 = float(tumor_xy[:, 0].max()), float(tumor_xy[:, 1].max())
        inset = np.minimum(np.minimum(tumor_xy[:, 0] - x0, x1 - tumor_xy[:, 0]), np.minimum(tumor_xy[:, 1] - y0, y1 - tumor_xy[:, 1]))
        tree_t = cKDTree(tumor_xy)
        for r, src_xy, key in (
            (None, xy[cd8], "cd8"),
            (None, xy[cyto], "cyto"),
        ):
            pass
        srcs = {
            "cd8": xy[cd8],
            "cyto": xy[cyto],
        }
        trees = {k: cKDTree(v) if len(v) else None for k, v in srcs.items()}
        for r in BALLS:
            inside = inset >= (0.5 * r)
            if int(inside.sum()) < 25:
                continue
            idx = np.flatnonzero(inside).astype(np.int32)
            counts = {}
            for name, tree in trees.items():
                if tree is None:
                    counts[name] = np.zeros(idx.size, dtype=np.float64)
                else:
                    counts[name] = np.array(
                        [len(hit) for hit in tree.query_ball_point(tumor_xy[idx], r)],
                        dtype=np.float64,
                    )
            # Leave-one-out mean CLDN4 inside the disk, as a sparse matvec.
            neigh = tree_t.query_ball_point(tumor_xy[idx], r)
            rows, cols, data = [], [], []
            for i, hits in enumerate(neigh):
                for j in hits:
                    if j == int(idx[i]):
                        continue
                    rows.append(i)
                    cols.append(int(j))
                    data.append(1.0)
            if rows:
                mat = sparse.csr_matrix((data, (rows, cols)), shape=(idx.size, cl.size))
            else:
                mat = sparse.csr_matrix((idx.size, cl.size))
            row_n = np.asarray(mat.sum(axis=1)).ravel()
            balls[r] = {"idx": idx, "cd8": counts["cd8"], "cyto": counts["cyto"], "mat": mat, "row_n": row_n}
        fovs.append(
            {
                "sample": sample,
                "fov": int(fov),
                "sec": SAMPLES.index(sample),
                "cl": cl,
                "iy": iy,
                "ix": ix,
                "ny": ny,
                "nx": nx,
                "coords": coords,
                "gauss": gauss,
                "balls": balls,
            }
        )
        if len(fovs) % 40 == 0:
            print(f"  built {len(fovs)} FOVs", flush=True)
    return fovs


def nw_from_cl(fov: dict, h: float, cl: np.ndarray) -> np.ndarray:
    g = fov["gauss"][h]
    grid = np.zeros((fov["ny"], fov["nx"]), dtype=np.float64)
    np.add.at(grid, (fov["iy"], fov["ix"]), cl)
    num = gaussian_filter(grid, sigma=g["sigma"], mode="constant", cval=0.0, truncate=4.0)
    nq = map_coordinates(num, fov["coords"], order=1, mode="nearest")
    keep = g["keep"]
    nq_loo = np.maximum(nq[keep] - g["k0"] * cl[keep], 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        nw = nq_loo / g["nn_loo"]
    nw[~np.isfinite(nw)] = np.nan
    return nw


def ball_nw(fov: dict, r: float, cl: np.ndarray) -> np.ndarray:
    b = fov["balls"][r]
    num = b["mat"].dot(cl)
    with np.errstate(divide="ignore", invalid="ignore"):
        nw = num / b["row_n"]
    nw[~np.isfinite(nw)] = np.nan
    return nw


def evaluate(fovs: list[dict], cl_of) -> dict[str, dict]:
    """cl_of(fov) -> CLDN4 vector aligned to fov['cl']."""
    out = {}
    # Gaussian specs
    for h in HS:
        buckets = {key: [] for key in ("sec", "own", "nw", "cd8", "cyto")}
        for fov in fovs:
            if h not in fov["gauss"]:
                continue
            g = fov["gauss"][h]
            keep = g["keep"]
            cl = cl_of(fov)
            nw = nw_from_cl(fov, h, cl)
            n = keep.size
            buckets["sec"].append(np.full(n, fov["sec"], dtype=np.int16))
            buckets["own"].append(cl[keep])
            buckets["nw"].append(nw)
            buckets["cd8"].append(g["cd8"])
            buckets["cyto"].append(g["cyto"])
        if not buckets["sec"]:
            continue
        sec = np.concatenate(buckets["sec"])
        own = np.concatenate(buckets["own"])
        nw = np.concatenate(buckets["nw"])
        n_by = np.array([(sec == s).sum() for s in range(len(SAMPLES))], dtype=float)
        for src in ("cd8", "cyto"):
            x = np.concatenate(buckets[src])
            for resp, y in (("own", own), ("nw", nw)):
                rho = section_rho(y, x, sec)
                out[f"gauss_{int(h)}_{src}_{resp}"] = pack_summaries(rho, n_by)
    for r in BALLS:
        buckets = {key: [] for key in ("sec", "own", "nw", "cd8", "cyto")}
        for fov in fovs:
            if r not in fov["balls"]:
                continue
            b = fov["balls"][r]
            cl = cl_of(fov)
            nw = ball_nw(fov, r, cl)
            n = b["idx"].size
            buckets["sec"].append(np.full(n, fov["sec"], dtype=np.int16))
            buckets["own"].append(cl[b["idx"]])
            buckets["nw"].append(nw)
            buckets["cd8"].append(b["cd8"])
            buckets["cyto"].append(b["cyto"])
        if not buckets["sec"]:
            continue
        sec = np.concatenate(buckets["sec"])
        own = np.concatenate(buckets["own"])
        nw = np.concatenate(buckets["nw"])
        n_by = np.array([(sec == s).sum() for s in range(len(SAMPLES))], dtype=float)
        for src in ("cd8", "cyto"):
            x = np.concatenate(buckets[src])
            for resp, y in (("own", own), ("nw", nw)):
                rho = section_rho(y, x, sec)
                out[f"disk_{int(r)}_{src}_{resp}"] = pack_summaries(rho, n_by)
    return out


def spec_table(observed: dict[str, dict], nulls: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    rows = []
    for key, pack in observed.items():
        kind, scale, src, resp = key.split("_")
        for weight in ("equal", "n", "patient"):
            stat = pack[weight]
            null = nulls[key][weight]
            finite = null[np.isfinite(null)]
            if np.isfinite(stat) and finite.size:
                p = float((1 + np.sum(np.abs(finite) >= abs(stat))) / (finite.size + 1))
            else:
                p = np.nan
            rows.append(
                {
                    "spec": key,
                    "density": kind,
                    "scale_um": int(scale),
                    "source": src,
                    "response": resp,
                    "weight": weight,
                    "stat": stat,
                    "p_two": p,
                    "n_neg": pack["n_neg"],
                    "n_finite": pack["n_finite"],
                    **{SAMPLES[i]: pack["rho"][i] for i in range(len(SAMPLES))},
                }
            )
    return pd.DataFrame(rows)


def choose_winner(table: pd.DataFrame) -> pd.Series:
    """Most negative summary with permutation p < 0.05. Pure CD8 wins ties."""
    ok = table[np.isfinite(table["stat"]) & (table["p_two"] < 0.05) & (table["stat"] < 0)].copy()
    if ok.empty:
        ok = table[np.isfinite(table["stat"])].copy()
    ok["cd8_first"] = ok["source"].eq("cd8").astype(int)
    ok = ok.sort_values(["stat", "cd8_first", "n_neg"], ascending=[True, False, False])
    return ok.iloc[0]


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


def plot_sweep(table: pd.DataFrame) -> None:
    sub = table[(table["weight"] == "equal") & (table["source"] == "cd8")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), sharey=True, constrained_layout=True)
    for ax, density, title in (
        (axes[0], "gauss", "Gaussian CD8 field"),
        (axes[1], "disk", "Disk CD8 count"),
    ):
        part = sub[sub["density"] == density]
        for resp, color, label in (
            ("own", "#4C78A8", "Cell CLDN4"),
            ("nw", "#E45756", "CLDN4 field"),
        ):
            p = part[part["response"] == resp].sort_values("scale_um")
            ax.plot(p["scale_um"], p["stat"], "-o", color=color, label=label)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.axhline(-0.042, color="0.7", lw=0.8, ls="--")
        ax.set_xlabel("Scale (µm)")
        ax.set_title(title)
    axes[0].set_ylabel("Unweighted mean of 8 section Spearman ρ")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Author CD8. Dashed line: pre-specified 50 µm cell-level mean (−0.042)", fontsize=11)
    savefig(fig, "sweep_unweighted_mean_rho")


def plot_winner(row: pd.Series) -> None:
    rhos = [row[s] for s in SAMPLES]
    colors = ["#E45756" if np.isfinite(v) and v < 0 else "#4C78A8" for v in rhos]
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    y = np.arange(len(SAMPLES))
    ax.axvline(0, color="0.35", lw=0.8)
    ax.barh(y, rhos, color=colors, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(SAMPLES)
    ax.invert_yaxis()
    ax.set_xlabel("Within-section Spearman ρ")
    ax.set_title(
        f"{row['density']} {int(row['scale_um'])} µm  {row['source']}  {row['response']}\n"
        f"{row['weight']} summary = {row['stat']:+.3f}   permutation p = {row['p_two']:.3f}"
    )
    savefig(fig, "sweep_winner_sections")


def main() -> None:
    apply_style()
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    print("loading", flush=True)
    df = pd.read_csv(DATA)
    df["x_um"] = df["x"].to_numpy(np.float64) * 1000.0
    df["y_um"] = df["y"].to_numpy(np.float64) * 1000.0
    print("building FOV kernels", flush=True)
    fovs = build_fovs(df)
    print(f"{len(fovs)} FOVs", flush=True)
    observed = evaluate(fovs, lambda fov: fov["cl"])
    cal = observed["gauss_50_cd8_own"]["equal"]
    print(f"calibration gauss 50 µm own-CLDN4 vs CD8 equal-mean ρ = {cal:+.4f}", flush=True)
    if not np.isfinite(cal) or abs(cal - (-0.042)) > 0.02:
        raise RuntimeError(f"calibration drifted from the published 50 µm mean: {cal}")
    keys = list(observed)
    nulls = {k: {w: np.full(N_PERM, np.nan) for w in ("equal", "n", "patient")} for k in keys}
    ids = {id(fov): i for i, fov in enumerate(fovs)}
    rng = np.random.default_rng(SEED)
    for b in range(N_PERM):
        shuffled = [rng.permutation(fov["cl"]) for fov in fovs]

        def cl_of(fov, _s=shuffled, _ids=ids):
            return _s[_ids[id(fov)]]

        stats = evaluate(fovs, cl_of)
        for k in keys:
            for w in ("equal", "n", "patient"):
                nulls[k][w][b] = stats[k][w]
        if (b + 1) % 25 == 0:
            print(f"  perm {b + 1}/{N_PERM}", flush=True)

    table = spec_table(observed, nulls)
    table.to_csv(TAB / "sweep_specs.csv", index=False)
    winner = choose_winner(table)
    # Selection-adjusted p: most negative summary on each permutation.
    obs_mat = []
    null_mat = []
    for k in keys:
        for w in ("equal", "n", "patient"):
            obs_mat.append(observed[k][w])
            null_mat.append(nulls[k][w])
    obs_min = float(np.nanmin(np.array(obs_mat, dtype=float)))
    null_min = np.nanmin(np.vstack(null_mat), axis=0)
    p_sel = float((1 + np.sum(null_min <= obs_min)) / (np.isfinite(null_min).sum() + 1))
    cd8_ok = table[(table["source"] == "cd8") & (table["p_two"] < 0.05) & (table["stat"] < 0)].sort_values("stat")
    if cd8_ok.empty:
        cd8_ok = table[(table["source"] == "cd8") & np.isfinite(table["stat"])].sort_values("stat")
    best_cd8 = cd8_ok.iloc[0]
    payload = {
        "n_perm": N_PERM,
        "seed": SEED,
        "n_specs": int(len(table)),
        "calibration_gauss50_cd8_own_equal": cal,
        "winner": {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in winner.to_dict().items()},
        "best_author_cd8": {
            k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in best_cd8.to_dict().items()
        },
        "selection_min_stat": obs_min,
        "selection_p": p_sel,
        "null_min_mean": float(np.nanmean(null_min)),
        "null_min_p05": float(np.nanpercentile(null_min, 5)),
    }
    # numpy types
    def _clean(obj):
        if isinstance(obj, dict):
            return {str(k): _clean(v) for k, v in obj.items()}
        if isinstance(obj, (np.floating, float)):
            v = float(obj)
            return None if not np.isfinite(v) else v
        if isinstance(obj, (np.integer, int)):
            return int(obj)
        return obj

    (TAB / "sweep_summary.json").write_text(json.dumps(_clean(payload), indent=2) + "\n")
    plot_sweep(table)
    plot_winner(winner)
    print("WINNER", winner["spec"], winner["weight"], round(float(winner["stat"]), 4), "p", winner["p_two"], flush=True)
    print("BEST CD8", best_cd8["spec"], best_cd8["weight"], round(float(best_cd8["stat"]), 4), "p", best_cd8["p_two"], flush=True)
    print("selection min", round(obs_min, 4), "p_sel", p_sel, flush=True)


if __name__ == "__main__":
    main()
