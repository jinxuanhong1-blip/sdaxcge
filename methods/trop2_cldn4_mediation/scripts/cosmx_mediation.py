#!/usr/bin/env python3
"""He 2022 CosMx: does the TACSTD2–cytotoxic-neighbor association shrink after CLDN4?

figshare 25976224, cosmx_human_nsclc_clustered.h5ad. Malignant cells are the
author cell_type matched to the section (tumor 5/6/9/12/13). Cytotoxic
neighbors are CD8 or NK cells inside 50 µm and inside 100 µm, the two radii
in the locked exclusion summary. Coordinates are global section pixels times
0.18 µm/px. Sections are never placed in one tree.

The independent unit is the donor (5), not the cell. Section coefficients (8)
are reported because the locked exclusion call is 8/8 sections and 5/5 donors.
Cell-level p-values are not used.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from stats import (  # noqa: E402
    acme_point,
    attenuation_percent,
    partial_spearman,
    spearman,
    verdict_row,
)

H5AD = Path("/tmp/cosmx_dl/cosmx_human_nsclc_clustered.h5ad")
OUT = ROOT / "results" / "tables"
EXPECTED_BYTES = 2_755_776_882
UM_PER_PX = 0.18
RADII = (50, 100)
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


def _decode(arr) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]


def _decode_cats(f, name: str):
    cats = _decode(f[f"obs/{name}/categories"][:])
    codes = np.asarray(f[f"obs/{name}/codes"][:])
    return np.asarray(cats, dtype=object)[codes], cats


def _extract_columns(f, gene_to_col: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
    indptr = f["layers/counts/indptr"]
    indices = f["layers/counts/indices"]
    data = f["layers/counts/data"]
    out = {g: np.zeros(n_cells, dtype=np.float32) for g in gene_to_col}
    chunk = 25000
    for start in range(0, n_cells, chunk):
        end = min(n_cells, start + chunk)
        i0 = int(indptr[start])
        i1 = int(indptr[end])
        idx = np.asarray(indices[i0:i1])
        dat = np.asarray(data[i0:i1], dtype=np.float32)
        local = np.asarray(indptr[start:end + 1]) - i0
        nnz = np.diff(local)
        rows = np.repeat(np.arange(start, end), nnz)
        for gene, col in gene_to_col.items():
            hit = idx == col
            if np.any(hit):
                out[gene][rows[hit]] = dat[hit]
        if start == 0 or end == n_cells or (start // chunk) % 10 == 0:
            print(f"  counts rows {end}/{n_cells}", flush=True)
    return out


def load_malignant_frame() -> pd.DataFrame:
    import h5py

    if not H5AD.is_file() or H5AD.stat().st_size != EXPECTED_BYTES:
        raise SystemExit(f"h5ad missing or wrong size at {H5AD}")
    t0 = time.time()
    f = h5py.File(H5AD, "r")
    genes = _decode(f["var/_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}
    for g in ("TACSTD2", "CLDN4"):
        if g not in gmap:
            raise SystemExit(f"{g} is absent from this CosMx object. Present example: {genes[:8]}")
    sample, sample_cats = _decode_cats(f, "sample")
    patient, _ = _decode_cats(f, "patient")
    cell_type, type_cats = _decode_cats(f, "cell_type")
    fov = np.asarray(f["obs/fov"][:], dtype=np.int32)
    n_counts = np.asarray(f["obs/n_counts"][:], dtype=np.float64)
    xy = np.asarray(f["obsm/spatial"][:], dtype=np.float64) * UM_PER_PX
    print(f"obs loaded in {time.time()-t0:.1f}s; extracting TACSTD2 and CLDN4", flush=True)
    raw = _extract_columns(f, {"TACSTD2": gmap["TACSTD2"], "CLDN4": gmap["CLDN4"]}, len(sample))
    f.close()

    missing = [s for s in SAMPLES if s not in set(sample_cats)]
    if missing:
        raise SystemExit(f"missing samples {missing}; found {sample_cats}")
    cd8_labels = [t for t in type_cats if "CD8" in t]
    nk_labels = [t for t in type_cats if t == "NK" or str(t).startswith("NK")]
    if not cd8_labels or not nk_labels:
        raise SystemExit(f"CD8/NK labels missing: {type_cats}")

    lib = n_counts.copy()
    lib[lib <= 0] = np.nan
    tac = np.log1p(raw["TACSTD2"] / lib * 1e4)
    cld = np.log1p(raw["CLDN4"] / lib * 1e4)
    tac = np.nan_to_num(tac, nan=0.0).astype(np.float32)
    cld = np.nan_to_num(cld, nan=0.0).astype(np.float32)

    malignant = np.array(
        [cell_type[i] == TUMOR_LABEL[sample[i]] for i in range(len(sample))],
        dtype=bool,
    )
    is_cd8 = np.isin(cell_type, cd8_labels)
    is_nk = np.isin(cell_type, nk_labels)
    is_cyto = is_cd8 | is_nk

    pieces = []
    inventory = []
    for s in SAMPLES:
        m = sample == s
        xy_s = xy[m]
        rng = np.random.default_rng(0)
        take = rng.choice(xy_s.shape[0], size=min(6000, xy_s.shape[0]), replace=False)
        tree_all = cKDTree(xy_s)
        dist, _ = tree_all.query(xy_s[take], k=2, workers=4)
        med_nn = float(np.median(dist[:, 1]))
        if not (3.0 <= med_nn <= 40.0):
            raise SystemExit(f"{s}: median NN {med_nn:.3f} µm is outside 3–40; check 0.18 µm/px")
        mal = malignant[m]
        cy = is_cyto[m]
        tum_xy = xy_s[mal]
        ref_xy = xy_s[cy]
        counts = np.zeros((int(mal.sum()), len(RADII)), dtype=np.float32)
        if len(ref_xy) and len(tum_xy):
            tree = cKDTree(ref_xy)
            for b, radius in enumerate(RADII):
                counts[:, b] = tree.query_ball_point(
                    tum_xy, r=float(radius), return_length=True, workers=4
                )
        pat = str(patient[m][0])
        piece = pd.DataFrame({
            "sample": s,
            "patient": pat,
            "tacstd2": tac[m][mal].astype(np.float64),
            "cldn4": cld[m][mal].astype(np.float64),
            "tacstd2_pos": raw["TACSTD2"][m][mal] > 0,
            "cldn4_pos": raw["CLDN4"][m][mal] > 0,
            "y50": counts[:, 0].astype(np.float64),
            "y100": counts[:, 1].astype(np.float64),
        })
        pieces.append(piece)
        inventory.append({
            "sample": s,
            "patient": pat,
            "tumor_label": TUMOR_LABEL[s],
            "n_cells": int(m.sum()),
            "n_malignant": int(mal.sum()),
            "n_cd8nk": int(cy.sum()),
            "frac_tacstd2_pos": float(piece.tacstd2_pos.mean()),
            "frac_cldn4_pos": float(piece.cldn4_pos.mean()),
            "mean_tacstd2_log1p_cp10k": float(piece.tacstd2.mean()),
            "mean_cldn4_log1p_cp10k": float(piece.cldn4.mean()),
            "mean_cd8nk_50": float(piece.y50.mean()),
            "mean_cd8nk_100": float(piece.y100.mean()),
            "median_nn_um": med_nn,
        })
        print(
            f"{s}: mal={int(mal.sum())} cd8nk={int(cy.sum())} "
            f"TACSTD2%={100*piece.tacstd2_pos.mean():.1f} CLDN4%={100*piece.cldn4_pos.mean():.1f} "
            f"NN={med_nn:.2f}",
            flush=True,
        )
    cells = pd.concat(pieces, ignore_index=True)
    meta = {
        "cd8_labels": cd8_labels,
        "nk_labels": nk_labels,
        "n_malignant": int(len(cells)),
        "inventory": inventory,
        "expression": "log1p(count / n_counts * 10000) in author-matched malignant cells",
        "neighbor": "count of CD8 or NK cells within the radius; trees are within section",
    }
    return cells, meta


def fit_unit(x, m, y) -> dict:
    raw_x, p_x = spearman(x, y)
    raw_m, p_m = spearman(m, y)
    px = partial_spearman(x, y, m)
    pm = partial_spearman(m, y, x)
    med = acme_point(x, m, y)
    med_log = acme_point(x, m, np.log1p(y))
    if med is None or med_log is None:
        raise SystemExit("ACME undefined (zero variance)")
    return {
        "n_cells": int(len(x)),
        "rho_tacstd2": raw_x,
        "p_tacstd2_cells": p_x,
        "rho_cldn4": raw_m,
        "p_cldn4_cells": p_m,
        "rho_tacstd2_given_cldn4": px["rho"],
        "p_tacstd2_given_cldn4_cells": px["p"],
        "rho_tacstd2_given_cldn4_rankresid": px["rank_residual_rho"],
        "rho_cldn4_given_tacstd2": pm["rho"],
        "p_cldn4_given_tacstd2_cells": pm["p"],
        "rho_tacstd2_cldn4": px["rho_xz"],
        "attenuation_pct": attenuation_percent(raw_x, px["rho"]),
        "acme": med["acme"],
        "ade": med["ade"],
        "total_linear": med["total"],
        "a_x_to_m": med["a"],
        "b_m_to_y": med["b"],
        "acme_log1p_y": med_log["acme"],
    }


def _mean_record(frame: pd.DataFrame, keys: list[str]) -> dict:
    out = {}
    for k in keys:
        out[k] = float(frame[k].mean())
    out["n_units"] = int(len(frame))
    out["n_negative_tacstd2"] = int((frame.rho_tacstd2 < 0).sum())
    out["n_negative_cldn4"] = int((frame.rho_cldn4 < 0).sum())
    out["n_negative_partial_tacstd2"] = int((frame.rho_tacstd2_given_cldn4 < 0).sum())
    out["n_negative_partial_cldn4"] = int((frame.rho_cldn4_given_tacstd2 < 0).sum())
    out["n_negative_acme"] = int((frame.acme < 0).sum())
    out["n_attenuated"] = int((
        (frame.rho_tacstd2 < 0) & (frame.rho_tacstd2_given_cldn4 > frame.rho_tacstd2)
    ).sum())
    return out


def cluster_bootstrap(units: pd.DataFrame, cluster: np.ndarray, n_boot: int = 4000, seed: int = 4639) -> dict:
    """Resample clusters, then average the unit rows that belong to them."""
    clusters = np.unique(cluster)
    groups = {c: units.index[cluster == c].to_numpy() for c in clusters}
    rng = np.random.default_rng(seed)
    keys = [
        "rho_tacstd2", "rho_cldn4",
        "rho_tacstd2_given_cldn4", "rho_cldn4_given_tacstd2", "acme",
    ]
    acc = {k: np.empty(n_boot) for k in keys}
    for i in range(n_boot):
        draw = rng.choice(clusters, size=len(clusters), replace=True)
        idx = np.concatenate([groups[c] for c in draw])
        taken = units.loc[idx]
        for k in keys:
            acc[k][i] = float(taken[k].mean())
    out = {"n_clusters": int(len(clusters)), "n_boot": n_boot}
    for k in keys:
        lo, hi = np.quantile(acc[k], [0.025, 0.975])
        est = acc[k]
        p = 2.0 * min(float(np.mean(est <= 0)), float(np.mean(est >= 0)))
        out[f"{k}_lo"] = float(lo)
        out[f"{k}_hi"] = float(hi)
        out[f"{k}_p"] = float(min(p, 1.0))
    return out


def between_unit(means: pd.DataFrame, ycol: str) -> dict:
    """Spearman across sections of section-mean expression vs section-mean neighbors."""
    rec = fit_unit(means.mean_tacstd2_log1p_cp10k, means.mean_cldn4_log1p_cp10k, means[ycol])
    # n is sections, so the cell-level p from fit_unit is the section-level p. Rename.
    rec["p_tacstd2"] = rec.pop("p_tacstd2_cells")
    rec["p_cldn4"] = rec.pop("p_cldn4_cells")
    rec["p_tacstd2_given_cldn4"] = rec.pop("p_tacstd2_given_cldn4_cells")
    rec["p_cldn4_given_tacstd2"] = rec.pop("p_cldn4_given_tacstd2_cells")
    rec["n"] = rec.pop("n_cells")
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells, meta = load_malignant_frame()
    detail_rows = []
    ymap = {"cd8nk_50um": "y50", "cd8nk_100um": "y100"}
    for sample, g in cells.groupby("sample", sort=False):
        for readout, ycol in ymap.items():
            rec = fit_unit(g.tacstd2, g.cldn4, g[ycol])
            rec.update({
                "sample": sample,
                "patient": g.patient.iloc[0],
                "readout": readout,
                "unit": "section",
            })
            detail_rows.append(rec)
    for patient, g in cells.groupby("patient", sort=False):
        for readout, ycol in ymap.items():
            rec = fit_unit(g.tacstd2, g.cldn4, g[ycol])
            rec.update({
                "sample": "",
                "patient": patient,
                "readout": readout,
                "unit": "donor",
            })
            detail_rows.append(rec)
    detail = pd.DataFrame(detail_rows)

    summary_rows = []
    keys = [
        "rho_tacstd2", "rho_cldn4", "rho_tacstd2_given_cldn4",
        "rho_cldn4_given_tacstd2", "attenuation_pct", "acme",
        "acme_log1p_y", "ade", "total_linear", "a_x_to_m", "b_m_to_y",
        "rho_tacstd2_cldn4",
    ]
    # Section rows are resampled by donor so the three LUAD-5 sections and
    # the two LUAD-9 sections are not treated as independent patients.
    for unit_name in ("section", "donor"):
        sub = detail[detail.unit == unit_name].copy()
        for readout in ymap:
            block = sub[sub.readout == readout].reset_index(drop=True)
            means = _mean_record(block, keys)
            # attenuation of the averaged coefficients, not the average of attenuations
            means["attenuation_of_mean_rho"] = attenuation_percent(
                means["rho_tacstd2"], means["rho_tacstd2_given_cldn4"]
            )
            clusters = block["patient"].to_numpy()
            boot = cluster_bootstrap(block, clusters)
            call = verdict_row(
                means["rho_tacstd2"], boot["rho_tacstd2_lo"], boot["rho_tacstd2_hi"],
                means["rho_tacstd2_given_cldn4"], boot["rho_tacstd2_given_cldn4_p"],
                means["acme"], boot["acme_lo"], boot["acme_hi"],
                means["rho_cldn4_given_tacstd2"],
            )
            summary_rows.append({
                "dataset": "He2022_CosMx",
                "unit": unit_name,
                "readout": readout,
                "verdict": call,
                **means,
                **boot,
            })

    inv = pd.DataFrame(meta["inventory"])
    between_rows = []
    for readout, ycol in (("cd8nk_50um", "mean_cd8nk_50"), ("cd8nk_100um", "mean_cd8nk_100")):
        rec = between_unit(inv, ycol)
        from stats import spearman_ci
        lo_x, hi_x = spearman_ci(rec["rho_tacstd2"], rec["n"])
        # One fit on 8 section means has no separate ACME interval. Pass a
        # placeholder interval that includes 0 so this sensitivity cannot
        # be called support on the product term alone.
        call = verdict_row(
            rec["rho_tacstd2"], lo_x, hi_x,
            rec["rho_tacstd2_given_cldn4"], rec["p_tacstd2_given_cldn4"],
            rec["acme"], min(float(rec["acme"]), -1.0), max(float(rec["acme"]), 1.0),
            rec["rho_cldn4_given_tacstd2"],
        )
        rec.update({"dataset": "He2022_CosMx", "unit": "between_section_means", "readout": readout, "verdict": call})
        between_rows.append(rec)

    sec50 = next(
        r for r in summary_rows
        if r["unit"] == "section" and r["readout"] == "cd8nk_50um"
    )
    # Prior continuous sweep (section mean Spearman of CLDN4 vs CD8+NK at 50 µm).
    if abs(sec50["rho_cldn4"] - (-0.0369809674423623)) > 1e-6:
        raise SystemExit(f"CLDN4 continuous check failed: {sec50['rho_cldn4']}")

    detail.to_csv(OUT / "cosmx_unit_coefficients.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT / "cosmx_summary.tsv", sep="\t", index=False)
    inv.to_csv(OUT / "cosmx_section_inventory.tsv", sep="\t", index=False)
    pd.DataFrame(between_rows).to_csv(OUT / "cosmx_between_section.tsv", sep="\t", index=False)
    meta_public = {k: v for k, v in meta.items() if k != "inventory"}
    (OUT / "cosmx_meta.json").write_text(json.dumps(meta_public, indent=2))
    print(pd.DataFrame(summary_rows)[[
        "unit", "readout", "rho_tacstd2", "rho_cldn4",
        "rho_tacstd2_given_cldn4", "rho_cldn4_given_tacstd2",
        "attenuation_of_mean_rho", "acme", "verdict",
        "n_negative_tacstd2", "n_negative_acme", "n_units",
    ]].to_string(index=False))
    print("malignant cells", meta["n_malignant"])


if __name__ == "__main__":
    main()
