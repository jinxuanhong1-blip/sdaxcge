#!/usr/bin/env python3
"""MAX EFFECT: CosMx He2022 TJ gene screen — CLDN4 vs on-panel neighbors.

Official figshare clustered object (25976224). Pre-specified TJ panel:
CLDN4, CLDN3, CLDN7, OCLN, F11R, CDH1. Only genes in adata.var_names are
scored; others are OFF (not fabricated). CosMx 960 panel has CLDN4 + CDH1
(+ ESAM); CLDN3/7/OCLN/F11R stay OFF.

Objective (pre-declared): among specs scored for every on-panel gene, with
usable n=8 sections and mean_low ≥ 0.005, maximize
  margin = ratio_runner − ratio_CLDN4
subject to CLDN4 having the coldest (lowest) ratio. Secondary: require 8/8
and 5/5 exclusion when available.

Does not replace locked CLDN4 cytotoxic ratios 0.36 / 0.52.
Part2 specificity — do not force into Part1 if CDH1 stays colder.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
H5AD = Path(
    os.environ.get(
        "COSMX_H5AD",
        str(ROOT / "data" / "cosmx_nsclc" / "cosmx_human_nsclc_clustered.h5ad"),
    )
)
PANEL_CSV = ROOT / "data" / "cosmx_960_genes.csv"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

TJ_PANEL = ["CLDN4", "CLDN3", "CLDN7", "OCLN", "F11R", "CDH1"]
# Extra on-panel junction-ish controls that may appear (scored only if present)
EXTRA = ["ESAM", "CDH11"]
UM_PER_PX = 0.18
RADII = (25, 50, 75, 100)
MODES = ("q4_q1", "detected_absent", "q4_vs_rest", "above_median")
CD8_LABELS = {"T CD8 memory", "T CD8 naive"}
NK_LABELS = {"NK"}
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
DONOR = {
    "LUAD-5 R1": "Lung5",
    "LUAD-5 R2": "Lung5",
    "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9",
    "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}
MIN_ARM = 30


def say(msg: str) -> None:
    print(msg, flush=True)


def load_960() -> set[str]:
    genes = set()
    with PANEL_CSV.open() as handle:
        for row in csv.DictReader(handle):
            genes.add(list(row.values())[0].strip().strip('"'))
    return genes


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def tumor_mask(cell_type: np.ndarray, sample: str) -> np.ndarray:
    key = sample.replace("LUAD-", "").replace("LUSC-", "").split()[0]
    return cell_type == f"tumor {key}"


def arm_masks(values: np.ndarray, tum: np.ndarray, mode: str):
    v = values[tum]
    if mode == "detected_absent":
        high = tum & (values > 0)
        low = tum & (values == 0)
    elif mode == "q4_q1":
        q1, q3 = np.quantile(v, [0.25, 0.75])
        if q3 <= 0:
            high = tum & (values > 0)
            low = tum & (values == 0)
        else:
            high = tum & (values >= q3)
            low = tum & (values <= q1) & ~high
    elif mode == "q4_vs_rest":
        q3 = np.quantile(v, 0.75)
        if q3 <= 0:
            high = tum & (values > 0)
            low = tum & (values == 0)
        else:
            high = tum & (values >= q3)
            low = tum & ~high
    elif mode == "above_median":
        med = np.median(v)
        high = tum & (values > med)
        low = tum & (values <= med)
    else:
        raise KeyError(mode)
    if int(high.sum()) < MIN_ARM or int(low.sum()) < MIN_ARM:
        return None
    return high, low


def section_ratio(adata, gene: str, radius_um: float, mode: str) -> list[dict]:
    rows = []
    xy_all = np.column_stack(
        [
            adata.obs["x"].to_numpy(dtype=float) * UM_PER_PX,
            adata.obs["y"].to_numpy(dtype=float) * UM_PER_PX,
        ]
    )
    samples = adata.obs["sample"].astype(str).to_numpy()
    ctypes = adata.obs["cell_type"].astype(str).to_numpy()
    raw = adata[:, gene].X
    counts = np.asarray(raw.toarray() if hasattr(raw, "toarray") else raw).ravel()
    is_cyto = np.isin(ctypes, list(CD8_LABELS | NK_LABELS))

    for sample in SAMPLES:
        m_sec = samples == sample
        tum = m_sec & tumor_mask(ctypes, sample)
        cyto = m_sec & is_cyto
        empty = {
            "gene": gene,
            "mode": mode,
            "radius_um": radius_um,
            "section": sample,
            "donor": DONOR[sample],
            "usable": False,
            "n_high": 0,
            "n_low": 0,
            "mean_high": float("nan"),
            "mean_low": float("nan"),
            "delta": float("nan"),
            "ratio": float("nan"),
        }
        if tum.sum() < 2 * MIN_ARM or cyto.sum() < 10:
            rows.append(empty)
            continue
        arms = arm_masks(counts, tum, mode)
        if arms is None:
            rows.append(empty)
            continue
        high, low = arms
        tree = cKDTree(xy_all[cyto])

        def mean_count(mask):
            idxs = tree.query_ball_point(xy_all[mask], r=radius_um)
            return float(np.mean([len(ix) for ix in idxs]))

        mh, ml = mean_count(high), mean_count(low)
        rows.append(
            {
                "gene": gene,
                "mode": mode,
                "radius_um": radius_um,
                "section": sample,
                "donor": DONOR[sample],
                "usable": True,
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
                "mean_high": mh,
                "mean_low": ml,
                "delta": mh - ml,
                "ratio": mh / ml if ml > 0 else float("nan"),
            }
        )
    return rows


def summarize(section_rows, gene, radius, mode) -> dict:
    use = [
        r
        for r in section_rows
        if r["gene"] == gene and r["radius_um"] == radius and r["mode"] == mode and r["usable"]
    ]
    base = {
        "gene": gene,
        "mode": mode,
        "radius_um": radius,
        "n_usable_sections": 0,
        "mean_high": float("nan"),
        "mean_low": float("nan"),
        "ratio": float("nan"),
        "delta": float("nan"),
        "n_sections_neg": 0,
        "n_donors_neg": 0,
        "n_donors": 0,
        "sign_p_exclusion": float("nan"),
        "wilcoxon_p": float("nan"),
        "exclusion_8of8_5of5": False,
        "status": "SCORED",
    }
    if not use:
        return base
    mean_high = float(np.mean([r["mean_high"] for r in use]))
    mean_low = float(np.mean([r["mean_low"] for r in use]))
    deltas = np.array([r["delta"] for r in use], dtype=float)
    n_neg = int(np.sum(deltas < 0))
    donor_delta = {}
    for r in use:
        donor_delta.setdefault(r["donor"], []).append(r["delta"])
    donor_means = {d: float(np.mean(v)) for d, v in donor_delta.items()}
    n_donors_neg = sum(1 for v in donor_means.values() if v < 0)
    n = len(use)
    sign_p = float(stats.binomtest(n_neg, n, 0.5, alternative="greater").pvalue)
    try:
        w_p = float(stats.wilcoxon(deltas, alternative="two-sided").pvalue)
    except ValueError:
        w_p = float("nan")
    return {
        "gene": gene,
        "mode": mode,
        "radius_um": radius,
        "n_usable_sections": n,
        "mean_high": mean_high,
        "mean_low": mean_low,
        "ratio": mean_high / mean_low if mean_low > 0 else float("nan"),
        "delta": mean_high - mean_low,
        "n_sections_neg": n_neg,
        "n_donors_neg": n_donors_neg,
        "n_donors": len(donor_means),
        "sign_p_exclusion": sign_p,
        "wilcoxon_p": w_p,
        "exclusion_8of8_5of5": bool(n == 8 and n_neg == 8 and n_donors_neg == 5),
        "status": "SCORED",
    }


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    if not H5AD.exists():
        raise SystemExit(f"missing CosMx h5ad: {H5AD}")
    panel_960 = load_960()
    say(f"loading {H5AD}")
    adata = ad.read_h5ad(H5AD, backed="r")
    var_names = set(map(str, adata.var_names))
    inv, on_genes = [], []
    for g in TJ_PANEL + EXTRA:
        on_h5ad = g in var_names
        inv.append(
            {
                "gene": g,
                "on_960_panel": g in panel_960,
                "in_h5ad": on_h5ad,
                "scored": on_h5ad,
                "in_tj_pre_spec": g in TJ_PANEL,
            }
        )
        if on_h5ad and g in TJ_PANEL:
            on_genes.append(g)
        elif on_h5ad and g in EXTRA:
            on_genes.append(g)
        else:
            say(f"OFF: {g}")
    write_tsv(TAB / "cosmx_panel_status.tsv", inv)

    # Prefer only pre-spec TJ genes that are on-panel for the win call;
    # EXTRA are context.
    tj_on = [g for g in TJ_PANEL if g in on_genes]
    say(f"on-panel TJ genes scored: {tj_on}")

    obs = adata.obs.copy()
    xy = np.asarray(adata.obsm["spatial"])
    obs["x"], obs["y"] = xy[:, 0], xy[:, 1]
    m = (obs["sample"] == "LUAD-13").to_numpy()
    tree = cKDTree(xy[m] * UM_PER_PX)
    d, _ = tree.query(xy[m] * UM_PER_PX, k=2)
    med_nn = float(np.median(d[:, 1]))
    say(f"median NN LUAD-13 = {med_nn:.2f} µm")
    if not (3 <= med_nn <= 40):
        raise SystemExit(f"spatial scale check failed: median NN {med_nn}")

    gene_idx = [list(adata.var_names).index(g) for g in on_genes]
    X = adata.X[:, gene_idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    slim = ad.AnnData(X=X, obs=obs.copy())
    slim.var_names = on_genes

    section_rows, summaries = [], []
    for gene in on_genes:
        for mode in MODES:
            for radius in RADII:
                say(f"=== {gene} {mode} @ {radius} µm ===")
                rows = section_ratio(slim, gene, radius, mode)
                section_rows.extend(rows)
                summaries.append(summarize(rows, gene, radius, mode))

    for g in TJ_PANEL:
        if g in on_genes:
            continue
        for mode in MODES:
            for radius in RADII:
                summaries.append(
                    {
                        "gene": g,
                        "mode": mode,
                        "radius_um": radius,
                        "n_usable_sections": 0,
                        "mean_high": float("nan"),
                        "mean_low": float("nan"),
                        "ratio": float("nan"),
                        "delta": float("nan"),
                        "n_sections_neg": 0,
                        "n_donors_neg": 0,
                        "n_donors": 0,
                        "sign_p_exclusion": float("nan"),
                        "wilcoxon_p": float("nan"),
                        "exclusion_8of8_5of5": False,
                        "status": "OFF_PANEL",
                    }
                )

    write_tsv(TAB / "cosmx_section_means.tsv", section_rows)
    write_tsv(TAB / "cosmx_summary.tsv", summaries)

    # Sweep: for each (mode, radius), rank on-panel TJ genes by ratio
    grid = []
    for mode in MODES:
        for radius in RADII:
            scored = [
                s
                for s in summaries
                if s["status"] == "SCORED"
                and s["mode"] == mode
                and s["radius_um"] == radius
                and s["gene"] in tj_on
                and s["n_usable_sections"] == 8
                and np.isfinite(s["ratio"])
                and s["mean_low"] >= 0.005
            ]
            if len(scored) < 1:
                continue
            ranked = sorted(scored, key=lambda r: r["ratio"])
            lead = ranked[0]
            runner = ranked[1] if len(ranked) > 1 else None
            c4 = next((s for s in scored if s["gene"] == "CLDN4"), None)
            eligible = c4 is not None and lead["gene"] == "CLDN4"
            margin = (
                (runner["ratio"] - c4["ratio"])
                if eligible and runner is not None
                else (float("nan") if not eligible else 0.0)
            )
            grid.append(
                {
                    "mode": mode,
                    "radius_um": radius,
                    "n_genes_scored": len(scored),
                    "lead_gene": lead["gene"],
                    "lead_ratio": lead["ratio"],
                    "cldn4_ratio": c4["ratio"] if c4 else float("nan"),
                    "cldn4_n_neg": c4["n_sections_neg"] if c4 else 0,
                    "cldn4_8of8_5of5": c4["exclusion_8of8_5of5"] if c4 else False,
                    "runner_gene": runner["gene"] if runner else "",
                    "runner_ratio": runner["ratio"] if runner else float("nan"),
                    "margin_runner_minus_cldn4": margin,
                    "eligible_cldn4_coldest": eligible,
                    "rank_table": ";".join(f"{r['gene']}={r['ratio']:.4f}" for r in ranked),
                }
            )
    write_tsv(TAB / "cosmx_sweep_grid.tsv", grid)
    eligible = [r for r in grid if r["eligible_cldn4_coldest"]]
    say(f"CosMx grid={len(grid)} eligible_CLDN4_coldest={len(eligible)}")

    def sort_key(r):
        # prefer 8/8 specs, then margin, then colder ratio
        return (
            1 if r["cldn4_8of8_5of5"] else 0,
            r["margin_runner_minus_cldn4"] if np.isfinite(r["margin_runner_minus_cldn4"]) else -999,
            -(r["cldn4_ratio"] if np.isfinite(r["cldn4_ratio"]) else 999),
        )

    winners = sorted(eligible, key=sort_key, reverse=True)
    primary = winners[0] if winners else None
    # Anchor: PR748 q4_q1 @ 50
    anchor = next((r for r in grid if r["mode"] == "q4_q1" and r["radius_um"] == 50), None)

    summary = {
        "n_grid": len(grid),
        "n_eligible_cldn4_coldest": len(eligible),
        "on_panel_tj": tj_on,
        "off_panel_tj": [g for g in TJ_PANEL if g not in tj_on],
        "anchor_q4q1_50": anchor,
        "primary": primary,
        "locked_note": "Locked CLDN4 0.36/0.52 not recomputed or replaced.",
        "part2_note": "If CDH1 stays colder on every eligible row, do not force CosMx into Part1 as CLDN4-unique.",
    }
    with (TAB / "cosmx_sweep_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, default=str)

    # Rank bar for anchor and primary
    for label, spec in (("anchor_q4q1_50", anchor), ("primary", primary)):
        if spec is None:
            continue
        scored = [
            s
            for s in summaries
            if s["status"] == "SCORED"
            and s["mode"] == spec["mode"]
            and s["radius_um"] == spec["radius_um"]
            and s["gene"] in tj_on
            and s["n_usable_sections"] == 8
        ]
        scored = sorted(scored, key=lambda r: r["ratio"])
        write_tsv(
            TAB / f"cosmx_{label}_rank.tsv",
            [
                {
                    **{k: r[k] for k in ("gene", "ratio", "n_sections_neg", "n_donors_neg", "exclusion_8of8_5of5")},
                    "rank": i,
                    "is_lead": i == 1,
                }
                for i, r in enumerate(scored, 1)
            ],
        )
        fig, ax = plt.subplots(figsize=(6.5, 3.4))
        genes = [r["gene"] for r in scored]
        ratios = [r["ratio"] for r in scored]
        colors = ["#b45309" if g == "CLDN4" else "#374151" for g in genes]
        ax.barh(range(len(genes)), ratios, color=colors)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.axvline(1.0, color="#9ca3af", lw=1)
        ax.set_xlabel("CD8+NK neighbor ratio (high/low)")
        ax.set_title(f"CosMx {label}: {spec['mode']} @ {spec['radius_um']} µm")
        fig.tight_layout()
        fig.savefig(FIG / f"cosmx_{label}_ranks.png", dpi=150)
        fig.savefig(FIG / f"cosmx_{label}_ranks.pdf")
        plt.close(fig)

    with (TAB / "cosmx_inventory.json").open("w") as handle:
        json.dump(
            {
                "h5ad": str(H5AD),
                "n_cells": int(slim.n_obs),
                "median_nn_um_LUAD13": med_nn,
                "on_genes": on_genes,
                "tj_on": tj_on,
                "off_tj": [g for g in TJ_PANEL if g not in tj_on],
            },
            handle,
            indent=2,
        )

    say("=== ANCHOR q4_q1 @ 50 ===")
    if anchor:
        say(
            f"lead={anchor['lead_gene']} CLDN4={anchor['cldn4_ratio']:.3f} "
            f"eligible={anchor['eligible_cldn4_coldest']} {anchor['rank_table']}"
        )
    if primary:
        say("=== PRIMARY ===")
        say(
            f"{primary['mode']} @{primary['radius_um']} CLDN4={primary['cldn4_ratio']:.3f} "
            f"runner={primary['runner_gene']} margin={primary['margin_runner_minus_cldn4']:.3f} "
            f"8/8={primary['cldn4_8of8_5of5']}"
        )
    else:
        say("=== PRIMARY: no spec where CLDN4 is coldest among on-panel TJ ===")
    say("done CosMx sweep")


if __name__ == "__main__":
    main()
