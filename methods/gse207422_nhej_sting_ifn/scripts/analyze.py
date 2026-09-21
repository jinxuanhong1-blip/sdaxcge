#!/usr/bin/env python3
"""GSE207422 NHEJ / STING / IFN split on public UMI.

Gate: TACSTD2 (TROP2) and CLDN4 must both be rows in the matrix.
If either is absent, write a gate-stop table and do not invent scores.

Primary unit is the post-treatment patient. Within A3-malignant cells,
CLDN4 log1p(CP10k) Q4 vs Q1 is the splitter. TACSTD2 is a second splitter,
not a gate. NHEJ, STING, and IFN effectors are scored separately and do
not share genes.

A3-malignant-like = epithelial marker-argmax AND zero UMI for
SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. This is not Hu et al. CopyKAT.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (
    A3_NORMAL_LUNG,
    LINEAGES,
    MODULES,
    PRIMARY_MODULES,
    STING_ALIASES,
    requested_symbols,
)

HERE = Path(__file__).resolve().parents[1]

PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}
COLOR = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}

MIN_MAL = 20
MIN_TAIL = 8


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    m = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(m)
    out[order] = adj
    q[ok] = out
    return q


def log1p_cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    umi = np.asarray(umi, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    out = np.full(umi.shape, np.nan)
    ok = lib > 0
    out[ok] = np.log1p(umi[ok] / lib[ok] * 1e4)
    return out


def score_log1p(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float64)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score_log1p(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def resolve_symbol(expr: dict[str, np.ndarray], name: str) -> str | None:
    options = STING_ALIASES.get(name, (name,))
    present = [g for g in options if g in expr]
    if not present:
        return None
    return max(present, key=lambda g: int(expr[g].sum()))


def module_vector(
    expr: dict[str, np.ndarray], genes: list[str], lib: np.ndarray
) -> tuple[np.ndarray, list[str]]:
    used: list[str] = []
    cols = []
    for g in genes:
        sym = resolve_symbol(expr, g)
        if sym is None:
            continue
        used.append(sym)
        cols.append(log1p_cp10k(expr[sym], lib))
    if not cols:
        return np.full(len(lib), np.nan), used
    return np.mean(np.vstack(cols), axis=0), used


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def spearman_row(x, y, contrast: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"contrast": contrast, "n": n, "spearman_rho": np.nan, "spearman_p": np.nan, "note": "too_few_or_constant"}
    rho, p = stats.spearmanr(x, y)
    return {"contrast": contrast, "n": n, "spearman_rho": float(rho), "spearman_p": float(p), "note": ""}


def wilcoxon_paired(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    n = int(len(a))
    base = {
        "contrast": contrast,
        "n": n,
        "mean_high": float(np.mean(a)) if n else np.nan,
        "mean_low": float(np.mean(b)) if n else np.nan,
        "median_delta": np.nan,
        "mean_delta": np.nan,
        "wilcoxon_w": np.nan,
        "p_value": np.nan,
        "n_high_gt_low": int(np.sum(a > b)) if n else 0,
        "n_high_lt_low": int(np.sum(a < b)) if n else 0,
        "note": "too_few_samples" if n < 3 else "",
    }
    if n < 3:
        return base
    d = a - b
    base["median_delta"] = float(np.median(d))
    base["mean_delta"] = float(np.mean(d))
    if np.allclose(d, 0):
        base["wilcoxon_w"] = 0.0
        base["p_value"] = 1.0
        base["note"] = "all_deltas_zero"
        return base
    w, p = stats.wilcoxon(a, b, alternative="two-sided", zero_method="wilcox", method="exact")
    base["wilcoxon_w"] = float(w)
    base["p_value"] = float(p)
    return base


def partial_spearman(x, y, z) -> float:
    """Pearson correlation of rank-residuals. NaN if a vector is constant."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if len(x) < 10 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan

    def resid(a, cov):
        A = np.column_stack([np.ones(len(cov)), stats.rankdata(cov)])
        coef, *_ = np.linalg.lstsq(A, stats.rankdata(a), rcond=None)
        return stats.rankdata(a) - A @ coef

    rx, ry = resid(x, z), resid(y, z)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(stats.pearsonr(rx, ry)[0])


def savefig(fig, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_gate_stop(outdir: Path, n_genes: int, genes_present: list[str], missing: list[str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    row = {
        "gate": "TACSTD2_and_CLDN4",
        "TACSTD2_present": "TACSTD2" in genes_present,
        "CLDN4_present": "CLDN4" in genes_present,
        "n_genes_in_matrix": n_genes,
        "status": "stop",
        "reason": "TROP2/TACSTD2 or CLDN4 is not a row in the public UMI. NHEJ/STING/IFN split was not run.",
    }
    pd.DataFrame([row]).to_csv(outdir / "tables" / "gate.tsv", sep="\t", index=False)
    (outdir / "summary.json").write_text(json.dumps({"status": "gate_stop", **row, "missing_panel": missing}, indent=2))
    print("GATE STOP", row["reason"], flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    outdir = args.outdir
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_nhej_sting_ifn.npz", allow_pickle=True)
    cells = extracted["cells"]
    genes = [str(g) for g in extracted["genes"]]
    mat = extracted["mat"]
    total = extracted["total"].astype(np.int64)
    n_genes_matrix = int(extracted["n_genes_in_matrix"][0])
    missing = [str(x) for x in extracted["missing"]]
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)

    if "TACSTD2" not in expr or "CLDN4" not in expr:
        write_gate_stop(outdir, n_genes_matrix, genes, missing)
        return

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    is_epi = lineage == "epithelial"
    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_malig = is_epi & (a3_normal == 0)
    is_tnk = np.isin(lineage, ["T", "NK"])
    is_mye = lineage == "myeloid"

    cldn4 = log1p_cp10k(expr["CLDN4"], total)
    tacstd2 = log1p_cp10k(expr["TACSTD2"], total)

    module_scores: dict[str, np.ndarray] = {}
    module_genes: dict[str, list[str]] = {}
    for name, glist in MODULES.items():
        sc, used = module_vector(expr, glist, total)
        module_scores[name] = sc
        module_genes[name] = used

    # Primary modules need at least 3 genes. Otherwise the split is not interpretable.
    module_ok = {name: len(module_genes[name]) >= 3 for name in MODULES}

    meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta.to_csv(tabdir / "sample_metadata.tsv", sep="\t", index=False)
    meta_map = meta.set_index("Sample")

    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_epi": is_epi,
            "is_malig_a3": is_malig,
            "is_tnk": is_tnk,
            "is_myeloid": is_mye,
            "total_umi": total,
            "CLDN4": cldn4,
            "TACSTD2": tacstd2,
        }
    )
    for name, sc in module_scores.items():
        per_cell[name] = sc
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)
    per_cell["timing"] = per_cell["Sample"].map(lambda s: meta_map.loc[s, "timing"] if s in meta_map.index else np.nan)

    # Gene detection. Post cells, by compartment.
    det_rows = []
    post_mask = per_cell["timing"].eq("post").to_numpy()
    compartments = {
        "all_post": post_mask,
        "a3_post": post_mask & is_malig,
        "epi_post": post_mask & is_epi,
        "myeloid_post": post_mask & is_mye,
        "tnk_post": post_mask & is_tnk,
    }
    scored_symbols = []
    for name, used in module_genes.items():
        for sym in used:
            scored_symbols.append((name, sym))
    for sym in ("TACSTD2", "CLDN4"):
        scored_symbols.append(("focal", sym))
    for name, sym in scored_symbols:
        raw = expr[sym]
        for comp, mask in compartments.items():
            nn = int(mask.sum())
            pos = int((raw[mask] > 0).sum()) if nn else 0
            det_rows.append(
                {
                    "module": name,
                    "gene": sym,
                    "compartment": comp,
                    "n_cells": nn,
                    "n_detected": pos,
                    "pct_detected": (100.0 * pos / nn) if nn else np.nan,
                    "total_umi": int(raw[mask].sum()) if nn else 0,
                }
            )
    det = pd.DataFrame(det_rows)
    det.to_csv(tabdir / "gene_detection.tsv", sep="\t", index=False)

    present_rows = []
    for name, requested in MODULES.items():
        for g in requested:
            sym = resolve_symbol(expr, g)
            present_rows.append(
                {
                    "module": name,
                    "requested": g,
                    "used_symbol": sym if sym else "",
                    "present": sym is not None,
                    "n_umi_all": int(expr[sym].sum()) if sym else 0,
                    "pct_cells_all": float((expr[sym] > 0).mean() * 100) if sym else 0.0,
                }
            )
    pd.DataFrame(present_rows).to_csv(tabdir / "module_genes.tsv", sep="\t", index=False)

    splitters = ("CLDN4", "TACSTD2")
    compartments_split = (
        ("is_malig_a3", "a3"),
        ("is_epi", "epi"),
    )
    patient_rows = []
    paired_rows = []
    partial_rows = []

    for sid, g in per_cell.groupby("Sample", sort=True):
        group = PAPER_GROUP.get(sid, "NA")
        timing = str(g["timing"].iloc[0])
        rec = {
            "Sample": sid,
            "paper_group": group,
            "timing": timing,
            "n_cells": int(len(g)),
            "n_epithelial": int(g["is_epi"].sum()),
            "n_malig_a3": int(g["is_malig_a3"].sum()),
            "n_tnk": int(g["is_tnk"].sum()),
            "n_myeloid": int(g["is_myeloid"].sum()),
        }
        for flag, tag in (
            ("is_malig_a3", "a3"),
            ("is_epi", "epi"),
            ("is_tnk", "tnk"),
            ("is_myeloid", "mye"),
        ):
            sub = g.loc[g[flag]]
            rec[f"n_{tag}"] = int(len(sub))
            if len(sub) == 0:
                rec[f"mean_CLDN4_{tag}"] = np.nan
                rec[f"mean_TACSTD2_{tag}"] = np.nan
                rec[f"pct_CLDN4_{tag}"] = np.nan
                rec[f"pct_TACSTD2_{tag}"] = np.nan
                for name in MODULES:
                    rec[f"mean_{name}_{tag}"] = np.nan
                continue
            rec[f"mean_CLDN4_{tag}"] = float(sub["CLDN4"].mean())
            rec[f"mean_TACSTD2_{tag}"] = float(sub["TACSTD2"].mean())
            rec[f"pct_CLDN4_{tag}"] = float((sub["CLDN4"] > 0).mean() * 100)
            rec[f"pct_TACSTD2_{tag}"] = float((sub["TACSTD2"] > 0).mean() * 100)
            for name in MODULES:
                rec[f"mean_{name}_{tag}"] = float(sub[name].mean())
        patient_rows.append(rec)

        for flag, tag in compartments_split:
            sub = g.loc[g[flag]]
            for splitter in splitters:
                n_sub = int(len(sub))
                base = {
                    "Sample": sid,
                    "paper_group": group,
                    "timing": timing,
                    "compartment": tag,
                    "splitter": splitter,
                    "n_cells": n_sub,
                }
                if n_sub < MIN_MAL:
                    paired_rows.append({**base, "n_high": 0, "n_low": 0, "eligible": False, "reason": f"n_cells<{MIN_MAL}"})
                    continue
                q = sub[splitter]
                q75 = float(q.quantile(0.75))
                q25 = float(q.quantile(0.25))
                if not np.isfinite(q75) or q75 <= q25:
                    paired_rows.append({
                        **base, "n_high": 0, "n_low": 0, "eligible": False, "reason": "no_spread",
                        "q75": q75, "q25": q25,
                    })
                    continue
                hi = sub.loc[q >= q75]
                lo = sub.loc[q <= q25]
                if len(hi) < MIN_TAIL or len(lo) < MIN_TAIL:
                    paired_rows.append({
                        **base,
                        "n_high": int(len(hi)),
                        "n_low": int(len(lo)),
                        "eligible": False,
                        "reason": f"tail<{MIN_TAIL}",
                        "q75": q75,
                        "q25": q25,
                    })
                    continue
                prow = {
                    **base,
                    "n_high": int(len(hi)),
                    "n_low": int(len(lo)),
                    "eligible": True,
                    "reason": "",
                    "q75": q75,
                    "q25": q25,
                    "splitter_high": float(hi[splitter].mean()),
                    "splitter_low": float(lo[splitter].mean()),
                    "umi_high": float(hi["total_umi"].mean()),
                    "umi_low": float(lo["total_umi"].mean()),
                    "CLDN4_high": float(hi["CLDN4"].mean()),
                    "CLDN4_low": float(lo["CLDN4"].mean()),
                    "TACSTD2_high": float(hi["TACSTD2"].mean()),
                    "TACSTD2_low": float(lo["TACSTD2"].mean()),
                }
                for name in MODULES:
                    prow[f"{name}_high"] = float(hi[name].mean())
                    prow[f"{name}_low"] = float(lo[name].mean())
                    prow[f"{name}_delta"] = prow[f"{name}_high"] - prow[f"{name}_low"]
                # Same Q4/Q1 cells, after a within-patient linear residual on log1p(UMI).
                if tag == "a3":
                    x = np.log1p(sub["total_umi"].to_numpy(dtype=float))
                    A = np.column_stack([np.ones(len(sub)), x])
                    hi_pos = sub.index.isin(hi.index)
                    lo_pos = sub.index.isin(lo.index)
                    for name in list(PRIMARY_MODULES) + ["sting_no_irf3", "oxphos", "prolif"]:
                        y = sub[name].to_numpy(dtype=float)
                        if not np.isfinite(y).all():
                            prow[f"{name}_resid_delta"] = np.nan
                            continue
                        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
                        resid = y - A @ coef
                        prow[f"{name}_resid_delta"] = float(resid[hi_pos].mean() - resid[lo_pos].mean())
                paired_rows.append(prow)

                # Depth sensitivity: partial Spearman inside this patient's cells.
                if tag == "a3" and n_sub >= 30:
                    for name in PRIMARY_MODULES:
                        rho_raw = stats.spearmanr(sub[splitter], sub[name]).statistic
                        rho_part = partial_spearman(sub[splitter], sub[name], np.log1p(sub["total_umi"]))
                        partial_rows.append({
                            "Sample": sid,
                            "paper_group": group,
                            "timing": timing,
                            "splitter": splitter,
                            "module": name,
                            "n_cells": n_sub,
                            "spearman_raw": float(rho_raw) if np.isfinite(rho_raw) else np.nan,
                            "spearman_partial_lib": float(rho_part) if np.isfinite(rho_part) else np.nan,
                        })

    patients = pd.DataFrame(patient_rows)
    paired = pd.DataFrame(paired_rows)
    partial = pd.DataFrame(partial_rows)
    patients.to_csv(tabdir / "per_patient.tsv", sep="\t", index=False)
    paired.to_csv(tabdir / "paired_high_low.tsv", sep="\t", index=False)
    partial.to_csv(tabdir / "partial_spearman_cells.tsv", sep="\t", index=False)

    tests = []

    def add(kind: str, family: str, row: dict) -> None:
        tests.append({"kind": kind, "family": family, **row})

    # Paired tests.
    for splitter in splitters:
        for comp in ("a3", "epi"):
            sub = paired[
                (paired["splitter"] == splitter)
                & (paired["compartment"] == comp)
                & (paired["timing"] == "post")
                & (paired["eligible"])
            ]
            family = f"paired_{splitter}_{comp}"
            for name in list(MODULES) + ["umi"]:
                if name == "umi":
                    add("paired", family, wilcoxon_paired(sub["umi_high"], sub["umi_low"], f"{splitter}_{comp}_umi"))
                else:
                    if sub.empty or f"{name}_high" not in sub.columns:
                        add("paired", family, wilcoxon_paired([], [], f"{splitter}_{comp}_{name}"))
                    else:
                        add("paired", family, wilcoxon_paired(sub[f"{name}_high"], sub[f"{name}_low"], f"{splitter}_{comp}_{name}"))
            if comp == "a3" and len(sub):
                for name in list(PRIMARY_MODULES) + ["sting_no_irf3", "oxphos", "prolif"]:
                    col = f"{name}_resid_delta"
                    if col not in sub.columns:
                        continue
                    add(
                        "paired_resid",
                        family,
                        wilcoxon_paired(sub[col], np.zeros(len(sub)), f"{splitter}_{comp}_{name}_resid_umi"),
                    )
            # Cross-module delta correlations (the split question).
            if len(sub) >= 4:
                for a_name, b_name in (
                    ("nhej", "sting"),
                    ("sting", "ifn_effector"),
                    ("nhej", "ifn_effector"),
                    ("nhej", "hr"),
                    ("nhej", "prolif"),
                    ("ifn_effector", "ifn_isg40"),
                ):
                    add(
                        "delta_spearman",
                        family,
                        spearman_row(sub[f"{a_name}_delta"], sub[f"{b_name}_delta"], f"delta_{splitter}_{comp}_{a_name}_vs_{b_name}"),
                    )

    # Patient-mean Spearman. Post, A3 n>=20.
    post = patients[patients["timing"] == "post"].copy()
    post_a3 = post[post["n_malig_a3"] >= MIN_MAL].copy()
    for splitter in ("CLDN4", "TACSTD2"):
        for name in MODULES:
            add(
                "sample_spearman",
                f"a3_mean_{splitter}",
                spearman_row(post_a3[f"mean_{splitter}_a3"], post_a3[f"mean_{name}_a3"], f"post_a3_{splitter}_vs_{name}"),
            )
        add(
            "sample_spearman",
            f"a3_mean_{splitter}",
            spearman_row(post_a3[f"mean_CLDN4_a3"], post_a3["mean_TACSTD2_a3"], "post_a3_CLDN4_vs_TACSTD2"),
        )

    # Immune compartment secondary: malignant CLDN4/TACSTD2 vs immune-cell modules.
    for immune_tag, immune_n, min_immune in (("mye", "n_myeloid", 20), ("tnk", "n_tnk", 20)):
        sub = post[(post["n_malig_a3"] >= MIN_MAL) & (post[immune_n] >= min_immune)]
        for splitter in ("CLDN4", "TACSTD2"):
            for name in PRIMARY_MODULES:
                add(
                    "cross_compartment",
                    f"a3_{splitter}_vs_{immune_tag}",
                    spearman_row(
                        sub[f"mean_{splitter}_a3"],
                        sub[f"mean_{name}_{immune_tag}"],
                        f"post_a3_{splitter}_vs_{immune_tag}_{name}",
                    ),
                )

    # Depth sensitivity summary: Wilcoxon of within-patient partial rhos vs 0.
    if len(partial):
        post_ids = set(post["Sample"])
        for splitter in splitters:
            for name in PRIMARY_MODULES:
                chunk = partial[(partial["splitter"] == splitter) & (partial["module"] == name) & (partial["Sample"].isin(post_ids))]
                # Restrict to patients eligible on the paired primary rule.
                elig_ids = set(paired[
                    (paired["splitter"] == splitter)
                    & (paired["compartment"] == "a3")
                    & (paired["timing"] == "post")
                    & (paired["eligible"])
                ]["Sample"])
                chunk = chunk[chunk["Sample"].isin(elig_ids)]
                rhos = chunk["spearman_partial_lib"].to_numpy(dtype=float)
                raws = chunk["spearman_raw"].to_numpy(dtype=float)
                add(
                    "partial_vs_0",
                    f"partial_{splitter}",
                    wilcoxon_paired(rhos, np.zeros(len(rhos)), f"partial_{splitter}_{name}"),
                )
                add(
                    "cell_spearman_vs_0",
                    f"partial_{splitter}",
                    wilcoxon_paired(raws, np.zeros(len(raws)), f"cellrho_{splitter}_{name}"),
                )

    tests_df = pd.DataFrame(tests)
    # BH within each family of primary paired module tests (3 tests).
    tests_df["p_bh"] = np.nan
    for family, idx in tests_df.groupby("family").groups.items():
        block = tests_df.loc[idx]
        primary_idx = block.index[block["contrast"].apply(lambda c: any(c.endswith("_" + m) for m in PRIMARY_MODULES))]
        if len(primary_idx) == 0:
            continue
        tests_df.loc[primary_idx, "p_bh"] = bh(tests_df.loc[primary_idx, "p_value"].to_numpy())
    tests_df.to_csv(tabdir / "tests.tsv", sep="\t", index=False)

    # Sign pattern for the primary CLDN4 A3 split.
    prim = paired[
        (paired["splitter"] == "CLDN4")
        & (paired["compartment"] == "a3")
        & (paired["timing"] == "post")
        & (paired["eligible"])
    ].copy()
    sign_rows = []
    if len(prim):
        for _, r in prim.iterrows():
            signs = {}
            for name in PRIMARY_MODULES:
                d = r[f"{name}_delta"]
                signs[name] = "up" if d > 0 else ("down" if d < 0 else "tie")
            n_up = sum(v == "up" for v in signs.values())
            sign_rows.append({
                "Sample": r["Sample"],
                "paper_group": r["paper_group"],
                "n_cells": r["n_cells"],
                **{f"{k}_delta": r[f"{k}_delta"] for k in PRIMARY_MODULES},
                **{f"{k}_sign": signs[k] for k in PRIMARY_MODULES},
                "pattern": "all_up" if n_up == 3 else ("all_down" if n_up == 0 else "mixed"),
            })
    signs_df = pd.DataFrame(sign_rows)
    signs_df.to_csv(tabdir / "sign_pattern.tsv", sep="\t", index=False)

    # Lineage counts.
    lin = (
        per_cell.groupby(["Sample", "timing", "paper_group", "lineage"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    lin.to_csv(tabdir / "lineage_counts.tsv", sep="\t", index=False)

    gate = pd.DataFrame([{
        "gate": "TACSTD2_and_CLDN4",
        "TACSTD2_present": True,
        "CLDN4_present": True,
        "TACSTD2_umi": int(expr["TACSTD2"].sum()),
        "CLDN4_umi": int(expr["CLDN4"].sum()),
        "TACSTD2_pct_cells": float((expr["TACSTD2"] > 0).mean() * 100),
        "CLDN4_pct_cells": float((expr["CLDN4"] > 0).mean() * 100),
        "n_genes_in_matrix": n_genes_matrix,
        "n_cells": n,
        "status": "run",
        "modules_ready": json.dumps({k: module_ok[k] for k in PRIMARY_MODULES}),
        "module_genes": json.dumps(module_genes),
    }])
    gate.to_csv(tabdir / "gate.tsv", sep="\t", index=False)

    # Figures
    _figures(prim, patients, det, figdir)

    summary = {
        "status": "run",
        "n_genes_in_matrix": n_genes_matrix,
        "n_cells": n,
        "n_epithelial": int(is_epi.sum()),
        "n_malig_a3": int(is_malig.sum()),
        "n_post_samples": int((patients["timing"] == "post").sum()),
        "n_post_a3_eligible_cldn4": int(len(prim)),
        "eligible_patients_cldn4": prim["Sample"].tolist() if len(prim) else [],
        "missing_panel_genes": missing,
        "module_genes_used": module_genes,
        "module_ready": module_ok,
        "min_malignant": MIN_MAL,
        "min_tail": MIN_TAIL,
        "primary": "post A3-malignant CLDN4 Q4 vs Q1; NHEJ, STING, IFN effector scored separately",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str), flush=True)
    print(tests_df[tests_df["family"].isin(["paired_CLDN4_a3", "paired_TACSTD2_a3"])][
        ["contrast", "n", "mean_high", "mean_low", "median_delta", "n_high_gt_low", "p_value", "p_bh"]
    ].to_string(index=False), flush=True)


def _figures(prim: pd.DataFrame, patients: pd.DataFrame, det: pd.DataFrame, figdir: Path) -> None:
    if len(prim):
        fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.6), sharey=False)
        titles = {"nhej": "NHEJ", "sting": "STING", "ifn_effector": "IFN effectors"}
        for ax, name in zip(axes, PRIMARY_MODULES):
            for _, r in prim.iterrows():
                c = COLOR.get(r["paper_group"], "#444")
                ax.plot([0, 1], [r[f"{name}_low"], r[f"{name}_high"]], color=c, lw=1.2, alpha=0.9)
                ax.scatter([0, 1], [r[f"{name}_low"], r[f"{name}_high"]], color=c, s=22, zorder=3)
            ax.set_xticks([0, 1], ["CLDN4 Q1", "CLDN4 Q4"])
            ax.set_title(titles[name])
            ax.set_xlim(-0.25, 1.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        handles = [
            plt.Line2D([0], [0], color=COLOR["MPR"], lw=1.5, label="MPR"),
            plt.Line2D([0], [0], color=COLOR["NMPR"], lw=1.5, label="NMPR"),
        ]
        axes[2].legend(handles=handles, frameon=False, loc="best")
        axes[0].set_ylabel("mean log1p(CP10k)")
        fig.suptitle(f"Same malignant cells, paired patients n={len(prim)}", y=1.02, fontsize=11)
        fig.tight_layout()
        savefig(fig, figdir / "fig_paired_cldn4_split")

        fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.3))
        pairs = (("nhej", "sting"), ("sting", "ifn_effector"), ("nhej", "ifn_effector"))
        for ax, (a, b) in zip(axes, pairs):
            for _, r in prim.iterrows():
                c = COLOR.get(r["paper_group"], "#444")
                ax.scatter(r[f"{a}_delta"], r[f"{b}_delta"], color=c, s=36, zorder=3)
                ax.text(r[f"{a}_delta"], r[f"{b}_delta"], r["Sample"].replace("BD_immune", "P"), fontsize=7, color=c)
            ax.axhline(0, color="#bbb", lw=0.6)
            ax.axvline(0, color="#bbb", lw=0.6)
            ax.set_xlabel(f"Δ {a}")
            ax.set_ylabel(f"Δ {b}")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        fig.tight_layout()
        savefig(fig, figdir / "fig_delta_split")

        fig, ax = plt.subplots(figsize=(7.4, 3.6))
        labels = [s.replace("BD_immune", "P") for s in prim["Sample"]]
        x = np.arange(len(prim))
        w = 0.25
        for i, (col, lab, color) in enumerate((
            ("nhej_resid_delta", "NHEJ", "#6b705c"),
            ("sting_resid_delta", "STING", "#2c6eaf"),
            ("ifn_effector_resid_delta", "IFN", "#d1495b"),
        )):
            ax.bar(x + (i - 1) * w, prim[col], width=w, label=lab, color=color)
        ax.axhline(0, color="#222", lw=0.6)
        ax.set_xticks(x, labels)
        ax.set_ylabel("Q4 − Q1 after UMI residual")
        ax.set_title("Within-patient residual on log1p(UMI)")
        ax.legend(frameon=False, ncol=3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        savefig(fig, figdir / "fig_resid_umi")

    post = patients[patients["timing"] == "post"]
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    order = post.sort_values("Sample")
    colors = [COLOR.get(g, "#444") for g in order["paper_group"]]
    ax.bar(order["Sample"].str.replace("BD_immune", "P"), order["n_malig_a3"], color=colors)
    ax.axhline(MIN_MAL, color="#222", ls="--", lw=0.8, label=f"floor n={MIN_MAL}")
    ax.set_ylabel("A3-malignant cells")
    ax.set_title("Post-treatment occupancy")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    savefig(fig, figdir / "fig_honest_n")

    # Detection of the proximal genes in malignant vs myeloid.
    focus_genes = ["XRCC6", "XRCC5", "PRKDC", "LIG4", "CGAS", "MB21D1", "STING1", "TMEM173", "TBK1", "IRF3", "ISG15", "MX1", "IFIT1", "TACSTD2", "CLDN4"]
    sub = det[det["gene"].isin(focus_genes) & det["compartment"].isin(["a3_post", "myeloid_post"])]
    if len(sub):
        genes_ord = [g for g in focus_genes if g in set(sub["gene"])]
        fig, ax = plt.subplots(figsize=(8.2, 3.6))
        x = np.arange(len(genes_ord))
        w = 0.38
        for i, (comp, label, color) in enumerate((
            ("a3_post", "A3-malignant", "#b08968"),
            ("myeloid_post", "Myeloid", "#2c6eaf"),
        )):
            vals = []
            for g in genes_ord:
                hit = sub[(sub["gene"] == g) & (sub["compartment"] == comp)]
                vals.append(float(hit["pct_detected"].iloc[0]) if len(hit) else 0.0)
            ax.bar(x + (i - 0.5) * w, vals, width=w, label=label, color=color)
        ax.set_xticks(x, genes_ord, rotation=60, ha="right")
        ax.set_ylabel("% cells detected")
        ax.set_title("Detection in post-treatment cells")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        savefig(fig, figdir / "fig_detection")


if __name__ == "__main__":
    main()
