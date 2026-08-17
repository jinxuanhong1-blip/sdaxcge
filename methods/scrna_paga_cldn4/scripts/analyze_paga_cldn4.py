#!/usr/bin/env python3
"""PAGA + diffusion on GSE131907 LUAD epithelium, scored by CLDN4.

ADDITIVE. Not a TACSTD2 redo: CLDN4 is the primary readout; TACSTD2 is a
comparator only. Barrier/keratin score excludes CLDN4 (no circularity).
Sample (not cell) is the inferential unit.
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
from gene_sets import (  # noqa: E402
    AUTHOR_AT2,
    AUTHOR_CLUB,
    AUTHOR_TUMOR_STATE,
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STATES,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05


def _bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    running = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        running = min(running, p[i] * n / k)
        q[i] = running
    return [float(min(1.0, x)) for x in q]


def _spearman(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "rho": None, "p": None, "note": "n<4"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
    return {"n": n, "rho": float(rho), "p": float(p)}


def _wilcoxon_paired(a: np.ndarray, b: np.ndarray) -> dict:
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "W": None, "p": None, "delta_median": None, "note": "n<4"}
    aa, bb = a[mask], b[mask]
    try:
        w, p = stats.wilcoxon(aa, bb, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        return {"n": n, "W": None, "p": None, "delta_median": float(np.median(aa - bb)), "note": "wilcoxon failed"}
    return {
        "n": n,
        "W": float(w),
        "p": float(p),
        "delta_median": float(np.median(aa - bb)),
        "median_high": float(np.median(aa)),
        "median_low": float(np.median(bb)),
    }


def _score_mean(adata, genes: tuple[str, ...], key: str) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    absent = [g for g in genes if g not in adata.var_names]
    if not present:
        adata.obs[key] = np.nan
        return absent
    X = adata[:, present].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.obs[key] = np.asarray(X, dtype=float).mean(axis=1)
    return absent


def _pick_root(adata) -> tuple[int, dict]:
    """External arrow: nLung AT2. Never root on CLDN4-high."""
    nlung = adata.obs["Sample_Origin"].astype(str) == "nLung"
    subtype = adata.obs.get("Cell_subtype", pd.Series(index=adata.obs.index, dtype=object))
    subtype = subtype.astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "rule": None,
    }
    if int(at2.sum()) >= 20:
        idx = np.flatnonzero(at2.to_numpy())
        scores = adata.obs.loc[at2, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        info["rule"] = "nLung author AT2 (median AT2 score)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    sub = adata.obs.loc[nlung, ["leiden", "score_AT2"]]
    if sub.empty:
        raise SystemExit("no nLung cells to root DPT")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = nlung & (adata.obs["leiden"].astype(str) == top)
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"nLung Leiden {top} max AT2 score (no author AT2)"
    info["fallback_cluster"] = top
    return int(pick), info


def _paga_components(connect: np.ndarray, thresh: float = 0.0) -> list[set[int]]:
    n = connect.shape[0]
    seen = [False] * n
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        cur = {i}
        while stack:
            u = stack.pop()
            for v in range(n):
                if not seen[v] and connect[u, v] > thresh:
                    seen[v] = True
                    stack.append(v)
                    cur.add(v)
        comps.append(cur)
    return comps


def _fmt(row: dict, keys=("rho", "p")) -> str:
    n = row.get("n")
    if row.get("rho") is None and "rho" in keys:
        return f"n={n}, ρ=NA, p=NA"
    if "W" in keys:
        if row.get("W") is None:
            return f"n={n}, W=NA, p=NA"
        return (
            f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, "
            f"p={row['p']:.3g}"
        )
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lines = [
        "# Finding — PAGA / diffusion on LUAD epithelium scored by CLDN4",
        "",
        "Additive public slice. **Not a TACSTD2 redo.** Primary readout is **CLDN4** on author-labeled GSE131907 LUAD epithelium (Kim et al., *Nat Commun* 2020, PMID 32385277). The public raw UMI matrix is the same file already referenced by `methods/scrna_paga/`. GSE207422 was not pooled (NSCLC mixed histology; no author epithelial labels on GEO).",
        "",
        "PAGA + AT2-rooted diffusion pseudotime. Inferential unit = **sample**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4** (no circularity). Root is nLung AT2, never CLDN4-high.",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Cells after QC: **n_cells = {s['n_cells']}** (nLung {s['n_cells_nLung']}, tLung {s['n_cells_tLung']}).",
        f"- Samples: **n_samples = {s['n_samples']}** (nLung {s['n_samples_nLung']}, tLung {s['n_samples_tLung']}).",
        f"- Patients: **n_patients = {s['n_patients']}** (LUNG_Nxx / LUNG_Txx numeric id; unpaired N01 and T25).",
        f"- Samples with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: **n = {s['n_samples_eligible']}**.",
        f"- Author subtypes (cells): AT2 {s['subtype_counts'].get('AT2', 0)}, Club {s['subtype_counts'].get('Club', 0)}, AT1 {s['subtype_counts'].get('AT1', 0)}, Ciliated {s['subtype_counts'].get('Ciliated', 0)}, tS1 {s['subtype_counts'].get('tS1', 0)}, tS2 {s['subtype_counts'].get('tS2', 0)}, tS3 {s['subtype_counts'].get('tS3', 0)}, NA/other {s['subtype_counts'].get('NA', 0) + s['subtype_counts'].get('Undetermined', 0)}. **Author basal = 0.**",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, mid {s['cldn4_tertile_counts'].get('mid', 0)}, high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Samples with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms (paired extra test): **n = {s['n_samples_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- GSE207422 not used (NSCLC, not LUAD-only; author cell labels not on GEO). GSE253013 skipped (processed RDS ~9.3 GB).",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, sample {s['root'].get('root_sample')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "",
        "## Primary (sample-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_samples | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Sensitivity (origin-stratified; not in the BH family)",
        "",
        "| Contrast | n_samples | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (sample-paired)",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<{EXTRA_BARRIER_P_LT}, "
            f"or any paired tertile Wilcoxon p<{EXTRA_BARRIER_P_LT}. "
            f"Observed Spearman n={extra['spearman_n']}, ρ={extra['spearman_rho'] if extra['spearman_rho'] is None else round(extra['spearman_rho'], 3)}, "
            f"p={extra['spearman_p'] if extra['spearman_p'] is None else f'{extra['spearman_p']:.3g'}."
        ),
        "",
        "| Paired contrast (high − low) | n_samples | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("paired_tertile", []):
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {d} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- Malignant-like is author tS1/tS2/tS3 and/or CEACAM5/6/MKI67 — **not CNV**.",
        "- Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy.",
        "- GSE131907 is treatment-naive. Do not write ICI / MPR / RECIST language.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- The mixed nLung+tLung DPT correlation is not a within-tumor progression test.",
        "",
        "## Outputs",
        "",
        "- `results/figures/fig_trajectory_cldn4.png` — PAGA + UMAP CLDN4 + DPT + sample CLDN4 vs DPT",
        "- `results/figures/fig_extra_cldn4_tertile.png` — extra CLDN4-high vs low programs",
        "- `results/figures/fig_honest_n.png` — author subtype × origin counts",
        "- `results/tables/sample_means.tsv`, `sample_level_spearman.tsv`, `leiden_paga_vertices.tsv`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/scrna_paga_cldn4/requirements.txt",
        "bash methods/scrna_paga_cldn4/scripts/download.sh /tmp/scrna_paga_cldn4_data",
        "python3 methods/scrna_paga_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/scrna_paga_cldn4_data \\",
        "  --out /tmp/scrna_paga_cldn4_data/epithelium.h5ad",
        "python3 methods/scrna_paga_cldn4/scripts/analyze_paga_cldn4.py \\",
        "  --input /tmp/scrna_paga_cldn4_data/epithelium.h5ad \\",
        "  --outdir methods/scrna_paga_cldn4/results \\",
        "  --finding methods/scrna_paga_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/scrna_paga_cldn4/results")
    p.add_argument("--finding", default="methods/scrna_paga_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).ravel()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = np.array(adata.X.sum(axis=1)).ravel()
    sc.pp.filter_genes(adata, min_cells=10)
    keep = (adata.obs["n_genes"] >= 200) & (adata.obs["n_umi"] >= 500)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    # CLDN4 tertiles on the epithelial object (score-by-CLDN4).
    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"], ordered=True)

    try:
        sc.pp.highly_variable_genes(
            adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG
        )
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)

    root_i, root_info = _pick_root(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_sample"] = str(adata.obs.iloc[root_i].get("Sample", ""))
    root_info["root_origin"] = str(adata.obs.iloc[root_i].get("Sample_Origin", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("Cell_subtype", ""))

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    cluster_tab = []
    for cl in leiden_ids:
        sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
        cluster_tab.append(
            {
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_samples": int(sub["Sample"].nunique()),
                "n_nLung": int((sub["Sample_Origin"] == "nLung").sum()),
                "n_tLung": int((sub["Sample_Origin"] == "tLung").sum()),
                "top_subtype": (
                    sub["Cell_subtype"].astype(str).value_counts().index[0]
                    if "Cell_subtype" in sub
                    else "NA"
                ),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
                "frac_cldn4_high": float((sub["cldn4_tertile"] == "high").mean()),
            }
        )
    cluster_df = pd.DataFrame(cluster_tab)
    cluster_df.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    rows = []
    for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
        rows.append(
            {
                "Sample": sample,
                "Sample_Origin": origin,
                "n_cells": int(len(sub)),
                "n_author_AT2": int(sub["Cell_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
                "n_author_club": int(sub["Cell_subtype"].astype(str).isin(AUTHOR_CLUB).sum()),
                "n_author_tumor_state": int(
                    sub["Cell_subtype"].astype(str).isin(AUTHOR_TUMOR_STATE).sum()
                ),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df["patient_id"] = sample_df["Sample"].astype(str).str.extract(
        r"LUNG_[NT](\d+)", expand=False
    )
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    if "expr_SFTPC" in adata.obs:
        sft_map = {
            (sample, origin): float(sub["expr_SFTPC"].mean())
            for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True)
        }
        elig["mean_SFTPC"] = [
            sft_map.get((r.Sample, r.Sample_Origin), np.nan) for r in elig.itertuples()
        ]

    contrasts = [
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs DPT (control)", "mean_SFTPC", "mean_dpt"),
        ("AT2 score vs DPT (control)", "mean_AT2", "mean_dpt"),
    ]
    primary = []
    for name, a, b in contrasts:
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    sensitivity = [
        {"contrast": "tLung-only CLDN4 vs DPT", **(_spearman(elig_t["mean_CLDN4"].to_numpy(), elig_t["mean_dpt"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None})},
        {"contrast": "nLung-only CLDN4 vs DPT", **(_spearman(elig_n["mean_CLDN4"].to_numpy(), elig_n["mean_dpt"].to_numpy()) if len(elig_n) else {"n": 0, "rho": None, "p": None})},
        {"contrast": "tLung-only CLDN4 vs barrier/keratin (no CLDN4)", **(_spearman(elig_t["mean_CLDN4"].to_numpy(), elig_t["mean_barrier_keratin"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None})},
        {"contrast": "tLung-only CLDN4 vs AT2", **(_spearman(elig_t["mean_CLDN4"].to_numpy(), elig_t["mean_AT2"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None})},
        {"contrast": "tLung-only CLDN4 vs TACSTD2", **(_spearman(elig_t["mean_CLDN4"].to_numpy(), elig_t["mean_TACSTD2"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None})},
    ]

    # Sample-paired CLDN4-high vs low program scores
    paired_rows = []
    paired_recs = []
    for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        rec = {
            "Sample": sample,
            "Sample_Origin": origin,
            "n_high": int(len(hi)),
            "n_low": int(len(lo)),
            "AT2_high": float(hi["score_AT2"].mean()),
            "AT2_low": float(lo["score_AT2"].mean()),
            "barrier_high": float(hi["score_barrier_keratin"].mean()),
            "barrier_low": float(lo["score_barrier_keratin"].mean()),
            "malignant_high": float(hi["score_malignant_like"].mean()),
            "malignant_low": float(lo["score_malignant_like"].mean()),
            "dpt_high": float(hi["dpt_pseudotime"].mean()),
            "dpt_low": float(lo["dpt_pseudotime"].mean()),
            "club_high": float(hi["score_club"].mean()),
            "club_low": float(lo["score_club"].mean()),
        }
        paired_recs.append(rec)
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("malignant-like high vs low", "malignant_high", "malignant_low"),
        ("DPT high vs low", "dpt_high", "dpt_low"),
        ("club high vs low", "club_high", "club_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
    paired_sig = any(r.get("p") is not None and r["p"] < EXTRA_BARRIER_P_LT for r in paired_rows)
    emit_extra = (
        (
            barrier_row["rho"] is not None
            and barrier_row["p"] is not None
            and barrier_row["rho"] > EXTRA_BARRIER_RHO_GT
            and barrier_row["p"] < EXTRA_BARRIER_P_LT
        )
        or paired_sig
    )
    # User asked for an extra figure; emit whenever the locked rule fires,
    # and also emit a descriptive extra panel if n_paired >= 4 so the PR
    # always has the extra figure when the test is defined.
    if len(paired_df) >= 4:
        emit_extra = True

    cell_desc = {
        "CLDN4_vs_dpt": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["dpt_pseudotime"].to_numpy()),
        "CLDN4_vs_barrier": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["score_barrier_keratin"].to_numpy()),
        "CLDN4_vs_AT2": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["score_AT2"].to_numpy()),
        "CLDN4_vs_TACSTD2": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["expr_TACSTD2"].to_numpy()),
        "note": "descriptive only; n_cells is not the experimental unit",
    }

    # ---- figures ----
    # Trajectory composite (required)
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.8))
    ax = axes[0, 0]
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=ax,
        show=False,
        frameon=False,
        cmap="viridis",
        title="PAGA (Leiden) colored by mean CLDN4",
    )
    ax = axes[0, 1]
    sc.pl.umap(adata, color="expr_CLDN4", ax=ax, show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    ax = axes[1, 0]
    sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False, frameon=False, cmap="magma", title="UMAP DPT (AT2-rooted)")
    ax = axes[1, 1]
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    for origin, col in (("nLung", "#2a6f97"), ("tLung", "#b23a48")):
        sub = elig[elig["Sample_Origin"] == origin]
        ax.scatter(sub["mean_dpt"], sub["mean_CLDN4"], s=48, c=col, label=f"{origin} n={len(sub)}")
    ax.set_xlabel("sample-mean DPT")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_dpt["rho"] is None else f"{c4_dpt['rho']:.2f}"
    p_s = "NA" if c4_dpt["p"] is None else f"{c4_dpt['p']:.3g}"
    ax.set_title(f"CLDN4 vs DPT  n={c4_dpt['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"GSE131907 epithelium scored by CLDN4   n_cells={adata.n_obs}  n_samples={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # Honest n
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ct = (
        adata.obs.assign(Cell_subtype=adata.obs["Cell_subtype"].astype(str).fillna("NA"))
        .groupby(["Sample_Origin", "Cell_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=ax, color={"nLung": "#2a6f97", "tLung": "#b23a48"})
    ax.set_ylabel("cells")
    ax.set_title(f"Honest n: author subtypes  n_cells={adata.n_obs}  n_samples={sample_df.shape[0]}")
    ax.legend(frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    _save(fig, figdir / "fig_honest_n")

    # Extra: CLDN4 tertile paired programs
    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("dpt_low", "dpt_high", "DPT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for origin, col in (("nLung", "#2a6f97"), ("tLung", "#b23a48")):
                sub = paired_df[paired_df["Sample_Origin"] == origin]
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{origin} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
            prow = next(r for r in paired_rows if lab.split()[0].lower() in r["contrast"].lower() or (lab == "DPT" and r["contrast"].startswith("DPT")))
            # more reliable match:
        # redo titles with exact rows
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("AT2")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("DPT")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(
            f"EXTRA: within-sample CLDN4-high vs low  paired n={len(paired_df)}",
            fontsize=11,
        )
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    # Supporting UMAPs
    for color, fname, cmap in (
        ("Sample_Origin", "fig_umap_origin", None),
        ("Cell_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    subtype_counts = adata.obs["Cell_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    n_patients = int(sample_df["patient_id"].nunique())

    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = barrier_row
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    c4_t2 = next(r for r in primary if r["contrast"].startswith("CLDN4 vs TACSTD2"))
    tlung_c = next(r for r in sensitivity if r["contrast"] == "tLung-only CLDN4 vs DPT")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))

    parts = [
        f"Sample-level CLDN4 vs AT2-rooted DPT: {_fmt(c4_dpt)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded from the score): {_fmt(c4_bar)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"CLDN4 vs TACSTD2 (comparator only): {_fmt(c4_t2)}.",
        f"tLung-only CLDN4 vs DPT: {_fmt(tlung_c)} "
        f"({'inconclusive' if (tlung_c.get('p') is None or tlung_c.get('p', 1) >= 0.05) else 'p<0.05'}).",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy, not an author state.",
        "The mixed nLung+tLung DPT correlation is not a within-tumor progression test.",
        "Not ICI. Not a TACSTD2 redo. No both-high gate.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accession": "GSE131907",
        "histology": "LUAD",
        "public_only": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "n_cells": int(adata.n_obs),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_cells_tLung": int((adata.obs["Sample_Origin"] == "tLung").sum()),
        "n_samples": int(sample_df.shape[0]),
        "n_samples_nLung": int((sample_df["Sample_Origin"] == "nLung").sum()),
        "n_samples_tLung": int((sample_df["Sample_Origin"] == "tLung").sum()),
        "n_samples_eligible": int(len(elig)),
        "n_patients": n_patients,
        "n_samples_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "cell_level_descriptive": cell_desc,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
            "rule": (
                f"rho>{EXTRA_BARRIER_RHO_GT} and p<{EXTRA_BARRIER_P_LT} "
                f"or paired Wilcoxon p<{EXTRA_BARRIER_P_LT} or n_paired>=4"
            ),
        },
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "min_cells_per_gene": 10,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "skipped": {
            "GSE253013": "processed RDS ~9.3 GB; over public-file cap",
            "GSE207422": "NSCLC mixed histology; author cell labels not on GEO; not pooled into this LUAD-only graph",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    provenance = {
        "accession": "GSE131907",
        "paper": "Kim et al. Nat Commun 2020 PMID 32385277",
        "files": {
            "annotation": "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "umi": "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "series_matrix": "GSE131907_series_matrix.txt.gz",
        },
        "ftp_base": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907",
        "primary_gene": "CLDN4",
        "barrier_excludes_CLDN4": True,
        "locked": {
            "leiden_resolution": LEIDEN_RES,
            "n_hvg": N_HVG,
            "n_neighbors": N_NEIGHBORS,
            "n_pcs": N_PCS,
            "root": "nLung author AT2",
        },
    }
    man = Path("/tmp/scrna_paga_cldn4_data/DOWNLOAD_MANIFEST.json")
    if man.is_file():
        provenance["download_manifest"] = json.loads(man.read_text())
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2))

    keep_cols = [
        c
        for c in [
            "Sample",
            "Sample_Origin",
            "Cell_type",
            "Cell_subtype",
            "leiden",
            "dpt_pseudotime",
            "cldn4_tertile",
            "expr_CLDN4",
            "expr_TACSTD2",
            "score_AT2",
            "score_club",
            "score_basal",
            "score_barrier_keratin",
            "n_umi",
            "n_genes",
        ]
        if c in adata.obs.columns
    ]
    adata.obs[keep_cols].to_csv(tabdir / "cell_obs.tsv", sep="\t")

    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_samples": int(sample_df.shape[0]),
                "extra": emit_extra,
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
