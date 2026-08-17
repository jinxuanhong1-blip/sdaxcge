#!/usr/bin/env python3
"""PAGA + diffusion on GSE148071 epithelium, scored by CLDN4.

ADDITIVE. Not a TACSTD2 redo: CLDN4 is the primary readout; TACSTD2 is a
comparator only. Barrier/keratin score excludes CLDN4 (no circularity).
Patient/sample (not cell) is the inferential unit. GEO n=42 is not the test n.
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
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STATES,
    TISCH_LEFTOVER_EPI,
    TISCH_MALIGNANT,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_GENES = 200
MIN_UMI = 500
MIN_CELLS_PER_GENE = 10
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
MAX_CELLS_PER_PATIENT = 500
GRAPH_CAP_TRIGGER = 20000
RANDOM_SEED = 0
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
        return {
            "n": n,
            "W": None,
            "p": None,
            "delta_median": float(np.median(aa - bb)),
            "note": "wilcoxon failed",
        }
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
    """Alveolar / AT2-high root. Never root on CLDN4-high."""
    lin = adata.obs["lineage"].astype(str)
    alveolar = lin.eq("Alveolar")
    info = {
        "n_alveolar": int(alveolar.sum()),
        "rule": None,
    }
    if int(alveolar.sum()) >= 20:
        idx = np.flatnonzero(alveolar.to_numpy())
        scores = adata.obs.loc[alveolar, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        info["rule"] = "TISCH Alveolar (median AT2 score)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    means = (
        adata.obs.groupby("leiden", observed=True)["score_AT2"]
        .mean()
        .sort_values(ascending=False)
    )
    top = str(means.index[0])
    cand = adata.obs["leiden"].astype(str) == top
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"Leiden {top} max AT2 score (alveolar n<20)"
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


def _cap_per_patient(adata, max_n: int, seed: int):
    rng = np.random.default_rng(seed)
    keep = []
    for pat, idx in adata.obs.groupby("patient", observed=True).groups.items():
        idx = np.asarray(idx)
        if len(idx) <= max_n:
            keep.extend(idx.tolist())
        else:
            take = rng.choice(idx, size=max_n, replace=False)
            keep.extend(take.tolist())
    return adata[keep].copy()


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lines = [
        "# Finding — GSE148071 PAGA / diffusion scored by CLDN4",
        "",
        "Additive public slice. **Not a TACSTD2 redo.** Primary readout is **CLDN4** on TISCH-labeled epithelium from GSE148071 (Wu et al., *Nat Commun* 2021, PMID 33953163): 42 advanced NSCLC diagnostic biopsies (Singleron GEXSCOPE). Counts are the public GEO raw UMI matrices. Labels are TISCH2 major-lineage (Malignant + leftover Alveolar / Basal / Epithelial). GEO does not ship author cell types or histology on the series matrix.",
        "",
        "PAGA + Alveolar/AT2-rooted diffusion pseudotime. Inferential unit = **patient** (one sample each). Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4** (no circularity). Root is TISCH Alveolar, never CLDN4-high. **Do not write n=42 as the test n.**",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- GEO patients / samples: **42** (Wu 2021). This is the catalog n, not the test n.",
        f"- TISCH cells: **{s['n_tisch_cells']}**. TISCH epithelial (Malignant+Alveolar+Basal+Epithelial): **{s['n_tisch_epithelial']}**.",
        f"- GEO barcodes matched to TISCH epithelium: **{s['n_matched']}**.",
        f"- After QC (min {MIN_GENES} genes, min {MIN_UMI} UMI): **n_cells_qc = {s['n_cells_qc']}** in **{s['n_patients_qc']}** patients.",
        f"- Graph cap: max {MAX_CELLS_PER_PATIENT}/patient when n>{GRAPH_CAP_TRIGGER} (seed {RANDOM_SEED}). Analysis object: **n_cells = {s['n_cells']}** in **n_patients = {s['n_patients']}**.",
        f"- Lineage on the analysis object: {s['lineage_counts']}.",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} analysis cells (Spearman): **n = {s['n_samples_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} TISCH Malignant cells on the object: **n = {s['n_malignant_eligible']}**. Leftover-epithelium–dominant (malignant <10): **{s['leftover_dominant_patients']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms (paired extra): **n = {s['n_samples_paired_tertile']}**.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, mid {s['cldn4_tertile_counts'].get('mid', 0)}, high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- No ICI / RECIST / MPR labels. Histology is not on the GEO series matrix — no LUAD/LUSC split.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, patient {s['root'].get('root_patient')}).",
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
        "## Sensitivity (not in the BH family)",
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
            f"Rule: always emit (requested extra figure). "
            f"Sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) "
            f"n={extra['spearman_n']}, ρ={extra['spearman_rho'] if extra['spearman_rho'] is None else round(extra['spearman_rho'], 3)}, "
            f"p={extra['spearman_p'] if extra['spearman_p'] is None else f'{extra['spearman_p']:.3g}'}. "
            f"Paired tertile n={s['n_samples_paired_tertile']}."
        ),
        "",
        "| Contrast | n | W | Δmed (high−low) | p |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s.get("paired_tertile", []):
        w = "NA" if r.get("W") is None else f"{r['W']:.1f}"
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {w} | {d} | {pv} |")
    lines += [
        "",
        "## Caveats",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- TISCH Malignant is a processed label, **not CNV re-called here**.",
        "- Several patients are leftover-epithelium–dominant (almost no TISCH Malignant). They stay in the mixed graph and are dropped from the malignant-only sensitivity.",
        "- Patient batch is strong (42 tumors). Per-patient cap reduces one-sample domination; it is not Harmony.",
        "- Mixed advanced NSCLC. Do not write LUAD-only. Do not write ICI language.",
        "- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.",
        "- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/gse148071_paga_cldn4/requirements.txt",
        "bash methods/gse148071_paga_cldn4/scripts/download.sh /tmp/gse148071_paga_data",
        "python3 methods/gse148071_paga_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse148071_paga_data \\",
        "  --out /tmp/gse148071_paga_data/epithelium.h5ad",
        "python3 methods/gse148071_paga_cldn4/scripts/analyze_paga_cldn4.py \\",
        "  --input /tmp/gse148071_paga_data/epithelium.h5ad \\",
        "  --outdir methods/gse148071_paga_cldn4 \\",
        "  --finding methods/gse148071_paga_cldn4/FINDING.md",
        "```",
        "",
        f"Trajectory figure: `figures/fig_trajectory_cldn4.png`. Extra figure: `figures/fig_extra_cldn4_tertile.png`.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default="methods/gse148071_paga_cldn4")
    ap.add_argument("--finding", default="methods/gse148071_paga_cldn4/FINDING.md")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    np.random.seed(RANDOM_SEED)

    adata = sc.read_h5ad(args.input)
    extract = dict(adata.uns.get("extract", {}))
    n_tisch_cells = int(extract.get("n_tisch_cells", 82267))
    n_tisch_epithelial = int(extract.get("n_tisch_epithelial", adata.n_obs))
    n_matched = int(extract.get("n_matched", adata.n_obs))

    adata.obs["lineage"] = adata.obs["lineage"].astype(str).str.strip()
    adata.obs["patient"] = adata.obs["patient"].astype(str)
    adata.obs["n_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["n_genes"] >= MIN_GENES].copy()
    adata = adata[adata.obs["n_umi"] >= MIN_UMI].copy()
    n_cells_qc = int(adata.n_obs)
    n_patients_qc = int(adata.obs["patient"].nunique())
    print(f"after QC n_cells={n_cells_qc} n_patients={n_patients_qc}", flush=True)

    capped = False
    if adata.n_obs > GRAPH_CAP_TRIGGER:
        adata = _cap_per_patient(adata, MAX_CELLS_PER_PATIENT, RANDOM_SEED)
        capped = True
        print(
            f"capped to max {MAX_CELLS_PER_PATIENT}/patient → n_cells={adata.n_obs}",
            flush=True,
        )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            X = adata[:, g].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(X, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    q1, q2 = np.nanquantile(adata.obs["expr_CLDN4"].to_numpy(), [1 / 3, 2 / 3])
    tert = pd.Series("mid", index=adata.obs_names)
    tert[adata.obs["expr_CLDN4"] <= q1] = "low"
    tert[adata.obs["expr_CLDN4"] > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"])

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)
    sc.tl.paga(adata, groups="leiden")
    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = sorted(adata.obs["leiden"].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)

    root_idx, root_info = _pick_root(adata)
    root_info["index"] = root_idx
    root_info["root_patient"] = str(adata.obs.iloc[root_idx]["patient"])
    root_info["root_lineage"] = str(adata.obs.iloc[root_idx]["lineage"])
    root_info["root_cldn4"] = float(adata.obs.iloc[root_idx]["expr_CLDN4"])
    adata.uns["iroot"] = root_idx
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    # sample table
    rows = []
    for pat, sub in adata.obs.groupby("patient", observed=True):
        n_mal = int(sub["lineage"].isin(TISCH_MALIGNANT).sum())
        n_leftover = int(sub["lineage"].isin(TISCH_LEFTOVER_EPI).sum())
        rows.append(
            {
                "patient": pat,
                "n_cells": int(len(sub)),
                "n_malignant": n_mal,
                "n_leftover_epi": n_leftover,
                "leftover_dominant": n_mal < MIN_CELLS_PER_SAMPLE_FOR_MEAN,
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()) if "expr_TACSTD2" in sub else np.nan,
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()) if "expr_SFTPC" in sub else np.nan,
            }
        )
    sample_df = pd.DataFrame(rows).sort_values("patient")
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    malig = sample_df[sample_df["n_malignant"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    leftover_dom = sample_df.loc[sample_df["leftover_dominant"], "patient"].tolist()

    primary_specs = [
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin score", "mean_CLDN4", "mean_barrier"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("TACSTD2 vs DPT (comparator)", "mean_TACSTD2", "mean_dpt"),
        ("SFTPC vs DPT (control)", "mean_SFTPC", "mean_dpt"),
    ]
    primary = []
    for name, a, b in primary_specs:
        r = _spearman(elig[a].to_numpy(), elig[b].to_numpy())
        r["contrast"] = name
        primary.append(r)
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q

    sensitivity = []
    for name, df, a, b in (
        ("malignant-only CLDN4 vs DPT", malig, "mean_CLDN4", "mean_dpt"),
        ("malignant-only CLDN4 vs AT2 score", malig, "mean_CLDN4", "mean_AT2"),
        ("malignant-only CLDN4 vs barrier/keratin", malig, "mean_CLDN4", "mean_barrier"),
        ("drop leftover-dominant CLDN4 vs DPT", sample_df[~sample_df["leftover_dominant"]], "mean_CLDN4", "mean_dpt"),
    ):
        use = df[df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN]
        r = _spearman(use[a].to_numpy(), use[b].to_numpy())
        r["contrast"] = name
        sensitivity.append(r)

    # paired tertile extra
    paired_rows_data = []
    for pat, sub in adata.obs.groupby("patient", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_rows_data.append(
            {
                "patient": pat,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_rows_data)
    paired_rows = []
    if not paired_df.empty:
        for contrast, hi, lo in (
            ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
            ("AT2 score high vs low", "AT2_high", "AT2_low"),
            ("DPT high vs low", "dpt_high", "dpt_low"),
        ):
            r = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            r["contrast"] = contrast
            paired_rows.append(r)
    else:
        paired_rows = [
            {"contrast": "barrier/keratin (no CLDN4) high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
            {"contrast": "AT2 score high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
            {"contrast": "DPT high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
        ]

    barrier_row = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
    emit_extra = True  # requested extra figure

    # ----- figures -----
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    sc.pl.paga(
        adata,
        color="leiden",
        ax=axes[0],
        show=False,
        frameon=False,
        title="PAGA (Leiden)",
    )
    sc.pl.umap(
        adata,
        color="expr_CLDN4",
        ax=axes[1],
        show=False,
        frameon=False,
        cmap="viridis",
        title="UMAP CLDN4",
    )
    ax = axes[2]
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    leftover = elig["leftover_dominant"]
    ax.scatter(
        elig.loc[~leftover, "mean_dpt"],
        elig.loc[~leftover, "mean_CLDN4"],
        s=48,
        c="#b23a48",
        label=f"malignant-bearing n={int((~leftover).sum())}",
    )
    ax.scatter(
        elig.loc[leftover, "mean_dpt"],
        elig.loc[leftover, "mean_CLDN4"],
        s=48,
        c="#2a6f97",
        label=f"leftover-epi dominant n={int(leftover.sum())}",
    )
    ax.set_xlabel("sample-mean DPT")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_dpt["rho"] is None else f"{c4_dpt['rho']:.2f}"
    p_s = "NA" if c4_dpt["p"] is None else f"{c4_dpt['p']:.3g}"
    ax.set_title(f"CLDN4 vs DPT  n={c4_dpt['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(
        f"GSE148071 epithelium scored by CLDN4   n_cells={adata.n_obs}  n_patients={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    ct = (
        adata.obs.groupby(["patient", "lineage"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.85)
    ax.set_ylabel("cells in analysis object")
    ax.set_title(
        f"Honest n: TISCH lineages after QC/cap  n_cells={adata.n_obs}  "
        f"n_patients={sample_df.shape[0]} / GEO 42"
    )
    ax.legend(frameon=False, fontsize=7)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", fontsize=7)
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("dpt_low", "dpt_high", "DPT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            ax.scatter(paired_df[lo], paired_df[hi], s=44, c="#4c6a92")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot(
                [lims[0] - pad, lims[1] + pad],
                [lims[0] - pad, lims[1] + pad],
                ls="--",
                c="0.6",
                lw=1,
            )
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("AT2")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("DPT")), keys=("W", "p")))
        fig.suptitle(
            f"EXTRA: within-sample CLDN4-high vs low  paired n={len(paired_df)}",
            fontsize=11,
        )
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("lineage", "fig_umap_lineage", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("dpt_pseudotime", "fig_umap_dpt", "viridis"),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = barrier_row
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    c4_t2 = next(r for r in primary if r["contrast"].startswith("CLDN4 vs TACSTD2"))
    mal_c = next(r for r in sensitivity if r["contrast"] == "malignant-only CLDN4 vs DPT")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))

    parts = [
        f"Sample-level CLDN4 vs Alveolar-rooted DPT: {_fmt(c4_dpt)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded from the score): {_fmt(c4_bar)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"CLDN4 vs TACSTD2 (comparator only): {_fmt(c4_t2)}.",
        f"Malignant-only CLDN4 vs DPT: {_fmt(mal_c)} "
        f"({'inconclusive' if (mal_c.get('p') is None or mal_c.get('p', 1) >= 0.05) else 'p<0.05'}).",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        f"Leftover-epithelium–dominant patients (TISCH Malignant <10 on the object): {leftover_dom or 'none'}.",
        "GEO n=42 is the catalog, not the Spearman n. Mixed advanced NSCLC; no ICI labels. Not a TACSTD2 redo. No both-high gate.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accession": "GSE148071",
        "histology": "advanced NSCLC (mixed; not on GEO series matrix)",
        "public_only": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "n_tisch_cells": n_tisch_cells,
        "n_tisch_epithelial": n_tisch_epithelial,
        "n_matched": n_matched,
        "n_cells_qc": n_cells_qc,
        "n_patients_qc": n_patients_qc,
        "n_cells": int(adata.n_obs),
        "n_patients": int(sample_df.shape[0]),
        "n_samples_eligible": int(len(elig)),
        "n_malignant_eligible": int(len(malig)),
        "n_samples_paired_tertile": int(len(paired_df)),
        "leftover_dominant_patients": leftover_dom,
        "lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "graph_capped": capped,
        "max_cells_per_patient": MAX_CELLS_PER_PATIENT if capped else None,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
            "rule": "always emit extra figure",
        },
        "verdict": verdict,
    }
    (outdir / "tables" / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig.to_csv(tabdir / "sample_eligible.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "paired_tertile.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tabdir / "primary_spearman.tsv", sep="\t", index=False)
    write_finding(Path(args.finding), {"summary": summary})

    keep_cols = [
        c
        for c in [
            "patient",
            "gsm",
            "lineage",
            "malignancy",
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
                "n_cells": int(adata.n_obs),
                "n_patients": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "extra": emit_extra,
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
