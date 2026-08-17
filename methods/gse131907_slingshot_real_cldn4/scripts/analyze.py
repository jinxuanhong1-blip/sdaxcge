#!/usr/bin/env python3
"""REAL Slingshot on GSE131907 epithelium/malignant, scored by CLDN4.

ADDITIVE. Not a redo of PR #325 (DPT/PAGA only). Not the winning-pair
Slingshot/DPT fallback (PR #449). Primary clock is Street 2018 Slingshot
(R if present, otherwise the Python equivalent in slingshot_py.py).
Root = nLung author AT2, never CLDN4-high. Inferential unit = sample.
Barrier/keratin and IFN scores exclude CLDN4. No TACSTD2∩CLDN4 dual-high.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
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
from slingshot_py import run_slingshot  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
SLING_DIMS = 5
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
N_LINEAGE_BINS = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05


def slingshot_r_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {"available": False, "reason": "Rscript not on PATH"}
    try:
        proc = subprocess.run(
            [rscript, "-e", 'cat(as.character(packageVersion("slingshot")))'],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"available": False, "reason": f"Rscript probe failed: {exc}"}
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:300],
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


def try_r_slingshot(
    embedding: np.ndarray,
    labels: np.ndarray,
    cell_ids: list[str],
    start_cluster: str,
    workdir: Path,
) -> dict | None:
    status = slingshot_r_status()
    if not status.get("available"):
        return None
    rscript = shutil.which("Rscript")
    emb_path = workdir / "slingshot_embedding.tsv"
    out_path = workdir / "slingshot_r_pseudotime.tsv"
    dim_names = [f"dim{i+1}" for i in range(embedding.shape[1])]
    frame = pd.DataFrame(embedding, columns=dim_names)
    frame.insert(0, "cluster", np.asarray(labels).astype(str))
    frame.insert(0, "cell_id", cell_ids)
    frame.to_csv(emb_path, sep="\t", index=False)
    r_script = Path(__file__).resolve().parent / "run_slingshot.R"
    proc = subprocess.run(
        [rscript, str(r_script), str(emb_path), str(start_cluster), str(out_path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if proc.returncode != 0 or not out_path.is_file():
        return {
            "ok": False,
            "stderr": (proc.stderr or proc.stdout or "")[:500],
            "status": status,
        }
    tab = pd.read_csv(out_path, sep="\t")
    pt_cols = [c for c in tab.columns if c not in {"cell_id", "cluster"} and not str(c).startswith("weight_")]
    return {
        "ok": True,
        "status": status,
        "table": tab,
        "pt_cols": pt_cols,
        "stdout": (proc.stdout or "").strip()[:300],
    }


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
    """External arrow: nLung AT2. Never root on CLDN4-high."""
    nlung = adata.obs["Sample_Origin"].astype(str) == "nLung"
    subtype = adata.obs.get("author_subtype", adata.obs.get("Cell_subtype"))
    subtype = subtype.astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    high = adata.obs["cldn4_tertile"].astype(str) == "high"
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "n_nLung_AT2_not_cldn4_high": int((at2 & ~high).sum()),
        "rule": None,
        "rejected_cldn4_high_root": False,
    }
    pool = at2 & ~high
    if int(pool.sum()) < 20:
        pool = at2
        info["note"] = "CLDN4-high filter left <20 AT2; used all nLung AT2"
    if int(pool.sum()) >= 20:
        idx = np.flatnonzero(pool.to_numpy())
        scores = adata.obs.loc[pool, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        if bool(high.iloc[pick]):
            # should not happen if pool excluded high
            info["rejected_cldn4_high_root"] = True
            not_high_idx = np.flatnonzero((at2 & ~high).to_numpy())
            if not_high_idx.size:
                scores2 = adata.obs.iloc[not_high_idx]["score_AT2"].to_numpy()
                pick = not_high_idx[int(np.nanargmin(np.abs(scores2 - np.nanmedian(scores2))))]
        info["rule"] = "nLung author AT2 (median AT2 score); never CLDN4-high"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    sub = adata.obs.loc[nlung & ~high, ["leiden", "score_AT2"]]
    if sub.empty:
        sub = adata.obs.loc[nlung, ["leiden", "score_AT2"]]
    if sub.empty:
        raise SystemExit("no nLung cells to root Slingshot")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = nlung & (adata.obs["leiden"].astype(str) == top) & ~high
    if int(cand.sum()) == 0:
        cand = nlung & (adata.obs["leiden"].astype(str) == top)
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"nLung Leiden {top} max AT2 score (no author AT2); never CLDN4-high"
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
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    clock = s["clock"]["name"]
    lines = [
        "# Finding — REAL Slingshot on GSE131907 epithelium/malignant, CLDN4-only",
        "",
        "ADDITIVE. **CLDN4 only.** Kim et al., *Nat Commun* 2020, PMID 32385277 "
        "(GSE131907). Epithelium / malignant cells only. This folder does **not** "
        "redo PR #325 (DPT/PAGA on nLung+tLung). It does **not** redo PR #449 "
        "(winning-pair DPT fallback; Slingshot R missing). GSE207422 and GSE205335 "
        "are not added. No TACSTD2∩CLDN4 dual-high gate.",
        "",
        f"Primary clock: **{clock}**. Inferential unit = **sample**. "
        "Cell-level ρ is descriptive. Barrier/keratin and IFN scores **exclude CLDN4**. "
        "Root is nLung author AT2, never CLDN4-high.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Catalog epithelium/malignant (before cap): **n_cells_catalog = {s['n_cells_catalog']}** "
        f"across **n_samples_catalog = {s['n_samples_catalog']}** "
        f"(origins {s['catalog_by_origin']}).",
        f"- Analysis cells after QC (cap ≤{s['cap_per_sample']}/sample, nLung AT2 protected): "
        f"**n_cells = {s['n_cells']}** "
        f"(nLung {s['n_cells_nLung']}, tLung {s['n_cells_tLung']}, "
        f"other malignant sites {s['n_cells_other']}).",
        f"- Samples: **n_samples = {s['n_samples']}** "
        f"(nLung {s['n_samples_nLung']}, tLung {s['n_samples_tLung']}, "
        f"other {s['n_samples_other']}).",
        f"- Patients: **n_patients = {s['n_patients']}**.",
        f"- Samples with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} cells used for Spearman: "
        f"**n = {s['n_samples_eligible']}**.",
        f"- Author AT2 in the object: **{s['n_author_AT2']}**. Author basal cells: "
        f"**{s['subtype_counts'].get('Basal', 0)}**.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Samples with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and "
        f"CLDN4-low arms: **n = {s['n_samples_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot lineages: **{s['n_lineages']}** from start cluster "
        f"{s['root'].get('start_cluster')} (root sample {s['root'].get('root_sample')}).",
        f"- Slingshot engine: {s['clock']['engine']}.",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among "
        f"{s['n_leiden']} Leiden vertices.",
        "- GSE207422 / GSE205335 not used. PE unlabeled epithelium dropped.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Slingshot embedding: first {SLING_DIMS} PCs (UMAP is visualization only).",
        f"- Root: {s['root']['rule']} (root cell index {s['root']['index']}, "
        f"sample {s['root'].get('root_sample')}, subtype {s['root'].get('root_subtype')}, "
        f"CLDN4 tertile {s['root'].get('root_cldn4_tertile')}).",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN genes: compact ISG panel in `scripts/gene_sets.py` (**CLDN4 out**).",
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
            f"Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<{EXTRA_BARRIER_P_LT}, "
            f"or any paired tertile Wilcoxon p<{EXTRA_BARRIER_P_LT}, or n_paired≥4. "
            f"Observed Spearman n={extra['spearman_n']}, ρ={extra_rho}, p={extra_p}."
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
        "- This is not a redo of PR #325 (GSE131907 DPT/PAGA).",
        "- The pooled Slingshot Spearman mixes nLung, tLung, and metastatic sites "
        "and is **not** a within-tumor progression test.",
        "- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or "
        "CEACAM5/6/MKI67 — **not CNV**.",
        "- GSE131907 is treatment-naive. Do not write ICI / MPR / RECIST language.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because Slingshot/PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot is an ordering on an embedding, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/sample_level_spearman.tsv` — **done criterion**",
        "- `results/tables/sample_means.tsv`",
        "- `results/tables/slingshot_lineages.tsv`",
        "- `results/tables/lineage_bin_means.tsv` — CLDN4 / barrier / IFN along lineages",
        "- `results/tables/leiden_paga_vertices.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_extra_lineage_programs.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/gse131907_slingshot_real_cldn4/requirements.txt",
        "python3 methods/gse131907_slingshot_real_cldn4/scripts/download.py \\",
        "  --out /tmp/gse131907_slingshot_data",
        "python3 methods/gse131907_slingshot_real_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse131907_slingshot_data \\",
        "  --out /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad",
        "python3 methods/gse131907_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad \\",
        "  --outdir methods/gse131907_slingshot_real_cldn4/results \\",
        "  --finding methods/gse131907_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/gse131907_slingshot_real_cldn4/results")
    p.add_argument("--finding", default="methods/gse131907_slingshot_real_cldn4/FINDING.md")
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
    extract_inv = dict(adata.uns.get("extract_inventory", {}))
    inv_json = Path(args.input).with_suffix(".inventory.json")
    if inv_json.is_file():
        try:
            disk_inv = json.loads(inv_json.read_text())
            extract_inv = {**extract_inv, **disk_inv}
        except json.JSONDecodeError:
            pass
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    if "author_subtype" not in adata.obs.columns:
        adata.obs["author_subtype"] = adata.obs.get(
            "Cell_subtype", pd.Series("NA", index=adata.obs.index)
        ).astype(str)

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

    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(
        tert, categories=["low", "mid", "high"], ordered=True
    )

    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
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
    root_obs = adata.obs.iloc[root_i]
    if str(root_obs.get("cldn4_tertile", "")) == "high":
        raise SystemExit("root landed on CLDN4-high; refusing")
    adata.uns["iroot"] = root_i
    # DPT is a sensitivity clock only — not the primary Slingshot claim.
    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].replace(
        [np.inf, -np.inf], np.nan
    )
    start_cluster = str(root_obs["leiden"])
    root_info.update(
        {
            "index": int(root_i),
            "root_sample": str(root_obs.get("Sample", "")),
            "root_origin": str(root_obs.get("Sample_Origin", "")),
            "root_subtype": str(root_obs.get("author_subtype", "")),
            "root_cldn4_tertile": str(root_obs.get("cldn4_tertile", "")),
            "start_cluster": start_cluster,
            "root_CLDN4": float(root_obs.get("expr_CLDN4", np.nan)),
            "root_AT2": float(root_obs.get("score_AT2", np.nan)),
        }
    )

    embedding = np.asarray(adata.obsm["X_pca"][:, :SLING_DIMS], dtype=float)
    labels = adata.obs["leiden"].astype(str).to_numpy()
    py_sling = run_slingshot(embedding, labels, start_cluster)
    adata.obs["sling_pseudotime"] = py_sling["shared_pseudotime"]
    engine = {
        "name": "Slingshot (Street 2018 Python equivalent: MST + principal curves)",
        "engine": "python:slingshot_py",
        "r": slingshot_r_status(),
    }
    r_try = try_r_slingshot(
        embedding,
        labels,
        adata.obs_names.astype(str).tolist(),
        start_cluster,
        outdir,
    )
    if r_try and r_try.get("ok"):
        tab = r_try["table"].set_index("cell_id").reindex(adata.obs_names.astype(str))
        pt = tab[r_try["pt_cols"]].to_numpy(dtype=float)
        adata.obs["sling_pseudotime"] = np.nanmean(pt, axis=1)
        for j, col in enumerate(r_try["pt_cols"]):
            adata.obs[f"sling_{col}"] = pt[:, j]
        engine = {
            "name": "Slingshot (Street 2018; Bioconductor R)",
            "engine": "R:slingshot",
            "r": r_try["status"],
            "pt_cols": r_try["pt_cols"],
        }
        py_sling["r_used"] = True
        py_sling["r_pt_cols"] = r_try["pt_cols"]
    else:
        py_sling["r_used"] = False
        if r_try:
            engine["r_error"] = r_try.get("stderr")
        for j in range(py_sling["n_lineages"]):
            adata.obs[f"sling_Lineage{j+1}"] = py_sling["pseudotime"][:, j]

    # lineage tables
    lin_rows = []
    for j, path in enumerate(py_sling.get("lineages", [])):
        col = f"sling_Lineage{j+1}"
        assigned = adata.obs[col].notna() if col in adata.obs else pd.Series(False, index=adata.obs.index)
        lin_rows.append(
            {
                "lineage": f"Lineage{j+1}",
                "clusters": ">".join(path),
                "n_clusters": int(len(path)),
                "n_cells_assigned": int(assigned.sum()),
                "n_samples_assigned": int(adata.obs.loc[assigned, "Sample"].nunique())
                if assigned.any()
                else 0,
                "mean_CLDN4": float(adata.obs.loc[assigned, "expr_CLDN4"].mean())
                if assigned.any()
                else None,
                "mean_barrier": float(adata.obs.loc[assigned, "score_barrier_keratin"].mean())
                if assigned.any()
                else None,
                "mean_IFN": float(adata.obs.loc[assigned, "score_IFN"].mean())
                if assigned.any()
                else None,
                "start_cluster": path[0] if path else None,
                "end_cluster": path[-1] if path else None,
            }
        )
    lin_df = pd.DataFrame(lin_rows)
    lin_df.to_csv(tabdir / "slingshot_lineages.tsv", sep="\t", index=False)
    pd.DataFrame(py_sling.get("edges", [])).to_csv(
        tabdir / "slingshot_mst_edges.tsv", sep="\t", index=False
    )

    # lineage-binned CLDN4 / barrier / IFN
    bin_rows = []
    for j, path in enumerate(py_sling.get("lineages", [])):
        col = f"sling_Lineage{j+1}"
        if col not in adata.obs:
            continue
        t = adata.obs[col].to_numpy()
        ok = np.isfinite(t)
        if int(ok.sum()) < 20:
            continue
        bins = np.linspace(0, 1, N_LINEAGE_BINS + 1)
        bid = np.digitize(t[ok], bins[1:-1], right=False)
        sub = adata.obs.loc[ok]
        for b in range(N_LINEAGE_BINS):
            m = bid == b
            if int(m.sum()) < 5:
                continue
            piece = sub.iloc[np.flatnonzero(m)]
            bin_rows.append(
                {
                    "lineage": f"Lineage{j+1}",
                    "bin": int(b),
                    "bin_lo": float(bins[b]),
                    "bin_hi": float(bins[b + 1]),
                    "n_cells": int(m.sum()),
                    "n_samples": int(piece["Sample"].nunique()),
                    "mean_pseudotime": float(piece[col].mean()),
                    "mean_CLDN4": float(piece["expr_CLDN4"].mean()),
                    "mean_barrier_keratin": float(piece["score_barrier_keratin"].mean()),
                    "mean_IFN": float(piece["score_IFN"].mean()),
                    "mean_AT2": float(piece["score_AT2"].mean()),
                    "mean_TACSTD2": float(piece["expr_TACSTD2"].mean()),
                }
            )
    bin_df = pd.DataFrame(bin_rows)
    bin_df.to_csv(tabdir / "lineage_bin_means.tsv", sep="\t", index=False)

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
                "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_sling": float(sub["sling_pseudotime"].mean()),
                "frac_cldn4_high": float((sub["cldn4_tertile"] == "high").mean()),
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    rows = []
    for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
        rows.append(
            {
                "Sample": sample,
                "Sample_Origin": origin,
                "n_cells": int(len(sub)),
                "n_author_AT2": int(sub["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
                "n_author_club": int(sub["author_subtype"].astype(str).isin(AUTHOR_CLUB).sum()),
                "n_author_tumor_state": int(
                    sub["author_subtype"].astype(str).isin(AUTHOR_TUMOR_STATE).sum()
                ),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_sling": float(sub["sling_pseudotime"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df["patient_id"] = (
        sample_df["Sample"].astype(str).str.extract(r"(?:LUNG_[NT]|NS_|EBUS_|BRONCHO_|EFFUSION_)?(\d+)", expand=False)
    )
    # Kim LUNG_Nxx / LUNG_Txx share a numeric id; mets keep their own sample id
    lung = sample_df["Sample"].astype(str).str.extract(r"LUNG_[NT](\d+)", expand=False)
    sample_df.loc[lung.notna(), "patient_id"] = lung[lung.notna()]
    if "expr_SFTPC" in adata.obs:
        sft_map = {
            (sample, origin): float(sub["expr_SFTPC"].mean())
            for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True)
        }
        sample_df["mean_SFTPC"] = [
            sft_map.get((r.Sample, r.Sample_Origin), np.nan) for r in sample_df.itertuples()
        ]
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Slingshot", "mean_CLDN4", "mean_sling"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN (no CLDN4)", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Slingshot (control)", "mean_SFTPC", "mean_sling"),
        ("AT2 score vs Slingshot (control)", "mean_AT2", "mean_sling"),
        ("IFN vs Slingshot", "mean_IFN", "mean_sling"),
        ("barrier/keratin vs Slingshot", "mean_barrier_keratin", "mean_sling"),
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

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    elig_tumor = elig[elig["Sample_Origin"] != "nLung"]
    elig_nl_tl = elig[elig["Sample_Origin"].isin(["nLung", "tLung"])]
    sensitivity = [
        {"contrast": "tLung-only CLDN4 vs Slingshot", **sp(elig_t, "mean_CLDN4", "mean_sling")},
        {"contrast": "nLung-only CLDN4 vs Slingshot", **sp(elig_n, "mean_CLDN4", "mean_sling")},
        {"contrast": "tumor/met (drop nLung) CLDN4 vs Slingshot", **sp(elig_tumor, "mean_CLDN4", "mean_sling")},
        {"contrast": "nLung+tLung only CLDN4 vs Slingshot", **sp(elig_nl_tl, "mean_CLDN4", "mean_sling")},
        {"contrast": "tLung-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_t, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "tLung-only CLDN4 vs IFN", **sp(elig_t, "mean_CLDN4", "mean_IFN")},
        {"contrast": "tLung-only CLDN4 vs AT2", **sp(elig_t, "mean_CLDN4", "mean_AT2")},
        {"contrast": "CLDN4 vs DPT (sensitivity clock)", **sp(elig, "mean_CLDN4", "mean_dpt")},
        {"contrast": "nLung+tLung CLDN4 vs DPT (PR #325-like slice)", **sp(elig_nl_tl, "mean_CLDN4", "mean_dpt")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "Sample": sample,
                "Sample_Origin": origin,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "malignant_high": float(hi["score_malignant_like"].mean()),
                "malignant_low": float(lo["score_malignant_like"].mean()),
                "sling_high": float(hi["sling_pseudotime"].mean()),
                "sling_low": float(lo["sling_pseudotime"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN (no CLDN4) high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("malignant-like high vs low", "malignant_high", "malignant_low"),
        ("Slingshot high vs low", "sling_high", "sling_low"),
        ("DPT high vs low", "dpt_high", "dpt_low"),
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
        or len(paired_df) >= 4
    )

    cell_desc = {
        "CLDN4_vs_sling": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["sling_pseudotime"].to_numpy()),
        "CLDN4_vs_barrier": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["score_barrier_keratin"].to_numpy()),
        "CLDN4_vs_IFN": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["score_IFN"].to_numpy()),
        "CLDN4_vs_AT2": _spearman(adata.obs["expr_CLDN4"].to_numpy(), adata.obs["score_AT2"].to_numpy()),
        "note": "descriptive only; n_cells is not the experimental unit",
    }

    origin_colors = {
        "nLung": "#2a6f97",
        "tLung": "#b23a48",
        "tL/B": "#e09f3e",
        "mLN": "#6a4c93",
        "mBrain": "#335c67",
    }

    # Trajectory composite
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    ax = axes[0, 0]
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=ax,
        show=False,
        frameon=False,
        cmap="viridis",
        title="PAGA (Leiden) mean CLDN4",
    )
    ax = axes[0, 1]
    sc.pl.umap(adata, color="expr_CLDN4", ax=ax, show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    ax = axes[1, 0]
    sc.pl.umap(
        adata,
        color="sling_pseudotime",
        ax=ax,
        show=False,
        frameon=False,
        cmap="magma",
        title="UMAP Slingshot (AT2-rooted)",
    )
    ax = axes[1, 1]
    c4_sl = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot")
    for origin, col in origin_colors.items():
        sub = elig[elig["Sample_Origin"] == origin]
        if sub.empty:
            continue
        ax.scatter(sub["mean_sling"], sub["mean_CLDN4"], s=48, c=col, label=f"{origin} n={len(sub)}")
    ax.set_xlabel("sample-mean Slingshot")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_sl["rho"] is None else f"{c4_sl['rho']:.2f}"
    p_s = "NA" if c4_sl["p"] is None else f"{c4_sl['p']:.3g}"
    ax.set_title(f"CLDN4 vs Slingshot  n={c4_sl['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(
        f"GSE131907 epithelium/malignant  Slingshot  CLDN4   "
        f"n_cells={adata.n_obs}  n_samples={sample_df.shape[0]}  lineages={py_sling['n_lineages']}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # Honest n
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["Sample_Origin", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=ax, color=[origin_colors.get(c, "0.5") for c in ct.index])
    ax.set_ylabel("cells")
    ax.set_title(f"Honest n: author subtypes  n_cells={adata.n_obs}  n_samples={sample_df.shape[0]}")
    ax.legend(frameon=False, fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN", "IFN (no CLDN4)"),
            ("sling_low", "sling_high", "Slingshot", "Slingshot"),
        )
        for ax, (lo, hi, key, lab) in zip(axes, panels):
            for origin, col in origin_colors.items():
                sub = paired_df[paired_df["Sample_Origin"] == origin]
                if sub.empty:
                    continue
                ax.scatter(sub[lo], sub[hi], s=40, c=col, label=f"{origin} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
            prow = next(r for r in paired_rows if r["contrast"].startswith(key))
            ax.set_title(_fmt(prow, keys=("W", "p")))
        axes[0].legend(fontsize=6, frameon=False)
        fig.suptitle(f"EXTRA: within-sample CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    # Extra: CLDN4, barrier, IFN along Slingshot lineages
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), sharex=True)
    programs = (
        ("mean_CLDN4", "CLDN4"),
        ("mean_barrier_keratin", "barrier/keratin (no CLDN4)"),
        ("mean_IFN", "IFN (no CLDN4)"),
    )
    if not bin_df.empty:
        for ax, (col, lab) in zip(axes, programs):
            for lin, sub in bin_df.groupby("lineage"):
                ax.plot(sub["mean_pseudotime"], sub[col], marker="o", lw=1.4, label=lin)
            ax.set_xlabel("Slingshot (lineage bin mean)")
            ax.set_ylabel(lab)
            ax.set_title(lab)
        axes[0].legend(fontsize=7, frameon=False)
    else:
        for ax, (_, lab) in zip(axes, programs):
            ax.set_title(f"{lab} (no bins)")
    fig.suptitle(
        f"EXTRA: CLDN4 / barrier / IFN along Slingshot lineages  n_lineages={py_sling['n_lineages']}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_extra_lineage_programs")

    # Extra: sample-level CLDN4 vs AT2 / barrier / IFN
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    sample_panels = (
        ("mean_AT2", "AT2 score", next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")),
        ("mean_barrier_keratin", "barrier/keratin (no CLDN4)", barrier_row),
        ("mean_IFN", "IFN (no CLDN4)", next(r for r in primary if r["contrast"] == "CLDN4 vs IFN (no CLDN4)")),
    )
    for ax, (xcol, lab, row) in zip(axes, sample_panels):
        for origin, col in origin_colors.items():
            sub = elig[elig["Sample_Origin"] == origin]
            if sub.empty:
                continue
            ax.scatter(sub[xcol], sub["mean_CLDN4"], s=40, c=col, label=f"{origin} n={len(sub)}")
        ax.set_xlabel(f"sample-mean {lab}")
        ax.set_ylabel("sample-mean CLDN4")
        ax.set_title(_fmt(row))
    axes[0].legend(fontsize=6, frameon=False)
    fig.suptitle(f"EXTRA: sample-level CLDN4 vs programs  n={len(elig)}", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_programs")

    for color, fname, cmap in (
        ("Sample_Origin", "fig_umap_origin", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_IFN", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    n_patients = int(sample_df["patient_id"].nunique())
    other_origins = ~adata.obs["Sample_Origin"].isin(["nLung", "tLung"])

    c4_sl = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = barrier_row
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN (no CLDN4)")
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))
    tlung_c = next(r for r in sensitivity if r["contrast"] == "tLung-only CLDN4 vs Slingshot")

    holds = (
        f"**What holds (sample n={c4_sl['n']}).** "
        f"CLDN4 vs Slingshot {_fmt(c4_sl)}. "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded) {_fmt(c4_bar)}. "
        f"CLDN4 vs IFN {_fmt(c4_ifn)}. "
        f"Within-sample CLDN4-high vs low barrier {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN {_fmt(pair_ifn, keys=('W', 'p'))}."
    )
    parts = [
        f"Sample-level CLDN4 vs AT2-rooted Slingshot: {_fmt(c4_sl)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded from the score): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN (CLDN4 excluded): {_fmt(c4_ifn)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"tLung-only CLDN4 vs Slingshot: {_fmt(tlung_c)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"Slingshot has {py_sling['n_lineages']} lineage(s) from Leiden {start_cluster}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "The mixed nLung+tLung+met Slingshot correlation is not a within-tumor progression test.",
        "Not ICI. Not a TACSTD2 redo. No both-high gate. Not a redo of PR #325 DPT.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accession": "GSE131907",
        "histology": "LUAD",
        "public_only": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "not_a_pr325_dpt_redo": True,
        "n_cells": int(adata.n_obs),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_cells_tLung": int((adata.obs["Sample_Origin"] == "tLung").sum()),
        "n_cells_other": int(other_origins.sum()),
        "n_cells_catalog": int(extract_inv.get("n_epithelium_malignant_catalog") or adata.n_obs),
        "n_samples_catalog": int(extract_inv.get("n_samples_catalog") or sample_df.shape[0]),
        "catalog_by_origin": extract_inv.get("by_origin", {}),
        "cap_per_sample": int(extract_inv.get("cap_per_sample") or 400),
        "n_samples": int(sample_df.shape[0]),
        "n_samples_nLung": int((sample_df["Sample_Origin"] == "nLung").sum()),
        "n_samples_tLung": int((sample_df["Sample_Origin"] == "tLung").sum()),
        "n_samples_other": int((~sample_df["Sample_Origin"].isin(["nLung", "tLung"])).sum()),
        "n_samples_eligible": int(len(elig)),
        "n_patients": n_patients,
        "n_samples_paired_tertile": int(len(paired_df)),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "n_lineages": int(py_sling["n_lineages"]),
        "lineages": py_sling.get("lineages"),
        "leiden_resolution": LEIDEN_RES,
        "clock": engine,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "cell_level_descriptive": cell_desc,
        "what_holds": holds,
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
            "PR325": "DPT/PAGA nLung+tLung; this folder is REAL Slingshot + IFN + mets",
            "GSE205335": "winning-pair PR #449; not added",
            "GSE207422": "NSCLC mixed histology; not pooled",
            "PE": "unlabeled epithelium dropped",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    provenance = {
        "accession": "GSE131907",
        "paper": "Kim et al. Nat Commun 2020 PMID 32385277",
        "primary_gene": "CLDN4",
        "barrier_excludes_CLDN4": True,
        "ifn_excludes_CLDN4": True,
        "no_dual_high": True,
        "root": "nLung author AT2, never CLDN4-high",
        "clock": engine,
        "locked": {
            "leiden_resolution": LEIDEN_RES,
            "n_hvg": N_HVG,
            "n_neighbors": N_NEIGHBORS,
            "n_pcs": N_PCS,
            "sling_dims": SLING_DIMS,
        },
    }
    man = Path("/tmp/gse131907_slingshot_data/DOWNLOAD_MANIFEST.json")
    if man.is_file():
        provenance["download_manifest"] = json.loads(man.read_text())
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2, default=str))

    keep_cols = [
        c
        for c in [
            "Sample",
            "Sample_Origin",
            "Cell_type",
            "author_subtype",
            "leiden",
            "sling_pseudotime",
            "dpt_pseudotime",
            "cldn4_tertile",
            "expr_CLDN4",
            "expr_TACSTD2",
            "score_AT2",
            "score_barrier_keratin",
            "score_IFN",
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
                "n_lineages": py_sling["n_lineages"],
                "clock": engine,
                "extra": emit_extra,
                "finding": args.finding,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
