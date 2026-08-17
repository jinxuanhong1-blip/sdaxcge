#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE123902+GSE189357 marker epithelium, CLDN4 only.

Installs/uses Bioconductor slingshot (not a DPT fallback). PAGA is scanpy.
Root = GSE123902 NORMAL AT2-like cell, never CLDN4-high.
Inferential unit = the given PR #459 pair (n=22). Tails may be thin — said so.
No TACSTD2∩CLDN4 dual-high gate.
Done when results/tables/slingshot_lineages.tsv exists.
"""
from __future__ import annotations

import argparse
import json
import os
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
    COMPARATOR,
    CONTROLS,
    FOCAL,
    HALLMARK_INTERFERON_ALPHA_RESPONSE,
    QC_NEG,
    STATES,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e

HERE = Path(__file__).resolve().parents[1]
LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
CAP_PER_UNIT = 350


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


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def maybe_harmony(adata) -> dict:
    info = {"used": False, "reason": None}
    try:
        import harmonypy
    except ImportError:
        info["reason"] = "harmonypy not installed; neighbors on PCA"
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        return info
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :N_PCS],
        adata.obs,
        "dataset",
        max_iter_harmony=20,
    )
    z = np.asarray(ho.Z_corr)
    if z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z
    elif z.shape[1] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z.T
    else:
        raise ValueError(f"Harmony Z_corr shape {z.shape} vs n_obs={adata.n_obs}")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony")
    info["used"] = True
    info["reason"] = "harmonypy on PCA, batch=dataset"
    return info


def pick_root(adata) -> tuple[int, dict]:
    """GSE123902 NORMAL AT2-like cell. Never CLDN4-high."""
    normal = (adata.obs["dataset"].astype(str) == "GSE123902") & (
        adata.obs["tissue"].astype(str) == "NORMAL"
    )
    tert = adata.obs["cldn4_tertile"].astype(str)
    not_high = tert != "high"
    pool = normal & not_high
    info = {
        "n_normal": int(normal.sum()),
        "n_normal_not_cldn4_high": int(pool.sum()),
        "rule": None,
        "rejected_cldn4_high_root": True,
    }
    if int(pool.sum()) >= 10:
        idx = np.flatnonzero(pool.to_numpy())
        scores = adata.obs.loc[pool, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        info["rule"] = "GSE123902 NORMAL, not CLDN4-high, median AT2 score"
        return int(pick), info
    if int(normal.sum()) >= 5:
        idx = np.flatnonzero(normal.to_numpy())
        at2 = adata.obs.loc[normal, "score_AT2"].to_numpy()
        c4 = adata.obs.loc[normal, "expr_CLDN4"].to_numpy()
        # max AT2 among the lowest-CLDN4 half of NORMAL
        lo = c4 <= np.nanmedian(c4)
        cand = idx[lo] if lo.any() else idx
        scs = at2[lo] if lo.any() else at2
        pick = cand[int(np.nanargmax(scs))]
        info["rule"] = "GSE123902 NORMAL, lowest-CLDN4 half, max AT2 (thin NORMAL pool)"
        return int(pick), info
    # last resort: Leiden with max AT2, pick a CLDN4-low cell in that cluster
    means = (
        adata.obs.groupby("leiden", observed=True)["score_AT2"]
        .mean()
        .sort_values(ascending=False)
    )
    top = str(means.index[0])
    cand = adata.obs["leiden"].astype(str) == top
    cand = cand & not_high
    if not cand.any():
        cand = adata.obs["leiden"].astype(str) == top
    idx = np.flatnonzero(cand.to_numpy())
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"fallback Leiden {top} max AT2 among non-CLDN4-high (no NORMAL)"
    info["fallback_cluster"] = top
    return int(pick), info


def _r_env() -> dict[str, str]:
    env = dict(**os.environ)
    lib = Path.home() / "R" / "library"
    lib.mkdir(parents=True, exist_ok=True)
    env["R_LIBS_USER"] = str(lib)
    return env


def require_slingshot() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript not on PATH. Run scripts/install_r_slingshot.sh")
    proc = subprocess.run(
        [rscript, "-e", '.libPaths(Sys.getenv("R_LIBS_USER")); cat(as.character(packageVersion("slingshot")))'],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
        env=_r_env(),
    )
    if proc.returncode != 0:
        raise SystemExit(
            "Bioconductor slingshot is not installed. "
            "Run scripts/install_r_slingshot.sh\n"
            + (proc.stderr or proc.stdout or "")[:400]
        )
    return {"available": True, "version": (proc.stdout or "").strip(), "rscript": rscript}


def run_slingshot(adata, start_cluster: str, tabdir: Path, rscript: str) -> pd.DataFrame:
    scratch = tabdir / "_slingshot_input"
    scratch.mkdir(parents=True, exist_ok=True)
    rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
    pca = np.asarray(adata.obsm[rep][:, :N_PCS], dtype=float)
    cells = adata.obs_names.astype(str)
    pca_df = pd.DataFrame(pca, index=cells, columns=[f"PC{i+1}" for i in range(pca.shape[1])])
    pca_path = scratch / "pca.tsv"
    cl_path = scratch / "cluster.tsv"
    pca_df.to_csv(pca_path, sep="\t", index_label="cell_id")
    pd.DataFrame({"cell_id": cells, "cluster": adata.obs["leiden"].astype(str).to_numpy()}).to_csv(
        cl_path, sep="\t", index=False
    )
    r_script = Path(__file__).resolve().parent / "run_slingshot.R"
    cmd = [
        rscript,
        str(r_script),
        "--pca",
        str(pca_path),
        "--cluster",
        str(cl_path),
        "--start",
        str(start_cluster),
        "--outdir",
        str(tabdir),
    ]
    print("RUN", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=_r_env())
    sys.stdout.write(proc.stdout or "")
    sys.stderr.write(proc.stderr or "")
    if proc.returncode != 0:
        raise SystemExit(f"slingshot R failed with code {proc.returncode}")
    lin_path = tabdir / "slingshot_lineages.tsv"
    if not lin_path.exists():
        raise SystemExit("slingshot_lineages.tsv was not written")
    return pd.read_csv(lin_path, sep="\t")


def write_finding(path: Path, s: dict) -> None:
    def rho_s(r: dict) -> str:
        if r.get("rho") is None:
            return "NA"
        return f"{r['rho']:.3f}"

    def p_s(r: dict) -> str:
        if r.get("p") is None:
            return "NA"
        return f"{r['p']:.3g}"

    def q_s(r: dict) -> str:
        if r.get("q") is None:
            return "NA"
        return f"{r['q']:.3g}"

    lines = [
        "# Finding — pair GSE123902+GSE189357, CLDN4-only REAL Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high. No CellChat.",
        "",
        "The pair that already **differs** is **taken as given** and is not re-audited (PR #459):",
        "",
        "- **GSE123902 + GSE189357 %pos · n=22 · ρ=−0.638 · p=0.003 · I²=0%**",
        "- Members: GSE123902 marker-malignant donors n=13 + GSE189357 marker-malignant patients n=9",
        "- Between-patient Q4 vs Q1 on that same vector is **r=−1.000 on tails 7/5 — thin**. That is not the trajectory n.",
        "",
        "This folder asks a **different** question: where do **CLDN4**, a CLDN4-excluded **barrier/keratin** score, and an **IFN** score sit on a **real Slingshot** lineage (PAGA for geometry) of marker epithelium from this pair.",
        "",
        f"Primary clock: **Bioconductor slingshot {s['slingshot']['version']}** (Street et al. 2018). "
        "PAGA is scanpy (Wolf et al. 2019). Inferential unit = **patient/donor** on the given n=22. "
        "Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. "
        f"Root = {s['root']['rule']}. Root is **not** CLDN4-high.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Given combo (do not re-audit): **n=22** (13+9). Q4 vs Q1 tails **7/5 are thin**.",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/unit): **n_cells = {s['n_cells']}** "
        f"(GSE123902 {s['n_cells_gse123902']}, GSE189357 {s['n_cells_gse189357']}).",
        f"- Given tumor/met units in the object: **n_given_units = {s['n_given_units']}** "
        f"(GSE123902 {s['n_given_gse123902']}, GSE189357 {s['n_given_gse189357']}).",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} cells used for Spearman: **n = {s['n_units_eligible']}**.",
        f"- GSE123902 NORMAL root-pool cells / units: {s['n_cells_normal']} / {s['n_normal_units']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: available=True; version={s['slingshot']['version']}; "
        f"n_lineages={s['n_lineages']}; start_cluster={s['root']['start_cluster']}.",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Catalog before cap is in `results/tables/extract_inventory.json`. Thin libraries (e.g. LX699 46 malignant) stay in the given n=22 and are flagged, not dropped.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Slingshot / DPT root: {s['root']['rule']} (root cell {s['root'].get('root_unit')}, Leiden {s['root']['start_cluster']}). Never CLDN4-high.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN score: mean of Hallmark interferon-alpha response genes present (CLDN4 not in the set). Compact ISG panel is extra.",
        "- Marker-epithelial gate (same as PR #459): (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.",
        "",
        "## Lineage table (done criterion)",
        "",
        "Real Slingshot lineages. Start cluster is the root Leiden (NORMAL / not CLDN4-high). "
        "Cell-level ρ along each lineage is descriptive. Unit-level tests are in the Spearman table.",
        "",
        "| lineage | start | end | n_clusters | n_cells | path | ρ CLDN4 | ρ barrier | ρ IFN |",
        "| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for r in s["lineages"]:
        lines.append(
            f"| {r['lineage_id']} | {r['start_cluster']} | {r['end_cluster']} | "
            f"{r['n_clusters']} | {r['n_cells_finite_pt']} | `{r['path']}` | "
            f"{r.get('rho_CLDN4', 'NA')} | {r.get('rho_barrier', 'NA')} | {r.get('rho_IFN', 'NA')} |"
        )
    lines += [
        "",
        f"Table: `results/tables/slingshot_lineages.tsv`. Primary lineage for unit-level tests: **{s['primary_lineage']}**.",
        "",
        "## Primary (given-unit Spearman, BH inside this list)",
        "",
        "Patient/donor is the unit. Given n=22; Spearman n is units with ≥10 cells after QC/cap. "
        "Q4 vs Q1 tails 7/5 remain **thin** and are not this table.",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(f"| {r['contrast']} | {r['n']} | {rho_s(r)} | {p_s(r)} | {q_s(r)} |")
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_units | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        lines.append(f"| {r['contrast']} | {r['n']} | {rho_s(r)} | {p_s(r)} |")
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)",
        "",
        f"Paired tertile n={s['n_units_paired_tertile']}. Thin arms are possible on GSE123902 mets (LX699, LX701).",
        "",
        "| Paired contrast (high − low) | n_units | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("paired_tertile", []):
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        lines.append(f"| {r['contrast']} | {r['n']} | {d} | {p_s(r)} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- The given T/NK %pos Spearman (n=22, ρ=−0.638) is **not re-audited**.",
        "- Q4 vs Q1 tails **7/5 are thin**. Do not treat r=−1 as a trajectory result.",
        "- Marker-epithelial is **not** author malignant and **not** CNV.",
        "- GSE123902 NORMAL cells are a root pool, not extra inferential units.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot is an ordering on this object, not a developmental clock.",
        "- Not ICI / MPR / RECIST. Not CellChat.",
        "",
        "## Outputs",
        "",
        "- `results/tables/slingshot_lineages.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_extra_along_pseudotime.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/pair_123902_189357_slingshot_cldn4/requirements.txt",
        "bash methods/pair_123902_189357_slingshot_cldn4/scripts/install_r_slingshot.sh",
        "python3 methods/pair_123902_189357_slingshot_cldn4/scripts/download.py \\",
        "  --out /tmp/geo_pair_123902_189357",
        "python3 methods/pair_123902_189357_slingshot_cldn4/scripts/extract.py \\",
        "  --tars /tmp/geo_pair_123902_189357 \\",
        "  --out /tmp/geo_pair_123902_189357/epithelium.h5ad",
        "python3 methods/pair_123902_189357_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/geo_pair_123902_189357/epithelium.h5ad \\",
        "  --outdir methods/pair_123902_189357_slingshot_cldn4/results \\",
        "  --finding methods/pair_123902_189357_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default=str(HERE / "results"))
    p.add_argument("--finding", default=str(HERE / "FINDING.md"))
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling_meta = require_slingshot()
    print(json.dumps({"slingshot": sling_meta}, indent=2), flush=True)

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
    absent["IFN_hallmark_IFNa"] = _score_mean(
        adata, HALLMARK_INTERFERON_ALPHA_RESPONSE, "score_IFN"
    )
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
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"], ordered=True)

    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    harmony = maybe_harmony(adata)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)

    root_i, root_info = pick_root(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    root_info["index"] = root_i
    root_obs = adata.obs.iloc[root_i]
    root_info["root_unit"] = str(root_obs.get("unit_id", ""))
    root_info["root_tissue"] = str(root_obs.get("tissue", ""))
    root_info["root_dataset"] = str(root_obs.get("dataset", ""))
    root_info["start_cluster"] = str(root_obs["leiden"])
    root_info["root_CLDN4"] = float(root_obs["expr_CLDN4"])
    root_info["root_AT2"] = float(root_obs["score_AT2"])
    root_info["root_tertile"] = str(root_obs["cldn4_tertile"])
    if root_info["root_tertile"] == "high":
        raise SystemExit("root landed in CLDN4-high; refusing")

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
                "n_units": int(sub["unit_id"].nunique()),
                "n_GSE123902": int((sub["dataset"] == "GSE123902").sum()),
                "n_GSE189357": int((sub["dataset"] == "GSE189357").sum()),
                "n_NORMAL": int((sub["tissue"].astype(str) == "NORMAL").sum()),
                "frac_NORMAL": float((sub["tissue"].astype(str) == "NORMAL").mean()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    lin_df = run_slingshot(adata, root_info["start_cluster"], tabdir, sling_meta["rscript"])
    pt = pd.read_csv(tabdir / "slingshot_pseudotime.tsv", sep="\t")
    pt = pt.set_index("cell_id")
    lin_cols = [c for c in pt.columns if c != "cell_id"]
    for c in lin_cols:
        adata.obs[f"sling_{c}"] = pt.reindex(adata.obs_names)[c].to_numpy()
    # primary lineage = most cells with finite PT
    n_fin = {c: int(np.isfinite(adata.obs[f"sling_{c}"]).sum()) for c in lin_cols}
    primary_lin = max(n_fin, key=n_fin.get) if n_fin else None
    if primary_lin is None:
        raise SystemExit("slingshot returned no lineages")
    adata.obs["sling_pt"] = adata.obs[f"sling_{primary_lin}"]

    # enrich lineage table with CLDN4 / barrier / IFN along PT
    lin_rows = lin_df.to_dict("records")
    for r in lin_rows:
        lid = r["lineage_id"]
        col = f"sling_{lid}"
        if col not in adata.obs:
            r["rho_CLDN4"] = r["rho_barrier"] = r["rho_IFN"] = "NA"
            continue
        sub = adata.obs[np.isfinite(adata.obs[col])]
        r["mean_CLDN4"] = float(sub["expr_CLDN4"].mean()) if len(sub) else None
        r["mean_barrier"] = float(sub["score_barrier_keratin"].mean()) if len(sub) else None
        r["mean_IFN"] = float(sub["score_IFN"].mean()) if len(sub) else None
        r["mean_AT2"] = float(sub["score_AT2"].mean()) if len(sub) else None
        for key, src in (
            ("rho_CLDN4", "expr_CLDN4"),
            ("rho_barrier", "score_barrier_keratin"),
            ("rho_IFN", "score_IFN"),
            ("rho_AT2", "score_AT2"),
        ):
            sp = _spearman(sub[src].to_numpy(), sub[col].to_numpy())
            r[key] = None if sp["rho"] is None else round(sp["rho"], 3)
            r[f"p_{key}"] = sp["p"]
            r[f"n_{key}"] = sp["n"]
    lin_out = pd.DataFrame(lin_rows)
    lin_out.to_csv(tabdir / "slingshot_lineages.tsv", sep="\t", index=False)

    rows = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        rows.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "patient": str(sub["patient"].iloc[0]),
                "tissue": str(sub["tissue"].iloc[0]),
                "role": str(sub["role"].iloc[0]),
                "n_cells": int(len(sub)),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "thin_arm": bool(
                    min(
                        int((sub["cldn4_tertile"] == "low").sum()),
                        int((sub["cldn4_tertile"] == "high").sum()),
                    )
                    < 20
                    and str(sub["role"].iloc[0]) == "given"
                ),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_sling_pt": float(sub["sling_pt"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_ifn_isg": float(sub["score_ifn_isg"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    if "expr_SFTPC" in adata.obs:
        sft_map = {
            unit: float(sub["expr_SFTPC"].mean())
            for unit, sub in adata.obs.groupby("unit_id", observed=True)
        }
        sample_df["mean_SFTPC"] = sample_df["unit_id"].map(sft_map)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)

    given = sample_df[sample_df["role"] == "given"].copy()
    elig = given[given["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling_pt"),
        ("CLDN4 vs DPT (companion)", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN (Hallmark IFNα)", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs compact ISG (extra)", "mean_CLDN4", "mean_ifn_isg"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("barrier vs Slingshot PT", "mean_barrier_keratin", "mean_sling_pt"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_sling_pt"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling_pt"),
        ("AT2 score vs Slingshot PT (control)", "mean_AT2", "mean_sling_pt"),
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

    elig_123 = elig[elig["dataset"] == "GSE123902"]
    elig_189 = elig[elig["dataset"] == "GSE189357"]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE123902-only CLDN4 vs Slingshot PT", **sp(elig_123, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "GSE189357-only CLDN4 vs Slingshot PT", **sp(elig_189, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "GSE123902-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_123, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE189357-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_189, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE123902-only CLDN4 vs IFN", **sp(elig_123, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE189357-only CLDN4 vs IFN", **sp(elig_189, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE123902-only barrier vs Slingshot PT", **sp(elig_123, "mean_barrier_keratin", "mean_sling_pt")},
        {"contrast": "GSE189357-only barrier vs Slingshot PT", **sp(elig_189, "mean_barrier_keratin", "mean_sling_pt")},
        {"contrast": "GSE123902-only IFN vs Slingshot PT", **sp(elig_123, "mean_IFN", "mean_sling_pt")},
        {"contrast": "GSE189357-only IFN vs Slingshot PT", **sp(elig_189, "mean_IFN", "mean_sling_pt")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        if str(sub["role"].iloc[0]) != "given":
            continue
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "thin_arm": min(len(hi), len(lo)) < 20,
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "pt_high": float(hi["sling_pt"].mean()),
                "pt_low": float(lo["sling_pt"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN (Hallmark IFNα) high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("Slingshot PT high vs low", "pt_high", "pt_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    def grab(name: str) -> dict:
        return next(r for r in primary if r["contrast"] == name)

    c4_pt = grab("CLDN4 vs Slingshot PT")
    c4_bar = grab("CLDN4 vs barrier/keratin (no CLDN4)")
    c4_ifn = grab("CLDN4 vs IFN (Hallmark IFNα)")
    bar_pt = grab("barrier vs Slingshot PT")
    ifn_pt = grab("IFN vs Slingshot PT")
    at2_pt = grab("AT2 score vs Slingshot PT (control)")

    holds = []
    if c4_bar.get("p") is not None and c4_bar["p"] < 0.05:
        holds.append(
            f"CLDN4 tracks CLDN4-excluded barrier/keratin (n={c4_bar['n']}, ρ={c4_bar['rho']:.3f}, p={c4_bar['p']:.3g})."
        )
    if c4_ifn.get("p") is not None and c4_ifn["p"] < 0.05:
        holds.append(
            f"CLDN4 vs IFN is nonzero (n={c4_ifn['n']}, ρ={c4_ifn['rho']:.3f}, p={c4_ifn['p']:.3g})."
        )
    if at2_pt.get("p") is not None and at2_pt["p"] < 0.05:
        holds.append(
            f"AT2 vs Slingshot PT control is nonzero (n={at2_pt['n']}, ρ={at2_pt['rho']:.3f}, p={at2_pt['p']:.3g})."
        )
    fails = []
    if c4_pt.get("p") is None or c4_pt["p"] >= 0.05:
        rho_txt = "NA" if c4_pt.get("rho") is None else f"{c4_pt['rho']:.3f}"
        p_txt = "NA" if c4_pt.get("p") is None else f"{c4_pt['p']:.3g}"
        fails.append(
            f"Pooled given-unit CLDN4 vs Slingshot PT is null "
            f"(n={c4_pt['n']}, ρ={rho_txt}, p={p_txt})."
        )
    what_holds = (
        f"**What holds (given units, Spearman n={elig.shape[0]}).** " + " ".join(holds)
        if holds
        else f"**What holds (given units, Spearman n={elig.shape[0]}).** No primary contrast is significant at p<0.05."
    )
    what_fails = (
        "**What does not hold / is thin.** " + " ".join(fails)
        + " Q4 vs Q1 tails 7/5 from PR #459 stay thin and are not re-used as a trajectory test."
    )
    what_holds = what_holds + " " + what_fails

    def fmt(r: dict) -> str:
        if r.get("rho") is None:
            return f"n={r['n']}, ρ=NA, p=NA"
        return f"n={r['n']}, ρ={r['rho']:.3f}, p={r['p']:.3g}"

    verdict = (
        f"Real Slingshot ({sling_meta['version']}) produced {len(lin_rows)} lineage(s) "
        f"rooted on Leiden {root_info['start_cluster']} ({root_info['rule']}; "
        f"root tertile={root_info['root_tertile']}, CLDN4={root_info['root_CLDN4']:.3f}). "
        f"Primary lineage {primary_lin}. "
        f"Given-unit CLDN4 vs Slingshot PT: {fmt(c4_pt)}. "
        f"CLDN4 vs barrier/keratin (no CLDN4): {fmt(c4_bar)}. "
        f"CLDN4 vs IFN (Hallmark IFNα): {fmt(c4_ifn)}. "
        f"Barrier vs PT: {fmt(bar_pt)}. IFN vs PT: {fmt(ifn_pt)}. "
        f"AT2 vs PT (control): {fmt(at2_pt)}. "
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices. "
        f"Given combo n=22 is not re-audited. Q4 vs Q1 tails 7/5 are thin. "
        "No dual-high gate. Not a T/NK redo."
    )

    # ---- figures ----
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.2))
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
        color="sling_pt",
        ax=ax,
        show=False,
        frameon=False,
        cmap="magma",
        title=f"UMAP Slingshot {primary_lin} (root not CLDN4-high)",
    )
    ax = axes[1, 1]
    colors = {"GSE123902": "#2a6f97", "GSE189357": "#b23a48"}
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_sling_pt"], sub["mean_CLDN4"], s=52, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("unit-mean Slingshot PT")
    ax.set_ylabel("unit-mean CLDN4")
    ax.set_title(f"CLDN4 vs Slingshot PT  {fmt(c4_pt)}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Pair GSE123902+GSE189357  CLDN4-only Slingshot/PAGA   "
        f"n_cells={adata.n_obs}  given_units={given.shape[0]}  (PR #459 n=22)",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # extra: CLDN4 + barrier + IFN along PT (cell hex + unit)
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 7.6))
    cell_pt = adata.obs["sling_pt"].to_numpy()
    for j, (title, col, cmap) in enumerate(
        (
            ("CLDN4", "expr_CLDN4", "viridis"),
            ("barrier/keratin (no CLDN4)", "score_barrier_keratin", "cividis"),
            ("IFN Hallmark IFNα", "score_IFN", "plasma"),
        )
    ):
        ax = axes[0, j]
        sc.pl.umap(adata, color=col, ax=ax, show=False, frameon=False, cmap=cmap, title=title)
        ax = axes[1, j]
        y = adata.obs[col].to_numpy()
        m = np.isfinite(cell_pt) & np.isfinite(y)
        ax.hexbin(cell_pt[m], y[m], gridsize=40, cmap="Greys", mincnt=1)
        for ds, c in colors.items():
            sub = elig[elig["dataset"] == ds]
            xkey = "mean_sling_pt"
            ykey = {
                "expr_CLDN4": "mean_CLDN4",
                "score_barrier_keratin": "mean_barrier_keratin",
                "score_IFN": "mean_IFN",
            }[col]
            ax.scatter(sub[xkey], sub[ykey], s=40, c=c, edgecolor="white", linewidth=0.4, label=ds)
        ax.set_xlabel("Slingshot PT")
        ax.set_ylabel(title)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("CLDN4 + barrier + IFN along Slingshot PT (hex = cells; points = given units)", fontsize=11)
    _save(fig, figdir / "fig_extra_along_pseudotime")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.0))
    if not paired_df.empty:
        for ax, (lab, a, b) in zip(
            axes,
            (
                ("barrier (no CLDN4)", "barrier_high", "barrier_low"),
                ("IFN", "IFN_high", "IFN_low"),
                ("Slingshot PT", "pt_high", "pt_low"),
            ),
        ):
            ax.scatter(paired_df[b], paired_df[a], s=40, c="#3d405b")
            lim = [
                min(paired_df[a].min(), paired_df[b].min()),
                max(paired_df[a].max(), paired_df[b].max()),
            ]
            ax.plot(lim, lim, color="0.7", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
            ax.set_title(f"paired n={len(paired_df)}")
    else:
        for ax in axes:
            ax.text(0.5, 0.5, "no paired units", ha="center")
            ax.axis("off")
    fig.suptitle("Extra: within-unit CLDN4-high vs low (thin arms flagged in table)", fontsize=11)
    _save(fig, figdir / "fig_extra_cldn4_tertile")

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    ax = axes[0]
    labels = [
        "given combo\n(PR #459)",
        "Q4 vs Q1\ntails (thin)",
        "units in\nobject",
        f"Spearman\n(≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} cells)",
        "paired\ntertile",
    ]
    vals = [22, 12, int(given.shape[0]), int(elig.shape[0]), int(len(paired_df))]
    ax.bar(labels, vals, color=["#2a6f97", "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f"])
    ax.set_ylabel("n (patients/donors)")
    ax.set_title("Honest n — unit is patient/donor")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.3, str(v), ha="center", fontsize=9)
    ax.set_ylim(0, max(vals) + 4)
    ax = axes[1]
    by = (
        adata.obs.groupby(["dataset", "role", "tissue"], observed=True)
        .size()
        .reset_index(name="n")
    )
    ax.barh(
        [f"{r.dataset} {r.role} {r.tissue}" for r in by.itertuples()],
        by["n"],
        color="#4a4e69",
    )
    ax.set_xlabel("cells after cap/QC")
    ax.set_title(f"n_cells={adata.n_obs} (cap ≤{CAP_PER_UNIT}/unit)")
    fig.suptitle("Tails may be thin — Q4 vs Q1 is 7/5, not the Slingshot n", fontsize=11)
    _save(fig, figdir / "fig_honest_n")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    sc.pl.umap(adata, color="dataset", ax=axes[0], show=False, frameon=False, title="dataset")
    sc.pl.umap(adata, color="tissue", ax=axes[1], show=False, frameon=False, title="tissue")
    sc.pl.umap(adata, color="score_AT2", ax=axes[2], show=False, frameon=False, cmap="YlGn", title="AT2 score")
    _save(fig, figdir / "fig_extra_umap_batch_at2")

    summary = {
        "given_combo": {
            "pr": 459,
            "n": 22,
            "rho_pctpos": -0.638,
            "p": 0.003,
            "q4q1_tails": "7/5 thin",
        },
        "n_cells": int(adata.n_obs),
        "n_cells_gse123902": int((adata.obs["dataset"] == "GSE123902").sum()),
        "n_cells_gse189357": int((adata.obs["dataset"] == "GSE189357").sum()),
        "n_cells_normal": int((adata.obs["tissue"].astype(str) == "NORMAL").sum()),
        "n_normal_units": int(
            adata.obs.loc[adata.obs["tissue"].astype(str) == "NORMAL", "unit_id"].nunique()
        ),
        "n_given_units": int(given.shape[0]),
        "n_given_gse123902": int((given["dataset"] == "GSE123902").sum()),
        "n_given_gse189357": int((given["dataset"] == "GSE189357").sum()),
        "n_units_eligible": int(elig.shape[0]),
        "n_units_paired_tertile": int(len(paired_df)),
        "n_leiden": int(len(leiden_ids)),
        "n_paga_components": int(len(comps)),
        "n_lineages": int(len(lin_rows)),
        "primary_lineage": primary_lin,
        "cap_per_unit": CAP_PER_UNIT,
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].value_counts().to_dict(),
        "genes_absent": absent,
        "harmony": harmony,
        "root": root_info,
        "slingshot": sling_meta,
        "lineages": lin_rows,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "what_holds": what_holds,
        "verdict": verdict,
        "dual_high": False,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), summary)
    print(f"DONE lineage table {tabdir / 'slingshot_lineages.tsv'}", flush=True)
    print(verdict, flush=True)


if __name__ == "__main__":
    main()
