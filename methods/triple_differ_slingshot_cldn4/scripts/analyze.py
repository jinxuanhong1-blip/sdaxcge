#!/usr/bin/env python3
"""REAL Slingshot/PAGA on the triple that differs, CLDN4 only.

ADDITIVE. GSE123902 + GSE131907 + GSE205335 epithelium.
Do not add GSE148071. This is not the 7-pool. No TACSTD2∩CLDN4 dual-high.

Try Harmony on the triple. If Harmony dies, per-dataset graphs, then stack
patient/sample CLDN4 vs that graph's pseudotime.

Primary clock: Slingshot (Street 2018) if R/slingshot is installed; else
documented AT2-rooted scanpy DPT. Root is never CLDN4-high.
Inferential unit = GSE123902 donor-sample, GSE131907 Sample, GSE205335 patient.
Done when stacked_sample_cldn4_pseudotime.tsv exists.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
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
    raise SystemExit("scanpy is required — pip install -r requirements.txt") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05
COHORTS = ("GSE123902", "GSE131907", "GSE205335")
COHORT_COLORS = {
    "GSE123902": "#4C78A8",
    "GSE131907": "#F58518",
    "GSE205335": "#54A24B",
}
NORMAL_ORIGINS = {"nLung", "NORMAL"}
SLINGSHOT_R = HERE / "run_slingshot.R"


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "available": False,
            "reason": "Rscript not on PATH",
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    try:
        proc = subprocess.run(
            [rscript, "-e", 'packageVersion("slingshot")'],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "available": False,
            "reason": f"Rscript slingshot probe failed: {exc}",
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:400],
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


def palantir_status() -> dict:
    try:
        import palantir  # noqa: F401

        return {"available": True, "module": "palantir"}
    except ImportError:
        return {
            "available": False,
            "reason": "python package palantir not installed",
            "note": "optional companion; Slingshot/DPT remains the primary clock",
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


def _is_normal(obs: pd.DataFrame) -> pd.Series:
    origin = obs["Sample_Origin"].astype(str)
    flag = obs.get("is_normal_tissue", pd.Series("False", index=obs.index)).astype(str)
    return origin.isin(NORMAL_ORIGINS) | flag.eq("True")


def _is_cldn4_high(obs: pd.DataFrame) -> pd.Series:
    return obs["cldn4_tertile"].astype(str).eq("high")


def _pick_root(adata, dataset: str | None = None) -> tuple[int, dict]:
    """External arrow. Never root on CLDN4-high."""
    obs = adata.obs
    ds = obs["dataset"].astype(str)
    subtype = obs["author_subtype"].astype(str)
    normal = _is_normal(obs)
    not_high = ~_is_cldn4_high(obs)
    info: dict = {"rule": None, "dataset_scope": dataset or "joint"}

    nlung_at2 = (
        ds.eq("GSE131907")
        & obs["Sample_Origin"].astype(str).eq("nLung")
        & subtype.isin(AUTHOR_AT2)
        & not_high
    )
    info["n_nLung_author_AT2_not_CLDN4high"] = int(nlung_at2.sum())
    if int(nlung_at2.sum()) >= 10:
        idx = np.flatnonzero(nlung_at2.to_numpy())
        scores = obs.loc[nlung_at2, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        info["rule"] = "GSE131907 nLung author AT2 (median AT2; not CLDN4-high)"
        return int(pick), info

    # GSE123902 NORMAL AT2-like, not CLDN4-high
    n123 = ds.eq("GSE123902") & obs["Sample_Origin"].astype(str).eq("NORMAL") & not_high
    info["n_GSE123902_NORMAL_not_CLDN4high"] = int(n123.sum())
    if int(n123.sum()) >= 10:
        idx = np.flatnonzero(n123.to_numpy())
        scores = obs.loc[n123, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        info["rule"] = "GSE123902 NORMAL marker-epi (median AT2; not CLDN4-high)"
        return int(pick), info

    # per-dataset / last resort: highest AT2 among not-CLDN4-high
    if "leiden" not in obs:
        raise SystemExit("leiden missing before root pick")
    pool = not_high
    if dataset is not None:
        pool = pool & ds.eq(dataset)
    if int(pool.sum()) < 5:
        pool = ds.eq(dataset) if dataset is not None else pd.Series(True, index=obs.index)
    sub = obs.loc[pool, ["leiden", "score_AT2"]]
    if sub.empty:
        raise SystemExit("no cells left to root (after excluding CLDN4-high)")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = pool & (obs["leiden"].astype(str) == top)
    scores = obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"Leiden {top} max AT2 among not-CLDN4-high (fallback)"
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


def try_harmony(adata) -> dict:
    info = {"used": False, "reason": None, "died": False}
    try:
        import harmonypy
    except ImportError:
        info["died"] = True
        info["reason"] = "harmonypy not installed"
        return info
    try:
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
        if not np.isfinite(adata.obsm["X_pca_harmony"]).all():
            raise ValueError("Harmony produced non-finite values")
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony")
        info["used"] = True
        info["reason"] = "harmonypy on PCA, batch=dataset"
        return info
    except Exception as exc:  # Harmony death → per-dataset graphs
        info["died"] = True
        info["reason"] = f"Harmony died: {type(exc).__name__}: {exc}"
        print(json.dumps({"harmony_died": info["reason"]}), flush=True)
        return info


def _compute_hvg_pca(adata) -> None:
    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
    except Exception:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")


def _leiden_paga_umap(adata) -> None:
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)


def run_slingshot_on(adata, root_i: int, sling: dict) -> dict:
    """Write slingshot_pseudotime onto adata.obs if R/slingshot works."""
    out = {"ran": False, "reason": None}
    if not sling.get("available"):
        out["reason"] = sling.get("reason")
        adata.obs["slingshot_pseudotime"] = np.nan
        return out
    if not SLINGSHOT_R.is_file():
        out["reason"] = f"missing {SLINGSHOT_R}"
        adata.obs["slingshot_pseudotime"] = np.nan
        return out
    start = str(adata.obs.iloc[root_i]["leiden"])
    if "X_pca_harmony" in adata.obsm:
        emb = np.asarray(adata.obsm["X_pca_harmony"][:, :N_PCS], dtype=float)
    else:
        emb = np.asarray(adata.obsm["X_pca"][:, :N_PCS], dtype=float)
    with tempfile.TemporaryDirectory(prefix="slingshot-") as tmp:
        tmp_p = Path(tmp)
        emb_df = pd.DataFrame(emb, index=adata.obs_names.astype(str))
        emb_path = tmp_p / "embedding.csv"
        clu_path = tmp_p / "clusters.csv"
        out_path = tmp_p / "pt.csv"
        emb_df.to_csv(emb_path)
        pd.DataFrame(
            {"cell": adata.obs_names.astype(str), "cluster": adata.obs["leiden"].astype(str)}
        ).to_csv(clu_path, index=False)
        rscript = shutil.which("Rscript")
        try:
            proc = subprocess.run(
                [rscript, str(SLINGSHOT_R), str(emb_path), str(clu_path), start, str(out_path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=600,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            out["reason"] = f"slingshot run failed: {exc}"
            adata.obs["slingshot_pseudotime"] = np.nan
            return out
        if proc.returncode != 0 or not out_path.is_file():
            out["reason"] = (proc.stderr or proc.stdout or "slingshot failed").strip()[:500]
            adata.obs["slingshot_pseudotime"] = np.nan
            return out
        pt = pd.read_csv(out_path)
        mapped = pt.set_index("cell")["slingshot_pseudotime"].reindex(adata.obs_names.astype(str))
        adata.obs["slingshot_pseudotime"] = mapped.to_numpy(dtype=float)
        out["ran"] = True
        out["start_cluster"] = start
        out["n_finite"] = int(np.isfinite(adata.obs["slingshot_pseudotime"]).sum())
        out["stdout"] = (proc.stdout or "").strip()[:300]
        return out


def assign_clock(adata, sling_run: dict) -> str:
    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    if sling_run.get("ran") and int(np.isfinite(adata.obs["slingshot_pseudotime"]).sum()) >= 20:
        adata.obs["pseudotime"] = adata.obs["slingshot_pseudotime"]
        return "slingshot"
    adata.obs["pseudotime"] = adata.obs["dpt_pseudotime"]
    return "scanpy_dpt"


def graph_one(adata, *, try_h: bool, sling: dict, dataset: str | None) -> dict:
    """HVG→PCA→(Harmony?)→neighbors→Leiden→PAGA→DPT→optional Slingshot."""
    _compute_hvg_pca(adata)
    harmony = {"used": False, "died": False, "reason": "not requested (per-dataset graph)"}
    if try_h:
        harmony = try_harmony(adata)
        if harmony.get("died"):
            return {"ok": False, "harmony": harmony}
        # try_harmony already set neighbors
    else:
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    _leiden_paga_umap(adata)
    root_i, root_info = _pick_root(adata, dataset=dataset)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_unit"] = str(adata.obs.iloc[root_i].get("unit_id", ""))
    root_info["root_origin"] = str(adata.obs.iloc[root_i].get("Sample_Origin", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("author_subtype", ""))
    root_info["root_dataset"] = str(adata.obs.iloc[root_i].get("dataset", ""))
    root_info["root_cldn4_tertile"] = str(adata.obs.iloc[root_i].get("cldn4_tertile", ""))
    if root_info["root_cldn4_tertile"] == "high":
        raise SystemExit("root landed on CLDN4-high; refusing")
    sling_run = run_slingshot_on(adata, root_i, sling)
    clock = assign_clock(adata, sling_run)
    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    return {
        "ok": True,
        "harmony": harmony,
        "root": root_info,
        "slingshot_run": sling_run,
        "clock": clock,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(adata.obs["leiden"].nunique()),
        "connect": connect,
        "leiden_ids": [str(x) for x in adata.obs["leiden"].cat.categories],
    }


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

    sling = s["slingshot"]
    clock = (
        "Slingshot (Street 2018)"
        if s.get("clock") == "slingshot"
        else "documented AT2-rooted diffusion pseudotime (scanpy DPT; Slingshot R missing or failed)"
    )
    graph = s.get("graph_mode")
    lines = [
        "# Finding — triple that differs (GSE123902+GSE131907+GSE205335), CLDN4-only trajectory",
        "",
        "ADDITIVE. **CLDN4 only.** The triple that **differs** is taken as given from PR #459 "
        "(author %pos GSE123902+GSE131907+GSE205335, n=56, ρ=−0.522 vs T/NK). "
        "This folder does **not** redo that T/NK Spearman. "
        "Do **not** add GSE148071. This is **not** the 7-pool. "
        "No TACSTD2∩CLDN4 dual-high gate.",
        "",
        f"Primary clock: **{clock}**. Graph mode: **{graph}**. "
        "Inferential unit = **sample/patient/donor-sample**. "
        "Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        "Root is never CLDN4-high (prefer GSE131907 nLung author AT2).",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/unit): **n_cells = {s['n_cells']}** "
        f"(GSE123902 {s['n_cells_gse123902']}, GSE131907 {s['n_cells_gse131907']}, "
        f"GSE205335 {s['n_cells_gse205335']}).",
        f"- Units (GSE123902 donor-sample + GSE131907 Sample + GSE205335 patient): "
        f"**n_units = {s['n_units']}** "
        f"(GSE123902 {s['n_units_gse123902']}, GSE131907 {s['n_units_gse131907']}, "
        f"GSE205335 {s['n_units_gse205335']}).",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: {s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- GSE123902 NORMAL cells: {s['n_cells_gse123902_normal']}.",
        f"- Author/marker subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE148071 not used. 7-pool extras not used. GSE131907 PE unlabeled epithelium dropped. "
        "GSE205335 normal-tissue samples dropped. GSE123902 36.5 GB H5 skipped; marker epithelium used.",
        f"- Graph: {graph}. Harmony: {s['harmony'].get('reason')}.",
        f"- Slingshot: available={sling.get('available')}; "
        f"{sling.get('reason') or sling.get('version')}. ran={s.get('slingshot_run', {}).get('ran')}.",
        f"- Palantir: available={s['palantir'].get('available')}; "
        f"{s['palantir'].get('reason') or s['palantir'].get('module')}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Root: {s['root']['rule']} (root cell index {s['root']['index']}, unit {s['root'].get('root_unit')}, "
        f"CLDN4 tertile {s['root'].get('root_cldn4_tertile')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "",
        "## Primary (stacked sample-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_units | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<{EXTRA_BARRIER_P_LT}, "
            f"or any paired tertile Wilcoxon p<{EXTRA_BARRIER_P_LT}, or n_paired≥4. "
            f"Observed Spearman n={extra['spearman_n']}, ρ={extra_rho}, p={extra_p}."
        ),
        "",
        "| Paired contrast (high − low) | n_units | Δ median | p |",
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
        "- The PR #459 n=56 T/NK Spearman is given and is **not** re-audited here.",
        "- The stacked pseudotime Spearman mixes cohorts and normal vs tumor and is **not** a within-tumor progression test.",
        "- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.",
        "- GSE123902 epithelium is marker-based (36.5 GB author H5 skipped), not author cell types.",
        "- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- GSE148071 was not added. This is not the 7-pool.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Pseudotime is an ordering, not a clock, unless Slingshot actually ran.",
        "",
        "## Outputs",
        "",
        "- `results/tables/stacked_sample_cldn4_pseudotime.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_stacked_cldn4_pseudotime.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/triple_differ_slingshot_cldn4/requirements.txt",
        "python3 methods/triple_differ_slingshot_cldn4/scripts/download.py \\",
        "  --outdir /tmp/triple_differ_raw",
        "python3 methods/triple_differ_slingshot_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/triple_differ_raw \\",
        "  --out /tmp/triple_differ_raw/epithelium.h5ad",
        "python3 methods/triple_differ_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/triple_differ_raw/epithelium.h5ad \\",
        "  --outdir methods/triple_differ_slingshot_cldn4/results \\",
        "  --finding methods/triple_differ_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/triple_differ_slingshot_cldn4/results")
    p.add_argument("--finding", default="methods/triple_differ_slingshot_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling = slingshot_status()
    pal = palantir_status()
    print(json.dumps({"slingshot": sling, "palantir": pal}, indent=2), flush=True)

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

    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"], ordered=True)

    graph_mode = "harmony_joint"
    per_dataset_info: dict[str, dict] = {}

    joint = adata.copy()
    ginfo = graph_one(joint, try_h=True, sling=sling, dataset=None)
    if not ginfo["ok"]:
        print("Harmony died — falling back to per-dataset graphs", flush=True)
        graph_mode = "per_dataset"
        adata.obs["pseudotime"] = np.nan
        adata.obs["dpt_pseudotime"] = np.nan
        adata.obs["slingshot_pseudotime"] = np.nan
        adata.obs["leiden"] = "NA"
        clocks = []
        roots = []
        n_comp = 0
        n_leid = 0
        last_harmony = ginfo["harmony"]
        last_sling_run = {"ran": False, "reason": "per-dataset"}
        for ds in COHORTS:
            mask = adata.obs["dataset"].astype(str).eq(ds)
            if int(mask.sum()) < 50:
                print(f"skip {ds}: n={int(mask.sum())}", flush=True)
                continue
            sub = adata[mask].copy()
            # restore log-normalized X from .X (already log1p)
            info = graph_one(sub, try_h=False, sling=sling, dataset=ds)
            if not info["ok"]:
                raise SystemExit(f"per-dataset graph failed for {ds}: {info}")
            per_dataset_info[ds] = {
                k: info[k]
                for k in ("harmony", "root", "slingshot_run", "clock", "n_paga_components", "n_leiden")
            }
            clocks.append(info["clock"])
            roots.append(info["root"])
            n_comp += info["n_paga_components"]
            n_leid += info["n_leiden"]
            last_sling_run = info["slingshot_run"]
            # write back
            adata.obs.loc[sub.obs_names, "pseudotime"] = sub.obs["pseudotime"].to_numpy()
            adata.obs.loc[sub.obs_names, "dpt_pseudotime"] = sub.obs["dpt_pseudotime"].to_numpy()
            if "slingshot_pseudotime" in sub.obs:
                adata.obs.loc[sub.obs_names, "slingshot_pseudotime"] = sub.obs[
                    "slingshot_pseudotime"
                ].to_numpy()
            adata.obs.loc[sub.obs_names, "leiden"] = (ds + ":" + sub.obs["leiden"].astype(str)).to_numpy()
            # keep last UMAP coords in a per-dataset slot
            adata.obsm.setdefault(f"X_umap_{ds}", np.full((adata.n_obs, 2), np.nan))
            um = np.asarray(sub.obsm["X_umap"])
            idx = adata.obs_names.get_indexer(sub.obs_names)
            coords = np.full((adata.n_obs, 2), np.nan)
            coords[idx] = um
            adata.obsm[f"X_umap_{ds}"] = coords
            # figures per dataset
            fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
            sc.pl.paga(
                sub,
                color="expr_CLDN4",
                ax=axes[0],
                show=False,
                frameon=False,
                cmap="viridis",
                title=f"{ds} PAGA CLDN4",
            )
            sc.pl.umap(sub, color="expr_CLDN4", ax=axes[1], show=False, frameon=False, cmap="viridis", title=f"{ds} UMAP CLDN4")
            sc.pl.umap(
                sub,
                color="pseudotime",
                ax=axes[2],
                show=False,
                frameon=False,
                cmap="magma",
                title=f"{ds} {info['clock']}",
            )
            _save(fig, figdir / f"fig_perdataset_{ds}")
            pd.DataFrame(info["connect"], index=info["leiden_ids"], columns=info["leiden_ids"]).to_csv(
                tabdir / f"paga_connectivities_{ds}.tsv", sep="\t"
            )
        clock = "slingshot" if clocks and all(c == "slingshot" for c in clocks) else (
            "mixed" if any(c == "slingshot" for c in clocks) else "scanpy_dpt"
        )
        # representative root: first GSE131907 if present
        root_info = next((r for r in roots if r.get("root_dataset") == "GSE131907"), roots[0] if roots else {})
        harmony = last_harmony
        sling_run = last_sling_run
        n_paga_components = n_comp
        n_leiden = n_leid
        connect = None
        leiden_ids = sorted(adata.obs["leiden"].astype(str).unique().tolist())
        # dummy UMAP for joint-style plots: concatenate per-dataset UMAPs with offsets
        umap = np.zeros((adata.n_obs, 2))
        for i, ds in enumerate(COHORTS):
            key = f"X_umap_{ds}"
            if key in adata.obsm:
                block = adata.obsm[key]
                finite = np.isfinite(block[:, 0])
                umap[finite] = block[finite] + np.array([i * 12.0, 0.0])
        adata.obsm["X_umap"] = umap
        adata.obs["leiden"] = pd.Categorical(adata.obs["leiden"].astype(str))
    else:
        adata = joint
        clock = ginfo["clock"]
        root_info = ginfo["root"]
        harmony = ginfo["harmony"]
        sling_run = ginfo["slingshot_run"]
        n_paga_components = ginfo["n_paga_components"]
        n_leiden = ginfo["n_leiden"]
        connect = ginfo["connect"]
        leiden_ids = ginfo["leiden_ids"]

    if connect is not None:
        cluster_tab = []
        for cl in leiden_ids:
            sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
            cluster_tab.append(
                {
                    "leiden": cl,
                    "n_cells": int(len(sub)),
                    "n_units": int(sub["unit_id"].nunique()),
                    "n_GSE123902": int((sub["dataset"] == "GSE123902").sum()),
                    "n_GSE131907": int((sub["dataset"] == "GSE131907").sum()),
                    "n_GSE205335": int((sub["dataset"] == "GSE205335").sum()),
                    "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                    "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                    "mean_AT2": float(sub["score_AT2"].mean()),
                    "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                    "mean_pseudotime": float(sub["pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
                }
            )
        pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
        pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
            tabdir / "paga_connectivities.tsv", sep="\t"
        )

    rows = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        rows.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "Sample": str(sub["Sample"].iloc[0]),
                "Sample_Origin": str(sub["Sample_Origin"].iloc[0]),
                "patient_id": str(sub["patient_id"].iloc[0]),
                "histology": str(sub["histology"].iloc[0]),
                "is_normal_tissue": str(sub["is_normal_tissue"].iloc[0])
                if "is_normal_tissue" in sub
                else str(sub["Sample_Origin"].iloc[0] in NORMAL_ORIGINS),
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
                "mean_pseudotime": float(sub["pseudotime"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
                "clock": clock,
                "graph_mode": graph_mode,
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

    # DONE CRITERION
    stacked = sample_df[
        [
            "unit_id",
            "dataset",
            "Sample",
            "Sample_Origin",
            "patient_id",
            "n_cells",
            "mean_CLDN4",
            "mean_pseudotime",
            "mean_dpt",
            "mean_AT2",
            "mean_barrier_keratin",
            "clock",
            "graph_mode",
        ]
    ].copy()
    stacked_path = tabdir / "stacked_sample_cldn4_pseudotime.tsv"
    stacked.to_csv(stacked_path, sep="\t", index=False)
    print(f"WROTE done-criterion table {stacked_path} n={len(stacked)}", flush=True)

    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs pseudotime", "mean_CLDN4", "mean_pseudotime"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs pseudotime (control)", "mean_SFTPC", "mean_pseudotime"),
        ("AT2 score vs pseudotime (control)", "mean_AT2", "mean_pseudotime"),
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

    def subset(frame, pred):
        return frame.loc[pred].copy()

    elig_123 = elig[elig["dataset"] == "GSE123902"]
    elig_131 = elig[elig["dataset"] == "GSE131907"]
    elig_205 = elig[elig["dataset"] == "GSE205335"]
    elig_t = elig[elig["Sample_Origin"].isin(["tLung", "PRIMARY"])]
    elig_n = elig[elig["Sample_Origin"].isin(list(NORMAL_ORIGINS))]
    elig_tumor = elig[~elig["Sample_Origin"].isin(list(NORMAL_ORIGINS))]
    elig_adcsq = elig_205[elig_205["histology"].isin(["ADC", "SQ"])]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE123902-only CLDN4 vs pseudotime", **sp(elig_123, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "GSE131907-only CLDN4 vs pseudotime", **sp(elig_131, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "GSE205335-only CLDN4 vs pseudotime", **sp(elig_205, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "GSE123902-only CLDN4 vs AT2", **sp(elig_123, "mean_CLDN4", "mean_AT2")},
        {"contrast": "GSE131907-only CLDN4 vs AT2", **sp(elig_131, "mean_CLDN4", "mean_AT2")},
        {"contrast": "GSE205335-only CLDN4 vs AT2", **sp(elig_205, "mean_CLDN4", "mean_AT2")},
        {"contrast": "GSE123902-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_123, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE131907-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_131, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE205335-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_205, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "tLung/PRIMARY-only CLDN4 vs pseudotime", **sp(elig_t, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "nLung/NORMAL-only CLDN4 vs pseudotime", **sp(elig_n, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "tumor-only (drop nLung/NORMAL) CLDN4 vs pseudotime", **sp(elig_tumor, "mean_CLDN4", "mean_pseudotime")},
        {"contrast": "GSE205335 ADC+SQ CLDN4 vs pseudotime", **sp(elig_adcsq, "mean_CLDN4", "mean_pseudotime")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "Sample_Origin": str(sub["Sample_Origin"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "malignant_high": float(hi["score_malignant_like"].mean()),
                "malignant_low": float(lo["score_malignant_like"].mean()),
                "pt_high": float(hi["pseudotime"].mean()),
                "pt_low": float(lo["pseudotime"].mean()),
                "club_high": float(hi["score_club"].mean()),
                "club_low": float(lo["score_club"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("malignant-like high vs low", "malignant_high", "malignant_low"),
        ("pseudotime high vs low", "pt_high", "pt_low"),
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
        or len(paired_df) >= 4
    )

    # ---- figures ----
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    ax = axes[0, 0]
    if graph_mode == "harmony_joint":
        sc.pl.paga(
            adata,
            color="expr_CLDN4",
            ax=ax,
            show=False,
            frameon=False,
            cmap="viridis",
            title="PAGA (Leiden) colored by mean CLDN4",
        )
    else:
        ax.set_axis_off()
        ax.set_title("PAGA is per-dataset (Harmony died)")
    ax = axes[0, 1]
    sc.pl.umap(adata, color="expr_CLDN4", ax=ax, show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    ax = axes[1, 0]
    sc.pl.umap(
        adata,
        color="pseudotime",
        ax=ax,
        show=False,
        frameon=False,
        cmap="magma",
        title=f"UMAP {clock} (root not CLDN4-high)",
    )
    ax = axes[1, 1]
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs pseudotime")
    for ds, col in COHORT_COLORS.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_pseudotime"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("sample-mean pseudotime")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_pt["rho"] is None else f"{c4_pt['rho']:.2f}"
    p_s = "NA" if c4_pt["p"] is None else f"{c4_pt['p']:.3g}"
    ax.set_title(f"stacked CLDN4 vs PT  n={c4_pt['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Triple-differ epithelium scored by CLDN4   n_cells={adata.n_obs}  n_units={sample_df.shape[0]}  mode={graph_mode}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for ds, col in COHORT_COLORS.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_pseudotime"], sub["mean_CLDN4"], s=52, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("sample-mean pseudotime")
    ax.set_ylabel("sample-mean CLDN4")
    ax.set_title(f"STACKED sample-level CLDN4 vs PT  {_fmt(c4_pt)}")
    ax.legend(fontsize=8, frameon=False)
    _save(fig, figdir / "fig_stacked_cldn4_pseudotime")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color=COHORT_COLORS)
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author/marker subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = sample_df.groupby("dataset").size().reindex(list(COHORTS)).fillna(0)
    axes[1].bar(
        unit_ct.index.astype(str),
        unit_ct.to_numpy(),
        color=[COHORT_COLORS[c] for c in unit_ct.index],
    )
    axes[1].set_ylabel("units")
    axes[1].set_title(f"Honest n units={sample_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("pt_low", "pt_high", "pseudotime"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for ds, col in COHORT_COLORS.items():
                sub = paired_df[paired_df["dataset"] == ds]
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{ds} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("AT2")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("pseudotime")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-unit CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("Sample_Origin", "fig_umap_origin", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        try:
            sc.pl.umap(adata, **kw)
        except Exception:
            ax.set_axis_off()
            ax.set_title(f"{fname} skipped")
        _save(fig, figdir / fname)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_AT2", "mean_CLDN4", c4_at2, "sample-mean AT2", "sample-mean CLDN4"),
        (
            axes[1],
            "mean_barrier_keratin",
            "mean_CLDN4",
            barrier_row,
            "sample-mean barrier/keratin (no CLDN4)",
            "sample-mean CLDN4",
        ),
    ):
        for ds, col in COHORT_COLORS.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: sample-level CLDN4 vs AT2 / barrier (CLDN4 excluded)", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_programs")

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs pseudotime")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = barrier_row
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    c4_t2 = next(r for r in primary if r["contrast"].startswith("CLDN4 vs TACSTD2"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))
    s123 = next(r for r in sensitivity if r["contrast"].startswith("GSE123902-only CLDN4 vs pseudotime"))
    s131 = next(r for r in sensitivity if r["contrast"].startswith("GSE131907-only CLDN4 vs pseudotime"))
    s205 = next(r for r in sensitivity if r["contrast"].startswith("GSE205335-only CLDN4 vs pseudotime"))

    parts = [
        f"Stacked sample-level CLDN4 vs AT2-rooted {clock}: {_fmt(c4_pt)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded from the score): {_fmt(c4_bar)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"CLDN4 vs TACSTD2 (comparator only): {_fmt(c4_t2)}.",
        f"GSE123902-only CLDN4 vs PT: {_fmt(s123)}.",
        f"GSE131907-only CLDN4 vs PT: {_fmt(s131)}.",
        f"GSE205335-only CLDN4 vs PT: {_fmt(s205)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"PAGA has {n_paga_components} component(s) at connectivity>0 among {n_leiden} Leiden vertices.",
        f"Graph mode={graph_mode}.",
        "The stacked pseudotime correlation mixes cohorts and normal vs tumor and is not a within-tumor progression test.",
        "Not a TACSTD2 redo. No both-high gate. GSE148071 not added. Not the 7-pool.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} units).** CLDN4 vs CLDN4-excluded "
        f"barrier/keratin: {_fmt(c4_bar)}. CLDN4 vs malignant-like: {_fmt(c4_mal)}. "
        f"Within-unit CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}. "
        f"**What does not automatically hold.** Stacked CLDN4 vs {clock}: {_fmt(c4_pt)}. "
        f"CLDN4 vs AT2: {_fmt(c4_at2)}."
    )

    n_normal_123 = int(
        (
            (adata.obs["dataset"] == "GSE123902")
            & adata.obs["Sample_Origin"].astype(str).eq("NORMAL")
        ).sum()
    )
    summary = {
        "accessions": list(COHORTS),
        "triple_that_differs": True,
        "gse148071_added": False,
        "seven_pool": False,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": clock,
        "graph_mode": graph_mode,
        "slingshot": sling,
        "slingshot_run": sling_run,
        "palantir": pal,
        "harmony": harmony,
        "per_dataset": per_dataset_info,
        "cap_per_unit": 350,
        "n_cells": int(adata.n_obs),
        "n_cells_gse123902": int((adata.obs["dataset"] == "GSE123902").sum()),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_cells_gse123902_normal": n_normal_123,
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_units": int(sample_df.shape[0]),
        "n_units_gse123902": int((sample_df["dataset"] == "GSE123902").sum()),
        "n_units_gse131907": int((sample_df["dataset"] == "GSE131907").sum()),
        "n_units_gse205335": int((sample_df["dataset"] == "GSE205335").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(n_paga_components),
        "n_leiden": int(n_leiden),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
        },
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "what_holds": what_holds,
        "done_table": str(stacked_path),
        "skipped": {
            "GSE148071": "explicitly not added",
            "7-pool": "explicitly not this folder",
            "GSE131907 PE": "unlabeled epithelium",
            "GSE205335 normal tissues": "Normal Lung / LN / Brain dropped",
            "Laughney H5": "36.5 GB skipped; marker epithelium",
            "slingshot": sling.get("reason"),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": list(COHORTS),
                "papers": {
                    "GSE123902": "Laughney et al. Nat Med 2020 PMID 31959988",
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "dual_high": False,
                "gse148071_added": False,
                "seven_pool": False,
                "slingshot": sling,
                "graph_mode": graph_mode,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "never CLDN4-high; prefer GSE131907 nLung author AT2",
                    "cap_per_unit": 350,
                },
            },
            indent=2,
        )
    )
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_units": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "graph_mode": graph_mode,
                "clock": clock,
                "extra": emit_extra,
                "finding": args.finding,
                "done_table": str(stacked_path),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
