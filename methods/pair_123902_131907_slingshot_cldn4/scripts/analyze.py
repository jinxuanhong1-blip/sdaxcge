#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE123902+GSE131907 epithelium, CLDN4 only.

ADDITIVE. PR #459 pair that differs. No dual-high. No GSE148071.
Root = GSE131907 nLung author AT2, never CLDN4-high.
Inferential unit = GSE123902 donor + GSE131907 sample.
Primary clock = Slingshot (Street 2018). PAGA is geometry. DPT is companion.
CLDN4 + barrier (CLDN4 excluded) + IFN are scored along pseudotime.
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

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


HERE = Path(__file__).resolve().parent
LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05
SLINGSHOT_RSCRIPT = HERE / "run_slingshot.R"


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


def require_slingshot() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit(
            "Rscript is required for REAL Slingshot. "
            "Install R and slingshot (see scripts/install_r_slingshot.sh). "
            "Do not stop at R-missing by falling back to DPT-only."
        )
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "library"))
    proc = subprocess.run(
        [rscript, "-e", 'cat(as.character(packageVersion("slingshot")))'],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )
    if proc.returncode != 0 or not (proc.stdout or "").strip():
        raise SystemExit(
            "slingshot R package is required. "
            f"probe stderr={(proc.stderr or '')[:400]} "
            "Run scripts/install_r_slingshot.sh. Do not fall back to DPT-only."
        )
    return {
        "available": True,
        "version": (proc.stdout or "").strip(),
        "rscript": rscript,
        "fallback": None,
    }


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


def pick_root_and_start_cluster(adata) -> tuple[int, str, dict]:
    """External arrow: GSE131907 nLung author AT2. Never root on CLDN4-high."""
    nlung = (adata.obs["dataset"].astype(str) == "GSE131907") & (
        adata.obs["Sample_Origin"].astype(str) == "nLung"
    )
    subtype = adata.obs["author_subtype"].astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    clust_mean_cldn = adata.obs.groupby("leiden", observed=True)["expr_CLDN4"].mean()
    cldn_rank = clust_mean_cldn.rank(ascending=False, method="min")
    high_clus = set(cldn_rank.index[cldn_rank <= max(1, int(np.ceil(len(cldn_rank) * 0.15)))].astype(str))
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "cldn4_high_clusters": sorted(high_clus),
        "cluster_mean_CLDN4": {str(k): float(v) for k, v in clust_mean_cldn.items()},
        "rule": None,
    }
    if int(at2.sum()) < 20:
        raise SystemExit("need ≥20 GSE131907 nLung author AT2 cells to root")
    at2_obs = adata.obs.loc[at2]
    counts = at2_obs["leiden"].astype(str).value_counts()
    # Prefer the AT2-rich cluster that is NOT CLDN4-high.
    start = None
    for cl, n in counts.items():
        if str(cl) not in high_clus and n >= 20:
            start = str(cl)
            break
    if start is None:
        # lowest mean CLDN4 among AT2-bearing clusters
        cand = [c for c in counts.index.astype(str) if counts[c] >= 10]
        if not cand:
            raise SystemExit("no AT2 cluster large enough to root")
        start = str(min(cand, key=lambda c: float(clust_mean_cldn.get(c, np.inf))))
        info["note"] = "all large AT2 clusters overlapped CLDN4-high rank; used lowest-CLDN4 AT2 cluster"
    if start in high_clus:
        raise SystemExit(
            f"refusing CLDN4-high start cluster {start} "
            f"(mean CLDN4={float(clust_mean_cldn.get(start, np.nan)):.3f})"
        )
    in_start = at2 & (adata.obs["leiden"].astype(str) == start)
    idx = np.flatnonzero(in_start.to_numpy())
    scores = adata.obs.loc[in_start, "score_AT2"].to_numpy()
    cldn = adata.obs.loc[in_start, "expr_CLDN4"].to_numpy()
    # median AT2, then break ties toward lower CLDN4
    order = np.lexsort((cldn, np.abs(scores - np.nanmedian(scores))))
    pick = int(idx[int(order[0])])
    info["rule"] = (
        "GSE131907 nLung author AT2 (median AT2, not CLDN4-high Leiden) "
        f"start_cluster={start}"
    )
    info["start_cluster"] = start
    info["start_cluster_mean_CLDN4"] = float(clust_mean_cldn.get(start, np.nan))
    info["start_cluster_n_AT2"] = int(in_start.sum())
    return pick, start, info


def run_slingshot(adata, start_cluster: str, work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    if "X_pca_harmony" in adata.obsm:
        rd = np.asarray(adata.obsm["X_pca_harmony"][:, :N_PCS], dtype=float)
    else:
        rd = np.asarray(adata.obsm["X_pca"][:, :N_PCS], dtype=float)
    cells = adata.obs_names.astype(str).to_numpy()
    embed = pd.DataFrame(rd, index=cells, columns=[f"PC{i+1}" for i in range(rd.shape[1])])
    embed_path = work / "embedding.tsv"
    cl_path = work / "clusters.tsv"
    embed.to_csv(embed_path, sep="\t")
    pd.DataFrame({"cell": cells, "cluster": adata.obs["leiden"].astype(str).to_numpy()}).to_csv(
        cl_path, sep="\t", index=False
    )
    rscript = shutil.which("Rscript")
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "library"))
    cmd = [
        rscript,
        str(SLINGSHOT_RSCRIPT),
        "--embed",
        str(embed_path),
        "--clusters",
        str(cl_path),
        "--start",
        str(start_cluster),
        "--outdir",
        str(work),
    ]
    print("RUN", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
    (work / "slingshot_stdout.txt").write_text(proc.stdout or "")
    (work / "slingshot_stderr.txt").write_text(proc.stderr or "")
    if proc.returncode != 0:
        raise SystemExit(
            "REAL Slingshot failed.\n"
            f"stdout:\n{(proc.stdout or '')[-2000:]}\n"
            f"stderr:\n{(proc.stderr or '')[-2000:]}"
        )
    pt = pd.read_csv(work / "slingshot_pseudotime.tsv", sep="\t")
    w = pd.read_csv(work / "slingshot_weights.tsv", sep="\t")
    lin = pd.read_csv(work / "slingshot_lineages.tsv", sep="\t")
    pt = pt.set_index("cell").reindex(cells)
    w = w.set_index("cell").reindex(cells)
    pt_cols = [c for c in pt.columns if c != "cell"]
    w_cols = [c for c in w.columns if c != "cell"]
    if not pt_cols:
        raise SystemExit("Slingshot returned no lineage columns")
    # Primary lineage = Lineage1 / first column (starts at the AT2 cluster).
    primary_col = pt_cols[0]
    for cand in ("Lineage1", "curve1"):
        if cand in pt_cols:
            primary_col = cand
            break
    adata.obs["sling_lineage_primary"] = primary_col
    adata.obs["sling_pseudotime"] = pd.to_numeric(pt[primary_col], errors="coerce").to_numpy()
    # Weighted mean across lineages for cells on any curve.
    w_mat = w[w_cols].to_numpy(dtype=float) if w_cols else np.ones((len(cells), 1))
    pt_mat = pt[pt_cols].to_numpy(dtype=float)
    w_mat = np.where(np.isfinite(pt_mat), w_mat, 0.0)
    wsum = w_mat.sum(axis=1)
    weighted = np.divide((np.nan_to_num(pt_mat) * w_mat).sum(axis=1), wsum, out=np.full(len(cells), np.nan), where=wsum > 0)
    adata.obs["sling_pseudotime_weighted"] = weighted
    n_on = int(np.isfinite(adata.obs["sling_pseudotime"]).sum())
    if n_on < 100:
        raise SystemExit(f"too few cells on primary Slingshot lineage: {n_on}")
    return {
        "n_lineages": int(len(pt_cols)),
        "primary_lineage": primary_col,
        "pt_columns": pt_cols,
        "n_cells_on_primary": n_on,
        "lineages": lin.to_dict(orient="records"),
        "stdout_tail": (proc.stdout or "")[-500:],
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
    lines = [
        "# Finding — pair GSE123902+GSE131907 epithelium, CLDN4-only REAL Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** Pair that **differs** in PR #459 "
        "(malignant CLDN4 %pos GSE123902+GSE131907, n=34, ρ=−0.575 vs T/NK). "
        "That T/NK rho is **not re-audited**. This folder asks a different question: "
        "where do **CLDN4**, a CLDN4-excluded **barrier/keratin** score, and a compact **IFN** "
        "score sit on a real Slingshot lineage. "
        "Does **not** redo GSE131907-only PAGA (PR #325) or the winning-pair "
        "GSE131907+GSE205335 Slingshot/DPT (PR #449). "
        "GSE148071 is not added. No TACSTD2∩CLDN4 dual-high gate.",
        "",
        f"Primary clock: **Slingshot {sling.get('version')}** "
        f"({s['slingshot_run'].get('n_lineages')} lineage(s); primary "
        f"`{s['slingshot_run'].get('primary_lineage')}`). "
        "PAGA is geometry only. DPT is a companion ordering, not the claim clock. "
        "Inferential unit = **GSE123902 donor + GSE131907 sample**. "
        "Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        "Root is GSE131907 nLung author AT2, never CLDN4-high.",
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
        f"(GSE131907 {s['n_cells_gse131907']}, GSE123902 {s['n_cells_gse123902']}).",
        f"- Units (GSE131907 Sample + GSE123902 donor): **n_units = {s['n_units']}** "
        f"(GSE131907 {s['n_units_gse131907']}, GSE123902 {s['n_units_gse123902']}).",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: {s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- Author / marker subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Cells on primary Slingshot lineage: **{s['slingshot_run'].get('n_cells_on_primary')}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE148071 not used. GSE123902 NORMAL samples dropped. "
        "GSE131907 PE unlabeled epithelium dropped. Author 36.5 GB H5 skipped.",
        f"- Slingshot: available=True; version={sling.get('version')}.",
        f"- Start cluster: {s['root'].get('start_cluster')} "
        f"(mean CLDN4={s['root'].get('start_cluster_mean_CLDN4')}); "
        f"CLDN4-high clusters refused as root: {s['root'].get('cldn4_high_clusters')}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Slingshot / DPT root: {s['root']['rule']} "
        f"(root cell index {s['root']['index']}, unit {s['root'].get('root_unit')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN: compact ISG / IFNG-response core (STAT1, IRF1/7/9, ISG15, IFIT1-3, MX1/2, OAS*, CXCL9/10/11, IDO1, GBP*, IFI*, RSAD2, B2M, TAP1, PSMB8/9). **CLDN4 out**.",
        "",
        "## Lineages (done criterion)",
        "",
        f"Primary lineage `{s['slingshot_run'].get('primary_lineage')}` starts at Leiden "
        f"{s['root'].get('start_cluster')} (nLung AT2, not CLDN4-high). "
        f"{s['slingshot_run'].get('n_lineages')} Slingshot lineage(s).",
        "",
        "| lineage | start | end | n_clusters | n_cells | path |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]
    for rec in s["slingshot_run"].get("lineages", []):
        lines.append(
            f"| {rec.get('lineage')} | {rec.get('start_cluster')} | {rec.get('end_cluster')} | "
            f"{rec.get('n_clusters')} | {rec.get('n_cells_on_lineage')} | `{rec.get('path')}` |"
        )
    lines += [
        "",
        "Machine table: `results/tables/slingshot_lineages.tsv`.",
        "",
        "## Primary (sample-level Spearman, BH inside this list)",
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
        "- Do not write n_cells as the inferential n. GSE123902 is **donor**; GSE131907 is **sample**.",
        "- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (GSE131907+GSE205335 DPT).",
        "- The pooled Slingshot Spearman mixes cohorts and nLung vs tumor and is **not** a within-tumor progression test.",
        "- GSE123902 epithelium is a marker gate (EPCAM/KRT+, PTPRC−), not the skipped 36.5 GB author H5.",
        "- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.",
        "- PR #459 T/NK ρ is given and was not re-audited.",
        "- No TACSTD2∩CLDN4 both-high gate. GSE148071 not added.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot is an ordering of Harmony space, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/slingshot_lineages.tsv` — **done criterion (lineage)**",
        "- `results/tables/sample_level_spearman.tsv` — **done criterion (sample-level)**",
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
        "bash methods/pair_123902_131907_slingshot_cldn4/scripts/install_r_slingshot.sh",
        "pip install -r methods/pair_123902_131907_slingshot_cldn4/requirements.txt",
        "python3 methods/pair_123902_131907_slingshot_cldn4/scripts/download.py \\",
        "  --out /tmp/pair_123902_131907",
        "python3 methods/pair_123902_131907_slingshot_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/pair_123902_131907 \\",
        "  --out /tmp/pair_123902_131907/epithelium.h5ad",
        "python3 methods/pair_123902_131907_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/pair_123902_131907/epithelium.h5ad \\",
        "  --outdir methods/pair_123902_131907_slingshot_cldn4/results \\",
        "  --finding methods/pair_123902_131907_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/pair_123902_131907_slingshot_cldn4/results")
    p.add_argument("--finding", default="methods/pair_123902_131907_slingshot_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling = require_slingshot()
    print(json.dumps({"slingshot": sling}, indent=2), flush=True)

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

    root_i, start_cluster, root_info = pick_root_and_start_cluster(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_unit"] = str(adata.obs.iloc[root_i].get("unit_id", ""))
    root_info["root_origin"] = str(adata.obs.iloc[root_i].get("Sample_Origin", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("author_subtype", ""))
    root_info["root_dataset"] = str(adata.obs.iloc[root_i].get("dataset", ""))
    root_info["root_CLDN4"] = float(adata.obs.iloc[root_i]["expr_CLDN4"])

    sling_run = run_slingshot(adata, start_cluster, outdir / "slingshot_work")
    # copy lineage table to results/tables (done criterion)
    src_lin = outdir / "slingshot_work" / "slingshot_lineages.tsv"
    if src_lin.is_file():
        shutil.copy(src_lin, tabdir / "slingshot_lineages.tsv")
        shutil.copy(outdir / "slingshot_work" / "slingshot_pseudotime.tsv", tabdir / "slingshot_pseudotime.tsv")
        shutil.copy(outdir / "slingshot_work" / "slingshot_weights.tsv", tabdir / "slingshot_weights.tsv")

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
                "n_GSE131907": int((sub["dataset"] == "GSE131907").sum()),
                "n_GSE123902": int((sub["dataset"] == "GSE123902").sum()),
                "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_sling": float(sub["sling_pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
                "is_start": cl == start_cluster,
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    sling_pt = adata.obs["sling_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["sling_pseudotime"] = sling_pt

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
    if "expr_SFTPC" in adata.obs:
        sft_map = {
            unit: float(sub["expr_SFTPC"].mean())
            for unit, sub in adata.obs.groupby("unit_id", observed=True)
        }
        sample_df["mean_SFTPC"] = sample_df["unit_id"].map(sft_map)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling"),
        ("barrier/keratin (no CLDN4) vs Slingshot PT", "mean_barrier_keratin", "mean_sling"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_sling"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling"),
        ("AT2 score vs Slingshot PT (control)", "mean_AT2", "mean_sling"),
        ("CLDN4 vs DPT (companion)", "mean_CLDN4", "mean_dpt"),
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

    elig_131 = elig[elig["dataset"] == "GSE131907"]
    elig_123 = elig[elig["dataset"] == "GSE123902"]
    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    elig_tumor = elig[~elig["Sample_Origin"].isin(["nLung"])]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE131907-only CLDN4 vs Slingshot PT", **sp(elig_131, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE123902-only CLDN4 vs Slingshot PT", **sp(elig_123, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE131907-only barrier vs Slingshot PT", **sp(elig_131, "mean_barrier_keratin", "mean_sling")},
        {"contrast": "GSE123902-only barrier vs Slingshot PT", **sp(elig_123, "mean_barrier_keratin", "mean_sling")},
        {"contrast": "GSE131907-only IFN vs Slingshot PT", **sp(elig_131, "mean_IFN", "mean_sling")},
        {"contrast": "GSE123902-only IFN vs Slingshot PT", **sp(elig_123, "mean_IFN", "mean_sling")},
        {"contrast": "GSE131907-only CLDN4 vs IFN", **sp(elig_131, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE123902-only CLDN4 vs IFN", **sp(elig_123, "mean_CLDN4", "mean_IFN")},
        {"contrast": "tLung-only CLDN4 vs Slingshot PT", **sp(elig_t, "mean_CLDN4", "mean_sling")},
        {"contrast": "nLung-only CLDN4 vs Slingshot PT", **sp(elig_n, "mean_CLDN4", "mean_sling")},
        {"contrast": "tumor-only (drop nLung) CLDN4 vs Slingshot PT", **sp(elig_tumor, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE131907-only CLDN4 vs DPT (companion)", **sp(elig_131, "mean_CLDN4", "mean_dpt")},
        {"contrast": "GSE123902-only CLDN4 vs DPT (companion)", **sp(elig_123, "mean_CLDN4", "mean_dpt")},
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
        ("IFN high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("malignant-like high vs low", "malignant_high", "malignant_low"),
        ("Slingshot PT high vs low", "sling_high", "sling_low"),
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

    colors = {"GSE131907": "#2a6f97", "GSE123902": "#6a994e"}

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
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
    sc.pl.umap(
        adata,
        color="sling_pseudotime",
        ax=ax,
        show=False,
        frameon=False,
        cmap="magma",
        title="UMAP Slingshot PT (AT2-rooted)",
    )
    ax = axes[1, 1]
    c4_sling = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_sling"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("sample-mean Slingshot PT")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_sling["rho"] is None else f"{c4_sling['rho']:.2f}"
    p_s = "NA" if c4_sling["p"] is None else f"{c4_sling['p']:.3g}"
    ax.set_title(f"CLDN4 vs Slingshot  n={c4_sling['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Pair GSE123902+GSE131907 epithelium scored by CLDN4   n_cells={adata.n_obs}  n_units={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # Extra: CLDN4 + barrier + IFN along Slingshot PT (sample-level)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.9))
    panels = (
        ("mean_CLDN4", "CLDN4 vs Slingshot PT", "sample-mean CLDN4"),
        ("mean_barrier_keratin", "barrier/keratin (no CLDN4) vs Slingshot PT", "sample-mean barrier (no CLDN4)"),
        ("mean_IFN", "IFN vs Slingshot PT", "sample-mean IFN"),
    )
    for ax, (y, contrast, ylab) in zip(axes, panels):
        row = next(r for r in primary if r["contrast"] == contrast)
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub["mean_sling"], sub[y], s=44, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel("sample-mean Slingshot PT")
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: CLDN4 + barrier + IFN along Slingshot PT (unit means)", fontsize=11)
    _save(fig, figdir / "fig_extra_along_pseudotime")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color={"GSE131907": "#2a6f97", "GSE123902": "#6a994e"})
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author / marker subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = sample_df.groupby("dataset").size()
    axes[1].bar(
        unit_ct.index.astype(str),
        unit_ct.to_numpy(),
        color=[colors.get(str(i), "#444") for i in unit_ct.index],
    )
    axes[1].set_ylabel("units")
    axes[1].set_title(
        f"Honest n units={sample_df.shape[0]} (eligible {len(elig)}; donor+sample)"
    )
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN"),
            ("sling_low", "sling_high", "Slingshot PT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for ds, col in colors.items():
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
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("IFN")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("Slingshot")), keys=("W", "p")))
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
        ("score_IFN", "fig_umap_IFN", "viridis"),
        ("sling_pseudotime", "fig_umap_slingshot", "magma"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_barrier_keratin", "mean_CLDN4", barrier_row, "sample-mean barrier/keratin (no CLDN4)", "sample-mean CLDN4"),
        (axes[1], "mean_IFN", "mean_CLDN4", c4_ifn, "sample-mean IFN", "sample-mean CLDN4"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: sample-level CLDN4 vs barrier / IFN", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_programs")

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_sling = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    bar_sling = next(r for r in primary if r["contrast"].startswith("barrier/keratin"))
    ifn_sling = next(r for r in primary if r["contrast"] == "IFN vs Slingshot PT")
    c4_bar = barrier_row
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    s131 = next(r for r in sensitivity if r["contrast"] == "GSE131907-only CLDN4 vs Slingshot PT")
    s123 = next(r for r in sensitivity if r["contrast"] == "GSE123902-only CLDN4 vs Slingshot PT")

    parts = [
        f"Sample-level CLDN4 vs AT2-rooted Slingshot PT: {_fmt(c4_sling)}.",
        f"Barrier/keratin (no CLDN4) vs Slingshot PT: {_fmt(bar_sling)}.",
        f"IFN vs Slingshot PT: {_fmt(ifn_sling)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"GSE131907-only CLDN4 vs Slingshot: {_fmt(s131)}.",
        f"GSE123902-only CLDN4 vs Slingshot: {_fmt(s123)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        f"Slingshot produced {sling_run['n_lineages']} lineage(s); primary={sling_run['primary_lineage']}.",
        "The pooled Slingshot correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test.",
        "REAL Slingshot was run. Not a TACSTD2 redo. No both-high gate. GSE148071 not added.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_sling['n']} units).** "
        f"CLDN4 vs Slingshot PT is {_fmt(c4_sling)}. "
        f"Barrier vs Slingshot PT is {_fmt(bar_sling)}. "
        f"IFN vs Slingshot PT is {_fmt(ifn_sling)}. "
        f"CLDN4 tracks barrier/keratin ({_fmt(c4_bar)}) and IFN ({_fmt(c4_ifn)}). "
        f"Within-unit CLDN4-high vs low barrier: {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN: {_fmt(pair_ifn, keys=('W', 'p'))}."
    )

    summary = {
        "accessions": ["GSE123902", "GSE131907"],
        "pr459_pair_that_differs": True,
        "gse148071_added": False,
        "not_a_pr325_redo": True,
        "not_a_pr449_redo": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": "slingshot",
        "slingshot": sling,
        "slingshot_run": sling_run,
        "harmony": harmony,
        "cap_per_unit": 350,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse123902": int((adata.obs["dataset"] == "GSE123902").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_units": int(sample_df.shape[0]),
        "n_units_gse131907": int((sample_df["dataset"] == "GSE131907").sum()),
        "n_units_gse123902": int((sample_df["dataset"] == "GSE123902").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
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
        "skipped": {
            "GSE148071": "explicitly not added",
            "GSE123902 NORMAL": "dropped (PR #459)",
            "GSE123902 author H5": "36.5 GB skipped",
            "GSE131907 PE": "unlabeled epithelium",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE123902", "GSE131907"],
                "papers": {
                    "GSE123902": "Laughney et al. Nat Med 2020 (human arm of GSE123904)",
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_excludes_CLDN4": True,
                "slingshot": sling,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "GSE131907 nLung author AT2, never CLDN4-high",
                    "unit": "GSE123902 donor + GSE131907 sample",
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
                "n_lineages": sling_run["n_lineages"],
                "extra": emit_extra,
                "finding": args.finding,
                "lineage_table": str(tabdir / "slingshot_lineages.tsv"),
                "sample_table": str(tabdir / "sample_level_spearman.tsv"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
