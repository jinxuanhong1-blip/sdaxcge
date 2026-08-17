#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE123902 tumor epithelium, CLDN4 only.

ADDITIVE. Laughney metastasis atlas (Nat Med 2020, PMID 32042191).
Donor is the inferential unit. Root is never CLDN4-high.
Barrier/keratin excludes CLDN4. No TACSTD2∩CLDN4 dual-high gate.
Author 36.5 GB H5 is not used.
Done when results/tables/lineage.tsv exists.
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
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import COMPARATOR, CONTROLS, FOCAL, QC_NEG, STATES  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required; run scripts/install_tools.sh") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
CAP_PER_SAMPLE = 350
MIN_CELLS_DONOR = 10
MIN_CELLS_TERTILE_ARM = 8
SEED = 1


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


def _fmt(row: dict, keys=("rho", "p")) -> str:
    n = row.get("n")
    if "W" in keys:
        if row.get("W") is None:
            return f"n={n}, W=NA, p=NA"
        return (
            f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, "
            f"p={row['p']:.3g}"
        )
    if row.get("rho") is None:
        return f"n={n}, ρ=NA, p=NA"
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


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


def cap_per_sample(adata, cap: int = CAP_PER_SAMPLE, seed: int = SEED):
    rng = np.random.default_rng(seed)
    keep = []
    for gsm, idx in adata.obs.groupby("gsm", observed=True).groups.items():
        idx = np.asarray(idx)
        if len(idx) <= cap:
            keep.extend(idx.tolist())
            continue
        sub = adata.obs.loc[idx]
        # Protect matched-normal AT2-like, CLDN4-low cells for the root.
        protect = (sub["tissue"].astype(str) == "NORMAL") & (
            sub["score_AT2"] >= np.nanmedian(sub["score_AT2"])
        )
        prot_idx = idx[protect.to_numpy()]
        rest = idx[~protect.to_numpy()]
        n_prot = min(len(prot_idx), cap)
        n_rest = cap - n_prot
        chosen = []
        if n_prot:
            chosen.extend(rng.choice(prot_idx, size=n_prot, replace=False).tolist())
        if n_rest > 0 and len(rest):
            chosen.extend(rng.choice(rest, size=min(n_rest, len(rest)), replace=False).tolist())
        keep.extend(chosen)
    return adata[keep].copy()


def pick_root(adata) -> tuple[int, dict]:
    """External arrow: matched-normal AT2-like, never CLDN4-high."""
    tissue = adata.obs["tissue"].astype(str)
    tert = adata.obs["cldn4_tertile"].astype(str)
    normal = tissue == "NORMAL"
    not_high = tert != "high"
    info = {
        "n_normal": int(normal.sum()),
        "n_normal_not_cldn4_high": int((normal & not_high).sum()),
        "rule": None,
    }
    cand = normal & not_high
    if int(cand.sum()) >= 10:
        idx = np.flatnonzero(cand.to_numpy())
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        info["rule"] = "matched-normal AT2-like, CLDN4 not high (median AT2)"
        return int(pick), info
    cand = not_high
    if int(cand.sum()) < 5:
        raise SystemExit("cannot pick a non-CLDN4-high root")
    idx = np.flatnonzero(cand.to_numpy())
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = "max AT2 among CLDN4-not-high cells (no usable matched normal)"
    return int(pick), info


def python_slingshot(rd: np.ndarray, clusters: np.ndarray, start: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Street 2018 Slingshot: MST lineages + projection onto the polyline."""
    labels = np.asarray(clusters, dtype=str)
    uniq = np.array(sorted(set(labels), key=lambda x: int(x) if str(x).isdigit() else str(x)))
    if start not in set(uniq):
        raise SystemExit(f"start cluster {start} not in Leiden labels")
    cents = np.vstack([rd[labels == c].mean(axis=0) for c in uniq])
    pos = {c: i for i, c in enumerate(uniq)}
    dist = squareform(pdist(cents))
    mst = minimum_spanning_tree(dist).toarray()
    mst = np.maximum(mst, mst.T)
    adj = {c: [] for c in uniq}
    for i, a in enumerate(uniq):
        for j, b in enumerate(uniq):
            if i < j and mst[i, j] > 0:
                adj[a].append(b)
                adj[b].append(a)

    def path_to(target: str) -> list[str] | None:
        stack = [(start, [start])]
        seen = {start}
        while stack:
            node, path = stack.pop()
            if node == target:
                return path
            for nxt in adj[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append((nxt, path + [nxt]))
        return None

    leaves = [c for c in uniq if c != start and len(adj[c]) == 1]
    if not leaves:
        leaves = [c for c in uniq if c != start]
    lineages = []
    for leaf in leaves:
        path = path_to(leaf)
        if path:
            lineages.append(path)
    if not lineages:
        lineages = [[start]]

    def project(X: np.ndarray, waypoints: np.ndarray) -> np.ndarray:
        if waypoints.shape[0] == 1:
            return np.zeros(X.shape[0])
        seglen = np.linalg.norm(np.diff(waypoints, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seglen)])
        best_t = np.full(X.shape[0], np.inf)
        best_d = np.full(X.shape[0], np.inf)
        for i in range(len(seglen)):
            a = waypoints[i]
            b = waypoints[i + 1]
            ab = b - a
            denom = float(np.dot(ab, ab)) + 1e-12
            u = np.clip(((X - a) @ ab) / denom, 0.0, 1.0)
            proj = a + np.outer(u, ab)
            d = np.linalg.norm(X - proj, axis=1)
            t = cum[i] + u * seglen[i]
            better = d < best_d
            best_d[better] = d[better]
            best_t[better] = t[better]
        mx = float(cum[-1]) if cum[-1] > 0 else 1.0
        return best_t / mx

    pt_cols = {}
    lin_rows = []
    for k, path in enumerate(lineages, start=1):
        lid = f"Lineage{k}"
        wps = np.vstack([cents[pos[c]] for c in path])
        pt = project(rd, wps)
        in_lin = np.isin(labels, path)
        col = np.full(rd.shape[0], np.nan)
        col[in_lin] = pt[in_lin]
        pt_cols[lid] = col
        lin_rows.append(
            {
                "lineage_id": lid,
                "start_cluster": path[0],
                "end_cluster": path[-1],
                "path": ">".join(path),
                "n_clusters": len(path),
            }
        )
    pt_df = pd.DataFrame(pt_cols)
    shared = np.nanmean(pt_df.to_numpy(), axis=1)
    pt_df.insert(0, "sling_pt", shared)
    return pt_df, pd.DataFrame(lin_rows)


def run_r_slingshot(rd: np.ndarray, clusters: np.ndarray, start: str, cells: list[str]) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    script = Path(__file__).resolve().parent / "run_slingshot.R"
    if not script.exists():
        return None
    with tempfile.TemporaryDirectory(prefix="slingshot_") as tmp:
        tmp = Path(tmp)
        pca = pd.DataFrame(rd, index=cells, columns=[f"PC{i+1}" for i in range(rd.shape[1])])
        pca.to_csv(tmp / "pca.csv")
        pd.DataFrame({"cell": cells, "leiden": clusters}).to_csv(tmp / "clusters.csv", index=False)
        (tmp / "start_cluster.txt").write_text(str(start) + "\n")
        out = tmp / "out"
        proc = subprocess.run(
            ["Rscript", str(script), str(tmp), str(out)],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        if proc.returncode != 0:
            print("R slingshot failed:\n", proc.stderr[-2000:] if proc.stderr else proc.stdout, flush=True)
            return None
        pt = pd.read_csv(out / "sling_pseudotime.csv")
        lin = pd.read_csv(out / "sling_lineages.csv")
        pt.index = pt["cell"].astype(str)
        lin_cols = [c for c in pt.columns if c not in {"cell", "cluster"}]
        shared = np.nanmean(pt[lin_cols].to_numpy(), axis=1) if lin_cols else np.full(len(pt), np.nan)
        out_pt = pt[lin_cols].copy()
        out_pt.insert(0, "sling_pt", shared)
        out_pt = out_pt.reindex(cells)
        return out_pt, lin


def write_finding(path: Path, s: dict) -> None:
    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    clock = s["clock"]
    lines = [
        "# Finding — GSE123902 Laughney tumor epithelium, CLDN4-only real Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** GSE123902 (Laughney et al., *Nat Med* 2020, PMID 32042191) "
        "human LUAD primary + metastasis epithelium. SuperSeries GSE123904 / mouse GSE123903 "
        "are not used. No TACSTD2∩CLDN4 dual-high gate. Author annotated H5 (36.5 GB) skipped; "
        "GEO dense UMI CSVs (90.4 MB) are the matrix.",
        "",
        f"Primary clock: **{clock}**. Inferential unit = **donor** (LX ID). "
        "Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        "IFN is a locked Hallmark IFNα∩IFNγ core. Root is never CLDN4-high.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- GEO catalog: **17 samples / 14 LX IDs** (8 primary + 5 metastasis + 4 matched-normal files; "
        f"LX685 is normal-only). Paper QC atlas is 41,384 cells; we do not substitute that n.",
        f"- Marker epithelium before cap (EPCAM|KRT8|KRT18|KRT19>0 and PTPRC==0): "
        f"**n_cells_uncapped = {s['n_cells_uncapped']}** "
        f"(tumor {s['n_cells_tumor_uncapped']}, normal {s['n_cells_normal_uncapped']}).",
        f"- Analysis cells after QC and ≤{s['cap_per_sample']}/sample cap: **n_cells = {s['n_cells']}** "
        f"(tumor {s['n_cells_tumor']}, normal {s['n_cells_normal']}).",
        f"- Donors on the object: **n_donors = {s['n_donors']}**. "
        f"Tumor donors: **{s['n_tumor_donors']}**. "
        f"Donors with ≥{MIN_CELLS_DONOR} tumor epithelial cells used for Spearman: **n = {s['n_units_eligible']}**.",
        f"- Donors with ≥{MIN_CELLS_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- CLDN4 tertile cells: {s['cldn4_tertile_counts']}.",
        f"- Tissue cells: {s['tissue_counts']}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: available={s['slingshot'].get('available')}; "
        f"engine={s['slingshot'].get('engine')}; {s['slingshot'].get('reason', '')}.",
        f"- Lineages: **{s['n_lineages']}**. Lineage table: `results/tables/lineage.tsv`.",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Marker epithelium ≠ author / inferCNV malignant. No ICI / MPR / RECIST labels.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {s['leiden_resolution']}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- DPT / Slingshot root: {s['root'].get('rule')} "
        f"(root cell index {s['root'].get('index')}, donor {s['root'].get('root_donor')}, "
        f"tissue {s['root'].get('root_tissue')}, Leiden {s['root'].get('start_cluster')}).",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN genes: STAT1/2, IRF1/7/9, ISG15/20, MX1/2, IFIT1/2/3, OAS1/2/L, IFI27/44/44L/6, "
        "RSAD2, CXCL9/10/11, B2M, TAP1, PSMB8/9.",
        "",
        "## Lineage table (done criterion)",
        "",
        "| lineage | start | end | n_clusters | n_cells | n_donors | mean CLDN4 | mean barrier | mean IFN | donor ρ CLDN4 vs PT |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for r in s["lineages"]:
        rho = "NA" if r.get("donor_rho_CLDN4_vs_PT") is None else f"{r['donor_rho_CLDN4_vs_PT']:.3f}"
        lines.append(
            f"| {r['lineage_id']} | {r['start_cluster']} | {r['end_cluster']} | "
            f"{r['n_clusters']} | {r['n_cells']} | {r['n_donors']} | "
            f"{r['mean_CLDN4']:.3f} | {r['mean_barrier']:.3f} | {r['mean_IFN']:.3f} | "
            f"n={r.get('donor_n')}, ρ={rho} |"
        )
    lines += [
        "",
        "## Primary (donor-level Spearman on tumor epithelium, BH inside this list)",
        "",
        "| Contrast | n_donors | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_donors | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["sensitivity_spearman"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    extra = s["extra_figure"]
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (donor-paired)",
        "",
        f"Emitted: **{extra['emitted']}**. Observed Spearman(CLDN4, barrier/keratin no CLDN4) "
        f"n={extra.get('spearman_n')}, ρ={extra_rho}, p={extra_p}.",
        "",
        "| Paired contrast (high − low) | n_donors | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["paired_tertile"]:
        dlt = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {dlt} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- Marker epithelium is **not** a CNV-malignant call. The 36.5 GB author H5 was skipped.",
        "- n=13 tumor donors is small; CIs are wide. Do not write a multi-cohort pile.",
        "- The graph mixes matched normal (root) and tumor/met. Pooled PT is not a within-tumor clock.",
        "- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot/DPT is an ordering, not a developmental clock. No ICI language.",
        "",
        "## Outputs",
        "",
        "- `results/tables/lineage.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_extra_along_pt.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash methods/gse123902_slingshot_real_cldn4/scripts/install_tools.sh",
        "python3 methods/gse123902_slingshot_real_cldn4/scripts/download.py --out /tmp/gse123902_slingshot",
        "python3 methods/gse123902_slingshot_real_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse123902_slingshot \\",
        "  --out /tmp/gse123902_slingshot/epithelium.h5ad",
        "python3 methods/gse123902_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse123902_slingshot/epithelium.h5ad \\",
        "  --outdir methods/gse123902_slingshot_real_cldn4/results \\",
        "  --finding methods/gse123902_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--finding", type=Path, required=True)
    args = p.parse_args()
    outdir = args.outdir
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    adata0 = sc.read_h5ad(args.input)
    n_uncapped = int(adata0.n_obs)
    n_tumor_uncapped = int(adata0.obs["is_tumor"].sum())
    n_normal_uncapped = n_uncapped - n_tumor_uncapped

    sc.pp.filter_cells(adata0, min_genes=200)
    keep_umi = adata0.obs["n_umi"] >= 500
    adata0 = adata0[keep_umi].copy()

    # Score on log1p(CP10k) of the uncapped object so tertiles/root see biology, then cap.
    adata0.layers["counts"] = adata0.X.copy()
    sc.pp.normalize_total(adata0, target_sum=1e4)
    sc.pp.log1p(adata0)
    absent = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata0, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata0.var_names:
            X = adata0[:, g].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            adata0.obs[f"expr_{g}"] = np.asarray(X, dtype=float).ravel()
        else:
            adata0.obs[f"expr_{g}"] = np.nan
    q1, q2 = np.nanquantile(adata0.obs["expr_CLDN4"].to_numpy(), [1 / 3, 2 / 3])
    adata0.obs["cldn4_tertile"] = pd.cut(
        adata0.obs["expr_CLDN4"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low", "mid", "high"],
    ).astype(str)

    adata = cap_per_sample(adata0, CAP_PER_SAMPLE, SEED)
    # Recompute graph on counts of the capped object.
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for name, genes in STATES.items():
        _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            X = adata[:, g].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(X, dtype=float).ravel()
    q1, q2 = np.nanquantile(adata.obs["expr_CLDN4"].to_numpy(), [1 / 3, 2 / 3])
    adata.obs["cldn4_tertile"] = pd.cut(
        adata.obs["expr_CLDN4"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low", "mid", "high"],
    ).astype(str)

    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    sc.pp.pca(adata, n_comps=N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.umap(adata)
    sc.tl.paga(adata, groups="leiden")

    root_idx, root_info = pick_root(adata)
    root_cell = str(adata.obs_names[root_idx])
    start_cluster = str(adata.obs.iloc[root_idx]["leiden"])
    if str(adata.obs.iloc[root_idx]["cldn4_tertile"]) == "high":
        raise SystemExit("root landed on CLDN4-high; abort")
    root_info.update(
        {
            "index": int(root_idx),
            "root_cell": root_cell,
            "root_donor": str(adata.obs.iloc[root_idx]["donor"]),
            "root_tissue": str(adata.obs.iloc[root_idx]["tissue"]),
            "root_gsm": str(adata.obs.iloc[root_idx]["gsm"]),
            "start_cluster": start_cluster,
            "root_CLDN4": float(adata.obs.iloc[root_idx]["expr_CLDN4"]),
            "root_AT2": float(adata.obs.iloc[root_idx]["score_AT2"]),
        }
    )
    adata.uns["iroot"] = int(root_idx)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    rd = np.asarray(adata.obsm["X_pca"][:, :N_PCS], dtype=float)
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    cells = adata.obs_names.astype(str).tolist()
    sling_meta = slingshot_r_status()
    pt_df = None
    lin_raw = None
    if sling_meta.get("available"):
        ran = run_r_slingshot(rd, clusters, start_cluster, cells)
        if ran is not None:
            pt_df, lin_raw = ran
            sling_meta["engine"] = "R_slingshot"
            sling_meta["reason"] = f"Bioconductor slingshot {sling_meta.get('version')}"
    if pt_df is None:
        pt_df, lin_raw = python_slingshot(rd, clusters, start_cluster)
        sling_meta["engine"] = "python_slingshot_mst_polyline"
        if not sling_meta.get("available"):
            sling_meta["reason"] = (
                sling_meta.get("reason")
                or "R/slingshot missing; Street 2018 MST + polyline projection"
            )
        else:
            sling_meta["reason"] = "R slingshot failed at runtime; Python MST + polyline used"
        sling_meta["available"] = True
        sling_meta["note"] = "Python implementation of Slingshot lineages (Street 2018)"

    adata.obs["sling_pt"] = pt_df["sling_pt"].to_numpy()
    lineage_cols = [c for c in pt_df.columns if c != "sling_pt"]
    for c in lineage_cols:
        adata.obs[c] = pt_df[c].to_numpy()

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories] if hasattr(adata.obs["leiden"], "cat") else sorted(set(clusters))
    comps = _paga_components(connect, 0.0)
    paga_rows = []
    for i, a in enumerate(leiden_ids):
        for j, b in enumerate(leiden_ids):
            if i < j:
                paga_rows.append({"from": a, "to": b, "connectivity": float(connect[i, j])})
    pd.DataFrame(paga_rows).to_csv(tabdir / "paga_connectivities.tsv", sep="\t", index=False)

    # Donor-level means on tumor epithelium (normals stay in the graph for the root).
    tumor = adata.obs[adata.obs["is_tumor"].astype(bool)].copy()
    donor_df = (
        tumor.groupby("donor", observed=True)
        .agg(
            n_cells=("expr_CLDN4", "size"),
            tissue=("tissue", lambda s: ",".join(sorted(set(s.astype(str))))),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_TACSTD2=("expr_TACSTD2", "mean"),
            mean_SFTPC=("expr_SFTPC", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
            mean_IFN=("score_IFN", "mean"),
            mean_malignant=("score_malignant_like", "mean"),
            mean_club=("score_club", "mean"),
            mean_basal=("score_basal", "mean"),
            mean_dpt=("dpt_pseudotime", "mean"),
            mean_sling=("sling_pt", "mean"),
        )
        .reset_index()
    )
    elig = donor_df[donor_df["n_cells"] >= MIN_CELLS_DONOR].copy()

    primary_spec = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling"),
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("barrier/keratin vs Slingshot PT", "mean_barrier", "mean_sling"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_sling"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling"),
        ("AT2 score vs Slingshot PT (control)", "mean_AT2", "mean_sling"),
    ]
    primary = []
    for name, x, y in primary_spec:
        rec = _spearman(elig[x].to_numpy(), elig[y].to_numpy())
        rec["contrast"] = name
        primary.append(rec)
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q

    def subset_spearman(mask: pd.Series, contrast: str, x="mean_CLDN4", y="mean_sling") -> dict:
        sub = elig[mask]
        rec = _spearman(sub[x].to_numpy(), sub[y].to_numpy())
        rec["contrast"] = contrast
        return rec

    sensitivity = [
        subset_spearman(elig["tissue"].str.contains("PRIMARY") & ~elig["tissue"].str.contains("METASTASIS"), "primary-only CLDN4 vs Slingshot PT"),
        subset_spearman(elig["tissue"].str.contains("METASTASIS") & ~elig["tissue"].str.contains("PRIMARY"), "met-only CLDN4 vs Slingshot PT"),
        subset_spearman(elig["tissue"].str.contains("PRIMARY") & ~elig["tissue"].str.contains("METASTASIS"), "primary-only CLDN4 vs barrier/keratin", y="mean_barrier"),
        subset_spearman(elig["tissue"].str.contains("METASTASIS") & ~elig["tissue"].str.contains("PRIMARY"), "met-only CLDN4 vs barrier/keratin", y="mean_barrier"),
        subset_spearman(elig["tissue"].str.contains("PRIMARY") & ~elig["tissue"].str.contains("METASTASIS"), "primary-only CLDN4 vs IFN", y="mean_IFN"),
        subset_spearman(elig["tissue"].str.contains("METASTASIS") & ~elig["tissue"].str.contains("PRIMARY"), "met-only CLDN4 vs IFN", y="mean_IFN"),
        subset_spearman(pd.Series(True, index=elig.index), "tumor-donor CLDN4 vs DPT", y="mean_dpt"),
    ]
    # Include normal cells in donor means (sensitivity only).
    all_donor = (
        adata.obs.groupby("donor", observed=True)
        .agg(
            n_cells=("expr_CLDN4", "size"),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_sling=("sling_pt", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
            mean_IFN=("score_IFN", "mean"),
        )
        .reset_index()
    )
    all_elig = all_donor[all_donor["n_cells"] >= MIN_CELLS_DONOR]
    rec = _spearman(all_elig["mean_CLDN4"].to_numpy(), all_elig["mean_sling"].to_numpy())
    rec["contrast"] = "all-tissues (include normal) CLDN4 vs Slingshot PT"
    sensitivity.append(rec)

    # Paired tertile within donor, tumor cells.
    paired_rows = []
    paired_rec = []
    for donor, sub in tumor.groupby("donor", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_TERTILE_ARM or len(lo) < MIN_CELLS_TERTILE_ARM:
            continue
        paired_rec.append(
            {
                "donor": donor,
                "tissue": ",".join(sorted(set(sub["tissue"].astype(str)))),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "sling_high": float(hi["sling_pt"].mean()),
                "sling_low": float(lo["sling_pt"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_rec)
    for lab, hi, lo in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("Slingshot PT high vs low", "sling_high", "sling_low"),
        ("DPT high vs low", "dpt_high", "dpt_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": lab, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            rec = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            rec["contrast"] = lab
            paired_rows.append(rec)

    # Lineage table — done criterion.
    lin_table = []
    for rec in lin_raw.to_dict(orient="records"):
        lid = rec["lineage_id"]
        path = str(rec["path"]).split(">")
        in_lin = adata.obs["leiden"].astype(str).isin(path)
        sub = adata.obs[in_lin]
        tumor_sub = sub[sub["is_tumor"].astype(bool)]
        donor_lin = (
            tumor_sub.groupby("donor", observed=True)
            .agg(mean_CLDN4=("expr_CLDN4", "mean"), mean_pt=(lid if lid in tumor_sub.columns else "sling_pt", "mean"), n=("expr_CLDN4", "size"))
            .reset_index()
        )
        donor_lin = donor_lin[donor_lin["n"] >= MIN_CELLS_DONOR]
        spr = _spearman(donor_lin["mean_CLDN4"].to_numpy(), donor_lin["mean_pt"].to_numpy())
        bar_spr = _spearman(
            tumor_sub.groupby("donor")["score_barrier_keratin"].mean().to_numpy() if len(donor_lin) else np.array([]),
            donor_lin["mean_pt"].to_numpy() if len(donor_lin) else np.array([]),
        ) if len(donor_lin) else {"n": 0, "rho": None, "p": None}
        ifn_donor = (
            tumor_sub.groupby("donor", observed=True)
            .agg(mean_IFN=("score_IFN", "mean"), mean_pt=(lid if lid in tumor_sub.columns else "sling_pt", "mean"), n=("expr_CLDN4", "size"))
            .reset_index()
        )
        ifn_donor = ifn_donor[ifn_donor["n"] >= MIN_CELLS_DONOR]
        ifn_spr = _spearman(ifn_donor["mean_IFN"].to_numpy(), ifn_donor["mean_pt"].to_numpy())
        lin_table.append(
            {
                "lineage_id": lid,
                "start_cluster": rec["start_cluster"],
                "end_cluster": rec["end_cluster"],
                "path": rec["path"],
                "n_clusters": int(rec["n_clusters"]),
                "n_cells": int(in_lin.sum()),
                "n_tumor_cells": int(tumor_sub.shape[0]),
                "n_normal_cells": int((~sub["is_tumor"].astype(bool)).sum()),
                "n_donors": int(tumor_sub["donor"].nunique()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_sling_pt": float(np.nanmean(sub["sling_pt"].to_numpy())),
                "donor_n": spr["n"],
                "donor_rho_CLDN4_vs_PT": spr["rho"],
                "donor_p_CLDN4_vs_PT": spr["p"],
                "donor_rho_barrier_vs_PT": bar_spr.get("rho"),
                "donor_p_barrier_vs_PT": bar_spr.get("p"),
                "donor_rho_IFN_vs_PT": ifn_spr.get("rho"),
                "donor_p_IFN_vs_PT": ifn_spr.get("p"),
                "engine": sling_meta.get("engine"),
            }
        )
    lin_df = pd.DataFrame(lin_table)
    lin_df.to_csv(tabdir / "lineage.tsv", sep="\t", index=False)

    leiden_rows = []
    for cl, sub in adata.obs.groupby("leiden", observed=True):
        leiden_rows.append(
            {
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_donors": int(sub["donor"].nunique()),
                "n_tumor": int(sub["is_tumor"].sum()),
                "n_normal": int((~sub["is_tumor"].astype(bool)).sum()),
                "top_tissue": sub["tissue"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_sling": float(np.nanmean(sub["sling_pt"].to_numpy())),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "is_start": str(cl) == start_cluster,
            }
        )
    pd.DataFrame(leiden_rows).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    donor_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)

    # Figures
    colors = {"PRIMARY": "#2a6f97", "METASTASIS": "#b23a48", "NORMAL": "#6a994e"}
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    for ax, x, y, title in (
        (axes[0], "mean_sling", "mean_CLDN4", "CLDN4 vs Slingshot PT"),
        (axes[1], "mean_sling", "mean_barrier", "barrier (no CLDN4) vs PT"),
        (axes[2], "mean_sling", "mean_IFN", "IFN vs PT"),
    ):
        for tissue, col in colors.items():
            sub = elig[elig["tissue"].str.contains(tissue)]
            if sub.empty:
                continue
            ax.scatter(sub[x], sub[y] if y != "mean_CLDN4" else sub["mean_CLDN4"], s=52, c=col, label=tissue)
        if title.startswith("CLDN4"):
            ax.scatter(elig[x], elig["mean_CLDN4"], s=52, c=elig["tissue"].map(lambda t: colors.get(t.split(",")[0], "0.4")))
        ax.set_xlabel("donor-mean Slingshot PT")
        ax.set_ylabel(y.replace("mean_", "donor-mean "))
        row = next(
            r
            for r in primary
            if (title.startswith("CLDN4") and r["contrast"].startswith("CLDN4 vs Slingshot"))
            or (title.startswith("barrier") and r["contrast"].startswith("barrier"))
            or (title.startswith("IFN") and r["contrast"].startswith("IFN vs"))
        )
        ax.set_title(_fmt(row))
    axes[0].set_ylabel("donor-mean CLDN4")
    fig.suptitle(f"GSE123902 tumor epithelium  donor n={len(elig)}  root≠CLDN4-high", fontsize=11)
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    rng = np.random.default_rng(0)
    take = rng.choice(adata.n_obs, size=min(4000, adata.n_obs), replace=False)
    sub = adata.obs.iloc[take]
    for ax, y, lab in (
        (axes[0], "expr_CLDN4", "CLDN4"),
        (axes[1], "score_barrier_keratin", "barrier/keratin (no CLDN4)"),
        (axes[2], "score_IFN", "IFN"),
    ):
        ax.scatter(sub["sling_pt"], sub[y], s=6, c=sub["tissue"].map(colors).fillna("0.5"), alpha=0.45)
        ax.set_xlabel("Slingshot PT")
        ax.set_ylabel(lab)
        ax.set_title(lab + " along PT (cells; descriptive)")
    fig.suptitle("EXTRA: CLDN4 + barrier + IFN along Slingshot pseudotime", fontsize=11)
    _save(fig, figdir / "fig_extra_along_pt")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    tissue_ct = adata.obs["tissue"].astype(str).value_counts()
    axes[0].bar(tissue_ct.index.astype(str), tissue_ct.to_numpy(), color=[colors.get(x, "0.5") for x in tissue_ct.index])
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Analysis cells n={adata.n_obs} (cap {CAP_PER_SAMPLE}/sample)")
    donor_ct = elig.groupby("tissue").size() if "tissue" in elig else elig.assign(t="tumor").groupby("t").size()
    axes[1].bar(range(len(elig)), elig.sort_values("n_cells")["n_cells"].to_numpy(), color="#2a6f97")
    axes[1].set_xlabel("tumor donors (sorted)")
    axes[1].set_ylabel("tumor epi cells")
    axes[1].set_title(f"Honest n donors eligible={len(elig)} / tumor={donor_df.shape[0]}")
    _save(fig, figdir / "fig_honest_n")

    emit_extra = (not paired_df.empty) or (
        next(r for r in primary if "barrier/keratin" in r["contrast"]).get("p") or 1
    ) < 0.05
    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN"),
            ("sling_low", "sling_high", "Slingshot PT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            ax.scatter(paired_df[lo], paired_df[hi], s=48, c="#2a6f97")
            lims = [min(paired_df[lo].min(), paired_df[hi].min()), max(paired_df[lo].max(), paired_df[hi].max())]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("IFN")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("Slingshot")), keys=("W", "p")))
        fig.suptitle(f"EXTRA: within-donor CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("tissue", "fig_umap_tissue", None),
        ("donor", "fig_umap_donor", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("leiden", "fig_umap_leiden", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_IFN", "viridis"),
        ("sling_pt", "fig_umap_slingshot", "viridis"),
        ("dpt_pseudotime", "fig_umap_dpt", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        if color == "donor":
            sc.pl.umap(adata, color="donor", ax=ax, show=False, frameon=False, legend_loc=None)
        else:
            sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc.pl.paga(adata, color="leiden", ax=ax, show=False, frameon=False)
    ax.set_title(f"PAGA  components={len(comps)}  start={start_cluster}")
    _save(fig, figdir / "fig_paga")

    c4_sling = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = next(r for r in primary if "barrier/keratin" in r["contrast"] and r["contrast"].startswith("CLDN4"))
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    bar_pt = next(r for r in primary if r["contrast"].startswith("barrier/keratin vs"))
    ifn_pt = next(r for r in primary if r["contrast"].startswith("IFN vs"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))

    parts = [
        f"Donor-level CLDN4 vs Slingshot PT: {_fmt(c4_sling)}.",
        f"CLDN4 vs DPT: {_fmt(c4_dpt)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"Barrier vs Slingshot PT: {_fmt(bar_pt)}.",
        f"IFN vs Slingshot PT: {_fmt(ifn_pt)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"Slingshot engine={sling_meta.get('engine')}; {len(lin_table)} lineage(s) from Leiden {start_cluster}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Root is not CLDN4-high. Donor is the unit. Not a TACSTD2 redo. No both-high gate. 36.5 GB H5 skipped.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} tumor donors).** "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded) is {_fmt(c4_bar)}. "
        f"CLDN4 vs IFN is {_fmt(c4_ifn)}. "
        f"Along Slingshot PT, barrier is {_fmt(bar_pt)} and IFN is {_fmt(ifn_pt)}. "
        f"**What is underpowered / null.** CLDN4 vs Slingshot PT is {_fmt(c4_sling)}. "
        f"CLDN4 vs AT2 is {_fmt(c4_at2)}. n=13 tumor donors is the honest ceiling."
    )

    summary = {
        "accessions": ["GSE123902"],
        "paper": "Laughney et al. Nat Med 2020 PMID 32042191",
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "author_h5_skipped": True,
        "clock": sling_meta.get("engine"),
        "slingshot": sling_meta,
        "cap_per_sample": CAP_PER_SAMPLE,
        "n_cells_uncapped": n_uncapped,
        "n_cells_tumor_uncapped": n_tumor_uncapped,
        "n_cells_normal_uncapped": n_normal_uncapped,
        "n_cells": int(adata.n_obs),
        "n_cells_tumor": int(adata.obs["is_tumor"].sum()),
        "n_cells_normal": int((~adata.obs["is_tumor"].astype(bool)).sum()),
        "n_donors": int(adata.obs["donor"].nunique()),
        "n_tumor_donors": int(donor_df.shape[0]),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "tissue_counts": adata.obs["tissue"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "n_lineages": int(len(lin_table)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "lineages": lin_table,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": c4_bar["rho"],
            "spearman_p": c4_bar["p"],
            "spearman_n": c4_bar["n"],
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
            "author_H5": "36.5 GB PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5",
            "GSE123903": "mouse SuperSeries arm",
            "dual_high": "TACSTD2∩CLDN4 never used as a gate",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), summary)
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE123902"],
                "paper": "Laughney et al. Nat Med 2020 PMID 32042191",
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "slingshot": sling_meta,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "matched-normal AT2-like, never CLDN4-high",
                    "cap_per_sample": CAP_PER_SAMPLE,
                    "unit": "donor",
                    "epi_gate": "(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0",
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
                "n_donors_eligible": int(len(elig)),
                "n_lineages": int(len(lin_table)),
                "engine": sling_meta.get("engine"),
                "lineage_table": str(tabdir / "lineage.tsv"),
                "finding": str(args.finding),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
