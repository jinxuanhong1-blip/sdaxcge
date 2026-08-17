#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE205335 malignant cells, CLDN4 only.

ADDITIVE. Ahn/Lee ICI (RECIST). Author-malignant only. Patient is the unit.
Root is the lowest-CLDN4 Leiden cluster — never a CLDN4-high cluster.
Barrier/keratin excludes CLDN4. IFN is a locked compact ISG panel.
No TACSTD2∩CLDN4 dual-high gate. A real lineage plot is required.
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
from gene_sets import COMPARATOR, CONTROLS, FOCAL, QC_NEG, STATES  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_PATIENT = 10
MIN_CELLS_PER_TERTILE_ARM = 8
MIN_ROOT_CLUSTER_CELLS = 40
MIN_ROOT_CLUSTER_PATIENTS = 3
NSCLC = {"ADC", "SQ"}
RESPONSE_R = {"R"}
RESPONSE_NR = {"NR"}


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "available": False,
            "reason": "Rscript not on PATH",
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    env_lib = str(Path.home() / "R" / "library")
    try:
        proc = subprocess.run(
            [
                rscript,
                "-e",
                f'.libPaths(c("{env_lib}", .libPaths())); packageVersion("slingshot")',
            ],
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


def _mwu(a: np.ndarray, b: np.ndarray) -> dict:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    aa = aa[np.isfinite(aa)]
    bb = bb[np.isfinite(bb)]
    n1, n2 = int(aa.size), int(bb.size)
    if n1 < 2 or n2 < 2:
        return {
            "n_r": n1,
            "n_nr": n2,
            "median_r": None,
            "median_nr": None,
            "delta_median": None,
            "U": None,
            "p": None,
            "note": "n<2",
        }
    u, p = stats.mannwhitneyu(aa, bb, alternative="two-sided")
    return {
        "n_r": n1,
        "n_nr": n2,
        "median_r": float(np.median(aa)),
        "median_nr": float(np.median(bb)),
        "delta_median": float(np.median(aa) - np.median(bb)),
        "U": float(u),
        "p": float(p),
    }


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
    if "n_r" in row and "rho" not in keys:
        pr = "NA" if row.get("p") is None else f"{row['p']:.3g}"
        dr = "NA" if row.get("delta_median") is None else f"{row['delta_median']:.3f}"
        return f"n={row.get('n_r')} vs {row.get('n_nr')}, Δmed={dr}, p={pr}"
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
        "patient_id",
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
    info["reason"] = "harmonypy on PCA, batch=patient_id"
    return info


def pick_root_not_cldn4_high(adata) -> tuple[int, dict]:
    """External arrow: lowest-CLDN4 Leiden cluster. Never a CLDN4-high cluster."""
    obs = adata.obs
    means = (
        obs.groupby("leiden", observed=True)
        .agg(
            n_cells=("expr_CLDN4", "size"),
            n_patients=("patient_id", "nunique"),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
        )
        .sort_values("mean_CLDN4")
    )
    high_clus = str(means["mean_CLDN4"].idxmax())
    eligible = means[
        (means["n_cells"] >= MIN_ROOT_CLUSTER_CELLS)
        & (means["n_patients"] >= MIN_ROOT_CLUSTER_PATIENTS)
        & (means.index.astype(str) != high_clus)
    ]
    info = {
        "high_cldn4_cluster": high_clus,
        "high_cldn4_mean": float(means.loc[high_clus, "mean_CLDN4"]),
        "eligible_clusters": [str(x) for x in eligible.index],
        "cluster_means": {
            str(i): {
                "n_cells": int(r.n_cells),
                "n_patients": int(r.n_patients),
                "mean_CLDN4": float(r.mean_CLDN4),
            }
            for i, r in means.iterrows()
        },
    }
    if eligible.empty:
        raise SystemExit("no eligible non-CLDN4-high Leiden cluster for the root")
    root_clus = str(eligible.index[0])
    cand = obs["leiden"].astype(str) == root_clus
    scores = obs.loc[cand, "expr_CLDN4"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
    info["root_cluster"] = root_clus
    info["root_cluster_mean_CLDN4"] = float(means.loc[root_clus, "mean_CLDN4"])
    info["rule"] = (
        f"lowest-CLDN4 Leiden {root_clus} (median CLDN4 cell); "
        f"never CLDN4-high cluster {high_clus}"
    )
    if info["root_cluster"] == high_clus:
        raise SystemExit("root cluster equals CLDN4-high cluster — refused")
    return int(pick), info


def run_slingshot(adata, root_cluster: str, work: Path) -> dict:
    status = slingshot_status()
    if not status.get("available"):
        return {"ran": False, **status}
    work.mkdir(parents=True, exist_ok=True)
    rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
    pca = np.asarray(adata.obsm[rep][:, :N_PCS], dtype=float)
    umap = np.asarray(adata.obsm["X_umap"], dtype=float)
    barcodes = adata.obs_names.astype(str)
    pd.DataFrame(pca, index=barcodes, columns=[f"PC{i+1}" for i in range(pca.shape[1])]).to_csv(
        work / "pca.csv"
    )
    pd.DataFrame(umap, index=barcodes, columns=["UMAP1", "UMAP2"]).to_csv(work / "umap.csv")
    pd.DataFrame({"leiden": adata.obs["leiden"].astype(str).to_numpy()}, index=barcodes).to_csv(
        work / "cell_meta.csv"
    )
    (work / "start_cluster.txt").write_text(str(root_cluster) + "\n")
    out = work / "slingshot_out"
    out.mkdir(exist_ok=True)
    rscript = shutil.which("Rscript")
    script = Path(__file__).resolve().parent / "run_slingshot.R"
    env_lib = str(Path.home() / "R" / "library")
    proc = subprocess.run(
        [rscript, str(script), str(work), str(out)],
        check=False,
        capture_output=True,
        text=True,
        timeout=600,
        env={**dict(**{k: v for k, v in __import__("os").environ.items()}), "R_LIBS_USER": env_lib},
    )
    (out / "slingshot_stdout.txt").write_text(proc.stdout or "")
    (out / "slingshot_stderr.txt").write_text(proc.stderr or "")
    if proc.returncode != 0:
        return {
            "ran": False,
            "available": True,
            "reason": (proc.stderr or proc.stdout or "slingshot failed").strip()[:600],
            "returncode": proc.returncode,
        }
    pt = pd.read_csv(out / "slingshot_pseudotime.csv")
    pt = pt.set_index("barcode").reindex(barcodes)
    adata.obs["sling_avg_pseudotime"] = pt["sling_avg_pseudotime"].to_numpy()
    lin = pd.read_csv(out / "slingshot_lineages.csv")
    curves = pd.read_csv(out / "slingshot_umap_curves.csv")
    info_txt = (out / "slingshot_info.txt").read_text() if (out / "slingshot_info.txt").exists() else ""
    return {
        "ran": True,
        "available": True,
        "version": status.get("version"),
        "n_lineages": int(lin["lineage"].nunique()) if len(lin) else 0,
        "lineages": lin.to_dict(orient="records"),
        "n_curve_points": int(len(curves)),
        "info": info_txt,
        "stdout_tail": (proc.stdout or "")[-400:],
        "curves_path": str(out / "slingshot_umap_curves.csv"),
        "lineages_path": str(out / "slingshot_lineages.csv"),
    }


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]
    sling = s["slingshot"]
    clock = "Slingshot + scanpy DPT" if sling.get("ran") else (
        "scanpy DPT only (Slingshot did not run: "
        + str(sling.get("reason") or "unavailable")
        + ")"
    )

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    rec = s["recist"]
    lines = [
        "# Finding — GSE205335 malignant-only REAL Slingshot/PAGA, CLDN4 only",
        "",
        "ADDITIVE. **CLDN4 only.** Ahn / Lee ICI biopsy/effusion cohort, GEO "
        "[GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) "
        "(eLife 98366). Author-malignant cells only (`lineage.sub == Malignant cells`) "
        "on non-normal tissues. This folder does **not** redo the winning-pair "
        "epithelium trajectory and does **not** substitute a RECIST table for a lineage. "
        "No TACSTD2∩CLDN4 dual-high gate. GSE131907 / GSE207422 are not added.",
        "",
        f"Primary clocks: **{clock}**. Inferential unit = **patient**. "
        "Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        f"Root is {s['root']['rule']}.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- GEO patients / samples: **26 / 33**. Four normal-only patients "
        f"({', '.join(s['patients_dropped_normal_only']) or 'none'}) have 0 malignant cells and are out.",
        f"- Author-malignant cells on non-normal tissues (catalog): **{s['n_malignant_catalog']}**.",
        f"- Analysis cells after QC (capped ≤{s['cap_per_patient']}/patient): "
        f"**n_cells = {s['n_cells']}** in **n_patients = {s['n_patients']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_PATIENT} malignant cells used for Spearman: "
        f"**n = {s['n_patients_eligible']}**.",
        f"- RECIST R (PR) / NR (SD+PD) / NE among analysis patients: "
        f"**{s['n_R']} / {s['n_NR']} / {s['n_NE']}**. MPR/NMPR is unlabeled (n=0).",
        f"- NSCLC ADC+SQ patients: **{s['n_nsclc']}**. SCLC+NUT remain in the primary n.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and "
        f"CLDN4-low arms: **n = {s['n_patients_paired_tertile']}**.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- Histology (cells): {s['histology_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: ran={sling.get('ran')}; available={sling.get('available')}; "
        f"{sling.get('reason') or sling.get('version')}.",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among "
        f"{s['n_leiden']} Leiden vertices.",
        "- Do not cite n_cells as the inferential n.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Root: {s['root']['rule']} (root cell index {s['root']['index']}, "
        f"patient {s['root'].get('root_patient')}).",
        f"- Root cluster mean CLDN4 = {s['root']['root_cluster_mean_CLDN4']:.3f}; "
        f"CLDN4-high cluster {s['root']['high_cldn4_cluster']} mean = "
        f"{s['root']['high_cldn4_mean']:.3f}.",
        f"- Slingshot start.clus = {s['root']['root_cluster']}.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN: compact Hallmark-like ISG panel (STAT1, IRF1, ISG15, MX1, OAS1, IFIT1, "
        "IFIT3, IFI6, BST2, GBP1, CXCL9, CXCL10, IFI27, RSAD2, IFITM1, OAS2, IRF7, "
        "IFI44L, MX2, ISG20).",
        "",
        "## Lineage (required; not a RECIST-only table)",
        "",
        f"PAGA has **{s['n_paga_components']}** component(s) among {s['n_leiden']} Leiden vertices. "
        f"Slingshot lineages: **{s['slingshot'].get('n_lineages', 0)}**. "
        f"{s.get('lineage_text', '')}",
        "",
        "Figures: `results/figures/fig_lineage_paga_slingshot.png` (PAGA + Slingshot curves), "
        "`results/figures/fig_along_pseudotime.png` (CLDN4 + barrier + IFN vs DPT).",
        "",
        "## Primary (patient-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## DPT / Slingshot vs RECIST (patient unit; not in the BH family)",
        "",
        "| Contrast | n_R vs n_NR | median R | median NR | Δ (R−NR) | p |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in rec.get("tests", []):
        mr = "NA" if r.get("median_r") is None else f"{r['median_r']:.3f}"
        mn = "NA" if r.get("median_nr") is None else f"{r['median_nr']:.3f}"
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(
            f"| {r['contrast']} | {r.get('n_r')} vs {r.get('n_nr')} | {mr} | {mn} | {d} | {pv} |"
        )
    lines += [
        "",
        f"RECIST R/NR/NE = **{s['n_R']} / {s['n_NR']} / {s['n_NE']}**. "
        "NE is dropped from R vs NR. MPR is not labelled and is not substituted.",
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_patients | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Observed Spearman(CLDN4, barrier_keratin_no_CLDN4) n={extra['spearman_n']}, "
            f"ρ={extra_rho}, p={extra_p}."
        ),
        "",
        "| Paired contrast (high − low) | n_patients | Δ median | p |",
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
        "- This is not a winning-pair (GSE131907+GSE205335) redo and not a GSE131907 PAGA redo.",
        "- Author `Malignant cells` is **not CNV**. Residual AT2 score in malignant cells is not a normal AT2 root.",
        "- DPT / Slingshot are orderings, not clocks. Harmony on patient removes a main-effect batch; it does not prove a within-tumor differentiation axis.",
        "- RECIST 6 vs 10 is thin. MPR/NMPR is unlabeled (n=0). Do not write MPR on these figures.",
        "- SCLC / NUT stay in the primary n. ADC+SQ is a sensitivity.",
        "- No TACSTD2∩CLDN4 both-high gate. TACSTD2 is a comparator only.",
        "- Do not write “malignant cells differentiate because PAGA is connected.”",
        "- Do not write “CLDN4 marks the ICI-resistant terminal.”",
        "- A sibling RECIST-traj folder may exist; this folder still emits a real lineage plot.",
        "",
        "## Outputs",
        "",
        "- `results/tables/lineage_vertices.tsv` — **lineage table (done criterion)**",
        "- `results/tables/patient_means.tsv` — **patient table (done criterion)**",
        "- `results/tables/patient_level_spearman.tsv`",
        "- `results/tables/recist_vs_dpt.tsv`",
        "- `results/tables/paga_connectivities.tsv`",
        "- `results/tables/slingshot_lineages.tsv`",
        "- `results/figures/fig_lineage_paga_slingshot.png`",
        "- `results/figures/fig_along_pseudotime.png`",
        "- `results/figures/fig_dpt_vs_recist.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/gse205335_slingshot_real_cldn4/requirements.txt",
        "# R / Bioconductor slingshot (user library)",
        "Rscript -e 'dir.create(Sys.getenv(\"R_LIBS_USER\", \"~/R/library\"), recursive=TRUE);",
        "  .libPaths(Sys.getenv(\"R_LIBS_USER\", \"~/R/library\"));",
        "  install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\");",
        "  BiocManager::install(c(\"slingshot\",\"SingleCellExperiment\",\"DelayedMatrixStats\"), ask=FALSE, update=FALSE)'",
        "python3 methods/gse205335_slingshot_real_cldn4/scripts/download.py \\",
        "  --out /tmp/gse205335_slingshot_real",
        "python3 methods/gse205335_slingshot_real_cldn4/scripts/extract_malignant.py \\",
        "  --data /tmp/gse205335_slingshot_real \\",
        "  --out /tmp/gse205335_slingshot_real/malignant.h5ad",
        "python3 methods/gse205335_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse205335_slingshot_real/malignant.h5ad \\",
        "  --outdir methods/gse205335_slingshot_real_cldn4/results \\",
        "  --finding methods/gse205335_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/gse205335_slingshot_real_cldn4/results")
    p.add_argument("--finding", default="methods/gse205335_slingshot_real_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling_probe = slingshot_status()
    print(json.dumps({"slingshot_probe": sling_probe}, indent=2), flush=True)

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    catalog_n = int(adata.uns.get("n_malignant_catalog", adata.n_obs))
    inv_path = Path(args.input).with_suffix(".inventory.json")
    inventory = json.loads(inv_path.read_text()) if inv_path.is_file() else {}
    catalog_n = int(inventory.get("catalog", {}).get("n_malignant_non_normal", adata.n_obs))
    dropped = inventory.get("catalog", {}).get("patients_dropped_normal_only", [])

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

    root_i, root_info = pick_root_not_cldn4_high(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_patient"] = str(adata.obs.iloc[root_i].get("patient_id", ""))
    root_info["root_recist"] = str(adata.obs.iloc[root_i].get("recist", ""))
    root_info["root_histology"] = str(adata.obs.iloc[root_i].get("histology", ""))

    sling = run_slingshot(adata, root_info["root_cluster"], outdir / "slingshot_work")
    if not sling.get("ran"):
        adata.obs["sling_avg_pseudotime"] = np.nan
    print(json.dumps({k: sling[k] for k in sling if k != "lineages"}, indent=2), flush=True)

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt

    cluster_tab = []
    for cl in leiden_ids:
        sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
        cluster_tab.append(
            {
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_patients": int(sub["patient_id"].nunique()),
                "is_root": cl == root_info["root_cluster"],
                "is_cldn4_high_cluster": cl == root_info["high_cldn4_cluster"],
                "top_histology": sub["histology"].astype(str).value_counts().index[0],
                "top_recist": sub["recist"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_sling": float(sub["sling_avg_pseudotime"].mean()),
            }
        )
    lineage_df = pd.DataFrame(cluster_tab)
    lineage_df.to_csv(tabdir / "lineage_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )
    if sling.get("ran") and sling.get("lineages"):
        pd.DataFrame(sling["lineages"]).to_csv(tabdir / "slingshot_lineages.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["lineage", "order", "leiden"]).to_csv(
            tabdir / "slingshot_lineages.tsv", sep="\t", index=False
        )

    rows = []
    for unit, sub in adata.obs.groupby("patient_id", observed=True):
        rows.append(
            {
                "patient_id": unit,
                "recist": str(sub["recist"].iloc[0]),
                "response": str(sub["response"].iloc[0]),
                "histology": str(sub["histology"].iloc[0]),
                "tissue": str(sub["tissue"].iloc[0]) if "tissue" in sub.columns else "NA",
                "n_cells": int(len(sub)),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_sling": float(sub["sling_avg_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
            }
        )
    patient_df = pd.DataFrame(rows).sort_values("patient_id")
    patient_df.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    elig = patient_df[patient_df["n_cells"] >= MIN_CELLS_PER_PATIENT].copy()

    contrasts = [
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("barrier/keratin vs DPT", "mean_barrier_keratin", "mean_dpt"),
        ("IFN vs DPT", "mean_IFN", "mean_dpt"),
        ("AT2 residual vs DPT (control)", "mean_AT2", "mean_dpt"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
    ]
    primary = []
    for name, a, b in contrasts:
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        if name == "CLDN4 vs Slingshot PT" and not sling.get("ran"):
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None, "note": "slingshot not run"})
            continue
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "patient_level_spearman.tsv", sep="\t", index=False)

    r_arm = elig[elig["response"].isin(RESPONSE_R)]
    nr_arm = elig[elig["response"].isin(RESPONSE_NR)]
    recist_tests = [
        {"contrast": "DPT, RECIST R vs NR", **_mwu(r_arm["mean_dpt"], nr_arm["mean_dpt"])},
        {"contrast": "Slingshot PT, RECIST R vs NR", **_mwu(r_arm["mean_sling"], nr_arm["mean_sling"])},
        {"contrast": "CLDN4, RECIST R vs NR", **_mwu(r_arm["mean_CLDN4"], nr_arm["mean_CLDN4"])},
        {"contrast": "barrier/keratin, RECIST R vs NR", **_mwu(r_arm["mean_barrier_keratin"], nr_arm["mean_barrier_keratin"])},
        {"contrast": "IFN, RECIST R vs NR", **_mwu(r_arm["mean_IFN"], nr_arm["mean_IFN"])},
    ]
    pd.DataFrame(recist_tests).to_csv(tabdir / "recist_vs_dpt.tsv", sep="\t", index=False)

    elig_nsclc = elig[elig["histology"].isin(NSCLC)]
    elig_drop_small = elig[elig["n_cells"] >= 50]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "ADC+SQ CLDN4 vs DPT", **sp(elig_nsclc, "mean_CLDN4", "mean_dpt")},
        {"contrast": "ADC+SQ CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_nsclc, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "ADC+SQ CLDN4 vs IFN", **sp(elig_nsclc, "mean_CLDN4", "mean_IFN")},
        {"contrast": "ADC+SQ IFN vs DPT", **sp(elig_nsclc, "mean_IFN", "mean_dpt")},
        {"contrast": "drop n_cells<50 CLDN4 vs DPT", **sp(elig_drop_small, "mean_CLDN4", "mean_dpt")},
        {"contrast": "CLDN4 %pos vs DPT", **sp(elig, "pct_CLDN4_pos", "mean_dpt")},
    ]
    if sling.get("ran"):
        sensitivity.append({"contrast": "ADC+SQ CLDN4 vs Slingshot PT", **sp(elig_nsclc, "mean_CLDN4", "mean_sling")})
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("patient_id", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "patient_id": unit,
                "recist": str(sub["recist"].iloc[0]),
                "histology": str(sub["histology"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN high vs low", "IFN_high", "IFN_low"),
        ("DPT high vs low", "dpt_high", "dpt_low"),
        ("AT2 residual high vs low", "AT2_high", "AT2_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if "barrier/keratin (no CLDN4)" in r["contrast"] and r["contrast"].startswith("CLDN4"))
    emit_extra = True

    recist_color = {"PR": "#2a9d8f", "SD": "#e9c46a", "PD": "#e76f51", "NE": "#8d99ae"}
    hist_color = {"ADC": "#264653", "SQ": "#2a9d8f", "SCLC": "#e76f51", "NUT": "#9b5de5"}

    # ---- LINEAGE PLOT (required) ----
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 10.0))
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
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, frameon=False, title="UMAP Leiden + Slingshot curves", legend_loc="on data")
    curves_path = Path(sling["curves_path"]) if sling.get("curves_path") else None
    if curves_path and curves_path.is_file():
        curves = pd.read_csv(curves_path)
        for lin, sub in curves.groupby("lineage"):
            sub = sub.sort_values("step")
            ax.plot(sub["UMAP1"], sub["UMAP2"], lw=2.4, label=str(lin))
        if len(curves):
            ax.legend(fontsize=7, frameon=False, loc="best")
    root_xy = adata.obsm["X_umap"][root_i]
    ax.scatter([root_xy[0]], [root_xy[1]], s=80, c="black", marker="*", zorder=5, label="root")
    ax = axes[1, 0]
    sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False, frameon=False, cmap="magma", title="UMAP DPT (root ≠ CLDN4-high)")
    ax = axes[1, 1]
    if sling.get("ran"):
        sc.pl.umap(adata, color="sling_avg_pseudotime", ax=ax, show=False, frameon=False, cmap="magma", title="UMAP Slingshot avg PT")
    else:
        ax.text(0.5, 0.5, "Slingshot not run\n" + str(sling.get("reason", ""))[:180], ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    fig.suptitle(
        f"GSE205335 malignant REAL lineage   n_cells={adata.n_obs}  n_patients={patient_df.shape[0]}  "
        f"root=Leiden {root_info['root_cluster']}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_lineage_paga_slingshot")

    # ---- CLDN4 + barrier + IFN along pseudotime ----
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    rng = np.random.default_rng(0)
    take = rng.choice(adata.n_obs, size=min(4000, adata.n_obs), replace=False)
    x = adata.obs["dpt_pseudotime"].to_numpy()[take]
    for ax, ykey, ylab, cmap in (
        (axes[0], "expr_CLDN4", "CLDN4", "viridis"),
        (axes[1], "score_barrier_keratin", "barrier/keratin (no CLDN4)", "cividis"),
        (axes[2], "score_IFN", "IFN ISG score", "coolwarm"),
    ):
        y = adata.obs[ykey].to_numpy()[take]
        ax.scatter(x, y, s=6, c=y, cmap=cmap, alpha=0.35, linewidths=0)
        order = np.argsort(x)
        if np.isfinite(x).sum() > 50:
            xx, yy = x[order], y[order]
            mask = np.isfinite(xx) & np.isfinite(yy)
            xx, yy = xx[mask], yy[mask]
            if xx.size > 30:
                w = max(30, xx.size // 40)
                ker = np.ones(w) / w
                ax.plot(xx, np.convolve(yy, ker, mode="same"), c="black", lw=1.6)
        ax.set_xlabel("DPT")
        ax.set_ylabel(ylab)
    fig.suptitle("CLDN4 + barrier + IFN along DPT (cells; trend is descriptive)", fontsize=11)
    _save(fig, figdir / "fig_along_pseudotime")

    # patient-level along-PT companion
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    for ax, y, ylab, rowname in (
        (axes[0], "mean_CLDN4", "patient-mean CLDN4", "CLDN4 vs DPT"),
        (axes[1], "mean_barrier_keratin", "patient-mean barrier (no CLDN4)", "barrier/keratin vs DPT"),
        (axes[2], "mean_IFN", "patient-mean IFN", "IFN vs DPT"),
    ):
        for recist, col in recist_color.items():
            sub = elig[elig["recist"] == recist]
            ax.scatter(sub["mean_dpt"], sub[y], s=52, c=col, label=f"{recist} n={len(sub)}")
        row = next(r for r in primary if r["contrast"] == rowname)
        ax.set_xlabel("patient-mean DPT")
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: patient-level CLDN4 / barrier / IFN vs DPT, colored by RECIST", fontsize=11)
    _save(fig, figdir / "fig_extra_patient_along_dpt")

    # ---- DPT vs RECIST ----
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    rec_order = ["PR", "SD", "PD", "NE"]
    for ax, col, lab, test_name in (
        (axes[0], "mean_dpt", "patient-mean DPT", "DPT, RECIST R vs NR"),
        (axes[1], "mean_sling", "patient-mean Slingshot PT", "Slingshot PT, RECIST R vs NR"),
        (axes[2], "mean_CLDN4", "patient-mean CLDN4", "CLDN4, RECIST R vs NR"),
    ):
        vals = [elig.loc[elig["recist"] == r, col].to_numpy() for r in rec_order]
        bp = ax.boxplot(vals, tick_labels=rec_order, patch_artist=True)
        for patch, recist in zip(bp["boxes"], rec_order):
            patch.set_facecolor(recist_color[recist])
            patch.set_alpha(0.7)
        for recist in rec_order:
            sub = elig[elig["recist"] == recist]
            ax.scatter(np.full(len(sub), rec_order.index(recist) + 1), sub[col], s=28, c="black", zorder=3)
        row = next(r for r in recist_tests if r["contrast"] == test_name)
        ax.set_ylabel(lab)
        ax.set_title(_fmt(row, keys=("p",)))
    fig.suptitle(
        f"DPT / Slingshot / CLDN4 vs RECIST   R={len(r_arm)} vs NR={len(nr_arm)} (NE dropped from test)",
        fontsize=11,
    )
    _save(fig, figdir / "fig_dpt_vs_recist")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    hist_ct = adata.obs["histology"].astype(str).value_counts()
    axes[0].bar(hist_ct.index.astype(str), hist_ct.to_numpy(), color=[hist_color.get(h, "#888") for h in hist_ct.index])
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Malignant cells by histology  n_cells={adata.n_obs}")
    rec_ct = patient_df["recist"].value_counts().reindex(rec_order).fillna(0)
    axes[1].bar(rec_ct.index.astype(str), rec_ct.to_numpy(), color=[recist_color[r] for r in rec_ct.index])
    axes[1].set_ylabel("patients")
    axes[1].set_title(f"Honest n patients={patient_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN"),
            ("dpt_low", "dpt_high", "DPT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for recist, col in recist_color.items():
                sub = paired_df[paired_df["recist"] == recist]
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{recist} n={len(sub)}")
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
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("DPT")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("histology", "fig_umap_histology", None),
        ("recist", "fig_umap_recist", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("expr_CLDN4", "fig_umap_cldn4", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_ifn", "coolwarm"),
        ("patient_id", "fig_umap_patient", None),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        if color == "patient_id":
            kw["legend_loc"] = None
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    n_R = int((patient_df["response"] == "R").sum())
    n_NR = int((patient_df["response"] == "NR").sum())
    n_NE = int((patient_df["response"] == "NE").sum())
    n_nsclc = int(patient_df["histology"].isin(NSCLC).sum())
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    c4_sling = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    c4_bar = barrier_row
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    bar_dpt = next(r for r in primary if r["contrast"] == "barrier/keratin vs DPT")
    ifn_dpt = next(r for r in primary if r["contrast"] == "IFN vs DPT")
    rec_dpt = next(r for r in recist_tests if r["contrast"].startswith("DPT"))
    rec_sling = next(r for r in recist_tests if r["contrast"].startswith("Slingshot"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))

    lineage_text = (
        f"Root Leiden {root_info['root_cluster']} is not the CLDN4-high cluster "
        f"{root_info['high_cldn4_cluster']}. "
        + (
            f"Slingshot paths: {sling.get('info', '').strip()}"
            if sling.get("ran")
            else "Slingshot curves were not drawn; PAGA + DPT still constitute the lineage geometry."
        )
    )
    parts = [
        f"Patient-level CLDN4 vs DPT: {_fmt(c4_dpt)}.",
        f"CLDN4 vs Slingshot PT: {_fmt(c4_sling)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"barrier vs DPT: {_fmt(bar_dpt)}.",
        f"IFN vs DPT: {_fmt(ifn_dpt)}.",
        f"DPT RECIST R vs NR: {_fmt(rec_dpt, keys=('p',))}.",
        f"Slingshot PT RECIST R vs NR: {_fmt(rec_sling, keys=('p',))}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Slingshot was "
        + ("run on Harmony-PCA with the non-CLDN4-high start cluster." if sling.get("ran") else "not run; DPT is the documented clock."),
        "Not a TACSTD2 redo. No both-high gate. Malignant-only. Patient unit.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} patients).** "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded) {_fmt(c4_bar)}. "
        f"CLDN4 vs IFN {_fmt(c4_ifn)}. "
        f"Within-patient CLDN4-high vs low barrier {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN {_fmt(pair_ifn, keys=('W', 'p'))}. "
        f"**Trajectory vs programs.** CLDN4 vs DPT {_fmt(c4_dpt)}; "
        f"barrier vs DPT {_fmt(bar_dpt)}; IFN vs DPT {_fmt(ifn_dpt)}. "
        f"**RECIST.** DPT R vs NR {_fmt(rec_dpt, keys=('p',))}. "
        "n_R vs n_NR is thin; a null is inconclusive, not “CLDN4 is unrelated to ICI.”"
    )

    summary = {
        "accession": "GSE205335",
        "citation": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
        "compartment": "author lineage.sub == Malignant cells; non-normal tissues",
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "gse131907_added": False,
        "gse207422_added": False,
        "clock": "slingshot+dpt" if sling.get("ran") else "scanpy_dpt",
        "slingshot": sling,
        "harmony": harmony,
        "cap_per_patient": int(inventory.get("cap_per_patient", 400)),
        "n_malignant_catalog": catalog_n,
        "n_cells": int(adata.n_obs),
        "n_patients": int(patient_df.shape[0]),
        "n_patients_eligible": int(len(elig)),
        "n_patients_paired_tertile": int(len(paired_df)),
        "n_R": n_R,
        "n_NR": n_NR,
        "n_NE": n_NE,
        "n_nsclc": n_nsclc,
        "patients_dropped_normal_only": dropped,
        "subtype_counts": adata.obs["author_subtype"].astype(str).value_counts().to_dict(),
        "histology_counts": adata.obs["histology"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "recist": {"n_R": n_R, "n_NR": n_NR, "n_NE": n_NE, "tests": recist_tests},
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
        },
        "lineage_text": lineage_text,
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "what_holds": what_holds,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accession": "GSE205335",
                "paper": "Ahn / Lee et al. eLife 2024",
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "dual_high": False,
                "root_not_cldn4_high": True,
                "unit": "patient",
                "slingshot": {k: sling[k] for k in sling if k != "lineages"},
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "lowest-CLDN4 Leiden cluster; never CLDN4-high",
                    "cap_per_patient": int(inventory.get("cap_per_patient", 400)),
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
                "n_patients": int(patient_df.shape[0]),
                "n_eligible": int(len(elig)),
                "slingshot_ran": bool(sling.get("ran")),
                "lineage_table": str(tabdir / "lineage_vertices.tsv"),
                "patient_table": str(tabdir / "patient_means.tsv"),
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
