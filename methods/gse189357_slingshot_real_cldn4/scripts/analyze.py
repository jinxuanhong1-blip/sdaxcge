#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE189357 tumor epithelium, CLDN4 only.

ADDITIVE. Zhu et al. Exp Mol Med 2022, GSE189357, n=9 patients (AIS/MIA/IAC).
Installs/runs Bioconductor slingshot. Root is never the CLDN4-high cluster.
Patient is the unit. Barrier/keratin excludes CLDN4. No dual-high gate.
Done when results/tables/slingshot_lineages.tsv exists.
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
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STAGE_ORDER,
    STATES,
    load_ifn,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required — run scripts/install_tools.sh") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05
HERE = Path(__file__).resolve().parent
RSCRIPT = HERE / "run_slingshot.R"


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {"available": False, "reason": "Rscript not on PATH"}
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env.setdefault("R_LIBS_USER", "/tmp/r-libs")
    try:
        proc = subprocess.run(
            [rscript, "-e", '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); packageVersion("slingshot")'],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"available": False, "reason": f"Rscript slingshot probe failed: {exc}"}
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:400],
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
        return {"n": n, "rho": None, "p": None, "note": "n<4 — n may be 9"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
    return {"n": n, "rho": float(rho), "p": float(p)}


def _wilcoxon_paired(a: np.ndarray, b: np.ndarray) -> dict:
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "W": None, "p": None, "delta_median": None, "note": "n<4 — n may be 9"}
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


def maybe_harmony(adata) -> dict:
    info = {"used": False, "reason": None, "batch": "patient"}
    try:
        import harmonypy
    except ImportError:
        info["reason"] = "harmonypy not installed; neighbors on PCA"
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        return info
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :N_PCS],
        adata.obs,
        "patient",
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
    info["reason"] = (
        "harmonypy on PCA, batch=patient. Each patient is one stage "
        "(AIS/MIA/IAC), so Harmony removes between-patient/stage mean shifts; "
        "the remaining axis is within-epithelium programs (n=9)"
    )
    return info


def pick_root_cluster(adata) -> tuple[str, int, dict]:
    """Start cluster = highest AT2 among clusters that are not CLDN4-high.

    Never root on the CLDN4-highest Leiden cluster. Never pick a CLDN4-high
    tertile cell as the DPT iroot.
    """
    tab = (
        adata.obs.groupby("leiden", observed=True)
        .agg(
            n=("leiden", "size"),
            mean_AT2=("score_AT2", "mean"),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_SFTPC=("expr_SFTPC", "mean") if "expr_SFTPC" in adata.obs else ("score_AT2", "mean"),
            frac_AIS=("stage", lambda s: float((s.astype(str) == "AIS").mean())),
            frac_cldn4_high=("cldn4_tertile", lambda s: float((s.astype(str) == "high").mean())),
        )
        .reset_index()
    )
    tab["leiden"] = tab["leiden"].astype(str)
    cldn4_high_cluster = str(tab.loc[tab["mean_CLDN4"].idxmax(), "leiden"])
    cand = tab[tab["leiden"] != cldn4_high_cluster].copy()
    if cand.empty:
        raise SystemExit("all Leiden clusters are the CLDN4-high cluster")
    cand = cand.sort_values(["mean_AT2", "frac_AIS", "mean_SFTPC"], ascending=False)
    start = str(cand.iloc[0]["leiden"])
    mask = (
        (adata.obs["leiden"].astype(str) == start)
        & (adata.obs["cldn4_tertile"].astype(str) != "high")
    )
    if int(mask.sum()) < 5:
        mask = adata.obs["leiden"].astype(str) == start
    scores = adata.obs.loc[mask, "score_AT2"].to_numpy()
    idx = np.flatnonzero(mask.to_numpy())
    pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
    info = {
        "rule": (
            "Leiden with highest AT2 among clusters that are not the "
            "CLDN4-highest cluster; root cell = median AT2 in that cluster "
            "excluding CLDN4-high tertile"
        ),
        "start_cluster": start,
        "cldn4_high_cluster_excluded": cldn4_high_cluster,
        "index": int(pick),
        "root_patient": str(adata.obs.iloc[pick]["patient"]),
        "root_stage": str(adata.obs.iloc[pick]["stage"]),
        "root_cldn4_tertile": str(adata.obs.iloc[pick]["cldn4_tertile"]),
        "root_AT2": float(adata.obs.iloc[pick]["score_AT2"]),
        "root_CLDN4": float(adata.obs.iloc[pick]["expr_CLDN4"]),
        "cluster_table": tab.to_dict(orient="records"),
    }
    if str(adata.obs.iloc[pick]["cldn4_tertile"]) == "high":
        raise SystemExit("root cell is CLDN4-high — refused")
    if start == cldn4_high_cluster:
        raise SystemExit("start cluster is the CLDN4-high cluster — refused")
    return start, int(pick), info


def run_slingshot(adata, start: str, tmp: Path) -> dict:
    tmp.mkdir(parents=True, exist_ok=True)
    if "X_pca_harmony" in adata.obsm:
        rd = np.asarray(adata.obsm["X_pca_harmony"][:, :N_PCS], dtype=float)
    else:
        rd = np.asarray(adata.obsm["X_pca"][:, :N_PCS], dtype=float)
    rd_df = pd.DataFrame(rd, index=adata.obs_names.astype(str), columns=[f"PC{i+1}" for i in range(rd.shape[1])])
    cl_df = pd.DataFrame({"leiden": adata.obs["leiden"].astype(str).to_numpy()}, index=adata.obs_names.astype(str))
    rd_df.to_csv(tmp / "reduced_dim.csv")
    cl_df.to_csv(tmp / "clusters.csv")
    (tmp / "start_cluster.txt").write_text(start + "\n")
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env.setdefault("R_LIBS_USER", "/tmp/r-libs")
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript missing after install — REAL Slingshot required")
    proc = subprocess.run(
        [rscript, str(RSCRIPT), str(tmp)],
        check=False,
        capture_output=True,
        text=True,
        timeout=1800,
        env=env,
    )
    (tmp / "slingshot_stdout.txt").write_text(proc.stdout or "")
    (tmp / "slingshot_stderr.txt").write_text(proc.stderr or "")
    if proc.returncode != 0:
        raise SystemExit(
            "Slingshot R failed:\n" + (proc.stderr or proc.stdout or "no output")[-2000:]
        )
    lin = pd.read_csv(tmp / "slingshot_lineages.csv")
    pt = pd.read_csv(tmp / "slingshot_pseudotime.csv")
    pt = pt.set_index("cell")
    missing = [c for c in adata.obs_names.astype(str) if c not in pt.index]
    if missing:
        raise SystemExit(f"Slingshot dropped {len(missing)} cells")
    pt = pt.loc[adata.obs_names.astype(str)]
    lin_cols = [c for c in pt.columns if c != "cell"]
    for c in lin_cols:
        adata.obs[c] = pt[c].to_numpy()
    # Primary clock: cell-wise mean of finite lineage times (common Slingshot collapse).
    arr = pt[lin_cols].to_numpy(dtype=float)
    adata.obs["slingshot_pseudotime"] = np.nanmean(arr, axis=1)
    # CLDN4-ward lineage = lineage whose end cluster has the highest mean CLDN4.
    cl_means = adata.obs.groupby("leiden", observed=True)["expr_CLDN4"].mean()
    best_i = None
    best_end = -np.inf
    for _, row in lin.iterrows():
        end = str(row["end_cluster"])
        val = float(cl_means.get(end, -np.inf))
        if val > best_end:
            best_end = val
            best_i = str(row["lineage"])
    if best_i is not None and best_i in adata.obs:
        adata.obs["slingshot_cldn4ward"] = adata.obs[best_i]
    else:
        adata.obs["slingshot_cldn4ward"] = adata.obs["slingshot_pseudotime"]
    version = (tmp / "slingshot_version.txt").read_text().strip() if (tmp / "slingshot_version.txt").exists() else "unknown"
    return {
        "available": True,
        "version": version,
        "n_lineages": int(len(lin)),
        "lineages": lin.to_dict(orient="records"),
        "cldn4ward_lineage": best_i,
        "stdout": (proc.stdout or "").strip()[-500:],
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

    lines = [
        "# Finding — GSE189357 tumor epithelium, CLDN4-only REAL Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** Zhu et al., *Exp Mol Med* 2022 "
        "(DOI 10.1038/s12276-022-00896-9), GEO **GSE189357**: nine treatment-naïve "
        "resected LUAD lesions (TD1–TD9; AIS=3, MIA=3, IAC=3). This folder does **not** "
        "redo winning-pair Slingshot (GSE131907+GSE205335) or GSE131907-only PAGA. "
        "No TACSTD2∩CLDN4 dual-high gate. Spatial GSE189487 is not used.",
        "",
        "Primary clock: **REAL Slingshot** (Street et al. 2018; Bioconductor `slingshot`). "
        "PAGA is geometry only. Inferential unit = **patient** (one 10x sample each). "
        "**n may be 9 — say so.** Cell-level ρ is descriptive. Barrier/keratin score "
        "**excludes CLDN4**. IFN is Hallmark IFNα ∪ IFNγ. Root is the highest-AT2 Leiden "
        "cluster that is **not** the CLDN4-highest cluster; the root cell is never CLDN4-high.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- GEO catalog: **n_patients = 9** (TD1–TD9). This is the catalog n and, after "
        f"the ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN}-cell rule, the Spearman n may still be 9. **Say so.**",
        f"- Stage split: AIS={s['n_AIS']}, MIA={s['n_MIA']}, IAC={s['n_IAC']} (n=3/stage — not a stage test).",
        f"- Marker-malignant cells after QC: **n_cells = {s['n_cells']}** "
        f"(gate `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`).",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Cells per patient: {s['cells_per_patient']}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: available={s['slingshot'].get('available')}; "
        f"version={s['slingshot'].get('version')}; n_lineages={s['slingshot'].get('n_lineages')}.",
        f"- Harmony: {s['harmony']['reason']}.",
        "- No nLung / uninvolved AT2 exists in this accession. Root is AT2-high tumor epithelium, not CLDN4-high.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: Harmony on **patient** (n=9; each patient is one stage).",
        f"- Slingshot start cluster: **{s['root']['start_cluster']}** "
        f"(excluded CLDN4-high cluster {s['root']['cldn4_high_cluster_excluded']}).",
        f"- Root cell: patient {s['root']['root_patient']} stage {s['root']['root_stage']} "
        f"tertile {s['root']['root_cldn4_tertile']} (never high).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN: Hallmark IFNα ∪ IFNγ (224 genes; CLDN4 not in the set).",
        "",
        "## Slingshot lineages (done criterion)",
        "",
        "| lineage | path | n_clusters | start | end | n_cells | mean PT | mean CLDN4 | mean barrier | mean IFN |",
        "| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in s.get("lineage_table", []):
        lines.append(
            f"| {r['lineage']} | {r['clusters']} | {r['n_clusters']} | {r['start_cluster']} | "
            f"{r['end_cluster']} | {r['n_cells']} | {r.get('mean_pseudotime', float('nan')):.3f} | "
            f"{r.get('mean_CLDN4', float('nan')):.3f} | {r.get('mean_barrier_keratin', float('nan')):.3f} | "
            f"{r.get('mean_IFN', float('nan')):.3f} |"
        )
    lines += [
        "",
        f"Machine table: `results/tables/slingshot_lineages.tsv`. "
        f"CLDN4-ward lineage used for along-pseudotime plots: **{s['slingshot'].get('cldn4ward_lineage')}**.",
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
        "n may be 9. A significant p at n=9 is a small-n result, not a cohort.",
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
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired) and programs along Slingshot",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Rule: always emit extra figures on this accession. "
            f"Observed Spearman(CLDN4, barrier) n={extra['spearman_n']}, ρ={extra_rho}, p={extra_p}."
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
        "- **n may be 9.** This is not a 65-unit winning-pair result and not a stage-powered test (n=3/stage).",
        "- Harmony on patient removes between-stage mean shifts. Do not read Slingshot as AIS→IAC invasion time.",
        "- No uninvolved nLung AT2 exists here. The root is AT2-high **tumor** epithelium.",
        "- Malignant is a marker gate, **not CNV** (inferCNV was not run).",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into IAC because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot/DPT is an ordering, not a clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/slingshot_lineages.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_along_pseudotime.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash methods/gse189357_slingshot_real_cldn4/scripts/install_tools.sh",
        "python3 methods/gse189357_slingshot_real_cldn4/scripts/download.py \\",
        "  --out /tmp/gse189357",
        "python3 methods/gse189357_slingshot_real_cldn4/scripts/extract_epithelium.py \\",
        "  --tar /tmp/gse189357/GSE189357_RAW.tar \\",
        "  --out /tmp/gse189357/epithelium.h5ad",
        "python3 methods/gse189357_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse189357/epithelium.h5ad \\",
        "  --outdir methods/gse189357_slingshot_real_cldn4/results \\",
        "  --finding methods/gse189357_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/gse189357_slingshot_real_cldn4/results")
    p.add_argument("--finding", default="methods/gse189357_slingshot_real_cldn4/FINDING.md")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    tmp = outdir / "tmp"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling_probe = slingshot_status()
    print(json.dumps({"slingshot_probe": sling_probe}, indent=2), flush=True)
    if not sling_probe.get("available"):
        raise SystemExit(
            "REAL Slingshot is required. Run scripts/install_tools.sh. "
            f"Probe: {sling_probe}"
        )

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

    ifn = load_ifn()
    score_sets = dict(STATES)
    score_sets["IFN"] = ifn["IFN"]
    score_sets["IFN_core"] = ifn["IFN_core"]

    absent: dict[str, list[str]] = {}
    for name, genes in score_sets.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", [])
            if g not in absent["single_genes"]:
                absent["single_genes"].append(g)

    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"], ordered=True)
    adata.obs["stage"] = pd.Categorical(adata.obs["stage"].astype(str), categories=list(STAGE_ORDER), ordered=True)

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

    start, root_i, root_info = pick_root_cluster(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)

    sling = run_slingshot(adata, start, tmp)

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
                "n_patients": int(sub["patient"].nunique()),
                "top_stage": sub["stage"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_slingshot": float(sub["slingshot_pseudotime"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "is_start": cl == start,
                "is_cldn4_high_cluster": cl == root_info["cldn4_high_cluster_excluded"],
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    lin_df = pd.DataFrame(sling["lineages"])
    cl_map = {r["leiden"]: r for r in cluster_tab}
    extra_cols = []
    for _, row in lin_df.iterrows():
        end = str(row["end_cluster"])
        path = str(row["clusters"]).split("->")
        cells = adata.obs[adata.obs["leiden"].astype(str).isin(path)]
        extra_cols.append(
            {
                "mean_CLDN4": float(cells["expr_CLDN4"].mean()) if len(cells) else np.nan,
                "mean_barrier_keratin": float(cells["score_barrier_keratin"].mean()) if len(cells) else np.nan,
                "mean_IFN": float(cells["score_IFN"].mean()) if len(cells) else np.nan,
                "mean_AT2": float(cells["score_AT2"].mean()) if len(cells) else np.nan,
                "end_mean_CLDN4": float(cl_map.get(end, {}).get("mean_CLDN4", np.nan)),
                "end_mean_IFN": float(cl_map.get(end, {}).get("mean_IFN", np.nan)),
            }
        )
    lin_df = pd.concat([lin_df.reset_index(drop=True), pd.DataFrame(extra_cols)], axis=1)
    lin_df.to_csv(tabdir / "slingshot_lineages.tsv", sep="\t", index=False)

    rows = []
    for unit, sub in adata.obs.groupby("patient", observed=True):
        rows.append(
            {
                "patient": unit,
                "unit_id": unit,
                "stage": str(sub["stage"].iloc[0]),
                "n_cells": int(len(sub)),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()),
                "mean_slingshot": float(sub["slingshot_pseudotime"].mean()),
                "mean_slingshot_cldn4ward": float(sub["slingshot_cldn4ward"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_IFN_core": float(sub["score_IFN_core"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_slingshot"),
        ("CLDN4 vs Slingshot CLDN4-ward PT", "mean_CLDN4", "mean_slingshot_cldn4ward"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN (Hallmark α∪γ)", "mean_CLDN4", "mean_IFN"),
        ("barrier/keratin vs Slingshot PT", "mean_barrier_keratin", "mean_slingshot"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_slingshot"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_slingshot"),
        ("AT2 score vs Slingshot PT (control)", "mean_AT2", "mean_slingshot"),
    ]
    primary = []
    for name, a, b in contrasts:
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    elig_no_iac = elig[elig["stage"] != "IAC"]
    elig_no_ais = elig[elig["stage"] != "AIS"]
    sensitivity = [
        {"contrast": "drop IAC CLDN4 vs Slingshot PT", **sp(elig_no_iac, "mean_CLDN4", "mean_slingshot")},
        {"contrast": "drop AIS CLDN4 vs Slingshot PT", **sp(elig_no_ais, "mean_CLDN4", "mean_slingshot")},
        {"contrast": "CLDN4 vs DPT (companion)", **sp(elig, "mean_CLDN4", "mean_dpt")},
        {"contrast": "IFN core vs Slingshot PT", **sp(elig, "mean_IFN_core", "mean_slingshot")},
        {"contrast": "CLDN4 vs IFN core", **sp(elig, "mean_CLDN4", "mean_IFN_core")},
        {"contrast": "CLDN4 vs malignant-like", **sp(elig, "mean_CLDN4", "mean_malignant_like")},
        {"contrast": "IFN vs Slingshot CLDN4-ward PT", **sp(elig, "mean_IFN", "mean_slingshot_cldn4ward")},
        {"contrast": "barrier vs Slingshot CLDN4-ward PT", **sp(elig, "mean_barrier_keratin", "mean_slingshot_cldn4ward")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("patient", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "patient": unit,
                "stage": str(sub["stage"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "sling_high": float(hi["slingshot_pseudotime"].mean()),
                "sling_low": float(lo["slingshot_pseudotime"].mean()),
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
        ("Slingshot PT high vs low", "sling_high", "sling_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if "barrier/keratin (no CLDN4)" in r["contrast"] and r["contrast"].startswith("CLDN4"))
    emit_extra = True

    stage_colors = {"AIS": "#4c9f70", "MIA": "#e09f3e", "IAC": "#9b2226"}

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=axes[0, 0],
        show=False,
        frameon=False,
        cmap="viridis",
        title="PAGA (Leiden) colored by mean CLDN4",
    )
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[0, 1], show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    sc.pl.umap(
        adata,
        color="slingshot_pseudotime",
        ax=axes[1, 0],
        show=False,
        frameon=False,
        cmap="magma",
        title=f"UMAP Slingshot PT (start={start})",
    )
    ax = axes[1, 1]
    c4_sl = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    for st, col in stage_colors.items():
        sub = elig[elig["stage"] == st]
        ax.scatter(sub["mean_slingshot"], sub["mean_CLDN4"], s=70, c=col, label=f"{st} n={len(sub)}")
        for _, rec in sub.iterrows():
            ax.annotate(rec["patient"], (rec["mean_slingshot"], rec["mean_CLDN4"]), fontsize=7, xytext=(4, 2), textcoords="offset points")
    ax.set_xlabel("patient-mean Slingshot PT")
    ax.set_ylabel("patient-mean CLDN4")
    ax.set_title(f"CLDN4 vs Slingshot  {_fmt(c4_sl)}  (n may be 9)")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"GSE189357 tumor epithelium scored by CLDN4   n_cells={adata.n_obs}  n_patients={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # Extra: CLDN4 + barrier + IFN along Slingshot (cells + patient means)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    pt = adata.obs["slingshot_cldn4ward"].to_numpy()
    finite = np.isfinite(pt)
    bins = np.linspace(np.nanmin(pt[finite]), np.nanmax(pt[finite]), 12) if finite.any() else np.array([0, 1])
    for ax, key, lab in (
        (axes[0], "expr_CLDN4", "CLDN4"),
        (axes[1], "score_barrier_keratin", "barrier/keratin (no CLDN4)"),
        (axes[2], "score_IFN", "IFN (Hallmark α∪γ)"),
    ):
        y = adata.obs[key].to_numpy()
        ax.scatter(pt[finite], y[finite], s=3, c="0.75", alpha=0.25, rasterized=True)
        if finite.sum() > 20:
            idx = np.digitize(pt[finite], bins)
            xs, ys = [], []
            for b in range(1, len(bins)):
                m = idx == b
                if m.sum() >= 8:
                    xs.append(0.5 * (bins[b - 1] + bins[b]))
                    ys.append(float(np.nanmean(y[finite][m])))
            if xs:
                ax.plot(xs, ys, c="#1d3557", lw=2)
        for st, col in stage_colors.items():
            sub = elig[elig["stage"] == st]
            ax.scatter(
                sub["mean_slingshot_cldn4ward"],
                sub[{"expr_CLDN4": "mean_CLDN4", "score_barrier_keratin": "mean_barrier_keratin", "score_IFN": "mean_IFN"}[key]],
                s=55,
                c=col,
                edgecolors="k",
                linewidths=0.4,
                label=st,
                zorder=3,
            )
        ax.set_xlabel("Slingshot PT (CLDN4-ward)")
        ax.set_ylabel(lab)
        ax.set_title(lab)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: CLDN4 + barrier + IFN along Slingshot (cells grey; patients colored by stage)", fontsize=11)
    _save(fig, figdir / "fig_along_pseudotime")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = adata.obs.groupby(["stage", "patient"], observed=True).size().unstack(fill_value=0)
    ct.T.plot(kind="bar", ax=axes[0], color=[stage_colors[s] for s in ct.index])
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Cells / patient  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8, title="stage")
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    stage_n = sample_df.groupby("stage").size().reindex(list(STAGE_ORDER)).fillna(0)
    axes[1].bar(stage_n.index.astype(str), stage_n.to_numpy(), color=[stage_colors[s] for s in stage_n.index])
    axes[1].set_ylabel("patients")
    axes[1].set_title(f"Honest n patients={sample_df.shape[0]} (eligible {len(elig)}; n may be 9)")
    for i, v in enumerate(stage_n.to_numpy()):
        axes[1].text(i, v + 0.05, str(int(v)), ha="center", fontsize=9)
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)", "barrier"),
            ("IFN_low", "IFN_high", "IFN", "IFN"),
            ("sling_low", "sling_high", "Slingshot PT", "Slingshot"),
        )
        for ax, (lo, hi, lab, key) in zip(axes, panels):
            for st, col in stage_colors.items():
                sub = paired_df[paired_df["stage"] == st]
                ax.scatter(sub[lo], sub[hi], s=50, c=col, label=f"{st} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
            ax.set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith(key.split()[0] if key != "Slingshot" else "Slingshot")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)} (n may be 9)", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("stage", "fig_umap_stage", None),
        ("patient", "fig_umap_patient", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_IFN", "viridis"),
        ("slingshot_pseudotime", "fig_umap_slingshot", "magma"),
        ("dpt_pseudotime", "fig_umap_dpt", "magma"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_ifn = next(r for r in primary if r["contrast"].startswith("CLDN4 vs IFN"))
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_AT2", "mean_CLDN4", c4_at2, "patient-mean AT2", "patient-mean CLDN4"),
        (axes[1], "mean_barrier_keratin", "mean_CLDN4", barrier_row, "patient-mean barrier (no CLDN4)", "patient-mean CLDN4"),
        (axes[2], "mean_IFN", "mean_CLDN4", c4_ifn, "patient-mean IFN", "patient-mean CLDN4"),
    ):
        for st, col in stage_colors.items():
            sub = elig[elig["stage"] == st]
            ax.scatter(sub[x], sub[y], s=55, c=col, label=f"{st} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: patient-level CLDN4 vs AT2 / barrier / IFN  (n may be 9)", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_programs")

    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    cells_per_patient = sample_df.set_index("patient")["n_cells"].to_dict()
    c4_sl = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    c4_bar = barrier_row
    c4_ifn = next(r for r in primary if r["contrast"].startswith("CLDN4 vs IFN"))
    bar_sl = next(r for r in primary if r["contrast"].startswith("barrier/keratin vs Slingshot"))
    ifn_sl = next(r for r in primary if r["contrast"] == "IFN vs Slingshot PT")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    at2_sl = next(r for r in primary if r["contrast"].startswith("AT2 score vs Slingshot"))

    parts = [
        f"Patient-level CLDN4 vs Slingshot PT: {_fmt(c4_sl)} (n may be 9).",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN (Hallmark α∪γ): {_fmt(c4_ifn)}.",
        f"barrier vs Slingshot PT: {_fmt(bar_sl)}.",
        f"IFN vs Slingshot PT: {_fmt(ifn_sl)}.",
        f"AT2 vs Slingshot PT (control): {_fmt(at2_sl)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Slingshot lineages={sling['n_lineages']}; start={start}; excluded CLDN4-high cluster={root_info['cldn4_high_cluster_excluded']}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "REAL Slingshot was run. Patient is the unit. n may be 9. Not a TACSTD2 redo. No both-high gate.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} patients; n may be 9).** "
        f"CLDN4 vs Slingshot PT {_fmt(c4_sl)}. "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded) {_fmt(c4_bar)}. "
        f"CLDN4 vs IFN {_fmt(c4_ifn)}. "
        f"Within-patient CLDN4-high vs low barrier {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN {_fmt(pair_ifn, keys=('W', 'p'))}. "
        f"**Limits.** Catalog n=9 (3 AIS / 3 MIA / 3 IAC). Harmony on patient removes "
        f"between-stage mean shifts. No nLung AT2. Marker-malignant ≠ CNV."
    )

    summary = {
        "accession": "GSE189357",
        "paper": "Zhu et al. Exp Mol Med 2022 DOI 10.1038/s12276-022-00896-9",
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": "slingshot",
        "slingshot": sling,
        "harmony": harmony,
        "n_cells": int(adata.n_obs),
        "n_units": int(sample_df.shape[0]),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "n_AIS": int((sample_df["stage"] == "AIS").sum()),
        "n_MIA": int((sample_df["stage"] == "MIA").sum()),
        "n_IAC": int((sample_df["stage"] == "IAC").sum()),
        "n_may_be_9": True,
        "cells_per_patient": cells_per_patient,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": {k: v for k, v in root_info.items() if k != "cluster_table"},
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "lineage_table": lin_df.to_dict(orient="records"),
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
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accession": "GSE189357",
                "paper": "Zhu et al. Exp Mol Med 2022 DOI 10.1038/s12276-022-00896-9",
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn": "Hallmark IFNalpha union IFNgamma",
                "slingshot": sling,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "highest-AT2 Leiden that is not CLDN4-highest; never CLDN4-high cell",
                    "unit": "patient",
                    "n_catalog": 9,
                },
            },
            indent=2,
        )
    )
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": int(adata.n_obs),
                "n_patients": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "n_lineages": sling["n_lineages"],
                "lineage_table": str(tabdir / "slingshot_lineages.tsv"),
                "finding": args.finding,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
