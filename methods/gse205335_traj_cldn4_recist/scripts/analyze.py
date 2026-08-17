#!/usr/bin/env python3
"""PAGA + DPT on GSE205335 epithelium, CLDN4 and RECIST at the patient.

ADDITIVE. GSE205335 only. Public processed UMI. No GSE148071. No dual-high.
Slingshot R if available, else scanpy DPT. Root is biological (AT2-like or
author non-malignant / low-CLDN4), never CLDN4-high.
Inferential unit = patient. Cell-level ρ is descriptive.
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
    AUTHOR_MALIGNANT,
    AUTHOR_NONMALIGNANT,
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
MIN_CELLS_PER_PATIENT = 10
MIN_CELLS_PER_TERTILE_ARM = 8
MIN_ROOT_CELLS = 20


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
            timeout=20,
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
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:300],
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
            "n1": n1,
            "n2": n2,
            "n": n1 + n2,
            "U": None,
            "p": None,
            "r_rb": None,
            "delta_median": None,
            "median_1": None,
            "median_2": None,
            "note": "n<2 per arm",
        }
    u, p = stats.mannwhitneyu(aa, bb, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n1 * n2) - 1.0
    return {
        "n1": n1,
        "n2": n2,
        "n": n1 + n2,
        "U": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "delta_median": float(np.median(aa) - np.median(bb)),
        "median_1": float(np.median(aa)),
        "median_2": float(np.median(bb)),
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


def _fmt_spearman(row: dict) -> str:
    n = row.get("n")
    if row.get("rho") is None:
        return f"n={n}, ρ=NA, p=NA"
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _fmt_mwu(row: dict, arm1="PR", arm2="PD/SD") -> str:
    if row.get("p") is None:
        return f"n={row.get('n')}, {arm1} n={row.get('n1')}, {arm2} n={row.get('n2')}, p=NA"
    return (
        f"n={row['n']} ({arm1} {row['n1']} vs {arm2} {row['n2']}), "
        f"Δmed={row['delta_median']:+.3f}, r={row['r_rb']:+.3f}, p={row['p']:.3g}"
    )


def _fmt_paired(row: dict) -> str:
    if row.get("W") is None:
        return f"n={row.get('n')}, W=NA, p=NA"
    return f"n={row['n']}, W={row['W']:.1f}, Δmed={row['delta_median']:+.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _pick_root(adata) -> tuple[int, dict]:
    """Biological root. Never CLDN4-high.

    Preference:
    1. Author Non-malignant cells with AT2 score at the median, if that pool
       is not CLDN4-high relative to the epithelial object.
    2. Leiden cluster with the highest mean AT2 score among clusters whose
       mean CLDN4 is at or below the object median (low-CLDN4 / AT2-like).
    3. If every cluster is CLDN4-high vs the median (degenerate), use the
       Leiden cluster with the lowest mean CLDN4.
    """
    cldn4 = adata.obs["expr_CLDN4"].to_numpy(dtype=float)
    at2 = adata.obs["score_AT2"].to_numpy(dtype=float)
    cldn4_med = float(np.nanmedian(cldn4))
    subtype = adata.obs["author_subtype"].astype(str)
    nonmal = subtype.isin(AUTHOR_NONMALIGNANT)
    info = {
        "n_author_nonmalignant": int(nonmal.sum()),
        "n_author_malignant": int(subtype.isin(AUTHOR_MALIGNANT).sum()),
        "object_median_CLDN4": cldn4_med,
        "rule": None,
        "rejected": [],
    }

    if int(nonmal.sum()) >= MIN_ROOT_CELLS:
        pool_cldn4 = float(np.nanmean(cldn4[nonmal.to_numpy()]))
        if pool_cldn4 <= cldn4_med:
            idx = np.flatnonzero(nonmal.to_numpy())
            scores = at2[idx]
            pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
            info["rule"] = (
                "author Non-malignant cells (median AT2 score); "
                "pool mean CLDN4 at or below object median"
            )
            info["root_pool_mean_CLDN4"] = pool_cldn4
            info["root_pool_mean_AT2"] = float(np.nanmean(at2[idx]))
            return int(pick), info
        info["rejected"].append(
            f"author Non-malignant pool mean CLDN4={pool_cldn4:.3f} > object median {cldn4_med:.3f}"
        )

    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    leiden = adata.obs["leiden"].astype(str)
    rows = []
    for cl in sorted(leiden.unique(), key=lambda x: int(x) if str(x).isdigit() else str(x)):
        mask = leiden.eq(cl).to_numpy()
        rows.append(
            {
                "leiden": cl,
                "n": int(mask.sum()),
                "mean_CLDN4": float(np.nanmean(cldn4[mask])),
                "mean_AT2": float(np.nanmean(at2[mask])),
                "frac_nonmalignant": float(nonmal.to_numpy()[mask].mean()) if mask.any() else 0.0,
            }
        )
    clus = pd.DataFrame(rows)
    info["leiden_means"] = clus.to_dict(orient="records")
    low = clus.loc[clus["mean_CLDN4"] <= cldn4_med].copy()
    if not low.empty:
        low = low.sort_values(["mean_AT2", "mean_CLDN4"], ascending=[False, True])
        top = low.iloc[0]
        rule = (
            f"Leiden {top['leiden']} highest AT2 among clusters with mean CLDN4 "
            f"≤ object median (AT2-like / low-CLDN4)"
        )
    else:
        top = clus.sort_values("mean_CLDN4", ascending=True).iloc[0]
        rule = (
            f"Leiden {top['leiden']} lowest mean CLDN4 "
            "(all clusters were CLDN4-high vs object median)"
        )
    cand = leiden.eq(str(top["leiden"])).to_numpy()
    idx = np.flatnonzero(cand)
    scores = at2[idx]
    pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
    info["rule"] = rule
    info["root_leiden"] = str(top["leiden"])
    info["root_cluster_mean_CLDN4"] = float(top["mean_CLDN4"])
    info["root_cluster_mean_AT2"] = float(top["mean_AT2"])
    return int(pick), info


def write_finding(path: Path, s: dict) -> None:
    sling = s["slingshot"]
    clock = (
        "Slingshot"
        if sling.get("available")
        else "documented AT2-like / low-CLDN4 diffusion pseudotime (scanpy DPT; Slingshot R missing)"
    )

    def spear_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    def recist_md(r: dict) -> str:
        dmed = "NA" if r.get("delta_median") is None else f"{r['delta_median']:+.3f}"
        rrb = "NA" if r.get("r_rb") is None else f"{r['r_rb']:+.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        return (
            f"| {r['contrast']} | {r['n']} | {r.get('n1', 'NA')} | {r.get('n2', 'NA')} | "
            f"{dmed} | {rrb} | {pv} |"
        )

    lines = [
        "# Finding — GSE205335 epithelium trajectory, CLDN4 and RECIST",
        "",
        "ADDITIVE. **CLDN4 only.** GSE205335 advanced NSCLC ICI (Ahn / Lee, *eLife* 2024), "
        "public processed UMI. **No GSE148071. No dual-high.** Does **not** re-audit "
        "PR #320 T/NK ρ. Patient is the unit.",
        "",
        f"Primary clock: **{clock}**. Root is a biological AT2-like / low-CLDN4 "
        "state, **never CLDN4-high**. Barrier/keratin score **excludes CLDN4**. "
        "RECIST test is **PR vs PD/SD** (NE excluded). Cell-level ρ is descriptive.",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- GEO catalog epithelium on non-normal tissues: **{s['n_cells_catalog']}** cells / "
        f"**{s['n_patients_catalog']}** patients.",
        f"- Analysis cells after QC (capped ≤{s['cap_per_patient']}/patient): "
        f"**n_cells = {s['n_cells']}**.",
        f"- Patients on the object: **n_patients = {s['n_patients']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_PATIENT} epithelial cells (Spearman): "
        f"**n = {s['n_patients_eligible']}**.",
        f"- RECIST-evaluable (PR vs PD/SD, NE out): **n = {s['n_recist']}** "
        f"(PR {s['n_pr']} vs PD/SD {s['n_pd_sd']}; NE {s['n_ne']}).",
        f"- Author malignant / non-malignant cells on the object: "
        f"{s['n_author_malignant']} / {s['n_author_nonmalignant']}.",
        f"- Patients with ≥{MIN_CELLS_PER_PATIENT} author-malignant cells: "
        f"**n = {s['n_patients_malignant_eligible']}**.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- Histology (patients): {s['histology_patients']}.",
        f"- RECIST (patients): {s['recist_patients']}.",
        f"- CLDN4 tertile cells: {s['cldn4_tertile_counts']}.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and "
        f"CLDN4-low arms: **n = {s['n_patients_paired']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: available={sling.get('available')}; {sling.get('reason', sling.get('version', ''))}.",
        "- GSE148071 not used. GSE131907 not used. Dual-high not used. PR #320 T/NK not re-tested.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {s['leiden_resolution']}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        "- No Harmony (single cohort; Harmony on patient would erase the RECIST contrast).",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, "
        f"patient {s['root']['patient']}, subtype {s['root']['author_subtype']}, "
        f"CLDN4={s['root']['expr_CLDN4']:.3f}, AT2={s['root']['score_AT2']:.3f}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among "
        f"{s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- RECIST arms: PR vs PD+SD. NE is excluded from the RECIST table, kept in CLDN4–DPT.",
        "",
        "## Primary — sample-level CLDN4 vs DPT (Spearman, BH inside this list)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(spear_md(r))
    lines += [
        "",
        "## Primary — sample-level DPT vs RECIST (PR vs PD/SD)",
        "",
        "Mann–Whitney on patient means. Positive Δmed / r = PR higher than PD/SD. "
        "NE excluded. Honest n is patients, not cells.",
        "",
        "| Contrast | n | n_PR | n_PD_SD | Δmed (PR−PD/SD) | r_rb | p |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_recist"]:
        lines.append(recist_md(r))
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_patients | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["sensitivity_spearman"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    lines += [
        "",
        "| Contrast | n | n_PR | n_PD_SD | Δmed (PR−PD/SD) | r_rb | p |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in s["sensitivity_recist"]:
        lines.append(recist_md(r))
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)",
        "",
        f"Emitted: **True**. Paired tertile n={s['n_patients_paired']}.",
        "",
        "| Paired contrast (high − low) | n_patients | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["paired_tertile"]:
        dmed = "NA" if r.get("delta_median") is None else f"{r['delta_median']:+.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {dmed} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- If CLDN4 vs DPT is null, that is **not** lineage proof. The RECIST split is still reported.",
        "- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Author Malignant cells is the published label, **not CNV re-called here**.",
        "- Q4 mixes SCLC with ADC/SQ; histology is reported, not hidden.",
        "- This folder does **not** re-audit PR #320 T/NK ρ.",
        "- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.",
        "- Slingshot R was not the clock unless `available=True`. DPT is an ordering, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/sample_level_spearman.tsv` — CLDN4 vs DPT (done criterion)",
        "- `results/tables/sample_level_recist.tsv` — DPT vs RECIST (done criterion)",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_trajectory_recist.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_extra_sample_cldn4_dpt_recist.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/gse205335_traj_cldn4_recist/requirements.txt",
        "python3 methods/gse205335_traj_cldn4_recist/scripts/download.py \\",
        "  --out /tmp/gse205335_traj",
        "python3 methods/gse205335_traj_cldn4_recist/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse205335_traj \\",
        "  --out /tmp/gse205335_traj/epithelium.h5ad",
        "python3 methods/gse205335_traj_cldn4_recist/scripts/analyze.py \\",
        "  --input /tmp/gse205335_traj/epithelium.h5ad \\",
        "  --outdir methods/gse205335_traj_cldn4_recist/results \\",
        "  --finding methods/gse205335_traj_cldn4_recist/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path("/tmp/gse205335_traj/epithelium.h5ad"))
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("methods/gse205335_traj_cldn4_recist/results"),
    )
    p.add_argument(
        "--finding",
        type=Path,
        default=Path("methods/gse205335_traj_cldn4_recist/FINDING.md"),
    )
    p.add_argument("--inventory", type=Path, default=None)
    args = p.parse_args()

    outdir = args.outdir
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sling = slingshot_status()
    print(json.dumps({"slingshot": sling}, indent=2), flush=True)

    adata = sc.read_h5ad(args.input)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    catalog_cells = int(adata.n_obs)
    catalog_patients = int(adata.obs["patient_id"].nunique())
    inv_path = args.inventory or args.input.with_suffix(".inventory.json")
    catalog_from_inv = {}
    if inv_path.is_file():
        catalog_from_inv = json.loads(inv_path.read_text())
        catalog_cells = int(
            catalog_from_inv.get("catalog", {}).get("n_epithelial_non_normal", catalog_cells)
        )
        catalog_patients = int(catalog_from_inv.get("n_patients", catalog_patients))

    sc.pp.filter_cells(adata, min_genes=200)
    if "n_counts" not in adata.obs:
        adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata = adata[adata.obs["n_counts"] >= 500].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent = {}
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
            absent.setdefault("focal", []).append(g)

    q1, q2 = np.nanquantile(adata.obs["expr_CLDN4"].to_numpy(), [1 / 3, 2 / 3])
    adata.obs["cldn4_tertile"] = pd.cut(
        adata.obs["expr_CLDN4"],
        bins=[-np.inf, q1, q2, np.inf],
        labels=["low", "mid", "high"],
    ).astype(str)

    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts"
    )
    sc.pp.pca(adata, n_comps=N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph")
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.umap(adata)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata)
    root_idx, root_info = _pick_root(adata)
    adata.uns["iroot"] = int(root_idx)
    sc.tl.dpt(adata)
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].astype(float)

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    leiden_ids = [str(x) for x in adata.obs["leiden"].astype(str).unique()]
    comps = _paga_components(connect, thresh=0.0)

    root_cell = adata.obs.iloc[root_idx]
    root_info.update(
        {
            "index": int(root_idx),
            "patient": str(root_cell["patient_id"]),
            "author_subtype": str(root_cell["author_subtype"]),
            "leiden": str(root_cell["leiden"]),
            "expr_CLDN4": float(root_cell["expr_CLDN4"]),
            "score_AT2": float(root_cell["score_AT2"]),
            "recist": str(root_cell.get("recist", "NA")),
            "histology": str(root_cell.get("histology", "NA")),
        }
    )
    print(json.dumps({"root": {k: v for k, v in root_info.items() if k != "leiden_means"}}, indent=2), flush=True)

    def _patient_frame(mask=None) -> pd.DataFrame:
        obs = adata.obs if mask is None else adata.obs.loc[mask]
        g = obs.groupby("patient_id", observed=True)
        out = g.agg(
            n_cells=("expr_CLDN4", "size"),
            n_malignant=("is_author_malignant", lambda s: int((s.astype(str) == "True").sum())),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_TACSTD2=("expr_TACSTD2", "mean"),
            mean_DPT=("dpt_pseudotime", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_club=("score_club", "mean"),
            mean_basal=("score_basal", "mean"),
            mean_barrier_keratin=("score_barrier_keratin", "mean"),
            mean_malignant_like=("score_malignant_like", "mean"),
            mean_SFTPC=("expr_SFTPC", "mean"),
            recist=("recist", "first"),
            recist_arm=("recist_arm", "first"),
            histology=("histology", "first"),
            tissue=("tissue", "first"),
        ).reset_index()
        return out

    sample_df = _patient_frame()
    elig = sample_df.loc[sample_df["n_cells"] >= MIN_CELLS_PER_PATIENT].copy()
    recist_elig = elig.loc[elig["recist_arm"].isin(["PR", "PD_SD"])].copy()
    mal_mask = adata.obs["is_author_malignant"].astype(str).eq("True")
    mal_df = _patient_frame(mal_mask)
    mal_elig = mal_df.loc[mal_df["n_cells"] >= MIN_CELLS_PER_PATIENT].copy()
    mal_recist = mal_elig.loc[mal_elig["recist_arm"].isin(["PR", "PD_SD"])].copy()

    nsclc = elig.loc[elig["histology"].isin(["ADC", "SQ"])].copy()
    nsclc_recist = nsclc.loc[nsclc["recist_arm"].isin(["PR", "PD_SD"])].copy()
    drop_sclc = elig.loc[elig["histology"] != "SCLC"].copy()
    drop_sclc_recist = drop_sclc.loc[drop_sclc["recist_arm"].isin(["PR", "PD_SD"])].copy()
    recist_only_spear = recist_elig

    primary = [
        {"contrast": "CLDN4 vs DPT", **_spearman(elig["mean_CLDN4"], elig["mean_DPT"])},
        {"contrast": "CLDN4 vs AT2 score", **_spearman(elig["mean_CLDN4"], elig["mean_AT2"])},
        {"contrast": "CLDN4 vs club score", **_spearman(elig["mean_CLDN4"], elig["mean_club"])},
        {"contrast": "CLDN4 vs basal score", **_spearman(elig["mean_CLDN4"], elig["mean_basal"])},
        {
            "contrast": "CLDN4 vs barrier/keratin (no CLDN4)",
            **_spearman(elig["mean_CLDN4"], elig["mean_barrier_keratin"]),
        },
        {
            "contrast": "CLDN4 vs malignant-like score",
            **_spearman(elig["mean_CLDN4"], elig["mean_malignant_like"]),
        },
        {
            "contrast": "CLDN4 vs TACSTD2 (comparator)",
            **_spearman(elig["mean_CLDN4"], elig["mean_TACSTD2"]),
        },
        {"contrast": "SFTPC vs DPT (control)", **_spearman(elig["mean_SFTPC"], elig["mean_DPT"])},
        {"contrast": "AT2 score vs DPT (control)", **_spearman(elig["mean_AT2"], elig["mean_DPT"])},
    ]
    qvals = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qvals):
        r["q"] = None if r["p"] is None else q

    def _recist_row(frame: pd.DataFrame, col: str, label: str) -> dict:
        pr = frame.loc[frame["recist_arm"] == "PR", col].to_numpy()
        nrs = frame.loc[frame["recist_arm"] == "PD_SD", col].to_numpy()
        return {"contrast": label, **_mwu(pr, nrs)}

    primary_recist = [
        _recist_row(recist_elig, "mean_DPT", "DPT vs RECIST (PR vs PD/SD)"),
        _recist_row(recist_elig, "mean_CLDN4", "CLDN4 vs RECIST (PR vs PD/SD)"),
        _recist_row(recist_elig, "mean_AT2", "AT2 vs RECIST (PR vs PD/SD)"),
        _recist_row(
            recist_elig, "mean_barrier_keratin", "barrier/keratin vs RECIST (PR vs PD/SD)"
        ),
    ]

    sensitivity = [
        {
            "contrast": "RECIST-evaluable only CLDN4 vs DPT",
            **_spearman(recist_only_spear["mean_CLDN4"], recist_only_spear["mean_DPT"]),
        },
        {
            "contrast": "malignant-only CLDN4 vs DPT",
            **_spearman(mal_elig["mean_CLDN4"], mal_elig["mean_DPT"]),
        },
        {
            "contrast": "ADC+SQ CLDN4 vs DPT",
            **_spearman(nsclc["mean_CLDN4"], nsclc["mean_DPT"]),
        },
        {
            "contrast": "drop-SCLC CLDN4 vs DPT",
            **_spearman(drop_sclc["mean_CLDN4"], drop_sclc["mean_DPT"]),
        },
    ]
    sensitivity_recist = [
        _recist_row(mal_recist, "mean_DPT", "malignant-only DPT vs RECIST"),
        _recist_row(nsclc_recist, "mean_DPT", "ADC+SQ DPT vs RECIST"),
        _recist_row(drop_sclc_recist, "mean_DPT", "drop-SCLC DPT vs RECIST"),
        _recist_row(mal_recist, "mean_CLDN4", "malignant-only CLDN4 vs RECIST"),
        _recist_row(nsclc_recist, "mean_CLDN4", "ADC+SQ CLDN4 vs RECIST"),
    ]

    # paired tertiles within patient
    paired_rows = []
    paired_records = []
    for pid, sub in adata.obs.groupby("patient_id", observed=True):
        hi = sub.loc[sub["cldn4_tertile"] == "high"]
        lo = sub.loc[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        rec = {
            "patient_id": pid,
            "recist": str(sub["recist"].iloc[0]),
            "recist_arm": str(sub["recist_arm"].iloc[0]),
            "histology": str(sub["histology"].iloc[0]),
            "n_high": int(len(hi)),
            "n_low": int(len(lo)),
            "barrier_high": float(hi["score_barrier_keratin"].mean()),
            "barrier_low": float(lo["score_barrier_keratin"].mean()),
            "AT2_high": float(hi["score_AT2"].mean()),
            "AT2_low": float(lo["score_AT2"].mean()),
            "dpt_high": float(hi["dpt_pseudotime"].mean()),
            "dpt_low": float(lo["dpt_pseudotime"].mean()),
            "mal_high": float(hi["score_malignant_like"].mean()),
            "mal_low": float(lo["score_malignant_like"].mean()),
        }
        paired_records.append(rec)
    paired_df = pd.DataFrame(paired_records)
    if not paired_df.empty:
        paired_rows = [
            {
                "contrast": "barrier/keratin (no CLDN4) high vs low",
                **_wilcoxon_paired(paired_df["barrier_high"], paired_df["barrier_low"]),
            },
            {
                "contrast": "AT2 high vs low",
                **_wilcoxon_paired(paired_df["AT2_high"], paired_df["AT2_low"]),
            },
            {
                "contrast": "malignant-like high vs low",
                **_wilcoxon_paired(paired_df["mal_high"], paired_df["mal_low"]),
            },
            {
                "contrast": "DPT high vs low",
                **_wilcoxon_paired(paired_df["dpt_high"], paired_df["dpt_low"]),
            },
        ]
    else:
        paired_rows = [
            {"contrast": name, "n": 0, "W": None, "p": None, "delta_median": None}
            for name in (
                "barrier/keratin (no CLDN4) high vs low",
                "AT2 high vs low",
                "malignant-like high vs low",
                "DPT high vs low",
            )
        ]

    # tables
    spear_tab = pd.DataFrame(primary)
    spear_tab.to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)
    recist_tab = pd.DataFrame(primary_recist + sensitivity_recist)
    recist_tab.to_csv(tabdir / "sample_level_recist.tsv", sep="\t", index=False)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    leiden_means = (
        adata.obs.groupby("leiden", observed=True)
        .agg(
            n_cells=("expr_CLDN4", "size"),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_DPT=("dpt_pseudotime", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
            frac_PR=("recist_arm", lambda s: float((s == "PR").mean())),
            frac_PD_SD=("recist_arm", lambda s: float((s == "PD_SD").mean())),
            frac_malignant=(
                "is_author_malignant",
                lambda s: float((s.astype(str) == "True").mean()),
            ),
        )
        .reset_index()
    )
    leiden_means.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=range(connect.shape[0]), columns=range(connect.shape[1])).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    # figures
    recist_colors = {"PR": "#2a9d8f", "PD_SD": "#e76f51", "NE_or_other": "#8d99ae"}
    hist_colors = {"ADC": "#264653", "SQ": "#e9c46a", "SCLC": "#e76f51", "NUT": "#6d597a"}

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[0], show=False, frameon=False, cmap="viridis")
    axes[0].set_title("CLDN4")
    sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[1], show=False, frameon=False, cmap="magma")
    axes[1].set_title("DPT")
    sc.pl.paga(adata, ax=axes[2], show=False, frameon=False)
    axes[2].set_title(f"PAGA  {len(comps)} component(s)")
    fig.suptitle(
        f"GSE205335 epithelium  n_cells={adata.n_obs}  n_patients={sample_df.shape[0]}  "
        f"root={root_info['rule'][:60]}",
        fontsize=10,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    sc.pl.umap(adata, color="recist_arm", ax=axes[0], show=False, frameon=False)
    axes[0].set_title("RECIST arm (cells)")
    for arm, col in recist_colors.items():
        sub = recist_elig if arm != "NE_or_other" else elig[elig["recist_arm"] == "NE_or_other"]
        if arm != "NE_or_other":
            sub = recist_elig[recist_elig["recist_arm"] == arm]
        axes[1].scatter(sub["mean_DPT"], sub["mean_CLDN4"], s=56, c=col, label=f"{arm} n={len(sub)}")
    axes[1].set_xlabel("patient-mean DPT")
    axes[1].set_ylabel("patient-mean CLDN4")
    axes[1].set_title(_fmt_spearman(primary[0]))
    axes[1].legend(fontsize=7, frameon=False)
    # box of DPT by RECIST
    boxes, labels, cols = [], [], []
    for arm, lab in (("PR", "PR"), ("PD_SD", "PD/SD")):
        vals = recist_elig.loc[recist_elig["recist_arm"] == arm, "mean_DPT"].to_numpy()
        boxes.append(vals)
        labels.append(f"{lab}\nn={len(vals)}")
        cols.append(recist_colors[arm])
    bp = axes[2].boxplot(boxes, patch_artist=True, widths=0.55)
    axes[2].set_xticklabels(labels)
    for patch, col in zip(bp["boxes"], cols):
        patch.set_facecolor(col)
        patch.set_alpha(0.7)
    axes[2].set_ylabel("patient-mean DPT")
    axes[2].set_title(_fmt_mwu(primary_recist[0]))
    fig.suptitle("GSE205335 RECIST along DPT (patient unit)", fontsize=11)
    _save(fig, figdir / "fig_trajectory_recist")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ct = (
        adata.obs.groupby(["author_subtype", "recist_arm"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.plot(kind="bar", ax=axes[0], color=[recist_colors.get(c, "#999") for c in ct.columns])
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=30, ha="right")
    recist_n = elig.groupby("recist").size()
    axes[1].bar(recist_n.index.astype(str), recist_n.to_numpy(), color="#2a6f97")
    axes[1].set_ylabel("patients")
    axes[1].set_title(
        f"Honest n patients={elig.shape[0]}  RECIST-evaluable={len(recist_elig)}"
    )
    _save(fig, figdir / "fig_honest_n")

    if not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("dpt_low", "dpt_high", "DPT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for arm, col in recist_colors.items():
                sub = paired_df[paired_df["recist_arm"] == arm]
                if sub.empty:
                    continue
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{arm} n={len(sub)}")
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
        axes[0].set_title(_fmt_paired(paired_rows[0]))
        axes[1].set_title(_fmt_paired(paired_rows[1]))
        axes[2].set_title(_fmt_paired(paired_rows[3]))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(
            f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)}",
            fontsize=11,
        )
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for ax, x, y, row, xlab, ylab in (
        (
            axes[0],
            "mean_DPT",
            "mean_CLDN4",
            primary[0],
            "patient-mean DPT",
            "patient-mean CLDN4",
        ),
        (
            axes[1],
            "mean_barrier_keratin",
            "mean_CLDN4",
            primary[4],
            "patient-mean barrier/keratin (no CLDN4)",
            "patient-mean CLDN4",
        ),
    ):
        for hist, col in hist_colors.items():
            sub = elig[elig["histology"] == hist]
            if sub.empty:
                continue
            marker = {"PR": "o", "PD_SD": "s", "NE_or_other": "^"}
            for arm, mk in marker.items():
                ss = sub[sub["recist_arm"] == arm]
                if ss.empty:
                    continue
                ax.scatter(
                    ss[x],
                    ss[y],
                    s=52,
                    c=col,
                    marker=mk,
                    label=f"{hist}/{arm} n={len(ss)}",
                )
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt_spearman(row))
        ax.legend(fontsize=6, frameon=False)
    fig.suptitle("EXTRA: patient-level CLDN4 vs DPT / barrier (RECIST+histology)", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_dpt_recist")

    for color, fname, cmap in (
        ("histology", "fig_umap_histology", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("recist", "fig_umap_recist", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("patient_id", "fig_umap_patient", None),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        if color == "patient_id":
            kw["legend"] = False
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    subtype_counts = adata.obs["author_subtype"].astype(str).value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    histology_patients = elig["histology"].value_counts().to_dict()
    recist_patients = elig["recist"].value_counts().to_dict()
    n_pr = int((recist_elig["recist_arm"] == "PR").sum())
    n_pd_sd = int((recist_elig["recist_arm"] == "PD_SD").sum())
    n_ne = int((elig["recist_arm"] == "NE_or_other").sum())

    c4_dpt = primary[0]
    rec_dpt = primary_recist[0]
    rec_c4 = primary_recist[1]
    pair_bar = paired_rows[0]
    pair_dpt = paired_rows[3]
    parts = [
        f"Patient-level CLDN4 vs AT2-like/low-CLDN4 DPT: {_fmt_spearman(c4_dpt)}.",
        f"CLDN4 vs AT2 score: {_fmt_spearman(primary[1])}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt_spearman(primary[4])}.",
        f"CLDN4 vs malignant-like: {_fmt_spearman(primary[5])}.",
        f"DPT vs RECIST (PR vs PD/SD): {_fmt_mwu(rec_dpt)}.",
        f"CLDN4 vs RECIST (PR vs PD/SD): {_fmt_mwu(rec_c4)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt_paired(pair_bar)}.",
        f"Paired CLDN4-high vs low DPT: {_fmt_paired(pair_dpt)}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Slingshot R was not run unless available; DPT is the documented clock. "
        "Not a TACSTD2 redo. No both-high gate. GSE148071 not used. PR #320 T/NK not re-tested.",
    ]
    if c4_dpt.get("p") is None or (c4_dpt.get("p") is not None and c4_dpt["p"] >= 0.05):
        parts.insert(
            1,
            "CLDN4 vs DPT is null at the patient — this is not lineage proof; RECIST is still reported.",
        )
    verdict = " ".join(parts)

    dpt_null = c4_dpt.get("p") is None or (c4_dpt.get("p") is not None and c4_dpt["p"] >= 0.05)
    rec_sig = rec_dpt.get("p") is not None and rec_dpt["p"] < 0.05
    what_holds = (
        f"**What holds (patient n={c4_dpt['n']}; RECIST n={rec_dpt['n']}, "
        f"PR {n_pr} vs PD/SD {n_pd_sd}).** "
        f"CLDN4 vs DPT: {_fmt_spearman(c4_dpt)}. "
        f"DPT vs RECIST (PR vs PD/SD): {_fmt_mwu(rec_dpt)}. "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt_spearman(primary[4])}. "
    )
    if dpt_null:
        what_holds += (
            "**What does not hold.** Patient-level CLDN4 vs DPT is null, so this is "
            "**not** a lineage / differentiation claim. The RECIST split is still the deliverable. "
        )
    if not rec_sig:
        what_holds += (
            "DPT does not significantly separate PR vs PD/SD at this n; that is reported, not hidden. "
        )
    what_holds += (
        f"Root used: {root_info['rule']}. "
        "Do not re-read this as a T/NK result (PR #320 is not re-audited)."
    )

    summary = {
        "accession": "GSE205335",
        "gse148071_used": False,
        "gse131907_used": False,
        "dual_high": False,
        "reaudit_pr320_tnk": False,
        "primary_gene": "CLDN4",
        "clock": "scanpy_dpt" if not sling.get("available") else "slingshot",
        "slingshot": sling,
        "cap_per_patient": int(catalog_from_inv.get("cap_per_patient", 500)),
        "n_cells_catalog": catalog_cells,
        "n_patients_catalog": catalog_patients,
        "n_cells": int(adata.n_obs),
        "n_patients": int(sample_df.shape[0]),
        "n_patients_eligible": int(len(elig)),
        "n_patients_malignant_eligible": int(len(mal_elig)),
        "n_patients_paired": int(len(paired_df)),
        "n_recist": int(len(recist_elig)),
        "n_pr": n_pr,
        "n_pd_sd": n_pd_sd,
        "n_ne": n_ne,
        "n_author_malignant": int(mal_mask.sum()),
        "n_author_nonmalignant": int((~mal_mask).sum()),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "histology_patients": histology_patients,
        "recist_patients": recist_patients,
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "primary_recist": primary_recist,
        "sensitivity_spearman": sensitivity,
        "sensitivity_recist": sensitivity_recist,
        "paired_tertile": paired_rows,
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean())
            if "expr_EPCAM" in adata.obs
            else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean())
            if "expr_PTPRC" in adata.obs
            else None,
        },
        "verdict": verdict,
        "what_holds": what_holds,
        "skipped": {
            "GSE148071": "explicitly not used",
            "GSE131907": "explicitly not used",
            "dual_high": "no TACSTD2∩CLDN4 gate",
            "PR320_TNK": "not re-audited",
            "slingshot": sling.get("reason"),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), summary)
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accession": "GSE205335",
                "paper": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "recist_contrast": "PR vs PD/SD",
                "unit": "patient",
                "slingshot": sling,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": root_info["rule"],
                    "cap_per_patient": summary["cap_per_patient"],
                    "harmony": False,
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
                "n_patients": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "n_recist": int(len(recist_elig)),
                "root": root_info["rule"],
                "finding": str(args.finding),
                "spearman_table": str(tabdir / "sample_level_spearman.tsv"),
                "recist_table": str(tabdir / "sample_level_recist.tsv"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
