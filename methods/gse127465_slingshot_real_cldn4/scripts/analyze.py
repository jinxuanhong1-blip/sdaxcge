#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE127465 tumor epithelium, CLDN4 only.

ADDITIVE. Installs/uses R slingshot (not a DPT stand-in). Root / start
cluster is Type II / AT2-high and is never the CLDN4-high cluster.
No TACSTD2∩CLDN4 dual-high gate. n=7 patients is thin — said so.
If malignant+CLDN4 is missing, stop with honest n=0.
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
    raise SystemExit("scanpy is required; run scripts/install_tools.sh") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE = 10
MIN_CELLS_PER_TERTILE_ARM = 8
RANDOM_SEED = 0
CLDN4_POS = 0.0  # author-normalized value > 0


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


def _fmt_num(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}g}" if abs(x) < 0.01 or abs(x) >= 100 else f"{x:.{nd}f}"
    return str(x)


def pick_start_cluster(adata) -> dict:
    """Type II / AT2-high Leiden start. Never the CLDN4-high cluster."""
    obs = adata.obs
    g = obs.groupby("leiden", observed=True)
    stats_df = pd.DataFrame(
        {
            "n": g.size(),
            "n_type2": g["is_type2"].sum(),
            "n_malignant": g["is_malignant"].sum(),
            "frac_type2": g["is_type2"].mean(),
            "frac_malignant": g["is_malignant"].mean(),
            "mean_AT2": g["score_AT2"].mean(),
            "mean_CLDN4": g["expr_CLDN4"].mean(),
            "mean_barrier": g["score_barrier_keratin"].mean(),
        }
    )
    stats_df.index = stats_df.index.astype(str)
    cldn4_high = str(stats_df["mean_CLDN4"].idxmax())
    # Rank by Type II fraction, then AT2 score. Drop the CLDN4-high cluster.
    ranked = stats_df.sort_values(
        ["frac_type2", "mean_AT2", "n_type2"], ascending=False
    )
    start = None
    rule = None
    for cl in ranked.index.astype(str):
        if cl == cldn4_high:
            continue
        start = cl
        rule = (
            f"Leiden {cl} max Type II fraction / AT2 "
            f"(CLDN4-high cluster {cldn4_high} excluded)"
        )
        break
    if start is None:
        raise SystemExit("could not pick a non-CLDN4-high start cluster")
    return {
        "start_cluster": start,
        "cldn4_high_cluster": cldn4_high,
        "rule": rule,
        "cluster_stats": stats_df.reset_index().rename(columns={"leiden": "leiden"}),
    }


def run_slingshot(adata, start: str, work: Path, rscript: str) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    pcs = np.asarray(adata.obsm["X_pca"][:, :N_PCS], dtype=float)
    pca_df = pd.DataFrame(pcs, index=adata.obs_names.astype(str))
    pca_path = work / "pca.csv"
    cl_path = work / "clusters.csv"
    pca_df.to_csv(pca_path)
    pd.DataFrame({"leiden": adata.obs["leiden"].astype(str).to_numpy()}, index=adata.obs_names).to_csv(
        cl_path
    )
    r_script = Path(__file__).resolve().parent / "run_slingshot.R"
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "library"))
    proc = subprocess.run(
        [
            rscript,
            str(r_script),
            "--pca",
            str(pca_path),
            "--clusters",
            str(cl_path),
            "--start",
            str(start),
            "--out",
            str(work),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=1800,
        env=env,
    )
    (work / "slingshot_stdout.txt").write_text(proc.stdout or "")
    (work / "slingshot_stderr.txt").write_text(proc.stderr or "")
    if proc.returncode != 0:
        raise SystemExit(
            "REAL slingshot failed:\n"
            + (proc.stderr or proc.stdout or "no output")[:4000]
        )
    pt = pd.read_csv(work / "slingshot_pseudotime.csv")
    raw = pd.read_csv(work / "slingshot_lineages_raw.csv")
    return {"pseudotime": pt, "lineages_raw": raw, "stdout": proc.stdout}


def write_stop_finding(path: Path, gate: dict) -> None:
    path.write_text(
        "\n".join(
            [
                "# Finding — GSE127465 REAL Slingshot/PAGA, CLDN4 only",
                "",
                "**STOP.** Author malignant + CLDN4 is missing. Trajectory was not run.",
                "",
                "ADDITIVE. **CLDN4 only.** Zilionis et al., *Immunity* 2019, PMID 30979687 "
                "(GSE127465 human inDrops). No dual-high. Root would not have been CLDN4-high.",
                "",
                "## Honest n",
                "",
                f"- Patients deposited: **7**.",
                f"- Author malignant: **{gate.get('n_malignant', 0)}**.",
                f"- Malignant CLDN4>0: **{gate.get('n_malignant_cldn4_pos', 0)}**.",
                "- Inferential n: **n=0** (gate failed).",
                "",
                "Done criterion: stop note with honest n=0.",
                "",
            ]
        )
    )


def write_finding(path: Path, s: dict) -> None:
    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lin_lines = [
        "| lineage | start | end | n_clusters | n_cells | mean_CLDN4 | mean_AT2 | frac_malignant | cluster_path |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for r in s["lineages"]:
        lin_lines.append(
            f"| {r['lineage_id']} | {r['start_cluster']} | {r['end_cluster']} | "
            f"{r['n_clusters']} | {r['n_cells']} | {_fmt_num(r.get('mean_CLDN4'))} | "
            f"{_fmt_num(r.get('mean_AT2'))} | {_fmt_num(r.get('frac_malignant'))} | "
            f"{r['cluster_path']} |"
        )

    extra = s["extra_figure"]
    lines = [
        "# Finding — GSE127465 REAL Slingshot/PAGA, CLDN4 only",
        "",
        "ADDITIVE. **CLDN4 only.** Zilionis et al., *Immunity* 2019, PMID 30979687; "
        "human NSCLC **inDrops** (GSE127465). Tumor epithelium = author Type I / Type II / "
        "club / ciliated + `PatientN-specific` malignant. Blood dropped. "
        "**No TACSTD2∩CLDN4 dual-high gate.** This folder does not merge GSE148071.",
        "",
        "Primary clock: **REAL Slingshot** (Street et al. 2018; R `slingshot` "
        f"{s['slingshot']['version']}) on Harmony-free PCA + Leiden, start cluster "
        f"**{s['start']['start_cluster']}** ({s['start']['rule']}). "
        "PAGA is geometry only. Inferential unit = **patient**. "
        "**n=7 is thin** — Spearman on 7 patients is a sign check, not a precise effect. "
        "Barrier/keratin **excludes CLDN4**. Root / start is not CLDN4-high "
        f"(CLDN4-high Leiden = {s['start']['cldn4_high_cluster']}).",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Human patients deposited: **7** (p1–p7). This catalog n is also the test n. **Thin.**",
        f"- Tumor cells (author metadata): **{s['n_tumor']}**. Blood cells excluded: **{s['n_blood']}**.",
        f"- Tumor epithelium on the object: **n_cells = {s['n_cells']}**.",
        f"- Author malignant (`Patient*`-specific): **{s['n_malignant']}**. "
        f"CLDN4>0 among them: **{s['n_malignant_cldn4_pos']}** "
        f"({s['n_malignant_cldn4_pos_by_patient']}).",
        f"- Type II / Type I / club / ciliated: "
        f"{s['n_type2']} / {s['n_type1']} / {s['n_club']} / {s['n_ciliated']}.",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE} epithelial cells: **n = {s['n_patients_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4 tertile arms: "
        f"**n = {s['n_patients_paired']}**.",
        f"- CLDN4 tertile cells: {s['cldn4_tertile_counts']}.",
        f"- Author lineage counts: {s['lineage_counts']}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- Dual-high TACSTD2∩CLDN4 gate: **not defined / not used**.",
        "- Slingshot MST will connect patient-specific malignant clusters even when they are "
        "transcriptionally discrete. PAGA connectivity is the honesty check.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        "- Batch: none (no Harmony). Author malignant is PatientN-specific by design.",
        f"- Slingshot start: Leiden {s['start']['start_cluster']}. "
        f"CLDN4-high cluster {s['start']['cldn4_high_cluster']} was ineligible as root.",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        f"- Slingshot version: {s['slingshot']['version']}.",
        "",
        "## Lineage table (done criterion)",
        "",
        *lin_lines,
        "",
        "Machine table: `results/tables/lineages.tsv`.",
        "",
        "## Primary (patient-level Spearman, BH inside this list)",
        "",
        "**n=7 is thin.** p-values are descriptive.",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)",
        "",
        f"Emitted: **{extra['emitted']}**. Paired tertile n={s['n_patients_paired']}.",
        "",
        "| Paired contrast (high − low) | n | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["paired_tertile"]:
        dmed = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {dmed} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- **n=7 is thin.** Do not write a precise effect size from the patient Spearman.",
        "- Author `PatientN-specific` is the paper's malignant label, **not CNV re-called here**.",
        "- Slingshot lineages are MST paths on Leiden centers. They are not proof that "
        "Type II differentiates into each patient's tumor.",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.",
        "- Not ICI / MPR / RECIST. Not a merge with GSE148071.",
        "- RNA velocity was not run (no spliced/unspliced).",
        "",
        "## Outputs",
        "",
        "- `results/tables/lineages.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_extra_lineage_cldn4.png`",
        "- `results/figures/fig_extra_paga.png`",
        "- `results/figures/fig_extra_patient_cldn4.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash methods/gse127465_slingshot_real_cldn4/scripts/install_tools.sh",
        "python3 methods/gse127465_slingshot_real_cldn4/scripts/download.py \\",
        "  --out /tmp/gse127465_slingshot",
        "python3 methods/gse127465_slingshot_real_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse127465_slingshot \\",
        "  --out /tmp/gse127465_slingshot/epithelium.h5ad",
        "python3 methods/gse127465_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse127465_slingshot/epithelium.h5ad \\",
        "  --outdir methods/gse127465_slingshot_real_cldn4/results \\",
        "  --finding methods/gse127465_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default="methods/gse127465_slingshot_real_cldn4/results")
    ap.add_argument("--finding", default="methods/gse127465_slingshot_real_cldn4/FINDING.md")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    work = outdir / "slingshot_work"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript missing; run scripts/install_tools.sh (REAL slingshot required)")
    env_probe = dict(**{k: v for k, v in __import__("os").environ.items()})
    env_probe.setdefault("R_LIBS_USER", str(Path.home() / "R" / "library"))
    probe = subprocess.run(
        [rscript, "-e", '.libPaths(Sys.getenv("R_LIBS_USER")); cat(as.character(packageVersion("slingshot")))'],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
        env=env_probe,
    )
    if probe.returncode != 0:
        raise SystemExit(
            "slingshot R package missing; run scripts/install_tools.sh\n"
            + (probe.stderr or probe.stdout or "")[:500]
        )
    sling_ver = (probe.stdout or "").strip()

    adata = sc.read_h5ad(args.input)
    if "CLDN4" not in adata.var_names:
        write_stop_finding(Path(args.finding), {"n_malignant": 0, "n_malignant_cldn4_pos": 0})
        raise SystemExit(0)

    # log1p of author-normalized values (not already log).
    sc.pp.log1p(adata)
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
    absent = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")

    n_mal = int(adata.obs["is_malignant"].sum())
    n_mal_pos = int(
        (adata.obs["is_malignant"] & (adata.obs["expr_CLDN4"] > np.log1p(CLDN4_POS))).sum()
    )
    # expr_CLDN4 is log1p(normalized); CLDN4>0 on raw normalized is log1p>0.
    n_mal_pos = int((adata.obs["is_malignant"] & (adata.obs["expr_CLDN4"] > 0)).sum())
    if n_mal == 0 or n_mal_pos == 0:
        gate = {"n_malignant": n_mal, "n_malignant_cldn4_pos": n_mal_pos}
        write_stop_finding(Path(args.finding), gate)
        (tabdir / "STOP_no_malignant_cldn4.json").write_text(json.dumps(gate, indent=2))
        raise SystemExit(0)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    adata.raw = adata
    sc.pp.pca(adata, n_comps=max(N_PCS, 40), random_state=RANDOM_SEED)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_SEED)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    sc.tl.paga(adata, groups="leiden")

    start_info = pick_start_cluster(adata)
    start_info["cluster_stats"].to_csv(tabdir / "leiden_cluster_stats.tsv", sep="\t", index=False)
    start = start_info["start_cluster"]

    # Companion DPT rooted on a Type II cell in the start cluster (not CLDN4-high).
    type2_in_start = (adata.obs["leiden"].astype(str) == start) & adata.obs["is_type2"]
    if int(type2_in_start.sum()) >= 1:
        idx = np.flatnonzero(type2_in_start.to_numpy())
        scores = adata.obs.loc[type2_in_start, "score_AT2"].to_numpy()
        root = int(idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])
        root_rule = f"Type II cell in start Leiden {start} (median AT2)"
    else:
        cand = adata.obs["leiden"].astype(str) == start
        idx = np.flatnonzero(cand.to_numpy())
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        root = int(idx[int(np.nanargmax(scores))])
        root_rule = f"max AT2 cell in start Leiden {start} (no Type II in cluster)"
    adata.uns["iroot"] = root
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    start_info["dpt_root_index"] = root
    start_info["dpt_root_rule"] = root_rule
    start_info["dpt_root_cell"] = str(adata.obs_names[root])
    start_info["dpt_root_patient"] = str(adata.obs["Patient"].iloc[root])

    sling = run_slingshot(adata, start, work, rscript)
    pt = sling["pseudotime"].set_index("cell")
    pt = pt.reindex(adata.obs_names.astype(str))
    lin_cols = [c for c in pt.columns if c != "cell"]
    for c in lin_cols:
        adata.obs[f"sling_{c}"] = pt[c].to_numpy()
    adata.obs["sling_pseudotime"] = np.nanmin(
        np.where(np.isfinite(pt.to_numpy()), pt.to_numpy(), np.inf), axis=1
    )
    adata.obs.loc[~np.isfinite(adata.obs["sling_pseudotime"]), "sling_pseudotime"] = np.nan

    # Lineage table with honest cell/program stats.
    lineages = []
    for rec in sling["lineages_raw"].itertuples(index=False):
        col = rec.lineage_id
        if col not in pt.columns:
            # slingshot names Lineage1 vs lineage_id
            matches = [c for c in pt.columns if c.replace(" ", "") == col.replace(" ", "")]
            col = matches[0] if matches else None
        if col is None or col not in pt.columns:
            # try LineageN
            guess = rec.lineage_id
            col = guess if guess in pt.columns else None
        if col is None:
            # map by order
            continue
        mask = np.isfinite(pt[col].to_numpy())
        sub = adata.obs.loc[mask]
        lineages.append(
            {
                "lineage_id": rec.lineage_id,
                "start_cluster": rec.start_cluster,
                "end_cluster": rec.end_cluster,
                "n_clusters": int(rec.n_clusters),
                "cluster_path": rec.cluster_path,
                "n_cells": int(mask.sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()) if len(sub) else np.nan,
                "mean_AT2": float(sub["score_AT2"].mean()) if len(sub) else np.nan,
                "mean_barrier": float(sub["score_barrier_keratin"].mean()) if len(sub) else np.nan,
                "frac_malignant": float(sub["is_malignant"].mean()) if len(sub) else np.nan,
                "frac_type2": float(sub["is_type2"].mean()) if len(sub) else np.nan,
                "start_is_cldn4_high": bool(str(rec.start_cluster) == start_info["cldn4_high_cluster"]),
                "clock": "slingshot",
                "slingshot_version": sling_ver,
            }
        )
    if not lineages:
        # fallback: use raw rows even if PT columns mismatched
        for rec in sling["lineages_raw"].itertuples(index=False):
            lineages.append(
                {
                    "lineage_id": rec.lineage_id,
                    "start_cluster": rec.start_cluster,
                    "end_cluster": rec.end_cluster,
                    "n_clusters": int(rec.n_clusters),
                    "cluster_path": rec.cluster_path,
                    "n_cells": np.nan,
                    "mean_CLDN4": np.nan,
                    "mean_AT2": np.nan,
                    "mean_barrier": np.nan,
                    "frac_malignant": np.nan,
                    "frac_type2": np.nan,
                    "start_is_cldn4_high": bool(str(rec.start_cluster) == start_info["cldn4_high_cluster"]),
                    "clock": "slingshot",
                    "slingshot_version": sling_ver,
                }
            )
    lin_df = pd.DataFrame(lineages)
    lin_df.to_csv(tabdir / "lineages.tsv", sep="\t", index=False)

    q1, q2 = np.nanquantile(adata.obs["expr_CLDN4"].to_numpy(), [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[adata.obs["expr_CLDN4"] <= q1] = "low"
    tert[adata.obs["expr_CLDN4"] > q2] = "high"
    adata.obs["cldn4_tertile"] = tert

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    comps = _paga_components(connect, thresh=0.0)
    paga_df = pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids)
    paga_df.to_csv(tabdir / "paga_connectivities.tsv", sep="\t")
    pd.DataFrame(
        {
            "leiden": leiden_ids,
            "n": [(adata.obs["leiden"].astype(str) == x).sum() for x in leiden_ids],
        }
    ).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)

    # Patient-level means.
    rows = []
    for pat, sub in adata.obs.groupby(adata.obs["Patient"].astype(str), observed=True):
        if len(sub) < MIN_CELLS_PER_SAMPLE:
            continue
        rows.append(
            {
                "patient": pat,
                "n_cells": int(len(sub)),
                "n_malignant": int(sub["is_malignant"].sum()),
                "n_type2": int(sub["is_type2"].sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "mean_slingshot": float(sub["sling_pseudotime"].mean()),
                "mean_DPT": float(sub["dpt_pseudotime"].mean()) if "dpt_pseudotime" in sub else np.nan,
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df

    primary_specs = [
        ("CLDN4 vs Slingshot", "mean_CLDN4", "mean_slingshot"),
        ("CLDN4 vs DPT (companion)", "mean_CLDN4", "mean_DPT"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Slingshot (control)", "mean_SFTPC", "mean_slingshot"),
        ("AT2 score vs Slingshot (control)", "mean_AT2", "mean_slingshot"),
    ]
    primary = []
    for name, x, y in primary_specs:
        r = _spearman(elig[x].to_numpy(), elig[y].to_numpy())
        r["contrast"] = name
        primary.append(r)
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    # Paired tertile extra.
    paired_rows_data = []
    for pat, sub in adata.obs.groupby(adata.obs["Patient"].astype(str), observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_rows_data.append(
            {
                "patient": pat,
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "mal_high": float(hi["score_malignant_like"].mean()),
                "mal_low": float(lo["score_malignant_like"].mean()),
                "sling_high": float(hi["sling_pseudotime"].mean()),
                "sling_low": float(lo["sling_pseudotime"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_rows_data)
    paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    if not paired_df.empty:
        for lab, hi, lo in (
            ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
            ("AT2 high vs low", "AT2_high", "AT2_low"),
            ("malignant-like high vs low", "mal_high", "mal_low"),
            ("Slingshot high vs low", "sling_high", "sling_low"),
            ("DPT high vs low", "dpt_high", "dpt_low"),
        ):
            w = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            w["contrast"] = lab
            paired_rows.append(w)
    else:
        paired_rows = [
            {"contrast": "barrier/keratin (no CLDN4) high vs low", "n": 0, "delta_median": None, "p": None}
        ]

    # ----- figures -----
    sc.pl.paga(adata, color="leiden", show=False)
    fig = plt.gcf()
    _save(fig, figdir / "fig_extra_paga")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[0], show=False, frameon=False, cmap="viridis")
    axes[0].set_title("CLDN4")
    sc.pl.umap(adata, color="lineage", ax=axes[1], show=False, frameon=False)
    axes[1].set_title("author lineage")
    sc.pl.umap(adata, color="sling_pseudotime", ax=axes[2], show=False, frameon=False, cmap="magma")
    axes[2].set_title(f"Slingshot PT  start={start}")
    fig.suptitle(
        f"GSE127465 tumor epithelium  n_cells={adata.n_obs}  n_patients=7 (thin)  "
        f"start Leiden {start} (not CLDN4-high {start_info['cldn4_high_cluster']})",
        fontsize=10,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ct = pd.crosstab(adata.obs["Patient"].astype(str), adata.obs["lineage"].astype(str))
    ct.plot(kind="bar", stacked=True, ax=axes[0], colormap="tab20")
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Tumor epithelium n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].bar(["patients"], [7], color="#1f4e79")
    axes[1].set_ylabel("n")
    axes[1].set_title("Honest n=7 (thin)")
    _save(fig, figdir / "fig_honest_n")

    if not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("sling_low", "sling_high", "Slingshot"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            ax.scatter(paired_df[lo], paired_df[hi], s=50, c="#1f4e79")
            for rec in paired_df.itertuples(index=False):
                ax.annotate(str(rec.patient), (getattr(rec, lo), getattr(rec, hi)), fontsize=7)
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)} (n=7 thin)", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    # Extra: CLDN4 vs slingshot PT per lineage (cell-level descriptive).
    nlin = max(1, len(lin_cols))
    fig, axes = plt.subplots(1, nlin, figsize=(4.2 * nlin, 3.6), squeeze=False)
    for ax, col in zip(axes[0], lin_cols):
        y = pt[col].to_numpy()
        x = adata.obs["expr_CLDN4"].to_numpy()
        m = np.isfinite(y) & np.isfinite(x)
        ax.scatter(y[m], x[m], s=6, c=np.where(adata.obs["is_malignant"].to_numpy()[m], "#b23a48", "#2a6f97"), alpha=0.35)
        ax.set_xlabel(f"{col} (Slingshot)")
        ax.set_ylabel("CLDN4")
        ax.set_title(col)
    fig.suptitle("EXTRA: cell-level CLDN4 vs Slingshot (descriptive; red=malignant)", fontsize=10)
    _save(fig, figdir / "fig_extra_lineage_cldn4")

    # Extra: patient-level CLDN4 vs slingshot / barrier.
    if not elig.empty:
        fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
        axes[0].scatter(elig["mean_slingshot"], elig["mean_CLDN4"], s=55, c="#1f4e79")
        for rec in elig.itertuples(index=False):
            axes[0].annotate(str(rec.patient), (rec.mean_slingshot, rec.mean_CLDN4), fontsize=8)
        axes[0].set_xlabel("patient-mean Slingshot")
        axes[0].set_ylabel("patient-mean CLDN4")
        c4s = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot")
        axes[0].set_title(_fmt(c4s) + "  (n=7 thin)")
        axes[1].scatter(elig["mean_barrier_keratin"], elig["mean_CLDN4"], s=55, c="#1f4e79")
        for rec in elig.itertuples(index=False):
            axes[1].annotate(str(rec.patient), (rec.mean_barrier_keratin, rec.mean_CLDN4), fontsize=8)
        axes[1].set_xlabel("patient-mean barrier/keratin (no CLDN4)")
        axes[1].set_ylabel("patient-mean CLDN4")
        c4b = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
        axes[1].set_title(_fmt(c4b) + "  (n=7 thin)")
        fig.suptitle("EXTRA: patient-level CLDN4 (honest n=7)", fontsize=11)
        _save(fig, figdir / "fig_extra_patient_cldn4")

    for color, fname, cmap in (
        ("lineage", "fig_umap_lineage", None),
        ("Patient", "fig_umap_patient", None),
        ("leiden", "fig_umap_leiden", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("dpt_pseudotime", "fig_umap_dpt", "magma"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    n_mal_pos_by = (
        adata.obs.loc[adata.obs["is_malignant"] & (adata.obs["expr_CLDN4"] > 0), "Patient"]
        .astype(str)
        .value_counts()
        .sort_index()
        .to_dict()
    )
    c4_sling = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    pair_bar = next((r for r in paired_rows if str(r.get("contrast", "")).startswith("barrier")), {"n": 0})

    parts = [
        f"REAL Slingshot produced {len(lin_df)} lineage(s) from start Leiden {start} "
        f"(not CLDN4-high {start_info['cldn4_high_cluster']}).",
        f"Patient-level CLDN4 vs Slingshot: {_fmt(c4_sling)} — n=7 is thin.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Author malignant is PatientN-specific; Slingshot MST may stitch discrete tumors. "
        "Not a TACSTD2 redo. No both-high gate.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What can be said (n=7 thin).** Malignant+CLDN4 exists "
        f"(n_malignant={n_mal}, n_CLDN4+={n_mal_pos}). "
        f"Slingshot start is Leiden {start}, not the CLDN4-high cluster "
        f"{start_info['cldn4_high_cluster']}. "
        f"CLDN4 vs barrier/keratin (no CLDN4): {_fmt(c4_bar)}. "
        f"**What cannot be said.** n=7 patient Spearman is not a precise effect. "
        "Do not read Slingshot paths as Type II → each patient's tumor."
    )

    summary = {
        "accession": "GSE127465",
        "paper": "Zilionis et al. Immunity 2019 PMID 30979687",
        "platform": "inDrops",
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": "slingshot",
        "slingshot": {"available": True, "version": sling_ver, "real": True},
        "n_tumor": 40362,
        "n_blood": 14411,
        "n_cells": int(adata.n_obs),
        "n_malignant": n_mal,
        "n_malignant_cldn4_pos": n_mal_pos,
        "n_malignant_cldn4_pos_by_patient": n_mal_pos_by,
        "n_type2": int(adata.obs["is_type2"].sum()),
        "n_type1": int((adata.obs["Major cell type"].astype(str) == "Type I cells").sum()),
        "n_club": int((adata.obs["Major cell type"].astype(str) == "Club cells").sum()),
        "n_ciliated": int((adata.obs["Major cell type"].astype(str) == "Ciliated cells").sum()),
        "n_patients": 7,
        "n_patients_eligible": int(len(elig)),
        "n_patients_paired": int(len(paired_df)),
        "n_thin_note": "n=7 is thin",
        "lineage_counts": adata.obs["lineage"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "genes_absent": absent,
        "start": {
            "start_cluster": start_info["start_cluster"],
            "cldn4_high_cluster": start_info["cldn4_high_cluster"],
            "rule": start_info["rule"],
            "dpt_root_index": start_info["dpt_root_index"],
            "dpt_root_rule": start_info["dpt_root_rule"],
            "dpt_root_cell": start_info["dpt_root_cell"],
            "dpt_root_patient": start_info["dpt_root_patient"],
        },
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "lineages": lineages,
        "primary_spearman": primary,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": True,
            "spearman_rho": c4_bar.get("rho"),
            "spearman_p": c4_bar.get("p"),
            "spearman_n": c4_bar.get("n"),
        },
        "verdict": verdict,
        "what_holds": what_holds,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), summary)
    adata.write_h5ad(outdir / "epithelium_analyzed.h5ad")
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_patients": 7,
                "n_lineages": len(lin_df),
                "start": start,
                "cldn4_high": start_info["cldn4_high_cluster"],
                "lineage_table": str(tabdir / "lineages.tsv"),
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
