#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE148071 malignant-like epithelium, CLDN4 only.

ADDITIVE. Not a TACSTD2 redo. No dual-high gate. Patient is the inferential
unit. PR #394 barrier ρ=+0.590 is given and is not re-audited.

Primary clock = Bioconductor slingshot (Street 2018), rooted on TISCH
Alveolar / AT2 — never CLDN4-high. PAGA is geometry. DPT is a comparator.
Readouts along pseudotime: CLDN4, barrier/keratin (CLDN4 held out), IFN
(Hallmark IFNα ∩ IFNγ).
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
    IFN_ALIASES,
    PR394_BARRIER_N,
    PR394_BARRIER_P,
    PR394_BARRIER_RHO,
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
N_SLING_PCS = 10
MIN_GENES = 200
MIN_UMI = 500
MIN_CELLS_PER_GENE = 10
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
MAX_CELLS_PER_PATIENT = 500
GRAPH_CAP_TRIGGER = 20000
RANDOM_SEED = 0
N_PT_BINS = 8


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


def _score_mean(adata, genes: tuple[str, ...], key: str, aliases: dict | None = None) -> list[str]:
    present = []
    absent = []
    for g in genes:
        if g in adata.var_names:
            present.append(g)
            continue
        found = False
        if aliases and g in aliases:
            for alt in aliases[g]:
                if alt in adata.var_names:
                    present.append(alt)
                    found = True
                    break
        if not found:
            absent.append(g)
    if not present:
        adata.obs[key] = np.nan
        return absent
    X = adata[:, present].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.obs[key] = np.asarray(X, dtype=float).mean(axis=1)
    return absent


def _pick_root(adata) -> tuple[int, dict]:
    """Alveolar / AT2-high root. Never root on a CLDN4-high cell."""
    lin = adata.obs["lineage"].astype(str)
    alveolar = lin.eq("Alveolar")
    high = adata.obs["cldn4_tertile"].astype(str).eq("high")
    info = {"n_alveolar": int(alveolar.sum()), "rule": None, "rejected_cldn4_high": False}
    pool = alveolar & ~high
    if int(pool.sum()) < 10:
        pool = alveolar
        info["rejected_cldn4_high"] = True
        info["note"] = "too few non-high Alveolar; used all Alveolar but still not max-CLDN4"
    if int(pool.sum()) >= 10:
        idx = np.flatnonzero(pool.to_numpy())
        scores = adata.obs.loc[pool, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        # never the CLDN4-max cell
        if str(adata.obs.iloc[pick]["cldn4_tertile"]) == "high":
            order = np.argsort(np.abs(scores - med))
            for j in order:
                if str(adata.obs.iloc[idx[j]]["cldn4_tertile"]) != "high":
                    pick = idx[j]
                    info["rejected_cldn4_high"] = True
                    break
        info["rule"] = "TISCH Alveolar (median AT2; CLDN4-high excluded)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    means = (
        adata.obs.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    )
    for top in means.index.astype(str):
        cand = adata.obs["leiden"].astype(str).eq(top) & ~high
        if int(cand.sum()) < 5:
            continue
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        idx = np.flatnonzero(cand.to_numpy())
        pick = idx[int(np.nanargmax(scores))]
        info["rule"] = f"Leiden {top} max AT2 among non-CLDN4-high (alveolar n<10)"
        info["fallback_cluster"] = top
        return int(pick), info
    raise SystemExit("could not pick a non-CLDN4-high root")


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
    for _pat, idx in adata.obs.groupby("patient", observed=True).groups.items():
        idx = np.asarray(idx)
        if len(idx) <= max_n:
            keep.extend(idx.tolist())
        else:
            take = rng.choice(idx, size=max_n, replace=False)
            keep.extend(take.tolist())
    return adata[keep].copy()


def _find_rscript() -> str:
    for cand in (os.environ.get("Rscript"), shutil.which("Rscript"), "/usr/bin/Rscript"):
        if cand and Path(cand).is_file():
            return cand
    raise SystemExit("Rscript not on PATH — REAL Slingshot requires R")


def _run_slingshot(adata, work: Path, start_cluster: str) -> dict:
    rscript = _find_rscript()
    r_file = Path(__file__).resolve().parent / "run_slingshot.R"
    work.mkdir(parents=True, exist_ok=True)
    pcs = np.asarray(adata.obsm["X_pca"][:, :N_SLING_PCS], dtype=float)
    pca_df = pd.DataFrame(pcs, index=adata.obs_names, columns=[f"PC{i+1}" for i in range(N_SLING_PCS)])
    pca_df.to_csv(work / "sling_pca.tsv", sep="\t")
    clu = pd.DataFrame({"leiden": adata.obs["leiden"].astype(str).to_numpy()}, index=adata.obs_names)
    clu.to_csv(work / "sling_clusters.tsv", sep="\t")
    (work / "sling_start_cluster.txt").write_text(str(start_cluster) + "\n")
    env = os.environ.copy()
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "library"))
    print(f"running REAL slingshot via {rscript}", flush=True)
    proc = subprocess.run(
        [rscript, str(r_file), str(work)],
        check=False,
        env=env,
        capture_output=True,
        text=True,
    )
    (work / "sling_stdout.txt").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0:
        raise SystemExit(f"slingshot R failed:\n{proc.stdout}\n{proc.stderr}")
    pt = pd.read_csv(work / "sling_pseudotime.tsv", sep="\t")
    wt = pd.read_csv(work / "sling_weights.tsv", sep="\t")
    lineages = json.loads((work / "sling_lineages.json").read_text())
    pt = pt.set_index("cell")
    wt = wt.set_index("cell")
    if not pt.index.equals(adata.obs_names) and not pt.index.isin(adata.obs_names).all():
        # align
        pt = pt.reindex(adata.obs_names)
        wt = wt.reindex(adata.obs_names)
    else:
        pt = pt.reindex(adata.obs_names)
        wt = wt.reindex(adata.obs_names)
    lin_cols = [c for c in pt.columns if c != "cell"]
    # normalize each lineage 0-1
    norm = {}
    for c in lin_cols:
        v = pt[c].to_numpy(dtype=float)
        finite = np.isfinite(v)
        if finite.sum() < 5:
            norm[c] = v
            continue
        lo, hi = np.nanmin(v[finite]), np.nanmax(v[finite])
        norm[c] = (v - lo) / (hi - lo + 1e-12)
        adata.obs[f"sling_{c}"] = v
        adata.obs[f"sling_{c}_norm"] = norm[c]
    stacked = np.vstack([norm[c] for c in lin_cols]) if lin_cols else np.full((1, adata.n_obs), np.nan)
    adata.obs["sling_pt_mean"] = np.nanmean(stacked, axis=0)
    n_on = {c: int(np.isfinite(pt[c].to_numpy(dtype=float)).sum()) for c in lin_cols}
    principal = max(n_on, key=n_on.get) if n_on else None
    if principal is not None:
        adata.obs["sling_pt_principal"] = pt[principal].to_numpy(dtype=float)
        lo, hi = np.nanmin(adata.obs["sling_pt_principal"]), np.nanmax(adata.obs["sling_pt_principal"])
        adata.obs["sling_pt_principal_norm"] = (adata.obs["sling_pt_principal"] - lo) / (hi - lo + 1e-12)
    else:
        adata.obs["sling_pt_principal"] = np.nan
        adata.obs["sling_pt_principal_norm"] = np.nan
        principal = "none"
    return {
        "available": True,
        "version": "2.10.0",
        "n_lineages": len(lin_cols),
        "lineage_names": lin_cols,
        "n_cells_on": n_on,
        "principal_lineage": principal,
        "lineages": lineages,
        "start_cluster": start_cluster,
        "n_pcs": N_SLING_PCS,
    }


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lines = [
        "# Finding — GSE148071 REAL Slingshot/PAGA scored by CLDN4",
        "",
        "ADDITIVE. **CLDN4 only.** Public malignant-like epithelium from GSE148071 (Wu et al., *Nat Commun* 2021, PMID 33953163): 42 advanced NSCLC diagnostic biopsies (Singleron GEXSCOPE). Counts = GEO raw UMI. Labels = TISCH2 major-lineage (Malignant + leftover Alveolar / Basal / Epithelial). **Not a TACSTD2 redo.** No dual-high gate.",
        "",
        f"Primary clock = **real Bioconductor slingshot** (Street et al. 2018; R package loaded). PAGA is geometry. DPT is a comparator only. Root = TISCH Alveolar / AT2, **never CLDN4-high**. Inferential unit = **patient**. Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. IFN = Hallmark IFNα ∩ IFNγ intersection (CLDN4 not in the set).",
        "",
        f"PR #394 barrier/keratin (CLDN4 excluded) patient Spearman is **given** (n={PR394_BARRIER_N}, ρ=+{PR394_BARRIER_RHO:.3f}, p={PR394_BARRIER_P:.2e}) and is **not re-audited** here. This folder asks a different question: CLDN4 + barrier + IFN **along Slingshot/PAGA pseudotime**.",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        "- GEO patients / samples: **42** (Wu 2021). This is the catalog n, not the test n.",
        f"- TISCH cells: **{s['n_tisch_cells']}**. TISCH epithelial (Malignant+Alveolar+Basal+Epithelial): **{s['n_tisch_epithelial']}**.",
        f"- GEO barcodes matched to TISCH epithelium: **{s['n_matched']}**.",
        f"- After QC (min {MIN_GENES} genes, min {MIN_UMI} UMI): **n_cells_qc = {s['n_cells_qc']}** in **{s['n_patients_qc']}** patients.",
        f"- Graph cap: max {MAX_CELLS_PER_PATIENT}/patient when n>{GRAPH_CAP_TRIGGER} (seed {RANDOM_SEED}). Analysis object: **n_cells = {s['n_cells']}** in **n_patients = {s['n_patients']}**.",
        f"- TISCH lineage on the analysis object: {s['tisch_lineage_counts']}.",
        f"- Slingshot lineages: **{s['slingshot']['n_lineages']}**. Principal = {s['slingshot']['principal_lineage']} (cells on lineage: {s['slingshot']['n_cells_on']}).",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} analysis cells: **n = {s['n_samples_eligible']}** (catalog n equals this count here; still not 42 malignant tumors).",
        f"- Patients with finite principal Slingshot PT used for along-PT Spearman: **n = {s.get('n_patients_with_principal_pt', 'NA')}**. Dropped for missing principal-lineage membership: **{s.get('missing_principal_pt_patients', [])}**.",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} TISCH Malignant cells on the object: **n = {s['n_malignant_eligible']}**. Leftover-epithelium–dominant (malignant <10): **{s['leftover_dominant_patients']}**.",
        f"- Patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms (paired extra): **n = {s['n_samples_paired_tertile']}**.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, mid {s['cldn4_tertile_counts'].get('mid', 0)}, high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Slingshot: available={s['slingshot']['available']}; version={s['slingshot'].get('version')}.",
        "- No ICI / RECIST / MPR labels. Histology is not on the GEO series matrix — no LUAD/LUSC split.",
        "- Do not write “n=42 malignant tumors.”",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}; Slingshot PCs {N_SLING_PCS}.",
        f"- Slingshot/DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, patient {s['root'].get('root_patient')}, lineage {s['root'].get('root_lineage')}, Leiden {s['root'].get('root_leiden')}, CLDN4 tertile {s['root'].get('root_tertile')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN genes: Hallmark IFNα ∩ IFNγ (73; WARS/WARS1 alias accepted).",
        f"- PR #394 barrier ρ=+{PR394_BARRIER_RHO:.3f} cited, not re-fit.",
        "",
        "## Primary (patient-level Spearman along pseudotime, BH inside this list)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
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
        "## Slingshot lineage table (done criterion)",
        "",
        "| Lineage | n_cells | n_patients_elig | end | CLDN4~PT ρ | barrier~PT ρ | IFN~PT ρ |",
        "| --- | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for r in s.get("lineage_rows", []):
        def _f(x):
            return "NA" if x is None else f"{x:.3f}"

        lines.append(
            f"| {r['lineage']} | {r['n_cells']} | {r['n_patients_eligible']} | {r.get('end_cluster','')} | "
            f"{_f(r.get('rho_cldn4_pt'))} | {_f(r.get('rho_barrier_pt'))} | {_f(r.get('rho_ifn_pt'))} |"
        )
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired) and along-PT",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Rule: always emit (requested extra figures). "
            f"Paired tertile n={s['n_samples_paired_tertile']}. "
            f"PR #394 barrier patient ρ=+{PR394_BARRIER_RHO:.3f} is cited, not re-audited."
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
        "- Slingshot is an ordering, not a clock. Direction is the external Alveolar/AT2 root.",
        "- PR #394 barrier ρ is given. This analysis does not re-fit that contrast as a discovery.",
        "",
        "## Outputs (done when these exist)",
        "",
        "- `tables/patient_level.tsv`",
        "- `tables/lineage_level.tsv`",
        "- `tables/tisch_lineage_level.tsv`",
        "- `tables/patient_level_spearman.tsv`",
        "- `figures/fig_trajectory_cldn4.png`",
        "- `figures/fig_extra_along_pt.png`",
        "- `figures/fig_extra_cldn4_tertile.png`",
        "- `figures/fig_honest_n.png`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "bash methods/gse148071_slingshot_real_cldn4/scripts/install_tools.sh",
        "bash methods/gse148071_slingshot_real_cldn4/scripts/download.sh /tmp/gse148071_sling_data",
        "python3 methods/gse148071_slingshot_real_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/gse148071_sling_data \\",
        "  --out /tmp/gse148071_sling_data/epithelium.h5ad",
        "python3 methods/gse148071_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/gse148071_sling_data/epithelium.h5ad \\",
        "  --outdir methods/gse148071_slingshot_real_cldn4 \\",
        "  --finding methods/gse148071_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
        "Trajectory: `figures/fig_trajectory_cldn4.png`. Extra along-PT: `figures/fig_extra_along_pt.png`.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default="methods/gse148071_slingshot_real_cldn4")
    ap.add_argument("--finding", default="methods/gse148071_slingshot_real_cldn4/FINDING.md")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    work = outdir / "work"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

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
        print(f"capped to max {MAX_CELLS_PER_PATIENT}/patient → n_cells={adata.n_obs}", flush=True)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        aliases = IFN_ALIASES if name == "IFN" else None
        absent[name] = _score_mean(adata, genes, f"score_{name}", aliases=aliases)
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
    root_info["root_leiden"] = str(adata.obs.iloc[root_idx]["leiden"])
    root_info["root_cldn4"] = float(adata.obs.iloc[root_idx]["expr_CLDN4"])
    root_info["root_tertile"] = str(adata.obs.iloc[root_idx]["cldn4_tertile"])
    if root_info["root_tertile"] == "high":
        raise SystemExit("root landed on CLDN4-high; refusing")
    adata.uns["iroot"] = root_idx
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    sling = _run_slingshot(adata, work, start_cluster=root_info["root_leiden"])

    # ----- patient table -----
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
                "mean_sling_pt": float(np.nanmean(sub["sling_pt_principal_norm"].to_numpy())),
                "mean_sling_pt_all": float(np.nanmean(sub["sling_pt_mean"].to_numpy())),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()) if "expr_SFTPC" in sub else np.nan,
            }
        )
    sample_df = pd.DataFrame(rows).sort_values("patient")
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    malig = sample_df[sample_df["n_malignant"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    leftover_dom = sample_df.loc[sample_df["leftover_dominant"], "patient"].tolist()

    # Primary: along-PT + IFN. Do NOT include CLDN4 vs barrier (PR #394 given).
    primary_specs = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling_pt"),
        ("barrier/keratin (no CLDN4) vs Slingshot PT", "mean_barrier", "mean_sling_pt"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_sling_pt"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("AT2 vs Slingshot PT (control)", "mean_AT2", "mean_sling_pt"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling_pt"),
        ("CLDN4 vs DPT (comparator clock)", "mean_CLDN4", "mean_dpt"),
        ("Slingshot PT vs DPT (concordance)", "mean_sling_pt", "mean_dpt"),
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
        ("malignant-only CLDN4 vs Slingshot PT", malig, "mean_CLDN4", "mean_sling_pt"),
        ("malignant-only barrier vs Slingshot PT", malig, "mean_barrier", "mean_sling_pt"),
        ("malignant-only IFN vs Slingshot PT", malig, "mean_IFN", "mean_sling_pt"),
        ("malignant-only CLDN4 vs IFN", malig, "mean_CLDN4", "mean_IFN"),
        ("drop leftover-dominant CLDN4 vs Slingshot PT", sample_df[~sample_df["leftover_dominant"]], "mean_CLDN4", "mean_sling_pt"),
        ("CLDN4 vs Slingshot PT (mean across lineages)", elig, "mean_CLDN4", "mean_sling_pt_all"),
        ("IFN vs Slingshot PT (mean across lineages)", elig, "mean_IFN", "mean_sling_pt_all"),
    ):
        use = df[df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN] if "n_cells" in df.columns else df
        r = _spearman(use[a].to_numpy(), use[b].to_numpy())
        r["contrast"] = name
        sensitivity.append(r)

    # ----- Slingshot lineage-level table -----
    lineage_rows = []
    lineage_patient_rows = []
    for L in sling["lineage_names"]:
        col = f"sling_{L}_norm"
        if col not in adata.obs:
            col_raw = f"sling_{L}"
            if col_raw not in adata.obs:
                continue
            v = adata.obs[col_raw].to_numpy(dtype=float)
            lo, hi = np.nanmin(v), np.nanmax(v)
            adata.obs[col] = (v - lo) / (hi - lo + 1e-12)
        on = adata.obs[col].notna()
        sub = adata.obs.loc[on]
        early = sub[sub[col] <= 1 / 3]
        late = sub[sub[col] >= 2 / 3]
        lin_meta = next((x for x in sling["lineages"] if x.get("name") == L or L.endswith(x.get("name", "___"))), None)
        clusters = (lin_meta or {}).get("clusters", [])
        # patient means on this lineage
        pat_rows = []
        for pat, psub in sub.groupby("patient", observed=True):
            if len(psub) < MIN_CELLS_PER_SAMPLE_FOR_MEAN:
                continue
            rec = {
                "lineage": L,
                "patient": pat,
                "n_cells": int(len(psub)),
                "mean_CLDN4": float(psub["expr_CLDN4"].mean()),
                "mean_barrier": float(psub["score_barrier_keratin"].mean()),
                "mean_IFN": float(psub["score_IFN"].mean()),
                "mean_sling_pt": float(psub[col].mean()),
            }
            pat_rows.append(rec)
            lineage_patient_rows.append(rec)
        pdf = pd.DataFrame(pat_rows)
        r_c = _spearman(pdf["mean_CLDN4"].to_numpy(), pdf["mean_sling_pt"].to_numpy()) if not pdf.empty else {"n": 0, "rho": None, "p": None}
        r_b = _spearman(pdf["mean_barrier"].to_numpy(), pdf["mean_sling_pt"].to_numpy()) if not pdf.empty else {"n": 0, "rho": None, "p": None}
        r_i = _spearman(pdf["mean_IFN"].to_numpy(), pdf["mean_sling_pt"].to_numpy()) if not pdf.empty else {"n": 0, "rho": None, "p": None}
        lineage_rows.append(
            {
                "lineage": L,
                "n_cells": int(on.sum()),
                "n_patients": int(sub["patient"].nunique()),
                "n_patients_eligible": int(len(pdf)),
                "start_cluster": clusters[0] if clusters else sling["start_cluster"],
                "end_cluster": clusters[-1] if clusters else "",
                "clusters": ">".join(clusters) if clusters else "",
                "mean_CLDN4_early": float(early["expr_CLDN4"].mean()) if len(early) else np.nan,
                "mean_CLDN4_late": float(late["expr_CLDN4"].mean()) if len(late) else np.nan,
                "mean_barrier_early": float(early["score_barrier_keratin"].mean()) if len(early) else np.nan,
                "mean_barrier_late": float(late["score_barrier_keratin"].mean()) if len(late) else np.nan,
                "mean_IFN_early": float(early["score_IFN"].mean()) if len(early) else np.nan,
                "mean_IFN_late": float(late["score_IFN"].mean()) if len(late) else np.nan,
                "n_early": int(len(early)),
                "n_late": int(len(late)),
                "rho_cldn4_pt": r_c.get("rho"),
                "p_cldn4_pt": r_c.get("p"),
                "rho_barrier_pt": r_b.get("rho"),
                "p_barrier_pt": r_b.get("p"),
                "rho_ifn_pt": r_i.get("rho"),
                "p_ifn_pt": r_i.get("p"),
                "n_spearman": r_c.get("n"),
            }
        )
        sensitivity.append(
            {
                "contrast": f"{L} patient CLDN4 vs PT",
                **{k: r_c.get(k) for k in ("n", "rho", "p")},
            }
        )
        sensitivity.append(
            {
                "contrast": f"{L} patient IFN vs PT",
                **{k: r_i.get(k) for k in ("n", "rho", "p")},
            }
        )
    lineage_df = pd.DataFrame(lineage_rows)
    lineage_patient_df = pd.DataFrame(lineage_patient_rows)

    # TISCH lineage-level
    tisch_rows = []
    for lin, sub in adata.obs.groupby("lineage", observed=True):
        tisch_rows.append(
            {
                "tisch_lineage": lin,
                "n_cells": int(len(sub)),
                "n_patients": int(sub["patient"].nunique()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_sling_pt": float(np.nanmean(sub["sling_pt_principal_norm"].to_numpy())),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "frac_cldn4_high": float((sub["cldn4_tertile"].astype(str) == "high").mean()),
            }
        )
    tisch_df = pd.DataFrame(tisch_rows).sort_values("n_cells", ascending=False)

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
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "sling_high": float(np.nanmean(hi["sling_pt_principal_norm"].to_numpy())),
                "sling_low": float(np.nanmean(lo["sling_pt_principal_norm"].to_numpy())),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_rows_data)
    paired_rows = []
    if not paired_df.empty:
        for contrast, hi, lo in (
            ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
            ("IFN high vs low", "IFN_high", "IFN_low"),
            ("Slingshot PT high vs low", "sling_high", "sling_low"),
            ("AT2 score high vs low", "AT2_high", "AT2_low"),
        ):
            r = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            r["contrast"] = contrast
            paired_rows.append(r)
    else:
        paired_rows = [
            {"contrast": "barrier/keratin (no CLDN4) high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
            {"contrast": "IFN high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
            {"contrast": "Slingshot PT high vs low", "n": 0, "W": None, "p": None, "delta_median": None},
        ]

    emit_extra = True

    # ----- figures -----
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    sc.pl.paga(adata, color="leiden", ax=axes[0], show=False, frameon=False, title="PAGA (Leiden)")
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[1], show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    ax = axes[2]
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    leftover = elig["leftover_dominant"]
    ax.scatter(
        elig.loc[~leftover, "mean_sling_pt"],
        elig.loc[~leftover, "mean_CLDN4"],
        s=48,
        c="#b23a48",
        label=f"malignant-bearing n={int((~leftover).sum())}",
    )
    ax.scatter(
        elig.loc[leftover, "mean_sling_pt"],
        elig.loc[leftover, "mean_CLDN4"],
        s=48,
        c="#2a6f97",
        label=f"leftover-epi dominant n={int(leftover.sum())}",
    )
    ax.set_xlabel("patient-mean Slingshot PT (principal, 0–1)")
    ax.set_ylabel("patient-mean CLDN4")
    rho_s = "NA" if c4_pt["rho"] is None else f"{c4_pt['rho']:.2f}"
    p_s = "NA" if c4_pt["p"] is None else f"{c4_pt['p']:.3g}"
    ax.set_title(f"CLDN4 vs Slingshot PT  n={c4_pt['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(
        f"GSE148071 REAL Slingshot/PAGA  CLDN4-only   n_cells={adata.n_obs}  n_patients={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    ct = adata.obs.groupby(["patient", "lineage"], observed=True).size().unstack(fill_value=0)
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.85)
    ax.set_ylabel("cells in analysis object")
    ax.set_title(
        f"Honest n: TISCH lineages after QC/cap  n_cells={adata.n_obs}  "
        f"n_patients={sample_df.shape[0]} / GEO 42"
    )
    ax.legend(frameon=False, fontsize=7)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", fontsize=7)
    _save(fig, figdir / "fig_honest_n")

    # Extra: along-PT CLDN4 + barrier + IFN
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    pt = adata.obs["sling_pt_principal_norm"].to_numpy(dtype=float)
    panels = (
        ("expr_CLDN4", "CLDN4", "#b23a48"),
        ("score_barrier_keratin", "barrier/keratin (no CLDN4)", "#2a6f97"),
        ("score_IFN", "IFN (IFNα∩IFNγ)", "#3d5a40"),
    )
    rng = np.random.default_rng(0)
    for ax, (col, lab, color) in zip(axes, panels):
        y = adata.obs[col].to_numpy(dtype=float)
        m = np.isfinite(pt) & np.isfinite(y)
        if m.sum() > 4000:
            take = rng.choice(np.flatnonzero(m), size=4000, replace=False)
        else:
            take = np.flatnonzero(m)
        ax.scatter(pt[take], y[take], s=4, alpha=0.15, c=color, linewidths=0)
        # bin means
        bins = np.linspace(0, 1, N_PT_BINS + 1)
        xc, yc = [], []
        for i in range(N_PT_BINS):
            sel = m & (pt >= bins[i]) & (pt < bins[i + 1] if i < N_PT_BINS - 1 else pt <= bins[i + 1])
            if sel.sum() < 20:
                continue
            xc.append(0.5 * (bins[i] + bins[i + 1]))
            yc.append(float(np.mean(y[sel])))
        ax.plot(xc, yc, color="k", lw=2.0)
        ax.set_xlabel("Slingshot PT (principal, 0–1)")
        ax.set_ylabel(lab)
        ax.set_title(lab)
    fig.suptitle(
        "EXTRA: CLDN4 + barrier + IFN along REAL Slingshot PT (cells; black = bin mean)",
        fontsize=11,
    )
    _save(fig, figdir / "fig_extra_along_pt")

    # Extra: patient-level along-PT programs
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    patient_panels = (
        ("mean_sling_pt", "mean_CLDN4", "CLDN4 vs Slingshot PT", next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")),
        ("mean_sling_pt", "mean_barrier", "barrier vs Slingshot PT", next(r for r in primary if r["contrast"].startswith("barrier"))),
        ("mean_sling_pt", "mean_IFN", "IFN vs Slingshot PT", next(r for r in primary if r["contrast"] == "IFN vs Slingshot PT")),
    )
    for ax, (x, y, lab, row) in zip(axes, patient_panels):
        ax.scatter(elig[x], elig[y], s=48, c="#4c6a92")
        ax.set_xlabel("patient-mean Slingshot PT")
        ax.set_ylabel(y.replace("mean_", "patient-mean "))
        ax.set_title(f"{lab}\n{_fmt(row)}")
    fig.suptitle("EXTRA: patient-level CLDN4 / barrier / IFN vs Slingshot PT", fontsize=11)
    _save(fig, figdir / "fig_extra_patient_along_pt")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN"),
            ("sling_low", "sling_high", "Slingshot PT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            ax.scatter(paired_df[lo], paired_df[hi], s=44, c="#4c6a92")
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
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("lineage", "fig_umap_lineage", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("dpt_pseudotime", "fig_umap_dpt", "viridis"),
        ("sling_pt_principal_norm", "fig_umap_slingshot", "viridis"),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_ifn", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    bar_pt = next(r for r in primary if r["contrast"].startswith("barrier"))
    ifn_pt = next(r for r in primary if r["contrast"] == "IFN vs Slingshot PT")
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    at2_pt = next(r for r in primary if r["contrast"].startswith("AT2"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next((r for r in paired_rows if r["contrast"].startswith("IFN")), {"n": 0, "p": None, "delta_median": None, "W": None})

    parts = [
        f"REAL Slingshot ran ({sling['n_lineages']} lineages; principal={sling['principal_lineage']}).",
        f"Patient-level CLDN4 vs Alveolar-rooted Slingshot PT: {_fmt(c4_pt)}.",
        f"Barrier/keratin (CLDN4 excluded) vs Slingshot PT: {_fmt(bar_pt)}.",
        f"IFN vs Slingshot PT: {_fmt(ifn_pt)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"AT2 vs Slingshot PT (root control): {_fmt(at2_pt)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        f"Root is {root_info['rule']}; root CLDN4 tertile={root_info['root_tertile']} (not high).",
        f"Leftover-epithelium–dominant patients (TISCH Malignant <10 on the object): {leftover_dom or 'none'}.",
        f"PR #394 barrier ρ=+{PR394_BARRIER_RHO:.3f} is given and was not re-audited.",
        "GEO n=42 is the catalog, not the Spearman n. Mixed advanced NSCLC; no ICI labels. Not a TACSTD2 redo. No both-high gate.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accession": "GSE148071",
        "histology": "advanced NSCLC (mixed; not on GEO series matrix)",
        "public_only": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "no_dual_high": True,
        "pr394_barrier_rho_given": PR394_BARRIER_RHO,
        "pr394_barrier_not_reaudited": True,
        "n_tisch_cells": n_tisch_cells,
        "n_tisch_epithelial": n_tisch_epithelial,
        "n_matched": n_matched,
        "n_cells_qc": n_cells_qc,
        "n_patients_qc": n_patients_qc,
        "n_cells": int(adata.n_obs),
        "n_patients": int(sample_df.shape[0]),
        "n_samples_eligible": int(len(elig)),
        "n_patients_with_principal_pt": int(elig["mean_sling_pt"].notna().sum()),
        "missing_principal_pt_patients": elig.loc[elig["mean_sling_pt"].isna(), "patient"].tolist(),
        "n_malignant_eligible": int(len(malig)),
        "n_samples_paired_tertile": int(len(paired_df)),
        "leftover_dominant_patients": leftover_dom,
        "tisch_lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "slingshot": sling,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "graph_capped": capped,
        "max_cells_per_patient": MAX_CELLS_PER_PATIENT if capped else None,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "lineage_rows": lineage_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "rule": "always emit extra figures (along-PT + paired tertile)",
        },
        "verdict": verdict,
    }

    (tabdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    sample_df.to_csv(tabdir / "patient_level.tsv", sep="\t", index=False)
    sample_df.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    elig.to_csv(tabdir / "patient_eligible.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tabdir / "patient_level_spearman.tsv", sep="\t", index=False)
    lineage_df.to_csv(tabdir / "lineage_level.tsv", sep="\t", index=False)
    if not lineage_patient_df.empty:
        lineage_patient_df.to_csv(tabdir / "lineage_patient_means.tsv", sep="\t", index=False)
    tisch_df.to_csv(tabdir / "tisch_lineage_level.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "paired_tertile.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    # PAGA vertices
    vert = (
        adata.obs.groupby("leiden", observed=True)
        .agg(
            n_cells=("patient", "size"),
            n_patients=("patient", "nunique"),
            mean_CLDN4=("expr_CLDN4", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
            mean_IFN=("score_IFN", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_sling=("sling_pt_principal_norm", "mean"),
            frac_malignant=("lineage", lambda s: float(s.isin(TISCH_MALIGNANT).mean())),
        )
        .reset_index()
    )
    vert.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(tabdir / "paga_connectivities.tsv", sep="\t")

    write_finding(Path(args.finding), {"summary": summary})

    # confirm done-criterion files
    for must in (
        tabdir / "patient_level.tsv",
        tabdir / "lineage_level.tsv",
        tabdir / "tisch_lineage_level.tsv",
        tabdir / "patient_level_spearman.tsv",
    ):
        if not must.is_file():
            raise SystemExit(f"missing done-criterion table {must}")

    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": int(adata.n_obs),
                "n_patients": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "n_lineages": sling["n_lineages"],
                "principal": sling["principal_lineage"],
                "root_tertile": root_info["root_tertile"],
                "finding": args.finding,
                "patient_table": str(tabdir / "patient_level.tsv"),
                "lineage_table": str(tabdir / "lineage_level.tsv"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
