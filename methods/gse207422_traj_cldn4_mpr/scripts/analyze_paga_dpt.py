#!/usr/bin/env python3
"""PAGA + DPT on GSE207422 epithelium, scored by CLDN4 vs MPR.

ADDITIVE CLDN4-only. Not a TACSTD2 redo. No dual-high gate.
Does not re-audit the given GSE207422-only T/NK-flat result.

Question: among A3-malignant-like epithelial cells, does CLDN4 sit on a
trajectory that also separates MPR vs NMPR? Patient is the unit.

Graph is built on marker epithelium (needed for a leftover/AT2-like root
that is not CLDN4-high). Primary means are A3-malignant-like cells only.
R/Slingshot is not present — PAGA + DPT only.
"""

from __future__ import annotations

import argparse
import json
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
MIN_GENES = 200
MIN_UMI = 500
MIN_CELLS_PER_GENE = 10
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
MAX_CELLS_PER_PATIENT = 500
GRAPH_CAP_TRIGGER = 8000
RANDOM_SEED = 0
COLOR = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}


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
    return {"n": n, "rho": float(rho), "p": float(p), "note": ""}


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


def exact_mwu(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "n_nmpr": n_a,
            "n_mpr": n_b,
            "mean_nmpr": float(np.mean(a)) if n_a else None,
            "mean_mpr": float(np.mean(b)) if n_b else None,
            "median_nmpr": float(np.median(a)) if n_a else None,
            "median_mpr": float(np.median(b)) if n_b else None,
            "delta_mean_nmpr_minus_mpr": None,
            "mwu_u": None,
            "exact_p": None,
            "note": "too_few_samples",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    u_obs, p_exact = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    return {
        "n_nmpr": n_a,
        "n_mpr": n_b,
        "mean_nmpr": float(np.mean(a)),
        "mean_mpr": float(np.mean(b)),
        "median_nmpr": float(np.median(a)),
        "median_mpr": float(np.median(b)),
        "delta_mean_nmpr_minus_mpr": float(np.mean(a) - np.mean(b)),
        "mwu_u": float(u_obs),
        "exact_p": float(p_exact),
        "note": "",
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
    """Leftover / AT2-high root. Never root on a CLDN4-high cell."""
    leftover = adata.obs["is_leftover_epi"].astype(bool)
    not_high = adata.obs["cldn4_tertile"].astype(str) != "high"
    info = {
        "n_leftover": int(leftover.sum()),
        "n_leftover_not_cldn4_high": int((leftover & not_high).sum()),
        "rule": None,
    }
    cand = leftover & not_high
    if int(cand.sum()) >= 20:
        idx = np.flatnonzero(cand.to_numpy())
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        info["rule"] = "leftover epithelium, not CLDN4-high, median AT2 score"
        return int(pick), info
    if int(leftover.sum()) >= 10:
        idx = np.flatnonzero(leftover.to_numpy())
        scores = adata.obs.loc[leftover, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmax(scores))]
        info["rule"] = "leftover epithelium max AT2 (CLDN4-high filter too thin)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    means = (
        adata.obs.groupby("leiden", observed=True)["score_AT2"]
        .mean()
        .sort_values(ascending=False)
    )
    for top in means.index.astype(str):
        cluster = adata.obs["leiden"].astype(str) == top
        pool = cluster & not_high
        if int(pool.sum()) == 0:
            continue
        scores = adata.obs.loc[pool, "score_AT2"].to_numpy()
        idx = np.flatnonzero(pool.to_numpy())
        pick = idx[int(np.nanargmax(scores))]
        info["rule"] = f"Leiden {top} max AT2 among non-CLDN4-high (leftover n<10)"
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


def _fmt_rho(row: dict) -> str:
    n = row.get("n")
    if row.get("rho") is None:
        return f"n={n}, ρ=NA, p=NA"
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _fmt_mwu(row: dict) -> str:
    if row.get("note") == "too_few_samples" or row.get("exact_p") is None:
        return (
            f"n={row.get('n_nmpr')} vs {row.get('n_mpr')} (too few for exact MWU)"
        )
    return (
        f"mean {row['mean_nmpr']:.3f} vs {row['mean_mpr']:.3f} "
        f"(Δ={row['delta_mean_nmpr_minus_mpr']:+.3f}); "
        f"exact p={row['exact_p']:.3g}; n={row['n_nmpr']} vs {row['n_mpr']}"
    )


def _fmt_pair(row: dict) -> str:
    n = row.get("n")
    if row.get("W") is None:
        return f"n={n}, W=NA, p=NA"
    return f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, p={row['p']:.3g}"


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


def _patient_block(sub: pd.DataFrame, prefix: str, mask: pd.Series) -> dict:
    use = sub.loc[mask]
    n = int(len(use))
    out = {f"{prefix}_n": n}
    if n == 0:
        for k in ("mean_CLDN4", "mean_dpt", "mean_AT2", "mean_barrier", "mean_TACSTD2"):
            out[f"{prefix}_{k}"] = np.nan
        return out
    out[f"{prefix}_mean_CLDN4"] = float(use["expr_CLDN4"].mean())
    out[f"{prefix}_mean_dpt"] = float(use["dpt_pseudotime"].mean())
    out[f"{prefix}_mean_AT2"] = float(use["score_AT2"].mean())
    out[f"{prefix}_mean_barrier"] = float(use["score_barrier_keratin"].mean())
    out[f"{prefix}_mean_TACSTD2"] = (
        float(use["expr_TACSTD2"].mean()) if "expr_TACSTD2" in use else np.nan
    )
    return out


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]

    def rho_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    def mwu_md(r: dict) -> str:
        p = "NA" if r.get("exact_p") is None else f"{r['exact_p']:.3g}"
        d = (
            "NA"
            if r.get("delta_mean_nmpr_minus_mpr") is None
            else f"{r['delta_mean_nmpr_minus_mpr']:+.3f}"
        )
        mn = (
            "NA"
            if r.get("mean_nmpr") is None
            else f"{r['mean_nmpr']:.3f}"
        )
        mm = "NA" if r.get("mean_mpr") is None else f"{r['mean_mpr']:.3f}"
        note = r.get("note") or ""
        return (
            f"| {r['contrast']} | {r['n_nmpr']} | {r['n_mpr']} | {mn} | {mm} | "
            f"{d} | {p} | {note} |"
        )

    empty = ", ".join(s["empty_malig_post"]) or "none"
    noisy = ", ".join(s["noisy_malig_post"]) or "none"
    lines = [
        "# Finding — GSE207422 PAGA / DPT scored by CLDN4 vs MPR",
        "",
        "Additive public slice. **CLDN4-only. Not a TACSTD2 redo. No dual-high.** "
        "Primary readout is **CLDN4** on marker-defined A3-malignant-like epithelium "
        "from GSE207422 (Hu et al., *Genome Med* 2023, PMID 36915183): neoadjuvant "
        "PD-1 + platinum, MPR vs NMPR. Counts are the public GEO processed UMI matrix. "
        "Author CopyKAT / epithelium barcodes are **not** on GEO. "
        "GSE207422-only T/NK was already flat and is **not** re-audited here.",
        "",
        "PAGA + leftover/AT2-rooted diffusion pseudotime. Inferential unit = **patient**. "
        "Cell-level ρ is descriptive and is not the claim. "
        "Barrier/keratin score **excludes CLDN4**. Root is leftover (normal-lung-program) "
        "epithelium, never a CLDN4-high cell. R was not present; Slingshot was not run.",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        "- GEO samples: **15** (3 pre-treatment biopsies + 12 post-treatment resections). "
        "Catalog n is not the test n.",
        f"- Public UMI matrix: **{s['n_cells_matrix']}** cells × **{s['n_genes_matrix']}** genes.",
        f"- Marker epithelial cells: **{s['n_epithelial_all']}**. "
        f"A3-malignant-like: **{s['n_malig_a3_all']}**. "
        f"Leftover epithelium: **{s['n_leftover_all']}**.",
        f"- After QC (min {MIN_GENES} genes, min {MIN_UMI} UMI): "
        f"**n_cells_qc = {s['n_cells_qc']}** in **{s['n_patients_qc']}** samples.",
        f"- Graph cap: max {MAX_CELLS_PER_PATIENT}/sample when n>{GRAPH_CAP_TRIGGER} "
        f"(seed {RANDOM_SEED}). Analysis object: **n_cells = {s['n_cells']}** "
        f"in **{s['n_samples_on_graph']}** samples. Do not quote this as the inferential n.",
        f"- Lineage/compartment on the analysis object: {s['compartment_counts']}.",
        "- Inferential cohort: **12 post-treatment patients** (MPR n=4 including pCR P06; "
        "NMPR n=8). Pre-treatment biopsies stay on the graph for the manifold/root and "
        "are dropped from every test.",
        f"- Post patients with 0 A3-malignant cells on the analysis object "
        f"(NaN, dropped from malignant-like tests): **{empty}**.",
        f"- Post patients with 1–9 A3-malignant cells (kept only in the noisy sensitivity): "
        f"**{noisy}**.",
        f"- Post patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} A3-malignant cells "
        f"(primary malignant-like Spearman): **n = {s['n_malig_eligible_post']}**.",
        f"- Post patients with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells "
        f"(complete-case epithelium companion): **n = {s['n_epi_eligible_post']}**.",
        f"- Post patients with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and "
        f"CLDN4-low arms (paired extra): **n = {s['n_samples_paired_tertile']}**.",
        f"- CLDN4 tertile cells on the graph: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- This set is small. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. "
        "Malignant-like MPR residual tumors are often empty. Do not inflate cell count.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, "
        f"patient {s['root'].get('root_patient')}, compartment {s['root'].get('root_compartment')}, "
        f"CLDN4={s['root'].get('root_cldn4')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among "
        f"{s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
        "(same rule as the given A3 slice; **not** CopyKAT).",
        "- Slingshot: not run (R absent).",
        "",
        "## Primary — patient-level CLDN4 vs DPT (A3-malignant-like, post, ≥10 cells)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_cldn4_vs_dpt"]:
        lines.append(rho_md(r))
    lines += [
        "",
        "## Primary — patient-level DPT vs MPR (exact two-sided MWU)",
        "",
        "| Contrast | n_NMPR | n_MPR | mean NMPR | mean MPR | Δ (NMPR−MPR) | exact p | note |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for r in s["primary_dpt_vs_mpr"]:
        lines.append(mwu_md(r))
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
        "## Extra figure — CLDN4-high vs CLDN4-low (sample-paired, epithelium on graph)",
        "",
        (
            f"Emitted: **True**. Rule: always emit (requested extra figure). "
            f"Paired tertile n={s['n_samples_paired_tertile']}."
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
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- A3-malignant-like is a marker rule, **not CNV**. Residual unmarked epithelium can leak in.",
        "- Empty MPR residuals (normal-lung program) are reported, not patched by dual-high.",
        "- GSE207422-only T/NK was flat in the given slice and is not re-tested here.",
        "- No TACSTD2∩CLDN4 both-high gate. TACSTD2 is a companion Spearman only.",
        "- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Do not write n=12 as the malignant-like test n when empty MPR were dropped.",
        "- Mixed NSCLC (LUAD + LUSC). Not LUAD-only.",
        "",
        "## Outputs",
        "",
        "- `figures/fig_trajectory_cldn4.png` — PAGA + UMAP CLDN4 + patient CLDN4 vs DPT",
        "- `figures/fig_dpt_vs_mpr.png` — patient DPT vs MPR (malignant-like and epithelium)",
        "- `figures/fig_honest_n.png` — per-patient malignant vs leftover counts",
        "- `figures/fig_extra_cldn4_tertile.png` — within-sample CLDN4-high vs low",
        "- `tables/patient_means.tsv` — per-patient means (do not sum cells across patients)",
        "- `tables/patient_cldn4_vs_dpt.tsv` — patient-level CLDN4 vs DPT",
        "- `tables/patient_dpt_vs_mpr.tsv` — patient-level DPT vs MPR",
        "- `tables/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/gse207422_traj_cldn4_mpr/requirements.txt",
        "python3 methods/gse207422_traj_cldn4_mpr/scripts/download.py",
        "python3 methods/gse207422_traj_cldn4_mpr/scripts/extract_epithelium.py",
        "python3 methods/gse207422_traj_cldn4_mpr/scripts/analyze_paga_dpt.py",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="/tmp/gse207422_traj_data/epithelium.h5ad")
    ap.add_argument("--outdir", default="methods/gse207422_traj_cldn4_mpr")
    ap.add_argument("--finding", default="methods/gse207422_traj_cldn4_mpr/FINDING.md")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    np.random.seed(RANDOM_SEED)

    adata = sc.read_h5ad(args.input)
    extract = dict(adata.uns.get("extract", {}))
    n_cells_matrix = int(extract.get("n_cells_matrix", 92330))
    n_genes_matrix = int(extract.get("n_genes_matrix", adata.n_vars))
    n_epithelial_all = int(extract.get("n_epithelial", adata.n_obs))
    n_malig_a3_all = int(extract.get("n_malig_a3", int(adata.obs["is_malig_a3"].sum())))
    n_leftover_all = int(extract.get("n_leftover_epi", int(adata.obs["is_leftover_epi"].sum())))

    adata.obs["patient"] = adata.obs["patient"].astype(str)
    adata.obs["paper_group"] = adata.obs["paper_group"].astype(str)
    if "n_umi" not in adata.obs:
        adata.obs["n_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
    if "n_genes" not in adata.obs:
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["n_genes"] >= MIN_GENES].copy()
    adata = adata[adata.obs["n_umi"] >= MIN_UMI].copy()
    n_cells_qc = int(adata.n_obs)
    n_patients_qc = int(adata.obs["patient"].nunique())
    print(f"after QC n_cells={n_cells_qc} n_samples={n_patients_qc}", flush=True)

    capped = False
    if adata.n_obs > GRAPH_CAP_TRIGGER:
        adata = _cap_per_patient(adata, MAX_CELLS_PER_PATIENT, RANDOM_SEED)
        capped = True
        print(
            f"capped to max {MAX_CELLS_PER_PATIENT}/sample → n_cells={adata.n_obs}",
            flush=True,
        )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    absent: dict[str, list[str]] = {}
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
    sc.tl.leiden(
        adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False
    )
    sc.tl.umap(adata)
    sc.tl.paga(adata, groups="leiden")
    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = sorted(
        adata.obs["leiden"].astype(str).unique(),
        key=lambda x: int(x) if x.isdigit() else x,
    )

    root_idx, root_info = _pick_root(adata)
    root_obs = adata.obs.iloc[root_idx]
    root_info["index"] = root_idx
    root_info["root_patient"] = str(root_obs["patient"])
    root_info["root_compartment"] = str(root_obs["compartment"])
    root_info["root_cldn4"] = float(root_obs["expr_CLDN4"])
    root_info["root_tertile"] = str(root_obs["cldn4_tertile"])
    root_info["root_group"] = str(root_obs["paper_group"])
    if str(root_obs["cldn4_tertile"]) == "high":
        raise SystemExit("root is CLDN4-high; refuse circular DPT")
    adata.uns["iroot"] = root_idx
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    rows = []
    for pat, sub in adata.obs.groupby("patient", observed=True):
        group = str(sub["paper_group"].iloc[0])
        timing = "pre" if group == "TN" else "post"
        path = (
            str(sub["path_response"].iloc[0])
            if "path_response" in sub.columns
            else group
        )
        histo = (
            str(sub["Pathology"].iloc[0]) if "Pathology" in sub.columns else ""
        )
        recist = str(sub["RECIST"].iloc[0]) if "RECIST" in sub.columns else ""
        block = {
            "patient": pat,
            "Sample": str(sub["Sample"].iloc[0]),
            "paper_group": group,
            "path_response": path,
            "timing": timing,
            "Pathology": histo,
            "RECIST": recist,
            "n_cells": int(len(sub)),
            "n_malignant": int(sub["is_malig_a3"].sum()),
            "n_leftover_epi": int(sub["is_leftover_epi"].sum()),
        }
        block.update(_patient_block(sub, "mal", sub["is_malig_a3"].astype(bool)))
        block.update(_patient_block(sub, "epi", pd.Series(True, index=sub.index)))
        block.update(
            _patient_block(sub, "leftover", sub["is_leftover_epi"].astype(bool))
        )
        rows.append(block)
    sample_df = pd.DataFrame(rows).sort_values("patient")
    post = sample_df[sample_df["timing"] == "post"].copy()
    if len(post) != 12:
        raise SystemExit(f"expected 12 post patients, got {len(post)}")

    malig_elig = post[post["mal_n"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    epi_elig = post[post["epi_n"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    leftover_elig = post[post["leftover_n"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()
    malig_any = post[post["mal_n"] >= 1].copy()
    empty_malig_post = post.loc[post["mal_n"] == 0, "patient"].tolist()
    noisy_malig_post = post.loc[
        (post["mal_n"] >= 1) & (post["mal_n"] < MIN_CELLS_PER_SAMPLE_FOR_MEAN),
        "patient",
    ].tolist()

    primary_specs = [
        ("A3-malignant CLDN4 vs DPT", "mal_mean_CLDN4", "mal_mean_dpt"),
        ("A3-malignant CLDN4 vs AT2 score", "mal_mean_CLDN4", "mal_mean_AT2"),
        ("A3-malignant CLDN4 vs barrier/keratin (no CLDN4)", "mal_mean_CLDN4", "mal_mean_barrier"),
        ("A3-malignant CLDN4 vs TACSTD2 (companion)", "mal_mean_CLDN4", "mal_mean_TACSTD2"),
        ("A3-malignant AT2 vs DPT (control)", "mal_mean_AT2", "mal_mean_dpt"),
    ]
    primary = []
    for name, a, b in primary_specs:
        r = _spearman(malig_elig[a].to_numpy(), malig_elig[b].to_numpy())
        r["contrast"] = name
        r["compartment"] = "A3_malignant"
        r["cohort"] = "post ≥10 malignant-like"
        primary.append(r)
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q

    sensitivity = []
    for name, df, a, b, note in (
        (
            "epithelium CLDN4 vs DPT (post complete-case)",
            epi_elig,
            "epi_mean_CLDN4",
            "epi_mean_dpt",
            "companion; not malignant-only",
        ),
        (
            "leftover-epi CLDN4 vs DPT (post ≥10 leftover)",
            leftover_elig,
            "leftover_mean_CLDN4",
            "leftover_mean_dpt",
            "leftover only",
        ),
        (
            "A3-malignant CLDN4 vs DPT including 1–9 cell patients",
            malig_any,
            "mal_mean_CLDN4",
            "mal_mean_dpt",
            "noisy; includes 1–9 cell residuals",
        ),
        (
            "epithelium CLDN4 vs barrier/keratin (post)",
            epi_elig,
            "epi_mean_CLDN4",
            "epi_mean_barrier",
            "companion",
        ),
        (
            "epithelium CLDN4 vs TACSTD2 (companion, not a gate)",
            epi_elig,
            "epi_mean_CLDN4",
            "epi_mean_TACSTD2",
            "companion",
        ),
    ):
        r = _spearman(df[a].to_numpy(), df[b].to_numpy())
        r["contrast"] = name
        r["note"] = (r.get("note") or "") + (("; " + note) if note else "")
        sensitivity.append(r)

    mpr_rows = []
    for name, df, col, note in (
        (
            "A3-malignant DPT, post ≥10 malignant-like",
            malig_elig,
            "mal_mean_dpt",
            "primary malignant-like gate",
        ),
        (
            "A3-malignant DPT, post any malignant-like (≥1 cell)",
            malig_any,
            "mal_mean_dpt",
            "noisy; includes 1–9 cell residuals",
        ),
        (
            "epithelial DPT, post complete-case",
            epi_elig,
            "epi_mean_dpt",
            "companion; leftover+malignant on the graph",
        ),
        (
            "leftover-epi DPT, post ≥10 leftover",
            leftover_elig,
            "leftover_mean_dpt",
            "leftover only",
        ),
        (
            "A3-malignant CLDN4, post ≥10 malignant-like",
            malig_elig,
            "mal_mean_CLDN4",
            "CLDN4 vs MPR companion on the same n",
        ),
        (
            "epithelial CLDN4, post complete-case",
            epi_elig,
            "epi_mean_CLDN4",
            "CLDN4 vs MPR companion",
        ),
    ):
        nmpr = df[df["paper_group"] == "NMPR"][col]
        mpr = df[df["paper_group"] == "MPR"][col]
        r = exact_mwu(nmpr, mpr)
        r["contrast"] = f"NMPR vs MPR: {name}"
        r["score"] = col
        r["patients_nmpr"] = ",".join(
            df.loc[df["paper_group"] == "NMPR", "patient"].tolist()
        )
        r["patients_mpr"] = ",".join(
            df.loc[df["paper_group"] == "MPR", "patient"].tolist()
        )
        r["note"] = (r.get("note") or "") or note
        if r.get("note") == "" and note:
            r["note"] = note
        elif note and note not in (r.get("note") or ""):
            r["note"] = f"{r.get('note')}; {note}" if r.get("note") else note
        mpr_rows.append(r)

    paired_rows_data = []
    for pat, sub in adata.obs.groupby("patient", observed=True):
        if str(sub["paper_group"].iloc[0]) == "TN":
            continue
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_rows_data.append(
            {
                "patient": pat,
                "paper_group": str(sub["paper_group"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_rows_data)
    paired_rows = []
    if not paired_df.empty:
        for contrast, hi, lo in (
            ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
            ("AT2 score high vs low", "AT2_high", "AT2_low"),
            ("DPT high vs low", "dpt_high", "dpt_low"),
        ):
            r = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            r["contrast"] = contrast
            paired_rows.append(r)
    else:
        paired_rows = [
            {
                "contrast": "barrier/keratin (no CLDN4) high vs low",
                "n": 0,
                "W": None,
                "p": None,
                "delta_median": None,
            },
            {
                "contrast": "AT2 score high vs low",
                "n": 0,
                "W": None,
                "p": None,
                "delta_median": None,
            },
            {
                "contrast": "DPT high vs low",
                "n": 0,
                "W": None,
                "p": None,
                "delta_median": None,
            },
        ]

    # ----- figures -----
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1))
    sc.pl.paga(
        adata,
        color="leiden",
        ax=axes[0],
        show=False,
        frameon=False,
        title="PAGA (Leiden)",
    )
    sc.pl.umap(
        adata,
        color="expr_CLDN4",
        ax=axes[1],
        show=False,
        frameon=False,
        cmap="viridis",
        title="UMAP CLDN4",
    )
    ax = axes[2]
    c4_dpt = next(r for r in primary if r["contrast"] == "A3-malignant CLDN4 vs DPT")
    for grp, sub in malig_elig.groupby("paper_group"):
        ax.scatter(
            sub["mal_mean_dpt"],
            sub["mal_mean_CLDN4"],
            s=56,
            c=COLOR.get(grp, "0.4"),
            label=f"{grp} n={len(sub)}",
            zorder=3,
        )
        for _, row in sub.iterrows():
            ax.annotate(
                row["patient"],
                (row["mal_mean_dpt"], row["mal_mean_CLDN4"]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=7,
            )
    ax.set_xlabel("patient-mean DPT (A3-malignant-like)")
    ax.set_ylabel("patient-mean CLDN4 (A3-malignant-like)")
    ax.set_title(f"CLDN4 vs DPT  {_fmt_rho(c4_dpt)}")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(
        f"GSE207422 epithelium scored by CLDN4   "
        f"graph n_cells={adata.n_obs}  inferential post n=12 (malignant-like test n={len(malig_elig)})",
        fontsize=10,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    panels = (
        (
            axes[0],
            malig_elig,
            "mal_mean_dpt",
            "A3-malignant-like DPT",
            next(r for r in mpr_rows if "A3-malignant DPT, post ≥10" in r["contrast"]),
        ),
        (
            axes[1],
            epi_elig,
            "epi_mean_dpt",
            "epithelial DPT (companion)",
            next(r for r in mpr_rows if "epithelial DPT, post complete-case" in r["contrast"]),
        ),
    )
    for ax, df, col, ylab, stats_row in panels:
        for i, grp in enumerate(["NMPR", "MPR"]):
            vals = df.loc[df["paper_group"] == grp, col].to_numpy()
            ax.scatter(
                np.full(vals.size, i) + np.linspace(-0.08, 0.08, max(vals.size, 1))[: vals.size],
                vals,
                s=50,
                c=COLOR[grp],
                zorder=3,
            )
            for _, row in df[df["paper_group"] == grp].iterrows():
                ax.annotate(
                    row["patient"],
                    (i, row[col]),
                    textcoords="offset points",
                    xytext=(6, 2),
                    fontsize=7,
                )
        ax.set_xticks([0, 1], ["NMPR", "MPR"])
        ax.set_ylabel(ylab)
        ax.set_title(_fmt_mwu(stats_row), fontsize=8)
    fig.suptitle(
        "Patient-level DPT vs MPR. Left = malignant-like primary n; "
        "right = epithelium complete-case companion. Do not inflate cell n.",
        fontsize=10,
    )
    _save(fig, figdir / "fig_dpt_vs_mpr")

    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    plot_df = sample_df.set_index("patient")[["n_malignant", "n_leftover_epi"]]
    plot_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=["#b23a48", "#2a6f97"],
        width=0.85,
    )
    ax.set_ylabel("cells on analysis object")
    ax.set_title(
        f"Honest n: A3-malignant vs leftover epithelium after QC/cap  "
        f"n_cells={adata.n_obs}  GEO 15 / post 12"
    )
    ax.legend(["A3-malignant-like", "leftover epi"], frameon=False, fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", fontsize=8)
    _save(fig, figdir / "fig_honest_n")

    if not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("dpt_low", "dpt_high", "DPT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for grp, sub in paired_df.groupby("paper_group"):
                ax.scatter(sub[lo], sub[hi], s=44, c=COLOR.get(grp, "0.4"), label=grp)
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
        axes[0].set_title(
            _fmt_pair(next(r for r in paired_rows if r["contrast"].startswith("barrier"))),
            fontsize=8,
        )
        axes[1].set_title(
            _fmt_pair(next(r for r in paired_rows if r["contrast"].startswith("AT2"))),
            fontsize=8,
        )
        axes[2].set_title(
            _fmt_pair(next(r for r in paired_rows if r["contrast"].startswith("DPT"))),
            fontsize=8,
        )
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(
            f"EXTRA: within-sample CLDN4-high vs low  paired post n={len(paired_df)}",
            fontsize=11,
        )
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("compartment", "fig_umap_compartment", None),
        ("paper_group", "fig_umap_mpr", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("dpt_pseudotime", "fig_umap_dpt", "viridis"),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    c4_dpt = next(r for r in primary if r["contrast"] == "A3-malignant CLDN4 vs DPT")
    c4_at2 = next(r for r in primary if r["contrast"] == "A3-malignant CLDN4 vs AT2 score")
    c4_bar = next(r for r in primary if "barrier" in r["contrast"])
    dpt_mpr = next(r for r in mpr_rows if "A3-malignant DPT, post ≥10" in r["contrast"])
    dpt_mpr_epi = next(
        r for r in mpr_rows if "epithelial DPT, post complete-case" in r["contrast"]
    )
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_dpt = next(r for r in paired_rows if r["contrast"].startswith("DPT"))

    parts = [
        f"Patient-level A3-malignant CLDN4 vs leftover/AT2-rooted DPT: {_fmt_rho(c4_dpt)}.",
        f"CLDN4 vs AT2 score (malignant-like): {_fmt_rho(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt_rho(c4_bar)}.",
        f"A3-malignant DPT vs MPR: {_fmt_mwu(dpt_mpr)}.",
        f"Epithelial DPT vs MPR (complete-case companion, not the malignant-like n): {_fmt_mwu(dpt_mpr_epi)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt_pair(pair_bar)}.",
        f"Paired CLDN4-high vs low DPT: {_fmt_pair(pair_dpt)}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        f"Empty A3-malignant post patients: {empty_malig_post or 'none'}. "
        f"Noisy 1–9 cell post patients: {noisy_malig_post or 'none'}.",
        "Honest n is the patient n after the ≥10 malignant-like gate. "
        "This set is small. Not a TACSTD2 redo. No both-high gate. "
        "T/NK was not re-tested.",
    ]
    verdict = " ".join(parts)

    leiden_rows = []
    for lab in leiden_ids:
        sub = adata.obs[adata.obs["leiden"].astype(str) == lab]
        leiden_rows.append(
            {
                "leiden": lab,
                "n_cells": int(len(sub)),
                "n_malignant": int(sub["is_malig_a3"].sum()),
                "n_leftover": int(sub["is_leftover_epi"].sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "frac_MPR": float((sub["paper_group"] == "MPR").mean()),
                "frac_NMPR": float((sub["paper_group"] == "NMPR").mean()),
            }
        )
    leiden_df = pd.DataFrame(leiden_rows)
    paga_df = pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids)

    summary = {
        "accession": "GSE207422",
        "paper": "Hu et al. Genome Med 2023 PMID 36915183",
        "histology": "NSCLC mixed (LUAD+LUSC); neoadjuvant PD-1+platinum",
        "public_only": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "no_dual_high": True,
        "tnk_not_reaudited": True,
        "slingshot_run": False,
        "slingshot_reason": "R not present",
        "n_cells_matrix": n_cells_matrix,
        "n_genes_matrix": n_genes_matrix,
        "n_epithelial_all": n_epithelial_all,
        "n_malig_a3_all": n_malig_a3_all,
        "n_leftover_all": n_leftover_all,
        "n_cells_qc": n_cells_qc,
        "n_patients_qc": n_patients_qc,
        "n_cells": int(adata.n_obs),
        "n_samples_on_graph": int(sample_df.shape[0]),
        "n_post": 12,
        "n_malig_eligible_post": int(len(malig_elig)),
        "n_epi_eligible_post": int(len(epi_elig)),
        "n_samples_paired_tertile": int(len(paired_df)),
        "empty_malig_post": empty_malig_post,
        "noisy_malig_post": noisy_malig_post,
        "malig_eligible_patients": malig_elig["patient"].tolist(),
        "compartment_counts": adata.obs["compartment"].value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"]
        .astype(str)
        .value_counts()
        .to_dict(),
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "graph_capped": capped,
        "max_cells_per_patient": MAX_CELLS_PER_PATIENT if capped else None,
        "primary_cldn4_vs_dpt": primary,
        "primary_dpt_vs_mpr": mpr_rows,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "verdict": verdict,
    }

    sample_df.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tabdir / "patient_cldn4_vs_dpt.tsv", sep="\t", index=False)
    pd.DataFrame(mpr_rows).to_csv(tabdir / "patient_dpt_vs_mpr.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)
    malig_elig.to_csv(tabdir / "patient_malig_eligible.tsv", sep="\t", index=False)
    epi_elig.to_csv(tabdir / "patient_epi_eligible.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "paired_tertile.tsv", sep="\t", index=False)
    leiden_df.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    paga_df.to_csv(tabdir / "paga_connectivities.tsv", sep="\t")
    (tabdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})

    keep_cols = [
        c
        for c in [
            "patient",
            "Sample",
            "paper_group",
            "path_response",
            "timing",
            "Pathology",
            "compartment",
            "is_malig_a3",
            "is_leftover_epi",
            "leiden",
            "dpt_pseudotime",
            "cldn4_tertile",
            "expr_CLDN4",
            "expr_TACSTD2",
            "score_AT2",
            "score_barrier_keratin",
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
                "n_cells": int(adata.n_obs),
                "n_post": 12,
                "n_malig_eligible": int(len(malig_elig)),
                "root": root_info["rule"],
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
