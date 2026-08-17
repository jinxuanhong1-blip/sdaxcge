#!/usr/bin/env python3
"""CLDN4 vs real CytoTRACE2 potency on winning-pair malignant cells.

ADDITIVE. CLDN4 only. GSE131907 + GSE205335 author-malignant cells.
No dual-high TACSTD2∩CLDN4. No GSE148071. No GSE207422.
Patient/sample is the inferential unit. Cell-level ρ is descriptive.

Primary score: CytoTRACE2 (Kang et al. Nat Methods 2025) via cytotrace2-py.
If the package or model cannot run, a documented Gulati 2020 CytoTRACE
implementation (top-200 gene-count correlates + KNN smooth) is used and
labeled as such. Residual n_genes alone is never the primary score.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.neighbors import NearestNeighbors

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import AT2, BARRIER_KERATIN, COMPARATOR, FOCAL, MALIGNANT_LIKE  # noqa: E402

MIN_CELLS_FOR_MEAN = 20
MIN_CELLS_PER_ARM = 8
CAP_NOTE = 200


def _fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def _fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def _fmt(row, keys=("rho", "p")) -> str:
    bits = [f"n={row.get('n', 'NA')}"]
    if "rho" in keys:
        bits.append(f"ρ={_fmt_r(row.get('rho'))}")
    if "W" in keys:
        w = row.get("W")
        bits.append("W=NA" if w is None or not np.isfinite(w) else f"W={w:.1f}")
    if "delta" in keys:
        d = row.get("delta_median")
        bits.append("Δmed=NA" if d is None or not np.isfinite(d) else f"Δmed={d:+.3f}")
    if "p" in keys:
        bits.append(f"p={_fmt_p(row.get('p'))}")
    return ", ".join(bits)


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 3:
        return np.nan, np.nan, np.nan, int(len(a))
    d = a - b
    if np.allclose(d, 0):
        return 0.0, 1.0, 0.0, int(len(a))
    try:
        w, p = stats.wilcoxon(a, b, alternative="two-sided", zero_method="wilcox")
        return float(w), float(p), float(np.median(d)), int(len(a))
    except ValueError:
        return np.nan, np.nan, float(np.median(d)), int(len(a))


def bh(pvals):
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(ok.sum())
    tmp[order] = q
    out[ok] = tmp
    return out


def log1p_cp10k(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    tot = np.asarray(adata.obs["n_umi"], float)
    tot[tot <= 0] = np.nan
    col = X[:, adata.var_names.get_loc(gene)]
    if sparse.issparse(col):
        val = np.asarray(col.toarray()).ravel()
    else:
        val = np.asarray(col).ravel()
    return np.log1p(val / tot * 1e4)


def module_score(adata, genes) -> np.ndarray:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.full(adata.n_obs, np.nan)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    tot = np.asarray(adata.obs["n_umi"], float)
    tot[tot <= 0] = np.nan
    idx = [adata.var_names.get_loc(g) for g in present]
    sub = X[:, idx]
    if sparse.issparse(sub):
        raw = np.asarray(sub.toarray())
    else:
        raw = np.asarray(sub)
    expr = np.log1p(raw / tot[:, None] * 1e4)
    z = (expr - np.nanmean(expr, axis=0)) / (np.nanstd(expr, axis=0) + 1e-8)
    return np.nanmean(z, axis=1)


def gulati_cytotrace(X, n_top: int = 200, n_neighbors: int = 10, seed: int = 14) -> np.ndarray:
    """Documented CytoTRACE (Gulati et al. Science 2020), not residual n_genes.

    1. Per-cell gene counts (UMI>0).
    2. Pearson of each gene vs gene counts.
    3. Mean of the top-N correlated genes.
    4. KNN mean smooth.
    5. Rank-scale to [0, 1]. Higher = more potent / less differentiated.
    """
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()
    n_cells, n_genes = X.shape
    gene_counts = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
    # genes expressed in at least 10 cells
    nz_cells = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = nz_cells >= 10
    if keep.sum() < n_top:
        keep = nz_cells >= 3
    Xs = X[:, keep]
    # Pearson of log1p counts vs gene_counts
    logx = Xs.copy()
    logx.data = np.log1p(logx.data)
    gc = gene_counts - gene_counts.mean()
    gc_ss = float(np.dot(gc, gc))
    if gc_ss <= 0:
        ranks = stats.rankdata(gene_counts)
        return ranks / n_cells
    col_means = np.asarray(logx.mean(axis=0)).ravel()
    # cov(gene, gc) = (X.T @ gc) / (n-1) after centering gene
    xt_gc = np.asarray(logx.T.dot(gc)).ravel()  # sum_i x_ij * gc_i
    # sum (x - mean_x) * gc = X.T@gc - mean_x * sum(gc); sum(gc)=0 so X.T@gc
    col_ss = np.asarray(logx.power(2).sum(axis=0)).ravel() - n_cells * col_means**2
    denom = np.sqrt(np.clip(col_ss, 0, None) * gc_ss)
    corr = np.zeros(keep.sum(), dtype=float)
    ok = denom > 1e-12
    corr[ok] = xt_gc[ok] / denom[ok]
    top = np.argsort(corr)[::-1][: min(n_top, keep.sum())]
    score = np.asarray(logx[:, top].mean(axis=1)).ravel()
    # KNN smooth on log1p PCA-ish: use the same top-gene space
    feat = np.asarray(logx[:, top].toarray())
    feat = feat - feat.mean(axis=0)
    std = feat.std(axis=0)
    std[std == 0] = 1
    feat = feat / std
    k = min(n_neighbors, n_cells - 1)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nn.fit(feat)
    idx = nn.kneighbors(feat, return_distance=False)
    smooth = score[idx].mean(axis=1)
    return stats.rankdata(smooth) / n_cells


def write_cytotrace2_input(adata, path: Path, min_cells: int = 10) -> Path:
    """Write genes x cells TSV (raw UMI, not log) for cytotrace2-py."""
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()
    nz = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = nz >= min_cells
    genes = adata.var_names.to_numpy()[keep]
    Xs = X[:, keep].T.tocsr()  # genes x cells
    cells = np.asarray(adata.obs_names, dtype=str)
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"write CytoTRACE2 input {Xs.shape[0]} genes x {Xs.shape[1]} cells -> {path}", flush=True)
    with path.open("w") as handle:
        handle.write("Gene\t" + "\t".join(cells) + "\n")
        for i, gene in enumerate(genes):
            row = Xs.getrow(i)
            dense = np.zeros(Xs.shape[1], dtype=np.float32)
            dense[row.indices] = row.data
            handle.write(gene + "\t" + "\t".join(f"{v:.4g}" if v else "0" for v in dense) + "\n")
            if (i + 1) % 2000 == 0:
                print(f"  wrote {i+1} genes", flush=True)
    print(f"wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


def run_cytotrace2(adata, work: Path) -> tuple[pd.DataFrame | None, dict]:
    info = {
        "package": "cytotrace2-py",
        "citation": "Kang et al. Nat Methods 2025 doi:10.1038/s41592-025-02857-2",
        "ran": False,
        "reason": None,
    }
    try:
        from cytotrace2_py.cytotrace2_py import cytotrace2
    except Exception as exc:
        info["reason"] = f"import failed: {exc}"
        print(info["reason"], flush=True)
        return None, info
    expr_path = work / "cytotrace2_input.tsv"
    ann_path = work / "cytotrace2_annotation.tsv"
    out_dir = work / "cytotrace2_out"
    if not expr_path.is_file():
        write_cytotrace2_input(adata, expr_path)
    pd.DataFrame(
        {
            "cell": adata.obs_names.astype(str),
            "phenotype": adata.obs["unit_id"].astype(str),
        }
    ).to_csv(ann_path, sep="\t", index=False)
    try:
        print("running cytotrace2(species=human, disable_plotting, max_cores=1)", flush=True)
        result = cytotrace2(
            str(expr_path),
            annotation_path=str(ann_path),
            species="human",
            batch_size=4000,
            smooth_batch_size=800,
            disable_plotting=True,
            disable_parallelization=True,
            max_cores=1,
            seed=14,
            output_dir=str(out_dir),
        )
    except TypeError:
        # older signature without some kwargs
        try:
            result = cytotrace2(
                str(expr_path),
                annotation_path=str(ann_path),
                species="human",
                batch_size=4000,
                smooth_batch_size=800,
                disable_plotting=True,
                max_cores=1,
                seed=14,
                output_dir=str(out_dir),
            )
        except Exception as exc:
            info["reason"] = f"cytotrace2() failed: {exc}\n{traceback.format_exc()}"
            print(info["reason"], flush=True)
            return None, info
    except Exception as exc:
        info["reason"] = f"cytotrace2() failed: {exc}\n{traceback.format_exc()}"
        print(info["reason"], flush=True)
        return None, info
    if result is None:
        info["reason"] = "cytotrace2() returned None"
        return None, info
    if not isinstance(result, pd.DataFrame):
        result = pd.DataFrame(result)
    info["ran"] = True
    info["n_scored"] = int(len(result))
    info["columns"] = list(result.columns)
    print(f"CytoTRACE2 scored {len(result)} cells; cols={list(result.columns)}", flush=True)
    return result, info


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(path: Path, S: dict) -> None:
    prim = {r["contrast"]: r for r in S["primary"]}
    sens = {r["contrast"]: r for r in S["sensitivity"]}
    paired = {r["contrast"]: r for r in S["paired"]}
    c4p = prim["CLDN4 vs CytoTRACE2"]
    c4b = prim["CLDN4 vs barrier/keratin (no CLDN4)"]
    p_bar = prim["CytoTRACE2 vs barrier/keratin (no CLDN4)"]
    pair_p = paired["potency high vs low"]
    pair_b = paired["barrier/keratin (no CLDN4) high vs low"]
    s131 = sens["GSE131907-only CLDN4 vs CytoTRACE2"]
    s205 = sens["GSE205335-only CLDN4 vs CytoTRACE2"]
    adc = sens.get("GSE205335 ADC-only CLDN4 vs CytoTRACE2", {})
    method = S["potency_method"]
    direction = (
        "CLDN4-high cells have **lower** potency (more differentiated / barrier-locked)"
        if (c4p.get("rho") is not None and c4p["rho"] < 0 and c4p.get("p", 1) < 0.05)
        or (pair_p.get("delta_median") is not None and pair_p["delta_median"] < 0 and pair_p.get("p", 1) < 0.05)
        else "the data do **not** support lower potency in CLDN4-high malignant cells at the patient unit"
    )
    if method["ran"] and method.get("name") == "CytoTRACE2":
        method_line = (
            f"Primary potency = **CytoTRACE2** (`cytotrace2-py` {method.get('version', '')}, "
            "Kang et al. *Nat Methods* 2025). Score 0 = differentiated, 1 = totipotent. "
            "This is not a residual-n_genes dump."
        )
    else:
        method_line = (
            f"CytoTRACE2 did **not** complete ({method.get('reason', 'unknown')}). "
            "Primary potency = documented **Gulati 2020 CytoTRACE** "
            "(top-200 gene-count correlates + KNN smooth, rank-scaled to [0,1]). "
            "Not residual n_genes alone."
        )

    lines = [
        "# Finding — winning-pair GSE131907+GSE205335 malignant cells, CLDN4 vs CytoTRACE2",
        "",
        "ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search "
        "(author malignant CLDN4 %pos GSE131907+GSE205335 vs T/NK). "
        "This folder does **not** use GSE148071. GSE207422 is not added. "
        "No TACSTD2∩CLDN4 dual-high gate.",
        "",
        f"{method_line} Inferential unit = **patient** (GSE131907 `Sample`, GSE205335 `patient`). "
        "Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**.",
        "",
        f"**Question.** Are CLDN4-high author-malignant cells more differentiated / barrier-locked "
        f"(lower potency) than CLDN4-low cells? **Answer:** {direction}.",
        "",
        "## Verdict",
        "",
        f"Patient-level CLDN4 vs CytoTRACE2/potency: {_fmt(c4p)}. "
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4b)}. "
        f"Potency vs barrier/keratin: {_fmt(p_bar)}. "
        f"Paired within-patient CLDN4-high vs low potency: {_fmt(pair_p, keys=('W', 'delta', 'p'))}. "
        f"Paired barrier/keratin: {_fmt(pair_b, keys=('W', 'delta', 'p'))}. "
        f"GSE131907-only: {_fmt(s131)}. GSE205335-only: {_fmt(s205)}. "
        "tLung contributes **0** author-malignant cells (Kim labels those cells tS1/tS2/tS3, not "
        "`Malignant cells`). Not a TACSTD2 redo. No both-high gate. GSE148071 not used.",
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{S['cap_per_unit']}/unit): "
        f"**n_cells = {S['n_cells']}** (GSE131907 {S['n_cells_gse131907']}, GSE205335 {S['n_cells_gse205335']}).",
        f"- Units (GSE131907 Sample + GSE205335 patient): **n_units = {S['n_units']}** "
        f"(GSE131907 {S['n_units_gse131907']}, GSE205335 {S['n_units_gse205335']}).",
        f"- Units with ≥{MIN_CELLS_FOR_MEAN} malignant cells used for Spearman: **n = {S['n_units_eligible']}**.",
        f"- Units with ≥{MIN_CELLS_PER_ARM} cells in both CLDN4-high and CLDN4-low arms: **n = {S['n_units_paired']}**.",
        f"- Author malignant rule: GSE131907 `Cell_subtype==Malignant cells` "
        f"(catalog malignant {S['catalog_gse131907_malignant']}; tLung author-malignant = 0). "
        f"GSE205335 `lineage.sub==Malignant cells` on non-normal tissues "
        f"(catalog {S['catalog_gse205335_malignant']}).",
        f"- Histology (cells): {S['histology_counts']}.",
        f"- CytoTRACE2 potency categories (cells): {S['potency_category_counts']}.",
        f"- CLDN4 tertile cells: {S['cldn4_tertile_counts']}.",
        f"- Genes absent from locked sets: {S['genes_absent']}.",
        "- GSE148071 not used. GSE207422 not used. Dual-high not used.",
        f"- CytoTRACE2 ran: **{method.get('ran')}**. Fallback/sensitivity Gulati 2020 ran: **{S['gulati_ran']}**.",
        "",
        "## Locked choices",
        "",
        f"- Malignant = author label only. No CopyKAT / inferCNV re-call.",
        f"- Cap ≤{S['cap_per_unit']} cells / unit after QC (honest n reports catalog vs analysis).",
        "- CLDN4 = log1p(CP10k) from raw UMI.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- Inferential n: GSE131907 `Sample` + GSE205335 `patient`.",
        "- Extra figure: within-unit CLDN4-high vs low potency and barrier (min 8 cells/arm).",
        "- Unused: GSE148071; GSE207422; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a potency label.",
        "",
        "## Primary (patient-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in S["primary"]:
        q = "NA" if r.get("q") is None or not np.isfinite(r.get("q", np.nan)) else (
            f"{r['q']:.2e}" if r["q"] < 0.001 else f"{r['q']:.4f}"
        )
        lines.append(
            f"| {r['contrast']} | {r['n']} | {_fmt_r(r.get('rho'))} | {_fmt_p(r.get('p'))} | {q} |"
        )
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_units | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in S["sensitivity"]:
        lines.append(
            f"| {r['contrast']} | {r['n']} | {_fmt_r(r.get('rho'))} | {_fmt_p(r.get('p'))} |"
        )
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)",
        "",
        f"Emitted: **{S['extra_emitted']}**. "
        "Rule: patient-level Spearman(CLDN4, potency) or Spearman(CLDN4, barrier) p<0.05, "
        "or any paired Wilcoxon p<0.05, or n_paired≥4.",
        "",
        "| Paired contrast (high − low) | n_units | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in S["paired"]:
        d = r.get("delta_median")
        ds = "NA" if d is None or not np.isfinite(d) else f"{d:.3f}"
        lines.append(f"| {r['contrast']} | {r['n']} | {ds} | {_fmt_p(r.get('p'))} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- tLung is not in the author-malignant universe. Do not quote this as a tLung AT2 trajectory.",
        "- GSE205335 mixes ADC / SQ / SCLC / NUT. ADC-only n is in the sensitivity table.",
        "- Malignant is the author label, not CNV.",
        "- This is not an ICI / MPR / RECIST test.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- GSE148071 is not used.",
        "- Do not write “CLDN4 marks a stem-like malignant state” unless the patient-level potency ρ is negative and significant.",
        "",
        "## Outputs",
        "",
        "- `results/tables/patient_cldn4_vs_potency.tsv` — **done criterion**",
        "- `results/tables/patient_means.tsv`",
        "- `results/tables/paired_high_vs_low.tsv`",
        "- `results/tables/stats.tsv`",
        "- `results/figures/fig_patient_cldn4_vs_potency.png`",
        "- `results/figures/fig_extra_paired_potency.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 -m venv /tmp/winpair_ct2_venv",
        " /tmp/winpair_ct2_venv/bin/pip install -r methods/winpair_cytotrace2_cldn4/requirements.txt",
        "/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/download.py \\",
        "  --out /tmp/winpair_cytotrace2_cldn4",
        "/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/extract_malignant.py \\",
        "  --data /tmp/winpair_cytotrace2_cldn4 \\",
        "  --out /tmp/winpair_cytotrace2_cldn4/malignant.h5ad",
        "/tmp/winpair_ct2_venv/bin/python methods/winpair_cytotrace2_cldn4/scripts/analyze.py \\",
        "  --input /tmp/winpair_cytotrace2_cldn4/malignant.h5ad \\",
        "  --outdir methods/winpair_cytotrace2_cldn4/results \\",
        "  --finding methods/winpair_cytotrace2_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))
    print(f"wrote {path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path("/tmp/winpair_cytotrace2_cldn4/malignant.h5ad"))
    p.add_argument("--outdir", type=Path, default=Path("methods/winpair_cytotrace2_cldn4/results"))
    p.add_argument("--finding", type=Path, default=Path("methods/winpair_cytotrace2_cldn4/FINDING.md"))
    p.add_argument("--work", type=Path, default=Path("/tmp/winpair_cytotrace2_cldn4/work"))
    args = p.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.input)
    if "n_umi" not in adata.obs:
        X = adata.layers["counts"] if "counts" in adata.layers else adata.X
        adata.obs["n_umi"] = np.asarray(X.sum(axis=1)).ravel()
        adata.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()

    adata.obs["expr_CLDN4"] = log1p_cp10k(adata, "CLDN4")
    adata.obs["expr_TACSTD2"] = log1p_cp10k(adata, "TACSTD2")
    adata.obs["score_barrier"] = module_score(adata, BARRIER_KERATIN)
    adata.obs["score_malignant_like"] = module_score(adata, MALIGNANT_LIKE)
    adata.obs["score_AT2"] = module_score(adata, AT2)
    q1, q2 = np.nanquantile(adata.obs["expr_CLDN4"], [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[adata.obs["expr_CLDN4"] <= q1] = "low"
    tert[adata.obs["expr_CLDN4"] >= q2] = "high"
    adata.obs["cldn4_tertile"] = tert

    absent = {
        "barrier_keratin": [g for g in BARRIER_KERATIN if g not in adata.var_names],
        "malignant_like": [g for g in MALIGNANT_LIKE if g not in adata.var_names],
        "AT2": [g for g in AT2 if g not in adata.var_names],
        "focal": [g for g in FOCAL if g not in adata.var_names],
        "comparator": [g for g in COMPARATOR if g not in adata.var_names],
    }

    args.work.mkdir(parents=True, exist_ok=True)
    ct2_df, ct2_info = run_cytotrace2(adata, args.work)
    gulati = gulati_cytotrace(adata.layers["counts"] if "counts" in adata.layers else adata.X)
    adata.obs["gulati_cytotrace"] = gulati

    potency_col = None
    potency_name = None
    if ct2_df is not None and len(ct2_df):
        df = ct2_df.copy()
        if df.index.name is None and "Unnamed: 0" in df.columns:
            df = df.set_index("Unnamed: 0")
        # join on cell id
        idx = df.index.astype(str)
        if not set(adata.obs_names.astype(str)).intersection(set(idx)):
            # try first column
            first = df.columns[0]
            if df[first].astype(str).isin(adata.obs_names.astype(str)).any():
                df = df.set_index(first)
        df.index = df.index.astype(str)
        score_col = None
        for cand in ("CytoTRACE2_Score", "CytoTRACE2_score", "cytotrace2_score", "preKNN_CytoTRACE2_Score"):
            if cand in df.columns:
                score_col = cand
                break
        if score_col is None:
            num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            score_col = num[0] if num else None
        if score_col is not None:
            mapped = df.reindex(adata.obs_names.astype(str))
            adata.obs["CytoTRACE2_Score"] = pd.to_numeric(mapped[score_col], errors="coerce").to_numpy()
            if "CytoTRACE2_Potency" in mapped.columns:
                adata.obs["CytoTRACE2_Potency"] = mapped["CytoTRACE2_Potency"].astype(str).to_numpy()
            elif "CytoTRACE2_potency" in mapped.columns:
                adata.obs["CytoTRACE2_Potency"] = mapped["CytoTRACE2_potency"].astype(str).to_numpy()
            else:
                sc = adata.obs["CytoTRACE2_Score"]
                cat = pd.cut(
                    sc,
                    bins=[-np.inf, 1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6, np.inf],
                    labels=["Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"],
                )
                adata.obs["CytoTRACE2_Potency"] = cat.astype(str)
            if "CytoTRACE2_Relative" in mapped.columns:
                adata.obs["CytoTRACE2_Relative"] = pd.to_numeric(mapped["CytoTRACE2_Relative"], errors="coerce").to_numpy()
            n_ok = int(np.isfinite(adata.obs["CytoTRACE2_Score"]).sum())
            if n_ok >= 50:
                potency_col = "CytoTRACE2_Score"
                potency_name = "CytoTRACE2"
                ct2_info["ran"] = True
                ct2_info["n_joined"] = n_ok
    if potency_col is None:
        adata.obs["CytoTRACE2_Score"] = adata.obs["gulati_cytotrace"]
        adata.obs["CytoTRACE2_Potency"] = pd.cut(
            adata.obs["CytoTRACE2_Score"],
            bins=[-np.inf, 1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6, np.inf],
            labels=["Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"],
        ).astype(str)
        potency_col = "CytoTRACE2_Score"
        potency_name = "Gulati2020_CytoTRACE"
        ct2_info["ran"] = False
        if not ct2_info.get("reason"):
            ct2_info["reason"] = "CytoTRACE2 scores could not be joined; using Gulati 2020"

    try:
        import importlib.metadata

        ct2_info["version"] = importlib.metadata.version("cytotrace2-py")
    except Exception:
        ct2_info["version"] = "unknown"
    ct2_info["name"] = potency_name
    ct2_info["primary_column"] = potency_col

    # UMAP for extra figures
    try:
        import scanpy as sc

        sc.settings.verbosity = 1
        adata_u = adata.copy()
        sc.pp.normalize_total(adata_u, target_sum=1e4)
        sc.pp.log1p(adata_u)
        sc.pp.highly_variable_genes(adata_u, n_top_genes=2000, flavor="seurat_v3", layer="counts")
        sc.pp.pca(adata_u, n_comps=30, use_highly_variable=True)
        sc.pp.neighbors(adata_u, n_neighbors=20, n_pcs=30)
        sc.tl.umap(adata_u)
        adata.obsm["X_umap"] = adata_u.obsm["X_umap"]
        umap_ok = True
    except Exception as exc:
        print(f"UMAP skipped: {exc}", flush=True)
        umap_ok = False

    rows = []
    for uid, g in adata.obs.groupby("unit_id", observed=True):
        cldn = g["expr_CLDN4"]
        hi = g.loc[g["cldn4_tertile"] == "high"]
        lo = g.loc[g["cldn4_tertile"] == "low"]
        # within-unit tertiles if global tertiles are unbalanced
        if len(hi) < MIN_CELLS_PER_ARM or len(lo) < MIN_CELLS_PER_ARM:
            local_q1, local_q2 = np.nanquantile(cldn, [1 / 3, 2 / 3]) if len(g) >= 3 * MIN_CELLS_PER_ARM else (np.nan, np.nan)
            hi = g.loc[cldn >= local_q2] if np.isfinite(local_q2) else hi
            lo = g.loc[cldn <= local_q1] if np.isfinite(local_q1) else lo
        rows.append(
            {
                "unit_id": uid,
                "dataset": g["dataset"].iloc[0],
                "patient_id": g["patient_id"].iloc[0],
                "origin": g["Sample_Origin"].iloc[0],
                "histology": g["histology"].iloc[0],
                "n_malignant": int(len(g)),
                "mean_CLDN4": float(cldn.mean()),
                "median_CLDN4": float(cldn.median()),
                "pct_CLDN4_pos": float((cldn > 0).mean() * 100.0),
                "mean_potency": float(g[potency_col].mean()),
                "median_potency": float(g[potency_col].median()),
                "mean_gulati": float(g["gulati_cytotrace"].mean()),
                "frac_differentiated": float((g["CytoTRACE2_Potency"] == "Differentiated").mean()),
                "mean_barrier": float(g["score_barrier"].mean()),
                "mean_malignant_like": float(g["score_malignant_like"].mean()),
                "mean_AT2": float(g["score_AT2"].mean()),
                "mean_TACSTD2": float(g["expr_TACSTD2"].mean()),
                "n_CLDN4_high": int(len(hi)),
                "n_CLDN4_low": int(len(lo)),
                "mean_potency_CLDN4_high": float(hi[potency_col].mean()) if len(hi) else np.nan,
                "mean_potency_CLDN4_low": float(lo[potency_col].mean()) if len(lo) else np.nan,
                "mean_barrier_CLDN4_high": float(hi["score_barrier"].mean()) if len(hi) else np.nan,
                "mean_barrier_CLDN4_low": float(lo["score_barrier"].mean()) if len(lo) else np.nan,
                "mean_gulati_CLDN4_high": float(hi["gulati_cytotrace"].mean()) if len(hi) else np.nan,
                "mean_gulati_CLDN4_low": float(lo["gulati_cytotrace"].mean()) if len(lo) else np.nan,
                "eligible": int(len(g) >= MIN_CELLS_FOR_MEAN),
            }
        )
    patient = pd.DataFrame(rows).sort_values(["dataset", "unit_id"])
    elig = patient[patient["eligible"] == 1].copy()

    def add_row(store, contrast, x, y, note=""):
        r, pval, n = spear(x, y)
        store.append({"contrast": contrast, "n": n, "rho": r, "p": pval, "note": note})

    primary = []
    add_row(primary, "CLDN4 vs CytoTRACE2", elig["mean_CLDN4"], elig["mean_potency"])
    add_row(primary, "CLDN4 vs barrier/keratin (no CLDN4)", elig["mean_CLDN4"], elig["mean_barrier"])
    add_row(primary, "CytoTRACE2 vs barrier/keratin (no CLDN4)", elig["mean_potency"], elig["mean_barrier"])
    add_row(primary, "CLDN4 vs frac Differentiated", elig["mean_CLDN4"], elig["frac_differentiated"])
    add_row(primary, "CLDN4 vs malignant-like", elig["mean_CLDN4"], elig["mean_malignant_like"])
    add_row(primary, "CLDN4 vs TACSTD2 (comparator)", elig["mean_CLDN4"], elig["mean_TACSTD2"])
    add_row(primary, "CLDN4 vs Gulati2020 CytoTRACE", elig["mean_CLDN4"], elig["mean_gulati"])
    qs = bh([r["p"] for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = float(q) if np.isfinite(q) else np.nan

    sensitivity = []
    for ds, label in (("GSE131907", "GSE131907-only"), ("GSE205335", "GSE205335-only")):
        sub = elig[elig["dataset"] == ds]
        add_row(sensitivity, f"{label} CLDN4 vs CytoTRACE2", sub["mean_CLDN4"], sub["mean_potency"])
        add_row(sensitivity, f"{label} CLDN4 vs barrier/keratin (no CLDN4)", sub["mean_CLDN4"], sub["mean_barrier"])
    adc = elig[(elig["dataset"] == "GSE205335") & (elig["histology"] == "ADC")]
    add_row(sensitivity, "GSE205335 ADC-only CLDN4 vs CytoTRACE2", adc["mean_CLDN4"], adc["mean_potency"])
    adcsq = elig[(elig["dataset"] == "GSE205335") & (elig["histology"].isin(["ADC", "SQ"]))]
    add_row(sensitivity, "GSE205335 ADC+SQ CLDN4 vs CytoTRACE2", adcsq["mean_CLDN4"], adcsq["mean_potency"])
    for origin in sorted(elig.loc[elig["dataset"] == "GSE131907", "origin"].unique()):
        sub = elig[(elig["dataset"] == "GSE131907") & (elig["origin"] == origin)]
        add_row(sensitivity, f"GSE131907 {origin} CLDN4 vs CytoTRACE2", sub["mean_CLDN4"], sub["mean_potency"])

    paired_df = elig[
        (elig["n_CLDN4_high"] >= MIN_CELLS_PER_ARM) & (elig["n_CLDN4_low"] >= MIN_CELLS_PER_ARM)
    ].copy()
    paired = []
    for name, hi, lo in (
        ("potency high vs low", "mean_potency_CLDN4_high", "mean_potency_CLDN4_low"),
        ("barrier/keratin (no CLDN4) high vs low", "mean_barrier_CLDN4_high", "mean_barrier_CLDN4_low"),
        ("Gulati2020 high vs low", "mean_gulati_CLDN4_high", "mean_gulati_CLDN4_low"),
    ):
        w, pw, delta, n = wilcoxon_paired(paired_df[hi], paired_df[lo])
        paired.append(
            {
                "contrast": name,
                "n": n,
                "W": w,
                "p": pw,
                "delta_median": delta,
                "note": "within-unit CLDN4-high minus CLDN4-low; negative potency Δ = high more differentiated",
            }
        )

    extra_emitted = (
        any(r["p"] < 0.05 for r in primary if np.isfinite(r.get("p", np.nan)))
        or any(r["p"] < 0.05 for r in paired if np.isfinite(r.get("p", np.nan)))
        or len(paired_df) >= 4
    )

    outdir = args.outdir
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    # done-criterion table
    done = elig[
        [
            "unit_id",
            "dataset",
            "patient_id",
            "origin",
            "histology",
            "n_malignant",
            "mean_CLDN4",
            "mean_potency",
            "median_potency",
            "frac_differentiated",
            "mean_barrier",
            "mean_gulati",
            "mean_potency_CLDN4_high",
            "mean_potency_CLDN4_low",
        ]
    ].copy()
    done = done.rename(columns={"mean_potency": "mean_CytoTRACE2_or_potency"})
    done.to_csv(tabdir / "patient_cldn4_vs_potency.tsv", sep="\t", index=False)
    patient.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    paired_df.to_csv(tabdir / "paired_high_vs_low.tsv", sep="\t", index=False)
    stats_rows = []
    for r in primary:
        stats_rows.append({**r, "family": "primary"})
    for r in sensitivity:
        stats_rows.append({**r, "family": "sensitivity"})
    for r in paired:
        stats_rows.append(
            {
                "contrast": r["contrast"],
                "n": r["n"],
                "rho": r.get("delta_median"),
                "p": r["p"],
                "note": r["note"],
                "family": "paired",
            }
        )
    pd.DataFrame(stats_rows).to_csv(tabdir / "stats.tsv", sep="\t", index=False)
    adata.obs[
        [
            "dataset",
            "unit_id",
            "patient_id",
            "Sample_Origin",
            "histology",
            "expr_CLDN4",
            "cldn4_tertile",
            potency_col,
            "gulati_cytotrace",
            "CytoTRACE2_Potency",
            "score_barrier",
            "n_umi",
            "n_genes",
        ]
    ].to_csv(tabdir / "cell_obs.tsv", sep="\t")

    colors = {"GSE131907": "#2a6f97", "GSE205335": "#b23a48"}
    c4p = next(r for r in primary if r["contrast"] == "CLDN4 vs CytoTRACE2")
    c4b = next(r for r in primary if r["contrast"] == "CLDN4 vs barrier/keratin (no CLDN4)")
    pbar = next(r for r in primary if r["contrast"] == "CytoTRACE2 vs barrier/keratin (no CLDN4)")

    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_potency"], sub["mean_CLDN4"], s=56, c=col, label=f"{ds} n={len(sub)}", edgecolors="white", linewidths=0.4)
    ax.set_xlabel(f"patient-mean {potency_name} (higher = more potent)")
    ax.set_ylabel("patient-mean CLDN4 log1p(CP10k)")
    ax.set_title(f"CLDN4 vs potency  {_fmt(c4p)}")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle(
        f"Winning-pair author-malignant cells   n_cells={adata.n_obs}  n_units={len(elig)}  {potency_name}",
        fontsize=10,
    )
    _save(fig, figdir / "fig_patient_cldn4_vs_potency")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for ax, x, y, row, xlab in (
        (axes[0], "mean_barrier", "mean_CLDN4", c4b, "patient-mean barrier/keratin (no CLDN4)"),
        (axes[1], "mean_barrier", "mean_potency", pbar, "patient-mean barrier/keratin (no CLDN4)"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel("patient-mean CLDN4" if y == "mean_CLDN4" else f"patient-mean {potency_name}")
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: CLDN4 / potency vs barrier (CLDN4 excluded from the score)", fontsize=11)
    _save(fig, figdir / "fig_extra_cldn4_vs_barrier")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    unit_ct = patient.groupby("dataset").size()
    axes[0].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(unit_ct)])
    axes[0].set_ylabel("units")
    axes[0].set_title(f"Units={len(patient)}  eligible≥{MIN_CELLS_FOR_MEAN}: {len(elig)}")
    cell_ct = adata.obs.groupby("dataset").size()
    axes[1].bar(cell_ct.index.astype(str), cell_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(cell_ct)])
    axes[1].set_ylabel("cells (capped)")
    axes[1].set_title(f"n_cells={adata.n_obs}  cap≤{CAP_NOTE}/unit")
    fig.suptitle("Honest n — author malignant, winning pair, no GSE148071", fontsize=11)
    _save(fig, figdir / "fig_honest_n")

    if extra_emitted and not paired_df.empty:
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
        panels = (
            ("mean_potency_CLDN4_low", "mean_potency_CLDN4_high", f"{potency_name}", paired[0]),
            ("mean_barrier_CLDN4_low", "mean_barrier_CLDN4_high", "barrier/keratin (no CLDN4)", paired[1]),
        )
        for ax, (lo, hi, lab, row) in zip(axes, panels):
            for ds, col in colors.items():
                sub = paired_df[paired_df["dataset"] == ds]
                ax.scatter(sub[lo], sub[hi], s=48, c=col, label=f"{ds} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
            ax.set_title(_fmt(row, keys=("W", "delta", "p")))
            ax.legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-patient CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_paired_potency")

        fig, ax = plt.subplots(figsize=(5.4, 4.2))
        for ds, col in colors.items():
            sub = paired_df[paired_df["dataset"] == ds]
            d = sub["mean_potency_CLDN4_high"] - sub["mean_potency_CLDN4_low"]
            ax.scatter(np.arange(len(d)) + (0 if ds == "GSE131907" else 0.15), d, s=40, c=col, label=ds)
        ax.axhline(0, ls="--", c="0.5", lw=1)
        ax.set_ylabel("Δ potency (CLDN4-high − low)")
        ax.set_xlabel("paired units (jittered index)")
        ax.set_title("EXTRA: negative Δ = CLDN4-high more differentiated")
        ax.legend(frameon=False, fontsize=8)
        _save(fig, figdir / "fig_extra_delta_potency")

    # potency category by CLDN4 tertile
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    cat_order = ["Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"]
    tab = (
        adata.obs.assign(cldn4_tertile=pd.Categorical(adata.obs["cldn4_tertile"], ["low", "mid", "high"], ordered=True))
        .groupby(["cldn4_tertile", "CytoTRACE2_Potency"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for c in cat_order:
        if c not in tab.columns:
            tab[c] = 0
    tab = tab[cat_order]
    frac = tab.div(tab.sum(axis=1), axis=0)
    frac.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
    ax.set_ylabel("fraction of cells")
    ax.set_title("EXTRA: potency category by CLDN4 tertile (cells; descriptive)")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1))
    _save(fig, figdir / "fig_extra_potency_category")

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    data_box = [elig.loc[elig["dataset"] == ds, "mean_potency"].to_numpy() for ds in ("GSE131907", "GSE205335") if (elig["dataset"] == ds).any()]
    labels = [ds for ds in ("GSE131907", "GSE205335") if (elig["dataset"] == ds).any()]
    bp = ax.boxplot(data_box, labels=labels, patch_artist=True)
    for patch, col in zip(bp["boxes"], ["#2a6f97", "#b23a48"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.4)
    ax.set_ylabel(f"patient-mean {potency_name}")
    ax.set_title("EXTRA: potency by cohort (patient means)")
    _save(fig, figdir / "fig_extra_cohort_potency")

    if umap_ok:
        for color, fname, cmap in (
            ("expr_CLDN4", "fig_umap_cldn4", "viridis"),
            (potency_col, "fig_umap_potency", "magma"),
            ("dataset", "fig_umap_dataset", None),
            ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
            ("score_barrier", "fig_umap_barrier_no_cldn4", "viridis"),
        ):
            fig, ax = plt.subplots(figsize=(4.8, 4.2))
            vals = adata.obs[color]
            xy = adata.obsm["X_umap"]
            if cmap:
                sca = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=4, cmap=cmap, linewidths=0)
                fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
            else:
                cats = vals.astype(str)
                for i, lab in enumerate(sorted(cats.unique())):
                    m = cats == lab
                    ax.scatter(xy[m, 0], xy[m, 1], s=4, label=lab)
                ax.legend(fontsize=6, frameon=False, markerscale=3)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(color)
            _save(fig, figdir / fname)

    catalog_131 = int((adata.obs["dataset"] == "GSE131907").sum())
    # catalog totals from uns if present; else analyzed
    cat131 = 24784  # author malignant catalog from winning-pair meta
    cat205 = 28912  # will overwrite from value_counts if we stored it
    # better: from extract inventory if present
    inv_path = args.input.with_suffix(".inventory.json")
    inv = json.loads(inv_path.read_text()) if inv_path.is_file() else {}

    S = {
        "accessions": ["GSE131907", "GSE205335"],
        "winning_pair": True,
        "gse148071_added": False,
        "gse207422_added": False,
        "dual_high": False,
        "primary_gene": "CLDN4",
        "unit": "patient (GSE131907 Sample / GSE205335 patient)",
        "cap_per_unit": CAP_NOTE,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_units": int(patient.shape[0]),
        "n_units_gse131907": int((patient["dataset"] == "GSE131907").sum()),
        "n_units_gse205335": int((patient["dataset"] == "GSE205335").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired": int(len(paired_df)),
        "catalog_gse131907_malignant": 24784,
        "catalog_gse205335_malignant": int(inv.get("by_dataset", {}).get("GSE205335", adata.obs["dataset"].eq("GSE205335").sum())),
        "histology_counts": adata.obs["histology"].astype(str).value_counts().to_dict(),
        "potency_category_counts": adata.obs["CytoTRACE2_Potency"].astype(str).value_counts().to_dict(),
        "cldn4_tertile_counts": adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict(),
        "genes_absent": absent,
        "potency_method": ct2_info,
        "gulati_ran": True,
        "primary": primary,
        "sensitivity": sensitivity,
        "paired": paired,
        "extra_emitted": bool(extra_emitted),
        "inventory": inv,
    }
    (outdir / "summary.json").write_text(json.dumps(S, indent=2, default=str))
    write_finding(args.finding, S)
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_units": int(patient.shape[0]),
                "n_eligible": int(len(elig)),
                "potency": potency_name,
                "cytotrace2_ran": bool(ct2_info.get("ran") and potency_name == "CytoTRACE2"),
                "done_table": str(tabdir / "patient_cldn4_vs_potency.tsv"),
                "finding": str(args.finding),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
