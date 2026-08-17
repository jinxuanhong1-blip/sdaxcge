#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA vs malignant CLDN4 on GSE148071.

Additive to prior GSE207422 TACSTD2 / CLDN4 Milo folders. Public GEO
per-sample raw-count TSVs only (Wu et al., Nat Commun 2021). Author
barcode labels and CNA IDs are not deposited; lineages are reconstructed
from canonical markers. miloR / edgeR are not used.

Independent unit = patient biopsy (one GEO sample per patient).
Do not cite cell count as n.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from knn_nhood import (  # noqa: E402
    EPS,
    attach_fdr,
    count_matrix,
    greedy_independent,
    make_nhoods,
    nhood_composition,
    nhood_membership,
    pca_knn,
    refine_indices,
    sample_paired_tnk_by_tacstd2,
    sample_sizes,
    select_hvg,
    self_test,
    spearman_da,
    spearman_safe,
    welch_da,
    wilcoxon_paired,
    write_json,
)

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
NORMAL_LUNG = ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
ALWAYS = ["TACSTD2", "CLDN4", "PTPRC", "EPCAM", "CD3D", "NKG7"]
HEADER_GENE = {"", "gene", "Gene", "GENE", "index", "Index", "symbol", "Symbol", "SYMBOL"}


def wanted_markers() -> set[str]:
    genes = set(ALWAYS)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    return genes


def norm_gene(raw: str) -> str:
    gene = raw.strip().strip('"')
    if "|" in gene:
        gene = gene.split("|")[-1]
    if ";" in gene:
        gene = gene.split(";")[0]
    return gene.split(".")[0]


def patient_from_name(path: Path) -> str:
    # GSM4453576_P1_exp.txt.gz → P1
    stem = path.name.replace(".txt.gz", "").replace(".txt", "")
    parts = stem.split("_")
    for p in parts:
        if p.startswith("P") and p[1:].isdigit():
            return p
    return parts[1] if len(parts) > 1 else stem


def parse_header(line: str) -> list[str]:
    header = line.rstrip("\n").split("\t")
    if header and header[0] in HEADER_GENE:
        return header[1:]
    return header


def list_exp_files(datadir: Path) -> list[Path]:
    files = sorted(datadir.glob("*_exp.txt.gz"))
    if not files:
        files = sorted(datadir.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz under {datadir}")
    return files


def stream_file(path: Path, markers: set[str], keep: set[str] | None):
    """Yield (gene, array) for every gene, and return cell_ids.

    If keep is set, only those genes are returned in `found`; all genes still
    contribute to n_umi / sum / sumsq when keep is None (pass 1).
    """
    with gzip.open(path, "rt") as handle:
        cell_ids = parse_header(handle.readline())
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        found: dict[str, np.ndarray] = {}
        gene_sum: dict[str, float] = {}
        gene_sumsq: dict[str, float] = {}
        n_genes = 0
        for line in handle:
            gene_raw, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = norm_gene(gene_raw)
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            n_genes += 1
            if keep is None:
                gene_sum[gene] = gene_sum.get(gene, 0.0) + float(arr.sum())
                gene_sumsq[gene] = gene_sumsq.get(gene, 0.0) + float(np.dot(arr, arr))
                if gene in markers:
                    found[gene] = arr
            elif gene in keep:
                found[gene] = arr
        return cell_ids, n_umi, found, gene_sum, gene_sumsq, n_genes


def assign_lineage(log_cp: dict[str, np.ndarray], n: int) -> np.ndarray:
    scores = {}
    for name, genes in LINEAGE_MARKERS.items():
        mats = [log_cp[g] for g in genes if g in log_cp]
        scores[name] = np.mean(np.vstack(mats), axis=0) if mats else np.zeros(n, dtype=np.float32)
    names = list(scores)
    mat = np.vstack([scores[n_] for n_ in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(n)]
    labels = np.array(names, dtype=object)[best]
    cd3 = log_cp.get("CD3D", np.zeros(n, dtype=np.float32))
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels = labels.copy()
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def gene_log(log_cp: dict[str, np.ndarray], name: str, n: int) -> np.ndarray:
    return log_cp[name] if name in log_cp else np.zeros(n, dtype=np.float32)


def plot_volcano(df: pd.DataFrame, x: str, title: str, path: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    testable = df["testable"].to_numpy()
    ax.scatter(
        df.loc[~testable, x],
        -np.log10(np.clip(df.loc[~testable, "p"], EPS, 1)),
        s=6,
        c="#bbbbbb",
        alpha=0.4,
        label="not testable",
    )
    t = df[testable]
    sig = t["SpatialFDR"] < 0.1
    ax.scatter(
        t.loc[~sig, x],
        -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)),
        s=8,
        c="#4c72b0",
        alpha=0.55,
        label="testable SpatialFDR≥0.1",
    )
    ax.scatter(
        t.loc[sig, x],
        -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)),
        s=14,
        c="#c44e52",
        alpha=0.85,
        label="SpatialFDR<0.1",
    )
    ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"−log10 p (sample-level test)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(x, y, s=8, c="#4c72b0", alpha=0.45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]]
    return {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((t["p"] < 0.05).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()) if len(t) else 0,
        "min_p": float(t["p"].min()) if len(t) else None,
        "min_SpatialFDR": float(t["SpatialFDR"].min()) if len(t) else None,
        "min_BH_FDR": float(t["BH_FDR"].min()) if len(t) else None,
    }


def _fmt(x, nd=3, sci=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.2e}" if sci else f"{x:.{nd}f}"


def write_finding(summary: dict, sample_tab: pd.DataFrame, out_path: Path) -> None:
    da = summary["da_malignant_cldn4"]
    da_hl = summary["da_cldn4_high_vs_low"]
    g = summary["graph"]
    s = summary["samples"]
    c = summary["composition"]
    hn = summary["honest_n"]
    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = [f"{r['sample']}={int(r['n_malignant'])}" for _, r in dropped.iterrows()]
    drop_txt = ", ".join(drop_bits) if drop_bits else "none"
    iface = c["interface_nhoods"]["spearman"]
    paired = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]
    w = paired["wilcoxon_high_vs_low"]
    text = f"""# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE148071)

Additive neighbourhood DA on the public Wu et al. 2021 advanced-NSCLC
biopsies (GSE148071). Prior TACSTD2 / CLDN4 Milo on GSE207422 is taken as
given and is not re-run. The covariate here is **sample malignant CLDN4**.

The independent unit is the **patient biopsy**, not the cell and not the
overlapping neighbourhood. miloR is not used.

## n (sample is the unit)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant-like cells | **{s['n_with_malignant_ge10']}** of {s['n']} GEO samples | Mean log1p-CP10k CLDN4 in malignant-like cells. Dropped if <10 such cells ({drop_txt}). |
| CLDN4 high vs low (Welch, secondary) | high / low | **{s['n_cldn4_high']} vs {s['n_cldn4_low']}** | Median split of the {s['n_with_malignant_ge10']} samples that have a malignant CLDN4 score. |

Do not cite {g['n_cells']:,} cells or {g['n_nhoods']:,} neighbourhoods as *n*. Neighbourhoods overlap. SpatialFDR is overlap-aware; a disjoint subset has n={hn['disjoint_nhood_subset_n']}.

Graph: {g['n_cells']:,} cells in the public matrices ({s['n']} samples); {g['n_nhoods']:,} neighbourhoods of size {g['nhood_size']}; k={g['k']}, d={g['d']}, {g['n_hvg']} HVG; PCA per-sample mean centering (not Harmony). Paper QC n after filtering was 90,406 cells / 42 patients; this run uses every barcode in the GEO count files after dropping empty libraries (UMI=0).

Malignant-like = epithelial AND NOT (alveolar/club/ciliated log1p-CP10k ≥ 1). Author inferCNV / CNA IDs are not public. n malignant-like = {summary['n_malignant']:,}; n T/NK = {summary['n_tnk']:,}.

GEO SOFT has age/sex only. Histology (LUAD/LUSC) lives in the paper supplement and is **not** used as a DA covariate. This series is advanced diagnostic biopsies, not an ICI-response cohort.

## SpatialFDR (all testable neighbourhoods)

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 | {da['n_testable']} | {da['n_p_lt_0.05']} | {_fmt(da.get('min_p'), sci=True)} | {_fmt(da.get('min_SpatialFDR'))} | **{da['n_SpatialFDR_lt_0.1']}** | {da['n_SpatialFDR_lt_0.05']} | {da['n_BH_FDR_lt_0.1']} |
| CLDN4 high vs low | {da_hl['n_testable']} | {da_hl['n_p_lt_0.05']} | {_fmt(da_hl.get('min_p'), sci=True)} | {_fmt(da_hl.get('min_SpatialFDR'))} | **{da_hl['n_SpatialFDR_lt_0.1']}** | {da_hl['n_SpatialFDR_lt_0.05']} | {da_hl['n_BH_FDR_lt_0.1']} |

DA model = sample-level Spearman (CLDN4) or Welch t-test (high vs low) on neighbourhood proportions. Not edgeR QLF. SpatialFDR = miloR `graphSpatialFDR` k-distance weights.

## Composition (transcriptional kNN, not histology)

Interface neighbourhoods (≥3 malignant-like and ≥3 T/NK cells): n={c['interface_nhoods']['n']}; Spearman malignant CLDN4 vs T/NK fraction ρ={_fmt(iface.get('rho'))}, p={_fmt(iface.get('p'), sci=True)}, n={iface.get('n')}.

Disjoint interface subset: n={c['disjoint_interface']['n_disjoint_interface']} of {c['disjoint_interface']['n_disjoint_nhoods']} disjoint neighbourhoods; ρ={_fmt(c['disjoint_interface']['spearman'].get('rho'))}, p={_fmt(c['disjoint_interface']['spearman'].get('p'), sci=True)}.

Sample-paired T/NK fraction in CLDN4-high vs CLDN4-low neighbourhoods (median split of neighbourhoods with ≥5 malignant-like cells; unit = sample): n={w.get('n')}; median T/NK high={_fmt(w.get('median_a'))}, low={_fmt(w.get('median_b'))}; Wilcoxon p={_fmt(w.get('p'), sci=True)}. High nhoods={paired['n_high_nhoods']}, low nhoods={paired['n_low_nhoods']}.

Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry (epithelium sits with epithelium). The sample-paired test is the only composition test whose independent unit is the patient.

## What this records

- This is an additive GSE148071 cut. GSE207422 Milo is not re-run.
- Numbers above are the CLDN4 neighbourhood tests that finished on the public counts.
- miloR was not available; p-values are the documented sample-level fallback.
- Cell count is not *n*. *n* is the number of patient biopsies that enter each test.

See `results/GSE148071/summary.json`, `da_malignant_cldn4.tsv`, `nhoods.tsv`, and `sample_scores.tsv`.
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE148071/files"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/gse148071_milo_cldn4/results/GSE148071"))
    ap.add_argument("--finding", type=Path, default=Path("methods/gse148071_milo_cldn4/FINDING.md"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-umi", type=float, default=1.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    self_test()

    files = list_exp_files(args.datadir)
    markers = wanted_markers()
    print(f"pass1: {len(files)} files", flush=True)

    cell_ids_all: list[str] = []
    sample_all: list[str] = []
    n_umi_all: list[np.ndarray] = []
    marker_blocks: dict[str, list[np.ndarray]] = defaultdict(list)
    gene_sum: dict[str, float] = {}
    gene_sumsq: dict[str, float] = {}
    n_cells_total = 0
    file_n: list[tuple[str, int, int]] = []

    for fp in files:
        patient = patient_from_name(fp)
        cell_ids, n_umi, found, gsum, gss, n_genes = stream_file(fp, markers, keep=None)
        n = len(cell_ids)
        cell_ids_all.extend(cell_ids)
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        n_cells_total += n
        file_n.append((patient, n, n_genes))
        for g, arr in found.items():
            marker_blocks[g].append(arr)
        for g, val in gsum.items():
            gene_sum[g] = gene_sum.get(g, 0.0) + val
            gene_sumsq[g] = gene_sumsq.get(g, 0.0) + gss[g]
        print(f"  {fp.name} {patient} n={n} genes={n_genes} markers={sorted(found)}", flush=True)

    n_umi = np.concatenate(n_umi_all)
    sample = np.array(sample_all, dtype=object)
    keep_cell = n_umi >= args.min_umi
    n_drop_empty = int((~keep_cell).sum())
    print(f"cells={n_cells_total} empty_or_below_min_umi={n_drop_empty}", flush=True)

    genes = sorted(gene_sum)
    means = np.array([gene_sum[g] / n_cells_total for g in genes], dtype=float)
    varis = np.array([gene_sumsq[g] / n_cells_total - means[i] ** 2 for i, g in enumerate(genes)], dtype=float)
    varis = np.maximum(varis, 0.0)
    hvg_idx = select_hvg(means, varis, n_hvg=args.n_hvg)
    hvg_names = [genes[i] for i in hvg_idx]
    keep_genes = set(hvg_names) | markers
    print(f"HVG={len(hvg_names)} keep_genes={len(keep_genes)}", flush=True)

    print("pass2: load HVG + markers", flush=True)
    expr_blocks: dict[str, list[np.ndarray]] = defaultdict(list)
    for fp in files:
        patient = patient_from_name(fp)
        cell_ids, _, found, _, _, _ = stream_file(fp, markers, keep=keep_genes)
        n = len(cell_ids)
        present = set(found)
        for g in keep_genes:
            expr_blocks[g].append(found[g] if g in present else np.zeros(n, dtype=np.float32))
        print(f"  {fp.name} stored {len(present)}/{len(keep_genes)}", flush=True)

    expr = {g: np.concatenate(blocks) for g, blocks in expr_blocks.items()}
    n = int(n_umi.size)
    scale = 1e4 / np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    lineage = assign_lineage(log_cp, n)
    if "CLDN4" not in log_cp:
        raise RuntimeError("CLDN4 is absent from the public GSE148071 count files")
    cldn4 = log_cp["CLDN4"]
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]])
    club = np.maximum(gene_log(log_cp, "SCGB1A1", n), gene_log(log_cp, "SCGB3A2", n))
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["TPPP3", "FOXJ1", "CAPS"]])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_epi = lineage == "Epithelial"
    is_malig = is_epi & (~clear_normal)
    is_tnk = np.isin(lineage, ["T", "NK"])

    if keep_cell.sum() < n:
        sample = sample[keep_cell]
        lineage = lineage[keep_cell]
        is_malig = is_malig[keep_cell]
        is_tnk = is_tnk[keep_cell]
        cldn4 = cldn4[keep_cell]
        n_umi = n_umi[keep_cell]
        scale = scale[keep_cell]
        expr = {g: v[keep_cell] for g, v in expr.items()}
        n = int(keep_cell.sum())

    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g] * scale) for g in hvg_present]).T.astype(np.float32)
    print(f"HVG matrix {X.shape}; PCA+kNN k={args.k} d={args.d}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample, n_pcs=args.d, k=args.k, random_state=args.seed)
    indices = refine_indices(pcs, knn_idx, prop=args.prop, random_state=args.seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, n)
    print(f"neighbourhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample), key=lambda s: (len(s), s))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)

    mal_score = []
    tnk_frac_sample = []
    n_mal_s = []
    for s in sample_levels:
        m = sample == s
        n_mal = int((m & is_malig).sum())
        n_mal_s.append(n_mal)
        mal_score.append(float(cldn4[m & is_malig].mean()) if n_mal >= 10 else np.nan)
        tnk_frac_sample.append(float(is_tnk[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac_sample = np.asarray(tnk_frac_sample, dtype=float)

    scored = np.isfinite(mal_score)
    group = np.array(["unscored"] * len(sample_levels), dtype=object)
    if scored.sum() >= 2:
        med_s = float(np.median(mal_score[scored]))
        group[scored & (mal_score >= med_s)] = "high"
        group[scored & (mal_score < med_s)] = "low"
    else:
        med_s = float("nan")

    da_hl = welch_da(counts, sizes, group, group_a="low", group_b="high", min_per_group=2, min_samples=4)
    da_hl = attach_fdr(da_hl, nhoods.k_distance)
    da_cldn = spearman_da(counts, sizes, mal_score, min_samples=5)
    da_cldn = attach_fdr(da_cldn, nhoods.k_distance)

    comp = nhood_composition(members, is_tnk, is_malig, cldn4, lineage)
    comp = comp.rename(columns={"tacstd2_all_mean": "cldn4_all_mean", "tacstd2_malig_mean": "cldn4_malig_mean"})
    merged = comp.merge(da_hl, on="nhood", suffixes=("", "_hl")).merge(
        da_cldn[["nhood", "spearman_rho", "p", "testable", "BH_FDR", "SpatialFDR"]].rename(
            columns={
                "spearman_rho": "cldn4_rho",
                "p": "cldn4_p",
                "testable": "cldn4_testable",
                "BH_FDR": "cldn4_BH_FDR",
                "SpatialFDR": "cldn4_SpatialFDR",
            }
        ),
        on="nhood",
    )
    merged.to_csv(args.outdir / "nhoods.tsv", sep="\t", index=False)
    da_hl.to_csv(args.outdir / "da_cldn4_high_vs_low.tsv", sep="\t", index=False)
    da_cldn.to_csv(args.outdir / "da_malignant_cldn4.tsv", sep="\t", index=False)

    iface = (comp["n_malig"] >= args.min_interface) & (comp["n_tnk"] >= args.min_interface)
    all_sp = spearman_safe(comp["cldn4_malig_mean"], comp["frac_tnk"])
    iface_sp = spearman_safe(comp.loc[iface, "cldn4_malig_mean"], comp.loc[iface, "frac_tnk"])
    indep = greedy_independent(members, max_shared=0)
    indep_mask = np.zeros(len(members), dtype=bool)
    indep_mask[indep] = True
    indep_sp = spearman_safe(
        comp.loc[indep_mask & iface, "cldn4_malig_mean"],
        comp.loc[indep_mask & iface, "frac_tnk"],
    )
    has_mal = comp["n_malig"] >= 5
    med = float(np.nanmedian(comp.loc[has_mal, "cldn4_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() >= med)
    low = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() < med)
    paired = sample_paired_tnk_by_tacstd2(members, sample, is_tnk, comp["cldn4_malig_mean"].to_numpy(), high, low)
    paired = paired.merge(
        pd.DataFrame(
            {
                "sample": sample_levels,
                "cldn4_arm": group,
                "malignant_cldn4": mal_score,
                "sample_frac_tnk": tnk_frac_sample,
                "n_malignant": n_mal_s,
                "n_cells": sizes.astype(int),
            }
        ),
        on="sample",
    )
    paired.to_csv(args.outdir / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])

    sample_tab = pd.DataFrame(
        {
            "sample": sample_levels,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "malignant_cldn4": mal_score,
            "cldn4_arm": group,
            "frac_tnk": tnk_frac_sample,
        }
    )
    sample_tab.to_csv(args.outdir / "sample_scores.tsv", sep="\t", index=False)
    pd.DataFrame(file_n, columns=["sample", "n_barcodes_file", "n_genes_file"]).to_csv(
        args.outdir / "file_inventory.tsv", sep="\t", index=False
    )

    summary = {
        "dataset": "GSE148071",
        "citation": "Wu et al. Nature Communications 2021 PMID 33953163",
        "public_only": True,
        "author_barcode_labels": False,
        "lineage_source": "marker_argmax_canonical",
        "malignant_definition": "epithelial AND NOT (alveolar/club/ciliated log1p-CP10k >= 1); not inferCNV",
        "gene": "CLDN4",
        "additive_to": "methods/scrna_milo and methods/scrna_milo_cldn4 on GSE207422 (taken as given)",
        "miloR": False,
        "da_model": "sample-level Spearman (malignant CLDN4) or Welch t-test (median-split high vs low); not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "graph": {
            "n_cells": int(n),
            "n_cells_before_umi_filter": int(n_cells_total),
            "n_dropped_low_umi": int(n_drop_empty),
            "k": args.k,
            "d": min(args.d, X.shape[1], n - 1),
            "prop": args.prop,
            "n_hvg": len(hvg_present),
            "n_nhoods": int(len(members)),
            "nhood_size": args.k + 1,
            "median_nhood_size": float(np.median([m.size for m in members])),
            "batch": "PCA per-sample mean centering (not Harmony)",
        },
        "samples": {
            "n": int(len(sample_levels)),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
            "n_cldn4_high": int((group == "high").sum()),
            "n_cldn4_low": int((group == "low").sum()),
            "median_malignant_cldn4_cut": med_s,
        },
        "lineage_counts": {k: int(v) for k, v in pd.Series(lineage).value_counts().items()},
        "n_malignant": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "da_malignant_cldn4": fdr_counts(da_cldn),
        "da_cldn4_high_vs_low": fdr_counts(da_hl),
        "composition": {
            "note": "transcriptional kNN != spatial niche; unrestricted CLDN4 vs T/NK is partly lineage geometry",
            "all_nhoods_maligCLDN4_vs_fracTNK": all_sp,
            "interface_nhoods": {
                "definition": f">={args.min_interface} malignant AND >={args.min_interface} T/NK cells",
                "n": int(iface.sum()),
                "spearman": iface_sp,
            },
            "disjoint_interface": {
                "n_disjoint_nhoods": int(indep.size),
                "n_disjoint_interface": int((indep_mask & iface.to_numpy()).sum()),
                "spearman": indep_sp,
            },
            "sample_paired_TNK_in_CLDN4high_vs_low_nhoods": {
                "high_low_cut": "median malignant CLDN4 among nhoods with >=5 malignant cells",
                "median_cut": med,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
            },
        },
        "honest_n": {
            "independent_unit_DA": f"patient biopsies (n={len(sample_levels)} GEO samples; {int(np.isfinite(mal_score).sum())} with ≥10 malignant-like cells)",
            "n_samples_with_malignant_CLDN4_score": int(np.isfinite(mal_score).sum()),
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(n),
        },
        "file_inventory": [{"sample": a, "n_barcodes": b, "n_genes": c} for a, b, c in file_n],
    }
    write_json(args.outdir / "summary.json", summary)
    write_finding(summary, sample_tab, args.finding)

    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "logFC_B_minus_A"}),
        "logFC_B_minus_A",
        "GSE148071 neighbourhood DA vs malignant CLDN4",
        args.outdir / "fig_da_cldn4_volcano.png",
        "Spearman ρ (nhood abundance vs sample malignant CLDN4)",
    )
    plot_volcano(
        da_hl,
        "logFC_B_minus_A",
        "GSE148071 neighbourhood DA: CLDN4-high vs low samples",
        args.outdir / "fig_da_cldn4_high_vs_low_volcano.png",
        "log2 FC (high − low sample proportion)",
    )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        args.outdir / "fig_interface_cldn4_vs_tnk.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"Interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                color = "#c44e52" if r["cldn4_arm"] == "high" else "#4c72b0"
                ax.plot([0, 1], [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]], "-o", c=color, alpha=0.75, ms=5)
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title("GSE148071 sample-paired (unit = patient)")
        fig.tight_layout()
        fig.savefig(args.outdir / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)

    print(
        json.dumps(
            {
                "n_nhoods": summary["graph"]["n_nhoods"],
                "n_samples": summary["samples"],
                "da_cldn4": summary["da_malignant_cldn4"],
                "da_high_low": summary["da_cldn4_high_vs_low"],
                "iface": summary["composition"]["interface_nhoods"],
                "paired": summary["composition"]["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"],
                "finding": str(args.finding),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
