#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA vs malignant CLDN4 on GSE207422.

Additive to the prior TACSTD2 Milo folder (methods/scrna_milo). Same public
GEO UMI, same graph / SpatialFDR fallback, same marker lineages. The
sample-level covariate and composition tests are CLDN4, not TACSTD2.

Public only: GEO supplementary UMI + sample xlsx. Author barcode labels are
not deposited; lineages are reconstructed from Hu et al. 2023 canonical
markers. miloR / edgeR are not used.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
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


def wanted_markers() -> set[str]:
    genes = set(ALWAYS)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    return genes


def stream_pass1(matrix_path: Path, markers: set[str]):
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        names: list[str] = []
        means: list[float] = []
        varis: list[float] = []
        found: dict[str, np.ndarray] = {}
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} values, expected {n}")
            n_umi += arr
            names.append(gene)
            mu = float(arr.mean())
            means.append(mu)
            varis.append(float(arr.var()) if n else 0.0)
            if gene in markers:
                found[gene] = arr
            if i % 4000 == 0:
                print(f"  pass1 {i} genes, markers={len(found)}", flush=True)
    return cell_ids, n_umi, names, np.asarray(means), np.asarray(varis), found


def stream_pass2(matrix_path: Path, keep: set[str], n_cells: int) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if gene not in keep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            out[gene] = arr
            if len(out) % 200 == 0:
                print(f"  pass2 stored {len(out)}/{len(keep)}", flush=True)
    missing = keep - set(out)
    if missing:
        print(f"  missing {len(missing)} requested genes (ok if absent)", flush=True)
    return out


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
    rec = {
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
    return rec


def write_finding(summary: dict, sample_tab: pd.DataFrame, out_path: Path) -> None:
    """Honest FINDING.md from the numbers that were just computed."""
    da = summary["da_malignant_cldn4"]
    da_mpr = summary["da_mpr_vs_nmpr"]
    g = summary["graph"]
    s = summary["samples_post"]
    c = summary["composition"]
    hn = summary["honest_n"]
    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = []
    for _, r in dropped.iterrows():
        drop_bits.append(f"{r['sample']}={int(r['n_malignant'])}")
    drop_txt = ", ".join(drop_bits) if drop_bits else "none"
    iface = c["interface_nhoods"]["spearman"]
    paired = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]
    w = paired["wilcoxon_high_vs_low"]
    med_h = w.get("median_a")
    med_l = w.get("median_b")
    med_h_s = "NA" if med_h is None else f"{med_h:.3f}"
    med_l_s = "NA" if med_l is None else f"{med_l:.3f}"
    p_w = "NA" if w.get("p") is None else f"{w['p']:.3g}"
    rho_i = iface.get("rho")
    p_i = iface.get("p")
    rho_i_s = "NA" if rho_i is None else f"{rho_i:.3f}"
    p_i_s = "NA" if p_i is None else f"{p_i:.3g}"
    min_sf = da.get("min_SpatialFDR")
    min_p = da.get("min_p")
    min_sf_s = "NA" if min_sf is None else f"{min_sf:.3f}"
    min_p_s = "NA" if min_p is None else f"{min_p:.2e}"
    min_sf_mpr = da_mpr.get("min_SpatialFDR")
    min_sf_mpr_s = "NA" if min_sf_mpr is None else f"{min_sf_mpr:.3f}"

    text = f"""# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE207422)

Additive to the prior TACSTD2 Milo run. Same public Hu et al. 2023 UMI, same
kNN graph / SpatialFDR fallback (not miloR). The covariate here is **sample
malignant CLDN4**, not TACSTD2.

The independent unit is the **post-treatment sample**, not the cell and not
the overlapping neighbourhood.

## n (sample is the unit)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant-like cells | **{s['n_with_malignant_ge10']}** of {s['n']} post samples | Mean log1p-CP10k CLDN4 in malignant-like cells. Dropped if <10 such cells ({drop_txt}). |
| MPR vs NMPR (Welch, secondary) | MPR / NMPR | **{s['MPR']} vs {s['NMPR']}** | Post-treatment only. pCR counted as MPR. Three pre-treatment biopsies excluded. |

Do not cite {g['cells_post']:,} post-treatment cells or {g['n_nhoods']:,} neighbourhoods as *n*. Neighbourhoods overlap. SpatialFDR is overlap-aware; a disjoint subset has n={hn['disjoint_nhood_subset_n']}.

Graph: {g['cells_all']:,} cells in the public matrix; {g['cells_post']:,} post-treatment cells; {g['n_nhoods']:,} neighbourhoods of size {g['nhood_size']}; k={g['k']}, d={g['d']}, {g['n_hvg']} HVG; PCA per-sample mean centering (not Harmony).

Malignant-like = epithelial AND NOT (alveolar/club/ciliated log1p-CP10k ≥ 1). Author CopyKAT IDs are not public. n malignant-like post = {summary['n_malignant_post']:,}; n T/NK post = {summary['n_tnk_post']:,}.

## SpatialFDR (all testable neighbourhoods)

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 | {da['n_testable']} | {da['n_p_lt_0.05']} | {min_p_s} | {min_sf_s} | **{da['n_SpatialFDR_lt_0.1']}** | {da['n_SpatialFDR_lt_0.05']} | {da['n_BH_FDR_lt_0.1']} |
| MPR vs NMPR | {da_mpr['n_testable']} | {da_mpr['n_p_lt_0.05']} | {da_mpr.get('min_p')} | {min_sf_mpr_s} | **{da_mpr['n_SpatialFDR_lt_0.1']}** | {da_mpr['n_SpatialFDR_lt_0.05']} | {da_mpr['n_BH_FDR_lt_0.1']} |

DA model = sample-level Spearman (CLDN4) or Welch t-test (MPR vs NMPR) on neighbourhood proportions. Not edgeR QLF. SpatialFDR = miloR `graphSpatialFDR` k-distance weights.

## Composition (transcriptional kNN, not histology)

Interface neighbourhoods (≥3 malignant-like and ≥3 T/NK cells): n={c['interface_nhoods']['n']}; Spearman malignant CLDN4 vs T/NK fraction ρ={rho_i_s}, p={p_i_s}, n={iface.get('n')}.

Disjoint interface subset: n={c['disjoint_interface']['n_disjoint_interface']} of {c['disjoint_interface']['n_disjoint_nhoods']} disjoint neighbourhoods; ρ={c['disjoint_interface']['spearman'].get('rho')}, p={c['disjoint_interface']['spearman'].get('p')}.

Sample-paired T/NK fraction in CLDN4-high vs CLDN4-low neighbourhoods (median split of neighbourhoods with ≥5 malignant-like cells; unit = sample): n={w.get('n')}; median T/NK high={med_h_s}, low={med_l_s}; Wilcoxon p={p_w}. High nhoods={paired['n_high_nhoods']}, low nhoods={paired['n_low_nhoods']}.

Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry (epithelium sits with epithelium). The sample-paired test is the only composition test whose independent unit is the patient.

## What this records

- Prior TACSTD2 Milo on this matrix is taken as given and is not re-run here.
- Numbers above are the CLDN4 neighbourhood tests that finished on the public UMI.
- miloR was not available; p-values are the documented sample-level fallback.
- GSE241934 was not run.

See `results/GSE207422/summary.json`, `da_malignant_cldn4.tsv`, `nhoods.tsv`, and `sample_scores.tsv`.
"""
    out_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    ap.add_argument("--metadata", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/scrna_milo_cldn4/results/GSE207422"))
    ap.add_argument("--finding", type=Path, default=Path("methods/scrna_milo_cldn4/FINDING.md"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    self_test()

    meta = pd.read_excel(args.metadata).dropna(subset=["Sample", "Patient"]).copy()
    meta["response_paper"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    markers = wanted_markers()
    print("pass1: gene stats + markers", flush=True)
    cell_ids, n_umi, gene_names, gene_mean, gene_var, marker_expr = stream_pass1(args.matrix, markers)
    n = len(cell_ids)
    print(f"cells={n} genes={len(gene_names)} median_nUMI={np.median(n_umi):.0f}", flush=True)

    hvg_idx = select_hvg(gene_mean, gene_var, n_hvg=args.n_hvg)
    hvg_names = [gene_names[i] for i in hvg_idx]
    keep = set(hvg_names) | markers
    print(f"pass2: {len(hvg_names)} HVG + {len(markers)} markers", flush=True)
    expr = stream_pass2(args.matrix, keep, n)
    for g, arr in marker_expr.items():
        expr.setdefault(g, arr)

    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    lineage = assign_lineage(log_cp, n)
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]])
    club = np.maximum(gene_log(log_cp, "SCGB1A1", n), gene_log(log_cp, "SCGB3A2", n))
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["TPPP3", "FOXJ1", "CAPS"]])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_epi = lineage == "Epithelial"
    is_malig = is_epi & (~clear_normal)
    is_tnk = np.isin(lineage, ["T", "NK"])
    cldn4 = log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))
    if "CLDN4" not in log_cp:
        raise RuntimeError("CLDN4 is absent from the public UMI matrix")

    sample = np.array([c.rsplit("_", 1)[0] for c in cell_ids])
    smap = meta.set_index("Sample")
    response = np.array([smap.loc[s, "response_paper"] if s in smap.index else "NA" for s in sample], dtype=object)
    is_post = np.array([bool(smap.loc[s, "is_post"]) if s in smap.index else False for s in sample])

    keep_cells = is_post
    print(f"post-treatment cells: {keep_cells.sum()} / {n}", flush=True)
    sample_k = sample[keep_cells]
    lineage_k = lineage[keep_cells]
    is_malig_k = is_malig[keep_cells]
    is_tnk_k = is_tnk[keep_cells]
    cldn_k = cldn4[keep_cells]
    n_k = int(keep_cells.sum())

    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g][keep_cells] * scale[keep_cells]) for g in hvg_present]).T.astype(np.float32)
    print(f"HVG matrix {X.shape}; PCA+kNN k={args.k} d={args.d}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample_k, n_pcs=args.d, k=args.k, random_state=args.seed)
    indices = refine_indices(pcs, knn_idx, prop=args.prop, random_state=args.seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, n_k)
    print(f"neighbourhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample_k))
    counts = count_matrix(members, sample_k, sample_levels)
    sizes = sample_sizes(sample_k, sample_levels)
    group = np.array([str(smap.loc[s, "response_paper"]) for s in sample_levels], dtype=object)

    mal_score = []
    tnk_frac_sample = []
    n_mal_s = []
    for s in sample_levels:
        m = sample_k == s
        n_mal = int((m & is_malig_k).sum())
        n_mal_s.append(n_mal)
        mal_score.append(float(cldn_k[m & is_malig_k].mean()) if n_mal >= 10 else np.nan)
        tnk_frac_sample.append(float(is_tnk_k[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac_sample = np.asarray(tnk_frac_sample, dtype=float)

    da_mpr = welch_da(counts, sizes, group, group_a="MPR", group_b="NMPR", min_per_group=2, min_samples=4)
    da_mpr = attach_fdr(da_mpr, nhoods.k_distance)
    da_cldn = spearman_da(counts, sizes, mal_score, min_samples=5)
    da_cldn = attach_fdr(da_cldn, nhoods.k_distance)

    # knn_nhood.nhood_composition names the gene columns tacstd2_*; rename to cldn4_*.
    comp = nhood_composition(members, is_tnk_k, is_malig_k, cldn_k, lineage_k)
    comp = comp.rename(columns={"tacstd2_all_mean": "cldn4_all_mean", "tacstd2_malig_mean": "cldn4_malig_mean"})
    merged = comp.merge(da_mpr, on="nhood", suffixes=("", "_mpr")).merge(
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
    da_mpr.to_csv(args.outdir / "da_mpr_vs_nmpr.tsv", sep="\t", index=False)
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
    paired = sample_paired_tnk_by_tacstd2(members, sample_k, is_tnk_k, comp["cldn4_malig_mean"].to_numpy(), high, low)
    paired = paired.merge(
        pd.DataFrame(
            {
                "sample": sample_levels,
                "response": group,
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

    cldn_ok = merged["cldn4_testable"].fillna(False).to_numpy()
    da_vs_comp = spearman_safe(merged.loc[cldn_ok, "cldn4_rho"], merged.loc[cldn_ok, "frac_tnk"])

    sample_tab = pd.DataFrame(
        {
            "sample": sample_levels,
            "patient": [str(smap.loc[s, "Patient"]) for s in sample_levels],
            "response": group,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "malignant_cldn4": mal_score,
            "frac_tnk": tnk_frac_sample,
        }
    )
    sample_tab.to_csv(args.outdir / "sample_scores.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "public_only": True,
        "author_barcode_labels": False,
        "lineage_source": "marker_argmax_Hu_canonical",
        "malignant_definition": "epithelial AND NOT (alveolar/club/ciliated log1p-CP10k >= 1); not CopyKAT",
        "gene": "CLDN4",
        "additive_to": "methods/scrna_milo TACSTD2 Milo (taken as given)",
        "miloR": False,
        "da_model": "sample-level Welch t-test (MPR vs NMPR) or Spearman (malignant CLDN4); not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "graph": {
            "cells_all": int(n),
            "cells_post": int(n_k),
            "k": args.k,
            "d": min(args.d, X.shape[1], n_k - 1),
            "prop": args.prop,
            "n_hvg": len(hvg_present),
            "n_nhoods": int(len(members)),
            "nhood_size": args.k + 1,
            "batch": "PCA per-sample mean centering (not Harmony)",
        },
        "samples_post": {
            "n": int(len(sample_levels)),
            "MPR": int((group == "MPR").sum()),
            "NMPR": int((group == "NMPR").sum()),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
        },
        "lineage_counts_post": {k: int(v) for k, v in pd.Series(lineage_k).value_counts().items()},
        "n_malignant_post": int(is_malig_k.sum()),
        "n_tnk_post": int(is_tnk_k.sum()),
        "da_mpr_vs_nmpr": fdr_counts(da_mpr),
        "da_malignant_cldn4": fdr_counts(da_cldn),
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
            "nhood_cldn4_DA_rho_vs_fracTNK": da_vs_comp,
        },
        "honest_n": {
            "independent_unit_DA": f"post-treatment samples (n={len(sample_levels)}; MPR={int((group == 'MPR').sum())} including pCR, NMPR={int((group == 'NMPR').sum())})",
            "n_samples_with_malignant_CLDN4_score": int(np.isfinite(mal_score).sum()),
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(n_k),
        },
    }
    write_json(args.outdir / "summary.json", summary)
    write_finding(summary, sample_tab, args.finding)

    plot_volcano(
        da_mpr,
        "logFC_B_minus_A",
        "GSE207422 neighbourhood DA: NMPR vs MPR",
        args.outdir / "fig_da_mpr_volcano.png",
        "log2 FC (NMPR − MPR sample proportion)",
    )
    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "logFC_B_minus_A"}),
        "logFC_B_minus_A",
        "GSE207422 neighbourhood DA vs malignant CLDN4",
        args.outdir / "fig_da_cldn4_volcano.png",
        "Spearman ρ (nhood abundance vs sample malignant CLDN4)",
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
                color = "#c44e52" if r["response"] == "NMPR" else "#4c72b0"
                ax.plot([0, 1], [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]], "-o", c=color, alpha=0.75, ms=5)
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title("GSE207422 sample-paired (unit = patient)")
        fig.tight_layout()
        fig.savefig(args.outdir / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)

    print(
        json.dumps(
            {
                "n_nhoods": summary["graph"]["n_nhoods"],
                "da_mpr": summary["da_mpr_vs_nmpr"],
                "da_cldn4": summary["da_malignant_cldn4"],
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
