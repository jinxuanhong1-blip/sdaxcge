#!/usr/bin/env python3
"""Tumor (malignant) vs immune-cell expression of TACSTD2 (TROP2) and CLDN4
in open TISCH2 lung cancer scRNA-seq objects.

TISCH2 ships each dataset as
  <DS>_expression.h5           MAESTRO/10x-style sparse matrix, genes x cells,
                               values are log2(TPM/10 + 1)-normalised expression
  <DS>_CellMetainfo_table.tsv  uniform annotations incl. "Celltype (major-lineage)"
                               with an explicit "Malignant" label

For every gene of interest we compute, per dataset:
  * per-lineage n / mean / median / %-expressing
  * cell-level Mann-Whitney U (malignant vs pooled immune, and vs each immune
    lineage), with AUROC effect size and log2 fold change of means
  * per-patient pseudobulk means + paired Wilcoxon signed-rank test
    (malignant vs immune within the same patient), which is robust to the
    pseudoreplication problem of cell-level tests
  * figures: per-lineage violins and per-patient paired pseudobulk plot

Controls: EPCAM (epithelial/tumor positive control), PTPRC/CD45 (immune
positive control -> confirms the direction of the pipeline is not an artifact).

Usage:
  python 02_trop2_cldn4_tumor_vs_immune.py --dataset NSCLC_GSE131907 \
      --data-dir data_fable_tisch --out-root results/fable_tisch
"""

from __future__ import annotations

import argparse
import json
import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

GENES = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC"]
TARGET_GENES = ["TACSTD2", "CLDN4"]

# TISCH2 major-lineage vocabulary -> compartment
IMMUNE_LINEAGES = {
    "CD8T", "CD8Tex", "CD4Tconv", "Treg", "Tprolif", "TMKI67", "NK",
    "B", "Plasma", "Mono/Macro", "Monocyte", "Macrophage", "DC", "pDC",
    "Mast", "Neutrophils", "ILC", "Myeloid", "Tcell", "NKT", "Erythrocytes",
}
NON_IMMUNE_NON_MALIGNANT = {
    "Fibroblasts", "Myofibroblasts", "SMC", "Endothelial", "Epithelial",
    "AT1", "AT2", "Alveolar", "Oligodendrocyte", "OPC", "Others",
    "Melanocytes", "Hepatocytes", "Ciliated", "Club", "Pericytes",
    "Basal",
}
# Residual non-malignant lung epithelium (used for the epithelial-restriction
# contrast: malignant vs leftover epithelial vs immune).
NONMALIGNANT_EPITHELIAL = {
    "Epithelial", "Alveolar", "AT1", "AT2", "Basal", "Ciliated", "Club",
}
MIN_CELLS_PER_LINEAGE = 30      # drop tiny lineages from per-lineage stats
MIN_CELLS_PER_PATIENT = 20      # per compartment, for pseudobulk pairing


def load_tisch_h5(path: str) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Return (cells x genes CSR matrix, gene names, barcodes) from a TISCH h5."""
    with h5py.File(path, "r") as f:
        # find the group holding the sparse matrix
        grp = None
        if all(k in f for k in ("data", "indices", "indptr")):
            grp = f
        else:
            for key in f.keys():
                g = f[key]
                if isinstance(g, h5py.Group) and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix group found in {path}: {list(f.keys())}")

        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(grp["shape"][:]) if "shape" in grp else None

        def _names(candidates):
            for c in candidates:
                node = grp.get(c) or f.get(c)
                if node is not None:
                    arr = node[:] if isinstance(node, h5py.Dataset) else None
                    if arr is None and isinstance(node, h5py.Group):
                        for sub in ("name", "id", "gene_names"):
                            if sub in node:
                                arr = node[sub][:]
                                break
                    if arr is not None:
                        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
            return None

        genes = _names(["features", "gene_names", "genes", "rownames"])
        barcodes = _names(["barcodes", "cell_names", "colnames"])

    if shape is None:
        raise ValueError("h5 matrix has no shape dataset")

    # 10x/MAESTRO convention: CSC with shape (n_genes, n_cells); a CSC
    # genes x cells matrix is the same buffer as a CSR cells x genes matrix.
    n_genes, n_cells = int(shape[0]), int(shape[1])
    if genes is not None and len(genes) == n_cells and len(barcodes) == n_genes:
        # transposed convention
        n_genes, n_cells = n_cells, n_genes
    mat = sp.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))
    if genes is None or len(genes) != n_genes or barcodes is None or len(barcodes) != n_cells:
        raise ValueError(
            f"gene/barcode dims do not match matrix: genes={None if genes is None else len(genes)}, "
            f"barcodes={None if barcodes is None else len(barcodes)}, shape=({n_genes},{n_cells})"
        )
    return mat, genes, barcodes


def auroc_from_u(u: float, n1: int, n2: int) -> float:
    return float(u) / (float(n1) * float(n2))


def gene_vector(mat: sp.csr_matrix, genes: np.ndarray, gene: str) -> np.ndarray | None:
    idx = np.where(genes == gene)[0]
    if len(idx) == 0:
        return None
    return np.asarray(mat[:, idx[0]].todense()).ravel()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--data-dir", default="data_fable_tisch")
    ap.add_argument("--out-root", default="results/fable_tisch")
    args = ap.parse_args()

    ds = args.dataset
    out_dir = os.path.join(args.out_root, ds)
    os.makedirs(out_dir, exist_ok=True)

    # ---------------- load ----------------
    mat, genes, barcodes = load_tisch_h5(os.path.join(args.data_dir, f"{ds}_expression.h5"))
    meta = pd.read_csv(os.path.join(args.data_dir, f"{ds}_CellMetainfo_table.tsv"), sep="\t")
    meta = meta.set_index(meta.columns[0])
    lineage_col = [c for c in meta.columns if "major-lineage" in c][0]
    patient_col = next((c for c in meta.columns if c.lower() == "patient"), None)

    # align metadata to matrix barcodes
    meta = meta.reindex(barcodes)
    keep = meta[lineage_col].notna().to_numpy()
    mat = mat[keep]
    meta = meta.loc[keep]
    lineage = meta[lineage_col].astype(str).str.strip().to_numpy()

    # -- define the malignant compartment --------------------------------
    # Preferred: the explicit TISCH2 "Malignant" major lineage.
    # Fallback (e.g. NSCLC_GSE131907, where TISCH2 provides no malignancy
    # call): restrict the WHOLE analysis to tumor-tissue cells and use
    # tumor-tissue epithelial cells as a malignant proxy. In LUAD tumor
    # tissue the epithelial compartment is dominated by cancer cells, but we
    # label it transparently as a proxy.
    source_col = next((c for c in meta.columns if c.lower() in ("source", "tissue")), None)

    # Tumor-tissue vs normal-tissue epithelial (same lineage, different Source).
    # Computed on the un-restricted object so GSE131907 keeps its Normal lung
    # epithelial cells as a true non-malignant epithelial comparator.
    tumor_vs_normal_epithelial = {}
    if source_col is not None and "Epithelial" in set(lineage):
        src0 = meta[source_col].astype(str).str.lower().to_numpy()
        in_tumor0 = np.char.find(src0.astype(str), "tumor") >= 0
        in_normal0 = np.char.find(src0.astype(str), "normal") >= 0
        epi0 = lineage == "Epithelial"
        if (epi0 & in_tumor0).sum() >= MIN_CELLS_PER_LINEAGE and (epi0 & in_normal0).sum() >= MIN_CELLS_PER_LINEAGE:
            for g in GENES:
                v = gene_vector(mat, genes, g)
                if v is None:
                    continue
                a, b = v[epi0 & in_tumor0], v[epi0 & in_normal0]
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                eps = 1e-9
                tumor_vs_normal_epithelial[g] = {
                    "n_tumor_epithelial": int(len(a)),
                    "n_normal_epithelial": int(len(b)),
                    "mean_tumor_epithelial": float(a.mean()),
                    "mean_normal_epithelial": float(b.mean()),
                    "pct_pos_tumor_epithelial": float((a > 0).mean() * 100),
                    "pct_pos_normal_epithelial": float((b > 0).mean() * 100),
                    "log2FC_of_means": float(np.log2((a.mean() + eps) / (b.mean() + eps))),
                    "auroc": auroc_from_u(u, len(a), len(b)),
                    "p_value": float(p),
                }

    if "Malignant" in set(lineage):
        malignant_name = "Malignant"
        malignant_definition = 'TISCH2 major-lineage == "Malignant"'
    else:
        if source_col is None or "Epithelial" not in set(lineage):
            raise SystemExit(f"{ds}: no Malignant lineage and no usable proxy")
        src = meta[source_col].astype(str).str.lower().to_numpy()
        in_tumor = np.char.find(src.astype(str), "tumor") >= 0
        print(f"[info] no 'Malignant' lineage in {ds}; restricting to "
              f"{int(in_tumor.sum())} tumor-tissue cells ({source_col}) and using "
              f"tumor-tissue Epithelial as malignant proxy")
        mat = mat[in_tumor]
        meta = meta.loc[in_tumor]
        lineage = lineage[in_tumor]
        malignant_name = "Epithelial(tumor tissue)"
        lineage = np.where(lineage == "Epithelial", malignant_name, lineage)
        malignant_definition = (
            f'proxy: major-lineage == "Epithelial" restricted to tumor-tissue '
            f'cells ({source_col} contains "tumor"); no TISCH2 malignancy call '
            f"available for this dataset"
        )

    compartment = np.full(len(lineage), "other", dtype=object)
    compartment[np.isin(lineage, list(IMMUNE_LINEAGES))] = "immune"
    compartment[lineage == malignant_name] = "malignant"

    expr = {}
    for g in GENES:
        v = gene_vector(mat, genes, g)
        if v is None:
            print(f"[warn] {g} not found in {ds}")
        else:
            expr[g] = v

    unknown = sorted(set(lineage) - IMMUNE_LINEAGES - NON_IMMUNE_NON_MALIGNANT - {"Malignant", malignant_name})
    if unknown:
        print(f"[warn] unclassified lineages treated as 'other': {unknown}")

    n_mal = int((compartment == "malignant").sum())
    n_imm = int((compartment == "immune").sum())
    print(f"{ds}: {len(lineage)} cells, malignant={n_mal}, immune={n_imm}")

    # ---------------- per-lineage summary ----------------
    rows = []
    for g, v in expr.items():
        for lin in sorted(set(lineage)):
            m = lineage == lin
            if m.sum() < MIN_CELLS_PER_LINEAGE:
                continue
            comp = "malignant" if lin == malignant_name else ("immune" if lin in IMMUNE_LINEAGES else "other")
            rows.append({
                "gene": g, "lineage": lin, "compartment": comp, "n_cells": int(m.sum()),
                "mean_lognorm": float(v[m].mean()),
                "median_lognorm": float(np.median(v[m])),
                "pct_expressing": float((v[m] > 0).mean() * 100),
            })
    lineage_summary = pd.DataFrame(rows).sort_values(["gene", "mean_lognorm"], ascending=[True, False])
    lineage_summary.to_csv(os.path.join(out_dir, "celltype_summary.csv"), index=False)

    # ---------------- cell-level tests ----------------
    test_rows = []
    mal_mask = compartment == "malignant"
    imm_mask = compartment == "immune"
    for g, v in expr.items():
        comparisons = [("all_immune", imm_mask)]
        comparisons += [
            (lin, lineage == lin)
            for lin in sorted(set(lineage[imm_mask]))
            if (lineage == lin).sum() >= MIN_CELLS_PER_LINEAGE
        ]
        for name, mask in comparisons:
            a, b = v[mal_mask], v[mask]
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            eps = 1e-9
            test_rows.append({
                "gene": g, "comparison": f"malignant_vs_{name}",
                "n_malignant": len(a), "n_other": len(b),
                "mean_malignant": float(a.mean()), "mean_other": float(b.mean()),
                "pct_pos_malignant": float((a > 0).mean() * 100),
                "pct_pos_other": float((b > 0).mean() * 100),
                "log2FC_of_means": float(np.log2((a.mean() + eps) / (b.mean() + eps))),
                "auroc": auroc_from_u(u, len(a), len(b)),
                "mannwhitney_U": float(u), "p_value": float(p),
            })
    tests = pd.DataFrame(test_rows)
    # BH-FDR within dataset
    m = len(tests)
    order = np.argsort(tests["p_value"].to_numpy())
    ranked = tests["p_value"].to_numpy()[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty(m)
    fdr[order] = np.minimum(ranked, 1.0)
    tests["fdr_bh"] = fdr
    tests.to_csv(os.path.join(out_dir, "stats_malignant_vs_immune.csv"), index=False)

    # ---------------- per-patient pseudobulk ----------------
    pseudo_rows, paired_stats = [], {}
    if patient_col is not None:
        patients = meta[patient_col].astype(str).to_numpy()
        for g, v in expr.items():
            pairs = []
            for pt in sorted(set(patients)):
                pm, pi = (patients == pt) & mal_mask, (patients == pt) & imm_mask
                if pm.sum() >= MIN_CELLS_PER_PATIENT and pi.sum() >= MIN_CELLS_PER_PATIENT:
                    mm, mi = float(v[pm].mean()), float(v[pi].mean())
                    pairs.append((pt, int(pm.sum()), int(pi.sum()), mm, mi))
                    pseudo_rows.append({
                        "gene": g, "patient": pt,
                        "n_malignant": int(pm.sum()), "n_immune": int(pi.sum()),
                        "mean_malignant": mm, "mean_immune": mi,
                    })
            if len(pairs) >= 5:
                a = np.array([p[3] for p in pairs])
                b = np.array([p[4] for p in pairs])
                try:
                    w, pw = stats.wilcoxon(a, b, alternative="two-sided")
                except ValueError:
                    w, pw = np.nan, np.nan
                paired_stats[g] = {
                    "n_patients": len(pairs),
                    "median_mean_malignant": float(np.median(a)),
                    "median_mean_immune": float(np.median(b)),
                    "n_patients_malignant_higher": int((a > b).sum()),
                    "wilcoxon_W": float(w), "p_value": float(pw),
                }
    if pseudo_rows:
        pd.DataFrame(pseudo_rows).to_csv(os.path.join(out_dir, "per_patient_pseudobulk.csv"), index=False)

    # ---------------- epithelial restriction ----------------
    # Three-way contrast: malignant / residual non-malignant epithelium / immune.
    # Also: among TACSTD2+ (or CLDN4+) cells, what fraction sit in each
    # compartment? That is the "restriction" of the positive population.
    epi_mask = np.isin(lineage, list(NONMALIGNANT_EPITHELIAL))
    # if we already remapped Epithelial -> malignant_name, epi_mask is empty
    # unless Alveolar/Basal/etc. remain.
    n_epi = int(epi_mask.sum())
    epi_restriction = {"n_nonmalignant_epithelial": n_epi,
                       "nonmalignant_epithelial_lineages": sorted(set(lineage[epi_mask]))}
    if n_epi >= MIN_CELLS_PER_LINEAGE:
        for g, v in expr.items():
            a, e, b = v[mal_mask], v[epi_mask], v[imm_mask]
            u_me, p_me = stats.mannwhitneyu(a, e, alternative="two-sided")
            u_ei, p_ei = stats.mannwhitneyu(e, b, alternative="two-sided")
            eps = 1e-9
            epi_restriction[g] = {
                "n_malignant": int(len(a)),
                "n_nonmalignant_epithelial": int(len(e)),
                "n_immune": int(len(b)),
                "mean_malignant": float(a.mean()),
                "mean_nonmalignant_epithelial": float(e.mean()),
                "mean_immune": float(b.mean()),
                "pct_pos_malignant": float((a > 0).mean() * 100),
                "pct_pos_nonmalignant_epithelial": float((e > 0).mean() * 100),
                "pct_pos_immune": float((b > 0).mean() * 100),
                "log2FC_malignant_vs_nonmalig_epi": float(np.log2((a.mean() + eps) / (e.mean() + eps))),
                "log2FC_nonmalig_epi_vs_immune": float(np.log2((e.mean() + eps) / (b.mean() + eps))),
                "auroc_malignant_vs_nonmalig_epi": auroc_from_u(u_me, len(a), len(e)),
                "p_malignant_vs_nonmalig_epi": float(p_me),
                "auroc_nonmalig_epi_vs_immune": auroc_from_u(u_ei, len(e), len(b)),
                "p_nonmalig_epi_vs_immune": float(p_ei),
            }

    # composition of the gene-positive population (epithelial restriction)
    pos_composition = {}
    for g, v in expr.items():
        pos = v > 0
        n_pos = int(pos.sum())
        if n_pos == 0:
            continue
        pos_composition[g] = {
            "n_positive": n_pos,
            "pct_of_positive_in_malignant": float((pos & mal_mask).sum() / n_pos * 100),
            "pct_of_positive_in_nonmalignant_epithelial": float((pos & epi_mask).sum() / n_pos * 100),
            "pct_of_positive_in_immune": float((pos & imm_mask).sum() / n_pos * 100),
            "pct_of_positive_in_other": float((pos & ~(mal_mask | epi_mask | imm_mask)).sum() / n_pos * 100),
        }

    # ---------------- per-patient immune fraction ----------------
    # For each patient with enough malignant + immune cells, record the
    # immune infiltrate fraction and the malignant-compartment gene means,
    # then Spearman-correlate the two. This asks whether TACSTD2/CLDN4
    # tumor expression tracks immune content (it should not, if the signal
    # is epithelial-restricted rather than an infiltrate artifact).
    immune_fraction_rows, immune_fraction_corr = [], {}
    if patient_col is not None:
        patients = meta[patient_col].astype(str).to_numpy()
        for pt in sorted(set(patients)):
            pm = (patients == pt) & mal_mask
            pi = (patients == pt) & imm_mask
            p_all = patients == pt
            if pm.sum() < MIN_CELLS_PER_PATIENT or pi.sum() < MIN_CELLS_PER_PATIENT:
                continue
            row = {
                "patient": pt,
                "n_cells": int(p_all.sum()),
                "n_malignant": int(pm.sum()),
                "n_immune": int(pi.sum()),
                "immune_fraction": float(pi.sum() / p_all.sum()),
                "malignant_fraction": float(pm.sum() / p_all.sum()),
            }
            for g, v in expr.items():
                row[f"{g}_mean_malignant"] = float(v[pm].mean())
                row[f"{g}_pct_pos_malignant"] = float((v[pm] > 0).mean() * 100)
            immune_fraction_rows.append(row)
        if immune_fraction_rows:
            if_df = pd.DataFrame(immune_fraction_rows)
            if_df.to_csv(os.path.join(out_dir, "per_patient_immune_fraction.csv"), index=False)
            if len(if_df) >= 5:
                for g in expr:
                    rho, pr = stats.spearmanr(if_df["immune_fraction"], if_df[f"{g}_mean_malignant"])
                    immune_fraction_corr[g] = {
                        "n_patients": int(len(if_df)),
                        "spearman_rho_immune_frac_vs_malignant_mean": float(rho),
                        "p_value": float(pr),
                    }

    # ---------------- sensitivity: tumor-tissue cells only ----------------
    # GSE131907 also contains normal lung / lymph node / effusion / metastasis
    # samples; check that the malignant-vs-immune contrast holds when we
    # restrict to cells from tumor tissue.
    tumor_tissue_sensitivity = {}
    source_col = next((c for c in meta.columns if c.lower() in ("source", "tissue")), None)
    if source_col is not None:
        src = meta[source_col].astype(str).str.lower().to_numpy()
        in_tumor = np.char.find(src.astype(str), "tumor") >= 0
        if in_tumor.sum() > 0 and (mal_mask & in_tumor).sum() >= 50 and (imm_mask & in_tumor).sum() >= 50:
            for g, v in expr.items():
                a, b = v[mal_mask & in_tumor], v[imm_mask & in_tumor]
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                eps = 1e-9
                tumor_tissue_sensitivity[g] = {
                    "source_column": source_col,
                    "n_malignant": int(len(a)), "n_immune": int(len(b)),
                    "mean_malignant": float(a.mean()), "mean_immune": float(b.mean()),
                    "pct_pos_malignant": float((a > 0).mean() * 100),
                    "pct_pos_immune": float((b > 0).mean() * 100),
                    "log2FC_of_means": float(np.log2((a.mean() + eps) / (b.mean() + eps))),
                    "auroc": auroc_from_u(u, len(a), len(b)),
                    "p_value": float(p),
                }

    # ---------------- summary json ----------------
    summary = {
        "dataset": ds,
        "malignant_definition": malignant_definition,
        "n_cells_total": int(len(lineage)),
        "n_malignant": n_mal,
        "n_immune": n_imm,
        "immune_lineages_present": sorted(set(lineage[imm_mask])),
        "genes_found": sorted(expr.keys()),
        "cell_level_malignant_vs_all_immune": {
            r["gene"]: {k: r[k] for k in (
                "mean_malignant", "mean_other", "pct_pos_malignant", "pct_pos_other",
                "log2FC_of_means", "auroc", "p_value")}
            for _, r in tests[tests.comparison == "malignant_vs_all_immune"].iterrows()
        },
        "per_patient_paired_wilcoxon": paired_stats,
        "tumor_tissue_only_sensitivity": tumor_tissue_sensitivity,
        "epithelial_restriction": epi_restriction,
        "positive_population_composition": pos_composition,
        "immune_fraction_spearman": immune_fraction_corr,
        "tumor_vs_normal_epithelial": tumor_vs_normal_epithelial,
    }
    with open(os.path.join(out_dir, "stats_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---------------- figures ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_lins = [l for l in lineage_summary[lineage_summary.gene == GENES[0]]
                 .sort_values("compartment").lineage.unique()]
    # order: malignant first, then immune by TACSTD2 mean desc, then other
    def _key(l):
        comp = 0 if l == malignant_name else (1 if l in IMMUNE_LINEAGES else 2)
        mu = lineage_summary[(lineage_summary.gene == "TACSTD2") & (lineage_summary.lineage == l)]["mean_lognorm"]
        return (comp, -(float(mu.iloc[0]) if len(mu) else 0.0))
    plot_lins = sorted(plot_lins, key=_key)

    fig, axes = plt.subplots(len(GENES), 1, figsize=(max(8, 0.6 * len(plot_lins)), 3.0 * len(GENES)), sharex=True)
    rng = np.random.default_rng(0)
    for ax, g in zip(np.atleast_1d(axes), GENES):
        if g not in expr:
            ax.set_visible(False)
            continue
        v = expr[g]
        data = [v[lineage == l] for l in plot_lins]
        parts = ax.violinplot([d if len(d) else np.zeros(1) for d in data],
                              showmeans=False, showextrema=False, widths=0.85)
        for i, (pc, l) in enumerate(zip(parts["bodies"], plot_lins)):
            pc.set_facecolor("#c0392b" if l == malignant_name else ("#2980b9" if l in IMMUNE_LINEAGES else "#95a5a6"))
            pc.set_alpha(0.7)
        means = [float(d.mean()) if len(d) else 0.0 for d in data]
        ax.scatter(np.arange(1, len(plot_lins) + 1), means, s=14, c="k", zorder=3, label="mean")
        ax.set_ylabel(f"{g}\nlog(TPM/10+1)")
        ax.legend(loc="upper right", fontsize=8, frameon=False)
    np.atleast_1d(axes)[-1].set_xticks(np.arange(1, len(plot_lins) + 1))
    np.atleast_1d(axes)[-1].set_xticklabels(plot_lins, rotation=60, ha="right")
    fig.suptitle(f"{ds}: expression by major lineage (red=malignant, blue=immune)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "violin_by_lineage.png"), dpi=160)
    plt.close(fig)

    if pseudo_rows:
        pb = pd.DataFrame(pseudo_rows)
        fig, axes = plt.subplots(1, len(TARGET_GENES), figsize=(4.2 * len(TARGET_GENES), 4), sharey=False)
        for ax, g in zip(np.atleast_1d(axes), TARGET_GENES):
            sub = pb[pb.gene == g]
            if sub.empty:
                ax.set_visible(False)
                continue
            for _, r in sub.iterrows():
                ax.plot([0, 1], [r.mean_malignant, r.mean_immune], "-o", color="#7f8c8d",
                        alpha=0.6, ms=4, lw=1)
            st = paired_stats.get(g, {})
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["malignant", "immune"])
            ax.set_ylabel("patient mean log(TPM/10+1)")
            title = f"{g} (n={st.get('n_patients', len(sub))} patients"
            if "p_value" in st and np.isfinite(st["p_value"]):
                title += f", Wilcoxon p={st['p_value']:.2e}"
            ax.set_title(title + ")", fontsize=10)
        fig.suptitle(f"{ds}: per-patient pseudobulk, malignant vs immune")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "per_patient_pseudobulk.png"), dpi=160)
        plt.close(fig)

    # dot-plot style summary: pct expressing vs mean
    fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(plot_lins))))
    ycoords = {l: i for i, l in enumerate(plot_lins[::-1])}
    colors = {"TACSTD2": "#c0392b", "CLDN4": "#e67e22", "EPCAM": "#27ae60", "PTPRC": "#2980b9"}
    for g in GENES:
        sub = lineage_summary[lineage_summary.gene == g]
        for _, r in sub.iterrows():
            if r.lineage not in ycoords:
                continue
            ax.scatter(r.mean_lognorm, ycoords[r.lineage] + {"TACSTD2": -0.25, "CLDN4": -0.08, "EPCAM": 0.08, "PTPRC": 0.25}[g],
                       s=8 + 1.6 * r.pct_expressing, color=colors[g], alpha=0.75,
                       label=g if r.lineage == plot_lins[0] else None)
    ax.set_yticks(range(len(plot_lins)))
    ax.set_yticklabels(plot_lins[::-1])
    ax.set_xlabel("mean log(TPM/10+1)  (dot size = % expressing)")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_title(f"{ds}: mean expression / % expressing by lineage")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "dotplot_summary.png"), dpi=160)
    plt.close(fig)

    # three-way epithelial-restriction violins
    if n_epi >= MIN_CELLS_PER_LINEAGE:
        fig, axes = plt.subplots(1, len(TARGET_GENES), figsize=(8, 4), sharey=False)
        labels = ["malignant", "non-malig.\nepithelial", "immune"]
        colors = ["#c0392b", "#27ae60", "#2980b9"]
        for ax, g in zip(np.atleast_1d(axes), TARGET_GENES):
            if g not in expr:
                ax.set_visible(False)
                continue
            v = expr[g]
            data = [v[mal_mask], v[epi_mask], v[imm_mask]]
            parts = ax.violinplot(data, showmeans=False, showextrema=False, widths=0.8)
            for pc, c in zip(parts["bodies"], colors):
                pc.set_facecolor(c)
                pc.set_alpha(0.7)
            ax.scatter([1, 2, 3], [d.mean() for d in data], s=16, c="k", zorder=3)
            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels(labels)
            ax.set_ylabel(f"{g} log(TPM/10+1)")
            st = epi_restriction.get(g, {})
            ax.set_title(
                f"{g}\nmalig vs epi p={st.get('p_malignant_vs_nonmalig_epi', float('nan')):.1e}\n"
                f"epi vs imm p={st.get('p_nonmalig_epi_vs_immune', float('nan')):.1e}",
                fontsize=9,
            )
        fig.suptitle(f"{ds}: epithelial restriction (malignant / leftover epithelium / immune)")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "epithelial_restriction.png"), dpi=160)
        plt.close(fig)

    if immune_fraction_rows and len(immune_fraction_rows) >= 5:
        if_df = pd.DataFrame(immune_fraction_rows)
        fig, axes = plt.subplots(1, len(TARGET_GENES), figsize=(8, 4), sharey=False)
        for ax, g in zip(np.atleast_1d(axes), TARGET_GENES):
            ax.scatter(if_df["immune_fraction"], if_df[f"{g}_mean_malignant"],
                       s=28, c="#c0392b", alpha=0.8)
            st = immune_fraction_corr.get(g, {})
            ax.set_xlabel("patient immune fraction")
            ax.set_ylabel(f"malignant {g} mean")
            ax.set_title(
                f"{g}  Spearman ρ={st.get('spearman_rho_immune_frac_vs_malignant_mean', float('nan')):.2f}"
                f"  p={st.get('p_value', float('nan')):.2e}",
                fontsize=9,
            )
        fig.suptitle(f"{ds}: malignant gene mean vs immune infiltrate fraction")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "immune_fraction_scatter.png"), dpi=160)
        plt.close(fig)

    print(f"done -> {out_dir}")


if __name__ == "__main__":
    main()
