"""Claim A10 / P1-P3: human bulk lung co-expression.

Datasets: GTEx v8 lung (normal, n~655) and TCGA-LUAD (tumour + adjacent normal),
both from recount3 (uniform reprocessing, so the two are directly comparable).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

N_BG_PAIRS = 200_000
N_MATCHED_SETS = 2000


def load_dataset(name: str, counts_file: str, annot: pd.DataFrame):
    counts = C.read_recount3_counts(os.path.join(C.RAW, counts_file))
    logcpm = C.to_log_cpm(counts)
    mat = C.collapse_to_symbols(logcpm, annot)
    uni = C.expressed_genes(mat)
    print(f"[{name}] {mat.shape[0]} symbols x {mat.shape[1]} samples; "
          f"{len(uni)} expressed genes", flush=True)
    return mat, uni


def detection_report(mat: pd.DataFrame, name: str) -> pd.DataFrame:
    genes = C.CLAIM_GENES_H + list(C.MARKERS_H)
    rows = []
    for g in genes:
        if g in mat.index:
            v = mat.loc[g]
            rows.append({"dataset": name, "gene": g, "detected": True,
                         "mean_log2cpm": float(v.mean()),
                         "sd_log2cpm": float(v.std()),
                         "frac_samples_gt1_log2cpm": float((v > 1).mean())})
        else:
            rows.append({"dataset": name, "gene": g, "detected": False,
                         "mean_log2cpm": np.nan, "sd_log2cpm": np.nan,
                         "frac_samples_gt1_log2cpm": np.nan})
    return pd.DataFrame(rows)


def pairwise_with_background(mat, uni, name, genes):
    bg = C.background_pair_distribution(mat, uni, n_pairs=N_BG_PAIRS)
    df = C.spearman_pairs(mat, genes)
    df.insert(0, "dataset", name)
    df["bg_percentile"] = [C.signed_percentile(r, bg) for r in df["spearman_rho"]]
    ci = [C.fisher_ci(r, n) for r, n in zip(df["spearman_rho"], df["n"])]
    df["ci_low"], df["ci_high"] = [c[0] for c in ci], [c[1] for c in ci]
    return df, bg


def module_coherence(mat, uni, name, module_genes, bg):
    """P1: is the module more internally coherent than expression-matched random sets?"""
    bins = C.expression_bins(mat, uni)
    obs = C.mean_pairwise_rho(mat, module_genes)
    sets = C.matched_random_sets(bins, module_genes, N_MATCHED_SETS)
    null = np.array([C.mean_pairwise_rho(mat, s) for s in sets])
    p = float((null >= obs).mean())
    return {
        "dataset": name,
        "module": ",".join(module_genes),
        "mean_pairwise_rho": obs,
        "matched_null_mean": float(np.nanmean(null)),
        "matched_null_p95": float(np.nanpercentile(null, 95)),
        "empirical_p_vs_matched_null": p,
        "percentile_vs_random_pairs": C.signed_percentile(obs, bg),
    }


def composition_controlled(mat, name):
    """P6-in-bulk: do claim associations survive adjusting for cell-composition proxies?

    Covariates: EPCAM (epithelial content) alone, then EPCAM plus airway/alveolar
    markers, then all lineage markers. If the association is only 'this biopsy has
    more airway epithelium', it should collapse.
    """
    covsets = {
        "none": [],
        "EPCAM": ["EPCAM"],
        "EPCAM+SFTPC+SCGB1A1": ["EPCAM", "SFTPC", "SCGB1A1"],
        "all_lineage_markers": list(C.MARKERS_H),
    }
    pairs = [(a, b) for a in C.MODULE_H for b in C.TARGET_H]
    pairs += [(a, C.NKX_H) for a in C.MODULE_H + C.TARGET_H]
    rows = []
    for cname, cgenes in covsets.items():
        cg = [g for g in cgenes if g in mat.index]
        covars = mat.loc[cg].T.to_numpy(dtype=float) if cg else None
        for a, b in pairs:
            if a not in mat.index or b not in mat.index:
                continue
            x = mat.loc[a].to_numpy(dtype=float)
            y = mat.loc[b].to_numpy(dtype=float)
            if covars is None:
                rho, p = stats.spearmanr(x, y)
            else:
                rho, p = C.partial_spearman(x, y, covars)
            rows.append({"dataset": name, "covariates": cname, "gene_a": a,
                         "gene_b": b, "n": len(x), "partial_rho": float(rho),
                         "p_value": float(p)})
    df = pd.DataFrame(rows)
    df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df


def module_vs_target(mat, name):
    """Module score against each target and against NKX2-1."""
    score = C.module_score(mat, C.MODULE_H)
    rows = []
    for g in C.TARGET_H + [C.NKX_H] + list(C.MARKERS_H):
        if g not in mat.index:
            continue
        rho, p = stats.spearmanr(score.to_numpy(), mat.loc[g].to_numpy(dtype=float))
        lo, hi = C.fisher_ci(rho, len(score))
        rows.append({"dataset": name, "x": "MODULE_score(ELF3,GRHL1,KLF4,TFAP2A)",
                     "y": g, "n": len(score), "spearman_rho": float(rho),
                     "ci_low": lo, "ci_high": hi, "p_value": float(p)})
    df = pd.DataFrame(rows)
    df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df, score


def _align_tcga_meta(mat, meta):
    """recount3 TCGA matrices are keyed by GDC file UUID, not rail_id."""
    for col in ("gdc_file_id", "external_id", "tcga_barcode"):
        if col not in meta.columns:
            continue
        m = meta.copy()
        m.index = m[col].astype(str)
        overlap = mat.columns.intersection(m.index)
        if len(overlap) >= 0.8 * mat.shape[1]:
            return m.loc[overlap]
    raise RuntimeError("could not align TCGA metadata to matrix columns")


def tcga_nkx_strata(mat, meta_path):
    """Natural perturbation: NKX2-1-low vs NKX2-1-high tumours.

    The claim predicts the module and its targets are UP where NKX2-1 is DOWN.
    """
    meta = C.read_recount3_metadata(meta_path)
    m = _align_tcga_meta(mat, meta)
    samp_type = None
    for cand in ["cgc_sample_sample_type", "gdc_cases.samples.sample_type",
                 "tcga.cgc_sample_sample_type"]:
        if cand in m.columns:
            samp_type = cand
            break
    types = m[samp_type].astype(str) if samp_type else pd.Series("unknown", index=m.index)
    is_tumour = types.str.contains("Tumor", case=False, na=False)
    tumours = types.index[is_tumour].tolist()
    normals = types.index[~is_tumour].tolist()
    print(f"[TCGA-LUAD] aligned {len(m)} samples; {len(tumours)} tumour, "
          f"{len(normals)} non-tumour", flush=True)

    nkx = mat.loc[C.NKX_H, tumours]
    lo_cut, hi_cut = nkx.quantile(0.25), nkx.quantile(0.75)
    lo = nkx.index[nkx <= lo_cut]
    hi = nkx.index[nkx >= hi_cut]
    rows = []
    for g in C.MODULE_H + C.TARGET_H + list(C.MARKERS_H):
        if g not in mat.index:
            continue
        a = mat.loc[g, lo].to_numpy(dtype=float)
        b = mat.loc[g, hi].to_numpy(dtype=float)
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"dataset": "TCGA-LUAD", "comparison": "NKX2-1-low vs NKX2-1-high tumours",
                     "gene": g, "n_low": len(a), "n_high": len(b),
                     "mean_log2cpm_low": float(a.mean()), "mean_log2cpm_high": float(b.mean()),
                     "log2fc_low_minus_high": float(a.mean() - b.mean()),
                     "cliffs_delta": C.cliffs_delta(a, b), "p_value": float(p)})
    df = pd.DataFrame(rows)
    df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df, tumours, normals


def main():
    annot = C.read_recount3_annotation(os.path.join(C.RAW, "human.gene_sums.G026.gtf.gz"))
    all_pairs, all_coh, all_comp, all_mod, all_det = [], [], [], [], []

    specs = [("GTEx-lung", "gtex.gene_sums.LUNG.G026.gz"),
             ("TCGA-LUAD", "tcga.gene_sums.LUAD.G026.gz")]
    mats = {}
    for name, f in specs:
        mat, uni = load_dataset(name, f, annot)
        mats[name] = (mat, uni)
        all_det.append(detection_report(mat, name))
        genes = C.CLAIM_GENES_H + list(C.MARKERS_H)
        df, bg = pairwise_with_background(mat, uni, name, genes)
        all_pairs.append(df)
        all_coh.append(module_coherence(mat, uni, name, C.MODULE_H, bg))
        all_comp.append(composition_controlled(mat, name))
        mv, score = module_vs_target(mat, name)
        all_mod.append(mv)
        score.to_frame("module_score").to_csv(
            os.path.join(C.INTERIM, f"module_score_{name}.csv"))
        np.save(os.path.join(C.INTERIM, f"bg_{name}.npy"), bg)

    # TCGA tumour strata + tumour vs normal
    mat_t, _ = mats["TCGA-LUAD"]
    strata, tumours, normals = tcga_nkx_strata(
        mat_t, os.path.join(C.RAW, "tcga.tcga.LUAD.MD.gz"))

    tn_rows = []
    for g in C.CLAIM_GENES_H + list(C.MARKERS_H):
        if g not in mat_t.index or not normals:
            continue
        a = mat_t.loc[g, tumours].to_numpy(dtype=float)
        b = mat_t.loc[g, normals].to_numpy(dtype=float)
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        tn_rows.append({"dataset": "TCGA-LUAD", "comparison": "tumour vs adjacent normal",
                        "gene": g, "n_tumour": len(a), "n_normal": len(b),
                        "log2fc_tumour_minus_normal": float(a.mean() - b.mean()),
                        "cliffs_delta": C.cliffs_delta(a, b), "p_value": float(p)})
    tn = pd.DataFrame(tn_rows)
    if len(tn):
        tn["fdr"] = C.bh_fdr(tn["p_value"].to_numpy())

    C.write_table(pd.concat(all_det, ignore_index=True), "human_bulk_gene_detection.csv")
    C.write_table(pd.concat(all_pairs, ignore_index=True), "human_bulk_pairwise_spearman.csv")
    C.write_table(pd.DataFrame(all_coh), "human_bulk_module_coherence.csv")
    C.write_table(pd.concat(all_comp, ignore_index=True), "human_bulk_partial_correlation.csv")
    C.write_table(pd.concat(all_mod, ignore_index=True), "human_bulk_module_score_assoc.csv")
    C.write_table(strata, "human_tcga_nkx2-1_low_vs_high.csv")
    C.write_table(tn, "human_tcga_tumour_vs_normal.csv")
    print("done", flush=True)


if __name__ == "__main__":
    main()
