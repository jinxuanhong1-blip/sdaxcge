#!/usr/bin/env python3
"""Score epithelial Cldn4 vs T/NK and exclusion in leftover mouse ICI scRNA."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    DATA,
    EXCLUSION_T,
    OUT,
    T_SCORE_GENES,
    classify,
    dump_json,
    gene_index,
    high_low_cut,
    load_panel,
    log1p_mat,
    mean_of,
    mtx_shape,
    mwu_or_nan,
    read_features,
    spearman,
    stream_mtx_panel,
    welch_or_nan,
    wilcoxon_signed,
)

SAMPLES = [
    dict(
        gse="GSE297632",
        gsm="GSM8995911",
        sample="LLC_Control",
        arm="control",
        genotype="LLC",
        icb="control",
        path="GSE297632/GSM8995911_Control_GEX_matrix.mtx.gz",
        features="GSE297632/GSM8995911_Control_GEX_features.tsv.gz",
        min_umi=200,
        note="Subcutaneous LLC; untreated control; n=1 library",
    ),
    dict(
        gse="GSE297632",
        gsm="GSM8995912",
        sample="LLC_aPD1_residual",
        arm="aPD1_residual",
        genotype="LLC",
        icb="icb",
        path="GSE297632/GSM8995912_PD_1_treatment_GEX_matrix.mtx.gz",
        features="GSE297632/GSM8995912_PD_1_treatment_GEX_features.tsv.gz",
        min_umi=200,
        note="Subcutaneous LLC; residual/tolerant cells 7d after anti-PD-1; n=1 library",
    ),
    dict(
        gse="GSE285606",
        gsm="GSM8705591",
        sample="WT_IgG_T1",
        arm="WT_IgG",
        genotype="344SQ_WT",
        icb="igg",
        path="GSE285606/GSM8705591_WT-IgG-T1.matrix.mtx.gz",
        features="GSE285606/GSM8705591_WT-IgG-T1.features.tsv.gz",
        min_umi=200,
        note="344SQ parental; IgG week 4; not ICB vs control",
    ),
    dict(
        gse="GSE285606",
        gsm="GSM8705592",
        sample="WT_IgG_T2",
        arm="WT_IgG",
        genotype="344SQ_WT",
        icb="igg",
        path="GSE285606/GSM8705592_WT-IgG-T2.matrix.mtx.gz",
        features="GSE285606/GSM8705592_WT-IgG-T2.features.tsv.gz",
        min_umi=200,
        note="344SQ parental; IgG week 6",
    ),
    dict(
        gse="GSE285606",
        gsm="GSM8705589",
        sample="PD1R1_IgG_T1",
        arm="PD1R1_IgG",
        genotype="344SQ_PD1R1",
        icb="igg",
        path="GSE285606/GSM8705589_140P-IgG-T1.matrix.mtx.gz",
        features="GSE285606/GSM8705589_140P-IgG-T1.features.tsv.gz",
        min_umi=200,
        note="344SQ PD1R1 (140P) acquired aPD-1-resistant; IgG week 4",
    ),
    dict(
        gse="GSE285606",
        gsm="GSM8705590",
        sample="PD1R1_IgG_T2",
        arm="PD1R1_IgG",
        genotype="344SQ_PD1R1",
        icb="igg",
        path="GSE285606/GSM8705590_140P-IgG-T2.matrix.mtx.gz",
        features="GSE285606/GSM8705590_140P-IgG-T2.features.tsv.gz",
        min_umi=200,
        note="344SQ PD1R1 (140P); IgG week 6",
    ),
]


def score_one(cfg, wanted):
    path = DATA / cfg["path"]
    feat = DATA / cfg["features"]
    rec = {k: cfg.get(k) for k in ["gse", "gsm", "sample", "arm", "genotype", "icb", "note"]}
    if not path.exists() or not feat.exists():
        rec["status"] = "missing_file"
        rec["detail"] = f"{path.exists()=} {feat.exists()=}"
        return rec, None
    genes = read_features(feat)
    ng, nc, nz = mtx_shape(path)
    rec["mtx_shape"] = f"{ng}x{nc}"
    rec["n_genes_features"] = len(genes)
    if ng != len(genes):
        rec["status"] = "features_mismatch"
        rec["detail"] = f"mtx {ng} genes vs features {len(genes)}"
        return rec, None
    idx = gene_index(genes, wanted)
    rec["genes_found"] = ",".join(sorted(idx))
    rec["Cldn4_in_matrix"] = "Cldn4" in idx
    rec["Tacstd2_in_matrix"] = "Tacstd2" in idx
    mat, umi, n_genes, n_use, n_drop, _keep = stream_mtx_panel(path, idx, nc, cfg.get("min_umi"))
    rec["n_cells_raw"] = nc
    rec["n_cells_used"] = n_use
    rec["n_cells_dropped"] = n_drop
    rec["median_umi"] = float(np.median(umi)) if n_use else np.nan
    labels = classify(mat)
    rec["n_epithelial"] = int((labels == "epithelial").sum())
    rec["n_tnk"] = int((labels == "tnk").sum())
    rec["n_other"] = int((labels == "other").sum())
    rec["frac_epithelial"] = rec["n_epithelial"] / n_use if n_use else np.nan
    rec["frac_tnk"] = rec["n_tnk"] / n_use if n_use else np.nan
    logmat = log1p_mat(mat)
    tscore = mean_of(logmat, T_SCORE_GENES)
    excl_t = mean_of(logmat, EXCLUSION_T)
    cldn = logmat.get("Cldn4")
    tac = logmat.get("Tacstd2")
    rec["status"] = "ok"

    # compartment means
    for comp in ["epithelial", "tnk", "other"]:
        mask = labels == comp
        k = int(mask.sum())
        rec[f"{comp}_Cldn4_mean"] = float(cldn[mask].mean()) if k and cldn is not None else np.nan
        rec[f"{comp}_Cldn4_pct"] = float((mat["Cldn4"][mask] > 0).mean() * 100) if k and "Cldn4" in mat else np.nan
        rec[f"{comp}_Tacstd2_mean"] = float(tac[mask].mean()) if k and tac is not None else np.nan
        rec[f"{comp}_Tscore_mean"] = float(tscore[mask].mean()) if k and tscore is not None else np.nan

    # same-library epithelial Cldn4 vs T/NK Cldn4 / T fraction
    rec["epi_minus_tnk_Cldn4"] = rec["epithelial_Cldn4_mean"] - rec["tnk_Cldn4_mean"]
    rec["exclusion_Cldn4_minus_Tfrac"] = rec["epithelial_Cldn4_mean"] - rec["frac_tnk"]
    rec["exclusion_Cldn4_minus_Tscore"] = rec["epithelial_Cldn4_mean"] - rec["tnk_Tscore_mean"]

    # cell-level exploratory (not a sample-level claim)
    if cldn is not None and tscore is not None:
        rec["cell_rho_Cldn4_Tscore"] = spearman(cldn, tscore)
        epi = labels == "epithelial"
        if epi.sum() >= 20:
            rec["epi_cell_rho_Cldn4_Tscore"] = spearman(cldn[epi], tscore[epi])
        for how in ("median", "quartile"):
            mask, meta = high_low_cut(cldn, tscore, how=how)
            rec[f"hiCldn4_loT_{how}"] = meta
            if epi.sum() >= 20:
                emask, emeta = high_low_cut(cldn[epi], tscore[epi], how=how)
                rec[f"epi_hiCldn4_loT_{how}"] = emeta

    cells = pd.DataFrame(
        {
            "gse": cfg["gse"],
            "sample": cfg["sample"],
            "arm": cfg["arm"],
            "label": labels,
            "umi": umi,
            "Cldn4": cldn if cldn is not None else np.nan,
            "Tacstd2": tac if tac is not None else np.nan,
            "Tscore": tscore if tscore is not None else np.nan,
            "excl_T": excl_t if excl_t is not None else np.nan,
        }
    )
    return rec, cells


def add_contrast(rows, contrasts, gse, name, a_samples, b_samples, field, extra):
    sub = [r for r in rows if r.get("gse") == gse and r.get("status") == "ok"]
    a = [r[field] for r in sub if r["sample"] in a_samples]
    b = [r[field] for r in sub if r["sample"] in b_samples]
    if not a or not b:
        return
    na, nb = len(a), len(b)
    mean_a = float(np.nanmean(a))
    mean_b = float(np.nanmean(b))
    delta = mean_b - mean_a
    _, p_t = welch_or_nan(b, a)
    _, p_u = mwu_or_nan(b, a)
    contrasts.append(
        {
            "gse": gse,
            "contrast": name,
            "field": field,
            "n_a": na,
            "n_b": nb,
            "mean_a": mean_a,
            "mean_b": mean_b,
            "delta": delta,
            "direction": "up" if delta > 0 else "down" if delta < 0 else "tie",
            "welch_p": p_t,
            "mwu_p": p_u,
            "note": extra if (na >= 2 and nb >= 2) else f"{extra}; n={na} vs {nb} — p not a sample-level claim" if min(na, nb) < 2 else extra,
        }
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wanted, _ = load_panel()
    sample_rows = []
    cell_frames = []
    for cfg in SAMPLES:
        print("SCORE", cfg["gse"], cfg["sample"], flush=True)
        rec, cells = score_one(cfg, wanted)
        sample_rows.append(rec)
        if cells is not None:
            cell_frames.append(cells)
        print(
            " ",
            rec.get("status"),
            "epi",
            rec.get("n_epithelial"),
            "tnk",
            rec.get("n_tnk"),
            "Cldn4",
            rec.get("Cldn4_in_matrix"),
            flush=True,
        )

    samp = pd.DataFrame(sample_rows)
    samp.to_csv(OUT / "scrna_sample_inventory.tsv", sep="\t", index=False)
    if cell_frames:
        cells = pd.concat(cell_frames, ignore_index=True)
        cells.to_csv(OUT / "scrna_cell_scores.tsv.gz", sep="\t", index=False)
    else:
        cells = pd.DataFrame()

    contrasts = []
    add_contrast(
        sample_rows,
        contrasts,
        "GSE297632",
        "aPD1_residual_vs_control",
        ["LLC_Control"],
        ["LLC_aPD1_residual"],
        "epithelial_Cldn4_mean",
        "subcutaneous LLC; n=1 vs 1",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE297632",
        "aPD1_residual_vs_control",
        ["LLC_Control"],
        ["LLC_aPD1_residual"],
        "frac_tnk",
        "subcutaneous LLC; n=1 vs 1",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE297632",
        "aPD1_residual_vs_control",
        ["LLC_Control"],
        ["LLC_aPD1_residual"],
        "exclusion_Cldn4_minus_Tscore",
        "subcutaneous LLC; n=1 vs 1",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE285606",
        "PD1R1_vs_WT_IgG",
        ["WT_IgG_T1", "WT_IgG_T2"],
        ["PD1R1_IgG_T1", "PD1R1_IgG_T2"],
        "epithelial_Cldn4_mean",
        "both arms IgG; acquired-resistant line vs parental; not ICB vs control",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE285606",
        "PD1R1_vs_WT_IgG",
        ["WT_IgG_T1", "WT_IgG_T2"],
        ["PD1R1_IgG_T1", "PD1R1_IgG_T2"],
        "frac_tnk",
        "both arms IgG; not ICB vs control",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE285606",
        "PD1R1_vs_WT_IgG",
        ["WT_IgG_T1", "WT_IgG_T2"],
        ["PD1R1_IgG_T1", "PD1R1_IgG_T2"],
        "tnk_Tscore_mean",
        "both arms IgG; not ICB vs control",
    )
    add_contrast(
        sample_rows,
        contrasts,
        "GSE285606",
        "PD1R1_vs_WT_IgG",
        ["WT_IgG_T1", "WT_IgG_T2"],
        ["PD1R1_IgG_T1", "PD1R1_IgG_T2"],
        "exclusion_Cldn4_minus_Tscore",
        "both arms IgG; not ICB vs control",
    )

    # paired epi vs T/NK Cldn4 within library
    ok = [r for r in sample_rows if r.get("status") == "ok" and r.get("n_epithelial", 0) >= 20 and r.get("n_tnk", 0) >= 20]
    if ok:
        epi = np.array([r["epithelial_Cldn4_mean"] for r in ok], float)
        tnk = np.array([r["tnk_Cldn4_mean"] for r in ok], float)
        w = wilcoxon_signed(epi, tnk)
        contrasts.append(
            {
                "gse": "POOL",
                "contrast": "epithelial_vs_tnk_Cldn4_paired",
                "field": "Cldn4_mean",
                "n_a": w["n"],
                "n_b": w["n"],
                "mean_a": float(np.nanmean(epi)),
                "mean_b": float(np.nanmean(tnk)),
                "delta": float(np.nanmean(epi - tnk)),
                "direction": f"{w['n_a_gt_b']}/{w['n']} epi>T/NK",
                "welch_p": np.nan,
                "mwu_p": w["p"],
                "note": "Wilcoxon signed-rank on libraries with ≥20 epi and ≥20 T/NK; not an ICB contrast",
            }
        )

    pd.DataFrame(contrasts).to_csv(OUT / "scrna_contrasts.tsv", sep="\t", index=False)
    dump_json(OUT / "scrna_sample_records.json", sample_rows)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
