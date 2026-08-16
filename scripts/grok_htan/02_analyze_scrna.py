#!/usr/bin/env python3
"""TACSTD2/CLDN4 vs immune in HTAN/HCA scRNA objects (<2GB)."""

from __future__ import annotations

import json
import traceback
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from lib_htan import (
    CLDN4_IDS,
    TACSTD2_IDS,
    bh_fdr,
    classify_compartment,
    ensure_dirs,
    extract_vector,
    is_cd8,
    is_tcell,
    mwu_safe,
    pick_celltype_col,
    pick_sample_col,
    resolve_gene_index,
    signature_scores,
    spearman_safe,
    write_json,
)

SCRNA_KEYS = [
    "hca_travaglini_10x",
    "hca_madissoon_lung",
    "htan_sclc_combined",
    "luad_histology_scrna",
]


def load_adata(path: Path):
    # backed read then subset genes of interest into memory
    a = ad.read_h5ad(path, backed="r")
    return a


def gene_pair(adata):
    t_idx, t_name = resolve_gene_index(adata.var, TACSTD2_IDS)
    c_idx, c_name = resolve_gene_index(adata.var, CLDN4_IDS)
    return t_idx, t_name, c_idx, c_name


def analyze_one(key: str, path: Path, out_tables: Path, out_fig: Path):
    adata = load_adata(path)
    ct_col = pick_celltype_col(adata.obs)
    samp_col = pick_sample_col(adata.obs)
    t_idx, t_name, c_idx, c_name = gene_pair(adata)
    info = {
        "key": key,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "celltype_col": ct_col,
        "sample_col": samp_col,
        "tacstd2_name": t_name,
        "cldn4_name": c_name,
        "obs_columns": list(map(str, adata.obs.columns)),
    }
    if t_idx is None or c_idx is None:
        info["error"] = "TACSTD2 or CLDN4 not found"
        write_json(out_tables / f"{key}_info.json", info)
        return info, None, None

    tac = extract_vector(adata, t_idx)
    cld = extract_vector(adata, c_idx)
    scores, coverage = signature_scores(adata)
    info["signature_coverage"] = coverage

    obs = adata.obs.copy()
    obs["TACSTD2"] = tac
    obs["CLDN4"] = cld
    for col in scores.columns:
        obs[col] = scores[col].to_numpy()
    if ct_col:
        obs["_compartment"] = obs[ct_col].map(classify_compartment)
        obs["_is_cd8"] = obs[ct_col].map(is_cd8)
        obs["_is_t"] = obs[ct_col].map(is_tcell)
    else:
        obs["_compartment"] = "unknown"
        obs["_is_cd8"] = False
        obs["_is_t"] = False

    # cell-type summary
    if ct_col:
        g = obs.groupby(ct_col, observed=True)
        ct_tab = g.agg(
            n=("TACSTD2", "size"),
            TACSTD2_mean=("TACSTD2", "mean"),
            TACSTD2_pct=( "TACSTD2", lambda s: float((s > 0).mean())),
            CLDN4_mean=("CLDN4", "mean"),
            CLDN4_pct=("CLDN4", lambda s: float((s > 0).mean())),
        ).reset_index()
        ct_tab["dataset"] = key
        ct_tab.to_csv(out_tables / f"{key}_celltype_expression.tsv", sep="\t", index=False)
    else:
        ct_tab = pd.DataFrame()

    # co-expression
    co_rows = []
    for label, mask in {
        "all": np.ones(len(obs), dtype=bool),
        "epithelial": (obs["_compartment"] == "epithelial").to_numpy(),
        "immune": (obs["_compartment"] == "immune").to_numpy(),
    }.items():
        r = spearman_safe(obs.loc[mask, "TACSTD2"], obs.loc[mask, "CLDN4"])
        r.update({"dataset": key, "subset": label, "pair": "TACSTD2_vs_CLDN4"})
        co_rows.append(r)

    # cell-level gene vs immune signatures (epithelial only)
    epi = obs["_compartment"] == "epithelial"
    sig_rows = []
    for gene in ("TACSTD2", "CLDN4"):
        for sig in scores.columns:
            r = spearman_safe(obs.loc[epi, gene], obs.loc[epi, sig])
            r.update({"dataset": key, "gene": gene, "signature": sig, "subset": "epithelial"})
            sig_rows.append(r)

    # sample-level: epithelial mean genes vs immune fractions
    sample_rows = []
    frac_rows = []
    if samp_col:
        for samp, sub in obs.groupby(samp_col, observed=True):
            n = len(sub)
            if n < 30:
                continue
            epi_s = sub["_compartment"] == "epithelial"
            rec = {
                "dataset": key,
                "sample": str(samp),
                "n": n,
                "n_epithelial": int(epi_s.sum()),
                "frac_epithelial": float(epi_s.mean()),
                "frac_immune": float((sub["_compartment"] == "immune").mean()),
                "frac_cd8": float(sub["_is_cd8"].mean()),
                "frac_t": float(sub["_is_t"].mean()),
                "TACSTD2_all": float(sub["TACSTD2"].mean()),
                "CLDN4_all": float(sub["CLDN4"].mean()),
                "TACSTD2_epi": float(sub.loc[epi_s, "TACSTD2"].mean()) if epi_s.any() else np.nan,
                "CLDN4_epi": float(sub.loc[epi_s, "CLDN4"].mean()) if epi_s.any() else np.nan,
            }
            for sig in scores.columns:
                rec[f"{sig}_all"] = float(sub[sig].mean())
            sample_rows.append(rec)
        samp_df = pd.DataFrame(sample_rows)
        if len(samp_df) >= 4:
            samp_df.to_csv(out_tables / f"{key}_sample_summary.tsv", sep="\t", index=False)
            for gene in ("TACSTD2_epi", "CLDN4_epi", "TACSTD2_all", "CLDN4_all"):
                for feat in ("frac_immune", "frac_cd8", "frac_t", "CYT_all", "CD8_Tcell_all",
                             "IFNG_6gene_all", "TLS_12chemokine_all", "TGFB_exclusion_all"):
                    if feat not in samp_df.columns:
                        continue
                    r = spearman_safe(samp_df[gene], samp_df[feat], min_n=4)
                    r.update({"dataset": key, "x": gene, "y": feat, "level": "sample"})
                    frac_rows.append(r)

    assoc = pd.DataFrame(co_rows + sig_rows + frac_rows)
    if len(assoc):
        assoc["q"] = bh_fdr(assoc["p"].to_numpy())
        assoc.to_csv(out_tables / f"{key}_associations.tsv", sep="\t", index=False)

    # figures
    sns.set_theme(style="whitegrid", font_scale=0.9)
    if ct_col and len(ct_tab):
        top = ct_tab.sort_values("n", ascending=False).head(18)
        fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
        for ax, gene in zip(axes, ("TACSTD2_mean", "CLDN4_mean")):
            sns.barplot(data=top, y=ct_col, x=gene, ax=ax, color="#3b6d9a")
            ax.set_xlabel(gene)
            ax.set_ylabel("")
        fig.suptitle(f"{key}: mean expression by cell type (top 18 by n)")
        fig.tight_layout()
        fig.savefig(out_fig / f"{key}_celltype_means.png", dpi=140)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    if epi.any() and epi.sum() > 2000:
        rng = np.random.default_rng(0)
        take = rng.choice(np.flatnonzero(epi.to_numpy()), size=2000, replace=False)
        x = obs.iloc[take]["TACSTD2"]
        y = obs.iloc[take]["CLDN4"]
    else:
        x = obs.loc[epi, "TACSTD2"] if epi.any() else obs["TACSTD2"]
        y = obs.loc[epi, "CLDN4"] if epi.any() else obs["CLDN4"]
    ax.scatter(x, y, s=6, alpha=0.25, c="#2c5f8a")
    r = spearman_safe(obs.loc[epi, "TACSTD2"] if epi.any() else obs["TACSTD2"],
                      obs.loc[epi, "CLDN4"] if epi.any() else obs["CLDN4"])
    ax.set_xlabel("TACSTD2")
    ax.set_ylabel("CLDN4")
    ax.set_title(f"{key} epithelial co-expression\nSpearman ρ={r['rho']:.3f} p={r['p']:.2e} n={r['n']}")
    fig.tight_layout()
    fig.savefig(out_fig / f"{key}_coexpression.png", dpi=140)
    plt.close(fig)

    write_json(out_tables / f"{key}_info.json", info)
    # close backed file
    if hasattr(adata, "file") and adata.file is not None:
        adata.file.close()
    return info, assoc, ct_tab


def main():
    paths = ensure_dirs()
    tables, figs = paths["tables"], paths["figures"]
    all_assoc = []
    all_ct = []
    runlog = []
    for key in SCRNA_KEYS:
        p = paths["data"] / f"{key}.h5ad"
        if not p.exists():
            runlog.append({"key": key, "status": "missing"})
            continue
        print(f"analyze {key}", flush=True)
        try:
            info, assoc, ct = analyze_one(key, p, tables, figs)
            runlog.append({"key": key, "status": "ok", "n_obs": info.get("n_obs"),
                           "tacstd2": info.get("tacstd2_name"), "cldn4": info.get("cldn4_name")})
            if assoc is not None and len(assoc):
                all_assoc.append(assoc)
            if ct is not None and len(ct):
                all_ct.append(ct)
        except Exception as e:
            runlog.append({"key": key, "status": "error", "error": str(e),
                           "trace": traceback.format_exc()[-800:]})
            print("ERROR", key, e, flush=True)
    if all_assoc:
        cat = pd.concat(all_assoc, ignore_index=True)
        cat["q"] = bh_fdr(cat["p"].to_numpy())
        cat.to_csv(tables / "scrna_associations_all.tsv", sep="\t", index=False)
    if all_ct:
        pd.concat(all_ct, ignore_index=True).to_csv(tables / "scrna_celltype_all.tsv", sep="\t", index=False)
    write_json(paths["notes"] / "scrna_runlog.json", runlog)
    print(json.dumps(runlog, indent=2))


if __name__ == "__main__":
    main()
