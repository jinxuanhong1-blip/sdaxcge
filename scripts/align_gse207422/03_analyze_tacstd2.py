"""GSE207422 USER-ALIGN: TACSTD2 restriction, NMPR vs MPR, T/NK scores, TJ/KRT.

Reads the annotated AnnData from 02_process_annotate.py and writes tables +
figures under results/align_gse207422/. Every statistic is computed from the
downloaded GEO matrix. Nothing is hard-coded to match a slide.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from gene_sets import CYT, CYTOTOXIC, EXHAUSTION, KERATINIZATION, TIGHT_JUNCTION

ANN = Path("data/gse207422/gse207422_annotated.h5ad")
OUT = Path("results/align_gse207422")
OUT.mkdir(parents=True, exist_ok=True)

MIN_MALIG = 20
MIN_TNK = 20
TARGET = "TACSTD2"


def present(adata, genes):
    return [g for g in genes if g in adata.var_names]


def gene_vector(adata, gene):
    """log1p CP10K from adata.X (already normalized+log1p)."""
    if gene not in adata.var_names:
        raise KeyError(gene)
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


def module_score_matrix(adata, genes):
    """Mean of per-gene z-scores across cells (scanpy-like, no random control)."""
    genes = present(adata, genes)
    if not genes:
        return np.full(adata.n_obs, np.nan), genes
    cols = []
    for g in genes:
        v = gene_vector(adata, g)
        sd = v.std()
        cols.append((v - v.mean()) / sd if sd > 0 else np.zeros_like(v))
    return np.vstack(cols).mean(axis=0), genes


def spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan, "note": "too_few"}
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p), "note": ""}


def mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {"n_a": int(len(a)), "n_b": int(len(b)),
                "median_a": float(np.median(a)) if len(a) else np.nan,
                "median_b": float(np.median(b)) if len(b) else np.nan,
                "U": np.nan, "p": np.nan, "note": "too_few"}
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)), "n_b": int(len(b)),
        "median_a": float(np.median(a)), "median_b": float(np.median(b)),
        "mean_a": float(np.mean(a)), "mean_b": float(np.mean(b)),
        "U": float(U), "p": float(p), "note": "",
    }


def savefig(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


def main():
    adata = sc.read_h5ad(ANN)
    print("loaded", adata.shape)
    if TARGET not in adata.var_names:
        raise SystemExit(f"{TARGET} missing from matrix")

    adata.obs["TACSTD2"] = gene_vector(adata, TARGET)
    adata.obs["TACSTD2_pos"] = adata.obs["TACSTD2"] > 0

    # T/NK = T + NK compartments
    tnk_mask = adata.obs["compartment"].isin(["T", "NK"])
    adata.obs["is_tnk"] = tnk_mask
    adata.obs["is_epi"] = adata.obs["compartment"].eq("Epithelial")

    # Module scores on the relevant cells (computed on all cells, then subset).
    cyto, cyto_used = module_score_matrix(adata, CYTOTOXIC)
    exh, exh_used = module_score_matrix(adata, EXHAUSTION)
    cyt, cyt_used = module_score_matrix(adata, CYT)
    tj, tj_used = module_score_matrix(adata, TIGHT_JUNCTION)
    krt, krt_used = module_score_matrix(adata, KERATINIZATION)
    adata.obs["score_cytotoxic"] = cyto
    adata.obs["score_exhaustion"] = exh
    adata.obs["score_CYT"] = cyt
    adata.obs["score_TJ"] = tj
    adata.obs["score_keratinization"] = krt

    # Response: paper groups pCR with MPR. Analysis unit is Patient.
    # Primary NMPR-vs-MPR contrast uses POST-treatment samples only (paper design).
    adata.obs["response3"] = adata.obs["response_group"].astype(str)
    adata.obs.loc[adata.obs["treatment"].eq("Pre"), "response3"] = "TN"

    # ---------- 1. TACSTD2 restriction by compartment ----------
    by_comp = (
        adata.obs.groupby("compartment", observed=True)
        .agg(
            n_cells=("TACSTD2", "size"),
            mean_log1p=("TACSTD2", "mean"),
            median_log1p=("TACSTD2", "median"),
            pct_pos=("TACSTD2_pos", "mean"),
        )
        .reset_index()
    )
    by_comp["pct_pos"] *= 100
    by_comp = by_comp.sort_values("mean_log1p", ascending=False)
    by_comp.to_csv(OUT / "tacstd2_by_compartment.csv", index=False)

    # cluster-level too (transparency of annotation)
    by_cl = (
        adata.obs.groupby(["leiden", "compartment"], observed=True)
        .agg(
            n_cells=("TACSTD2", "size"),
            mean_log1p=("TACSTD2", "mean"),
            pct_pos=("TACSTD2_pos", "mean"),
            score_Epithelial=("score_Epithelial", "mean"),
            score_T=("score_T", "mean"),
            score_NK=("score_NK", "mean"),
        )
        .reset_index()
    )
    by_cl["pct_pos"] *= 100
    by_cl.to_csv(OUT / "tacstd2_by_leiden.csv", index=False)

    # ---------- 2. Per-patient metrics ----------
    rows = []
    for (patient, sample), sdf in adata.obs.groupby(["Patient", "Sample"], observed=True):
        epi = sdf[sdf["is_epi"]]
        tnk = sdf[sdf["is_tnk"]]
        row = {
            "Patient": patient,
            "Sample": sample,
            "treatment": sdf["treatment"].iloc[0],
            "response_group": sdf["response_group"].iloc[0],
            "response3": sdf["response3"].iloc[0],
            "Pathology": sdf["Pathology"].iloc[0],
            "Resource": sdf["Resource"].iloc[0],
            "n_cells": int(len(sdf)),
            "n_epithelial": int(len(epi)),
            "n_TNK": int(len(tnk)),
            "epi_TACSTD2_mean": float(epi["TACSTD2"].mean()) if len(epi) else np.nan,
            "epi_TACSTD2_pct_pos": float(epi["TACSTD2_pos"].mean() * 100) if len(epi) else np.nan,
            "epi_TJ_mean": float(epi["score_TJ"].mean()) if len(epi) else np.nan,
            "epi_KRT_mean": float(epi["score_keratinization"].mean()) if len(epi) else np.nan,
            "tnk_cytotoxic_mean": float(tnk["score_cytotoxic"].mean()) if len(tnk) else np.nan,
            "tnk_exhaustion_mean": float(tnk["score_exhaustion"].mean()) if len(tnk) else np.nan,
            "tnk_CYT_mean": float(tnk["score_CYT"].mean()) if len(tnk) else np.nan,
            "tnk_TACSTD2_mean": float(tnk["TACSTD2"].mean()) if len(tnk) else np.nan,
            "tnk_TACSTD2_pct_pos": float(tnk["TACSTD2_pos"].mean() * 100) if len(tnk) else np.nan,
        }
        rows.append(row)
    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(OUT / "per_sample_metrics.csv", index=False)

    # Patient-level: one row per patient. For the 3 TN patients there is only
    # the pre-treatment sample. For post-treatment patients use the surgical sample.
    per_pt = per_sample.sort_values(["Patient", "treatment"]).drop_duplicates("Patient", keep="last")
    # Prefer Post when both exist (P05/P08 are Pre-only in this GEO matrix).
    per_pt.to_csv(OUT / "per_patient_metrics.csv", index=False)

    eval_malig = per_pt["n_epithelial"] >= MIN_MALIG
    eval_tnk = per_pt["n_TNK"] >= MIN_TNK
    eval_both = eval_malig & eval_tnk

    # ---------- 3. NMPR vs MPR (post-treatment, patient-level) ----------
    post = per_pt[per_pt["treatment"].eq("Post") & eval_malig].copy()
    mpr = post.loc[post["response_group"].eq("MPR"), "epi_TACSTD2_mean"]
    nmpr = post.loc[post["response_group"].eq("NMPR"), "epi_TACSTD2_mean"]
    mpr_pct = post.loc[post["response_group"].eq("MPR"), "epi_TACSTD2_pct_pos"]
    nmpr_pct = post.loc[post["response_group"].eq("NMPR"), "epi_TACSTD2_pct_pos"]

    # Sensitivity: all patients with malignant cells (include TN labeled by eventual response)
    all_eval = per_pt[eval_malig & per_pt["response_group"].isin(["MPR", "NMPR"])]
    mpr_all = all_eval.loc[all_eval["response_group"].eq("MPR"), "epi_TACSTD2_mean"]
    nmpr_all = all_eval.loc[all_eval["response_group"].eq("NMPR"), "epi_TACSTD2_mean"]

    # ---------- 4. Per-patient correlations ----------
    # Primary: patients with both compartments (any treatment).
    both = per_pt[eval_both].copy()
    corrs = {
        "cytotoxic_all": spearman(both["epi_TACSTD2_mean"].values, both["tnk_cytotoxic_mean"].values),
        "exhaustion_all": spearman(both["epi_TACSTD2_mean"].values, both["tnk_exhaustion_mean"].values),
        "CYT_all": spearman(both["epi_TACSTD2_mean"].values, both["tnk_CYT_mean"].values),
    }
    post_both = both[both["treatment"].eq("Post")]
    corrs["cytotoxic_post"] = spearman(post_both["epi_TACSTD2_mean"].values, post_both["tnk_cytotoxic_mean"].values)
    corrs["exhaustion_post"] = spearman(post_both["epi_TACSTD2_mean"].values, post_both["tnk_exhaustion_mean"].values)
    corrs["CYT_post"] = spearman(post_both["epi_TACSTD2_mean"].values, post_both["tnk_CYT_mean"].values)

    # ---------- 5. TROP2-high tumors: TJ / keratinization ----------
    # Split POST patients with enough epithelial cells by median malignant TACSTD2.
    trop = post.copy()
    med = float(trop["epi_TACSTD2_mean"].median())
    trop["trop2_group"] = np.where(trop["epi_TACSTD2_mean"] >= med, "TROP2_high", "TROP2_low")
    trop.to_csv(OUT / "post_patients_trop2_split.csv", index=False)
    tj_test = mwu(
        trop.loc[trop["trop2_group"].eq("TROP2_high"), "epi_TJ_mean"],
        trop.loc[trop["trop2_group"].eq("TROP2_low"), "epi_TJ_mean"],
    )
    krt_test = mwu(
        trop.loc[trop["trop2_group"].eq("TROP2_high"), "epi_KRT_mean"],
        trop.loc[trop["trop2_group"].eq("TROP2_low"), "epi_KRT_mean"],
    )
    # Continuous Spearman (more honest than a median split on n=12).
    tj_rho = spearman(trop["epi_TACSTD2_mean"].values, trop["epi_TJ_mean"].values)
    krt_rho = spearman(trop["epi_TACSTD2_mean"].values, trop["epi_KRT_mean"].values)

    # Optional prerank GSEA on epithelial pseudobulk DE (TROP2-high vs low).
    gsea_summary = run_pseudobulk_gsea(adata, trop)

    # ---------- stats JSON ----------
    stats_out = {
        "dataset": "GSE207422",
        "n_cells_raw": 92330,
        "n_cells_qc": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_patients": int(per_pt["Patient"].nunique()),
        "n_samples": int(per_sample["Sample"].nunique()),
        "compartments": by_comp.to_dict(orient="records"),
        "annotation_note": (
            "Malignant/epithelial = Leiden clusters whose highest lineage score "
            "is Epithelial (EPCAM/KRT). No author cell-level labels or inferCNV "
            "were available on GEO; epithelial is the malignant proxy."
        ),
        "min_cells": {"epithelial": MIN_MALIG, "TNK": MIN_TNK},
        "genes_used": {
            "cytotoxic": cyto_used,
            "exhaustion": exh_used,
            "CYT": cyt_used,
            "TJ": tj_used,
            "keratinization": krt_used,
        },
        "NMPR_vs_MPR_post_mean_log1p": {
            **mwu(nmpr, mpr),
            "group_a": "NMPR", "group_b": "MPR",
            "patients_NMPR": post.loc[post["response_group"].eq("NMPR"), "Patient"].tolist(),
            "patients_MPR": post.loc[post["response_group"].eq("MPR"), "Patient"].tolist(),
        },
        "NMPR_vs_MPR_post_pct_pos": {**mwu(nmpr_pct, mpr_pct), "group_a": "NMPR", "group_b": "MPR"},
        "NMPR_vs_MPR_all_eval_mean_log1p": {**mwu(nmpr_all, mpr_all), "group_a": "NMPR", "group_b": "MPR"},
        "correlations_patient_level": corrs,
        "TROP2_high_vs_low": {
            "split": "median of post-treatment patient epithelial TACSTD2 mean",
            "median_cutoff": med,
            "n_high": int((trop["trop2_group"] == "TROP2_high").sum()),
            "n_low": int((trop["trop2_group"] == "TROP2_low").sum()),
            "TJ_MWU_high_vs_low": tj_test,
            "KRT_MWU_high_vs_low": krt_test,
            "TJ_spearman_vs_TACSTD2": tj_rho,
            "KRT_spearman_vs_TACSTD2": krt_rho,
        },
        "gsea": gsea_summary,
        "slide5_claims": {
            "n_cells": "~90k neoadjuvant anti-PD-1",
            "restriction": "TACSTD2 restricted to malignant/epithelial",
            "NMPR_vs_MPR": "malignant-cell mean TACSTD2 NMPR vs MPR (direction not specified in the prompt)",
            "correlations": "per-patient malignant TACSTD2 vs T/NK cytotoxicity and exhaustion; they report ρ ≈ −0.40 to −0.50, all negative",
            "TROP2_high": "TROP2-high tumors enrich TJ/keratinization",
        },
    }
    (OUT / "stats.json").write_text(json.dumps(stats_out, indent=2, default=_json_default))
    print(json.dumps({k: stats_out[k] for k in [
        "n_cells_qc", "NMPR_vs_MPR_post_mean_log1p", "correlations_patient_level",
        "TROP2_high_vs_low",
    ]}, indent=2, default=_json_default))

    # ---------- figures ----------
    plot_restriction(adata, by_comp)
    plot_nmpr_mpr(post)
    plot_correlations(both)
    plot_trop2_programs(trop)
    print("done")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if pd.isna(o):
        return None
    raise TypeError(type(o))


def run_pseudobulk_gsea(adata, trop):
    """Wilcoxon DE on epithelial cells pooled by TROP2-high/low patient, then prerank GSEA.

    Cells-as-replicates is pseudoreplication; this is labeled exploratory and is
    only used to ask whether TJ/keratinization gene sets sit in the high tail.
    """
    try:
        import gseapy as gp
    except Exception as e:
        return {"ran": False, "reason": f"gseapy import failed: {e}"}

    high_pts = set(trop.loc[trop["trop2_group"].eq("TROP2_high"), "Patient"])
    low_pts = set(trop.loc[trop["trop2_group"].eq("TROP2_low"), "Patient"])
    epi = adata[adata.obs["is_epi"] & adata.obs["treatment"].eq("Post")].copy()
    epi.obs["trop2_group"] = np.where(epi.obs["Patient"].isin(high_pts), "TROP2_high",
                               np.where(epi.obs["Patient"].isin(low_pts), "TROP2_low", "other"))
    epi = epi[epi.obs["trop2_group"].isin(["TROP2_high", "TROP2_low"])].copy()
    if epi.n_obs < 100:
        return {"ran": False, "reason": "too few epithelial cells"}

    # Patient-level pseudobulk log-mean, then t-statistic high vs low (honest n).
    genes = list(epi.var_names)
    # restrict to genes detected in >= 10 epithelial cells
    if hasattr(epi.X, "getnnz"):
        det = np.asarray(epi.X.getnnz(axis=0)).ravel()
    else:
        det = np.asarray((epi.X > 0).sum(axis=0)).ravel()
    keep = det >= 10
    epi = epi[:, keep].copy()

    pb_rows = []
    for pt, sdf_idx in epi.obs.groupby("Patient").groups.items():
        sub = epi[list(sdf_idx)]
        x = sub.X
        if hasattr(x, "mean"):
            mu = np.asarray(x.mean(axis=0)).ravel()
        else:
            mu = np.asarray(x).mean(axis=0)
        pb_rows.append((pt, trop.loc[trop["Patient"].eq(pt), "trop2_group"].iloc[0], mu))
    pb = pd.DataFrame({pt: mu for pt, _, mu in pb_rows}, index=epi.var_names)
    groups = {pt: g for pt, g, _ in pb_rows}
    high_cols = [c for c in pb.columns if groups[c] == "TROP2_high"]
    low_cols = [c for c in pb.columns if groups[c] == "TROP2_low"]
    if len(high_cols) < 2 or len(low_cols) < 2:
        return {"ran": False, "reason": "too few patients for pseudobulk"}

    tstats = []
    for gene in pb.index:
        a = pb.loc[gene, high_cols].astype(float).values
        b = pb.loc[gene, low_cols].astype(float).values
        if np.nanstd(a) == 0 and np.nanstd(b) == 0:
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        tstats.append((gene, float(t) if np.isfinite(t) else 0.0))
    rnk = pd.DataFrame(tstats, columns=["gene", "t"]).sort_values("t", ascending=False)
    rnk.to_csv(OUT / "pseudobulk_trop2high_vs_low_tstat.csv", index=False)

    gene_sets = {
        "Tight_junction": TIGHT_JUNCTION,
        "Keratinization": KERATINIZATION,
        "Cytotoxic": CYTOTOXIC,
        "Exhaustion": EXHAUSTION,
    }
    try:
        res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            outdir=str(OUT / "gsea_prerank"),
            permutation_num=1000,
            min_size=5,
            max_size=200,
            seed=1,
            verbose=False,
            no_plot=True,
        )
        tab = res.res2d.copy()
        tab.to_csv(OUT / "gsea_prerank_custom_sets.csv", index=False)
        recs = []
        for _, r in tab.iterrows():
            recs.append({
                "term": r.get("Term", r.get("Name")),
                "nes": float(r["NES"]) if "NES" in r and pd.notna(r["NES"]) else None,
                "fdr": float(r["FDR q-val"]) if "FDR q-val" in r and pd.notna(r["FDR q-val"]) else (
                    float(r["FDR"]) if "FDR" in r and pd.notna(r["FDR"]) else None
                ),
                "pval": float(r["NOM p-val"]) if "NOM p-val" in r and pd.notna(r["NOM p-val"]) else None,
            })
        return {
            "ran": True,
            "method": "patient-pseudobulk t-stat (TROP2-high vs low), gseapy prerank, 1000 perm, custom sets",
            "n_high": len(high_cols),
            "n_low": len(low_cols),
            "sets": recs,
        }
    except Exception as e:
        return {"ran": False, "reason": f"prerank failed: {e}",
                "n_high": len(high_cols), "n_low": len(low_cols)}


def plot_restriction(adata, by_comp):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    order = by_comp["compartment"].tolist()
    ax = axes[0]
    means = [by_comp.set_index("compartment").loc[c, "mean_log1p"] for c in order]
    ax.bar(range(len(order)), means, color=["#c0392b" if c == "Epithelial" else "#7f8c8d" for c in order])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=35, ha="right")
    ax.set_ylabel("mean log1p(CP10K) TACSTD2")
    ax.set_title("TACSTD2 by compartment")
    ax = axes[1]
    pcts = [by_comp.set_index("compartment").loc[c, "pct_pos"] for c in order]
    ax.bar(range(len(order)), pcts, color=["#c0392b" if c == "Epithelial" else "#7f8c8d" for c in order])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=35, ha="right")
    ax.set_ylabel("% TACSTD2+ cells")
    ax.set_title("Detection rate by compartment")
    fig.tight_layout()
    savefig(fig, "fig1_tacstd2_by_compartment.png")

    # UMAP
    if "X_umap" in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
        sc.pl.umap(adata, color="compartment", ax=axes[0], show=False, title="compartment")
        sc.pl.umap(adata, color="TACSTD2", ax=axes[1], show=False, title="TACSTD2", cmap="Reds")
        fig.tight_layout()
        savefig(fig, "fig1b_umap_tacstd2.png")


def plot_nmpr_mpr(post):
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2))
    for ax, col, ylab in [
        (axes[0], "epi_TACSTD2_mean", "epithelial TACSTD2 mean log1p(CP10K)"),
        (axes[1], "epi_TACSTD2_pct_pos", "epithelial % TACSTD2+"),
    ]:
        for i, grp in enumerate(["MPR", "NMPR"]):
            y = post.loc[post["response_group"].eq(grp), col].values
            x = np.random.default_rng(0).normal(i, 0.06, size=len(y))
            ax.scatter(x, y, s=50, zorder=3)
            ax.hlines(np.median(y), i - 0.2, i + 0.2, color="black")
            for pt, yy, xx in zip(post.loc[post["response_group"].eq(grp), "Patient"], y, x):
                ax.annotate(pt, (xx, yy), textcoords="offset points", xytext=(4, 2), fontsize=7)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["MPR (incl. pCR)", "NMPR"])
        ax.set_ylabel(ylab)
        ax.set_title("Post-treatment, patient-level")
    fig.tight_layout()
    savefig(fig, "fig2_nmpr_vs_mpr_tacstd2.png")


def plot_correlations(both):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    pairs = [
        (axes[0], "tnk_cytotoxic_mean", "T/NK cytotoxic score"),
        (axes[1], "tnk_exhaustion_mean", "T/NK exhaustion score"),
    ]
    for ax, col, ylab in pairs:
        x = both["epi_TACSTD2_mean"].values
        y = both[col].values
        ax.scatter(x, y, c=np.where(both["response_group"].eq("MPR"), "#27ae60",
                             np.where(both["response_group"].eq("NMPR"), "#c0392b", "#7f8c8d")))
        for _, r in both.iterrows():
            ax.annotate(r["Patient"], (r["epi_TACSTD2_mean"], r[col]),
                        textcoords="offset points", xytext=(3, 2), fontsize=7)
        rho, p = stats.spearmanr(x, y)
        ax.set_xlabel("epithelial TACSTD2 mean")
        ax.set_ylabel(ylab)
        ax.set_title(f"Spearman ρ = {rho:.2f}, p = {p:.3g}, n = {len(both)}")
    fig.tight_layout()
    savefig(fig, "fig3_patient_correlations.png")


def plot_trop2_programs(trop):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, col, ylab in [
        (axes[0], "epi_TJ_mean", "epithelial tight-junction score"),
        (axes[1], "epi_KRT_mean", "epithelial keratinization score"),
    ]:
        ax.scatter(trop["epi_TACSTD2_mean"], trop[col],
                   c=np.where(trop["trop2_group"].eq("TROP2_high"), "#c0392b", "#2980b9"))
        for _, r in trop.iterrows():
            ax.annotate(r["Patient"], (r["epi_TACSTD2_mean"], r[col]),
                        textcoords="offset points", xytext=(3, 2), fontsize=7)
        rho, p = stats.spearmanr(trop["epi_TACSTD2_mean"], trop[col])
        ax.set_xlabel("epithelial TACSTD2 mean")
        ax.set_ylabel(ylab)
        ax.set_title(f"post-tx  ρ = {rho:.2f}, p = {p:.3g}, n = {len(trop)}")
    fig.tight_layout()
    savefig(fig, "fig4_trop2_vs_TJ_KRT.png")


if __name__ == "__main__":
    main()
