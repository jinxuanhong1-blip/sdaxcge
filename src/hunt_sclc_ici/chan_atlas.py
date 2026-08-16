"""Chan et al. 2021 SCLC single-cell atlas: TACSTD2 / CLDN4 vs immune.

Source: cellxgene collection 62e8f058-9c37-48bc-9200-e767f318a8ec
        combined object (147,137 cells; SCLC + LUAD + normal).
This is PUBLIC. Results here are REAL.

Honesty notes:
- Chan annotates tumor cells as SCLC-A / SCLC-N / SCLC-P (plus a small
  'SCLC' leftover). There is NO SCLC-I tumor-cell label. SCLC-I is a
  bulk / TME concept (Gay 2021). We therefore treat 'immune' as
  compartment + patient-level immune fraction, not as a fourth tumor
  subtype.
- Expression is the cellxgene-normalized matrix (log-scale, not raw UMI).
  Detection = value > 0 on that scale. We do not re-normalize.
- Cells are not independent. Primary inference is at the *patient*
  (donor) level; cell-level rates are descriptive only.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from . import bulk_george as bg
from . import signatures as sig

DATA = bg.DATA
TABLES = bg.TABLES
LOGS = bg.LOGS
H5AD = DATA / "chan_combined.h5ad"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler(LOGS / "chan_atlas.log", mode="w"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("chan_atlas")

# Genes pulled from the combined object. Keep this list short: we extract
# one column at a time from a backed h5ad (~1.5 GB).
GENES = (
    sig.GENES_OF_INTEREST
    + sig.CONTEXT_SURFACE_TARGETS
    + ["ASCL1", "NEUROD1", "POU2F3", "YAP1", "CD8A", "CD8B",
       "GZMA", "PRF1", "CD274", "PDCD1", "CXCL9", "CXCL10"]
)

# Minimum tumor epithelial cells to keep a donor in patient-level stats
MIN_TUMOR_CELLS = 30
MIN_TOTAL_CELLS = 80


def _symbol_to_ensembl(adata: ad.AnnData) -> dict[str, str]:
    names = adata.var["feature_name"].astype(str)
    # feature_name can be a Categorical; values are gene symbols
    mapping = {}
    for ens, sym in zip(adata.var_names, names):
        mapping.setdefault(str(sym), str(ens))
    return mapping


def extract_gene_matrix(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    """Return cells x genes from a backed h5ad (one column at a time)."""
    sym2ens = _symbol_to_ensembl(adata)
    missing = [g for g in genes if g not in sym2ens]
    if missing:
        log.warning("Genes not in Chan object: %s", missing)
    present = [g for g in genes if g in sym2ens]
    cols = {}
    var_index = list(adata.var_names)
    ens_to_i = {e: i for i, e in enumerate(var_index)}
    for g in present:
        ens = sym2ens[g]
        i = ens_to_i[ens]
        col = adata.X[:, i]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        cols[g] = col.astype(np.float32)
        log.info("Extracted %s: nonzero %.1f%%, max %.2f",
                 g, 100 * (cols[g] > 0).mean(), float(cols[g].max()))
    return pd.DataFrame(cols, index=adata.obs.index)


def _bh(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def compartment_table(obs: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """Detection rate + mean (all cells, and among detected) by compartment.

    Restricted to SCLC samples (histo == 'SCLC').
    """
    sclc = obs["histo"].astype(str) == "SCLC"
    rows = []
    for comp, mask_comp in [
        ("Epithelial", obs["cell_type_coarse"].astype(str) == "Epithelial"),
        ("Lymphoid", obs["cell_type_coarse"].astype(str) == "Lymphoid"),
        ("Myeloid", obs["cell_type_coarse"].astype(str) == "Myeloid"),
        ("Mesenchymal", obs["cell_type_coarse"].astype(str) == "Mesenchymal"),
    ]:
        m = sclc & mask_comp
        n = int(m.sum())
        rec = {"compartment": comp, "n_cells": n}
        for g in sig.GENES_OF_INTEREST + sig.CONTEXT_SURFACE_TARGETS:
            if g not in expr.columns:
                continue
            v = expr.loc[m, g].to_numpy()
            rec[f"{g}_detect_frac"] = float((v > 0).mean()) if n else np.nan
            rec[f"{g}_mean"] = float(v.mean()) if n else np.nan
            rec[f"{g}_mean_detected"] = float(v[v > 0].mean()) if (v > 0).any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def tumor_subtype_table(obs: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """TACSTD2/CLDN4 in SCLC tumor epithelial cells by Chan subtype label."""
    sclc_epi = (
        (obs["histo"].astype(str) == "SCLC")
        & (obs["cell_type_coarse"].astype(str) == "Epithelial")
    )
    rows = []
    for lab in ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC"]:
        m = sclc_epi & (obs["cell_type_med"].astype(str) == lab)
        n = int(m.sum())
        rec = {"chan_tumor_label": lab, "n_cells": n,
               "n_donors": int(obs.loc[m, "donor_id"].nunique()) if n else 0}
        for g in sig.GENES_OF_INTEREST + ["ASCL1", "NEUROD1", "POU2F3", "YAP1"]:
            if g not in expr.columns:
                continue
            v = expr.loc[m, g].to_numpy()
            rec[f"{g}_detect_frac"] = float((v > 0).mean()) if n else np.nan
            rec[f"{g}_mean"] = float(v.mean()) if n else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def patient_table(obs: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """One row per SCLC donor: tumor-cell ADC-target means vs immune fraction."""
    sclc = obs["histo"].astype(str) == "SCLC"
    rows = []
    for donor, idx in obs.loc[sclc].groupby("donor_id", observed=True).groups.items():
        o = obs.loc[idx]
        e = expr.loc[idx]
        n = len(o)
        tumor = o["cell_type_coarse"].astype(str) == "Epithelial"
        immune = o["cell_type_general"].astype(str) == "Immune"
        tcell = o["cell_type_med"].astype(str) == "T cell"
        n_tumor = int(tumor.sum())
        n_immune = int(immune.sum())
        rec = {
            "donor_id": str(donor),
            "n_cells": n,
            "n_tumor": n_tumor,
            "n_immune": n_immune,
            "n_tcell": int(tcell.sum()),
            "frac_immune": n_immune / n if n else np.nan,
            "frac_tcell": int(tcell.sum()) / n if n else np.nan,
            "treatment": str(o["treatment"].mode().iloc[0]) if len(o) else "",
            "procedure": str(o["procedure"].mode().iloc[0]) if len(o) else "",
        }
        # dominant Chan tumor-cell subtype (among epithelial)
        if n_tumor:
            rec["dominant_tumor_label"] = str(
                o.loc[tumor, "cell_type_med"].astype(str).value_counts().idxmax())
        else:
            rec["dominant_tumor_label"] = "no_tumor_cells"
        ici = "Immunotherapy" in rec["treatment"]
        rec["ever_immunotherapy"] = bool(ici)
        rec["treatment_naive"] = rec["treatment"] == "Naive"
        for g in sig.GENES_OF_INTEREST + ["ASCL1", "NEUROD1", "POU2F3", "YAP1", "CD8A"]:
            if g not in e.columns:
                continue
            if n_tumor:
                tv = e.loc[tumor, g].to_numpy()
                rec[f"tumor_{g}_mean"] = float(tv.mean())
                rec[f"tumor_{g}_detect"] = float((tv > 0).mean())
            else:
                rec[f"tumor_{g}_mean"] = np.nan
                rec[f"tumor_{g}_detect"] = np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["keep_for_stats"] = (
        (df["n_tumor"] >= MIN_TUMOR_CELLS) & (df["n_cells"] >= MIN_TOTAL_CELLS)
    )
    return df


def patient_correlations(pat: pd.DataFrame) -> pd.DataFrame:
    """Spearman: tumor-cell ADC-target vs patient immune fraction."""
    use = pat[pat["keep_for_stats"]].copy()
    rows = []
    axes = ["frac_immune", "frac_tcell"]
    # both mean (includes zeros) and detection fraction — TACSTD2 is sparse
    gene_cols = []
    for c in use.columns:
        if c.startswith("tumor_") and (c.endswith("_mean") or c.endswith("_detect")):
            gene_cols.append(c)
    for gcol in gene_cols:
        metric = "mean" if gcol.endswith("_mean") else "detect_frac"
        gene = gcol.replace("tumor_", "").replace("_mean", "").replace("_detect", "")
        for ax in axes:
            x = use[gcol].astype(float)
            y = use[ax].astype(float)
            m = x.notna() & y.notna()
            if m.sum() < 6:
                continue
            rho, p = stats.spearmanr(x[m], y[m])
            rows.append({
                "tumor_gene": gene,
                "metric": metric,
                "axis": ax,
                "n_donors": int(m.sum()),
                "spearman_rho": float(rho),
                "p_value": float(p),
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value_BH"] = _bh(df["p_value"].tolist())
        df = df.sort_values(["tumor_gene", "p_value"]).reset_index(drop=True)
    return df


def subtype_vs_target(pat: pd.DataFrame) -> pd.DataFrame:
    """Patient-level TACSTD2/CLDN4 by dominant Chan tumor label (A/N/P)."""
    use = pat[pat["keep_for_stats"]].copy()
    use = use[use["dominant_tumor_label"].isin(["SCLC-A", "SCLC-N", "SCLC-P"])]
    rows = []
    for g in sig.GENES_OF_INTEREST:
        col = f"tumor_{g}_mean"
        if col not in use.columns:
            continue
        rec = {"gene": g}
        groups = []
        for lab in ["SCLC-A", "SCLC-N", "SCLC-P"]:
            v = use.loc[use["dominant_tumor_label"] == lab, col].astype(float).dropna()
            rec[f"n_{lab}"] = int(v.size)
            rec[f"median_{lab}"] = float(v.median()) if v.size else np.nan
            if v.size:
                groups.append(v.to_numpy())
        if len(groups) >= 2 and all(len(v) >= 2 for v in groups):
            try:
                rec["kruskal_H"], rec["kruskal_p"] = stats.kruskal(*groups)
            except ValueError:
                rec["kruskal_H"], rec["kruskal_p"] = np.nan, np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def run() -> dict:
    if not H5AD.exists():
        raise FileNotFoundError(
            f"{H5AD} missing. Run results/hunt_sclc_ici/data/download_data.sh")
    log.info("Opening %s (backed)", H5AD)
    adata = ad.read_h5ad(H5AD, backed="r")
    log.info("AnnData %s cells x %s genes", adata.n_obs, adata.n_vars)

    expr = extract_gene_matrix(adata, GENES)
    obs = adata.obs.copy()
    # align
    expr.index = obs.index
    adata.file.close()

    comp = compartment_table(obs, expr)
    comp.to_csv(TABLES / "chan_compartment.tsv", sep="\t", index=False)
    log.info("Compartment table:\n%s", comp.to_string(index=False))

    subt = tumor_subtype_table(obs, expr)
    subt.to_csv(TABLES / "chan_tumor_subtype_cells.tsv", sep="\t", index=False)
    log.info("Tumor-cell subtype table:\n%s", subt.to_string(index=False))

    pat = patient_table(obs, expr)
    pat.to_csv(TABLES / "chan_patients.tsv", sep="\t", index=False)
    log.info("SCLC donors: %d total, %d kept for stats (min tumor cells=%d)",
             len(pat), int(pat["keep_for_stats"].sum()), MIN_TUMOR_CELLS)

    corr = patient_correlations(pat)
    corr.to_csv(TABLES / "chan_patient_correlations.tsv", sep="\t", index=False)
    log.info("Patient-level correlations:\n%s", corr.to_string(index=False))

    # Sensitivity: the TACSTD2–immune association is inflated by donors whose
    # dominant epithelial label is 'Epithelial Stroma' (not SCLC-A/N/P).
    sens_rows = []
    for label, mask in [
        ("all_kept", pat["keep_for_stats"]),
        ("ANP_tumor_dominant",
         pat["keep_for_stats"]
         & pat["dominant_tumor_label"].isin(["SCLC-A", "SCLC-N", "SCLC-P"])),
        ("treatment_naive", pat["keep_for_stats"] & pat["treatment_naive"]),
    ]:
        sub = patient_correlations(pat.assign(keep_for_stats=mask))
        if not sub.empty:
            sub.insert(0, "subset", label)
            sens_rows.append(sub)
    if sens_rows:
        sens = pd.concat(sens_rows, ignore_index=True)
        sens.to_csv(TABLES / "chan_patient_correlations_sensitivity.tsv",
                    sep="\t", index=False)
        log.info("Sensitivity correlations:\n%s",
                 sens[sens["tumor_gene"].isin(["TACSTD2", "CLDN4"])].to_string(index=False))

    grp = subtype_vs_target(pat)
    grp.to_csv(TABLES / "chan_patient_subtype_vs_target.tsv", sep="\t", index=False)
    log.info("Patient-level subtype vs target:\n%s", grp.to_string(index=False))

    summary = {
        "source": "Chan et al. 2021 Cancer Cell; cellxgene combined object",
        "n_cells_total": int(len(obs)),
        "n_cells_SCLC": int((obs["histo"].astype(str) == "SCLC").sum()),
        "n_SCLC_donors": int(pat.shape[0]),
        "n_SCLC_donors_kept": int(pat["keep_for_stats"].sum()),
        "note": "No SCLC-I tumor-cell label in Chan; immune is TME fraction.",
    }
    (TABLES / "chan_summary.json").write_text(json.dumps(summary, indent=2))
    return {"obs": obs, "expr": expr, "pat": pat, "comp": comp,
            "subt": subt, "corr": corr, "grp": grp}


if __name__ == "__main__":
    run()
