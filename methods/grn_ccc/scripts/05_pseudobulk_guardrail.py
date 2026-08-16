# %% [markdown]
# # 05 - Pseudobulk guardrail: the test that actually supports the claim
#
# **Why this is the primary analysis, not an afterthought.** LIANA/CellChat/
# NicheNet/SCENIC are hypothesis generators built on co-expression and generic
# priors. The clean, hard-to-overclaim question is direct: *within epithelium,
# do TACSTD2/CLDN4-high cells express fewer T-cell-recruitment chemokines than
# low cells, at the PATIENT level?* Pseudobulk + a patient-level model answers
# that without treating cells as independent replicates.
#
# Aggregate counts to (patient x epi_state), then compare chemokine expression
# with a paired/mixed design. Cells are not replicates; patients are.

# %%
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
P = CFG["params"]
IN_H5AD = "data/lung_ici_preprocessed.h5ad"
OUTDIR = Path("results/pseudobulk"); OUTDIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
epi = adata[adata.obs[P["celltype_key"]].astype(str).str.contains("pithel")].copy()
epi.X = epi.layers["counts"] if "counts" in epi.layers else epi.X

recr = CFG["t_cell_recruitment"]
chemokines = sorted({g for ax in recr.values() for g in ax["ligands"]})
chemokines = [g for g in chemokines if g in epi.var_names]

# %% [markdown]
# ## 1. Sum counts to pseudobulk (patient x state), keep library sizes

# %%
grp = epi.obs.groupby([P["sample_key"], "epi_state"], observed=True)
pb_rows, meta_rows = [], []
X = epi.X.tocsr() if hasattr(epi.X, "tocsr") else np.asarray(epi.X)
for (pat, state), idx in grp.indices.items():
    if state == "non_epithelial":
        continue
    sub = X[idx]
    total = np.asarray(sub.sum(axis=0)).ravel()
    pb_rows.append(total)
    meta_rows.append({"patient": pat, "state": state,
                      "n_cells": len(idx), "libsize": total.sum()})
pb = pd.DataFrame(pb_rows, columns=epi.var_names)
meta = pd.DataFrame(meta_rows)
# Drop tiny pseudobulk samples (unstable): require >= 20 cells.
keep = meta["n_cells"] >= 20
pb, meta = pb[keep].reset_index(drop=True), meta[keep].reset_index(drop=True)

# %% [markdown]
# ## 2. CPM-normalize and test high vs low per chemokine (paired by patient)
# We use CPM + log for a quick, transparent readout, and a paired Wilcoxon on
# patients that have BOTH states. For a publication use limma-voom / DESeq2 with
# `~ patient + state` (see the R snippet in the playbook) -- same idea, better
# variance modeling.

# %%
cpm = pb.div(meta["libsize"].values, axis=0) * 1e6
logcpm = np.log1p(cpm)
tbl = pd.concat([meta[["patient", "state", "n_cells"]], logcpm[chemokines]], axis=1)

from scipy.stats import wilcoxon
rows = []
for ck in chemokines:
    wide = tbl.pivot_table(index="patient", columns="state", values=ck)
    both = wide.dropna(subset=["TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low"])
    if len(both) >= P["min_patients_per_group"]:
        hi, lo = both["TACSTD2_CLDN4_high"], both["TACSTD2_CLDN4_low"]
        try:
            _, pval = wilcoxon(hi, lo)
        except ValueError:
            pval = np.nan
        rows.append((ck, len(both), float(hi.mean()), float(lo.mean()),
                     float((hi - lo).mean()), pval))
out = pd.DataFrame(rows, columns=["chemokine", "n_paired", "mean_logcpm_high",
                                  "mean_logcpm_low", "mean_delta_high_minus_low",
                                  "pval"])
if len(out):
    from statsmodels.stats.multitest import multipletests
    out["padj"] = multipletests(out["pval"].fillna(1), method="fdr_bh")[1]
out = out.sort_values("mean_delta_high_minus_low")
out.to_csv(OUTDIR / "chemokine_high_vs_low_pseudobulk.csv", index=False)
pb.assign(**meta.to_dict("series")).to_csv(OUTDIR / "pseudobulk_counts.csv", index=False)
print("Negative delta = LOWER recruitment chemokine in TACSTD2/CLDN4-high "
      "epithelium (supports the hypothesis):")
print(out.to_string(index=False))

# %% [markdown]
# ## 3. Fraction expressing (dropout sanity check)
# "Reduced signal" must not be pure dropout. Report the fraction of cells
# expressing each chemokine per state; a real difference should show in both
# level (above) and fraction (here).

# %%
frac_rows = []
for state in ["TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low"]:
    sub = epi[epi.obs["epi_state"] == state]
    Xs = sub.X.tocsr() if hasattr(sub.X, "tocsr") else np.asarray(sub.X)
    for ck in chemokines:
        j = sub.var_names.get_loc(ck)
        col = Xs[:, j]
        frac = float((col > 0).sum()) / sub.n_obs if sub.n_obs else np.nan
        frac_rows.append((state, ck, frac))
frac = (pd.DataFrame(frac_rows, columns=["state", "chemokine", "frac_expr"])
        .pivot(index="chemokine", columns="state", values="frac_expr"))
frac.to_csv(OUTDIR / "chemokine_fraction_expressing.csv")
print(frac)
print(f"[done] pseudobulk guardrail in {OUTDIR}")
