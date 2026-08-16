# %% [markdown]
# # 04b - SCENIC downstream: state-defining regulons & the chemokine question
#
# Two jobs:
#  (0) export a loom of epithelium (subsampled) for `04_pyscenic.sh`;
#  (1) after AUCell, find regulons whose activity separates TACSTD2/CLDN4-high
#      vs -low epithelium, take a CONSENSUS across seeds, and check where the
#      T-cell-recruitment chemokines sit relative to those regulons.
#
# **Overclaim guard.** A regulon activity difference is *association*. "TF X
# represses CXCL10 in the high state" requires: (a) CXCL10 anticorrelated with
# the regulon, (b) an X motif in the CXCL10 locus (cisTarget/ATAC), and (c) ideally
# perturbation/orthogonal evidence. This script produces (a)+(b) as hypotheses
# only. See playbook §SCENIC overclaims.

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
OUTDIR = Path("results/scenic"); OUTDIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Step 0. Export epithelium loom for GRN inference (run once, then 04_pyscenic.sh)
# Subsample to keep GRNBoost2 tractable; SCENIC needs raw/normalized counts.

# %%
def export_loom(max_cells=20000):
    import loompy  # noqa
    adata = sc.read_h5ad(IN_H5AD)
    epi = adata[adata.obs[P["celltype_key"]].astype(str).str.contains("pithel")].copy()
    if epi.n_obs > max_cells:
        idx = np.random.RandomState(P["random_seed"]).choice(
            epi.n_obs, max_cells, replace=False)
        epi = epi[idx].copy()
    # SCENIC expects genes x cells not required; loom via scanpy write_loom is cells x genes.
    epi.X = epi.layers["counts"] if "counts" in epi.layers else epi.X
    epi.write_loom("data/epithelium_for_scenic.loom", write_obsm_varm=False)
    print(f"[loom] {epi.n_obs} cells x {epi.n_vars} genes -> data/epithelium_for_scenic.loom")

# export_loom()   # <- uncomment for the first run

# %% [markdown]
# ## Step 1. Load AUCell results from >=3 seeds and build a consensus
# Keep only regulons recovered by a majority of seeds (robustness against
# GRNBoost2 stochasticity).

# %%
def load_auc(loom_path):
    import loompy
    with loompy.connect(loom_path, mode="r") as ds:
        auc = pd.DataFrame(ds.ca["RegulonsAUC"]) if "RegulonsAUC" in ds.ca else None
        if auc is None:  # pyscenic stores AUC as a structured array
            arr = ds.ca["RegulonsAUC"]
            auc = pd.DataFrame(arr.tolist(), columns=arr.dtype.names)
        auc.index = ds.ca["CellID"]
    return auc

seed_looms = sorted(OUTDIR.glob("aucell.seed*.loom"))
if not seed_looms:
    print("[info] no AUCell looms yet; run 04_pyscenic.sh for >=3 seeds first.")
else:
    aucs = {p.stem: load_auc(p) for p in seed_looms}
    common = set.intersection(*[set(a.columns) for a in aucs.values()])
    print(f"{len(common)} regulons common to all {len(aucs)} seeds.")
    auc = pd.concat([a[list(common)] for a in aucs.values()]).groupby(level=0).mean()

    adata = sc.read_h5ad(IN_H5AD)
    obs = adata.obs.loc[auc.index]
    hi = obs["epi_state"] == "TACSTD2_CLDN4_high"
    lo = obs["epi_state"] == "TACSTD2_CLDN4_low"

    # %% [markdown]
    # ## Step 2. Which regulons define the high state? (patient-aware)
    # Use patient-level means to avoid pseudoreplication, then Wilcoxon.
    # %%
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    pat = obs[P["sample_key"]].values
    rows = []
    for reg in common:
        v = auc[reg].values
        pm = pd.DataFrame({"pat": pat, "state": np.where(hi, "hi", np.where(lo, "lo", "x")),
                           "auc": v})
        pm = pm[pm.state != "x"].groupby(["pat", "state"])["auc"].mean().unstack()
        if {"hi", "lo"}.issubset(pm.columns):
            a, b = pm["hi"].dropna(), pm["lo"].dropna()
            if len(a) >= P["min_patients_per_group"] and len(b) >= P["min_patients_per_group"]:
                _, pval = mannwhitneyu(a, b)
                rows.append((reg, a.median(), b.median(), a.median() - b.median(), pval))
    reg_de = pd.DataFrame(rows, columns=["regulon", "median_hi", "median_lo",
                                         "delta", "pval"]).sort_values("delta")
    if len(reg_de):
        reg_de["padj"] = multipletests(reg_de["pval"], method="fdr_bh")[1]
    reg_de.to_csv(OUTDIR / "regulons_high_vs_low_patientlevel.csv", index=False)
    print(reg_de.head(15).to_string(index=False))

    # %% [markdown]
    # ## Step 3. Are recruitment chemokines targets of the high-state regulons,
    # and are they anticorrelated with them? (hypothesis for repression)
    # %%
    regulons = pd.read_csv(sorted(OUTDIR.glob("regulons.seed*.csv"))[0],
                           header=[0, 1], index_col=[0, 1], skipinitialspace=True) \
        if list(OUTDIR.glob("regulons.seed*.csv")) else None
    chemokines = sorted({g for ax in CFG["t_cell_recruitment"].values()
                         for g in ax["ligands"]})
    # Correlate each state-up regulon's AUC with chemokine expression in epithelium.
    epi = adata[adata.obs_names.isin(auc.index)]
    up_regs = reg_de.sort_values("delta", ascending=False).head(10)["regulon"].tolist()
    corr_rows = []
    for reg in up_regs:
        a = auc.loc[epi.obs_names, reg].values
        for ck in [c for c in chemokines if c in epi.var_names]:
            x = np.asarray(epi[:, ck].X.todense()).ravel() if hasattr(epi.X, "todense") \
                else np.asarray(epi[:, ck].X).ravel()
            if x.std() > 0:
                r = np.corrcoef(a, x)[0, 1]
                corr_rows.append((reg, ck, r))
    corr = pd.DataFrame(corr_rows, columns=["regulon", "chemokine", "pearson_r"])
    corr.to_csv(OUTDIR / "regulon_chemokine_correlation.csv", index=False)
    print("Negative r = chemokine LOW where the high-state regulon is active "
          "(hypothesis of repression, NOT proof):")
    print(corr.sort_values("pearson_r").head(15).to_string(index=False))

print(f"[done] SCENIC downstream in {OUTDIR}")
