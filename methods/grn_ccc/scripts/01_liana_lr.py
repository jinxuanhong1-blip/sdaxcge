# %% [markdown]
# # 01 - LIANA: ligand-receptor communication (consensus)
#
# **Question this addresses.** Do TACSTD2/CLDN4-high epithelial *senders* show
# weaker inferred communication toward T/NK cells through T-cell-recruitment
# axes (CXCL9/10/11->CXCR3, CCL5->CCR5, CXCL16->CXCR6) than TACSTD2/CLDN4-low
# epithelium?
#
# **What LIANA does / does not do.** LIANA aggregates several LR methods
# (CellPhoneDB, NATMI, Connectome, logFC, SingleCellSignalR, geometric-mean)
# into a consensus rank. It scores **co-expression** of a ligand in the sender
# cluster and a receptor in the receiver cluster. It does **not** observe
# secretion, diffusion, spatial proximity, or downstream signaling, and its
# permutation p-values test *cluster-label specificity within one sample*, not
# reproducibility across patients. Treat output as ranked hypotheses. See the
# playbook §LIANA overclaims.
#
# Input: ../scripts/00_preprocess.py output (log-norm in .X, `epi_state` set).

# %%
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
P = CFG["params"]

IN_H5AD = "data/lung_ici_preprocessed.h5ad"
OUTDIR = Path("results/liana"); OUTDIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)

# %% [markdown]
# ## 1. Build the grouping used as sender/receiver identities
# Split epithelium into high/low and keep other lineages as-is, so LIANA can
# contrast `Epi_high -> T_cell` vs `Epi_low -> T_cell` directly.

# %%
ct = adata.obs[P["celltype_key"]].astype(str)
grp = ct.copy()
grp[adata.obs["epi_state"] == "TACSTD2_CLDN4_high"] = "Epi_TACSTD2high"
grp[adata.obs["epi_state"] == "TACSTD2_CLDN4_low"] = "Epi_TACSTD2low"
adata.obs["cc_group"] = pd.Categorical(grp)
print(adata.obs["cc_group"].value_counts())

# %% [markdown]
# ## 2. Run the consensus rank_aggregate
# `use_raw=False` -> use log-norm .X (NOT batch-corrected values). `expr_prop`
# drops LR pairs where <10% of a cluster expresses the gene: sparse chemokines
# otherwise produce noise. This is a WITHIN-sample analysis; for cross-patient
# claims see §5.

# %%
li.mt.rank_aggregate(
    adata,
    groupby="cc_group",
    expr_prop=P["min_expr_fraction"],
    use_raw=False,
    verbose=True,
    key_added="liana_res",
)
res = adata.uns["liana_res"].copy()
res.to_csv(OUTDIR / "liana_all_interactions.csv", index=False)
print(res.head())

# %% [markdown]
# ## 3. Filter to the hypothesis: epithelium -> T/NK on recruitment axes

# %%
recr = CFG["t_cell_recruitment"]
focus_ligands = sorted({g for ax in recr.values() for g in ax["ligands"]})
senders = ["Epi_TACSTD2high", "Epi_TACSTD2low"]
receivers = [c for c in adata.obs["cc_group"].cat.categories
             if ("T_cell" in c) or ("NK" in c) or ("Tcell" in c)]

def has_any(complex_str, genes):
    # ligand/receptor may be a complex "A_B"; match if any subunit is in genes.
    return any(s in genes for s in str(complex_str).split("_"))

mask = (
    res["source"].isin(senders)
    & res["target"].isin(receivers)
    & res["ligand_complex"].apply(lambda x: has_any(x, set(focus_ligands)))
)
focus = res[mask].sort_values("magnitude_rank")
focus.to_csv(OUTDIR / "liana_epi_to_T_recruitment.csv", index=False)
print(focus[["source", "target", "ligand_complex", "receptor_complex",
             "magnitude_rank", "specificity_rank"]].to_string(index=False))

# %% [markdown]
# ## 4. High vs low: is the recruitment signal weaker in TACSTD2/CLDN4-high?
# Compare the consensus magnitude_rank (lower = stronger) of the SAME LR pair
# from the high vs the low sender. A pair that is present/stronger from `low`
# but absent/weaker from `high` supports "reduced recruitment signal".
# This is descriptive; the statistic that supports a *claim* is §5.

# %%
key = ["ligand_complex", "receptor_complex", "target"]
piv = (
    focus.assign(sender=np.where(focus["source"] == "Epi_TACSTD2high", "high", "low"))
    .pivot_table(index=key, columns="sender", values="magnitude_rank", aggfunc="min")
    .reset_index()
)
piv["stronger_in"] = np.where(
    piv.get("high", np.nan).fillna(np.inf) < piv.get("low", np.nan).fillna(np.inf),
    "high", "low")
piv.to_csv(OUTDIR / "liana_high_vs_low.csv", index=False)
print(piv.to_string(index=False))

# %% [markdown]
# ## 5. Cross-patient version (the one you can actually defend)
# Run LIANA per sample, then compare the per-sample magnitude of focus pairs
# between conditions with a patient-level test. Cells are NOT replicates;
# patients are. Requires >= `min_patients_per_group` per group.

# %%
try:
    li.mt.rank_aggregate.by_sample(
        adata,
        sample_key=P["sample_key"],
        groupby="cc_group",
        expr_prop=P["min_expr_fraction"],
        use_raw=False,
        verbose=True,
        key_added="liana_by_sample",
    )
    bs = adata.uns["liana_by_sample"].copy()
    bs = bs[
        bs["source"].isin(senders)
        & bs["target"].isin(receivers)
        & bs["ligand_complex"].apply(lambda x: has_any(x, set(focus_ligands)))
    ]
    # Attach condition per sample and test high vs low WITHIN each condition,
    # or Epi_high across conditions -- with a nonparametric patient-level test.
    smeta = (adata.obs[[P["sample_key"], P["condition_key"]]]
             .drop_duplicates().set_index(P["sample_key"]))
    bs = bs.merge(smeta, left_on=P["sample_key"], right_index=True, how="left")
    bs.to_csv(OUTDIR / "liana_by_sample_focus.csv", index=False)

    from scipy.stats import mannwhitneyu
    rows = []
    for (lig, rec, tgt), sub in bs.groupby(["ligand_complex", "receptor_complex", "target"]):
        hi = sub.loc[sub.source == "Epi_TACSTD2high", "magnitude_rank"].dropna()
        lo = sub.loc[sub.source == "Epi_TACSTD2low", "magnitude_rank"].dropna()
        if len(hi) >= P["min_patients_per_group"] and len(lo) >= P["min_patients_per_group"]:
            u, p = mannwhitneyu(hi, lo, alternative="two-sided")
            rows.append((lig, rec, tgt, len(hi), len(lo),
                         float(hi.median()), float(lo.median()), p))
    stat = pd.DataFrame(rows, columns=["ligand", "receptor", "target",
                                       "n_hi", "n_lo", "median_rank_hi",
                                       "median_rank_lo", "pval"])
    if len(stat):
        from statsmodels.stats.multitest import multipletests
        stat["padj"] = multipletests(stat["pval"], method="fdr_bh")[1]
    stat.to_csv(OUTDIR / "liana_high_vs_low_patientlevel.csv", index=False)
    print(stat.to_string(index=False))
except Exception as e:
    print(f"[note] per-sample step skipped: {e}\n"
          "Need >=3 patients/group and liana>=1.0 with rank_aggregate.by_sample.")

# %% [markdown]
# ## 6. Visualize (spotlight only the focus interactions)

# %%
try:
    dp = li.pl.dotplot(
        adata,
        uns_key="liana_res",
        colour="magnitude_rank",
        size="specificity_rank",
        source_labels=senders,
        target_labels=receivers,
        top_n=25,
        orderby="magnitude_rank",
        orderby_ascending=True,   # lower magnitude_rank = stronger
        figure_size=(8, 7),
        return_fig=True,
    )
    dp.save(str(OUTDIR / "liana_dotplot.png"), dpi=150, verbose=False)
except Exception as e:
    print(f"[note] plotting skipped: {e}")

print(f"[done] LIANA outputs in {OUTDIR}")
