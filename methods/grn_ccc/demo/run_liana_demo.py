# %% [markdown]
# # DEMO - LIANA end-to-end on a tiny dataset (mechanics + readout check)
#
# **Data policy.** This repo does not ship or download a large dataset. Per the
# playbook, a real demo is run "only if a small scRNA <2GB is available":
#   * If you set `DEMO_H5AD` to a small `.h5ad` with columns `cell_type`,
#     `epi_state`, `patient_id`, the demo uses it.
#   * Otherwise it SYNTHESIZES a small, labeled AnnData (no network, <2 MB) whose
#     ground truth is "TACSTD2/CLDN4-high epithelium expresses fewer T-cell
#     recruitment chemokines". This validates the LIANA pipeline AND shows the
#     expected readout. It is a mechanics/logic check, not biological evidence.
#
# Run: `python run_liana_demo.py`  (needs scanpy + liana; see ../env/requirements.txt)

# %%
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import liana as li

OUTDIR = Path(__file__).resolve().parent / "output"
OUTDIR.mkdir(exist_ok=True)
RNG = np.random.default_rng(0)
DEMO_H5AD = os.environ.get("DEMO_H5AD", "")

# Genes: recruitment chemokines + receptors + lineage markers.
CHEMOKINES = ["CXCL9", "CXCL10", "CXCL11", "CCL5", "CXCL16"]
RECEPTORS = ["CXCR3", "CCR5", "CXCR6"]
EPI = ["EPCAM", "KRT8", "KRT18", "TACSTD2", "CLDN4"]
TCELL = ["CD3D", "CD3E", "CD8A", "IL7R"]
MYELOID = ["LYZ", "CD68"]
GENES = EPI + TCELL + MYELOID + CHEMOKINES + RECEPTORS

# LIANA refuses to run unless >=2% of its resource genes are present. For the
# synthetic dataset we therefore draw background gene names from LIANA's own
# consensus resource (real symbols) so the coverage check passes.
_resource = li.rs.select_resource("consensus")
_res_genes = pd.unique(_resource[["ligand", "receptor"]].values.ravel())
NOISE_GENES = [g for g in _res_genes if isinstance(g, str)
               and "_" not in g and g not in GENES][:400]
ALLGENES = GENES + NOISE_GENES


def synth():
    """Small labeled AnnData with a planted 'reduced recruitment in high' signal."""
    specs = [  # (cell_type, epi_state, n_cells)
        ("Epithelial", "TACSTD2_CLDN4_high", 400),
        ("Epithelial", "TACSTD2_CLDN4_low", 400),
        ("T_cell", "non_epithelial", 500),
        ("Myeloid", "non_epithelial", 300),
    ]
    blocks, obs = [], []
    n_pat = 6
    for ct, state, n in specs:
        base = {g: 0.15 for g in ALLGENES}          # low ambient baseline
        if ct == "Epithelial":
            for g in EPI:
                base[g] = 6.0
            if state == "TACSTD2_CLDN4_high":
                for g in ["TACSTD2", "CLDN4"]:
                    base[g] = 10.0
                # planted effect: high state has ~5x LOWER chemokine ligands
                for g in CHEMOKINES:
                    base[g] = 0.6
            else:
                for g in CHEMOKINES:
                    base[g] = 3.0
        if ct == "T_cell":
            for g in TCELL:
                base[g] = 6.0
            for g in RECEPTORS:
                base[g] = 4.0
        if ct == "Myeloid":
            for g in MYELOID:
                base[g] = 6.0
            base["CCL5"] = 2.0                        # myeloid also makes some CCL5
        lam = np.array([base[g] for g in ALLGENES])
        blocks.append(RNG.poisson(lam=lam, size=(n, len(ALLGENES))))
        pats = RNG.integers(0, n_pat, size=n)
        for k in range(n):
            obs.append({"cell_type": ct, "epi_state": state,
                        "patient_id": f"P{pats[k]:02d}",
                        "ici_response": "R" if pats[k] % 2 == 0 else "NR"})
    X = np.vstack(blocks).astype(np.float32)
    a = ad.AnnData(X=X, obs=pd.DataFrame(obs), var=pd.DataFrame(index=ALLGENES))
    a.obs_names = [f"cell{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = a.obs["cell_type"].astype("category")
    return a


# %% ---- load or synthesize ---------------------------------------------------
if DEMO_H5AD and Path(DEMO_H5AD).exists():
    adata = sc.read_h5ad(DEMO_H5AD)
    print(f"[demo] using provided dataset {DEMO_H5AD}: {adata.shape}")
    source = "provided"
else:
    adata = synth()
    print(f"[demo] synthesized labeled AnnData: {adata.shape} (mechanics check)")
    source = "synthetic"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %% ---- build sender/receiver groups (Epi split high/low) --------------------
ct = adata.obs["cell_type"].astype(str)
grp = ct.copy()
grp[adata.obs["epi_state"] == "TACSTD2_CLDN4_high"] = "Epi_TACSTD2high"
grp[adata.obs["epi_state"] == "TACSTD2_CLDN4_low"] = "Epi_TACSTD2low"
adata.obs["cc_group"] = pd.Categorical(grp)
print(adata.obs["cc_group"].value_counts())

# %% ---- run LIANA consensus --------------------------------------------------
li.mt.rank_aggregate(adata, groupby="cc_group", expr_prop=0.1,
                     use_raw=False, verbose=False, key_added="liana_res")
res = adata.uns["liana_res"]
res.to_csv(OUTDIR / "demo_liana_all.csv", index=False)

# %% ---- extract the hypothesis: epithelium -> T cell on recruitment axes -----
focus_lig = set(CHEMOKINES)
senders = ["Epi_TACSTD2high", "Epi_TACSTD2low"]


def has_any(x, genes):
    return any(s in genes for s in str(x).split("_"))


m = (res["source"].isin(senders) & (res["target"] == "T_cell")
     & res["ligand_complex"].apply(lambda x: has_any(x, focus_lig)))
focus = res[m].copy()
focus.to_csv(OUTDIR / "demo_liana_epi_to_T.csv", index=False)

piv = (focus.assign(sender=np.where(focus.source == "Epi_TACSTD2high", "high", "low"))
       .pivot_table(index=["ligand_complex", "receptor_complex"],
                    columns="sender", values="magnitude_rank", aggfunc="min"))
print("\n=== Epithelium -> T cell recruitment interactions "
      "(magnitude_rank; LOWER = stronger) ===")
print(piv.to_string())

# %% ---- readout --------------------------------------------------------------
if {"high", "low"}.issubset(piv.columns):
    weaker = (piv["high"].fillna(np.inf) > piv["low"].fillna(np.inf)).mean()
    print(f"\nFraction of recruitment pairs weaker (worse rank) from the HIGH "
          f"state: {weaker:.0%}")
    verdict = ("as planted: recruitment signal is weaker from TACSTD2/CLDN4-high"
               if weaker >= 0.5 else "not weaker in this run")
    print(f"[demo readout | {source}] {verdict}.")
print(f"[done] demo outputs in {OUTDIR}")
