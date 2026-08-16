# %% [markdown]
# # 06 - Load a public lung ICI scRNA object into the playbook schema
#
# Maps a GEO processed matrix onto the columns the rest of the pipeline
# expects (`cell_type`, `patient_id`, `ici_response`, then `epi_state` via
# the same score as `00_preprocess.py`). Default = **GSE207422** (175 MB,
# <2 GB demo cap). See `../config/public_datasets.yaml`.
#
# This is a *loader*, not a re-analysis of the paper. Authors already QC'd
# and removed doublets; ambient RNA cannot be re-corrected from this matrix.
# Record that in `.uns['ambient_correction']`.

# %%
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
PUB = yaml.safe_load((HERE.parent / "config" / "public_datasets.yaml").read_text())
P = CFG["params"]

# ---- EDIT -------------------------------------------------------------------
ACCESSION = "GSE207422"
CACHE = Path("data/public") / ACCESSION
OUT_H5AD = Path("data/lung_ici_preprocessed.h5ad")
# -----------------------------------------------------------------------------

ds = PUB["datasets"][ACCESSION]


def download(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[cache] {dest}")
        return dest
    import urllib.request
    print(f"[get] {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


# %% [markdown]
# ## GSE207422: gene x cell UMI matrix + xlsx metadata
# Matrix is genes (rows) x cells (columns), gzipped TSV. Metadata xlsx has
# one row per cell (or per sample — inspect and join). Column names vary
# slightly across GEO dumps; the aliases below cover the published table.

# %%
def load_gse207422(cache: Path):
    files = {f["name"]: download(f["url"], cache / f["name"]) for f in ds["processed"]}
    mtx = files["GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"]
    meta_p = files["GSE207422_NSCLC_scRNAseq_metadata.xlsx"]

    # pandas read_csv on a genes x cells TSV is RAM-heavy (~1–2 GB). Use
    # scanpy's reader when possible; fall back to pandas.
    print("[read] UMI matrix (this can take a few minutes)…")
    try:
        adata = sc.read_text(mtx).T          # cells x genes
    except Exception:
        df = pd.read_csv(mtx, sep="\t", index_col=0, compression="gzip")
        import anndata as ad
        adata = ad.AnnData(df.T)

    meta = pd.read_excel(meta_p)
    # GSE207422 xlsx is SAMPLE-level (one row per BD_immuneXX), not cell-level.
    # Cell barcodes look like BD_immune01_612637. Join on the sample prefix.
    if "Sample" in meta.columns and meta["Sample"].astype(str).str.startswith("BD_immune").any():
        meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
        meta = meta.set_index("Sample")
        sample = adata.obs_names.to_series().str.extract(r"(BD_immune\d+)", expand=False)
        adata.obs["sample"] = sample.values
        colmap = {"Patient": P["sample_key"],
                  "Pathologic Response": P["condition_key"],
                  "RECIST": "recist", "Pathology": "pathology"}
        for src, dest in colmap.items():
            if src in meta.columns:
                adata.obs[dest] = sample.map(meta[src]).values
        print("[join] sample-level metadata via BD_immune prefix")
    else:
        id_col = next((c for c in meta.columns
                       if str(c).lower() in {"cell", "cellid", "cell_id", "barcode",
                                             "cell.name", "cell_name"}),
                      meta.columns[0])
        meta = meta.set_index(id_col)
        meta.index = meta.index.astype(str)
        adata.obs_names = adata.obs_names.astype(str)
        shared = adata.obs_names.intersection(meta.index)
        if len(shared) < 0.5 * adata.n_obs:
            print(f"[warn] only {len(shared)}/{adata.n_obs} barcodes match metadata.")
            adata.obs = adata.obs.join(meta, how="left")
        else:
            adata = adata[shared].copy()
            adata.obs = meta.loc[shared]

    # Map published columns -> playbook schema (aliases, case-insensitive).
    # GSE207422 does NOT deposit a cell-type column; marker-gate in that case.
    aliases = {
        P["celltype_key"]: ["celltype", "cell_type", "cell.type", "major_celltype",
                            "celltype_major", "annotation", "cluster_celltype"],
        P["sample_key"]:   ["patient", "patient_id", "patient.id", "orig.ident",
                            "sample", "sample_id", "donor"],
        P["condition_key"]: ["pathologic_response", "path_response", "response",
                             "MPR", "mpr", "recist", "group"],
    }
    lower = {c.lower(): c for c in adata.obs.columns}

    def pick(cands):
        for a in cands:
            if a.lower() in lower:
                return lower[a.lower()]
        return None

    for dest, cands in aliases.items():
        src = pick(cands)
        if src and dest not in adata.obs:
            adata.obs[dest] = adata.obs[src].astype(str)
        elif dest not in adata.obs:
            print(f"[action needed] no column matching {cands}; set adata.obs['{dest}']")

    if P["celltype_key"] not in adata.obs:
        print("[note] no cell_type column deposited; marker-gating lineages "
              "(see demo/run_gse207422_demo.py for the same rule).")
        scores = {}
        Xlog = adata.X
        # expect already log-norm later; use raw means as a coarse gate
        for name, genes in CFG["lineage_markers"].items():
            g = [x for x in genes if x in adata.var_names]
            scores[name] = (np.asarray(adata[:, g].X.mean(axis=1)).ravel()
                            if g else np.zeros(adata.n_obs))
        S = pd.DataFrame(scores, index=adata.obs_names)
        # map lineage_markers keys -> playbook labels
        remap = {"epithelial": "Epithelial", "t_cell": "T_cell",
                 "nk_cell": "NK", "myeloid": "Myeloid", "b_cell": "B_cell",
                 "fibroblast": "Fibroblast", "endothelial": "Endothelial"}
        S = S.rename(columns=remap)
        adata.obs[P["celltype_key"]] = S.idxmax(axis=1)
        adata.obs.loc[S.max(axis=1) < 0.05, P["celltype_key"]] = "Unknown"

    adata.uns["ambient_correction"] = (
        "NOT re-run; GSE207422 authors QC + Scrublet (Hu et al. 2023). "
        "Ambient RNA cannot be re-estimated from this processed matrix."
    )
    adata.uns["public_dataset"] = {"accession": ACCESSION, **{k: ds[k] for k in
        ("title", "platform", "ici_setting", "n_patients_scrna") if k in ds}}
    return adata


# %% [markdown]
# ## Define epi_state on the loaded object (same rule as 00_preprocess.py)

# %%
def define_epi_state(ad):
    ct = ad.obs.get(P["celltype_key"], pd.Series(index=ad.obs_names, dtype=str))
    epi_mask = ct.astype(str).str.contains("pithel|tumor|malignan|cancer", case=False)
    if epi_mask.sum() < 50:
        print("[warn] few epithelial cells by name; scoring ALL cells then gating "
              "on EPCAM/KRT later. Check cell_type labels.")
        epi_mask = pd.Series(True, index=ad.obs_names)
    epi = ad[epi_mask].copy()
    pos = [g for g in CFG["state_markers"]["positive"] if g in epi.var_names]
    if len(pos) < 2:
        raise SystemExit(f"TACSTD2/CLDN4 not in var_names: {pos}")
    sc.tl.score_genes(epi, pos, score_name="tacstd2_cldn4_score",
                      random_state=P["random_seed"])
    s = epi.obs["tacstd2_cldn4_score"].to_numpy()
    thr = float(np.quantile(s, P["state_high_quantile"]))
    label = np.where(s >= thr, "TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low")
    ad.obs["epi_state"] = "non_epithelial"
    ad.obs.loc[epi.obs_names, "epi_state"] = label
    ad.obs["epi_state"] = ad.obs["epi_state"].astype("category")
    ad.uns["epi_state_def"] = {"method": "quantile", "threshold": thr,
                               "markers": pos, "accession": ACCESSION}
    return ad


# %%
if ACCESSION == "GSE207422":
    adata = load_gse207422(CACHE)
else:
    raise SystemExit(f"loader for {ACCESSION} not implemented; add a function "
                     "mirroring load_gse207422() and map columns to the schema.")

# Keep raw counts; log-norm for CCC tools.
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
if "log1p" not in adata.uns:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

adata = define_epi_state(adata)
print(adata)
print(adata.obs["epi_state"].value_counts())
print("ambient:", adata.uns.get("ambient_correction"))

OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUT_H5AD)
print(f"[done] {OUT_H5AD}")
