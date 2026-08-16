# %% [markdown]
# # DEMO — GSE207422 (public NSCLC ICI scRNA, 175 MB, <2 GB)
#
# Hu et al., *Genome Med* 2023 (PMID 36869384). Neoadjuvant PD-1 + chemo,
# 15 patients, ~92k cells, BD Rhapsody (not 10x). Authors' processed UMI
# matrix + sample metadata. Ambient/doublets already handled upstream.
#
# **What this demo claims / does not claim.** It runs the playbook's
# *measurable* tests on a real public object:
#   1. marker-gate epithelium vs T/NK (GEO did not deposit cell-type labels);
#   2. split epithelium into TACSTD2/CLDN4-high vs -low (top quartile);
#   3. patient-level pseudobulk of recruitment chemokines (PRIMARY);
#   4. LIANA consensus on the focus axes (hypothesis only).
# It does **not** claim "reduced T-cell recruitment" unless the patient-level
# test supports it. n=15, MPR n=4 — underpowered for response splits.
#
# Run: `python run_gse207422_demo.py`  (downloads ~176 MB on first run)

# %%
from __future__ import annotations
import gzip
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import liana as li
import yaml

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
P = CFG["params"]
OUT = HERE / "gse207422_example"
OUT.mkdir(exist_ok=True)
CACHE = Path("/tmp/gse207422")
CACHE.mkdir(exist_ok=True)

URL_MTX = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
           "suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
URL_META = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/"
            "suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx")

CHEMOKINES = sorted({g for ax in CFG["t_cell_recruitment"].values()
                     for g in ax["ligands"]})
RECEPTORS = sorted({r for ax in CFG["t_cell_recruitment"].values()
                    for r in ax["receptor"]})
FOCUS_REC = {"CXCR3", "CCR5", "CXCR6", "CCR1"}


def fetch(url, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}")
    urllib.request.urlretrieve(url, dest)
    return dest


# %% [markdown]
# ## 1. Stream only the genes we need (keeps RAM well under 2 GB)
# Full dense 24k × 92k would be ~9 GB. We keep lineage + state + recruitment
# genes plus LIANA consensus symbols so the resource-coverage check passes.

# %%
def wanted_genes() -> set[str]:
    w = set()
    for block in CFG["state_markers"].values():
        w.update(block)
    for block in CFG["lineage_markers"].values():
        w.update(block)
    w.update(CHEMOKINES)
    w.update(RECEPTORS)
    for ax in CFG["immunosuppressive_recruitment"].values():
        w.update(ax["ligands"]); w.update(ax["receptor"])
    res = li.rs.select_resource("consensus")
    w.update(pd.unique(res[["ligand", "receptor"]].values.ravel()))
    return {g for g in w if isinstance(g, str) and "_" not in g}


def stream_matrix(path: Path, keep: set[str]) -> ad.AnnData:
    print(f"[read] streaming {path.name} for {len(keep)} symbols…")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        genes, rows = [], []
        for line in f:
            gene, rest = line.split("\t", 1)
            if gene not in keep:
                continue
            genes.append(gene)
            rows.append(np.fromstring(rest, sep="\t", dtype=np.float32))
    X = np.vstack(rows).T                         # cells x genes
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=cells),
                   var=pd.DataFrame(index=genes))
    a.obs_names = cells
    a.var_names = genes
    print(f"[matrix] {a.n_obs} cells x {a.n_vars} genes (streamed)")
    return a


# %% [markdown]
# ## 2. Attach sample metadata (xlsx is sample-level, not cell-level)

# %%
def attach_sample_meta(a: ad.AnnData, xlsx: Path) -> ad.AnnData:
    meta = pd.read_excel(xlsx)
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.set_index("Sample")
    sample = a.obs_names.to_series().str.extract(r"(BD_immune\d+)", expand=False)
    a.obs["sample"] = sample.values
    a.obs["patient_id"] = sample.map(meta["Patient"]).values
    a.obs["ici_response"] = sample.map(meta["Pathologic Response"]).values
    a.obs["recist"] = sample.map(meta["RECIST"]).values
    a.obs["pathology"] = sample.map(meta["Pathology"]).values
    a.obs["resource"] = sample.map(meta["Resource"]).values
    print(a.obs["patient_id"].value_counts().sort_index().to_string())
    return a


# %% [markdown]
# ## 3. Marker-gate lineages (GEO deposited no cell-type column)

# %%
def _mean_present(a, genes):
    g = [x for x in genes if x in a.var_names]
    if not g:
        return np.zeros(a.n_obs)
    X = a[:, g].X
    return np.asarray(X.mean(axis=1)).ravel()


def marker_gate(a: ad.AnnData) -> ad.AnnData:
    # Use log-norm for scores, keep counts in .X until then.
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    scores = {
        "Epithelial": _mean_present(a, CFG["lineage_markers"]["epithelial"]),
        "T_cell":     _mean_present(a, CFG["lineage_markers"]["t_cell"]),
        "NK":         _mean_present(a, CFG["lineage_markers"]["nk_cell"]),
        "Myeloid":    _mean_present(a, CFG["lineage_markers"]["myeloid"]),
        "B_cell":     _mean_present(a, CFG["lineage_markers"]["b_cell"]),
        "Fibroblast": _mean_present(a, CFG["lineage_markers"]["fibroblast"]),
        "Endothelial": _mean_present(a, CFG["lineage_markers"]["endothelial"]),
    }
    S = pd.DataFrame(scores, index=a.obs_names)
    a.obs["cell_type"] = S.idxmax(axis=1)
    # Require the winning score to beat a floor so ambient crumbs aren't typed.
    a.obs.loc[S.max(axis=1) < 0.4, "cell_type"] = "Unknown"
    print(a.obs["cell_type"].value_counts().to_string())
    return a


def define_epi_state(a: ad.AnnData) -> ad.AnnData:
    epi = a[a.obs["cell_type"] == "Epithelial"].copy()
    pos = [g for g in CFG["state_markers"]["positive"] if g in epi.var_names]
    sc.tl.score_genes(epi, pos, score_name="tacstd2_cldn4_score",
                      random_state=P["random_seed"])
    thr = float(np.quantile(epi.obs["tacstd2_cldn4_score"],
                            P["state_high_quantile"]))
    lab = np.where(epi.obs["tacstd2_cldn4_score"] >= thr,
                   "TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low")
    a.obs["epi_state"] = "non_epithelial"
    a.obs.loc[epi.obs_names, "epi_state"] = lab
    a.obs["epi_state"] = a.obs["epi_state"].astype("category")
    a.uns["epi_state_def"] = {"method": "quantile", "threshold": thr,
                              "markers": pos, "accession": "GSE207422"}
    print(a.obs["epi_state"].value_counts().to_string())
    print("epi_state_def", a.uns["epi_state_def"])
    return a


# %% [markdown]
# ## 4. PRIMARY: patient-level chemokine pseudobulk (epithelium high vs low)

# %%
def pseudobulk_guardrail(a: ad.AnnData) -> pd.DataFrame:
    epi = a[a.obs["cell_type"] == "Epithelial"].copy()
    Xc = epi.layers["counts"]
    ck = [g for g in CHEMOKINES if g in epi.var_names]
    rows = []
    for (pat, state), idx in epi.obs.groupby(["patient_id", "epi_state"],
                                             observed=True).indices.items():
        if state == "non_epithelial" or len(idx) < 20:
            continue
        sub = np.asarray(Xc[idx])
        lib = sub.sum()
        if lib <= 0:
            continue
        cpm = sub.sum(axis=0) / lib * 1e6
        rec = {"patient": pat, "state": state, "n_cells": len(idx)}
        for g in ck:
            rec[g] = float(np.log1p(cpm[epi.var_names.get_loc(g)]))
        rows.append(rec)
    tbl = pd.DataFrame(rows)
    from scipy.stats import wilcoxon
    from statsmodels.stats.multitest import multipletests
    out = []
    for g in ck:
        wide = tbl.pivot_table(index="patient", columns="state", values=g)
        both = wide.dropna(subset=["TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low"])
        if len(both) < P["min_patients_per_group"]:
            out.append((g, len(both), np.nan, np.nan, np.nan, np.nan))
            continue
        hi, lo = both["TACSTD2_CLDN4_high"], both["TACSTD2_CLDN4_low"]
        try:
            _, p = wilcoxon(hi, lo)
        except ValueError:
            p = np.nan
        out.append((g, len(both), float(hi.mean()), float(lo.mean()),
                    float((hi - lo).mean()), p))
    res = pd.DataFrame(out, columns=["chemokine", "n_paired", "mean_logcpm_high",
                                     "mean_logcpm_low", "delta_high_minus_low",
                                     "pval"])
    if res["pval"].notna().any():
        mask = res["pval"].notna()
        res.loc[mask, "padj"] = multipletests(res.loc[mask, "pval"],
                                              method="fdr_bh")[1]
    res = res.sort_values("delta_high_minus_low")
    res.to_csv(OUT / "gse207422_chemokine_pseudobulk.csv", index=False)
    tbl.to_csv(OUT / "gse207422_pseudobulk_logcpm.csv", index=False)
    print("\n=== PRIMARY pseudobulk (neg. delta = lower in HIGH state) ===")
    print(res.to_string(index=False))
    return res


# %% [markdown]
# ## 5. LIANA on focus axes (hypothesis generator, not the claim)

# %%
def run_liana(a: ad.AnnData) -> pd.DataFrame:
    keep_types = {"Epithelial", "T_cell", "NK"}
    sub = a[a.obs["cell_type"].isin(keep_types)].copy()
    grp = sub.obs["cell_type"].astype(str)
    grp = grp.mask(sub.obs["epi_state"] == "TACSTD2_CLDN4_high", "Epi_TACSTD2high")
    grp = grp.mask(sub.obs["epi_state"] == "TACSTD2_CLDN4_low", "Epi_TACSTD2low")
    sub.obs["cc_group"] = pd.Categorical(grp)
    print(sub.obs["cc_group"].value_counts().to_string())
    li.mt.rank_aggregate(sub, groupby="cc_group", expr_prop=0.10,
                         use_raw=False, verbose=False, key_added="liana_res")
    res = sub.uns["liana_res"]
    senders = ["Epi_TACSTD2high", "Epi_TACSTD2low"]
    rec_ok = res["receptor_complex"].apply(
        lambda x: any(r in str(x).split("_") for r in FOCUS_REC))
    lig_ok = res["ligand_complex"].apply(
        lambda x: any(c in str(x).split("_") for c in CHEMOKINES))
    focus = res[res["source"].isin(senders)
                & res["target"].isin(["T_cell", "NK"])
                & lig_ok & rec_ok].copy()
    focus.to_csv(OUT / "gse207422_liana_focus.csv", index=False)
    piv = (focus.assign(sender=np.where(focus.source == "Epi_TACSTD2high",
                                        "high", "low"))
           .pivot_table(index=["ligand_complex", "receptor_complex", "target"],
                        columns="sender", values="magnitude_rank",
                        aggfunc="min"))
    piv.to_csv(OUT / "gse207422_liana_high_vs_low.csv")
    print("\n=== LIANA focus (magnitude_rank; LOWER = stronger) ===")
    print(piv.to_string() if len(piv) else "[none passed expr_prop / filter]")
    return focus


# %%
def main():
    mtx = fetch(URL_MTX, CACHE / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
    xlsx = fetch(URL_META, CACHE / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    keep = wanted_genes()
    a = stream_matrix(mtx, keep)
    a = attach_sample_meta(a, xlsx)
    a = marker_gate(a)
    a = define_epi_state(a)
    a.uns["ambient_correction"] = (
        "NOT re-run; GSE207422 authors QC + Scrublet (Hu et al. 2023). "
        "Ambient cannot be re-estimated from this processed matrix."
    )
    a.uns["public_dataset"] = {
        "accession": "GSE207422",
        "platform": "BD Rhapsody WTA (not 10x)",
        "n_patients": 15,
        "note": "marker-gated lineages; GEO deposited no cell-type column",
    }
    pb = pseudobulk_guardrail(a)
    try:
        run_liana(a)
    except Exception as e:
        print(f"[note] LIANA skipped: {e}")

    n_lower = int((pb["delta_high_minus_low"] < 0).sum()) if len(pb) else 0
    n_sig = int(((pb.get("padj", pd.Series(dtype=float)) < 0.05)
                 & (pb["delta_high_minus_low"] < 0)).sum()) if "padj" in pb else 0
    print("\n=== DEMO READOUT (GSE207422) ===")
    print(f"chemokines lower in TACSTD2/CLDN4-high: {n_lower}/{len(pb)}; "
          f"FDR<0.05: {n_sig}")
    print("This is a patient-level description of THIS cohort, not a "
          "general claim about T-cell recruitment.")
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
