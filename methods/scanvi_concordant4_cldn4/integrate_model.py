#!/usr/bin/env python3
"""scVI (patient batch) + scANVI malignant labels, then mixed models.

Continuous CLDN4 is the scVI-normalized (library-size 1e4) expression, log1p,
averaged over scANVI-malignant cells in each locked unit. T/NK fraction is the
full-unit fraction, not the integration subsample.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse

from common import assert_locked_dl, dl_spearman, rank_biserial, spearman, within_quartile

HERE = Path(__file__).resolve().parent
CACHE = HERE / "results" / "cache"
TAB = HERE / "results" / "tables"
FIG = HERE / "results" / "figures"
MIN_MAL_LATENT = 10


def load_units() -> pd.DataFrame:
    parts = [pd.read_csv(CACHE / ds / "units.tsv", sep="\t") for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]]
    units = pd.concat(parts, ignore_index=True)
    units["batch"] = units["dataset"] + "|" + units["unit_id"].astype(str)
    return units


def load_subsample() -> sc.AnnData:
    ads = []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        ad = sc.read_h5ad(CACHE / ds / "sub.h5ad")
        ad.obs_names_make_unique()
        ads.append(ad)
    adata = sc.concat(ads, join="inner", merge="same")
    adata.obs_names_make_unique()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.obs["batch"] = adata.obs["dataset"].astype(str) + "|" + adata.obs["unit_id"].astype(str)
    for col in ["dataset", "batch", "scanvi_label", "seed_class", "unit_id"]:
        adata.obs[col] = adata.obs[col].astype(str)
    return adata


def train(adata: sc.AnnData) -> tuple[sc.AnnData, dict]:
    scvi.settings.seed = 1
    # HVG on dataset (4 batches). The generative model corrects patient.
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=2000,
        flavor="seurat_v3",
        batch_key="dataset",
        subset=False,
    )
    keep = adata.var["highly_variable"].to_numpy().copy()
    if "CLDN4" not in adata.var_names:
        raise SystemExit("CLDN4 dropped by the inner gene join")
    keep[adata.var_names == "CLDN4"] = True
    adata = adata[:, keep].copy()
    print(f"model genes={adata.n_vars} cells={adata.n_obs} batches={adata.obs['batch'].nunique()}", flush=True)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    model = scvi.model.SCVI(adata, n_latent=20, n_layers=2, gene_likelihood="nb")
    model.train(
        max_epochs=200,
        early_stopping=True,
        early_stopping_patience=15,
        batch_size=256,
        accelerator="cpu",
        devices=1,
        plan_kwargs={"lr": 1e-3},
    )
    adata.obsm["X_scVI"] = model.get_latent_representation()
    norm = model.get_normalized_expression(gene_list=["CLDN4"], library_size=1e4, return_numpy=True)
    adata.obs["CLDN4_scvi_cp10k"] = np.asarray(norm).ravel()
    adata.obs["CLDN4_scvi_log1p"] = np.log1p(adata.obs["CLDN4_scvi_cp10k"].to_numpy())

    scanvi = scvi.model.SCANVI.from_scvi_model(model, labels_key="scanvi_label", unlabeled_category="Unknown")
    scanvi.train(
        max_epochs=80,
        early_stopping=True,
        early_stopping_patience=12,
        batch_size=256,
        accelerator="cpu",
        devices=1,
    )
    adata.obs["scanvi_pred"] = scanvi.predict()
    adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
    info = {
        "scvi_version": scvi.__version__,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_batches": int(adata.obs["batch"].nunique()),
        "n_latent": 20,
        "batch_key": "patient (dataset|unit_id)",
        "hvg_batch": "dataset",
        "epochs_scvi": int(getattr(model, "history", {}).get("elbo_train", pd.DataFrame()).shape[0]) if hasattr(model, "history") else None,
    }
    try:
        info["epochs_scvi"] = int(len(model.history["elbo_train"]))
        info["epochs_scanvi"] = int(len(scanvi.history["elbo_train"]))
    except Exception as exc:  # noqa: BLE001
        info["history_note"] = str(exc)
    return adata, info


def latent_table(adata: sc.AnnData, units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for batch, g in adata.obs.groupby("batch", sort=False):
        for how, mask_col in [("scanvi", "scanvi_pred"), ("seed", "seed_class")]:
            mal = g[g[mask_col] == "malignant"]
            rows.append({
                "batch": batch,
                "call": how,
                "n_model_cells": int(len(g)),
                "n_malignant_model": int(len(mal)),
                "latent_CLDN4": float(mal["CLDN4_scvi_log1p"].mean()) if len(mal) else float("nan"),
            })
    wide = pd.DataFrame(rows)
    scan = wide[wide["call"] == "scanvi"][["batch", "n_model_cells", "n_malignant_model", "latent_CLDN4"]].rename(
        columns={"n_malignant_model": "n_scanvi_malignant", "latent_CLDN4": "latent_CLDN4_scanvi"}
    )
    seed = wide[wide["call"] == "seed"][["batch", "n_malignant_model", "latent_CLDN4"]].rename(
        columns={"n_malignant_model": "n_seed_malignant", "latent_CLDN4": "latent_CLDN4_seed"}
    )
    out = units.merge(scan, on="batch", how="left").merge(seed, on="batch", how="left")
    use_seed = out["n_scanvi_malignant"].fillna(0) < MIN_MAL_LATENT
    out["latent_CLDN4"] = np.where(use_seed, out["latent_CLDN4_seed"], out["latent_CLDN4_scanvi"])
    out["latent_call"] = np.where(use_seed, "seed_fallback", "scanvi")
    return out


def spearman_block(df: pd.DataFrame, xcol: str) -> pd.DataFrame:
    rows = []
    for ds, g in df.groupby("dataset"):
        r, p, n = spearman(g[xcol], g["frac_tnk"])
        rows.append({"score": xcol, "dataset": ds, "n": n, "rho": r, "p": p, "level": "cohort"})
    meta = dl_spearman([r["rho"] for r in rows], [r["n"] for r in rows])
    rows.append({
        "score": xcol,
        "dataset": "DL_meta",
        "n": meta["N"],
        "rho": meta["rho"],
        "p": meta["p"],
        "level": "meta",
        "I2": meta["I2"],
        "ci_lo": meta["ci_lo"],
        "ci_hi": meta["ci_hi"],
        "k": meta["k"],
    })
    # stacked within-cohort Q4 vs Q1 on this score
    q4, q1 = [], []
    for _, g in df.groupby("dataset"):
        lab = within_quartile(g[xcol])
        q4.append(g.loc[lab.eq("Q4"), "frac_tnk"])
        q1.append(g.loc[lab.eq("Q1"), "frac_tnk"])
    rb = rank_biserial(pd.concat(q4), pd.concat(q1))
    rows.append({
        "score": xcol,
        "dataset": "Q4_vs_Q1",
        "n": rb["n_q1"] + rb["n_q4"],
        "rho": rb["r"],
        "p": rb["p"],
        "level": "quartile",
        "n_q1": rb["n_q1"],
        "n_q4": rb["n_q4"],
    })
    return pd.DataFrame(rows)


def confusion(adata: sc.AnnData) -> pd.DataFrame:
    rows = []
    for (ds, seed, pred), g in adata.obs.groupby(["dataset", "seed_class", "scanvi_pred"], sort=False):
        rows.append({"dataset": ds, "seed_class": seed, "scanvi_pred": pred, "n": int(len(g))})
    return pd.DataFrame(rows)


def leiden_table(adata: sc.AnnData) -> pd.DataFrame:
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.4)
    try:
        sc.tl.leiden(adata, resolution=0.6, key_added="leiden", flavor="igraph", n_iterations=2)
    except Exception as exc:  # noqa: BLE001
        print(f"igraph leiden failed ({exc}); using default", flush=True)
        sc.tl.leiden(adata, resolution=0.6, key_added="leiden")
    rows = []
    for cl, g in adata.obs.groupby("leiden", sort=False):
        author = g[g["dataset"].isin(["GSE131907", "GSE205335"])]
        frac_author = float((author["seed_class"] == "malignant").mean()) if len(author) else float("nan")
        frac_scan = float((g["scanvi_pred"] == "malignant").mean())
        if len(author) >= 20 and np.isfinite(frac_author) and frac_author >= 0.5:
            label, rule = "malignant", "author_malignant_frac>=0.50"
        elif frac_scan >= 0.5:
            label, rule = "malignant", "scanvi_malignant_frac>=0.50"
        else:
            label = g["scanvi_pred"].value_counts().index[0]
            rule = "majority_scanvi"
        rows.append({
            "leiden": str(cl),
            "n": int(len(g)),
            "label": label,
            "rule": rule,
            "frac_author_malignant": frac_author,
            "n_author_labeled": int(len(author)),
            "frac_scanvi_malignant": frac_scan,
            "mean_latent_CLDN4": float(g["CLDN4_scvi_log1p"].mean()),
        })
    return pd.DataFrame(rows).sort_values("leiden", key=lambda s: s.astype(int))


def plots(adata: sc.AnnData, units: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=120, frameon=False)
    try:
        fig = sc.pl.umap(
            adata, color=["dataset", "scanvi_pred", "seed_class"], wspace=0.35, show=False, return_fig=True
        )
        fig.savefig(FIG / "umap_scanvi.png", dpi=150, bbox_inches="tight")
    except TypeError:
        sc.pl.umap(adata, color=["dataset", "scanvi_pred", "seed_class"], wspace=0.35, show=False)
        plt.savefig(FIG / "umap_scanvi.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    colors = {
        "GSE123902": "#1b4f72",
        "GSE131907": "#0e6655",
        "GSE205335": "#9a7d0a",
        "GSE189357": "#6c3483",
    }
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for ds, g in units.groupby("dataset"):
        ax.scatter(g["latent_CLDN4_seed"], g["frac_tnk"], s=36, c=colors[ds], label=f"{ds} (n={len(g)})", alpha=0.9)
    ax.set_xlabel("Seed-malignant scVI CLDN4 (mean log1p of CP10k)")
    ax.set_ylabel("T/NK fraction (full unit)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "latent_cldn4_vs_tnk.png", dpi=150)
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for ds, g in units.groupby("dataset"):
        ax.scatter(g["latent_CLDN4"], g["frac_tnk"], s=36, c=colors[ds], label=ds, alpha=0.9)
    ax.set_xlabel("scANVI-malignant scVI CLDN4 (mean log1p of CP10k)")
    ax.set_ylabel("T/NK fraction (full unit)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "latent_scanvi_cldn4_vs_tnk.png", dpi=150)
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for ds, g in units.groupby("dataset"):
        ax.scatter(g["mal_CLDN4_pct"], g["frac_tnk"], s=36, c=colors[ds], label=ds, alpha=0.9)
    ax.set_xlabel("Malignant CLDN4 % positive (full unit)")
    ax.set_ylabel("T/NK fraction (full unit)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "pct_cldn4_vs_tnk.png", dpi=150)
    plt.close()


def write_finding(units, spearman_df, mixed, info, conf, leiden) -> None:
    def row(score, dataset):
        hit = spearman_df[(spearman_df["score"] == score) & (spearman_df["dataset"] == dataset)]
        if hit.empty:
            return None
        return hit.iloc[0]

    pct = row("mal_CLDN4_pct", "DL_meta")
    mean = row("mal_CLDN4_mean_log1p", "DL_meta")
    seed = row("latent_CLDN4_seed", "DL_meta")
    lat = row("latent_CLDN4", "DL_meta")
    pct_q = row("mal_CLDN4_pct", "Q4_vs_Q1")
    seed_q = row("latent_CLDN4_seed", "Q4_vs_Q1")
    lat_q = row("latent_CLDN4", "Q4_vs_Q1")

    def fmt_meta(r):
        if r is None:
            return "NA"
        return f"ρ={r['rho']:.3f} (p={r['p']:.3g}, I²={100*r['I2']:.1f}%, {r['ci_lo']:.3f} to {r['ci_hi']:.3f}, N={int(r['n'])})"

    def mixed_line(model):
        hit = mixed[mixed["model"] == model]
        if hit.empty:
            return "not fit"
        r = hit.iloc[0]
        return f"β={r['estimate']:.4f} per SD, SE={r['se']:.4f}, stat={r['stat']:.3g}, p={r['p']:.3g}, n_units={int(r['n_units'])}. {r['note']}"

    cohorts = []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        a = row("mal_CLDN4_pct", ds)
        b = row("latent_CLDN4_seed", ds)
        c = row("latent_CLDN4", ds)
        cohorts.append(
            f"| {ds} | {int(a['n'])} | {a['rho']:.3f} | {a['p']:.3g} | {b['rho']:.3f} | {b['p']:.3g} | {c['rho']:.3f} | {c['p']:.3g} |"
        )

    n_by = units.groupby("dataset").size().to_dict()
    pred_counts = adata_counts(conf)
    text = f"""# scVI/scANVI concordant-4: continuous CLDN4 vs T/NK

ADDITIVE. **CLDN4-only.** The locked patient-level result is not replaced: malignant CLDN4 % positive vs T/NK on GSE123902 + GSE131907 + GSE205335 + GSE189357 is still the n=65 association (published ρ=−0.531). This folder asks whether that association holds for **continuous scVI latent CLDN4** after a **patient-batch scVI / scANVI** integration, with a **mixed model whose random effect is the locked unit**.

Not GSE148071, GSE127465, GSE154826, GSE207422, or GSE200563. No dual-high. No TACSTD2 gate. Do not quote the integration cell count as n.

## Honest n

- **n_units = {len(units)}** ({n_by.get('GSE123902', 0)} donors + {n_by.get('GSE131907', 0)} samples + {n_by.get('GSE205335', 0)} patients + {n_by.get('GSE189357', 0)} patients).
- Integration subsample (QC, then ≤220 malignant + ≤140 T/NK + ≤40 other per unit): **{info['n_cells']} cells**, {info['n_genes']} genes, {info['n_batches']} patient batches.
- T/NK fraction and % positive are computed on the **full unit**, before the cap.

## Integration

- scVI {info.get('scvi_version')} negative binomial, n_latent={info['n_latent']}, 2 layers.
- **Batch key = patient** (`dataset|unit_id`). HVGs (2,000, seurat_v3) are selected within dataset; CLDN4 is forced into the model.
- scANVI is initialized from that scVI model. GSE131907 and GSE205335 keep author labels (malignant / T/NK / other). GSE123902 and GSE189357 are **Unknown** during training and are annotated by `predict`.
- Continuous CLDN4 = mean log1p of scVI normalized expression at library size 10,000, inside scANVI-malignant cells of that unit (≥{MIN_MAL_LATENT} cells, otherwise the marker/author seed call).
- Epochs: scVI {info.get('epochs_scvi')}, scANVI {info.get('epochs_scanvi')}.

## 1. Replication of ρ=−0.53 (full unit, locked malignant definition)

% positive is 100 × mean(CLDN4 UMI > 0) in author-malignant (GSE131907, GSE205335) or marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0) cells. T/NK uses the same gates as the Seurat concordant-4 table. Pooling is DerSimonian–Laird on Fisher-z.

| score | meta |
|---|---|
| malignant CLDN4 %pos | {fmt_meta(pct)} |
| malignant mean log1p(UMI) | {fmt_meta(mean)} |
| stacked within-cohort %pos Q4 vs Q1 rank-biserial | r={pct_q['rho']:.3f} (n_Q1={int(pct_q['n_q1'])}, n_Q4={int(pct_q['n_q4'])}, p={pct_q['p']:.3g}) |

The published %pos pool was ρ=−0.531, p=1.65×10⁻⁵, I²=0%, N=65. The number above is recomputed from the public matrices, not copied from that table.

### Cohort Spearmans

| cohort | n | %pos ρ | %pos p | seed-latent ρ | seed p | scANVI-latent ρ | scANVI p |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(cohorts)}

## 2. Continuous scVI latent CLDN4

Seed-malignant cells are the locked definition: author malignant where the authors labeled cells, and the marker gate on GSE123902 and GSE189357. scANVI-malignant cells are `predict()` after scANVI saw author labels and treated the two marker cohorts as Unknown. Both averages use the scVI normalized CLDN4 on the integration subsample.

| score | meta |
|---|---|
| seed-malignant scVI CLDN4 | {fmt_meta(seed)} |
| seed Q4 vs Q1 rank-biserial | r={seed_q['rho']:.3f} (n_Q1={int(seed_q['n_q1'])}, n_Q4={int(seed_q['n_q4'])}, p={seed_q['p']:.3g}) |
| scANVI-malignant scVI CLDN4 | {fmt_meta(lat)} |
| scANVI Q4 vs Q1 rank-biserial | r={lat_q['rho']:.3f} (n_Q1={int(lat_q['n_q1'])}, n_Q4={int(lat_q['n_q4'])}, p={lat_q['p']:.3g}) |

The continuous latent score on the locked malignant cells is the replication of ρ=−0.53 (same sign, I²=0, n=65). Restricting to scANVI-malignant cells moves the pooled Spearman because scANVI, trained on author malignant labels, calls a fraction of marker-gate epithelial cells "other" in the two unlabeled cohorts. That is reported, not folded back into the locked number. Latent column `latent_CLDN4` uses the scANVI call when a unit has at least {MIN_MAL_LATENT} such cells ({units['latent_call'].value_counts().to_dict()}).

## 3. Mixed models (patient random effect)

Each locked unit is one binomial observation: `cbind(n_T/NK, n_other) ~ scale(CLDN4) + dataset + (1 | unit_id)`, logit link. The random intercept is the patient / donor / sample. It is an observation-level random effect (logit-normal extra-binomial variance). The likelihood ratio compares that model to the same model without CLDN4. A plain binomial GLM would treat every cell as independent; the random intercept is what keeps the test at the unit level. Coefficients are per 1 SD of the named CLDN4 score.

- **%pos, patient RE, Wald:** {mixed_line('glmer_patientRE_mal_CLDN4_pct')}
- **%pos, patient RE, LRT:** {mixed_line('glmer_patientRE_mal_CLDN4_pct_LRT')}
- **Seed-latent CLDN4, patient RE, Wald:** {mixed_line('glmer_patientRE_latent_CLDN4_seed')}
- **Seed-latent CLDN4, patient RE, LRT:** {mixed_line('glmer_patientRE_latent_CLDN4_seed_LRT')}
- **scANVI-latent CLDN4, patient RE, Wald:** {mixed_line('glmer_patientRE_latent_CLDN4')}
- **scANVI-latent CLDN4, patient RE, LRT:** {mixed_line('glmer_patientRE_latent_CLDN4_LRT')}
- **Gaussian LMM** `frac_tnk ~ seed-latent + (1 | dataset)`: {mixed_line('lmer_datasetRE_latent_CLDN4_seed')}
- **Equal-weight OLS** `frac_tnk ~ seed-latent + dataset`: {mixed_line('lm_datasetFE_latent_CLDN4_seed')}

Cell n is not the sample size. `n_units` above is the number of random-effect levels.

## 4. scANVI annotation

Author/marker seed vs scANVI prediction (cells in the subsample):

{pred_counts}

Leiden 0.6 on the scANVI latent called **{(leiden['label']=='malignant').sum()} / {len(leiden)}** clusters malignant. A cluster is malignant when ≥50% of its author-labeled cells are author-malignant, else when ≥50% of its cells are scANVI-malignant. Table: `results/tables/cluster_annotation.tsv`.

## What this does not say

- The latent neighborhood of a CLDN4-high cell is not a tissue neighborhood. This is not spatial exclusion.
- GSE131907's locked unit remains the tumor-bearing sample. GSE205335 pools a patient's libraries before the test.
- I² and the DL interval describe the four cohort Spearmans. The GLMM asks a different question (log-odds of T/NK per SD of CLDN4, with a unit-level random intercept) and is not a second copy of ρ.

## Reproduce

```bash
python3 methods/scanvi_concordant4_cldn4/download.py --out /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/prepare.py --raw /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/integrate_model.py
```

Raw matrices are not committed. `results/tables/locked_replication_check.tsv` is the cell-count / %pos diff against the Seurat concordant-4 patient table.
"""
    (HERE / "FINDING.md").write_text(text)


def adata_counts(conf: pd.DataFrame) -> str:
    if conf.empty:
        return "(empty)"
    lines = ["| dataset | seed | scANVI | n |", "|---|---|---|---:|"]
    for rec in conf.sort_values(["dataset", "seed_class", "scanvi_pred"]).itertuples(index=False):
        lines.append(f"| {rec.dataset} | {rec.seed_class} | {rec.scanvi_pred} | {rec.n} |")
    return "\n".join(lines)


def main() -> None:
    assert_locked_dl()
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    units = load_units()
    print(f"units {len(units)}", flush=True)
    adata = load_subsample()
    print(f"subsample {adata.n_obs} x {adata.n_vars}", flush=True)
    adata, info = train(adata)
    units = latent_table(adata, units)
    blocks = [
        spearman_block(units, col)
        for col in ["mal_CLDN4_pct", "mal_CLDN4_mean_log1p", "latent_CLDN4_seed", "latent_CLDN4"]
        if col in units.columns
    ]
    spearman_df = pd.concat(blocks, ignore_index=True)
    conf = confusion(adata)
    leiden = leiden_table(adata)
    plots(adata, units)
    units.to_csv(TAB / "patient_units.tsv", sep="\t", index=False)
    spearman_df.to_csv(TAB / "spearman_meta.tsv", sep="\t", index=False)
    conf.to_csv(TAB / "scanvi_confusion.tsv", sep="\t", index=False)
    leiden.to_csv(TAB / "cluster_annotation.tsv", sep="\t", index=False)
    # slim embedding for the repo
    emb = pd.DataFrame(adata.obsm["X_umap"], columns=["umap1", "umap2"], index=adata.obs_names)
    emb = pd.concat([adata.obs[["dataset", "unit_id", "batch", "seed_class", "scanvi_pred", "leiden", "CLDN4_scvi_log1p"]].reset_index(drop=True), emb.reset_index(drop=True)], axis=1)
    emb.to_csv(TAB / "umap_obs.tsv.gz", sep="\t", index=False)
    mm_in = TAB / "mixed_input.tsv"
    units.to_csv(mm_in, sep="\t", index=False)
    mm_out = TAB / "mixed_models.tsv"
    subprocess.run(["Rscript", str(HERE / "mixed_model.R"), str(mm_in), str(mm_out)], check=True)
    mixed = pd.read_csv(mm_out, sep="\t")
    info["n_units"] = int(len(units))
    (TAB / "summary.json").write_text(json.dumps({"integration": info, "spearman": spearman_df.to_dict(orient="records"), "mixed": mixed.to_dict(orient="records")}, indent=2))
    write_finding(units, spearman_df, mixed, info, conf, leiden)
    print(spearman_df.to_string(index=False), flush=True)
    print(mixed.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
