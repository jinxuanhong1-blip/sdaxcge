#!/usr/bin/env python3
"""Harmony joint object: epithelial TACSTD2/CLDN4 vs T/NK and B/TLS-like.

Unit is the tumor donor (Kim tumor site, Xiang patient, Wu patient,
Zilionis patient). This is extra atlas n, not an ICI test.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import harmonypy as hm

from config import (
    EXCLUDED,
    EXTRACTED,
    GENE_PANEL,
    HARMONY_MAX_PER_DONOR,
    HARMONY_NPC,
    LINEAGES,
    MIN_B,
    MIN_EPI,
    MIN_TNK,
    NORMAL_LUNG,
    RANDOM_SEED,
    RESULTS,
    TARGETS,
    TLS_CHEMOKINE,
)

DATASETS = ["GSE131907", "GSE253013", "GSE148071", "GSE127465"]

AUTHOR_EPI = {
    "epithelial cells",
    "epithelial",
    "malignant cells",
    "type i cells",
    "type ii cells",
    "ciliated cells",
    "club cells",
    "at1",
    "at2",
    "ciliated",
    "club",
}
AUTHOR_T = {"t lymphocytes", "t cells", "tt cells", "cd8 t", "cd4 t"}
AUTHOR_NK = {"nk cells", "tnk cells"}
AUTHOR_B = {"b lymphocytes", "b cells", "tb cells", "plasma cells", "tplasma cells"}


def log1p_mean(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return float("nan")
    return float(np.mean(np.log1p(np.clip(x, 0, None))))


def pct_pos(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return float("nan")
    return float(np.mean(x > 0) * 100.0)


def marker_score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(np.clip(expr[g], 0, None)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0).astype(np.float32)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = [k for k in LINEAGES if k != "plasma"]
    scores = np.vstack([marker_score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.12) | ((best_val - second) < 0.04)] = "other"
    return assigned


def map_author(label: str) -> str:
    s = str(label).strip().lower()
    if not s:
        return ""
    if "patient" in s and "specific" in s:
        return "epithelial"
    if s in AUTHOR_EPI or s.startswith("epithelial"):
        return "epithelial"
    if s in AUTHOR_T or s.startswith("t cell") or s.startswith("tt cell"):
        return "T"
    if s in AUTHOR_NK or "nk" == s or s.startswith("tnk"):
        return "NK"
    if s in AUTHOR_B or "plasma" in s or s.startswith("tb cell") or s.startswith("b cell"):
        return "B"
    if "myeloid" in s or "mac" in s or "mono" in s or "neutrophil" in s or "mast" in s:
        return "myeloid"
    if "fibro" in s or "smooth muscle" in s:
        return "fibroblast"
    if "endoth" in s:
        return "endothelial"
    return ""


def spearman_row(x, y, contrast: str, extra: dict | None = None) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    rec = {
        "contrast": contrast,
        "n": int(len(x)),
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "note": "",
    }
    if extra:
        rec.update(extra)
    if len(x) < 5:
        rec["note"] = "too_few_donors"
        return rec
    rho, p = stats.spearmanr(x, y)
    rec["spearman_rho"] = float(rho)
    rec["spearman_p"] = float(p)
    return rec


def residualize(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Subtract per-group mean rank so pooled Spearman is not a batch effect."""
    out = np.full(values.shape, np.nan, dtype=float)
    for g in pd.unique(groups):
        idx = groups == g
        v = values[idx]
        ok = np.isfinite(v)
        if ok.sum() < 3:
            continue
        ranks = stats.rankdata(v[ok])
        out_idx = np.where(idx)[0][ok]
        out[out_idx] = ranks - ranks.mean()
    return out


def load_dataset(name: str) -> tuple[dict[str, np.ndarray], pd.DataFrame] | None:
    d = EXTRACTED / name
    expr_p = d / "gene_panel.npz"
    meta_p = d / "cell_metadata.tsv"
    if not expr_p.exists() or not meta_p.exists():
        print(f"skip {name}: not extracted", flush=True)
        return None
    packed = np.load(expr_p)
    expr = {k: packed[k] for k in packed.files if k != "total"}
    meta = pd.read_csv(meta_p, sep="\t")
    if "dataset" not in meta.columns:
        meta["dataset"] = name
    if "donor" not in meta.columns:
        meta["donor"] = meta["patient"] if "patient" in meta.columns else "unknown"
    if "tissue" not in meta.columns:
        meta["tissue"] = "Tumor"
    n = len(meta)
    for g, v in list(expr.items()):
        if len(v) != n:
            raise SystemExit(f"{name} {g} length {len(v)} != {n}")
    info = {}
    gi = d / "gene_index.json"
    if gi.exists():
        info = json.loads(gi.read_text())
    print(f"loaded {name} n={n} TACSTD2={info.get('TACSTD2')} CLDN4={info.get('CLDN4')}", flush=True)
    return expr, meta


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(ok.sum())
    tmp[order] = q
    out[ok] = tmp
    return out.tolist()


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    loaded = []
    catalog = []
    for name in DATASETS:
        got = load_dataset(name)
        gi = EXTRACTED / name / "gene_index.json"
        rec = {"dataset": name, "included": got is not None}
        if gi.exists():
            rec.update(json.loads(gi.read_text()))
        if got is None:
            rec["reason"] = "extraction_missing"
            catalog.append(rec)
            continue
        expr, meta = got
        rec["n_cells_extracted"] = int(len(meta))
        rec["TACSTD2"] = "TACSTD2" in expr
        rec["CLDN4"] = "CLDN4" in expr
        if not rec["TACSTD2"] or not rec["CLDN4"]:
            rec["included"] = False
            rec["reason"] = "missing_TACSTD2_or_CLDN4"
            catalog.append(rec)
            continue
        loaded.append((name, expr, meta))
        catalog.append(rec)

    for acc, why in EXCLUDED.items():
        catalog.append({"dataset": acc, "included": False, "reason": why})

    if not loaded:
        raise SystemExit("no datasets extracted")

    # Build joint tables
    metas = []
    expr_joint: dict[str, list[np.ndarray]] = {g: [] for g in GENE_PANEL}
    for name, expr, meta in loaded:
        n = len(meta)
        tumor = meta["tissue"].astype(str).str.upper().isin({"TUMOR", "T", "TLUNG"})
        # GSE253013 uses Tumor / ANT
        if name == "GSE253013":
            tumor = meta["tissue"].astype(str).eq("Tumor")
        if name == "GSE131907":
            tumor = np.ones(n, dtype=bool)  # already tumor-site filtered
        meta = meta.loc[tumor].copy()
        expr = {g: v[tumor.to_numpy()] for g, v in expr.items()}
        n = len(meta)
        lineage = assign_lineage(expr, n)
        meta["lineage"] = lineage
        author = meta["author_cell_type"].map(map_author) if "author_cell_type" in meta.columns else ""
        meta["author_lineage"] = author if isinstance(author, pd.Series) else ""
        normal = marker_score(expr, [g for g in NORMAL_LUNG if g in expr], n)
        meta["malignant_like"] = (meta["lineage"] == "epithelial") & (normal <= 0.05)
        if name == "GSE131907" and "author_subtype" in meta.columns:
            malig_sub = meta["author_subtype"].isin(["Malignant cells", "tS1", "tS2", "tS3"])
            meta["author_malignant"] = malig_sub
        elif name == "GSE127465" and "author_cell_type" in meta.columns:
            meta["author_malignant"] = meta["author_cell_type"].astype(str).str.contains(
                "Patient", case=False, na=False
            ) & meta["author_cell_type"].astype(str).str.contains("specific", case=False, na=False)
        else:
            meta["author_malignant"] = meta["malignant_like"]
        meta["tnk"] = meta["lineage"].isin(["T", "NK"])
        meta["b_tls"] = meta["lineage"].isin(["B", "plasma"])
        # TLS-like: B/plasma OR high CXCL13/CCL19/CCL21 among non-epithelial
        tls = marker_score(expr, [g for g in TLS_CHEMOKINE if g in expr], n)
        meta["tls_chemokine"] = tls
        meta["b_tls_like"] = meta["b_tls"] | ((meta["lineage"] != "epithelial") & (tls >= 0.8))
        for g in TARGETS:
            meta[g] = expr[g] if g in expr else np.nan
        for g in GENE_PANEL:
            expr_joint[g].append(expr[g] if g in expr else np.zeros(n, dtype=np.float32))
        metas.append(meta.reset_index(drop=True))

    joint = pd.concat(metas, ignore_index=True)
    for g in list(expr_joint):
        expr_joint[g] = np.concatenate(expr_joint[g]) if expr_joint[g] else np.zeros(0)

    # Per-donor metrics
    rows = []
    for (dataset, donor), pdf in joint.groupby(["dataset", "donor"], sort=False):
        n_cells = len(pdf)
        n_epi = int((pdf["lineage"] == "epithelial").sum())
        n_mal = int(pdf["malignant_like"].sum())
        n_tnk = int(pdf["tnk"].sum())
        n_b = int(pdf["b_tls"].sum())
        n_btls = int(pdf["b_tls_like"].sum())
        rec = {
            "dataset": dataset,
            "donor": str(donor),
            "n_cells": n_cells,
            "n_epithelial": n_epi,
            "n_malignant_like": n_mal,
            "n_tnk": n_tnk,
            "n_b": n_b,
            "n_b_tls_like": n_btls,
            "tnk_fraction": n_tnk / n_cells if n_cells else np.nan,
            "b_fraction": n_b / n_cells if n_cells else np.nan,
            "b_tls_like_fraction": n_btls / n_cells if n_cells else np.nan,
            "eligible": n_epi >= MIN_EPI and n_tnk >= MIN_TNK,
            "eligible_b": n_epi >= MIN_EPI and n_b >= MIN_B,
        }
        for gene in TARGETS:
            rec[f"epi_{gene}_mean_log1p"] = log1p_mean(pdf.loc[pdf["lineage"] == "epithelial", gene])
            rec[f"epi_{gene}_pct_pos"] = pct_pos(pdf.loc[pdf["lineage"] == "epithelial", gene])
            rec[f"mal_{gene}_mean_log1p"] = log1p_mean(pdf.loc[pdf["malignant_like"], gene])
            rec[f"mal_{gene}_pct_pos"] = pct_pos(pdf.loc[pdf["malignant_like"], gene])
        rows.append(rec)
    per = pd.DataFrame(rows)
    per.to_csv(RESULTS / "per_donor_metrics.tsv", sep="\t", index=False)

    tests = []
    elig = per[per["eligible"]].copy()
    elig_b = per[per["eligible_b"]].copy()

    primary_specs = [
        ("epi_TACSTD2_mean_log1p", "tnk_fraction", "epithelial TACSTD2 mean log1p vs T/NK fraction", elig),
        ("epi_CLDN4_mean_log1p", "tnk_fraction", "epithelial CLDN4 mean log1p vs T/NK fraction", elig),
        ("mal_TACSTD2_mean_log1p", "tnk_fraction", "malignant-like TACSTD2 mean log1p vs T/NK fraction", elig),
        ("mal_CLDN4_mean_log1p", "tnk_fraction", "malignant-like CLDN4 mean log1p vs T/NK fraction", elig),
        ("epi_TACSTD2_mean_log1p", "b_tls_like_fraction", "epithelial TACSTD2 mean log1p vs B/TLS-like fraction", elig_b),
        ("epi_CLDN4_mean_log1p", "b_tls_like_fraction", "epithelial CLDN4 mean log1p vs B/TLS-like fraction", elig_b),
        ("mal_TACSTD2_mean_log1p", "b_tls_like_fraction", "malignant-like TACSTD2 mean log1p vs B/TLS-like fraction", elig_b),
        ("mal_CLDN4_mean_log1p", "b_tls_like_fraction", "malignant-like CLDN4 mean log1p vs B/TLS-like fraction", elig_b),
        ("epi_TACSTD2_pct_pos", "tnk_fraction", "epithelial TACSTD2 %pos vs T/NK fraction", elig),
        ("epi_CLDN4_pct_pos", "tnk_fraction", "epithelial CLDN4 %pos vs T/NK fraction", elig),
    ]

    for xcol, ycol, label, df in primary_specs:
        tests.append(spearman_row(df[xcol], df[ycol], label + " (joint, all eligible donors)"))
        # dataset-residualized ranks
        xr = residualize(df[xcol].to_numpy(), df["dataset"].to_numpy())
        yr = residualize(df[ycol].to_numpy(), df["dataset"].to_numpy())
        tests.append(spearman_row(xr, yr, label + " (dataset-residualized ranks)"))
        for ds, sub in df.groupby("dataset"):
            tests.append(spearman_row(sub[xcol], sub[ycol], label + f" [{ds}]", {"dataset": ds}))

    # Kim tLung-only sensitivity
    kim = elig[elig["dataset"] == "GSE131907"]
    tests.append(
        spearman_row(
            kim["epi_TACSTD2_mean_log1p"],
            kim["tnk_fraction"],
            "epithelial TACSTD2 vs T/NK [GSE131907 tumor sites only]",
        )
    )

    tests_df = pd.DataFrame(tests)
    primary_mask = tests_df["contrast"].str.contains(r"\(joint, all eligible", regex=True)
    tests_df.loc[primary_mask, "bh_fdr"] = bh_fdr(tests_df.loc[primary_mask, "spearman_p"].tolist())
    tests_df.to_csv(RESULTS / "association_statistics.tsv", sep="\t", index=False)

    # Harmony on a balanced subsample
    rng = np.random.default_rng(RANDOM_SEED)
    take = []
    for (dataset, donor), idx in joint.groupby(["dataset", "donor"]).groups.items():
        idx = np.asarray(list(idx))
        if len(idx) > HARMONY_MAX_PER_DONOR:
            idx = rng.choice(idx, size=HARMONY_MAX_PER_DONOR, replace=False)
        take.append(idx)
    take = np.concatenate(take)
    take.sort()
    genes_ok = [g for g in GENE_PANEL if g in expr_joint and float(np.max(expr_joint[g])) > 0]
    X = np.vstack([np.log1p(np.clip(expr_joint[g][take], 0, None)) for g in genes_ok]).T
    # z-score within dataset so Harmony is not just protocol scale
    sub = joint.iloc[take].reset_index(drop=True)
    Xs = np.zeros_like(X)
    for ds, idx in sub.groupby("dataset").groups.items():
        idx = np.asarray(list(idx))
        block = X[idx]
        scaler = StandardScaler()
        Xs[idx] = scaler.fit_transform(block)
    n_pc = min(HARMONY_NPC, Xs.shape[1], Xs.shape[0] - 1)
    pcs = PCA(n_components=n_pc, random_state=RANDOM_SEED).fit_transform(Xs)
    ho = hm.run_harmony(pcs, sub, vars_use=["dataset"], max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr).T
    sub["harmony_1"] = Z[:, 0]
    sub["harmony_2"] = Z[:, 1]
    sub.to_csv(RESULTS / "harmony_subsample.tsv", sep="\t", index=False)

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0), constrained_layout=True)
    pairs = [
        ("epi_TACSTD2_mean_log1p", "tnk_fraction", "Epithelial TACSTD2 vs T/NK"),
        ("epi_CLDN4_mean_log1p", "tnk_fraction", "Epithelial CLDN4 vs T/NK"),
        ("epi_TACSTD2_mean_log1p", "b_tls_like_fraction", "Epithelial TACSTD2 vs B/TLS-like"),
        ("epi_CLDN4_mean_log1p", "b_tls_like_fraction", "Epithelial CLDN4 vs B/TLS-like"),
    ]
    colors = {
        "GSE131907": "#1f4e79",
        "GSE253013": "#b85c38",
        "GSE148071": "#2e7d4f",
        "GSE127465": "#6b4c9a",
    }
    for ax, (yc, xc, title) in zip(axes.ravel(), pairs):
        df = elig if "tnk" in xc else elig_b
        for ds, subdf in df.groupby("dataset"):
            ax.scatter(
                subdf[xc],
                subdf[yc],
                s=36,
                c=colors.get(ds, "gray"),
                label=ds,
                alpha=0.85,
                edgecolors="none",
            )
        blk = spearman_row(df[yc], df[xc], title)
        ax.set_xlabel(xc.replace("_", " "))
        ax.set_ylabel(yc.replace("_", " "))
        ax.set_title(f"{title}\nn={blk['n']}  ρ={blk['spearman_rho']:.2f}  p={blk['spearman_p']:.2g}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Joint naive LUAD/NSCLC tumor scRNA (Harmony object) · not an ICI test", fontsize=11)
    fig.savefig(RESULTS / "fig_joint_tacstd2_cldn4_vs_immune.png", dpi=180, bbox_inches="tight")
    fig.savefig(RESULTS / "fig_joint_tacstd2_cldn4_vs_immune.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    for ax, color_by, cmap in (
        (axes[0], "dataset", None),
        (axes[1], "lineage", None),
    ):
        cats = sub[color_by].astype(str)
        uniq = sorted(cats.unique())
        pal = plt.cm.tab10(np.linspace(0, 1, max(len(uniq), 1)))
        for i, u in enumerate(uniq):
            m = cats == u
            ax.scatter(sub.loc[m, "harmony_1"], sub.loc[m, "harmony_2"], s=4, alpha=0.45, label=u, c=[pal[i]])
        ax.set_xlabel("Harmony 1")
        ax.set_ylabel("Harmony 2")
        ax.set_title(color_by)
        ax.legend(markerscale=3, fontsize=7, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.savefig(RESULTS / "fig_harmony_joint.png", dpi=160)
    fig.savefig(RESULTS / "fig_harmony_joint.pdf")
    plt.close(fig)

    primary = tests_df[primary_mask].copy()
    summary = {
        "scope": "extra_atlas_n_not_ici",
        "unit": "tumor_donor",
        "harmony": {
            "batch": "dataset",
            "n_cells_subsampled": int(len(sub)),
            "n_genes": len(genes_ok),
            "n_pcs": n_pc,
            "max_per_donor": HARMONY_MAX_PER_DONOR,
        },
        "n_joint_tumor_cells": int(len(joint)),
        "n_donors_all": int(len(per)),
        "n_eligible_epi_tnk": int(elig.shape[0]),
        "n_eligible_epi_b": int(elig_b.shape[0]),
        "n_by_dataset": per.groupby("dataset").size().to_dict(),
        "n_eligible_by_dataset": elig.groupby("dataset").size().to_dict(),
        "primary": primary.to_dict(orient="records"),
        "excluded": EXCLUDED,
        "catalog": catalog,
        "eligibility": {"min_epithelial": MIN_EPI, "min_tnk": MIN_TNK, "min_b": MIN_B},
        "note": (
            "Scores are joint-object marker lineages (argmax of mean log1p "
            "lineage genes). Malignant-like = epithelial with near-zero normal "
            "lung markers. B/TLS-like = B/plasma or non-epithelial cells with "
            "high CXCL13/CCL19/CCL21. Not histology. Not ICI."
        ),
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame(catalog).to_csv(RESULTS / "dataset_catalog.tsv", sep="\t", index=False)
    print(json.dumps({k: summary[k] for k in ("n_eligible_epi_tnk", "n_eligible_by_dataset", "primary")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
