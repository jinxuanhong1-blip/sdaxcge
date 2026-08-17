#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA vs malignant CLDN4 on GSE131907.

Additive. GSE207422 Milo is a different agent and is not run here.
Author cell labels are public (Kim et al. 2020). miloR / edgeR are not used.

Primary graph = tLung (same tissue). Secondary same-site graph = mBrain.
Independent unit = sample. Do not cite cell count as n.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from knn_nhood import (  # noqa: E402
    EPS,
    attach_fdr,
    count_matrix,
    greedy_independent,
    make_nhoods,
    nhood_composition,
    nhood_membership,
    pca_knn,
    refine_indices,
    sample_paired_tnk_by_gene,
    sample_sizes,
    select_hvg,
    self_test,
    spearman_da,
    spearman_safe,
    welch_da,
    wilcoxon_paired,
    write_json,
)

MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
IMMUNE_TYPES = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
}
TNK_TYPES = {"T lymphocytes", "NK cells"}
KEEP_TYPES = IMMUNE_TYPES | {"Epithelial cells"}
MARKERS = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "NKG7",
    "GNLY",
    "CD79A",
    "MS4A1",
    "LYZ",
    "CD68",
    "PECAM1",
    "COL1A1",
]
COHORTS = ("tLung", "mBrain")


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    geo = pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})
    geo = geo.rename(columns={"title": "Sample"})
    return geo


def origin_from_sample(name: str) -> str:
    if name.startswith("LUNG_N"):
        return "nLung"
    if name.startswith("LUNG_T"):
        return "tLung"
    if name.startswith("LN_"):
        return "nLN"
    if name.startswith("EFFUSION_"):
        return "PE"
    if name.startswith("NS_"):
        return "mBrain"
    if name.startswith(("EBUS_", "BRONCHO_")):
        return "unknown_lung_or_ln"
    return "unknown"


def load_annotation(ann_path: Path, series_path: Path, cell_ids: list[str]) -> pd.DataFrame:
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index: {list(ann.columns)}")
    ann = ann.set_index("Index").reindex(cell_ids).reset_index()
    if ann["Sample"].isna().any():
        n_miss = int(ann["Sample"].isna().sum())
        raise SystemExit(f"{n_miss} matrix cell IDs missing from annotation")
    geo = parse_series_matrix(series_path)
    geo_small = geo[["Sample"]].copy()
    for col in ("patient_id", "tumor_stage", "source_name_ch1", "lung_cancer_subtype"):
        if col in geo.columns:
            geo_small[col] = geo[col]
    # GEO origin is often in source_name or a characteristics field
    if "sample_origin" in geo.columns:
        geo_small["Sample_Origin_geo"] = geo["sample_origin"]
    elif "source_name_ch1" in geo.columns:
        geo_small["Sample_Origin_geo"] = geo["source_name_ch1"]
    ann = ann.merge(geo_small, on="Sample", how="left")
    if "Sample_Origin" in ann.columns:
        origin = ann["Sample_Origin"].fillna("")
    else:
        origin = pd.Series("", index=ann.index)
    mapped = ann["Sample"].map(origin_from_sample)
    # Prefer author Sample_Origin when it is a known site code
    known = {"tLung", "nLung", "tL/B", "mLN", "nLN", "PE", "mBrain"}
    use_author = origin.isin(known)
    ann["origin"] = np.where(use_author, origin, mapped)
    # EBUS/BRONCHO without author origin: leave as unknown; author file usually has it
    return ann


def stream_pass1(
    matrix_path: Path,
    markers: set[str],
    cohort_masks: dict[str, np.ndarray],
):
    t0 = time.time()
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[0] not in {"Index", "Gene", "gene"}:
            # Kim matrix starts with Index
            pass
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        names: list[str] = []
        stats = {c: {"mean": [], "var": []} for c in cohort_masks}
        found: dict[str, np.ndarray] = {}
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} values, expected {n}")
            n_umi += arr
            names.append(gene)
            for c, mask in cohort_masks.items():
                sub = arr[mask]
                stats[c]["mean"].append(float(sub.mean()) if sub.size else 0.0)
                stats[c]["var"].append(float(sub.var()) if sub.size else 0.0)
            if gene in markers:
                found[gene] = arr
            if i % 2000 == 0:
                print(
                    f"  pass1 {i} genes, markers={len(found)} [{time.time() - t0:.0f}s]",
                    flush=True,
                )
    out_stats = {
        c: (np.asarray(stats[c]["mean"]), np.asarray(stats[c]["var"])) for c in stats
    }
    return cell_ids, n_umi, names, out_stats, found


def stream_pass2(matrix_path: Path, keep: set[str], keep_idx: np.ndarray, n_cells: int) -> dict[str, np.ndarray]:
    t0 = time.time()
    out: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene not in keep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            out[gene] = arr[keep_idx].astype(np.float32, copy=False)
            if len(out) % 200 == 0:
                print(f"  pass2 stored {len(out)}/{len(keep)} [{time.time() - t0:.0f}s]", flush=True)
    missing = keep - set(out)
    if missing:
        print(f"  missing {sorted(missing)[:12]} ({len(missing)} genes)", flush=True)
    return out


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]]
    rec = {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((t["p"] < 0.05).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()) if len(t) else 0,
        "min_p": float(t["p"].min()) if len(t) else None,
        "min_SpatialFDR": float(t["SpatialFDR"].min()) if len(t) else None,
        "min_BH_FDR": float(t["BH_FDR"].min()) if len(t) else None,
    }
    return rec


def fmt(x, nd=3, sci=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if sci:
        return f"{x:.2e}"
    return f"{x:.{nd}f}"


def plot_volcano(df: pd.DataFrame, x: str, title: str, path: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    testable = df["testable"].to_numpy()
    ax.scatter(
        df.loc[~testable, x],
        -np.log10(np.clip(df.loc[~testable, "p"], EPS, 1)),
        s=6,
        c="#bbbbbb",
        alpha=0.4,
        label="not testable",
    )
    t = df[testable]
    sig = t["SpatialFDR"] < 0.1
    ax.scatter(
        t.loc[~sig, x],
        -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)),
        s=8,
        c="#4c72b0",
        alpha=0.55,
        label="testable SpatialFDR≥0.1",
    )
    ax.scatter(
        t.loc[sig, x],
        -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)),
        s=14,
        c="#c44e52",
        alpha=0.85,
        label="SpatialFDR<0.1",
    )
    ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"−log10 p (sample-level test)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(x, y, s=8, c="#4c72b0", alpha=0.45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run_cohort(
    name: str,
    local: np.ndarray,
    sample_keep: np.ndarray,
    lineage_keep: np.ndarray,
    is_malig_keep: np.ndarray,
    is_tnk_keep: np.ndarray,
    cldn_keep: np.ndarray,
    n_umi_keep: np.ndarray,
    expr: dict[str, np.ndarray],
    hvg_names: list[str],
    outdir: Path,
    k: int,
    d: int,
    prop: float,
    n_hvg_requested: int,
    min_interface: int,
    seed: int,
    min_malig: int,
) -> dict:
    sample = sample_keep[local]
    lineage = lineage_keep[local]
    is_malig = is_malig_keep[local]
    is_tnk = is_tnk_keep[local]
    cldn = cldn_keep[local]
    scale = np.where(n_umi_keep[local] > 0, 1e4 / n_umi_keep[local], 0.0)
    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g][local] * scale) for g in hvg_present]).T.astype(np.float32)
    print(f"{name}: cells={local.sum()} HVG={X.shape[1]} samples={pd.unique(sample).size}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample, n_pcs=d, k=k, random_state=seed)
    indices = refine_indices(pcs, knn_idx, prop=prop, random_state=seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods)
    print(f"{name}: nhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)

    mal_score = []
    tnk_frac = []
    n_mal_s = []
    n_tnk_s = []
    for s in sample_levels:
        m = sample == s
        n_mal = int((m & is_malig).sum())
        n_mal_s.append(n_mal)
        n_tnk_s.append(int((m & is_tnk).sum()))
        mal_score.append(float(cldn[m & is_malig].mean()) if n_mal >= min_malig else np.nan)
        tnk_frac.append(float(is_tnk[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac = np.asarray(tnk_frac, dtype=float)

    da_cldn = attach_fdr(spearman_da(counts, sizes, mal_score, min_samples=5), nhoods.k_distance)
    scored = np.isfinite(mal_score)
    group = np.full(len(sample_levels), "unscored", dtype=object)
    if scored.sum() >= 4:
        med = float(np.median(mal_score[scored]))
        group[scored & (mal_score > med)] = "high"
        group[scored & (mal_score < med)] = "low"
        # ties at the median stay unscored
    n_high = int((group == "high").sum())
    n_low = int((group == "low").sum())
    if n_high >= 2 and n_low >= 2:
        da_split = attach_fdr(
            welch_da(counts, sizes, group, group_a="low", group_b="high", min_per_group=2, min_samples=4),
            nhoods.k_distance,
        )
    else:
        da_split = pd.DataFrame(
            {
                "nhood": np.arange(len(members)),
                "p": np.nan,
                "testable": False,
                "logFC_B_minus_A": np.nan,
                "BH_FDR": np.nan,
                "SpatialFDR": np.nan,
            }
        )

    comp = nhood_composition(members, is_tnk, is_malig, cldn, lineage, gene_prefix="cldn4")
    merged = comp.merge(
        da_cldn.rename(
            columns={
                "spearman_rho": "cldn4_rho",
                "p": "cldn4_p",
                "testable": "cldn4_testable",
                "BH_FDR": "cldn4_BH_FDR",
                "SpatialFDR": "cldn4_SpatialFDR",
                "n_cells": "n_cells_da",
            }
        ),
        on="nhood",
        how="left",
    )
    dest = outdir / name
    dest.mkdir(parents=True, exist_ok=True)
    merged.to_csv(dest / "nhoods.tsv", sep="\t", index=False)
    da_cldn.to_csv(dest / "da_malignant_cldn4.tsv", sep="\t", index=False)
    da_split.to_csv(dest / "da_cldn4_median_split.tsv", sep="\t", index=False)

    iface = (comp["n_malig"] >= min_interface) & (comp["n_tnk"] >= min_interface)
    all_sp = spearman_safe(comp["cldn4_malig_mean"], comp["frac_tnk"])
    iface_sp = spearman_safe(comp.loc[iface, "cldn4_malig_mean"], comp.loc[iface, "frac_tnk"])
    indep = greedy_independent(members, max_shared=0)
    indep_mask = np.zeros(len(members), dtype=bool)
    indep_mask[indep] = True
    indep_sp = spearman_safe(
        comp.loc[indep_mask & iface.to_numpy(), "cldn4_malig_mean"],
        comp.loc[indep_mask & iface.to_numpy(), "frac_tnk"],
    )
    has_mal = comp["n_malig"] >= 5
    med_n = float(np.nanmedian(comp.loc[has_mal, "cldn4_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() >= med_n)
    low = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() < med_n)
    paired = sample_paired_tnk_by_gene(members, sample, is_tnk, high, low)
    sample_tab = pd.DataFrame(
        {
            "sample": sample_levels,
            "origin": name,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "n_tnk": n_tnk_s,
            "malignant_cldn4": mal_score,
            "frac_tnk": tnk_frac,
            "cldn4_arm": group,
        }
    )
    paired = paired.merge(sample_tab, on="sample")
    paired.to_csv(dest / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    sample_tab.to_csv(dest / "sample_scores.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    sample_level = spearman_safe(mal_score, tnk_frac)

    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "stat"}),
        "stat",
        f"GSE131907 {name}: nhood abundance vs malignant CLDN4",
        dest / "fig_da_cldn4_volcano.png",
        "Spearman ρ (nhood proportion vs sample malignant CLDN4)",
    )
    if da_split["testable"].any():
        plot_volcano(
            da_split,
            "logFC_B_minus_A",
            f"GSE131907 {name}: nhood DA CLDN4-high vs low",
            dest / "fig_da_cldn4_split_volcano.png",
            "log2 FC (high − low sample proportion)",
        )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        dest / "fig_interface_cldn4_vs_tnk.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"{name} interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                ax.plot(
                    [0, 1],
                    [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]],
                    "-o",
                    c="#4c72b0",
                    alpha=0.75,
                    ms=5,
                )
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title(f"GSE131907 {name} sample-paired (unit = sample)")
        fig.tight_layout()
        fig.savefig(dest / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)

    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = [f"{r.sample}={int(r.n_malignant)}" for r in dropped.itertuples()]
    summary = {
        "cohort": name,
        "n_cells_in_graph": int(local.sum()),
        "n_malignant": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "lineage_counts": {k: int(v) for k, v in pd.Series(lineage).value_counts().items()},
        "graph": {
            "k": k,
            "d": min(d, X.shape[1], int(local.sum()) - 1),
            "prop": prop,
            "n_hvg": len(hvg_present),
            "n_hvg_requested": n_hvg_requested,
            "n_nhoods": int(len(members)),
            "nhood_size": k + 1,
            "median_nhood_size": float(np.median([m.size for m in members])),
            "batch": "PCA per-sample mean centering (not Harmony)",
        },
        "samples": {
            "n": int(len(sample_levels)),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
            "n_cldn4_high": n_high,
            "n_cldn4_low": n_low,
            "dropped_lt10_malignant": drop_bits,
            "median_split_cut": float(np.median(mal_score[scored])) if scored.any() else None,
        },
        "sample_level_malignant_CLDN4_vs_TNK_fraction": sample_level,
        "da_malignant_cldn4": fdr_counts(da_cldn),
        "da_cldn4_median_split": fdr_counts(da_split),
        "composition": {
            "note": "transcriptional kNN != spatial niche",
            "all_nhoods_maligCLDN4_vs_fracTNK": all_sp,
            "interface_nhoods": {
                "definition": f">={min_interface} malignant AND >={min_interface} T/NK cells",
                "n": int(iface.sum()),
                "spearman": iface_sp,
            },
            "disjoint_interface": {
                "n_disjoint_nhoods": int(indep.size),
                "n_disjoint_interface": int((indep_mask & iface.to_numpy()).sum()),
                "spearman": indep_sp,
            },
            "sample_paired_TNK_in_CLDN4high_vs_low_nhoods": {
                "high_low_cut": "median malignant CLDN4 among nhoods with >=5 malignant cells",
                "median_cut": med_n,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
            },
        },
        "honest_n": {
            "independent_unit_DA": f"{name} samples (n={len(sample_levels)}; scored={int(np.isfinite(mal_score).sum())})",
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(local.sum()),
        },
    }
    write_json(dest / "summary.json", summary)
    return summary


def write_finding(overall: dict, out_path: Path) -> None:
    t = overall["cohorts"]["tLung"]
    b = overall["cohorts"]["mBrain"]

    def block(tag: str, s: dict) -> str:
        da = s["da_malignant_cldn4"]
        sp = s["da_cldn4_median_split"]
        c = s["composition"]
        iface = c["interface_nhoods"]["spearman"]
        paired = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]
        w = paired["wilcoxon_high_vs_low"]
        sl = s["sample_level_malignant_CLDN4_vs_TNK_fraction"]
        drop = ", ".join(s["samples"]["dropped_lt10_malignant"]) or "none"
        return f"""### {tag}

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 author-malignant cells | **{s['samples']['n_with_malignant_ge10']}** of {s['samples']['n']} | Mean log1p-CP10k CLDN4 in author malignant cells (tS1/tS2/tS3 or Malignant cells). Dropped if <10 ({drop}). |
| Median split (Welch, secondary) | high / low | **{s['samples']['n_cldn4_high']} vs {s['samples']['n_cldn4_low']}** | Median of the scored samples. Ties at the median unlabeled. |

Do not cite {s['n_cells_in_graph']:,} graph cells or {s['graph']['n_nhoods']:,} neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n={s['honest_n']['disjoint_nhood_subset_n']}.

Graph: {s['n_cells_in_graph']:,} epithelium+immune cells; {s['n_malignant']:,} author-malignant; {s['n_tnk']:,} T/NK; {s['graph']['n_nhoods']:,} neighbourhoods of size {s['graph']['nhood_size']}; k={s['graph']['k']}, d={s['graph']['d']}, {s['graph']['n_hvg']} HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | {da['n_testable']} | {da['n_p_lt_0.05']} | {fmt(da.get('min_p'), sci=True)} | {fmt(da.get('min_SpatialFDR'))} | **{da['n_SpatialFDR_lt_0.1']}** | {da['n_SpatialFDR_lt_0.05']} | {da['n_BH_FDR_lt_0.1']} |
| Median split high vs low | {sp['n_testable']} | {sp['n_p_lt_0.05']} | {fmt(sp.get('min_p'), sci=True)} | {fmt(sp.get('min_SpatialFDR'))} | **{sp['n_SpatialFDR_lt_0.1']}** | {sp['n_SpatialFDR_lt_0.05']} | {sp['n_BH_FDR_lt_0.1']} |

Sample-level malignant CLDN4 vs T/NK fraction (no neighbourhoods): ρ={fmt(sl.get('rho'))}, p={fmt(sl.get('p'), sci=True)}, n={sl.get('n')}.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n={c['interface_nhoods']['n']}; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ={fmt(iface.get('rho'))}, p={fmt(iface.get('p'), sci=True)}, n={iface.get('n')}.

Disjoint interface subset: n={c['disjoint_interface']['n_disjoint_interface']} of {c['disjoint_interface']['n_disjoint_nhoods']} disjoint neighbourhoods; ρ={fmt(c['disjoint_interface']['spearman'].get('rho'))}, p={fmt(c['disjoint_interface']['spearman'].get('p'), sci=True)}.

Sample-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample): n={w.get('n')}; median T/NK high={fmt(w.get('median_a'))}, low={fmt(w.get('median_b'))}; Wilcoxon p={fmt(w.get('p'), sci=True)}. High nhoods={paired['n_high_nhoods']}, low nhoods={paired['n_low_nhoods']}.
"""

    text = f"""# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE131907)

Additive. Prior GSE131907 epithelial TACSTD2 work is taken as given. This
folder asks whether **transcriptional neighbourhoods** are differentially
abundant with **sample-level malignant CLDN4** on Kim et al. 2020
(GSE131907; PMID 32385277). GSE207422 Milo is a different agent and was
not run.

The independent unit is the **sample**, not the cell and not the overlapping
neighbourhood. This atlas is treatment-naive: there is **no ICI / MPR /
RECIST** label. Site is a confounder, so tLung and mBrain are separate
graphs. Do not pool them into one *n*.

**Not miloR.** DA = sample-level Spearman (primary) or Welch t-test
(median split) on neighbourhood proportions. SpatialFDR = miloR
`graphSpatialFDR` k-distance weights. Author `Cell_type` / `Cell_subtype`
are used. Malignant = `tS1` / `tS2` / `tS3` / `Malignant cells`.

Public matrix: {overall['n_cells_matrix']:,} cells, {overall['n_genes_matrix']:,} genes
(raw UMI text). Graph cells are epithelium + immune at one site.

## n and SpatialFDR

{block("tLung (primary, same tissue)", t)}

{block("mBrain (secondary, same-site metastasis)", b)}

## What this does not say

- It does not invent an ICI or MPR contrast. GEO has none.
- It does not treat 200k cells as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry
  (epithelium sits with epithelium). The sample-paired test is the only
  composition test whose unit is the patient.
- miloR / edgeR QLF numbers are not claimed.
- The 3 GB log2TPM text and EGA FASTQ were not used; UMI is the Milo input.

See `results/tLung/` and `results/mBrain/` (`summary.json`,
`da_malignant_cldn4.tsv`, `sample_scores.tsv`).
"""
    out_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE131907"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/gse131907_milo_cldn4/results"))
    ap.add_argument("--finding", type=Path, default=Path("methods/gse131907_milo_cldn4/FINDING.md"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--min-malignant", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    self_test()

    ann_path = args.datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = args.datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    series_path = args.datadir / "GSE131907_series_matrix.txt.gz"
    for p in (ann_path, umi_path, series_path):
        if not p.exists():
            raise SystemExit(f"missing {p}; run scripts/download.py")

    # Need cell IDs from the matrix header before annotation align.
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    print(f"matrix header cells={len(cell_ids)}", flush=True)
    ann = load_annotation(ann_path, series_path, cell_ids)
    cell_type = ann["Cell_type"].to_numpy()
    subtype = ann["Cell_subtype"].fillna("").to_numpy() if "Cell_subtype" in ann.columns else np.array([""] * len(ann))
    origin = ann["origin"].to_numpy()
    sample = ann["Sample"].to_numpy()
    keep_type = np.isin(cell_type, list(KEEP_TYPES))
    is_malig_all = np.isin(subtype, list(MALIGNANT_SUBTYPES))
    is_tnk_all = np.isin(cell_type, list(TNK_TYPES))
    cohort_masks = {c: (origin == c) & keep_type for c in COHORTS}
    for c, m in cohort_masks.items():
        print(
            f"  {c}: graph_cells={int(m.sum())} malignant={int((m & is_malig_all).sum())} "
            f"samples={pd.unique(sample[m]).size}",
            flush=True,
        )
        if m.sum() < 200:
            raise SystemExit(f"{c}: too few epithelium+immune cells ({int(m.sum())})")

    print("pass1: gene stats + markers", flush=True)
    cell_ids2, n_umi, gene_names, cohort_stats, marker_expr = stream_pass1(
        umi_path, set(MARKERS), cohort_masks
    )
    if cell_ids2 != cell_ids:
        raise SystemExit("matrix header changed between reads")
    n = len(cell_ids)
    print(f"cells={n} genes={len(gene_names)} median_nUMI={np.median(n_umi):.0f}", flush=True)
    if "CLDN4" not in marker_expr:
        raise SystemExit("CLDN4 is absent from the public UMI matrix")

    hvg_names = {}
    keep_genes = set(MARKERS)
    for c in COHORTS:
        mu, va = cohort_stats[c]
        idx = select_hvg(mu, va, n_hvg=args.n_hvg)
        hvg_names[c] = [gene_names[i] for i in idx]
        keep_genes.update(hvg_names[c])
        print(f"  {c} HVG={len(hvg_names[c])}", flush=True)

    keep_any = cohort_masks["tLung"] | cohort_masks["mBrain"]
    keep_idx = np.flatnonzero(keep_any)
    print(f"pass2: {len(keep_genes)} genes × {keep_idx.size} kept cells", flush=True)
    expr = stream_pass2(umi_path, keep_genes, keep_idx, n)
    for g, arr in marker_expr.items():
        expr.setdefault(g, arr[keep_idx].astype(np.float32, copy=False))

    sample_k = sample[keep_idx]
    lineage_k = cell_type[keep_idx]
    is_malig_k = is_malig_all[keep_idx]
    is_tnk_k = is_tnk_all[keep_idx]
    n_umi_k = n_umi[keep_idx]
    scale_k = np.where(n_umi_k > 0, 1e4 / n_umi_k, 0.0)
    cldn_k = np.log1p(expr["CLDN4"] * scale_k).astype(np.float32)
    origin_k = origin[keep_idx]

    # Inventory of author labels on the full matrix (honest, not used as n).
    inventory = {
        "n_cells_matrix": int(n),
        "n_genes_matrix": int(len(gene_names)),
        "cell_type_counts": {k: int(v) for k, v in pd.Series(cell_type).value_counts().items()},
        "origin_counts": {k: int(v) for k, v in pd.Series(origin).value_counts().items()},
        "n_author_malignant": int(is_malig_all.sum()),
        "umi_bytes": int(umi_path.stat().st_size),
        "skipped": [
            "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
            "EGAD00001005054",
            "GSE207422",
        ],
    }
    write_json(args.outdir / "inventory.json", inventory)
    ann[["Index", "Sample", "Cell_type", "Cell_subtype", "origin"]].assign(
        keep_in_graph=keep_type, author_malignant=is_malig_all
    ).groupby(["origin", "Cell_type"], dropna=False).size().reset_index(name="n").to_csv(
        args.outdir / "label_counts.tsv", sep="\t", index=False
    )

    cohort_summaries = {}
    for c in COHORTS:
        local = origin_k == c
        cohort_summaries[c] = run_cohort(
            name=c,
            local=local,
            sample_keep=sample_k,
            lineage_keep=lineage_k,
            is_malig_keep=is_malig_k,
            is_tnk_keep=is_tnk_k,
            cldn_keep=cldn_k,
            n_umi_keep=n_umi_k,
            expr=expr,
            hvg_names=hvg_names[c],
            outdir=args.outdir,
            k=args.k,
            d=args.d,
            prop=args.prop,
            n_hvg_requested=args.n_hvg,
            min_interface=args.min_interface,
            seed=args.seed,
            min_malig=args.min_malignant,
        )

    overall = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        "public_only": True,
        "author_barcode_labels": True,
        "malignant_definition": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
        "gene": "CLDN4",
        "additive_to": "prior GSE131907 TACSTD2 epithelial work; GSE207422 Milo is a different agent",
        "miloR": False,
        "da_model": "sample-level Spearman (malignant CLDN4) or Welch t-test (median split); not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "ici_labels": False,
        "n_cells_matrix": inventory["n_cells_matrix"],
        "n_genes_matrix": inventory["n_genes_matrix"],
        "cohorts": cohort_summaries,
        "inventory": inventory,
    }
    write_json(args.outdir / "summary.json", overall)
    write_finding(overall, args.finding)
    print(json.dumps({"finding": str(args.finding), "tLung": cohort_summaries["tLung"]["da_malignant_cldn4"], "mBrain": cohort_summaries["mBrain"]["da_malignant_cldn4"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
