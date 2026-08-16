#!/usr/bin/env python3
"""
Hunt: TACSTD2/CLDN4-high tumor spots vs T-cell spots in E-MTAB-13530 Visium.

Question
--------
In human NSCLC 10x Visium sections (E-MTAB-13530), are spots with high
tight-junction / epithelial-tumor signal (TACSTD2 = TROP2, CLDN4) spatially
SEGREGATED from T-cell spots? I.e. is there evidence of immune exclusion
(anti-colocalization) around TROP2/CLDN4-high tumor regions?

Design (honest, multi-sample)
-----------------------------
* 20 tumor ("T") sections across 8 patients.
* Per section we compute three complementary, permutation-based statistics:
    1. Spearman correlation between the TJ score and the T-cell score across
       spots (descriptive; spot-level p-values are anti-conservative because
       Visium spots are spatially autocorrelated -- reported but NOT used for
       the primary inference).
    2. Nearest-neighbour distance (microns): for each T-cell-high spot, the
       distance to the nearest TJ-high spot. Null = permute the TJ-high labels
       among spots (fixed count). Anti-colocalization => observed median NN
       distance LARGER than null. Reported as nn_ratio and nn_z.
    3. Neighbourhood-enrichment z-score (primary): count edges in the hex
       spatial graph that connect TJ-high(only) spots to T-cell-high(only)
       spots. Null = permute the categorical spot labels. NEGATIVE z =
       fewer contacts than expected = segregation / anti-colocalization.
       (Same idea as squidpy.gr.nhood_enrichment, implemented explicitly.)
* Primary inference is done ACROSS samples, not across spots: per-section
  effect sizes are averaged to one value per patient, and a Wilcoxon
  signed-rank test asks whether the patient-level effect differs from 0.
  This avoids pseudoreplication from autocorrelated spots.

Everything is written to results/hunt_visium_tj/.
"""

import io
import json
import os
import sys
import tarfile
import warnings
from glob import glob

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
DATA_DIR = "data/E-MTAB-13530"
OUT_DIR = "results/hunt_visium_tj"
RNG_SEED = 0
N_PERM = 1000
HIGH_QUANTILE = 0.80          # "high" = top 20% of the signature within a section
MIN_GENES_PER_SPOT = 200
MIN_COUNTS_PER_SPOT = 500
MIN_CELLS_PER_GENE = 3
SPOT_UM = 55.0                # Visium spot diameter in microns
RADIUS_FACTOR = 1.4          # neighbour radius = factor * median NN spot spacing

TJ_GENES = ["TACSTD2", "CLDN4"]
TCELL_GENES = ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "CD8A", "CD8B"]
# sanity-check reference signatures
EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19"]
IMMUNE_GENES = ["PTPRC"]

rng = np.random.default_rng(RNG_SEED)


def log(msg):
    print(msg, flush=True)


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def read_positions_from_tar(tar_path):
    """Return DataFrame indexed by barcode with pixel row/col from a spatial.tar."""
    with tarfile.open(tar_path) as tf:
        pos_member = next(m for m in tf.getmembers()
                          if m.name.endswith("tissue_positions_list.csv")
                          or m.name.endswith("tissue_positions.csv"))
        raw = tf.extractfile(pos_member).read().decode()
        scale_member = next(m for m in tf.getmembers()
                            if m.name.endswith("scalefactors_json.json"))
        scalefactors = json.loads(tf.extractfile(scale_member).read().decode())

    first = raw.splitlines()[0].split(",")
    has_header = "barcode" in first[0].lower()
    df = pd.read_csv(io.StringIO(raw), header=0 if has_header else None)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col",
                  "pxl_row_in_fullres", "pxl_col_in_fullres"]
    df = df.set_index("barcode")
    return df, scalefactors


def load_section(sample):
    h5 = os.path.join(DATA_DIR, f"{sample}-filtered_feature_bc_matrix.h5")
    tar = os.path.join(DATA_DIR, f"{sample}-spatial.tar")
    adata = sc.read_10x_h5(h5)
    adata.var_names_make_unique()
    pos, scalefactors = read_positions_from_tar(tar)
    common = adata.obs_names.intersection(pos.index)
    adata = adata[common].copy()
    pos = pos.loc[common]
    # pixel coordinates: (x=col, y=row)
    adata.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    adata.obs["array_row"] = pos["array_row"].to_numpy()
    adata.obs["array_col"] = pos["array_col"].to_numpy()
    um_per_px = SPOT_UM / float(scalefactors["spot_diameter_fullres"])
    adata.uns["um_per_px"] = um_per_px
    adata.obs["sample"] = sample
    adata.obs["patient"] = sample.split("_")[0]
    return adata


# ----------------------------------------------------------------------------
# Scoring
# ----------------------------------------------------------------------------
def preprocess_and_score(adata):
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_SPOT)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS_PER_SPOT)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    def score(name, genes):
        present = [g for g in genes if g in adata.var_names]
        if not present:
            adata.obs[name] = np.nan
            return present
        sc.tl.score_genes(adata, present, score_name=name, random_state=RNG_SEED)
        return present

    used = {
        "tj_score": score("tj_score", TJ_GENES),
        "tcell_score": score("tcell_score", TCELL_GENES),
        "epi_score": score("epi_score", EPI_GENES),
        "immune_score": score("immune_score", IMMUNE_GENES),
    }
    return adata, used


def build_graph(coords):
    """Hex neighbour graph from pixel coords; returns edge arrays (i, j) i<j and spacing."""
    tree = cKDTree(coords)
    d, _ = tree.query(coords, k=2)
    spacing = np.median(d[:, 1])
    radius = RADIUS_FACTOR * spacing
    pairs = tree.query_pairs(radius, output_type="ndarray")
    return pairs, spacing


def pair_stats(coords, um, n, a_only, b_only, both, pairs, prefix):
    """Anti-colocalization stats between group A-only and group B-only spots.

    Returns dict with <prefix>_nn_ratio, <prefix>_nn_z, <prefix>_nn_p_farther,
    <prefix>_enrich_z, <prefix>_enrich_p (+ raw counts).
    NEGATIVE enrich_z / nn_ratio>1 => segregation (anti-colocalization).
    """
    out = {}
    if a_only.sum() < 5 or b_only.sum() < 5 or len(pairs) == 0:
        for k in ["nn_obs_um", "nn_null_mean_um", "nn_ratio", "nn_z", "nn_p_farther",
                  "cross_edges_obs", "cross_edges_exp", "enrich_z", "enrich_p"]:
            out[f"{prefix}_{k}"] = np.nan
        return out

    # nearest-neighbour distance (B-only -> nearest A-only), microns
    b_coords = coords[b_only]
    obs_nn = np.median(cKDTree(coords[a_only]).query(b_coords, k=1)[0]) * um
    idx = np.arange(n)
    k_a = int(a_only.sum())
    null = np.empty(N_PERM)
    for t in range(N_PERM):
        perm = rng.choice(idx, size=k_a, replace=False)
        null[t] = np.median(cKDTree(coords[perm]).query(b_coords, k=1)[0]) * um
    nm, nsd = null.mean(), null.std(ddof=1)
    out[f"{prefix}_nn_obs_um"] = float(obs_nn)
    out[f"{prefix}_nn_null_mean_um"] = float(nm)
    out[f"{prefix}_nn_ratio"] = float(obs_nn / nm)
    out[f"{prefix}_nn_z"] = float((obs_nn - nm) / nsd) if nsd > 0 else np.nan
    out[f"{prefix}_nn_p_farther"] = float((1 + np.sum(null >= obs_nn)) / (N_PERM + 1))

    # neighbourhood enrichment on the hex graph
    i, j = pairs[:, 0], pairs[:, 1]
    lab = np.zeros(n, dtype=np.int8)
    lab[a_only] = 1
    lab[b_only] = 2
    lab[both] = 3

    def cross(labels):
        li, lj = labels[i], labels[j]
        return np.sum(((li == 1) & (lj == 2)) | ((li == 2) & (lj == 1)))

    obs_c = cross(lab)
    null = np.empty(N_PERM)
    for t in range(N_PERM):
        null[t] = cross(rng.permutation(lab))
    nm, nsd = null.mean(), null.std(ddof=1)
    out[f"{prefix}_cross_edges_obs"] = int(obs_c)
    out[f"{prefix}_cross_edges_exp"] = float(nm)
    out[f"{prefix}_enrich_z"] = float((obs_c - nm) / nsd) if nsd > 0 else np.nan
    centered = np.abs(null - nm)
    out[f"{prefix}_enrich_p"] = float((1 + np.sum(centered >= abs(obs_c - nm))) / (N_PERM + 1))
    return out


# ----------------------------------------------------------------------------
# Per-section statistics
# ----------------------------------------------------------------------------
def analyse_section(adata):
    n = adata.n_obs
    tj = adata.obs["tj_score"].to_numpy()
    tc = adata.obs["tcell_score"].to_numpy()
    epi = adata.obs["epi_score"].to_numpy()
    coords = adata.obsm["spatial"]
    um = adata.uns["um_per_px"]

    tc_high = tc >= np.quantile(tc, HIGH_QUANTILE)
    tj_high = tj >= np.quantile(tj, HIGH_QUANTILE)
    epi_high = epi >= np.quantile(epi, HIGH_QUANTILE)

    # primary groups: TROP2/CLDN4-high vs T-cell
    tj_a = tj_high & ~tc_high
    tj_b = tc_high & ~tj_high
    tj_both = tj_high & tc_high
    # epithelial control groups: EPCAM/KRT-high vs T-cell
    epi_a = epi_high & ~tc_high
    epi_b = tc_high & ~epi_high
    epi_both = epi_high & tc_high

    res = dict(
        sample=adata.obs["sample"].iloc[0],
        patient=adata.obs["patient"].iloc[0],
        n_spots=int(n),
        n_tj_high=int(tj_high.sum()),
        n_tcell_high=int(tc_high.sum()),
        n_epi_high=int(epi_high.sum()),
        n_tj_tcell_both=int(tj_both.sum()),
        um_per_px=float(um),
    )

    # descriptive spot-level correlations
    res["spearman_rho"] = float(spearmanr(tj, tc)[0])
    res["spearman_p_naive"] = float(spearmanr(tj, tc)[1])
    res["spearman_rho_epi_tcell"] = float(spearmanr(epi, tc)[0])
    # sanity checks that signatures capture the intended biology
    if adata.obs["epi_score"].notna().all():
        res["rho_tj_vs_epithelial"] = float(spearmanr(tj, epi)[0])
    if adata.obs["immune_score"].notna().all():
        res["rho_tcell_vs_ptprc"] = float(spearmanr(tc, adata.obs["immune_score"])[0])

    pairs, spacing = build_graph(coords)
    res["graph_spacing_um"] = float(spacing * um)
    res["n_edges"] = int(len(pairs))

    res.update(pair_stats(coords, um, n, tj_a, tj_b, tj_both, pairs, "tj"))
    res.update(pair_stats(coords, um, n, epi_a, epi_b, epi_both, pairs, "epi"))
    # specificity: does TROP2/CLDN4 segregate MORE than generic epithelium?
    res["enrich_z_tj_minus_epi"] = res["tj_enrich_z"] - res["epi_enrich_z"]

    masks = dict(tj_high=tj_high, tc_high=tc_high,
                 a_only=tj_a, b_only=tj_b, both=tj_both, coords=coords)
    return res, masks


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
def plot_section_map(sample, masks, out_path):
    coords = masks["coords"]
    fig, ax = plt.subplots(figsize=(6, 6))
    other = ~(masks["a_only"] | masks["b_only"] | masks["both"])
    ax.scatter(coords[other, 0], coords[other, 1], s=6, c="#dddddd", label="other", linewidths=0)
    ax.scatter(coords[masks["a_only"], 0], coords[masks["a_only"], 1], s=10, c="#d62728",
               label="TACSTD2/CLDN4-high", linewidths=0)
    ax.scatter(coords[masks["b_only"], 0], coords[masks["b_only"], 1], s=10, c="#1f77b4",
               label="T-cell-high", linewidths=0)
    ax.scatter(coords[masks["both"], 0], coords[masks["both"], 1], s=12, c="#9467bd",
               label="both", linewidths=0)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(f"{sample}: TJ-high vs T-cell-high spots")
    ax.axis("off")
    ax.legend(loc="upper right", fontsize=8, markerscale=1.5, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_summary(df, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    order = df.sort_values("tj_enrich_z")
    colors = ["#d62728" if z < 0 else "#7f7f7f" for z in order["tj_enrich_z"]]
    axes[0].barh(order["sample"], order["tj_enrich_z"], color=colors)
    axes[0].axvline(0, color="k", lw=1)
    axes[0].set_xlabel("neighbourhood-enrichment z\n(negative = segregation)")
    axes[0].set_title("TROP2/CLDN4-high vs T-cell (primary)")

    # specificity: TROP2/CLDN4 vs generic epithelial control, per section
    axes[1].scatter(df["epi_enrich_z"], df["tj_enrich_z"], c="#333")
    lo = min(df["epi_enrich_z"].min(), df["tj_enrich_z"].min())
    axes[1].plot([lo, 0], [lo, 0], "k--", lw=0.8, label="y = x")
    axes[1].set_xlabel("EPCAM/KRT-high vs T-cell  (enrich z)")
    axes[1].set_ylabel("TROP2/CLDN4-high vs T-cell  (enrich z)")
    axes[1].set_title("Specificity vs generic epithelium")
    axes[1].legend(fontsize=8)

    axes[2].hist(df["spearman_rho"], bins=12, color="#1f77b4", edgecolor="k")
    axes[2].axvline(0, color="k", lw=1)
    axes[2].set_xlabel("Spearman rho (TROP2/CLDN4 vs T-cell score, per spot)")
    axes[2].set_title("Per-section spot-level correlation")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "maps"), exist_ok=True)

    samples = sorted({os.path.basename(p).split("-filtered")[0]
                      for p in glob(os.path.join(DATA_DIR, "*-filtered_feature_bc_matrix.h5"))})
    log(f"Found {len(samples)} sections: {samples}")

    rows = []
    example_masks = {}
    for s in samples:
        try:
            adata = load_section(s)
            adata, used = preprocess_and_score(adata)
            if adata.n_obs < 100:
                log(f"[skip] {s}: only {adata.n_obs} spots after QC")
                continue
            res, masks = analyse_section(adata)
            rows.append(res)
            example_masks[s] = masks
            log(f"[ok] {s}: n={res['n_spots']} tj_high={res['n_tj_high']} "
                f"tcell_high={res['n_tcell_high']} tj_enrich_z={res['tj_enrich_z']:.2f} "
                f"epi_enrich_z={res['epi_enrich_z']:.2f} "
                f"tj_nn_ratio={res['tj_nn_ratio']:.3f} rho={res['spearman_rho']:.3f}")
        except Exception as e:
            log(f"[error] {s}: {type(e).__name__}: {e}")

    df = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
    df.to_csv(os.path.join(OUT_DIR, "per_section_stats.csv"), index=False)

    # ---- aggregate to patient level (avoid pseudoreplication) ----
    patient = df.groupby("patient").agg(
        n_sections=("sample", "size"),
        tj_enrich_z=("tj_enrich_z", "mean"),
        epi_enrich_z=("epi_enrich_z", "mean"),
        enrich_z_tj_minus_epi=("enrich_z_tj_minus_epi", "mean"),
        tj_nn_ratio=("tj_nn_ratio", "mean"),
        spearman_rho=("spearman_rho", "mean"),
    ).reset_index()
    patient.to_csv(os.path.join(OUT_DIR, "per_patient_stats.csv"), index=False)

    # ---- across-sample inference ----
    def signed_test(x, alternative="less"):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        if len(x) < 3 or np.allclose(x, 0):
            return dict(n=len(x), median=float(np.median(x)) if len(x) else np.nan,
                        W=np.nan, p=np.nan, alternative=alternative)
        try:
            W, p = wilcoxon(x, alternative=alternative)
        except ValueError:
            W, p = np.nan, np.nan
        return dict(n=int(len(x)), median=float(np.median(x)), W=float(W),
                    p=float(p), alternative=alternative)

    summary = {
        "n_sections": int(len(df)),
        "n_patients": int(patient.shape[0]),
        "config": dict(high_quantile=HIGH_QUANTILE, n_perm=N_PERM,
                       tj_genes=TJ_GENES, tcell_genes=TCELL_GENES,
                       epi_genes=EPI_GENES, radius_factor=RADIUS_FACTOR, seed=RNG_SEED),
        # section-level descriptive counts (TROP2/CLDN4 vs T-cell)
        "sections_tj_enrich_z_negative": int((df["tj_enrich_z"] < 0).sum()),
        "sections_tj_enrich_z_positive": int((df["tj_enrich_z"] > 0).sum()),
        "sections_tj_nn_ratio_gt1": int((df["tj_nn_ratio"] > 1).sum()),
        "sections_spearman_negative": int((df["spearman_rho"] < 0).sum()),
        # PRIMARY across-patient tests (H1: segregation)
        "patient_tj_enrich_z_vs0_less": signed_test(patient["tj_enrich_z"], alternative="less"),
        "patient_tj_nn_ratio_vs1_greater": signed_test(patient["tj_nn_ratio"].to_numpy() - 1.0, alternative="greater"),
        "patient_spearman_vs0_less": signed_test(patient["spearman_rho"], alternative="less"),
        # SPECIFICITY: is TROP2/CLDN4 segregation stronger than generic epithelium?
        "sections_tj_more_segregated_than_epi": int((df["enrich_z_tj_minus_epi"] < 0).sum()),
        "patient_tj_vs_epi_enrich_z_less": signed_test(patient["enrich_z_tj_minus_epi"], alternative="less"),
        "patient_tj_vs_epi_enrich_z_twosided": signed_test(patient["enrich_z_tj_minus_epi"], alternative="two-sided"),
        # section-level test (secondary; pseudoreplicated)
        "section_tj_enrich_z_vs0_less": signed_test(df["tj_enrich_z"], alternative="less"),
        # sanity checks
        "median_rho_tj_vs_epithelial": float(df.get("rho_tj_vs_epithelial", pd.Series(dtype=float)).median()),
        "median_rho_tcell_vs_ptprc": float(df.get("rho_tcell_vs_ptprc", pd.Series(dtype=float)).median()),
        "median_epi_enrich_z": float(df["epi_enrich_z"].median()),
        "median_tj_enrich_z": float(df["tj_enrich_z"].median()),
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    log("\n=== SUMMARY ===")
    log(json.dumps(summary, indent=2))

    # ---- plots ----
    plot_summary(df, os.path.join(OUT_DIR, "summary_overview.png"))
    for s, masks in example_masks.items():
        plot_section_map(s, masks, os.path.join(OUT_DIR, "maps", f"{s}.png"))

    log("\nWrote outputs to " + OUT_DIR)


if __name__ == "__main__":
    main()
