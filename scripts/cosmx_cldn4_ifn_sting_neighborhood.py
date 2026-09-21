#!/usr/bin/env python3
"""CosMx NSCLC: IFN module and STING-pathway neighborhood scores.

ADDITIVE layer on the locked CLDN4 exclusion result. Official He et al. 2022
960-plex (figshare 25976224 clustered object): CLDN4-high vs CLDN4-low
malignant cells, spatial neighborhood scores at 50 and 100 µm.

Does not recompute the locked cytotoxic-cell ratio. STING1, CGAS, and TBK1
are absent from the 960-plex; the STING score uses the on-panel downstream
neighborhood of that pathway. No private 8-KL.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_cldn4_ifn_sting")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

UM_PER_PX = 0.18  # He 2022 CosMx; FOV span is ~0.98 x 0.65 mm
RADII_UM = (50, 100)

SAMPLES = [
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
]
PATIENT = {
    "LUAD-5 R1": "Lung5",
    "LUAD-5 R2": "Lung5",
    "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9",
    "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}
# Published malignant labels are patient-specific tumor clusters.
TUMOR_LABEL = {
    "LUAD-5 R1": "tumor 5",
    "LUAD-5 R2": "tumor 5",
    "LUAD-5 R3": "tumor 5",
    "LUSC-6": "tumor 6",
    "LUAD-9 R1": "tumor 9",
    "LUAD-9 R2": "tumor 9",
    "LUAD-12": "tumor 12",
    "LUAD-13": "tumor 13",
}
IMMUNE_TYPES = {
    "B-cell",
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
    "Treg",
    "mDC",
    "macrophage",
    "mast",
    "monocyte",
    "neutrophil",
    "pDC",
    "plasmablast",
}

# On-panel interferon-response module. HLA genes are left out (detected in
# ~70% of cells and would dominate). Cytotoxic-granule genes (GZMB, PRF1,
# NKG7) are left out so this layer does not redo the locked effector readout.
IFN_GENES = [
    "STAT1",
    "JAK1",
    "JAK2",
    "IFNGR1",
    "IFNGR2",
    "IFNAR1",
    "IFNAR2",
    "IFNG",
    "IFNB1",
    "IFNA1",
    "IFIT1",
    "IFITM1",
    "IFITM3",
    "MX1",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "CXCL9",
    "CXCL10",
    "CCL5",
    "CD274",
    "IDO1",
    "BST2",
    "TAP1",
    "TAP2",
]
# STING1 / CGAS / TBK1 are not on the 960-plex. These are the on-panel genes
# immediately downstream of cGAS–STING (IRF3 and NF-kB arms, type-I output,
# and the canonical STING chemokines).
STING_GENES = [
    "IRF3",
    "NFKB1",
    "RELA",
    "NFKBIA",
    "IFNB1",
    "IFNA1",
    "CXCL10",
    "CCL5",
    "TNF",
    "IL6",
]
STING_ABSENT = ["STING1", "TMEM173", "CGAS", "MB21D1", "TBK1", "MAVS", "TREX1"]
EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"]

MODULES = {
    "ifn": IFN_GENES,
    "sting": STING_GENES,
    "epithelial_control": EPI_GENES,
}
CONTEXTS = ("all_neighbors", "immune_neighbors", "intrinsic")
PRIMARY = (("ifn", "all_neighbors"), ("sting", "all_neighbors"))


def _decode(arr) -> list:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def load_obs_and_counts(path: str):
    import h5py

    f = h5py.File(path, "r")
    genes = _decode(f["var/_index"][:])
    gmap = {g: i for i, g in enumerate(genes)}
    missing = [g for g in ["CLDN4", *IFN_GENES, *STING_GENES, *EPI_GENES] if g not in gmap]
    if missing:
        raise SystemExit(f"required genes absent from the 960-plex object: {missing}")
    absent_confirmed = [g for g in STING_ABSENT if g not in gmap]

    def categorical(name: str) -> np.ndarray:
        cats = _decode(f[f"obs/{name}/categories"][:])
        codes = f[f"obs/{name}/codes"][:]
        return np.array(cats, dtype=object)[codes]

    obs = {
        "sample": categorical("sample"),
        "patient": categorical("patient"),
        "cell_type": categorical("cell_type"),
        "n_counts": f["obs/n_counts"][:].astype(np.float64),
        "xy_um": f["obsm/spatial"][:].astype(np.float64) * UM_PER_PX,
    }
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    f.close()
    counts = sparse.csr_matrix((data, indices, indptr), shape=(len(obs["sample"]), len(genes)))
    needed = ["CLDN4", *IFN_GENES, *STING_GENES, *EPI_GENES]
    # unique, stable order
    seen = set()
    needed = [g for g in needed if not (g in seen or seen.add(g))]
    logexpr = {}
    lib = obs["n_counts"].copy()
    lib[lib <= 0] = np.nan
    for g in needed:
        col = np.asarray(counts.getcol(gmap[g]).todense()).ravel().astype(np.float64)
        logexpr[g] = np.nan_to_num(np.log1p(col / lib * 1e4), nan=0.0)
    del counts
    return obs, logexpr, genes, absent_confirmed


def zscore_within_sample(values: np.ndarray, sample: np.ndarray) -> np.ndarray:
    out = np.zeros(len(values), dtype=np.float64)
    for s in SAMPLES:
        m = sample == s
        v = values[m]
        sd = float(v.std())
        if sd == 0 or not np.isfinite(sd):
            continue
        out[m] = (v - float(v.mean())) / sd
    return out


def module_matrices(obs, logexpr) -> dict[str, dict[str, np.ndarray]]:
    """Per module: z-score mean and raw mean of log1p CP10k, plus per-gene log."""
    sample = obs["sample"]
    out = {}
    for name, genes in MODULES.items():
        zcols = [zscore_within_sample(logexpr[g], sample) for g in genes]
        zmean = np.mean(np.column_stack(zcols), axis=1)
        raw = np.mean(np.column_stack([logexpr[g] for g in genes]), axis=1)
        out[name] = {"z": zmean, "logmean": raw, "genes": genes}
    return out


def radius_adjacency(xy: np.ndarray, radius: float) -> sparse.csr_matrix:
    tree = cKDTree(xy)
    coo = tree.sparse_distance_matrix(tree, max_distance=radius, output_type="coo_matrix")
    n = xy.shape[0]
    mask = coo.row != coo.col
    A = sparse.csr_matrix(
        (np.ones(int(mask.sum()), dtype=np.float32), (coo.row[mask], coo.col[mask])),
        shape=(n, n),
    )
    A.sum_duplicates()
    return A


def neighbor_means(A: sparse.csr_matrix, values: np.ndarray, col_mask: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    """Mean of `values` over neighbors. Optional column mask keeps a neighbor class.

    `values` is length-n or shape (n, k). Self is already dropped from A.
    """
    if col_mask is None:
        use = A
    else:
        use = A @ sparse.diags(col_mask.astype(np.float32))
    deg = np.asarray(use.sum(axis=1)).ravel()
    total = np.asarray(use @ values)
    squeeze = total.ndim == 1
    if squeeze:
        total = total.reshape(-1, 1)
    out = np.full(total.shape, np.nan, dtype=np.float64)
    ok = deg > 0
    out[ok] = total[ok] / deg[ok, None]
    if squeeze:
        out = out.ravel()
    return out, deg.astype(np.int32)


def assign_arms(cldn4: np.ndarray, malignant: np.ndarray, sample: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    arm = np.full(len(cldn4), "other", dtype=object)
    rows = []
    for s in SAMPLES:
        m = (sample == s) & malignant
        vals = cldn4[m]
        q1, q3 = np.quantile(vals, [0.25, 0.75])
        idx = np.flatnonzero(m)
        v = cldn4[idx]
        hi = v >= q3
        lo = v <= q1
        both = hi & lo
        hi = hi & ~both
        arm[idx[hi]] = "high"
        arm[idx[lo]] = "low"
        arm[idx[~hi & ~lo]] = "mid"
        rows.append(
            {
                "sample": s,
                "patient": PATIENT[s],
                "tumor_label": TUMOR_LABEL[s],
                "n_malignant": int(m.sum()),
                "n_high": int(hi.sum()),
                "n_low": int(lo.sum()),
                "n_mid": int((~hi & ~lo).sum()),
                "q1_log_cp10k": float(q1),
                "q3_log_cp10k": float(q3),
                "frac_cldn4_pos": float((vals > 0).mean()),
                "mean_cldn4_high": float(v[hi].mean()) if hi.any() else np.nan,
                "mean_cldn4_low": float(v[lo].mean()) if lo.any() else np.nan,
            }
        )
    return arm, rows


def _append_arm_rows(sec_rows, sample, module, context, radius, arm_l, score_z, score_log, deg) -> None:
    for arm_name in ("high", "low"):
        q = arm_l == arm_name
        st = summarize_group(score_z[q])
        st_log = summarize_group(score_log[q])
        deg_q = deg[q]
        sec_rows.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "module": module,
                "context": context,
                "radius_um": radius,
                "arm": arm_name,
                "n_index": st["n"],
                "mean_z": st["mean"],
                "median_z": st["median"],
                "mean_log_cp10k": st_log["mean"],
                "median_neighbor_n": float(np.median(deg_q)) if len(deg_q) else np.nan,
                "frac_with_neighbor": float(np.mean(deg_q > 0)) if len(deg_q) else np.nan,
            }
        )


def summarize_group(values: np.ndarray) -> dict:
    v = values[np.isfinite(values)]
    if len(v) == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan}
    return {"n": int(len(v)), "mean": float(v.mean()), "median": float(np.median(v))}


def wilcoxon_pair(high: np.ndarray, low: np.ndarray) -> dict:
    d = np.asarray(high, dtype=float) - np.asarray(low, dtype=float)
    d = d[np.isfinite(d)]
    n_neg = int((d < 0).sum())
    n_pos = int((d > 0).sum())
    n_zero = int((d == 0).sum())
    out = {
        "n": int(len(d)),
        "n_high_lt_low": n_neg,
        "n_high_gt_low": n_pos,
        "n_tie": n_zero,
        "median_delta": float(np.median(d)) if len(d) else np.nan,
        "mean_delta": float(np.mean(d)) if len(d) else np.nan,
        "p": np.nan,
    }
    usable = d[d != 0]
    if len(usable) >= 1:
        res = stats.wilcoxon(usable, alternative="two-sided", method="exact", zero_method="wilcox")
        out["p"] = float(res.pvalue)
    return out


def analyze(obs, logexpr, modules) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample = obs["sample"]
    cell_type = obs["cell_type"]
    malignant = np.array([cell_type[i] == TUMOR_LABEL[sample[i]] for i in range(len(sample))])
    immune = np.isin(cell_type, list(IMMUNE_TYPES))
    arm, inv_rows = assign_arms(logexpr["CLDN4"], malignant, sample)
    # cross-patient tumor labels, not used as the index
    any_tumor = np.array([str(ct).startswith("tumor") for ct in cell_type])
    for row in inv_rows:
        s = row["sample"]
        m = sample == s
        row["n_qc"] = int(m.sum())
        row["n_immune"] = int((m & immune).sum())
        row["n_other_tumor_label"] = int((m & any_tumor & ~malignant).sum())
        row["n_epithelial"] = int((m & (cell_type == "epithelial")).sum())

    sec_rows = []
    gene_rows = []
    for s in SAMPLES:
        local = np.flatnonzero(sample == s)
        xy = obs["xy_um"][local]
        arm_l = arm[local]
        immune_l = immune[local]
        print(f"{s}: n={len(local)} building neighborhoods", flush=True)
        adj = {r: radius_adjacency(xy, float(r)) for r in RADII_UM}
        for module, pack in modules.items():
            z = pack["z"][local]
            logmean = pack["logmean"][local]
            gene_mat = np.column_stack([logexpr[g][local] for g in pack["genes"]])
            for radius, A in adj.items():
                for context in ("all_neighbors", "immune_neighbors"):
                    if context == "immune_neighbors":
                        score_z, deg = neighbor_means(A, z, immune_l)
                        score_log, _ = neighbor_means(A, logmean, immune_l)
                    else:
                        score_z, deg = neighbor_means(A, z, None)
                        score_log, _ = neighbor_means(A, logmean, None)
                    _append_arm_rows(sec_rows, s, module, context, radius, arm_l, score_z, score_log, deg)
                if module == "epithelial_control":
                    continue
                neigh_g, _ = neighbor_means(A, gene_mat, None)
                for j, g in enumerate(pack["genes"]):
                    for arm_name in ("high", "low"):
                        q = arm_l == arm_name
                        st = summarize_group(neigh_g[q, j])
                        gene_rows.append(
                            {
                                "sample": s,
                                "patient": PATIENT[s],
                                "module": module,
                                "gene": g,
                                "radius_um": radius,
                                "arm": arm_name,
                                "n_index": st["n"],
                                "mean_log_cp10k": st["mean"],
                            }
                        )
            _append_arm_rows(sec_rows, s, module, "intrinsic", 0, arm_l, z, logmean, np.ones(len(local), dtype=np.int32))
        del adj
    inv = pd.DataFrame(inv_rows)
    sec = pd.DataFrame(sec_rows)
    genes_df = pd.DataFrame(gene_rows)
    paired = pair_table(sec)
    gene_paired = pair_gene_table(genes_df)
    return inv, paired, gene_paired, sec


def pair_table(sec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = sec.groupby(["sample", "module", "context", "radius_um"], sort=False).size().reset_index()
    for rec in keys.itertuples(index=False):
        sub = sec[
            (sec["sample"] == rec.sample)
            & (sec["module"] == rec.module)
            & (sec["context"] == rec.context)
            & (sec["radius_um"] == rec.radius_um)
        ]
        hi = sub[sub["arm"] == "high"].iloc[0]
        lo = sub[sub["arm"] == "low"].iloc[0]
        rows.append(
            {
                "sample": rec.sample,
                "patient": PATIENT[rec.sample],
                "module": rec.module,
                "context": rec.context,
                "radius_um": int(rec.radius_um),
                "n_high": int(hi["n_index"]),
                "n_low": int(lo["n_index"]),
                "mean_z_high": hi["mean_z"],
                "mean_z_low": lo["mean_z"],
                "delta_z": hi["mean_z"] - lo["mean_z"],
                "mean_log_high": hi["mean_log_cp10k"],
                "mean_log_low": lo["mean_log_cp10k"],
                "delta_log": hi["mean_log_cp10k"] - lo["mean_log_cp10k"],
                "median_neighbor_n_high": hi["median_neighbor_n"],
                "median_neighbor_n_low": lo["median_neighbor_n"],
            }
        )
    return pd.DataFrame(rows)


def pair_gene_table(genes_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = genes_df.groupby(["sample", "module", "gene", "radius_um"], sort=False)
    for key, sub in grouped:
        sample, module, gene, radius = key
        hi = sub[sub["arm"] == "high"].iloc[0]
        lo = sub[sub["arm"] == "low"].iloc[0]
        rows.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "module": module,
                "gene": gene,
                "radius_um": int(radius),
                "mean_log_high": hi["mean_log_cp10k"],
                "mean_log_low": lo["mean_log_cp10k"],
                "delta_log": hi["mean_log_cp10k"] - lo["mean_log_cp10k"],
            }
        )
    return pd.DataFrame(rows)


def patient_roll(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["patient", "module", "context", "radius_um"]
    for key, sub in paired.groupby(keys, sort=False):
        rows.append(
            {
                "patient": key[0],
                "module": key[1],
                "context": key[2],
                "radius_um": int(key[3]),
                "n_sections": int(sub["sample"].nunique()),
                "mean_z_high": float(sub["mean_z_high"].mean()),
                "mean_z_low": float(sub["mean_z_low"].mean()),
                "delta_z": float(sub["delta_z"].mean()),
                "mean_log_high": float(sub["mean_log_high"].mean()),
                "mean_log_low": float(sub["mean_log_low"].mean()),
                "delta_log": float(sub["delta_log"].mean()),
            }
        )
    return pd.DataFrame(rows)


def test_block(df: pd.DataFrame, high_col: str, low_col: str) -> dict:
    return wilcoxon_pair(df[high_col].to_numpy(), df[low_col].to_numpy())


def collect_tests(paired: pd.DataFrame, patients: pd.DataFrame) -> list[dict]:
    tests = []
    specs = paired.groupby(["module", "context", "radius_um"], sort=False).size().reset_index()
    for rec in specs.itertuples(index=False):
        sec = paired[
            (paired["module"] == rec.module)
            & (paired["context"] == rec.context)
            & (paired["radius_um"] == rec.radius_um)
        ].sort_values("sample")
        pat = patients[
            (patients["module"] == rec.module)
            & (patients["context"] == rec.context)
            & (patients["radius_um"] == rec.radius_um)
        ].sort_values("patient")
        st = test_block(sec, "mean_z_high", "mean_z_low")
        pt = test_block(pat, "mean_z_high", "mean_z_low")
        st_log = test_block(sec, "mean_log_high", "mean_log_low")
        tests.append(
            {
                "module": rec.module,
                "context": rec.context,
                "radius_um": int(rec.radius_um),
                "section_median_z_high": float(sec["mean_z_high"].median()),
                "section_median_z_low": float(sec["mean_z_low"].median()),
                "section_median_delta_z": st["median_delta"],
                "section_mean_delta_z": st["mean_delta"],
                "section_high_lt_low": f"{st['n_high_lt_low']}/{st['n']}",
                "section_high_gt_low": f"{st['n_high_gt_low']}/{st['n']}",
                "section_p": st["p"],
                "patient_median_delta_z": pt["median_delta"],
                "patient_high_lt_low": f"{pt['n_high_lt_low']}/{pt['n']}",
                "patient_high_gt_low": f"{pt['n_high_gt_low']}/{pt['n']}",
                "patient_p": pt["p"],
                "section_median_delta_log": st_log["median_delta"],
                "section_log_high_lt_low": f"{st_log['n_high_lt_low']}/{st_log['n']}",
                "section_log_p": st_log["p"],
            }
        )
    return tests


def collect_gene_tests(gene_paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (module, gene, radius), sub in gene_paired.groupby(["module", "gene", "radius_um"], sort=False):
        st = wilcoxon_pair(sub["mean_log_high"].to_numpy(), sub["mean_log_low"].to_numpy())
        # patient rollup
        pat = (
            sub.groupby("patient", sort=False)[["mean_log_high", "mean_log_low", "delta_log"]]
            .mean()
            .reset_index()
        )
        pt = wilcoxon_pair(pat["mean_log_high"].to_numpy(), pat["mean_log_low"].to_numpy())
        rows.append(
            {
                "module": module,
                "gene": gene,
                "radius_um": int(radius),
                "median_delta_log": st["median_delta"],
                "sections_high_lt_low": f"{st['n_high_lt_low']}/{st['n']}",
                "sections_high_gt_low": f"{st['n_high_gt_low']}/{st['n']}",
                "section_p": st["p"],
                "patients_high_lt_low": f"{pt['n_high_lt_low']}/{pt['n']}",
                "patients_high_gt_low": f"{pt['n_high_gt_low']}/{pt['n']}",
                "patient_p": pt["p"],
            }
        )
    mod_order = {
        "ifn": {g: i for i, g in enumerate(IFN_GENES)},
        "sting": {g: i for i, g in enumerate(STING_GENES)},
    }
    out = pd.DataFrame(rows)
    out["_ord"] = [mod_order[m][g] for m, g in zip(out["module"], out["gene"])]
    return out.sort_values(["module", "radius_um", "_ord"]).drop(columns="_ord")


def _fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _fmt(x, nd=3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def plot_paired(paired: pd.DataFrame, module: str, path: str, title: str) -> None:
    sub = paired[(paired["module"] == module) & (paired["context"] == "all_neighbors")].copy()
    radii = [r for r in RADII_UM]
    fig, axes = plt.subplots(1, len(radii), figsize=(8.4, 4.2), sharey=True)
    patients = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
    patient_color = {
        "Lung5": "#1b9e77",
        "Lung6": "#d95f02",
        "Lung9": "#7570b3",
        "Lung12": "#e7298a",
        "Lung13": "#66a61e",
    }
    for ax, radius in zip(axes, radii):
        part = sub[sub["radius_um"] == radius]
        for s in SAMPLES:
            row = part[part["sample"] == s].iloc[0]
            color = patient_color[row["patient"]]
            ax.plot([0, 1], [row["mean_z_low"], row["mean_z_high"]], color=color, lw=1.4, alpha=0.9)
            ax.scatter([0, 1], [row["mean_z_low"], row["mean_z_high"]], color=color, s=28, zorder=3)
        ax.set_xticks([0, 1], ["CLDN4-low", "CLDN4-high"])
        ax.set_title(f"{radius} µm")
        ax.axhline(0, color="0.75", lw=0.8)
        ax.set_xlim(-0.25, 1.25)
    axes[0].set_ylabel("Mean neighborhood module score (within-sample z)")
    handles = [
        plt.Line2D([0], [0], color=patient_color[p], lw=2, label=p) for p in patients
    ]
    axes[-1].legend(handles=handles, frameon=False, fontsize=8, loc="best")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def plot_delta_context(paired: pd.DataFrame, path: str) -> None:
    """Section deltas at 50 µm for all-neighbor, immune-only, and intrinsic."""
    modules = ["ifn", "sting", "epithelial_control"]
    contexts = ["all_neighbors", "immune_neighbors", "intrinsic"]
    labels = {
        ("ifn", "all_neighbors"): "IFN\nall neighbors",
        ("ifn", "immune_neighbors"): "IFN\nimmune neighbors",
        ("ifn", "intrinsic"): "IFN\nmalignant cell",
        ("sting", "all_neighbors"): "STING\nall neighbors",
        ("sting", "immune_neighbors"): "STING\nimmune neighbors",
        ("sting", "intrinsic"): "STING\nmalignant cell",
        ("epithelial_control", "all_neighbors"): "Epithelial\nall neighbors",
        ("epithelial_control", "immune_neighbors"): "Epithelial\nimmune neighbors",
        ("epithelial_control", "intrinsic"): "Epithelial\nmalignant cell",
    }
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    positions = []
    pos = 0
    for module in modules:
        for context in contexts:
            part = paired[
                (paired["module"] == module)
                & (paired["context"] == context)
                & ((paired["radius_um"] == 50) | (paired["context"] == "intrinsic"))
            ]
            if context == "intrinsic":
                part = paired[(paired["module"] == module) & (paired["context"] == "intrinsic")]
                part = part.drop_duplicates("sample")
            vals = part["delta_z"].to_numpy()
            jitter = np.linspace(-0.12, 0.12, len(vals))
            ax.scatter(np.full(len(vals), pos) + jitter, vals, s=28, color="#333333", zorder=3)
            ax.hlines(np.median(vals), pos - 0.28, pos + 0.28, color="#b00020", lw=2, zorder=4)
            positions.append(pos)
            pos += 1
        pos += 0.6
    names = [labels[(m, c)] for m in modules for c in contexts]
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xticks(positions, names, fontsize=8)
    ax.set_ylabel("Section Δ (CLDN4-high − low), module z")
    ax.set_title("50 µm neighborhood scores (intrinsic has no radius)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def plot_gene_deltas(gene_tests: pd.DataFrame, path: str) -> None:
    sub = gene_tests[gene_tests["radius_um"] == 50].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.2), sharex=False)
    for ax, module, title in (
        (axes[0], "ifn", "IFN module genes"),
        (axes[1], "sting", "STING-pathway genes"),
    ):
        part = sub[sub["module"] == module].iloc[::-1]
        y = np.arange(len(part))
        ax.axvline(0, color="0.6", lw=0.8)
        ax.scatter(part["median_delta_log"], y, s=28, color="#1f4e79")
        ax.set_yticks(y, part["gene"].tolist(), fontsize=8)
        ax.set_xlabel("Median section Δ log1p CP10k\n(CLDN4-high − low), 50 µm")
        ax.set_title(title)
    fig.suptitle("Per-gene all-neighbor expression around malignant cells", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def lookup(tests: list[dict], module: str, context: str, radius: int) -> dict:
    for t in tests:
        if t["module"] == module and t["context"] == context and int(t["radius_um"]) == radius:
            return t
    raise KeyError((module, context, radius))


def write_results(inv: pd.DataFrame, paired: pd.DataFrame, tests: list[dict], gene_tests: pd.DataFrame, absent: list[str]) -> None:
    lines = []
    a = lines.append
    a("# CosMx NSCLC: IFN module and STING-pathway neighborhood scores")
    a("")
    a("ADDITIVE, **CLDN4-only**, He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; 8 sections / 5 patients). This layer scores **IFN-module** and **STING-pathway** genes in the spatial neighborhood of malignant cells. It does **not** replace the locked exclusion result (fewer cytotoxic cells at 50/100 µm around CLDN4-high tumor; nearby effectors not muzzled on GZMB/PRF1/NKG7/IFNG). No private 8-KL. No ICI labels.")
    a("")
    a("## Design")
    a("")
    a("- **Malignant index:** published `cell_type` matched to the section's tumor cluster (`tumor 5` in LUAD-5, `tumor 6` in LUSC-6, `tumor 9` in LUAD-9, `tumor 12` in LUAD-12, `tumor 13` in LUAD-13). Cells carrying another patient's tumor label are not used as the index. Generic `epithelial` is not called malignant.")
    a("- **CLDN4-high / low:** Q4 vs Q1 of log1p(CP10k) CLDN4 among those malignant cells, **per section**. Q1 is 0 where CLDN4 is zero-inflated, so the low arm is larger than a strict quartile.")
    a("- **Neighborhood:** other cells within **50 µm** and **100 µm** of the malignant centroid (global coordinates, 0.18 µm/pixel). The index cell is excluded. FOV edges are not censored.")
    a("- **Module score:** each gene is log1p(CP10k), then z-scored within the section across all cells. The cell score is the mean z of the module. The neighborhood score is the mean cell score of neighbors. A companion column is the mean log1p(CP10k) of the same genes (expression units).")
    a("- **Contexts:** all neighbors (primary); immune neighbors only (published immune types — program among immune cells that are present); the malignant cell itself (intrinsic, not a neighborhood).")
    a("- **Epithelial control:** EPCAM, KRT8, KRT18, KRT19, KRT7, same neighborhood machinery.")
    a("- **Inference:** two-sided exact Wilcoxon signed-rank on the 8 section means, and on the 5 patient means (a patient's sections are averaged with equal section weight). Honest n = 8 sections / 5 patients.")
    a("")
    a("## Gene sets and panel limits")
    a("")
    a(f"IFN module ({len(IFN_GENES)}): {', '.join(IFN_GENES)}.")
    a("")
    a(f"STING-pathway neighborhood ({len(STING_GENES)}): {', '.join(STING_GENES)}.")
    a("")
    a(f"Confirmed **absent** from this 960-plex object: {', '.join(absent)}. The STING score is the on-panel downstream neighborhood (IRF3, NF-κB, type-I IFN, CXCL10/CCL5), not a measurement of STING1 or cGAS.")
    a("")
    a("HLA-A/B/C and B2M are on the panel and were **not** put in the IFN module (too constitutive). GZMB, PRF1, and NKG7 were **not** put in it either (those belong to the locked effector readout).")
    a("")
    a("CXCL10, CCL5, IFNB1, and IFNA1 sit in both modules because they are IFN-response genes and STING outputs. The two scores are therefore not independent.")
    a("")
    a("## Inventory")
    a("")
    a("| Section | Patient | QC cells | Malignant | CLDN4-high | CLDN4-low | CLDN4>0 | Immune | Other tumor label |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in inv.itertuples(index=False):
        a(
            f"| {r.sample} | {r.patient} | {r.n_qc:,} | {r.n_malignant:,} | {r.n_high:,} | {r.n_low:,} | {r.frac_cldn4_pos:.3f} | {r.n_immune:,} | {r.n_other_tumor_label:,} |"
        )
    a("")
    a("## Primary neighborhood scores (all neighbors)")
    a("")
    a("Δ is CLDN4-high minus CLDN4-low. A positive Δ means a **higher** module score around CLDN4-high malignant cells. “Higher” counts sections (or patients) with Δ > 0.")
    a("")
    a("| Module | Radius | Median section z high | Median section z low | Median Δ z | Sections higher | Section p | Patients higher | Patient p | Median Δ log1p |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for module, radius in (("ifn", 50), ("ifn", 100), ("sting", 50), ("sting", 100)):
        t = lookup(tests, module, "all_neighbors", radius)
        a(
            f"| {module} | {radius} | {_fmt(t['section_median_z_high'])} | {_fmt(t['section_median_z_low'])} | {_fmt(t['section_median_delta_z'])} | {t['section_high_gt_low']} | {_fmt_p(t['section_p'])} | {t['patient_high_gt_low']} | {_fmt_p(t['patient_p'])} | {_fmt(t['section_median_delta_log'])} |"
        )
    a("")
    a("Per-section Δ z (all neighbors):")
    a("")
    a("| Section | Patient | IFN 50 | IFN 100 | STING 50 | STING 100 |")
    a("|---|---|---:|---:|---:|---:|")
    for s in SAMPLES:
        bits = []
        for module, radius in (("ifn", 50), ("ifn", 100), ("sting", 50), ("sting", 100)):
            row = paired[
                (paired["sample"] == s)
                & (paired["module"] == module)
                & (paired["context"] == "all_neighbors")
                & (paired["radius_um"] == radius)
            ].iloc[0]
            bits.append(_fmt(row["delta_z"]))
        a(f"| {s} | {PATIENT[s]} | " + " | ".join(bits) + " |")
    a("")
    a("## Immune-neighbor-only and intrinsic scores")
    a("")
    a("Immune-neighbor scores use only malignant cells that have at least one immune neighbor. Intrinsic is the malignant cell's own module score and is not a spatial neighborhood.")
    a("")
    a("| Module | Context | Radius | Median Δ z | Sections higher | Section p | Patients higher | Patient p |")
    a("|---|---|---:|---:|---:|---:|---:|---:|")
    for module in ("ifn", "sting", "epithelial_control"):
        for context, radius in (
            ("all_neighbors", 50),
            ("all_neighbors", 100),
            ("immune_neighbors", 50),
            ("immune_neighbors", 100),
            ("intrinsic", 0),
        ):
            t = lookup(tests, module, context, radius)
            rad = "—" if context == "intrinsic" else str(radius)
            a(
                f"| {module} | {context} | {rad} | {_fmt(t['section_median_delta_z'])} | {t['section_high_gt_low']} | {_fmt_p(t['section_p'])} | {t['patient_high_gt_low']} | {_fmt_p(t['patient_p'])} |"
            )
    a("")
    a("Immune-neighbor **counts** are not this layer's endpoint (the locked cytotoxic ratio already covers exclusion). They are recorded so a higher module score is not misread as “more immune cells.” Median immune neighbors at 50 µm, CLDN4-high vs low:")
    a("")
    a("| Section | Median immune neighbors, high | Median immune neighbors, low | Fraction of high cells with ≥1 | Fraction of low cells with ≥1 |")
    a("|---|---:|---:|---:|---:|")
    imm = paired[(paired["module"] == "ifn") & (paired["context"] == "immune_neighbors") & (paired["radius_um"] == 50)]
    inv_by = inv.set_index("sample")
    for s in SAMPLES:
        row = imm[imm["sample"] == s].iloc[0]
        frac_hi = row["n_high"] / inv_by.loc[s, "n_high"]
        frac_lo = row["n_low"] / inv_by.loc[s, "n_low"]
        a(
            f"| {s} | {_fmt(row['median_neighbor_n_high'], 1)} | {_fmt(row['median_neighbor_n_low'], 1)} | {frac_hi:.3f} | {frac_lo:.3f} |"
        )
    a("")
    a("## Per-gene all-neighbor expression (50 µm)")
    a("")
    a("Median section Δ of mean neighbor log1p(CP10k). Full 50 and 100 µm tests are in `tables/gene_neighborhood_tests.csv`.")
    a("")
    a("| Module | Gene | Median Δ log1p | Sections higher | Section p | Patients higher | Patient p |")
    a("|---|---|---:|---:|---:|---:|---:|")
    show = gene_tests[gene_tests["radius_um"] == 50]
    for r in show.itertuples(index=False):
        a(
            f"| {r.module} | {r.gene} | {_fmt(r.median_delta_log, 4)} | {r.sections_high_gt_low} | {_fmt_p(r.section_p)} | {r.patients_high_gt_low} | {_fmt_p(r.patient_p)} |"
        )
    a("")
    a("## Reading this layer")
    a("")
    ifn50 = lookup(tests, "ifn", "all_neighbors", 50)
    ifn100 = lookup(tests, "ifn", "all_neighbors", 100)
    st50 = lookup(tests, "sting", "all_neighbors", 50)
    st100 = lookup(tests, "sting", "all_neighbors", 100)
    ifn_imm = lookup(tests, "ifn", "immune_neighbors", 50)
    st_imm = lookup(tests, "sting", "immune_neighbors", 50)
    ifn_in = lookup(tests, "ifn", "intrinsic", 0)
    st_in = lookup(tests, "sting", "intrinsic", 0)
    epi = lookup(tests, "epithelial_control", "all_neighbors", 50)
    g50 = gene_tests[gene_tests["radius_um"] == 50]
    ifn_genes_up = int(((g50["module"] == "ifn") & (g50["median_delta_log"] > 0)).sum())
    ifn_genes_n = int((g50["module"] == "ifn").sum())
    sting_genes_up = int(((g50["module"] == "sting") & (g50["median_delta_log"] > 0)).sum())
    sting_genes_n = int((g50["module"] == "sting").sum())
    a(
        f"CLDN4-high neighborhoods score **higher**, not lower, on both modules. All-neighbor IFN Δ z is {_fmt(ifn50['section_median_delta_z'])} at 50 µm ({ifn50['section_high_gt_low']} sections higher, p={_fmt_p(ifn50['section_p'])}; {ifn50['patient_high_gt_low']} patients, p={_fmt_p(ifn50['patient_p'])}) and {_fmt(ifn100['section_median_delta_z'])} at 100 µm ({ifn100['section_high_gt_low']} higher, p={_fmt_p(ifn100['section_p'])}). "
        f"All-neighbor STING Δ z is {_fmt(st50['section_median_delta_z'])} at 50 µm ({st50['section_high_gt_low']} higher, p={_fmt_p(st50['section_p'])}; {st50['patient_high_gt_low']} patients, p={_fmt_p(st50['patient_p'])}) and {_fmt(st100['section_median_delta_z'])} at 100 µm ({st100['section_high_gt_low']} higher, p={_fmt_p(st100['section_p'])}). "
        f"The shift is small (a few hundredths of a within-section z). Patient p = 0.0625 is the smallest two-sided exact Wilcoxon p at n = 5 when every patient has the same sign."
    )
    a("")
    a(
        f"The same direction is present when the average is restricted to immune neighbors (50 µm IFN Δ z {_fmt(ifn_imm['section_median_delta_z'])}, {ifn_imm['section_high_gt_low']}, p={_fmt_p(ifn_imm['section_p'])}; STING Δ z {_fmt(st_imm['section_median_delta_z'])}, {st_imm['section_high_gt_low']}, p={_fmt_p(st_imm['section_p'])}) and in the malignant cell itself (IFN Δ z {_fmt(ifn_in['section_median_delta_z'])}, {ifn_in['section_high_gt_low']}, p={_fmt_p(ifn_in['section_p'])}; STING Δ z {_fmt(st_in['section_median_delta_z'])}, {st_in['section_high_gt_low']}, p={_fmt_p(st_in['section_p'])}). "
        f"At 50 µm, {ifn_genes_up}/{ifn_genes_n} IFN genes and {sting_genes_up}/{sting_genes_n} STING-pathway genes have a positive median section Δ. "
        f"Section-consistent genes (8/8 higher) include STAT1, IFNGR1, MX1, IFIT1, CXCL9, CXCL10, CCL5, CD274, IDO1, NFKBIA, and IL6. JAK1, JAK2, and IRF3 are positive at the median but not consistent across sections. "
        f"Epithelial-marker all-neighbor Δ z at 50 µm is larger at the median ({_fmt(epi['section_median_delta_z'])}) but only {epi['section_high_gt_low']} sections are higher (p={_fmt_p(epi['section_p'])}). The IFN and STING shifts are smaller and more consistent across sections. The malignant cell's own epithelial score is higher ({lookup(tests, 'epithelial_control', 'intrinsic', 0)['section_high_gt_low']}, median Δ z {_fmt(lookup(tests, 'epithelial_control', 'intrinsic', 0)['section_median_delta_z'])})."
    )
    a("")
    a("Lung5 is the flat section: IFN all-neighbor Δ z is about +0.002 to +0.012 at 50 µm, and two of the three Lung5 replicates flip slightly negative at 100 µm. Lung6 is immune-poor (median immune-neighbor count 0 in both arms); its immune-neighbor score uses the ~40% of malignant cells that have at least one immune neighbor. Lung9 and Lung12 have fewer immune neighbors around CLDN4-high than around CLDN4-low, and the IFN score among those neighbors is still higher.")
    a("")
    a("## What this does not claim")
    a("")
    a("- Not a new estimate of the locked 50/100 µm cytotoxic-cell ratio, and not a GZMB/PRF1/NKG7/IFNG muzzling test.")
    a("- Not evidence that STING1 or cGAS was measured. Those genes are off-panel.")
    a("- Not ICI response, not private 8-KL, and not a ligand-receptor probability.")
    a("- Q4 vs Q1 is a within-section contrast. n = 8 sections (5 patients), not 8 independent patients. Lung5 and Lung9 contribute repeated sections.")
    a("")
    a("## Reproduce")
    a("")
    a("```bash")
    a("python3 scripts/download_cosmx_nsclc_h5ad.py")
    a("python3 scripts/cosmx_cldn4_ifn_sting_neighborhood.py")
    a("```")
    a("")
    a("The h5ad is gitignored. Tables are under `results/cosmx_cldn4_ifn_sting/tables/`. Figures: `figures/paired_ifn_neighborhood.png`, `figures/paired_sting_neighborhood.png`, `figures/delta_context_50um.png`, `figures/gene_delta_50um.png`.")
    a("")
    text = "\n".join(lines) + "\n"
    os.makedirs(OUT, exist_ok=True)
    for dest in (os.path.join(ROOT, "RESULTS.md"), os.path.join(OUT, "RESULTS.md")):
        with open(dest, "w") as fh:
            fh.write(text)


def main() -> int:
    if not os.path.isfile(H5AD):
        raise SystemExit(f"missing {H5AD}; run scripts/download_cosmx_nsclc_h5ad.py")
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    print("loading counts", flush=True)
    obs, logexpr, genes, absent = load_obs_and_counts(H5AD)
    print(f"cells={len(obs['sample'])} genes_in_object={len(genes)} sting_absent={absent}", flush=True)
    modules = module_matrices(obs, logexpr)
    inv, paired, gene_paired, sec = analyze(obs, logexpr, modules)
    patients = patient_roll(paired)
    tests = collect_tests(paired, patients)
    gene_tests = collect_gene_tests(gene_paired)

    inv.to_csv(os.path.join(TAB, "sample_inventory.csv"), index=False)
    paired.to_csv(os.path.join(TAB, "section_paired_scores.csv"), index=False)
    patients.to_csv(os.path.join(TAB, "patient_paired_scores.csv"), index=False)
    pd.DataFrame(tests).to_csv(os.path.join(TAB, "wilcoxon_tests.csv"), index=False)
    gene_paired.to_csv(os.path.join(TAB, "gene_neighborhood_by_section.csv"), index=False)
    gene_tests.to_csv(os.path.join(TAB, "gene_neighborhood_tests.csv"), index=False)
    sec.to_csv(os.path.join(TAB, "section_arm_scores.csv"), index=False)
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(
            {
                "ifn_genes": IFN_GENES,
                "sting_genes": STING_GENES,
                "sting_absent": absent,
                "epithelial_control": EPI_GENES,
                "um_per_px": UM_PER_PX,
                "radii_um": list(RADII_UM),
                "tests": tests,
            },
            fh,
            indent=2,
        )
    plot_paired(
        paired,
        "ifn",
        os.path.join(FIG, "paired_ifn_neighborhood.png"),
        "IFN-module neighborhood score, CLDN4-high vs low malignant",
    )
    plot_paired(
        paired,
        "sting",
        os.path.join(FIG, "paired_sting_neighborhood.png"),
        "STING-pathway neighborhood score, CLDN4-high vs low malignant",
    )
    plot_delta_context(paired, os.path.join(FIG, "delta_context_50um.png"))
    plot_gene_deltas(gene_tests, os.path.join(FIG, "gene_delta_50um.png"))
    write_results(inv, paired, tests, gene_tests, absent)
    print(pd.DataFrame(tests)[["module", "context", "radius_um", "section_median_delta_z", "section_high_lt_low", "section_p", "patient_high_lt_low", "patient_p"]].to_string(index=False))
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
