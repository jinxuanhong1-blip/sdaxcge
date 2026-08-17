#!/usr/bin/env python3
"""CellChat-style LR inference on public GSE207422 — CLDN4 only.

Compare CLDN4-high vs CLDN4-low epithelial/malignant cells for
outgoing signals to T/NK and incoming signals from T/NK.

Prior TACSTD2 CellChat (methods/scrna_cellchat) is taken as given.
This folder does not re-split on TACSTD2.

Probability is the published CellChat Hill / mass-action on truncated
means (CellChatDB v2) plus a high/low label permutation. The CellChat
R package is not run. Networks are not drawn for non-significant pairs.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
NBOOT = 100
SEED = 1
MIN_GROUP = 25
MARKER = "CLDN4"

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "IFNG", "TNF"]

# Curated symbols used only to score which split matches the immune-cold / barrier trend.
BARRIER_LIGANDS = {
    "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "JAM2", "JAM3",
    "CEACAM1", "CEACAM5", "CEACAM6", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4",
    "PVR", "EPCAM", "DSG2", "DSC2", "CADM1",
}
INHIB_LIGANDS = {
    "CD274", "PDCD1LG2", "LGALS9", "HLA-E", "HLA-G", "HLA-F", "TGFB1", "TGFB2",
    "TGFB3", "CD80", "CD86", "CD276", "VSIR", "PVR", "NECTIN2", "CD47", "CDH1",
}
RECRUIT_LIGANDS = {
    "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL3", "CCL4", "IL15",
    "IL2", "IL18", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
}
ATTACK_LIGANDS = {"IFNG", "TNF", "FASLG", "TNFSF10", "LTA"}


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"})
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(rec.ligand_symbol if rec.ligand_symbol == rec.ligand_symbol else rec.ligand)
        recp = parse_symbols(rec.receptor_symbol if rec.receptor_symbol == rec.receptor_symbol else rec.receptor)
        if not lig:
            lig = parse_symbols(rec.ligand)
        if not recp:
            recp = parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        if any(g not in matrix_genes for g in lig + recp):
            continue
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


def wanted_genes(lr: pd.DataFrame) -> set[str]:
    genes = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    for rec in lr.itertuples(index=False):
        genes.update(rec.ligand_genes)
        genes.update(rec.receptor_genes)
    return genes


def sample_from_barcode(bc: str) -> str:
    parts = bc.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else bc


def stream_matrix(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def geom_mean_rows(mat: np.ndarray) -> np.ndarray:
    """Geometric mean over genes (rows). Zero if any gene mean is 0."""
    if mat.ndim == 1:
        return mat
    if mat.shape[0] == 1:
        return mat[0]
    out = np.exp(np.mean(np.log(np.clip(mat, 1e-12, None)), axis=0))
    out[np.any(mat <= 0, axis=0)] = 0.0
    return out


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def trim_mean_fast(sub: np.ndarray, proportiontocut: float = TRIM) -> np.ndarray:
    """Same definition as scipy.stats.trim_mean along axis=1."""
    n = sub.shape[1]
    if n == 0:
        return np.zeros(sub.shape[0], dtype=np.float64)
    if n == 1:
        return sub[:, 0].astype(np.float64)
    k = int(n * proportiontocut)
    if k == 0:
        return sub.mean(axis=1)
    s = np.sort(sub, axis=1)
    return s[:, k : n - k].mean(axis=1)


def group_trim_means(expr: np.ndarray, pos: np.ndarray, labels: np.ndarray, groups: list[str]):
    n_g = expr.shape[0]
    means = np.zeros((n_g, len(groups)), dtype=np.float64)
    props = np.zeros((n_g, len(groups)), dtype=np.float64)
    counts = {}
    for j, g in enumerate(groups):
        idx = np.flatnonzero(labels == g)
        counts[g] = int(idx.size)
        if idx.size == 0:
            continue
        sub = expr[:, idx]
        means[:, j] = trim_mean_fast(sub)
        props[:, j] = pos[:, idx].mean(axis=1)
    return means, props, counts


def complex_value(gene_means: np.ndarray, gene_props: np.ndarray, gene_index: dict[str, int], subunits: tuple[str, ...]):
    ix = [gene_index[g] for g in subunits]
    mu = geom_mean_rows(gene_means[ix])
    pr = gene_props[ix].min(axis=0) if len(ix) > 1 else gene_props[ix[0]]
    return mu, pr


def pair_specs(lr: pd.DataFrame, gene_index: dict[str, int]):
    specs = []
    for rec in lr.itertuples(index=False):
        specs.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": rec.ligand_genes,
                "receptor_genes": rec.receptor_genes,
                "lig_ix": [gene_index[g] for g in rec.ligand_genes],
                "rec_ix": [gene_index[g] for g in rec.receptor_genes],
            }
        )
    return specs


def score_pairs(specs, gene_means, gene_props, groups, counts, src_tgt):
    gpos = {g: i for i, g in enumerate(groups)}
    cache = {}

    def cached(ix_key, ix):
        hit = cache.get(ix_key)
        if hit is None:
            mu = geom_mean_rows(gene_means[ix])
            pr = gene_props[ix].min(axis=0) if len(ix) > 1 else gene_props[ix[0]]
            hit = (mu, pr)
            cache[ix_key] = hit
        return hit

    rows = []
    probs = []
    for spec in specs:
        lig_mu, lig_pr = cached(("L", tuple(spec["lig_ix"])), spec["lig_ix"])
        rec_mu, rec_pr = cached(("R", tuple(spec["rec_ix"])), spec["rec_ix"])
        for src, tgt in src_tgt:
            if src not in gpos or tgt not in gpos:
                continue
            i, j = gpos[src], gpos[tgt]
            if counts.get(src, 0) < MIN_GROUP or counts.get(tgt, 0) < MIN_GROUP:
                continue
            detected = (float(lig_pr[i]) >= EXPR_PROP) and (float(rec_pr[j]) >= EXPR_PROP)
            prob = hill_prob(float(lig_mu[i]), float(rec_mu[j])) if detected else 0.0
            rows.append(
                {
                    "interaction_name": spec["interaction_name"],
                    "pathway_name": spec["pathway_name"],
                    "annotation": spec["annotation"],
                    "ligand": spec["ligand"],
                    "receptor": spec["receptor"],
                    "ligand_genes": "|".join(spec["ligand_genes"]),
                    "receptor_genes": "|".join(spec["receptor_genes"]),
                    "source": src,
                    "target": tgt,
                    "n_source": counts[src],
                    "n_target": counts[tgt],
                    "ligand_mean": float(lig_mu[i]),
                    "receptor_mean": float(rec_mu[j]),
                    "ligand_prop": float(lig_pr[i]),
                    "receptor_prop": float(rec_pr[j]),
                    "detected": bool(detected),
                    "prob": prob,
                }
            )
            probs.append(prob)
    return pd.DataFrame(rows), np.asarray(probs, dtype=np.float64)


def permute_mal_high_low(expr, pos, specs, labels, groups, src_tgt, mal_mask, nboot, rng):
    """Shuffle CLDN4-high/low among malignant cells; T/NK labels stay fixed."""
    obs_means, obs_props, counts = group_trim_means(expr, pos, labels, groups)
    obs, obs_prob = score_pairs(specs, obs_means, obs_props, groups, counts, src_tgt)
    ge = np.zeros(len(obs), dtype=np.int32)

    mal_idx = np.flatnonzero(mal_mask)
    mal_labs = labels[mal_idx].copy()
    work = labels.copy()
    for b in range(nboot):
        work[mal_idx] = rng.permutation(mal_labs)
        means, props, counts_b = group_trim_means(expr, pos, work, groups)
        _, perm_prob = score_pairs(specs, means, props, groups, counts_b, src_tgt)
        ge += (perm_prob >= obs_prob).astype(np.int32)
        if (b + 1) % 25 == 0:
            print(f"    perm {b + 1}/{nboot}", flush=True)
    pval = (ge + 1) / (nboot + 1)
    obs = obs.copy()
    obs["pval"] = pval
    obs["significant"] = (obs["pval"] < 0.05) & (obs["prob"] > 0) & obs["detected"]
    return obs


def ligand_class(genes: str) -> str:
    parts = set(genes.split("|"))
    tags = []
    if parts & BARRIER_LIGANDS:
        tags.append("barrier")
    if parts & INHIB_LIGANDS:
        tags.append("inhibitory")
    if parts & RECRUIT_LIGANDS:
        tags.append("recruit")
    if parts & ATTACK_LIGANDS:
        tags.append("attack")
    return "|".join(tags) if tags else "other"


def contrast_high_low(df: pd.DataFrame, high: str, low: str, partner: str, direction: str) -> pd.DataFrame:
    if direction == "outgoing":
        a = df[(df.source == high) & (df.target == partner)].set_index("interaction_name")
        b = df[(df.source == low) & (df.target == partner)].set_index("interaction_name")
    else:
        a = df[(df.source == partner) & (df.target == high)].set_index("interaction_name")
        b = df[(df.source == partner) & (df.target == low)].set_index("interaction_name")
    common = a.index.intersection(b.index)
    if len(common) == 0:
        return pd.DataFrame()
    out = a.loc[common, ["pathway_name", "annotation", "ligand", "receptor", "ligand_genes", "receptor_genes"]].copy()
    out["direction"] = direction
    out["partner"] = partner
    out["prob_high"] = a.loc[common, "prob"].to_numpy()
    out["prob_low"] = b.loc[common, "prob"].to_numpy()
    out["pval_high"] = a.loc[common, "pval"].to_numpy()
    out["pval_low"] = b.loc[common, "pval"].to_numpy()
    out["sig_high"] = a.loc[common, "significant"].to_numpy()
    out["sig_low"] = b.loc[common, "significant"].to_numpy()
    out["delta_prob"] = out["prob_high"] - out["prob_low"]
    out["ligand_class"] = [ligand_class(x) for x in out["ligand_genes"]]
    out["sig_either"] = out["sig_high"] | out["sig_low"]
    out["sig_diff"] = out["sig_either"] & (out["delta_prob"] != 0)
    return out.reset_index()


def cold_barrier_score(contrasts: list[pd.DataFrame]) -> dict:
    n_bar_up = n_inh_up = n_rec_down = n_atk_down = 0
    n_bar_down = n_rec_up = n_sig = 0
    for c in contrasts:
        if c.empty:
            continue
        sig = c[c["sig_diff"]]
        n_sig += int(len(sig))
        for rec in sig.itertuples(index=False):
            klass = rec.ligand_class
            d = rec.delta_prob
            if rec.direction == "outgoing":
                if "barrier" in klass:
                    n_bar_up += int(d > 0)
                    n_bar_down += int(d < 0)
                if "inhibitory" in klass:
                    n_inh_up += int(d > 0)
                if "recruit" in klass:
                    n_rec_up += int(d > 0)
                    n_rec_down += int(d < 0)
            else:
                if "attack" in klass:
                    n_atk_down += int(d < 0)
    score = n_bar_up + n_inh_up + n_rec_down + n_atk_down - n_rec_up - n_bar_down
    return {
        "score": int(score),
        "n_sig_diff": int(n_sig),
        "n_barrier_up_in_high": int(n_bar_up),
        "n_barrier_down_in_high": int(n_bar_down),
        "n_inhib_up_in_high": int(n_inh_up),
        "n_recruit_up_in_high": int(n_rec_up),
        "n_recruit_down_in_high": int(n_rec_down),
        "n_attack_in_down_in_high": int(n_atk_down),
    }


def load_sample_meta(path: Path) -> pd.DataFrame:
    meta = pd.read_excel(path)
    meta = meta.dropna(subset=["Sample"]).copy()
    meta["Sample"] = meta["Sample"].astype(str)
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post-treatment", case=False)
    resp = meta["Pathologic Response"].astype(str)
    meta["response"] = np.where(resp.isin(["MPR", "pCR"]), "MPR", np.where(resp.eq("NMPR"), "NMPR", "NE"))
    return meta


def make_split(lineage, marker, sample, response, is_post, mode: str):
    """Return labels array and malignant-mask for permutation. Split is CLDN4 only."""
    n = len(lineage)
    labels = np.array(["drop"] * n, dtype=object)
    tnk = np.isin(lineage, ["T", "NK"])
    epi = lineage == "Epithelial"

    if mode == "tertile_post":
        keep_s = is_post
        mal = epi & keep_s
        vals = marker[mal]
        if vals.size < 30:
            return labels, np.zeros(n, bool), {}
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        labels[tnk & keep_s] = "TNK"
        labels[mal & (marker <= q1)] = "Mal_low"
        labels[mal & (marker >= q2)] = "Mal_high"
        info = {"q_low": float(q1), "q_high": float(q2), "split": "tertile", "samples": "post", "marker": MARKER}
    elif mode == "tertile_nmpr":
        keep_s = is_post & (response == "NMPR")
        mal = epi & keep_s
        vals = marker[mal]
        if vals.size < 30:
            return labels, np.zeros(n, bool), {}
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        labels[tnk & keep_s] = "TNK"
        labels[mal & (marker <= q1)] = "Mal_low"
        labels[mal & (marker >= q2)] = "Mal_high"
        info = {"q_low": float(q1), "q_high": float(q2), "split": "tertile", "samples": "NMPR_post", "marker": MARKER}
    elif mode == "tertile_mpr":
        keep_s = is_post & (response == "MPR")
        mal = epi & keep_s
        vals = marker[mal]
        if vals.size < 30:
            return labels, np.zeros(n, bool), {}
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        labels[tnk & keep_s] = "TNK"
        labels[mal & (marker <= q1)] = "Mal_low"
        labels[mal & (marker >= q2)] = "Mal_high"
        info = {"q_low": float(q1), "q_high": float(q2), "split": "tertile", "samples": "MPR_post", "marker": MARKER}
    elif mode == "median_post":
        keep_s = is_post
        mal = epi & keep_s
        vals = marker[mal]
        if vals.size < 30:
            return labels, np.zeros(n, bool), {}
        med = float(np.median(vals))
        labels[tnk & keep_s] = "TNK"
        labels[mal & (marker <= med)] = "Mal_low"
        labels[mal & (marker > med)] = "Mal_high"
        info = {"median": med, "split": "median", "samples": "post", "marker": MARKER}
    elif mode == "combo_mpr":
        keep_s = is_post & np.isin(response, ["MPR", "NMPR"])
        mal = epi & keep_s
        vals = marker[mal]
        if vals.size < 30:
            return labels, np.zeros(n, bool), {}
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        for resp_name in ("MPR", "NMPR"):
            rm = keep_s & (response == resp_name)
            labels[tnk & rm] = f"TNK_{resp_name}"
            labels[mal & rm & (marker <= q1)] = f"Mal_low_{resp_name}"
            labels[mal & rm & (marker >= q2)] = f"Mal_high_{resp_name}"
        info = {"q_low": float(q1), "q_high": float(q2), "split": "tertile_x_MPR", "samples": "post", "marker": MARKER}
    else:
        raise ValueError(mode)

    mal_mask = np.isin(labels, [x for x in np.unique(labels) if str(x).startswith("Mal_")])
    return labels, mal_mask, info


def src_tgt_for(mode: str) -> list[tuple[str, str]]:
    if mode != "combo_mpr":
        return [
            ("Mal_high", "TNK"),
            ("Mal_low", "TNK"),
            ("TNK", "Mal_high"),
            ("TNK", "Mal_low"),
        ]
    pairs = []
    for resp in ("NMPR", "MPR"):
        pairs.extend(
            [
                (f"Mal_high_{resp}", f"TNK_{resp}"),
                (f"Mal_low_{resp}", f"TNK_{resp}"),
                (f"TNK_{resp}", f"Mal_high_{resp}"),
                (f"TNK_{resp}", f"Mal_low_{resp}"),
            ]
        )
    return pairs


def n_table(labels, sample, response, lineage, marker) -> pd.DataFrame:
    rows = []
    for g in sorted(set(labels) - {"drop"}):
        idx = labels == g
        rows.append(
            {
                "group": g,
                "n_cells": int(idx.sum()),
                "n_samples": int(pd.Series(sample[idx]).nunique()),
                "n_MPR_samples": int(pd.Series(sample[idx][response[idx] == "MPR"]).nunique()),
                "n_NMPR_samples": int(pd.Series(sample[idx][response[idx] == "NMPR"]).nunique()),
                "mean_CLDN4_log1p_cp10k": float(np.mean(marker[idx])) if idx.any() else None,
            }
        )
    return pd.DataFrame(rows)


def plot_counts(n_df: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.bar(n_df["group"], n_df["n_cells"], color="#4C72B0")
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=25)
    for i, r in n_df.iterrows():
        ax.text(i, r.n_cells, f"n={int(r.n_cells)}\n{int(r.n_samples)} pts", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_nsig(contrasts: dict[str, pd.DataFrame], path: Path, title: str) -> None:
    names, high, low = [], [], []
    for name, c in contrasts.items():
        names.append(name)
        if c.empty:
            high.append(0)
            low.append(0)
        else:
            high.append(int(c["sig_high"].sum()))
            low.append(int(c["sig_low"].sum()))
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.bar(x - 0.18, high, 0.36, label="CLDN4-high arm p<0.05", color="#C44E52")
    ax.bar(x + 0.18, low, 0.36, label="CLDN4-low arm p<0.05", color="#4C72B0")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel("significant LR pairs")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_top_delta(c: pd.DataFrame, path: Path, title: str, k: int = 20) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    if c is None or c.empty or not c["sig_diff"].any():
        ax.text(0.5, 0.5, "No significant LR pairs (p<0.05, detected)", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(path, dpi=140)
        plt.close(fig)
        return
    s = c[c["sig_diff"]].copy()
    s["absd"] = s["delta_prob"].abs()
    s = s.sort_values("absd", ascending=False).head(k)
    s = s.sort_values("delta_prob")
    colors = ["#C44E52" if d > 0 else "#4C72B0" for d in s["delta_prob"]]
    labels = [f"{a} ({b})" for a, b in zip(s["interaction_name"], s["direction"])]
    ax.barh(range(len(s)), s["delta_prob"], color=colors)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("prob(high) − prob(low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_per_sample_n(df: pd.DataFrame, path: Path, title: str) -> None:
    """Extra honesty figure: per-sample epithelial / TNK counts."""
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    x = np.arange(len(df))
    ax.bar(x - 0.18, df["n_epithelial"], 0.36, label="epithelial (malignant proxy)", color="#8C6D31")
    ax.bar(x + 0.18, df["n_TNK"], 0.36, label="T+NK", color="#4C72B0")
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample"], rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_extra_ligand_table(table: pd.DataFrame, path: Path, title: str, k: int = 18) -> None:
    """Extra figure: ligand table as a drawn table (kept split, sig only)."""
    fig, ax = plt.subplots(figsize=(11.2, 6.4))
    ax.set_axis_off()
    if table is None or table.empty:
        ax.text(0.5, 0.5, "No significant differential LR pairs", ha="center", va="center")
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    s = table.copy()
    s["absd"] = s["delta_prob"].abs()
    s = s.sort_values(["direction", "absd"], ascending=[True, False])
    # keep a balanced outgoing/incoming slice
    out = s[s["direction"] == "outgoing"].head(k // 2 + 2)
    inn = s[s["direction"] == "incoming"].head(max(4, k - len(out)))
    s = pd.concat([out, inn], ignore_index=True)
    cells = []
    for rec in s.itertuples(index=False):
        cells.append(
            [
                rec.direction,
                rec.interaction_name,
                rec.pathway_name,
                rec.ligand_class,
                f"{rec.prob_high:.3f}",
                f"{rec.prob_low:.3f}",
                f"{rec.delta_prob:+.3f}",
                f"{min(rec.pval_high, rec.pval_low):.3f}",
            ]
        )
    col_labels = ["dir", "pair", "pathway", "class", "P high", "P low", "ΔP", "p min"]
    tbl = ax.table(cellText=cells, colLabels=col_labels, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.0, 1.25)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#333333")
            cell.set_text_props(color="white", weight="bold")
        elif r > 0 and c == 6:
            try:
                d = float(cells[r - 1][6])
            except ValueError:
                d = 0.0
            cell.set_facecolor("#F4C2C2" if d > 0 else "#C2D4F4")
    ax.set_title(title, pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_ligand_table(contrasts: dict[str, pd.DataFrame], path: Path) -> pd.DataFrame:
    frames = []
    for name, c in contrasts.items():
        if c is None or c.empty:
            continue
        sub = c[c["sig_diff"]].copy()
        if sub.empty:
            continue
        sub["contrast"] = name
        frames.append(sub)
    if not frames:
        empty = pd.DataFrame(
            columns=[
                "interaction_name", "pathway_name", "direction", "ligand", "receptor",
                "ligand_class", "prob_high", "prob_low", "delta_prob", "pval_high",
                "pval_low", "sig_high", "sig_low", "contrast",
            ]
        )
        empty.to_csv(path, sep="\t", index=False)
        return empty
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["direction", "delta_prob"], ascending=[True, False])
    keep = [
        "interaction_name", "pathway_name", "annotation", "direction", "ligand", "receptor",
        "ligand_genes", "receptor_genes", "ligand_class", "prob_high", "prob_low",
        "delta_prob", "pval_high", "pval_low", "sig_high", "sig_low", "contrast",
    ]
    out = out[[c for c in keep if c in out.columns]]
    out.to_csv(path, sep="\t", index=False)
    return out


def run_mode(mode, expr, pos, specs, lineage, marker, sample, response, is_post, rng, out_dir: Path):
    print(f"== {mode} ==", flush=True)
    labels, mal_mask, info = make_split(lineage, marker, sample, response, is_post, mode)
    groups = sorted(set(labels) - {"drop"})
    if not groups or mal_mask.sum() < MIN_GROUP * 2:
        print(f"  skip {mode}: insufficient cells", flush=True)
        return None
    n_df = n_table(labels, sample, response, lineage, marker)
    n_df.to_csv(out_dir / f"n_cells_{mode}.tsv", sep="\t", index=False)
    plot_counts(n_df, f"GSE207422 CLDN4 {mode} cell counts", out_dir / f"fig_n_{mode}.png")
    src_tgt = src_tgt_for(mode)
    used = np.isin(labels, groups)
    pairs = permute_mal_high_low(
        expr[:, used],
        pos[:, used],
        specs,
        labels[used],
        groups,
        src_tgt,
        mal_mask[used],
        NBOOT,
        rng,
    )
    pairs.to_csv(out_dir / f"lr_pairs_{mode}.tsv", sep="\t", index=False)

    contrasts = {}
    if mode == "combo_mpr":
        for resp in ("NMPR", "MPR"):
            contrasts[f"out_{resp}"] = contrast_high_low(
                pairs, f"Mal_high_{resp}", f"Mal_low_{resp}", f"TNK_{resp}", "outgoing"
            )
            contrasts[f"in_{resp}"] = contrast_high_low(
                pairs, f"Mal_high_{resp}", f"Mal_low_{resp}", f"TNK_{resp}", "incoming"
            )
    else:
        contrasts["outgoing"] = contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "outgoing")
        contrasts["incoming"] = contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "incoming")

    for name, c in contrasts.items():
        if not c.empty:
            c = c.copy()
            c["contrast"] = name
            c.to_csv(out_dir / f"contrast_{mode}_{name}.tsv", sep="\t", index=False)
    score = cold_barrier_score(list(contrasts.values()))
    plot_nsig(contrasts, out_dir / f"fig_nsig_{mode}.png", f"Significant LR pairs — CLDN4 {mode}")
    if mode != "combo_mpr" and "outgoing" in contrasts:
        plot_top_delta(
            contrasts["outgoing"],
            out_dir / f"fig_top_outgoing_{mode}.png",
            f"Outgoing Mal→T/NK Δprob (sig only) — CLDN4 {mode}",
        )
        plot_top_delta(
            contrasts["incoming"],
            out_dir / f"fig_top_incoming_{mode}.png",
            f"Incoming T/NK→Mal Δprob (sig only) — CLDN4 {mode}",
        )
    else:
        if "out_NMPR" in contrasts:
            plot_top_delta(
                contrasts["out_NMPR"],
                out_dir / f"fig_top_outgoing_{mode}_NMPR.png",
                "Outgoing Mal→T/NK Δprob (sig only) — CLDN4 combo NMPR",
            )
            plot_top_delta(
                contrasts["in_NMPR"],
                out_dir / f"fig_top_incoming_{mode}_NMPR.png",
                "Incoming T/NK→Mal Δprob (sig only) — CLDN4 combo NMPR",
            )

    n_sig_out = int(pairs.loc[pairs.source.astype(str).str.startswith("Mal_") & pairs.significant].shape[0])
    n_sig_in = int(pairs.loc[pairs.target.astype(str).str.startswith("Mal_") & pairs.significant].shape[0])
    summary = {
        "mode": mode,
        "marker": MARKER,
        "info": info,
        "n_cells": n_df.to_dict(orient="records"),
        "n_lr_tested": int(len(pairs)),
        "n_detected": int(pairs["detected"].sum()),
        "n_significant": int(pairs["significant"].sum()),
        "n_sig_outgoing_mal_to_tnk": n_sig_out,
        "n_sig_incoming_tnk_to_mal": n_sig_in,
        "cold_barrier": score,
        "nboot": NBOOT,
        "expr_prop": EXPR_PROP,
        "trim": TRIM,
        "kh": KH,
    }
    (out_dir / f"summary_{mode}.json").write_text(json.dumps(summary, indent=2))
    print(
        f"  {mode}: tested={len(pairs)} detected={int(pairs.detected.sum())} "
        f"sig={int(pairs.significant.sum())} score={score['score']}",
        flush=True,
    )
    return summary, pairs, contrasts, labels


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT / "results")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    gene_list = []
    with gzip.open(args.matrix, "rt") as handle:
        handle.readline()
        for line in handle:
            gene_list.append(line.split("\t", 1)[0])
    matrix_genes = set(gene_list)
    print(f"matrix genes={len(matrix_genes)}", flush=True)

    lr_all = load_lr(args.db, matrix_genes)
    print(f"LR pairs with all subunits in matrix: {len(lr_all)}", flush=True)
    wanted = wanted_genes(lr_all)
    print(f"streaming {len(wanted)} genes", flush=True)
    cell_ids, found, n_umi, n_streamed = stream_matrix(args.matrix, wanted)
    n = len(cell_ids)
    if MARKER not in found:
        raise SystemExit(f"{MARKER} not in UMI matrix")

    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
    lineage = assign_lineage(scores, cd3)
    marker = log_cp[MARKER]
    sample = np.array([sample_from_barcode(bc) for bc in cell_ids])
    meta = load_sample_meta(args.meta)
    meta.to_csv(args.out / "sample_metadata.tsv", sep="\t", index=False)
    mmap = meta.set_index("Sample")
    response = np.array([mmap.loc[s, "response"] if s in mmap.index else "NA" for s in sample], dtype=object)
    is_post = np.array([bool(mmap.loc[s, "is_post"]) if s in mmap.index else False for s in sample])

    lr_genes = sorted({g for rec in lr_all.itertuples(index=False) for g in rec.ligand_genes + rec.receptor_genes})
    gene_index = {g: i for i, g in enumerate(lr_genes)}
    expr = np.vstack([log_cp[g] for g in lr_genes])
    pos = np.vstack([(found[g] > 0).astype(np.float32) for g in lr_genes])
    specs = pair_specs(lr_all, gene_index)

    lineage_counts = (
        pd.DataFrame({"sample": sample, "lineage": lineage, "response": response, "is_post": is_post})
        .groupby(["sample", "lineage", "response", "is_post"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    lineage_counts.to_csv(args.out / "lineage_counts.tsv", sep="\t", index=False)

    # Honest per-sample n (post-treatment only).
    post_rows = []
    for s in sorted(set(sample[is_post])):
        m = sample == s
        epi_n = int(((lineage == "Epithelial") & m).sum())
        tnk_n = int((np.isin(lineage, ["T", "NK"]) & m).sum())
        resp = str(response[m][0]) if m.any() else "NA"
        cldn = marker[(lineage == "Epithelial") & m]
        post_rows.append(
            {
                "sample": s,
                "response": resp,
                "n_epithelial": epi_n,
                "n_TNK": tnk_n,
                "n_total": int(m.sum()),
                "mean_CLDN4_epithelial": float(cldn.mean()) if cldn.size else None,
                "frac_CLDN4_pos_epithelial": float((cldn > 0).mean()) if cldn.size else None,
            }
        )
    per_sample = pd.DataFrame(post_rows)
    per_sample.to_csv(args.out / "per_sample_post.tsv", sep="\t", index=False)
    plot_per_sample_n(
        per_sample,
        args.out / "fig_n_per_sample.png",
        "GSE207422 post-treatment epithelial and T/NK counts (honest n)",
    )

    rng = np.random.default_rng(SEED)
    modes = ["tertile_post", "median_post", "tertile_nmpr", "tertile_mpr", "combo_mpr"]
    run_summaries = []
    kept_pack = None
    all_mode_contrasts = {}
    for mode in modes:
        got = run_mode(mode, expr, pos, specs, lineage, marker, sample, response, is_post, rng, args.out)
        if got is not None:
            run_summaries.append(got[0])
            all_mode_contrasts[mode] = got[2]

    comparable = [s for s in run_summaries if s["mode"] != "combo_mpr"]
    if comparable:
        nmpr = next((s for s in comparable if s["mode"] == "tertile_nmpr"), None)
        mpr = next((s for s in comparable if s["mode"] == "tertile_mpr"), None)
        if (
            nmpr
            and mpr
            and nmpr["cold_barrier"]["n_recruit_up_in_high"] == 0
            and mpr["cold_barrier"]["n_recruit_up_in_high"] > 0
            and nmpr["cold_barrier"]["score"] >= mpr["cold_barrier"]["score"]
        ):
            kept = nmpr
        else:
            kept = max(
                comparable,
                key=lambda s: (
                    s["cold_barrier"]["score"],
                    -s["cold_barrier"]["n_recruit_up_in_high"],
                    sum(r["n_samples"] for r in s["n_cells"] if str(r["group"]).startswith("Mal_")),
                ),
            )
    else:
        kept = None

    if kept is not None:
        kept_contrasts = all_mode_contrasts[kept["mode"]]
        ligand_tbl = write_ligand_table(kept_contrasts, args.out / "ligand_table_kept.tsv")
        # also write a compact named table for the PR
        ligand_tbl.to_csv(args.out / "ligand_table.tsv", sep="\t", index=False)
        plot_extra_ligand_table(
            ligand_tbl,
            args.out / "fig_extra_ligand_table.png",
            f"Extra figure — CLDN4 kept split ligand table ({kept['mode']}; sig ΔP only)",
        )
        # key pairs convenience file
        if not ligand_tbl.empty:
            key = ligand_tbl.copy()
            key["absd"] = key["delta_prob"].abs()
            key.sort_values("absd", ascending=False).head(40).to_csv(
                args.out / "key_pairs_kept.tsv", sep="\t", index=False
            )

    header = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023, PMID 36869384",
        "marker": MARKER,
        "additive_to": "methods/scrna_cellchat (TACSTD2 CellChat-style; taken as given)",
        "matrix_cells": n,
        "matrix_genes": n_streamed,
        "lr_pairs_in_matrix": int(len(lr_all)),
        "lineage_rule": "argmax canonical marker score on log1p(CP10k); T vs NK broken by CD3E",
        "malignant_compartment": "Epithelial lineage (author CopyKAT IDs are not public; A3 taken as given)",
        "tnk_compartment": "T or NK lineage, merged as TNK for LR tests",
        "primary_samples": "12 post-treatment resections (MPR includes pCR). Pre-treatment biopsies excluded from splits.",
        "algorithm": "CellChat-like: 10% truncated mean, Hill Kh=0.5, CellChatDB v2 protein pairs, expr_prop>=0.10, nboot=100. High/low labels permuted among malignant cells. Split gene = CLDN4 only.",
        "not_run": "CellChat R package; LIANA; TACSTD2 split; author Seurat/CopyKAT objects",
        "kept_split": None if kept is None else kept["mode"],
        "kept_reason": "Single-contrast only. After MPR split: keep NMPR tertile if it has recruit-up=0 and MPR does not; otherwise highest cold/barrier score with fewest recruit-up pairs.",
        "runs": run_summaries,
        "kept": kept,
        "per_sample_post": per_sample.to_dict(orient="records"),
    }
    (args.out / "summary.json").write_text(json.dumps(header, indent=2))
    print(
        json.dumps(
            {
                "kept_split": header["kept_split"],
                "runs": [
                    (s["mode"], s["cold_barrier"]["score"], s["n_significant"]) for s in run_summaries
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
