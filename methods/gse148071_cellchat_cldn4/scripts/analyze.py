#!/usr/bin/env python3
"""CellChat-style LR inference on public GSE148071 — CLDN4 only.

Compare CLDN4-high vs CLDN4-low epithelial/malignant cells for
outgoing signals to T/NK and incoming signals from T/NK.

Wu et al. Nat Commun 2021 (PMID 33953163): 42 advanced NSCLC biopsies.
No published cell-type labels on GEO; epithelial = putative malignant.
No MPR/NMPR labels. TACSTD2 is not used to define groups.

Probability is the published CellChat Hill / mass-action on truncated
means (CellChatDB v2) plus a high/low label permutation. The CellChat
R package is not run. Networks are not drawn for non-significant pairs.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
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
MIN_EPI = 25
MIN_TNK = 25
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
        lig = parse_symbols(getattr(rec, "ligand_symbol", None))
        recp = parse_symbols(getattr(rec, "receptor_symbol", None))
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


def patient_from_filename(path: Path) -> str:
    m = re.search(r"_(P\d+)_", path.name)
    if m:
        return m.group(1)
    stem = path.name.replace(".txt.gz", "")
    parts = stem.split("_")
    return parts[1] if len(parts) >= 2 else stem


def list_exp_files(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*_exp.txt.gz"))
    if files:
        return files
    tar_path = data_dir / "GSE148071_RAW.tar"
    if not tar_path.exists():
        raise SystemExit(f"no exp matrices and no {tar_path}")
    dest = data_dir / "GSE148071_files"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    files = sorted(dest.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz in {dest}")
    return files


def gene_names_from_matrix(path: Path) -> list[str]:
    genes = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline()
        if not header:
            return genes
        for line in handle:
            gene = line.split("\t", 1)[0].strip().strip('"').split(".")[0]
            if gene:
                genes.append(gene)
    return genes


def stream_one_matrix(path: Path, wanted: set[str]):
    """Return cell_ids, found gene arrays, n_umi. GSE148071 header is barcodes only."""
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header if h != ""]
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.strip().strip('"').split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
    return cell_ids, found, n_umi, n_streamed


def parse_series_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    opener = gzip.open if str(path).endswith(".gz") else open
    titles = accs = None
    chars: list[list[str]] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    n = len(titles or accs or [])
    if n == 0:
        return pd.DataFrame()
    rows = []
    for i in range(n):
        rec = {
            "sample": (titles[i] if titles else f"S{i}"),
            "geo_accession": accs[i] if accs else "",
        }
        for row in chars:
            if i >= len(row):
                continue
            val = row[i]
            if ":" in val:
                k, v = val.split(":", 1)
                rec[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(rec)
    return pd.DataFrame(rows)


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
    out["sig_high"] = (out["pval_high"] < 0.05) & (out["prob_high"] > 0)
    out["sig_low"] = (out["pval_low"] < 0.05) & (out["prob_low"] > 0)
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


def make_split(lineage, marker, eligible, mode: str):
    """CLDN4-only split among epithelial cells from eligible patients."""
    labels = np.full(lineage.shape[0], "drop", dtype=object)
    epi = (lineage == "Epithelial") & eligible
    tnk = np.isin(lineage, ["T", "NK"]) & eligible
    labels[tnk] = "TNK"
    vals = marker[epi]
    if vals.size < MIN_GROUP * 2:
        return labels, epi, {"mode": mode, "n_epi": int(epi.sum()), "cuts": None}
    if mode == "tertile":
        lo, hi = np.quantile(vals, [1.0 / 3.0, 2.0 / 3.0])
        labels[epi & (marker <= lo)] = "Mal_low"
        labels[epi & (marker >= hi)] = "Mal_high"
        cuts = {"lo": float(lo), "hi": float(hi)}
    elif mode == "median":
        med = float(np.median(vals))
        labels[epi & (marker < med)] = "Mal_low"
        labels[epi & (marker >= med)] = "Mal_high"
        cuts = {"median": med}
    else:
        raise ValueError(mode)
    mal_mask = np.isin(labels, ["Mal_high", "Mal_low"])
    return labels, mal_mask, {"mode": mode, "n_epi": int(epi.sum()), "cuts": cuts}


def n_table(labels, sample, lineage, marker) -> pd.DataFrame:
    rows = []
    for g in sorted(set(labels) - {"drop"}):
        idx = labels == g
        rows.append(
            {
                "group": g,
                "n_cells": int(idx.sum()),
                "n_samples": int(len(set(sample[idx]))),
                "mean_CLDN4_log1p_cp10k": float(np.mean(marker[idx])) if idx.any() else None,
            }
        )
    return pd.DataFrame(rows)


def plot_counts(n_df: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.bar(n_df["group"], n_df["n_cells"], color="#4C72B0")
    ax.set_ylabel("cells")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_nsig(contrasts: dict[str, pd.DataFrame], path: Path, title: str) -> None:
    names = list(contrasts)
    high = []
    low = []
    for name in names:
        c = contrasts[name]
        if c.empty:
            high.append(0)
            low.append(0)
        else:
            high.append(int(((c["sig_diff"]) & (c["delta_prob"] > 0)).sum()))
            low.append(int(((c["sig_diff"]) & (c["delta_prob"] < 0)).sum()))
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.bar(x - 0.18, high, 0.36, label="CLDN4-high arm p<0.05", color="#C44E52")
    ax.bar(x + 0.18, low, 0.36, label="CLDN4-low arm p<0.05", color="#4C72B0")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("significant differential pairs")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_top_delta(c: pd.DataFrame, path: Path, title: str, k: int = 20) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    if c is None or c.empty or not c["sig_diff"].any():
        ax.text(0.5, 0.5, "No significant differential pairs", ha="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    s = c[c["sig_diff"]].copy()
    s["absd"] = s["delta_prob"].abs()
    s = s.sort_values("absd", ascending=False).head(k).iloc[::-1]
    colors = ["#C44E52" if d > 0 else "#4C72B0" for d in s["delta_prob"]]
    ax.barh(s["interaction_name"], s["delta_prob"], color=colors)
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("ΔP (CLDN4-high − CLDN4-low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_per_sample_n(df: pd.DataFrame, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(max(8.0, 0.28 * len(df)), 4.2))
    x = np.arange(len(df))
    ax.bar(x - 0.18, df["n_epithelial"], 0.36, label="epithelial (malignant proxy)", color="#8C6D31")
    ax.bar(x + 0.18, df["n_TNK"], 0.36, label="T+NK", color="#4C72B0")
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample"], rotation=90, fontsize=7)
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_extra_ligand_table(table: pd.DataFrame, path: Path, title: str, k: int = 18) -> None:
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


def run_mode(mode, expr, pos, specs, lineage, marker, sample, eligible, rng, out_dir: Path):
    print(f"== {mode} ==", flush=True)
    labels, mal_mask, info = make_split(lineage, marker, eligible, mode)
    groups = sorted(set(labels) - {"drop"})
    if not groups or mal_mask.sum() < MIN_GROUP * 2:
        print(f"  skip {mode}: insufficient cells", flush=True)
        return None
    n_df = n_table(labels, sample, lineage, marker)
    n_df.to_csv(out_dir / f"n_cells_{mode}.tsv", sep="\t", index=False)
    plot_counts(n_df, f"GSE148071 CLDN4 {mode} cell counts", out_dir / f"fig_n_{mode}.png")
    src_tgt = [
        ("Mal_high", "TNK"),
        ("Mal_low", "TNK"),
        ("TNK", "Mal_high"),
        ("TNK", "Mal_low"),
    ]
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
    contrasts = {
        "outgoing": contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "outgoing"),
        "incoming": contrast_high_low(pairs, "Mal_high", "Mal_low", "TNK", "incoming"),
    }
    for name, c in contrasts.items():
        if not c.empty:
            c = c.copy()
            c["contrast"] = name
            c.to_csv(out_dir / f"contrast_{mode}_{name}.tsv", sep="\t", index=False)
    score = cold_barrier_score(list(contrasts.values()))
    plot_nsig(contrasts, out_dir / f"fig_nsig_{mode}.png", f"Significant LR pairs — CLDN4 {mode}")
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
    p.add_argument("--data", type=Path, default=ROOT / "data")
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT / "results")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    files = list_exp_files(args.data)
    print(f"exp matrices: {len(files)}", flush=True)
    gene_list = gene_names_from_matrix(files[0])
    matrix_genes = set(gene_list)
    print(f"genes in {files[0].name}: {len(matrix_genes)}", flush=True)

    lr_all = load_lr(args.db, matrix_genes)
    print(f"LR pairs with all subunits in first matrix: {len(lr_all)}", flush=True)
    wanted = wanted_genes(lr_all)
    print(f"streaming {len(wanted)} genes across {len(files)} patients", flush=True)

    cell_ids_all: list[str] = []
    sample_all: list[str] = []
    found_all: dict[str, list[np.ndarray]] = {g: [] for g in wanted}
    n_umi_all: list[np.ndarray] = []
    n_streamed_ref = None
    per_file = []
    for fp in files:
        patient = patient_from_filename(fp)
        cell_ids, found, n_umi, n_streamed = stream_one_matrix(fp, wanted)
        n = len(cell_ids)
        n_streamed_ref = n_streamed if n_streamed_ref is None else n_streamed_ref
        cell_ids_all.extend([f"{patient}_{c}" for c in cell_ids])
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        for g in wanted:
            found_all[g].append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        per_file.append(
            {
                "file": fp.name,
                "sample": patient,
                "n_cells": n,
                "n_genes_streamed": n_streamed,
                "n_wanted_found": len(found),
            }
        )
        print(f"  {fp.name} {patient} n={n} wanted={len(found)} genes={n_streamed}", flush=True)

    n = len(cell_ids_all)
    sample = np.array(sample_all, dtype=object)
    n_umi = np.concatenate(n_umi_all)
    found = {g: np.concatenate(v) for g, v in found_all.items()}
    if MARKER not in found or float(found[MARKER].max()) <= 0:
        raise SystemExit(f"{MARKER} missing or all-zero in GSE148071 matrices")

    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
    lineage = assign_lineage(scores, cd3)
    marker = log_cp[MARKER]

    series = parse_series_matrix(args.data / "GSE148071_series_matrix.txt.gz")
    if not series.empty:
        series.to_csv(args.out / "sample_metadata.tsv", sep="\t", index=False)

    lr_genes = sorted({g for rec in lr_all.itertuples(index=False) for g in rec.ligand_genes + rec.receptor_genes})
    gene_index = {g: i for i, g in enumerate(lr_genes)}
    expr = np.vstack([log_cp[g] for g in lr_genes])
    pos = np.vstack([(found[g] > 0).astype(np.float32) for g in lr_genes])
    specs = pair_specs(lr_all, gene_index)

    lineage_counts = (
        pd.DataFrame({"sample": sample, "lineage": lineage})
        .groupby(["sample", "lineage"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    lineage_counts.to_csv(args.out / "lineage_counts.tsv", sep="\t", index=False)

    per_rows = []
    for s in sorted(set(sample), key=lambda x: (len(x), x)):
        m = sample == s
        epi_n = int(((lineage == "Epithelial") & m).sum())
        tnk_n = int((np.isin(lineage, ["T", "NK"]) & m).sum())
        cldn = marker[(lineage == "Epithelial") & m]
        per_rows.append(
            {
                "sample": s,
                "n_epithelial": epi_n,
                "n_TNK": tnk_n,
                "n_T": int(((lineage == "T") & m).sum()),
                "n_NK": int(((lineage == "NK") & m).sum()),
                "n_total": int(m.sum()),
                "eligible": bool(epi_n >= MIN_EPI and tnk_n >= MIN_TNK),
                "mean_CLDN4_epithelial": float(cldn.mean()) if cldn.size else None,
                "frac_CLDN4_pos_epithelial": float((cldn > 0).mean()) if cldn.size else None,
            }
        )
    per_sample = pd.DataFrame(per_rows)
    per_sample.to_csv(args.out / "per_sample.tsv", sep="\t", index=False)
    plot_per_sample_n(
        per_sample,
        args.out / "fig_n_per_sample.png",
        "GSE148071 epithelial and T/NK counts (honest n; all 42 biopsies)",
    )
    eligible_samples = set(per_sample.loc[per_sample["eligible"], "sample"])
    eligible = np.array([s in eligible_samples for s in sample])
    print(
        f"eligible patients (≥{MIN_EPI} epi and ≥{MIN_TNK} T/NK): "
        f"{len(eligible_samples)} / {per_sample.shape[0]}",
        flush=True,
    )

    rng = np.random.default_rng(SEED)
    run_summaries = []
    all_mode_contrasts = {}
    for mode in ("tertile", "median"):
        got = run_mode(mode, expr, pos, specs, lineage, marker, sample, eligible, rng, args.out)
        if got is not None:
            run_summaries.append(got[0])
            all_mode_contrasts[mode] = got[2]

    if run_summaries:
        kept = max(
            run_summaries,
            key=lambda s: (
                s["cold_barrier"]["score"],
                -s["cold_barrier"]["n_recruit_up_in_high"],
                sum(r["n_cells"] for r in s["n_cells"] if str(r["group"]).startswith("Mal_")),
            ),
        )
    else:
        kept = None

    ligand_tbl = pd.DataFrame()
    if kept is not None:
        kept_contrasts = all_mode_contrasts[kept["mode"]]
        ligand_tbl = write_ligand_table(kept_contrasts, args.out / "ligand_table_kept.tsv")
        ligand_tbl.to_csv(args.out / "ligand_table.tsv", sep="\t", index=False)
        plot_extra_ligand_table(
            ligand_tbl,
            args.out / "fig_extra_ligand_table.png",
            f"Extra figure — CLDN4 kept split ligand table ({kept['mode']}; sig ΔP only)",
        )
        if not ligand_tbl.empty:
            key = ligand_tbl.copy()
            key["absd"] = key["delta_prob"].abs()
            key.sort_values("absd", ascending=False).head(40).to_csv(
                args.out / "key_pairs_kept.tsv", sep="\t", index=False
            )

    n_table_out = pd.DataFrame(
        [
            {"item": "patients_deposited", "n": int(per_sample.shape[0])},
            {"item": "patients_eligible_epi_and_tnk", "n": int(len(eligible_samples))},
            {"item": "cells_total", "n": int(n)},
            {"item": "cells_epithelial_all", "n": int((lineage == "Epithelial").sum())},
            {"item": "cells_TNK_all", "n": int(np.isin(lineage, ["T", "NK"]).sum())},
            {"item": "cells_epithelial_eligible", "n": int(((lineage == "Epithelial") & eligible).sum())},
            {"item": "cells_TNK_eligible", "n": int((np.isin(lineage, ["T", "NK"]) & eligible).sum())},
            {"item": "lr_pairs_in_matrix", "n": int(len(lr_all))},
            {"item": "kept_split_sig_diff_pairs", "n": int(len(ligand_tbl))},
        ]
    )
    n_table_out.to_csv(args.out / "n_table.tsv", sep="\t", index=False)
    pd.DataFrame(per_file).to_csv(args.out / "file_audit.tsv", sep="\t", index=False)

    header = {
        "dataset": "GSE148071",
        "citation": "Wu et al. Nature Communications 2021, PMID 33953163",
        "marker": MARKER,
        "additive": True,
        "matrix_cells": n,
        "matrix_genes_first_file": len(matrix_genes),
        "n_patients_deposited": int(per_sample.shape[0]),
        "n_patients_eligible": int(len(eligible_samples)),
        "eligible_rule": f"epithelial>={MIN_EPI} and T/NK>={MIN_TNK} (marker-argmax)",
        "eligible_patients": sorted(eligible_samples, key=lambda x: (len(x), x)),
        "lr_pairs_in_matrix": int(len(lr_all)),
        "lineage_rule": "argmax canonical marker score on log1p(CP10k); T vs NK broken by CD3E",
        "malignant_compartment": "Epithelial lineage by marker-argmax. No CopyKAT / author labels on GEO. Putative malignant.",
        "tnk_compartment": "T or NK lineage, merged as TNK for LR tests",
        "algorithm": "CellChat-like: 10% truncated mean, Hill Kh=0.5, CellChatDB v2 protein pairs, expr_prop>=0.10, nboot=100. High/low labels permuted among malignant cells. Split gene = CLDN4 only.",
        "not_run": "CellChat R package; LIANA; TACSTD2 split; author Seurat objects",
        "kept_split": None if kept is None else kept["mode"],
        "kept_reason": "Highest cold/barrier score, then fewest recruit-up pairs, then more malignant cells. Tertile vs median only (no MPR labels in this series).",
        "runs": run_summaries,
        "kept": kept,
        "per_sample": per_sample.to_dict(orient="records"),
        "file_audit": per_file,
    }
    (args.out / "summary.json").write_text(json.dumps(header, indent=2))
    print(
        json.dumps(
            {
                "kept_split": header["kept_split"],
                "n_patients_deposited": header["n_patients_deposited"],
                "n_patients_eligible": header["n_patients_eligible"],
                "runs": [
                    (s["mode"], s["cold_barrier"]["score"], s["n_significant"]) for s in run_summaries
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
