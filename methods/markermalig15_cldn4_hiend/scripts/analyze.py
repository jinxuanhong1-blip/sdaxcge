#!/usr/bin/env python3
"""CLDN4-only CellChat-style + NicheNet-style on the n=15 marker-malignant combo.

The given combo (GSE253013 tumor n=9 + GSE291670 n=6, CLDN4 %pos vs T/NK
ρ=−0.714) is taken from existing notes and is not re-audited.

Sender split is CLDN4-only among marker-malignant cells. TACSTD2 is scored
as a companion and is never used to call high/low. No dual-high gate.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import (  # noqa: E402
    ATTACK_LIGANDS,
    BARRIER_LIGANDS,
    CYTOTOXICITY,
    EXHAUSTION,
    EXTRA,
    IFN,
    INHIB_LIGANDS,
    LINEAGE_MARKERS,
    NORMAL_LUNG,
    RECRUIT_LIGANDS,
)

KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
NBOOT = 100
SEED = 1
MIN_GROUP = 25
DETECT_FRAC = 0.10
GIVEN_RHO = -0.714
GIVEN_P = 0.0217
GIVEN_I2 = 23
GIVEN_N = 15

# GSE253013 tumors with very few marker-malignant cells in the combo notes
THIN_MALIG = {"MRC004", "MRC007"}


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_cellchat_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    lig_col = "ligand.symbol" if "ligand.symbol" in inter.columns else "ligand"
    rec_col = "receptor.symbol" if "receptor.symbol" in inter.columns else "receptor"
    rows = []
    for rec in inter.itertuples(index=False):
        d = rec._asdict()
        lig = parse_symbols(d.get(lig_col) or d.get("ligand.symbol") or d.get("ligand"))
        recp = parse_symbols(d.get(rec_col) or d.get("receptor.symbol") or d.get("receptor"))
        if not lig or not recp:
            continue
        if any(g not in matrix_genes for g in lig + recp):
            continue
        rows.append(
            {
                "interaction_name": d["interaction_name"],
                "pathway_name": d["pathway_name"],
                "annotation": d["annotation"],
                "ligand": d["ligand"],
                "receptor": d["receptor"],
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


def module(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
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
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
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


def group_trim_means(expr, pos, labels, groups):
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

    rows, probs = [], []
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
    obs_means, obs_props, counts = group_trim_means(expr, pos, labels, groups)
    obs, obs_prob = score_pairs(specs, obs_means, obs_props, groups, counts, src_tgt)
    if obs.empty:
        return obs
    ge = np.zeros(len(obs), dtype=np.int32)
    mal_idx = np.flatnonzero(mal_mask)
    mal_labs = labels[mal_idx].copy()
    work = labels.copy()
    for b in range(nboot):
        work[mal_idx] = rng.permutation(mal_labs)
        means, props, counts_b = group_trim_means(expr, pos, work, groups)
        _, perm_prob = score_pairs(specs, means, props, groups, counts_b, src_tgt)
        if perm_prob.size == obs_prob.size:
            ge += (perm_prob >= obs_prob).astype(np.int32)
        if (b + 1) % 25 == 0:
            print(f"    perm {b + 1}/{nboot}", flush=True)
    obs = obs.copy()
    obs["pval"] = (ge + 1) / (nboot + 1)
    obs["significant"] = (obs["pval"] < 0.05) & (obs["prob"] > 0) & obs["detected"]
    return obs


def ligand_class(genes: str) -> str:
    parts = set(str(genes).split("|"))
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
    if df.empty:
        return pd.DataFrame()
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


def gene_index_10x(genes: np.ndarray) -> dict[str, int]:
    idx: dict[str, int] = {}
    for i, g in enumerate(genes):
        if g not in idx:
            idx[g] = i
    return idx


def read_10x(prefix: Path):
    mtx = scipy.io.mmread(str(prefix) + "_matrix.mtx.gz").tocsr().astype(np.float32)
    with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh]
    genes = []
    with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
        mtx = mtx.T.tocsr()
    if mtx.shape[0] != len(genes) or mtx.shape[1] != len(barcodes):
        raise ValueError(f"{prefix}: mtx {mtx.shape} genes {len(genes)} barcodes {len(barcodes)}")
    return mtx, np.array(genes), np.array(barcodes)


def extract_dense(mtx, gidx: dict[str, int], names: list[str]) -> dict[str, np.ndarray]:
    out = {}
    for g in names:
        if g in gidx:
            out[g] = np.asarray(mtx[gidx[g], :].todense()).ravel()
    return out


def wanted_panel(panel_path: Path) -> list[str]:
    genes = [ln.strip() for ln in panel_path.read_text().splitlines() if ln.strip()]
    extra = list(EXTRA) + list(NORMAL_LUNG) + list(CYTOTOXICITY) + list(IFN) + list(EXHAUSTION)
    for vs in LINEAGE_MARKERS.values():
        extra.extend(vs)
    seen, out = set(), []
    for g in extra + genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def _cells_from_logcp(
    log_cp: dict[str, np.ndarray],
    keep_idx: np.ndarray,
    meta: dict,
) -> pd.DataFrame:
    genes = list(log_cp)
    mat = np.column_stack([log_cp[g][keep_idx] for g in genes]).astype(np.float32, copy=False)
    base = pd.DataFrame(meta)
    gene_df = pd.DataFrame(mat, columns=[f"g_{g}" for g in genes])
    return pd.concat([base.reset_index(drop=True), gene_df], axis=1)


def process_gse291670(datadir: Path, panel: list[str], samples: pd.DataFrame) -> pd.DataFrame:
    keep = set(panel)
    parts = []
    for rec in samples[samples.cohort == "GSE291670"].itertuples(index=False):
        prefix = datadir / rec.sample
        print(f"GSE291670 {rec.sample}", flush=True)
        mtx, genes, barcodes = read_10x(prefix)
        gidx = gene_index_10x(genes)
        n_umi = np.asarray(mtx.sum(axis=0)).ravel()
        n_genes = np.asarray((mtx > 0).sum(axis=0)).ravel()
        mt_idx = [i for i, g in enumerate(genes) if str(g).startswith("MT-")]
        mt_umi = np.asarray(mtx[mt_idx, :].sum(axis=0)).ravel() if mt_idx else np.zeros(mtx.shape[1])
        pct_mt = np.where(n_umi > 0, 100.0 * mt_umi / n_umi, 0.0)
        qc = (n_genes >= 200) & (pct_mt < 20) & (n_genes < 8000)
        mtx = mtx[:, qc]
        barcodes = barcodes[qc]
        n_umi = n_umi[qc]
        n = mtx.shape[1]
        expr = extract_dense(mtx, gidx, [g for g in keep if g in gidx])
        scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
        log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
        scores = {name: module(log_cp, gs, n) for name, gs in LINEAGE_MARKERS.items()}
        lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n))))
        normal = module(log_cp, NORMAL_LUNG, n)
        is_epi = lineage == "Epithelial"
        nl_cut = float(np.quantile(normal[is_epi], 0.75)) if is_epi.sum() >= 20 else 0.3
        malig = is_epi & (normal <= nl_cut)
        tnk = np.isin(lineage, ["T", "NK"])
        keep_idx = np.flatnonzero(malig | tnk)
        print(f"  keep malig+T/NK {keep_idx.size}/{n}", flush=True)
        meta = {
            "cohort": ["GSE291670"] * keep_idx.size,
            "sample": [rec.sample] * keep_idx.size,
            "patient": [rec.patient] * keep_idx.size,
            "response": [rec.response] * keep_idx.size,
            "barcode": barcodes[keep_idx],
            "lineage": lineage[keep_idx],
            "marker_malignant": malig[keep_idx],
            "tnk": tnk[keep_idx],
            "n_umi": n_umi[keep_idx],
            "CLDN4": log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))[keep_idx],
            "TACSTD2": log_cp.get("TACSTD2", np.zeros(n, dtype=np.float32))[keep_idx],
        }
        parts.append(_cells_from_logcp(log_cp, keep_idx, meta))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def process_gse253013(extracted: Path, samples: pd.DataFrame) -> pd.DataFrame:
    packed = np.load(extracted / "gene_panel.npz")
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(extracted / "cell_metadata.tsv", sep="\t")
    n = len(meta)
    log_cp: dict[str, np.ndarray] = {}
    if "nCount_RNA" in meta.columns:
        lib = meta["nCount_RNA"].to_numpy(dtype=float)
    elif "total" in expr:
        lib = np.asarray(expr["total"], dtype=float)
    else:
        lib = np.sum([expr[g] for g in expr if g != "total"], axis=0)
    scale = np.where(lib > 0, 1e4 / lib, 0.0)
    for g, umi in expr.items():
        if g == "total":
            continue
        log_cp[g] = np.log1p(np.asarray(umi, dtype=float) * scale).astype(np.float32)
    # Combo-note lineage (PR #256 / #290): log1p(UMI) modules, not CP10k.
    # epithelial = EPCAM/KRT8/KRT18/KRT19; malignant-like = epi & normal<=0.05
    umi = {g: np.asarray(expr[g], dtype=np.float32) for g in expr if g != "total"}
    note_lin = {
        "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
        "T": ["CD3D", "CD3E", "CD2"],
        "NK": ["NKG7", "GNLY", "FGFBP2"],
        "B": ["CD79A", "MS4A1"],
        "myeloid": ["LYZ", "CD68", "CD14"],
        "fibroblast": ["COL1A1", "DCN"],
        "endothelial": ["VWF", "PECAM1"],
    }
    names = list(note_lin)
    sc = np.vstack([module({g: np.log1p(umi[g]) for g in umi}, note_lin[k], n) for k in names])
    best = sc.argmax(axis=0)
    best_val = sc.max(axis=0)
    second = np.partition(sc, -2, axis=0)[-2]
    lineage = np.array(names, dtype=object)[best]
    lineage[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    lineage = np.where(lineage == "epithelial", "Epithelial", lineage)
    lineage = np.where(lineage == "T", "T", lineage)
    lineage = np.where(lineage == "NK", "NK", lineage)
    normal = module({g: np.log1p(umi[g]) for g in umi}, [g for g in ("SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3") if g in umi], n)
    malig = (lineage == "Epithelial") & (normal <= 0.05)
    tnk = np.isin(lineage, ["T", "NK"])
    if "patient" not in meta.columns:
        raise SystemExit("GSE253013 metadata missing patient")
    if "tissue" not in meta.columns:
        raise SystemExit("GSE253013 metadata missing tissue")
    keep_patients = set(samples.loc[samples.cohort == "GSE253013", "patient"])
    tumor = meta["tissue"].astype(str).str.upper().isin(["TUMOR", "T"])
    use = tumor & meta["patient"].isin(keep_patients) & (malig | tnk)
    idx = np.flatnonzero(np.asarray(use))
    print(f"GSE253013 tumor malig+T/NK cells in n=15: {idx.size}", flush=True)
    meta_d = {
        "cohort": ["GSE253013"] * idx.size,
        "sample": [f"{p}_Tumor" for p in meta["patient"].iloc[idx].astype(str)],
        "patient": meta["patient"].iloc[idx].astype(str).to_numpy(),
        "response": [""] * idx.size,
        "barcode": idx.astype(str),
        "lineage": lineage[idx],
        "marker_malignant": malig[idx],
        "tnk": tnk[idx],
        "n_umi": lib[idx],
        "CLDN4": log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))[idx],
        "TACSTD2": log_cp.get("TACSTD2", np.zeros(n, dtype=np.float32))[idx],
    }
    return _cells_from_logcp(log_cp, idx, meta_d)


def gene_cols(df: pd.DataFrame) -> list[str]:
    return [c[2:] for c in df.columns if c.startswith("g_")]


def split_cldn4_only(cells: pd.DataFrame) -> pd.DataFrame:
    """CLDN4-only high/low among marker-malignant cells, within cohort.

    Median log1p(CP10k) when the median is > 0. If CLDN4 is zero-inflated
    (median 0, as in GSE291670 %pos scoring), high = detected (CLDN4 > 0).
    TACSTD2 is never used to call the split.
    """
    out = cells.copy()
    out["cldn4_high"] = False
    out["cldn4_low"] = False
    out["dual_high_companion"] = False
    medians: dict[str, float] = {}
    for cohort, gidx in out.groupby("cohort").groups.items():
        mal = out.loc[gidx, "marker_malignant"].to_numpy()
        if mal.sum() == 0:
            continue
        cld = out.loc[gidx, "CLDN4"].to_numpy()
        tac = out.loc[gidx, "TACSTD2"].to_numpy() if "TACSTD2" in out.columns else np.zeros(len(gidx))
        med_c = float(np.median(cld[mal]))
        med_t = float(np.median(tac[mal]))
        if med_c <= 0:
            high = mal & (cld > 0)
            low = mal & (cld <= 0)
        else:
            high = mal & (cld >= med_c)
            low = mal & (cld < med_c)
        dual = high & (tac >= med_t) & (med_t > 0)
        if med_t <= 0:
            dual = high & (tac > 0)
        out.loc[gidx, "cldn4_high"] = high
        out.loc[gidx, "cldn4_low"] = low
        out.loc[gidx, "dual_high_companion"] = dual
        medians[str(cohort)] = med_c
    out.attrs["cldn4_median"] = medians
    return out


def cellchat_on(cells: pd.DataFrame, lr: pd.DataFrame, tag: str, rng: np.random.Generator):
    genes = gene_cols(cells)
    present = set(genes)
    lr_use = lr[lr["ligand_genes"].map(lambda xs: all(g in present for g in xs)) & lr["receptor_genes"].map(lambda xs: all(g in present for g in xs))].copy()
    if lr_use.empty:
        return pd.DataFrame(), pd.DataFrame(), {}
    gidx = {g: i for i, g in enumerate(genes)}
    expr = np.vstack([cells[f"g_{g}"].to_numpy(dtype=np.float32) for g in genes])
    pos = expr > 0
    labels = np.full(len(cells), "drop", dtype=object)
    labels[cells["cldn4_high"].to_numpy()] = "Mal_high"
    labels[cells["cldn4_low"].to_numpy()] = "Mal_low"
    labels[cells["tnk"].to_numpy()] = "TNK"
    groups = ["Mal_high", "Mal_low", "TNK"]
    src_tgt = [("Mal_high", "TNK"), ("Mal_low", "TNK"), ("TNK", "Mal_high"), ("TNK", "Mal_low")]
    specs = pair_specs(lr_use, gidx)
    mal_mask = np.isin(labels, ["Mal_high", "Mal_low"])
    print(f"CellChat {tag}: high={int((labels=='Mal_high').sum())} low={int((labels=='Mal_low').sum())} tnk={int((labels=='TNK').sum())} pairs={len(specs)}", flush=True)
    scored = permute_mal_high_low(expr, pos, specs, labels, groups, src_tgt, mal_mask, NBOOT, rng)
    out_c = contrast_high_low(scored, "Mal_high", "Mal_low", "TNK", "outgoing")
    in_c = contrast_high_low(scored, "Mal_high", "Mal_low", "TNK", "incoming")
    contrast = pd.concat([out_c, in_c], ignore_index=True) if not out_c.empty or not in_c.empty else pd.DataFrame()
    if not contrast.empty:
        contrast.insert(0, "split", tag)
    ninfo = {
        "tag": tag,
        "n_mal_high": int((labels == "Mal_high").sum()),
        "n_mal_low": int((labels == "Mal_low").sum()),
        "n_tnk": int((labels == "TNK").sum()),
        "n_pairs_scored": int(len(scored)) if not scored.empty else 0,
        "n_detected": int(scored["detected"].sum()) if not scored.empty else 0,
        "n_sig_diff": int(contrast["sig_diff"].sum()) if not contrast.empty else 0,
        "cldn4_median": cells.attrs.get("cldn4_median", {}),
    }
    return scored, contrast, ninfo


def ligand_activity(lt: pd.DataFrame, geneset: set[str], background: list[str], ligands: list[str]) -> pd.DataFrame:
    genes = [g for g in background if g in lt.index]
    if len(genes) < 20 or not ligands:
        return pd.DataFrame()
    y = np.array([1 if g in geneset else 0 for g in genes], dtype=int)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    X = lt.loc[genes, use].to_numpy(dtype=float)
    rows = []
    for j, L in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s, y)
        try:
            auroc = float(roc_auc_score(y, s))
        except ValueError:
            auroc = np.nan
        try:
            aupr = float(average_precision_score(y, s))
        except ValueError:
            aupr = np.nan
        rows.append(
            {
                "ligand": L,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "auroc": auroc,
                "aupr": aupr,
                "n_background": len(genes),
                "n_geneset": int(y.sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def expressed_in(sub: pd.DataFrame, genes: list[str], frac: float = DETECT_FRAC) -> list[str]:
    keep = []
    if sub.empty:
        return keep
    for g in genes:
        col = f"g_{g}"
        if col not in sub.columns:
            continue
        if float((sub[col] > 0).mean()) >= frac:
            keep.append(g)
    return keep


def patient_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cohort, patient), g in cells.groupby(["cohort", "patient"]):
        mal = g[g["marker_malignant"]]
        high = mal[mal["cldn4_high"]]
        low = mal[mal["cldn4_low"]]
        tnk = g[g["tnk"]]
        rec = {
            "cohort": cohort,
            "patient": patient,
            "sample": g["sample"].iloc[0],
            "response": g["response"].iloc[0] if "response" in g else "",
            "n_cells": int(len(g)),
            "n_marker_malignant": int(len(mal)),
            "n_cldn4_high": int(len(high)),
            "n_cldn4_low": int(len(low)),
            "n_tnk": int(len(tnk)),
            "tnk_fraction": float(g["tnk"].mean()),
            "mal_CLDN4_mean": float(mal["CLDN4"].mean()) if len(mal) else np.nan,
            "mal_CLDN4_pct_pos": float((mal["CLDN4"] > 0).mean()) if len(mal) else np.nan,
            "mal_TACSTD2_companion": float(mal["TACSTD2"].mean()) if len(mal) else np.nan,
            "n_dual_high_companion": int(mal["dual_high_companion"].sum()) if len(mal) else 0,
            "thin_malig_note": patient in THIN_MALIG,
        }
        if len(tnk):
            rec["tnk_cyto"] = float(tnk[[f"g_{x}" for x in CYTOTOXICITY if f"g_{x}" in tnk.columns]].mean().mean()) if any(f"g_{x}" in tnk.columns for x in CYTOTOXICITY) else np.nan
            rec["tnk_ifn"] = float(tnk[[f"g_{x}" for x in IFN if f"g_{x}" in tnk.columns]].mean().mean()) if any(f"g_{x}" in tnk.columns for x in IFN) else np.nan
            rec["tnk_exh"] = float(tnk[[f"g_{x}" for x in EXHAUSTION if f"g_{x}" in tnk.columns]].mean().mean()) if any(f"g_{x}" in tnk.columns for x in EXHAUSTION) else np.nan
        else:
            rec["tnk_cyto"] = rec["tnk_ifn"] = rec["tnk_exh"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["cohort", "patient"])


def plot_given_combo(combo: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    colors = {"GSE253013": "#4C72B0", "GSE291670": "#C44E52"}
    for cohort, g in combo.groupby("cohort"):
        ax.scatter(
            g["cldn4_pct_pos"],
            g["tnk_fraction"],
            s=55,
            c=colors.get(cohort, "gray"),
            label=f"{cohort} n={len(g)}",
            edgecolors="k",
            linewidths=0.4,
        )
        for rec in g.itertuples(index=False):
            ax.annotate(rec.patient, (rec.cldn4_pct_pos, rec.tnk_fraction), fontsize=7, alpha=0.8)
    ax.set_xlabel("marker-malignant CLDN4 %pos (combo notes)")
    ax.set_ylabel("T/NK fraction (combo notes)")
    ax.set_title(f"Given combo n={GIVEN_N}  ρ={GIVEN_RHO} (not re-audited)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_n_honesty(pat: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    axes[0].bar(pat["patient"], pat["n_marker_malignant"], color="#4C72B0")
    axes[0].axhline(MIN_GROUP, color="k", ls="--", lw=0.8, label=f"CellChat min {MIN_GROUP}")
    axes[0].set_ylabel("marker-malignant cells")
    axes[0].set_title(f"Honest n={len(pat)} samples (do not inflate)")
    axes[0].tick_params(axis="x", rotation=60, labelsize=7)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].bar(pat["patient"], pat["n_tnk"], color="#55A868")
    axes[1].axhline(MIN_GROUP, color="k", ls="--", lw=0.8)
    axes[1].set_ylabel("T/NK cells")
    axes[1].set_title("T/NK per sample")
    axes[1].tick_params(axis="x", rotation=60, labelsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_lr_table(contrast: pd.DataFrame, path: Path, title: str, top: int = 16) -> None:
    if contrast.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "no differential LR pairs", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    sig = contrast[contrast["sig_diff"]].copy()
    if sig.empty:
        sig = contrast.copy()
    sig = sig.reindex(sig["delta_prob"].abs().sort_values(ascending=False).index).head(top)
    fig, ax = plt.subplots(figsize=(8.2, max(3.2, 0.32 * len(sig) + 1.2)))
    y = np.arange(len(sig))
    colors = ["#C44E52" if d > 0 else "#4C72B0" for d in sig["delta_prob"]]
    ax.barh(y, sig["delta_prob"], color=colors)
    labels = [f"{a}–{b} ({c})" for a, b, c in zip(sig["ligand"], sig["receptor"], sig["direction"])]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("ΔP (CLDN4-high − CLDN4-low)")
    ax.set_title(title)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_nichenet(act: pd.DataFrame, path: Path, title: str, top: int = 15) -> None:
    if act.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "no NicheNet ligand activity", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    topd = act.head(top)
    fig, ax = plt.subplots(figsize=(6.6, max(3.0, 0.32 * len(topd) + 1.0)))
    y = np.arange(len(topd))
    ax.barh(y, topd["pearson"], color="#8172B3")
    ax.set_yticks(y)
    ax.set_yticklabels(topd["ligand"], fontsize=8)
    ax.set_xlabel("NicheNet ligand activity (Pearson vs target set)")
    ax.set_title(title)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def md_table(df: pd.DataFrame, cols: list[str], fmt: dict | None = None, max_rows: int = 20) -> str:
    fmt = fmt or {}
    use = df.loc[:, [c for c in cols if c in df.columns]].head(max_rows)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" if not c.startswith("n_") and c not in {"prob_high", "prob_low", "delta_prob", "pearson", "auroc"} else "---:" for c in cols) + " |"
    lines = [header, sep]
    for rec in use.itertuples(index=False):
        vals = []
        for c, v in zip(use.columns, rec):
            if c in fmt and isinstance(v, float) and np.isfinite(v):
                vals.append(fmt[c].format(v))
            elif isinstance(v, float) and np.isfinite(v):
                vals.append(f"{v:.3g}")
            else:
                vals.append("" if pd.isna(v) else str(v))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append(f"| … | {len(df) - max_rows} more rows in TSV |")
    return "\n".join(lines)


def write_finding(out: Path, pat: pd.DataFrame, contrast: pd.DataFrame, activity: pd.DataFrame, summary: dict) -> None:
    sig = contrast[contrast["sig_diff"]].copy() if not contrast.empty else pd.DataFrame()
    if not sig.empty:
        sig = sig.reindex(sig["delta_prob"].abs().sort_values(ascending=False).index)
    outg = sig[sig["direction"] == "outgoing"] if not sig.empty else pd.DataFrame()
    top_act = activity.head(12) if not activity.empty else pd.DataFrame()
    n_thin = int(pat["thin_malig_note"].sum()) if "thin_malig_note" in pat.columns else 0
    text = f"""# FINDING — CLDN4-only high-end on the marker-malignant combo (n=15)

**Additive. CLDN4 only. No dual-high.** The marker-malignant combo that already differs is taken as given from PR #290 / #279 and is **not re-audited**:

| combo | k | N | CLDN4 ρ (p, I²) |
|---|---:|---:|---|
| Marker-malignant %pos vs T/NK: GSE253013 + GSE291670 | 2 | **15** | **−0.714 (0.0217, 23%)** |

**n=15 is small. Do not inflate it.** Adjacent-normal GSE253013 libraries are not added. Extra GEO samples are not added. Cell-level p-values are not a 15-patient mixed model.

Sender definition is **CLDN4-only** among marker-malignant cells (median `log1p(CP10k)`). TACSTD2 is a companion column. Dual-high (TACSTD2-high ∩ CLDN4-high) is **not** the sender set.

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| Given combo (not re-audited) | **15** | 9 GSE253013 tumors + 6 GSE291670 |
| GSE253013 tumor patients | **9** | MRC001–004, MRC006–010 (no MRC005 in the notes) |
| GSE291670 tumors | **6** | MPR-1/2/3 vs Non-MPR-1/2/3 |
| Thin marker-malignant tumors in the notes | **{n_thin}** | MRC004 n_mal=16, MRC007 n_mal=18 (combo notes) |
| Marker-malignant cells used here | **{summary.get("n_marker_malignant", "?")}** | CLDN4-high {summary.get("n_cldn4_high", "?")} / low {summary.get("n_cldn4_low", "?")} |
| T/NK cells used here | **{summary.get("n_tnk", "?")}** | receiver |
| Dual-high companion cells (not used as senders) | {summary.get("n_dual_high_companion", "?")} | TACSTD2 scored only |
| CellChatDB v2 pairs in the extracted panel | {summary.get("n_lr_in_panel", "?")} | protein pairs |
| Significant differential directed tests | **{summary.get("n_sig_diff", 0)}** | 100 high/low permutations; smallest p = 1/101 = 0.0099 |

Patient/sample table: `results/patient_table.tsv`. Combo membership: `data/combo/samples_n15.tsv`.

### Ligand–receptor table (CellChat-style, kept median, Mal CLDN4-high vs low → T/NK)

Hill / mass-action probability on 10% truncated means (CellChatDB v2). The CellChat R package is not run. Full table: `results/lr_table.tsv`.

"""
    if outg.empty:
        text += "_No outgoing pair passed permutation p<0.05 with ΔP ≠ 0 on these n=15 samples._\n\n"
    else:
        text += md_table(
            outg,
            ["interaction_name", "pathway_name", "ligand_class", "prob_high", "prob_low", "delta_prob", "pval_high"],
            {"prob_high": "{:.3f}", "prob_low": "{:.3f}", "delta_prob": "{:+.3f}", "pval_high": "{:.3g}"},
            max_rows=18,
        )
        text += "\n\n"
    text += "### NicheNet-style ligand activity (CLDN4-high malignant → T/NK IFN / cytotoxicity)\n\n"
    text += (
        "Unsigned NicheNet-v2 regulatory potential. Target set is the a-priori T/NK IFN + cytotoxicity list, "
        "not a data-mined DE list. Ligands are those expressed in ≥10% of CLDN4-high marker-malignant cells. "
        f"Background is T/NK-expressed genes. Honest n remains {GIVEN_N}.\n\n"
    )
    if top_act.empty:
        text += "_NicheNet prior was not scored (missing prior or too few expressed ligands)._\n\n"
    else:
        text += md_table(
            top_act,
            ["ligand", "pearson", "pearson_p", "auroc", "aupr", "n_geneset"],
            {"pearson": "{:.3f}", "pearson_p": "{:.3g}", "auroc": "{:.3f}", "aupr": "{:.3f}"},
            max_rows=12,
        )
        text += "\n\n"
    text += f"""### What is not supported

- Re-auditing ρ=−0.714. That number is given.
- **n > 15.** ANT libraries, extra cohorts, or cell-as-N inflation are out of scope.
- Dual-high (TACSTD2-high and CLDN4-high) as the sender definition.
- A confirmatory ICI-response claim. GSE253013 has **no** public MPR/R labels. GSE291670 is 3 vs 3.
- Patient mixed-model communication. CellChat here is cell-pooled truncated means with a high/low label permutation.

### Files

- `results/lr_table.tsv` — CellChat-style differential LR table (the PR done criterion)
- `results/nichenet_ligand_activity.tsv` — NicheNet-style ligand activity
- `results/patient_table.tsv` — honest n=15 sample counts
- `results/fig_given_combo_n15.png` — given combo scatter (notes, not re-audited)
- `results/fig_n_honesty.png` — per-sample malignant / T/NK counts
- `results/fig_lr_delta.png` — extra LR ΔP figure
- `results/fig_nichenet_activity.png` — extra NicheNet figure
- `results/summary.json`

## 中文

在已经给出的 marker-malignant 组合（GSE253013 肿瘤 9 例 + GSE291670 6 例，**n=15**，CLDN4 %pos vs T/NK **ρ=−0.714**，不重审）上，做 **CLDN4-only** 高端：CellChat 风格的恶性 CLDN4-high vs low → T/NK，以及 NicheNet 风格的配体活性。**禁止 dual-high。** n=15 很小，不把 ANT 或其它队列加进来把 n 做大。完整 LR 表见 `results/lr_table.tsv`。
"""
    (out / "FINDING.md").write_text(text)
    (ROOT / "FINDING.md").write_text(text)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse291670", type=Path, default=Path("/tmp/markermalig15_data/GSE291670"))
    ap.add_argument("--gse253013-extracted", type=Path, default=Path("/tmp/markermalig15_data/GSE253013/extracted"))
    ap.add_argument("--nichenet-lt", type=Path, default=Path("/tmp/nichenet_priors/ligand_target.parquet"))
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    out = args.out
    figdir = ROOT / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    combo = pd.read_csv(ROOT / "data" / "combo" / "samples_n15.tsv", sep="\t")
    if len(combo) != GIVEN_N:
        raise SystemExit(f"combo notes must list {GIVEN_N} samples, got {len(combo)}")
    panel = wanted_panel(ROOT / "data" / "lr_gene_panel.txt")
    db_genes = set(panel)
    lr = load_cellchat_lr(ROOT / "db", db_genes)
    print(f"CellChat pairs in panel: {len(lr)}", flush=True)

    frames = []
    gse291 = process_gse291670(args.gse291670, panel, combo)
    frames.append(gse291)
    if (args.gse253013_extracted / "gene_panel.npz").exists() and (args.gse253013_extracted / "cell_metadata.tsv").exists():
        frames.append(process_gse253013(args.gse253013_extracted, combo))
    else:
        print("WARNING: GSE253013 extract missing; CellChat/NicheNet will use GSE291670 only but n stays 15 in the given combo", flush=True)

    cells = pd.concat(frames, ignore_index=True)
    # restrict to the 15 named samples
    cells = cells[cells["sample"].isin(set(combo["sample"])) | cells["patient"].isin(set(combo["patient"]))].copy()
    cells = split_cldn4_only(cells)
    print(
        f"cells={len(cells)} malig={int(cells.marker_malignant.sum())} "
        f"high={int(cells.cldn4_high.sum())} low={int(cells.cldn4_low.sum())} "
        f"tnk={int(cells.tnk.sum())} dual_companion={int(cells.dual_high_companion.sum())}",
        flush=True,
    )

    pat = patient_table(cells)
    pat.to_csv(out / "patient_table.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    contrasts = []
    scored_all = []
    ninfos = []
    for cohort, sub in cells.groupby("cohort"):
        sub = split_cldn4_only(sub)
        scored, contrast, ninfo = cellchat_on(sub, lr, cohort, rng)
        ninfo["cohort"] = cohort
        ninfos.append(ninfo)
        if not scored.empty:
            scored.insert(0, "cohort", cohort)
            scored_all.append(scored)
        if not contrast.empty:
            contrast.insert(0, "cohort", cohort)
            contrasts.append(contrast)
    # No cross-cohort cell pool: GSE253013 and GSE291670 CLDN4 scales differ
    # (%pos tens of percent vs <20%). Pooling cells would inflate n and mix platforms.

    contrast = pd.concat(contrasts, ignore_index=True) if contrasts else pd.DataFrame()
    scored = pd.concat(scored_all, ignore_index=True) if scored_all else pd.DataFrame()
    if not scored.empty:
        scored.to_csv(out / "lr_pairs_scored.tsv", sep="\t", index=False)

    # Primary LR table = per-cohort outgoing+incoming differential pairs
    if not contrast.empty:
        primary = contrast.copy()
        primary = primary.sort_values(["cohort", "sig_diff", "delta_prob"], ascending=[True, False, False])
        primary.to_csv(out / "lr_table.tsv", sep="\t", index=False)
        contrast.to_csv(out / "lr_contrasts_all.tsv", sep="\t", index=False)
    else:
        primary = pd.DataFrame()
        pd.DataFrame(columns=["interaction_name", "pathway_name", "direction", "prob_high", "prob_low", "delta_prob"]).to_csv(
            out / "lr_table.tsv", sep="\t", index=False
        )

    # NicheNet
    activity = pd.DataFrame()
    if args.nichenet_lt.exists():
        lt = pd.read_parquet(args.nichenet_lt)
        if lt.shape[0] < lt.shape[1]:
            print(f"warning: lt shape {lt.shape}", flush=True)
        high = cells[cells["cldn4_high"]]
        tnk = cells[cells["tnk"]]
        ligands = expressed_in(high, list(lt.columns) if hasattr(lt, "columns") else [])
        panel_in_lt = [g for g in gene_cols(cells) if g in lt.index]
        tnk_genes = expressed_in(tnk, panel_in_lt)
        if len(tnk_genes) < 30:
            tnk_genes = panel_in_lt
        target = set(IFN) | set(CYTOTOXICITY)
        activity = ligand_activity(lt, target, tnk_genes, ligands)
        if not activity.empty:
            activity.insert(0, "target_set", "a_priori_ifn_plus_cytotoxicity")
            activity.to_csv(out / "nichenet_ligand_activity.tsv", sep="\t", index=False)
        print(f"NicheNet ligands expressed in CLDN4-high: {len(ligands)}; activity rows={len(activity)}", flush=True)
    else:
        print("NicheNet prior parquet missing", flush=True)
        pd.DataFrame(columns=["ligand", "pearson", "pearson_p", "auroc"]).to_csv(
            out / "nichenet_ligand_activity.tsv", sep="\t", index=False
        )

    # extra patient-level ligand vs T/NK (honest n<=15)
    lig_cols = [c for c in cells.columns if c.startswith("g_") and c[2:] in set(RECRUIT_LIGANDS | INHIB_LIGANDS | BARRIER_LIGANDS)]
    lig_rows = []
    for (cohort, patient), g in cells.groupby(["cohort", "patient"]):
        high = g[g["cldn4_high"]]
        rec = {"cohort": cohort, "patient": patient, "n_cldn4_high": int(len(high))}
        for c in lig_cols:
            rec[c[2:]] = float(high[c].mean()) if len(high) else np.nan
        lig_rows.append(rec)
    lig_pat = pd.DataFrame(lig_rows)
    if not lig_pat.empty:
        lig_pat.to_csv(out / "patient_cldn4high_ligand_means.tsv", sep="\t", index=False)

    plot_given_combo(combo, out / "fig_given_combo_n15.png")
    plot_given_combo(combo, figdir / "fig_given_combo_n15.png")
    plot_n_honesty(pat, out / "fig_n_honesty.png")
    plot_n_honesty(pat, figdir / "fig_n_honesty.png")
    plot_lr_table(primary, out / "fig_lr_delta.png", "CellChat-style ΔP  CLDN4-high vs low → T/NK  (n=15 samples, per cohort)")
    plot_lr_table(primary, figdir / "fig_lr_delta.png", "CellChat-style ΔP  CLDN4-high vs low → T/NK  (n=15 samples, per cohort)")
    plot_nichenet(activity, out / "fig_nichenet_activity.png", "NicheNet-style ligand activity (IFN+cyto targets)")
    plot_nichenet(activity, figdir / "fig_nichenet_activity.png", "NicheNet-style ligand activity (IFN+cyto targets)")

    # extra: CLDN4-high vs T/NK IFN scatter at patient level
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    for cohort, g in pat.groupby("cohort"):
        ax.scatter(g["mal_CLDN4_mean"], g["tnk_ifn"], s=50, label=f"{cohort} n={len(g)}")
        for rec in g.itertuples(index=False):
            ax.annotate(rec.patient, (rec.mal_CLDN4_mean, rec.tnk_ifn), fontsize=7)
    ax.set_xlabel("marker-malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("T/NK IFN program")
    ax.set_title(f"Patient-level (n={len(pat)}; small)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig_patient_cldn4_vs_tnk_ifn.png", dpi=160)
    fig.savefig(figdir / "fig_patient_cldn4_vs_tnk_ifn.png", dpi=160)
    plt.close(fig)

    summary = {
        "given_combo": {
            "cohorts": "GSE253013+GSE291670",
            "n": GIVEN_N,
            "rho_CLDN4_pct_vs_tnk": GIVEN_RHO,
            "p": GIVEN_P,
            "I2": GIVEN_I2,
            "re_audited": False,
            "source": "PR #290 / #279 highlighted_combos marker/tnk/pct",
        },
        "cldn4_only": True,
        "dual_high_used_as_sender": False,
        "n_samples_in_cells": int(pat["patient"].nunique()),
        "n_marker_malignant": int(cells["marker_malignant"].sum()),
        "n_cldn4_high": int(cells["cldn4_high"].sum()),
        "n_cldn4_low": int(cells["cldn4_low"].sum()),
        "n_tnk": int(cells["tnk"].sum()),
        "n_dual_high_companion": int(cells["dual_high_companion"].sum()),
        "n_lr_in_panel": int(len(lr)),
        "n_sig_diff": int(primary["sig_diff"].sum()) if not primary.empty else 0,
        "cellchat_ninfo": ninfos,
        "nichenet_n_ligands_scored": int(len(activity)),
        "honest_n_note": "n=15 is small; do not inflate",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(out, pat, primary, activity, summary)
    print(json.dumps({k: summary[k] for k in summary if k != "cellchat_ninfo"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
