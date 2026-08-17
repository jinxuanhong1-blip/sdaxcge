#!/usr/bin/env python3
"""Shared CellChat-style Hill probability, lineage, and patient-delta meta."""
from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_CELLS_ARM = 10
MIN_CELLS_COMP = 20
MIN_MAL_Q4 = 40

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
NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER",
    "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1",
]
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "IFNG", "TNF", "MKI67"]


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir, matrix_genes: set[str]) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(
        columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"}
    )
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


def wanted_genes(lr: pd.DataFrame | None = None) -> set[str]:
    genes = set(EXTRA) | set(NORMAL_LUNG)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    if lr is not None:
        for rec in lr.itertuples(index=False):
            genes.update(rec.ligand_genes)
            genes.update(rec.receptor_genes)
    return genes


def ligand_class(genes) -> str:
    parts = set(genes) if not isinstance(genes, str) else set(str(genes).split("|"))
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


def trim_mean_1d(values: np.ndarray, proportiontocut: float = TRIM) -> float:
    n = int(values.size)
    if n == 0:
        return 0.0
    if n == 1:
        return float(values[0])
    k = int(n * proportiontocut)
    if k == 0:
        return float(values.mean())
    s = np.sort(values)
    return float(s[k : n - k].mean())


def geom_mean(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    if np.any(arr <= 0):
        return 0.0
    if arr.size == 1:
        return float(arr[0])
    return float(np.exp(np.mean(np.log(arr))))


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


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
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def marker_lineage(log_cp: dict[str, np.ndarray], n: int) -> np.ndarray:
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    cd3 = log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n)))
    return assign_lineage(scores, cd3)


def marker_malignant(lineage: np.ndarray, log_cp: dict[str, np.ndarray]) -> np.ndarray:
    n = lineage.size
    is_epi = lineage == "Epithelial"
    normal = module_score(log_cp, NORMAL_LUNG, n)
    if int(is_epi.sum()) >= 20:
        nl_cut = float(np.quantile(normal[is_epi], 0.75))
    else:
        nl_cut = 0.3
    return is_epi & (normal < nl_cut)


def compartment_gene_stats(log_cp, pos, idx: np.ndarray):
    means, props = {}, {}
    if idx.size == 0:
        return means, props
    for gene, vec in log_cp.items():
        means[gene] = trim_mean_1d(vec[idx])
        props[gene] = float(pos[gene][idx].mean())
    return means, props


def complex_from_maps(means: dict[str, float], props: dict[str, float], subunits: tuple[str, ...]):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    return geom_mean(vals), float(min(prs)) if prs else 0.0


def score_outgoing(lr: pd.DataFrame, log_cp, pos, mal_idx, tnk_idx) -> list[dict]:
    rows = []
    if mal_idx.size < MIN_CELLS_ARM or tnk_idx.size < MIN_CELLS_ARM:
        return rows
    mal_mu, mal_pr = compartment_gene_stats(log_cp, pos, mal_idx)
    tnk_mu, tnk_pr = compartment_gene_stats(log_cp, pos, tnk_idx)
    for rec in lr.itertuples(index=False):
        lig_mal, lig_mal_p = complex_from_maps(mal_mu, mal_pr, rec.ligand_genes)
        rec_tnk, rec_tnk_p = complex_from_maps(tnk_mu, tnk_pr, rec.receptor_genes)
        detected = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "lr_class": ligand_class(rec.ligand_genes),
                "prob": hill_prob(lig_mal, rec_tnk) if detected else 0.0,
                "detected": bool(detected),
                "ligand_mean": lig_mal,
                "receptor_mean": rec_tnk,
                "ligand_prop": lig_mal_p,
                "receptor_prop": rec_tnk_p,
            }
        )
    return rows


def highend_split(cldn4: np.ndarray, idx: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (high_idx, low_idx) among malignant barcodes. CLDN4-only.

    Zero-inflated CLDN4 collapses quantile cuts (Q1=Q3=0). High-end vs
    low-end is therefore the top vs bottom slice after a stable sort.
    """
    if idx.size == 0:
        return None
    vals = cldn4[idx]
    order = np.argsort(vals, kind="mergesort")
    n = int(idx.size)
    if mode == "median":
        if n < 2 * MIN_CELLS_ARM:
            return None
        k = n // 2
    elif mode == "tertile":
        if n < 3 * MIN_CELLS_ARM:
            return None
        k = n // 3
    elif mode == "q4q1":
        if n < MIN_MAL_Q4:
            return None
        k = max(MIN_CELLS_ARM, n // 4)
    else:
        raise ValueError(mode)
    if k < MIN_CELLS_ARM or (n - k) < MIN_CELLS_ARM:
        return None
    lo = idx[order[:k]]
    hi = idx[order[-k:]]
    if hi.size < MIN_CELLS_ARM or lo.size < MIN_CELLS_ARM:
        return None
    # refuse a no-contrast split (all ties)
    if float(vals[order[-k:]].max()) <= float(vals[order[:k]].max()) and float(vals[order[-k:]].mean()) == float(vals[order[:k]].mean()):
        return None
    return hi, lo


def assign_patient_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def dersimonian_laird(effects: np.ndarray, variances: np.ndarray) -> dict:
    """Random-effects mean of patient/unit deltas. Honest n = number of effects."""
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    mask = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[mask], v[mask]
    k = int(y.size)
    if k == 0:
        return {"k": 0, "n": 0, "mean": np.nan, "se": np.nan, "p": np.nan, "i2": np.nan, "tau2": np.nan}
    if k == 1:
        se = float(np.sqrt(v[0]))
        z = float(y[0] / se) if se > 0 else np.nan
        p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
        return {"k": 1, "n": 1, "mean": float(y[0]), "se": se, "p": p, "i2": 0.0, "tau2": 0.0}
    w = 1.0 / v
    fe = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (v + tau2)
    mean = float(np.sum(w_star * y) / np.sum(w_star))
    se = float(np.sqrt(1.0 / np.sum(w_star)))
    z = mean / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    return {"k": k, "n": k, "mean": mean, "se": se, "p": p, "i2": i2, "tau2": tau2}


def unit_delta_stats(deltas: Iterable[float]) -> dict:
    arr = np.asarray(list(deltas), dtype=float)
    arr = arr[np.isfinite(arr)]
    n = int(arr.size)
    if n < 2:
        return {
            "n": n,
            "mean": float(arr[0]) if n == 1 else np.nan,
            "sd": np.nan,
            "se": np.nan,
            "median": float(arr[0]) if n == 1 else np.nan,
            "wilcoxon_p": np.nan,
        }
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1))
    se = sd / np.sqrt(n) if n else np.nan
    try:
        p = float(stats.wilcoxon(arr, alternative="two-sided").pvalue)
    except ValueError:
        p = np.nan
    return {"n": n, "mean": mean, "sd": sd, "se": se, "median": float(np.median(arr)), "wilcoxon_p": p}


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}" if value != 0 else f"{value:.{digits}f}"
