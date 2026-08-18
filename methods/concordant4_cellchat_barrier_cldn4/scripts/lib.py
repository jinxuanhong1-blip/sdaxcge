#!/usr/bin/env python3
"""Pre-specified thesis families + CellChat-style Hill scoring.

Thesis (already correct; this folder tests it on the concordant four):
  CLDN4-high → more barrier/inhibitory outgoing to T/NK
    (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
  CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing
    (CXCL9/10–CXCR3, CCL5, HLA–CD8; often NOT detected)

Barrier-up-in-high is ON-thesis. Do not bury it as a recruit-up skip.
"""
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

# Public-matrix aliases (GSE131907 stores NECTIN2 as PVRL2).
GENE_ALIASES = {
    "PVRL2": "NECTIN2",
    "JAM1": "F11R",
}

# Pre-specified families. Interaction names are CellChatDB v2 protein pairs.
# Direction: thesis_expect = "high>low" means ΔP = P_high − P_low should be > 0.
FAMILY_BARRIER_INHIB = [
    {
        "interaction_name": "JAM1_ITGAL_ITGB2",
        "ligand": "F11R",
        "receptor": "ITGAL_ITGB2",
        "axis": "F11R",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "NECTIN2_TIGIT",
        "ligand": "NECTIN2",
        "receptor": "TIGIT",
        "axis": "NECTIN2-TIGIT",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "CDH1_ITGAE_ITGB7",
        "ligand": "CDH1",
        "receptor": "ITGAE_ITGB7",
        "axis": "CDH1",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "CDH1_KLRG1",
        "ligand": "CDH1",
        "receptor": "KLRG1",
        "axis": "CDH1",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "LGALS9_HAVCR2",
        "ligand": "LGALS9",
        "receptor": "HAVCR2",
        "axis": "LGALS9",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "LGALS9_CD44",
        "ligand": "LGALS9",
        "receptor": "CD44",
        "axis": "LGALS9",
        "thesis_expect": "high>low",
    },
    {
        "interaction_name": "LGALS9_CD45",
        "ligand": "LGALS9",
        "receptor": "PTPRC",
        "axis": "LGALS9",
        "thesis_expect": "high>low",
    },
]

FAMILY_IFN_RECRUIT_MHCI = [
    {
        "interaction_name": "CXCL9_CXCR3",
        "ligand": "CXCL9",
        "receptor": "CXCR3",
        "axis": "CXCL9/10-CXCR3",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "CXCL10_CXCR3",
        "ligand": "CXCL10",
        "receptor": "CXCR3",
        "axis": "CXCL9/10-CXCR3",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "CCL5_CCR5",
        "ligand": "CCL5",
        "receptor": "CCR5",
        "axis": "CCL5",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "CCL5_CCR1",
        "ligand": "CCL5",
        "receptor": "CCR1",
        "axis": "CCL5",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "HLA-A_CD8A",
        "ligand": "HLA-A",
        "receptor": "CD8A",
        "axis": "HLA-CD8",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "HLA-B_CD8A",
        "ligand": "HLA-B",
        "receptor": "CD8A",
        "axis": "HLA-CD8",
        "thesis_expect": "low>high",
    },
    {
        "interaction_name": "HLA-C_CD8A",
        "ligand": "HLA-C",
        "receptor": "CD8A",
        "axis": "HLA-CD8",
        "thesis_expect": "low>high",
    },
]

FAMILIES = {
    "barrier_inhibitory": {
        "label": "barrier/inhibitory (CLDN4-high → more outgoing)",
        "thesis_expect": "high>low",
        "pairs": FAMILY_BARRIER_INHIB,
    },
    "ifn_recruit_mhci": {
        "label": "IFN / T-recruit / MHC-I (CLDN4-low / KD-like → more outgoing)",
        "thesis_expect": "low>high",
        "pairs": FAMILY_IFN_RECRUIT_MHCI,
    },
}

EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")

# PR #424 / #474 comparison rows (not re-audited).
COMPARISON_ROWS = [
    {
        "source": "PR #424",
        "combo": "GSE131907+GSE205335",
        "row": "NECTIN2–TIGIT",
        "split": "median",
        "n": 39,
        "n_note": "26/13",
        "delta": 0.094,
        "p": 3.6e-12,
        "metric": "median ΔP",
        "note": "given; not re-audited",
    },
    {
        "source": "PR #474",
        "combo": "GSE123902+GSE205335",
        "row": "FAMILY barrier",
        "split": "median",
        "n": 32,
        "n_note": "11/21",
        "delta": 0.0868,
        "p": 2.00e-08,
        "metric": "mean ΔP",
        "note": "given; not re-audited",
    },
]


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir, names: set[str], matrix_genes: set[str] | None = None) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(
        columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"}
    )
    rows = []
    for rec in inter.itertuples(index=False):
        if rec.interaction_name not in names:
            continue
        lig = parse_symbols(getattr(rec, "ligand_symbol", None))
        recp = parse_symbols(getattr(rec, "receptor_symbol", None))
        if not lig:
            lig = parse_symbols(rec.ligand)
        if not recp:
            recp = parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        if matrix_genes is not None and any(g not in matrix_genes for g in lig + recp):
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


def thesis_genes() -> set[str]:
    genes = set(EPI_MARKERS) | set(TNK_MARKERS) | {
        "CLDN4",
        "TACSTD2",
        "PTPRC",
        "CD3G",
        "CD4",
        "NCAM1",
        "PVRL2",
        "JAM1",
    }
    for fam in FAMILIES.values():
        for rec in fam["pairs"]:
            genes.add(rec["ligand"])
            genes.update(parse_symbols(rec["receptor"]))
    genes.update(
        {
            "F11R",
            "ITGAL",
            "ITGB2",
            "NECTIN2",
            "TIGIT",
            "CDH1",
            "ITGAE",
            "ITGB7",
            "KLRG1",
            "LGALS9",
            "HAVCR2",
            "CD44",
            "CXCL9",
            "CXCL10",
            "CXCR3",
            "CCL5",
            "CCR5",
            "CCR1",
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "CD8A",
            "CD8B",
        }
    )
    return genes


def apply_aliases(extracted: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Map public-matrix aliases onto CellChatDB / HGNC names."""
    out = dict(extracted)
    for src, dest in GENE_ALIASES.items():
        if dest not in out and src in out:
            out[dest] = out[src]
    return out


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
    if idx.size == 0:
        return None
    vals = cldn4[idx]
    if mode == "median":
        if idx.size < 2 * MIN_CELLS_ARM:
            return None
        cut = float(np.median(vals))
        hi = idx[vals > cut]
        lo = idx[vals <= cut]
    elif mode == "tertile":
        if idx.size < 3 * MIN_CELLS_ARM:
            return None
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        hi = idx[vals >= q2]
        lo = idx[vals <= q1]
    elif mode == "q4q1":
        if idx.size < MIN_MAL_Q4:
            return None
        q1, q3 = np.quantile(vals, [0.25, 0.75])
        hi = idx[vals >= q3]
        lo = idx[vals <= q1]
    elif mode == "pctpos":
        # Cell-level %pos gate: CLDN4>0 vs CLDN4==0 when both arms meet the floor.
        if idx.size < 2 * MIN_CELLS_ARM:
            return None
        hi = idx[vals > 0]
        lo = idx[vals <= 0]
    else:
        raise ValueError(mode)
    if hi.size < MIN_CELLS_ARM or lo.size < MIN_CELLS_ARM:
        return None
    return hi, lo


def dersimonian_laird(effects: np.ndarray, variances: np.ndarray) -> dict:
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
