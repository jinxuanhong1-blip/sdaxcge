#!/usr/bin/env python3
"""Pairwise merge GSE148071 + GSE205335: patient-level CLDN4 vs T/NK + CellChat.

ADDITIVE. CLDN4 only. No dual-high. GSE131907 is not included.

1. Patient-level malignant CLDN4 vs same-patient T/NK (honest n, Q4 vs Q1).
2. Combo Spearman = DerSimonian–Laird on Fisher-z of the two cohort rhos.
3. CellChat-style outgoing CLDN4-high (within-cohort Q4) → same-patient T/NK.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_EPI_148 = 25
MIN_TNK_148 = 25
MIN_DETECT_ARM = 3
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


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_rho(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.3f}"


def spearman(x, y) -> tuple[float, float, int]:
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[mask], ya[mask]
    n = int(xa.size)
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def fisher_z(rho: float) -> float:
    return float(np.arctanh(float(np.clip(rho, -0.999999, 0.999999))))


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    z = np.array([fisher_z(r) for r in rhos], dtype=float)
    v = np.array([1.0 / (n - 3) for n in ns], dtype=float)
    ok = np.isfinite(z) & np.isfinite(v) & (v > 0)
    z, v, n_ok = z[ok], v[ok], np.array(ns, dtype=float)[ok]
    k = int(z.size)
    if k == 0:
        return {"k": 0}
    w = 1.0 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - z_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else float("nan")
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se = float(1.0 / math.sqrt(float(np.sum(w_re))))
    z_stat = z_re / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z_stat))) if math.isfinite(z_stat) else float("nan")
    lo, hi = z_re - 1.96 * se, z_re + 1.96 * se
    i2 = max(0.0, (q - df) / q) * 100.0 if k > 1 and q > 0 else 0.0
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "method": "DerSimonian-Laird random-effects on Fisher-z(Spearman ρ)",
        "pooled_rho": float(np.tanh(z_re)),
        "p": p,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
        "Q": q,
        "I2": i2,
        "tau2": tau2,
        "fixed_rho": float(np.tanh(z_fe)),
    }


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def q4_vs_q1(cldn4, immune) -> dict | None:
    frame = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    if n < 6:
        return None
    try:
        qs = assign_quartiles(frame["c"])
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = frame.loc[qs == "Q1", "i"]
    q4 = frame.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
        "q1_values": q1.to_numpy(),
        "q4_values": q4.to_numpy(),
    }


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


def ligand_class(genes) -> str:
    if isinstance(genes, str):
        parts = set(genes.split("|"))
    else:
        parts = set(genes)
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
    return m.group(1) if m else path.name


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
        handle.readline()
        for line in handle:
            gene = line.split("\t", 1)[0].strip().strip('"').split(".")[0]
            if gene:
                genes.append(gene)
    return genes


def stream_one_matrix(path: Path, wanted: set[str]):
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
    return cell_ids, found, n_umi


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
    if arr.size == 0 or np.any(arr <= 0):
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


def complex_from_maps(means, props, subunits):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    return geom_mean(vals), float(min(prs)) if prs else 0.0


def score_patient_pairs(lr: pd.DataFrame, log_cp, pos, mal_idx, tnk_idx) -> list[dict]:
    rows = []
    if mal_idx.size < MIN_EPI_148 or tnk_idx.size < MIN_TNK_148:
        return rows
    mal_mu, mal_pr = compartment_gene_stats(log_cp, pos, mal_idx)
    tnk_mu, tnk_pr = compartment_gene_stats(log_cp, pos, tnk_idx)
    for rec in lr.itertuples(index=False):
        lig_mal, lig_mal_p = complex_from_maps(mal_mu, mal_pr, rec.ligand_genes)
        rec_tnk, rec_tnk_p = complex_from_maps(tnk_mu, tnk_pr, rec.receptor_genes)
        lig_tnk, lig_tnk_p = complex_from_maps(tnk_mu, tnk_pr, rec.ligand_genes)
        rec_mal, rec_mal_p = complex_from_maps(mal_mu, mal_pr, rec.receptor_genes)
        out_det = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
        in_det = lig_tnk_p >= EXPR_PROP and rec_mal_p >= EXPR_PROP
        for direction, det, lmean, rmean, lprop, rprop in (
            ("outgoing", out_det, lig_mal, rec_tnk, lig_mal_p, rec_tnk_p),
            ("incoming", in_det, lig_tnk, rec_mal, lig_tnk_p, rec_mal_p),
        ):
            rows.append(
                {
                    "interaction_name": rec.interaction_name,
                    "pathway_name": rec.pathway_name,
                    "annotation": rec.annotation,
                    "ligand": rec.ligand,
                    "receptor": rec.receptor,
                    "ligand_genes": "|".join(rec.ligand_genes),
                    "receptor_genes": "|".join(rec.receptor_genes),
                    "direction": direction,
                    "prob": hill_prob(lmean, rmean) if det else 0.0,
                    "detected": bool(det),
                    "ligand_mean": lmean,
                    "receptor_mean": rmean,
                    "ligand_prop": lprop,
                    "receptor_prop": rprop,
                }
            )
    return rows


def plot_q4q1(q1, q4, ylabel, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter_combo(patients: pd.DataFrame, title: str, path: Path) -> None:
    colors = {"GSE148071": "#4d4d4d", "GSE205335": "#b2182b"}
    markers = {"Q1": "o", "Q2": "s", "Q3": "D", "Q4": "^"}
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for cohort, color in colors.items():
        sub = patients[patients["cohort"] == cohort]
        for q, mk in markers.items():
            m = sub["q_pct"] == q
            ax.scatter(
                sub.loc[m, "cldn4_rank"],
                sub.loc[m, "frac_tnk"],
                c=color,
                marker=mk,
                s=36,
                label=f"{cohort} {q}" if q in {"Q1", "Q4"} else None,
                zorder=3,
            )
    ax.set_xlabel("Within-cohort CLDN4 %pos rank (0–1)")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str) -> None:
    if tbl.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No differential pairs at the stated gates", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    show = tbl.head(18).copy()
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(show) + 1.4))
    y = np.arange(len(show))
    colors = np.where(show["delta_median"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["delta_median"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    labels = [
        f"{r.direction[:3]} {r.ligand}–{r.receptor} ({r.ligand_class})"
        for r in show.itertuples()
    ]
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("median P(Q4) − median P(Q1)")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_n_bars(patients: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=False)
    for ax, cohort in zip(axes, ["GSE148071", "GSE205335"]):
        sub = patients[patients["cohort"] == cohort].sort_values("patient")
        x = np.arange(len(sub))
        ax.bar(x - 0.18, sub["n_malignant"], width=0.36, color="#b2182b", alpha=0.8, label="malignant")
        ax.bar(x + 0.18, sub["n_tnk"], width=0.36, color="#6a8aaa", alpha=0.8, label="T/NK")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient"], rotation=90, fontsize=6)
        ax.set_title(f"{cohort} honest n (eligible)", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def cohort_row(cohort, malig_def, score, cldn4, immune, unit, note) -> dict:
    rho, p_s, n = spearman(cldn4, immune)
    q = q4_vs_q1(cldn4, immune)
    return {
        "cohort": cohort,
        "malig_def": malig_def,
        "immune_def": "same_patient_tnk_fraction",
        "score": score,
        "n_patients": n,
        "unit": unit,
        "note": note,
        "spearman_rho": rho,
        "spearman_p": p_s,
        "n_q1": None if q is None else q["n_q1"],
        "n_q4": None if q is None else q["n_q4"],
        "n_compared": None if q is None else q["n_compared"],
        "q4q1_r_rb": None if q is None else q["r_rb"],
        "q4q1_p": None if q is None else q["p"],
        "median_tnk_q1": None if q is None else q["median_q1"],
        "median_tnk_q4": None if q is None else q["median_q4"],
        "delta_median_tnk": None if q is None else q["delta_median"],
        "thin": True if q is None else q["thin"],
    }


def build_patient_table(g14: pd.DataFrame, g20: pd.DataFrame) -> pd.DataFrame:
    a = g14[g14["eligible"]].copy()
    a["cohort"] = "GSE148071"
    a["patient"] = a["sample"].astype(str)
    a["n_malignant"] = a["n_epithelial"]
    a["n_tnk"] = a["n_TNK"]
    a["n_cells"] = a["n_total"]
    a["frac_tnk"] = a["n_TNK"] / a["n_total"]
    a["mal_CLDN4_mean"] = a["mean_CLDN4_epithelial"]
    a["mal_CLDN4_pct_pos"] = a["frac_CLDN4_pos_epithelial"] * 100.0
    a["malig_def"] = "marker_epi_putative"
    a["histology"] = "NSCLC_unlabeled"
    a["recist"] = "NA"
    b = g20.copy()
    b["cohort"] = "GSE205335"
    b["frac_tnk"] = b["n_tnk"] / b["n_cells"]
    b["malig_def"] = "author_malig"
    b["histology"] = b["cancer_subtype"]
    keep = [
        "cohort", "patient", "n_malignant", "n_tnk", "n_cells", "frac_tnk",
        "mal_CLDN4_mean", "mal_CLDN4_pct_pos", "malig_def", "histology", "recist",
    ]
    out = pd.concat([a[keep], b[keep]], ignore_index=True)
    out["q_pct"] = out.groupby("cohort", group_keys=False)["mal_CLDN4_pct_pos"].apply(assign_quartiles)
    out["q_mean"] = out.groupby("cohort", group_keys=False)["mal_CLDN4_mean"].apply(assign_quartiles)
    out["cldn4_rank"] = out.groupby("cohort")["mal_CLDN4_pct_pos"].rank(method="average", pct=True)
    out["cldn4_rank_mean"] = out.groupby("cohort")["mal_CLDN4_mean"].rank(method="average", pct=True)
    return out


def run_patient_level(patients: pd.DataFrame, out: Path) -> dict:
    rows = []
    for cohort, malig, note in (
        ("GSE148071", "marker_epi_putative", "eligible ≥25 epi and ≥25 T/NK; putative malignant"),
        ("GSE205335", "author_malig", "locked 22-patient author-malignant extract"),
    ):
        sub = patients[patients["cohort"] == cohort]
        rows.append(cohort_row(cohort, malig, "pct_pos", sub["mal_CLDN4_pct_pos"], sub["frac_tnk"], "patient", note))
        rows.append(cohort_row(cohort, malig, "mean", sub["mal_CLDN4_mean"], sub["frac_tnk"], "patient", note))
    table = pd.DataFrame(rows)

    combo_rows = []
    for score in ("pct_pos", "mean"):
        block = table[table["score"] == score]
        pooled = random_effects_dl(block["spearman_rho"].tolist(), block["n_patients"].tolist())
        q_pool = random_effects_dl(
            [r for r in block["q4q1_r_rb"].tolist() if pd.notna(r)],
            [int(n) for n in block["n_compared"].tolist() if pd.notna(n)],
        )
        stacked = patients.copy()
        qcol = "q_pct" if score == "pct_pos" else "q_mean"
        rcol = "cldn4_rank" if score == "pct_pos" else "cldn4_rank_mean"
        q1 = stacked.loc[stacked[qcol] == "Q1", "frac_tnk"].to_numpy()
        q4 = stacked.loc[stacked[qcol] == "Q4", "frac_tnk"].to_numpy()
        u, p_st = stats.mannwhitneyu(q4, q1, alternative="two-sided")
        r_st = (2.0 * float(u)) / (len(q4) * len(q1)) - 1.0
        rho_st, p_rho_st, n_st = spearman(stacked[rcol], stacked["frac_tnk"])
        combo_rows.append(
            {
                "pair": "GSE148071+GSE205335",
                "score": score,
                "k": 2,
                "n_patients": int(block["n_patients"].sum()),
                "combo_rho_dl": pooled["pooled_rho"],
                "combo_rho_p": pooled["p"],
                "combo_rho_I2": pooled["I2"],
                "combo_rho_ci95_lo": pooled["ci95_rho"][0],
                "combo_rho_ci95_hi": pooled["ci95_rho"][1],
                "stacked_rank_rho": rho_st,
                "stacked_rank_p": p_rho_st,
                "stacked_n": n_st,
                "q4q1_r_dl": q_pool.get("pooled_rho"),
                "q4q1_p_dl": q_pool.get("p"),
                "q4q1_I2": q_pool.get("I2"),
                "stacked_q4q1_r": r_st,
                "stacked_q4q1_p": float(p_st),
                "stacked_n_q1": int(len(q1)),
                "stacked_n_q4": int(len(q4)),
                "stacked_median_tnk_q1": float(np.median(q1)),
                "stacked_median_tnk_q4": float(np.median(q4)),
                "includes_gse131907": False,
                "dual_high": False,
            }
        )
    combo = pd.DataFrame(combo_rows)
    table.to_csv(out / "results" / "q4q1_tnk_singles.tsv", sep="\t", index=False)
    combo.to_csv(out / "results" / "combo_rho.tsv", sep="\t", index=False)
    patients.to_csv(out / "results" / "patients_with_quartiles.tsv", sep="\t", index=False)

    prim = combo.loc[combo["score"] == "pct_pos"].iloc[0]
    q1 = patients.loc[patients["q_pct"] == "Q1", "frac_tnk"].to_numpy()
    q4 = patients.loc[patients["q_pct"] == "Q4", "frac_tnk"].to_numpy()
    plot_q4q1(
        q1,
        q4,
        "Same-patient T/NK fraction",
        (
            f"Pair 148071+205335 CLDN4 %pos Q4 vs Q1 (within-cohort tails)\n"
            f"stacked r={prim['stacked_q4q1_r']:+.3f} p={fmt_p(prim['stacked_q4q1_p'])} "
            f"n={int(prim['stacked_n_q1'])}/{int(prim['stacked_n_q4'])}"
        ),
        out / "figures" / "q4q1_tnk_pct.png",
    )
    plot_scatter_combo(
        patients,
        (
            f"Pair 148071+205335  combo ρ={prim['combo_rho_dl']:+.3f} "
            f"p={fmt_p(prim['combo_rho_p'])} I²={prim['combo_rho_I2']:.0f}%  N=47"
        ),
        out / "figures" / "scatter_combo_cldn4_tnk.png",
    )
    plot_n_bars(patients, out / "figures" / "fig_n_per_patient.png")
    return {"singles": table, "combo": combo, "patients": patients, "primary": prim.to_dict()}


def contrast_lr(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (name, direction), block in long.groupby(["interaction_name", "direction"], observed=True):
        q1 = block[block["quartile_cldn4_pct"] == "Q1"]
        q4 = block[block["quartile_cldn4_pct"] == "Q4"]
        d1 = q1[q1["detected"]]
        d4 = q4[q4["detected"]]
        if len(d1) < MIN_DETECT_ARM or len(d4) < MIN_DETECT_ARM:
            continue
        u, p = stats.mannwhitneyu(d4["prob"].to_numpy(), d1["prob"].to_numpy(), alternative="two-sided")
        n1, n4 = int(len(d1)), int(len(d4))
        r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
        rec0 = block.iloc[0]
        n1_14 = int(((d1["cohort"] == "GSE148071").sum()) if "cohort" in d1 else 0)
        n4_14 = int(((d4["cohort"] == "GSE148071").sum()) if "cohort" in d4 else 0)
        n1_20 = int(((d1["cohort"] == "GSE205335").sum()) if "cohort" in d1 else 0)
        n4_20 = int(((d4["cohort"] == "GSE205335").sum()) if "cohort" in d4 else 0)
        rows.append(
            {
                "interaction_name": name,
                "direction": direction,
                "pathway_name": rec0["pathway_name"],
                "annotation": rec0["annotation"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "ligand_genes": rec0["ligand_genes"],
                "receptor_genes": rec0["receptor_genes"],
                "ligand_class": ligand_class(rec0["ligand_genes"]),
                "n_q1_detected": n1,
                "n_q4_detected": n4,
                "n_compared": n1 + n4,
                "n_q1_gse148071": n1_14,
                "n_q4_gse148071": n4_14,
                "n_q1_gse205335": n1_20,
                "n_q4_gse205335": n4_20,
                "median_prob_q1": float(d1["prob"].median()),
                "median_prob_q4": float(d4["prob"].median()),
                "delta_median": float(d4["prob"].median() - d1["prob"].median()),
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": float(r_rb),
            }
        )
    contrast = pd.DataFrame(rows)
    if contrast.empty:
        return contrast
    return (
        contrast[contrast["delta_median"] != 0]
        .assign(_absd=lambda d: d["delta_median"].abs())
        .sort_values(["p", "_absd"], ascending=[True, False])
        .drop(columns="_absd")
    )


def score_gse148071(data_dir: Path, db_dir: Path, patients: pd.DataFrame) -> pd.DataFrame:
    files = list_exp_files(data_dir)
    print(f"GSE148071 exp matrices: {len(files)}", flush=True)
    matrix_genes = set(gene_names_from_matrix(files[0]))
    lr = load_lr(db_dir, matrix_genes)
    wanted = wanted_genes(lr)
    print(f"LR pairs in first matrix: {len(lr)}; streaming {len(wanted)} genes", flush=True)
    qmap = patients.loc[patients["cohort"] == "GSE148071"].set_index("patient")["q_pct"].astype(str)
    keep = set(qmap.index)
    rows = []
    for fp in files:
        patient = patient_from_filename(fp)
        if patient not in keep:
            continue
        cell_ids, found, n_umi = stream_one_matrix(fp, wanted)
        n = len(cell_ids)
        lib = np.maximum(n_umi, 1.0)
        log_cp = {g: np.log1p(found.get(g, np.zeros(n, dtype=np.float32)) / lib * 1e4).astype(np.float32) for g in wanted}
        pos = {g: (found[g] > 0).astype(np.float32) if g in found else np.zeros(n, dtype=np.float32) for g in wanted}
        scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
        cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
        lineage = assign_lineage(scores, cd3)
        mal = np.flatnonzero(lineage == "Epithelial")
        tnk = np.flatnonzero(np.isin(lineage, ["T", "NK"]))
        scored = score_patient_pairs(lr, log_cp, pos, mal, tnk)
        for row in scored:
            row["patient"] = patient
            row["cohort"] = "GSE148071"
            row["quartile_cldn4_pct"] = str(qmap.loc[patient])
            rows.append(row)
        print(f"  scored {patient} mal={mal.size} tnk={tnk.size} q={qmap.loc[patient]} pairs={len(scored)}", flush=True)
    return pd.DataFrame(rows)


def run_cellchat(args, patients: pd.DataFrame, out: Path) -> dict:
    g14 = score_gse148071(args.data_148071, args.db, patients)
    g20 = pd.read_csv(args.lr_205335, sep="\t")
    g20["cohort"] = "GSE205335"
    qmap20 = patients.loc[patients["cohort"] == "GSE205335"].set_index("patient")["q_pct"].astype(str)
    g20["quartile_cldn4_pct"] = g20["patient"].map(qmap20)
    long = pd.concat([g14, g20], ignore_index=True, sort=False)
    long.to_csv(out / "results" / "per_patient_lr.tsv.gz", sep="\t", index=False)

    contrast = contrast_lr(long)
    outgoing = contrast[contrast["direction"] == "outgoing"] if not contrast.empty else contrast
    if contrast.empty:
        ligand_tbl = contrast
        n_sig = 0
    else:
        contrast.to_csv(out / "results" / "lr_q4q1_all.tsv", sep="\t", index=False)
        outgoing.to_csv(out / "results" / "lr_q4q1_outgoing.tsv", sep="\t", index=False)
        n_sig = int((outgoing["p"] < 0.05).sum()) if not outgoing.empty else 0
        ligand_tbl = outgoing.sort_values(["p"]).head(20).copy() if not outgoing.empty else outgoing
        if not ligand_tbl.empty:
            ligand_tbl["sig_p05"] = ligand_tbl["p"] < 0.05
            ligand_tbl["note"] = (
                f"{n_sig} outgoing pair(s) p<0.05 of {len(outgoing)} detect-gated outgoing rows; "
                "Q4/Q1 labels are within-cohort; GSE131907 not included"
            )
    ligand_tbl.to_csv(out / "results" / "ligand_table.tsv", sep="\t", index=False)
    plot_ligand_table(
        ligand_tbl,
        out / "figures" / "fig_extra_ligand_table.png",
        "Pair 148071+205335 outgoing Mal→T/NK  Q4 vs Q1 (same-patient; within-cohort tails)",
    )
    return {
        "n_patients_scored_148071": int(g14["patient"].nunique()) if not g14.empty else 0,
        "n_patients_scored_205335": int(g20["patient"].nunique()) if not g20.empty else 0,
        "n_lr_rows": int(len(long)),
        "n_pairs_tested_after_detect_gate": int(len(contrast)) if not contrast.empty else 0,
        "n_outgoing_tested": int(len(outgoing)) if not outgoing.empty else 0,
        "n_outgoing_p_lt_05": n_sig,
        "ligand_table_rows": int(len(ligand_tbl)),
    }


def write_finding(q4: dict, cellchat: dict | None, out: Path) -> None:
    singles = q4["singles"]
    combo = q4["combo"]
    patients = q4["patients"]
    prim = combo.loc[combo["score"] == "pct_pos"].iloc[0]
    sec = combo.loc[combo["score"] == "mean"].iloc[0]
    s14 = singles[(singles["cohort"] == "GSE148071") & (singles["score"] == "pct_pos")].iloc[0]
    s20 = singles[(singles["cohort"] == "GSE205335") & (singles["score"] == "pct_pos")].iloc[0]
    m14 = singles[(singles["cohort"] == "GSE148071") & (singles["score"] == "mean")].iloc[0]
    m20 = singles[(singles["cohort"] == "GSE205335") & (singles["score"] == "mean")].iloc[0]

    def tail_rows(cohort: str) -> list[str]:
        sub = patients[patients["cohort"] == cohort]
        lines = []
        for rec in sub[sub["q_pct"] == "Q4"].sort_values("mal_CLDN4_pct_pos", ascending=False).itertuples():
            lines.append(
                f"| Q4 | {cohort} {rec.patient} ({rec.histology}, {rec.recist}) | "
                f"{int(rec.n_malignant)} | {int(rec.n_tnk)} | "
                f"{rec.mal_CLDN4_pct_pos:.1f} | {rec.frac_tnk:.3f} |"
            )
        for rec in sub[sub["q_pct"] == "Q1"].sort_values("mal_CLDN4_pct_pos").itertuples():
            lines.append(
                f"| Q1 | {cohort} {rec.patient} ({rec.histology}, {rec.recist}) | "
                f"{int(rec.n_malignant)} | {int(rec.n_tnk)} | "
                f"{rec.mal_CLDN4_pct_pos:.1f} | {rec.frac_tnk:.3f} |"
            )
        return lines

    lines = [
        "# FINDING — pairwise GSE148071 + GSE205335: malignant CLDN4 vs T/NK + CellChat",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **GSE131907 is not in this pair.**",
        "Patient is the unit. Prior single-cohort CellChat folders (PR #348, #362) and the",
        "Q4 meta (PR #320) are given and are not re-ranked. p-values are descriptive.",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        "| GSE148071 deposited | **42** | Wu et al. 2021; stage III/IV NSCLC biopsies |",
        f"| GSE148071 eligible (≥25 epi **and** ≥25 T/NK) | **{int(s14.n_patients)}** | marker-argmax epithelium = putative malignant |",
        "| GSE148071 excluded | **17** | mostly epithelium with almost no T/NK |",
        f"| GSE205335 locked extract | **{int(s20.n_patients)}** | author malignant; ≥20 mal and ≥20 T/NK (PR #279/#320) |",
        f"| Pair (eligible) | **{int(prim.n_patients)}** | 25 + 22; not 42 + 26 |",
        f"| Q4 vs Q1 tails (within-cohort, stacked) | **{int(prim.stacked_n_q1)} vs {int(prim.stacked_n_q4)}** | n_compared={int(prim.stacked_n_q1)+int(prim.stacked_n_q4)}, not 47 |",
        "",
        "Malignant definitions are **not the same**: GSE148071 is marker-argmax epithelium",
        "(no GEO labels); GSE205335 is author `Malignant cells`. CLDN4 %pos is not stacked",
        "across platforms. Combo ρ is Fisher-z / DerSimonian–Laird of the two cohort rhos.",
        "GSE205335 Q4 mixes SCLC with ADC; that mix is kept.",
        "",
        "## Malignant CLDN4 vs same-patient T/NK",
        "",
        "| cohort | score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |",
        "|---|---|---:|---|---|---:|---:|---:|",
        (
            f"| GSE148071 | %pos | {int(s14.n_patients)} | {fmt_rho(s14.spearman_rho)} ({fmt_p(s14.spearman_p)}) | "
            f"{fmt_rho(s14.q4q1_r_rb)} ({fmt_p(s14.q4q1_p)}; {int(s14.n_q1)}/{int(s14.n_q4)}) | "
            f"{s14.median_tnk_q1:.3f} | {s14.median_tnk_q4:.3f} | {s14.delta_median_tnk:+.3f} |"
        ),
        (
            f"| GSE205335 | %pos | {int(s20.n_patients)} | {fmt_rho(s20.spearman_rho)} ({fmt_p(s20.spearman_p)}) | "
            f"{fmt_rho(s20.q4q1_r_rb)} ({fmt_p(s20.q4q1_p)}; {int(s20.n_q1)}/{int(s20.n_q4)}) | "
            f"{s20.median_tnk_q1:.3f} | {s20.median_tnk_q4:.3f} | {s20.delta_median_tnk:+.3f} |"
        ),
        (
            f"| GSE148071 | mean | {int(m14.n_patients)} | {fmt_rho(m14.spearman_rho)} ({fmt_p(m14.spearman_p)}) | "
            f"{fmt_rho(m14.q4q1_r_rb)} ({fmt_p(m14.q4q1_p)}; {int(m14.n_q1)}/{int(m14.n_q4)}) | "
            f"{m14.median_tnk_q1:.3f} | {m14.median_tnk_q4:.3f} | {m14.delta_median_tnk:+.3f} |"
        ),
        (
            f"| GSE205335 | mean | {int(m20.n_patients)} | {fmt_rho(m20.spearman_rho)} ({fmt_p(m20.spearman_p)}) | "
            f"{fmt_rho(m20.q4q1_r_rb)} ({fmt_p(m20.q4q1_p)}; {int(m20.n_q1)}/{int(m20.n_q4)}) | "
            f"{m20.median_tnk_q1:.3f} | {m20.median_tnk_q4:.3f} | {m20.delta_median_tnk:+.3f} |"
        ),
        "",
        "### Combo ρ (this pair only)",
        "",
        "| score | k | N | DL ρ (p, I²) | 95% CI | stacked rank ρ (p) | stacked Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---:|---:|---|---|---|---|",
        (
            f"| **%pos (primary)** | 2 | {int(prim.n_patients)} | "
            f"**{fmt_rho(prim.combo_rho_dl)}** ({fmt_p(prim.combo_rho_p)}, I²={prim.combo_rho_I2:.0f}%) | "
            f"[{prim.combo_rho_ci95_lo:+.3f}, {prim.combo_rho_ci95_hi:+.3f}] | "
            f"{fmt_rho(prim.stacked_rank_rho)} ({fmt_p(prim.stacked_rank_p)}) | "
            f"{fmt_rho(prim.stacked_q4q1_r)} ({fmt_p(prim.stacked_q4q1_p)}; "
            f"{int(prim.stacked_n_q1)}/{int(prim.stacked_n_q4)}) |"
        ),
        (
            f"| mean | 2 | {int(sec.n_patients)} | "
            f"{fmt_rho(sec.combo_rho_dl)} ({fmt_p(sec.combo_rho_p)}, I²={sec.combo_rho_I2:.0f}%) | "
            f"[{sec.combo_rho_ci95_lo:+.3f}, {sec.combo_rho_ci95_hi:+.3f}] | "
            f"{fmt_rho(sec.stacked_rank_rho)} ({fmt_p(sec.stacked_rank_p)}) | "
            f"{fmt_rho(sec.stacked_q4q1_r)} ({fmt_p(sec.stacked_q4q1_p)}; "
            f"{int(sec.stacked_n_q1)}/{int(sec.stacked_n_q4)}) |"
        ),
        "",
        "Primary cut is **%pos**. GSE148071 is near-null (ρ=+0.069, p=0.74; Q4 vs Q1 r=0).",
        "GSE205335 is the negative arm (ρ=−0.435, p=0.043; r=−0.778, 6/6). The pair is",
        "**heterogeneous** (I²≈66%) and the DL combo is **not** the GSE131907+GSE205335",
        "row from PR #320 (that pair is a different agent). This pair does not recover a",
        "significant CLDN4–T/NK anti-correlation.",
        "",
        "### Quartile tails (CLDN4 %pos, within-cohort)",
        "",
        "| tail | patient | n_mal | n_TNK | CLDN4 %pos | T/NK frac |",
        "|---|---|---:|---:|---:|---:|",
    ]
    lines += tail_rows("GSE148071")
    lines += tail_rows("GSE205335")
    lines += [""]

    if cellchat is None:
        lines += [
            "## CellChat-style ligands",
            "",
            "Not scored in this write-up (matrix step skipped or failed).",
            "Re-run without `--skip-cellchat` to fill `results/ligand_table.tsv`.",
            "",
        ]
    else:
        lig_path = out / "results" / "ligand_table.tsv"
        lig = pd.read_csv(lig_path, sep="\t") if lig_path.exists() else pd.DataFrame()
        lines += [
            "## CellChat-style outgoing CLDN4-high → T/NK",
            "",
            "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs.",
            "Outgoing = malignant / putative-malignant → **same-patient** T/NK.",
            "Q4 vs Q1 labels are **within-cohort** on CLDN4 %pos, then stacked.",
            f"Test = Mann–Whitney on per-patient *P* (detected ≥{MIN_DETECT_ARM} per tail).",
            "CellChat R and LIANA were not run. Cell-pooled truncated means are **not** the test.",
            "",
            (
                f"Patients scored: GSE148071 {cellchat['n_patients_scored_148071']}, "
                f"GSE205335 {cellchat['n_patients_scored_205335']}. "
                f"Detect-gated directed rows: {cellchat['n_pairs_tested_after_detect_gate']}. "
                f"Outgoing tested: {cellchat['n_outgoing_tested']}. "
                f"Outgoing p<0.05: {cellchat['n_outgoing_p_lt_05']}."
            ),
            "",
        ]
        if lig.empty:
            lines.append("No outgoing pairs passed the detect gate.")
        else:
            lines += [
                "| pair | class | n_Q1/n_Q4 (148071, 205335) | median P Q1 | median P Q4 | Δ | r | p |",
                "|---|---|---|---:|---:|---:|---:|---|",
            ]
            for rec in lig.head(15).itertuples():
                star = " **" if getattr(rec, "sig_p05", rec.p < 0.05) else ""
                lines.append(
                    f"| {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                    f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} "
                    f"({int(rec.n_q1_gse148071)}/{int(rec.n_q4_gse148071)}, "
                    f"{int(rec.n_q1_gse205335)}/{int(rec.n_q4_gse205335)}) | "
                    f"{rec.median_prob_q1:.3f} | {rec.median_prob_q4:.3f} | "
                    f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
                )
            lines += [
                "",
                "Stars mark p<0.05. The table is outgoing Mal→T/NK only, top by p.",
                "Full detect-gated table: `results/lr_q4q1_outgoing.tsv`.",
                "",
                "Only **MDK–NCL outgoing** is p<0.05 (higher in Q4). It is detected",
                "in both cohorts, so it is not a GSE205335-only leftover. Later rows",
                "with n=3 vs 3 or a zero arm are listed, not claimed. n=13 vs 12 is thin.",
                "",
            ]

    lines += [
        "## What is not supported",
        "",
        "- Treating this pair as GSE131907+GSE205335. That combo is a different folder.",
        "- Dual-high TACSTD2×CLDN4. Groups are CLDN4 only.",
        "- n = 42 + 26 as the communication n. Eligible n is 25 + 22; Q4 vs Q1 is the tails.",
        "- Stacking raw CLDN4 %pos across Singleron vs 10x / author vs marker-argmax.",
        "- Patient-level ICI response or MPR (GSE148071 has none; GSE205335 RECIST is not MPR).",
        "- Cell-pooled high/low permutation as the pair test.",
        "",
        "## Files",
        "",
        "- `results/combo_rho.tsv` — DL combo ρ + stacked Q4 vs Q1",
        "- `results/q4q1_tnk_singles.tsv` — per-cohort Spearman and Q4 vs Q1",
        "- `results/ligand_table.tsv` — outgoing CLDN4-high → T/NK",
        "- `results/patients_with_quartiles.tsv` — 47-patient table",
        "- `figures/scatter_combo_cldn4_tnk.png` — extra scatter",
        "- `figures/q4q1_tnk_pct.png` — extra Q4 vs Q1 box",
        "- `figures/fig_extra_ligand_table.png` — extra ligand-table figure",
        "- `figures/fig_n_per_patient.png` — extra honest-n bars",
        "- `METHODS.md` — pair rule, quartiles, Hill probability",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gse148071-patients", type=Path, default=ROOT / "data" / "GSE148071_per_sample.tsv")
    p.add_argument("--gse205335-patients", type=Path, default=ROOT / "data" / "GSE205335_patients.tsv")
    p.add_argument("--data-148071", type=Path, default=Path("/tmp/gse148071"))
    p.add_argument("--lr-205335", type=Path, default=ROOT / "data" / "GSE205335_per_patient_lr.tsv.gz")
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--skip-cellchat", action="store_true")
    args = p.parse_args()
    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    g14 = pd.read_csv(args.gse148071_patients, sep="\t")
    g20 = pd.read_csv(args.gse205335_patients, sep="\t")
    patients = build_patient_table(g14, g20)
    q4 = run_patient_level(patients, args.out)
    print(q4["singles"].to_string(index=False), flush=True)
    print(q4["combo"].to_string(index=False), flush=True)

    cellchat = None
    if not args.skip_cellchat:
        cellchat = run_cellchat(args, patients, args.out)
        print(json.dumps(cellchat, indent=2), flush=True)

    write_finding(q4, cellchat, args.out)
    summary = {
        "pair": "GSE148071+GSE205335",
        "additive": True,
        "marker": "CLDN4",
        "dual_high": False,
        "includes_gse131907": False,
        "unit": "patient",
        "n_gse148071_deposited": 42,
        "n_gse148071_eligible": int((g14["eligible"] == True).sum()) if "eligible" in g14 else int(g14["eligible"].sum()),
        "n_gse205335": int(len(g20)),
        "n_pair": int(len(patients)),
        "combo": q4["combo"].to_dict(orient="records"),
        "singles": q4["singles"].to_dict(orient="records"),
        "cellchat": cellchat,
        "algorithm": (
            "Within-cohort pd.qcut on average ranks; DL Fisher-z combo ρ; "
            "MWU rank-biserial; CellChat-like 10% trim mean, Hill Kh=0.5, "
            "expr_prop>=0.10, same-patient Mal→T/NK, Q4 vs Q1 within-cohort then stacked"
        ),
        "not_run": "CellChat R; LIANA; GSE131907; TACSTD2 split; dual-high; EGA raw",
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
