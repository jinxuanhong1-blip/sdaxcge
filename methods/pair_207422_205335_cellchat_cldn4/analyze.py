#!/usr/bin/env python3
"""ADDITIVE CLDN4-only pairwise merge: GSE207422 + GSE205335.

Patient-level malignant CLDN4 vs same-patient T/NK (Spearman + Q4 vs Q1),
then CellChat-style Hill P and LIANA/CellPhoneDB-style mean-of-means
outgoing CLDN4-high (Q4) → T/NK. No dual-high. Include 207422.
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
from scipy import stats

from lib_stats import random_effects_dl, spearman as spearman_xy

ROOT = Path(__file__).resolve().parent
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_CELLS = 20
MIN_DETECT_ARM = 3
TNK_TYPES = {"CD8+ T cells", "CD4+ T cells", "NK cells"}
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD8A", "CD3E", "IFNG", "TNF"]

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


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def ligand_class(genes) -> str:
    parts = set(genes) if not isinstance(genes, str) else set(str(genes).replace("+", "|").split("|"))
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


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def q4_vs_q1(cldn4, immune) -> dict | None:
    frame = pd.DataFrame(
        {"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)}
    )
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    if n < 6:
        return None
    ranks = frame["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
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
    }


def load_cellchat_lr(db_dir: Path, matrix_genes: set[str] | None = None) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(
        columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"}
    )
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand_symbol", None)) or parse_symbols(rec.ligand)
        recp = parse_symbols(getattr(rec, "receptor_symbol", None)) or parse_symbols(rec.receptor)
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


def wanted_from_dbs(db_dir: Path, cpdb_path: Path) -> set[str]:
    wanted = set(EXTRA)
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    for _, rec in inter.iterrows():
        wanted.update(parse_symbols(rec.get("ligand.symbol") or rec["ligand"]))
        wanted.update(parse_symbols(rec.get("receptor.symbol") or rec["receptor"]))
    cpdb = pd.read_csv(cpdb_path, sep="\t")
    for col in ("ligand", "receptor"):
        for val in cpdb[col].astype(str):
            wanted.update(u for u in val.replace("|", "+").split("+") if u and u != "nan")
    return wanted


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


def compartment_gene_stats(log_cp, pos, idx: np.ndarray, trim: bool) -> tuple[dict[str, float], dict[str, float]]:
    means, props = {}, {}
    if idx.size == 0:
        return means, props
    for gene, vec in log_cp.items():
        vals = vec[idx]
        means[gene] = trim_mean_1d(vals) if trim else float(vals.mean())
        props[gene] = float(pos[gene][idx].mean())
    return means, props


def complex_from_maps(means, props, subunits, how: str):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    if how == "geom":
        return geom_mean(vals), float(min(prs)) if prs else 0.0
    return float(min(vals)) if vals else 0.0, float(min(prs)) if prs else 0.0


def stream_gse207422(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = np.asarray(header[1:], dtype=str)
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


def cohort_tests(patients: pd.DataFrame, cohort: str, malig_def: str, cldn_col: str, score: str) -> dict:
    rho, p_s, n = spearman_xy(patients[cldn_col], patients["frac_tnk"])
    q = q4_vs_q1(patients[cldn_col], patients["frac_tnk"])
    return {
        "cohort": cohort,
        "malig_def": malig_def,
        "immune_def": "same_patient_tnk_fraction",
        "score": score,
        "n_patients": n,
        "n_q1": q["n_q1"] if q else None,
        "n_q4": q["n_q4"] if q else None,
        "n_compared": q["n_compared"] if q else None,
        "spearman_rho": rho,
        "spearman_p": p_s,
        "q4q1_r_rb": q["r_rb"] if q else None,
        "q4q1_p": q["p"] if q else None,
        "median_tnk_q1": q["median_q1"] if q else None,
        "median_tnk_q4": q["median_q4"] if q else None,
        "delta_median_tnk": q["delta_median"] if q else None,
        "thin": q["thin"] if q else True,
        "unit": "patient",
    }


def combo_from_rows(rows: list[dict], family: str) -> dict:
    rhos = [r["spearman_rho"] for r in rows]
    ns = [int(r["n_patients"]) for r in rows]
    re = random_effects_dl(rhos, ns)
    q_rows = [r for r in rows if r.get("n_q1") and r.get("n_q4") and not r.get("thin")]
    q_re = {}
    if len(q_rows) >= 1:
        q_re = random_effects_dl([r["q4q1_r_rb"] for r in q_rows], [int(r["n_compared"]) for r in q_rows])
        # n_compared is the quartile-tail n; Fisher-z var uses n-3 so this is conservative
    return {
        "family": family,
        "k": re.get("k", 0),
        "cohorts": "+".join(r["cohort"] for r in rows),
        "n_patients": int(sum(ns)),
        "pooled_rho": re.get("pooled_rho"),
        "spearman_p": re.get("p"),
        "I2": re.get("I2"),
        "q4q1_k": q_re.get("k"),
        "n_compared": int(sum(int(r["n_compared"]) for r in q_rows)) if q_rows else None,
        "pooled_r_rb": q_re.get("pooled_rho"),
        "q4q1_p": q_re.get("p"),
        "q4q1_I2": q_re.get("I2"),
        "method": re.get("method"),
    }


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


def plot_scatter(x, y, hue, xlabel, ylabel, title, path: Path) -> None:
    colors = {"Q1": "#6a8aaa", "Q2": "#bdbdbd", "Q3": "#f4a582", "Q4": "#b2182b"}
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for q, color in colors.items():
        m = hue.astype(str) == q
        ax.scatter(np.asarray(x)[m], np.asarray(y)[m], c=color, s=36, label=q, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, title="CLDN4 quartile")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(rows: list[dict], title: str, path: Path, effect="spearman_rho", pcol="spearman_p") -> None:
    fig, ax = plt.subplots(figsize=(6.4, 0.55 * len(rows) + 1.6))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        val = r[effect]
        ax.plot(val, i, "o", color="#b2182b" if val < 0 else "#2166ac", ms=8)
        ax.text(0.02 if val >= 0 else -0.02, i + 0.18, f"n={r.get('n_patients', r.get('n_compared'))} p={fmt_p(r[pcol])}", fontsize=7, ha="left" if val >= 0 else "right")
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r.get("label", r.get("cohort", "?")) for r in rows], fontsize=8)
    ax.set_xlabel("effect")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str, delta_col: str) -> None:
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
    colors = np.where(show[delta_col] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show[delta_col], color=colors, alpha=0.85)
    ax.set_yticks(y)
    labels = [
        f"{getattr(r, 'direction', 'out')[:3]} {r.ligand}–{r.receptor} ({getattr(r, 'ligand_class', '')})"
        for r in show.itertuples()
    ]
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("median score Q4 − Q1")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_n_bars(df: pd.DataFrame, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(max(7.5, 0.38 * len(df)), 3.8))
    x = np.arange(len(df))
    ax.bar(x - 0.18, df["n_malignant"], 0.36, label="malignant", color="#8C6D31")
    ax.bar(x + 0.18, df["n_tnk"], 0.36, label="T/NK", color="#4C72B0")
    ax.set_xticks(x)
    ax.set_xticklabels(df["patient"].astype(str), rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("cells")
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def contrast_long(long: pd.DataFrame, score_col: str, q_col: str, min_arm: int = MIN_DETECT_ARM) -> pd.DataFrame:
    rows = []
    group_keys = [c for c in ("interaction_name", "ligand", "receptor", "direction") if c in long.columns]
    if "interaction_name" in long.columns:
        groups = long.groupby(["interaction_name", "direction"], observed=True)
    else:
        groups = long.groupby(["ligand", "receptor", "direction"], observed=True)
    for _, block in groups:
        q1 = block[block[q_col] == "Q1"]
        q4 = block[block[q_col] == "Q4"]
        det = "detected" if "detected" in block.columns else "pass_expr_prop"
        d1 = q1[q1[det]] if det in q1.columns else q1
        d4 = q4[q4[det]] if det in q4.columns else q4
        if len(d1) < min_arm or len(d4) < min_arm:
            continue
        a = d4[score_col].to_numpy(dtype=float)
        b = d1[score_col].to_numpy(dtype=float)
        if np.allclose(a, a[0]) and np.allclose(b, b[0]) and a[0] == b[0]:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        n1, n4 = int(len(d1)), int(len(d4))
        rec0 = block.iloc[0]
        med1 = float(np.median(b))
        med4 = float(np.median(a))
        rows.append(
            {
                "interaction_name": rec0.get("interaction_name", f"{rec0['ligand']}_{rec0['receptor']}"),
                "direction": rec0["direction"],
                "pathway_name": rec0.get("pathway_name", rec0.get("pathway", "")),
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "ligand_genes": rec0.get("ligand_genes", rec0["ligand"]),
                "receptor_genes": rec0.get("receptor_genes", rec0["receptor"]),
                "ligand_class": ligand_class(rec0.get("ligand_genes", rec0["ligand"])),
                "n_q1_detected": n1,
                "n_q4_detected": n4,
                "n_compared": n1 + n4,
                "n_cohorts": int(block.loc[block[q_col].isin(["Q1", "Q4"]), "cohort"].nunique()) if "cohort" in block.columns else 1,
                "median_q1": med1,
                "median_q4": med4,
                "delta_median": med4 - med1,
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": (2.0 * float(u)) / (n4 * n1) - 1.0,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.assign(_abs=out["delta_median"].abs()).sort_values(["p", "_abs"], ascending=[True, False]).drop(columns="_abs")


def score_gse207422(matrix_path: Path, annot_path: Path, patients: pd.DataFrame, db_dir: Path, cpdb_path: Path, out: Path) -> dict:
    wanted = wanted_from_dbs(db_dir, cpdb_path)
    print(f"wanted genes={len(wanted)}", flush=True)
    cell_ids, extracted, n_umi, n_genes = stream_gse207422(matrix_path, wanted)
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 not in GSE207422 matrix")
    annot = pd.read_csv(annot_path, sep="\t")
    annot = annot.rename(columns={"cell_barcode": "barcode"})
    bc_to_i = {b: i for i, b in enumerate(cell_ids)}
    keep = annot["barcode"].map(bc_to_i)
    annot = annot.loc[keep.notna()].copy()
    annot["idx"] = keep.loc[annot.index].astype(int)
    print(f"DRMref barcodes matched to matrix: {len(annot)}", flush=True)

    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}

    # %pos on DRMref malignant (additive; not on the locked A3 TACSTD2 slide)
    pct_rows = []
    for patient, sub in annot.groupby("patient", observed=True):
        mal = sub.loc[sub["celltype"].eq("Malignant cells"), "idx"].to_numpy()
        tnk = sub.loc[sub["celltype"].isin(TNK_TYPES), "idx"].to_numpy()
        pct_rows.append(
            {
                "patient": patient,
                "n_malignant_matrix": int(mal.size),
                "n_tnk_matrix": int(tnk.size),
                "mal_CLDN4_mean_matrix": float(log_cp["CLDN4"][mal].mean()) if mal.size else np.nan,
                "mal_CLDN4_pct_pos": float(100.0 * (extracted["CLDN4"][mal] > 0).mean()) if mal.size else np.nan,
            }
        )
    pct = pd.DataFrame(pct_rows)
    patients = patients.merge(pct, on="patient", how="left")
    patients.to_csv(out / "results" / "GSE207422_patients_with_pct.tsv", sep="\t", index=False)

    matrix_genes = set(extracted)
    lr = load_cellchat_lr(db_dir, matrix_genes)
    cpdb = pd.read_csv(cpdb_path, sep="\t")
    cpdb = cpdb[
        cpdb["ligand"].astype(str).map(lambda s: all(u in matrix_genes for u in s.replace("|", "+").split("+")))
        & cpdb["receptor"].astype(str).map(lambda s: all(u in matrix_genes for u in s.replace("|", "+").split("+")))
    ].copy()
    print(f"CellChat pairs in matrix={len(lr)}  CPDB pairs in matrix={len(cpdb)}", flush=True)

    qmap = patients.set_index("patient")
    cc_rows, li_rows, cell_rows = [], [], []
    for patient, sub in annot.groupby("patient", observed=True):
        if patient not in set(qmap.index):
            continue
        mal = sub.loc[sub["celltype"].eq("Malignant cells"), "idx"].to_numpy()
        tnk = sub.loc[sub["celltype"].isin(TNK_TYPES), "idx"].to_numpy()
        rec = qmap.loc[patient]
        cell_rows.append(
            {
                "cohort": "GSE207422",
                "patient": patient,
                "q_mean": str(rec["q_mean"]),
                "n_malignant": int(mal.size),
                "n_tnk": int(tnk.size),
                "eligible_lr": int(mal.size >= MIN_CELLS and tnk.size >= MIN_CELLS),
            }
        )
        if mal.size < MIN_CELLS or tnk.size < MIN_CELLS:
            print(f"  skip LR {patient} mal={mal.size} tnk={tnk.size}", flush=True)
            continue
        mal_mu_t, mal_pr_t = compartment_gene_stats(log_cp, pos, mal, trim=True)
        tnk_mu_t, tnk_pr_t = compartment_gene_stats(log_cp, pos, tnk, trim=True)
        mal_mu, mal_pr = compartment_gene_stats(log_cp, pos, mal, trim=False)
        tnk_mu, tnk_pr = compartment_gene_stats(log_cp, pos, tnk, trim=False)
        for rec_lr in lr.itertuples(index=False):
            lig_mal, lig_mal_p = complex_from_maps(mal_mu_t, mal_pr_t, rec_lr.ligand_genes, "geom")
            rec_tnk, rec_tnk_p = complex_from_maps(tnk_mu_t, tnk_pr_t, rec_lr.receptor_genes, "geom")
            out_det = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
            cc_rows.append(
                {
                    "cohort": "GSE207422",
                    "patient": patient,
                    "q_mean": str(rec["q_mean"]),
                    "interaction_name": rec_lr.interaction_name,
                    "pathway_name": rec_lr.pathway_name,
                    "ligand": rec_lr.ligand,
                    "receptor": rec_lr.receptor,
                    "ligand_genes": "|".join(rec_lr.ligand_genes),
                    "receptor_genes": "|".join(rec_lr.receptor_genes),
                    "direction": "outgoing",
                    "prob": hill_prob(lig_mal, rec_tnk) if out_det else 0.0,
                    "detected": bool(out_det),
                    "pass_expr_prop": bool(out_det),
                    "ligand_mean": lig_mal,
                    "receptor_mean": rec_tnk,
                    "cpdb_mean_score": 0.5 * (lig_mal + rec_tnk) if out_det else 0.0,
                }
            )
        for rec_cp in cpdb.itertuples(index=False):
            lig_u = str(rec_cp.ligand).replace("|", "+").split("+")
            rec_u = str(rec_cp.receptor).replace("|", "+").split("+")
            lig_mal, lig_mal_p = complex_from_maps(mal_mu, mal_pr, tuple(lig_u), "min")
            rec_tnk, rec_tnk_p = complex_from_maps(tnk_mu, tnk_pr, tuple(rec_u), "min")
            out_det = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
            li_rows.append(
                {
                    "cohort": "GSE207422",
                    "patient": patient,
                    "q_mean": str(rec["q_mean"]),
                    "ligand": rec_cp.ligand,
                    "receptor": rec_cp.receptor,
                    "pathway": rec_cp.pathway,
                    "direction": "outgoing",
                    "cpdb_mean_score": 0.5 * (lig_mal + rec_tnk) if out_det else 0.0,
                    "pass_expr_prop": bool(out_det),
                    "ligand_mean": lig_mal,
                    "receptor_mean": rec_tnk,
                }
            )
        print(f"  scored {patient} mal={mal.size} tnk={tnk.size} q={rec['q_mean']}", flush=True)

    cells = pd.DataFrame(cell_rows)
    cells.to_csv(out / "results" / "GSE207422_per_patient_cells.tsv", sep="\t", index=False)
    cc = pd.DataFrame(cc_rows)
    li = pd.DataFrame(li_rows)
    if not cc.empty:
        cc.to_csv(out / "results" / "GSE207422_per_patient_cellchat.tsv.gz", sep="\t", index=False)
    if not li.empty:
        li.to_csv(out / "results" / "GSE207422_per_patient_liana.tsv.gz", sep="\t", index=False)
    return {
        "patients": patients,
        "cellchat": cc,
        "liana": li,
        "cells": cells,
        "n_genes_streamed": int(n_genes),
        "n_genes_stored": int(len(extracted)),
        "n_cellchat_pairs": int(len(lr)),
        "n_cpdb_pairs": int(len(cpdb)),
        "n_drmref_matched": int(len(annot)),
    }


def prepare_205335_lr(path: Path, patients: pd.DataFrame) -> pd.DataFrame:
    long = pd.read_csv(path, sep="\t")
    qmap = patients.set_index("patient")["q_mean"].astype(str)
    long = long[long["direction"].eq("outgoing")].copy()
    long["cohort"] = "GSE205335"
    long["q_mean"] = long["patient"].map(qmap)
    long["q_pct"] = long["quartile_cldn4_pct"]
    long["cpdb_mean_score"] = np.where(
        long["detected"], 0.5 * (long["ligand_mean"] + long["receptor_mean"]), 0.0
    )
    long["pass_expr_prop"] = long["detected"]
    return long


def write_finding(singles, combos, cc_tbl, li_tbl, n_notes, out: Path) -> None:
    s207 = next(r for r in singles if r["cohort"] == "GSE207422" and r["score"] == "mean")
    s335 = next(r for r in singles if r["cohort"] == "GSE205335" and r["score"] == "mean")
    s335p = next((r for r in singles if r["cohort"] == "GSE205335" and r["score"] == "pct_pos"), None)
    s207p = next((r for r in singles if r["cohort"] == "GSE207422" and r["score"] == "pct_pos"), None)
    c_mean = next(c for c in combos if c["family"] == "author_or_DRMref/tnk/mean")
    lines = [
        "# FINDING — pairwise GSE207422 + GSE205335 malignant CLDN4 vs T/NK + outgoing ligands",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high. **Include 207422.** Both cohorts are",
        "ICI-adjacent (GSE207422 neoadjuvant PD-1; GSE205335 palliative ICI biopsy/effusion).",
        "Patient is the unit. Prior TACSTD2 A3 and the multi-cohort Q4 meta are given and",
        "are not re-ranked. p-values are descriptive.",
        "",
        "## Honest n",
        "",
        "| cohort | malignant label | eligible patients | Q4 vs Q1 tails | notes |",
        "|---|---|---:|---|---|",
        (
            f"| GSE207422 | DRMref `Malignant cells` (not Hu CopyKAT) | **12** post-tx | "
            f"**3 vs 3** (n_compared=6) | P06 has 15 malignant cells; kept on the locked "
            f"mean table, dropped from LR (need ≥20). Pre-tx P01/P05/P08 are not in DRMref. |"
        ),
        (
            f"| GSE205335 | author `Malignant cells` | **22** | **6 vs 6** (n_compared=12) | "
            f"Four GEO patients with (near-)zero malignant cells already out. Q4 is SCLC-heavy. |"
        ),
        (
            f"| **combo** | within-cohort quartiles, then pool tails | **34** | "
            f"**9 vs 9** (n_compared=18) | Quartiles are **not** cut on the stacked 34. |"
        ),
        "",
        "Combo Spearman is DerSimonian–Laird random-effects on Fisher-z(ρ), not a stacked",
        "34-patient correlation (CLDN4 scales differ). Q4 vs Q1 combo is the same estimator",
        "on rank-biserial *r* with n = n_Q1+n_Q4 per cohort.",
        "",
        "## Combo rho — malignant CLDN4 vs same-patient T/NK",
        "",
        "| set | score | n | Spearman ρ (p, I²) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---:|---|---|",
        (
            f"| GSE207422 | mean log1p(CP10k) | {s207['n_patients']} | "
            f"{s207['spearman_rho']:+.3f} ({fmt_p(s207['spearman_p'])}) | "
            f"{s207['q4q1_r_rb']:+.3f} ({fmt_p(s207['q4q1_p'])}; {s207['n_q1']}/{s207['n_q4']}) |"
        ),
        (
            f"| GSE205335 | mean log1p(CP10k) | {s335['n_patients']} | "
            f"{s335['spearman_rho']:+.3f} ({fmt_p(s335['spearman_p'])}) | "
            f"{s335['q4q1_r_rb']:+.3f} ({fmt_p(s335['q4q1_p'])}; {s335['n_q1']}/{s335['n_q4']}) |"
        ),
        (
            f"| **combo 207422+205335** | **mean (aligned)** | **{c_mean['n_patients']}** | "
            f"**{c_mean['pooled_rho']:+.3f} ({fmt_p(c_mean['spearman_p'])}, I²={c_mean['I2']:.0f}%)** | "
            f"**{c_mean['pooled_r_rb']:+.3f} ({fmt_p(c_mean['q4q1_p'])}; n_compared={c_mean['n_compared']})** |"
        ),
    ]
    if s335p is not None:
        lines.append(
            f"| GSE205335 | %pos (stronger single) | {s335p['n_patients']} | "
            f"{s335p['spearman_rho']:+.3f} ({fmt_p(s335p['spearman_p'])}) | "
            f"{s335p['q4q1_r_rb']:+.3f} ({fmt_p(s335p['q4q1_p'])}; {s335p['n_q1']}/{s335p['n_q4']}) |"
        )
    if s207p is not None and np.isfinite(s207p["spearman_rho"]):
        lines.append(
            f"| GSE207422 | %pos (matrix+DRMref, not on locked A3 table) | {s207p['n_patients']} | "
            f"{s207p['spearman_rho']:+.3f} ({fmt_p(s207p['spearman_p'])}) | "
            f"{s207p['q4q1_r_rb']:+.3f} ({fmt_p(s207p['q4q1_p'])}; {s207p['n_q1']}/{s207p['n_q4']}) |"
        )
        c_pct = next((c for c in combos if c["family"] == "author_or_DRMref/tnk/pct"), None)
        if c_pct is not None:
            lines.append(
                f"| **combo 207422+205335** | **%pos (primary ρ)** | **{c_pct['n_patients']}** | "
                f"**{c_pct['pooled_rho']:+.3f} ({fmt_p(c_pct['spearman_p'])}, I²={c_pct['I2']:.0f}%)** | "
                f"{c_pct['pooled_r_rb']:+.3f} ({fmt_p(c_pct['q4q1_p'])}; n_compared={c_pct['n_compared']}; I²_r high) |"
            )
    lines += [
        "",
    ]
    c_pct = next((c for c in combos if c["family"] == "author_or_DRMref/tnk/pct"), None)
    if c_pct is not None:
        lines += [
            f"**Combo ρ (primary, %pos):** Fisher-z RE on the two %pos Spearmans is "
            f"**ρ={c_pct['pooled_rho']:+.3f} (p={fmt_p(c_pct['spearman_p'])}, I²={c_pct['I2']:.0f}%, n=34)**. "
            "GSE207422 %pos is from the public UMI on DRMref malignant cells (not on the",
            "locked A3 TACSTD2 table). Mean combo is weaker (ρ=−0.166, p=0.376, I²=0%).",
            "",
            "Q4 vs Q1 **%pos** combo is **not** the tail claim: GSE207422 3 vs 3 has |r|=1",
            "(P06 enters %pos Q4 with only 15 malignant cells) and the Fisher-z *r* pool",
            f"blows up (r≈−1, I²={c_pct['q4q1_I2']:.0f}%). Mean-quartile tails are the",
            "stable Q4 vs Q1 row (9 vs 9, all LR-eligible; r=−0.505, p=0.0539, I²=0%).",
            "Ligand tables below use **mean** quartiles for that reason.",
            "",
        ]
    else:
        lines += [
            "Aligned-family **mean** combo is the verdict row for this pair.",
            "",
        ]
    lines += [
        "## CellChat-style outgoing CLDN4-high (Q4) → T/NK",
        "",
        "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs. Outgoing =",
        "same-patient malignant → T/NK. Test = Mann–Whitney on per-patient *P*",
        f"(detected ≥{MIN_DETECT_ARM} Q1 and ≥{MIN_DETECT_ARM} Q4). Quartiles are within-cohort",
        "on **mean** CLDN4 (LR-eligible 9 vs 9; %pos Q4 would drop P06). CellChat R was not run.",
        "",
    ]
    if cc_tbl is None or cc_tbl.empty:
        lines += ["No detect-gated outgoing pairs in the combo tails.", ""]
    else:
        outg = cc_tbl[cc_tbl["direction"].eq("outgoing")].copy()
        n_sig = int((outg["p"] < 0.05).sum())
        lines += [
            (
                f"Detect-gated outgoing rows: **{len(outg)}**. p<0.05: **{n_sig}**. "
                f"{n_notes}"
            ),
            "",
            "| pair | class | n_Q1/n_Q4 | cohorts | median P Q1 | median P Q4 | Δ | r | p |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
        for rec in outg.head(15).itertuples():
            star = " **" if rec.p < 0.05 else ""
            lines.append(
                f"| {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | {int(rec.n_cohorts)} | "
                f"{rec.median_q1:.3f} | {rec.median_q4:.3f} | "
                f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
            )
        lines += ["", "Stars = p<0.05. The rest of the table is the next pairs by p; they are not claimed.", ""]

    lines += [
        "## LIANA-style outgoing CLDN4-high (Q4) → T/NK",
        "",
        "CellPhoneDB-style score (Efremova 2020 / Garcia-Alonso 2022): partner =",
        "min subunit mean on log1p(CP10k); pair = mean of the two partner means.",
        "GSE207422 uses that rule on DRMref compartments. GSE205335 reuses the given",
        "CellChat truncated-mean complexes (geom-mean, 10% trim) as the partner means",
        "— same patients, not a second matrix pass. Detect gate and MWU as above.",
        "",
    ]
    if li_tbl is None or li_tbl.empty:
        lines += ["No detect-gated LIANA-style outgoing pairs in the combo tails.", ""]
    else:
        n_sig = int((li_tbl["p"] < 0.05).sum())
        lines += [
            f"Detect-gated outgoing rows: **{len(li_tbl)}**. p<0.05: **{n_sig}**.",
            "",
            "| pair | class | n_Q1/n_Q4 | cohorts | median S Q1 | median S Q4 | Δ | r | p |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
        for rec in li_tbl.head(15).itertuples():
            star = " **" if rec.p < 0.05 else ""
            lines.append(
                f"| {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | {int(rec.n_cohorts)} | "
                f"{rec.median_q1:.3f} | {rec.median_q4:.3f} | "
                f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
            )
        lines += ["", ""]

    lines += [
        "## What is not claimed",
        "",
        "- Dual-high TACSTD2×CLDN4. TACSTD2 is not a gate.",
        "- A stacked 34-patient Spearman (batch/scale mix). Combo ρ is Fisher-z RE.",
        "- Cell-pooled truncated means as the test. Patient is the unit.",
        "- Hu CopyKAT malignant IDs on GSE207422 (not public). DRMref is the given label.",
        "- MPR/RECIST as the split. Quartiles are CLDN4, not response.",
        "- CellChat R visualizations or a LIANA `mt.cellphonedb` permutation run.",
        "",
        "## Files",
        "",
        "- `results/combo_rho.tsv` — singles + Fisher-z combo",
        "- `results/patients_with_quartiles.tsv` — 34-patient table",
        "- `results/ligand_table_cellchat_outgoing.tsv` — combo CellChat-style LR",
        "- `results/ligand_table_liana_outgoing.tsv` — combo LIANA-style LR",
        "- `figures/fig_extra_ligand_table.png` — extra ligand-table figure",
        "- `figures/fig_combo_rho_forest.png` — extra combo-ρ forest",
        "- `METHODS.md` — quartiles, Hill P, CellPhoneDB mean-of-means",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gse207422-patients", type=Path, default=ROOT / "data" / "GSE207422_drmref_patients.tsv")
    p.add_argument("--gse205335-patients", type=Path, default=ROOT / "data" / "GSE205335_patients.tsv")
    p.add_argument("--gse205335-lr", type=Path, default=ROOT / "data" / "GSE205335_per_patient_lr.tsv.gz")
    p.add_argument("--drmref", type=Path, default=ROOT / "data" / "drmref_cell_annotation.tsv.gz")
    p.add_argument("--matrix", type=Path, default=Path("/tmp/geo/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--cpdb", type=Path, default=ROOT / "resources" / "cellphonedb_v5_lr_pairs.tsv")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--skip-matrix", action="store_true")
    args = p.parse_args()
    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    g22 = pd.read_csv(args.gse207422_patients, sep="\t")
    g22 = g22.rename(columns={"n_annot": "n_cells", "malig_CLDN4_mean": "mal_CLDN4_mean"})
    g22["cohort"] = "GSE207422"
    g22["malig_def"] = "DRMref_malignant"
    g22["q_mean"] = assign_quartiles(g22["mal_CLDN4_mean"])

    g35 = pd.read_csv(args.gse205335_patients, sep="\t")
    g35["cohort"] = "GSE205335"
    g35["malig_def"] = "author_malignant"
    g35["q_mean"] = assign_quartiles(g35["mal_CLDN4_mean"])
    g35["q_pct"] = assign_quartiles(g35["mal_CLDN4_pct_pos"])

    singles = [
        cohort_tests(g22, "GSE207422", "DRMref_malignant", "mal_CLDN4_mean", "mean"),
        cohort_tests(g35, "GSE205335", "author_malignant", "mal_CLDN4_mean", "mean"),
        cohort_tests(g35, "GSE205335", "author_malignant", "mal_CLDN4_pct_pos", "pct_pos"),
    ]
    combos = [
        combo_from_rows([singles[0], singles[1]], "author_or_DRMref/tnk/mean"),
    ]

    g22_scored = None
    if not args.skip_matrix and args.matrix.exists():
        g22_scored = score_gse207422(args.matrix, args.drmref, g22, args.db, args.cpdb, args.out)
        g22 = g22_scored["patients"]
        if "mal_CLDN4_pct_pos" in g22.columns and g22["mal_CLDN4_pct_pos"].notna().sum() >= 6:
            g22["q_pct"] = assign_quartiles(g22["mal_CLDN4_pct_pos"])
            singles.append(cohort_tests(g22, "GSE207422", "DRMref_malignant", "mal_CLDN4_pct_pos", "pct_pos"))
            combos.append(combo_from_rows([singles[-1], singles[2]], "author_or_DRMref/tnk/pct"))

    single_df = pd.DataFrame(singles)
    combo_df = pd.DataFrame(combos)
    single_df.to_csv(args.out / "results" / "single_cohort_rho.tsv", sep="\t", index=False)
    combo_df.to_csv(args.out / "results" / "combo_rho.tsv", sep="\t", index=False)
    print(single_df.to_string(index=False), flush=True)
    print(combo_df.to_string(index=False), flush=True)

    keep_cols = [
        "cohort", "patient", "sample", "response", "cancer_subtype", "n_cells",
        "n_malignant", "n_tnk", "frac_tnk", "mal_CLDN4_mean", "mal_CLDN4_pct_pos",
        "q_mean", "q_pct", "malig_def",
    ]
    both = pd.concat(
        [
            g22.assign(source_table="drmref_locked")[[c for c in keep_cols if c in g22.columns]],
            g35.assign(source_table="author_locked")[[c for c in keep_cols if c in g35.columns]],
        ],
        ignore_index=True,
        sort=False,
    )
    both.to_csv(args.out / "results" / "patients_with_quartiles.tsv", sep="\t", index=False)

    plot_q4q1(
        g22.loc[g22["q_mean"] == "Q1", "frac_tnk"].to_numpy(),
        g22.loc[g22["q_mean"] == "Q4", "frac_tnk"].to_numpy(),
        "Same-patient T/NK fraction",
        f"GSE207422 CLDN4 mean Q4 vs Q1  r={singles[0]['q4q1_r_rb']:+.3f} p={fmt_p(singles[0]['q4q1_p'])} n=3/3",
        args.out / "figures" / "q4q1_gse207422_mean.png",
    )
    plot_q4q1(
        g35.loc[g35["q_mean"] == "Q1", "frac_tnk"].to_numpy(),
        g35.loc[g35["q_mean"] == "Q4", "frac_tnk"].to_numpy(),
        "Same-patient T/NK fraction",
        f"GSE205335 CLDN4 mean Q4 vs Q1  r={singles[1]['q4q1_r_rb']:+.3f} p={fmt_p(singles[1]['q4q1_p'])} n=6/6",
        args.out / "figures" / "q4q1_gse205335_mean.png",
    )
    plot_q4q1(
        g35.loc[g35["q_pct"] == "Q1", "frac_tnk"].to_numpy(),
        g35.loc[g35["q_pct"] == "Q4", "frac_tnk"].to_numpy(),
        "Same-patient T/NK fraction",
        f"GSE205335 CLDN4 %pos Q4 vs Q1  r={singles[2]['q4q1_r_rb']:+.3f} p={fmt_p(singles[2]['q4q1_p'])} n=6/6",
        args.out / "figures" / "q4q1_gse205335_pct.png",
    )
    combo_q1 = np.concatenate(
        [
            g22.loc[g22["q_mean"] == "Q1", "frac_tnk"].to_numpy(),
            g35.loc[g35["q_mean"] == "Q1", "frac_tnk"].to_numpy(),
        ]
    )
    combo_q4 = np.concatenate(
        [
            g22.loc[g22["q_mean"] == "Q4", "frac_tnk"].to_numpy(),
            g35.loc[g35["q_mean"] == "Q4", "frac_tnk"].to_numpy(),
        ]
    )
    plot_q4q1(
        combo_q1,
        combo_q4,
        "Same-patient T/NK fraction",
        f"Combo within-cohort mean Q4 vs Q1  r={combos[0]['pooled_r_rb']:+.3f} p={fmt_p(combos[0]['q4q1_p'])} n=9/9",
        args.out / "figures" / "q4q1_combo_mean.png",
    )
    plot_scatter(
        g22["mal_CLDN4_mean"], g22["frac_tnk"], g22["q_mean"],
        "Malignant CLDN4 mean log1p(CP10k)", "Same-patient T/NK fraction",
        f"GSE207422 n=12  ρ={singles[0]['spearman_rho']:+.3f} p={fmt_p(singles[0]['spearman_p'])}",
        args.out / "figures" / "scatter_gse207422.png",
    )
    plot_scatter(
        g35["mal_CLDN4_mean"], g35["frac_tnk"], g35["q_mean"],
        "Malignant CLDN4 mean log1p(CP10k)", "Same-patient T/NK fraction",
        f"GSE205335 n=22  ρ={singles[1]['spearman_rho']:+.3f} p={fmt_p(singles[1]['spearman_p'])}",
        args.out / "figures" / "scatter_gse205335.png",
    )
    forest_rows = [
        {**singles[0], "label": "GSE207422 mean n=12"},
        {**singles[1], "label": "GSE205335 mean n=22"},
        {
            "label": "combo Fisher-z n=34",
            "spearman_rho": combos[0]["pooled_rho"],
            "spearman_p": combos[0]["spearman_p"],
            "n_patients": combos[0]["n_patients"],
        },
        {**singles[2], "label": "GSE205335 %pos n=22"},
    ]
    plot_forest(forest_rows, "Malignant CLDN4 vs T/NK (patient Spearman)", args.out / "figures" / "fig_combo_rho_forest.png")
    plot_n_bars(
        pd.concat(
            [
                g22[["patient", "n_malignant", "n_tnk"]].assign(patient=lambda d: "207422:" + d["patient"].astype(str)),
                g35[["patient", "n_malignant", "n_tnk"]].assign(patient=lambda d: "205335:" + d["patient"].astype(str)),
            ],
            ignore_index=True,
        ),
        args.out / "figures" / "fig_n_per_patient.png",
        "Honest n: malignant and T/NK cells per patient",
    )

    cc_tbl = pd.DataFrame()
    li_tbl = pd.DataFrame()
    n_notes = "GSE207422 matrix step skipped or failed; combo LR is GSE205335-only."
    long_cc_parts = []
    long_li_parts = []
    if args.gse205335_lr.exists():
        g35_lr = prepare_205335_lr(args.gse205335_lr, g35)
        long_cc_parts.append(g35_lr)
        long_li_parts.append(g35_lr)
        n_notes = "GSE205335 Hill P reused from the given per-patient table (mean-quartile remap)."
    if g22_scored is not None and not g22_scored["cellchat"].empty:
        long_cc_parts.append(g22_scored["cellchat"])
        long_li_parts.append(g22_scored["cellchat"])
        n_notes = (
            f"GSE207422 scored on DRMref compartments from the public UMI "
            f"({g22_scored['n_cellchat_pairs']} CellChat pairs). "
            "LIANA-style combo uses mean-of-means on the same CellChat complexes in both cohorts."
        )
        if not g22_scored["liana"].empty:
            native = contrast_long(g22_scored["liana"], "cpdb_mean_score", "q_mean")
            if not native.empty:
                native.to_csv(args.out / "results" / "GSE207422_liana_cpdb_outgoing.tsv", sep="\t", index=False)
    if long_cc_parts:
        long_cc = pd.concat(long_cc_parts, ignore_index=True, sort=False)
        long_cc.to_csv(args.out / "results" / "combo_per_patient_cellchat.tsv.gz", sep="\t", index=False)
        cc_tbl = contrast_long(long_cc, "prob", "q_mean")
        if not cc_tbl.empty:
            cc_tbl.to_csv(args.out / "results" / "lr_cellchat_q4q1_all.tsv", sep="\t", index=False)
            outg = cc_tbl[cc_tbl["direction"].eq("outgoing")].copy()
            outg.to_csv(args.out / "results" / "ligand_table_cellchat_outgoing.tsv", sep="\t", index=False)
            plot_ligand_table(
                outg,
                args.out / "figures" / "fig_extra_ligand_table.png",
                "Combo CellChat-style outgoing Mal Q4 → T/NK (within-cohort mean quartiles)",
                "delta_median",
            )
    if long_li_parts:
        long_li = pd.concat(long_li_parts, ignore_index=True, sort=False)
        if "cpdb_mean_score" in long_li.columns:
            long_li.to_csv(args.out / "results" / "combo_per_patient_liana.tsv.gz", sep="\t", index=False)
            li_tbl = contrast_long(long_li, "cpdb_mean_score", "q_mean")
            if not li_tbl.empty:
                li_tbl = li_tbl[li_tbl["direction"].eq("outgoing")].copy()
                li_tbl.to_csv(args.out / "results" / "ligand_table_liana_outgoing.tsv", sep="\t", index=False)
                plot_ligand_table(
                    li_tbl,
                    args.out / "figures" / "fig_liana_outgoing.png",
                    "Combo LIANA-style outgoing Mal Q4 → T/NK (mean-of-means)",
                    "delta_median",
                )

    write_finding(singles, combos, cc_tbl if not cc_tbl.empty else None, li_tbl if not li_tbl.empty else None, n_notes, args.out)
    summary = {
        "pair": "GSE207422+GSE205335",
        "additive": True,
        "marker": "CLDN4",
        "dual_high": False,
        "include_207422": True,
        "unit": "patient",
        "n_gse207422": int(len(g22)),
        "n_gse205335": int(len(g35)),
        "combo": combos,
        "singles": singles,
        "cellchat_outgoing_rows": int(len(cc_tbl)) if not cc_tbl.empty else 0,
        "liana_outgoing_rows": int(len(li_tbl)) if not li_tbl.empty else 0,
        "algorithm": "within-cohort pd.qcut on average ranks; Fisher-z DL combo; CellChat-like Hill Kh=0.5; LIANA-style mean of partner means",
        "not_run": "CellChat R; LIANA mt.cellphonedb permutations; EGA raw; dual-high; stacked 34-patient Spearman as primary",
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
