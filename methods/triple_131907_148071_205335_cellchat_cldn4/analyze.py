#!/usr/bin/env python3
"""Triple GSE131907+GSE148071+GSE205335: CLDN4-only patient table, then CellChat.

ADDITIVE. Patient is the unit. Locked PR #320 and GSE148071 n=25 partial
rhos are not re-audited. CellChat R is not required (Jin 2021 Hill).
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from lib_stats import q4_vs_q1, random_effects_dl, spearman, spearman_p_from_rho, stouffer

ROOT = Path(__file__).resolve().parent
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_CELLS_131907 = 20
MIN_CELLS_148071 = 25
MIN_CELLS_205335 = 20
MIN_DETECT_ARM = 3

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
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD8A", "CD3E", "IFNG", "TNF", "CD4", "NCAM1"]
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
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_rho(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:+.3f}"


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


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


def load_lr(db_dir: Path, matrix_genes: set[str]) -> pd.DataFrame:
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


def wanted_genes(db_dir: Path) -> set[str]:
    genes = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    for _, rec in inter.iterrows():
        genes.update(parse_symbols(rec.get("ligand.symbol") or rec["ligand"]))
        genes.update(parse_symbols(rec.get("receptor.symbol") or rec["receptor"]))
    return genes


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
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def compartment_gene_stats(log_cp, pos, idx: np.ndarray):
    means, props = {}, {}
    if idx.size == 0:
        return means, props
    for gene, vec in log_cp.items():
        means[gene] = trim_mean_1d(vec[idx])
        props[gene] = float(pos[gene][idx].mean())
    return means, props


def complex_from_maps(means, props, subunits: tuple[str, ...]):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    return geom_mean(vals), float(min(prs)) if prs else 0.0


def score_patient_pairs(lr: pd.DataFrame, log_cp, pos, mal_idx, tnk_idx, min_cells: int) -> list[dict]:
    rows = []
    if mal_idx.size < min_cells or tnk_idx.size < min_cells:
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
        for direction, det, lm, rm, lp, rp in (
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
                    "prob": hill_prob(lm, rm) if det else 0.0,
                    "detected": bool(det),
                    "ligand_mean": lm,
                    "receptor_mean": rm,
                    "ligand_prop": lp,
                    "receptor_prop": rp,
                }
            )
    return rows


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def stream_matrix_keep(matrix_path: Path, wanted: set[str], keep_idx: np.ndarray | None = None):
    """Stream a genes×cells TSV.gz. Optionally keep a column subset."""
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        if keep_idx is None:
            keep_idx = np.arange(n)
        n_keep = int(keep_idx.size)
        n_umi = np.zeros(n_keep, dtype=np.float64)
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
                raise ValueError(f"{matrix_path.name} {gene}: {arr.size} != {n}")
            kept = arr[keep_idx]
            n_umi += kept
            if gene in wanted and gene not in found:
                found[gene] = kept
            n_streamed += 1
            if n_streamed % 5000 == 0:
                print(f"  {matrix_path.name} genes={n_streamed} stored={len(found)}", flush=True)
    cells = [cell_ids[i] for i in keep_idx]
    return cells, found, n_umi, n_streamed


def load_rds_genes(path: Path, wanted: set[str]):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read RDS", flush=True)
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"extracted {len(extracted)} / {len(wanted)} genes from {tuple(obj.Dim)}", flush=True)
    return extracted, library_umi, barcodes, set(genes.tolist())


# ---------------------------------------------------------------------------
# Phase 1 — triple patient table (no matrix required)
# ---------------------------------------------------------------------------

def build_patient_table(data_dir: Path) -> pd.DataFrame:
    given = json.loads((data_dir / "given_rhos.json").read_text())
    rows = []

    g13 = pd.read_csv(data_dir / "GSE131907_samples.tsv", sep="\t")
    mal = g13[g13["n_malignant"] > 0].copy()
    for rec in mal.itertuples():
        rows.append(
            {
                "cohort": "GSE131907",
                "unit_id": rec.sample,
                "unit": "sample",
                "malig_def": "author_malig",
                "n_cells": int(rec.n_cells),
                "n_malignant": int(rec.n_malignant),
                "n_tnk": int(rec.n_tnk),
                "frac_tnk": float(rec.frac_tnk),
                "cldn4_mean": float(rec.mal_CLDN4_mean),
                "cldn4_pct": float(rec.mal_CLDN4_pct),
                "pass_floor": int(rec.n_malignant) >= MIN_CELLS_131907
                and int(rec.n_tnk) >= MIN_CELLS_131907,
                "floor_rule": "author_malig>0 locked PR#320; CellChat ≥20/≥20",
                "rho_source": given["GSE131907"]["source"],
            }
        )

    g20 = pd.read_csv(data_dir / "GSE205335_patients.tsv", sep="\t")
    for rec in g20.itertuples():
        rows.append(
            {
                "cohort": "GSE205335",
                "unit_id": rec.patient,
                "unit": "patient",
                "malig_def": "author_malig",
                "n_cells": int(rec.n_cells),
                "n_malignant": int(rec.n_malignant),
                "n_tnk": int(rec.n_tnk),
                "frac_tnk": float(rec.frac_tnk),
                "cldn4_mean": float(rec.mal_CLDN4_mean),
                "cldn4_pct": float(rec.mal_CLDN4_pct_pos),
                "pass_floor": int(rec.n_malignant) >= MIN_CELLS_205335
                and int(rec.n_tnk) >= MIN_CELLS_205335,
                "floor_rule": "≥20 author-mal and ≥20 T/NK locked PR#320",
                "rho_source": given["GSE205335"]["source"],
            }
        )

    g14 = pd.read_csv(data_dir / "GSE148071_per_sample.tsv", sep="\t")
    el = g14[g14["eligible"] == True].copy()  # noqa: E712
    for rec in el.itertuples():
        n_cells = int(rec.n_total)
        rows.append(
            {
                "cohort": "GSE148071",
                "unit_id": rec.sample,
                "unit": "patient",
                "malig_def": "putative_epithelial",
                "n_cells": n_cells,
                "n_malignant": int(rec.n_epithelial),
                "n_tnk": int(rec.n_TNK),
                "frac_tnk": float(rec.n_TNK) / float(n_cells) if n_cells else np.nan,
                "cldn4_mean": float(rec.mean_CLDN4_epithelial),
                "cldn4_pct": float(rec.frac_CLDN4_pos_epithelial) * 100.0,
                "pass_floor": bool(rec.eligible),
                "floor_rule": "≥25 putative epi and ≥25 T/NK (given n=25)",
                "rho_source": given["GSE148071"]["source"],
            }
        )

    table = pd.DataFrame(rows)
    table["q_pct"] = ""
    table["q_mean"] = ""
    for cohort, sub in table.groupby("cohort", sort=False):
        idx = sub.index
        table.loc[idx, "q_pct"] = assign_quartiles(sub["cldn4_pct"]).astype(str).to_numpy()
        table.loc[idx, "q_mean"] = assign_quartiles(sub["cldn4_mean"]).astype(str).to_numpy()
    return table


def combo_from_given(given: dict, patient_table: pd.DataFrame) -> pd.DataFrame:
    """Combo rho table. Locked singles are used as-is."""
    rows = []
    members = []
    for cohort in ("GSE131907", "GSE205335", "GSE148071"):
        rec = given[cohort]
        rho = float(rec["rho"])
        n = int(rec["n"])
        p = rec.get("p")
        if p is None:
            p = spearman_p_from_rho(rho, n)
        members.append((cohort, rho, float(p), n, rec.get("q4q1_r"), rec.get("n_q1"), rec.get("n_q4"), rec.get("q4q1_p")))
        rows.append(
            {
                "analysis": "spearman_single_given",
                "combo": cohort,
                "k": 1,
                "n": n,
                "n_q1": rec.get("n_q1"),
                "n_q4": rec.get("n_q4"),
                "effect": rho,
                "p": p,
                "I2": 0.0,
                "unit": rec["unit"],
                "score": rec["score"],
                "source": rec["source"],
                "re_audited": False,
            }
        )
        if rec.get("q4q1_r") is not None:
            rows.append(
                {
                    "analysis": "q4q1_single_given",
                    "combo": cohort,
                    "k": 1,
                    "n": int(rec["n_q1"]) + int(rec["n_q4"]),
                    "n_q1": rec["n_q1"],
                    "n_q4": rec["n_q4"],
                    "effect": rec["q4q1_r"],
                    "p": rec["q4q1_p"],
                    "I2": 0.0,
                    "unit": rec["unit"],
                    "score": rec["score"],
                    "source": rec["source"],
                    "re_audited": False,
                }
            )

    pair = given["pair_131907_205335_q4q1"]
    rows.append(
        {
            "analysis": "q4q1_pair_given",
            "combo": "GSE131907+GSE205335",
            "k": 2,
            "n": pair["n_compared"],
            "n_q1": pair["n_q1"],
            "n_q4": pair["n_q4"],
            "effect": pair["r"],
            "p": pair["p"],
            "I2": pair["I2"],
            "unit": "sample+patient",
            "score": "pct",
            "source": pair["source"],
            "re_audited": False,
        }
    )

    rhos = [m[1] for m in members]
    ps = [m[2] for m in members]
    ns = [m[3] for m in members]
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    rows.append(
        {
            "analysis": "spearman_triple_combo",
            "combo": "GSE131907+GSE148071+GSE205335",
            "k": re["k"],
            "n": re["n_patients_total"],
            "n_q1": None,
            "n_q4": None,
            "effect": re["pooled_rho"],
            "p": re["p"],
            "I2": re["I2"],
            "stouffer_z": st.get("z"),
            "stouffer_p": st.get("p"),
            "ci95_lo": re["ci95_rho"][0],
            "ci95_hi": re["ci95_rho"][1],
            "unit": "sample+patient",
            "score": "pct + given partial",
            "source": "NEW triple Fisher-z of locked singles",
            "re_audited": False,
        }
    )

    # Pair Spearman of the two PR #320 locked pct rows (not the Q4 pair).
    re_pair = random_effects_dl([members[0][1], members[1][1]], [members[0][3], members[1][3]])
    rows.append(
        {
            "analysis": "spearman_pair_given",
            "combo": "GSE131907+GSE205335",
            "k": 2,
            "n": members[0][3] + members[1][3],
            "n_q1": None,
            "n_q4": None,
            "effect": re_pair["pooled_rho"],
            "p": re_pair["p"],
            "I2": re_pair["I2"],
            "unit": "sample+patient",
            "score": "pct",
            "source": "PR #320 author/tnk/pct reconstructed from locked singles (not a re-audit of Q4 r=−0.705)",
            "re_audited": False,
        }
    )

    # New Q4 vs Q1 on the 25 GSE148071 eligible rows (T/NK fraction vs CLDN4 %pos).
    # This is a new cut for the triple table, not an audit of the given −0.49.
    sub = patient_table[patient_table["cohort"] == "GSE148071"]
    q = q4_vs_q1(sub["cldn4_pct"], sub["frac_tnk"])
    if q is not None:
        rows.append(
            {
                "analysis": "q4q1_148071_table_extra",
                "combo": "GSE148071",
                "k": 1,
                "n": q["n_compared"],
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "effect": q["r_rb"],
                "p": q["p"],
                "I2": 0.0,
                "unit": "patient",
                "score": "pct vs tnk_frac (table extra; not the given −0.49)",
                "source": "NEW Q4 on the n=25 eligibility table",
                "re_audited": False,
                "thin": q["thin"],
            }
        )
        # Do not Fisher-z pool this extra Q4 with the locked PR #320 pair.
        # The given −0.49 is a different cut; T/NK-fraction Q4 is reported only.

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Phase 2 — CellChat-style per-patient Hill P
# ---------------------------------------------------------------------------

def contrast_q4q1(long: pd.DataFrame) -> pd.DataFrame:
    contrast_rows = []
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
        contrast_rows.append(
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
                "n_cohorts": int(block.loc[block["detected"], "cohort"].nunique()),
                "median_prob_q1": float(d1["prob"].median()),
                "median_prob_q4": float(d4["prob"].median()),
                "delta_median": float(d4["prob"].median() - d1["prob"].median()),
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": float(r_rb),
            }
        )
    contrast = pd.DataFrame(contrast_rows)
    if contrast.empty:
        return contrast
    contrast = contrast[contrast["delta_median"] != 0].copy()
    return contrast.assign(_absd=contrast["delta_median"].abs()).sort_values(
        ["p", "_absd"], ascending=[True, False]
    ).drop(columns="_absd")


def run_cellchat_131907(args, patients: pd.DataFrame, wanted: set[str]) -> pd.DataFrame:
    keep = set(patients.loc[patients["cohort"] == "GSE131907", "unit_id"])
    ann = pd.read_csv(args.gse131907_ann, sep="\t")
    ann = ann[ann["Sample"].isin(keep)].copy()
    if ann.empty:
        raise SystemExit("GSE131907 annotation has no locked samples")
    print(f"GSE131907 keep samples={len(keep)} cells={len(ann)}", flush=True)
    with gzip.open(args.gse131907_matrix, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:] if header[0] in {"", "Index", "index", "gene", "Gene"} else header
    id_to_i = {c: i for i, c in enumerate(cell_ids)}
    # Annotation Index is usually BARCODE_SAMPLE
    keys = ann["Index"].astype(str)
    keep_idx = np.array([id_to_i[k] for k in keys if k in id_to_i], dtype=int)
    if keep_idx.size < 100:
        # try Barcode_Sample
        keys = ann["Barcode"].astype(str) + "_" + ann["Sample"].astype(str)
        keep_idx = np.array([id_to_i[k] for k in keys if k in id_to_i], dtype=int)
    if keep_idx.size < 100:
        raise SystemExit(f"GSE131907 barcode match failed ({keep_idx.size})")
    print(f"GSE131907 matched cells={keep_idx.size}", flush=True)
    cells, found, n_umi, _ = stream_matrix_keep(args.gse131907_matrix, wanted, keep_idx)
    key_order = list(cells)
    ann = ann.set_index("Index")
    if key_order[0] not in ann.index:
        ann = ann.reset_index()
        ann["Index"] = ann["Barcode"].astype(str) + "_" + ann["Sample"].astype(str)
        ann = ann.set_index("Index")
    meta = ann.loc[key_order].reset_index(drop=True)
    mal = (
        meta["Cell_subtype"].eq("Malignant cells")
        | meta["Cell_type.refined"].eq("Malignant cells")
    ).to_numpy()
    tnk = (
        meta["Cell_type"].isin(["T lymphocytes", "NK cells"])
        | meta["Cell_type.refined"].isin(["T lymphocytes", "NK cells"])
        | meta["Cell_type"].astype(str).str.contains(r"T lymph|NK", case=False, regex=True)
    ).to_numpy()
    print(f"GSE131907 labels mal={int(mal.sum())} tnk={int(tnk.sum())} types={meta['Cell_type'].value_counts().to_dict()}", flush=True)
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    pos = {g: (found[g] > 0).astype(np.float32) for g in found}
    lr = load_lr(args.db, set(found))
    qmap = patients.set_index("unit_id")["q_pct"]
    out = []
    for sample, sub in meta.groupby("Sample", observed=True):
        if sample not in keep:
            continue
        idx = sub.index.to_numpy()
        mal_idx = idx[mal[idx]]
        tnk_idx = idx[tnk[idx]]
        scored = score_patient_pairs(lr, log_cp, pos, mal_idx, tnk_idx, MIN_CELLS_131907)
        for row in scored:
            row["cohort"] = "GSE131907"
            row["unit_id"] = sample
            row["quartile_cldn4_pct"] = str(qmap.loc[sample])
            row["n_malignant"] = int(mal_idx.size)
            row["n_tnk"] = int(tnk_idx.size)
            out.append(row)
        print(f"  GSE131907 {sample} mal={mal_idx.size} tnk={tnk_idx.size}", flush=True)
    return pd.DataFrame(out)


def run_cellchat_148071(args, patients: pd.DataFrame, wanted: set[str]) -> pd.DataFrame:
    keep = set(patients.loc[patients["cohort"] == "GSE148071", "unit_id"])
    files_dir = args.gse148071_dir
    if (files_dir / "GSE148071_RAW.tar").exists() and not any(files_dir.glob("*_exp.txt.gz")):
        extract = files_dir / "GSE148071_files"
        extract.mkdir(exist_ok=True)
        if not any(extract.glob("*_exp.txt.gz")):
            with tarfile.open(files_dir / "GSE148071_RAW.tar", "r") as tf:
                tf.extractall(extract)
        files_dir = extract
    matrices = sorted(files_dir.glob("*_exp.txt.gz"))
    if not matrices:
        matrices = sorted(files_dir.rglob("*_exp.txt.gz"))
    if not matrices:
        raise SystemExit(f"no GSE148071 matrices in {files_dir}")
    # gene universe from first file
    with gzip.open(matrices[0], "rt") as handle:
        handle.readline()
        genes0 = [line.split("\t", 1)[0].strip().strip('"').split(".")[0] for line in handle]
    lr = load_lr(args.db, set(genes0))
    qmap = patients.set_index("unit_id")["q_pct"]
    out = []
    for path in matrices:
        # GSM4453576_P1_exp.txt.gz
        stem = path.name.replace("_exp.txt.gz", "")
        patient = stem.split("_", 1)[-1]
        if patient not in keep:
            continue
        cells, found, n_umi, _ = stream_matrix_keep(path, wanted)
        n = len(cells)
        lib = np.maximum(n_umi, 1.0)
        log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
        pos = {g: (found[g] > 0).astype(np.float32) for g in found}
        scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
        cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
        labels = assign_lineage(scores, cd3)
        mal_idx = np.flatnonzero(labels == "Epithelial")
        tnk_idx = np.flatnonzero(np.isin(labels, ["T", "NK"]))
        scored = score_patient_pairs(lr, log_cp, pos, mal_idx, tnk_idx, MIN_CELLS_148071)
        for row in scored:
            row["cohort"] = "GSE148071"
            row["unit_id"] = patient
            row["quartile_cldn4_pct"] = str(qmap.loc[patient])
            row["n_malignant"] = int(mal_idx.size)
            row["n_tnk"] = int(tnk_idx.size)
            out.append(row)
        print(f"  GSE148071 {patient} epi={mal_idx.size} tnk={tnk_idx.size}", flush=True)
    return pd.DataFrame(out)


def run_cellchat_205335(args, patients: pd.DataFrame, wanted: set[str]) -> pd.DataFrame:
    keep = set(patients.loc[patients["cohort"] == "GSE205335", "unit_id"])
    identities = pd.read_csv(args.gse205335_id, sep="\t")
    meta = pd.read_csv(args.gse205335_map)
    extracted, library_umi, barcodes, matrix_genes = load_rds_genes(args.gse205335_matrix, wanted)
    lr = load_lr(args.db, matrix_genes)
    indexed = identities.set_index("barcode")
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    cells = cells.merge(
        meta[["orig.ident", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    cells["compartment"] = np.select(
        [cells["lineage.sub"].eq("Malignant cells"), cells["lineage.total"].eq("T/NK cells")],
        ["Malignant", "T/NK"],
        default="Other",
    )
    lib = np.maximum(cells["total_umi"].to_numpy(), 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    qmap = patients.set_index("unit_id")["q_pct"]
    out = []
    for patient, sub in cells.groupby("patient", observed=True):
        if patient not in keep:
            continue
        idx = sub.index.to_numpy()
        mal = idx[sub["compartment"].to_numpy() == "Malignant"]
        tnk = idx[sub["compartment"].to_numpy() == "T/NK"]
        scored = score_patient_pairs(lr, log_cp, pos, mal, tnk, MIN_CELLS_205335)
        for row in scored:
            row["cohort"] = "GSE205335"
            row["unit_id"] = patient
            row["quartile_cldn4_pct"] = str(qmap.loc[patient])
            row["n_malignant"] = int(mal.size)
            row["n_tnk"] = int(tnk.size)
            out.append(row)
        print(f"  GSE205335 {patient} mal={mal.size} tnk={tnk.size}", flush=True)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Figures + FINDING
# ---------------------------------------------------------------------------

def plot_scatter_grid(table: pd.DataFrame, path: Path) -> None:
    cohorts = ["GSE131907", "GSE148071", "GSE205335"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharey=False)
    colors = {"Q1": "#6a8aaa", "Q2": "#bdbdbd", "Q3": "#f4a582", "Q4": "#b2182b"}
    for ax, cohort in zip(axes, cohorts):
        sub = table[table["cohort"] == cohort]
        for q, color in colors.items():
            m = sub["q_pct"] == q
            ax.scatter(sub.loc[m, "cldn4_pct"], sub.loc[m, "frac_tnk"], c=color, s=28, label=q)
        ax.set_title(f"{cohort} n={len(sub)}", fontsize=9)
        ax.set_xlabel("Malignant CLDN4 %pos")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Same-unit T/NK fraction")
    axes[2].legend(frameon=False, fontsize=7, title="Q")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(combo: pd.DataFrame, path: Path) -> None:
    show = combo[combo["analysis"].isin(["spearman_single_given", "spearman_triple_combo"])].copy()
    show = show.sort_values("analysis", ascending=False)
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    y = np.arange(len(show))
    ax.axvline(0, color="0.4", lw=0.8)
    ax.errorbar(show["effect"], y, fmt="o", color="#b2182b", ms=6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.combo} (n={int(r.n)})" for r in show.itertuples()], fontsize=8)
    ax.set_xlabel("Spearman ρ (CLDN4 vs T/NK)")
    ax.set_title("Triple combo — locked singles + new pool", fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_q4_boxes(table: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, cohort in zip(axes, ["GSE131907", "GSE148071", "GSE205335"]):
        sub = table[table["cohort"] == cohort]
        q1 = sub.loc[sub["q_pct"] == "Q1", "frac_tnk"].to_numpy()
        q4 = sub.loc[sub["q_pct"] == "Q4", "frac_tnk"].to_numpy()
        bp = ax.boxplot([q1, q4], tick_labels=[f"Q1 n={len(q1)}", f"Q4 n={len(q4)}"], patch_artist=True, widths=0.55)
        for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(0)
        for i, vals in enumerate((q1, q4), start=1):
            ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=14, zorder=3)
        ax.set_title(cohort, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("T/NK fraction")
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
    show = tbl.head(16).copy()
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(show) + 1.4))
    y = np.arange(len(show))
    colors = np.where(show["delta_median"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["delta_median"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.direction[:3]} {r.ligand}–{r.receptor} ({r.ligand_class})" for r in show.itertuples()],
        fontsize=8,
    )
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


def plot_n_units(table: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for cohort, sub in table.groupby("cohort"):
        ax.scatter(sub["n_malignant"], sub["n_tnk"], s=28, label=f"{cohort} n={len(sub)}")
    ax.axhline(20, color="0.6", ls="--", lw=0.7)
    ax.axvline(20, color="0.6", ls="--", lw=0.7)
    ax.set_xlabel("Malignant / putative-malignant cells")
    ax.set_ylabel("T/NK cells")
    ax.set_title("Honest n: cells per unit (floors)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(combo: pd.DataFrame, patients: pd.DataFrame, ligand: pd.DataFrame | None, out: Path) -> None:
    given_note = (
        "ADDITIVE. **CLDN4 only.** No dual-high. No GSE207422. The "
        "131907+205335-only CellChat is a different agent and is not re-run. "
        "PR #320 pair Q4 *r*=−0.705 and GSE148071 n=25 partial ρ=−0.49 are "
        "**given** and are not re-audited. Patient is the unit "
        "(GSE131907 is sample-level)."
    )
    trip = combo[combo["analysis"] == "spearman_triple_combo"].iloc[0]
    pair_s = combo[combo["analysis"] == "spearman_pair_given"].iloc[0]
    pair_q = combo[combo["analysis"] == "q4q1_pair_given"].iloc[0]
    singles = combo[combo["analysis"] == "spearman_single_given"]

    n13 = int((patients["cohort"] == "GSE131907").sum())
    n14 = int((patients["cohort"] == "GSE148071").sum())
    n20 = int((patients["cohort"] == "GSE205335").sum())
    n_floor = int(patients["pass_floor"].sum())

    lines = [
        "# FINDING — triple GSE131907 + GSE148071 + GSE205335, CLDN4-only CellChat",
        "",
        given_note,
        "",
        "p-values are descriptive.",
        "",
        "## Honest n (triple patient table first)",
        "",
        "| cohort | unit | malignant | floor | n | CellChat floor pass |",
        "|---|---|---|---|---:|---:|",
        f"| GSE131907 | sample | author Malignant cells | PR #320 extract | **{n13}** | {int(((patients.cohort=='GSE131907')&patients.pass_floor).sum())} |",
        f"| GSE148071 | patient | putative epithelium | ≥25 / ≥25 (given partial) | **{n14}** | {int(((patients.cohort=='GSE148071')&patients.pass_floor).sum())} |",
        f"| GSE205335 | patient | author malignant | ≥20 / ≥20 | **{n20}** | {int(((patients.cohort=='GSE205335')&patients.pass_floor).sum())} |",
        f"| **triple** | mixed | — | — | **{n13+n14+n20}** | **{n_floor}** |",
        "",
        "Do not write n=44 (GSE131907 series) or n=42 (GSE148071 deposited) or "
        "n=26 (GSE205335 GEO patients). The computable units are above.",
        "",
        "## Combo rho table (locked singles → new triple)",
        "",
        "Primary family is malignant CLDN4 vs same-unit T/NK. Spearman pool is "
        "DerSimonian–Laird on Fisher-z of the **given** singles (not re-audited). "
        "The PR #320 pair Q4 row is listed as given.",
        "",
        "| analysis | combo | k | N | effect (p) | source |",
        "|---|---|---:|---:|---|---|",
    ]
    for rec in combo.itertuples():
        nq = ""
        if pd.notna(getattr(rec, "n_q1", None)) and pd.notna(getattr(rec, "n_q4", None)):
            nq = f" · Q1/Q4={int(rec.n_q1)}/{int(rec.n_q4)}"
        i2 = f", I²={rec.I2:.0f}%" if pd.notna(rec.I2) and rec.k > 1 else ""
        lines.append(
            f"| {rec.analysis} | {rec.combo} | {int(rec.k)} | {int(rec.n)}{nq} | "
            f"{fmt_rho(rec.effect)} ({fmt_p(rec.p)}{i2}) | {rec.source} |"
        )

    lines += [
        "",
        f"**New triple Spearman** (k=3, N={int(trip.n)}): ρ={fmt_rho(trip.effect)} "
        f"p={fmt_p(trip.p)} I²={trip.I2:.0f}%. "
        f"Locked pair Spearman GSE131907+GSE205335 is ρ={fmt_rho(pair_s.effect)} "
        f"(PR #320 Q4 *r*={fmt_rho(pair_q.effect)} is given and not re-ranked).",
        "",
        "GSE148071 n=25 partial ρ=−0.49 is the locked single. The eligibility "
        "table’s T/NK-fraction Spearman is a different cut and is not used to "
        "replace that given ρ.",
        "",
        "### Locked singles (not re-audited)",
        "",
        "| cohort | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---:|---|---|---|",
    ]
    for rec in singles.itertuples():
        q = combo[
            (combo["analysis"] == "q4q1_single_given") & (combo["combo"] == rec.combo)
        ]
        qtxt = "—"
        if len(q):
            qq = q.iloc[0]
            qtxt = f"{fmt_rho(qq.effect)} ({fmt_p(qq.p)}; {int(qq.n_q1)}/{int(qq.n_q4)})"
        lines.append(
            f"| {rec.combo} | {int(rec.n)} | {rec.unit} | {fmt_rho(rec.effect)} ({fmt_p(rec.p)}) | {qtxt} |"
        )
    lines.append("| GSE148071 | 25 | patient | −0.490 (given partial) | — |")

    # fix duplicate 148071 if already in singles
    # singles includes 148071 from given_rhos — the extra line is redundant if so.
    # Keep one 148071 row from singles only. Remove the hardcoded extra if present in singles.
    if (singles["combo"] == "GSE148071").any():
        lines = [ln for ln in lines if ln != "| GSE148071 | 25 | patient | −0.490 (given partial) | — |"]

    lines += [
        "",
        "Full combo table: [`results/combo_rho_table.tsv`](results/combo_rho_table.tsv). "
        "Patient table: [`results/triple_patient_table.tsv`](results/triple_patient_table.tsv).",
        "",
    ]

    if ligand is None or ligand.empty:
        lines += [
            "## Ligand table (outgoing CLDN4-high malignant → T/NK)",
            "",
            "Not scored in this write-up (matrix step skipped or no pair passed the detect gate). "
            "Re-run without `--skip-cellchat`.",
            "",
        ]
    else:
        outg = ligand[ligand["direction"] == "outgoing"] if "direction" in ligand.columns else ligand
        n_sig = int((outg["p"] < 0.05).sum()) if len(outg) else 0
        lines += [
            "## Ligand table (outgoing CLDN4-high malignant → same-patient T/NK)",
            "",
            "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs. "
            "Outgoing = malignant → **same-patient** T/NK. Within-cohort CLDN4 %pos "
            "quartiles (scales are not mixed). Test = Mann–Whitney on per-patient *P* "
            f"(detected arm n≥{MIN_DETECT_ARM}). CellChat R was not run. "
            "Cell-pooled means are not the test.",
            "",
            f"Detect-gated outgoing rows: {len(outg)}. p<0.05: {n_sig}.",
            "",
            "| pair | class | n_Q1/n_Q4 | cohorts | median P Q1 | median P Q4 | Δ | r | p |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
        show = outg.head(18) if len(outg) else ligand.head(18)
        for rec in show.itertuples():
            star = " **" if rec.p < 0.05 else ""
            lines.append(
                f"| {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | "
                f"{int(getattr(rec, 'n_cohorts', 0))} | "
                f"{rec.median_prob_q1:.3f} | {rec.median_prob_q4:.3f} | "
                f"{rec.delta_median:+.3f} | {fmt_rho(rec.r_rb)} | {fmt_p(rec.p)} |"
            )
        lines += [
            "",
            "Full LR table: [`results/ligand_table.tsv`](results/ligand_table.tsv) "
            "and [`results/lr_q4q1_all.tsv`](results/lr_q4q1_all.tsv).",
            "",
            "CD274–PDCD1 / NECTIN2–TIGIT / CXCL9–CXCR3 are **not invented** if they "
            "fail `expr_prop ≥ 0.10` on either side.",
            "",
        ]

    lines += [
        "## What was not done",
        "",
        "- No dual-high TACSTD2×CLDN4 score.",
        "- No GSE207422-only CellChat.",
        "- No 131907+205335-only CellChat redo.",
        "- GSE131907 tS1–tS3 / tLung-only CellChat is the other agent’s slice.",
        "- CellChat R and LIANA were not run.",
        "",
        "## Files",
        "",
        "- `results/combo_rho_table.tsv` — locked singles + new triple pool",
        "- `results/triple_patient_table.tsv` — honest n / quartiles",
        "- `results/ligand_table.tsv` — CellChat-style outgoing table",
        "- `figures/scatter_triple_cldn4_tnk.png`",
        "- `figures/forest_triple_spearman.png`",
        "- `figures/q4q1_tnk_boxes.png`",
        "- `figures/fig_extra_ligand_table.png`",
        "- `figures/n_cells_per_unit.png`",
        "",
        "Reproduce: `python3 methods/triple_131907_148071_205335_cellchat_cldn4/analyze.py`",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=ROOT / "data")
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--gse131907-matrix", type=Path, default=Path("/tmp/geo/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"))
    p.add_argument("--gse131907-ann", type=Path, default=Path("/tmp/geo/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"))
    p.add_argument("--gse148071-dir", type=Path, default=Path("/tmp/geo/GSE148071/files"))
    p.add_argument("--gse205335-matrix", type=Path, default=Path("/tmp/geo/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz"))
    p.add_argument("--gse205335-id", type=Path, default=Path("/tmp/geo/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz"))
    p.add_argument("--gse205335-map", type=Path, default=ROOT / "data" / "gse205335_gsm_sample_metadata.csv")
    p.add_argument("--skip-cellchat", action="store_true")
    args = p.parse_args()
    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    given = json.loads((args.data / "given_rhos.json").read_text())
    patients = build_patient_table(args.data)
    patients.to_csv(args.out / "results" / "triple_patient_table.tsv", sep="\t", index=False)
    combo = combo_from_given(given, patients)
    combo.to_csv(args.out / "results" / "combo_rho_table.tsv", sep="\t", index=False)
    print(combo.to_string(index=False), flush=True)
    print(patients.groupby("cohort").size().to_string(), flush=True)

    plot_scatter_grid(patients, args.out / "figures" / "scatter_triple_cldn4_tnk.png")
    plot_forest(combo, args.out / "figures" / "forest_triple_spearman.png")
    plot_q4_boxes(patients, args.out / "figures" / "q4q1_tnk_boxes.png")
    plot_n_units(patients, args.out / "figures" / "n_cells_per_unit.png")

    ligand = None
    long = None
    if not args.skip_cellchat:
        wanted = wanted_genes(args.db)
        parts = []
        if args.gse131907_matrix.exists() and args.gse131907_ann.exists():
            parts.append(run_cellchat_131907(args, patients, wanted))
        else:
            print("skip GSE131907 CellChat (matrix missing)", flush=True)
        if args.gse148071_dir.exists():
            parts.append(run_cellchat_148071(args, patients, wanted))
        else:
            print("skip GSE148071 CellChat (matrices missing)", flush=True)
        if args.gse205335_matrix.exists() and args.gse205335_id.exists():
            parts.append(run_cellchat_205335(args, patients, wanted))
        else:
            print("skip GSE205335 CellChat (matrix missing)", flush=True)
        parts = [x for x in parts if x is not None and len(x)]
        if parts:
            long = pd.concat(parts, ignore_index=True)
            long.to_csv(args.out / "results" / "per_patient_lr.tsv.gz", sep="\t", index=False)
            contrast = contrast_q4q1(long)
            contrast.to_csv(args.out / "results" / "lr_q4q1_all.tsv", sep="\t", index=False)
            ligand = contrast[contrast["direction"] == "outgoing"].copy() if len(contrast) else contrast
            if len(ligand):
                ligand = ligand.sort_values("p")
                ligand["sig_p05"] = ligand["p"] < 0.05
            ligand.to_csv(args.out / "results" / "ligand_table.tsv", sep="\t", index=False)
            plot_ligand_table(
                ligand if ligand is not None else pd.DataFrame(),
                args.out / "figures" / "fig_extra_ligand_table.png",
                "Triple CellChat-style outgoing Mal→T/NK  Q4 vs Q1 (same-patient)",
            )

    write_finding(combo, patients, ligand, args.out)
    summary = {
        "additive": True,
        "marker": "CLDN4",
        "cohorts": ["GSE131907", "GSE148071", "GSE205335"],
        "n_units": int(len(patients)),
        "n_pass_floor": int(patients["pass_floor"].sum()),
        "combo": combo.to_dict(orient="records"),
        "cellchat": None
        if ligand is None
        else {
            "n_lr_rows": int(len(ligand)),
            "n_p_lt_05": int((ligand["p"] < 0.05).sum()) if len(ligand) else 0,
        },
        "not_run": "CellChat R; LIANA; GSE207422; dual-high; 131907+205335-only CellChat redo",
        "algorithm": "locked singles Fisher-z DL; Jin 2021 Hill Kh=0.5, 10% trim, expr_prop>=0.10, same-patient Mal→T/NK, within-cohort Q4 vs Q1",
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
