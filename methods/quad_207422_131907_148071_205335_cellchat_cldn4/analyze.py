#!/usr/bin/env python3
"""QUAD merge: GSE207422 + GSE131907 + GSE148071 + GSE205335.

ADDITIVE. CLDN4 only. Patient is the unit.
1) Patient-level malignant CLDN4 vs same-patient T/NK (combo rho + LOO + Q4 vs Q1).
2) CellChat-style outgoing CLDN4-high malignant → T/NK (Jin 2021 Hill).
GSE207422 is IN the merge. PR #320 and GSE207422-only T/NK are given, not re-audited.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import json
import re
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

ROOT = Path(__file__).resolve().parent
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_CELLS = 20
MIN_DETECT_ARM = 3
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
COHORTS = ("GSE207422", "GSE131907", "GSE148071", "GSE205335")

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
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD8A", "CD3E", "CD4", "IFNG", "TNF"]
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
COHORT_COLORS = {
    "GSE207422": "#1b9e77",
    "GSE131907": "#d95f02",
    "GSE148071": "#7570b3",
    "GSE205335": "#e7298a",
}


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str] | None = None) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"})
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


def wanted_genes(lr: pd.DataFrame, extra: bool = True) -> set[str]:
    genes = set(EXTRA) if extra else set()
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
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


def stream_gene_matrix(path: Path, wanted: set[str], gene_dotstrip: bool = False):
    with gzip.open(path, "rt") as handle:
        header = [h.strip().strip('"') for h in handle.readline().rstrip("\n").split("\t") if h != ""]
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        found: dict[str, np.ndarray] = {}
        n_streamed = 0
        all_genes: list[str] = []
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.strip().strip('"')
            if gene_dotstrip:
                gene = gene.split(".")[0]
            all_genes.append(gene)
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
            if n_streamed % 4000 == 0:
                print(f"  stream {path.name} genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done {path.name} genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed, all_genes


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


def marker_compartments(log_cp: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
    lineage = assign_lineage(scores, cd3)
    mal = lineage == "Epithelial"
    tnk = np.isin(lineage, ["T", "NK"])
    return mal, tnk


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
    if mal_idx.size < MIN_CELLS or tnk_idx.size < MIN_CELLS:
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
        base = {
            "interaction_name": rec.interaction_name,
            "pathway_name": rec.pathway_name,
            "annotation": rec.annotation,
            "ligand": rec.ligand,
            "receptor": rec.receptor,
            "ligand_genes": "|".join(rec.ligand_genes),
            "receptor_genes": "|".join(rec.receptor_genes),
        }
        rows.append(
            {
                **base,
                "direction": "outgoing",
                "prob": hill_prob(lig_mal, rec_tnk) if out_det else 0.0,
                "detected": bool(out_det),
                "ligand_mean": lig_mal,
                "receptor_mean": rec_tnk,
                "ligand_prop": lig_mal_p,
                "receptor_prop": rec_tnk_p,
            }
        )
        rows.append(
            {
                **base,
                "direction": "incoming",
                "prob": hill_prob(lig_tnk, rec_mal) if in_det else 0.0,
                "detected": bool(in_det),
                "ligand_mean": lig_tnk,
                "receptor_mean": rec_mal,
                "ligand_prop": lig_tnk_p,
                "receptor_prop": rec_mal_p,
            }
        )
    return rows


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


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def q4_vs_q1(cldn4, immune) -> dict | None:
    frame = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
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


def fisher_z_pool(rhos: list[float], ns: list[int]) -> dict:
    pairs = [(r, n) for r, n in zip(rhos, ns) if np.isfinite(r) and n > 3]
    if not pairs:
        return {"k": 0, "N": 0, "rho": float("nan"), "p": float("nan"), "I2": float("nan")}
    zs = np.array([np.arctanh(np.clip(r, -0.999999, 0.999999)) for r, _ in pairs])
    ws = np.array([n - 3 for _, n in pairs], dtype=float)
    zbar = float(np.sum(ws * zs) / np.sum(ws))
    se = float(1.0 / np.sqrt(np.sum(ws)))
    p = float(2 * stats.norm.sf(abs(zbar) / se))
    q = float(np.sum(ws * (zs - zbar) ** 2))
    df = len(pairs) - 1
    i2 = float(max(0.0, (q - df) / q)) if q > 0 and df > 0 else 0.0
    return {
        "k": len(pairs),
        "N": int(sum(n for _, n in pairs)),
        "rho": float(np.tanh(zbar)),
        "p": p,
        "I2": i2,
        "Q": q,
    }


def mwu_two_groups(q4, q1) -> dict | None:
    a = np.asarray(q4, dtype=float)
    b = np.asarray(q1, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n4, n1 = int(a.size), int(b.size)
    if n4 < 2 or n1 < 2:
        return None
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(np.median(b)),
        "median_q4": float(np.median(a)),
        "delta_median": float(np.median(a) - np.median(b)),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n1 < 3 or n4 < 3 or (n1 + n4) < 8,
    }


def parse_series_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    if "title" not in fields:
        return pd.DataFrame()
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


def load_selected_genes_rds(path: Path, wanted: set[str]):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    return extracted, library_umi, barcodes, genes.tolist()


def patient_rows_from_masks(cohort, malig_def, patient, mal, tnk, n_umi, cldn_log, cldn_raw, extra=None):
    rows = []
    for pid in pd.unique(patient):
        m = patient == pid
        mal_idx = np.flatnonzero(m & mal)
        tnk_idx = np.flatnonzero(m & tnk)
        n_mal = int(mal_idx.size)
        n_tnk = int(tnk_idx.size)
        n_cells = int(m.sum())
        rec = {
            "cohort": cohort,
            "patient": str(pid),
            "malig_def": malig_def,
            "n_cells": n_cells,
            "n_malignant": n_mal,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n_cells if n_cells else float("nan"),
            "mal_CLDN4_pct_pos": float(100.0 * (cldn_raw[mal_idx] > 0).mean()) if n_mal else float("nan"),
            "mal_CLDN4_mean": float(cldn_log[mal_idx].mean()) if n_mal else float("nan"),
            "eligible": bool(n_mal >= MIN_CELLS and n_tnk >= MIN_CELLS),
        }
        if extra:
            rec.update(extra.get(str(pid), {}))
        rows.append(rec)
    return pd.DataFrame(rows)


def score_cohort_lr(lr, log_cp, pos, patient, mal, tnk, keep_patients) -> pd.DataFrame:
    rows = []
    for pid in keep_patients:
        m = patient == pid
        mal_idx = np.flatnonzero(m & mal)
        tnk_idx = np.flatnonzero(m & tnk)
        scored = score_patient_pairs(lr, log_cp, pos, mal_idx, tnk_idx)
        for row in scored:
            row["patient"] = str(pid)
            rows.append(row)
        print(f"    scored {pid} mal={mal_idx.size} tnk={tnk_idx.size}", flush=True)
    return pd.DataFrame(rows)


def process_gse207422(data_dir: Path, db_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("== GSE207422 ==", flush=True)
    matrix = data_dir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_path = data_dir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    dummy = load_lr(db_dir)
    wanted = wanted_genes(dummy)
    cell_ids, found, n_umi, _, all_genes = stream_gene_matrix(matrix, wanted)
    lr = load_lr(db_dir, set(all_genes))
    n = len(cell_ids)
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    pos = {g: (found[g] > 0).astype(np.float32) for g in found}
    mal, tnk = marker_compartments(log_cp, n)
    sample = np.array(["_".join(bc.split("_")[:2]) if "_" in bc else bc for bc in cell_ids], dtype=object)
    meta = pd.read_excel(meta_path).dropna(subset=["Sample"]).copy()
    meta["Sample"] = meta["Sample"].astype(str)
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post-treatment", case=False)
    resp = meta["Pathologic Response"].astype(str)
    meta["response"] = np.where(resp.isin(["MPR", "pCR"]), "MPR", np.where(resp.eq("NMPR"), "NMPR", "NE"))
    mmap = meta.set_index("Sample")
    is_post = np.array([bool(mmap.loc[s, "is_post"]) if s in mmap.index else False for s in sample])
    response = np.array([mmap.loc[s, "response"] if s in mmap.index else "NA" for s in sample], dtype=object)
    patient_from_sample = np.array(
        [str(mmap.loc[s, "Patient"]) if s in mmap.index and "Patient" in mmap.columns else s for s in sample],
        dtype=object,
    )
    keep = is_post
    patient = patient_from_sample
    extra = {}
    for pid in pd.unique(patient[keep]):
        m = (patient == pid) & keep
        extra[str(pid)] = {
            "response": str(response[m][0]) if m.any() else "NA",
            "sample": ",".join(sorted(pd.unique(sample[m]))),
            "note": "post-treatment; marker epithelial as malignant proxy",
        }
    patients = patient_rows_from_masks(
        "GSE207422",
        "marker_epithelial_post",
        patient[keep],
        mal[keep],
        tnk[keep],
        n_umi[keep],
        log_cp["CLDN4"][keep],
        found["CLDN4"][keep],
        extra,
    )
    keep_ids = patients.loc[patients["eligible"], "patient"].tolist()
    log_keep = {g: v[keep] for g, v in log_cp.items()}
    pos_keep = {g: v[keep] for g, v in pos.items()}
    long = score_cohort_lr(lr, log_keep, pos_keep, patient[keep], mal[keep], tnk[keep], keep_ids)
    long["cohort"] = "GSE207422"
    return patients, long


def process_gse131907(data_dir: Path, db_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("== GSE131907 ==", flush=True)
    matrix = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann_path = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    series_path = data_dir / "GSE131907_series_matrix.txt.gz"
    dummy = load_lr(db_dir)
    wanted = wanted_genes(dummy, extra=True)
    cell_ids, found, n_umi, _, all_genes = stream_gene_matrix(matrix, wanted, gene_dotstrip=True)
    lr = load_lr(db_dir, set(all_genes))
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    per = ann.set_index("Index").reindex(cell_ids).reset_index()
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix barcodes do not align with annotation Index")
    series = parse_series_matrix(series_path)
    sample_meta = series.rename(columns={"title": "Sample"})
    pmap = (
        sample_meta.set_index("Sample")["patient_id"]
        if "patient_id" in sample_meta.columns
        else pd.Series(dtype=str)
    )
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    pos = {g: (found[g] > 0).astype(np.float32) for g in found}
    origin = per["Sample_Origin"].to_numpy()
    sample = per["Sample"].to_numpy()
    patient = np.array([pmap[s] if s in pmap.index else s for s in sample], dtype=object)
    mal = (
        (per["Cell_type"] == "Epithelial cells") & per["Cell_subtype"].isin(MALIG_SUBTYPES)
    ).to_numpy()
    tnk = per["Cell_type"].isin(["T lymphocytes", "NK cells"]).to_numpy()
    keep = np.isin(origin, TUMOR_ORIGINS)
    extra = {}
    for pid in pd.unique(patient[keep]):
        m = (patient == pid) & keep
        extra[str(pid)] = {
            "n_samples": int(pd.Series(sample[m]).nunique()),
            "origins": ",".join(sorted(pd.Series(origin[m]).dropna().unique())),
            "note": "author malignant/tS*; tumor-origin cells pooled per patient",
        }
    patients = patient_rows_from_masks(
        "GSE131907",
        "author_malig_tS_tumor_origin",
        patient[keep],
        mal[keep],
        tnk[keep],
        n_umi[keep],
        log_cp["CLDN4"][keep],
        found["CLDN4"][keep],
        extra,
    )
    keep_ids = patients.loc[patients["eligible"], "patient"].tolist()
    log_keep = {g: v[keep] for g, v in log_cp.items()}
    pos_keep = {g: v[keep] for g, v in pos.items()}
    long = score_cohort_lr(lr, log_keep, pos_keep, patient[keep], mal[keep], tnk[keep], keep_ids)
    long["cohort"] = "GSE131907"
    return patients, long


def list_gse148071(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*_exp.txt.gz"))
    if files:
        return files
    tar_path = data_dir / "GSE148071_RAW.tar"
    dest = data_dir / "GSE148071_files"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    return sorted(dest.rglob("*_exp.txt.gz"))


def process_gse148071(data_dir: Path, db_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("== GSE148071 ==", flush=True)
    files = list_gse148071(data_dir)
    dummy = load_lr(db_dir)
    wanted = wanted_genes(dummy)
    cell_ids_all: list[str] = []
    sample_all: list[str] = []
    found_all: dict[str, list[np.ndarray]] = {g: [] for g in wanted}
    n_umi_all: list[np.ndarray] = []
    gene_universe: set[str] | None = None
    for fp in files:
        m = re.search(r"_(P\d+)_", fp.name)
        patient = m.group(1) if m else fp.name.split("_")[1]
        cell_ids, found, n_umi, _, all_genes = stream_gene_matrix(fp, wanted, gene_dotstrip=True)
        if gene_universe is None:
            gene_universe = set(all_genes)
        else:
            gene_universe &= set(all_genes)
        n = len(cell_ids)
        cell_ids_all.extend([f"{patient}_{c}" for c in cell_ids])
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        for g in wanted:
            found_all[g].append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        print(f"  {fp.name} {patient} n={n}", flush=True)
    n_umi = np.concatenate(n_umi_all)
    found = {g: np.concatenate(v) for g, v in found_all.items() if v}
    patient = np.array(sample_all, dtype=object)
    lr = load_lr(db_dir, gene_universe or set())
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    pos = {g: (found[g] > 0).astype(np.float32) for g in found}
    mal, tnk = marker_compartments(log_cp, int(patient.size))
    extra = {str(p): {"note": "marker epithelial as malignant proxy; one biopsy / patient"} for p in pd.unique(patient)}
    patients = patient_rows_from_masks(
        "GSE148071",
        "marker_epithelial",
        patient,
        mal,
        tnk,
        n_umi,
        log_cp["CLDN4"],
        found["CLDN4"],
        extra,
    )
    keep_ids = patients.loc[patients["eligible"], "patient"].tolist()
    long = score_cohort_lr(lr, log_cp, pos, patient, mal, tnk, keep_ids)
    long["cohort"] = "GSE148071"
    return patients, long


def process_gse205335(data_dir: Path, db_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("== GSE205335 ==", flush=True)
    matrix = data_dir / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ident_path = data_dir / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data_dir / "GSE205335_family.soft.gz"
    dummy = load_lr(db_dir)
    wanted = wanted_genes(dummy)
    extracted, library_umi, barcodes, all_genes = load_selected_genes_rds(matrix, wanted)
    lr = load_lr(db_dir, set(all_genes))
    identities = pd.read_csv(ident_path, sep="\t")
    indexed = identities.set_index("barcode")
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    metadata = parse_geo_soft(soft_path)
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise SystemExit("GSE205335 identity/SOFT patient merge failed")
    mal = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    tnk = cells["lineage.total"].eq("T/NK cells").to_numpy()
    patient = cells["patient"].astype(str).to_numpy()
    lib = np.maximum(cells["total_umi"].to_numpy(), 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    extra = {}
    for rec in cells.drop_duplicates("patient").itertuples():
        extra[str(rec.patient)] = {
            "response": str(rec.recist),
            "histology": str(rec.cancer_subtype),
            "tissue": str(rec.tissue),
            "note": "author malignant; author T/NK",
        }
    patients = patient_rows_from_masks(
        "GSE205335",
        "author_malig",
        patient,
        mal,
        tnk,
        cells["total_umi"].to_numpy(),
        log_cp["CLDN4"],
        extracted["CLDN4"],
        extra,
    )
    keep_ids = patients.loc[patients["eligible"], "patient"].tolist()
    long = score_cohort_lr(lr, log_cp, pos, patient, mal, tnk, keep_ids)
    long["cohort"] = "GSE205335"
    return patients, long


def add_quartiles(patients: pd.DataFrame) -> pd.DataFrame:
    out = patients.copy()
    out["q_within"] = pd.Series(index=out.index, dtype=object)
    for cohort, sub in out.groupby("cohort", observed=True):
        out.loc[sub.index, "q_within"] = assign_quartiles(sub["mal_CLDN4_pct_pos"]).astype(str)
    out["q_global"] = assign_quartiles(out["mal_CLDN4_pct_pos"]).astype(str)
    return out


def rho_row(label, subset: pd.DataFrame, score="mal_CLDN4_pct_pos") -> dict:
    rho, p, n = spearman(subset[score], subset["frac_tnk"])
    q = q4_vs_q1(subset[score], subset["frac_tnk"])
    pooled = None
    if "q_within" in subset.columns:
        pooled = mwu_two_groups(
            subset.loc[subset["q_within"] == "Q4", "frac_tnk"],
            subset.loc[subset["q_within"] == "Q1", "frac_tnk"],
        )
    return {
        "analysis": label,
        "k": int(subset["cohort"].nunique()) if "cohort" in subset.columns else 1,
        "n_patients": n,
        "cohorts": ",".join(sorted(subset["cohort"].unique())) if "cohort" in subset.columns else "",
        "spearman_rho": rho,
        "spearman_p": p,
        "q4q1_within_r": pooled["r_rb"] if pooled else (q["r_rb"] if q else float("nan")),
        "q4q1_within_p": pooled["p"] if pooled else (q["p"] if q else float("nan")),
        "n_q1": pooled["n_q1"] if pooled else (q["n_q1"] if q else None),
        "n_q4": pooled["n_q4"] if pooled else (q["n_q4"] if q else None),
        "n_compared": pooled["n_compared"] if pooled else (q["n_compared"] if q else None),
        "median_tnk_q1": pooled["median_q1"] if pooled else (q["median_q1"] if q else float("nan")),
        "median_tnk_q4": pooled["median_q4"] if pooled else (q["median_q4"] if q else float("nan")),
        "delta_median_tnk": pooled["delta_median"] if pooled else (q["delta_median"] if q else float("nan")),
        "thin": pooled["thin"] if pooled else (q["thin"] if q else True),
    }


def contrast_lr(long: pd.DataFrame, patients: pd.DataFrame, qcol: str) -> pd.DataFrame:
    qmap = patients.set_index(["cohort", "patient"])[qcol].astype(str)
    work = long.copy()
    work["quartile"] = [
        qmap.loc[(r.cohort, r.patient)] if (r.cohort, r.patient) in qmap.index else "NA"
        for r in work.itertuples()
    ]
    rows = []
    for (name, direction), block in work.groupby(["interaction_name", "direction"], observed=True):
        q1 = block[(block["quartile"] == "Q1") & block["detected"]]
        q4 = block[(block["quartile"] == "Q4") & block["detected"]]
        if len(q1) < MIN_DETECT_ARM or len(q4) < MIN_DETECT_ARM:
            continue
        u, p = stats.mannwhitneyu(q4["prob"].to_numpy(), q1["prob"].to_numpy(), alternative="two-sided")
        n1, n4 = int(len(q1)), int(len(q4))
        r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
        rec0 = block.iloc[0]
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
                "n_q1_cohorts": int(q1["cohort"].nunique()),
                "n_q4_cohorts": int(q4["cohort"].nunique()),
                "median_prob_q1": float(q1["prob"].median()),
                "median_prob_q4": float(q4["prob"].median()),
                "delta_median": float(q4["prob"].median() - q1["prob"].median()),
                "p": float(p),
                "r_rb": float(r_rb),
            }
        )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out = out[out["delta_median"] != 0].copy()
    return out.assign(_abs=out["delta_median"].abs()).sort_values(["p", "_abs"], ascending=[True, False]).drop(columns="_abs")


def plot_scatter(patients, path: Path, rho, p, n) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for cohort, sub in patients.groupby("cohort", observed=True):
        ax.scatter(
            sub["mal_CLDN4_pct_pos"],
            sub["frac_tnk"],
            c=COHORT_COLORS[cohort],
            s=38,
            label=f"{cohort} n={len(sub)}",
            zorder=3,
        )
    ax.set_xlabel("Malignant CLDN4 % positive")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(f"QUAD CLDN4 vs T/NK  n={n}  ρ={rho:+.3f} p={fmt_p(p)}", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_q4q1(q1, q4, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], patch_artist=True, widths=0.55)
    for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(rows: pd.DataFrame, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 0.45 * len(rows) + 1.6))
    y = np.arange(len(rows))
    ax.axvline(0, color="0.4", lw=0.8)
    ax.scatter(rows["spearman_rho"], y, c="#b2182b", s=36, zorder=3)
    ax.set_yticks(y)
    labels = [f"{r.analysis} (n={int(r.n_patients)})" for r in rows.itertuples()]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (CLDN4 %pos vs T/NK)")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
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


def plot_n_bars(patients: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    counts = patients.groupby("cohort").size().reindex(COHORTS)
    ax.bar(range(len(counts)), counts.values, color=[COHORT_COLORS[c] for c in counts.index])
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=20, ha="right")
    ax.set_ylabel("Eligible patients")
    ax.set_title(f"QUAD honest n (floor ≥{MIN_CELLS} mal + ≥{MIN_CELLS} T/NK)", fontsize=9)
    for i, v in enumerate(counts.values):
        ax.text(i, v + 0.3, str(int(v)), ha="center", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(patients, combo, loo, fisher, ligand, out: Path) -> None:
    prim = combo.loc[combo["analysis"] == "QUAD pooled patients"].iloc[0]
    q4 = patients[patients["q_within"] == "Q4"]
    q1 = patients[patients["q_within"] == "Q1"]
    lines = [
        "# FINDING — QUAD GSE207422+GSE131907+GSE148071+GSE205335 malignant CLDN4 vs T/NK + CellChat",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high. Patient is the unit.",
        "**GSE207422 is IN the merge** (not excluded). GSE207422-only CLDN4 vs T/NK",
        "(flat) and PR #320 (GSE131907+GSE205335 Q4 r=−0.705) are given and are",
        "not re-audited or re-ranked here. p-values are descriptive.",
        "",
        "## Honest n",
        "",
        f"Eligible = ≥{MIN_CELLS} malignant and ≥{MIN_CELLS} T/NK cells in the same patient.",
        "Malignant definitions stay cohort-native: author malignant/tS* (GSE131907",
        "tumor-origin, pooled per patient), author malignant (GSE205335), marker",
        "epithelial (GSE148071; GSE207422 post-treatment; CopyKAT IDs are not public).",
        "",
        "| cohort | eligible n | malig def | dropped for floor |",
        "|---|---:|---|---:|",
    ]
    for cohort in COHORTS:
        sub = patients[patients["cohort"] == cohort]
        # patients table already filtered to eligible; n_dropped written in note via inventory later
        lines.append(
            f"| {cohort} | {len(sub)} | {sub['malig_def'].iloc[0] if len(sub) else 'NA'} | see `patients_all.tsv` |"
        )
    lines += [
        "",
        f"QUAD pooled **n={int(prim.n_patients)}** patients. Q4 vs Q1 uses",
        "**within-cohort** CLDN4 %pos quartiles, then pools the tails",
        f"(n_Q1={int(prim.n_q1)}, n_Q4={int(prim.n_q4)}, n_compared={int(prim.n_compared)}), not the mid quartiles.",
        "",
        "## Combo rho (primary = pooled patients)",
        "",
        "| analysis | k | N | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | Δ T/NK |",
        "|---|---:|---:|---|---|---:|",
    ]
    show = pd.concat([combo, loo], ignore_index=True)
    for rec in show.itertuples():
        lines.append(
            f"| {rec.analysis} | {int(rec.k)} | {int(rec.n_patients)} | "
            f"{rec.spearman_rho:+.3f} ({fmt_p(rec.spearman_p)}) | "
            f"{rec.q4q1_within_r:+.3f} ({fmt_p(rec.q4q1_within_p)}; "
            f"{int(rec.n_q1)}/{int(rec.n_q4)}) | {rec.delta_median_tnk:+.3f} |"
        )
    lines += [
        "",
        (
            f"Fisher-z meta of the four within-cohort Spearmans: "
            f"k={fisher['k']} N={fisher['N']} ρ={fisher['rho']:+.3f} "
            f"p={fmt_p(fisher['p'])} I²={100 * fisher['I2']:.0f}%."
        ),
        "",
        "Leave-one-cohort-out is the same pooled-patient Spearman after dropping",
        "that cohort. Including GSE207422 is the point of this merge.",
        "",
        "### Within-cohort tails (CLDN4 %pos, Q4 / Q1)",
        "",
        "| cohort | tail | n | median CLDN4 %pos | median T/NK |",
        "|---|---|---:|---:|---:|",
    ]
    for cohort in COHORTS:
        sub = patients[patients["cohort"] == cohort]
        for tail in ("Q4", "Q1"):
            t = sub[sub["q_within"] == tail]
            if t.empty:
                continue
            lines.append(
                f"| {cohort} | {tail} | {len(t)} | {t['mal_CLDN4_pct_pos'].median():.1f} | {t['frac_tnk'].median():.3f} |"
            )
    lines += [
        "",
        f"Pooled within-cohort tails: Q4 n={len(q4)} median T/NK={q4['frac_tnk'].median():.3f}; "
        f"Q1 n={len(q1)} median T/NK={q1['frac_tnk'].median():.3f}.",
        "",
        "## CellChat-style ligands (patients that pass floors)",
        "",
        "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs",
        "(10% truncated mean, \\(K_h=0.5\\), `expr_prop ≥ 0.10`). Outgoing =",
        "malignant → **same-patient** T/NK. Incoming = T/NK → malignant.",
        "Test = Mann–Whitney on per-patient *P* among within-cohort Q4 vs Q1",
        f"tails (detected in ≥{MIN_DETECT_ARM} Q1 and ≥{MIN_DETECT_ARM} Q4).",
        "CellChat R and LIANA were not run. Cell-pooled stacked means are not the test.",
        "",
    ]
    if ligand.empty:
        lines.append("No pairs passed the detect gate.")
    else:
        n_sig = int((ligand["p"] < 0.05).sum())
        n_out = int(((ligand["direction"] == "outgoing") & (ligand["p"] < 0.05)).sum())
        lines += [
            (
                f"Detect-gated direction×pair rows: {len(ligand)}. "
                f"p<0.05: {n_sig} (outgoing {n_out})."
            ),
            "",
            "| direction | pair | class | n_Q1/n_Q4 | median P Q1 | median P Q4 | Δ | r | p |",
            "|---|---|---|---|---:|---:|---:|---:|---|",
        ]
        show_l = ligand.head(15)
        for rec in show_l.itertuples():
            star = " **" if rec.p < 0.05 else ""
            lines.append(
                f"| {rec.direction} | {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | "
                f"{rec.median_prob_q1:.3f} | {rec.median_prob_q4:.3f} | "
                f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
            )
        lines += [
            "",
            "Stars mark p<0.05. The rest of the table is the next detect-gated",
            "pairs by p; they are not claimed. Full table:",
            "`tables/ligand_table.tsv` / `results/ligand_table.tsv`.",
            "",
        ]
    lines += [
        "## Files",
        "",
        "- `tables/combo_rho.tsv` — QUAD / per-cohort / LOO Spearman + Q4 vs Q1",
        "- `tables/ligand_table.tsv` — CellChat-style differential pairs",
        "- `results/patients_eligible.tsv` — patient table + within/global quartiles",
        "- `figures/scatter_quad_cldn4_tnk.png` — extra scatter",
        "- `figures/q4q1_tnk_within.png` — extra Q4 vs Q1 box",
        "- `figures/loo_rho_forest.png` — extra LOO forest",
        "- `figures/fig_extra_ligand_table.png` — extra ligand-table figure",
        "- `METHODS.md` — floors, labels, Hill probability",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("/tmp/quad_cldn4_geo"))
    parser.add_argument("--db", type=Path, default=ROOT / "db")
    parser.add_argument("--out", type=Path, default=ROOT)
    parser.add_argument("--skip-cellchat", action="store_true")
    args = parser.parse_args()
    results = args.out / "results"
    figures = args.out / "figures"
    tables = args.out / "tables"
    for d in (results, figures, tables):
        d.mkdir(parents=True, exist_ok=True)

    processors = {
        "GSE207422": process_gse207422,
        "GSE131907": process_gse131907,
        "GSE148071": process_gse148071,
        "GSE205335": process_gse205335,
    }
    patient_parts = []
    lr_parts = []
    inventory = {}
    for cohort, fn in processors.items():
        patients, long = fn(args.data, args.db)
        patients.to_csv(results / f"patients_all_{cohort}.tsv", sep="\t", index=False)
        elig = patients[patients["eligible"]].copy()
        elig.to_csv(results / f"patients_eligible_{cohort}.tsv", sep="\t", index=False)
        if not args.skip_cellchat and not long.empty:
            long.to_csv(results / f"per_patient_lr_{cohort}.tsv.gz", sep="\t", index=False)
        inventory[cohort] = {
            "n_scored_or_seen": int(len(patients)),
            "n_eligible": int(elig.shape[0]),
            "n_dropped_floor": int((~patients["eligible"]).sum()),
            "malig_def": str(patients["malig_def"].iloc[0]) if len(patients) else None,
        }
        print(json.dumps({cohort: inventory[cohort]}, indent=2), flush=True)
        patient_parts.append(elig)
        if not args.skip_cellchat:
            lr_parts.append(long)
        del patients, long
        gc.collect()

    patients = add_quartiles(pd.concat(patient_parts, ignore_index=True))
    patients.to_csv(results / "patients_eligible.tsv", sep="\t", index=False)
    patients.to_csv(tables / "patients_eligible.tsv", sep="\t", index=False)

    combo_rows = [rho_row("QUAD pooled patients", patients)]
    for cohort in COHORTS:
        sub = patients[patients["cohort"] == cohort]
        if len(sub):
            combo_rows.append(rho_row(f"{cohort} only", sub))
    combo = pd.DataFrame(combo_rows)
    loo_rows = []
    for cohort in COHORTS:
        sub = patients[patients["cohort"] != cohort]
        loo_rows.append(rho_row(f"LOO drop {cohort}", sub))
    loo = pd.DataFrame(loo_rows)
    per_cohort = combo[combo["analysis"].str.endswith("only")]
    fisher = fisher_z_pool(per_cohort["spearman_rho"].tolist(), per_cohort["n_patients"].tolist())
    combo_out = pd.concat([combo, loo], ignore_index=True)
    combo_out["fisher_z_rho"] = fisher["rho"]
    combo_out["fisher_z_p"] = fisher["p"]
    combo_out["fisher_z_I2"] = fisher["I2"]
    combo_out.to_csv(results / "combo_rho.tsv", sep="\t", index=False)
    combo_out.to_csv(tables / "combo_rho.tsv", sep="\t", index=False)
    print(combo_out.to_string(index=False), flush=True)
    print("fisher-z", fisher, flush=True)

    ligand = pd.DataFrame()
    if not args.skip_cellchat and lr_parts:
        long = pd.concat(lr_parts, ignore_index=True)
        ligand = contrast_lr(long, patients, "q_within")
        if not ligand.empty:
            ligand["sig_p05"] = ligand["p"] < 0.05
            ligand.to_csv(results / "lr_q4q1_all.tsv", sep="\t", index=False)
            ligand.to_csv(tables / "lr_q4q1_all.tsv", sep="\t", index=False)
        ligand.to_csv(results / "ligand_table.tsv", sep="\t", index=False)
        ligand.to_csv(tables / "ligand_table.tsv", sep="\t", index=False)
        outgoing = ligand[ligand["direction"] == "outgoing"] if not ligand.empty else ligand
        outgoing.to_csv(results / "ligand_table_outgoing.tsv", sep="\t", index=False)
        outgoing.to_csv(tables / "ligand_table_outgoing.tsv", sep="\t", index=False)

    prim = combo.loc[combo["analysis"] == "QUAD pooled patients"].iloc[0]
    plot_scatter(
        patients,
        figures / "scatter_quad_cldn4_tnk.png",
        prim.spearman_rho,
        prim.spearman_p,
        int(prim.n_patients),
    )
    plot_q4q1(
        patients.loc[patients["q_within"] == "Q1", "frac_tnk"].to_numpy(),
        patients.loc[patients["q_within"] == "Q4", "frac_tnk"].to_numpy(),
        f"QUAD within-cohort CLDN4 Q4 vs Q1 T/NK\nr={prim.q4q1_within_r:+.3f} p={fmt_p(prim.q4q1_within_p)} "
        f"n={int(prim.n_q1)}/{int(prim.n_q4)}",
        figures / "q4q1_tnk_within.png",
    )
    plot_forest(pd.concat([combo, loo], ignore_index=True), figures / "loo_rho_forest.png", "QUAD + leave-one-cohort-out Spearman")
    plot_forest(combo, figures / "cohort_rho_forest.png", "QUAD and per-cohort Spearman")
    plot_n_bars(patients, figures / "fig_n_eligible.png")
    plot_ligand_table(
        ligand if not ligand.empty else pd.DataFrame(),
        figures / "fig_extra_ligand_table.png",
        "QUAD CellChat-style Mal↔T/NK  within-cohort Q4 vs Q1 (patient unit)",
    )
    if not ligand.empty:
        plot_ligand_table(
            ligand[ligand["direction"] == "outgoing"],
            figures / "fig_extra_outgoing_lr.png",
            "QUAD outgoing Mal→T/NK  within-cohort Q4 vs Q1",
        )

    write_finding(patients, combo, loo, fisher, ligand if not ligand.empty else pd.DataFrame(), args.out)
    summary = {
        "merge": "GSE207422+GSE131907+GSE148071+GSE205335",
        "gse207422_in_merge": True,
        "additive": True,
        "marker": "CLDN4",
        "unit": "patient",
        "min_cells": MIN_CELLS,
        "inventory": inventory,
        "combo": combo_out.to_dict(orient="records"),
        "fisher_z": fisher,
        "n_lr_rows": int(len(ligand)),
        "n_lr_p_lt_05": int((ligand["p"] < 0.05).sum()) if not ligand.empty else 0,
        "algorithm": "patient-level Spearman + within-cohort qcut Q4 vs Q1 MWU; CellChat-like 10% trim mean, Hill Kh=0.5, expr_prop>=0.10, same-patient Mal→T/NK",
        "not_run": "CellChat R; LIANA; EGA raw; dual-high; re-audit of GSE207422-only T/NK or PR #320",
    }
    (results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
