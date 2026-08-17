#!/usr/bin/env python3
"""CLDN4-only high-end LR on the given B-frac combo (GSE131907 + GSE241934 IIT).

Patient is the unit. Dual-high is not used. GSE207422 is not loaded.

CellChat-style: Jin 2021 Hill probability (10% truncated mean, Kh=0.5,
CellChatDB v2 protein pairs, expr_prop >= 0.10).

LIANA-style: CellPhoneDB mean-of-means on log1p(CP10k) (Efremova 2020),
CellPhoneDB v5 pairs, expr_prop >= 0.10. The LIANA R/Python package is
optional and is recorded as not-run if it does not import.

Outgoing only: CLDN4-high malignant → B, → TLS-like (B + CXCL13+ T/NK), → T/NK.
The given combo Spearman (n=32, ρ=−0.513) is cited, not re-audited.
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
from scipy.stats import trim_mean, wilcoxon
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_MAL = 10
MIN_B = 10
MIN_TLS = 10
MIN_TNK = 20
SEED = 1

# Given combo members (PR #290). Do not re-audit ρ.
GSE131907_AUTHOR_MALIG = [
    "BRONCHO_11", "BRONCHO_58", "EBUS_06", "EBUS_10", "EBUS_12", "EBUS_13",
    "EBUS_15", "EBUS_19", "EBUS_28", "EBUS_49", "EBUS_51",
    "NS_02", "NS_03", "NS_04", "NS_06", "NS_07", "NS_12", "NS_13",
    "NS_16", "NS_17", "NS_19",
]
IIT_PATIENTS = [
    "P343", "P438", "P519", "P529", "P531", "P547", "P549", "P586", "P589", "P590", "P591",
]

MALIG_SUBTYPES = {"Malignant cells"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
B_TYPES = {"B lymphocytes"}

BARRIER = {
    "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "JAM2", "JAM3",
    "CEACAM1", "CEACAM5", "CEACAM6", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4",
    "PVR", "EPCAM", "DSG2", "DSC2", "CADM1",
}
INHIB = {
    "CD274", "PDCD1LG2", "LGALS9", "HLA-E", "HLA-G", "HLA-F", "TGFB1", "TGFB2",
    "TGFB3", "CD80", "CD86", "CD276", "VSIR", "PVR", "NECTIN2", "CD47", "CDH1",
}
RECRUIT = {
    "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL3", "CCL4", "IL15",
    "IL2", "IL18", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "CXCL13",
    "CCL19", "CCL21", "CXCL12",
}
ATTACK = {"IFNG", "TNF", "FASLG", "TNFSF10", "LTA"}
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "NKG7", "MS4A1", "CD79A", "CXCL13"]


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def ligand_class(genes) -> str:
    if isinstance(genes, str):
        parts = set(genes.replace("+", "|").split("|"))
    else:
        parts = set(genes)
    tags = []
    if parts & BARRIER:
        tags.append("barrier")
    if parts & INHIB:
        tags.append("inhibitory")
    if parts & RECRUIT:
        tags.append("recruit")
    if parts & ATTACK:
        tags.append("attack")
    return "|".join(tags) if tags else "other"


def load_cellchat_lr(db_dir: Path) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand.symbol", None)) or parse_symbols(rec.ligand)
        recp = parse_symbols(getattr(rec, "receptor.symbol", None)) or parse_symbols(rec.receptor)
        if not lig or not recp:
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
                "pair_origin": "cellchatdb_v2",
            }
        )
    return pd.DataFrame(rows)


def load_cpdb_lr(res_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(res_dir / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    rows = []
    for rec in df.itertuples(index=False):
        lig = tuple(str(rec.ligand).split("+"))
        recp = tuple(str(rec.receptor).split("+"))
        if not lig or not recp:
            continue
        rows.append(
            {
                "interaction_name": f"{rec.ligand}_{rec.receptor}",
                "pathway_name": rec.pathway if hasattr(rec, "pathway") else "",
                "annotation": rec.classification if hasattr(rec, "classification") else "",
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": lig,
                "receptor_genes": recp,
                "pair_origin": getattr(rec, "pair_origin", "cellphonedb_v5"),
            }
        )
    return pd.DataFrame(rows)


def wanted_genes(cc: pd.DataFrame, cp: pd.DataFrame) -> set[str]:
    genes = set(EXTRA)
    for df in (cc, cp):
        for rec in df.itertuples(index=False):
            genes.update(rec.ligand_genes)
            genes.update(rec.receptor_genes)
    return genes


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def geom_mean(vals: list[float]) -> float:
    arr = np.asarray(vals, dtype=float)
    if np.any(arr <= 0):
        return 0.0
    return float(np.exp(np.mean(np.log(arr))))


def group_stats(logx: dict[str, np.ndarray], pos: dict[str, np.ndarray], mask: np.ndarray) -> dict:
    n = int(mask.sum())
    out = {"n": n, "mean": {}, "trim": {}, "prop": {}}
    if n == 0:
        return out
    for g, arr in logx.items():
        v = arr[mask]
        out["mean"][g] = float(np.mean(v))
        out["trim"][g] = float(v[0] if n == 1 else trim_mean(v, TRIM))
        out["prop"][g] = float(np.mean(pos[g][mask]))
    return out


def partner_value(stats: dict, subunits: tuple[str, ...], how: str) -> tuple[float, float]:
    if any(g not in stats["mean"] for g in subunits):
        return np.nan, np.nan
    props = [stats["prop"][g] for g in subunits]
    if how == "trim":
        mus = [stats["trim"][g] for g in subunits]
        mu = geom_mean(mus) if len(mus) > 1 else mus[0]
    else:
        mus = [stats["mean"][g] for g in subunits]
        mu = float(min(mus))
    return float(mu), float(min(props))


def score_pairs(pairs: pd.DataFrame, sender: dict, receiver: dict, method: str) -> pd.DataFrame:
    how = "trim" if method == "cellchat" else "mean"
    rows = []
    for rec in pairs.itertuples(index=False):
        l_mu, l_pr = partner_value(sender, rec.ligand_genes, how)
        r_mu, r_pr = partner_value(receiver, rec.receptor_genes, how)
        if not np.isfinite(l_mu) or not np.isfinite(r_mu):
            continue
        detected = (l_pr >= EXPR_PROP) and (r_pr >= EXPR_PROP)
        if method == "cellchat":
            score = hill_prob(l_mu, r_mu) if detected else 0.0
        else:
            score = 0.5 * (l_mu + r_mu) if detected else 0.0
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "ligand_class": ligand_class(rec.ligand_genes),
                "pair_origin": rec.pair_origin,
                "method": method,
                "n_sender": sender["n"],
                "n_receiver": receiver["n"],
                "ligand_mean": l_mu,
                "receptor_mean": r_mu,
                "ligand_prop": l_pr,
                "receptor_prop": r_pr,
                "detected": bool(detected),
                "score": float(score),
            }
        )
    return pd.DataFrame(rows)


def parse_series_matrix(path: Path) -> pd.DataFrame:
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
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def stream_gene_cell_tsv(path: Path, genes: set[str]):
    log(f"[stream-tsv] {path.name} keep={len(genes)}")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n = len(cells)
        lib = np.zeros(n, dtype=np.float64)
        keep: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                continue
            lib += vals
            n_genes += 1
            if gene in genes and gene not in keep:
                keep[gene] = vals
            if n_genes % 4000 == 0:
                log(f"  genes={n_genes} stored={len(keep)}")
    log(f"[stream-tsv] cells={n} genes={n_genes} stored={len(keep)}")
    return np.array(cells, dtype=object), keep, lib.astype(np.float32), n_genes


def stream_mtx(mtx: Path, features: Path, barcodes: Path, genes: set[str], lib_meta=None):
    feats = pd.read_csv(features, sep="\t", header=None, dtype=str)
    if feats.shape[1] >= 2:
        symbols = feats[1].astype(str).str.split(".").str[0]
        alt = feats[0].astype(str).str.split(".").str[0]
        symbols = np.where(symbols.isin(genes), symbols, alt)
    else:
        symbols = feats[0].astype(str).str.split(".").str[0].to_numpy()
    symbols = np.asarray(symbols)
    want_rows = {i for i, s in enumerate(symbols) if s in genes}
    bcs = pd.read_csv(barcodes, sep="\t", header=None, dtype=str)[0].to_numpy()
    n = len(bcs)
    keep = {s: np.zeros(n, dtype=np.float32) for s in set(symbols[list(want_rows)])} if want_rows else {}
    lib = None if lib_meta is not None else np.zeros(n, dtype=np.float64)
    log(f"[stream-mtx] {mtx.name} cells={n} want_rows={len(want_rows)}")
    with gzip.open(mtx, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            break
        n_lines = 0
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            v = float(parts[2])
            if lib is not None:
                lib[c] += v
            if r in want_rows:
                keep[symbols[r]][c] += v
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                log(f"  mtx entries={n_lines}")
    if lib is None:
        lib = np.asarray(lib_meta, dtype=np.float32)
    else:
        lib = lib.astype(np.float32)
    log(f"[stream-mtx] stored={len(keep)}")
    return bcs, keep, lib


def log1p_cp10k(umi: dict[str, np.ndarray], lib: np.ndarray) -> dict[str, np.ndarray]:
    scale = np.where(lib > 0, 1e4 / lib, 0.0).astype(np.float32)
    return {g: np.log1p(arr * scale).astype(np.float32) for g, arr in umi.items()}


def load_gse131907(data: Path, genes: set[str]) -> dict:
    ann = pd.read_csv(data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    series = parse_series_matrix(data / "GSE131907" / "GSE131907_series_matrix.txt.gz")
    series = series.rename(columns={"title": "Sample"})
    pmap = {}
    if "patient_id" in series.columns:
        pmap = series.drop_duplicates("Sample").set_index("Sample")["patient_id"].to_dict()
    cells, umi, lib, n_genes = stream_gene_cell_tsv(
        data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", genes
    )
    per = ann.set_index("Index").reindex(cells)
    sample = per["Sample"].fillna("NA").astype(str).to_numpy()
    keep_s = np.isin(sample, GSE131907_AUTHOR_MALIG)
    idx = np.flatnonzero(keep_s)
    cells = cells[idx]
    sample = sample[idx]
    lib = lib[idx]
    umi = {g: v[idx] for g, v in umi.items()}
    per = per.iloc[idx]
    origin = per["Sample_Origin"].fillna("").astype(str).to_numpy()
    ctype = per["Cell_type"].fillna("").astype(str).to_numpy()
    subtype = per["Cell_subtype"].fillna("").astype(str).to_numpy()
    patient = np.array([pmap.get(s, s) for s in sample], dtype=object)
    is_mal = (ctype == "Epithelial cells") & np.isin(subtype, list(MALIG_SUBTYPES))
    is_b = np.isin(ctype, list(B_TYPES))
    is_t = ctype == "T lymphocytes"
    is_nk = ctype == "NK cells"
    is_tnk = is_t | is_nk
    logx = log1p_cp10k(umi, lib)
    pos = {g: (umi[g] > 0).astype(np.float32) for g in umi}
    cxcl13 = umi["CXCL13"] > 0 if "CXCL13" in umi else np.zeros(len(cells), dtype=bool)
    is_tls = is_b | (is_tnk & cxcl13)
    log(
        f"[GSE131907] kept-combo-samples cells={len(cells)} mal={int(is_mal.sum())} "
        f"B={int(is_b.sum())} TNK={int(is_tnk.sum())} TLS={int(is_tls.sum())} "
        f"patients={pd.Series(patient).nunique()}"
    )
    return {
        "cohort": "GSE131907",
        "cells": cells,
        "patient": patient,
        "sample": sample,
        "origin": origin,
        "is_mal": is_mal,
        "is_b": is_b,
        "is_tnk": is_tnk,
        "is_tls": is_tls,
        "logx": logx,
        "pos": pos,
        "cldn4": logx.get("CLDN4", np.zeros(len(cells), dtype=np.float32)),
        "n_genes_matrix": n_genes,
        "n_cells_matrix": 208506,
        "mal_def": "author Cell_subtype == Malignant cells (combo author_malig n=21)",
    }


def _map_iit_lineage(major: str, fine: str = "") -> str:
    fine_l = str(fine).lower()
    if "plasma" in fine_l:
        return "plasma"
    t = str(major).strip().lower()
    if t in {"epi", "epithelial", "cancer", "tumor", "malig", "malignant"}:
        return "epithelial"
    if t in {"b", "bcell", "b_cell", "b-cell"} or t.startswith("b"):
        return "B"
    if t in {"t", "tcell", "t_cell"} or t.startswith("t"):
        return "T"
    if t in {"nk", "nkcell"}:
        return "NK"
    return "other"


def load_gse241934_iit(data: Path, genes: set[str]) -> dict:
    d = data / "GSE241934_IIT"
    meta = pd.read_csv(d / "GSE241934_IIT_Meta.txt.gz", sep="\t")
    log(f"[IIT meta] columns={list(meta.columns)[:30]} n={len(meta)}")
    key = "cellID" if "cellID" in meta.columns else ("barcode" if "barcode" in meta.columns else meta.columns[0])
    sid_col = "sampleID" if "sampleID" in meta.columns else ("orig.ident" if "orig.ident" in meta.columns else key)
    major_col = None
    for c in ("major_cell_type", "major.cell.type", "cell.type", "cell_type", "celltype", "CellType"):
        if c in meta.columns:
            major_col = c
            break
    fine_col = "cell.type" if "cell.type" in meta.columns and major_col != "cell.type" else None
    bcs_file = pd.read_csv(d / "GSE241934_IIT_barcodes.tsv.gz", sep="\t", header=None, dtype=str)[0]
    lib_meta = None
    if "nCount_RNA" in meta.columns:
        lookup = dict(zip(meta[key].astype(str), meta["nCount_RNA"].astype(float)))
        lib_meta = np.array([lookup.get(str(b), np.nan) for b in bcs_file], dtype=np.float64)
        if not np.isfinite(lib_meta).all():
            # try suffix match
            suffix = {str(k).split("_")[-1]: v for k, v in lookup.items()}
            lib_meta = np.array(
                [lookup.get(str(b), suffix.get(str(b).split("_")[-1], np.nan)) for b in bcs_file],
                dtype=np.float64,
            )
        if not np.isfinite(lib_meta).all():
            lib_meta = None
    cells, umi, lib = stream_mtx(
        d / "GSE241934_IIT_Matrix.mtx.gz",
        d / "GSE241934_IIT_features.tsv.gz",
        d / "GSE241934_IIT_barcodes.tsv.gz",
        genes,
        lib_meta=lib_meta,
    )
    meta_idx = meta.set_index(meta[key].astype(str))
    patient = []
    major_vals = []
    lineage = []
    for cell in cells:
        row = None
        if cell in meta_idx.index:
            row = meta_idx.loc[cell]
        else:
            tok = str(cell).split("_")[-1]
            if tok in meta_idx.index:
                row = meta_idx.loc[tok]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        if row is None:
            patient.append("unknown")
            major_vals.append("")
            lineage.append("other")
            continue
        patient.append(str(row[sid_col]))
        maj = str(row[major_col]) if major_col else ""
        fine = str(row[fine_col]) if fine_col else ""
        major_vals.append(maj)
        lineage.append(_map_iit_lineage(maj, fine))
    patient = np.array(patient, dtype=object)
    lineage = np.array(lineage, dtype=object)
    # keep combo IIT patients when present; otherwise keep any patient with epi
    keep_p = np.isin(patient, IIT_PATIENTS)
    if keep_p.sum() == 0:
        log("[IIT] listed patient IDs not found; using all meta patients with epithelial cells")
        keep_p = np.ones(len(patient), dtype=bool)
    idx = np.flatnonzero(keep_p)
    cells = cells[idx]
    patient = patient[idx]
    lineage = lineage[idx]
    lib = lib[idx]
    umi = {g: v[idx] for g, v in umi.items()}
    is_mal = lineage == "epithelial"
    is_b = lineage == "B"
    is_tnk = np.isin(lineage, ["T", "NK"])
    logx = log1p_cp10k(umi, lib)
    pos = {g: (umi[g] > 0).astype(np.float32) for g in umi}
    cxcl13 = umi["CXCL13"] > 0 if "CXCL13" in umi else np.zeros(len(cells), dtype=bool)
    is_tls = is_b | (is_tnk & cxcl13)
    log(
        f"[GSE241934_IIT] cells={len(cells)} mal/epi={int(is_mal.sum())} "
        f"B={int(is_b.sum())} TNK={int(is_tnk.sum())} TLS={int(is_tls.sum())} "
        f"patients={pd.Series(patient).nunique()} majors={sorted(set(major_vals))[:12]}"
    )
    return {
        "cohort": "GSE241934_IIT",
        "cells": cells,
        "patient": patient,
        "sample": patient,
        "origin": np.array(["IIT"] * len(cells), dtype=object),
        "is_mal": is_mal,
        "is_b": is_b,
        "is_tnk": is_tnk,
        "is_tls": is_tls,
        "logx": logx,
        "pos": pos,
        "cldn4": logx.get("CLDN4", np.zeros(len(cells), dtype=np.float32)),
        "n_genes_matrix": None,
        "n_cells_matrix": int(len(keep_p)),
        "mal_def": "author Epi (combo author_epi n=11)",
        "author_majors": sorted({m for m in major_vals if m}),
    }


def split_high_low(cldn4: np.ndarray, mal: np.ndarray, rule: str) -> tuple[np.ndarray, np.ndarray, dict]:
    vals = cldn4[mal]
    hi = np.zeros(len(cldn4), dtype=bool)
    lo = np.zeros(len(cldn4), dtype=bool)
    info = {"rule": rule, "n_mal": int(mal.sum())}
    if vals.size < MIN_MAL * 2:
        return hi, lo, info
    if rule == "median":
        thr = float(np.median(vals))
        hi[mal & (cldn4 >= thr)] = True
        lo[mal & (cldn4 < thr)] = True
        info["threshold"] = thr
    else:
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        hi[mal & (cldn4 >= q2)] = True
        lo[mal & (cldn4 <= q1)] = True
        info["q_low"] = float(q1)
        info["q_high"] = float(q2)
    info["n_high"] = int(hi.sum())
    info["n_low"] = int(lo.sum())
    return hi, lo, info


def score_patient(obj: dict, pid: str, pairs_cc: pd.DataFrame, pairs_cp: pd.DataFrame, rule: str) -> tuple[list[dict], dict]:
    m = obj["patient"] == pid
    mal = obj["is_mal"] & m
    hi, lo, info = split_high_low(obj["cldn4"], mal, rule)
    dests = {
        "B": obj["is_b"] & m,
        "TLS": obj["is_tls"] & m,
        "TNK": obj["is_tnk"] & m,
    }
    mins = {"B": MIN_B, "TLS": MIN_TLS, "TNK": MIN_TNK}
    rec = {
        "cohort": obj["cohort"],
        "patient": pid,
        "rule": rule,
        "n_mal": int(mal.sum()),
        "n_high": int(hi.sum()),
        "n_low": int(lo.sum()),
        "n_B": int(dests["B"].sum()),
        "n_TLS": int(dests["TLS"].sum()),
        "n_TNK": int(dests["TNK"].sum()),
        "mean_CLDN4_high": float(obj["cldn4"][hi].mean()) if hi.any() else np.nan,
        "mean_CLDN4_low": float(obj["cldn4"][lo].mean()) if lo.any() else np.nan,
        "eligible_mal_split": int(hi.sum()) >= MIN_MAL and int(lo.sum()) >= MIN_MAL,
    }
    rec.update({f"eligible_{k}": rec["eligible_mal_split"] and int(dests[k].sum()) >= mins[k] for k in dests})
    rows = []
    if not rec["eligible_mal_split"]:
        return rows, rec
    hi_s = group_stats(obj["logx"], obj["pos"], hi)
    lo_s = group_stats(obj["logx"], obj["pos"], lo)
    for dest, mask in dests.items():
        if not rec[f"eligible_{dest}"]:
            continue
        recv = group_stats(obj["logx"], obj["pos"], mask)
        for method, pairs in (("cellchat", pairs_cc), ("liana_cpdb", pairs_cp)):
            a = score_pairs(pairs, hi_s, recv, "cellchat" if method == "cellchat" else "liana")
            b = score_pairs(pairs, lo_s, recv, "cellchat" if method == "cellchat" else "liana")
            if a.empty or b.empty:
                continue
            merged = a.merge(
                b,
                on=["interaction_name", "pathway_name", "annotation", "ligand", "receptor",
                    "ligand_genes", "receptor_genes", "ligand_class", "pair_origin", "method"],
                suffixes=("_high", "_low"),
            )
            merged["delta"] = merged["score_high"] - merged["score_low"]
            merged["pass_either"] = merged["detected_high"] | merged["detected_low"]
            merged["method"] = method
            merged["cohort"] = obj["cohort"]
            merged["patient"] = pid
            merged["destination"] = dest
            merged["rule"] = rule
            rows.append(merged)
    return rows, rec


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 6 or np.allclose(a, b, equal_nan=False):
        return np.nan
    try:
        return float(wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def aggregate_patients(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return raw
    rows = []
    keys = ["method", "destination", "rule", "interaction_name", "pathway_name",
            "ligand", "receptor", "ligand_genes", "receptor_genes", "ligand_class", "pair_origin"]
    for key, sub in raw.groupby(keys, observed=True, dropna=False):
        keep = sub[sub["pass_either"]]
        if len(keep) < 3:
            continue
        p = wilcoxon_safe(keep["score_high"], keep["score_low"])
        rec = dict(zip(keys, key))
        rec.update(
            {
                "n_patients": int(keep["patient"].nunique()),
                "n_cohorts": int(keep["cohort"].nunique()),
                "cohorts": ",".join(sorted(keep["cohort"].unique())),
                "median_delta": float(np.median(keep["delta"])),
                "mean_delta": float(np.mean(keep["delta"])),
                "mean_score_high": float(keep["score_high"].mean()),
                "mean_score_low": float(keep["score_low"].mean()),
                "frac_detected_high": float(keep["detected_high"].mean()),
                "frac_detected_low": float(keep["detected_low"].mean()),
                "median_n_sender_high": float(keep["n_sender_high"].median()),
                "median_n_receiver": float(keep["n_receiver_high"].median()),
                "pval": p,
            }
        )
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["padj"] = np.nan
    for (method, dest, rule), idx in out.groupby(["method", "destination", "rule"]).groups.items():
        mask = out.index.isin(idx) & out["pval"].notna()
        if mask.sum() == 0:
            continue
        out.loc[mask, "padj"] = multipletests(out.loc[mask, "pval"], method="fdr_bh")[1]
    out = out.sort_values(["method", "destination", "rule", "median_delta"])
    return out


def plot_n(n_df: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), sharey=False)
    for ax, cohort in zip(axes, ["GSE131907", "GSE241934_IIT"]):
        sub = n_df[n_df.cohort == cohort].sort_values("patient")
        if sub.empty:
            ax.set_title(cohort + " (empty)")
            continue
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["n_high"], width=0.2, label="CLDN4-high mal", color="#C44E52")
        ax.bar(x, sub["n_low"], width=0.2, label="CLDN4-low mal", color="#4C72B0")
        ax.bar(x + 0.2, sub["n_B"], width=0.2, label="B", color="#55A868")
        ax.bar(x + 0.4, sub["n_TNK"], width=0.2, label="T/NK", color="#DD8452")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient"], rotation=80, ha="right", fontsize=7)
        ax.set_title(f"{cohort} (patient unit)")
        ax.set_ylabel("cells")
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_top(df: pd.DataFrame, path: Path, title: str, k: int = 20) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    if df is None or df.empty:
        ax.text(0.5, 0.5, "no pairs", ha="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    s = df.copy()
    s["absd"] = s["median_delta"].abs()
    s = s.sort_values("absd", ascending=False).head(k).sort_values("median_delta")
    colors = ["#C44E52" if d > 0 else "#4C72B0" for d in s["median_delta"]]
    labels = [f"{a}–{b} ({c})" for a, b, c in zip(s["ligand"], s["receptor"], s["ligand_class"])]
    ax.barh(range(len(s)), s["median_delta"], color=colors)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("median patient Δ (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_focus(df: pd.DataFrame, path: Path, title: str) -> None:
    want = df[df["ligand_class"].str.contains("recruit|inhibitory|barrier|attack", na=False)].copy()
    plot_top(want, path, title, k=24)


def write_finding(out: Path, given: dict, n_df: pd.DataFrame, lr: pd.DataFrame, summary: dict) -> None:
    def md(df: pd.DataFrame, cols: list[str], n=12) -> list[str]:
        if df is None or df.empty:
            return ["_(empty)_", ""]
        show = df.head(n)
        lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
        for rec in show.itertuples(index=False):
            row = []
            for c in cols:
                v = getattr(rec, c)
                if isinstance(v, float) and pd.notna(v):
                    row.append(f"{v:.3g}")
                else:
                    row.append(str(v))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
        return lines

    n_combo = int(given["n"])
    lines = [
        "# FINDING — CLDN4-only high-end LR on the given B-frac combo",
        "",
        "Additive. **CLDN4 only** (no dual-high TACSTD2∩CLDN4). Patient is the unit.",
        "The B-fraction combo that already differs is **taken as given** and is not re-audited:",
        "",
        f"- PR #290 `methods/scrna_cldn4_combo`, family `author/b/mean`",
        f"- **GSE131907 + GSE241934 IIT · k=2 · n={n_combo} · ρ=−0.513 · p=0.00388 · I²=0%**",
        "- Members: GSE131907 author-malignant n=21 (ρ=−0.465) + GSE241934 IIT author-epi n=11 (ρ=−0.609)",
        "- GSE207422-only B-frac (PR #400) is NS and is **not** re-run.",
        "",
        "This slice asks a different question: which **outgoing ligand–receptor** pairs from",
        "**CLDN4-high vs CLDN4-low malignant cells** go to **B / TLS-like** and to **T/NK**",
        "on those same two public cohorts.",
        "",
        "## Honest n (LR, not the given Spearman)",
        "",
        "The given n=32 is the B-fraction Spearman unit. LR n is the number of patients",
        "who pass the within-patient CLDN4 split and destination floors. Cells are not n.",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| Given combo patients (do not re-audit) | {n_combo} | 21 + 11 |",
        f"| GSE131907 patients in matrix (author-malig samples) | {summary['n_patients_gse131907']} | sample≈patient on the 21 mets |",
        f"| GSE241934 IIT patients in matrix | {summary['n_patients_iit']} | author Epi |",
        f"| Patients with ≥{MIN_MAL} high and ≥{MIN_MAL} low malignant (median) | {summary['n_split_median']} | CLDN4-only median |",
        f"| Eligible vs B (median) | {summary['n_B_median']} | ≥{MIN_B} B cells |",
        f"| Eligible vs TLS-like (median) | {summary['n_TLS_median']} | B + CXCL13+ T/NK; ≥{MIN_TLS} |",
        f"| Eligible vs T/NK (median) | {summary['n_TNK_median']} | ≥{MIN_TNK} T/NK |",
        f"| Patients with ≥{MIN_MAL} high/low (tertile extra) | {summary['n_split_tertile']} | high-end extra |",
        "",
        "TLS-like is a dissociated cellular proxy (B + CXCL13+ T/NK), not a histologic follicle.",
        "",
        "### Per-patient floors (median split)",
        "",
    ]
    show_n = n_df[n_df.rule == "median"][
        ["cohort", "patient", "n_mal", "n_high", "n_low", "n_B", "n_TLS", "n_TNK",
         "eligible_B", "eligible_TLS", "eligible_TNK"]
    ]
    lines += md(show_n, list(show_n.columns), n=40)

    def block(method: str, dest: str, rule: str, title: str) -> None:
        sub = lr[(lr.method == method) & (lr.destination == dest) & (lr.rule == rule)].copy()
        lines.append(f"## {title}")
        lines.append("")
        if sub.empty:
            lines.append("No pairs with ≥3 patients passing `expr_prop` on either arm.")
            lines.append("")
            return
        n_pat = int(sub["n_patients"].max()) if len(sub) else 0
        n_sig = int(((sub["padj"] < 0.05) & sub["padj"].notna()).sum())
        n_neg = int((sub["median_delta"] < 0).sum())
        n_pos = int((sub["median_delta"] > 0).sum())
        lines.append(
            f"Honest n = **{n_pat} patients**. Pairs scored = {len(sub)}. "
            f"median Δ<0: {n_neg}; Δ>0: {n_pos}; BH-FDR<0.05: {n_sig}."
        )
        lines.append("")
        sig = sub[sub["padj"].notna() & (sub["padj"] < 0.05)].sort_values("median_delta")
        if len(sig):
            lines.append("FDR < 0.05 (sorted by median Δ):")
            lines.append("")
            lines.extend(md(sig, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=20))
        else:
            lines.append("No pair at BH-FDR < 0.05. Ranked by |median Δ| (descriptive):")
            lines.append("")
            top = sub.assign(absd=sub["median_delta"].abs()).sort_values("absd", ascending=False)
            lines.extend(md(top, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=15))
        foc = sub[sub["ligand_class"].str.contains("recruit|inhibitory|barrier|attack", na=False)]
        if len(foc):
            lines.append("Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):")
            lines.append("")
            foc = foc.assign(absd=foc["median_delta"].abs()).sort_values("absd", ascending=False)
            lines.extend(md(foc, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=16))

    block("cellchat", "B", "median", "CellChat-style outgoing CLDN4-high → B (median, patient unit)")
    block("cellchat", "TLS", "median", "CellChat-style outgoing CLDN4-high → TLS-like (median)")
    block("cellchat", "TNK", "median", "CellChat-style outgoing CLDN4-high → T/NK (median)")
    block("liana_cpdb", "B", "median", "LIANA/CellPhoneDB-style outgoing CLDN4-high → B (median)")
    block("liana_cpdb", "TLS", "median", "LIANA/CellPhoneDB-style outgoing CLDN4-high → TLS-like (median)")
    block("liana_cpdb", "TNK", "median", "LIANA/CellPhoneDB-style outgoing CLDN4-high → T/NK (median)")
    block("cellchat", "B", "tertile", "Extra: CellChat tertile (high-end vs low-end) → B")
    block("cellchat", "TNK", "tertile", "Extra: CellChat tertile → T/NK")

    lines += [
        "## Methods (short)",
        "",
        "- Public processed GEO only. Matrices are **not** concatenated across accessions.",
        "- CLDN4-high / low is a **within-patient** split of malignant (or IIT author-epi) `log1p(CP10k)` CLDN4.",
        "  Median is primary. Tertile is the high-end extra (middle third dropped).",
        "- CellChat-style: 10% truncated mean, Hill \(K_h=0.5\), CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`.",
        "- LIANA-style: CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`.",
        f"  LIANA package status: `{summary.get('liana_package', 'not_run')}`.",
        "- Inference: paired Wilcoxon on patient scores (high vs low). BH-FDR within method × destination × rule.",
        "- TACSTD2 is a companion gene, never a gate.",
        "",
        "## What is not claimed",
        "",
        "- The given B-fraction ρ=−0.513 is not re-computed here.",
        "- This is not GSE207422 and not an ICI/MPR contrast on GSE131907 (treatment-naive atlas).",
        "- TLS-like is not a follicle. CopyKAT was not re-run.",
        "- Cell-pooled permutation p-values are not the inferential unit; the patient is.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/bfrac32_cldn4_hiend/scripts/download.py --outdir /tmp/bfrac32_cldn4_hiend",
        "python3 methods/bfrac32_cldn4_hiend/scripts/analyze.py --data /tmp/bfrac32_cldn4_hiend --out methods/bfrac32_cldn4_hiend",
        "```",
        "",
        "Primary table: [`results/lr_table.tsv`](results/lr_table.tsv).",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("/tmp/bfrac32_cldn4_hiend"))
    ap.add_argument("--out", type=Path, default=HERE)
    args = ap.parse_args()
    res = args.out / "results"
    fig = args.out / "figures"
    res.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    given = json.loads((HERE / "given_combo.json").read_text())
    cc = load_cellchat_lr(HERE / "db")
    cp = load_cpdb_lr(HERE / "resources")
    genes = wanted_genes(cc, cp)
    log(f"CellChat pairs={len(cc)} CPDB pairs={len(cp)} genes={len(genes)}")

    obj1 = load_gse131907(args.data, genes)
    obj2 = load_gse241934_iit(args.data, genes)

    n_rows = []
    raw_chunks = []
    for obj in (obj1, obj2):
        for pid in sorted(pd.Series(obj["patient"]).unique()):
            if pid in {"unknown", "NA", ""}:
                continue
            for rule in ("median", "tertile"):
                chunks, rec = score_patient(obj, pid, cc, cp, rule)
                n_rows.append(rec)
                raw_chunks.extend(chunks)
                log(
                    f"  {obj['cohort']} {pid} {rule} high={rec['n_high']} low={rec['n_low']} "
                    f"B={rec['n_B']} TLS={rec['n_TLS']} TNK={rec['n_TNK']} "
                    f"eligB={rec['eligible_B']} eligTNK={rec['eligible_TNK']}"
                )

    n_df = pd.DataFrame(n_rows)
    n_df.to_csv(res / "patient_n.tsv", sep="\t", index=False)
    raw = pd.concat(raw_chunks, ignore_index=True) if raw_chunks else pd.DataFrame()
    if len(raw):
        # keep a compact patient-level dump of detected pairs only
        slim = raw[raw["pass_either"]].copy()
        slim.to_csv(res / "patient_outgoing.tsv", sep="\t", index=False)
    else:
        slim = raw
        raw.to_csv(res / "patient_outgoing.tsv", sep="\t", index=False)

    lr = aggregate_patients(raw)
    lr.to_csv(res / "lr_table.tsv", sep="\t", index=False)
    if len(lr):
        lr[lr.method == "cellchat"].to_csv(res / "lr_table_cellchat.tsv", sep="\t", index=False)
        lr[lr.method == "liana_cpdb"].to_csv(res / "lr_table_liana.tsv", sep="\t", index=False)
        for method, dest, rule in (
            ("cellchat", "B", "median"),
            ("cellchat", "TLS", "median"),
            ("cellchat", "TNK", "median"),
            ("liana_cpdb", "B", "median"),
            ("liana_cpdb", "TLS", "median"),
            ("liana_cpdb", "TNK", "median"),
            ("cellchat", "B", "tertile"),
            ("cellchat", "TNK", "tertile"),
        ):
            sub = lr[(lr.method == method) & (lr.destination == dest) & (lr.rule == rule)]
            tag = f"{method}_{dest}_{rule}"
            plot_top(sub, fig / f"fig_top_{tag}.png", f"{method} {dest} {rule}: median patient Δ")
            plot_focus(sub, fig / f"fig_focus_{tag}.png", f"{method} {dest} {rule}: barrier/inhib/recruit/attack")

    plot_n(n_df[n_df.rule == "median"], fig / "fig_n_per_patient.png")

    def n_elig(rule, col):
        sub = n_df[n_df.rule == rule]
        return int(sub[col].sum()) if len(sub) and col in sub.columns else 0

    liana_pkg = "not_run"
    try:
        import liana  # noqa: F401

        liana_pkg = "imported_but_not_used_patient_cpdb_is_primary"
    except Exception as exc:
        liana_pkg = f"not_run ({type(exc).__name__})"

    summary = {
        "given_combo_n": given["n"],
        "given_combo_rho": given["rho"],
        "n_patients_gse131907": int(n_df[n_df.cohort == "GSE131907"].patient.nunique()) if len(n_df) else 0,
        "n_patients_iit": int(n_df[n_df.cohort == "GSE241934_IIT"].patient.nunique()) if len(n_df) else 0,
        "n_split_median": n_elig("median", "eligible_mal_split"),
        "n_split_tertile": n_elig("tertile", "eligible_mal_split"),
        "n_B_median": n_elig("median", "eligible_B"),
        "n_TLS_median": n_elig("median", "eligible_TLS"),
        "n_TNK_median": n_elig("median", "eligible_TNK"),
        "n_lr_rows": int(len(lr)),
        "n_patient_pair_rows": int(len(slim)),
        "cellchat_pairs": int(len(cc)),
        "cpdb_pairs": int(len(cp)),
        "liana_package": liana_pkg,
        "unit": "patient",
        "cldn4_only": True,
        "dual_high": False,
        "gse207422_included": False,
        "gse131907_mal_def": obj1["mal_def"],
        "iit_mal_def": obj2["mal_def"],
        "iit_author_majors": obj2.get("author_majors", []),
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(args.out, given, n_df, lr, summary)
    # also copy FINDING to results for the "PR has LR table" readers
    (res / "FINDING.md").write_text((args.out / "FINDING.md").read_text())
    log(json.dumps(summary, indent=2))
    log(f"[done] lr_table rows={len(lr)} -> {res / 'lr_table.tsv'}")


if __name__ == "__main__":
    main()
