#!/usr/bin/env python3
"""CLDN4-only high-end LR on the given differing pair GSE123902 + GSE189357.

The PR #459 %pos Spearman (n=22, ρ=−0.638) is cited, not re-audited.
No dual-high. Patient/donor is the unit.

CellChat-style: Jin 2021 Hill probability (10% truncated mean, Kh=0.5,
CellChatDB v2 protein pairs, expr_prop >= 0.10).

LIANA-style extra: CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5
pairs, expr_prop >= 0.10.

Outgoing only: CLDN4-high marker-malignant → same-unit T/NK.
Honest n may be < 22. Between-patient Q4 tails (7/5) are thin.
"""
from __future__ import annotations

import argparse
import gzip
import json
import tarfile
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
MIN_TNK = 20
MIN_MAL_Q4 = 40
THIN_ARM = 20

EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]

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
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "NKG7"] + EPI + TNK


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
    return {g.upper() for g in genes}


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


def log1p_cp10k(umi: dict[str, np.ndarray], lib: np.ndarray) -> dict[str, np.ndarray]:
    scale = np.where(lib > 0, 1e4 / lib, 0.0).astype(np.float32)
    return {g: np.log1p(arr * scale).astype(np.float32) for g, arr in umi.items()}


def marker_masks(umi: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def col(g: str) -> np.ndarray:
        return umi[g] if g in umi else np.zeros(n, dtype=np.float32)

    epi = np.zeros(n, dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    ptprc = col("PTPRC")
    mal = epi & (ptprc == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in TNK:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    return mal, tnk


def _tissue_from_name(name: str) -> str:
    if "NORMAL" in name:
        return "NORMAL"
    if "METASTASIS" in name:
        return "METASTASIS"
    if "PRIMARY" in name:
        return "PRIMARY"
    return "OTHER"


def stream_dense_csv_gz(handle, genes: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """cells × genes dense CSV. col0 = cell id. Returns barcodes, umi, lib."""
    header = handle.readline().decode().strip().split(",")
    symbols = [g.split(".")[0].upper() for g in header[1:]]
    keep_idx = [i for i, g in enumerate(symbols) if g in genes]
    keep_names = [symbols[i] for i in keep_idx]
    cells: list[str] = []
    libs: list[float] = []
    cols: dict[str, list[float]] = {g: [] for g in keep_names}
    n_line = 0
    for raw in handle:
        line = raw.decode().strip()
        if not line:
            continue
        cell, _, rest = line.partition(",")
        vals = np.fromstring(rest, sep=",", dtype=np.float32)
        if vals.size != len(symbols):
            continue
        cells.append(cell)
        libs.append(float(vals.sum()))
        for g, i in zip(keep_names, keep_idx):
            cols[g].append(float(vals[i]))
        n_line += 1
        if n_line % 2000 == 0:
            log(f"    dense rows={n_line}")
    n = len(cells)
    umi = {g: np.asarray(v, dtype=np.float32) for g, v in cols.items()}
    # if duplicate symbols, sum
    merged: dict[str, np.ndarray] = {}
    for g, arr in umi.items():
        merged[g] = arr if g not in merged else merged[g] + arr
    return np.asarray(cells, dtype=object), merged, np.asarray(libs, dtype=np.float32)


def load_gse123902(tar_path: Path, genes: set[str]) -> dict:
    """Tumor/met donors only. Prefer PRIMARY over METASTASIS. Drop NORMAL."""
    rows = []
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".csv.gz")]
        log(f"[GSE123902] members={len(members)}")
        for m in members:
            name = Path(m.name).name
            parts = name.split("_")
            gsm = parts[0]
            donor = parts[2] if len(parts) > 2 else name
            tissue = _tissue_from_name(name)
            if tissue == "NORMAL":
                log(f"  skip NORMAL {name}")
                continue
            log(f"  {name} donor={donor} tissue={tissue}")
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            cells, umi, lib = stream_dense_csv_gz(raw, genes)
            mal, tnk = marker_masks(umi, len(cells))
            rows.append(
                {
                    "donor": donor,
                    "gsm": gsm,
                    "tissue": tissue,
                    "file": name,
                    "cells": cells,
                    "umi": umi,
                    "lib": lib,
                    "is_mal": mal,
                    "is_tnk": tnk,
                }
            )
            log(
                f"    n={len(cells)} mal={int(mal.sum())} tnk={int(tnk.sum())} "
                f"CLDN4={'yes' if 'CLDN4' in umi else 'NO'}"
            )

    # one donor: PRIMARY over METASTASIS
    by_donor: dict[str, list] = {}
    for rec in rows:
        by_donor.setdefault(rec["donor"], []).append(rec)
    keep = []
    for donor, recs in sorted(by_donor.items()):
        recs = sorted(recs, key=lambda r: 0 if r["tissue"] == "PRIMARY" else 1)
        keep.append(recs[0])

    # concatenate kept libraries
    cells = np.concatenate([r["cells"] for r in keep]) if keep else np.array([], dtype=object)
    patient = np.concatenate([np.array([r["donor"]] * len(r["cells"]), dtype=object) for r in keep]) if keep else np.array([], dtype=object)
    tissue = np.concatenate([np.array([r["tissue"]] * len(r["cells"]), dtype=object) for r in keep]) if keep else np.array([], dtype=object)
    lib = np.concatenate([r["lib"] for r in keep]) if keep else np.array([], dtype=np.float32)
    is_mal = np.concatenate([r["is_mal"] for r in keep]) if keep else np.array([], dtype=bool)
    is_tnk = np.concatenate([r["is_tnk"] for r in keep]) if keep else np.array([], dtype=bool)
    all_genes = set()
    for r in keep:
        all_genes.update(r["umi"])
    umi = {}
    for g in all_genes:
        chunks = []
        for r in keep:
            if g in r["umi"]:
                chunks.append(r["umi"][g])
            else:
                chunks.append(np.zeros(len(r["cells"]), dtype=np.float32))
        umi[g] = np.concatenate(chunks)
    logx = log1p_cp10k(umi, lib)
    pos = {g: (umi[g] > 0).astype(np.float32) for g in umi}
    log(
        f"[GSE123902] donors={len(keep)} cells={len(cells)} mal={int(is_mal.sum())} "
        f"TNK={int(is_tnk.sum())}"
    )
    return {
        "cohort": "GSE123902",
        "cells": cells,
        "patient": patient,
        "tissue": tissue,
        "is_mal": is_mal,
        "is_tnk": is_tnk,
        "logx": logx,
        "pos": pos,
        "cldn4": logx.get("CLDN4", np.zeros(len(cells), dtype=np.float32)),
        "mal_def": "marker_malig (EPCAM|KRT8|KRT18|KRT19)>0 & PTPRC==0; tumor/met donors",
        "unit": "donor",
        "n_units_matrix": int(pd.Series(patient).nunique()) if len(patient) else 0,
    }


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_features" in n or f"_{sample}_genes" in n)
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            p = line.decode().strip().split("\t")
            genes.append((p[1] if len(p) > 1 else p[0]).split(".")[0].upper())
    return genes


def _read_10x_barcodes(tf: tarfile.TarFile, members: dict, sample: str) -> np.ndarray:
    name = next(n for n in members if f"_{sample}_barcodes" in n)
    bcs = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            bcs.append(line.decode().strip().split("\t")[0])
    return np.asarray(bcs, dtype=object)


def stream_sample_mtx(tf, members, sample: str, genes: set[str]):
    symbols = _read_10x_features(tf, members, sample)
    bcs = _read_10x_barcodes(tf, members, sample)
    n = len(bcs)
    want_rows = {i for i, s in enumerate(symbols) if s in genes}
    keep: dict[str, np.ndarray] = {}
    for i in want_rows:
        keep.setdefault(symbols[i], np.zeros(n, dtype=np.float32))
    lib = np.zeros(n, dtype=np.float64)
    mtx_name = next(n_ for n_ in members if f"_{sample}_matrix.mtx" in n_)
    log(f"  {sample} cells={n} want_rows={len(want_rows)}")
    with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as f:
        for line in f:
            if line.decode().startswith("%"):
                continue
            break
        n_lines = 0
        for line in f:
            parts = line.decode().split()
            if len(parts) < 3:
                continue
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            v = float(parts[2])
            lib[c] += v
            if r in want_rows:
                keep[symbols[r]][c] += v
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                log(f"    mtx entries={n_lines}")
    return bcs, keep, lib.astype(np.float32)


def load_gse189357(tar_path: Path, genes: set[str]) -> dict:
    chunks = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        samples = [f"TD{i}" for i in range(1, 10)]
        for sample in samples:
            cells, umi, lib = stream_sample_mtx(tf, members, sample, genes)
            mal, tnk = marker_masks(umi, len(cells))
            chunks.append(
                {
                    "patient": sample,
                    "cells": cells,
                    "umi": umi,
                    "lib": lib,
                    "is_mal": mal,
                    "is_tnk": tnk,
                }
            )
            log(f"    {sample} n={len(cells)} mal={int(mal.sum())} tnk={int(tnk.sum())}")

    cells = np.concatenate([c["cells"] for c in chunks])
    patient = np.concatenate([np.array([c["patient"]] * len(c["cells"]), dtype=object) for c in chunks])
    lib = np.concatenate([c["lib"] for c in chunks])
    is_mal = np.concatenate([c["is_mal"] for c in chunks])
    is_tnk = np.concatenate([c["is_tnk"] for c in chunks])
    all_genes = set()
    for c in chunks:
        all_genes.update(c["umi"])
    umi = {}
    for g in all_genes:
        parts = []
        for c in chunks:
            parts.append(c["umi"][g] if g in c["umi"] else np.zeros(len(c["cells"]), dtype=np.float32))
        umi[g] = np.concatenate(parts)
    logx = log1p_cp10k(umi, lib)
    pos = {g: (umi[g] > 0).astype(np.float32) for g in umi}
    log(f"[GSE189357] patients=9 cells={len(cells)} mal={int(is_mal.sum())} TNK={int(is_tnk.sum())}")
    return {
        "cohort": "GSE189357",
        "cells": cells,
        "patient": patient,
        "tissue": np.array(["TUMOR"] * len(cells), dtype=object),
        "is_mal": is_mal,
        "is_tnk": is_tnk,
        "logx": logx,
        "pos": pos,
        "cldn4": logx.get("CLDN4", np.zeros(len(cells), dtype=np.float32)),
        "mal_def": "marker_malig (EPCAM|KRT8|KRT18|KRT19)>0 & PTPRC==0; TD1–TD9",
        "unit": "patient",
        "n_units_matrix": 9,
    }


def split_high_low(cldn4: np.ndarray, mal: np.ndarray, rule: str) -> tuple[np.ndarray, np.ndarray, dict]:
    vals = cldn4[mal]
    hi = np.zeros(len(cldn4), dtype=bool)
    lo = np.zeros(len(cldn4), dtype=bool)
    info = {"rule": rule, "n_mal": int(mal.sum())}
    if rule == "median":
        if vals.size < MIN_MAL * 2:
            return hi, lo, info
        thr = float(np.median(vals))
        hi[mal & (cldn4 > thr)] = True
        lo[mal & (cldn4 <= thr)] = True
        info["threshold"] = thr
    elif rule == "tertile":
        if vals.size < MIN_MAL * 3:
            return hi, lo, info
        q1, q2 = np.quantile(vals, [1 / 3, 2 / 3])
        hi[mal & (cldn4 >= q2)] = True
        lo[mal & (cldn4 <= q1)] = True
        info["q_low"] = float(q1)
        info["q_high"] = float(q2)
    elif rule == "q4q1":
        if vals.size < MIN_MAL_Q4:
            return hi, lo, info
        q1, q3 = np.quantile(vals, [0.25, 0.75])
        hi[mal & (cldn4 >= q3)] = True
        lo[mal & (cldn4 <= q1)] = True
        info["q_low"] = float(q1)
        info["q_high"] = float(q3)
    else:
        raise ValueError(rule)
    info["n_high"] = int(hi.sum())
    info["n_low"] = int(lo.sum())
    return hi, lo, info


def score_patient(obj: dict, pid: str, pairs_cc: pd.DataFrame, pairs_cp: pd.DataFrame, rule: str):
    m = obj["patient"] == pid
    mal = obj["is_mal"] & m
    tnk = obj["is_tnk"] & m
    hi, lo, info = split_high_low(obj["cldn4"], mal, rule)
    n_high, n_low, n_tnk = int(hi.sum()), int(lo.sum()), int(tnk.sum())
    thin = (n_high < THIN_ARM) or (n_low < THIN_ARM)
    rec = {
        "cohort": obj["cohort"],
        "patient": pid,
        "unit": obj["unit"],
        "rule": rule,
        "n_mal": int(mal.sum()),
        "n_high": n_high,
        "n_low": n_low,
        "n_TNK": n_tnk,
        "mean_CLDN4_high": float(obj["cldn4"][hi].mean()) if hi.any() else np.nan,
        "mean_CLDN4_low": float(obj["cldn4"][lo].mean()) if lo.any() else np.nan,
        "eligible_mal_split": n_high >= MIN_MAL and n_low >= MIN_MAL,
        "eligible_TNK": n_high >= MIN_MAL and n_low >= MIN_MAL and n_tnk >= MIN_TNK,
        "thin_tail": bool(thin and n_high >= MIN_MAL and n_low >= MIN_MAL),
        "given_combo_member": True,
    }
    rows = []
    if not rec["eligible_TNK"]:
        return rows, rec
    hi_s = group_stats(obj["logx"], obj["pos"], hi)
    lo_s = group_stats(obj["logx"], obj["pos"], lo)
    recv = group_stats(obj["logx"], obj["pos"], tnk)
    for method, pairs in (("cellchat", pairs_cc), ("liana_cpdb", pairs_cp)):
        a = score_pairs(pairs, hi_s, recv, "cellchat" if method == "cellchat" else "liana")
        b = score_pairs(pairs, lo_s, recv, "cellchat" if method == "cellchat" else "liana")
        if a.empty or b.empty:
            continue
        merged = a.merge(
            b,
            on=[
                "interaction_name", "pathway_name", "annotation", "ligand", "receptor",
                "ligand_genes", "receptor_genes", "ligand_class", "pair_origin", "method",
            ],
            suffixes=("_high", "_low"),
        )
        merged["delta"] = merged["score_high"] - merged["score_low"]
        merged["pass_either"] = merged["detected_high"] | merged["detected_low"]
        merged["method"] = method
        merged["cohort"] = obj["cohort"]
        merged["patient"] = pid
        merged["destination"] = "TNK"
        merged["rule"] = rule
        merged["thin_tail"] = rec["thin_tail"]
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
    keys = [
        "method", "destination", "rule", "interaction_name", "pathway_name",
        "ligand", "receptor", "ligand_genes", "receptor_genes", "ligand_class", "pair_origin",
    ]
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
                "n_thin_tail": int(keep.loc[keep["thin_tail"] == True, "patient"].nunique())
                if "thin_tail" in keep.columns
                else 0,
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
    for (_, _, _), idx in out.groupby(["method", "destination", "rule"]).groups.items():
        mask = out.index.isin(idx) & out["pval"].notna()
        if mask.sum() == 0:
            continue
        out.loc[mask, "padj"] = multipletests(out.loc[mask, "pval"], method="fdr_bh")[1]
    return out.sort_values(["method", "destination", "rule", "median_delta"])


def plot_n(n_df: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), sharey=False)
    for ax, cohort in zip(axes, ["GSE123902", "GSE189357"]):
        sub = n_df[n_df.cohort == cohort].sort_values("patient")
        if sub.empty:
            ax.set_title(cohort + " (empty)")
            continue
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["n_high"], width=0.2, label="CLDN4-high mal", color="#C44E52")
        ax.bar(x, sub["n_low"], width=0.2, label="CLDN4-low mal", color="#4C72B0")
        ax.bar(x + 0.2, sub["n_TNK"], width=0.2, label="T/NK", color="#DD8452")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient"], rotation=80, ha="right", fontsize=7)
        ax.set_title(f"{cohort} (patient/donor unit)")
        ax.set_ylabel("cells")
        ax.legend(frameon=False, fontsize=7)
        for i, thin in enumerate(sub["thin_tail"].tolist()):
            if thin:
                ax.text(i, max(sub["n_TNK"].max(), sub["n_high"].max()) * 0.95, "thin", ha="center", fontsize=6, color="#7a2d0b")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_honest_n(n_df: pd.DataFrame, given: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    labels = [
        "Given combo\n(do not re-audit)",
        "Given Q4 vs Q1\n(7/5 thin)",
        "Median split\n→ T/NK",
        "Tertile\n→ T/NK",
        "Q4 vs Q1 within\n→ T/NK",
        "Median thin-tail\nunits",
    ]
    med = n_df[n_df.rule == "median"]
    ter = n_df[n_df.rule == "tertile"]
    q4 = n_df[n_df.rule == "q4q1"]
    vals = [
        int(given["n"]),
        int(given.get("n_q1", 7) + given.get("n_q4", 5)),
        int(med["eligible_TNK"].sum()) if len(med) else 0,
        int(ter["eligible_TNK"].sum()) if len(ter) else 0,
        int(q4["eligible_TNK"].sum()) if len(q4) else 0,
        int(med["thin_tail"].sum()) if len(med) else 0,
    ]
    colors = ["#4a4a4a", "#C44E52", "#4C72B0", "#55A868", "#8172B3", "#CCB974"]
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("n (patients / donors)")
    ax.set_title("Honest n — given Spearman vs CellChat floors (tails may be thin)")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.15, str(v), ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_top(df: pd.DataFrame, path: Path, title: str, k: int = 20) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.6))
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
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_focus(df: pd.DataFrame, path: Path, title: str) -> None:
    want = df[df["ligand_class"].str.contains("recruit|inhibitory|barrier|attack", na=False)].copy()
    plot_top(want, path, title, k=24)


def md_table(df: pd.DataFrame, cols: list[str], n: int = 12) -> list[str]:
    if df is None or df.empty:
        return ["_(empty)_", ""]
    show = df.head(n)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, rec in show.iterrows():
        row = []
        for c in cols:
            v = rec[c]
            if isinstance(v, float) and pd.notna(v):
                row.append(f"{v:.3g}" if abs(v) < 1 or abs(v) >= 10 else f"{v:.3f}")
            else:
                row.append(str(v))
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    return lines


def write_finding(out: Path, given: dict, n_df: pd.DataFrame, lr: pd.DataFrame, summary: dict) -> None:
    med = n_df[n_df.rule == "median"] if len(n_df) else n_df
    ter = n_df[n_df.rule == "tertile"] if len(n_df) else n_df
    q4 = n_df[n_df.rule == "q4q1"] if len(n_df) else n_df
    n_med = int(med["eligible_TNK"].sum()) if len(med) else 0
    n_ter = int(ter["eligible_TNK"].sum()) if len(ter) else 0
    n_q4 = int(q4["eligible_TNK"].sum()) if len(q4) else 0
    n_thin = int(med["thin_tail"].sum()) if len(med) else 0
    n_thin_q4 = int(q4["thin_tail"].sum()) if len(q4) else 0

    cc_med = lr[(lr.method == "cellchat") & (lr.rule == "median") & (lr.destination == "TNK")] if len(lr) else lr
    n_sig = int(((cc_med["padj"] < 0.05) & cc_med["padj"].notna()).sum()) if len(cc_med) else 0
    n_neg = int((cc_med["median_delta"] < 0).sum()) if len(cc_med) else 0
    n_pos = int((cc_med["median_delta"] > 0).sum()) if len(cc_med) else 0

    lines = [
        "# FINDING — CLDN4-only high-end CellChat on GSE123902 + GSE189357",
        "",
        "ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit.",
        "",
        "The pair that already **differs** is **taken as given** and is not re-audited:",
        "",
        f"- PR #459 `methods/scrna_cldn4_combo_enum`, pair **GSE123902+GSE189357 %pos**",
        f"- **n=22 · ρ=−0.638 · p=0.003 · I²=0%**",
        f"- Members: GSE123902 marker-malignant donors n=13 + GSE189357 marker-malignant patients n=9",
        f"- Between-patient Q4 vs Q1 on that same %pos vector is **r=−1.000 p=0.003 on tails 7/5 — thin**. That is not the CellChat n.",
        "",
        "This slice asks a different question: which **outgoing ligand–receptor** pairs from",
        "**CLDN4-high vs CLDN4-low marker-malignant cells** go to **same-unit T/NK**.",
        "",
        "## Honest n (LR, not the given Spearman)",
        "",
        "The given n=22 is the %pos Spearman unit. LR n is the number of patients/donors",
        "who pass the within-unit CLDN4 split and T/NK floors. Cells are not n.",
        "Tails may be thin — said here, not hidden.",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| Given combo (do not re-audit) | {given['n']} | 13 + 9; %pos ρ=−0.638 |",
        f"| Given between-patient Q4 vs Q1 | {given.get('n_q1', 7)}/{given.get('n_q4', 5)} | **thin tails**; r=−1 is not CellChat n |",
        f"| Eligible vs T/NK (median) | **{n_med}** | ≥10 high/low mal + ≥20 T/NK |",
        f"| Eligible vs T/NK (tertile extra) | {n_ter} | middle third dropped |",
        f"| Eligible vs T/NK (Q4 vs Q1 extra) | {n_q4} | ≥40 mal; high-end |",
        f"| Median units with thin arm (<20 high or low) | {n_thin} | kept and flagged |",
        f"| Q4 vs Q1 units with thin arm | {n_thin_q4} | expected on GSE123902 mets |",
        "",
        "GSE123902 LX699 (46 malignant) and LX701 (90) are the thinnest tumor/met donors.",
        "Q4 vs Q1 within those units is a thin-tail extra, not the primary n.",
        "",
        "### Per-patient floors (median split)",
        "",
    ]
    show_n = med[["cohort", "patient", "n_mal", "n_high", "n_low", "n_TNK", "eligible_TNK", "thin_tail"]] if len(med) else med
    lines += md_table(show_n, list(show_n.columns), n=40) if len(show_n) else ["_(empty)_", ""]

    def block(method: str, rule: str, title: str) -> None:
        sub = lr[(lr.method == method) & (lr.destination == "TNK") & (lr.rule == rule)].copy() if len(lr) else lr
        lines.append(f"## {title}")
        lines.append("")
        if sub is None or sub.empty:
            lines.append("No pairs with ≥3 patients passing `expr_prop` on either arm.")
            lines.append("")
            return
        n_pat = int(sub["n_patients"].max()) if len(sub) else 0
        n_sig_b = int(((sub["padj"] < 0.05) & sub["padj"].notna()).sum())
        n_neg_b = int((sub["median_delta"] < 0).sum())
        n_pos_b = int((sub["median_delta"] > 0).sum())
        lines.append(
            f"Honest n = **{n_pat} patients/donors**. Pairs scored = {len(sub)}. "
            f"median Δ<0: {n_neg_b}; Δ>0: {n_pos_b}; BH-FDR<0.05: {n_sig_b}."
        )
        lines.append("")
        sig = sub[sub["padj"].notna() & (sub["padj"] < 0.05)].sort_values("median_delta")
        if len(sig):
            lines.append("FDR < 0.05 (sorted by median Δ):")
            lines.append("")
            lines.extend(md_table(sig, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=20))
        else:
            lines.append("No pair at BH-FDR < 0.05. Ranked by |median Δ| (descriptive):")
            lines.append("")
            top = sub.assign(absd=sub["median_delta"].abs()).sort_values("absd", ascending=False)
            lines.extend(md_table(top, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=15))
        foc = sub[sub["ligand_class"].str.contains("recruit|inhibitory|barrier|attack", na=False)]
        if len(foc):
            lines.append("Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):")
            lines.append("")
            foc = foc.assign(absd=foc["median_delta"].abs()).sort_values("absd", ascending=False)
            lines.extend(md_table(foc, ["ligand", "receptor", "ligand_class", "n_patients", "median_delta", "pval", "padj"], n=16))

    block("cellchat", "median", "CellChat-style outgoing CLDN4-high → T/NK (median, patient unit)")
    block("liana_cpdb", "median", "Extra: LIANA/CellPhoneDB-style outgoing CLDN4-high → T/NK (median)")
    block("cellchat", "tertile", "Extra: CellChat tertile (high-end vs low-end) → T/NK")
    block("cellchat", "q4q1", "Extra: CellChat within-unit Q4 vs Q1 → T/NK (thin tails flagged)")

    lines += [
        "## Methods (short)",
        "",
        "- Public processed GEO only. Matrices are **not** concatenated across accessions.",
        "- Marker-malignant and T/NK gates match PR #459 (not author labels).",
        "- CLDN4-high / low is a **within-unit** split of marker-malignant `log1p(CP10k)` CLDN4.",
        "  Median is primary. Tertile and Q4 vs Q1 are high-end extras.",
        "- CellChat-style: 10% truncated mean, Hill Kh=0.5, CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`.",
        "- LIANA-style extra: CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`.",
        f"  LIANA package status: `{summary.get('liana_package', 'not_run')}`.",
        "- Inference: paired Wilcoxon on patient scores (high vs low). BH-FDR within method × rule.",
        "- TACSTD2 is a companion gene, never a gate.",
        "",
        "## What is not claimed",
        "",
        "- The given %pos ρ=−0.638 is not re-computed here.",
        "- Between-patient Q4 vs Q1 r=−1 on 7/5 is thin and is not the LR n.",
        "- This is not a dual-high TACSTD2∩CLDN4 gate and not the 7-cohort pool.",
        "- Cell-pooled permutation p-values are not the inferential unit; the patient/donor is.",
        "",
        "## Extra figures",
        "",
        "- `figures/fig_n_per_patient.png` — per-unit high/low malignant and T/NK floors",
        "- `figures/fig_honest_n.png` — given n vs CellChat n (thin tails called out)",
        "- `figures/fig_top_cellchat_TNK_median.png` / `fig_focus_cellchat_TNK_median.png`",
        "- tertile and Q4 vs Q1 extras under `figures/`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/pair_123902_189357_hiend_cldn4/scripts/download.py",
        "python3 methods/pair_123902_189357_hiend_cldn4/scripts/analyze.py",
        "```",
        "",
        "Primary table: [`results/lr_table.tsv`](results/lr_table.tsv).",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("/tmp/geo_pair_123902_189357"))
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

    p123 = args.data / "GSE123902_RAW.tar"
    p189 = args.data / "GSE189357_RAW.tar"
    if not p123.exists() or not p189.exists():
        raise SystemExit(f"missing tars under {args.data}; run scripts/download.py")

    obj1 = load_gse123902(p123, genes)
    obj2 = load_gse189357(p189, genes)

    n_rows = []
    raw_chunks = []
    for obj in (obj1, obj2):
        for pid in sorted(pd.Series(obj["patient"]).unique()):
            if pid in {"unknown", "NA", ""}:
                continue
            for rule in ("median", "tertile", "q4q1"):
                chunks, rec = score_patient(obj, pid, cc, cp, rule)
                n_rows.append(rec)
                raw_chunks.extend(chunks)
                log(
                    f"  {obj['cohort']} {pid} {rule} high={rec['n_high']} low={rec['n_low']} "
                    f"TNK={rec['n_TNK']} elig={rec['eligible_TNK']} thin={rec['thin_tail']}"
                )

    n_df = pd.DataFrame(n_rows)
    n_df.to_csv(res / "patient_n.tsv", sep="\t", index=False)
    raw = pd.concat(raw_chunks, ignore_index=True) if raw_chunks else pd.DataFrame()
    if len(raw):
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
        focus = lr[lr["ligand_class"].str.contains("recruit|inhibitory|barrier|attack", na=False)]
        focus.to_csv(res / "lr_table_focus.tsv", sep="\t", index=False)
        for method, rule in (
            ("cellchat", "median"),
            ("cellchat", "tertile"),
            ("cellchat", "q4q1"),
            ("liana_cpdb", "median"),
        ):
            sub = lr[(lr.method == method) & (lr.destination == "TNK") & (lr.rule == rule)]
            tag = f"{method}_TNK_{rule}"
            plot_top(sub, fig / f"fig_top_{tag}.png", f"{method} T/NK {rule}: median patient Δ")
            plot_focus(sub, fig / f"fig_focus_{tag}.png", f"{method} T/NK {rule}: barrier/inhib/recruit/attack")

    plot_n(n_df[n_df.rule == "median"] if len(n_df) else n_df, fig / "fig_n_per_patient.png")
    plot_honest_n(n_df, given, fig / "fig_honest_n.png")

    liana_pkg = "not_run"
    try:
        import liana  # noqa: F401

        liana_pkg = "imported_but_not_used_patient_cpdb_is_primary"
    except Exception as exc:  # noqa: BLE001
        liana_pkg = f"not_run ({type(exc).__name__})"

    def n_elig(rule, col):
        sub = n_df[n_df.rule == rule]
        return int(sub[col].sum()) if len(sub) and col in sub.columns else 0

    summary = {
        "given_combo_n": given["n"],
        "given_combo_rho": given["rho"],
        "given_q4q1_thin": True,
        "given_q4q1_n": f"{given.get('n_q1', 7)}/{given.get('n_q4', 5)}",
        "n_patients_gse123902": int(n_df[n_df.cohort == "GSE123902"].patient.nunique()) if len(n_df) else 0,
        "n_patients_gse189357": int(n_df[n_df.cohort == "GSE189357"].patient.nunique()) if len(n_df) else 0,
        "n_TNK_median": n_elig("median", "eligible_TNK"),
        "n_TNK_tertile": n_elig("tertile", "eligible_TNK"),
        "n_TNK_q4q1": n_elig("q4q1", "eligible_TNK"),
        "n_thin_median": n_elig("median", "thin_tail"),
        "n_lr_rows": int(len(lr)),
        "n_patient_pair_rows": int(len(slim)),
        "cellchat_pairs": int(len(cc)),
        "cpdb_pairs": int(len(cp)),
        "liana_package": liana_pkg,
        "unit": "patient/donor",
        "cldn4_only": True,
        "dual_high": False,
        "destination": "TNK",
        "gse123902_mal_def": obj1["mal_def"],
        "gse189357_mal_def": obj2["mal_def"],
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(args.out, given, n_df, lr, summary)
    (res / "FINDING.md").write_text((args.out / "FINDING.md").read_text())
    log(json.dumps(summary, indent=2))
    log(f"[done] lr_table rows={len(lr)} -> {res / 'lr_table.tsv'}")


if __name__ == "__main__":
    main()
