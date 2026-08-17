#!/usr/bin/env python3
"""Pair GSE189357 + GSE205335: CellChat-style CLDN4-high malignant → T/NK.

ADDITIVE. CLDN4 only. No dual-high. Patient is the unit.
PR #459 combo (%pos n=31 ρ=−0.478, Q4 r=−0.750) is given and is not re-ranked.

Primary test: within-patient split of malignant cells (median and Q4 vs Q1),
outgoing to same-patient T/NK, Wilcoxon signed-rank on per-patient P.
Honest paired n requires both CLDN4 bins (≥10 cells) and T/NK ≥20.

CellChat R / LIANA are not run. Probability is Jin et al. 2021
(10% truncated mean, Hill K_h=0.5, CellChatDB v2 protein pairs).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
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

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_BIN = 10
MIN_TNK = 20
MIN_DETECT_PAIRED = 6
MIN_DETECT_ARM = 3
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "CD8B", "NKG7", "IFNG"] + EPI + TNK_MARKERS

CLASSICAL_MHC1 = {"HLA-A", "HLA-B", "HLA-C"}
RECRUIT_LIGANDS = {"CXCL9", "CXCL10", "CXCL16", "CCL4", "CCL5"}
KEY_NAMES = {
    "CD274_PDCD1",
    "NECTIN2_TIGIT",
    "CXCL9_CXCR3",
    "CXCL10_CXCR3",
    "CXCL16_CXCR6",
    "CCL4_CCR5",
    "CCL5_CCR5",
    "CCL5_CCR1",
    "HLA-A_CD8A",
    "HLA-B_CD8A",
    "HLA-C_CD8A",
    "HLA-A_CD8B",
    "HLA-B_CD8B",
    "HLA-C_CD8B",
    "HLA-E_CD8A",
    "HLA-E_KLRK1",
}


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


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


def pair_focus(name: str, ligand_genes: tuple[str, ...] | str, pathway: str) -> str:
    lig = set(ligand_genes.split("|") if isinstance(ligand_genes, str) else ligand_genes)
    if name == "CD274_PDCD1":
        return "CD274–PDCD1"
    if name == "NECTIN2_TIGIT":
        return "NECTIN2–TIGIT"
    if lig & RECRUIT_LIGANDS:
        return "T-recruit"
    if pathway == "MHC-I" or lig & CLASSICAL_MHC1:
        return "MHC-I"
    return "other"


def ligand_class(genes) -> str:
    parts = set(genes.split("|") if isinstance(genes, str) else genes)
    tags = []
    if parts & {"CDH1", "CLDN4", "F11R", "NECTIN2", "PVR", "EPCAM", "CEACAM5"}:
        tags.append("barrier")
    if parts & {"CD274", "LGALS9", "HLA-E", "HLA-G", "NECTIN2", "PVR", "TGFB1"}:
        tags.append("inhibitory")
    if parts & RECRUIT_LIGANDS:
        tags.append("recruit")
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


def complex_from_maps(means: dict[str, float], props: dict[str, float], subunits: tuple[str, ...]):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    return geom_mean(vals), float(min(prs)) if prs else 0.0


def compartment_gene_stats(log_cp, pos, idx: np.ndarray):
    means, props = {}, {}
    if idx.size == 0:
        return means, props
    for gene, vec in log_cp.items():
        means[gene] = trim_mean_1d(vec[idx])
        props[gene] = float(pos[gene][idx].mean())
    return means, props


def score_lr(lr: pd.DataFrame, log_cp, pos, src_idx, tgt_idx) -> list[dict]:
    rows = []
    if src_idx.size < MIN_BIN or tgt_idx.size < MIN_TNK:
        return rows
    src_mu, src_pr = compartment_gene_stats(log_cp, pos, src_idx)
    tgt_mu, tgt_pr = compartment_gene_stats(log_cp, pos, tgt_idx)
    for rec in lr.itertuples(index=False):
        lig, lig_p = complex_from_maps(src_mu, src_pr, rec.ligand_genes)
        recp, rec_p = complex_from_maps(tgt_mu, tgt_pr, rec.receptor_genes)
        detected = lig_p >= EXPR_PROP and rec_p >= EXPR_PROP
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "prob": hill_prob(lig, recp) if detected else 0.0,
                "detected": bool(detected),
                "ligand_mean": lig,
                "receptor_mean": recp,
                "ligand_prop": lig_p,
                "receptor_prop": rec_p,
            }
        )
    return rows


def split_malignant(cldn4: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray] | None:
    n = int(cldn4.size)
    if n < MIN_BIN * 2:
        return None
    if mode == "median":
        med = float(np.median(cldn4))
        if med == 0.0:
            high = np.flatnonzero(cldn4 > 0)
            low = np.flatnonzero(cldn4 == 0)
        else:
            high = np.flatnonzero(cldn4 >= med)
            low = np.flatnonzero(cldn4 < med)
    elif mode == "q4q1":
        ranks = pd.Series(cldn4).rank(method="average")
        try:
            qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
        except ValueError:
            return None
        if qs.nunique() < 4:
            return None
        high = np.flatnonzero(qs.to_numpy() == "Q4")
        low = np.flatnonzero(qs.to_numpy() == "Q1")
    else:
        raise ValueError(mode)
    if high.size < MIN_BIN or low.size < MIN_BIN:
        return None
    return high, low


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
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
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_rds_genes(path: Path, wanted: set[str]):
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read GSE205335 RDS", flush=True)
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
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
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"GSE205335 extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    return extracted, library_umi, barcodes, set(genes.tolist())


def _group_index(labels: np.ndarray):
    labels = np.asarray(labels)
    order = np.argsort(labels, kind="mergesort")
    sorted_lab = labels[order]
    cuts = np.flatnonzero(sorted_lab[1:] != sorted_lab[:-1]) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(order)]
    return [(sorted_lab[s], order[s:e]) for s, e in zip(starts, ends)]


def _member(members: dict[str, tarfile.TarInfo], sample: str, kind: str) -> str:
    for name in members:
        if f"_{sample}_{kind}" in Path(name).name:
            return name
    raise KeyError(f"no {kind} for {sample}")


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = _member(members, sample, "features")
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for line in handle:
            parts = line.decode().strip().split("\t")
            genes.append((parts[1] if len(parts) > 1 else parts[0]).upper())
    return genes


def _read_10x_nbarcodes(tf: tarfile.TarFile, members: dict, sample: str) -> int:
    name = _member(members, sample, "barcodes")
    n = 0
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for _ in handle:
            n += 1
    return n


def _stream_10x_wanted(tf, members, sample, genes, wanted) -> tuple[dict[str, np.ndarray], np.ndarray]:
    n_cells = _read_10x_nbarcodes(tf, members, sample)
    keep_idx = {i for i, g in enumerate(genes) if g in wanted}
    extracted = {genes[i]: np.zeros(n_cells, dtype=np.float32) for i in keep_idx}
    library = np.zeros(n_cells, dtype=np.float64)
    mtx_name = _member(members, sample, "matrix.mtx")
    with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as handle:
        for line in handle:
            if line.startswith(b"%"):
                continue
            break
        for line in handle:
            a, b, v = line.decode().split()
            gi = int(a) - 1
            ci = int(b) - 1
            val = float(v)
            library[ci] += val
            if gi in keep_idx:
                extracted[genes[gi]][ci] = val
    return extracted, library


def marker_mal_tnk(umi: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def col(g: str) -> np.ndarray:
        return umi[g] if g in umi else np.zeros(n, dtype=np.float32)

    epi = np.zeros(n, dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    mal = epi & (col("PTPRC") == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in TNK_MARKERS:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    return mal, tnk


def load_gse189357(tar_path: Path, wanted: set[str]) -> dict:
    print("== GSE189357 ==", flush=True)
    patients = []
    log_parts: dict[str, list[np.ndarray]] = {}
    pos_parts: dict[str, list[np.ndarray]] = {}
    mal_parts = []
    tnk_parts = []
    patient_parts = []
    gene_universe: set[str] | None = None
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            genes = _read_10x_features(tf, members, sample)
            if gene_universe is None:
                gene_universe = set(genes)
            else:
                gene_universe &= set(genes)
            extracted, library = _stream_10x_wanted(tf, members, sample, genes, wanted)
            n = int(library.size)
            lib = np.maximum(library, 1.0)
            for g in wanted:
                vec = extracted.get(g, np.zeros(n, dtype=np.float32))
                log_parts.setdefault(g, []).append(np.log1p(vec / lib * 1e4).astype(np.float32))
                pos_parts.setdefault(g, []).append((vec > 0).astype(np.float32))
            mal, tnk = marker_mal_tnk(extracted, n)
            mal_parts.append(mal)
            tnk_parts.append(tnk)
            patient_parts.append(np.array([sample] * n, dtype=object))
            patients.append(sample)
            print(
                f"  {sample}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())} "
                f"genes={len(extracted)}",
                flush=True,
            )
    log_cp = {g: np.concatenate(chunks) for g, chunks in log_parts.items() if len(chunks) == 9}
    pos = {g: np.concatenate(chunks) for g, chunks in pos_parts.items() if len(chunks) == 9}
    return {
        "cohort": "GSE189357",
        "log_cp": log_cp,
        "pos": pos,
        "patient": np.concatenate(patient_parts),
        "sample": np.concatenate(patient_parts),
        "mal": np.concatenate(mal_parts),
        "tnk": np.concatenate(tnk_parts),
        "genes": gene_universe or set(),
        "n_cells": int(sum(p.size for p in patient_parts)),
        "malig_def": "marker_malig",
    }


def load_gse205335(args, wanted: set[str]) -> dict:
    print("== GSE205335 ==", flush=True)
    identities = pd.read_csv(args.gse205335_identities, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("GSE205335 identity barcodes are not unique")
    metadata = parse_geo_soft(args.gse205335_soft)
    extracted, library_umi, barcodes, genes = load_rds_genes(args.gse205335_matrix, wanted)
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 missing from GSE205335 UMI")
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(
            f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("GSE205335 identity samples did not match GEO metadata")
    lib = np.maximum(cells["total_umi"].to_numpy(), 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    mal = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    tnk = cells["lineage.total"].eq("T/NK cells").to_numpy()
    return {
        "cohort": "GSE205335",
        "log_cp": log_cp,
        "pos": pos,
        "patient": cells["patient"].to_numpy(),
        "sample": cells["orig.ident"].to_numpy(),
        "mal": mal,
        "tnk": tnk,
        "genes": genes,
        "n_cells": int(len(barcodes)),
        "histology": cells["cancer_subtype"].to_numpy(),
        "recist": cells["recist"].to_numpy(),
        "malig_def": "author_malig",
    }


def score_paired(ds: dict, lr: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    elig_rows = []
    pair_rows = []
    for patient, idx in _group_index(ds["patient"]):
        mal = idx[ds["mal"][idx]]
        tnk = idx[ds["tnk"][idx]]
        if tnk.size < MIN_TNK or mal.size < MIN_BIN * 2:
            elig_rows.append(
                {
                    "cohort": ds["cohort"],
                    "patient": patient,
                    "split": mode,
                    "eligible": False,
                    "reason": "floor",
                    "n_malignant": int(mal.size),
                    "n_tnk": int(tnk.size),
                    "n_high": 0,
                    "n_low": 0,
                }
            )
            continue
        cldn4 = ds["log_cp"]["CLDN4"][mal]
        split = split_malignant(cldn4, mode)
        if split is None:
            elig_rows.append(
                {
                    "cohort": ds["cohort"],
                    "patient": patient,
                    "split": mode,
                    "eligible": False,
                    "reason": "thin_bins",
                    "n_malignant": int(mal.size),
                    "n_tnk": int(tnk.size),
                    "n_high": 0,
                    "n_low": 0,
                }
            )
            continue
        high_rel, low_rel = split
        high = mal[high_rel]
        low = mal[low_rel]
        elig_rows.append(
            {
                "cohort": ds["cohort"],
                "patient": patient,
                "split": mode,
                "eligible": True,
                "reason": "ok",
                "n_malignant": int(mal.size),
                "n_tnk": int(tnk.size),
                "n_high": int(high.size),
                "n_low": int(low.size),
                "mean_CLDN4_high": float(ds["log_cp"]["CLDN4"][high].mean()),
                "mean_CLDN4_low": float(ds["log_cp"]["CLDN4"][low].mean()),
            }
        )
        high_rows = score_lr(lr, ds["log_cp"], ds["pos"], high, tnk)
        low_rows = score_lr(lr, ds["log_cp"], ds["pos"], low, tnk)
        low_map = {r["interaction_name"]: r for r in low_rows}
        for hr in high_rows:
            lr_low = low_map[hr["interaction_name"]]
            pair_rows.append(
                {
                    "cohort": ds["cohort"],
                    "patient": patient,
                    "split": mode,
                    "direction": "outgoing",
                    "interaction_name": hr["interaction_name"],
                    "pathway_name": hr["pathway_name"],
                    "annotation": hr["annotation"],
                    "ligand": hr["ligand"],
                    "receptor": hr["receptor"],
                    "ligand_genes": hr["ligand_genes"],
                    "receptor_genes": hr["receptor_genes"],
                    "focus": pair_focus(hr["interaction_name"], hr["ligand_genes"], hr["pathway_name"]),
                    "ligand_class": ligand_class(hr["ligand_genes"]),
                    "prob_high": hr["prob"],
                    "prob_low": lr_low["prob"],
                    "detected_high": hr["detected"],
                    "detected_low": lr_low["detected"],
                    "detected_both": bool(hr["detected"] and lr_low["detected"]),
                    "n_high": int(high.size),
                    "n_low": int(low.size),
                    "n_tnk": int(tnk.size),
                }
            )
        print(
            f"  paired {ds['cohort']} {patient} {mode} high={high.size} low={low.size} tnk={tnk.size}",
            flush=True,
        )
    return pd.DataFrame(elig_rows), pd.DataFrame(pair_rows)


def wilcoxon_safe(high: np.ndarray, low: np.ndarray) -> tuple[float, float]:
    delta = high - low
    if delta.size < 3 or np.allclose(delta, 0):
        return float("nan"), 1.0
    try:
        res = stats.wilcoxon(high, low, alternative="two-sided", zero_method="wilcox")
        return float(res.statistic), float(res.pvalue)
    except ValueError:
        return float("nan"), 1.0


def contrast_paired(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if long.empty:
        return pd.DataFrame()
    for (split, name), block in long.groupby(["split", "interaction_name"], observed=True):
        both = block[block["detected_both"]]
        n = int(len(both))
        if n < MIN_DETECT_PAIRED:
            continue
        rec0 = block.iloc[0]
        w, p = wilcoxon_safe(both["prob_high"].to_numpy(), both["prob_low"].to_numpy())
        med_h = float(both["prob_high"].median())
        med_l = float(both["prob_low"].median())
        paired_delta = both["prob_high"] - both["prob_low"]
        rows.append(
            {
                "split": split,
                "direction": "outgoing",
                "interaction_name": name,
                "pathway_name": rec0["pathway_name"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "ligand_genes": rec0["ligand_genes"],
                "receptor_genes": rec0["receptor_genes"],
                "focus": rec0["focus"],
                "ligand_class": rec0["ligand_class"],
                "n_paired": n,
                "n_GSE189357": int((both["cohort"] == "GSE189357").sum()),
                "n_GSE205335": int((both["cohort"] == "GSE205335").sum()),
                "median_prob_high": med_h,
                "median_prob_low": med_l,
                "delta_median": med_h - med_l,
                "median_paired_delta": float(paired_delta.median()),
                "wilcoxon_w": w,
                "p": p,
                "thin": n < 8,
                "test": "wilcoxon_signed_rank_paired_patient",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return (
        out.assign(_abs=out["delta_median"].abs())
        .sort_values(["p", "_abs"], ascending=[True, False])
        .drop(columns="_abs")
    )


def score_between_unit(ds: dict, lr: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keep = set(units["unit_id"].astype(str))
    for unit, idx in _group_index(ds["patient"].astype(str)):
        if str(unit) not in keep:
            continue
        mal = idx[ds["mal"][idx]]
        tnk = idx[ds["tnk"][idx]]
        scored = score_lr(lr, ds["log_cp"], ds["pos"], mal, tnk)
        for rec in scored:
            rec["cohort"] = ds["cohort"]
            rec["unit_id"] = str(unit)
            rec["n_malignant"] = int(mal.size)
            rec["n_tnk"] = int(tnk.size)
            rec["direction"] = "outgoing"
            rec["focus"] = pair_focus(rec["interaction_name"], rec["ligand_genes"], rec["pathway_name"])
            rec["ligand_class"] = ligand_class(rec["ligand_genes"])
            rows.append(rec)
        print(f"  between {ds['cohort']} {unit} mal={mal.size} tnk={tnk.size}", flush=True)
    return pd.DataFrame(rows)


def assign_groups(values: pd.Series, mode: str) -> pd.Series:
    if mode == "median":
        med = float(values.median())
        return pd.Series(np.where(values >= med, "high", "low"), index=values.index)
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def contrast_between(long: pd.DataFrame, units: pd.DataFrame, mode: str) -> pd.DataFrame:
    frame = long.merge(
        units[["cohort", "unit_id", "cldn4_pct", "frac_tnk"]],
        on=["cohort", "unit_id"],
        how="left",
    )
    grouped = []
    for cohort, sub in frame.groupby("cohort", observed=True):
        labs = assign_groups(sub.drop_duplicates("unit_id").set_index("unit_id")["cldn4_pct"], mode)
        sub = sub.copy()
        sub["arm"] = sub["unit_id"].map(labs.astype(str))
        grouped.append(sub)
    frame = pd.concat(grouped, ignore_index=True)
    high_lab, low_lab = ("Q4", "Q1") if mode == "q4q1" else ("high", "low")
    rows = []
    for name, block in frame.groupby("interaction_name", observed=True):
        h = block[(block["arm"] == high_lab) & block["detected"]]
        l = block[(block["arm"] == low_lab) & block["detected"]]
        if len(h) < MIN_DETECT_ARM or len(l) < MIN_DETECT_ARM:
            continue
        u, p = stats.mannwhitneyu(h["prob"].to_numpy(), l["prob"].to_numpy(), alternative="two-sided")
        n_h, n_l = int(len(h)), int(len(l))
        rec0 = block.iloc[0]
        rows.append(
            {
                "split": f"between_{mode}",
                "direction": "outgoing",
                "interaction_name": name,
                "pathway_name": rec0["pathway_name"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "focus": rec0["focus"],
                "ligand_class": rec0["ligand_class"],
                "n_high": n_h,
                "n_low": n_l,
                "n_compared": n_h + n_l,
                "n_GSE189357": int(((h["cohort"] == "GSE189357").sum()) + (l["cohort"] == "GSE189357").sum()),
                "n_GSE205335": int(((h["cohort"] == "GSE205335").sum()) + (l["cohort"] == "GSE205335").sum()),
                "median_prob_high": float(h["prob"].median()),
                "median_prob_low": float(l["prob"].median()),
                "delta_median": float(h["prob"].median() - l["prob"].median()),
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": (2.0 * float(u)) / (n_h * n_l) - 1.0,
                "test": "mannwhitney_between_patient_given_combo",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["p", "delta_median"], ascending=[True, False])


def plot_ligand_bars(tbl: pd.DataFrame, path: Path, title: str, value="delta_median") -> None:
    fig, ax = plt.subplots(figsize=(8.6, max(2.2, 0.38 * max(len(tbl), 1) + 1.3)))
    if tbl.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No pairs at the stated gates", ha="center")
    else:
        show = tbl.head(18).copy()
        y = np.arange(len(show))
        colors = np.where(show[value] >= 0, "#b2182b", "#2166ac")
        ax.barh(y, show[value], color=colors, alpha=0.85)
        ax.set_yticks(y)
        labels = [
            f"{r.ligand}–{r.receptor}  n={int(getattr(r, 'n_paired', getattr(r, 'n_compared', 0)))}  p={fmt_p(r.p)}"
            for r in show.itertuples()
        ]
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="0.3", lw=0.8)
        ax.set_xlabel("ΔP (CLDN4-high − CLDN4-low)")
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_honest_n(elig: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    for ax, mode in zip(axes, ["median", "q4q1"]):
        sub = elig[elig["split"] == mode]
        ok = sub[sub["eligible"]]
        colors = np.where(ok["cohort"] == "GSE189357", "#4C72B0", "#C44E52")
        ax.scatter(ok["n_high"], ok["n_low"], c=colors, s=28)
        ax.axvline(MIN_BIN, color="0.5", ls="--", lw=0.7)
        ax.axhline(MIN_BIN, color="0.5", ls="--", lw=0.7)
        ax.set_xlabel("n malignant CLDN4-high")
        ax.set_title(f"{mode}: {int(ok.shape[0])} paired patients", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("n malignant CLDN4-low")
    fig.suptitle("Honest paired n (both bins ≥10; T/NK ≥20 already applied)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_n_bars(elig: pd.DataFrame, path: Path) -> None:
    sub = elig[(elig["split"] == "median")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), sharey=False)
    for ax, cohort in zip(axes, ["GSE189357", "GSE205335"]):
        block = sub[sub["cohort"] == cohort].sort_values("patient")
        x = np.arange(len(block))
        ax.bar(x - 0.18, block["n_malignant"], width=0.36, color="#b2182b", alpha=0.8, label="malignant")
        ax.bar(x + 0.18, block["n_tnk"], width=0.36, color="#6a8aaa", alpha=0.8, label="T/NK")
        ax.set_xticks(x)
        ax.set_xticklabels(block["patient"], rotation=90, fontsize=6)
        ax.set_title(f"{cohort} cells / patient", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_key_paired(long: pd.DataFrame, path: Path, mode: str) -> None:
    keys = [
        "CD274_PDCD1",
        "NECTIN2_TIGIT",
        "CXCL16_CXCR6",
        "CXCL9_CXCR3",
        "CXCL10_CXCR3",
        "CCL5_CCR5",
        "HLA-A_CD8A",
        "HLA-B_CD8A",
    ]
    sub = long[(long["split"] == mode) & long["interaction_name"].isin(keys) & long["detected_both"]]
    names = [k for k in keys if k in set(sub["interaction_name"])]
    fig, axes = plt.subplots(2, 4, figsize=(11.2, 5.6))
    axes = axes.ravel()
    for i, name in enumerate(names + [None] * (8 - len(names))):
        ax = axes[i]
        if name is None:
            ax.axis("off")
            continue
        b = sub[sub["interaction_name"] == name]
        for row in b.itertuples():
            color = "#4C72B0" if row.cohort == "GSE189357" else "#C44E52"
            ax.plot([0, 1], [row.prob_low, row.prob_high], color="0.65", lw=0.7)
            ax.scatter([0, 1], [row.prob_low, row.prob_high], c=color, s=14, zorder=3)
        ax.set_xticks([0, 1], ["low", "high"])
        ax.set_title(f"{name}\nn={len(b)}", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"Paired patient P (outgoing Mal→T/NK) — {mode}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_given_combo(units: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for cohort, color in (("GSE189357", "#4C72B0"), ("GSE205335", "#C44E52")):
        s = units[units["cohort"] == cohort]
        ax.scatter(s["cldn4_rank"], s["frac_tnk"], c=color, s=32, label=cohort, zorder=3)
    ax.set_xlabel("Within-cohort CLDN4 %pos rank (0–1)")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(
        "PR #459 given combo  n=31  ρ=−0.478  Q4 r=−0.750\n(not re-ranked)",
        fontsize=9,
    )
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def md_table(df: pd.DataFrame, n: int = 12) -> list[str]:
    if df.empty:
        return ["No pairs passed the detect / n gates.", ""]
    lines = [
        "| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |",
        "|---|---|---:|---|---:|---:|---:|---|",
    ]
    for rec in df.head(n).itertuples():
        star = " **" if rec.p < 0.05 else ""
        n189 = getattr(rec, "n_GSE189357", "")
        n205 = getattr(rec, "n_GSE205335", "")
        n_txt = f"{int(n189)}/{int(n205)}" if n189 != "" else f"{int(rec.n_high)}/{int(rec.n_low)}"
        n_use = int(getattr(rec, "n_paired", getattr(rec, "n_compared", 0)))
        lines.append(
            f"| {rec.focus} | {rec.ligand}–{rec.receptor}{star} | {n_use} | {n_txt} | "
            f"{rec.median_prob_high:.3f} | {rec.median_prob_low:.3f} | "
            f"{rec.delta_median:+.3f} | {fmt_p(rec.p)} |"
        )
    lines.append("")
    return lines


def locked_combo_units(root: Path) -> pd.DataFrame:
    a = pd.read_csv(root / "data" / "GSE189357_marker_units.tsv", sep="\t")
    a = a[a["eligible"]].copy()
    a["cohort"] = "GSE189357"
    a["unit_id"] = a["patient"].astype(str)
    a["cldn4_pct"] = a["mal_CLDN4_pct"] * 100.0
    a["frac_tnk"] = a["frac_tnk"]
    b = pd.read_csv(root / "data" / "GSE205335_patients.tsv", sep="\t")
    b["cohort"] = "GSE205335"
    b["unit_id"] = b["patient"].astype(str)
    b["cldn4_pct"] = b["mal_CLDN4_pct_pos"]
    keep = ["cohort", "unit_id", "cldn4_pct", "frac_tnk", "n_malignant", "n_tnk"]
    out = pd.concat([a[keep], b[keep]], ignore_index=True)
    out["cldn4_rank"] = out.groupby("cohort")["cldn4_pct"].rank(method="average", pct=True)
    return out


def force_key_table(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if long.empty:
        return pd.DataFrame()
    for (split, name), block in long.groupby(["split", "interaction_name"], observed=True):
        rec0 = block.iloc[0]
        if rec0["focus"] == "other" and name not in KEY_NAMES:
            continue
        both = block[block["detected_both"]]
        w, pval = (
            wilcoxon_safe(both["prob_high"].to_numpy(), both["prob_low"].to_numpy())
            if len(both)
            else (float("nan"), float("nan"))
        )
        rows.append(
            {
                "split": split,
                "direction": "outgoing",
                "interaction_name": name,
                "pathway_name": rec0["pathway_name"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "focus": rec0["focus"],
                "ligand_class": rec0["ligand_class"],
                "n_paired": int(len(both)),
                "n_GSE189357": int((both["cohort"] == "GSE189357").sum()) if len(both) else 0,
                "n_GSE205335": int((both["cohort"] == "GSE205335").sum()) if len(both) else 0,
                "median_prob_high": float(both["prob_high"].median()) if len(both) else np.nan,
                "median_prob_low": float(both["prob_low"].median()) if len(both) else np.nan,
                "delta_median": (
                    float(both["prob_high"].median() - both["prob_low"].median()) if len(both) else np.nan
                ),
                "median_paired_delta": float((both["prob_high"] - both["prob_low"]).median()) if len(both) else np.nan,
                "p": pval,
                "thin": len(both) < 8,
                "test": "wilcoxon_signed_rank_paired_patient",
                "note": "" if len(both) >= MIN_DETECT_PAIRED else f"n_paired={len(both)} below gate {MIN_DETECT_PAIRED}",
            }
        )
    key = pd.DataFrame(rows)
    if key.empty:
        return key
    return key.sort_values(["split", "focus", "p"])


def write_finding(
    out: Path,
    paired: pd.DataFrame,
    between: pd.DataFrame,
    elig: pd.DataFrame,
    key: pd.DataFrame,
) -> None:
    med_elig = elig[(elig["split"] == "median") & elig["eligible"]]
    q_elig = elig[(elig["split"] == "q4q1") & elig["eligible"]]
    med_tab = paired[paired["split"] == "median"] if not paired.empty else paired
    q_tab = paired[paired["split"] == "q4q1"] if not paired.empty else paired
    key_med = key[key["split"] == "median"] if not key.empty else key
    key_q = key[key["split"] == "q4q1"] if not key.empty else key
    bet_med = between[between["split"] == "between_median"] if not between.empty else between
    bet_q = between[between["split"] == "between_q4q1"] if not between.empty else between

    def n_sig(df):
        return int((df["p"] < 0.05).sum()) if not df.empty else 0

    dropped_med = int(((elig.split == "median") & ~elig.eligible).sum())
    dropped_q = int(((elig.split == "q4q1") & ~elig.eligible).sum())
    lines = [
        "# FINDING — pair GSE189357 + GSE205335: CLDN4-only high-end Mal→T/NK",
        "",
        "ADDITIVE **CLDN4 only** on the pairwise combo that already **differs**",
        "(PR #459: malignant CLDN4 %pos vs same-patient T/NK **n=31 ρ=−0.478** p=0.009 I²=0%;",
        "Q4 vs Q1 **r=−0.750** p=0.010, 8/8). **No TACSTD2∩CLDN4 dual-high.**",
        "The given Spearman / Q4 row is not re-ranked. p-values are descriptive.",
        "CellChat R and LIANA were not run.",
        "",
        "## Given combo (PR #459; not re-audited)",
        "",
        "| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---:|---|---|",
        "| GSE189357+GSE205335 | %pos | 31 | −0.478 (0.009, 0%) | −0.750 (0.010; 8/8) |",
        "",
        "Members: GSE189357 marker-malignant n=9 + GSE205335 author-malignant n=22.",
        "Malignant definitions differ and are stated on every row.",
        "",
        "## Method (one paragraph)",
        "",
        "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs: 10% truncated",
        "mean of `log1p(CP10k)`, complexes = geometric mean of subunits, detected if",
        f"each complex has expressing-cell fraction ≥ {EXPR_PROP:.2f},",
        f"P = (L·R)/({KH}+L·R). Outgoing = malignant → **same-patient** T/NK.",
        "Malignant cells are split **within each patient** by CLDN4 (median and Q4 vs Q1).",
        f"A patient enters only with both bins ≥{MIN_BIN} cells and T/NK ≥{MIN_TNK}.",
        "The test is Wilcoxon signed-rank on per-patient P (pair detected on both arms,",
        f"n≥{MIN_DETECT_PAIRED}). Unit = patient. GSE189357 malignant = marker-malignant",
        "(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0; GSE205335 malignant = author",
        "`Malignant cells`.",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        "| Given combo (PR #459 %pos) | **31** | 9 + 22; do not re-audit |",
        "| Given Q4 vs Q1 tails | **8 vs 8** | within-cohort, stacked |",
        (
            f"| High-end median paired | **{int(len(med_elig))}** | "
            f"GSE189357 {int((med_elig.cohort=='GSE189357').sum())} + "
            f"GSE205335 {int((med_elig.cohort=='GSE205335').sum())}; "
            f"dropped {dropped_med} |"
        ),
        (
            f"| High-end Q4 vs Q1 paired | **{int(len(q_elig))}** | "
            f"GSE189357 {int((q_elig.cohort=='GSE189357').sum())} + "
            f"GSE205335 {int((q_elig.cohort=='GSE205335').sum())}; "
            f"dropped {dropped_q} |"
        ),
        "",
        "Cells are not n. The given n=31 is the combo Spearman n, not the communication n.",
        "Patients that lack a high bin, a low bin, or T/NK are out of the paired test.",
        "",
        "## Primary ligand table — within-patient paired (outgoing)",
        "",
        f"Detect-gated pairs (median): {int(len(med_tab))} "
        f"({n_sig(med_tab)} with p<0.05). "
        f"Q4 vs Q1: {int(len(q_tab))} ({n_sig(q_tab)} with p<0.05).",
        "",
        "### Median split",
        "",
    ]
    lines += md_table(med_tab)
    lines += ["### Q4 vs Q1 split", ""]
    lines += md_table(q_tab)
    lines += [
        "## Focused pairs (always shown if subunits exist)",
        "",
        "MHC-I (HLA-A/B/C–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1, NECTIN2–TIGIT.",
        "These rows are the same test; they are not a second discovery pass.",
        "",
        "### Median",
        "",
    ]
    lines += md_table(key_med, n=20)
    lines += ["### Q4 vs Q1", ""]
    lines += md_table(key_q, n=20)
    lines += [
        "## Companion — between-patient on the given 31-patient combo",
        "",
        "All-malignant → same-patient T/NK. Quartiles / median are **within cohort**,",
        "then stacked. This is the slice whose T/NK association is already given",
        "(not re-ranked). Test = Mann–Whitney on per-patient Mal→T/NK P",
        f"(detected ≥{MIN_DETECT_ARM} per arm).",
        "",
        "### Between-patient median",
        "",
    ]
    lines += md_table(bet_med)
    lines += ["### Between-patient Q4 vs Q1", ""]
    lines += md_table(bet_q)
    lines += [
        "## What this does not say",
        "",
        "- It does not build a TACSTD2∩CLDN4 dual-high score.",
        "- It does not re-audit the given PR #459 Spearman / Q4 row.",
        "- It does not treat this pair as GSE131907+GSE205335 (that is a different folder).",
        "- Cell-pooled stacked means are not the test.",
        "- GSE205335 Q4 remains SCLC-heavy on the given between-patient tails.",
        "- No FASTQ. No CellChat R. No LIANA.",
        "",
        "## Files",
        "",
        "- `results/ligand_table.tsv` — patient-level paired LR (n / Δ / p)",
        "- `results/ligand_table_key.tsv` — MHC-I / recruit / CD274–PDCD1 / NECTIN2–TIGIT",
        "- `results/ligand_table_between.tsv` — given-combo between-patient companion",
        "- `results/eligibility.tsv` — honest paired n",
        "- `figures/` — extra figures",
        "- `METHODS.md` — labels, Hill probability, floors",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def wanted_genes(db: Path) -> set[str]:
    wanted = set(EXTRA)
    inter = pd.read_csv(db / "interaction_cellchatdb_v2_protein.csv")
    for rec in inter.itertuples(index=False):
        wanted.update(parse_symbols(getattr(rec, "ligand_symbol", rec.ligand)))
        wanted.update(parse_symbols(getattr(rec, "receptor_symbol", rec.receptor)))
        wanted.update(parse_symbols(rec.ligand))
        wanted.update(parse_symbols(rec.receptor))
    return wanted


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/geo_pair_189357_205335"))
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--db", type=Path, default=ROOT / "db")
    args = p.parse_args()
    args.gse189357_tar = args.data / "gse189357" / "GSE189357_RAW.tar"
    args.gse205335_matrix = args.data / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    args.gse205335_identities = args.data / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    args.gse205335_soft = args.data / "gse205335" / "GSE205335_family.soft.gz"
    for path in (
        args.gse189357_tar,
        args.gse205335_matrix,
        args.gse205335_identities,
        args.gse205335_soft,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run scripts/download.py")

    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    wanted = wanted_genes(args.db)
    print(f"wanted genes={len(wanted)}", flush=True)
    ds189 = load_gse189357(args.gse189357_tar, wanted)
    ds205 = load_gse205335(args, wanted)
    both_genes = set(ds189["log_cp"]) & set(ds205["log_cp"])
    lr = load_lr(args.db, both_genes)
    print(f"LR pairs with all subunits in both matrices: {len(lr)}", flush=True)
    if "CLDN4" not in ds189["log_cp"]:
        raise SystemExit("CLDN4 missing from GSE189357")

    elig_parts, long_parts = [], []
    for ds in (ds189, ds205):
        for mode in ("median", "q4q1"):
            elig, long = score_paired(ds, lr, mode)
            elig_parts.append(elig)
            long_parts.append(long)
    elig = pd.concat(elig_parts, ignore_index=True)
    long = pd.concat([x for x in long_parts if not x.empty], ignore_index=True)
    elig.to_csv(args.out / "results" / "eligibility.tsv", sep="\t", index=False)
    long.to_csv(args.out / "results" / "per_patient_lr_paired.tsv.gz", sep="\t", index=False)

    paired = contrast_paired(long)
    paired.to_csv(args.out / "results" / "ligand_table.tsv", sep="\t", index=False)
    key = force_key_table(long)
    key.to_csv(args.out / "results" / "ligand_table_key.tsv", sep="\t", index=False)

    units = locked_combo_units(ROOT)
    units.to_csv(args.out / "results" / "given_combo_units.tsv", sep="\t", index=False)
    between_long = pd.concat(
        [
            score_between_unit(ds189, lr, units[units.cohort == "GSE189357"]),
            score_between_unit(ds205, lr, units[units.cohort == "GSE205335"]),
        ],
        ignore_index=True,
    )
    between_long.to_csv(args.out / "results" / "per_unit_lr_between.tsv.gz", sep="\t", index=False)
    between = pd.concat(
        [contrast_between(between_long, units, "median"), contrast_between(between_long, units, "q4q1")],
        ignore_index=True,
    )
    between.to_csv(args.out / "results" / "ligand_table_between.tsv", sep="\t", index=False)

    plot_ligand_bars(
        paired[paired["split"] == "median"] if not paired.empty else paired,
        args.out / "figures" / "fig_extra_ligand_table.png",
        "Pair 189357+205335 outgoing ΔP  (within-patient median; paired)",
    )
    plot_ligand_bars(
        paired[paired["split"] == "q4q1"] if not paired.empty else paired,
        args.out / "figures" / "fig_q4q1_ligand_table.png",
        "Pair 189357+205335 outgoing ΔP  (within-patient Q4 vs Q1; paired)",
    )
    plot_ligand_bars(
        key[key["split"] == "median"] if not key.empty else key,
        args.out / "figures" / "fig_key_pairs_median.png",
        "Focused pairs  MHC-I / T-recruit / CD274–PDCD1 / NECTIN2–TIGIT  (median)",
    )
    plot_honest_n(elig, args.out / "figures" / "fig_honest_n.png")
    plot_n_bars(elig, args.out / "figures" / "fig_n_per_patient.png")
    if not long.empty:
        plot_key_paired(long, args.out / "figures" / "fig_key_pairs_paired_median.png", "median")
        plot_key_paired(long, args.out / "figures" / "fig_key_pairs_paired_q4q1.png", "q4q1")
    plot_given_combo(units, args.out / "figures" / "fig_given_combo_cldn4_tnk.png")
    if not between.empty:
        plot_ligand_bars(
            between[between["split"] == "between_q4q1"],
            args.out / "figures" / "fig_between_q4q1.png",
            "Companion: between-patient Q4 vs Q1 on the given n=31 combo",
        )

    given = json.loads((ROOT / "given_combo.json").read_text())
    summary = {
        "datasets": ["GSE189357", "GSE205335"],
        "marker": "CLDN4",
        "additive": True,
        "dual_high": False,
        "algorithm": "Jin 2021 Hill Kh=0.5, 10% trim mean, expr_prop>=0.10, CellChatDB v2 protein",
        "cellchat_r": False,
        "liana": False,
        "unit": "patient",
        "floors": {"min_bin": MIN_BIN, "min_tnk": MIN_TNK, "min_detect_paired": MIN_DETECT_PAIRED},
        "n_lr_both_matrices": int(len(lr)),
        "n_paired_median": int(((elig.split == "median") & elig.eligible).sum()),
        "n_paired_q4q1": int(((elig.split == "q4q1") & elig.eligible).sum()),
        "n_ligand_table": int(len(paired)),
        "n_ligand_p_lt_05_median": int(((paired.split == "median") & (paired.p < 0.05)).sum())
        if not paired.empty
        else 0,
        "n_ligand_p_lt_05_q4q1": int(((paired.split == "q4q1") & (paired.p < 0.05)).sum())
        if not paired.empty
        else 0,
        "pr459_given": given,
        "not_used": "FASTQ; TACSTD2 groups; dual-high; CellChat R; LIANA; GSE131907",
        "matrix_sha256": {
            "GSE189357": sha256(args.gse189357_tar),
            "GSE205335": sha256(args.gse205335_matrix),
        },
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(args.out, paired, between, elig, key)
    slim = {k: summary[k] for k in summary if k != "matrix_sha256"}
    print(json.dumps(slim, indent=2), flush=True)
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
