#!/usr/bin/env python3
"""Merged GSE131907 + GSE205335 CellChat-style CLDN4-high vs low → T/NK.

ADDITIVE. CLDN4 only. Patient is the unit.
CellChat R / LIANA are not run. Probability is Jin et al. 2021
(10% truncated mean, Hill K_h=0.5, CellChatDB v2 protein pairs).

Primary test: within-patient split of malignant cells (median and Q4 vs Q1),
outgoing to same-patient T/NK, Wilcoxon signed-rank on per-patient P.
Honest paired n requires both CLDN4 bins and a T/NK floor.

Companion: between-unit Mal→T/NK on the locked PR #320 slice
(GSE131907 author-Malignant-cells samples n=21 + GSE205335 patients n=22).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
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
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "CD8B", "NKG7", "IFNG", "NECTIN2", "TIGIT"]
# Kim et al. GSE131907 still uses pre-HGNC nectin symbols.
GENE_ALIAS = {"PVRL1": "NECTIN1", "PVRL2": "NECTIN2", "PVRL3": "NECTIN3", "PVRL4": "NECTIN4"}

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
    """Return (high_idx_into_cldn4, low_idx_into_cldn4) or None if floors fail."""
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


def stream_matrix(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = GENE_ALIAS.get(gene.split(".")[0], gene.split(".")[0])
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  GSE131907 stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"GSE131907 stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


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


def load_gse131907(args, wanted: set[str]) -> dict:
    print("== GSE131907 ==", flush=True)
    cell_ids, found, n_umi, n_streamed = stream_matrix(args.gse131907_matrix, wanted)
    if "CLDN4" not in found:
        raise SystemExit("CLDN4 missing from GSE131907 UMI")
    ann = pd.read_csv(args.gse131907_ann, sep="\t", dtype=str)
    per = ann.set_index("Index").reindex(cell_ids).reset_index()
    if per["Sample"].isna().any():
        raise SystemExit("GSE131907 matrix/annotation mismatch")
    series = parse_series_matrix(args.gse131907_series)
    sample_meta = series.rename(columns={"title": "Sample"})
    keep_cols = [
        c
        for c in [
            "Sample",
            "geo_accession",
            "patient_id",
            "tumor_stage",
            "tissue_origin_abbrevation",
            "source_name_ch1",
        ]
        if c in sample_meta.columns
    ]
    sample_meta = sample_meta[keep_cols].drop_duplicates("Sample")
    pmap = sample_meta.set_index("Sample")["patient_id"]
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    pos = {g: (found[g] > 0).astype(np.float32) for g in found}
    epi = per["Cell_type"].eq("Epithelial cells")
    mal_ts = epi & per["Cell_subtype"].isin(MALIG_SUBTYPES)
    mal_author = epi & per["Cell_subtype"].eq("Malignant cells")
    tnk = per["Cell_type"].isin(["T lymphocytes", "NK cells"])
    tumor = per["Sample_Origin"].isin(TUMOR_ORIGINS)
    patient = np.array([pmap[s] if s in pmap.index else s for s in per["Sample"]], dtype=object)
    return {
        "cohort": "GSE131907",
        "log_cp": log_cp,
        "pos": pos,
        "umi": found,
        "patient": patient,
        "sample": per["Sample"].to_numpy(),
        "origin": per["Sample_Origin"].to_numpy(),
        "mal_ts": (mal_ts & tumor).to_numpy(),
        "mal_author": (mal_author & tumor).to_numpy(),
        "tnk": (tnk & tumor).to_numpy(),
        "genes": set(found) | set(wanted),
        "n_streamed": n_streamed,
        "n_cells": int(len(cell_ids)),
        "sample_meta": sample_meta,
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
        raise ValueError(f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
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
        "umi": extracted,
        "patient": cells["patient"].to_numpy(),
        "sample": cells["orig.ident"].to_numpy(),
        "origin": cells["tissue"].fillna("").to_numpy(),
        "mal_ts": mal,
        "mal_author": mal,
        "tnk": tnk,
        "genes": genes,
        "n_cells": int(len(barcodes)),
        "histology": cells["cancer_subtype"].to_numpy(),
        "recist": cells["recist"].to_numpy(),
    }


def patient_inventory(ds: dict, mal_key: str) -> pd.DataFrame:
    rows = []
    for patient, idx in _group_index(ds["patient"]):
        mal = idx[ds[mal_key][idx]]
        tnk = idx[ds["tnk"][idx]]
        samples = sorted(set(ds["sample"][idx[ds[mal_key][idx] | ds["tnk"][idx]]].tolist()))
        origins = sorted(set(ds["origin"][mal].tolist())) if mal.size else []
        cldn = ds["log_cp"]["CLDN4"][mal] if mal.size else np.array([])
        umi = ds["umi"]["CLDN4"][mal] if mal.size else np.array([])
        rows.append(
            {
                "cohort": ds["cohort"],
                "patient": patient,
                "n_malignant": int(mal.size),
                "n_tnk": int(tnk.size),
                "n_samples": int(len(samples)),
                "samples": ",".join(samples),
                "origins": ",".join(str(x) for x in origins),
                "mal_CLDN4_mean": float(cldn.mean()) if mal.size else np.nan,
                "mal_CLDN4_pct_pos": float(100 * (umi > 0).mean()) if mal.size else np.nan,
                "frac_tnk": float(tnk.size / max(mal.size + tnk.size, 1)),
                "eligible_floors": bool(mal.size >= MIN_BIN * 2 and tnk.size >= MIN_TNK),
            }
        )
    return pd.DataFrame(rows)


def _group_index(labels: np.ndarray) -> list[tuple[object, np.ndarray]]:
    order = np.argsort(labels, kind="mergesort")
    sorted_lab = labels[order]
    cuts = np.flatnonzero(sorted_lab[1:] != sorted_lab[:-1]) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(order)]
    return [(sorted_lab[s], order[s:e]) for s, e in zip(starts, ends)]


def score_paired(ds: dict, lr: pd.DataFrame, mal_key: str, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    elig_rows = []
    pair_rows = []
    for patient, idx in _group_index(ds["patient"]):
        mal = idx[ds[mal_key][idx]]
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
                "n_GSE131907": int((both["cohort"] == "GSE131907").sum()),
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
    return out.assign(_abs=out["delta_median"].abs()).sort_values(
        ["p", "_abs"], ascending=[True, False]
    ).drop(columns="_abs")


def score_between_unit(ds: dict, lr: pd.DataFrame, units: pd.DataFrame, mal_key: str) -> pd.DataFrame:
    """All-malignant → T/NK per locked unit (sample or patient)."""
    rows = []
    keep = set(units["unit_id"].astype(str))
    labels = ds["sample"] if ds["cohort"] == "GSE131907" else ds["patient"]
    for unit, idx in _group_index(labels.astype(str)):
        if str(unit) not in keep:
            continue
        mal = idx[ds[mal_key][idx]]
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
    frame = long.merge(units[["cohort", "unit_id", "cldn4_pct", "frac_tnk"]], on=["cohort", "unit_id"], how="left")
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
                "median_prob_high": float(h["prob"].median()),
                "median_prob_low": float(l["prob"].median()),
                "delta_median": float(h["prob"].median() - l["prob"].median()),
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": (2.0 * float(u)) / (n_h * n_l) - 1.0,
                "test": "mannwhitney_between_unit_PR320_slice",
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
        ax.scatter(ok["n_high"], ok["n_low"], c=np.where(ok["cohort"] == "GSE131907", "#4C72B0", "#C44E52"), s=28)
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


def plot_key_paired(long: pd.DataFrame, path: Path, mode: str) -> None:
    keys = ["CD274_PDCD1", "NECTIN2_TIGIT", "CXCL16_CXCR6", "CXCL9_CXCR3", "CXCL10_CXCR3", "CCL5_CCR5", "HLA-A_CD8A", "HLA-B_CD8A"]
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
        for j, row in enumerate(b.itertuples()):
            ax.plot([0, 1], [row.prob_low, row.prob_high], color="0.65", lw=0.7)
            ax.scatter([0, 1], [row.prob_low, row.prob_high], c="#4C72B0" if row.cohort == "GSE131907" else "#C44E52", s=14, zorder=3)
        ax.set_xticks([0, 1], ["low", "high"])
        ax.set_title(f"{name}\nn={len(b)}", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"Paired patient P (outgoing Mal→T/NK) — {mode}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_given_slice(units: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for cohort, color in (("GSE131907", "#4C72B0"), ("GSE205335", "#C44E52")):
        s = units[units["cohort"] == cohort]
        ax.scatter(s["cldn4_pct"], s["frac_tnk"], c=color, s=32, label=cohort, zorder=3)
    rho, p = stats.spearmanr(units["cldn4_pct"], units["frac_tnk"])
    ax.set_xlabel("Malignant CLDN4 %pos (locked extract)")
    ax.set_ylabel("Same-unit T/NK fraction")
    ax.set_title(f"PR #320 slice n={len(units)}  ρ={rho:+.3f} p={fmt_p(p)}\n(given; not re-ranked)", fontsize=9)
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
        "| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |",
        "|---|---|---:|---|---:|---:|---:|---|",
    ]
    for rec in df.head(n).itertuples():
        star = " **" if rec.p < 0.05 else ""
        n131 = getattr(rec, "n_GSE131907", "")
        n205 = getattr(rec, "n_GSE205335", "")
        n_txt = f"{int(n131)}/{int(n205)}" if n131 != "" else f"{int(rec.n_high)}/{int(rec.n_low)}"
        n_use = int(getattr(rec, "n_paired", getattr(rec, "n_compared", 0)))
        lines.append(
            f"| {rec.focus} | {rec.ligand}–{rec.receptor}{star} | {n_use} | {n_txt} | "
            f"{rec.median_prob_high:.3f} | {rec.median_prob_low:.3f} | "
            f"{rec.delta_median:+.3f} | {fmt_p(rec.p)} |"
        )
    lines.append("")
    return lines


def write_finding(out: Path, summary: dict, paired: pd.DataFrame, between: pd.DataFrame, elig: pd.DataFrame, key: pd.DataFrame) -> None:
    def _considered(frame: pd.DataFrame) -> pd.DataFrame:
        return frame[(frame["n_malignant"] > 0) & (frame["n_tnk"] > 0)]

    med_all = _considered(elig[elig["split"] == "median"])
    q_all = _considered(elig[elig["split"] == "q4q1"])
    med_elig = med_all[med_all["eligible"]]
    q_elig = q_all[q_all["eligible"]]
    med_tab = paired[paired["split"] == "median"] if not paired.empty else paired
    q_tab = paired[paired["split"] == "q4q1"] if not paired.empty else paired
    key_med = key[key["split"] == "median"] if not key.empty else key
    key_q = key[key["split"] == "q4q1"] if not key.empty else key
    bet_med = between[between["split"] == "between_median"] if not between.empty else between
    bet_q = between[between["split"] == "between_q4q1"] if not between.empty else between

    def n_sig(df):
        return int((df["p"] < 0.05).sum()) if not df.empty else 0

    lines = [
        "# FINDING — merged GSE131907 + GSE205335 CellChat-style CLDN4-high → T/NK",
        "",
        "ADDITIVE **CLDN4 only** on the merged public UMI slice that already differs",
        "(PR #320: malignant CLDN4 Q4 vs Q1 vs same-patient T/NK n=23 r=−0.705 p=0.0003;",
        "continuous n=43 ρ=−0.479). **No TACSTD2∩CLDN4 dual-high.** GSE207422 was not run.",
        "p-values are descriptive. CellChat R was not run.",
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
        f"n≥{MIN_DETECT_PAIRED}). Unit = GEO patient (GSE131907 tumor origins pooled;",
        "PE / nLung / nLN dropped). GSE131907 malignant = author `Malignant cells` + tS1/tS2/tS3;",
        "GSE205335 malignant = author `Malignant cells`.",
        "",
        "## Honest paired n",
        "",
        f"| split | considered (mal>0 and T/NK>0) | eligible | GSE131907 | GSE205335 | out (floor / thin bins) |",
        f"|---|---:|---:|---:|---:|---:|",
        (
            f"| median | {int(len(med_all))} | {int(len(med_elig))} | "
            f"{int((med_elig.cohort=='GSE131907').sum())} | "
            f"{int((med_elig.cohort=='GSE205335').sum())} | "
            f"{int((~med_all.eligible).sum())} |"
        ),
        (
            f"| Q4 vs Q1 | {int(len(q_all))} | {int(len(q_elig))} | "
            f"{int((q_elig.cohort=='GSE131907').sum())} | "
            f"{int((q_elig.cohort=='GSE205335').sum())} | "
            f"{int((~q_all.eligible).sum())} |"
        ),
        "",
        "Cells are not n. nLung / nLN / PE-only patients (0 tumor malignant) are not",
        "in the considered column. Q4 vs Q1 loses patients when `qcut` cannot form",
        "four ranks (CLDN4 ties, often a large zero mass).",
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
    def key_block(df: pd.DataFrame, focuses: list[str]) -> pd.DataFrame:
        if df.empty:
            return df
        show = df[df["focus"].isin(focuses) & (df["n_paired"] > 0)].copy()
        order = {k: i for i, k in enumerate(focuses)}
        show["_o"] = show["focus"].map(order)
        return show.sort_values(["_o", "p"]).drop(columns="_o")

    focus_order = ["MHC-I", "T-recruit", "CD274–PDCD1", "NECTIN2–TIGIT"]
    lines += [
        "## Focused pairs (same test; not a second discovery pass)",
        "",
        "MHC-I (HLA-A/B/C/E–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1,",
        "NECTIN2–TIGIT. Rows with n_paired=0 (never detected on both arms) are omitted",
        "here and kept in `ligand_table_key.tsv`. CXCL9–CXCR3 was not detected.",
        "",
        "### Median",
        "",
    ]
    for split_name, kdf in (("Median", key_med), ("Q4 vs Q1", key_q)):
        lines += [f"### {split_name}", ""]
        for foc in focus_order:
            sub = key_block(kdf, [foc])
            if sub.empty:
                lines += [f"**{foc}.** Not detected on both arms in ≥1 patient.", ""]
                continue
            lines += [f"**{foc}**", ""]
            lines += md_table(sub, n=12)
    lines += [
        "## Companion — between-unit on the locked PR #320 slice",
        "",
        "Same 21 GSE131907 author-`Malignant cells` samples (n_mal≥20) + 22 GSE205335",
        "patients. Quartiles / median are **within cohort**, then stacked. This is the",
        "slice whose T/NK association is already given (not re-ranked). Test =",
        "Mann–Whitney on per-unit Mal→T/NK P (detected ≥3 per arm).",
        "",
        "### Between-unit median",
        "",
    ]
    lines += md_table(bet_med)
    lines += ["### Between-unit Q4 vs Q1", ""]
    lines += md_table(bet_q)
    lines += [
        "## What this does not say",
        "",
        "- It does not build a TACSTD2∩CLDN4 dual-high score.",
        "- It does not call the thesis a failure. The given T/NK association on this",
        "  merged slice stands (PR #320).",
        "- Cell-pooled stacked means are not the test.",
        "- GSE207422 patient-level CLDN4 vs T/NK is flat and was not run.",
        "- No FASTQ. The 2.86 GB GSE131907 log2TPM text matrix was skipped (UMI exists).",
        "",
        "## Files",
        "",
        "- `results/ligand_table.tsv` — patient-level paired LR (n / Δ / p)",
        "- `results/ligand_table_key.tsv` — MHC-I / recruit / CD274–PDCD1 / NECTIN2–TIGIT",
        "- `results/ligand_table_between.tsv` — PR #320-slice between-unit companion",
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


def locked_pr320_units(root: Path) -> pd.DataFrame:
    s131 = pd.read_csv(root / "data" / "GSE131907_samples.tsv", sep="\t")
    s131 = s131[s131["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    s131 = s131[s131["n_malignant"] >= 20].copy()
    s131["cohort"] = "GSE131907"
    s131["unit_id"] = s131["sample"].astype(str)
    s131["cldn4_pct"] = s131["mal_CLDN4_pct"]
    s131["frac_tnk"] = s131["frac_tnk"]
    s205 = pd.read_csv(root / "data" / "GSE205335_patients.tsv", sep="\t")
    s205["cohort"] = "GSE205335"
    s205["unit_id"] = s205["patient"].astype(str)
    s205["cldn4_pct"] = s205["mal_CLDN4_pct_pos"]
    keep = ["cohort", "unit_id", "cldn4_pct", "frac_tnk", "n_malignant", "n_tnk"]
    return pd.concat([s131[keep], s205[keep]], ignore_index=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/geo_merge_131907_205335"))
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--db", type=Path, default=ROOT / "db")
    args = p.parse_args()
    args.gse131907_matrix = args.data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    args.gse131907_ann = args.data / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    args.gse131907_series = args.data / "gse131907" / "GSE131907_series_matrix.txt.gz"
    args.gse205335_matrix = args.data / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    args.gse205335_identities = args.data / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    args.gse205335_soft = args.data / "gse205335" / "GSE205335_family.soft.gz"
    for path in (
        args.gse131907_matrix,
        args.gse131907_ann,
        args.gse131907_series,
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
    ds131 = load_gse131907(args, wanted)
    ds205 = load_gse205335(args, wanted)
    genes = (ds131["genes"] & set(ds131["log_cp"])) | (ds205["genes"] & set(ds205["log_cp"]))
    # A pair is kept if every subunit is in BOTH matrices (merged scoring).
    both_genes = set(ds131["log_cp"]) & set(ds205["log_cp"])
    lr = load_lr(args.db, both_genes)
    print(f"LR pairs with all subunits in both matrices: {len(lr)}", flush=True)

    inv131 = patient_inventory(ds131, "mal_ts")
    inv205 = patient_inventory(ds205, "mal_author")
    inv = pd.concat([inv131, inv205], ignore_index=True)
    inv.to_csv(args.out / "results" / "patient_inventory.tsv", sep="\t", index=False)

    elig_parts, long_parts = [], []
    for ds, mal_key in ((ds131, "mal_ts"), (ds205, "mal_author")):
        for mode in ("median", "q4q1"):
            elig, long = score_paired(ds, lr, mal_key, mode)
            elig_parts.append(elig)
            long_parts.append(long)
    elig = pd.concat(elig_parts, ignore_index=True)
    long = pd.concat([x for x in long_parts if not x.empty], ignore_index=True)
    elig.to_csv(args.out / "results" / "eligibility.tsv", sep="\t", index=False)
    long.to_csv(args.out / "results" / "per_patient_lr_paired.tsv.gz", sep="\t", index=False)

    paired = contrast_paired(long)
    paired.to_csv(args.out / "results" / "ligand_table.tsv", sep="\t", index=False)
    key = paired[paired["focus"] != "other"].copy() if not paired.empty else paired
    # Always attach key names even if they missed the n gate, from long.
    key_force_rows = []
    if not long.empty:
        for (split, name), block in long.groupby(["split", "interaction_name"], observed=True):
            if name not in KEY_NAMES and pair_focus(name, block.iloc[0]["ligand_genes"], block.iloc[0]["pathway_name"]) == "other":
                continue
            both = block[block["detected_both"]]
            rec0 = block.iloc[0]
            if rec0["focus"] == "other" and name not in KEY_NAMES:
                continue
            w, pval = wilcoxon_safe(both["prob_high"].to_numpy(), both["prob_low"].to_numpy()) if len(both) else (float("nan"), float("nan"))
            key_force_rows.append(
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
                    "n_GSE131907": int((both["cohort"] == "GSE131907").sum()) if len(both) else 0,
                    "n_GSE205335": int((both["cohort"] == "GSE205335").sum()) if len(both) else 0,
                    "median_prob_high": float(both["prob_high"].median()) if len(both) else np.nan,
                    "median_prob_low": float(both["prob_low"].median()) if len(both) else np.nan,
                    "delta_median": float(both["prob_high"].median() - both["prob_low"].median()) if len(both) else np.nan,
                    "median_paired_delta": float((both["prob_high"] - both["prob_low"]).median()) if len(both) else np.nan,
                    "p": pval,
                    "thin": len(both) < 8,
                    "test": "wilcoxon_signed_rank_paired_patient",
                    "note": "" if len(both) >= MIN_DETECT_PAIRED else f"n_paired={len(both)} below gate {MIN_DETECT_PAIRED}",
                }
            )
    key = pd.DataFrame(key_force_rows)
    if not key.empty:
        key = key.sort_values(["split", "focus", "p"])
    key.to_csv(args.out / "results" / "ligand_table_key.tsv", sep="\t", index=False)

    units = locked_pr320_units(ROOT)
    units.to_csv(args.out / "results" / "pr320_locked_units.tsv", sep="\t", index=False)
    between_long = pd.concat(
        [
            score_between_unit(ds131, lr, units[units.cohort == "GSE131907"], "mal_author"),
            score_between_unit(ds205, lr, units[units.cohort == "GSE205335"], "mal_author"),
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
        "Merged CellChat-style outgoing ΔP  (within-patient median; paired)",
    )
    plot_ligand_bars(
        paired[paired["split"] == "q4q1"] if not paired.empty else paired,
        args.out / "figures" / "fig_q4q1_ligand_table.png",
        "Merged CellChat-style outgoing ΔP  (within-patient Q4 vs Q1; paired)",
    )
    plot_ligand_bars(
        key[key["split"] == "median"] if not key.empty else key,
        args.out / "figures" / "fig_key_pairs_median.png",
        "Focused pairs  MHC-I / T-recruit / CD274–PDCD1 / NECTIN2–TIGIT  (median)",
    )
    plot_honest_n(elig, args.out / "figures" / "fig_honest_n.png")
    if not long.empty:
        plot_key_paired(long, args.out / "figures" / "fig_key_pairs_paired_median.png", "median")
        plot_key_paired(long, args.out / "figures" / "fig_key_pairs_paired_q4q1.png", "q4q1")
    plot_given_slice(units, args.out / "figures" / "fig_pr320_slice_cldn4_tnk.png")
    if not between.empty:
        plot_ligand_bars(
            between[between["split"] == "between_q4q1"],
            args.out / "figures" / "fig_between_q4q1.png",
            "Companion: between-unit Q4 vs Q1 on the PR #320 slice",
        )

    summary = {
        "datasets": ["GSE131907", "GSE205335"],
        "marker": "CLDN4",
        "additive": True,
        "dual_high": False,
        "gse207422_run": False,
        "algorithm": "Jin 2021 Hill Kh=0.5, 10% trim mean, expr_prop>=0.10, CellChatDB v2 protein",
        "cellchat_r": False,
        "liana": False,
        "unit": "patient",
        "floors": {"min_bin": MIN_BIN, "min_tnk": MIN_TNK, "min_detect_paired": MIN_DETECT_PAIRED},
        "n_lr_both_matrices": int(len(lr)),
        "n_paired_median": int(((elig.split == "median") & elig.eligible).sum()),
        "n_paired_q4q1": int(((elig.split == "q4q1") & elig.eligible).sum()),
        "n_ligand_table": int(len(paired)),
        "n_ligand_p_lt_05_median": int(((paired.split == "median") & (paired.p < 0.05)).sum()) if not paired.empty else 0,
        "n_ligand_p_lt_05_q4q1": int(((paired.split == "q4q1") & (paired.p < 0.05)).sum()) if not paired.empty else 0,
        "pr320_given": {
            "q4q1_n": 23,
            "q4q1_r": -0.705,
            "q4q1_p": 0.0003,
            "spearman_n": 43,
            "spearman_rho": -0.479,
        },
        "not_used": "FASTQ; GSE131907 2.86GB log2TPM; GSE207422; TACSTD2 groups",
        "matrix_sha256": {
            "GSE131907": sha256(args.gse131907_matrix),
            "GSE205335": sha256(args.gse205335_matrix),
        },
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(args.out, summary, paired, between, elig, key)
    print(json.dumps({k: summary[k] for k in summary if k != "matrix_sha256"}, indent=2), flush=True)
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
