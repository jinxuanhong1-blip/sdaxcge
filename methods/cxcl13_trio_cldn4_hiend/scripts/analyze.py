#!/usr/bin/env python3
"""Patient-level CellChat-style LR on the CXCL13+ trio (CLDN4-only).

Trio (given, not re-audited): GSE148071 + GSE207422 + GSE253013
n=60 ρ=−0.425 vs CXCL13+ T fraction (PR #290 / methods/scrna_cldn4_combo).

Outgoing only: CLDN4-high malignant → CXCL13+ T and → T/NK.
Gates are CLDN4-only (no TACSTD2, no dual-high). Not GSE207422-only.
CellChat R and LIANA are not run.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_MAL_ARM = 10
MIN_TNK = 15
MIN_CXCL13T = 8
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
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "NCAM1", "IFNG", "TNF", "CXCL13"]

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

GIVEN_TRIO = {
    "cohorts": ["GSE148071", "GSE207422", "GSE253013"],
    "n": 60,
    "rho": -0.425,
    "p": 0.00121,
    "source": "methods/scrna_cldn4_combo (PR #290); not re-audited",
}


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


def wanted_genes(lr: pd.DataFrame) -> set[str]:
    genes = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    for rec in lr.itertuples(index=False):
        genes.update(rec.ligand_genes)
        genes.update(rec.receptor_genes)
    return genes


def write_gene_panel(path: Path, lr: pd.DataFrame) -> list[str]:
    genes = sorted(wanted_genes(lr))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(genes) + "\n")
    return genes


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


def geom_mean_rows(mat: np.ndarray) -> np.ndarray:
    if mat.ndim == 1:
        return mat
    if mat.shape[0] == 1:
        return mat[0]
    out = np.exp(np.mean(np.log(np.clip(mat, 1e-12, None)), axis=0))
    out[np.any(mat <= 0, axis=0)] = 0.0
    return out


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def trim_mean_fast(sub: np.ndarray, proportiontocut: float = TRIM) -> np.ndarray:
    n = sub.shape[1]
    if n == 0:
        return np.zeros(sub.shape[0], dtype=np.float64)
    if n == 1:
        return sub[:, 0].astype(np.float64)
    k = int(n * proportiontocut)
    if k == 0:
        return sub.mean(axis=1)
    s = np.sort(sub, axis=1)
    return s[:, k : n - k].mean(axis=1)


def group_trim_means(expr: np.ndarray, pos: np.ndarray, labels: np.ndarray, groups: list[str]):
    n_g = expr.shape[0]
    means = np.zeros((n_g, len(groups)), dtype=np.float64)
    props = np.zeros((n_g, len(groups)), dtype=np.float64)
    counts = {}
    for j, g in enumerate(groups):
        idx = np.flatnonzero(labels == g)
        counts[g] = int(idx.size)
        if idx.size == 0:
            continue
        sub = expr[:, idx]
        means[:, j] = trim_mean_fast(sub)
        props[:, j] = pos[:, idx].mean(axis=1)
    return means, props, counts


def pair_specs(lr: pd.DataFrame, gene_index: dict[str, int]):
    specs = []
    for rec in lr.itertuples(index=False):
        if any(g not in gene_index for g in rec.ligand_genes + rec.receptor_genes):
            continue
        specs.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": rec.ligand_genes,
                "receptor_genes": rec.receptor_genes,
                "lig_ix": [gene_index[g] for g in rec.ligand_genes],
                "rec_ix": [gene_index[g] for g in rec.receptor_genes],
            }
        )
    return specs


def score_pairs(specs, gene_means, gene_props, groups, counts, src_tgt, min_src=MIN_MAL_ARM, min_tgt=MIN_TNK):
    gpos = {g: i for i, g in enumerate(groups)}
    cache = {}

    def cached(ix_key, ix):
        hit = cache.get(ix_key)
        if hit is None:
            mu = geom_mean_rows(gene_means[ix])
            pr = gene_props[ix].min(axis=0) if len(ix) > 1 else gene_props[ix[0]]
            hit = (mu, pr)
            cache[ix_key] = hit
        return hit

    rows = []
    for spec in specs:
        lig_mu, lig_pr = cached(("L", tuple(spec["lig_ix"])), spec["lig_ix"])
        rec_mu, rec_pr = cached(("R", tuple(spec["rec_ix"])), spec["rec_ix"])
        for src, tgt in src_tgt:
            if src not in gpos or tgt not in gpos:
                continue
            i, j = gpos[src], gpos[tgt]
            if counts.get(src, 0) < min_src or counts.get(tgt, 0) < min_tgt:
                continue
            detected = (float(lig_pr[i]) >= EXPR_PROP) and (float(rec_pr[j]) >= EXPR_PROP)
            prob = hill_prob(float(lig_mu[i]), float(rec_mu[j])) if detected else 0.0
            rows.append(
                {
                    "interaction_name": spec["interaction_name"],
                    "pathway_name": spec["pathway_name"],
                    "annotation": spec["annotation"],
                    "ligand": spec["ligand"],
                    "receptor": spec["receptor"],
                    "ligand_genes": "|".join(spec["ligand_genes"]),
                    "receptor_genes": "|".join(spec["receptor_genes"]),
                    "source": src,
                    "target": tgt,
                    "n_source": counts[src],
                    "n_target": counts[tgt],
                    "ligand_mean": float(lig_mu[i]),
                    "receptor_mean": float(rec_mu[j]),
                    "ligand_prop": float(lig_pr[i]),
                    "receptor_prop": float(rec_pr[j]),
                    "detected": bool(detected),
                    "prob": prob,
                }
            )
    return pd.DataFrame(rows)


def ligand_class(genes: str) -> str:
    parts = set(str(genes).split("|"))
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


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    try:
        out[ok] = false_discovery_control(p[ok], method="bh")
    except Exception:
        # fallback
        order = np.argsort(p[ok])
        m = ok.sum()
        ranked = p[ok][order]
        q = ranked * m / (np.arange(1, m + 1))
        q = np.minimum.accumulate(q[::-1])[::-1]
        tmp = np.empty(m)
        tmp[order] = np.clip(q, 0, 1)
        out[ok] = tmp
    return out


# ---------------------------------------------------------------------------
# Cohort loaders
# ---------------------------------------------------------------------------

def patient_from_filename(path: Path) -> str:
    m = re.search(r"_(P\d+)_", path.name)
    if m:
        return m.group(1)
    stem = path.name.replace(".txt.gz", "")
    parts = stem.split("_")
    return parts[1] if len(parts) >= 2 else stem


def list_exp_files(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*_exp.txt.gz"))
    if files:
        return files
    tar_path = data_dir / "GSE148071_RAW.tar"
    if not tar_path.exists():
        raise SystemExit(f"no exp matrices and no {tar_path}")
    dest = data_dir / "files"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    files = sorted(dest.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz in {dest}")
    return files


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
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
    return cell_ids, found, n_umi, n_streamed


def stream_gse207422(matrix_path: Path, wanted: set[str]):
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
            gene = gene.strip().strip('"').split(".")[0]
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


def finish_cohort(cohort: str, sample: np.ndarray, found: dict[str, np.ndarray], n_umi: np.ndarray, extra_meta=None):
    n = len(sample)
    lib = np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
    lineage = assign_lineage(scores, cd3)
    cxcl = found.get("CXCL13", np.zeros(n, dtype=np.float32))
    is_t = lineage == "T"
    is_cxcl13t = is_t & (cxcl > 0)
    is_tnk = np.isin(lineage, ["T", "NK"])
    is_mal = lineage == "Epithelial"
    rec = {
        "cohort": cohort,
        "sample": sample,
        "found": found,
        "n_umi": n_umi,
        "log_cp": log_cp,
        "lineage": lineage,
        "is_mal": is_mal,
        "is_tnk": is_tnk,
        "is_cxcl13t": is_cxcl13t,
        "n": n,
    }
    if extra_meta:
        rec.update(extra_meta)
    return rec


def load_gse148071(data_dir: Path, wanted: set[str]):
    files = list_exp_files(data_dir)
    print(f"GSE148071 exp matrices: {len(files)}", flush=True)
    cell_ids_all: list[str] = []
    sample_all: list[str] = []
    found_all: dict[str, list[np.ndarray]] = {g: [] for g in wanted}
    n_umi_all: list[np.ndarray] = []
    for fp in files:
        patient = patient_from_filename(fp)
        cell_ids, found, n_umi, n_streamed = stream_one_matrix(fp, wanted)
        n = len(cell_ids)
        cell_ids_all.extend([f"{patient}_{c}" for c in cell_ids])
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        for g in wanted:
            found_all[g].append(found[g] if g in found else np.zeros(n, dtype=np.float32))
        print(f"  {fp.name} {patient} n={n} wanted={len(found)} genes={n_streamed}", flush=True)
    sample = np.array(sample_all, dtype=object)
    n_umi = np.concatenate(n_umi_all)
    found = {g: np.concatenate(v) for g, v in found_all.items()}
    return finish_cohort("GSE148071", sample, found, n_umi)


def load_gse207422(data_dir: Path, wanted: set[str]):
    matrix = data_dir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_path = data_dir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    print(f"GSE207422 streaming {matrix}", flush=True)
    cell_ids, found, n_umi, n_streamed = stream_gse207422(matrix, wanted)
    sample = np.array(["_".join(bc.split("_")[:2]) for bc in cell_ids], dtype=object)
    extra = {}
    if meta_path.exists():
        meta = pd.read_excel(meta_path).dropna(subset=["Sample"]).copy()
        meta["Sample"] = meta["Sample"].astype(str)
        smap = dict(zip(meta["Sample"], meta["Patient"].astype(str)))
        rmap = dict(zip(meta["Sample"], meta["Resource"].astype(str)))
        pmap = dict(zip(meta["Sample"], meta["Pathologic Response"].astype(str)))
        patient = np.array([smap.get(s, s) for s in sample], dtype=object)
        extra = {
            "patient_id": patient,
            "resource": np.array([rmap.get(s, "") for s in sample], dtype=object),
            "response": np.array([pmap.get(s, "") for s in sample], dtype=object),
        }
        sample = patient
    rec = finish_cohort("GSE207422", sample, found, n_umi, extra)
    rec["n_streamed"] = n_streamed
    return rec


def load_gse253013(extracted: Path, wanted: set[str]):
    npz = extracted / "gene_panel.npz"
    meta_path = extracted / "cell_metadata.tsv"
    if not npz.exists():
        raise FileNotFoundError(f"missing {npz}; run extract + assemble first")
    print(f"GSE253013 loading {npz}", flush=True)
    data = np.load(npz)
    n = int(data["total"].shape[0]) if "total" in data.files else int(data[data.files[0]].shape[0])
    found = {}
    for g in wanted:
        if g in data.files:
            found[g] = np.asarray(data[g], dtype=np.float32)
        else:
            found[g] = np.zeros(n, dtype=np.float32)
    n_umi = np.asarray(data["total"], dtype=np.float64) if "total" in data.files else np.ones(n)
    meta = pd.read_csv(meta_path, sep="\t") if meta_path.exists() else pd.DataFrame()
    if "patient" in meta.columns:
        sample = meta["patient"].astype(str).to_numpy()
    elif "orig.ident" in meta.columns:
        sample = meta["orig.ident"].astype(str).to_numpy()
    else:
        sample = np.array(["MRC_unknown"] * n, dtype=object)
    tissue = meta["tissue"].astype(str).to_numpy() if "tissue" in meta.columns else np.array(["Tumor"] * n)
    rec = finish_cohort("GSE253013", sample, found, n_umi, {"tissue": tissue})
    # Tumor only — adjacent lung is not the CXCL13+ trio unit
    if "tissue" in rec:
        keep = rec["tissue"] == "Tumor"
        if keep.sum() < rec["n"] and keep.sum() > 1000:
            print(f"  GSE253013 keep Tumor {int(keep.sum())}/{rec['n']}", flush=True)
            rec = subset_rec(rec, keep)
    # Prefer author epithelial/malignant if present
    if "cell_type" in meta.columns:
        ct = meta["cell_type"].astype(str).to_numpy()
        if "tissue" in rec and rec["n"] != len(ct):
            ct = ct[rec.get("_keep_idx", np.arange(len(ct)))]
        if len(ct) == rec["n"]:
            rec["author_cell_type"] = ct
            epi_like = np.array(
                [re.search(r"epithel|malign|tumor|cancer|aluad|at2|club", str(x), re.I) is not None for x in ct]
            )
            if epi_like.sum() >= 50:
                rec["is_mal"] = epi_like
                print(f"  GSE253013 author-like epithelial/malignant n={int(epi_like.sum())}", flush=True)
    return rec


def subset_rec(rec: dict, keep: np.ndarray) -> dict:
    out = dict(rec)
    out["n"] = int(keep.sum())
    out["sample"] = rec["sample"][keep]
    out["n_umi"] = rec["n_umi"][keep]
    out["lineage"] = rec["lineage"][keep]
    out["is_mal"] = rec["is_mal"][keep]
    out["is_tnk"] = rec["is_tnk"][keep]
    out["is_cxcl13t"] = rec["is_cxcl13t"][keep]
    out["found"] = {g: v[keep] for g, v in rec["found"].items()}
    out["log_cp"] = {g: v[keep] for g, v in rec["log_cp"].items()}
    for k in ("patient_id", "resource", "response", "tissue", "author_cell_type"):
        if k in rec:
            out[k] = rec[k][keep]
    out["_keep_idx"] = np.flatnonzero(keep)
    return out


# ---------------------------------------------------------------------------
# Patient-level scoring
# ---------------------------------------------------------------------------

def score_patient(expr, pos, specs, labels, groups, src_tgt, min_src, min_tgt):
    means, props, counts = group_trim_means(expr, pos, labels, groups)
    return score_pairs(specs, means, props, groups, counts, src_tgt, min_src=min_src, min_tgt=min_tgt)


def patient_rows_from_pairs(pairs: pd.DataFrame, target: str) -> pd.DataFrame:
    if pairs.empty:
        return pd.DataFrame()
    high = pairs[(pairs.source == "Mal_high") & (pairs.target == target)].set_index("interaction_name")
    low = pairs[(pairs.source == "Mal_low") & (pairs.target == target)].set_index("interaction_name")
    common = high.index.intersection(low.index)
    if len(common) == 0:
        return pd.DataFrame()
    out = high.loc[common, ["pathway_name", "annotation", "ligand", "receptor", "ligand_genes", "receptor_genes"]].copy()
    out["target"] = target
    out["prob_high"] = high.loc[common, "prob"].to_numpy()
    out["prob_low"] = low.loc[common, "prob"].to_numpy()
    out["delta_prob"] = out["prob_high"] - out["prob_low"]
    out["detected_high"] = high.loc[common, "detected"].to_numpy()
    out["detected_low"] = low.loc[common, "detected"].to_numpy()
    out["n_mal_high"] = high.loc[common, "n_source"].to_numpy()
    out["n_mal_low"] = low.loc[common, "n_source"].to_numpy()
    out["n_target"] = high.loc[common, "n_target"].to_numpy()
    out["ligand_class"] = [ligand_class(x) for x in out["ligand_genes"]]
    return out.reset_index()


def run_patient_level(rec: dict, specs, lr_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    log_cp = rec["log_cp"]
    found = rec["found"]
    sample = rec["sample"]
    gene_index = {g: i for i, g in enumerate(lr_genes)}
    expr = np.vstack([log_cp.get(g, np.zeros(rec["n"], dtype=np.float32)) for g in lr_genes])
    pos = np.vstack([(found.get(g, np.zeros(rec["n"], dtype=np.float32)) > 0).astype(np.float32) for g in lr_genes])
    marker = log_cp.get(MARKER, np.zeros(rec["n"], dtype=np.float32))

    n_rows = []
    pair_rows = []
    for pat in sorted(set(sample.tolist()), key=lambda x: (len(str(x)), str(x))):
        m = sample == pat
        mal = rec["is_mal"] & m
        tnk = rec["is_tnk"] & m
        cx = rec["is_cxcl13t"] & m
        n_mal = int(mal.sum())
        n_tnk = int(tnk.sum())
        n_cx = int(cx.sum())
        n_t = int(((rec["lineage"] == "T") & m).sum())
        cldn = marker[mal]
        med = float(np.median(cldn)) if n_mal else float("nan")
        high = mal & (marker > med) if n_mal else mal
        low = mal & (marker <= med) if n_mal else mal
        n_hi = int(high.sum())
        n_lo = int(low.sum())
        elig_tnk = n_hi >= MIN_MAL_ARM and n_lo >= MIN_MAL_ARM and n_tnk >= MIN_TNK
        elig_cx = n_hi >= MIN_MAL_ARM and n_lo >= MIN_MAL_ARM and n_cx >= MIN_CXCL13T
        n_rows.append(
            {
                "cohort": rec["cohort"],
                "patient": pat,
                "n_total": int(m.sum()),
                "n_epithelial": n_mal,
                "n_mal_high": n_hi,
                "n_mal_low": n_lo,
                "n_TNK": n_tnk,
                "n_T": n_t,
                "n_CXCL13pos_T": n_cx,
                "cldn4_median": med,
                "mean_CLDN4_epithelial": float(cldn.mean()) if n_mal else None,
                "eligible_TNK": elig_tnk,
                "eligible_CXCL13T": elig_cx,
                "gate": "within-patient median CLDN4 log1p(CP10k); CLDN4-only",
            }
        )
        if not (elig_tnk or elig_cx):
            continue
        labels = np.array(["drop"] * rec["n"], dtype=object)
        labels[high] = "Mal_high"
        labels[low] = "Mal_low"
        src_tgt = []
        min_tgt = {}
        if elig_tnk:
            labels[tnk] = "TNK"
            src_tgt += [("Mal_high", "TNK"), ("Mal_low", "TNK")]
        if elig_cx:
            labels[cx] = "CXCL13T"
            src_tgt += [("Mal_high", "CXCL13T"), ("Mal_low", "CXCL13T")]
        used = labels != "drop"
        groups = sorted(set(labels[used]))
        pairs = score_patient(
            expr[:, used],
            pos[:, used],
            specs,
            labels[used],
            groups,
            src_tgt,
            MIN_MAL_ARM,
            min(MIN_TNK, MIN_CXCL13T),
        )
        if pairs.empty:
            continue
        for target, elig in (("TNK", elig_tnk), ("CXCL13T", elig_cx)):
            if not elig:
                continue
            wide = patient_rows_from_pairs(pairs, target)
            if wide.empty:
                continue
            wide.insert(0, "cohort", rec["cohort"])
            wide.insert(1, "patient", pat)
            pair_rows.append(wide)
        print(
            f"  {rec['cohort']} {pat} mal={n_mal} hi/lo={n_hi}/{n_lo} "
            f"TNK={n_tnk} CXCL13T={n_cx} tnk={elig_tnk} cx={elig_cx}",
            flush=True,
        )
    n_df = pd.DataFrame(n_rows)
    p_df = pd.concat(pair_rows, ignore_index=True) if pair_rows else pd.DataFrame()
    return n_df, p_df


def pool_pairs(patient_pairs: pd.DataFrame) -> pd.DataFrame:
    if patient_pairs.empty:
        return pd.DataFrame()
    rows = []
    keys = ["target", "interaction_name", "pathway_name", "annotation", "ligand", "receptor", "ligand_genes", "receptor_genes", "ligand_class"]
    for key, g in patient_pairs.groupby(keys, dropna=False):
        rec = dict(zip(keys, key))
        d = g["delta_prob"].to_numpy(dtype=float)
        # patients where the pair is detected on at least one arm
        det = (g["detected_high"] | g["detected_low"]).to_numpy()
        d_det = d[det]
        rec["n_patients_scored"] = int(len(g))
        rec["n_patients_detected"] = int(det.sum())
        rec["n_cohorts"] = int(g["cohort"].nunique())
        rec["cohorts"] = "+".join(sorted(g["cohort"].unique()))
        rec["median_prob_high"] = float(np.median(g["prob_high"]))
        rec["median_prob_low"] = float(np.median(g["prob_low"]))
        rec["median_delta"] = float(np.median(d))
        rec["mean_delta"] = float(np.mean(d))
        rec["n_delta_pos"] = int((d > 0).sum())
        rec["n_delta_neg"] = int((d < 0).sum())
        rec["n_delta_zero"] = int((d == 0).sum())
        rec["frac_delta_pos"] = float((d > 0).mean())
        p_w = np.nan
        if d_det.size >= 6 and np.any(d_det != 0):
            try:
                p_w = float(wilcoxon(d_det, zero_method="wilcox", alternative="two-sided").pvalue)
            except ValueError:
                p_w = np.nan
        rec["wilcoxon_p"] = p_w
        rec["n_wilcoxon"] = int(d_det.size)
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr_bh"] = np.nan
    for tgt, idx in out.groupby("target").groups.items():
        out.loc[list(idx), "fdr_bh"] = bh_fdr(out.loc[list(idx), "wilcoxon_p"].to_numpy())
    out = out.sort_values(["target", "median_delta"], ascending=[True, False])
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_n(n_df: pd.DataFrame, path: Path, title: str) -> None:
    if n_df.empty:
        return
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=False)
    cohorts = list(n_df["cohort"].unique())
    colors = {"GSE148071": "#4C78A8", "GSE207422": "#F58518", "GSE253013": "#54A24B"}
    for ax, col, ylab in zip(
        axes,
        ["n_epithelial", "n_TNK", "n_CXCL13pos_T"],
        ["epithelial (putative mal.)", "T/NK", "CXCL13+ T"],
    ):
        x = 0
        ticks = []
        labs = []
        for c in cohorts:
            sub = n_df[n_df.cohort == c].reset_index(drop=True)
            xs = np.arange(len(sub)) + x
            ax.bar(xs, sub[col], color=colors.get(c, "grey"), width=0.8)
            ticks.extend(xs.tolist())
            labs.extend(sub["patient"].astype(str).tolist())
            x += len(sub) + 1
        ax.set_ylabel(ylab)
        ax.axhline(MIN_MAL_ARM if col == "n_epithelial" else (MIN_TNK if col == "n_TNK" else MIN_CXCL13T), color="k", ls="--", lw=0.8)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labs, rotation=90, fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_eligible(n_df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    rows = []
    for c, g in n_df.groupby("cohort"):
        rows.append((c, "deposited", len(g)))
        rows.append((c, "T/NK eligible", int(g.eligible_TNK.sum())))
        rows.append((c, "CXCL13+ T eligible", int(g.eligible_CXCL13T.sum())))
    df = pd.DataFrame(rows, columns=["cohort", "set", "n"])
    sets = ["deposited", "T/NK eligible", "CXCL13+ T eligible"]
    cohorts = list(n_df["cohort"].unique())
    x = np.arange(len(cohorts))
    w = 0.25
    for i, s in enumerate(sets):
        vals = [int(df[(df.cohort == c) & (df.set == s)]["n"].sum()) for c in cohorts]
        ax.bar(x + (i - 1) * w, vals, w, label=s)
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts)
    ax.set_ylabel("patients")
    ax.set_title("Honest n — CLDN4-only within-patient median gate")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_top_delta(table: pd.DataFrame, path: Path, title: str, k: int = 18) -> None:
    if table.empty:
        return
    sig = table[table["wilcoxon_p"].notna() & (table["n_wilcoxon"] >= 6)].copy()
    if sig.empty:
        sig = table.copy()
    sig = sig.reindex(sig["median_delta"].abs().sort_values(ascending=False).index).head(k)
    sig = sig.sort_values("median_delta")
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.32 * len(sig) + 1.4)))
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in sig["median_delta"]]
    ax.barh(np.arange(len(sig)), sig["median_delta"], color=colors)
    ax.set_yticks(np.arange(len(sig)))
    ax.set_yticklabels([f"{a}–{b}" for a, b in zip(sig["ligand"], sig["receptor"])], fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("median ΔP (CLDN4-high − CLDN4-low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(table: pd.DataFrame, path: Path, title: str, k: int = 16) -> None:
    if table.empty:
        return
    show = table[table["n_wilcoxon"] >= 6].copy()
    if show.empty:
        show = table.copy()
    show = show.reindex(show["median_delta"].abs().sort_values(ascending=False).index).head(k)
    fig, ax = plt.subplots(figsize=(11.5, 0.42 * len(show) + 1.8))
    ax.axis("off")
    cols = ["target", "ligand", "receptor", "ligand_class", "n_wilcoxon", "median_delta", "wilcoxon_p", "n_delta_pos"]
    cell = show[cols].copy()
    cell["median_delta"] = cell["median_delta"].map(lambda x: f"{x:+.3f}")
    cell["wilcoxon_p"] = cell["wilcoxon_p"].map(lambda x: "NA" if pd.isna(x) else f"{x:.3g}")
    ax.table(cellText=cell.values, colLabels=cols, loc="center", cellLoc="center")
    ax.set_title(title, pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest_key(patient_pairs: pd.DataFrame, names: list[str], target: str, path: Path, title: str) -> None:
    sub = patient_pairs[(patient_pairs.target == target) & (patient_pairs.interaction_name.isin(names))]
    if sub.empty:
        return
    fig, axes = plt.subplots(len(names), 1, figsize=(8.5, 2.2 * len(names)), sharex=True)
    if len(names) == 1:
        axes = [axes]
    colors = {"GSE148071": "#4C78A8", "GSE207422": "#F58518", "GSE253013": "#54A24B"}
    for ax, name in zip(axes, names):
        g = sub[sub.interaction_name == name]
        if g.empty:
            ax.set_title(name + " (absent)")
            continue
        y = 0
        yt, yl = [], []
        for c, gg in g.groupby("cohort"):
            xs = gg["delta_prob"].to_numpy()
            jitter = np.random.default_rng(1).uniform(-0.12, 0.12, size=len(xs))
            ax.scatter(xs, np.full(len(xs), y) + jitter, s=18, color=colors.get(c, "grey"), alpha=0.8)
            ax.plot([np.median(xs)], [y], "k|", ms=14)
            yt.append(y)
            yl.append(f"{c} n={len(gg)}")
            y += 1
        ax.axvline(0, color="k", lw=0.7)
        ax.set_yticks(yt)
        ax.set_yticklabels(yl, fontsize=8)
        ax.set_title(name, fontsize=10)
    axes[-1].set_xlabel("per-patient ΔP (high − low)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# FINDING
# ---------------------------------------------------------------------------

def md_table(df: pd.DataFrame, cols: list[str], fmt: dict | None = None) -> str:
    fmt = fmt or {}
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for rec in df.itertuples(index=False):
        row = []
        for c in cols:
            v = getattr(rec, c)
            if c in fmt and pd.notna(v):
                row.append(fmt[c].format(v))
            elif pd.isna(v):
                row.append("NA")
            else:
                row.append(str(v))
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def write_finding(n_df: pd.DataFrame, pooled: pd.DataFrame, patient_pairs: pd.DataFrame, out_dir: Path, missing: list[str]) -> None:
    n_dep = n_df.groupby("cohort").size().to_dict()
    n_tnk = n_df.groupby("cohort")["eligible_TNK"].sum().to_dict()
    n_cx = n_df.groupby("cohort")["eligible_CXCL13T"].sum().to_dict()
    n_tnk_tot = int(n_df["eligible_TNK"].sum())
    n_cx_tot = int(n_df["eligible_CXCL13T"].sum())
    cohorts_done = sorted(n_df["cohort"].unique())

    def top_block(target: str, sign: str, k: int = 12) -> pd.DataFrame:
        g = pooled[pooled.target == target].copy()
        if g.empty:
            return g
        g = g[g["n_wilcoxon"] >= 6]
        if sign == "pos":
            g = g[g["median_delta"] > 0]
            return g.sort_values(["wilcoxon_p", "median_delta"], ascending=[True, False]).head(k)
        g = g[g["median_delta"] < 0]
        return g.sort_values(["wilcoxon_p", "median_delta"], ascending=[True, True]).head(k)

    cols = ["interaction_name", "pathway_name", "ligand_class", "n_wilcoxon", "n_cohorts", "median_delta", "wilcoxon_p", "n_delta_pos"]
    fmt = {"median_delta": "{:+.3f}", "wilcoxon_p": "{:.3g}", "n_wilcoxon": "{:.0f}", "n_cohorts": "{:.0f}", "n_delta_pos": "{:.0f}"}

    tnk_up = top_block("TNK", "pos")
    tnk_dn = top_block("TNK", "neg")
    cx_up = top_block("CXCL13T", "pos")
    cx_dn = top_block("CXCL13T", "neg")

    missing_txt = ", ".join(missing) if missing else "none"
    text = f"""# FINDING — CellChat-style outgoing CLDN4-high malignant on the CXCL13+ trio

**CLDN4 only. No dual-high. Not GSE207422-only.** Patient is the unit.

The CXCL13+ trio that already differs is taken as given and is **not re-audited**:
GSE148071 + GSE207422 + GSE253013, **n=60, ρ=−0.425, p=0.00121**
(`methods/scrna_cldn4_combo`, PR #290, family `tls/cxcl13pos/mean`).

This folder adds **outgoing CellChat-style probability** from **CLDN4-high vs CLDN4-low**
putative malignant epithelium toward **CXCL13+ T** and toward **T/NK**, scored
**within each patient** (median CLDN4 gate), then tested across patients
(Wilcoxon on ΔP). CellChat R and LIANA were not run.

---

## Honest n (communication, not the given Spearman n=60)

Given combo n=60 is mean-CLDN4 vs CXCL13+ T *fraction*. Communication n is the
number of patients with enough cells on **both** CLDN4 arms **and** the receiver.

| Cohort | Deposited (this extract) | Eligible Mal→T/NK | Eligible Mal→CXCL13+ T |
| --- | ---: | ---: | ---: |
| GSE148071 | {n_dep.get("GSE148071", 0)} | {int(n_tnk.get("GSE148071", 0))} | {int(n_cx.get("GSE148071", 0))} |
| GSE207422 | {n_dep.get("GSE207422", 0)} | {int(n_tnk.get("GSE207422", 0))} | {int(n_cx.get("GSE207422", 0))} |
| GSE253013 | {n_dep.get("GSE253013", 0)} | {int(n_tnk.get("GSE253013", 0))} | {int(n_cx.get("GSE253013", 0))} |
| **Trio** | **{int(n_df.shape[0])}** | **{n_tnk_tot}** | **{n_cx_tot}** |

Floors: ≥{MIN_MAL_ARM} CLDN4-high and ≥{MIN_MAL_ARM} CLDN4-low epithelial cells
(within-patient median of `log1p(CP10k)` CLDN4); ≥{MIN_TNK} T/NK; ≥{MIN_CXCL13T} CXCL13+ T
(T lineage and CXCL13 UMI > 0). Missing this run: **{missing_txt}**.
Cohorts scored: {", ".join(cohorts_done) if cohorts_done else "none"}.

Do **not** read the LR table as n=60. Patients below the floor stay in
`results/patient_n.tsv` and are not scored.

## Gate

- Marker: **CLDN4 only**. TACSTD2 is recorded in the stream and is not used.
- No dual-high TACSTD2×CLDN4 quadrant.
- Split: **within-patient median** among epithelial cells (putative malignant;
  marker-argmax; GSE253013 uses author epithelial/malignant when present).
- Probability: Jin et al. 2021 CellChat Hill / mass-action on 10% truncated
  means, CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`, Kh=0.5.
- Test: Wilcoxon signed-rank on per-patient ΔP = P(high→receiver) − P(low→receiver),
  among patients where the pair is detected on at least one arm (n≥6).
  BH-FDR within receiver. p-values are descriptive.

## Ligand table — outgoing CLDN4-high → CXCL13+ T

Full table: `results/ligand_table.tsv`.

**Higher from CLDN4-high (median ΔP > 0):**

{md_table(cx_up, cols, fmt) if not cx_up.empty else "_no pair with n_wilcoxon≥6 and median ΔP>0_"}

**Higher from CLDN4-low (median ΔP < 0):**

{md_table(cx_dn, cols, fmt) if not cx_dn.empty else "_no pair with n_wilcoxon≥6 and median ΔP<0_"}

## Ligand table — outgoing CLDN4-high → T/NK

**Higher from CLDN4-high (median ΔP > 0):**

{md_table(tnk_up, cols, fmt) if not tnk_up.empty else "_no pair with n_wilcoxon≥6 and median ΔP>0_"}

**Higher from CLDN4-low (median ΔP < 0):**

{md_table(tnk_dn, cols, fmt) if not tnk_dn.empty else "_no pair with n_wilcoxon≥6 and median ΔP<0_"}

## What is not supported

- Re-using **n=60** as the communication n.
- A GSE207422-only CellChat (already in `methods/scrna_cellchat_cldn4`).
- Dual-high TACSTD2×CLDN4 gates.
- Cell-pooled permutation as a patient claim.
- ICI response or histology as the unit of this table.

## Files

`results/ligand_table.tsv` (pooled patient-level LR),
`results/lr_pairs_patient.tsv`, `results/patient_n.tsv`,
`results/fig_n_per_patient.png`, `results/fig_n_eligible.png`,
`results/fig_extra_ligand_table.png`, `results/fig_top_outgoing_CXCL13T.png`,
`results/fig_top_outgoing_TNK.png`, `results/summary.json`.
"""
    (out_dir.parent / "FINDING.md").write_text(text)
    (out_dir / "FINDING.md").write_text(text)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/cxcl13_trio_data"))
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT / "results")
    p.add_argument("--cohorts", default="GSE148071,GSE207422,GSE253013")
    p.add_argument("--gse253013-extracted", type=Path, default=Path("/tmp/cxcl13_trio_data/GSE253013/extracted"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (ROOT / "figures").mkdir(parents=True, exist_ok=True)

    want_cohorts = [c.strip() for c in args.cohorts.split(",") if c.strip()]
    lr_all = load_lr(args.db, None)
    wanted = wanted_genes(lr_all)
    write_gene_panel(ROOT / "db" / "gene_panel.txt", lr_all)
    print(f"CellChatDB protein pairs={len(lr_all)} wanted_genes={len(wanted)}", flush=True)

    n_frames = []
    p_frames = []
    missing = []
    loaded = {}

    loaders = {
        "GSE148071": lambda: load_gse148071(args.data / "GSE148071", wanted),
        "GSE207422": lambda: load_gse207422(args.data / "GSE207422", wanted),
        "GSE253013": lambda: load_gse253013(args.gse253013_extracted, wanted),
    }

    for cohort in want_cohorts:
        try:
            print(f"== load {cohort} ==", flush=True)
            rec = loaders[cohort]()
            loaded[cohort] = rec
            print(
                f"  {cohort} cells={rec['n']} epi={int(rec['is_mal'].sum())} "
                f"TNK={int(rec['is_tnk'].sum())} CXCL13T={int(rec['is_cxcl13t'].sum())}",
                flush=True,
            )
        except Exception as exc:
            print(f"SKIP {cohort}: {exc}", flush=True)
            missing.append(f"{cohort} ({exc})")

    # Shared LR universe = genes present in at least one loaded cohort
    present = set()
    for rec in loaded.values():
        present.update(g for g, v in rec["found"].items() if float(np.max(v)) > 0)
    lr = load_lr(args.db, present)
    lr_genes = sorted({g for rec in lr.itertuples(index=False) for g in rec.ligand_genes + rec.receptor_genes})
    specs = pair_specs(lr, {g: i for i, g in enumerate(lr_genes)})
    print(f"LR pairs with subunits in union={len(lr)} genes={len(lr_genes)} specs={len(specs)}", flush=True)

    for cohort, rec in loaded.items():
        print(f"== patient-level {cohort} ==", flush=True)
        n_df, p_df = run_patient_level(rec, specs, lr_genes)
        n_df.to_csv(args.out / f"patient_n_{cohort}.tsv", sep="\t", index=False)
        if not p_df.empty:
            p_df.to_csv(args.out / f"lr_pairs_patient_{cohort}.tsv", sep="\t", index=False)
        n_frames.append(n_df)
        if not p_df.empty:
            p_frames.append(p_df)
        # free big arrays
        del rec
        loaded[cohort] = None

    n_all = pd.concat(n_frames, ignore_index=True) if n_frames else pd.DataFrame()
    p_all = pd.concat(p_frames, ignore_index=True) if p_frames else pd.DataFrame()
    n_all.to_csv(args.out / "patient_n.tsv", sep="\t", index=False)
    if not p_all.empty:
        p_all.to_csv(args.out / "lr_pairs_patient.tsv", sep="\t", index=False)

    pooled = pool_pairs(p_all)
    if not pooled.empty:
        pooled.to_csv(args.out / "ligand_table.tsv", sep="\t", index=False)
        pooled.to_csv(ROOT / "figures" / "ligand_table.tsv", sep="\t", index=False)

    if not n_all.empty:
        plot_n(n_all, args.out / "fig_n_per_patient.png", "Trio cell counts by patient (honest n)")
        plot_eligible(n_all, args.out / "fig_n_eligible.png")
        for src, dst in (
            (args.out / "fig_n_per_patient.png", ROOT / "figures" / "fig_n_per_patient.png"),
            (args.out / "fig_n_eligible.png", ROOT / "figures" / "fig_n_eligible.png"),
        ):
            if src.exists():
                dst.write_bytes(src.read_bytes())

    if not pooled.empty:
        for tgt in ("CXCL13T", "TNK"):
            sub = pooled[pooled.target == tgt]
            plot_top_delta(
                sub,
                args.out / f"fig_top_outgoing_{tgt}.png",
                f"Outgoing Mal CLDN4-high → {tgt} median ΔP (patient-level)",
            )
        plot_ligand_table(
            pooled,
            args.out / "fig_extra_ligand_table.png",
            "Extra: patient-level outgoing LR (top |median ΔP|)",
        )
        for fn in ("fig_top_outgoing_CXCL13T.png", "fig_top_outgoing_TNK.png", "fig_extra_ligand_table.png"):
            src = args.out / fn
            if src.exists():
                (ROOT / "figures" / fn).write_bytes(src.read_bytes())

    if not p_all.empty:
        key = [
            "NECTIN2_TIGIT",
            "LGALS9_PTPRC",
            "CXCL16_CXCR6",
            "HLA-E_CD8A",
            "CD274_PDCD1",
            "F11R_ITGAL_ITGB2",
            "MDK_NCL",
        ]
        have = [x for x in key if x in set(p_all.interaction_name)]
        if have:
            plot_forest_key(p_all, have[:5], "TNK", args.out / "fig_keypairs_TNK.png", "Key pairs Mal→T/NK ΔP")
            plot_forest_key(p_all, have[:5], "CXCL13T", args.out / "fig_keypairs_CXCL13T.png", "Key pairs Mal→CXCL13+ T ΔP")

    write_finding(n_all, pooled, p_all, args.out, missing)

    summary = {
        "given_trio": GIVEN_TRIO,
        "gate": "CLDN4-only within-patient median of log1p(CP10k) among epithelial cells",
        "no_dual_high": True,
        "not_gse207422_only": True,
        "min_mal_arm": MIN_MAL_ARM,
        "min_tnk": MIN_TNK,
        "min_cxcl13t": MIN_CXCL13T,
        "kh": KH,
        "trim": TRIM,
        "expr_prop": EXPR_PROP,
        "n_deposited": n_all.groupby("cohort").size().to_dict() if not n_all.empty else {},
        "n_eligible_TNK": int(n_all["eligible_TNK"].sum()) if not n_all.empty else 0,
        "n_eligible_CXCL13T": int(n_all["eligible_CXCL13T"].sum()) if not n_all.empty else 0,
        "n_lr_rows_patient": int(len(p_all)),
        "n_lr_rows_pooled": int(len(pooled)),
        "missing": missing,
        "not_run": "CellChat R; LIANA; TACSTD2 split; dual-high; GSE207422-only re-run; given n=60 Spearman re-audit",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("wrote", args.out / "ligand_table.tsv" if (args.out / "ligand_table.tsv").exists() else "NO ligand table")


if __name__ == "__main__":
    main()
