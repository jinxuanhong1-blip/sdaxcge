#!/usr/bin/env python3
"""Patient-level extract for GSE123902, GSE131907, and GSE189357.

GSE205335 is extracted in R (RDS). Discordant accessions are not read.
Unit definitions match the locked concordant-4 tables:
  GSE123902 donor, primary preferred over metastasis, normal dropped
  GSE131907 sample (21 malignant-bearing mLN / tL/B / mBrain)
  GSE189357 patient, marker malignant
"""

from __future__ import annotations

import gzip
import os
import tarfile
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

RAW = os.environ.get("CONCORDANT4_RAW", "/tmp/concordant4_raw")
HERE = os.environ.get(
    "NICHENET_HERE", "/workspace/methods/nichenet_concordant4_cldn4"
)
OUT = os.environ.get("NICHENET_EXTRACT", "/tmp/nichenet_work/extract")
PRIOR = os.environ.get("NICHENET_PRIOR", "/tmp/nichenet_prior")
os.makedirs(OUT, exist_ok=True)

ALIASES = {"PVRL2": "NECTIN2", "JAM1": "F11R", "PVRL1": "NECTIN1", "PVRL3": "NECTIN3"}
EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")
TNK = ("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
MIN_MAL = 40
MIN_TNK = 20
MALIG_SUB = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def canon(g: str) -> str:
    g = str(g).upper()
    return ALIASES.get(g, g)


def load_ligands() -> list[str]:
    with open(os.path.join(PRIOR, "ligands.txt")) as f:
        return [canon(x.strip()) for x in f if x.strip()]


def quartile_masks(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    n = x.size
    if n < 4:
        z = np.zeros(n, dtype=bool)
        return z, z, False
    r = rankdata(x, method="ordinal")  # R rank(..., ties.method="first")
    q1 = int(np.floor(n * 0.25))
    q4 = int(np.ceil(n * 0.75))
    if q1 < 1 or q4 >= n:
        z = np.zeros(n, dtype=bool)
        return z, z, False
    return r > q4, r <= q1, True


def group_stats(counts: np.ndarray, lib: np.ndarray) -> tuple[float, float]:
    """Mean log1p(CP10k) and fraction UMI>0. counts and lib aligned."""
    if counts.size == 0:
        return np.nan, np.nan
    ok = lib > 0
    if not np.any(ok):
        return np.nan, np.nan
    c = counts[ok]
    ell = lib[ok]
    return float(np.log1p(1e4 * c / ell).mean()), float(np.mean(c > 0))


def splits_from_mal(cldn4_umi: np.ndarray, lib: np.ndarray) -> dict:
    logc = np.log1p(1e4 * cldn4_umi / np.maximum(lib, 1.0))
    hi, lo, ok = quartile_masks(logc)
    med = float(np.median(logc)) if logc.size else np.nan
    return {
        "log": logc,
        "q_hi": hi,
        "q_lo": lo,
        "q_ok": ok,
        "med_hi": logc > med if logc.size else hi,
        "med_lo": logc <= med if logc.size else lo,
        "pos": cldn4_umi > 0,
        "neg": cldn4_umi == 0,
        "pct": float(np.mean(cldn4_umi > 0)) if cldn4_umi.size else np.nan,
        "mean": float(logc.mean()) if logc.size else np.nan,
        "mean_hi": float(logc[hi].mean()) if np.any(hi) else np.nan,
        "mean_lo": float(logc[lo].mean()) if np.any(lo) else np.nan,
    }


def ligand_rows(cohort, unit, ligands, gene_index, counts_by_gene, mal_lib, sp) -> list[dict]:
    """counts_by_gene: dict gene -> malignant-cell count vector aligned to mal_lib.

    Ligands absent from the matrix are omitted (not filled with zeros).
    """
    rows = []
    n = mal_lib.size
    for lig in ligands:
        vec = counts_by_gene.get(lig)
        if vec is None:
            continue
        m_q4, _ = group_stats(vec[sp["q_hi"]], mal_lib[sp["q_hi"]])
        m_q1, _ = group_stats(vec[sp["q_lo"]], mal_lib[sp["q_lo"]])
        m_mh, _ = group_stats(vec[sp["med_hi"]], mal_lib[sp["med_hi"]])
        m_ml, _ = group_stats(vec[sp["med_lo"]], mal_lib[sp["med_lo"]])
        m_pos, _ = group_stats(vec[sp["pos"]], mal_lib[sp["pos"]]) if np.any(sp["pos"]) else (np.nan, np.nan)
        m_neg, _ = group_stats(vec[sp["neg"]], mal_lib[sp["neg"]]) if np.any(sp["neg"]) else (np.nan, np.nan)
        mean_mal, pct_mal = group_stats(vec, mal_lib)
        rows.append(
            {
                "cohort": cohort,
                "unit_id": unit,
                "ligand": lig,
                "delta_q4q1": m_q4 - m_q1,
                "mean_q4": m_q4,
                "mean_q1": m_q1,
                "delta_median": m_mh - m_ml,
                "delta_pctpos": m_pos - m_neg,
                "pct_mal": pct_mal,
                "mean_mal": mean_mal,
            }
        )
    return rows


def marker_masks(symbol_to_row: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def any_pos(genes):
        hit = False
        acc = np.zeros(n, dtype=bool)
        for g in genes:
            v = symbol_to_row.get(g)
            if v is None:
                continue
            acc |= v > 0
            hit = True
        return acc if hit else np.zeros(n, dtype=bool)

    ptprc = symbol_to_row.get("PTPRC", np.zeros(n, dtype=np.float32))
    mal = any_pos(EPI) & (ptprc == 0)
    tnk = any_pos(TNK) & ~mal
    return mal, tnk


def finalize_unit(cohort, unit, n_cells, n_mal, n_tnk, sp) -> dict:
    eligible = n_mal >= MIN_MAL and n_tnk >= MIN_TNK and sp["q_ok"]
    separated = bool(
        eligible
        and np.isfinite(sp["mean_hi"])
        and np.isfinite(sp["mean_lo"])
        and sp["mean_hi"] > sp["mean_lo"]
    )
    return {
        "cohort": cohort,
        "unit_id": unit,
        "n_cells": int(n_cells),
        "n_mal": int(n_mal),
        "n_tnk": int(n_tnk),
        "n_high": int(np.sum(sp["q_hi"])),
        "n_low": int(np.sum(sp["q_lo"])),
        "cldn4_pct": sp["pct"],
        "cldn4_mean": sp["mean"],
        "cldn4_mean_high": sp["mean_hi"],
        "cldn4_mean_low": sp["mean_lo"],
        "eligible_q4": eligible,
        "cldn4_separated": separated,
    }


def dense_patient(cohort, unit, symbols: list[str], X: np.ndarray, ligands: list[str]):
    """X is genes x cells float32."""
    n_genes, n_cells = X.shape
    lib = X.sum(axis=0).astype(np.float64)
    sym = [canon(s) for s in symbols]
    # first occurrence wins
    index = {}
    for i, g in enumerate(sym):
        if g not in index:
            index[g] = i
    rows = {g: X[i] for g, i in index.items()}
    mal, tnk = marker_masks(rows, n_cells)
    mal_ix = np.flatnonzero(mal)
    tnk_ix = np.flatnonzero(tnk)
    cldn4 = rows.get("CLDN4", np.zeros(n_cells, dtype=np.float32))[mal_ix]
    mal_lib = lib[mal_ix]
    sp = splits_from_mal(cldn4.astype(np.float64), mal_lib)
    unit_row = finalize_unit(cohort, unit, n_cells, mal_ix.size, tnk_ix.size, sp)
    tnk_lib = lib[tnk_ix]
    # T/NK means for every gene we have (later intersected with the prior)
    genes = list(index.keys())
    mean = np.empty(len(genes), dtype=np.float32)
    pct = np.empty(len(genes), dtype=np.float32)
    if tnk_ix.size:
        sub = X[:, tnk_ix]
        logged = np.log1p(1e4 * sub / np.maximum(tnk_lib, 1.0))
        # duplicate symbols: first row only, matching index
        for j, g in enumerate(genes):
            i = index[g]
            mean[j] = logged[i].mean()
            pct[j] = np.mean(sub[i] > 0)
    else:
        mean[:] = np.nan
        pct[:] = np.nan
    counts_by_gene = {}
    if mal_ix.size:
        for lig in ligands:
            i = index.get(lig)
            if i is None:
                continue
            counts_by_gene[lig] = X[i, mal_ix]
    lig_rows = []
    if unit_row["eligible_q4"]:
        lig_rows = ligand_rows(cohort, unit, ligands, index, counts_by_gene, mal_lib, sp)
    return unit_row, genes, mean, pct, lig_rows


def extract_gse123902(ligands: list[str]):
    meta = pd.read_csv(os.path.join(HERE, "data", "GSE123902_marker_units.tsv"), sep="\t")
    tumor = meta[meta.tissue.isin(["PRIMARY", "METASTASIS"])].copy()
    tumor["pri"] = (tumor.tissue != "PRIMARY").astype(int)
    tumor = tumor.sort_values(["patient", "pri"]).drop_duplicates("patient")
    csv_dir = os.path.join(RAW, "GSE123902", "csv")
    os.makedirs(csv_dir, exist_ok=True)
    if not any(n.endswith("_dense.csv.gz") for n in os.listdir(csv_dir)):
        log("untar GSE123902")
        with tarfile.open(os.path.join(RAW, "GSE123902", "GSE123902_RAW.tar")) as tar:
            tar.extractall(csv_dir)
    units, ligs = [], []
    gene_store = {}
    for rec in tumor.itertuples(index=False):
        fp = os.path.join(csv_dir, rec.file)
        if not os.path.exists(fp):
            log(f"missing {rec.file}")
            continue
        log(f"GSE123902 {rec.patient} {rec.tissue}")
        df = pd.read_csv(fp)
        symbols = [str(c) for c in df.columns[1:]]
        X = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False).T
        unit_row, genes, mean, pct, lig_rows = dense_patient(
            "GSE123902", str(rec.patient), symbols, X, ligands
        )
        units.append(unit_row)
        ligs.extend(lig_rows)
        gene_store[str(rec.patient)] = (genes, mean, pct)
        log(
            f"  cells {unit_row['n_cells']} mal {unit_row['n_mal']} tnk {unit_row['n_tnk']} "
            f"eligible {unit_row['eligible_q4']} sep {unit_row['cldn4_separated']}"
        )
        del df, X
    _dump("GSE123902", units, ligs, gene_store)


def _dump(cohort, units, ligs, gene_store):
    pd.DataFrame(units).to_csv(os.path.join(OUT, f"units_{cohort}.tsv"), sep="\t", index=False)
    pd.DataFrame(ligs).to_csv(os.path.join(OUT, f"ligands_{cohort}.tsv"), sep="\t", index=False)
    # align genes across units: union, missing = nan
    all_genes = sorted(set().union(*[set(g) for g, _, _ in gene_store.values()])) if gene_store else []
    gix = {g: i for i, g in enumerate(all_genes)}
    cols = list(gene_store.keys())
    mean = np.full((len(all_genes), len(cols)), np.nan, dtype=np.float32)
    pct = np.full_like(mean, np.nan)
    for j, u in enumerate(cols):
        genes, m, p = gene_store[u]
        for g, mv, pv in zip(genes, m, p):
            mean[gix[g], j] = mv
            pct[gix[g], j] = pv
    np.savez_compressed(
        os.path.join(OUT, f"tnk_{cohort}.npz"),
        genes=np.array(all_genes, dtype=object),
        units=np.array(cols, dtype=object),
        mean=mean,
        pct=pct,
    )
    log(f"wrote {cohort} units {len(units)} ligand-rows {len(ligs)} tnk genes {len(all_genes)}")


def extract_gse189357(ligands: list[str]):
    meta = pd.read_csv(os.path.join(HERE, "data", "GSE189357_sample_metadata.tsv"), sep="\t")
    ex = os.path.join(RAW, "GSE189357", "raw")
    os.makedirs(ex, exist_ok=True)
    need = f"{meta.gsm.iloc[0]}_{meta.patient.iloc[0]}_matrix.mtx.gz"
    if not os.path.exists(os.path.join(ex, need)):
        log("untar GSE189357")
        with tarfile.open(os.path.join(RAW, "GSE189357", "GSE189357_RAW.tar")) as tar:
            tar.extractall(ex)
    units, ligs = [], []
    gene_store = {}
    for rec in meta.itertuples(index=False):
        prefix = os.path.join(ex, f"{rec.gsm}_{rec.patient}")
        mtx = prefix + "_matrix.mtx.gz"
        feat = prefix + "_features.tsv.gz"
        log(f"GSE189357 {rec.patient}")
        symbols = pd.read_csv(feat, sep="\t", header=None)[1].astype(str).tolist()
        sp = mmread(mtx).tocsr()
        # 10x MTX is genes x cells
        if sp.shape[0] != len(symbols):
            raise SystemExit(f"feature/matrix mismatch {rec.patient} {sp.shape} vs {len(symbols)}")
        X = sp.astype(np.float32)
        # densify only if reasonable
        n_cells = X.shape[1]
        if n_cells > 40000:
            raise SystemExit("unexpected cell count")
        Xd = X.toarray()
        unit_row, genes, mean, pct, lig_rows = dense_patient(
            "GSE189357", str(rec.patient), symbols, Xd, ligands
        )
        units.append(unit_row)
        ligs.extend(lig_rows)
        gene_store[str(rec.patient)] = (genes, mean, pct)
        log(
            f"  cells {unit_row['n_cells']} mal {unit_row['n_mal']} tnk {unit_row['n_tnk']} "
            f"eligible {unit_row['eligible_q4']}"
        )
        del X, Xd, sp
    _dump("GSE189357", units, ligs, gene_store)


def _fromstring_line(rest: str, n: int) -> np.ndarray:
    v = np.fromstring(rest, sep="\t", dtype=np.float32)
    if v.size == n:
        return v
    if v.size > n:
        return v[:n]
    out = np.zeros(n, dtype=np.float32)
    out[: v.size] = v
    return out


def extract_gse131907(ligands: list[str]):
    samples = pd.read_csv(os.path.join(HERE, "data", "GSE131907_samples.tsv"), sep="\t")
    locked = samples.loc[samples.n_malignant > 0, "sample"].astype(str).tolist()
    locked_set = set(locked)
    log(f"GSE131907 locked samples {len(locked)}")
    ann = pd.read_csv(
        os.path.join(RAW, "GSE131907", "GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
        sep="\t",
    )
    ann = ann[ann["Sample"].astype(str).isin(locked_set)]
    sid = {s: i for i, s in enumerate(locked)}
    path = os.path.join(RAW, "GSE131907", "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")[1:]
    n = len(header)
    log(f"  matrix cells {n} annotated locked {len(ann)}")
    code = np.full(n, -1, dtype=np.int16)
    is_mal = np.zeros(n, dtype=bool)
    is_tnk = np.zeros(n, dtype=bool)
    ann_map = {}
    for rec in ann.itertuples(index=False):
        ann_map[str(rec.Index)] = (
            sid[str(rec.Sample)],
            str(rec.Cell_subtype) in MALIG_SUB,
            str(rec.Cell_type) in TNK_TYPES,
        )
    for i, cell in enumerate(header):
        hit = ann_map.get(cell)
        if hit is None:
            continue
        code[i] = hit[0]
        is_mal[i] = hit[1]
        is_tnk[i] = hit[2] and not hit[1]
    log(f"  matched mal {int(is_mal.sum())} tnk {int(is_tnk.sum())}")

    mal_idx = [np.flatnonzero((code == i) & is_mal) for i in range(len(locked))]
    tnk_idx = [np.flatnonzero((code == i) & is_tnk) for i in range(len(locked))]

    log("  pass 1 library size + CLDN4")
    lib = np.zeros(n, dtype=np.float64)
    cldn4 = np.zeros(n, dtype=np.float32)
    t0 = time.time()
    with gzip.open(path, "rt") as f:
        f.readline()
        k = 0
        for line in f:
            tab = line.find("\t")
            gene = canon(line[:tab])
            v = _fromstring_line(line[tab + 1 :], n)
            lib += v
            if gene == "CLDN4":
                cldn4 = v
            k += 1
            if k % 4000 == 0:
                log(f"    pass1 {k} genes {time.time()-t0:.0f}s")
    log(f"  pass1 done {k} genes {time.time()-t0:.0f}s lib median {np.median(lib):.0f}")

    splits = []
    units = []
    for i, sample in enumerate(locked):
        ix = mal_idx[i]
        c_umi = cldn4[ix].astype(np.float64)
        c_lib = lib[ix]
        sp = splits_from_mal(c_umi, c_lib)
        # masks are local to malignant cells; also keep global index arrays for arms
        sp["g_q_hi"] = ix[sp["q_hi"]] if ix.size else ix
        sp["g_q_lo"] = ix[sp["q_lo"]] if ix.size else ix
        sp["g_med_hi"] = ix[sp["med_hi"]] if ix.size else ix
        sp["g_med_lo"] = ix[sp["med_lo"]] if ix.size else ix
        sp["g_pos"] = ix[sp["pos"]] if ix.size else ix
        sp["g_neg"] = ix[sp["neg"]] if ix.size else ix
        sp["g_all"] = ix
        splits.append(sp)
        units.append(
            finalize_unit("GSE131907", sample, int((code == i).sum()), ix.size, tnk_idx[i].size, sp)
        )
        log(
            f"  {sample} mal {ix.size} tnk {tnk_idx[i].size} "
            f"pct {sp['pct']:.3f} eligible {units[-1]['eligible_q4']} sep {units[-1]['cldn4_separated']}"
        )

    # pass 2: T/NK profiles for prior targets, ligand deltas for eligible units
    target_set = set()
    with open(os.path.join(PRIOR, "targets.txt")) as f:
        for line in f:
            target_set.add(canon(line.strip()))
    ligand_set = set(ligands)
    want = target_set | ligand_set | {"CLDN4"}
    log(f"  pass 2 wanted genes {len(want)}")
    # store tnk mean/pct only for genes we see
    tnk_mean = defaultdict(lambda: np.full(len(locked), np.nan, dtype=np.float32))
    tnk_pct = defaultdict(lambda: np.full(len(locked), np.nan, dtype=np.float32))
    # ligand malignant counts are computed on the fly into rows
    lig_rows = []
    # prebuild per-ligand accumulators: we'll append when gene is a ligand
    t1 = time.time()
    seen = 0
    with gzip.open(path, "rt") as f:
        f.readline()
        k = 0
        for line in f:
            tab = line.find("\t")
            gene = canon(line[:tab])
            k += 1
            if gene not in want:
                if k % 4000 == 0:
                    log(f"    pass2 scanned {k} kept {seen} {time.time()-t1:.0f}s")
                continue
            v = _fromstring_line(line[tab + 1 :], n)
            seen += 1
            for i in range(len(locked)):
                ix = tnk_idx[i]
                if ix.size == 0:
                    continue
                m, p = group_stats(v[ix], lib[ix])
                tnk_mean[gene][i] = m
                tnk_pct[gene][i] = p
            if gene in ligand_set:
                for i, sample in enumerate(locked):
                    if not units[i]["eligible_q4"]:
                        continue
                    sp = splits[i]
                    m_q4, _ = group_stats(v[sp["g_q_hi"]], lib[sp["g_q_hi"]])
                    m_q1, _ = group_stats(v[sp["g_q_lo"]], lib[sp["g_q_lo"]])
                    m_mh, _ = group_stats(v[sp["g_med_hi"]], lib[sp["g_med_hi"]])
                    m_ml, _ = group_stats(v[sp["g_med_lo"]], lib[sp["g_med_lo"]])
                    m_pos, _ = (
                        group_stats(v[sp["g_pos"]], lib[sp["g_pos"]])
                        if sp["g_pos"].size
                        else (np.nan, np.nan)
                    )
                    m_neg, _ = (
                        group_stats(v[sp["g_neg"]], lib[sp["g_neg"]])
                        if sp["g_neg"].size
                        else (np.nan, np.nan)
                    )
                    mean_mal, pct_mal = group_stats(v[sp["g_all"]], lib[sp["g_all"]])
                    lig_rows.append(
                        {
                            "cohort": "GSE131907",
                            "unit_id": sample,
                            "ligand": gene,
                            "delta_q4q1": m_q4 - m_q1,
                            "mean_q4": m_q4,
                            "mean_q1": m_q1,
                            "delta_median": m_mh - m_ml,
                            "delta_pctpos": m_pos - m_neg,
                            "pct_mal": pct_mal,
                            "mean_mal": mean_mal,
                        }
                    )
            if seen % 2000 == 0:
                log(f"    pass2 kept {seen} scanned {k} {time.time()-t1:.0f}s")
    genes = sorted(tnk_mean.keys())
    mean = np.vstack([tnk_mean[g] for g in genes]) if genes else np.zeros((0, len(locked)))
    pct = np.vstack([tnk_pct[g] for g in genes]) if genes else mean.copy()
    pd.DataFrame(units).to_csv(os.path.join(OUT, "units_GSE131907.tsv"), sep="\t", index=False)
    pd.DataFrame(lig_rows).to_csv(os.path.join(OUT, "ligands_GSE131907.tsv"), sep="\t", index=False)
    np.savez_compressed(
        os.path.join(OUT, "tnk_GSE131907.npz"),
        genes=np.array(genes, dtype=object),
        units=np.array(locked, dtype=object),
        mean=mean.astype(np.float32),
        pct=pct.astype(np.float32),
    )
    log(f"wrote GSE131907 units {len(units)} ligand-rows {len(lig_rows)} tnk genes {len(genes)}")


def main():
    ligands = load_ligands()
    which = os.environ.get("EXTRACT_WHICH", "123902,189357,131907")
    if "123902" in which:
        extract_gse123902(ligands)
    if "189357" in which:
        extract_gse189357(ligands)
    if "131907" in which:
        extract_gse131907(ligands)


if __name__ == "__main__":
    main()
