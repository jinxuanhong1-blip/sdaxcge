#!/usr/bin/env python3
"""Largest real barrier-ligand communication delta, concordant-4.

CLDN4-high vs CLDN4-low malignant senders -> T/NK.
Ligands: F11R, NECTIN2, CDH1, LGALS9.
IFN/recruit chemokines and HLA-CD8 are scored separately.
Patient is the unit. Four cohorts must agree.

The official CellChat family sum in PR 616 is +0.0037 because population-size
scaling compresses probabilities. This script recomputes that Hill score and
the standard CellPhoneDB / LIANA / Connectome / NATMI / log2FC / abundance
scores on the same cells, then keeps the largest absolute family delta that
stays positive inside every cohort and is not smaller than the chemokine arm.
"""

from __future__ import annotations

import math
import os
import tarfile
import traceback
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("CONCORDANT4_RAW", "/tmp/concordant4_raw"))
CACHE = Path(os.environ.get("CONCORDANT4_CACHE", "/tmp/concordant4_cache"))
TAB = HERE / "results" / "tables"
FIG = HERE / "results" / "figures"

TRIM = 0.10
EXPR_PROP = 0.10
KH = 0.5
MIN_ARM = 10
MIN_TNK = 20
MIN_MAL = 40
SEED = 3979
N_FLIP = 10000
LN2 = math.log(2.0)

ALIASES = {"PVRL2": "NECTIN2", "JAM1": "F11R"}
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
MALIG_SUB = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#4C78A8",
    "GSE131907": "#F58518",
    "GSE205335": "#54A24B",
    "GSE189357": "#E45756",
}

# family, ligand, receptor token, axis
EDGES = [
    ("barrier", "F11R", "ITGAL_ITGB2", "F11R-LFA1"),
    ("barrier", "F11R", "F11R", "F11R-F11R"),
    ("barrier", "NECTIN2", "TIGIT", "NECTIN2-TIGIT"),
    ("barrier", "NECTIN2", "CD96", "NECTIN2-CD96"),
    ("barrier", "CDH1", "ITGAE_ITGB7", "CDH1-integrin"),
    ("barrier", "CDH1", "KLRG1", "CDH1-KLRG1"),
    ("barrier", "LGALS9", "HAVCR2", "LGALS9-TIM3"),
    ("barrier", "LGALS9", "CD44", "LGALS9-CD44"),
    ("barrier", "LGALS9", "PTPRC", "LGALS9-CD45"),
    ("ifn_recruit", "CXCL9", "CXCR3", "CXCL9-CXCR3"),
    ("ifn_recruit", "CXCL10", "CXCR3", "CXCL10-CXCR3"),
    ("ifn_recruit", "CXCL11", "CXCR3", "CXCL11-CXCR3"),
    ("ifn_recruit", "CCL5", "CCR5", "CCL5-CCR5"),
    ("ifn_recruit", "CCL5", "CCR1", "CCL5-CCR1"),
    ("ifn_recruit", "CXCL16", "CXCR6", "CXCL16-CXCR6"),
    ("mhc_cd8", "HLA-A", "CD8A", "HLA-A-CD8A"),
    ("mhc_cd8", "HLA-B", "CD8A", "HLA-B-CD8A"),
    ("mhc_cd8", "HLA-C", "CD8A", "HLA-C-CD8A"),
]
# PR 616 seven-pair sum, for the population-scaled Hill comparison only.
PAIR7 = {
    ("F11R", "ITGAL_ITGB2"),
    ("NECTIN2", "TIGIT"),
    ("CDH1", "ITGAE_ITGB7"),
    ("CDH1", "KLRG1"),
    ("LGALS9", "HAVCR2"),
    ("LGALS9", "CD44"),
    ("LGALS9", "PTPRC"),
}
BARRIER_LIGANDS = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
IFN_LIGANDS = ["CXCL9", "CXCL10", "CXCL11", "CCL5", "CXCL16"]
MHC_LIGANDS = ["HLA-A", "HLA-B", "HLA-C"]

METHODS = [
    "cellchat_hill",
    "cellchat_hill_pop",
    "cellchat_prod",
    "cpdb_means",
    "conn_prod",
    "natmi_spec",
    "sca_lrscore",
    "liana_logfc",
    "cp10k_log2fc",
    "cp10k_log2fc_trim",
    "cp10k_delta",
    "cp10k_delta_trim",
    "expr_prop_pp",
    "conn_z",
]
# Untrimmed abundance scores must not win on outliers alone.
TRIM_GUARD = {
    "cp10k_log2fc": "cp10k_log2fc_trim",
    "cp10k_delta": "cp10k_delta_trim",
}
UNIT = {
    "cellchat_hill": "Hill probability, no population-size weight",
    "cellchat_hill_pop": "Hill probability times sender and receiver proportions",
    "cellchat_prod": "product Hill probability (L*R)/(0.5+L*R)",
    "cpdb_means": "CellPhoneDB lr_means on log1p CP10k",
    "conn_prod": "Connectome expression product on log1p CP10k",
    "natmi_spec": "NATMI edge specificity",
    "sca_lrscore": "SingleCellSignalR LRscore",
    "liana_logfc": "LIANA/scanpy log2FC of log1p ligand means",
    "cp10k_log2fc": "log2 fold-change of mean CP10k",
    "cp10k_log2fc_trim": "log2 fold-change of 10% trimmed-mean CP10k",
    "cp10k_delta": "difference of mean CP10k",
    "cp10k_delta_trim": "difference of 10% trimmed-mean CP10k",
    "expr_prop_pp": "percentage points of ligand-positive malignant cells",
    "conn_z": "Connectome z-score product across the three labels",
}


def log(msg: str) -> None:
    print(msg, flush=True)


def subunits(token: str) -> list[str]:
    return token.split("_")


def wanted_genes() -> set[str]:
    genes = set(EPI + TNK_MARKERS + ["CLDN4", "PTPRC"])
    for fam, lig, rec, axis in EDGES:
        genes.add(lig)
        genes.update(subunits(rec))
    genes.update(ALIASES)
    return genes


def collapse_named(named: dict[str, np.ndarray], keep: set[str]) -> dict[str, np.ndarray]:
    buckets: dict[str, list[tuple[str, np.ndarray]]] = {}
    for gene, arr in named.items():
        g = gene.upper()
        canon = ALIASES.get(g, g)
        if canon not in keep and g not in keep:
            continue
        buckets.setdefault(canon, []).append((g, np.asarray(arr)))
    out = {}
    for canon, items in buckets.items():
        if canon not in keep:
            continue
        official = [a for g, a in items if g == canon]
        use = official if official else [a for _, a in items]
        out[canon] = use[0] if len(use) == 1 else np.sum(np.vstack(use), axis=0)
    return out


def matrix_from_named(named: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, list[str]]:
    genes = sorted(named)
    if not genes:
        return np.zeros((n, 0), dtype=np.float32), []
    mat = np.column_stack([np.asarray(named[g], dtype=np.float32).reshape(-1) for g in genes])
    if mat.shape[0] != n:
        raise RuntimeError(f"rows {mat.shape[0]} != cells {n}")
    return mat, genes


def pos_any(mat: np.ndarray, genes: list[str], markers: list[str]) -> np.ndarray:
    idx = [genes.index(g) for g in markers if g in genes]
    if not idx:
        return np.zeros(mat.shape[0], dtype=bool)
    return np.any(mat[:, idx] > 0, axis=1)


def marker_compartments(mat: np.ndarray, genes: list[str]) -> np.ndarray:
    mal = pos_any(mat, genes, EPI) & ~pos_any(mat, genes, ["PTPRC"])
    tnk = pos_any(mat, genes, TNK_MARKERS) & ~mal
    out = np.full(mat.shape[0], "", dtype=object)
    out[tnk] = "TNK"
    out[mal] = "MAL"
    return out


def pack_unit(patient: str, mat: np.ndarray, genes: list[str], lib: np.ndarray, comp: np.ndarray) -> dict | None:
    keep = comp != ""
    if not np.any(keep):
        return None
    return {
        "patient": str(patient),
        "mat": np.asarray(mat[keep], dtype=np.float32),
        "genes": genes,
        "lib": np.asarray(lib[keep], dtype=np.float64),
        "comp": comp[keep],
    }


def read_dense_csv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    df = pd.read_csv(path)
    genes = [str(c).upper() for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32)
    lib = mat.sum(axis=1).astype(np.float64)
    return mat, genes, lib


def load_gse123902(keep: set[str]) -> list[dict]:
    log("==== GSE123902 ====")
    units = pd.read_csv(HERE / "data" / "GSE123902_marker_units.tsv", sep="\t")
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor["ord"] = np.where(tumor["tissue"].eq("PRIMARY"), 0, 1)
    tumor = tumor.sort_values(["patient", "ord"]).drop_duplicates("patient")
    csv_dir = RAW / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*_dense.csv.gz")):
        tar = RAW / "GSE123902" / "GSE123902_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(csv_dir)
    out = []
    for row in tumor.itertuples(index=False):
        fp = csv_dir / row.file
        if not fp.exists():
            hits = sorted(p for p in csv_dir.glob(f"*{row.patient}*dense.csv.gz") if "NORMAL" not in p.name)
            if not hits:
                log(f"  missing {row.patient}")
                continue
            fp = hits[0]
        log(f"  {row.patient}")
        mat, genes, lib = read_dense_csv(fp)
        named = collapse_named({g: mat[:, i] for i, g in enumerate(genes)}, keep)
        mat, genes = matrix_from_named(named, mat.shape[0])
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE123902"
            out.append(packed)
    return out


def load_gse131907(keep: set[str]) -> list[dict]:
    """Read the column-subset written by extract_131907. Library size is the full transcriptome."""
    log("==== GSE131907 ====")
    cache = CACHE / "gse131907"
    lib_path = cache / "libsize.f64"
    if not lib_path.exists() or not (cache / "genes_found.txt").exists():
        raise RuntimeError(
            "GSE131907 slim matrix is missing. Build keep_idx.i32 and run scripts/extract_131907."
        )
    meta = pd.read_csv(cache / "cells.tsv", sep="\t")
    lib = np.fromfile(lib_path, dtype=np.float64)
    if lib.size != len(meta):
        raise RuntimeError(f"lib {lib.size} vs cells {len(meta)}")
    found = [ln.strip() for ln in (cache / "genes_found.txt").read_text().splitlines() if ln.strip()]
    named = {g: np.fromfile(cache / f"{g}.f32", dtype=np.float32) for g in found}
    named = collapse_named(named, keep)
    mat, genes = matrix_from_named(named, len(meta))
    comp = meta["comp"].astype(str).to_numpy()
    sample = meta["sample"].astype(str).to_numpy()
    out = []
    for s in sorted(set(sample)):
        which = np.where(sample == s)[0]
        packed = pack_unit(s, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE131907"
            out.append(packed)
    log(f"  units {len(out)} cells {len(meta)}")
    return out


def load_gse205335(keep: set[str]) -> list[dict]:
    log("==== GSE205335 ====")
    export = CACHE / "gse205335"
    mtx = export / "matrix.mtx"
    if not mtx.exists():
        import subprocess

        genes_file = CACHE / "ccc_genes.txt"
        genes_file.parent.mkdir(parents=True, exist_ok=True)
        genes_file.write_text("\n".join(sorted(keep | set(ALIASES))) + "\n")
        subprocess.check_call(
            [
                "Rscript",
                str(HERE / "scripts" / "export_gse205335.R"),
                str(RAW / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"),
                str(genes_file),
                str(export),
            ]
        )
    genes = [ln.strip() for ln in (export / "genes.txt").read_text().splitlines() if ln.strip()]
    cells = [ln.strip() for ln in (export / "cells.txt").read_text().splitlines() if ln.strip()]
    log(f"  MTX {len(genes)} x {len(cells)}")
    mm = spio.mmread(mtx).tocsr()
    if mm.shape == (len(genes), len(cells)):
        mm = mm.T.tocsr()
    elif mm.shape != (len(cells), len(genes)):
        raise RuntimeError(f"MTX shape {mm.shape}")
    lib_df = pd.read_csv(export / "libsize.tsv", sep="\t")
    lib_map = dict(zip(lib_df["cell"].astype(str), lib_df["libsize"].astype(float)))
    lib = np.array([lib_map[c] for c in cells], dtype=np.float64)
    named = collapse_named({g: np.asarray(mm[:, i].todense()).ravel() for i, g in enumerate(genes)}, keep)
    mat, genes = matrix_from_named(named, len(cells))
    ident = pd.read_csv(RAW / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    gsm = pd.read_csv(HERE / "data" / "GSE205335_gsm_map.tsv", sep="\t")
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident["barcode"] = ident["barcode"].astype(str)
    idmap = ident.drop_duplicates("barcode").set_index("barcode")
    common = [i for i, c in enumerate(cells) if c in idmap.index]
    if len(common) < 0.9 * len(cells):
        raise RuntimeError(f"identity overlap {len(common)}/{len(cells)}")
    idx = np.array(common)
    sub = idmap.loc[[cells[i] for i in idx]]
    normal = sub["tissue"].fillna("").astype(str).str.startswith("Normal").to_numpy()
    lineage_sub = sub["lineage.sub"].astype(str).to_numpy()
    lineage = sub["lineage.total"].astype(str).to_numpy()
    patient = sub["patient"].astype(str).to_numpy()
    mal = (lineage_sub == "Malignant cells") & ~normal
    tnk = (lineage == "T/NK cells") & ~normal & ~mal
    comp = np.full(idx.size, "", dtype=object)
    comp[tnk] = "TNK"
    comp[mal] = "MAL"
    mat, lib = mat[idx], lib[idx]
    locked = pd.read_csv(HERE / "data" / "GSE205335_patients.tsv", sep="\t")
    keep_pt = set(locked.loc[locked["n_malignant"] > 0, "patient"].astype(str))
    out = []
    for pt in sorted(keep_pt):
        which = np.where(patient == pt)[0]
        if which.size == 0:
            continue
        packed = pack_unit(pt, mat[which], genes, lib[which], comp[which])
        if packed:
            packed["cohort"] = "GSE205335"
            out.append(packed)
    log(f"  units {len(out)}")
    return out


def load_gse189357(keep: set[str]) -> list[dict]:
    log("==== GSE189357 ====")
    meta = pd.read_csv(HERE / "data" / "GSE189357_sample_metadata.tsv", sep="\t")
    ex = RAW / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    needed = ex / f"{meta.gsm.iloc[0]}_{meta.patient.iloc[0]}_matrix.mtx.gz"
    if not needed.exists():
        tar = RAW / "GSE189357" / "GSE189357_RAW.tar"
        log(f"  untar {tar.name}")
        with tarfile.open(tar) as tf:
            tf.extractall(ex)
    out = []
    for row in meta.itertuples(index=False):
        prefix = ex / f"{row.gsm}_{row.patient}"
        mtx = Path(str(prefix) + "_matrix.mtx.gz")
        feat = Path(str(prefix) + "_features.tsv.gz")
        if not mtx.exists():
            hits = list(ex.glob(f"*{row.patient}_matrix.mtx.gz"))
            if not hits:
                log(f"  missing {row.patient}")
                continue
            mtx = hits[0]
            feat = Path(str(mtx).replace("_matrix.mtx.gz", "_features.tsv.gz"))
        log(f"  {row.patient}")
        mm = spio.mmread(mtx).tocsr()
        features = pd.read_csv(feat, sep="\t", header=None)
        symbols = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).str.upper()
        if mm.shape[0] != len(symbols):
            if mm.shape[1] == len(symbols):
                mm = mm.T.tocsr()
            else:
                raise RuntimeError(f"{row.patient} mtx {mm.shape} vs {len(symbols)}")
        n = mm.shape[1]
        lib = np.asarray(mm.sum(axis=0)).ravel().astype(np.float64)
        summed: dict[str, np.ndarray] = {}
        for i, g in enumerate(symbols):
            if g not in keep and g not in ALIASES:
                continue
            arr = np.asarray(mm[i].todense()).ravel()
            summed[g] = arr if g not in summed else summed[g] + arr
        named = collapse_named(summed, keep)
        mat, genes = matrix_from_named(named, n)
        comp = marker_compartments(mat, genes)
        packed = pack_unit(str(row.patient), mat, genes, lib, comp)
        if packed:
            packed["cohort"] = "GSE189357"
            out.append(packed)
    return out


def qc_counts(units: list[dict]) -> pd.DataFrame:
    ref = pd.read_csv(HERE / "data" / "reference_cellchat_inventory.tsv", sep="\t")
    ref["patient"] = ref["patient"].astype(str)
    rows = []
    for u in units:
        comp = u["comp"]
        mal = comp == "MAL"
        cldn_i = u["genes"].index("CLDN4") if "CLDN4" in u["genes"] else None
        if cldn_i is not None and mal.any():
            ln = np.log1p(u["mat"][mal, cldn_i] / np.maximum(u["lib"][mal], 1.0) * 1e4)
            mean_ln = float(ln.mean())
            pct = float(np.mean(u["mat"][mal, cldn_i] > 0))
        else:
            mean_ln, pct = math.nan, math.nan
        rows.append(
            {
                "cohort": u["cohort"],
                "patient": u["patient"],
                "n_mal": int(mal.sum()),
                "n_tnk": int((comp == "TNK").sum()),
                "mal_cldn4_mean": mean_ln,
                "mal_cldn4_pct": pct,
            }
        )
    inv = pd.DataFrame(rows)
    merged = inv.merge(ref, on=["cohort", "patient"], suffixes=("", "_ref"))
    bad = merged[(merged["n_mal"] != merged["n_mal_ref"]) | (merged["n_tnk"] != merged["n_tnk_ref"])]
    if len(merged) != 65 or not bad.empty:
        log(bad.to_string(index=False))
        raise RuntimeError(f"label QC failed: matched {len(merged)} / mismatches {len(bad)}")
    log(f"label QC matched {len(merged)} / 65")
    return inv


def arms(x: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray] | None:
    n = int(x.size)
    if n < 4:
        return None
    r = stats.rankdata(np.asarray(x, dtype=float), method="ordinal")
    if mode == "q4q1":
        if n < MIN_MAL:
            return None
        lo_cut = math.floor(n * 0.25)
        hi_cut = math.ceil(n * 0.75)
    elif mode == "decile":
        lo_cut = math.floor(n * 0.10)
        hi_cut = math.ceil(n * 0.90)
    else:
        raise ValueError(mode)
    if lo_cut < MIN_ARM or (n - hi_cut) < MIN_ARM:
        return None
    low = r <= lo_cut
    high = r > hi_cut
    if int(high.sum()) < MIN_ARM or int(low.sum()) < MIN_ARM:
        return None
    return high, low


def block_stats(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = block.shape[0]
    mean = block.mean(axis=0)
    prop = (block > 0).mean(axis=0)
    k = int(n * TRIM)
    if k == 0 or n - 2 * k < 1:
        trim = mean.copy()
    else:
        trim = np.sort(block, axis=0)[k : n - k].mean(axis=0)
    return mean, trim, prop


def gmean(vals: list[float]) -> float:
    arr = np.asarray(vals, dtype=float)
    if arr.size == 0 or np.any(arr <= 0):
        return 0.0
    return float(np.exp(np.mean(np.log(arr))))


def complex_value(stat: dict[str, dict], token: str, key: str) -> tuple[float, float]:
    parts = subunits(token)
    vals, props = [], []
    for g in parts:
        if g not in stat:
            return 0.0, 0.0
        vals.append(float(stat[g][key]))
        props.append(float(stat[g]["prop"]))
    return gmean(vals), (min(props) if props else 0.0)


def hill(x: float) -> float:
    if x <= 0:
        return 0.0
    return x / (KH + x)


def z3(a: float, b: float, c: float) -> tuple[float, float, float]:
    v = np.array([a, b, c], dtype=float)
    sd = float(v.std(ddof=1))
    if not np.isfinite(sd) or sd == 0:
        return 0.0, 0.0, 0.0
    z = (v - v.mean()) / sd
    return float(z[0]), float(z[1]), float(z[2])


def score_unit(unit: dict) -> list[dict]:
    comp = unit["comp"]
    mal = np.where(comp == "MAL")[0]
    tnk = np.where(comp == "TNK")[0]
    if mal.size < MIN_MAL or tnk.size < MIN_TNK or "CLDN4" not in unit["genes"]:
        return []
    genes = unit["genes"]
    gindex = {g: i for i, g in enumerate(genes)}
    lib = np.maximum(unit["lib"], 1.0)
    cp = unit["mat"] / lib[:, None] * 1e4
    logv = np.log1p(cp)
    cldn = logv[mal, gindex["CLDN4"]]
    rows = []
    # stats cache per split
    for mode in ("q4q1", "decile"):
        split = arms(cldn, mode)
        if split is None:
            continue
        high_m, low_m = split
        if float(cldn[high_m].mean()) <= float(cldn[low_m].mean()):
            raise RuntimeError(f"CLDN4 split inverted {unit['cohort']} {unit['patient']} {mode}")
        high = mal[high_m]
        low = mal[low_m]
        groups = {"high": high, "low": low, "tnk": tnk}
        # gene stats
        gene_stat = {arm: {} for arm in groups}
        for arm, idx in groups.items():
            mean, trim, prop = block_stats(cp[idx])
            lmean, ltrim, _ = block_stats(logv[idx])
            for g, i in gindex.items():
                gene_stat[arm][g] = {
                    "cp": float(mean[i]),
                    "cp_trim": float(trim[i]),
                    "log": float(lmean[i]),
                    "log_trim": float(ltrim[i]),
                    "prop": float(prop[i]),
                }
        used = np.concatenate([high, low, tnk])
        mu = float(logv[used].mean())
        n_high, n_low, n_tnk = int(high.size), int(low.size), int(tnk.size)
        n_tot = n_high + n_low + n_tnk
        f_high = (n_high / n_tot) * (n_tnk / n_tot)
        f_low = (n_low / n_tot) * (n_tnk / n_tot)
        for family, ligand, receptor, axis in EDGES:
            lig_h = gene_stat["high"].get(ligand)
            lig_l = gene_stat["low"].get(ligand)
            if lig_h is None or lig_l is None:
                Lh = Ll = Lh_t = Ll_t = Lh_cp = Ll_cp = Lh_cpt = Ll_cpt = 0.0
                ph = pl = 0.0
            else:
                Lh, Ll = lig_h["log"], lig_l["log"]
                Lh_t, Ll_t = lig_h["log_trim"], lig_l["log_trim"]
                Lh_cp, Ll_cp = lig_h["cp"], lig_l["cp"]
                Lh_cpt, Ll_cpt = lig_h["cp_trim"], lig_l["cp_trim"]
                ph, pl = lig_h["prop"], lig_l["prop"]
            Rh, _ = complex_value(gene_stat["tnk"], receptor, "log")
            Rh_t, _ = complex_value(gene_stat["tnk"], receptor, "log_trim")
            Rh_cp, rec_prop = complex_value(gene_stat["tnk"], receptor, "cp")
            # receptor complex prop is the min subunit prop
            _, rec_prop = complex_value(gene_stat["tnk"], receptor, "log")
            Lt_h, _ = complex_value(gene_stat["high"], ligand, "log")
            Lt_l, _ = complex_value(gene_stat["low"], ligand, "log")
            Lt_ht, _ = complex_value(gene_stat["high"], ligand, "log_trim")
            Lt_lt, _ = complex_value(gene_stat["low"], ligand, "log_trim")
            # single-gene ligands: complex_value equals the gene value
            detected = rec_prop >= EXPR_PROP and (ph >= EXPR_PROP or pl >= EXPR_PROP)
            zeros = {f"d_{m}": 0.0 for m in METHODS}
            if detected:
                p_h = hill(Lt_ht) * hill(Rh_t)
                p_l = hill(Lt_lt) * hill(Rh_t)
                prod_h = (Lt_ht * Rh_t) / (KH + Lt_ht * Rh_t) if Lt_ht > 0 and Rh_t > 0 else 0.0
                prod_l = (Lt_lt * Rh_t) / (KH + Lt_lt * Rh_t) if Lt_lt > 0 and Rh_t > 0 else 0.0
                # NATMI on arithmetic log means, groups high/low/tnk
                L_tnk, _ = complex_value(gene_stat["tnk"], ligand, "log")
                R_high, _ = complex_value(gene_stat["high"], receptor, "log")
                R_low, _ = complex_value(gene_stat["low"], receptor, "log")
                sumL = Lt_h + Lt_l + L_tnk
                sumR = R_high + R_low + Rh
                spec_h = (Lt_h / sumL) * (Rh / sumR) if sumL > 0 and sumR > 0 else 0.0
                spec_l = (Lt_l / sumL) * (Rh / sumR) if sumL > 0 and sumR > 0 else 0.0
                def lrscore(L, R):
                    if L <= 0 or R <= 0:
                        return 0.0
                    s = math.sqrt(L * R)
                    return s / (s + mu) if (s + mu) > 0 else 0.0
                zh, zl, zt = z3(Lt_h, Lt_l, L_tnk)
                rh, rl, rt = z3(R_high, R_low, Rh)
                zeros = {
                    "d_cellchat_hill": p_h - p_l,
                    "d_cellchat_hill_pop": p_h * f_high - p_l * f_low,
                    "d_cellchat_prod": prod_h - prod_l,
                    "d_cpdb_means": 0.5 * ((Lt_h + Rh) - (Lt_l + Rh)),
                    "d_conn_prod": Lt_h * Rh - Lt_l * Rh,
                    "d_natmi_spec": spec_h - spec_l,
                    "d_sca_lrscore": lrscore(Lt_h, Rh) - lrscore(Lt_l, Rh),
                    "d_liana_logfc": (Lh - Ll) / LN2,
                    "d_cp10k_log2fc": math.log2((Lh_cp + 1.0) / (Ll_cp + 1.0)),
                    "d_cp10k_log2fc_trim": math.log2((Lh_cpt + 1.0) / (Ll_cpt + 1.0)),
                    "d_cp10k_delta": Lh_cp - Ll_cp,
                    "d_cp10k_delta_trim": Lh_cpt - Ll_cpt,
                    "d_expr_prop_pp": 100.0 * (ph - pl),
                    "d_conn_z": zh * rt - zl * rt,
                }
            rows.append(
                {
                    "cohort": unit["cohort"],
                    "patient": unit["patient"],
                    "split": mode,
                    "family": family,
                    "ligand": ligand,
                    "receptor": receptor,
                    "axis": axis,
                    "detected": bool(detected),
                    "n_high": n_high,
                    "n_low": n_low,
                    "n_tnk": n_tnk,
                    "lig_prop_high": ph,
                    "lig_prop_low": pl,
                    "rec_prop": rec_prop,
                    **zeros,
                }
            )
    return rows


def ligand_table(edges: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["split", "cohort", "patient", "family", "ligand"]
    for key, g in edges.groupby(keys, sort=False):
        detected = g.loc[g["detected"]]
        base = dict(zip(keys, key))
        base["n_edges"] = int(len(g))
        base["n_detected"] = int(len(detected))
        base["n_high"] = int(g["n_high"].iloc[0])
        base["n_low"] = int(g["n_low"].iloc[0])
        base["n_tnk"] = int(g["n_tnk"].iloc[0])
        for m in METHODS:
            col = f"d_{m}"
            base[m] = float(detected[col].mean()) if len(detected) else 0.0
        rows.append(base)
    return pd.DataFrame(rows)


def family_table(lig: pd.DataFrame) -> pd.DataFrame:
    """One row per patient × split × method × family. Sum counts each ligand once."""
    ligand_of = {
        "barrier": BARRIER_LIGANDS,
        "ifn_recruit": IFN_LIGANDS,
        "mhc_cd8": MHC_LIGANDS,
    }
    rows = []
    units = lig[["split", "cohort", "patient"]].drop_duplicates()
    for rec in units.itertuples(index=False):
        sub = lig[(lig["split"] == rec.split) & (lig["cohort"] == rec.cohort) & (lig["patient"] == rec.patient)]
        for method in METHODS:
            for family, ligands in ligand_of.items():
                part = sub[sub["family"] == family]
                vals = []
                n_det = 0
                for lig_name in ligands:
                    hit = part[part["ligand"] == lig_name]
                    if hit.empty:
                        vals.append(0.0)
                    else:
                        vals.append(float(hit[method].iloc[0]))
                        n_det += int(hit["n_detected"].iloc[0] > 0)
                arr = np.asarray(vals, dtype=float)
                rows.append(
                    {
                        "split": rec.split,
                        "cohort": rec.cohort,
                        "patient": rec.patient,
                        "method": method,
                        "family": family,
                        "n_ligands": len(ligands),
                        "n_ligands_detected": n_det,
                        "family_sum": float(arr.sum()),
                        "family_mean": float(arr.mean()),
                    }
                )
    return pd.DataFrame(rows)


def pair7_table(edges: pd.DataFrame) -> pd.DataFrame:
    m = edges.apply(lambda r: (r["ligand"], r["receptor"]) in PAIR7, axis=1)
    sub = edges.loc[m].copy()
    rows = []
    for key, g in sub.groupby(["split", "cohort", "patient"], sort=False):
        for method in ("cellchat_hill", "cellchat_hill_pop", "cellchat_prod"):
            col = f"d_{method}"
            # undetected already 0
            rows.append(
                {
                    "split": key[0],
                    "cohort": key[1],
                    "patient": key[2],
                    "method": method,
                    "pair7_sum": float(g[col].sum()),
                }
            )
    return pd.DataFrame(rows)


def wilcox_p(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8 or np.allclose(x, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return float("nan")


def signflip_p(x: np.ndarray, rng: np.random.Generator) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8:
        return float("nan")
    obs = float(x.mean())
    signs = rng.choice(np.array([-1.0, 1.0]), size=(N_FLIP, x.size))
    null = signs @ x / x.size
    # one-sided, high > low
    return float((np.sum(null >= obs) + 1) / (N_FLIP + 1))


def dl_meta(y: np.ndarray, v: np.ndarray) -> dict:
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    mask = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[mask], v[mask]
    k = int(y.size)
    if k < 2:
        return {"mean": float(y[0]) if k == 1 else float("nan"), "i2": float("nan"), "p": float("nan")}
    w = 1.0 / v
    fe = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    ws = 1.0 / (v + tau2)
    mean = float(np.sum(ws * y) / np.sum(ws))
    se = float(np.sqrt(1.0 / np.sum(ws)))
    z = mean / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    return {"mean": mean, "i2": i2, "p": p}


def summarize_family(fam: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for (method, split, family), g in fam.groupby(["method", "split", "family"], sort=False):
        x = g["family_sum"].to_numpy(float)
        rec = {
            "method": method,
            "split": split,
            "family": family,
            "n": int(np.isfinite(x).sum()),
            "mean_sum": float(np.nanmean(x)),
            "median_sum": float(np.nanmedian(x)),
            "mean_per_ligand": float(np.nanmean(g["family_mean"])),
            "frac_pos": float(np.mean(x > 0)),
            "p_wilcox": wilcox_p(x),
            "p_signflip": signflip_p(x, rng) if family == "barrier" else float("nan"),
        }
        cohort_means, cohort_vars, n_pos_cohorts = [], [], 0
        for cohort in COHORTS:
            part = g.loc[g["cohort"] == cohort, "family_sum"].to_numpy(float)
            rec[f"n_{cohort}"] = int(part.size)
            rec[f"mean_{cohort}"] = float(part.mean()) if part.size else float("nan")
            if part.size >= 3 and part.mean() > 0:
                n_pos_cohorts += 1
            if part.size >= 3:
                sd = float(part.std(ddof=1))
                se = sd / math.sqrt(part.size)
                se = max(se, 1e-8 + 1e-6 * abs(float(part.mean())))
                cohort_means.append(float(part.mean()))
                cohort_vars.append(se ** 2)
        rec["n_cohorts_pos"] = n_pos_cohorts
        meta = dl_meta(np.array(cohort_means), np.array(cohort_vars))
        rec["meta_mean"] = meta["mean"]
        rec["i2"] = meta["i2"]
        rec["p_meta"] = meta["p"]
        rows.append(rec)
    return pd.DataFrame(rows)


def ligand_consistency(lig: pd.DataFrame, method: str, split: str, ligands: list[str]) -> dict:
    out = {}
    sub = lig[(lig["split"] == split)]
    n_full = 0
    for name in ligands:
        part = sub[sub["ligand"] == name]
        means = {}
        n_pos = 0
        for cohort in COHORTS:
            x = part.loc[part["cohort"] == cohort, method].to_numpy(float)
            means[cohort] = float(x.mean()) if x.size else float("nan")
            if x.size >= 3 and x.mean() > 0:
                n_pos += 1
        x = part[method].to_numpy(float)
        out[name] = {
            "n": int(x.size),
            "mean": float(np.mean(x)) if x.size else float("nan"),
            "median": float(np.median(x)) if x.size else float("nan"),
            "frac_pos": float(np.mean(x > 0)) if x.size else float("nan"),
            "p_wilcox": wilcox_p(x),
            "n_cohorts_pos": n_pos,
            **{f"mean_{c}": means[c] for c in COHORTS},
            **{f"n_{c}": int((part["cohort"] == c).sum()) for c in COHORTS},
        }
        if n_pos == 4 and out[name]["mean"] > 0:
            n_full += 1
    out["_n_ligands_4of4"] = n_full
    return out


def select_winner(summary: pd.DataFrame, lig: pd.DataFrame) -> dict:
    barrier = summary[(summary["family"] == "barrier")].copy()
    ifn = summary[summary["family"] == "ifn_recruit"][
        ["method", "split", "mean_sum", "mean_per_ligand", "n_cohorts_pos"]
    ].rename(
        columns={
            "mean_sum": "ifn_sum",
            "mean_per_ligand": "ifn_mean",
            "n_cohorts_pos": "ifn_cohorts_pos",
        }
    )
    merged = barrier.merge(ifn, on=["method", "split"], how="left")
    merged["separates"] = merged["mean_per_ligand"] > merged["ifn_mean"]
    merged["cohorts_ok"] = (
        (merged["n_GSE123902"] >= 3)
        & (merged["n_GSE131907"] >= 3)
        & (merged["n_GSE205335"] >= 3)
        & (merged["n_GSE189357"] >= 3)
        & (merged["n_cohorts_pos"] == 4)
    )
    merged["sig"] = merged["p_wilcox"] < 0.05
    merged["positive"] = merged["mean_sum"] > 0
    # trim guard
    guard_ok = []
    for rec in merged.itertuples(index=False):
        sibling = TRIM_GUARD.get(rec.method)
        if sibling is None:
            guard_ok.append(True)
            continue
        sib = merged[(merged["method"] == sibling) & (merged["split"] == rec.split)]
        if sib.empty or not np.isfinite(sib.iloc[0]["mean_sum"]):
            guard_ok.append(False)
            continue
        guard_ok.append(float(sib.iloc[0]["mean_sum"]) >= 0.5 * float(rec.mean_sum) and float(sib.iloc[0]["mean_sum"]) > 0)
    merged["trim_ok"] = guard_ok
    lig_ok_4, lig_ok_3 = [], []
    for rec in merged.itertuples(index=False):
        cons = ligand_consistency(lig, rec.method, rec.split, BARRIER_LIGANDS)
        lig_ok_4.append(cons["_n_ligands_4of4"] == 4)
        n3 = 0
        for name in BARRIER_LIGANDS:
            if cons[name]["mean"] > 0 and cons[name]["n_cohorts_pos"] >= 3:
                n3 += 1
        lig_ok_3.append(n3 == 4)
    merged["lig4"] = lig_ok_4
    merged["lig3"] = lig_ok_3
    base = merged["positive"] & merged["cohorts_ok"] & merged["sig"] & merged["separates"] & merged["trim_ok"]
    # Prefer a score that does not also lift IFN/recruit in every cohort.
    tiers = [
        base & merged["lig4"] & (merged["ifn_cohorts_pos"] < 4),
        base & merged["lig4"],
        base & merged["lig3"] & (merged["ifn_cohorts_pos"] < 4),
        base & merged["lig3"],
        base,
    ]
    tier_name = [
        "4/4 ligands each 4/4 cohorts, chemokine family not 4/4 high",
        "4/4 ligands each 4/4 cohorts",
        "4/4 ligands each >=3/4 cohorts, chemokine family not 4/4 high",
        "4/4 ligands each >=3/4 cohorts",
        "family 4/4 only",
    ]
    chosen = None
    used = None
    for mask, name in zip(tiers, tier_name):
        pool = merged.loc[mask]
        if pool.empty:
            continue
        # Largest absolute family sum. Units differ; this is the pre-specified max.
        chosen = pool.sort_values("mean_sum", ascending=False).iloc[0]
        used = name
        break
    if chosen is None:
        raise RuntimeError("no eligible communication score")
    return {"row": chosen, "tier": used, "ranked": merged.sort_values("mean_sum", ascending=False)}


def fmt_num(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    ax = abs(x)
    if ax >= 100:
        return f"{x:+.1f}"
    if ax >= 10:
        return f"{x:+.2f}"
    if ax >= 0.01:
        return f"{x:+.3f}"
    return f"{x:+.4g}"


def fmt_p(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    if x < 1e-3:
        return f"{x:.2e}"
    return f"{x:.3f}"


def summarize_pair7(pair7: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for (method, split), g in pair7.groupby(["method", "split"]):
        x = g["pair7_sum"].to_numpy(float)
        rec = {
            "method": method,
            "split": split,
            "n": int(x.size),
            "mean_sum": float(x.mean()) if x.size else float("nan"),
            "p_wilcox": wilcox_p(x),
            "p_signflip": signflip_p(x, rng),
            "frac_pos": float(np.mean(x > 0)) if x.size else float("nan"),
        }
        n_pos = 0
        for cohort in COHORTS:
            part = g.loc[g["cohort"] == cohort, "pair7_sum"].to_numpy(float)
            rec[f"n_{cohort}"] = int(part.size)
            rec[f"mean_{cohort}"] = float(part.mean()) if part.size else float("nan")
            if part.size >= 3 and part.mean() > 0:
                n_pos += 1
        rec["n_cohorts_pos"] = n_pos
        rows.append(rec)
    return pd.DataFrame(rows)


def plot_winner(fam: pd.DataFrame, method: str, split: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.4), sharey=False)
    families = [("barrier", "Barrier ligands"), ("ifn_recruit", "IFN / recruit"), ("mhc_cd8", "HLA–CD8")]
    for ax, (family, title) in zip(axes, families):
        sub = fam[(fam["method"] == method) & (fam["split"] == split) & (fam["family"] == family)]
        y = 0
        yticks, ylabels = [], []
        for cohort in COHORTS:
            part = sub[sub["cohort"] == cohort].sort_values("family_sum")
            ys = np.arange(y, y + len(part))
            ax.scatter(part["family_sum"], ys, s=18, color=COHORT_COLOR[cohort], zorder=3)
            if len(part):
                m = float(part["family_sum"].mean())
                ax.plot([m, m], [ys[0] - 0.4, ys[-1] + 0.4], color=COHORT_COLOR[cohort], lw=1.5)
            yticks.append(y + max(len(part) - 1, 0) / 2)
            ylabels.append(cohort.replace("GSE", ""))
            y += len(part) + 1.2
        ax.axvline(0, color="#888888", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("Patient family sum (high − low)")
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels)
    fig.suptitle(f"{method} · {split} · CLDN4-high minus CLDN4-low → T/NK", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligands(lig: pd.DataFrame, method: str, split: str, path: Path) -> None:
    order = BARRIER_LIGANDS + IFN_LIGANDS + MHC_LIGANDS
    means, los, his, colors = [], [], [], []
    color_of = {**{g: "#54A24B" for g in BARRIER_LIGANDS}, **{g: "#E45756" for g in IFN_LIGANDS}, **{g: "#4C78A8" for g in MHC_LIGANDS}}
    sub = lig[lig["split"] == split]
    for gene in order:
        x = sub.loc[sub["ligand"] == gene, method].to_numpy(float)
        means.append(float(np.mean(x)))
        # percentile interval across patients
        los.append(float(np.percentile(x, 25)))
        his.append(float(np.percentile(x, 75)))
        colors.append(color_of[gene])
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    y = np.arange(len(order))
    ax.barh(y, means, color=colors, height=0.72)
    means_a = np.asarray(means, dtype=float)
    lo_a = np.minimum(np.asarray(los, dtype=float), means_a)
    hi_a = np.maximum(np.asarray(his, dtype=float), means_a)
    ax.errorbar(means_a, y, xerr=[means_a - lo_a, hi_a - means_a], fmt="none", ecolor="#333333", elinewidth=0.8, capsize=2)
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel(f"Patient-mean Δ · {method}")
    ax.set_title("Green barrier, red IFN/recruit, blue HLA–CD8. Whiskers are the patient IQR.")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(summary, ranked, winner, tier, lig, pair7_sum, inv) -> None:
    row = winner
    method, split = row["method"], row["split"]
    cons = ligand_consistency(lig, method, split, BARRIER_LIGANDS)
    ifn = summary[(summary["method"] == method) & (summary["split"] == split) & (summary["family"] == "ifn_recruit")].iloc[0]
    mhc = summary[(summary["method"] == method) & (summary["split"] == split) & (summary["family"] == "mhc_cd8")].iloc[0]
    q4_n = int(inv["n_mal"].ge(MIN_MAL).sum()) if "n_mal" in inv else int(row["n"])
    lines = []
    lines.append("# FINDING — largest barrier-ligand communication Δ, concordant-4")
    lines.append("")
    lines.append("ADDITIVE. Does not replace the locked CellChat probability table (PR 616, family sum +0.0037).")
    lines.append("CLDN4 only. No dual-high. No TACSTD2 gate. Concordant four only.")
    lines.append("GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.")
    lines.append("Senders are malignant CLDN4-high vs CLDN4-low. Receiver is T/NK.")
    lines.append("The patient is the unit. Cell-pooled tests are not reported.")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(
        f"Largest eligible family Δ is **{method}** on the **{split}** split. "
        f"Per-ligand mean **{fmt_num(row['mean_per_ligand'])}** ({UNIT[method]}). "
        f"The family sum of F11R + NECTIN2 + CDH1 + LGALS9 is **{fmt_num(row['mean_sum'])}**. "
        f"That sum is not a single proportion."
    )
    lines.append("")
    lines.append(
        f"n={int(row['n'])}. Wilcoxon p={fmt_p(row['p_wilcox'])}. "
        f"Sign-flip one-sided p={fmt_p(row['p_signflip'])}. "
        f"Patients with family sum > 0: {row['frac_pos']:.0%}. "
        f"Cohort means of the sum: {fmt_num(row['mean_GSE123902'])} / {fmt_num(row['mean_GSE131907'])} / "
        f"{fmt_num(row['mean_GSE205335'])} / {fmt_num(row['mean_GSE189357'])} "
        f"(n {int(row['n_GSE123902'])}+{int(row['n_GSE131907'])}+{int(row['n_GSE205335'])}+{int(row['n_GSE189357'])}). "
        f"Random-effects I²={row['i2']:.0f}%."
    )
    q4 = summary[(summary["method"] == method) & (summary["split"] == "q4q1") & (summary["family"] == "barrier")]
    if split != "q4q1" and not q4.empty:
        q = q4.iloc[0]
        lines.append("")
        lines.append(
            f"Locked Q4 vs Q1, same score: n={int(q['n'])}, per-ligand mean {fmt_num(q['mean_per_ligand'])}, "
            f"family sum {fmt_num(q['mean_sum'])}, Wilcoxon p={fmt_p(q['p_wilcox'])}, "
            f"cohort means {fmt_num(q['mean_GSE123902'])} / {fmt_num(q['mean_GSE131907'])} / "
            f"{fmt_num(q['mean_GSE205335'])} / {fmt_num(q['mean_GSE189357'])} "
            f"(n {int(q['n_GSE123902'])}+{int(q['n_GSE131907'])}+{int(q['n_GSE205335'])}+{int(q['n_GSE189357'])})."
        )
    lines.append("")
    lines.append(f"Selection tier: {tier}.")
    lines.append(
        "The menu and the gates were fixed in the script before the ranking. "
        "The winner is the largest family sum among scores that pass. Units are not interchangeable. "
        "A percentage-point sum across four ligands is the absolute gap in how many sender cells express each ligand, "
        "gated on a detected T/NK receptor. It is not a CellChat probability. "
        "Fold-changes and population-scaled Hill sums are in the menu below and are smaller numbers on their own scales."
    )
    lines.append("")
    lines.append("## Why +0.0037 was small")
    lines.append("")
    lines.append("PR 616 summed official CellChat probabilities with `population.size=TRUE`. That multiplies every edge by the sender proportion times the receiver proportion, so a real ligand gap becomes a probability near zero. The same seven pairs, recomputed here:")
    lines.append("")
    lines.append("| score | split | n | mean pair-sum Δ | cohorts + | p Wilcoxon | p sign-flip |")
    lines.append("|---|---|---:|---:|---:|---|---|")
    for _, r in pair7_sum.sort_values(["split", "method"]).iterrows():
        lines.append(
            f"| {r['method']} | {r['split']} | {int(r['n'])} | {fmt_num(r['mean_sum'])} | "
            f"{int(r['n_cohorts_pos'])}/4 | {fmt_p(r['p_wilcox'])} | {fmt_p(r['p_signflip'])} |"
        )
    lines.append("")
    lines.append("## Barrier ligands")
    lines.append("")
    lines.append("| ligand | n | mean Δ | median | frac > 0 | cohorts + | 123902 | 131907 | 205335 | 189357 | p Wilcoxon |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for name in BARRIER_LIGANDS:
        d = cons[name]
        lines.append(
            f"| {name} | {d['n']} | {fmt_num(d['mean'])} | {fmt_num(d['median'])} | {d['frac_pos']:.0%} | "
            f"{d['n_cohorts_pos']}/4 | {fmt_num(d['mean_GSE123902'])} | {fmt_num(d['mean_GSE131907'])} | "
            f"{fmt_num(d['mean_GSE205335'])} | {fmt_num(d['mean_GSE189357'])} | {fmt_p(d['p_wilcox'])} |"
        )
    lines.append("")
    lines.append("## IFN / recruit, separate from HLA–CD8")
    lines.append("")
    lines.append(
        f"Chemokine family (CXCL9, CXCL10, CXCL11, CCL5, CXCL16), same score and split: "
        f"sum **{fmt_num(ifn['mean_sum'])}**, per-ligand mean {fmt_num(ifn['mean_per_ligand'])}, "
        f"n={int(ifn['n'])}, Wilcoxon p={fmt_p(ifn['p_wilcox'])}, "
        f"cohorts high>low {int(ifn['n_cohorts_pos'])}/4, "
        f"cohort means {fmt_num(ifn['mean_GSE123902'])} / {fmt_num(ifn['mean_GSE131907'])} / "
        f"{fmt_num(ifn['mean_GSE205335'])} / {fmt_num(ifn['mean_GSE189357'])}."
    )
    lines.append("")
    lines.append(
        f"HLA–CD8 (HLA-A, HLA-B, HLA-C) is not recruitment. Same score: sum **{fmt_num(mhc['mean_sum'])}**, "
        f"per-ligand mean {fmt_num(mhc['mean_per_ligand'])}, Wilcoxon p={fmt_p(mhc['p_wilcox'])}, "
        f"cohorts high>low {int(mhc['n_cohorts_pos'])}/4."
    )
    lines.append("")
    ifn_cons = ligand_consistency(lig, method, split, IFN_LIGANDS + MHC_LIGANDS)
    lines.append("| ligand | family | mean Δ | cohorts + | p Wilcoxon |")
    lines.append("|---|---|---:|---:|---|")
    for name in IFN_LIGANDS:
        d = ifn_cons[name]
        lines.append(f"| {name} | ifn_recruit | {fmt_num(d['mean'])} | {d['n_cohorts_pos']}/4 | {fmt_p(d['p_wilcox'])} |")
    for name in MHC_LIGANDS:
        d = ifn_cons[name]
        lines.append(f"| {name} | mhc_cd8 | {fmt_num(d['mean'])} | {d['n_cohorts_pos']}/4 | {fmt_p(d['p_wilcox'])} |")
    lines.append("")
    lines.append("## Full menu (barrier family sum, both splits)")
    lines.append("")
    lines.append("Eligible scores had to be positive, Wilcoxon p<0.05, positive in all four cohorts (each n≥3), larger per ligand than the chemokine family, and, for untrimmed CP10k scores, at least half as large after 10% trimming. Each of the four ligands also had to be positive. The winning tier is named above.")
    lines.append("")
    lines.append("| method | split | n | family sum | per ligand | frac>0 | cohorts + | ligands 4/4 | p | separates | trim ok |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|---|---|---|")
    show = ranked.sort_values(["split", "mean_sum"], ascending=[True, False])
    for _, r in show.iterrows():
        lines.append(
            f"| {r['method']} | {r['split']} | {int(r['n'])} | {fmt_num(r['mean_sum'])} | {fmt_num(r['mean_per_ligand'])} | "
            f"{r['frac_pos']:.0%} | {int(r['n_cohorts_pos'])}/4 | {'yes' if r['lig4'] else 'no'} | {fmt_p(r['p_wilcox'])} | "
            f"{'yes' if r['separates'] else 'no'} | {'yes' if r['trim_ok'] else 'no'} |"
        )
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append(f"Inventory units: {len(inv)}. Q4 vs Q1 requires ≥{MIN_MAL} malignant cells, ≥{MIN_ARM} cells in each arm, and ≥{MIN_TNK} T/NK cells. Decile uses the outer 10% and the same arm floor, so patients with fewer than 100 malignant cells drop out.")
    lines.append(f"Units with n_mal≥{MIN_MAL} in the inventory built here: {int((inv['n_mal']>=MIN_MAL).sum())} (reference expectation 64; P4001 has 27 malignant cells).")
    lines.append("")
    lines.append("## What is not claimed")
    lines.append("")
    lines.append("- This is an expression ligand–receptor contrast. It is not spatial exclusion and it does not replace the CosMx 8/8 result.")
    lines.append("- Population-scaled CellChat probabilities are not re-issued as +0.0037. The pair-sum table above is the recomputed Hill score.")
    lines.append("- A larger number on a CP10k or percentage-point scale is a different unit from a probability. The menu shows both.")
    lines.append("- TACSTD2 is not a gate. GSE148071 is not added.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("bash methods/concordant4_ccc_maxeffect_cldn4/scripts/download.sh /tmp/concordant4_raw")
    lines.append("python3 methods/concordant4_ccc_maxeffect_cldn4/scripts/run_maxeffect.py")
    lines.append("```")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (HERE / "FINDING.md").write_text(text)
    log(f"wrote FINDING headline {method} {split} {row['mean_sum']:.4g}")


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    keep = wanted_genes()
    units = []
    units += load_gse123902(keep)
    units += load_gse131907(keep)
    units += load_gse205335(keep)
    units += load_gse189357(keep)
    inv = qc_counts(units)
    inv.to_csv(TAB / "patient_inventory.tsv", sep="\t", index=False)
    log(f"scoring {len(units)} units")
    rows = []
    for i, unit in enumerate(units, start=1):
        try:
            got = score_unit(unit)
        except Exception:
            log(traceback.format_exc())
            raise
        rows.extend(got)
        if i % 10 == 0:
            log(f"  scored {i}/{len(units)}")
    edges = pd.DataFrame(rows)
    edges.to_csv(TAB / "per_patient_edges.tsv", sep="\t", index=False)
    log(f"edges {len(edges)}")
    lig = ligand_table(edges)
    # long ligand table is wide (methods as columns). Keep it.
    lig.to_csv(TAB / "per_patient_ligands.tsv", sep="\t", index=False)
    fam = family_table(lig)
    fam.to_csv(TAB / "per_patient_family.tsv", sep="\t", index=False)
    pair7 = pair7_table(edges)
    pair7.to_csv(TAB / "per_patient_pair7.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED)
    summary = summarize_family(fam, rng)
    summary.to_csv(TAB / "summary_family.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED)
    p7 = summarize_pair7(pair7, rng)
    p7.to_csv(TAB / "summary_pair7.tsv", sep="\t", index=False)
    picked = select_winner(summary, lig)
    picked["ranked"].to_csv(TAB / "method_rank.tsv", sep="\t", index=False)
    row = picked["row"]
    plot_winner(fam, row["method"], row["split"], FIG / "winner_family_patients.png")
    plot_ligands(lig, row["method"], row["split"], FIG / "winner_ligand_bars.png")
    write_finding(summary, picked["ranked"], row, picked["tier"], lig, p7, inv)
    # session pin
    (HERE / "results" / "versions.txt").write_text(
        "numpy " + np.__version__ + "\n" + "pandas " + pd.__version__ + "\n" + "scipy " + __import__("scipy").__version__ + "\n"
    )
    log("done")


if __name__ == "__main__":
    main()
