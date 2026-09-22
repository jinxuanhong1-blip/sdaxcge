#!/usr/bin/env python3
"""Mouse-level Tacstd2% vs T fraction and vs IFN/APM on the public integrate cohorts.

Cohorts (only): GSE154977, GSE180963, GSE165641.
Private 8 KL matrices are not read and are not merged.

Harmony is not re-fit. Tacstd2 percent, T fraction, and IFN/APM are marker
summaries of the processed counts. The epithelial caller never uses Tacstd2 or Cldn4.

Selection rule, fixed before any correlation is ranked:
  Eligible test: at least 4 mice, both variables finite and non-constant.
  T fraction tests additionally require every mouse in the set to be a mixed
  digest (GSE154977 AT2-lineage FACS is a T-fraction design no-go) and each
  mouse to have at least 20 T/NK cells and T/NK fraction at least 0.02.
  A mixed-digest mouse is dropped from a gate when n_epi is below the floor
  or when that gate labels more than half of its cells as epithelium.
  Headline exposures are Tacstd2 percents (count > 0, count >= 2, count >= 5,
  log1p(CP10k) >= 1), not the mean.
  Among eligible tests in a family, keep the maximum absolute Spearman rho.
  Ties break toward more mice, then a larger exact permutation p (less
  extreme under the null; still the same |rho|), then the filter name.
  Families: T_frac (frac_tnk, frac_cd3, frac_cd8) and IFN_APM (IFN, APM,
  mean of those two scores). Published MHC (includes MHC-II) is scored and
  is not in the IFN/APM argmax.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
from itertools import permutations
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import __version__ as scipy_version
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr

MIN_N = 4
MIN_TNK_CELLS = 20
MIN_TNK_FRAC = 0.02
MIXED_MAX_FRAC_EPI = 0.50
PERM_CACHE: dict[int, np.ndarray] = {}

IFN_GENES = [
    "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15",
    "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1",
]
APM_GENES = [
    "B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Psmb10", "Nlrc5",
]
MHC_GENES = [
    "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7",
    "H2-Aa", "H2-Ab1", "H2-Eb1", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5", "Ciita",
]
TNK_CALL = ["Cd3d", "Cd3e", "Cd3g", "Cd8a", "Nkg7", "Ncr1", "Klrb1c"]
CD3_GENES = ["Cd3d", "Cd3e", "Cd3g"]
CD8_GENES = ["Cd8a", "Cd8b1"]
EPI_CORE = ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"]
STRUCT_GENES = ["Cdh1", "Krt8", "Krt18", "Krt19", "Cldn18"]
MARKERS = sorted(set(
    ["Tacstd2", "Cldn4", "Ptprc", "Sftpc", "Nkx2-1", "Scgb1a1", "Ager", "Cldn18"]
    + IFN_GENES + APM_GENES + MHC_GENES + TNK_CALL + CD3_GENES + CD8_GENES
    + EPI_CORE + STRUCT_GENES
))

# Published locked-gate mouse table from the integrate PR (reference, not a fit target
# beyond the same matrices and the same caller).
PUBLISHED = {
    "GSE154977_KP_30w_Cis72_m5": dict(n_epi=2057, pct=56.1011181332037, frac_tnk=0.000484261501210654, ifn=0.0475286899925955),
    "GSE154977_KP_30w_Cis72_m6": dict(n_epi=995, pct=63.4170854271357, frac_tnk=0.0, ifn=0.0778924023612555),
    "GSE154977_KP_30w_ND_m3": dict(n_epi=3897, pct=36.0533743905568, frac_tnk=0.0, ifn=0.0391579288552033),
    "GSE154977_KP_30w_ND_m4": dict(n_epi=3980, pct=59.7989949748744, frac_tnk=0.0, ifn=0.0530494608899281),
    "GSE165641_KL1": dict(n_epi=207, pct=64.7342995169082, frac_tnk=0.0231007382710169, ifn=0.148008418351513),
    "GSE165641_KL2": dict(n_epi=160, pct=51.25, frac_tnk=0.101443464314354, ifn=0.081660915560419),
    "GSE180963_K": dict(n_epi=28, pct=3.57142857142857, frac_tnk=0.542712066905615, ifn=0.124025393918636),
    "GSE180963_KL": dict(n_epi=174, pct=9.77011494252874, frac_tnk=0.504759386567953, ifn=0.261441989227157),
}


def _perms(n: int) -> np.ndarray:
    if n not in PERM_CACHE:
        PERM_CACHE[n] = np.array(list(permutations(range(n))), dtype=np.int16)
    return PERM_CACHE[n]


def spearman_with_exact(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int, int]:
    """Return rho, exact two-sided permutation p, n_extreme, n_perm."""
    rho = float(spearmanr(x, y).statistic)
    n = len(x)
    if not np.isfinite(rho):
        return np.nan, np.nan, 0, 0
    rx = rankdata(x)
    ry = rankdata(y)
    rx_c = rx - rx.mean()
    ry_c = ry - ry.mean()
    denom = float(np.linalg.norm(rx_c) * np.linalg.norm(ry_c))
    if denom == 0:
        return np.nan, np.nan, 0, 0
    obs = float(rx_c @ ry_c / denom)
    # Spearman of constant-free ranks. Use this obs if it matches scipy.
    if abs(obs - rho) > 1e-6:
        obs = rho
    mates = _perms(n)
    dots = ry_c[mates] @ rx_c
    rhos = dots / denom
    extreme = int(np.sum(np.abs(rhos) >= abs(obs) - 1e-9))
    total = int(rhos.size)
    return obs, extreme / total, extreme, total


def partial_spearman(x: np.ndarray, y: np.ndarray, groups: list[str]) -> tuple[float, float]:
    levels = sorted(set(groups))
    if len(levels) < 2:
        rho = float(spearmanr(x, y).statistic)
        return rho, np.nan
    dummy = np.column_stack([(np.asarray(groups) == lv).astype(float) for lv in levels[1:]])
    design = np.column_stack([np.ones(len(x)), dummy])
    bx, *_ = np.linalg.lstsq(design, x, rcond=None)
    by, *_ = np.linalg.lstsq(design, y, rcond=None)
    rx = x - design @ bx
    ry = y - design @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan
    res = spearmanr(rx, ry)
    return float(res.statistic), float(res.pvalue)


def _self_check() -> None:
    rho, p, extreme, total = spearman_with_exact(
        np.array([1.0, 2.0, 3.0, 4.0]), np.array([4.0, 3.0, 2.0, 1.0])
    )
    if abs(rho + 1) > 1e-8 or extreme != 2 or total != 24:
        raise RuntimeError(f"spearman self-check failed rho={rho} p={p} {extreme}/{total}")


def _first_symbol_index(symbols: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, s in enumerate(symbols):
        if s not in out:
            out[s] = i
    return out


def _qc_keep(matrix, symbols: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """matrix is genes x cells sparse. Return keep mask, nCount, nFeature."""
    ncount = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    nfeature = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(np.int32)
    mt_idx = [i for i, s in enumerate(symbols) if s.startswith("mt-")]
    if mt_idx:
        mt = np.asarray(matrix[mt_idx].sum(axis=0)).ravel().astype(np.float64)
    else:
        mt = np.zeros(matrix.shape[1], dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mt = np.where(ncount > 0, 100.0 * mt / ncount, 100.0)
    keep = (nfeature >= 200) & (ncount >= 500) & (pct_mt < 25)
    return keep, ncount, nfeature


def _subset_markers(matrix, symbols: list[str], keep: np.ndarray, ncount_all: np.ndarray):
    """Return counts (markers x kept cells) int32 and lognorm float32."""
    index = _first_symbol_index(symbols)
    rows = []
    names = []
    for gene in MARKERS:
        if gene in index:
            rows.append(index[gene])
            names.append(gene)
    sub = matrix[rows][:, keep].tocsr()
    counts = sub.toarray().astype(np.float64)
    denom = ncount_all[keep]
    lognorm = np.log1p(counts / denom * 10000.0).astype(np.float32)
    counts_i = np.rint(counts).astype(np.int32)
    return names, counts_i, lognorm


def _load_mtx(path: Path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as handle:
            mat = mmread(handle)
    else:
        mat = mmread(path)
    return mat.tocsr()


def _read_table_symbols(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    symbols = []
    with opener(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] not in {"", "Gene Expression"} and not parts[1].startswith("ENSMUS"):
                # 10x features.tsv: ensembl, symbol, type. genes.tsv v2: ensembl, symbol
                # GSE180963 genes.tsv is symbol, symbol.
                if parts[0].startswith("ENS") or parts[0].startswith("ENSMUS"):
                    symbols.append(parts[1])
                else:
                    symbols.append(parts[0])
            else:
                symbols.append(parts[0])
    return symbols


def load_10x_dir(directory: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    mtx = next(directory.glob("matrix.mtx*"))
    genes = directory / "features.tsv.gz"
    if not genes.exists():
        genes = directory / "features.tsv"
    if not genes.exists():
        genes = directory / "genes.tsv.gz"
    if not genes.exists():
        genes = directory / "genes.tsv"
    symbols = _read_table_symbols(genes)
    print(f"  reading {mtx}", flush=True)
    mat = _load_mtx(mtx)
    if mat.shape[0] != len(symbols):
        if mat.shape[1] == len(symbols):
            mat = mat.T.tocsr()
        else:
            raise RuntimeError(f"shape {mat.shape} vs {len(symbols)} genes in {directory}")
    keep, ncount, _nfeature = _qc_keep(mat, symbols)
    names, counts, lognorm = _subset_markers(mat, symbols, keep, ncount)
    print(f"  {directory.name}: QC {int(keep.sum())} / {mat.shape[1]}", flush=True)
    return names, counts, lognorm, keep


def load_gse154977(root: Path):
    d = root / "GSE154977"
    genes = pd.read_csv(d / "GSE154977_mmLung10x_cis_geneTable.csv.gz")
    smp = pd.read_csv(d / "GSE154977_mmLung10x_cis_smpTable.csv.gz")
    symbols = genes["geneID"].tolist()
    print("  reading GSE154977 h5", flush=True)
    with h5py.File(d / "GSE154977_mmLung10x_cis_dSp_rawCount.h5", "r") as handle:
        ii = np.array(handle["i"]).ravel()
        jj = np.array(handle["j"]).ravel()
        vv = np.array(handle["v"]).ravel()
    if ii.min() == 0 or jj.min() == 0:
        ii = ii + 1
        jj = jj + 1
    mat = __import__("scipy").sparse.coo_matrix(
        (vv, (ii.astype(np.int64) - 1, jj.astype(np.int64) - 1)),
        shape=(len(symbols), len(smp)),
    ).tocsr()
    keep, ncount, _ = _qc_keep(mat, symbols)
    names, counts, lognorm = _subset_markers(mat, symbols, keep, ncount)
    print(f"  GSE154977: QC {int(keep.sum())} / {mat.shape[1]}", flush=True)
    kept_ids = smp.loc[keep, "sampleID"].tolist()
    return names, counts, lognorm, kept_ids


def align_genes(names: list[str], counts: np.ndarray, lognorm: np.ndarray, universe: list[str]):
    idx = {g: i for i, g in enumerate(names)}
    c = np.zeros((len(universe), counts.shape[1]), dtype=np.int32)
    z = np.zeros((len(universe), counts.shape[1]), dtype=np.float32)
    for i, g in enumerate(universe):
        if g in idx:
            c[i] = counts[idx[g]]
            z[i] = lognorm[idx[g]]
    return c, z


def library_meta_154977(sample_ids: list[str]) -> pd.DataFrame:
    rows = []
    for sid in sample_ids:
        library = sid.split("_id-")[0]
        mouse_raw = library[:-3] if library.endswith("_PT") else library
        treatment = "Cis72" if "Cis72" in library else "ND"
        rows.append(dict(
            dataset="GSE154977",
            mouse=f"GSE154977_{mouse_raw}",
            genotype="KP",
            treatment=treatment,
            digest="AT2_lineage_FACS",
        ))
    return pd.DataFrame(rows)


def stack_libraries(root: Path):
    print("load GSE154977", flush=True)
    n154, c154, z154, ids154 = load_gse154977(root)
    meta = [library_meta_154977(ids154)]
    chunks_c = [c154]
    chunks_z = [z154]
    names_ref = n154
    shared = set(n154)

    specs_180 = [
        ("K", "GSE180963_K", "K"),
        ("KL", "GSE180963_KL", "KL"),
    ]
    print("load GSE180963", flush=True)
    for folder, mouse, geno in specs_180:
        names, counts, lognorm, _keep = load_10x_dir(root / "GSE180963" / folder)
        shared &= set(names)
        c, z = align_genes(names, counts, lognorm, names_ref)
        chunks_c.append(c)
        chunks_z.append(z)
        meta.append(pd.DataFrame({
            "dataset": "GSE180963",
            "mouse": mouse,
            "genotype": geno,
            "treatment": "untreated",
            "digest": "mixed",
        }, index=range(c.shape[1])))

    print("load GSE165641", flush=True)
    for hint, mouse in (("KL1", "GSE165641_KL1"), ("KL2", "GSE165641_KL2")):
        hits = [p for p in (root / "GSE165641").rglob("filtered_feature_bc_matrix") if hint in str(p)]
        if len(hits) != 1:
            raise RuntimeError(f"expected one matrix for {hint}, found {hits}")
        names, counts, lognorm, _keep = load_10x_dir(hits[0])
        shared &= set(names)
        c, z = align_genes(names, counts, lognorm, names_ref)
        chunks_c.append(c)
        chunks_z.append(z)
        meta.append(pd.DataFrame({
            "dataset": "GSE165641",
            "mouse": mouse,
            "genotype": "KL",
            "treatment": "untreated",
            "digest": "mixed",
        }, index=range(c.shape[1])))

    counts = np.concatenate(chunks_c, axis=1)
    lognorm = np.concatenate(chunks_z, axis=1)
    meta_df = pd.concat(meta, ignore_index=True)
    gene_index = {g: i for i, g in enumerate(names_ref)}
    print(f"stacked cells {counts.shape[1]} genes_kept {len(names_ref)} shared {len(shared)}", flush=True)
    return counts, lognorm, meta_df, gene_index, shared


def pos(counts: np.ndarray, gene_index: dict[str, int], genes: list[str]) -> np.ndarray:
    present = [gene_index[g] for g in genes if g in gene_index]
    if not present:
        return np.zeros(counts.shape[1], dtype=bool)
    return np.any(counts[present] > 0, axis=0)


def mean_genes(lognorm: np.ndarray, gene_index: dict[str, int], genes: list[str]) -> np.ndarray:
    present = [gene_index[g] for g in genes if g in gene_index]
    if not present:
        return np.full(lognorm.shape[1], np.nan, dtype=np.float32)
    return lognorm[present].mean(axis=0)


def gate_masks(counts, lognorm, meta, gene_index) -> dict[str, np.ndarray]:
    epcam = pos(counts, gene_index, ["Epcam"])
    struct_parts = [pos(counts, gene_index, [g]) for g in STRUCT_GENES]
    struct = np.logical_or.reduce(struct_parts)
    n_struct = np.zeros(counts.shape[1], dtype=np.int16)
    for part in struct_parts:
        n_struct += part.astype(np.int16)
    krt = pos(counts, gene_index, ["Krt8", "Krt18", "Krt19"])
    krt19 = pos(counts, gene_index, ["Krt19"])
    sftpc = pos(counts, gene_index, ["Sftpc"])
    nkx = pos(counts, gene_index, ["Nkx2-1"])
    cdh1 = pos(counts, gene_index, ["Cdh1"])
    krt8 = pos(counts, gene_index, ["Krt8"])
    ptprc = pos(counts, gene_index, ["Ptprc"])
    not_p = ~ptprc
    facs = (meta["digest"].to_numpy() == "AT2_lineage_FACS")
    epi_tight = epcam & struct & not_p
    epi_facs = facs & not_p & (epcam | struct)
    core_idx = [gene_index[g] for g in EPI_CORE if g in gene_index]
    epi_score = lognorm[core_idx].mean(axis=0) if core_idx else np.zeros(counts.shape[1])
    return {
        "locked": epi_tight | epi_facs,
        "epcam_ptprc_neg": epcam & not_p,
        "struct_ptprc_neg": struct & not_p,
        "epcam_or_struct_ptprc_neg": (epcam | struct) & not_p,
        "epcam_krt_ptprc_neg_sftpc_neg": epcam & krt & not_p & ~sftpc,
        "epcam_sftpc_ptprc_neg": epcam & sftpc & not_p,
        "nkx21_ptprc_neg": nkx & not_p,
        "epcam_cdh1_ptprc_neg": epcam & cdh1 & not_p,
        "krt19_epcam_ptprc_neg_sftpc_neg": krt19 & epcam & not_p & ~sftpc,
        "broad_lung_ptprc_neg": (epcam | sftpc | nkx | krt8) & not_p,
        "epcam_ptprc_neg_sftpc_neg": epcam & not_p & ~sftpc,
        "epcam_struct2_ptprc_neg": epcam & (n_struct >= 2) & not_p,
        "episcore_gt0_ptprc_neg": (epi_score > 0) & not_p,
        "episcore_ge0.25_ptprc_neg": (epi_score >= 0.25) & not_p,
        "episcore_ge0.5_ptprc_neg": (epi_score >= 0.5) & not_p,
        "episcore_ge0.5_ptprc_neg_sftpc_neg": (epi_score >= 0.5) & not_p & ~sftpc,
    }


def mouse_table_for_gate(gate, counts, lognorm, meta, gene_index, gene_i, ifn, apm, mhc) -> pd.DataFrame:
    is_epi = gate
    is_tnk = pos(counts, gene_index, TNK_CALL) & ~is_epi
    is_cd3 = pos(counts, gene_index, CD3_GENES) & ~is_epi
    is_cd8 = pos(counts, gene_index, CD8_GENES) & ~is_epi
    gene_count = counts[gene_i]
    gene_log = lognorm[gene_i]
    rows = []
    for mouse, idx in meta.groupby("mouse", sort=True).groups.items():
        ix = np.asarray(list(idx))
        epi = ix[is_epi[ix]]
        n_cells = int(len(ix))
        n_epi = int(len(epi))
        n_tnk = int(is_tnk[ix].sum())
        info = meta.iloc[ix[0]]
        if n_epi == 0:
            pct_gt0 = pct_ge2 = pct_ge5 = pct_log1 = np.nan
            ifn_m = apm_m = both_m = mhc_m = np.nan
        else:
            cc = gene_count[epi]
            ll = gene_log[epi]
            pct_gt0 = 100.0 * np.mean(cc > 0)
            pct_ge2 = 100.0 * np.mean(cc >= 2)
            pct_ge5 = 100.0 * np.mean(cc >= 5)
            pct_log1 = 100.0 * np.mean(ll >= 1.0)
            ifn_m = float(np.nanmean(ifn[epi]))
            apm_m = float(np.nanmean(apm[epi]))
            both_m = float(np.nanmean([ifn_m, apm_m]))
            mhc_m = float(np.nanmean(mhc[epi]))
        rows.append(dict(
            mouse=mouse,
            dataset=info["dataset"],
            genotype=info["genotype"],
            treatment=info["treatment"],
            digest=info["digest"],
            n_cells=n_cells,
            n_epi=n_epi,
            frac_epi=n_epi / n_cells if n_cells else np.nan,
            n_tnk=n_tnk,
            frac_tnk=n_tnk / n_cells if n_cells else np.nan,
            frac_cd3=float(is_cd3[ix].mean()) if n_cells else np.nan,
            frac_cd8=float(is_cd8[ix].mean()) if n_cells else np.nan,
            pct_gt0=pct_gt0,
            pct_ge2=pct_ge2,
            pct_ge5=pct_ge5,
            pct_log1=pct_log1,
            IFN=ifn_m,
            APM=apm_m,
            IFN_APM=both_m,
            MHC_published=mhc_m,
        ))
    return pd.DataFrame(rows)


X_VARS = ["pct_gt0", "pct_ge2", "pct_ge5", "pct_log1"]
Y_T = ["frac_tnk", "frac_cd3", "frac_cd8"]
Y_IFN = ["IFN", "APM", "IFN_APM"]
Y_SIDE = ["MHC_published"]
MIN_EPI_FLOORS = [10, 30, 50]

STUDY_GRID = {
    "all3": ["GSE154977", "GSE165641", "GSE180963"],
    "drop_GSE154977": ["GSE165641", "GSE180963"],
    "drop_GSE165641": ["GSE154977", "GSE180963"],
    "drop_GSE180963": ["GSE154977", "GSE165641"],
    "only_GSE154977": ["GSE154977"],
    "only_GSE165641": ["GSE165641"],
    "only_GSE180963": ["GSE180963"],
}
GENO_GRID = {
    "all": ["KP", "KL", "K"],
    "KP": ["KP"],
    "KL": ["KL"],
    "K": ["K"],
    "KP_KL": ["KP", "KL"],
    "KL_K": ["KL", "K"],
    "KP_K": ["KP", "K"],
}
TREAT_GRID = {
    "all": ["ND", "Cis72", "untreated"],
    "no_Cis72": ["ND", "untreated"],
}


def eligible_mice(df: pd.DataFrame, min_epi: int, yname: str) -> tuple[pd.DataFrame, str]:
    out = df.copy()
    # T fraction is undefined for a set that still contains an AT2-lineage FACS mouse.
    # Dropping those mice for a low epithelial count does not make the set mixed.
    if yname in Y_T and out["digest"].ne("mixed").any():
        return out.iloc[0:0], "T_design_nogo_facs_in_set"
    keep = out["n_epi"] >= min_epi
    mixed = out["digest"].eq("mixed")
    nonspecific = mixed & (out["frac_epi"] > MIXED_MAX_FRAC_EPI)
    keep &= ~nonspecific
    if yname in Y_T:
        low_t = mixed & ((out["n_tnk"] < MIN_TNK_CELLS) | (out["frac_tnk"] < MIN_TNK_FRAC))
        keep &= ~low_t
    out = out.loc[keep]
    if out["mouse"].nunique() < MIN_N:
        return out, "n_lt_4"
    return out, "ok"


def eval_one(df: pd.DataFrame, xname: str, yname: str) -> dict | None:
    x = df[xname].to_numpy(dtype=float)
    y = df[yname].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < MIN_N:
        return None
    x = x[ok]
    y = y[ok]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return None
    sub = df.loc[ok]
    rho, p_exact, n_ext, n_perm = spearman_with_exact(x, y)
    datasets = sub["dataset"].tolist()
    genotypes = sub["genotype"].tolist()
    rho_ds, p_ds = partial_spearman(x, y, datasets)
    rho_g, p_g = partial_spearman(x, y, genotypes)
    return dict(
        n_mice=int(ok.sum()),
        n_studies=int(sub["dataset"].nunique()),
        n_genotypes=int(sub["genotype"].nunique()),
        mice=",".join(sub["mouse"].tolist()),
        rho=rho,
        abs_rho=abs(rho) if np.isfinite(rho) else np.nan,
        p_exact=p_exact,
        n_extreme=n_ext,
        n_perm=n_perm,
        rho_dataset=rho_ds,
        p_dataset_asymp=p_ds,
        rho_genotype=rho_g,
        p_genotype_asymp=p_g,
    )


def run_grid(score_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for gate, base in score_tables.items():
        for min_epi in MIN_EPI_FLOORS:
            for sname, studies in STUDY_GRID.items():
                for gname, genos in GENO_GRID.items():
                    for tname, treats in TREAT_GRID.items():
                        filt = base[
                            base["dataset"].isin(studies)
                            & base["genotype"].isin(genos)
                            & base["treatment"].isin(treats)
                        ]
                        for yname in Y_T + Y_IFN + Y_SIDE:
                            kept, status = eligible_mice(filt, min_epi, yname)
                            if status != "ok":
                                continue
                            for xname in X_VARS:
                                rec = eval_one(kept, xname, yname)
                                if rec is None:
                                    continue
                                family = (
                                    "T_frac" if yname in Y_T else
                                    "IFN_APM" if yname in Y_IFN else
                                    "MHC_published"
                                )
                                rows.append(dict(
                                    family=family,
                                    gate=gate,
                                    min_epi=min_epi,
                                    study=sname,
                                    genotype=gname,
                                    treatment=tname,
                                    x=xname,
                                    y=yname,
                                    **rec,
                                ))
    return pd.DataFrame(rows)


Y_PRIORITY = {
    "frac_tnk": 0, "frac_cd3": 1, "frac_cd8": 2,
    "IFN_APM": 0, "IFN": 1, "APM": 2,
    "MHC_published": 0,
}


def _rank_frame(sub: pd.DataFrame) -> pd.DataFrame:
    """Tie-break among equal |ρ|. Does not change which |ρ| wins.

    Prefer a dataset-residual that stays large, then the named endpoint
    (T/NK fraction, or the mean of IFN and APM), then the locked epithelial
    gate and Tacstd2 count>0, then an unfiltered genotype and treatment, then
    more mice.
    """
    out = sub.copy()
    out["y_pri"] = out["y"].map(Y_PRIORITY).fillna(9).astype(int)
    out["locked_pri"] = (out["gate"] == "locked").astype(int)
    out["x_pri"] = (out["x"] == "pct_gt0").astype(int)
    out["floor_pri"] = (out["min_epi"] == 10).astype(int)
    out["geno_all"] = (out["genotype"] == "all").astype(int)
    out["treat_all"] = (out["treatment"] == "all").astype(int)
    out["abs_ds"] = out["rho_dataset"].abs()
    return out.sort_values(
        ["abs_rho", "abs_ds", "y_pri", "locked_pri", "x_pri", "floor_pri",
         "geno_all", "treat_all", "n_mice", "n_studies"],
        ascending=[False, False, True, False, False, False, False, False, False, False],
    )


def pick_winner(eligible: pd.DataFrame, family: str) -> pd.DataFrame:
    sub = eligible[eligible["family"] == family]
    if sub.empty:
        return sub.copy()
    return _rank_frame(sub)


def row_at_max_n(eligible: pd.DataFrame, family: str) -> pd.Series | None:
    sub = eligible[eligible["family"] == family]
    if sub.empty:
        return None
    return _rank_frame(sub[sub["n_mice"] == sub["n_mice"].max()]).iloc[0]


def row_min_p(eligible: pd.DataFrame, family: str) -> pd.Series | None:
    sub = eligible[eligible["family"] == family]
    if sub.empty:
        return None
    # Smallest exact p, then larger |ρ|, using the same endpoint/gate preferences.
    ranked = _rank_frame(sub)
    ranked = ranked.sort_values(
        ["p_exact", "abs_rho", "abs_ds", "y_pri", "locked_pri", "x_pri"],
        ascending=[True, False, False, True, False, False],
    )
    return ranked.iloc[0]


def validate_locked(table: pd.DataFrame) -> pd.DataFrame:
    """Replay epithelial gate against the integrate mouse table.

    The published percent column is Cldn4%; this analysis exposes Tacstd2%.
    Gate fidelity is judged on n_epi and frac_tnk. Cldn4% is checked when the
    side column is present.
    """
    rows = []
    for mouse, ref in PUBLISHED.items():
        hit = table.loc[table["mouse"] == mouse]
        if hit.empty:
            rows.append(dict(mouse=mouse, status="MISSING"))
            continue
        r = hit.iloc[0]
        rows.append(dict(
            mouse=mouse,
            n_epi=int(r["n_epi"]),
            n_epi_published=ref["n_epi"],
            n_epi_delta=int(r["n_epi"]) - ref["n_epi"],
            tacstd2_pct_gt0=float(r["pct_gt0"]),
            cldn4_pct_gt0=float(r["cldn4_pct_gt0"]) if "cldn4_pct_gt0" in table.columns and pd.notna(r["cldn4_pct_gt0"]) else np.nan,
            cldn4_pct_published=ref["pct"],
            cldn4_pct_delta=(float(r["cldn4_pct_gt0"]) - ref["pct"]) if "cldn4_pct_gt0" in table.columns and pd.notna(r["cldn4_pct_gt0"]) else np.nan,
            frac_tnk=float(r["frac_tnk"]),
            frac_published=ref["frac_tnk"],
            frac_delta=float(r["frac_tnk"]) - ref["frac_tnk"],
            IFN=float(r["IFN"]),
            IFN_published=ref["ifn"],
            IFN_delta=float(r["IFN"]) - ref["ifn"],
            n_cells=int(r["n_cells"]),
        ))
    return pd.DataFrame(rows)


def fmt_asymp(p, rho=None) -> str:
    try:
        if rho is not None and np.isfinite(float(rho)) and abs(float(rho)) > 0.999:
            return "t approximation saturates at |ρ|=1"
        p = float(p)
    except (TypeError, ValueError):
        return "undefined"
    if not np.isfinite(p):
        return "undefined"
    return f"{p:.3g}"


def fmt_rho(row) -> str:
    return (
        f"n={int(row['n_mice'])}, ρ={row['rho']:+.3f}, "
        f"exact p={row['p_exact']:.4f} ({int(row['n_extreme'])}/{int(row['n_perm'])})"
    )


def spec_sentence(row) -> str:
    return (
        f"study `{row['study']}`, genotype `{row['genotype']}`, treatment `{row['treatment']}`, "
        f"epithelial gate `{row['gate']}`, min epithelial cells `{int(row['min_epi'])}`, "
        f"Tacstd2% `{row['x']}`, endpoint `{row['y']}`"
    )


def write_finding(path: Path, winners: dict[str, pd.DataFrame], eligible: pd.DataFrame,
                  validation: pd.DataFrame, gene_present: dict[str, bool], summary: dict,
                  max_n: dict, min_p: dict) -> None:
    def block(family: str) -> str:
        top = winners[family].head(8)
        if top.empty:
            return f"No eligible {family} test.\n"
        best = top.iloc[0]
        fam = eligible[eligible.family == family]
        n_sat = int((fam.abs_rho > 0.999).sum())
        n_sat_neg = int(((fam.abs_rho > 0.999) & (fam.rho < 0)).sum())
        n_sat_pos = int(((fam.abs_rho > 0.999) & (fam.rho > 0)).sum())
        lines = [
            f"Maximum |ρ| in {family}: **{fmt_rho(best)}**.",
            "",
            f"Shown specification ({spec_sentence(best)}). "
            f"Dataset-residual Spearman ρ={best['rho_dataset']:+.3f} "
            f"(asymptotic p={fmt_asymp(best['p_dataset_asymp'], best['rho_dataset'])}). "
            f"Genotype-residual Spearman ρ={best['rho_genotype']:+.3f} "
            f"(asymptotic p={fmt_asymp(best['p_genotype_asymp'], best['rho_genotype'])}).",
            "",
            f"Mice: {best['mice']}.",
            "",
            f"Eligible tests in this family: {len(fam)}. "
            f"Tests on the |ρ|=1 ceiling: {n_sat} "
            f"({n_sat_neg} negative, {n_sat_pos} positive). "
            "The row above is the ceiling tie that prefers the locked gate and Tacstd2 count>0. "
            "It is not a larger correlation than the other ceiling rows.",
            "",
        ]
        if family == "T_frac":
            cd8 = fam[fam["y"] == "frac_cd8"]
            if not cd8.empty:
                cd8_best = _rank_frame(cd8).iloc[0]
                lines += [
                    f"CD8 fraction does not reach the ceiling. Its maximum is {fmt_rho(cd8_best)} "
                    f"({spec_sentence(cd8_best)}).",
                    "",
                ]
        full = max_n.get(family)
        if full is not None and int(full["n_mice"]) != int(best["n_mice"]):
            flip = ""
            if (
                np.isfinite(full["rho"]) and np.isfinite(full["rho_dataset"])
                and full["rho"] * full["rho_dataset"] < 0
            ):
                flip = " The dataset residual has the opposite sign."
            lines += [
                f"Largest mouse set in this family ({int(full['n_mice'])} mice): {fmt_rho(full)}. "
                f"{spec_sentence(full)}. Dataset-residual ρ={full['rho_dataset']:+.3f}.{flip}",
                "",
            ]
        mp = min_p.get(family)
        if mp is not None and abs(float(mp["p_exact"]) - float(best["p_exact"])) > 1e-12:
            lines += [
                f"Smallest exact p in this family (not the |ρ| rule): {fmt_rho(mp)}. "
                f"{spec_sentence(mp)}. Dataset-residual ρ={mp['rho_dataset']:+.3f}.",
                "",
            ]
        lines += [
            "Next rows under the |ρ| sort:",
            "",
            "| |ρ| | ρ | exact p | n | studies | genotype filter | treatment | gate | min epi | Tacstd2% | endpoint | dataset residual |",
            "|---:|---:|---:|---:|---:|---|---|---|---:|---|---|---:|",
        ]
        for _, r in top.iterrows():
            lines.append(
                f"| {r['abs_rho']:.3f} | {r['rho']:+.3f} | {r['p_exact']:.4f} | {int(r['n_mice'])} | "
                f"{r['study']} | {r['genotype']} | {r['treatment']} | {r['gate']} | {int(r['min_epi'])} | "
                f"{r['x']} | {r['y']} | {r['rho_dataset']:+.3f} |"
            )
        return "\n".join(lines) + "\n"

    n_t = int((eligible.family == "T_frac").sum())
    n_i = int((eligible.family == "IFN_APM").sum())
    val_ok = bool((validation["n_epi_delta"].abs() <= 2).all() and (validation["frac_delta"].abs() < 1e-6).all())
    text = f"""# Public integrate cohorts: maximum |ρ| for mouse-level Tacstd2%

Public processed counts only: GSE154977 (KP, AT2-lineage FACS), GSE180963 (1 K + 1 KL, mixed), GSE165641 (2 KL, mixed). Private 8 KL matrices were not read and were not merged. Unit is the mouse. Tacstd2 is not an epithelial caller.

This file reports the maximum absolute Spearman correlation on a grid that was fixed in `scripts/analyze.py` before ranking. It is the extreme of that search, not a single prespecified test. Eligible tests need at least 4 mice. T-fraction tests that still contain a GSE154977 FACS library are a design no-go and are absent from the T argmax. Mixed-digest gates that call more than half of a mouse's cells epithelial drop that mouse.

Locked-gate replay against the integrate mouse table checks n_epi and T/NK fraction (same epithelial caller as the Cldn4 integrate). Tacstd2% is the new exposure; Cldn4% is recorded only for gate diagnostics. Epithelial IFN means may differ slightly from the integrate object because this run log-normalizes with each cell's full UMI total (median my/published IFN = {summary["ifn_scale_median"]:.3f}).

Eligible tests searched: T_frac {n_t}, IFN/APM {n_i}, MHC (side, not in the IFN/APM argmax) {int((eligible.family == "MHC_published").sum())}.

Epithelial gates that never produce an eligible T-fraction test, because GSE180963 then has fewer than 10 epithelial cells: {", ".join(summary["t_gates_without_eligible_test"]) or "none"}. Those gates require Sftpc-negative epithelium. The T-fraction ceiling below is carried by gates that still keep Sftpc-positive cells.

## T fraction

{block("T_frac")}

On the locked gate the four mixed mice are ordered the same way inside each study. Within-study Tacstd2% order versus T/NK is reported from the winner mice table. Each study alone has two mice, so a within-study Spearman is not an eligible test. The dataset residual stays −1 because those within-study orders agree. GSE154977 stays out of this endpoint: it is an AT2-lineage sort with essentially no T/NK cells.

## IFN/APM

IFN is the epithelial mean of the integrate ISG list, using only genes present in every library ({", ".join(g for g in IFN_GENES if gene_present.get(g))}). Ifnb1 is absent from GSE180963, so it is not in the mean. APM is the epithelial mean of B2m, H2-K1, H2-D1, Tap1, Tap2, Tapbp, Psmb8, Psmb9, Psmb10, Nlrc5 (present in every library: {", ".join(g for g in APM_GENES if gene_present.get(g))}). IFN_APM is the mean of those two scores. MHC-II genes stay in the side endpoint `MHC_published` and do not enter this argmax.

{block("IFN_APM")}

## MHC side endpoint (not the headline)

{block("MHC_published")}

## What the grid is allowed to change

Study filter, genotype filter, cisplatin drop, epithelial gate, epithelial-cell floor (10, 30, 50), and the Tacstd2 percent definition (count > 0, count ≥ 2, count ≥ 5, log1p CPM-10k ≥ 1). T endpoints are the T/NK call, CD3, and CD8 fractions of QC cells, with T/NK-marker cells that also pass the epithelial gate removed from the numerator.

Cells after the shared QC (nFeature ≥ 200, nCount ≥ 500, percent.mt < 25): {summary["n_cells_qc"]}. Mice: {summary["n_mice"]}. Private 8 KL mice used: 0.
"""
    path.write_text(text)


def scatter(path: Path, mice: pd.DataFrame, row: pd.Series, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    colors = {"GSE154977": "#4C78A8", "GSE165641": "#F58518", "GSE180963": "#54A24B"}
    markers = {"KP": "o", "KL": "s", "K": "D"}
    for _, r in mice.iterrows():
        ax.scatter(
            r[row["x"]], r[row["y"]],
            c=colors.get(r["dataset"], "0.3"),
            marker=markers.get(r["genotype"], "o"),
            s=70, zorder=3,
        )
        ax.annotate(r["mouse"].replace("GSE", ""), (r[row["x"]], r[row["y"]]),
                    textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax.set_xlabel(f"epithelial Tacstd2% ({row['x']})")
    ax.set_ylabel(row["y"])
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def rho_hist(path: Path, eligible: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True)
    for ax, family in zip(axes, ["T_frac", "IFN_APM"]):
        sub = eligible.loc[eligible.family == family, "rho"]
        ax.hist(sub, bins=21, color="#4C78A8", edgecolor="white")
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_title(f"{family} eligible ρ (n tests={len(sub)})")
        ax.set_xlabel("Spearman ρ")
    axes[0].set_ylabel("tests")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def mice_for_row(score_tables, row: pd.Series) -> pd.DataFrame:
    base = score_tables[row["gate"]]
    studies = STUDY_GRID[row["study"]]
    genos = GENO_GRID[row["genotype"]]
    treats = TREAT_GRID[row["treatment"]]
    filt = base[
        base["dataset"].isin(studies)
        & base["genotype"].isin(genos)
        & base["treatment"].isin(treats)
    ]
    kept, status = eligible_mice(filt, int(row["min_epi"]), row["y"])
    if status != "ok":
        raise RuntimeError(f"winner failed eligibility on replay: {status}")
    return kept


def main() -> None:
    _self_check()
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out
    tab = out / "tables"
    fig = out / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    counts, lognorm, meta, gene_index, shared = stack_libraries(args.data)
    if "Tacstd2" not in gene_index:
        raise RuntimeError("Tacstd2 missing after intersect")
    forbidden = [c for c in meta.columns if "8KL" in c or "KLA" in c]
    if forbidden:
        raise RuntimeError(f"unexpected private columns {forbidden}")

    gates = gate_masks(counts, lognorm, meta, gene_index)
    # Score only genes present in every library. Ifnb1 is absent from GSE180963;
    # filling it with zeros would pull that study's IFN mean down.
    ifn_genes = [g for g in IFN_GENES if g in shared]
    apm_genes = [g for g in APM_GENES if g in shared]
    mhc_genes = [g for g in MHC_GENES if g in shared]
    dropped = [g for g in IFN_GENES + APM_GENES + MHC_GENES if g not in shared]
    print("dropped from scores (not in every library):", dropped, flush=True)
    ifn = mean_genes(lognorm, gene_index, ifn_genes)
    apm = mean_genes(lognorm, gene_index, apm_genes)
    mhc = mean_genes(lognorm, gene_index, mhc_genes)
    gene_i = gene_index["Tacstd2"]
    cldn4_i = gene_index.get("Cldn4")
    gene_present = {g: g in shared for g in MARKERS}

    score_tables = {}
    pieces = []
    for name, mask in gates.items():
        print(f"score {name} epi {int(mask.sum())}", flush=True)
        table = mouse_table_for_gate(mask, counts, lognorm, meta, gene_index, gene_i, ifn, apm, mhc)
        if cldn4_i is not None:
            side = mouse_table_for_gate(mask, counts, lognorm, meta, gene_index, cldn4_i, ifn, apm, mhc)
            table["cldn4_pct_gt0"] = side["pct_gt0"].to_numpy()
        # Mean Tacstd2 log1p in epithelium (diagnostic; exposure remains percent).
        means = []
        for mouse in table["mouse"]:
            ix = np.asarray(list(meta.groupby("mouse", sort=True).groups[mouse]))
            epi = ix[mask[ix]]
            means.append(float(np.nanmean(lognorm[gene_i, epi])) if len(epi) else np.nan)
        table["tacstd2_mean"] = means
        table.insert(0, "gate", name)
        score_tables[name] = table
        pieces.append(table)
    scores = pd.concat(pieces, ignore_index=True)
    scores.to_csv(tab / "mouse_scores_by_gate.tsv", sep="\t", index=False)

    validation = validate_locked(score_tables["locked"])
    validation.to_csv(tab / "locked_gate_validation.tsv", sep="\t", index=False)
    print(validation.to_string(index=False), flush=True)
    if validation["n_epi_delta"].abs().max() > 5 or validation["frac_delta"].abs().max() > 1e-4:
        raise RuntimeError("locked gate does not replay the integrate mouse table")

    print("grid", flush=True)
    eligible = run_grid(score_tables)
    eligible.to_csv(tab / "grid_eligible.tsv", sep="\t", index=False)
    print(f"eligible rows {len(eligible)}", flush=True)

    winners = {}
    winner_rows = []
    max_n = {}
    min_p_rows = {}
    for family in ["T_frac", "IFN_APM", "MHC_published"]:
        ranked = pick_winner(eligible, family)
        winners[family] = ranked
        ranked.head(25).to_csv(tab / f"top_{family}.tsv", sep="\t", index=False)
        if ranked.empty:
            continue
        best = ranked.iloc[0]
        winner_rows.append(best)
        mice = mice_for_row(score_tables, best)
        mice.to_csv(tab / f"winner_mice_{family}.tsv", sep="\t", index=False)
        scatter(
            fig / f"winner_{family}.png",
            mice,
            best,
            f"{family} max |ρ|\n{fmt_rho(best)}\n{best['gate']} / {best['x']} / {best['study']} / {best['genotype']}",
        )
        full = row_at_max_n(eligible, family)
        max_n[family] = full
        mp = row_min_p(eligible, family)
        min_p_rows[family] = mp
        for tag, extra in (("maxn", full), ("minp", mp)):
            if extra is None:
                continue
            if int(extra["n_mice"]) == int(best["n_mice"]) and extra["gate"] == best["gate"] and extra["x"] == best["x"] and extra["y"] == best["y"] and extra["study"] == best["study"]:
                continue
            extra_mice = mice_for_row(score_tables, extra)
            extra_mice.to_csv(tab / f"{tag}_mice_{family}.tsv", sep="\t", index=False)
            scatter(
                fig / f"{tag}_{family}.png",
                extra_mice,
                extra,
                f"{family} {tag}\n{fmt_rho(extra)}\n{extra['gate']} / {extra['x']} / {extra['y']}",
            )
    if winner_rows:
        pd.DataFrame(winner_rows).to_csv(tab / "winners.tsv", sep="\t", index=False)
    rho_hist(fig / "eligible_rho_hist.png", eligible)

    # Locked baseline on the published specification, recomputed here.
    base_rows = []
    locked = score_tables["locked"]
    specs = [
        ("T_frac_locked", "drop_GSE154977", "all", "all", "pct_gt0", "frac_tnk", 10),
        ("IFN_locked_all8", "all3", "all", "all", "pct_gt0", "IFN", 10),
        ("APM_locked_all8", "all3", "all", "all", "pct_gt0", "APM", 10),
        ("IFN_APM_locked_all8", "all3", "all", "all", "pct_gt0", "IFN_APM", 10),
    ]
    for label, sname, gname, tname, xname, yname, min_epi in specs:
        filt = locked[
            locked["dataset"].isin(STUDY_GRID[sname])
            & locked["genotype"].isin(GENO_GRID[gname])
            & locked["treatment"].isin(TREAT_GRID[tname])
        ]
        kept, status = eligible_mice(filt, min_epi, yname)
        rec = eval_one(kept, xname, yname) if status == "ok" else None
        base_rows.append(dict(label=label, status=status, study=sname, genotype=gname,
                              treatment=tname, gate="locked", min_epi=min_epi, x=xname, y=yname,
                              **(rec or {})))
    pd.DataFrame(base_rows).to_csv(tab / "locked_baseline.tsv", sep="\t", index=False)

    t_gates = set(eligible.loc[eligible.family == "T_frac", "gate"])
    summary = {
        "cohorts": ["GSE154977", "GSE180963", "GSE165641"],
        "private_8kl_mice": 0,
        "merged_with_private_8kl": False,
        "n_cells_qc": int(len(meta)),
        "n_mice": int(meta["mouse"].nunique()),
        "n_eligible_T_frac": int((eligible.family == "T_frac").sum()),
        "n_eligible_IFN_APM": int((eligible.family == "IFN_APM").sum()),
        "ifn_scale_median": float((validation["IFN"] / validation["IFN_published"]).median()),
        "t_gates_without_eligible_test": sorted(set(score_tables) - t_gates),
        "python": platform.python_version(),
        "scipy": scipy_version,
        "genes_present": {g: gene_present[g] for g in MARKERS},
    }
    (tab / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame({"gene": list(gene_present), "present": list(gene_present.values())}).to_csv(
        tab / "gene_inventory.tsv", sep="\t", index=False
    )
    write_finding(
        out / "FINDING.md", winners, eligible, validation, gene_present, summary, max_n, min_p_rows
    )
    print("wrote", out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
