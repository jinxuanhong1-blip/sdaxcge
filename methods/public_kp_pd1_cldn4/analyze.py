#!/usr/bin/env python3
"""Score public mouse KP / Kras / ICI-sensitive lung tumors: Cldn4 vs T/NK and PD-1 vs control.

Additive public MOUSE. Cldn4-only. No private 8-KL. No mega-merge.
Processed GEO matrices only. Series are scored separately.
"""

from __future__ import annotations

import argparse
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

from gene_sets import CLDN4_ENS, CLDN4_SYMBOL, SCRNA_TNK_GATE, SCRNA_TSCORE, TNK_SYMBOLS

HERE = Path(__file__).resolve().parent
DEFAULT_DL = Path("/tmp/geo_dl")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

GSE114601_META = {
    "s1795": {"treatment": "anti-PD1", "arm": "PD-1"},
    "s2521": {"treatment": "anti-PD1", "arm": "PD-1"},
    "s2596": {"treatment": "Vehicle", "arm": "control"},
    "s2617": {"treatment": "Vehicle", "arm": "control"},
    "s2503": {"treatment": "JQ1", "arm": "JQ1"},
    "s2625": {"treatment": "JQ1", "arm": "JQ1"},
    "s2467": {"treatment": "anti-PD1-JQ1", "arm": "PD-1+JQ1"},
    "s2669": {"treatment": "anti-PD1-JQ1", "arm": "PD-1+JQ1"},
}

# Titles on GEO; each library is 3–4 mice pooled.
GSE157880_META = {
    "0-1_S13": {"title": "IgG_0Gy_1", "treatment": "IgG", "rt_gy": 0, "arm": "control"},
    "0-2_S14": {"title": "IgG_0Gy_2", "treatment": "IgG", "rt_gy": 0, "arm": "control"},
    "0-3_S15": {"title": "IgG_0Gy_3", "treatment": "IgG", "rt_gy": 0, "arm": "control"},
    "4_2_S26": {"title": "IgG_4Gy_1", "treatment": "IgG", "rt_gy": 4, "arm": "control"},
    "4_3_S27": {"title": "IgG_4Gy_2", "treatment": "IgG", "rt_gy": 4, "arm": "control"},
    "8_1_S31": {"title": "IgG_8Gy_1", "treatment": "IgG", "rt_gy": 8, "arm": "control"},
    "8_2_S32": {"title": "IgG_8Gy_2", "treatment": "IgG", "rt_gy": 8, "arm": "control"},
    "8_3_S33": {"title": "IgG_8Gy_3", "treatment": "IgG", "rt_gy": 8, "arm": "control"},
    "0-4_S16": {"title": "PD-1_0Gy_1", "treatment": "PD-1", "rt_gy": 0, "arm": "PD-1"},
    "0-5_S17": {"title": "PD-1_0Gy_2", "treatment": "PD-1", "rt_gy": 0, "arm": "PD-1"},
    "4_4_S28": {"title": "PD-1_4Gy_1", "treatment": "PD-1", "rt_gy": 4, "arm": "PD-1"},
    "4_5_S29": {"title": "PD-1_4Gy_2", "treatment": "PD-1", "rt_gy": 4, "arm": "PD-1"},
    "4_6_S30": {"title": "PD-1_4Gy_3", "treatment": "PD-1", "rt_gy": 4, "arm": "PD-1"},
    "8_4_S34": {"title": "PD-1_8Gy_1", "treatment": "PD-1", "rt_gy": 8, "arm": "PD-1"},
    "8_5_S35": {"title": "PD-1_8Gy_2", "treatment": "PD-1", "rt_gy": 8, "arm": "PD-1"},
    "8_6_S36": {"title": "PD-1_8Gy_3", "treatment": "PD-1", "rt_gy": 8, "arm": "PD-1"},
}

# total viable-cell samples only (not TAM / CD4 sorts).
GSE169194_TOTAL = {
    "X1740R_159_04": {"replicate": 1, "treatment": "IgG", "arm": "control"},
    "X1740R_159_08": {"replicate": 2, "treatment": "IgG", "arm": "control"},
    "X1740R_159_12": {"replicate": 3, "treatment": "IgG", "arm": "control"},
    "X1740R_159_16": {"replicate": 1, "treatment": "A2V", "arm": "A2V"},
    "X1740R_159_20": {"replicate": 2, "treatment": "A2V", "arm": "A2V"},
    "X1740R_159_24": {"replicate": 3, "treatment": "A2V", "arm": "A2V"},
    "X1740R_159_32": {"replicate": 2, "treatment": "A2V+aPD-1", "arm": "A2V+PD-1"},
    "X1740R_159_36": {"replicate": 3, "treatment": "A2V+aPD-1", "arm": "A2V+PD-1"},
    "X1740R_159_40": {"replicate": 4, "treatment": "A2V+aPD-1", "arm": "A2V+PD-1"},
}

GSE246922_KP = {
    "KP_1": {"group": "parental", "arm": "parental", "rep": 1},
    "KP_2": {"group": "parental", "arm": "parental", "rep": 2},
    "KP_3": {"group": "parental", "arm": "parental", "rep": 3},
    "KPy_1": {"group": "chronic_IFNG", "arm": "chronic_IFNG", "rep": 1},
    "KPy_2": {"group": "chronic_IFNG", "arm": "chronic_IFNG", "rep": 2},
    "KPy_3": {"group": "chronic_IFNG", "arm": "chronic_IFNG", "rep": 3},
    "ResKP2_1": {"group": "1st_relapse", "arm": "ICB_relapse", "rep": 1},
    "ResKP2_2": {"group": "1st_relapse", "arm": "ICB_relapse", "rep": 2},
    "ResKP2_3": {"group": "1st_relapse", "arm": "ICB_relapse", "rep": 3},
    "ResResKP_1": {"group": "2nd_relapse", "arm": "ICB_relapse", "rep": 1},
    "ResResKP_2": {"group": "2nd_relapse", "arm": "ICB_relapse", "rep": 2},
    "ResResKP_3": {"group": "2nd_relapse", "arm": "ICB_relapse", "rep": 3},
    "ResKPlate_1": {"group": "late_relapse", "arm": "ICB_relapse", "rep": 1},
    "ResKPlate_2": {"group": "late_relapse", "arm": "ICB_relapse", "rep": 2},
    "ResKPlate_3": {"group": "late_relapse", "arm": "ICB_relapse", "rep": 3},
}

GSE133604_LIBS = [
    ("GSM3912860_1_Ctrl_matrix.mtx.gz", "Ctrl", "control", "KP WT"),
    ("GSM3912861_2_CtrlplusPD1_matrix.mtx.gz", "Ctrl+PD-1", "PD-1", "KP WT"),
    ("GSM3912862_3_ko_matrix.mtx.gz", "Asf1a-KO", "control", "KP Asf1a-KO"),
    ("GSM3912863_4_KOplusPD1_matrix.mtx.gz", "Asf1a-KO+PD-1", "PD-1", "KP Asf1a-KO"),
]


def log2p1(x: pd.Series | np.ndarray) -> pd.Series | np.ndarray:
    return np.log2(np.asarray(x, dtype=float) + 1.0)


def mean_z(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    present = [c for c in cols if c in df.columns]
    if not present:
        return pd.Series(np.nan, index=df.index)
    z = df[present].apply(
        lambda s: (s - s.mean()) / s.std(ddof=0) if float(s.std(ddof=0) or 0) else 0.0,
        axis=0,
    )
    return z.mean(axis=1)


def spearman(x: pd.Series, y: pd.Series) -> dict:
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r = stats.spearmanr(x[mask], y[mask])
    return {"n": n, "rho": float(r.statistic), "p": float(r.pvalue)}


def welch_mwu(ctrl: np.ndarray, treat: np.ndarray) -> dict:
    ctrl = np.asarray(ctrl, dtype=float)
    treat = np.asarray(treat, dtype=float)
    ctrl = ctrl[np.isfinite(ctrl)]
    treat = treat[np.isfinite(treat)]
    out = {
        "n_ctrl": int(ctrl.size),
        "n_treat": int(treat.size),
        "mean_ctrl": float(np.mean(ctrl)) if ctrl.size else np.nan,
        "mean_treat": float(np.mean(treat)) if treat.size else np.nan,
        "delta": np.nan,
        "welch_t": np.nan,
        "welch_p": np.nan,
        "mwu_u": np.nan,
        "mwu_p": np.nan,
    }
    if ctrl.size and treat.size:
        out["delta"] = float(np.mean(treat) - np.mean(ctrl))
    if ctrl.size < 2 or treat.size < 2:
        return out
    t_res = stats.ttest_ind(treat, ctrl, equal_var=False, alternative="two-sided")
    u_res = stats.mannwhitneyu(treat, ctrl, alternative="two-sided")
    out["welch_t"] = float(t_res.statistic)
    out["welch_p"] = float(t_res.pvalue)
    out["mwu_u"] = float(u_res.statistic)
    out["mwu_p"] = float(u_res.pvalue)
    return out


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.2f}" if p >= 0.01 else f"{p:.4f}"


def fmt_num(x: float, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    return f"{x:.{nd}f}"


def load_ens_map(dl: Path) -> dict[str, str]:
    genes = pd.read_csv(dl / "GSE133604_genes.tsv.gz", sep="\t", header=None)
    genes.columns = ["ens", "symbol", "type"][: genes.shape[1]]
    return dict(zip(genes["ens"].astype(str), genes["symbol"].astype(str)))


def score_bulk_matrix(
    expr: pd.DataFrame,
    meta: dict[str, dict],
    series: str,
    unit: str,
    gene_id_kind: str,
    ens_map: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """expr: genes x samples. Returns per-unit table."""
    if gene_id_kind == "ensembl":
        assert ens_map is not None
        symbols = expr.index.astype(str).map(lambda g: ens_map.get(g, g))
        expr = expr.copy()
        expr.index = symbols
        expr = expr.groupby(expr.index).max()
    expr.index = expr.index.astype(str)

    present_tnk = [g for g in TNK_SYMBOLS if g in expr.index]
    if CLDN4_SYMBOL not in expr.index:
        raise RuntimeError(f"{series}: Cldn4 not in features")

    samples = [c for c in expr.columns if c in meta]
    rows = []
    tnk_log = {}
    for sid in samples:
        cldn4 = float(expr.loc[CLDN4_SYMBOL, sid])
        tnk_vals = {g: float(expr.loc[g, sid]) for g in present_tnk}
        tnk_mean = float(np.mean(list(tnk_vals.values()))) if tnk_vals else np.nan
        tnk_log[sid] = float(np.mean(log2p1(np.array(list(tnk_vals.values()))))) if tnk_vals else np.nan
        rec = {
            "series": series,
            "unit_id": sid,
            "unit": unit,
            "cldn4": cldn4,
            "cldn4_log2p1": float(np.log2(cldn4 + 1.0)),
            "tnk_mean": tnk_mean,
            "tnk_mean_log2p1": tnk_log[sid],
            "n_tnk_genes": len(present_tnk),
            **meta[sid],
        }
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["tnk_z"] = mean_z(out, ["tnk_mean_log2p1"])
    out["cldn4_z"] = mean_z(out, ["cldn4_log2p1"])
    return out, present_tnk


def contrast_rows(df: pd.DataFrame, series: str, model: str, ctrl_arm: str, treat_arm: str, label: str) -> list[dict]:
    a = df[df["arm"] == ctrl_arm]
    b = df[df["arm"] == treat_arm]
    rows = []
    for metric, col in [("Cldn4", "cldn4"), ("Cldn4_log2p1", "cldn4_log2p1"), ("T/NK", "tnk_mean_log2p1")]:
        stats_d = welch_mwu(a[col].to_numpy(), b[col].to_numpy())
        rows.append(
            {
                "series": series,
                "model": model,
                "contrast": label,
                "metric": metric,
                **stats_d,
            }
        )
    return rows


def score_gse114601(dl: Path, tables: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    expr = pd.read_csv(dl / "GSE114601_counts.normalized.csv.gz", index_col=0)
    units, genes = score_bulk_matrix(expr, GSE114601_META, "GSE114601", "mouse_nodule", "symbol")
    units.to_csv(tables / "GSE114601_mouse_level.tsv", sep="\t", index=False)
    # primary: Vehicle vs anti-PD1
    prim = units[units["arm"].isin(["control", "PD-1"])]
    sp_all = spearman(units["cldn4_log2p1"], units["tnk_mean_log2p1"])
    sp_prim = spearman(prim["cldn4_log2p1"], prim["tnk_mean_log2p1"])
    contrasts = contrast_rows(units, "GSE114601", "KP GEMM Kras/Trp53 nodules", "control", "PD-1", "anti-PD1 vs Vehicle")
    contrasts += contrast_rows(units, "GSE114601", "KP GEMM Kras/Trp53 nodules", "control", "JQ1", "JQ1 vs Vehicle")
    contrasts += contrast_rows(units, "GSE114601", "KP GEMM Kras/Trp53 nodules", "control", "PD-1+JQ1", "anti-PD1+JQ1 vs Vehicle")
    summary = {
        "series": "GSE114601",
        "model": "KP GEMM (KrasTrp53) lung tumor nodules",
        "unit": "biological mouse nodule (n=2 / arm)",
        "cldn4_present": True,
        "tnk_genes": genes,
        "spearman_all_n8": sp_all,
        "spearman_vehicle_pd1_n4": sp_prim,
        "primary_contrast": "anti-PD1 vs Vehicle",
        "n_ctrl": 2,
        "n_pd1": 2,
    }
    return units, contrasts, summary


def score_gse157880(dl: Path, tables: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    raw = pd.read_csv(dl / "GSE157880_Bulk048.txt.gz", sep="\t")
    gene_col = "Gene Symbol" if "Gene Symbol" in raw.columns else "gene_name"
    expr = raw.set_index(gene_col)
    sample_cols = [c for c in expr.columns if c in GSE157880_META]
    expr = expr[sample_cols]
    expr = expr.groupby(expr.index).max()
    units, genes = score_bulk_matrix(expr, GSE157880_META, "GSE157880", "pooled_library", "symbol")
    units.to_csv(tables / "GSE157880_sample_level.tsv", sep="\t", index=False)
    contrasts = []
    for gy in (0, 4, 8):
        sub = units[units["rt_gy"] == gy]
        contrasts += contrast_rows(
            sub,
            "GSE157880",
            f"HKP1 orthotopic lung, {gy} Gy",
            "control",
            "PD-1",
            f"PD-1 vs IgG at {gy} Gy (pooled 3-4 mice/library)",
        )
    # also all libraries, ignoring RT (not the claim)
    sp0 = spearman(units.loc[units["rt_gy"] == 0, "cldn4_log2p1"], units.loc[units["rt_gy"] == 0, "tnk_mean_log2p1"])
    sp_all = spearman(units["cldn4_log2p1"], units["tnk_mean_log2p1"])
    summary = {
        "series": "GSE157880",
        "model": "HKP1 (Kras/p53) orthotopic lung",
        "unit": "pooled library (3-4 mice each); not one-mouse",
        "cldn4_present": True,
        "tnk_genes": genes,
        "spearman_0Gy": sp0,
        "spearman_all_doses": sp_all,
        "primary_contrast": "PD-1 vs IgG at 0 Gy",
        "n_ctrl_0Gy": int((units.arm.eq("control") & units.rt_gy.eq(0)).sum()),
        "n_pd1_0Gy": int((units.arm.eq("PD-1") & units.rt_gy.eq(0)).sum()),
    }
    return units, contrasts, summary


def score_gse169194(dl: Path, tables: Path, ens_map: dict[str, str]) -> tuple[pd.DataFrame, list[dict], dict]:
    raw = pd.read_csv(dl / "GSE169194_annotated_log2_norm_count.txt.gz", sep="\t")
    raw = raw.rename(columns={raw.columns[0]: "gene"})
    expr = raw.set_index("gene")
    # already log2-norm; do not log again for Cldn4
    if CLDN4_ENS in expr.index:
        expr = expr.rename(index={CLDN4_ENS: CLDN4_SYMBOL})
    units, genes = score_bulk_matrix(expr, GSE169194_TOTAL, "GSE169194", "mouse_total_viable", "ensembl", ens_map)
    # deposited values are already log2; overwrite log2p1 with deposited Cldn4
    units["cldn4_log2p1"] = units["cldn4"]
    units["tnk_mean_log2p1"] = units["tnk_mean"]
    units["tnk_z"] = mean_z(units, ["tnk_mean"])
    units["cldn4_z"] = mean_z(units, ["cldn4"])
    units.to_csv(tables / "GSE169194_total_mouse_level.tsv", sep="\t", index=False)
    contrasts = contrast_rows(
        units,
        "GSE169194",
        "KPM GEMM (Kras/p53/Msh2) total viable tumor cells",
        "control",
        "A2V+PD-1",
        "A2V+aPD-1 vs IgG (no PD-1 monotherapy)",
    )
    contrasts += contrast_rows(
        units,
        "GSE169194",
        "KPM GEMM (Kras/p53/Msh2) total viable tumor cells",
        "control",
        "A2V",
        "A2V vs IgG",
    )
    igg = units[units["arm"] == "control"]
    sp_igg = spearman(igg["cldn4"], igg["tnk_mean"])
    sp_all = spearman(units["cldn4"], units["tnk_mean"])
    summary = {
        "series": "GSE169194",
        "model": "KPM GEMM Kras/p53/Msh2 (ICI-sensitized vs KP)",
        "unit": "total viable cells from one mouse tumor",
        "cldn4_present": True,
        "tnk_genes": genes,
        "spearman_IgG": sp_igg,
        "spearman_all_total": sp_all,
        "primary_contrast": "A2V+aPD-1 vs IgG (combo; no PD-1 mono)",
        "n_igg": 3,
        "n_a2v": 3,
        "n_a2v_pd1": 3,
        "note": "total_A2V_aPD1_1 not deposited",
    }
    return units, contrasts, summary


def score_gse246922(dl: Path, tables: Path, ens_map: dict[str, str]) -> tuple[pd.DataFrame, list[dict], dict]:
    expr = pd.read_csv(dl / "GSE246922_KP_RNA_counts_vst.csv.gz", index_col=0)
    units, genes = score_bulk_matrix(expr, GSE246922_KP, "GSE246922", "CD45neg_tumor", "ensembl", ens_map)
    units["tnk_note"] = "CD45- sort; T/NK genes are leak, not a T/NK fraction"
    units.to_csv(tables / "GSE246922_KP_CD45neg.tsv", sep="\t", index=False)
    # parental vs any ICB relapse (n=3 vs 9) and vs 1st relapse only
    first = units[units["group"].isin(["parental", "1st_relapse"])].copy()
    first["arm"] = first["group"].map({"parental": "control", "1st_relapse": "PD-1"})
    contrasts = contrast_rows(
        first,
        "GSE246922",
        "KP CD45- cancer cells (parental vs 1st ICB relapse)",
        "control",
        "PD-1",
        "1st ICB relapse vs parental (not on-treatment PD-1 vs control)",
    )
    rel = units[units["arm"].isin(["parental", "ICB_relapse"])].copy()
    rel["arm"] = rel["arm"].map({"parental": "control", "ICB_relapse": "PD-1"})
    contrasts += contrast_rows(
        rel,
        "GSE246922",
        "KP CD45- cancer cells (parental vs any ICB relapse)",
        "control",
        "PD-1",
        "any ICB relapse vs parental (not on-treatment)",
    )
    par = units[units["arm"] == "parental"]
    sp = spearman(units["cldn4_log2p1"], units["tnk_mean_log2p1"])
    summary = {
        "series": "GSE246922",
        "model": "KP CD45- cancer cells; parental / IFNG / ICB-relapse",
        "unit": "CD45- sorted tumor (n=3 / group)",
        "cldn4_present": True,
        "tnk_usable": False,
        "tnk_genes_present_but_leak": genes,
        "spearman_all_not_a_TNK_test": sp,
        "n_parental": 3,
        "n_1st_relapse": 3,
        "mean_cldn4_parental": float(par["cldn4"].mean()) if len(par) else np.nan,
    }
    return units, contrasts, summary


def _gene_idx(symbols: pd.Series, names: list[str]) -> dict[str, int]:
    out = {}
    for i, s in enumerate(symbols.astype(str)):
        if s in names and s not in out:
            out[s] = i
    return out


def score_gse133604(dl: Path, tables: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    raw_dir = dl / "GSE133604"
    if not raw_dir.exists():
        raw_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(dl / "GSE133604_RAW.tar") as tar:
            tar.extractall(raw_dir)
    genes = pd.read_csv(dl / "GSE133604_genes.tsv.gz", sep="\t", header=None)
    symbols = genes.iloc[:, 1].astype(str)
    need = [CLDN4_SYMBOL, "Epcam", "Cdh1", "Krt8", *SCRNA_TNK_GATE, *SCRNA_TSCORE]
    gix = _gene_idx(symbols, need)
    if CLDN4_SYMBOL not in gix:
        raise RuntimeError("GSE133604: Cldn4 missing")

    rows = []
    for fname, lib, arm, genotype in GSE133604_LIBS:
        path = raw_dir / fname
        print(f"  MTX {fname}", flush=True)
        mat = mmread(path).tocsr()  # genes x cells
        n_genes, n_cells = mat.shape
        umi = np.asarray(mat.sum(axis=0)).ravel()
        keep = umi >= 200
        mat = mat[:, keep]
        umi = umi[keep]
        n_keep = int(mat.shape[1])

        def vec(sym: str) -> np.ndarray:
            if sym not in gix:
                return np.zeros(n_keep, dtype=float)
            return np.asarray(mat[gix[sym], :].todense()).ravel().astype(float)

        epcam = vec("Epcam")
        cdh1 = vec("Cdh1")
        krt8 = vec("Krt8")
        epi = (epcam > 0) | ((cdh1 > 0) & (krt8 > 0))
        t_any = np.zeros(n_keep, dtype=bool)
        for g in SCRNA_TNK_GATE:
            t_any |= vec(g) > 0
        # Epcam wins over ambient T
        tnk = t_any & ~epi
        cldn4 = vec(CLDN4_SYMBOL)
        tscore_genes = [g for g in SCRNA_TSCORE if g in gix]
        if tscore_genes:
            tscore = np.mean(np.log1p(np.vstack([vec(g) for g in tscore_genes])), axis=0)
        else:
            tscore = np.zeros(n_keep)
        epi_idx = np.where(epi)[0]
        rows.append(
            {
                "series": "GSE133604",
                "unit_id": lib,
                "unit": "10x_library",
                "arm": arm,
                "genotype": genotype,
                "n_cells_deposited": int(n_cells),
                "n_cells_umi200": n_keep,
                "n_epi": int(epi.sum()),
                "n_tnk": int(tnk.sum()),
                "tnk_fraction": float(tnk.sum() / n_keep) if n_keep else np.nan,
                "cldn4": float(np.mean(np.log1p(cldn4[epi_idx]))) if epi_idx.size else np.nan,
                "cldn4_log2p1": float(np.mean(np.log2(cldn4[epi_idx] + 1.0))) if epi_idx.size else np.nan,
                "epi_cldn4_pct_pos": float(np.mean(cldn4[epi_idx] > 0)) if epi_idx.size else np.nan,
                "tnk_mean_log2p1": float(np.mean(tscore)),
                "epi_minus_tnk_cldn4": (
                    float(np.mean(np.log1p(cldn4[epi_idx])) - np.mean(np.log1p(cldn4[tnk])))
                    if epi_idx.size and tnk.sum()
                    else np.nan
                ),
            }
        )
    units = pd.DataFrame(rows)
    units.to_csv(tables / "GSE133604_library_level.tsv", sep="\t", index=False)
    # n=1 vs 1: still emit contrast rows with p=NA
    wt = units[units["genotype"] == "KP WT"]
    ko = units[units["genotype"] == "KP Asf1a-KO"]
    contrasts = contrast_rows(wt, "GSE133604", "KP GEMM lung 10x (WT)", "control", "PD-1", "Ctrl+PD-1 vs Ctrl (n=1 vs 1)")
    contrasts += contrast_rows(
        ko, "GSE133604", "KP GEMM lung 10x (Asf1a-KO)", "control", "PD-1", "KO+PD-1 vs KO (n=1 vs 1)"
    )
    summary = {
        "series": "GSE133604",
        "model": "KP GEMM KrasG12D;Trp53-/- lung 10x",
        "unit": "one 10x library per arm (honest n=1 vs 1)",
        "cldn4_present": True,
        "n_ctrl_wt": 1,
        "n_pd1_wt": 1,
        "direction_only": True,
    }
    return units, contrasts, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dl", type=Path, default=DEFAULT_DL)
    parser.add_argument("--out", type=Path, default=HERE / "tables")
    args = parser.parse_args()
    tables = args.out
    tables.mkdir(parents=True, exist_ok=True)
    dl = args.dl

    ens_map = load_ens_map(dl)
    all_contrasts: list[dict] = []
    summaries: dict[str, dict] = {}
    spearman_rows: list[dict] = []

    print("GSE114601", flush=True)
    u114, c114, s114 = score_gse114601(dl, tables)
    all_contrasts += c114
    summaries["GSE114601"] = s114
    spearman_rows.append({"series": "GSE114601", "subset": "all_8_nodules", **s114["spearman_all_n8"]})
    spearman_rows.append({"series": "GSE114601", "subset": "Vehicle+PD1_n4", **s114["spearman_vehicle_pd1_n4"]})

    print("GSE157880", flush=True)
    u157, c157, s157 = score_gse157880(dl, tables)
    all_contrasts += c157
    summaries["GSE157880"] = s157
    spearman_rows.append({"series": "GSE157880", "subset": "0Gy_only", **s157["spearman_0Gy"]})
    spearman_rows.append({"series": "GSE157880", "subset": "all_doses", **s157["spearman_all_doses"]})

    print("GSE169194", flush=True)
    u169, c169, s169 = score_gse169194(dl, tables, ens_map)
    all_contrasts += c169
    summaries["GSE169194"] = s169
    spearman_rows.append({"series": "GSE169194", "subset": "IgG_total", **s169["spearman_IgG"]})
    spearman_rows.append({"series": "GSE169194", "subset": "all_total", **s169["spearman_all_total"]})

    print("GSE246922", flush=True)
    u246, c246, s246 = score_gse246922(dl, tables, ens_map)
    all_contrasts += c246
    summaries["GSE246922"] = s246
    spearman_rows.append({"series": "GSE246922", "subset": "all_CD45neg_NOT_TNK", **s246["spearman_all_not_a_TNK_test"]})

    print("GSE133604", flush=True)
    u133, c133, s133 = score_gse133604(dl, tables)
    all_contrasts += c133
    summaries["GSE133604"] = s133

    pd.DataFrame(all_contrasts).to_csv(tables / "contrasts.tsv", sep="\t", index=False)
    pd.DataFrame(spearman_rows).to_csv(tables / "spearman.tsv", sep="\t", index=False)

    # combinatorial (per-series, not merged expression)
    combo = []
    for row in all_contrasts:
        if row["metric"] not in ("Cldn4_log2p1", "T/NK"):
            continue
        if "vs Vehicle" in row["contrast"] or "0 Gy" in row["contrast"] or "Ctrl+PD-1 vs Ctrl" in row["contrast"]:
            combo.append(row)
        if row["series"] == "GSE169194" and "A2V+aPD-1 vs IgG" in row["contrast"]:
            combo.append(row)
    pd.DataFrame(combo).to_csv(tables / "combinatorial_pd1_vs_control.tsv", sep="\t", index=False)

    (tables / "summary.json").write_text(json.dumps(summaries, indent=2, default=str), encoding="utf-8")
    print("wrote", tables)
    for s, rec in summaries.items():
        print(s, rec.get("primary_contrast"), "n", rec.get("n_ctrl", rec.get("n_ctrl_0Gy", rec.get("n_igg"))))


if __name__ == "__main__":
    main()
