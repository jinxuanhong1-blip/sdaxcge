#!/usr/bin/env python3
"""GSE207422: TACSTD2 and CLDN4, scored separately, plus dual-high.

Slide 5–6 describe a TROP2-high group with a higher tumor fraction and lower
CD8 T/NK infiltration, and lower TROP2 in MPR than in NMPR. This script
re-downloads nothing itself; it reads the GEO UMI matrix and asks the same
three questions of TACSTD2, of CLDN4, and of the dual-high state.

CLDN4 is not gated on TACSTD2. A flat CLDN4 result is reported as flat.

Primary rules (fixed before looking at p-values):
- Unit: one sample per patient. Primary cohort = post-treatment samples with
  at least 10 A3-malignant cells.
- A3-malignant: epithelial-lineage call, and zero UMI for
  SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3. TACSTD2 and CLDN4 are never markers.
- Gene level: mean log1p(CP10k) inside A3-malignant cells.
- High: at or above the median of that gene in the eligible cohort.
- Dual-high: at or above the median for both genes. Compared with everyone else.
- Tumor fraction: A3-malignant cells / all cells in the sample.
- CD8/NK fraction: (CD3+ CD8A+ cells, or CD3-negative NK-marker cells) / all cells.
- MPR: pathologic response with pCR folded into MPR, as in Hu et al.
- Tests: two-sided exact Mann-Whitney when SciPy allows it; Spearman for the
  continuous gene; two-sided Fisher exact for high vs MPR/NMPR.

Whole-sample (all-cell) averages are reported only as a labeled companion.
Both genes are epithelial-restricted, so an all-cell average moves with tumor
fraction and is not an independent test of that endpoint.
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
from scipy import stats

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
# A3 rule used by the prior CLDN4-only GSE207422 run. Not TACSTD2. Not CLDN4.
NORMAL_LUNG_ZERO = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
GENES_OF_INTEREST = ["TACSTD2", "CLDN4"]
EXTRA = ["PTPRC", "CD8A", "CD8B", "CD4", "CD3D", "CD3E", "CD3G", "NCR1", "KLRF1", "GNLY", "NKG7"]
MIN_MALIGNANT = 10

ENDPOINTS = [
    ("frac_tumor", "A3-malignant fraction"),
    ("frac_cd8nk", "CD8 T + NK fraction"),
    ("frac_tnk", "all T + NK fraction"),
]


def marker_universe() -> set[str]:
    genes = set(GENES_OF_INTEREST) | set(NORMAL_LUNG_ZERO) | set(EXTRA)
    for group in LINEAGE_MARKERS.values():
        genes.update(group)
    return genes


def stream_matrix(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        n_umi = np.zeros(n_cells, dtype=np.float64)
        n_genes = np.zeros(n_cells, dtype=np.int32)
        n_rows = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            n_umi += arr
            n_genes += arr > 0
            if gene in wanted:
                found[gene] = arr
            n_rows += 1
            if n_rows % 4000 == 0:
                print(f"  streamed {n_rows} genes, stored {len(found)}", flush=True)
    print(f"  streamed {n_rows} genes, stored {len(found)}", flush=True)
    return cell_ids, found, n_umi, n_genes, n_rows


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n_cells: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n_cells, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0).astype(np.float32)


def assign_lineage(scores: dict[str, np.ndarray]) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[name] for name in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    second = np.partition(mat, -2, axis=0)[-2]
    labels = np.array(names, dtype=object)[best].copy()
    labels[(top < 0.12) | ((top - second) < 0.02)] = "Unassigned"
    return labels


def sample_of(barcode: str, sample_ids: list[str]) -> str:
    for sample_id in sample_ids:
        if barcode.startswith(sample_id + "_"):
            return sample_id
    return barcode.rsplit("_", 1)[0]


def _finite(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def mannwhitney(a, b) -> dict:
    a = _finite(a)
    b = _finite(b)
    rec = {
        "n_high": int(len(a)),
        "n_low": int(len(b)),
        "median_high": float(np.median(a)) if len(a) else None,
        "median_low": float(np.median(b)) if len(b) else None,
        "U": None,
        "p": None,
        "method": "skipped",
    }
    if len(a) < 1 or len(b) < 1:
        return rec
    try:
        result = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        method = "exact"
    except ValueError:
        result = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
        method = "asymptotic"
    rec.update(U=float(result.statistic), p=float(result.pvalue), method=method)
    return rec


def spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return {"n": int(len(x)), "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(len(x)), "rho": float(rho), "p": float(p)}


def fisher_response(high_mpr: np.ndarray, low_mpr: np.ndarray) -> dict:
    """high_mpr / low_mpr are boolean arrays: True = MPR (pCR folded in)."""
    table = np.array(
        [
            [int(np.sum(high_mpr)), int(np.sum(~high_mpr))],
            [int(np.sum(low_mpr)), int(np.sum(~low_mpr))],
        ],
        dtype=int,
    )
    odds, p = stats.fisher_exact(table, alternative="two-sided")
    return {
        "table_high_mpr": int(table[0, 0]),
        "table_high_nmpr": int(table[0, 1]),
        "table_low_mpr": int(table[1, 0]),
        "table_low_nmpr": int(table[1, 1]),
        "odds_ratio_high_vs_mpr": float(odds),
        "p": float(p),
    }


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{digits}f}"


def median_split(values: pd.Series) -> tuple[pd.Series, float]:
    med = float(np.nanmedian(values.to_numpy(dtype=float)))
    high = values >= med
    return high, med


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_excel(path)
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["Sample"] = meta["Sample"].astype(str).str.strip()
    meta["Patient"] = meta["Patient"].astype(str).str.strip()
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta["is_pre"] = meta["Resource"].astype(str).str.contains("Pre", case=False)
    meta["response_paper"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["residual_tumor"] = pd.to_numeric(meta["Residual Tumor"], errors="coerce")
    return meta


def build_cells(cell_ids, expr, n_umi, n_genes, meta: pd.DataFrame) -> pd.DataFrame:
    n = len(cell_ids)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0).astype(np.float32)
    log_cp = {gene: np.log1p(expr[gene] * scale).astype(np.float32) for gene in expr}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    lineage = assign_lineage(scores)
    sample_ids = sorted(meta["Sample"].unique(), key=len, reverse=True)
    samples = np.array([sample_of(bc, sample_ids) for bc in cell_ids])

    def umi(gene: str) -> np.ndarray:
        if gene not in expr:
            return np.zeros(n, dtype=np.float32)
        return expr[gene]

    normal = np.zeros(n, dtype=np.float32)
    for gene in NORMAL_LUNG_ZERO:
        normal += umi(gene)
    cd3 = (umi("CD3D") > 0) | (umi("CD3E") > 0) | (umi("CD3G") > 0)
    cd8 = cd3 & (umi("CD8A") > 0)
    nk = (~cd3) & ((umi("NCR1") > 0) | (umi("KLRF1") > 0) | ((umi("GNLY") > 0) & (umi("NKG7") > 0)))
    epithelial = lineage == "Epithelial"
    malignant = epithelial & (normal == 0)
    # Sensitivity: drop PTPRC/CD3 doublets from the malignant call.
    malignant_strict = malignant & (umi("PTPRC") == 0) & (~cd3)

    frame = pd.DataFrame(
        {
            "barcode": cell_ids,
            "Sample": samples,
            "n_umi": n_umi,
            "n_genes": n_genes,
            "lineage": lineage,
            "is_epithelial": epithelial,
            "is_malignant": malignant,
            "is_malignant_strict": malignant_strict,
            "is_cd8": cd8,
            "is_nk": nk,
            "is_cd8nk": cd8 | nk,
            "is_tnk_lineage": np.isin(lineage, ["T", "NK"]),
            "TACSTD2": umi("TACSTD2"),
            "CLDN4": umi("CLDN4"),
            "tacstd2_log": log_cp.get("TACSTD2", np.zeros(n, dtype=np.float32)),
            "cldn4_log": log_cp.get("CLDN4", np.zeros(n, dtype=np.float32)),
            "PTPRC": umi("PTPRC"),
            "EPCAM": umi("EPCAM"),
        }
    )
    return frame.merge(meta, on="Sample", how="left")


def _mean(series: pd.Series) -> float:
    if len(series) == 0:
        return np.nan
    return float(series.mean())


def _pct(series: pd.Series) -> float:
    if len(series) == 0:
        return np.nan
    return float(100.0 * (series > 0).mean())


def per_sample_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, sdf in cells.groupby("Sample", sort=True):
        mal = sdf[sdf["is_malignant"]]
        epi = sdf[sdf["is_epithelial"]]
        strict = sdf[sdf["is_malignant_strict"]]
        n = len(sdf)
        row = {
            "Sample": sample,
            "Patient": sdf["Patient"].iloc[0],
            "Resource": sdf["Resource"].iloc[0],
            "is_post": bool(sdf["is_post"].iloc[0]),
            "is_pre": bool(sdf["is_pre"].iloc[0]),
            "Pathology": sdf["Pathology"].iloc[0],
            "Pathologic_Response": sdf["Pathologic Response"].iloc[0],
            "response_paper": sdf["response_paper"].iloc[0],
            "residual_tumor": sdf["residual_tumor"].iloc[0],
            "RECIST": sdf["RECIST"].iloc[0],
            "n_cells": int(n),
            "n_epithelial": int(len(epi)),
            "n_malignant": int(len(mal)),
            "n_malignant_strict": int(len(strict)),
            "n_cd8": int(sdf["is_cd8"].sum()),
            "n_nk": int(sdf["is_nk"].sum()),
            "n_cd8nk": int(sdf["is_cd8nk"].sum()),
            "n_tnk_lineage": int(sdf["is_tnk_lineage"].sum()),
            "frac_tumor": float(len(mal) / n),
            "frac_epithelial": float(len(epi) / n),
            "frac_tumor_strict": float(len(strict) / n),
            "frac_cd8nk": float(sdf["is_cd8nk"].mean()),
            "frac_tnk": float(sdf["is_tnk_lineage"].mean()),
            "tacstd2_mal_mean": _mean(mal["tacstd2_log"]),
            "cldn4_mal_mean": _mean(mal["cldn4_log"]),
            "tacstd2_mal_pct": _pct(mal["TACSTD2"]),
            "cldn4_mal_pct": _pct(mal["CLDN4"]),
            "tacstd2_epi_mean": _mean(epi["tacstd2_log"]),
            "cldn4_epi_mean": _mean(epi["cldn4_log"]),
            "tacstd2_strict_mean": _mean(strict["tacstd2_log"]),
            "cldn4_strict_mean": _mean(strict["cldn4_log"]),
            "tacstd2_all_mean": _mean(sdf["tacstd2_log"]),
            "cldn4_all_mean": _mean(sdf["cldn4_log"]),
            "tacstd2_mal_pct_pos_vs_cd8nk": _pct(sdf.loc[sdf["is_cd8nk"], "TACSTD2"]),
            "cldn4_cd8nk_pct": _pct(sdf.loc[sdf["is_cd8nk"], "CLDN4"]),
            "tacstd2_cd8nk_pct": _pct(sdf.loc[sdf["is_cd8nk"], "TACSTD2"]),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values("Patient").reset_index(drop=True)


def annotate_splits(df: pd.DataFrame, tac_col: str, cldn_col: str) -> pd.DataFrame:
    out = df.copy()
    tac_high, tac_med = median_split(out[tac_col])
    cldn_high, cldn_med = median_split(out[cldn_col])
    out["tacstd2_group"] = np.where(tac_high, "high", "low")
    out["cldn4_group"] = np.where(cldn_high, "high", "low")
    out["dual_high"] = tac_high & cldn_high
    out["tacstd2_median_cut"] = tac_med
    out["cldn4_median_cut"] = cldn_med
    return out


def group_rows(df: pd.DataFrame, group_col: str, high_label: str, gene: str, cohort: str, expr_def: str) -> list[dict]:
    high = df[df[group_col] == high_label] if group_col != "dual_high" else df[df["dual_high"]]
    low = df[df[group_col] != high_label] if group_col != "dual_high" else df[~df["dual_high"]]
    if group_col != "dual_high":
        high = df[df[group_col] == "high"]
        low = df[df[group_col] == "low"]
    rows = []
    for col, label in ENDPOINTS:
        rec = mannwhitney(high[col], low[col])
        rec.update(
            {
                "gene": gene,
                "cohort": cohort,
                "expression": expr_def,
                "test": "high_vs_low_composition",
                "endpoint": label,
                "endpoint_col": col,
                "delta_median_high_minus_low": (
                    None
                    if rec["median_high"] is None or rec["median_low"] is None
                    else float(rec["median_high"] - rec["median_low"])
                ),
            }
        )
        rows.append(rec)
    mpr_df = df[df["response_paper"].isin(["MPR", "NMPR"])]
    high_m = mpr_df[mpr_df[group_col] == "high"] if group_col != "dual_high" else mpr_df[mpr_df["dual_high"]]
    low_m = mpr_df[mpr_df[group_col] == "low"] if group_col != "dual_high" else mpr_df[~mpr_df["dual_high"]]
    if group_col == "dual_high":
        high_m = mpr_df[mpr_df["dual_high"]]
        low_m = mpr_df[~mpr_df["dual_high"]]
    else:
        high_m = mpr_df[mpr_df[group_col] == "high"]
        low_m = mpr_df[mpr_df[group_col] == "low"]
    if len(high_m) and len(low_m):
        fish = fisher_response(
            (high_m["response_paper"] == "MPR").to_numpy(),
            (low_m["response_paper"] == "MPR").to_numpy(),
        )
        rows.append(
            {
                "gene": gene,
                "cohort": cohort,
                "expression": expr_def,
                "test": "high_vs_low_MPR",
                "endpoint": "MPR among high vs low",
                "endpoint_col": "response_paper",
                "n_high": int(len(high_m)),
                "n_low": int(len(low_m)),
                "median_high": fish["table_high_mpr"] / len(high_m),
                "median_low": fish["table_low_mpr"] / len(low_m),
                "delta_median_high_minus_low": fish["table_high_mpr"] / len(high_m) - fish["table_low_mpr"] / len(low_m),
                "U": fish["odds_ratio_high_vs_mpr"],
                "p": fish["p"],
                "method": "fisher_exact",
                "n_high_mpr": fish["table_high_mpr"],
                "n_high_nmpr": fish["table_high_nmpr"],
                "n_low_mpr": fish["table_low_mpr"],
                "n_low_nmpr": fish["table_low_nmpr"],
            }
        )
    return rows


def expression_mpr_rows(df: pd.DataFrame, gene: str, col: str, cohort: str, expr_def: str) -> list[dict]:
    use = df[df["response_paper"].isin(["MPR", "NMPR"])].copy()
    nmpr = use.loc[use["response_paper"] == "NMPR", col]
    mpr = use.loc[use["response_paper"] == "MPR", col]
    # mannwhitney() labels the first array "high". Here group A is NMPR.
    rec = mannwhitney(nmpr, mpr)
    rec.update(
        {
            "gene": gene,
            "cohort": cohort,
            "expression": expr_def,
            "test": "NMPR_vs_MPR_expression",
            "endpoint": "malignant/epithelial mean log1p(CP10k)",
            "endpoint_col": col,
            "n_nmpr": rec.pop("n_high"),
            "n_mpr": rec.pop("n_low"),
            "median_nmpr": rec.pop("median_high"),
            "median_mpr": rec.pop("median_low"),
        }
    )
    rec["delta_median_nmpr_minus_mpr"] = (
        None
        if rec["median_nmpr"] is None or rec["median_mpr"] is None
        else float(rec["median_nmpr"] - rec["median_mpr"])
    )
    return [rec]


def continuous_rows(df: pd.DataFrame, gene: str, expr_col: str, cohort: str, expr_def: str) -> list[dict]:
    rows = []
    for col, label in ENDPOINTS:
        rec = spearman(df[expr_col], df[col])
        rows.append(
            {
                "gene": gene,
                "cohort": cohort,
                "expression": expr_def,
                "test": "spearman",
                "endpoint": label,
                "endpoint_col": col,
                "n": rec["n"],
                "rho": rec["rho"],
                "p": rec["p"],
            }
        )
    use = df[df["response_paper"].isin(["MPR", "NMPR"])]
    rec = spearman(use[expr_col], use["residual_tumor"])
    rows.append(
        {
            "gene": gene,
            "cohort": cohort,
            "expression": expr_def,
            "test": "spearman_residual_tumor_companion",
            "endpoint": "residual tumor fraction",
            "endpoint_col": "residual_tumor",
            "n": rec["n"],
            "rho": rec["rho"],
            "p": rec["p"],
        }
    )
    return rows


def run_block(
    df: pd.DataFrame,
    cohort: str,
    expr_def: str,
    tac_col: str,
    cldn_col: str,
    tumor_col: str,
    cd8_col: str,
    tnk_col: str,
    tumor_label: str = "A3-malignant fraction",
) -> tuple[pd.DataFrame, list[dict], list[dict]]:
    """Score one cohort. Tumor/immune columns may be swapped for a sensitivity."""
    use = df.dropna(subset=[tac_col, cldn_col]).copy()
    use["frac_tumor"] = use[tumor_col]
    use["frac_cd8nk"] = use[cd8_col]
    use["frac_tnk"] = use[tnk_col]
    use["expr_tac"] = use[tac_col]
    use["expr_cldn"] = use[cldn_col]
    use = annotate_splits(use, tac_col, cldn_col)
    tests: list[dict] = []
    tests += group_rows(use, "tacstd2_group", "high", "TACSTD2", cohort, expr_def)
    tests += group_rows(use, "cldn4_group", "high", "CLDN4", cohort, expr_def)
    tests += group_rows(use, "dual_high", "high", "dual-high", cohort, expr_def)
    tests += expression_mpr_rows(use, "TACSTD2", tac_col, cohort, expr_def)
    tests += expression_mpr_rows(use, "CLDN4", cldn_col, cohort, expr_def)
    cont = []
    cont += continuous_rows(use, "TACSTD2", tac_col, cohort, expr_def)
    cont += continuous_rows(use, "CLDN4", cldn_col, cohort, expr_def)
    for bucket in (tests, cont):
        for rec in bucket:
            if rec.get("endpoint") == "A3-malignant fraction":
                rec["endpoint"] = tumor_label
    return use, tests, cont


def balance_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, col in [("TACSTD2", "tacstd2_group"), ("CLDN4", "cldn4_group"), ("dual-high", "dual_high")]:
        for key, sub in df.groupby(col):
            label = "high" if key is True or key == "high" else ("low" if key is False or key == "low" else str(key))
            if col == "dual_high":
                label = "dual-high" if key else "not dual-high"
            rows.append(
                {
                    "gene": gene,
                    "group": label,
                    "n": int(len(sub)),
                    "n_MPR": int((sub["response_paper"] == "MPR").sum()),
                    "n_NMPR": int((sub["response_paper"] == "NMPR").sum()),
                    "n_Adeno": int((sub["Pathology"] == "Adeno").sum()),
                    "n_Squamous": int((sub["Pathology"] == "Squamous").sum()),
                    "n_pre": int(sub["is_pre"].sum()) if "is_pre" in sub else 0,
                    "patients": ",".join(sub["Patient"].astype(str)),
                }
            )
    return pd.DataFrame(rows)


def restriction_table(cells: pd.DataFrame) -> pd.DataFrame:
    """Per-sample percent positive, then the median across post-treatment samples."""
    post_samples = cells.loc[cells["is_post"], "Sample"].unique()
    rows = []
    for gene, umi_col in [("TACSTD2", "TACSTD2"), ("CLDN4", "CLDN4")]:
        for compartment, mask_col in [
            ("A3-malignant", "is_malignant"),
            ("epithelial", "is_epithelial"),
            ("CD8/NK", "is_cd8nk"),
            ("T/NK lineage", "is_tnk_lineage"),
        ]:
            pcts = []
            for sample in post_samples:
                sub = cells[(cells["Sample"] == sample) & cells[mask_col]]
                if len(sub) < 5:
                    continue
                pcts.append(100.0 * (sub[umi_col] > 0).mean())
            rows.append(
                {
                    "gene": gene,
                    "compartment": compartment,
                    "n_post_samples_ge5cells": int(len(pcts)),
                    "median_pct_pos": float(np.median(pcts)) if pcts else None,
                }
            )
    return pd.DataFrame(rows)


def _scatter_panel(ax, df: pd.DataFrame, xcol: str, ycol: str, title: str, xlabel: str, ylabel: str) -> None:
    colors = {"MPR": "#2166ac", "NMPR": "#b2182b", "NE": "#888888"}
    for _, row in df.iterrows():
        ax.scatter(row[xcol], row[ycol], c=colors.get(row["response_paper"], "#555555"), s=42, zorder=3)
        ax.annotate(row["Patient"], (row[xcol], row[ycol]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    rec = spearman(df[xcol], df[ycol])
    ax.set_title(f"{title}\nSpearman ρ={fmt_num(rec['rho'], 2)}  P={fmt_p(rec['p'])}  n={rec['n']}", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def _group_panel(ax, df: pd.DataFrame, group_col: str, ycol: str, title: str, ylabel: str, dual: bool = False) -> None:
    colors = {"MPR": "#2166ac", "NMPR": "#b2182b", "NE": "#888888"}
    if dual:
        groups = [(~df["dual_high"], "not dual"), (df["dual_high"], "dual-high")]
    else:
        groups = [(df[group_col] == "low", "low"), (df[group_col] == "high", "high")]
    rng = np.random.default_rng(1)
    for i, (mask, name) in enumerate(groups):
        sub = df[mask]
        y = sub[ycol].to_numpy(dtype=float)
        if len(y) == 0:
            continue
        x = rng.normal(i, 0.04, size=len(y))
        ax.scatter(
            x,
            y,
            c=[colors.get(v, "#555") for v in sub["response_paper"]],
            s=36,
            zorder=3,
        )
        ax.hlines(np.median(y), i - 0.22, i + 0.22, color="black", lw=1.6, zorder=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[1] for g in groups])
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel)


def _save_composition(df: pd.DataFrame, outdir: Path, stem: str, tumor_ylabel: str, suptitle: str) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), constrained_layout=True)
    specs = [
        (0, "tacstd2_group", False, "TACSTD2"),
        (1, "cldn4_group", False, "CLDN4"),
        (2, "dual_high", True, "dual-high"),
    ]
    for col, group_col, dual, name in specs:
        _group_panel(axes[0, col], df, group_col, "frac_tumor", f"{name}: tumor fraction", tumor_ylabel, dual)
        _group_panel(axes[1, col], df, group_col, "frac_cd8nk", f"{name}: CD8 T + NK", "CD8 T + NK fraction", dual)
    fig.suptitle(suptitle, fontsize=11)
    fig.savefig(outdir / f"{stem}.png", dpi=180)
    fig.savefig(outdir / f"{stem}.pdf")
    plt.close(fig)


def _save_mpr(df: pd.DataFrame, outdir: Path, stem: str, tac_col: str, cldn_col: str, ylabel: str, suptitle: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.2), constrained_layout=True)
    order = ["MPR", "NMPR"]
    colors = {"MPR": "#2166ac", "NMPR": "#b2182b"}
    rng = np.random.default_rng(2)
    for ax, col, name in [(axes[0], tac_col, "TACSTD2"), (axes[1], cldn_col, "CLDN4")]:
        for i, grp in enumerate(order):
            sub = df[df["response_paper"] == grp]
            y = sub[col].to_numpy(dtype=float)
            x = rng.normal(i, 0.04, size=len(y))
            ax.scatter(x, y, c=colors[grp], s=42, zorder=3)
            for _, row in sub.iterrows():
                ax.annotate(row["Patient"], (i, row[col]), fontsize=7, xytext=(4, 0), textcoords="offset points")
            if len(y):
                ax.hlines(np.median(y), i - 0.22, i + 0.22, color="black", lw=1.6)
        nmpr = df.loc[df["response_paper"] == "NMPR", col]
        mpr = df.loc[df["response_paper"] == "MPR", col]
        rec = mannwhitney(nmpr, mpr)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"MPR\nn={len(mpr)}", f"NMPR\nn={len(nmpr)}"])
        ax.set_ylabel(ylabel)
        ax.set_title(f"{name}: NMPR vs MPR\nexact MW P={fmt_p(rec['p'])}", fontsize=10)
    fig.suptitle(suptitle, fontsize=11)
    fig.savefig(outdir / f"{stem}.png", dpi=180)
    fig.savefig(outdir / f"{stem}.pdf")
    plt.close(fig)


def _save_scatter(df: pd.DataFrame, outdir: Path, stem: str, tac_col: str, cldn_col: str, tumor_label: str, expr_label: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2), constrained_layout=True)
    _scatter_panel(axes[0, 0], df, tac_col, "frac_tumor", f"TACSTD2 vs {tumor_label}", expr_label, tumor_label)
    _scatter_panel(axes[0, 1], df, tac_col, "frac_cd8nk", "TACSTD2 vs CD8/NK", expr_label, "CD8 T + NK fraction")
    _scatter_panel(axes[1, 0], df, cldn_col, "frac_tumor", f"CLDN4 vs {tumor_label}", expr_label.replace("TACSTD2", "CLDN4"), tumor_label)
    _scatter_panel(axes[1, 1], df, cldn_col, "frac_cd8nk", "CLDN4 vs CD8/NK", expr_label.replace("TACSTD2", "CLDN4"), "CD8 T + NK fraction")
    fig.savefig(outdir / f"{stem}.png", dpi=180)
    fig.savefig(outdir / f"{stem}.pdf")
    plt.close(fig)


def make_figures(epi: pd.DataFrame, a3: pd.DataFrame, outdir: Path) -> None:
    _save_composition(
        epi, outdir, "fig_epithelial_high_low",
        "epithelial fraction",
        "All 12 post-treatment samples. Expression = epithelial mean. Blue = MPR (pCR included). Red = NMPR.",
    )
    _save_mpr(
        epi, outdir, "fig_epithelial_mpr",
        "expr_tac", "expr_cldn",
        "mean log1p(CP10k) in epithelial cells",
        "All 12 post-treatment samples. pCR is counted as MPR.",
    )
    _save_scatter(
        epi, outdir, "fig_epithelial_spearman",
        "expr_tac", "expr_cldn", "epithelial fraction", "epithelial mean log1p(CP10k)",
    )
    _save_composition(
        a3, outdir, "fig_a3_high_low",
        "A3-malignant fraction",
        "A3 subset only (samples with ≥10 A3-malignant cells). Three of four MPR samples are absent.",
    )
    _save_mpr(
        a3, outdir, "fig_a3_mpr",
        "expr_tac", "expr_cldn",
        "mean log1p(CP10k) in A3-malignant cells",
        "A3 subset. Only one MPR sample still has ≥10 A3-malignant cells.",
    )


def lookup(tests: pd.DataFrame, **kwargs) -> pd.Series | None:
    hit = tests
    for key, value in kwargs.items():
        hit = hit[hit[key] == value]
    if hit.empty:
        return None
    return hit.iloc[0]


def prose_composition(tests: pd.DataFrame, gene: str, cohort: str, expr: str) -> str:
    bits = []
    for endpoint in [
        "A3-malignant fraction",
        "epithelial fraction",
        "PTPRC-negative malignant fraction",
        "CD8 T + NK fraction",
        "all T + NK fraction",
    ]:
        row = lookup(tests, gene=gene, cohort=cohort, expression=expr, test="high_vs_low_composition", endpoint=endpoint)
        if row is None:
            continue
        bits.append(
            f"{endpoint}: high median {fmt_num(row['median_high'])} (n={int(row['n_high'])}) vs "
            f"low median {fmt_num(row['median_low'])} (n={int(row['n_low'])}), "
            f"Δ={fmt_num(row['delta_median_high_minus_low'])}, {row['method']} P={fmt_p(row['p'])}"
        )
    mpr = lookup(tests, gene=gene, cohort=cohort, expression=expr, test="high_vs_low_MPR")
    if mpr is not None:
        bits.append(
            f"MPR rate: high {int(mpr['n_high_mpr'])}/{int(mpr['n_high'])} vs "
            f"low {int(mpr['n_low_mpr'])}/{int(mpr['n_low'])}, Fisher P={fmt_p(mpr['p'])}"
        )
    return "; ".join(bits)


def prose_mpr_expr(tests: pd.DataFrame, gene: str, cohort: str, expr: str) -> str:
    row = lookup(tests, gene=gene, cohort=cohort, expression=expr, test="NMPR_vs_MPR_expression")
    if row is None:
        return "not tested"
    return (
        f"NMPR median {fmt_num(row['median_nmpr'])} (n={int(row['n_nmpr'])}) vs "
        f"MPR median {fmt_num(row['median_mpr'])} (n={int(row['n_mpr'])}), "
        f"Δ(NMPR−MPR)={fmt_num(row['delta_median_nmpr_minus_mpr'])}, {row['method']} P={fmt_p(row['p'])}"
    )


def prose_spearman(cont: pd.DataFrame, gene: str, cohort: str, expr: str) -> str:
    bits = []
    sub = cont[(cont["gene"] == gene) & (cont["cohort"] == cohort) & (cont["expression"] == expr)]
    for _, row in sub.iterrows():
        bits.append(f"{row['endpoint']} ρ={fmt_num(row['rho'], 2)} P={fmt_p(row['p'])} (n={int(row['n'])})")
    return "; ".join(bits)


def _membership_lines(balance: pd.DataFrame) -> list[str]:
    lines = []
    for _, row in balance.iterrows():
        lines.append(
            f"- {row['gene']} {row['group']}: n={int(row['n'])}, MPR={int(row['n_MPR'])}, NMPR={int(row['n_NMPR'])}, "
            f"Adeno={int(row['n_Adeno'])}, Squamous={int(row['n_Squamous'])} ({row['patients']})"
        )
    return lines


def write_finding(
    path: Path,
    summary: dict,
    tests: pd.DataFrame,
    cont: pd.DataFrame,
    balance_epi: pd.DataFrame,
    balance_a3: pd.DataFrame,
    restrict: pd.DataFrame,
) -> None:
    epi_cohort = "post_epithelial_n12"
    epi_expr = "epithelial_mean_log1p_cp10k"
    a3_cohort = "post_ge10_A3"
    a3_expr = "A3_malignant_mean_log1p_cp10k"
    lines = [
        "# GSE207422: TACSTD2 and CLDN4, scored separately",
        "",
        "Public neoadjuvant PD-1 plus chemotherapy NSCLC scRNA-seq (Hu et al., Genome Medicine 2023, PMID 36869384).",
        "Slide 5–6 use this cohort for a TROP2-high state: higher tumor fraction (stated P=0.047), lower CD8 T/NK (stated P=0.015), and lower TROP2 in MPR than in NMPR.",
        "This note scores TACSTD2, CLDN4, and dual-high with the same rules. CLDN4 is not required to match TACSTD2.",
        "",
        "## Data",
        "",
        f"- GEO processed UMI matrix re-downloaded ({summary['matrix_bytes']:,} bytes, sha256 `{summary['matrix_sha256']}`).",
        f"- {summary['n_genes_matrix']:,} genes × {summary['n_cells']:,} cells. The slide text says 90,652 cells; the deposited matrix has {summary['n_cells']:,}, which matches the paper. No cells were dropped to force the slide count.",
        f"- Median genes/cell = {summary['median_genes_per_cell']:.0f} (paper 1,256).",
        "- Sample metadata: 15 patients, one sample each. 12 post-treatment surgery samples (MPR n=4, including pCR P06; NMPR n=8) and 3 pre-treatment biopsies.",
        "- GEO does not deposit author barcode labels or CopyKAT calls. Lineages are reconstructed from canonical markers. TACSTD2 and CLDN4 are not used to call a cell.",
        "- Epithelial: highest lineage score is epithelial. A3-malignant: epithelial and zero UMI for SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3.",
        f"- A3 leaves fewer than 10 malignant cells in {len(summary['dropped_post'])} of 12 post-treatment samples ({', '.join(summary['dropped_post'])}), including 3 of 4 MPR samples (P06 pCR, P11, P14). An MPR contrast on the A3 subset is 1 MPR vs 6 NMPR. The cohort that still contains every post-treatment sample is the epithelial mean (n=12).",
        "",
        "## Where each gene is expressed (post-treatment)",
        "",
        "Median percent of cells with UMI > 0, across post-treatment samples that have at least 5 cells in the compartment.",
        "",
    ]
    for _, row in restrict.iterrows():
        lines.append(
            f"- {row['gene']} in {row['compartment']}: {fmt_num(row['median_pct_pos'], 1)}% "
            f"(n samples={int(row['n_post_samples_ge5cells'])})"
        )
    lines += [
        "",
        "Both genes are epithelial-restricted (about 1% of CD8/NK cells positive, under 1% of T/NK-lineage cells). A whole-sample average therefore moves with epithelial fraction. That all-cell contrast is reported below and is not an independent test of tumor content. The expression values in the main tests are means inside epithelial cells, or inside A3-malignant cells.",
        "",
        "## Post-treatment samples that still contain the MPR cases (epithelial mean, n=12)",
        "",
        f"Median cuts on the epithelial mean log1p(CP10k): TACSTD2 = {fmt_num(summary['epi_tac_median'])}, CLDN4 = {fmt_num(summary['epi_cldn_median'])}. High = at or above that median. Dual-high = both.",
        "",
        f"- **TACSTD2 high vs low.** {prose_composition(tests, 'TACSTD2', epi_cohort, epi_expr)}.",
        f"- **CLDN4 high vs low.** {prose_composition(tests, 'CLDN4', epi_cohort, epi_expr)}.",
        f"- **Dual-high vs not** (n dual-high = {summary['n_dual_epi']}). {prose_composition(tests, 'dual-high', epi_cohort, epi_expr)}.",
        "",
        "Continuous Spearman, same 12 samples:",
        "",
        f"- TACSTD2: {prose_spearman(cont, 'TACSTD2', epi_cohort, epi_expr)}.",
        f"- CLDN4: {prose_spearman(cont, 'CLDN4', epi_cohort, epi_expr)}.",
        "",
        "Expression itself, NMPR versus MPR (pCR counted as MPR). A positive NMPR−MPR difference is higher expression in NMPR, which is the direction stated for TROP2 on the slide.",
        "",
        f"- TACSTD2: {prose_mpr_expr(tests, 'TACSTD2', epi_cohort, epi_expr)}.",
        f"- CLDN4: {prose_mpr_expr(tests, 'CLDN4', epi_cohort, epi_expr)}.",
        "",
        "The TACSTD2 MPR values sit inside the NMPR range (see `figures/fig_epithelial_mpr.png`). The median split puts all 4 MPR samples on the TACSTD2-low side; the continuous test still overlaps. CLDN4 MPR values cover the same span as NMPR.",
        "",
        "Group membership:",
        "",
    ]
    lines += _membership_lines(balance_epi)
    allcell_bits = []
    for gene in ("TACSTD2", "CLDN4"):
        hit = cont[
            (cont["gene"] == gene)
            & (cont["cohort"] == "post_all12_allcell_mean")
            & (cont["endpoint"] == "epithelial fraction")
            & (cont["test"] == "spearman")
        ]
        if not hit.empty:
            row = hit.iloc[0]
            allcell_bits.append(f"{gene} ρ={fmt_num(row['rho'], 2)} (P={fmt_p(row['p'])})")
    lines += [
        "",
        "Whole-sample mean versus epithelial fraction, same 12 samples: " + "; ".join(allcell_bits) + ". "
        "That correlation is the epithelial restriction above. It is not an independent test of tumor content, and it is not used as the tumor-fraction result.",
        "",
        "## A3-malignant subset (not an MPR test)",
        "",
        f"Samples with ≥10 A3-malignant cells: n={summary['n_primary']}. Cuts: TACSTD2 = {fmt_num(summary['tac_median'])}, CLDN4 = {fmt_num(summary['cldn_median'])}. One MPR sample remains (P03).",
        "",
        f"- **TACSTD2.** {prose_composition(tests, 'TACSTD2', a3_cohort, a3_expr)}.",
        f"- **CLDN4.** {prose_composition(tests, 'CLDN4', a3_cohort, a3_expr)}.",
        f"- **Dual-high** (n={summary['n_dual']}). {prose_composition(tests, 'dual-high', a3_cohort, a3_expr)}.",
        "",
        f"- TACSTD2 expression: {prose_mpr_expr(tests, 'TACSTD2', a3_cohort, a3_expr)}.",
        f"- CLDN4 expression: {prose_mpr_expr(tests, 'CLDN4', a3_cohort, a3_expr)}.",
        "",
        "Spearman on this subset:",
        "",
        f"- TACSTD2: {prose_spearman(cont, 'TACSTD2', a3_cohort, a3_expr)}.",
        f"- CLDN4: {prose_spearman(cont, 'CLDN4', a3_cohort, a3_expr)}.",
        "",
        "Membership:",
        "",
    ]
    lines += _membership_lines(balance_a3)
    lines += [
        "",
        "## How the two genes compare",
        "",
        summary["reading"],
        "",
        "## Other cuts",
        "",
        "Same tests under other pre-specified definitions. Not used to pick a CLDN4 result.",
        "",
    ]
    shown_cohorts = {"post_epithelial_n12", "post_ge10_A3"}
    sens = tests[~tests["cohort"].isin(shown_cohorts)]
    keep_tests = {"high_vs_low_composition", "high_vs_low_MPR", "NMPR_vs_MPR_expression"}
    shown = sens[sens["test"].isin(keep_tests)]
    for _, row in shown.iterrows():
        if row["test"] == "NMPR_vs_MPR_expression":
            lines.append(
                f"- {row['cohort']} | {row['gene']} | {row['expression']} | NMPR vs MPR expression: "
                f"medians {fmt_num(row['median_nmpr'])} vs {fmt_num(row['median_mpr'])} "
                f"(n={int(row['n_nmpr'])} vs {int(row['n_mpr'])}), P={fmt_p(row['p'])}"
            )
        elif row["test"] == "high_vs_low_MPR":
            lines.append(
                f"- {row['cohort']} | {row['gene']} | {row['expression']} | MPR rate high "
                f"{int(row['n_high_mpr'])}/{int(row['n_high'])} vs low {int(row['n_low_mpr'])}/{int(row['n_low'])}, "
                f"Fisher P={fmt_p(row['p'])}"
            )
        elif row["endpoint"] in (
            "A3-malignant fraction",
            "epithelial fraction",
            "PTPRC-negative malignant fraction",
            "CD8 T + NK fraction",
            "all T + NK fraction",
        ):
            lines.append(
                f"- {row['cohort']} | {row['gene']} | {row['expression']} | {row['endpoint']}: "
                f"high {fmt_num(row['median_high'])} vs low {fmt_num(row['median_low'])} "
                f"(n={int(row['n_high'])} vs {int(row['n_low'])}), Δ={fmt_num(row['delta_median_high_minus_low'])}, "
                f"P={fmt_p(row['p'])}"
            )
    lines += [
        "",
        "All-cell means versus epithelial fraction are in `tables/spearman.tsv` under cohort `post_all12_allcell_mean`. Full test rows: `tables/tests.tsv`.",
        "",
        "## Caveats",
        "",
        "- Barcode labels are reconstructed. This is not a cell-for-cell replay of the slide, and the slide's stated P values are not copied forward as results.",
        "- n=12 post-treatment samples. A median split is one cut; the Spearman column does not depend on that cut. With 4 MPR samples, a Fisher P near 0.06 is a small-sample count.",
        "- P07 (squamous, NMPR) is an epithelial-fraction outlier (about half the cells). Medians, not means, are the high-vs-low summaries.",
        "- CD8/NK is marker-defined: CD3+ and CD8A UMI > 0, or CD3-negative with NCR1, KLRF1, or GNLY together with NKG7. The all-T/NK column is the lineage-score compartment used in the earlier CLDN4-only run.",
        "- Adenocarcinoma versus squamous is unbalanced inside some splits. Counts are listed; there is no multivariable model at this n.",
        "- Dual-high is the intersection of the two median groups. A dual-high count is not a CLDN4 result.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def reading_paragraph(tests: pd.DataFrame, cont: pd.DataFrame) -> str:
    """Factual comparison. Does not require CLDN4 to follow TACSTD2."""
    cohort = "post_epithelial_n12"
    expr = "epithelial_mean_log1p_cp10k"

    def delta(gene, endpoint):
        row = lookup(tests, gene=gene, cohort=cohort, expression=expr, test="high_vs_low_composition", endpoint=endpoint)
        if row is None or row["delta_median_high_minus_low"] is None:
            return None, None
        return float(row["delta_median_high_minus_low"]), float(row["p"]) if pd.notna(row["p"]) else None

    def expr_delta(gene):
        row = lookup(tests, gene=gene, cohort=cohort, expression=expr, test="NMPR_vs_MPR_expression")
        if row is None:
            return None, None
        return float(row["delta_median_nmpr_minus_mpr"]), float(row["p"]) if pd.notna(row["p"]) else None

    def rho(gene, endpoint):
        hit = cont[
            (cont["gene"] == gene)
            & (cont["cohort"] == cohort)
            & (cont["expression"] == expr)
            & (cont["endpoint"] == endpoint)
            & (cont["test"] == "spearman")
        ]
        if hit.empty:
            return None, None
        return hit.iloc[0]["rho"], hit.iloc[0]["p"]

    t_tumor, t_tumor_p = delta("TACSTD2", "epithelial fraction")
    c_tumor, c_tumor_p = delta("CLDN4", "epithelial fraction")
    t_cd8, t_cd8_p = delta("TACSTD2", "CD8 T + NK fraction")
    c_cd8, c_cd8_p = delta("CLDN4", "CD8 T + NK fraction")
    t_mpr, t_mpr_p = expr_delta("TACSTD2")
    c_mpr, c_mpr_p = expr_delta("CLDN4")
    tr, tp = rho("TACSTD2", "CD8 T + NK fraction")
    cr, cp = rho("CLDN4", "CD8 T + NK fraction")
    tt, ttp = rho("TACSTD2", "all T + NK fraction")
    ct, ctp = rho("CLDN4", "all T + NK fraction")
    return (
        f"On all 12 post-treatment samples, a higher epithelial TACSTD2 does not bring a higher epithelial fraction "
        f"(high−low Δ={fmt_num(t_tumor)}, P={fmt_p(t_tumor_p)}) or a lower CD8 T/NK fraction "
        f"(Δ={fmt_num(t_cd8)}, P={fmt_p(t_cd8_p)}; Spearman ρ={fmt_num(tr, 2)}, P={fmt_p(tp)}). "
        f"The same is true for CLDN4 (tumor-fraction Δ={fmt_num(c_tumor)}, P={fmt_p(c_tumor_p)}; "
        f"CD8/NK Δ={fmt_num(c_cd8)}, P={fmt_p(c_cd8_p)}; Spearman ρ={fmt_num(cr, 2)}, P={fmt_p(cp)}). "
        f"CLDN4 versus the lineage T/NK fraction is ρ={fmt_num(ct, 2)} (P={fmt_p(ctp)}); "
        f"TACSTD2 versus that fraction is ρ={fmt_num(tt, 2)} (P={fmt_p(ttp)}). "
        f"The MPR contrast is where the genes differ. Epithelial TACSTD2 is higher in NMPR than in MPR "
        f"(Δ={fmt_num(t_mpr)}, P={fmt_p(t_mpr_p)}), and the TACSTD2-high half contains no MPR sample. "
        f"Epithelial CLDN4 does not shift between NMPR and MPR (Δ={fmt_num(c_mpr)}, P={fmt_p(c_mpr_p)}), "
        f"and the CLDN4 median split is 2 MPR versus 2 MPR. "
        "Dual-high is reported as its own row; it is not used to rewrite the CLDN4 estimate."
    )


def sanity(cells: pd.DataFrame) -> dict:
    mal = cells[cells["is_malignant"]]
    tnk = cells[cells["is_tnk_lineage"]]
    epi = cells[cells["is_epithelial"]]
    return {
        "n_epithelial": int(len(epi)),
        "n_malignant": int(len(mal)),
        "n_malignant_strict": int(cells["is_malignant_strict"].sum()),
        "n_cd8nk": int(cells["is_cd8nk"].sum()),
        "n_tnk_lineage": int(len(tnk)),
        "lineage_counts": {k: int(v) for k, v in cells["lineage"].value_counts().items()},
        "EPCAM_pct_malignant": float(100 * (mal["EPCAM"] > 0).mean()) if len(mal) else None,
        "EPCAM_pct_tnk": float(100 * (tnk["EPCAM"] > 0).mean()) if len(tnk) else None,
        "PTPRC_pct_malignant": float(100 * (mal["PTPRC"] > 0).mean()) if len(mal) else None,
        "PTPRC_pct_tnk": float(100 * (tnk["PTPRC"] > 0).mean()) if len(tnk) else None,
        "unmatched_samples": sorted(set(cells.loc[cells["Patient"].isna(), "Sample"].unique())),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path)
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--per-sample", type=Path, help="Skip the matrix stream and reuse per_sample_all.tsv")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    figdir = args.outdir / "figures"
    tabdir = args.outdir / "tables"
    figdir.mkdir(exist_ok=True)
    tabdir.mkdir(exist_ok=True)

    prior_summary = {}
    summary_path = args.outdir / "summary.json"
    if summary_path.exists():
        prior_summary = json.loads(summary_path.read_text())

    if args.per_sample:
        per = pd.read_csv(args.per_sample, sep="\t")

        def as_bool(series: pd.Series) -> pd.Series:
            if series.dtype == bool:
                return series
            return series.astype(str).str.strip().str.lower().isin(["true", "1"])

        per["is_post"] = as_bool(per["is_post"])
        per["is_pre"] = as_bool(per["is_pre"])
        restrict = pd.read_csv(tabdir / "expression_by_compartment.tsv", sep="\t")
        checks = prior_summary.get("sanity", {})
        missing = prior_summary.get("missing_markers", [])
        n_gene_rows = int(prior_summary["n_genes_matrix"])
        n_cells = int(prior_summary["n_cells"])
        median_genes = float(prior_summary["median_genes_per_cell"])
        matrix_bytes = int(prior_summary["matrix_bytes"])
        matrix_sha = prior_summary["matrix_sha256"]
        cells = None
    else:
        if args.matrix is None or args.metadata is None:
            raise SystemExit("pass --matrix and --metadata, or --per-sample")
        meta = load_metadata(args.metadata)
        meta.to_csv(tabdir / "sample_metadata.tsv", sep="\t", index=False)
        wanted = marker_universe()
        print(f"streaming matrix; extracting {len(wanted)} genes", flush=True)
        cell_ids, expr, n_umi, n_genes, n_gene_rows = stream_matrix(args.matrix, wanted)
        missing = sorted(wanted - set(expr))
        print(f"cells={len(cell_ids)} missing_markers={missing}", flush=True)
        cells = build_cells(cell_ids, expr, n_umi, n_genes, meta)
        checks = sanity(cells)
        print(json.dumps({k: checks[k] for k in ["n_epithelial", "n_malignant", "n_cd8nk", "n_tnk_lineage", "PTPRC_pct_malignant", "unmatched_samples"]}), flush=True)
        if checks["unmatched_samples"]:
            raise SystemExit(f"barcodes did not match metadata: {checks['unmatched_samples']}")
        per = per_sample_table(cells)
        per.to_csv(tabdir / "per_sample_all.tsv", sep="\t", index=False)
        restrict = restriction_table(cells)
        restrict.to_csv(tabdir / "expression_by_compartment.tsv", sep="\t", index=False)
        pd.crosstab(cells["Sample"], cells["lineage"]).to_csv(tabdir / "lineage_by_sample.tsv", sep="\t")
        n_cells = int(len(cells))
        median_genes = float(np.median(n_genes))
        matrix_bytes = args.matrix.stat().st_size
        matrix_sha = __import__("hashlib").sha256(args.matrix.read_bytes()).hexdigest()

    post = per[(per["is_post"]) & (per["n_malignant"] >= MIN_MALIGNANT)].copy()
    dropped_post = sorted(per.loc[per["is_post"] & (per["n_malignant"] < MIN_MALIGNANT), "Patient"].astype(str))
    # Response-known, including pre-treatment biopsies that carry an MPR/NMPR label.
    known = per[(per["response_paper"].isin(["MPR", "NMPR"])) & (per["n_malignant"] >= MIN_MALIGNANT)].copy()
    epi_post = per[(per["is_post"]) & (per["n_epithelial"] >= MIN_MALIGNANT)].copy()
    strict_post = per[(per["is_post"]) & (per["n_malignant_strict"] >= MIN_MALIGNANT)].copy()

    primary, t1, c1 = run_block(
        post, "post_ge10_A3", "A3_malignant_mean_log1p_cp10k",
        "tacstd2_mal_mean", "cldn4_mal_mean", "frac_tumor", "frac_cd8nk", "frac_tnk",
    )
    _, t2, c2 = run_block(
        known, "response_known_ge10_A3", "A3_malignant_mean_log1p_cp10k",
        "tacstd2_mal_mean", "cldn4_mal_mean", "frac_tumor", "frac_cd8nk", "frac_tnk",
    )
    epi, t3, c3 = run_block(
        epi_post, "post_epithelial_n12", "epithelial_mean_log1p_cp10k",
        "tacstd2_epi_mean", "cldn4_epi_mean", "frac_epithelial", "frac_cd8nk", "frac_tnk",
        tumor_label="epithelial fraction",
    )
    _, t4, c4 = run_block(
        strict_post, "post_ge10_A3_PTPRCneg", "A3_PTPRCneg_mean_log1p_cp10k",
        "tacstd2_strict_mean", "cldn4_strict_mean", "frac_tumor_strict", "frac_cd8nk", "frac_tnk",
        tumor_label="PTPRC-negative malignant fraction",
    )
    # All-cell mean on every post-treatment sample. Tumor-fraction contrasts from
    # this average are not independent of epithelial content.
    post_all = per[per["is_post"]].copy()
    _, t5, c5 = run_block(
        post_all, "post_all12_allcell_mean", "all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction",
        "tacstd2_all_mean", "cldn4_all_mean", "frac_tumor", "frac_cd8nk", "frac_tnk",
    )
    for gene, col in [("TACSTD2", "tacstd2_all_mean"), ("CLDN4", "cldn4_all_mean")]:
        rec = spearman(post_all[col], post_all["frac_epithelial"])
        c5.append(
            {
                "gene": gene,
                "cohort": "post_all12_allcell_mean",
                "expression": "all_cell_mean_log1p_cp10k_NOT_independent_of_tumor_fraction",
                "test": "spearman",
                "endpoint": "epithelial fraction",
                "endpoint_col": "frac_epithelial",
                "n": rec["n"],
                "rho": rec["rho"],
                "p": rec["p"],
            }
        )

    tests = pd.DataFrame(t1 + t2 + t3 + t4 + t5)
    cont = pd.DataFrame(c1 + c2 + c3 + c4 + c5)
    tests.to_csv(tabdir / "tests.tsv", sep="\t", index=False)
    cont.to_csv(tabdir / "spearman.tsv", sep="\t", index=False)
    epi.to_csv(tabdir / "per_sample_epithelial.tsv", sep="\t", index=False)
    primary.to_csv(tabdir / "per_sample_a3.tsv", sep="\t", index=False)
    balance_epi = balance_table(epi)
    balance_a3 = balance_table(primary)
    balance_epi.insert(0, "cohort", "post_epithelial_n12")
    balance_a3.insert(0, "cohort", "post_ge10_A3")
    pd.concat([balance_epi, balance_a3], ignore_index=True).to_csv(tabdir / "group_balance.tsv", sep="\t", index=False)

    make_figures(epi, primary, figdir)

    summary = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "matrix_bytes": matrix_bytes,
        "matrix_sha256": matrix_sha,
        "n_cells": n_cells,
        "n_genes_matrix": int(n_gene_rows),
        "median_genes_per_cell": median_genes,
        "missing_markers": missing,
        "n_primary": int(len(primary)),
        "n_epithelial_post": int(len(epi)),
        "dropped_post": dropped_post,
        "n_dual": int(primary["dual_high"].sum()),
        "n_dual_epi": int(epi["dual_high"].sum()),
        "tac_median": float(primary["tacstd2_median_cut"].iloc[0]),
        "cldn_median": float(primary["cldn4_median_cut"].iloc[0]),
        "epi_tac_median": float(epi["tacstd2_median_cut"].iloc[0]),
        "epi_cldn_median": float(epi["cldn4_median_cut"].iloc[0]),
        "sanity": checks,
        "primary_patients": primary["Patient"].tolist(),
        "epithelial_patients": epi["Patient"].tolist(),
    }
    summary["reading"] = reading_paragraph(tests, cont)
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    write_finding(args.outdir / "FINDING.md", summary, tests, cont, balance_epi, balance_a3, restrict)
    print(summary["reading"])
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
