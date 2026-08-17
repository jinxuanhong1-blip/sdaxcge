#!/usr/bin/env python3
"""Score Cldn4 and IFN/MHC on the public GSE297630 processed matrix.

Reproduces the given Tacstd2 log2FC +2.05 (anti-PD-1 tolerant P vs control C,
n=3 vs 3) from GSE297630_processed_data.xlsx Signal columns, then scores
CLDN4 and pre-specified IFN/MHC sets on the same table.

Usage:
  python3 methods/gse297630_cldn4/score_gse297630.py
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import SETS, TARGETS

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
CACHE = HERE / "_cache"
XLSX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297630/"
    "suppl/GSE297630_processed_data.xlsx"
)
XLSX_NAME = "GSE297630_processed_data.xlsx"

CTRL = [
    "C-1_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
    "C-2_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
    "C-3_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
]
TRT = [
    "P-1_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
    "P-2_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
    "P-3_(Clariom_S_Mouse).sst-rma-gene-full.chp Signal",
]
SAMPLE_LABELS = ["C1", "C2", "C3", "P1", "P2", "P3"]


def fetch_xlsx() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / XLSX_NAME
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    print(f"downloading {XLSX_URL}")
    with urllib.request.urlopen(XLSX_URL, timeout=180) as resp:
        path.write_bytes(resp.read())
    return path


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Expression", skiprows=4)
    d = df.dropna(subset=["Gene Symbol"]).copy()
    d["sym"] = d["Gene Symbol"].astype(str).str.upper().str.strip()
    d = d[~d["sym"].str.match(r"^\d")]
    mat = d.groupby("sym")[CTRL + TRT].mean()
    mat.columns = SAMPLE_LABELS
    return mat


def welch_mwu(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    tw = stats.ttest_ind(a, b, equal_var=False)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(tw.statistic), float(tw.pvalue), float(mw.pvalue)


def gene_row(mat: pd.DataFrame, gene: str) -> dict:
    c = mat.loc[gene, ["C1", "C2", "C3"]].astype(float).to_numpy()
    p = mat.loc[gene, ["P1", "P2", "P3"]].astype(float).to_numpy()
    t, tp, up = welch_mwu(p, c)
    return {
        "gene": gene,
        "n_control": 3,
        "n_tolerant": 3,
        "mean_control": round(float(c.mean()), 4),
        "mean_tolerant": round(float(p.mean()), 4),
        "sd_control": round(float(c.std(ddof=1)), 4),
        "sd_tolerant": round(float(p.std(ddof=1)), 4),
        "log2FC_P_minus_C": round(float(p.mean() - c.mean()), 4),
        "welch_t": round(t, 4),
        "welch_p": tp,
        "mwu_p": up,
        "C1": float(c[0]),
        "C2": float(c[1]),
        "C3": float(c[2]),
        "P1": float(p[0]),
        "P2": float(p[1]),
        "P3": float(p[2]),
    }


def zscore_rows(m: pd.DataFrame) -> pd.DataFrame:
    return m.sub(m.mean(axis=1), axis=0).div(
        m.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )


def set_row(mat: pd.DataFrame, name: str, genes: list[str]) -> dict:
    used = [g for g in genes if g in mat.index]
    missing = [g for g in genes if g not in mat.index]
    z = zscore_rows(mat.loc[used])
    sc = z.mean(axis=0)
    c = sc[["C1", "C2", "C3"]].astype(float).to_numpy()
    p = sc[["P1", "P2", "P3"]].astype(float).to_numpy()
    t, tp, up = welch_mwu(p, c)
    raw = mat.loc[used].mean(axis=0)
    return {
        "set": name,
        "n_requested": len(genes),
        "n_present": len(used),
        "n_missing": len(missing),
        "missing": ";".join(missing),
        "n_control": 3,
        "n_tolerant": 3,
        "mean_z_control": round(float(c.mean()), 4),
        "mean_z_tolerant": round(float(p.mean()), 4),
        "delta_z_P_minus_C": round(float(p.mean() - c.mean()), 4),
        "raw_mean_log2_control": round(float(raw[["C1", "C2", "C3"]].mean()), 4),
        "raw_mean_log2_tolerant": round(float(raw[["P1", "P2", "P3"]].mean()), 4),
        "raw_log2FC": round(
            float(raw[["P1", "P2", "P3"]].mean() - raw[["C1", "C2", "C3"]].mean()),
            4,
        ),
        "welch_t": round(t, 4),
        "welch_p": tp,
        "mwu_p": up,
        "z_C1": round(float(c[0]), 4),
        "z_C2": round(float(c[1]), 4),
        "z_C3": round(float(c[2]), 4),
        "z_P1": round(float(p[0]), 4),
        "z_P2": round(float(p[1]), 4),
        "z_P3": round(float(p[2]), 4),
        "genes": ";".join(used),
    }


def spearman_pair(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    rho, p = stats.spearmanr(a.astype(float), b.astype(float))
    return float(rho), float(p)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mat = load_matrix(fetch_xlsx())

    target_rows = [gene_row(mat, g) for g in TARGETS]
    pd.DataFrame(target_rows).to_csv(OUT / "targets.tsv", sep="\t", index=False)

    set_rows = [set_row(mat, name, genes) for name, genes in SETS.items()]
    pd.DataFrame(set_rows).to_csv(OUT / "set_scores.tsv", sep="\t", index=False)

    union = []
    seen = set()
    for genes in SETS.values():
        for g in genes:
            if g not in seen and g in mat.index:
                seen.add(g)
                union.append(gene_row(mat, g))
    gene_df = pd.DataFrame(union).sort_values("log2FC_P_minus_C", ascending=False)
    gene_df.to_csv(OUT / "ifn_mhc_genes.tsv", sep="\t", index=False)

    per_sample = pd.DataFrame(
        {
            "sample": SAMPLE_LABELS,
            "group": ["control"] * 3 + ["antiPD1_tolerant"] * 3,
            "TACSTD2": mat.loc["TACSTD2", SAMPLE_LABELS].astype(float).to_numpy(),
            "CLDN4": mat.loc["CLDN4", SAMPLE_LABELS].astype(float).to_numpy(),
        }
    )
    for name, genes in SETS.items():
        used = [g for g in genes if g in mat.index]
        z = zscore_rows(mat.loc[used]).mean(axis=0)
        per_sample[f"z_{name}"] = z[SAMPLE_LABELS].astype(float).to_numpy()
    per_sample.to_csv(OUT / "per_sample.tsv", sep="\t", index=False)

    corr_rows = []
    for left, right in [
        ("TACSTD2", "CLDN4"),
        ("TACSTD2", "z_CORE6"),
        ("TACSTD2", "z_IFN_ISG"),
        ("TACSTD2", "z_MHC_I"),
        ("TACSTD2", "z_MHC_APM"),
        ("TACSTD2", "z_IFN_MHC"),
        ("CLDN4", "z_CORE6"),
        ("CLDN4", "z_IFN_ISG"),
        ("CLDN4", "z_MHC_I"),
        ("CLDN4", "z_MHC_APM"),
        ("CLDN4", "z_IFN_MHC"),
    ]:
        rho, p = spearman_pair(per_sample[left], per_sample[right])
        corr_rows.append(
            {
                "x": left,
                "y": right,
                "n": 6,
                "spearman_rho": round(rho, 4),
                "spearman_p": p,
                "note": "n=6; treatment-driven; one-or-two-point fragile",
            }
        )
    pd.DataFrame(corr_rows).to_csv(OUT / "correlations.tsv", sep="\t", index=False)

    ifn_fc = gene_df[gene_df["gene"].isin(SETS["IFN_ISG"])]["log2FC_P_minus_C"]
    mhci_fc = gene_df[gene_df["gene"].isin(SETS["MHC_I"])]["log2FC_P_minus_C"]
    tac = next(r for r in target_rows if r["gene"] == "TACSTD2")
    cld = next(r for r in target_rows if r["gene"] == "CLDN4")
    by_set = {r["set"]: r for r in set_rows}

    summary = {
        "accession": "GSE297630",
        "contrast": "anti-PD-1 tolerant (P, n=3) vs control (C, n=3)",
        "matrix": "GSE297630_processed_data.xlsx Expression Signal columns",
        "honest_n": "3 vs 3; two-sided MWU cannot go below 0.1",
        "tacstd2_given_log2FC": 2.05,
        "tacstd2_reproduced_log2FC": tac["log2FC_P_minus_C"],
        "tacstd2_welch_p": tac["welch_p"],
        "tacstd2_mwu_p": tac["mwu_p"],
        "cldn4_log2FC": cld["log2FC_P_minus_C"],
        "cldn4_welch_p": cld["welch_p"],
        "cldn4_mwu_p": cld["mwu_p"],
        "ifn_isg_delta_z": by_set["IFN_ISG"]["delta_z_P_minus_C"],
        "ifn_isg_welch_p": by_set["IFN_ISG"]["welch_p"],
        "ifn_isg_up_down": f"{int((ifn_fc > 0).sum())} up / {int((ifn_fc < 0).sum())} down",
        "mhc_i_delta_z": by_set["MHC_I"]["delta_z_P_minus_C"],
        "mhc_i_welch_p": by_set["MHC_I"]["welch_p"],
        "mhc_i_all_up": bool((mhci_fc > 0).all()),
        "ifn_mhc_delta_z": by_set["IFN_MHC"]["delta_z_P_minus_C"],
        "ifn_mhc_welch_p": by_set["IFN_MHC"]["welch_p"],
        "n_genes_matrix": int(mat.shape[0]),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
