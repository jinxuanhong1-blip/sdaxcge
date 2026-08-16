"""Verification / sanity checks (positive controls) for the analysis.

These do not test the hypothesis; they confirm that genes, compartments and
normalisation behave as biology dictates, so the main results can be trusted:
  V1  TACSTD2/CLDN4 are strongly enriched in epithelial vs immune cells (scRNA).
  V2  Marker specificity: EPCAM up in epithelial, PTPRC up in immune.
  V3  TACSTD2 and CLDN4 positively correlate across bulk tumors (co-epithelial).
  V4  Group sizes match the deposited sample metadata.
  V5  Pipeline is deterministic: cached re-run reproduces pseudobulk values.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import config as C


def main() -> None:
    checks = []

    # ---- scRNA-based checks ----
    cells = pd.read_csv(C.TABLES_DIR / "gse207422_scrna_cells_annotated.csv.gz")
    epi = cells[cells.is_epithelial]
    imm = cells[(cells.PTPRC_raw > 0) & (cells.EPCAM_raw == 0)]

    for gene in C.TARGET_GENES:
        me, mi = epi[gene].mean(), imm[gene].mean()
        u, p = stats.mannwhitneyu(epi[gene], imm[gene], alternative="greater")
        checks.append({"check": f"V1 {gene} epithelial>immune",
                       "epi_mean": round(me, 3), "immune_mean": round(mi, 3),
                       "p_value": p, "pass": bool(me > mi and p < 1e-6)})

    # V2 marker specificity
    for mk, hi_epi in [("EPCAM_raw", True), ("PTPRC_raw", False)]:
        me = (epi[mk] > 0).mean() * 100
        mi = (imm[mk] > 0).mean() * 100
        ok = (me > mi) if hi_epi else (mi > me)
        checks.append({"check": f"V2 {mk} specificity",
                       "epi_pctpos": round(me, 1), "immune_pctpos": round(mi, 1),
                       "pass": bool(ok)})

    # ---- V3 TACSTD2/CLDN4 co-expression across bulk tumors ----
    bulk = pd.read_csv(C.raw_path("GSE207422_bulk_expr"), sep="\t", index_col=0)
    t = bulk.loc["TACSTD2"].astype(float)
    c = bulk.loc["CLDN4"].astype(float)
    r, p = stats.spearmanr(t, c)
    checks.append({"check": "V3 TACSTD2~CLDN4 corr (GSE207422 bulk)",
                   "spearman_r": round(float(r), 3), "p_value": float(p),
                   "pass": bool(r > 0)})

    # ---- V4 group sizes ----
    pseudo = pd.read_csv(C.TABLES_DIR / "gse207422_scrna_pseudobulk_epithelial.csv")
    n_pre = int((pseudo.timepoint == "pre").sum())
    n_post = int((pseudo.timepoint == "post").sum())
    checks.append({"check": "V4 scRNA sample counts (pre/post)",
                   "n_pre": n_pre, "n_post": n_post,
                   "pass": bool(n_pre == 3 and n_post == 12)})

    # ---- V5 determinism ----
    epi_grp = epi.groupby("sample")["TACSTD2"].mean()
    ref = pseudo.set_index("sample")["TACSTD2_mean"]
    common = ref.index.intersection(epi_grp.index)
    max_abs = float(np.max(np.abs(ref.loc[common].values - epi_grp.loc[common].values)))
    checks.append({"check": "V5 pseudobulk reproducible",
                   "max_abs_diff": max_abs, "pass": bool(max_abs < 1e-9)})

    # ---- V6 GSE248249 pair inventory matches the series record ----
    pairs = pd.read_csv(C.TABLES_DIR / "gse248249_patient_pairs.csv")
    n_pat = pairs[pairs.gene == "TACSTD2"]["patient"].nunique()
    checks.append({"check": "V6 GSE248249 13 same-patient pairs",
                   "n_pairs": int(n_pat), "pass": bool(n_pat == 13)})

    # ---- V7 Clariom probes present and finite ----
    long = pd.read_csv(C.TABLES_DIR / "gse248249_paired_long.csv")
    finite = bool(np.isfinite(long["TACSTD2"]).all() and np.isfinite(long["CLDN4"]).all())
    checks.append({"check": "V7 GSE248249 TACSTD2/CLDN4 finite",
                   "n_samples": int(len(long)), "pass": finite and len(long) == 42})

    report = pd.DataFrame(checks)
    report.to_csv(C.TABLES_DIR / "verification_report.csv", index=False)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(report.to_string())
    n_pass = int(report["pass"].sum())
    print(f"\nVERIFICATION: {n_pass}/{len(report)} checks passed.")
    if n_pass != len(report):
        raise SystemExit("One or more verification checks FAILED")


if __name__ == "__main__":
    main()
