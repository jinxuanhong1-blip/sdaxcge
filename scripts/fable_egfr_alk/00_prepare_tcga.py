"""Prepare a compact, analysis-ready TCGA table for one cohort (LUAD or LUSC).

Reads the large Xena/GDC STAR-TPM matrix in chunks, keeps only the genes we need,
derives EGFR / ALK driver status from the somatic-mutation file (plus an EML4-ALK
expression-outlier proxy), merges overall survival, and writes a small per-sample
table under results/fable_egfr_alk/processed/.

Usage:
    python 00_prepare_tcga.py LUAD
    python 00_prepare_tcga.py LUSC
"""
from __future__ import annotations

import sys
import gzip
import numpy as np
import pandas as pd

import common as C

NONSILENT = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant",
    "protein_altering_variant", "coding_sequence_variant",
}


def wanted_ensembl() -> dict:
    """symbol -> Ensembl(unversioned) for every gene we extract."""
    m = {}
    m.update(C.TARGETS)
    m.update(C.DRIVERS)
    m.update(C.IMMUNE_GENES)
    return m


def extract_expression(tpm_path: str, ens2sym: dict) -> pd.DataFrame:
    """Return samples x gene-symbol expression (log2 TPM+1) for wanted genes."""
    keep_rows = []
    header = None
    with gzip.open(tpm_path, "rt") as fh:
        for chunk in pd.read_csv(fh, sep="\t", chunksize=5000, dtype=str):
            if header is None:
                header = list(chunk.columns)
            ens_nover = chunk["Ensembl_ID"].str.replace(r"\.\d+$", "", regex=True)
            mask = ens_nover.isin(ens2sym.keys())
            if mask.any():
                sub = chunk[mask].copy()
                sub["ens"] = ens_nover[mask].values
                keep_rows.append(sub)
    if not keep_rows:
        raise RuntimeError("No target genes found in expression matrix")
    expr = pd.concat(keep_rows, ignore_index=True)
    expr["symbol"] = expr["ens"].map(ens2sym)
    expr = expr.drop(columns=["Ensembl_ID", "ens"]).set_index("symbol")
    expr = expr.astype(float)
    # samples as rows, genes as columns
    mat = expr.T
    mat.index.name = "sample"
    return mat


def driver_status(mut_path: str) -> pd.DataFrame:
    mut = pd.read_csv(mut_path, sep="\t")
    mut = mut[mut["effect"].isin(NONSILENT)]
    egfr = set(mut.loc[mut["gene"] == "EGFR", "sample"])
    alk = set(mut.loc[mut["gene"] == "ALK", "sample"])
    all_samples = set(mut["sample"])
    df = pd.DataFrame({"sample": sorted(all_samples)})
    df["EGFR_mut"] = df["sample"].isin(egfr).astype(int)
    df["ALK_mut"] = df["sample"].isin(alk).astype(int)
    return df


def main(cohort: str) -> None:
    tpm = f"{C.DATA}/TCGA-{cohort}.star_tpm.tsv.gz"
    mut = f"{C.DATA}/TCGA-{cohort}.somaticmutation_wxs.tsv.gz"
    surv = f"{C.DATA}/TCGA-{cohort}.survival.tsv.gz"

    ens2sym = {v: k for k, v in wanted_ensembl().items()}
    expr = extract_expression(tpm, ens2sym)

    # Restrict to primary tumours (barcode sample-type code 01) to avoid mixing
    # normal / metastatic tissue.
    sample_type = expr.index.to_series().str[13:15]
    expr = expr[sample_type == "01"].copy()

    # Driver status (mutation barcodes carry the -01A suffix; align on 15-char id)
    drv = driver_status(mut)
    drv["sample15"] = drv["sample"].str[:15]
    drv_g = drv.groupby("sample15")[["EGFR_mut", "ALK_mut"]].max().reset_index()

    tab = expr.reset_index()
    tab["sample15"] = tab["sample"].str[:15]
    tab = tab.merge(drv_g, on="sample15", how="left")
    # Samples not present in the mutation file are not WXS-profiled -> unknown.
    tab["mut_profiled"] = tab["sample15"].isin(set(drv["sample15"])).astype(int)
    tab["EGFR_mut"] = tab["EGFR_mut"].fillna(0).astype(int)
    tab["ALK_mut"] = tab["ALK_mut"].fillna(0).astype(int)

    # EML4-ALK fusion proxy: aberrant high ALK expression.  ALK is essentially
    # silent in normal lung / most LUAD; fusion drives strong 3'-kinase-domain
    # expression.  Flag samples above (median + 3*MAD) and in the top 5%.
    alk = tab["ALK"].astype(float)
    mad = (alk - alk.median()).abs().median() * 1.4826
    thr_mad = alk.median() + 3 * (mad if mad > 0 else alk.std())
    thr_pct = alk.quantile(0.95)
    tab["ALK_expr_outlier"] = ((alk > thr_mad) & (alk > thr_pct)).astype(int)

    # ALK-driven == EML4-ALK fusion proxy (expression outlier).  ALK point
    # mutations are NOT used as a driver definition: ALK is a large gene and
    # most non-silent SNVs in this smoking-enriched cohort are passengers.  The
    # ALK_mut column is retained separately for transparency only.
    tab["ALK_driven"] = tab["ALK_expr_outlier"].astype(int)

    # Driver group label (EGFR takes precedence when co-occurring)
    def grp(r):
        if r["EGFR_mut"] == 1:
            return "EGFR_mut"
        if r["ALK_driven"] == 1:
            return "ALK_driven"
        if r["mut_profiled"] == 1:
            return "EGFR_ALK_wt"
        return "unknown"

    tab["driver_group"] = tab.apply(grp, axis=1)

    # Survival
    sv = pd.read_csv(surv, sep="\t")
    sv["sample15"] = sv["sample"].str[:15]
    sv_g = sv.groupby("sample15")[["OS", "OS.time"]].first().reset_index()
    tab = tab.merge(sv_g, on="sample15", how="left")

    tab["cohort"] = cohort
    out = f"{C.PROCESSED}/tcga_{cohort.lower()}_samples.csv"
    tab.to_csv(out, index=False)

    # brief console summary
    print(f"[{cohort}] primary tumours: {len(tab)}")
    print(tab["driver_group"].value_counts().to_string())
    print(f"ALK expr-outlier: {int(tab['ALK_expr_outlier'].sum())} | "
          f"ALK point-mut: {int(tab['ALK_mut'].sum())} | "
          f"EGFR-mut: {int(tab['EGFR_mut'].sum())}")
    print(f"wrote {out}  shape={tab.shape}")


if __name__ == "__main__":
    cohort = (sys.argv[1] if len(sys.argv) > 1 else "LUAD").upper()
    main(cohort)
