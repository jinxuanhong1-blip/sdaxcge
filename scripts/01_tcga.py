#!/usr/bin/env python3
"""TCGA LUAD / LUSC bulk RNA-seq (UCSC Xena GDC hub, STAR log2(tpm+1)).

Primary solid tumours only (barcode sample-type code 01), one sample per patient.
Usage: python3 scripts/01_tcga.py LUAD
       python3 scripts/01_tcga.py LUSC
"""
import os
import sys
import pandas as pd
import lib_analysis as L

RAW = "data/raw"


def load_tcga(cohort):
    path = f"{RAW}/TCGA-{cohort}.star_tpm.tsv.gz"
    print(f"Loading {path} ...")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    print("  matrix:", expr.shape)

    pm = pd.read_csv(f"{RAW}/gencode.v36.gene.probemap", sep="\t")
    id2sym = dict(zip(pm["id"], pm["gene"]))
    expr.index = [id2sym.get(i, i) for i in expr.index]
    expr = expr.groupby(level=0).max()

    def sample_code(bc):
        parts = bc.split("-")
        return parts[3][:2] if len(parts) > 3 else "NA"

    tumors = [c for c in expr.columns if sample_code(c) == "01"]
    expr = expr[tumors]
    keep, seen = [], set()
    for c in expr.columns:
        patient = "-".join(c.split("-")[:3])
        if patient not in seen:
            seen.add(patient)
            keep.append(c)
    expr = expr[keep]
    print("  primary tumours (unique patients):", expr.shape[1])
    return expr


def run_cohort(cohort):
    out = f"results/claim_A8A9/tcga_{cohort.lower()}"
    os.makedirs(out, exist_ok=True)
    expr = load_tcga(cohort)
    if L.TROP2 not in expr.index:
        raise SystemExit(f"{L.TROP2} missing from {cohort}")
    trop2 = expr.loc[L.TROP2]
    stats = L.per_gene_stats(expr, trop2, min_frac_expressed=0.10)
    stats.to_csv(f"{out}/per_gene_stats.csv")
    print("  genes tested:", stats.shape[0])
    panel = stats.reindex(L.PANEL)
    panel.to_csv(f"{out}/panel_stats.csv")
    print(panel[["log2FC_high_vs_low", "spearman_rho", "welch_fdr", "up_in_trop2_high"]])
    gsea = L.run_prerank_gsea(stats, out)
    print(gsea[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))
    print(f"TCGA-{cohort} done.")


def main():
    cohorts = sys.argv[1:] or ["LUAD", "LUSC"]
    for c in cohorts:
        run_cohort(c.upper())


if __name__ == "__main__":
    main()
