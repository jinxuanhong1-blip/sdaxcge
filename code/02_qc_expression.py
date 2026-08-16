#!/usr/bin/env python3
"""Sanity checks on the ICB compendium before any CLDN4 analysis.

Checks
  1. expression units / value range (are these log2 TPM?)
  2. sample ids in the expression matrix match the metadata patient ids
  3. CLDN4 is present, and how detectable it is per study (epithelial gene:
     expected near-zero in melanoma and other non-epithelial tumours)
  4. positive controls: GAPDH (ubiquitous), PTPRC/CD45 (immune), EPCAM
     (epithelial), MLANA (melanocytic), GFAP (glial)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DER = ROOT / "data" / "derived"
OUT = ROOT / "results" / "claim_B5"
OUT.mkdir(parents=True, exist_ok=True)

CONTROLS = {
    "CLDN4": "ENSG00000189143",
    "CLDN3": "ENSG00000165215",
    "CLDN7": "ENSG00000181885",
    "GAPDH": "ENSG00000111640",
    "ACTB": "ENSG00000075624",
    "PTPRC": "ENSG00000081237",
    "EPCAM": "ENSG00000119888",
    "MLANA": "ENSG00000120215",
    "GFAP": "ENSG00000131095",
    "CD8A": "ENSG00000153563",
    "PDCD1": "ENSG00000188389",
    "CD274": "ENSG00000120217",
}


def load():
    expr = np.load(DER / "expr.npy")
    genes = (DER / "expr_genes.txt").read_text().split()
    samples = (DER / "expr_samples.txt").read_text().split("\n")
    samples = [s for s in samples if s]
    meta = pd.read_csv(ROOT / "data" / "raw" / "merged_metadata.tsv", sep="\t")
    return expr, genes, samples, meta


def main() -> None:
    expr, genes, samples, meta = load()
    gene_nov = [g.split(".")[0] for g in genes]
    idx = {g: i for i, g in enumerate(gene_nov)}

    lines: list[str] = []
    lines.append(f"expr matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
    lines.append(f"NaN fraction: {np.isnan(expr).mean():.6f}")
    lines.append(
        "value range: min %.3f  1%% %.3f  median %.3f  99%% %.3f  max %.3f"
        % (
            np.nanmin(expr),
            np.nanpercentile(expr, 1),
            np.nanmedian(expr),
            np.nanpercentile(expr, 99),
            np.nanmax(expr),
        )
    )
    # If values are log2(TPM + eps), then 2**x summed over genes ~ 1e6 per sample.
    colsum = np.nansum(np.power(2.0, expr[:, :5].astype(np.float64)), axis=0)
    lines.append("sum(2^x) for first 5 samples (1e6 => log2 TPM): " + ", ".join(f"{v:.3g}" for v in colsum))

    # sample id alignment
    meta_ids = set(meta["patientid"])
    overlap = len(meta_ids & set(samples))
    lines.append(f"metadata rows: {len(meta)}; expression columns: {len(samples)}; id overlap: {overlap}")

    lines.append("")
    lines.append("control gene availability and pan-compendium distribution (log2 units):")
    for sym, ens in CONTROLS.items():
        if ens not in idx:
            lines.append(f"  {sym:6s} {ens}  NOT FOUND")
            continue
        v = expr[idx[ens]]
        lines.append(
            f"  {sym:6s} {ens}  median {np.nanmedian(v):7.3f}  p10 {np.nanpercentile(v, 10):7.3f}"
            f"  p90 {np.nanpercentile(v, 90):7.3f}  max {np.nanmax(v):7.3f}"
        )

    # CLDN4 per study
    cldn4 = expr[idx[CONTROLS["CLDN4"]]]
    ptprc = expr[idx[CONTROLS["PTPRC"]]]
    df = pd.DataFrame({"patientid": samples, "CLDN4": cldn4, "PTPRC": ptprc}).merge(
        meta[["patientid", "study", "cancer_type", "response", "recist", "rna"]], on="patientid", how="left"
    )
    df.to_csv(DER / "cldn4_per_sample.tsv", sep="\t", index=False)

    floor = np.nanmin(expr)
    lines.append("")
    lines.append(f"expression floor (matrix minimum) = {floor:.4f}; treated as 'undetected'")
    lines.append("")
    lines.append("CLDN4 by study (log2 units; frac_at_floor = fraction of samples at the floor):")
    g = (
        df.groupby(["study", "rna"])
        .agg(
            n=("CLDN4", "size"),
            cldn4_median=("CLDN4", "median"),
            cldn4_p25=("CLDN4", lambda s: s.quantile(0.25)),
            cldn4_p75=("CLDN4", lambda s: s.quantile(0.75)),
            cldn4_iqr=("CLDN4", lambda s: s.quantile(0.75) - s.quantile(0.25)),
            frac_at_floor=("CLDN4", lambda s: float(np.mean(np.isclose(s, floor, atol=0.05)))),
            ptprc_median=("PTPRC", "median"),
        )
        .reset_index()
    )
    lines.append(g.to_string(index=False))

    text = "\n".join(lines) + "\n"
    (OUT / "qc_expression.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
