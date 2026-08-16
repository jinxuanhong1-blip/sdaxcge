"""Claim A10 / P4: does experimental NKX2-1 loss raise the ELF3 module and TACSTD2/CLDN4?

Contrasts are hard-coded from GEO sample titles / author column names — never
inferred from the claim genes themselves.

GSE188435 is FoxA1/2 deletion in an NKX2-1-positive background. It does not
manipulate NKX2-1 and is therefore not a P4 test; it is recorded as inspected
and excluded.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

GEO = os.path.join(C.RAW, "geo")


def mouse_id_to_symbol():
    annot = C.read_recount3_annotation(os.path.join(C.RAW, "mouse.gene_sums.M023.gtf.gz"))
    # strip version
    annot["ens"] = annot["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    return dict(zip(annot["ens"], annot["gene_name"]))


def collapse_symbol(df, symbol_col, value_cols) -> pd.DataFrame:
    sub = df[[symbol_col] + list(value_cols)].copy()
    sub[symbol_col] = sub[symbol_col].astype(str).str.split(",").str[0].str.strip()
    sub = sub[sub[symbol_col].ne("") & ~sub[symbol_col].isin(["nan", "-", "."])]
    for c in value_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    return sub.groupby(symbol_col, sort=False).mean(numeric_only=True)


def remap_ensembl(mat: pd.DataFrame, id2sym: dict) -> pd.DataFrame:
    idx = mat.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    names = idx.map(lambda x: id2sym.get(x, x))
    mat = mat.copy()
    mat.index = names
    order = mat.mean(axis=1).sort_values(ascending=False).index
    mat = mat.loc[order]
    return mat[~mat.index.duplicated(keep="first")]


def logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.divide(lib, axis=1) * 1e6 + 1.0)


def contrast(mat, lo, hi, dataset, comparison, genes):
    rows = []
    for g in genes:
        if g not in mat.index:
            rows.append({"dataset": dataset, "comparison": comparison, "gene": g,
                         "n_nkx_low": len(lo), "n_nkx_high": len(hi),
                         "detected": False, "log2fc_low_minus_high": np.nan,
                         "cliffs_delta": np.nan, "p_value": np.nan})
            continue
        a = mat.loc[g, lo].to_numpy(dtype=float)
        b = mat.loc[g, hi].to_numpy(dtype=float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) >= 2 and len(b) >= 2:
            try:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                p = np.nan
            d = C.cliffs_delta(a, b)
        else:
            p, d = np.nan, np.nan
        rows.append({"dataset": dataset, "comparison": comparison, "gene": g,
                     "n_nkx_low": int(len(a)), "n_nkx_high": int(len(b)),
                     "detected": True,
                     "mean_low": float(np.nanmean(a)) if len(a) else np.nan,
                     "mean_high": float(np.nanmean(b)) if len(b) else np.nan,
                     "log2fc_low_minus_high": float(np.nanmean(a) - np.nanmean(b)) if len(a) and len(b) else np.nan,
                     "cliffs_delta": d, "p_value": p})
    df = pd.DataFrame(rows)
    df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df


def gse129340():
    frames = []
    for label, fname in [
        ("GSE129340_H441", "GSE129340_H441_TTF1KD_genes.fpkm_tracking.gz"),
        ("GSE129340_H209", "GSE129340_H209_TTF1KD_genes.fpkm_tracking.gz"),
    ]:
        df = pd.read_csv(os.path.join(GEO, fname), sep="\t")
        cols = ["siNC_FPKM", "siTTF1#1_FPKM", "siTTF1#2_FPKM"]
        mat = collapse_symbol(df, "gene_short_name", cols)
        mat = np.log2(mat.clip(lower=0) + 1.0)
        frames.append(contrast(mat, ["siTTF1#1_FPKM", "siTTF1#2_FPKM"], ["siNC_FPKM"],
                               label, "siTTF-1 vs siNC (n=2 vs 1; MWU undefined)",
                               C.CLAIM_GENES_H))
    return frames


def gse229541():
    df = pd.read_csv(os.path.join(GEO, "GSE229541_NKX2-1_RNAseq_TPM_Final.csv.gz"))
    df = df.rename(columns={df.columns[0]: "gene"})
    num = [c for c in df.columns if c != "gene"]
    mat = collapse_symbol(df, "gene", num)
    mat = np.log2(mat.clip(lower=0) + 1.0)
    frames = []
    # Loss-of-function: low = sg/sh NKX2-1, high = non-targeting / shLuc
    specs = [
        ("GSE229541_H358_sg",
         [c for c in num if c.startswith("NCI-H358_sgNKX2-1_")],
         [c for c in num if c.startswith("NCI-H358_sgNT_")]),
        ("GSE229541_H2087_sg",
         [c for c in num if c.startswith("NCI-H2087_sgNKX2-1_")],
         [c for c in num if c.startswith("NCI-H2087_sgNT_")]),
        ("GSE229541_H441_sg",
         [c for c in num if c.startswith("NCI-H441_sgNKX2-1_")],
         [c for c in num if c.startswith("NCI-H441_sgNT_")]),
        ("GSE229541_H441_sh",
         [c for c in num if c.startswith("NCI-H441_shNKX2-1") and "Set2" not in c],
         [c for c in num if c.startswith("NCI-H441_shLuc") and "Set2" not in c]),
        ("GSE229541_H2087_sh",
         [c for c in num if c.startswith("NCI-H2087_shNKX2-1")],
         [c for c in num if c.startswith("NCI-H2087_shLuc")]),
        ("GSE229541_H441_sh_set2",
         [c for c in num if "H441_shNKX2-1" in c and "Set2" in c],
         [c for c in num if c.startswith("NCI-H441_shLuc_Set2")]),
        ("GSE229541_H2087_CRISPRi",
         [c for c in num if c.startswith("NCI-H2087_dCas9-KRAB-MeCP2_sgE764")],
         [c for c in num if c.startswith("NCI-H2087_dCas9-KRAB-MeCP2_sgNT")]),
        ("GSE229541_H358_CRISPRi",
         [c for c in num if c.startswith("NCI-H358_dCas9-KRAB-MeCP2_sgE764")],
         [c for c in num if c.startswith("NCI-H358_dCas9-KRAB-MeCP2_sgNT")]),
        ("GSE229541_HCC78_CRISPRi",
         [c for c in num if c.startswith("HCC78_dCas9-KRAB-MeCP2_sgE764")],
         [c for c in num if c.startswith("HCC78_dCas9-KRAB-MeCP2_sgNT")]),
    ]
    for name, lo, hi in specs:
        if lo and hi:
            frames.append(contrast(mat, lo, hi, name, "NKX2-1 loss vs control", C.CLAIM_GENES_H))
    # Gain-of-function: claim predicts MODULE/TARGET go DOWN when NKX2-1 is raised,
    # so we still report low=GFP (endogenous) minus high=OE. Positive log2FC = claim.
    oe = [
        ("GSE229541_PC9_OE",
         [c for c in num if c.startswith("PC9_EF1a-GFP_Rep")],
         [c for c in num if c.startswith("PC9_EF1a-NKX2-1_Rep") or c.startswith("PC9_hPGK-NKX2-1_Rep")]),
        ("GSE229541_H1975_OE",
         [c for c in num if c.startswith("NCI-H1975_EF1a-GFP_Rep")],
         [c for c in num if c.startswith("NCI-H1975_EF1a-NKX2-1_Rep") or c.startswith("NCI-H1975_hPGK-NKX2-1_Rep")]),
    ]
    for name, lo, hi in oe:
        if lo and hi:
            frames.append(contrast(mat, lo, hi, name,
                                   "parental/GFP vs NKX2-1 overexpression",
                                   C.CLAIM_GENES_H))
    return frames


def gse129583():
    df = pd.ExcelFile(os.path.join(GEO, "GSE129583_processeddata_matrix_Nkx2-1_RNAseq.xlsx")).parse("Sheet1")
    frames = []
    blocks = [
        ("GSE129583_AT1_P5", "gene_short_name",
         ["Nkx2-1;Aqp5 P5 AT1 Mutant_A", "Nkx2-1;Aqp5 P5 AT1 Mutant_B", "Nkx2-1;Aqp5 P5 AT1 Mutant_C"],
         ["Nkx2-1;Aqp5 P5 AT1 Control_A", "Nkx2-1;Aqp5 P5 AT1 Control_B", "Nkx2-1;Aqp5 P5 AT1 Control_C"]),
        ("GSE129583_AT2_P8P9", "gene_short_name",
         ["Nkx2-1; Sftpc P8 AT2 Mutant_A", "Nkx2-1; Sftpc P8 AT2 Mutant_B", "Nkx2-1; Sftpc P9 AT2 Mutant_C"],
         ["Nkx2-1; Sftpc P8 AT2 Control_A", "Nkx2-1; Sftpc P8 AT2 Control_B", "Nkx2-1; Sftpc P9 AT2 Control_C"]),
    ]
    for name, sym, lo, hi in blocks:
        mat = collapse_symbol(df, sym, lo + hi)
        # symbols are mixed case; normalise to mouse claim names
        mat.index = mat.index.astype(str)
        lower = {i.lower(): i for i in mat.index}
        keep = {}
        for g in C.CLAIM_GENES_M + list(C.MARKERS_M):
            if g.lower() in lower:
                keep[g] = lower[g.lower()]
        mat = mat.loc[list(keep.values())]
        mat.index = list(keep.keys())
        mat = np.log2(mat.clip(lower=0) + 1.0)
        frames.append(contrast(mat, lo, hi, name, "Nkx2-1 mutant vs control",
                               C.CLAIM_GENES_M))
    return frames


def gse115899(id2sym):
    df = pd.read_csv(os.path.join(GEO, "GSE115899_counts_13438R.txt.gz"), sep="\t")
    df = df.rename(columns={df.columns[0]: "gene"})
    # drop unnamed empty columns
    keep_cols = [c for c in df.columns if not str(c).startswith("Unnamed")]
    df = df[keep_cols]
    num = [c for c in df.columns if c != "gene"]
    mat = collapse_symbol(df, "gene", num)
    mat = remap_ensembl(mat, id2sym)
    mat = logcpm(mat)
    lo = ["13438X5", "13438X6", "13438X7"]  # NKX2-1-negative
    hi = ["13438X1", "13438X2", "13438X3"]  # NKX2-1-positive
    return [contrast(mat, lo, hi, "GSE115899",
                     "Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3)",
                     C.CLAIM_GENES_M)]


def gse145152(id2sym):
    df = pd.read_csv(os.path.join(GEO, "GSE145152_counts_14489R.txt.gz"), sep="\t")
    df = df.rename(columns={df.columns[0]: "gene"})
    num = [c for c in df.columns if c != "gene"]
    mat = collapse_symbol(df, "gene", num)
    mat = remap_ensembl(mat, id2sym)
    mat = logcpm(mat)
    # control-chow only: BPN (Nkx2-1 f/f) vs BP (Nkx2-1 f/+)
    lo = ["14489X6", "14489X7", "14489X8", "14489X9"]
    hi = ["14489X1", "14489X2", "14489X3", "14489X4", "14489X5"]
    return [contrast(mat, lo, hi, "GSE145152",
                     "Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow",
                     C.CLAIM_GENES_M)]


def main():
    id2sym = mouse_id_to_symbol()
    frames = []
    frames += gse129340()
    frames += gse229541()
    frames += gse129583()
    frames += gse115899(id2sym)
    frames += gse145152(id2sym)
    out = pd.concat(frames, ignore_index=True)
    C.write_table(out, "perturbation_contrasts.csv")
    notes = [
        "GSE188435 inspected and EXCLUDED: FoxA1/2 deletion in NKX2-1-positive LUAD; NKX2-1 is not the perturbed factor.",
        "GSE145152 MAPK-inhibitor arms excluded so the contrast is genotype-only.",
        "GSE229541 osimertinib arms excluded (drug confound).",
        "GSE129340 has n=1 control; log2FC is reported, MWU p-values are undefined.",
    ]
    with open(os.path.join(C.LOGS, "perturbation_notes.txt"), "w") as fh:
        fh.write("\n".join(notes) + "\n")
        fh.write(f"contrasts written: {out.dataset.nunique()}\n")
    print(out.groupby("dataset").size().to_string())
    print("done", flush=True)


if __name__ == "__main__":
    main()
